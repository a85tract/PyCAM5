module stratiform

!-------------------------------------------------------------------------------------------------------
!
! Provides the CAM interface to the Rasch and Kristjansson (RK) 
! prognostic cloud microphysics, and the cam3/4 macrophysics.
!
!-------------------------------------------------------------------------------------------------------

use shr_kind_mod,   only: r8=>shr_kind_r8
use ppgrid,         only: pcols, pver, pverp
use physconst,      only: gravit, latvap, latice
use phys_control,   only: phys_getopts
use constituents,   only: pcnst

use cam_logfile,    only: iulog
use cam_abortutils, only: endrun
use perf_mod
use spmd_utils,     only: masterproc

implicit none
private
save

public :: stratiform_register, stratiform_init_cnst, stratiform_implements_cnst
public :: stratiform_init
public :: stratiform_tend

! Physics buffer indices 
integer  ::  qcwat_idx          = 0 
integer  ::  lcwat_idx          = 0 
integer  ::  tcwat_idx          = 0 

integer  ::  cld_idx            = 0 
integer  ::  ast_idx            = 0 
integer  ::  concld_idx         = 0 
integer  ::  fice_idx           = 0 

integer  ::  qme_idx            = 0 
integer  ::  prain_idx          = 0 
integer  ::  nevapr_idx         = 0 

integer  ::  wsedl_idx          = 0

integer  ::  rei_idx            = 0 
integer  ::  rel_idx            = 0 

integer  ::  shfrc_idx          = 0

integer  ::  prec_str_idx       = 0
integer  ::  snow_str_idx       = 0
integer  ::  prec_sed_idx       = 0
integer  ::  snow_sed_idx       = 0
integer  ::  prec_pcw_idx       = 0
integer  ::  snow_pcw_idx       = 0

integer  ::  ls_flxprc_idx      = 0
integer  ::  ls_flxsnw_idx      = 0

integer, parameter :: ncnst = 2                    ! Number of constituents
character(len=8), dimension(ncnst), parameter &    ! Constituent names
                   :: cnst_names = (/'CLDLIQ', 'CLDICE'/)
logical            :: use_shfrc                       ! Local copy of flag from convect_shallow_use_shfrc

logical            :: do_cnst = .false. ! True when this module has registered constituents.
logical            :: use_native_impl = .false.
logical            :: impl_selected = .false.
logical            :: stratiform_entered_logged = .false.
logical            :: use_native_microphys_shell_impl = .false.
logical            :: microphys_shell_impl_selected = .false.
logical            :: microphys_shell_entered_logged = .false.
logical            :: detrain_assign_entered_logged = .false.
logical            :: sedimentation_diag_entered_logged = .false.
logical            :: rhdfda_entered_logged = .false.
logical            :: repartition_entered_logged = .false.
logical            :: forcing_prep_entered_logged = .false.
logical            :: microphys_tend_diag_entered_logged = .false.
logical            :: cloud_mixing_diag_entered_logged = .false.
logical            :: postcloud_diag_entered_logged = .false.
logical            :: store_oldcloud_entered_logged = .false.
integer            :: branch_mask = 0
logical            :: branch_selected = .false.

integer :: &
   ixcldliq,     &! cloud liquid amount index
   ixcldice       ! cloud ice amount index

!===============================================================================
contains
!===============================================================================

subroutine stratiform_register

   !---------------------------------------------------------------------- !
   !                                                                       !
   ! Register the constituents (cloud liquid and cloud ice) and the fields !
   ! in the physics buffer.                                                !
   !                                                                       !
   !---------------------------------------------------------------------- !

   use constituents, only: cnst_add, pcnst
   use physconst,    only: mwh2o, cpair
    
   use physics_buffer, only : pbuf_add_field, dtype_r8, dyn_time_lvls

   !-----------------------------------------------------------------------

   ! Take note of the fact that we are registering constituents.
   do_cnst = .true.

   ! Register cloud water and save indices.
   call cnst_add(cnst_names(1), mwh2o, cpair, 0._r8, ixcldliq, &
      longname='Grid box averaged cloud liquid amount', is_convtran1=.true.)
   call cnst_add(cnst_names(2), mwh2o, cpair, 0._r8, ixcldice, &
      longname='Grid box averaged cloud ice amount', is_convtran1=.true.)

   call pbuf_add_field('QCWAT',  'global', dtype_r8, (/pcols,pver,dyn_time_lvls/), qcwat_idx)
   call pbuf_add_field('LCWAT',  'global', dtype_r8, (/pcols,pver,dyn_time_lvls/), lcwat_idx)
   call pbuf_add_field('TCWAT',  'global', dtype_r8, (/pcols,pver,dyn_time_lvls/), tcwat_idx)

   call pbuf_add_field('CLD',    'global', dtype_r8, (/pcols,pver,dyn_time_lvls/), cld_idx)
   call pbuf_add_field('AST',    'global', dtype_r8, (/pcols,pver,dyn_time_lvls/), ast_idx)
   call pbuf_add_field('CONCLD', 'global', dtype_r8, (/pcols,pver,dyn_time_lvls/), concld_idx)

   call pbuf_add_field('FICE',   'physpkg', dtype_r8, (/pcols,pver/), fice_idx) 

   call pbuf_add_field('QME',       'physpkg', dtype_r8, (/pcols,pver/), qme_idx)
   call pbuf_add_field('PRAIN',     'physpkg', dtype_r8, (/pcols,pver/), prain_idx)
   call pbuf_add_field('NEVAPR',    'physpkg', dtype_r8, (/pcols,pver/), nevapr_idx)

   call pbuf_add_field('WSEDL',     'physpkg', dtype_r8, (/pcols,pver/), wsedl_idx)

   call pbuf_add_field('REI',       'physpkg', dtype_r8, (/pcols,pver/), rei_idx)
   call pbuf_add_field('REL',       'physpkg', dtype_r8, (/pcols,pver/), rel_idx)

   call pbuf_add_field('LS_FLXPRC', 'physpkg', dtype_r8, (/pcols,pverp/), ls_flxprc_idx)
   call pbuf_add_field('LS_FLXSNW', 'physpkg', dtype_r8, (/pcols,pverp/), ls_flxsnw_idx)

end subroutine stratiform_register

!===============================================================================

function stratiform_implements_cnst(name)

  !----------------------------------------------------------------------------- ! 
  !                                                                              !    
  ! Return true if specified constituent is implemented by this package          !
  !                                                                              !
  !----------------------------------------------------------------------------- !

   character(len=*), intent(in) :: name      ! constituent name
   logical :: stratiform_implements_cnst     ! return value

   !-----------------------------------------------------------------------

   stratiform_implements_cnst = (do_cnst .and. any(name == cnst_names))

end function stratiform_implements_cnst

!===============================================================================

subroutine stratiform_init_cnst(name, q, gcid)

   !----------------------------------------------------------------------- !
   !                                                                        !
   ! Initialize the cloud water mixing ratios (liquid and ice), if they are !
   ! not read from the initial file                                         ! 
   !                                                                        !
   !----------------------------------------------------------------------- !

   character(len=*), intent(in)  :: name     ! constituent name
   real(r8),         intent(out) :: q(:,:)   ! mass mixing ratio (gcol, plev)
   integer,          intent(in)  :: gcid(:)  ! global column id
   !-----------------------------------------------------------------------

   if (any(name == cnst_names)) q = 0.0_r8

end subroutine stratiform_init_cnst

!===============================================================================

subroutine stratiform_init()

   !-------------------------------------------- !
   !                                             !
   ! Initialize the cloud water parameterization !
   !                                             ! 
   !-------------------------------------------- !

   use physics_buffer,  only: physics_buffer_desc, pbuf_get_index
   use constituents,    only: cnst_get_ind, cnst_name, cnst_longname, sflxnam, apcnst, bpcnst
   use cam_history,     only: addfld, add_default, phys_decomp
   use convect_shallow, only: convect_shallow_use_shfrc
   use phys_control,    only: cam_physpkg_is
   use physconst,       only: tmelt, rhodair, rh2o
   use cldwat,          only: inimc
    
   integer :: m, mm
   logical :: history_amwg         ! output the variables used by the AMWG diag package
   logical :: history_aerosol      ! Output the MAM aerosol tendencies
   logical :: history_budget       ! Output tendencies and state variables for CAM4
                                   ! temperature, water vapor, cloud ice and cloud
                                   ! liquid budgets.
   integer :: history_budget_histfile_num ! output history file number for budget fields
   !-----------------------------------------------------------------------

   call phys_getopts( history_aerosol_out        = history_aerosol      , &
                      history_amwg_out   = history_amwg                 , & 
                      history_budget_out         = history_budget       , &
                      history_budget_histfile_num_out = history_budget_histfile_num)

   ! Find out whether shfrc from convect_shallow will be used in cldfrc
   if( convect_shallow_use_shfrc() ) then
      use_shfrc = .true.
      shfrc_idx = pbuf_get_index('shfrc')
   else 
      use_shfrc = .false.
   endif

   ! Register history variables

   do m = 1, ncnst
      call cnst_get_ind( cnst_names(m), mm )
      call addfld( cnst_name(mm), 'kg/kg   ', pver, 'A', cnst_longname(mm)                   , phys_decomp )
      call addfld( sflxnam  (mm), 'kg/m2/s ',    1, 'A', trim(cnst_name(mm))//' surface flux', phys_decomp )
      if (history_amwg) then
         call add_default( cnst_name(mm), 1, ' ' )
         call add_default( sflxnam  (mm), 1, ' ' )
      endif
   enddo

   call addfld (apcnst(ixcldliq), 'kg/kg   ', pver, 'A', trim(cnst_name(ixcldliq))//' after physics'  , phys_decomp)
   call addfld (apcnst(ixcldice), 'kg/kg   ', pver, 'A', trim(cnst_name(ixcldice))//' after physics'  , phys_decomp)
   call addfld (bpcnst(ixcldliq), 'kg/kg   ', pver, 'A', trim(cnst_name(ixcldliq))//' before physics' , phys_decomp)
   call addfld (bpcnst(ixcldice), 'kg/kg   ', pver, 'A', trim(cnst_name(ixcldice))//' before physics' , phys_decomp)

   if( history_budget) then
      call add_default (cnst_name(ixcldliq), history_budget_histfile_num, ' ')
      call add_default (cnst_name(ixcldice), history_budget_histfile_num, ' ')
      call add_default (apcnst   (ixcldliq), history_budget_histfile_num, ' ')
      call add_default (apcnst   (ixcldice), history_budget_histfile_num, ' ')
      call add_default (bpcnst   (ixcldliq), history_budget_histfile_num, ' ')
      call add_default (bpcnst   (ixcldice), history_budget_histfile_num, ' ')
   end if

   call addfld ('FWAUT    ', 'fraction', pver, 'A', 'Relative importance of liquid autoconversion'            ,phys_decomp)
   call addfld ('FSAUT    ', 'fraction', pver, 'A', 'Relative importance of ice autoconversion'               ,phys_decomp)
   call addfld ('FRACW    ', 'fraction', pver, 'A', 'Relative importance of rain accreting liquid'            ,phys_decomp)
   call addfld ('FSACW    ', 'fraction', pver, 'A', 'Relative importance of snow accreting liquid'            ,phys_decomp)
   call addfld ('FSACI    ', 'fraction', pver, 'A', 'Relative importance of snow accreting ice'               ,phys_decomp)
   call addfld ('CME      ', 'kg/kg/s ', pver, 'A', 'Rate of cond-evap within the cloud'                      ,phys_decomp)
   call addfld ('CMEICE   ', 'kg/kg/s ', pver, 'A', 'Rate of cond-evap of ice within the cloud'               ,phys_decomp)
   call addfld ('CMELIQ   ', 'kg/kg/s ', pver, 'A', 'Rate of cond-evap of liq within the cloud'               ,phys_decomp)
   call addfld ('ICE2PR   ', 'kg/kg/s ', pver, 'A', 'Rate of conversion of ice to precip'                     ,phys_decomp)
   call addfld ('LIQ2PR   ', 'kg/kg/s ', pver, 'A', 'Rate of conversion of liq to precip'                     ,phys_decomp)
   call addfld ('ZMDLF    ', 'kg/kg/s ', pver, 'A', 'Detrained liquid water from ZM convection'               ,phys_decomp)
   call addfld ('SHDLF    ', 'kg/kg/s ', pver, 'A', 'Detrained liquid water from shallow convection'          ,phys_decomp)

   call addfld ('PRODPREC ', 'kg/kg/s ', pver, 'A', 'Rate of conversion of condensate to precip'              ,phys_decomp)
   call addfld ('EVAPPREC ', 'kg/kg/s ', pver, 'A', 'Rate of evaporation of falling precip'                   ,phys_decomp)
   call addfld ('EVAPSNOW ', 'kg/kg/s ', pver, 'A', 'Rate of evaporation of falling snow'                     ,phys_decomp)
   call addfld ('HPROGCLD ', 'W/kg'    , pver, 'A', 'Heating from prognostic clouds'                          ,phys_decomp)
   call addfld ('HCME     ', 'W/kg'    , pver, 'A', 'Heating from cond-evap within the cloud'                 ,phys_decomp)
   call addfld ('HEVAP    ', 'W/kg'    , pver, 'A', 'Heating from evaporation of falling precip'              ,phys_decomp)
   call addfld ('HFREEZ   ', 'W/kg'    , pver, 'A', 'Heating rate due to freezing of precip'                  ,phys_decomp)
   call addfld ('HMELT    ', 'W/kg'    , pver, 'A', 'Heating from snow melt'                                  ,phys_decomp)
   call addfld ('HREPART  ', 'W/kg'    , pver, 'A', 'Heating from cloud ice/liquid repartitioning'            ,phys_decomp)
   call addfld ('REPARTICE', 'kg/kg/s' , pver, 'A', 'Cloud ice tendency from cloud ice/liquid repartitioning' ,phys_decomp)
   call addfld ('REPARTLIQ', 'kg/kg/s' , pver, 'A', 'Cloud liq tendency from cloud ice/liquid repartitioning' ,phys_decomp)
   call addfld ('FICE     ', 'fraction', pver, 'A', 'Fractional ice content within cloud'                     ,phys_decomp)
   call addfld ('ICWMR    ', 'kg/kg   ', pver, 'A', 'Prognostic in-cloud water mixing ratio'                  ,phys_decomp)
   call addfld ('ICIMR    ', 'kg/kg   ', pver, 'A', 'Prognostic in-cloud ice mixing ratio'                    ,phys_decomp)
   call addfld ('PCSNOW   ', 'm/s     ', 1   , 'A', 'Snow fall from prognostic clouds'                        ,phys_decomp)
 
   call addfld ('DQSED    ', 'kg/kg/s ', pver, 'A', 'Water vapor tendency from cloud sedimentation'           ,phys_decomp)
   call addfld ('DLSED    ', 'kg/kg/s ', pver, 'A', 'Cloud liquid tendency from sedimentation'                ,phys_decomp)
   call addfld ('DISED    ', 'kg/kg/s ', pver, 'A', 'Cloud ice tendency from sedimentation'                   ,phys_decomp)
   call addfld ('HSED     ', 'W/kg    ', pver, 'A', 'Heating from cloud sediment evaporation'                 ,phys_decomp)
   call addfld ('SNOWSED  ', 'm/s     ', 1   , 'A', 'Snow from cloud ice sedimentation'                       ,phys_decomp)
   call addfld ('RAINSED  ', 'm/s     ', 1   , 'A', 'Rain from cloud liquid sedimentation'                    ,phys_decomp)
   call addfld ('PRECSED  ', 'm/s     ', 1   , 'A', 'Precipitation from cloud sedimentation'                  ,phys_decomp)


   call addfld ('CNVCLD   ', 'fraction', 1,    'A', 'Vertically integrated convective cloud amount'           ,phys_decomp)
   call addfld ('CLDST    ', 'fraction', pver, 'A', 'Stratus cloud fraction'                                  ,phys_decomp)
   call addfld ('CONCLD   ', 'fraction', pver, 'A', 'Convective cloud cover'                                  ,phys_decomp)
	
   call addfld ('AST'      ,'fraction' , pver, 'A', 'Stratus cloud fraction'                                  ,phys_decomp)
   call addfld ('LIQCLDF  ', 'fraction', pver, 'A', 'Stratus Liquid cloud fraction'                           ,phys_decomp)
   call addfld ('ICECLDF  ', 'fraction', pver, 'A', 'Stratus ICE cloud fraction'                              ,phys_decomp)
   call addfld ('IWC      ', 'kg/m3   ', pver, 'A', 'Grid box average ice water content'                      ,phys_decomp)
   call addfld ('LWC      ', 'kg/m3   ', pver, 'A', 'Grid box average liquid water content'                   ,phys_decomp)
   call addfld ('ICWNC    ', 'm-3     ', pver, 'A', 'Prognostic in-cloud water number conc'                   ,phys_decomp)
   call addfld ('ICINC    ', 'm-3     ', pver, 'A', 'Prognostic in-cloud ice number conc'                     ,phys_decomp)
   call addfld ('EFFLIQ   ', 'Micron  ', pver, 'A', 'Prognostic droplet effective radius'                     ,phys_decomp)
   call addfld ('EFFLIQ_IND','Micron  ', pver, 'A', 'Prognostic droplet effective radius (indirect effect)'   ,phys_decomp)
   call addfld ('EFFICE   ', 'Micron  ', pver, 'A', 'Prognostic ice effective radius'                         ,phys_decomp)
   call addfld ('REL',       'micron',   pver, 'A', 'effective liquid drop radius'                            ,phys_decomp)
   call addfld ('REI',       'micron',   pver, 'A', 'effective ice particle radius'                           ,phys_decomp)

   if ( history_budget ) then

      call add_default ('EVAPSNOW ', history_budget_histfile_num, ' ')
      call add_default ('EVAPPREC ', history_budget_histfile_num, ' ')
      call add_default ('CMELIQ   ', history_budget_histfile_num, ' ')

      if( cam_physpkg_is('cam3') .or. cam_physpkg_is('cam4') ) then

         call add_default ('ZMDLF    ', history_budget_histfile_num, ' ')
         call add_default ('CME      ', history_budget_histfile_num, ' ')
         call add_default ('DQSED    ', history_budget_histfile_num, ' ')
         call add_default ('DISED    ', history_budget_histfile_num, ' ')
         call add_default ('DLSED    ', history_budget_histfile_num, ' ')
         call add_default ('HSED     ', history_budget_histfile_num, ' ')
         call add_default ('CMEICE   ', history_budget_histfile_num, ' ')
         call add_default ('LIQ2PR   ', history_budget_histfile_num, ' ')
         call add_default ('ICE2PR   ', history_budget_histfile_num, ' ')
         call add_default ('HCME     ', history_budget_histfile_num, ' ')
         call add_default ('HEVAP    ', history_budget_histfile_num, ' ')
         call add_default ('HFREEZ   ', history_budget_histfile_num, ' ')
         call add_default ('HMELT    ', history_budget_histfile_num, ' ')
         call add_default ('HREPART  ', history_budget_histfile_num, ' ')
         call add_default ('HPROGCLD ', history_budget_histfile_num, ' ')
         call add_default ('REPARTLIQ', history_budget_histfile_num, ' ')
         call add_default ('REPARTICE', history_budget_histfile_num, ' ')

      end if

   end if

   if (history_amwg) then
      call add_default ('ICWMR', 1, ' ')
      call add_default ('ICIMR', 1, ' ')
      call add_default ('CONCLD  ', 1, ' ')
      call add_default ('FICE    ', 1, ' ')
   endif

   ! History Variables for COSP/CFMIP
   call addfld ('LS_FLXPRC', 'kg/m2/s', pverp, 'A', 'ls stratiform gbm interface rain+snow flux', phys_decomp)
   call addfld ('LS_FLXSNW', 'kg/m2/s', pverp, 'A', 'ls stratiform gbm interface snow flux', phys_decomp)
   call addfld ('PRACWO', '1/s', pver, 'A', 'Accretion of cloud water by rain', phys_decomp)
   call addfld ('PSACWO', '1/s', pver, 'A', 'Accretion of cloud water by snow', phys_decomp)
   call addfld ('PSACIO', '1/s', pver, 'A', 'Accretion of cloud ice by snow', phys_decomp)

   call addfld ('CLDLIQSTR   ', 'kg/kg', pver, 'A', 'Stratiform CLDLIQ', phys_decomp)
   call addfld ('CLDICESTR   ', 'kg/kg', pver, 'A', 'Stratiform CLDICE', phys_decomp)
   call addfld ('CLDLIQCON   ', 'kg/kg', pver, 'A', 'Convective CLDLIQ', phys_decomp)
   call addfld ('CLDICECON   ', 'kg/kg', pver, 'A', 'Convective CLDICE', phys_decomp)

   prec_str_idx = pbuf_get_index('PREC_STR')
   snow_str_idx = pbuf_get_index('SNOW_STR')
   prec_pcw_idx = pbuf_get_index('PREC_PCW')
   snow_pcw_idx = pbuf_get_index('SNOW_PCW')
   prec_sed_idx = pbuf_get_index('PREC_SED')
   snow_sed_idx = pbuf_get_index('SNOW_SED')

   ! Initialize cldwat with constants.
   call inimc(tmelt, rhodair/1000.0_r8, gravit, rh2o)

end subroutine stratiform_init

!===============================================================================

subroutine stratiform_tend( &
   state, ptend_all, pbuf, dtime, icefrac, &
   landfrac, ocnfrac, landm, snowh, dlf,   &
   dlf2, rliq, cmfmc, cmfmc2, ts,          &
   sst, zdu)

   !-------------------------------------------------------- !  
   !                                                         ! 
   ! Interface to sedimentation, detrain, cloud fraction and !
   !        cloud macro - microphysics subroutines           !
   !                                                         ! 
   !-------------------------------------------------------- !

   use cloud_fraction,   only: cldfrc, cldfrc_fice
   use physics_types,    only: physics_state, physics_ptend
   use physics_types,    only: physics_ptend_init, physics_update
   use physics_types,    only: physics_ptend_sum,  physics_state_copy
   use physics_types,    only: physics_state_dealloc
   use cam_history,      only: outfld
   use physics_buffer,   only: physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx
   use pkg_cld_sediment, only: cld_sediment_vel, cld_sediment_tend
   use cldwat,           only: pcond
   use pkg_cldoptics,    only: cldefr
   use phys_control,     only: cam_physpkg_is

   ! Arguments
   type(physics_state), intent(in)    :: state       ! State variables
   type(physics_ptend), intent(out)   :: ptend_all   ! Package tendencies
   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8), intent(in)  :: dtime                    ! Timestep
   real(r8), intent(in)  :: icefrac (pcols)          ! Sea ice fraction (fraction)
   real(r8), intent(in)  :: landfrac(pcols)          ! Land fraction (fraction)
   real(r8), intent(in)  :: ocnfrac (pcols)          ! Ocean fraction (fraction)
   real(r8), intent(in)  :: landm(pcols)             ! Land fraction ramped over water
   real(r8), intent(in)  :: snowh(pcols)             ! Snow depth over land, water equivalent (m)

   real(r8), intent(in)  :: dlf(pcols,pver)          ! Detrained water from convection schemes
   real(r8), intent(in)  :: dlf2(pcols,pver)         ! Detrained water from shallow convection scheme
   real(r8), intent(in)  :: rliq(pcols)              ! Vertical integral of liquid not yet in q(ixcldliq)
   real(r8), intent(in)  :: cmfmc(pcols,pverp)       ! Deep + Shallow Convective mass flux [ kg /s/m^2 ]
   real(r8), intent(in)  :: cmfmc2(pcols,pverp)      ! Shallow convective mass flux [ kg/s/m^2 ]

   real(r8), intent(in)  :: ts(pcols)                ! Surface temperature
   real(r8), intent(in)  :: sst(pcols)               ! Sea surface temperature
   real(r8), intent(in)  :: zdu(pcols,pver)          ! Detrainment rate from deep convection

  ! Local variables

   type(physics_state)   :: state1                   ! Local copy of the state variable
   type(physics_ptend)   :: ptend_loc                ! Package tendencies

   integer :: i, k, m
   integer :: lchnk                                  ! Chunk identifier
   integer :: ncol                                   ! Number of atmospheric columns
   integer :: itim_old

   ! Physics buffer fields
   real(r8), pointer :: prec_str(:)          ! [Total] Sfc flux of precip from stratiform [ m/s ] 
   real(r8), pointer :: snow_str(:)          ! [Total] Sfc flux of snow from stratiform   [ m/s ]
   real(r8), pointer :: prec_sed(:)          ! Surface flux of total cloud water from sedimentation
   real(r8), pointer :: snow_sed(:)          ! Surface flux of cloud ice from sedimentation
   real(r8), pointer :: prec_pcw(:)          ! Sfc flux of precip from microphysics [ m/s ]
   real(r8), pointer :: snow_pcw(:)          ! Sfc flux of snow from microphysics [ m/s ]

   real(r8), pointer, dimension(:,:) :: qcwat        ! Cloud water old q
   real(r8), pointer, dimension(:,:) :: tcwat        ! Cloud water old temperature
   real(r8), pointer, dimension(:,:) :: lcwat        ! Cloud liquid water old q
   real(r8), pointer, dimension(:,:) :: cld          ! Total cloud fraction
   real(r8), pointer, dimension(:,:) :: fice         ! Cloud ice/water partitioning ratio.
   real(r8), pointer, dimension(:,:) :: ast          ! Relative humidity cloud fraction
   real(r8), pointer, dimension(:,:) :: concld       ! Convective cloud fraction
   real(r8), pointer, dimension(:,:) :: qme          ! rate of cond-evap of condensate (1/s)
   real(r8), pointer, dimension(:,:) :: prain        ! Total precipitation (rain + snow)
   real(r8), pointer, dimension(:,:) :: nevapr       ! Evaporation of total precipitation (rain + snow)
   real(r8), pointer, dimension(:,:) :: rel          ! Liquid effective drop radius (microns)
   real(r8), pointer, dimension(:,:) :: rei          ! Ice effective drop size (microns)
   real(r8), pointer, dimension(:,:) :: wsedl        ! Sedimentation velocity of liquid stratus cloud droplet [ m/s ]
   real(r8), pointer, dimension(:,:) :: shfrc        ! Cloud fraction from shallow convection scheme

   real(r8), target :: shfrc_local(pcols,pver)

   ! physics buffer fields for COSP simulator (RK only)
   real(r8), pointer, dimension(:,:) :: rkflxprc     ! RK grid-box mean flux_large_scale_cloud_rain+snow at interfaces (kg/m2/s)
   real(r8), pointer, dimension(:,:) :: rkflxsnw     ! RK grid-box mean flux_large_scale_cloud_snow at interfaces (kg/m2/s)

   ! Local variables for stratiform_sediment
   real(r8) :: rain(pcols)                           ! Surface flux of cloud liquid
   real(r8) :: pvliq(pcols,pverp)                    ! Vertical velocity of cloud liquid drops (Pa/s)
   real(r8) :: pvice(pcols,pverp)                    ! Vertical velocity of cloud ice particles (Pa/s)

   ! Local variables for cldfrc

   real(r8) :: cldst(pcols,pver)                       ! Stratus cloud fraction
   real(r8) :: rhcloud(pcols,pver)                     ! Relative humidity cloud (last timestep)
   real(r8) :: rhcloud2(pcols,pver)                    ! Relative humidity cloud (perturbation)
   real(r8) :: clc(pcols)                              ! Column convective cloud amount
   real(r8) :: relhum(pcols,pver)                      ! RH, output to determine drh/da
   real(r8) :: rhu00(pcols,pver)
   real(r8) :: rhu002(pcols,pver)                      ! Same as rhu00 but for perturbed rh 
   real(r8) :: rhdfda(pcols,pver)
   real(r8) :: cld2(pcols,pver)                        ! Same as cld but for perturbed rh
   real(r8) :: concld2(pcols,pver)                     ! Same as concld but for perturbed rh 
   real(r8) :: cldst2(pcols,pver)                      ! Same as cldst but for perturbed rh 
   real(r8) :: relhum2(pcols,pver)                     ! RH after  perturbation            
   real(r8) :: icecldf(pcols,pver)                     ! Ice cloud fraction
   real(r8) :: liqcldf(pcols,pver)                     ! Liquid cloud fraction (combined into cloud)
   real(r8) :: icecldf_out(pcols,pver)                 ! Ice cloud fraction
   real(r8) :: liqcldf_out(pcols,pver)                 ! Liquid cloud fraction (combined into cloud)
   real(r8) :: icecldf2(pcols,pver)                    ! Ice cloud fraction
   real(r8) :: liqcldf2(pcols,pver)                    ! Liquid cloud fraction (combined into cloud)

   ! Local variables for microphysics

   real(r8) :: rdtime                                  ! 1./dtime
   real(r8) :: qtend(pcols,pver)                       ! Moisture tendencies
   real(r8) :: ttend(pcols,pver)                       ! Temperature tendencies
   real(r8) :: ltend(pcols,pver)                       ! Cloud liquid water tendencies
   real(r8) :: evapheat(pcols,pver)                    ! Heating rate due to evaporation of precip
   real(r8) :: evapsnow(pcols,pver)                    ! Local evaporation of snow
   real(r8) :: prfzheat(pcols,pver)                    ! Heating rate due to freezing of precip (W/kg)
   real(r8) :: meltheat(pcols,pver)                    ! Heating rate due to phase change of precip
   real(r8) :: cmeheat (pcols,pver)                    ! Heating rate due to phase change of precip
   real(r8) :: prodsnow(pcols,pver)                    ! Local production of snow
   real(r8) :: totcw(pcols,pver)                       ! Total cloud water mixing ratio
   real(r8) :: fsnow(pcols,pver)                       ! Fractional snow production
   real(r8) :: repartht(pcols,pver)                    ! Heating rate due to phase repartition of input precip
   real(r8) :: icimr(pcols,pver)                       ! In cloud ice mixing ratio
   real(r8) :: icwmr(pcols,pver)                       ! In cloud water mixing ratio
   real(r8) :: fwaut(pcols,pver)              
   real(r8) :: fsaut(pcols,pver)              
   real(r8) :: fracw(pcols,pver)              
   real(r8) :: fsacw(pcols,pver)              
   real(r8) :: fsaci(pcols,pver)              
   real(r8) :: cmeice(pcols,pver)                      ! Rate of cond-evap of ice within the cloud
   real(r8) :: cmeliq(pcols,pver)                      ! Rate of cond-evap of liq within the cloud
   real(r8) :: ice2pr(pcols,pver)                      ! Rate of conversion of ice to precip
   real(r8) :: liq2pr(pcols,pver)                      ! Rate of conversion of liquid to precip
   real(r8) :: liq2snow(pcols,pver)                    ! Rate of conversion of liquid to snow

   ! Local variables for CFMIP calculations
   real(r8) :: mr_lsliq(pcols,pver)                     ! mixing_ratio_large_scale_cloud_liquid (kg/kg)
   real(r8) :: mr_lsice(pcols,pver)                     ! mixing_ratio_large_scale_cloud_ice (kg/kg)
   real(r8) :: mr_ccliq(pcols,pver)                     ! mixing_ratio_convective_cloud_liquid (kg/kg)
   real(r8) :: mr_ccice(pcols,pver)                     ! mixing_ratio_convective_cloud_ice (kg/kg)

   real(r8) :: pracwo(pcols,pver)                       ! RK accretion of cloud water by rain (1/s)
   real(r8) :: psacwo(pcols,pver)                       ! RK accretion of cloud water by snow (1/s)
   real(r8) :: psacio(pcols,pver)                       ! RK accretion of cloud ice by snow (1/s)

   real(r8) :: iwc(pcols,pver)                         ! Grid box average ice water content
   real(r8) :: lwc(pcols,pver)                         ! Grid box average liquid water content  
   
   logical  :: lq(pcnst)
   logical  :: use_shfrc_local
   logical  :: cam3_local
   ! ======================================================================

   lchnk = state%lchnk
   ncol  = state%ncol

   call stratiform_select_impl()
   if (.not. use_native_impl) then
      call stratiform_log_entered_once(stratiform_entered_logged, 'STRATIFORM_PROOF_FILE', &
           'stratiform_tend entered (active helper shells direct = codon; pcond/cldfrc = native)')
      call stratiform_select_branches(use_shfrc, cam_physpkg_is('cam3'))
      use_shfrc_local = iand(branch_mask, 1) /= 0
      cam3_local = iand(branch_mask, 2) /= 0
   else
      use_shfrc_local = use_shfrc
      cam3_local = cam_physpkg_is('cam3')
   end if

   call physics_state_copy(state,state1)             ! Copy state to local state1.

   ! Associate pointers with physics buffer fields

   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, qcwat_idx,   qcwat,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
   call pbuf_get_field(pbuf, tcwat_idx,   tcwat,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
   call pbuf_get_field(pbuf, lcwat_idx,   lcwat,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

   call pbuf_get_field(pbuf, cld_idx,     cld,     start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
   call pbuf_get_field(pbuf, concld_idx,  concld,  start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
   call pbuf_get_field(pbuf, ast_idx,     ast,     start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

   call pbuf_get_field(pbuf, fice_idx,    fice)
 
   call pbuf_get_field(pbuf, prec_str_idx, prec_str)
   call pbuf_get_field(pbuf, snow_str_idx, snow_str)
   call pbuf_get_field(pbuf, prec_sed_idx, prec_sed)
   call pbuf_get_field(pbuf, snow_sed_idx, snow_sed)
   call pbuf_get_field(pbuf, prec_pcw_idx, prec_pcw)
   call pbuf_get_field(pbuf, snow_pcw_idx, snow_pcw)

   call pbuf_get_field(pbuf, qme_idx,    qme )
   call pbuf_get_field(pbuf, prain_idx,  prain)
   call pbuf_get_field(pbuf, nevapr_idx, nevapr)

   call pbuf_get_field(pbuf, rel_idx,    rel)
   call pbuf_get_field(pbuf, rei_idx,    rei)

   call pbuf_get_field(pbuf, wsedl_idx,  wsedl)

   ! ------------- !
   ! Sedimentation !
   ! ------------- !

   ! Allow the cloud liquid drops and ice particles to sediment.
   ! This is done before adding convectively detrained cloud water, 
   ! because the phase of the detrained water is unknown.

   call t_startf('stratiform_sediment')

   call cld_sediment_vel( ncol,                                                           &
                          icefrac, landfrac, ocnfrac, state1%pmid, state1%pdel, state1%t, &
                          cld, state1%q(:,:,ixcldliq), state1%q(:,:,ixcldice),            & 
                          pvliq, pvice, landm, snowh )

   if (use_native_impl) then
      wsedl(:ncol,:pver) = pvliq(:ncol,:pver)/gravit/(state1%pmid(:ncol,:pver)/(287.15_r8*state1%t(:ncol,:pver)))
   end if

   lq(:)        = .FALSE.
   lq(1)        = .TRUE.
   lq(ixcldice) = .TRUE.
   lq(ixcldliq) = .TRUE.
   call physics_ptend_init(ptend_loc, state%psetcols, 'pcwsediment', ls=.true., lq=lq)! Initialize local ptend type

   call cld_sediment_tend( ncol, dtime ,                                                             &
                           state1%pint, state1%pmid, state1%pdel, state1%t,                          &
                           cld, state1%q(:,:,ixcldliq), state1%q(:,:,ixcldice), pvliq, pvice,        &
                           ptend_loc%q(:,:,ixcldliq), ptend_loc%q(:,:,ixcldice), ptend_loc%q(:,:,1), &
                           ptend_loc%s, rain, snow_sed )

   ! Convert rain and snow fluxes at the surface from [kg/m2/s] to [m/s]
   ! Compute total precipitation flux at the surface in [m/s]

   if (use_native_impl) then
      snow_sed(:ncol) = snow_sed(:ncol)/1000._r8
      rain(:ncol)     = rain(:ncol)/1000._r8
      prec_sed(:ncol) = rain(:ncol) + snow_sed(:ncol)
   else
      call stratiform_sedimentation_diag(ncol, state1%pmid, state1%t, pvliq, wsedl, rain, snow_sed, prec_sed)
   end if

   ! Record history variables
   call outfld( 'DQSED'   ,ptend_loc%q(:,:,1)       , pcols,lchnk )
   call outfld( 'DISED'   ,ptend_loc%q(:,:,ixcldice), pcols,lchnk )
   call outfld( 'DLSED'   ,ptend_loc%q(:,:,ixcldliq), pcols,lchnk )
   call outfld( 'HSED'    ,ptend_loc%s              , pcols,lchnk )
   call outfld( 'PRECSED' ,prec_sed                 , pcols,lchnk )
   call outfld( 'SNOWSED' ,snow_sed                 , pcols,lchnk )
   call outfld( 'RAINSED' ,rain                     , pcols,lchnk )

   ! Add tendency from this process to tend from other processes here
   call physics_ptend_init(ptend_all, state%psetcols, 'stratiform')
   call physics_ptend_sum( ptend_loc, ptend_all, ncol )

   ! Update physics state type state1 with ptend_loc 
   call physics_update( state1, ptend_loc, dtime )

   call t_stopf('stratiform_sediment')

   ! Accumulate prec and snow flux at the surface [ m/s ]
   prec_str(:ncol) = prec_sed(:ncol)
   snow_str(:ncol) = snow_sed(:ncol)

   ! ----------------------------------------------------------------------------- !
   ! Detrainment of convective condensate into the environment or stratiform cloud !
   ! ----------------------------------------------------------------------------- !

   ! Put all of the detraining cloud water from convection into the large scale cloud.
   ! It all goes in liquid for the moment.
   ! Strictly speaking, this approach is detraining all the cconvective water into 
   ! the environment, not the large-scale cloud.

   lq(:)        = .FALSE.
   lq(ixcldliq) = .TRUE.
   call physics_ptend_init( ptend_loc, state1%psetcols, 'pcwdetrain', lq=lq)
   
   call stratiform_detrain_assign(ncol, dlf, ptend_loc%q)

   call outfld( 'ZMDLF', dlf, pcols, lchnk )
   call outfld( 'SHDLF', dlf2, pcols, lchnk )

   ! Add hie detrainment tendency to tend from the other prior processes

   call physics_ptend_sum( ptend_loc, ptend_all, ncol )
   call physics_update( state1, ptend_loc, dtime )

   ! Accumulate prec and snow, reserved liquid has now been used.

   prec_str(:ncol) = prec_str(:ncol) - rliq(:ncol)  ! ( snow contribution is zero )

   ! -------------------------------------- !
   ! Computation of Various Cloud Fractions !
   ! -------------------------------------- !

   ! ----------------------------------------------------------------------------- !
   ! Treatment of cloud fraction in CAM4 and CAM5 differs                          !  
   ! (1) CAM4                                                                      !
   !     . Cumulus AMT = Deep    Cumulus AMT ( empirical fcn of mass flux ) +      !
   !                     Shallow Cumulus AMT ( empirical fcn of mass flux )        !
   !     . Stratus AMT = max( RH stratus AMT, Stability Stratus AMT )              !
   !     . Cumulus and Stratus are 'minimally' overlapped without hierarchy.       !
   !     . Cumulus LWC,IWC is assumed to be the same as Stratus LWC,IWC            !
   ! (2) CAM5                                                                      !
   !     . Cumulus AMT = Deep    Cumulus AMT ( empirical fcn of mass flux ) +      !
   !                     Shallow Cumulus AMT ( internally fcn of mass flux and w ) !
   !     . Stratus AMT = fcn of environmental-mean RH ( no Stability Stratus )     !
   !     . Cumulus and Stratus are non-overlapped with higher priority on Cumulus  !
   !     . Cumulus ( both Deep and Shallow ) has its own LWC and IWC.              !
   ! ----------------------------------------------------------------------------- ! 

   if( use_shfrc_local ) then
       call pbuf_get_field(pbuf, shfrc_idx, shfrc )
   else
       shfrc=>shfrc_local
       shfrc(:,:) = 0._r8
   endif

   ! Stratus ('ast' = max(alst,aist)) and total cloud fraction ('cld = ast + concld')
   ! will be computed using this updated 'concld' in the stratiform macrophysics 
   ! scheme (mmacro_pcond) later below. 

   call t_startf("cldfrc")
   call cldfrc( lchnk, ncol, pbuf,                                  &
                state1%pmid, state1%t, state1%q(:,:,1), state1%omega, state1%phis, &
                shfrc, use_shfrc_local,                                            &
                cld, rhcloud, clc, state1%pdel,                                    &
                cmfmc, cmfmc2, landfrac,snowh, concld, cldst,                      &
                ts, sst, state1%pint(:,pverp), zdu, ocnfrac, rhu00,                &
                state1%q(:,:,ixcldice), icecldf, liqcldf,                          &
                relhum, 0 )    

   ! Re-calculate cloud with perturbed rh add call cldfrc to estimate rhdfda.

   call cldfrc( lchnk, ncol, pbuf,                                  &
                state1%pmid, state1%t, state1%q(:,:,1), state1%omega, state1%phis, &
                shfrc, use_shfrc_local,                                            &
                cld2, rhcloud2, clc, state1%pdel,                                  &
                cmfmc, cmfmc2, landfrac, snowh, concld2, cldst2,                   &
                ts, sst, state1%pint(:,pverp), zdu, ocnfrac, rhu002,               &
                state1%q(:,:,ixcldice), icecldf2, liqcldf2,                        &
                relhum2, 1 )              

   call t_stopf("cldfrc")

   call stratiform_rhdfda(ncol, relhum, rhu00, cld, cld2, rhdfda)

   ! ---------------------------------------------- !
   ! Stratiform Cloud Macrophysics and Microphysics !
   ! ---------------------------------------------- !

   call t_startf('stratiform_microphys')

   rdtime = 1._r8/dtime

   ! Define fractional amount of stratus condensate and precipitation in ice phase.
   ! This uses a ramp ( -30 ~ -10 for fice, -5 ~ 0 for fsnow ). 
   ! The ramp within convective cloud may be different

   call cldfrc_fice(ncol, state1%t, fice, fsnow)

   ! Perform repartitioning of stratiform condensate.    
   ! Corresponding heating tendency will be added later. 

   lq(:)        = .FALSE.
   lq(ixcldice) = .true.
   lq(ixcldliq) = .true.
   call physics_ptend_init( ptend_loc, state1%psetcols, 'cldwat-repartition', lq=lq )

   call stratiform_repartition(ncol, rdtime, fice, state1%q, totcw, repartht, ptend_loc%q)

   call outfld( 'REPARTICE', ptend_loc%q(:,:,ixcldice), pcols, lchnk )
   call outfld( 'REPARTLIQ', ptend_loc%q(:,:,ixcldliq), pcols, lchnk )

   call physics_ptend_sum( ptend_loc, ptend_all, ncol )
   call physics_update( state1, ptend_loc, dtime )

   ! Determine repartition heating from change in cloud ice.

   ! Non-micro and non-macrophysical external advective forcings to compute net condensation rate. 
   ! Note that advective forcing of condensate is aggregated into liquid phase.
   call stratiform_forcing_prep(ncol, dtime, state1%q, state1%t, qcwat, tcwat, lcwat, totcw, repartht, qtend, ttend, ltend)

   ! Compute Stratiform Macro-Microphysical Tendencies

   ! Add rain and snow fluxes as output variables from pcond, and into physics buffer
   call pbuf_get_field(pbuf, ls_flxprc_idx, rkflxprc)
   call pbuf_get_field(pbuf, ls_flxsnw_idx, rkflxsnw)

   call t_startf('pcond')
   call pcond( lchnk, ncol,                                                &
               state1%t, ttend, state1%q(1,1,1), qtend, state1%omega,      &
               totcw, state1%pmid , state1%pdel, cld, fice, fsnow,         &
               qme, prain, prodsnow, nevapr, evapsnow, evapheat, prfzheat, &
               meltheat, prec_pcw, snow_pcw, dtime, fwaut,                 &
               fsaut, fracw, fsacw, fsaci, ltend,                          &
               rhdfda, rhu00, icefrac, state1%zi, ice2pr, liq2pr,          &
               liq2snow, snowh, rkflxprc, rkflxsnw, pracwo, psacwo, psacio )
   call t_stopf('pcond')

   lq(:)        = .FALSE.
   lq(1)        = .true.
   lq(ixcldice) = .true.
   lq(ixcldliq) = .true.
   call physics_ptend_init( ptend_loc, state1%psetcols, 'cldwat', ls=.true., lq=lq)

   call stratiform_select_microphys_shell_impl()

   if (use_native_microphys_shell_impl) then
      call stratiform_microphys_tend_diag(ncol, dtime, cld, fice, qme, nevapr, evapheat, prfzheat, meltheat, ice2pr, liq2pr, &
           state1%q, prec_pcw, snow_pcw, repartht, ptend_loc%q, ptend_loc%s, ast, icimr, icwmr, cmeheat, cmeice, cmeliq)

      call stratiform_cloud_mixing_diag(ncol, state%q, cld, concld, mr_ccliq, mr_ccice, mr_lsliq, mr_lsice)
   else
      call stratiform_microphys_shell_codon_wrap(ncol, dtime, cld, concld, fice, qme, nevapr, evapheat, prfzheat, meltheat, &
           ice2pr, liq2pr, state1%q, prec_pcw, snow_pcw, repartht, ptend_loc%q, ptend_loc%s, ast, icimr, icwmr, cmeheat, &
           cmeice, cmeliq, mr_ccliq, mr_ccice, mr_lsliq, mr_lsice)
   end if

   ! Record history variables

   call outfld( 'FWAUT'   , fwaut,       pcols, lchnk )
   call outfld( 'FSAUT'   , fsaut,       pcols, lchnk )
   call outfld( 'FRACW'   , fracw,       pcols, lchnk )
   call outfld( 'FSACW'   , fsacw,       pcols, lchnk )
   call outfld( 'FSACI'   , fsaci,       pcols, lchnk )

   call outfld( 'PCSNOW'  , snow_pcw,    pcols, lchnk )
   call outfld( 'FICE'    , fice,        pcols, lchnk )
   call outfld( 'CMEICE'  , cmeice,      pcols, lchnk )
   call outfld( 'CMELIQ'  , cmeliq,      pcols, lchnk )
   call outfld( 'ICE2PR'  , ice2pr,      pcols, lchnk )
   call outfld( 'LIQ2PR'  , liq2pr,      pcols, lchnk )
   call outfld( 'HPROGCLD', ptend_loc%s, pcols, lchnk )
   call outfld( 'HEVAP   ', evapheat,    pcols, lchnk )
   call outfld( 'HMELT'   , meltheat,    pcols, lchnk )
   call outfld( 'HCME'    , cmeheat ,    pcols, lchnk )
   call outfld( 'HFREEZ'  , prfzheat,    pcols, lchnk )
   call outfld( 'HREPART' , repartht,    pcols, lchnk )
   call outfld('LS_FLXPRC', rkflxprc,    pcols, lchnk )
   call outfld('LS_FLXSNW', rkflxsnw,    pcols, lchnk )
   call outfld('PRACWO'   , pracwo,      pcols, lchnk )
   call outfld('PSACWO'   , psacwo,      pcols, lchnk )
   call outfld('PSACIO'   , psacio,      pcols, lchnk )

   call outfld( 'CLDLIQSTR  ', mr_lsliq,    pcols, lchnk )
   call outfld( 'CLDICESTR  ', mr_lsice,    pcols, lchnk )
   call outfld( 'CLDLIQCON  ', mr_ccliq,    pcols, lchnk )
   call outfld( 'CLDICECON  ', mr_ccice,    pcols, lchnk )

   ! ------------------------------- !
   ! Update microphysical tendencies !
   ! ------------------------------- !

   call physics_ptend_sum( ptend_loc, ptend_all, ncol )
   call physics_update( state1, ptend_loc, dtime )

   if (.not. cam3_local) then

      call t_startf("cldfrc")
      call cldfrc( lchnk, ncol, pbuf,                                  &
                   state1%pmid, state1%t, state1%q(:,:,1), state1%omega, state1%phis, &
                   shfrc, use_shfrc_local,                                            &
                   cld, rhcloud, clc, state1%pdel,                                    &
                   cmfmc, cmfmc2, landfrac, snowh, concld, cldst,                     &
                   ts, sst, state1%pint(:,pverp), zdu, ocnfrac, rhu00,                &
                   state1%q(:,:,ixcldice), icecldf, liqcldf,                          &
                   relhum, 0 )    
      call t_stopf("cldfrc")

   endif

   call outfld( 'CONCLD  ', concld, pcols, lchnk )
   call outfld( 'CLDST   ', cldst,  pcols, lchnk )
   call outfld( 'CNVCLD  ', clc,    pcols, lchnk )
   call outfld( 'AST',      ast,    pcols, lchnk )   

   call stratiform_postcloud_diag(ncol, state1%q, state1%pmid, state1%t, rhcloud, iwc, lwc, icimr, icwmr)

   call outfld( 'IWC'      , iwc,         pcols, lchnk )
   call outfld( 'LWC'      , lwc,         pcols, lchnk )
   call outfld( 'ICIMR'    , icimr,       pcols, lchnk )
   call outfld( 'ICWMR'    , icwmr,       pcols, lchnk )
   call outfld( 'CME'      , qme,         pcols, lchnk )
   call outfld( 'PRODPREC' , prain,       pcols, lchnk )
   call outfld( 'EVAPPREC' , nevapr,      pcols, lchnk )
   call outfld( 'EVAPSNOW' , evapsnow,    pcols, lchnk )

   call t_stopf('stratiform_microphys')

   prec_str(:ncol) = prec_str(:ncol) + prec_pcw(:ncol)
   snow_str(:ncol) = snow_str(:ncol) + snow_pcw(:ncol)

   ! Save variables for use in the macrophysics at the next time step

   call stratiform_store_oldcloud(ncol, state1%q, state1%t, qcwat, tcwat, lcwat)
  
   ! Cloud water and ice particle sizes, saved in physics buffer for radiation

   call cldefr( lchnk, ncol, landfrac, state1%t, rel, rei, state1%ps, state1%pmid, landm, icefrac, snowh )

   call outfld('REL', rel, pcols, lchnk)
   call outfld('REI', rei, pcols, lchnk)

   call physics_state_dealloc(state1)

end subroutine stratiform_tend

!===============================================================================

subroutine stratiform_detrain_assign(ncol_local, dlf_local, ptend_q_local)

   use iso_c_binding, only: c_int64_t, c_loc, c_ptr

   implicit none

   integer, intent(in) :: ncol_local
   real(r8), target, intent(in) :: dlf_local(pcols,pver)
   real(r8), target, intent(inout) :: ptend_q_local(pcols,pver,pcnst)

   interface
      subroutine stratiform_detrain_assign_codon(ncol_c, pcols_c, pver_c, pcnst_c, ixcldliq_c, dlf_p, ptend_q_p) &
           bind(c, name="stratiform_detrain_assign_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pcnst_c, ixcldliq_c
         type(c_ptr), value :: dlf_p, ptend_q_p
      end subroutine stratiform_detrain_assign_codon
   end interface

   if (use_native_impl) then
      call stratiform_detrain_assign_native(ncol_local, dlf_local, ptend_q_local)
      return
   end if

   call stratiform_log_entered_once(detrain_assign_entered_logged, 'STRATIFORM_PROOF_FILE', &
        'stratiform_detrain_assign entered (detrain cloud liquid assign direct = codon)')

   call stratiform_detrain_assign_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
        int(pcnst, c_int64_t), int(ixcldliq, c_int64_t), c_loc(dlf_local), c_loc(ptend_q_local))

end subroutine stratiform_detrain_assign

!===============================================================================

subroutine stratiform_detrain_assign_native(ncol_local, dlf_local, ptend_q_local)

   implicit none

   integer, intent(in) :: ncol_local
   real(r8), intent(in) :: dlf_local(pcols,pver)
   real(r8), intent(inout) :: ptend_q_local(pcols,pver,pcnst)

   integer :: i, k

   do k = 1, pver
      do i = 1, ncol_local
         ptend_q_local(i,k,ixcldliq) = dlf_local(i,k)
      end do
   end do

end subroutine stratiform_detrain_assign_native

!===============================================================================

subroutine stratiform_append_impl_proof(env_name, proof_line)

   character(len=*), intent(in) :: env_name, proof_line
   character(len=512) :: proof_path
   integer :: status, n, unit_id

   if (.not. masterproc) return

   proof_path = ''
   call get_environment_variable(env_name, value=proof_path, length=n, status=status)
   if (status /= 0 .or. n <= 0) return

   open(newunit=unit_id, file=trim(adjustl(proof_path(:n))), status='unknown', action='write', &
        position='append', iostat=status)
   if (status /= 0) return

   write(unit_id,'(A)') trim(proof_line)
   close(unit_id)

end subroutine stratiform_append_impl_proof

!===============================================================================

subroutine stratiform_log_entered_once(entered_logged, env_name, proof_line)

   logical, intent(inout) :: entered_logged
   character(len=*), intent(in) :: env_name, proof_line

   if (entered_logged) return
   entered_logged = .true.

   if (masterproc) then
      write(iulog,*) trim(proof_line)
      call stratiform_append_impl_proof(env_name, proof_line)
      call flush(iulog)
   end if

end subroutine stratiform_log_entered_once

!===============================================================================

subroutine stratiform_select_microphys_shell_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (microphys_shell_impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('STRATIFORM_MICROPHYS_SHELL_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_microphys_shell_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_microphys_shell_impl = .false.
   end if

   microphys_shell_impl_selected = .true.

   if (masterproc) then
      if (use_native_microphys_shell_impl) then
         write(iulog,*) 'stratiform_microphys_shell implementation = native'
         call stratiform_append_impl_proof('STRATIFORM_MICROPHYS_SHELL_PROOF_FILE', &
              'stratiform_microphys_shell implementation = native')
      else
         write(iulog,*) 'stratiform_microphys_shell implementation = codon'
         call stratiform_append_impl_proof('STRATIFORM_MICROPHYS_SHELL_PROOF_FILE', &
              'stratiform_microphys_shell implementation = codon')
      end if
      call flush(iulog)
   end if

end subroutine stratiform_select_microphys_shell_impl

!===============================================================================

subroutine stratiform_microphys_shell_codon_wrap(ncol_local, dtime_local, cld_local, concld_local, fice_local, qme_local, &
     nevapr_local, evapheat_local, prfzheat_local, meltheat_local, ice2pr_local, liq2pr_local, state_q_local, prec_pcw_local, &
     snow_pcw_local, repartht_local, ptend_q_local, ptend_s_local, ast_local, icimr_local, icwmr_local, cmeheat_local, &
     cmeice_local, cmeliq_local, mr_ccliq_local, mr_ccice_local, mr_lsliq_local, mr_lsice_local)

   use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

   integer, intent(in) :: ncol_local
   real(r8), intent(in) :: dtime_local
   real(r8), target, intent(in) :: cld_local(pcols,pver), concld_local(pcols,pver), fice_local(pcols,pver)
   real(r8), target, intent(in) :: qme_local(pcols,pver), nevapr_local(pcols,pver), evapheat_local(pcols,pver)
   real(r8), target, intent(in) :: prfzheat_local(pcols,pver), meltheat_local(pcols,pver), ice2pr_local(pcols,pver)
   real(r8), target, intent(in) :: liq2pr_local(pcols,pver), state_q_local(pcols,pver,pcnst), repartht_local(pcols,pver)
   real(r8), target, intent(inout) :: prec_pcw_local(pcols), snow_pcw_local(pcols), ptend_q_local(pcols,pver,pcnst)
   real(r8), target, intent(inout) :: ptend_s_local(pcols,pver), ast_local(pcols,pver), icimr_local(pcols,pver)
   real(r8), target, intent(inout) :: icwmr_local(pcols,pver), cmeheat_local(pcols,pver), cmeice_local(pcols,pver)
   real(r8), target, intent(inout) :: cmeliq_local(pcols,pver), mr_ccliq_local(pcols,pver), mr_ccice_local(pcols,pver)
   real(r8), target, intent(inout) :: mr_lsliq_local(pcols,pver), mr_lsice_local(pcols,pver)

   interface
      subroutine stratiform_microphys_shell_codon(ncol_c, pcols_c, pver_c, pcnst_c, ixcldliq_c, ixcldice_c, dtime_c, latvap_c, &
           latice_c, cld_p, concld_p, fice_p, qme_p, nevapr_p, evapheat_p, prfzheat_p, meltheat_p, ice2pr_p, liq2pr_p, &
           state_q_p, prec_pcw_p, snow_pcw_p, repartht_p, ptend_q_p, ptend_s_p, ast_p, icimr_p, icwmr_p, cmeheat_p, &
           cmeice_p, cmeliq_p, mr_ccliq_p, mr_ccice_p, mr_lsliq_p, mr_lsice_p) bind(c, name="stratiform_microphys_shell_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pcnst_c, ixcldliq_c, ixcldice_c
         real(c_double), value :: dtime_c, latvap_c, latice_c
         type(c_ptr), value :: cld_p, concld_p, fice_p, qme_p, nevapr_p, evapheat_p, prfzheat_p, meltheat_p
         type(c_ptr), value :: ice2pr_p, liq2pr_p, state_q_p, prec_pcw_p, snow_pcw_p, repartht_p, ptend_q_p, ptend_s_p
         type(c_ptr), value :: ast_p, icimr_p, icwmr_p, cmeheat_p, cmeice_p, cmeliq_p, mr_ccliq_p, mr_ccice_p
         type(c_ptr), value :: mr_lsliq_p, mr_lsice_p
      end subroutine stratiform_microphys_shell_codon
   end interface

   call stratiform_log_entered_once(microphys_shell_entered_logged, 'STRATIFORM_MICROPHYS_SHELL_PROOF_FILE', &
        'stratiform_microphys_shell entered (microphysics tendency/mixing diagnostics direct = codon)')

   call stratiform_microphys_shell_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
        int(pcnst, c_int64_t), int(ixcldliq, c_int64_t), int(ixcldice, c_int64_t), dtime_local, latvap, latice, &
        c_loc(cld_local), c_loc(concld_local), c_loc(fice_local), c_loc(qme_local), c_loc(nevapr_local), c_loc(evapheat_local), &
        c_loc(prfzheat_local), c_loc(meltheat_local), c_loc(ice2pr_local), c_loc(liq2pr_local), c_loc(state_q_local), &
        c_loc(prec_pcw_local), c_loc(snow_pcw_local), c_loc(repartht_local), c_loc(ptend_q_local), c_loc(ptend_s_local), &
        c_loc(ast_local), c_loc(icimr_local), c_loc(icwmr_local), c_loc(cmeheat_local), c_loc(cmeice_local), c_loc(cmeliq_local), &
        c_loc(mr_ccliq_local), c_loc(mr_ccice_local), c_loc(mr_lsliq_local), c_loc(mr_lsice_local))

end subroutine stratiform_microphys_shell_codon_wrap

!===============================================================================

subroutine stratiform_sedimentation_diag(ncol_local, state_pmid_local, state_t_local, pvliq_local, wsedl_local, rain_local, &
     snow_sed_local, prec_sed_local)

   use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

   implicit none

   integer, intent(in) :: ncol_local
   real(r8), target, intent(in) :: state_pmid_local(pcols,pver), state_t_local(pcols,pver), pvliq_local(pcols,pverp)
   real(r8), target, intent(inout) :: wsedl_local(pcols,pver), prec_sed_local(pcols)
   real(r8), target, intent(inout) :: rain_local(pcols), snow_sed_local(pcols)

   interface
      subroutine stratiform_sedimentation_diag_codon(ncol_c, pcols_c, pver_c, gravit_c, state_pmid_p, state_t_p, pvliq_p, &
           wsedl_p, rain_p, snow_sed_p, prec_sed_p) bind(c, name="stratiform_sedimentation_diag_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         real(c_double), value :: gravit_c
         type(c_ptr), value :: state_pmid_p, state_t_p, pvliq_p, wsedl_p, rain_p, snow_sed_p, prec_sed_p
      end subroutine stratiform_sedimentation_diag_codon
   end interface

   if (use_native_impl) then
      call stratiform_sedimentation_diag_native(ncol_local, state_pmid_local, state_t_local, pvliq_local, wsedl_local, &
           rain_local, snow_sed_local, prec_sed_local)
      return
   end if

   call stratiform_log_entered_once(sedimentation_diag_entered_logged, 'STRATIFORM_PROOF_FILE', &
        'stratiform_sedimentation_diag entered (sedimentation diagnostics direct = codon)')

   call stratiform_sedimentation_diag_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
        gravit, c_loc(state_pmid_local), c_loc(state_t_local), c_loc(pvliq_local), c_loc(wsedl_local), c_loc(rain_local), &
        c_loc(snow_sed_local), c_loc(prec_sed_local))

end subroutine stratiform_sedimentation_diag

!===============================================================================

subroutine stratiform_sedimentation_diag_native(ncol_local, state_pmid_local, state_t_local, pvliq_local, wsedl_local, &
     rain_local, snow_sed_local, prec_sed_local)

   implicit none

   integer, intent(in) :: ncol_local
   real(r8), intent(in) :: state_pmid_local(pcols,pver), state_t_local(pcols,pver), pvliq_local(pcols,pverp)
   real(r8), intent(inout) :: wsedl_local(pcols,pver), prec_sed_local(pcols)
   real(r8), intent(inout) :: rain_local(pcols), snow_sed_local(pcols)

   integer :: i, k

   do k = 1, pver
      do i = 1, ncol_local
         wsedl_local(i,k) = pvliq_local(i,k)/gravit/(state_pmid_local(i,k)/(287.15_r8*state_t_local(i,k)))
      end do
   end do

   snow_sed_local(:ncol_local) = snow_sed_local(:ncol_local)/1000._r8
   rain_local(:ncol_local) = rain_local(:ncol_local)/1000._r8
   prec_sed_local(:ncol_local) = rain_local(:ncol_local) + snow_sed_local(:ncol_local)

end subroutine stratiform_sedimentation_diag_native

!===============================================================================

subroutine stratiform_rhdfda(ncol_local, relhum_local, rhu00_local, cld_local, cld2_local, rhdfda_local)

   use iso_c_binding, only: c_int64_t, c_loc, c_ptr

   implicit none

   integer, intent(in) :: ncol_local
   real(r8), target, intent(in) :: relhum_local(pcols,pver), cld_local(pcols,pver), cld2_local(pcols,pver)
   real(r8), target, intent(inout) :: rhu00_local(pcols,pver)
   real(r8), target, intent(inout) :: rhdfda_local(pcols,pver)

   interface
      subroutine stratiform_rhdfda_codon(ncol_c, pcols_c, pver_c, relhum_p, rhu00_p, cld_p, cld2_p, rhdfda_p) &
           bind(c, name="stratiform_rhdfda_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         type(c_ptr), value :: relhum_p, rhu00_p, cld_p, cld2_p, rhdfda_p
      end subroutine stratiform_rhdfda_codon
   end interface

   if (use_native_impl) then
      call stratiform_rhdfda_native(ncol_local, relhum_local, rhu00_local, cld_local, cld2_local, rhdfda_local)
      return
   end if

   call stratiform_log_entered_once(rhdfda_entered_logged, 'STRATIFORM_PROOF_FILE', &
        'stratiform_rhdfda entered (RH perturbation diagnostic direct = codon)')

   call stratiform_rhdfda_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), c_loc(relhum_local), &
        c_loc(rhu00_local), c_loc(cld_local), c_loc(cld2_local), c_loc(rhdfda_local))

end subroutine stratiform_rhdfda

!===============================================================================

subroutine stratiform_rhdfda_native(ncol_local, relhum_local, rhu00_local, cld_local, cld2_local, rhdfda_local)

   implicit none

   integer, intent(in) :: ncol_local
   real(r8), intent(in) :: relhum_local(pcols,pver), cld_local(pcols,pver), cld2_local(pcols,pver)
   real(r8), intent(inout) :: rhu00_local(pcols,pver)
   real(r8), intent(inout) :: rhdfda_local(pcols,pver)

   integer :: i, k

   rhu00_local(:ncol_local,1) = 2.0_r8
   do k = 1, pver
      do i = 1, ncol_local
         if( relhum_local(i,k) < rhu00_local(i,k) ) then
            rhdfda_local(i,k) = 0.0_r8
         elseif( relhum_local(i,k) >= 1.0_r8 ) then
            rhdfda_local(i,k) = 0.0_r8
         else
            if( ( cld2_local(i,k) - cld_local(i,k) ) < 1.e-4_r8 ) then
               rhdfda_local(i,k) = 0.01_r8*relhum_local(i,k)*1.e+4_r8
            else
               rhdfda_local(i,k) = 0.01_r8*relhum_local(i,k)/(cld2_local(i,k)-cld_local(i,k))
            endif
         endif
      enddo
   enddo

end subroutine stratiform_rhdfda_native

!===============================================================================

subroutine stratiform_repartition(ncol_local, rdtime_local, fice_local, state_q_local, totcw_local, repartht_local, ptend_q_local)

   use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

   implicit none

   integer, intent(in) :: ncol_local
   real(r8), intent(in) :: rdtime_local
   real(r8), target, intent(in) :: fice_local(pcols,pver), state_q_local(pcols,pver,pcnst)
   real(r8), target, intent(inout) :: totcw_local(pcols,pver), repartht_local(pcols,pver)
   real(r8), target, intent(inout) :: ptend_q_local(pcols,pver,pcnst)

   interface
      subroutine stratiform_repartition_codon(ncol_c, pcols_c, pver_c, pcnst_c, ixcldliq_c, ixcldice_c, rdtime_c, fice_p, &
           state_q_p, totcw_p, repartht_p, ptend_q_p) bind(c, name="stratiform_repartition_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pcnst_c, ixcldliq_c, ixcldice_c
         real(c_double), value :: rdtime_c
         type(c_ptr), value :: fice_p, state_q_p, totcw_p, repartht_p, ptend_q_p
      end subroutine stratiform_repartition_codon
   end interface

   if (use_native_impl) then
      call stratiform_repartition_native(ncol_local, rdtime_local, fice_local, state_q_local, totcw_local, repartht_local, &
           ptend_q_local)
      return
   end if

   call stratiform_log_entered_once(repartition_entered_logged, 'STRATIFORM_PROOF_FILE', &
        'stratiform_repartition entered (cloud condensate repartition direct = codon)')

   call stratiform_repartition_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
        int(pcnst, c_int64_t), int(ixcldliq, c_int64_t), int(ixcldice, c_int64_t), rdtime_local, c_loc(fice_local), &
        c_loc(state_q_local), c_loc(totcw_local), c_loc(repartht_local), c_loc(ptend_q_local))

end subroutine stratiform_repartition

!===============================================================================

subroutine stratiform_repartition_native(ncol_local, rdtime_local, fice_local, state_q_local, totcw_local, repartht_local, &
     ptend_q_local)

   implicit none

   integer, intent(in) :: ncol_local
   real(r8), intent(in) :: rdtime_local
   real(r8), intent(in) :: fice_local(pcols,pver), state_q_local(pcols,pver,pcnst)
   real(r8), intent(inout) :: totcw_local(pcols,pver), repartht_local(pcols,pver)
   real(r8), intent(inout) :: ptend_q_local(pcols,pver,pcnst)

   totcw_local(:ncol_local,:pver) = state_q_local(:ncol_local,:pver,ixcldice) + state_q_local(:ncol_local,:pver,ixcldliq)
   repartht_local(:ncol_local,:pver) = state_q_local(:ncol_local,:pver,ixcldice)
   ptend_q_local(:ncol_local,:pver,ixcldice) = rdtime_local * (totcw_local(:ncol_local,:pver)*fice_local(:ncol_local,:pver) - &
        state_q_local(:ncol_local,:pver,ixcldice))
   ptend_q_local(:ncol_local,:pver,ixcldliq) = rdtime_local * ( &
        totcw_local(:ncol_local,:pver)*(1.0_r8-fice_local(:ncol_local,:pver)) - &
        state_q_local(:ncol_local,:pver,ixcldliq))

end subroutine stratiform_repartition_native

!===============================================================================

subroutine stratiform_forcing_prep(ncol_local, dtime_local, state_q_local, state_t_local, qcwat_local, tcwat_local, lcwat_local, &
     totcw_local, repartht_local, qtend_local, ttend_local, ltend_local)

   use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

   implicit none

   integer, intent(in) :: ncol_local
   real(r8), intent(in) :: dtime_local
   real(r8), target, intent(in) :: state_q_local(pcols,pver,pcnst), state_t_local(pcols,pver), qcwat_local(pcols,pver), &
        tcwat_local(pcols,pver), lcwat_local(pcols,pver), totcw_local(pcols,pver)
   real(r8), target, intent(inout) :: repartht_local(pcols,pver)
   real(r8), target, intent(inout) :: qtend_local(pcols,pver), ttend_local(pcols,pver), ltend_local(pcols,pver)

   interface
      subroutine stratiform_forcing_prep_codon(ncol_c, pcols_c, pver_c, pcnst_c, ixcldliq_c, ixcldice_c, dtime_c, latice_c, &
           state_q_p, state_t_p, qcwat_p, tcwat_p, lcwat_p, totcw_p, repartht_p, qtend_p, ttend_p, ltend_p) &
           bind(c, name="stratiform_forcing_prep_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pcnst_c, ixcldliq_c, ixcldice_c
         real(c_double), value :: dtime_c, latice_c
         type(c_ptr), value :: state_q_p, state_t_p, qcwat_p, tcwat_p, lcwat_p, totcw_p, repartht_p, qtend_p, ttend_p, ltend_p
      end subroutine stratiform_forcing_prep_codon
   end interface

   if (use_native_impl) then
      call stratiform_forcing_prep_native(ncol_local, dtime_local, state_q_local, state_t_local, qcwat_local, tcwat_local, &
           lcwat_local, totcw_local, repartht_local, qtend_local, ttend_local, ltend_local)
      return
   end if

   call stratiform_log_entered_once(forcing_prep_entered_logged, 'STRATIFORM_PROOF_FILE', &
        'stratiform_forcing_prep entered (pcond forcing prep direct = codon)')

   call stratiform_forcing_prep_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
        int(pcnst, c_int64_t), int(ixcldliq, c_int64_t), int(ixcldice, c_int64_t), dtime_local, latice, &
        c_loc(state_q_local), c_loc(state_t_local), c_loc(qcwat_local), c_loc(tcwat_local), c_loc(lcwat_local), &
        c_loc(totcw_local), c_loc(repartht_local), c_loc(qtend_local), c_loc(ttend_local), c_loc(ltend_local))

end subroutine stratiform_forcing_prep

!===============================================================================

subroutine stratiform_forcing_prep_native(ncol_local, dtime_local, state_q_local, state_t_local, qcwat_local, tcwat_local, &
     lcwat_local, totcw_local, repartht_local, qtend_local, ttend_local, ltend_local)

   implicit none

   integer, intent(in) :: ncol_local
   real(r8), intent(in) :: dtime_local
   real(r8), intent(in) :: state_q_local(pcols,pver,pcnst), state_t_local(pcols,pver), qcwat_local(pcols,pver), &
        tcwat_local(pcols,pver), lcwat_local(pcols,pver), totcw_local(pcols,pver)
   real(r8), intent(inout) :: repartht_local(pcols,pver)
   real(r8), intent(inout) :: qtend_local(pcols,pver), ttend_local(pcols,pver), ltend_local(pcols,pver)

   real(r8) :: rdtime_local

   rdtime_local = 1._r8/dtime_local
   repartht_local(:ncol_local,:pver) = (latice/dtime_local) * (state_q_local(:ncol_local,:pver,ixcldice) - &
        repartht_local(:ncol_local,:pver))
   qtend_local(:ncol_local,:pver) = (state_q_local(:ncol_local,:pver,1) - qcwat_local(:ncol_local,:pver)) * rdtime_local
   ttend_local(:ncol_local,:pver) = (state_t_local(:ncol_local,:pver) - tcwat_local(:ncol_local,:pver)) * rdtime_local
   ltend_local(:ncol_local,:pver) = (totcw_local(:ncol_local,:pver) - lcwat_local(:ncol_local,:pver)) * rdtime_local

end subroutine stratiform_forcing_prep_native

!===============================================================================

subroutine stratiform_microphys_tend_diag(ncol_local, dtime_local, cld_local, fice_local, qme_local, nevapr_local, &
     evapheat_local, prfzheat_local, meltheat_local, ice2pr_local, liq2pr_local, state_q_local, prec_pcw_local, snow_pcw_local, &
     repartht_local, ptend_q_local, ptend_s_local, ast_local, icimr_local, icwmr_local, cmeheat_local, cmeice_local, cmeliq_local)

   use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

   implicit none

   integer, intent(in) :: ncol_local
   real(r8), intent(in) :: dtime_local
   real(r8), target, intent(in) :: cld_local(pcols,pver), fice_local(pcols,pver), qme_local(pcols,pver), nevapr_local(pcols,pver), &
        evapheat_local(pcols,pver), prfzheat_local(pcols,pver), meltheat_local(pcols,pver), ice2pr_local(pcols,pver), &
        liq2pr_local(pcols,pver), state_q_local(pcols,pver,pcnst), repartht_local(pcols,pver)
   real(r8), target, intent(inout) :: prec_pcw_local(pcols), snow_pcw_local(pcols), ptend_q_local(pcols,pver,pcnst)
   real(r8), target, intent(inout) :: ptend_s_local(pcols,pver), ast_local(pcols,pver), icimr_local(pcols,pver), &
        icwmr_local(pcols,pver), cmeheat_local(pcols,pver), cmeice_local(pcols,pver), cmeliq_local(pcols,pver)

   interface
      subroutine stratiform_microphys_tend_diag_codon(ncol_c, pcols_c, pver_c, pcnst_c, ixcldliq_c, ixcldice_c, dtime_c, &
           latvap_c, latice_c, cld_p, fice_p, qme_p, nevapr_p, evapheat_p, prfzheat_p, meltheat_p, ice2pr_p, liq2pr_p, &
           state_q_p, prec_pcw_p, snow_pcw_p, repartht_p, ptend_q_p, ptend_s_p, ast_p, icimr_p, icwmr_p, cmeheat_p, &
           cmeice_p, cmeliq_p) bind(c, name="stratiform_microphys_tend_diag_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pcnst_c, ixcldliq_c, ixcldice_c
         real(c_double), value :: dtime_c, latvap_c, latice_c
         type(c_ptr), value :: cld_p, fice_p, qme_p, nevapr_p, evapheat_p, prfzheat_p, meltheat_p, ice2pr_p, liq2pr_p
         type(c_ptr), value :: state_q_p, prec_pcw_p, snow_pcw_p, repartht_p, ptend_q_p, ptend_s_p, ast_p, icimr_p, icwmr_p
         type(c_ptr), value :: cmeheat_p, cmeice_p, cmeliq_p
      end subroutine stratiform_microphys_tend_diag_codon
   end interface

   if (use_native_microphys_shell_impl) then
      call stratiform_microphys_tend_diag_native(ncol_local, dtime_local, cld_local, fice_local, qme_local, nevapr_local, &
           evapheat_local, prfzheat_local, meltheat_local, ice2pr_local, liq2pr_local, state_q_local, prec_pcw_local, &
           snow_pcw_local, repartht_local, ptend_q_local, ptend_s_local, ast_local, icimr_local, icwmr_local, cmeheat_local, &
           cmeice_local, cmeliq_local)
      return
   end if

   call stratiform_microphys_tend_diag_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
        int(pcnst, c_int64_t), int(ixcldliq, c_int64_t), int(ixcldice, c_int64_t), dtime_local, latvap, latice, &
        c_loc(cld_local), c_loc(fice_local), c_loc(qme_local), c_loc(nevapr_local), c_loc(evapheat_local), c_loc(prfzheat_local), &
        c_loc(meltheat_local), c_loc(ice2pr_local), c_loc(liq2pr_local), c_loc(state_q_local), c_loc(prec_pcw_local), &
        c_loc(snow_pcw_local), c_loc(repartht_local), c_loc(ptend_q_local), c_loc(ptend_s_local), c_loc(ast_local), &
        c_loc(icimr_local), c_loc(icwmr_local), c_loc(cmeheat_local), c_loc(cmeice_local), c_loc(cmeliq_local))

end subroutine stratiform_microphys_tend_diag

!===============================================================================

subroutine stratiform_microphys_tend_diag_native(ncol_local, dtime_local, cld_local, fice_local, qme_local, nevapr_local, &
     evapheat_local, prfzheat_local, meltheat_local, ice2pr_local, liq2pr_local, state_q_local, prec_pcw_local, snow_pcw_local, &
     repartht_local, ptend_q_local, ptend_s_local, ast_local, icimr_local, icwmr_local, cmeheat_local, cmeice_local, cmeliq_local)

   implicit none

   integer, intent(in) :: ncol_local
   real(r8), intent(in) :: dtime_local
   real(r8), intent(in) :: cld_local(pcols,pver), fice_local(pcols,pver), qme_local(pcols,pver), nevapr_local(pcols,pver), &
        evapheat_local(pcols,pver), prfzheat_local(pcols,pver), meltheat_local(pcols,pver), ice2pr_local(pcols,pver), &
        liq2pr_local(pcols,pver), state_q_local(pcols,pver,pcnst), repartht_local(pcols,pver)
   real(r8), intent(inout) :: prec_pcw_local(pcols), snow_pcw_local(pcols), ptend_q_local(pcols,pver,pcnst)
   real(r8), intent(inout) :: ptend_s_local(pcols,pver), ast_local(pcols,pver), icimr_local(pcols,pver), icwmr_local(pcols,pver), &
        cmeheat_local(pcols,pver), cmeice_local(pcols,pver), cmeliq_local(pcols,pver)

   integer :: i, k

   do k = 1, pver
      do i = 1, ncol_local
         ptend_s_local(i,k) = qme_local(i,k)*(latvap + latice*fice_local(i,k)) + evapheat_local(i,k) + prfzheat_local(i,k) + &
              meltheat_local(i,k) + repartht_local(i,k)
         ptend_q_local(i,k,1) = -qme_local(i,k) + nevapr_local(i,k)
         ptend_q_local(i,k,ixcldice) = qme_local(i,k)*fice_local(i,k) - ice2pr_local(i,k)
         ptend_q_local(i,k,ixcldliq) = qme_local(i,k)*(1._r8-fice_local(i,k)) - liq2pr_local(i,k)
      end do
   end do

   do k = 1, pver
      do i = 1, ncol_local
         ast_local(i,k) = cld_local(i,k)
         icimr_local(i,k) = (state_q_local(i,k,ixcldice) + dtime_local*ptend_q_local(i,k,ixcldice)) / max(0.01_r8, ast_local(i,k))
         icwmr_local(i,k) = (state_q_local(i,k,ixcldliq) + dtime_local*ptend_q_local(i,k,ixcldliq)) / max(0.01_r8, ast_local(i,k))
      end do
   end do

   snow_pcw_local(:ncol_local) = snow_pcw_local(:ncol_local)/1000._r8
   prec_pcw_local(:ncol_local) = prec_pcw_local(:ncol_local)/1000._r8

   do k = 1, pver
      do i = 1, ncol_local
         cmeheat_local(i,k) = qme_local(i,k) * (latvap + latice*fice_local(i,k))
         cmeice_local(i,k) = qme_local(i,k) * fice_local(i,k)
         cmeliq_local(i,k) = qme_local(i,k) * (1._r8 - fice_local(i,k))
      end do
   end do

end subroutine stratiform_microphys_tend_diag_native

!===============================================================================

subroutine stratiform_cloud_mixing_diag(ncol_local, state_q_local, cld_local, concld_local, mr_ccliq_local, mr_ccice_local, &
     mr_lsliq_local, mr_lsice_local)

   use iso_c_binding, only: c_int64_t, c_loc, c_ptr

   implicit none

   integer, intent(in) :: ncol_local
   real(r8), target, intent(in) :: state_q_local(pcols,pver,pcnst), cld_local(pcols,pver), concld_local(pcols,pver)
   real(r8), target, intent(inout) :: mr_ccliq_local(pcols,pver), mr_ccice_local(pcols,pver), mr_lsliq_local(pcols,pver), &
        mr_lsice_local(pcols,pver)

   interface
      subroutine stratiform_cloud_mixing_diag_codon(ncol_c, pcols_c, pver_c, pcnst_c, ixcldliq_c, ixcldice_c, state_q_p, cld_p, &
           concld_p, mr_ccliq_p, mr_ccice_p, mr_lsliq_p, mr_lsice_p) bind(c, name="stratiform_cloud_mixing_diag_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pcnst_c, ixcldliq_c, ixcldice_c
         type(c_ptr), value :: state_q_p, cld_p, concld_p, mr_ccliq_p, mr_ccice_p, mr_lsliq_p, mr_lsice_p
      end subroutine stratiform_cloud_mixing_diag_codon
   end interface

   if (use_native_microphys_shell_impl) then
      call stratiform_cloud_mixing_diag_native(ncol_local, state_q_local, cld_local, concld_local, mr_ccliq_local, &
           mr_ccice_local, mr_lsliq_local, mr_lsice_local)
      return
   end if

   call stratiform_cloud_mixing_diag_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
        int(pcnst, c_int64_t), int(ixcldliq, c_int64_t), int(ixcldice, c_int64_t), c_loc(state_q_local), c_loc(cld_local), &
        c_loc(concld_local), c_loc(mr_ccliq_local), c_loc(mr_ccice_local), c_loc(mr_lsliq_local), c_loc(mr_lsice_local))

end subroutine stratiform_cloud_mixing_diag

!===============================================================================

subroutine stratiform_cloud_mixing_diag_native(ncol_local, state_q_local, cld_local, concld_local, mr_ccliq_local, mr_ccice_local, &
     mr_lsliq_local, mr_lsice_local)

   implicit none

   integer, intent(in) :: ncol_local
   real(r8), intent(in) :: state_q_local(pcols,pver,pcnst), cld_local(pcols,pver), concld_local(pcols,pver)
   real(r8), intent(inout) :: mr_ccliq_local(pcols,pver), mr_ccice_local(pcols,pver), mr_lsliq_local(pcols,pver), &
        mr_lsice_local(pcols,pver)

   integer :: i, k

   do k = 1, pver
      do i = 1, ncol_local
         if (cld_local(i,k) .gt. 0._r8) then
            mr_ccliq_local(i,k) = (state_q_local(i,k,ixcldliq)/cld_local(i,k))*concld_local(i,k)
            mr_ccice_local(i,k) = (state_q_local(i,k,ixcldice)/cld_local(i,k))*concld_local(i,k)
            mr_lsliq_local(i,k) = (state_q_local(i,k,ixcldliq)/cld_local(i,k))*(cld_local(i,k)-concld_local(i,k))
            mr_lsice_local(i,k) = (state_q_local(i,k,ixcldice)/cld_local(i,k))*(cld_local(i,k)-concld_local(i,k))
         else
            mr_ccliq_local(i,k) = 0._r8
            mr_ccice_local(i,k) = 0._r8
            mr_lsliq_local(i,k) = 0._r8
            mr_lsice_local(i,k) = 0._r8
         end if
      end do
   end do

end subroutine stratiform_cloud_mixing_diag_native

!===============================================================================

subroutine stratiform_postcloud_diag(ncol_local, state_q_local, state_pmid_local, state_t_local, rhcloud_local, iwc_local, &
     lwc_local, icimr_local, icwmr_local)

   use iso_c_binding, only: c_int64_t, c_loc, c_ptr

   implicit none

   integer, intent(in) :: ncol_local
   real(r8), target, intent(in) :: state_q_local(pcols,pver,pcnst), state_pmid_local(pcols,pver), state_t_local(pcols,pver), &
        rhcloud_local(pcols,pver)
   real(r8), target, intent(inout) :: iwc_local(pcols,pver), lwc_local(pcols,pver), icimr_local(pcols,pver), icwmr_local(pcols,pver)

   interface
      subroutine stratiform_postcloud_diag_codon(ncol_c, pcols_c, pver_c, pcnst_c, ixcldliq_c, ixcldice_c, state_q_p, &
           state_pmid_p, state_t_p, rhcloud_p, iwc_p, lwc_p, icimr_p, icwmr_p) bind(c, name="stratiform_postcloud_diag_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pcnst_c, ixcldliq_c, ixcldice_c
         type(c_ptr), value :: state_q_p, state_pmid_p, state_t_p, rhcloud_p, iwc_p, lwc_p, icimr_p, icwmr_p
      end subroutine stratiform_postcloud_diag_codon
   end interface

   if (use_native_impl) then
      call stratiform_postcloud_diag_native(ncol_local, state_q_local, state_pmid_local, state_t_local, rhcloud_local, iwc_local, &
           lwc_local, icimr_local, icwmr_local)
      return
   end if

   call stratiform_postcloud_diag_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
        int(pcnst, c_int64_t), int(ixcldliq, c_int64_t), int(ixcldice, c_int64_t), c_loc(state_q_local), &
        c_loc(state_pmid_local), c_loc(state_t_local), c_loc(rhcloud_local), c_loc(iwc_local), c_loc(lwc_local), &
        c_loc(icimr_local), c_loc(icwmr_local))

end subroutine stratiform_postcloud_diag

!===============================================================================

subroutine stratiform_postcloud_diag_native(ncol_local, state_q_local, state_pmid_local, state_t_local, rhcloud_local, iwc_local, &
     lwc_local, icimr_local, icwmr_local)

   implicit none

   integer, intent(in) :: ncol_local
   real(r8), intent(in) :: state_q_local(pcols,pver,pcnst), state_pmid_local(pcols,pver), state_t_local(pcols,pver), &
        rhcloud_local(pcols,pver)
   real(r8), intent(inout) :: iwc_local(pcols,pver), lwc_local(pcols,pver), icimr_local(pcols,pver), icwmr_local(pcols,pver)

   integer :: i, k

   do k = 1, pver
      do i = 1, ncol_local
         iwc_local(i,k) = state_q_local(i,k,ixcldice)*state_pmid_local(i,k)/(287.15_r8*state_t_local(i,k))
         lwc_local(i,k) = state_q_local(i,k,ixcldliq)*state_pmid_local(i,k)/(287.15_r8*state_t_local(i,k))
         icimr_local(i,k) = state_q_local(i,k,ixcldice) / max(0.01_r8,rhcloud_local(i,k))
         icwmr_local(i,k) = state_q_local(i,k,ixcldliq) / max(0.01_r8,rhcloud_local(i,k))
      end do
   end do

end subroutine stratiform_postcloud_diag_native

!===============================================================================

subroutine stratiform_store_oldcloud(ncol_local, state_q_local, state_t_local, qcwat_local, tcwat_local, lcwat_local)

   use iso_c_binding, only: c_int64_t, c_loc, c_ptr

   implicit none

   integer, intent(in) :: ncol_local
   real(r8), target, intent(in) :: state_q_local(pcols,pver,pcnst), state_t_local(pcols,pver)
   real(r8), target, intent(inout) :: qcwat_local(pcols,pver), tcwat_local(pcols,pver), lcwat_local(pcols,pver)

   interface
      subroutine stratiform_store_oldcloud_codon(ncol_c, pcols_c, pver_c, pcnst_c, ixcldliq_c, ixcldice_c, state_q_p, &
           state_t_p, qcwat_p, tcwat_p, lcwat_p) bind(c, name="stratiform_store_oldcloud_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pcnst_c, ixcldliq_c, ixcldice_c
         type(c_ptr), value :: state_q_p, state_t_p, qcwat_p, tcwat_p, lcwat_p
      end subroutine stratiform_store_oldcloud_codon
   end interface

   if (use_native_impl) then
      call stratiform_store_oldcloud_native(ncol_local, state_q_local, state_t_local, qcwat_local, tcwat_local, lcwat_local)
      return
   end if

   call stratiform_store_oldcloud_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
        int(pcnst, c_int64_t), int(ixcldliq, c_int64_t), int(ixcldice, c_int64_t), c_loc(state_q_local), c_loc(state_t_local), &
        c_loc(qcwat_local), c_loc(tcwat_local), c_loc(lcwat_local))

end subroutine stratiform_store_oldcloud

!===============================================================================

subroutine stratiform_store_oldcloud_native(ncol_local, state_q_local, state_t_local, qcwat_local, tcwat_local, lcwat_local)

   implicit none

   integer, intent(in) :: ncol_local
   real(r8), intent(in) :: state_q_local(pcols,pver,pcnst), state_t_local(pcols,pver)
   real(r8), intent(inout) :: qcwat_local(pcols,pver), tcwat_local(pcols,pver), lcwat_local(pcols,pver)

   integer :: k

   do k = 1, pver
      qcwat_local(:ncol_local,k) = state_q_local(:ncol_local,k,1)
      tcwat_local(:ncol_local,k) = state_t_local(:ncol_local,k)
      lcwat_local(:ncol_local,k) = state_q_local(:ncol_local,k,ixcldice) + state_q_local(:ncol_local,k,ixcldliq)
   end do

end subroutine stratiform_store_oldcloud_native

!===============================================================================

subroutine stratiform_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('STRATIFORM_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_impl = .false.
   end if

   impl_selected = .true.

   if (masterproc) then
      if (use_native_impl) then
         write(iulog,*) 'stratiform_tend implementation = native'
         write(*,*) 'stratiform_tend implementation = native'
         call stratiform_append_impl_proof('STRATIFORM_PROOF_FILE', &
              'stratiform_tend implementation = native')
      else
         write(iulog,*) 'stratiform_tend implementation = codon'
         write(*,*) 'stratiform_tend implementation = codon'
         call stratiform_append_impl_proof('STRATIFORM_PROOF_FILE', &
              'stratiform_tend implementation = codon')
      end if
      call flush(iulog)
   end if

end subroutine stratiform_select_impl

!===============================================================================

subroutine stratiform_select_branches(use_shfrc_in, cam3_in)

   use iso_c_binding, only: c_int64_t, c_loc, c_ptr

   logical, intent(in) :: use_shfrc_in
   logical, intent(in) :: cam3_in

   integer(c_int64_t), target :: branch_mask_c

   interface
      subroutine stratiform_select_branches_codon(use_shfrc_c, cam3_c, branch_mask_p) &
           bind(c, name="stratiform_select_branches_codon")
        use iso_c_binding, only: c_int64_t, c_ptr
        integer(c_int64_t), value :: use_shfrc_c, cam3_c
        type(c_ptr), value :: branch_mask_p
      end subroutine stratiform_select_branches_codon
   end interface

   if (branch_selected) return

   branch_mask_c = 0_c_int64_t
   call stratiform_select_branches_codon( &
        merge(1_c_int64_t, 0_c_int64_t, use_shfrc_in), &
        merge(1_c_int64_t, 0_c_int64_t, cam3_in), &
        c_loc(branch_mask_c) &
   )

   branch_mask = int(branch_mask_c)
   branch_selected = .true.

end subroutine stratiform_select_branches

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

   subroutine debug_microphys_1(state1,ptend,i,k, &
        dtime,qme,fice,snow_pcw,prec_pcw, &
        prain,nevapr,prodsnow, evapsnow, &
        ice2pr,liq2pr,liq2snow)

     use physics_types, only: physics_state, physics_ptend
     use physconst,     only: tmelt

     implicit none
     
     integer, intent(in) :: i,k
     type(physics_state), intent(in) :: state1   ! local copy of the state variable
     type(physics_ptend), intent(in) :: ptend  ! local copy of the ptend variable
     real(r8), intent(in)  :: dtime                ! timestep
     real(r8), intent(in) :: qme(pcols,pver)          ! local condensation - evaporation of cloud water

     real(r8), intent(in) :: prain(pcols,pver)          ! local production of precipitation
     real(r8), intent(in) :: nevapr(pcols,pver)          ! local evaporation of precipitation
     real(r8), intent(in) :: prodsnow(pcols,pver)          ! local production of snow
     real(r8), intent(in) :: evapsnow(pcols,pver)          ! local evaporation of snow
     real(r8), intent(in) :: ice2pr(pcols,pver)   ! rate of conversion of ice to precip
     real(r8), intent(in) :: liq2pr(pcols,pver)   ! rate of conversion of liquid to precip
     real(r8), intent(in) :: liq2snow(pcols,pver)   ! rate of conversion of liquid to snow
     real(r8), intent(in) :: fice    (pcols,pver)          ! Fractional ice content within cloud
     real(r8), intent(in) :: snow_pcw(pcols)
     real(r8), intent(in) :: prec_pcw(pcols)

     real(r8) hs1, qv1, ql1, qi1, qs1, qr1, fice2, pr1, w1, w2, w3, fliq, res
     real(r8) w4, wl, wv, wi, wlf, wvf, wif, qif, qlf, qvf

     pr1 = 0
     hs1 = 0
     qv1 = 0
     ql1 = 0
     qi1 = 0
     qs1 = 0
     qr1 = 0
     w1 = 0
     wl = 0
     wv = 0
     wi = 0
     wlf = 0
     wvf = 0 
     wif = 0


     write(iulog,*) 
     write(iulog,*) ' input state, t, q, l, i ', k, state1%t(i,k), state1%q(i,k,1), state1%q(i,k,ixcldliq),  state1%q(i,k,ixcldice)
     write(iulog,*) ' rain, snow, total from components before accumulation ', qr1, qs1, qr1+qs1
     write(iulog,*) ' total precip before accumulation                      ', k, pr1

     wv = wv + state1%q(i,k,1       )*state1%pdel(i,k)/gravit
     wl = wl + state1%q(i,k,ixcldliq)*state1%pdel(i,k)/gravit
     wi = wi + state1%q(i,k,ixcldice)*state1%pdel(i,k)/gravit

     qvf = state1%q(i,k,1) + ptend%q(i,k,1)*dtime
     qlf = state1%q(i,k,ixcldliq) + ptend%q(i,k,ixcldliq)*dtime
     qif = state1%q(i,k,ixcldice) + ptend%q(i,k,ixcldice)*dtime

     if (qvf.lt.0._r8) then
        write(iulog,*) ' qvf is negative *******', qvf
     endif
     if (qlf.lt.0._r8) then
        write(iulog,*) ' qlf is negative *******', qlf
     endif
     if (qif.lt.0._r8) then
        write(iulog,*) ' qif is negative *******', qif
     endif
     write(iulog,*) ' qvf, qlf, qif ', qvf, qlf, qif

     wvf = wvf + qvf*state1%pdel(i,k)/gravit
     wlf = wlf + qlf*state1%pdel(i,k)/gravit
     wif = wif + qif*state1%pdel(i,k)/gravit

     hs1 = hs1 + ptend%s(i,k)*state1%pdel(i,k)/gravit
     pr1 = pr1 + state1%pdel(i,k)/gravit*(prain(i,k)-nevapr(i,k))
     qv1 = qv1 - (qme(i,k)-nevapr(i,k))*state1%pdel(i,k)/gravit    ! vdot
     w1  = w1  + (qme(i,k)-prain(i,k))*state1%pdel(i,k)/gravit    ! cdot
     qi1 = qi1 + ((qme(i,k))*fice(i,k)        -ice2pr(i,k) )*state1%pdel(i,k)/gravit   ! idot
     ql1 = ql1 + ((qme(i,k))*(1._r8-fice(i,k))-liq2pr(i,k) )*state1%pdel(i,k)/gravit   ! ldot

     qr1 = qr1 &
          + ( liq2pr(i,k)-liq2snow(i,k)   &     ! production of rain
          -(nevapr(i,k)-evapsnow(i,k)) &     ! rain evaporation
          )*state1%pdel(i,k)/gravit
     qs1 = qs1 &
          + ( ice2pr(i,k) + liq2snow(i,k) &     ! production of snow.Note last term has phase change
          -evapsnow(i,k)               &     ! snow evaporation
          )*state1%pdel(i,k)/gravit

     if (state1%t(i,k).gt.tmelt) then
        qr1 = qr1 + qs1
        qs1 = 0._r8
     endif
     write(iulog,*) ' rain, snow, total after accumulation ', qr1, qs1, qr1+qs1
     write(iulog,*) ' total precip after accumulation      ', k, pr1
     write(iulog,*)
     write(iulog,*) ' layer prain, nevapr, pdel ', prain(i,k), nevapr(i,k), state1%pdel(i,k)
     write(iulog,*) ' layer prodsnow, ice2pr+liq2snow ', prodsnow(i,k), ice2pr(i,k)+liq2snow(i,k)
     write(iulog,*) ' layer prain-prodsnow, liq2pr-liq2snow ', prain(i,k)-prodsnow(i,k), liq2pr(i,k)-liq2snow(i,k)
     write(iulog,*) ' layer evapsnow, evaprain ', k, evapsnow(i,k), nevapr(i,k)-evapsnow(i,k)
     write(iulog,*) ' layer ice2pr, liq2pr, liq2snow ', ice2pr(i,k), liq2pr(i,k), liq2snow(i,k)
     write(iulog,*) ' layer ice2pr+liq2pr, prain ', ice2pr(i,k)+liq2pr(i,k), prain(i,k)
     write(iulog,*)
     write(iulog,*) ' qv1 vapor removed from col after accum  (vdot)   ', k, qv1
     write(iulog,*) ' - (precip produced - vapor removed) after accum  ', k, -pr1-qv1
     write(iulog,*) ' condensate produce after accum                   ', k, w1
     write(iulog,*) ' liq+ice tends accum                              ', k, ql1+qi1
     write(iulog,*) ' change in total water after accum                ', k, qv1+ql1+qi1
     write(iulog,*) ' imbalance in colum after accum                   ', k, qs1+qr1+qv1+ql1+qi1
     write(iulog,*) ' fice at this lev ', fice(i,k)
     write(iulog,*)

     res = abs((qs1+qr1+qv1+ql1+qi1)/max(abs(qv1),abs(ql1),abs(qi1),abs(qs1),abs(qr1),1.e-36_r8))
     write(iulog,*) ' relative residual in column method 1             ', k, res

     write(iulog,*) ' relative residual in column method 2             ',&
	 k, abs((qs1+qr1+qv1+ql1+qi1)/max(abs(qv1+ql1+qi1),1.e-36_r8))
     !            if (abs((qs1+qr1+qv1+ql1+qi1)/(qs1+qr1+1.e-36)).gt.1.e-14) then
     if (res.gt.1.e-14_r8) then
        call endrun ('STRATIFORM_TEND')
     endif

     !             w3  = qme(i,k) * (latvap + latice*fice(i,k)) &
     !               + evapheat(i,k) + prfzheat(i,k) + meltheat(i,k)

     res = qs1+qr1-pr1
     w4 = max(abs(qs1),abs(qr1),abs(pr1)) 
     if (w4.gt.0._r8)  then
        if (res/w4.gt.1.e-14_r8) then
           write(iulog,*) ' imbalance in precips calculated two ways '
           write(iulog,*) ' res/w4, pr1, qr1, qs1, qr1+qs1 ', &
                res/w4, pr1, qr1, qs1, qr1+qs1
           !                   call endrun()
        endif
     endif
     if (k.eq.pver) then
        write(iulog,*) ' pcond returned precip, rain and snow rates ', prec_pcw(i), prec_pcw(i)-snow_pcw(i), snow_pcw(i)
        write(iulog,*) ' I calculate ', pr1, qr1, qs1
        !               call endrun
        write(iulog,*) ' byrons water check ', wv+wl+wi-pr1*dtime, wvf+wlf+wif
     endif
     write(iulog,*)


   end subroutine debug_microphys_1

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

   subroutine debug_microphys_2(state1,&
        snow_pcw,fsaut,fsacw ,fsaci, meltheat)

     use ppgrid,        only: pver
     use physconst,     only: tmelt
     use physics_types, only: physics_state
     
     implicit none

     type(physics_state), intent(in) :: state1   ! local copy of the state variable
     real(r8), intent(in) :: snow_pcw(pcols)
     real(r8), intent(in) :: fsaut(pcols,pver)              
     real(r8), intent(in) :: fsacw(pcols,pver)              
     real(r8), intent(in) :: fsaci(pcols,pver)              
     real(r8), intent(in) :: meltheat(pcols,pver)          ! heating rate due to phase change of precip


     integer  i,ncol,lchnk


     ncol = state1%ncol
     lchnk = state1%lchnk
     
     do i = 1,ncol
        if (snow_pcw(i) .gt. 0.01_r8/8.64e4_r8  .and.  state1%t(i,pver) .gt. tmelt) then
           write(iulog,*) ' stratiform: snow, temp, ', i, lchnk, &
                snow_pcw(i), state1%t(i,pver)
           write(iulog,*) ' t ', state1%t(i,:)
           write(iulog,*) ' fsaut ', fsaut(i,:)
           write(iulog,*) ' fsacw ', fsacw(i,:)
           write(iulog,*) ' fsaci ', fsaci(i,:)
           write(iulog,*) ' meltheat ', meltheat(i,:)
           call endrun ('STRATIFORM_TEND')
        endif
        
        if (snow_pcw(i)*8.64e4_r8 .lt. -1.e-5_r8) then
           write(iulog,*) ' neg snow ', snow_pcw(i)*8.64e4_r8
           write(iulog,*) ' stratiform: snow_pcw, temp, ', i, lchnk, &
                snow_pcw(i), state1%t(i,pver)
           write(iulog,*) ' t ', state1%t(i,:)
           write(iulog,*) ' fsaut ', fsaut(i,:)
           write(iulog,*) ' fsacw ', fsacw(i,:)
           write(iulog,*) ' fsaci ', fsaci(i,:)
           write(iulog,*) ' meltheat ', meltheat(i,:)
           call endrun ('STRATIFORM_TEND')
        endif
     end do
     
   end subroutine debug_microphys_2

  end module stratiform
