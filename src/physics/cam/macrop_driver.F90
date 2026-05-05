  module macrop_driver

  !-------------------------------------------------------------------------------------------------------
  ! Purpose:
  !
  ! Provides the CAM interface to the prognostic cloud macrophysics
  !
  ! Author: Andrew Gettelman, Cheryl Craig October 2010
  ! Origin: modified from stratiform.F90 elements 
  !    (Boville 2002, Coleman 2004, Park 2009, Kay 2010)
  !-------------------------------------------------------------------------------------------------------

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use spmd_utils,    only: masterproc
  use ppgrid,        only: pcols, pver, pverp
  use physconst,     only: latice, latvap
  use phys_control,  only: phys_getopts
  use constituents,  only: cnst_get_ind, pcnst
  use physics_buffer,    only: physics_buffer_desc, pbuf_set_field, pbuf_get_field, pbuf_old_tim_idx
  use time_manager,      only: is_first_step
  use cldwat2m_macro,    only: ini_macro
  use perf_mod,          only: t_startf, t_stopf
  use cam_logfile,       only: iulog
  use cam_abortutils,    only: endrun

  implicit none
  private
  save

  public :: macrop_driver_readnl
  public :: macrop_driver_register
  public :: macrop_driver_init
  public :: macrop_driver_tend
  public :: ice_macro_tend

  logical, public :: do_cldice             ! .true., park macrophysics is prognosing cldice
  logical, public :: do_cldliq             ! .true., park macrophysics is prognosing cldliq
  logical, public :: do_detrain            ! .true., park macrophysics is detraining ice into stratiform

  ! ------------------------- !
  ! Private Module Parameters !
  ! ------------------------- !

  ! 'cu_det_st' : If .true. (.false.), detrain cumulus liquid condensate into the pre-existing liquid stratus 
  !               (environment) without (with) macrophysical evaporation. If there is no pre-esisting stratus, 
  !               evaporate cumulus liquid condensate. This option only influences the treatment of cumulus
  !               liquid condensate, not cumulus ice condensate.

  logical, parameter :: cu_det_st  = .false.  

  logical :: micro_do_icesupersat

  ! Parameters used for selecting generalized critical RH for liquid and ice stratus
  integer :: rhminl_opt = 0
  integer :: rhmini_opt = 0


  character(len=16) :: shallow_scheme
  logical           :: use_shfrc                       ! Local copy of flag from convect_shallow_use_shfrc

  integer :: &
    ixcldliq,     &! cloud liquid amount index
    ixcldice,     &! cloud ice amount index
    ixnumliq,     &! cloud liquid number index
    ixnumice,     &! cloud ice water index
    qcwat_idx,    &! qcwat index in physics buffer
    lcwat_idx,    &! lcwat index in physics buffer
    iccwat_idx,   &! iccwat index in physics buffer
    nlwat_idx,    &! nlwat index in physics buffer
    niwat_idx,    &! niwat index in physics buffer
    tcwat_idx,    &! tcwat index in physics buffer
    CC_T_idx,     &!
    CC_qv_idx,    &!
    CC_ql_idx,    &!
    CC_qi_idx,    &!
    CC_nl_idx,    &!
    CC_ni_idx,    &!
    CC_qlst_idx,  &!
    cld_idx,      &! cld index in physics buffer
    ast_idx,      &! stratiform cloud fraction index in physics buffer
    aist_idx,     &! ice stratiform cloud fraction index in physics buffer
    alst_idx,     &! liquid stratiform cloud fraction index in physics buffer
    qist_idx,     &! ice stratiform in-cloud IWC 
    qlst_idx,     &! liquid stratiform in-cloud LWC  
    concld_idx,   &! concld index in physics buffer
    fice_idx,     &  
    cmeliq_idx,   &  
    shfrc_idx,    &
    naai_idx 

  integer :: &
    tke_idx = -1,       &! tke defined at the model interfaces
    qtl_flx_idx = -1,   &! overbar(w'qtl' where qtl = qv + ql) from the PBL scheme
    qti_flx_idx = -1,   &! overbar(w'qti' where qti = qv + qi) from the PBL scheme
    cmfr_det_idx = -1,  &! detrained convective mass flux from UNICON
    qlr_det_idx = -1,   &! detrained convective ql from UNICON  
    qir_det_idx = -1     ! detrained convective qi from UNICON  

  logical :: use_native_impl = .false.
  logical :: impl_selected = .false.
  integer :: branch_mask = 0
  logical :: branch_selected = .false.
  logical :: wtrc_detrain_impl_logged = .false.
  logical :: use_native_wtrc_shell_impl = .false.
  logical :: wtrc_shell_impl_selected = .false.
  logical :: use_native_wtrc_process_impl = .false.
  logical :: wtrc_process_impl_selected = .false.
  logical :: mmacro_input_shell_logged = .false.
  logical :: mmacro_post_fields_shell_logged = .false.
  logical :: cfmip_diag_shell_logged = .false.
  logical :: detrain_init_shell_logged = .false.
  logical :: detrain_post_shell_logged = .false.
  logical :: mmacro_config_check_logged = .false.
  logical :: detrain_core_logged = .false.
  logical :: clr_old_diag_logged = .false.
  logical :: forcing_prep_logged = .false.
  logical :: ptend_assign_logged = .false.
  logical :: store_state_logged = .false.

  contains

  ! ===============================================================================
  subroutine macrop_driver_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Namelist variables
   logical  :: macro_park_do_cldice  = .true.   ! do_cldice = .true., park macrophysics is prognosing cldice
   logical  :: macro_park_do_cldliq  = .true.   ! do_cldliq = .true., park macrophysics is prognosing cldliq
   logical  :: macro_park_do_detrain = .true.   ! do_detrain = .true., park macrophysics is detraining ice into stratiform

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'macrop_driver_readnl'

   namelist /macro_park_nl/ macro_park_do_cldice, macro_park_do_cldliq, macro_park_do_detrain
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'macro_park_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, macro_park_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)

      ! set local variables

      do_cldice  = macro_park_do_cldice
      do_cldliq  = macro_park_do_cldliq
      do_detrain = macro_park_do_detrain

   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(do_cldice,             1, mpilog, 0, mpicom)
   call mpibcast(do_cldliq,             1, mpilog, 0, mpicom)
   call mpibcast(do_detrain,            1, mpilog, 0, mpicom)
#endif

end subroutine macrop_driver_readnl

  !================================================================================================

  subroutine macrop_driver_register

  !---------------------------------------------------------------------- !
  !                                                                       !
  ! Register the constituents (cloud liquid and cloud ice) and the fields !
  ! in the physics buffer.                                                !
  !                                                                       !
  !---------------------------------------------------------------------- !

   
   use physics_buffer, only : pbuf_add_field, dtype_r8, dyn_time_lvls

  !-----------------------------------------------------------------------

    call phys_getopts(shallow_scheme_out=shallow_scheme)

    call pbuf_add_field('AST',      'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), ast_idx)
    call pbuf_add_field('AIST',     'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), aist_idx)
    call pbuf_add_field('ALST',     'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), alst_idx)
    call pbuf_add_field('QIST',     'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), qist_idx)
    call pbuf_add_field('QLST',     'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), qlst_idx)
    call pbuf_add_field('CLD',      'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cld_idx)
    call pbuf_add_field('CONCLD',   'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), concld_idx)

    call pbuf_add_field('QCWAT',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), qcwat_idx)
    call pbuf_add_field('LCWAT',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), lcwat_idx)
    call pbuf_add_field('ICCWAT',   'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), iccwat_idx)
    call pbuf_add_field('NLWAT',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), nlwat_idx)
    call pbuf_add_field('NIWAT',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), niwat_idx)
    call pbuf_add_field('TCWAT',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), tcwat_idx)

    call pbuf_add_field('FICE',     'physpkg', dtype_r8, (/pcols,pver/), fice_idx)

    call pbuf_add_field('CMELIQ',   'physpkg', dtype_r8, (/pcols,pver/), cmeliq_idx)

  end subroutine macrop_driver_register

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine macrop_driver_init(pbuf2d)

  !-------------------------------------------- !
  !                                             !
  ! Initialize the cloud water parameterization !
  !                                             ! 
  !-------------------------------------------- !
    use physics_buffer, only : pbuf_get_index
    use cam_history,     only: addfld, add_default, phys_decomp
    use convect_shallow, only: convect_shallow_use_shfrc
    
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    logical              :: history_aerosol      ! Output the MAM aerosol tendencies
    logical              :: history_budget       ! Output tendencies and state variables for CAM4
                                                 ! temperature, water vapor, cloud ice and cloud
                                                 ! liquid budgets.
    integer              :: history_budget_histfile_num ! output history file number for budget fields
    integer :: istat
    character(len=*), parameter :: subname = 'macrop_driver_init'
    !-----------------------------------------------------------------------

    ! Initialization routine for cloud macrophysics
    if (shallow_scheme .eq. 'UNICON') rhminl_opt = 1
    call ini_macro(rhminl_opt, rhmini_opt)

    call phys_getopts(history_aerosol_out              = history_aerosol      , &
                      history_budget_out               = history_budget       , &
                      history_budget_histfile_num_out  = history_budget_histfile_num, &
                      micro_do_icesupersat_out         = micro_do_icesupersat)

  ! Find out whether shfrc from convect_shallow will be used in cldfrc

    if( convect_shallow_use_shfrc() ) then
        use_shfrc = .true.
        shfrc_idx = pbuf_get_index('shfrc')
    else 
        use_shfrc = .false.
    endif

    call addfld ('DPDLFLIQ ', 'kg/kg/s ', pver, 'A', 'Detrained liquid water from deep convection'             ,phys_decomp)
    call addfld ('DPDLFICE ', 'kg/kg/s ', pver, 'A', 'Detrained ice from deep convection'                      ,phys_decomp)
    call addfld ('SHDLFLIQ ', 'kg/kg/s ', pver, 'A', 'Detrained liquid water from shallow convection'          ,phys_decomp)
    call addfld ('SHDLFICE ', 'kg/kg/s ', pver, 'A', 'Detrained ice from shallow convection'                   ,phys_decomp)
    call addfld ('DPDLFT   ', 'K/s     ', pver, 'A', 'T-tendency due to deep convective detrainment'           ,phys_decomp)
    call addfld ('SHDLFT   ', 'K/s     ', pver, 'A', 'T-tendency due to shallow convective detrainment'        ,phys_decomp)

    call addfld ('ZMDLF    ', 'kg/kg/s ', pver, 'A', 'Detrained liquid water from ZM convection'               ,phys_decomp)

    call addfld ('MACPDT   ', 'W/kg    ', pver, 'A', 'Heating tendency - Revised  macrophysics'                ,phys_decomp)
    call addfld ('MACPDQ   ', 'kg/kg/s ', pver, 'A', 'Q tendency - Revised macrophysics'                       ,phys_decomp)
    call addfld ('MACPDLIQ ', 'kg/kg/s ', pver, 'A', 'CLDLIQ tendency - Revised macrophysics'                  ,phys_decomp)
    call addfld ('MACPDICE ', 'kg/kg/s ', pver, 'A', 'CLDICE tendency - Revised macrophysics'                  ,phys_decomp)

    call addfld ('CLDVAPADJ', 'kg/kg/s ', pver, 'A', &
         'Q tendency associated with liq/ice adjustment - Revised macrophysics' ,phys_decomp)
    call addfld ('CLDLIQADJ', 'kg/kg/s ', pver, 'A', 'CLDLIQ adjustment tendency - Revised macrophysics'       ,phys_decomp)
    call addfld ('CLDICEADJ', 'kg/kg/s ', pver, 'A', 'CLDICE adjustment tendency - Revised macrophysics'       ,phys_decomp)
    call addfld ('CLDLIQDET', 'kg/kg/s ', pver, 'A', &
         'Detrainment of conv cld liq into envrionment  - Revised macrophysics' ,phys_decomp)
    call addfld ('CLDICEDET', 'kg/kg/s ', pver, 'A', &
         'Detrainment of conv cld ice into envrionment  - Revised macrophysics' ,phys_decomp)
    call addfld ('CLDLIQLIM', 'kg/kg/s ', pver, 'A', 'CLDLIQ limiting tendency - Revised macrophysics'         ,phys_decomp)
    call addfld ('CLDICELIM', 'kg/kg/s ', pver, 'A', 'CLDICE limiting tendency - Revised macrophysics'         ,phys_decomp)

    call addfld ('AST',       '1',        pver, 'A', 'Stratus cloud fraction',                                  phys_decomp)
    call addfld ('LIQCLDF',   '1',        pver, 'A', 'Stratus Liquid cloud fraction',                           phys_decomp)
    call addfld ('ICECLDF',   '1',        pver, 'A', 'Stratus ICE cloud fraction',                              phys_decomp)

    call addfld ('CLDST    ', 'fraction', pver, 'A', 'Stratus cloud fraction'                                  ,phys_decomp)
    call addfld ('CONCLD   ', 'fraction', pver, 'A', 'Convective cloud cover'                                  ,phys_decomp)
 
    call addfld ('CLR_LIQ',   'fraction', pver, 'A', 'Clear sky fraction for liquid stratus'  , phys_decomp)
    call addfld ('CLR_ICE',   'fraction', pver, 'A', 'Clear sky fraction for ice stratus'     , phys_decomp)

    call addfld ('CLDLIQSTR   ', 'kg/kg', pver, 'A', 'Stratiform CLDLIQ'                                  ,phys_decomp)
    call addfld ('CLDICESTR   ', 'kg/kg', pver, 'A', 'Stratiform CLDICE'                                  ,phys_decomp)
    call addfld ('CLDLIQCON   ', 'kg/kg', pver, 'A', 'Convective CLDLIQ'                                  ,phys_decomp)
    call addfld ('CLDICECON   ', 'kg/kg', pver, 'A', 'Convective CLDICE'                                  ,phys_decomp)

    call addfld ('CLDSICE  ', 'kg/kg   ', pver, 'A', 'CloudSat equivalent ice mass mixing ratio'               ,phys_decomp)
    call addfld ('CMELIQ   ', 'kg/kg/s ', pver, 'A', 'Rate of cond-evap of liq within the cloud'               ,phys_decomp)

    call addfld ('TTENDICE',      'K/s ', pver, 'A', 'T tendency from Ice Saturation Adjustment'       ,phys_decomp)
    call addfld ('QVTENDICE', 'kg/kg/s ', pver, 'A', 'Q tendency from Ice Saturation Adjustment'       ,phys_decomp)
    call addfld ('QITENDICE', 'kg/kg/s ', pver, 'A', 'CLDICE tendency from Ice Saturation Adjustment'       ,phys_decomp)
    call addfld ('NITENDICE', 'kg/kg/s ', pver, 'A', 'NUMICE tendency from Ice Saturation Adjustment'       ,phys_decomp)
    if ( history_budget ) then

          call add_default ('DPDLFLIQ ', history_budget_histfile_num, ' ')
          call add_default ('DPDLFICE ', history_budget_histfile_num, ' ')
          call add_default ('SHDLFLIQ ', history_budget_histfile_num, ' ')
          call add_default ('SHDLFICE ', history_budget_histfile_num, ' ')
          call add_default ('DPDLFT   ', history_budget_histfile_num, ' ')
          call add_default ('SHDLFT   ', history_budget_histfile_num, ' ')
          call add_default ('ZMDLF    ', history_budget_histfile_num, ' ')

          call add_default ('MACPDT   ', history_budget_histfile_num, ' ')
          call add_default ('MACPDQ   ', history_budget_histfile_num, ' ')
          call add_default ('MACPDLIQ ', history_budget_histfile_num, ' ')
          call add_default ('MACPDICE ', history_budget_histfile_num, ' ')
 
          call add_default ('CLDVAPADJ', history_budget_histfile_num, ' ')
          call add_default ('CLDLIQLIM', history_budget_histfile_num, ' ')
          call add_default ('CLDLIQDET', history_budget_histfile_num, ' ')
          call add_default ('CLDLIQADJ', history_budget_histfile_num, ' ')
          call add_default ('CLDICELIM', history_budget_histfile_num, ' ')
          call add_default ('CLDICEDET', history_budget_histfile_num, ' ')
          call add_default ('CLDICEADJ', history_budget_histfile_num, ' ')

          call add_default ('CMELIQ   ', history_budget_histfile_num, ' ')

    end if

    ! Get constituent indices
    call cnst_get_ind('CLDLIQ', ixcldliq)
    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('NUMLIQ', ixnumliq)
    call cnst_get_ind('NUMICE', ixnumice)

    ! Get physics buffer indices
    CC_T_idx    = pbuf_get_index('CC_T')
    CC_qv_idx   = pbuf_get_index('CC_qv')
    CC_ql_idx   = pbuf_get_index('CC_ql')
    CC_qi_idx   = pbuf_get_index('CC_qi')
    CC_nl_idx   = pbuf_get_index('CC_nl')
    CC_ni_idx   = pbuf_get_index('CC_ni')
    CC_qlst_idx = pbuf_get_index('CC_qlst')

    if (micro_do_icesupersat) then 
       naai_idx      = pbuf_get_index('NAAI')
    endif

    if (rhminl_opt > 0 .or. rhmini_opt > 0) then
       cmfr_det_idx = pbuf_get_index('cmfr_det', istat)
       if (istat < 0) call endrun(subname//': macrop option requires cmfr_det in pbuf')
       if (rhminl_opt > 0) then
          qlr_det_idx  = pbuf_get_index('qlr_det', istat)
          if (istat < 0) call endrun(subname//': macrop option requires qlr_det in pbuf')
       end if
       if (rhmini_opt > 0) then
          qir_det_idx  = pbuf_get_index('qir_det', istat)
          if (istat < 0) call endrun(subname//': macrop option requires qir_det in pbuf')
       end if
    end if

    if (rhminl_opt == 2 .or. rhmini_opt == 2) then
       tke_idx = pbuf_get_index('tke')
       if (rhminl_opt == 2) then
          qtl_flx_idx = pbuf_get_index('qtl_flx', istat)
          if (istat < 0) call endrun(subname//': macrop option requires qtl_flx in pbuf')
       end if
       if (rhmini_opt == 2) then
          qti_flx_idx = pbuf_get_index('qti_flx', istat)
          if (istat < 0) call endrun(subname//': macrop option requires qti_flx in pbuf')
       end if
    end if

    ! Init pbuf fields.  Note that the fields CLD, CONCLD, QCWAT, LCWAT, 
    ! ICCWAT, and TCWAT are initialized in phys_inidat.
    if (is_first_step()) then
       call pbuf_set_field(pbuf2d, ast_idx,    0._r8)
       call pbuf_set_field(pbuf2d, aist_idx,   0._r8)
       call pbuf_set_field(pbuf2d, alst_idx,   0._r8)
       call pbuf_set_field(pbuf2d, qist_idx,   0._r8)
       call pbuf_set_field(pbuf2d, qlst_idx,   0._r8)
       call pbuf_set_field(pbuf2d, nlwat_idx,  0._r8)
       call pbuf_set_field(pbuf2d, niwat_idx,  0._r8)
       call pbuf_set_field(pbuf2d, fice_idx,   0._r8)
       call pbuf_set_field(pbuf2d, cmeliq_idx, 0._r8)
    end if

  end subroutine macrop_driver_init

  !============================================================================ !
  !                                                                             !
  !============================================================================ !


  subroutine macrop_driver_tend(                             &
             state, ptend, dtime, landfrac,  &
             ocnfrac,  snowh,                       &
             dlf, dlf2,wtdlf, cmfmc, cmfmc2, ts,          &
             sst, zdu,       &
             pbuf, &
             det_s, det_ice)

  !-------------------------------------------------------- !  
  !                                                         ! 
  ! Purpose:                                                !
  !                                                         !
  ! Interface to detrain, cloud fraction and                !
  !     cloud macrophysics subroutines                      !
  !                                                         ! 
  ! Author: A. Gettelman, C. Craig, Oct 2010                !
  ! based on stratiform_tend by D.B. Coleman 4/2010         !
  !                                                         !
  !-------------------------------------------------------- !

  use cloud_fraction,   only: cldfrc, cldfrc_fice
  use physics_types,    only: physics_state, physics_ptend
  use physics_types,    only: physics_ptend_init, physics_update
  use physics_types,    only: physics_ptend_sum,  physics_state_copy
  use physics_types,    only: physics_state_dealloc
  use cam_history,      only: outfld
  use constituents,     only: cnst_get_ind, pcnst
  use cldwat2m_macro,   only: mmacro_pcond
  use physconst,        only: cpair, tmelt, gravit
  use time_manager,     only: get_nstep
  use water_tracer_vars,only: trace_water, wtrc_detrain_in_macrop, wtrc_nwset,&
                              wtrc_iatype, iwspec, wtrc_indices, wtrc_ncnst
  use water_tracers,    only: wtrc_init_rates, wtrc_add_rates, wtrc_apply_rates
  use water_types,      only: pwtype, iwtvap, iwtliq, iwtice

  use ref_pres,         only: top_lev => trop_cloud_top_lev

  !
  ! Input arguments
  !

  type(physics_state), intent(in)    :: state       ! State variables
  type(physics_ptend), intent(out)   :: ptend       ! macrophysics parameterization tendencies
  type(physics_buffer_desc), pointer :: pbuf(:)     ! Physics buffer

  real(r8), intent(in)  :: dtime                    ! Timestep
  real(r8), intent(in)  :: landfrac(pcols)          ! Land fraction (fraction)
  real(r8), intent(in)  :: ocnfrac (pcols)          ! Ocean fraction (fraction)
  real(r8), intent(in)  :: snowh(pcols)             ! Snow depth over land, water equivalent (m)
  real(r8), intent(in)  :: dlf(pcols,pver)          ! Detrained water from convection schemes
  real(r8), intent(in)  :: dlf2(pcols,pver)         ! Detrained water from shallow convection scheme
  real(r8), intent(in)  :: cmfmc(pcols,pverp)       ! Deep + Shallow Convective mass flux [ kg /s/m^2 ]
  real(r8), intent(in)  :: cmfmc2(pcols,pverp)      ! Shallow convective mass flux [ kg/s/m^2 ]

  real(r8), intent(in)  :: ts(pcols)                ! Surface temperature
  real(r8), intent(in)  :: sst(pcols)               ! Sea surface temperature
  real(r8), intent(in)  :: zdu(pcols,pver)          ! Detrainment rate from deep convection

  ! water tracers:
  real(r8), intent(in)  :: wtdlf(pcols,pver,wtrc_nwset) !detrained water tracers from convection [kg/kg/s]

  ! These two variables are needed for energy check    
  real(r8), intent(out) :: det_s(pcols)             ! Integral of detrained static energy from ice
  real(r8), intent(out) :: det_ice(pcols)           ! Integral of detrained ice for energy check

  !
  ! Local variables
  !

  type(physics_state) :: state_loc                  ! Local copy of the state variable
  type(physics_ptend) :: ptend_loc                  ! Local parameterization tendencies

  integer :: lchnk                                  ! Chunk identifier
  integer :: ncol                                   ! Number of atmospheric columns

  ! Physics buffer fields

  integer itim_old
  real(r8), pointer, dimension(:,:) :: qcwat        ! Cloud water old q
  real(r8), pointer, dimension(:,:) :: tcwat        ! Cloud water old temperature
  real(r8), pointer, dimension(:,:) :: lcwat        ! Cloud liquid water old q
  real(r8), pointer, dimension(:,:) :: iccwat       ! Cloud ice water old q
  real(r8), pointer, dimension(:,:) :: nlwat        ! Cloud liquid droplet number condentration. old.
  real(r8), pointer, dimension(:,:) :: niwat        ! Cloud ice    droplet number condentration. old.
  real(r8), pointer, dimension(:,:) :: CC_T         ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_qv        ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_ql        ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_qi        ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_nl        ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_ni        ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_qlst      ! In-liquid stratus microphysical tendency
  real(r8), pointer, dimension(:,:) :: cld          ! Total cloud fraction
  real(r8), pointer, dimension(:,:) :: ast          ! Relative humidity cloud fraction
  real(r8), pointer, dimension(:,:) :: aist         ! Physical ice stratus fraction
  real(r8), pointer, dimension(:,:) :: alst         ! Physical liquid stratus fraction
  real(r8), pointer, dimension(:,:) :: qist         ! Physical in-cloud IWC
  real(r8), pointer, dimension(:,:) :: qlst         ! Physical in-cloud LWC
  real(r8), pointer, dimension(:,:) :: concld       ! Convective cloud fraction

  real(r8), pointer, dimension(:,:) :: shfrc        ! Cloud fraction from shallow convection scheme

  real(r8), pointer, dimension(:,:) :: cmeliq

  real(r8), pointer, dimension(:,:) :: tke
  real(r8), pointer, dimension(:,:) :: qtl_flx
  real(r8), pointer, dimension(:,:) :: qti_flx
  real(r8), pointer, dimension(:,:) :: cmfr_det
  real(r8), pointer, dimension(:,:) :: qlr_det
  real(r8), pointer, dimension(:,:) :: qir_det

  ! Convective cloud to the physics buffer for purposes of ql contrib. to radn.

  real(r8), pointer, dimension(:,:) :: fice_ql      ! Cloud ice/water partitioning ratio.

  real(r8), pointer, dimension(:,:) :: naai         ! Number concentration of activated ice nuclei
 
  real(r8) :: latsub

  ! tendencies for ice saturation adjustment
  real(r8)  :: stend(pcols,pver)
  real(r8)  :: qvtend(pcols,pver)
  real(r8)  :: qitend(pcols,pver)
  real(r8)  :: initend(pcols,pver)

  ! Local variables for cldfrc

  real(r8)  cldst(pcols,pver)                       ! Stratus cloud fraction
  real(r8)  rhcloud(pcols,pver)                     ! Relative humidity cloud (last timestep)
  real(r8)  clc(pcols)                              ! Column convective cloud amount
  real(r8)  rhu00(pcols,pver)                       ! RH threshold for cloud
  real(r8)  icecldf(pcols,pver)                     ! Ice cloud fraction
  real(r8)  liqcldf(pcols,pver)                     ! Liquid cloud fraction (combined into cloud)
  real(r8)  relhum(pcols,pver)                      ! RH, output to determine drh/da

  ! Local variables for macrophysics

  real(r8)  rdtime                                  ! 1./dtime
  real(r8)  qtend(pcols,pver)                       ! Moisture tendencies
  real(r8)  ttend(pcols,pver)                       ! Temperature tendencies
  real(r8)  ltend(pcols,pver)                       ! Cloud liquid water tendencies
  real(r8)  fice(pcols,pver)                        ! Fractional ice content within cloud
  real(r8)  fsnow(pcols,pver)                       ! Fractional snow production
  real(r8)  homoo(pcols,pver)  
  real(r8)  qcreso(pcols,pver)  
  real(r8)  prcio(pcols,pver)  
  real(r8)  praio(pcols,pver)  
  real(r8)  qireso(pcols,pver)
  real(r8)  ftem(pcols,pver)
  real(r8)  pracso (pcols,pver) 
  real(r8)  dpdlfliq(pcols,pver)
  real(r8)  dpdlfice(pcols,pver)
  real(r8)  shdlfliq(pcols,pver)
  real(r8)  shdlfice(pcols,pver)
  real(r8)  dpdlft  (pcols,pver)
  real(r8)  shdlft  (pcols,pver)

  real(r8)  qc(pcols,pver)
  real(r8)  qi(pcols,pver)
  real(r8)  nc(pcols,pver)
  real(r8)  ni(pcols,pver)

  logical   lq(pcnst)

  ! Output from mmacro_pcond

  real(r8)  tlat(pcols,pver)
  real(r8)  qvlat(pcols,pver)
  real(r8)  qcten(pcols,pver)
  real(r8)  qiten(pcols,pver)
  real(r8)  ncten(pcols,pver)
  real(r8)  niten(pcols,pver)

  ! Output from mmacro_pcond

  real(r8)  qvadj(pcols,pver)                       ! Macro-physics adjustment tendency from "positive_moisture" call (vapor)
  real(r8)  qladj(pcols,pver)                       ! Macro-physics adjustment tendency from "positive_moisture" call (liquid)
  real(r8)  qiadj(pcols,pver)                       ! Macro-physics adjustment tendency from "positive_moisture" call (ice)
  real(r8)  qllim(pcols,pver)                       ! Macro-physics tendency from "instratus_condensate" call (liquid)
  real(r8)  qilim(pcols,pver)                       ! Macro-physics tendency from "instratus_condensate" call (ice)

  ! For revised macophysics, mmacro_pcond

  real(r8)  itend(pcols,pver)
  real(r8)  lmitend(pcols,pver)
  real(r8)  zeros(pcols,pver)
  real(r8)  t_inout(pcols,pver)
  real(r8)  qv_inout(pcols,pver)
  real(r8)  ql_inout(pcols,pver)
  real(r8)  qi_inout(pcols,pver)
  real(r8)  concld_old(pcols,pver)

  ! Note that below 'clr_old' is defined using 'alst_old' not 'ast_old' for full consistency with the 
  ! liquid condensation process which is using 'alst' not 'ast'. 
  ! For microconsistency use 'concld_old', since 'alst_old' was computed using 'concld_old'.
  ! Since convective updraft fractional area is small, it does not matter whether 'concld' or 'concld_old' is used.
  ! Note also that 'clri_old' is defined using 'ast_old' since current microphysics is operating on 'ast_old' 
  real(r8)  clrw_old(pcols,pver) ! (1 - concld_old - alst_old)
  real(r8)  clri_old(pcols,pver) ! (1 - concld_old -  ast_old)

  real(r8)  nl_inout(pcols,pver)
  real(r8)  ni_inout(pcols,pver)

  real(r8)  nltend(pcols,pver)
  real(r8)  nitend(pcols,pver)


  ! For detraining cumulus condensate into the 'stratus' without evaporation
  ! This is for use in mmacro_pcond

  real(r8)  dlf_T(pcols,pver)
  real(r8)  dlf_qv(pcols,pver)
  real(r8)  dlf_ql(pcols,pver)
  real(r8)  dlf_qi(pcols,pver)
  real(r8)  dlf_nl(pcols,pver)
  real(r8)  dlf_ni(pcols,pver)

  ! Local variables for CFMIP calculations
  real(r8) :: mr_lsliq(pcols,pver)  ! mixing_ratio_large_scale_cloud_liquid (kg/kg)
  real(r8) :: mr_lsice(pcols,pver)  ! mixing_ratio_large_scale_cloud_ice (kg/kg)
  real(r8) :: mr_ccliq(pcols,pver)  ! mixing_ratio_convective_cloud_liquid (kg/kg)
  real(r8) :: mr_ccice(pcols,pver)  ! mixing_ratio_convective_cloud_ice (kg/kg)

  ! CloudSat equivalent ice mass mixing ratio (kg/kg)
  real(r8) :: cldsice(pcols,pver)

  ! Local variables for water tracers/isotopes
  real(r8)              :: process_rates(pcols,pver,pwtype,pwtype,pwtype) ! Process rates (kg/kg/sec)
  integer               :: m                                              ! water set index
  logical               :: isOk
  integer               :: iwtype
  logical               :: micro_do_icesupersat_local
  logical               :: trace_water_local
  logical               :: wtrc_detrain_in_macrop_local
  logical               :: cu_det_st_local
  logical               :: use_shfrc_local
  logical               :: do_cldice_local
  logical               :: do_cldliq_local
  logical               :: do_detrain_local

  real(r8) pqctn(pcols,pver)
  real(r8) nqctn(pcols,pver)
  real(r8) pqitn(pcols,pver)
  real(r8) nqitn(pcols,pver)
  integer :: mmacro_config_mask

  ! ======================================================================

  call macrop_driver_select_impl()
  if (.not. use_native_impl) then
     call macrop_driver_select_branches(micro_do_icesupersat, trace_water, wtrc_detrain_in_macrop, &
          cu_det_st, use_shfrc, do_cldice, do_cldliq, do_detrain)
     micro_do_icesupersat_local = iand(branch_mask, 1) /= 0
     trace_water_local = iand(branch_mask, 2) /= 0
     wtrc_detrain_in_macrop_local = iand(branch_mask, 4) /= 0
     cu_det_st_local = iand(branch_mask, 8) /= 0
     use_shfrc_local = iand(branch_mask, 16) /= 0
     do_cldice_local = iand(branch_mask, 32) /= 0
     do_cldliq_local = iand(branch_mask, 64) /= 0
     do_detrain_local = iand(branch_mask, 128) /= 0
  else
     micro_do_icesupersat_local = micro_do_icesupersat
     trace_water_local = trace_water
     wtrc_detrain_in_macrop_local = wtrc_detrain_in_macrop
     cu_det_st_local = cu_det_st
     use_shfrc_local = use_shfrc
     do_cldice_local = do_cldice
     do_cldliq_local = do_cldliq
     do_detrain_local = do_detrain
  end if

  if (micro_do_icesupersat_local) then 
     call pbuf_get_field(pbuf, naai_idx, naai)
  endif

  lchnk = state%lchnk
  ncol  = state%ncol

  call physics_state_copy(state, state_loc)            ! Copy state to local state_loc.

  ! Associate pointers with physics buffer fields

  itim_old = pbuf_old_tim_idx()

  call pbuf_get_field(pbuf, qcwat_idx,   qcwat,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, tcwat_idx,   tcwat,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, lcwat_idx,   lcwat,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, iccwat_idx,  iccwat,  start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, nlwat_idx,   nlwat,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, niwat_idx,   niwat,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

  call pbuf_get_field(pbuf, cc_t_idx,    cc_t,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, cc_qv_idx,   cc_qv,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, cc_ql_idx,   cc_ql,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, cc_qi_idx,   cc_qi,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, cc_nl_idx,   cc_nl,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, cc_ni_idx,   cc_ni,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, cc_qlst_idx, cc_qlst, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

  call pbuf_get_field(pbuf, cld_idx,     cld,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, concld_idx,  concld, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, ast_idx,     ast,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, aist_idx,    aist,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, alst_idx,    alst,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, qist_idx,    qist,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, qlst_idx,    qlst,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

  call pbuf_get_field(pbuf, cmeliq_idx,  cmeliq)

! For purposes of convective ql.

  call pbuf_get_field(pbuf, fice_idx,     fice_ql )


  ! Initialize convective detrainment tendency

  call macrop_driver_detrain_init_shell(ncol, dlf_T, dlf_qv, dlf_ql, dlf_qi, dlf_nl, dlf_ni, det_s, det_ice, &
       dpdlfliq, dpdlfice, shdlfliq, shdlfice, dpdlft, shdlft)

   ! ------------------------------------- !
   ! From here, process computation begins ! 
   ! ------------------------------------- !

   ! ----------------------------------------------------------------------------- !
   ! Detrainment of convective condensate into the environment or stratiform cloud !
   ! ----------------------------------------------------------------------------- !

   lq(:)        = .FALSE.
   lq(ixcldliq) = .TRUE.
   lq(ixcldice) = .TRUE.
   lq(ixnumliq) = .TRUE.
   lq(ixnumice) = .TRUE.

   !allow water tracers to change:
   if ((trace_water_local) .and. (wtrc_detrain_in_macrop_local)) then
     do m=1,wtrc_nwset
       lq(wtrc_iatype(m,iwtliq)) = .TRUE.
       lq(wtrc_iatype(m,iwtice)) = .TRUE.
     end do  
   end if

   call physics_ptend_init(ptend_loc, state%psetcols, 'pcwdetrain', ls=.true., lq=lq)   ! Initialize local physics_ptend object

     ! Procedures :
     ! (1) Partition detrained convective cloud water into liquid and ice based on T.
     !     This also involves heating.
     !     If convection scheme can handle this internally, this step is not necssary.
     ! (2) Assuming a certain effective droplet radius, computes number concentration
     !     of detrained convective cloud liquid and ice.
     ! (3) If 'cu_det_st = .true' ('false'), detrain convective cloud 'liquid' into 
     !     the pre-existing 'liquid' stratus ( mean environment ).  The former does
     !     not involve any macrophysical evaporation while the latter does. This is
     !     a kind of 'targetted' deposition. Then, force in-stratus LWC to be bounded 
     !     by qcst_min and qcst_max in mmacro_pcond.
     ! (4) In contrast to liquid, convective ice is detrained into the environment 
     !     and involved in the sublimation. Similar bounds as liquid stratus are imposed.
     ! This is the key procesure generating upper-level cirrus clouds.
     ! The unit of dlf : [ kg/kg/s ]

   call macrop_driver_detrain_core(ncol, do_detrain_local, cu_det_st_local, state_loc%t, state_loc%pdel, dlf, dlf2, ptend_loc%q(:,:,ixcldliq), &
        ptend_loc%q(:,:,ixcldice), ptend_loc%q(:,:,ixnumliq), ptend_loc%q(:,:,ixnumice), ptend_loc%s, det_s, det_ice, &
        dlf_t, dlf_qv, dlf_ql, dlf_qi, dlf_nl, dlf_ni, dpdlfliq, dpdlfice, shdlfliq, shdlfice, dpdlft, shdlft)

   if ((trace_water_local) .and. (wtrc_detrain_in_macrop_local)) then
      call macrop_driver_wtrc_detrain(ncol, state_loc%t, wtdlf, ptend_loc%q)
   end if

   call outfld( 'DPDLFLIQ ', dpdlfliq, pcols, lchnk )
   call outfld( 'DPDLFICE ', dpdlfice, pcols, lchnk )
   call outfld( 'SHDLFLIQ ', shdlfliq, pcols, lchnk )
   call outfld( 'SHDLFICE ', shdlfice, pcols, lchnk )
   call outfld( 'DPDLFT   ', dpdlft  , pcols, lchnk )
   call outfld( 'SHDLFT   ', shdlft  , pcols, lchnk )

   call outfld( 'ZMDLF',     dlf     , pcols, state_loc%lchnk )

   call macrop_driver_detrain_post_shell(ncol, det_ice)  ! divide by density of water

   ! Add the detrainment tendency to the output tendency
   call physics_ptend_init(ptend, state%psetcols, 'macrop')
   call physics_ptend_sum(ptend_loc, ptend, ncol)

   ! update local copy of state with the detrainment tendency
   ! ptend_loc is reset to zero by this call
   call physics_update(state_loc, ptend_loc, dtime)

   if (micro_do_icesupersat_local) then 

      ! -------------------------------------- !
      ! Ice Saturation Adjustment Computation  !
      ! -------------------------------------- !

      lq(:)        = .FALSE.

      lq(1)        = .true.
      lq(ixcldice) = .true.
      lq(ixnumice) = .true.

      latsub = latvap + latice

      call physics_ptend_init(ptend_loc, state%psetcols, 'iceadj', ls=.true., lq=lq)

      stend(:ncol,:)=0._r8
      qvtend(:ncol,:)=0._r8
      qitend(:ncol,:)=0._r8
      initend(:ncol,:)=0._r8

      call ice_macro_tend(naai(:ncol,top_lev:pver),state%t(:ncol,top_lev:pver), &
           state%pmid(:ncol,top_lev:pver),state%q(:ncol,top_lev:pver,1),state%q(:ncol,top_lev:pver,ixcldice),&
           state%q(:ncol,top_lev:pver,ixnumice),latsub,dtime,&
           stend(:ncol,top_lev:pver),qvtend(:ncol,top_lev:pver),qitend(:ncol,top_lev:pver),&
           initend(:ncol,top_lev:pver))

      ! update local copy of state with the tendencies
      ptend_loc%q(:ncol,top_lev:pver,1)=qvtend(:ncol,top_lev:pver)
      ptend_loc%q(:ncol,top_lev:pver,ixcldice)=qitend(:ncol,top_lev:pver)  
      ptend_loc%q(:ncol,top_lev:pver,ixnumice)=initend(:ncol,top_lev:pver)
      ptend_loc%s(:ncol,top_lev:pver)=stend(:ncol,top_lev:pver) 

      ! Add the ice tendency to the output tendency
      call physics_ptend_sum(ptend_loc, ptend, ncol)
 
      ! ptend_loc is reset to zero by this call
      call physics_update(state_loc, ptend_loc, dtime)

      ! Write output for tendencies:
      call outfld( 'TTENDICE',  stend/cpair, pcols, lchnk )
      call outfld( 'QVTENDICE', qvtend, pcols, lchnk )
      call outfld( 'QITENDICE', qitend, pcols, lchnk )
      call outfld( 'NITENDICE', initend, pcols, lchnk )

   endif

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

   concld_old(:ncol,top_lev:pver) = concld(:ncol,top_lev:pver)

   nullify(tke, qtl_flx, qti_flx, cmfr_det, qlr_det, qir_det)
   if (tke_idx      > 0) call pbuf_get_field(pbuf, tke_idx, tke)
   if (qtl_flx_idx  > 0) call pbuf_get_field(pbuf, qtl_flx_idx,  qtl_flx)
   if (qti_flx_idx  > 0) call pbuf_get_field(pbuf, qti_flx_idx,  qti_flx)
   if (cmfr_det_idx > 0) call pbuf_get_field(pbuf, cmfr_det_idx, cmfr_det)
   if (qlr_det_idx  > 0) call pbuf_get_field(pbuf, qlr_det_idx,  qlr_det)
   if (qir_det_idx  > 0) call pbuf_get_field(pbuf, qir_det_idx,  qir_det)

   clrw_old(:ncol,:top_lev-1) = 0._r8
   clri_old(:ncol,:top_lev-1) = 0._r8
   call macrop_driver_clr_old_diag(ncol, concld, alst, ast, clrw_old, clri_old)

   if( use_shfrc_local ) then
       call pbuf_get_field(pbuf, shfrc_idx, shfrc )
   else 
       allocate(shfrc(pcols,pver))
       shfrc(:,:) = 0._r8
   endif

   ! CAM5 only uses 'concld' output from the below subroutine. 
   ! Stratus ('ast' = max(alst,aist)) and total cloud fraction ('cld = ast + concld')
   ! will be computed using this updated 'concld' in the stratiform macrophysics 
   ! scheme (mmacro_pcond) later below. 

   call t_startf("cldfrc")

   call cldfrc( lchnk, ncol, pbuf,                                                 &
                state_loc%pmid, state_loc%t, state_loc%q(:,:,1), state_loc%omega,  &
                state_loc%phis, shfrc, use_shfrc_local,                            &
                cld, rhcloud, clc, state_loc%pdel,                                 &
                cmfmc, cmfmc2, landfrac,snowh, concld, cldst,                      &
                ts, sst, state_loc%pint(:,pverp), zdu, ocnfrac, rhu00,             &
                state_loc%q(:,:,ixcldice), icecldf, liqcldf,                       &
                relhum, 0 )

   call t_stopf("cldfrc")

   ! ---------------------------------------------- !
   ! Stratiform Cloud Macrophysics and Microphysics !
   ! ---------------------------------------------- !

   lchnk  = state_loc%lchnk
   ncol   = state_loc%ncol
   rdtime = 1._r8/dtime

 ! Define fractional amount of stratus condensate and precipitation in ice phase.
 ! This uses a ramp ( -30 ~ -10 for fice, -5 ~ 0 for fsnow ). 
 ! The ramp within convective cloud may be different

   call cldfrc_fice( ncol, state_loc%t, fice, fsnow )


   lq(:)        = .FALSE.

   lq(1)        = .true.
   lq(ixcldice) = .true.
   lq(ixcldliq) = .true.

   lq(ixnumliq) = .true.
   lq(ixnumice) = .true.

   !Water tracers:
   do m=1,wtrc_ncnst
     lq(wtrc_indices(m)) = .true.
   end do

   ! Initialize local physics_ptend object again
   call physics_ptend_init(ptend_loc, state%psetcols, 'macro_park', &
        ls=.true., lq=lq )  

 ! --------------------------------- !
 ! Liquid Macrop_Driver Macrophysics !
 ! --------------------------------- !

   call t_startf('mmacro_pcond')

   call macrop_driver_mmacro_input_shell(ncol, state_loc%q, zeros, qc, qi, nc, ni)

 ! In CAM5, 'microphysical forcing' ( CC_... ) and 'the other advective forcings' ( ttend, ... ) 
 ! are separately provided into the prognostic microp_driver macrophysics scheme. This is an
 ! attempt to resolve in-cloud and out-cloud forcings. 

   call macrop_driver_forcing_prep(ncol, get_nstep(), rdtime, state_loc%t, state_loc%q(:,:,1), qc, qi, nc, ni, tcwat, &
        qcwat, lcwat, iccwat, nlwat, niwat, cc_t, cc_qv, cc_ql, cc_qi, cc_nl, cc_ni, cc_qlst, ttend, qtend, ltend, &
        itend, nltend, nitend, lmitend, t_inout, qv_inout, ql_inout, qi_inout, nl_inout, ni_inout)

 ! Liquid Microp_Driver Macrophysics.
 ! The main roles of this subroutines are
 ! (1) compute net condensation rate of stratiform liquid ( cmeliq )
 ! (2) compute liquid stratus and ice stratus fractions. 
 ! Note 'ttend...' are advective tendencies except microphysical process while
 !      'CC...'    are microphysical tendencies. 

   call mmacro_pcond( lchnk, ncol, dtime, state_loc%pmid, state_loc%pdel,        &
                      t_inout, qv_inout, ql_inout, qi_inout, nl_inout, ni_inout, &                  
                      ttend, qtend, lmitend, itend, nltend, nitend,              &
                      CC_T, CC_qv, CC_ql, CC_qi, CC_nl, CC_ni, CC_qlst,          & 
                      dlf_T, dlf_qv, dlf_ql, dlf_qi, dlf_nl, dlf_ni,             &
                      concld_old, concld, clrw_old, clri_old, landfrac, snowh,   &
                      tke, qtl_flx, qti_flx, cmfr_det, qlr_det, qir_det,         &
                      tlat, qvlat, qcten, qiten, ncten, niten,                   &
                      cmeliq, qvadj, qladj, qiadj, qllim, qilim,                 &
                      cld, alst, aist, qlst, qist, do_cldice_local ) 

 ! Copy of concld/fice to put in physics buffer
 ! Below are used only for convective cloud.

   call macrop_driver_mmacro_post_fields_shell(ncol, fice, alst, aist, fice_ql, ast)

 ! Compute net stratus fraction using maximum over-lapping assumption

   call t_stopf('mmacro_pcond')

   call macrop_driver_ptend_assign(ncol, tlat, qvlat, qcten, qiten, ncten, niten, ptend_loc%s, ptend_loc%q(:,:,1), &
        ptend_loc%q(:,:,ixcldliq), ptend_loc%q(:,:,ixcldice), ptend_loc%q(:,:,ixnumliq), ptend_loc%q(:,:,ixnumice))

   mmacro_config_mask = macrop_driver_mmacro_config_check(ncol, do_cldice_local, do_cldliq_local, qiten, niten, qcten, ncten)
   if (iand(mmacro_config_mask, 1) /= 0) then
      call endrun("macrop_driver:ERROR - "// &
           "Cldwat is configured not to prognose cloud ice, but mmacro_pcond has ice mass tendencies.")
   end if
   if (iand(mmacro_config_mask, 2) /= 0) then
      call endrun("macrop_driver:ERROR -"// &
           " Cldwat is configured not to prognose cloud ice, but mmacro_pcond has ice number tendencies.")
   end if
   if (iand(mmacro_config_mask, 4) /= 0) then
      call endrun("macrop_driver:ERROR - "// &
           "Cldwat is configured not to prognose cloud liquid, but mmacro_pcond has liquid mass tendencies.")
   end if
   if (iand(mmacro_config_mask, 8) /= 0) then
      call endrun("macrop_driver:ERROR - "// &
           "Cldwat is configured not to prognose cloud liquid, but mmacro_pcond has liquid number tendencies.")
   end if

!-------------------------------------------------
!water tracers
!-------------------------------------------------

   ! If doing water isotopes, then apply these processes to the isotopic water
   ! species.
   !
   ! NOTE: The detrained ice and water vapor should be handled in the convection
   ! routines, but for now a generic tendency may be applied here based upon the.
   !
   ! Assume that the liquid and ice tendencies are from vapor<->liquid and
   ! vapor<->ice processes, not ice<->liquid. If this is wrong, then individual
   ! rates will need to be determined in mmacro_pcond.
   
   if (trace_water_local) then

     call macrop_driver_select_wtrc_shell_impl()

     if (use_native_wtrc_shell_impl) then
        call macrop_driver_select_wtrc_process_impl()

        ! Setup the process rate matrix using the pre-sedimentation state and
        ! the calculated process rates.
        !
        ! NOTE: The reverse of the process is filled in automatically.
        if (use_native_wtrc_process_impl) then
           call wtrc_init_rates(top_lev, process_rates)

           ! Processes that consume water vapor:

           !initalize variables - JN
           pqctn(:,top_lev:) = 0._r8
           nqctn(:,top_lev:) = 0._r8
           pqitn(:,top_lev:) = 0._r8
           nqitn(:,top_lev:) = 0._r8

           !split into positive and negative tendencies - JN
           call macrop_driver_wtrc_split_tend(ncol, qcten, qiten, pqctn, nqctn, pqitn, nqitn)

           call wtrc_add_rates(process_rates, ncol, top_lev, iwtvap, iwtvap, iwtvap, qvlat + qcten + qiten)

           call wtrc_add_rates(process_rates, ncol, top_lev, iwtvap, iwtliq, iwtvap, pqctn)
           call wtrc_add_rates(process_rates, ncol, top_lev, iwtvap, iwtliq, iwtliq, nqctn)
           call wtrc_add_rates(process_rates, ncol, top_lev, iwtvap, iwtice, iwtvap, pqitn)
           call wtrc_add_rates(process_rates, ncol, top_lev, iwtvap, iwtice, iwtice, nqitn)
        else
           call macrop_driver_wtrc_process_rates_codon_wrap(ncol, qvlat, qcten, qiten, process_rates)
        end if

       !debugging:
!        call wtrc_add_rates(process_rates, ncol, top_lev, iwtvap, iwtvap, qvlat, do_reverse=.false. )
!        call wtrc_add_rates(process_rates, ncol, top_lev, iwtliq, iwtliq, qcten, do_reverse=.false. )
!        call wtrc_add_rates(process_rates, ncol, top_lev, iwtice, iwtice, qiten, do_reverse=.false. )

        ! Apply these rates.
        call wtrc_apply_rates(state_loc, ptend_loc, pbuf, top_lev, dtime, .false., pre_rates=process_rates, &
                              prelat=tlat)
     else
        call macrop_driver_wtrc_shell_codon_wrap(ncol, dtime, state_loc%q, state_loc%t, state_loc%pmid, ptend_loc%q, &
             qvlat, qcten, qiten, tlat, process_rates)
     end if

   end if !water tracers

!-----------------------------------------------------

   ! update the output tendencies with the mmacro_pcond tendencies
   call physics_ptend_sum(ptend_loc, ptend, ncol)

   ! state_loc is the equlibrium state after macrophysics
   call physics_update(state_loc, ptend_loc, dtime)

   call outfld('CLR_LIQ', clrw_old,  pcols, lchnk)
   call outfld('CLR_ICE', clri_old,  pcols, lchnk)

   call outfld( 'MACPDT   ', tlat ,  pcols, lchnk )
   call outfld( 'MACPDQ   ', qvlat,  pcols, lchnk )
   call outfld( 'MACPDLIQ ', qcten,  pcols, lchnk )
   call outfld( 'MACPDICE ', qiten,  pcols, lchnk )
   call outfld( 'CLDVAPADJ', qvadj,  pcols, lchnk )
   call outfld( 'CLDLIQADJ', qladj,  pcols, lchnk )
   call outfld( 'CLDICEADJ', qiadj,  pcols, lchnk )
   call outfld( 'CLDLIQDET', dlf_ql, pcols, lchnk )
   call outfld( 'CLDICEDET', dlf_qi, pcols, lchnk )
   call outfld( 'CLDLIQLIM', qllim,  pcols, lchnk )
   call outfld( 'CLDICELIM', qilim,  pcols, lchnk )

   call outfld( 'ICECLDF ', aist,   pcols, lchnk )
   call outfld( 'LIQCLDF ', alst,   pcols, lchnk )
   call outfld( 'AST',      ast,    pcols, lchnk )

   call outfld( 'CONCLD  ', concld, pcols, lchnk )
   call outfld( 'CLDST   ', cldst,  pcols, lchnk )

   call outfld( 'CMELIQ'  , cmeliq, pcols, lchnk )


   ! calculations and outfld calls for CLDLIQSTR, CLDICESTR, CLDLIQCON, CLDICECON for CFMIP

   call macrop_driver_cfmip_diag_shell(ncol, cld, state_loc%q(:,:,ixcldliq), state_loc%q(:,:,ixcldice), &
        mr_ccliq, mr_ccice, mr_lsliq, mr_lsice)

   call outfld( 'CLDLIQSTR  ', mr_lsliq,    pcols, lchnk )
   call outfld( 'CLDICESTR  ', mr_lsice,    pcols, lchnk )
   call outfld( 'CLDLIQCON  ', mr_ccliq,    pcols, lchnk )
   call outfld( 'CLDICECON  ', mr_ccice,    pcols, lchnk )

   ! ------------------------------------------------- !
   ! Save equilibrium state variables for macrophysics !        
   ! at the next time step                             !
   ! ------------------------------------------------- !
   cldsice = 0._r8
   call macrop_driver_store_state(ncol, state_loc%t, state_loc%q(:,:,1), state_loc%q(:,:,ixcldliq), &
        state_loc%q(:,:,ixcldice), state_loc%q(:,:,ixnumliq), state_loc%q(:,:,ixnumice), tcwat, qcwat, lcwat, iccwat, &
        nlwat, niwat, cldsice)

   call outfld( 'CLDSICE'    , cldsice,   pcols, lchnk )

   ! ptend_loc is deallocated in physics_update above
   call physics_state_dealloc(state_loc)

end subroutine macrop_driver_tend

!============================================================================ !
!                                                                             !
!============================================================================ !

subroutine macrop_driver_wtrc_detrain(ncol_local, state_t_local, wtdlf_local, ptend_q_local)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use water_tracer_vars, only: wtrc_nwset, wtrc_iatype
  use water_types, only: iwtliq, iwtice
  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local
  real(r8), target, intent(in) :: state_t_local(pcols,pver), wtdlf_local(pcols,pver,wtrc_nwset)
  real(r8), target, intent(inout) :: ptend_q_local(pcols,pver,pcnst)

  integer(c_int64_t), target :: liq_type_c(wtrc_nwset), ice_type_c(wtrc_nwset)
  integer :: m

  interface
     subroutine macrop_driver_wtrc_detrain_codon(ncol_c, pcols_c, pver_c, top_lev_c, wtrc_nwset_c, state_t_p, wtdlf_p, &
          liq_type_p, ice_type_p, ptend_q_p) bind(c, name="macrop_driver_wtrc_detrain_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c, wtrc_nwset_c
       type(c_ptr), value :: state_t_p, wtdlf_p, liq_type_p, ice_type_p, ptend_q_p
     end subroutine macrop_driver_wtrc_detrain_codon
  end interface

  if (masterproc .and. .not. wtrc_detrain_impl_logged) then
     if (use_native_impl) then
        write(iulog,*) 'macrop_driver_wtrc_detrain implementation = native'
     else
        write(iulog,*) 'macrop_driver_wtrc_detrain implementation = codon'
     end if
     call flush(iulog)
     wtrc_detrain_impl_logged = .true.
  end if

  if (use_native_impl) then
     call macrop_driver_wtrc_detrain_native(ncol_local, state_t_local, wtdlf_local, wtrc_iatype(:,iwtliq), wtrc_iatype(:,iwtice), &
          ptend_q_local)
     return
  end if

  do m = 1, wtrc_nwset
     liq_type_c(m) = int(wtrc_iatype(m,iwtliq), c_int64_t)
     ice_type_c(m) = int(wtrc_iatype(m,iwtice), c_int64_t)
  end do

  call macrop_driver_wtrc_detrain_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
       int(top_lev, c_int64_t), int(wtrc_nwset, c_int64_t), c_loc(state_t_local), c_loc(wtdlf_local), c_loc(liq_type_c), &
       c_loc(ice_type_c), c_loc(ptend_q_local))

end subroutine macrop_driver_wtrc_detrain

!============================================================================ !

subroutine macrop_driver_wtrc_detrain_native(ncol_local, state_t_local, wtdlf_local, liq_type_local, ice_type_local, ptend_q_local)

  use water_tracer_vars, only: wtrc_nwset
  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local
  real(r8), intent(in) :: state_t_local(pcols,pver), wtdlf_local(pcols,pver,wtrc_nwset)
  integer, intent(in) :: liq_type_local(wtrc_nwset), ice_type_local(wtrc_nwset)
  real(r8), intent(inout) :: ptend_q_local(pcols,pver,pcnst)

  real(r8) :: dum1_local
  integer :: i, k, m

  do k = top_lev, pver
     do i = 1, ncol_local
        if( state_t_local(i,k) > 268.15_r8 ) then
           dum1_local = 0.0_r8
        elseif( state_t_local(i,k) < 238.15_r8 ) then
           dum1_local = 1.0_r8
        else
           dum1_local = ( 268.15_r8 - state_t_local(i,k) ) / 30._r8
        endif
        do m = 1, wtrc_nwset
           ptend_q_local(i,k,liq_type_local(m)) = wtdlf_local(i,k,m) * (1._r8 - dum1_local)
           ptend_q_local(i,k,ice_type_local(m)) = wtdlf_local(i,k,m) * dum1_local
        end do
     end do
  end do

end subroutine macrop_driver_wtrc_detrain_native

!============================================================================ !

subroutine macrop_driver_detrain_core(ncol_local, do_detrain_local, cu_det_st_local, state_t_local, state_pdel_local, &
     dlf_local, dlf2_local, ptend_ql_local, ptend_qi_local, ptend_nl_local, ptend_ni_local, ptend_s_local, det_s_local, &
     det_ice_local, dlf_t_local, dlf_qv_local, dlf_ql_local, dlf_qi_local, dlf_nl_local, dlf_ni_local, dpdlfliq_local, &
     dpdlfice_local, shdlfliq_local, shdlfice_local, dpdlft_local, shdlft_local)

  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
  use physconst, only: cpair, gravit
  use ref_pres,  only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local
  logical, intent(in) :: do_detrain_local, cu_det_st_local
  real(r8), target, intent(in) :: state_t_local(pcols,pver), state_pdel_local(pcols,pver)
  real(r8), target, intent(in) :: dlf_local(pcols,pver), dlf2_local(pcols,pver)
  real(r8), target, intent(inout) :: ptend_ql_local(pcols,pver), ptend_qi_local(pcols,pver)
  real(r8), target, intent(inout) :: ptend_nl_local(pcols,pver), ptend_ni_local(pcols,pver)
  real(r8), target, intent(inout) :: ptend_s_local(pcols,pver), det_s_local(pcols), det_ice_local(pcols)
  real(r8), target, intent(inout) :: dlf_t_local(pcols,pver), dlf_qv_local(pcols,pver), dlf_ql_local(pcols,pver)
  real(r8), target, intent(inout) :: dlf_qi_local(pcols,pver), dlf_nl_local(pcols,pver), dlf_ni_local(pcols,pver)
  real(r8), target, intent(inout) :: dpdlfliq_local(pcols,pver), dpdlfice_local(pcols,pver)
  real(r8), target, intent(inout) :: shdlfliq_local(pcols,pver), shdlfice_local(pcols,pver)
  real(r8), target, intent(inout) :: dpdlft_local(pcols,pver), shdlft_local(pcols,pver)
  real(r8) :: nl_denom_a_local, nl_denom_b_local, ni_denom_a_local, ni_denom_b_local
  logical, save :: detrain_debug_checked = .false.
  logical, save :: detrain_debug_enabled = .false.
  logical, save :: detrain_debug_reported = .false.
  character(len=32) :: debug_name
  integer :: debug_status, debug_len
  real(r8) :: ptend_ql_ref(pcols,pver), ptend_qi_ref(pcols,pver), ptend_nl_ref(pcols,pver), ptend_ni_ref(pcols,pver)
  real(r8) :: ptend_s_ref(pcols,pver), det_s_ref(pcols), det_ice_ref(pcols)
  real(r8) :: dlf_t_ref(pcols,pver), dlf_qv_ref(pcols,pver), dlf_ql_ref(pcols,pver), dlf_qi_ref(pcols,pver)
  real(r8) :: dlf_nl_ref(pcols,pver), dlf_ni_ref(pcols,pver)
  real(r8) :: dpdlfliq_ref(pcols,pver), dpdlfice_ref(pcols,pver), shdlfliq_ref(pcols,pver), shdlfice_ref(pcols,pver)
  real(r8) :: dpdlft_ref(pcols,pver), shdlft_ref(pcols,pver)

  interface
     subroutine macrop_driver_detrain_core_codon(ncol_c, pcols_c, pver_c, top_lev_c, do_detrain_c, cu_det_st_c, cpair_c, &
          gravit_c, latice_c, nl_denom_a_c, nl_denom_b_c, ni_denom_a_c, ni_denom_b_c, state_t_p, state_pdel_p, dlf_p, &
          dlf2_p, ptend_ql_p, ptend_qi_p, ptend_nl_p, ptend_ni_p, ptend_s_p, det_s_p, det_ice_p, dlf_t_p, dlf_qv_p, &
          dlf_ql_p, dlf_qi_p, dlf_nl_p, dlf_ni_p, dpdlfliq_p, &
          dpdlfice_p, shdlfliq_p, shdlfice_p, dpdlft_p, shdlft_p) bind(c, name="macrop_driver_detrain_core_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c, do_detrain_c, cu_det_st_c
       real(c_double), value :: cpair_c, gravit_c, latice_c, nl_denom_a_c, nl_denom_b_c, ni_denom_a_c, ni_denom_b_c
       type(c_ptr), value :: state_t_p, state_pdel_p, dlf_p, dlf2_p, ptend_ql_p, ptend_qi_p, ptend_nl_p, ptend_ni_p
       type(c_ptr), value :: ptend_s_p, det_s_p, det_ice_p, dlf_t_p, dlf_qv_p, dlf_ql_p, dlf_qi_p, dlf_nl_p, dlf_ni_p
       type(c_ptr), value :: dpdlfliq_p, dpdlfice_p, shdlfliq_p, shdlfice_p, dpdlft_p, shdlft_p
     end subroutine macrop_driver_detrain_core_codon
  end interface

  if (use_native_impl) then
     call macrop_driver_detrain_core_native(ncol_local, do_detrain_local, cu_det_st_local, state_t_local, state_pdel_local, &
          dlf_local, dlf2_local, ptend_ql_local, ptend_qi_local, ptend_nl_local, ptend_ni_local, ptend_s_local, det_s_local, &
          det_ice_local, dlf_t_local, dlf_qv_local, dlf_ql_local, dlf_qi_local, dlf_nl_local, dlf_ni_local, dpdlfliq_local, &
          dpdlfice_local, shdlfliq_local, shdlfice_local, dpdlft_local, shdlft_local)
     return
  end if

  if (.not. detrain_debug_checked) then
     debug_name = ''
     call get_environment_variable('MACROP_DRIVER_DETRAIN_DEBUG', value=debug_name, length=debug_len, status=debug_status)
     detrain_debug_enabled = debug_status == 0 .and. debug_len > 0 .and. trim(adjustl(debug_name(:debug_len))) /= '0'
     detrain_debug_checked = .true.
  end if

  nl_denom_a_local = 4._r8*3.14_r8*8.e-6_r8**3*997._r8
  nl_denom_b_local = 4._r8*3.14_r8*10.e-6_r8**3*997._r8
  ni_denom_a_local = 4._r8*3.14_r8*25.e-6_r8**3*500._r8
  ni_denom_b_local = 4._r8*3.14_r8*50.e-6_r8**3*500._r8

  if (detrain_debug_enabled .and. .not. detrain_debug_reported) then
     ptend_ql_ref = ptend_ql_local
     ptend_qi_ref = ptend_qi_local
     ptend_nl_ref = ptend_nl_local
     ptend_ni_ref = ptend_ni_local
     ptend_s_ref = ptend_s_local
     det_s_ref = det_s_local
     det_ice_ref = det_ice_local
     dlf_t_ref = dlf_t_local
     dlf_qv_ref = dlf_qv_local
     dlf_ql_ref = dlf_ql_local
     dlf_qi_ref = dlf_qi_local
     dlf_nl_ref = dlf_nl_local
     dlf_ni_ref = dlf_ni_local
     dpdlfliq_ref = dpdlfliq_local
     dpdlfice_ref = dpdlfice_local
     shdlfliq_ref = shdlfliq_local
     shdlfice_ref = shdlfice_local
     dpdlft_ref = dpdlft_local
     shdlft_ref = shdlft_local

     call macrop_driver_detrain_core_native(ncol_local, do_detrain_local, cu_det_st_local, state_t_local, state_pdel_local, &
          dlf_local, dlf2_local, ptend_ql_ref, ptend_qi_ref, ptend_nl_ref, ptend_ni_ref, ptend_s_ref, det_s_ref, det_ice_ref, &
          dlf_t_ref, dlf_qv_ref, dlf_ql_ref, dlf_qi_ref, dlf_nl_ref, dlf_ni_ref, dpdlfliq_ref, dpdlfice_ref, shdlfliq_ref, &
          shdlfice_ref, dpdlft_ref, shdlft_ref)
  end if

  if (masterproc .and. .not. detrain_core_logged) then
     write(iulog,*) 'macrop_driver detrain core entered = codon'
     call macrop_driver_append_impl_proof('MACROP_DRIVER_DETRAIN_SHELL_PROOF_FILE', &
          'macrop_driver detrain core entered = codon')
     call flush(iulog)
     detrain_core_logged = .true.
  end if

  call macrop_driver_detrain_core_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
       int(top_lev, c_int64_t), merge(1_c_int64_t, 0_c_int64_t, do_detrain_local), &
       merge(1_c_int64_t, 0_c_int64_t, cu_det_st_local), cpair, gravit, latice, nl_denom_a_local, nl_denom_b_local, &
       ni_denom_a_local, ni_denom_b_local, c_loc(state_t_local), c_loc(state_pdel_local), c_loc(dlf_local), c_loc(dlf2_local), &
       c_loc(ptend_ql_local), c_loc(ptend_qi_local), c_loc(ptend_nl_local), c_loc(ptend_ni_local), c_loc(ptend_s_local), &
       c_loc(det_s_local), c_loc(det_ice_local), c_loc(dlf_t_local), c_loc(dlf_qv_local), c_loc(dlf_ql_local), &
       c_loc(dlf_qi_local), c_loc(dlf_nl_local), c_loc(dlf_ni_local), c_loc(dpdlfliq_local), c_loc(dpdlfice_local), &
       c_loc(shdlfliq_local), c_loc(shdlfice_local), c_loc(dpdlft_local), c_loc(shdlft_local))

  if (detrain_debug_enabled .and. .not. detrain_debug_reported) then
     write(iulog,*) 'macrop_driver_detrain_core debug maxdiff ptend_ql =', &
          maxval(abs(ptend_ql_local(1:ncol_local,top_lev:pver) - ptend_ql_ref(1:ncol_local,top_lev:pver)))
     write(iulog,*) 'macrop_driver_detrain_core debug maxdiff ptend_qi =', &
          maxval(abs(ptend_qi_local(1:ncol_local,top_lev:pver) - ptend_qi_ref(1:ncol_local,top_lev:pver)))
     write(iulog,*) 'macrop_driver_detrain_core debug maxdiff ptend_nl =', &
          maxval(abs(ptend_nl_local(1:ncol_local,top_lev:pver) - ptend_nl_ref(1:ncol_local,top_lev:pver)))
     write(iulog,*) 'macrop_driver_detrain_core debug maxdiff ptend_ni =', &
          maxval(abs(ptend_ni_local(1:ncol_local,top_lev:pver) - ptend_ni_ref(1:ncol_local,top_lev:pver)))
     write(iulog,*) 'macrop_driver_detrain_core debug maxdiff ptend_s =', &
          maxval(abs(ptend_s_local(1:ncol_local,top_lev:pver) - ptend_s_ref(1:ncol_local,top_lev:pver)))
     write(iulog,*) 'macrop_driver_detrain_core debug maxdiff det_s =', &
          maxval(abs(det_s_local(1:ncol_local) - det_s_ref(1:ncol_local)))
     write(iulog,*) 'macrop_driver_detrain_core debug maxdiff det_ice =', &
          maxval(abs(det_ice_local(1:ncol_local) - det_ice_ref(1:ncol_local)))
     write(iulog,*) 'macrop_driver_detrain_core debug maxdiff dpdlfliq =', &
          maxval(abs(dpdlfliq_local(1:ncol_local,top_lev:pver) - dpdlfliq_ref(1:ncol_local,top_lev:pver)))
     write(iulog,*) 'macrop_driver_detrain_core debug maxdiff dpdlfice =', &
          maxval(abs(dpdlfice_local(1:ncol_local,top_lev:pver) - dpdlfice_ref(1:ncol_local,top_lev:pver)))
     write(iulog,*) 'macrop_driver_detrain_core debug maxdiff shdlfliq =', &
          maxval(abs(shdlfliq_local(1:ncol_local,top_lev:pver) - shdlfliq_ref(1:ncol_local,top_lev:pver)))
     write(iulog,*) 'macrop_driver_detrain_core debug maxdiff shdlfice =', &
          maxval(abs(shdlfice_local(1:ncol_local,top_lev:pver) - shdlfice_ref(1:ncol_local,top_lev:pver)))
     write(iulog,*) 'macrop_driver_detrain_core debug maxdiff dpdlft =', &
          maxval(abs(dpdlft_local(1:ncol_local,top_lev:pver) - dpdlft_ref(1:ncol_local,top_lev:pver)))
     write(iulog,*) 'macrop_driver_detrain_core debug maxdiff shdlft =', &
          maxval(abs(shdlft_local(1:ncol_local,top_lev:pver) - shdlft_ref(1:ncol_local,top_lev:pver)))
     call flush(iulog)
     detrain_debug_reported = .true.
  end if

end subroutine macrop_driver_detrain_core

!============================================================================ !

subroutine macrop_driver_detrain_core_native(ncol_local, do_detrain_local, cu_det_st_local, state_t_local, state_pdel_local, &
     dlf_local, dlf2_local, ptend_ql_local, ptend_qi_local, ptend_nl_local, ptend_ni_local, ptend_s_local, det_s_local, &
     det_ice_local, dlf_t_local, dlf_qv_local, dlf_ql_local, dlf_qi_local, dlf_nl_local, dlf_ni_local, dpdlfliq_local, &
     dpdlfice_local, shdlfliq_local, shdlfice_local, dpdlft_local, shdlft_local)

  use physconst, only: cpair, gravit
  use ref_pres,  only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local
  logical, intent(in) :: do_detrain_local, cu_det_st_local
  real(r8), intent(in) :: state_t_local(pcols,pver), state_pdel_local(pcols,pver)
  real(r8), intent(in) :: dlf_local(pcols,pver), dlf2_local(pcols,pver)
  real(r8), intent(inout) :: ptend_ql_local(pcols,pver), ptend_qi_local(pcols,pver)
  real(r8), intent(inout) :: ptend_nl_local(pcols,pver), ptend_ni_local(pcols,pver)
  real(r8), intent(inout) :: ptend_s_local(pcols,pver), det_s_local(pcols), det_ice_local(pcols)
  real(r8), intent(inout) :: dlf_t_local(pcols,pver), dlf_qv_local(pcols,pver), dlf_ql_local(pcols,pver)
  real(r8), intent(inout) :: dlf_qi_local(pcols,pver), dlf_nl_local(pcols,pver), dlf_ni_local(pcols,pver)
  real(r8), intent(inout) :: dpdlfliq_local(pcols,pver), dpdlfice_local(pcols,pver)
  real(r8), intent(inout) :: shdlfliq_local(pcols,pver), shdlfice_local(pcols,pver)
  real(r8), intent(inout) :: dpdlft_local(pcols,pver), shdlft_local(pcols,pver)
  integer :: i, k
  real(r8) :: dum1_local

  do k = top_lev, pver
     do i = 1, ncol_local
        if( state_t_local(i,k) > 268.15_r8 ) then
            dum1_local = 0.0_r8
        elseif( state_t_local(i,k) < 238.15_r8 ) then
            dum1_local = 1.0_r8
        else
            dum1_local = ( 268.15_r8 - state_t_local(i,k) ) / 30._r8
        endif

        if (do_detrain_local) then
           ptend_ql_local(i,k) = dlf_local(i,k) * ( 1._r8 - dum1_local )
           ptend_qi_local(i,k) = dlf_local(i,k) * dum1_local
           ptend_nl_local(i,k) = 3._r8 * ( max(0._r8, ( dlf_local(i,k) - dlf2_local(i,k) )) * ( 1._r8 - dum1_local ) ) / &
                (4._r8*3.14_r8* 8.e-6_r8**3*997._r8) + &
                3._r8 * (                         dlf2_local(i,k)    * ( 1._r8 - dum1_local ) ) / &
                (4._r8*3.14_r8*10.e-6_r8**3*997._r8)
           ptend_ni_local(i,k) = 3._r8 * ( max(0._r8, ( dlf_local(i,k) - dlf2_local(i,k) )) *  dum1_local ) / &
                (4._r8*3.14_r8*25.e-6_r8**3*500._r8) + &
                3._r8 * (                         dlf2_local(i,k)    *  dum1_local ) / &
                (4._r8*3.14_r8*50.e-6_r8**3*500._r8)
           ptend_s_local(i,k) = dlf_local(i,k) * dum1_local * latice
        else
           ptend_ql_local(i,k) = 0._r8
           ptend_qi_local(i,k) = 0._r8
           ptend_nl_local(i,k) = 0._r8
           ptend_ni_local(i,k) = 0._r8
           ptend_s_local(i,k) = 0._r8
        end if

        det_s_local(i) = det_s_local(i) + ptend_s_local(i,k) * state_pdel_local(i,k) / gravit
        det_ice_local(i) = det_ice_local(i) - ptend_qi_local(i,k) * state_pdel_local(i,k) / gravit

        if( cu_det_st_local ) then
            dlf_t_local(i,k) = ptend_s_local(i,k)/cpair
            dlf_qv_local(i,k) = 0._r8
            dlf_ql_local(i,k) = ptend_ql_local(i,k)
            dlf_qi_local(i,k) = ptend_qi_local(i,k)
            dlf_nl_local(i,k) = ptend_nl_local(i,k)
            dlf_ni_local(i,k) = ptend_ni_local(i,k)
            ptend_ql_local(i,k) = 0._r8
            ptend_qi_local(i,k) = 0._r8
            ptend_nl_local(i,k) = 0._r8
            ptend_ni_local(i,k) = 0._r8
            ptend_s_local(i,k) = 0._r8
            dpdlfliq_local(i,k) = 0._r8
            dpdlfice_local(i,k) = 0._r8
            shdlfliq_local(i,k) = 0._r8
            shdlfice_local(i,k) = 0._r8
            dpdlft_local(i,k) = 0._r8
            shdlft_local(i,k) = 0._r8
         else
            dpdlfliq_local(i,k) = ( dlf_local(i,k) - dlf2_local(i,k) ) * ( 1._r8 - dum1_local )
            dpdlfice_local(i,k) = ( dlf_local(i,k) - dlf2_local(i,k) ) * ( dum1_local )
            shdlfliq_local(i,k) = dlf2_local(i,k) * ( 1._r8 - dum1_local )
            shdlfice_local(i,k) = dlf2_local(i,k) * ( dum1_local )
            dpdlft_local(i,k) = ( dlf_local(i,k) - dlf2_local(i,k) ) * dum1_local * latice/cpair
            shdlft_local(i,k) = dlf2_local(i,k) * dum1_local * latice/cpair
        endif
     end do
  end do

end subroutine macrop_driver_detrain_core_native

!============================================================================ !

subroutine macrop_driver_detrain_init_shell(ncol_local, dlf_T_local, dlf_qv_local, dlf_ql_local, dlf_qi_local, &
     dlf_nl_local, dlf_ni_local, det_s_local, det_ice_local, dpdlfliq_local, dpdlfice_local, shdlfliq_local, &
     shdlfice_local, dpdlft_local, shdlft_local)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr

  integer, intent(in) :: ncol_local
  real(r8), target, intent(inout) :: dlf_T_local(pcols,pver), dlf_qv_local(pcols,pver), dlf_ql_local(pcols,pver)
  real(r8), target, intent(inout) :: dlf_qi_local(pcols,pver), dlf_nl_local(pcols,pver), dlf_ni_local(pcols,pver)
  real(r8), target, intent(inout) :: det_s_local(pcols), det_ice_local(pcols)
  real(r8), target, intent(inout) :: dpdlfliq_local(pcols,pver), dpdlfice_local(pcols,pver)
  real(r8), target, intent(inout) :: shdlfliq_local(pcols,pver), shdlfice_local(pcols,pver)
  real(r8), target, intent(inout) :: dpdlft_local(pcols,pver), shdlft_local(pcols,pver)

  interface
     subroutine macrop_driver_detrain_init_shell_codon(ncol_c, pcols_c, pver_c, dlf_T_p, dlf_qv_p, dlf_ql_p, &
          dlf_qi_p, dlf_nl_p, dlf_ni_p, det_s_p, det_ice_p, dpdlfliq_p, dpdlfice_p, shdlfliq_p, shdlfice_p, &
          dpdlft_p, shdlft_p) bind(c, name="macrop_driver_detrain_init_shell_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
       type(c_ptr), value :: dlf_T_p, dlf_qv_p, dlf_ql_p, dlf_qi_p, dlf_nl_p, dlf_ni_p, det_s_p, det_ice_p
       type(c_ptr), value :: dpdlfliq_p, dpdlfice_p, shdlfliq_p, shdlfice_p, dpdlft_p, shdlft_p
     end subroutine macrop_driver_detrain_init_shell_codon
  end interface

  if (use_native_impl) then
     call macrop_driver_detrain_init_shell_native(ncol_local, dlf_T_local, dlf_qv_local, dlf_ql_local, dlf_qi_local, &
          dlf_nl_local, dlf_ni_local, det_s_local, det_ice_local, dpdlfliq_local, dpdlfice_local, shdlfliq_local, &
          shdlfice_local, dpdlft_local, shdlft_local)
     return
  end if

  if (masterproc .and. .not. detrain_init_shell_logged) then
     write(iulog,*) 'macrop_driver detrain init shell entered = codon'
     call macrop_driver_append_impl_proof('MACROP_DRIVER_DETRAIN_SHELL_PROOF_FILE', &
          'macrop_driver detrain init shell entered = codon')
     call flush(iulog)
     detrain_init_shell_logged = .true.
  end if

  call macrop_driver_detrain_init_shell_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
       c_loc(dlf_T_local), c_loc(dlf_qv_local), c_loc(dlf_ql_local), c_loc(dlf_qi_local), c_loc(dlf_nl_local), &
       c_loc(dlf_ni_local), c_loc(det_s_local), c_loc(det_ice_local), c_loc(dpdlfliq_local), c_loc(dpdlfice_local), &
       c_loc(shdlfliq_local), c_loc(shdlfice_local), c_loc(dpdlft_local), c_loc(shdlft_local))

end subroutine macrop_driver_detrain_init_shell

!============================================================================ !

subroutine macrop_driver_detrain_init_shell_native(ncol_local, dlf_T_local, dlf_qv_local, dlf_ql_local, dlf_qi_local, &
     dlf_nl_local, dlf_ni_local, det_s_local, det_ice_local, dpdlfliq_local, dpdlfice_local, shdlfliq_local, &
     shdlfice_local, dpdlft_local, shdlft_local)

  integer, intent(in) :: ncol_local
  real(r8), intent(inout) :: dlf_T_local(pcols,pver), dlf_qv_local(pcols,pver), dlf_ql_local(pcols,pver)
  real(r8), intent(inout) :: dlf_qi_local(pcols,pver), dlf_nl_local(pcols,pver), dlf_ni_local(pcols,pver)
  real(r8), intent(inout) :: det_s_local(pcols), det_ice_local(pcols)
  real(r8), intent(inout) :: dpdlfliq_local(pcols,pver), dpdlfice_local(pcols,pver)
  real(r8), intent(inout) :: shdlfliq_local(pcols,pver), shdlfice_local(pcols,pver)
  real(r8), intent(inout) :: dpdlft_local(pcols,pver), shdlft_local(pcols,pver)

  dlf_T_local(:,:) = 0._r8
  dlf_qv_local(:,:) = 0._r8
  dlf_ql_local(:,:) = 0._r8
  dlf_qi_local(:,:) = 0._r8
  dlf_nl_local(:,:) = 0._r8
  dlf_ni_local(:,:) = 0._r8
  det_s_local(:) = 0._r8
  det_ice_local(:) = 0._r8
  dpdlfliq_local(:,:) = 0._r8
  dpdlfice_local(:,:) = 0._r8
  shdlfliq_local(:,:) = 0._r8
  shdlfice_local(:,:) = 0._r8
  dpdlft_local(:,:) = 0._r8
  shdlft_local(:,:) = 0._r8

end subroutine macrop_driver_detrain_init_shell_native

!============================================================================ !

subroutine macrop_driver_detrain_post_shell(ncol_local, det_ice_local)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr

  integer, intent(in) :: ncol_local
  real(r8), target, intent(inout) :: det_ice_local(pcols)

  interface
     subroutine macrop_driver_detrain_post_shell_codon(ncol_c, pcols_c, det_ice_p) &
          bind(c, name="macrop_driver_detrain_post_shell_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c
       type(c_ptr), value :: det_ice_p
     end subroutine macrop_driver_detrain_post_shell_codon
  end interface

  if (use_native_impl) then
     call macrop_driver_detrain_post_shell_native(ncol_local, det_ice_local)
     return
  end if

  if (masterproc .and. .not. detrain_post_shell_logged) then
     write(iulog,*) 'macrop_driver detrain post shell entered = codon'
     call macrop_driver_append_impl_proof('MACROP_DRIVER_DETRAIN_SHELL_PROOF_FILE', &
          'macrop_driver detrain post shell entered = codon')
     call flush(iulog)
     detrain_post_shell_logged = .true.
  end if

  call macrop_driver_detrain_post_shell_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), c_loc(det_ice_local))

end subroutine macrop_driver_detrain_post_shell

!============================================================================ !

subroutine macrop_driver_detrain_post_shell_native(ncol_local, det_ice_local)

  integer, intent(in) :: ncol_local
  real(r8), intent(inout) :: det_ice_local(pcols)

  det_ice_local(:ncol_local) = det_ice_local(:ncol_local) / 1000._r8

end subroutine macrop_driver_detrain_post_shell_native

!============================================================================ !

subroutine macrop_driver_mmacro_input_shell(ncol_local, state_q_local, zeros_local, qc_local, qi_local, nc_local, ni_local)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use constituents, only: pcnst
  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local
  real(r8), target, intent(in) :: state_q_local(pcols,pver,pcnst)
  real(r8), target, intent(inout) :: zeros_local(pcols,pver), qc_local(pcols,pver), qi_local(pcols,pver)
  real(r8), target, intent(inout) :: nc_local(pcols,pver), ni_local(pcols,pver)

  interface
     subroutine macrop_driver_mmacro_input_shell_codon(ncol_c, pcols_c, pver_c, pcnst_c, top_lev_c, ixcldliq_c, &
          ixcldice_c, ixnumliq_c, ixnumice_c, state_q_p, zeros_p, qc_p, qi_p, nc_p, ni_p) &
          bind(c, name="macrop_driver_mmacro_input_shell_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pcnst_c, top_lev_c
       integer(c_int64_t), value :: ixcldliq_c, ixcldice_c, ixnumliq_c, ixnumice_c
       type(c_ptr), value :: state_q_p, zeros_p, qc_p, qi_p, nc_p, ni_p
     end subroutine macrop_driver_mmacro_input_shell_codon
  end interface

  if (use_native_impl) then
     call macrop_driver_mmacro_input_shell_native(ncol_local, state_q_local, zeros_local, qc_local, qi_local, nc_local, ni_local)
     return
  end if

  if (masterproc .and. .not. mmacro_input_shell_logged) then
     write(iulog,*) 'macrop_driver mmacro input shell entered = codon'
     call macrop_driver_append_impl_proof('MACROP_DRIVER_MMPCOND_SHELL_PROOF_FILE', &
          'macrop_driver mmacro input shell entered = codon')
     call flush(iulog)
     mmacro_input_shell_logged = .true.
  end if

  call macrop_driver_mmacro_input_shell_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), &
       int(pver, c_int64_t), int(pcnst, c_int64_t), int(top_lev, c_int64_t), int(ixcldliq, c_int64_t), &
       int(ixcldice, c_int64_t), int(ixnumliq, c_int64_t), int(ixnumice, c_int64_t), c_loc(state_q_local), &
       c_loc(zeros_local), c_loc(qc_local), c_loc(qi_local), c_loc(nc_local), c_loc(ni_local))

end subroutine macrop_driver_mmacro_input_shell

!============================================================================ !

subroutine macrop_driver_mmacro_input_shell_native(ncol_local, state_q_local, zeros_local, qc_local, qi_local, nc_local, ni_local)

  use constituents, only: pcnst
  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local
  real(r8), intent(in) :: state_q_local(pcols,pver,pcnst)
  real(r8), intent(inout) :: zeros_local(pcols,pver), qc_local(pcols,pver), qi_local(pcols,pver)
  real(r8), intent(inout) :: nc_local(pcols,pver), ni_local(pcols,pver)

  zeros_local(:ncol_local,top_lev:pver) = 0._r8
  qc_local(:ncol_local,top_lev:pver) = state_q_local(:ncol_local,top_lev:pver,ixcldliq)
  qi_local(:ncol_local,top_lev:pver) = state_q_local(:ncol_local,top_lev:pver,ixcldice)
  nc_local(:ncol_local,top_lev:pver) = state_q_local(:ncol_local,top_lev:pver,ixnumliq)
  ni_local(:ncol_local,top_lev:pver) = state_q_local(:ncol_local,top_lev:pver,ixnumice)

end subroutine macrop_driver_mmacro_input_shell_native

!============================================================================ !

subroutine macrop_driver_mmacro_post_fields_shell(ncol_local, fice_local, alst_local, aist_local, fice_ql_local, ast_local)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local
  real(r8), target, intent(in) :: fice_local(pcols,pver), alst_local(pcols,pver), aist_local(pcols,pver)
  real(r8), target, intent(inout) :: fice_ql_local(pcols,pver), ast_local(pcols,pver)

  interface
     subroutine macrop_driver_mmacro_post_fields_shell_codon(ncol_c, pcols_c, pver_c, top_lev_c, fice_p, &
          alst_p, aist_p, fice_ql_p, ast_p) bind(c, name="macrop_driver_mmacro_post_fields_shell_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c
       type(c_ptr), value :: fice_p, alst_p, aist_p, fice_ql_p, ast_p
     end subroutine macrop_driver_mmacro_post_fields_shell_codon
  end interface

  if (use_native_impl) then
     call macrop_driver_mmacro_post_fields_shell_native(ncol_local, fice_local, alst_local, aist_local, fice_ql_local, ast_local)
     return
  end if

  if (masterproc .and. .not. mmacro_post_fields_shell_logged) then
     write(iulog,*) 'macrop_driver mmacro post fields shell entered = codon'
     call macrop_driver_append_impl_proof('MACROP_DRIVER_MMPCOND_SHELL_PROOF_FILE', &
          'macrop_driver mmacro post fields shell entered = codon')
     call flush(iulog)
     mmacro_post_fields_shell_logged = .true.
  end if

  call macrop_driver_mmacro_post_fields_shell_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), &
       int(pver, c_int64_t), int(top_lev, c_int64_t), c_loc(fice_local), c_loc(alst_local), c_loc(aist_local), &
       c_loc(fice_ql_local), c_loc(ast_local))

end subroutine macrop_driver_mmacro_post_fields_shell

!============================================================================ !

subroutine macrop_driver_mmacro_post_fields_shell_native(ncol_local, fice_local, alst_local, aist_local, fice_ql_local, ast_local)

  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local
  real(r8), intent(in) :: fice_local(pcols,pver), alst_local(pcols,pver), aist_local(pcols,pver)
  real(r8), intent(inout) :: fice_ql_local(pcols,pver), ast_local(pcols,pver)

  fice_ql_local(:ncol_local,:top_lev-1) = 0._r8
  fice_ql_local(:ncol_local,top_lev:pver) = fice_local(:ncol_local,top_lev:pver)
  ast_local(:ncol_local,:top_lev-1) = 0._r8
  ast_local(:ncol_local,top_lev:pver) = max(alst_local(:ncol_local,top_lev:pver), aist_local(:ncol_local,top_lev:pver))

end subroutine macrop_driver_mmacro_post_fields_shell_native

!============================================================================ !

integer function macrop_driver_mmacro_config_check(ncol_local, do_cldice_local, do_cldliq_local, qiten_local, niten_local, &
     qcten_local, ncten_local) result(mask)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local
  logical, intent(in) :: do_cldice_local, do_cldliq_local
  real(r8), target, intent(in) :: qiten_local(pcols,pver), niten_local(pcols,pver)
  real(r8), target, intent(in) :: qcten_local(pcols,pver), ncten_local(pcols,pver)
  integer(c_int64_t), target :: mask64

  interface
     subroutine macrop_driver_mmacro_config_check_codon(ncol_c, pcols_c, pver_c, top_lev_c, do_cldice_c, &
          do_cldliq_c, qiten_p, niten_p, qcten_p, ncten_p, mask_p) &
          bind(c, name="macrop_driver_mmacro_config_check_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c, do_cldice_c, do_cldliq_c
       type(c_ptr), value :: qiten_p, niten_p, qcten_p, ncten_p, mask_p
     end subroutine macrop_driver_mmacro_config_check_codon
  end interface

  if (use_native_impl) then
     mask = macrop_driver_mmacro_config_check_native(ncol_local, do_cldice_local, do_cldliq_local, qiten_local, niten_local, &
          qcten_local, ncten_local)
     return
  end if

  if (masterproc .and. .not. mmacro_config_check_logged) then
     write(iulog,*) 'macrop_driver mmacro config check entered = codon'
     call macrop_driver_append_impl_proof('MACROP_DRIVER_MMPCOND_SHELL_PROOF_FILE', &
          'macrop_driver mmacro config check entered = codon')
     call flush(iulog)
     mmacro_config_check_logged = .true.
  end if

  mask64 = 0_c_int64_t
  call macrop_driver_mmacro_config_check_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
       int(top_lev, c_int64_t), merge(1_c_int64_t, 0_c_int64_t, do_cldice_local), &
       merge(1_c_int64_t, 0_c_int64_t, do_cldliq_local), c_loc(qiten_local), c_loc(niten_local), c_loc(qcten_local), &
       c_loc(ncten_local), c_loc(mask64))
  mask = int(mask64)

end function macrop_driver_mmacro_config_check

!============================================================================ !

integer function macrop_driver_mmacro_config_check_native(ncol_local, do_cldice_local, do_cldliq_local, qiten_local, niten_local, &
     qcten_local, ncten_local) result(mask)

  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local
  logical, intent(in) :: do_cldice_local, do_cldliq_local
  real(r8), intent(in) :: qiten_local(pcols,pver), niten_local(pcols,pver)
  real(r8), intent(in) :: qcten_local(pcols,pver), ncten_local(pcols,pver)

  mask = 0
  if ((.not. do_cldice_local) .and. any(qiten_local(:ncol_local,top_lev:pver) /= 0.0_r8)) mask = ior(mask, 1)
  if ((.not. do_cldice_local) .and. any(niten_local(:ncol_local,top_lev:pver) /= 0.0_r8)) mask = ior(mask, 2)
  if ((.not. do_cldliq_local) .and. any(qcten_local(:ncol_local,top_lev:pver) /= 0.0_r8)) mask = ior(mask, 4)
  if ((.not. do_cldliq_local) .and. any(ncten_local(:ncol_local,top_lev:pver) /= 0.0_r8)) mask = ior(mask, 8)

end function macrop_driver_mmacro_config_check_native

!============================================================================ !

subroutine macrop_driver_cfmip_diag_shell(ncol_local, cld_local, state_ql_local, state_qi_local, mr_ccliq_local, &
     mr_ccice_local, mr_lsliq_local, mr_lsice_local)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local
  real(r8), target, intent(in) :: cld_local(pcols,pver), state_ql_local(pcols,pver), state_qi_local(pcols,pver)
  real(r8), target, intent(inout) :: mr_ccliq_local(pcols,pver), mr_ccice_local(pcols,pver)
  real(r8), target, intent(inout) :: mr_lsliq_local(pcols,pver), mr_lsice_local(pcols,pver)

  interface
     subroutine macrop_driver_cfmip_diag_shell_codon(ncol_c, pcols_c, pver_c, top_lev_c, cld_p, state_ql_p, &
          state_qi_p, mr_ccliq_p, mr_ccice_p, mr_lsliq_p, mr_lsice_p) bind(c, name="macrop_driver_cfmip_diag_shell_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c
       type(c_ptr), value :: cld_p, state_ql_p, state_qi_p, mr_ccliq_p, mr_ccice_p, mr_lsliq_p, mr_lsice_p
     end subroutine macrop_driver_cfmip_diag_shell_codon
  end interface

  if (use_native_impl) then
     call macrop_driver_cfmip_diag_shell_native(ncol_local, cld_local, state_ql_local, state_qi_local, mr_ccliq_local, &
          mr_ccice_local, mr_lsliq_local, mr_lsice_local)
     return
  end if

  if (masterproc .and. .not. cfmip_diag_shell_logged) then
     write(iulog,*) 'macrop_driver cfmip diag shell entered = codon'
     call macrop_driver_append_impl_proof('MACROP_DRIVER_CFMIP_DIAG_SHELL_PROOF_FILE', &
          'macrop_driver cfmip diag shell entered = codon')
     call flush(iulog)
     cfmip_diag_shell_logged = .true.
  end if

  call macrop_driver_cfmip_diag_shell_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
       int(top_lev, c_int64_t), c_loc(cld_local), c_loc(state_ql_local), c_loc(state_qi_local), c_loc(mr_ccliq_local), &
       c_loc(mr_ccice_local), c_loc(mr_lsliq_local), c_loc(mr_lsice_local))

end subroutine macrop_driver_cfmip_diag_shell

!============================================================================ !

subroutine macrop_driver_cfmip_diag_shell_native(ncol_local, cld_local, state_ql_local, state_qi_local, mr_ccliq_local, &
     mr_ccice_local, mr_lsliq_local, mr_lsice_local)

  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local
  real(r8), intent(in) :: cld_local(pcols,pver), state_ql_local(pcols,pver), state_qi_local(pcols,pver)
  real(r8), intent(inout) :: mr_ccliq_local(pcols,pver), mr_ccice_local(pcols,pver)
  real(r8), intent(inout) :: mr_lsliq_local(pcols,pver), mr_lsice_local(pcols,pver)
  integer :: i, k

  mr_ccliq_local = 0._r8
  mr_ccice_local = 0._r8
  mr_lsliq_local = 0._r8
  mr_lsice_local = 0._r8

  do k = top_lev, pver
     do i = 1, ncol_local
        if (cld_local(i,k) .gt. 0._r8) then
           mr_lsliq_local(i,k) = state_ql_local(i,k)
           mr_lsice_local(i,k) = state_qi_local(i,k)
        else
           mr_lsliq_local(i,k) = 0._r8
           mr_lsice_local(i,k) = 0._r8
        end if
     end do
  end do

end subroutine macrop_driver_cfmip_diag_shell_native

!============================================================================ !

subroutine macrop_driver_clr_old_diag(ncol_local, concld_local, alst_local, ast_local, clrw_old_local, clri_old_local)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local
  real(r8), target, intent(in) :: concld_local(pcols,pver), alst_local(pcols,pver), ast_local(pcols,pver)
  real(r8), target, intent(inout) :: clrw_old_local(pcols,pver), clri_old_local(pcols,pver)

  interface
     subroutine macrop_driver_clr_old_diag_codon(ncol_c, pcols_c, pver_c, top_lev_c, concld_p, alst_p, ast_p, &
          clrw_old_p, clri_old_p) bind(c, name="macrop_driver_clr_old_diag_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c
       type(c_ptr), value :: concld_p, alst_p, ast_p, clrw_old_p, clri_old_p
     end subroutine macrop_driver_clr_old_diag_codon
  end interface

  if (use_native_impl) then
     call macrop_driver_clr_old_diag_native(ncol_local, concld_local, alst_local, ast_local, clrw_old_local, clri_old_local)
     return
  end if

  if (masterproc .and. .not. clr_old_diag_logged) then
     write(iulog,*) 'macrop_driver clr old diag entered = codon'
     call macrop_driver_append_impl_proof('MACROP_DRIVER_MMPCOND_SHELL_PROOF_FILE', &
          'macrop_driver clr old diag entered = codon')
     call flush(iulog)
     clr_old_diag_logged = .true.
  end if

  call macrop_driver_clr_old_diag_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
       int(top_lev, c_int64_t), c_loc(concld_local), c_loc(alst_local), c_loc(ast_local), c_loc(clrw_old_local), &
       c_loc(clri_old_local))

end subroutine macrop_driver_clr_old_diag

!============================================================================ !

subroutine macrop_driver_clr_old_diag_native(ncol_local, concld_local, alst_local, ast_local, clrw_old_local, clri_old_local)

  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local
  real(r8), intent(in) :: concld_local(pcols,pver), alst_local(pcols,pver), ast_local(pcols,pver)
  real(r8), intent(inout) :: clrw_old_local(pcols,pver), clri_old_local(pcols,pver)
  integer :: i, k

  do k = top_lev, pver
     do i = 1, ncol_local
        clrw_old_local(i,k) = max( 0._r8, min( 1._r8, 1._r8 - concld_local(i,k) - alst_local(i,k) ) )
        clri_old_local(i,k) = max( 0._r8, min( 1._r8, 1._r8 - concld_local(i,k) -  ast_local(i,k) ) )
     end do
  end do

end subroutine macrop_driver_clr_old_diag_native

!============================================================================ !

subroutine macrop_driver_forcing_prep(ncol_local, nstep_local, rdtime_local, state_t_local, state_qv_local, qc_local, qi_local, &
     nc_local, ni_local, tcwat_local, qcwat_local, lcwat_local, iccwat_local, nlwat_local, niwat_local, cc_t_local, &
     cc_qv_local, cc_ql_local, cc_qi_local, cc_nl_local, cc_ni_local, cc_qlst_local, ttend_local, qtend_local, ltend_local, &
     itend_local, nltend_local, nitend_local, lmitend_local, t_inout_local, qv_inout_local, ql_inout_local, qi_inout_local, &
     nl_inout_local, ni_inout_local)

  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local, nstep_local
  real(r8), intent(in) :: rdtime_local
  real(r8), target, intent(in) :: state_t_local(pcols,pver), state_qv_local(pcols,pver)
  real(r8), target, intent(in) :: qc_local(pcols,pver), qi_local(pcols,pver), nc_local(pcols,pver), ni_local(pcols,pver)
  real(r8), target, intent(inout) :: tcwat_local(pcols,pver), qcwat_local(pcols,pver), lcwat_local(pcols,pver)
  real(r8), target, intent(inout) :: iccwat_local(pcols,pver), nlwat_local(pcols,pver), niwat_local(pcols,pver)
  real(r8), target, intent(inout) :: cc_t_local(pcols,pver), cc_qv_local(pcols,pver), cc_ql_local(pcols,pver)
  real(r8), target, intent(inout) :: cc_qi_local(pcols,pver), cc_nl_local(pcols,pver), cc_ni_local(pcols,pver)
  real(r8), target, intent(inout) :: cc_qlst_local(pcols,pver)
  real(r8), target, intent(inout) :: ttend_local(pcols,pver), qtend_local(pcols,pver), ltend_local(pcols,pver)
  real(r8), target, intent(inout) :: itend_local(pcols,pver), nltend_local(pcols,pver), nitend_local(pcols,pver)
  real(r8), target, intent(inout) :: lmitend_local(pcols,pver), t_inout_local(pcols,pver), qv_inout_local(pcols,pver)
  real(r8), target, intent(inout) :: ql_inout_local(pcols,pver), qi_inout_local(pcols,pver)
  real(r8), target, intent(inout) :: nl_inout_local(pcols,pver), ni_inout_local(pcols,pver)

  interface
     subroutine macrop_driver_forcing_prep_codon(ncol_c, pcols_c, pver_c, top_lev_c, nstep_c, rdtime_c, state_t_p, &
          state_qv_p, qc_p, qi_p, nc_p, ni_p, tcwat_p, qcwat_p, lcwat_p, iccwat_p, nlwat_p, niwat_p, cc_t_p, cc_qv_p, &
          cc_ql_p, cc_qi_p, cc_nl_p, cc_ni_p, cc_qlst_p, ttend_p, qtend_p, ltend_p, itend_p, nltend_p, nitend_p, &
          lmitend_p, t_inout_p, qv_inout_p, ql_inout_p, qi_inout_p, nl_inout_p, ni_inout_p) bind(c, name="macrop_driver_forcing_prep_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c, nstep_c
       real(c_double), value :: rdtime_c
       type(c_ptr), value :: state_t_p, state_qv_p, qc_p, qi_p, nc_p, ni_p, tcwat_p, qcwat_p, lcwat_p, iccwat_p
       type(c_ptr), value :: nlwat_p, niwat_p, cc_t_p, cc_qv_p, cc_ql_p, cc_qi_p, cc_nl_p, cc_ni_p, cc_qlst_p
       type(c_ptr), value :: ttend_p, qtend_p, ltend_p, itend_p, nltend_p, nitend_p, lmitend_p, t_inout_p
       type(c_ptr), value :: qv_inout_p, ql_inout_p, qi_inout_p, nl_inout_p, ni_inout_p
     end subroutine macrop_driver_forcing_prep_codon
  end interface

  if (use_native_impl) then
     call macrop_driver_forcing_prep_native(ncol_local, nstep_local, rdtime_local, state_t_local, state_qv_local, qc_local, &
          qi_local, nc_local, ni_local, tcwat_local, qcwat_local, lcwat_local, iccwat_local, nlwat_local, niwat_local, &
          cc_t_local, cc_qv_local, cc_ql_local, cc_qi_local, cc_nl_local, cc_ni_local, cc_qlst_local, ttend_local, &
          qtend_local, ltend_local, itend_local, nltend_local, nitend_local, lmitend_local, t_inout_local, qv_inout_local, &
          ql_inout_local, qi_inout_local, nl_inout_local, ni_inout_local)
     return
  end if

  if (masterproc .and. .not. forcing_prep_logged) then
     write(iulog,*) 'macrop_driver forcing prep entered = codon'
     call macrop_driver_append_impl_proof('MACROP_DRIVER_MMPCOND_SHELL_PROOF_FILE', &
          'macrop_driver forcing prep entered = codon')
     call flush(iulog)
     forcing_prep_logged = .true.
  end if

  call macrop_driver_forcing_prep_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
       int(top_lev, c_int64_t), int(nstep_local, c_int64_t), rdtime_local, c_loc(state_t_local), c_loc(state_qv_local), &
       c_loc(qc_local), c_loc(qi_local), c_loc(nc_local), c_loc(ni_local), c_loc(tcwat_local), c_loc(qcwat_local), &
       c_loc(lcwat_local), c_loc(iccwat_local), c_loc(nlwat_local), c_loc(niwat_local), c_loc(cc_t_local), c_loc(cc_qv_local), &
       c_loc(cc_ql_local), c_loc(cc_qi_local), c_loc(cc_nl_local), c_loc(cc_ni_local), c_loc(cc_qlst_local), c_loc(ttend_local), &
       c_loc(qtend_local), c_loc(ltend_local), c_loc(itend_local), c_loc(nltend_local), c_loc(nitend_local), c_loc(lmitend_local), &
       c_loc(t_inout_local), c_loc(qv_inout_local), c_loc(ql_inout_local), c_loc(qi_inout_local), c_loc(nl_inout_local), &
       c_loc(ni_inout_local))

end subroutine macrop_driver_forcing_prep

!============================================================================ !

subroutine macrop_driver_forcing_prep_native(ncol_local, nstep_local, rdtime_local, state_t_local, state_qv_local, qc_local, &
     qi_local, nc_local, ni_local, tcwat_local, qcwat_local, lcwat_local, iccwat_local, nlwat_local, niwat_local, cc_t_local, &
     cc_qv_local, cc_ql_local, cc_qi_local, cc_nl_local, cc_ni_local, cc_qlst_local, ttend_local, qtend_local, ltend_local, &
     itend_local, nltend_local, nitend_local, lmitend_local, t_inout_local, qv_inout_local, ql_inout_local, qi_inout_local, &
     nl_inout_local, ni_inout_local)

  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local, nstep_local
  real(r8), intent(in) :: rdtime_local
  real(r8), intent(in) :: state_t_local(pcols,pver), state_qv_local(pcols,pver)
  real(r8), intent(in) :: qc_local(pcols,pver), qi_local(pcols,pver), nc_local(pcols,pver), ni_local(pcols,pver)
  real(r8), intent(inout) :: tcwat_local(pcols,pver), qcwat_local(pcols,pver), lcwat_local(pcols,pver)
  real(r8), intent(inout) :: iccwat_local(pcols,pver), nlwat_local(pcols,pver), niwat_local(pcols,pver)
  real(r8), intent(inout) :: cc_t_local(pcols,pver), cc_qv_local(pcols,pver), cc_ql_local(pcols,pver), cc_qi_local(pcols,pver)
  real(r8), intent(inout) :: cc_nl_local(pcols,pver), cc_ni_local(pcols,pver), cc_qlst_local(pcols,pver)
  real(r8), intent(inout) :: ttend_local(pcols,pver), qtend_local(pcols,pver), ltend_local(pcols,pver), itend_local(pcols,pver)
  real(r8), intent(inout) :: nltend_local(pcols,pver), nitend_local(pcols,pver), lmitend_local(pcols,pver)
  real(r8), intent(inout) :: t_inout_local(pcols,pver), qv_inout_local(pcols,pver), ql_inout_local(pcols,pver)
  real(r8), intent(inout) :: qi_inout_local(pcols,pver), nl_inout_local(pcols,pver), ni_inout_local(pcols,pver)
  integer :: i, k

  if( nstep_local .le. 1 ) then
     tcwat_local(:ncol_local,:) = state_t_local(:ncol_local,:)
     qcwat_local(:ncol_local,:) = state_qv_local(:ncol_local,:)
     lcwat_local(:ncol_local,:) = qc_local(:ncol_local,:) + qi_local(:ncol_local,:)
     iccwat_local(:ncol_local,:) = qi_local(:ncol_local,:)
     nlwat_local(:ncol_local,:) = nc_local(:ncol_local,:)
     niwat_local(:ncol_local,:) = ni_local(:ncol_local,:)
     ttend_local(:ncol_local,:) = 0._r8
     qtend_local(:ncol_local,:) = 0._r8
     ltend_local(:ncol_local,:) = 0._r8
     itend_local(:ncol_local,:) = 0._r8
     nltend_local(:ncol_local,:) = 0._r8
     nitend_local(:ncol_local,:) = 0._r8
     cc_t_local(:ncol_local,:) = 0._r8
     cc_qv_local(:ncol_local,:) = 0._r8
     cc_ql_local(:ncol_local,:) = 0._r8
     cc_qi_local(:ncol_local,:) = 0._r8
     cc_nl_local(:ncol_local,:) = 0._r8
     cc_ni_local(:ncol_local,:) = 0._r8
     cc_qlst_local(:ncol_local,:) = 0._r8
  else
     ttend_local(:ncol_local,top_lev:pver) = ( state_t_local(:ncol_local,top_lev:pver) - tcwat_local(:ncol_local,top_lev:pver) ) * rdtime_local &
          - cc_t_local(:ncol_local,top_lev:pver)
     qtend_local(:ncol_local,top_lev:pver) = ( state_qv_local(:ncol_local,top_lev:pver) - qcwat_local(:ncol_local,top_lev:pver) ) * rdtime_local &
          - cc_qv_local(:ncol_local,top_lev:pver)
     ltend_local(:ncol_local,top_lev:pver) = ( qc_local(:ncol_local,top_lev:pver) + qi_local(:ncol_local,top_lev:pver) - lcwat_local(:ncol_local,top_lev:pver) ) * rdtime_local &
          - (cc_ql_local(:ncol_local,top_lev:pver) + cc_qi_local(:ncol_local,top_lev:pver))
     itend_local(:ncol_local,top_lev:pver) = ( qi_local(:ncol_local,top_lev:pver) - iccwat_local(:ncol_local,top_lev:pver) ) * rdtime_local &
          - cc_qi_local(:ncol_local,top_lev:pver)
     nltend_local(:ncol_local,top_lev:pver) = ( nc_local(:ncol_local,top_lev:pver) - nlwat_local(:ncol_local,top_lev:pver) ) * rdtime_local &
          - cc_nl_local(:ncol_local,top_lev:pver)
     nitend_local(:ncol_local,top_lev:pver) = ( ni_local(:ncol_local,top_lev:pver) - niwat_local(:ncol_local,top_lev:pver) ) * rdtime_local &
          - cc_ni_local(:ncol_local,top_lev:pver)
  endif

  do k = top_lev, pver
     do i = 1, ncol_local
        lmitend_local(i,k) = ltend_local(i,k) - itend_local(i,k)
        t_inout_local(i,k) = tcwat_local(i,k)
        qv_inout_local(i,k) = qcwat_local(i,k)
        ql_inout_local(i,k) = lcwat_local(i,k) - iccwat_local(i,k)
        qi_inout_local(i,k) = iccwat_local(i,k)
        nl_inout_local(i,k) = nlwat_local(i,k)
        ni_inout_local(i,k) = niwat_local(i,k)
     end do
  end do

end subroutine macrop_driver_forcing_prep_native

!============================================================================ !

subroutine macrop_driver_ptend_assign(ncol_local, tlat_local, qvlat_local, qcten_local, qiten_local, ncten_local, niten_local, &
     ptend_s_local, ptend_qv_local, ptend_ql_local, ptend_qi_local, ptend_nl_local, ptend_ni_local)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local
  real(r8), target, intent(in) :: tlat_local(pcols,pver), qvlat_local(pcols,pver), qcten_local(pcols,pver)
  real(r8), target, intent(in) :: qiten_local(pcols,pver), ncten_local(pcols,pver), niten_local(pcols,pver)
  real(r8), target, intent(inout) :: ptend_s_local(pcols,pver), ptend_qv_local(pcols,pver), ptend_ql_local(pcols,pver)
  real(r8), target, intent(inout) :: ptend_qi_local(pcols,pver), ptend_nl_local(pcols,pver), ptend_ni_local(pcols,pver)

  interface
     subroutine macrop_driver_ptend_assign_codon(ncol_c, pcols_c, pver_c, top_lev_c, tlat_p, qvlat_p, qcten_p, qiten_p, &
          ncten_p, niten_p, ptend_s_p, ptend_qv_p, ptend_ql_p, ptend_qi_p, ptend_nl_p, ptend_ni_p) bind(c, name="macrop_driver_ptend_assign_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c
       type(c_ptr), value :: tlat_p, qvlat_p, qcten_p, qiten_p, ncten_p, niten_p, ptend_s_p, ptend_qv_p
       type(c_ptr), value :: ptend_ql_p, ptend_qi_p, ptend_nl_p, ptend_ni_p
     end subroutine macrop_driver_ptend_assign_codon
  end interface

  if (use_native_impl) then
     call macrop_driver_ptend_assign_native(ncol_local, tlat_local, qvlat_local, qcten_local, qiten_local, ncten_local, &
          niten_local, ptend_s_local, ptend_qv_local, ptend_ql_local, ptend_qi_local, ptend_nl_local, ptend_ni_local)
     return
  end if

  if (masterproc .and. .not. ptend_assign_logged) then
     write(iulog,*) 'macrop_driver ptend assign entered = codon'
     call macrop_driver_append_impl_proof('MACROP_DRIVER_MMPCOND_SHELL_PROOF_FILE', &
          'macrop_driver ptend assign entered = codon')
     call flush(iulog)
     ptend_assign_logged = .true.
  end if

  call macrop_driver_ptend_assign_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
       int(top_lev, c_int64_t), c_loc(tlat_local), c_loc(qvlat_local), c_loc(qcten_local), c_loc(qiten_local), &
       c_loc(ncten_local), c_loc(niten_local), c_loc(ptend_s_local), c_loc(ptend_qv_local), c_loc(ptend_ql_local), &
       c_loc(ptend_qi_local), c_loc(ptend_nl_local), c_loc(ptend_ni_local))

end subroutine macrop_driver_ptend_assign

!============================================================================ !

subroutine macrop_driver_ptend_assign_native(ncol_local, tlat_local, qvlat_local, qcten_local, qiten_local, ncten_local, &
     niten_local, ptend_s_local, ptend_qv_local, ptend_ql_local, ptend_qi_local, ptend_nl_local, ptend_ni_local)

  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local
  real(r8), intent(in) :: tlat_local(pcols,pver), qvlat_local(pcols,pver), qcten_local(pcols,pver), qiten_local(pcols,pver)
  real(r8), intent(in) :: ncten_local(pcols,pver), niten_local(pcols,pver)
  real(r8), intent(inout) :: ptend_s_local(pcols,pver), ptend_qv_local(pcols,pver), ptend_ql_local(pcols,pver)
  real(r8), intent(inout) :: ptend_qi_local(pcols,pver), ptend_nl_local(pcols,pver), ptend_ni_local(pcols,pver)
  integer :: i, k

  do k = top_lev, pver
     do i = 1, ncol_local
        ptend_s_local(i,k) = tlat_local(i,k)
        ptend_qv_local(i,k) = qvlat_local(i,k)
        ptend_ql_local(i,k) = qcten_local(i,k)
        ptend_qi_local(i,k) = qiten_local(i,k)
        ptend_nl_local(i,k) = ncten_local(i,k)
        ptend_ni_local(i,k) = niten_local(i,k)
     end do
  end do

end subroutine macrop_driver_ptend_assign_native

!============================================================================ !

subroutine macrop_driver_wtrc_split_tend(ncol_local, qcten_local, qiten_local, pqctn_local, nqctn_local, pqitn_local, nqitn_local)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local
  real(r8), target, intent(in) :: qcten_local(pcols,pver), qiten_local(pcols,pver)
  real(r8), target, intent(inout) :: pqctn_local(pcols,pver), nqctn_local(pcols,pver), pqitn_local(pcols,pver), nqitn_local(pcols,pver)

  interface
     subroutine macrop_driver_wtrc_split_tend_codon(ncol_c, pcols_c, pver_c, top_lev_c, qcten_p, qiten_p, pqctn_p, nqctn_p, &
          pqitn_p, nqitn_p) bind(c, name="macrop_driver_wtrc_split_tend_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c
       type(c_ptr), value :: qcten_p, qiten_p, pqctn_p, nqctn_p, pqitn_p, nqitn_p
     end subroutine macrop_driver_wtrc_split_tend_codon
  end interface

  if (use_native_impl) then
     call macrop_driver_wtrc_split_tend_native(ncol_local, qcten_local, qiten_local, pqctn_local, nqctn_local, pqitn_local, nqitn_local)
     return
  end if

  call macrop_driver_wtrc_split_tend_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
       int(top_lev, c_int64_t), c_loc(qcten_local), c_loc(qiten_local), c_loc(pqctn_local), c_loc(nqctn_local), &
       c_loc(pqitn_local), c_loc(nqitn_local))

end subroutine macrop_driver_wtrc_split_tend

!============================================================================ !

subroutine macrop_driver_wtrc_split_tend_native(ncol_local, qcten_local, qiten_local, pqctn_local, nqctn_local, pqitn_local, nqitn_local)

  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local
  real(r8), intent(in) :: qcten_local(pcols,pver), qiten_local(pcols,pver)
  real(r8), intent(inout) :: pqctn_local(pcols,pver), nqctn_local(pcols,pver), pqitn_local(pcols,pver), nqitn_local(pcols,pver)
  integer :: i, k

  do i = 1, ncol_local
     do k = top_lev, pver
        if(qcten_local(i,k) .lt. 0._r8) then
          nqctn_local(i,k) = qcten_local(i,k)
        else
          pqctn_local(i,k) = qcten_local(i,k)
        end if
        if(qiten_local(i,k) .lt. 0._r8) then
          nqitn_local(i,k) = qiten_local(i,k)
        else
          pqitn_local(i,k) = qiten_local(i,k)
        end if
     end do
  end do

end subroutine macrop_driver_wtrc_split_tend_native

!============================================================================ !

subroutine macrop_driver_cloud_mixing_diag(ncol_local, cld_local, state_ql_local, state_qi_local, mr_lsliq_local, mr_lsice_local)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local
  real(r8), target, intent(in) :: cld_local(pcols,pver), state_ql_local(pcols,pver), state_qi_local(pcols,pver)
  real(r8), target, intent(inout) :: mr_lsliq_local(pcols,pver), mr_lsice_local(pcols,pver)

  interface
     subroutine macrop_driver_cloud_mixing_diag_codon(ncol_c, pcols_c, pver_c, top_lev_c, cld_p, state_ql_p, state_qi_p, &
          mr_lsliq_p, mr_lsice_p) bind(c, name="macrop_driver_cloud_mixing_diag_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c
       type(c_ptr), value :: cld_p, state_ql_p, state_qi_p, mr_lsliq_p, mr_lsice_p
     end subroutine macrop_driver_cloud_mixing_diag_codon
  end interface

  if (use_native_impl) then
     call macrop_driver_cloud_mixing_diag_native(ncol_local, cld_local, state_ql_local, state_qi_local, mr_lsliq_local, mr_lsice_local)
     return
  end if

  call macrop_driver_cloud_mixing_diag_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
       int(top_lev, c_int64_t), c_loc(cld_local), c_loc(state_ql_local), c_loc(state_qi_local), c_loc(mr_lsliq_local), &
       c_loc(mr_lsice_local))

end subroutine macrop_driver_cloud_mixing_diag

!============================================================================ !

subroutine macrop_driver_cloud_mixing_diag_native(ncol_local, cld_local, state_ql_local, state_qi_local, mr_lsliq_local, mr_lsice_local)

  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local
  real(r8), intent(in) :: cld_local(pcols,pver), state_ql_local(pcols,pver), state_qi_local(pcols,pver)
  real(r8), intent(inout) :: mr_lsliq_local(pcols,pver), mr_lsice_local(pcols,pver)
  integer :: i, k

  do k = top_lev, pver
     do i = 1, ncol_local
        if (cld_local(i,k) .gt. 0._r8) then
           mr_lsliq_local(i,k) = state_ql_local(i,k)
           mr_lsice_local(i,k) = state_qi_local(i,k)
        else
           mr_lsliq_local(i,k) = 0._r8
           mr_lsice_local(i,k) = 0._r8
        end if
     end do
  end do

end subroutine macrop_driver_cloud_mixing_diag_native

!============================================================================ !

subroutine macrop_driver_store_state(ncol_local, state_t_local, state_qv_local, state_ql_local, state_qi_local, state_nl_local, &
     state_ni_local, tcwat_local, qcwat_local, lcwat_local, iccwat_local, nlwat_local, niwat_local, cldsice_local)

  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
  use physconst, only: tmelt
  use ref_pres,  only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local
  real(r8), target, intent(in) :: state_t_local(pcols,pver), state_qv_local(pcols,pver), state_ql_local(pcols,pver)
  real(r8), target, intent(in) :: state_qi_local(pcols,pver), state_nl_local(pcols,pver), state_ni_local(pcols,pver)
  real(r8), target, intent(inout) :: tcwat_local(pcols,pver), qcwat_local(pcols,pver), lcwat_local(pcols,pver)
  real(r8), target, intent(inout) :: iccwat_local(pcols,pver), nlwat_local(pcols,pver), niwat_local(pcols,pver)
  real(r8), target, intent(inout) :: cldsice_local(pcols,pver)

  interface
     subroutine macrop_driver_store_state_codon(ncol_c, pcols_c, pver_c, top_lev_c, tmelt_c, state_t_p, state_qv_p, &
          state_ql_p, state_qi_p, state_nl_p, state_ni_p, tcwat_p, qcwat_p, lcwat_p, iccwat_p, nlwat_p, niwat_p, &
          cldsice_p) bind(c, name="macrop_driver_store_state_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c
       real(c_double), value :: tmelt_c
       type(c_ptr), value :: state_t_p, state_qv_p, state_ql_p, state_qi_p, state_nl_p, state_ni_p, tcwat_p, qcwat_p
       type(c_ptr), value :: lcwat_p, iccwat_p, nlwat_p, niwat_p, cldsice_p
     end subroutine macrop_driver_store_state_codon
  end interface

  if (use_native_impl) then
     call macrop_driver_store_state_native(ncol_local, state_t_local, state_qv_local, state_ql_local, state_qi_local, &
          state_nl_local, state_ni_local, tcwat_local, qcwat_local, lcwat_local, iccwat_local, nlwat_local, niwat_local, &
          cldsice_local)
     return
  end if

  if (masterproc .and. .not. store_state_logged) then
     write(iulog,*) 'macrop_driver store state entered = codon'
     call macrop_driver_append_impl_proof('MACROP_DRIVER_MMPCOND_SHELL_PROOF_FILE', &
          'macrop_driver store state entered = codon')
     call flush(iulog)
     store_state_logged = .true.
  end if

  call macrop_driver_store_state_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
       int(top_lev, c_int64_t), tmelt, c_loc(state_t_local), c_loc(state_qv_local), c_loc(state_ql_local), c_loc(state_qi_local), &
       c_loc(state_nl_local), c_loc(state_ni_local), c_loc(tcwat_local), c_loc(qcwat_local), c_loc(lcwat_local), &
       c_loc(iccwat_local), c_loc(nlwat_local), c_loc(niwat_local), c_loc(cldsice_local))

end subroutine macrop_driver_store_state

!============================================================================ !

subroutine macrop_driver_store_state_native(ncol_local, state_t_local, state_qv_local, state_ql_local, state_qi_local, state_nl_local, &
     state_ni_local, tcwat_local, qcwat_local, lcwat_local, iccwat_local, nlwat_local, niwat_local, cldsice_local)

  use physconst, only: tmelt
  use ref_pres,  only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local
  real(r8), intent(in) :: state_t_local(pcols,pver), state_qv_local(pcols,pver), state_ql_local(pcols,pver)
  real(r8), intent(in) :: state_qi_local(pcols,pver), state_nl_local(pcols,pver), state_ni_local(pcols,pver)
  real(r8), intent(inout) :: tcwat_local(pcols,pver), qcwat_local(pcols,pver), lcwat_local(pcols,pver), iccwat_local(pcols,pver)
  real(r8), intent(inout) :: nlwat_local(pcols,pver), niwat_local(pcols,pver), cldsice_local(pcols,pver)
  integer :: i, k

  do k = top_lev, pver
     do i = 1, ncol_local
        tcwat_local(i,k) = state_t_local(i,k)
        qcwat_local(i,k) = state_qv_local(i,k)
        lcwat_local(i,k) = state_ql_local(i,k) + state_qi_local(i,k)
        iccwat_local(i,k) = state_qi_local(i,k)
        nlwat_local(i,k) = state_nl_local(i,k)
        niwat_local(i,k) = state_ni_local(i,k)
        cldsice_local(i,k) = lcwat_local(i,k) * min(1.0_r8, max(0.0_r8, (tmelt - tcwat_local(i,k)) / 20._r8))
     end do
  end do

end subroutine macrop_driver_store_state_native

!============================================================================ !

subroutine macrop_driver_append_impl_proof(env_name, proof_line)

  character(len=*), intent(in) :: env_name, proof_line
  character(len=512) :: proof_path
  integer :: status, n, unit_id

  proof_path = ''
  call get_environment_variable(env_name, value=proof_path, length=n, status=status)
  if (status /= 0 .or. n <= 0) return

  open(newunit=unit_id, file=trim(adjustl(proof_path(:n))), status='unknown', action='write', &
       position='append', iostat=status)
  if (status /= 0) return

  write(unit_id,'(A)') trim(proof_line)
  close(unit_id)

end subroutine macrop_driver_append_impl_proof

!============================================================================ !

subroutine macrop_driver_select_wtrc_shell_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (wtrc_shell_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('MACROP_DRIVER_WTRC_SHELL_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_wtrc_shell_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_wtrc_shell_impl = .false.
  end if

  wtrc_shell_impl_selected = .true.

  if (masterproc) then
     if (use_native_wtrc_shell_impl) then
        write(iulog,*) 'macrop_driver_wtrc_shell implementation = native'
        call macrop_driver_append_impl_proof('MACROP_DRIVER_WTRC_SHELL_PROOF_FILE', &
             'macrop_driver_wtrc_shell implementation = native')
     else
        write(iulog,*) 'macrop_driver_wtrc_shell implementation = codon'
        call macrop_driver_append_impl_proof('MACROP_DRIVER_WTRC_SHELL_PROOF_FILE', &
             'macrop_driver_wtrc_shell implementation = codon')
     end if
     call flush(iulog)
  end if

end subroutine macrop_driver_select_wtrc_shell_impl

!============================================================================ !

subroutine macrop_driver_select_wtrc_process_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (wtrc_process_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('MACROP_DRIVER_WTRC_PROCESS_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_wtrc_process_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_wtrc_process_impl = .false.
  end if

  wtrc_process_impl_selected = .true.

  if (masterproc) then
     if (use_native_wtrc_process_impl) then
        write(iulog,*) 'macrop_driver_wtrc_process implementation = native'
        call macrop_driver_append_impl_proof('MACROP_DRIVER_WTRC_PROCESS_PROOF_FILE', &
             'macrop_driver_wtrc_process implementation = native')
     else
        write(iulog,*) 'macrop_driver_wtrc_process implementation = codon'
        call macrop_driver_append_impl_proof('MACROP_DRIVER_WTRC_PROCESS_PROOF_FILE', &
             'macrop_driver_wtrc_process implementation = codon')
     end if
     call flush(iulog)
  end if

end subroutine macrop_driver_select_wtrc_process_impl

!============================================================================ !

subroutine macrop_driver_wtrc_shell_codon_wrap(ncol_local, dtime_local, state_q_local, state_t_local, state_pmid_local, &
     ptend_q_local, qvlat_local, qcten_local, qiten_local, prelat_local, process_rates_local)

  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
  use constituents, only: pcnst
  use physconst, only: cpair, epsilo
  use ref_pres, only: top_lev => trop_cloud_top_lev
  use water_isotopes, only: pwtspec
  use water_tracer_vars, only: wisotope, wtrc_iatype, wtrc_iawset, wtrc_bulk_indices, wtrc_indices, &
       wtrc_ncnst, wtrc_niter, wtrc_nwset, iwspec, wtrc_qmin
  use water_tracers, only: wtrc_get_rstd
  use water_types, only: pwtype, iwtvap, iwtliq, iwtice

  integer, intent(in) :: ncol_local
  real(r8), intent(in) :: dtime_local
  real(r8), target, intent(in) :: state_q_local(pcols,pver,pcnst), state_t_local(pcols,pver), state_pmid_local(pcols,pver)
  real(r8), target, intent(inout) :: ptend_q_local(pcols,pver,pcnst)
  real(r8), target, intent(in) :: qvlat_local(pcols,pver), qcten_local(pcols,pver), qiten_local(pcols,pver)
  real(r8), target, intent(in) :: prelat_local(pcols,pver)
  real(r8), target, intent(inout) :: process_rates_local(pcols,pver,pwtype,pwtype,pwtype)

  integer :: ispec
  integer(c_int64_t), target :: wisotope_c
  integer(c_int64_t), target :: wtrc_iawset64(pwtype,wtrc_nwset)
  integer(c_int64_t), target :: wtrc_iatype64(wtrc_nwset,pwtype)
  integer(c_int64_t), target :: wtrc_bulk_indices64(pwtype)
  integer(c_int64_t), target :: wtrc_indices64(wtrc_ncnst)
  integer(c_int64_t), target :: iwspec64(pcnst)
  real(c_double), target :: rstd(pwtspec)
  real(r8), target :: qloc_local(pcols,pver,pcnst)
  real(r8), target :: qloc0_local(pcols,pver,pcnst)
  real(r8), target :: tloc_local(pcols,pver)
  real(r8), target :: diff_local(pcols,pver,pwtype)

  interface
     subroutine macrop_driver_wtrc_shell_codon(ncol_c, pcols_c, pver_c, pcnst_c, pwtype_c, top_lev_c, wtrc_niter_c, &
          wtrc_ncnst_c, wtrc_nwset_c, wisotope_c, iwtvap_c, iwtliq_c, iwtice_c, cpair_c, dtime_c, wtrc_qmin_c, epsilo_c, &
          state_q_p, state_t_p, state_pmid_p, ptend_q_p, qvlat_p, qcten_p, qiten_p, prelat_p, process_rates_p, qloc_p, &
          qloc0_p, tloc_p, diff_p, wtrc_iawset_p, wtrc_iatype_p, wtrc_bulk_indices_p, wtrc_indices_p, iwspec_p, rstd_p) &
          bind(c, name="macrop_driver_wtrc_shell_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pcnst_c, pwtype_c, top_lev_c, wtrc_niter_c
       integer(c_int64_t), value :: wtrc_ncnst_c, wtrc_nwset_c, wisotope_c, iwtvap_c, iwtliq_c, iwtice_c
       real(c_double), value :: cpair_c, dtime_c, wtrc_qmin_c, epsilo_c
       type(c_ptr), value :: state_q_p, state_t_p, state_pmid_p, ptend_q_p, qvlat_p, qcten_p, qiten_p, prelat_p
       type(c_ptr), value :: process_rates_p, qloc_p, qloc0_p, tloc_p, diff_p
       type(c_ptr), value :: wtrc_iawset_p, wtrc_iatype_p, wtrc_bulk_indices_p, wtrc_indices_p, iwspec_p, rstd_p
     end subroutine macrop_driver_wtrc_shell_codon
  end interface

  do ispec = 1, pwtspec
     rstd(ispec) = real(wtrc_get_rstd(ispec), c_double)
  end do
  do ispec = 1, wtrc_nwset
     wtrc_iawset64(:,ispec) = int(wtrc_iawset(:,ispec), c_int64_t)
  end do
  do ispec = 1, pwtype
     wtrc_iatype64(:,ispec) = int(wtrc_iatype(1:wtrc_nwset,ispec), c_int64_t)
     wtrc_bulk_indices64(ispec) = int(wtrc_bulk_indices(ispec), c_int64_t)
  end do
  do ispec = 1, wtrc_ncnst
     wtrc_indices64(ispec) = int(wtrc_indices(ispec), c_int64_t)
  end do
  do ispec = 1, pcnst
     iwspec64(ispec) = int(iwspec(ispec), c_int64_t)
  end do

  wisotope_c = merge(1_c_int64_t, 0_c_int64_t, wisotope)

  if (ncol_local > 0 .and. top_lev <= pver) then
     call macrop_driver_wtrc_shell_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
          int(pcnst, c_int64_t), int(pwtype, c_int64_t), int(top_lev, c_int64_t), int(wtrc_niter, c_int64_t), &
          int(wtrc_ncnst, c_int64_t), int(wtrc_nwset, c_int64_t), wisotope_c, int(iwtvap, c_int64_t), &
          int(iwtliq, c_int64_t), int(iwtice, c_int64_t), real(cpair, c_double), real(dtime_local, c_double), &
          real(wtrc_qmin, c_double), real(epsilo, c_double), c_loc(state_q_local), c_loc(state_t_local), c_loc(state_pmid_local), &
          c_loc(ptend_q_local), c_loc(qvlat_local), c_loc(qcten_local), c_loc(qiten_local), c_loc(prelat_local), &
          c_loc(process_rates_local), c_loc(qloc_local), c_loc(qloc0_local), c_loc(tloc_local), c_loc(diff_local), &
          c_loc(wtrc_iawset64), c_loc(wtrc_iatype64), c_loc(wtrc_bulk_indices64), c_loc(wtrc_indices64), c_loc(iwspec64), &
          c_loc(rstd))
  end if

end subroutine macrop_driver_wtrc_shell_codon_wrap

!============================================================================ !

subroutine macrop_driver_wtrc_process_rates_codon_wrap(ncol_local, qvlat_local, qcten_local, qiten_local, process_rates_local)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use water_types, only: pwtype, iwtvap, iwtliq, iwtice
  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local
  real(r8), target, intent(in) :: qvlat_local(pcols,pver), qcten_local(pcols,pver), qiten_local(pcols,pver)
  real(r8), target, intent(inout) :: process_rates_local(pcols,pver,pwtype,pwtype,pwtype)

  interface
     subroutine macrop_driver_wtrc_process_rates_codon(ncol_c, pcols_c, pver_c, pwtype_c, top_lev_c, iwtvap_c, iwtliq_c, &
          iwtice_c, qvlat_p, qcten_p, qiten_p, process_rates_p) bind(c, name="macrop_driver_wtrc_process_rates_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pwtype_c, top_lev_c
       integer(c_int64_t), value :: iwtvap_c, iwtliq_c, iwtice_c
       type(c_ptr), value :: qvlat_p, qcten_p, qiten_p, process_rates_p
     end subroutine macrop_driver_wtrc_process_rates_codon
  end interface

  if (ncol_local > 0 .and. top_lev <= pver) then
     call macrop_driver_wtrc_process_rates_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
          int(pwtype, c_int64_t), int(top_lev, c_int64_t), int(iwtvap, c_int64_t), int(iwtliq, c_int64_t), &
          int(iwtice, c_int64_t), c_loc(qvlat_local), c_loc(qcten_local), c_loc(qiten_local), c_loc(process_rates_local))
  end if

end subroutine macrop_driver_wtrc_process_rates_codon_wrap

!============================================================================ !

subroutine macrop_driver_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('MACROP_DRIVER_IMPL', value=impl_name, length=n, status=status)

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
        write(iulog,*) 'macrop_driver_tend implementation = native'
     else
        write(iulog,*) 'macrop_driver_tend implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine macrop_driver_select_impl

!============================================================================ !
!                                                                             !
!============================================================================ !

subroutine macrop_driver_select_branches(micro_do_icesupersat_in, trace_water_in, &
     wtrc_detrain_in_macrop_in, cu_det_st_in, use_shfrc_in, do_cldice_in, do_cldliq_in, do_detrain_in)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr

  logical, intent(in) :: micro_do_icesupersat_in
  logical, intent(in) :: trace_water_in
  logical, intent(in) :: wtrc_detrain_in_macrop_in
  logical, intent(in) :: cu_det_st_in
  logical, intent(in) :: use_shfrc_in
  logical, intent(in) :: do_cldice_in
  logical, intent(in) :: do_cldliq_in
  logical, intent(in) :: do_detrain_in

  integer(c_int64_t), target :: branch_mask_c

  interface
     subroutine macrop_driver_select_branches_codon(micro_do_icesupersat_c, trace_water_c, &
          wtrc_detrain_in_macrop_c, cu_det_st_c, use_shfrc_c, do_cldice_c, do_cldliq_c, do_detrain_c, &
          branch_mask_p) bind(c, name="macrop_driver_select_branches_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: micro_do_icesupersat_c, trace_water_c
       integer(c_int64_t), value :: wtrc_detrain_in_macrop_c, cu_det_st_c, use_shfrc_c
       integer(c_int64_t), value :: do_cldice_c, do_cldliq_c, do_detrain_c
       type(c_ptr), value :: branch_mask_p
     end subroutine macrop_driver_select_branches_codon
  end interface

  if (branch_selected) return

  branch_mask_c = 0_c_int64_t
  call macrop_driver_select_branches_codon( &
       merge(1_c_int64_t, 0_c_int64_t, micro_do_icesupersat_in), &
       merge(1_c_int64_t, 0_c_int64_t, trace_water_in), &
       merge(1_c_int64_t, 0_c_int64_t, wtrc_detrain_in_macrop_in), &
       merge(1_c_int64_t, 0_c_int64_t, cu_det_st_in), &
       merge(1_c_int64_t, 0_c_int64_t, use_shfrc_in), &
       merge(1_c_int64_t, 0_c_int64_t, do_cldice_in), &
       merge(1_c_int64_t, 0_c_int64_t, do_cldliq_in), &
       merge(1_c_int64_t, 0_c_int64_t, do_detrain_in), &
       c_loc(branch_mask_c) &
  )

  branch_mask = int(branch_mask_c)
  branch_selected = .true.

end subroutine macrop_driver_select_branches

! Saturation adjustment for ice
! Add ice mass if supersaturated
elemental subroutine ice_macro_tend(naai,t,p,qv,qi,ni,xxls,deltat,stend,qvtend,qitend,nitend) 

  use wv_sat_methods, only: wv_sat_qsat_ice

  real(r8), intent(in)  :: naai   !Activated number of ice nuclei 
  real(r8), intent(in)  :: t      !temperature (k)
  real(r8), intent(in)  :: p      !pressure (pa0
  real(r8), intent(in)  :: qv     !water vapor mixing ratio
  real(r8), intent(in)  :: qi     !ice mixing ratio
  real(r8), intent(in)  :: ni     !ice number concentration
  real(r8), intent(in)  :: xxls   !latent heat of sublimation
  real(r8), intent(in)  :: deltat !timestep
  real(r8), intent(out) :: stend  ! 'temperature' tendency 
  real(r8), intent(out) :: qvtend !vapor tendency
  real(r8), intent(out) :: qitend !ice mass tendency
  real(r8), intent(out) :: nitend !ice number tendency  
 
  real(r8) :: ESI
  real(r8) :: QSI
  real(r8) :: tau
  logical  :: tau_constant

  tau_constant = .true.

  stend = 0._r8
  qvtend = 0._r8
  qitend = 0._r8
  nitend = 0._r8

  ! calculate qsati from t,p,q

  call wv_sat_qsat_ice(t, p, ESI, QSI)

  if (naai.gt.1.e-18_r8.and.qv.gt.QSI) then

     !optional timescale on condensation
     !tau in sections. Try 300. or tau = f(T): 300s  t> 268, 1800s for t<238
     !     
     if (.not. tau_constant) then
        if( t.gt. 268.15_r8 ) then
           tau = 300.0_r8
        elseif(t.lt.238.15_r8 ) then
           tau = 1800._r8
        else
           tau = 300._r8 + (1800._r8 - 300._r8) * ( 268.15_r8 - t ) / 30._r8
        endif
     else
         tau = 300._r8
     end if

     qitend = (qv-QSI)/deltat !* exp(-tau/deltat)
     qvtend = 0._r8 - qitend
     stend  = qitend * xxls    ! moist static energy tend...[J/kg/s] !

     ! kg(h2o)/kg(air)/s * J/kg(h2o)  = J/kg(air)/s (=W/kg)
     ! if ice exists (more than 1 L-1) and there is condensation, do not add to number (= growth), else, add 10um ice

     if (ni.lt.1.e3_r8.and.(qi+qitend*deltat).gt.1e-18_r8) then
        nitend = nitend + 3._r8 * qitend/(4._r8*3.14_r8* 10.e-6_r8**3*997._r8)
     endif

  endif

end subroutine ice_macro_tend

end module macrop_driver
