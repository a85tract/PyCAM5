module zm_conv_intr
!---------------------------------------------------------------------------------
! Purpose:
!
! CAM interface to the Zhang-McFarlane deep convection scheme
!
! Author: D.B. Coleman
! January 2010 modified by J. Kay to add COSP simulator fields to physics buffer
!
!Water isotopes/tracers added by J. Nusbaumer
!
!---------------------------------------------------------------------------------
   use shr_kind_mod, only: r8=>shr_kind_r8
   use physconst,    only: cpair                              
   use ppgrid,       only: pver, pcols, pverp, begchunk, endchunk
   use zm_conv,      only: zm_conv_evap, zm_convr, convtran, momtran
   use cam_history,  only: outfld, addfld, phys_decomp
   use perf_mod
   use cam_logfile,  only: iulog
   use constituents, only: cnst_add
   use spmd_utils,   only: masterproc
   
   implicit none
   private
   save

   ! Public methods

   public ::&
      zm_conv_register,           &! register fields in physics buffer
      zm_conv_readnl,             &! read namelist
      zm_conv_init,               &! initialize donner_deep module
      zm_conv_tend,               &! return tendencies
      zm_conv_tend_2               ! return tendencies

   ! Private module data

   real(r8), allocatable, target, dimension(:,:,:) :: mu  !(pcols,pver,begchunk:endchunk)
   real(r8), allocatable, target, dimension(:,:,:) :: eu  !(pcols,pver,begchunk:endchunk)
   real(r8), allocatable, target, dimension(:,:,:) :: du  !(pcols,pver,begchunk:endchunk)
   real(r8), allocatable, target, dimension(:,:,:) :: md  !(pcols,pver,begchunk:endchunk)
   real(r8), allocatable, target, dimension(:,:,:) :: ed  !(pcols,pver,begchunk:endchunk)
   real(r8), allocatable, target, dimension(:,:,:) :: dp  !(pcols,pver,begchunk:endchunk)
        ! wg layer thickness in mbs (between upper/lower interface).
   real(r8), allocatable, dimension(:,:)   :: dsubcld  !(pcols,begchunk:endchunk)
        ! wg layer thickness in mbs between lcl and maxi.

   integer, allocatable, dimension(:,:) :: jt   !(pcols,begchunk:endchunk)
        ! wg top  level index of deep cumulus convection.
   integer, allocatable, dimension(:,:) :: maxg !(pcols,begchunk:endchunk)
        ! wg gathered values of maxi.
   integer, allocatable, dimension(:,:) :: ideep !(pcols,begchunk:endchunk)               
        ! w holds position of gathered points vs longitude index

   integer, allocatable, dimension(:) :: lengath !(begchunk:endchunk)

   integer ::& ! indices for fields in the physics buffer
      dp_flxprc_idx, &
      dp_flxsnw_idx, &
      dp_cldliq_idx, &
      ixorg,       &
      dp_cldice_idx, &
      prec_dp_idx,   &
      snow_dp_idx

   real(r8), parameter :: unset_r8 = huge(1.0_r8)
   real(r8) :: zmconv_c0_lnd = unset_r8
   real(r8) :: zmconv_c0_ocn = unset_r8
   real(r8) :: zmconv_ke     = unset_r8
   real(r8) :: zmconv_ke_lnd = unset_r8
   logical  :: zmconv_org                !  Parameterization for sub-grid scale convective organization for the ZM deep 
                                         !  convective scheme based on Mapes and Neale (2011)

   !Needed for Water tracers:
   integer  limcnv           !top interface level limit for convection

!  indices for fields in the physics buffer
   integer  ::    cld_idx          = 0    
   integer  ::    icwmrdp_idx      = 0     
   integer  ::    rprddp_idx       = 0    
   integer  ::    fracis_idx       = 0   
   integer  ::    nevapr_dpcu_idx  = 0    

   logical :: use_native_zm_post_shell = .false.
   logical :: zm_post_shell_selected = .false.
   logical :: zm_convr_post_logged = .false.
   logical :: zm_conv_evap_prep_logged = .false.
   logical :: zm_conv_evap_post_logged = .false.
   logical :: zm_conv_evap_hist_logged = .false.
   logical :: zm_momtran_post_logged = .false.
   logical :: zm_convtran1_ratio_logged = .false.
   logical :: zm_convtran2_dpdry_logged = .false.


!=========================================================================================
contains
!=========================================================================================

subroutine zm_conv_register

!----------------------------------------
! Purpose: register fields with the physics buffer
!----------------------------------------

  use physics_buffer, only : pbuf_add_field, dtype_r8

  implicit none

  integer idx

! Flux of precipitation from deep convection (kg/m2/s)
   call pbuf_add_field('DP_FLXPRC','global',dtype_r8,(/pcols,pverp/),dp_flxprc_idx) 

! Flux of snow from deep convection (kg/m2/s) 
   call pbuf_add_field('DP_FLXSNW','global',dtype_r8,(/pcols,pverp/),dp_flxsnw_idx) 

! deep gbm cloud liquid water (kg/kg)
   call pbuf_add_field('DP_CLDLIQ','global',dtype_r8,(/pcols,pver/), dp_cldliq_idx)  

! deep gbm cloud liquid water (kg/kg)    
   call pbuf_add_field('DP_CLDICE','global',dtype_r8,(/pcols,pver/), dp_cldice_idx)  

   if (zmconv_org) then
      call cnst_add('ZM_ORG',0._r8,0._r8,0._r8,ixorg,longname='organization parameter')
   endif

end subroutine zm_conv_register

!=========================================================================================

subroutine zm_conv_readnl(nlfile)

   use cam_abortutils,  only: endrun
   use spmd_utils,      only: masterproc
   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'zm_conv_readnl'

   namelist /zmconv_nl/ zmconv_c0_lnd, zmconv_c0_ocn, zmconv_ke, zmconv_ke_lnd, zmconv_org
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'zmconv_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, zmconv_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)

   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(zmconv_c0_lnd,            1, mpir8,  0, mpicom)
   call mpibcast(zmconv_c0_ocn,            1, mpir8,  0, mpicom)
   call mpibcast(zmconv_ke,                1, mpir8,  0, mpicom)
   call mpibcast(zmconv_ke_lnd,            1, mpir8,  0, mpicom)
   call mpibcast(zmconv_org,               1, mpilog, 0, mpicom)
#endif

end subroutine zm_conv_readnl

!=========================================================================================

subroutine zm_conv_init(pref_edge)

!----------------------------------------
! Purpose:  declare output fields, initialize variables needed by convection
!----------------------------------------

  use cam_history,    only: outfld, addfld, add_default, phys_decomp
  use ppgrid,         only: pcols, pver
  use zm_conv,        only: zm_convi
  use pmgrid,         only: plev,plevp
  use spmd_utils,     only: masterproc
  use error_messages, only: alloc_err
  use phys_control,   only: phys_deepconv_pbl, phys_getopts, cam_physpkg_is
  use physics_buffer, only: pbuf_get_index

  implicit none

  real(r8),intent(in) :: pref_edge(plevp)        ! reference pressures at interfaces


  logical :: no_deep_pbl    ! if true, no deep convection in PBL
!  integer  limcnv           ! top interface level limit for convection
  integer k, istat
  logical :: history_budget ! output tendencies and state variables for CAM4
                            ! temperature, water vapor, cloud ice and cloud
                            ! liquid budgets.
  integer :: history_budget_histfile_num ! output history file number for budget fields

!
! Allocate space for arrays private to this module
!
     allocate( mu(pcols,pver,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'zm_conv_tend', 'mu', &
                      pcols*pver*((endchunk-begchunk)+1) )
     allocate( eu(pcols,pver,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'zm_conv_tend', 'eu', &
                      pcols*pver*((endchunk-begchunk)+1) )
     allocate( du(pcols,pver,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'zm_conv_tend', 'du', &
                      pcols*pver*((endchunk-begchunk)+1) )
     allocate( md(pcols,pver,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'zm_conv_tend', 'md', &
                      pcols*pver*((endchunk-begchunk)+1) )
     allocate( ed(pcols,pver,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'zm_conv_tend', 'ed', &
                      pcols*pver*((endchunk-begchunk)+1) )
     allocate( dp(pcols,pver,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'zm_conv_tend', 'dp', &
                      pcols*pver*((endchunk-begchunk)+1) )
     allocate( dsubcld(pcols,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'zm_conv_tend', 'dsubcld', &
                      pcols*((endchunk-begchunk)+1) )
     allocate( jt(pcols,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'zm_conv_tend', 'jt', &
                      pcols*((endchunk-begchunk)+1) )
     allocate( maxg(pcols,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'zm_conv_tend', 'maxg', &
                      pcols*((endchunk-begchunk)+1) )
     allocate( ideep(pcols,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'zm_conv_tend', 'ideep', &
                      pcols*((endchunk-begchunk)+1) )
     allocate( lengath(begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'zm_conv_tend', 'lengath', &
                      ((endchunk-begchunk)+1) )


! 
! Register fields with the output buffer
!

    if (zmconv_org) then
       call addfld ('ZM_ORG     ','-       ',pver, 'A','Organization parameter',                        phys_decomp)
       call addfld ('ZM_ORG2D   ','-       ',pver, 'A','Organization parameter 2D',                     phys_decomp)
    endif
    call addfld ('PRECZ   ','m/s     ',1,    'A','total precipitation from ZM convection',        phys_decomp)
    call addfld ('ZMDT    ','K/s     ',pver, 'A','T tendency - Zhang-McFarlane moist convection', phys_decomp)
    call addfld ('ZMDQ    ','kg/kg/s ',pver, 'A','Q tendency - Zhang-McFarlane moist convection', phys_decomp)
    call addfld ('ZMDICE ','kg/kg/s ',pver, 'A','Cloud ice tendency - Zhang-McFarlane convection',phys_decomp)
    call addfld ('ZMDLIQ ','kg/kg/s ',pver, 'A','Cloud liq tendency - Zhang-McFarlane convection',phys_decomp)
    call addfld ('EVAPTZM ','K/s     ',pver, 'A','T tendency - Evaporation/snow prod from Zhang convection',phys_decomp)
    call addfld ('FZSNTZM ','K/s     ',pver, 'A','T tendency - Rain to snow conversion from Zhang convection',phys_decomp)
    call addfld ('EVSNTZM ','K/s     ',pver, 'A','T tendency - Snow to rain prod from Zhang convection',phys_decomp)
    call addfld ('EVAPQZM ','kg/kg/s ',pver, 'A','Q tendency - Evaporation from Zhang-McFarlane moist convection',phys_decomp)
    
    call addfld ('ZMFLXPRC','kg/m2/s ',pverp, 'A','Flux of precipitation from ZM convection'       ,phys_decomp)
    call addfld ('ZMFLXSNW','kg/m2/s ',pverp, 'A','Flux of snow from ZM convection'                ,phys_decomp)
    call addfld ('ZMNTPRPD','kg/kg/s ',pver , 'A','Net precipitation production from ZM convection',phys_decomp)
    call addfld ('ZMNTSNPD','kg/kg/s ',pver , 'A','Net snow production from ZM convection'         ,phys_decomp)
    call addfld ('ZMEIHEAT','W/kg'    ,pver , 'A','Heating by ice and evaporation in ZM convection',phys_decomp)
    
    call addfld ('CMFMCDZM','kg/m2/s ',pverp,'A','Convection mass flux from ZM deep ',phys_decomp)
    call addfld ('PRECCDZM','m/s     ',1,    'A','Convective precipitation rate from ZM deep',phys_decomp)

    call addfld ('PCONVB','Pa'    ,1 , 'A','convection base pressure',phys_decomp)
    call addfld ('PCONVT','Pa'    ,1 , 'A','convection top  pressure',phys_decomp)

    call addfld ('CAPE',   'J/kg',       1, 'A', 'Convectively available potential energy', phys_decomp)
    call addfld ('FREQZM ','fraction  ',1  ,'A', 'Fractional occurance of ZM convection',phys_decomp) 

    call addfld ('ZMMTT ', 'K/s',     pver, 'A', 'T tendency - ZM convective momentum transport',phys_decomp)
    call addfld ('ZMMTU',  'm/s2',    pver, 'A', 'U tendency - ZM convective momentum transport',  phys_decomp)
    call addfld ('ZMMTV',  'm/s2',    pver, 'A', 'V tendency - ZM convective momentum transport',  phys_decomp)

    call addfld ('ZMMU',   'kg/m2/s', pver, 'A', 'ZM convection updraft mass flux',   phys_decomp)
    call addfld ('ZMMD',   'kg/m2/s', pver, 'A', 'ZM convection downdraft mass flux', phys_decomp)

    call addfld ('ZMUPGU', 'm/s2',    pver, 'A', 'zonal force from ZM updraft pressure gradient term',       phys_decomp)
    call addfld ('ZMUPGD', 'm/s2',    pver, 'A', 'zonal force from ZM downdraft pressure gradient term',     phys_decomp)
    call addfld ('ZMVPGU', 'm/s2',    pver, 'A', 'meridional force from ZM updraft pressure gradient term',  phys_decomp)
    call addfld ('ZMVPGD', 'm/s2',    pver, 'A', 'merdional force from ZM downdraft pressure gradient term', phys_decomp)

    call addfld ('ZMICUU', 'm/s',     pver, 'A', 'ZM in-cloud U updrafts',      phys_decomp)
    call addfld ('ZMICUD', 'm/s',     pver, 'A', 'ZM in-cloud U downdrafts',    phys_decomp)
    call addfld ('ZMICVU', 'm/s',     pver, 'A', 'ZM in-cloud V updrafts',      phys_decomp)
    call addfld ('ZMICVD', 'm/s',     pver, 'A', 'ZM in-cloud V downdrafts',    phys_decomp)
    
    call phys_getopts( history_budget_out = history_budget, &
                       history_budget_histfile_num_out = history_budget_histfile_num)

    if (zmconv_org) then
       call add_default('ZM_ORG', 1, ' ')
       call add_default('ZM_ORG2D', 1, ' ')
    endif
    if ( history_budget ) then
       call add_default('EVAPTZM  ', history_budget_histfile_num, ' ')
       call add_default('EVAPQZM  ', history_budget_histfile_num, ' ')
       call add_default('ZMDT     ', history_budget_histfile_num, ' ')
       call add_default('ZMDQ     ', history_budget_histfile_num, ' ')
       call add_default('ZMDLIQ   ', history_budget_histfile_num, ' ')
       call add_default('ZMDICE   ', history_budget_histfile_num, ' ')

       if( cam_physpkg_is('cam4') .or. cam_physpkg_is('cam5') ) then
          call add_default('ZMMTT    ', history_budget_histfile_num, ' ')
       end if

    end if
!
! Limit deep convection to regions below 40 mb
! Note this calculation is repeated in the shallow convection interface
!
    limcnv = 0   ! null value to check against below
    if (pref_edge(1) >= 4.e3_r8) then
       limcnv = 1
    else
       do k=1,plev
          if (pref_edge(k) < 4.e3_r8 .and. pref_edge(k+1) >= 4.e3_r8) then
             limcnv = k
             exit
          end if
       end do
       if ( limcnv == 0 ) limcnv = plevp
    end if
    
    if (masterproc) then
       write(iulog,*)'ZM_CONV_INIT: Deep convection will be capped at intfc ',limcnv, &
            ' which is ',pref_edge(limcnv),' pascals'
    end if
        
    no_deep_pbl = phys_deepconv_pbl()
    call zm_convi(limcnv,zmconv_c0_lnd, zmconv_c0_ocn, zmconv_ke, zmconv_ke_lnd, &
                  zmconv_org, no_deep_pbl_in = no_deep_pbl)

    cld_idx         = pbuf_get_index('CLD')
    icwmrdp_idx     = pbuf_get_index('ICWMRDP')
    rprddp_idx      = pbuf_get_index('RPRDDP')
    fracis_idx      = pbuf_get_index('FRACIS')
    nevapr_dpcu_idx = pbuf_get_index('NEVAPR_DPCU')
    prec_dp_idx     = pbuf_get_index('PREC_DP')
    snow_dp_idx     = pbuf_get_index('SNOW_DP')

end subroutine zm_conv_init
!=========================================================================================
!subroutine zm_conv_tend(state, ptend, tdt)

subroutine zm_conv_tend(pblh    ,mcon    ,cme     , &
     tpert   ,dlf     ,pflx    ,zdu      , &
     rliq    , &
     ztodt   , &
     jctop   ,jcbot , &
     state   ,ptend_all   ,landfrac   ,pbuf , &
     wtdlf   ) 
  

   use cam_history,   only: outfld
   use physics_types, only: physics_state, physics_ptend
   use physics_types, only: physics_ptend_init, physics_update
   use physics_types, only: physics_state_copy, physics_state_dealloc
   use physics_types, only: physics_ptend_sum, physics_ptend_dealloc

   use phys_grid,     only: get_lat_p, get_lon_p
   use time_manager,  only: get_nstep, is_first_step
   use physics_buffer, only: pbuf_get_field, physics_buffer_desc, pbuf_old_tim_idx,&
                             pbuf_get_index
   use constituents,  only: pcnst, cnst_get_ind, cnst_is_convtran1
   use check_energy,  only: check_energy_chng
   use physconst,     only: gravit
   use phys_control,  only: cam_physpkg_is
   use water_types,   only: iwtvap, iwtundef, iwtliq, iwtice, iwtcvrain, iwtcvsnow
   use water_tracer_vars, only: trace_water, wtrc_iatype, wtrc_nwset, wtrc_iawset, &
                                wtrc_srfpcp_indices, wtrc_bulk_indices
   use water_tracers, only: wtrc_check_h2o, wtrc_precip_evap, wtrc_q1q2_pjr
   use constituents,  only: cnst_name
   use iso_c_binding, only: c_int64_t

   ! Arguments

   type(physics_state), intent(in),target   :: state          ! Physics state variables
   type(physics_ptend), intent(out)         :: ptend_all      ! individual parameterization tendencies
   type(physics_buffer_desc), pointer       :: pbuf(:)

   real(r8), intent(in) :: ztodt                       ! 2 delta t (model time increment)
   real(r8), intent(in) :: pblh(pcols)                 ! Planetary boundary layer height
   real(r8), intent(in) :: tpert(pcols)                ! Thermal temperature excess
   real(r8), intent(in) :: landfrac(pcols)             ! RBN - Landfrac 

   real(r8), target, intent(out) :: mcon(pcols,pverp)  ! Convective mass flux--m sub c
   real(r8), intent(out) :: dlf(pcols,pver)    ! scattrd version of the detraining cld h2o tend
   real(r8), intent(out) :: pflx(pcols,pverp)  ! scattered precip flux at each level
   real(r8), intent(out) :: cme(pcols,pver)    ! cmf condensation - evaporation
   real(r8), intent(out) :: zdu(pcols,pver)    ! detraining mass flux

   real(r8), intent(out) :: rliq(pcols) ! reserved liquid (not yet in cldliq) for energy integrals

   !Water tracers:
   real(r8), intent(out) :: wtdlf(pcols,pver,wtrc_nwset) ! Detraining water tracer cld from convection

   ! Local variables

   integer :: i,k,m
   integer :: ilon                      ! global longitude index of a column
   integer :: ilat                      ! global latitude index of a column
   integer :: nstep
   integer :: ixcldice, ixcldliq      ! constituent indices for cloud liquid and ice water.
   integer :: lchnk                   ! chunk identifier
   integer :: ncol                    ! number of atmospheric columns
   integer :: itim_old                ! for physics buffer fields

   real(r8), target :: ftem(pcols,pver)              ! Temporary workspace for outfld variables
   real(r8) :: ntprprd(pcols,pver)    ! evap outfld: net precip production in layer
   real(r8) :: ntsnprd(pcols,pver)    ! evap outfld: net snow production in layer
   real(r8), target :: tend_s_snwprd  (pcols,pver) ! Heating rate of snow production
   real(r8), target :: tend_s_snwevmlt(pcols,pver) ! Heating rate of evap/melting of snow
   real(r8) :: fake_dpdry(pcols,pver) ! used in convtran call

   !----------------------
   !Local water tracer variables (should probably be cleaned up)
   !---------------------- 
!   real(r8) wtdlft(pcols,pver,pcnst)      ! Temporary Detraining water tracer cld from convection
   real(r8) wtprect(pcols,pcnst)          ! Temporary tracer total precipitation from ZM convection
   real(r8) wtsnowt(pcols,pcnst)          ! Temproary tracer snow from ZM convection
!   real(r8) wtcmet(pcols,pver,pcnst)      ! Temporary condensation - evaporation 
!   real(r8) wtql(pcols,pver,pcnst)        ! updraft cloud liquid water ZM scheme
!   real(r8) wtqi(pcols,pver,pcnst)        ! updraft cloud ice water ZM scheme
!   real(r8) wtqu(pcols,pver,pcnst)        ! updraft water vapour ZM scheme
   real(r8) qu(pcols,pver)                ! updraft water vapour ZM scheme 
   real(r8) qd(pcols,pver)                ! downdraft water vapour ZM scheme 
   real(r8) dz(pcols,pver)                ! layer thickness
   real(r8) eps0(pcols)                   ! Not sure...
   real(r8) cu(pcols,pver)                ! condensation rate
   real(r8) evp(pcols,pver)               ! evaporation rate
   real(r8) tu(pcols,pver)                ! temperature of updraft
   real(r8) td(pcols,pver)                ! temperature of downdraft   
!   real(r8) wtqd(pcols,pver,pcnst)        ! downdraft water vapour ZM scheme
   real(r8) wtrprd(pcols,pver,pcnst)      ! rain production in ZM convection
   real(r8) rppe(pcols,pver)              ! ZM rain production pre-evaporation
   real(r8) evpstore(pcols,pver)          ! bulk precipiation evaporation rate
   real(r8) substore(pcols,pver)          ! bulk snow sublimation rate

   real(r8) qst(pcols,pver)               ! bulk environmental saturation mixing ratio
   real(r8) qsthat(pcols,pver)            ! qst at interfaces
   real(r8) hmn(pcols,pver)               ! bulk environmental moist static energy
   real(r8) hsat(pcols,pver)              ! bulk env. saturation moist static energy
   real(r8) hsthat(pcols,pver)            ! hsat at interfaces
   real(r8) wteu(pcols,pver)              ! eu in units of 1/m
   real(r8) wted(pcols,pver)              ! ed in units of 1/m
   real(r8) wtdu(pcols,pver)              ! du in units of 1/m
   real(r8) wtmu(pcols,pver)              ! mu pre-unit conversion
   real(r8) wtmd(pcols,pver)              ! md pre-unit conversion
   real(r8) wtcu(pcols,pver)              ! cu pre-unit conversion
   real(r8) wtevp(pcols,pver)             ! evp pre-unit conversion
   real(r8) rpdpc(pcols,pver)             ! rprd pre-unit conversion
   real(r8) qds(pcols,pver)               ! bulk water saturation mixing ratio in downdraft

   real(r8), target :: Rwt(pcols,pver,wtrc_nwset,2)  !Convtran water tracer ratio

   logical done(pcols,pver)               ! end of convr updraft loop
   real(r8) c0mask(pcols)                 ! autoconversion rates
   integer jd(pcols)                      ! downdraft initation level
   integer lel(pcols)                     ! Equilibrium level (EL)
   integer lcl(pcols)                     ! Lifted condensation level (LCL)
   integer wtpcidx                        ! Physics Buffer index
   integer wtsnidx                        ! Physics Buffer index
!   integer ijctop(pcols)                  ! top deep-con indices
!   integer ijcbot(pcols)                  ! cloud base indices
   logical :: isOk                        ! Used to check water tracer mass conservation

   !---------------------------------------------

   ! physics types
   type(physics_state) :: state1        ! locally modify for evaporation to use, not returned
   type(physics_ptend),target :: ptend_loc     ! package tendencies
   type(physics_state) :: state2        ! debugging

   ! physics buffer fields
   real(r8), pointer, dimension(:)   :: prec         ! total precipitation
   real(r8), pointer, dimension(:)   :: snow         ! snow from ZM convection 
   real(r8), pointer, dimension(:,:) :: cld
   real(r8), pointer, dimension(:,:) :: ql           ! wg grid slice of cloud liquid water.
   real(r8), pointer, dimension(:,:) :: rprd         ! rain production rate
   real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble
   real(r8), pointer, dimension(:,:) :: evapcdp      ! Evaporation of deep convective precipitation
   real(r8), pointer, dimension(:,:) :: flxprec      ! Convective-scale flux of precip at interfaces (kg/m2/s)
   real(r8), pointer, dimension(:,:) :: flxsnow      ! Convective-scale flux of snow   at interfaces (kg/m2/s)
   real(r8), pointer, dimension(:,:) :: dp_cldliq
   real(r8), pointer, dimension(:,:) :: dp_cldice

   !Water tracers
    real(r8), pointer,dimension(:) :: wtprec  !tracer total precipitation
    real(r8), pointer,dimension(:) :: wtsnow

   real(r8) :: jctop(pcols)  ! o row of top-of-deep-convection indices passed out.
   real(r8) :: jcbot(pcols)  ! o row of base of cloud indices passed out.

   real(r8), target :: pcont(pcols), pconb(pcols), freqzm(pcols)

   ! history output fields
   real(r8) :: cape(pcols)        ! w  convective available potential energy.
   real(r8), target :: mu_out(pcols,pver)
   real(r8), target :: md_out(pcols,pver)

   ! used in momentum transport calculation
   real(r8) :: winds(pcols, pver, 2)
   real(r8), target :: wind_tends(pcols, pver, 2)
   real(r8) :: pguall(pcols, pver, 2)
   real(r8) :: pgdall(pcols, pver, 2)
   real(r8) :: icwu(pcols,pver, 2)
   real(r8) :: icwd(pcols,pver, 2)
   real(r8), target :: seten(pcols, pver)
   logical  :: l_windt(2)
   real(r8) :: tfinal1, tfinal2
   integer  :: ii
   integer(c_int64_t), target :: ideep64(pcols), jt64(pcols), maxg64(pcols)
   integer(c_int64_t), target :: liq_type64(wtrc_nwset), ice_type64(wtrc_nwset)
   
   real(r8),pointer :: zm_org2d(:,:)
   real(r8),pointer :: orgt(:,:), org(:,:)

   logical  :: lq(pcnst)

   !----------------------------------------------------------------------

   ! initialize
   lchnk = state%lchnk
   ncol  = state%ncol
   nstep = get_nstep()

   ftem = 0._r8   
   mu_out(:,:) = 0._r8
   md_out(:,:) = 0._r8
   wind_tends(:ncol,:pver,:) = 0.0_r8

   call physics_state_copy(state,state1)             ! copy state to local state1.

   if(trace_water) then
     isOk = wtrc_check_h2o("before-deep", state1, state1%q, ztodt)
   end if

   lq(:) = .FALSE.
   lq(1) = .TRUE.
   if (zmconv_org) then
      lq(ixorg) = .TRUE.
   endif

   !water tracers:
   do ii=1,wtrc_nwset
     lq(wtrc_iatype(ii,iwtvap)) = .TRUE.
   end do

   call physics_ptend_init(ptend_loc, state%psetcols, 'zm_convr', ls=.true., lq=lq)! initialize local ptend type

!
! Associate pointers with physics buffer fields
!
   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, cld_idx,         cld,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

   call pbuf_get_field(pbuf, icwmrdp_idx,     ql )
   call pbuf_get_field(pbuf, rprddp_idx,      rprd )
   call pbuf_get_field(pbuf, fracis_idx,      fracis, start=(/1,1,1/),    kount=(/pcols, pver, pcnst/) )
   call pbuf_get_field(pbuf, nevapr_dpcu_idx, evapcdp )
   call pbuf_get_field(pbuf, prec_dp_idx,     prec )
   call pbuf_get_field(pbuf, snow_dp_idx,     snow )

!
! Begin with Zhang-McFarlane (1996) convection parameterization
!
   call t_startf ('zm_convr')

   if (zmconv_org) then
      allocate(zm_org2d(pcols,pver))
      org => state%q(:,:,ixorg) 
      orgt => ptend_loc%q(:,:,ixorg)
   endif

   call zm_convr(   lchnk   ,ncol    , &
                    state%t       ,state%q(:,:,1),      prec    ,jctop   ,jcbot   , &
                    pblh    ,state%zm      ,state%phis    ,state%zi      ,ptend_loc%q(:,:,1)    , &
                    ptend_loc%s    , state%pmid     ,state%pint    ,state%pdel     , &
                    .5_r8*ztodt    ,mcon    ,cme     , cape,      &
                    tpert   ,dlf     ,pflx    ,zdu     ,rprd    , &
                    mu(:,:,lchnk),md(:,:,lchnk),du(:,:,lchnk),eu(:,:,lchnk),ed(:,:,lchnk)      , &
                    dp(:,:,lchnk) ,dsubcld(:,lchnk) ,jt(:,lchnk),maxg(:,lchnk),ideep(:,lchnk)   , &
                    lengath(lchnk) ,ql      ,rliq  ,landfrac,  &
                    org, orgt, zm_org2d,                       &
                    qu  ,qd  ,dz ,rppe   ,eps0 ,cu , &
                    evp ,tu  ,td  ,jd  ,done  ,lel      ,lcl ,qst ,qsthat     ,&
                    hmn ,hsat    ,hsthat      ,wteu     ,wted,wtdu,wtmu       ,wtmd ,wtcu , &
                    c0mask, rpdpc, qds ,wtevp ) 
                    

   call outfld('CAPE', cape, pcols, lchnk)        ! RBN - CAPE output
!
! Output fractional occurance of ZM convection
!
   do i = 1, lengath(lchnk)
      ideep64(i) = int(ideep(i,lchnk), c_int64_t)
      jt64(i) = int(jt(i,lchnk), c_int64_t)
      maxg64(i) = int(maxg(i,lchnk), c_int64_t)
   end do

   call zm_convr_post_shell(ncol, lengath(lchnk), gravit, cpair, &
        mcon, mu(1,1,lchnk), md(1,1,lchnk), ideep64, jt64, maxg64, &
        ptend_loc%s, state%ps, state%pmid, freqzm, mu_out, md_out, ftem, pcont, pconb)

   call outfld('FREQZM  ',freqzm          ,pcols   ,lchnk   )
!
! Convert mass flux from reported mb/s to kg/m^2/s
!
   call outfld('ZMMU', mu_out(1,1), pcols, lchnk)
   call outfld('ZMMD', md_out(1,1), pcols, lchnk)

   call outfld('ZMDT    ',ftem           ,pcols   ,lchnk   )
   call outfld('ZMDQ    ',ptend_loc%q(1,1,1) ,pcols   ,lchnk   )
   call t_stopf ('zm_convr')

!    do i = 1,pcols
!    do i = 1,nco
    call outfld('PCONVT  ',pcont          ,pcols   ,lchnk   )
    call outfld('PCONVB  ',pconb          ,pcols   ,lchnk   )

  call physics_ptend_init(ptend_all, state%psetcols, 'zm_conv_tend')
!------------------------------------------
!Water tracers
!------------------------------------------

   if (trace_water) then
     call t_startf ('wtrc_zm_convr')

!
! initialize tracer variables
!       wtdlft(:ncol,:,:) = 0._r8
 !      wtcmet(:ncol,:,:) = 0._r8
 !      wtql(:ncol,:,:) = 0._r8
 !      wtqi(:ncol,:,:) = 0._r8
 !      wtqu(:ncol,:,:) = 0._r8
 !      wtqd(:ncol,:,:) = 0._r8
       wtdlf(:ncol,:,:) = 0._r8
       wtrprd(:ncol,:,:) = 0._r8
       wtprect(:ncol,:) = 0._r8
!
! assign totals 
!       wtcmet(:ncol,:,1) = cme(:ncol,:)
!       wtdlft(:ncol,:,1) = dlf(:ncol,:)
!       wtql(:ncol,:,1) = ql(:ncol,:) ! partitioned to liquid on output
!       wtqi(:ncol,:,1) = 0.0         ! partitioned to ice on output
!       wtqu(:ncol,:,1) = qu(:ncol,:)
!       wtqd(:ncol,:,1) = qd(:ncol,:)
       wtrprd(:ncol,:,1) = rprd(:ncol,:)
       wtprect(:ncol,1) = prec(:ncol)

!       isOk = wtrc_check_h2o("before-deep", state, state%q, ztodt,ptend=ptend_loc) !check water tracers before

!New (simplified) version three:

        call wtrc_q1q2_pjr( ptend_loc%q       ,ideep(:,lchnk), lengath(lchnk), &
                    state%q          ,qu      ,done  ,lel ,lcl , &
                    wtdu    ,du(:,:,lchnk)    ,dp(:,:,lchnk)   ,qst ,qsthat  , &
                    hmn     ,hsat  ,hsthat    ,wteu            , &
                    mu(:,:,lchnk)    ,md(:,:,lchnk)  ,qd ,tu  ,td , &
                    dsubcld(:,lchnk) ,maxg(:,lchnk)  ,jt(:,lchnk) ,ncol      , &
                    limcnv-1         ,wted    ,rppe  ,wtmu  ,wtmd , &
                    wtevp   ,evp     ,wtcu    ,cu    ,state%pmid  , &
                    dz      ,qds     ,rpdpc   ,c0mask,eps0        , &
                    jd      ,wtrprd  ,wtdlf   ,ztodt )

     call t_stopf ('wtrc_zm_convr')
   end if
!--------------------------------------------------

  ! This name triggers a special case in physics_types.F90:physics_update()
  call physics_ptend_init(ptend_all, state%psetcols, 'convect_deep')

  ! add tendency from this process to tendencies from other processes
  call physics_ptend_sum(ptend_loc,ptend_all, ncol)

  ! update physics state type state1 with ptend_loc 
  call physics_update(state1, ptend_loc, ztodt)

!-----------------------
!Check water tracer mass
!-----------------------
if (trace_water) then
  isOk = wtrc_check_h2o("after-deep zm", state1, state1%q, ztodt)
end if
!-----------------------

  ! initialize ptend for next process
  lq(:) = .FALSE.
  lq(1) = .TRUE.
  if (zmconv_org) then
     lq(ixorg) = .TRUE.
  endif
 
  !Water tracers:  
  do ii=1,wtrc_nwset
    lq(wtrc_iatype(ii,iwtvap)) = .TRUE.
  end do

  call physics_ptend_init(ptend_loc, state1%psetcols, 'zm_conv_evap', ls=.true., lq=lq)

  call t_startf ('zm_conv_evap')
!
! Determine the phase of the precipitation produced and add latent heat of fusion
! Evaporate some of the precip directly into the environment (Sundqvist)
! Allow this to use the updated state1 and the fresh ptend_loc type
! heating and specific humidity tendencies produced
!

    call pbuf_get_field(pbuf, dp_flxprc_idx, flxprec    )
    call pbuf_get_field(pbuf, dp_flxsnw_idx, flxsnow    )
    call pbuf_get_field(pbuf, dp_cldliq_idx, dp_cldliq  )
    call pbuf_get_field(pbuf, dp_cldice_idx, dp_cldice  )
    call zm_conv_evap_prep_shell(ncol, dp_cldliq, dp_cldice)

    call zm_conv_evap(state1%ncol,state1%lchnk, &
         state1%t,state1%pmid,state1%pdel,state1%q(:pcols,:pver,1), &
         landfrac, &
         ptend_loc%s, tend_s_snwprd, tend_s_snwevmlt, ptend_loc%q(:pcols,:pver,1), &
         rprd, cld, ztodt, prec, snow, &
         evpstore, substore, ntprprd, ntsnprd, flxprec, flxsnow)

    call zm_conv_evap_post_shell(ncol, ztodt, zmconv_org, ixorg, evapcdp, ptend_loc%q, state%q)
    

!----------------------------------------------------------
! Water tracers: evaporate falling precip, and assign state
!----------------------------------------------------------
   if (trace_water) then
     !Check incoming water balance
!     isOk = wtrc_check_h2o("before-evap", state, state%q, ztodt,ptend=ptend_loc)

     !Calculate water tracer rain evaporation (old version):
!     call wtrc_zm_conv_evap(state1  , ptend_loc , wtrprd, cld , ztodt , &
!                            flxprec, flxsnow, fevprec, fevsnow, ntprprd, ntsnprd, &
!                            tend_s_snwprd, tend_s_snwevmlt, wtprect, wtsnowt, .false. )

      call wtrc_precip_evap(state1, wtrprd, ztodt, evpstore, substore, ptend_loc%q, wtprect, wtsnowt)

     !assign values to physics buffer variables
     do i=1,wtrc_nwset
       call pbuf_get_field(pbuf, wtrc_srfpcp_indices(iwtcvrain,i), wtprec)
       call pbuf_get_field(pbuf, wtrc_srfpcp_indices(iwtcvsnow,i), wtsnow)

       wtprec(:) = wtprect(:,wtrc_iatype(i,iwtvap)) - wtsnowt(:,wtrc_iatype(i,iwtvap)) !assign values (should be rain only)
       wtsnow(:) = wtsnowt(:,wtrc_iatype(i,iwtvap))                                    !(snow only)
     end do
   end if
!----------------------------------------------------------

!
! Write out variables from zm_conv_evap
!
   call zm_conv_evap_hist_shell(1, ncol, cpair, ptend_loc%s, tend_s_snwprd, tend_s_snwevmlt, ftem)
   call outfld('EVAPTZM ',ftem           ,pcols   ,lchnk   )
   call zm_conv_evap_hist_shell(2, ncol, cpair, ptend_loc%s, tend_s_snwprd, tend_s_snwevmlt, ftem)
   call outfld('FZSNTZM ',ftem           ,pcols   ,lchnk   )
   call zm_conv_evap_hist_shell(3, ncol, cpair, ptend_loc%s, tend_s_snwprd, tend_s_snwevmlt, ftem)
   call outfld('EVSNTZM ',ftem           ,pcols   ,lchnk   )
   call outfld('EVAPQZM ',ptend_loc%q(1,1,1) ,pcols   ,lchnk   )
   call outfld('ZMFLXPRC', flxprec, pcols, lchnk)
   call outfld('ZMFLXSNW', flxsnow, pcols, lchnk)
   call outfld('ZMNTPRPD', ntprprd, pcols, lchnk)
   call outfld('ZMNTSNPD', ntsnprd, pcols, lchnk)
   call outfld('ZMEIHEAT', ptend_loc%s, pcols, lchnk)
   call outfld('CMFMCDZM   ',mcon ,  pcols   ,lchnk   )
   call outfld('PRECCDZM   ',prec,  pcols   ,lchnk   )


   call t_stopf ('zm_conv_evap')

   call outfld('PRECZ   ', prec   , pcols, lchnk)

  ! add tendency from this process to tend from other processes here
  call physics_ptend_sum(ptend_loc,ptend_all, ncol)

  ! update physics state type state1 with ptend_loc 
  call physics_update(state1, ptend_loc, ztodt)

  !-----------------------
  !Check water tracer mass 
  !-----------------------
  if (trace_water) then
    isOk = wtrc_check_h2o("after-evap", state1, state1%q, ztodt)
  end if
  !-----------------------


  ! Momentum Transport (non-cam3 physics)

  if ( .not. cam_physpkg_is('cam3')) then

     call physics_ptend_init(ptend_loc, state1%psetcols, 'momtran', ls=.true., lu=.true., lv=.true.)

     winds(:ncol,:pver,1) = state1%u(:ncol,:pver)
     winds(:ncol,:pver,2) = state1%v(:ncol,:pver)
   
     l_windt(1) = .true.
     l_windt(2) = .true.

     call t_startf ('momtran')
     call momtran (lchnk, ncol,                                        &
                   l_windt,winds, 2,  mu(1,1,lchnk), md(1,1,lchnk),   &
                   du(1,1,lchnk), eu(1,1,lchnk), ed(1,1,lchnk), dp(1,1,lchnk), dsubcld(1,lchnk),  &
                   jt(1,lchnk),maxg(1,lchnk), ideep(1,lchnk), 1, lengath(lchnk),  &
                   nstep,  wind_tends, pguall, pgdall, icwu, icwd, ztodt, seten )  
     call t_stopf ('momtran')

     call zm_momtran_post_shell(ncol, cpair, wind_tends, seten, ptend_loc%u, ptend_loc%v, ptend_loc%s, ftem)

     call physics_ptend_sum(ptend_loc,ptend_all, ncol)

     ! update physics state type state1 with ptend_loc 
     call physics_update(state1, ptend_loc, ztodt)

     if (zmconv_org) then
        call outfld('ZM_ORG', state%q(:,:,ixorg), pcols, lchnk)
        call outfld('ZM_ORG2D', zm_org2d, pcols, lchnk)
     endif
     call outfld('ZMMTT', ftem             , pcols, lchnk)
     call outfld('ZMMTU', wind_tends(1,1,1), pcols, lchnk)
     call outfld('ZMMTV', wind_tends(1,1,2), pcols, lchnk)
   
     ! Output apparent force from  pressure gradient
     call outfld('ZMUPGU', pguall(1,1,1), pcols, lchnk)
     call outfld('ZMUPGD', pgdall(1,1,1), pcols, lchnk)
     call outfld('ZMVPGU', pguall(1,1,2), pcols, lchnk)
     call outfld('ZMVPGD', pgdall(1,1,2), pcols, lchnk)

     ! Output in-cloud winds
     call outfld('ZMICUU', icwu(1,1,1), pcols, lchnk)
     call outfld('ZMICUD', icwd(1,1,1), pcols, lchnk)
     call outfld('ZMICVU', icwu(1,1,2), pcols, lchnk)
     call outfld('ZMICVD', icwd(1,1,2), pcols, lchnk)

   end if

   ! Transport cloud water and ice only
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)

   lq(:)  = .FALSE.
   lq(2:) = cnst_is_convtran1(2:)

   !transport water tracer condensate (NOTE:  May be better to use is_convtran1 logical
   !instead -JN):
   if(trace_water) then
     do ii=1,wtrc_nwset
       lq(wtrc_iatype(ii,iwtliq)) = .TRUE.
       lq(wtrc_iatype(ii,iwtice)) = .TRUE.
     end do
   end if

   call physics_ptend_init(ptend_loc, state1%psetcols, 'convtran1', lq=lq)

   ! dpdry is not used in this call to convtran since the cloud liquid and ice mixing
   ! ratios are moist
   fake_dpdry(:,:) = 0._r8

   call t_startf ('convtran1')
   call convtran (lchnk,                                        &
                  ptend_loc%lq,state1%q, pcnst,  mu(:,:,lchnk), md(:,:,lchnk),   &
                  du(:,:,lchnk), eu(:,:,lchnk), ed(:,:,lchnk), dp(:,:,lchnk), dsubcld(:,lchnk),  &
                  jt(:,lchnk),maxg(:,lchnk), ideep(:,lchnk), 1, lengath(lchnk),  &
                  nstep,   fracis,  ptend_loc%q, fake_dpdry, Rwt)
 
   !Modify water tendencies via ratios:
   do m = 1,wtrc_nwset
      liq_type64(m) = int(wtrc_iatype(m,iwtliq), c_int64_t)
      ice_type64(m) = int(wtrc_iatype(m,iwtice), c_int64_t)
   end do
   call zm_convtran1_ratio_shell(ncol, wtrc_nwset, Rwt, ptend_loc%q, liq_type64, ice_type64)

   call t_stopf ('convtran1')

   call outfld('ZMDICE ',ptend_loc%q(1,1,ixcldice) ,pcols   ,lchnk   )
   call outfld('ZMDLIQ ',ptend_loc%q(1,1,ixcldliq) ,pcols   ,lchnk   )

   ! add tendency from this process to tend from other processes here
   call physics_ptend_sum(ptend_loc,ptend_all, ncol)

   call physics_state_dealloc(state1)
   call physics_ptend_dealloc(ptend_loc)

   if (zmconv_org) then
      deallocate(zm_org2d)
   end if

end subroutine zm_conv_tend
!=========================================================================================

subroutine zm_conv_append_post_shell_proof(proof_line)

   character(len=*), intent(in) :: proof_line
   character(len=512) :: proof_path
   integer :: status, n, unit_id

   proof_path = ''
   call get_environment_variable('ZM_CONV_POST_SHELL_PROOF_FILE', value=proof_path, length=n, status=status)
   if (status /= 0 .or. n <= 0) return

   open(newunit=unit_id, file=trim(adjustl(proof_path(:n))), status='unknown', action='write', &
        position='append', iostat=status)
   if (status /= 0) return

   write(unit_id,'(A)') trim(proof_line)
   close(unit_id)

end subroutine zm_conv_append_post_shell_proof

!=========================================================================================

subroutine zm_conv_select_post_shell_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (zm_post_shell_selected) return

   impl_name = 'codon'
   call get_environment_variable('ZM_CONV_POST_SHELL_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_zm_post_shell = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_zm_post_shell = .false.
   end if

   zm_post_shell_selected = .true.

   if (masterproc) then
      if (use_native_zm_post_shell) then
         write(iulog,*) 'zm_conv_post_shell implementation = native'
         call zm_conv_append_post_shell_proof('zm_conv_post_shell implementation = native')
      else
         write(iulog,*) 'zm_conv_post_shell implementation = codon'
         call zm_conv_append_post_shell_proof('zm_conv_post_shell implementation = codon')
      end if
      call flush(iulog)
   end if

end subroutine zm_conv_select_post_shell_impl

!=========================================================================================

subroutine zm_convr_post_shell(ncol_local, lengath_local, gravit_local, cpair_local, &
     mcon_local, mu_local, md_local, ideep64_local, jt64_local, maxg64_local, ptend_s_local, &
     state_ps_local, state_pmid_local, freqzm_local, mu_out_local, md_out_local, ftem_local, pcont_local, pconb_local)

   use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

   integer, intent(in) :: ncol_local, lengath_local
   real(r8), intent(in) :: gravit_local, cpair_local
   real(r8), target, intent(inout) :: mcon_local(pcols,pverp)
   real(r8), target, intent(in) :: mu_local(pcols,pver), md_local(pcols,pver)
   integer(c_int64_t), target, intent(in) :: ideep64_local(pcols), jt64_local(pcols), maxg64_local(pcols)
   real(r8), target, intent(in) :: ptend_s_local(pcols,pver), state_ps_local(pcols), state_pmid_local(pcols,pver)
   real(r8), target, intent(inout) :: freqzm_local(pcols), mu_out_local(pcols,pver), md_out_local(pcols,pver)
   real(r8), target, intent(inout) :: ftem_local(pcols,pver), pcont_local(pcols), pconb_local(pcols)

   interface
      subroutine zm_convr_post_shell_codon(ncol_c, pcols_c, pver_c, pverp_c, lengath_c, gravit_c, cpair_c, &
           mcon_p, mu_p, md_p, ideep_p, jt_p, maxg_p, ptend_s_p, state_ps_p, state_pmid_p, freqzm_p, &
           mu_out_p, md_out_p, ftem_p, pcont_p, pconb_p) bind(c, name="zm_convr_post_shell_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pverp_c, lengath_c
         real(c_double), value :: gravit_c, cpair_c
         type(c_ptr), value :: mcon_p, mu_p, md_p, ideep_p, jt_p, maxg_p, ptend_s_p, state_ps_p, state_pmid_p
         type(c_ptr), value :: freqzm_p, mu_out_p, md_out_p, ftem_p, pcont_p, pconb_p
      end subroutine zm_convr_post_shell_codon
   end interface

   call zm_conv_select_post_shell_impl()

   if (use_native_zm_post_shell) then
      call zm_convr_post_shell_native(ncol_local, lengath_local, gravit_local, cpair_local, &
           mcon_local, mu_local, md_local, ideep64_local, jt64_local, maxg64_local, ptend_s_local, &
           state_ps_local, state_pmid_local, freqzm_local, mu_out_local, md_out_local, ftem_local, pcont_local, pconb_local)
      return
   end if

   if (masterproc .and. .not. zm_convr_post_logged) then
      write(iulog,*) 'zm_convr post shell entered (freq/massflux/history/base-top direct = codon)'
      call zm_conv_append_post_shell_proof('zm_convr post shell entered (freq/massflux/history/base-top direct = codon)')
      call flush(iulog)
      zm_convr_post_logged = .true.
   end if

   call zm_convr_post_shell_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
        int(pverp, c_int64_t), int(lengath_local, c_int64_t), real(gravit_local, c_double), real(cpair_local, c_double), &
        c_loc(mcon_local), c_loc(mu_local), c_loc(md_local), c_loc(ideep64_local), c_loc(jt64_local), c_loc(maxg64_local), &
        c_loc(ptend_s_local), c_loc(state_ps_local), c_loc(state_pmid_local), c_loc(freqzm_local), &
        c_loc(mu_out_local), c_loc(md_out_local), c_loc(ftem_local), c_loc(pcont_local), c_loc(pconb_local))

end subroutine zm_convr_post_shell

!=========================================================================================

subroutine zm_convr_post_shell_native(ncol_local, lengath_local, gravit_local, cpair_local, &
     mcon_local, mu_local, md_local, ideep64_local, jt64_local, maxg64_local, ptend_s_local, &
     state_ps_local, state_pmid_local, freqzm_local, mu_out_local, md_out_local, ftem_local, pcont_local, pconb_local)

   use iso_c_binding, only: c_int64_t

   integer, intent(in) :: ncol_local, lengath_local
   real(r8), intent(in) :: gravit_local, cpair_local
   real(r8), intent(inout) :: mcon_local(pcols,pverp)
   real(r8), intent(in) :: mu_local(pcols,pver), md_local(pcols,pver)
   integer(c_int64_t), intent(in) :: ideep64_local(pcols), jt64_local(pcols), maxg64_local(pcols)
   real(r8), intent(in) :: ptend_s_local(pcols,pver), state_ps_local(pcols), state_pmid_local(pcols,pver)
   real(r8), intent(inout) :: freqzm_local(pcols), mu_out_local(pcols,pver), md_out_local(pcols,pver)
   real(r8), intent(inout) :: ftem_local(pcols,pver), pcont_local(pcols), pconb_local(pcols)
   integer :: i, k, ii

   freqzm_local(:) = 0._r8
   do i = 1, lengath_local
      freqzm_local(int(ideep64_local(i))) = 1.0_r8
   end do

   mcon_local(:ncol_local,:pver) = mcon_local(:ncol_local,:pver) * 100._r8/gravit_local

   mu_out_local(:,:) = 0._r8
   md_out_local(:,:) = 0._r8
   do i = 1, lengath_local
      do k = 1, pver
         ii = int(ideep64_local(i))
         mu_out_local(ii,k) = mu_local(i,k) * 100._r8/gravit_local
         md_out_local(ii,k) = md_local(i,k) * 100._r8/gravit_local
      end do
   end do

   ftem_local(:ncol_local,:pver) = ptend_s_local(:ncol_local,:pver)/cpair_local

   pcont_local(:ncol_local) = state_ps_local(:ncol_local)
   pconb_local(:ncol_local) = state_ps_local(:ncol_local)
   do i = 1, lengath_local
      if (maxg64_local(i) .gt. jt64_local(i)) then
         ii = int(ideep64_local(i))
         pcont_local(ii) = state_pmid_local(ii,int(jt64_local(i)))
         pconb_local(ii) = state_pmid_local(ii,int(maxg64_local(i)))
      endif
   end do

end subroutine zm_convr_post_shell_native

!=========================================================================================

subroutine zm_conv_evap_prep_shell(ncol_local, dp_cldliq_local, dp_cldice_local)

   use iso_c_binding, only: c_int64_t, c_loc, c_ptr

   integer, intent(in) :: ncol_local
   real(r8), pointer, intent(inout) :: dp_cldliq_local(:,:), dp_cldice_local(:,:)

   interface
      subroutine zm_conv_evap_prep_shell_codon(ncol_c, pcols_c, pver_c, dp_cldliq_p, dp_cldice_p) &
           bind(c, name="zm_conv_evap_prep_shell_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         type(c_ptr), value :: dp_cldliq_p, dp_cldice_p
      end subroutine zm_conv_evap_prep_shell_codon
   end interface

   call zm_conv_select_post_shell_impl()

   if (use_native_zm_post_shell) then
      dp_cldliq_local(:ncol_local,:) = 0._r8
      dp_cldice_local(:ncol_local,:) = 0._r8
      return
   end if

   if (masterproc .and. .not. zm_conv_evap_prep_logged) then
      write(iulog,*) 'zm_conv_evap prep shell entered (cloud buffers zero direct = codon)'
      call zm_conv_append_post_shell_proof('zm_conv_evap prep shell entered (cloud buffers zero direct = codon)')
      call flush(iulog)
      zm_conv_evap_prep_logged = .true.
   end if

   call zm_conv_evap_prep_shell_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
        c_loc(dp_cldliq_local), c_loc(dp_cldice_local))

end subroutine zm_conv_evap_prep_shell

!=========================================================================================

subroutine zm_conv_evap_post_shell(ncol_local, ztodt_local, do_org_local, ixorg_local, evapcdp_local, ptend_q_local, state_q_local)

   use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

   integer, intent(in) :: ncol_local, ixorg_local
   real(r8), intent(in) :: ztodt_local
   logical, intent(in) :: do_org_local
   real(r8), pointer, intent(inout) :: evapcdp_local(:,:)
   real(r8), target, intent(inout) :: ptend_q_local(pcols,pver,*)
   real(r8), target, intent(in) :: state_q_local(pcols,pver,*)
   integer(c_int64_t) :: do_org_c

   interface
      subroutine zm_conv_evap_post_shell_codon(ncol_c, pcols_c, pver_c, ixorg_c, do_org_c, ztodt_c, &
           evapcdp_p, ptend_q_p, state_q_p) bind(c, name="zm_conv_evap_post_shell_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, ixorg_c, do_org_c
         real(c_double), value :: ztodt_c
         type(c_ptr), value :: evapcdp_p, ptend_q_p, state_q_p
      end subroutine zm_conv_evap_post_shell_codon
   end interface

   call zm_conv_select_post_shell_impl()

   if (use_native_zm_post_shell) then
      evapcdp_local(:ncol_local,:pver) = ptend_q_local(:ncol_local,:pver,1)
      if (do_org_local) then
         ptend_q_local(:ncol_local,:pver,ixorg_local) = &
              min(1._r8,max(0._r8,(50._r8*1000._r8*1000._r8*abs(evapcdp_local(:ncol_local,:pver))) &
              -(state_q_local(:ncol_local,:pver,ixorg_local)/10800._r8)))
         ptend_q_local(:ncol_local,:pver,ixorg_local) = &
              (ptend_q_local(:ncol_local,:pver,ixorg_local) - state_q_local(:ncol_local,:pver,ixorg_local))/ztodt_local
      endif
      return
   end if

   if (do_org_local) then
      do_org_c = 1_c_int64_t
   else
      do_org_c = 0_c_int64_t
   end if

   if (masterproc .and. .not. zm_conv_evap_post_logged) then
      write(iulog,*) 'zm_conv_evap post shell entered (evap/org tendencies direct = codon)'
      call zm_conv_append_post_shell_proof('zm_conv_evap post shell entered (evap/org tendencies direct = codon)')
      call flush(iulog)
      zm_conv_evap_post_logged = .true.
   end if

   call zm_conv_evap_post_shell_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
        int(ixorg_local, c_int64_t), do_org_c, real(ztodt_local, c_double), c_loc(evapcdp_local), &
        c_loc(ptend_q_local), c_loc(state_q_local))

end subroutine zm_conv_evap_post_shell

!=========================================================================================

subroutine zm_conv_evap_hist_shell(mode, ncol_local, cpair_local, ptend_s_local, tend_s_snwprd_local, tend_s_snwevmlt_local, &
     ftem_local)

   use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

   integer, intent(in) :: mode, ncol_local
   real(r8), intent(in) :: cpair_local
   real(r8), target, intent(in) :: ptend_s_local(pcols,pver), tend_s_snwprd_local(pcols,pver), tend_s_snwevmlt_local(pcols,pver)
   real(r8), target, intent(inout) :: ftem_local(pcols,pver)

   interface
      subroutine zm_conv_evap_hist_shell_codon(mode_c, ncol_c, pcols_c, pver_c, cpair_c, ptend_s_p, &
           tend_s_snwprd_p, tend_s_snwevmlt_p, ftem_p) bind(c, name="zm_conv_evap_hist_shell_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: mode_c, ncol_c, pcols_c, pver_c
         real(c_double), value :: cpair_c
         type(c_ptr), value :: ptend_s_p, tend_s_snwprd_p, tend_s_snwevmlt_p, ftem_p
      end subroutine zm_conv_evap_hist_shell_codon
   end interface

   call zm_conv_select_post_shell_impl()

   if (use_native_zm_post_shell) then
      select case (mode)
      case (1)
         ftem_local(:ncol_local,:pver) = ptend_s_local(:ncol_local,:pver)/cpair_local
      case (2)
         ftem_local(:ncol_local,:pver) = tend_s_snwprd_local(:ncol_local,:pver)/cpair_local
      case (3)
         ftem_local(:ncol_local,:pver) = tend_s_snwevmlt_local(:ncol_local,:pver)/cpair_local
      end select
      return
   end if

   if (masterproc .and. .not. zm_conv_evap_hist_logged) then
      write(iulog,*) 'zm_conv_evap hist shell entered (evap heating prep direct = codon)'
      call zm_conv_append_post_shell_proof('zm_conv_evap hist shell entered (evap heating prep direct = codon)')
      call flush(iulog)
      zm_conv_evap_hist_logged = .true.
   end if

   call zm_conv_evap_hist_shell_codon(int(mode, c_int64_t), int(ncol_local, c_int64_t), &
        int(pcols, c_int64_t), int(pver, c_int64_t), real(cpair_local, c_double), c_loc(ptend_s_local), &
        c_loc(tend_s_snwprd_local), c_loc(tend_s_snwevmlt_local), c_loc(ftem_local))

end subroutine zm_conv_evap_hist_shell

!=========================================================================================

subroutine zm_momtran_post_shell(ncol_local, cpair_local, wind_tends_local, seten_local, ptend_u_local, ptend_v_local, &
     ptend_s_local, ftem_local)

   use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

   integer, intent(in) :: ncol_local
   real(r8), intent(in) :: cpair_local
   real(r8), target, intent(in) :: wind_tends_local(pcols,pver,2), seten_local(pcols,pver)
   real(r8), target, intent(inout) :: ptend_u_local(pcols,pver), ptend_v_local(pcols,pver), ptend_s_local(pcols,pver)
   real(r8), target, intent(inout) :: ftem_local(pcols,pver)

   interface
      subroutine zm_momtran_post_shell_codon(ncol_c, pcols_c, pver_c, cpair_c, wind_tends_p, seten_p, &
           ptend_u_p, ptend_v_p, ptend_s_p, ftem_p) bind(c, name="zm_momtran_post_shell_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         real(c_double), value :: cpair_c
         type(c_ptr), value :: wind_tends_p, seten_p, ptend_u_p, ptend_v_p, ptend_s_p, ftem_p
      end subroutine zm_momtran_post_shell_codon
   end interface

   call zm_conv_select_post_shell_impl()

   if (use_native_zm_post_shell) then
      ptend_u_local(:ncol_local,:pver) = wind_tends_local(:ncol_local,:pver,1)
      ptend_v_local(:ncol_local,:pver) = wind_tends_local(:ncol_local,:pver,2)
      ptend_s_local(:ncol_local,:pver) = seten_local(:ncol_local,:pver)
      ftem_local(:ncol_local,:pver) = seten_local(:ncol_local,:pver)/cpair_local
      return
   end if

   if (masterproc .and. .not. zm_momtran_post_logged) then
      write(iulog,*) 'zm_momtran post shell entered (wind/static tendencies direct = codon)'
      call zm_conv_append_post_shell_proof('zm_momtran post shell entered (wind/static tendencies direct = codon)')
      call flush(iulog)
      zm_momtran_post_logged = .true.
   end if

   call zm_momtran_post_shell_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
        real(cpair_local, c_double), c_loc(wind_tends_local), c_loc(seten_local), c_loc(ptend_u_local), &
        c_loc(ptend_v_local), c_loc(ptend_s_local), c_loc(ftem_local))

end subroutine zm_momtran_post_shell

!=========================================================================================

subroutine zm_convtran1_ratio_shell(ncol_local, wtrc_nwset_local, Rwt_local, ptend_q_local, liq_type64_local, ice_type64_local)

   use iso_c_binding, only: c_int64_t, c_loc, c_ptr

   integer, intent(in) :: ncol_local, wtrc_nwset_local
   real(r8), target, intent(in) :: Rwt_local(pcols,pver,wtrc_nwset_local,2)
   real(r8), target, intent(inout) :: ptend_q_local(pcols,pver,*)
   integer(c_int64_t), target, intent(in) :: liq_type64_local(wtrc_nwset_local), ice_type64_local(wtrc_nwset_local)

   interface
      subroutine zm_convtran1_ratio_shell_codon(ncol_c, pcols_c, pver_c, wtrc_nwset_c, Rwt_p, ptend_q_p, &
           liq_type_p, ice_type_p) bind(c, name="zm_convtran1_ratio_shell_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, wtrc_nwset_c
         type(c_ptr), value :: Rwt_p, ptend_q_p, liq_type_p, ice_type_p
      end subroutine zm_convtran1_ratio_shell_codon
   end interface

   call zm_conv_select_post_shell_impl()

   if (use_native_zm_post_shell) then
      call zm_convtran1_ratio_shell_native(ncol_local, wtrc_nwset_local, Rwt_local, ptend_q_local, &
           liq_type64_local, ice_type64_local)
      return
   end if

   if (masterproc .and. .not. zm_convtran1_ratio_logged) then
      write(iulog,*) 'zm_convtran1 ratio shell entered (water tracer ratios direct = codon)'
      call zm_conv_append_post_shell_proof('zm_convtran1 ratio shell entered (water tracer ratios direct = codon)')
      call flush(iulog)
      zm_convtran1_ratio_logged = .true.
   end if

   call zm_convtran1_ratio_shell_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
        int(wtrc_nwset_local, c_int64_t), c_loc(Rwt_local), c_loc(ptend_q_local), c_loc(liq_type64_local), &
        c_loc(ice_type64_local))

end subroutine zm_convtran1_ratio_shell

!=========================================================================================

subroutine zm_convtran1_ratio_shell_native(ncol_local, wtrc_nwset_local, Rwt_local, ptend_q_local, &
     liq_type64_local, ice_type64_local)

   use iso_c_binding, only: c_int64_t

   integer, intent(in) :: ncol_local, wtrc_nwset_local
   real(r8), intent(in) :: Rwt_local(pcols,pver,wtrc_nwset_local,2)
   real(r8), intent(inout) :: ptend_q_local(pcols,pver,*)
   integer(c_int64_t), intent(in) :: liq_type64_local(wtrc_nwset_local), ice_type64_local(wtrc_nwset_local)
   integer :: i, k, m

   do m = 2,wtrc_nwset_local
     do k=1,pver
       do i=1,ncol_local
         ptend_q_local(i,k,int(liq_type64_local(m))) = Rwt_local(i,k,m,1)*ptend_q_local(i,k,int(liq_type64_local(1)))
         ptend_q_local(i,k,int(ice_type64_local(m))) = Rwt_local(i,k,m,2)*ptend_q_local(i,k,int(ice_type64_local(1)))
       end do
     end do
   end do

end subroutine zm_convtran1_ratio_shell_native

!=========================================================================================

subroutine zm_convtran2_dpdry_shell(lengath_local, state_pdeldry_local, ideep64_local, dpdry_local)

   use iso_c_binding, only: c_int64_t, c_loc, c_ptr

   integer, intent(in) :: lengath_local
   real(r8), target, intent(in) :: state_pdeldry_local(pcols,pver)
   integer(c_int64_t), target, intent(in) :: ideep64_local(pcols)
   real(r8), target, intent(inout) :: dpdry_local(pcols,pver)

   interface
      subroutine zm_convtran2_dpdry_shell_codon(pcols_c, pver_c, lengath_c, state_pdeldry_p, ideep_p, dpdry_p) &
           bind(c, name="zm_convtran2_dpdry_shell_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: pcols_c, pver_c, lengath_c
         type(c_ptr), value :: state_pdeldry_p, ideep_p, dpdry_p
      end subroutine zm_convtran2_dpdry_shell_codon
   end interface

   call zm_conv_select_post_shell_impl()

   if (use_native_zm_post_shell) then
      call zm_convtran2_dpdry_shell_native(lengath_local, state_pdeldry_local, ideep64_local, dpdry_local)
      return
   end if

   if (masterproc .and. .not. zm_convtran2_dpdry_logged) then
      write(iulog,*) 'zm_convtran2 dpdry shell entered (dry pressure prep direct = codon)'
      call zm_conv_append_post_shell_proof('zm_convtran2 dpdry shell entered (dry pressure prep direct = codon)')
      call flush(iulog)
      zm_convtran2_dpdry_logged = .true.
   end if

   call zm_convtran2_dpdry_shell_codon(int(pcols, c_int64_t), int(pver, c_int64_t), int(lengath_local, c_int64_t), &
        c_loc(state_pdeldry_local), c_loc(ideep64_local), c_loc(dpdry_local))

end subroutine zm_convtran2_dpdry_shell

!=========================================================================================

subroutine zm_convtran2_dpdry_shell_native(lengath_local, state_pdeldry_local, ideep64_local, dpdry_local)

   use iso_c_binding, only: c_int64_t

   integer, intent(in) :: lengath_local
   real(r8), intent(in) :: state_pdeldry_local(pcols,pver)
   integer(c_int64_t), intent(in) :: ideep64_local(pcols)
   real(r8), intent(inout) :: dpdry_local(pcols,pver)
   integer :: i

   dpdry_local = 0._r8
   do i = 1,lengath_local
      dpdry_local(i,:) = state_pdeldry_local(int(ideep64_local(i)),:)/100._r8
   end do

end subroutine zm_convtran2_dpdry_shell_native

!=========================================================================================

subroutine zm_conv_tend_2( state,  ptend,  ztodt, pbuf)

   use physics_types, only: physics_state, physics_ptend, physics_ptend_init
   use time_manager,  only: get_nstep
   use physics_buffer, only: pbuf_get_index, pbuf_get_field, physics_buffer_desc
   use constituents,  only: pcnst, cnst_get_ind, cnst_is_convtran1
   use error_messages, only: alloc_err
   use water_types,   only: iwtvap, iwtcvsnow 
   use water_tracer_vars, only: trace_water, wtrc_iatype, wtrc_nwset
   use iso_c_binding, only: c_int64_t
 
! Arguments
   type(physics_state), intent(in )   :: state          ! Physics state variables
   type(physics_ptend), intent(out)   :: ptend          ! indivdual parameterization tendencies
   
   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)

! Local variables
   integer :: i, lchnk, istat
   integer :: nstep
   real(r8), target, dimension(pcols,pver) :: dpdry

   !water tracers (not necessary, could make variable optional - JN):
   real(r8) :: Rwt(pcols,pver,wtrc_nwset,2) !water tracer ratio
   integer  :: m                            !loop variable 
   integer(c_int64_t), target :: ideep64(pcols)

! physics buffer fields 
   integer ifld
   real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble
   logical   :: lq(pcnst)

!
! Initialize
!
  lq(:) = .FALSE.
  lq(:) = .not. cnst_is_convtran1(:)

 !transport water tracer condensate (NOTE:  May be better to use is_convtran1 logical
 !instead -JN):
  if(trace_water) then
    do m=iwtvap,iwtcvsnow
      do i=1,wtrc_nwset
        lq(wtrc_iatype(i,m)) = .FALSE.
      end do
    end do
  end if

  call physics_ptend_init(ptend, state%psetcols, 'convtran2', lq=lq )

!
! Associate pointers with physics buffer fields
!
   ifld = pbuf_get_index('FRACIS')
   call pbuf_get_field(pbuf, fracis_idx, fracis, start=(/1,1,1/), kount=(/pcols, pver, pcnst/) )

!
! Transport all constituents except cloud water and ice
!

  lchnk = state%lchnk

   nstep = get_nstep()

   if (any(ptend%lq(:))) then
      ! initialize dpdry for call to convtran
      ! it is used for tracers of dry mixing ratio type
      do i = 1,lengath(lchnk)
         ideep64(i) = int(ideep(i,lchnk), c_int64_t)
      end do
      call zm_convtran2_dpdry_shell(lengath(lchnk), state%pdeldry, ideep64, dpdry)

      call t_startf ('convtran2')
      call convtran (lchnk,                                        &
                     ptend%lq,state%q, pcnst,  mu(:,:,lchnk), md(:,:,lchnk),   &
                     du(:,:,lchnk), eu(:,:,lchnk), ed(:,:,lchnk), dp(:,:,lchnk), dsubcld(:,lchnk),  &
                     jt(:,lchnk),maxg(:,lchnk),ideep(:,lchnk), 1, lengath(lchnk),  &
                     nstep,   fracis,  ptend%q, dpdry, Rwt)
      call t_stopf ('convtran2')
   end if

end subroutine zm_conv_tend_2

!=========================================================================================



end module zm_conv_intr
