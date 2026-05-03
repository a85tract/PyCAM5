module vertical_diffusion

  !----------------------------------------------------------------------------------------------------- !
  ! Module to compute vertical diffusion of momentum,  moisture, trace constituents                      !
  ! and static energy. Separate modules compute                                                          !
  !   1. stresses associated with turbulent flow over orography                                          !
  !      ( turbulent mountain stress )                                                                   !
  !   2. eddy diffusivities, including nonlocal tranport terms                                           !
  !   3. molecular diffusivities                                                                         !
  ! Lastly, a implicit diffusion solver is called, and tendencies retrieved by                           !
  ! differencing the diffused and initial states.                                                        !
  !                                                                                                      !
  ! Calling sequence:                                                                                    !
  !                                                                                                      !
  !  vertical_diffusion_init      Initializes vertical diffustion constants and modules                  !
  !        init_molec_diff        Initializes molecular diffusivity module                               !
  !        init_eddy_diff         Initializes eddy diffusivity module (includes PBL)                     !
  !        init_tms               Initializes turbulent mountain stress module                           !
  !        init_vdiff             Initializes diffusion solver module                                    !
  !  vertical_diffusion_ts_init   Time step initialization (only used for upper boundary condition)      !
  !  vertical_diffusion_tend      Computes vertical diffusion tendencies                                 !
  !        compute_tms            Computes turbulent mountain stresses                                   !
  !        compute_eddy_diff      Computes eddy diffusivities and countergradient terms                  !
  !        compute_vdiff          Solves vertical diffusion equations, including molecular diffusivities !
  !                                                                                                      !
  !---------------------------Code history-------------------------------------------------------------- !
  ! J. Rosinski : Jun. 1992                                                                              !
  ! J. McCaa    : Sep. 2004                                                                              !
  ! S. Park     : Aug. 2006, Dec. 2008. Jan. 2010                                                        ! 
  !----------------------------------------------------------------------------------------------------- !

  use iso_c_binding,    only : c_int64_t
  use shr_kind_mod,     only : r8 => shr_kind_r8, i4=> shr_kind_i4
  use ppgrid,           only : pcols, pver, pverp
  use constituents,     only : pcnst, qmin, cnst_get_ind
  use diffusion_solver, only : vdiff_selector
  use cam_abortutils,   only : endrun
  use error_messages,   only : handle_errmsg
  use physconst,        only :          &
                               cpair  , &     ! Specific heat of dry air
                               gravit , &     ! Acceleration due to gravity
                               rair   , &     ! Gas constant for dry air
                               zvir   , &     ! rh2o/rair - 1
                               latvap , &     ! Latent heat of vaporization
                               latice , &     ! Latent heat of fusion
                               karman , &     ! von Karman constant
                               mwdry  , &     ! Molecular weight of dry air
                               avogad , &     ! Avogadro's number
                               boltz  , &     ! Boltzman's constant
                               tms_orocnst,&  ! turbulent mountain stress parameter
                               tms_z0fac      ! Factor determining z_0 from orographic standard deviation [no unit]
  use cam_history,      only : fieldname_len
  use perf_mod
  use cam_logfile,      only : iulog
  use spmd_utils,       only : masterproc
  use ref_pres,         only : do_molec_diff
  use phys_control,     only : phys_getopts, waccmx_is
  use time_manager,     only : is_first_step

  implicit none
  private      
  save
  
  ! ----------------- !
  ! Public interfaces !
  ! ----------------- !

  public vd_readnl
  public vd_register                                   ! Register multi-time-level variables with physics buffer
  public vertical_diffusion_init                       ! Initialization
  public vertical_diffusion_ts_init                    ! Time step initialization (only used for upper boundary condition)
  public vertical_diffusion_tend                       ! Full vertical diffusion routine

  ! ------------ !
  ! Private data !
  ! ------------ !

  character(len=16)    :: eddy_scheme                  ! Default set in phys_control.F90, use namelist to change
                                                       !     'HB'       = Holtslag and Boville (default)
                                                       !     'HBR'      = Holtslag and Boville and Rash 
                                                       !     'diag_TKE' = Bretherton and Park ( UW Moist Turbulence Scheme )
  integer, parameter   :: nturb = 5                    ! Number of iterations for solution ( when 'diag_TKE' scheme is selected )
  logical, parameter   :: wstarent = .true.            ! Use wstar (.true.) or TKE (.false.) entrainment closure
                                                       ! ( when 'diag_TKE' scheme is selected )
  logical              :: do_pseudocon_diff = .false.  ! If .true., do pseudo-conservative variables diffusion

  character(len=16)    :: shallow_scheme               ! Shallow convection scheme
  character(len=16)    :: microp_scheme                ! Microphysics scheme

  type(vdiff_selector) :: fieldlist_wet                ! Logical switches for moist mixing ratio diffusion
  type(vdiff_selector) :: fieldlist_dry                ! Logical switches for dry mixing ratio diffusion
  type(vdiff_selector) :: fieldlist_molec              ! Logical switches for molecular diffusion
  integer              :: ntop                         ! Top interface level to which vertical diffusion is applied ( = 1 ).
  integer              :: nbot                         ! Bottom interface level to which vertical diffusion is applied ( = pver ).
  integer              :: tke_idx, kvh_idx, kvm_idx    ! TKE and eddy diffusivity indices for fields in the physics buffer
  integer              :: kvt_idx                      ! Index for kinematic molecular conductivity
  integer              :: turbtype_idx, smaw_idx       ! Turbulence type and instability functions
  integer              :: tauresx_idx, tauresy_idx     ! Redisual stress for implicit surface stress

  character(len=fieldname_len) :: vdiffnam(pcnst)      ! Names of vertical diffusion tendencies
  integer              :: ixcldice, ixcldliq           ! Constituent indices for cloud liquid and ice water
  integer              :: ixnumice, ixnumliq


  logical              :: history_amwg                 ! output the variables used by the AMWG diag package
  logical              :: history_eddy                 ! output the eddy variables
  logical              :: history_budget               ! Output tendencies and state variables for CAM4 T, qv, ql, qi
  integer              :: history_budget_histfile_num  ! output history file number for budget fields
  logical              :: history_waccm                ! output variables of interest for WACCM runs

  integer              :: qrl_idx    = 0               ! pbuf index 
  integer              :: wsedl_idx  = 0               ! pbuf index

  integer              :: pblh_idx, tpert_idx, qpert_idx

  ! pbuf fields for unicon
  integer              :: bprod_idx    = -1
  integer              :: ipbl_idx     = -1
  integer              :: kpblh_idx    = -1
  integer              :: wstarPBL_idx = -1
  integer              :: tkes_idx     = -1
  integer              :: went_idx     = -1
  integer              :: qtl_flx_idx  = -1            ! for use in cloud macrophysics when UNICON is on
  integer              :: qti_flx_idx  = -1            ! for use in cloud macrophysics when UNICON is on

  real(r8), parameter  :: unset_r8 = huge(1._r8)
  real(r8)             :: kv_top_pressure              ! Pressure defining the bottom of the upper atmosphere for kvh scaling (Pa)
  real(r8)             :: kv_top_scale                 ! Eddy diffusivity scale factor for upper atmosphere
  real(r8)             :: kv_freetrop_scale            ! Eddy diffusivity scale factor for the free troposphere
  real(r8)             :: eddy_lbulk_max               ! Maximum master length for diag_TKE
  real(r8)             :: eddy_leng_max                ! Maximum dissipation length for diag_TKE
  real(r8)             :: eddy_max_bot_pressure        ! Bottom pressure level (hPa) for eddy_leng_max
  real(r8)             :: eddy_moist_entrain_a2l = unset_r8 ! Moist entrainment enhancement param
  logical              :: diff_cnsrv_mass_check        ! do mass conservation check
  logical              :: do_tms                       ! switch for turbulent mountain stress
  logical              :: do_iss                       ! switch for implicit turbulent surface stress
  logical              :: prog_modal_aero = .false.    ! set true if prognostic modal aerosols are present
  logical              :: use_native_ts_init_impl = .false.
  logical              :: ts_init_impl_selected = .false.
  logical              :: use_native_tend_impl = .false.
  logical              :: tend_impl_selected = .false.
  logical              :: use_native_flux_diag_impl = .false.
  logical              :: flux_diag_impl_selected = .false.
  logical              :: use_native_ptend_core_impl = .false.
  logical              :: ptend_core_impl_selected = .false.
  logical              :: use_native_pre_pbl_diag_impl = .false.
  logical              :: pre_pbl_diag_impl_selected = .false.
  logical              :: use_native_post_pbl_state_impl = .false.
  logical              :: post_pbl_state_impl_selected = .false.
  logical              :: use_native_modal_aero_flux_impl = .false.
  logical              :: modal_aero_flux_impl_selected = .false.
  logical              :: use_native_obklen_diag_impl = .false.
  logical              :: obklen_diag_impl_selected = .false.
  logical              :: use_native_pre_qsat_rh_impl = .false.
  logical              :: pre_qsat_rh_impl_selected = .false.
  logical              :: use_native_post_qsat_diag_impl = .false.
  logical              :: post_qsat_diag_impl_selected = .false.
  logical              :: use_native_diag_batch_impl = .false.
  logical              :: diag_batch_impl_selected = .false.
  logical              :: diag_batch_entered_logged = .false.
  logical              :: use_native_core_batch_impl = .false.
  logical              :: core_batch_impl_selected = .false.
  logical              :: core_batch_entered_logged = .false.
  integer              :: tend_branch_mask = 0
  logical              :: tend_branch_selected = .false.
  integer              :: pmam_ncnst = 0               ! number of prognostic modal aerosol constituents
  integer, allocatable, target :: pmam_cnst_idx(:)     ! constituent indices of prognostic modal aerosols
  integer(c_int64_t), allocatable, target :: pmam_cnst_idx_c(:) ! 64-bit indices for Codon interop

contains

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !
  subroutine vd_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand
    use spmd_utils,      only: masterproc
  
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
  
    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'vd_readnl'
  
    namelist /vert_diff_nl/ kv_top_pressure, kv_top_scale, kv_freetrop_scale, eddy_lbulk_max, eddy_leng_max, &
         eddy_max_bot_pressure, eddy_moist_entrain_a2l, diff_cnsrv_mass_check, do_iss
    !-----------------------------------------------------------------------------
  
    if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'vert_diff_nl', status=ierr)
      if (ierr == 0) then
        read(unitn, vert_diff_nl, iostat=ierr)
        if (ierr /= 0) then
          call endrun(subname // ':: ERROR reading namelist')
        end if
      end if
      close(unitn)
      call freeunit(unitn)
    end if
  
#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(kv_top_pressure,                 1 , mpir8,   0, mpicom)
    call mpibcast(kv_top_scale,                    1 , mpir8,   0, mpicom)
    call mpibcast(kv_freetrop_scale,               1 , mpir8,   0, mpicom)
    call mpibcast(eddy_lbulk_max,                  1 , mpir8,   0, mpicom)
    call mpibcast(eddy_leng_max,                   1 , mpir8,   0, mpicom)
    call mpibcast(eddy_max_bot_pressure,           1 , mpir8,   0, mpicom)
    call mpibcast(eddy_moist_entrain_a2l,          1 , mpir8,   0, mpicom)
    call mpibcast(diff_cnsrv_mass_check,           1 , mpilog,  0, mpicom)
    call mpibcast(do_iss,                          1 , mpilog,  0, mpicom)
#endif

  end subroutine vd_readnl

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vd_register()

    !------------------------------------------------ !
    ! Register physics buffer fields and constituents !
    !------------------------------------------------ !

    use physics_buffer,      only : pbuf_add_field, dtype_r8, dtype_i4

    ! Get eddy_scheme setting from phys_control.F90

    call phys_getopts( eddy_scheme_out          =          eddy_scheme, & 
                       shallow_scheme_out       =       shallow_scheme, &
                       microp_scheme_out        =        microp_scheme, &
                       do_tms_out               =               do_tms)

    ! Add fields to physics buffer

    ! kvt is used by gw_drag.  only needs physpkg scope.
    call pbuf_add_field('kvt', 'physpkg', dtype_r8, (/pcols,pverp/), kvt_idx) 


    call pbuf_add_field('pblh',     'global', dtype_r8, (/pcols/),        pblh_idx)
    call pbuf_add_field('tke',      'global', dtype_r8, (/pcols, pverp/), tke_idx) 
    call pbuf_add_field('kvh',      'global', dtype_r8, (/pcols, pverp/), kvh_idx) 
    call pbuf_add_field('kvm',      'global', dtype_r8, (/pcols, pverp/), kvm_idx ) 
    call pbuf_add_field('turbtype', 'global', dtype_i4, (/pcols, pverp/), turbtype_idx)
    call pbuf_add_field('smaw',     'global', dtype_r8, (/pcols, pverp/), smaw_idx) 

    call pbuf_add_field('tauresx',  'global', dtype_r8, (/pcols/),        tauresx_idx)
    call pbuf_add_field('tauresy',  'global', dtype_r8, (/pcols/),        tauresy_idx)

    call pbuf_add_field('tpert', 'global', dtype_r8, (/pcols/),                       tpert_idx)
    call pbuf_add_field('qpert', 'global', dtype_r8, (/pcols,pcnst/),                 qpert_idx)

    if (trim(shallow_scheme) == 'UNICON') then
       call pbuf_add_field('bprod',    'global', dtype_r8, (/pcols,pverp/), bprod_idx)
       call pbuf_add_field('ipbl',     'global', dtype_i4, (/pcols/),       ipbl_idx)
       call pbuf_add_field('kpblh',    'global', dtype_i4, (/pcols/),       kpblh_idx)
       call pbuf_add_field('wstarPBL', 'global', dtype_r8, (/pcols/),       wstarPBL_idx)
       call pbuf_add_field('tkes',     'global', dtype_r8, (/pcols/),       tkes_idx)
       call pbuf_add_field('went',     'global', dtype_r8, (/pcols/),       went_idx)
       call pbuf_add_field('qtl_flx',  'global', dtype_r8, (/pcols, pverp/), qtl_flx_idx)
       call pbuf_add_field('qti_flx',  'global', dtype_r8, (/pcols, pverp/), qti_flx_idx) 
    end if

  end subroutine vd_register

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_init(pbuf2d)

    !------------------------------------------------------------------!
    ! Initialization of time independent fields for vertical diffusion !
    ! Calls initialization routines for subsidiary modules             !
    !----------------------------------------------------------------- !

    use cam_history,       only : addfld, add_default, phys_decomp
    use eddy_diff,         only : init_eddy_diff
    use hb_diff,           only : init_hb_diff
    use molec_diff,        only : init_molec_diff
    use trb_mtn_stress,    only : init_tms
    use diffusion_solver,  only : init_vdiff, new_fieldlist_vdiff, vdiff_select
    use constituents,      only : cnst_get_ind, cnst_get_type_byind, cnst_name, cnst_get_molec_byind
    use spmd_utils,        only : masterproc
    use ref_pres,          only : ntop_molec, nbot_molec, press_lim_idx, pref_mid
    use physics_buffer,    only : pbuf_set_field, pbuf_get_index, physics_buffer_desc
    use rad_constituents,  only : rad_cnst_get_info, rad_cnst_get_mode_num_idx, &
                                  rad_cnst_get_mam_mmr_idx

    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    character(128) :: errstring   ! Error status for init_vdiff
    integer        :: ntop_eddy   ! Top    interface level to which eddy vertical diffusion is applied ( = 1 )
    integer        :: nbot_eddy   ! Bottom interface level to which eddy vertical diffusion is applied ( = pver )
    integer        :: k           ! Vertical loop index

    real(r8), parameter :: ntop_eddy_pres = 1.e-5_r8 ! Pressure below which eddy diffusion is not done in WACCM-X. (Pa)

    integer :: im, l, m, nmodes, nspec

    ! ----------------------------------------------------------------- !

    if (masterproc) then
       write(iulog,*)'Initializing vertical diffusion (vertical_diffusion_init)'
    end if

    ! ----------------------------------------------------------------- !
    ! Get indices of cloud liquid and ice within the constituents array !
    ! ----------------------------------------------------------------- !

    call cnst_get_ind( 'CLDLIQ', ixcldliq )
    call cnst_get_ind( 'CLDICE', ixcldice )
    if( microp_scheme == 'MG' ) then
        call cnst_get_ind( 'NUMLIQ', ixnumliq )
        call cnst_get_ind( 'NUMICE', ixnumice )
    endif

    ! prog_modal_aero determines whether prognostic modal aerosols are present in the run.
    call phys_getopts(prog_modal_aero_out=prog_modal_aero)
    if (prog_modal_aero) then

       ! Get the constituent indices of the number and mass mixing ratios of the modal 
       ! aerosols.
       !
       ! N.B. - This implementation assumes that the prognostic modal aerosols are 
       !        impacting the climate calculation (i.e., can get info from list 0).
       ! 

       ! First need total number of mam constituents
       call rad_cnst_get_info(0, nmodes=nmodes)
       do m = 1, nmodes
          call rad_cnst_get_info(0, m, nspec=nspec)
          pmam_ncnst = pmam_ncnst + 1 + nspec
       end do

       allocate(pmam_cnst_idx(pmam_ncnst))
       allocate(pmam_cnst_idx_c(pmam_ncnst))

       ! Get the constituent indicies
       im = 1
       do m = 1, nmodes
          call rad_cnst_get_mode_num_idx(m, pmam_cnst_idx(im))
          im = im + 1
          call rad_cnst_get_info(0, m, nspec=nspec)
          do l = 1, nspec
             call rad_cnst_get_mam_mmr_idx(m, l, pmam_cnst_idx(im))
             im = im + 1
          end do
       end do

       pmam_cnst_idx_c(:) = int(pmam_cnst_idx(:), c_int64_t)
    end if

    ! ---------------------------------------------------------------------------------------- !
    ! Initialize molecular diffusivity module                                                  !
    ! Note that computing molecular diffusivities is a trivial expense, but constituent        !
    ! diffusivities depend on their molecular weights. Decomposing the diffusion matric        !
    ! for each constituent is a needless expense unless the diffusivity is significant.        !
    ! ---------------------------------------------------------------------------------------- !

    !----------------------------------------------------------------------------------------
    ! Initialize molecular diffusion and get top and bottom molecular diffusion limits
    !----------------------------------------------------------------------------------------

    if( do_molec_diff ) then
       call init_molec_diff( r8, pcnst, rair, mwdry, avogad, gravit, &
            cpair, boltz, errstring)

       call handle_errmsg(errstring, subname="init_molec_diff")

       call addfld( 'TTPXMLC', 'K/S', 1, 'A', 'Top interf. temp. flux: molec. viscosity', phys_decomp )
       if( masterproc ) write(iulog,fmt='(a,i3,5x,a,i3)') 'NTOP_MOLEC =', ntop_molec, 'NBOT_MOLEC =', nbot_molec
    end if

    ! ---------------------------------- !    
    ! Initialize eddy diffusivity module !
    ! ---------------------------------- !
 
    ! ntop_eddy must be 1 or <= nbot_molec
    ! Currently, it is always 1 except for WACCM-X.
    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
       ntop_eddy  = press_lim_idx(ntop_eddy_pres, top=.true.)
    else
       ntop_eddy = 1
    end if
    nbot_eddy  = pver

    if (masterproc) write(iulog, fmt='(a,i3,5x,a,i3)') 'NTOP_EDDY  =', ntop_eddy, 'NBOT_EDDY  =', nbot_eddy

    select case ( eddy_scheme )
    case ( 'diag_TKE' ) 
        if( masterproc ) write(iulog,*) &
             'vertical_diffusion_init: eddy_diffusivity scheme: UW Moist Turbulence Scheme by Bretherton and Park'
        call init_eddy_diff( r8, pver, gravit, cpair, rair, zvir, latvap, latice, &
                             ntop_eddy, nbot_eddy, karman, eddy_lbulk_max, eddy_leng_max, &
                             eddy_max_bot_pressure, eddy_moist_entrain_a2l)
        if( masterproc ) write(iulog,*) 'vertical_diffusion: nturb, ntop_eddy, nbot_eddy ', nturb, ntop_eddy, nbot_eddy
    case ( 'HB', 'HBR')
        if( masterproc ) write(iulog,*) 'vertical_diffusion_init: eddy_diffusivity scheme:  Holtslag and Boville'
        call init_hb_diff(gravit, cpair, ntop_eddy, nbot_eddy, pref_mid, &
                          karman, eddy_scheme)
        call addfld('HB_ri', 'no',      pver,  'A',  'Richardson Number (HB Scheme), I',  phys_decomp )
    end select
    
    ! The vertical diffusion solver must operate 
    ! over the full range of molecular and eddy diffusion

    ntop = min(ntop_molec,ntop_eddy)
    nbot = max(nbot_molec,nbot_eddy)
    
    ! ------------------------------------------- !
    ! Initialize turbulent mountain stress module !
    ! ------------------------------------------- !

    if( do_tms ) then
       call init_tms( r8, tms_orocnst, tms_z0fac, karman, gravit, rair, errstring )

       call handle_errmsg(errstring, subname="init_tms")

       call addfld( 'TAUTMSX' ,'N/m2  ',  1,  'A',  'Zonal      turbulent mountain surface stress',  phys_decomp )
       call addfld( 'TAUTMSY' ,'N/m2  ',  1,  'A',  'Meridional turbulent mountain surface stress',  phys_decomp )

       if (history_amwg) then
          call add_default('TAUTMSX ', 1, ' ')
          call add_default('TAUTMSY ', 1, ' ')
       end if

       if (masterproc) then
          write(iulog,*)'Using turbulent mountain stress module'
          write(iulog,*)'  tms_orocnst = ',tms_orocnst
          write(iulog,*)'  tms_z0fac = ',tms_z0fac
       end if
    endif
    
    ! ---------------------------------- !
    ! Initialize diffusion solver module !
    ! ---------------------------------- !

    call init_vdiff( r8, iulog, rair, gravit, do_iss, errstring )
    call handle_errmsg(errstring, subname="init_vdiff")

    ! Use fieldlist_wet to select the fields which will be diffused using moist mixing ratios ( all by default )
    ! Use fieldlist_dry to select the fields which will be diffused using dry   mixing ratios.

    fieldlist_wet = new_fieldlist_vdiff( pcnst)
    fieldlist_dry = new_fieldlist_vdiff( pcnst)
    fieldlist_molec = new_fieldlist_vdiff( pcnst)

    if( vdiff_select( fieldlist_wet, 'u' ) .ne. '' ) call endrun( vdiff_select( fieldlist_wet, 'u' ) )
    if( vdiff_select( fieldlist_wet, 'v' ) .ne. '' ) call endrun( vdiff_select( fieldlist_wet, 'v' ) )
    if( vdiff_select( fieldlist_wet, 's' ) .ne. '' ) call endrun( vdiff_select( fieldlist_wet, 's' ) )

    constit_loop: do k = 1, pcnst

       if (prog_modal_aero) then
     ! Do not diffuse droplet number - treated in dropmixnuc
          if (k == ixnumliq) cycle constit_loop
          ! Don't diffuse modal aerosol - treated in dropmixnuc
          do m = 1, pmam_ncnst
             if (k == pmam_cnst_idx(m)) cycle constit_loop
       enddo
       end if

       if( cnst_get_type_byind(k) .eq. 'wet' ) then
          if( vdiff_select( fieldlist_wet, 'q', k ) .ne. '' ) call endrun( vdiff_select( fieldlist_wet, 'q', k ) )
       else
          if( vdiff_select( fieldlist_dry, 'q', k ) .ne. '' ) call endrun( vdiff_select( fieldlist_dry, 'q', k ) )
       endif
    
       ! ----------------------------------------------- !
       ! Select constituents for molecular diffusion     !
       ! ----------------------------------------------- !
       if ( cnst_get_molec_byind(k) .eq. 'minor' ) then
         if( vdiff_select(fieldlist_molec,'q',k) .ne. '' ) call endrun( vdiff_select( fieldlist_molec,'q',k ) )
       endif

    end do constit_loop
    
    ! ------------------------ !
    ! Diagnostic output fields !
    ! ------------------------ !

    do k = 1, pcnst
       vdiffnam(k) = 'VD'//cnst_name(k)
       if( k == 1 ) vdiffnam(k) = 'VD01'    !**** compatibility with old code ****
       call addfld( vdiffnam(k), 'kg/kg/s ', pver, 'A', 'Vertical diffusion of '//cnst_name(k), phys_decomp )
    end do

    call addfld( 'TKE'         , 'm2/s2'  , pverp  , 'A', 'Turbulent Kinetic Energy'                          , phys_decomp )
    call addfld( 'PBLH'        , 'm'      , 1      , 'A', 'PBL height'                                        , phys_decomp )
    call addfld( 'TPERT'       , 'K'      , 1      , 'A', 'Perturbation temperature (eddies in PBL)'          , phys_decomp )
    call addfld( 'QPERT'       , 'kg/kg'  , 1      , 'A', 'Perturbation specific humidity (eddies in PBL)'    , phys_decomp )
    call addfld( 'USTAR'       , 'm/s'    , 1      , 'A', 'Surface friction velocity'                         , phys_decomp )
    call addfld( 'KVH'         , 'm2/s'   , pverp  , 'A', 'Vertical diffusion diffusivities (heat/moisture)'  , phys_decomp )
    call addfld( 'KVM'         , 'm2/s'   , pverp  , 'A', 'Vertical diffusion diffusivities (momentum)'       , phys_decomp )
    call addfld( 'KVT'         , 'm2/s'   , pverp  , 'A', 'Vertical diffusion kinematic molecular conductivity', phys_decomp )
    call addfld( 'CGS'         , 's/m2'   , pverp  , 'A', 'Counter-gradient coeff on surface kinematic fluxes', phys_decomp )
    call addfld( 'DTVKE'       , 'K/s'    , pver   , 'A', 'dT/dt vertical diffusion KE dissipation'           , phys_decomp )
    call addfld( 'DTV'         , 'K/s'    , pver   , 'A', 'T vertical diffusion'                              , phys_decomp )
    call addfld( 'DUV'         , 'm/s2'   , pver   , 'A', 'U vertical diffusion'                              , phys_decomp )
    call addfld( 'DVV'         , 'm/s2'   , pver   , 'A', 'V vertical diffusion'                              , phys_decomp )
    call addfld( 'QT'          , 'kg/kg'  , pver   , 'A', 'Total water mixing ratio'                          , phys_decomp )
    call addfld( 'SL'          , 'J/kg'   , pver   , 'A', 'Liquid water static energy'                        , phys_decomp )
    call addfld( 'SLV'         , 'J/kg'   , pver   , 'A', 'Liq wat virtual static energy'                     , phys_decomp )
    call addfld( 'SLFLX'       , 'W/m2'   , pverp  , 'A', 'Liquid static energy flux'                         , phys_decomp ) 
    call addfld( 'QTFLX'       , 'W/m2'   , pverp  , 'A', 'Total water flux'                                  , phys_decomp ) 
    call addfld( 'UFLX'        , 'W/m2'   , pverp  , 'A', 'Zonal momentum flux'                               , phys_decomp ) 
    call addfld( 'VFLX'        , 'W/m2'   , pverp  , 'A', 'Meridional momentm flux'                           , phys_decomp ) 
    call addfld( 'WGUSTD'      , 'm/s'    , 1      , 'A', 'wind gusts from turbulence'                        , phys_decomp )

    ! ---------------------------------------------------------------------------- !
    ! Below ( with '_PBL') are for detailed analysis of UW Moist Turbulence Scheme !
    ! ---------------------------------------------------------------------------- !

    call addfld( 'qt_pre_PBL  ', 'kg/kg'  , pver   , 'A', 'qt_prePBL'                                         , phys_decomp )
    call addfld( 'sl_pre_PBL  ', 'J/kg'   , pver   , 'A', 'sl_prePBL'                                         , phys_decomp )
    call addfld( 'slv_pre_PBL ', 'J/kg'   , pver   , 'A', 'slv_prePBL'                                        , phys_decomp )
    call addfld( 'u_pre_PBL   ', 'm/s'    , pver   , 'A', 'u_prePBL'                                          , phys_decomp )
    call addfld( 'v_pre_PBL   ', 'm/s'    , pver   , 'A', 'v_prePBL'                                          , phys_decomp )
    call addfld( 'qv_pre_PBL  ', 'kg/kg'  , pver   , 'A', 'qv_prePBL'                                         , phys_decomp )
    call addfld( 'ql_pre_PBL  ', 'kg/kg'  , pver   , 'A', 'ql_prePBL'                                         , phys_decomp )
    call addfld( 'qi_pre_PBL  ', 'kg/kg'  , pver   , 'A', 'qi_prePBL'                                         , phys_decomp )
    call addfld( 't_pre_PBL   ', 'K'      , pver   , 'A', 't_prePBL'                                          , phys_decomp )
    call addfld( 'rh_pre_PBL  ', '%'      , pver   , 'A', 'rh_prePBL'                                         , phys_decomp )

    call addfld( 'qt_aft_PBL  ', 'kg/kg'  , pver   , 'A', 'qt_afterPBL'                                       , phys_decomp )
    call addfld( 'sl_aft_PBL  ', 'J/kg'   , pver   , 'A', 'sl_afterPBL'                                       , phys_decomp )
    call addfld( 'slv_aft_PBL ', 'J/kg'   , pver   , 'A', 'slv_afterPBL'                                      , phys_decomp )
    call addfld( 'u_aft_PBL   ', 'm/s'    , pver   , 'A', 'u_afterPBL'                                        , phys_decomp )
    call addfld( 'v_aft_PBL   ', 'm/s'    , pver   , 'A', 'v_afterPBL'                                        , phys_decomp )
    call addfld( 'qv_aft_PBL  ', 'kg/kg'  , pver   , 'A', 'qv_afterPBL'                                       , phys_decomp )
    call addfld( 'ql_aft_PBL  ', 'kg/kg'  , pver   , 'A', 'ql_afterPBL'                                       , phys_decomp )
    call addfld( 'qi_aft_PBL  ', 'kg/kg'  , pver   , 'A', 'qi_afterPBL'                                       , phys_decomp )
    call addfld( 't_aft_PBL   ', 'K'      , pver   , 'A', 't_afterPBL'                                        , phys_decomp )
    call addfld( 'rh_aft_PBL  ', '%'      , pver   , 'A', 'rh_afterPBL'                                       , phys_decomp )

    call addfld( 'slflx_PBL   ', 'J/m2/s' , pverp  , 'A', 'sl flux by PBL'                                    , phys_decomp ) 
    call addfld( 'qtflx_PBL   ', 'kg/m2/s', pverp  , 'A', 'qt flux by PBL'                                    , phys_decomp ) 
    call addfld( 'uflx_PBL    ', 'kg/m/s2', pverp  , 'A', 'u flux by PBL'                                     , phys_decomp ) 
    call addfld( 'vflx_PBL    ', 'kg/m/s2', pverp  , 'A', 'v flux by PBL'                                     , phys_decomp ) 

    call addfld( 'slflx_cg_PBL', 'J/m2/s' , pverp  , 'A', 'sl_cg flux by PBL'                                 , phys_decomp ) 
    call addfld( 'qtflx_cg_PBL', 'kg/m2/s', pverp  , 'A', 'qt_cg flux by PBL'                                 , phys_decomp ) 
    call addfld( 'uflx_cg_PBL ', 'kg/m/s2', pverp  , 'A', 'u_cg flux by PBL'                                  , phys_decomp ) 
    call addfld( 'vflx_cg_PBL ', 'kg/m/s2', pverp  , 'A', 'v_cg flux by PBL'                                  , phys_decomp ) 

    call addfld( 'qtten_PBL   ', 'kg/kg/s', pver   , 'A', 'qt tendency by PBL'                                , phys_decomp )
    call addfld( 'slten_PBL   ', 'J/kg/s' , pver   , 'A', 'sl tendency by PBL'                                , phys_decomp )
    call addfld( 'uten_PBL    ', 'm/s2'   , pver   , 'A', 'u tendency by PBL'                                 , phys_decomp )
    call addfld( 'vten_PBL    ', 'm/s2'   , pver   , 'A', 'v tendency by PBL'                                 , phys_decomp )
    call addfld( 'qvten_PBL   ', 'kg/kg/s', pver   , 'A', 'qv tendency by PBL'                                , phys_decomp )
    call addfld( 'qlten_PBL   ', 'kg/kg/s', pver   , 'A', 'ql tendency by PBL'                                , phys_decomp )
    call addfld( 'qiten_PBL   ', 'kg/kg/s', pver   , 'A', 'qi tendency by PBL'                                , phys_decomp )
    call addfld( 'tten_PBL    ', 'K/s'    , pver   , 'A', 'T tendency by PBL'                                 , phys_decomp )
    call addfld( 'rhten_PBL   ', '%/s'    , pver   , 'A', 'RH tendency by PBL'                                , phys_decomp )

    call addfld ('ustar',     ' ',1, 'A',' ',phys_decomp)
    call addfld ('obklen',    ' ',1, 'A',' ',phys_decomp)

    if( eddy_scheme .eq. 'diag_TKE' ) then    
       call addfld( 'BPROD   ',  'M2/S3   ',pverp,   'A', 'Buoyancy Production'                               ,phys_decomp)
       call addfld( 'SFI     ',  'FRACTION',pverp,   'A', 'Interface-layer sat frac'                          ,phys_decomp)    
       call addfld( 'SPROD   ',  'M2/S3   ',pverp,   'A', 'Shear Production'                                  ,phys_decomp)   
    endif
 
    ! ----------------------------
    ! determine default variables
    ! ----------------------------
 
    call phys_getopts( history_amwg_out = history_amwg, &
                       history_eddy_out = history_eddy, &
                       history_budget_out = history_budget, &
                       history_budget_histfile_num_out = history_budget_histfile_num, &
                       history_waccm_out = history_waccm)

    if (history_amwg) then
       call add_default(  vdiffnam(1), 1, ' ' )
       call add_default( 'DTV'       , 1, ' ' )  
       call add_default( 'PBLH'      , 1, ' ' )
       if( eddy_scheme .eq. 'diag_TKE' ) then    
          call add_default( 'WGUSTD  ', 1, ' ' )
       endif
    endif
 
    if (history_eddy) then
       if( eddy_scheme .eq. 'diag_TKE' ) then    
          call add_default( 'UFLX    ', 1, ' ' )
          call add_default( 'VFLX    ', 1, ' ' )
       endif
    endif

    if( history_budget ) then
        call add_default( vdiffnam(ixcldliq), history_budget_histfile_num, ' ' )
        call add_default( vdiffnam(ixcldice), history_budget_histfile_num, ' ' )
        if( history_budget_histfile_num > 1 ) then
           call add_default(  vdiffnam(1), history_budget_histfile_num, ' ' )
           call add_default( 'DTV'       , history_budget_histfile_num, ' ' )
        end if
    end if

    if ( history_waccm ) then
       if (do_molec_diff) then
          call add_default ( 'TTPXMLC', 1, ' ' )
       end if
       call add_default( 'DUV'     , 1, ' ' )
       call add_default( 'DVV'     , 1, ' ' )
    end if
     ! ----------------------------
   

     qrl_idx   = pbuf_get_index('QRL')
     wsedl_idx = pbuf_get_index('WSEDL')


     ! Initialization of some pbuf fields
     if (is_first_step()) then
        ! Initialization of pbuf fields tke, kvh, kvm are done in phys_inidat
        call pbuf_set_field(pbuf2d, turbtype_idx, 0    )
        call pbuf_set_field(pbuf2d, smaw_idx,     0.0_r8)
        call pbuf_set_field(pbuf2d, tauresx_idx,  0.0_r8)
        call pbuf_set_field(pbuf2d, tauresy_idx,  0.0_r8)
        if (trim(shallow_scheme) == 'UNICON') then
           call pbuf_set_field(pbuf2d, bprod_idx,    1.0e-5_r8)
           call pbuf_set_field(pbuf2d, ipbl_idx,     0    )
           call pbuf_set_field(pbuf2d, kpblh_idx,    1    )
           call pbuf_set_field(pbuf2d, wstarPBL_idx, 0.0_r8)
           call pbuf_set_field(pbuf2d, tkes_idx,     0.0_r8)
           call pbuf_set_field(pbuf2d, went_idx,     0.0_r8)
           call pbuf_set_field(pbuf2d, qtl_flx_idx,  0.0_r8)
           call pbuf_set_field(pbuf2d, qti_flx_idx,  0.0_r8)   
        end if
     end if

  end subroutine vertical_diffusion_init

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_ts_init( pbuf2d, state )

    !-------------------------------------------------------------- !
    ! Timestep dependent setting,                                   !
    ! At present only invokes upper bc code for molecular diffusion !
    !-------------------------------------------------------------- !
    use molec_diff    , only : init_timestep_molec_diff
    use physics_types , only : physics_state
    use ppgrid        , only : begchunk, endchunk
    
    use physics_buffer, only : physics_buffer_desc

    type(physics_state), intent(in) :: state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    interface
       subroutine vertical_diffusion_ts_init_codon() bind(c, name="vertical_diffusion_ts_init_codon")
       end subroutine vertical_diffusion_ts_init_codon
    end interface

    call vertical_diffusion_ts_init_select_impl()

    if (use_native_ts_init_impl .or. do_molec_diff) then
       call vertical_diffusion_ts_init_native(pbuf2d, state)
       return
    end if

    call vertical_diffusion_ts_init_codon()

  end subroutine vertical_diffusion_ts_init

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_ts_init_native( pbuf2d, state )

    !-------------------------------------------------------------- !
    ! Timestep dependent setting,                                   !
    ! At present only invokes upper bc code for molecular diffusion !
    !-------------------------------------------------------------- !
    use molec_diff    , only : init_timestep_molec_diff
    use physics_types , only : physics_state
    use ppgrid        , only : begchunk, endchunk
    
    use physics_buffer, only : physics_buffer_desc

    type(physics_state), intent(in) :: state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    if (do_molec_diff) call init_timestep_molec_diff(pbuf2d, state )

  end subroutine vertical_diffusion_ts_init_native

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_ts_init_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (ts_init_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('VERTICAL_DIFFUSION_TS_INIT_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_ts_init_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_ts_init_impl = .false.
    end if

    ts_init_impl_selected = .true.

    if (masterproc) then
       if (use_native_ts_init_impl) then
          write(iulog,*) 'vertical_diffusion_ts_init implementation = native'
       else
          write(iulog,*) 'vertical_diffusion_ts_init implementation = codon'
       end if
    end if

  end subroutine vertical_diffusion_ts_init_select_impl

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_tend_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (tend_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('VERTICAL_DIFFUSION_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_tend_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_tend_impl = .false.
    end if

    tend_impl_selected = .true.

    if (masterproc) then
       if (use_native_tend_impl) then
          write(iulog,*) 'vertical_diffusion_tend implementation = native'
       else
          write(iulog,*) 'vertical_diffusion_tend implementation = codon'
       end if
    end if

  end subroutine vertical_diffusion_tend_select_impl

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_tend_select_branches(do_tms_in, do_molec_diff_in, use_diag_tke_in, &
       use_hb_family_in, shallow_unicon_in, prog_modal_aero_in, do_pseudocon_diff_in, &
       diff_cnsrv_mass_check_in, waccmx_special_in)

    use iso_c_binding, only: c_int64_t, c_loc, c_ptr

    logical, intent(in) :: do_tms_in
    logical, intent(in) :: do_molec_diff_in
    logical, intent(in) :: use_diag_tke_in
    logical, intent(in) :: use_hb_family_in
    logical, intent(in) :: shallow_unicon_in
    logical, intent(in) :: prog_modal_aero_in
    logical, intent(in) :: do_pseudocon_diff_in
    logical, intent(in) :: diff_cnsrv_mass_check_in
    logical, intent(in) :: waccmx_special_in

    integer(c_int64_t), target :: branch_mask_c

    interface
       subroutine vertical_diffusion_tend_select_branches_codon(do_tms_c, do_molec_diff_c, use_diag_tke_c, &
            use_hb_family_c, shallow_unicon_c, prog_modal_aero_c, do_pseudocon_diff_c, &
            diff_cnsrv_mass_check_c, waccmx_special_c, branch_mask_p) &
            bind(c, name="vertical_diffusion_tend_select_branches_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: do_tms_c, do_molec_diff_c, use_diag_tke_c, use_hb_family_c
         integer(c_int64_t), value :: shallow_unicon_c, prog_modal_aero_c, do_pseudocon_diff_c
         integer(c_int64_t), value :: diff_cnsrv_mass_check_c, waccmx_special_c
         type(c_ptr), value :: branch_mask_p
       end subroutine vertical_diffusion_tend_select_branches_codon
    end interface

    if (tend_branch_selected) return

    branch_mask_c = 0_c_int64_t
    call vertical_diffusion_tend_select_branches_codon( &
         merge(1_c_int64_t, 0_c_int64_t, do_tms_in), &
         merge(1_c_int64_t, 0_c_int64_t, do_molec_diff_in), &
         merge(1_c_int64_t, 0_c_int64_t, use_diag_tke_in), &
         merge(1_c_int64_t, 0_c_int64_t, use_hb_family_in), &
         merge(1_c_int64_t, 0_c_int64_t, shallow_unicon_in), &
         merge(1_c_int64_t, 0_c_int64_t, prog_modal_aero_in), &
         merge(1_c_int64_t, 0_c_int64_t, do_pseudocon_diff_in), &
         merge(1_c_int64_t, 0_c_int64_t, diff_cnsrv_mass_check_in), &
         merge(1_c_int64_t, 0_c_int64_t, waccmx_special_in), &
         c_loc(branch_mask_c) &
    )

    tend_branch_mask = int(branch_mask_c)
    tend_branch_selected = .true.

  end subroutine vertical_diffusion_tend_select_branches

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_diag_batch_append_proof(proof_line)

    character(len=*), intent(in) :: proof_line
    character(len=512) :: proof_file
    integer :: status, n, unitno

    proof_file = ''
    call get_environment_variable('VERTICAL_DIFFUSION_DIAG_BATCH_PROOF_FILE', value=proof_file, length=n, status=status)
    if (status == 0 .and. n > 0) then
       open(newunit=unitno, file=trim(proof_file(:n)), status='unknown', position='append', action='write')
       write(unitno,'(A)') trim(proof_line)
       close(unitno)
    end if

  end subroutine vertical_diffusion_diag_batch_append_proof

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_diag_batch_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (diag_batch_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('VERTICAL_DIFFUSION_DIAG_BATCH_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_diag_batch_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_diag_batch_impl = .false.
    end if

    diag_batch_impl_selected = .true.

    if (masterproc) then
       if (use_native_diag_batch_impl) then
          write(iulog,*) 'vertical_diffusion_diag_batch implementation = native'
          call vertical_diffusion_diag_batch_append_proof('vertical_diffusion_diag_batch selector entered implementation = native')
       else
          write(iulog,*) 'vertical_diffusion_diag_batch implementation = codon'
          call vertical_diffusion_diag_batch_append_proof('vertical_diffusion_diag_batch selector entered implementation = codon')
       end if
       call flush(iulog)
    end if

  end subroutine vertical_diffusion_diag_batch_select_impl

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_diag_batch_log_entered()

    if (diag_batch_entered_logged) return
    diag_batch_entered_logged = .true.

    if (masterproc) then
       write(iulog,*) 'vertical_diffusion_diag_batch entered (pre/post PBL diagnostics direct = codon)'
       call vertical_diffusion_diag_batch_append_proof('vertical_diffusion_diag_batch entered (pre/post PBL diagnostics direct = codon)')
       call flush(iulog)
    end if

  end subroutine vertical_diffusion_diag_batch_log_entered

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_diag_batch_call(stage, ncol, psetcols_local, rztodt_local, ztodt_local, &
       state_q_local, state_s_local, state_u_local, state_v_local, state_t_local, state_zm_local, &
       q_tmp_local, s_tmp_local, u_tmp_local, v_tmp_local, sl_prePBL_local, qt_prePBL_local, slv_prePBL_local, &
       ftem_local, ftem_prePBL_local, ptend_q_local, ptend_s_local, ptend_u_local, ptend_v_local, &
       qv_aft_PBL_local, ql_aft_PBL_local, qi_aft_PBL_local, s_aft_PBL_local, t_aftPBL_local, &
       u_aft_PBL_local, v_aft_PBL_local, ftem_aftPBL_local, tten_local, rhten_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_null_ptr, c_ptr

    integer, intent(in) :: stage
    integer, intent(in) :: ncol
    integer, optional, intent(in) :: psetcols_local
    real(r8), optional, intent(in) :: rztodt_local, ztodt_local
    real(r8), target, optional, intent(in) :: state_q_local(pcols,pver,pcnst)
    real(r8), target, optional, intent(in) :: state_s_local(pcols,pver)
    real(r8), target, optional, intent(in) :: state_u_local(pcols,pver)
    real(r8), target, optional, intent(in) :: state_v_local(pcols,pver)
    real(r8), target, optional, intent(in) :: state_t_local(pcols,pver)
    real(r8), target, optional, intent(in) :: state_zm_local(pcols,pver)
    real(r8), target, optional, intent(in) :: q_tmp_local(pcols,pver,pcnst)
    real(r8), target, optional, intent(in) :: s_tmp_local(pcols,pver)
    real(r8), target, optional, intent(in) :: u_tmp_local(pcols,pver)
    real(r8), target, optional, intent(in) :: v_tmp_local(pcols,pver)
    real(r8), target, optional, intent(in) :: sl_prePBL_local(pcols,pver)
    real(r8), target, optional, intent(in) :: qt_prePBL_local(pcols,pver)
    real(r8), target, optional, intent(in) :: slv_prePBL_local(pcols,pver)
    real(r8), target, optional, intent(in) :: ftem_local(pcols,pver)
    real(r8), target, optional, intent(in) :: ftem_prePBL_local(pcols,pver)
    real(r8), target, optional, intent(in) :: ptend_q_local(:,:,:)
    real(r8), target, optional, intent(in) :: ptend_s_local(:,:)
    real(r8), target, optional, intent(in) :: ptend_u_local(:,:)
    real(r8), target, optional, intent(in) :: ptend_v_local(:,:)
    real(r8), target, optional, intent(in) :: qv_aft_PBL_local(pcols,pver)
    real(r8), target, optional, intent(in) :: ql_aft_PBL_local(pcols,pver)
    real(r8), target, optional, intent(in) :: qi_aft_PBL_local(pcols,pver)
    real(r8), target, optional, intent(in) :: s_aft_PBL_local(pcols,pver)
    real(r8), target, optional, intent(in) :: t_aftPBL_local(pcols,pver)
    real(r8), target, optional, intent(in) :: u_aft_PBL_local(pcols,pver)
    real(r8), target, optional, intent(in) :: v_aft_PBL_local(pcols,pver)
    real(r8), target, optional, intent(in) :: ftem_aftPBL_local(pcols,pver)
    real(r8), target, optional, intent(in) :: tten_local(pcols,pver)
    real(r8), target, optional, intent(in) :: rhten_local(pcols,pver)

    integer(c_int64_t) :: psetcols_c
    real(c_double) :: rztodt_c, ztodt_c
    type(c_ptr) :: state_q_p, state_s_p, state_u_p, state_v_p, state_t_p, state_zm_p
    type(c_ptr) :: q_tmp_p, s_tmp_p, u_tmp_p, v_tmp_p, sl_prePBL_p, qt_prePBL_p, slv_prePBL_p
    type(c_ptr) :: ftem_p, ftem_prePBL_p, ptend_q_p, ptend_s_p, ptend_u_p, ptend_v_p
    type(c_ptr) :: qv_aft_PBL_p, ql_aft_PBL_p, qi_aft_PBL_p, s_aft_PBL_p, t_aftPBL_p
    type(c_ptr) :: u_aft_PBL_p, v_aft_PBL_p, ftem_aftPBL_p, tten_p, rhten_p

    interface
       subroutine vertical_diffusion_diag_batch_codon(stage_c, ncol_c, pcols_c, pver_c, pcnst_c, psetcols_c, &
            ixcldliq_c, ixcldice_c, latvap_c, latice_c, zvir_c, rztodt_c, ztodt_c, gravit_c, cpair_c, &
            state_q_p, state_s_p, state_u_p, state_v_p, state_t_p, state_zm_p, q_tmp_p, s_tmp_p, u_tmp_p, v_tmp_p, &
            sl_prePBL_p, qt_prePBL_p, slv_prePBL_p, ftem_p, ftem_prePBL_p, ptend_q_p, ptend_s_p, ptend_u_p, &
            ptend_v_p, qv_aft_PBL_p, ql_aft_PBL_p, qi_aft_PBL_p, s_aft_PBL_p, t_aftPBL_p, u_aft_PBL_p, &
            v_aft_PBL_p, ftem_aftPBL_p, tten_p, rhten_p) bind(c, name="vertical_diffusion_diag_batch_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: stage_c, ncol_c, pcols_c, pver_c, pcnst_c, psetcols_c, ixcldliq_c, ixcldice_c
         real(c_double), value :: latvap_c, latice_c, zvir_c, rztodt_c, ztodt_c, gravit_c, cpair_c
         type(c_ptr), value :: state_q_p, state_s_p, state_u_p, state_v_p, state_t_p, state_zm_p
         type(c_ptr), value :: q_tmp_p, s_tmp_p, u_tmp_p, v_tmp_p, sl_prePBL_p, qt_prePBL_p, slv_prePBL_p
         type(c_ptr), value :: ftem_p, ftem_prePBL_p, ptend_q_p, ptend_s_p, ptend_u_p, ptend_v_p
         type(c_ptr), value :: qv_aft_PBL_p, ql_aft_PBL_p, qi_aft_PBL_p, s_aft_PBL_p, t_aftPBL_p
         type(c_ptr), value :: u_aft_PBL_p, v_aft_PBL_p, ftem_aftPBL_p, tten_p, rhten_p
       end subroutine vertical_diffusion_diag_batch_codon
    end interface

    psetcols_c = int(pcols, c_int64_t)
    if (present(psetcols_local)) psetcols_c = int(psetcols_local, c_int64_t)
    rztodt_c = 0._c_double
    ztodt_c = 0._c_double
    if (present(rztodt_local)) rztodt_c = real(rztodt_local, c_double)
    if (present(ztodt_local)) ztodt_c = real(ztodt_local, c_double)

    state_q_p = c_null_ptr; if (present(state_q_local)) state_q_p = c_loc(state_q_local)
    state_s_p = c_null_ptr; if (present(state_s_local)) state_s_p = c_loc(state_s_local)
    state_u_p = c_null_ptr; if (present(state_u_local)) state_u_p = c_loc(state_u_local)
    state_v_p = c_null_ptr; if (present(state_v_local)) state_v_p = c_loc(state_v_local)
    state_t_p = c_null_ptr; if (present(state_t_local)) state_t_p = c_loc(state_t_local)
    state_zm_p = c_null_ptr; if (present(state_zm_local)) state_zm_p = c_loc(state_zm_local)
    q_tmp_p = c_null_ptr; if (present(q_tmp_local)) q_tmp_p = c_loc(q_tmp_local)
    s_tmp_p = c_null_ptr; if (present(s_tmp_local)) s_tmp_p = c_loc(s_tmp_local)
    u_tmp_p = c_null_ptr; if (present(u_tmp_local)) u_tmp_p = c_loc(u_tmp_local)
    v_tmp_p = c_null_ptr; if (present(v_tmp_local)) v_tmp_p = c_loc(v_tmp_local)
    sl_prePBL_p = c_null_ptr; if (present(sl_prePBL_local)) sl_prePBL_p = c_loc(sl_prePBL_local)
    qt_prePBL_p = c_null_ptr; if (present(qt_prePBL_local)) qt_prePBL_p = c_loc(qt_prePBL_local)
    slv_prePBL_p = c_null_ptr; if (present(slv_prePBL_local)) slv_prePBL_p = c_loc(slv_prePBL_local)
    ftem_p = c_null_ptr; if (present(ftem_local)) ftem_p = c_loc(ftem_local)
    ftem_prePBL_p = c_null_ptr; if (present(ftem_prePBL_local)) ftem_prePBL_p = c_loc(ftem_prePBL_local)
    ptend_q_p = c_null_ptr; if (present(ptend_q_local)) ptend_q_p = c_loc(ptend_q_local)
    ptend_s_p = c_null_ptr; if (present(ptend_s_local)) ptend_s_p = c_loc(ptend_s_local)
    ptend_u_p = c_null_ptr; if (present(ptend_u_local)) ptend_u_p = c_loc(ptend_u_local)
    ptend_v_p = c_null_ptr; if (present(ptend_v_local)) ptend_v_p = c_loc(ptend_v_local)
    qv_aft_PBL_p = c_null_ptr; if (present(qv_aft_PBL_local)) qv_aft_PBL_p = c_loc(qv_aft_PBL_local)
    ql_aft_PBL_p = c_null_ptr; if (present(ql_aft_PBL_local)) ql_aft_PBL_p = c_loc(ql_aft_PBL_local)
    qi_aft_PBL_p = c_null_ptr; if (present(qi_aft_PBL_local)) qi_aft_PBL_p = c_loc(qi_aft_PBL_local)
    s_aft_PBL_p = c_null_ptr; if (present(s_aft_PBL_local)) s_aft_PBL_p = c_loc(s_aft_PBL_local)
    t_aftPBL_p = c_null_ptr; if (present(t_aftPBL_local)) t_aftPBL_p = c_loc(t_aftPBL_local)
    u_aft_PBL_p = c_null_ptr; if (present(u_aft_PBL_local)) u_aft_PBL_p = c_loc(u_aft_PBL_local)
    v_aft_PBL_p = c_null_ptr; if (present(v_aft_PBL_local)) v_aft_PBL_p = c_loc(v_aft_PBL_local)
    ftem_aftPBL_p = c_null_ptr; if (present(ftem_aftPBL_local)) ftem_aftPBL_p = c_loc(ftem_aftPBL_local)
    tten_p = c_null_ptr; if (present(tten_local)) tten_p = c_loc(tten_local)
    rhten_p = c_null_ptr; if (present(rhten_local)) rhten_p = c_loc(rhten_local)

    call vertical_diffusion_diag_batch_log_entered()
    call vertical_diffusion_diag_batch_codon( &
         int(stage, c_int64_t), int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
         int(pcnst, c_int64_t), psetcols_c, int(ixcldliq, c_int64_t), int(ixcldice, c_int64_t), &
         real(latvap, c_double), real(latice, c_double), real(zvir, c_double), rztodt_c, ztodt_c, &
         real(gravit, c_double), real(cpair, c_double), state_q_p, state_s_p, state_u_p, state_v_p, state_t_p, &
         state_zm_p, q_tmp_p, s_tmp_p, u_tmp_p, v_tmp_p, sl_prePBL_p, qt_prePBL_p, slv_prePBL_p, ftem_p, &
         ftem_prePBL_p, ptend_q_p, ptend_s_p, ptend_u_p, ptend_v_p, qv_aft_PBL_p, ql_aft_PBL_p, qi_aft_PBL_p, &
         s_aft_PBL_p, t_aftPBL_p, u_aft_PBL_p, v_aft_PBL_p, ftem_aftPBL_p, tten_p, rhten_p)

  end subroutine vertical_diffusion_diag_batch_call

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_core_batch_append_proof(proof_line)

    character(len=*), intent(in) :: proof_line
    character(len=512) :: proof_file
    integer :: status, n, unitno

    proof_file = ''
    call get_environment_variable('VERTICAL_DIFFUSION_CORE_BATCH_PROOF_FILE', value=proof_file, length=n, status=status)
    if (status == 0 .and. n > 0) then
       open(newunit=unitno, file=trim(proof_file(:n)), status='unknown', position='append', action='write')
       write(unitno,'(A)') trim(proof_line)
       close(unitno)
    end if

  end subroutine vertical_diffusion_core_batch_append_proof

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_core_batch_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (core_batch_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('VERTICAL_DIFFUSION_CORE_BATCH_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_core_batch_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_core_batch_impl = .false.
    end if

    core_batch_impl_selected = .true.

    if (masterproc) then
       if (use_native_core_batch_impl) then
          write(iulog,*) 'vertical_diffusion_core_batch implementation = native'
          call vertical_diffusion_core_batch_append_proof('vertical_diffusion_core_batch selector entered implementation = native')
       else
          write(iulog,*) 'vertical_diffusion_core_batch implementation = codon'
          call vertical_diffusion_core_batch_append_proof('vertical_diffusion_core_batch selector entered implementation = codon')
       end if
       call flush(iulog)
    end if

  end subroutine vertical_diffusion_core_batch_select_impl

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_core_batch_log_entered()

    if (core_batch_entered_logged) return
    core_batch_entered_logged = .true.

    if (masterproc) then
       write(iulog,*) 'vertical_diffusion_core_batch entered (modal/flux/ptend direct = codon)'
       call vertical_diffusion_core_batch_append_proof('vertical_diffusion_core_batch entered (modal/flux/ptend direct = codon)')
       call flush(iulog)
    end if

  end subroutine vertical_diffusion_core_batch_log_entered

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_core_batch_call(stage, ncol, psetcols_local, rztodt_local, ztodt_local, &
       q_tmp_local, s_tmp_local, u_tmp_local, v_tmp_local, state_q_local, state_s_local, state_u_local, &
       state_v_local, state_rpdel_local, pint_local, zi_local, zm_local, cflx_local, kvh_local, kvm_local, &
       cgs_local, cgh_local, shflx_local, tautotx_local, tautoty_local, sl_local, qt_local, slv_local, &
       sl_prePBL_local, qt_prePBL_local, slflx_local, qtflx_local, uflx_local, vflx_local, slflx_cg_local, &
       qtflx_cg_local, uflx_cg_local, vflx_cg_local, ptend_q_local, ptend_s_local, ptend_u_local, &
       ptend_v_local, slten_local, qtten_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_null_ptr, c_ptr

    integer, intent(in) :: stage
    integer, intent(in) :: ncol
    integer, optional, intent(in) :: psetcols_local
    real(r8), optional, intent(in) :: rztodt_local, ztodt_local
    real(r8), target, optional, intent(in) :: q_tmp_local(pcols,pver,pcnst)
    real(r8), target, optional, intent(in) :: s_tmp_local(pcols,pver)
    real(r8), target, optional, intent(in) :: u_tmp_local(pcols,pver)
    real(r8), target, optional, intent(in) :: v_tmp_local(pcols,pver)
    real(r8), target, optional, intent(in) :: state_q_local(pcols,pver,pcnst)
    real(r8), target, optional, intent(in) :: state_s_local(pcols,pver)
    real(r8), target, optional, intent(in) :: state_u_local(pcols,pver)
    real(r8), target, optional, intent(in) :: state_v_local(pcols,pver)
    real(r8), target, optional, intent(in) :: state_rpdel_local(pcols,pver)
    real(r8), target, optional, intent(in) :: pint_local(pcols,pverp)
    real(r8), target, optional, intent(in) :: zi_local(pcols,pverp)
    real(r8), target, optional, intent(in) :: zm_local(pcols,pver)
    real(r8), target, optional, intent(in) :: cflx_local(pcols,pcnst)
    real(r8), target, optional, intent(in) :: kvh_local(pcols,pverp)
    real(r8), target, optional, intent(in) :: kvm_local(pcols,pverp)
    real(r8), target, optional, intent(in) :: cgs_local(pcols,pverp)
    real(r8), target, optional, intent(in) :: cgh_local(pcols,pverp)
    real(r8), target, optional, intent(in) :: shflx_local(pcols)
    real(r8), target, optional, intent(in) :: tautotx_local(pcols)
    real(r8), target, optional, intent(in) :: tautoty_local(pcols)
    real(r8), target, optional, intent(in) :: sl_local(pcols,pver)
    real(r8), target, optional, intent(in) :: qt_local(pcols,pver)
    real(r8), target, optional, intent(in) :: slv_local(pcols,pver)
    real(r8), target, optional, intent(in) :: sl_prePBL_local(pcols,pver)
    real(r8), target, optional, intent(in) :: qt_prePBL_local(pcols,pver)
    real(r8), target, optional, intent(in) :: slflx_local(pcols,pverp)
    real(r8), target, optional, intent(in) :: qtflx_local(pcols,pverp)
    real(r8), target, optional, intent(in) :: uflx_local(pcols,pverp)
    real(r8), target, optional, intent(in) :: vflx_local(pcols,pverp)
    real(r8), target, optional, intent(in) :: slflx_cg_local(pcols,pverp)
    real(r8), target, optional, intent(in) :: qtflx_cg_local(pcols,pverp)
    real(r8), target, optional, intent(in) :: uflx_cg_local(pcols,pverp)
    real(r8), target, optional, intent(in) :: vflx_cg_local(pcols,pverp)
    real(r8), target, optional, intent(in) :: ptend_q_local(:,:,:)
    real(r8), target, optional, intent(in) :: ptend_s_local(:,:)
    real(r8), target, optional, intent(in) :: ptend_u_local(:,:)
    real(r8), target, optional, intent(in) :: ptend_v_local(:,:)
    real(r8), target, optional, intent(in) :: slten_local(pcols,pver)
    real(r8), target, optional, intent(in) :: qtten_local(pcols,pver)

    integer(c_int64_t) :: psetcols_c
    real(c_double) :: rztodt_c, ztodt_c
    type(c_ptr) :: q_tmp_p, s_tmp_p, u_tmp_p, v_tmp_p, state_q_p, state_s_p, state_u_p, state_v_p
    type(c_ptr) :: state_rpdel_p, pint_p, zi_p, zm_p, cflx_p, kvh_p, kvm_p, cgs_p, cgh_p
    type(c_ptr) :: shflx_p, tautotx_p, tautoty_p, sl_p, qt_p, slv_p, sl_prePBL_p, qt_prePBL_p
    type(c_ptr) :: slflx_p, qtflx_p, uflx_p, vflx_p, slflx_cg_p, qtflx_cg_p, uflx_cg_p, vflx_cg_p
    type(c_ptr) :: ptend_q_p, ptend_s_p, ptend_u_p, ptend_v_p, slten_p, qtten_p, pmam_cnst_idx_p

    interface
       subroutine vertical_diffusion_core_batch_codon(stage_c, ncol_c, pcols_c, pver_c, pverp_c, pcnst_c, &
            psetcols_c, pmam_ncnst_c, ixcldliq_c, ixcldice_c, latvap_c, latice_c, zvir_c, rair_c, gravit_c, &
            cpair_c, rztodt_c, ztodt_c, q_tmp_p, s_tmp_p, u_tmp_p, v_tmp_p, state_q_p, state_s_p, state_u_p, &
            state_v_p, state_rpdel_p, pint_p, zi_p, zm_p, cflx_p, kvh_p, kvm_p, cgs_p, cgh_p, shflx_p, &
            tautotx_p, tautoty_p, sl_p, qt_p, slv_p, sl_prePBL_p, qt_prePBL_p, slflx_p, qtflx_p, uflx_p, &
            vflx_p, slflx_cg_p, qtflx_cg_p, uflx_cg_p, vflx_cg_p, pmam_cnst_idx_p, ptend_q_p, ptend_s_p, &
            ptend_u_p, ptend_v_p, slten_p, qtten_p) bind(c, name="vertical_diffusion_core_batch_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: stage_c, ncol_c, pcols_c, pver_c, pverp_c, pcnst_c, psetcols_c, pmam_ncnst_c
         integer(c_int64_t), value :: ixcldliq_c, ixcldice_c
         real(c_double), value :: latvap_c, latice_c, zvir_c, rair_c, gravit_c, cpair_c, rztodt_c, ztodt_c
         type(c_ptr), value :: q_tmp_p, s_tmp_p, u_tmp_p, v_tmp_p, state_q_p, state_s_p, state_u_p, state_v_p
         type(c_ptr), value :: state_rpdel_p, pint_p, zi_p, zm_p, cflx_p, kvh_p, kvm_p, cgs_p, cgh_p
         type(c_ptr), value :: shflx_p, tautotx_p, tautoty_p, sl_p, qt_p, slv_p, sl_prePBL_p, qt_prePBL_p
         type(c_ptr), value :: slflx_p, qtflx_p, uflx_p, vflx_p, slflx_cg_p, qtflx_cg_p, uflx_cg_p, vflx_cg_p
         type(c_ptr), value :: pmam_cnst_idx_p, ptend_q_p, ptend_s_p, ptend_u_p, ptend_v_p, slten_p, qtten_p
       end subroutine vertical_diffusion_core_batch_codon
    end interface

    psetcols_c = int(pcols, c_int64_t)
    if (present(psetcols_local)) psetcols_c = int(psetcols_local, c_int64_t)
    rztodt_c = 0._c_double
    ztodt_c = 0._c_double
    if (present(rztodt_local)) rztodt_c = real(rztodt_local, c_double)
    if (present(ztodt_local)) ztodt_c = real(ztodt_local, c_double)

    q_tmp_p = c_null_ptr; if (present(q_tmp_local)) q_tmp_p = c_loc(q_tmp_local)
    s_tmp_p = c_null_ptr; if (present(s_tmp_local)) s_tmp_p = c_loc(s_tmp_local)
    u_tmp_p = c_null_ptr; if (present(u_tmp_local)) u_tmp_p = c_loc(u_tmp_local)
    v_tmp_p = c_null_ptr; if (present(v_tmp_local)) v_tmp_p = c_loc(v_tmp_local)
    state_q_p = c_null_ptr; if (present(state_q_local)) state_q_p = c_loc(state_q_local)
    state_s_p = c_null_ptr; if (present(state_s_local)) state_s_p = c_loc(state_s_local)
    state_u_p = c_null_ptr; if (present(state_u_local)) state_u_p = c_loc(state_u_local)
    state_v_p = c_null_ptr; if (present(state_v_local)) state_v_p = c_loc(state_v_local)
    state_rpdel_p = c_null_ptr; if (present(state_rpdel_local)) state_rpdel_p = c_loc(state_rpdel_local)
    pint_p = c_null_ptr; if (present(pint_local)) pint_p = c_loc(pint_local)
    zi_p = c_null_ptr; if (present(zi_local)) zi_p = c_loc(zi_local)
    zm_p = c_null_ptr; if (present(zm_local)) zm_p = c_loc(zm_local)
    cflx_p = c_null_ptr; if (present(cflx_local)) cflx_p = c_loc(cflx_local)
    kvh_p = c_null_ptr; if (present(kvh_local)) kvh_p = c_loc(kvh_local)
    kvm_p = c_null_ptr; if (present(kvm_local)) kvm_p = c_loc(kvm_local)
    cgs_p = c_null_ptr; if (present(cgs_local)) cgs_p = c_loc(cgs_local)
    cgh_p = c_null_ptr; if (present(cgh_local)) cgh_p = c_loc(cgh_local)
    shflx_p = c_null_ptr; if (present(shflx_local)) shflx_p = c_loc(shflx_local)
    tautotx_p = c_null_ptr; if (present(tautotx_local)) tautotx_p = c_loc(tautotx_local)
    tautoty_p = c_null_ptr; if (present(tautoty_local)) tautoty_p = c_loc(tautoty_local)
    sl_p = c_null_ptr; if (present(sl_local)) sl_p = c_loc(sl_local)
    qt_p = c_null_ptr; if (present(qt_local)) qt_p = c_loc(qt_local)
    slv_p = c_null_ptr; if (present(slv_local)) slv_p = c_loc(slv_local)
    sl_prePBL_p = c_null_ptr; if (present(sl_prePBL_local)) sl_prePBL_p = c_loc(sl_prePBL_local)
    qt_prePBL_p = c_null_ptr; if (present(qt_prePBL_local)) qt_prePBL_p = c_loc(qt_prePBL_local)
    slflx_p = c_null_ptr; if (present(slflx_local)) slflx_p = c_loc(slflx_local)
    qtflx_p = c_null_ptr; if (present(qtflx_local)) qtflx_p = c_loc(qtflx_local)
    uflx_p = c_null_ptr; if (present(uflx_local)) uflx_p = c_loc(uflx_local)
    vflx_p = c_null_ptr; if (present(vflx_local)) vflx_p = c_loc(vflx_local)
    slflx_cg_p = c_null_ptr; if (present(slflx_cg_local)) slflx_cg_p = c_loc(slflx_cg_local)
    qtflx_cg_p = c_null_ptr; if (present(qtflx_cg_local)) qtflx_cg_p = c_loc(qtflx_cg_local)
    uflx_cg_p = c_null_ptr; if (present(uflx_cg_local)) uflx_cg_p = c_loc(uflx_cg_local)
    vflx_cg_p = c_null_ptr; if (present(vflx_cg_local)) vflx_cg_p = c_loc(vflx_cg_local)
    ptend_q_p = c_null_ptr; if (present(ptend_q_local)) ptend_q_p = c_loc(ptend_q_local)
    ptend_s_p = c_null_ptr; if (present(ptend_s_local)) ptend_s_p = c_loc(ptend_s_local)
    ptend_u_p = c_null_ptr; if (present(ptend_u_local)) ptend_u_p = c_loc(ptend_u_local)
    ptend_v_p = c_null_ptr; if (present(ptend_v_local)) ptend_v_p = c_loc(ptend_v_local)
    slten_p = c_null_ptr; if (present(slten_local)) slten_p = c_loc(slten_local)
    qtten_p = c_null_ptr; if (present(qtten_local)) qtten_p = c_loc(qtten_local)
    pmam_cnst_idx_p = c_null_ptr
    if (allocated(pmam_cnst_idx_c)) pmam_cnst_idx_p = c_loc(pmam_cnst_idx_c)

    call vertical_diffusion_core_batch_log_entered()
    call vertical_diffusion_core_batch_codon( &
         int(stage, c_int64_t), int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
         int(pverp, c_int64_t), int(pcnst, c_int64_t), psetcols_c, int(pmam_ncnst, c_int64_t), &
         int(ixcldliq, c_int64_t), int(ixcldice, c_int64_t), real(latvap, c_double), real(latice, c_double), &
         real(zvir, c_double), real(rair, c_double), real(gravit, c_double), real(cpair, c_double), &
         rztodt_c, ztodt_c, q_tmp_p, s_tmp_p, u_tmp_p, v_tmp_p, state_q_p, state_s_p, state_u_p, state_v_p, &
         state_rpdel_p, pint_p, zi_p, zm_p, cflx_p, kvh_p, kvm_p, cgs_p, cgh_p, shflx_p, tautotx_p, tautoty_p, &
         sl_p, qt_p, slv_p, sl_prePBL_p, qt_prePBL_p, slflx_p, qtflx_p, uflx_p, vflx_p, slflx_cg_p, qtflx_cg_p, &
         uflx_cg_p, vflx_cg_p, pmam_cnst_idx_p, ptend_q_p, ptend_s_p, ptend_u_p, ptend_v_p, slten_p, qtten_p)

  end subroutine vertical_diffusion_core_batch_call

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_flux_diag_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (flux_diag_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('VERTICAL_DIFFUSION_FLUX_DIAG_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_flux_diag_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_flux_diag_impl = .false.
    end if

    flux_diag_impl_selected = .true.

    if (masterproc) then
       if (use_native_flux_diag_impl) then
          write(iulog,*) 'vertical_diffusion_flux_diag implementation = native'
       else
          write(iulog,*) 'vertical_diffusion_flux_diag implementation = codon'
       end if
    end if

  end subroutine vertical_diffusion_flux_diag_select_impl

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_flux_diag(ncol, q_tmp_local, s_tmp_local, u_tmp_local, v_tmp_local, &
       pint_local, zi_local, zm_local, cflx_local, kvh_local, kvm_local, cgs_local, cgh_local, &
       shflx_local, tautotx_local, tautoty_local, sl_local, qt_local, slv_local, slflx_local, qtflx_local, &
       uflx_local, vflx_local, slflx_cg_local, qtflx_cg_local, uflx_cg_local, vflx_cg_local)

    use iso_c_binding, only: c_int64_t, c_double, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: q_tmp_local(pcols,pver,pcnst)
    real(r8), target, intent(in) :: s_tmp_local(pcols,pver)
    real(r8), target, intent(in) :: u_tmp_local(pcols,pver)
    real(r8), target, intent(in) :: v_tmp_local(pcols,pver)
    real(r8), target, intent(in) :: pint_local(pcols,pverp)
    real(r8), target, intent(in) :: zi_local(pcols,pverp)
    real(r8), target, intent(in) :: zm_local(pcols,pver)
    real(r8), target, intent(in) :: cflx_local(pcols,pcnst)
    real(r8), target, intent(in) :: kvh_local(pcols,pverp)
    real(r8), target, intent(in) :: kvm_local(pcols,pverp)
    real(r8), target, intent(in) :: cgs_local(pcols,pverp)
    real(r8), target, intent(in) :: cgh_local(pcols,pverp)
    real(r8), target, intent(in) :: shflx_local(pcols)
    real(r8), target, intent(in) :: tautotx_local(pcols)
    real(r8), target, intent(in) :: tautoty_local(pcols)
    real(r8), target, intent(inout) :: sl_local(pcols,pver)
    real(r8), target, intent(inout) :: qt_local(pcols,pver)
    real(r8), target, intent(inout) :: slv_local(pcols,pver)
    real(r8), target, intent(inout) :: slflx_local(pcols,pverp)
    real(r8), target, intent(inout) :: qtflx_local(pcols,pverp)
    real(r8), target, intent(inout) :: uflx_local(pcols,pverp)
    real(r8), target, intent(inout) :: vflx_local(pcols,pverp)
    real(r8), target, intent(inout) :: slflx_cg_local(pcols,pverp)
    real(r8), target, intent(inout) :: qtflx_cg_local(pcols,pverp)
    real(r8), target, intent(inout) :: uflx_cg_local(pcols,pverp)
    real(r8), target, intent(inout) :: vflx_cg_local(pcols,pverp)

    interface
       subroutine vertical_diffusion_flux_diag_codon(ncol_c, pcols_c, pver_c, pverp_c, ixcldliq_c, ixcldice_c, &
            latvap_c, latice_c, zvir_c, rair_c, gravit_c, cpair_c, q_tmp_p, s_tmp_p, u_tmp_p, v_tmp_p, pint_p, &
            zi_p, zm_p, cflx_p, kvh_p, kvm_p, cgs_p, cgh_p, shflx_p, tautotx_p, tautoty_p, sl_p, qt_p, slv_p, &
            slflx_p, qtflx_p, uflx_p, vflx_p, slflx_cg_p, qtflx_cg_p, uflx_cg_p, vflx_cg_p) &
            bind(c, name="vertical_diffusion_flux_diag_codon")
         use iso_c_binding, only: c_int64_t, c_double, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pverp_c, ixcldliq_c, ixcldice_c
         real(c_double), value :: latvap_c, latice_c, zvir_c, rair_c, gravit_c, cpair_c
         type(c_ptr), value :: q_tmp_p, s_tmp_p, u_tmp_p, v_tmp_p, pint_p, zi_p, zm_p, cflx_p
         type(c_ptr), value :: kvh_p, kvm_p, cgs_p, cgh_p, shflx_p, tautotx_p, tautoty_p
         type(c_ptr), value :: sl_p, qt_p, slv_p, slflx_p, qtflx_p, uflx_p, vflx_p
         type(c_ptr), value :: slflx_cg_p, qtflx_cg_p, uflx_cg_p, vflx_cg_p
       end subroutine vertical_diffusion_flux_diag_codon
    end interface

    call vertical_diffusion_flux_diag_select_impl()

    if (use_native_flux_diag_impl) then
       call vertical_diffusion_flux_diag_native(ncol, q_tmp_local, s_tmp_local, u_tmp_local, v_tmp_local, &
            pint_local, zi_local, zm_local, cflx_local, kvh_local, kvm_local, cgs_local, cgh_local, shflx_local, &
            tautotx_local, tautoty_local, sl_local, qt_local, slv_local, slflx_local, qtflx_local, uflx_local, &
            vflx_local, slflx_cg_local, qtflx_cg_local, uflx_cg_local, vflx_cg_local)
       return
    end if

    call vertical_diffusion_flux_diag_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(pverp, c_int64_t), &
         int(ixcldliq, c_int64_t), int(ixcldice, c_int64_t), real(latvap, c_double), real(latice, c_double), &
         real(zvir, c_double), real(rair, c_double), real(gravit, c_double), real(cpair, c_double), &
         c_loc(q_tmp_local), c_loc(s_tmp_local), c_loc(u_tmp_local), c_loc(v_tmp_local), c_loc(pint_local), &
         c_loc(zi_local), c_loc(zm_local), c_loc(cflx_local), c_loc(kvh_local), c_loc(kvm_local), c_loc(cgs_local), &
         c_loc(cgh_local), c_loc(shflx_local), c_loc(tautotx_local), c_loc(tautoty_local), c_loc(sl_local), &
         c_loc(qt_local), c_loc(slv_local), c_loc(slflx_local), c_loc(qtflx_local), c_loc(uflx_local), &
         c_loc(vflx_local), c_loc(slflx_cg_local), c_loc(qtflx_cg_local), c_loc(uflx_cg_local), c_loc(vflx_cg_local) &
    )

  end subroutine vertical_diffusion_flux_diag

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_flux_diag_native(ncol, q_tmp_local, s_tmp_local, u_tmp_local, v_tmp_local, &
       pint_local, zi_local, zm_local, cflx_local, kvh_local, kvm_local, cgs_local, cgh_local, shflx_local, &
       tautotx_local, tautoty_local, sl_local, qt_local, slv_local, slflx_local, qtflx_local, uflx_local, &
       vflx_local, slflx_cg_local, qtflx_cg_local, uflx_cg_local, vflx_cg_local)

    integer, intent(in) :: ncol
    real(r8), intent(in) :: q_tmp_local(pcols,pver,pcnst)
    real(r8), intent(in) :: s_tmp_local(pcols,pver)
    real(r8), intent(in) :: u_tmp_local(pcols,pver)
    real(r8), intent(in) :: v_tmp_local(pcols,pver)
    real(r8), intent(in) :: pint_local(pcols,pverp)
    real(r8), intent(in) :: zi_local(pcols,pverp)
    real(r8), intent(in) :: zm_local(pcols,pver)
    real(r8), intent(in) :: cflx_local(pcols,pcnst)
    real(r8), intent(in) :: kvh_local(pcols,pverp)
    real(r8), intent(in) :: kvm_local(pcols,pverp)
    real(r8), intent(in) :: cgs_local(pcols,pverp)
    real(r8), intent(in) :: cgh_local(pcols,pverp)
    real(r8), intent(in) :: shflx_local(pcols)
    real(r8), intent(in) :: tautotx_local(pcols)
    real(r8), intent(in) :: tautoty_local(pcols)
    real(r8), intent(inout) :: sl_local(pcols,pver)
    real(r8), intent(inout) :: qt_local(pcols,pver)
    real(r8), intent(inout) :: slv_local(pcols,pver)
    real(r8), intent(inout) :: slflx_local(pcols,pverp)
    real(r8), intent(inout) :: qtflx_local(pcols,pverp)
    real(r8), intent(inout) :: uflx_local(pcols,pverp)
    real(r8), intent(inout) :: vflx_local(pcols,pverp)
    real(r8), intent(inout) :: slflx_cg_local(pcols,pverp)
    real(r8), intent(inout) :: qtflx_cg_local(pcols,pverp)
    real(r8), intent(inout) :: uflx_cg_local(pcols,pverp)
    real(r8), intent(inout) :: vflx_cg_local(pcols,pverp)

    integer :: i, k
    real(r8) :: rhoair_local

    sl_local(:ncol,:pver)  = s_tmp_local(:ncol,:) -   latvap           * q_tmp_local(:ncol,:,ixcldliq) &
                                               - ( latvap + latice) * q_tmp_local(:ncol,:,ixcldice)
    qt_local(:ncol,:pver)  = q_tmp_local(:ncol,:,1) + q_tmp_local(:ncol,:,ixcldliq) &
                                              + q_tmp_local(:ncol,:,ixcldice)
    slv_local(:ncol,:pver) = sl_local(:ncol,:pver) * ( 1._r8 + zvir*qt_local(:ncol,:pver) )

    slflx_local(:ncol,1) = 0._r8
    qtflx_local(:ncol,1) = 0._r8
    uflx_local(:ncol,1)  = 0._r8
    vflx_local(:ncol,1)  = 0._r8

    slflx_cg_local(:ncol,1) = 0._r8
    qtflx_cg_local(:ncol,1) = 0._r8
    uflx_cg_local(:ncol,1)  = 0._r8
    vflx_cg_local(:ncol,1)  = 0._r8

    do k = 2, pver
       do i = 1, ncol
          rhoair_local     = pint_local(i,k) / &
               ( rair * ( ( 0.5_r8*(slv_local(i,k)+slv_local(i,k-1)) - gravit*zi_local(i,k))/cpair ) )
          slflx_local(i,k) = kvh_local(i,k) * &
                               ( - rhoair_local*(sl_local(i,k-1)-sl_local(i,k))/(zm_local(i,k-1)-zm_local(i,k)) &
                                 + cgh_local(i,k) )
          qtflx_local(i,k) = kvh_local(i,k) * &
                               ( - rhoair_local*(qt_local(i,k-1)-qt_local(i,k))/(zm_local(i,k-1)-zm_local(i,k)) &
                                 + rhoair_local*(cflx_local(i,1)+cflx_local(i,ixcldliq)+cflx_local(i,ixcldice))*cgs_local(i,k) )
          uflx_local(i,k)  = kvm_local(i,k) * &
                               ( - rhoair_local*(u_tmp_local(i,k-1)-u_tmp_local(i,k))/(zm_local(i,k-1)-zm_local(i,k)))
          vflx_local(i,k)  = kvm_local(i,k) * &
                               ( - rhoair_local*(v_tmp_local(i,k-1)-v_tmp_local(i,k))/(zm_local(i,k-1)-zm_local(i,k)))
          slflx_cg_local(i,k) = kvh_local(i,k) * cgh_local(i,k)
          qtflx_cg_local(i,k) = kvh_local(i,k) * rhoair_local * &
               ( cflx_local(i,1) + cflx_local(i,ixcldliq) + cflx_local(i,ixcldice) ) * cgs_local(i,k)
          uflx_cg_local(i,k)  = 0._r8
          vflx_cg_local(i,k)  = 0._r8
       end do
    end do

    slflx_local(:ncol,pverp) = shflx_local(:ncol)
    qtflx_local(:ncol,pverp) = cflx_local(:ncol,1)
    uflx_local(:ncol,pverp)  = tautotx_local(:ncol)
    vflx_local(:ncol,pverp)  = tautoty_local(:ncol)

    slflx_cg_local(:ncol,pverp) = 0._r8
    qtflx_cg_local(:ncol,pverp) = 0._r8
    uflx_cg_local(:ncol,pverp)  = 0._r8
    vflx_cg_local(:ncol,pverp)  = 0._r8

  end subroutine vertical_diffusion_flux_diag_native

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_ptend_core_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (ptend_core_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('VERTICAL_DIFFUSION_PTEND_CORE_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_ptend_core_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_ptend_core_impl = .false.
    end if

    ptend_core_impl_selected = .true.

    if (masterproc) then
       if (use_native_ptend_core_impl) then
          write(iulog,*) 'vertical_diffusion_ptend_core implementation = native'
       else
          write(iulog,*) 'vertical_diffusion_ptend_core implementation = codon'
       end if
    end if

  end subroutine vertical_diffusion_ptend_core_select_impl

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_ptend_core(ncol, psetcols_local, q_tmp_local, s_tmp_local, u_tmp_local, v_tmp_local, &
       state_q_local, state_s_local, state_u_local, state_v_local, sl_local, qt_local, sl_prePBL_local, &
       qt_prePBL_local, rztodt_local, ptend_q_local, ptend_s_local, ptend_u_local, ptend_v_local, slten_local, &
       qtten_local)

    use iso_c_binding, only: c_int64_t, c_double, c_loc, c_ptr

    integer, intent(in) :: ncol
    integer, intent(in) :: psetcols_local
    real(r8), target, intent(in) :: q_tmp_local(pcols,pver,pcnst)
    real(r8), target, intent(in) :: s_tmp_local(pcols,pver)
    real(r8), target, intent(in) :: u_tmp_local(pcols,pver)
    real(r8), target, intent(in) :: v_tmp_local(pcols,pver)
    real(r8), target, intent(in) :: state_q_local(pcols,pver,pcnst)
    real(r8), target, intent(in) :: state_s_local(pcols,pver)
    real(r8), target, intent(in) :: state_u_local(pcols,pver)
    real(r8), target, intent(in) :: state_v_local(pcols,pver)
    real(r8), target, intent(in) :: sl_local(pcols,pver)
    real(r8), target, intent(in) :: qt_local(pcols,pver)
    real(r8), target, intent(in) :: sl_prePBL_local(pcols,pver)
    real(r8), target, intent(in) :: qt_prePBL_local(pcols,pver)
    real(r8), intent(in) :: rztodt_local
    real(r8), target, intent(inout) :: ptend_q_local(psetcols_local,pver,pcnst)
    real(r8), target, intent(inout) :: ptend_s_local(psetcols_local,pver)
    real(r8), target, intent(inout) :: ptend_u_local(psetcols_local,pver)
    real(r8), target, intent(inout) :: ptend_v_local(psetcols_local,pver)
    real(r8), target, intent(inout) :: slten_local(pcols,pver)
    real(r8), target, intent(inout) :: qtten_local(pcols,pver)

    interface
       subroutine vertical_diffusion_ptend_core_codon(ncol_c, pcols_c, pver_c, pcnst_c, psetcols_c, rztodt_c, &
            q_tmp_p, s_tmp_p, u_tmp_p, v_tmp_p, state_q_p, state_s_p, state_u_p, state_v_p, sl_p, qt_p, &
            sl_prePBL_p, qt_prePBL_p, ptend_q_p, ptend_s_p, ptend_u_p, ptend_v_p, slten_p, qtten_p) &
            bind(c, name="vertical_diffusion_ptend_core_codon")
         use iso_c_binding, only: c_int64_t, c_double, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pcnst_c, psetcols_c
         real(c_double), value :: rztodt_c
         type(c_ptr), value :: q_tmp_p, s_tmp_p, u_tmp_p, v_tmp_p, state_q_p, state_s_p, state_u_p, state_v_p
         type(c_ptr), value :: sl_p, qt_p, sl_prePBL_p, qt_prePBL_p, ptend_q_p, ptend_s_p, ptend_u_p, ptend_v_p
         type(c_ptr), value :: slten_p, qtten_p
       end subroutine vertical_diffusion_ptend_core_codon
    end interface

    call vertical_diffusion_ptend_core_select_impl()

    if (use_native_ptend_core_impl) then
       call vertical_diffusion_ptend_core_native(ncol, psetcols_local, q_tmp_local, s_tmp_local, u_tmp_local, &
            v_tmp_local, state_q_local, state_s_local, state_u_local, state_v_local, sl_local, qt_local, &
            sl_prePBL_local, qt_prePBL_local, rztodt_local, ptend_q_local, ptend_s_local, ptend_u_local, &
            ptend_v_local, slten_local, qtten_local)
       return
    end if

    call vertical_diffusion_ptend_core_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(pcnst, c_int64_t), &
         int(psetcols_local, c_int64_t), real(rztodt_local, c_double), c_loc(q_tmp_local), c_loc(s_tmp_local), &
         c_loc(u_tmp_local), c_loc(v_tmp_local), c_loc(state_q_local), c_loc(state_s_local), c_loc(state_u_local), &
         c_loc(state_v_local), c_loc(sl_local), c_loc(qt_local), c_loc(sl_prePBL_local), c_loc(qt_prePBL_local), &
         c_loc(ptend_q_local), c_loc(ptend_s_local), c_loc(ptend_u_local), c_loc(ptend_v_local), c_loc(slten_local), &
         c_loc(qtten_local) &
    )

  end subroutine vertical_diffusion_ptend_core

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_ptend_core_native(ncol, psetcols_local, q_tmp_local, s_tmp_local, u_tmp_local, &
       v_tmp_local, state_q_local, state_s_local, state_u_local, state_v_local, sl_local, qt_local, &
       sl_prePBL_local, qt_prePBL_local, rztodt_local, ptend_q_local, ptend_s_local, ptend_u_local, ptend_v_local, &
       slten_local, qtten_local)

    integer, intent(in) :: ncol
    integer, intent(in) :: psetcols_local
    real(r8), intent(in) :: q_tmp_local(pcols,pver,pcnst)
    real(r8), intent(in) :: s_tmp_local(pcols,pver)
    real(r8), intent(in) :: u_tmp_local(pcols,pver)
    real(r8), intent(in) :: v_tmp_local(pcols,pver)
    real(r8), intent(in) :: state_q_local(pcols,pver,pcnst)
    real(r8), intent(in) :: state_s_local(pcols,pver)
    real(r8), intent(in) :: state_u_local(pcols,pver)
    real(r8), intent(in) :: state_v_local(pcols,pver)
    real(r8), intent(in) :: sl_local(pcols,pver)
    real(r8), intent(in) :: qt_local(pcols,pver)
    real(r8), intent(in) :: sl_prePBL_local(pcols,pver)
    real(r8), intent(in) :: qt_prePBL_local(pcols,pver)
    real(r8), intent(in) :: rztodt_local
    real(r8), intent(inout) :: ptend_q_local(psetcols_local,pver,pcnst)
    real(r8), intent(inout) :: ptend_s_local(psetcols_local,pver)
    real(r8), intent(inout) :: ptend_u_local(psetcols_local,pver)
    real(r8), intent(inout) :: ptend_v_local(psetcols_local,pver)
    real(r8), intent(inout) :: slten_local(pcols,pver)
    real(r8), intent(inout) :: qtten_local(pcols,pver)

    ptend_s_local(:ncol,:)       = ( s_tmp_local(:ncol,:) - state_s_local(:ncol,:) ) * rztodt_local
    ptend_u_local(:ncol,:)       = ( u_tmp_local(:ncol,:) - state_u_local(:ncol,:) ) * rztodt_local
    ptend_v_local(:ncol,:)       = ( v_tmp_local(:ncol,:) - state_v_local(:ncol,:) ) * rztodt_local
    ptend_q_local(:ncol,:pver,:) = ( q_tmp_local(:ncol,:pver,:) - state_q_local(:ncol,:pver,:) ) * rztodt_local
    slten_local(:ncol,:)         = ( sl_local(:ncol,:) - sl_prePBL_local(:ncol,:) ) * rztodt_local
    qtten_local(:ncol,:)         = ( qt_local(:ncol,:) - qt_prePBL_local(:ncol,:) ) * rztodt_local

  end subroutine vertical_diffusion_ptend_core_native

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_pre_pbl_diag_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (pre_pbl_diag_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('VERTICAL_DIFFUSION_PRE_PBL_DIAG_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_pre_pbl_diag_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_pre_pbl_diag_impl = .false.
    end if

    pre_pbl_diag_impl_selected = .true.

    if (masterproc) then
       if (use_native_pre_pbl_diag_impl) then
          write(iulog,*) 'vertical_diffusion_pre_pbl_diag implementation = native'
       else
          write(iulog,*) 'vertical_diffusion_pre_pbl_diag implementation = codon'
       end if
    end if

  end subroutine vertical_diffusion_pre_pbl_diag_select_impl

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_pre_pbl_diag(ncol, state_q_local, state_s_local, state_u_local, state_v_local, &
       q_tmp_local, s_tmp_local, u_tmp_local, v_tmp_local, sl_prePBL_local, qt_prePBL_local, slv_prePBL_local)

    use iso_c_binding, only: c_int64_t, c_double, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: state_q_local(pcols,pver,pcnst)
    real(r8), target, intent(in) :: state_s_local(pcols,pver)
    real(r8), target, intent(in) :: state_u_local(pcols,pver)
    real(r8), target, intent(in) :: state_v_local(pcols,pver)
    real(r8), target, intent(inout) :: q_tmp_local(pcols,pver,pcnst)
    real(r8), target, intent(inout) :: s_tmp_local(pcols,pver)
    real(r8), target, intent(inout) :: u_tmp_local(pcols,pver)
    real(r8), target, intent(inout) :: v_tmp_local(pcols,pver)
    real(r8), target, intent(inout) :: sl_prePBL_local(pcols,pver)
    real(r8), target, intent(inout) :: qt_prePBL_local(pcols,pver)
    real(r8), target, intent(inout) :: slv_prePBL_local(pcols,pver)

    interface
       subroutine vertical_diffusion_pre_pbl_diag_codon(ncol_c, pcols_c, pver_c, pcnst_c, ixcldliq_c, ixcldice_c, &
            latvap_c, latice_c, zvir_c, state_q_p, state_s_p, state_u_p, state_v_p, q_tmp_p, s_tmp_p, u_tmp_p, &
            v_tmp_p, sl_prePBL_p, qt_prePBL_p, slv_prePBL_p) bind(c, name="vertical_diffusion_pre_pbl_diag_codon")
         use iso_c_binding, only: c_int64_t, c_double, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pcnst_c, ixcldliq_c, ixcldice_c
         real(c_double), value :: latvap_c, latice_c, zvir_c
         type(c_ptr), value :: state_q_p, state_s_p, state_u_p, state_v_p, q_tmp_p, s_tmp_p, u_tmp_p, v_tmp_p
         type(c_ptr), value :: sl_prePBL_p, qt_prePBL_p, slv_prePBL_p
       end subroutine vertical_diffusion_pre_pbl_diag_codon
    end interface

    call vertical_diffusion_diag_batch_select_impl()

    if (use_native_diag_batch_impl) then
       call vertical_diffusion_pre_pbl_diag_native(ncol, state_q_local, state_s_local, state_u_local, state_v_local, &
            q_tmp_local, s_tmp_local, u_tmp_local, v_tmp_local, sl_prePBL_local, qt_prePBL_local, slv_prePBL_local)
       return
    end if

    call vertical_diffusion_diag_batch_call(1, ncol, state_q_local=state_q_local, state_s_local=state_s_local, &
         state_u_local=state_u_local, state_v_local=state_v_local, q_tmp_local=q_tmp_local, s_tmp_local=s_tmp_local, &
         u_tmp_local=u_tmp_local, v_tmp_local=v_tmp_local, sl_prePBL_local=sl_prePBL_local, &
         qt_prePBL_local=qt_prePBL_local, slv_prePBL_local=slv_prePBL_local)

  end subroutine vertical_diffusion_pre_pbl_diag

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_pre_pbl_diag_native(ncol, state_q_local, state_s_local, state_u_local, state_v_local, &
       q_tmp_local, s_tmp_local, u_tmp_local, v_tmp_local, sl_prePBL_local, qt_prePBL_local, slv_prePBL_local)

    integer, intent(in) :: ncol
    real(r8), intent(in) :: state_q_local(pcols,pver,pcnst)
    real(r8), intent(in) :: state_s_local(pcols,pver)
    real(r8), intent(in) :: state_u_local(pcols,pver)
    real(r8), intent(in) :: state_v_local(pcols,pver)
    real(r8), intent(inout) :: q_tmp_local(pcols,pver,pcnst)
    real(r8), intent(inout) :: s_tmp_local(pcols,pver)
    real(r8), intent(inout) :: u_tmp_local(pcols,pver)
    real(r8), intent(inout) :: v_tmp_local(pcols,pver)
    real(r8), intent(inout) :: sl_prePBL_local(pcols,pver)
    real(r8), intent(inout) :: qt_prePBL_local(pcols,pver)
    real(r8), intent(inout) :: slv_prePBL_local(pcols,pver)

    q_tmp_local(:ncol,:,:) = state_q_local(:ncol,:,:)
    s_tmp_local(:ncol,:) = state_s_local(:ncol,:)
    u_tmp_local(:ncol,:) = state_u_local(:ncol,:)
    v_tmp_local(:ncol,:) = state_v_local(:ncol,:)

    sl_prePBL_local(:ncol,:pver)  = s_tmp_local(:ncol,:) -   latvap * q_tmp_local(:ncol,:,ixcldliq) &
                                                   - ( latvap + latice) * q_tmp_local(:ncol,:,ixcldice)
    qt_prePBL_local(:ncol,:pver)  = q_tmp_local(:ncol,:,1) + q_tmp_local(:ncol,:,ixcldliq) &
                                                      + q_tmp_local(:ncol,:,ixcldice)
    slv_prePBL_local(:ncol,:pver) = sl_prePBL_local(:ncol,:pver) * ( 1._r8 + zvir*qt_prePBL_local(:ncol,:pver) )

  end subroutine vertical_diffusion_pre_pbl_diag_native

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_post_pbl_state_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (post_pbl_state_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('VERTICAL_DIFFUSION_POST_PBL_STATE_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_post_pbl_state_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_post_pbl_state_impl = .false.
    end if

    post_pbl_state_impl_selected = .true.

    if (masterproc) then
       if (use_native_post_pbl_state_impl) then
          write(iulog,*) 'vertical_diffusion_post_pbl_state implementation = native'
       else
          write(iulog,*) 'vertical_diffusion_post_pbl_state implementation = codon'
       end if
    end if

  end subroutine vertical_diffusion_post_pbl_state_select_impl

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_post_pbl_state(ncol, psetcols_local, state_q_local, state_s_local, state_u_local, &
       state_v_local, state_zm_local, ptend_q_local, ptend_s_local, ptend_u_local, ptend_v_local, ztodt_local, &
       qv_aft_PBL_local, ql_aft_PBL_local, qi_aft_PBL_local, s_aft_PBL_local, t_aftPBL_local, u_aft_PBL_local, &
       v_aft_PBL_local)

    use iso_c_binding, only: c_int64_t, c_double, c_loc, c_ptr

    integer, intent(in) :: ncol
    integer, intent(in) :: psetcols_local
    real(r8), target, intent(in) :: state_q_local(pcols,pver,pcnst)
    real(r8), target, intent(in) :: state_s_local(pcols,pver)
    real(r8), target, intent(in) :: state_u_local(pcols,pver)
    real(r8), target, intent(in) :: state_v_local(pcols,pver)
    real(r8), target, intent(in) :: state_zm_local(pcols,pver)
    real(r8), target, intent(in) :: ptend_q_local(psetcols_local,pver,pcnst)
    real(r8), target, intent(in) :: ptend_s_local(psetcols_local,pver)
    real(r8), target, intent(in) :: ptend_u_local(psetcols_local,pver)
    real(r8), target, intent(in) :: ptend_v_local(psetcols_local,pver)
    real(r8), intent(in) :: ztodt_local
    real(r8), target, intent(inout) :: qv_aft_PBL_local(pcols,pver)
    real(r8), target, intent(inout) :: ql_aft_PBL_local(pcols,pver)
    real(r8), target, intent(inout) :: qi_aft_PBL_local(pcols,pver)
    real(r8), target, intent(inout) :: s_aft_PBL_local(pcols,pver)
    real(r8), target, intent(inout) :: t_aftPBL_local(pcols,pver)
    real(r8), target, intent(inout) :: u_aft_PBL_local(pcols,pver)
    real(r8), target, intent(inout) :: v_aft_PBL_local(pcols,pver)

    interface
       subroutine vertical_diffusion_post_pbl_state_codon(ncol_c, pcols_c, pver_c, pcnst_c, psetcols_c, &
            ixcldliq_c, ixcldice_c, ztodt_c, gravit_c, cpair_c, state_q_p, state_s_p, state_u_p, state_v_p, &
            state_zm_p, ptend_q_p, ptend_s_p, ptend_u_p, ptend_v_p, qv_aft_PBL_p, ql_aft_PBL_p, qi_aft_PBL_p, &
            s_aft_PBL_p, t_aftPBL_p, u_aft_PBL_p, v_aft_PBL_p) bind(c, name="vertical_diffusion_post_pbl_state_codon")
         use iso_c_binding, only: c_int64_t, c_double, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pcnst_c, psetcols_c, ixcldliq_c, ixcldice_c
         real(c_double), value :: ztodt_c, gravit_c, cpair_c
         type(c_ptr), value :: state_q_p, state_s_p, state_u_p, state_v_p, state_zm_p
         type(c_ptr), value :: ptend_q_p, ptend_s_p, ptend_u_p, ptend_v_p
         type(c_ptr), value :: qv_aft_PBL_p, ql_aft_PBL_p, qi_aft_PBL_p, s_aft_PBL_p, t_aftPBL_p
         type(c_ptr), value :: u_aft_PBL_p, v_aft_PBL_p
       end subroutine vertical_diffusion_post_pbl_state_codon
    end interface

    call vertical_diffusion_diag_batch_select_impl()

    if (use_native_diag_batch_impl) then
      call vertical_diffusion_post_pbl_state_native(ncol, psetcols_local, state_q_local, state_s_local, state_u_local, &
           state_v_local, state_zm_local, ptend_q_local, ptend_s_local, ptend_u_local, ptend_v_local, ztodt_local, &
           qv_aft_PBL_local, ql_aft_PBL_local, qi_aft_PBL_local, s_aft_PBL_local, t_aftPBL_local, u_aft_PBL_local, &
           v_aft_PBL_local)
      return
    end if

    call vertical_diffusion_diag_batch_call(3, ncol, psetcols_local=psetcols_local, ztodt_local=ztodt_local, &
         state_q_local=state_q_local, state_s_local=state_s_local, state_u_local=state_u_local, &
         state_v_local=state_v_local, state_zm_local=state_zm_local, ptend_q_local=ptend_q_local, &
         ptend_s_local=ptend_s_local, ptend_u_local=ptend_u_local, ptend_v_local=ptend_v_local, &
         qv_aft_PBL_local=qv_aft_PBL_local, ql_aft_PBL_local=ql_aft_PBL_local, qi_aft_PBL_local=qi_aft_PBL_local, &
         s_aft_PBL_local=s_aft_PBL_local, t_aftPBL_local=t_aftPBL_local, u_aft_PBL_local=u_aft_PBL_local, &
         v_aft_PBL_local=v_aft_PBL_local)

  end subroutine vertical_diffusion_post_pbl_state

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_post_pbl_state_native(ncol, psetcols_local, state_q_local, state_s_local, &
       state_u_local, state_v_local, state_zm_local, ptend_q_local, ptend_s_local, ptend_u_local, ptend_v_local, &
       ztodt_local, qv_aft_PBL_local, ql_aft_PBL_local, qi_aft_PBL_local, s_aft_PBL_local, t_aftPBL_local, &
       u_aft_PBL_local, v_aft_PBL_local)

    integer, intent(in) :: ncol
    integer, intent(in) :: psetcols_local
    real(r8), intent(in) :: state_q_local(pcols,pver,pcnst)
    real(r8), intent(in) :: state_s_local(pcols,pver)
    real(r8), intent(in) :: state_u_local(pcols,pver)
    real(r8), intent(in) :: state_v_local(pcols,pver)
    real(r8), intent(in) :: state_zm_local(pcols,pver)
    real(r8), intent(in) :: ptend_q_local(psetcols_local,pver,pcnst)
    real(r8), intent(in) :: ptend_s_local(psetcols_local,pver)
    real(r8), intent(in) :: ptend_u_local(psetcols_local,pver)
    real(r8), intent(in) :: ptend_v_local(psetcols_local,pver)
    real(r8), intent(in) :: ztodt_local
    real(r8), intent(inout) :: qv_aft_PBL_local(pcols,pver)
    real(r8), intent(inout) :: ql_aft_PBL_local(pcols,pver)
    real(r8), intent(inout) :: qi_aft_PBL_local(pcols,pver)
    real(r8), intent(inout) :: s_aft_PBL_local(pcols,pver)
    real(r8), intent(inout) :: t_aftPBL_local(pcols,pver)
    real(r8), intent(inout) :: u_aft_PBL_local(pcols,pver)
    real(r8), intent(inout) :: v_aft_PBL_local(pcols,pver)

    qv_aft_PBL_local(:ncol,:pver) = state_q_local(:ncol,:pver,1)         + ptend_q_local(:ncol,:pver,1)        * ztodt_local
    ql_aft_PBL_local(:ncol,:pver) = state_q_local(:ncol,:pver,ixcldliq)  + ptend_q_local(:ncol,:pver,ixcldliq) * ztodt_local
    qi_aft_PBL_local(:ncol,:pver) = state_q_local(:ncol,:pver,ixcldice)  + ptend_q_local(:ncol,:pver,ixcldice) * ztodt_local
    s_aft_PBL_local(:ncol,:pver)  = state_s_local(:ncol,:pver)           + ptend_s_local(:ncol,:pver)          * ztodt_local
    t_aftPBL_local(:ncol,:pver)   = ( s_aft_PBL_local(:ncol,:pver) - gravit*state_zm_local(:ncol,:pver) ) / cpair
    u_aft_PBL_local(:ncol,:pver)  = state_u_local(:ncol,:pver)           + ptend_u_local(:ncol,:pver)          * ztodt_local
    v_aft_PBL_local(:ncol,:pver)  = state_v_local(:ncol,:pver)           + ptend_v_local(:ncol,:pver)          * ztodt_local

  end subroutine vertical_diffusion_post_pbl_state_native

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_modal_aero_flux_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (modal_aero_flux_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('VERTICAL_DIFFUSION_MODAL_AERO_FLUX_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_modal_aero_flux_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_modal_aero_flux_impl = .false.
    end if

    modal_aero_flux_impl_selected = .true.

    if (masterproc) then
       if (use_native_modal_aero_flux_impl) then
          write(iulog,*) 'vertical_diffusion_modal_aero_flux implementation = native'
       else
          write(iulog,*) 'vertical_diffusion_modal_aero_flux implementation = codon'
       end if
    end if

  end subroutine vertical_diffusion_modal_aero_flux_select_impl

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_modal_aero_flux(ncol, state_rpdel_local, cflx_local, ztodt_local, q_tmp_local)

    use iso_c_binding, only: c_int64_t, c_double, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: state_rpdel_local(pcols,pver)
    real(r8), target, intent(in) :: cflx_local(pcols,pcnst)
    real(r8), intent(in) :: ztodt_local
    real(r8), target, intent(inout) :: q_tmp_local(pcols,pver,pcnst)

    interface
       subroutine vertical_diffusion_modal_aero_flux_codon(ncol_c, pcols_c, pver_c, pcnst_c, pmam_ncnst_c, ztodt_c, &
            gravit_c, state_rpdel_p, pmam_cnst_idx_p, cflx_p, q_tmp_p) bind(c, name="vertical_diffusion_modal_aero_flux_codon")
         use iso_c_binding, only: c_int64_t, c_double, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pcnst_c, pmam_ncnst_c
         real(c_double), value :: ztodt_c, gravit_c
         type(c_ptr), value :: state_rpdel_p, pmam_cnst_idx_p, cflx_p, q_tmp_p
       end subroutine vertical_diffusion_modal_aero_flux_codon
    end interface

    call vertical_diffusion_modal_aero_flux_select_impl()

    if (use_native_modal_aero_flux_impl) then
      call vertical_diffusion_modal_aero_flux_native(ncol, state_rpdel_local, cflx_local, ztodt_local, q_tmp_local)
      return
    end if

    call vertical_diffusion_modal_aero_flux_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(pcnst, c_int64_t), &
         int(pmam_ncnst, c_int64_t), real(ztodt_local, c_double), real(gravit, c_double), c_loc(state_rpdel_local), &
         c_loc(pmam_cnst_idx_c), c_loc(cflx_local), c_loc(q_tmp_local) &
    )

  end subroutine vertical_diffusion_modal_aero_flux

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_modal_aero_flux_native(ncol, state_rpdel_local, cflx_local, ztodt_local, q_tmp_local)

    integer, intent(in) :: ncol
    real(r8), intent(in) :: state_rpdel_local(pcols,pver)
    real(r8), intent(in) :: cflx_local(pcols,pcnst)
    real(r8), intent(in) :: ztodt_local
    real(r8), intent(inout) :: q_tmp_local(pcols,pver,pcnst)

    integer :: l, m
    real(r8) :: tmp1_local(pcols)

    tmp1_local(:ncol) = ztodt_local * gravit * state_rpdel_local(:ncol,pver)
    do m = 1, pmam_ncnst
       l = pmam_cnst_idx(m)
       q_tmp_local(:ncol,pver,l) = q_tmp_local(:ncol,pver,l) + tmp1_local(:ncol) * cflx_local(:ncol,l)
    enddo

  end subroutine vertical_diffusion_modal_aero_flux_native

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_obklen_diag_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (obklen_diag_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('VERTICAL_DIFFUSION_OBKLEN_DIAG_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_obklen_diag_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_obklen_diag_impl = .false.
    end if

    obklen_diag_impl_selected = .true.

    if (masterproc) then
       if (use_native_obklen_diag_impl) then
          write(iulog,*) 'vertical_diffusion_obklen_diag implementation = native'
       else
          write(iulog,*) 'vertical_diffusion_obklen_diag implementation = codon'
       end if
    end if

  end subroutine vertical_diffusion_obklen_diag_select_impl

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_obklen_diag(ncol, state_t_local, state_exner_local, state_q_local, cflx_local, &
       shflx_local, rrho_local, ustar_local, th_local, thvs_local, khfs_local, kqfs_local, kbfs_local, obklen_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: state_t_local(pcols,pver), state_exner_local(pcols,pver)
    real(r8), target, intent(in) :: state_q_local(pcols,pver,pcnst), cflx_local(pcols,pcnst), shflx_local(pcols)
    real(r8), target, intent(in) :: rrho_local(pcols), ustar_local(pcols)
    real(r8), target, intent(inout) :: th_local(pcols,pver), thvs_local(pcols)
    real(r8), target, intent(inout) :: khfs_local(pcols), kqfs_local(pcols), kbfs_local(pcols), obklen_local(pcols)

    interface
       subroutine vertical_diffusion_obklen_diag_codon(ncol_c, pcols_c, pver_c, cpair_c, zvir_c, gravit_c, karman_c, &
            state_t_p, state_exner_p, state_q_p, cflx_p, shflx_p, rrho_p, ustar_p, th_p, thvs_p, khfs_p, kqfs_p, &
            kbfs_p, obklen_p) bind(c, name="vertical_diffusion_obklen_diag_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         real(c_double), value :: cpair_c, zvir_c, gravit_c, karman_c
         type(c_ptr), value :: state_t_p, state_exner_p, state_q_p, cflx_p, shflx_p, rrho_p, ustar_p
         type(c_ptr), value :: th_p, thvs_p, khfs_p, kqfs_p, kbfs_p, obklen_p
       end subroutine vertical_diffusion_obklen_diag_codon
    end interface

    call vertical_diffusion_obklen_diag_select_impl()

    if (use_native_obklen_diag_impl) then
       call vertical_diffusion_obklen_diag_native(ncol, state_t_local, state_exner_local, state_q_local, cflx_local, &
            shflx_local, rrho_local, ustar_local, th_local, thvs_local, khfs_local, kqfs_local, kbfs_local, &
            obklen_local)
       return
    end if

    call vertical_diffusion_obklen_diag_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), real(cpair, c_double), &
         real(zvir, c_double), real(gravit, c_double), real(karman, c_double), c_loc(state_t_local), &
         c_loc(state_exner_local), c_loc(state_q_local), c_loc(cflx_local), c_loc(shflx_local), c_loc(rrho_local), &
         c_loc(ustar_local), c_loc(th_local), c_loc(thvs_local), c_loc(khfs_local), c_loc(kqfs_local), c_loc(kbfs_local), &
         c_loc(obklen_local) &
    )

  end subroutine vertical_diffusion_obklen_diag

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_obklen_diag_native(ncol, state_t_local, state_exner_local, state_q_local, cflx_local, &
       shflx_local, rrho_local, ustar_local, th_local, thvs_local, khfs_local, kqfs_local, kbfs_local, obklen_local)

    use pbl_utils, only: calc_obklen, virtem

    integer, intent(in) :: ncol
    real(r8), intent(in) :: state_t_local(pcols,pver), state_exner_local(pcols,pver)
    real(r8), intent(in) :: state_q_local(pcols,pver,pcnst), cflx_local(pcols,pcnst), shflx_local(pcols)
    real(r8), intent(in) :: rrho_local(pcols), ustar_local(pcols)
    real(r8), intent(inout) :: th_local(pcols,pver), thvs_local(pcols)
    real(r8), intent(inout) :: khfs_local(pcols), kqfs_local(pcols), kbfs_local(pcols), obklen_local(pcols)

    th_local(:ncol,pver) = state_t_local(:ncol,pver) * state_exner_local(:ncol,pver)
    thvs_local(:ncol) = virtem(th_local(:ncol,pver), state_q_local(:ncol,pver,1))
    call calc_obklen(th_local(:ncol,pver), thvs_local(:ncol), cflx_local(:ncol,1), shflx_local(:ncol), rrho_local(:ncol), &
         ustar_local(:ncol), khfs_local(:ncol), kqfs_local(:ncol), kbfs_local(:ncol), obklen_local(:ncol))

  end subroutine vertical_diffusion_obklen_diag_native

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_pre_qsat_rh_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (pre_qsat_rh_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('VERTICAL_DIFFUSION_PRE_QSAT_RH_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_pre_qsat_rh_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_pre_qsat_rh_impl = .false.
    end if

    pre_qsat_rh_impl_selected = .true.

    if (masterproc) then
       if (use_native_pre_qsat_rh_impl) then
          write(iulog,*) 'vertical_diffusion_pre_qsat_rh implementation = native'
       else
          write(iulog,*) 'vertical_diffusion_pre_qsat_rh implementation = codon'
       end if
    end if

  end subroutine vertical_diffusion_pre_qsat_rh_select_impl

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_pre_qsat_rh(ncol, state_q_local, ftem_local, ftem_prePBL_local)

    use iso_c_binding, only: c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: state_q_local(pcols,pver,pcnst)
    real(r8), target, intent(in) :: ftem_local(pcols,pver)
    real(r8), target, intent(inout) :: ftem_prePBL_local(pcols,pver)

    interface
       subroutine vertical_diffusion_pre_qsat_rh_codon(ncol_c, pcols_c, pver_c, state_q_p, ftem_p, ftem_prePBL_p) &
            bind(c, name="vertical_diffusion_pre_qsat_rh_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         type(c_ptr), value :: state_q_p, ftem_p, ftem_prePBL_p
       end subroutine vertical_diffusion_pre_qsat_rh_codon
    end interface

    call vertical_diffusion_diag_batch_select_impl()

    if (use_native_diag_batch_impl) then
       call vertical_diffusion_pre_qsat_rh_native(ncol, state_q_local, ftem_local, ftem_prePBL_local)
       return
    end if

    call vertical_diffusion_diag_batch_call(2, ncol, state_q_local=state_q_local, ftem_local=ftem_local, &
         ftem_prePBL_local=ftem_prePBL_local)

  end subroutine vertical_diffusion_pre_qsat_rh

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_pre_qsat_rh_native(ncol, state_q_local, ftem_local, ftem_prePBL_local)

    integer, intent(in) :: ncol
    real(r8), intent(in) :: state_q_local(pcols,pver,pcnst)
    real(r8), intent(in) :: ftem_local(pcols,pver)
    real(r8), intent(inout) :: ftem_prePBL_local(pcols,pver)

    ftem_prePBL_local(:ncol,:pver) = state_q_local(:ncol,:pver,1) / ftem_local(:ncol,:pver) * 100._r8

  end subroutine vertical_diffusion_pre_qsat_rh_native

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_post_qsat_diag_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (post_qsat_diag_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('VERTICAL_DIFFUSION_POST_QSAT_DIAG_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_post_qsat_diag_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_post_qsat_diag_impl = .false.
    end if

    post_qsat_diag_impl_selected = .true.

    if (masterproc) then
       if (use_native_post_qsat_diag_impl) then
          write(iulog,*) 'vertical_diffusion_post_qsat_diag implementation = native'
       else
          write(iulog,*) 'vertical_diffusion_post_qsat_diag implementation = codon'
       end if
    end if

  end subroutine vertical_diffusion_post_qsat_diag_select_impl

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_post_qsat_diag(ncol, state_t_local, qv_aft_PBL_local, ftem_prePBL_local, &
       t_aftPBL_local, ftem_local, rztodt_local, ftem_aftPBL_local, tten_local, rhten_local)

    use iso_c_binding, only: c_int64_t, c_double, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: state_t_local(pcols,pver)
    real(r8), target, intent(in) :: qv_aft_PBL_local(pcols,pver)
    real(r8), target, intent(in) :: ftem_prePBL_local(pcols,pver)
    real(r8), target, intent(in) :: t_aftPBL_local(pcols,pver)
    real(r8), target, intent(in) :: ftem_local(pcols,pver)
    real(r8), intent(in) :: rztodt_local
    real(r8), target, intent(inout) :: ftem_aftPBL_local(pcols,pver)
    real(r8), target, intent(inout) :: tten_local(pcols,pver)
    real(r8), target, intent(inout) :: rhten_local(pcols,pver)

    interface
       subroutine vertical_diffusion_post_qsat_diag_codon(ncol_c, pcols_c, pver_c, rztodt_c, state_t_p, &
            qv_aft_PBL_p, ftem_prePBL_p, t_aftPBL_p, ftem_p, ftem_aftPBL_p, tten_p, rhten_p) &
            bind(c, name="vertical_diffusion_post_qsat_diag_codon")
         use iso_c_binding, only: c_int64_t, c_double, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         real(c_double), value :: rztodt_c
         type(c_ptr), value :: state_t_p, qv_aft_PBL_p, ftem_prePBL_p, t_aftPBL_p, ftem_p
         type(c_ptr), value :: ftem_aftPBL_p, tten_p, rhten_p
       end subroutine vertical_diffusion_post_qsat_diag_codon
    end interface

    call vertical_diffusion_diag_batch_select_impl()

    if (use_native_diag_batch_impl) then
       call vertical_diffusion_post_qsat_diag_native(ncol, state_t_local, qv_aft_PBL_local, ftem_prePBL_local, &
            t_aftPBL_local, ftem_local, rztodt_local, ftem_aftPBL_local, tten_local, rhten_local)
       return
    end if

    call vertical_diffusion_diag_batch_call(4, ncol, rztodt_local=rztodt_local, state_t_local=state_t_local, &
         qv_aft_PBL_local=qv_aft_PBL_local, ftem_prePBL_local=ftem_prePBL_local, t_aftPBL_local=t_aftPBL_local, &
         ftem_local=ftem_local, ftem_aftPBL_local=ftem_aftPBL_local, tten_local=tten_local, rhten_local=rhten_local)

  end subroutine vertical_diffusion_post_qsat_diag

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_post_qsat_diag_native(ncol, state_t_local, qv_aft_PBL_local, ftem_prePBL_local, &
       t_aftPBL_local, ftem_local, rztodt_local, ftem_aftPBL_local, tten_local, rhten_local)

    integer, intent(in) :: ncol
    real(r8), intent(in) :: state_t_local(pcols,pver)
    real(r8), intent(in) :: qv_aft_PBL_local(pcols,pver)
    real(r8), intent(in) :: ftem_prePBL_local(pcols,pver)
    real(r8), intent(in) :: t_aftPBL_local(pcols,pver)
    real(r8), intent(in) :: ftem_local(pcols,pver)
    real(r8), intent(in) :: rztodt_local
    real(r8), intent(inout) :: ftem_aftPBL_local(pcols,pver)
    real(r8), intent(inout) :: tten_local(pcols,pver)
    real(r8), intent(inout) :: rhten_local(pcols,pver)

    ftem_aftPBL_local(:ncol,:pver) = qv_aft_PBL_local(:ncol,:pver) / ftem_local(:ncol,:pver) * 100._r8
    tten_local(:ncol,:pver)        = ( t_aftPBL_local(:ncol,:pver)    - state_t_local(:ncol,:pver) )     * rztodt_local
    rhten_local(:ncol,:pver)       = ( ftem_aftPBL_local(:ncol,:pver) - ftem_prePBL_local(:ncol,:pver) ) * rztodt_local

  end subroutine vertical_diffusion_post_qsat_diag_native

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vertical_diffusion_tend( &
                                      ztodt    , state    ,                  &
                                      taux     , tauy     , shflx    , cflx, &
                                      ustar    , obklen   , ptend    , &
                                      cldn     , ocnfrac  , landfrac , sgh      , pbuf) 
    !---------------------------------------------------- !
    ! This is an interface routine for vertical diffusion !
    !---------------------------------------------------- !
    use physics_buffer,     only : physics_buffer_desc, pbuf_get_field, pbuf_set_field
    use physics_types,      only : physics_state, physics_ptend, physics_ptend_init
    use cam_history,        only : outfld
    
    use trb_mtn_stress,     only : compute_tms
    use eddy_diff,          only : compute_eddy_diff
    use hb_diff,            only : compute_hb_diff
    use wv_saturation,      only : qsat
    use molec_diff,         only : compute_molec_diff, vd_lu_qdecomp
    use constituents,       only : qmincg, qmin
    use diffusion_solver,   only : compute_vdiff, any, operator(.not.)
    use physconst,          only : cpairv, rairv !Needed for calculation of upward H flux
    use time_manager,       only : get_nstep
    use constituents,       only : cnst_get_type_byind, cnst_name, cnst_fixed_ubc, cnst_fixed_ubflx
    use physconst,          only : pi
    use pbl_utils,          only : virtem, calc_obklen

    ! --------------- !
    ! Input Arguments !
    ! --------------- !

    type(physics_state), intent(in)    :: state                     ! Physics state variables

    real(r8),            intent(in)    :: taux(pcols)               ! x surface stress  [ N/m2 ]
    real(r8),            intent(in)    :: tauy(pcols)               ! y surface stress  [ N/m2 ]
    real(r8),            intent(in)    :: shflx(pcols)              ! Surface sensible heat flux  [ w/m2 ]
    real(r8),            intent(in)    :: cflx(pcols,pcnst)         ! Surface constituent flux [ kg/m2/s ]
    real(r8),            intent(in)    :: ztodt                     ! 2 delta-t [ s ]
    real(r8),            intent(in)    :: cldn(pcols,pver)          ! New stratus fraction [ fraction ]
    real(r8),            intent(in)    :: ocnfrac(pcols)            ! Ocean fraction
    real(r8),            intent(in)    :: landfrac(pcols)           ! Land fraction
    real(r8),            intent(in)    :: sgh(pcols)                ! Standard deviation of orography [ unit ? ]

    ! ---------------------- !
    ! Input-Output Arguments !
    ! ---------------------- !
    
    type(physics_ptend), intent(out) :: ptend                       ! Individual parameterization tendencies
    type(physics_buffer_desc), pointer :: pbuf(:)

    ! ---------------- !
    ! Output Arguments !
    ! ---------------- !

    real(r8),            intent(out)   :: ustar(pcols)              ! Surface friction velocity [ m/s ]
    real(r8),            intent(out)   :: obklen(pcols)             ! Obukhov length [ m ]

    ! --------------- !
    ! Local Variables !
    ! --------------- !

    logical  :: kvinit                                              ! Tell compute_eddy_diff/ caleddy
                                                                    ! to initialize kvh, kvm (uses kvf)
    character(128) :: errstring                                     ! Error status for compute_vdiff
    real(r8), pointer, dimension(:,:) :: qrl                        ! LW radiative cooling rate
    real(r8), pointer, dimension(:,:) :: wsedl                      ! Sedimentation velocity
                                                                    ! of stratiform liquid cloud droplet [ m/s ]

    integer  :: lchnk                                               ! Chunk identifier
    integer  :: ncol                                                ! Number of atmospheric columns
    integer  :: i, k, l, m                                          ! column, level, constituent indices
    integer  :: ierr                                                ! status for allocate/deallocate

    real(r8) :: dtk(pcols,pver)                                     ! T tendency from KE dissipation
    real(r8), pointer   :: tke(:,:)                                 ! Turbulent kinetic energy [ m2/s2 ]
    integer(i4),pointer :: turbtype(:,:)                            ! Turbulent interface types [ no unit ]
    real(r8), pointer   :: smaw(:,:)                                ! Normalized Galperin instability function
                                                                    ! ( 0<= <=4.964 and 1 at neutral )

    real(r8), pointer   :: qtl_flx(:,:)                             ! overbar(w'qtl') where qtl = qv + ql
    real(r8), pointer   :: qti_flx(:,:)                             ! overbar(w'qti') where qti = qv + qi

    real(r8) :: cgs(pcols,pverp)                                    ! Counter-gradient star  [ cg/flux ]
    real(r8) :: cgh(pcols,pverp)                                    ! Counter-gradient term for heat
    real(r8) :: rztodt                                              ! 1./ztodt [ 1/s ]
    real(r8) :: ksrftms(pcols)                                      ! Turbulent mountain stress surface drag coefficient [ kg/s/m2 ]
    real(r8) :: tautmsx(pcols)                                      ! U component of turbulent mountain stress [ N/m2 ]
    real(r8) :: tautmsy(pcols)                                      ! V component of turbulent mountain stress [ N/m2 ]
    real(r8) :: tautotx(pcols)                                      ! U component of total surface stress [ N/m2 ]
    real(r8) :: tautoty(pcols)                                      ! V component of total surface stress [ N/m2 ]

    real(r8), pointer :: kvh_in(:,:)                                ! kvh from previous timestep [ m2/s ]
    real(r8), pointer :: kvm_in(:,:)                                ! kvm from previous timestep [ m2/s ]
    real(r8), pointer :: kvt(:,:)                                   ! Molecular kinematic conductivity for temperature [  ]
    real(r8) :: kvq(pcols,pverp)                                    ! Eddy diffusivity for constituents [ m2/s ]
    real(r8) :: kvh(pcols,pverp)                                    ! Eddy diffusivity for heat [ m2/s ]
    real(r8) :: kvm(pcols,pverp)                                    ! Eddy diffusivity for momentum [ m2/s ]
    real(r8), pointer :: bprod(:,:)                                 ! Buoyancy production of tke [ m2/s3 ]
    real(r8) :: sprod(pcols,pverp)                                  ! Shear production of tke [ m2/s3 ]
    real(r8) :: sfi(pcols,pverp)                                    ! Saturation fraction at interfaces [ fraction ]
    real(r8) :: sl(pcols,pver)
    real(r8) :: qt(pcols,pver)
    real(r8) :: slv(pcols,pver)
    real(r8) :: sl_prePBL(pcols,pver)
    real(r8) :: qt_prePBL(pcols,pver)
    real(r8) :: slv_prePBL(pcols,pver)
    real(r8) :: slten(pcols,pver)
    real(r8) :: qtten(pcols,pver)
    real(r8) :: slvten(pcols,pver)
    real(r8) :: slflx(pcols,pverp)
    real(r8) :: qtflx(pcols,pverp)
    real(r8) :: uflx(pcols,pverp)
    real(r8) :: vflx(pcols,pverp)
    real(r8) :: slflx_cg(pcols,pverp)
    real(r8) :: qtflx_cg(pcols,pverp)
    real(r8) :: uflx_cg(pcols,pverp)
    real(r8) :: vflx_cg(pcols,pverp)
    real(r8) :: th(pcols,pver)                                      ! Potential temperature
    real(r8) :: topflx(pcols)                                       ! Molecular heat flux at top interface
    real(r8) :: wpert(pcols)                                        ! Turbulent wind gusts
    real(r8) :: rhoair

    real(r8) :: ri(pcols,pver)                                      ! richardson number (HB output)
    
    ! for obklen calculation outside HB
    real(r8) :: thvs(pcols)                                         ! Virtual potential temperature at surface
    real(r8) :: rrho(pcols)                                         ! Reciprocal of density at surface
    real(r8) :: khfs(pcols)                                         ! sfc kinematic heat flux [mK/s]
    real(r8) :: kqfs(pcols)                                         ! sfc kinematic water vapor flux [m/s]
    real(r8) :: kbfs(pcols)                                         ! sfc kinematic buoyancy flux [m^2/s^3]

    real(r8) :: ftem(pcols,pver)                                    ! Saturation vapor pressure before PBL
    real(r8) :: ftem_prePBL(pcols,pver)                             ! Saturation vapor pressure before PBL
    real(r8) :: ftem_aftPBL(pcols,pver)                             ! Saturation vapor pressure after PBL
    real(r8) :: tem2(pcols,pver)                                    ! Saturation specific humidity and RH
    real(r8) :: t_aftPBL(pcols,pver)                                ! Temperature after PBL diffusion
    real(r8) :: tten(pcols,pver)                                    ! Temperature tendency by PBL diffusion
    real(r8) :: rhten(pcols,pver)                                   ! RH tendency by PBL diffusion
    real(r8) :: qv_aft_PBL(pcols,pver)                              ! qv after PBL diffusion
    real(r8) :: ql_aft_PBL(pcols,pver)                              ! ql after PBL diffusion
    real(r8) :: qi_aft_PBL(pcols,pver)                              ! qi after PBL diffusion
    real(r8) :: s_aft_PBL(pcols,pver)                               ! s after PBL diffusion
    real(r8) :: u_aft_PBL(pcols,pver)                               ! u after PBL diffusion
    real(r8) :: v_aft_PBL(pcols,pver)                               ! v after PBL diffusion
    real(r8) :: qv_pro(pcols,pver) 
    real(r8) :: ql_pro(pcols,pver)
    real(r8) :: qi_pro(pcols,pver)
    real(r8) :: s_pro(pcols,pver)
    real(r8) :: t_pro(pcols,pver)
    real(r8), pointer :: tauresx(:)                                      ! Residual stress to be added in vdiff to correct
    real(r8), pointer :: tauresy(:)                                      ! for turb stress mismatch between sfc and atm accumulated.
    integer(i4), pointer :: ipbl(:)
    integer(i4), pointer :: kpblh(:)
    real(r8), pointer :: wstarPBL(:)
    real(r8), pointer :: tkes(:)
    real(r8), pointer :: went(:)
    real(r8) :: tpertPBL(pcols)
    real(r8) :: qpertPBL(pcols)
    real(r8) :: rairi(pcols,pver+1)                                 ! interface gas constant needed for compute_vdiff

    real(r8), pointer :: tpert(:)
    real(r8), pointer :: qpert(:)
    real(r8), pointer :: pblh(:)

    real(r8) :: tmp1(pcols)                                         ! Temporary storage

    integer  :: nstep
    real(r8) :: sum1, sum2, sum3, pdelx 
    real(r8) :: sflx

    ! Copy state so we can pass to intent(inout) routines that return
    ! new state instead of a tendency.
    real(r8) :: s_tmp(pcols,pver)
    real(r8) :: u_tmp(pcols,pver)
    real(r8) :: v_tmp(pcols,pver)
    real(r8) :: q_tmp(pcols,pver,pcnst)

    logical  :: lq(pcnst)
    logical  :: do_tms_local
    logical  :: do_molec_diff_local
    logical  :: use_diag_tke_local
    logical  :: use_hb_family_local
    logical  :: shallow_unicon_local
    logical  :: prog_modal_aero_local
    logical  :: do_pseudocon_diff_local
    logical  :: diff_cnsrv_mass_check_local
    logical  :: waccmx_special_local

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !

    rztodt = 1._r8 / ztodt
    lchnk  = state%lchnk
    ncol   = state%ncol

    call vertical_diffusion_tend_select_impl()
    if (.not. use_native_tend_impl) then
       call vertical_diffusion_tend_select_branches( &
            do_tms, do_molec_diff, trim(eddy_scheme) == 'diag_TKE', &
            trim(eddy_scheme) == 'HB' .or. trim(eddy_scheme) == 'HBR', &
            trim(shallow_scheme) == 'UNICON', prog_modal_aero, do_pseudocon_diff, &
            diff_cnsrv_mass_check, waccmx_is('ionosphere') .or. waccmx_is('neutral') )
       do_tms_local = iand(tend_branch_mask, 1) /= 0
       do_molec_diff_local = iand(tend_branch_mask, 2) /= 0
       use_diag_tke_local = iand(tend_branch_mask, 4) /= 0
       use_hb_family_local = iand(tend_branch_mask, 8) /= 0
       shallow_unicon_local = iand(tend_branch_mask, 16) /= 0
       prog_modal_aero_local = iand(tend_branch_mask, 32) /= 0
       do_pseudocon_diff_local = iand(tend_branch_mask, 64) /= 0
       diff_cnsrv_mass_check_local = iand(tend_branch_mask, 128) /= 0
       waccmx_special_local = iand(tend_branch_mask, 256) /= 0
    else
       do_tms_local = do_tms
       do_molec_diff_local = do_molec_diff
       use_diag_tke_local = trim(eddy_scheme) == 'diag_TKE'
       use_hb_family_local = trim(eddy_scheme) == 'HB' .or. trim(eddy_scheme) == 'HBR'
       shallow_unicon_local = trim(shallow_scheme) == 'UNICON'
       prog_modal_aero_local = prog_modal_aero
       do_pseudocon_diff_local = do_pseudocon_diff
       diff_cnsrv_mass_check_local = diff_cnsrv_mass_check
       waccmx_special_local = waccmx_is('ionosphere') .or. waccmx_is('neutral')
    end if

    call pbuf_get_field(pbuf, tauresx_idx,  tauresx)
    call pbuf_get_field(pbuf, tauresy_idx,  tauresy)
    call pbuf_get_field(pbuf, tpert_idx,    tpert)
    call pbuf_get_field(pbuf, qpert_idx,    qpert)
    call pbuf_get_field(pbuf, pblh_idx,     pblh)
    call pbuf_get_field(pbuf, turbtype_idx, turbtype)

    ! ---------------------------------------- !
    ! Computation of turbulent mountain stress !
    ! ---------------------------------------- !
   
    ! Consistent with the computation of 'normal' drag coefficient, we are using 
    ! the raw input (u,v) to compute 'ksrftms', not the provisionally-marched 'u,v' 
    ! within the iteration loop of the PBL scheme. 

    if( do_tms_local ) then
        call compute_tms( pcols      , pver     , ncol    ,              &
                          state%u    , state%v  , state%t , state%pmid , & 
                          state%exner, state%zm , sgh     , ksrftms    , & 
                          tautmsx    , tautmsy  , landfrac )
      ! Here, both 'taux, tautmsx' are explicit surface stresses.        
      ! Note that this 'tautotx, tautoty' are different from the total stress
      ! that has been actually added into the atmosphere. This is because both
      ! taux and tautmsx are fully implicitly treated within compute_vdiff.
      ! However, 'tautotx, tautoty' are not used in the actual numerical
      ! computation in this module.   
        tautotx(:ncol) = taux(:ncol) + tautmsx(:ncol)
        tautoty(:ncol) = tauy(:ncol) + tautmsy(:ncol)
    else
        ksrftms(:ncol) = 0._r8
        tautotx(:ncol) = taux(:ncol)
        tautoty(:ncol) = tauy(:ncol)
    endif

    !----------------------------------------------------------------------- !
    !   Computation of eddy diffusivities - Select appropriate PBL scheme    !
    !----------------------------------------------------------------------- !
    call pbuf_get_field(pbuf, kvm_idx,  kvm_in)
    call pbuf_get_field(pbuf, kvh_idx,  kvh_in)
    call pbuf_get_field(pbuf, smaw_idx, smaw)
    call pbuf_get_field(pbuf, tke_idx,  tke)

    if (use_diag_tke_local) then

       ! ---------------------------------------------------------------- !
       ! At first time step, have eddy_diff.F90:caleddy() use kvh=kvm=kvf !
       ! This has to be done in compute_eddy_diff after kvf is calculated !
       ! ---------------------------------------------------------------- !

       if( is_first_step() ) then
           kvinit = .true.
       else
           kvinit = .false.
       endif

       ! ---------------------------------------------- !
       ! Get LW radiative heating out of physics buffer !
       ! ---------------------------------------------- !
       call pbuf_get_field(pbuf, qrl_idx,   qrl)

       call pbuf_get_field(pbuf, wsedl_idx, wsedl)

       ! These fields are put into the pbuf for UNICON only.
       if (shallow_unicon_local) then
          call pbuf_get_field(pbuf, bprod_idx,    bprod)
          call pbuf_get_field(pbuf, ipbl_idx,     ipbl)
          call pbuf_get_field(pbuf, kpblh_idx,    kpblh)
          call pbuf_get_field(pbuf, wstarPBL_idx, wstarPBL)
          call pbuf_get_field(pbuf, tkes_idx,     tkes)
          call pbuf_get_field(pbuf, went_idx,     went)
       else
          allocate(bprod(pcols,pverp), ipbl(pcols), kpblh(pcols), wstarPBL(pcols), tkes(pcols), went(pcols))
       end if

       call compute_eddy_diff( lchnk    ,                                                                    &
                               pcols    , pver        , ncol       , state%t    , state%q(:,:,1) , ztodt   , &
                               state%q(:,:,ixcldliq)  , state%q(:,:,ixcldice)   ,                            &
                               state%s  , state%pdel  , state%rpdel, cldn       , qrl            , wsedl   , &
                               state%zm , state%zi    , state%pmid , state%pint , state%u        , state%v , &
                               taux     , tauy        , shflx      , cflx(:,1)  , wstarent       , nturb   , &
                               rrho     , ustar       , pblh       , kvm_in     , kvh_in         , kvm     , &
                               kvh      , kvq         , cgh        ,                                         &
                               cgs      , tpert       , qpert      , wpert      , tke            , bprod   , &
                               sprod    , sfi         , kvinit     ,                                         &
                               tauresx  , tauresy     , ksrftms    ,                                         &
                               ipbl     , kpblh       , wstarPBL   , tkes       , went           , turbtype, &
                               smaw )

       ! The diag_TKE scheme does not calculate the Monin-Obukhov length, which is used in dry deposition calculations.
       ! Use the routines from pbl_utils to accomplish this. Assumes ustar and rrho have been set.
       call vertical_diffusion_obklen_diag(ncol, state%t, state%exner, state%q, cflx, shflx, rrho, ustar, th, thvs, &
            khfs, kqfs, kbfs, obklen)

       ! ----------------------------------------------- !       
       ! Store TKE in pbuf for use by shallow convection !
       ! ----------------------------------------------- !   

       tpertPBL(:ncol) = tpert(:ncol)
       qpertPBL(:ncol) = qpert(:ncol)


       ! The diffusivities from diag_TKE can be much larger than from HB in the free
       ! troposphere and upper atmosphere. These seem to be larger than observations,
       ! and in WACCM the gw_drag code is already applying an eddy diffusivity in the
       ! upper atmosphere. Optionally, adjust the diffusivities in the free troposphere
       ! or the upper atmosphere.
       !
       ! NOTE: Further investigation should be done as to why the diffusivities are
       ! larger in diag_TKE.
       if ((kv_freetrop_scale /= 1._r8) .or. ((kv_top_scale /= 1._r8) .and. (kv_top_pressure > 0._r8))) then
         do i = 1, ncol
           do k = 1, pverp
           
             ! Outside of the boundary layer?
             if (state%zi(i,k) > pblh(i)) then

               ! In the upper atmosphere?
               if (state%pint(i,k) <= kv_top_pressure) then
                 kvh(i,k) = kvh(i,k) * kv_top_scale
                 kvm(i,k) = kvm(i,k) * kv_top_scale
                 kvq(i,k) = kvq(i,k) * kv_top_scale
               else
                 kvh(i,k) = kvh(i,k) * kv_freetrop_scale
                 kvm(i,k) = kvm(i,k) * kv_freetrop_scale
                 kvq(i,k) = kvq(i,k) * kv_freetrop_scale
               end if
             else
               exit
             end if
           end do
         end do
       end if

       ! Write out fields that are only used by this scheme

       call outfld( 'BPROD   ', bprod, pcols, lchnk )
       call outfld( 'SPROD   ', sprod, pcols, lchnk )
       call outfld( 'SFI     ', sfi,   pcols, lchnk )

       if (.not. shallow_unicon_local) then
          deallocate(bprod, ipbl, kpblh, wstarPBL, tkes, went)
       end if

    else if (use_hb_family_local) then

       ! Modification : We may need to use 'taux' instead of 'tautotx' here, for
       !                consistency with the previous HB scheme.

       th(:ncol,:pver) = state%t(:ncol,:pver) * state%exner(:ncol,:pver)

       call compute_hb_diff( lchnk     , ncol     ,                                &
                             th        , state%t  , state%q , state%zm , state%zi, &
                             state%pmid, state%u  , state%v , tautotx  , tautoty , &
                             shflx     , cflx(:,1), obklen  , ustar    , pblh    , &
                             kvm       , kvh      , kvq     , cgh      , cgs     , &
                             tpert     , qpert    , cldn    , ocnfrac  , tke     , &
                             ri        , &
                             eddy_scheme )
    
       call outfld( 'HB_ri',          ri,         pcols,   lchnk )

       wpert = 0._r8  

    end if

    call outfld( 'ustar',   ustar(:), pcols, lchnk )
    call outfld( 'obklen', obklen(:), pcols, lchnk )

    ! kvh (in pbuf) is used by other physics parameterizations, and as an initial guess in compute_eddy_diff
    ! on the next timestep.  It is not updated by the compute_vdiff call below.
    call pbuf_set_field(pbuf, kvh_idx, kvh)

    ! kvm (in pbuf) is only used as an initial guess in compute_eddy_diff on the next timestep.
    ! The contributions for molecular diffusion made to kvm by the call to compute_vdiff below 
    ! are not included in the pbuf as these are not needed in the initial guess by compute_eddy_diff.
    call pbuf_set_field(pbuf, kvm_idx, kvm)

    call outfld( 'WGUSTD' , wpert, pcols, lchnk )

    !------------------------------------ ! 
    !    Application of diffusivities     !
    !------------------------------------ !

    !------------------------------------------------------ !
    ! Write profile output before applying diffusion scheme !
    !------------------------------------------------------ !

    call vertical_diffusion_pre_pbl_diag(ncol, state%q, state%s, state%u, state%v, q_tmp, s_tmp, u_tmp, v_tmp, &
         sl_prePBL, qt_prePBL, slv_prePBL)

    call qsat(state%t(:ncol,:), state%pmid(:ncol,:), &
         tem2(:ncol,:), ftem(:ncol,:))
    call vertical_diffusion_pre_qsat_rh(ncol, state%q, ftem, ftem_prePBL)

    call outfld( 'qt_pre_PBL   ', qt_prePBL,                 pcols, lchnk )
    call outfld( 'sl_pre_PBL   ', sl_prePBL,                 pcols, lchnk )
    call outfld( 'slv_pre_PBL  ', slv_prePBL,                pcols, lchnk )
    call outfld( 'u_pre_PBL    ', state%u,                   pcols, lchnk )
    call outfld( 'v_pre_PBL    ', state%v,                   pcols, lchnk )
    call outfld( 'qv_pre_PBL   ', state%q(:ncol,:,1),        pcols, lchnk )
    call outfld( 'ql_pre_PBL   ', state%q(:ncol,:,ixcldliq), pcols, lchnk )
    call outfld( 'qi_pre_PBL   ', state%q(:ncol,:,ixcldice), pcols, lchnk )
    call outfld( 't_pre_PBL    ', state%t,                   pcols, lchnk )
    call outfld( 'rh_pre_PBL   ', ftem_prePBL,               pcols, lchnk )

    ! --------------------------------------------------------------------------------- !
    ! Call the diffusivity solver and solve diffusion equation                          !
    ! The final two arguments are optional function references to                       !
    ! constituent-independent and constituent-dependent moleculuar diffusivity routines !
    ! --------------------------------------------------------------------------------- !

    !------------------------------------------------------------------------ 
    !  Check to see if constituent dependent gas constant needed (WACCM-X) 
    !------------------------------------------------------------------------
    if (waccmx_special_local) then 
      rairi(:ncol,1) = rairv(:ncol,1,lchnk)
      do k = 2, pver
        do i = 1, ncol
          rairi(i,k) = 0.5_r8 * (rairv(i,k,lchnk)+rairv(i,k-1,lchnk))
        end do
      end do
    else
      rairi(:ncol,:pver+1) = rair 
    endif

    ! Modification : We may need to output 'tautotx_im,tautoty_im' from below 'compute_vdiff' and
    !                separately print out as diagnostic output, because these are different from
    !                the explicit 'tautotx, tautoty' computed above. 
    ! Note that the output 'tauresx,tauresy' from below subroutines are fully implicit ones.

    call pbuf_get_field(pbuf, kvt_idx, kvt)

    if( any(fieldlist_wet) ) then

        call compute_vdiff( state%lchnk   ,                                                                     &
                            pcols         , pver               , pcnst        , ncol          , state%pmid    , &
                            state%pint    , state%pdel         , &
                            state%rpdel        , state%t      , ztodt         , taux          , &
                            tauy          , shflx              , cflx         , ntop          , nbot          , &
                            kvh           , kvm                , kvq          , cgs           , cgh           , &
                            state%zi      , ksrftms            , qmincg       , fieldlist_wet , fieldlist_molec,&
                            u_tmp         , v_tmp              , q_tmp        , s_tmp         ,                 &
                            tautmsx       , tautmsy            , dtk          , topflx        , errstring     , &
                            tauresx       , tauresy            , 1            , cpairv(:,:,state%lchnk), rairi, &
                            do_molec_diff_local , compute_molec_diff , vd_lu_qdecomp, kvt )

        call handle_errmsg(errstring, subname="compute_vdiff", &
             extra_msg="Error in fieldlist_wet call from vertical_diffusion.")

    end if
 
    if( any( fieldlist_dry ) ) then

        if( do_molec_diff_local ) then
            errstring = "Design flaw: dry vdiff not currently supported with molecular diffusion"
            call endrun(errstring)
        end if

        call compute_vdiff( state%lchnk   ,                                                                     &
                            pcols         , pver               , pcnst        , ncol          , state%pmiddry , &
                            state%pintdry , state%pdeldry      , &
                            state%rpdeldry     , state%t      , ztodt         , taux          , &       
                            tauy          , shflx              , cflx         , ntop          , nbot          , &       
                            kvh           , kvm                , kvq          , cgs           , cgh           , &   
                            state%zi      , ksrftms            , qmincg       , fieldlist_dry , fieldlist_molec,&
                            u_tmp         , v_tmp              , q_tmp        , s_tmp         ,                 &
                            tautmsx       , tautmsy            , dtk          , topflx        , errstring     , &
                            tauresx       , tauresy            , 1            , cpairv(:,:,state%lchnk), rairi, &
                            do_molec_diff_local , compute_molec_diff , vd_lu_qdecomp )

        call handle_errmsg(errstring, subname="compute_vdiff", &
             extra_msg="Error in fieldlist_dry call from vertical_diffusion.")

    end if

    if (prog_modal_aero_local) then

       ! Modal aerosol species not diffused, so just add the explicit surface fluxes to the
       ! lowest layer

       call vertical_diffusion_core_batch_select_impl()
       if (use_native_core_batch_impl) then
          call vertical_diffusion_modal_aero_flux(ncol, state%rpdel, cflx, ztodt, q_tmp)
       else
          call vertical_diffusion_core_batch_call(1, ncol, ztodt_local=ztodt, state_rpdel_local=state%rpdel, &
               cflx_local=cflx, q_tmp_local=q_tmp)
       end if
    end if

    ! -------------------------------------------------------- !
    ! Diagnostics and output writing after applying PBL scheme !
    ! -------------------------------------------------------- !

    call vertical_diffusion_core_batch_select_impl()
    if (use_native_core_batch_impl) then
       call vertical_diffusion_flux_diag(ncol, q_tmp, s_tmp, u_tmp, v_tmp, state%pint, state%zi, state%zm, &
            cflx, kvh, kvm, cgs, cgh, shflx, tautotx, tautoty, sl, qt, slv, slflx, qtflx, uflx, vflx, &
            slflx_cg, qtflx_cg, uflx_cg, vflx_cg)
    else
       call vertical_diffusion_core_batch_call(2, ncol, q_tmp_local=q_tmp, s_tmp_local=s_tmp, u_tmp_local=u_tmp, &
            v_tmp_local=v_tmp, pint_local=state%pint, zi_local=state%zi, zm_local=state%zm, cflx_local=cflx, &
            kvh_local=kvh, kvm_local=kvm, cgs_local=cgs, cgh_local=cgh, shflx_local=shflx, &
            tautotx_local=tautotx, tautoty_local=tautoty, sl_local=sl, qt_local=qt, slv_local=slv, &
            slflx_local=slflx, qtflx_local=qtflx, uflx_local=uflx, vflx_local=vflx, &
            slflx_cg_local=slflx_cg, qtflx_cg_local=qtflx_cg, uflx_cg_local=uflx_cg, vflx_cg_local=vflx_cg)
    end if

    if (shallow_unicon_local) then
       call pbuf_get_field(pbuf, qtl_flx_idx,  qtl_flx)
       call pbuf_get_field(pbuf, qti_flx_idx,  qti_flx)
       qtl_flx(:ncol,1) = 0._r8
       qti_flx(:ncol,1) = 0._r8
       do k = 2, pver
          do i = 1, ncol
             ! For use in the cloud macrophysics
             ! Note that density is not addd here. Also, only consider local transport term.
             qtl_flx(i,k) = - kvh(i,k)*(q_tmp(i,k-1,1)-q_tmp(i,k,1)+q_tmp(i,k-1,ixcldliq)-q_tmp(i,k,ixcldliq))/&
                                       (state%zm(i,k-1)-state%zm(i,k))
             qti_flx(i,k) = - kvh(i,k)*(q_tmp(i,k-1,1)-q_tmp(i,k,1)+q_tmp(i,k-1,ixcldice)-q_tmp(i,k,ixcldice))/&
                                       (state%zm(i,k-1)-state%zm(i,k))
          end do
       end do
       do i = 1, ncol
          rhoair = state%pint(i,pverp)/(rair*((slv(i,pver)-gravit*state%zi(i,pverp))/cpair))
          qtl_flx(i,pverp) = cflx(i,1)/rhoair
          qti_flx(i,pverp) = cflx(i,1)/rhoair
       end do
    end if

    ! --------------------------------------------------------------- !
    ! Convert the new profiles into vertical diffusion tendencies.    !
    ! Convert KE dissipative heat change into "temperature" tendency. !
    ! --------------------------------------------------------------- !

    ! All variables are modified by vertical diffusion

    lq(:) = .TRUE.
    call physics_ptend_init(ptend,state%psetcols, "vertical diffusion", &
         ls=.true., lu=.true., lv=.true., lq=lq)

    call vertical_diffusion_core_batch_select_impl()
    if (use_native_core_batch_impl) then
       call vertical_diffusion_ptend_core(ncol, state%psetcols, q_tmp, s_tmp, u_tmp, v_tmp, state%q, state%s, &
            state%u, state%v, sl, qt, sl_prePBL, qt_prePBL, rztodt, ptend%q, ptend%s, ptend%u, ptend%v, slten, &
            qtten)
    else
       call vertical_diffusion_core_batch_call(3, ncol, psetcols_local=state%psetcols, rztodt_local=rztodt, &
            q_tmp_local=q_tmp, s_tmp_local=s_tmp, u_tmp_local=u_tmp, v_tmp_local=v_tmp, state_q_local=state%q, &
            state_s_local=state%s, state_u_local=state%u, state_v_local=state%v, sl_local=sl, qt_local=qt, &
            sl_prePBL_local=sl_prePBL, qt_prePBL_local=qt_prePBL, ptend_q_local=ptend%q, ptend_s_local=ptend%s, &
            ptend_u_local=ptend%u, ptend_v_local=ptend%v, slten_local=slten, qtten_local=qtten)
    end if

    ! ----------------------------------------------------------- !
    ! In order to perform 'pseudo-conservative varible diffusion' !
    ! perform the following two stages:                           !
    !                                                             !
    ! I.  Re-set (1) 'qvten' by 'qtten', and 'qlten = qiten = 0'  !
    !            (2) 'sten'  by 'slten', and                      !
    !            (3) 'qlten = qiten = 0'                          !
    !                                                             !
    ! II. Apply 'positive_moisture'                               !
    !                                                             !
    ! ----------------------------------------------------------- !

    if( use_diag_tke_local .and. do_pseudocon_diff_local ) then

         ptend%q(:ncol,:pver,1) = qtten(:ncol,:pver)
         ptend%s(:ncol,:pver)   = slten(:ncol,:pver)
         ptend%q(:ncol,:pver,ixcldliq) = 0._r8         
         ptend%q(:ncol,:pver,ixcldice) = 0._r8         
         ptend%q(:ncol,:pver,ixnumliq) = 0._r8         
         ptend%q(:ncol,:pver,ixnumice) = 0._r8         

         do i = 1, ncol
            do k = 1, pver
               qv_pro(i,k) = state%q(i,k,1)        + ptend%q(i,k,1)             * ztodt       
               ql_pro(i,k) = state%q(i,k,ixcldliq) + ptend%q(i,k,ixcldliq)      * ztodt
               qi_pro(i,k) = state%q(i,k,ixcldice) + ptend%q(i,k,ixcldice)      * ztodt              
               s_pro(i,k)  = state%s(i,k)          + ptend%s(i,k)               * ztodt
               t_pro(i,k)  = state%t(i,k)          + (1._r8/cpair)*ptend%s(i,k) * ztodt
            end do 
         end do 

         call positive_moisture( cpair, latvap, latvap+latice, ncol, pver, ztodt, qmin(1), qmin(ixcldliq), qmin(ixcldice),    &
                                 state%pdel(:ncol,pver:1:-1), qv_pro(:ncol,pver:1:-1), ql_pro(:ncol,pver:1:-1), &
                                 qi_pro(:ncol,pver:1:-1), t_pro(:ncol,pver:1:-1), s_pro(:ncol,pver:1:-1),       &
                                 ptend%q(:ncol,pver:1:-1,1), ptend%q(:ncol,pver:1:-1,ixcldliq),                 &
                                 ptend%q(:ncol,pver:1:-1,ixcldice), ptend%s(:ncol,pver:1:-1) )

    end if

    ! ----------------------------------------------------------------- !
    ! Re-calculate diagnostic output variables after vertical diffusion !
    ! ----------------------------------------------------------------- !
 
    call vertical_diffusion_post_pbl_state(ncol, state%psetcols, state%q, state%s, state%u, state%v, state%zm, &
         ptend%q, ptend%s, ptend%u, ptend%v, ztodt, qv_aft_PBL, ql_aft_PBL, qi_aft_PBL, s_aft_PBL, t_aftPBL, &
         u_aft_PBL, v_aft_PBL)

    call qsat(t_aftPBL(:ncol,:pver), state%pmid(:ncol,:pver), &
         tem2(:ncol,:pver), ftem(:ncol,:pver))
    call vertical_diffusion_post_qsat_diag(ncol, state%t, qv_aft_PBL, ftem_prePBL, t_aftPBL, ftem, rztodt, &
         ftem_aftPBL, tten, rhten)

    ! -------------------------------------------------------------- !
    ! mass conservation check.........
    ! -------------------------------------------------------------- !
    if (diff_cnsrv_mass_check_local) then

       ! Conservation check
       nstep = get_nstep()
       
       do m = 1, pcnst
          fixed_ubc: if ((.not.cnst_fixed_ubc(m)).and.(.not.cnst_fixed_ubflx(m))) then
             col_loop: do i = 1, ncol
                sum1 = 0._r8
                sum2 = 0._r8
                sum3 = 0._r8
                do k = 1, pver
                   if(cnst_get_type_byind(m).eq.'wet') then
                      pdelx = state%pdel(i,k)
                   else
                      pdelx = state%pdeldry(i,k)
                   endif
                   sum1 = sum1 + state%q(i,k,m)*pdelx/gravit                          ! total column
                   sum2 = sum2 +(state%q(i,k,m)+ptend%q(i,k,m)*ztodt)*pdelx/ gravit   ! total column after tendancy is applied
                   sum3 = sum3 +(               ptend%q(i,k,m)*ztodt)*pdelx/ gravit   ! rate of change in column
                enddo
                sum1 = sum1 + (cflx(i,m) * ztodt) ! add in surface flux (kg/m2)
                sflx = (cflx(i,m) * ztodt) 
                if (sum1>1.e-36_r8) then
                   if( abs((sum2-sum1)/sum1) .gt. 1.e-12_r8  ) then
                      write(iulog,'(a,a8,a,I4,2f8.3,5e25.16)') &
                           'MASSCHECK vert diff : nstep,lon,lat,mass1,mass2,sum3,sflx,rel-diff : ', &
                           trim(cnst_name(m)), ' : ', nstep, state%lon(i)*180._r8/pi, state%lat(i)*180._r8/pi, &
                           sum1, sum2, sum3, sflx, abs(sum2-sum1)/sum1
                      call endrun('vertical_diffusion_tend : mass not conserved' )
                   endif
                endif
             enddo col_loop
          endif fixed_ubc
       enddo
    endif

    ! -------------------------------------------------------------- !
    ! Writing state variables after PBL scheme for detailed analysis !
    ! -------------------------------------------------------------- !

    call outfld( 'sl_aft_PBL'   , sl,                        pcols, lchnk )
    call outfld( 'qt_aft_PBL'   , qt,                        pcols, lchnk )
    call outfld( 'slv_aft_PBL'  , slv,                       pcols, lchnk )
    call outfld( 'u_aft_PBL'    , u_aft_PBL,                 pcols, lchnk )
    call outfld( 'v_aft_PBL'    , v_aft_PBL,                 pcols, lchnk )
    call outfld( 'qv_aft_PBL'   , qv_aft_PBL,                pcols, lchnk )
    call outfld( 'ql_aft_PBL'   , ql_aft_PBL,                pcols, lchnk )
    call outfld( 'qi_aft_PBL'   , qi_aft_PBL,                pcols, lchnk )
    call outfld( 't_aft_PBL '   , t_aftPBL,                  pcols, lchnk )
    call outfld( 'rh_aft_PBL'   , ftem_aftPBL,               pcols, lchnk )
    call outfld( 'slflx_PBL'    , slflx,                     pcols, lchnk )
    call outfld( 'qtflx_PBL'    , qtflx,                     pcols, lchnk )
    call outfld( 'uflx_PBL'     , uflx,                      pcols, lchnk )
    call outfld( 'vflx_PBL'     , vflx,                      pcols, lchnk )
    call outfld( 'slflx_cg_PBL' , slflx_cg,                  pcols, lchnk )
    call outfld( 'qtflx_cg_PBL' , qtflx_cg,                  pcols, lchnk )
    call outfld( 'uflx_cg_PBL'  , uflx_cg,                   pcols, lchnk )
    call outfld( 'vflx_cg_PBL'  , vflx_cg,                   pcols, lchnk )
    call outfld( 'slten_PBL'    , slten,                     pcols, lchnk )
    call outfld( 'qtten_PBL'    , qtten,                     pcols, lchnk )
    call outfld( 'uten_PBL'     , ptend%u(:ncol,:),          pcols, lchnk )
    call outfld( 'vten_PBL'     , ptend%v(:ncol,:),          pcols, lchnk )
    call outfld( 'qvten_PBL'    , ptend%q(:ncol,:,1),        pcols, lchnk )
    call outfld( 'qlten_PBL'    , ptend%q(:ncol,:,ixcldliq), pcols, lchnk )
    call outfld( 'qiten_PBL'    , ptend%q(:ncol,:,ixcldice), pcols, lchnk )
    call outfld( 'tten_PBL'     , tten,                      pcols, lchnk )
    call outfld( 'rhten_PBL'    , rhten,                     pcols, lchnk )

    ! ------------------------------------------- !
    ! Writing the other standard output variables !
    ! ------------------------------------------- !

    call outfld( 'QT'           , qt,                        pcols, lchnk )
    call outfld( 'SL'           , sl,                        pcols, lchnk )
    call outfld( 'SLV'          , slv,                       pcols, lchnk )
    call outfld( 'SLFLX'        , slflx,                     pcols, lchnk )
    call outfld( 'QTFLX'        , qtflx,                     pcols, lchnk )
    call outfld( 'UFLX'         , uflx,                      pcols, lchnk )
    call outfld( 'VFLX'         , vflx,                      pcols, lchnk )
    call outfld( 'TKE'          , tke,                       pcols, lchnk )

    call outfld( 'PBLH'         , pblh,                      pcols, lchnk )
    call outfld( 'TPERT'        , tpert,                     pcols, lchnk )
    call outfld( 'QPERT'        , qpert,                     pcols, lchnk )
    call outfld( 'USTAR'        , ustar,                     pcols, lchnk )
    call outfld( 'KVH'          , kvh,                       pcols, lchnk )
    call outfld( 'KVT'          , kvt,                       pcols, lchnk )
    call outfld( 'KVM'          , kvm,                       pcols, lchnk )
    call outfld( 'CGS'          , cgs,                       pcols, lchnk )
    dtk(:ncol,:) = dtk(:ncol,:) / cpair              ! Normalize heating for history
    call outfld( 'DTVKE'        , dtk,                       pcols, lchnk )
    dtk(:ncol,:) = ptend%s(:ncol,:) / cpair          ! Normalize heating for history using dtk
    call outfld( 'DTV'          , dtk,                       pcols, lchnk ) 
    call outfld( 'DUV'          , ptend%u,                   pcols, lchnk )
    call outfld( 'DVV'          , ptend%v,                   pcols, lchnk )
    do m = 1, pcnst
       call outfld( vdiffnam(m) , ptend%q(1,1,m),            pcols, lchnk )
    end do
    if( do_tms_local ) then
      ! Here, 'tautmsx,tautmsy' are implicit 'tms' that have been actually
      ! added into the atmosphere.
        call outfld( 'TAUTMSX'  , tautmsx,                   pcols, lchnk )
        call outfld( 'TAUTMSY'  , tautmsy,                   pcols, lchnk )
    end if
    if( do_molec_diff_local ) then
        call outfld( 'TTPXMLC'  , topflx,                    pcols, lchnk )
    end if

    return
  end subroutine vertical_diffusion_tend

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine positive_moisture( cp, xlv, xls, ncol, mkx, dt, qvmin, qlmin, qimin, & 
                                dp, qv, ql, qi, t, s, qvten, qlten, qiten, sten )
  ! ------------------------------------------------------------------------------- !
  ! If any 'ql < qlmin, qi < qimin, qv < qvmin' are developed in any layer,         !
  ! force them to be larger than minimum value by (1) condensating water vapor      !
  ! into liquid or ice, and (2) by transporting water vapor from the very lower     !
  ! layer. '2._r8' is multiplied to the minimum values for safety.                  !
  ! Update final state variables and tendencies associated with this correction.    !
  ! If any condensation happens, update (s,t) too.                                  !
  ! Note that (qv,ql,qi,t,s) are final state variables after applying corresponding !
  ! input tendencies.                                                               !
  ! Be careful the order of k : '1': near-surface layer, 'mkx' : top layer          ! 
  ! ------------------------------------------------------------------------------- !
    implicit none
    integer,  intent(in)     :: ncol, mkx
    real(r8), intent(in)     :: cp, xlv, xls
    real(r8), intent(in)     :: dt, qvmin, qlmin, qimin
    real(r8), intent(in)     :: dp(ncol,mkx)
    real(r8), intent(inout)  :: qv(ncol,mkx), ql(ncol,mkx), qi(ncol,mkx), t(ncol,mkx), s(ncol,mkx)
    real(r8), intent(inout)  :: qvten(ncol,mkx), qlten(ncol,mkx), qiten(ncol,mkx), sten(ncol,mkx)
    integer   i, k
    real(r8)  dql, dqi, dqv, sum, aa, dum 

  ! Modification : I should check whether this is exactly same as the one used in
  !                shallow convection and cloud macrophysics.

    do i = 1, ncol
       do k = mkx, 1, -1    ! From the top to the 1st (lowest) layer from the surface
          dql        = max(0._r8,1._r8*qlmin-ql(i,k))
          dqi        = max(0._r8,1._r8*qimin-qi(i,k))
          qlten(i,k) = qlten(i,k) +  dql/dt
          qiten(i,k) = qiten(i,k) +  dqi/dt
          qvten(i,k) = qvten(i,k) - (dql+dqi)/dt
          sten(i,k)  = sten(i,k)  + xlv * (dql/dt) + xls * (dqi/dt)
          ql(i,k)    = ql(i,k) +  dql
          qi(i,k)    = qi(i,k) +  dqi
          qv(i,k)    = qv(i,k) -  dql - dqi
          s(i,k)     = s(i,k)  +  xlv * dql + xls * dqi
          t(i,k)     = t(i,k)  + (xlv * dql + xls * dqi)/cp
          dqv        = max(0._r8,1._r8*qvmin-qv(i,k))
          qvten(i,k) = qvten(i,k) + dqv/dt
          qv(i,k)    = qv(i,k)    + dqv
          if( k .ne. 1 ) then 
              qv(i,k-1)    = qv(i,k-1)    - dqv*dp(i,k)/dp(i,k-1)
              qvten(i,k-1) = qvten(i,k-1) - dqv*dp(i,k)/dp(i,k-1)/dt
          endif
          qv(i,k) = max(qv(i,k),qvmin)
          ql(i,k) = max(ql(i,k),qlmin)
          qi(i,k) = max(qi(i,k),qimin)
       end do
       ! Extra moisture used to satisfy 'qv(i,1)=qvmin' is proportionally 
       ! extracted from all the layers that has 'qv > 2*qvmin'. This fully
       ! preserves column moisture. 
       if( dqv .gt. 1.e-20_r8 ) then
           sum = 0._r8
           do k = 1, mkx
              if( qv(i,k) .gt. 2._r8*qvmin ) sum = sum + qv(i,k)*dp(i,k)
           enddo
           aa = dqv*dp(i,1)/max(1.e-20_r8,sum)
           if( aa .lt. 0.5_r8 ) then
               do k = 1, mkx
                  if( qv(i,k) .gt. 2._r8*qvmin ) then
                      dum        = aa*qv(i,k)
                      qv(i,k)    = qv(i,k) - dum
                      qvten(i,k) = qvten(i,k) - dum/dt
                  endif
               enddo 
           else 
               write(iulog,*) 'Full positive_moisture is impossible in vertical_diffusion'
           endif
       endif 
    end do
    return

  end subroutine positive_moisture

end module vertical_diffusion
