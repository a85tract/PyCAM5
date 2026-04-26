  module eddy_diff

  !--------------------------------------------------------------------------------- !
  !                                                                                  !
  ! The University of Washington Moist Turbulence Scheme to compute eddy diffusion   ! 
  ! coefficients associated with dry and moist turbulences in the whole              !
  ! atmospheric layers.                                                              !
  !                                                                                  !
  ! For detailed description of the code and its performances, see                   !
  !                                                                                  !
  ! 1.'A new moist turbulence parametrization in the Community Atmosphere Model'     !
  !    by Christopher S. Bretherton and Sungsu Park. J. Climate. 2009. 22. 3422-3448 !
  ! 2.'The University of Washington shallow convection and moist turbulence schemes  !
  !    and their impact on climate simulations with the Community Atmosphere Model'  !
  !    by Sungsu Park and Christopher S. Bretherton. J. Climate. 2009. 22. 3449-3469 !
  !                                                                                  !
  ! For questions on the scheme and code, send an email to                           !
  !     Sungsu Park      at sungsup@ucar.edu (tel: 303-497-1375)                     !
  !     Chris Bretherton at breth@washington.edu                                     !
  !                                                                                  ! 
  ! Developed by Chris Bretherton at the University of Washington, Seattle, WA.      !
  !              Sungsu Park      at the CGD/NCAR, Boulder, CO.                      !
  ! Last coded on May.2006, Dec.2009 by Sungsu Park.                                 !
  !                                                                                  !  
  !--------------------------------------------------------------------------------- !

  use diffusion_solver, only: vdiff_selector
  use cam_history,      only: outfld, addfld, phys_decomp
  use cam_logfile,      only: iulog
  use ppgrid,           only: pver  
  use cam_abortutils,   only: endrun
  use spmd_utils,       only: masterproc, iam
  use wv_saturation,    only: qsat, estblf, svp_to_qsat

  implicit none
  private
  save

  public init_eddy_diff
  public compute_eddy_diff

  type(vdiff_selector)        :: fieldlist_wet                  ! Logical switches for moist mixing ratio diffusion
  type(vdiff_selector)        :: fieldlist_molec                ! Logical switches for molecular diffusion
  integer,          parameter :: r8 = selected_real_kind(12)    ! 8 byte real
  integer,          parameter :: i4 = selected_int_kind( 6)     ! 4 byte integer
  ! --------------------------------- !
  ! PBL Parameters used in the UW PBL !
  ! --------------------------------- !

  character,        parameter :: sftype         = 'l'           ! Method for calculating saturation fraction

  character(len=4), parameter :: choice_evhc    = 'maxi'        ! 'orig', 'ramp', 'maxi' : recommended to be used with choice_radf 
  character(len=6), parameter :: choice_radf    = 'maxi'        ! 'orig', 'ramp', 'maxi' : recommended to be used with choice_evhc 
  character(len=6), parameter :: choice_SRCL    = 'nonamb'      ! 'origin', 'remove', 'nonamb'
 
  character(len=6), parameter :: choice_tunl    = 'rampcl'      ! 'origin', 'rampsl'(Sungsu), 'rampcl'(Chris)
  real(r8),         parameter :: ctunl          =  2._r8        !  Maximum asympt leng = ctunl*tunl when choice_tunl = 'rampsl(cl)'
                                                                ! [ no unit ]
  character(len=6), parameter :: choice_leng    = 'origin'      ! 'origin', 'takemn'
  real(r8),         parameter :: cleng          =  3._r8        !  Order of 'leng' when choice_leng = 'origin' [ no unit ]
  character(len=6), parameter :: choice_tkes    = 'ibprod'      ! 'ibprod' (include tkes in computing bprod), 'ebprod'(exclude)

  real(r8)                    :: lbulk_max      =  40.e3_r8     ! Maximum master length scale designed to address issues in the
                                                                ! upper atmosphere where vertical model resolution is coarse [ m ].
                                                                ! In order not to disturb turbulence characteristics in the lower
                                                                ! troposphere, this should be set at least larger than ~ a few km.  
  real(r8)                    :: leng_max(pver) =  40.e3_r8     ! Maximum length scale designed to address issues in the upper
                                                                ! atmosphere.

  ! Parameters for 'sedimentation-entrainment feedback' for liquid stratus
  ! If .false.,  no sedimentation entrainment feedback ( i.e., use default evhc )

  logical,          parameter :: id_sedfact     = .false.
  real(r8),         parameter :: ased           =  9._r8        !  Valid only when id_sedfact = .true.

  ! --------------------------------------------------------------------------------------------------- !
  ! Parameters governing entrainment efficiency A = a1l(i)*evhc, evhc = 1 + a2l * a3l * L * ql / jt2slv !
  ! Here, 'ql' is cloud-top LWC and 'jt2slv' is the jump in 'slv' across                                !
  ! the cloud-top entrainment zone ( across two grid layers to consider full mixture )                  !
  ! --------------------------------------------------------------------------------------------------- !

  real(r8),         parameter :: a1l            =   0.10_r8     ! Dry entrainment efficiency for TKE closure
                                                                ! a1l = 0.2*tunl*erat^-1.5,
                                                                ! where erat = <e>/wstar^2 for dry CBL =  0.3.

  real(r8),         parameter :: a1i            =   0.2_r8      ! Dry entrainment efficiency for wstar closure
  real(r8),         parameter :: ccrit          =   0.5_r8      ! Minimum allowable sqrt(tke)/wstar.
                                                                ! Used in solving cubic equation for 'ebrk'
  real(r8),         parameter :: wstar3factcrit =   0.5_r8      ! 1/wstar3factcrit is the maximally allowed enhancement of
                                                                ! 'wstar3' due to entrainment.

  real(r8)                    :: a2l                            ! Moist entrainment enhancement param (recommended range : 10~30 )
  real(r8),         parameter :: a3l            =   0.8_r8      ! Approximation to a complicated thermodynamic parameters

  real(r8),         parameter :: jbumin         =   .001_r8     ! Minimum buoyancy jump at an entrainment jump, [m/s2]
  real(r8),         parameter :: evhcmax        =   10._r8      ! Upper limit of evaporative enhancement factor

  real(r8),         parameter :: onet           =   1._r8/3._r8 ! 1/3 power in wind gradient expression [ no unit ]
  integer,          parameter :: ncvmax         =   pver        ! Max numbers of CLs (good to set to 'pver')
  real(r8),         parameter :: qmin           =   1.e-5_r8    ! Minimum grid-mean LWC counted as clouds [kg/kg]
  real(r8),         parameter :: ntzero         =   1.e-12_r8   ! Not zero (small positive number used in 's2')
  real(r8),         parameter :: b1             =   5.8_r8      ! TKE dissipation D = e^3/(b1*leng), e = b1*W.
  real(r8)                    :: b123                           ! b1**(2/3)
  real(r8),         parameter :: tunl           =   0.085_r8    ! Asympt leng = tunl*(turb lay depth)
  real(r8),         parameter :: alph1          =   0.5562_r8   ! alph1~alph5 : Galperin instability function parameters
  real(r8),         parameter :: alph2          =  -4.3640_r8   !               These coefficients are used to calculate 
  real(r8),         parameter :: alph3          = -34.6764_r8   !               'sh' and 'sm' from 'gh'.
  real(r8),         parameter :: alph4          =  -6.1272_r8   !
  real(r8),         parameter :: alph5          =   0.6986_r8   !
  real(r8),         parameter :: ricrit         =   0.19_r8     ! Critical Richardson number for turbulence.
                                                                ! Can be any value >= 0.19.
  real(r8),         parameter :: ae             =   1._r8       ! TKE transport efficiency [no unit]
  real(r8),         parameter :: rinc           =  -0.04_r8     ! Minimum W/<W> used for CL merging test 
  real(r8),         parameter :: wpertmin       =   1.e-6_r8    ! Minimum PBL eddy vertical velocity perturbation
  real(r8),         parameter :: wfac           =   1._r8       ! Ratio of 'wpert' to sqrt(tke) for CL.
  real(r8),         parameter :: tfac           =   1._r8       ! Ratio of 'tpert' to (w't')/wpert for CL.
                                                                ! Same ratio also used for q
  real(r8),         parameter :: fak            =   8.5_r8      ! Constant in surface temperature excess for stable STL.
                                                                ! [ no unit ]
  real(r8),         parameter :: rcapmin        =   0.1_r8      ! Minimum allowable e/<e> in a CL
  real(r8),         parameter :: rcapmax        =   2.0_r8      ! Maximum allowable e/<e> in a CL
  real(r8),         parameter :: tkemax         =  20._r8       ! TKE is capped at tkemax [m2/s2]
  real(r8),         parameter :: lambda         =   0.5_r8      ! Under-relaxation factor ( 0 < lambda =< 1 )

  logical,          parameter :: use_kvf        =  .false.      ! .true. (.false.) : initialize kvh/kvm =  kvf ( 0. )
  logical,          parameter :: use_dw_surf    =  .true.       ! Used in 'zisocl'. Default is 'true'
                                                                ! If 'true', surface interfacial energy does not contribute
                                                                ! to the CL mean stability functions after finishing merging.
                                                                ! For this case, 'dl2n2_surf' is only used for a merging test
                                                                ! based on 'l2n2'
                                                                ! If 'false',surface interfacial enery explicitly contribute to
                                                                ! CL mean stability functions after finishing merging.
                                                                ! For this case, 'dl2n2_surf' and 'dl2s2_surf' are directly used
                                                                ! for calculating surface interfacial layer energetics

  logical,          parameter :: set_qrlzero    =  .false.      ! .true. ( .false.) : turning-off ( on) radiative-turbulence
                                                                ! interaction by setting qrl = 0.

  ! ------------------------------------- !
  ! PBL Parameters not used in the UW PBL !
  ! ------------------------------------- !

  real(r8),         parameter :: pblmaxp        =  4.e4_r8      ! PBL max depth in pressure units. 
  real(r8),         parameter :: zkmin          =  0.01_r8      ! Minimum kneutral*f(ri). 
  real(r8),         parameter :: betam          = 15.0_r8       ! Constant in wind gradient expression.
  real(r8),         parameter :: betas          =  5.0_r8       ! Constant in surface layer gradient expression.
  real(r8),         parameter :: betah          = 15.0_r8       ! Constant in temperature gradient expression.
  real(r8),         parameter :: fakn           =  7.2_r8       ! Constant in turbulent prandtl number.
  real(r8),         parameter :: ricr           =  0.3_r8       ! Critical richardson number.
  real(r8),         parameter :: sffrac         =  0.1_r8       ! Surface layer fraction of boundary layer
  real(r8),         parameter :: binm           =  betam*sffrac ! betam * sffrac
  real(r8),         parameter :: binh           =  betah*sffrac ! betah * sffrac

  ! ------------------------------------------------------- !
  ! PBL constants set using values from other parts of code !
  ! ------------------------------------------------------- !

  real(r8)                    :: cpair                          ! Specific heat of dry air
  real(r8)                    :: rair                           ! Gas const for dry air
  real(r8)                    :: zvir                           ! rh2o/rair - 1
  real(r8)                    :: latvap                         ! Latent heat of vaporization
  real(r8)                    :: latice                         ! Latent heat of fusion
  real(r8)                    :: latsub                         ! Latent heat of sublimation
  real(r8)                    :: g                              ! Gravitational acceleration
  real(r8)                    :: vk                             ! Von Karman's constant
  real(r8)                    :: ccon                           ! fak * sffrac * vk

  integer                     :: ntop_turb                      ! Top interface level to which turbulent vertical diffusion
                                                                ! is applied ( = 1 )
  integer                     :: nbot_turb                      ! Bottom interface level to which turbulent vertical diff
                                                                ! is applied ( = pver )

  real(r8), allocatable, target :: ml2(:)                       ! Mixing lengths squared. Not used in the UW PBL.
                                                                ! Used for computing free air diffusivity.
  logical                     :: use_native_surface_stress_diag_impl = .false.
  logical                     :: surface_stress_diag_impl_selected = .false.
  logical                     :: use_native_austausch_atm_impl = .false.
  logical                     :: austausch_atm_impl_selected = .false.
  logical                     :: use_native_kv_init_impl = .false.
  logical                     :: kv_init_impl_selected = .false.
  logical                     :: use_native_error_pbl_impl = .false.
  logical                     :: error_pbl_impl_selected = .false.
  logical                     :: use_native_kv_relax_impl = .false.
  logical                     :: kv_relax_impl_selected = .false.
  logical                     :: use_native_init_fields_impl = .false.
  logical                     :: init_fields_impl_selected = .false.
  logical                     :: use_native_rebuild_thermo_impl = .false.
  logical                     :: rebuild_thermo_impl_selected = .false.
  logical                     :: use_native_trbintd_midpoint_impl = .false.
  logical                     :: trbintd_midpoint_impl_selected = .false.
  logical                     :: use_native_trbintd_slopes_impl = .false.
  logical                     :: trbintd_slopes_impl_selected = .false.
  logical                     :: use_native_trbintd_sfdiag_interface_impl = .false.
  logical                     :: trbintd_sfdiag_interface_impl_selected = .false.
  logical                     :: use_native_trbintd_core_impl = .false.
  logical                     :: trbintd_core_impl_selected = .false.
  logical                     :: use_native_zero_nonlocal_impl = .false.
  logical                     :: zero_nonlocal_impl_selected = .false.
  logical                     :: use_native_restore_fields_impl = .false.
  logical                     :: restore_fields_impl_selected = .false.
  logical                     :: use_native_wstar_pbl_impl = .false.
  logical                     :: wstar_pbl_impl_selected = .false.
  logical                     :: use_native_exacol_impl = .false.
  logical                     :: exacol_impl_selected = .false.
  logical                     :: use_native_compute_radf_impl = .false.
  logical                     :: compute_radf_impl_selected = .false.
  logical                     :: use_native_caleddy_init_impl = .false.
  logical                     :: caleddy_init_impl_selected = .false.
  logical                     :: use_native_caleddy_diaginit_impl = .false.
  logical                     :: caleddy_diaginit_impl_selected = .false.
  logical                     :: use_native_caleddy_regime_diag_impl = .false.
  logical                     :: caleddy_regime_diag_impl_selected = .false.
  logical                     :: use_native_caleddy_stable_config_impl = .false.
  logical                     :: caleddy_stable_config_impl_selected = .false.
  logical                     :: use_native_caleddy_surface_tke_impl = .false.
  logical                     :: caleddy_surface_tke_impl_selected = .false.
  logical                     :: use_native_zisocl_surface_energy_impl = .false.
  logical                     :: zisocl_surface_energy_impl_selected = .false.
  logical                     :: use_native_zisocl_surface_state_impl = .false.
  logical                     :: zisocl_surface_state_impl_selected = .false.
  logical                     :: use_native_zisocl_surface_extend_impl = .false.
  logical                     :: zisocl_surface_extend_impl_selected = .false.
  logical                     :: use_native_zisocl_sbcl_state_impl = .false.
  logical                     :: zisocl_sbcl_state_impl_selected = .false.
  logical                     :: use_native_zisocl_initial_state_impl = .false.
  logical                     :: zisocl_initial_state_impl_selected = .false.
  logical                     :: use_native_zisocl_extended_state_impl = .false.
  logical                     :: zisocl_extended_state_impl_selected = .false.
  logical                     :: use_native_zisocl_non_sbcl_state_impl = .false.
  logical                     :: zisocl_non_sbcl_state_impl_selected = .false.
  logical                     :: use_native_zisocl_upward_state_impl = .false.
  logical                     :: zisocl_upward_state_impl_selected = .false.
  logical                     :: use_native_zisocl_downward_state_impl = .false.
  logical                     :: zisocl_downward_state_impl_selected = .false.
  logical                     :: use_native_zisocl_stability_impl = .false.
  logical                     :: zisocl_stability_impl_selected = .false.
  logical                     :: use_native_zisocl_layer_energy_impl = .false.
  logical                     :: zisocl_layer_energy_impl_selected = .false.
  logical                     :: use_native_zisocl_interface_energy_impl = .false.
  logical                     :: zisocl_interface_energy_impl_selected = .false.
  logical                     :: use_native_caleddy_clprep_impl = .false.
  logical                     :: caleddy_clprep_impl_selected = .false.
  logical                     :: use_native_caleddy_closure_impl = .false.
  logical                     :: caleddy_closure_impl_selected = .false.
  logical                     :: use_native_caleddy_srcl_impl = .false.
  logical                     :: caleddy_srcl_impl_selected = .false.
  logical                     :: use_native_caleddy_stl_impl = .false.
  logical                     :: caleddy_stl_impl_selected = .false.
  logical                     :: use_native_caleddy_diag_impl = .false.
  logical                     :: caleddy_diag_impl_selected = .false.

  CONTAINS

  !============================================================================ !
  !                                                                             !
  !============================================================================ !
  
  subroutine init_eddy_diff( kind, pver, gravx, cpairx, rairx, zvirx, & 
                             latvapx, laticex, ntop_eddy, nbot_eddy, vkx, &
                             eddy_lbulk_max, eddy_leng_max, eddy_max_bot_pressure, &
                             eddy_moist_entrain_a2l)
    !---------------------------------------------------------------- ! 
    ! Purpose:                                                        !
    ! Initialize time independent constants/variables of PBL package. !
    !---------------------------------------------------------------- !
    use diffusion_solver, only: new_fieldlist_vdiff, vdiff_select
    use cam_history,      only: outfld, addfld, phys_decomp
    use ref_pres,         only: pref_mid
    
    ! --------- !
    ! Arguments !
    ! --------- !
    integer,  intent(in) :: kind       ! Kind of reals being passed in
    integer,  intent(in) :: pver       ! Number of vertical layers
    integer,  intent(in) :: ntop_eddy  ! Top interface level to which eddy vertical diffusivity is applied ( = 1 )
    integer,  intent(in) :: nbot_eddy  ! Bottom interface level to which eddy vertical diffusivity is applied ( = pver )
    real(r8), intent(in) :: gravx      ! Acceleration of gravity
    real(r8), intent(in) :: cpairx     ! Specific heat of dry air
    real(r8), intent(in) :: rairx      ! Gas constant for dry air
    real(r8), intent(in) :: zvirx      ! rh2o/rair - 1
    real(r8), intent(in) :: latvapx    ! Latent heat of vaporization
    real(r8), intent(in) :: laticex    ! Latent heat of fusion
    real(r8), intent(in) :: vkx        ! Von Karman's constant
    real(r8), intent(in) :: eddy_lbulk_max ! Maximum master length scale
    real(r8), intent(in) :: eddy_leng_max  ! Maximum dissipation length scale
    real(r8), intent(in) :: eddy_max_bot_pressure  ! Bottom pressure level (hPa) at which namelist leng_max and lbulk_max
                                                   ! are applied
    real(r8), intent(in) :: eddy_moist_entrain_a2l ! Moist entrainment enhancement param

    integer              :: k          ! Vertical loop index

    if( kind .ne. r8 ) then
        write(iulog,*) 'wrong KIND of reals passed to init_diffusvity -- exiting.'
        call endrun('init_eddy_diff: wrong KIND of reals passed to init_diffusvity')
    endif

    ! --------------- !
    ! Basic constants !
    ! --------------- !

    cpair     = cpairx
    rair      = rairx
    g         = gravx
    zvir      = zvirx
    latvap    = latvapx
    latice    = laticex
    latsub    = latvap + latice
    vk        = vkx
    ccon      = fak*sffrac*vk
    ntop_turb = ntop_eddy
    nbot_turb = nbot_eddy
    b123      = b1**(2._r8/3._r8)
    a2l       = eddy_moist_entrain_a2l
    
    lbulk_max = eddy_lbulk_max
    do k = 1,pver
      if ( pref_mid(k) .le. eddy_max_bot_pressure*1.D2 ) leng_max(k)  = eddy_leng_max
    end do

    if (masterproc) then
       write(iulog,*)'init_eddy_diff: eddy_leng_max=',eddy_leng_max,' lbulk_max=',lbulk_max
       do k = 1,pver
          write(iulog,*)'init_eddy_diff:',k,pref_mid(k),'leng_max=',leng_max(k)
       end do
    end if

    ! Set the square of the mixing lengths. Only for CAM3 HB PBL scheme.
    ! Not used for UW moist PBL. Used for free air eddy diffusivity.

    allocate(ml2(pver+1))
    ml2(1:ntop_turb) = 0._r8
    do k = ntop_turb + 1, nbot_turb
       ml2(k) = 30.0_r8**2
    end do
    ml2(nbot_turb+1:pver+1) = 0._r8
    
    ! Get fieldlists to pass to diffusion solver.
    fieldlist_wet   = new_fieldlist_vdiff(1)
    fieldlist_molec = new_fieldlist_vdiff(1)

    ! Select the fields which will be diffused 

    if(vdiff_select(fieldlist_wet,'s').ne.'')   call endrun( vdiff_select(fieldlist_wet,'s') )
    if(vdiff_select(fieldlist_wet,'q',1).ne.'') call endrun( vdiff_select(fieldlist_wet,'q',1) )
    if(vdiff_select(fieldlist_wet,'u').ne.'')   call endrun( vdiff_select(fieldlist_wet,'u') )
    if(vdiff_select(fieldlist_wet,'v').ne.'')   call endrun( vdiff_select(fieldlist_wet,'v') )

    ! ------------------------------------------------------------------- !
    ! Writing outputs for detailed analysis of UW moist turbulence scheme !
    ! ------------------------------------------------------------------- !

    call addfld('UW_errorPBL',      'm2/s',    1,      'A',  'Error function of UW PBL',                              phys_decomp )
    call addfld('UW_n2',            's-2',     pver,   'A',  'Buoyancy Frequency, LI',                                phys_decomp )
    call addfld('UW_s2',            's-2',     pver,   'A',  'Shear Frequency, LI',                                   phys_decomp )
    call addfld('UW_ri',            'no',      pver,   'A',  'Interface Richardson Number, I',                        phys_decomp )
    call addfld('UW_sfuh',          'no',      pver,   'A',  'Upper-Half Saturation Fraction, L',                     phys_decomp )
    call addfld('UW_sflh',          'no',      pver,   'A',  'Lower-Half Saturation Fraction, L',                     phys_decomp )
    call addfld('UW_sfi',           'no',      pver+1, 'A',  'Interface Saturation Fraction, I',                      phys_decomp )
    call addfld('UW_cldn',          'no',      pver,   'A',  'Cloud Fraction, L',                                     phys_decomp )
    call addfld('UW_qrl',           'g*W/m2',  pver,   'A',  'LW cooling rate, L',                                    phys_decomp )
    call addfld('UW_ql',            'kg/kg',   pver,   'A',  'ql(LWC), L',                                            phys_decomp )
    call addfld('UW_chu',           'g*kg/J',  pver+1, 'A',  'Buoyancy Coefficient, chu, I',                          phys_decomp )
    call addfld('UW_chs',           'g*kg/J',  pver+1, 'A',  'Buoyancy Coefficient, chs, I',                          phys_decomp )
    call addfld('UW_cmu',           'g/kg/kg', pver+1, 'A',  'Buoyancy Coefficient, cmu, I',                          phys_decomp )
    call addfld('UW_cms',           'g/kg/kg', pver+1, 'A',  'Buoyancy Coefficient, cms, I',                          phys_decomp )
    call addfld('UW_tke',           'm2/s2',   pver+1, 'A',  'TKE, I',                                                phys_decomp )
    call addfld('UW_wcap',          'm2/s2',   pver+1, 'A',  'Wcap, I',                                               phys_decomp )
    call addfld('UW_bprod',         'm2/s3',   pver+1, 'A',  'Buoyancy production, I',                                phys_decomp )
    call addfld('UW_sprod',         'm2/s3',   pver+1, 'A',  'Shear production, I',                                   phys_decomp )
    call addfld('UW_kvh',           'm2/s',    pver+1, 'A',  'Eddy diffusivity of heat, I',                           phys_decomp )
    call addfld('UW_kvm',           'm2/s',    pver+1, 'A',  'Eddy diffusivity of uv, I',                             phys_decomp )
    call addfld('UW_pblh',          'm',       1,      'A',  'PBLH, 1',                                               phys_decomp )
    call addfld('UW_pblhp',         'Pa',      1,      'A',  'PBLH pressure, 1',                                      phys_decomp )
    call addfld('UW_tpert',         'K',       1,      'A',  'Convective T excess, 1',                                phys_decomp )
    call addfld('UW_qpert',         'kg/kg',   1,      'A',  'Convective qt excess, I',                               phys_decomp )
    call addfld('UW_wpert',         'm/s',     1,      'A',  'Convective W excess, I',                                phys_decomp )
    call addfld('UW_ustar',         'm/s',     1,      'A',  'Surface Frictional Velocity, 1',                        phys_decomp )
    call addfld('UW_tkes',          'm2/s2',   1,      'A',  'Surface TKE, 1',                                        phys_decomp )
    call addfld('UW_minpblh',       'm',       1,      'A',  'Minimum PBLH, 1',                                       phys_decomp )
    call addfld('UW_turbtype',      'no',      pver+1, 'A',  'Interface Turbulence Type, I',                          phys_decomp )
    call addfld('UW_kbase_o',       'no',      ncvmax, 'A',  'Initial CL Base Exterbal Interface Index, CL',          phys_decomp )
    call addfld('UW_ktop_o',        'no',      ncvmax, 'A',  'Initial Top Exterbal Interface Index, CL',              phys_decomp )
    call addfld('UW_ncvfin_o',      '#',       1,      'A',  'Initial Total Number of CL regimes, CL',                phys_decomp )
    call addfld('UW_kbase_mg',      'no',      ncvmax, 'A',  'kbase after merging, CL',                               phys_decomp )
    call addfld('UW_ktop_mg',       'no',      ncvmax, 'A',  'ktop after merging, CL',                                phys_decomp )
    call addfld('UW_ncvfin_mg',     '#',       1,      'A',  'ncvfin after merging, CL',                              phys_decomp )
    call addfld('UW_kbase_f',       'no',      ncvmax, 'A',  'Final kbase with SRCL, CL',                             phys_decomp )
    call addfld('UW_ktop_f',        'no',      ncvmax, 'A',  'Final ktop with SRCL, CL',                              phys_decomp )
    call addfld('UW_ncvfin_f',      '#',       1,      'A',  'Final ncvfin with SRCL, CL',                            phys_decomp )
    call addfld('UW_wet',           'm/s',     ncvmax, 'A',  'Entrainment rate at CL top, CL',                        phys_decomp )
    call addfld('UW_web',           'm/s',     ncvmax, 'A',  'Entrainment rate at CL base, CL',                       phys_decomp )
    call addfld('UW_jtbu',          'm/s2',    ncvmax, 'A',  'Buoyancy jump across CL top, CL',                       phys_decomp )
    call addfld('UW_jbbu',          'm/s2',    ncvmax, 'A',  'Buoyancy jump across CL base, CL',                      phys_decomp )
    call addfld('UW_evhc',          'no',      ncvmax, 'A',  'Evaporative enhancement factor, CL',                    phys_decomp )
    call addfld('UW_jt2slv',        'J/kg',    ncvmax, 'A',  'slv jump for evhc, CL',                                 phys_decomp )
    call addfld('UW_n2ht',          's-2',     ncvmax, 'A',  'n2 at just below CL top interface, CL',                 phys_decomp )
    call addfld('UW_n2hb',          's-2',     ncvmax, 'A',  'n2 at just above CL base interface',                    phys_decomp )
    call addfld('UW_lwp',           'kg/m2',   ncvmax, 'A',  'LWP in the CL top layer, CL',                           phys_decomp )
    call addfld('UW_optdepth',      'no',      ncvmax, 'A',  'Optical depth of the CL top layer, CL',                 phys_decomp )
    call addfld('UW_radfrac',       'no',      ncvmax, 'A',  'Fraction of radiative cooling confined in the CL top',  phys_decomp )
    call addfld('UW_radf',          'm2/s3',   ncvmax, 'A',  'Buoyancy production at the CL top by radf, I',          phys_decomp )
    call addfld('UW_wstar',         'm/s',     ncvmax, 'A',  'Convective velocity, Wstar, CL',                        phys_decomp )
    call addfld('UW_wstar3fact',    'no',      ncvmax, 'A',  'Enhancement of wstar3 due to entrainment, CL',          phys_decomp )
    call addfld('UW_ebrk',          'm2/s2',   ncvmax, 'A',  'CL-averaged TKE, CL',                                   phys_decomp )
    call addfld('UW_wbrk',          'm2/s2',   ncvmax, 'A',  'CL-averaged W, CL',                                     phys_decomp )
    call addfld('UW_lbrk',          'm',       ncvmax, 'A',  'CL internal thickness, CL',                             phys_decomp )
    call addfld('UW_ricl',          'no',      ncvmax, 'A',  'CL-averaged Ri, CL',                                    phys_decomp )
    call addfld('UW_ghcl',          'no',      ncvmax, 'A',  'CL-averaged gh, CL',                                    phys_decomp )
    call addfld('UW_shcl',          'no',      ncvmax, 'A',  'CL-averaged sh, CL',                                    phys_decomp )
    call addfld('UW_smcl',          'no',      ncvmax, 'A',  'CL-averaged sm, CL',                                    phys_decomp )
    call addfld('UW_gh',            'no',      pver+1, 'A',  'gh at all interfaces, I',                               phys_decomp )
    call addfld('UW_sh',            'no',      pver+1, 'A',  'sh at all interfaces, I',                               phys_decomp )
    call addfld('UW_sm',            'no',      pver+1, 'A',  'sm at all interfaces, I',                               phys_decomp )
    call addfld('UW_ria',           'no',      pver+1, 'A',  'ri at all interfaces, I',                               phys_decomp )
    call addfld('UW_leng',          'm/s',     pver+1, 'A',  'Turbulence length scale, I',                            phys_decomp )
    ! For sedimentation-entrainment feedback analysis
    call addfld('UW_wsed',          'm/s',     ncvmax, 'A',  'Sedimentation velocity at CL top, CL',                  phys_decomp )

  return

  end subroutine init_eddy_diff

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !
  
  subroutine compute_eddy_diff( lchnk  ,                                                            &
                                pcols  , pver   , ncol     , t       , qv       , ztodt   ,         &
                                ql     , qi     , s        , pdel, rpdel   , cldn     , qrl     , wsedl , &
                                z      , zi     , pmid     , pi      , u        , v       ,         &
                                taux   , tauy   , shflx    , qflx    , wstarent , nturb   , rrho  , &
                                ustar  , pblh   , kvm_in   , kvh_in  , kvm_out  , kvh_out , kvq   , & 
                                cgh    , cgs    , tpert    , qpert   , wpert    , tke     , bprod , &
                                sprod  , sfi    , kvinit   ,                                        &
                                tauresx, tauresy, ksrftms  ,                                        &
                                ipbl   , kpblh  , wstarPBL , tkes    , went     ,turbtype, sm_aw )
       
    !-------------------------------------------------------------------- ! 
    ! Purpose: Interface to compute eddy diffusivities.                   !
    !          Eddy diffusivities are calculated in a fully implicit way  !
    !          through iteration process.                                 !   
    ! Author:  Sungsu Park. August. 2006.                                 !
    !                       May.    2008.                                 !
    !-------------------------------------------------------------------- !

    use diffusion_solver, only: compute_vdiff
    use cam_history,      only: outfld, addfld, phys_decomp
  ! use physics_types,    only: physics_state
    use phys_debug_util,  only: phys_debug_col
    use time_manager,     only: is_first_step, get_nstep
    use physconst,        only: cpairv, rairv, rair !Needed for call to compute_vdiff 
    use phys_control,     only: waccmx_is
    use pbl_utils,        only: calc_ustar
    use error_messages,   only: handle_errmsg

    implicit none

  ! type(physics_state)     :: state                     ! Physics state variables

    ! --------------- !
    ! Input Variables !
    ! --------------- ! 

    integer,  intent(in)    :: lchnk   
    integer,  intent(in)    :: pcols                     ! Number of atmospheric columns [ # ]
    integer,  intent(in)    :: pver                      ! Number of atmospheric layers  [ # ]
    integer,  intent(in)    :: ncol                      ! Number of atmospheric columns [ # ]
    integer,  intent(in)    :: nturb                     ! Number of iteration steps for calculating eddy diffusivity [ # ]
    logical,  intent(in)    :: wstarent                  ! .true. means use the 'wstar' entrainment closure. 
    logical,  intent(in)    :: kvinit                    ! 'true' means time step = 1 : used for initializing kvh, kvm
                                                         ! (uses kvf or zero)
    real(r8), intent(in)    :: ztodt                     ! Physics integration time step 2 delta-t [ s ]
    real(r8), intent(in)    :: t(pcols,pver)             ! Temperature [K]
    real(r8), intent(in)    :: qv(pcols,pver)            ! Water vapor  specific humidity [ kg/kg ]
    real(r8), intent(in)    :: ql(pcols,pver)            ! Liquid water specific humidity [ kg/kg ]
    real(r8), intent(in)    :: qi(pcols,pver)            ! Ice specific humidity [ kg/kg ]
    real(r8), intent(in)    :: s(pcols,pver)             ! Dry static energy [ J/kg ]
    real(r8), intent(in)    :: pdel(pcols,pver)          ! thickness of the layer [ Pa ]
    real(r8), intent(in)    :: rpdel(pcols,pver)         ! 1./pdel where 'pdel' is thickness of the layer [ 1/Pa ]
    real(r8), intent(in)    :: cldn(pcols,pver)          ! Stratiform cloud fraction [ fraction ]
    real(r8), intent(in)    :: qrl(pcols,pver)           ! LW cooling rate
    real(r8), intent(in)    :: wsedl(pcols,pver)         ! Sedimentation velocity of liquid stratus cloud droplet [ m/s ]
    real(r8), intent(in)    :: z(pcols,pver)             ! Layer mid-point height above surface [ m ]
    real(r8), intent(in)    :: zi(pcols,pver+1)          ! Interface height above surface [ m ]
    real(r8), intent(in)    :: pmid(pcols,pver)          ! Layer mid-point pressure [ Pa ]
    real(r8), intent(in)    :: pi(pcols,pver+1)          ! Interface pressure [ Pa ]
    real(r8), intent(in)    :: u(pcols,pver)             ! Zonal velocity [ m/s ]
    real(r8), intent(in)    :: v(pcols,pver)             ! Meridional velocity [ m/s ]
    real(r8), intent(in)    :: taux(pcols)               ! Zonal wind stress at surface [ N/m2 ]
    real(r8), intent(in)    :: tauy(pcols)               ! Meridional wind stress at surface [ N/m2 ]
    real(r8), intent(in)    :: shflx(pcols)              ! Sensible heat flux at surface [ unit ? ]
    real(r8), intent(in)    :: qflx(pcols)               ! Water vapor flux at surface [ unit ? ]
    real(r8), intent(in)    :: kvm_in(pcols,pver+1)      ! kvm saved from last timestep [ m2/s ]
    real(r8), intent(in)    :: kvh_in(pcols,pver+1)      ! kvh saved from last timestep [ m2/s ]
    real(r8), intent(in)    :: ksrftms(pcols)            ! Surface drag coefficient of turbulent mountain stress [ unit ? ]

    ! ---------------- !
    ! Output Variables !
    ! ---------------- ! 

    real(r8), intent(out)   :: kvm_out(pcols,pver+1)     ! Eddy diffusivity for momentum [ m2/s ]
    real(r8), intent(out)   :: kvh_out(pcols,pver+1)     ! Eddy diffusivity for heat [ m2/s ]
    real(r8), intent(out)   :: kvq(pcols,pver+1)         ! Eddy diffusivity for constituents, moisture and tracers [ m2/s ]
                                                         ! (note not having '_out')
    real(r8), intent(out)   :: rrho(pcols)               ! Reciprocal of density at the lowest layer
    real(r8), intent(out)   :: ustar(pcols)              ! Surface friction velocity [ m/s ]
    real(r8), intent(out)   :: pblh(pcols)               ! PBL top height [ m ]
    real(r8), intent(out)   :: cgh(pcols,pver+1)         ! Counter-gradient term for heat [ J/kg/m ]
    real(r8), intent(out)   :: cgs(pcols,pver+1)         ! Counter-gradient star [ cg/flux ]
    real(r8), intent(out)   :: tpert(pcols)              ! Convective temperature excess [ K ]
    real(r8), intent(out)   :: qpert(pcols)              ! Convective humidity excess [ kg/kg ]
    real(r8), intent(out)   :: wpert(pcols)              ! Turbulent velocity excess [ m/s ]
    real(r8), intent(out)   :: tke(pcols,pver+1)         ! Turbulent kinetic energy [ m2/s2 ]
    real(r8), intent(out)   :: bprod(pcols,pver+1)       ! Buoyancy production [ m2/s3 ] 
    real(r8), intent(out)   :: sprod(pcols,pver+1)       ! Shear production [ m2/s3 ] 
    real(r8), intent(out)   :: sfi(pcols,pver+1)         ! Interfacial layer saturation fraction [ fraction ]
    integer(i4), intent(out)   :: turbtype(pcols,pver+1)    ! Turbulence type identifier at all interfaces [ no unit ]
    real(r8), intent(out)   :: sm_aw(pcols,pver+1)       ! Normalized Galperin instability function for momentum [ no unit ]
                                                         ! This is 1 when neutral condition (Ri=0),
                                                         ! 4.964 for maximum unstable case, and 0 when Ri > Ricrit=0.19. 
    integer(i4), intent(out) :: ipbl(pcols)              ! If 1, PBL is CL, while if 0, PBL is STL.
    integer(i4), intent(out) :: kpblh(pcols)             ! Layer index containing PBL top within or at the base interface
    real(r8), intent(out)   :: wstarPBL(pcols)           ! Convective velocity within PBL [ m/s ]
    real(r8), intent(out)   :: tkes(pcols)               ! TKE at surface interface [ m2/s2 ]
    real(r8), intent(out)   :: went(pcols)               ! Entrainment rate at the PBL top interface [ m/s ]

    ! ---------------------- !
    ! Input-Output Variables !
    ! ---------------------- ! 

    real(r8), intent(inout) :: tauresx(pcols)            ! Residual stress to be added in vdiff to correct for turb
    real(r8), intent(inout) :: tauresy(pcols)            ! Stress mismatch between sfc and atm accumulated in prior timesteps

    ! --------------- !
    ! Local Variables !
    ! --------------- !

    integer                    icol
    integer                    i, k, iturb, status

    character(128)          :: errstring                 ! Error status for compute_vdiff

    real(r8)                :: kvf(pcols,pver+1)         ! Free atmospheric eddy diffusivity [ m2/s ]
    real(r8)                :: kvm(pcols,pver+1)         ! Eddy diffusivity for momentum [ m2/s ]
    real(r8)                :: kvh(pcols,pver+1)         ! Eddy diffusivity for heat [ m2/s ]
    real(r8)                :: kvm_preo(pcols,pver+1)    ! Eddy diffusivity for momentum [ m2/s ]
    real(r8)                :: kvh_preo(pcols,pver+1)    ! Eddy diffusivity for heat [ m2/s ]
    real(r8)                :: kvm_pre(pcols,pver+1)     ! Eddy diffusivity for momentum [ m2/s ]
    real(r8)                :: kvh_pre(pcols,pver+1)     ! Eddy diffusivity for heat [ m2/s ]
    real(r8)                :: errorPBL(pcols)           ! Error function showing whether PBL produced convergent solution or not.
                                                         ! [ unit ? ]
    real(r8)                :: s2(pcols,pver)            ! Shear squared, defined at interfaces except surface [ s-2 ]
    real(r8)                :: n2(pcols,pver)            ! Buoyancy frequency, defined at interfaces except surface [ s-2 ]
    real(r8)                :: ri(pcols,pver)            ! Richardson number, 'n2/s2', defined at interfaces except surface [ s-2 ]
    real(r8)                :: pblhp(pcols)              ! PBL top pressure [ Pa ]
    real(r8)                :: minpblh(pcols)            ! Minimum PBL height based on surface stress

    real(r8)                :: qt(pcols,pver)            ! Total specific humidity [ kg/kg ]
    real(r8)                :: sfuh(pcols,pver)          ! Saturation fraction in upper half-layer [ fraction ]
    real(r8)                :: sflh(pcols,pver)          ! Saturation fraction in lower half-layer [ fraction ]
    real(r8)                :: sl(pcols,pver)            ! Liquid water static energy [ J/kg ]
    real(r8)                :: slv(pcols,pver)           ! Liquid water virtual static energy [ J/kg ]
    real(r8)                :: slslope(pcols,pver)       ! Slope of 'sl' in each layer
    real(r8)                :: qtslope(pcols,pver)       ! Slope of 'qt' in each layer
    real(r8)                :: qvfd(pcols,pver)          ! Specific humidity for diffusion [ kg/kg ]
    real(r8)                :: tfd(pcols,pver)           ! Temperature for diffusion [ K ]
    real(r8)                :: slfd(pcols,pver)          ! Liquid static energy [ J/kg ]
    real(r8)                :: qtfd(pcols,pver)          ! Total specific humidity [ kg/kg ] 
    real(r8)                :: qlfd(pcols,pver)          ! Liquid water specific humidity for diffusion [ kg/kg ]
    real(r8)                :: ufd(pcols,pver)           ! U-wind for diffusion [ m/s ]
    real(r8)                :: vfd(pcols,pver)           ! V-wind for diffusion [ m/s ]

    ! Buoyancy coefficients : w'b' = ch * w'sl' + cm * w'qt'

    real(r8)                :: chu(pcols,pver+1)         ! Heat buoyancy coef for dry states, defined at each interface, finally.
    real(r8)                :: chs(pcols,pver+1)         ! Heat buoyancy coef for sat states, defined at each interface, finally. 
    real(r8)                :: cmu(pcols,pver+1)         ! Moisture buoyancy coef for dry states,
                                                         ! defined at each interface, finally.
    real(r8)                :: cms(pcols,pver+1)         ! Moisture buoyancy coef for sat states,
                                                         ! defined at each interface, finally. 

    real(r8)                :: jnk1d(pcols)
    real(r8)                :: jnk2d(pcols,pver+1)  
    real(r8)                :: zero(pcols)
    real(r8)                :: zero2d(pcols,pver+1)
    real(r8)                :: es                     ! Saturation vapor pressure
    real(r8)                :: qs                     ! Saturation specific humidity
    real(r8)                :: ep2, templ, temps

    ! ------------------------------- !
    ! Variables for diagnostic output !
    ! ------------------------------- !

    real(r8)                :: kbase_o(pcols,ncvmax)     ! Original external base interface index of CL from 'exacol'
    real(r8)                :: ktop_o(pcols,ncvmax)      ! Original external top  interface index of CL from 'exacol'
    real(r8)                :: ncvfin_o(pcols)           ! Original number of CLs from 'exacol'
    real(r8)                :: kbase_mg(pcols,ncvmax)    ! 'kbase' after extending-merging from 'zisocl'
    real(r8)                :: ktop_mg(pcols,ncvmax)     ! 'ktop' after extending-merging from 'zisocl'
    real(r8)                :: ncvfin_mg(pcols)          ! 'ncvfin' after extending-merging from 'zisocl'
    real(r8)                :: kbase_f(pcols,ncvmax)     ! Final 'kbase' after extending-merging & including SRCL
    real(r8)                :: ktop_f(pcols,ncvmax)      ! Final 'ktop' after extending-merging & including SRCL
    real(r8)                :: ncvfin_f(pcols)           ! Final 'ncvfin' after extending-merging & including SRCL
    real(r8)                :: wet(pcols,ncvmax)         ! Entrainment rate at the CL top  [ m/s ] 
    real(r8)                :: web(pcols,ncvmax)         ! Entrainment rate at the CL base [ m/s ].
                                                         ! Set to zero if CL is based at surface.
    real(r8)                :: jtbu(pcols,ncvmax)        ! Buoyancy jump across the CL top  [ m/s2 ]  
    real(r8)                :: jbbu(pcols,ncvmax)        ! Buoyancy jump across the CL base [ m/s2 ]  
    real(r8)                :: evhc(pcols,ncvmax)        ! Evaporative enhancement factor at the CL top
    real(r8)                :: jt2slv(pcols,ncvmax)      ! Jump of slv ( across two layers ) at CL top used only for evhc [ J/kg ]
    real(r8)                :: n2ht(pcols,ncvmax)        ! n2 defined at the CL top  interface but using
                                                         ! sfuh(kt)   instead of sfi(kt) [ s-2 ] 
    real(r8)                :: n2hb(pcols,ncvmax)        ! n2 defined at the CL base interface but using
                                                         ! sflh(kb-1) instead of sfi(kb) [ s-2 ]
    real(r8)                :: lwp(pcols,ncvmax)         ! LWP in the CL top layer [ kg/m2 ]
    real(r8)                :: opt_depth(pcols,ncvmax)   ! Optical depth of the CL top layer
    real(r8)                :: radinvfrac(pcols,ncvmax)  ! Fraction of radiative cooling confined in the top portion of CL top layer
    real(r8)                :: radf(pcols,ncvmax)        ! Buoyancy production at the CL top due to LW radiative cooling [ m2/s3 ]
    real(r8)                :: wstar(pcols,ncvmax)       ! Convective velocity in each CL [ m/s ]
    real(r8)                :: wstar3fact(pcols,ncvmax)  ! Enhancement of 'wstar3' due to entrainment (inverse) [ no unit ]
    real(r8)                :: ebrk(pcols,ncvmax)        ! Net mean TKE of CL including entrainment effect [ m2/s2 ]
    real(r8)                :: wbrk(pcols,ncvmax)        ! Net mean normalized TKE (W) of CL,
                                                         ! 'ebrk/b1' including entrainment effect [ m2/s2 ]
    real(r8)                :: lbrk(pcols,ncvmax)        ! Energetic internal thickness of CL [m]
    real(r8)                :: ricl(pcols,ncvmax)        ! CL internal mean Richardson number
    real(r8)                :: ghcl(pcols,ncvmax)        ! Half of normalized buoyancy production of CL
    real(r8)                :: shcl(pcols,ncvmax)        ! Galperin instability function of heat-moisture of CL
    real(r8)                :: smcl(pcols,ncvmax)        ! Galperin instability function of mementum of CL
    real(r8)                :: ghi(pcols,pver+1)         ! Half of normalized buoyancy production at all interfaces
    real(r8)                :: shi(pcols,pver+1)         ! Galperin instability function of heat-moisture at all interfaces
    real(r8)                :: smi(pcols,pver+1)         ! Galperin instability function of heat-moisture at all interfaces
    real(r8)                :: rii(pcols,pver+1)         ! Interfacial Richardson number defined at all interfaces
    real(r8)                :: lengi(pcols,pver+1)       ! Turbulence length scale at all interfaces [ m ]
    real(r8)                :: wcap(pcols,pver+1)        ! Normalized TKE at all interfaces [ m2/s2 ]
    real(r8)                :: rairi(pcols,pver+1)       ! interface gas constant needed for compute_vdiff
    ! For sedimentation-entrainment feedback
    real(r8)                :: wsed(pcols,ncvmax)        ! Sedimentation velocity at the top of each CL [ m/s ]

    ! ---------- !
    ! Initialize !
    ! ---------- !

    call eddy_diff_init_fields(ncol, pcols, pver, u, v, t, qv, ql, zero, zero2d, ufd, vfd, tfd, qvfd, qlfd)

    ! ----------------------- !
    ! Main Computation Begins ! 
    ! ----------------------- !

    ufd(:ncol,:)  = u(:ncol,:)
    vfd(:ncol,:)  = v(:ncol,:)
    tfd(:ncol,:)  = t(:ncol,:)
    qvfd(:ncol,:) = qv(:ncol,:)
    qlfd(:ncol,:) = ql(:ncol,:)
    
    do iturb = 1, nturb

     ! Total stress includes 'tms'.
     ! Here, in computing 'tms', we can use either iteratively changed 'ufd,vfd' or the
     ! initially given 'u,v' to the PBL scheme. Note that normal stress, 'taux, tauy'
     ! are not changed by iteration. In order to treat 'tms' in a fully implicit way,
     ! I am using updated wind, here.

     ! Compute ustar
       call eddy_diff_surface_stress_diag(ncol, pcols, pver, tfd, pmid, taux, tauy, ksrftms, ufd, vfd, rrho, ustar, &
            minpblh)

     ! Calculate (qt,sl,n2,s2,ri) from a given set of (t,qv,ql,qi,u,v)

       call trbintd( &
                     pcols    , pver    , ncol  , z       , ufd     , vfd     , tfd   , pmid    , &
                     s2       , n2      , ri    , zi      , pi      , cldn    , qtfd  , qvfd    , &
                     qlfd     , qi      , sfi   , sfuh    , sflh    , slfd    , slv   , slslope , &
                     qtslope  , chs     , chu   , cms     , cmu     )

     ! Save initial (i.e., before iterative diffusion) profile of (qt,sl) at each iteration.         
     ! Only necessary for (qt,sl) not (u,v) because (qt,sl) are newly calculated variables. 

       if( iturb .eq. 1 ) then
           qt(:ncol,:) = qtfd(:ncol,:)
           sl(:ncol,:) = slfd(:ncol,:)
       endif

     ! Get free atmosphere exchange coefficients. This 'kvf' is not used in UW moist PBL scheme

       call austausch_atm( pcols, pver, ncol, ri, s2, kvf )

     ! Initialize kvh/kvm to send to caleddy, depending on model timestep and iteration number
     ! This is necessary for 'wstar-based' entrainment closure.

       call eddy_diff_kv_init(ncol, pcols, pver, iturb, kvinit, kvf, kvh_in, kvm_in, kvh_out, kvm_out, kvh, kvm)

     ! Calculate eddy diffusivity (kvh_out,kvm_out) and (tke,bprod,sprod) using
     ! a given (kvh,kvm) which are used only for initializing (bprod,sprod)  at
     ! the first part of caleddy. (bprod,sprod) are fully updated at the end of
     ! caleddy after calculating (kvh_out,kvm_out) 

       call caleddy( pcols     , pver      , ncol      ,                     &
                     slfd      , qtfd      , qlfd      , slv      ,ufd     , &
                     vfd       , pi        , z         , zi       ,          &
                     qflx      , shflx     , slslope   , qtslope  ,          &
                     chu       , chs       , cmu       , cms      ,sfuh    , &
                     sflh      , n2        , s2        , ri       ,rrho    , &
                     pblh      , ustar     ,                                 &
                     kvh       , kvm       , kvh_out   , kvm_out  ,          &
                     tpert     , qpert     , qrl       , kvf      , tke    , &
                     wstarent  , bprod     , sprod     , minpblh  , wpert  , &
                     tkes      , went      , turbtype  , sm_aw    ,          & 
                     kbase_o   , ktop_o    , ncvfin_o  ,                     &
                     kbase_mg  , ktop_mg   , ncvfin_mg ,                     &                  
                     kbase_f   , ktop_f    , ncvfin_f  ,                     &                  
                     wet       , web       , jtbu      , jbbu     ,          &
                     evhc      , jt2slv    , n2ht      , n2hb     ,          & 
                     lwp       , opt_depth , radinvfrac, radf     ,          &
                     wstar     , wstar3fact,                                 &
                     ebrk      , wbrk      , lbrk      , ricl     , ghcl   , & 
                     shcl      , smcl      , ghi       , shi      , smi    , &
                     rii       , lengi     , wcap      , pblhp    , cldn   , &
                     ipbl      , kpblh     , wsedl     , wsed)

     ! Calculate errorPBL to check whether PBL produced convergent solutions or not.

       if( iturb .eq. nturb ) then
           call eddy_diff_error_pbl(ncol, pcols, pver, kvh, kvh_out, errorPBL)
       end if

     ! Eddy diffusivities which will be used for the initialization of (bprod,
     ! sprod) in 'caleddy' at the next iteration step.

       if( iturb .gt. 1 .and. iturb .lt. nturb ) then
           call eddy_diff_kv_relax(ncol, pcols, pver, lambda, kvm, kvh, kvm_out, kvh_out)
       endif

     ! Set nonlocal terms to zero for flux diagnostics, since not used by caleddy.

       call eddy_diff_zero_nonlocal(ncol, pcols, pver, cgh, cgs)

       if( iturb .lt. nturb ) then

         ! Each time we diffuse the original state

           call eddy_diff_restore_fields(ncol, pcols, pver, sl, qt, u, v, slfd, qtfd, ufd, vfd)

         !------------------------------------------------------------------------ 
         !  Check to see if constituent dependent gas constant needed (WACCM-X)
         !------------------------------------------------------------------------
         if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
           rairi(:ncol,1) = rairv(:ncol,1,lchnk)
           do k = 2, pver
             do i = 1, ncol
               rairi(i,k) = 0.5_r8 * (rairv(i,k,lchnk)+rairv(i,k-1,lchnk))
             end do
           end do      
         else
           rairi(:ncol,:pver+1) = rair          
         endif

         ! Diffuse initial profile of each time step using a given (kvh_out,kvm_out)
         ! In the below 'compute_vdiff', (slfd,qtfd,ufd,vfd) are 'inout' variables.

         call compute_vdiff( lchnk   ,                                                  &
                             pcols   , pver     , 1        , ncol         , pmid      , &
                             pi      , pdel, rpdel    , t        , ztodt        , taux      , &
                             tauy    , shflx    , qflx     , ntop_turb    , nbot_turb , &
                             kvh_out , kvm_out  , kvh_out  , cgs          , cgh       , &
                             zi      , ksrftms  , zero     , fieldlist_wet, fieldlist_molec, &
                             ufd     , vfd      , qtfd     , slfd         ,             &
                             jnk1d   , jnk1d    , jnk2d    , jnk1d        , errstring , &
                             tauresx , tauresy  , 0        , cpairv(:,:,lchnk), rairi , .false. )

         call handle_errmsg(errstring, subname="compute_vdiff", &
              extra_msg="compute_vdiff called from eddy_diff")

         ! Retrieve (tfd,qvfd,qlfd) from (slfd,qtfd) in order to 
         ! use 'trbintd' at the next iteration.
          
          call eddy_diff_rebuild_thermo(ncol, pcols, pver, cpair, latvap, latsub, g, rair, slfd, qtfd, qi, z, pmid, &
               qlfd, qvfd, tfd)
       endif

     ! Debug 
     ! icol = phys_debug_col(lchnk) 
     ! if( icol > 0 .and. get_nstep() .ge. 1 ) then
     !     write(iulog,*) ' '
     !     write(iulog,*) 'eddy_diff debug at the end of iteration' 
     !     write(iulog,*) 't,     qv,     ql,     cld,     u,     v'
     !     do k = pver-3, pver
     !        write (iulog,*) k, tfd(icol,k), qvfd(icol,k), qlfd(icol,k), cldn(icol,k), ufd(icol,k), vfd(icol,k)
     !     end do
     ! endif
     ! Debug

    end do  ! End of 'iturb' iteration

    kvq(:ncol,:) = kvh_out(:ncol,:)

  ! Compute 'wstar' within the PBL for use in the future convection scheme.

    call eddy_diff_wstar_pbl(ncol, pcols, ncvmax, ipbl, wstar, wstarPBL)

    ! --------------------------------------------------------------- !
    ! Writing for detailed diagnostic analysis of UW moist PBL scheme !
    ! --------------------------------------------------------------- !

    call outfld( 'UW_errorPBL',    errorPBL,   pcols,   lchnk )

    call outfld( 'UW_n2',          n2,         pcols,   lchnk )
    call outfld( 'UW_s2',          s2,         pcols,   lchnk )
    call outfld( 'UW_ri',          ri,         pcols,   lchnk )

    call outfld( 'UW_sfuh',        sfuh,       pcols,   lchnk )
    call outfld( 'UW_sflh',        sflh,       pcols,   lchnk )
    call outfld( 'UW_sfi',         sfi,        pcols,   lchnk )

    call outfld( 'UW_cldn',        cldn,       pcols,   lchnk )
    call outfld( 'UW_qrl',         qrl,        pcols,   lchnk )
    call outfld( 'UW_ql',          qlfd,       pcols,   lchnk )

    call outfld( 'UW_chu',         chu,        pcols,   lchnk )
    call outfld( 'UW_chs',         chs,        pcols,   lchnk )
    call outfld( 'UW_cmu',         cmu,        pcols,   lchnk )
    call outfld( 'UW_cms',         cms,        pcols,   lchnk )

    call outfld( 'UW_tke',         tke,        pcols,   lchnk )
    call outfld( 'UW_wcap',        wcap,       pcols,   lchnk )
    call outfld( 'UW_bprod',       bprod,      pcols,   lchnk )
    call outfld( 'UW_sprod',       sprod,      pcols,   lchnk )

    call outfld( 'UW_kvh',         kvh_out,    pcols,   lchnk )
    call outfld( 'UW_kvm',         kvm_out,    pcols,   lchnk )

    call outfld( 'UW_pblh',        pblh,       pcols,   lchnk )
    call outfld( 'UW_pblhp',       pblhp,      pcols,   lchnk )
    call outfld( 'UW_tpert',       tpert,      pcols,   lchnk )
    call outfld( 'UW_qpert',       qpert,      pcols,   lchnk )
    call outfld( 'UW_wpert',       wpert,      pcols,   lchnk )

    call outfld( 'UW_ustar',       ustar,      pcols,   lchnk )
    call outfld( 'UW_tkes',        tkes,       pcols,   lchnk )
    call outfld( 'UW_minpblh',     minpblh,    pcols,   lchnk )

    call outfld( 'UW_turbtype',    real(turbtype,r8),   pcols,   lchnk )

    call outfld( 'UW_kbase_o',     kbase_o,    pcols,   lchnk )
    call outfld( 'UW_ktop_o',      ktop_o,     pcols,   lchnk )
    call outfld( 'UW_ncvfin_o',    ncvfin_o,   pcols,   lchnk )

    call outfld( 'UW_kbase_mg',    kbase_mg,   pcols,   lchnk )
    call outfld( 'UW_ktop_mg',     ktop_mg,    pcols,   lchnk )
    call outfld( 'UW_ncvfin_mg',   ncvfin_mg,  pcols,   lchnk )

    call outfld( 'UW_kbase_f',     kbase_f,    pcols,   lchnk )
    call outfld( 'UW_ktop_f',      ktop_f,     pcols,   lchnk )
    call outfld( 'UW_ncvfin_f',    ncvfin_f,   pcols,   lchnk ) 

    call outfld( 'UW_wet',         wet,        pcols,   lchnk )
    call outfld( 'UW_web',         web,        pcols,   lchnk )
    call outfld( 'UW_jtbu',        jtbu,       pcols,   lchnk )
    call outfld( 'UW_jbbu',        jbbu,       pcols,   lchnk )
    call outfld( 'UW_evhc',        evhc,       pcols,   lchnk )
    call outfld( 'UW_jt2slv',      jt2slv,     pcols,   lchnk )
    call outfld( 'UW_n2ht',        n2ht,       pcols,   lchnk )
    call outfld( 'UW_n2hb',        n2hb,       pcols,   lchnk )
    call outfld( 'UW_lwp',         lwp,        pcols,   lchnk )
    call outfld( 'UW_optdepth',    opt_depth,  pcols,   lchnk )
    call outfld( 'UW_radfrac',     radinvfrac, pcols,   lchnk )
    call outfld( 'UW_radf',        radf,       pcols,   lchnk )
    call outfld( 'UW_wstar',       wstar,      pcols,   lchnk )
    call outfld( 'UW_wstar3fact',  wstar3fact, pcols,   lchnk )
    call outfld( 'UW_ebrk',        ebrk,       pcols,   lchnk )
    call outfld( 'UW_wbrk',        wbrk,       pcols,   lchnk )
    call outfld( 'UW_lbrk',        lbrk,       pcols,   lchnk )
    call outfld( 'UW_ricl',        ricl,       pcols,   lchnk )
    call outfld( 'UW_ghcl',        ghcl,       pcols,   lchnk )
    call outfld( 'UW_shcl',        shcl,       pcols,   lchnk )
    call outfld( 'UW_smcl',        smcl,       pcols,   lchnk )

    call outfld( 'UW_gh',          ghi,        pcols,   lchnk )
    call outfld( 'UW_sh',          shi,        pcols,   lchnk )
    call outfld( 'UW_sm',          smi,        pcols,   lchnk )
    call outfld( 'UW_ria',         rii,        pcols,   lchnk )
    call outfld( 'UW_leng',        lengi,      pcols,   lchnk )

    call outfld( 'UW_wsed',        wsed,       pcols,   lchnk )

    return
    
  end subroutine compute_eddy_diff

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !
  
  subroutine sfdiag( pcols   , pver    , ncol    , qt      , ql      , sl      , &
                     pi      , pm      , zi      , cld     , sfi     , sfuh    , &
                     sflh    , slslope , qtslope )
    !----------------------------------------------------------------------- ! 
    !                                                                        !
    ! Purpose: Interface for calculating saturation fractions  at upper and  ! 
    !          lower-half layers, & interfaces for use by turbulence scheme  !
    !                                                                        !
    ! Method : Various but 'l' should be chosen for consistency.             !
    !                                                                        ! 
    ! Author : B. Stevens and C. Bretherton (August 2000)                    !
    !          Sungsu Park. August 2006.                                     !
    !                       May.   2008.                                     ! 
    !                                                                        !  
    ! S.Park : The computed saturation fractions are repeatedly              !
    !          used to compute buoyancy coefficients in'trbintd' & 'caleddy'.!  
    !----------------------------------------------------------------------- !

    implicit none       

    ! --------------- !
    ! Input arguments !
    ! --------------- !

    integer,  intent(in)  :: pcols               ! Number of atmospheric columns   
    integer,  intent(in)  :: pver                ! Number of atmospheric layers   
    integer,  intent(in)  :: ncol                ! Number of atmospheric columns   

    real(r8), intent(in)  :: sl(pcols,pver)      ! Liquid water static energy [ J/kg ]
    real(r8), intent(in)  :: qt(pcols,pver)      ! Total water specific humidity [ kg/kg ]
    real(r8), intent(in)  :: ql(pcols,pver)      ! Liquid water specific humidity [ kg/kg ]
    real(r8), intent(in)  :: pi(pcols,pver+1)    ! Interface pressures [ Pa ]
    real(r8), intent(in)  :: pm(pcols,pver)      ! Layer mid-point pressures [ Pa ]
    real(r8), intent(in)  :: zi(pcols,pver+1)    ! Interface heights [ m ]
    real(r8), intent(in)  :: cld(pcols,pver)     ! Stratiform cloud fraction [ fraction ]
    real(r8), intent(in)  :: slslope(pcols,pver) ! Slope of 'sl' in each layer
    real(r8), intent(in)  :: qtslope(pcols,pver) ! Slope of 'qt' in each layer

    ! ---------------- !
    ! Output arguments !
    ! ---------------- !

    real(r8), intent(out) :: sfi(pcols,pver+1)   ! Interfacial layer saturation fraction [ fraction ]
    real(r8), intent(out) :: sfuh(pcols,pver)    ! Saturation fraction in upper half-layer [ fraction ]
    real(r8), intent(out) :: sflh(pcols,pver)    ! Saturation fraction in lower half-layer [ fraction ]

    ! --------------- !
    ! Local Variables !
    ! --------------- !

    integer               :: i                   ! Longitude index
    integer               :: k                   ! Vertical index
    integer               :: km1                 ! k-1
    integer               :: status              ! Status returned by function calls
    real(r8)              :: sltop, slbot        ! sl at top/bot of grid layer
    real(r8)              :: qttop, qtbot        ! qt at top/bot of grid layer
    real(r8)              :: tltop, tlbot  ! Liquid water temperature at top/bot of grid layer
    real(r8)              :: qxtop, qxbot        ! Sat excess at top/bot of grid layer
    real(r8)              :: qxm                 ! Sat excess at midpoint
    real(r8)              :: es               ! Saturation vapor pressure
    real(r8)              :: qs               ! Saturation spec. humidity
    real(r8)              :: cldeff(pcols,pver)  ! Effective Cloud Fraction [ fraction ]

    ! ----------------------- !
    ! Main Computation Begins ! 
    ! ----------------------- !

    sfi(1:ncol,:)    = 0._r8
    sfuh(1:ncol,:)   = 0._r8
    sflh(1:ncol,:)   = 0._r8
    cldeff(1:ncol,:) = 0._r8

    select case (sftype)
    case ('d')
       ! ----------------------------------------------------------------------- !
       ! Simply use the given stratus fraction ('horizontal' cloud partitioning) !
       ! ----------------------------------------------------------------------- !
       do k = ntop_turb + 1, nbot_turb
          km1 = k - 1
          do i = 1, ncol
             sfuh(i,k) = cld(i,k)
             sflh(i,k) = cld(i,k)
             sfi(i,k)  = 0.5_r8 * ( sflh(i,km1) + min( sflh(i,km1), sfuh(i,k) ) )
          end do
       end do
       do i = 1, ncol
          sfi(i,pver+1) = sflh(i,pver) 
       end do
    case ('l')
       ! ------------------------------------------ !
       ! Use modified stratus fraction partitioning !
       ! ------------------------------------------ !
       do k = ntop_turb + 1, nbot_turb
          km1 = k - 1
          do i = 1, ncol
             cldeff(i,k) = cld(i,k)
             sfuh(i,k)   = cld(i,k)
             sflh(i,k)   = cld(i,k)
             if( ql(i,k) .lt. qmin ) then
                 sfuh(i,k) = 0._r8
                 sflh(i,k) = 0._r8
             end if
           ! Modification : The contribution of ice should be carefully considered.
             if( choice_evhc .eq. 'ramp' .or. choice_radf .eq. 'ramp' ) then 
                 cldeff(i,k) = cld(i,k) * min( ql(i,k) / qmin, 1._r8 )
                 sfuh(i,k)   = cldeff(i,k)
                 sflh(i,k)   = cldeff(i,k)
             elseif( choice_evhc .eq. 'maxi' .or. choice_radf .eq. 'maxi' ) then 
                 cldeff(i,k) = cld(i,k)
                 sfuh(i,k)   = cldeff(i,k)
                 sflh(i,k)   = cldeff(i,k)
             endif
           ! At the stratus top, take the minimum interfacial saturation fraction
             sfi(i,k) = 0.5_r8 * ( sflh(i,km1) + min( sfuh(i,k), sflh(i,km1) ) )
           ! Modification : Currently sfi at the top and surface interfaces are set to be zero.
           !                Also, sfuh and sflh in the top model layer is set to be zero.
           !                However, I may need to set 
           !                         do i = 1, ncol
           !                            sfi(i,pver+1) = sflh(i,pver) 
           !                         end do
           !                for treating surface-based fog. 
           ! OK. I added below block similar to the other cases.
          end do
       end do
       do i = 1, ncol
          sfi(i,pver+1) = sflh(i,pver)
       end do
    case ('u')
       ! ------------------------------------------------------------------------- !
       ! Use unsaturated buoyancy - since sfi, sfuh, sflh have already been zeroed !
       ! nothing more need be done for this case.                                  !
       ! ------------------------------------------------------------------------- !
    case ('z')
       ! ------------------------------------------------------------------------- !
       ! Calculate saturation fraction based on whether the air just above or just !
       ! below the interface is saturated, i.e. with vertical cloud partitioning.  !
       ! The saturation fraction of the interfacial layer between mid-points k and !
       ! k+1 is computed by averaging the saturation fraction   of the half-layers !
       ! above and below the interface,  with a special provision   for cloud tops !
       ! (more cloud in the half-layer below than in the half-layer above).In each !
       ! half-layer, vertical partitioning of  cloud based on the slopes diagnosed !
       ! above is used.     Loop down through the layers, computing the saturation !
       ! fraction in each half-layer (sfuh for upper half, sflh for lower half).   !
       ! Once sfuh(i,k) is computed, use with sflh(i,k-1) to determine  saturation !
       ! fraction sfi(i,k) for interfacial layer k-0.5.                            !
       ! This is 'not' chosen for full consistent treatment of stratus fraction in !
       ! all physics schemes.                                                      !
       ! ------------------------------------------------------------------------- !
       do k = ntop_turb + 1, nbot_turb
          km1 = k - 1
          do i = 1, ncol
           ! Compute saturation excess at the mid-point of layer k
             sltop    = sl(i,k) + slslope(i,k) * ( pi(i,k) - pm(i,k) )      
             qttop    = qt(i,k) + qtslope(i,k) * ( pi(i,k) - pm(i,k) )
             tltop = ( sltop - g * zi(i,k) ) / cpair 
             call qsat( tltop, pi(i,k), es, qs)
             qxtop    = qttop - qs
             slbot    = sl(i,k) + slslope(i,k) * ( pi(i,k+1) - pm(i,k) )      
             qtbot    = qt(i,k) + qtslope(i,k) * ( pi(i,k+1) - pm(i,k) )
             tlbot = ( slbot - g * zi(i,k+1) ) / cpair 
             call qsat( tlbot, pi(i,k+1), es, qs)
             qxbot    = qtbot - qs
             qxm      = qxtop + ( qxbot - qxtop ) * ( pm(i,k) - pi(i,k) ) / ( pi(i,k+1) - pi(i,k) )
           ! Find the saturation fraction sfuh(i,k) of the upper half of layer k.
             if( ( qxtop .lt. 0._r8 ) .and. ( qxm .lt. 0._r8 ) ) then
                   sfuh(i,k) = 0._r8 
             else if( ( qxtop .gt. 0._r8 ) .and. ( qxm .gt. 0._r8 ) ) then
                   sfuh(i,k) = 1._r8  
             else ! Either qxm < 0 and qxtop > 0 or vice versa
                   sfuh(i,k) = max( qxtop, qxm ) / abs( qxtop - qxm )
             end if
           ! Combine with sflh(i) (still for layer k-1) to get interfac layer saturation fraction
             sfi(i,k) = 0.5_r8 * ( sflh(i,k-1) + min( sflh(i,k-1), sfuh(i,k) ) )
           ! Update sflh to be for the lower half of layer k.             
             if( ( qxbot .lt. 0._r8 ) .and. ( qxm .lt. 0._r8 ) ) then
                   sflh(i,k) = 0._r8 
             else if( ( qxbot .gt. 0._r8 ) .and. ( qxm .gt. 0._r8 ) ) then
                   sflh(i,k) = 1._r8 
             else ! Either qxm < 0 and qxbot > 0 or vice versa
                   sflh(i,k) = max( qxbot, qxm ) / abs( qxbot - qxm )
             end if
          end do  ! i
       end do ! k
       do i = 1, ncol
          sfi(i,pver+1) = sflh(i,pver)  ! Saturation fraction in the lowest half-layer. 
       end do
    end select

  return
  end subroutine sfdiag
  
  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !
 
  subroutine trbintd( pcols   , pver    , ncol    ,                               &
                      z       , u       , v       ,                               &
                      t       , pmid    ,                                         &
                      s2      , n2      , ri      ,                               &
                      zi      , pi      , cld     ,                               &
                      qt      , qv      , ql      , qi      , sfi     , sfuh    , &
                      sflh    , sl      , slv     , slslope , qtslope ,           &
                      chs     , chu     , cms     , cmu     )
    !----------------------------------------------------------------------- !
    ! Purpose: Calculate buoyancy coefficients at all interfaces including   !
    !          surface. Also, computes the profiles of ( sl,qt,n2,s2,ri ).   !
    !          Note that (n2,s2,ri) are defined at each interfaces except    !
    !          surface.                                                      !
    !                                                                        !
    ! Author: B. Stevens  ( Extracted from pbldiff, August, 2000 )           !
    !         Sungsu Park ( August 2006, May. 2008 )                         !
    !----------------------------------------------------------------------- !

    implicit none

    ! --------------- !
    ! Input arguments !
    ! --------------- !

    integer,  intent(in)  :: pcols                            ! Number of atmospheric columns   
    integer,  intent(in)  :: pver                             ! Number of atmospheric layers   
    integer,  intent(in)  :: ncol                             ! Number of atmospheric columns
    real(r8), target, intent(in)  :: z(pcols,pver)            ! Layer mid-point height above surface [ m ]
    real(r8), target, intent(in)  :: u(pcols,pver)            ! Layer mid-point u [ m/s ]
    real(r8), target, intent(in)  :: v(pcols,pver)            ! Layer mid-point v [ m/s ]
    real(r8), target, intent(in)  :: t(pcols,pver)            ! Layer mid-point temperature [ K ]
    real(r8), target, intent(in)  :: pmid(pcols,pver)         ! Layer mid-point pressure [ Pa ]
    real(r8), intent(in)  :: zi(pcols,pver+1)                 ! Interface height [ m ]
    real(r8), intent(in)  :: pi(pcols,pver+1)                 ! Interface pressure [ Pa ]
    real(r8), target, intent(in)  :: cld(pcols,pver)          ! Stratus fraction
    real(r8), target, intent(in)  :: qv(pcols,pver)           ! Water vapor specific humidity [ kg/kg ]
    real(r8), target, intent(in)  :: ql(pcols,pver)           ! Liquid water specific humidity [ kg/kg ]
    real(r8), target, intent(in)  :: qi(pcols,pver)           ! Ice water specific humidity [ kg/kg ]

    ! ---------------- !
    ! Output arguments !
    ! ---------------- !

    real(r8), target, intent(out) :: s2(pcols,pver)           ! Interfacial ( except surface ) shear squared [ s-2 ]
    real(r8), target, intent(out) :: n2(pcols,pver)           ! Interfacial ( except surface ) buoyancy frequency [ s-2 ]
    real(r8), target, intent(out) :: ri(pcols,pver)           ! Interfacial ( except surface ) Richardson number, 'n2/s2'
 
    real(r8), target, intent(out) :: qt(pcols,pver)           ! Total specific humidity [ kg/kg ]
    real(r8), target, intent(out) :: sfi(pcols,pver+1)        ! Interfacial layer saturation fraction [ fraction ]
    real(r8), target, intent(out) :: sfuh(pcols,pver)         ! Saturation fraction in upper half-layer [ fraction ]
    real(r8), target, intent(out) :: sflh(pcols,pver)         ! Saturation fraction in lower half-layer [ fraction ]
    real(r8), target, intent(out) :: sl(pcols,pver)           ! Liquid water static energy [ J/kg ] 
    real(r8), target, intent(out) :: slv(pcols,pver)          ! Liquid water virtual static energy [ J/kg ]
   
    real(r8), target, intent(out) :: chu(pcols,pver+1)        ! Heat buoyancy coef for dry states at all interfaces, finally.
                                                              ! [ unit ? ]
    real(r8), target, intent(out) :: chs(pcols,pver+1)        ! heat buoyancy coef for sat states at all interfaces, finally.
                                                              ! [ unit ? ]
    real(r8), target, intent(out) :: cmu(pcols,pver+1)        ! Moisture buoyancy coef for dry states at all interfaces, finally.
                                                              ! [ unit ? ]
    real(r8), target, intent(out) :: cms(pcols,pver+1)        ! Moisture buoyancy coef for sat states at all interfaces, finally.
                                                              ! [ unit ? ]
    real(r8), target, intent(out) :: slslope(pcols,pver)      ! Slope of 'sl' in each layer
    real(r8), target, intent(out) :: qtslope(pcols,pver)      ! Slope of 'qt' in each layer
 
    ! --------------- !
    ! Local Variables !
    ! --------------- ! 

    integer               :: k                                ! Level index
    real(r8), target      :: qs(pcols,pver)                   ! Saturation specific humidity
    real(r8)              :: es(pcols,pver)                   ! Saturation vapor pressure
    real(r8), target      :: gam(pcols,pver)                  ! (l/cp)*(d(qs)/dT)
    real(r8), target      :: dsldp_b(pcols), dqtdp_b(pcols)   ! Slopes across interface below

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !

    ! Calculate conservative scalars (qt,sl,slv) and buoyancy coefficients at the layer mid-points.
    ! Note that 'ntop_turb = 1', 'nbot_turb = pver'

    do k = ntop_turb, nbot_turb
       call qsat( t(:ncol,k), pmid(:ncol,k), es(:ncol,k), qs(:ncol,k), gam=gam(:ncol,k))
    end do

    call eddy_diff_trbintd_core(ncol, pcols, pver, t, z, u, v, qv, ql, qi, gam, pmid, pi, zi, cld, qt, sl, slv, &
         slslope, qtslope, dsldp_b, dqtdp_b, chu, chs, cmu, cms, sfi, sfuh, sflh, n2, s2, ri)

  return

  end subroutine trbintd

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_trbintd_core_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (trbintd_core_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_TRBINTD_CORE_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_trbintd_core_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_trbintd_core_impl = .false.
    end if

    trbintd_core_impl_selected = .true.

    if (masterproc) then
       if (use_native_trbintd_core_impl) then
          write(iulog,*) 'eddy_diff_trbintd_core implementation = native'
       else
          write(iulog,*) 'eddy_diff_trbintd_core implementation = codon'
       end if
    end if

  end subroutine eddy_diff_trbintd_core_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_trbintd_core(ncol, pcols, pver, t_local, z_local, u_local, v_local, qv_local, ql_local, qi_local, &
       gam_local, pmid_local, pi_local, zi_local, cld_local, qt_local, sl_local, slv_local, slslope_local, qtslope_local, &
       dsldp_b_local, dqtdp_b_local, chu_local, chs_local, cmu_local, cms_local, sfi_local, sfuh_local, sflh_local, &
       n2_local, s2_local, ri_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: ncol, pcols, pver
    real(r8), target, intent(in) :: t_local(pcols,pver), z_local(pcols,pver), u_local(pcols,pver), v_local(pcols,pver)
    real(r8), target, intent(in) :: qv_local(pcols,pver), ql_local(pcols,pver), qi_local(pcols,pver), gam_local(pcols,pver)
    real(r8), target, intent(in) :: pmid_local(pcols,pver), cld_local(pcols,pver)
    real(r8), intent(in) :: pi_local(pcols,pver+1), zi_local(pcols,pver+1)
    real(r8), target, intent(out) :: qt_local(pcols,pver), sl_local(pcols,pver), slv_local(pcols,pver)
    real(r8), target, intent(out) :: slslope_local(pcols,pver), qtslope_local(pcols,pver)
    real(r8), target, intent(out) :: dsldp_b_local(pcols), dqtdp_b_local(pcols)
    real(r8), target, intent(out) :: chu_local(pcols,pver+1), chs_local(pcols,pver+1), cmu_local(pcols,pver+1), &
         cms_local(pcols,pver+1)
    real(r8), target, intent(out) :: sfi_local(pcols,pver+1), sfuh_local(pcols,pver), sflh_local(pcols,pver)
    real(r8), target, intent(out) :: n2_local(pcols,pver), s2_local(pcols,pver), ri_local(pcols,pver)

    interface
       subroutine eddy_diff_trbintd_core_codon(ncol_c, pcols_c, pver_c, cpair_c, latvap_c, latsub_c, g_c, zvir_c, ntzero_c, &
            t_p, z_p, u_p, v_p, qv_p, ql_p, qi_p, gam_p, pmid_p, cld_p, qt_p, sl_p, slv_p, slslope_p, qtslope_p, &
            dsldp_b_p, dqtdp_b_p, chu_p, chs_p, cmu_p, cms_p, sfi_p, sfuh_p, sflh_p, n2_p, s2_p, ri_p) &
            bind(c, name="eddy_diff_trbintd_core_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         real(c_double), value :: cpair_c, latvap_c, latsub_c, g_c, zvir_c, ntzero_c
         type(c_ptr), value :: t_p, z_p, u_p, v_p, qv_p, ql_p, qi_p, gam_p, pmid_p, cld_p, qt_p, sl_p, slv_p, &
              slslope_p, qtslope_p, dsldp_b_p, dqtdp_b_p, chu_p, chs_p, cmu_p, cms_p, sfi_p, sfuh_p, sflh_p, n2_p, s2_p, &
              ri_p
       end subroutine eddy_diff_trbintd_core_codon
    end interface

    call eddy_diff_trbintd_core_select_impl()

    if (use_native_trbintd_core_impl) then
       call eddy_diff_trbintd_core_native(ncol, pcols, pver, t_local, z_local, u_local, v_local, qv_local, ql_local, qi_local, &
            gam_local, pmid_local, pi_local, zi_local, cld_local, qt_local, sl_local, slv_local, slslope_local, qtslope_local, &
            dsldp_b_local, dqtdp_b_local, chu_local, chs_local, cmu_local, cms_local, sfi_local, sfuh_local, sflh_local, &
            n2_local, s2_local, ri_local)
       return
    end if

    call eddy_diff_trbintd_core_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), real(cpair, c_double), &
         real(latvap, c_double), real(latsub, c_double), real(g, c_double), real(zvir, c_double), real(ntzero, c_double), &
         c_loc(t_local), c_loc(z_local), c_loc(u_local), c_loc(v_local), c_loc(qv_local), c_loc(ql_local), c_loc(qi_local), &
         c_loc(gam_local), c_loc(pmid_local), c_loc(cld_local), c_loc(qt_local), c_loc(sl_local), c_loc(slv_local), &
         c_loc(slslope_local), c_loc(qtslope_local), c_loc(dsldp_b_local), c_loc(dqtdp_b_local), c_loc(chu_local), &
         c_loc(chs_local), c_loc(cmu_local), c_loc(cms_local), c_loc(sfi_local), c_loc(sfuh_local), c_loc(sflh_local), &
         c_loc(n2_local), c_loc(s2_local), c_loc(ri_local))

  end subroutine eddy_diff_trbintd_core

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_trbintd_core_native(ncol, pcols, pver, t_local, z_local, u_local, v_local, qv_local, ql_local, qi_local, &
       gam_local, pmid_local, pi_local, zi_local, cld_local, qt_local, sl_local, slv_local, slslope_local, qtslope_local, &
       dsldp_b_local, dqtdp_b_local, chu_local, chs_local, cmu_local, cms_local, sfi_local, sfuh_local, sflh_local, &
       n2_local, s2_local, ri_local)

    implicit none

    integer, intent(in) :: ncol, pcols, pver
    real(r8), intent(in) :: t_local(pcols,pver), z_local(pcols,pver), u_local(pcols,pver), v_local(pcols,pver)
    real(r8), intent(in) :: qv_local(pcols,pver), ql_local(pcols,pver), qi_local(pcols,pver), gam_local(pcols,pver)
    real(r8), intent(in) :: pmid_local(pcols,pver), pi_local(pcols,pver+1), zi_local(pcols,pver+1), cld_local(pcols,pver)
    real(r8), intent(out) :: qt_local(pcols,pver), sl_local(pcols,pver), slv_local(pcols,pver)
    real(r8), intent(out) :: slslope_local(pcols,pver), qtslope_local(pcols,pver)
    real(r8), intent(out) :: dsldp_b_local(pcols), dqtdp_b_local(pcols)
    real(r8), intent(out) :: chu_local(pcols,pver+1), chs_local(pcols,pver+1), cmu_local(pcols,pver+1), cms_local(pcols,pver+1)
    real(r8), intent(out) :: sfi_local(pcols,pver+1), sfuh_local(pcols,pver), sflh_local(pcols,pver)
    real(r8), intent(out) :: n2_local(pcols,pver), s2_local(pcols,pver), ri_local(pcols,pver)

    call eddy_diff_trbintd_midpoint_native(ncol, pcols, pver, t_local, z_local, qv_local, ql_local, qi_local, gam_local, &
         qt_local, sl_local, slv_local, chu_local, chs_local, cmu_local, cms_local)

    call eddy_diff_trbintd_slopes_native(ncol, pcols, pver, pmid_local, sl_local, qt_local, slslope_local, qtslope_local, &
         dsldp_b_local, dqtdp_b_local)

    call eddy_diff_trbintd_sfdiag_interface_native(ncol, pcols, pver, ql_local, sl_local, qt_local, pi_local, pmid_local, &
         zi_local, cld_local, slslope_local, qtslope_local, u_local, v_local, z_local, chu_local, chs_local, cmu_local, &
         cms_local, sfi_local, sfuh_local, sflh_local, n2_local, s2_local, ri_local)

  end subroutine eddy_diff_trbintd_core_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_trbintd_midpoint_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (trbintd_midpoint_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_TRBINTD_MIDPOINT_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_trbintd_midpoint_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_trbintd_midpoint_impl = .false.
    end if

    trbintd_midpoint_impl_selected = .true.

    if (masterproc) then
       if (use_native_trbintd_midpoint_impl) then
          write(iulog,*) 'eddy_diff_trbintd_midpoint implementation = native'
       else
          write(iulog,*) 'eddy_diff_trbintd_midpoint implementation = codon'
       end if
    end if

  end subroutine eddy_diff_trbintd_midpoint_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_trbintd_midpoint(ncol, pcols, pver, t_local, z_local, qv_local, ql_local, qi_local, gam_local, &
       qt_local, sl_local, slv_local, chu_local, chs_local, cmu_local, cms_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: ncol, pcols, pver
    real(r8), target, intent(in) :: t_local(pcols,pver), z_local(pcols,pver), qv_local(pcols,pver), ql_local(pcols,pver), &
         qi_local(pcols,pver), gam_local(pcols,pver)
    real(r8), target, intent(out) :: qt_local(pcols,pver), sl_local(pcols,pver), slv_local(pcols,pver)
    real(r8), target, intent(out) :: chu_local(pcols,pver+1), chs_local(pcols,pver+1), cmu_local(pcols,pver+1), &
         cms_local(pcols,pver+1)

    interface
       subroutine eddy_diff_trbintd_midpoint_codon(ncol_c, pcols_c, pver_c, cpair_c, latvap_c, latsub_c, g_c, zvir_c, &
            t_p, z_p, qv_p, ql_p, qi_p, gam_p, qt_p, sl_p, slv_p, chu_p, chs_p, cmu_p, cms_p) &
            bind(c, name="eddy_diff_trbintd_midpoint_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         real(c_double), value :: cpair_c, latvap_c, latsub_c, g_c, zvir_c
         type(c_ptr), value :: t_p, z_p, qv_p, ql_p, qi_p, gam_p, qt_p, sl_p, slv_p, chu_p, chs_p, cmu_p, cms_p
       end subroutine eddy_diff_trbintd_midpoint_codon
    end interface

    call eddy_diff_trbintd_midpoint_select_impl()

    if (use_native_trbintd_midpoint_impl) then
       call eddy_diff_trbintd_midpoint_native(ncol, pcols, pver, t_local, z_local, qv_local, ql_local, qi_local, gam_local, &
            qt_local, sl_local, slv_local, chu_local, chs_local, cmu_local, cms_local)
       return
    end if

    call eddy_diff_trbintd_midpoint_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), real(cpair, c_double), real(latvap, c_double), &
         real(latsub, c_double), real(g, c_double), real(zvir, c_double), c_loc(t_local), c_loc(z_local), c_loc(qv_local), &
         c_loc(ql_local), c_loc(qi_local), c_loc(gam_local), c_loc(qt_local), c_loc(sl_local), c_loc(slv_local), &
         c_loc(chu_local), c_loc(chs_local), c_loc(cmu_local), c_loc(cms_local) &
    )

  end subroutine eddy_diff_trbintd_midpoint

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_trbintd_midpoint_native(ncol, pcols, pver, t_local, z_local, qv_local, ql_local, qi_local, gam_local, &
       qt_local, sl_local, slv_local, chu_local, chs_local, cmu_local, cms_local)

    implicit none

    integer, intent(in) :: ncol, pcols, pver
    real(r8), intent(in) :: t_local(pcols,pver), z_local(pcols,pver), qv_local(pcols,pver), ql_local(pcols,pver), &
         qi_local(pcols,pver), gam_local(pcols,pver)
    real(r8), intent(out) :: qt_local(pcols,pver), sl_local(pcols,pver), slv_local(pcols,pver)
    real(r8), intent(out) :: chu_local(pcols,pver+1), chs_local(pcols,pver+1), cmu_local(pcols,pver+1), cms_local(pcols,pver+1)

    integer :: i, k
    real(r8) :: bfact_local

    do k = 1, pver
       do i = 1, ncol
          qt_local(i,k)  = qv_local(i,k) + ql_local(i,k) + qi_local(i,k)
          sl_local(i,k)  = cpair * t_local(i,k) + g * z_local(i,k) - latvap * ql_local(i,k) - latsub * qi_local(i,k)
          slv_local(i,k) = sl_local(i,k) * ( 1._r8 + zvir * qt_local(i,k) )
          bfact_local    = g / ( t_local(i,k) * ( 1._r8 + zvir * qv_local(i,k) - ql_local(i,k) - qi_local(i,k) ) )
          chu_local(i,k) = ( 1._r8 + zvir * qt_local(i,k) ) * bfact_local / cpair
          chs_local(i,k) = ( ( 1._r8 + ( 1._r8 + zvir ) * gam_local(i,k) * cpair * t_local(i,k) / latvap ) / &
               ( 1._r8 + gam_local(i,k) ) ) * bfact_local / cpair
          cmu_local(i,k) = zvir * bfact_local * t_local(i,k)
          cms_local(i,k) = latvap * chs_local(i,k)  -  bfact_local * t_local(i,k)
       end do
    end do

    do i = 1, ncol
       chu_local(i,pver+1) = chu_local(i,pver)
       chs_local(i,pver+1) = chs_local(i,pver)
       cmu_local(i,pver+1) = cmu_local(i,pver)
       cms_local(i,pver+1) = cms_local(i,pver)
    end do

  end subroutine eddy_diff_trbintd_midpoint_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_trbintd_slopes_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (trbintd_slopes_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_TRBINTD_SLOPES_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_trbintd_slopes_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_trbintd_slopes_impl = .false.
    end if

    trbintd_slopes_impl_selected = .true.

    if (masterproc) then
       if (use_native_trbintd_slopes_impl) then
          write(iulog,*) 'eddy_diff_trbintd_slopes implementation = native'
       else
          write(iulog,*) 'eddy_diff_trbintd_slopes implementation = codon'
       end if
    end if

  end subroutine eddy_diff_trbintd_slopes_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_trbintd_slopes(ncol, pcols, pver, pmid_local, sl_local, qt_local, slslope_local, qtslope_local, &
       dsldp_b_local, dqtdp_b_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: ncol, pcols, pver
    real(r8), target, intent(in) :: pmid_local(pcols,pver), sl_local(pcols,pver), qt_local(pcols,pver)
    real(r8), target, intent(out) :: slslope_local(pcols,pver), qtslope_local(pcols,pver)
    real(r8), target, intent(out) :: dsldp_b_local(pcols), dqtdp_b_local(pcols)

    interface
       subroutine eddy_diff_trbintd_slopes_codon(ncol_c, pcols_c, pver_c, pmid_p, sl_p, qt_p, slslope_p, qtslope_p, &
            dsldp_b_p, dqtdp_b_p) bind(c, name="eddy_diff_trbintd_slopes_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         type(c_ptr), value :: pmid_p, sl_p, qt_p, slslope_p, qtslope_p, dsldp_b_p, dqtdp_b_p
       end subroutine eddy_diff_trbintd_slopes_codon
    end interface

    call eddy_diff_trbintd_slopes_select_impl()

    if (use_native_trbintd_slopes_impl) then
       call eddy_diff_trbintd_slopes_native(ncol, pcols, pver, pmid_local, sl_local, qt_local, slslope_local, qtslope_local, &
            dsldp_b_local, dqtdp_b_local)
       return
    end if

    call eddy_diff_trbintd_slopes_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), c_loc(pmid_local), &
         c_loc(sl_local), c_loc(qt_local), c_loc(slslope_local), c_loc(qtslope_local), c_loc(dsldp_b_local), &
         c_loc(dqtdp_b_local))

  end subroutine eddy_diff_trbintd_slopes

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_trbintd_slopes_native(ncol, pcols, pver, pmid_local, sl_local, qt_local, slslope_local, qtslope_local, &
       dsldp_b_local, dqtdp_b_local)

    implicit none

    integer, intent(in) :: ncol, pcols, pver
    real(r8), intent(in) :: pmid_local(pcols,pver), sl_local(pcols,pver), qt_local(pcols,pver)
    real(r8), intent(out) :: slslope_local(pcols,pver), qtslope_local(pcols,pver)
    real(r8), intent(out) :: dsldp_b_local(pcols), dqtdp_b_local(pcols)

    integer :: i, k
    real(r8) :: product_local, dsldp_a_local, dqtdp_a_local

    do i = 1, ncol
       slslope_local(i,pver) = ( sl_local(i,pver) - sl_local(i,pver-1) ) / ( pmid_local(i,pver) - pmid_local(i,pver-1) )
       qtslope_local(i,pver) = ( qt_local(i,pver) - qt_local(i,pver-1) ) / ( pmid_local(i,pver) - pmid_local(i,pver-1) )
       slslope_local(i,1)    = ( sl_local(i,2) - sl_local(i,1) ) / ( pmid_local(i,2) - pmid_local(i,1) )
       qtslope_local(i,1)    = ( qt_local(i,2) - qt_local(i,1) ) / ( pmid_local(i,2) - pmid_local(i,1) )
       dsldp_b_local(i)      = slslope_local(i,1)
       dqtdp_b_local(i)      = qtslope_local(i,1)
    end do

    do k = 2, pver - 1
       do i = 1, ncol
          dsldp_a_local    = dsldp_b_local(i)
          dqtdp_a_local    = dqtdp_b_local(i)
          dsldp_b_local(i) = ( sl_local(i,k+1) - sl_local(i,k) ) / ( pmid_local(i,k+1) - pmid_local(i,k) )
          dqtdp_b_local(i) = ( qt_local(i,k+1) - qt_local(i,k) ) / ( pmid_local(i,k+1) - pmid_local(i,k) )
          product_local    = dsldp_a_local * dsldp_b_local(i)
          if( product_local .le. 0._r8 ) then
              slslope_local(i,k) = 0._r8
          else if( product_local .gt. 0._r8 .and. dsldp_a_local .lt. 0._r8 ) then
              slslope_local(i,k) = max( dsldp_a_local, dsldp_b_local(i) )
          else if( product_local .gt. 0._r8 .and. dsldp_a_local .gt. 0._r8 ) then
              slslope_local(i,k) = min( dsldp_a_local, dsldp_b_local(i) )
          end if
          product_local = dqtdp_a_local * dqtdp_b_local(i)
          if( product_local .le. 0._r8 ) then
              qtslope_local(i,k) = 0._r8
          else if( product_local .gt. 0._r8 .and. dqtdp_a_local .lt. 0._r8 ) then
              qtslope_local(i,k) = max( dqtdp_a_local, dqtdp_b_local(i) )
          else if( product_local .gt. 0._r8 .and. dqtdp_a_local .gt. 0._r8 ) then
              qtslope_local(i,k) = min( dqtdp_a_local, dqtdp_b_local(i) )
          end if
       end do
    end do

  end subroutine eddy_diff_trbintd_slopes_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_trbintd_sfdiag_interface_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (trbintd_sfdiag_interface_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_TRBINTD_SFDIAG_INTERFACE_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_trbintd_sfdiag_interface_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_trbintd_sfdiag_interface_impl = .false.
    end if

    trbintd_sfdiag_interface_impl_selected = .true.

    if (masterproc) then
       if (use_native_trbintd_sfdiag_interface_impl) then
          write(iulog,*) 'eddy_diff_trbintd_sfdiag_interface implementation = native'
       else
          write(iulog,*) 'eddy_diff_trbintd_sfdiag_interface implementation = codon'
       end if
    end if

  end subroutine eddy_diff_trbintd_sfdiag_interface_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_trbintd_sfdiag_interface(ncol, pcols, pver, ql_local, sl_local, qt_local, pi_local, pmid_local, &
       zi_local, cld_local, slslope_local, qtslope_local, u_local, v_local, z_local, chu_local, chs_local, cmu_local, &
       cms_local, sfi_local, sfuh_local, sflh_local, n2_local, s2_local, ri_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: ncol, pcols, pver
    real(r8), intent(in) :: ql_local(pcols,pver), pi_local(pcols,pver+1), pmid_local(pcols,pver), zi_local(pcols,pver+1)
    real(r8), intent(in) :: slslope_local(pcols,pver), qtslope_local(pcols,pver)
    real(r8), target, intent(in) :: sl_local(pcols,pver), qt_local(pcols,pver), cld_local(pcols,pver), u_local(pcols,pver), &
         v_local(pcols,pver), z_local(pcols,pver)
    real(r8), target, intent(inout) :: chu_local(pcols,pver+1), chs_local(pcols,pver+1), cmu_local(pcols,pver+1), &
         cms_local(pcols,pver+1)
    real(r8), target, intent(out) :: sfi_local(pcols,pver+1), sfuh_local(pcols,pver), sflh_local(pcols,pver)
    real(r8), target, intent(out) :: n2_local(pcols,pver), s2_local(pcols,pver), ri_local(pcols,pver)

    interface
       subroutine eddy_diff_trbintd_sfdiag_interface_codon(ncol_c, pcols_c, pver_c, ntzero_c, cld_p, u_p, v_p, z_p, sl_p, &
            qt_p, chu_p, chs_p, cmu_p, cms_p, sfi_p, sfuh_p, sflh_p, n2_p, s2_p, ri_p) &
            bind(c, name="eddy_diff_trbintd_sfdiag_interface_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         real(c_double), value :: ntzero_c
         type(c_ptr), value :: cld_p, u_p, v_p, z_p, sl_p, qt_p, chu_p, chs_p, cmu_p, cms_p, sfi_p, sfuh_p, sflh_p, &
              n2_p, s2_p, ri_p
       end subroutine eddy_diff_trbintd_sfdiag_interface_codon
    end interface

    call eddy_diff_trbintd_sfdiag_interface_select_impl()

    if (use_native_trbintd_sfdiag_interface_impl) then
       call eddy_diff_trbintd_sfdiag_interface_native(ncol, pcols, pver, ql_local, sl_local, qt_local, pi_local, pmid_local, &
            zi_local, cld_local, slslope_local, qtslope_local, u_local, v_local, z_local, chu_local, chs_local, cmu_local, &
            cms_local, sfi_local, sfuh_local, sflh_local, n2_local, s2_local, ri_local)
       return
    end if

    call eddy_diff_trbintd_sfdiag_interface_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
         real(ntzero, c_double), c_loc(cld_local), c_loc(u_local), c_loc(v_local), c_loc(z_local), c_loc(sl_local), &
         c_loc(qt_local), c_loc(chu_local), c_loc(chs_local), c_loc(cmu_local), c_loc(cms_local), c_loc(sfi_local), &
         c_loc(sfuh_local), c_loc(sflh_local), c_loc(n2_local), c_loc(s2_local), c_loc(ri_local))

  end subroutine eddy_diff_trbintd_sfdiag_interface

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_trbintd_sfdiag_interface_native(ncol, pcols, pver, ql_local, sl_local, qt_local, pi_local, pmid_local, &
       zi_local, cld_local, slslope_local, qtslope_local, u_local, v_local, z_local, chu_local, chs_local, cmu_local, &
       cms_local, sfi_local, sfuh_local, sflh_local, n2_local, s2_local, ri_local)

    implicit none

    integer, intent(in) :: ncol, pcols, pver
    real(r8), intent(in) :: ql_local(pcols,pver), sl_local(pcols,pver), qt_local(pcols,pver)
    real(r8), intent(in) :: pi_local(pcols,pver+1), pmid_local(pcols,pver), zi_local(pcols,pver+1), cld_local(pcols,pver)
    real(r8), intent(in) :: slslope_local(pcols,pver), qtslope_local(pcols,pver), u_local(pcols,pver), v_local(pcols,pver)
    real(r8), intent(in) :: z_local(pcols,pver)
    real(r8), intent(inout) :: chu_local(pcols,pver+1), chs_local(pcols,pver+1), cmu_local(pcols,pver+1), cms_local(pcols,pver+1)
    real(r8), intent(out) :: sfi_local(pcols,pver+1), sfuh_local(pcols,pver), sflh_local(pcols,pver)
    real(r8), intent(out) :: n2_local(pcols,pver), s2_local(pcols,pver), ri_local(pcols,pver)

    integer :: i, k, km1
    real(r8) :: rdz_local, dsldz_local, dqtdz_local, ch_local, cm_local

    call sfdiag(pcols, pver, ncol, qt_local, ql_local, sl_local, pi_local, pmid_local, zi_local, cld_local, sfi_local, &
         sfuh_local, sflh_local, slslope_local, qtslope_local)

    do k = nbot_turb, ntop_turb + 1, -1
       km1 = k - 1
       do i = 1, ncol
          rdz_local      = 1._r8 / ( z_local(i,km1) - z_local(i,k) )
          dsldz_local    = ( sl_local(i,km1) - sl_local(i,k) ) * rdz_local
          dqtdz_local    = ( qt_local(i,km1) - qt_local(i,k) ) * rdz_local
          chu_local(i,k) = ( chu_local(i,km1) + chu_local(i,k) ) * 0.5_r8
          chs_local(i,k) = ( chs_local(i,km1) + chs_local(i,k) ) * 0.5_r8
          cmu_local(i,k) = ( cmu_local(i,km1) + cmu_local(i,k) ) * 0.5_r8
          cms_local(i,k) = ( cms_local(i,km1) + cms_local(i,k) ) * 0.5_r8
          ch_local       = chu_local(i,k) * ( 1._r8 - sfi_local(i,k) ) + chs_local(i,k) * sfi_local(i,k)
          cm_local       = cmu_local(i,k) * ( 1._r8 - sfi_local(i,k) ) + cms_local(i,k) * sfi_local(i,k)
          n2_local(i,k)  = ch_local * dsldz_local + cm_local * dqtdz_local
          s2_local(i,k)  = ( ( u_local(i,km1) - u_local(i,k) )**2 + ( v_local(i,km1) - v_local(i,k) )**2 ) * rdz_local**2
          s2_local(i,k)  = max( ntzero, s2_local(i,k) )
          ri_local(i,k)  = n2_local(i,k) / s2_local(i,k)
       end do
    end do

    do i = 1, ncol
       n2_local(i,1) = n2_local(i,2)
       s2_local(i,1) = s2_local(i,2)
       ri_local(i,1) = ri_local(i,2)
    end do

  end subroutine eddy_diff_trbintd_sfdiag_interface_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine surface_stress_diag_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (surface_stress_diag_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_SURFACE_STRESS_DIAG_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_surface_stress_diag_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_surface_stress_diag_impl = .false.
    end if

    surface_stress_diag_impl_selected = .true.

    if (masterproc) then
       if (use_native_surface_stress_diag_impl) then
          write(iulog,*) 'eddy_diff_surface_stress_diag implementation = native'
       else
          write(iulog,*) 'eddy_diff_surface_stress_diag implementation = codon'
       end if
    end if

  end subroutine surface_stress_diag_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_surface_stress_diag(ncol, pcols, pver, tfd_local, pmid_local, taux_local, tauy_local, &
       ksrftms_local, ufd_local, vfd_local, rrho_local, ustar_local, minpblh_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: ncol, pcols, pver
    real(r8), target, intent(in) :: tfd_local(pcols,pver), pmid_local(pcols,pver)
    real(r8), target, intent(in) :: taux_local(pcols), tauy_local(pcols), ksrftms_local(pcols)
    real(r8), target, intent(in) :: ufd_local(pcols,pver), vfd_local(pcols,pver)
    real(r8), target, intent(inout) :: rrho_local(pcols), ustar_local(pcols), minpblh_local(pcols)

    interface
       subroutine eddy_diff_surface_stress_diag_codon(ncol_c, pcols_c, pver_c, rair_c, ustar_min_c, tfd_p, &
            pmid_p, taux_p, tauy_p, ksrftms_p, ufd_p, vfd_p, rrho_p, ustar_p, minpblh_p) &
            bind(c, name="eddy_diff_surface_stress_diag_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         real(c_double), value :: rair_c, ustar_min_c
         type(c_ptr), value :: tfd_p, pmid_p, taux_p, tauy_p, ksrftms_p, ufd_p, vfd_p, rrho_p, ustar_p, minpblh_p
       end subroutine eddy_diff_surface_stress_diag_codon
    end interface

    call surface_stress_diag_select_impl()

    if (use_native_surface_stress_diag_impl) then
       call eddy_diff_surface_stress_diag_native(ncol, pcols, pver, tfd_local, pmid_local, taux_local, tauy_local, &
            ksrftms_local, ufd_local, vfd_local, rrho_local, ustar_local, minpblh_local)
       return
    end if

    call eddy_diff_surface_stress_diag_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), real(rair, c_double), &
         real(0.01_r8, c_double), c_loc(tfd_local), c_loc(pmid_local), c_loc(taux_local), c_loc(tauy_local), &
         c_loc(ksrftms_local), c_loc(ufd_local), c_loc(vfd_local), c_loc(rrho_local), c_loc(ustar_local), &
         c_loc(minpblh_local) &
    )

  end subroutine eddy_diff_surface_stress_diag

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_surface_stress_diag_native(ncol, pcols, pver, tfd_local, pmid_local, taux_local, tauy_local, &
       ksrftms_local, ufd_local, vfd_local, rrho_local, ustar_local, minpblh_local)

    use pbl_utils, only: calc_ustar

    implicit none

    integer, intent(in) :: ncol, pcols, pver
    real(r8), intent(in) :: tfd_local(pcols,pver), pmid_local(pcols,pver)
    real(r8), intent(in) :: taux_local(pcols), tauy_local(pcols), ksrftms_local(pcols)
    real(r8), intent(in) :: ufd_local(pcols,pver), vfd_local(pcols,pver)
    real(r8), intent(inout) :: rrho_local(pcols), ustar_local(pcols), minpblh_local(pcols)

    call calc_ustar( tfd_local(:ncol,pver), pmid_local(:ncol,pver), &
         taux_local(:ncol) - ksrftms_local(:ncol) * ufd_local(:ncol,pver), & ! Zonal wind stress
         tauy_local(:ncol) - ksrftms_local(:ncol) * vfd_local(:ncol,pver), & ! Meridional wind stress
         rrho_local(:ncol), ustar_local(:ncol))
    minpblh_local(:ncol) = 100.0_r8 * ustar_local(:ncol)

  end subroutine eddy_diff_surface_stress_diag_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !
  
  subroutine austausch_atm_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (austausch_atm_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('AUSTAUSCH_ATM_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_austausch_atm_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_austausch_atm_impl = .false.
    end if

    austausch_atm_impl_selected = .true.

    if (masterproc) then
       if (use_native_austausch_atm_impl) then
          write(iulog,*) 'austausch_atm implementation = native'
       else
          write(iulog,*) 'austausch_atm implementation = codon'
       end if
    end if

  end subroutine austausch_atm_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine austausch_atm( pcols, pver, ncol, ri, s2, kvf )

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    integer,  intent(in)  :: pcols                ! Number of atmospheric columns
    integer,  intent(in)  :: pver                 ! Number of atmospheric layers
    integer,  intent(in)  :: ncol                 ! Number of atmospheric columns

    real(r8), target, intent(in)  :: s2(pcols,pver)       ! Shear squared
    real(r8), target, intent(in)  :: ri(pcols,pver)       ! Richardson no
    real(r8), target, intent(out) :: kvf(pcols,pver+1)    ! Eddy diffusivity for heat and tracers

    interface
       subroutine austausch_atm_codon(ncol_c, pcols_c, pver_c, ntop_turb_c, nbot_turb_c, zkmin_c, &
            ri_p, s2_p, ml2_p, kvf_p) bind(c, name="austausch_atm_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, ntop_turb_c, nbot_turb_c
         real(c_double), value :: zkmin_c
         type(c_ptr), value :: ri_p, s2_p, ml2_p, kvf_p
       end subroutine austausch_atm_codon
    end interface

    call austausch_atm_select_impl()

    if (use_native_austausch_atm_impl) then
       call austausch_atm_native(pcols, pver, ncol, ri, s2, kvf)
       return
    end if

    call austausch_atm_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(ntop_turb, c_int64_t), &
         int(nbot_turb, c_int64_t), real(zkmin, c_double), c_loc(ri), c_loc(s2), c_loc(ml2), c_loc(kvf) &
    )

  end subroutine austausch_atm

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine austausch_atm_native( pcols, pver, ncol, ri, s2, kvf )

    !---------------------------------------------------------------------- ! 
    !                                                                       !
    ! Purpose: Computes exchange coefficients for free turbulent flows.     !
    !          This is not used in the UW moist turbulence scheme.          !
    !                                                                       !
    ! Method:                                                               !
    !                                                                       !
    ! The free atmosphere diffusivities are based on standard mixing length !
    ! forms for the neutral diffusivity multiplied by functns of Richardson !
    ! number. K = l^2 * |dV/dz| * f(Ri). The same functions are used for    !
    ! momentum, potential temperature, and constitutents.                   !
    !                                                                       !
    ! The stable Richardson num function (Ri>0) is taken from Holtslag and  !
    ! Beljaars (1989), ECMWF proceedings. f = 1 / (1 + 10*Ri*(1 + 8*Ri))    !
    ! The unstable Richardson number function (Ri<0) is taken from  CCM1.   !
    ! f = sqrt(1 - 18*Ri)                                                   !
    !                                                                       !
    ! Author: B. Stevens (rewrite, August 2000)                             !
    !                                                                       !
    !---------------------------------------------------------------------- !
    implicit none
    
    ! --------------- ! 
    ! Input arguments !
    ! --------------- !

    integer,  intent(in)  :: pcols                ! Number of atmospheric columns   
    integer,  intent(in)  :: pver                 ! Number of atmospheric layers   
    integer,  intent(in)  :: ncol                 ! Number of atmospheric columns

    real(r8), intent(in)  :: s2(pcols,pver)       ! Shear squared
    real(r8), intent(in)  :: ri(pcols,pver)       ! Richardson no

    ! ---------------- !
    ! Output arguments !
    ! ---------------- !

    real(r8), intent(out) :: kvf(pcols,pver+1)    ! Eddy diffusivity for heat and tracers

    ! --------------- !
    ! Local Variables !
    ! --------------- !

    real(r8)              :: fofri                ! f(ri)
    real(r8)              :: kvn                  ! Neutral Kv

    integer               :: i                    ! Longitude index
    integer               :: k                    ! Vertical index

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !

    kvf(:ncol,:)           = 0.0_r8
    kvf(:ncol,pver+1)      = 0.0_r8
    kvf(:ncol,1:ntop_turb) = 0.0_r8

    ! Compute the free atmosphere vertical diffusion coefficients: kvh = kvq = kvm. 

    do k = ntop_turb + 1, nbot_turb
       do i = 1, ncol
          if( ri(i,k) < 0.0_r8 ) then
              fofri = sqrt( max( 1._r8 - 18._r8 * ri(i,k), 0._r8 ) )
          else 
              fofri = 1.0_r8 / ( 1.0_r8 + 10.0_r8 * ri(i,k) * ( 1.0_r8 + 8.0_r8 * ri(i,k) ) )    
          end if
          kvn = ml2(k) * sqrt(s2(i,k))
          kvf(i,k) = max( zkmin, kvn * fofri )
       end do
    end do

    return

    end subroutine austausch_atm_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_kv_init_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (kv_init_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_KV_INIT_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_kv_init_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_kv_init_impl = .false.
    end if

    kv_init_impl_selected = .true.

    if (masterproc) then
       if (use_native_kv_init_impl) then
          write(iulog,*) 'eddy_diff_kv_init implementation = native'
       else
          write(iulog,*) 'eddy_diff_kv_init implementation = codon'
       end if
    end if

  end subroutine eddy_diff_kv_init_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_kv_init(ncol, pcols, pver, iturb, kvinit_local, kvf_local, kvh_in_local, kvm_in_local, &
       kvh_out_local, kvm_out_local, kvh_local, kvm_local)

    use iso_c_binding, only: c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: ncol, pcols, pver, iturb
    logical, intent(in) :: kvinit_local
    real(r8), target, intent(in) :: kvf_local(pcols,pver+1), kvh_in_local(pcols,pver+1), kvm_in_local(pcols,pver+1)
    real(r8), target, intent(in) :: kvh_out_local(pcols,pver+1), kvm_out_local(pcols,pver+1)
    real(r8), target, intent(inout) :: kvh_local(pcols,pver+1), kvm_local(pcols,pver+1)

    interface
       subroutine eddy_diff_kv_init_codon(ncol_c, pcols_c, pver_c, iturb_c, kvinit_c, use_kvf_c, kvf_p, kvh_in_p, &
            kvm_in_p, kvh_out_p, kvm_out_p, kvh_p, kvm_p) bind(c, name="eddy_diff_kv_init_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, iturb_c, kvinit_c, use_kvf_c
         type(c_ptr), value :: kvf_p, kvh_in_p, kvm_in_p, kvh_out_p, kvm_out_p, kvh_p, kvm_p
       end subroutine eddy_diff_kv_init_codon
    end interface

    call eddy_diff_kv_init_select_impl()

    if (use_native_kv_init_impl) then
       call eddy_diff_kv_init_native(ncol, pcols, pver, iturb, kvinit_local, kvf_local, kvh_in_local, kvm_in_local, &
            kvh_out_local, kvm_out_local, kvh_local, kvm_local)
       return
    end if

    call eddy_diff_kv_init_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(iturb, c_int64_t), &
         int(merge(1, 0, kvinit_local), c_int64_t), int(merge(1, 0, use_kvf), c_int64_t), c_loc(kvf_local), &
         c_loc(kvh_in_local), c_loc(kvm_in_local), c_loc(kvh_out_local), c_loc(kvm_out_local), c_loc(kvh_local), &
         c_loc(kvm_local) &
    )

  end subroutine eddy_diff_kv_init

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_kv_init_native(ncol, pcols, pver, iturb, kvinit_local, kvf_local, kvh_in_local, kvm_in_local, &
       kvh_out_local, kvm_out_local, kvh_local, kvm_local)

    implicit none

    integer, intent(in) :: ncol, pcols, pver, iturb
    logical, intent(in) :: kvinit_local
    real(r8), intent(in) :: kvf_local(pcols,pver+1), kvh_in_local(pcols,pver+1), kvm_in_local(pcols,pver+1)
    real(r8), intent(in) :: kvh_out_local(pcols,pver+1), kvm_out_local(pcols,pver+1)
    real(r8), intent(inout) :: kvh_local(pcols,pver+1), kvm_local(pcols,pver+1)

    if( iturb .eq. 1 ) then
       if( kvinit_local ) then
       ! First iteration of first model timestep : Use free tropospheric value or zero.
         if( use_kvf ) then
             kvh_local(:ncol,:) = kvf_local(:ncol,:)
             kvm_local(:ncol,:) = kvf_local(:ncol,:)
         else
             kvh_local(:ncol,:) = 0._r8
             kvm_local(:ncol,:) = 0._r8
         endif
       else
       ! First iteration on any model timestep except the first : Use value from previous timestep
         kvh_local(:ncol,:) = kvh_in_local(:ncol,:)
         kvm_local(:ncol,:) = kvm_in_local(:ncol,:)
       endif
    else
     ! Not the first iteration : Use from previous iteration
      kvh_local(:ncol,:) = kvh_out_local(:ncol,:)
      kvm_local(:ncol,:) = kvm_out_local(:ncol,:)
    endif

  end subroutine eddy_diff_kv_init_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_error_pbl_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (error_pbl_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_ERROR_PBL_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_error_pbl_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_error_pbl_impl = .false.
    end if

    error_pbl_impl_selected = .true.

    if (masterproc) then
       if (use_native_error_pbl_impl) then
          write(iulog,*) 'eddy_diff_error_pbl implementation = native'
       else
          write(iulog,*) 'eddy_diff_error_pbl implementation = codon'
       end if
    end if

  end subroutine eddy_diff_error_pbl_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_error_pbl(ncol, pcols, pver, kvh_local, kvh_out_local, errorPBL_local)

    use iso_c_binding, only: c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: ncol, pcols, pver
    real(r8), target, intent(in) :: kvh_local(pcols,pver+1), kvh_out_local(pcols,pver+1)
    real(r8), target, intent(inout) :: errorPBL_local(pcols)

    interface
       subroutine eddy_diff_error_pbl_codon(ncol_c, pcols_c, pver_c, kvh_p, kvh_out_p, errorPBL_p) &
            bind(c, name="eddy_diff_error_pbl_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         type(c_ptr), value :: kvh_p, kvh_out_p, errorPBL_p
       end subroutine eddy_diff_error_pbl_codon
    end interface

    call eddy_diff_error_pbl_select_impl()

    if (use_native_error_pbl_impl) then
       call eddy_diff_error_pbl_native(ncol, pcols, pver, kvh_local, kvh_out_local, errorPBL_local)
       return
    end if

    call eddy_diff_error_pbl_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), c_loc(kvh_local), c_loc(kvh_out_local), &
         c_loc(errorPBL_local) &
    )

  end subroutine eddy_diff_error_pbl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_error_pbl_native(ncol, pcols, pver, kvh_local, kvh_out_local, errorPBL_local)

    implicit none

    integer, intent(in) :: ncol, pcols, pver
    real(r8), intent(in) :: kvh_local(pcols,pver+1), kvh_out_local(pcols,pver+1)
    real(r8), intent(inout) :: errorPBL_local(pcols)

    integer :: i, k

    do i = 1, ncol
       errorPBL_local(i) = 0._r8
       do k = 1, pver
          errorPBL_local(i) = errorPBL_local(i) + ( kvh_local(i,k) - kvh_out_local(i,k) )**2
       end do
       errorPBL_local(i) = sqrt(errorPBL_local(i)/pver)
    end do

  end subroutine eddy_diff_error_pbl_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_kv_relax_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (kv_relax_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_KV_RELAX_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_kv_relax_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_kv_relax_impl = .false.
    end if

    kv_relax_impl_selected = .true.

    if (masterproc) then
       if (use_native_kv_relax_impl) then
          write(iulog,*) 'eddy_diff_kv_relax implementation = native'
       else
          write(iulog,*) 'eddy_diff_kv_relax implementation = codon'
       end if
    end if

  end subroutine eddy_diff_kv_relax_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_kv_relax(ncol, pcols, pver, lambda_local, kvm_local, kvh_local, kvm_out_local, kvh_out_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: ncol, pcols, pver
    real(r8), intent(in) :: lambda_local
    real(r8), target, intent(in) :: kvm_local(pcols,pver+1), kvh_local(pcols,pver+1)
    real(r8), target, intent(inout) :: kvm_out_local(pcols,pver+1), kvh_out_local(pcols,pver+1)

    interface
       subroutine eddy_diff_kv_relax_codon(ncol_c, pcols_c, pver_c, lambda_c, kvm_p, kvh_p, kvm_out_p, kvh_out_p) &
            bind(c, name="eddy_diff_kv_relax_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         real(c_double), value :: lambda_c
         type(c_ptr), value :: kvm_p, kvh_p, kvm_out_p, kvh_out_p
       end subroutine eddy_diff_kv_relax_codon
    end interface

    call eddy_diff_kv_relax_select_impl()

    if (use_native_kv_relax_impl) then
       call eddy_diff_kv_relax_native(ncol, pcols, pver, lambda_local, kvm_local, kvh_local, kvm_out_local, kvh_out_local)
       return
    end if

    call eddy_diff_kv_relax_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), real(lambda_local, c_double), &
         c_loc(kvm_local), c_loc(kvh_local), c_loc(kvm_out_local), c_loc(kvh_out_local) &
    )

  end subroutine eddy_diff_kv_relax

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_kv_relax_native(ncol, pcols, pver, lambda_local, kvm_local, kvh_local, kvm_out_local, kvh_out_local)

    implicit none

    integer, intent(in) :: ncol, pcols, pver
    real(r8), intent(in) :: lambda_local
    real(r8), intent(in) :: kvm_local(pcols,pver+1), kvh_local(pcols,pver+1)
    real(r8), intent(inout) :: kvm_out_local(pcols,pver+1), kvh_out_local(pcols,pver+1)

    kvm_out_local(:ncol,:) = lambda_local * kvm_out_local(:ncol,:) + ( 1._r8 - lambda_local ) * kvm_local(:ncol,:)
    kvh_out_local(:ncol,:) = lambda_local * kvh_out_local(:ncol,:) + ( 1._r8 - lambda_local ) * kvh_local(:ncol,:)

  end subroutine eddy_diff_kv_relax_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_init_fields_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (init_fields_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_INIT_FIELDS_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_init_fields_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_init_fields_impl = .false.
    end if

    init_fields_impl_selected = .true.

    if (masterproc) then
       if (use_native_init_fields_impl) then
          write(iulog,*) 'eddy_diff_init_fields implementation = native'
       else
          write(iulog,*) 'eddy_diff_init_fields implementation = codon'
       end if
    end if

  end subroutine eddy_diff_init_fields_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_init_fields(ncol, pcols, pver, u_local, v_local, t_local, qv_local, ql_local, zero_local, &
       zero2d_local, ufd_local, vfd_local, tfd_local, qvfd_local, qlfd_local)

    use iso_c_binding, only: c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: ncol, pcols, pver
    real(r8), target, intent(in) :: u_local(pcols,pver), v_local(pcols,pver), t_local(pcols,pver), qv_local(pcols,pver), &
         ql_local(pcols,pver)
    real(r8), target, intent(inout) :: zero_local(pcols), zero2d_local(pcols,pver+1), ufd_local(pcols,pver), &
         vfd_local(pcols,pver), tfd_local(pcols,pver), qvfd_local(pcols,pver), qlfd_local(pcols,pver)

    interface
       subroutine eddy_diff_init_fields_codon(ncol_c, pcols_c, pver_c, u_p, v_p, t_p, qv_p, ql_p, zero_p, zero2d_p, &
            ufd_p, vfd_p, tfd_p, qvfd_p, qlfd_p) bind(c, name="eddy_diff_init_fields_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         type(c_ptr), value :: u_p, v_p, t_p, qv_p, ql_p, zero_p, zero2d_p, ufd_p, vfd_p, tfd_p, qvfd_p, qlfd_p
       end subroutine eddy_diff_init_fields_codon
    end interface

    call eddy_diff_init_fields_select_impl()

    if (use_native_init_fields_impl) then
       call eddy_diff_init_fields_native(ncol, pcols, pver, u_local, v_local, t_local, qv_local, ql_local, zero_local, &
            zero2d_local, ufd_local, vfd_local, tfd_local, qvfd_local, qlfd_local)
       return
    end if

    call eddy_diff_init_fields_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), c_loc(u_local), c_loc(v_local), &
         c_loc(t_local), c_loc(qv_local), c_loc(ql_local), c_loc(zero_local), c_loc(zero2d_local), c_loc(ufd_local), &
         c_loc(vfd_local), c_loc(tfd_local), c_loc(qvfd_local), c_loc(qlfd_local) &
    )

  end subroutine eddy_diff_init_fields

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_init_fields_native(ncol, pcols, pver, u_local, v_local, t_local, qv_local, ql_local, zero_local, &
       zero2d_local, ufd_local, vfd_local, tfd_local, qvfd_local, qlfd_local)

    implicit none

    integer, intent(in) :: ncol, pcols, pver
    real(r8), intent(in) :: u_local(pcols,pver), v_local(pcols,pver), t_local(pcols,pver), qv_local(pcols,pver), &
         ql_local(pcols,pver)
    real(r8), intent(inout) :: zero_local(pcols), zero2d_local(pcols,pver+1), ufd_local(pcols,pver), vfd_local(pcols,pver), &
         tfd_local(pcols,pver), qvfd_local(pcols,pver), qlfd_local(pcols,pver)

    zero_local(:) = 0._r8
    zero2d_local(:,:) = 0._r8
    ufd_local(:ncol,:) = u_local(:ncol,:)
    vfd_local(:ncol,:) = v_local(:ncol,:)
    tfd_local(:ncol,:) = t_local(:ncol,:)
    qvfd_local(:ncol,:) = qv_local(:ncol,:)
    qlfd_local(:ncol,:) = ql_local(:ncol,:)

  end subroutine eddy_diff_init_fields_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  function eddy_diff_estblf_cb(t_local) bind(C, name="eddy_diff_estblf_cb") result(es_local)

    use iso_c_binding, only: c_double

    implicit none

    real(c_double), value, intent(in) :: t_local
    real(c_double) :: es_local

    es_local = estblf(real(t_local, r8))

  end function eddy_diff_estblf_cb

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  function eddy_diff_svp_to_qsat_cb(es_local, p_local) bind(C, name="eddy_diff_svp_to_qsat_cb") result(qs_local)

    use iso_c_binding, only: c_double

    implicit none

    real(c_double), value, intent(in) :: es_local, p_local
    real(c_double) :: qs_local

    qs_local = svp_to_qsat(real(es_local, r8), real(p_local, r8))

  end function eddy_diff_svp_to_qsat_cb

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_rebuild_thermo_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (rebuild_thermo_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_REBUILD_THERMO_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_rebuild_thermo_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_rebuild_thermo_impl = .false.
    end if

    rebuild_thermo_impl_selected = .true.

    if (masterproc) then
       if (use_native_rebuild_thermo_impl) then
          write(iulog,*) 'eddy_diff_rebuild_thermo implementation = native'
       else
          write(iulog,*) 'eddy_diff_rebuild_thermo implementation = codon'
       end if
    end if

  end subroutine eddy_diff_rebuild_thermo_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_rebuild_thermo(ncol, pcols, pver, cpair_local, latvap_local, latsub_local, g_local, rair_local, &
       slfd_local, qtfd_local, qi_local, z_local, pmid_local, qlfd_local, qvfd_local, tfd_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: ncol, pcols, pver
    real(r8), intent(in) :: cpair_local, latvap_local, latsub_local, g_local, rair_local
    real(r8), target, intent(in) :: slfd_local(pcols,pver), qtfd_local(pcols,pver), qi_local(pcols,pver), z_local(pcols,pver), &
         pmid_local(pcols,pver)
    real(r8), target, intent(inout) :: qlfd_local(pcols,pver), qvfd_local(pcols,pver), tfd_local(pcols,pver)

    interface
       subroutine eddy_diff_rebuild_thermo_codon(ncol_c, pcols_c, pver_c, cpair_c, latvap_c, latsub_c, g_c, rair_c, &
            slfd_p, qtfd_p, qi_p, z_p, pmid_p, qlfd_p, qvfd_p, tfd_p) bind(c, name="eddy_diff_rebuild_thermo_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         real(c_double), value :: cpair_c, latvap_c, latsub_c, g_c, rair_c
         type(c_ptr), value :: slfd_p, qtfd_p, qi_p, z_p, pmid_p, qlfd_p, qvfd_p, tfd_p
       end subroutine eddy_diff_rebuild_thermo_codon
    end interface

    call eddy_diff_rebuild_thermo_select_impl()

    if (use_native_rebuild_thermo_impl) then
       call eddy_diff_rebuild_thermo_native(ncol, pcols, pver, cpair_local, latvap_local, latsub_local, g_local, &
            rair_local, slfd_local, qtfd_local, qi_local, z_local, pmid_local, qlfd_local, qvfd_local, tfd_local)
       return
    end if

    call eddy_diff_rebuild_thermo_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), real(cpair_local, c_double), &
         real(latvap_local, c_double), real(latsub_local, c_double), real(g_local, c_double), real(rair_local, c_double), &
         c_loc(slfd_local), c_loc(qtfd_local), c_loc(qi_local), c_loc(z_local), c_loc(pmid_local), c_loc(qlfd_local), &
         c_loc(qvfd_local), c_loc(tfd_local) &
    )

  end subroutine eddy_diff_rebuild_thermo

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_rebuild_thermo_native(ncol, pcols, pver, cpair_local, latvap_local, latsub_local, g_local, &
       rair_local, slfd_local, qtfd_local, qi_local, z_local, pmid_local, qlfd_local, qvfd_local, tfd_local)

    implicit none

    integer, intent(in) :: ncol, pcols, pver
    real(r8), intent(in) :: cpair_local, latvap_local, latsub_local, g_local, rair_local
    real(r8), intent(in) :: slfd_local(pcols,pver), qtfd_local(pcols,pver), qi_local(pcols,pver), z_local(pcols,pver), &
         pmid_local(pcols,pver)
    real(r8), intent(inout) :: qlfd_local(pcols,pver), qvfd_local(pcols,pver), tfd_local(pcols,pver)

    integer :: i, k
    real(r8) :: es_local, qs_local, templ_local, temps_local, ep2_local

    do k = 1, pver
       do i = 1, ncol
          templ_local = ( slfd_local(i,k) - g_local*z_local(i,k) ) / cpair_local
          call qsat( templ_local, pmid_local(i,k), es_local, qs_local)
          ep2_local = .622_r8
          temps_local = templ_local + ( qtfd_local(i,k) - qs_local ) / ( cpair_local / latvap_local + &
               latvap_local * qs_local / ( rair_local * templ_local**2 ) )
          call qsat( temps_local, pmid_local(i,k), es_local, qs_local)
          qlfd_local(i,k) = max( qtfd_local(i,k) - qi_local(i,k) - qs_local ,0._r8 )
          qvfd_local(i,k) = max( 0._r8, qtfd_local(i,k) - qi_local(i,k) - qlfd_local(i,k) )
          tfd_local(i,k) = ( slfd_local(i,k) + latvap_local * qlfd_local(i,k) + latsub_local * qi_local(i,k) - &
               g_local*z_local(i,k)) / cpair_local
       end do
    end do

  end subroutine eddy_diff_rebuild_thermo_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zero_nonlocal_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (zero_nonlocal_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_ZERO_NONLOCAL_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_zero_nonlocal_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_zero_nonlocal_impl = .false.
    end if

    zero_nonlocal_impl_selected = .true.

    if (masterproc) then
       if (use_native_zero_nonlocal_impl) then
          write(iulog,*) 'eddy_diff_zero_nonlocal implementation = native'
       else
          write(iulog,*) 'eddy_diff_zero_nonlocal implementation = codon'
       end if
    end if

  end subroutine eddy_diff_zero_nonlocal_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zero_nonlocal(ncol, pcols, pver, cgh_local, cgs_local)

    use iso_c_binding, only: c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: ncol, pcols, pver
    real(r8), target, intent(inout) :: cgh_local(pcols,pver+1), cgs_local(pcols,pver+1)

    interface
       subroutine eddy_diff_zero_nonlocal_codon(ncol_c, pcols_c, pver_c, cgh_p, cgs_p) &
            bind(c, name="eddy_diff_zero_nonlocal_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         type(c_ptr), value :: cgh_p, cgs_p
       end subroutine eddy_diff_zero_nonlocal_codon
    end interface

    call eddy_diff_zero_nonlocal_select_impl()

    if (use_native_zero_nonlocal_impl) then
       call eddy_diff_zero_nonlocal_native(ncol, pcols, pver, cgh_local, cgs_local)
       return
    end if

    call eddy_diff_zero_nonlocal_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), c_loc(cgh_local), c_loc(cgs_local) &
    )

  end subroutine eddy_diff_zero_nonlocal

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zero_nonlocal_native(ncol, pcols, pver, cgh_local, cgs_local)

    implicit none

    integer, intent(in) :: ncol, pcols, pver
    real(r8), intent(inout) :: cgh_local(pcols,pver+1), cgs_local(pcols,pver+1)

    cgh_local(:ncol,:) = 0._r8
    cgs_local(:ncol,:) = 0._r8

  end subroutine eddy_diff_zero_nonlocal_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_restore_fields_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (restore_fields_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_RESTORE_FIELDS_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_restore_fields_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_restore_fields_impl = .false.
    end if

    restore_fields_impl_selected = .true.

    if (masterproc) then
       if (use_native_restore_fields_impl) then
          write(iulog,*) 'eddy_diff_restore_fields implementation = native'
       else
          write(iulog,*) 'eddy_diff_restore_fields implementation = codon'
       end if
    end if

  end subroutine eddy_diff_restore_fields_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_restore_fields(ncol, pcols, pver, sl_local, qt_local, u_local, v_local, slfd_local, qtfd_local, &
       ufd_local, vfd_local)

    use iso_c_binding, only: c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: ncol, pcols, pver
    real(r8), target, intent(in) :: sl_local(pcols,pver), qt_local(pcols,pver), u_local(pcols,pver), v_local(pcols,pver)
    real(r8), target, intent(inout) :: slfd_local(pcols,pver), qtfd_local(pcols,pver), ufd_local(pcols,pver), &
         vfd_local(pcols,pver)

    interface
       subroutine eddy_diff_restore_fields_codon(ncol_c, pcols_c, pver_c, sl_p, qt_p, u_p, v_p, slfd_p, qtfd_p, ufd_p, &
            vfd_p) bind(c, name="eddy_diff_restore_fields_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         type(c_ptr), value :: sl_p, qt_p, u_p, v_p, slfd_p, qtfd_p, ufd_p, vfd_p
       end subroutine eddy_diff_restore_fields_codon
    end interface

    call eddy_diff_restore_fields_select_impl()

    if (use_native_restore_fields_impl) then
       call eddy_diff_restore_fields_native(ncol, pcols, pver, sl_local, qt_local, u_local, v_local, slfd_local, &
            qtfd_local, ufd_local, vfd_local)
       return
    end if

    call eddy_diff_restore_fields_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), c_loc(sl_local), c_loc(qt_local), &
         c_loc(u_local), c_loc(v_local), c_loc(slfd_local), c_loc(qtfd_local), c_loc(ufd_local), c_loc(vfd_local) &
    )

  end subroutine eddy_diff_restore_fields

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_restore_fields_native(ncol, pcols, pver, sl_local, qt_local, u_local, v_local, slfd_local, &
       qtfd_local, ufd_local, vfd_local)

    implicit none

    integer, intent(in) :: ncol, pcols, pver
    real(r8), intent(in) :: sl_local(pcols,pver), qt_local(pcols,pver), u_local(pcols,pver), v_local(pcols,pver)
    real(r8), intent(inout) :: slfd_local(pcols,pver), qtfd_local(pcols,pver), ufd_local(pcols,pver), vfd_local(pcols,pver)

    slfd_local(:ncol,:) = sl_local(:ncol,:)
    qtfd_local(:ncol,:) = qt_local(:ncol,:)
    ufd_local(:ncol,:) = u_local(:ncol,:)
    vfd_local(:ncol,:) = v_local(:ncol,:)

  end subroutine eddy_diff_restore_fields_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_wstar_pbl_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (wstar_pbl_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_WSTAR_PBL_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_wstar_pbl_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_wstar_pbl_impl = .false.
    end if

    wstar_pbl_impl_selected = .true.

    if (masterproc) then
       if (use_native_wstar_pbl_impl) then
          write(iulog,*) 'eddy_diff_wstar_pbl implementation = native'
       else
          write(iulog,*) 'eddy_diff_wstar_pbl implementation = codon'
       end if
    end if

  end subroutine eddy_diff_wstar_pbl_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_wstar_pbl(ncol, pcols, ncvmax_local, ipbl_local, wstar_local, wstarPBL_local)

    use iso_c_binding, only: c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: ncol, pcols, ncvmax_local
    integer(i4), target, intent(in) :: ipbl_local(pcols)
    real(r8), target, intent(in) :: wstar_local(pcols,ncvmax_local)
    real(r8), target, intent(inout) :: wstarPBL_local(pcols)

    interface
       subroutine eddy_diff_wstar_pbl_codon(ncol_c, pcols_c, ncvmax_c, ipbl_p, wstar_p, wstarPBL_p) &
            bind(c, name="eddy_diff_wstar_pbl_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, ncvmax_c
         type(c_ptr), value :: ipbl_p, wstar_p, wstarPBL_p
       end subroutine eddy_diff_wstar_pbl_codon
    end interface

    call eddy_diff_wstar_pbl_select_impl()

    if (use_native_wstar_pbl_impl) then
       call eddy_diff_wstar_pbl_native(ncol, pcols, ncvmax_local, ipbl_local, wstar_local, wstarPBL_local)
       return
    end if

    call eddy_diff_wstar_pbl_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(ncvmax_local, c_int64_t), c_loc(ipbl_local), c_loc(wstar_local), &
         c_loc(wstarPBL_local) &
    )

  end subroutine eddy_diff_wstar_pbl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_wstar_pbl_native(ncol, pcols, ncvmax_local, ipbl_local, wstar_local, wstarPBL_local)

    implicit none

    integer, intent(in) :: ncol, pcols, ncvmax_local
    integer(i4), intent(in) :: ipbl_local(pcols)
    real(r8), intent(in) :: wstar_local(pcols,ncvmax_local)
    real(r8), intent(inout) :: wstarPBL_local(pcols)

    integer :: i

    do i = 1, ncol
       if(ipbl_local(i) .eq. 1) then
           wstarPBL_local(i) = max( 0._r8, wstar_local(i,1) )
       else
           wstarPBL_local(i) = 0._r8
       endif
    end do

  end subroutine eddy_diff_wstar_pbl_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_init_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (caleddy_init_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_CALEDDY_INIT_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_caleddy_init_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_caleddy_init_impl = .false.
    end if

    caleddy_init_impl_selected = .true.

    if (masterproc) then
       if (use_native_caleddy_init_impl) then
          write(iulog,*) 'eddy_diff_caleddy_init implementation = native'
       else
          write(iulog,*) 'eddy_diff_caleddy_init implementation = codon'
       end if
    end if

  end subroutine eddy_diff_caleddy_init_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_init(ncol_local, pcols_local, pver_local, qrlzero_mode_local, cldeff_mode_local, &
       tkes_mode_local, use_kvf_mode_local, qmin_local, vk_local, ql_local, qrlin_local, cld_local, kvf_local, &
       kvh_in_local, kvm_in_local, n2_local, s2_local, shflx_local, qflx_local, rrho_local, ustar_local, z_local, &
       chu_local, chs_local, cmu_local, cms_local, sflh_local, qrlw_local, cldeff_local, kvh_local, kvm_local, bflxs_local, &
       bprod_local, sprod_local, wcap_local, leng_local, tke_local, turbtype_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: ncol_local, pcols_local, pver_local
    integer, intent(in) :: qrlzero_mode_local, cldeff_mode_local, tkes_mode_local, use_kvf_mode_local
    real(r8), intent(in) :: qmin_local, vk_local
    real(r8), target, intent(in) :: ql_local(pcols_local,pver_local), qrlin_local(pcols_local,pver_local)
    real(r8), target, intent(in) :: cld_local(pcols_local,pver_local), kvf_local(pcols_local,pver_local+1)
    real(r8), target, intent(in) :: kvh_in_local(pcols_local,pver_local+1), kvm_in_local(pcols_local,pver_local+1)
    real(r8), target, intent(in) :: n2_local(pcols_local,pver_local), s2_local(pcols_local,pver_local)
    real(r8), target, intent(in) :: shflx_local(pcols_local), qflx_local(pcols_local), rrho_local(pcols_local)
    real(r8), target, intent(in) :: ustar_local(pcols_local), z_local(pcols_local,pver_local)
    real(r8), target, intent(in) :: chu_local(pcols_local,pver_local+1), chs_local(pcols_local,pver_local+1)
    real(r8), target, intent(in) :: cmu_local(pcols_local,pver_local+1), cms_local(pcols_local,pver_local+1)
    real(r8), target, intent(in) :: sflh_local(pcols_local,pver_local)
    real(r8), target, intent(inout) :: qrlw_local(pcols_local,pver_local), cldeff_local(pcols_local,pver_local)
    real(r8), target, intent(inout) :: kvh_local(pcols_local,pver_local+1), kvm_local(pcols_local,pver_local+1)
    real(r8), target, intent(inout) :: bflxs_local(pcols_local), bprod_local(pcols_local,pver_local+1)
    real(r8), target, intent(inout) :: sprod_local(pcols_local,pver_local+1), wcap_local(pcols_local,pver_local+1)
    real(r8), target, intent(inout) :: leng_local(pcols_local,pver_local+1), tke_local(pcols_local,pver_local+1)
    integer(i4), target, intent(inout) :: turbtype_local(pcols_local,pver_local+1)

    interface
       subroutine eddy_diff_caleddy_init_codon(ncol_c, pcols_c, pver_c, qrlzero_mode_c, cldeff_mode_c, tkes_mode_c, &
            use_kvf_mode_c, qmin_c, vk_c, ql_p, qrlin_p, cld_p, kvf_p, kvh_in_p, kvm_in_p, n2_p, s2_p, shflx_p, qflx_p, &
            rrho_p, ustar_p, z_p, chu_p, chs_p, cmu_p, cms_p, sflh_p, qrlw_p, cldeff_p, kvh_p, kvm_p, bflxs_p, bprod_p, &
            sprod_p, wcap_p, leng_p, tke_p, turbtype_p) bind(c, name="eddy_diff_caleddy_init_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         integer(c_int64_t), value :: qrlzero_mode_c, cldeff_mode_c, tkes_mode_c, use_kvf_mode_c
         real(c_double), value :: qmin_c, vk_c
         type(c_ptr), value :: ql_p, qrlin_p, cld_p, kvf_p, kvh_in_p, kvm_in_p, n2_p, s2_p, shflx_p, qflx_p, rrho_p
         type(c_ptr), value :: ustar_p, z_p, chu_p, chs_p, cmu_p, cms_p, sflh_p, qrlw_p, cldeff_p, kvh_p, kvm_p
         type(c_ptr), value :: bflxs_p, bprod_p, sprod_p, wcap_p, leng_p, tke_p, turbtype_p
       end subroutine eddy_diff_caleddy_init_codon
    end interface

    call eddy_diff_caleddy_init_select_impl()

    if (use_native_caleddy_init_impl) then
       call eddy_diff_caleddy_init_native(ncol_local, pcols_local, pver_local, qrlzero_mode_local, cldeff_mode_local, &
            tkes_mode_local, use_kvf_mode_local, qmin_local, vk_local, ql_local, qrlin_local, cld_local, kvf_local, &
            kvh_in_local, kvm_in_local, n2_local, s2_local, shflx_local, qflx_local, rrho_local, ustar_local, z_local, &
            chu_local, chs_local, cmu_local, cms_local, sflh_local, qrlw_local, cldeff_local, kvh_local, kvm_local, &
            bflxs_local, bprod_local, sprod_local, wcap_local, leng_local, tke_local, turbtype_local)
       return
    end if

    call eddy_diff_caleddy_init_codon(int(ncol_local, c_int64_t), int(pcols_local, c_int64_t), int(pver_local, c_int64_t), &
         int(qrlzero_mode_local, c_int64_t), int(cldeff_mode_local, c_int64_t), int(tkes_mode_local, c_int64_t), &
         int(use_kvf_mode_local, c_int64_t), real(qmin_local, c_double), real(vk_local, c_double), c_loc(ql_local), &
         c_loc(qrlin_local), c_loc(cld_local), c_loc(kvf_local), c_loc(kvh_in_local), c_loc(kvm_in_local), c_loc(n2_local), &
         c_loc(s2_local), c_loc(shflx_local), c_loc(qflx_local), c_loc(rrho_local), c_loc(ustar_local), c_loc(z_local), &
         c_loc(chu_local), c_loc(chs_local), c_loc(cmu_local), c_loc(cms_local), c_loc(sflh_local), c_loc(qrlw_local), &
         c_loc(cldeff_local), c_loc(kvh_local), c_loc(kvm_local), c_loc(bflxs_local), c_loc(bprod_local), c_loc(sprod_local), &
         c_loc(wcap_local), c_loc(leng_local), c_loc(tke_local), c_loc(turbtype_local))

  end subroutine eddy_diff_caleddy_init

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_init_native(ncol_local, pcols_local, pver_local, qrlzero_mode_local, cldeff_mode_local, &
       tkes_mode_local, use_kvf_mode_local, qmin_local, vk_local, ql_local, qrlin_local, cld_local, kvf_local, kvh_in_local, &
       kvm_in_local, n2_local, s2_local, shflx_local, qflx_local, rrho_local, ustar_local, z_local, chu_local, chs_local, &
       cmu_local, cms_local, sflh_local, qrlw_local, cldeff_local, kvh_local, kvm_local, bflxs_local, bprod_local, sprod_local, &
       wcap_local, leng_local, tke_local, turbtype_local)

    implicit none

    integer, intent(in) :: ncol_local, pcols_local, pver_local
    integer, intent(in) :: qrlzero_mode_local, cldeff_mode_local, tkes_mode_local, use_kvf_mode_local
    real(r8), intent(in) :: qmin_local, vk_local
    real(r8), intent(in) :: ql_local(pcols_local,pver_local), qrlin_local(pcols_local,pver_local)
    real(r8), intent(in) :: cld_local(pcols_local,pver_local), kvf_local(pcols_local,pver_local+1)
    real(r8), intent(in) :: kvh_in_local(pcols_local,pver_local+1), kvm_in_local(pcols_local,pver_local+1)
    real(r8), intent(in) :: n2_local(pcols_local,pver_local), s2_local(pcols_local,pver_local)
    real(r8), intent(in) :: shflx_local(pcols_local), qflx_local(pcols_local), rrho_local(pcols_local)
    real(r8), intent(in) :: ustar_local(pcols_local), z_local(pcols_local,pver_local)
    real(r8), intent(in) :: chu_local(pcols_local,pver_local+1), chs_local(pcols_local,pver_local+1)
    real(r8), intent(in) :: cmu_local(pcols_local,pver_local+1), cms_local(pcols_local,pver_local+1)
    real(r8), intent(in) :: sflh_local(pcols_local,pver_local)
    real(r8), intent(inout) :: qrlw_local(pcols_local,pver_local), cldeff_local(pcols_local,pver_local)
    real(r8), intent(inout) :: kvh_local(pcols_local,pver_local+1), kvm_local(pcols_local,pver_local+1)
    real(r8), intent(inout) :: bflxs_local(pcols_local), bprod_local(pcols_local,pver_local+1)
    real(r8), intent(inout) :: sprod_local(pcols_local,pver_local+1), wcap_local(pcols_local,pver_local+1)
    real(r8), intent(inout) :: leng_local(pcols_local,pver_local+1), tke_local(pcols_local,pver_local+1)
    integer(i4), intent(inout) :: turbtype_local(pcols_local,pver_local+1)

    integer :: i, k
    real(r8) :: ch, cm

    if (qrlzero_mode_local /= 0) then
       qrlw_local(:,:) = 0._r8
    else
       qrlw_local(:ncol_local,:pver_local) = qrlin_local(:ncol_local,:pver_local)
    endif

    do k = 1, pver_local
       do i = 1, ncol_local
          if (cldeff_mode_local /= 0) then
             cldeff_local(i,k) = cld_local(i,k) * min(ql_local(i,k) / qmin_local, 1._r8)
          else
             cldeff_local(i,k) = cld_local(i,k)
          endif
       end do
    end do

    if (use_kvf_mode_local /= 0) then
       kvh_local(:,:) = kvf_local(:,:)
       kvm_local(:,:) = kvf_local(:,:)
    else
       kvh_local(:,:) = 0._r8
       kvm_local(:,:) = 0._r8
    endif

    wcap_local(:,:) = 0._r8
    leng_local(:,:) = 0._r8
    tke_local(:,:)  = 0._r8
    turbtype_local(:,:) = 0

    do k = 2, pver_local
       do i = 1, ncol_local
          bprod_local(i,k) = -kvh_in_local(i,k) * n2_local(i,k)
          sprod_local(i,k) =  kvm_in_local(i,k) * s2_local(i,k)
       end do
    end do

    do i = 1, ncol_local
       bprod_local(i,1) = 0._r8
       sprod_local(i,1) = 0._r8
       ch = chu_local(i,pver_local+1) * ( 1._r8 - sflh_local(i,pver_local) ) + chs_local(i,pver_local+1) * sflh_local(i,pver_local)
       cm = cmu_local(i,pver_local+1) * ( 1._r8 - sflh_local(i,pver_local) ) + cms_local(i,pver_local+1) * sflh_local(i,pver_local)
       bflxs_local(i) = ch * shflx_local(i) * rrho_local(i) + cm * qflx_local(i) * rrho_local(i)
       if (tkes_mode_local /= 0) then
          bprod_local(i,pver_local+1) = bflxs_local(i)
       else
          bprod_local(i,pver_local+1) = 0._r8
       endif
       sprod_local(i,pver_local+1) = (ustar_local(i)**3)/(vk_local*z_local(i,pver_local))
    end do

  end subroutine eddy_diff_caleddy_init_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_diaginit_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (caleddy_diaginit_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_CALEDDY_DIAGINIT_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_caleddy_diaginit_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_caleddy_diaginit_impl = .false.
    end if

    caleddy_diaginit_impl_selected = .true.

    if (masterproc) then
       if (use_native_caleddy_diaginit_impl) then
          write(iulog,*) 'eddy_diff_caleddy_diaginit implementation = native'
       else
          write(iulog,*) 'eddy_diff_caleddy_diaginit implementation = codon'
       end if
    end if

  end subroutine eddy_diff_caleddy_diaginit_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_diaginit(ncol_local, pcols_local, pver_local, ncvmax_local, went_local, wet_CL_local, &
       web_CL_local, jtbu_CL_local, jbbu_CL_local, evhc_CL_local, jt2slv_CL_local, n2ht_CL_local, n2hb_CL_local, lwp_CL_local, &
       opt_depth_CL_local, radinvfrac_CL_local, radf_CL_local, wstar_CL_local, wstar3fact_CL_local, ricl_local, ghcl_local, &
       shcl_local, smcl_local, ebrk_local, wbrk_local, lbrk_local, gh_a_local, sh_a_local, sm_a_local, ri_a_local, &
       sm_aw_local, ipbl_local, kpblh_local, wsed_CL_local)

    use iso_c_binding, only: c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: ncol_local, pcols_local, pver_local, ncvmax_local
    real(r8), target, intent(inout) :: went_local(pcols_local)
    real(r8), target, intent(inout) :: wet_CL_local(pcols_local,ncvmax_local), web_CL_local(pcols_local,ncvmax_local)
    real(r8), target, intent(inout) :: jtbu_CL_local(pcols_local,ncvmax_local), jbbu_CL_local(pcols_local,ncvmax_local)
    real(r8), target, intent(inout) :: evhc_CL_local(pcols_local,ncvmax_local), jt2slv_CL_local(pcols_local,ncvmax_local)
    real(r8), target, intent(inout) :: n2ht_CL_local(pcols_local,ncvmax_local), n2hb_CL_local(pcols_local,ncvmax_local)
    real(r8), target, intent(inout) :: lwp_CL_local(pcols_local,ncvmax_local), opt_depth_CL_local(pcols_local,ncvmax_local)
    real(r8), target, intent(inout) :: radinvfrac_CL_local(pcols_local,ncvmax_local), radf_CL_local(pcols_local,ncvmax_local)
    real(r8), target, intent(inout) :: wstar_CL_local(pcols_local,ncvmax_local), wstar3fact_CL_local(pcols_local,ncvmax_local)
    real(r8), target, intent(inout) :: ricl_local(pcols_local,ncvmax_local), ghcl_local(pcols_local,ncvmax_local)
    real(r8), target, intent(inout) :: shcl_local(pcols_local,ncvmax_local), smcl_local(pcols_local,ncvmax_local)
    real(r8), target, intent(inout) :: ebrk_local(pcols_local,ncvmax_local), wbrk_local(pcols_local,ncvmax_local)
    real(r8), target, intent(inout) :: lbrk_local(pcols_local,ncvmax_local), wsed_CL_local(pcols_local,ncvmax_local)
    real(r8), target, intent(inout) :: gh_a_local(pcols_local,pver_local+1), sh_a_local(pcols_local,pver_local+1)
    real(r8), target, intent(inout) :: sm_a_local(pcols_local,pver_local+1), ri_a_local(pcols_local,pver_local+1)
    real(r8), target, intent(inout) :: sm_aw_local(pcols_local,pver_local+1)
    integer(i4), target, intent(inout) :: ipbl_local(pcols_local), kpblh_local(pcols_local)

    interface
       subroutine eddy_diff_caleddy_diaginit_codon(ncol_c, pcols_c, pver_c, ncvmax_c, went_p, wet_CL_p, web_CL_p, jtbu_CL_p, &
            jbbu_CL_p, evhc_CL_p, jt2slv_CL_p, n2ht_CL_p, n2hb_CL_p, lwp_CL_p, opt_depth_CL_p, radinvfrac_CL_p, radf_CL_p, &
            wstar_CL_p, wstar3fact_CL_p, ricl_p, ghcl_p, shcl_p, smcl_p, ebrk_p, wbrk_p, lbrk_p, gh_a_p, sh_a_p, sm_a_p, &
            ri_a_p, sm_aw_p, ipbl_p, kpblh_p, wsed_CL_p) bind(c, name="eddy_diff_caleddy_diaginit_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, ncvmax_c
         type(c_ptr), value :: went_p, wet_CL_p, web_CL_p, jtbu_CL_p, jbbu_CL_p, evhc_CL_p, jt2slv_CL_p, n2ht_CL_p
         type(c_ptr), value :: n2hb_CL_p, lwp_CL_p, opt_depth_CL_p, radinvfrac_CL_p, radf_CL_p, wstar_CL_p
         type(c_ptr), value :: wstar3fact_CL_p, ricl_p, ghcl_p, shcl_p, smcl_p, ebrk_p, wbrk_p, lbrk_p, gh_a_p
         type(c_ptr), value :: sh_a_p, sm_a_p, ri_a_p, sm_aw_p, ipbl_p, kpblh_p, wsed_CL_p
       end subroutine eddy_diff_caleddy_diaginit_codon
    end interface

    call eddy_diff_caleddy_diaginit_select_impl()

    if (use_native_caleddy_diaginit_impl) then
       call eddy_diff_caleddy_diaginit_native(ncol_local, pcols_local, pver_local, ncvmax_local, went_local, wet_CL_local, &
            web_CL_local, jtbu_CL_local, jbbu_CL_local, evhc_CL_local, jt2slv_CL_local, n2ht_CL_local, n2hb_CL_local, &
            lwp_CL_local, opt_depth_CL_local, radinvfrac_CL_local, radf_CL_local, wstar_CL_local, wstar3fact_CL_local, &
            ricl_local, ghcl_local, shcl_local, smcl_local, ebrk_local, wbrk_local, lbrk_local, gh_a_local, sh_a_local, &
            sm_a_local, ri_a_local, sm_aw_local, ipbl_local, kpblh_local, wsed_CL_local)
       return
    end if

    call eddy_diff_caleddy_diaginit_codon(int(ncol_local, c_int64_t), int(pcols_local, c_int64_t), int(pver_local, c_int64_t), &
         int(ncvmax_local, c_int64_t), c_loc(went_local), c_loc(wet_CL_local), c_loc(web_CL_local), c_loc(jtbu_CL_local), &
         c_loc(jbbu_CL_local), c_loc(evhc_CL_local), c_loc(jt2slv_CL_local), c_loc(n2ht_CL_local), c_loc(n2hb_CL_local), &
         c_loc(lwp_CL_local), c_loc(opt_depth_CL_local), c_loc(radinvfrac_CL_local), c_loc(radf_CL_local), c_loc(wstar_CL_local), &
         c_loc(wstar3fact_CL_local), c_loc(ricl_local), c_loc(ghcl_local), c_loc(shcl_local), c_loc(smcl_local), &
         c_loc(ebrk_local), c_loc(wbrk_local), c_loc(lbrk_local), c_loc(gh_a_local), c_loc(sh_a_local), c_loc(sm_a_local), &
         c_loc(ri_a_local), c_loc(sm_aw_local), c_loc(ipbl_local), c_loc(kpblh_local), c_loc(wsed_CL_local))

  end subroutine eddy_diff_caleddy_diaginit

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_diaginit_native(ncol_local, pcols_local, pver_local, ncvmax_local, went_local, wet_CL_local, &
       web_CL_local, jtbu_CL_local, jbbu_CL_local, evhc_CL_local, jt2slv_CL_local, n2ht_CL_local, n2hb_CL_local, lwp_CL_local, &
       opt_depth_CL_local, radinvfrac_CL_local, radf_CL_local, wstar_CL_local, wstar3fact_CL_local, ricl_local, ghcl_local, &
       shcl_local, smcl_local, ebrk_local, wbrk_local, lbrk_local, gh_a_local, sh_a_local, sm_a_local, ri_a_local, &
       sm_aw_local, ipbl_local, kpblh_local, wsed_CL_local)

    implicit none

    integer, intent(in) :: ncol_local, pcols_local, pver_local, ncvmax_local
    real(r8), intent(inout) :: went_local(pcols_local)
    real(r8), intent(inout) :: wet_CL_local(pcols_local,ncvmax_local), web_CL_local(pcols_local,ncvmax_local)
    real(r8), intent(inout) :: jtbu_CL_local(pcols_local,ncvmax_local), jbbu_CL_local(pcols_local,ncvmax_local)
    real(r8), intent(inout) :: evhc_CL_local(pcols_local,ncvmax_local), jt2slv_CL_local(pcols_local,ncvmax_local)
    real(r8), intent(inout) :: n2ht_CL_local(pcols_local,ncvmax_local), n2hb_CL_local(pcols_local,ncvmax_local)
    real(r8), intent(inout) :: lwp_CL_local(pcols_local,ncvmax_local), opt_depth_CL_local(pcols_local,ncvmax_local)
    real(r8), intent(inout) :: radinvfrac_CL_local(pcols_local,ncvmax_local), radf_CL_local(pcols_local,ncvmax_local)
    real(r8), intent(inout) :: wstar_CL_local(pcols_local,ncvmax_local), wstar3fact_CL_local(pcols_local,ncvmax_local)
    real(r8), intent(inout) :: ricl_local(pcols_local,ncvmax_local), ghcl_local(pcols_local,ncvmax_local)
    real(r8), intent(inout) :: shcl_local(pcols_local,ncvmax_local), smcl_local(pcols_local,ncvmax_local)
    real(r8), intent(inout) :: ebrk_local(pcols_local,ncvmax_local), wbrk_local(pcols_local,ncvmax_local)
    real(r8), intent(inout) :: lbrk_local(pcols_local,ncvmax_local), wsed_CL_local(pcols_local,ncvmax_local)
    real(r8), intent(inout) :: gh_a_local(pcols_local,pver_local+1), sh_a_local(pcols_local,pver_local+1)
    real(r8), intent(inout) :: sm_a_local(pcols_local,pver_local+1), ri_a_local(pcols_local,pver_local+1)
    real(r8), intent(inout) :: sm_aw_local(pcols_local,pver_local+1)
    integer(i4), intent(inout) :: ipbl_local(pcols_local), kpblh_local(pcols_local)

    integer :: i

    do i = 1, ncol_local
       went_local(i)                  = 0._r8
       wet_CL_local(i,:ncvmax_local)        = 0._r8
       web_CL_local(i,:ncvmax_local)        = 0._r8
       jtbu_CL_local(i,:ncvmax_local)       = 0._r8
       jbbu_CL_local(i,:ncvmax_local)       = 0._r8
       evhc_CL_local(i,:ncvmax_local)       = 0._r8
       jt2slv_CL_local(i,:ncvmax_local)     = 0._r8
       n2ht_CL_local(i,:ncvmax_local)       = 0._r8
       n2hb_CL_local(i,:ncvmax_local)       = 0._r8
       lwp_CL_local(i,:ncvmax_local)        = 0._r8
       opt_depth_CL_local(i,:ncvmax_local)  = 0._r8
       radinvfrac_CL_local(i,:ncvmax_local) = 0._r8
       radf_CL_local(i,:ncvmax_local)       = 0._r8
       wstar_CL_local(i,:ncvmax_local)      = 0._r8
       wstar3fact_CL_local(i,:ncvmax_local) = 0._r8
       ricl_local(i,:ncvmax_local)          = 0._r8
       ghcl_local(i,:ncvmax_local)          = 0._r8
       shcl_local(i,:ncvmax_local)          = 0._r8
       smcl_local(i,:ncvmax_local)          = 0._r8
       ebrk_local(i,:ncvmax_local)          = 0._r8
       wbrk_local(i,:ncvmax_local)          = 0._r8
       lbrk_local(i,:ncvmax_local)          = 0._r8
       gh_a_local(i,:pver_local+1)          = 0._r8
       sh_a_local(i,:pver_local+1)          = 0._r8
       sm_a_local(i,:pver_local+1)          = 0._r8
       ri_a_local(i,:pver_local+1)          = 0._r8
       sm_aw_local(i,:pver_local+1)         = 0._r8
       ipbl_local(i)                        = 0
       kpblh_local(i)                       = pver_local
       wsed_CL_local(i,:ncvmax_local)       = 0._r8
    end do

  end subroutine eddy_diff_caleddy_diaginit_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_regime_diag_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (caleddy_regime_diag_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_CALEDDY_REGIME_DIAG_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_caleddy_regime_diag_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_caleddy_regime_diag_impl = .false.
    end if

    caleddy_regime_diag_impl_selected = .true.

    if (masterproc) then
       if (use_native_caleddy_regime_diag_impl) then
          write(iulog,*) 'eddy_diff_caleddy_regime_diag implementation = native'
       else
          write(iulog,*) 'eddy_diff_caleddy_regime_diag implementation = codon'
       end if
    end if

  end subroutine eddy_diff_caleddy_regime_diag_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_regime_diag(i_local, pcols_local, ncvmax_local, kbase_local, ktop_local, ncvfin_local, &
       kbase_diag_local, ktop_diag_local, ncvfin_diag_local)

    use iso_c_binding, only: c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: i_local, pcols_local, ncvmax_local
    integer(i4), target, intent(in) :: kbase_local(pcols_local,ncvmax_local), ktop_local(pcols_local,ncvmax_local)
    integer(i4), target, intent(in) :: ncvfin_local(pcols_local)
    real(r8), target, intent(inout) :: kbase_diag_local(pcols_local,ncvmax_local), ktop_diag_local(pcols_local,ncvmax_local)
    real(r8), target, intent(inout) :: ncvfin_diag_local(pcols_local)

    interface
       subroutine eddy_diff_caleddy_regime_diag_codon(i_c, pcols_c, ncvmax_c, kbase_p, ktop_p, ncvfin_p, kbase_diag_p, &
            ktop_diag_p, ncvfin_diag_p) bind(c, name="eddy_diff_caleddy_regime_diag_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: i_c, pcols_c, ncvmax_c
         type(c_ptr), value :: kbase_p, ktop_p, ncvfin_p, kbase_diag_p, ktop_diag_p, ncvfin_diag_p
       end subroutine eddy_diff_caleddy_regime_diag_codon
    end interface

    call eddy_diff_caleddy_regime_diag_select_impl()

    if (use_native_caleddy_regime_diag_impl) then
       call eddy_diff_caleddy_regime_diag_native(i_local, pcols_local, ncvmax_local, kbase_local, ktop_local, ncvfin_local, &
            kbase_diag_local, ktop_diag_local, ncvfin_diag_local)
       return
    end if

    call eddy_diff_caleddy_regime_diag_codon(int(i_local, c_int64_t), int(pcols_local, c_int64_t), int(ncvmax_local, c_int64_t), &
         c_loc(kbase_local), c_loc(ktop_local), c_loc(ncvfin_local), c_loc(kbase_diag_local), c_loc(ktop_diag_local), &
         c_loc(ncvfin_diag_local))

  end subroutine eddy_diff_caleddy_regime_diag

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_regime_diag_native(i_local, pcols_local, ncvmax_local, kbase_local, ktop_local, ncvfin_local, &
       kbase_diag_local, ktop_diag_local, ncvfin_diag_local)

    implicit none

    integer, intent(in) :: i_local, pcols_local, ncvmax_local
    integer(i4), intent(in) :: kbase_local(pcols_local,ncvmax_local), ktop_local(pcols_local,ncvmax_local)
    integer(i4), intent(in) :: ncvfin_local(pcols_local)
    real(r8), intent(inout) :: kbase_diag_local(pcols_local,ncvmax_local), ktop_diag_local(pcols_local,ncvmax_local)
    real(r8), intent(inout) :: ncvfin_diag_local(pcols_local)

    integer :: k

    do k = 1, ncvmax_local
       kbase_diag_local(i_local,k) = real(kbase_local(i_local,k),r8)
       ktop_diag_local(i_local,k)  = real(ktop_local(i_local,k),r8)
    end do
    ncvfin_diag_local(i_local) = real(ncvfin_local(i_local),r8)

  end subroutine eddy_diff_caleddy_regime_diag_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_stable_config_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (caleddy_stable_config_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_CALEDDY_STABLE_CONFIG_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_caleddy_stable_config_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_caleddy_stable_config_impl = .false.
    end if

    caleddy_stable_config_impl_selected = .true.

    if (masterproc) then
       if (use_native_caleddy_stable_config_impl) then
          write(iulog,*) 'eddy_diff_caleddy_stable_config implementation = native'
       else
          write(iulog,*) 'eddy_diff_caleddy_stable_config implementation = codon'
       end if
    end if

  end subroutine eddy_diff_caleddy_stable_config_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_stable_config(ricrit_local, b1_local, alph2_local, alph3_local, alph4_local, alph5_local, &
       alph4exs_local, ghmin_local)

    use iso_c_binding, only: c_double, c_loc, c_ptr

    implicit none

    real(r8), intent(in) :: ricrit_local, b1_local, alph2_local, alph3_local, alph4_local, alph5_local
    real(r8), target, intent(out) :: alph4exs_local, ghmin_local
    integer(i4), target :: stable_config_status_local

    interface
       subroutine eddy_diff_caleddy_stable_config_codon(ricrit_c, b1_c, alph2_c, alph3_c, alph4_c, alph5_c, alph4exs_p, &
            ghmin_p, status_p) bind(c, name="eddy_diff_caleddy_stable_config_codon")
         use iso_c_binding, only: c_double, c_ptr
         real(c_double), value :: ricrit_c, b1_c, alph2_c, alph3_c, alph4_c, alph5_c
         type(c_ptr), value :: alph4exs_p, ghmin_p, status_p
       end subroutine eddy_diff_caleddy_stable_config_codon
    end interface

    stable_config_status_local = 0_i4

    call eddy_diff_caleddy_stable_config_select_impl()

    if (use_native_caleddy_stable_config_impl) then
       call eddy_diff_caleddy_stable_config_native(ricrit_local, b1_local, alph2_local, alph3_local, alph4_local, &
            alph5_local, alph4exs_local, ghmin_local, stable_config_status_local)
    else
       call eddy_diff_caleddy_stable_config_codon(ricrit_local, b1_local, alph2_local, alph3_local, alph4_local, alph5_local, &
            c_loc(alph4exs_local), c_loc(ghmin_local), c_loc(stable_config_status_local))
    end if

    if (stable_config_status_local .ne. 0_i4) then
       write(iulog,*) 'Error : ricrit should be larger than 0.19 in UW PBL'
       call endrun('CALEDDY Error: ricrit should be larger than 0.19 in UW PBL')
    end if

  end subroutine eddy_diff_caleddy_stable_config

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_stable_config_native(ricrit_local, b1_local, alph2_local, alph3_local, alph4_local, alph5_local, &
       alph4exs_local, ghmin_local, stable_config_status_local)

    implicit none

    real(r8), intent(in) :: ricrit_local, b1_local, alph2_local, alph3_local, alph4_local, alph5_local
    real(r8), intent(out) :: alph4exs_local, ghmin_local
    integer(i4), intent(out) :: stable_config_status_local

    stable_config_status_local = 0_i4

    if( ricrit_local .eq. 0.19_r8 ) then
        alph4exs_local = alph4_local
        ghmin_local    = -3.5334_r8
    elseif( ricrit_local .gt. 0.19_r8 ) then
        alph4exs_local = -2._r8 * b1_local * alph2_local / ( alph3_local - 2._r8 * b1_local * alph5_local ) / ricrit_local
        ghmin_local    = -1.e10_r8
    else
        stable_config_status_local = 1_i4
    endif

  end subroutine eddy_diff_caleddy_stable_config_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_surface_tke_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (caleddy_surface_tke_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_CALEDDY_SURFACE_TKE_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_caleddy_surface_tke_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_caleddy_surface_tke_impl = .false.
    end if

    caleddy_surface_tke_impl_selected = .true.

    if (masterproc) then
       if (use_native_caleddy_surface_tke_impl) then
          write(iulog,*) 'eddy_diff_caleddy_surface_tke implementation = native'
       else
          write(iulog,*) 'eddy_diff_caleddy_surface_tke implementation = codon'
       end if
    end if

  end subroutine eddy_diff_caleddy_surface_tke_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_surface_tke(i_local, pcols_local, pver_local, b1_local, vk_local, tkemax_local, z_local, &
       bprod_local, sprod_local, tkes_local, tke_local, wcap_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: i_local, pcols_local, pver_local
    real(r8), intent(in) :: b1_local, vk_local, tkemax_local
    real(r8), target, intent(in) :: z_local(pcols_local,pver_local), bprod_local(pcols_local,pver_local+1)
    real(r8), target, intent(in) :: sprod_local(pcols_local,pver_local+1)
    real(r8), target, intent(inout) :: tkes_local(pcols_local), tke_local(pcols_local,pver_local+1)
    real(r8), target, intent(inout) :: wcap_local(pcols_local,pver_local+1)

    interface
       subroutine eddy_diff_caleddy_surface_tke_codon(i_c, pcols_c, pver_c, b1_c, vk_c, tkemax_c, z_p, bprod_p, sprod_p, &
            tkes_p, tke_p, wcap_p) bind(c, name="eddy_diff_caleddy_surface_tke_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: i_c, pcols_c, pver_c
         real(c_double), value :: b1_c, vk_c, tkemax_c
         type(c_ptr), value :: z_p, bprod_p, sprod_p, tkes_p, tke_p, wcap_p
       end subroutine eddy_diff_caleddy_surface_tke_codon
    end interface

    call eddy_diff_caleddy_surface_tke_select_impl()

    if (use_native_caleddy_surface_tke_impl) then
       call eddy_diff_caleddy_surface_tke_native(i_local, pcols_local, pver_local, b1_local, vk_local, tkemax_local, z_local, &
            bprod_local, sprod_local, tkes_local, tke_local, wcap_local)
       return
    end if

    call eddy_diff_caleddy_surface_tke_codon(int(i_local, c_int64_t), int(pcols_local, c_int64_t), int(pver_local, c_int64_t), &
         b1_local, vk_local, tkemax_local, c_loc(z_local), c_loc(bprod_local), c_loc(sprod_local), c_loc(tkes_local), &
         c_loc(tke_local), c_loc(wcap_local))

  end subroutine eddy_diff_caleddy_surface_tke

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_surface_tke_native(i_local, pcols_local, pver_local, b1_local, vk_local, tkemax_local, z_local, &
       bprod_local, sprod_local, tkes_local, tke_local, wcap_local)

    implicit none

    integer, intent(in) :: i_local, pcols_local, pver_local
    real(r8), intent(in) :: b1_local, vk_local, tkemax_local
    real(r8), intent(in) :: z_local(pcols_local,pver_local), bprod_local(pcols_local,pver_local+1)
    real(r8), intent(in) :: sprod_local(pcols_local,pver_local+1)
    real(r8), intent(inout) :: tkes_local(pcols_local), tke_local(pcols_local,pver_local+1)
    real(r8), intent(inout) :: wcap_local(pcols_local,pver_local+1)

    tkes_local(i_local) = max(b1_local*vk_local*z_local(i_local,pver_local)*(bprod_local(i_local,pver_local+1)+ &
         sprod_local(i_local,pver_local+1)), 1.e-7_r8)**(2._r8/3._r8)
    tkes_local(i_local) = min(tkes_local(i_local), tkemax_local)
    tke_local(i_local,pver_local+1)  = tkes_local(i_local)
    wcap_local(i_local,pver_local+1) = tkes_local(i_local)/b1_local

  end subroutine eddy_diff_caleddy_surface_tke_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_surface_energy_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (zisocl_surface_energy_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_ZISOCL_SURFACE_ENERGY_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_zisocl_surface_energy_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_zisocl_surface_energy_impl = .false.
    end if

    zisocl_surface_energy_impl_selected = .true.

    if (masterproc) then
       if (use_native_zisocl_surface_energy_impl) then
          write(iulog,*) 'eddy_diff_zisocl_surface_energy implementation = native'
       else
          write(iulog,*) 'eddy_diff_zisocl_surface_energy implementation = codon'
       end if
    end if

  end subroutine eddy_diff_zisocl_surface_energy_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_surface_energy(z_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local, gh_local, &
       sh_local, sm_local, dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, dw_surf_local)

    use iso_c_binding, only: c_double, c_loc, c_ptr

    implicit none

    real(r8), target, intent(in) :: z_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local
    real(r8), target, intent(out) :: gh_local, sh_local, sm_local
    real(r8), target, intent(out) :: dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, dw_surf_local

    interface
       subroutine eddy_diff_zisocl_surface_energy_codon(alph1_c, alph2_c, alph3_c, alph4_c, alph5_c, b1_c, vk_c, &
            z_surf_c, bprod_surf_c, sprod_surf_c, tkes_surf_c, gh_p, sh_p, sm_p, dlint_surf_p, dl2n2_surf_p, &
            dl2s2_surf_p, dw_surf_p) bind(c, name="eddy_diff_zisocl_surface_energy_codon")
         use iso_c_binding, only: c_double, c_ptr
         real(c_double), value :: alph1_c, alph2_c, alph3_c, alph4_c, alph5_c, b1_c, vk_c
         real(c_double), value :: z_surf_c, bprod_surf_c, sprod_surf_c, tkes_surf_c
         type(c_ptr), value :: gh_p, sh_p, sm_p, dlint_surf_p, dl2n2_surf_p, dl2s2_surf_p, dw_surf_p
       end subroutine eddy_diff_zisocl_surface_energy_codon
    end interface

    call eddy_diff_zisocl_surface_energy_select_impl()

    if (use_native_zisocl_surface_energy_impl) then
       call eddy_diff_zisocl_surface_energy_native(z_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local, gh_local, &
            sh_local, sm_local, dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, dw_surf_local)
       return
    end if

    call eddy_diff_zisocl_surface_energy_codon(alph1, alph2, alph3, alph4, alph5, b1, vk, z_surf_local, bprod_surf_local, &
         sprod_surf_local, tkes_surf_local, c_loc(gh_local), c_loc(sh_local), c_loc(sm_local), c_loc(dlint_surf_local), &
         c_loc(dl2n2_surf_local), c_loc(dl2s2_surf_local), c_loc(dw_surf_local))

  end subroutine eddy_diff_zisocl_surface_energy

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_surface_energy_native(z_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local, gh_local, &
       sh_local, sm_local, dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, dw_surf_local)

    implicit none

    real(r8), intent(in) :: z_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local
    real(r8), intent(out) :: gh_local, sh_local, sm_local
    real(r8), intent(out) :: dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, dw_surf_local

    real(r8) :: gg_local

    gg_local = 0.5_r8*vk*z_surf_local*bprod_surf_local/(tkes_surf_local**(3._r8/2._r8))
    gh_local = gg_local/(alph5-gg_local*alph3)
    gh_local = min(max(gh_local,-3.5334_r8),0.0233_r8)
    sh_local = alph5/(1._r8+alph3*gh_local)
    sm_local = (alph1 + alph2*gh_local)/(1._r8+alph3*gh_local)/(1._r8+alph4*gh_local)
    dlint_surf_local = z_surf_local
    dl2n2_surf_local = -vk*(z_surf_local**2._r8)*bprod_surf_local/(sh_local*sqrt(tkes_surf_local))
    dl2s2_surf_local =  vk*(z_surf_local**2._r8)*sprod_surf_local/(sm_local*sqrt(tkes_surf_local))
    dw_surf_local = (tkes_surf_local/b1)*z_surf_local

  end subroutine eddy_diff_zisocl_surface_energy_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_surface_state_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (zisocl_surface_state_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_ZISOCL_SURFACE_STATE_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_zisocl_surface_state_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_zisocl_surface_state_impl = .false.
    end if

    zisocl_surface_state_impl_selected = .true.

    if (masterproc) then
       if (use_native_zisocl_surface_state_impl) then
          write(iulog,*) 'eddy_diff_zisocl_surface_state implementation = native'
       else
          write(iulog,*) 'eddy_diff_zisocl_surface_state implementation = codon'
       end if
    end if

  end subroutine eddy_diff_zisocl_surface_state_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_surface_state(kb_is_surface_local, use_dw_surf_local, zi_top_local, zi_base_local, z_surf_local, &
       bflxs_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local, lbulk_local, gh_local, sh_local, sm_local, &
       dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, dw_surf_local, lint_local, l2n2_local, l2s2_local, wint_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: kb_is_surface_local, use_dw_surf_local
    real(r8), intent(in) :: zi_top_local, zi_base_local, z_surf_local
    real(r8), intent(in) :: bflxs_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local
    real(r8), target, intent(out) :: lbulk_local, dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, dw_surf_local
    real(r8), target, intent(out) :: lint_local, l2n2_local, l2s2_local, wint_local
    real(r8), target, intent(inout) :: gh_local, sh_local, sm_local

    interface
       subroutine eddy_diff_zisocl_surface_state_codon(alph1_c, alph2_c, alph3_c, alph4_c, alph5_c, b1_c, vk_c, lbulk_max_c, &
            kb_is_surface_c, use_dw_surf_c, zi_top_c, zi_base_c, z_surf_c, bflxs_surf_c, bprod_surf_c, sprod_surf_c, &
            tkes_surf_c, lbulk_p, gh_p, sh_p, sm_p, dlint_surf_p, dl2n2_surf_p, dl2s2_surf_p, dw_surf_p, lint_p, l2n2_p, &
            l2s2_p, wint_p) bind(c, name="eddy_diff_zisocl_surface_state_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         real(c_double), value :: alph1_c, alph2_c, alph3_c, alph4_c, alph5_c, b1_c, vk_c, lbulk_max_c
         integer(c_int64_t), value :: kb_is_surface_c, use_dw_surf_c
         real(c_double), value :: zi_top_c, zi_base_c, z_surf_c, bflxs_surf_c, bprod_surf_c, sprod_surf_c, tkes_surf_c
         type(c_ptr), value :: lbulk_p, gh_p, sh_p, sm_p, dlint_surf_p, dl2n2_surf_p, dl2s2_surf_p, dw_surf_p
         type(c_ptr), value :: lint_p, l2n2_p, l2s2_p, wint_p
       end subroutine eddy_diff_zisocl_surface_state_codon
    end interface

    call eddy_diff_zisocl_surface_state_select_impl()

    if (use_native_zisocl_surface_state_impl) then
       call eddy_diff_zisocl_surface_state_native(kb_is_surface_local, use_dw_surf_local, zi_top_local, zi_base_local, z_surf_local, &
            bflxs_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local, lbulk_local, gh_local, sh_local, sm_local, &
            dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, dw_surf_local, lint_local, l2n2_local, l2s2_local, wint_local)
       return
    end if

    call eddy_diff_zisocl_surface_state_codon(alph1, alph2, alph3, alph4, alph5, b1, vk, lbulk_max, int(kb_is_surface_local, c_int64_t), &
         int(use_dw_surf_local, c_int64_t), zi_top_local, zi_base_local, z_surf_local, bflxs_surf_local, bprod_surf_local, &
         sprod_surf_local, tkes_surf_local, c_loc(lbulk_local), c_loc(gh_local), c_loc(sh_local), c_loc(sm_local), &
         c_loc(dlint_surf_local), c_loc(dl2n2_surf_local), c_loc(dl2s2_surf_local), c_loc(dw_surf_local), c_loc(lint_local), &
         c_loc(l2n2_local), c_loc(l2s2_local), c_loc(wint_local))

  end subroutine eddy_diff_zisocl_surface_state

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_surface_state_native(kb_is_surface_local, use_dw_surf_local, zi_top_local, zi_base_local, z_surf_local, &
       bflxs_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local, lbulk_local, gh_local, sh_local, sm_local, &
       dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, dw_surf_local, lint_local, l2n2_local, l2s2_local, wint_local)

    implicit none

    integer, intent(in) :: kb_is_surface_local, use_dw_surf_local
    real(r8), intent(in) :: zi_top_local, zi_base_local, z_surf_local
    real(r8), intent(in) :: bflxs_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local
    real(r8), intent(out) :: lbulk_local, dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, dw_surf_local
    real(r8), intent(out) :: lint_local, l2n2_local, l2s2_local, wint_local
    real(r8), intent(inout) :: gh_local, sh_local, sm_local

    lbulk_local = zi_top_local - zi_base_local
    lbulk_local = min(lbulk_local, lbulk_max)
    dlint_surf_local = 0._r8
    dl2n2_surf_local = 0._r8
    dl2s2_surf_local = 0._r8
    dw_surf_local = 0._r8

    if (kb_is_surface_local .ne. 0) then
       if (bflxs_surf_local .gt. 0._r8) then
          call eddy_diff_zisocl_surface_energy(z_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local, gh_local, &
               sh_local, sm_local, dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, dw_surf_local)
       else
          lbulk_local = zi_top_local - z_surf_local
          lbulk_local = min(lbulk_local, lbulk_max)
       end if
    end if

    lint_local = dlint_surf_local
    l2n2_local = dl2n2_surf_local
    l2s2_local = dl2s2_surf_local
    wint_local = dw_surf_local

    if (use_dw_surf_local .ne. 0) then
       l2n2_local = 0._r8
       l2s2_local = 0._r8
    else
       wint_local = 0._r8
    end if

  end subroutine eddy_diff_zisocl_surface_state_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_surface_extend_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (zisocl_surface_extend_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_ZISOCL_SURFACE_EXTEND_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_zisocl_surface_extend_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_zisocl_surface_extend_impl = .false.
    end if

    zisocl_surface_extend_impl_selected = .true.

    if (masterproc) then
       if (use_native_zisocl_surface_extend_impl) then
          write(iulog,*) 'eddy_diff_zisocl_surface_extend implementation = native'
       else
          write(iulog,*) 'eddy_diff_zisocl_surface_extend implementation = codon'
       end if
    end if

  end subroutine eddy_diff_zisocl_surface_extend_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_surface_extend(bflxs_surf_local, z_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local, &
       sh_local, dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, dw_surf_local, lint_local, l2n2_local, l2s2_local, &
       wint_local)

    use iso_c_binding, only: c_double, c_loc, c_ptr

    implicit none

    real(r8), target, intent(in) :: bflxs_surf_local, z_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local, sh_local
    real(r8), target, intent(out) :: dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, dw_surf_local
    real(r8), target, intent(inout) :: lint_local, l2n2_local, l2s2_local, wint_local

    interface
       subroutine eddy_diff_zisocl_surface_extend_codon(alph1_c, alph2_c, alph3_c, alph4_c, alph5_c, b1_c, vk_c, tkemax_c, &
            bflxs_surf_c, z_surf_c, bprod_surf_c, sprod_surf_c, tkes_surf_c, sh_c, dlint_surf_p, dl2n2_surf_p, dl2s2_surf_p, &
            dw_surf_p, lint_p, l2n2_p, l2s2_p, wint_p) bind(c, name="eddy_diff_zisocl_surface_extend_codon")
         use iso_c_binding, only: c_double, c_ptr
         real(c_double), value :: alph1_c, alph2_c, alph3_c, alph4_c, alph5_c, b1_c, vk_c, tkemax_c
         real(c_double), value :: bflxs_surf_c, z_surf_c, bprod_surf_c, sprod_surf_c, tkes_surf_c, sh_c
         type(c_ptr), value :: dlint_surf_p, dl2n2_surf_p, dl2s2_surf_p, dw_surf_p, lint_p, l2n2_p, l2s2_p, wint_p
       end subroutine eddy_diff_zisocl_surface_extend_codon
    end interface

    call eddy_diff_zisocl_surface_extend_select_impl()

    if (use_native_zisocl_surface_extend_impl) then
       call eddy_diff_zisocl_surface_extend_native(bflxs_surf_local, z_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local, &
            sh_local, dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, dw_surf_local, lint_local, l2n2_local, l2s2_local, &
            wint_local)
       return
    end if

    call eddy_diff_zisocl_surface_extend_codon(alph1, alph2, alph3, alph4, alph5, b1, vk, tkemax, bflxs_surf_local, z_surf_local, &
         bprod_surf_local, sprod_surf_local, tkes_surf_local, sh_local, c_loc(dlint_surf_local), c_loc(dl2n2_surf_local), &
         c_loc(dl2s2_surf_local), c_loc(dw_surf_local), c_loc(lint_local), c_loc(l2n2_local), c_loc(l2s2_local), c_loc(wint_local))

  end subroutine eddy_diff_zisocl_surface_extend

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_surface_extend_native(bflxs_surf_local, z_surf_local, bprod_surf_local, sprod_surf_local, &
       tkes_surf_local, sh_local, dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, dw_surf_local, lint_local, &
       l2n2_local, l2s2_local, wint_local)

    implicit none

    real(r8), intent(in) :: bflxs_surf_local, z_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local, sh_local
    real(r8), intent(out) :: dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, dw_surf_local
    real(r8), intent(inout) :: lint_local, l2n2_local, l2s2_local, wint_local

    real(r8) :: gh_surf_local, sh_surf_local, sm_surf_local

    if (bflxs_surf_local .gt. 0._r8) then
       call eddy_diff_zisocl_surface_energy(z_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local, gh_surf_local, &
            sh_surf_local, sm_surf_local, dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, dw_surf_local)
    else
       dlint_surf_local = 0._r8
       dl2n2_surf_local = 0._r8
       dl2s2_surf_local = 0._r8
       dw_surf_local = 0._r8
    end if

    lint_local = lint_local + dlint_surf_local
    l2n2_local = l2n2_local + dl2n2_surf_local
    l2n2_local = -min(-l2n2_local, tkemax*lint_local/(b1*sh_local))
    l2s2_local = l2s2_local + dl2s2_surf_local
    wint_local = wint_local + dw_surf_local

  end subroutine eddy_diff_zisocl_surface_extend_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_sbcl_state_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (zisocl_sbcl_state_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_ZISOCL_SBCL_STATE_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_zisocl_sbcl_state_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_zisocl_sbcl_state_impl = .false.
    end if

    zisocl_sbcl_state_impl_selected = .true.

    if (masterproc) then
       if (use_native_zisocl_sbcl_state_impl) then
          write(iulog,*) 'eddy_diff_zisocl_sbcl_state implementation = native'
       else
          write(iulog,*) 'eddy_diff_zisocl_sbcl_state implementation = codon'
       end if
    end if

  end subroutine eddy_diff_zisocl_sbcl_state_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_sbcl_state(choice_tkes_ebprod_local, sh_local, dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, &
       dw_surf_local, lint_local, l2n2_local, l2s2_local, wint_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: choice_tkes_ebprod_local
    real(r8), target, intent(in) :: sh_local, dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, dw_surf_local
    real(r8), target, intent(out) :: lint_local, l2n2_local, l2s2_local, wint_local

    interface
       subroutine eddy_diff_zisocl_sbcl_state_codon(choice_tkes_ebprod_c, sh_c, dlint_surf_c, dl2n2_surf_c, dl2s2_surf_c, &
            dw_surf_c, lint_p, l2n2_p, l2s2_p, wint_p) bind(c, name="eddy_diff_zisocl_sbcl_state_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: choice_tkes_ebprod_c
         real(c_double), value :: sh_c, dlint_surf_c, dl2n2_surf_c, dl2s2_surf_c, dw_surf_c
         type(c_ptr), value :: lint_p, l2n2_p, l2s2_p, wint_p
       end subroutine eddy_diff_zisocl_sbcl_state_codon
    end interface

    call eddy_diff_zisocl_sbcl_state_select_impl()

    if (use_native_zisocl_sbcl_state_impl) then
       call eddy_diff_zisocl_sbcl_state_native(choice_tkes_ebprod_local, sh_local, dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, &
            dw_surf_local, lint_local, l2n2_local, l2s2_local, wint_local)
       return
    end if

    call eddy_diff_zisocl_sbcl_state_codon(int(choice_tkes_ebprod_local, c_int64_t), sh_local, dlint_surf_local, dl2n2_surf_local, &
         dl2s2_surf_local, dw_surf_local, c_loc(lint_local), c_loc(l2n2_local), c_loc(l2s2_local), c_loc(wint_local))

  end subroutine eddy_diff_zisocl_sbcl_state

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_sbcl_state_native(choice_tkes_ebprod_local, sh_local, dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, &
       dw_surf_local, lint_local, l2n2_local, l2s2_local, wint_local)

    implicit none

    integer, intent(in) :: choice_tkes_ebprod_local
    real(r8), intent(in) :: sh_local, dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, dw_surf_local
    real(r8), intent(out) :: lint_local, l2n2_local, l2s2_local, wint_local

    lint_local = dlint_surf_local
    l2n2_local = dl2n2_surf_local
    l2s2_local = dl2s2_surf_local
    wint_local = dw_surf_local

    if (choice_tkes_ebprod_local .ne. 0) then
       l2n2_local = -wint_local / sh_local
    end if

  end subroutine eddy_diff_zisocl_sbcl_state_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_initial_state_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (zisocl_initial_state_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_ZISOCL_INITIAL_STATE_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_zisocl_initial_state_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_zisocl_initial_state_impl = .false.
    end if

    zisocl_initial_state_impl_selected = .true.

    if (masterproc) then
       if (use_native_zisocl_initial_state_impl) then
          write(iulog,*) 'eddy_diff_zisocl_initial_state implementation = native'
       else
          write(iulog,*) 'eddy_diff_zisocl_initial_state implementation = codon'
       end if
    end if

  end subroutine eddy_diff_zisocl_initial_state_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_initial_state(i_local, kt_local, kb_local, pcols_local, pver_local, use_dw_surf_local, &
       choice_tkes_ebprod_local, z_surf_local, bflxs_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local, z_local, &
       zi_local, n2_local, s2_local, leng_max_local, lbulk_local, gh_local, sh_local, sm_local, dlint_surf_local, &
       dl2n2_surf_local, dl2s2_surf_local, dw_surf_local, lint_local, l2n2_local, l2s2_local, wint_local, ricll_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: i_local, kt_local, kb_local, pcols_local, pver_local, use_dw_surf_local, choice_tkes_ebprod_local
    real(r8), intent(in) :: z_surf_local, bflxs_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local
    real(r8), target, intent(in) :: z_local(pcols_local,pver_local), zi_local(pcols_local,pver_local+1)
    real(r8), target, intent(in) :: n2_local(pcols_local,pver_local), s2_local(pcols_local,pver_local)
    real(r8), target, intent(in) :: leng_max_local(pver_local)
    real(r8), target, intent(inout) :: lbulk_local, gh_local, sh_local, sm_local
    real(r8), target, intent(inout) :: dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, dw_surf_local
    real(r8), target, intent(inout) :: lint_local, l2n2_local, l2s2_local, wint_local, ricll_local

    integer :: kb_is_surface_local, tunl_mode_local, leng_mode_local

    interface
       subroutine eddy_diff_zisocl_initial_state_codon(i_c, kt_c, kb_c, pcols_c, pver_c, kb_is_surface_c, use_dw_surf_c, &
            choice_tkes_ebprod_c, tunl_mode_c, leng_mode_c, alph1_c, alph2_c, alph3_c, alph4_c, alph5_c, b1_c, vk_c, ntzero_c, &
            ricrit_c, lbulk_max_c, tunl_c, ctunl_c, cleng_c, tkemax_c, z_surf_c, bflxs_surf_c, bprod_surf_c, sprod_surf_c, &
            tkes_surf_c, z_p, zi_p, n2_p, s2_p, leng_max_p, lbulk_p, gh_p, sh_p, sm_p, dlint_surf_p, dl2n2_surf_p, &
            dl2s2_surf_p, dw_surf_p, lint_p, l2n2_p, l2s2_p, wint_p, ricll_p) bind(c, &
            name="eddy_diff_zisocl_initial_state_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: i_c, kt_c, kb_c, pcols_c, pver_c, kb_is_surface_c, use_dw_surf_c
         integer(c_int64_t), value :: choice_tkes_ebprod_c, tunl_mode_c, leng_mode_c
         real(c_double), value :: alph1_c, alph2_c, alph3_c, alph4_c, alph5_c, b1_c, vk_c, ntzero_c, ricrit_c
         real(c_double), value :: lbulk_max_c, tunl_c, ctunl_c, cleng_c, tkemax_c
         real(c_double), value :: z_surf_c, bflxs_surf_c, bprod_surf_c, sprod_surf_c, tkes_surf_c
         type(c_ptr), value :: z_p, zi_p, n2_p, s2_p, leng_max_p, lbulk_p, gh_p, sh_p, sm_p, dlint_surf_p, dl2n2_surf_p
         type(c_ptr), value :: dl2s2_surf_p, dw_surf_p, lint_p, l2n2_p, l2s2_p, wint_p, ricll_p
       end subroutine eddy_diff_zisocl_initial_state_codon
    end interface

    kb_is_surface_local = 0
    if( kb_local .eq. pver_local + 1 ) kb_is_surface_local = 1

    tunl_mode_local = 0
    if( choice_tunl .eq. 'rampcl' ) then
        tunl_mode_local = 1
    elseif( choice_tunl .eq. 'rampsl' ) then
        tunl_mode_local = 2
    end if

    leng_mode_local = 1
    if( choice_leng .eq. 'origin' ) then
        leng_mode_local = 0
    end if

    call eddy_diff_zisocl_initial_state_select_impl()

    if (use_native_zisocl_initial_state_impl) then
       call eddy_diff_zisocl_initial_state_native(i_local, kt_local, kb_local, pcols_local, pver_local, use_dw_surf_local, &
            choice_tkes_ebprod_local, z_surf_local, bflxs_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local, &
            z_local, zi_local, n2_local, s2_local, leng_max_local, lbulk_local, gh_local, sh_local, sm_local, dlint_surf_local, &
            dl2n2_surf_local, dl2s2_surf_local, dw_surf_local, lint_local, l2n2_local, l2s2_local, wint_local, ricll_local)
       return
    end if

    call eddy_diff_zisocl_initial_state_codon(int(i_local, c_int64_t), int(kt_local, c_int64_t), int(kb_local, c_int64_t), &
         int(pcols_local, c_int64_t), int(pver_local, c_int64_t), int(kb_is_surface_local, c_int64_t), &
         int(use_dw_surf_local, c_int64_t), int(choice_tkes_ebprod_local, c_int64_t), int(tunl_mode_local, c_int64_t), &
         int(leng_mode_local, c_int64_t), alph1, alph2, alph3, alph4, alph5, b1, vk, ntzero, ricrit, lbulk_max, tunl, ctunl, &
         cleng, tkemax, z_surf_local, bflxs_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local, c_loc(z_local), &
         c_loc(zi_local), c_loc(n2_local), c_loc(s2_local), c_loc(leng_max_local), c_loc(lbulk_local), c_loc(gh_local), &
         c_loc(sh_local), c_loc(sm_local), c_loc(dlint_surf_local), c_loc(dl2n2_surf_local), c_loc(dl2s2_surf_local), &
         c_loc(dw_surf_local), c_loc(lint_local), c_loc(l2n2_local), c_loc(l2s2_local), c_loc(wint_local), c_loc(ricll_local))

  end subroutine eddy_diff_zisocl_initial_state

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_initial_state_native(i_local, kt_local, kb_local, pcols_local, pver_local, use_dw_surf_local, &
       choice_tkes_ebprod_local, z_surf_local, bflxs_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local, z_local, &
       zi_local, n2_local, s2_local, leng_max_local, lbulk_local, gh_local, sh_local, sm_local, dlint_surf_local, &
       dl2n2_surf_local, dl2s2_surf_local, dw_surf_local, lint_local, l2n2_local, l2s2_local, wint_local, ricll_local)

    implicit none

    integer, intent(in) :: i_local, kt_local, kb_local, pcols_local, pver_local, use_dw_surf_local, choice_tkes_ebprod_local
    real(r8), intent(in) :: z_surf_local, bflxs_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local
    real(r8), intent(in) :: z_local(pcols_local,pver_local), zi_local(pcols_local,pver_local+1)
    real(r8), intent(in) :: n2_local(pcols_local,pver_local), s2_local(pcols_local,pver_local)
    real(r8), intent(in) :: leng_max_local(pver_local)
    real(r8), intent(inout) :: lbulk_local, gh_local, sh_local, sm_local
    real(r8), intent(inout) :: dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, dw_surf_local
    real(r8), intent(inout) :: lint_local, l2n2_local, l2s2_local, wint_local, ricll_local

    integer :: kb_is_surface_local

    kb_is_surface_local = 0
    if( kb_local .eq. pver_local + 1 ) kb_is_surface_local = 1

    call eddy_diff_zisocl_surface_state(kb_is_surface_local, use_dw_surf_local, zi_local(i_local,kt_local), zi_local(i_local,kb_local), &
         z_surf_local, bflxs_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local, lbulk_local, gh_local, sh_local, &
         sm_local, dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, dw_surf_local, lint_local, l2n2_local, l2s2_local, &
         wint_local)
    if( kb_local .eq. pver_local+1 .and. bflxs_surf_local .gt. 0._r8 ) then
        ricll_local = min(-(sm_local/sh_local)*(bprod_surf_local/sprod_surf_local),ricrit)
    end if

    if( kt_local .lt. kb_local - 1 ) then
        call eddy_diff_zisocl_non_sbcl_state_native(i_local, kt_local, kb_local, pcols_local, pver_local, lbulk_local, z_local, &
             zi_local, n2_local, s2_local, leng_max_local, lint_local, l2n2_local, l2s2_local, wint_local, ricll_local, gh_local, &
             sh_local, sm_local)
    else
        call eddy_diff_zisocl_sbcl_state(choice_tkes_ebprod_local, sh_local, dlint_surf_local, dl2n2_surf_local, &
             dl2s2_surf_local, dw_surf_local, lint_local, l2n2_local, l2s2_local, wint_local)
    end if

    l2n2_local = -min(-l2n2_local, tkemax*lint_local/(b1*sh_local))
    l2s2_local =  min( l2s2_local, tkemax*lint_local/(b1*sm_local))

  end subroutine eddy_diff_zisocl_initial_state_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_extended_state_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (zisocl_extended_state_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_ZISOCL_EXTENDED_STATE_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_zisocl_extended_state_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_zisocl_extended_state_impl = .false.
    end if

    zisocl_extended_state_impl_selected = .true.

    if (masterproc) then
       if (use_native_zisocl_extended_state_impl) then
          write(iulog,*) 'eddy_diff_zisocl_extended_state implementation = native'
       else
          write(iulog,*) 'eddy_diff_zisocl_extended_state implementation = codon'
       end if
    end if

  end subroutine eddy_diff_zisocl_extended_state_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_extended_state(i_local, kt_local, kb_local, pcols_local, pver_local, use_dw_surf_local, &
       zi_top_local, zi_base_local, z_surf_local, bflxs_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local, &
       z_local, zi_local, n2_local, s2_local, leng_max_local, gh_local, sh_local, sm_local, lint_local, wint_local, ricll_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: i_local, kt_local, kb_local, pcols_local, pver_local, use_dw_surf_local
    real(r8), intent(in) :: zi_top_local, zi_base_local, z_surf_local
    real(r8), intent(in) :: bflxs_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local
    real(r8), target, intent(in) :: z_local(pcols_local,pver_local), zi_local(pcols_local,pver_local+1)
    real(r8), target, intent(in) :: n2_local(pcols_local,pver_local), s2_local(pcols_local,pver_local)
    real(r8), target, intent(in) :: leng_max_local(pver_local)
    real(r8), target, intent(out) :: gh_local, sh_local, sm_local, lint_local, wint_local, ricll_local

    integer :: tunl_mode_local, leng_mode_local

    interface
       subroutine eddy_diff_zisocl_extended_state_codon(i_c, kt_c, kb_c, pcols_c, pver_c, use_dw_surf_c, tunl_mode_c, leng_mode_c, &
            alph1_c, alph2_c, alph3_c, alph4_c, alph5_c, b1_c, vk_c, ntzero_c, ricrit_c, lbulk_max_c, tunl_c, ctunl_c, cleng_c, &
            tkemax_c, zi_top_c, zi_base_c, z_surf_c, bflxs_surf_c, bprod_surf_c, sprod_surf_c, tkes_surf_c, z_p, zi_p, n2_p, s2_p, &
            leng_max_p, gh_p, sh_p, sm_p, lint_p, wint_p, ricll_p) bind(c, name="eddy_diff_zisocl_extended_state_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: i_c, kt_c, kb_c, pcols_c, pver_c, use_dw_surf_c, tunl_mode_c, leng_mode_c
         real(c_double), value :: alph1_c, alph2_c, alph3_c, alph4_c, alph5_c, b1_c, vk_c, ntzero_c, ricrit_c, lbulk_max_c
         real(c_double), value :: tunl_c, ctunl_c, cleng_c, tkemax_c
         real(c_double), value :: zi_top_c, zi_base_c, z_surf_c, bflxs_surf_c, bprod_surf_c, sprod_surf_c, tkes_surf_c
         type(c_ptr), value :: z_p, zi_p, n2_p, s2_p, leng_max_p, gh_p, sh_p, sm_p, lint_p, wint_p, ricll_p
       end subroutine eddy_diff_zisocl_extended_state_codon
    end interface

    tunl_mode_local = 0
    if( choice_tunl .eq. 'rampcl' ) then
        tunl_mode_local = 1
    elseif( choice_tunl .eq. 'rampsl' ) then
        tunl_mode_local = 2
    end if

    leng_mode_local = 1
    if( choice_leng .eq. 'origin' ) then
        leng_mode_local = 0
    end if

    call eddy_diff_zisocl_extended_state_select_impl()

    if (use_native_zisocl_extended_state_impl) then
       call eddy_diff_zisocl_extended_state_native(i_local, kt_local, kb_local, pcols_local, pver_local, use_dw_surf_local, &
            zi_top_local, zi_base_local, z_surf_local, bflxs_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local, &
            z_local, zi_local, n2_local, s2_local, leng_max_local, gh_local, sh_local, sm_local, lint_local, wint_local, ricll_local)
       return
    end if

    call eddy_diff_zisocl_extended_state_codon(int(i_local, c_int64_t), int(kt_local, c_int64_t), int(kb_local, c_int64_t), &
         int(pcols_local, c_int64_t), int(pver_local, c_int64_t), int(use_dw_surf_local, c_int64_t), int(tunl_mode_local, c_int64_t), &
         int(leng_mode_local, c_int64_t), alph1, alph2, alph3, alph4, alph5, b1, vk, ntzero, ricrit, lbulk_max, tunl, ctunl, cleng, &
         tkemax, zi_top_local, zi_base_local, z_surf_local, bflxs_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local, &
         c_loc(z_local), c_loc(zi_local), c_loc(n2_local), c_loc(s2_local), c_loc(leng_max_local), c_loc(gh_local), c_loc(sh_local), &
         c_loc(sm_local), c_loc(lint_local), c_loc(wint_local), c_loc(ricll_local))

  end subroutine eddy_diff_zisocl_extended_state

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_extended_state_native(i_local, kt_local, kb_local, pcols_local, pver_local, use_dw_surf_local, &
       zi_top_local, zi_base_local, z_surf_local, bflxs_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local, &
       z_local, zi_local, n2_local, s2_local, leng_max_local, gh_local, sh_local, sm_local, lint_local, wint_local, ricll_local)

    implicit none

    integer, intent(in) :: i_local, kt_local, kb_local, pcols_local, pver_local, use_dw_surf_local
    real(r8), intent(in) :: zi_top_local, zi_base_local, z_surf_local
    real(r8), intent(in) :: bflxs_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local
    real(r8), intent(in) :: z_local(pcols_local,pver_local), zi_local(pcols_local,pver_local+1)
    real(r8), intent(in) :: n2_local(pcols_local,pver_local), s2_local(pcols_local,pver_local)
    real(r8), intent(in) :: leng_max_local(pver_local)
    real(r8), intent(out) :: gh_local, sh_local, sm_local, lint_local, wint_local, ricll_local

    integer :: k_local, kb_is_surface_mode_local
    real(r8) :: lbulk_local
    real(r8) :: dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, dw_surf_local
    real(r8) :: dzinc_local, dl2n2_local, dl2s2_local, l2n2_local, l2s2_local

    kb_is_surface_mode_local = 0
    if( kb_local .eq. pver_local + 1 ) kb_is_surface_mode_local = 1

    call eddy_diff_zisocl_surface_state(kb_is_surface_mode_local, use_dw_surf_local, zi_top_local, zi_base_local, z_surf_local, &
         bflxs_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local, lbulk_local, gh_local, sh_local, sm_local, &
         dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, dw_surf_local, lint_local, l2n2_local, l2s2_local, wint_local)

    do k_local = kt_local + 1, kb_local - 1
       call eddy_diff_zisocl_layer_energy_native(i_local, k_local, pcols_local, pver_local, lbulk_local, z_local, zi_local, n2_local, &
            s2_local, leng_max_local, dzinc_local, dl2n2_local, dl2s2_local)
       lint_local = lint_local + dzinc_local
       l2n2_local = l2n2_local + dl2n2_local
       l2s2_local = l2s2_local + dl2s2_local
    end do

    call eddy_diff_zisocl_stability_native(l2n2_local, l2s2_local, ricll_local, gh_local, sh_local, sm_local)
    wint_local = max( wint_local - sh_local*l2n2_local + sm_local*l2s2_local, 0.01_r8 )

  end subroutine eddy_diff_zisocl_extended_state_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_non_sbcl_state_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (zisocl_non_sbcl_state_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_ZISOCL_NON_SBCL_STATE_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_zisocl_non_sbcl_state_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_zisocl_non_sbcl_state_impl = .false.
    end if

    zisocl_non_sbcl_state_impl_selected = .true.

    if (masterproc) then
       if (use_native_zisocl_non_sbcl_state_impl) then
          write(iulog,*) 'eddy_diff_zisocl_non_sbcl_state implementation = native'
       else
          write(iulog,*) 'eddy_diff_zisocl_non_sbcl_state implementation = codon'
       end if
    end if

  end subroutine eddy_diff_zisocl_non_sbcl_state_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_non_sbcl_state(i_local, kt_local, kb_local, pcols_local, pver_local, lbulk_local, z_local, zi_local, &
       n2_local, s2_local, leng_max_local, lint_local, l2n2_local, l2s2_local, wint_local, ricll_local, gh_local, sh_local, sm_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: i_local, kt_local, kb_local, pcols_local, pver_local
    real(r8), intent(in) :: lbulk_local
    real(r8), target, intent(in) :: z_local(pcols_local,pver_local), zi_local(pcols_local,pver_local+1)
    real(r8), target, intent(in) :: n2_local(pcols_local,pver_local), s2_local(pcols_local,pver_local)
    real(r8), target, intent(in) :: leng_max_local(pver_local)
    real(r8), target, intent(inout) :: lint_local, l2n2_local, l2s2_local, wint_local
    real(r8), target, intent(out) :: ricll_local, gh_local, sh_local, sm_local

    integer :: tunl_mode_local, leng_mode_local

    interface
       subroutine eddy_diff_zisocl_non_sbcl_state_codon(i_c, kt_c, kb_c, pcols_c, pver_c, tunl_mode_c, leng_mode_c, alph1_c, &
            alph2_c, alph3_c, alph4_c, alph5_c, b1_c, vk_c, ntzero_c, ricrit_c, tunl_c, ctunl_c, cleng_c, lbulk_c, z_p, zi_p, &
            n2_p, s2_p, leng_max_p, lint_p, l2n2_p, l2s2_p, wint_p, ricll_p, gh_p, sh_p, sm_p) &
            bind(c, name="eddy_diff_zisocl_non_sbcl_state_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: i_c, kt_c, kb_c, pcols_c, pver_c, tunl_mode_c, leng_mode_c
         real(c_double), value :: alph1_c, alph2_c, alph3_c, alph4_c, alph5_c, b1_c, vk_c, ntzero_c, ricrit_c
         real(c_double), value :: tunl_c, ctunl_c, cleng_c, lbulk_c
         type(c_ptr), value :: z_p, zi_p, n2_p, s2_p, leng_max_p, lint_p, l2n2_p, l2s2_p, wint_p, ricll_p, gh_p, sh_p, sm_p
       end subroutine eddy_diff_zisocl_non_sbcl_state_codon
    end interface

    tunl_mode_local = 0
    if( choice_tunl .eq. 'rampcl' ) then
        tunl_mode_local = 1
    elseif( choice_tunl .eq. 'rampsl' ) then
        tunl_mode_local = 2
    end if

    leng_mode_local = 1
    if( choice_leng .eq. 'origin' ) then
        leng_mode_local = 0
    end if

    call eddy_diff_zisocl_non_sbcl_state_select_impl()

    if (use_native_zisocl_non_sbcl_state_impl) then
       call eddy_diff_zisocl_non_sbcl_state_native(i_local, kt_local, kb_local, pcols_local, pver_local, lbulk_local, z_local, zi_local, &
            n2_local, s2_local, leng_max_local, lint_local, l2n2_local, l2s2_local, wint_local, ricll_local, gh_local, sh_local, &
            sm_local)
       return
    end if

    call eddy_diff_zisocl_non_sbcl_state_codon(int(i_local, c_int64_t), int(kt_local, c_int64_t), int(kb_local, c_int64_t), &
         int(pcols_local, c_int64_t), int(pver_local, c_int64_t), int(tunl_mode_local, c_int64_t), int(leng_mode_local, c_int64_t), &
         alph1, alph2, alph3, alph4, alph5, b1, vk, ntzero, ricrit, tunl, ctunl, cleng, lbulk_local, c_loc(z_local), c_loc(zi_local), &
         c_loc(n2_local), c_loc(s2_local), c_loc(leng_max_local), c_loc(lint_local), c_loc(l2n2_local), c_loc(l2s2_local), &
         c_loc(wint_local), c_loc(ricll_local), c_loc(gh_local), c_loc(sh_local), c_loc(sm_local))

  end subroutine eddy_diff_zisocl_non_sbcl_state

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_non_sbcl_state_native(i_local, kt_local, kb_local, pcols_local, pver_local, lbulk_local, z_local, zi_local, &
       n2_local, s2_local, leng_max_local, lint_local, l2n2_local, l2s2_local, wint_local, ricll_local, gh_local, sh_local, sm_local)

    implicit none

    integer, intent(in) :: i_local, kt_local, kb_local, pcols_local, pver_local
    real(r8), intent(in) :: lbulk_local
    real(r8), intent(in) :: z_local(pcols_local,pver_local), zi_local(pcols_local,pver_local+1)
    real(r8), intent(in) :: n2_local(pcols_local,pver_local), s2_local(pcols_local,pver_local)
    real(r8), intent(in) :: leng_max_local(pver_local)
    real(r8), intent(inout) :: lint_local, l2n2_local, l2s2_local, wint_local
    real(r8), intent(out) :: ricll_local, gh_local, sh_local, sm_local

    integer :: k_local
    real(r8) :: dzinc_local, dl2n2_local, dl2s2_local

    do k_local = kb_local - 1, kt_local + 1, -1
       call eddy_diff_zisocl_layer_energy_native(i_local, k_local, pcols_local, pver_local, lbulk_local, z_local, zi_local, n2_local, &
            s2_local, leng_max_local, dzinc_local, dl2n2_local, dl2s2_local)
       l2n2_local = l2n2_local + dl2n2_local
       l2s2_local = l2s2_local + dl2s2_local
       lint_local = lint_local + dzinc_local
    end do

    call eddy_diff_zisocl_stability_native(l2n2_local, l2s2_local, ricll_local, gh_local, sh_local, sm_local)
    wint_local = wint_local - sh_local*l2n2_local + sm_local*l2s2_local

  end subroutine eddy_diff_zisocl_non_sbcl_state_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_upward_state_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (zisocl_upward_state_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_ZISOCL_UPWARD_STATE_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_zisocl_upward_state_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_zisocl_upward_state_impl = .false.
    end if

    zisocl_upward_state_impl_selected = .true.

    if (masterproc) then
       if (use_native_zisocl_upward_state_impl) then
          write(iulog,*) 'eddy_diff_zisocl_upward_state implementation = native'
       else
          write(iulog,*) 'eddy_diff_zisocl_upward_state implementation = codon'
       end if
    end if

  end subroutine eddy_diff_zisocl_upward_state_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_upward_state(i_local, kt_local, pcols_local, pver_local, ntop_turb_local, ncv_local, &
       lbulk_local, sh_local, sm_local, z_local, zi_local, n2_local, s2_local, leng_max_local, ncvfin_local, kbase_local, &
       ktop_local, lint_local, l2n2_local, l2s2_local, wint_local, status_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: i_local, pcols_local, pver_local, ntop_turb_local, ncv_local
    integer, intent(inout) :: kt_local
    real(r8), intent(in) :: lbulk_local, sh_local, sm_local
    real(r8), target, intent(in) :: z_local(pcols_local,pver_local), zi_local(pcols_local,pver_local+1)
    real(r8), target, intent(in) :: n2_local(pcols_local,pver_local), s2_local(pcols_local,pver_local)
    real(r8), target, intent(in) :: leng_max_local(pver_local)
    integer(i4), target, intent(inout) :: ncvfin_local(pcols_local), kbase_local(pcols_local,ncvmax)
    integer(i4), target, intent(inout) :: ktop_local(pcols_local,ncvmax)
    real(r8), target, intent(inout) :: lint_local, l2n2_local, l2s2_local, wint_local
    integer, intent(out) :: status_local

    integer :: tunl_mode_local, leng_mode_local
    integer(c_int64_t), target :: kt_codon, status_codon

    interface
       subroutine eddy_diff_zisocl_upward_state_codon(i_c, pcols_c, pver_c, ntop_turb_c, ncv_c, tunl_mode_c, leng_mode_c, &
            tunl_c, ctunl_c, cleng_c, vk_c, rinc_c, tkemax_c, b1_c, lbulk_c, sh_c, sm_c, z_p, zi_p, n2_p, s2_p, leng_max_p, &
            ncvfin_p, kbase_p, ktop_p, kt_p, lint_p, l2n2_p, l2s2_p, wint_p, status_p) bind(c, &
            name="eddy_diff_zisocl_upward_state_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: i_c, pcols_c, pver_c, ntop_turb_c, ncv_c, tunl_mode_c, leng_mode_c
         real(c_double), value :: tunl_c, ctunl_c, cleng_c, vk_c, rinc_c, tkemax_c, b1_c, lbulk_c, sh_c, sm_c
         type(c_ptr), value :: z_p, zi_p, n2_p, s2_p, leng_max_p, ncvfin_p, kbase_p, ktop_p, kt_p, lint_p, l2n2_p, l2s2_p
         type(c_ptr), value :: wint_p, status_p
       end subroutine eddy_diff_zisocl_upward_state_codon
    end interface

    tunl_mode_local = 0
    if( choice_tunl .eq. 'rampcl' ) then
        tunl_mode_local = 1
    elseif( choice_tunl .eq. 'rampsl' ) then
        tunl_mode_local = 2
    end if

    leng_mode_local = 1
    if( choice_leng .eq. 'origin' ) then
        leng_mode_local = 0
    end if

    call eddy_diff_zisocl_upward_state_select_impl()

    status_local = 0

    if (use_native_zisocl_upward_state_impl) then
       call eddy_diff_zisocl_upward_state_native(i_local, kt_local, pcols_local, pver_local, ntop_turb_local, ncv_local, &
            lbulk_local, sh_local, sm_local, z_local, zi_local, n2_local, s2_local, leng_max_local, ncvfin_local, kbase_local, &
            ktop_local, lint_local, l2n2_local, l2s2_local, wint_local, status_local)
       return
    end if

    kt_codon = int(kt_local, c_int64_t)
    status_codon = 0_c_int64_t

    call eddy_diff_zisocl_upward_state_codon(int(i_local, c_int64_t), int(pcols_local, c_int64_t), int(pver_local, c_int64_t), &
         int(ntop_turb_local, c_int64_t), int(ncv_local, c_int64_t), int(tunl_mode_local, c_int64_t), &
         int(leng_mode_local, c_int64_t), tunl, ctunl, cleng, vk, rinc, tkemax, b1, lbulk_local, sh_local, sm_local, &
         c_loc(z_local), c_loc(zi_local), c_loc(n2_local), c_loc(s2_local), c_loc(leng_max_local), c_loc(ncvfin_local), &
         c_loc(kbase_local), c_loc(ktop_local), c_loc(kt_codon), c_loc(lint_local), c_loc(l2n2_local), c_loc(l2s2_local), &
         c_loc(wint_local), c_loc(status_codon))

    kt_local = int(kt_codon)
    status_local = int(status_codon)

  end subroutine eddy_diff_zisocl_upward_state

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_upward_state_native(i_local, kt_local, pcols_local, pver_local, ntop_turb_local, ncv_local, &
       lbulk_local, sh_local, sm_local, z_local, zi_local, n2_local, s2_local, leng_max_local, ncvfin_local, kbase_local, &
       ktop_local, lint_local, l2n2_local, l2s2_local, wint_local, status_local)

    implicit none

    integer, intent(in) :: i_local, pcols_local, pver_local, ntop_turb_local, ncv_local
    integer, intent(inout) :: kt_local
    real(r8), intent(in) :: lbulk_local, sh_local, sm_local
    real(r8), intent(in) :: z_local(pcols_local,pver_local), zi_local(pcols_local,pver_local+1)
    real(r8), intent(in) :: n2_local(pcols_local,pver_local), s2_local(pcols_local,pver_local)
    real(r8), intent(in) :: leng_max_local(pver_local)
    integer(i4), intent(inout) :: ncvfin_local(pcols_local), kbase_local(pcols_local,ncvmax)
    integer(i4), intent(inout) :: ktop_local(pcols_local,ncvmax)
    real(r8), intent(inout) :: lint_local, l2n2_local, l2s2_local, wint_local
    integer, intent(out) :: status_local

    integer :: cntu_local, k_local, incv_local, ktinc_local
    real(r8) :: dzinc_local, dl2n2_local, dl2s2_local, dwinc_local

    status_local = 0
    cntu_local = 0

    call eddy_diff_zisocl_interface_energy_native(i_local, kt_local, pcols_local, pver_local, lbulk_local, sh_local, sm_local, &
         z_local, zi_local, n2_local, s2_local, leng_max_local, dzinc_local, dl2n2_local, dl2s2_local, dwinc_local)

    do while ( -dl2n2_local .gt. (-rinc*l2n2_local/(1._r8-rinc)) .and. (kt_local > ntop_turb_local + 2 .or. &
         z_local(i_local,kt_local) < 50000._r8) )

       lint_local = lint_local + dzinc_local
       l2n2_local = l2n2_local + dl2n2_local
       l2n2_local = -min(-l2n2_local, tkemax*lint_local/(b1*sh_local))
       l2s2_local = l2s2_local + dl2s2_local
       wint_local = wint_local + dwinc_local

       kt_local = kt_local - 1
       if( kt_local .eq. ntop_turb_local ) then
           status_local = 1
           return
       end if

       ktinc_local = kbase_local(i_local,ncv_local+cntu_local+1) - 1

       if( kt_local .eq. ktinc_local ) then

           do k_local = kbase_local(i_local,ncv_local+cntu_local+1) - 1, ktop_local(i_local,ncv_local+cntu_local+1) + 1, -1

              call eddy_diff_zisocl_interface_energy_native(i_local, k_local, pcols_local, pver_local, lbulk_local, sh_local, &
                   sm_local, z_local, zi_local, n2_local, s2_local, leng_max_local, dzinc_local, dl2n2_local, dl2s2_local, &
                   dwinc_local)

              lint_local = lint_local + dzinc_local
              l2n2_local = l2n2_local + dl2n2_local
              l2n2_local = -min(-l2n2_local, tkemax*lint_local/(b1*sh_local))
              l2s2_local = l2s2_local + dl2s2_local
              wint_local = wint_local + dwinc_local

           end do

           kt_local = ktop_local(i_local,ncv_local+cntu_local+1)
           ncvfin_local(i_local) = ncvfin_local(i_local) - 1
           cntu_local = cntu_local + 1

       end if

       call eddy_diff_zisocl_interface_energy_native(i_local, kt_local, pcols_local, pver_local, lbulk_local, sh_local, sm_local, &
            z_local, zi_local, n2_local, s2_local, leng_max_local, dzinc_local, dl2n2_local, dl2s2_local, dwinc_local)

    end do

    if( cntu_local .gt. 0 ) then
        do incv_local = 1, ncvfin_local(i_local) - ncv_local
           kbase_local(i_local,ncv_local+incv_local) = kbase_local(i_local,ncv_local+cntu_local+incv_local)
           ktop_local(i_local,ncv_local+incv_local)  = ktop_local(i_local,ncv_local+cntu_local+incv_local)
        end do
    end if

  end subroutine eddy_diff_zisocl_upward_state_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_downward_state_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (zisocl_downward_state_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_ZISOCL_DOWNWARD_STATE_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_zisocl_downward_state_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_zisocl_downward_state_impl = .false.
    end if

    zisocl_downward_state_impl_selected = .true.

    if (masterproc) then
       if (use_native_zisocl_downward_state_impl) then
          write(iulog,*) 'eddy_diff_zisocl_downward_state implementation = native'
       else
          write(iulog,*) 'eddy_diff_zisocl_downward_state implementation = codon'
       end if
    end if

  end subroutine eddy_diff_zisocl_downward_state_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_downward_state(i_local, kb_local, ncv_local, ncvinit_local, pcols_local, pver_local, &
       lbulk_local, sh_local, sm_local, z_surf_local, bflxs_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local, &
       z_local, zi_local, n2_local, s2_local, leng_max_local, dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, &
       dw_surf_local, ncvfin_local, kbase_local, ktop_local, lint_local, l2n2_local, l2s2_local, wint_local, status_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: i_local, ncvinit_local, pcols_local, pver_local
    integer, intent(inout) :: kb_local, ncv_local
    real(r8), intent(in) :: lbulk_local, sh_local, sm_local
    real(r8), intent(in) :: z_surf_local, bflxs_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local
    real(r8), target, intent(in) :: z_local(pcols_local,pver_local), zi_local(pcols_local,pver_local+1)
    real(r8), target, intent(in) :: n2_local(pcols_local,pver_local), s2_local(pcols_local,pver_local)
    real(r8), target, intent(in) :: leng_max_local(pver_local)
    real(r8), target, intent(inout) :: dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, dw_surf_local
    integer(i4), target, intent(inout) :: ncvfin_local(pcols_local), kbase_local(pcols_local,ncvmax)
    integer(i4), target, intent(inout) :: ktop_local(pcols_local,ncvmax)
    real(r8), target, intent(inout) :: lint_local, l2n2_local, l2s2_local, wint_local
    integer, intent(out) :: status_local

    integer :: tunl_mode_local, leng_mode_local
    integer(c_int64_t), target :: kb_codon, ncv_codon, status_codon

    interface
       subroutine eddy_diff_zisocl_downward_state_codon(i_c, pcols_c, pver_c, ncvinit_c, tunl_mode_c, leng_mode_c, &
            tunl_c, ctunl_c, cleng_c, vk_c, rinc_c, tkemax_c, b1_c, alph1_c, alph2_c, alph3_c, alph4_c, alph5_c, &
            lbulk_c, sh_c, sm_c, z_surf_c, bflxs_surf_c, bprod_surf_c, sprod_surf_c, tkes_surf_c, z_p, zi_p, n2_p, s2_p, &
            leng_max_p, dlint_surf_p, dl2n2_surf_p, dl2s2_surf_p, dw_surf_p, ncvfin_p, kbase_p, ktop_p, kb_p, ncv_p, lint_p, &
            l2n2_p, l2s2_p, wint_p, status_p) bind(c, name="eddy_diff_zisocl_downward_state_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: i_c, pcols_c, pver_c, ncvinit_c, tunl_mode_c, leng_mode_c
         real(c_double), value :: tunl_c, ctunl_c, cleng_c, vk_c, rinc_c, tkemax_c, b1_c
         real(c_double), value :: alph1_c, alph2_c, alph3_c, alph4_c, alph5_c
         real(c_double), value :: lbulk_c, sh_c, sm_c, z_surf_c, bflxs_surf_c, bprod_surf_c, sprod_surf_c, tkes_surf_c
         type(c_ptr), value :: z_p, zi_p, n2_p, s2_p, leng_max_p, dlint_surf_p, dl2n2_surf_p, dl2s2_surf_p, dw_surf_p
         type(c_ptr), value :: ncvfin_p, kbase_p, ktop_p, kb_p, ncv_p, lint_p, l2n2_p, l2s2_p, wint_p, status_p
       end subroutine eddy_diff_zisocl_downward_state_codon
    end interface

    tunl_mode_local = 0
    if( choice_tunl .eq. 'rampcl' ) then
        tunl_mode_local = 1
    elseif( choice_tunl .eq. 'rampsl' ) then
        tunl_mode_local = 2
    end if

    leng_mode_local = 1
    if( choice_leng .eq. 'origin' ) then
        leng_mode_local = 0
    end if

    call eddy_diff_zisocl_downward_state_select_impl()

    status_local = 0

    if (use_native_zisocl_downward_state_impl) then
       call eddy_diff_zisocl_downward_state_native(i_local, kb_local, ncv_local, ncvinit_local, pcols_local, pver_local, &
            lbulk_local, sh_local, sm_local, z_surf_local, bflxs_surf_local, bprod_surf_local, sprod_surf_local, &
            tkes_surf_local, z_local, zi_local, n2_local, s2_local, leng_max_local, dlint_surf_local, dl2n2_surf_local, &
            dl2s2_surf_local, dw_surf_local, ncvfin_local, kbase_local, ktop_local, lint_local, l2n2_local, l2s2_local, &
            wint_local, status_local)
       return
    end if

    kb_codon = int(kb_local, c_int64_t)
    ncv_codon = int(ncv_local, c_int64_t)
    status_codon = 0_c_int64_t

    call eddy_diff_zisocl_downward_state_codon(int(i_local, c_int64_t), int(pcols_local, c_int64_t), int(pver_local, c_int64_t), &
         int(ncvinit_local, c_int64_t), int(tunl_mode_local, c_int64_t), int(leng_mode_local, c_int64_t), tunl, ctunl, cleng, &
         vk, rinc, tkemax, b1, alph1, alph2, alph3, alph4, alph5, lbulk_local, sh_local, sm_local, z_surf_local, &
         bflxs_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local, c_loc(z_local), c_loc(zi_local), c_loc(n2_local), &
         c_loc(s2_local), c_loc(leng_max_local), c_loc(dlint_surf_local), c_loc(dl2n2_surf_local), c_loc(dl2s2_surf_local), &
         c_loc(dw_surf_local), c_loc(ncvfin_local), c_loc(kbase_local), c_loc(ktop_local), c_loc(kb_codon), c_loc(ncv_codon), &
         c_loc(lint_local), c_loc(l2n2_local), c_loc(l2s2_local), c_loc(wint_local), c_loc(status_codon))

    kb_local = int(kb_codon)
    ncv_local = int(ncv_codon)
    status_local = int(status_codon)

  end subroutine eddy_diff_zisocl_downward_state

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_downward_state_native(i_local, kb_local, ncv_local, ncvinit_local, pcols_local, pver_local, &
       lbulk_local, sh_local, sm_local, z_surf_local, bflxs_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local, &
       z_local, zi_local, n2_local, s2_local, leng_max_local, dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, &
       dw_surf_local, ncvfin_local, kbase_local, ktop_local, lint_local, l2n2_local, l2s2_local, wint_local, status_local)

    implicit none

    integer, intent(in) :: i_local, ncvinit_local, pcols_local, pver_local
    integer, intent(inout) :: kb_local, ncv_local
    real(r8), intent(in) :: lbulk_local, sh_local, sm_local
    real(r8), intent(in) :: z_surf_local, bflxs_surf_local, bprod_surf_local, sprod_surf_local, tkes_surf_local
    real(r8), intent(in) :: z_local(pcols_local,pver_local), zi_local(pcols_local,pver_local+1)
    real(r8), intent(in) :: n2_local(pcols_local,pver_local), s2_local(pcols_local,pver_local)
    real(r8), intent(in) :: leng_max_local(pver_local)
    real(r8), intent(inout) :: dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, dw_surf_local
    integer(i4), intent(inout) :: ncvfin_local(pcols_local), kbase_local(pcols_local,ncvmax)
    integer(i4), intent(inout) :: ktop_local(pcols_local,ncvmax)
    real(r8), intent(inout) :: lint_local, l2n2_local, l2s2_local, wint_local
    integer, intent(out) :: status_local

    integer :: cntd_local, k_local, incv_local, kbinc_local
    real(r8) :: dzinc_local, dl2n2_local, dl2s2_local, dwinc_local

    status_local = 0
    cntd_local = 0

    if( kb_local .eq. pver_local + 1 ) return

    call eddy_diff_zisocl_interface_energy_native(i_local, kb_local, pcols_local, pver_local, lbulk_local, sh_local, sm_local, &
         z_local, zi_local, n2_local, s2_local, leng_max_local, dzinc_local, dl2n2_local, dl2s2_local, dwinc_local)

    do while( ( -dl2n2_local .gt. (-rinc*l2n2_local/(1._r8-rinc)) ) .and. (kb_local .ne. pver_local+1) )

       lint_local = lint_local + dzinc_local
       l2n2_local = l2n2_local + dl2n2_local
       l2n2_local = -min(-l2n2_local, tkemax*lint_local/(b1*sh_local))
       l2s2_local = l2s2_local + dl2s2_local
       wint_local = wint_local + dwinc_local

       kb_local = kb_local + 1

       kbinc_local = 0
       if( ncv_local .gt. 1 ) kbinc_local = ktop_local(i_local,ncv_local-1) + 1
       if( kb_local .eq. kbinc_local ) then

           do k_local = ktop_local(i_local,ncv_local-1) + 1, kbase_local(i_local,ncv_local-1) - 1

              call eddy_diff_zisocl_interface_energy_native(i_local, k_local, pcols_local, pver_local, lbulk_local, sh_local, &
                   sm_local, z_local, zi_local, n2_local, s2_local, leng_max_local, dzinc_local, dl2n2_local, dl2s2_local, &
                   dwinc_local)

              lint_local = lint_local + dzinc_local
              l2n2_local = l2n2_local + dl2n2_local
              l2n2_local = -min(-l2n2_local, tkemax*lint_local/(b1*sh_local))
              l2s2_local = l2s2_local + dl2s2_local
              wint_local = wint_local + dwinc_local

           end do

           kb_local = kbase_local(i_local,ncv_local-1)
           ncv_local = ncv_local - 1
           ncvfin_local(i_local) = ncvfin_local(i_local) - 1
           cntd_local = cntd_local + 1

       end if

       if( kb_local .eq. pver_local + 1 ) then

           call eddy_diff_zisocl_surface_extend(bflxs_surf_local, z_surf_local, bprod_surf_local, sprod_surf_local, &
                tkes_surf_local, sh_local, dlint_surf_local, dl2n2_surf_local, dl2s2_surf_local, dw_surf_local, lint_local, &
                l2n2_local, l2s2_local, wint_local)

       else

           call eddy_diff_zisocl_interface_energy_native(i_local, kb_local, pcols_local, pver_local, lbulk_local, sh_local, sm_local, &
                z_local, zi_local, n2_local, s2_local, leng_max_local, dzinc_local, dl2n2_local, dl2s2_local, dwinc_local)

       end if

    end do

    if( (kb_local.eq.pver_local+1) .and. (ncv_local.ne.1) ) then
        status_local = 1
        return
    end if

    if( cntd_local .gt. 0 ) then
        do incv_local = 1, ncvfin_local(i_local) - ncv_local
           kbase_local(i_local,ncv_local+incv_local) = kbase_local(i_local,ncvinit_local+incv_local)
           ktop_local(i_local,ncv_local+incv_local)  = ktop_local(i_local,ncvinit_local+incv_local)
        end do
    end if

  end subroutine eddy_diff_zisocl_downward_state_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_stability_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (zisocl_stability_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_ZISOCL_STABILITY_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_zisocl_stability_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_zisocl_stability_impl = .false.
    end if

    zisocl_stability_impl_selected = .true.

    if (masterproc) then
       if (use_native_zisocl_stability_impl) then
          write(iulog,*) 'eddy_diff_zisocl_stability implementation = native'
       else
          write(iulog,*) 'eddy_diff_zisocl_stability implementation = codon'
       end if
    end if

  end subroutine eddy_diff_zisocl_stability_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_stability(l2n2_local, l2s2_local, ricll_local, gh_local, sh_local, sm_local)

    use iso_c_binding, only: c_double, c_loc, c_ptr

    implicit none

    real(r8), target, intent(in) :: l2n2_local, l2s2_local
    real(r8), target, intent(out) :: ricll_local, gh_local, sh_local, sm_local

    interface
       subroutine eddy_diff_zisocl_stability_codon(alph1_c, alph2_c, alph3_c, alph4_c, alph5_c, b1_c, ntzero_c, ricrit_c, &
            l2n2_c, l2s2_c, ricll_p, gh_p, sh_p, sm_p) bind(c, name="eddy_diff_zisocl_stability_codon")
         use iso_c_binding, only: c_double, c_ptr
         real(c_double), value :: alph1_c, alph2_c, alph3_c, alph4_c, alph5_c, b1_c, ntzero_c, ricrit_c
         real(c_double), value :: l2n2_c, l2s2_c
         type(c_ptr), value :: ricll_p, gh_p, sh_p, sm_p
       end subroutine eddy_diff_zisocl_stability_codon
    end interface

    call eddy_diff_zisocl_stability_select_impl()

    if (use_native_zisocl_stability_impl) then
       call eddy_diff_zisocl_stability_native(l2n2_local, l2s2_local, ricll_local, gh_local, sh_local, sm_local)
       return
    end if

    call eddy_diff_zisocl_stability_codon(alph1, alph2, alph3, alph4, alph5, b1, ntzero, ricrit, l2n2_local, l2s2_local, &
         c_loc(ricll_local), c_loc(gh_local), c_loc(sh_local), c_loc(sm_local))

  end subroutine eddy_diff_zisocl_stability

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_stability_native(l2n2_local, l2s2_local, ricll_local, gh_local, sh_local, sm_local)

    implicit none

    real(r8), intent(in) :: l2n2_local, l2s2_local
    real(r8), intent(out) :: ricll_local, gh_local, sh_local, sm_local

    real(r8) :: trma_local, trmb_local, trmc_local, det_local

    ricll_local = min(l2n2_local/max(l2s2_local,ntzero),ricrit)
    trma_local = alph3*alph4*ricll_local+2._r8*b1*(alph2-alph4*alph5*ricll_local)
    trmb_local = ricll_local*(alph3+alph4)+2._r8*b1*(-alph5*ricll_local+alph1)
    trmc_local = ricll_local
    det_local = max(trmb_local*trmb_local-4._r8*trma_local*trmc_local,0._r8)
    gh_local = (-trmb_local + sqrt(det_local))/2._r8/trma_local
    gh_local = min(max(gh_local,-3.5334_r8),0.0233_r8)
    sh_local = alph5/(1._r8+alph3*gh_local)
    sm_local = (alph1 + alph2*gh_local)/(1._r8+alph3*gh_local)/(1._r8+alph4*gh_local)

  end subroutine eddy_diff_zisocl_stability_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_layer_energy_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (zisocl_layer_energy_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_ZISOCL_LAYER_ENERGY_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_zisocl_layer_energy_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_zisocl_layer_energy_impl = .false.
    end if

    zisocl_layer_energy_impl_selected = .true.

    if (masterproc) then
       if (use_native_zisocl_layer_energy_impl) then
          write(iulog,*) 'eddy_diff_zisocl_layer_energy implementation = native'
       else
          write(iulog,*) 'eddy_diff_zisocl_layer_energy implementation = codon'
       end if
    end if

  end subroutine eddy_diff_zisocl_layer_energy_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_layer_energy(i_local, k_local, pcols_local, pver_local, lbulk_local, z_local, zi_local, &
       n2_local, s2_local, leng_max_local, dzinc_local, dl2n2_local, dl2s2_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: i_local, k_local, pcols_local, pver_local
    real(r8), intent(in) :: lbulk_local
    real(r8), target, intent(in) :: z_local(pcols_local,pver_local), zi_local(pcols_local,pver_local+1)
    real(r8), target, intent(in) :: n2_local(pcols_local,pver_local), s2_local(pcols_local,pver_local)
    real(r8), target, intent(in) :: leng_max_local(pver_local)
    real(r8), target, intent(out) :: dzinc_local, dl2n2_local, dl2s2_local

    integer :: tunl_mode_local, leng_mode_local

    interface
       subroutine eddy_diff_zisocl_layer_energy_codon(i_c, k_c, pcols_c, pver_c, tunl_mode_c, leng_mode_c, tunl_c, &
            ctunl_c, cleng_c, lbulk_c, vk_c, z_p, zi_p, n2_p, s2_p, leng_max_p, dzinc_p, dl2n2_p, dl2s2_p) &
            bind(c, name="eddy_diff_zisocl_layer_energy_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: i_c, k_c, pcols_c, pver_c, tunl_mode_c, leng_mode_c
         real(c_double), value :: tunl_c, ctunl_c, cleng_c, lbulk_c, vk_c
         type(c_ptr), value :: z_p, zi_p, n2_p, s2_p, leng_max_p, dzinc_p, dl2n2_p, dl2s2_p
       end subroutine eddy_diff_zisocl_layer_energy_codon
    end interface

    tunl_mode_local = 0
    if( choice_tunl .eq. 'rampcl' ) then
        tunl_mode_local = 1
    elseif( choice_tunl .eq. 'rampsl' ) then
        tunl_mode_local = 2
    end if

    leng_mode_local = 1
    if( choice_leng .eq. 'origin' ) then
        leng_mode_local = 0
    end if

    call eddy_diff_zisocl_layer_energy_select_impl()

    if (use_native_zisocl_layer_energy_impl) then
       call eddy_diff_zisocl_layer_energy_native(i_local, k_local, pcols_local, pver_local, lbulk_local, z_local, zi_local, &
            n2_local, s2_local, leng_max_local, dzinc_local, dl2n2_local, dl2s2_local)
       return
    end if

    call eddy_diff_zisocl_layer_energy_codon(int(i_local, c_int64_t), int(k_local, c_int64_t), int(pcols_local, c_int64_t), &
         int(pver_local, c_int64_t), int(tunl_mode_local, c_int64_t), int(leng_mode_local, c_int64_t), tunl, ctunl, cleng, &
         lbulk_local, vk, c_loc(z_local), c_loc(zi_local), c_loc(n2_local), c_loc(s2_local), c_loc(leng_max_local), &
         c_loc(dzinc_local), c_loc(dl2n2_local), c_loc(dl2s2_local))

  end subroutine eddy_diff_zisocl_layer_energy

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_layer_energy_native(i_local, k_local, pcols_local, pver_local, lbulk_local, z_local, zi_local, &
       n2_local, s2_local, leng_max_local, dzinc_local, dl2n2_local, dl2s2_local)

    implicit none

    integer, intent(in) :: i_local, k_local, pcols_local, pver_local
    real(r8), intent(in) :: lbulk_local
    real(r8), intent(in) :: z_local(pcols_local,pver_local), zi_local(pcols_local,pver_local+1)
    real(r8), intent(in) :: n2_local(pcols_local,pver_local), s2_local(pcols_local,pver_local)
    real(r8), intent(in) :: leng_max_local(pver_local)
    real(r8), intent(out) :: dzinc_local, dl2n2_local, dl2s2_local

    real(r8) :: tunlramp_local, lz_local

    if( choice_tunl .eq. 'rampcl' ) then
        tunlramp_local = 0.5_r8*(1._r8+ctunl)*tunl
    elseif( choice_tunl .eq. 'rampsl' ) then
        tunlramp_local = ctunl*tunl
    else
        tunlramp_local = tunl
    endif
    if( choice_leng .eq. 'origin' ) then
        lz_local = ( (vk*zi_local(i_local,k_local))**(-cleng) + (tunlramp_local*lbulk_local)**(-cleng) )**(-1._r8/cleng)
    else
        lz_local = min( vk*zi_local(i_local,k_local), tunlramp_local*lbulk_local )
    endif
    lz_local = min(leng_max_local(k_local), lz_local)

    dzinc_local = z_local(i_local,k_local-1) - z_local(i_local,k_local)
    dl2n2_local = lz_local*lz_local*n2_local(i_local,k_local)*dzinc_local
    dl2s2_local = lz_local*lz_local*s2_local(i_local,k_local)*dzinc_local

  end subroutine eddy_diff_zisocl_layer_energy_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_interface_energy_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (zisocl_interface_energy_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_ZISOCL_INTERFACE_ENERGY_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_zisocl_interface_energy_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_zisocl_interface_energy_impl = .false.
    end if

    zisocl_interface_energy_impl_selected = .true.

    if (masterproc) then
       if (use_native_zisocl_interface_energy_impl) then
          write(iulog,*) 'eddy_diff_zisocl_interface_energy implementation = native'
       else
          write(iulog,*) 'eddy_diff_zisocl_interface_energy implementation = codon'
       end if
    end if

  end subroutine eddy_diff_zisocl_interface_energy_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_interface_energy(i_local, k_local, pcols_local, pver_local, lbulk_local, sh_local, sm_local, &
       z_local, zi_local, n2_local, s2_local, leng_max_local, dzinc_local, dl2n2_local, dl2s2_local, dwinc_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: i_local, k_local, pcols_local, pver_local
    real(r8), intent(in) :: lbulk_local, sh_local, sm_local
    real(r8), target, intent(in) :: z_local(pcols_local,pver_local), zi_local(pcols_local,pver_local+1)
    real(r8), target, intent(in) :: n2_local(pcols_local,pver_local), s2_local(pcols_local,pver_local)
    real(r8), target, intent(in) :: leng_max_local(pver_local)
    real(r8), target, intent(out) :: dzinc_local, dl2n2_local, dl2s2_local, dwinc_local

    integer :: tunl_mode_local, leng_mode_local

    interface
       subroutine eddy_diff_zisocl_interface_energy_codon(i_c, k_c, pcols_c, pver_c, tunl_mode_c, leng_mode_c, tunl_c, &
            ctunl_c, cleng_c, lbulk_c, sh_c, sm_c, vk_c, z_p, zi_p, n2_p, s2_p, leng_max_p, dzinc_p, dl2n2_p, dl2s2_p, &
            dwinc_p) bind(c, name="eddy_diff_zisocl_interface_energy_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: i_c, k_c, pcols_c, pver_c, tunl_mode_c, leng_mode_c
         real(c_double), value :: tunl_c, ctunl_c, cleng_c, lbulk_c, sh_c, sm_c, vk_c
         type(c_ptr), value :: z_p, zi_p, n2_p, s2_p, leng_max_p, dzinc_p, dl2n2_p, dl2s2_p, dwinc_p
       end subroutine eddy_diff_zisocl_interface_energy_codon
    end interface

    tunl_mode_local = 0
    if( choice_tunl .eq. 'rampcl' ) then
        tunl_mode_local = 1
    elseif( choice_tunl .eq. 'rampsl' ) then
        tunl_mode_local = 2
    end if

    leng_mode_local = 1
    if( choice_leng .eq. 'origin' ) then
        leng_mode_local = 0
    end if

    call eddy_diff_zisocl_interface_energy_select_impl()

    if (use_native_zisocl_interface_energy_impl) then
       call eddy_diff_zisocl_interface_energy_native(i_local, k_local, pcols_local, pver_local, lbulk_local, sh_local, sm_local, &
            z_local, zi_local, n2_local, s2_local, leng_max_local, dzinc_local, dl2n2_local, dl2s2_local, dwinc_local)
       return
    end if

    call eddy_diff_zisocl_interface_energy_codon(int(i_local, c_int64_t), int(k_local, c_int64_t), int(pcols_local, c_int64_t), &
         int(pver_local, c_int64_t), int(tunl_mode_local, c_int64_t), int(leng_mode_local, c_int64_t), tunl, ctunl, cleng, &
         lbulk_local, sh_local, sm_local, vk, c_loc(z_local), c_loc(zi_local), c_loc(n2_local), c_loc(s2_local), &
         c_loc(leng_max_local), c_loc(dzinc_local), c_loc(dl2n2_local), c_loc(dl2s2_local), c_loc(dwinc_local))

  end subroutine eddy_diff_zisocl_interface_energy

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_zisocl_interface_energy_native(i_local, k_local, pcols_local, pver_local, lbulk_local, sh_local, sm_local, &
       z_local, zi_local, n2_local, s2_local, leng_max_local, dzinc_local, dl2n2_local, dl2s2_local, dwinc_local)

    implicit none

    integer, intent(in) :: i_local, k_local, pcols_local, pver_local
    real(r8), intent(in) :: lbulk_local, sh_local, sm_local
    real(r8), intent(in) :: z_local(pcols_local,pver_local), zi_local(pcols_local,pver_local+1)
    real(r8), intent(in) :: n2_local(pcols_local,pver_local), s2_local(pcols_local,pver_local)
    real(r8), intent(in) :: leng_max_local(pver_local)
    real(r8), intent(out) :: dzinc_local, dl2n2_local, dl2s2_local, dwinc_local

    real(r8) :: tunlramp_local, lz_local

    if( choice_tunl .eq. 'rampcl' ) then
        tunlramp_local = 0.5_r8*(1._r8+ctunl)*tunl
    elseif( choice_tunl .eq. 'rampsl' ) then
        tunlramp_local = ctunl*tunl
    else
        tunlramp_local = tunl
    endif
    if( choice_leng .eq. 'origin' ) then
        lz_local = ( (vk*zi_local(i_local,k_local))**(-cleng) + (tunlramp_local*lbulk_local)**(-cleng) )**(-1._r8/cleng)
    else
        lz_local = min( vk*zi_local(i_local,k_local), tunlramp_local*lbulk_local )
    endif
    lz_local = min(leng_max_local(k_local), lz_local)

    dzinc_local = z_local(i_local,k_local-1) - z_local(i_local,k_local)
    dl2n2_local = lz_local*lz_local*n2_local(i_local,k_local)*dzinc_local
    dl2s2_local = lz_local*lz_local*s2_local(i_local,k_local)*dzinc_local
    dwinc_local = -sh_local*dl2n2_local + sm_local*dl2s2_local

  end subroutine eddy_diff_zisocl_interface_energy_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_compute_radf_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (compute_radf_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_COMPUTE_RADF_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_compute_radf_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_compute_radf_impl = .false.
    end if

    compute_radf_impl_selected = .true.

    if (masterproc) then
       if (use_native_compute_radf_impl) then
          write(iulog,*) 'eddy_diff_compute_radf implementation = native'
       else
          write(iulog,*) 'eddy_diff_compute_radf implementation = codon'
       end if
    end if

  end subroutine eddy_diff_compute_radf_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_compute_radf(choice_radf_local, i_local, pcols_local, pver_local, ncvmax_local, ncvfin_local, &
       ktop_local, qmin_local, ql_local, pi_local, qrlw_local, g_local, cldeff_local, zi_local, chs_local, lwp_CL_local, &
       opt_depth_CL_local, radinvfrac_CL_local, radf_CL_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    character(len=*), intent(in) :: choice_radf_local
    integer, intent(in) :: i_local, pcols_local, pver_local, ncvmax_local
    integer(i4), target, intent(in) :: ncvfin_local(pcols_local), ktop_local(pcols_local,ncvmax_local)
    real(r8), intent(in) :: qmin_local, g_local
    real(r8), target, intent(in) :: ql_local(pcols_local,pver_local), pi_local(pcols_local,pver_local+1)
    real(r8), target, intent(in) :: qrlw_local(pcols_local,pver_local), cldeff_local(pcols_local,pver_local)
    real(r8), target, intent(in) :: zi_local(pcols_local,pver_local+1), chs_local(pcols_local,pver_local+1)
    real(r8), target, intent(inout) :: lwp_CL_local(pcols_local,ncvmax_local), opt_depth_CL_local(pcols_local,ncvmax_local)
    real(r8), target, intent(inout) :: radinvfrac_CL_local(pcols_local,ncvmax_local), radf_CL_local(pcols_local,ncvmax_local)

    integer :: radf_mode

    interface
       subroutine eddy_diff_compute_radf_codon(i_c, pcols_c, pver_c, ncvmax_c, radf_mode_c, qmin_c, g_c, ncvfin_p, &
            ktop_p, ql_p, pi_p, qrlw_p, cldeff_p, zi_p, chs_p, lwp_CL_p, opt_depth_CL_p, radinvfrac_CL_p, radf_CL_p) &
            bind(c, name="eddy_diff_compute_radf_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: i_c, pcols_c, pver_c, ncvmax_c, radf_mode_c
         real(c_double), value :: qmin_c, g_c
         type(c_ptr), value :: ncvfin_p, ktop_p, ql_p, pi_p, qrlw_p, cldeff_p, zi_p, chs_p, lwp_CL_p, opt_depth_CL_p, &
              radinvfrac_CL_p, radf_CL_p
       end subroutine eddy_diff_compute_radf_codon
    end interface

    call eddy_diff_compute_radf_select_impl()

    if (use_native_compute_radf_impl) then
       call eddy_diff_compute_radf_native(choice_radf_local, i_local, pcols_local, pver_local, ncvmax_local, ncvfin_local, &
            ktop_local, qmin_local, ql_local, pi_local, qrlw_local, g_local, cldeff_local, zi_local, chs_local, lwp_CL_local, &
            opt_depth_CL_local, radinvfrac_CL_local, radf_CL_local)
       return
    end if

    radf_mode = 2
    if (choice_radf_local .eq. 'orig') then
       radf_mode = 0
    elseif (choice_radf_local .eq. 'ramp') then
       radf_mode = 1
    endif

    call eddy_diff_compute_radf_codon(int(i_local, c_int64_t), int(pcols_local, c_int64_t), int(pver_local, c_int64_t), &
         int(ncvmax_local, c_int64_t), int(radf_mode, c_int64_t), real(qmin_local, c_double), real(g_local, c_double), &
         c_loc(ncvfin_local), c_loc(ktop_local), c_loc(ql_local), c_loc(pi_local), c_loc(qrlw_local), c_loc(cldeff_local), &
         c_loc(zi_local), c_loc(chs_local), c_loc(lwp_CL_local), c_loc(opt_depth_CL_local), c_loc(radinvfrac_CL_local), &
         c_loc(radf_CL_local))

  end subroutine eddy_diff_compute_radf

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_compute_radf_native(choice_radf_local, i_local, pcols_local, pver_local, ncvmax_local, ncvfin_local, &
       ktop_local, qmin_local, ql_local, pi_local, qrlw_local, g_local, cldeff_local, zi_local, chs_local, lwp_CL_local, &
       opt_depth_CL_local, radinvfrac_CL_local, radf_CL_local)

    use pbl_utils, only: compute_radf

    implicit none

    character(len=*), intent(in) :: choice_radf_local
    integer, intent(in) :: i_local, pcols_local, pver_local, ncvmax_local
    integer(i4), intent(in) :: ncvfin_local(pcols_local), ktop_local(pcols_local,ncvmax_local)
    real(r8), intent(in) :: qmin_local, g_local
    real(r8), intent(in) :: ql_local(pcols_local,pver_local), pi_local(pcols_local,pver_local+1)
    real(r8), intent(in) :: qrlw_local(pcols_local,pver_local), cldeff_local(pcols_local,pver_local)
    real(r8), intent(in) :: zi_local(pcols_local,pver_local+1), chs_local(pcols_local,pver_local+1)
    real(r8), intent(inout) :: lwp_CL_local(pcols_local,ncvmax_local), opt_depth_CL_local(pcols_local,ncvmax_local)
    real(r8), intent(inout) :: radinvfrac_CL_local(pcols_local,ncvmax_local), radf_CL_local(pcols_local,ncvmax_local)

    call compute_radf(choice_radf_local, i_local, pcols_local, pver_local, ncvmax_local, ncvfin_local, ktop_local, &
         qmin_local, ql_local, pi_local, qrlw_local, g_local, cldeff_local, zi_local, chs_local, lwp_CL_local(i_local,:), &
         opt_depth_CL_local(i_local,:), radinvfrac_CL_local(i_local,:), radf_CL_local(i_local,:))

  end subroutine eddy_diff_compute_radf_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_srcl_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (caleddy_srcl_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_CALEDDY_SRCL_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_caleddy_srcl_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_caleddy_srcl_impl = .false.
    end if

    caleddy_srcl_impl_selected = .true.

    if (masterproc) then
       if (use_native_caleddy_srcl_impl) then
          write(iulog,*) 'eddy_diff_caleddy_srcl implementation = native'
       else
          write(iulog,*) 'eddy_diff_caleddy_srcl implementation = codon'
       end if
    end if

  end subroutine eddy_diff_caleddy_srcl_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_srcl(choice_srcl_local, i_local, pcols_local, pver_local, ncvmax_local, ntop_turb_local, &
       nbot_turb_local, qmin_local, ricrit_local, b1_local, vk_local, alph1_local, alph2_local, alph3_local, &
       alph4exs_local, alph5_local, ghmin_local, ql_local, qrlw_local, ri_local, sfuh_local, chu_local, chs_local, &
       cmu_local, cms_local, slslope_local, qtslope_local, z_local, bflxs_local, tkes_local, bprod_local, sprod_local, &
       ncvfin_local, kbase_local, ktop_local, belongcv_local, ricl_local, ghcl_local, shcl_local, smcl_local, lbrk_local, &
       wbrk_local, ebrk_local, ncvsurf_local, srcl_status_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    character(len=*), intent(in) :: choice_srcl_local
    integer, intent(in) :: i_local, pcols_local, pver_local, ncvmax_local, ntop_turb_local, nbot_turb_local
    real(r8), intent(in) :: qmin_local, ricrit_local, b1_local, vk_local
    real(r8), intent(in) :: alph1_local, alph2_local, alph3_local, alph4exs_local, alph5_local, ghmin_local
    real(r8), target, intent(in) :: ql_local(pcols_local,pver_local), qrlw_local(pcols_local,pver_local)
    real(r8), target, intent(in) :: ri_local(pcols_local,pver_local), sfuh_local(pcols_local,pver_local)
    real(r8), target, intent(in) :: chu_local(pcols_local,pver_local+1), chs_local(pcols_local,pver_local+1)
    real(r8), target, intent(in) :: cmu_local(pcols_local,pver_local+1), cms_local(pcols_local,pver_local+1)
    real(r8), target, intent(in) :: slslope_local(pcols_local,pver_local), qtslope_local(pcols_local,pver_local)
    real(r8), target, intent(in) :: z_local(pcols_local,pver_local), bflxs_local(pcols_local), tkes_local(pcols_local)
    real(r8), target, intent(in) :: bprod_local(pcols_local,pver_local+1), sprod_local(pcols_local,pver_local+1)
    integer(i4), target, intent(inout) :: ncvfin_local(pcols_local), kbase_local(pcols_local,ncvmax_local)
    integer(i4), target, intent(inout) :: ktop_local(pcols_local,ncvmax_local)
    logical, intent(inout) :: belongcv_local(pcols_local,pver_local+1)
    real(r8), target, intent(inout) :: ricl_local(pcols_local,ncvmax_local), ghcl_local(pcols_local,ncvmax_local)
    real(r8), target, intent(inout) :: shcl_local(pcols_local,ncvmax_local), smcl_local(pcols_local,ncvmax_local)
    real(r8), target, intent(inout) :: lbrk_local(pcols_local,ncvmax_local), wbrk_local(pcols_local,ncvmax_local)
    real(r8), target, intent(inout) :: ebrk_local(pcols_local,ncvmax_local)
    integer(i4), intent(inout) :: ncvsurf_local
    integer(i4), intent(out) :: srcl_status_local

    integer :: srcl_mode, k, ncv
    integer(i4), target :: belongcv_mask_local(pver_local+1)
    integer(i4), target :: ncvsurf_codon, srcl_status_codon

    interface
       subroutine eddy_diff_caleddy_srcl_codon(i_c, pcols_c, pver_c, ncvmax_c, ntop_turb_c, nbot_turb_c, srcl_mode_c, &
            qmin_c, ricrit_c, b1_c, vk_c, alph1_c, alph2_c, alph3_c, alph4exs_c, alph5_c, ghmin_c, ql_p, qrlw_p, ri_p, &
            sfuh_p, chu_p, chs_p, cmu_p, cms_p, slslope_p, qtslope_p, z_p, bflxs_p, tkes_p, bprod_p, sprod_p, ncvfin_p, &
            kbase_p, ktop_p, ricl_p, ghcl_p, shcl_p, smcl_p, lbrk_p, wbrk_p, ebrk_p, belong_mask_p, ncvsurf_p, &
            srcl_status_p) bind(c, name="eddy_diff_caleddy_srcl_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: i_c, pcols_c, pver_c, ncvmax_c, ntop_turb_c, nbot_turb_c, srcl_mode_c
         real(c_double), value :: qmin_c, ricrit_c, b1_c, vk_c, alph1_c, alph2_c, alph3_c, alph4exs_c, alph5_c, ghmin_c
         type(c_ptr), value :: ql_p, qrlw_p, ri_p, sfuh_p, chu_p, chs_p, cmu_p, cms_p, slslope_p, qtslope_p, z_p, bflxs_p
         type(c_ptr), value :: tkes_p, bprod_p, sprod_p, ncvfin_p, kbase_p, ktop_p, ricl_p, ghcl_p, shcl_p, smcl_p
         type(c_ptr), value :: lbrk_p, wbrk_p, ebrk_p, belong_mask_p, ncvsurf_p, srcl_status_p
       end subroutine eddy_diff_caleddy_srcl_codon
    end interface

    call eddy_diff_caleddy_srcl_select_impl()

    srcl_status_local = 0_i4

    if (use_native_caleddy_srcl_impl) then
       call eddy_diff_caleddy_srcl_native(choice_srcl_local, i_local, pcols_local, pver_local, ncvmax_local, ntop_turb_local, &
            nbot_turb_local, qmin_local, ricrit_local, b1_local, vk_local, alph1_local, alph2_local, alph3_local, &
            alph4exs_local, alph5_local, ghmin_local, ql_local, qrlw_local, ri_local, sfuh_local, chu_local, chs_local, &
            cmu_local, cms_local, slslope_local, qtslope_local, z_local, bflxs_local, tkes_local, bprod_local, sprod_local, &
            ncvfin_local, kbase_local, ktop_local, belongcv_local, ricl_local, ghcl_local, shcl_local, smcl_local, lbrk_local, &
            wbrk_local, ebrk_local, ncvsurf_local, srcl_status_local)
       return
    end if

    srcl_mode = 1
    if (choice_srcl_local .eq. 'remove') then
       srcl_mode = 0
    else if (choice_srcl_local .eq. 'nonamb') then
       srcl_mode = 2
    end if

    ncvsurf_codon = ncvsurf_local
    srcl_status_codon = 0_i4

    call eddy_diff_caleddy_srcl_codon(int(i_local, c_int64_t), int(pcols_local, c_int64_t), int(pver_local, c_int64_t), &
         int(ncvmax_local, c_int64_t), int(ntop_turb_local, c_int64_t), int(nbot_turb_local, c_int64_t), &
         int(srcl_mode, c_int64_t), qmin_local, ricrit_local, b1_local, vk_local, alph1_local, alph2_local, alph3_local, &
         alph4exs_local, alph5_local, ghmin_local, c_loc(ql_local), c_loc(qrlw_local), c_loc(ri_local), c_loc(sfuh_local), &
         c_loc(chu_local), c_loc(chs_local), c_loc(cmu_local), c_loc(cms_local), c_loc(slslope_local), c_loc(qtslope_local), &
         c_loc(z_local), c_loc(bflxs_local), c_loc(tkes_local), c_loc(bprod_local), c_loc(sprod_local), c_loc(ncvfin_local), &
         c_loc(kbase_local), c_loc(ktop_local), c_loc(ricl_local), c_loc(ghcl_local), c_loc(shcl_local), c_loc(smcl_local), &
         c_loc(lbrk_local), c_loc(wbrk_local), c_loc(ebrk_local), c_loc(belongcv_mask_local), c_loc(ncvsurf_codon), &
         c_loc(srcl_status_codon))

    ncvsurf_local = ncvsurf_codon
    srcl_status_local = srcl_status_codon

    do k = 1, pver_local + 1
       belongcv_local(i_local,k) = .false.
    end do

    do ncv = 1, ncvfin_local(i_local)
       do k = ktop_local(i_local,ncv), kbase_local(i_local,ncv)
          belongcv_local(i_local,k) = .true.
       end do
    end do

  end subroutine eddy_diff_caleddy_srcl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_srcl_native(choice_srcl_local, i_local, pcols_local, pver_local, ncvmax_local, ntop_turb_local, &
       nbot_turb_local, qmin_local, ricrit_local, b1_local, vk_local, alph1_local, alph2_local, alph3_local, &
       alph4exs_local, alph5_local, ghmin_local, ql_local, qrlw_local, ri_local, sfuh_local, chu_local, chs_local, &
       cmu_local, cms_local, slslope_local, qtslope_local, z_local, bflxs_local, tkes_local, bprod_local, sprod_local, &
       ncvfin_local, kbase_local, ktop_local, belongcv_local, ricl_local, ghcl_local, shcl_local, smcl_local, lbrk_local, &
       wbrk_local, ebrk_local, ncvsurf_local, srcl_status_local)

    implicit none

    character(len=*), intent(in) :: choice_srcl_local
    integer, intent(in) :: i_local, pcols_local, pver_local, ncvmax_local, ntop_turb_local, nbot_turb_local
    real(r8), intent(in) :: qmin_local, ricrit_local, b1_local, vk_local
    real(r8), intent(in) :: alph1_local, alph2_local, alph3_local, alph4exs_local, alph5_local, ghmin_local
    real(r8), intent(in) :: ql_local(pcols_local,pver_local), qrlw_local(pcols_local,pver_local)
    real(r8), intent(in) :: ri_local(pcols_local,pver_local), sfuh_local(pcols_local,pver_local)
    real(r8), intent(in) :: chu_local(pcols_local,pver_local+1), chs_local(pcols_local,pver_local+1)
    real(r8), intent(in) :: cmu_local(pcols_local,pver_local+1), cms_local(pcols_local,pver_local+1)
    real(r8), intent(in) :: slslope_local(pcols_local,pver_local), qtslope_local(pcols_local,pver_local)
    real(r8), intent(in) :: z_local(pcols_local,pver_local), bflxs_local(pcols_local), tkes_local(pcols_local)
    real(r8), intent(in) :: bprod_local(pcols_local,pver_local+1), sprod_local(pcols_local,pver_local+1)
    integer(i4), intent(inout) :: ncvfin_local(pcols_local), kbase_local(pcols_local,ncvmax_local)
    integer(i4), intent(inout) :: ktop_local(pcols_local,ncvmax_local)
    logical, intent(inout) :: belongcv_local(pcols_local,pver_local+1)
    real(r8), intent(inout) :: ricl_local(pcols_local,ncvmax_local), ghcl_local(pcols_local,ncvmax_local)
    real(r8), intent(inout) :: shcl_local(pcols_local,ncvmax_local), smcl_local(pcols_local,ncvmax_local)
    real(r8), intent(inout) :: lbrk_local(pcols_local,ncvmax_local), wbrk_local(pcols_local,ncvmax_local)
    real(r8), intent(inout) :: ebrk_local(pcols_local,ncvmax_local)
    integer(i4), intent(inout) :: ncvsurf_local
    integer(i4), intent(out) :: srcl_status_local

    logical :: in_cl_local
    integer :: k, ncv, ncvf, ncvnew
    real(r8) :: ch, cm, n2htSRCL, gg, gh

    srcl_status_local = 0_i4

    ncv  = 1
    ncvf = ncvfin_local(i_local)

    if( choice_srcl_local .eq. 'remove' ) goto 222

    do k = nbot_turb_local, ntop_turb_local + 1, -1

       if( ql_local(i_local,k) .gt. qmin_local .and. ql_local(i_local,k-1) .lt. qmin_local .and. &
            qrlw_local(i_local,k) .lt. 0._r8 .and. ri_local(i_local,k) .ge. ricrit_local ) then

           if( choice_srcl_local .eq. 'nonamb' .and. belongcv_local(i_local,k+1) ) then
               go to 220
           endif

           ch = ( 1._r8 - sfuh_local(i_local,k) ) * chu_local(i_local,k) + sfuh_local(i_local,k) * chs_local(i_local,k)
           cm = ( 1._r8 - sfuh_local(i_local,k) ) * cmu_local(i_local,k) + sfuh_local(i_local,k) * cms_local(i_local,k)

           n2htSRCL = ch * slslope_local(i_local,k) + cm * qtslope_local(i_local,k)

           if( n2htSRCL .le. 0._r8 ) then

               in_cl_local = .false.

               do while ( ncv .le. ncvf )
                  if( ktop_local(i_local,ncv) .le. k ) then
                     if( kbase_local(i_local,ncv) .gt. k ) then
                        in_cl_local = .true.
                     endif
                     exit
                  else
                     ncv = ncv + 1
                  end if
               end do

               if( .not. in_cl_local ) then

                  ncvfin_local(i_local)       =  ncvfin_local(i_local) + 1
                  ncvnew                      =  ncvfin_local(i_local)
                  ktop_local(i_local,ncvnew)  =  k
                  kbase_local(i_local,ncvnew) =  k+1
                  belongcv_local(i_local,k)   = .true.
                  belongcv_local(i_local,k+1) = .true.

                  if( k .lt. pver_local ) then

                      wbrk_local(i_local,ncvnew) = 0._r8
                      ebrk_local(i_local,ncvnew) = 0._r8
                      lbrk_local(i_local,ncvnew) = 0._r8
                      ghcl_local(i_local,ncvnew) = 0._r8
                      shcl_local(i_local,ncvnew) = 0._r8
                      smcl_local(i_local,ncvnew) = 0._r8
                      ricl_local(i_local,ncvnew) = 0._r8

                  else

                      if( bflxs_local(i_local) .gt. 0._r8 ) then
                          ebrk_local(i_local,ncvnew) = tkes_local(i_local)
                          lbrk_local(i_local,ncvnew) = z_local(i_local,pver_local)
                          wbrk_local(i_local,ncvnew) = tkes_local(i_local) / b1_local
                          srcl_status_local = 1_i4
                          return

                      else

                          ebrk_local(i_local,ncvnew) = 0._r8
                          lbrk_local(i_local,ncvnew) = 0._r8
                          wbrk_local(i_local,ncvnew) = 0._r8

                      endif

                      gg = 0.5_r8 * vk_local * z_local(i_local,pver_local) * bprod_local(i_local,pver_local+1) / &
                           ( tkes_local(i_local)**(3._r8/2._r8) )
                      if( abs(alph5_local-gg*alph3_local) .le. 1.e-7_r8 ) then
                         gh = ghmin_local
                      else
                         gh = gg / ( alph5_local - gg * alph3_local )
                      end if
                      gh = min(max(gh,ghmin_local),0.0233_r8)
                      ghcl_local(i_local,ncvnew) =  gh
                      shcl_local(i_local,ncvnew) =  max(0._r8,alph5_local/(1._r8+alph3_local*gh))
                      smcl_local(i_local,ncvnew) =  max(0._r8,(alph1_local + alph2_local*gh)/(1._r8+alph3_local*gh)/ &
                           (1._r8+alph4exs_local*gh))
                      ricl_local(i_local,ncvnew) = -(smcl_local(i_local,ncvnew)/shcl_local(i_local,ncvnew))* &
                           (bprod_local(i_local,pver_local+1)/sprod_local(i_local,pver_local+1))

                      ncvsurf_local = ncvnew

                   end if

               end if

           end if

       end if

220 continue

    end do

222 continue

  end subroutine eddy_diff_caleddy_srcl_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_stl_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (caleddy_stl_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_CALEDDY_STL_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_caleddy_stl_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_caleddy_stl_impl = .false.
    end if

    caleddy_stl_impl_selected = .true.

    if (masterproc) then
       if (use_native_caleddy_stl_impl) then
          write(iulog,*) 'eddy_diff_caleddy_stl implementation = native'
       else
          write(iulog,*) 'eddy_diff_caleddy_stl implementation = codon'
       end if
    end if

  end subroutine eddy_diff_caleddy_stl_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_stl(i_local, pcols_local, pver_local, ncvmax_local, tunl_mode_local, leng_mode_local, &
       ricrit_local, tunl_local, ctunl_local, cleng_local, lbulk_max_local, tkemax_local, b1_local, ae_local, alph1_local, &
       alph2_local, alph3_local, alph4exs_local, alph5_local, ghmin_local, vk_local, fak_local, cpair_local, ri_local, &
       z_local, zi_local, pi_local, n2_local, s2_local, shflx_local, qflx_local, rrho_local, ustar_local, leng_max_local, &
       ncvfin_local, ktop_local, kbase_local, kvh_local, kvm_local, leng_local, tke_local, wcap_local, bprod_local, &
       sprod_local, turbtype_local, sm_aw_local, pblh_local, pblhp_local, wpert_local, tpert_local, qpert_local, &
       ipbl_local, kpblh_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: i_local, pcols_local, pver_local, ncvmax_local
    integer, intent(in) :: tunl_mode_local, leng_mode_local
    real(r8), intent(in) :: ricrit_local, tunl_local, ctunl_local, cleng_local, lbulk_max_local, tkemax_local
    real(r8), intent(in) :: b1_local, ae_local, alph1_local, alph2_local, alph3_local, alph4exs_local, alph5_local
    real(r8), intent(in) :: ghmin_local, vk_local, fak_local, cpair_local
    real(r8), target, intent(in) :: ri_local(pcols_local,pver_local), z_local(pcols_local,pver_local)
    real(r8), target, intent(in) :: zi_local(pcols_local,pver_local+1), pi_local(pcols_local,pver_local+1)
    real(r8), target, intent(in) :: n2_local(pcols_local,pver_local), s2_local(pcols_local,pver_local)
    real(r8), target, intent(in) :: shflx_local(pcols_local), qflx_local(pcols_local), rrho_local(pcols_local)
    real(r8), target, intent(in) :: ustar_local(pcols_local), leng_max_local(pver_local+1)
    integer(i4), target, intent(in) :: ncvfin_local(pcols_local), ktop_local(pcols_local,ncvmax_local)
    integer(i4), target, intent(in) :: kbase_local(pcols_local,ncvmax_local)
    real(r8), target, intent(inout) :: kvh_local(pcols_local,pver_local+1), kvm_local(pcols_local,pver_local+1)
    real(r8), target, intent(inout) :: leng_local(pcols_local,pver_local+1), tke_local(pcols_local,pver_local+1)
    real(r8), target, intent(inout) :: wcap_local(pcols_local,pver_local+1), bprod_local(pcols_local,pver_local+1)
    real(r8), target, intent(inout) :: sprod_local(pcols_local,pver_local+1), sm_aw_local(pcols_local,pver_local+1)
    integer(i4), target, intent(inout) :: turbtype_local(pcols_local,pver_local+1), ipbl_local(pcols_local), kpblh_local(pcols_local)
    real(r8), target, intent(inout) :: pblh_local(pcols_local), pblhp_local(pcols_local), wpert_local(pcols_local)
    real(r8), target, intent(inout) :: tpert_local(pcols_local), qpert_local(pcols_local)

    integer(i4), target :: clmask_local(pver_local+1), stlmask_local(pver_local+1)

    interface
       subroutine eddy_diff_caleddy_stl_codon(i_c, pcols_c, pver_c, ncvmax_c, tunl_mode_c, leng_mode_c, ricrit_c, tunl_c, &
            ctunl_c, cleng_c, lbulk_max_c, tkemax_c, b1_c, ae_c, alph1_c, alph2_c, alph3_c, alph4exs_c, alph5_c, ghmin_c, &
            vk_c, fak_c, cpair_c, ri_p, z_p, zi_p, pi_p, n2_p, s2_p, shflx_p, qflx_p, rrho_p, ustar_p, leng_max_p, &
            ncvfin_p, ktop_p, kbase_p, kvh_p, kvm_p, leng_p, tke_p, wcap_p, bprod_p, sprod_p, turbtype_p, sm_aw_p, pblh_p, &
            pblhp_p, wpert_p, tpert_p, qpert_p, ipbl_p, kpblh_p, clmask_p, stlmask_p) bind(c, name="eddy_diff_caleddy_stl_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: i_c, pcols_c, pver_c, ncvmax_c, tunl_mode_c, leng_mode_c
         real(c_double), value :: ricrit_c, tunl_c, ctunl_c, cleng_c, lbulk_max_c, tkemax_c, b1_c, ae_c, alph1_c, alph2_c
         real(c_double), value :: alph3_c, alph4exs_c, alph5_c, ghmin_c, vk_c, fak_c, cpair_c
         type(c_ptr), value :: ri_p, z_p, zi_p, pi_p, n2_p, s2_p, shflx_p, qflx_p, rrho_p, ustar_p, leng_max_p
         type(c_ptr), value :: ncvfin_p, ktop_p, kbase_p, kvh_p, kvm_p, leng_p, tke_p, wcap_p, bprod_p, sprod_p
         type(c_ptr), value :: turbtype_p, sm_aw_p, pblh_p, pblhp_p, wpert_p, tpert_p, qpert_p, ipbl_p, kpblh_p
         type(c_ptr), value :: clmask_p, stlmask_p
       end subroutine eddy_diff_caleddy_stl_codon
    end interface

    call eddy_diff_caleddy_stl_select_impl()

    if (use_native_caleddy_stl_impl) then
       call eddy_diff_caleddy_stl_native(i_local, pcols_local, pver_local, ncvmax_local, tunl_mode_local, leng_mode_local, &
            ricrit_local, tunl_local, ctunl_local, cleng_local, lbulk_max_local, tkemax_local, b1_local, ae_local, &
            alph1_local, alph2_local, alph3_local, alph4exs_local, alph5_local, ghmin_local, vk_local, fak_local, &
            cpair_local, ri_local, z_local, zi_local, pi_local, n2_local, s2_local, shflx_local, qflx_local, rrho_local, &
            ustar_local, leng_max_local, ncvfin_local, ktop_local, kbase_local, kvh_local, kvm_local, leng_local, tke_local, &
            wcap_local, bprod_local, sprod_local, turbtype_local, sm_aw_local, pblh_local, pblhp_local, wpert_local, &
            tpert_local, qpert_local, ipbl_local, kpblh_local)
       return
    end if

    call eddy_diff_caleddy_stl_codon(int(i_local, c_int64_t), int(pcols_local, c_int64_t), int(pver_local, c_int64_t), &
         int(ncvmax_local, c_int64_t), int(tunl_mode_local, c_int64_t), int(leng_mode_local, c_int64_t), ricrit_local, &
         tunl_local, ctunl_local, cleng_local, lbulk_max_local, tkemax_local, b1_local, ae_local, alph1_local, alph2_local, &
         alph3_local, alph4exs_local, alph5_local, ghmin_local, vk_local, fak_local, cpair_local, c_loc(ri_local), &
         c_loc(z_local), c_loc(zi_local), c_loc(pi_local), c_loc(n2_local), c_loc(s2_local), c_loc(shflx_local), &
         c_loc(qflx_local), c_loc(rrho_local), c_loc(ustar_local), c_loc(leng_max_local), c_loc(ncvfin_local), c_loc(ktop_local), &
         c_loc(kbase_local), c_loc(kvh_local), c_loc(kvm_local), c_loc(leng_local), c_loc(tke_local), c_loc(wcap_local), &
         c_loc(bprod_local), c_loc(sprod_local), c_loc(turbtype_local), c_loc(sm_aw_local), c_loc(pblh_local), &
         c_loc(pblhp_local), c_loc(wpert_local), c_loc(tpert_local), c_loc(qpert_local), c_loc(ipbl_local), c_loc(kpblh_local), &
         c_loc(clmask_local), c_loc(stlmask_local))

  end subroutine eddy_diff_caleddy_stl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_stl_native(i_local, pcols_local, pver_local, ncvmax_local, tunl_mode_local, leng_mode_local, &
       ricrit_local, tunl_local, ctunl_local, cleng_local, lbulk_max_local, tkemax_local, b1_local, ae_local, alph1_local, &
       alph2_local, alph3_local, alph4exs_local, alph5_local, ghmin_local, vk_local, fak_local, cpair_local, ri_local, &
       z_local, zi_local, pi_local, n2_local, s2_local, shflx_local, qflx_local, rrho_local, ustar_local, leng_max_local, &
       ncvfin_local, ktop_local, kbase_local, kvh_local, kvm_local, leng_local, tke_local, wcap_local, bprod_local, &
       sprod_local, turbtype_local, sm_aw_local, pblh_local, pblhp_local, wpert_local, tpert_local, qpert_local, &
       ipbl_local, kpblh_local)

    implicit none

    integer, intent(in) :: i_local, pcols_local, pver_local, ncvmax_local
    integer, intent(in) :: tunl_mode_local, leng_mode_local
    real(r8), intent(in) :: ricrit_local, tunl_local, ctunl_local, cleng_local, lbulk_max_local, tkemax_local
    real(r8), intent(in) :: b1_local, ae_local, alph1_local, alph2_local, alph3_local, alph4exs_local, alph5_local
    real(r8), intent(in) :: ghmin_local, vk_local, fak_local, cpair_local
    real(r8), intent(in) :: ri_local(pcols_local,pver_local), z_local(pcols_local,pver_local)
    real(r8), intent(in) :: zi_local(pcols_local,pver_local+1), pi_local(pcols_local,pver_local+1)
    real(r8), intent(in) :: n2_local(pcols_local,pver_local), s2_local(pcols_local,pver_local)
    real(r8), intent(in) :: shflx_local(pcols_local), qflx_local(pcols_local), rrho_local(pcols_local)
    real(r8), intent(in) :: ustar_local(pcols_local), leng_max_local(pver_local+1)
    integer(i4), intent(in) :: ncvfin_local(pcols_local), ktop_local(pcols_local,ncvmax_local)
    integer(i4), intent(in) :: kbase_local(pcols_local,ncvmax_local)
    real(r8), intent(inout) :: kvh_local(pcols_local,pver_local+1), kvm_local(pcols_local,pver_local+1)
    real(r8), intent(inout) :: leng_local(pcols_local,pver_local+1), tke_local(pcols_local,pver_local+1)
    real(r8), intent(inout) :: wcap_local(pcols_local,pver_local+1), bprod_local(pcols_local,pver_local+1)
    real(r8), intent(inout) :: sprod_local(pcols_local,pver_local+1), sm_aw_local(pcols_local,pver_local+1)
    integer(i4), intent(inout) :: turbtype_local(pcols_local,pver_local+1), ipbl_local(pcols_local), kpblh_local(pcols_local)
    real(r8), intent(inout) :: pblh_local(pcols_local), pblhp_local(pcols_local), wpert_local(pcols_local)
    real(r8), intent(inout) :: tpert_local(pcols_local), qpert_local(pcols_local)

    logical :: belongcv_local(pver_local+1), belongst_local(pver_local+1)
    integer :: k, ks, ncv, kt, kb, ktopbl_local
    real(r8) :: tunlramp, lbulk, gh, sh, sm
    real(r8) :: trma, trmb, trmc, det
    real(r8) :: leng_imsi, tke_imsi, kvh_imsi, kvm_imsi

    belongcv_local(:) = .false.
    belongst_local(:) = .false.

    do ncv = 1, ncvfin_local(i_local)
       do k = ktop_local(i_local,ncv), kbase_local(i_local,ncv)
          belongcv_local(k) = .true.
       end do
    end do

    belongst_local(1) = .false.
    do k = 2, pver_local
       belongst_local(k) = ( ri_local(i_local,k) .lt. ricrit_local ) .and. ( .not. belongcv_local(k) )
       if( belongst_local(k) .and. ( .not. belongst_local(k-1) ) ) then
           kt = k
       elseif( .not. belongst_local(k) .and. belongst_local(k-1) ) then
           kb = k - 1
           lbulk = z_local(i_local,kt-1) - z_local(i_local,kb)
           lbulk = min( lbulk, lbulk_max_local )
           do ks = kt, kb
              if( tunl_mode_local .eq. 2 ) then
                  tunlramp = max(1.e-3_r8,ctunl_local*tunl_local*exp(-log(ctunl_local)*ri_local(i_local,ks)/ricrit_local))
              else
                  tunlramp = tunl_local
              endif
              if( leng_mode_local .eq. 0 ) then
                  leng_local(i_local,ks) = ( (vk_local*zi_local(i_local,ks))**(-cleng_local) + &
                       (tunlramp*lbulk)**(-cleng_local) )**(-1._r8/cleng_local)
              else
                  leng_local(i_local,ks) = min( vk_local*zi_local(i_local,ks), tunlramp*lbulk )
              endif
              leng_local(i_local,ks) = min(leng_max_local(ks), leng_local(i_local,ks))
           end do
       end if
    end do

    belongst_local(pver_local+1) = .not. belongcv_local(pver_local+1)

    if( belongst_local(pver_local+1) ) then

        turbtype_local(i_local,pver_local+1) = 1

        if( belongst_local(pver_local) ) then
            lbulk = z_local(i_local,kt-1)
        else
            kt = pver_local+1
            lbulk = z_local(i_local,kt-1)
        end if
        lbulk = min( lbulk, lbulk_max_local )

        ktopbl_local = kt - 1
        pblh_local(i_local) = z_local(i_local,ktopbl_local)
        pblhp_local(i_local) = 0.5_r8 * ( pi_local(i_local,ktopbl_local) + pi_local(i_local,ktopbl_local+1) )

        do ks = kt, pver_local
           if( tunl_mode_local .eq. 2 ) then
               tunlramp = max(1.e-3_r8,ctunl_local*tunl_local*exp(-log(ctunl_local)*ri_local(i_local,ks)/ricrit_local))
           else
               tunlramp = tunl_local
           endif
           if( leng_mode_local .eq. 0 ) then
               leng_local(i_local,ks) = ( (vk_local*zi_local(i_local,ks))**(-cleng_local) + &
                    (tunlramp*lbulk)**(-cleng_local) )**(-1._r8/cleng_local)
           else
               leng_local(i_local,ks) = min( vk_local*zi_local(i_local,ks), tunlramp*lbulk )
           endif
           leng_local(i_local,ks) = min(leng_max_local(ks), leng_local(i_local,ks))
        end do

        wpert_local(i_local) = 0._r8
        tpert_local(i_local) = max(shflx_local(i_local)*rrho_local(i_local)/cpair_local*fak_local/ustar_local(i_local),0._r8)
        qpert_local(i_local) = max(qflx_local(i_local)*rrho_local(i_local)*fak_local/ustar_local(i_local),0._r8)

        ipbl_local(i_local) = 0
        kpblh_local(i_local) = ktopbl_local

    end if

    do k = 2, pver_local

       if( belongst_local(k) ) then

           turbtype_local(i_local,k) = 1
           trma = alph3_local*alph4exs_local*ri_local(i_local,k) + 2._r8*b1_local*(alph2_local-alph4exs_local*alph5_local*ri_local(i_local,k))
           trmb = (alph3_local+alph4exs_local)*ri_local(i_local,k) + 2._r8*b1_local*(-alph5_local*ri_local(i_local,k)+alph1_local)
           trmc = ri_local(i_local,k)
           det = max(trmb*trmb-4._r8*trma*trmc,0._r8)
           gh = (-trmb + sqrt(det))/(2._r8*trma)
           gh = min(max(gh,ghmin_local),0.0233_r8)
           sh = max(0._r8,alph5_local/(1._r8+alph3_local*gh))
           sm = max(0._r8,(alph1_local + alph2_local*gh)/(1._r8+alph3_local*gh)/(1._r8+alph4exs_local*gh))

           tke_local(i_local,k) = b1_local*(leng_local(i_local,k)**2)*(-sh*n2_local(i_local,k)+sm*s2_local(i_local,k))
           tke_local(i_local,k) = min(tke_local(i_local,k),tkemax_local)
           wcap_local(i_local,k) = tke_local(i_local,k)/b1_local
           kvh_local(i_local,k) = leng_local(i_local,k) * sqrt(tke_local(i_local,k)) * sh
           kvm_local(i_local,k) = leng_local(i_local,k) * sqrt(tke_local(i_local,k)) * sm
           bprod_local(i_local,k) = -kvh_local(i_local,k) * n2_local(i_local,k)
           sprod_local(i_local,k) = kvm_local(i_local,k) * s2_local(i_local,k)

           sm_aw_local(i_local,k) = sm/alph1_local

       end if

    end do

    do k = 2, pver_local

       if( ( turbtype_local(i_local,k) .eq. 3 ) .or. ( turbtype_local(i_local,k) .eq. 4 ) .or. &
           ( turbtype_local(i_local,k) .eq. 5 ) ) then

           trma = alph3_local*alph4exs_local*ri_local(i_local,k) + 2._r8*b1_local*(alph2_local-alph4exs_local*alph5_local*ri_local(i_local,k))
           trmb = (alph3_local+alph4exs_local)*ri_local(i_local,k) + 2._r8*b1_local*(-alph5_local*ri_local(i_local,k)+alph1_local)
           trmc = ri_local(i_local,k)
           det  = max(trmb*trmb-4._r8*trma*trmc,0._r8)
           gh   = (-trmb + sqrt(det))/(2._r8*trma)
           gh   = min(max(gh,ghmin_local),0.0233_r8)
           sh   = max(0._r8,alph5_local/(1._r8+alph3_local*gh))
           sm   = max(0._r8,(alph1_local + alph2_local*gh)/(1._r8+alph3_local*gh)/(1._r8+alph4exs_local*gh))

           lbulk = z_local(i_local,k-1) - z_local(i_local,k)
           lbulk = min( lbulk, lbulk_max_local )

           if( tunl_mode_local .eq. 2 ) then
               tunlramp = max(1.e-3_r8,ctunl_local*tunl_local*exp(-log(ctunl_local)*ri_local(i_local,k)/ricrit_local))
           else
               tunlramp = tunl_local
           endif
           if( leng_mode_local .eq. 0 ) then
               leng_imsi = ( (vk_local*zi_local(i_local,k))**(-cleng_local) + (tunlramp*lbulk)**(-cleng_local) )**(-1._r8/cleng_local)
           else
               leng_imsi = min( vk_local*zi_local(i_local,k), tunlramp*lbulk )
           endif
           leng_imsi = min(leng_max_local(k), leng_imsi)

           tke_imsi = b1_local*(leng_imsi**2)*(-sh*n2_local(i_local,k)+sm*s2_local(i_local,k))
           tke_imsi = min(max(tke_imsi,0._r8),tkemax_local)
           kvh_imsi = leng_imsi * sqrt(tke_imsi) * sh
           kvm_imsi = leng_imsi * sqrt(tke_imsi) * sm

           if( kvh_local(i_local,k) .lt. kvh_imsi ) then
               kvh_local(i_local,k) = kvh_imsi
               kvm_local(i_local,k) = kvm_imsi
               leng_local(i_local,k) = leng_imsi
               tke_local(i_local,k) = tke_imsi
               wcap_local(i_local,k) = tke_imsi / b1_local
               bprod_local(i_local,k) = -kvh_imsi * n2_local(i_local,k)
               sprod_local(i_local,k) = kvm_imsi * s2_local(i_local,k)
               sm_aw_local(i_local,k) = sm/alph1_local
               turbtype_local(i_local,k) = 1
           endif

       end if

    end do

  end subroutine eddy_diff_caleddy_stl_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_diag_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (caleddy_diag_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_CALEDDY_DIAG_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_caleddy_diag_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_caleddy_diag_impl = .false.
    end if

    caleddy_diag_impl_selected = .true.

    if (masterproc) then
       if (use_native_caleddy_diag_impl) then
          write(iulog,*) 'eddy_diff_caleddy_diag implementation = native'
       else
          write(iulog,*) 'eddy_diff_caleddy_diag implementation = codon'
       end if
    end if

  end subroutine eddy_diff_caleddy_diag_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_diag(i_local, pcols_local, pver_local, ricrit_local, tkes_local, b1_local, alph1_local, &
       alph2_local, alph3_local, alph4_local, alph4exs_local, alph5_local, ghmin_local, vk_local, z_local, ri_local, &
       bflxs_local, bprod_local, sprod_local, gh_a_local, sh_a_local, sm_a_local, ri_a_local, sm_aw_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: i_local, pcols_local, pver_local
    real(r8), intent(in) :: ricrit_local, b1_local, alph1_local, alph2_local, alph3_local, alph4_local
    real(r8), intent(in) :: alph4exs_local, alph5_local, ghmin_local, vk_local
    real(r8), target, intent(in) :: tkes_local(pcols_local), z_local(pcols_local,pver_local), ri_local(pcols_local,pver_local)
    real(r8), target, intent(in) :: bflxs_local(pcols_local), sprod_local(pcols_local,pver_local+1)
    real(r8), target, intent(inout) :: bprod_local(pcols_local,pver_local+1)
    real(r8), target, intent(inout) :: gh_a_local(pcols_local,pver_local+1), sh_a_local(pcols_local,pver_local+1)
    real(r8), target, intent(inout) :: sm_a_local(pcols_local,pver_local+1), ri_a_local(pcols_local,pver_local+1)
    real(r8), target, intent(inout) :: sm_aw_local(pcols_local,pver_local+1)

    interface
       subroutine eddy_diff_caleddy_diag_codon(i_c, pcols_c, pver_c, ricrit_c, b1_c, alph1_c, alph2_c, alph3_c, alph4_c, &
            alph4exs_c, alph5_c, ghmin_c, vk_c, tkes_p, z_p, ri_p, bflxs_p, bprod_p, sprod_p, gh_a_p, sh_a_p, sm_a_p, &
            ri_a_p, sm_aw_p) bind(c, name="eddy_diff_caleddy_diag_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: i_c, pcols_c, pver_c
         real(c_double), value :: ricrit_c, b1_c, alph1_c, alph2_c, alph3_c, alph4_c, alph4exs_c, alph5_c, ghmin_c, vk_c
         type(c_ptr), value :: tkes_p, z_p, ri_p, bflxs_p, bprod_p, sprod_p, gh_a_p, sh_a_p, sm_a_p, ri_a_p, sm_aw_p
       end subroutine eddy_diff_caleddy_diag_codon
    end interface

    call eddy_diff_caleddy_diag_select_impl()

    if (use_native_caleddy_diag_impl) then
       call eddy_diff_caleddy_diag_native(i_local, pcols_local, pver_local, ricrit_local, tkes_local, b1_local, alph1_local, &
            alph2_local, alph3_local, alph4_local, alph4exs_local, alph5_local, ghmin_local, vk_local, z_local, ri_local, &
            bflxs_local, bprod_local, sprod_local, gh_a_local, sh_a_local, sm_a_local, ri_a_local, sm_aw_local)
       return
    end if

    call eddy_diff_caleddy_diag_codon(int(i_local, c_int64_t), int(pcols_local, c_int64_t), int(pver_local, c_int64_t), &
         ricrit_local, b1_local, alph1_local, alph2_local, alph3_local, alph4_local, alph4exs_local, alph5_local, &
         ghmin_local, vk_local, c_loc(tkes_local), c_loc(z_local), c_loc(ri_local), c_loc(bflxs_local), c_loc(bprod_local), &
         c_loc(sprod_local), c_loc(gh_a_local), c_loc(sh_a_local), c_loc(sm_a_local), c_loc(ri_a_local), c_loc(sm_aw_local))

  end subroutine eddy_diff_caleddy_diag

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_diag_native(i_local, pcols_local, pver_local, ricrit_local, tkes_local, b1_local, alph1_local, &
       alph2_local, alph3_local, alph4_local, alph4exs_local, alph5_local, ghmin_local, vk_local, z_local, ri_local, &
       bflxs_local, bprod_local, sprod_local, gh_a_local, sh_a_local, sm_a_local, ri_a_local, sm_aw_local)

    implicit none

    integer, intent(in) :: i_local, pcols_local, pver_local
    real(r8), intent(in) :: ricrit_local, b1_local, alph1_local, alph2_local, alph3_local, alph4_local
    real(r8), intent(in) :: alph4exs_local, alph5_local, ghmin_local, vk_local
    real(r8), intent(in) :: tkes_local(pcols_local), z_local(pcols_local,pver_local), ri_local(pcols_local,pver_local)
    real(r8), intent(in) :: bflxs_local(pcols_local), sprod_local(pcols_local,pver_local+1)
    real(r8), intent(inout) :: bprod_local(pcols_local,pver_local+1)
    real(r8), intent(inout) :: gh_a_local(pcols_local,pver_local+1), sh_a_local(pcols_local,pver_local+1)
    real(r8), intent(inout) :: sm_a_local(pcols_local,pver_local+1), ri_a_local(pcols_local,pver_local+1)
    real(r8), intent(inout) :: sm_aw_local(pcols_local,pver_local+1)

    integer :: k
    real(r8) :: gg, gh, trma, trmb, trmc, det

    bprod_local(i_local,pver_local+1) = bflxs_local(i_local)

    gg = 0.5_r8*vk_local*z_local(i_local,pver_local)*bprod_local(i_local,pver_local+1)/(tkes_local(i_local)**(3._r8/2._r8))
    if( abs(alph5_local-gg*alph3_local) .le. 1.e-7_r8 ) then
        if( bprod_local(i_local,pver_local+1) .gt. 0._r8 ) then
            gh = -3.5334_r8
        else
            gh = ghmin_local
        endif
    else
        gh = gg/(alph5_local-gg*alph3_local)
    end if

    if( bprod_local(i_local,pver_local+1) .gt. 0._r8 ) then
        gh = min(max(gh,-3.5334_r8),0.0233_r8)
    else
        gh = min(max(gh,ghmin_local),0.0233_r8)
    endif

    gh_a_local(i_local,pver_local+1) = gh
    sh_a_local(i_local,pver_local+1) = max(0._r8,alph5_local/(1._r8+alph3_local*gh))
    if( bprod_local(i_local,pver_local+1) .gt. 0._r8 ) then
        sm_a_local(i_local,pver_local+1) = max(0._r8,(alph1_local+alph2_local*gh)/(1._r8+alph3_local*gh)/(1._r8+alph4_local*gh))
    else
        sm_a_local(i_local,pver_local+1) = max(0._r8,(alph1_local+alph2_local*gh)/(1._r8+alph3_local*gh)/(1._r8+alph4exs_local*gh))
    endif
    sm_aw_local(i_local,pver_local+1) = sm_a_local(i_local,pver_local+1)/alph1_local
    ri_a_local(i_local,pver_local+1)  = -(sm_a_local(i_local,pver_local+1)/sh_a_local(i_local,pver_local+1))* &
         (bprod_local(i_local,pver_local+1)/sprod_local(i_local,pver_local+1))

    do k = 1, pver_local
       if( ri_local(i_local,k) .lt. 0._r8 ) then
           trma = alph3_local*alph4_local*ri_local(i_local,k) + 2._r8*b1_local*(alph2_local-alph4_local*alph5_local*ri_local(i_local,k))
           trmb = (alph3_local+alph4_local)*ri_local(i_local,k) + 2._r8*b1_local*(-alph5_local*ri_local(i_local,k)+alph1_local)
           trmc = ri_local(i_local,k)
           det  = max(trmb*trmb-4._r8*trma*trmc,0._r8)
           gh   = (-trmb + sqrt(det))/(2._r8*trma)
           gh   = min(max(gh,-3.5334_r8),0.0233_r8)
           gh_a_local(i_local,k) = gh
           sh_a_local(i_local,k) = max(0._r8,alph5_local/(1._r8+alph3_local*gh))
           sm_a_local(i_local,k) = max(0._r8,(alph1_local+alph2_local*gh)/(1._r8+alph3_local*gh)/(1._r8+alph4_local*gh))
           ri_a_local(i_local,k) = ri_local(i_local,k)
       else
           if( ri_local(i_local,k) .gt. ricrit_local ) then
               gh_a_local(i_local,k) = ghmin_local
               sh_a_local(i_local,k) = 0._r8
               sm_a_local(i_local,k) = 0._r8
               ri_a_local(i_local,k) = ri_local(i_local,k)
           else
               trma = alph3_local*alph4exs_local*ri_local(i_local,k) + 2._r8*b1_local*(alph2_local-alph4exs_local*alph5_local*ri_local(i_local,k))
               trmb = (alph3_local+alph4exs_local)*ri_local(i_local,k) + 2._r8*b1_local*(-alph5_local*ri_local(i_local,k)+alph1_local)
               trmc = ri_local(i_local,k)
               det  = max(trmb*trmb-4._r8*trma*trmc,0._r8)
               gh   = (-trmb + sqrt(det))/(2._r8*trma)
               gh   = min(max(gh,ghmin_local),0.0233_r8)
               gh_a_local(i_local,k) = gh
               sh_a_local(i_local,k) = max(0._r8,alph5_local/(1._r8+alph3_local*gh))
               sm_a_local(i_local,k) = max(0._r8,(alph1_local+alph2_local*gh)/(1._r8+alph3_local*gh)/(1._r8+alph4exs_local*gh))
               ri_a_local(i_local,k) = ri_local(i_local,k)
           endif
       endif
    end do

  end subroutine eddy_diff_caleddy_diag_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_clprep_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (caleddy_clprep_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_CALEDDY_CLPREP_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_caleddy_clprep_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_caleddy_clprep_impl = .false.
    end if

    caleddy_clprep_impl_selected = .true.

    if (masterproc) then
       if (use_native_caleddy_clprep_impl) then
          write(iulog,*) 'eddy_diff_caleddy_clprep implementation = native'
       else
          write(iulog,*) 'eddy_diff_caleddy_clprep implementation = codon'
       end if
    end if

  end subroutine eddy_diff_caleddy_clprep_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_clprep(i_local, ncv_local, pcols_local, pver_local, ncvmax_local, tunl_mode_local, &
       leng_mode_local, evhc_mode_local, tunl_local, ctunl_local, cleng_local, lbulk_max_local, qmin_local, g_local, &
       vk_local, latvap_local, a2l_local, a3l_local, jbumin_local, evhcmax_local, ql_local, slv_local, sl_local, qt_local, &
       u_local, v_local, zi_local, z_local, n2_local, s2_local, sfuh_local, sflh_local, chu_local, chs_local, cmu_local, &
       cms_local, cldeff_local, bflxs_local, bprod_local, kbase_local, ktop_local, ricl_local, shcl_local, smcl_local, &
       radf_local, leng_max_local, leng_local, wcap_local, lbulk_local, jbzm_local, jbbu_local, n2hb_local, vyb_local, &
       vub_local, jtzm_local, jtbu_local, jt2slv_local, n2ht_local, vyt_local, vut_local, evhc_local, dzht_local, &
       dzhb_local, wstar3_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: i_local, ncv_local, pcols_local, pver_local, ncvmax_local
    integer, intent(in) :: tunl_mode_local, leng_mode_local, evhc_mode_local
    real(r8), intent(in) :: tunl_local, ctunl_local, cleng_local, lbulk_max_local
    real(r8), intent(in) :: qmin_local, g_local, vk_local, latvap_local
    real(r8), intent(in) :: a2l_local, a3l_local, jbumin_local, evhcmax_local, radf_local
    real(r8), target, intent(in) :: ql_local(pcols_local,pver_local), slv_local(pcols_local,pver_local)
    real(r8), target, intent(in) :: sl_local(pcols_local,pver_local), qt_local(pcols_local,pver_local)
    real(r8), target, intent(in) :: u_local(pcols_local,pver_local), v_local(pcols_local,pver_local)
    real(r8), target, intent(in) :: zi_local(pcols_local,pver_local+1), z_local(pcols_local,pver_local)
    real(r8), target, intent(in) :: n2_local(pcols_local,pver_local), s2_local(pcols_local,pver_local)
    real(r8), target, intent(in) :: sfuh_local(pcols_local,pver_local), sflh_local(pcols_local,pver_local)
    real(r8), target, intent(in) :: chu_local(pcols_local,pver_local+1), chs_local(pcols_local,pver_local+1)
    real(r8), target, intent(in) :: cmu_local(pcols_local,pver_local+1), cms_local(pcols_local,pver_local+1)
    real(r8), target, intent(in) :: cldeff_local(pcols_local,pver_local), bflxs_local(pcols_local)
    real(r8), target, intent(in) :: bprod_local(pcols_local,pver_local+1), leng_max_local(pver_local)
    integer(i4), target, intent(in) :: kbase_local(pcols_local,ncvmax_local), ktop_local(pcols_local,ncvmax_local)
    real(r8), target, intent(in) :: ricl_local(pcols_local,ncvmax_local), shcl_local(pcols_local,ncvmax_local)
    real(r8), target, intent(in) :: smcl_local(pcols_local,ncvmax_local)
    real(r8), target, intent(inout) :: leng_local(pcols_local,pver_local+1), wcap_local(pcols_local,pver_local+1)
    real(r8), intent(out) :: lbulk_local, jbzm_local, jbbu_local, n2hb_local, vyb_local, vub_local
    real(r8), intent(out) :: jtzm_local, jtbu_local, jt2slv_local, n2ht_local, vyt_local, vut_local
    real(r8), intent(out) :: evhc_local, dzht_local, dzhb_local, wstar3_local

    real(r8), target :: clprep_state_local(16)

    interface
       subroutine eddy_diff_caleddy_clprep_codon(i_c, ncv_c, pcols_c, pver_c, ncvmax_c, tunl_mode_c, leng_mode_c, &
            evhc_mode_c, tunl_c, ctunl_c, cleng_c, lbulk_max_c, qmin_c, g_c, vk_c, latvap_c, a2l_c, a3l_c, jbumin_c, &
            evhcmax_c, ql_p, slv_p, sl_p, qt_p, u_p, v_p, zi_p, z_p, n2_p, s2_p, sfuh_p, sflh_p, chu_p, chs_p, cmu_p, &
            cms_p, cldeff_p, bflxs_p, bprod_p, kbase_p, ktop_p, ricl_p, shcl_p, smcl_p, radf_c, leng_max_p, leng_p, &
            wcap_p, clprep_state_p) bind(c, name="eddy_diff_caleddy_clprep_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: i_c, ncv_c, pcols_c, pver_c, ncvmax_c, tunl_mode_c, leng_mode_c, evhc_mode_c
         real(c_double), value :: tunl_c, ctunl_c, cleng_c, lbulk_max_c, qmin_c, g_c, vk_c, latvap_c, a2l_c, a3l_c
         real(c_double), value :: jbumin_c, evhcmax_c, radf_c
         type(c_ptr), value :: ql_p, slv_p, sl_p, qt_p, u_p, v_p, zi_p, z_p, n2_p, s2_p, sfuh_p, sflh_p
         type(c_ptr), value :: chu_p, chs_p, cmu_p, cms_p, cldeff_p, bflxs_p, bprod_p, kbase_p, ktop_p
         type(c_ptr), value :: ricl_p, shcl_p, smcl_p, leng_max_p, leng_p, wcap_p, clprep_state_p
       end subroutine eddy_diff_caleddy_clprep_codon
    end interface

    call eddy_diff_caleddy_clprep_select_impl()

    if (use_native_caleddy_clprep_impl) then
       call eddy_diff_caleddy_clprep_native(i_local, ncv_local, pcols_local, pver_local, ncvmax_local, tunl_mode_local, &
            leng_mode_local, evhc_mode_local, tunl_local, ctunl_local, cleng_local, lbulk_max_local, qmin_local, g_local, &
            vk_local, latvap_local, a2l_local, a3l_local, jbumin_local, evhcmax_local, ql_local, slv_local, sl_local, &
            qt_local, u_local, v_local, zi_local, z_local, n2_local, s2_local, sfuh_local, sflh_local, chu_local, &
            chs_local, cmu_local, cms_local, cldeff_local, bflxs_local, bprod_local, kbase_local, ktop_local, ricl_local, &
            shcl_local, smcl_local, radf_local, leng_max_local, leng_local, wcap_local, lbulk_local, jbzm_local, jbbu_local, &
            n2hb_local, vyb_local, vub_local, jtzm_local, jtbu_local, jt2slv_local, n2ht_local, vyt_local, vut_local, &
            evhc_local, dzht_local, dzhb_local, wstar3_local)
       return
    end if

    call eddy_diff_caleddy_clprep_codon(int(i_local, c_int64_t), int(ncv_local, c_int64_t), int(pcols_local, c_int64_t), &
         int(pver_local, c_int64_t), int(ncvmax_local, c_int64_t), int(tunl_mode_local, c_int64_t), &
         int(leng_mode_local, c_int64_t), int(evhc_mode_local, c_int64_t), tunl_local, ctunl_local, cleng_local, &
         lbulk_max_local, qmin_local, g_local, vk_local, latvap_local, a2l_local, a3l_local, jbumin_local, evhcmax_local, &
         c_loc(ql_local), c_loc(slv_local), c_loc(sl_local), c_loc(qt_local), c_loc(u_local), c_loc(v_local), c_loc(zi_local), &
         c_loc(z_local), c_loc(n2_local), c_loc(s2_local), c_loc(sfuh_local), c_loc(sflh_local), c_loc(chu_local), &
         c_loc(chs_local), c_loc(cmu_local), c_loc(cms_local), c_loc(cldeff_local), c_loc(bflxs_local), c_loc(bprod_local), &
         c_loc(kbase_local), c_loc(ktop_local), c_loc(ricl_local), c_loc(shcl_local), c_loc(smcl_local), radf_local, &
         c_loc(leng_max_local), c_loc(leng_local), c_loc(wcap_local), c_loc(clprep_state_local))

    lbulk_local  = clprep_state_local(1)
    jbzm_local   = clprep_state_local(2)
    jbbu_local   = clprep_state_local(3)
    n2hb_local   = clprep_state_local(4)
    vyb_local    = clprep_state_local(5)
    vub_local    = clprep_state_local(6)
    jtzm_local   = clprep_state_local(7)
    jtbu_local   = clprep_state_local(8)
    jt2slv_local = clprep_state_local(9)
    n2ht_local   = clprep_state_local(10)
    vyt_local    = clprep_state_local(11)
    vut_local    = clprep_state_local(12)
    evhc_local   = clprep_state_local(13)
    dzht_local   = clprep_state_local(14)
    dzhb_local   = clprep_state_local(15)
    wstar3_local = clprep_state_local(16)

  end subroutine eddy_diff_caleddy_clprep

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_clprep_native(i_local, ncv_local, pcols_local, pver_local, ncvmax_local, tunl_mode_local, &
       leng_mode_local, evhc_mode_local, tunl_local, ctunl_local, cleng_local, lbulk_max_local, qmin_local, g_local, &
       vk_local, latvap_local, a2l_local, a3l_local, jbumin_local, evhcmax_local, ql_local, slv_local, sl_local, qt_local, &
       u_local, v_local, zi_local, z_local, n2_local, s2_local, sfuh_local, sflh_local, chu_local, chs_local, cmu_local, &
       cms_local, cldeff_local, bflxs_local, bprod_local, kbase_local, ktop_local, ricl_local, shcl_local, smcl_local, &
       radf_local, leng_max_local, leng_local, wcap_local, lbulk_local, jbzm_local, jbbu_local, n2hb_local, vyb_local, &
       vub_local, jtzm_local, jtbu_local, jt2slv_local, n2ht_local, vyt_local, vut_local, evhc_local, dzht_local, &
       dzhb_local, wstar3_local)

    implicit none

    integer, intent(in) :: i_local, ncv_local, pcols_local, pver_local, ncvmax_local
    integer, intent(in) :: tunl_mode_local, leng_mode_local, evhc_mode_local
    real(r8), intent(in) :: tunl_local, ctunl_local, cleng_local, lbulk_max_local
    real(r8), intent(in) :: qmin_local, g_local, vk_local, latvap_local
    real(r8), intent(in) :: a2l_local, a3l_local, jbumin_local, evhcmax_local, radf_local
    real(r8), intent(in) :: ql_local(pcols_local,pver_local), slv_local(pcols_local,pver_local)
    real(r8), intent(in) :: sl_local(pcols_local,pver_local), qt_local(pcols_local,pver_local)
    real(r8), intent(in) :: u_local(pcols_local,pver_local), v_local(pcols_local,pver_local)
    real(r8), intent(in) :: zi_local(pcols_local,pver_local+1), z_local(pcols_local,pver_local)
    real(r8), intent(in) :: n2_local(pcols_local,pver_local), s2_local(pcols_local,pver_local)
    real(r8), intent(in) :: sfuh_local(pcols_local,pver_local), sflh_local(pcols_local,pver_local)
    real(r8), intent(in) :: chu_local(pcols_local,pver_local+1), chs_local(pcols_local,pver_local+1)
    real(r8), intent(in) :: cmu_local(pcols_local,pver_local+1), cms_local(pcols_local,pver_local+1)
    real(r8), intent(in) :: cldeff_local(pcols_local,pver_local), bflxs_local(pcols_local)
    real(r8), intent(in) :: bprod_local(pcols_local,pver_local+1), leng_max_local(pver_local)
    integer(i4), intent(in) :: kbase_local(pcols_local,ncvmax_local), ktop_local(pcols_local,ncvmax_local)
    real(r8), intent(in) :: ricl_local(pcols_local,ncvmax_local), shcl_local(pcols_local,ncvmax_local)
    real(r8), intent(in) :: smcl_local(pcols_local,ncvmax_local)
    real(r8), intent(inout) :: leng_local(pcols_local,pver_local+1), wcap_local(pcols_local,pver_local+1)
    real(r8), intent(out) :: lbulk_local, jbzm_local, jbbu_local, n2hb_local, vyb_local, vub_local
    real(r8), intent(out) :: jtzm_local, jtbu_local, jt2slv_local, n2ht_local, vyt_local, vut_local
    real(r8), intent(out) :: evhc_local, dzht_local, dzhb_local, wstar3_local

    integer :: k, kb, kt
    real(r8) :: tunlramp, qleff
    real(r8) :: jbsl, jbqt, jbu, jbv, jtsl, jtqt, jtu, jtv
    real(r8) :: ch, cm

    kt = ktop_local(i_local,ncv_local)
    kb = kbase_local(i_local,ncv_local)

    if( kb .eq. (pver_local+1) .and. bflxs_local(i_local) .le. 0._r8 ) then
        lbulk_local = zi_local(i_local,kt) - z_local(i_local,pver_local)
    else
        lbulk_local = zi_local(i_local,kt) - zi_local(i_local,kb)
    end if
    lbulk_local = min( lbulk_local, lbulk_max_local )

    do k = min(kb,pver_local), kt, -1
       if( tunl_mode_local .eq. 1 ) then
           tunlramp = ctunl_local*tunl_local*(1._r8-(1._r8-1._r8/ctunl_local)*exp(min(0._r8,ricl_local(i_local,ncv_local))))
           tunlramp = min(max(tunlramp,tunl_local),ctunl_local*tunl_local)
       elseif( tunl_mode_local .eq. 2 ) then
           tunlramp = ctunl_local*tunl_local
       else
           tunlramp = tunl_local
       endif
       if( leng_mode_local .eq. 0 ) then
           leng_local(i_local,k) = ( (vk_local*zi_local(i_local,k))**(-cleng_local) + (tunlramp*lbulk_local)**(-cleng_local) )**(-1._r8/cleng_local)
       else
           leng_local(i_local,k) = min( vk_local*zi_local(i_local,k), tunlramp*lbulk_local )
       endif
       leng_local(i_local,k) = min(leng_max_local(k), leng_local(i_local,k))
       wcap_local(i_local,k) = (leng_local(i_local,k)**2) * (-shcl_local(i_local,ncv_local)*n2_local(i_local,k)+smcl_local(i_local,ncv_local)*s2_local(i_local,k))
    end do

    if( kb .lt. pver_local+1 ) then

        jbzm_local = z_local(i_local,kb-1) - z_local(i_local,kb)
        jbsl = sl_local(i_local,kb-1) - sl_local(i_local,kb)
        jbqt = qt_local(i_local,kb-1) - qt_local(i_local,kb)
        jbbu_local = n2_local(i_local,kb) * jbzm_local
        jbbu_local = max(jbbu_local,jbumin_local)
        jbu  = u_local(i_local,kb-1) - u_local(i_local,kb)
        jbv  = v_local(i_local,kb-1) - v_local(i_local,kb)
        ch   = (1._r8 -sflh_local(i_local,kb-1))*chu_local(i_local,kb) + sflh_local(i_local,kb-1)*chs_local(i_local,kb)
        cm   = (1._r8 -sflh_local(i_local,kb-1))*cmu_local(i_local,kb) + sflh_local(i_local,kb-1)*cms_local(i_local,kb)
        n2hb_local = (ch*jbsl + cm*jbqt)/jbzm_local
        vyb_local  = n2hb_local*jbzm_local/jbbu_local
        vub_local  = min(1._r8,(jbu**2+jbv**2)/(jbbu_local*jbzm_local))

    else

        jbzm_local = 0._r8
        jbbu_local = 0._r8
        n2hb_local = 0._r8
        vyb_local  = 0._r8
        vub_local  = 0._r8

    end if

    jtzm_local = z_local(i_local,kt-1) - z_local(i_local,kt)
    jtsl = sl_local(i_local,kt-1) - sl_local(i_local,kt)
    jtqt = qt_local(i_local,kt-1) - qt_local(i_local,kt)
    jtbu_local = n2_local(i_local,kt)*jtzm_local
    jtbu_local = max(jtbu_local,jbumin_local)
    jtu  = u_local(i_local,kt-1) - u_local(i_local,kt)
    jtv  = v_local(i_local,kt-1) - v_local(i_local,kt)
    ch   = (1._r8 -sfuh_local(i_local,kt))*chu_local(i_local,kt) + sfuh_local(i_local,kt)*chs_local(i_local,kt)
    cm   = (1._r8 -sfuh_local(i_local,kt))*cmu_local(i_local,kt) + sfuh_local(i_local,kt)*cms_local(i_local,kt)
    n2ht_local = (ch*jtsl + cm*jtqt)/jtzm_local
    vyt_local  = n2ht_local*jtzm_local/jtbu_local
    vut_local  = min(1._r8,(jtu**2+jtv**2)/(jtbu_local*jtzm_local))

    evhc_local   = 1._r8
    jt2slv_local = 0._r8

    if( evhc_mode_local .eq. 0 ) then

        if( ql_local(i_local,kt) .gt. qmin_local .and. ql_local(i_local,kt-1) .lt. qmin_local ) then
            jt2slv_local = slv_local(i_local,max(kt-2,1)) - slv_local(i_local,kt)
            jt2slv_local = max( jt2slv_local, jbumin_local*slv_local(i_local,kt-1)/g_local )
            evhc_local   = 1._r8 + a2l_local * a3l_local * latvap_local * ql_local(i_local,kt) / jt2slv_local
            evhc_local   = min( evhc_local, evhcmax_local )
        end if

    elseif( evhc_mode_local .eq. 1 ) then

        jt2slv_local = slv_local(i_local,max(kt-2,1)) - slv_local(i_local,kt)
        jt2slv_local = max( jt2slv_local, jbumin_local*slv_local(i_local,kt-1)/g_local )
        evhc_local   = 1._r8 + max(cldeff_local(i_local,kt)-cldeff_local(i_local,kt-1),0._r8) * a2l_local * a3l_local * &
             latvap_local * ql_local(i_local,kt) / jt2slv_local
        evhc_local   = min( evhc_local, evhcmax_local )

    else

        qleff        = max( ql_local(i_local,kt-1), ql_local(i_local,kt) )
        jt2slv_local = slv_local(i_local,max(kt-2,1)) - slv_local(i_local,kt)
        jt2slv_local = max( jt2slv_local, jbumin_local*slv_local(i_local,kt-1)/g_local )
        evhc_local   = 1._r8 + a2l_local * a3l_local * latvap_local * qleff / jt2slv_local
        evhc_local   = min( evhc_local, evhcmax_local )

    endif

    dzht_local   = zi_local(i_local,kt)  - z_local(i_local,kt)
    dzhb_local   = z_local(i_local,kb-1) - zi_local(i_local,kb)
    wstar3_local = radf_local * dzht_local
    do k = kt + 1, kb - 1
         wstar3_local =  wstar3_local + bprod_local(i_local,k) * ( z_local(i_local,k-1) - z_local(i_local,k) )
    end do
    if( kb .eq. (pver_local+1) .and. bflxs_local(i_local) .gt. 0._r8 ) then
       wstar3_local = wstar3_local + bflxs_local(i_local) * dzhb_local
    end if
    wstar3_local = max( 2.5_r8 * wstar3_local, 0._r8 )

  end subroutine eddy_diff_caleddy_clprep_native

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_closure_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (caleddy_closure_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_CALEDDY_CLOSURE_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_caleddy_closure_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_caleddy_closure_impl = .false.
    end if

    caleddy_closure_impl_selected = .true.

    if (masterproc) then
       if (use_native_caleddy_closure_impl) then
          write(iulog,*) 'eddy_diff_caleddy_closure implementation = native'
       else
          write(iulog,*) 'eddy_diff_caleddy_closure implementation = codon'
       end if
    end if

  end subroutine eddy_diff_caleddy_closure_select_impl

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_closure(i_local, pcols_local, pver_local, ncvmax_local, tunl_mode_local, leng_mode_local, &
       evhc_mode_local, wstarent_mode_local, sedfact_mode_local, tunl_local, ctunl_local, cleng_local, lbulk_max_local, &
       tkemax_local, b1_local, ae_local, alph1_local, a1l_local, a1i_local, ccrit_local, wstar3factcrit_local, &
       ntzero_local, onet_local, rcapmin_local, rcapmax_local, wfac_local, wpertmin_local, tfac_local, qmin_local, &
       g_local, vk_local, cpair_local, latvap_local, a2l_local, a3l_local, jbumin_local, evhcmax_local, ased_local, &
       ql_local, slv_local, sl_local, qt_local, u_local, v_local, pi_local, zi_local, z_local, n2_local, s2_local, &
       shflx_local, qflx_local, rrho_local, sfuh_local, sflh_local, chu_local, chs_local, cmu_local, cms_local, &
       cldeff_local, bflxs_local, bprod_local, sprod_local, wsedl_local, ncvfin_local, kbase_local, ktop_local, &
       belongcv_local, lbrk_local, ebrk_local, wbrk_local, ricl_local, shcl_local, smcl_local, radf_CL_local, &
       wsed_CL_local, leng_max_local, wet_CL_local, web_CL_local, jtbu_CL_local, jbbu_CL_local, evhc_CL_local, &
       jt2slv_CL_local, n2ht_CL_local, n2hb_CL_local, wstar_CL_local, wstar3fact_CL_local, leng_local, wcap_local, &
       tke_local, kvh_local, kvm_local, turbtype_local, sm_aw_local, pblh_local, pblhp_local, wpert_local, tpert_local, &
       qpert_local, ipbl_local, kpblh_local, went_local, ncvsurf_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: i_local, pcols_local, pver_local, ncvmax_local
    integer, intent(in) :: tunl_mode_local, leng_mode_local, evhc_mode_local
    integer, intent(in) :: wstarent_mode_local, sedfact_mode_local
    real(r8), intent(in) :: tunl_local, ctunl_local, cleng_local, lbulk_max_local, tkemax_local
    real(r8), intent(in) :: b1_local, ae_local, alph1_local, a1l_local, a1i_local, ccrit_local
    real(r8), intent(in) :: wstar3factcrit_local, ntzero_local, onet_local, rcapmin_local, rcapmax_local
    real(r8), intent(in) :: wfac_local, wpertmin_local, tfac_local, qmin_local, g_local, vk_local, cpair_local
    real(r8), intent(in) :: latvap_local, a2l_local, a3l_local, jbumin_local, evhcmax_local, ased_local
    real(r8), target, intent(in) :: ql_local(pcols_local,pver_local), slv_local(pcols_local,pver_local)
    real(r8), target, intent(in) :: sl_local(pcols_local,pver_local), qt_local(pcols_local,pver_local)
    real(r8), target, intent(in) :: u_local(pcols_local,pver_local), v_local(pcols_local,pver_local)
    real(r8), target, intent(in) :: pi_local(pcols_local,pver_local+1), zi_local(pcols_local,pver_local+1)
    real(r8), target, intent(in) :: z_local(pcols_local,pver_local), n2_local(pcols_local,pver_local)
    real(r8), target, intent(in) :: s2_local(pcols_local,pver_local), shflx_local(pcols_local), qflx_local(pcols_local)
    real(r8), target, intent(in) :: rrho_local(pcols_local), sfuh_local(pcols_local,pver_local), sflh_local(pcols_local,pver_local)
    real(r8), target, intent(in) :: chu_local(pcols_local,pver_local+1), chs_local(pcols_local,pver_local+1)
    real(r8), target, intent(in) :: cmu_local(pcols_local,pver_local+1), cms_local(pcols_local,pver_local+1)
    real(r8), target, intent(in) :: cldeff_local(pcols_local,pver_local), bflxs_local(pcols_local)
    real(r8), target, intent(inout) :: bprod_local(pcols_local,pver_local+1), sprod_local(pcols_local,pver_local+1)
    real(r8), target, intent(in) :: wsedl_local(pcols_local,pver_local)
    integer(i4), target, intent(in) :: ncvfin_local(pcols_local), kbase_local(pcols_local,ncvmax_local)
    integer(i4), target, intent(in) :: ktop_local(pcols_local,ncvmax_local)
    logical, intent(inout) :: belongcv_local(pcols_local,pver_local+1)
    real(r8), target, intent(in) :: lbrk_local(pcols_local,ncvmax_local)
    real(r8), target, intent(inout) :: ebrk_local(pcols_local,ncvmax_local), wbrk_local(pcols_local,ncvmax_local)
    real(r8), target, intent(in) :: ricl_local(pcols_local,ncvmax_local), shcl_local(pcols_local,ncvmax_local)
    real(r8), target, intent(in) :: smcl_local(pcols_local,ncvmax_local), radf_CL_local(pcols_local,ncvmax_local)
    real(r8), target, intent(inout) :: wsed_CL_local(pcols_local,ncvmax_local), wet_CL_local(pcols_local,ncvmax_local)
    real(r8), target, intent(inout) :: web_CL_local(pcols_local,ncvmax_local), jtbu_CL_local(pcols_local,ncvmax_local)
    real(r8), target, intent(inout) :: jbbu_CL_local(pcols_local,ncvmax_local), evhc_CL_local(pcols_local,ncvmax_local)
    real(r8), target, intent(inout) :: jt2slv_CL_local(pcols_local,ncvmax_local), n2ht_CL_local(pcols_local,ncvmax_local)
    real(r8), target, intent(inout) :: n2hb_CL_local(pcols_local,ncvmax_local), wstar_CL_local(pcols_local,ncvmax_local)
    real(r8), target, intent(inout) :: wstar3fact_CL_local(pcols_local,ncvmax_local)
    real(r8), target, intent(in) :: leng_max_local(pver_local)
    real(r8), target, intent(inout) :: leng_local(pcols_local,pver_local+1), wcap_local(pcols_local,pver_local+1)
    real(r8), target, intent(inout) :: tke_local(pcols_local,pver_local+1), kvh_local(pcols_local,pver_local+1)
    real(r8), target, intent(inout) :: kvm_local(pcols_local,pver_local+1), sm_aw_local(pcols_local,pver_local+1)
    integer(i4), target, intent(inout) :: turbtype_local(pcols_local,pver_local+1)
    real(r8), target, intent(inout) :: pblh_local(pcols_local), pblhp_local(pcols_local), wpert_local(pcols_local)
    real(r8), target, intent(inout) :: tpert_local(pcols_local), qpert_local(pcols_local), went_local(pcols_local)
    integer(i4), target, intent(inout) :: ipbl_local(pcols_local), kpblh_local(pcols_local)
    integer(i4), intent(in) :: ncvsurf_local

    integer(i4), target :: zero_tke_mask_local(pver_local+1), closure_status_local(3)
    integer :: k

    interface
       subroutine eddy_diff_caleddy_closure_codon(i_c, pcols_c, pver_c, ncvmax_c, tunl_mode_c, leng_mode_c, evhc_mode_c, &
            wstarent_mode_c, sedfact_mode_c, ncvsurf_c, tunl_c, ctunl_c, cleng_c, lbulk_max_c, tkemax_c, b1_c, ae_c, &
            alph1_c, a1l_c, a1i_c, ccrit_c, wstar3factcrit_c, ntzero_c, onet_c, rcapmin_c, rcapmax_c, wfac_c, &
            wpertmin_c, tfac_c, qmin_c, g_c, vk_c, cpair_c, latvap_c, a2l_c, a3l_c, jbumin_c, evhcmax_c, ased_c, &
            ql_p, slv_p, sl_p, qt_p, u_p, v_p, pi_p, zi_p, z_p, n2_p, s2_p, shflx_p, qflx_p, rrho_p, sfuh_p, sflh_p, &
            chu_p, chs_p, cmu_p, cms_p, cldeff_p, bflxs_p, bprod_p, sprod_p, wsedl_p, ncvfin_p, kbase_p, ktop_p, lbrk_p, &
            ebrk_p, wbrk_p, ricl_p, shcl_p, smcl_p, radf_CL_p, wsed_CL_p, leng_max_p, wet_CL_p, web_CL_p, jtbu_CL_p, &
            jbbu_CL_p, evhc_CL_p, jt2slv_CL_p, n2ht_CL_p, n2hb_CL_p, wstar_CL_p, wstar3fact_CL_p, leng_p, wcap_p, tke_p, &
            kvh_p, kvm_p, turbtype_p, sm_aw_p, pblh_p, pblhp_p, wpert_p, tpert_p, qpert_p, ipbl_p, kpblh_p, went_p, &
            zero_tke_mask_p, closure_status_p) bind(c, name="eddy_diff_caleddy_closure_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: i_c, pcols_c, pver_c, ncvmax_c, tunl_mode_c, leng_mode_c, evhc_mode_c
         integer(c_int64_t), value :: wstarent_mode_c, sedfact_mode_c, ncvsurf_c
         real(c_double), value :: tunl_c, ctunl_c, cleng_c, lbulk_max_c, tkemax_c, b1_c, ae_c, alph1_c, a1l_c, a1i_c
         real(c_double), value :: ccrit_c, wstar3factcrit_c, ntzero_c, onet_c, rcapmin_c, rcapmax_c, wfac_c, wpertmin_c
         real(c_double), value :: tfac_c, qmin_c, g_c, vk_c, cpair_c, latvap_c, a2l_c, a3l_c, jbumin_c, evhcmax_c, ased_c
         type(c_ptr), value :: ql_p, slv_p, sl_p, qt_p, u_p, v_p, pi_p, zi_p, z_p, n2_p, s2_p, shflx_p, qflx_p, rrho_p
         type(c_ptr), value :: sfuh_p, sflh_p, chu_p, chs_p, cmu_p, cms_p, cldeff_p, bflxs_p, bprod_p, sprod_p, wsedl_p
         type(c_ptr), value :: ncvfin_p, kbase_p, ktop_p, lbrk_p, ebrk_p, wbrk_p, ricl_p, shcl_p, smcl_p, radf_CL_p
         type(c_ptr), value :: wsed_CL_p, leng_max_p, wet_CL_p, web_CL_p, jtbu_CL_p, jbbu_CL_p, evhc_CL_p, jt2slv_CL_p
         type(c_ptr), value :: n2ht_CL_p, n2hb_CL_p, wstar_CL_p, wstar3fact_CL_p, leng_p, wcap_p, tke_p, kvh_p, kvm_p
         type(c_ptr), value :: turbtype_p, sm_aw_p, pblh_p, pblhp_p, wpert_p, tpert_p, qpert_p, ipbl_p, kpblh_p, went_p
         type(c_ptr), value :: zero_tke_mask_p, closure_status_p
       end subroutine eddy_diff_caleddy_closure_codon
    end interface

    call eddy_diff_caleddy_closure_select_impl()

    if (use_native_caleddy_closure_impl) then
       call eddy_diff_caleddy_closure_native(i_local, pcols_local, pver_local, ncvmax_local, tunl_mode_local, leng_mode_local, &
            evhc_mode_local, wstarent_mode_local, sedfact_mode_local, tunl_local, ctunl_local, cleng_local, lbulk_max_local, &
            tkemax_local, b1_local, ae_local, alph1_local, a1l_local, a1i_local, ccrit_local, wstar3factcrit_local, &
            ntzero_local, onet_local, rcapmin_local, rcapmax_local, wfac_local, wpertmin_local, tfac_local, qmin_local, &
            g_local, vk_local, cpair_local, latvap_local, a2l_local, a3l_local, jbumin_local, evhcmax_local, ased_local, &
            ql_local, slv_local, sl_local, qt_local, u_local, v_local, pi_local, zi_local, z_local, n2_local, s2_local, &
            shflx_local, qflx_local, rrho_local, sfuh_local, sflh_local, chu_local, chs_local, cmu_local, cms_local, &
            cldeff_local, bflxs_local, bprod_local, sprod_local, wsedl_local, ncvfin_local, kbase_local, ktop_local, &
            belongcv_local, lbrk_local, ebrk_local, wbrk_local, ricl_local, shcl_local, smcl_local, radf_CL_local, &
            wsed_CL_local, leng_max_local, wet_CL_local, web_CL_local, jtbu_CL_local, jbbu_CL_local, evhc_CL_local, &
            jt2slv_CL_local, n2ht_CL_local, n2hb_CL_local, wstar_CL_local, wstar3fact_CL_local, leng_local, wcap_local, &
            tke_local, kvh_local, kvm_local, turbtype_local, sm_aw_local, pblh_local, pblhp_local, wpert_local, tpert_local, &
            qpert_local, ipbl_local, kpblh_local, went_local, ncvsurf_local)
       return
    end if

    call eddy_diff_caleddy_closure_codon(int(i_local, c_int64_t), int(pcols_local, c_int64_t), int(pver_local, c_int64_t), &
         int(ncvmax_local, c_int64_t), int(tunl_mode_local, c_int64_t), int(leng_mode_local, c_int64_t), &
         int(evhc_mode_local, c_int64_t), int(wstarent_mode_local, c_int64_t), int(sedfact_mode_local, c_int64_t), &
         int(ncvsurf_local, c_int64_t), tunl_local, ctunl_local, cleng_local, lbulk_max_local, tkemax_local, b1_local, &
         ae_local, alph1_local, a1l_local, a1i_local, ccrit_local, wstar3factcrit_local, ntzero_local, onet_local, &
         rcapmin_local, rcapmax_local, wfac_local, wpertmin_local, tfac_local, qmin_local, g_local, vk_local, cpair_local, &
         latvap_local, a2l_local, a3l_local, jbumin_local, evhcmax_local, ased_local, c_loc(ql_local), c_loc(slv_local), &
         c_loc(sl_local), c_loc(qt_local), c_loc(u_local), c_loc(v_local), c_loc(pi_local), c_loc(zi_local), c_loc(z_local), &
         c_loc(n2_local), c_loc(s2_local), c_loc(shflx_local), c_loc(qflx_local), c_loc(rrho_local), c_loc(sfuh_local), &
         c_loc(sflh_local), c_loc(chu_local), c_loc(chs_local), c_loc(cmu_local), c_loc(cms_local), c_loc(cldeff_local), &
         c_loc(bflxs_local), c_loc(bprod_local), c_loc(sprod_local), c_loc(wsedl_local), c_loc(ncvfin_local), c_loc(kbase_local), &
         c_loc(ktop_local), c_loc(lbrk_local), c_loc(ebrk_local), c_loc(wbrk_local), c_loc(ricl_local), c_loc(shcl_local), &
         c_loc(smcl_local), c_loc(radf_CL_local), c_loc(wsed_CL_local), c_loc(leng_max_local), c_loc(wet_CL_local), &
         c_loc(web_CL_local), c_loc(jtbu_CL_local), c_loc(jbbu_CL_local), c_loc(evhc_CL_local), c_loc(jt2slv_CL_local), &
         c_loc(n2ht_CL_local), c_loc(n2hb_CL_local), c_loc(wstar_CL_local), c_loc(wstar3fact_CL_local), c_loc(leng_local), &
         c_loc(wcap_local), c_loc(tke_local), c_loc(kvh_local), c_loc(kvm_local), c_loc(turbtype_local), c_loc(sm_aw_local), &
         c_loc(pblh_local), c_loc(pblhp_local), c_loc(wpert_local), c_loc(tpert_local), c_loc(qpert_local), c_loc(ipbl_local), &
         c_loc(kpblh_local), c_loc(went_local), c_loc(zero_tke_mask_local), c_loc(closure_status_local))

    if (closure_status_local(1) .ne. 0_i4) then
       write(iulog,*) 'CALEDDY: Warning, CL with zero TKE, i, kt, kb ', i_local, closure_status_local(2), closure_status_local(3)
       do k = 1, pver_local + 1
          if (zero_tke_mask_local(k) .ne. 0_i4) belongcv_local(i_local,k) = .false.
       end do
    end if

  end subroutine eddy_diff_caleddy_closure

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine eddy_diff_caleddy_closure_native(i_local, pcols_local, pver_local, ncvmax_local, tunl_mode_local, leng_mode_local, &
       evhc_mode_local, wstarent_mode_local, sedfact_mode_local, tunl_local, ctunl_local, cleng_local, lbulk_max_local, &
       tkemax_local, b1_local, ae_local, alph1_local, a1l_local, a1i_local, ccrit_local, wstar3factcrit_local, &
       ntzero_local, onet_local, rcapmin_local, rcapmax_local, wfac_local, wpertmin_local, tfac_local, qmin_local, &
       g_local, vk_local, cpair_local, latvap_local, a2l_local, a3l_local, jbumin_local, evhcmax_local, ased_local, &
       ql_local, slv_local, sl_local, qt_local, u_local, v_local, pi_local, zi_local, z_local, n2_local, s2_local, &
       shflx_local, qflx_local, rrho_local, sfuh_local, sflh_local, chu_local, chs_local, cmu_local, cms_local, &
       cldeff_local, bflxs_local, bprod_local, sprod_local, wsedl_local, ncvfin_local, kbase_local, ktop_local, &
       belongcv_local, lbrk_local, ebrk_local, wbrk_local, ricl_local, shcl_local, smcl_local, radf_CL_local, &
       wsed_CL_local, leng_max_local, wet_CL_local, web_CL_local, jtbu_CL_local, jbbu_CL_local, evhc_CL_local, &
       jt2slv_CL_local, n2ht_CL_local, n2hb_CL_local, wstar_CL_local, wstar3fact_CL_local, leng_local, wcap_local, &
       tke_local, kvh_local, kvm_local, turbtype_local, sm_aw_local, pblh_local, pblhp_local, wpert_local, tpert_local, &
       qpert_local, ipbl_local, kpblh_local, went_local, ncvsurf_local)

    implicit none

    integer, intent(in) :: i_local, pcols_local, pver_local, ncvmax_local
    integer, intent(in) :: tunl_mode_local, leng_mode_local, evhc_mode_local
    integer, intent(in) :: wstarent_mode_local, sedfact_mode_local
    real(r8), intent(in) :: tunl_local, ctunl_local, cleng_local, lbulk_max_local, tkemax_local
    real(r8), intent(in) :: b1_local, ae_local, alph1_local, a1l_local, a1i_local, ccrit_local
    real(r8), intent(in) :: wstar3factcrit_local, ntzero_local, onet_local, rcapmin_local, rcapmax_local
    real(r8), intent(in) :: wfac_local, wpertmin_local, tfac_local, qmin_local, g_local, vk_local, cpair_local
    real(r8), intent(in) :: latvap_local, a2l_local, a3l_local, jbumin_local, evhcmax_local, ased_local
    real(r8), intent(in) :: ql_local(pcols_local,pver_local), slv_local(pcols_local,pver_local)
    real(r8), intent(in) :: sl_local(pcols_local,pver_local), qt_local(pcols_local,pver_local)
    real(r8), intent(in) :: u_local(pcols_local,pver_local), v_local(pcols_local,pver_local)
    real(r8), intent(in) :: pi_local(pcols_local,pver_local+1), zi_local(pcols_local,pver_local+1)
    real(r8), intent(in) :: z_local(pcols_local,pver_local), n2_local(pcols_local,pver_local), s2_local(pcols_local,pver_local)
    real(r8), intent(in) :: shflx_local(pcols_local), qflx_local(pcols_local), rrho_local(pcols_local)
    real(r8), intent(in) :: sfuh_local(pcols_local,pver_local), sflh_local(pcols_local,pver_local)
    real(r8), intent(in) :: chu_local(pcols_local,pver_local+1), chs_local(pcols_local,pver_local+1)
    real(r8), intent(in) :: cmu_local(pcols_local,pver_local+1), cms_local(pcols_local,pver_local+1)
    real(r8), intent(in) :: cldeff_local(pcols_local,pver_local), bflxs_local(pcols_local)
    real(r8), intent(inout) :: bprod_local(pcols_local,pver_local+1), sprod_local(pcols_local,pver_local+1)
    real(r8), intent(in) :: wsedl_local(pcols_local,pver_local)
    integer(i4), intent(in) :: ncvfin_local(pcols_local), kbase_local(pcols_local,ncvmax_local), ktop_local(pcols_local,ncvmax_local)
    logical, intent(inout) :: belongcv_local(pcols_local,pver_local+1)
    real(r8), intent(in) :: lbrk_local(pcols_local,ncvmax_local)
    real(r8), intent(inout) :: ebrk_local(pcols_local,ncvmax_local), wbrk_local(pcols_local,ncvmax_local)
    real(r8), intent(in) :: ricl_local(pcols_local,ncvmax_local), shcl_local(pcols_local,ncvmax_local)
    real(r8), intent(in) :: smcl_local(pcols_local,ncvmax_local), radf_CL_local(pcols_local,ncvmax_local)
    real(r8), intent(inout) :: wsed_CL_local(pcols_local,ncvmax_local), wet_CL_local(pcols_local,ncvmax_local)
    real(r8), intent(inout) :: web_CL_local(pcols_local,ncvmax_local), jtbu_CL_local(pcols_local,ncvmax_local)
    real(r8), intent(inout) :: jbbu_CL_local(pcols_local,ncvmax_local), evhc_CL_local(pcols_local,ncvmax_local)
    real(r8), intent(inout) :: jt2slv_CL_local(pcols_local,ncvmax_local), n2ht_CL_local(pcols_local,ncvmax_local)
    real(r8), intent(inout) :: n2hb_CL_local(pcols_local,ncvmax_local), wstar_CL_local(pcols_local,ncvmax_local)
    real(r8), intent(inout) :: wstar3fact_CL_local(pcols_local,ncvmax_local)
    real(r8), intent(in) :: leng_max_local(pver_local)
    real(r8), intent(inout) :: leng_local(pcols_local,pver_local+1), wcap_local(pcols_local,pver_local+1)
    real(r8), intent(inout) :: tke_local(pcols_local,pver_local+1), kvh_local(pcols_local,pver_local+1)
    real(r8), intent(inout) :: kvm_local(pcols_local,pver_local+1), sm_aw_local(pcols_local,pver_local+1)
    integer(i4), intent(inout) :: turbtype_local(pcols_local,pver_local+1)
    real(r8), intent(inout) :: pblh_local(pcols_local), pblhp_local(pcols_local), wpert_local(pcols_local)
    real(r8), intent(inout) :: tpert_local(pcols_local), qpert_local(pcols_local), went_local(pcols_local)
    integer(i4), intent(inout) :: ipbl_local(pcols_local), kpblh_local(pcols_local)
    integer(i4), intent(in) :: ncvsurf_local

    integer :: ncv, k, ktblw, kb, kt, ktopbl_local
    real(r8) :: lbulk, jbzm, jbbu, n2hb, vyb, vub, jtzm, jtbu, jt2slv, n2ht, vyt, vut, evhc
    real(r8) :: dzht, dzhb, wstar3, radf, web, wet, sedfact, qleff, cet, ceb, wstar, wstar3fact
    real(r8) :: fact, trma, trmp, trmq, qq, rmin, fmin, rcrit, fcrit, rootp, rcap, kentr, dzhb5, dzht5, tke_imsi
    logical  :: noroot

    ktblw = 0

    do ncv = 1, ncvfin_local(i_local)
       kt = ktop_local(i_local,ncv)
       kb = kbase_local(i_local,ncv)
       radf = radf_CL_local(i_local,ncv)

       call eddy_diff_caleddy_clprep_native(i_local, ncv, pcols_local, pver_local, ncvmax_local, tunl_mode_local, &
            leng_mode_local, evhc_mode_local, tunl_local, ctunl_local, cleng_local, lbulk_max_local, qmin_local, g_local, &
            vk_local, latvap_local, a2l_local, a3l_local, jbumin_local, evhcmax_local, ql_local, slv_local, sl_local, &
            qt_local, u_local, v_local, zi_local, z_local, n2_local, s2_local, sfuh_local, sflh_local, chu_local, &
            chs_local, cmu_local, cms_local, cldeff_local, bflxs_local, bprod_local, kbase_local, ktop_local, ricl_local, &
            shcl_local, smcl_local, radf, leng_max_local, leng_local, wcap_local, lbulk, jbzm, jbbu, n2hb, vyb, vub, &
            jtzm, jtbu, jt2slv, n2ht, vyt, vut, evhc, dzht, dzhb, wstar3)

       web = 0._r8
       wstar = 0._r8

       if( sedfact_mode_local .ne. 0 ) then
           sedfact = exp(-ased_local*wsedl_local(i_local,kt)/(wstar3**(1._r8/3._r8)+1.e-6_r8))
           wsed_CL_local(i_local,ncv) = wsedl_local(i_local,kt)
           if( evhc_mode_local .eq. 0 ) then
               if (ql_local(i_local,kt).gt.qmin_local .and. ql_local(i_local,kt-1).lt.qmin_local) then
                   jt2slv = slv_local(i_local,max(kt-2,1)) - slv_local(i_local,kt)
                   jt2slv = max(jt2slv, jbumin_local*slv_local(i_local,kt-1)/g_local)
                   evhc = 1._r8+sedfact*a2l_local*a3l_local*latvap_local*ql_local(i_local,kt) / jt2slv
                   evhc = min(evhc,evhcmax_local)
               end if
           elseif( evhc_mode_local .eq. 1 ) then
               jt2slv = slv_local(i_local,max(kt-2,1)) - slv_local(i_local,kt)
               jt2slv = max(jt2slv, jbumin_local*slv_local(i_local,kt-1)/g_local)
               evhc = 1._r8+max(cldeff_local(i_local,kt)-cldeff_local(i_local,kt-1),0._r8)*sedfact*a2l_local*a3l_local* &
                    latvap_local*ql_local(i_local,kt) / jt2slv
               evhc = min(evhc,evhcmax_local)
           else
               qleff  = max(ql_local(i_local,kt-1),ql_local(i_local,kt))
               jt2slv = slv_local(i_local,max(kt-2,1)) - slv_local(i_local,kt)
               jt2slv = max(jt2slv, jbumin_local*slv_local(i_local,kt-1)/g_local)
               evhc = 1._r8+sedfact*a2l_local*a3l_local*latvap_local*qleff / jt2slv
               evhc = min(evhc,evhcmax_local)
           endif
       end if

       if( wstar3 .gt. 0._r8 ) then
           cet = a1i_local * evhc / ( jtbu * lbulk )
           if( kb .eq. pver_local + 1 ) then
               wstar3fact = max( 1._r8 + 2.5_r8 * cet * n2ht * jtzm * dzht, wstar3factcrit_local )
           else
               ceb = a1i_local / ( jbbu * lbulk )
               wstar3fact = max( 1._r8 + 2.5_r8 * cet * n2ht * jtzm * dzht + 2.5_r8 * ceb * n2hb * jbzm * dzhb, &
                    wstar3factcrit_local )
           end if
           wstar3 = wstar3 / wstar3fact
       else
           wstar3fact = 0._r8
           cet        = 0._r8
           ceb        = 0._r8
       end if

       fact = ( evhc * ( -vyt + vut ) * dzht + ( -vyb + vub ) * dzhb * leng_local(i_local,kb) / leng_local(i_local,kt) ) / lbulk

       if( wstarent_mode_local .ne. 0 ) then
           trma = 1._r8
           trmp = ebrk_local(i_local,ncv) * ( lbrk_local(i_local,ncv) / lbulk ) / 3._r8 + ntzero_local
           trmq = 0.5_r8 * b1_local * ( leng_local(i_local,kt) / lbulk ) * ( radf * dzht + a1i_local * fact * wstar3 )

           rmin  = sqrt(trmp)
           fmin  = rmin * ( rmin * rmin - 3._r8 * trmp ) - 2._r8 * trmq
           wstar = wstar3**onet_local
           rcrit = ccrit_local * wstar
           fcrit = rcrit * ( rcrit * rcrit - 3._r8 * trmp ) - 2._r8 * trmq
           noroot = ( ( rmin .lt. rcrit ) .and. ( fcrit .gt. 0._r8 ) ) .or. ( ( rmin .ge. rcrit ) .and. ( fmin  .gt. 0._r8 ) )
           if( noroot ) then
               trma = 1._r8 - b1_local * ( leng_local(i_local,kt) / lbulk ) * a1i_local * fact / ccrit_local**3
               trma = max( trma, 0.5_r8 )
               trmp = trmp / trma
               trmq = 0.5_r8 * b1_local * ( leng_local(i_local,kt) / lbulk ) * radf * dzht / trma
           end if

           qq = trmq**2 - trmp**3
           if( qq .ge. 0._r8 ) then
               rootp = ( trmq + sqrt(qq) )**(1._r8/3._r8) + ( max( trmq - sqrt(qq), 0._r8 ) )**(1._r8/3._r8)
           else
               rootp = 2._r8 * sqrt(trmp) * cos( acos( trmq / sqrt(trmp**3) ) / 3._r8 )
           end if

           if( noroot )  wstar3 = ( rootp / ccrit_local )**3
           wet = cet * wstar3
           if( kb .lt. pver_local + 1 ) web = ceb * wstar3
       else
           trma = 1._r8 - b1_local * a1l_local * fact
           trma = max( trma, 0.5_r8 )
           trmp = ebrk_local(i_local,ncv) * ( lbrk_local(i_local,ncv) / lbulk ) / ( 3._r8 * trma )
           trmq = 0.5_r8 * b1_local * ( leng_local(i_local,kt) / lbulk ) * radf * dzht / trma

           qq = trmq**2 - trmp**3
           if( qq .ge. 0._r8 ) then
               rootp = ( trmq + sqrt(qq) )**(1._r8/3._r8) + ( max( trmq - sqrt(qq), 0._r8 ) )**(1._r8/3._r8)
           else
               rootp = 2._r8 * sqrt(trmp) * cos( acos( trmq / sqrt(trmp**3) ) / 3._r8 )
           end if

           wet = a1l_local * rootp * min( evhc * rootp**2 / ( leng_local(i_local,kt) * jtbu ), 1._r8 )
           if( kb .lt. pver_local + 1 ) web = a1l_local * rootp * min( evhc * rootp**2 / ( leng_local(i_local,kb) * jbbu ), 1._r8 )
       end if

       ebrk_local(i_local,ncv) = rootp**2
       ebrk_local(i_local,ncv) = min(ebrk_local(i_local,ncv),tkemax_local)
       wbrk_local(i_local,ncv) = ebrk_local(i_local,ncv)/b1_local

       if( ebrk_local(i_local,ncv) .le. 0._r8 ) then
           write(iulog,*) 'CALEDDY: Warning, CL with zero TKE, i, kt, kb ', i_local, kt, kb
           belongcv_local(i_local,kt) = .false.
           belongcv_local(i_local,kb) = .false.
       end if

       do k = kb - 1, kt + 1, -1
          rcap = ( b1_local * ae_local + wcap_local(i_local,k) / wbrk_local(i_local,ncv) ) / ( b1_local * ae_local + 1._r8 )
          rcap = min( max(rcap,rcapmin_local), rcapmax_local )
          tke_local(i_local,k) = ebrk_local(i_local,ncv) * rcap
          tke_local(i_local,k) = min( tke_local(i_local,k), tkemax_local )
          kvh_local(i_local,k) = leng_local(i_local,k) * sqrt(tke_local(i_local,k)) * shcl_local(i_local,ncv)
          kvm_local(i_local,k) = leng_local(i_local,k) * sqrt(tke_local(i_local,k)) * smcl_local(i_local,ncv)
          bprod_local(i_local,k) = -kvh_local(i_local,k) * n2_local(i_local,k)
          sprod_local(i_local,k) =  kvm_local(i_local,k) * s2_local(i_local,k)
          turbtype_local(i_local,k) = 2
          sm_aw_local(i_local,k) = smcl_local(i_local,ncv)/alph1_local
       end do

       kentr = wet * jtzm
       kvh_local(i_local,kt) = kentr
       kvm_local(i_local,kt) = kentr
       bprod_local(i_local,kt) = -kentr * n2ht + radf
       sprod_local(i_local,kt) =  kentr * s2_local(i_local,kt)
       turbtype_local(i_local,kt) = 4
       trmp = -b1_local * ae_local / ( 1._r8 + b1_local * ae_local )
       trmq = -(bprod_local(i_local,kt)+sprod_local(i_local,kt))*b1_local*leng_local(i_local,kt) / &
            (1._r8+b1_local*ae_local)/(ebrk_local(i_local,ncv)**(3._r8/2._r8))
       rcap = compute_cubic(0._r8,trmp,trmq)**2._r8
       rcap = min( max(rcap,rcapmin_local), rcapmax_local )
       tke_local(i_local,kt)  = ebrk_local(i_local,ncv) * rcap
       tke_local(i_local,kt)  = min( tke_local(i_local,kt), tkemax_local )
       sm_aw_local(i_local,kt) = smcl_local(i_local,ncv) / alph1_local

       if( kb .lt. pver_local + 1 ) then
           kentr = web * jbzm
           if( kb .ne. ktblw ) then
               kvh_local(i_local,kb) = kentr
               kvm_local(i_local,kb) = kentr
               bprod_local(i_local,kb) = -kvh_local(i_local,kb)*n2hb
               sprod_local(i_local,kb) =  kvm_local(i_local,kb)*s2_local(i_local,kb)
               turbtype_local(i_local,kb) = 3
               trmp = -b1_local*ae_local/(1._r8+b1_local*ae_local)
               trmq = -(bprod_local(i_local,kb)+sprod_local(i_local,kb))*b1_local*leng_local(i_local,kb) / &
                    (1._r8+b1_local*ae_local)/(ebrk_local(i_local,ncv)**(3._r8/2._r8))
               rcap = compute_cubic(0._r8,trmp,trmq)**2._r8
               rcap = min( max(rcap,rcapmin_local), rcapmax_local )
               tke_local(i_local,kb)  = ebrk_local(i_local,ncv) * rcap
               tke_local(i_local,kb)  = min( tke_local(i_local,kb),tkemax_local )
           else
               kvh_local(i_local,kb) = kvh_local(i_local,kb) + kentr
               kvm_local(i_local,kb) = kvm_local(i_local,kb) + kentr
               dzhb5 = z_local(i_local,kb-1) - zi_local(i_local,kb)
               dzht5 = zi_local(i_local,kb) - z_local(i_local,kb)
               bprod_local(i_local,kb) = ( dzht5*bprod_local(i_local,kb) - dzhb5*kentr*n2hb ) / ( dzhb5 + dzht5 )
               sprod_local(i_local,kb) = ( dzht5*sprod_local(i_local,kb) + dzhb5*kentr*s2_local(i_local,kb) ) / ( dzhb5 + dzht5 )
               trmp = -b1_local*ae_local/(1._r8+b1_local*ae_local)
               trmq = -kentr*(s2_local(i_local,kb)-n2hb)*b1_local*leng_local(i_local,kb) / &
                    (1._r8+b1_local*ae_local)/(ebrk_local(i_local,ncv)**(3._r8/2._r8))
               rcap = compute_cubic(0._r8,trmp,trmq)**2._r8
               rcap = min( max(rcap,rcapmin_local), rcapmax_local )
               tke_imsi = ebrk_local(i_local,ncv) * rcap
               tke_imsi = min( tke_imsi, tkemax_local )
               tke_local(i_local,kb)  = ( dzht5*tke_local(i_local,kb) + dzhb5*tke_imsi ) / ( dzhb5 + dzht5 )
               tke_local(i_local,kb)  = min(tke_local(i_local,kb),tkemax_local)
               turbtype_local(i_local,kb) = 5
           end if
       else
           rcap = (b1_local*ae_local + wcap_local(i_local,kb)/wbrk_local(i_local,ncv))/(b1_local*ae_local + 1._r8)
           rcap = min( max(rcap,rcapmin_local), rcapmax_local )
           tke_local(i_local,kb) = ebrk_local(i_local,ncv) * rcap
           tke_local(i_local,kb) = min( tke_local(i_local,kb),tkemax_local )
       end if

       sm_aw_local(i_local,kb) = smcl_local(i_local,ncv)/alph1_local

       wcap_local(i_local,kt) = (bprod_local(i_local,kt)+sprod_local(i_local,kt))*leng_local(i_local,kt)/sqrt(max(tke_local(i_local,kt),1.e-6_r8))
       if( kb .lt. pver_local + 1 ) then
           wcap_local(i_local,kb) = (bprod_local(i_local,kb)+sprod_local(i_local,kb))*leng_local(i_local,kb)/sqrt(max(tke_local(i_local,kb),1.e-6_r8))
       end if

       ktblw = kt

       wet_CL_local(i_local,ncv)        = wet
       web_CL_local(i_local,ncv)        = web
       jtbu_CL_local(i_local,ncv)       = jtbu
       jbbu_CL_local(i_local,ncv)       = jbbu
       evhc_CL_local(i_local,ncv)       = evhc
       jt2slv_CL_local(i_local,ncv)     = jt2slv
       n2ht_CL_local(i_local,ncv)       = n2ht
       n2hb_CL_local(i_local,ncv)       = n2hb
       wstar_CL_local(i_local,ncv)      = wstar
       wstar3fact_CL_local(i_local,ncv) = wstar3fact
    end do

    if( ncvsurf_local .gt. 0 ) then
        ktopbl_local = ktop_local(i_local,ncvsurf_local)
        pblh_local(i_local)   = zi_local(i_local, ktopbl_local)
        pblhp_local(i_local)  = pi_local(i_local, ktopbl_local)
        wpert_local(i_local)  = max(wfac_local*sqrt(ebrk_local(i_local,ncvsurf_local)),wpertmin_local)
        tpert_local(i_local)  = max(abs(shflx_local(i_local)*rrho_local(i_local)/cpair_local)*tfac_local/wpert_local(i_local),0._r8)
        qpert_local(i_local)  = max(abs(qflx_local(i_local)*rrho_local(i_local))*tfac_local/wpert_local(i_local),0._r8)
        if( bflxs_local(i_local) .gt. 0._r8 ) then
            turbtype_local(i_local,pver_local+1) = 2
        else
            turbtype_local(i_local,pver_local+1) = 3
        endif
        ipbl_local(i_local)  = 1
        kpblh_local(i_local) = max(ktopbl_local-1, 1)
        went_local(i_local)  = wet_CL_local(i_local,ncvsurf_local)
    end if

  end subroutine eddy_diff_caleddy_closure_native

    ! ---------------------------------------------------------------------------- !
    !                                                                              !
    ! The University of Washington Moist Turbulence Scheme                         !
    !                                                                              !
    ! Authors : Chris Bretherton at the University of Washington, Seattle, WA      ! 
    !           Sungsu Park at the CGD/NCAR, Boulder, CO                           !
    !                                                                              !
    ! ---------------------------------------------------------------------------- !

    subroutine caleddy( pcols        , pver         , ncol        ,                             &
                        sl           , qt           , ql          , slv        , u            , &
                        v            , pi           , z           , zi         ,                &
                        qflx         , shflx        , slslope     , qtslope    ,                &
                        chu          , chs          , cmu         , cms        , sfuh         , &
                        sflh         , n2           , s2          , ri         , rrho         , &
                        pblh         , ustar        ,                                           &
                        kvh_in       , kvm_in       , kvh         , kvm        ,                &
                        tpert        , qpert        , qrlin       , kvf        , tke          , & 
                        wstarent     , bprod        , sprod       , minpblh    , wpert        , &
                        tkes         , went         , turbtype    , sm_aw      ,                &
                        kbase_o      , ktop_o       , ncvfin_o    ,                             & 
                        kbase_mg     , ktop_mg      , ncvfin_mg   ,                             & 
                        kbase_f      , ktop_f       , ncvfin_f    ,                             & 
                        wet_CL       , web_CL       , jtbu_CL     , jbbu_CL    ,                &
                        evhc_CL      , jt2slv_CL    , n2ht_CL     , n2hb_CL    , lwp_CL       , &
                        opt_depth_CL , radinvfrac_CL, radf_CL     , wstar_CL   , wstar3fact_CL, &
                        ebrk         , wbrk         , lbrk        , ricl       , ghcl         , & 
                        shcl         , smcl         ,                                           &
                        gh_a         , sh_a         , sm_a        , ri_a       , leng         , & 
                        wcap         , pblhp        , cld         , ipbl       , kpblh        , &
                        wsedl        , wsed_CL )

    !--------------------------------------------------------------------------------- !
    !                                                                                  !
    ! Purpose : This is a driver routine to compute eddy diffusion coefficients        !
    !           for heat (sl), momentum (u, v), moisture (qt), and other  trace        !
    !           constituents.   This scheme uses first order closure for stable        !
    !           turbulent layers (STL). For convective layers (CL), entrainment        !
    !           closure is used at the CL external interfaces, which is coupled        !
    !           to the diagnosis of a CL regime mean TKE from the instantaneous        !
    !           thermodynamic and velocity profiles.   The CLs are diagnosed by        !
    !           extending original CL layers of moist static instability   into        !
    !           adjacent weakly stably stratified interfaces,   stopping if the        !
    !           stability is too strong.   This allows a realistic depiction of        !
    !           dry convective boundary layers with a downgradient approach.           !
    !                                                                                  !   
    ! NOTE:     This routine currently assumes ntop_turb = 1, nbot_turb = pver         !
    !           ( turbulent diffusivities computed at all interior interfaces )        !
    !           and will require modification to handle a different ntop_turb.         ! 
    !                                                                                  !
    ! Authors:  Sungsu Park and Chris Bretherton. 08/2006, 05/2008.                    !
    !                                                                                  ! 
    ! For details, see                                                                 !
    !                                                                                  !
    ! 1. 'A new moist turbulence parametrization in the Community Atmosphere Model'    !
    !     by Christopher S. Bretherton & Sungsu Park. J. Climate. 22. 3422-3448. 2009. !
    !                                                                                  !
    ! 2. 'The University of Washington shallow convection and moist turbulence schemes !
    !     and their impact on climate simulations with the Community Atmosphere Model' !
    !     by Sungsu Park & Christopher S. Bretherton. J. Climate. 22. 3449-3469. 2009. !
    !                                                                                  !
    ! For questions on the scheme and code, send an email to                           !
    !     sungsup@ucar.edu or breth@washington.edu                                     !
    !                                                                                  !
    !--------------------------------------------------------------------------------- !

    ! ---------------- !
    ! Inputs variables !
    ! ---------------- !

    implicit none
    integer,  intent(in) :: pcols                     ! Number of atmospheric columns   
    integer,  intent(in) :: pver                      ! Number of atmospheric layers   
    integer,  intent(in) :: ncol                      ! Number of atmospheric columns   
    real(r8), intent(in) :: u(pcols,pver)             ! U wind [ m/s ]
    real(r8), intent(in) :: v(pcols,pver)             ! V wind [ m/s ]
    real(r8), intent(in) :: sl(pcols,pver)            ! Liquid water static energy, cp * T + g * z - Lv * ql - Ls * qi [ J/kg ]
    real(r8), intent(in) :: slv(pcols,pver)           ! Liquid water virtual static energy, sl * ( 1 + 0.608 * qt ) [ J/kg ]
    real(r8), intent(in) :: qt(pcols,pver)            ! Total speccific humidity  qv + ql + qi [ kg/kg ] 
    real(r8), intent(in) :: ql(pcols,pver)            ! Liquid water specific humidity [ kg/kg ]
    real(r8), intent(in) :: pi(pcols,pver+1)          ! Interface pressures [ Pa ]
    real(r8), intent(in) :: z(pcols,pver)             ! Layer midpoint height above surface [ m ]
    real(r8), intent(in) :: zi(pcols,pver+1)          ! Interface height above surface, i.e., zi(pver+1) = 0 all over the globe
                                                      ! [ m ]
    real(r8), intent(in) :: chu(pcols,pver+1)         ! Buoyancy coeffi. unsaturated sl (heat) coef. at all interfaces.
                                                      ! [ unit ? ]
    real(r8), intent(in) :: chs(pcols,pver+1)         ! Buoyancy coeffi. saturated sl (heat) coef. at all interfaces.
                                                      ! [ unit ? ]
    real(r8), intent(in) :: cmu(pcols,pver+1)         ! Buoyancy coeffi. unsaturated qt (moisture) coef. at all interfaces
                                                      ! [ unit ? ]
    real(r8), intent(in) :: cms(pcols,pver+1)         ! Buoyancy coeffi. saturated qt (moisture) coef. at all interfaces
                                                      ! [ unit ? ]
    real(r8), intent(in) :: sfuh(pcols,pver)          ! Saturation fraction in upper half-layer [ fraction ]
    real(r8), intent(in) :: sflh(pcols,pver)          ! Saturation fraction in lower half-layer [ fraction ]
    real(r8), intent(in) :: n2(pcols,pver)            ! Interfacial (except surface) moist buoyancy frequency [ s-2 ]
    real(r8), intent(in) :: s2(pcols,pver)            ! Interfacial (except surface) shear frequency [ s-2 ]
    real(r8), intent(in) :: ri(pcols,pver)            ! Interfacial (except surface) Richardson number
    real(r8), intent(in) :: qflx(pcols)               ! Kinematic surface constituent ( water vapor ) flux [ kg/m2/s ]
    real(r8), intent(in) :: shflx(pcols)              ! Kinematic surface heat flux [ unit ? ] 
    real(r8), intent(in) :: slslope(pcols,pver)       ! Slope of 'sl' in each layer [ J/kg/Pa ]
    real(r8), intent(in) :: qtslope(pcols,pver)       ! Slope of 'qt' in each layer [ kg/kg/Pa ]
    real(r8), intent(in) :: qrlin(pcols,pver)         ! Input grid-mean LW heating rate : [ K/s ] * cpair * dp = [ W/kg*Pa ]
    real(r8), intent(in) :: wsedl(pcols,pver)         ! Sedimentation velocity of liquid stratus cloud droplet [ m/s ]
    real(r8), intent(in) :: ustar(pcols)              ! Surface friction velocity [ m/s ]
    real(r8), intent(in) :: rrho(pcols)               ! 1./bottom mid-point density. Specific volume [ m3/kg ]
    real(r8), intent(in) :: kvf(pcols,pver+1)         ! Free atmosphere eddy diffusivity [ m2/s ]
    logical,  intent(in) :: wstarent                  ! Switch for choosing wstar3 entrainment parameterization
    real(r8), intent(in) :: minpblh(pcols)            ! Minimum PBL height based on surface stress [ m ]
    real(r8), intent(in) :: kvh_in(pcols,pver+1)      ! kvh saved from last timestep or last iterative step [ m2/s ] 
    real(r8), intent(in) :: kvm_in(pcols,pver+1)      ! kvm saved from last timestep or last iterative step [ m2/s ]
    real(r8), intent(in) :: cld(pcols,pver)           ! Stratus Cloud Fraction [ fraction ]

    ! ---------------- !
    ! Output variables !
    ! ---------------- !

    real(r8), intent(out) :: kvh(pcols,pver+1)        ! Eddy diffusivity for heat, moisture, and tracers [ m2/s ]
    real(r8), intent(out) :: kvm(pcols,pver+1)        ! Eddy diffusivity for momentum [ m2/s ]
    real(r8), intent(out) :: pblh(pcols)              ! PBL top height [ m ]
    real(r8), intent(out) :: pblhp(pcols)             ! PBL top height pressure [ Pa ]
    real(r8), intent(out) :: tpert(pcols)             ! Convective temperature excess [ K ]
    real(r8), intent(out) :: qpert(pcols)             ! Convective humidity excess [ kg/kg ]
    real(r8), intent(out) :: wpert(pcols)             ! Turbulent velocity excess [ m/s ]
    real(r8), intent(out) :: tkes(pcols)              ! TKE at surface [ m2/s2 ] 
    real(r8), intent(out) :: went(pcols)              ! Entrainment rate at the PBL top interface [ m/s ] 
    real(r8), intent(out) :: tke(pcols,pver+1)        ! Turbulent kinetic energy [ m2/s2 ], 'tkes' at surface, pver+1.
    real(r8), intent(out) :: bprod(pcols,pver+1)      ! Buoyancy production [ m2/s3 ],     'bflxs' at surface, pver+1.
    real(r8), intent(out) :: sprod(pcols,pver+1)      ! Shear production [ m2/s3 ], (ustar(i)**3)/(vk*z(i,pver))
                                                      ! at surface, pver+1.
    integer(i4), intent(out) :: turbtype(pcols,pver+1) ! Turbulence type at each interface:
                                                      ! 0. = Non turbulence interface
                                                      ! 1. = Stable turbulence interface
                                                      ! 2. = CL interior interface ( if bflxs > 0, surface is this )
                                                      ! 3. = Bottom external interface of CL
                                                      ! 4. = Top external interface of CL.
                                                      ! 5. = Double entraining CL external interface 
    real(r8), intent(out) :: sm_aw(pcols,pver+1)      ! Galperin instability function of momentum for use in the microphysics
                                                      ! [ no unit ]
    integer(i4), intent(out) :: ipbl(pcols)           ! If 1, PBL is CL, while if 0, PBL is STL.
    integer(i4), intent(out) :: kpblh(pcols)          ! Layer index containing PBL within or at the base interface
    real(r8), intent(out) :: wsed_CL(pcols,ncvmax)    ! Sedimentation velocity at the top of each CL [ m/s ]

    ! --------------------------- !
    ! Diagnostic output variables !
    ! --------------------------- !

    real(r8) :: kbase_o(pcols,ncvmax)                 ! Original external base interface index of CL just after 'exacol'
    real(r8) :: ktop_o(pcols,ncvmax)                  ! Original external top  interface index of CL just after 'exacol'
    real(r8) :: ncvfin_o(pcols)                       ! Original number of CLs just after 'exacol'
    real(r8) :: kbase_mg(pcols,ncvmax)                ! kbase  just after extending-merging (after 'zisocl') but without SRCL
    real(r8) :: ktop_mg(pcols,ncvmax)                 ! ktop   just after extending-merging (after 'zisocl') but without SRCL
    real(r8) :: ncvfin_mg(pcols)                      ! ncvfin just after extending-merging (after 'zisocl') but without SRCL
    real(r8) :: kbase_f(pcols,ncvmax)                 ! Final kbase  after adding SRCL
    real(r8) :: ktop_f(pcols,ncvmax)                  ! Final ktop   after adding SRCL
    real(r8) :: ncvfin_f(pcols)                       ! Final ncvfin after adding SRCL
    real(r8) :: wet_CL(pcols,ncvmax)                  ! Entrainment rate at the CL top [ m/s ] 
    real(r8) :: web_CL(pcols,ncvmax)                  ! Entrainment rate at the CL base [ m/s ]
    real(r8) :: jtbu_CL(pcols,ncvmax)                 ! Buoyancy jump across the CL top [ m/s2 ]  
    real(r8) :: jbbu_CL(pcols,ncvmax)                 ! Buoyancy jump across the CL base [ m/s2 ]  
    real(r8) :: evhc_CL(pcols,ncvmax)                 ! Evaporative enhancement factor at the CL top
    real(r8) :: jt2slv_CL(pcols,ncvmax)               ! Jump of slv ( across two layers ) at CL top for use only in evhc [ J/kg ]
    real(r8) :: n2ht_CL(pcols,ncvmax)                 ! n2 defined at the CL top  interface
                                                      ! but using sfuh(kt)   instead of sfi(kt) [ s-2 ]
    real(r8) :: n2hb_CL(pcols,ncvmax)                 ! n2 defined at the CL base interface
                                                      ! but using sflh(kb-1) instead of sfi(kb) [ s-2 ]
    real(r8) :: lwp_CL(pcols,ncvmax)                  ! LWP in the CL top layer [ kg/m2 ]
    real(r8) :: opt_depth_CL(pcols,ncvmax)            ! Optical depth of the CL top layer
    real(r8) :: radinvfrac_CL(pcols,ncvmax)           ! Fraction of LW radiative cooling confined in the top portion of CL
    real(r8) :: radf_CL(pcols,ncvmax)                 ! Buoyancy production at the CL top due to radiative cooling [ m2/s3 ]
    real(r8) :: wstar_CL(pcols,ncvmax)                ! Convective velocity of CL including entrainment contribution finally [ m/s ]
    real(r8) :: wstar3fact_CL(pcols,ncvmax)           ! "wstar3fact" of CL. Entrainment enhancement of wstar3 (inverse)

    real(r8) :: gh_a(pcols,pver+1)                    ! Half of normalized buoyancy production, -l2n2/2e. [ no unit ]
    real(r8) :: sh_a(pcols,pver+1)                    ! Galperin instability function of heat-moisture at all interfaces [ no unit ]
    real(r8) :: sm_a(pcols,pver+1)                    ! Galperin instability function of momentum      at all interfaces [ no unit ]
    real(r8) :: ri_a(pcols,pver+1)                    ! Interfacial Richardson number                  at all interfaces [ no unit ]

    real(r8) :: ebrk(pcols,ncvmax)                    ! Net CL mean TKE [ m2/s2 ]
    real(r8) :: wbrk(pcols,ncvmax)                    ! Net CL mean normalized TKE [ m2/s2 ]
    real(r8) :: lbrk(pcols,ncvmax)                    ! Net energetic integral thickness of CL [ m ]
    real(r8) :: ricl(pcols,ncvmax)                    ! Mean Richardson number of CL ( l2n2/l2s2 )
    real(r8) :: ghcl(pcols,ncvmax)                    ! Half of normalized buoyancy production of CL                 
    real(r8) :: shcl(pcols,ncvmax)                    ! Instability function of heat and moisture of CL
    real(r8) :: smcl(pcols,ncvmax)                    ! Instability function of momentum of CL

    real(r8) :: leng(pcols,pver+1)                    ! Turbulent length scale [ m ], 0 at the surface.
    real(r8) :: wcap(pcols,pver+1)                    ! Normalized TKE [m2/s2], 'tkes/b1' at the surface and 'tke/b1' at
                                                      ! the top/bottom entrainment interfaces of CL assuming no transport.
    ! ------------------------ !
    ! Local Internal Variables !
    ! ------------------------ !

    logical :: belongcv(pcols,pver+1)                 ! True for interfaces in a CL (both interior and exterior are included)
    logical :: in_CL                                  ! True if interfaces k,k+1 both in same CL.
    logical :: extend                                 ! True when CL is extended in zisocl
    logical :: extend_up                              ! True when CL is extended upward in zisocl
    logical :: extend_dn                              ! True when CL is extended downward in zisocl

    integer :: i                                      ! Longitude index
    integer :: k                                      ! Vertical index
    integer :: ks                                     ! Vertical index
    integer :: ncvfin(pcols)                          ! Total number of CL in column
    integer :: ncvf                                   ! Total number of CL in column prior to adding SRCL
    integer :: ncv                                    ! Index of current CL
    integer :: ncvnew                                 ! Index of added SRCL appended after regular CLs from 'zisocl'
    integer(i4) :: ncvsurf                            ! If nonzero, CL index based on surface
                                                      ! (usually 1, but can be > 1 when SRCL is based at sfc)
    integer(i4) :: srcl_status
    integer :: qrlzero_mode                          ! Encoded set_qrlzero for Codon helper
    integer :: cldeff_mode                           ! Encoded cldeff ramp choice for Codon helper
    integer :: tkes_mode                             ! Encoded choice_tkes for Codon helper
    integer :: use_kvf_mode                          ! Encoded use_kvf for Codon helper
    integer :: tunl_mode                              ! Encoded choice_tunl for Codon helper
    integer :: leng_mode                              ! Encoded choice_leng for Codon helper
    integer :: evhc_mode                              ! Encoded choice_evhc for Codon helper
    integer :: wstarent_mode                          ! Encoded wstarent for Codon helper
    integer :: sedfact_mode                           ! Encoded id_sedfact for Codon helper
    integer :: kbase(pcols,ncvmax)                    ! Vertical index of CL base interface
    integer :: ktop(pcols,ncvmax)                     ! Vertical index of CL top interface
    integer :: kb, kt                                 ! kbase and ktop for current CL
    integer :: ktblw                                  ! ktop of the CL located at just below the current CL

    integer  :: ktopbl(pcols)                         ! PBL top height or interface index 
    real(r8) :: bflxs(pcols)                          ! Surface buoyancy flux [ m2/s3 ]
    real(r8) :: rcap                                  ! 'tke/ebrk' at all interfaces of CL.
                                                      ! Set to 1 at the CL entrainment interfaces
    real(r8) :: jtzm                                  ! Interface layer thickness of CL top interface [ m ]
    real(r8) :: jtsl                                  ! Jump of s_l across CL top interface [ J/kg ]
    real(r8) :: jtqt                                  ! Jump of q_t across CL top interface [ kg/kg ]
    real(r8) :: jtbu                                  ! Jump of buoyancy across CL top interface [ m/s2 ]
    real(r8) :: jtu                                   ! Jump of u across CL top interface [ m/s ]
    real(r8) :: jtv                                   ! Jump of v across CL top interface [ m/s ]
    real(r8) :: jt2slv                                ! Jump of slv ( across two layers ) at CL top for use only in evhc [ J/kg ]
    real(r8) :: radf                                  ! Buoyancy production at the CL top due to radiative cooling [ m2/s3 ]
    real(r8) :: jbzm                                  ! Interface layer thickness of CL base interface [ m ]
    real(r8) :: jbsl                                  ! Jump of s_l across CL base interface [ J/kg ]
    real(r8) :: jbqt                                  ! Jump of q_t across CL top interface [ kg/kg ]
    real(r8) :: jbbu                                  ! Jump of buoyancy across CL base interface [ m/s2 ]
    real(r8) :: jbu                                   ! Jump of u across CL base interface [ m/s ]
    real(r8) :: jbv                                   ! Jump of v across CL base interface [ m/s ]
    real(r8) :: ch                                    ! Buoyancy coefficients defined at the CL top and base interfaces
                                                      ! using CL internal
    real(r8) :: cm                                    ! sfuh(kt) and sflh(kb-1) instead of sfi(kt) and sfi(kb), respectively.
                                                      ! These are used for entrainment calculation at CL external interfaces
                                                      ! and SRCL identification.
    real(r8) :: n2ht                                  ! n2 defined at the CL top  interface
                                                      ! but using sfuh(kt)   instead of sfi(kt) [ s-2 ]
    real(r8) :: n2hb                                  ! n2 defined at the CL base interface
                                                      ! but using sflh(kb-1) instead of sfi(kb) [ s-2 ]
    real(r8) :: n2htSRCL                              ! n2 defined at the upper-half layer of SRCL.
                                                      ! This is used only for identifying SRCL.
                                                      ! n2htSRCL use SRCL internal slope sl and qt
                                                      ! as well as sfuh(kt) instead of sfi(kt) [ s-2 ]
    real(r8) :: gh                                    ! Half of normalized buoyancy production ( -l2n2/2e ) [ no unit ]
    real(r8) :: sh                                    ! Galperin instability function for heat and moisture
    real(r8) :: sm                                    ! Galperin instability function for momentum
    real(r8) :: lbulk                                 ! Depth of turbulent layer, Master length scale (not energetic length)
    real(r8) :: dzht                                  ! Thickness of top    half-layer [ m ]
    real(r8) :: dzhb                                  ! Thickness of bottom half-layer [ m ]
    real(r8) :: rootp                                 ! Sqrt(net CL-mean TKE including entrainment contribution) [ m/s ]     
    real(r8) :: evhc                                  ! Evaporative enhancement factor: (1+E)
                                                      ! with E = evap. cool. efficiency [ no unit ]
    real(r8) :: kentr                                 ! Effective entrainment diffusivity 'wet*dz', 'web*dz' [ m2/s ]
    real(r8) :: lwp                                   ! Liquid water path in the layer kt [ kg/m2 ]
    real(r8) :: opt_depth                             ! Optical depth of the layer kt [ no unit ]
    real(r8) :: radinvfrac                            ! Fraction of LW cooling in the layer kt
                                                      ! concentrated at the CL top [ no unit ]
    real(r8) :: wet                                   ! CL top entrainment rate [ m/s ]
    real(r8) :: web                                   ! CL bot entrainment rate [ m/s ]. Set to zero if CL is based at surface.
    real(r8) :: vyt                                   ! n2ht/n2 at the CL top  interface
    real(r8) :: vyb                                   ! n2hb/n2 at the CL base interface
    real(r8) :: vut                                   ! Inverse Ri (=s2/n2) at the CL top  interface
    real(r8) :: vub                                   ! Inverse Ri (=s2/n2) at the CL base interface
    real(r8) :: fact                                  ! Factor relating TKE generation to entrainment [ no unit ]
    real(r8) :: trma                                  ! Intermediate variables used for solving quadratic ( for gh from ri )
    real(r8) :: trmb                                  ! and cubic equations ( for ebrk: the net CL mean TKE )
    real(r8) :: trmc                                  !
    real(r8) :: trmp                                  !
    real(r8) :: trmq                                  !
    real(r8) :: qq                                    ! 
    real(r8) :: det                                   !
    real(r8) :: gg                                    ! Intermediate variable used for calculating stability functions of
                                                      ! SRCL or SBCL based at the surface with bflxs > 0.
    real(r8) :: dzhb5                                 ! Half thickness of the bottom-most layer of current CL regime
    real(r8) :: dzht5                                 ! Half thickness of the top-most layer of adjacent CL regime
                                                      ! just below current CL
    real(r8) :: qrlw(pcols,pver)                      ! Local grid-mean LW heating rate : [K/s] * cpair * dp = [ W/kg*Pa ]

    real(r8) :: cldeff(pcols,pver)                    ! Effective stratus fraction
    real(r8) :: qleff                                 ! Used for computing evhc
    real(r8) :: tunlramp                              ! Ramping tunl
    real(r8) :: leng_imsi                             ! For Kv = max(Kv_STL, Kv_entrain)
    real(r8) :: tke_imsi                              !
    real(r8) :: kvh_imsi                              !
    real(r8) :: kvm_imsi                              !
    real(r8) :: alph4exs                              ! For extended stability function in the stable regime
    real(r8) :: ghmin                                 !   

    real(r8) :: sedfact                               ! For 'sedimentation-entrainment feedback' 

    ! Local variables specific for 'wstar' entrainment closure

    real(r8) :: cet                                   ! Proportionality coefficient between wet and wstar3
    real(r8) :: ceb                                   ! Proportionality coefficient between web and wstar3
    real(r8) :: wstar                                 ! Convective velocity for CL [ m/s ]
    real(r8) :: wstar3                                ! Cubed convective velocity for CL [ m3/s3 ]
    real(r8) :: wstar3fact                            ! 1/(relative change of wstar^3 by entrainment)
    real(r8) :: rmin                                  ! sqrt(p)
    real(r8) :: fmin                                  ! f(rmin), where f(r) = r^3 - 3*p*r - 2q
    real(r8) :: rcrit                                 ! ccrit*wstar
    real(r8) :: fcrit                                 ! f(rcrit)
    logical     noroot                                ! True if f(r) has no root r > rcrit

    !-----------------------!
    ! Start of Main Program !
    !-----------------------!
    
    ! Option: Turn-off LW radiative-turbulence interaction in PBL scheme
    !         by setting qrlw = 0.  Logical parameter 'set_qrlzero'  was
    !         defined in the first part of 'eddy_diff.F90' module. 

    ! For an extended stability function in the stable regime, re-define
    ! alph4exe and ghmin. This is for future work.

    call eddy_diff_caleddy_stable_config(ricrit, b1, alph2, alph3, alph4, alph5, alph4exs, ghmin)

    tunl_mode = 0
    if( choice_tunl .eq. 'rampcl' ) then
        tunl_mode = 1
    elseif( choice_tunl .eq. 'rampsl' ) then
        tunl_mode = 2
    end if

    leng_mode = 1
    if( choice_leng .eq. 'origin' ) then
        leng_mode = 0
    end if

    evhc_mode = 2
    if( choice_evhc .eq. 'orig' ) then
        evhc_mode = 0
    elseif( choice_evhc .eq. 'ramp' ) then
        evhc_mode = 1
    end if

    wstarent_mode = 0
    if (wstarent) then
        wstarent_mode = 1
    end if

    sedfact_mode = 0
    if (id_sedfact) then
        sedfact_mode = 1
    end if

    qrlzero_mode = 0
    if (set_qrlzero) then
        qrlzero_mode = 1
    end if

    cldeff_mode = 0
    if (choice_evhc .eq. 'ramp' .or. choice_radf .eq. 'ramp') then
        cldeff_mode = 1
    end if

    tkes_mode = 0
    if (choice_tkes .eq. 'ibprod') then
        tkes_mode = 1
    end if

    use_kvf_mode = 0
    if (use_kvf) then
        use_kvf_mode = 1
    end if

    call eddy_diff_caleddy_init(ncol, pcols, pver, qrlzero_mode, cldeff_mode, tkes_mode, use_kvf_mode, qmin, vk, ql, qrlin, &
         cld, kvf, kvh_in, kvm_in, n2, s2, shflx, qflx, rrho, ustar, z, chu, chs, cmu, cms, sflh, qrlw, cldeff, kvh, kvm, &
         bflxs, bprod, sprod, wcap, leng, tke, turbtype)

    !
    ! Initialization of Diagnostic Output
    !

    call eddy_diff_caleddy_diaginit(ncol, pcols, pver, ncvmax, went, wet_CL, web_CL, jtbu_CL, jbbu_CL, evhc_CL, &
         jt2slv_CL, n2ht_CL, n2hb_CL, lwp_CL, opt_depth_CL, radinvfrac_CL, radf_CL, wstar_CL, wstar3fact_CL, ricl, ghcl, &
         shcl, smcl, ebrk, wbrk, lbrk, gh_a, sh_a, sm_a, ri_a, sm_aw, ipbl, kpblh, wsed_CL)

    ! Initially identify CL regimes in 'exacol'
    !    ktop  : Interface index of the CL top  external interface
    !    kbase : Interface index of the CL base external interface
    !    ncvfin: Number of total CLs
    ! Note that if surface buoyancy flux is positive ( bflxs = bprod(i,pver+1) > 0 ),
    ! surface interface is identified as an internal interface of CL. However, even
    ! though bflxs <= 0, if 'pver' interface is a CL internal interface (ri(pver)<0),
    ! surface interface is identified as an external interface of CL. If bflxs =< 0 
    ! and ri(pver) >= 0, then surface interface is identified as a stable turbulent
    ! intereface (STL) as shown at the end of 'caleddy'. Even though a 'minpblh' is
    ! passed into 'exacol', it is not used in the 'exacol'.

    call exacol( pcols, pver, ncol, ri, bflxs, minpblh, zi, ktop, kbase, ncvfin )

    ! Diagnostic output of CL interface indices before performing 'extending-merging'
    ! of CL regimes in 'zisocl'
    do i = 1, ncol
       call eddy_diff_caleddy_regime_diag(i, pcols, ncvmax, kbase, ktop, ncvfin, kbase_o, ktop_o, ncvfin_o)
    end do

    ! ----------------------------------- !
    ! Perform calculation for each column !
    ! ----------------------------------- !

    do i = 1, ncol

       ! Define Surface Interfacial Layer TKE, 'tkes'.
       ! In the current code, 'tkes' is used as representing TKE of surface interfacial
       ! layer (low half-layer of surface-based grid layer). In the code, when bflxs>0,
       ! surface interfacial layer is assumed to be energetically  coupled to the other
       ! parts of the CL regime based at the surface. In this sense, it is conceptually
       ! more reasonable to include both 'bprod' and 'sprod' in the definition of 'tkes'.
       ! Since 'tkes' cannot be negative, it is lower bounded by small positive number. 
       ! Note that inclusion of 'bprod' in the definition of 'tkes' may increase 'ebrk'
       ! and 'wstar3', and eventually, 'wet' at the CL top, especially when 'bflxs>0'.
       ! This might help to solve the problem of too shallow PBLH over the overcast Sc
       ! regime. If I want to exclude 'bprod(i,pver+1)' in calculating 'tkes' even when
       ! bflxs > 0, all I should to do is to set 'bprod(i,pver+1) = 0' in the above 
       ! initialization 'do' loop (explained above), NOT changing the formulation of
       ! tkes(i) in the below block. This is because for consistent treatment in the 
       ! other parts of the code also.
  
     ! tkes(i) = (b1*vk*z(i,pver)*sprod(i,pver+1))**(2._r8/3._r8)
       call eddy_diff_caleddy_surface_tke(i, pcols, pver, b1, vk, tkemax, z, bprod, sprod, tkes, tke, wcap)

       ! Extend and merge the initially identified CLs, relabel the CLs, and calculate
       ! CL internal mean energetics and stability functions in 'zisocl'. 
       ! The CL nearest to the surface is CL(1) and the CL index, ncv, increases 
       ! with height. The following outputs are from 'zisocl'. Here, the dimension
       ! of below outputs are (pcols,ncvmax) (except the 'ncvfin(pcols)' and 
       ! 'belongcv(pcols,pver+1)) and 'ncv' goes from 1 to 'ncvfin'. 
       ! For 'ncv = ncvfin+1, ncvmax', below output are already initialized to be zero. 
       !      ncvfin       : Total number of CLs
       !      kbase(ncv)   : Base external interface index of CL
       !      ktop         : Top  external interface index of CL
       !      belongcv     : True if the interface (either internal or external) is CL  
       !      ricl         : Mean Richardson number of internal CL
       !      ghcl         : Normalized buoyancy production '-l2n2/2e' [no unit] of internal CL
       !      shcl         : Galperin instability function of heat-moisture of internal CL
       !      smcl         : Galperin instability function of momentum of internal CL
       !      lbrk, <l>int : Thickness of (energetically) internal CL (lint, [m])
       !      wbrk, <W>int : Mean normalized TKE of internal CL  ([m2/s2])
       !      ebrk, <e>int : Mean TKE of internal CL (b1*wbrk,[m2/s2])
       ! The ncvsurf is an identifier saying which CL regime is based at the surface.
       ! If 'ncvsurf=1', then the first CL regime is based at the surface. If surface
       ! interface is not a part of CL (neither internal nor external), 'ncvsurf = 0'.
       ! After identifying and including SRCLs into the normal CL regimes (where newly
       ! identified SRCLs are simply appended to the normal CL regimes using regime 
       ! indices of 'ncvfin+1','ncvfin+2' (as will be shown in the below SRCL part),..
       ! where 'ncvfin' is the final CL regime index produced after extending-merging 
       ! in 'zisocl' but before adding SRCLs), if any newly identified SRCL (e.g., 
       ! 'ncvfin+1') is based at surface, then 'ncvsurf = ncvfin+1'. Thus 'ncvsurf' can
       ! be 0, 1, or >1. 'ncvsurf' can be a useful diagnostic output.   

       ncvsurf = 0
       if( ncvfin(i) .gt. 0 ) then 
           call zisocl( pcols  , pver     , i        ,           &
                        z      , zi       , n2       , s2      , & 
                        bprod  , sprod    , bflxs    , tkes    , &
                        ncvfin , kbase    , ktop     , belongcv, &
                        ricl   , ghcl     , shcl     , smcl    , & 
                        lbrk   , wbrk     , ebrk     ,           & 
                        extend , extend_up, extend_dn )
           if( kbase(i,1) .eq. pver + 1 ) ncvsurf = 1
       else
           belongcv(i,:) = .false.
       endif

       ! Diagnostic output after finishing extending-merging process in 'zisocl'
       ! Since we are adding SRCL additionally, we need to print out these here.

       call eddy_diff_caleddy_regime_diag(i, pcols, ncvmax, kbase, ktop, ncvfin, kbase_mg, ktop_mg, ncvfin_mg)

       ! ----------------------- !
       ! Identification of SRCLs !
       ! ----------------------- !

     ! Modification : This cannot identify the 'cirrus' layer due to the condition of
     !                ql(i,k) .gt. qmin. This should be modified in future to identify
     !                a single thin cirrus layer.  
     !                Instead of ql, we may use cldn in future, including ice 
     !                contribution.

       ! ------------------------------------------------------------------------------ !
       ! Find single-layer radiatively-driven cloud-topped convective layers (SRCLs).   !
       ! SRCLs extend through a single model layer k, with entrainment at the top and   !
       ! bottom interfaces, unless bottom interface is the surface.                     !
       ! The conditions for an SRCL is identified are:                                  ! 
       !                                                                                !
       !   1. Cloud in the layer, k : ql(i,k) .gt. qmin = 1.e-5 [ kg/kg ]               !
       !   2. No cloud in the above layer (else assuming that some fraction of the LW   !
       !      flux divergence in layer k is concentrated at just below top interface    !
       !      of layer k is invalid). Then, this condition might be sensitive to the    !
       !      vertical resolution of grid.                                              !
       !   3. LW radiative cooling (SW heating is assumed uniformly distributed through !
       !      layer k, so not relevant to buoyancy production) in the layer k. However, !
       !      SW production might also contribute, which may be considered in a future. !
       !   4. Internal stratification 'n2ht' of upper-half layer should be unstable.    !
       !      The 'n2ht' is pure internal stratification of upper half layer, obtained  !
       !      using internal slopes of sl, qt in layer k (in contrast to conventional   !
       !      interfacial slope) and saturation fraction in the upper-half layer,       !
       !      sfuh(k) (in contrast to sfi(k)).                                          !
       !   5. Top and bottom interfaces not both in the same existing convective layer. !
       !      If SRCL is within the previouisly identified CL regimes, we don't define  !
       !      a new SRCL.                                                               !
       !   6. k >= ntop_turb + 1 = 2                                                    !
       !   7. Ri at the top interface > ricrit = 0.19 (otherwise turbulent mixing will  !
       !      broadly distribute the cloud top in the vertical, preventing localized    !
       !      radiative destabilization at the top interface).                          !
       !                                                                                !
       ! Note if 'k = pver', it identifies a surface-based single fog layer, possibly,  !
       ! warm advection fog. Note also the CL regime index of SRCLs itself increases    !
       ! with height similar to the regular CLs indices identified from 'zisocl'.       !
       ! ------------------------------------------------------------------------------ !

       call eddy_diff_caleddy_srcl(choice_SRCL, i, pcols, pver, ncvmax, ntop_turb, nbot_turb, qmin, ricrit, b1, vk, alph1, &
            alph2, alph3, alph4exs, alph5, ghmin, ql, qrlw, ri, sfuh, chu, chs, cmu, cms, slslope, qtslope, z, bflxs, &
            tkes, bprod, sprod, ncvfin, kbase, ktop, belongcv, ricl, ghcl, shcl, smcl, lbrk, wbrk, ebrk, ncvsurf, srcl_status)

       if (srcl_status .ne. 0_i4) then
           write(iulog,*) 'Major mistake in SRCL: bflxs > 0 for surface-based SRCL'
           write(iulog,*) 'bflxs = ', bflxs(i)
           write(iulog,*) 'ncvfin_o = ', ncvfin_o(i)
           write(iulog,*) 'ncvfin_mg = ', ncvfin_mg(i)
           do ks = 1, ncvmax
              write(iulog,*) 'ncv =', ks, ' ', kbase_o(i,ks), ktop_o(i,ks), kbase_mg(i,ks), ktop_mg(i,ks)
           end do
           call endrun('CALEDDY: Major mistake in SRCL: bflxs > 0 for surface-based SRCL')
       end if

       ! -------------------------------------------------------------------------- !
       ! Up to this point, we identified all kinds of CL regimes :                  !
       !   1. A SBCL. By construction, 'bflxs > 0' for SBCL.                        !
       !   2. Surface-based CL with multiple layers and 'bflxs =< 0'                !
       !   3. Surface-based CL with multiple layers and 'bflxs > 0'                 !
       !   4. Regular elevated CL with two entraining interfaces                    ! 
       !   5. SRCLs. If SRCL is based at surface, it will be bflxs < 0.             !
       ! '1-4' were identified from 'zisocl' while '5' were identified separately   !
       ! after performing 'zisocl'. CL regime index of '1-4' increases with height  !
       ! ( e.g., CL = 1 is the CL regime nearest to the surface ) while CL regime   !
       ! index of SRCL is simply appended after the final index of CL regimes from  !
       ! 'zisocl'. However, CL regime indices of SRCLs itself increases with height !
       ! when there are multiple SRCLs, similar to the regular CLs from 'zisocl'.   !
       ! -------------------------------------------------------------------------- !

       ! Diagnostic output of final CL regimes indices
       
       call eddy_diff_caleddy_regime_diag(i, pcols, ncvmax, kbase, ktop, ncvfin, kbase_f, ktop_f, ncvfin_f)

       ! --------------------------------------------------------------------- !
       ! Compute radf for each CL in column by calling subroutine compute_radf !
       ! --------------------------------------------------------------------- !
       call eddy_diff_compute_radf(choice_radf, i, pcols, pver, ncvmax, ncvfin, ktop, qmin, ql, pi, qrlw, g, cldeff, zi, &
            chs, lwp_CL, opt_depth_CL, radinvfrac_CL, radf_CL)

       ! ---------------------------------------- !
       ! Perform do loop for individual CL regime !
       ! ---------------------------------------- ! -------------------------------- !
       ! For individual CLs, compute                                                 !
       !   1. Entrainment rates at the CL top and (if any) base interfaces using     !
       !      appropriate entrainment closure (current code use 'wstar' closure).    !
       !   2. Net CL mean (i.e., including entrainment contribution) TKE (ebrk)      !
       !      and normalized TKE (wbrk).                                             ! 
       !   3. TKE (tke) and normalized TKE (wcap) profiles at all CL interfaces.     !
       !   4. ( kvm, kvh ) profiles at all CL interfaces.                            !
       !   5. ( bprod, sprod ) profiles at all CL interfaces.                        !
       ! Also calculate                                                              !
       !   1. PBL height as the top external interface of surface-based CL, if any.  !
       !   2. Characteristic excesses of convective 'updraft velocity (wpert)',      !
       !      'temperature (tpert)', and 'moisture (qpert)' in the surface-based CL, !
       !      if any, for use in the separate convection scheme.                     ! 
       ! If there is no surface-based CL, 'PBL height' and 'convective excesses' are !
       ! calculated later from surface-based STL (Stable Turbulent Layer) properties.!
       ! --------------------------------------------------------------------------- !

       call eddy_diff_caleddy_closure(i, pcols, pver, ncvmax, tunl_mode, leng_mode, evhc_mode, wstarent_mode, sedfact_mode, &
            tunl, ctunl, cleng, lbulk_max, tkemax, b1, ae, alph1, a1l, a1i, ccrit, wstar3factcrit, ntzero, onet, rcapmin, &
            rcapmax, wfac, wpertmin, tfac, qmin, g, vk, cpair, latvap, a2l, a3l, jbumin, evhcmax, ased, ql, slv, sl, qt, u, &
            v, pi, zi, z, n2, s2, shflx, qflx, rrho, sfuh, sflh, chu, chs, cmu, cms, cldeff, bflxs, bprod, sprod, wsedl, &
            ncvfin, kbase, ktop, belongcv, lbrk, ebrk, wbrk, ricl, shcl, smcl, radf_CL, wsed_CL, leng_max, wet_CL, web_CL, &
            jtbu_CL, jbbu_CL, evhc_CL, jt2slv_CL, n2ht_CL, n2hb_CL, wstar_CL, wstar3fact_CL, leng, wcap, tke, kvh, kvm, &
            turbtype, sm_aw, pblh, pblhp, wpert, tpert, qpert, ipbl, kpblh, went, ncvsurf)

       call eddy_diff_caleddy_stl(i, pcols, pver, ncvmax, tunl_mode, leng_mode, ricrit, tunl, ctunl, cleng, lbulk_max, &
            tkemax, b1, ae, alph1, alph2, alph3, alph4exs, alph5, ghmin, vk, fak, cpair, ri, z, zi, pi, n2, s2, shflx, &
            qflx, rrho, ustar, leng_max, ncvfin, ktop, kbase, kvh, kvm, leng, tke, wcap, bprod, sprod, turbtype, sm_aw, &
            pblh, pblhp, wpert, tpert, qpert, ipbl, kpblh)

       ! As an option, we can impose a certain minimum back-ground diffusivity.

       ! do k = 1, pver+1
       !    kvh(i,k) = max(0.01_r8,kvh(i,k))
       !    kvm(i,k) = max(0.01_r8,kvm(i,k))
       ! enddo
 
       ! --------------------------------------------------------------------- !
       ! Diagnostic Output                                                     !
       ! Just for diagnostic purpose, calculate stability functions at  each   !
       ! interface including surface. Instead of assuming neutral stability,   !
       ! explicitly calculate stability functions using an reverse procedure   !
       ! starting from tkes(i) similar to the case of SRCL and SBCL in zisocl. !
       ! Note that it is possible to calculate stability functions even when   !
       ! bflxs < 0. Note that this inverse method allows us to define Ri even  !
       ! at the surface. Note also tkes(i) and sprod(i,pver+1) are always      !
       ! positive values by limiters (e.g., ustar_min = 0.01).                 !
       ! Dec.12.2006 : Also just for diagnostic output, re-set                 !
       ! 'bprod(i,pver+1)= bflxs(i)' here. Note that this setting does not     !
       ! influence numerical calculation at all - it is just for diagnostic    !
       ! output.                                                               !
       ! --------------------------------------------------------------------- !

       call eddy_diff_caleddy_diag(i, pcols, pver, ricrit, tkes, b1, alph1, alph2, alph3, alph4, alph4exs, alph5, ghmin, &
            vk, z, ri, bflxs, bprod, sprod, gh_a, sh_a, sm_a, ri_a, sm_aw)

    end do   ! End of column index loop, i

    return

    end subroutine caleddy

    !============================================================================== !
    !                                                                               !
    !============================================================================== !

  subroutine exacol_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (exacol_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('EDDY_DIFF_EXACOL_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_exacol_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_exacol_impl = .false.
    end if

    exacol_impl_selected = .true.

    if (masterproc) then
       if (use_native_exacol_impl) then
          write(iulog,*) 'eddy_diff_exacol implementation = native'
       else
          write(iulog,*) 'eddy_diff_exacol implementation = codon'
       end if
    end if

  end subroutine exacol_select_impl

  !============================================================================== !
  !                                                                               !
  !============================================================================== !

    subroutine exacol( pcols, pver, ncol, ri, bflxs, minpblh, zi, ktop, kbase, ncvfin )

    use iso_c_binding, only: c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: pcols
    integer, intent(in) :: pver
    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: ri(pcols,pver)
    real(r8), target, intent(in) :: bflxs(pcols)
    real(r8), intent(in) :: minpblh(pcols)
    real(r8), intent(in) :: zi(pcols,pver+1)
    integer(i4), target, intent(out) :: ktop(pcols,ncvmax)
    integer(i4), target, intent(out) :: kbase(pcols,ncvmax)
    integer(i4), target, intent(out) :: ncvfin(pcols)

    interface
       subroutine eddy_diff_exacol_codon(ncol_c, pcols_c, pver_c, ncvmax_c, ntop_turb_c, ri_p, bflxs_p, ktop_p, kbase_p, &
            ncvfin_p) bind(c, name="eddy_diff_exacol_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, ncvmax_c, ntop_turb_c
         type(c_ptr), value :: ri_p, bflxs_p, ktop_p, kbase_p, ncvfin_p
       end subroutine eddy_diff_exacol_codon
    end interface

    call exacol_select_impl()

    if (use_native_exacol_impl) then
       call exacol_native(pcols, pver, ncol, ri, bflxs, minpblh, zi, ktop, kbase, ncvfin)
       return
    end if

    call eddy_diff_exacol_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(ncvmax, c_int64_t), &
         int(ntop_turb, c_int64_t), c_loc(ri), c_loc(bflxs), c_loc(ktop), c_loc(kbase), c_loc(ncvfin))

    return

    end subroutine exacol

    !============================================================================== !
    !                                                                               !
    !============================================================================== !

    subroutine exacol_native( pcols, pver, ncol, ri, bflxs, minpblh, zi, ktop, kbase, ncvfin )

    ! ---------------------------------------------------------------------------- !
    ! Object : Find unstable CL regimes and determine the indices                  !
    !          kbase, ktop which delimit these unstable layers :                   !
    !          ri(kbase) > 0 and ri(ktop) > 0, but ri(k) < 0 for ktop < k < kbase. ! 
    ! Author : Chris  Bretherton 08/2000,                                          !
    !          Sungsu Park       08/2006, 11/2008                                  !
    !----------------------------------------------------------------------------- !

    implicit none

    ! --------------- !
    ! Input variables !
    ! --------------- !

    integer,  intent(in) :: pcols                  ! Number of atmospheric columns   
    integer,  intent(in) :: pver                   ! Number of atmospheric vertical layers   
    integer,  intent(in) :: ncol                   ! Number of atmospheric columns   

    real(r8), intent(in) :: ri(pcols,pver)         ! Moist gradient Richardson no.
    real(r8), intent(in) :: bflxs(pcols)           ! Buoyancy flux at surface
    real(r8), intent(in) :: minpblh(pcols)         ! Minimum PBL height based on surface stress
    real(r8), intent(in) :: zi(pcols,pver+1)       ! Interface heights

    ! ---------------- !
    ! Output variables !      
    ! ---------------- !

    integer(i4), intent(out) :: kbase(pcols,ncvmax)    ! External interface index of CL base
    integer(i4), intent(out) :: ktop(pcols,ncvmax)     ! External interface index of CL top
    integer(i4), intent(out) :: ncvfin(pcols)          ! Total number of CLs

    ! --------------- !
    ! Local variables !
    ! --------------- !

    integer              :: i
    integer              :: k
    integer              :: ncv
    real(r8)             :: rimaxentr
    real(r8)             :: riex(pver+1)           ! Column Ri profile extended to surface

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !

    do i = 1, ncol
       ncvfin(i) = 0
       do ncv = 1, ncvmax
          ktop(i,ncv)  = 0
          kbase(i,ncv) = 0
       end do
    end do

    ! ------------------------------------------------------ !
    ! Find CL regimes starting from the surface going upward !
    ! ------------------------------------------------------ !
    
    rimaxentr = 0._r8   
    
    do i = 1, ncol

       riex(2:pver) = ri(i,2:pver)

       ! Below allows consistent treatment of surface and other interfaces.
       ! Simply, if surface buoyancy flux is positive, Ri of surface is set to be negative.

       riex(pver+1) = rimaxentr - bflxs(i) 

       ncv = 0
       k   = pver + 1 ! Work upward from surface interface

       do while ( k .gt. ntop_turb + 1 )

        ! Below means that if 'bflxs > 0' (do not contain '=' sign), surface
        ! interface is energetically interior surface. 
       
          if( riex(k) .lt. rimaxentr ) then 

              ! Identify a new CL

              ncv = ncv + 1

              ! First define 'kbase' as the first interface below the lower-most unstable interface
              ! Thus, Richardson number at 'kbase' is positive.

              kbase(i,ncv) = min(k+1,pver+1)

              ! Decrement k until top unstable level

              do while( riex(k) .lt. rimaxentr .and. k .gt. ntop_turb + 1 )
                 k = k - 1
              end do

              ! ktop is the first interface above upper-most unstable interface
              ! Thus, Richardson number at 'ktop' is positive. 

              ktop(i,ncv) = k
             
          else

              ! Search upward for a CL.

              k = k - 1

          end if

       end do ! End of CL regime finding for each atmospheric column

       ncvfin(i) = ncv    

    end do  ! End of atmospheric column do loop

    return 

    end subroutine exacol_native

    !============================================================================== !
    !                                                                               !
    !============================================================================== !
    
    subroutine zisocl( pcols  , pver  , long ,                                 & 
                       z      , zi    , n2   ,  s2      ,                      & 
                       bprod  , sprod , bflxs,  tkes    ,                      & 
                       ncvfin , kbase , ktop ,  belongcv,                      & 
                       ricl   , ghcl  , shcl ,  smcl    ,                      &
                       lbrk   , wbrk  , ebrk ,  extend  , extend_up, extend_dn )

    !------------------------------------------------------------------------ !
    ! Object : This 'zisocl' vertically extends original CLs identified from  !
    !          'exacol' using a merging test based on either 'wint' or 'l2n2' !
    !          and identify new CL regimes. Similar to the case of 'exacol',  !
    !          CL regime index increases with height.  After identifying new  !
    !          CL regimes ( kbase, ktop, ncvfin ),calculate CL internal mean  !
    !          energetics (lbrk : energetic thickness integral, wbrk, ebrk )  !
    !          and stability functions (ricl, ghcl, shcl, smcl) by including  !
    !          surface interfacial layer contribution when bflxs > 0.   Note  !
    !          that there are two options in the treatment of the energetics  !
    !          of surface interfacial layer (use_dw_surf= 'true' or 'false')  !
    ! Author : Sungsu Park 08/2006, 11/2008                                   !
    !------------------------------------------------------------------------ !

    implicit none

    ! --------------- !    
    ! Input variables !
    ! --------------- !

    integer,  intent(in)   :: long                    ! Longitude of the column
    integer,  intent(in)   :: pcols                   ! Number of atmospheric columns   
    integer,  intent(in)   :: pver                    ! Number of atmospheric vertical layers   
    real(r8), intent(in)   :: z(pcols, pver)          ! Layer mid-point height [ m ]
    real(r8), intent(in)   :: zi(pcols, pver+1)       ! Interface height [ m ]
    real(r8), intent(in)   :: n2(pcols, pver)         ! Buoyancy frequency at interfaces except surface [ s-2 ]
    real(r8), intent(in)   :: s2(pcols, pver)         ! Shear frequency at interfaces except surface [ s-2 ]
    real(r8), intent(in)   :: bprod(pcols,pver+1)     ! Buoyancy production [ m2/s3 ]. bprod(i,pver+1) = bflxs 
    real(r8), intent(in)   :: sprod(pcols,pver+1)     ! Shear production [ m2/s3 ]. sprod(i,pver+1) = usta**3/(vk*z(i,pver))
    real(r8), intent(in)   :: bflxs(pcols)            ! Surface buoyancy flux [ m2/s3 ]. bprod(i,pver+1) = bflxs 
    real(r8), intent(in)   :: tkes(pcols)             ! TKE at the surface [ s2/s2 ]

    ! ---------------------- !
    ! Input/output variables !
    ! ---------------------- !

    integer, intent(inout) :: kbase(pcols,ncvmax)     ! Base external interface index of CL
    integer, intent(inout) :: ktop(pcols,ncvmax)      ! Top external interface index of CL
    integer, intent(inout) :: ncvfin(pcols)           ! Total number of CLs

    ! ---------------- !
    ! Output variables !
    ! ---------------- !

    logical,  intent(out) :: belongcv(pcols,pver+1)   ! True if interface is in a CL ( either internal or external )
    real(r8), intent(out) :: ricl(pcols,ncvmax)       ! Mean Richardson number of internal CL
    real(r8), intent(out) :: ghcl(pcols,ncvmax)       ! Half of normalized buoyancy production of internal CL
    real(r8), intent(out) :: shcl(pcols,ncvmax)       ! Galperin instability function of heat-moisture of internal CL
    real(r8), intent(out) :: smcl(pcols,ncvmax)       ! Galperin instability function of momentum of internal CL
    real(r8), intent(out) :: lbrk(pcols,ncvmax)       ! Thickness of (energetically) internal CL ( lint, [m] )
    real(r8), intent(out) :: wbrk(pcols,ncvmax)       ! Mean normalized TKE of internal CL  [ m2/s2 ]
    real(r8), intent(out) :: ebrk(pcols,ncvmax)       ! Mean TKE of internal CL ( b1*wbrk, [m2/s2] )

    ! ------------------ !
    ! Internal variables !
    ! ------------------ !

    logical               :: extend                   ! True when CL is extended in zisocl
    logical               :: extend_up                ! True when CL is extended upward in zisocl
    logical               :: extend_dn                ! True when CL is extended downward in zisocl
    logical               :: bottom                   ! True when CL base is at surface ( kb = pver + 1 )

    integer               :: i                        ! Local index for the longitude
    integer               :: ncv                      ! CL Index increasing with height
    integer               :: incv
    integer               :: k
    integer               :: kb                       ! Local index for kbase
    integer               :: kt                       ! Local index for ktop
    integer               :: kb_before                ! Previous value of kb before helper update
    integer               :: kt_before                ! Previous value of kt before helper update
    integer               :: ncvinit                  ! Value of ncv at routine entrance 
    integer               :: cntu                     ! Number of merged CLs during upward   extension of individual CL
    integer               :: cntd                     ! Number of merged CLs during downward extension of individual CL
    integer               :: kbinc                    ! Index for incorporating underlying CL
    integer               :: ktinc                    ! Index for incorporating  overlying CL
    integer               :: kb_is_surface_mode
    integer               :: use_dw_surf_mode
    integer               :: choice_tkes_ebprod_mode
    integer               :: upward_state_status
    integer               :: downward_state_status

    real(r8)              :: wint                     ! Normalized TKE of internal CL
    real(r8)              :: dwinc                    ! Normalized TKE of CL external interfaces
    real(r8)              :: dw_surf                  ! Normalized TKE of surface interfacial layer
    real(r8)              :: dzinc
    real(r8)              :: gh
    real(r8)              :: sh
    real(r8)              :: sm
    real(r8)              :: l2n2                     ! Vertical integral of 'l^2N^2' over CL. Include thickness product
    real(r8)              :: l2s2                     ! Vertical integral of 'l^2S^2' over CL. Include thickness product
    real(r8)              :: dl2n2                    ! Vertical integration of 'l^2*N^2' of CL external interfaces
    real(r8)              :: dl2s2                    ! Vertical integration of 'l^2*S^2' of CL external interfaces
    real(r8)              :: dl2n2_surf               ! 'dl2n2' defined in the surface interfacial layer
    real(r8)              :: dl2s2_surf               ! 'dl2s2' defined in the surface interfacial layer  
    real(r8)              :: lint                     ! Thickness of (energetically) internal CL
    real(r8)              :: dlint                    ! Interfacial layer thickness of CL external interfaces
    real(r8)              :: dlint_surf               ! Surface interfacial layer thickness 
    real(r8)              :: lbulk                    ! Master Length Scale : Whole CL thickness from top to base external interface
    real(r8)              :: ricll                    ! Mean Richardson number of internal CL 
    real(r8)              :: zbot                     ! Height of CL base
    real(r8)              :: l2rat                    ! Square of ratio of actual to initial CL (not used)

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- ! 

    i = long

    kb_is_surface_mode = 0
    use_dw_surf_mode = 0
    choice_tkes_ebprod_mode = 0
    if( use_dw_surf ) use_dw_surf_mode = 1
    if( choice_tkes .eq. 'ebprod' ) choice_tkes_ebprod_mode = 1

    ! Initialize main output variables
    
    do k = 1, ncvmax
       ricl(i,k) = 0._r8
       ghcl(i,k) = 0._r8
       shcl(i,k) = 0._r8
       smcl(i,k) = 0._r8
       lbrk(i,k) = 0._r8
       wbrk(i,k) = 0._r8
       ebrk(i,k) = 0._r8
    end do
    extend    = .false.
    extend_up = .false.
    extend_dn = .false.

    ! ----------------------------------------------------------- !
    ! Loop over each CL to see if any of them need to be extended !
    ! ----------------------------------------------------------- !

    ncv = 1

    do while( ncv .le. ncvfin(i) )

       ncvinit = ncv
       cntu    = 0
       cntd    = 0
       kb      = kbase(i,ncv) 
       kt      = ktop(i,ncv)
       
       ! ---------------------------------------------------------------------------- !
       ! Calculation of CL interior energetics including surface before extension     !
       ! ---------------------------------------------------------------------------- !
       ! Note that the contribution of interior interfaces (not surface) to 'wint' is !
       ! accounted by using '-sh*l2n2 + sm*l2s2' while the contribution of surface is !
       ! accounted by using 'dwsurf = tkes/b1' when bflxs > 0. This approach is fully !
       ! reasonable. Another possible alternative,  which seems to be also consistent !
       ! is to calculate 'dl2n2_surf'  and  'dl2s2_surf' of surface interfacial layer !
       ! separately, and this contribution is explicitly added by initializing 'l2n2' !
       ! 'l2s2' not by zero, but by 'dl2n2_surf' and 'ds2n2_surf' below.  At the same !
       ! time, 'dwsurf' should be excluded in 'wint' calculation below. The only diff.!
       ! between two approaches is that in case of the latter approach, contributions !
       ! of surface interfacial layer to the CL mean stability function (ri,gh,sh,sm) !
       ! are explicitly included while the first approach is not. In this sense,  the !
       ! second approach seems to be more conceptually consistent,   but currently, I !
       ! (Sungsu) will keep the first default approach. There is a switch             !
       ! 'use_dw_surf' at the first part of eddy_diff.F90 chosing one of              !
       ! these two options.                                                           !
       ! ---------------------------------------------------------------------------- !
       
       ! ------------------------------------------------------ !   
       ! Step 0: Calculate surface interfacial layer energetics !
       ! ------------------------------------------------------ !

       call eddy_diff_zisocl_initial_state(i, kt, kb, pcols, pver, use_dw_surf_mode, choice_tkes_ebprod_mode, z(i,pver), bflxs(i), &
            bprod(i,pver+1), sprod(i,pver+1), tkes(i), z, zi, n2, s2, leng_max, lbulk, gh, sh, sm, dlint_surf, dl2n2_surf, &
            dl2s2_surf, dw_surf, lint, l2n2, l2s2, wint, ricll)
       
       ! Note that at this stage, ( gh, sh, sm )  are the values of surface
       ! interfacial layer if there is no pure internal interface, while if
       ! there is pure internal interface, ( gh, sh, sm ) are the values of
       ! pure CL interfaces or the values that include both the CL internal
       ! interfaces and surface interfaces, depending on the 'use_dw_surf'.       
       
       ! ----------------------------------------------------------------------- !
       ! Perform vertical extension-merging process                              !
       ! ----------------------------------------------------------------------- !
       ! During the merging process, we assumed ( lbulk, sh, sm ) of CL external !
       ! interfaces are the same as the ones of the original merging CL. This is !
       ! an inevitable approximation since we don't know  ( sh, sm ) of external !
       ! interfaces at this stage.     Note that current default merging test is !
       ! purely based on buoyancy production without including shear production, !
       ! since we used 'l2n2' instead of 'wint' as a merging parameter. However, !
       ! merging test based on 'wint' maybe conceptually more attractable.       !
       ! Downward CL merging process is identical to the upward merging process, !
       ! but when the base of extended CL reaches to the surface, surface inter  !
       ! facial layer contribution to the energetic of extended CL must be done  !
       ! carefully depending on the sign of surface buoyancy flux. The contribu  !
       ! tion of surface interfacial layer energetic is included to the internal !
       ! energetics of merging CL only when bflxs > 0.                           !
       ! ----------------------------------------------------------------------- !
       
       ! ---------------------------- !
       ! Step 1. Extend the CL upward !
       ! ---------------------------- !
       
       extend = .false.    ! This will become .true. if CL top or base is extended
       kt_before = kt
       call eddy_diff_zisocl_upward_state(i, kt, pcols, pver, ntop_turb, ncv, lbulk, sh, sm, z, zi, n2, s2, leng_max, ncvfin, &
            kbase, ktop, lint, l2n2, l2s2, wint, upward_state_status)
       if( kt .ne. kt_before ) then
           extend    = .true.
           extend_up = .true.
       end if
       if( upward_state_status .ne. 0 ) then
           write(iulog,*) 'zisocl: Error: Tried to extend CL to the model top'
           call endrun('zisocl: Error: Tried to extend CL to the model top')
       end if

       ! ------------------------------ !
       ! Step 2. Extend the CL downward !
       ! ------------------------------ !
       
       kb_before = kb
       call eddy_diff_zisocl_downward_state(i, kb, ncv, ncvinit, pcols, pver, lbulk, sh, sm, z(i,pver), bflxs(i), &
            bprod(i,pver+1), sprod(i,pver+1), tkes(i), z, zi, n2, s2, leng_max, dlint_surf, dl2n2_surf, dl2s2_surf, dw_surf, &
            ncvfin, kbase, ktop, lint, l2n2, l2s2, wint, downward_state_status)
       if( kb .ne. kb_before ) then
           extend    = .true.
           extend_dn = .true.
       end if
       if( downward_state_status .ne. 0 ) then
           write(iulog,*) 'Major mistake zisocl: the CL based at surface is not indexed 1'
           call endrun('Major mistake zisocl: the CL based at surface is not indexed 1')
       end if

       ! Sanity check for positive wint.

       if( wint .lt. 0.01_r8 ) then
           wint = 0.01_r8
       end if

       ! -------------------------------------------------------------------------- !
       ! Finally update CL mean internal energetics including surface contribution  !
       ! after finishing all the CL extension-merging process.  As mentioned above, !
       ! there are two possible ways in the treatment of surface interfacial layer, !
       ! either through 'dw_surf' or 'dl2n2_surf and dl2s2_surf' by setting logical !
       ! variable 'use_dw_surf' =.true. or .false.    In any cases, we should avoid !
       ! double counting of surface interfacial layer and one single consistent way !
       ! should be used throughout the program.                                     !
       ! -------------------------------------------------------------------------- !

       if( extend ) then

           ktop(i,ncv)  = kt
           kbase(i,ncv) = kb
           call eddy_diff_zisocl_extended_state(i, kt, kb, pcols, pver, use_dw_surf_mode, zi(i,kt), zi(i,kb), z(i,pver), bflxs(i), &
                bprod(i,pver+1), sprod(i,pver+1), tkes(i), z, zi, n2, s2, leng_max, gh, sh, sm, lint, wint, ricll)

       end if

       ! ---------------------------------------------------------------------- !
       ! Calculate final output variables of each CL (either has merged or not) !
       ! ---------------------------------------------------------------------- !

       lbrk(i,ncv) = lint
       wbrk(i,ncv) = wint/lint
       ebrk(i,ncv) = b1*wbrk(i,ncv)
       ebrk(i,ncv) = min(ebrk(i,ncv),tkemax)
       ricl(i,ncv) = ricll 
       ghcl(i,ncv) = gh 
       shcl(i,ncv) = sh
       smcl(i,ncv) = sm

       ! Increment counter for next CL. I should check if the increament of 'ncv'
       ! below is reasonable or not, since whenever CL is merged during downward
       ! extension process, 'ncv' is lowered down continuously within 'do' loop.
       ! But it seems that below 'ncv = ncv + 1' is perfectly correct.

       ncv = ncv + 1

    end do                   ! End of loop over each CL regime, ncv.

    ! ---------------------------------------------------------- !
    ! Re-initialize external interface indices which are not CLs !
    ! ---------------------------------------------------------- !

    do ncv = ncvfin(i) + 1, ncvmax
       ktop(i,ncv)  = 0
       kbase(i,ncv) = 0
    end do

    ! ------------------------------------------------ !
    ! Update CL interface identifiers, 'belongcv'      !
    ! CL external interfaces are also identified as CL !
    ! ------------------------------------------------ !

    do k = 1, pver + 1
       belongcv(i,k) = .false.
    end do

    do ncv = 1, ncvfin(i)
       do k = ktop(i,ncv), kbase(i,ncv)
          belongcv(i,k) = .true.
       end do
    end do

    return

    end subroutine zisocl

    real(r8) function compute_cubic(a,b,c)
    ! ------------------------------------------------------------------------- !
    ! Solve canonical cubic : x^3 + a*x^2 + b*x + c = 0,  x = sqrt(e)/sqrt(<e>) !
    ! Set x = max(xmin,x) at the end                                            ! 
    ! ------------------------------------------------------------------------- !
    implicit none
    real(r8), intent(in)     :: a, b, c
    real(r8)  qq, rr, dd, theta, aa, bb, x1, x2, x3
    real(r8), parameter      :: xmin = 1.e-2_r8
    
    qq = (a**2-3._r8*b)/9._r8 
    rr = (2._r8*a**3 - 9._r8*a*b + 27._r8*c)/54._r8
    
    dd = rr**2 - qq**3
    if( dd .le. 0._r8 ) then
        theta = acos(rr/qq**(3._r8/2._r8))
        x1 = -2._r8*sqrt(qq)*cos(theta/3._r8) - a/3._r8
        x2 = -2._r8*sqrt(qq)*cos((theta+2._r8*3.141592_r8)/3._r8) - a/3._r8
        x3 = -2._r8*sqrt(qq)*cos((theta-2._r8*3.141592_r8)/3._r8) - a/3._r8
        compute_cubic = max(max(max(x1,x2),x3),xmin)        
        return
    else
        if( rr .ge. 0._r8 ) then
            aa = -(sqrt(rr**2-qq**3)+rr)**(1._r8/3._r8)
        else
            aa =  (sqrt(rr**2-qq**3)-rr)**(1._r8/3._r8)
        endif
        if( aa .eq. 0._r8 ) then
            bb = 0._r8
        else
            bb = qq/aa
        endif
        compute_cubic = max((aa+bb)-a/3._r8,xmin) 
        return
    endif

    return
    end function compute_cubic

END MODULE eddy_diff
