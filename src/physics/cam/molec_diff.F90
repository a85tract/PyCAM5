module molec_diff

  !------------------------------------------------------------------------------------------------- !
  ! Module to compute molecular diffusivity for various constituents                                 !
  !                                                                                                  !
  ! Public interfaces :                                                                              !
  !                                                                                                  !
  !    init_molec_diff           Initializes time independent coefficients                           !
  !    init_timestep_molec_diff  Time-step initialization for molecular diffusivity                  !
  !    compute_molec_diff        Computes constituent-independent terms for moleculuar diffusivity   !
  !    vd_lu_qdecomp             Computes constituent-dependent terms for moleculuar diffusivity and !
  !                              updates terms in the triadiagonal matrix used for the implicit      !
  !                              solution of the diffusion equation                                  !
  !                                                                                                  !
  !---------------------------Code history---------------------------------------------------------- !
  ! Modularized     :  J. McCaa, September 2004                                                      !
  ! Lastly Arranged :  S. Park,  January.  2010                                                      !
  !                    M. Mills, November  2011
  !------------------------------------------------------------------------------------------------- !

  use perf_mod
  use physconst,    only : mbarv
  use constituents, only : pcnst
  use phys_control, only : waccmx_is             !WACCM-X runtime switch
  use ref_pres,     only : nbot_molec, ntop_molec

  implicit none
  private
  save

  public init_molec_diff
  public init_timestep_molec_diff
  public compute_molec_diff
  public vd_lu_qdecomp
  public molec_diff_get_inputs

  ! ---------- !
  ! Parameters !
  ! ---------- !

  integer,  parameter   :: r8 = selected_real_kind(12) ! 8 byte real

  real(r8), parameter   :: km_fac = 3.55E-7_r8         ! Molecular viscosity constant [ unit ? ]
  real(r8), parameter   :: pr_num = 1._r8              ! Prandtl number [ no unit ]
  real(r8), parameter   :: pwr    = 2._r8/3._r8        ! Exponentiation factor [ no unit ]
  real(r8), parameter   :: d0     = 1.52E20_r8         ! Diffusion factor [ m-1 s-1 ] molec sqrt(kg/kmol/K) [ unit ? ]
                                                       ! Aerononmy, Part B, Banks and Kockarts (1973), p39
                                                       ! Note text cites 1.52E18 cm-1 ...

  real(r8)              :: rair                        ! Gas constant for dry air
  real(r8)              :: mw_dry                      ! Molecular weight of dry air
  real(r8)              :: n_avog                      ! Avogadro's number [ molec/kmol ]
  real(r8)              :: gravit
  real(r8)              :: cpair
  real(r8)              :: kbtz                        ! Boltzman constant

  real(r8), allocatable :: mw_fac(:)                   ! sqrt(1/M_q + 1/M_d) in constituent diffusivity [  unit ? ]
  real(r8), allocatable :: alphath(:)                  ! Thermal diffusion factor, -0.38 for H, 0 for others

  logical :: waccmx_mode = .false.

contains

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine init_molec_diff( kind, ncnst, rair_in, mw_dry_in, n_avog_in, gravit_in, &
                              cpair_in, kbtz_in, errstring)

    use constituents,     only : cnst_mw, cnst_get_ind
    use upper_bc,         only : ubc_init
    use physics_buffer,   only : physics_buffer_desc

    integer,  intent(in)  :: kind           ! Kind of reals being passed in
    integer,  intent(in)  :: ncnst          ! Number of constituents
    real(r8), intent(in)  :: rair_in
    real(r8), intent(in)  :: mw_dry_in      ! Molecular weight of dry air
    real(r8), intent(in)  :: n_avog_in      ! Avogadro's number [ molec/kmol ]
    real(r8), intent(in)  :: gravit_in
    real(r8), intent(in)  :: cpair_in
    real(r8), intent(in)  :: kbtz_in        ! Boltzman constant

    character(len=*), intent(out) :: errstring

    ! Local

    integer               :: k              ! Level index
    integer               :: m              ! Constituent index
    integer               :: indx_H         ! Constituent index for H
    integer               :: ierr           ! Allocate error check

    errstring = ' '

    rair       = rair_in
    mw_dry     = mw_dry_in
    n_avog     = n_avog_in
    gravit     = gravit_in
    cpair      = cpair_in
    kbtz       = kbtz_in

    if( kind /= r8 ) then
       errstring = 'inconsistent KIND of reals passed to init_molec_diff'
       return
    end if

    ! Determine whether WACCM-X is on.
    waccmx_mode = waccmx_is('ionosphere') .or. waccmx_is('neutral')

  ! Initialize upper boundary condition variables

    call ubc_init()

  ! Molecular weight factor in constitutent diffusivity
  ! ***** FAKE THIS FOR NOW USING MOLECULAR WEIGHT OF DRY AIR FOR ALL TRACERS ****

    allocate(mw_fac(ncnst))
    do m = 1, ncnst
       mw_fac(m) = d0 * mw_dry * sqrt(1._r8/mw_dry + 1._r8/cnst_mw(m)) / n_avog
    end do

    !--------------------------------------------------------------------------------------------
    ! For WACCM-X, get H data index and initialize thermal diffusion coefficient
    !--------------------------------------------------------------------------------------------
    if ( waccmx_mode ) then

      call cnst_get_ind('H',  indx_H)

      allocate(alphath(ncnst), stat=ierr)
      if ( ierr /= 0 ) then
         errstring = 'allocate failed in init_molec_diff'
         return
      end if
      alphath(:ncnst) = 0._r8
      alphath(indx_H) = -0.38_r8

    endif

  end subroutine init_molec_diff

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine init_timestep_molec_diff(pbuf2d, state)
    !--------------------------- !
    ! Timestep dependent setting !
    !--------------------------- !
    use upper_bc,     only : ubc_timestep_init
    use physics_types,only: physics_state
    use ppgrid,       only: begchunk, endchunk
    use physics_buffer, only : physics_buffer_desc

    type(physics_state), intent(in) :: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    call ubc_timestep_init( pbuf2d, state)

  end subroutine init_timestep_molec_diff

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine molec_diff_get_inputs( &
       lchnk, pcols_in, pver_in, ncnst, ncol, pint, zi, &
       waccmx_mode_out, ntop_molec_out, nbot_molec_out, &
       d0_out, km_fac_out, pr_num_out, pwr_out, n_avog_out, mw_dry_out, &
       mbarv_out, rairv_out, kmvis_out, kmcnd_out, cnst_mw_out, &
       cnst_fixed_ubc_out, cnst_fixed_ubflx_out, mw_fac_out, &
       alphath_out, ubc_t_out, ubc_mmr_out, ubc_flux_out)

    use upper_bc, only: ubc_get_vals
    use constituents, only: cnst_mw, cnst_fixed_ubc, cnst_fixed_ubflx
    use physconst, only: rairv, kmvis, kmcnd

    integer, intent(in) :: lchnk, pcols_in, pver_in, ncnst, ncol
    real(r8), intent(in) :: pint(pcols_in,pver_in+1)
    real(r8), intent(in) :: zi(pcols_in,pver_in+1)
    logical, intent(out) :: waccmx_mode_out
    integer, intent(out) :: ntop_molec_out, nbot_molec_out
    real(r8), intent(out) :: d0_out, km_fac_out, pr_num_out, pwr_out
    real(r8), intent(out) :: n_avog_out, mw_dry_out
    real(r8), intent(out) :: mbarv_out(pcols_in,pver_in)
    real(r8), intent(out) :: rairv_out(pcols_in,pver_in)
    real(r8), intent(out) :: kmvis_out(pcols_in,pver_in+1)
    real(r8), intent(out) :: kmcnd_out(pcols_in,pver_in+1)
    real(r8), intent(out) :: cnst_mw_out(ncnst)
    logical, intent(out) :: cnst_fixed_ubc_out(ncnst)
    logical, intent(out) :: cnst_fixed_ubflx_out(ncnst)
    real(r8), intent(out) :: mw_fac_out(ncnst)
    real(r8), intent(out) :: alphath_out(ncnst)
    real(r8), intent(out) :: ubc_t_out(pcols_in)
    real(r8), intent(out) :: ubc_mmr_out(pcols_in,ncnst)
    real(r8), intent(out) :: ubc_flux_out(ncnst)

    call ubc_get_vals(lchnk, ncol, ntop_molec, pint, zi, ubc_t_out, &
         ubc_mmr_out, ubc_flux_out)

    waccmx_mode_out = waccmx_mode
    ntop_molec_out = ntop_molec
    nbot_molec_out = nbot_molec
    d0_out = d0
    km_fac_out = km_fac
    pr_num_out = pr_num
    pwr_out = pwr
    n_avog_out = n_avog
    mw_dry_out = mw_dry

    mbarv_out = mbarv(:,:,lchnk)
    rairv_out = rairv(:,:,lchnk)
    cnst_mw_out = cnst_mw(:ncnst)
    cnst_fixed_ubc_out = cnst_fixed_ubc(:ncnst)
    cnst_fixed_ubflx_out = cnst_fixed_ubflx(:ncnst)
    mw_fac_out = mw_fac(:ncnst)

    alphath_out = 0._r8
    kmvis_out = 0._r8
    kmcnd_out = 0._r8
    if (waccmx_mode) then
       alphath_out = alphath(:ncnst)
       kmvis_out = kmvis(:,:,lchnk)
       kmcnd_out = kmcnd(:,:,lchnk)
    end if

  end subroutine molec_diff_get_inputs

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  integer function compute_molec_diff( lchnk             ,                                          &
       pcols             , pver                , ncnst     , ncol     , t      , pmid   , pint   ,  &
       zi                , ztodt               , kvm       , kvt      , tint   , rhoi   , tmpi2  ,  &
       kq_scal           , ubc_t               , ubc_mmr   , ubc_flux , dse_top, cc_top ,           &
       cnst_mw_out       , cnst_fixed_ubc_out  , cnst_fixed_ubflx_out , mw_fac_out      ,           &
       ntop_molec_out    , nbot_molec_out      , kvt_returned )

    use molec_diff_kernel, only: compute_molec_diff_kernel
    use physconst, only: cpairv

    ! --------------------- !
    ! Input-Output Argument !
    ! --------------------- !

    integer,  intent(in)    :: pcols
    integer,  intent(in)    :: pver
    integer,  intent(in)    :: ncnst
    integer,  intent(in)    :: ncol                      ! Number of atmospheric columns
    integer,  intent(in)    :: lchnk                     ! Chunk identifier
    real(r8), intent(in)    :: t(pcols,pver)             ! Temperature input
    real(r8), intent(in)    :: pmid(pcols,pver)          ! Midpoint pressures
    real(r8), intent(in)    :: pint(pcols,pver+1)        ! Interface pressures
    real(r8), intent(in)    :: zi(pcols,pver+1)          ! Interface heights
    real(r8), intent(in)    :: ztodt                     ! 2 delta-t

    real(r8), intent(inout) :: kvm(pcols,pver+1)         ! Viscosity ( diffusivity for momentum )
    real(r8), intent(out)   :: kvt(pcols,pver+1)         ! Kinematic molecular conductivity
    real(r8), intent(inout) :: tint(pcols,pver+1)        ! Interface temperature
    real(r8), intent(inout) :: rhoi(pcols,pver+1)        ! Density ( rho ) at interfaces
    real(r8), intent(inout) :: tmpi2(pcols,pver+1)       ! dt*(g*rho)**2/dp at interfaces

    real(r8), intent(out)   :: kq_scal(pcols,pver+1)     ! kq_fac*sqrt(T)*m_d/rho for molecular diffusivity
    real(r8), intent(out)   :: ubc_t(pcols)              ! Upper boundary temperature (K)
    real(r8), intent(out)   :: ubc_mmr(pcols,ncnst)      ! Upper boundary mixing ratios [ kg/kg ]
    real(r8), intent(out)   :: ubc_flux(ncnst)           ! Upper boundary flux [ kg/s/m^2 ]
    real(r8), intent(out)   :: cnst_mw_out(ncnst)
    logical,  intent(out)   :: cnst_fixed_ubc_out(ncnst)
    logical,  intent(out)   :: cnst_fixed_ubflx_out(ncnst)
    real(r8), intent(out)   :: mw_fac_out(pcols,pver+1,ncnst) ! composition dependent mw_fac on interface level
    real(r8), intent(out)   :: dse_top(pcols)            ! dse on top boundary
    real(r8), intent(out)   :: cc_top(pcols)             ! Lower diagonal at top interface
    integer,  intent(out)   :: ntop_molec_out
    integer,  intent(out)   :: nbot_molec_out
    logical,  intent(out)   :: kvt_returned              ! Whether we actually returned kvt (vs. kvh).

    ! --------------- !
    ! Local variables !
    ! --------------- !

    logical :: waccmx_mode_in
    integer :: ntop_molec_in, nbot_molec_in
    real(r8) :: d0_in, km_fac_in, pr_num_in, pwr_in
    real(r8) :: n_avog_in, mw_dry_in
    real(r8) :: mbarv_in(pcols,pver)
    real(r8) :: rairv_in(pcols,pver)
    real(r8) :: kmvis_in(pcols,pver+1)
    real(r8) :: kmcnd_in(pcols,pver+1)
    real(r8) :: cnst_mw_in(ncnst)
    logical :: cnst_fixed_ubc_in(ncnst)
    logical :: cnst_fixed_ubflx_in(ncnst)
    real(r8) :: mw_fac_in(ncnst)
    real(r8) :: alphath_in(ncnst)
    real(r8) :: ubc_t_in(pcols)
    real(r8) :: ubc_mmr_in(pcols,ncnst)
    real(r8) :: ubc_flux_in(ncnst)

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !

    call molec_diff_get_inputs( &
         lchnk, pcols, pver, ncnst, ncol, pint, zi, &
         waccmx_mode_in, ntop_molec_in, nbot_molec_in, &
         d0_in, km_fac_in, pr_num_in, pwr_in, n_avog_in, mw_dry_in, &
         mbarv_in, rairv_in, kmvis_in, kmcnd_in, cnst_mw_in, &
         cnst_fixed_ubc_in, cnst_fixed_ubflx_in, mw_fac_in, &
         alphath_in, ubc_t_in, ubc_mmr_in, ubc_flux_in)

    ntop_molec_out = ntop_molec_in
    nbot_molec_out = nbot_molec_in

    call compute_molec_diff_kernel( &
         pcols, pver, ncnst, ncol, t, pmid, pint, zi, ztodt, &
         waccmx_mode_in, ntop_molec_in, nbot_molec_in, gravit, &
         d0_in, km_fac_in, pr_num_in, pwr_in, n_avog_in, mw_dry_in, &
         mbarv_in, rairv_in, kmvis_in, kmcnd_in, cpairv(:,:,lchnk), &
         cnst_mw_in, cnst_fixed_ubc_in, cnst_fixed_ubflx_in, &
         mw_fac_in, ubc_t_in, ubc_mmr_in, ubc_flux_in, kvm, kvt, &
         tint, rhoi, tmpi2, kq_scal, ubc_t, ubc_mmr, ubc_flux, &
         dse_top, cc_top, cnst_mw_out, cnst_fixed_ubc_out, &
         cnst_fixed_ubflx_out, mw_fac_out, kvt_returned)

    compute_molec_diff = 1
    return
  end function compute_molec_diff

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  function vd_lu_qdecomp( &
       pcols , pver   , ncol       , fixed_ubc  , mw     , &
       kv    , kq_scal, mw_facm    , dpidz_sq   , p      , &
       interface_boundary, molec_boundary, rhoi   ,        &
       tint  , ztodt  , ntop_molec , nbot_molec , nbot   , &
       lchnk , t          , m      , no_molec_decomp)      result(decomp)

    use coords_1d, only: Coords1D
    use linear_1d_operators, only: BoundaryType, TriDiagDecomp
    use molec_diff_kernel, only: vd_lu_qdecomp_kernel

    !------------------------------------------------------------------------------ !
    ! Add the molecular diffusivity to the turbulent diffusivity for a consitutent. !
    ! Update the superdiagonal (ca(k)), diagonal (cb(k)) and subdiagonal (cc(k))    !
    ! coefficients of the tridiagonal diffusion matrix, also ze and denominator.    !
    !------------------------------------------------------------------------------ !

    ! ---------------------- !
    ! Input-Output Arguments !
    ! ---------------------- !

    integer,  intent(in)    :: pcols
    integer,  intent(in)    :: pver
    integer,  intent(in)    :: ncol                  ! Number of atmospheric columns

    integer,  intent(in)    :: ntop_molec
    integer,  intent(in)    :: nbot_molec
    integer,  intent(in)    :: nbot

    logical,  intent(in)    :: fixed_ubc             ! Fixed upper boundary condition flag
    real(r8), intent(in)    :: kv(pcols,pver+1)      ! Eddy diffusivity
    real(r8), intent(in)    :: kq_scal(pcols,pver+1) ! Molecular diffusivity ( kq_fac*sqrt(T)*m_d/rho )
    real(r8), intent(in)    :: mw                    ! Molecular weight for this constituent
    real(r8), intent(in)    :: mw_facm(pcols,pver+1) ! composition dependent sqrt(1/M_q + 1/M_d) for this constituent
    real(r8), intent(in)    :: dpidz_sq(ncol,pver+1) ! (g*rho)**2 (square of vertical derivative of pint)
    type(Coords1D), intent(in) :: p                  ! Pressure coordinates
    type(BoundaryType), intent(in) :: interface_boundary ! Boundary on grid edge.
    type(BoundaryType), intent(in) :: molec_boundary ! Boundary at edge of molec_diff region.
    real(r8), intent(in)    :: rhoi(pcols,pver+1)    ! Density at interfaces [ kg/m3 ]
    real(r8), intent(in)    :: tint(pcols,pver+1)    ! Interface temperature [ K ]
    real(r8), intent(in)    :: ztodt                 ! 2 delta-t [ s ]

    integer,  intent(in)    :: lchnk		    ! Chunk number
    real(r8), intent(in)    :: t(pcols,pver)	    ! temperature
    integer,  intent(in)    :: m 		    ! cnst index

    ! Decomposition covering levels without vertical diffusion.
    type(TriDiagDecomp), intent(in) :: no_molec_decomp

    ! LU decomposition information for solver.
    type(TriDiagDecomp) :: decomp

    ! --------------- !
    ! Local Variables !
    ! --------------- !

    real(r8) :: alphath_m

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !

    ! --------------------------------------------------------------------- !
    ! Determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs of the !
    ! tridiagonal diffusion matrix. The diagonal elements  (cb=1+ca+cc) are !
    ! a combination of ca and cc; they are not required by the solver.      !
    !---------------------------------------------------------------------- !

    call t_startf('vd_lu_qdecomp')

    alphath_m = 0._r8
    if (waccmx_mode) alphath_m = alphath(m)

    decomp = vd_lu_qdecomp_kernel( &
         pcols, pver, ncol, fixed_ubc, mw, kv, kq_scal, mw_facm, &
         dpidz_sq, p, interface_boundary, molec_boundary, rhoi, tint, &
         ztodt, ntop_molec, nbot_molec, nbot, t, waccmx_mode, &
         mw_dry, mbarv(:,:,lchnk), alphath_m, no_molec_decomp)

    call t_stopf('vd_lu_qdecomp')

  end function vd_lu_qdecomp

end module molec_diff
