module gw_common

!
! This module contains code common to different gravity wave
! parameterizations.
!
use gw_utils, only: r8
use coords_1d, only: Coords1D
use spmd_utils, only: masterproc
use cam_logfile, only: iulog
use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr


implicit none
private
save

! Public interface.

public :: GWBand

public :: gw_common_init
public :: gw_prof
public :: gw_drag_prof
public :: calc_taucd, momentum_flux, momentum_fixer
public :: energy_change, energy_fixer
public :: coriolis_speed, adjust_inertial

public :: pver
public :: west, east, north, south
public :: pi
public :: gravit
public :: rair

! Number of levels in the atmosphere.
integer, target, protected :: pver = 0

! Whether or not to enforce an upper boundary condition of tau = 0.
logical, target :: tau_0_ubc = .false.

! Index the cardinal directions.
integer, parameter :: west = 1
integer, parameter :: east = 2
integer, parameter :: south = 3
integer, parameter :: north = 4

! 3.14159...
real(r8), parameter :: pi = acos(-1._r8)

! Acceleration due to gravity.
real(r8), target, protected :: gravit = huge(1._r8)

! Gas constant for dry air.
real(r8), target, protected :: rair = huge(1._r8)

!
! Private variables
!

! Interface levels for gravity wave sources.
integer, target :: ktop = huge(1)

! Background diffusivity.
real(r8), parameter :: dback = 0.05_r8

! rair/gravit
real(r8), target :: rog = huge(1._r8)

! Newtonian cooling coefficients.
real(r8), allocatable, target :: alpha(:)

!
! Limits to keep values reasonable.
!

! Minimum non-zero stress.
real(r8), parameter :: taumin = 1.e-10_r8
! Maximum wind tendency from stress divergence (before efficiency applied).
! 400 m/s/day
real(r8), parameter :: tndmax = 400._r8 / 86400._r8
! Maximum allowed change in u-c (before efficiency applied).
real(r8), parameter :: umcfac = 0.5_r8
! Minimum value of (u-c)**2.
real(r8), parameter :: ubmc2mn = 0.01_r8

logical :: use_native_gw_prof_impl = .false.
logical :: gw_prof_impl_selected = .false.
logical :: gw_prof_entered_logged = .false.
logical :: use_native_energy_change_impl = .false.
logical :: energy_change_impl_selected = .false.
logical :: energy_change_entered_logged = .false.
logical :: use_native_gw_drag_prof_core_impl = .false.
logical :: gw_drag_prof_core_impl_selected = .false.
logical :: gw_drag_prof_core_entered_logged = .false.
logical :: use_native_gw_diff_solver_impl = .false.
logical :: gw_diff_solver_impl_selected = .false.
logical :: gw_diff_solver_entered_logged = .false.
logical :: use_native_gw_common_init_impl = .false.
logical :: gw_common_init_impl_selected = .false.
logical :: gw_common_init_direct_logged = .false.
logical :: use_native_new_gwband_impl = .false.
logical :: new_gwband_impl_selected = .false.
logical :: new_gwband_direct_logged = .false.

interface
   subroutine gw_common_new_gwband_codon(ngwv_c, dc_c, fcrit2_c, wavelength_c, pi_c, &
        ngwv_p, dc_p, fcrit2_p, cref_p, kwv_p, effkwv_p) bind(c, name="gw_common_new_gwband_codon")
     use iso_c_binding, only: c_double, c_int64_t, c_ptr
     integer(c_int64_t), value :: ngwv_c
     real(c_double), value :: dc_c, fcrit2_c, wavelength_c, pi_c
     type(c_ptr), value :: ngwv_p, dc_p, fcrit2_p, cref_p, kwv_p, effkwv_p
   end subroutine gw_common_new_gwband_codon
   subroutine gw_common_init_scalars_codon(pver_in_c, ktop_in_c, gravit_in_c, rair_in_c, &
        tau_0_ubc_in_c, pver_p, tau_0_ubc_p, ktop_p, gravit_p, rair_p, rog_p) &
        bind(c, name="gw_common_init_scalars_codon")
     use iso_c_binding, only: c_double, c_int64_t, c_ptr
     integer(c_int64_t), value :: pver_in_c, ktop_in_c, tau_0_ubc_in_c
     real(c_double), value :: gravit_in_c, rair_in_c
     type(c_ptr), value :: pver_p, tau_0_ubc_p, ktop_p, gravit_p, rair_p, rog_p
   end subroutine gw_common_init_scalars_codon
   subroutine gw_common_init_alpha_codon(n_c, alpha_in_p, alpha_p) bind(c, name="gw_common_init_alpha_codon")
     use iso_c_binding, only: c_int64_t, c_ptr
     integer(c_int64_t), value :: n_c
     type(c_ptr), value :: alpha_in_p, alpha_p
   end subroutine gw_common_init_alpha_codon
end interface

! Type describing a band of wavelengths into which gravity waves can be
! emitted.
! Currently this has to have uniform spacing (i.e. adjacent elements of
! cref are exactly dc apart).
type :: GWBand
   ! Dimension of the spectrum.
   integer :: ngwv
   ! Delta between nearest phase speeds [m/s].
   real(r8) :: dc
   ! Reference speeds [m/s].
   real(r8), allocatable :: cref(:)
   ! Critical Froude number, squared (usually 1, but CAM3 used 0.5).
   real(r8) :: fcrit2
   ! Horizontal wave number [1/m].
   real(r8) :: kwv
   ! Effective horizontal wave number [1/m] (fcrit2*kwv).
   real(r8) :: effkwv
end type GWBand

interface GWBand
   module procedure new_GWBand
end interface

contains

!==========================================================================

! Constructor for a GWBand that calculates derived components.
function new_GWBand(ngwv, dc, fcrit2, wavelength) result(band)
  ! Used directly to set the type's components.
  integer, intent(in) :: ngwv
  real(r8), intent(in) :: dc
  real(r8), intent(in) :: fcrit2
  ! Wavelength in meters.
  real(r8), intent(in) :: wavelength

  ! Output.
  type(GWBand), target :: band

  ! Wavenumber index.
  integer :: l

  call new_gwband_select_impl()

  ! Uniform phase speed reference grid.
  allocate(band%cref(-ngwv:ngwv))

  if (use_native_new_gwband_impl) then
     ! Simple assignments.
     band%ngwv = ngwv
     band%dc = dc
     band%fcrit2 = fcrit2

     ! Uniform phase speed reference grid.
     band%cref = [( dc * l, l = -ngwv, ngwv )]

     ! Wavenumber and effective wavenumber come from the wavelength.
     band%kwv = 2._r8*pi / wavelength
     band%effkwv = band%fcrit2 * band%kwv
  else
     call gw_common_new_gwband_codon(int(ngwv, c_int64_t), real(dc, c_double), &
          real(fcrit2, c_double), real(wavelength, c_double), real(pi, c_double), &
          c_loc(band%ngwv), c_loc(band%dc), c_loc(band%fcrit2), c_loc(band%cref(-ngwv)), &
          c_loc(band%kwv), c_loc(band%effkwv))
     call new_gwband_note_direct()
  end if

end function new_GWBand

!==========================================================================

subroutine new_gwband_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (new_gwband_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('GW_COMMON_NEW_GWBAND_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_new_gwband_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_new_gwband_impl = .false.
  end if

  new_gwband_impl_selected = .true.

  if (masterproc) then
     if (use_native_new_gwband_impl) then
        write(iulog,*) 'new_GWBand implementation = native'
     else
        write(iulog,*) 'new_GWBand implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine new_gwband_select_impl

!==========================================================================

subroutine new_gwband_note_direct()

  if (new_gwband_direct_logged) return
  new_gwband_direct_logged = .true.

  if (masterproc) then
     write(iulog,*) 'new_GWBand direct = codon'
     call flush(iulog)
  end if

end subroutine new_gwband_note_direct

!==========================================================================

subroutine gw_common_init(pver_in, &
     tau_0_ubc_in, ktop_in, gravit_in, rair_in, alpha_in, errstring)

  integer,  intent(in) :: pver_in
  logical,  intent(in) :: tau_0_ubc_in
  integer,  intent(in) :: ktop_in
  real(r8), intent(in) :: gravit_in
  real(r8), intent(in) :: rair_in
  real(r8), target, intent(in) :: alpha_in(:)
  ! Report any errors from this routine.
  character(len=*), intent(out) :: errstring

  integer :: ierr

  call gw_common_init_select_impl()

  errstring = ""

  allocate(alpha(pver_in+1), stat=ierr, errmsg=errstring)
  if (ierr /= 0) return

  if (use_native_gw_common_init_impl) then
     pver = pver_in
     tau_0_ubc = tau_0_ubc_in
     ktop = ktop_in
     gravit = gravit_in
     rair = rair_in
     alpha = alpha_in
     rog = rair/gravit
  else
     call gw_common_init_scalars_codon(int(pver_in, c_int64_t), int(ktop_in, c_int64_t), &
          real(gravit_in, c_double), real(rair_in, c_double), &
          merge(1_c_int64_t, 0_c_int64_t, tau_0_ubc_in), &
          c_loc(pver), c_loc(tau_0_ubc), c_loc(ktop), c_loc(gravit), c_loc(rair), c_loc(rog))
     call gw_common_init_alpha_codon(int(pver_in + 1, c_int64_t), c_loc(alpha_in(1)), c_loc(alpha(1)))
     call gw_common_init_note_direct()
  end if

end subroutine gw_common_init

!==========================================================================

subroutine gw_common_init_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (gw_common_init_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('GW_COMMON_INIT_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_gw_common_init_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_gw_common_init_impl = .false.
  end if

  gw_common_init_impl_selected = .true.

  if (masterproc) then
     if (use_native_gw_common_init_impl) then
        write(iulog,*) 'gw_common_init implementation = native'
     else
        write(iulog,*) 'gw_common_init implementation = codon'
     end if
     call flush(iulog)
  end if

end subroutine gw_common_init_select_impl

!==========================================================================

subroutine gw_common_init_note_direct()

  if (gw_common_init_direct_logged) return
  gw_common_init_direct_logged = .true.

  if (masterproc) then
     write(iulog,*) 'gw_common_init direct = codon'
     call flush(iulog)
  end if

end subroutine gw_common_init_note_direct

!==========================================================================

subroutine gw_prof (ncol, p, cpair, t, rhoi, nm, ni)
  !-----------------------------------------------------------------------
  ! Selectable wrapper for background gravity-wave profile calculations.
  !-----------------------------------------------------------------------
  integer, intent(in) :: ncol
  type(Coords1D), intent(in) :: p
  real(r8), intent(in) :: cpair
  real(r8), intent(in) :: t(ncol,pver)
  real(r8), intent(out) :: rhoi(ncol,pver+1)
  real(r8), intent(out) :: nm(ncol,pver), ni(ncol,pver+1)

  real(r8), target :: ti(ncol,pver+1)

  call gw_prof_select_impl()

  if (use_native_gw_prof_impl) then
     call gw_prof_native(ncol, p, cpair, t, rhoi, nm, ni)
  else
     call gw_prof_note_entered()
     call gw_prof_codon_wrap(ncol, pver, cpair, rair, gravit, p%ifc, p%rdst, &
          t, rhoi, nm, ni, ti)
  end if

end subroutine gw_prof

!==========================================================================

subroutine gw_prof_append_proof(proof_line)

  character(len=*), intent(in) :: proof_line

  character(len=512) :: proof_file
  integer :: status, n, unitno

  proof_file = ''
  call get_environment_variable('GW_PROF_PROOF_FILE', value=proof_file, length=n, status=status)
  if (status == 0 .and. n > 0) then
     open(newunit=unitno, file=trim(proof_file(:n)), status='unknown', position='append', action='write')
     write(unitno,'(A)') trim(proof_line)
     close(unitno)
  end if

end subroutine gw_prof_append_proof

!==========================================================================

subroutine gw_prof_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (gw_prof_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('GW_PROF_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_gw_prof_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_gw_prof_impl = .false.
  end if

  gw_prof_impl_selected = .true.

  if (masterproc) then
     if (use_native_gw_prof_impl) then
        write(iulog,*) 'gw_prof implementation = native'
        call gw_prof_append_proof('gw_prof selector entered implementation = native')
     else
        write(iulog,*) 'gw_prof implementation = codon'
        call gw_prof_append_proof('gw_prof selector entered implementation = codon')
     end if
     call flush(iulog)
  end if

end subroutine gw_prof_select_impl

!==========================================================================

subroutine gw_prof_note_entered()

  if (gw_prof_entered_logged) return
  gw_prof_entered_logged = .true.

  if (masterproc) then
     write(iulog,*) 'gw_prof entered (unified background-profile stage dispatch = codon)'
     call gw_prof_append_proof('gw_prof entered (unified background-profile stage dispatch = codon)')
     call flush(iulog)
  end if

end subroutine gw_prof_note_entered

!==========================================================================

subroutine gw_prof_codon_wrap(ncol_local, pver_local, cpair_local, rair_local, gravit_local, &
     p_ifc_local, p_rdst_local, t_local, rhoi_local, nm_local, ni_local, ti_local)

  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

  integer, intent(in) :: ncol_local, pver_local
  real(r8), intent(in) :: cpair_local, rair_local, gravit_local
  real(r8), target, intent(in) :: p_ifc_local(ncol_local,pver_local+1)
  real(r8), target, intent(in) :: p_rdst_local(ncol_local,pver_local-1)
  real(r8), target, intent(in) :: t_local(ncol_local,pver_local)
  real(r8), target, intent(inout) :: rhoi_local(ncol_local,pver_local+1)
  real(r8), target, intent(inout) :: nm_local(ncol_local,pver_local)
  real(r8), target, intent(inout) :: ni_local(ncol_local,pver_local+1)
  real(r8), target, intent(inout) :: ti_local(ncol_local,pver_local+1)

  interface
     subroutine gw_prof_stage_dispatch_codon(ncol_c, pver_c, cpair_c, rair_c, gravit_c, &
          p_ifc_p, p_rdst_p, t_p, rhoi_p, nm_p, ni_p, ti_p) &
          bind(c, name="gw_prof_stage_dispatch_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pver_c
       real(c_double), value :: cpair_c, rair_c, gravit_c
       type(c_ptr), value :: p_ifc_p, p_rdst_p, t_p, rhoi_p, nm_p, ni_p, ti_p
     end subroutine gw_prof_stage_dispatch_codon
  end interface

  call gw_prof_stage_dispatch_codon(int(ncol_local, c_int64_t), int(pver_local, c_int64_t), &
       real(cpair_local, c_double), real(rair_local, c_double), real(gravit_local, c_double), &
       c_loc(p_ifc_local), c_loc(p_rdst_local), c_loc(t_local), c_loc(rhoi_local), &
       c_loc(nm_local), c_loc(ni_local), c_loc(ti_local))

end subroutine gw_prof_codon_wrap

!==========================================================================

subroutine gw_prof_native (ncol, p, cpair, t, rhoi, nm, ni)
  !-----------------------------------------------------------------------
  ! Compute profiles of background state quantities for the multiple
  ! gravity wave drag parameterization.
  !
  ! The parameterization is assumed to operate only where water vapor
  ! concentrations are negligible in determining the density.
  !-----------------------------------------------------------------------
  use gw_utils, only: midpoint_interp
  !------------------------------Arguments--------------------------------
  ! Column dimension.
  integer, intent(in) :: ncol
  ! Pressure coordinates.
  type(Coords1D), intent(in) :: p

  ! Specific heat of dry air, constant pressure.
  real(r8), intent(in) :: cpair
  ! Midpoint temperatures.
  real(r8), intent(in) :: t(ncol,pver)

  ! Interface density.
  real(r8), intent(out) :: rhoi(ncol,pver+1)
  ! Midpoint and interface Brunt-Vaisalla frequencies.
  real(r8), intent(out) :: nm(ncol,pver), ni(ncol,pver+1)

  !---------------------------Local Storage-------------------------------
  ! Column and level indices.
  integer :: i,k

  ! dt/dp
  real(r8) :: dtdp
  ! Brunt-Vaisalla frequency squared.
  real(r8) :: n2

  ! Interface temperature.
  real(r8) :: ti(ncol,pver+1)

  ! Minimum value of Brunt-Vaisalla frequency squared.
  real(r8), parameter :: n2min = 5.e-5_r8

  !------------------------------------------------------------------------
  ! Determine the interface densities and Brunt-Vaisala frequencies.
  !------------------------------------------------------------------------

  ! The top interface values are calculated assuming an isothermal
  ! atmosphere above the top level.
  k = 1
  do i = 1, ncol
     ti(i,k) = t(i,k)
     rhoi(i,k) = p%ifc(i,k) / (rair*ti(i,k))
     ni(i,k) = sqrt(gravit*gravit / (cpair*ti(i,k)))
  end do

  ! Interior points use centered differences.
  ti(:,2:pver) = midpoint_interp(t)
  do k = 2, pver
     do i = 1, ncol
        rhoi(i,k) = p%ifc(i,k) / (rair*ti(i,k))
        dtdp = (t(i,k)-t(i,k-1)) * p%rdst(i,k-1)
        n2 = gravit*gravit/ti(i,k) * (1._r8/cpair - rhoi(i,k)*dtdp)
        ni(i,k) = sqrt(max(n2min, n2))
     end do
  end do

  ! Bottom interface uses bottom level temperature, density; next interface
  ! B-V frequency.
  k = pver+1
  do i = 1, ncol
     ti(i,k) = t(i,k-1)
     rhoi(i,k) = p%ifc(i,k) / (rair*ti(i,k))
     ni(i,k) = ni(i,k-1)
  end do

  !------------------------------------------------------------------------
  ! Determine the midpoint Brunt-Vaisala frequencies.
  !------------------------------------------------------------------------
  nm = midpoint_interp(ni)

end subroutine gw_prof_native

!==========================================================================

subroutine gw_drag_prof(ncol, band, p, src_level, tend_level, dt, &
     t,    &
     piln, rhoi,    nm,   ni,  ubm,  ubi,  xv,    yv,   &
     effgw,      c, kvtt, q,   dse,  tau,  utgw,  vtgw, &
     ttgw, qtgw, egwdffi,   gwut, dttdf, dttke, ro_adjust)

  !-----------------------------------------------------------------------
  ! Solve for the drag profile from the multiple gravity wave drag
  ! parameterization.
  ! 1. scan up from the wave source to determine the stress profile
  ! 2. scan down the stress profile to determine the tendencies
  !     => apply bounds to the tendency
  !          a. from wkb solution
  !          b. from computational stability constraints
  !     => adjust stress on interface below to reflect actual bounded
  !        tendency
  !-----------------------------------------------------------------------

  use gw_diffusion, only: gw_ediff, gw_diff_tend
  use linear_1d_operators, only: TriDiagDecomp

  !------------------------------Arguments--------------------------------
  ! Column dimension.
  integer, intent(in) :: ncol
  ! Wavelengths.
  type(GWBand), intent(in) :: band
  ! Pressure coordinates.
  type(Coords1D), intent(in) :: p
  ! Level from which gravity waves are propagated upward.
  integer, intent(in) :: src_level(ncol)
  ! Lowest level where wind tendencies are calculated.
  integer, intent(in) :: tend_level(ncol)
  ! Using tend_level > src_level allows the orographic waves to prescribe
  ! wave propagation up to a certain level, but then allow wind tendencies
  ! and adjustments to tau below that level.

  ! Time step.
  real(r8), intent(in) :: dt

  ! Midpoint and interface temperatures.
  real(r8), intent(in) :: t(ncol,pver)
  ! Log of interface pressures.
  real(r8), intent(in) :: piln(ncol,pver+1)
  ! Interface densities.
  real(r8), intent(in) :: rhoi(ncol,pver+1)
  ! Midpoint and interface Brunt-Vaisalla frequencies.
  real(r8), intent(in) :: nm(ncol,pver), ni(ncol,pver+1)
  ! Projection of wind at midpoints and interfaces.
  real(r8), intent(in) :: ubm(ncol,pver), ubi(ncol,pver+1)
  ! Unit vectors of source wind (zonal and meridional components).
  real(r8), intent(in) :: xv(ncol), yv(ncol)
  ! Tendency efficiency.
  real(r8), intent(in) :: effgw(ncol)
  ! Wave phase speeds for each column.
  real(r8), intent(in) :: c(ncol,-band%ngwv:band%ngwv)
  ! Molecular thermal diffusivity.
  real(r8), intent(in) :: kvtt(ncol,pver+1)
  ! Constituent array.
  real(r8), intent(in) :: q(:,:,:)
  ! Dry static energy.
  real(r8), intent(in) :: dse(ncol,pver)

  ! Wave Reynolds stress.
  real(r8), intent(inout) :: tau(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Zonal/meridional wind tendencies.
  real(r8), intent(out) :: utgw(ncol,pver), vtgw(ncol,pver)
  ! Gravity wave heating tendency.
  real(r8), intent(out) :: ttgw(ncol,pver)
  ! Gravity wave constituent tendency.
  real(r8), intent(out) :: qtgw(:,:,:)

  ! Effective gravity wave diffusivity at interfaces.
  real(r8), intent(out) :: egwdffi(ncol,pver+1)

  ! Gravity wave wind tendency for each wave.
  real(r8), intent(out) :: gwut(ncol,pver,-band%ngwv:band%ngwv)

  ! Temperature tendencies from diffusion and kinetic energy.
  real(r8), intent(out) :: dttdf(ncol,pver)
  real(r8), intent(out) :: dttke(ncol,pver)

  ! Adjustment parameter for IGWs.
  real(r8), intent(in), optional :: &
       ro_adjust(ncol,-band%ngwv:band%ngwv,pver+1)

  !---------------------------Local storage-------------------------------

  ! Level, wavenumber, and constituent loop indices.
  integer :: k, l, m

  ! Lowest tendency and source levels.
  integer :: kbot_tend, kbot_src

  ! "Total" and saturation diffusivity.
  real(r8) :: d(ncol)
  ! Imaginary part of vertical wavenumber.
  real(r8) :: mi(ncol)
  ! Stress after damping.
  real(r8) :: taudmp(ncol)
  ! Saturation stress.
  real(r8) :: tausat(ncol)
  ! (ub-c) and (ub-c)**2
  real(r8) :: ubmc(ncol), ubmc2(ncol)
  ! Temporary ubar tendencies (overall, and at wave l).
  real(r8) :: ubt(ncol,pver), ubtl(ncol)
  real(r8) :: wrk(ncol)
  ! Ratio used for ubt tndmax limiting.
  real(r8) :: ubt_lim_ratio(ncol)
  ! Workspaces for Codon diffusion solver path.
  real(r8), target :: egwdffm(ncol,pver)
  real(r8), target :: egwdff_lev(ncol)
  real(r8), target :: dpidz_sq(ncol,pver+1)
  real(r8), target :: gw_diff_coef(ncol,pver+1)
  real(r8), target :: gw_diff_qnew(ncol,pver)
  real(r8), target :: gw_diff_spr(ncol,pver)
  real(r8), target :: gw_diff_sub(ncol,pver)
  real(r8), target :: gw_diff_diag(ncol,pver)
  real(r8), target :: gw_diff_ca(ncol,pver)
  real(r8), target :: gw_diff_ze(ncol,pver)
  real(r8), target :: gw_diff_dnom(ncol,pver)
  real(r8), target :: gw_diff_zf(ncol,pver)
  ! Whether to keep this core in native Fortran.
  logical :: use_native_core
  ! Whether to keep the diffusion solver cluster in native Fortran.
  logical :: use_native_diff_solver

  ! LU decomposition.
  type(TriDiagDecomp) :: decomp

  !------------------------------------------------------------------------

  ! Lowest levels that loops need to iterate over.
  kbot_tend = maxval(tend_level)
  kbot_src = maxval(src_level)

  call gw_drag_prof_core_select_impl()
  use_native_core = use_native_gw_drag_prof_core_impl .or. present(ro_adjust)
  call gw_diff_solver_select_impl()
  use_native_diff_solver = use_native_gw_diff_solver_impl .or. present(ro_adjust)

  if (use_native_core) then

  ! Initialize gravity wave drag tendencies to zero.

  utgw = 0._r8
  vtgw = 0._r8

  gwut = 0._r8

  dttke = 0._r8
  ttgw = 0._r8

  ! Workaround floating point exception issues on Intel by initializing
  ! everything that's first set in a where block.
  mi = 0._r8
  taudmp = 0._r8
  tausat = 0._r8
  ubmc = 0._r8
  ubmc2 = 0._r8
  wrk = 0._r8

  !------------------------------------------------------------------------
  ! Compute the stress profiles and diffusivities
  !------------------------------------------------------------------------

  ! Loop from bottom to top to get stress profiles.
  do k = kbot_src, ktop, -1

     ! Determine the diffusivity for each column.

     d = dback + kvtt(:,k)

     do l = -band%ngwv, band%ngwv

        ! Determine the absolute value of the saturation stress.
        ! Define critical levels where the sign of (u-c) changes between
        ! interfaces.
        ubmc = ubi(:,k) - c(:,l)

        tausat = 0.0_r8
        where (src_level >= k)
           ! Test to see if u-c has the same sign here as the level below.
           where (ubmc > 0.0_r8 .eqv. ubi(:,k+1) > c(:,l))
              tausat = abs(band%effkwv * rhoi(:,k) * ubmc**3 / &
                   (2._r8*ni(:,k)))
           end where
        end where

        if (present(ro_adjust)) then
           where (src_level >= k)
              tausat = tausat * sqrt(ro_adjust(:,l,k))
           end where
        end if

        where (src_level >= k)

           ! Compute stress for each wave. The stress at this level is the
           ! min of the saturation stress and the stress at the level below
           ! reduced by damping. The sign of the stress must be the same as
           ! at the level below.

           ubmc2 = max(ubmc**2, ubmc2mn)
           mi = ni(:,k) / (2._r8 * band%kwv * ubmc2) * &
                (alpha(k) + ni(:,k)**2/ubmc2 * d)
           wrk = -2._r8*mi*rog*t(:,k)*(piln(:,k+1) - piln(:,k))

           taudmp = tau(:,l,k+1) * exp(wrk)

           ! For some reason, PGI 14.1 loses bit-for-bit reproducibility if
           ! we limit tau, so instead limit the arrays used to set it.
           where (tausat <= taumin) tausat = 0._r8
           where (taudmp <= taumin) taudmp = 0._r8

           tau(:,l,k) = min(taudmp, tausat)

        end where
     end do

  end do

  ! Force tau at the top of the model to zero, if requested.
  if (tau_0_ubc) tau(:,:,ktop) = 0._r8

  ! Apply efficiency to completed stress profile.
  do k = ktop, kbot_tend+1
     do l = -band%ngwv, band%ngwv
        where (k-1 <= tend_level)
           tau(:,l,k) = tau(:,l,k) * effgw
        end where
     end do
  end do

  !------------------------------------------------------------------------
  ! Compute the tendencies from the stress divergence.
  !------------------------------------------------------------------------

  ! Loop over levels from top to bottom
  do k = ktop, kbot_tend

     ! Accumulate the mean wind tendency over wavenumber.
     ubt(:,k) = 0.0_r8

     do l = -band%ngwv, band%ngwv    ! loop over wave

        ! Determine the wind tendency, including excess stress carried down
        ! from above.
        ubtl = gravit * (tau(:,l,k+1)-tau(:,l,k)) * p%rdel(:,k)

        ! Apply first tendency limit to maintain numerical stability.
        ! Enforce du/dt < |c-u|/dt  so u-c cannot change sign
        !    (u^n+1 = u^n + du/dt * dt)
        ! The limiter is somewhat stricter, so that we don't come anywhere
        ! near reversing c-u.
        ubtl = min(ubtl, umcfac * abs(c(:,l)-ubm(:,k)) / dt)

        where (k <= tend_level)

           ! Save tendency for each wave (for later computation of kzz):
           gwut(:,k,l) = sign(ubtl, c(:,l)-ubm(:,k))
           ubt(:,k) = ubt(:,k) + gwut(:,k,l)

        end where

     end do

     ! Apply second tendency limit to maintain numerical stability.
     ! Enforce du/dt < tndmax so that ridicuously large tendencies are not
     ! permitted.
     ! This can only happen above tend_level, so don't bother checking the
     ! level explicitly.
     where (abs(ubt(:,k)) > tndmax)
        ubt_lim_ratio = tndmax/abs(ubt(:,k))
        ubt(:,k) = ubt_lim_ratio * ubt(:,k)
     elsewhere
        ubt_lim_ratio = 1._r8
     end where

     do l = -band%ngwv, band%ngwv
        gwut(:,k,l) = ubt_lim_ratio*gwut(:,k,l)
        ! Redetermine the effective stress on the interface below from the
        ! wind tendency. If the wind tendency was limited above, then the
        ! new stress will be smaller than the old stress, causing stress
        ! divergence in the next layer down. This smoothes large stress
        ! divergences downward while conserving total stress.
        where (k <= tend_level)
           tau(:,l,k+1) = tau(:,l,k) + &
                abs(gwut(:,k,l)) * p%del(:,k) / gravit
        end where
     end do

     ! Project the mean wind tendency onto the components.
     where (k <= tend_level)
        utgw(:,k) = ubt(:,k) * xv
        vtgw(:,k) = ubt(:,k) * yv
     end where

     ! End of level loop.
  end do

  else

     call gw_drag_prof_core_note_entered()
     call gw_drag_prof_core_codon_wrap(1, ncol, pver, pver+1, band%ngwv, ktop, kbot_tend, kbot_src, &
          merge(1, 0, tau_0_ubc), dback, taumin, tndmax, umcfac, ubmc2mn, &
          band%effkwv, band%kwv, gravit, rog, dt, alpha, p%del, p%rdel, &
          t, piln, rhoi, ni, ubm, ubi, xv, yv, effgw, c, kvtt, src_level, tend_level, tau, &
          utgw, vtgw, ttgw, gwut, dttdf, dttke, d, mi, taudmp, tausat, ubmc, ubmc2, &
          ubt, ubtl, wrk, ubt_lim_ratio)

  end if

  if (use_native_diff_solver) then

     ! Calculate effective diffusivity and LU decomposition for the
     ! vertical diffusion solver.
     call gw_ediff (ncol, pver, band%ngwv, kbot_tend, ktop, tend_level, &
          gwut, ubm, nm, rhoi, dt, gravit, p, c, &
          egwdffi, decomp, ro_adjust=ro_adjust)

     ! Calculate tendency on each constituent.
     do m = 1, size(q,3)

        call gw_diff_tend(ncol, pver, kbot_tend, ktop, q(:,:,m), &
             dt, decomp, qtgw(:,:,m))

     enddo

     ! Calculate tendency from diffusing dry static energy (dttdf).
     call gw_diff_tend(ncol, pver, kbot_tend, ktop, dse, dt, decomp, dttdf)

  else

     call gw_diff_solver_note_entered()
     call gw_diff_solver_codon_wrap(1, ncol, pver, pver+1, size(q,3), band%ngwv, &
          kbot_tend, ktop, dt, gravit, gwut, ubm, nm, rhoi, c, tend_level, &
          p%del, p%rdel, p%rdst, q, dse, egwdffi, qtgw, dttdf, &
          egwdffm, egwdff_lev, dpidz_sq, gw_diff_coef, gw_diff_qnew, &
          gw_diff_spr, gw_diff_sub, gw_diff_diag, gw_diff_ca, gw_diff_ze, &
          gw_diff_dnom, gw_diff_zf)

     call gw_diff_solver_codon_wrap(2, ncol, pver, pver+1, size(q,3), band%ngwv, &
          kbot_tend, ktop, dt, gravit, gwut, ubm, nm, rhoi, c, tend_level, &
          p%del, p%rdel, p%rdst, q, dse, egwdffi, qtgw, dttdf, &
          egwdffm, egwdff_lev, dpidz_sq, gw_diff_coef, gw_diff_qnew, &
          gw_diff_spr, gw_diff_sub, gw_diff_diag, gw_diff_ca, gw_diff_ze, &
          gw_diff_dnom, gw_diff_zf)

  end if

  ! Evaluate second temperature tendency term: Conversion of kinetic
  ! energy into thermal.
  if (use_native_core) then
     do l = -band%ngwv, band%ngwv
        do k = ktop, kbot_tend
           dttke(:,k) = dttke(:,k) - (ubm(:,k) - c(:,l)) * gwut(:,k,l)
        end do
     end do

     ttgw = dttke + dttdf
  else
     call gw_drag_prof_core_codon_wrap(2, ncol, pver, pver+1, band%ngwv, ktop, kbot_tend, kbot_src, &
          merge(1, 0, tau_0_ubc), dback, taumin, tndmax, umcfac, ubmc2mn, &
          band%effkwv, band%kwv, gravit, rog, dt, alpha, p%del, p%rdel, &
          t, piln, rhoi, ni, ubm, ubi, xv, yv, effgw, c, kvtt, src_level, tend_level, tau, &
          utgw, vtgw, ttgw, gwut, dttdf, dttke, d, mi, taudmp, tausat, ubmc, ubmc2, &
          ubt, ubtl, wrk, ubt_lim_ratio)
  end if

  ! Deallocate decomp.
  if (use_native_diff_solver) call decomp%finalize()

end subroutine gw_drag_prof

!==========================================================================

subroutine gw_drag_prof_core_append_proof(proof_line)

  character(len=*), intent(in) :: proof_line

  character(len=512) :: proof_file
  integer :: status, n, unitno

  proof_file = ''
  call get_environment_variable('GW_DRAG_PROF_CORE_PROOF_FILE', value=proof_file, length=n, status=status)
  if (status == 0 .and. n > 0) then
     open(newunit=unitno, file=trim(proof_file(:n)), status='unknown', position='append', action='write')
     write(unitno,'(A)') trim(proof_line)
     close(unitno)
  end if

end subroutine gw_drag_prof_core_append_proof

!==========================================================================

subroutine gw_drag_prof_core_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (gw_drag_prof_core_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('GW_DRAG_PROF_CORE_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_gw_drag_prof_core_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_gw_drag_prof_core_impl = .false.
  end if

  gw_drag_prof_core_impl_selected = .true.

  if (masterproc) then
     if (use_native_gw_drag_prof_core_impl) then
        write(iulog,*) 'gw_drag_prof_core implementation = native'
        call gw_drag_prof_core_append_proof('gw_drag_prof_core selector entered implementation = native')
     else
        write(iulog,*) 'gw_drag_prof_core implementation = codon'
        call gw_drag_prof_core_append_proof('gw_drag_prof_core selector entered implementation = codon')
     end if
     call flush(iulog)
  end if

end subroutine gw_drag_prof_core_select_impl

!==========================================================================

subroutine gw_drag_prof_core_note_entered()

  if (gw_drag_prof_core_entered_logged) return
  gw_drag_prof_core_entered_logged = .true.

  if (masterproc) then
     write(iulog,*) 'gw_drag_prof_core entered (unified gw-common dispatch core loops direct = codon)'
     call gw_drag_prof_core_append_proof('gw_drag_prof_core entered (unified gw-common dispatch core loops direct = codon)')
     call flush(iulog)
  end if

end subroutine gw_drag_prof_core_note_entered

!==========================================================================

subroutine gw_drag_prof_core_codon_wrap(stage, ncol_local, pver_local, pverp_local, ngwv_local, &
     ktop_local, kbot_tend_local, kbot_src_local, tau_0_ubc_local, dback_local, taumin_local, &
     tndmax_local, umcfac_local, ubmc2mn_local, effkwv_local, kwv_local, gravit_local, rog_local, &
     dt_local, alpha_local, p_del_local, p_rdel_local, t_local, piln_local, rhoi_local, ni_local, &
     ubm_local, ubi_local, xv_local, yv_local, effgw_local, c_local, kvtt_local, src_level_local, &
     tend_level_local, tau_local, utgw_local, vtgw_local, ttgw_local, gwut_local, dttdf_local, &
     dttke_local, d_local, mi_local, taudmp_local, tausat_local, ubmc_local, ubmc2_local, &
     ubt_local, ubtl_local, wrk_local, ubt_lim_ratio_local)

  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

  integer, intent(in) :: stage
  integer, intent(in) :: ncol_local, pver_local, pverp_local, ngwv_local
  integer, intent(in) :: ktop_local, kbot_tend_local, kbot_src_local, tau_0_ubc_local
  real(r8), intent(in) :: dback_local, taumin_local, tndmax_local, umcfac_local, ubmc2mn_local
  real(r8), intent(in) :: effkwv_local, kwv_local, gravit_local, rog_local, dt_local
  real(r8), target, intent(in) :: alpha_local(pverp_local)
  real(r8), target, intent(in) :: p_del_local(ncol_local,pver_local), p_rdel_local(ncol_local,pver_local)
  real(r8), target, intent(in) :: t_local(ncol_local,pver_local), piln_local(ncol_local,pverp_local)
  real(r8), target, intent(in) :: rhoi_local(ncol_local,pverp_local), ni_local(ncol_local,pverp_local)
  real(r8), target, intent(in) :: ubm_local(ncol_local,pver_local), ubi_local(ncol_local,pverp_local)
  real(r8), target, intent(in) :: xv_local(ncol_local), yv_local(ncol_local), effgw_local(ncol_local)
  real(r8), target, intent(in) :: c_local(ncol_local,-ngwv_local:ngwv_local)
  real(r8), target, intent(in) :: kvtt_local(ncol_local,pverp_local)
  integer, intent(in) :: src_level_local(ncol_local), tend_level_local(ncol_local)
  real(r8), target, intent(inout) :: tau_local(ncol_local,-ngwv_local:ngwv_local,pverp_local)
  real(r8), target, intent(inout) :: utgw_local(ncol_local,pver_local), vtgw_local(ncol_local,pver_local)
  real(r8), target, intent(inout) :: ttgw_local(ncol_local,pver_local)
  real(r8), target, intent(inout) :: gwut_local(ncol_local,pver_local,-ngwv_local:ngwv_local)
  real(r8), target, intent(in) :: dttdf_local(ncol_local,pver_local)
  real(r8), target, intent(inout) :: dttke_local(ncol_local,pver_local)
  real(r8), target, intent(inout) :: d_local(ncol_local), mi_local(ncol_local), taudmp_local(ncol_local)
  real(r8), target, intent(inout) :: tausat_local(ncol_local), ubmc_local(ncol_local), ubmc2_local(ncol_local)
  real(r8), target, intent(inout) :: ubt_local(ncol_local,pver_local), ubtl_local(ncol_local)
  real(r8), target, intent(inout) :: wrk_local(ncol_local), ubt_lim_ratio_local(ncol_local)

  integer(c_int64_t), target :: src_level_i8(ncol_local), tend_level_i8(ncol_local)
  integer :: i

  interface
     subroutine gw_drag_prof_core_codon(stage_c, ncol_c, pver_c, pverp_c, ngwv_c, &
          ktop_c, kbot_tend_c, kbot_src_c, tau_0_ubc_c, dback_c, taumin_c, tndmax_c, &
          umcfac_c, ubmc2mn_c, effkwv_c, kwv_c, gravit_c, rog_c, dt_c, alpha_p, &
          p_del_p, p_rdel_p, t_p, piln_p, rhoi_p, ni_p, ubm_p, ubi_p, xv_p, yv_p, &
          effgw_p, c_p, kvtt_p, src_level_p, tend_level_p, tau_p, utgw_p, vtgw_p, &
          ttgw_p, gwut_p, dttdf_p, dttke_p, d_p, mi_p, taudmp_p, tausat_p, ubmc_p, &
          ubmc2_p, ubt_p, ubtl_p, wrk_p, ubt_lim_ratio_p) &
          bind(c, name="gw_drag_prof_core_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: stage_c, ncol_c, pver_c, pverp_c, ngwv_c
       integer(c_int64_t), value :: ktop_c, kbot_tend_c, kbot_src_c, tau_0_ubc_c
       real(c_double), value :: dback_c, taumin_c, tndmax_c, umcfac_c, ubmc2mn_c
       real(c_double), value :: effkwv_c, kwv_c, gravit_c, rog_c, dt_c
       type(c_ptr), value :: alpha_p, p_del_p, p_rdel_p, t_p, piln_p, rhoi_p, ni_p
       type(c_ptr), value :: ubm_p, ubi_p, xv_p, yv_p, effgw_p, c_p, kvtt_p
       type(c_ptr), value :: src_level_p, tend_level_p, tau_p, utgw_p, vtgw_p, ttgw_p, gwut_p
       type(c_ptr), value :: dttdf_p, dttke_p, d_p, mi_p, taudmp_p, tausat_p, ubmc_p, ubmc2_p
       type(c_ptr), value :: ubt_p, ubtl_p, wrk_p, ubt_lim_ratio_p
     end subroutine gw_drag_prof_core_codon
  end interface

  do i = 1, ncol_local
     src_level_i8(i) = int(src_level_local(i), c_int64_t)
     tend_level_i8(i) = int(tend_level_local(i), c_int64_t)
  end do

  call gw_drag_prof_core_codon(int(stage, c_int64_t), int(ncol_local, c_int64_t), &
       int(pver_local, c_int64_t), int(pverp_local, c_int64_t), int(ngwv_local, c_int64_t), &
       int(ktop_local, c_int64_t), int(kbot_tend_local, c_int64_t), int(kbot_src_local, c_int64_t), &
       int(tau_0_ubc_local, c_int64_t), real(dback_local, c_double), real(taumin_local, c_double), &
       real(tndmax_local, c_double), real(umcfac_local, c_double), real(ubmc2mn_local, c_double), &
       real(effkwv_local, c_double), real(kwv_local, c_double), real(gravit_local, c_double), &
       real(rog_local, c_double), real(dt_local, c_double), c_loc(alpha_local), c_loc(p_del_local), &
       c_loc(p_rdel_local), c_loc(t_local), c_loc(piln_local), c_loc(rhoi_local), c_loc(ni_local), &
       c_loc(ubm_local), c_loc(ubi_local), c_loc(xv_local), c_loc(yv_local), c_loc(effgw_local), &
       c_loc(c_local), c_loc(kvtt_local), c_loc(src_level_i8), c_loc(tend_level_i8), c_loc(tau_local), &
       c_loc(utgw_local), c_loc(vtgw_local), c_loc(ttgw_local), c_loc(gwut_local), c_loc(dttdf_local), &
       c_loc(dttke_local), c_loc(d_local), c_loc(mi_local), c_loc(taudmp_local), c_loc(tausat_local), &
       c_loc(ubmc_local), c_loc(ubmc2_local), c_loc(ubt_local), c_loc(ubtl_local), c_loc(wrk_local), &
       c_loc(ubt_lim_ratio_local))

end subroutine gw_drag_prof_core_codon_wrap

!==========================================================================

subroutine gw_diff_solver_append_proof(proof_line)

  character(len=*), intent(in) :: proof_line

  character(len=512) :: proof_file
  integer :: status, n, unitno

  proof_file = ''
  call get_environment_variable('GW_DIFF_SOLVER_PROOF_FILE', value=proof_file, length=n, status=status)
  if (status == 0 .and. n > 0) then
     open(newunit=unitno, file=trim(proof_file(:n)), status='unknown', position='append', action='write')
     write(unitno,'(A)') trim(proof_line)
     close(unitno)
  end if

end subroutine gw_diff_solver_append_proof

!==========================================================================

subroutine gw_diff_solver_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (gw_diff_solver_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('GW_DIFF_SOLVER_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_gw_diff_solver_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_gw_diff_solver_impl = .false.
  end if

  gw_diff_solver_impl_selected = .true.

  if (masterproc) then
     if (use_native_gw_diff_solver_impl) then
        write(iulog,*) 'gw_diff_solver implementation = native'
        call gw_diff_solver_append_proof('gw_diff_solver selector entered implementation = native')
     else
        write(iulog,*) 'gw_diff_solver implementation = codon'
        call gw_diff_solver_append_proof('gw_diff_solver selector entered implementation = codon')
     end if
     call flush(iulog)
  end if

end subroutine gw_diff_solver_select_impl

!==========================================================================

subroutine gw_diff_solver_note_entered()

  if (gw_diff_solver_entered_logged) return
  gw_diff_solver_entered_logged = .true.

  if (masterproc) then
     write(iulog,*) 'gw_diff_solver entered (unified gw-common dispatch diff solver direct = codon)'
     call gw_diff_solver_append_proof('gw_diff_solver entered (unified gw-common dispatch diff solver direct = codon)')
     call flush(iulog)
  end if

end subroutine gw_diff_solver_note_entered

!==========================================================================

subroutine gw_diff_solver_codon_wrap(stage, ncol_local, pver_local, pverp_local, pcnst_local, ngwv_local, &
     kbot_local, ktop_local, dt_local, gravit_local, gwut_local, ubm_local, nm_local, rho_local, c_local, &
     tend_level_local, p_del_local, p_rdel_local, p_rdst_local, q_local, dse_local, egwdffi_local, &
     qtgw_local, dttdf_local, egwdffm_local, egwdff_lev_local, dpidz_sq_local, coef_q_diff_local, &
     qnew_local, spr_local, sub_local, diag_local, ca_local, ze_local, dnom_local, zf_local)

  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

  integer, intent(in) :: stage
  integer, intent(in) :: ncol_local, pver_local, pverp_local, pcnst_local, ngwv_local
  integer, intent(in) :: kbot_local, ktop_local
  real(r8), intent(in) :: dt_local, gravit_local
  real(r8), target, intent(in) :: gwut_local(ncol_local,pver_local,-ngwv_local:ngwv_local)
  real(r8), target, intent(in) :: ubm_local(ncol_local,pver_local), nm_local(ncol_local,pver_local)
  real(r8), target, intent(in) :: rho_local(ncol_local,pverp_local)
  real(r8), target, intent(in) :: c_local(ncol_local,-ngwv_local:ngwv_local)
  integer, intent(in) :: tend_level_local(ncol_local)
  real(r8), target, intent(in) :: p_del_local(ncol_local,pver_local)
  real(r8), target, intent(in) :: p_rdel_local(ncol_local,pver_local)
  real(r8), target, intent(in) :: p_rdst_local(ncol_local,pver_local-1)
  real(r8), target, intent(in) :: q_local(ncol_local,pver_local,pcnst_local)
  real(r8), target, intent(in) :: dse_local(ncol_local,pver_local)
  real(r8), target, intent(inout) :: egwdffi_local(ncol_local,pverp_local)
  real(r8), target, intent(inout) :: qtgw_local(ncol_local,pver_local,pcnst_local)
  real(r8), target, intent(inout) :: dttdf_local(ncol_local,pver_local)
  real(r8), target, intent(inout) :: egwdffm_local(ncol_local,pver_local)
  real(r8), target, intent(inout) :: egwdff_lev_local(ncol_local)
  real(r8), target, intent(inout) :: dpidz_sq_local(ncol_local,pverp_local)
  real(r8), target, intent(inout) :: coef_q_diff_local(ncol_local,pverp_local)
  real(r8), target, intent(inout) :: qnew_local(ncol_local,pver_local)
  real(r8), target, intent(inout) :: spr_local(ncol_local,pver_local)
  real(r8), target, intent(inout) :: sub_local(ncol_local,pver_local)
  real(r8), target, intent(inout) :: diag_local(ncol_local,pver_local)
  real(r8), target, intent(inout) :: ca_local(ncol_local,pver_local)
  real(r8), target, intent(inout) :: ze_local(ncol_local,pver_local)
  real(r8), target, intent(inout) :: dnom_local(ncol_local,pver_local)
  real(r8), target, intent(inout) :: zf_local(ncol_local,pver_local)

  integer(c_int64_t), target :: tend_level_i8(ncol_local)
  integer :: i

  interface
     subroutine gw_diff_solver_codon(stage_c, ncol_c, pver_c, pverp_c, pcnst_c, ngwv_c, &
          kbot_c, ktop_c, dt_c, gravit_c, gwut_p, ubm_p, nm_p, rho_p, c_p, tend_level_p, &
          p_del_p, p_rdel_p, p_rdst_p, q_p, dse_p, egwdffi_p, qtgw_p, dttdf_p, &
          egwdffm_p, egwdff_lev_p, dpidz_sq_p, coef_q_diff_p, qnew_p, spr_p, sub_p, &
          diag_p, ca_p, ze_p, dnom_p, zf_p) &
          bind(c, name="gw_diff_solver_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: stage_c, ncol_c, pver_c, pverp_c, pcnst_c, ngwv_c
       integer(c_int64_t), value :: kbot_c, ktop_c
       real(c_double), value :: dt_c, gravit_c
       type(c_ptr), value :: gwut_p, ubm_p, nm_p, rho_p, c_p, tend_level_p
       type(c_ptr), value :: p_del_p, p_rdel_p, p_rdst_p, q_p, dse_p, egwdffi_p
       type(c_ptr), value :: qtgw_p, dttdf_p, egwdffm_p, egwdff_lev_p, dpidz_sq_p
       type(c_ptr), value :: coef_q_diff_p, qnew_p, spr_p, sub_p, diag_p, ca_p, ze_p, dnom_p, zf_p
     end subroutine gw_diff_solver_codon
  end interface

  do i = 1, ncol_local
     tend_level_i8(i) = int(tend_level_local(i), c_int64_t)
  end do

  call gw_diff_solver_codon(int(stage, c_int64_t), int(ncol_local, c_int64_t), &
       int(pver_local, c_int64_t), int(pverp_local, c_int64_t), int(pcnst_local, c_int64_t), &
       int(ngwv_local, c_int64_t), int(kbot_local, c_int64_t), int(ktop_local, c_int64_t), &
       real(dt_local, c_double), real(gravit_local, c_double), c_loc(gwut_local), c_loc(ubm_local), &
       c_loc(nm_local), c_loc(rho_local), c_loc(c_local), c_loc(tend_level_i8), c_loc(p_del_local), &
       c_loc(p_rdel_local), c_loc(p_rdst_local), c_loc(q_local), c_loc(dse_local), c_loc(egwdffi_local), &
       c_loc(qtgw_local), c_loc(dttdf_local), c_loc(egwdffm_local), c_loc(egwdff_lev_local), &
       c_loc(dpidz_sq_local), c_loc(coef_q_diff_local), c_loc(qnew_local), c_loc(spr_local), &
       c_loc(sub_local), c_loc(diag_local), c_loc(ca_local), c_loc(ze_local), c_loc(dnom_local), &
       c_loc(zf_local))

end subroutine gw_diff_solver_codon_wrap

!==========================================================================

! Calculate Reynolds stress for waves propagating in each cardinal
! direction.

function calc_taucd(ncol, ngwv, tend_level, tau, c, xv, yv, ubi) &
     result(taucd)

  ! Column and gravity wave wavenumber dimensions.
  integer, intent(in) :: ncol, ngwv
  ! Lowest level where wind tendencies are calculated.
  integer, intent(in) :: tend_level(:)
  ! Wave Reynolds stress.
  real(r8), intent(in) :: tau(:,-ngwv:,:)
  ! Wave phase speeds for each column.
  real(r8), intent(in) :: c(:,-ngwv:)
  ! Unit vectors of source wind (zonal and meridional components).
  real(r8), intent(in) :: xv(:), yv(:)
  ! Projection of wind at interfaces.
  real(r8), intent(in) :: ubi(:,:)

  real(r8) :: taucd(ncol,pver+1,4)

  ! Indices.
  integer :: i, k, l

  ! ubi at tend_level.
  real(r8) :: ubi_tend(ncol)

  ! Signed wave Reynolds stress.
  real(r8) :: tausg(ncol)

  ! Reynolds stress for waves propagating behind and forward of the wind.
  real(r8) :: taub(ncol)
  real(r8) :: tauf(ncol)

  taucd = 0._r8
  tausg = 0._r8

  ubi_tend = (/ (ubi(i,tend_level(i)+1), i = 1, ncol) /)

  do k = ktop, maxval(tend_level)+1

     taub = 0._r8
     tauf = 0._r8

     do l = -ngwv, ngwv
        where (k-1 <= tend_level)

           tausg = sign(tau(:,l,k), c(:,l)-ubi(:,k))

           where ( c(:,l) < ubi_tend )
              taub = taub + tausg
           elsewhere
              tauf = tauf + tausg
           end where

        end where
     end do

     where (k-1 <= tend_level)
        where (xv > 0._r8)
           taucd(:,k,east) = tauf * xv
           taucd(:,k,west) = taub * xv
        elsewhere
           taucd(:,k,east) = taub * xv
           taucd(:,k,west) = tauf * xv
        end where

        where ( yv > 0._r8)
           taucd(:,k,north) = tauf * yv
           taucd(:,k,south) = taub * yv
        elsewhere
           taucd(:,k,north) = taub * yv
           taucd(:,k,south) = tauf * yv
        end where
     end where

  end do

end function calc_taucd

!==========================================================================

! Calculate the amount of momentum conveyed from below the gravity wave
! region, to the region where gravity waves are calculated.
subroutine momentum_flux(tend_level, taucd, um_flux, vm_flux)

  ! Bottom stress level.
  integer, intent(in) :: tend_level(:)
  ! Projected stresses.
  real(r8), intent(in) :: taucd(:,:,:)
  ! Components of momentum change sourced from the bottom.
  real(r8), intent(out) :: um_flux(:), vm_flux(:)

  integer :: i

  ! Tendency for U & V below source level.
  do i = 1, size(tend_level)
     um_flux(i) = taucd(i,tend_level(i)+1, east) + &
                  taucd(i,tend_level(i)+1, west)
     vm_flux(i) = taucd(i,tend_level(i)+1,north) + &
                  taucd(i,tend_level(i)+1,south)
  end do

end subroutine momentum_flux

!==========================================================================

! Subtracts a change in momentum in the gravity wave levels from wind
! tendencies in lower levels, ensuring momentum conservation.
subroutine momentum_fixer(tend_level, p, um_flux, vm_flux, utgw, vtgw)

  ! Bottom stress level.
  integer, intent(in) :: tend_level(:)
  ! Pressure coordinates.
  type(Coords1D), intent(in) :: p
  ! Components of momentum change sourced from the bottom.
  real(r8), intent(in) :: um_flux(:), vm_flux(:)
  ! Wind tendencies.
  real(r8), intent(inout) :: utgw(:,:), vtgw(:,:)

  ! Indices.
  integer :: i, k
  ! Reciprocal of total mass.
  real(r8) :: rdm(size(tend_level))
  ! Average changes in velocity from momentum change being spread over
  ! total mass.
  real(r8) :: du(size(tend_level)), dv(size(tend_level))

  ! Total mass from ground to source level: rho*dz = dp/gravit
  do i = 1, size(tend_level)
     rdm(i) = gravit/(p%ifc(i,pver+1)-p%ifc(i,tend_level(i)+1))
  end do

  ! Average velocity changes.
  du = -um_flux*rdm
  dv = -vm_flux*rdm

  do k = minval(tend_level)+1, pver
     where (k > tend_level)
        utgw(:,k) = utgw(:,k) + du
        vtgw(:,k) = vtgw(:,k) + dv
     end where
  end do
  
end subroutine momentum_fixer

!==========================================================================

! Calculate the change in total energy from tendencies up to this point.
subroutine energy_change(dt, p, u, v, dudt, dvdt, dsdt, de)

  ! Selectable wrapper for total-column gravity-wave energy diagnostics.
  ! The Codon path receives explicit-shape Fortran dummies so array sections
  ! such as ptend%u(:ncol,:) are materialized with the original ncol stride.

  ! Time step.
  real(r8), intent(in) :: dt
  ! Pressure coordinates.
  type(Coords1D), intent(in) :: p
  ! Winds at start of time step.
  real(r8), intent(in) :: u(:,:), v(:,:)
  ! Wind tendencies.
  real(r8), intent(in) :: dudt(:,:), dvdt(:,:)
  ! Heating tendency.
  real(r8), intent(in) :: dsdt(:,:)
  ! Change in energy.
  real(r8), intent(out) :: de(:)

  call energy_change_select_impl()

  if (use_native_energy_change_impl) then
     call energy_change_native(dt, p, u, v, dudt, dvdt, dsdt, de)
  else
     call energy_change_note_entered()
     call energy_change_codon_wrap(size(de), pver, dt, gravit, p%del, u, v, dudt, dvdt, dsdt, de)
  end if

end subroutine energy_change

!==========================================================================

subroutine energy_change_append_proof(proof_line)

  character(len=*), intent(in) :: proof_line

  character(len=512) :: proof_file
  integer :: status, n, unitno

  proof_file = ''
  call get_environment_variable('GW_ENERGY_CHANGE_PROOF_FILE', value=proof_file, length=n, status=status)
  if (status == 0 .and. n > 0) then
     open(newunit=unitno, file=trim(proof_file(:n)), status='unknown', position='append', action='write')
     write(unitno,'(A)') trim(proof_line)
     close(unitno)
  end if

end subroutine energy_change_append_proof

!==========================================================================

subroutine energy_change_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (energy_change_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('GW_ENERGY_CHANGE_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_energy_change_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_energy_change_impl = .false.
  end if

  energy_change_impl_selected = .true.

  if (masterproc) then
     if (use_native_energy_change_impl) then
        write(iulog,*) 'gw_energy_change implementation = native'
        call energy_change_append_proof('gw_energy_change selector entered implementation = native')
     else
        write(iulog,*) 'gw_energy_change implementation = codon'
        call energy_change_append_proof('gw_energy_change selector entered implementation = codon')
     end if
     call flush(iulog)
  end if

end subroutine energy_change_select_impl

!==========================================================================

subroutine energy_change_note_entered()

  if (energy_change_entered_logged) return
  energy_change_entered_logged = .true.

  if (masterproc) then
     write(iulog,*) 'gw_energy_change entered (unified column-energy stage dispatch = codon)'
     call energy_change_append_proof('gw_energy_change entered (unified column-energy stage dispatch = codon)')
     call flush(iulog)
  end if

end subroutine energy_change_note_entered

!==========================================================================

subroutine energy_change_codon_wrap(ncol_local, pver_local, dt_local, gravit_local, &
     p_del_local, u_local, v_local, dudt_local, dvdt_local, dsdt_local, de_local)

  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

  integer, intent(in) :: ncol_local, pver_local
  real(r8), intent(in) :: dt_local, gravit_local
  real(r8), target, intent(in) :: p_del_local(ncol_local,pver_local)
  real(r8), target, intent(in) :: u_local(ncol_local,pver_local), v_local(ncol_local,pver_local)
  real(r8), target, intent(in) :: dudt_local(ncol_local,pver_local), dvdt_local(ncol_local,pver_local)
  real(r8), target, intent(in) :: dsdt_local(ncol_local,pver_local)
  real(r8), target, intent(inout) :: de_local(ncol_local)

  interface
     subroutine gw_energy_change_stage_dispatch_codon(ncol_c, pver_c, dt_c, gravit_c, &
          p_del_p, u_p, v_p, dudt_p, dvdt_p, dsdt_p, de_p) &
          bind(c, name="gw_energy_change_stage_dispatch_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pver_c
       real(c_double), value :: dt_c, gravit_c
       type(c_ptr), value :: p_del_p, u_p, v_p, dudt_p, dvdt_p, dsdt_p, de_p
     end subroutine gw_energy_change_stage_dispatch_codon
  end interface

  call gw_energy_change_stage_dispatch_codon(int(ncol_local, c_int64_t), int(pver_local, c_int64_t), &
       real(dt_local, c_double), real(gravit_local, c_double), &
       c_loc(p_del_local), c_loc(u_local), c_loc(v_local), c_loc(dudt_local), &
       c_loc(dvdt_local), c_loc(dsdt_local), c_loc(de_local))

end subroutine energy_change_codon_wrap

!==========================================================================

subroutine energy_change_native(dt, p, u, v, dudt, dvdt, dsdt, de)

  ! Time step.
  real(r8), intent(in) :: dt
  ! Pressure coordinates.
  type(Coords1D), intent(in) :: p
  ! Winds at start of time step.
  real(r8), intent(in) :: u(:,:), v(:,:)
  ! Wind tendencies.
  real(r8), intent(in) :: dudt(:,:), dvdt(:,:)
  ! Heating tendency.
  real(r8), intent(in) :: dsdt(:,:)
  ! Change in energy.
  real(r8), intent(out) :: de(:)

  ! Level index.
  integer :: k

  ! Net gain/loss of total energy in the column.
  de = 0.0_r8
  do k = 1, pver
     de = de + p%del(:,k)/gravit * (dsdt(:,k) + &
          dudt(:,k)*(u(:,k)+dudt(:,k)*0.5_r8*dt) + &
          dvdt(:,k)*(v(:,k)+dvdt(:,k)*0.5_r8*dt) )
  end do

end subroutine energy_change_native

!==========================================================================

! Subtract change in energy from the heating tendency in the levels below
! the gravity wave region.
subroutine energy_fixer(tend_level, p, de, ttgw)

  ! Bottom stress level.
  integer, intent(in) :: tend_level(:)
  ! Pressure coordinates.
  type(Coords1D), intent(in) :: p
  ! Change in energy.
  real(r8), intent(in) :: de(:)
  ! Heating tendency.
  real(r8), intent(inout) :: ttgw(:,:)

  ! Column/level indices.
  integer :: i, k
  ! Energy change to apply divided by all the mass it is spread across.
  real(r8) :: de_dm(size(tend_level))

  do i = 1, size(tend_level)
     de_dm(i) = -de(i)*gravit/(p%ifc(i,pver+1)-p%ifc(i,tend_level(i)+1))
  end do

  ! Subtract net gain/loss of total energy below tend_level.
  do k = minval(tend_level)+1, pver
     where (k > tend_level)
        ttgw(:,k) = ttgw(:,k) + de_dm
     end where
  end do

end subroutine energy_fixer

!==========================================================================

! Calculates absolute value of the local Coriolis frequency divided by the
! spatial frequency kwv, which gives a characteristic speed in m/s.
function coriolis_speed(band, lat)
  ! Inertial gravity wave lengths.
  type(GWBand), intent(in) :: band
  ! Latitude in radians.
  real(r8), intent(in) :: lat(:)

  real(r8) :: coriolis_speed(size(lat))

  ! 24*3600 = 86400 seconds in a day.
  real(r8), parameter :: omega_earth = 2._r8*pi/86400._r8

  coriolis_speed = abs(sin(lat) * 2._r8 * omega_earth / band%kwv)

end function coriolis_speed

!==========================================================================

subroutine adjust_inertial(band, tend_level, &
     u_coriolis, c, ubi, tau, ro_adjust)
  ! Inertial gravity wave lengths.
  type(GWBand), intent(in) :: band
  ! Levels above which tau is calculated.
  integer, intent(in) :: tend_level(:)
  ! Absolute value of the Coriolis frequency for each column,
  ! divided by kwv [m/s].
  real(r8), intent(in) :: u_coriolis(:)
  ! Wave propagation speed.
  real(r8), intent(in) :: c(:,-band%ngwv:)
  ! Wind speed in the direction of wave propagation.
  real(r8), intent(in) :: ubi(:,:)

  ! Tau will be adjusted by blocking wave propagation through cells where
  ! the Coriolis effect prevents it.
  real(r8), intent(inout) :: tau(:,-band%ngwv:,:)
  ! Dimensionless Coriolis term used to reduce gravity wave strength.
  ! Equal to max(0, 1 - (1/ro)^2), where ro is the Rossby number of the
  ! wind with respect to inertial waves.
  real(r8), intent(out) :: ro_adjust(:,-band%ngwv:,:)

  ! Column/level/wavenumber indices.
  integer :: i, k, l

  ! For each column and wavenumber, are we clear of levels that block
  ! upward propagation?
  logical :: unblocked_mask(size(tend_level),-band%ngwv:band%ngwv)

  unblocked_mask = .true.
  ro_adjust = 0._r8

  ! Iterate from the bottom up, through every interface level where tau is
  ! set.
  do k = maxval(tend_level)+1, ktop, -1
     do l = -band%ngwv, band%ngwv
        do i = 1, size(tend_level)
           ! Only operate on valid levels for this column.
           if (k <= tend_level(i) + 1) then
              ! Block waves if Coriolis is too strong.
              ! By setting the mask in this way, we avoid division by zero.
              unblocked_mask(i,l) = unblocked_mask(i,l) .and. &
                   (abs(ubi(i,k) - c(i,l)) > u_coriolis(i))
              if (unblocked_mask(i,l)) then
                 ro_adjust(i,l,k) = &
                      1._r8 - (u_coriolis(i)/(ubi(i,k)-c(i,l)))**2
              else
                 tau(i,l,k) = 0._r8
              end if
           end if
        end do
     end do
  end do

end subroutine adjust_inertial

end module gw_common
