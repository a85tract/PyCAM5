module ap_gw_prof_scheme

  use gw_utils, only: r8
  use coords_1d, only: Coords1D

  implicit none
  private

  public :: gw_prof_run

contains

  !> \section arg_table_gw_prof_run Argument Table
  !! \htmlinclude gw_prof_run.html
  !!
subroutine gw_prof_run(pver, pverp, rair, gravit, ncol, p, cpair, t, rhoi, nm, ni)
  !-----------------------------------------------------------------------
  ! Compute profiles of background state quantities for the multiple
  ! gravity wave drag parameterization.
  !
  ! The parameterization is assumed to operate only where water vapor
  ! concentrations are negligible in determining the density.
  !-----------------------------------------------------------------------
  use gw_utils, only: midpoint_interp
  !------------------------------Arguments--------------------------------
  integer, intent(in) :: pver, pverp
  real(r8), intent(in) :: rair, gravit
  ! Column dimension.
  integer, intent(in) :: ncol
  ! Pressure coordinates.
  type(Coords1D), intent(in) :: p

  ! Specific heat of dry air, constant pressure.
  real(r8), intent(in) :: cpair
  ! Midpoint temperatures.
  real(r8), intent(in) :: t(ncol,pver)

  ! Interface density.
  real(r8), intent(out) :: rhoi(ncol,pverp)
  ! Midpoint and interface Brunt-Vaisalla frequencies.
  real(r8), intent(out) :: nm(ncol,pver), ni(ncol,pverp)

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

end subroutine gw_prof_run

end module ap_gw_prof_scheme
