module gw_oro

!
! This module handles gravity waves from orographic sources, and was
! extracted from gw_drag in May 2013.
!
use gw_utils, only: r8
use coords_1d, only: Coords1D

implicit none
private
save

! Public interface
public :: gw_oro_src

contains

!==========================================================================

subroutine gw_oro_src(ncol, band, p, &
     u, v, t, sgh, zm, nm, &
     src_level, tend_level, tau, ubm, ubi, xv, yv, c)
  use gw_common, only: GWBand, pver, rair
  use ap_gw_oro_src_scheme, only: gw_oro_src_run
  use perf_mod, only: t_startf, t_stopf
  !-----------------------------------------------------------------------
  ! Orographic source for multiple gravity wave drag parameterization.
  !
  ! The stress is returned for a single wave with c=0, over orography.
  ! For points where the orographic variance is small (including ocean),
  ! the returned stress is zero.
  !------------------------------Arguments--------------------------------
  ! Column dimension.
  integer, intent(in) :: ncol
  ! Band to emit orographic waves in.
  ! Regardless, we will only ever emit into l = 0.
  type(GWBand), intent(in) :: band
  ! Pressure coordinates.
  type(Coords1D), intent(in) :: p

  ! Midpoint zonal/meridional winds.
  real(r8), intent(in) :: u(ncol,pver), v(ncol,pver)
  ! Midpoint temperatures.
  real(r8), intent(in) :: t(ncol,pver)
  ! Standard deviation of orography.
  real(r8), intent(in) :: sgh(ncol)
  ! Midpoint altitudes.
  real(r8), intent(in) :: zm(ncol,pver)
  ! Midpoint Brunt-Vaisalla frequencies.
  real(r8), intent(in) :: nm(ncol,pver)

  ! Indices of top gravity wave source level and lowest level where wind
  ! tendencies are allowed.
  integer, intent(out) :: src_level(ncol)
  integer, intent(out) :: tend_level(ncol)

  ! Wave Reynolds stress.
  real(r8), intent(out) :: tau(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Projection of wind at midpoints and interfaces.
  real(r8), intent(out) :: ubm(ncol,pver), ubi(ncol,pver+1)
  ! Unit vectors of source wind (zonal and meridional components).
  real(r8), intent(out) :: xv(ncol), yv(ncol)
  ! Phase speeds.
  real(r8), intent(out) :: c(ncol,-band%ngwv:band%ngwv)

  call t_startf('ap_gw_oro_src_run')
  call gw_oro_src_run(pver, pver+1, rair, band%ngwv, 2*band%ngwv+1, &
       band%fcrit2, band%kwv, ncol, p, u, v, t, sgh, zm, nm, &
       src_level, tend_level, tau, ubm, ubi, xv, yv, c)
  call t_stopf('ap_gw_oro_src_run')

end subroutine gw_oro_src

end module gw_oro
