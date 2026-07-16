module gw_diffusion

!
! This module contains code computing the effective diffusion of
! constituents and dry static energy due to gravity wave breaking.
!

use gw_utils, only: r8
use linear_1d_operators, only: TriDiagDecomp

implicit none
private
save

public :: gw_ediff
public :: gw_ediff_fields
public :: gw_diff_tend

contains

!==========================================================================

subroutine gw_ediff(ncol, pver, ngwv, kbot, ktop, tend_level, &
     gwut, ubm, nm, rho, dt, gravit, p, c, &
     egwdffi, decomp, ro_adjust)
!
! Calculate effective diffusivity associated with GW forcing.
!
! Author: F. Sassi, Jan 31, 2001
!
  use coords_1d, only: Coords1D

!-------------------------------Input Arguments----------------------------

  ! Column, level, and gravity wave spectrum dimensions.
  integer, intent(in) :: ncol, pver, ngwv
  ! Bottom and top levels to operate on.
  integer, intent(in) :: kbot, ktop
  ! Per-column bottom index where tendencies are applied.
  integer, intent(in) :: tend_level(ncol)
  ! GW zonal wind tendencies at midpoint.
  real(r8), intent(in) :: gwut(ncol,pver,-ngwv:ngwv)
  ! Projection of wind at midpoints.
  real(r8), intent(in) :: ubm(ncol,pver)
  ! Brunt-Vaisalla frequency.
  real(r8), intent(in) :: nm(ncol,pver)

  ! Density at interfaces.
  real(r8), intent(in) :: rho(ncol,pver+1)
  ! Time step.
  real(r8), intent(in) :: dt
  ! Acceleration due to gravity.
  real(r8), intent(in) :: gravit
  ! Pressure coordinates.
  type(Coords1D), intent(in) :: p
  ! Wave phase speeds for each column.
  real(r8), intent(in) :: c(ncol,-ngwv:ngwv)

  ! Adjustment parameter for IGWs.
  real(r8), intent(in), optional :: &
       ro_adjust(ncol,-ngwv:ngwv,pver+1)

!-----------------------------Output Arguments-----------------------------
  ! Effective gw diffusivity at interfaces.
  real(r8), intent(out) :: egwdffi(ncol,pver+1)
  ! LU decomposition.
  type(TriDiagDecomp), intent(out) :: decomp

  call gw_ediff_fields(ncol, pver, ngwv, kbot, ktop, tend_level, &
       gwut, ubm, nm, rho, dt, gravit, p%rdel, p%rdst, c, &
       egwdffi, decomp, ro_adjust)

end subroutine gw_ediff

!==========================================================================

subroutine gw_ediff_fields(ncol, pver, ngwv, kbot, ktop, tend_level, &
     gwut, ubm, nm, rho, dt, gravit, p_rdel, p_rdst, c, &
     egwdffi, decomp, ro_adjust)
!
! Calculate effective diffusivity associated with GW forcing from
! intrinsic pressure-coordinate fields.
!
  use gw_utils, only: midpoint_interp
  use coords_1d, only: Coords1D
  use vdiff_lu_solver, only: fin_vol_lu_decomp

!-------------------------------Input Arguments----------------------------

  ! Column, level, and gravity wave spectrum dimensions.
  integer, intent(in) :: ncol, pver, ngwv
  ! Bottom and top levels to operate on.
  integer, intent(in) :: kbot, ktop
  ! Per-column bottom index where tendencies are applied.
  integer, intent(in) :: tend_level(ncol)
  ! GW zonal wind tendencies at midpoint.
  real(r8), intent(in) :: gwut(ncol,pver,-ngwv:ngwv)
  ! Projection of wind at midpoints.
  real(r8), intent(in) :: ubm(ncol,pver)
  ! Brunt-Vaisalla frequency.
  real(r8), intent(in) :: nm(ncol,pver)
  ! Density at interfaces.
  real(r8), intent(in) :: rho(ncol,pver+1)
  ! Time step and acceleration due to gravity.
  real(r8), intent(in) :: dt, gravit
  ! Inverse pressure-layer thicknesses and midpoint distances.
  real(r8), intent(in) :: p_rdel(ncol,pver)
  real(r8), intent(in) :: p_rdst(ncol,pver-1)
  ! Wave phase speeds for each column.
  real(r8), intent(in) :: c(ncol,-ngwv:ngwv)

  ! Adjustment parameter for IGWs.
  real(r8), intent(in), optional :: &
       ro_adjust(ncol,-ngwv:ngwv,pver+1)

!-----------------------------Output Arguments-----------------------------

  ! Effective gw diffusivity at interfaces.
  real(r8), intent(out) :: egwdffi(ncol,pver+1)
  ! LU decomposition.
  type(TriDiagDecomp), intent(out) :: decomp

!-----------------------------Local Workspace------------------------------

  ! Effective gw diffusivity at midpoints.
  real(r8) :: egwdffm(ncol,pver)
  ! Temporary used to hold gw_diffusivity for one level and wavenumber.
  real(r8) :: egwdff_lev(ncol)
  ! (dp/dz)^2 == (gravit*rho)^2
  real(r8) :: dpidz_sq(ncol,pver+1)
  ! Pressure-coordinate section needed by the legacy LU implementation.
  type(Coords1D) :: p_section
  ! Level and wave indices.
  integer :: k, l
  ! Inverse Prandtl number.
  real(r8), parameter :: prndl=0.25_r8

!--------------------------------------------------------------------------

  egwdffi = 0._r8
  egwdffm = 0._r8

  ! Calculate effective diffusivity at midpoints.
  do l = -ngwv, ngwv
     do k = ktop, kbot

        egwdff_lev = &
             prndl * 0.5_r8 * gwut(:,k,l) * (c(:,l)-ubm(:,k)) / nm(:,k)**2

        ! IGWs have a different Prandtl number, and need ro_adjust factor.
        if (present(ro_adjust)) then
           egwdff_lev = egwdff_lev * 4._r8 * ro_adjust(:,l,k)**2
        end if

        egwdffm(:,k) = egwdffm(:,k) + egwdff_lev

     end do
  end do

  ! Interpolate effective diffusivity to interfaces.
  ! Assume zero at top and bottom interfaces.
  egwdffi(:,ktop+1:kbot) = midpoint_interp(egwdffm(:,ktop:kbot))

  ! Do not calculate diffusivities below level where tendencies are
  ! actually allowed.
  do k = ktop+1, kbot
     where (k > tend_level) egwdffi(:,k) = 0.0_r8
  enddo

  ! Calculate (dp/dz)^2.
  dpidz_sq = rho*gravit
  dpidz_sq = dpidz_sq*dpidz_sq

  ! Construct only the pressure fields consumed by diffusion_operator.
  ! This reproduces the slices copied by Coords1D%section without exposing
  ! Coords1D at the scheme boundary.
  p_section%n = ncol
  p_section%d = kbot-ktop+1
  allocate(p_section%rdel(ncol,p_section%d))
  allocate(p_section%rdst(ncol,p_section%d-1))
  p_section%rdel = p_rdel(:,ktop:kbot)
  p_section%rdst = p_rdst(:,ktop:kbot-1)

  ! Decompose the diffusion matrix.
  decomp = fin_vol_lu_decomp(dt, p_section, &
       coef_q_diff=egwdffi(:,ktop:kbot+1)*dpidz_sq(:,ktop:kbot+1))

end subroutine gw_ediff_fields

!==========================================================================

subroutine gw_diff_tend(ncol, pver, kbot, ktop, q, dt, decomp, dq)

!
! Calculates tendencies from effective diffusion due to gravity wave
! breaking.
!
! Method:
! A constituent flux on interfaces is given by:
!
!              rho * (w'q') = rho * Deff qz
!
! where (all evaluated on interfaces):
!
!        rho   = density
!        qz    = constituent vertical gradient
!        Deff  = effective diffusivity
!
! An effective diffusivity is calculated by adding up the diffusivities
! from all waves (see gw_ediff). The tendency is calculated by invoking LU
! decomposition and solving as for a regular diffusion equation.
!
! Author: Sassi - Jan 2001
!--------------------------------------------------------------------------

!---------------------------Input Arguments--------------------------------

  ! Column and level dimensions.
  integer, intent(in) :: ncol, pver
  ! Bottom and top levels to operate on.
  integer, intent(in) :: kbot, ktop

  ! Constituent to diffuse.
  real(r8), intent(in) :: q(ncol,pver)
  ! Time step.
  real(r8), intent(in) :: dt

  ! LU decomposition.
  type(TriDiagDecomp), intent(in) :: decomp

!--------------------------Output Arguments--------------------------------

  ! Constituent tendencies.
  real(r8), intent(out) :: dq(ncol,pver)

!--------------------------Local Workspace---------------------------------

  ! Temporary storage for constituent.
  real(r8) :: qnew(ncol,pver)

!--------------------------------------------------------------------------

  dq   = 0.0_r8
  qnew = q

  call decomp%left_div(qnew(:,ktop:kbot))

  ! Evaluate tendency to be reported back.
  dq = (qnew-q) / dt

end subroutine gw_diff_tend

end module gw_diffusion
