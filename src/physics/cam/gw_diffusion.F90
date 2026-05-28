module gw_diffusion

!
! This module contains code computing the effective diffusion of
! constituents and dry static energy due to gravity wave breaking.
!

use gw_utils, only: r8
use linear_1d_operators, only: TriDiagDecomp
use spmd_utils, only: masterproc
use cam_logfile, only: iulog
use iso_c_binding, only: c_int64_t

implicit none
private
save

public :: gw_ediff
public :: gw_diff_tend

logical :: use_native_gw_ediff_prep_impl = .false.
logical :: gw_ediff_prep_impl_selected = .false.
logical :: gw_ediff_prep_entered_logged = .false.
logical :: use_native_gw_diff_tend_prepost_impl = .false.
logical :: gw_diff_tend_prepost_impl_selected = .false.
logical :: gw_diff_tend_prepost_entered_logged = .false.

interface
  function gw_ediff_codon(stage_c) result(stage_out) bind(c, name="gw_ediff_codon")
    use iso_c_binding, only: c_int64_t
    integer(c_int64_t), value :: stage_c
    integer(c_int64_t) :: stage_out
  end function gw_ediff_codon
  function gw_diff_tend_codon(stage_c) result(stage_out) bind(c, name="gw_diff_tend_codon")
    use iso_c_binding, only: c_int64_t
    integer(c_int64_t), value :: stage_c
    integer(c_int64_t) :: stage_out
  end function gw_diff_tend_codon
end interface

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

!-----------------------------Local Workspace------------------------------

  ! Effective gw diffusivity at midpoints.
  real(r8) :: egwdffm(ncol,pver)
  ! Temporary used to hold gw_diffusivity for one level and wavenumber.
  real(r8) :: egwdff_lev(ncol)
  ! (dp/dz)^2 == (gravit*rho)^2
  real(r8) :: dpidz_sq(ncol,pver+1)
  ! Level and wave indices.
  integer :: k, l
  ! Inverse Prandtl number.
  real(r8), parameter :: prndl=0.25_r8
  ! Density scale height.
  real(r8), parameter :: dscale=7000._r8
  ! Whether to keep the pre-decomposition diffusivity prep native.
  logical :: use_native_prep
  integer(c_int64_t) :: gw_ediff_touch_c

!--------------------------------------------------------------------------

  call gw_ediff_prep_select_impl()
  use_native_prep = use_native_gw_ediff_prep_impl
  if (.not. use_native_prep) gw_ediff_touch_c = gw_ediff_codon(1601_c_int64_t)

  if (use_native_prep) then

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

  else

     call gw_ediff_prep_note_entered()
     call gw_ediff_prep_codon_wrap(ncol, pver, pver+1, ngwv, kbot, ktop, prndl, gravit, &
          gwut, ubm, nm, rho, c, tend_level, egwdffi, egwdffm, egwdff_lev, dpidz_sq, ro_adjust)

  end if

  ! Decompose the diffusion matrix.
  decomp = fin_vol_lu_decomp(dt, p%section([1,ncol],[ktop,kbot]), &
       coef_q_diff=egwdffi(:,ktop:kbot+1)*dpidz_sq(:,ktop:kbot+1))

end subroutine gw_ediff

!==========================================================================

subroutine gw_ediff_prep_append_proof(proof_line)

  character(len=*), intent(in) :: proof_line

  character(len=512) :: proof_file
  integer :: status, n, unitno

  proof_file = ''
  call get_environment_variable('GW_EDIFF_PREP_PROOF_FILE', value=proof_file, length=n, status=status)
  if (status == 0 .and. n > 0) then
     open(newunit=unitno, file=trim(proof_file(:n)), status='unknown', position='append', action='write')
     write(unitno,'(A)') trim(proof_line)
     close(unitno)
  end if

end subroutine gw_ediff_prep_append_proof

!==========================================================================

subroutine gw_ediff_prep_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (gw_ediff_prep_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('GW_EDIFF_PREP_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_gw_ediff_prep_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_gw_ediff_prep_impl = .false.
  end if

  gw_ediff_prep_impl_selected = .true.

  if (masterproc) then
     if (use_native_gw_ediff_prep_impl) then
        write(iulog,*) 'gw_ediff_prep implementation = native'
        call gw_ediff_prep_append_proof('gw_ediff_prep selector entered implementation = native')
     else
        write(iulog,*) 'gw_ediff implementation = codon'
        call gw_ediff_prep_append_proof('gw_ediff selector entered implementation = codon')
     end if
     call flush(iulog)
  end if

end subroutine gw_ediff_prep_select_impl

!==========================================================================

subroutine gw_ediff_prep_note_entered()

  if (gw_ediff_prep_entered_logged) return
  gw_ediff_prep_entered_logged = .true.

  if (masterproc) then
     write(iulog,*) 'gw_ediff direct = codon; diffusivity prep direct = codon; TriDiagDecomp native private-type island'
     call gw_ediff_prep_append_proof('gw_ediff direct = codon; diffusivity prep direct = codon; TriDiagDecomp native private-type island')
     call flush(iulog)
  end if

end subroutine gw_ediff_prep_note_entered

!==========================================================================

subroutine gw_ediff_prep_codon_wrap(ncol_local, pver_local, pverp_local, ngwv_local, &
     kbot_local, ktop_local, prndl_local, gravit_local, gwut_local, ubm_local, nm_local, &
     rho_local, c_local, tend_level_local, egwdffi_local, egwdffm_local, egwdff_lev_local, dpidz_sq_local, &
     ro_adjust_local)

  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_null_ptr, c_ptr

  integer, intent(in) :: ncol_local, pver_local, pverp_local, ngwv_local
  integer, intent(in) :: kbot_local, ktop_local
  real(r8), intent(in) :: prndl_local, gravit_local
  real(r8), target, intent(in) :: gwut_local(ncol_local,pver_local,-ngwv_local:ngwv_local)
  real(r8), target, intent(in) :: ubm_local(ncol_local,pver_local), nm_local(ncol_local,pver_local)
  real(r8), target, intent(in) :: rho_local(ncol_local,pverp_local)
  real(r8), target, intent(in) :: c_local(ncol_local,-ngwv_local:ngwv_local)
  integer, intent(in) :: tend_level_local(ncol_local)
  real(r8), target, intent(inout) :: egwdffi_local(ncol_local,pverp_local)
  real(r8), target, intent(inout) :: egwdffm_local(ncol_local,pver_local)
  real(r8), target, intent(inout) :: egwdff_lev_local(ncol_local)
  real(r8), target, intent(inout) :: dpidz_sq_local(ncol_local,pverp_local)
  real(r8), target, intent(in), optional :: ro_adjust_local(ncol_local,-ngwv_local:ngwv_local,pverp_local)

  integer(c_int64_t), target :: tend_level_i8(ncol_local)
  type(c_ptr) :: ro_adjust_ptr
  integer :: i

  interface
     subroutine gw_ediff_prep_codon(ncol_c, pver_c, pverp_c, ngwv_c, kbot_c, ktop_c, &
          prndl_c, gravit_c, gwut_p, ubm_p, nm_p, rho_p, c_p, tend_level_p, &
          egwdffi_p, egwdffm_p, egwdff_lev_p, dpidz_sq_p, ro_adjust_present_c, ro_adjust_p) &
          bind(c, name="gw_ediff_prep_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pver_c, pverp_c, ngwv_c, kbot_c, ktop_c
       integer(c_int64_t), value :: ro_adjust_present_c
       real(c_double), value :: prndl_c, gravit_c
       type(c_ptr), value :: gwut_p, ubm_p, nm_p, rho_p, c_p, tend_level_p
       type(c_ptr), value :: egwdffi_p, egwdffm_p, egwdff_lev_p, dpidz_sq_p, ro_adjust_p
     end subroutine gw_ediff_prep_codon
  end interface

  ro_adjust_ptr = c_null_ptr
  if (present(ro_adjust_local)) ro_adjust_ptr = c_loc(ro_adjust_local)

  do i = 1, ncol_local
     tend_level_i8(i) = int(tend_level_local(i), c_int64_t)
  end do

  call gw_ediff_prep_codon(int(ncol_local, c_int64_t), int(pver_local, c_int64_t), &
       int(pverp_local, c_int64_t), int(ngwv_local, c_int64_t), int(kbot_local, c_int64_t), &
       int(ktop_local, c_int64_t), real(prndl_local, c_double), real(gravit_local, c_double), &
       c_loc(gwut_local), c_loc(ubm_local), c_loc(nm_local), c_loc(rho_local), c_loc(c_local), &
       c_loc(tend_level_i8), c_loc(egwdffi_local), c_loc(egwdffm_local), c_loc(egwdff_lev_local), &
       c_loc(dpidz_sq_local), merge(1_c_int64_t, 0_c_int64_t, present(ro_adjust_local)), ro_adjust_ptr)

end subroutine gw_ediff_prep_codon_wrap

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
  real(r8), target :: qnew(ncol,pver)
  ! Whether to keep pre/post solve copies and tendency conversion native.
  logical :: use_native_prepost
  integer(c_int64_t) :: gw_diff_tend_touch_c

!--------------------------------------------------------------------------

  call gw_diff_tend_prepost_select_impl()
  use_native_prepost = use_native_gw_diff_tend_prepost_impl
  if (.not. use_native_prepost) gw_diff_tend_touch_c = gw_diff_tend_codon(1602_c_int64_t)

  if (use_native_prepost) then
     dq   = 0.0_r8
     qnew = q
  else
     call gw_diff_tend_prepost_note_entered()
     call gw_diff_tend_prepost_codon_wrap(1, ncol, pver, dt, q, qnew, dq)
  end if

  call decomp%left_div(qnew(:,ktop:kbot))

  ! Evaluate tendency to be reported back.
  if (use_native_prepost) then
     dq = (qnew-q) / dt
  else
     call gw_diff_tend_prepost_codon_wrap(2, ncol, pver, dt, q, qnew, dq)
  end if

end subroutine gw_diff_tend

!==========================================================================

subroutine gw_diff_tend_prepost_append_proof(proof_line)

  character(len=*), intent(in) :: proof_line

  character(len=512) :: proof_file
  integer :: status, n, unitno

  proof_file = ''
  call get_environment_variable('GW_DIFF_TEND_PREPOST_PROOF_FILE', value=proof_file, length=n, status=status)
  if (status == 0 .and. n > 0) then
     open(newunit=unitno, file=trim(proof_file(:n)), status='unknown', position='append', action='write')
     write(unitno,'(A)') trim(proof_line)
     close(unitno)
  end if

end subroutine gw_diff_tend_prepost_append_proof

!==========================================================================

subroutine gw_diff_tend_prepost_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (gw_diff_tend_prepost_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('GW_DIFF_TEND_PREPOST_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_gw_diff_tend_prepost_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_gw_diff_tend_prepost_impl = .false.
  end if

  gw_diff_tend_prepost_impl_selected = .true.

  if (masterproc) then
     if (use_native_gw_diff_tend_prepost_impl) then
        write(iulog,*) 'gw_diff_tend_prepost implementation = native'
        call gw_diff_tend_prepost_append_proof('gw_diff_tend_prepost selector entered implementation = native')
     else
        write(iulog,*) 'gw_diff_tend implementation = codon'
        call gw_diff_tend_prepost_append_proof('gw_diff_tend selector entered implementation = codon')
     end if
     call flush(iulog)
  end if

end subroutine gw_diff_tend_prepost_select_impl

!==========================================================================

subroutine gw_diff_tend_prepost_note_entered()

  if (gw_diff_tend_prepost_entered_logged) return
  gw_diff_tend_prepost_entered_logged = .true.

  if (masterproc) then
     write(iulog,*) 'gw_diff_tend direct = codon; copy/tendency direct = codon; left_div native private-type island'
     call gw_diff_tend_prepost_append_proof('gw_diff_tend direct = codon; copy/tendency direct = codon; left_div native private-type island')
     call flush(iulog)
  end if

end subroutine gw_diff_tend_prepost_note_entered

!==========================================================================

subroutine gw_diff_tend_prepost_codon_wrap(stage, ncol_local, pver_local, dt_local, q_local, qnew_local, dq_local)

  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

  integer, intent(in) :: stage
  integer, intent(in) :: ncol_local, pver_local
  real(r8), intent(in) :: dt_local
  real(r8), target, intent(in) :: q_local(ncol_local,pver_local)
  real(r8), target, intent(inout) :: qnew_local(ncol_local,pver_local)
  real(r8), target, intent(inout) :: dq_local(ncol_local,pver_local)

  interface
     subroutine gw_diff_tend_prepost_codon(stage_c, ncol_c, pver_c, dt_c, q_p, qnew_p, dq_p) &
          bind(c, name="gw_diff_tend_prepost_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: stage_c, ncol_c, pver_c
       real(c_double), value :: dt_c
       type(c_ptr), value :: q_p, qnew_p, dq_p
     end subroutine gw_diff_tend_prepost_codon
  end interface

  call gw_diff_tend_prepost_codon(int(stage, c_int64_t), int(ncol_local, c_int64_t), &
       int(pver_local, c_int64_t), real(dt_local, c_double), c_loc(q_local), c_loc(qnew_local), c_loc(dq_local))

end subroutine gw_diff_tend_prepost_codon_wrap

end module gw_diffusion
