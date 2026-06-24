module gw_utils

!
! This module contains utility code for the gravity wave modules.
!

use iso_c_binding, only: c_int64_t, c_loc

implicit none
private
save

! Real kind for gravity wave parameterization.
integer, public, parameter :: r8 = selected_real_kind(12)

! Public interface
public :: get_unit_vector
public :: dot_2d
public :: midpoint_interp

logical :: use_native_gw_utils_impl = .false.
logical :: gw_utils_impl_selected = .false.
logical :: gw_utils_proof_written = .false.
logical :: gw_utils_get_unit_vector_logged = .false.
logical :: gw_utils_dot_2d_logged = .false.
logical :: gw_utils_midpoint_interp_logged = .false.

interface
   subroutine get_unit_vector_codon(u_p, v_p, u_n_p, v_n_p, mag_p, n_c) &
        bind(c, name="get_unit_vector_codon")
     use iso_c_binding, only: c_int64_t, c_ptr
     type(c_ptr), value :: u_p, v_p, u_n_p, v_n_p, mag_p
     integer(c_int64_t), value :: n_c
   end subroutine get_unit_vector_codon

   subroutine dot_2d_codon(u1_p, v1_p, u2_p, v2_p, out_p, n_c) &
        bind(c, name="dot_2d_codon")
     use iso_c_binding, only: c_int64_t, c_ptr
     type(c_ptr), value :: u1_p, v1_p, u2_p, v2_p, out_p
     integer(c_int64_t), value :: n_c
   end subroutine dot_2d_codon

   subroutine midpoint_interp_codon(arr_p, interp_p, n1_c, n2_c) &
        bind(c, name="midpoint_interp_codon")
     use iso_c_binding, only: c_int64_t, c_ptr
     type(c_ptr), value :: arr_p, interp_p
     integer(c_int64_t), value :: n1_c, n2_c
   end subroutine midpoint_interp_codon
end interface

contains

subroutine gw_utils_select_impl()

  use cam_logfile, only: iulog
  use spmd_utils, only: masterproc

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (gw_utils_impl_selected) return

  impl_name = 'codon'
  call cam_codon_get_impl('GW_UTILS_IMPL', impl_name, n, status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_gw_utils_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_gw_utils_impl = .false.
  end if

  gw_utils_impl_selected = .true.

  if (masterproc) then
     if (use_native_gw_utils_impl) then
        write(iulog,*) 'gw_utils implementation = native'
     else
        write(iulog,*) 'gw_utils implementation = codon'
     end if
  end if

end subroutine gw_utils_select_impl

subroutine gw_utils_proof_once()

  use cam_logfile, only: iulog
  use spmd_utils, only: masterproc

  if (gw_utils_proof_written) return
  gw_utils_proof_written = .true.

  if (masterproc) then
     write(iulog,'(A)') 'gw_utils entered (vector and midpoint helpers = codon)'
  end if

end subroutine gw_utils_proof_once

subroutine gw_utils_log_direct(logged, proof_line)

  use cam_logfile, only: iulog
  use spmd_utils, only: masterproc

  logical, intent(inout) :: logged
  character(len=*), intent(in) :: proof_line

  if (logged) return
  logged = .true.

  if (masterproc) then
     write(iulog,'(A)') trim(proof_line)
  end if

end subroutine gw_utils_log_direct

! Take two components of a vector, and find the unit vector components and
! total magnitude.
subroutine get_unit_vector(u, v, u_n, v_n, mag)
  real(r8), contiguous, target, intent(in) :: u(:)
  real(r8), contiguous, target, intent(in) :: v(:)
  real(r8), contiguous, target, intent(out) :: u_n(:)
  real(r8), contiguous, target, intent(out) :: v_n(:)
  real(r8), contiguous, target, intent(out) :: mag(:)

  integer :: i

  call gw_utils_select_impl()

  if (.not. use_native_gw_utils_impl) then
     call gw_utils_proof_once()
     call get_unit_vector_codon(c_loc(u), c_loc(v), c_loc(u_n), c_loc(v_n), c_loc(mag), &
          int(size(mag), c_int64_t))
     call gw_utils_log_direct(gw_utils_get_unit_vector_logged, 'get_unit_vector direct = codon')
     return
  end if

  mag = sqrt(u*u + v*v)

  ! Has to be a loop/if instead of a where, because floating point
  ! exceptions can trigger even on a masked divide-by-zero operation
  ! (especially on Intel).
  do i = 1, size(mag)
     if (mag(i) > 0._r8) then
        u_n(i) = u(i)/mag(i)
        v_n(i) = v(i)/mag(i)
     else
        u_n(i) = 0._r8
        v_n(i) = 0._r8
     end if
  end do

end subroutine get_unit_vector

! Vectorized version of a 2D dot product (since the intrinsic dot_product
! is more suitable for arrays of contiguous vectors).
function dot_2d(u1, v1, u2, v2)
  real(r8), contiguous, target, intent(in) :: u1(:), v1(:)
  real(r8), contiguous, target, intent(in) :: u2(:), v2(:)

  real(r8), target :: dot_2d(size(u1))

  call gw_utils_select_impl()

  if (.not. use_native_gw_utils_impl) then
     call gw_utils_proof_once()
     call dot_2d_codon(c_loc(u1), c_loc(v1), c_loc(u2), c_loc(v2), c_loc(dot_2d), &
          int(size(u1), c_int64_t))
     call gw_utils_log_direct(gw_utils_dot_2d_logged, 'dot_2d direct = codon')
  else
     dot_2d = u1*u2 + v1*v2
  end if

end function dot_2d

! Pure function that interpolates the values of the input array along
! dimension 2. This is obviously not a very generic routine, unlike, say,
! CAM's lininterp. But it's used often enough that it seems worth providing
! here.
function midpoint_interp(arr) result(interp)
  real(r8), contiguous, target, intent(in) :: arr(:,:)
  real(r8), target :: interp(size(arr,1),size(arr,2)-1)

  integer :: i

  call gw_utils_select_impl()

  if (.not. use_native_gw_utils_impl) then
     call gw_utils_proof_once()
     call midpoint_interp_codon(c_loc(arr), c_loc(interp), &
          int(size(arr,1), c_int64_t), int(size(arr,2), c_int64_t))
     call gw_utils_log_direct(gw_utils_midpoint_interp_logged, 'midpoint_interp direct = codon')
  else
     do i = 1, size(interp,2)
        interp(:,i) = 0.5_r8 * (arr(:,i)+arr(:,i+1))
     end do
  end if

end function midpoint_interp

end module gw_utils
