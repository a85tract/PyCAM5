




      module mo_lu_solve

      private
      public :: lu_slv

      logical :: lu_slv_use_native_impl = .false.
      logical :: lu_slv_impl_selected = .false.

      contains

      subroutine lu_slv01( lu, b )


      use shr_kind_mod, only : r8 => shr_kind_r8
      use iso_c_binding, only : c_loc, c_ptr

      implicit none

!-----------------------------------------------------------------------
! ... Dummy args
!-----------------------------------------------------------------------
      real(r8), target, intent(in) :: lu(:)
      real(r8), target, intent(inout) :: b(:)

!-----------------------------------------------------------------------
! ... Local variables
!-----------------------------------------------------------------------

      interface
         subroutine lu_slv_codon(lu_p, b_p) bind(c, name="lu_slv_codon")
            use iso_c_binding, only : c_ptr
            type(c_ptr), value :: lu_p, b_p
         end subroutine lu_slv_codon
      end interface

      call lu_slv_select_impl()

      if (.not. lu_slv_use_native_impl) then
         call lu_slv_codon(c_loc(lu), c_loc(b))
         return
      end if

!-----------------------------------------------------------------------
! ... solve L * y = b
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! ... Solve U * x = y
!-----------------------------------------------------------------------
         b(20) = b(20) * lu(22)
         b(19) = b(19) * lu(21)
         b(18) = b(18) * lu(20)
         b(17) = b(17) * lu(19)
         b(16) = b(16) * lu(18)
         b(15) = b(15) * lu(17)
         b(14) = b(14) * lu(16)
         b(13) = b(13) * lu(15)
         b(12) = b(12) * lu(14)
         b(11) = b(11) * lu(13)
         b(10) = b(10) * lu(12)
         b(9) = b(9) * lu(11)
         b(8) = b(8) * lu(10)
         b(7) = b(7) * lu(9)
         b(6) = b(6) * lu(8)
         b(5) = b(5) * lu(7)
         b(4) = b(4) * lu(6)
         b(3) = b(3) - lu(5) * b(4)
         b(3) = b(3) * lu(4)
         b(2) = b(2) - lu(3) * b(3)
         b(2) = b(2) * lu(2)
         b(1) = b(1) * lu(1)
      end subroutine lu_slv01
      subroutine lu_slv( lu, b )
      use shr_kind_mod, only : r8 => shr_kind_r8
      implicit none
!-----------------------------------------------------------------------
! ... Dummy args
!-----------------------------------------------------------------------
      real(r8), intent(in) :: lu(:)
      real(r8), intent(inout) :: b(:)
      call lu_slv01( lu, b )
      end subroutine lu_slv

      subroutine lu_slv_select_impl()

      use cam_logfile, only : iulog
      use spmd_utils, only : masterproc

      implicit none

      character(len=32) :: impl_name
      integer :: status, n, i, code

      if (lu_slv_impl_selected) return

      impl_name = 'codon'
      call cam_codon_get_impl('LU_SLV_IMPL', impl_name, n, status)

      if (status == 0 .and. n > 0) then
         do i = 1, n
            code = iachar(impl_name(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
               impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
         end do
         lu_slv_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
      else
         lu_slv_use_native_impl = .false.
      end if

      lu_slv_impl_selected = .true.

      if (masterproc) then
         if (lu_slv_use_native_impl) then
            write(iulog,*) 'lu_slv implementation = native'
         else
            write(iulog,*) 'lu_slv implementation = codon'
         end if
         call flush(iulog)
      end if

      end subroutine lu_slv_select_impl

      end module mo_lu_solve
