




      module mo_lu_factor

      private
      public :: lu_fac

      logical :: lu_fac_use_native_impl = .false.
      logical :: lu_fac_impl_selected = .false.
      logical :: lu_fac01_logged = .false.

      contains

      subroutine lu_fac01( lu )


      use shr_kind_mod, only : r8 => shr_kind_r8
      use iso_c_binding, only : c_loc, c_ptr
      use cam_logfile, only : iulog
      use spmd_utils, only : masterproc

      implicit none

!-----------------------------------------------------------------------
! ... dummy args
!-----------------------------------------------------------------------
      real(r8), target, intent(inout) :: lu(:)

      interface
         subroutine lu_fac01_codon(lu_p) bind(c, name="lu_fac01_codon")
            use iso_c_binding, only : c_ptr
            type(c_ptr), value :: lu_p
         end subroutine lu_fac01_codon
      end interface

      call lu_fac_select_impl()

      if (.not. lu_fac_use_native_impl) then
         call lu_fac01_codon(c_loc(lu))
         if (masterproc .and. .not. lu_fac01_logged) then
            write(iulog,*) 'lu_fac01 implementation = codon'
            lu_fac01_logged = .true.
            call flush(iulog)
         end if
         return
      end if

         lu(1) = 1._r8 / lu(1)

         lu(2) = 1._r8 / lu(2)

         lu(4) = 1._r8 / lu(4)

         lu(6) = 1._r8 / lu(6)

         lu(7) = 1._r8 / lu(7)

         lu(8) = 1._r8 / lu(8)

         lu(9) = 1._r8 / lu(9)

         lu(10) = 1._r8 / lu(10)

         lu(11) = 1._r8 / lu(11)

         lu(12) = 1._r8 / lu(12)

         lu(13) = 1._r8 / lu(13)

         lu(14) = 1._r8 / lu(14)

         lu(15) = 1._r8 / lu(15)

         lu(16) = 1._r8 / lu(16)

         lu(17) = 1._r8 / lu(17)

         lu(18) = 1._r8 / lu(18)

         lu(19) = 1._r8 / lu(19)

         lu(20) = 1._r8 / lu(20)

         lu(21) = 1._r8 / lu(21)

         lu(22) = 1._r8 / lu(22)


      end subroutine lu_fac01

      subroutine lu_fac( lu )


      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!-----------------------------------------------------------------------
! ... dummy args
!-----------------------------------------------------------------------
      real(r8), intent(inout) :: lu(:)

      call lu_fac01( lu )

      end subroutine lu_fac

      subroutine lu_fac_select_impl()

      use cam_logfile, only : iulog
      use spmd_utils, only : masterproc

      implicit none

      character(len=32) :: impl_name
      integer :: status, n, i, code

      if (lu_fac_impl_selected) return

      impl_name = 'codon'
      call cam_codon_get_impl('LU_FAC_IMPL', impl_name, n, status)

      if (status == 0 .and. n > 0) then
         do i = 1, n
            code = iachar(impl_name(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
               impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
         end do
         lu_fac_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
      else
         lu_fac_use_native_impl = .false.
      end if

      lu_fac_impl_selected = .true.

      if (masterproc) then
         if (lu_fac_use_native_impl) then
            write(iulog,*) 'lu_fac implementation = native'
         else
            write(iulog,*) 'lu_fac implementation = codon'
         end if
         call flush(iulog)
      end if

      end subroutine lu_fac_select_impl

      end module mo_lu_factor
