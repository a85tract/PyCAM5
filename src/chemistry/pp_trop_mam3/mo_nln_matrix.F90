




      module mo_nln_matrix

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: nlnmat

      logical :: nlnmat_use_native_impl = .false.
      logical :: nlnmat_impl_selected = .false.

      contains

      subroutine nlnmat( mat, y, rxt, lmat, dti )

      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use iso_c_binding, only : c_double, c_int64_t, c_loc, c_ptr

      implicit none

!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: dti
      real(r8), target, intent(in) :: lmat(nzcnt)
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), target, intent(inout) :: mat(nzcnt)

      interface
         subroutine nlnmat_codon(mat_p, lmat_p, dti_c) bind(c, name="nlnmat_codon")
            use iso_c_binding, only : c_double, c_ptr
            type(c_ptr), value :: mat_p, lmat_p
            real(c_double), value :: dti_c
         end subroutine nlnmat_codon
      end interface

      call nlnmat_select_impl()

      if (.not. nlnmat_use_native_impl) then
         call nlnmat_codon(c_loc(mat), c_loc(lmat), real(dti, c_double))
         return
      end if

      call nlnmat_finit( mat, lmat, dti )
      end subroutine nlnmat
      subroutine nlnmat_finit( mat, lmat, dti )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: dti
      real(r8), intent(in) :: lmat(nzcnt)
      real(r8), intent(inout) :: mat(nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
         mat( 1) = lmat( 1)
         mat( 2) = lmat( 2)
         mat( 3) = lmat( 3)
         mat( 4) = lmat( 4)
         mat( 5) = lmat( 5)
         mat( 6) = lmat( 6)
         mat( 7) = lmat( 7)
         mat( 8) = lmat( 8)
         mat( 9) = lmat( 9)
         mat( 10) = lmat( 10)
         mat( 11) = lmat( 11)
         mat( 12) = lmat( 12)
         mat( 13) = lmat( 13)
         mat( 14) = lmat( 14)
         mat( 15) = lmat( 15)
         mat( 16) = lmat( 16)
         mat( 17) = lmat( 17)
         mat( 18) = lmat( 18)
         mat( 19) = lmat( 19)
         mat( 20) = lmat( 20)
         mat( 21) = lmat( 21)
         mat( 22) = lmat( 22)
         mat( 1) = mat( 1) - dti
         mat( 2) = mat( 2) - dti
         mat( 4) = mat( 4) - dti
         mat( 6) = mat( 6) - dti
         mat( 7) = mat( 7) - dti
         mat( 8) = mat( 8) - dti
         mat( 9) = mat( 9) - dti
         mat( 10) = mat( 10) - dti
         mat( 11) = mat( 11) - dti
         mat( 12) = mat( 12) - dti
         mat( 13) = mat( 13) - dti
         mat( 14) = mat( 14) - dti
         mat( 15) = mat( 15) - dti
         mat( 16) = mat( 16) - dti
         mat( 17) = mat( 17) - dti
         mat( 18) = mat( 18) - dti
         mat( 19) = mat( 19) - dti
         mat( 20) = mat( 20) - dti
         mat( 21) = mat( 21) - dti
         mat( 22) = mat( 22) - dti
      end subroutine nlnmat_finit

      subroutine nlnmat_select_impl()

      use cam_logfile, only : iulog
      use spmd_utils, only : masterproc

      implicit none

      character(len=32) :: impl_name
      integer :: status, n, i, code

      if (nlnmat_impl_selected) return

      impl_name = 'codon'
      call get_environment_variable('NLNMAT_IMPL', value=impl_name, length=n, status=status)

      if (status == 0 .and. n > 0) then
         do i = 1, n
            code = iachar(impl_name(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
               impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
         end do
         nlnmat_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
      else
         nlnmat_use_native_impl = .false.
      end if

      nlnmat_impl_selected = .true.

      if (masterproc) then
         if (nlnmat_use_native_impl) then
            write(iulog,*) 'nlnmat implementation = native'
         else
            write(iulog,*) 'nlnmat implementation = codon'
         end if
         call flush(iulog)
      end if

      end subroutine nlnmat_select_impl

      end module mo_nln_matrix
