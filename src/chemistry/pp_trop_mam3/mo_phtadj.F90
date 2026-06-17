




      module mo_phtadj

      use cam_logfile, only : iulog
      use spmd_utils, only : masterproc

      private
      public :: phtadj

      logical, save :: phtadj_codon_logged = .false.

      contains

      logical function phtadj_use_native()
      character(len=32) :: impl_name
      integer :: n, status, i, code

      impl_name = 'codon'
      call cam_codon_get_impl('PHTADJ_IMPL', impl_name, n, status)
      if (status == 0 .and. n > 0) then
         do i = 1, n
            code = iachar(impl_name(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
               impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
         end do
         phtadj_use_native = trim(adjustl(impl_name(:n))) == 'native'
      else
         phtadj_use_native = .false.
      end if
      end function phtadj_use_native

      subroutine phtadj( p_rate, inv, m, ncol )

      use chem_mods, only : nfs, phtcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
      use ppgrid, only : pver
      use iso_c_binding, only : c_int64_t

      implicit none

!--------------------------------------------------------------------
! ... dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: ncol
      real(r8), intent(in) :: inv(:,:,:)
      real(r8), intent(in) :: m(:,:)
      real(r8), intent(inout) :: p_rate(:,:,:)

!--------------------------------------------------------------------
! ... local variables
!--------------------------------------------------------------------
      interface
         function phtadj_codon() result(ok) bind(c, name="phtadj_codon")
           use iso_c_binding, only : c_int64_t
           integer(c_int64_t) :: ok
         end function phtadj_codon
      end interface

      if (masterproc .and. .not. phtadj_codon_logged) then
         if (phtadj_use_native()) then
            write(iulog,*) 'phtadj implementation = native'
         else
            write(iulog,*) 'phtadj implementation = codon'
         end if
         phtadj_codon_logged = .true.
         call flush(iulog)
      end if
      if (phtadj_use_native()) return
      if (phtadj_codon() /= 1_c_int64_t) then
         stop 2
      end if

      end subroutine phtadj

      end module mo_phtadj
