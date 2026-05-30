




      module mo_phtadj

      use cam_logfile, only : iulog
      use spmd_utils, only : masterproc

      private
      public :: phtadj

      logical, save :: phtadj_codon_logged = .false.

      contains

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

      if (phtadj_codon() /= 1_c_int64_t) then
         stop 2
      end if
      if (masterproc .and. .not. phtadj_codon_logged) then
         write(iulog,*) 'phtadj implementation = codon'
         phtadj_codon_logged = .true.
         call flush(iulog)
      end if

      end subroutine phtadj

      end module mo_phtadj
