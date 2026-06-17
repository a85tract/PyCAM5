




      module mo_indprd

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: indprd

      logical :: indprd_use_native_impl = .false.
      logical :: indprd_impl_selected = .false.

      contains

      subroutine indprd( class, prod, nprod, y, extfrc, rxt, ncol )

      use chem_mods, only : gas_pcnst, extcnt, rxntot
      use iso_c_binding, only : c_int64_t, c_loc, c_ptr
      use ppgrid, only : pver

      implicit none

!--------------------------------------------------------------------
! ... dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: class
      integer, intent(in) :: ncol
      integer, intent(in) :: nprod
      real(r8), intent(in) :: y(ncol,pver,gas_pcnst)
      real(r8), target, intent(in) :: rxt(ncol,pver,rxntot)
      real(r8), target, intent(in) :: extfrc(ncol,pver,extcnt)
      real(r8), target, intent(inout) :: prod(ncol,pver,nprod)

      interface
         subroutine indprd_codon(class_c, ncol_c, pver_c, nprod_c, rxt_p, extfrc_p, prod_p) bind(c, name="indprd_codon")
            use iso_c_binding, only : c_int64_t, c_ptr
            integer(c_int64_t), value :: class_c, ncol_c, pver_c, nprod_c
            type(c_ptr), value :: rxt_p, extfrc_p, prod_p
         end subroutine indprd_codon
      end interface

      call indprd_select_impl()

      if (.not. indprd_use_native_impl) then
         call indprd_codon(int(class, c_int64_t), int(ncol, c_int64_t), int(pver, c_int64_t), &
              int(nprod, c_int64_t), c_loc(rxt), c_loc(extfrc), c_loc(prod))
         return
      end if

      call indprd_native( class, prod, nprod, y, extfrc, rxt, ncol )
      end subroutine indprd

      subroutine indprd_native( class, prod, nprod, y, extfrc, rxt, ncol )

      use chem_mods, only : gas_pcnst, extcnt, rxntot
      use ppgrid, only : pver

      implicit none

!--------------------------------------------------------------------
! ... dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: class
      integer, intent(in) :: ncol
      integer, intent(in) :: nprod
      real(r8), intent(in) :: y(ncol,pver,gas_pcnst)
      real(r8), intent(in) :: rxt(ncol,pver,rxntot)
      real(r8), intent(in) :: extfrc(ncol,pver,extcnt)
      real(r8), intent(inout) :: prod(ncol,pver,nprod)

!--------------------------------------------------------------------
! ... "independent" production for Implicit species
!--------------------------------------------------------------------
      if( class == 4 ) then
         prod(:,:,1) =rxt(:,:,2)

         prod(:,:,2) = 0._r8

         prod(:,:,3) = + extfrc(:,:,1)

         prod(:,:,4) = 0._r8

         prod(:,:,5) = 0._r8

         prod(:,:,6) = + extfrc(:,:,2)

         prod(:,:,7) = + extfrc(:,:,4)

         prod(:,:,8) = 0._r8

         prod(:,:,9) = + extfrc(:,:,5)

         prod(:,:,10) = 0._r8

         prod(:,:,11) = 0._r8

         prod(:,:,12) = + extfrc(:,:,6)

         prod(:,:,13) = + extfrc(:,:,3)

         prod(:,:,14) = 0._r8

         prod(:,:,15) = 0._r8

         prod(:,:,16) = + extfrc(:,:,7)

         prod(:,:,17) = 0._r8

         prod(:,:,18) = 0._r8

         prod(:,:,19) = 0._r8

         prod(:,:,20) = 0._r8

      end if

      end subroutine indprd_native

      subroutine indprd_select_impl()

      use cam_logfile, only : iulog
      use spmd_utils, only : masterproc

      implicit none

      character(len=32) :: impl_name
      integer :: status, n, i, code

      if (indprd_impl_selected) return

      impl_name = 'codon'
      call cam_codon_get_impl('INDPRD_IMPL', impl_name, n, status)

      if (status == 0 .and. n > 0) then
         do i = 1, n
            code = iachar(impl_name(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
               impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
         end do
         indprd_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
      else
         indprd_use_native_impl = .false.
      end if

      indprd_impl_selected = .true.

      if (masterproc) then
         if (indprd_use_native_impl) then
            write(iulog,*) 'indprd implementation = native'
         else
            write(iulog,*) 'indprd implementation = codon'
         end if
         call flush(iulog)
      end if

      end subroutine indprd_select_impl

      end module mo_indprd
