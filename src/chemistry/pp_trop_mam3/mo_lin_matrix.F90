




      module mo_lin_matrix

      private
      public :: linmat

      logical :: linmat_use_native_impl = .false.
      logical :: linmat_impl_selected = .false.

      contains

      subroutine linmat01( mat, y, rxt, het_rates )
!----------------------------------------------
! ... linear matrix entries for implicit species
!----------------------------------------------

      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8
      use iso_c_binding, only : c_loc, c_ptr

      implicit none

!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), target, intent(in) :: rxt(rxntot)
      real(r8), target, intent(in) :: het_rates(max(1,gas_pcnst))
      real(r8), target, intent(inout) :: mat(nzcnt)

      interface
         subroutine linmat_codon(mat_p, rxt_p, het_rates_p) bind(c, name="linmat_codon")
            use iso_c_binding, only : c_ptr
            type(c_ptr), value :: mat_p, rxt_p, het_rates_p
         end subroutine linmat_codon
      end interface

      call linmat_select_impl()

      if (.not. linmat_use_native_impl) then
         call linmat_codon(c_loc(mat), c_loc(rxt), c_loc(het_rates))
         return
      end if

         mat(1) = -( rxt(1) + rxt(3) + het_rates(1) )

         mat(2) = -( het_rates(2) )
         mat(3) = rxt(4)

         mat(4) = -( rxt(4) + het_rates(3) )
         mat(5) = rxt(5) + .500_r8*rxt(6) + rxt(7)

         mat(6) = -( rxt(5) + rxt(6) + rxt(7) + het_rates(4) )

         mat(7) = -( het_rates(5) )

         mat(8) = -( het_rates(6) )

         mat(9) = -( het_rates(7) )

         mat(10) = -( het_rates(8) )

         mat(11) = -( het_rates(9) )

         mat(12) = -( het_rates(10) )

         mat(13) = -( het_rates(11) )

         mat(14) = -( het_rates(12) )

         mat(15) = -( het_rates(13) )

         mat(16) = -( het_rates(14) )

         mat(17) = -( het_rates(15) )

         mat(18) = -( het_rates(16) )

         mat(19) = -( het_rates(17) )

         mat(20) = -( het_rates(18) )

         mat(21) = -( het_rates(19) )

         mat(22) = -( het_rates(20) )


      end subroutine linmat01

      subroutine linmat( mat, y, rxt, het_rates )
!----------------------------------------------
! ... linear matrix entries for implicit species
!----------------------------------------------

      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      use shr_kind_mod, only : r8 => shr_kind_r8

      implicit none

!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(in) :: het_rates(max(1,gas_pcnst))
      real(r8), intent(inout) :: mat(nzcnt)

      call linmat01( mat, y, rxt, het_rates )

      end subroutine linmat

      subroutine linmat_select_impl()

      use cam_logfile, only : iulog
      use spmd_utils, only : masterproc

      implicit none

      character(len=32) :: impl_name
      integer :: status, n, i, code

      if (linmat_impl_selected) return

      impl_name = 'codon'
      call get_environment_variable('LINMAT_IMPL', value=impl_name, length=n, status=status)

      if (status == 0 .and. n > 0) then
         do i = 1, n
            code = iachar(impl_name(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
               impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
         end do
         linmat_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
      else
         linmat_use_native_impl = .false.
      end if

      linmat_impl_selected = .true.

      if (masterproc) then
         if (linmat_use_native_impl) then
            write(iulog,*) 'linmat implementation = native'
         else
            write(iulog,*) 'linmat implementation = codon'
         end if
         call flush(iulog)
      end if

      end subroutine linmat_select_impl

      end module mo_lin_matrix
