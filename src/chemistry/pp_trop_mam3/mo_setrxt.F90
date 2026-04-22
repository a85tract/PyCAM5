
      module mo_setrxt

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: setrxt
      public :: setrxt_hrates

      logical :: setrxt_use_native_impl = .false.
      logical :: setrxt_impl_selected = .false.

      contains

      subroutine setrxt( rate, temp, m, ncol )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl
      use iso_c_binding, only : c_int64_t, c_loc, c_ptr

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      real(r8), target, intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), target, intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer  ::  n
      real(r8)  ::  itemp(ncol,pver)
      real(r8)  ::  exp_fac(ncol,pver)

      interface
         subroutine setrxt_codon(ncol_c, pcols_c, pver_c, temp_p, rate_p) bind(c, name="setrxt_codon")
            use iso_c_binding, only : c_int64_t, c_ptr
            integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
            type(c_ptr), value :: temp_p, rate_p
         end subroutine setrxt_codon
      end interface

      call setrxt_select_impl()

      if (.not. setrxt_use_native_impl) then
         call setrxt_codon( &
              int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), c_loc(temp), c_loc(rate) &
         )
         return
      end if

      itemp(:ncol,:) = 1._r8 / temp(:ncol,:)
      n = ncol*pver
      rate(:,:,3) = 2.9e-12_r8 * exp( -160._r8 * itemp(:,:) )
      rate(:,:,5) = 9.6e-12_r8 * exp( -234._r8 * itemp(:,:) )
      rate(:,:,7) = 1.9e-13_r8 * exp( 520._r8 * itemp(:,:) )

      end subroutine setrxt

      subroutine setrxt_select_impl()

      use cam_logfile, only : iulog
      use spmd_utils, only : masterproc

      implicit none

      character(len=32) :: impl_name
      integer :: status, n, i, code

      if (setrxt_impl_selected) return

      impl_name = 'codon'
      call get_environment_variable('SETRXT_IMPL', value=impl_name, length=n, status=status)

      if (status == 0 .and. n > 0) then
         do i = 1, n
            code = iachar(impl_name(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
               impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
         end do
         setrxt_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
      else
         setrxt_use_native_impl = .false.
      end if

      setrxt_impl_selected = .true.

      if (masterproc) then
         if (setrxt_use_native_impl) then
            write(iulog,*) 'setrxt implementation = native'
         else
            write(iulog,*) 'setrxt implementation = codon'
         end if
         call flush(iulog)
      end if

      end subroutine setrxt_select_impl


      subroutine setrxt_hrates( rate, temp, m, ncol, kbot )

      use ppgrid,       only : pver, pcols
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : rxntot
      use mo_jpl,    only : jpl

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      integer, intent(in) :: ncol
      integer, intent(in) :: kbot
      real(r8), intent(in)    :: temp(pcols,pver)
      real(r8), intent(in)    :: m(ncol,pver)
      real(r8), intent(inout) :: rate(ncol,pver,rxntot)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      integer  ::  n
      real(r8)  ::  itemp(ncol,kbot)
      real(r8)  ::  exp_fac(ncol,kbot)


      end subroutine setrxt_hrates

      end module mo_setrxt
