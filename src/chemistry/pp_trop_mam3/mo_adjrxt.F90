




      module mo_adjrxt

      private
      public :: adjrxt

      logical :: adjrxt_use_native_impl = .false.
      logical :: adjrxt_impl_selected = .false.

      contains

      subroutine adjrxt( rate, inv, m, ncol )

      use ppgrid, only : pver
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : nfs, rxntot
      use iso_c_binding, only : c_int64_t, c_loc, c_ptr

      implicit none

!--------------------------------------------------------------------
! ... dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: ncol
      real(r8), target, intent(in) :: inv(ncol,pver,nfs)
      real(r8), target, intent(in) :: m(ncol,pver)
      real(r8), target, intent(inout) :: rate(ncol,pver,rxntot)

!--------------------------------------------------------------------
! ... local variables
!--------------------------------------------------------------------
      real(r8) :: im(ncol,pver)

      interface
         subroutine adjrxt_codon(ncol_c, pver_c, rate_p, inv_p, m_p) bind(c, name="adjrxt_codon")
            use iso_c_binding, only : c_int64_t, c_ptr
            integer(c_int64_t), value :: ncol_c, pver_c
            type(c_ptr), value :: rate_p, inv_p, m_p
         end subroutine adjrxt_codon
      end interface

      call adjrxt_select_impl()

      if (.not. adjrxt_use_native_impl) then
         call adjrxt_codon( &
              int(ncol, c_int64_t), int(pver, c_int64_t), c_loc(rate), c_loc(inv), c_loc(m) &
         )
         return
      end if


         rate(:,:, 3) = rate(:,:, 3) * inv(:,:, 6)
         rate(:,:, 4) = rate(:,:, 4) * inv(:,:, 6)
         rate(:,:, 5) = rate(:,:, 5) * inv(:,:, 6)
         rate(:,:, 6) = rate(:,:, 6) * inv(:,:, 6)
         rate(:,:, 7) = rate(:,:, 7) * inv(:,:, 7)
         im(:,:) = 1._r8 / m(:,:)
         rate(:,:, 2) = rate(:,:, 2) * inv(:,:, 8) * inv(:,:, 8) * im(:,:)

      end subroutine adjrxt

      subroutine adjrxt_select_impl()

      use cam_logfile, only : iulog
      use spmd_utils, only : masterproc

      implicit none

      character(len=32) :: impl_name
      integer :: status, n, i, code

      if (adjrxt_impl_selected) return

      impl_name = 'codon'
      call get_environment_variable('ADJRXT_IMPL', value=impl_name, length=n, status=status)

      if (status == 0 .and. n > 0) then
         do i = 1, n
            code = iachar(impl_name(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
               impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
         end do
         adjrxt_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
      else
         adjrxt_use_native_impl = .false.
      end if

      adjrxt_impl_selected = .true.

      if (masterproc) then
         if (adjrxt_use_native_impl) then
            write(iulog,*) 'adjrxt implementation = native'
         else
            write(iulog,*) 'adjrxt implementation = codon'
         end if
         call flush(iulog)
      end if

      end subroutine adjrxt_select_impl

      end module mo_adjrxt
