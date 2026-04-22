
      module mo_negtrc

      private
      public :: negtrc

      logical :: negtrc_use_native_impl = .false.
      logical :: negtrc_impl_selected = .false.

      contains

      subroutine negtrc( header, fld, ncol )
!-----------------------------------------------------------------------
!  	... Check for negative constituent values and
!	    replace with zero value
!-----------------------------------------------------------------------

      use shr_kind_mod, only: r8 => shr_kind_r8
      use chem_mods,   only : gas_pcnst
      use iso_c_binding, only : c_int64_t, c_loc, c_ptr
      use ppgrid,      only : pver

      implicit none

!-----------------------------------------------------------------------
!  	... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in)          :: ncol
      character(len=*), intent(in) :: header
      real(r8), target, intent(inout) :: fld(ncol,pver,gas_pcnst) ! field to check

      interface
         subroutine negtrc_codon(ncol_c, pver_c, gas_pcnst_c, fld_p) bind(c, name="negtrc_codon")
            use iso_c_binding, only : c_int64_t, c_ptr
            integer(c_int64_t), value :: ncol_c, pver_c, gas_pcnst_c
            type(c_ptr), value :: fld_p
         end subroutine negtrc_codon
      end interface

      call negtrc_select_impl()

      if (.not. negtrc_use_native_impl) then
         call negtrc_codon(int(ncol, c_int64_t), int(pver, c_int64_t), int(gas_pcnst, c_int64_t), c_loc(fld))
         return
      end if

      call negtrc_native( header, fld, ncol )
      end subroutine negtrc

      subroutine negtrc_native( header, fld, ncol )
!-----------------------------------------------------------------------
!  	... Check for negative constituent values and
!	    replace with zero value
!-----------------------------------------------------------------------

      use shr_kind_mod, only: r8 => shr_kind_r8
      use chem_mods,   only : gas_pcnst
      use ppgrid,      only : pver

      implicit none

!-----------------------------------------------------------------------
!  	... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in)          :: ncol
      character(len=*), intent(in) :: header
      real(r8), intent(inout)      :: fld(ncol,pver,gas_pcnst) ! field to check

!-----------------------------------------------------------------------
!  	... Local variables
!-----------------------------------------------------------------------
      integer :: m
      integer :: nneg                       ! flag counter

      do m  = 1,gas_pcnst
         nneg = count( fld(:,:,m) < 0._r8 )
	 if( nneg > 0 ) then
            where( fld(:,:,m) < 0._r8 )
	       fld(:,:,m) = 0._r8
	    endwhere
	 end if
      end do

      end subroutine negtrc_native

      subroutine negtrc_select_impl()

      use cam_logfile, only : iulog
      use spmd_utils, only : masterproc

      implicit none

      character(len=32) :: impl_name
      integer :: status, n, i, code

      if (negtrc_impl_selected) return

      impl_name = 'codon'
      call get_environment_variable('NEGTRC_IMPL', value=impl_name, length=n, status=status)

      if (status == 0 .and. n > 0) then
         do i = 1, n
            code = iachar(impl_name(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
               impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
         end do
         negtrc_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
      else
         negtrc_use_native_impl = .false.
      end if

      negtrc_impl_selected = .true.

      if (masterproc) then
         if (negtrc_use_native_impl) then
            write(iulog,*) 'negtrc implementation = native'
         else
            write(iulog,*) 'negtrc implementation = codon'
         end if
         call flush(iulog)
      end if

      end subroutine negtrc_select_impl

      end module mo_negtrc
