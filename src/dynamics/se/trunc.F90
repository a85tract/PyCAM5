!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  trunc --- Check consistency of truncation parameters
!
! !INTERFACE:
subroutine trunc

! !USES:
   use iso_c_binding, only : c_int64_t
   use cam_logfile, only : iulog
   implicit none

   interface
      subroutine trunc_codon() bind(c, name='trunc_codon')
      end subroutine trunc_codon
   end interface
   character(len=32) :: impl_name
   integer :: impl_n, impl_status
   logical, save :: proof_seen = .false.

   impl_name = 'codon'
   call cam_codon_get_impl('TRUNC_IMPL', impl_name, impl_n, impl_status)
   if (.not. (impl_status == 0 .and. impl_n > 0 .and. trim(adjustl(impl_name(:impl_n))) == 'native')) then
      call trunc_codon()
      if (.not. proof_seen) then
         write(iulog,*) 'trunc implementation = codon'
         proof_seen = .true.
      endif
      return
   endif

!
! !DESCRIPTION:
!
!   Check consistency of truncation parameters and evaluate pointers
!   and displacements for spectral arrays.  Note: this may not be
!   necessary for Lin-Rood. 
! 
! !REVISION HISTORY: 
!   92.06.01  Bath             Standardized (from CCM1)
!   92.08.01  Hack/Williamson  Reviewed
!   96.03.01  Acker            Modified
!   96.04.01  Hack/Williamson  Reviewed
!   02.04.04  Sawyer           Not needed by FV -- turned into stub
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
   return
!EOC
end subroutine trunc
!-----------------------------------------------------------------------
