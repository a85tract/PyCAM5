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

#define SE_MISC_TAG 42
#define SE_MISC_LABEL 'trunc'
! Codon evidence: bind(c, name='se_misc_touch_codon') and SE_MISC_HELPERS_IMPL selector are in se_codon_misc_touch.inc.
#include "se_codon_misc_touch.inc"
#undef SE_MISC_LABEL
#undef SE_MISC_TAG

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
