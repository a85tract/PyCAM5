
module units

   use cam_abortutils, only: endrun
   use shr_file_mod,   only: shr_file_getUnit, shr_file_freeUnit
   use iso_c_binding,  only: c_int64_t
   use cam_logfile,    only: iulog

implicit none

PRIVATE

   public :: getunit, freeunit, units_misc_touch

CONTAINS

   subroutine units_misc_touch()
#define CAM_MISC_TAG 232
#define CAM_MISC_LABEL 'units'
! Codon evidence: bind(c, name='cam_misc_touch_codon') and CAM_MISC_HELPERS_IMPL selector are in cam_misc_codon_touch.inc.
#include "cam_misc_codon_touch.inc"
#undef CAM_MISC_LABEL
#undef CAM_MISC_TAG
   end subroutine units_misc_touch

   integer function getunit (iu)
!
! Arguments
!
   integer, intent(in), optional :: iu   ! desired unit number
!
! Local workspace
!

#define CAM_MISC_TAG 233
#define CAM_MISC_LABEL 'getunit'
! Codon evidence: bind(c, name='cam_misc_touch_codon') and CAM_MISC_HELPERS_IMPL selector are in cam_misc_codon_touch.inc.
#include "cam_misc_codon_touch.inc"
#undef CAM_MISC_LABEL
#undef CAM_MISC_TAG

     getunit = shr_file_getUnit( iu )

   end function getunit

!#######################################################################

   subroutine freeunit (iu)
!
! Arguments
!
   integer, intent(in) :: iu       ! unit number to be freed

#define CAM_MISC_TAG 234
#define CAM_MISC_LABEL 'freeunit'
! Codon evidence: bind(c, name='cam_misc_touch_codon') and CAM_MISC_HELPERS_IMPL selector are in cam_misc_codon_touch.inc.
#include "cam_misc_codon_touch.inc"
#undef CAM_MISC_LABEL
#undef CAM_MISC_TAG

   call shr_file_freeUnit( iu )

   return
   end subroutine freeunit

end module units
