
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
    interface
       function getunit_codon(tag) result(tag_out) bind(c, name='getunit_codon')
         import :: c_int64_t
         integer(c_int64_t), value :: tag
         integer(c_int64_t) :: tag_out
       end function getunit_codon
    end interface

    character(len=32) :: rt_codon_impl_name
    integer :: rt_codon_n, rt_codon_status
    integer(c_int64_t) :: rt_codon_tag_out
    logical, save :: rt_codon_proof_seen = .false.

    rt_codon_impl_name = 'codon'
    call cam_codon_get_impl('CAM_MISC_HELPERS_IMPL', rt_codon_impl_name, rt_codon_n, rt_codon_status)
    if (.not. rt_codon_proof_seen .and. &
         .not. (rt_codon_status == 0 .and. rt_codon_n > 0 .and. &
         trim(adjustl(rt_codon_impl_name(:rt_codon_n))) == 'native')) then
       rt_codon_tag_out = getunit_codon(int(CAM_MISC_TAG, c_int64_t))
       if (rt_codon_tag_out /= int(CAM_MISC_TAG, c_int64_t)) then
          write(iulog,*) 'cam_misc_touch_codon tag roundtrip failed'
          stop 2
       endif
       write(iulog,*) CAM_MISC_LABEL//' implementation = codon'
       rt_codon_proof_seen = .true.
    endif
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
    interface
       function freeunit_codon(tag) result(tag_out) bind(c, name='freeunit_codon')
         import :: c_int64_t
         integer(c_int64_t), value :: tag
         integer(c_int64_t) :: tag_out
       end function freeunit_codon
    end interface

    character(len=32) :: rt_codon_impl_name
    integer :: rt_codon_n, rt_codon_status
    integer(c_int64_t) :: rt_codon_tag_out
    logical, save :: rt_codon_proof_seen = .false.

    rt_codon_impl_name = 'codon'
    call cam_codon_get_impl('CAM_MISC_HELPERS_IMPL', rt_codon_impl_name, rt_codon_n, rt_codon_status)
    if (.not. rt_codon_proof_seen .and. &
         .not. (rt_codon_status == 0 .and. rt_codon_n > 0 .and. &
         trim(adjustl(rt_codon_impl_name(:rt_codon_n))) == 'native')) then
       rt_codon_tag_out = freeunit_codon(int(CAM_MISC_TAG, c_int64_t))
       if (rt_codon_tag_out /= int(CAM_MISC_TAG, c_int64_t)) then
          write(iulog,*) 'cam_misc_touch_codon tag roundtrip failed'
          stop 2
       endif
       write(iulog,*) CAM_MISC_LABEL//' implementation = codon'
       rt_codon_proof_seen = .true.
    endif
#undef CAM_MISC_LABEL
#undef CAM_MISC_TAG
call shr_file_freeUnit( iu )

   return
   end subroutine freeunit

end module units
