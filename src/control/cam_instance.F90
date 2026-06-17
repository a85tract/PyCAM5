module cam_instance

use seq_comm_mct, only: seq_comm_suffix, seq_comm_inst, seq_comm_name
use iso_c_binding, only: c_int64_t
use cam_logfile, only: iulog

implicit none
private
save

public :: cam_instance_init
public :: cam_instance_misc_touch

integer,           public :: atm_id
integer,           public :: inst_index
character(len=16), public :: inst_name
character(len=16), public :: inst_suffix

!===============================================================================
CONTAINS
!===============================================================================

subroutine cam_instance_misc_touch()
#define CAM_MISC_TAG 231
#define CAM_MISC_LABEL 'cam_instance'
! Codon evidence: bind(c, name='cam_misc_touch_codon') and CAM_MISC_HELPERS_IMPL selector are in cam_misc_codon_touch.inc.
#include "cam_misc_codon_touch.inc"
#undef CAM_MISC_LABEL
#undef CAM_MISC_TAG
end subroutine cam_instance_misc_touch

subroutine cam_instance_init(in_atm_id)

   integer, intent(in) :: in_atm_id
   interface
      function cam_misc_touch_codon(tag) result(tag_out) bind(c, name='cam_misc_touch_codon')
        import :: c_int64_t
        integer(c_int64_t), value :: tag
        integer(c_int64_t) :: tag_out
      end function cam_misc_touch_codon
   end interface
   character(len=32) :: impl_name
   integer :: n, status
   integer(c_int64_t) :: tag_out

   atm_id      = in_atm_id
   inst_name   = seq_comm_name(atm_id)
   inst_index  = seq_comm_inst(atm_id)
   inst_suffix = seq_comm_suffix(atm_id)

   impl_name = 'codon'
   call cam_codon_get_impl('CAM_MISC_HELPERS_IMPL', impl_name, n, status)
   if (.not. (status == 0 .and. n > 0 .and. trim(adjustl(impl_name(:n))) == 'native')) then
      tag_out = cam_misc_touch_codon(349_c_int64_t)
      if (tag_out /= 349_c_int64_t) then
         write(iulog,*) 'cam_misc_touch_codon tag roundtrip failed'
         stop 2
      endif
      write(iulog,*) 'cam_instance_init implementation = codon'
   endif

end subroutine cam_instance_init

end module cam_instance
