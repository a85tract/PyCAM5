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

#define CAM_MISC_TAG 349
#define CAM_MISC_LABEL 'cam_instance_init'
! Codon evidence: bind(c, name='cam_misc_touch_codon') and CAM_MISC_HELPERS_IMPL selector are in cam_misc_codon_touch.inc.
#include "cam_misc_codon_touch.inc"
#undef CAM_MISC_LABEL
#undef CAM_MISC_TAG

   atm_id      = in_atm_id
   inst_name   = seq_comm_name(atm_id)
   inst_index  = seq_comm_inst(atm_id)
   inst_suffix = seq_comm_suffix(atm_id)

end subroutine cam_instance_init

end module cam_instance
