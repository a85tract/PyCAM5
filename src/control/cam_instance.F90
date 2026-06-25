module cam_instance

use seq_comm_mct, only: seq_comm_suffix, seq_comm_inst, seq_comm_name
use iso_c_binding, only: c_int64_t, c_loc, c_ptr
use cam_logfile, only: iulog

implicit none
private
save

public :: cam_instance_init
public :: cam_instance_misc_touch

integer,           public, target :: atm_id
integer,           public, target :: inst_index
character(len=16), public, target :: inst_name
character(len=16), public, target :: inst_suffix

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
      subroutine cam_instance_init_codon(in_atm_id_c, inst_index_in_c, name_in_p, suffix_in_p, &
           atm_id_p, inst_index_p, inst_name_p, inst_suffix_p, name_len_c) &
           bind(c, name='cam_instance_init_codon')
        import :: c_int64_t, c_ptr
        integer(c_int64_t), value :: in_atm_id_c, inst_index_in_c, name_len_c
        type(c_ptr), value :: name_in_p, suffix_in_p, atm_id_p, inst_index_p, inst_name_p, inst_suffix_p
      end subroutine cam_instance_init_codon
   end interface
   character(len=32) :: impl_name
   integer :: n, status
   character(len=16), target :: inst_name_in, inst_suffix_in
   integer :: inst_index_in

   impl_name = 'codon'
   call cam_codon_get_impl('CAM_INSTANCE_INIT_IMPL', impl_name, n, status)
   if (.not. (status == 0 .and. n > 0 .and. trim(adjustl(impl_name(:n))) == 'native')) then
      inst_name_in = seq_comm_name(in_atm_id)
      inst_index_in = seq_comm_inst(in_atm_id)
      inst_suffix_in = seq_comm_suffix(in_atm_id)
      call cam_instance_init_codon(int(in_atm_id, c_int64_t), int(inst_index_in, c_int64_t), &
           c_loc(inst_name_in(1:1)), c_loc(inst_suffix_in(1:1)), c_loc(atm_id), c_loc(inst_index), &
           c_loc(inst_name(1:1)), c_loc(inst_suffix(1:1)), int(len(inst_name), c_int64_t))
      write(iulog,*) 'cam_instance_init implementation = codon'
      return
   endif

   atm_id      = in_atm_id
   inst_name   = seq_comm_name(atm_id)
   inst_index  = seq_comm_inst(atm_id)
   inst_suffix = seq_comm_suffix(atm_id)

end subroutine cam_instance_init

end module cam_instance
