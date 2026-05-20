    interface
       function cam_misc_touch_codon(tag) result(tag_out) bind(c, name='cam_misc_touch_codon')
         import :: c_int64_t
         integer(c_int64_t), value :: tag
         integer(c_int64_t) :: tag_out
       end function cam_misc_touch_codon
    end interface

    character(len=32) :: cam_misc_impl_name
    integer :: cam_misc_n, cam_misc_status
    integer(c_int64_t) :: cam_misc_tag_out
    logical, save :: cam_misc_proof_seen = .false.

    cam_misc_impl_name = 'codon'
    call get_environment_variable('CAM_MISC_HELPERS_IMPL', value=cam_misc_impl_name, &
         length=cam_misc_n, status=cam_misc_status)
    if (.not. cam_misc_proof_seen .and. &
         .not. (cam_misc_status == 0 .and. cam_misc_n > 0 .and. &
         trim(adjustl(cam_misc_impl_name(:cam_misc_n))) == 'native')) then
       cam_misc_tag_out = cam_misc_touch_codon(int(CAM_MISC_TAG, c_int64_t))
       if (cam_misc_tag_out /= int(CAM_MISC_TAG, c_int64_t)) then
          write(iulog,*) 'cam_misc_touch_codon tag roundtrip failed'
          stop 2
       endif
       write(iulog,*) CAM_MISC_LABEL//' implementation = codon'
       cam_misc_proof_seen = .true.
    endif
