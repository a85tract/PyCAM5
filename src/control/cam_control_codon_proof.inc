    character(len=32) :: cam_control_impl_name
    integer :: cam_control_n, cam_control_status
    integer(c_int64_t) :: cam_control_tag_out
    logical, save :: cam_control_proof_seen = .false.

    cam_control_impl_name = 'codon'
    call get_environment_variable('CAM_CONTROL_PROOFS_IMPL', value=cam_control_impl_name, &
         length=cam_control_n, status=cam_control_status)
    if (.not. cam_control_proof_seen .and. &
         .not. (cam_control_status == 0 .and. cam_control_n > 0 .and. &
         trim(adjustl(cam_control_impl_name(:cam_control_n))) == 'native')) then
