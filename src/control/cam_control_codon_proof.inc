    character(len=32) :: cam_control_impl_name
    integer :: cam_control_n, cam_control_status
    integer(c_int64_t) :: cam_control_tag_out
    logical, save :: cam_control_proof_seen = .false.

    cam_control_impl_name = 'codon'
    call cam_codon_get_impl('CAM_CONTROL_PROOFS_IMPL', cam_control_impl_name, cam_control_n, cam_control_status)
    if (.not. cam_control_proof_seen .and. &
         .not. (cam_control_status == 0 .and. cam_control_n > 0 .and. &
         trim(adjustl(cam_control_impl_name(:cam_control_n))) == 'native')) then
