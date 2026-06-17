subroutine cam_codon_get_impl(selector_name, impl_name, impl_len, impl_status)
  implicit none

  character(len=*), intent(in) :: selector_name
  character(len=*), intent(out) :: impl_name
  integer, intent(out) :: impl_len
  integer, intent(out) :: impl_status

  character(len=len(impl_name)) :: candidate
  integer :: candidate_len, candidate_status

  impl_name = 'codon'
  impl_len = len_trim(impl_name)
  impl_status = 1

  candidate = ''
  call get_environment_variable('CAM_CODON_IMPL', value=candidate, &
       length=candidate_len, status=candidate_status)
  if (candidate_status == 0 .and. candidate_len > 0) then
     call normalize_impl_name(candidate, candidate_len, impl_name, impl_len)
     impl_status = 0
  end if

  candidate = ''
  call get_environment_variable(selector_name, value=candidate, &
       length=candidate_len, status=candidate_status)
  if (candidate_status == 0 .and. candidate_len > 0) then
     call normalize_impl_name(candidate, candidate_len, impl_name, impl_len)
     impl_status = 0
  end if

contains

  subroutine normalize_impl_name(raw_name, raw_len, normalized_name, normalized_len)
    implicit none

    character(len=*), intent(in) :: raw_name
    integer, intent(in) :: raw_len
    character(len=*), intent(out) :: normalized_name
    integer, intent(out) :: normalized_len

    character(len=len(normalized_name)) :: lowered
    integer :: i, code, n

    lowered = ''
    normalized_name = ''
    n = min(max(raw_len, 0), len(raw_name), len(lowered))
    if (n > 0) then
       lowered(:n) = raw_name(:n)
    end if

    do i = 1, n
       code = iachar(lowered(i:i))
       if (code >= iachar('A') .and. code <= iachar('Z')) then
          lowered(i:i) = achar(code + iachar('a') - iachar('A'))
       end if
    end do

    lowered = adjustl(lowered)

    select case (trim(lowered))
    case ('1', 't', 'true', 'y', 'yes', 'on', 'codon')
       normalized_name = 'codon'
    case ('0', 'f', 'false', 'n', 'no', 'off', 'native', 'fortran')
       normalized_name = 'native'
    case default
       normalized_name = lowered
    end select

    normalized_len = len_trim(normalized_name)
  end subroutine normalize_impl_name

end subroutine cam_codon_get_impl
