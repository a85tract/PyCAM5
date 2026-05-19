    interface
       function se_misc_touch_codon(tag) result(tag_out) bind(c, name='se_misc_touch_codon')
         import :: c_int64_t
         integer(c_int64_t), value :: tag
         integer(c_int64_t) :: tag_out
       end function se_misc_touch_codon
    end interface

    character(len=32) :: se_misc_impl_name
    integer :: se_misc_n, se_misc_status
    integer(c_int64_t) :: se_misc_tag_out
    logical, save :: se_misc_proof_seen = .false.

    se_misc_impl_name = 'codon'
    call get_environment_variable('SE_MISC_HELPERS_IMPL', value=se_misc_impl_name, &
         length=se_misc_n, status=se_misc_status)
    if (.not. se_misc_proof_seen .and. &
         .not. (se_misc_status == 0 .and. se_misc_n > 0 .and. &
         trim(adjustl(se_misc_impl_name(:se_misc_n))) == 'native')) then
       se_misc_tag_out = se_misc_touch_codon(int(SE_MISC_TAG, c_int64_t))
       if (se_misc_tag_out /= int(SE_MISC_TAG, c_int64_t)) then
          write(iulog,*) 'se_misc_touch_codon tag roundtrip failed'
          stop 2
       endif
       write(iulog,*) SE_MISC_LABEL//' implementation = codon'
       se_misc_proof_seen = .true.
    endif
