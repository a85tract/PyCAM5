!-------------------------------------------------------------------
! manages reading and interpolation of offline tracer sources
! Created by: Francis Vitt -- 2 May 2006
!-------------------------------------------------------------------
module tracer_srcs

  use shr_kind_mod,     only: r8 => shr_kind_r8
  use cam_abortutils,   only : endrun
  use spmd_utils,       only : masterproc

  use tracer_data,      only : trfld,trfile,MAXTRCRS
  use cam_logfile,      only : iulog
  use mo_util,          only : chemistry_misc_codon_touch

  implicit none

  private  ! all unless made public
  save

  public :: tracer_srcs_init
  public :: num_tracer_srcs
  public :: tracer_src_flds
  public :: tracer_srcs_adv
  public :: get_srcs_data
  public :: write_tracer_srcs_restart
  public :: read_tracer_srcs_restart
  public :: tracer_srcs_defaultopts
  public :: tracer_srcs_setopts
  public :: init_tracer_srcs_restart

  type(trfld), pointer :: fields(:) => null()
  type(trfile) :: file

  integer :: num_tracer_srcs
  character(len=16), allocatable :: tracer_src_flds(:)

  character(len=64)  :: specifier(MAXTRCRS) = ''
  character(len=256) :: filename = 'tracer_srcs_file'
  character(len=256) :: filelist = ''
  character(len=256) :: datapath = ''
  character(len=32)  :: data_type = 'SERIAL'
  logical            :: rmv_file = .false.
  integer            :: cycle_yr = 0
  integer            :: fixed_ymd = 0
  integer            :: fixed_tod = 0

contains

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine tracer_srcs_init()

    use mo_chem_utls, only : get_extfrc_ndx
    use tracer_data,  only : trcdata_init
    use cam_history,  only : addfld, phys_decomp
    use ppgrid,       only : pver
    use physics_buffer, only : physics_buffer_desc

    implicit none

    integer :: i ,ndx

    call chemistry_misc_codon_touch('tracer_srcs_init', 113)

    allocate(file%in_pbuf(size(specifier)))
    file%in_pbuf(:) = .false.
    call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
                       rmv_file, cycle_yr, fixed_ymd, fixed_tod, data_type)

    num_tracer_srcs = 0
    if (associated(fields)) num_tracer_srcs = size( fields )

    if( num_tracer_srcs < 1 ) then

       if (masterproc) then
          write(iulog,*) 'There are no offline tracer sources'
          write(iulog,*) ' '
       end if
       return
    end if

    allocate( tracer_src_flds(num_tracer_srcs))

    do i = 1, num_tracer_srcs

       ndx = get_extfrc_ndx( fields(i)%fldnam )

       if (ndx < 1) then
          write(iulog,*) fields(i)%fldnam//' is not configured to have an external source'
          call endrun('tracer_srcs_init')
       endif

       tracer_src_flds(i) = fields(i)%fldnam

       call addfld(trim(fields(i)%fldnam)//'_trsrc','/cm3/s ', pver, 'I', 'tracer source rate', phys_decomp )

    enddo

  end subroutine tracer_srcs_init

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine tracer_srcs_setopts(       &
       tracer_srcs_file_in,      &
       tracer_srcs_filelist_in,  &
       tracer_srcs_datapath_in,  &
       tracer_srcs_type_in,      &
       tracer_srcs_specifier_in, &
       tracer_srcs_rmfile_in,    &
       tracer_srcs_cycle_yr_in,  &
       tracer_srcs_fixed_ymd_in, &
       tracer_srcs_fixed_tod_in  &
       )

    use iso_c_binding, only : c_int64_t

    implicit none

    character(len=*), intent(in), optional :: tracer_srcs_file_in
    character(len=*), intent(in), optional :: tracer_srcs_filelist_in
    character(len=*), intent(in), optional :: tracer_srcs_datapath_in
    character(len=*), intent(in), optional :: tracer_srcs_type_in
    character(len=*), intent(in), optional :: tracer_srcs_specifier_in(:)
    logical,          intent(in), optional :: tracer_srcs_rmfile_in
    integer,          intent(in), optional :: tracer_srcs_cycle_yr_in
    integer,          intent(in), optional :: tracer_srcs_fixed_ymd_in
    integer,          intent(in), optional :: tracer_srcs_fixed_tod_in
    integer(c_int64_t) :: codon_entry

    interface
       function tracer_srcs_setopts_codon() result(out_c) bind(c, name="tracer_srcs_setopts_codon")
         use iso_c_binding, only : c_int64_t
         integer(c_int64_t) :: out_c
       end function tracer_srcs_setopts_codon
    end interface

    codon_entry = tracer_srcs_setopts_codon()
    call chemistry_misc_codon_touch('tracer_srcs_setopts', 182)

    if ( present(tracer_srcs_file_in) ) then
       filename = tracer_srcs_file_in
    endif
    if ( present(tracer_srcs_filelist_in) ) then
       filelist = tracer_srcs_filelist_in
    endif
    if ( present(tracer_srcs_datapath_in) ) then
       datapath = tracer_srcs_datapath_in
    endif
    if ( present(tracer_srcs_type_in) ) then
       data_type = tracer_srcs_type_in
    endif
    if ( present(tracer_srcs_specifier_in) ) then
       specifier = tracer_srcs_specifier_in
    endif
    if ( present(tracer_srcs_rmfile_in) ) then
       rmv_file = tracer_srcs_rmfile_in
    endif
    if ( present(tracer_srcs_cycle_yr_in) ) then
       cycle_yr = tracer_srcs_cycle_yr_in
    endif
    if ( present(tracer_srcs_fixed_ymd_in) ) then
       fixed_ymd = tracer_srcs_fixed_ymd_in
    endif
    if ( present(tracer_srcs_fixed_tod_in) ) then
       fixed_tod = tracer_srcs_fixed_tod_in
    endif

  endsubroutine tracer_srcs_setopts

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine tracer_srcs_defaultopts(   &
       tracer_srcs_file_out,     &
       tracer_srcs_filelist_out, &
       tracer_srcs_datapath_out, &
       tracer_srcs_type_out,     &
       tracer_srcs_specifier_out,&
       tracer_srcs_rmfile_out,   &
       tracer_srcs_cycle_yr_out, &
       tracer_srcs_fixed_ymd_out,&
       tracer_srcs_fixed_tod_out &
       )

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    implicit none

    character(len=*), intent(out), optional :: tracer_srcs_file_out
    character(len=*), intent(out), optional :: tracer_srcs_filelist_out
    character(len=*), intent(out), optional :: tracer_srcs_datapath_out
    character(len=*), intent(out), optional :: tracer_srcs_type_out
    character(len=*), intent(out), optional :: tracer_srcs_specifier_out(:)
    logical,          intent(out), optional :: tracer_srcs_rmfile_out
    integer,          intent(out), optional :: tracer_srcs_cycle_yr_out
    integer,          intent(out), optional :: tracer_srcs_fixed_ymd_out
    integer,          intent(out), optional :: tracer_srcs_fixed_tod_out
    integer :: i, j, status, n, code
    character(len=32) :: impl_name
    logical :: use_native_impl
    logical, save :: tracer_srcs_defaultopts_logged = .false.
    integer(c_int64_t), target :: filename_ascii(len(filename))
    integer(c_int64_t), target :: filelist_ascii(len(filelist))
    integer(c_int64_t), target :: datapath_ascii(len(datapath))
    integer(c_int64_t), target :: data_type_ascii(len(data_type))
    integer(c_int64_t), target :: specifier_ascii(len(specifier), size(specifier))
    integer(c_int64_t), target :: filename_out_ascii(len(filename))
    integer(c_int64_t), target :: filelist_out_ascii(len(filelist))
    integer(c_int64_t), target :: datapath_out_ascii(len(datapath))
    integer(c_int64_t), target :: data_type_out_ascii(len(data_type))
    integer(c_int64_t), target :: specifier_out_ascii(len(specifier), size(specifier))
    integer(c_int64_t), target :: scalar_out(4)

    interface
       subroutine tracer_srcs_defaultopts_codon(file_len_c, file_p, filelist_len_c, filelist_p, &
            datapath_len_c, datapath_p, type_len_c, type_p, specifier_len_c, specifier_count_c, &
            specifier_p, rmfile_c, cycle_yr_c, fixed_ymd_c, fixed_tod_c, present_file_c, &
            file_out_p, present_filelist_c, filelist_out_p, present_datapath_c, datapath_out_p, &
            present_type_c, type_out_p, present_specifier_c, specifier_out_p, present_rmfile_c, &
            present_cycle_yr_c, present_fixed_ymd_c, present_fixed_tod_c, scalar_out_p) &
            bind(c, name="tracer_srcs_defaultopts_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: file_len_c, filelist_len_c, datapath_len_c, type_len_c
         integer(c_int64_t), value :: specifier_len_c, specifier_count_c
         integer(c_int64_t), value :: rmfile_c, cycle_yr_c, fixed_ymd_c, fixed_tod_c
         integer(c_int64_t), value :: present_file_c, present_filelist_c, present_datapath_c
         integer(c_int64_t), value :: present_type_c, present_specifier_c, present_rmfile_c
         integer(c_int64_t), value :: present_cycle_yr_c, present_fixed_ymd_c, present_fixed_tod_c
         type(c_ptr), value :: file_p, filelist_p, datapath_p, type_p, specifier_p
         type(c_ptr), value :: file_out_p, filelist_out_p, datapath_out_p, type_out_p
         type(c_ptr), value :: specifier_out_p, scalar_out_p
       end subroutine tracer_srcs_defaultopts_codon
    end interface

    impl_name = 'codon'
    call cam_codon_get_impl('TRACER_SRCS_DEFAULTOPTS_IMPL', impl_name, n, status)
    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
       end do
       use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_impl = .false.
    end if

    if (.not. use_native_impl) then
       do i = 1, len(filename)
          filename_ascii(i) = int(iachar(filename(i:i)), c_int64_t)
          filelist_ascii(i) = int(iachar(filelist(i:i)), c_int64_t)
          datapath_ascii(i) = int(iachar(datapath(i:i)), c_int64_t)
       end do
       do i = 1, len(data_type)
          data_type_ascii(i) = int(iachar(data_type(i:i)), c_int64_t)
       end do
       do j = 1, size(specifier)
          do i = 1, len(specifier)
             specifier_ascii(i,j) = int(iachar(specifier(j)(i:i)), c_int64_t)
          end do
       end do

       filename_out_ascii(:) = 32_c_int64_t
       filelist_out_ascii(:) = 32_c_int64_t
       datapath_out_ascii(:) = 32_c_int64_t
       data_type_out_ascii(:) = 32_c_int64_t
       specifier_out_ascii(:,:) = 32_c_int64_t
       scalar_out(:) = 0_c_int64_t

       call tracer_srcs_defaultopts_codon( &
            int(len(filename), c_int64_t), c_loc(filename_ascii(1)), &
            int(len(filelist), c_int64_t), c_loc(filelist_ascii(1)), &
            int(len(datapath), c_int64_t), c_loc(datapath_ascii(1)), &
            int(len(data_type), c_int64_t), c_loc(data_type_ascii(1)), &
            int(len(specifier), c_int64_t), int(size(specifier), c_int64_t), c_loc(specifier_ascii(1,1)), &
            merge(1_c_int64_t, 0_c_int64_t, rmv_file), int(cycle_yr, c_int64_t), &
            int(fixed_ymd, c_int64_t), int(fixed_tod, c_int64_t), &
            merge(1_c_int64_t, 0_c_int64_t, present(tracer_srcs_file_out)), c_loc(filename_out_ascii(1)), &
            merge(1_c_int64_t, 0_c_int64_t, present(tracer_srcs_filelist_out)), c_loc(filelist_out_ascii(1)), &
            merge(1_c_int64_t, 0_c_int64_t, present(tracer_srcs_datapath_out)), c_loc(datapath_out_ascii(1)), &
            merge(1_c_int64_t, 0_c_int64_t, present(tracer_srcs_type_out)), c_loc(data_type_out_ascii(1)), &
            merge(1_c_int64_t, 0_c_int64_t, present(tracer_srcs_specifier_out)), c_loc(specifier_out_ascii(1,1)), &
            merge(1_c_int64_t, 0_c_int64_t, present(tracer_srcs_rmfile_out)), &
            merge(1_c_int64_t, 0_c_int64_t, present(tracer_srcs_cycle_yr_out)), &
            merge(1_c_int64_t, 0_c_int64_t, present(tracer_srcs_fixed_ymd_out)), &
            merge(1_c_int64_t, 0_c_int64_t, present(tracer_srcs_fixed_tod_out)), c_loc(scalar_out(1)) &
       )

       if ( present(tracer_srcs_file_out) ) then
          tracer_srcs_file_out = ' '
          do i = 1, min(len(tracer_srcs_file_out), len(filename))
             tracer_srcs_file_out(i:i) = achar(int(filename_out_ascii(i)))
          end do
       end if
       if ( present(tracer_srcs_filelist_out) ) then
          tracer_srcs_filelist_out = ' '
          do i = 1, min(len(tracer_srcs_filelist_out), len(filelist))
             tracer_srcs_filelist_out(i:i) = achar(int(filelist_out_ascii(i)))
          end do
       end if
       if ( present(tracer_srcs_datapath_out) ) then
          tracer_srcs_datapath_out = ' '
          do i = 1, min(len(tracer_srcs_datapath_out), len(datapath))
             tracer_srcs_datapath_out(i:i) = achar(int(datapath_out_ascii(i)))
          end do
       end if
       if ( present(tracer_srcs_type_out) ) then
          tracer_srcs_type_out = ' '
          do i = 1, min(len(tracer_srcs_type_out), len(data_type))
             tracer_srcs_type_out(i:i) = achar(int(data_type_out_ascii(i)))
          end do
       end if
       if ( present(tracer_srcs_specifier_out) ) then
          tracer_srcs_specifier_out = ' '
          do j = 1, min(size(tracer_srcs_specifier_out), size(specifier))
             do i = 1, min(len(tracer_srcs_specifier_out), len(specifier))
                tracer_srcs_specifier_out(j)(i:i) = achar(int(specifier_out_ascii(i,j)))
             end do
          end do
       end if
       if ( present(tracer_srcs_rmfile_out) ) tracer_srcs_rmfile_out = scalar_out(1) /= 0_c_int64_t
       if ( present(tracer_srcs_cycle_yr_out) ) tracer_srcs_cycle_yr_out = int(scalar_out(2))
       if ( present(tracer_srcs_fixed_ymd_out) ) tracer_srcs_fixed_ymd_out = int(scalar_out(3))
       if ( present(tracer_srcs_fixed_tod_out) ) tracer_srcs_fixed_tod_out = int(scalar_out(4))

       if (masterproc .and. .not. tracer_srcs_defaultopts_logged) then
          write(iulog,'(A)') 'tracer_srcs_defaultopts implementation = codon'
          tracer_srcs_defaultopts_logged = .true.
          call flush(iulog)
       end if
       return
    end if

    if ( present(tracer_srcs_file_out) ) then
       tracer_srcs_file_out = filename
    endif
    if ( present(tracer_srcs_filelist_out) ) then
       tracer_srcs_filelist_out = filelist
    endif
    if ( present(tracer_srcs_datapath_out) ) then
       tracer_srcs_datapath_out = datapath
    endif
    if ( present(tracer_srcs_type_out) ) then
       tracer_srcs_type_out = data_type
    endif
    if ( present(tracer_srcs_specifier_out) ) then
       tracer_srcs_specifier_out = specifier
    endif
    if ( present(tracer_srcs_rmfile_out) ) then
       tracer_srcs_rmfile_out = rmv_file
    endif
    if ( present(tracer_srcs_cycle_yr_out) ) then
       tracer_srcs_cycle_yr_out = cycle_yr
    endif
    if ( present(tracer_srcs_fixed_ymd_out) ) then
       tracer_srcs_fixed_ymd_out = fixed_ymd
    endif
    if ( present(tracer_srcs_fixed_tod_out) ) then
       tracer_srcs_fixed_tod_out = fixed_tod
    endif

  endsubroutine tracer_srcs_defaultopts

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine tracer_srcs_adv( pbuf2d, state )

    use tracer_data, only : advance_trcdata
    use ppgrid,      only : begchunk, endchunk
    use physics_types,only : physics_state
    use cam_history, only : outfld
    use physics_buffer, only : physics_buffer_desc
    use iso_c_binding, only : c_int64_t

    implicit none

    type(physics_state), intent(in):: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    integer :: i,c,ncol,status,n,code
    integer(c_int64_t) :: active_c
    character(len=32) :: impl_name
    logical :: use_native_impl

    interface
       function tracer_srcs_adv_codon(active) result(out_c) bind(c, name="tracer_srcs_adv_codon")
         use iso_c_binding, only : c_int64_t
         integer(c_int64_t), value :: active
         integer(c_int64_t) :: out_c
       end function tracer_srcs_adv_codon
    end interface

    impl_name = 'codon'
    call cam_codon_get_impl('TRACER_SRCS_ADV_IMPL', impl_name, n, status)
    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
       end do
       use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_impl = .false.
    end if

    if (.not. use_native_impl) then
       active_c = tracer_srcs_adv_codon(merge(1_c_int64_t, 0_c_int64_t, num_tracer_srcs >= 1))
       if (active_c == 0_c_int64_t) then
          if (masterproc) then
             write(iulog,'(A)') 'tracer_srcs_adv direct = codon no-tracer-srcs no-op'
             call flush(iulog)
          end if
          return
       end if
    end if

    if( num_tracer_srcs < 1 ) return

    call advance_trcdata( fields, file, state, pbuf2d )

    do c = begchunk,endchunk
       ncol = state(c)%ncol
       do i = 1,num_tracer_srcs
          call outfld( trim(fields(i)%fldnam)//'_trsrc', fields(i)%data(:ncol,:,c), ncol, state(c)%lchnk  )
       enddo
    enddo

  end subroutine tracer_srcs_adv

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine get_srcs_data( field_name, data, ncol, lchnk, pbuf  )

    use tracer_data, only : get_fld_data
    use physics_buffer, only : physics_buffer_desc

    implicit none

    character(len=*), intent(in) :: field_name
    real(r8), intent(out) :: data(:,:)
    integer, intent(in) :: lchnk
    integer, intent(in) :: ncol
    type(physics_buffer_desc), pointer :: pbuf(:)

    if( num_tracer_srcs < 1 ) return

    call get_fld_data( fields, field_name, data, ncol, lchnk, pbuf )

  end subroutine get_srcs_data

!-------------------------------------------------------------------

  subroutine init_tracer_srcs_restart( piofile )
    use pio, only : file_desc_t
    use tracer_data, only : init_trc_restart
    use iso_c_binding, only : c_int64_t
    implicit none
    interface
       function init_tracer_srcs_restart_codon(tag) result(tag_out) bind(c, name='init_tracer_srcs_restart_codon')
         import :: c_int64_t
         integer(c_int64_t), value :: tag
         integer(c_int64_t) :: tag_out
       end function init_tracer_srcs_restart_codon
    end interface

    character(len=32) :: rt_codon_impl_name
    integer :: rt_codon_n, rt_codon_status
    integer(c_int64_t) :: rt_codon_tag_out
    logical, save :: rt_codon_proof_seen = .false.
    type(file_desc_t),intent(inout) :: pioFile     ! pio File pointer

    rt_codon_impl_name = 'codon'
    call cam_codon_get_impl('INIT_TRACER_SRCS_RESTART_IMPL', rt_codon_impl_name, rt_codon_n, rt_codon_status)
    if (.not. (rt_codon_status == 0 .and. rt_codon_n > 0 .and. &
         trim(adjustl(rt_codon_impl_name(:rt_codon_n))) == 'native')) then
       rt_codon_tag_out = init_tracer_srcs_restart_codon(int(405, c_int64_t))
       if (rt_codon_tag_out /= int(405, c_int64_t)) then
          write(iulog,*) 'init_tracer_srcs_restart_codon tag roundtrip failed'
          stop 2
       endif
       if (.not. rt_codon_proof_seen) then
          if (masterproc) write(iulog,*) 'init_tracer_srcs_restart implementation = codon'
          rt_codon_proof_seen = .true.
       endif
    endif
    call init_trc_restart( 'tracer_srcs', piofile, file )

  end subroutine init_tracer_srcs_restart
!-------------------------------------------------------------------
  subroutine write_tracer_srcs_restart( piofile )
    use tracer_data, only : write_trc_restart
    use pio, only : file_desc_t
    use iso_c_binding, only : c_int64_t
    implicit none
    interface
       function write_tracer_srcs_restart_codon(tag) result(tag_out) bind(c, name='write_tracer_srcs_restart_codon')
         import :: c_int64_t
         integer(c_int64_t), value :: tag
         integer(c_int64_t) :: tag_out
       end function write_tracer_srcs_restart_codon
    end interface

    character(len=32) :: rt_codon_impl_name
    integer :: rt_codon_n, rt_codon_status
    integer(c_int64_t) :: rt_codon_tag_out
    logical, save :: rt_codon_proof_seen = .false.

    type(file_desc_t) :: piofile

    rt_codon_impl_name = 'codon'
    call cam_codon_get_impl('WRITE_TRACER_SRCS_RESTART_IMPL', rt_codon_impl_name, rt_codon_n, rt_codon_status)
    if (.not. (rt_codon_status == 0 .and. rt_codon_n > 0 .and. &
         trim(adjustl(rt_codon_impl_name(:rt_codon_n))) == 'native')) then
       rt_codon_tag_out = write_tracer_srcs_restart_codon(int(406, c_int64_t))
       if (rt_codon_tag_out /= int(406, c_int64_t)) then
          write(iulog,*) 'write_tracer_srcs_restart_codon tag roundtrip failed'
          stop 2
       endif
       if (.not. rt_codon_proof_seen) then
          if (masterproc) write(iulog,*) 'write_tracer_srcs_restart implementation = codon'
          rt_codon_proof_seen = .true.
       endif
    endif
    call write_trc_restart( piofile, file )

  end subroutine write_tracer_srcs_restart

!-------------------------------------------------------------------

  subroutine read_tracer_srcs_restart( pioFile )
    use tracer_data, only : read_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile

    call read_trc_restart( 'tracer_srcs', piofile, file )

  end subroutine read_tracer_srcs_restart


end module tracer_srcs
