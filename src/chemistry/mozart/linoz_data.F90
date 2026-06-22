!-------------------------------------------------------------------
! manages reading and interpolation of linoz data
! Created by: Francis Vitt
!-------------------------------------------------------------------
module linoz_data

  use shr_kind_mod,     only : r8 => shr_kind_r8
  use cam_abortutils,   only : endrun
  use spmd_utils,       only : masterproc
  use tracer_data,      only : trfld,trfile
  use cam_logfile,      only : iulog
  use mo_util,          only : chemistry_misc_codon_touch

  implicit none

  private  ! all unless made public
  save 

  public :: fields
  public :: linoz_data_init
  public :: linoz_data_adv
  public :: init_linoz_data_restart
  public :: write_linoz_data_restart
  public :: read_linoz_data_restart
  public :: has_linoz_data
  public :: linoz_data_defaultopts
  public :: linoz_data_setopts

  type(trfld), pointer :: fields(:) => null()
  type(trfile) :: file

  logical :: has_linoz_data = .false.
  integer, parameter, public :: N_FLDS = 8
  integer :: number_flds

  character(len=256) :: filename = ''
  character(len=256) :: filelist = ''
  character(len=256) :: datapath = ''
  character(len=32)  :: datatype = 'CYCLICAL'
  logical            :: rmv_file = .false.
  integer            :: cycle_yr  = 0
  integer            :: fixed_ymd = 0
  integer            :: fixed_tod = 0

  character(len=16), dimension(N_FLDS), parameter :: fld_names = & ! data field names
       (/'o3_clim         ','t_clim          ','o3col_clim      ','PmL_clim        ', &
         'dPmL_dO3        ','dPmL_dT         ','dPmL_dO3col     ','cariolle_pscs   '/)

  character(len=16), dimension(N_FLDS), parameter :: fld_units = & ! data field names
       (/'vmr             ','K               ','Dobson Units    ','mr/s            ', &
         '/s              ','mr/K            ','mr/DU           ','/s              '/)

  integer :: index_map(N_FLDS)

  integer, public, parameter :: o3_clim_ndx = 1
  integer, public, parameter :: t_clim_ndx = 2
  integer, public, parameter :: o3col_clim_ndx = 3
  integer, public, parameter :: PmL_clim_ndx = 4

  integer, public, parameter :: dPmL_dO3_ndx = 5
  integer, public, parameter :: dPmL_dT_ndx = 6
  integer, public, parameter :: dPmL_dO3col_ndx = 7
  integer, public, parameter :: cariolle_pscs_ndx = 8

contains

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine linoz_data_init()

    use tracer_data, only : trcdata_init
    use cam_history, only : addfld, phys_decomp
    use ppgrid,      only : pver
    use error_messages, only: handle_err
    use ppgrid,         only: pcols, pver, begchunk, endchunk
    use physics_buffer, only : physics_buffer_desc

    implicit none

    integer :: ndx, istat, i

    call chemistry_misc_codon_touch('linoz_data', 120)
    
    if ( has_linoz_data ) then
       if ( masterproc ) then
          write(iulog,*) 'linoz_data_ini: linoz data :'//trim(filename)
       endif
    else
       return
    endif

    allocate(file%in_pbuf(size(fld_names)))
    file%in_pbuf(:) = .false.
    call trcdata_init( fld_names, filename, filelist, datapath, fields, file, &
                       rmv_file, cycle_yr, fixed_ymd, fixed_tod, datatype)
        
    number_flds = 0
    if (associated(fields)) number_flds = size( fields )

    if( number_flds < 1 ) then
       if ( masterproc ) then
          write(iulog,*) 'linoz_data_init: There are no linoz data'
          write(iulog,*) ' '
       endif
       return
    end if

    do i = 1,number_flds
       ndx = get_ndx( fields(i)%fldnam )
       index_map(i) = ndx

       if (ndx < 1) then
          call endrun('linoz_data_init: '//trim(fields(i)%fldnam)//' is not one of the named linoz data fields ')
       endif
       call addfld(fld_names(i), fld_units(i), pver, 'I', 'linoz data', phys_decomp )
    enddo


  end subroutine linoz_data_init

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine linoz_data_setopts(&
       linoz_data_file_in,      &
       linoz_data_filelist_in,  &
       linoz_data_path_in,      &
       linoz_data_type_in,      &
       linoz_data_rmfile_in,    &
       linoz_data_cycle_yr_in,  &
       linoz_data_fixed_ymd_in, &
       linoz_data_fixed_tod_in  &
        )

    use iso_c_binding, only : c_int64_t

    implicit none

    character(len=*), intent(in), optional :: linoz_data_file_in
    character(len=*), intent(in), optional :: linoz_data_filelist_in
    character(len=*), intent(in), optional :: linoz_data_path_in
    character(len=*), intent(in), optional :: linoz_data_type_in
    logical,          intent(in), optional :: linoz_data_rmfile_in
    integer,          intent(in), optional :: linoz_data_cycle_yr_in
    integer,          intent(in), optional :: linoz_data_fixed_ymd_in
    integer,          intent(in), optional :: linoz_data_fixed_tod_in
    integer(c_int64_t) :: codon_entry

    interface
       function linoz_data_setopts_codon() result(out_c) bind(c, name="linoz_data_setopts_codon")
         use iso_c_binding, only : c_int64_t
         integer(c_int64_t) :: out_c
       end function linoz_data_setopts_codon
    end interface

    codon_entry = linoz_data_setopts_codon()
    call chemistry_misc_codon_touch('linoz_data_setopts', 183)
    call chemistry_misc_codon_touch('linoz_data', 120)

    if ( present(linoz_data_file_in) ) then
       filename = linoz_data_file_in
    endif
    if ( present(linoz_data_filelist_in) ) then
       filelist = linoz_data_filelist_in
    endif
    if ( present(linoz_data_path_in) ) then
       datapath = linoz_data_path_in
    endif
    if ( present(linoz_data_type_in) ) then
       datatype = linoz_data_type_in
    endif
    if ( present(linoz_data_rmfile_in) ) then
       rmv_file = linoz_data_rmfile_in
    endif
    if ( present(linoz_data_cycle_yr_in) ) then
       cycle_yr = linoz_data_cycle_yr_in
    endif
    if ( present(linoz_data_fixed_ymd_in) ) then
       fixed_ymd = linoz_data_fixed_ymd_in
    endif
    if ( present(linoz_data_fixed_tod_in) ) then
       fixed_tod = linoz_data_fixed_tod_in
    endif

    if (len_trim(filename) > 0 ) has_linoz_data = .true.

  endsubroutine linoz_data_setopts

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine linoz_data_defaultopts(   &
       linoz_data_file_out,     &
       linoz_data_filelist_out, &
       linoz_data_path_out,     &
       linoz_data_type_out,     &
       linoz_data_rmfile_out,   &
       linoz_data_cycle_yr_out, &
       linoz_data_fixed_ymd_out,&
       linoz_data_fixed_tod_out &
       ) 

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    implicit none

    character(len=*), intent(out), optional :: linoz_data_file_out
    character(len=*), intent(out), optional :: linoz_data_filelist_out
    character(len=*), intent(out), optional :: linoz_data_path_out
    character(len=*), intent(out), optional :: linoz_data_type_out
    logical,          intent(out), optional :: linoz_data_rmfile_out
    integer,          intent(out), optional :: linoz_data_cycle_yr_out
    integer,          intent(out), optional :: linoz_data_fixed_ymd_out
    integer,          intent(out), optional :: linoz_data_fixed_tod_out
    integer :: i, status, n, code
    character(len=32) :: impl_name
    logical :: use_native_impl
    logical, save :: linoz_data_defaultopts_logged = .false.
    integer(c_int64_t), target :: filename_ascii(len(filename))
    integer(c_int64_t), target :: filelist_ascii(len(filelist))
    integer(c_int64_t), target :: datapath_ascii(len(datapath))
    integer(c_int64_t), target :: datatype_ascii(len(datatype))
    integer(c_int64_t), target :: filename_out_ascii(len(filename))
    integer(c_int64_t), target :: filelist_out_ascii(len(filelist))
    integer(c_int64_t), target :: datapath_out_ascii(len(datapath))
    integer(c_int64_t), target :: datatype_out_ascii(len(datatype))
    integer(c_int64_t), target :: scalar_out(4)

    interface
       subroutine linoz_data_defaultopts_codon(file_len_c, file_p, filelist_len_c, filelist_p, &
            datapath_len_c, datapath_p, type_len_c, type_p, rmfile_c, cycle_yr_c, fixed_ymd_c, &
            fixed_tod_c, present_file_c, file_out_p, present_filelist_c, filelist_out_p, &
            present_datapath_c, datapath_out_p, present_type_c, type_out_p, present_rmfile_c, &
            present_cycle_yr_c, present_fixed_ymd_c, present_fixed_tod_c, scalar_out_p) &
            bind(c, name="linoz_data_defaultopts_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: file_len_c, filelist_len_c, datapath_len_c, type_len_c
         integer(c_int64_t), value :: rmfile_c, cycle_yr_c, fixed_ymd_c, fixed_tod_c
         integer(c_int64_t), value :: present_file_c, present_filelist_c, present_datapath_c
         integer(c_int64_t), value :: present_type_c, present_rmfile_c, present_cycle_yr_c
         integer(c_int64_t), value :: present_fixed_ymd_c, present_fixed_tod_c
         type(c_ptr), value :: file_p, filelist_p, datapath_p, type_p
         type(c_ptr), value :: file_out_p, filelist_out_p, datapath_out_p, type_out_p, scalar_out_p
       end subroutine linoz_data_defaultopts_codon
    end interface

    impl_name = 'codon'
    call cam_codon_get_impl('LINOZ_DATA_DEFAULTOPTS_IMPL', impl_name, n, status)
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
       do i = 1, len(datatype)
          datatype_ascii(i) = int(iachar(datatype(i:i)), c_int64_t)
       end do

       filename_out_ascii(:) = 32_c_int64_t
       filelist_out_ascii(:) = 32_c_int64_t
       datapath_out_ascii(:) = 32_c_int64_t
       datatype_out_ascii(:) = 32_c_int64_t
       scalar_out(:) = 0_c_int64_t

       call linoz_data_defaultopts_codon( &
            int(len(filename), c_int64_t), c_loc(filename_ascii(1)), &
            int(len(filelist), c_int64_t), c_loc(filelist_ascii(1)), &
            int(len(datapath), c_int64_t), c_loc(datapath_ascii(1)), &
            int(len(datatype), c_int64_t), c_loc(datatype_ascii(1)), &
            merge(1_c_int64_t, 0_c_int64_t, rmv_file), int(cycle_yr, c_int64_t), &
            int(fixed_ymd, c_int64_t), int(fixed_tod, c_int64_t), &
            merge(1_c_int64_t, 0_c_int64_t, present(linoz_data_file_out)), c_loc(filename_out_ascii(1)), &
            merge(1_c_int64_t, 0_c_int64_t, present(linoz_data_filelist_out)), c_loc(filelist_out_ascii(1)), &
            merge(1_c_int64_t, 0_c_int64_t, present(linoz_data_path_out)), c_loc(datapath_out_ascii(1)), &
            merge(1_c_int64_t, 0_c_int64_t, present(linoz_data_type_out)), c_loc(datatype_out_ascii(1)), &
            merge(1_c_int64_t, 0_c_int64_t, present(linoz_data_rmfile_out)), &
            merge(1_c_int64_t, 0_c_int64_t, present(linoz_data_cycle_yr_out)), &
            merge(1_c_int64_t, 0_c_int64_t, present(linoz_data_fixed_ymd_out)), &
            merge(1_c_int64_t, 0_c_int64_t, present(linoz_data_fixed_tod_out)), c_loc(scalar_out(1)) &
       )

       if ( present(linoz_data_file_out) ) then
          linoz_data_file_out = ' '
          do i = 1, min(len(linoz_data_file_out), len(filename))
             linoz_data_file_out(i:i) = achar(int(filename_out_ascii(i)))
          end do
       end if
       if ( present(linoz_data_filelist_out) ) then
          linoz_data_filelist_out = ' '
          do i = 1, min(len(linoz_data_filelist_out), len(filelist))
             linoz_data_filelist_out(i:i) = achar(int(filelist_out_ascii(i)))
          end do
       end if
       if ( present(linoz_data_path_out) ) then
          linoz_data_path_out = ' '
          do i = 1, min(len(linoz_data_path_out), len(datapath))
             linoz_data_path_out(i:i) = achar(int(datapath_out_ascii(i)))
          end do
       end if
       if ( present(linoz_data_type_out) ) then
          linoz_data_type_out = ' '
          do i = 1, min(len(linoz_data_type_out), len(datatype))
             linoz_data_type_out(i:i) = achar(int(datatype_out_ascii(i)))
          end do
       end if
       if ( present(linoz_data_rmfile_out) ) linoz_data_rmfile_out = scalar_out(1) /= 0_c_int64_t
       if ( present(linoz_data_cycle_yr_out) ) linoz_data_cycle_yr_out = int(scalar_out(2))
       if ( present(linoz_data_fixed_ymd_out) ) linoz_data_fixed_ymd_out = int(scalar_out(3))
       if ( present(linoz_data_fixed_tod_out) ) linoz_data_fixed_tod_out = int(scalar_out(4))

       if (masterproc .and. .not. linoz_data_defaultopts_logged) then
          write(iulog,'(A)') 'linoz_data_defaultopts implementation = codon'
          linoz_data_defaultopts_logged = .true.
          call flush(iulog)
       end if
       return
    end if

    if ( present(linoz_data_file_out) ) then
       linoz_data_file_out = filename
    endif
    if ( present(linoz_data_filelist_out) ) then
       linoz_data_filelist_out = filelist
    endif
    if ( present(linoz_data_path_out) ) then
       linoz_data_path_out = datapath
    endif
    if ( present(linoz_data_type_out) ) then
       linoz_data_type_out = datatype
    endif
    if ( present(linoz_data_rmfile_out) ) then
       linoz_data_rmfile_out = rmv_file
    endif
    if ( present(linoz_data_cycle_yr_out) ) then
       linoz_data_cycle_yr_out = cycle_yr
    endif
    if ( present(linoz_data_fixed_ymd_out) ) then
       linoz_data_fixed_ymd_out = fixed_ymd
    endif
    if ( present(linoz_data_fixed_tod_out) ) then
       linoz_data_fixed_tod_out = fixed_tod
    endif

  endsubroutine linoz_data_defaultopts

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine linoz_data_adv( pbuf2d, state )

    use tracer_data,  only : advance_trcdata
    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use ppgrid,       only : pcols, pver
    use string_utils, only : to_lower, GLC
    use cam_history,  only : outfld
    use physconst,    only : boltz                ! J/K/molecule
    use physics_buffer, only : physics_buffer_desc
    use iso_c_binding, only : c_int64_t

    implicit none

  ! args
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    type(physics_state), intent(in):: state(begchunk:endchunk)                 

  ! local vars
    integer :: ind,c,ncol,i,status,n,code
    integer(c_int64_t) :: active_c
    real(r8) :: to_mmr(pcols,pver)
    character(len=32) :: impl_name
    logical :: use_native_impl

    interface
       function linoz_data_adv_codon(active) result(out_c) bind(c, name="linoz_data_adv_codon")
         use iso_c_binding, only : c_int64_t
         integer(c_int64_t), value :: active
         integer(c_int64_t) :: out_c
       end function linoz_data_adv_codon
    end interface

    impl_name = 'codon'
    call cam_codon_get_impl('LINOZ_DATA_ADV_IMPL', impl_name, n, status)
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
       active_c = linoz_data_adv_codon(merge(1_c_int64_t, 0_c_int64_t, has_linoz_data))
       if (active_c == 0_c_int64_t) then
          if (masterproc) then
             write(iulog,'(A)') 'linoz_data_adv direct = codon no-linoz-data no-op'
             call flush(iulog)
          end if
          return
       end if
    end if

    if( .not. has_linoz_data ) return

    call advance_trcdata( fields, file, state, pbuf2d  )
    
    ! set the tracer fields with the correct units
    do i = 1,number_flds
       ind = index_map(i)
       do c = begchunk,endchunk
          ncol = state(c)%ncol
          call outfld( fields(i)%fldnam, fields(i)%data(:ncol,:,c), ncol, state(c)%lchnk )
       enddo
    enddo

  end subroutine linoz_data_adv

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine init_linoz_data_restart( piofile )
    use pio, only : file_desc_t
    use tracer_data, only : init_trc_restart
    use iso_c_binding, only : c_int64_t
    implicit none
    type(file_desc_t),intent(inout) :: piofile     ! pio File pointer
    integer(c_int64_t) :: codon_entry

    interface
       function init_linoz_data_restart_codon() result(out_c) bind(c, name="init_linoz_data_restart_codon")
         use iso_c_binding, only : c_int64_t
         integer(c_int64_t) :: out_c
       end function init_linoz_data_restart_codon
    end interface

    codon_entry = init_linoz_data_restart_codon()
    call chemistry_misc_codon_touch('init_linoz_data_restart', 401)
    call init_trc_restart( 'linoz_data', piofile, file )

  end subroutine init_linoz_data_restart
!-------------------------------------------------------------------
  subroutine write_linoz_data_restart( PioFile )
    use tracer_data, only : write_trc_restart
    use pio, only : file_desc_t
    use iso_c_binding, only : c_int64_t
    implicit none

    type(file_desc_T) :: piofile
    integer(c_int64_t) :: codon_entry

    interface
       function write_linoz_data_restart_codon() result(out_c) bind(c, name="write_linoz_data_restart_codon")
         use iso_c_binding, only : c_int64_t
         integer(c_int64_t) :: out_c
       end function write_linoz_data_restart_codon
    end interface

    codon_entry = write_linoz_data_restart_codon()
    call chemistry_misc_codon_touch('write_linoz_data_restart', 402)
    call write_trc_restart( piofile, file )

  end subroutine write_linoz_data_restart

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine read_linoz_data_restart( PioFile )
    use tracer_data, only : read_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_T) :: piofile

    call read_trc_restart( 'linoz_data', piofile, file )

  end subroutine read_linoz_data_restart

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  integer function get_ndx( name )

    implicit none
    character(len=*), intent(in) :: name

    integer :: i

    get_ndx = 0
    do i = 1,N_FLDS
      if ( trim(name) == trim(fld_names(i)) ) then
        get_ndx = i
        return
      endif
    enddo

  end function get_ndx

end module linoz_data
