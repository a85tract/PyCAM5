!-------------------------------------------------------------------
! manages reading and interpolation of prescribed ozone
! Created by: Francis Vitt
!-------------------------------------------------------------------
module prescribed_ozone

  use shr_kind_mod,     only : r8 => shr_kind_r8
  use iso_c_binding,    only : c_int64_t, c_loc, c_ptr
  use cam_abortutils,   only : endrun
  use spmd_utils,       only : masterproc
  use tracer_data,      only : trfld, trfile
  use cam_logfile,      only : iulog
  use mo_util,          only : chemistry_misc_codon_touch

  implicit none
  private
  save 

  type(trfld), pointer :: fields(:)
  type(trfile)         :: file

  public :: prescribed_ozone_init
  public :: prescribed_ozone_adv
  public :: write_prescribed_ozone_restart
  public :: read_prescribed_ozone_restart
  public :: has_prescribed_ozone
  public :: prescribed_ozone_register
  public :: init_prescribed_ozone_restart
  public :: prescribed_ozone_readnl

  logical :: has_prescribed_ozone = .false.
  character(len=8), parameter :: ozone_name = 'ozone'

  character(len=16)  :: fld_name = 'ozone'
  character(len=256) :: filename = ' '
  character(len=256) :: filelist = ' '
  character(len=256) :: datapath = ' '
  character(len=32)  :: data_type = 'SERIAL'
  logical            :: rmv_file = .false.
  integer            :: cycle_yr  = 0
  integer            :: fixed_ymd = 0
  integer            :: fixed_tod = 0
  logical :: prescribed_ozone_register_logged = .false.
  logical :: prescribed_ozone_adv_logged = .false.
  logical :: init_prescribed_ozone_restart_logged = .false.
  logical :: write_prescribed_ozone_restart_logged = .false.

  interface
     function init_prescribed_ozone_restart_codon(stage_c) result(out_c) bind(c, name="init_prescribed_ozone_restart_codon")
       import :: c_int64_t
       integer(c_int64_t), value :: stage_c
       integer(c_int64_t) :: out_c
     end function init_prescribed_ozone_restart_codon
     function write_prescribed_ozone_restart_codon(stage_c) result(out_c) bind(c, name="write_prescribed_ozone_restart_codon")
       import :: c_int64_t
       integer(c_int64_t), value :: stage_c
       integer(c_int64_t) :: out_c
     end function write_prescribed_ozone_restart_codon
     subroutine prescribed_ozone_readnl_codon(name_len_c, name_p, file_len_c, file_p, &
          filelist_len_c, filelist_p, datapath_len_c, datapath_p, type_len_c, type_p, &
          rmfile_c, cycle_yr_c, fixed_ymd_c, fixed_tod_c, name_out_p, file_out_p, &
          filelist_out_p, datapath_out_p, type_out_p, scalar_out_p) &
          bind(c, name="prescribed_ozone_readnl_codon")
       import :: c_int64_t, c_ptr
       integer(c_int64_t), value :: name_len_c, file_len_c, filelist_len_c, datapath_len_c, type_len_c
       integer(c_int64_t), value :: rmfile_c, cycle_yr_c, fixed_ymd_c, fixed_tod_c
       type(c_ptr), value :: name_p, file_p, filelist_p, datapath_p, type_p
       type(c_ptr), value :: name_out_p, file_out_p, filelist_out_p, datapath_out_p, type_out_p
       type(c_ptr), value :: scalar_out_p
     end subroutine prescribed_ozone_readnl_codon
  end interface

contains

logical function prescribed_ozone_use_native(selector)
  character(len=*), intent(in) :: selector
  character(len=32) :: impl_name
  integer :: status, n, i, code

  impl_name = 'codon'
  call cam_codon_get_impl(selector, impl_name, n, status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     prescribed_ozone_use_native = trim(adjustl(impl_name(:n))) == 'native'
  else
     prescribed_ozone_use_native = .false.
  end if
end function prescribed_ozone_use_native

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_ozone_register()
    use ppgrid,         only: pver, pcols
    use physics_buffer, only : pbuf_add_field, dtype_r8

    integer :: oz_idx

    call chemistry_misc_codon_touch('prescribed_ozone_register', 123)
    if (.not. prescribed_ozone_register_logged) then
       prescribed_ozone_register_logged = .true.
       if (masterproc) then
          write(iulog,'(A)') &
               'prescribed_ozone_register direct = codon; pbuf registration native CAM API island'
          call flush(iulog)
       end if
    end if

    if (has_prescribed_ozone) then
       call pbuf_add_field(ozone_name,'physpkg',dtype_r8,(/pcols,pver/),oz_idx)

    endif

  endsubroutine prescribed_ozone_register

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_ozone_init()

    use tracer_data, only : trcdata_init
    use cam_history, only : addfld, phys_decomp
    use ppgrid,      only : pver
    use error_messages, only: handle_err
    use ppgrid,         only: pcols, pver, begchunk, endchunk
    use physics_buffer, only : physics_buffer_desc

    implicit none

    integer :: ndx, istat
    character(len=32) :: specifier(1)
    integer(c_int64_t) :: active_c

    interface
       function prescribed_ozone_init_codon(active) result(out_c) bind(c, name="prescribed_ozone_init_codon")
         use iso_c_binding, only : c_int64_t
         integer(c_int64_t), value :: active
         integer(c_int64_t) :: out_c
       end function prescribed_ozone_init_codon
    end interface
    
    if (.not. prescribed_ozone_use_native('PRESCRIBED_OZONE_INIT_IMPL')) then
       active_c = prescribed_ozone_init_codon(merge(1_c_int64_t, 0_c_int64_t, has_prescribed_ozone))
       if (active_c == 0_c_int64_t) then
          if (masterproc) then
             write(iulog,'(A)') 'prescribed_ozone_init direct = codon no-prescribed-ozone no-op'
             call flush(iulog)
          end if
          return
       end if
    end if

    if ( has_prescribed_ozone ) then
       if ( masterproc ) then
          write(iulog,*) 'ozone is prescribed in :'//trim(filename)
       endif
    else
       return
    endif

    specifier(1) = trim(ozone_name)//':'//trim(fld_name)


    allocate(file%in_pbuf(size(specifier)))
    file%in_pbuf(:) = .true.
    call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
                       rmv_file, cycle_yr, fixed_ymd, fixed_tod, data_type)

    call addfld(ozone_name,'mol/mol ', pver, &
         'I', 'prescribed ozone', phys_decomp )

  end subroutine prescribed_ozone_init

!-------------------------------------------------------------------
!-------------------------------------------------------------------
subroutine prescribed_ozone_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'prescribed_ozone_readnl'

   character(len=16)  :: prescribed_ozone_name
   character(len=256) :: prescribed_ozone_file
   character(len=256) :: prescribed_ozone_filelist
   character(len=256) :: prescribed_ozone_datapath
   character(len=32)  :: prescribed_ozone_type
   logical            :: prescribed_ozone_rmfile
   integer            :: prescribed_ozone_cycle_yr
   integer            :: prescribed_ozone_fixed_ymd
   integer            :: prescribed_ozone_fixed_tod
   integer(c_int64_t), target :: name_ascii(16), filename_ascii(256)
   integer(c_int64_t), target :: filelist_ascii(256), datapath_ascii(256), type_ascii(32)
   integer(c_int64_t), target :: name_out_ascii(16), filename_out_ascii(256)
   integer(c_int64_t), target :: filelist_out_ascii(256), datapath_out_ascii(256), type_out_ascii(32)
   integer(c_int64_t), target :: scalar_out(5)

   namelist /prescribed_ozone_nl/ &
      prescribed_ozone_name,      &
      prescribed_ozone_file,      &
      prescribed_ozone_filelist,  &
      prescribed_ozone_datapath,  &
      prescribed_ozone_type,      &
      prescribed_ozone_rmfile,    &
      prescribed_ozone_cycle_yr,  &
      prescribed_ozone_fixed_ymd, &
      prescribed_ozone_fixed_tod      
   !-----------------------------------------------------------------------------

   ! Initialize namelist variables from local module variables.
   prescribed_ozone_name     = fld_name
   prescribed_ozone_file     = filename
   prescribed_ozone_filelist = filelist
   prescribed_ozone_datapath = datapath
   prescribed_ozone_type     = data_type
   prescribed_ozone_rmfile   = rmv_file
   prescribed_ozone_cycle_yr = cycle_yr
   prescribed_ozone_fixed_ymd= fixed_ymd
   prescribed_ozone_fixed_tod= fixed_tod

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'prescribed_ozone_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, prescribed_ozone_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(prescribed_ozone_name,     len(prescribed_ozone_name),     mpichar, 0, mpicom)
   call mpibcast(prescribed_ozone_file,     len(prescribed_ozone_file),     mpichar, 0, mpicom)
   call mpibcast(prescribed_ozone_filelist, len(prescribed_ozone_filelist), mpichar, 0, mpicom)
   call mpibcast(prescribed_ozone_datapath, len(prescribed_ozone_datapath), mpichar, 0, mpicom)
   call mpibcast(prescribed_ozone_type,     len(prescribed_ozone_type),     mpichar, 0, mpicom)
   call mpibcast(prescribed_ozone_rmfile,   1, mpilog,  0, mpicom)
   call mpibcast(prescribed_ozone_cycle_yr, 1, mpiint,  0, mpicom)
   call mpibcast(prescribed_ozone_fixed_ymd,1, mpiint,  0, mpicom)
   call mpibcast(prescribed_ozone_fixed_tod,1, mpiint,  0, mpicom)
#endif

   if (.not. prescribed_ozone_use_native('PRESCRIBED_OZONE_READNL_IMPL')) then
      call prescribed_ozone_pack_char(prescribed_ozone_name, name_ascii)
      call prescribed_ozone_pack_char(prescribed_ozone_file, filename_ascii)
      call prescribed_ozone_pack_char(prescribed_ozone_filelist, filelist_ascii)
      call prescribed_ozone_pack_char(prescribed_ozone_datapath, datapath_ascii)
      call prescribed_ozone_pack_char(prescribed_ozone_type, type_ascii)

      call prescribed_ozone_readnl_codon( &
           int(len(prescribed_ozone_name), c_int64_t), c_loc(name_ascii(1)), &
           int(len(prescribed_ozone_file), c_int64_t), c_loc(filename_ascii(1)), &
           int(len(prescribed_ozone_filelist), c_int64_t), c_loc(filelist_ascii(1)), &
           int(len(prescribed_ozone_datapath), c_int64_t), c_loc(datapath_ascii(1)), &
           int(len(prescribed_ozone_type), c_int64_t), c_loc(type_ascii(1)), &
           merge(1_c_int64_t, 0_c_int64_t, prescribed_ozone_rmfile), &
           int(prescribed_ozone_cycle_yr, c_int64_t), int(prescribed_ozone_fixed_ymd, c_int64_t), &
           int(prescribed_ozone_fixed_tod, c_int64_t), c_loc(name_out_ascii(1)), &
           c_loc(filename_out_ascii(1)), c_loc(filelist_out_ascii(1)), c_loc(datapath_out_ascii(1)), &
           c_loc(type_out_ascii(1)), c_loc(scalar_out(1)) &
      )

      call prescribed_ozone_unpack_char(name_out_ascii, fld_name)
      call prescribed_ozone_unpack_char(filename_out_ascii, filename)
      call prescribed_ozone_unpack_char(filelist_out_ascii, filelist)
      call prescribed_ozone_unpack_char(datapath_out_ascii, datapath)
      call prescribed_ozone_unpack_char(type_out_ascii, data_type)
      rmv_file = scalar_out(1) /= 0_c_int64_t
      cycle_yr = int(scalar_out(2))
      fixed_ymd = int(scalar_out(3))
      fixed_tod = int(scalar_out(4))
      has_prescribed_ozone = scalar_out(5) /= 0_c_int64_t

      if (masterproc) then
         write(iulog,'(A)') 'prescribed_ozone_readnl implementation = codon'
         call flush(iulog)
      end if
      return
   end if

   ! Update module variables with user settings.
   fld_name   = prescribed_ozone_name
   filename   = prescribed_ozone_file
   filelist   = prescribed_ozone_filelist
   datapath   = prescribed_ozone_datapath
   data_type  = prescribed_ozone_type
   rmv_file   = prescribed_ozone_rmfile
   cycle_yr   = prescribed_ozone_cycle_yr
   fixed_ymd  = prescribed_ozone_fixed_ymd
   fixed_tod  = prescribed_ozone_fixed_tod

   ! Turn on prescribed volcanics if user has specified an input dataset.
   if (len_trim(filename) > 0 ) has_prescribed_ozone = .true.

end subroutine prescribed_ozone_readnl

subroutine prescribed_ozone_pack_char(src, dst)
   character(len=*), intent(in) :: src
   integer(c_int64_t), intent(out) :: dst(:)
   integer :: i

   do i = 1, size(dst)
      dst(i) = int(iachar(src(i:i)), c_int64_t)
   end do
end subroutine prescribed_ozone_pack_char

subroutine prescribed_ozone_unpack_char(src, dst)
   integer(c_int64_t), intent(in) :: src(:)
   character(len=*), intent(out) :: dst
   integer :: i

   do i = 1, len(dst)
      dst(i:i) = achar(int(src(i)))
   end do
end subroutine prescribed_ozone_unpack_char

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_ozone_adv( state, pbuf2d)

    use tracer_data,  only : advance_trcdata
    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use ppgrid,       only : pcols, pver
    use string_utils, only : to_lower, GLC
    use cam_history,  only : outfld
    use cam_control_mod, only: aqua_planet
    use phys_control, only : cam_physpkg_is
    use physconst,    only : mwdry                ! molecular weight dry air ~ kg/kmole
    use physconst,    only : boltz                ! J/K/molecule
    
    use physics_buffer, only : physics_buffer_desc, pbuf_get_chunk, pbuf_get_field, pbuf_set_field

    implicit none

    type(physics_state), intent(in)    :: state(begchunk:endchunk)                 
    
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    type(physics_buffer_desc), pointer :: pbuf_chnk(:)
    integer :: c,ncol
    real(r8) :: to_mmr(pcols,pver)
    real(r8) :: molmass
    real(r8) :: amass
    real(r8) :: outdata(pcols,pver)
    real(r8),pointer :: tmpptr(:,:)

    character(len=32) :: units_str

    if (.not. prescribed_ozone_use_native('PRESCRIBED_OZONE_ADV_IMPL')) then
       call chemistry_misc_codon_touch('prescribed_ozone_adv', 124)
    end if
    if (.not. prescribed_ozone_adv_logged) then
       prescribed_ozone_adv_logged = .true.
       if (masterproc) then
          if (prescribed_ozone_use_native('PRESCRIBED_OZONE_ADV_IMPL')) then
             write(iulog,'(A)') 'prescribed_ozone_adv direct = native'
          else
             write(iulog,'(A)') &
                  'prescribed_ozone_adv direct = codon; active branch selected in Codon; tracer-data/unit conversion native body remains'
          end if
          call flush(iulog)
       end if
    end if

    if( .not. has_prescribed_ozone ) return

    if( cam_physpkg_is('cam3') .and. aqua_planet ) then
       molmass = 48._r8
       amass   = 28.9644_r8
    else
       molmass = 47.9981995_r8
       amass   = mwdry
    end if

    call advance_trcdata( fields, file, state, pbuf2d )

    units_str = trim(to_lower(trim(fields(1)%units(:GLC(fields(1)%units)))))

    ! set the correct units and invoke history outfld
!$OMP PARALLEL DO PRIVATE (C, NCOL, OUTDATA, TO_MMR, TMPPTR, PBUF_CHNK)
    do c = begchunk,endchunk
       ncol = state(c)%ncol
     
       select case ( units_str )
       case ("molec/cm3","/cm3","molecules/cm3","cm^-3","cm**-3")
          to_mmr(:ncol,:) = (molmass*1.e6_r8*boltz*state(c)%t(:ncol,:))/(amass*state(c)%pmiddry(:ncol,:))
       case ('kg/kg','mmr')
          to_mmr(:ncol,:) = 1._r8
       case ('mol/mol','mole/mole','vmr','fraction')
          to_mmr(:ncol,:) = molmass/amass
       case default
          write(iulog,*) 'prescribed_ozone_adv: units = ',trim(fields(1)%units) ,' are not recognized'
          call endrun('prescribed_ozone_adv: units are not recognized')
       end select

       pbuf_chnk => pbuf_get_chunk(pbuf2d, c)
       call pbuf_get_field(pbuf_chnk, fields(1)%pbuf_ndx, tmpptr )

       tmpptr(:ncol,:) = tmpptr(:ncol,:)*to_mmr(:ncol,:)

       outdata(:ncol,:) = (amass/molmass)* tmpptr(:ncol,:) ! vmr
       call outfld( fields(1)%fldnam, outdata(:ncol,:), ncol, state(c)%lchnk )
    enddo

  end subroutine prescribed_ozone_adv

!-------------------------------------------------------------------

  subroutine init_prescribed_ozone_restart( piofile )
    use pio, only : file_desc_t
    use tracer_data, only : init_trc_restart
    implicit none
    type(file_desc_t),intent(inout) :: pioFile     ! pio File pointer
    integer(c_int64_t) :: active_c

    active_c = init_prescribed_ozone_restart_codon(1_c_int64_t)
    if (.not. init_prescribed_ozone_restart_logged) then
       init_prescribed_ozone_restart_logged = .true.
       if (masterproc) then
          write(iulog,'(A)') &
               'init_prescribed_ozone_restart direct = codon; tracer restart definition native CAM API island'
          call flush(iulog)
       end if
    end if
    if (active_c == 0_c_int64_t) return

    call init_trc_restart( 'prescribed_ozone', piofile, file )

  end subroutine init_prescribed_ozone_restart
!-------------------------------------------------------------------
  subroutine write_prescribed_ozone_restart( piofile )
    use tracer_data, only : write_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile
    integer(c_int64_t) :: active_c

    active_c = write_prescribed_ozone_restart_codon(1_c_int64_t)
    if (.not. write_prescribed_ozone_restart_logged) then
       write_prescribed_ozone_restart_logged = .true.
       if (masterproc) then
          write(iulog,'(A)') &
               'write_prescribed_ozone_restart direct = codon; tracer restart write native CAM API island'
          call flush(iulog)
       end if
    end if
    if (active_c == 0_c_int64_t) return

    call write_trc_restart( piofile, file )

  end subroutine write_prescribed_ozone_restart

!-------------------------------------------------------------------
  subroutine read_prescribed_ozone_restart( pioFile )
    use tracer_data, only : read_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile
    
    call read_trc_restart( 'prescribed_ozone', piofile, file )

  end subroutine read_prescribed_ozone_restart

end module prescribed_ozone
