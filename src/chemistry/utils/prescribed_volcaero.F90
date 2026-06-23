!-------------------------------------------------------------------
! manages reading and interpolation of prescribed volcanic aerosol
! Created by: Francis Vitt
!-------------------------------------------------------------------
module prescribed_volcaero

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

  public :: prescribed_volcaero_readnl
  public :: prescribed_volcaero_register
  public :: prescribed_volcaero_init
  public :: prescribed_volcaero_adv
  public :: write_prescribed_volcaero_restart
  public :: read_prescribed_volcaero_restart
  public :: has_prescribed_volcaero
  public :: init_prescribed_volcaero_restart


  logical :: has_prescribed_volcaero = .false.
  character(len=8), parameter :: volcaero_name = 'VOLC_MMR'
  character(len=13), parameter :: volcrad_name = 'VOLC_RAD_GEOM'
  character(len=9), parameter :: volcmass_name = 'VOLC_MASS'
  character(len=11), parameter :: volcmass_column_name = 'VOLC_MASS_C'

  ! These variables are settable via the namelist (with longer names)
  character(len=16)  :: fld_name = 'MMRVOLC'
  character(len=256) :: filename = ''
  character(len=256) :: filelist = ''
  character(len=256) :: datapath = ''
  character(len=32)  :: data_type = 'SERIAL'
  logical            :: rmv_file = .false.
  integer            :: cycle_yr  = 0
  integer            :: fixed_ymd = 0
  integer            :: fixed_tod = 0
  integer            :: radius_ndx
  logical :: prescribed_volcaero_adv_logged = .false.
  logical :: init_prescribed_volcaero_restart_logged = .false.
  logical :: write_prescribed_volcaero_restart_logged = .false.
  logical :: prescribed_volcaero_register_logged = .false.

  interface
     function prescribed_volcaero_adv_codon(active_c) result(out_c) bind(c, name="prescribed_volcaero_adv_codon")
       import :: c_int64_t
       integer(c_int64_t), value :: active_c
       integer(c_int64_t) :: out_c
     end function prescribed_volcaero_adv_codon
     function init_prescribed_volcaero_restart_codon(stage_c) result(out_c) bind(c, name="init_prescribed_volcaero_restart_codon")
       import :: c_int64_t
       integer(c_int64_t), value :: stage_c
       integer(c_int64_t) :: out_c
     end function init_prescribed_volcaero_restart_codon
     function write_prescribed_volcaero_restart_codon(stage_c) result(out_c) bind(c, name="write_prescribed_volcaero_restart_codon")
       import :: c_int64_t
       integer(c_int64_t), value :: stage_c
       integer(c_int64_t) :: out_c
     end function write_prescribed_volcaero_restart_codon
     subroutine prescribed_volcaero_readnl_codon(name_len_c, name_p, file_len_c, file_p, &
          filelist_len_c, filelist_p, datapath_len_c, datapath_p, type_len_c, type_p, &
          rmfile_c, cycle_yr_c, fixed_ymd_c, fixed_tod_c, name_out_p, file_out_p, &
          filelist_out_p, datapath_out_p, type_out_p, scalar_out_p) &
          bind(c, name="prescribed_volcaero_readnl_codon")
       import :: c_int64_t, c_ptr
       integer(c_int64_t), value :: name_len_c, file_len_c, filelist_len_c, datapath_len_c, type_len_c
       integer(c_int64_t), value :: rmfile_c, cycle_yr_c, fixed_ymd_c, fixed_tod_c
       type(c_ptr), value :: name_p, file_p, filelist_p, datapath_p, type_p
       type(c_ptr), value :: name_out_p, file_out_p, filelist_out_p, datapath_out_p, type_out_p
       type(c_ptr), value :: scalar_out_p
     end subroutine prescribed_volcaero_readnl_codon
  end interface

contains

logical function prescribed_volcaero_use_native(selector)
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
     prescribed_volcaero_use_native = trim(adjustl(impl_name(:n))) == 'native'
  else
     prescribed_volcaero_use_native = .false.
  end if
end function prescribed_volcaero_use_native

!-------------------------------------------------------------------
!-------------------------------------------------------------------
subroutine prescribed_volcaero_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'prescribed_volcaero_readnl'

   character(len=16)  :: prescribed_volcaero_name
   character(len=256) :: prescribed_volcaero_file
   character(len=256) :: prescribed_volcaero_filelist
   character(len=256) :: prescribed_volcaero_datapath
   character(len=32)  :: prescribed_volcaero_type
   logical            :: prescribed_volcaero_rmfile
   integer            :: prescribed_volcaero_cycle_yr
   integer            :: prescribed_volcaero_fixed_ymd
   integer            :: prescribed_volcaero_fixed_tod
   integer(c_int64_t), target :: name_ascii(16), filename_ascii(256)
   integer(c_int64_t), target :: filelist_ascii(256), datapath_ascii(256), type_ascii(32)
   integer(c_int64_t), target :: name_out_ascii(16), filename_out_ascii(256)
   integer(c_int64_t), target :: filelist_out_ascii(256), datapath_out_ascii(256), type_out_ascii(32)
   integer(c_int64_t), target :: scalar_out(5)

   namelist /prescribed_volcaero_nl/ &
      prescribed_volcaero_name,      &
      prescribed_volcaero_file,      &
      prescribed_volcaero_filelist,  &
      prescribed_volcaero_datapath,  &
      prescribed_volcaero_type,      &
      prescribed_volcaero_rmfile,    &
      prescribed_volcaero_cycle_yr,  &
      prescribed_volcaero_fixed_ymd, &
      prescribed_volcaero_fixed_tod
   !-----------------------------------------------------------------------------

   ! Initialize namelist variables from local module variables.
   prescribed_volcaero_name     = fld_name
   prescribed_volcaero_file     = filename
   prescribed_volcaero_filelist = filelist
   prescribed_volcaero_datapath = datapath
   prescribed_volcaero_type     = data_type
   prescribed_volcaero_rmfile   = rmv_file
   prescribed_volcaero_cycle_yr = cycle_yr
   prescribed_volcaero_fixed_ymd= fixed_ymd
   prescribed_volcaero_fixed_tod= fixed_tod

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'prescribed_volcaero_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, prescribed_volcaero_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(prescribed_volcaero_name,     len(prescribed_volcaero_name),     mpichar, 0, mpicom)
   call mpibcast(prescribed_volcaero_file,     len(prescribed_volcaero_file),     mpichar, 0, mpicom)
   call mpibcast(prescribed_volcaero_filelist, len(prescribed_volcaero_filelist), mpichar, 0, mpicom)
   call mpibcast(prescribed_volcaero_datapath, len(prescribed_volcaero_datapath), mpichar, 0, mpicom)
   call mpibcast(prescribed_volcaero_type,     len(prescribed_volcaero_type),     mpichar, 0, mpicom)
   call mpibcast(prescribed_volcaero_rmfile,   1, mpilog,  0, mpicom)
   call mpibcast(prescribed_volcaero_cycle_yr, 1, mpiint,  0, mpicom)
   call mpibcast(prescribed_volcaero_fixed_ymd,1, mpiint,  0, mpicom)
   call mpibcast(prescribed_volcaero_fixed_tod,1, mpiint,  0, mpicom)
#endif

   if (.not. prescribed_volcaero_use_native('PRESCRIBED_VOLCAERO_READNL_IMPL')) then
      call prescribed_volcaero_pack_char(prescribed_volcaero_name, name_ascii)
      call prescribed_volcaero_pack_char(prescribed_volcaero_file, filename_ascii)
      call prescribed_volcaero_pack_char(prescribed_volcaero_filelist, filelist_ascii)
      call prescribed_volcaero_pack_char(prescribed_volcaero_datapath, datapath_ascii)
      call prescribed_volcaero_pack_char(prescribed_volcaero_type, type_ascii)

      call prescribed_volcaero_readnl_codon( &
           int(len(prescribed_volcaero_name), c_int64_t), c_loc(name_ascii(1)), &
           int(len(prescribed_volcaero_file), c_int64_t), c_loc(filename_ascii(1)), &
           int(len(prescribed_volcaero_filelist), c_int64_t), c_loc(filelist_ascii(1)), &
           int(len(prescribed_volcaero_datapath), c_int64_t), c_loc(datapath_ascii(1)), &
           int(len(prescribed_volcaero_type), c_int64_t), c_loc(type_ascii(1)), &
           merge(1_c_int64_t, 0_c_int64_t, prescribed_volcaero_rmfile), &
           int(prescribed_volcaero_cycle_yr, c_int64_t), int(prescribed_volcaero_fixed_ymd, c_int64_t), &
           int(prescribed_volcaero_fixed_tod, c_int64_t), c_loc(name_out_ascii(1)), &
           c_loc(filename_out_ascii(1)), c_loc(filelist_out_ascii(1)), c_loc(datapath_out_ascii(1)), &
           c_loc(type_out_ascii(1)), c_loc(scalar_out(1)) &
      )

      call prescribed_volcaero_unpack_char(name_out_ascii, fld_name)
      call prescribed_volcaero_unpack_char(filename_out_ascii, filename)
      call prescribed_volcaero_unpack_char(filelist_out_ascii, filelist)
      call prescribed_volcaero_unpack_char(datapath_out_ascii, datapath)
      call prescribed_volcaero_unpack_char(type_out_ascii, data_type)
      rmv_file = scalar_out(1) /= 0_c_int64_t
      cycle_yr = int(scalar_out(2))
      fixed_ymd = int(scalar_out(3))
      fixed_tod = int(scalar_out(4))
      has_prescribed_volcaero = scalar_out(5) /= 0_c_int64_t

      if (masterproc) then
         write(iulog,'(A)') 'prescribed_volcaero_readnl implementation = codon'
         call flush(iulog)
      end if
      return
   end if

   ! Update module variables with user settings.
   fld_name   = prescribed_volcaero_name
   filename   = prescribed_volcaero_file
   filelist   = prescribed_volcaero_filelist
   datapath   = prescribed_volcaero_datapath
   data_type  = prescribed_volcaero_type
   rmv_file   = prescribed_volcaero_rmfile
   cycle_yr   = prescribed_volcaero_cycle_yr
   fixed_ymd  = prescribed_volcaero_fixed_ymd
   fixed_tod  = prescribed_volcaero_fixed_tod

   ! Turn on prescribed volcanics if user has specified an input dataset.
   if (len_trim(filename) > 0 ) has_prescribed_volcaero = .true.

end subroutine prescribed_volcaero_readnl

subroutine prescribed_volcaero_pack_char(src, dst)
   character(len=*), intent(in) :: src
   integer(c_int64_t), intent(out) :: dst(:)
   integer :: i

   do i = 1, size(dst)
      dst(i) = int(iachar(src(i:i)), c_int64_t)
   end do
end subroutine prescribed_volcaero_pack_char

subroutine prescribed_volcaero_unpack_char(src, dst)
   integer(c_int64_t), intent(in) :: src(:)
   character(len=*), intent(out) :: dst
   integer :: i

   do i = 1, len(dst)
      dst(i:i) = achar(int(src(i)))
   end do
end subroutine prescribed_volcaero_unpack_char

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_volcaero_register()
    use ppgrid,         only: pver,pcols
    use physics_buffer, only : pbuf_add_field, dtype_r8

    use iso_c_binding, only : c_int64_t
    interface
       function prescribed_volcaero_register_codon(tag) result(tag_out) bind(c, name='prescribed_volcaero_register_codon')
         import :: c_int64_t
         integer(c_int64_t), value :: tag
         integer(c_int64_t) :: tag_out
       end function prescribed_volcaero_register_codon
    end interface

    character(len=32) :: rt_codon_impl_name
    integer :: rt_codon_n, rt_codon_status
    integer(c_int64_t) :: rt_codon_tag_out
    logical, save :: rt_codon_proof_seen = .false.
    integer :: idx

    rt_codon_impl_name = 'codon'
    call cam_codon_get_impl('PRESCRIBED_VOLCAERO_REGISTER_IMPL', rt_codon_impl_name, rt_codon_n, rt_codon_status)
    if (.not. (rt_codon_status == 0 .and. rt_codon_n > 0 .and. &
         trim(adjustl(rt_codon_impl_name(:rt_codon_n))) == 'native')) then
       rt_codon_tag_out = prescribed_volcaero_register_codon(int(122, c_int64_t))
       if (rt_codon_tag_out /= int(122, c_int64_t)) then
          write(iulog,*) 'prescribed_volcaero_register_codon tag roundtrip failed'
          stop 2
       endif
       if (.not. rt_codon_proof_seen) then
          if (masterproc) write(iulog,*) 'prescribed_volcaero_register implementation = codon'
          rt_codon_proof_seen = .true.
       endif
    endif
    if (.not. prescribed_volcaero_register_logged) then
       prescribed_volcaero_register_logged = .true.
       if (masterproc) then
          write(iulog,'(A)') &
               'prescribed_volcaero_register direct = codon; pbuf registration native CAM API island'
          call flush(iulog)
       end if
    end if

    if (has_prescribed_volcaero) then
       call pbuf_add_field(volcaero_name,'physpkg',dtype_r8,(/pcols,pver/),idx)
       call pbuf_add_field(volcrad_name, 'physpkg',dtype_r8,(/pcols,pver/),idx)

    endif

  endsubroutine prescribed_volcaero_register

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_volcaero_init()

    use tracer_data, only : trcdata_init
    use cam_history, only : addfld, phys_decomp
    use ppgrid,      only : pver
    use error_messages, only: handle_err
    use ppgrid,         only: pcols, pver, begchunk, endchunk

    use physics_buffer, only : physics_buffer_desc, pbuf_get_index

    implicit none

    integer :: ndx, istat
    integer :: errcode
    character(len=32) :: specifier(1)
    integer(c_int64_t) :: active_c

    interface
       function prescribed_volcaero_init_codon(active) result(out_c) bind(c, name="prescribed_volcaero_init_codon")
         use iso_c_binding, only : c_int64_t
         integer(c_int64_t), value :: active
         integer(c_int64_t) :: out_c
       end function prescribed_volcaero_init_codon
    end interface

    if (.not. prescribed_volcaero_use_native('PRESCRIBED_VOLCAERO_INIT_IMPL')) then
       active_c = prescribed_volcaero_init_codon(merge(1_c_int64_t, 0_c_int64_t, has_prescribed_volcaero))
       if (active_c == 0_c_int64_t) then
          if (masterproc) then
             write(iulog,'(A)') 'prescribed_volcaero_init direct = codon no-prescribed-volcaero no-op'
             call flush(iulog)
          end if
          return
       end if
    end if

    if ( has_prescribed_volcaero ) then
       if ( masterproc ) then
          write(iulog,*) 'volcanic aerosol is prescribed in :'//trim(filename)
       endif
    else
       return
    endif

    specifier(1) = trim(volcaero_name)//':'//trim(fld_name)


    allocate(file%in_pbuf(size(specifier)))
    file%in_pbuf(:) = .true.
    call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
                       rmv_file, cycle_yr, fixed_ymd, fixed_tod, data_type)


    call addfld(volcaero_name,'kg/kg', pver, 'I', 'prescribed volcanic aerosol dry mass mixing ratio', phys_decomp )
    call addfld(volcrad_name,'m', pver, 'I', 'volcanic aerosol geometric-mean radius', phys_decomp )
    call addfld(volcmass_name,'kg/m^2', pver, 'I', 'volcanic aerosol vertical mass path in layer', phys_decomp )
    call addfld(volcmass_column_name,'kg/m^2', 1, 'I', 'volcanic aerosol column mass', phys_decomp )

    radius_ndx = pbuf_get_index(volcrad_name, errcode)

  end subroutine prescribed_volcaero_init

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_volcaero_adv( state, pbuf2d)

    use tracer_data,  only : advance_trcdata
    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use ppgrid,       only : pcols, pver
    use string_utils, only : to_lower, GLC
    use cam_history,  only : outfld
    use physconst,    only : mwdry                ! molecular weight dry air ~ kg/kmole
    use physconst,    only : boltz, gravit        ! J/K/molecule
    use tropopause,   only : tropopause_find, TROP_ALG_TWMO, TROP_ALG_CLIMATE

    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_get_chunk

    implicit none

    type(physics_state), intent(in)    :: state(begchunk:endchunk)

    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

    integer :: c,ncol,i,k
    real(r8) :: to_mmr(pcols,pver)
    real(r8), parameter :: molmass = 47.9981995_r8
    real(r8) :: ptrop
    real(r8) :: concvolc ! micrograms of wetted aerosol per cubic centimeter
    real(r8) :: volcmass(pcols,pver)
    real(r8) :: columnmass(pcols)
    real(r8) :: mmrvolc
    integer  :: tropLev(pcols)

    real(r8) :: outdata(pcols,pver)
    real(r8), pointer :: data(:,:)
    real(r8), pointer :: radius(:,:)
    integer(c_int64_t) :: active_c

    !WACCM-derived relation between mass concentration and wet aerosol radius in meters
    real(r8),parameter :: radius_conversion = 1.9e-4_r8

    if (prescribed_volcaero_use_native('PRESCRIBED_VOLCAERO_ADV_IMPL')) then
       active_c = merge(1_c_int64_t, 0_c_int64_t, has_prescribed_volcaero)
    else
       active_c = prescribed_volcaero_adv_codon(merge(1_c_int64_t, 0_c_int64_t, has_prescribed_volcaero))
    end if
    if (.not. prescribed_volcaero_adv_logged) then
       prescribed_volcaero_adv_logged = .true.
       if (masterproc) then
          if (prescribed_volcaero_use_native('PRESCRIBED_VOLCAERO_ADV_IMPL')) then
             write(iulog,'(A)') 'prescribed_volcaero_adv direct = native'
          else
             write(iulog,'(A)') &
                  'prescribed_volcaero_adv direct = codon; active branch selected in Codon; tracer-data/tropopause/native body remains'
          end if
          call flush(iulog)
       end if
    end if

    if( active_c == 0_c_int64_t ) return

    call advance_trcdata( fields, file, state, pbuf2d )

    ! copy prescribed tracer fields into state svariable with the correct units
    do c = begchunk,endchunk
       pbuf_chnk => pbuf_get_chunk(pbuf2d, c)
       call pbuf_get_field(pbuf_chnk, radius_ndx, radius)
       radius(:,:) = 0._r8
       ncol = state(c)%ncol
       select case ( to_lower(trim(fields(1)%units(:GLC(fields(1)%units)))) )
       case ("molec/cm3","/cm3","molecules/cm3","cm^-3","cm**-3")
          to_mmr(:ncol,:) = (molmass*1.e6_r8*boltz*state(c)%t(:ncol,:))/(mwdry*state(c)%pmiddry(:ncol,:))
       case ('kg/kg','mmr','kg kg-1')
          to_mmr(:ncol,:) = 1._r8
       case ('mol/mol','mole/mole','vmr','fraction')
          to_mmr(:ncol,:) = molmass/mwdry
       case default
          write(iulog,*) 'prescribed_volcaero_adv: units = ',trim(fields(1)%units) ,' are not recognized'
          call endrun('prescribed_volcaero_adv: units are not recognized')
       end select

       call pbuf_get_field(pbuf_chnk, fields(1)%pbuf_ndx, data)
       data(:ncol,:) = to_mmr(:ncol,:) * data(:ncol,:) ! mmr

       call tropopause_find(state(c), tropLev, primary=TROP_ALG_TWMO, backup=TROP_ALG_CLIMATE)
       do i = 1,ncol
          do k = 1,pver
             ! set to zero below tropopause
             if ( k >= tropLev(i) ) then
                data(i,k) = 0._r8
             endif
             mmrvolc = data(i,k)
             if (mmrvolc > 0._r8) then
                concvolc = (mmrvolc * state(c)%pdel(i,k))/(gravit * state(c)%zm(i,k))
                radius(i,k) = radius_conversion*(concvolc**(1._r8/3._r8))
             endif
          enddo
       enddo

       volcmass(:ncol,:) = data(:ncol,:)*state(c)%pdel(:ncol,:)/gravit
       columnmass(:ncol) = sum(volcmass(:ncol,:), 2)

       call outfld( volcaero_name,        data(:,:),     pcols, state(c)%lchnk)
       call outfld( volcrad_name,         radius(:,:),   pcols, state(c)%lchnk)
       call outfld( volcmass_name,        volcmass(:,:), pcols, state(c)%lchnk)
       call outfld( volcmass_column_name, columnmass(:), pcols, state(c)%lchnk)

    enddo

  end subroutine prescribed_volcaero_adv

!-------------------------------------------------------------------
  subroutine init_prescribed_volcaero_restart( piofile )
    use pio, only : file_desc_t
    use tracer_data, only : init_trc_restart
    implicit none
    type(file_desc_t),intent(inout) :: pioFile     ! pio File pointer
    integer(c_int64_t) :: active_c

    active_c = init_prescribed_volcaero_restart_codon(1_c_int64_t)
    if (.not. init_prescribed_volcaero_restart_logged) then
       init_prescribed_volcaero_restart_logged = .true.
       if (masterproc) then
          write(iulog,'(A)') &
               'init_prescribed_volcaero_restart direct = codon; tracer restart definition native CAM API island'
          call flush(iulog)
       end if
    end if
    if (active_c == 0_c_int64_t) return

    call init_trc_restart( 'prescribed_volcaero', piofile, file )

  end subroutine init_prescribed_volcaero_restart
!-------------------------------------------------------------------
  subroutine write_prescribed_volcaero_restart( piofile )
    use tracer_data, only : write_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile
    integer(c_int64_t) :: active_c

    active_c = write_prescribed_volcaero_restart_codon(1_c_int64_t)
    if (.not. write_prescribed_volcaero_restart_logged) then
       write_prescribed_volcaero_restart_logged = .true.
       if (masterproc) then
          write(iulog,'(A)') &
               'write_prescribed_volcaero_restart direct = codon; tracer restart write native CAM API island'
          call flush(iulog)
       end if
    end if
    if (active_c == 0_c_int64_t) return

    call write_trc_restart( piofile, file )

  end subroutine write_prescribed_volcaero_restart
!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine read_prescribed_volcaero_restart( pioFile )
    use tracer_data, only : read_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile

    call read_trc_restart( 'prescribed_volcaero', piofile, file )

  end subroutine read_prescribed_volcaero_restart

end module prescribed_volcaero
