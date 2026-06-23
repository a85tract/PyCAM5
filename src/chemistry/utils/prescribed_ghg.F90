!-------------------------------------------------------------------
! manages reading and interpolation of prescribed ghg tracers
! Created by: Francis Vitt
!-------------------------------------------------------------------
module prescribed_ghg

  use shr_kind_mod,     only : r8 => shr_kind_r8
  use iso_c_binding,    only : c_int64_t
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

  public :: prescribed_ghg_init
  public :: prescribed_ghg_adv
  public :: write_prescribed_ghg_restart
  public :: read_prescribed_ghg_restart
  public :: has_prescribed_ghg
  public :: prescribed_ghg_register
  public :: init_prescribed_ghg_restart
  public :: prescribed_ghg_readnl

  logical :: has_prescribed_ghg = .false.
  integer, parameter, public :: N_GHG = 5
  integer :: number_flds

  character(len=256) :: filename = ''
  character(len=256) :: filelist = ''
  character(len=256) :: datapath = ''
  character(len=32)  :: datatype = 'SERIAL'
  logical            :: rmv_file = .false.
  integer            :: cycle_yr  = 0
  integer            :: fixed_ymd = 0
  integer            :: fixed_tod = 0
  character(len=16)  :: specifier(N_GHG) = ''

  character(len=8)    :: ghg_names(N_GHG) = (/ 'prsd_co2',  'prsd_ch4',  'prsd_n2o',  'prsd_f11',  'prsd_f12'  /)
  real(r8), parameter :: molmass(N_GHG)   = (/ 44.00980_r8, 16.04060_r8, 44.01288_r8, 137.3675_r8, 120.9132_r8 /)

  integer :: index_map(N_GHG)
  logical :: prescribed_ghg_register_logged = .false.
  logical :: prescribed_ghg_adv_logged = .false.
  logical :: init_prescribed_ghg_restart_logged = .false.
  logical :: write_prescribed_ghg_restart_logged = .false.

  interface
     function prescribed_ghg_register_codon(active_c) result(out_c) bind(c, name="prescribed_ghg_register_codon")
       import :: c_int64_t
       integer(c_int64_t), value :: active_c
       integer(c_int64_t) :: out_c
     end function prescribed_ghg_register_codon
     function prescribed_ghg_adv_codon(active_c) result(out_c) bind(c, name="prescribed_ghg_adv_codon")
       import :: c_int64_t
       integer(c_int64_t), value :: active_c
       integer(c_int64_t) :: out_c
     end function prescribed_ghg_adv_codon
     function init_prescribed_ghg_restart_codon(stage_c) result(out_c) bind(c, name="init_prescribed_ghg_restart_codon")
       import :: c_int64_t
       integer(c_int64_t), value :: stage_c
       integer(c_int64_t) :: out_c
     end function init_prescribed_ghg_restart_codon
     function write_prescribed_ghg_restart_codon(stage_c) result(out_c) bind(c, name="write_prescribed_ghg_restart_codon")
       import :: c_int64_t
       integer(c_int64_t), value :: stage_c
       integer(c_int64_t) :: out_c
     end function write_prescribed_ghg_restart_codon
  end interface

contains

logical function prescribed_ghg_use_native(selector)
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
     prescribed_ghg_use_native = trim(adjustl(impl_name(:n))) == 'native'
  else
     prescribed_ghg_use_native = .false.
  end if
end function prescribed_ghg_use_native

!-------------------------------------------------------------------
!-------------------------------------------------------------------
subroutine prescribed_ghg_readnl(nlfile)

   use iso_c_binding, only : c_int64_t
   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

    interface
       function prescribed_ghg_readnl_codon(tag) result(tag_out) bind(c, name='prescribed_ghg_readnl_codon')
         import :: c_int64_t
         integer(c_int64_t), value :: tag
         integer(c_int64_t) :: tag_out
       end function prescribed_ghg_readnl_codon
    end interface

    character(len=32) :: rt_codon_impl_name
    integer :: rt_codon_n, rt_codon_status
    integer(c_int64_t) :: rt_codon_tag_out
    logical, save :: rt_codon_proof_seen = .false.


   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'prescribed_ghg_readnl'

   character(len=16)  :: prescribed_ghg_specifier(N_GHG)
   character(len=256) :: prescribed_ghg_file
   character(len=256) :: prescribed_ghg_filelist
   character(len=256) :: prescribed_ghg_datapath
   character(len=32)  :: prescribed_ghg_type
   logical            :: prescribed_ghg_rmfile
   integer            :: prescribed_ghg_cycle_yr
   integer            :: prescribed_ghg_fixed_ymd
   integer            :: prescribed_ghg_fixed_tod

   namelist /prescribed_ghg_nl/ &
      prescribed_ghg_specifier, &
      prescribed_ghg_file,      &
      prescribed_ghg_filelist,  &
      prescribed_ghg_datapath,  &
      prescribed_ghg_type,      &
      prescribed_ghg_rmfile,    &
      prescribed_ghg_cycle_yr,  &
      prescribed_ghg_fixed_ymd, &
      prescribed_ghg_fixed_tod      
   !-----------------------------------------------------------------------------

    rt_codon_impl_name = 'codon'
    call cam_codon_get_impl('PRESCRIBED_GHG_READNL_IMPL', rt_codon_impl_name, rt_codon_n, rt_codon_status)
    if (.not. (rt_codon_status == 0 .and. rt_codon_n > 0 .and. &
         trim(adjustl(rt_codon_impl_name(:rt_codon_n))) == 'native')) then
       rt_codon_tag_out = prescribed_ghg_readnl_codon(int(117, c_int64_t))
       if (rt_codon_tag_out /= int(117, c_int64_t)) then
          write(iulog,*) 'prescribed_ghg_readnl_codon tag roundtrip failed'
          stop 2
       endif
       if (.not. rt_codon_proof_seen) then
          write(iulog,*) 'prescribed_ghg_readnl implementation = codon'
          rt_codon_proof_seen = .true.
       endif
    endif

   ! Initialize namelist variables from local module variables.
   prescribed_ghg_specifier= specifier
   prescribed_ghg_file     = filename
   prescribed_ghg_filelist = filelist
   prescribed_ghg_datapath = datapath
   prescribed_ghg_type     = datatype
   prescribed_ghg_rmfile   = rmv_file
   prescribed_ghg_cycle_yr = cycle_yr
   prescribed_ghg_fixed_ymd= fixed_ymd
   prescribed_ghg_fixed_tod= fixed_tod

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'prescribed_ghg_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, prescribed_ghg_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(prescribed_ghg_specifier,len(prescribed_ghg_specifier(1))*N_GHG,     mpichar, 0, mpicom)
   call mpibcast(prescribed_ghg_file,     len(prescribed_ghg_file),     mpichar, 0, mpicom)
   call mpibcast(prescribed_ghg_filelist, len(prescribed_ghg_filelist), mpichar, 0, mpicom)
   call mpibcast(prescribed_ghg_datapath, len(prescribed_ghg_datapath), mpichar, 0, mpicom)
   call mpibcast(prescribed_ghg_type,     len(prescribed_ghg_type),     mpichar, 0, mpicom)
   call mpibcast(prescribed_ghg_rmfile,   1, mpilog,  0, mpicom)
   call mpibcast(prescribed_ghg_cycle_yr, 1, mpiint,  0, mpicom)
   call mpibcast(prescribed_ghg_fixed_ymd,1, mpiint,  0, mpicom)
   call mpibcast(prescribed_ghg_fixed_tod,1, mpiint,  0, mpicom)
#endif

   ! Update module variables with user settings.
   specifier  = prescribed_ghg_specifier
   filename   = prescribed_ghg_file
   filelist   = prescribed_ghg_filelist
   datapath   = prescribed_ghg_datapath
   datatype   = prescribed_ghg_type
   rmv_file   = prescribed_ghg_rmfile
   cycle_yr   = prescribed_ghg_cycle_yr
   fixed_ymd  = prescribed_ghg_fixed_ymd
   fixed_tod  = prescribed_ghg_fixed_tod

   ! Turn on prescribed volcanics if user has specified an input dataset.
   if (len_trim(filename) > 0 ) has_prescribed_ghg = .true.

end subroutine prescribed_ghg_readnl

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_ghg_register()
    use ppgrid,         only: pver, pcols
    use physics_buffer, only : pbuf_add_field, dtype_r8

    integer :: i,idx
    integer(c_int64_t) :: active_c

    active_c = prescribed_ghg_register_codon(merge(1_c_int64_t, 0_c_int64_t, has_prescribed_ghg))
    if (.not. prescribed_ghg_register_logged) then
       prescribed_ghg_register_logged = .true.
       if (masterproc) then
          write(iulog,'(A)') &
               'prescribed_ghg_register direct = codon; pbuf registration native CAM API island'
          call flush(iulog)
       end if
    end if

    if (active_c /= 0_c_int64_t) then
       do i = 1,N_GHG
          call pbuf_add_field(ghg_names(i),'physpkg',dtype_r8,(/pcols,pver/),idx)
       enddo
    endif

  endsubroutine prescribed_ghg_register
!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_ghg_init()

    use tracer_data, only : trcdata_init
    use cam_history, only : addfld, phys_decomp
    use ppgrid,      only : pver
    use error_messages, only: handle_err
    use ppgrid,         only: pcols, pver, begchunk, endchunk
    use physics_buffer, only : physics_buffer_desc

    implicit none

    integer :: ndx, istat, i
    integer(c_int64_t) :: active_c

    interface
       function prescribed_ghg_init_codon(active) result(out_c) bind(c, name="prescribed_ghg_init_codon")
         use iso_c_binding, only : c_int64_t
         integer(c_int64_t), value :: active
         integer(c_int64_t) :: out_c
       end function prescribed_ghg_init_codon
    end interface

    if (.not. prescribed_ghg_use_native('PRESCRIBED_GHG_INIT_IMPL')) then
       active_c = prescribed_ghg_init_codon(merge(1_c_int64_t, 0_c_int64_t, has_prescribed_ghg))
       if (active_c == 0_c_int64_t) then
          if (masterproc) then
             write(iulog,'(A)') 'prescribed_ghg_init direct = codon no-prescribed-ghg no-op'
             call flush(iulog)
          end if
          return
       end if
    end if
    
    if ( has_prescribed_ghg ) then
       if ( masterproc ) then
          write(iulog,*) 'ghg is prescribed in :'//trim(filename)
       endif
    else
       return
    endif

    allocate(file%in_pbuf(size(specifier)))
    file%in_pbuf(:) = .true.
    call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
                       rmv_file, cycle_yr, fixed_ymd, fixed_tod, datatype)
        
    number_flds = 0
    if (associated(fields)) number_flds = size( fields )

    if( number_flds < 1 ) then
       if ( masterproc ) then
          write(iulog,*) 'There are no prescribed ghg tracers'
          write(iulog,*) ' '
       endif
       return
    end if

    do i = 1,number_flds
       ndx = get_ndx( fields(i)%fldnam )
       index_map(i) = ndx

       if (ndx < 1) then
          call endrun('prescribed_ghg_init: '//trim(fields(i)%fldnam)//' is not one of the named ghg fields in pbuf2d')
       endif
       call addfld( fields(i)%fldnam,'kg/kg', pver, 'I', 'prescribed ghg', phys_decomp )
    enddo

  end subroutine prescribed_ghg_init

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_ghg_adv( state, pbuf2d)

    use tracer_data,  only : advance_trcdata
    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use ppgrid,       only : pcols, pver
    use string_utils, only : to_lower, GLC
    use cam_history,  only : outfld
    use physconst,    only : mwdry                ! molecular weight dry air ~ kg/kmole
    use physconst,    only : boltz                ! J/K/molecule
    
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_set_field, pbuf_get_chunk

    implicit none

    type(physics_state), intent(in)    :: state(begchunk:endchunk)                 
    
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    type(physics_buffer_desc), pointer :: pbuf_chnk(:)
    integer :: ind,c,ncol,i
    real(r8) :: to_mmr(pcols,pver)
    real(r8) :: outdata(pcols,pver)
    real(r8),pointer :: tmpptr(:,:)

    character(len=32) :: units_str
    integer(c_int64_t) :: active_c

    if (prescribed_ghg_use_native('PRESCRIBED_GHG_ADV_IMPL')) then
       active_c = merge(1_c_int64_t, 0_c_int64_t, has_prescribed_ghg)
    else
       active_c = prescribed_ghg_adv_codon(merge(1_c_int64_t, 0_c_int64_t, has_prescribed_ghg))
    end if
    if (.not. prescribed_ghg_adv_logged) then
       prescribed_ghg_adv_logged = .true.
       if (masterproc) then
          if (prescribed_ghg_use_native('PRESCRIBED_GHG_ADV_IMPL')) then
             write(iulog,'(A)') 'prescribed_ghg_adv direct = native'
          else
             write(iulog,'(A)') &
                  'prescribed_ghg_adv direct = codon; active branch selected in Codon; tracer-data/unit conversion native body remains'
          end if
          call flush(iulog)
       end if
    end if

    if( active_c == 0_c_int64_t ) return

    call advance_trcdata( fields, file, state, pbuf2d )
    
    ! set the correct units and invoke history outfld
    do i = 1,number_flds
       ind = index_map(i)

       units_str = trim(to_lower(trim(fields(i)%units(:GLC(fields(i)%units)))))

!$OMP PARALLEL DO PRIVATE (C, NCOL, OUTDATA, TO_MMR, tmpptr, pbuf_chnk)
       do c = begchunk,endchunk
          ncol = state(c)%ncol

          select case ( units_str )
          case ("molec/cm3","/cm3","molecules/cm3","cm^-3","cm**-3")
             to_mmr(:ncol,:) = (molmass(ind)*1.e6_r8*boltz*state(c)%t(:ncol,:))/(mwdry*state(c)%pmiddry(:ncol,:))
          case ('kg/kg','mmr')
             to_mmr(:ncol,:) = 1._r8
          case ('mol/mol','mole/mole','vmr','fraction')
             to_mmr(:ncol,:) = molmass(ind)/mwdry
          case default
             print*, 'prescribed_ghg_adv: units = ',trim(fields(i)%units) ,' are not recognized'
             call endrun('prescribed_ghg_adv: units are not recognized')
          end select

          pbuf_chnk => pbuf_get_chunk(pbuf2d, c)
          call pbuf_get_field(pbuf_chnk, fields(i)%pbuf_ndx, tmpptr )

          tmpptr(:ncol,:) = tmpptr(:ncol,:)*to_mmr(:ncol,:)

          outdata(:ncol,:) = tmpptr(:ncol,:) 
          call outfld( fields(1)%fldnam, outdata(:ncol,:), ncol, state(c)%lchnk )

       enddo
    enddo

  end subroutine prescribed_ghg_adv

!-------------------------------------------------------------------

!-------------------------------------------------------------------
  subroutine init_prescribed_ghg_restart( piofile )
    use pio, only : file_desc_t
    use tracer_data, only : init_trc_restart
    implicit none
    type(file_desc_t),intent(inout) :: pioFile     ! pio File pointer
    integer(c_int64_t) :: active_c

    active_c = init_prescribed_ghg_restart_codon(1_c_int64_t)
    if (.not. init_prescribed_ghg_restart_logged) then
       init_prescribed_ghg_restart_logged = .true.
       if (masterproc) then
          write(iulog,'(A)') &
               'init_prescribed_ghg_restart direct = codon; tracer restart definition native CAM API island'
          call flush(iulog)
       end if
    end if
    if (active_c == 0_c_int64_t) return

    call init_trc_restart( 'prescribed_ghg', piofile, file )

  end subroutine init_prescribed_ghg_restart
!-------------------------------------------------------------------
  subroutine write_prescribed_ghg_restart( piofile )
    use tracer_data, only : write_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile
    integer(c_int64_t) :: active_c

    active_c = write_prescribed_ghg_restart_codon(1_c_int64_t)
    if (.not. write_prescribed_ghg_restart_logged) then
       write_prescribed_ghg_restart_logged = .true.
       if (masterproc) then
          write(iulog,'(A)') &
               'write_prescribed_ghg_restart direct = codon; tracer restart write native CAM API island'
          call flush(iulog)
       end if
    end if
    if (active_c == 0_c_int64_t) return

    call write_trc_restart( piofile, file )

  end subroutine write_prescribed_ghg_restart

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine read_prescribed_ghg_restart( pioFile )
    use tracer_data, only : read_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile

    call read_trc_restart( 'prescribed_ghg', piofile, file )

  end subroutine read_prescribed_ghg_restart
!-------------------------------------------------------------------
  integer function get_ndx( name )

    implicit none
    character(len=*), intent(in) :: name

    integer :: i

    get_ndx = 0
    do i = 1,N_GHG
      if ( trim(name) == trim(ghg_names(i)) ) then
        get_ndx = i
        return
      endif
    enddo

  end function get_ndx

end module prescribed_ghg
