module gcr_ionization

  use shr_kind_mod,   only : r8 => shr_kind_r8
  use cam_abortutils, only : endrun
  use spmd_utils,     only : masterproc
  use tracer_data,    only : trfld,trfile
  use cam_logfile,    only : iulog
  use physics_buffer, only : physics_buffer_desc
  use cam_history,    only : outfld, addfld, phys_decomp
  use physics_types,  only : physics_state
  use ppgrid,         only : begchunk, endchunk
  use ppgrid,         only : pcols, pver
  use tracer_data,    only : trcdata_init, advance_trcdata
  use mo_util,        only : chemistry_misc_codon_touch
  use iso_c_binding,  only : c_int64_t

  implicit none
  private 
  public :: gcr_ionization_readnl
  public :: gcr_ionization_init
  public :: gcr_ionization_adv
  public :: gcr_ionization_noxhox

  type(trfld), pointer :: fields(:)
  type(trfile), save :: file

  character(len=32)  :: specifier(1) = 'prod'
  character(len=256) :: filename = ''
  character(len=256) :: filelist = ''
  character(len=256) :: datapath = ''
  character(len=32)  :: datatype = 'SERIAL'
  logical            :: rmv_file = .false.
  integer            :: cycle_yr  = 0
  integer            :: fixed_ymd = 0
  integer            :: fixed_tod = 0

  logical :: has_gcr_ionization = .false.
  logical :: gcr_ionization_init_use_native_impl = .false.
  logical :: gcr_ionization_init_impl_selected = .false.
  logical :: gcr_ionization_init_proof_written = .false.
  logical :: gcr_ionization_adv_use_native_impl = .false.
  logical :: gcr_ionization_adv_impl_selected = .false.
  logical :: gcr_ionization_adv_proof_written = .false.

  interface
    function gcr_ionization_init_codon() result(out_c) bind(c, name="gcr_ionization_init_codon")
      use iso_c_binding, only : c_int64_t
      integer(c_int64_t) :: out_c
    end function gcr_ionization_init_codon

    function gcr_ionization_adv_codon() result(out_c) bind(c, name="gcr_ionization_adv_codon")
      use iso_c_binding, only : c_int64_t
      integer(c_int64_t) :: out_c
    end function gcr_ionization_adv_codon
  end interface

contains
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  subroutine gcr_ionization_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'gcr_ionization_readnl'

    character(len=16)  ::  gcr_ionization_fldname
    character(len=256) ::  gcr_ionization_filename
    character(len=256) ::  gcr_ionization_datapath
    character(len=256) ::  gcr_ionization_filelist
    character(len=32)  ::  gcr_ionization_datatype
    integer            ::  gcr_ionization_cycle_yr
    integer            ::  gcr_ionization_fixed_ymd
    integer            ::  gcr_ionization_fixed_tod

    namelist /gcr_ionization_nl/ &
         gcr_ionization_fldname, &
         gcr_ionization_filename, &
         gcr_ionization_datapath, &
         gcr_ionization_filelist, &
         gcr_ionization_datatype, &
         gcr_ionization_cycle_yr, &
         gcr_ionization_fixed_ymd, &
         gcr_ionization_fixed_tod

    gcr_ionization_fldname = specifier(1)
    gcr_ionization_filename = filename
    gcr_ionization_datapath = datapath
    gcr_ionization_filelist = filelist
    gcr_ionization_datatype = datatype
    gcr_ionization_cycle_yr = cycle_yr
    gcr_ionization_fixed_ymd = fixed_ymd
    gcr_ionization_fixed_tod = fixed_tod

    ! Read namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'gcr_ionization_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, gcr_ionization_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(gcr_ionization_fldname,  len(gcr_ionization_fldname),  mpichar, 0, mpicom)
    call mpibcast(gcr_ionization_filename, len(gcr_ionization_filename), mpichar, 0, mpicom)
    call mpibcast(gcr_ionization_filelist, len(gcr_ionization_filelist), mpichar, 0, mpicom)
    call mpibcast(gcr_ionization_datapath, len(gcr_ionization_datapath), mpichar, 0, mpicom)
    call mpibcast(gcr_ionization_datatype, len(gcr_ionization_datatype), mpichar, 0, mpicom)
    call mpibcast(gcr_ionization_cycle_yr, 1, mpiint,  0, mpicom)
    call mpibcast(gcr_ionization_fixed_ymd,1, mpiint,  0, mpicom)
    call mpibcast(gcr_ionization_fixed_tod,1, mpiint,  0, mpicom)
#endif

    ! Update module variables with user settings.
    specifier(1) = gcr_ionization_fldname
    filename  = gcr_ionization_filename
    filelist  = gcr_ionization_filelist
    datapath  = gcr_ionization_datapath
    datatype  = gcr_ionization_datatype
    cycle_yr  = gcr_ionization_cycle_yr
    fixed_ymd = gcr_ionization_fixed_ymd
    fixed_tod = gcr_ionization_fixed_tod

    ! Turn on galactic cosmic rays if user has specified an input dataset.
    if (len_trim(filename) > 0 ) has_gcr_ionization = .true.
    call chemistry_misc_codon_touch('gcr_ionization', 124)

  end subroutine gcr_ionization_readnl

  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  subroutine gcr_ionization_init()

    integer(c_int64_t) :: out_c

    call gcr_ionization_init_select_impl()

    if (.not.gcr_ionization_init_use_native_impl) then
       out_c = gcr_ionization_init_codon()
       call gcr_ionization_init_proof_once()
    end if

    if (.not.has_gcr_ionization) return
    
    allocate(file%in_pbuf(size(specifier)))
    file%in_pbuf(:) = .false.
    call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
                       rmv_file, cycle_yr, fixed_ymd, fixed_tod, datatype )

    call addfld('GCRION','/cm3/sec', pver, 'I', 'GCR ionizaton rate', phys_decomp )

  end subroutine gcr_ionization_init

  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  subroutine gcr_ionization_adv( pbuf2d, state )
    type(physics_state), intent(in):: state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    integer(c_int64_t) :: out_c

    call gcr_ionization_adv_select_impl()

    if (.not.gcr_ionization_adv_use_native_impl) then
       out_c = gcr_ionization_adv_codon()
       call gcr_ionization_adv_proof_once()
    end if

    if (.not.has_gcr_ionization) return

    call advance_trcdata( fields, file, state, pbuf2d )

  end subroutine gcr_ionization_adv

  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  subroutine gcr_ionization_init_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (gcr_ionization_init_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('GCR_IONIZATION_INIT_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       gcr_ionization_init_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       gcr_ionization_init_use_native_impl = .false.
    end if

    gcr_ionization_init_impl_selected = .true.

    if (masterproc) then
       if (gcr_ionization_init_use_native_impl) then
          write(iulog,*) 'gcr_ionization_init implementation = native'
       else
          write(iulog,*) 'gcr_ionization_init implementation = codon'
       end if
       call flush(iulog)
    end if

  end subroutine gcr_ionization_init_select_impl

  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  subroutine gcr_ionization_init_proof_once()

    if (gcr_ionization_init_proof_written) return
    gcr_ionization_init_proof_written = .true.

    if (masterproc) then
       write(iulog,'(A)') 'gcr_ionization_init entered (availability gate = codon; active data setup = native when enabled)'
       call flush(iulog)
    end if

  end subroutine gcr_ionization_init_proof_once

  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  subroutine gcr_ionization_adv_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (gcr_ionization_adv_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('GCR_IONIZATION_ADV_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       gcr_ionization_adv_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       gcr_ionization_adv_use_native_impl = .false.
    end if

    gcr_ionization_adv_impl_selected = .true.

    if (masterproc) then
       if (gcr_ionization_adv_use_native_impl) then
          write(iulog,*) 'gcr_ionization_adv implementation = native'
       else
          write(iulog,*) 'gcr_ionization_adv implementation = codon'
       end if
       call flush(iulog)
    end if

  end subroutine gcr_ionization_adv_select_impl

  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  subroutine gcr_ionization_adv_proof_once()

    if (gcr_ionization_adv_proof_written) return
    gcr_ionization_adv_proof_written = .true.

    if (masterproc) then
       write(iulog,'(A)') 'gcr_ionization_adv entered (availability gate = codon; active tracer advance = native when enabled)'
       call flush(iulog)
    end if

  end subroutine gcr_ionization_adv_proof_once

  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  subroutine gcr_ionization_get( ncol, lchnk, ionpairs )

    integer, intent(in) :: lchnk
    integer, intent(in) :: ncol
    real(r8), intent(out) :: ionpairs(:,:)

    ionpairs(:,:) = 0._r8

    if (.not.has_gcr_ionization) return

    ionpairs(:ncol,:) = fields(1)%data(:ncol,:,lchnk)
    call outfld( 'GCRION', ionpairs(:ncol,:), ncol, lchnk )

  end subroutine gcr_ionization_get

  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  subroutine gcr_ionization_noxhox( ncol, lchnk, zmid, gcr_nox, gcr_hox )
    use spehox,  only : hox_prod_factor

    integer, intent(in) :: ncol, lchnk
    real(r8), intent(in) :: zmid(:,:)

    real(r8), intent(out) :: gcr_nox(:,:)
    real(r8), intent(out) :: gcr_hox(:,:)

    real(r8) :: hoxprod_factor(pver)
    real(r8) :: ionpairs(pcols,pver)
    integer :: i

    gcr_nox(:,:) = 0._r8
    gcr_hox(:,:) = 0._r8

    if (.not.has_gcr_ionization) return
    call  gcr_ionization_get( ncol, lchnk, ionpairs )

    gcr_nox(:ncol,:) = ionpairs(:ncol,:)

    do i = 1,ncol
       hoxprod_factor(:pver) = hox_prod_factor( ionpairs(i,:pver), zmid(i,:pver) )
       gcr_hox(i,:pver) = hoxprod_factor(:pver) * ionpairs(i,:pver)
    end do

  end subroutine gcr_ionization_noxhox

end module gcr_ionization
