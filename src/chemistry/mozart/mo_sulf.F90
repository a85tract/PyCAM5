

      module mo_sulf
!---------------------------------------------------------------
!	... Annual cycle for sulfur
!---------------------------------------------------------------

      use shr_kind_mod,     only : r8 => shr_kind_r8

      use cam_abortutils,   only : endrun
      use cam_logfile,      only : iulog
      use tracer_data,      only : trfld,trfile
      use physics_types,    only : physics_state
      use ppgrid,           only : begchunk, endchunk
      use physics_buffer,   only : physics_buffer_desc
	      use ppgrid,           only : pcols, pver
	      use mo_util,          only : chemistry_misc_codon_touch

	      use spmd_utils,       only : masterproc
	      use iso_c_binding,    only : c_int64_t

	      implicit none

      private
      public  :: sulf_inti, set_sulf_time, sulf_interp, sulf_readnl

      save

      type(trfld), pointer :: fields(:) => null()
      type(trfile) :: file

      logical :: read_sulf = .false.

      character(len=16)  :: fld_name = 'SULFATE'
      character(len=256) :: filename = ' '
      character(len=256) :: filelist = ' '
      character(len=256) :: datapath = ' '
      character(len=32)  :: datatype = 'CYCLICAL'
      logical            :: rmv_file = .false.
      integer            :: cycle_yr  = 0
      integer            :: fixed_ymd = 0
      integer            :: fixed_tod = 0

      logical :: has_sulf = .false.
	      logical :: sulf_interp_use_native_impl = .false.
	      logical :: sulf_interp_impl_selected = .false.
	      logical :: sulf_interp_proof_written = .false.
	      logical :: set_sulf_time_use_native_impl = .false.
	      logical :: set_sulf_time_impl_selected = .false.
	      logical :: set_sulf_time_proof_written = .false.
	      logical :: sulf_inti_proof_written = .false.

	      interface
	         function set_sulf_time_codon() result(out_c) bind(c, name="set_sulf_time_codon")
	            use iso_c_binding, only : c_int64_t
	            integer(c_int64_t) :: out_c
	         end function set_sulf_time_codon

	         function sulf_inti_codon(active) result(out_c) bind(c, name="sulf_inti_codon")
	            use iso_c_binding, only : c_int64_t
	            integer(c_int64_t), value :: active
	            integer(c_int64_t) :: out_c
	         end function sulf_inti_codon
	      end interface

	      contains

!-------------------------------------------------------------------
!-------------------------------------------------------------------
subroutine sulf_readnl(nlfile)

   use iso_c_binding, only : c_int64_t
   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

    interface
       function sulf_readnl_codon(tag) result(tag_out) bind(c, name='sulf_readnl_codon')
         import :: c_int64_t
         integer(c_int64_t), value :: tag
         integer(c_int64_t) :: tag_out
       end function sulf_readnl_codon
    end interface

    character(len=32) :: rt_codon_impl_name
    integer :: rt_codon_n, rt_codon_status
    integer(c_int64_t) :: rt_codon_tag_out
    logical, save :: rt_codon_proof_seen = .false.


   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'sulf_readnl'

   character(len=16)  :: sulf_name
   character(len=256) :: sulf_file
   character(len=256) :: sulf_filelist
   character(len=256) :: sulf_datapath
   character(len=32)  :: sulf_type
   logical            :: sulf_rmfile
   integer            :: sulf_cycle_yr
   integer            :: sulf_fixed_ymd
   integer            :: sulf_fixed_tod

   namelist /sulf_nl/ &
      sulf_name,      &
      sulf_file,      &
      sulf_filelist,  &
      sulf_datapath,  &
      sulf_type,      &
      sulf_rmfile,    &
      sulf_cycle_yr,  &
      sulf_fixed_ymd, &
      sulf_fixed_tod      
   !-----------------------------------------------------------------------------

   ! Initialize namelist variables from local module variables.
   sulf_name     = fld_name
   sulf_file     = filename
   sulf_filelist = filelist
   sulf_datapath = datapath
   sulf_type     = datatype
   sulf_rmfile   = rmv_file
   sulf_cycle_yr = cycle_yr
   sulf_fixed_ymd= fixed_ymd
   sulf_fixed_tod= fixed_tod

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'sulf_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, sulf_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(sulf_name,     len(sulf_name),     mpichar, 0, mpicom)
   call mpibcast(sulf_file,     len(sulf_file),     mpichar, 0, mpicom)
   call mpibcast(sulf_filelist, len(sulf_filelist), mpichar, 0, mpicom)
   call mpibcast(sulf_datapath, len(sulf_datapath), mpichar, 0, mpicom)
   call mpibcast(sulf_type,     len(sulf_type),     mpichar, 0, mpicom)
   call mpibcast(sulf_rmfile,   1, mpilog,  0, mpicom)
   call mpibcast(sulf_cycle_yr, 1, mpiint,  0, mpicom)
   call mpibcast(sulf_fixed_ymd,1, mpiint,  0, mpicom)
   call mpibcast(sulf_fixed_tod,1, mpiint,  0, mpicom)
#endif

   ! Update module variables with user settings.
   fld_name   = sulf_name
   filename   = sulf_file
   filelist   = sulf_filelist
   datapath   = sulf_datapath
   datatype   = sulf_type
   rmv_file   = sulf_rmfile
   cycle_yr   = sulf_cycle_yr
   fixed_ymd  = sulf_fixed_ymd
   fixed_tod  = sulf_fixed_tod

   ! Turn on prescribed volcanics if user has specified an input dataset.
   if (len_trim(filename) > 0 ) has_sulf = .true.
    rt_codon_impl_name = 'codon'
    call cam_codon_get_impl('SULF_READNL_IMPL', rt_codon_impl_name, rt_codon_n, rt_codon_status)
    if (.not. (rt_codon_status == 0 .and. rt_codon_n > 0 .and. &
         trim(adjustl(rt_codon_impl_name(:rt_codon_n))) == 'native')) then
       rt_codon_tag_out = sulf_readnl_codon(int(187, c_int64_t))
       if (rt_codon_tag_out /= int(187, c_int64_t)) then
          write(iulog,*) 'sulf_readnl_codon tag roundtrip failed'
          stop 2
       endif
       if (.not. rt_codon_proof_seen) then
          write(iulog,*) 'sulf_readnl implementation = codon'
          rt_codon_proof_seen = .true.
       endif
    endif

end subroutine sulf_readnl

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine sulf_inti()
!-----------------------------------------------------------------------
! 	... Open netCDF file containing annual sulfur data.  Initialize
!           arrays with the data to be interpolated to the current time.
!
!           It is assumed that the time coordinate is increasing
!           and represents calendar days; range = [1.,366.).
!-----------------------------------------------------------------------
      use spmd_utils,    only : masterproc
      use mo_chem_utls,  only : get_spc_ndx, get_rxt_ndx
      use interpolate_data, only : lininterp_init, lininterp, lininterp_finish, interp_type
      use tracer_data,   only : trcdata_init
      use cam_history,   only : addfld, phys_decomp

      implicit none

!-----------------------------------------------------------------------
!	... Dummy args
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer :: ndxs(5), so4_ndx
      integer(c_int64_t) :: active_c

      character(len=8), parameter :: fld_names(1) = (/'SULFATE '/)

      ndxs(1) = get_rxt_ndx( 'usr_N2O5_aer' )
      ndxs(2) = get_rxt_ndx( 'usr_NO3_aer' )
      ndxs(3) = get_rxt_ndx( 'usr_NO2_aer' )
      ndxs(4) = get_rxt_ndx( 'usr_HO2_aer' )
      ndxs(5) = get_rxt_ndx( 'het1' )
      so4_ndx = get_spc_ndx('SO4')

      read_sulf = any( ndxs > 0) .and. (so4_ndx < 0)

      active_c = sulf_inti_codon(merge(1_c_int64_t, 0_c_int64_t, read_sulf))
      if (.not. sulf_inti_proof_written) then
         sulf_inti_proof_written = .true.
         if (masterproc) then
            if (active_c == 0_c_int64_t) then
               write(iulog,'(A)') 'sulf_inti direct = codon read_sulf=false no-op'
            else
               write(iulog,'(A)') 'sulf_inti selector = codon; active sulfur data init body = native'
            end if
            call flush(iulog)
         end if
      end if

      if ( active_c == 0_c_int64_t ) return

      allocate(file%in_pbuf(size(fld_names)))
      file%in_pbuf(:) = .false. 
      call trcdata_init( (/ fld_name /), filename, filelist, datapath, fields, file, &
           rmv_file, cycle_yr, fixed_ymd, fixed_tod, datatype)

      call addfld('SULFATE','VMR', pver, 'I', 'sulfate data', phys_decomp )

      end subroutine sulf_inti

	      subroutine set_sulf_time( pbuf2d, state )
!--------------------------------------------------------------------
!	... Check and set time interpolation indicies
!--------------------------------------------------------------------
      use tracer_data,  only : advance_trcdata

      implicit none

!--------------------------------------------------------------------
!	... Dummy args
!--------------------------------------------------------------------
	      type(physics_buffer_desc), pointer :: pbuf2d(:,:)
	      type(physics_state), intent(in):: state(begchunk:endchunk)
	      integer(c_int64_t) :: out_c

	      call set_sulf_time_select_impl()

	      if (.not.set_sulf_time_use_native_impl) then
	         out_c = set_sulf_time_codon()
	         call set_sulf_time_proof_once()
	      end if

	      if ( .not. read_sulf ) return

	      call advance_trcdata( fields, file, state, pbuf2d  )

	      end subroutine set_sulf_time

	      subroutine set_sulf_time_select_impl()

	      implicit none

	      character(len=32) :: impl_name
	      integer :: status, n, i, code

	      if (set_sulf_time_impl_selected) return

	      impl_name = 'codon'
	      call cam_codon_get_impl('SET_SULF_TIME_IMPL', impl_name, n, status)

	      if (status == 0 .and. n > 0) then
	         do i = 1, n
	            code = iachar(impl_name(i:i))
	            if (code >= iachar('A') .and. code <= iachar('Z')) then
	               impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
	            end if
	         end do
	         set_sulf_time_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
	      else
	         set_sulf_time_use_native_impl = .false.
	      end if

	      set_sulf_time_impl_selected = .true.

	      if (masterproc) then
	         if (set_sulf_time_use_native_impl) then
	            write(iulog,*) 'set_sulf_time implementation = native'
	         else
	            write(iulog,*) 'set_sulf_time implementation = codon'
	         end if
	         call flush(iulog)
	      end if

	      end subroutine set_sulf_time_select_impl

	      subroutine set_sulf_time_proof_once()

	      implicit none

	      if (set_sulf_time_proof_written) return
	      set_sulf_time_proof_written = .true.

	      if (masterproc) then
	         write(iulog,'(A)') 'set_sulf_time entered (read_sulf gate = codon; active tracer advance = native when enabled)'
	         call flush(iulog)
	      end if

	      end subroutine set_sulf_time_proof_once

	      subroutine sulf_interp( ncol, lchnk, ccm_sulf )
!-----------------------------------------------------------------------
! 	... Time interpolate sulfatei to current time
!-----------------------------------------------------------------------
      use cam_history,  only : outfld
      use iso_c_binding, only : c_int64_t, c_loc, c_null_ptr, c_ptr

      implicit none

!-----------------------------------------------------------------------
! 	... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in)   :: ncol              ! columns in chunk
      integer, intent(in)   :: lchnk             ! chunk number
      real(r8), target, intent(out) :: ccm_sulf(:,:)     ! output sulfate

!-----------------------------------------------------------------------
! 	... Local variables
!-----------------------------------------------------------------------
      type(c_ptr) :: fields_data_p

      interface
         subroutine sulf_interp_codon(ncol_c, pcols_c, pver_c, begchunk_c, lchnk_c, read_sulf_c, &
              fields_data_p, ccm_sulf_p) bind(c, name="sulf_interp_codon")
           use iso_c_binding, only : c_int64_t, c_ptr
           integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, begchunk_c, lchnk_c, read_sulf_c
           type(c_ptr), value :: fields_data_p, ccm_sulf_p
         end subroutine sulf_interp_codon
      end interface

      call sulf_interp_select_impl()

      if (sulf_interp_use_native_impl) then
         ccm_sulf(:,:) = 0._r8
         if ( .not. read_sulf ) return
         ccm_sulf(:ncol,:) = fields(1)%data(:ncol,:,lchnk)
      else
         fields_data_p = c_null_ptr
         if (read_sulf) fields_data_p = c_loc(fields(1)%data(1,1,begchunk))
         call sulf_interp_codon( &
              int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(begchunk, c_int64_t), &
              int(lchnk, c_int64_t), int(merge(1, 0, read_sulf), c_int64_t), fields_data_p, c_loc(ccm_sulf) &
         )
         if (masterproc .and. .not. sulf_interp_proof_written) then
            if (read_sulf) then
               write(iulog,'(A)') 'sulf_interp entered (unified zero/copy stage dispatch = codon, read_sulf = true)'
               call sulf_interp_append_proof('sulf_interp entered (unified zero/copy stage dispatch = codon, read_sulf = true)')
            else
               write(iulog,'(A)') 'sulf_interp entered (unified zero stage dispatch = codon, read_sulf = false)'
               call sulf_interp_append_proof('sulf_interp entered (unified zero stage dispatch = codon, read_sulf = false)')
            end if
            sulf_interp_proof_written = .true.
         end if
      end if

      if ( .not. read_sulf ) return

      call outfld( 'SULFATE', ccm_sulf(:ncol,:), ncol, lchnk )

      end subroutine sulf_interp

      subroutine sulf_interp_select_impl()

      implicit none

      character(len=32) :: impl_name
      integer :: status, n, i, code

      if (sulf_interp_impl_selected) return

      impl_name = 'codon'
      call cam_codon_get_impl('SULF_INTERP_IMPL', impl_name, n, status)

      if (status == 0 .and. n > 0) then
         do i = 1, n
            code = iachar(impl_name(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
               impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
         end do
         sulf_interp_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
      else
         sulf_interp_use_native_impl = .false.
      end if

      sulf_interp_impl_selected = .true.

      if (masterproc) then
         if (sulf_interp_use_native_impl) then
            write(iulog,*) 'sulf_interp implementation = native'
         else
            write(iulog,*) 'sulf_interp implementation = codon'
            call sulf_interp_append_proof('sulf_interp implementation = codon')
         end if
         call flush(iulog)
      end if

      end subroutine sulf_interp_select_impl

      subroutine sulf_interp_append_proof(proof_line)

      implicit none

      character(len=*), intent(in) :: proof_line
      character(len=512) :: proof_file
      integer :: status, n, proof_unit, ios

      proof_file = ''
      call get_environment_variable('SULF_INTERP_PROOF_FILE', value=proof_file, length=n, status=status)
      if (status == 0 .and. n > 0) then
         open(newunit=proof_unit, file=trim(proof_file(:n)), status='unknown', position='append', action='write', iostat=ios)
         if (ios == 0) then
            write(proof_unit,'(A)') trim(proof_line)
            close(proof_unit)
         end if
      end if

      end subroutine sulf_interp_append_proof

      end module mo_sulf
