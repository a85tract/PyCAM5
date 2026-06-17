module cam3_ozone_data

!----------------------------------------------------------------------- 
! Purpose:
!
! Interpolates zonal ozone datasets used by CAM3 and puts the field 'O3' into
! the physics buffer.
! 
! Revision history:
! 2004-07-31  B. Eaton       Assemble module from comozp.F90, oznini.F90, oznint.F90, radozn.F90
! 2004-08-19  B. Eaton       Modify ozone_data_vert_interp to return mass mixing ratio.
! 2004-08-30  B. Eaton       Add ozone_data_get_cnst method.
! 2008 June   B. Eaton       Change name to cam3_ozone_data to support backwards compatibility
!                            for reading the CAM3 ozone data.  Add *_readnl method so module
!                            reads its own namelist.  Add cam3_ozone_data_on variable to
!                            turn the module on from the namelist.  By default it's off.
!-----------------------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8
use iso_c_binding,  only: c_int64_t
use spmd_utils,     only: masterproc
use ppgrid,         only: begchunk, endchunk, pcols, pver
use cam_abortutils, only: endrun
use cam_logfile,    only: iulog
use physics_types,  only: physics_state
use boundarydata,   only: boundarydata_type, boundarydata_init, boundarydata_update, &
                          boundarydata_vert_interp
use mpishorthand

implicit none
private
save

! Public methods
public ::&
   cam3_ozone_data_readnl,        &! get namelist input
   cam3_ozone_data_register,      &! register ozone with physics buffer
   cam3_ozone_data_init,          &! open dataset and spatially interpolate data bounding initial time
   cam3_ozone_data_timestep_init   ! interpolate to current time

! Namelist variables
logical, public    :: cam3_ozone_data_on = .false. ! switch to turn module on/off
logical            :: ozncyc = .true.  ! .true. => assume annual cycle ozone data
character(len=256) :: bndtvo = ' '     ! full pathname for time-variant ozone dataset

! Local
integer            :: oz_idx           ! index into phys_buffer for ozone

type(boundarydata_type) :: ozonedata
character(len=6), parameter, dimension(1) :: nc_name = (/'OZONE '/) ! constituent names

logical :: use_native_cam3_ozone_data_impl = .false.
logical :: cam3_ozone_data_impl_selected = .false.
logical :: cam3_ozone_data_proof_written = .false.
logical :: cam3_ozone_data_readnl_logged = .false.

interface
   function cam3_ozone_data_readnl_codon(flag_c, cyc_c) result(out_c) &
        bind(c, name="cam3_ozone_data_readnl_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: flag_c, cyc_c
      integer(c_int64_t) :: out_c
   end function cam3_ozone_data_readnl_codon
   function cam3_ozone_data_flag_codon(flag_c) result(out_c) bind(c, name="cam3_ozone_data_flag_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: flag_c
      integer(c_int64_t) :: out_c
   end function cam3_ozone_data_flag_codon
end interface

!================================================================================================
contains
!================================================================================================

subroutine cam3_ozone_data_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (cam3_ozone_data_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('CAM3_OZONE_DATA_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_cam3_ozone_data_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_cam3_ozone_data_impl = .false.
   end if

   cam3_ozone_data_impl_selected = .true.

   if (masterproc) then
      if (use_native_cam3_ozone_data_impl) then
         write(iulog,*) 'cam3_ozone_data implementation = native'
      else
         write(iulog,*) 'cam3_ozone_data implementation = codon'
      end if
   end if

end subroutine cam3_ozone_data_select_impl

!================================================================================================

subroutine cam3_ozone_data_proof_once()

   if (cam3_ozone_data_proof_written) return
   cam3_ozone_data_proof_written = .true.

   if (masterproc) then
      write(iulog,'(A)') 'cam3_ozone_data entered (runtime namelist flag helpers = codon)'
   end if

end subroutine cam3_ozone_data_proof_once

!================================================================================================

subroutine cam3_ozone_data_log_direct(logged, proof_line)

   logical, intent(inout) :: logged
   character(len=*), intent(in) :: proof_line

   if (logged) return
   logged = .true.

   if (masterproc) then
      write(iulog,'(A)') trim(proof_line)
   end if

end subroutine cam3_ozone_data_log_direct

!================================================================================================

logical function cam3_ozone_data_flag(flag)

   logical, intent(in) :: flag
   integer(c_int64_t) :: flag_c, out_c

   call cam3_ozone_data_select_impl()

   if (use_native_cam3_ozone_data_impl) then
      cam3_ozone_data_flag = flag
      return
   end if

   call cam3_ozone_data_proof_once()
   if (flag) then
      flag_c = 1_c_int64_t
   else
      flag_c = 0_c_int64_t
   end if
   out_c = cam3_ozone_data_flag_codon(flag_c)
   cam3_ozone_data_flag = out_c /= 0_c_int64_t

end function cam3_ozone_data_flag

!================================================================================================

subroutine cam3_ozone_data_finalize_flags()

   cam3_ozone_data_on = cam3_ozone_data_flag(cam3_ozone_data_on)
   ozncyc = cam3_ozone_data_flag(ozncyc)

end subroutine cam3_ozone_data_finalize_flags

!================================================================================================

subroutine cam3_ozone_data_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   integer(c_int64_t) :: out_c
   logical :: group_found
   character(len=*), parameter :: subname = 'cam3_ozone_data_readnl'

   namelist /cam3_ozone_data_nl/ cam3_ozone_data_on, bndtvo, ozncyc
   !-----------------------------------------------------------------------------

   call cam3_ozone_data_select_impl()
   if (.not. use_native_cam3_ozone_data_impl) then
      group_found = .false.
      if (masterproc) then
         unitn = getunit()
         open( unitn, file=trim(nlfile), status='old' )
         call find_group_name(unitn, 'cam3_ozone_data_nl', status=ierr)
         group_found = ierr == 0
         close(unitn)
         call freeunit(unitn)
      end if
#ifdef SPMD
      call mpibcast(group_found, 1, mpilog, 0, mpicom)
#endif
      if (.not. group_found) then
         call cam3_ozone_data_proof_once()
         out_c = cam3_ozone_data_readnl_codon( &
              merge(1_c_int64_t, 0_c_int64_t, cam3_ozone_data_on), &
              merge(1_c_int64_t, 0_c_int64_t, ozncyc))
         cam3_ozone_data_on = mod(out_c, 2_c_int64_t) /= 0_c_int64_t
         ozncyc = mod(out_c / 2_c_int64_t, 2_c_int64_t) /= 0_c_int64_t
         call cam3_ozone_data_log_direct(cam3_ozone_data_readnl_logged, &
              'cam3_ozone_data_readnl direct = codon')
         return
      end if
   end if

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'cam3_ozone_data_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, cam3_ozone_data_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(cam3_ozone_data_on, 1, mpilog, 0, mpicom)
   call mpibcast(bndtvo, len(bndtvo), mpichar, 0, mpicom)
   call mpibcast(ozncyc, 1, mpilog, 0, mpicom)
#endif

   call cam3_ozone_data_finalize_flags()

end subroutine cam3_ozone_data_readnl

!================================================================================================

subroutine cam3_ozone_data_register()
   use physics_buffer, only : pbuf_add_field, dtype_r8

   call pbuf_add_field('O3','physpkg',dtype_r8,(/pcols,pver/),oz_idx)

end subroutine cam3_ozone_data_register

!================================================================================================

subroutine cam3_ozone_data_init(phys_state)
!----------------------------------------------------------------------- 
! 
! Purpose: Do initial read of time-variant ozone boundary dataset, containing
!          ozone mixing ratios as a function of latitude and pressure.  Read two
!          consecutive months between which the current date lies.  Routine
!          RADOZ2 then evaluates the two path length integrals (with and without
!          pressure weighting) from zero to the interfaces between the input
!          levels.  It also stores the contribution to the integral from each
!          layer.
! 
! Method: Call appropriate netcdf wrapper routines and interpolate to model grid
! 
! Author: CCM Core Group
! Modified: P. Worley, August 2003, for chunking and performance optimization
!           J. Edwards, Dec 2005, functionality now performed by zonalbndrydata
!-----------------------------------------------------------------------

   use cam_history,      only: addfld, phys_decomp

   type(physics_state), intent(in) :: phys_state(begchunk:endchunk) 
   !-----------------------------------------------------------------------
    
   call addfld ('O3VMR', 'm3/m3', pver, 'A', 'Ozone volume mixing ratio', phys_decomp, sampling_seq='rad_lwsw')


   ! Initialize for one field (arg_4=1) and do not vertically interpolate (arg_6=3)
   call boundarydata_init(bndtvo, phys_state, nc_name, 1, ozonedata, 3)

   if (masterproc) then
      write(iulog,*)'cam3_ozone_data_init: Initializing CAM3 prescribed ozone'
      write(iulog,*)'Time-variant boundary dataset (ozone) is: ', trim(bndtvo)
      if (ozncyc) then
         write(iulog,*)'OZONE dataset will be reused for each model year'
      else
         write(iulog,*)'OZONE dataset will not be cycled'
      end if
   end if

end subroutine cam3_ozone_data_init

!================================================================================================

subroutine cam3_ozone_data_timestep_init(pbuf2d,  phys_state)
!----------------------------------------------------------------------- 
! 
! Purpose: Interpolate ozone mixing ratios to current time, reading in new monthly
!          data if necessary, and spatially interpolating it.
! 
! Method: Find next month of ozone data to interpolate.  Linearly interpolate 
!         vertically and horizontally
! 
!-----------------------------------------------------------------------

   
   use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_get_chunk

   
   type(physics_state), intent(in) :: phys_state(begchunk:endchunk) 
   type(physics_buffer_desc), pointer :: pbuf2d(:,:)
   real(r8),pointer :: tmpptr(:,:)

   integer lchnk
    
   call boundarydata_update(phys_state, ozonedata)

   do lchnk = begchunk, endchunk
      call pbuf_get_field(pbuf_get_chunk(pbuf2d, lchnk), oz_idx, tmpptr)
      call ozone_data_get_cnst(phys_state(lchnk), tmpptr)
   enddo

end subroutine cam3_ozone_data_timestep_init

!================================================================================================

subroutine ozone_data_get_cnst(state, q)

   use cam_history, only: outfld
   use physconst,   only: mwo3

   type(physics_state),  intent(in) :: state
   real(r8)                         :: q(:,:)     ! constituent mass mixing ratio

   ! local variables
   integer :: lchnk            ! chunk identifier
   integer :: i, k
   real(r8) :: ozmixin(pcols,ozonedata%levsiz)
   ! *** N.B. this hardwired mw of dry air needs to be changed to the share value
   real(r8), parameter :: mwdry = 28.9644_r8  ! Effective molecular weight of dry air (g/mol)
   real(r8), parameter :: mwr =  mwo3/mwdry   ! convert from the dataset values of vmr to mmr
   !-------------------------------------------------------------------------------

   lchnk = state%lchnk

   ozmixin=0._r8
   do k=1,ozonedata%levsiz
      do i=1,state%ncol
         ozmixin(i,k) = ozonedata%datainst(state%latmapback(i),k,lchnk,1)
      end do
   end do
   call boundarydata_vert_interp(lchnk, state%ncol, ozonedata%levsiz, &
                                 1, ozonedata%pin, state%pmid, ozmixin , q)

   call outfld('O3VMR', q, pcols, lchnk)

   do k=1,pver
      do i=1,state%ncol
         q(i,k) = mwr*q(i,k)
      end do
   end do
    
end subroutine ozone_data_get_cnst

!================================================================================================

end module cam3_ozone_data
