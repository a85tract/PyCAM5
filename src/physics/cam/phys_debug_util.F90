module phys_debug_util

!----------------------------------------------------------------------------------------

! Module to facilitate debugging of physics parameterizations.
!
! The user requests a location for debugging in lat/lon coordinates
! (degrees).  The initialization routine does a global search to find the
! column in the physics grid closest to the requested location.  The local
! indices of that column in the physics decomposition are stored as module
! data.  The user code then passes the local chunk index of the chunked
! data into the subroutine that will write diagnostic information for the
! column.  The function phys_debug_col returns the local column index if
! the column of interest is contained in the chunk, and zero otherwise.
! Printing is done only if a column index >0 is returned.
!
! Phil Rasch, B. Eaton, Feb 2008
!----------------------------------------------------------------------------------------

use shr_kind_mod,    only: r8 => shr_kind_r8
use iso_c_binding,   only: c_double, c_int64_t
use phys_grid,       only: phys_grid_find_col, get_rlat_p, get_rlon_p
use spmd_utils,      only: masterproc, iam
use cam_logfile,     only: iulog
use cam_abortutils,  only: endrun

implicit none
private
save

real(r8), parameter :: uninit_r8 = huge(1._r8)

! Public methods
public phys_debug_readnl  ! read namelist input
public phys_debug_init    ! initialize the method to a chunk and column
public phys_debug_col     ! return local column index in debug chunk

! Namelist variables
real(r8) :: phys_debug_lat = uninit_r8 ! latitude of requested debug column location in degrees
real(r8) :: phys_debug_lon = uninit_r8 ! longitude of requested debug column location in degrees


integer :: debchunk = -999            ! local index of the chuck we will debug
integer :: debcol   = -999            ! the column within the chunk we will debug

logical :: use_native_phys_debug_util_impl = .false.
logical :: phys_debug_util_impl_selected = .false.
logical :: phys_debug_util_proof_written = .false.

interface
   function phys_debug_value_codon(value_c) result(out_c) bind(c, name="phys_debug_value_codon")
      use iso_c_binding, only: c_double
      real(c_double), value :: value_c
      real(c_double) :: out_c
   end function phys_debug_value_codon

   function phys_debug_has_location_codon(lat_set_c, lon_set_c) result(out_c) &
        bind(c, name="phys_debug_has_location_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: lat_set_c, lon_set_c
      integer(c_int64_t) :: out_c
   end function phys_debug_has_location_codon

   function phys_debug_init_codon(lat_set_c, lon_set_c) result(out_c) bind(c, name="phys_debug_init_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: lat_set_c, lon_set_c
      integer(c_int64_t) :: out_c
   end function phys_debug_init_codon

   function phys_debug_col_codon(chunk_c, debchunk_c, debcol_c) result(out_c) &
        bind(c, name="phys_debug_col_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: chunk_c, debchunk_c, debcol_c
      integer(c_int64_t) :: out_c
   end function phys_debug_col_codon
end interface

!================================================================================
contains
!================================================================================

subroutine phys_debug_util_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (phys_debug_util_impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('PHYS_DEBUG_UTIL_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_phys_debug_util_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_phys_debug_util_impl = .false.
   end if

   phys_debug_util_impl_selected = .true.

   if (masterproc) then
      if (use_native_phys_debug_util_impl) then
         write(iulog,*) 'phys_debug_util implementation = native'
      else
         write(iulog,*) 'phys_debug_util implementation = codon'
      end if
   end if

end subroutine phys_debug_util_select_impl

!================================================================================

subroutine phys_debug_util_proof_once()

   if (phys_debug_util_proof_written) return
   phys_debug_util_proof_written = .true.

   if (masterproc) then
      write(iulog,'(A)') 'phys_debug_util entered (debug namelist/location helpers = codon)'
   end if

end subroutine phys_debug_util_proof_once

!================================================================================

real(r8) function phys_debug_value(value_in)

   real(r8), intent(in) :: value_in
   real(c_double) :: out_c

   call phys_debug_util_select_impl()

   if (use_native_phys_debug_util_impl) then
      phys_debug_value = value_in
      return
   end if

   call phys_debug_util_proof_once()
   out_c = phys_debug_value_codon(real(value_in, c_double))
   phys_debug_value = real(out_c, r8)

end function phys_debug_value

!================================================================================

logical function phys_debug_has_location(lat, lon)

   real(r8), intent(in) :: lat, lon
   integer(c_int64_t) :: lat_set_c, lon_set_c, out_c

   call phys_debug_util_select_impl()

   if (use_native_phys_debug_util_impl) then
      phys_debug_has_location = .not. (lat == uninit_r8 .or. lon == uninit_r8)
      return
   end if

   call phys_debug_util_proof_once()
   if (lat == uninit_r8) then
      lat_set_c = 0_c_int64_t
   else
      lat_set_c = 1_c_int64_t
   end if
   if (lon == uninit_r8) then
      lon_set_c = 0_c_int64_t
   else
      lon_set_c = 1_c_int64_t
   end if
   out_c = phys_debug_has_location_codon(lat_set_c, lon_set_c)
   phys_debug_has_location = out_c /= 0_c_int64_t

end function phys_debug_has_location

!================================================================================

subroutine phys_debug_finalize_values()

   phys_debug_lat = phys_debug_value(phys_debug_lat)
   phys_debug_lon = phys_debug_value(phys_debug_lon)

end subroutine phys_debug_finalize_values

!================================================================================

subroutine phys_debug_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'phys_debug_readnl'

   namelist /phys_debug_nl/ phys_debug_lat, phys_debug_lon
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'phys_debug_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, phys_debug_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(phys_debug_lat, 1, mpir8, 0, mpicom)
   call mpibcast(phys_debug_lon, 1, mpir8, 0, mpicom)
#endif

   call phys_debug_finalize_values()

end subroutine phys_debug_readnl

!================================================================================

subroutine phys_debug_init()

   integer  :: owner, lchunk, icol
   real(r8) :: deblat, deblon
   integer(c_int64_t) :: lat_set_c, lon_set_c, active_c
   !-----------------------------------------------------------------------------

   ! If no debug column specified then do nothing
   call phys_debug_util_select_impl()
   lat_set_c = merge(1_c_int64_t, 0_c_int64_t, phys_debug_lat /= uninit_r8)
   lon_set_c = merge(1_c_int64_t, 0_c_int64_t, phys_debug_lon /= uninit_r8)
   if (use_native_phys_debug_util_impl) then
      if (.not. (lat_set_c /= 0_c_int64_t .and. lon_set_c /= 0_c_int64_t)) return
   else
      call phys_debug_util_proof_once()
      active_c = phys_debug_init_codon(lat_set_c, lon_set_c)
      if (active_c == 0_c_int64_t) return
   end if

   ! User has specified a column location for debugging.  Find the closest
   ! column in the physics grid.
   call phys_grid_find_col(phys_debug_lat, phys_debug_lon, owner, lchunk, icol)

   ! If the column is owned by this process then save its local indices
   if (iam == owner) then
      debchunk         = lchunk
      debcol           = icol
      deblat           = get_rlat_p(lchunk, icol)*57.296_r8  ! approximate conversion for log output only
      deblon           = get_rlon_p(lchunk, icol)*57.296_r8
      write(iulog,*) 'phys_debug_init: debugging column at lat=', deblat, '  lon=', deblon
   end if

end subroutine phys_debug_init

!================================================================================

integer function phys_debug_col(chunk)

   integer,  intent(in) :: chunk
   integer(c_int64_t) :: out_c
   !-----------------------------------------------------------------------------

   call phys_debug_util_select_impl()
   if (use_native_phys_debug_util_impl) then
      if (chunk == debchunk) then
         phys_debug_col = debcol
      else
         phys_debug_col = 0
      endif
   else
      call phys_debug_util_proof_once()
      out_c = phys_debug_col_codon(int(chunk, c_int64_t), int(debchunk, c_int64_t), int(debcol, c_int64_t))
      phys_debug_col = int(out_c)
   end if

end function phys_debug_col

!================================================================================

end module phys_debug_util
