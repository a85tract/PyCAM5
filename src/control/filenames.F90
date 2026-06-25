module filenames

! Module and methods to handle filenames needed for the model. This 
! includes input filenames, and most output filenames that the model
! uses. All filenames that the model uses will use methods or data
! constructed by this module. In some cases (such as the cam_history module)
! other modules or routines will store the actual filenames used, but
! this module is used to determine the names.

use time_manager,     only: get_curr_date, get_prev_date
use shr_kind_mod,     only: shr_kind_cm, shr_kind_cl
use cam_abortutils,   only: endrun
use cam_logfile,      only: iulog
use spmd_utils,       only: masterproc
use iso_c_binding,    only: c_int64_t, c_loc

implicit none
private
save

public get_dir                                  ! Get the directory name from a full path
public interpret_filename_spec                  ! Interpret a filename specifier

character(shr_kind_cl), public :: ncdata = 'ncdata'       ! full pathname for initial dataset
character(shr_kind_cl), public :: bnd_topo = 'bnd_topo'   ! full pathname for topography dataset

character(shr_kind_cl), public :: absems_data = 'absems_data' ! full pathname for time-invariant absorption dataset

character(shr_kind_cm), public :: caseid = ' '  ! Case identifier
logical, public :: brnch_retain_casename = .false.

integer, parameter :: nlen = shr_kind_cl                ! String length

logical :: use_native_interpret_filename_spec_impl = .false.
logical :: interpret_filename_spec_impl_selected = .false.
logical :: interpret_filename_spec_logged = .false.

interface
   subroutine interpret_filename_spec_codon(spec_len_c, spec_ascii_p, has_number_c, number_c, &
        year_c, month_c, day_c, ncsec_c, case_len_c, case_ascii_p, out_len_c, out_ascii_p, &
        status_p) bind(c, name="interpret_filename_spec_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: spec_len_c, has_number_c, number_c
      integer(c_int64_t), value :: year_c, month_c, day_c, ncsec_c
      integer(c_int64_t), value :: case_len_c, out_len_c
      type(c_ptr), value :: spec_ascii_p, case_ascii_p, out_ascii_p, status_p
   end subroutine interpret_filename_spec_codon
end interface

!===============================================================================
CONTAINS
!===============================================================================

subroutine interpret_filename_spec_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (interpret_filename_spec_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('INTERPRET_FILENAME_SPEC_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_interpret_filename_spec_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_interpret_filename_spec_impl = .false.
   end if

   interpret_filename_spec_impl_selected = .true.

   if (masterproc) then
      if (use_native_interpret_filename_spec_impl) then
         write(iulog,*) 'interpret_filename_spec implementation = native'
      else
         write(iulog,*) 'interpret_filename_spec implementation = codon'
      end if
   end if

end subroutine interpret_filename_spec_select_impl

!===============================================================================

subroutine interpret_filename_spec_log_direct()

   if (interpret_filename_spec_logged) return
   interpret_filename_spec_logged = .true.

   if (masterproc) then
      write(iulog,'(A)') 'interpret_filename_spec direct = codon; time-manager/endrun native boundary'
      call flush(iulog)
   end if

end subroutine interpret_filename_spec_log_direct

!===============================================================================

subroutine filenames_pack_ascii(text, ascii, n)

   character(len=*), intent(in) :: text
   integer(c_int64_t), intent(out) :: ascii(:)
   integer, intent(in) :: n
   integer :: i

   do i = 1, n
      ascii(i) = int(iachar(text(i:i)), c_int64_t)
   end do

end subroutine filenames_pack_ascii

!===============================================================================

character(len=nlen) function get_dir( filepath )

! Return the directory from a filename with a full path

   ! arguments
   character(len=*), intent(in) :: filepath ! Full path for a filename

   ! local variables
   integer :: filenameposition   ! Character position for last character of directory
   !-----------------------------------------------------------------------------

   ! Get the directory name of the input dataset
   filenameposition = index( filepath, '/', back=.true. )
   if ( filenameposition == 0 )then
      get_dir  = './'
   else
      get_dir  = filepath(1:filenameposition)
   end if

end function get_dir

!===============================================================================

character(len=nlen) function interpret_filename_spec( filename_spec, number, prev, case, &
   yr_spec, mon_spec, day_spec, sec_spec )

! Create a filename from a filename specifier. The 
! filename specifyer includes codes for setting things such as the
! year, month, day, seconds in day, caseid, and tape number. This
! routine is private to filenames.F90
!
! Interpret filename specifyer string with: 
!
!      %c for case, 
!      %t for optional number argument sent into function
!      %y for year
!      %m for month
!      %d for day
!      %s for second
!      %% for the "%" character
!
! If the filename specifyer has spaces " ", they will be trimmed out
! of the resulting filename.

   ! arguments
   character(len=*), intent(in)           :: filename_spec   ! Filename specifier to use
   integer         , intent(in), optional :: number          ! Number to use for %t field
   logical         , intent(in), optional :: prev            ! If should label with previous time-step
   character(len=*), intent(in), optional :: case            ! Optional casename
   integer         , intent(in), optional :: yr_spec         ! Simulation year
   integer         , intent(in), optional :: mon_spec        ! Simulation month
   integer         , intent(in), optional :: day_spec        ! Simulation day
   integer         , intent(in), optional :: sec_spec        ! Seconds into current simulation day

   ! Local variables
   integer :: year  ! Simulation year
   integer :: month ! Simulation month
   integer :: day   ! Simulation day
   integer :: ncsec ! Seconds into current simulation day
   character(len=nlen) :: string    ! Temporary character string 
   character(len=nlen) :: format    ! Format character string 
   integer :: i, n  ! Loop variables
   logical :: previous              ! If should label with previous time-step
   logical :: done
   integer :: spec_n, case_n
   integer(c_int64_t), target :: spec_codes(len(filename_spec)), case_codes(nlen)
   integer(c_int64_t), target :: out_codes(nlen), status_c
   integer(c_int64_t) :: has_number_c, number_c
   character(len=nlen) :: case_work
   !-----------------------------------------------------------------------------

   call interpret_filename_spec_select_impl()

   if ( len_trim(filename_spec) == 0 )then
      call endrun ('INTERPRET_FILENAME_SPEC: filename specifier is empty')
   end if
   if ( index(trim(filename_spec)," ") /= 0 )then
      call endrun ('INTERPRET_FILENAME_SPEC: filename specifier can not contain a space:'//trim(filename_spec))
   end if
   !
   ! Determine year, month, day and sec to put in filename
   !
   if (present(yr_spec) .and. present(mon_spec) .and. present(day_spec) .and. present(sec_spec)) then
      year  = yr_spec
      month = mon_spec
      day   = day_spec
      ncsec = sec_spec
   else
      if ( .not. present(prev) ) then
         previous = .false.
      else
         previous = prev
      end if
      if ( previous ) then
         call get_prev_date(year, month, day, ncsec)
      else
         call get_curr_date(year, month, day, ncsec)
      end if
   end if

   if (.not. use_native_interpret_filename_spec_impl) then
      spec_n = len_trim(filename_spec)
      case_work = ''
      if (present(case)) then
         case_work = case
      else
         case_work = caseid
      end if
      case_n = len_trim(case_work)
      call filenames_pack_ascii(filename_spec, spec_codes, spec_n)
      if (case_n > 0) then
         call filenames_pack_ascii(case_work, case_codes, case_n)
      else
         case_codes(1) = int(iachar(' '), c_int64_t)
      end if
      if (present(number)) then
         has_number_c = 1_c_int64_t
         number_c = int(number, c_int64_t)
      else
         has_number_c = 0_c_int64_t
         number_c = 0_c_int64_t
      end if
      status_c = -1_c_int64_t
      call interpret_filename_spec_log_direct()
      call interpret_filename_spec_codon(int(spec_n, c_int64_t), c_loc(spec_codes(1)), &
         has_number_c, number_c, int(year, c_int64_t), int(month, c_int64_t), &
         int(day, c_int64_t), int(ncsec, c_int64_t), int(case_n, c_int64_t), &
         c_loc(case_codes(1)), int(nlen, c_int64_t), c_loc(out_codes(1)), c_loc(status_c))
      select case (int(status_c))
      case (0)
         interpret_filename_spec = ''
         do i = 1, nlen
            interpret_filename_spec(i:i) = achar(int(out_codes(i)))
         end do
         return
      case (1)
         call endrun ('INTERPRET_FILENAME_SPEC: number needed in filename_spec')
      case (2)
         call endrun ('INTERPRET_FILENAME_SPEC: number is too large')
      case (3)
         call endrun ('INTERPRET_FILENAME_SPEC: Invalid expansion character')
      case (4)
         call endrun ('INTERPRET_FILENAME_SPEC: Resultant filename too long')
      case (5)
         call endrun ('INTERPRET_FILENAME_SPEC: Resulting filename is empty')
      case default
         call endrun ('INTERPRET_FILENAME_SPEC: Codon expansion failed')
      end select
   end if
   !
   ! Go through each character in the filename specifyer and interpret if special string
   !
   i = 1
   interpret_filename_spec = ''
   do while ( i <= len_trim(filename_spec) )
      !
      ! If following is an expansion string
      !
      if ( filename_spec(i:i) == "%" )then
         i = i + 1
         select case( filename_spec(i:i) )
         case( 'c' )   ! caseid
            if ( present(case) )then
               string = trim(case)
            else
               string = trim(caseid)
            end if
         case( 't' )   ! number
            if ( .not. present(number) )then
               write(iulog,*) 'INTERPRET_FILENAME_SPEC: number needed in filename_spec' &
                  , ', but not provided to subroutine'
               write(iulog,*) 'filename_spec = ', filename_spec
               call endrun
            end if
            if (      number > 999 ) then
               format = '(i4.4)'
               if ( number > 9999 ) then
                  write(iulog,*) 'INTERPRET_FILENAME_SPEC: number is too large: ', number
                  call endrun
               end if
            else if ( number > 99  ) then
               format = '(i3.3)'
            else if ( number > 9   ) then
               format = '(i2.2)'
            else
               format = '(i1.1)'
            end if
            write(string,format) number
         case( 'y' )   ! year
            if ( year > 99999   ) then
               format = '(i6.6)'
            else if ( year > 9999    ) then
               format = '(i5.5)'
            else
               format = '(i4.4)'
            end if
            write(string,format) year
         case( 'm' )   ! month
            write(string,'(i2.2)') month
         case( 'd' )   ! day
            write(string,'(i2.2)') day
         case( 's' )   ! second
            write(string,'(i5.5)') ncsec
         case( '%' )   ! percent character
            string = "%"
         case default
            call endrun ('INTERPRET_FILENAME_SPEC: Invalid expansion character: '//filename_spec(i:i))
         end select
         !
         ! Otherwise take normal text up to the next "%" character
         !
      else
         n = index( filename_spec(i:), "%" )
         if ( n == 0 ) n = len_trim( filename_spec(i:) ) + 1
         if ( n == 0 ) exit 
         string = filename_spec(i:n+i-2)
         i = n + i - 2
      end if
      if ( len_trim(interpret_filename_spec) == 0 )then
         interpret_filename_spec = trim(string)
      else
         if ( (len_trim(interpret_filename_spec)+len_trim(string)) >= nlen )then
            call endrun ('INTERPRET_FILENAME_SPEC: Resultant filename too long')
         end if
         interpret_filename_spec = trim(interpret_filename_spec) // trim(string)
      end if
      i = i + 1

   end do
   if ( len_trim(interpret_filename_spec) == 0 )then
      call endrun ('INTERPRET_FILENAME_SPEC: Resulting filename is empty')
   end if

end function interpret_filename_spec

end module filenames
