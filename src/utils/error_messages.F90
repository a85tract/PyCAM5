module error_messages

   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: 
   ! General purpose routines for issuing error messages.
   ! 
   ! Author: B. Eaton
   ! 
   !----------------------------------------------------------------------- 
   use cam_abortutils, only: endrun
   use cam_logfile,    only: iulog
   use iso_c_binding,  only: c_int64_t

   implicit none
   save
   private
   public :: &
      alloc_err,      &! Issue error message after non-zero return from an allocate statement.
      handle_err,     &! Issue error message after non-zero return from anything
      handle_ncerr     ! Handle error returns from netCDF library procedures.

   ! If an error message string is not empty, abort with that string as the
   ! error message.
   public :: handle_errmsg
   public :: error_messages_misc_touch

!##############################################################################
contains
!##############################################################################

   logical function error_messages_use_native(selector)
      character(len=*), intent(in) :: selector
      character(len=32) :: impl_name
      integer :: n, status, i, code

      error_messages_use_native = .false.
      impl_name = 'codon'
      call cam_codon_get_impl(selector, impl_name, n, status)
      if (status == 0 .and. n > 0) then
         do i = 1, n
            code = iachar(impl_name(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
               impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
         end do
         error_messages_use_native = trim(adjustl(impl_name(:n))) == 'native'
      end if
   end function error_messages_use_native

!##############################################################################

   subroutine error_messages_misc_touch()
#define CAM_MISC_TAG 225
#define CAM_MISC_LABEL 'error_messages'
! Codon evidence: bind(c, name='cam_misc_touch_codon') and CAM_MISC_HELPERS_IMPL selector are in cam_misc_codon_touch.inc.
#include "cam_misc_codon_touch.inc"
#undef CAM_MISC_LABEL
#undef CAM_MISC_TAG
   end subroutine error_messages_misc_touch

!##############################################################################

   subroutine alloc_err( istat, routine, name, nelem )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Issue error message after non-zero return from an allocate statement.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      integer, intent(in) ::&
         istat           ! status from allocate statement
      character(len=*), intent(in) ::&
         routine,       &! routine that called allocate
         name            ! name of array
      integer, intent(in) ::&
         nelem           ! number of elements attempted to allocate
      interface
         function alloc_err_codon(istat) result(is_error) bind(c, name='alloc_err_codon')
           import :: c_int64_t
           integer(c_int64_t), value :: istat
           integer(c_int64_t) :: is_error
         end function alloc_err_codon
      end interface
      logical, save :: alloc_err_codon_logged = .false.
      logical, save :: alloc_err_native_logged = .false.
      !-----------------------------------------------------------------------

      if (error_messages_use_native('ALLOC_ERR_IMPL')) then
         if (.not. alloc_err_native_logged) then
            write(iulog,*) 'alloc_err implementation = native'
            alloc_err_native_logged = .true.
         end if
         if ( istat .ne. 0 ) then
            write(iulog,*)'ERROR trying to allocate memory in routine: ' &
                      //trim(routine)
            write(iulog,*)'  Variable name: '//trim(name)
            write(iulog,*)'  Number of elements: ',nelem
            call endrun ('ALLOC_ERR')
         end if
         return
      end if

      if (.not. alloc_err_codon_logged) then
         write(iulog,*) 'alloc_err implementation = codon'
         alloc_err_codon_logged = .true.
      end if
      if (alloc_err_codon(int(istat, c_int64_t)) /= 0_c_int64_t) then
         write(iulog,*)'ERROR trying to allocate memory in routine: ' &
                   //trim(routine)
         write(iulog,*)'  Variable name: '//trim(name)
         write(iulog,*)'  Number of elements: ',nelem
         call endrun ('ALLOC_ERR')
      end if

      return

   end subroutine alloc_err

!##############################################################################

   subroutine handle_err( istat, msg )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Issue error message after non-zero return from anything.
      !
      ! Author: T. Henderson
      !----------------------------------------------------------------------- 

      integer,          intent(in) :: istat  ! status, zero = "no error"
      character(len=*), intent(in) :: msg    ! error message to print
      interface
         function handle_err_codon(istat) result(is_error) bind(c, name='handle_err_codon')
           import :: c_int64_t
           integer(c_int64_t), value :: istat
           integer(c_int64_t) :: is_error
         end function handle_err_codon
      end interface
      logical, save :: handle_err_codon_logged = .false.
      logical, save :: handle_err_native_logged = .false.
      !-----------------------------------------------------------------------

      if (error_messages_use_native('HANDLE_ERR_IMPL')) then
         if (.not. handle_err_native_logged) then
            write(iulog,*) 'handle_err implementation = native'
            handle_err_native_logged = .true.
         end if
         if ( istat .ne. 0 ) then
            call endrun (trim(msg))
         end if
         return
      end if

      if (.not. handle_err_codon_logged) then
         write(iulog,*) 'handle_err implementation = codon'
         handle_err_codon_logged = .true.
      end if
      if (handle_err_codon(int(istat, c_int64_t)) /= 0_c_int64_t) then
         call endrun (trim(msg))
      end if

      return

   end subroutine handle_err

!##############################################################################

   subroutine handle_ncerr( ret, mes, line )
      
      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Check netCDF library function return code.  If error detected 
      ! issue error message then abort.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

!-----------------------------------------------------------------------
     use netcdf
!-----------------------------------------------------------------------

      integer, intent(in) ::&
         ret                 ! return code from netCDF library routine
      character(len=*), intent(in) ::&
         mes                 ! message to be printed if error detected
      integer, intent(in), optional :: line
      interface
         function handle_ncerr_codon(ret, noerr) result(is_error) &
              bind(c, name='handle_ncerr_codon')
           import :: c_int64_t
           integer(c_int64_t), value :: ret
           integer(c_int64_t), value :: noerr
           integer(c_int64_t) :: is_error
         end function handle_ncerr_codon
      end interface
      logical, save :: handle_ncerr_codon_logged = .false.
      logical, save :: handle_ncerr_native_logged = .false.
      !-----------------------------------------------------------------------

      if (error_messages_use_native('HANDLE_NCERR_IMPL')) then
         if (.not. handle_ncerr_native_logged) then
            write(iulog,*) 'handle_ncerr implementation = native'
            handle_ncerr_native_logged = .true.
         end if
         if ( ret .ne. NF90_NOERR ) then
            if(present(line)) then
               write(iulog,*) mes, line
            else
               write(iulog,*) mes
            end if
            write(iulog,*) nf90_strerror( ret )
            call endrun ('HANDLE_NCERR')
         endif
         return
      end if

      if (.not. handle_ncerr_codon_logged) then
         write(iulog,*) 'handle_ncerr implementation = codon'
         handle_ncerr_codon_logged = .true.
      end if
      if (handle_ncerr_codon(int(ret, c_int64_t), int(NF90_NOERR, c_int64_t)) /= 0_c_int64_t) then
         if(present(line)) then
            write(iulog,*) mes, line
         else	
            write(iulog,*) mes
         end if
         write(iulog,*) nf90_strerror( ret )
         call endrun ('HANDLE_NCERR')
      endif

      return

   end subroutine handle_ncerr

!##############################################################################

   subroutine handle_errmsg(errmsg, subname, extra_msg)

     ! String that is asserted to be null.
     character(len=*), intent(in)           :: errmsg
     ! Name of procedure generating the message.
     character(len=*), intent(in), optional :: subname
     ! Additional message from the procedure calling this one.
     character(len=*), intent(in), optional :: extra_msg
     integer(c_int64_t) :: active_c
     logical, save :: proof_logged = .false.

     interface
        function handle_errmsg_codon(active) result(out_c) bind(c, name='handle_errmsg_codon')
          import :: c_int64_t
          integer(c_int64_t), value :: active
          integer(c_int64_t) :: out_c
        end function handle_errmsg_codon
     end interface

     if (.not. error_messages_use_native('HANDLE_ERRMSG_IMPL')) then
        active_c = handle_errmsg_codon(merge(1_c_int64_t, 0_c_int64_t, len_trim(errmsg) /= 0))
        if (active_c == 0_c_int64_t) then
           if (.not. proof_logged) then
              write(iulog,'(A)') 'handle_errmsg direct = codon empty-message no-op'
              call flush(iulog)
              proof_logged = .true.
           end if
           return
        end if
     end if

     if (trim(errmsg) /= "") then

        if (present(extra_msg)) &
             write(iulog,*) "handle_errmsg: &
             &Message from caller: ",trim(extra_msg)

        if (present(subname)) then
           call endrun("ERROR: handle_errmsg: "// &
                trim(subname)//": "//trim(errmsg))
        else
           call endrun("ERROR: handle_errmsg: "// &
                "Error message received from routine: "//trim(errmsg))
        end if

     end if

   end subroutine handle_errmsg

!##############################################################################

end module error_messages
