module string_utils

   use iso_c_binding, only: c_int64_t, c_loc, c_ptr
   use cam_logfile,   only: iulog

   implicit none
   private

! Public interface methods

   public ::&
      to_upper, &   ! Convert character string to upper case
      to_lower, &   ! Convert character string to lower case
      INCSTR, &     ! increments a string
      GLC, &        ! Position of last significant character in string
      string_utils_misc_touch

contains

subroutine string_utils_misc_touch()
#define CAM_MISC_TAG 223
#define CAM_MISC_LABEL 'string_utils'
! Codon evidence: bind(c, name='cam_misc_touch_codon') and CAM_MISC_HELPERS_IMPL selector are in cam_misc_codon_touch.inc.
#include "cam_misc_codon_touch.inc"
#undef CAM_MISC_LABEL
#undef CAM_MISC_TAG
end subroutine string_utils_misc_touch

logical function string_utils_use_native(selector)
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
      string_utils_use_native = trim(adjustl(impl_name(:n))) == 'native'
   else
      string_utils_use_native = .false.
   end if
end function string_utils_use_native

function to_upper(str)

!----------------------------------------------------------------------- 
! Purpose: 
! Convert character string to upper case.
! 
! Method: 
! Use achar and iachar intrinsics to ensure use of ascii collating sequence.
!
! Author:  B. Eaton, July 2001
!     
! $Id$
!----------------------------------------------------------------------- 
   implicit none

   character(len=*), intent(in) :: str      ! String to convert to upper case
   character(len=len(str))      :: to_upper

! Local variables

   integer :: i                ! Index
   integer :: aseq             ! ascii collating sequence
   integer :: lower_to_upper   ! integer to convert case
   character(len=1) :: ctmp    ! Character temporary
   integer(c_int64_t), target :: in_codes(len(str)), out_codes(len(str))
   logical, save :: to_upper_codon_logged = .false.
   logical, save :: to_upper_native_logged = .false.
   interface
      subroutine to_upper_codon(in_p, out_p, n) bind(c, name='to_upper_codon')
         import :: c_int64_t, c_ptr
         type(c_ptr), value :: in_p, out_p
         integer(c_int64_t), value :: n
      end subroutine to_upper_codon
   end interface
!-----------------------------------------------------------------------
   if (.not. string_utils_use_native('TO_UPPER_IMPL')) then
      do i = 1, len(str)
         in_codes(i) = int(iachar(str(i:i)), c_int64_t)
      end do
      call to_upper_codon(c_loc(in_codes(1)), c_loc(out_codes(1)), int(len(str), c_int64_t))
      do i = 1, len(str)
         to_upper(i:i) = achar(int(out_codes(i)))
      end do
      if (.not. to_upper_codon_logged) then
         write(iulog,*) 'to_upper implementation = codon'
         to_upper_codon_logged = .true.
      endif
      return
   endif

   if (.not. to_upper_native_logged) then
      write(iulog,*) 'to_upper implementation = native'
      to_upper_native_logged = .true.
   endif

   lower_to_upper = iachar("A") - iachar("a")

   do i = 1, len(str)
      ctmp = str(i:i)
      aseq = iachar(ctmp)
      if ( aseq >= iachar("a") .and. aseq <= iachar("z") ) &
           ctmp = achar(aseq + lower_to_upper)
      to_upper(i:i) = ctmp
   end do

end function to_upper

function to_lower(str)

!----------------------------------------------------------------------- 
! Purpose: 
! Convert character string to lower case.
! 
! Method: 
! Use achar and iachar intrinsics to ensure use of ascii collating sequence.
!
! Author:  B. Eaton, July 2001
!     
! $Id$
!----------------------------------------------------------------------- 
   implicit none

   character(len=*), intent(in) :: str      ! String to convert to lower case
   character(len=len(str))      :: to_lower

! Local variables

   integer :: i                ! Index
   integer :: aseq             ! ascii collating sequence
   integer :: upper_to_lower   ! integer to convert case
   character(len=1) :: ctmp    ! Character temporary
   integer(c_int64_t), target :: in_codes(len(str)), out_codes(len(str))
   logical, save :: to_lower_codon_logged = .false.
   logical, save :: to_lower_native_logged = .false.
   interface
      subroutine to_lower_codon(in_p, out_p, n) bind(c, name='to_lower_codon')
         import :: c_int64_t, c_ptr
         type(c_ptr), value :: in_p, out_p
         integer(c_int64_t), value :: n
      end subroutine to_lower_codon
   end interface
!-----------------------------------------------------------------------
   if (.not. string_utils_use_native('TO_LOWER_IMPL')) then
      do i = 1, len(str)
         in_codes(i) = int(iachar(str(i:i)), c_int64_t)
      end do
      call to_lower_codon(c_loc(in_codes(1)), c_loc(out_codes(1)), int(len(str), c_int64_t))
      do i = 1, len(str)
         to_lower(i:i) = achar(int(out_codes(i)))
      end do
      if (.not. to_lower_codon_logged) then
         write(iulog,*) 'to_lower implementation = codon'
         to_lower_codon_logged = .true.
      endif
      return
   endif

   if (.not. to_lower_native_logged) then
      write(iulog,*) 'to_lower implementation = native'
      to_lower_native_logged = .true.
   endif

   upper_to_lower = iachar("a") - iachar("A")

   do i = 1, len(str)
      ctmp = str(i:i)
      aseq = iachar(ctmp)
      if ( aseq >= iachar("A") .and. aseq <= iachar("Z") ) &
           ctmp = achar(aseq + upper_to_lower)
      to_lower(i:i) = ctmp
   end do

end function to_lower

integer function INCSTR( s, inc )
  !-----------------------------------------------------------------------
  ! 	... Increment a string whose ending characters are digits.
  !           The incremented integer must be in the range [0 - (10**n)-1]
  !           where n is the number of trailing digits.
  !           Return values:
  !
  !            0 success
  !           -1 error: no trailing digits in string
  !           -2 error: incremented integer is out of range
  !-----------------------------------------------------------------------

  implicit none

  !-----------------------------------------------------------------------
  ! 	... Dummy variables
  !-----------------------------------------------------------------------
  integer, intent(in) :: &
       inc                                       ! value to increment string (may be negative)
  character(len=*), intent(inout) :: &
       s                                         ! string with trailing digits


  !-----------------------------------------------------------------------
  ! 	... Local variables
  !-----------------------------------------------------------------------
  integer :: &
       i, &                          ! index
       lstr, &                       ! number of significant characters in string
       lnd, &                        ! position of last non-digit
       ndigit, &                     ! number of trailing digits
       ival, &                       ! integer value of trailing digits
       pow, &                        ! power of ten
       digit                         ! integer value of a single digit

  lstr   = GLC( s )
  lnd    = LASTND( s )
  ndigit = lstr - lnd

  if( ndigit == 0 ) then
     INCSTR = -1
     return
  end if

  !-----------------------------------------------------------------------
  !     	... Calculate integer corresponding to trailing digits.
  !-----------------------------------------------------------------------
  ival = 0
  pow  = 0
  do i = lstr,lnd+1,-1
     digit = ICHAR(s(i:i)) - ICHAR('0')
     ival  = ival + digit * 10**pow
     pow   = pow + 1
  end do

  !-----------------------------------------------------------------------
  !     	... Increment the integer
  !-----------------------------------------------------------------------
  ival = ival + inc
  if( ival < 0 .or. ival > 10**ndigit-1 ) then
     INCSTR = -2
     return
  end if

  !-----------------------------------------------------------------------
  !     	... Overwrite trailing digits
  !-----------------------------------------------------------------------
  pow = ndigit
  do i = lnd+1,lstr
     digit  = MOD( ival,10**pow ) / 10**(pow-1)
     s(i:i) = CHAR( ICHAR('0') + digit )
     pow    = pow - 1
  end do

  INCSTR = 0

end function INCSTR

integer function LASTND( cs )
  !-----------------------------------------------------------------------
  ! 	... Position of last non-digit in the first input token.
  ! 	    Return values:
  !     	    > 0  => position of last non-digit
  !     	    = 0  => token is all digits (or empty)
  !-----------------------------------------------------------------------

  implicit none

  !-----------------------------------------------------------------------
  ! 	... Dummy arguments
  !-----------------------------------------------------------------------
  character(len=*), intent(in) :: cs       !  Input character string

  !-----------------------------------------------------------------------
  ! 	... Local variables
  !-----------------------------------------------------------------------
  integer :: n, nn, digit

  n = GLC( cs )
  if( n == 0 ) then     ! empty string
     LASTND = 0
     return
  end if

  do nn = n,1,-1
     digit = ICHAR( cs(nn:nn) ) - ICHAR('0')
     if( digit < 0 .or. digit > 9 ) then
        LASTND = nn
        return
     end if
  end do

  LASTND = 0    ! all characters are digits

end function LASTND

integer function GLC( cs )
  !-----------------------------------------------------------------------
  ! 	... Position of last significant character in string. 
  !           Here significant means non-blank or non-null.
  !           Return values:
  !               > 0  => position of last significant character
  !               = 0  => no significant characters in string
  !-----------------------------------------------------------------------

  implicit none

  !-----------------------------------------------------------------------
  ! 	... Dummy arguments
  !-----------------------------------------------------------------------
  character(len=*), intent(in) :: cs       !  Input character string

  !-----------------------------------------------------------------------
  ! 	... Local variables
  !-----------------------------------------------------------------------
  integer :: l, n
  integer(c_int64_t), target :: cs_codes(len(cs))
  logical, save :: glc_codon_logged = .false.
  logical, save :: glc_native_logged = .false.
  interface
     function glc_codon(chars_p, n) result(last_pos) bind(c, name='glc_codon')
        import :: c_int64_t, c_ptr
        type(c_ptr), value :: chars_p
        integer(c_int64_t), value :: n
        integer(c_int64_t) :: last_pos
     end function glc_codon
  end interface

  if (.not. string_utils_use_native('GLC_IMPL')) then
     do n = 1, len(cs)
        cs_codes(n) = int(iachar(cs(n:n)), c_int64_t)
     end do
     GLC = int(glc_codon(c_loc(cs_codes(1)), int(len(cs), c_int64_t)))
     if (.not. glc_codon_logged) then
        write(iulog,*) 'glc implementation = codon'
        glc_codon_logged = .true.
     endif
     return
  endif

  if (.not. glc_native_logged) then
     write(iulog,*) 'glc implementation = native'
     glc_native_logged = .true.
  endif

  l = LEN( cs )
  if( l == 0 ) then
     GLC = 0
     return
  end if

  do n = l,1,-1
     if( cs(n:n) /= ' ' .and. cs(n:n) /= CHAR(0) ) then
        exit
     end if
  end do
  GLC = n

end function GLC

end module string_utils
