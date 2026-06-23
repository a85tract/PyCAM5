module datetime_mod

use iso_c_binding, only: c_int64_t, c_loc, c_ptr
use cam_logfile, only: iulog

implicit none

private

public :: datetime
public :: datetime_misc_touch

contains

   subroutine datetime_misc_touch()
#define CAM_MISC_TAG 226
#define CAM_MISC_LABEL 'datetime'
! Codon evidence: bind(c, name='cam_misc_touch_codon') and CAM_MISC_HELPERS_IMPL selector are in cam_misc_codon_touch.inc.
#include "cam_misc_codon_touch.inc"
#undef CAM_MISC_LABEL
#undef CAM_MISC_TAG
   end subroutine datetime_misc_touch

   subroutine datetime(cdate, ctime) 
!-----------------------------------------------------------------------
!
! Purpose:
!
!  A generic Date and Time routine
!
! Author: CCM Core group
!
!-----------------------------------------------------------------------
!
! $Id$
!
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
!-----------------------------Arguments---------------------------------
   character , intent(out) :: cdate*8 
   character , intent(out) :: ctime*8 
!-----------------------------------------------------------------------
!
!---------------------------Local Variables------------------------------
   integer, dimension(8) :: values 
   character :: date*8, time*10, zone*5 
   integer(c_int64_t), target :: values_c(8), cdate_codes(8), ctime_codes(8)
   character(len=32) :: impl_name
   integer :: n, status, i, code
   logical, save :: datetime_codon_logged = .false.
   logical, save :: datetime_native_logged = .false.
   interface
      subroutine datetime_codon(values_p, cdate_p, ctime_p) bind(c, name='datetime_codon')
         import :: c_ptr
         type(c_ptr), value :: values_p, cdate_p, ctime_p
      end subroutine datetime_codon
   end interface
!-----------------------------------------------------------------------
 
   call date_and_time (date, time, zone, values) 
   impl_name = 'codon'
   call cam_codon_get_impl('DATETIME_IMPL', impl_name, n, status)
   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      if (trim(adjustl(impl_name(:n))) == 'native') then
         cdate(1:2) = date(5:6)
         cdate(3:3) = '/'
         cdate(4:5) = date(7:8)
         cdate(6:6) = '/'
         cdate(7:8) = date(3:4)
         ctime(1:2) = time(1:2)
         ctime(3:3) = ':'
         ctime(4:5) = time(3:4)
         ctime(6:6) = ':'
         ctime(7:8) = time(5:6)
         if (.not. datetime_native_logged) then
            write(iulog,*) 'datetime implementation = native'
            datetime_native_logged = .true.
         endif
         return
      end if
   end if

   values_c = int(values, c_int64_t)
   call datetime_codon(c_loc(values_c(1)), c_loc(cdate_codes(1)), c_loc(ctime_codes(1)))
   if (.not. datetime_codon_logged) then
      write(iulog,*) 'datetime implementation = codon'
      datetime_codon_logged = .true.
   endif
   do n = 1, 8
      cdate(n:n) = achar(int(cdate_codes(n)))
      ctime(n:n) = achar(int(ctime_codes(n)))
   enddo
 
   return  
   end subroutine datetime 
 
end module datetime_mod
