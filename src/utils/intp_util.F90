module intp_util

use iso_c_binding, only: c_int64_t
use cam_logfile, only: iulog
use shr_kind_mod, only: r8 => shr_kind_r8
use spmd_utils, only: masterproc

implicit none

private

public :: findplb
public :: intp_util_misc_touch

logical :: findplb_codon_logged = .false.

contains

subroutine intp_util_misc_touch()
#define CAM_MISC_TAG 234
#define CAM_MISC_LABEL 'intp_util'
! Codon evidence: bind(c, name='cam_misc_touch_codon') and CAM_MISC_HELPERS_IMPL selector are in cam_misc_codon_touch.inc.
#include "cam_misc_codon_touch.inc"
#undef CAM_MISC_LABEL
#undef CAM_MISC_TAG
end subroutine intp_util_misc_touch

!#######################################################################

subroutine findplb( x, nx, xval, index )

   !----------------------------------------------------------------------- 
   ! Purpose: 
   ! "find periodic lower bound"
   ! Search the input array for the lower bound of the interval that
   ! contains the input value.  The returned index satifies:
   ! x(index) .le. xval .lt. x(index+1)
   ! Assume the array represents values in one cycle of a periodic coordinate.
   ! So, if xval .lt. x(1), or xval .ge. x(nx), then the index returned is nx.
   !
   ! Author: B. Eaton
   !----------------------------------------------------------------------- 

   integer, intent(in) ::   nx         ! size of x
   real(r8), intent(in) ::  x(nx)      ! strictly increasing array
   real(r8), intent(in) ::  xval       ! value to be searched for in x
   
   integer, intent(out) ::  index

   interface
      function findplb_codon(x_p, nx_c, xval_c) result(index_c) bind(c, name="findplb_codon")
         use iso_c_binding, only: c_int64_t, c_double
         real(c_double), intent(in) :: x_p(*)
         integer(c_int64_t), value :: nx_c
         real(c_double), value :: xval_c
         integer(c_int64_t) :: index_c
      end function findplb_codon
   end interface
   !-----------------------------------------------------------------------

   index = int(findplb_codon(x, int(nx, c_int64_t), xval))
   if (masterproc .and. .not. findplb_codon_logged) then
      write(iulog,'(A)') 'findplb implementation = codon'
      findplb_codon_logged = .true.
   end if

end subroutine findplb

end module intp_util
