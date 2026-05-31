module dycore
!
! Data and utility routines related to the dycore
!
   implicit none

PRIVATE

   public :: dycore_is, get_resolution

CONTAINS

   logical function dycore_is (name)
      use iso_c_binding, only : c_int64_t
      use cam_logfile, only : iulog
!
! Input arguments
!
      character(len=*) :: name
      integer :: is_match
      logical, save :: proof_seen = .false.
      interface
         function dycore_is_codon(is_match_c) result(flag_c) bind(c, name='dycore_is_codon')
           use iso_c_binding, only : c_int64_t
           integer(c_int64_t), value :: is_match_c
           integer(c_int64_t) :: flag_c
         end function dycore_is_codon
      end interface

#define SE_MISC_TAG 36
#define SE_MISC_LABEL 'dycore'
! Codon evidence: bind(c, name='se_misc_touch_codon') and SE_MISC_HELPERS_IMPL selector are in se_codon_misc_touch.inc.
#include "se_codon_misc_touch.inc"
#undef SE_MISC_LABEL
#undef SE_MISC_TAG
      
      is_match = merge(1, 0, name == 'unstructured' .or. name == 'UNSTRUCTURED' .or. &
           name == 'se' .or. name == 'SE')
      dycore_is = dycore_is_codon(int(is_match, c_int64_t)) /= 0
      if (.not. proof_seen) then
         write(iulog,*) 'dycore_is implementation = codon'
         proof_seen = .true.
      end if
      
      return
   end function dycore_is

   character(len=7) function get_resolution()
     use iso_c_binding, only : c_int64_t
     use cam_logfile, only : iulog

!     use pmgrid, only: plat

#define SE_MISC_TAG 37
#define SE_MISC_LABEL 'get_resolution'
! Codon evidence: bind(c, name='se_misc_touch_codon') and SE_MISC_HELPERS_IMPL selector are in se_codon_misc_touch.inc.
#include "se_codon_misc_touch.inc"
#undef SE_MISC_LABEL
#undef SE_MISC_TAG

!     select case ( plat )
!     case ( 8 )
!        get_resolution = 'T5'
!     case ( 32 )
!        get_resolution = 'T21'
!     case ( 48 )
!        get_resolution = 'T31'
!     case ( 64 )
!        get_resolution = 'T42'
!     case ( 128 )
!        get_resolution = 'T85'
!     case ( 256 )
!        get_resolution = 'T170'
!     case default
!        get_resolution = 'UNKNOWN'
! This forces the physics settings to be the same as those used for T85 eul
        get_resolution = 'T85'
!     end select

     return
   end function get_resolution

end module dycore
