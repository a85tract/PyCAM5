subroutine zenith(calday  ,clat    , clon   ,coszrs  ,ncol    )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute cosine of solar zenith angle for albedo and radiation
!   computations.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J. Kiehl
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use shr_const_mod, only: orb_pi => shr_const_pi
   use shr_orb_mod
   use cam_control_mod, only: lambm0, obliqr, eccen, mvelpp
   use cam_logfile, only: iulog
   use spmd_utils, only: masterproc
   use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
   implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: ncol                 ! number of positions
   real(r8), intent(in) :: calday              ! Calendar day, including fraction
   real(r8), target, intent(in) :: clat(ncol)  ! Current centered latitude (radians)
   real(r8), target, intent(in) :: clon(ncol)  ! Centered longitude (radians)
!
! Output arguments
!
   real(r8), target, intent(out) :: coszrs(ncol) ! Cosine solar zenith angle
!
!---------------------------Local variables-----------------------------
!
   integer i         ! Position loop index
   real(r8) delta    ! Solar declination angle  in radians
   real(r8) eccf     ! Earth orbit eccentricity factor
   logical, save :: zenith_use_native_impl = .false.
   logical, save :: zenith_impl_selected = .false.

   interface
      subroutine zenith_codon(ncol_c, calday_c, pi_c, delta_c, clat_p, clon_p, coszrs_p) bind(c, name="zenith_codon")
         use iso_c_binding, only : c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c
         real(c_double), value :: calday_c, pi_c, delta_c
         type(c_ptr), value :: clat_p, clon_p, coszrs_p
      end subroutine zenith_codon
   end interface
!
!-----------------------------------------------------------------------
!
   call zenith_select_impl()

   call shr_orb_decl (calday  ,eccen     ,mvelpp  ,lambm0  ,obliqr  , &
                      delta   ,eccf      )
!
! Compute local cosine solar zenith angle,
!
   if (.not. zenith_use_native_impl) then
      call zenith_codon(int(ncol, c_int64_t), calday, orb_pi, delta, c_loc(clat), c_loc(clon), c_loc(coszrs))
      return
   end if

   do i=1,ncol
      coszrs(i) = shr_orb_cosz( calday, clat(i), clon(i), delta )
   end do

contains

subroutine zenith_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (zenith_impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('ZENITH_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      zenith_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      zenith_use_native_impl = .false.
   end if

   zenith_impl_selected = .true.

   if (masterproc) then
      if (zenith_use_native_impl) then
         write(iulog,*) 'zenith implementation = native'
      else
         write(iulog,*) 'zenith implementation = codon'
      end if
      call flush(iulog)
   end if

end subroutine zenith_select_impl

end subroutine zenith
