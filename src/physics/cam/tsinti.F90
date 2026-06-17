
subroutine tsinti (tmeltx, latvapx, rairx, stebolx, laticex)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize surface temperature calculation constants
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: L. Buja
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use cam_control_mod, only: rair, stebol, snwedp, snwedp, latice, tmelt, latvap
   use cam_logfile, only: iulog
   use iso_c_binding, only: c_double, c_loc, c_ptr
   use spmd_utils, only: masterproc
   implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   real(r8) tmeltx         ! Melting temperature of snow and ice
   real(r8) latvapx        ! Latent heat of vaporization
   real(r8) rairx          ! Gas constant for dry air
   real(r8) stebolx        ! Stefan-Boltzmann constant
   real(r8) laticex        ! latent heat of fusion
!
   logical, save :: use_native_tsinti_impl = .false.
   logical, save :: tsinti_impl_selected = .false.
   logical, save :: tsinti_proof_written = .false.
   real(c_double), target :: tsinti_values(6)

   interface
      subroutine tsinti_codon(tmelt_c, latvap_c, rair_c, stebol_c, latice_c, values_p) &
           bind(c, name="tsinti_codon")
         use iso_c_binding, only: c_double, c_ptr
         real(c_double), value :: tmelt_c, latvap_c, rair_c, stebol_c, latice_c
         type(c_ptr), value :: values_p
      end subroutine tsinti_codon
      function tsinti_param_codon(value_c) result(out_c) bind(c, name="tsinti_param_codon")
         use iso_c_binding, only: c_double
         real(c_double), value :: value_c
         real(c_double) :: out_c
      end function tsinti_param_codon
   end interface
!
!-----------------------------------------------------------------------
!
   call tsinti_select_impl()
   if (use_native_tsinti_impl) then
      latice = laticex    ! Latent heat of fusion at 0'C = 3.336e5 J/Kg
      tmelt  = tmeltx
      latvap = latvapx
      rair   = rairx
      stebol = stebolx
      snwedp = 10.0_r8       ! 10:1 Snow:water equivalent depth factor
      return
   end if

   call tsinti_proof_once()
   call tsinti_codon(real(tmeltx, c_double), real(latvapx, c_double), real(rairx, c_double), &
        real(stebolx, c_double), real(laticex, c_double), c_loc(tsinti_values))
   latice = real(tsinti_values(1), r8)    ! Latent heat of fusion at 0'C = 3.336e5 J/Kg
   tmelt  = real(tsinti_values(2), r8)
   latvap = real(tsinti_values(3), r8)
   rair   = real(tsinti_values(4), r8)
   stebol = real(tsinti_values(5), r8)
   snwedp = real(tsinti_values(6), r8)       ! 10:1 Snow:water equivalent depth factor
!
   return

contains

subroutine tsinti_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (tsinti_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('TSINTI_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_tsinti_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_tsinti_impl = .false.
   end if

   tsinti_impl_selected = .true.

   if (masterproc) then
      if (use_native_tsinti_impl) then
         write(iulog,*) 'tsinti implementation = native'
      else
         write(iulog,*) 'tsinti implementation = codon'
      end if
   end if

end subroutine tsinti_select_impl

subroutine tsinti_proof_once()

   if (tsinti_proof_written) return
   tsinti_proof_written = .true.

   if (masterproc) then
      write(iulog,'(A)') 'tsinti entered (surface constants scalar helpers = codon)'
   end if

end subroutine tsinti_proof_once

real(r8) function tsinti_param(value_in)

   real(r8), intent(in) :: value_in
   real(c_double) :: out_c

   call tsinti_select_impl()

   if (use_native_tsinti_impl) then
      tsinti_param = value_in
      return
   end if

   call tsinti_proof_once()
   out_c = tsinti_param_codon(real(value_in, c_double))
   tsinti_param = real(out_c, r8)

end function tsinti_param

end subroutine tsinti
