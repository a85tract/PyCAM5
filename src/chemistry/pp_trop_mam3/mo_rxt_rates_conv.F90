module mo_rxt_rates_conv
  use shr_kind_mod, only : r8 => shr_kind_r8
  use cam_logfile,  only : iulog
  use spmd_utils,   only : masterproc
  use iso_c_binding, only : c_int64_t, c_loc, c_ptr
  implicit none
  private
   public :: set_rates
   logical, save :: set_rates_logged = .false.
contains
   logical function set_rates_use_codon()
      character(len=32) :: impl_name
      integer :: status, n, i, code

      impl_name = 'native'
      call cam_codon_get_impl('SET_RATES_IMPL', impl_name, n, status)
      if (status == 0 .and. n > 0) then
         do i = 1, n
            code = iachar(impl_name(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
               impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
         end do
         set_rates_use_codon = trim(adjustl(impl_name(:n))) == 'codon'
      else
         set_rates_use_codon = .false.
      end if
   end function set_rates_use_codon

   subroutine set_rates( rxt_rates, sol, ncol )
      real(r8), target, intent(inout) :: rxt_rates(:,:,:)
      real(r8), target, intent(in) :: sol(:,:,:)
      integer, intent(in) :: ncol
      interface
         subroutine set_rates_codon(ncol_c, rxt_d1_c, rxt_d2_c, sol_d1_c, sol_d2_c, rxt_rates_p, sol_p) &
              bind(c, name="set_rates_codon")
            use iso_c_binding, only : c_int64_t, c_ptr
            integer(c_int64_t), value :: ncol_c, rxt_d1_c, rxt_d2_c, sol_d1_c, sol_d2_c
            type(c_ptr), value :: rxt_rates_p, sol_p
         end subroutine set_rates_codon
      end interface

      if (set_rates_use_codon()) then
         call set_rates_codon(int(ncol, c_int64_t), int(size(rxt_rates,1), c_int64_t), &
              int(size(rxt_rates,2), c_int64_t), int(size(sol,1), c_int64_t), int(size(sol,2), c_int64_t), &
              c_loc(rxt_rates(1,1,1)), c_loc(sol(1,1,1)))
         if (masterproc .and. .not. set_rates_logged) then
            write(iulog,'(A)') 'set_rates direct = codon'
            call flush(iulog)
            set_rates_logged = .true.
         end if
         return
      end if

      rxt_rates(:ncol,:,     1) = rxt_rates(:ncol,:,     1)*sol(:ncol,:,     1)                                                ! rate_const*H2O2
                                                                                                                               ! rate_const
      rxt_rates(:ncol,:,     3) = rxt_rates(:ncol,:,     3)*sol(:ncol,:,     1)                                                ! rate_const*OH*H2O2
      rxt_rates(:ncol,:,     4) = rxt_rates(:ncol,:,     4)*sol(:ncol,:,     3)                                                ! rate_const*OH*SO2
      rxt_rates(:ncol,:,     5) = rxt_rates(:ncol,:,     5)*sol(:ncol,:,     4)                                                ! rate_const*OH*DMS
      rxt_rates(:ncol,:,     6) = rxt_rates(:ncol,:,     6)*sol(:ncol,:,     4)                                                ! rate_const*OH*DMS
      rxt_rates(:ncol,:,     7) = rxt_rates(:ncol,:,     7)*sol(:ncol,:,     4)                                                ! rate_const*NO3*DMS
  end subroutine set_rates
end module mo_rxt_rates_conv
