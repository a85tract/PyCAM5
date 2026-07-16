! Runtime-only pp_trop_mam3 mechanism operations.  Generated mechanism
! dimensions and tables are supplied explicitly by the host facade.
module ap_trop_mam3_gas_phase_kernels

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private

  public :: trop_mam3_set_thermal_rates
  public :: trop_mam3_adjust_rates
  public :: trop_mam3_adjust_photolysis
  public :: trop_mam3_apply_o3s
  public :: trop_mam3_h2so4_before
  public :: trop_mam3_h2so4_after

contains

  subroutine trop_mam3_set_thermal_rates(ncol, pver, rxntot, temp, rate, ierr)
    integer, intent(in) :: ncol, pver, rxntot
    real(r8), intent(in) :: temp(ncol,pver)
    real(r8), intent(inout) :: rate(ncol,pver,rxntot)
    integer, intent(out) :: ierr
    real(r8) :: inverse_temperature(ncol,pver)

    ierr = 0
    if (rxntot < 7) then
       ierr = 1
       return
    end if
    inverse_temperature = 1._r8 / temp
    rate(:,:,3) = 2.9e-12_r8 * exp(-160._r8 * inverse_temperature)
    rate(:,:,5) = 9.6e-12_r8 * exp(-234._r8 * inverse_temperature)
    rate(:,:,7) = 1.9e-13_r8 * exp(520._r8 * inverse_temperature)
  end subroutine trop_mam3_set_thermal_rates

  subroutine trop_mam3_adjust_rates(ncol, pver, nfs, rxntot, rate, inv, &
       total_density, ierr)
    integer, intent(in) :: ncol, pver, nfs, rxntot
    real(r8), intent(inout) :: rate(ncol,pver,rxntot)
    real(r8), intent(in) :: inv(ncol,pver,nfs)
    real(r8), intent(in) :: total_density(ncol,pver)
    integer, intent(out) :: ierr
    real(r8) :: inverse_density(ncol,pver)

    ierr = 0
    if (nfs < 8 .or. rxntot < 7) then
       ierr = 1
       return
    end if
    rate(:,:,3) = rate(:,:,3) * inv(:,:,6)
    rate(:,:,4) = rate(:,:,4) * inv(:,:,6)
    rate(:,:,5) = rate(:,:,5) * inv(:,:,6)
    rate(:,:,6) = rate(:,:,6) * inv(:,:,6)
    rate(:,:,7) = rate(:,:,7) * inv(:,:,7)
    inverse_density = 1._r8 / total_density
    rate(:,:,2) = rate(:,:,2) * inv(:,:,8) * inv(:,:,8) * inverse_density
  end subroutine trop_mam3_adjust_rates

  subroutine trop_mam3_adjust_photolysis(ncol, pver, nfs, phtcnt,      &
       photolysis_rate, inv, total_density)
    integer, intent(in) :: ncol, pver, nfs, phtcnt
    real(r8), intent(inout) :: photolysis_rate(ncol,pver,phtcnt)
    real(r8), intent(in) :: inv(ncol,pver,nfs)
    real(r8), intent(in) :: total_density(ncol,pver)

    ! The generated pp_trop_mam3 phtadj body is intentionally empty.
    if (ncol < 0 .or. pver < 0 .or. nfs < 0 .or. phtcnt < 0) return
    if (size(photolysis_rate) < 0 .or. size(inv) < 0 .or. &
        size(total_density) < 0) return
  end subroutine trop_mam3_adjust_photolysis

  subroutine trop_mam3_apply_o3s(ncol, pver, gas_pcnst, o3_index,    &
       o3s_index, troplev, timestep, o3s_loss, vmr)
    integer, intent(in) :: ncol, pver, gas_pcnst, o3_index, o3s_index
    integer, intent(in) :: troplev(ncol)
    real(r8), intent(in) :: timestep
    real(r8), intent(in) :: o3s_loss(ncol,pver)
    real(r8), intent(inout) :: vmr(ncol,pver,gas_pcnst)
    integer :: i

    if (o3_index > 0 .and. o3s_index > 0) then
       do i = 1, ncol
          vmr(i,1:troplev(i),o3s_index) = vmr(i,1:troplev(i),o3_index)
          vmr(i,troplev(i)+1:pver,o3s_index) =                       &
               vmr(i,troplev(i)+1:pver,o3s_index) *                 &
               exp(-timestep*o3s_loss(i,troplev(i)+1:pver))
       end do
    end if
  end subroutine trop_mam3_apply_o3s

  subroutine trop_mam3_h2so4_before(ncol, pver, gas_pcnst,           &
       h2so4_index, vmr, increment)
    integer, intent(in) :: ncol, pver, gas_pcnst, h2so4_index
    real(r8), intent(in) :: vmr(ncol,pver,gas_pcnst)
    real(r8), intent(out) :: increment(ncol,pver)

    if (h2so4_index > 0) then
       increment = vmr(:,:,h2so4_index)
    else
       increment = 0._r8
    end if
  end subroutine trop_mam3_h2so4_before

  subroutine trop_mam3_h2so4_after(ncol, pver, gas_pcnst,            &
       h2so4_index, vmr, increment)
    integer, intent(in) :: ncol, pver, gas_pcnst, h2so4_index
    real(r8), intent(in) :: vmr(ncol,pver,gas_pcnst)
    real(r8), intent(inout) :: increment(ncol,pver)

    if (h2so4_index > 0) increment = vmr(:,:,h2so4_index) - increment
  end subroutine trop_mam3_h2so4_after

end module ap_trop_mam3_gas_phase_kernels
