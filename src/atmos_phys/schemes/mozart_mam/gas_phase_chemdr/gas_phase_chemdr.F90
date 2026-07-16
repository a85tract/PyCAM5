! Standalone pp_trop_mam3 gas-phase solver orchestration.  CAM registration,
! state extraction, history, timers, and fatal handling remain in
! mo_gas_phase_chemdr.F90.
module gas_phase_chemdr_scheme

  use shr_kind_mod, only : r8 => shr_kind_r8
  use ap_exp_sol_scheme, only : exp_sol_run
  use ap_imp_sol_trop_mam3_scheme, only : imp_sol_run
  use ap_trop_mam3_gas_phase_kernels, only : &
       trop_mam3_set_thermal_rates, trop_mam3_adjust_rates, &
       trop_mam3_adjust_photolysis, trop_mam3_apply_o3s,   &
       trop_mam3_h2so4_before, trop_mam3_h2so4_after

  implicit none
  private

  public :: gas_phase_chemdr_run
  public :: gas_phase_set_thermal_rates_run
  public :: gas_phase_adjust_rates_run
  public :: gas_phase_adjust_photolysis_run

contains

  !> \section arg_table_gas_phase_chemdr_run Argument Table
  !! \htmlinclude gas_phase_chemdr_run.html
  subroutine gas_phase_chemdr_run(                                &
       phase, ncol, pver, gas_pcnst, rxntot, extcnt, clscnt1, clscnt4, &
       nzcnt, itermax, o3_index, o3s_index, h2so4_index,        &
       timestep, total_density, ltrop, clsmap4, permute4,       &
       class_independent_reaction_count, epsilon, factor, small,&
       reaction_rates, het_rates, extfrc, vmr, o3s_loss,       &
       del_h2so4_gasprod, prod_out, loss_out,                  &
       explicit_entry_count, implicit_entry_count,             &
       solver_failure_count, last_failure_dt, errmsg, errflg)

    integer, intent(in) :: phase, ncol, pver, gas_pcnst, rxntot, extcnt
    integer, intent(in) :: clscnt1, clscnt4, nzcnt, itermax
    integer, intent(in) :: o3_index, o3s_index, h2so4_index
    real(r8), intent(in) :: timestep
    real(r8), intent(in) :: total_density(ncol,pver)
    integer, intent(in) :: ltrop(ncol)
    integer, intent(in) :: clsmap4(clscnt4), permute4(clscnt4)
    integer, intent(in) :: class_independent_reaction_count
    real(r8), intent(in) :: epsilon(clscnt4)
    logical, intent(in) :: factor(itermax)
    real(r8), intent(in) :: small
    real(r8), intent(in) :: reaction_rates(ncol,pver,rxntot)
    real(r8), intent(in) :: het_rates(ncol,pver,gas_pcnst)
    real(r8), intent(in) :: extfrc(ncol,pver,extcnt)
    real(r8), intent(inout) :: vmr(ncol,pver,gas_pcnst)
    real(r8), intent(out) :: o3s_loss(ncol,pver)
    real(r8), intent(inout) :: del_h2so4_gasprod(ncol,pver)
    real(r8), intent(out) :: prod_out(ncol,pver,clscnt4)
    real(r8), intent(out) :: loss_out(ncol,pver,clscnt4)
    integer, intent(inout) :: explicit_entry_count
    integer, intent(inout) :: implicit_entry_count
    integer, intent(out) :: solver_failure_count
    real(r8), intent(out) :: last_failure_dt
    character(len=512), intent(out) :: errmsg
    integer, intent(out) :: errflg

    character(len=512) :: child_errmsg
    integer :: child_errflg

    errmsg = ''
    errflg = 0
    o3s_loss = 0._r8
    prod_out = 0._r8
    loss_out = 0._r8
    solver_failure_count = 0
    last_failure_dt = 0._r8

    if (phase == 1) then
       call trop_mam3_h2so4_before(ncol, pver, gas_pcnst, h2so4_index, &
            vmr, del_h2so4_gasprod)
       call exp_sol_run(clscnt1, explicit_entry_count, child_errmsg, &
            child_errflg)
       if (child_errflg /= 0) then
          errmsg = child_errmsg
          errflg = child_errflg
       end if
       return
    else if (phase /= 2) then
       errmsg = 'gas_phase_chemdr_run phase must be 1 or 2'
       errflg = 1
       return
    end if

    call imp_sol_run(ncol, pver, gas_pcnst, rxntot, extcnt, clscnt4, &
         nzcnt, itermax, vmr, reaction_rates, het_rates, extfrc,    &
         timestep, total_density, ltrop, clsmap4, permute4,        &
         class_independent_reaction_count, epsilon, factor, small, &
         prod_out, loss_out, solver_failure_count, last_failure_dt,&
         implicit_entry_count, child_errmsg, child_errflg)
    if (child_errflg /= 0) then
       errmsg = child_errmsg
       errflg = child_errflg
       prod_out = 0._r8
       loss_out = 0._r8
       return
    end if

    call trop_mam3_apply_o3s(ncol, pver, gas_pcnst, o3_index, o3s_index, &
         ltrop, timestep, o3s_loss, vmr)
    call trop_mam3_h2so4_after(ncol, pver, gas_pcnst, h2so4_index, vmr, &
         del_h2so4_gasprod)

  end subroutine gas_phase_chemdr_run

  subroutine gas_phase_set_thermal_rates_run(ncol, pver, rxntot,  &
       temperature, reaction_rates, ierr)
    integer, intent(in) :: ncol, pver, rxntot
    real(r8), intent(in) :: temperature(ncol,pver)
    real(r8), intent(inout) :: reaction_rates(ncol,pver,rxntot)
    integer, intent(out) :: ierr

    call trop_mam3_set_thermal_rates(ncol, pver, rxntot,          &
         temperature, reaction_rates, ierr)
  end subroutine gas_phase_set_thermal_rates_run

  subroutine gas_phase_adjust_rates_run(ncol, pver, nfs, rxntot, &
       reaction_rates, invariants, total_density, ierr)
    integer, intent(in) :: ncol, pver, nfs, rxntot
    real(r8), intent(inout) :: reaction_rates(ncol,pver,rxntot)
    real(r8), intent(in) :: invariants(ncol,pver,nfs)
    real(r8), intent(in) :: total_density(ncol,pver)
    integer, intent(out) :: ierr

    call trop_mam3_adjust_rates(ncol, pver, nfs, rxntot,         &
         reaction_rates, invariants, total_density, ierr)
  end subroutine gas_phase_adjust_rates_run

  subroutine gas_phase_adjust_photolysis_run(ncol, pver, nfs,    &
       phtcnt, photolysis_rate, invariants, total_density)
    integer, intent(in) :: ncol, pver, nfs, phtcnt
    real(r8), intent(inout) :: photolysis_rate(ncol,pver,phtcnt)
    real(r8), intent(in) :: invariants(ncol,pver,nfs)
    real(r8), intent(in) :: total_density(ncol,pver)

    call trop_mam3_adjust_photolysis(ncol, pver, nfs, phtcnt,    &
         photolysis_rate, invariants, total_density)
  end subroutine gas_phase_adjust_photolysis_run

end module gas_phase_chemdr_scheme
