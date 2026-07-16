module ap_modal_aero_gasaerexch_scheme

  use shr_kind_mod, only : r8 => shr_kind_r8
  use ap_modal_aero_gasaerexch_kernels, only : &
       modal_aero_gasaerexch_kernel

  implicit none
  private

  public :: modal_aero_gasaerexch_sub_run

contains

  !> \section arg_table_modal_aero_gasaerexch_sub_run Argument Table
  !! \htmlinclude modal_aero_gasaerexch_sub_run.html
  subroutine modal_aero_gasaerexch_sub_run(                    &
       lchnk, ncol, nstep, pcols_in, pver_in, pcnstxx_in,      &
       ntot_amode_in, pcnst_in, top_lev_in, maxspec_in,       &
       ntot_aspectype_in, loffset, deltat,                    &
       l_so4g_in, l_nh4g_in, l_msag_in, l_soag_in,           &
       modefrm_pcage_in, modetoo_pcage_in,                   &
       nspecfrm_pcage_in, lspecfrm_pcage_in,                 &
       lspectoo_pcage_in, modeptr_pcarbon_in,                &
       numptr_amode_in, nspec_amode_in,                      &
       lmassptr_amode_in, lspectype_amode_in,                &
       lptr_so4_a_amode_in, lptr_nh4_a_amode_in,            &
       lptr_soa_a_amode_in, lptr_pom_a_amode_in,            &
       sigmag_amode_in, alnsg_amode_in,                     &
       specmw_amode_in, specdens_amode_in,                  &
       specmw_so4_amode_in, specmw_nh4_amode_in,            &
       specmw_soa_amode_in, specdens_so4_amode_in,          &
       specdens_nh4_amode_in, specdens_soa_amode_in,        &
       dr_so4_monolayers_pcage_in, soa_equivso4_factor_in,  &
       gravity_in, mwdry_in, rair_in, adv_mass_in,          &
       t, pmid, pdel, qh2o, troplev, q, qqcw,               &
       dqdt_other, dqqcwdt_other, dgncur_a, dgncur_awet,    &
       has_sulfeq, sulfeq, qsrflx, qqcwsrflx,               &
       dotend, dotendqqcw, dotendrn, dotendqqcwrn,          &
       ierr_out, rename_error_out)

    integer, intent(in) :: lchnk, ncol, nstep
    integer, intent(in) :: pcols_in, pver_in, pcnstxx_in
    integer, intent(in) :: ntot_amode_in, pcnst_in, top_lev_in
    integer, intent(in) :: maxspec_in, ntot_aspectype_in, loffset
    integer, intent(in) :: l_so4g_in, l_nh4g_in, l_msag_in, l_soag_in
    integer, intent(in) :: modefrm_pcage_in, modetoo_pcage_in
    integer, intent(in) :: nspecfrm_pcage_in, modeptr_pcarbon_in
    integer, intent(in) :: lspecfrm_pcage_in(maxspec_in)
    integer, intent(in) :: lspectoo_pcage_in(maxspec_in)
    integer, intent(in) :: numptr_amode_in(ntot_amode_in)
    integer, intent(in) :: nspec_amode_in(ntot_amode_in)
    integer, intent(in) :: lmassptr_amode_in(maxspec_in,ntot_amode_in)
    integer, intent(in) :: lspectype_amode_in(maxspec_in,ntot_amode_in)
    integer, intent(in) :: lptr_so4_a_amode_in(ntot_amode_in)
    integer, intent(in) :: lptr_nh4_a_amode_in(ntot_amode_in)
    integer, intent(in) :: lptr_soa_a_amode_in(ntot_amode_in)
    integer, intent(in) :: lptr_pom_a_amode_in(ntot_amode_in)
    real(r8), intent(in) :: deltat
    real(r8), intent(in) :: sigmag_amode_in(ntot_amode_in)
    real(r8), intent(in) :: alnsg_amode_in(ntot_amode_in)
    real(r8), intent(in) :: specmw_amode_in(ntot_aspectype_in)
    real(r8), intent(in) :: specdens_amode_in(ntot_aspectype_in)
    real(r8), intent(in) :: specmw_so4_amode_in
    real(r8), intent(in) :: specmw_nh4_amode_in
    real(r8), intent(in) :: specmw_soa_amode_in
    real(r8), intent(in) :: specdens_so4_amode_in
    real(r8), intent(in) :: specdens_nh4_amode_in
    real(r8), intent(in) :: specdens_soa_amode_in
    real(r8), intent(in) :: dr_so4_monolayers_pcage_in
    real(r8), intent(in) :: soa_equivso4_factor_in
    real(r8), intent(in) :: gravity_in, mwdry_in, rair_in
    real(r8), intent(in) :: adv_mass_in(pcnstxx_in)
    real(r8), intent(in) :: t(pcols_in,pver_in)
    real(r8), intent(in) :: pmid(pcols_in,pver_in)
    real(r8), intent(in) :: pdel(pcols_in,pver_in)
    real(r8), intent(in) :: qh2o(pcols_in,pver_in)
    integer, intent(in) :: troplev(pcols_in)
    real(r8), intent(inout) :: q(ncol,pver_in,pcnstxx_in)
    real(r8), intent(inout) :: qqcw(ncol,pver_in,pcnstxx_in)
    real(r8), intent(in) :: dqdt_other(ncol,pver_in,pcnstxx_in)
    real(r8), intent(in) :: dqqcwdt_other(ncol,pver_in,pcnstxx_in)
    real(r8), intent(in) :: dgncur_a(pcols_in,pver_in,ntot_amode_in)
    real(r8), intent(in) :: dgncur_awet(pcols_in,pver_in,ntot_amode_in)
    logical, intent(in) :: has_sulfeq
    real(r8), intent(in) :: sulfeq(pcols_in,pver_in,ntot_amode_in)
    real(r8), intent(out) :: qsrflx(pcols_in,pcnstxx_in,2)
    real(r8), intent(out) :: qqcwsrflx(pcols_in,pcnstxx_in,2)
    logical, intent(out) :: dotend(pcnstxx_in)
    logical, intent(out) :: dotendqqcw(pcnstxx_in)
    logical, intent(out) :: dotendrn(pcnstxx_in)
    logical, intent(out) :: dotendqqcwrn(pcnstxx_in)
    integer, intent(out) :: ierr_out, rename_error_out

    call modal_aero_gasaerexch_kernel(                         &
         lchnk, ncol, nstep, pcols_in, pver_in, pcnstxx_in,   &
         ntot_amode_in, pcnst_in, top_lev_in, maxspec_in,    &
         ntot_aspectype_in, loffset, deltat,                 &
         l_so4g_in, l_nh4g_in, l_msag_in, l_soag_in,        &
         modefrm_pcage_in, modetoo_pcage_in,                &
         nspecfrm_pcage_in, lspecfrm_pcage_in,              &
         lspectoo_pcage_in, modeptr_pcarbon_in,             &
         numptr_amode_in, nspec_amode_in,                   &
         lmassptr_amode_in, lspectype_amode_in,             &
         lptr_so4_a_amode_in, lptr_nh4_a_amode_in,         &
         lptr_soa_a_amode_in, lptr_pom_a_amode_in,         &
         sigmag_amode_in, alnsg_amode_in,                  &
         specmw_amode_in, specdens_amode_in,               &
         specmw_so4_amode_in, specmw_nh4_amode_in,         &
         specmw_soa_amode_in, specdens_so4_amode_in,       &
         specdens_nh4_amode_in, specdens_soa_amode_in,     &
         dr_so4_monolayers_pcage_in, soa_equivso4_factor_in, &
         gravity_in, mwdry_in, rair_in, adv_mass_in,       &
         t, pmid, pdel, qh2o, troplev, q, qqcw,            &
         dqdt_other, dqqcwdt_other, dgncur_a, dgncur_awet, &
         has_sulfeq, sulfeq, qsrflx, qqcwsrflx,            &
         dotend, dotendqqcw, dotendrn, dotendqqcwrn,       &
         ierr_out, rename_error_out)

  end subroutine modal_aero_gasaerexch_sub_run

end module ap_modal_aero_gasaerexch_scheme
