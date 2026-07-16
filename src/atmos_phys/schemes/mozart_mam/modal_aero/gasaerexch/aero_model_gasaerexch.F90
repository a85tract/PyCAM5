! Standalone orchestration for modal aerosol gas exchange, nucleation, and
! coagulation.  CAM state extraction, setsox, history, timers, and errors are
! intentionally retained by the host facade in aero_model.F90.
module ap_aero_model_gasaerexch_scheme

  use shr_kind_mod, only : r8 => shr_kind_r8
  use ap_modal_aero_gasaerexch_scheme, only : &
       modal_aero_gasaerexch_sub_run
  use ap_modal_aero_newnuc_scheme, only : modal_aero_newnuc_sub_run
  use ap_modal_aero_coag_scheme, only : modal_aero_coag_sub_run

  implicit none
  private

  public :: aero_model_gasaerexch_run

contains

  !> \section arg_table_aero_model_gasaerexch_run Argument Table
  !! \htmlinclude aero_model_gasaerexch_run.html
  subroutine aero_model_gasaerexch_run(                            &
       lchnk, ncol, nstep, pcols_in, pver_in, pcnstxx_in,         &
       pcnst_in, ntot_amode_in, top_lev_in, maxspec_in,          &
       ntot_aspectype_in, maxpair_in, loffset, deltat,           &
       l_so4g_in, l_nh4g_in, l_msag_in, l_soag_in,              &
       modefrm_pcage_in, modetoo_pcage_in, nspecfrm_pcage_in,   &
       lspecfrm_pcage_in, lspectoo_pcage_in,                    &
       modeptr_accum_in, modeptr_aitken_in,                     &
       modeptr_pcarbon_in, numptr_amode_in,                     &
       mprognum_amode_in, nspec_amode_in,                       &
       lmassptr_amode_in, lspectype_amode_in,                   &
       lptr_so4_a_amode_in, lptr_nh4_a_amode_in,               &
       lptr_soa_a_amode_in, lptr_pom_a_amode_in,               &
       sigmag_amode_in, alnsg_amode_in, specmw_amode_in,       &
       specdens_amode_in, specmw_so4_amode_in,                 &
       specmw_nh4_amode_in, specmw_soa_amode_in,               &
       specdens_so4_amode_in, specdens_nh4_amode_in,           &
       specdens_soa_amode_in, dr_so4_monolayers_pcage_in,      &
       n_so4_monolayers_pcage_in, soa_equivso4_factor_in,      &
       pair_option_acoag_in, npair_acoag_in,                   &
       modefrm_acoag_in, modetoo_acoag_in,                     &
       nspecfrm_acoag_in, lspecfrm_acoag_in,                   &
       lspectoo_acoag_in, l_h2so4_sv_in, l_nh3_sv_in,         &
       lnumait_sv_in, lnh4ait_sv_in, lso4ait_sv_in,           &
       lptr_nh4_aitken_in, dgnumlo_aitken_in,                 &
       dgnum_aitken_in, dgnumhi_aitken_in, gravity_in,        &
       mwdry_in, rair_in, tmelt_in, adv_mass_in,              &
       t, pmid, pdel, mbar, zm, pblh, qh2o, qv_sat, cld,     &
       troplev, vmr0, vmr_before_setsox,                      &
       vmrcw_before_setsox, vmr, vmrcw,                       &
       del_h2so4_gasprod, dgncur_a, dgncur_awet, wetdens_a,  &
       has_sulfeq, sulfeq, gs_flux, aq_flux,                  &
       ga_qsrflx, ga_qqcwsrflx, ga_dotend, ga_dotendqqcw,     &
       ga_dotendrn, ga_dotendqqcwrn, nn_qsrflx, nn_dotend,   &
       co_qsrflx, co_dotend, ga_ierr, rename_ierr, co_ierr)

    integer, intent(in) :: lchnk, ncol, nstep
    integer, intent(in) :: pcols_in, pver_in, pcnstxx_in
    integer, intent(in) :: pcnst_in, ntot_amode_in, top_lev_in
    integer, intent(in) :: maxspec_in, ntot_aspectype_in, maxpair_in
    integer, intent(in) :: loffset
    integer, intent(in) :: l_so4g_in, l_nh4g_in, l_msag_in, l_soag_in
    integer, intent(in) :: modefrm_pcage_in, modetoo_pcage_in
    integer, intent(in) :: nspecfrm_pcage_in
    integer, intent(in) :: lspecfrm_pcage_in(maxspec_in)
    integer, intent(in) :: lspectoo_pcage_in(maxspec_in)
    integer, intent(in) :: modeptr_accum_in, modeptr_aitken_in
    integer, intent(in) :: modeptr_pcarbon_in
    integer, intent(in) :: numptr_amode_in(ntot_amode_in)
    integer, intent(in) :: mprognum_amode_in(ntot_amode_in)
    integer, intent(in) :: nspec_amode_in(ntot_amode_in)
    integer, intent(in) :: lmassptr_amode_in(maxspec_in,ntot_amode_in)
    integer, intent(in) :: lspectype_amode_in(maxspec_in,ntot_amode_in)
    integer, intent(in) :: lptr_so4_a_amode_in(ntot_amode_in)
    integer, intent(in) :: lptr_nh4_a_amode_in(ntot_amode_in)
    integer, intent(in) :: lptr_soa_a_amode_in(ntot_amode_in)
    integer, intent(in) :: lptr_pom_a_amode_in(ntot_amode_in)
    integer, intent(in) :: pair_option_acoag_in, npair_acoag_in
    integer, intent(in) :: modefrm_acoag_in(maxpair_in)
    integer, intent(in) :: modetoo_acoag_in(maxpair_in)
    integer, intent(in) :: nspecfrm_acoag_in(maxpair_in)
    integer, intent(in) :: lspecfrm_acoag_in(maxspec_in,maxpair_in)
    integer, intent(in) :: lspectoo_acoag_in(maxspec_in,maxpair_in)
    integer, intent(in) :: l_h2so4_sv_in, l_nh3_sv_in
    integer, intent(in) :: lnumait_sv_in, lnh4ait_sv_in, lso4ait_sv_in
    integer, intent(in) :: lptr_nh4_aitken_in

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
    real(r8), intent(in) :: n_so4_monolayers_pcage_in
    real(r8), intent(in) :: soa_equivso4_factor_in
    real(r8), intent(in) :: dgnumlo_aitken_in
    real(r8), intent(in) :: dgnum_aitken_in
    real(r8), intent(in) :: dgnumhi_aitken_in
    real(r8), intent(in) :: gravity_in, mwdry_in, rair_in, tmelt_in
    real(r8), intent(in) :: adv_mass_in(pcnstxx_in)

    real(r8), intent(in) :: t(pcols_in,pver_in)
    real(r8), intent(in) :: pmid(pcols_in,pver_in)
    real(r8), intent(in) :: pdel(pcols_in,pver_in)
    real(r8), intent(in) :: mbar(pcols_in,pver_in)
    real(r8), intent(in) :: zm(pcols_in,pver_in)
    real(r8), intent(in) :: pblh(pcols_in)
    real(r8), intent(in) :: qh2o(pcols_in,pver_in)
    real(r8), intent(in) :: qv_sat(pcols_in,pver_in)
    real(r8), intent(in) :: cld(ncol,pver_in)
    integer, intent(in) :: troplev(pcols_in)
    real(r8), intent(in) :: vmr0(ncol,pver_in,pcnstxx_in)
    real(r8), intent(in) :: vmr_before_setsox(ncol,pver_in,pcnstxx_in)
    real(r8), intent(in) :: vmrcw_before_setsox(ncol,pver_in,pcnstxx_in)
    real(r8), intent(inout) :: vmr(ncol,pver_in,pcnstxx_in)
    real(r8), intent(inout) :: vmrcw(ncol,pver_in,pcnstxx_in)
    real(r8), intent(in) :: del_h2so4_gasprod(ncol,pver_in)
    real(r8), intent(in) :: dgncur_a(pcols_in,pver_in,ntot_amode_in)
    real(r8), intent(in) :: dgncur_awet(pcols_in,pver_in,ntot_amode_in)
    real(r8), intent(in) :: wetdens_a(pcols_in,pver_in,ntot_amode_in)
    logical, intent(in) :: has_sulfeq
    real(r8), intent(in) :: sulfeq(pcols_in,pver_in,ntot_amode_in)

    real(r8), intent(out) :: gs_flux(pcols_in,pcnstxx_in)
    real(r8), intent(out) :: aq_flux(pcols_in,pcnstxx_in)
    real(r8), intent(out) :: ga_qsrflx(pcols_in,pcnstxx_in,2)
    real(r8), intent(out) :: ga_qqcwsrflx(pcols_in,pcnstxx_in,2)
    logical, intent(out) :: ga_dotend(pcnstxx_in)
    logical, intent(out) :: ga_dotendqqcw(pcnstxx_in)
    logical, intent(out) :: ga_dotendrn(pcnstxx_in)
    logical, intent(out) :: ga_dotendqqcwrn(pcnstxx_in)
    real(r8), intent(out) :: nn_qsrflx(pcols_in,pcnst_in)
    logical, intent(out) :: nn_dotend(pcnst_in)
    real(r8), intent(out) :: co_qsrflx(pcols_in,pcnst_in)
    logical, intent(out) :: co_dotend(pcnst_in)
    integer, intent(out) :: ga_ierr, rename_ierr, co_ierr

    integer :: k, m
    real(r8) :: gs_dqdt(ncol,pver_in,pcnstxx_in)
    real(r8) :: dqdt_other(ncol,pver_in,pcnstxx_in)
    real(r8) :: dqqcwdt_other(ncol,pver_in,pcnstxx_in)
    real(r8) :: del_h2so4_aeruptk(ncol,pver_in)

    gs_flux = 0.0_r8
    aq_flux = 0.0_r8
    gs_dqdt = (vmr_before_setsox - vmr0) / deltat
    dqdt_other = (vmr - vmr_before_setsox) / deltat
    dqqcwdt_other = (vmrcw - vmrcw_before_setsox) / deltat

    do m = 1, pcnstxx_in
       do k = 1, pver_in
          gs_flux(1:ncol,m) = gs_flux(1:ncol,m) +                &
               gs_dqdt(:,k,m) * adv_mass_in(m) /                &
               mbar(1:ncol,k) *                                 &
               pdel(1:ncol,k) / gravity_in
          aq_flux(1:ncol,m) = aq_flux(1:ncol,m) +               &
               dqdt_other(:,k,m) * adv_mass_in(m) /             &
               mbar(1:ncol,k) * pdel(1:ncol,k) / gravity_in
       end do
    end do

    if (l_so4g_in > 0 .and. l_so4g_in <= pcnstxx_in) then
       del_h2so4_aeruptk = vmr(:,:,l_so4g_in)
    else
       del_h2so4_aeruptk = 0.0_r8
    end if

    call modal_aero_gasaerexch_sub_run(                           &
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
         t, pmid, pdel, qh2o, troplev, vmr, vmrcw,            &
         dqdt_other, dqqcwdt_other, dgncur_a, dgncur_awet,    &
         has_sulfeq, sulfeq, ga_qsrflx, ga_qqcwsrflx,         &
         ga_dotend, ga_dotendqqcw, ga_dotendrn,               &
         ga_dotendqqcwrn, ga_ierr, rename_ierr)

    if (l_so4g_in > 0 .and. l_so4g_in <= pcnstxx_in) then
       del_h2so4_aeruptk = vmr(:,:,l_so4g_in) - del_h2so4_aeruptk
    end if

    call modal_aero_newnuc_sub_run(                               &
         lchnk, ncol, nstep, pcols_in, pver_in, pcnstxx_in,      &
         pcnst_in, top_lev_in, loffset, deltat,                 &
         l_h2so4_sv_in, l_nh3_sv_in, lnumait_sv_in,             &
         lnh4ait_sv_in, lso4ait_sv_in, lptr_nh4_aitken_in,     &
         dgnumlo_aitken_in, dgnum_aitken_in,                   &
         dgnumhi_aitken_in, specdens_so4_amode_in,             &
         specmw_so4_amode_in, specmw_nh4_amode_in,             &
         gravity_in, t, pmid, pdel, zm, pblh, qh2o, qv_sat,   &
         cld, vmr, del_h2so4_gasprod, del_h2so4_aeruptk,      &
         nn_qsrflx, nn_dotend)

    call modal_aero_coag_sub_run(                                 &
         lchnk, ncol, nstep, pcols_in, pver_in, pcnstxx_in,      &
         ntot_amode_in, pcnst_in, top_lev_in, maxspec_in,       &
         ntot_aspectype_in, maxpair_in, loffset, deltat,        &
         pair_option_acoag_in, npair_acoag_in,                  &
         modefrm_acoag_in, modetoo_acoag_in,                   &
         nspecfrm_acoag_in, lspecfrm_acoag_in,                 &
         lspectoo_acoag_in, modeptr_accum_in,                  &
         modeptr_aitken_in, modeptr_pcarbon_in,                &
         numptr_amode_in, mprognum_amode_in,                   &
         nspec_amode_in, lmassptr_amode_in,                    &
         lspectype_amode_in, sigmag_amode_in,                  &
         alnsg_amode_in, lptr_so4_a_amode_in,                  &
         lptr_nh4_a_amode_in, lptr_soa_a_amode_in,            &
         specmw_amode_in, specdens_amode_in,                   &
         specmw_so4_amode_in, specdens_so4_amode_in,          &
         specmw_nh4_amode_in, specdens_nh4_amode_in,          &
         specmw_soa_amode_in, specdens_soa_amode_in,          &
         n_so4_monolayers_pcage_in, soa_equivso4_factor_in,   &
         tmelt_in, t, pmid, pdel, vmr, dgncur_a,              &
         dgncur_awet, wetdens_a, co_qsrflx, co_dotend, co_ierr)

  end subroutine aero_model_gasaerexch_run

end module ap_aero_model_gasaerexch_scheme
