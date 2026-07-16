module ap_modal_aero_coag_scheme

  use shr_kind_mod, only : r8 => shr_kind_r8
  use ap_modal_aero_coag_kernels, only : modal_aero_coag_kernel

  implicit none
  private

  public :: modal_aero_coag_sub_run

contains

  !> \section arg_table_modal_aero_coag_sub_run Argument Table
  !! \htmlinclude modal_aero_coag_sub_run.html
  !!
  subroutine modal_aero_coag_sub_run(                             &
       lchnk, ncol, nstep, pcols_in, pver_in, pcnstxx_in,        &
       ntot_amode_in, pcnst_in, top_lev_in, maxspec_in,         &
       ntot_aspectype_in, maxpair_in, loffset, deltat_main,     &
       pair_option_in,                                          &
       npair_in, modefrm_in, modetoo_in, nspecfrm_in,           &
       lspecfrm_in, lspectoo_in, modeptr_accum_in,              &
       modeptr_aitken_in, modeptr_pcarbon_in, numptr_amode_in,  &
       mprognum_amode_in, nspec_amode_in, lmassptr_amode_in,    &
       lspectype_amode_in, sigmag_amode_in, alnsg_amode_in,     &
       lptr_so4_a_amode_in, lptr_nh4_a_amode_in,                &
       lptr_soa_a_amode_in, specmw_amode_in, specdens_amode_in, &
       specmw_so4_amode_in, specdens_so4_amode_in,              &
       specmw_nh4_amode_in, specdens_nh4_amode_in,              &
       specmw_soa_amode_in, specdens_soa_amode_in,              &
       n_so4_monolayers_pcage_in, soa_equivso4_factor_in,       &
       tmelt_in, t, pmid, pdel, q, dgncur_a, dgncur_awet,      &
       wetdens_a,                                               &
       qsrflx_out, dotend_out, ierr_out)

    integer, intent(in) :: lchnk, ncol, nstep
    integer, intent(in) :: pcols_in, pver_in, pcnstxx_in
    integer, intent(in) :: ntot_amode_in, pcnst_in, top_lev_in
    integer, intent(in) :: maxspec_in, ntot_aspectype_in
    integer, intent(in) :: maxpair_in, loffset
    integer, intent(in) :: pair_option_in, npair_in
    integer, intent(in) :: modefrm_in(maxpair_in)
    integer, intent(in) :: modetoo_in(maxpair_in)
    integer, intent(in) :: nspecfrm_in(maxpair_in)
    integer, intent(in) :: lspecfrm_in(maxspec_in,maxpair_in)
    integer, intent(in) :: lspectoo_in(maxspec_in,maxpair_in)
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

    real(r8), intent(in) :: deltat_main
    real(r8), intent(in) :: sigmag_amode_in(ntot_amode_in)
    real(r8), intent(in) :: alnsg_amode_in(ntot_amode_in)
    real(r8), intent(in) :: specmw_amode_in(ntot_aspectype_in)
    real(r8), intent(in) :: specdens_amode_in(ntot_aspectype_in)
    real(r8), intent(in) :: specmw_so4_amode_in
    real(r8), intent(in) :: specdens_so4_amode_in
    real(r8), intent(in) :: specmw_nh4_amode_in
    real(r8), intent(in) :: specdens_nh4_amode_in
    real(r8), intent(in) :: specmw_soa_amode_in
    real(r8), intent(in) :: specdens_soa_amode_in
    real(r8), intent(in) :: n_so4_monolayers_pcage_in
    real(r8), intent(in) :: soa_equivso4_factor_in
    real(r8), intent(in) :: tmelt_in
    real(r8), intent(in) :: t(pcols_in,pver_in)
    real(r8), intent(in) :: pmid(pcols_in,pver_in)
    real(r8), intent(in) :: pdel(pcols_in,pver_in)
    real(r8), intent(inout) :: q(ncol,pver_in,pcnstxx_in)
    real(r8), intent(in) :: dgncur_a(pcols_in,pver_in,ntot_amode_in)
    real(r8), intent(in) :: dgncur_awet(pcols_in,pver_in,ntot_amode_in)
    real(r8), intent(in) :: wetdens_a(pcols_in,pver_in,ntot_amode_in)
    real(r8), intent(out) :: qsrflx_out(pcols_in,pcnst_in)
    logical, intent(out) :: dotend_out(pcnst_in)
    integer, intent(out) :: ierr_out

    call modal_aero_coag_kernel(                                  &
         lchnk, ncol, nstep, pcols_in, pver_in, pcnstxx_in,      &
         ntot_amode_in, pcnst_in, top_lev_in, maxspec_in,       &
         ntot_aspectype_in, maxpair_in, loffset, deltat_main,   &
         pair_option_in,                                        &
         npair_in, modefrm_in, modetoo_in, nspecfrm_in,         &
         lspecfrm_in, lspectoo_in, modeptr_accum_in,            &
         modeptr_aitken_in, modeptr_pcarbon_in, numptr_amode_in,&
         mprognum_amode_in, nspec_amode_in, lmassptr_amode_in,  &
         lspectype_amode_in, sigmag_amode_in, alnsg_amode_in,   &
         lptr_so4_a_amode_in, lptr_nh4_a_amode_in,              &
         lptr_soa_a_amode_in, specmw_amode_in,                  &
         specdens_amode_in, specmw_so4_amode_in,                &
         specdens_so4_amode_in, specmw_nh4_amode_in,            &
         specdens_nh4_amode_in, specmw_soa_amode_in,            &
         specdens_soa_amode_in, n_so4_monolayers_pcage_in,     &
         soa_equivso4_factor_in, tmelt_in, t, pmid, pdel, q,   &
         dgncur_a,                                              &
         dgncur_awet, wetdens_a, qsrflx_out, dotend_out,        &
         ierr_out)

  end subroutine modal_aero_coag_sub_run

end module ap_modal_aero_coag_scheme
