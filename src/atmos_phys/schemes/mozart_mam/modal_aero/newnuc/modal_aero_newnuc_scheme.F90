module ap_modal_aero_newnuc_scheme

  use shr_kind_mod, only : r8 => shr_kind_r8
  use ap_modal_aero_newnuc_kernels, only : modal_aero_newnuc_kernel

  implicit none
  private

  public :: modal_aero_newnuc_sub_run

contains

  !> \section arg_table_modal_aero_newnuc_sub_run Argument Table
  !! \htmlinclude modal_aero_newnuc_sub_run.html
  !!
  subroutine modal_aero_newnuc_sub_run(                           &
       lchnk, ncol, nstep, pcols_in, pver_in, pcnstxx_in,        &
       pcnst_in, top_lev_in, loffset, deltat,                    &
       l_h2so4_sv_in, l_nh3_sv_in, lnumait_sv_in,                &
       lnh4ait_sv_in, lso4ait_sv_in, lptr_nh4_aitken_in,        &
       dgnumlo_aitken_in, dgnum_aitken_in, dgnumhi_aitken_in,   &
       specdens_so4_in, specmw_so4_in, specmw_nh4_in,           &
       gravity_in, t, pmid, pdel, zm, pblh, qv, qv_sat, cld, q, &
       del_h2so4_gasprod, del_h2so4_aeruptk,                    &
       qsrflx_out, dotend_out)

    integer, intent(in) :: lchnk
    integer, intent(in) :: ncol
    integer, intent(in) :: nstep
    integer, intent(in) :: pcols_in
    integer, intent(in) :: pver_in
    integer, intent(in) :: pcnstxx_in
    integer, intent(in) :: pcnst_in
    integer, intent(in) :: top_lev_in
    integer, intent(in) :: loffset
    integer, intent(in) :: l_h2so4_sv_in
    integer, intent(in) :: l_nh3_sv_in
    integer, intent(in) :: lnumait_sv_in
    integer, intent(in) :: lnh4ait_sv_in
    integer, intent(in) :: lso4ait_sv_in
    integer, intent(in) :: lptr_nh4_aitken_in

    real(r8), intent(in) :: deltat
    real(r8), intent(in) :: dgnumlo_aitken_in
    real(r8), intent(in) :: dgnum_aitken_in
    real(r8), intent(in) :: dgnumhi_aitken_in
    real(r8), intent(in) :: specdens_so4_in
    real(r8), intent(in) :: specmw_so4_in
    real(r8), intent(in) :: specmw_nh4_in
    real(r8), intent(in) :: gravity_in

    real(r8), intent(in) :: t(pcols_in,pver_in)
    real(r8), intent(in) :: pmid(pcols_in,pver_in)
    real(r8), intent(in) :: pdel(pcols_in,pver_in)
    real(r8), intent(in) :: zm(pcols_in,pver_in)
    real(r8), intent(in) :: pblh(pcols_in)
    real(r8), intent(in) :: qv(pcols_in,pver_in)
    real(r8), intent(in) :: qv_sat(pcols_in,pver_in)
    real(r8), intent(in) :: cld(ncol,pver_in)
    real(r8), intent(inout) :: q(ncol,pver_in,pcnstxx_in)
    real(r8), intent(in) :: del_h2so4_gasprod(ncol,pver_in)
    real(r8), intent(in) :: del_h2so4_aeruptk(ncol,pver_in)
    real(r8), intent(out) :: qsrflx_out(pcols_in,pcnst_in)
    logical, intent(out) :: dotend_out(pcnst_in)

    call modal_aero_newnuc_kernel(                                &
         lchnk, ncol, nstep, pcols_in, pver_in, pcnstxx_in,      &
         pcnst_in, top_lev_in, loffset, deltat,                  &
         l_h2so4_sv_in, l_nh3_sv_in, lnumait_sv_in,              &
         lnh4ait_sv_in, lso4ait_sv_in, lptr_nh4_aitken_in,      &
         dgnumlo_aitken_in, dgnum_aitken_in, dgnumhi_aitken_in, &
         specdens_so4_in, specmw_so4_in, specmw_nh4_in,         &
         gravity_in, t, pmid, pdel, zm, pblh, qv, qv_sat, cld,  &
         q, del_h2so4_gasprod, del_h2so4_aeruptk,               &
         qsrflx_out, dotend_out)

  end subroutine modal_aero_newnuc_sub_run

end module ap_modal_aero_newnuc_scheme
