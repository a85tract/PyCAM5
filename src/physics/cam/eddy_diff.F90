module eddy_diff

  use shr_kind_mod, only : r8 => shr_kind_r8, i4 => shr_kind_i4
  use ap_compute_eddy_diff_scheme, only : init_eddy_diff, &
       compute_eddy_diff_run
  use perf_mod, only : t_startf, t_stopf

  implicit none
  private

  public :: init_eddy_diff
  public :: compute_eddy_diff

contains

  subroutine compute_eddy_diff( lchnk, pcols, pver, ncol, t, qv, ztodt, &
       ql, qi, s, pdel, rpdel, cldn, qrl, wsedl, z, zi, pmid, pi, u, v, &
       taux, tauy, shflx, qflx, wstarent, nturb, rrho, ustar, pblh, &
       kvm_in, kvh_in, kvm_out, kvh_out, kvq, cgh, cgs, tpert, qpert, &
       wpert, tke, bprod, sprod, sfi, kvinit, tauresx, tauresy, ksrftms, &
       ipbl, kpblh, wstarPBL, tkes, went, turbtype, sm_aw )

    integer, intent(in) :: lchnk
    integer, intent(in) :: pcols
    integer, intent(in) :: pver
    integer, intent(in) :: ncol
    integer, intent(in) :: nturb
    logical, intent(in) :: wstarent
    logical, intent(in) :: kvinit
    real(r8), intent(in) :: ztodt
    real(r8), intent(in) :: t(pcols,pver)
    real(r8), intent(in) :: qv(pcols,pver)
    real(r8), intent(in) :: ql(pcols,pver)
    real(r8), intent(in) :: qi(pcols,pver)
    real(r8), intent(in) :: s(pcols,pver)
    real(r8), intent(in) :: pdel(pcols,pver)
    real(r8), intent(in) :: rpdel(pcols,pver)
    real(r8), intent(in) :: cldn(pcols,pver)
    real(r8), intent(in) :: qrl(pcols,pver)
    real(r8), intent(in) :: wsedl(pcols,pver)
    real(r8), intent(in) :: z(pcols,pver)
    real(r8), intent(in) :: zi(pcols,pver+1)
    real(r8), intent(in) :: pmid(pcols,pver)
    real(r8), intent(in) :: pi(pcols,pver+1)
    real(r8), intent(in) :: u(pcols,pver)
    real(r8), intent(in) :: v(pcols,pver)
    real(r8), intent(in) :: taux(pcols)
    real(r8), intent(in) :: tauy(pcols)
    real(r8), intent(in) :: shflx(pcols)
    real(r8), intent(in) :: qflx(pcols)
    real(r8), intent(in) :: kvm_in(pcols,pver+1)
    real(r8), intent(in) :: kvh_in(pcols,pver+1)
    real(r8), intent(in) :: ksrftms(pcols)

    real(r8), intent(out) :: kvm_out(pcols,pver+1)
    real(r8), intent(out) :: kvh_out(pcols,pver+1)
    real(r8), intent(out) :: kvq(pcols,pver+1)
    real(r8), intent(out) :: rrho(pcols)
    real(r8), intent(out) :: ustar(pcols)
    real(r8), intent(out) :: pblh(pcols)
    real(r8), intent(out) :: cgh(pcols,pver+1)
    real(r8), intent(out) :: cgs(pcols,pver+1)
    real(r8), intent(out) :: tpert(pcols)
    real(r8), intent(out) :: qpert(pcols)
    real(r8), intent(out) :: wpert(pcols)
    real(r8), intent(out) :: tke(pcols,pver+1)
    real(r8), intent(out) :: bprod(pcols,pver+1)
    real(r8), intent(out) :: sprod(pcols,pver+1)
    real(r8), intent(out) :: sfi(pcols,pver+1)
    integer(i4), intent(out) :: turbtype(pcols,pver+1)
    real(r8), intent(out) :: sm_aw(pcols,pver+1)
    integer(i4), intent(out) :: ipbl(pcols)
    integer(i4), intent(out) :: kpblh(pcols)
    real(r8), intent(out) :: wstarPBL(pcols)
    real(r8), intent(out) :: tkes(pcols)
    real(r8), intent(out) :: went(pcols)

    real(r8), intent(inout) :: tauresx(pcols)
    real(r8), intent(inout) :: tauresy(pcols)

    character(len=512) :: errmsg
    integer :: errflg

    call t_startf('ap_compute_eddy_diff_run')
    call compute_eddy_diff_run(lchnk, pcols, pver, pver+1, ncol, t, qv, &
         ztodt, ql, qi, s, pdel, rpdel, cldn, qrl, wsedl, z, zi, pmid, &
         pi, u, v, taux, tauy, shflx, qflx, wstarent, nturb, rrho, &
         ustar, pblh, kvm_in, kvh_in, kvm_out, kvh_out, kvq, cgh, cgs, &
         tpert, qpert, wpert, tke, bprod, sprod, sfi, kvinit, tauresx, &
         tauresy, ksrftms, ipbl, kpblh, wstarPBL, tkes, went, turbtype, &
         sm_aw, errmsg, errflg)
    call t_stopf('ap_compute_eddy_diff_run')

  end subroutine compute_eddy_diff

end module eddy_diff
