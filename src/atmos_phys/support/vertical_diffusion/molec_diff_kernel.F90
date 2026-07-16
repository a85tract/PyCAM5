module molec_diff_kernel

  use shr_kind_mod, only: r8 => shr_kind_r8
  use coords_1d, only: Coords1D
  use linear_1d_operators, only: BoundaryType, TriDiagDecomp

  implicit none
  private

  public :: compute_molec_diff_kernel
  public :: vd_lu_qdecomp_kernel

contains

  subroutine compute_molec_diff_kernel( &
       pcols, pver, ncnst, ncol, t, pmid, pint, zi, ztodt, &
       waccmx_mode, ntop_molec, nbot_molec, gravit, d0, km_fac, &
       pr_num, pwr, n_avog, mw_dry, mbarv, rairv, kmvis, kmcnd, &
       cpairv, cnst_mw_in, cnst_fixed_ubc_in, cnst_fixed_ubflx_in, &
       mw_fac_in, ubc_t_in, ubc_mmr_in, ubc_flux_in, kvm, kvt, &
       tint, rhoi, tmpi2, kq_scal, ubc_t, ubc_mmr, ubc_flux, &
       dse_top, cc_top, cnst_mw, cnst_fixed_ubc, cnst_fixed_ubflx, &
       mw_fac_out, kvt_returned)

    integer, intent(in) :: pcols, pver, ncnst, ncol
    integer, intent(in) :: ntop_molec, nbot_molec
    logical, intent(in) :: waccmx_mode
    real(r8), intent(in) :: t(pcols,pver)
    real(r8), intent(in) :: pmid(pcols,pver)
    real(r8), intent(in) :: pint(pcols,pver+1)
    real(r8), intent(in) :: zi(pcols,pver+1)
    real(r8), intent(in) :: ztodt
    real(r8), intent(in) :: gravit, d0, km_fac, pr_num, pwr
    real(r8), intent(in) :: n_avog, mw_dry
    real(r8), intent(in) :: mbarv(pcols,pver)
    real(r8), intent(in) :: rairv(pcols,pver)
    real(r8), intent(in) :: kmvis(pcols,pver+1)
    real(r8), intent(in) :: kmcnd(pcols,pver+1)
    real(r8), intent(in) :: cpairv(pcols,pver)
    real(r8), intent(in) :: cnst_mw_in(ncnst)
    logical, intent(in) :: cnst_fixed_ubc_in(ncnst)
    logical, intent(in) :: cnst_fixed_ubflx_in(ncnst)
    real(r8), intent(in) :: mw_fac_in(ncnst)
    real(r8), intent(in) :: ubc_t_in(pcols)
    real(r8), intent(in) :: ubc_mmr_in(pcols,ncnst)
    real(r8), intent(in) :: ubc_flux_in(ncnst)

    real(r8), intent(inout) :: kvm(pcols,pver+1)
    real(r8), intent(out) :: kvt(pcols,pver+1)
    real(r8), intent(inout) :: tint(pcols,pver+1)
    real(r8), intent(inout) :: rhoi(pcols,pver+1)
    real(r8), intent(inout) :: tmpi2(pcols,pver+1)
    real(r8), intent(out) :: kq_scal(pcols,pver+1)
    real(r8), intent(out) :: ubc_t(pcols)
    real(r8), intent(out) :: ubc_mmr(pcols,ncnst)
    real(r8), intent(out) :: ubc_flux(ncnst)
    real(r8), intent(out) :: dse_top(pcols)
    real(r8), intent(out) :: cc_top(pcols)
    real(r8), intent(out) :: cnst_mw(ncnst)
    logical, intent(out) :: cnst_fixed_ubc(ncnst)
    logical, intent(out) :: cnst_fixed_ubflx(ncnst)
    real(r8), intent(out) :: mw_fac_out(pcols,pver+1,ncnst)
    logical, intent(out) :: kvt_returned

    integer :: i, k, m
    real(r8) :: mbarvi
    real(r8) :: mkvisc

    kvt_returned = waccmx_mode

    ubc_t = ubc_t_in
    ubc_mmr = ubc_mmr_in
    ubc_flux = ubc_flux_in
    cnst_mw = cnst_mw_in
    cnst_fixed_ubc = cnst_fixed_ubc_in
    cnst_fixed_ubflx = cnst_fixed_ubflx_in

    do m = 1, ncnst
       if (.not. cnst_fixed_ubc(m)) then
          ubc_mmr(:,m) = 0._r8
       end if
    end do

    if (kvt_returned) then
       do m = 1, ncnst
          do k = ntop_molec+1, nbot_molec
             do i = 1, ncol
                mbarvi = 0.5_r8 * (mbarv(i,k-1)+mbarv(i,k))
                mw_fac_out(i,k,m) = d0 * mbarvi * &
                     sqrt(1._r8/mbarvi + 1._r8/cnst_mw(m)) / n_avog
             end do
          end do
          mw_fac_out(:ncol,ntop_molec,m) = &
               1.5_r8*mw_fac_out(:ncol,ntop_molec+1,m) - &
               .5_r8*mw_fac_out(:ncol,ntop_molec+2,m)
          do k = nbot_molec+1, pver+1
             mw_fac_out(:ncol,k,m) = mw_fac_out(:ncol,nbot_molec,m)
          end do
       end do
    else
       do k = 1, pver+1
          do i = 1, ncol
             mw_fac_out(i,k,:ncnst) = mw_fac_in(:ncnst)
          end do
       end do
    end if

    if (kvt_returned) then
       tint(:ncol,ntop_molec) = &
            1.5_r8*tint(:ncol,ntop_molec+1) - &
            .5_r8*tint(:ncol,ntop_molec+2)
    else
       tint(:ncol,ntop_molec) = ubc_t(:ncol)
    end if

    rhoi(:ncol,ntop_molec) = pint(:ncol,ntop_molec) / &
         (rairv(:ncol,ntop_molec) * tint(:ncol,ntop_molec))
    tmpi2(:ncol,ntop_molec) = ztodt * &
         (gravit * rhoi(:ncol,ntop_molec))**2 / &
         (pmid(:ncol,ntop_molec) - pint(:ncol,ntop_molec))

    kvt = 0._r8
    kq_scal = 0._r8
    if (kvt_returned) then
       do k = ntop_molec, nbot_molec
          do i = 1, ncol
             mkvisc = kmvis(i,k) / rhoi(i,k)
             kvm(i,k) = kvm(i,k) + mkvisc
             mkvisc = kmcnd(i,k) / rhoi(i,k)
             kvt(i,k) = mkvisc
             kq_scal(i,k) = sqrt(tint(i,k)) / rhoi(i,k)
          end do
       end do
    else
       do k = ntop_molec, nbot_molec
          do i = 1, ncol
             mkvisc = km_fac * tint(i,k)**pwr / rhoi(i,k)
             kvm(i,k) = kvm(i,k) + mkvisc
             kvt(i,k) = mkvisc * pr_num * cpairv(i,k)
             kq_scal(i,k) = sqrt(tint(i,k)) / rhoi(i,k)
          end do
       end do
    end if

    dse_top(:ncol) = cpairv(:ncol,ntop_molec) * &
         tint(:ncol,ntop_molec) + gravit * zi(:ncol,ntop_molec)

    if (kvt_returned) then
       do i = 1, ncol
          cc_top(i) = ztodt * gravit**2 * rhoi(i,ntop_molec) * &
               km_fac * ubc_t(i)**pwr / &
               (pmid(i,1) - pint(i,1))
       end do
    else
       cc_top = 0._r8
    end if

  end subroutine compute_molec_diff_kernel

  function vd_lu_qdecomp_kernel( &
       pcols, pver, ncol, fixed_ubc, mw, kv, kq_scal, mw_facm, &
       dpidz_sq, p, interface_boundary, molec_boundary, rhoi, tint, &
       ztodt, ntop_molec, nbot_molec, nbot, t, waccmx_mode, mw_dry, &
       mbarv, alphath_m, no_molec_decomp) result(decomp)

    use vdiff_lu_solver, only: fin_vol_lu_decomp

    integer, intent(in) :: pcols, pver, ncol
    integer, intent(in) :: ntop_molec, nbot_molec, nbot
    logical, intent(in) :: fixed_ubc, waccmx_mode
    real(r8), intent(in) :: mw
    real(r8), intent(in) :: kv(pcols,pver+1)
    real(r8), intent(in) :: kq_scal(pcols,pver+1)
    real(r8), intent(in) :: mw_facm(pcols,pver+1)
    real(r8), intent(in) :: dpidz_sq(ncol,pver+1)
    type(Coords1D), intent(in) :: p
    type(BoundaryType), intent(in) :: interface_boundary
    type(BoundaryType), intent(in) :: molec_boundary
    real(r8), intent(in) :: rhoi(pcols,pver+1)
    real(r8), intent(in) :: tint(pcols,pver+1)
    real(r8), intent(in) :: ztodt
    real(r8), intent(in) :: t(pcols,pver)
    real(r8), intent(in) :: mw_dry
    real(r8), intent(in) :: mbarv(pcols,pver)
    real(r8), intent(in) :: alphath_m
    type(TriDiagDecomp), intent(in) :: no_molec_decomp
    type(TriDiagDecomp) :: decomp

    integer :: k
    real(r8) :: kmq(ncol,nbot_molec+1)
    real(r8) :: mw_term(ncol,nbot_molec+1)
    real(r8) :: diff_coef(ncol,nbot_molec+1)
    real(r8) :: advect_v(ncol,nbot_molec+1)
    real(r8) :: gradm(ncol,nbot_molec+1)
    real(r8) :: gradt(ncol,nbot_molec+1)
    real(r8) :: mbarvi(ncol)

    kmq = 0._r8
    mw_term = 0._r8
    gradm = 0._r8
    gradt = 0._r8

    if (waccmx_mode) then
       k = ntop_molec
       mbarvi = .75_r8*mbarv(:ncol,k) + &
            0.5_r8*mbarv(:ncol,k+1) - .25_r8*mbarv(:ncol,k+2)
       mw_term(:,k) = (mw/mbarvi - 1._r8) / p%ifc(:,k)
       gradm(:,k) = (mbarv(:ncol,k)-mbarvi) / &
            (p%mid(:,k)-p%ifc(:,k)) / &
            (mbarv(:ncol,k)+mbarvi)*2._r8

       if (alphath_m /= 0._r8) then
          gradt(:,k) = alphath_m*(t(:ncol,k)-tint(:ncol,k)) / &
               (p%mid(:ncol,k)-p%ifc(:ncol,k)) / &
               (t(:ncol,k)+tint(:ncol,k))*2._r8
       end if

       do k = ntop_molec+1, nbot_molec
          mbarvi = 0.5_r8 * (mbarv(:ncol,k-1)+mbarv(:ncol,k))
          mw_term(:,k) = (mw/mbarvi - 1._r8) / p%ifc(:,k)
          gradm(:,k) = (mbarv(:ncol,k)-mbarv(:ncol,k-1)) * &
               p%rdst(:,k-1)/mbarvi
       end do

       if (alphath_m /= 0._r8) then
          do k = ntop_molec+1, nbot_molec
             gradt(:,k) = alphath_m*(t(:ncol,k)-t(:ncol,k-1)) * &
                  p%rdst(:,k-1)/tint(:ncol,k)
          end do
       end if
    else
       do k = ntop_molec, nbot_molec
          mw_term(:,k) = (mw/mw_dry - 1._r8) / p%ifc(:ncol,k)
       end do
    end if

    do k = ntop_molec, nbot_molec
       kmq(:,k) = kq_scal(:ncol,k) * mw_facm(:ncol,k)
    end do

    diff_coef = kv(:ncol,:nbot_molec+1) + kmq
    advect_v = kmq*mw_term
    if (waccmx_mode) then
       advect_v = advect_v - kmq*gradt - &
            (kv(:ncol,:nbot_molec+1) + kmq)*gradm
    end if

    diff_coef = dpidz_sq(:,:nbot_molec+1) * diff_coef
    advect_v = dpidz_sq(:,:nbot_molec+1) * advect_v

    if (fixed_ubc) then
       decomp = fin_vol_lu_decomp(ztodt, p, &
            coef_q_diff=diff_coef, coef_q_adv=advect_v, &
            upper_bndry=interface_boundary, &
            lower_bndry=molec_boundary, &
            graft_decomp=no_molec_decomp)
    else
       decomp = fin_vol_lu_decomp(ztodt, p, &
            coef_q_diff=diff_coef, coef_q_adv=advect_v, &
            lower_bndry=molec_boundary, &
            graft_decomp=no_molec_decomp)
    end if

  end function vd_lu_qdecomp_kernel

end module molec_diff_kernel
