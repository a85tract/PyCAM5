from math import exp, sqrt


@export
def gw_tend_select_branches_codon(
    do_molec_diff: int,
    use_gw_convect_dp: int,
    use_gw_convect_sh: int,
    use_gw_front: int,
    use_gw_front_igw: int,
    use_gw_oro: int,
    branch_mask_p: cobj,
):
    branch_mask = Ptr[int](branch_mask_p)

    mask = 0
    if do_molec_diff != 0:
        mask |= 1
    if use_gw_convect_dp != 0:
        mask |= 2
    if use_gw_convect_sh != 0:
        mask |= 4
    if use_gw_front != 0:
        mask |= 8
    if use_gw_front_igw != 0:
        mask |= 16
    if use_gw_oro != 0:
        mask |= 32

    branch_mask[0] = mask


@inline
def _idx2(i: int, k: int, ld1: int) -> int:
    """Fortran array(i,k) with first dimension ld1."""
    return (i - 1) + (k - 1) * ld1


@inline
def _idx3(i: int, k: int, m: int, ld1: int, ld2: int) -> int:
    """Fortran array(i,k,m) with dimensions (ld1,ld2,*)."""
    return (i - 1) + (k - 1) * ld1 + (m - 1) * ld1 * ld2


@inline
def _idx_tau(i: int, l: int, k: int, ncol: int, ngwv: int) -> int:
    """Fortran tau(i,l,k) with dimensions (ncol,-ngwv:ngwv,*)."""
    return (i - 1) + (l + ngwv) * ncol + (k - 1) * ncol * (2 * ngwv + 1)


@inline
def _idx_c(i: int, l: int, ncol: int, ngwv: int) -> int:
    """Fortran c(i,l) with dimensions (ncol,-ngwv:ngwv)."""
    return (i - 1) + (l + ngwv) * ncol


@inline
def _idx_gwut(i: int, k: int, l: int, ncol: int, pver: int, ngwv: int) -> int:
    """Fortran gwut(i,k,l) with dimensions (ncol,pver,-ngwv:ngwv)."""
    return (i - 1) + (k - 1) * ncol + (l + ngwv) * ncol * pver


@inline
def _sign_fortran(a: float, b: float) -> float:
    if b >= 0.0:
        return abs(a)
    return -abs(a)


@export
def gw_common_init_scalars_codon(
    pver_in: int,
    ktop_in: int,
    gravit_in: float,
    rair_in: float,
    tau_0_ubc_in: int,
    pver_p: cobj,
    tau_0_ubc_p: cobj,
    ktop_p: cobj,
    gravit_p: cobj,
    rair_p: cobj,
    rog_p: cobj,
):
    pver = Ptr[i32](pver_p)
    tau_0_ubc = Ptr[i32](tau_0_ubc_p)
    ktop = Ptr[i32](ktop_p)
    gravit = Ptr[float](gravit_p)
    rair = Ptr[float](rair_p)
    rog = Ptr[float](rog_p)

    pver[0] = i32(pver_in)
    tau_0_ubc[0] = i32(tau_0_ubc_in)
    ktop[0] = i32(ktop_in)
    gravit[0] = gravit_in
    rair[0] = rair_in
    rog[0] = rair[0] / gravit[0]


@export
def gw_common_new_gwband_codon(
    ngwv_in: int,
    dc_in: float,
    fcrit2_in: float,
    wavelength: float,
    pi_in: float,
    ngwv_p: cobj,
    dc_p: cobj,
    fcrit2_p: cobj,
    cref_p: cobj,
    kwv_p: cobj,
    effkwv_p: cobj,
):
    ngwv = Ptr[i32](ngwv_p)
    dc = Ptr[float](dc_p)
    fcrit2 = Ptr[float](fcrit2_p)
    cref = Ptr[float](cref_p)
    kwv = Ptr[float](kwv_p)
    effkwv = Ptr[float](effkwv_p)

    ngwv[0] = i32(ngwv_in)
    dc[0] = dc_in
    fcrit2[0] = fcrit2_in

    for l in range(-ngwv_in, ngwv_in + 1):
        cref[l + ngwv_in] = dc_in * float(l)

    kwv[0] = 2.0 * pi_in / wavelength
    effkwv[0] = fcrit2[0] * kwv[0]


@export
def gw_common_init_alpha_codon(n: int, alpha_in_p: cobj, alpha_p: cobj):
    alpha_in = Ptr[float](alpha_in_p)
    alpha = Ptr[float](alpha_p)

    for i in range(n):
        alpha[i] = alpha_in[i]


@export
def gw_prof_codon(
    ncol: int,
    pver: int,
    cpair: float,
    rair: float,
    gravit: float,
    p_ifc_p: cobj,
    p_rdst_p: cobj,
    t_p: cobj,
    rhoi_p: cobj,
    nm_p: cobj,
    ni_p: cobj,
    ti_p: cobj,
):
    p_ifc = Ptr[float](p_ifc_p)
    p_rdst = Ptr[float](p_rdst_p)
    t = Ptr[float](t_p)
    rhoi = Ptr[float](rhoi_p)
    nm = Ptr[float](nm_p)
    ni = Ptr[float](ni_p)
    ti = Ptr[float](ti_p)

    n2min = 5.0e-5

    k = 1
    for i in range(1, ncol + 1):
        ti[_idx2(i, k, ncol)] = t[_idx2(i, k, ncol)]
        rhoi[_idx2(i, k, ncol)] = p_ifc[_idx2(i, k, ncol)] / (rair * ti[_idx2(i, k, ncol)])
        ni[_idx2(i, k, ncol)] = sqrt((gravit * gravit) / (cpair * ti[_idx2(i, k, ncol)]))

    for k_mid in range(1, pver):
        for i in range(1, ncol + 1):
            ti[_idx2(i, k_mid + 1, ncol)] = 0.5 * (t[_idx2(i, k_mid, ncol)] + t[_idx2(i, k_mid + 1, ncol)])

    for k in range(2, pver + 1):
        for i in range(1, ncol + 1):
            rhoi[_idx2(i, k, ncol)] = p_ifc[_idx2(i, k, ncol)] / (rair * ti[_idx2(i, k, ncol)])
            dtdp = (t[_idx2(i, k, ncol)] - t[_idx2(i, k - 1, ncol)]) * p_rdst[_idx2(i, k - 1, ncol)]
            n2 = ((gravit * gravit) / ti[_idx2(i, k, ncol)]) * ((1.0 / cpair) - (rhoi[_idx2(i, k, ncol)] * dtdp))
            ni[_idx2(i, k, ncol)] = sqrt(max(n2min, n2))

    k = pver + 1
    for i in range(1, ncol + 1):
        ti[_idx2(i, k, ncol)] = t[_idx2(i, k - 1, ncol)]
        rhoi[_idx2(i, k, ncol)] = p_ifc[_idx2(i, k, ncol)] / (rair * ti[_idx2(i, k, ncol)])
        ni[_idx2(i, k, ncol)] = ni[_idx2(i, k - 1, ncol)]

    for k_mid in range(1, pver + 1):
        for i in range(1, ncol + 1):
            nm[_idx2(i, k_mid, ncol)] = 0.5 * (ni[_idx2(i, k_mid, ncol)] + ni[_idx2(i, k_mid + 1, ncol)])


@export
def gw_prof_stage_dispatch_codon(
    ncol: int,
    pver: int,
    cpair: float,
    rair: float,
    gravit: float,
    p_ifc_p: cobj,
    p_rdst_p: cobj,
    t_p: cobj,
    rhoi_p: cobj,
    nm_p: cobj,
    ni_p: cobj,
    ti_p: cobj,
):
    gw_prof_codon(
        ncol,
        pver,
        cpair,
        rair,
        gravit,
        p_ifc_p,
        p_rdst_p,
        t_p,
        rhoi_p,
        nm_p,
        ni_p,
        ti_p,
    )


@export
def gw_energy_change_codon(
    ncol: int,
    pver: int,
    dt: float,
    gravit: float,
    p_del_p: cobj,
    u_p: cobj,
    v_p: cobj,
    dudt_p: cobj,
    dvdt_p: cobj,
    dsdt_p: cobj,
    de_p: cobj,
):
    p_del = Ptr[float](p_del_p)
    u = Ptr[float](u_p)
    v = Ptr[float](v_p)
    dudt = Ptr[float](dudt_p)
    dvdt = Ptr[float](dvdt_p)
    dsdt = Ptr[float](dsdt_p)
    de = Ptr[float](de_p)

    for i in range(1, ncol + 1):
        de[i - 1] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, ncol)
            de[i - 1] = de[i - 1] + p_del[idx] / gravit * (
                dsdt[idx]
                + dudt[idx] * (u[idx] + dudt[idx] * 0.5 * dt)
                + dvdt[idx] * (v[idx] + dvdt[idx] * 0.5 * dt)
            )


@export
def gw_energy_change_stage_dispatch_codon(
    ncol: int,
    pver: int,
    dt: float,
    gravit: float,
    p_del_p: cobj,
    u_p: cobj,
    v_p: cobj,
    dudt_p: cobj,
    dvdt_p: cobj,
    dsdt_p: cobj,
    de_p: cobj,
):
    gw_energy_change_codon(
        ncol,
        pver,
        dt,
        gravit,
        p_del_p,
        u_p,
        v_p,
        dudt_p,
        dvdt_p,
        dsdt_p,
        de_p,
    )


def _gw_drag_prof_core_impl(
    stage: int,
    ncol: int,
    pver: int,
    pverp: int,
    ngwv: int,
    ktop: int,
    kbot_tend: int,
    kbot_src: int,
    tau_0_ubc: int,
    dback: float,
    taumin: float,
    tndmax: float,
    umcfac: float,
    ubmc2mn: float,
    effkwv: float,
    kwv: float,
    gravit: float,
    rog: float,
    dt: float,
    alpha_p: cobj,
    p_del_p: cobj,
    p_rdel_p: cobj,
    t_p: cobj,
    piln_p: cobj,
    rhoi_p: cobj,
    ni_p: cobj,
    ubm_p: cobj,
    ubi_p: cobj,
    xv_p: cobj,
    yv_p: cobj,
    effgw_p: cobj,
    c_p: cobj,
    kvtt_p: cobj,
    src_level_p: cobj,
    tend_level_p: cobj,
    tau_p: cobj,
    utgw_p: cobj,
    vtgw_p: cobj,
    ttgw_p: cobj,
    gwut_p: cobj,
    dttdf_p: cobj,
    dttke_p: cobj,
    d_p: cobj,
    mi_p: cobj,
    taudmp_p: cobj,
    tausat_p: cobj,
    ubmc_p: cobj,
    ubmc2_p: cobj,
    ubt_p: cobj,
    ubtl_p: cobj,
    wrk_p: cobj,
    ubt_lim_ratio_p: cobj,
):
    alpha = Ptr[float](alpha_p)
    p_del = Ptr[float](p_del_p)
    p_rdel = Ptr[float](p_rdel_p)
    t = Ptr[float](t_p)
    piln = Ptr[float](piln_p)
    rhoi = Ptr[float](rhoi_p)
    ni = Ptr[float](ni_p)
    ubm = Ptr[float](ubm_p)
    ubi = Ptr[float](ubi_p)
    xv = Ptr[float](xv_p)
    yv = Ptr[float](yv_p)
    effgw = Ptr[float](effgw_p)
    c = Ptr[float](c_p)
    kvtt = Ptr[float](kvtt_p)
    src_level = Ptr[int](src_level_p)
    tend_level = Ptr[int](tend_level_p)
    tau = Ptr[float](tau_p)
    utgw = Ptr[float](utgw_p)
    vtgw = Ptr[float](vtgw_p)
    ttgw = Ptr[float](ttgw_p)
    gwut = Ptr[float](gwut_p)
    dttdf = Ptr[float](dttdf_p)
    dttke = Ptr[float](dttke_p)
    d = Ptr[float](d_p)
    mi = Ptr[float](mi_p)
    taudmp = Ptr[float](taudmp_p)
    tausat = Ptr[float](tausat_p)
    ubmc = Ptr[float](ubmc_p)
    ubmc2 = Ptr[float](ubmc2_p)
    ubt = Ptr[float](ubt_p)
    ubtl = Ptr[float](ubtl_p)
    wrk = Ptr[float](wrk_p)
    ubt_lim_ratio = Ptr[float](ubt_lim_ratio_p)

    if stage == 1:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, ncol)
                utgw[idx] = 0.0
                vtgw[idx] = 0.0
                dttke[idx] = 0.0
                ttgw[idx] = 0.0

        for l in range(-ngwv, ngwv + 1):
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    gwut[_idx_gwut(i, k, l, ncol, pver, ngwv)] = 0.0

        for i in range(1, ncol + 1):
            mi[i - 1] = 0.0
            taudmp[i - 1] = 0.0
            tausat[i - 1] = 0.0
            ubmc[i - 1] = 0.0
            ubmc2[i - 1] = 0.0
            wrk[i - 1] = 0.0

        for k in range(kbot_src, ktop - 1, -1):
            for i in range(1, ncol + 1):
                d[i - 1] = dback + kvtt[_idx2(i, k, ncol)]

            for l in range(-ngwv, ngwv + 1):
                for i in range(1, ncol + 1):
                    ubmc[i - 1] = ubi[_idx2(i, k, ncol)] - c[_idx_c(i, l, ncol, ngwv)]

                for i in range(1, ncol + 1):
                    tausat[i - 1] = 0.0

                for i in range(1, ncol + 1):
                    c_val = c[_idx_c(i, l, ncol, ngwv)]
                    if src_level[i - 1] >= k:
                        if (ubmc[i - 1] > 0.0) == (ubi[_idx2(i, k + 1, ncol)] > c_val):
                            tausat[i - 1] = abs(
                                effkwv * rhoi[_idx2(i, k, ncol)] * ubmc[i - 1] ** 3
                                / (2.0 * ni[_idx2(i, k, ncol)])
                            )

                for i in range(1, ncol + 1):
                    if src_level[i - 1] >= k:
                        idx = _idx2(i, k, ncol)
                        c_idx = _idx_c(i, l, ncol, ngwv)
                        tau_idx_next = _idx_tau(i, l, k + 1, ncol, ngwv)
                        tau_idx = _idx_tau(i, l, k, ncol, ngwv)

                        ubmc2[i - 1] = max(ubmc[i - 1] ** 2, ubmc2mn)
                        mi[i - 1] = ni[idx] / (2.0 * kwv * ubmc2[i - 1]) * (
                            alpha[k - 1] + ni[idx] ** 2 / ubmc2[i - 1] * d[i - 1]
                        )
                        wrk[i - 1] = -2.0 * mi[i - 1] * rog * t[idx] * (
                            piln[_idx2(i, k + 1, ncol)] - piln[idx]
                        )

                        taudmp[i - 1] = tau[tau_idx_next] * exp(wrk[i - 1])

                        if tausat[i - 1] <= taumin:
                            tausat[i - 1] = 0.0
                        if taudmp[i - 1] <= taumin:
                            taudmp[i - 1] = 0.0

                        tau[tau_idx] = min(taudmp[i - 1], tausat[i - 1])

        if tau_0_ubc != 0:
            for l in range(-ngwv, ngwv + 1):
                for i in range(1, ncol + 1):
                    tau[_idx_tau(i, l, ktop, ncol, ngwv)] = 0.0

        for k in range(ktop, kbot_tend + 2):
            for l in range(-ngwv, ngwv + 1):
                for i in range(1, ncol + 1):
                    if k - 1 <= tend_level[i - 1]:
                        tau_idx = _idx_tau(i, l, k, ncol, ngwv)
                        tau[tau_idx] = tau[tau_idx] * effgw[i - 1]

        for k in range(ktop, kbot_tend + 1):
            for i in range(1, ncol + 1):
                ubt[_idx2(i, k, ncol)] = 0.0

            for l in range(-ngwv, ngwv + 1):
                for i in range(1, ncol + 1):
                    ubtl[i - 1] = (
                        gravit
                        * (tau[_idx_tau(i, l, k + 1, ncol, ngwv)] - tau[_idx_tau(i, l, k, ncol, ngwv)])
                        * p_rdel[_idx2(i, k, ncol)]
                    )

                for i in range(1, ncol + 1):
                    c_val = c[_idx_c(i, l, ncol, ngwv)]
                    ubtl[i - 1] = min(ubtl[i - 1], umcfac * abs(c_val - ubm[_idx2(i, k, ncol)]) / dt)

                for i in range(1, ncol + 1):
                    if k <= tend_level[i - 1]:
                        idx = _idx2(i, k, ncol)
                        gwut_idx = _idx_gwut(i, k, l, ncol, pver, ngwv)
                        gwut[gwut_idx] = _sign_fortran(ubtl[i - 1], c[_idx_c(i, l, ncol, ngwv)] - ubm[idx])
                        ubt[idx] = ubt[idx] + gwut[gwut_idx]

            for i in range(1, ncol + 1):
                idx = _idx2(i, k, ncol)
                if abs(ubt[idx]) > tndmax:
                    ubt_lim_ratio[i - 1] = tndmax / abs(ubt[idx])
                    ubt[idx] = ubt_lim_ratio[i - 1] * ubt[idx]
                else:
                    ubt_lim_ratio[i - 1] = 1.0

            for l in range(-ngwv, ngwv + 1):
                for i in range(1, ncol + 1):
                    gwut_idx = _idx_gwut(i, k, l, ncol, pver, ngwv)
                    gwut[gwut_idx] = ubt_lim_ratio[i - 1] * gwut[gwut_idx]

                for i in range(1, ncol + 1):
                    if k <= tend_level[i - 1]:
                        tau[_idx_tau(i, l, k + 1, ncol, ngwv)] = tau[_idx_tau(i, l, k, ncol, ngwv)] + (
                            abs(gwut[_idx_gwut(i, k, l, ncol, pver, ngwv)]) * p_del[_idx2(i, k, ncol)] / gravit
                        )

            for i in range(1, ncol + 1):
                if k <= tend_level[i - 1]:
                    idx = _idx2(i, k, ncol)
                    utgw[idx] = ubt[idx] * xv[i - 1]
                    vtgw[idx] = ubt[idx] * yv[i - 1]

    elif stage == 2:
        for l in range(-ngwv, ngwv + 1):
            for k in range(ktop, kbot_tend + 1):
                for i in range(1, ncol + 1):
                    idx = _idx2(i, k, ncol)
                    dttke[idx] = dttke[idx] - (
                        (ubm[idx] - c[_idx_c(i, l, ncol, ngwv)])
                        * gwut[_idx_gwut(i, k, l, ncol, pver, ngwv)]
                    )

        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, ncol)
                ttgw[idx] = dttke[idx] + dttdf[idx]


@export
def gw_ediff_prep_codon(
    ncol: int,
    pver: int,
    pverp: int,
    ngwv: int,
    kbot: int,
    ktop: int,
    prndl: float,
    gravit: float,
    gwut_p: cobj,
    ubm_p: cobj,
    nm_p: cobj,
    rho_p: cobj,
    c_p: cobj,
    tend_level_p: cobj,
    egwdffi_p: cobj,
    egwdffm_p: cobj,
    egwdff_lev_p: cobj,
    dpidz_sq_p: cobj,
):
    gwut = Ptr[float](gwut_p)
    ubm = Ptr[float](ubm_p)
    nm = Ptr[float](nm_p)
    rho = Ptr[float](rho_p)
    c = Ptr[float](c_p)
    tend_level = Ptr[int](tend_level_p)
    egwdffi = Ptr[float](egwdffi_p)
    egwdffm = Ptr[float](egwdffm_p)
    egwdff_lev = Ptr[float](egwdff_lev_p)
    dpidz_sq = Ptr[float](dpidz_sq_p)

    for k in range(1, pverp + 1):
        for i in range(1, ncol + 1):
            egwdffi[_idx2(i, k, ncol)] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            egwdffm[_idx2(i, k, ncol)] = 0.0

    for l in range(-ngwv, ngwv + 1):
        for k in range(ktop, kbot + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, ncol)
                egwdff_lev[i - 1] = (
                    prndl
                    * 0.5
                    * gwut[_idx_gwut(i, k, l, ncol, pver, ngwv)]
                    * (c[_idx_c(i, l, ncol, ngwv)] - ubm[idx])
                    / (nm[idx] ** 2)
                )

            for i in range(1, ncol + 1):
                idx = _idx2(i, k, ncol)
                egwdffm[idx] = egwdffm[idx] + egwdff_lev[i - 1]

    for k in range(ktop + 1, kbot + 1):
        for i in range(1, ncol + 1):
            egwdffi[_idx2(i, k, ncol)] = 0.5 * (egwdffm[_idx2(i, k - 1, ncol)] + egwdffm[_idx2(i, k, ncol)])

    for k in range(ktop + 1, kbot + 1):
        for i in range(1, ncol + 1):
            if k > tend_level[i - 1]:
                egwdffi[_idx2(i, k, ncol)] = 0.0

    for k in range(1, pverp + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, ncol)
            dpidz_sq[idx] = rho[idx] * gravit

    for k in range(1, pverp + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, ncol)
            dpidz_sq[idx] = dpidz_sq[idx] * dpidz_sq[idx]


@export
def gw_diff_tend_prepost_codon(
    stage: int,
    ncol: int,
    pver: int,
    dt: float,
    q_p: cobj,
    qnew_p: cobj,
    dq_p: cobj,
):
    q = Ptr[float](q_p)
    qnew = Ptr[float](qnew_p)
    dq = Ptr[float](dq_p)

    if stage == 1:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, ncol)
                dq[idx] = 0.0
                qnew[idx] = q[idx]

    elif stage == 2:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, ncol)
                dq[idx] = (qnew[idx] - q[idx]) / dt


def _gw_diff_solver_impl(
    stage: int,
    ncol: int,
    pver: int,
    pverp: int,
    pcnst: int,
    ngwv: int,
    kbot: int,
    ktop: int,
    dt: float,
    gravit: float,
    gwut_p: cobj,
    ubm_p: cobj,
    nm_p: cobj,
    rho_p: cobj,
    c_p: cobj,
    tend_level_p: cobj,
    p_del_p: cobj,
    p_rdel_p: cobj,
    p_rdst_p: cobj,
    q_p: cobj,
    dse_p: cobj,
    egwdffi_p: cobj,
    qtgw_p: cobj,
    dttdf_p: cobj,
    egwdffm_p: cobj,
    egwdff_lev_p: cobj,
    dpidz_sq_p: cobj,
    coef_q_diff_p: cobj,
    qnew_p: cobj,
    spr_p: cobj,
    sub_p: cobj,
    diag_p: cobj,
    ca_p: cobj,
    ze_p: cobj,
    dnom_p: cobj,
    zf_p: cobj,
):
    gwut = Ptr[float](gwut_p)
    ubm = Ptr[float](ubm_p)
    nm = Ptr[float](nm_p)
    rho = Ptr[float](rho_p)
    c = Ptr[float](c_p)
    tend_level = Ptr[int](tend_level_p)
    p_del = Ptr[float](p_del_p)
    p_rdel = Ptr[float](p_rdel_p)
    p_rdst = Ptr[float](p_rdst_p)
    q = Ptr[float](q_p)
    dse = Ptr[float](dse_p)
    egwdffi = Ptr[float](egwdffi_p)
    qtgw = Ptr[float](qtgw_p)
    dttdf = Ptr[float](dttdf_p)
    egwdffm = Ptr[float](egwdffm_p)
    egwdff_lev = Ptr[float](egwdff_lev_p)
    dpidz_sq = Ptr[float](dpidz_sq_p)
    coef_q_diff = Ptr[float](coef_q_diff_p)
    qnew = Ptr[float](qnew_p)
    spr = Ptr[float](spr_p)
    sub = Ptr[float](sub_p)
    diag = Ptr[float](diag_p)
    ca = Ptr[float](ca_p)
    ze = Ptr[float](ze_p)
    dnom = Ptr[float](dnom_p)
    zf = Ptr[float](zf_p)

    prndl = 0.25
    ncel = kbot - ktop + 1

    if stage == 1:
        for k in range(1, pverp + 1):
            for i in range(1, ncol + 1):
                egwdffi[_idx2(i, k, ncol)] = 0.0

        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                egwdffm[_idx2(i, k, ncol)] = 0.0

        for l in range(-ngwv, ngwv + 1):
            for k in range(ktop, kbot + 1):
                for i in range(1, ncol + 1):
                    idx = _idx2(i, k, ncol)
                    egwdff_lev[i - 1] = (
                        prndl
                        * 0.5
                        * gwut[_idx_gwut(i, k, l, ncol, pver, ngwv)]
                        * (c[_idx_c(i, l, ncol, ngwv)] - ubm[idx])
                        / (nm[idx] ** 2)
                    )

                for i in range(1, ncol + 1):
                    idx = _idx2(i, k, ncol)
                    egwdffm[idx] = egwdffm[idx] + egwdff_lev[i - 1]

        for k in range(ktop + 1, kbot + 1):
            for i in range(1, ncol + 1):
                egwdffi[_idx2(i, k, ncol)] = 0.5 * (
                    egwdffm[_idx2(i, k - 1, ncol)] + egwdffm[_idx2(i, k, ncol)]
                )

        for k in range(ktop + 1, kbot + 1):
            for i in range(1, ncol + 1):
                if k > tend_level[i - 1]:
                    egwdffi[_idx2(i, k, ncol)] = 0.0

        for k in range(1, pverp + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, ncol)
                dpidz_sq[idx] = rho[idx] * gravit

        for k in range(1, pverp + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, ncol)
                dpidz_sq[idx] = dpidz_sq[idx] * dpidz_sq[idx]

        for sk in range(1, ncel + 2):
            k = ktop + sk - 1
            for i in range(1, ncol + 1):
                coef_q_diff[_idx2(i, sk, ncol)] = (
                    egwdffi[_idx2(i, k, ncol)] * dpidz_sq[_idx2(i, k, ncol)]
                )

    elif stage == 2:
        for sk in range(1, ncel):
            k = ktop + sk - 1
            for i in range(1, ncol + 1):
                spr[_idx2(i, sk, ncol)] = (
                    coef_q_diff[_idx2(i, sk + 1, ncol)]
                    * p_rdst[_idx2(i, k, ncol)]
                    * p_rdel[_idx2(i, k, ncol)]
                )
                sub[_idx2(i, sk, ncol)] = (
                    coef_q_diff[_idx2(i, sk + 1, ncol)]
                    * p_rdst[_idx2(i, k, ncol)]
                    * p_rdel[_idx2(i, k + 1, ncol)]
                )

        for sk in range(1, ncel):
            for i in range(1, ncol + 1):
                diag[_idx2(i, sk, ncol)] = -spr[_idx2(i, sk, ncol)]

        for i in range(1, ncol + 1):
            diag[_idx2(i, ncel, ncol)] = -0.0

        for sk in range(2, ncel + 1):
            for i in range(1, ncol + 1):
                diag[_idx2(i, sk, ncol)] = diag[_idx2(i, sk, ncol)] - sub[_idx2(i, sk - 1, ncol)]

        for sk in range(1, ncel):
            for i in range(1, ncol + 1):
                spr[_idx2(i, sk, ncol)] = spr[_idx2(i, sk, ncol)] * (-dt)

        for sk in range(1, ncel):
            for i in range(1, ncol + 1):
                sub[_idx2(i, sk, ncol)] = sub[_idx2(i, sk, ncol)] * (-dt)

        for sk in range(1, ncel + 1):
            for i in range(1, ncol + 1):
                diag[_idx2(i, sk, ncol)] = diag[_idx2(i, sk, ncol)] * (-dt)

        for sk in range(1, ncel + 1):
            for i in range(1, ncol + 1):
                diag[_idx2(i, sk, ncol)] = diag[_idx2(i, sk, ncol)] + 1.0

        for sk in range(1, ncel):
            for i in range(1, ncol + 1):
                ca[_idx2(i, sk, ncol)] = -spr[_idx2(i, sk, ncol)]

        for i in range(1, ncol + 1):
            ca[_idx2(i, ncel, ncol)] = -0.0

        for i in range(1, ncol + 1):
            dnom[_idx2(i, ncel, ncol)] = 1.0 / diag[_idx2(i, ncel, ncol)]

        for sk in range(ncel - 1, 0, -1):
            for i in range(1, ncol + 1):
                ze[_idx2(i, sk + 1, ncol)] = -sub[_idx2(i, sk, ncol)] * dnom[_idx2(i, sk + 1, ncol)]

            for i in range(1, ncol + 1):
                dnom[_idx2(i, sk, ncol)] = 1.0 / (
                    diag[_idx2(i, sk, ncol)] - ca[_idx2(i, sk, ncol)] * ze[_idx2(i, sk + 1, ncol)]
                )

        for i in range(1, ncol + 1):
            ze[_idx2(i, 1, ncol)] = -0.0

        for m in range(1, pcnst + 1):
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    idx2 = _idx2(i, k, ncol)
                    idx3 = _idx3(i, k, m, ncol, pver)
                    qtgw[idx3] = 0.0
                    qnew[idx2] = q[idx3]

            for i in range(1, ncol + 1):
                zf[_idx2(i, ncel, ncol)] = qnew[_idx2(i, kbot, ncol)] * dnom[_idx2(i, ncel, ncol)]

            for sk in range(ncel - 1, 0, -1):
                k = ktop + sk - 1
                for i in range(1, ncol + 1):
                    zf[_idx2(i, sk, ncol)] = (
                        qnew[_idx2(i, k, ncol)] + ca[_idx2(i, sk, ncol)] * zf[_idx2(i, sk + 1, ncol)]
                    ) * dnom[_idx2(i, sk, ncol)]

            for i in range(1, ncol + 1):
                qnew[_idx2(i, ktop, ncol)] = zf[_idx2(i, 1, ncol)]

            for sk in range(2, ncel + 1):
                k = ktop + sk - 1
                for i in range(1, ncol + 1):
                    qnew[_idx2(i, k, ncol)] = (
                        zf[_idx2(i, sk, ncol)] + ze[_idx2(i, sk, ncol)] * qnew[_idx2(i, k - 1, ncol)]
                    )

            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    idx2 = _idx2(i, k, ncol)
                    idx3 = _idx3(i, k, m, ncol, pver)
                    qtgw[idx3] = (qnew[idx2] - q[idx3]) / dt

        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, ncol)
                dttdf[idx] = 0.0
                qnew[idx] = dse[idx]

        for i in range(1, ncol + 1):
            zf[_idx2(i, ncel, ncol)] = qnew[_idx2(i, kbot, ncol)] * dnom[_idx2(i, ncel, ncol)]

        for sk in range(ncel - 1, 0, -1):
            k = ktop + sk - 1
            for i in range(1, ncol + 1):
                zf[_idx2(i, sk, ncol)] = (
                    qnew[_idx2(i, k, ncol)] + ca[_idx2(i, sk, ncol)] * zf[_idx2(i, sk + 1, ncol)]
                ) * dnom[_idx2(i, sk, ncol)]

        for i in range(1, ncol + 1):
            qnew[_idx2(i, ktop, ncol)] = zf[_idx2(i, 1, ncol)]

        for sk in range(2, ncel + 1):
            k = ktop + sk - 1
            for i in range(1, ncol + 1):
                qnew[_idx2(i, k, ncol)] = (
                    zf[_idx2(i, sk, ncol)] + ze[_idx2(i, sk, ncol)] * qnew[_idx2(i, k - 1, ncol)]
                )

        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, ncol)
                dttdf[idx] = (qnew[idx] - dse[idx]) / dt


def _gw_common_driver_dispatch(
    group: int,
    stage: int,
    ncol: int,
    pver: int,
    pverp: int,
    pcnst: int,
    ngwv: int,
    ktop: int,
    kbot_tend: int,
    kbot_src: int,
    kbot: int,
    tau_0_ubc: int,
    dback: float,
    taumin: float,
    tndmax: float,
    umcfac: float,
    ubmc2mn: float,
    effkwv: float,
    kwv: float,
    gravit: float,
    rog: float,
    dt: float,
    alpha_p: cobj,
    p_del_p: cobj,
    p_rdel_p: cobj,
    p_rdst_p: cobj,
    t_p: cobj,
    piln_p: cobj,
    rhoi_p: cobj,
    ni_p: cobj,
    ubm_p: cobj,
    ubi_p: cobj,
    xv_p: cobj,
    yv_p: cobj,
    effgw_p: cobj,
    c_p: cobj,
    kvtt_p: cobj,
    src_level_p: cobj,
    tend_level_p: cobj,
    tau_p: cobj,
    utgw_p: cobj,
    vtgw_p: cobj,
    ttgw_p: cobj,
    gwut_p: cobj,
    dttdf_p: cobj,
    dttke_p: cobj,
    d_p: cobj,
    mi_p: cobj,
    taudmp_p: cobj,
    tausat_p: cobj,
    ubmc_p: cobj,
    ubmc2_p: cobj,
    ubt_p: cobj,
    ubtl_p: cobj,
    wrk_p: cobj,
    ubt_lim_ratio_p: cobj,
    nm_p: cobj,
    q_p: cobj,
    dse_p: cobj,
    egwdffi_p: cobj,
    qtgw_p: cobj,
    egwdffm_p: cobj,
    egwdff_lev_p: cobj,
    dpidz_sq_p: cobj,
    coef_q_diff_p: cobj,
    qnew_p: cobj,
    spr_p: cobj,
    sub_p: cobj,
    diag_p: cobj,
    ca_p: cobj,
    ze_p: cobj,
    dnom_p: cobj,
    zf_p: cobj,
):
    if group == 1:
        _gw_drag_prof_core_impl(
            stage, ncol, pver, pverp, ngwv, ktop, kbot_tend, kbot_src, tau_0_ubc,
            dback, taumin, tndmax, umcfac, ubmc2mn, effkwv, kwv, gravit, rog, dt,
            alpha_p, p_del_p, p_rdel_p, t_p, piln_p, rhoi_p, ni_p, ubm_p, ubi_p,
            xv_p, yv_p, effgw_p, c_p, kvtt_p, src_level_p, tend_level_p, tau_p,
            utgw_p, vtgw_p, ttgw_p, gwut_p, dttdf_p, dttke_p, d_p, mi_p,
            taudmp_p, tausat_p, ubmc_p, ubmc2_p, ubt_p, ubtl_p, wrk_p,
            ubt_lim_ratio_p,
        )
    elif group == 2:
        _gw_diff_solver_impl(
            stage, ncol, pver, pverp, pcnst, ngwv, kbot, ktop, dt, gravit, gwut_p,
            ubm_p, nm_p, rhoi_p, c_p, tend_level_p, p_del_p, p_rdel_p, p_rdst_p,
            q_p, dse_p, egwdffi_p, qtgw_p, dttdf_p, egwdffm_p, egwdff_lev_p,
            dpidz_sq_p, coef_q_diff_p, qnew_p, spr_p, sub_p, diag_p, ca_p, ze_p,
            dnom_p, zf_p,
        )


@export
def gw_drag_prof_core_codon(
    stage: int,
    ncol: int,
    pver: int,
    pverp: int,
    ngwv: int,
    ktop: int,
    kbot_tend: int,
    kbot_src: int,
    tau_0_ubc: int,
    dback: float,
    taumin: float,
    tndmax: float,
    umcfac: float,
    ubmc2mn: float,
    effkwv: float,
    kwv: float,
    gravit: float,
    rog: float,
    dt: float,
    alpha_p: cobj,
    p_del_p: cobj,
    p_rdel_p: cobj,
    t_p: cobj,
    piln_p: cobj,
    rhoi_p: cobj,
    ni_p: cobj,
    ubm_p: cobj,
    ubi_p: cobj,
    xv_p: cobj,
    yv_p: cobj,
    effgw_p: cobj,
    c_p: cobj,
    kvtt_p: cobj,
    src_level_p: cobj,
    tend_level_p: cobj,
    tau_p: cobj,
    utgw_p: cobj,
    vtgw_p: cobj,
    ttgw_p: cobj,
    gwut_p: cobj,
    dttdf_p: cobj,
    dttke_p: cobj,
    d_p: cobj,
    mi_p: cobj,
    taudmp_p: cobj,
    tausat_p: cobj,
    ubmc_p: cobj,
    ubmc2_p: cobj,
    ubt_p: cobj,
    ubtl_p: cobj,
    wrk_p: cobj,
    ubt_lim_ratio_p: cobj,
):
    _gw_common_driver_dispatch(
        1, stage, ncol, pver, pverp, 0, ngwv, ktop, kbot_tend, kbot_src, 0,
        tau_0_ubc, dback, taumin, tndmax, umcfac, ubmc2mn, effkwv, kwv, gravit,
        rog, dt, alpha_p, p_del_p, p_rdel_p, p_rdel_p, t_p, piln_p, rhoi_p, ni_p,
        ubm_p, ubi_p, xv_p, yv_p, effgw_p, c_p, kvtt_p, src_level_p, tend_level_p,
        tau_p, utgw_p, vtgw_p, ttgw_p, gwut_p, dttdf_p, dttke_p, d_p, mi_p,
        taudmp_p, tausat_p, ubmc_p, ubmc2_p, ubt_p, ubtl_p, wrk_p,
        ubt_lim_ratio_p, ni_p, dttdf_p, dttdf_p, dttdf_p, dttdf_p, dttdf_p, d_p,
        dttdf_p, dttdf_p, dttdf_p, dttdf_p, dttdf_p, dttdf_p, dttdf_p, dttdf_p,
        dttdf_p, dttdf_p,
    )


@export
def gw_diff_solver_codon(
    stage: int,
    ncol: int,
    pver: int,
    pverp: int,
    pcnst: int,
    ngwv: int,
    kbot: int,
    ktop: int,
    dt: float,
    gravit: float,
    gwut_p: cobj,
    ubm_p: cobj,
    nm_p: cobj,
    rho_p: cobj,
    c_p: cobj,
    tend_level_p: cobj,
    p_del_p: cobj,
    p_rdel_p: cobj,
    p_rdst_p: cobj,
    q_p: cobj,
    dse_p: cobj,
    egwdffi_p: cobj,
    qtgw_p: cobj,
    dttdf_p: cobj,
    egwdffm_p: cobj,
    egwdff_lev_p: cobj,
    dpidz_sq_p: cobj,
    coef_q_diff_p: cobj,
    qnew_p: cobj,
    spr_p: cobj,
    sub_p: cobj,
    diag_p: cobj,
    ca_p: cobj,
    ze_p: cobj,
    dnom_p: cobj,
    zf_p: cobj,
):
    _gw_common_driver_dispatch(
        2, stage, ncol, pver, pverp, pcnst, ngwv, ktop, 0, 0, kbot, 0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, gravit, 0.0, dt, dse_p, p_del_p,
        p_rdel_p, p_rdst_p, dse_p, rho_p, rho_p, nm_p, ubm_p, ubm_p, egwdff_lev_p,
        egwdff_lev_p, egwdff_lev_p, c_p, egwdffi_p, tend_level_p, tend_level_p,
        egwdffi_p, dttdf_p, dttdf_p, dttdf_p, gwut_p, dttdf_p, dttdf_p,
        egwdff_lev_p, egwdff_lev_p, egwdff_lev_p, egwdff_lev_p, egwdff_lev_p,
        egwdff_lev_p, qnew_p, egwdff_lev_p, egwdff_lev_p, egwdff_lev_p, nm_p,
        q_p, dse_p, egwdffi_p, qtgw_p, egwdffm_p, egwdff_lev_p, dpidz_sq_p,
        coef_q_diff_p, qnew_p, spr_p, sub_p, diag_p, ca_p, ze_p, dnom_p, zf_p,
    )


@export
def gw_oro_src_codon(
    ncol: int,
    pver: int,
    ngwv: int,
    fcrit2: float,
    kwv: float,
    rair: float,
    p_mid_p: cobj,
    p_del_p: cobj,
    p_ifc_p: cobj,
    u_p: cobj,
    v_p: cobj,
    t_p: cobj,
    sgh_p: cobj,
    zm_p: cobj,
    nm_p: cobj,
    src_level_p: cobj,
    tend_level_p: cobj,
    tau_p: cobj,
    ubm_p: cobj,
    ubi_p: cobj,
    xv_p: cobj,
    yv_p: cobj,
    c_p: cobj,
    hdsp_p: cobj,
    tauoro_p: cobj,
    nsrc_p: cobj,
    rsrc_p: cobj,
    usrc_p: cobj,
    vsrc_p: cobj,
    dpsrc_p: cobj,
):
    p_mid = Ptr[float](p_mid_p)
    p_del = Ptr[float](p_del_p)
    p_ifc = Ptr[float](p_ifc_p)
    u = Ptr[float](u_p)
    v = Ptr[float](v_p)
    t = Ptr[float](t_p)
    sgh = Ptr[float](sgh_p)
    zm = Ptr[float](zm_p)
    nm = Ptr[float](nm_p)
    src_level = Ptr[int](src_level_p)
    tend_level = Ptr[int](tend_level_p)
    tau = Ptr[float](tau_p)
    ubm = Ptr[float](ubm_p)
    ubi = Ptr[float](ubi_p)
    xv = Ptr[float](xv_p)
    yv = Ptr[float](yv_p)
    c = Ptr[float](c_p)
    hdsp = Ptr[float](hdsp_p)
    tauoro = Ptr[float](tauoro_p)
    nsrc = Ptr[float](nsrc_p)
    rsrc = Ptr[float](rsrc_p)
    usrc = Ptr[float](usrc_p)
    vsrc = Ptr[float](vsrc_p)
    dpsrc = Ptr[float](dpsrc_p)

    orohmin = 10.0
    orovmin = 2.0

    for i in range(1, ncol + 1):
        hdsp[i - 1] = 2.0 * sgh[i - 1]

    k = pver
    for i in range(1, ncol + 1):
        src_level[i - 1] = k - 1
        rsrc[i - 1] = (p_mid[_idx2(i, k, ncol)] / (rair * t[_idx2(i, k, ncol)])) * p_del[_idx2(i, k, ncol)]
        usrc[i - 1] = u[_idx2(i, k, ncol)] * p_del[_idx2(i, k, ncol)]
        vsrc[i - 1] = v[_idx2(i, k, ncol)] * p_del[_idx2(i, k, ncol)]
        nsrc[i - 1] = nm[_idx2(i, k, ncol)] * p_del[_idx2(i, k, ncol)]

    for k in range(pver - 1, 0, -1):
        for i in range(1, ncol + 1):
            if hdsp[i - 1] > sqrt(zm[_idx2(i, k, ncol)] * zm[_idx2(i, k + 1, ncol)]):
                src_level[i - 1] = k - 1
                rsrc[i - 1] = rsrc[i - 1] + (
                    p_mid[_idx2(i, k, ncol)] / (rair * t[_idx2(i, k, ncol)])
                ) * p_del[_idx2(i, k, ncol)]
                usrc[i - 1] = usrc[i - 1] + u[_idx2(i, k, ncol)] * p_del[_idx2(i, k, ncol)]
                vsrc[i - 1] = vsrc[i - 1] + v[_idx2(i, k, ncol)] * p_del[_idx2(i, k, ncol)]
                nsrc[i - 1] = nsrc[i - 1] + nm[_idx2(i, k, ncol)] * p_del[_idx2(i, k, ncol)]

        all_found = True
        for i in range(1, ncol + 1):
            if src_level[i - 1] < k:
                all_found = False
        if all_found:
            break

    for i in range(1, ncol + 1):
        dpsrc[i - 1] = p_ifc[_idx2(i, pver + 1, ncol)] - p_ifc[_idx2(i, src_level[i - 1] + 1, ncol)]

    for i in range(1, ncol + 1):
        rsrc[i - 1] = rsrc[i - 1] / dpsrc[i - 1]
        usrc[i - 1] = usrc[i - 1] / dpsrc[i - 1]
        vsrc[i - 1] = vsrc[i - 1] / dpsrc[i - 1]
        nsrc[i - 1] = nsrc[i - 1] / dpsrc[i - 1]

    for i in range(1, ncol + 1):
        ubi[_idx2(i, pver + 1, ncol)] = sqrt(usrc[i - 1] * usrc[i - 1] + vsrc[i - 1] * vsrc[i - 1])

    for i in range(1, ncol + 1):
        mag = ubi[_idx2(i, pver + 1, ncol)]
        if mag > 0.0:
            xv[i - 1] = usrc[i - 1] / mag
            yv[i - 1] = vsrc[i - 1] / mag
        else:
            xv[i - 1] = 0.0
            yv[i - 1] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            ubm[_idx2(i, k, ncol)] = u[_idx2(i, k, ncol)] * xv[i - 1] + v[_idx2(i, k, ncol)] * yv[i - 1]

    for i in range(1, ncol + 1):
        ubi[_idx2(i, 1, ncol)] = ubm[_idx2(i, 1, ncol)]

    for k in range(2, pver + 1):
        for i in range(1, ncol + 1):
            ubi[_idx2(i, k, ncol)] = 0.5 * (ubm[_idx2(i, k - 1, ncol)] + ubm[_idx2(i, k, ncol)])

    for i in range(1, ncol + 1):
        ubi_surface = ubi[_idx2(i, pver + 1, ncol)]
        if (ubi_surface > orovmin) and (hdsp[i - 1] > orohmin):
            sghmax = fcrit2 * ((ubi_surface / nsrc[i - 1]) ** 2)
            tauoro[i - 1] = 0.5 * kwv * min(hdsp[i - 1] ** 2, sghmax) * rsrc[i - 1] * nsrc[i - 1] * ubi_surface
        else:
            tauoro[i - 1] = 0.0
            src_level[i - 1] = pver

    for k in range(1, pver + 2):
        for l in range(-ngwv, ngwv + 1):
            for i in range(1, ncol + 1):
                tau[_idx_tau(i, l, k, ncol, ngwv)] = 0.0

    min_src_level = src_level[0]
    for i in range(2, ncol + 1):
        if src_level[i - 1] < min_src_level:
            min_src_level = src_level[i - 1]

    for k in range(pver, min_src_level - 1, -1):
        for i in range(1, ncol + 1):
            if src_level[i - 1] <= k:
                tau[_idx_tau(i, 0, k + 1, ncol, ngwv)] = tauoro[i - 1]

    for i in range(1, ncol + 1):
        tend_level[i - 1] = pver

    for l in range(-ngwv, ngwv + 1):
        for i in range(1, ncol + 1):
            c[_idx_c(i, l, ncol, ngwv)] = 0.0


@export
def gw_oro_src_stage_dispatch_codon(
    ncol: int,
    pver: int,
    ngwv: int,
    fcrit2: float,
    kwv: float,
    rair: float,
    p_mid_p: cobj,
    p_del_p: cobj,
    p_ifc_p: cobj,
    u_p: cobj,
    v_p: cobj,
    t_p: cobj,
    sgh_p: cobj,
    zm_p: cobj,
    nm_p: cobj,
    src_level_p: cobj,
    tend_level_p: cobj,
    tau_p: cobj,
    ubm_p: cobj,
    ubi_p: cobj,
    xv_p: cobj,
    yv_p: cobj,
    c_p: cobj,
    hdsp_p: cobj,
    tauoro_p: cobj,
    nsrc_p: cobj,
    rsrc_p: cobj,
    usrc_p: cobj,
    vsrc_p: cobj,
    dpsrc_p: cobj,
):
    gw_oro_src_codon(
        ncol,
        pver,
        ngwv,
        fcrit2,
        kwv,
        rair,
        p_mid_p,
        p_del_p,
        p_ifc_p,
        u_p,
        v_p,
        t_p,
        sgh_p,
        zm_p,
        nm_p,
        src_level_p,
        tend_level_p,
        tau_p,
        ubm_p,
        ubi_p,
        xv_p,
        yv_p,
        c_p,
        hdsp_p,
        tauoro_p,
        nsrc_p,
        rsrc_p,
        usrc_p,
        vsrc_p,
        dpsrc_p,
    )


def _gw_tend_prep_impl(
    stage: int,
    ncol: int,
    psetcols: int,
    pcols: int,
    pver: int,
    pverp: int,
    pcnst: int,
    effgw_oro: float,
    eps: float,
    state_s_p: cobj,
    state_t_p: cobj,
    state_u_p: cobj,
    state_v_p: cobj,
    state_q_p: cobj,
    state_lnpint_p: cobj,
    state_zm_p: cobj,
    dse_p: cobj,
    t_p: cobj,
    u_p: cobj,
    v_p: cobj,
    q_p: cobj,
    piln_p: cobj,
    zm_p: cobj,
    egwdffi_tot_p: cobj,
    flx_heat_p: cobj,
    landfrac_p: cobj,
    sgh_p: cobj,
    effgw_p: cobj,
    sgh_scaled_p: cobj,
):
    state_s = Ptr[float](state_s_p)
    state_t = Ptr[float](state_t_p)
    state_u = Ptr[float](state_u_p)
    state_v = Ptr[float](state_v_p)
    state_q = Ptr[float](state_q_p)
    state_lnpint = Ptr[float](state_lnpint_p)
    state_zm = Ptr[float](state_zm_p)
    dse = Ptr[float](dse_p)
    t = Ptr[float](t_p)
    u = Ptr[float](u_p)
    v = Ptr[float](v_p)
    q = Ptr[float](q_p)
    piln = Ptr[float](piln_p)
    zm = Ptr[float](zm_p)
    egwdffi_tot = Ptr[float](egwdffi_tot_p)
    flx_heat = Ptr[float](flx_heat_p)
    landfrac = Ptr[float](landfrac_p)
    sgh = Ptr[float](sgh_p)
    effgw = Ptr[float](effgw_p)
    sgh_scaled = Ptr[float](sgh_scaled_p)

    if stage == 1:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                src_idx = _idx2(i, k, psetcols)
                dst_idx = _idx2(i, k, ncol)
                dse[dst_idx] = state_s[src_idx]
                t[dst_idx] = state_t[src_idx]
                u[dst_idx] = state_u[src_idx]
                v[dst_idx] = state_v[src_idx]
                zm[dst_idx] = state_zm[src_idx]

        for k in range(1, pverp + 1):
            for i in range(1, ncol + 1):
                piln[_idx2(i, k, ncol)] = state_lnpint[_idx2(i, k, psetcols)]
                egwdffi_tot[_idx2(i, k, ncol)] = 0.0

        for m in range(1, pcnst + 1):
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    q[_idx3(i, k, m, ncol, pver)] = state_q[_idx3(i, k, m, psetcols, pver)]

        for i in range(1, pcols + 1):
            flx_heat[i - 1] = 0.0

    elif stage == 2:
        for i in range(1, ncol + 1):
            landfrac_i = landfrac[i - 1]
            if landfrac_i >= eps:
                effgw[i - 1] = effgw_oro * landfrac_i
                sgh_scaled[i - 1] = sgh[i - 1] / sqrt(landfrac_i)
            else:
                effgw[i - 1] = 0.0
                sgh_scaled[i - 1] = 0.0


def _gw_tend_history_prep_impl(
    stage: int,
    ncol: int,
    psetcols: int,
    pcols: int,
    pver: int,
    cpair: float,
    dttdf_p: cobj,
    dttke_p: cobj,
    ptend_s_p: cobj,
    cpairv_p: cobj,
    ttgwsdf_oro_p: cobj,
    ttgwske_oro_p: cobj,
    ttgw_total_p: cobj,
):
    dttdf = Ptr[float](dttdf_p)
    dttke = Ptr[float](dttke_p)
    ptend_s = Ptr[float](ptend_s_p)
    cpairv = Ptr[float](cpairv_p)
    ttgwsdf_oro = Ptr[float](ttgwsdf_oro_p)
    ttgwske_oro = Ptr[float](ttgwske_oro_p)
    ttgw_total = Ptr[float](ttgw_total_p)

    if stage == 1:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, ncol)
                ttgwsdf_oro[idx] = dttdf[idx] / cpair
                ttgwske_oro[idx] = dttke[idx] / cpair

    elif stage == 2:
        for k in range(1, pver + 1):
            for i in range(1, pcols + 1):
                ttgw_total[_idx2(i, k, pcols)] = ptend_s[_idx2(i, k, psetcols)] / cpairv[_idx2(i, k, pcols)]


def _gw_tend_oro_post_impl(
    stage: int,
    ncol: int,
    psetcols: int,
    pcols: int,
    pver: int,
    pcnst: int,
    pverp: int,
    egwdffi_tot_p: cobj,
    egwdffi_p: cobj,
    utgw_p: cobj,
    vtgw_p: cobj,
    ttgw_p: cobj,
    cpairv_p: cobj,
    ptend_u_p: cobj,
    ptend_v_p: cobj,
    ptend_s_p: cobj,
    qtgw_p: cobj,
    ptend_q_p: cobj,
    tau_p: cobj,
    xv_p: cobj,
    yv_p: cobj,
    tau0x_p: cobj,
    tau0y_p: cobj,
):
    egwdffi_tot = Ptr[float](egwdffi_tot_p)
    egwdffi = Ptr[float](egwdffi_p)
    utgw = Ptr[float](utgw_p)
    vtgw = Ptr[float](vtgw_p)
    ttgw = Ptr[float](ttgw_p)
    cpairv = Ptr[float](cpairv_p)
    ptend_u = Ptr[float](ptend_u_p)
    ptend_v = Ptr[float](ptend_v_p)
    ptend_s = Ptr[float](ptend_s_p)
    qtgw = Ptr[float](qtgw_p)
    ptend_q = Ptr[float](ptend_q_p)
    tau = Ptr[float](tau_p)
    xv = Ptr[float](xv_p)
    yv = Ptr[float](yv_p)
    tau0x = Ptr[float](tau0x_p)
    tau0y = Ptr[float](tau0y_p)

    if stage == 1:
        for k in range(1, pverp + 1):
            for i in range(1, ncol + 1):
                idx_local = _idx2(i, k, ncol)
                egwdffi_tot[idx_local] = egwdffi_tot[idx_local] + egwdffi[idx_local]

        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx_local = _idx2(i, k, ncol)
                idx_ptend = _idx2(i, k, psetcols)
                ptend_u[idx_ptend] = ptend_u[idx_ptend] + utgw[idx_local]
                ptend_v[idx_ptend] = ptend_v[idx_ptend] + vtgw[idx_local]
                ptend_s[idx_ptend] = ptend_s[idx_ptend] + ttgw[idx_local]
                ttgw[idx_local] = ttgw[idx_local] / cpairv[_idx2(i, k, pcols)]

    elif stage == 2:
        for m in range(1, pcnst + 1):
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    ptend_idx = _idx3(i, k, m, psetcols, pver)
                    local_idx = _idx3(i, k, m, ncol, pver)
                    ptend_q[ptend_idx] = ptend_q[ptend_idx] + qtgw[local_idx]

        tau_surface_offset = pver * ncol
        for i in range(1, ncol + 1):
            tau_surface = tau[(i - 1) + tau_surface_offset]
            tau0x[i - 1] = tau_surface * xv[i - 1]
            tau0y[i - 1] = tau_surface * yv[i - 1]


def _gw_tend_driver_dispatch(
    group: int,
    stage: int,
    ncol: int,
    psetcols: int,
    pcols: int,
    pver: int,
    pverp: int,
    pcnst: int,
    effgw_oro: float,
    eps: float,
    cpair: float,
    state_s_p: cobj,
    state_t_p: cobj,
    state_u_p: cobj,
    state_v_p: cobj,
    state_q_p: cobj,
    state_lnpint_p: cobj,
    state_zm_p: cobj,
    dse_p: cobj,
    t_p: cobj,
    u_p: cobj,
    v_p: cobj,
    q_p: cobj,
    piln_p: cobj,
    zm_p: cobj,
    egwdffi_tot_p: cobj,
    flx_heat_p: cobj,
    landfrac_p: cobj,
    sgh_p: cobj,
    effgw_p: cobj,
    sgh_scaled_p: cobj,
    dttdf_p: cobj,
    dttke_p: cobj,
    ptend_s_p: cobj,
    cpairv_p: cobj,
    ttgwsdf_oro_p: cobj,
    ttgwske_oro_p: cobj,
    ttgw_total_p: cobj,
    egwdffi_p: cobj,
    utgw_p: cobj,
    vtgw_p: cobj,
    ttgw_p: cobj,
    ptend_u_p: cobj,
    ptend_v_p: cobj,
    qtgw_p: cobj,
    ptend_q_p: cobj,
    tau_p: cobj,
    xv_p: cobj,
    yv_p: cobj,
    tau0x_p: cobj,
    tau0y_p: cobj,
):
    if group == 1:
        _gw_tend_prep_impl(
            stage, ncol, psetcols, pcols, pver, pverp, pcnst, effgw_oro, eps,
            state_s_p, state_t_p, state_u_p, state_v_p, state_q_p, state_lnpint_p,
            state_zm_p, dse_p, t_p, u_p, v_p, q_p, piln_p, zm_p, egwdffi_tot_p,
            flx_heat_p, landfrac_p, sgh_p, effgw_p, sgh_scaled_p,
        )
    elif group == 2:
        _gw_tend_history_prep_impl(
            stage, ncol, psetcols, pcols, pver, cpair, dttdf_p, dttke_p,
            ptend_s_p, cpairv_p, ttgwsdf_oro_p, ttgwske_oro_p, ttgw_total_p,
        )
    elif group == 3:
        _gw_tend_oro_post_impl(
            stage, ncol, psetcols, pcols, pver, pcnst, pverp, egwdffi_tot_p,
            egwdffi_p, utgw_p, vtgw_p, ttgw_p, cpairv_p, ptend_u_p, ptend_v_p,
            ptend_s_p, qtgw_p, ptend_q_p, tau_p, xv_p, yv_p, tau0x_p, tau0y_p,
        )


@export
def gw_tend_prep_codon(
    stage: int,
    ncol: int,
    psetcols: int,
    pcols: int,
    pver: int,
    pverp: int,
    pcnst: int,
    effgw_oro: float,
    eps: float,
    state_s_p: cobj,
    state_t_p: cobj,
    state_u_p: cobj,
    state_v_p: cobj,
    state_q_p: cobj,
    state_lnpint_p: cobj,
    state_zm_p: cobj,
    dse_p: cobj,
    t_p: cobj,
    u_p: cobj,
    v_p: cobj,
    q_p: cobj,
    piln_p: cobj,
    zm_p: cobj,
    egwdffi_tot_p: cobj,
    flx_heat_p: cobj,
    landfrac_p: cobj,
    sgh_p: cobj,
    effgw_p: cobj,
    sgh_scaled_p: cobj,
):
    _gw_tend_driver_dispatch(
        1, stage, ncol, psetcols, pcols, pver, pverp, pcnst, effgw_oro, eps, 0.0,
        state_s_p, state_t_p, state_u_p, state_v_p, state_q_p, state_lnpint_p,
        state_zm_p, dse_p, t_p, u_p, v_p, q_p, piln_p, zm_p, egwdffi_tot_p,
        flx_heat_p, landfrac_p, sgh_p, effgw_p, sgh_scaled_p, dse_p, dse_p,
        dse_p, dse_p, dse_p, dse_p, dse_p, egwdffi_tot_p, t_p, u_p, t_p,
        dse_p, dse_p, q_p, q_p, egwdffi_tot_p, effgw_p, sgh_scaled_p, effgw_p,
        sgh_scaled_p,
    )


@export
def gw_tend_history_prep_codon(
    stage: int,
    ncol: int,
    psetcols: int,
    pcols: int,
    pver: int,
    cpair: float,
    dttdf_p: cobj,
    dttke_p: cobj,
    ptend_s_p: cobj,
    cpairv_p: cobj,
    ttgwsdf_oro_p: cobj,
    ttgwske_oro_p: cobj,
    ttgw_total_p: cobj,
):
    _gw_tend_driver_dispatch(
        2, stage, ncol, psetcols, pcols, pver, 0, 0, 0.0, 0.0, cpair,
        dttdf_p, dttdf_p, dttdf_p, dttdf_p, dttdf_p, dttdf_p, dttdf_p, dttdf_p,
        dttdf_p, dttdf_p, dttdf_p, dttdf_p, dttdf_p, dttdf_p, dttdf_p, dttdf_p,
        dttdf_p, dttdf_p, dttdf_p, dttdf_p, dttdf_p, dttke_p, ptend_s_p,
        cpairv_p, ttgwsdf_oro_p, ttgwske_oro_p, ttgw_total_p, dttdf_p, dttdf_p,
        dttdf_p, dttdf_p, dttdf_p, dttdf_p, dttdf_p, dttdf_p, dttdf_p, dttdf_p,
        dttdf_p, dttdf_p, dttdf_p,
    )


@export
def gw_tend_oro_post_codon(
    stage: int,
    ncol: int,
    psetcols: int,
    pcols: int,
    pver: int,
    pcnst: int,
    pverp: int,
    egwdffi_tot_p: cobj,
    egwdffi_p: cobj,
    utgw_p: cobj,
    vtgw_p: cobj,
    ttgw_p: cobj,
    cpairv_p: cobj,
    ptend_u_p: cobj,
    ptend_v_p: cobj,
    ptend_s_p: cobj,
    qtgw_p: cobj,
    ptend_q_p: cobj,
    tau_p: cobj,
    xv_p: cobj,
    yv_p: cobj,
    tau0x_p: cobj,
    tau0y_p: cobj,
):
    _gw_tend_driver_dispatch(
        3, stage, ncol, psetcols, pcols, pver, pverp, pcnst, 0.0, 0.0, 0.0,
        egwdffi_tot_p, egwdffi_tot_p, egwdffi_tot_p, egwdffi_tot_p, qtgw_p,
        egwdffi_tot_p, egwdffi_tot_p, egwdffi_tot_p, ttgw_p, utgw_p, vtgw_p,
        qtgw_p, egwdffi_tot_p, egwdffi_tot_p, egwdffi_tot_p, egwdffi_tot_p,
        egwdffi_tot_p, egwdffi_tot_p, egwdffi_tot_p, egwdffi_tot_p, ttgw_p,
        ttgw_p, ptend_s_p, cpairv_p, ttgw_p, ttgw_p, ttgw_p, egwdffi_p, utgw_p,
        vtgw_p, ttgw_p, ptend_u_p, ptend_v_p, qtgw_p, ptend_q_p, tau_p, xv_p,
        yv_p, tau0x_p, tau0y_p,
    )
