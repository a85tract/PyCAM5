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
