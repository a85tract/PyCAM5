from math import exp, log, sqrt


@export
def rrtmg_init_real_passthrough_codon(value: float) -> float:
    return value


@export
def rrtmg_init_int_passthrough_codon(value: int) -> int:
    return value


@inline
def _idx2(i: int, k: int, ld1: int) -> int:
    """Fortran array declared as (ld1, n2), both bounds starting at 1."""
    return (i - 1) + (k - 1) * ld1


@inline
def _idx2_lb0(i: int, k0: int, ld1: int) -> int:
    """Fortran array declared as (ld1, 0:n2)."""
    return (i - 1) + k0 * ld1


@inline
def _idx2_dim1_lb0(k0: int, b: int, ub1: int) -> int:
    """Fortran array declared as (0:ub1, n2)."""
    return k0 + (b - 1) * (ub1 + 1)


@inline
def _idx2_dim2_lb(i: int, b: int, ld1: int, lb2: int) -> int:
    """Fortran array declared as (ld1, lb2:ub2)."""
    return (i - 1) + (b - lb2) * ld1


@inline
def _idx3(a: int, b: int, c: int, ld1: int, ld2: int) -> int:
    """Fortran array declared as (ld1, ld2, n3), bounds starting at 1."""
    return (a - 1) + (b - 1) * ld1 + (c - 1) * ld1 * ld2


@inline
def _lw_tau_major_3pt(
    specparm: float,
    speccomb: float,
    ind: int,
    ig: int,
    fac000: float,
    fac100: float,
    fac200: float,
    fac010: float,
    fac110: float,
    fac210: float,
    coeff: Ptr[float],
    ld1: int,
) -> float:
    if specparm < 0.125:
        return speccomb * (
            fac000 * coeff[_idx2(ind, ig, ld1)]
            + fac100 * coeff[_idx2(ind + 1, ig, ld1)]
            + fac200 * coeff[_idx2(ind + 2, ig, ld1)]
            + fac010 * coeff[_idx2(ind + 9, ig, ld1)]
            + fac110 * coeff[_idx2(ind + 10, ig, ld1)]
            + fac210 * coeff[_idx2(ind + 11, ig, ld1)]
        )
    if specparm > 0.875:
        return speccomb * (
            fac200 * coeff[_idx2(ind - 1, ig, ld1)]
            + fac100 * coeff[_idx2(ind, ig, ld1)]
            + fac000 * coeff[_idx2(ind + 1, ig, ld1)]
            + fac210 * coeff[_idx2(ind + 8, ig, ld1)]
            + fac110 * coeff[_idx2(ind + 9, ig, ld1)]
            + fac010 * coeff[_idx2(ind + 10, ig, ld1)]
        )
    return speccomb * (
        fac000 * coeff[_idx2(ind, ig, ld1)]
        + fac100 * coeff[_idx2(ind + 1, ig, ld1)]
        + fac010 * coeff[_idx2(ind + 9, ig, ld1)]
        + fac110 * coeff[_idx2(ind + 10, ig, ld1)]
    )


@inline
def _idx3_lb0_dim2(a: int, b0: int, c: int, ld1: int, ub2: int) -> int:
    """Fortran array declared as (ld1, 0:ub2, n3)."""
    return (a - 1) + b0 * ld1 + (c - 1) * ld1 * (ub2 + 1)


@inline
def _max_1em80(x: float) -> float:
    if x > 1.0e-80:
        return x
    return 1.0e-80


@inline
def _interp1_bndry_table(yin: Ptr[float], table: Ptr[float], nin: int, band: int, y: float) -> float:
    increasing = True
    if yin[0] > yin[1]:
        increasing = False

    if increasing:
        if y <= yin[0]:
            value = table[_idx2(1, band, nin)]
            return value * 1.0 + value * 0.0
        if y > yin[nin - 1]:
            value = table[_idx2(nin, band, nin)]
            return value * 1.0 + value * 0.0
        for jj in range(1, nin):
            if y > yin[jj - 1] and y <= yin[jj]:
                denom = yin[jj] - yin[jj - 1]
                wgts = (yin[jj] - y) / denom
                wgtn = (y - yin[jj - 1]) / denom
                return table[_idx2(jj, band, nin)] * wgts + table[_idx2(jj + 1, band, nin)] * wgtn
    else:
        if y > yin[0]:
            value = table[_idx2(1, band, nin)]
            return value * 1.0 + value * 0.0
        if y <= yin[nin - 1]:
            value = table[_idx2(nin, band, nin)]
            return value * 1.0 + value * 0.0
        for jj in range(1, nin):
            if y <= yin[jj - 1] and y > yin[jj]:
                denom = yin[jj] - yin[jj - 1]
                wgts = (yin[jj] - y) / denom
                wgtn = (y - yin[jj - 1]) / denom
                return table[_idx2(jj, band, nin)] * wgts + table[_idx2(jj + 1, band, nin)] * wgtn

    return table[0] * 0.0


@inline
def _interp1_bndry_jjm(yin: Ptr[float], nin: int, y: float) -> int:
    increasing = True
    if yin[0] > yin[1]:
        increasing = False

    if increasing:
        if y <= yin[0]:
            return 1
        if y > yin[nin - 1]:
            return nin
        for jj in range(1, nin):
            if y > yin[jj - 1] and y <= yin[jj]:
                return jj
    else:
        if y > yin[0]:
            return 1
        if y <= yin[nin - 1]:
            return nin
        for jj in range(1, nin):
            if y <= yin[jj - 1] and y > yin[jj]:
                return jj

    return 1


@inline
def _interp1_bndry_jjp(yin: Ptr[float], nin: int, y: float) -> int:
    jjm = _interp1_bndry_jjm(yin, nin, y)
    if jjm == 1:
        increasing = True
        if yin[0] > yin[1]:
            increasing = False
        if (increasing and y <= yin[0]) or ((not increasing) and y > yin[0]):
            return 1
    if jjm == nin:
        increasing = True
        if yin[0] > yin[1]:
            increasing = False
        if (increasing and y > yin[nin - 1]) or ((not increasing) and y <= yin[nin - 1]):
            return nin
    return jjm + 1


@inline
def _interp1_wgts(yin: Ptr[float], jjm: int, jjp: int, y: float) -> float:
    if jjm == jjp:
        return 1.0
    return (yin[jjp - 1] - y) / (yin[jjp - 1] - yin[jjm - 1])


@inline
def _interp1_wgtn(yin: Ptr[float], jjm: int, jjp: int, y: float) -> float:
    if jjm == jjp:
        return 0.0
    return (y - yin[jjm - 1]) / (yin[jjp - 1] - yin[jjm - 1])


@inline
def _liquid_lambda_grid_value(
    g_lambda: Ptr[float],
    nmu: int,
    ilambda: int,
    mu_jjm: int,
    mu_jjp: int,
    mu_wgts: float,
    mu_wgtn: float,
) -> float:
    return (
        g_lambda[_idx2(mu_jjm, ilambda, nmu)] * mu_wgts
        + g_lambda[_idx2(mu_jjp, ilambda, nmu)] * mu_wgtn
    )


@inline
def _liquid_lambda_jjm(
    g_lambda: Ptr[float],
    nmu: int,
    nlambda: int,
    lamc: float,
    mu_jjm: int,
    mu_jjp: int,
    mu_wgts: float,
    mu_wgtn: float,
) -> int:
    first = _liquid_lambda_grid_value(g_lambda, nmu, 1, mu_jjm, mu_jjp, mu_wgts, mu_wgtn)
    second = _liquid_lambda_grid_value(g_lambda, nmu, 2, mu_jjm, mu_jjp, mu_wgts, mu_wgtn)
    last = _liquid_lambda_grid_value(g_lambda, nmu, nlambda, mu_jjm, mu_jjp, mu_wgts, mu_wgtn)
    increasing = True
    if first > second:
        increasing = False

    if increasing:
        if lamc <= first:
            return 1
        if lamc > last:
            return nlambda
        for jj in range(1, nlambda):
            lo = _liquid_lambda_grid_value(g_lambda, nmu, jj, mu_jjm, mu_jjp, mu_wgts, mu_wgtn)
            hi = _liquid_lambda_grid_value(g_lambda, nmu, jj + 1, mu_jjm, mu_jjp, mu_wgts, mu_wgtn)
            if lamc > lo and lamc <= hi:
                return jj
    else:
        if lamc > first:
            return 1
        if lamc <= last:
            return nlambda
        for jj in range(1, nlambda):
            lo = _liquid_lambda_grid_value(g_lambda, nmu, jj, mu_jjm, mu_jjp, mu_wgts, mu_wgtn)
            hi = _liquid_lambda_grid_value(g_lambda, nmu, jj + 1, mu_jjm, mu_jjp, mu_wgts, mu_wgtn)
            if lamc <= lo and lamc > hi:
                return jj

    return 1


@inline
def _liquid_lambda_jjp(
    g_lambda: Ptr[float],
    nmu: int,
    nlambda: int,
    lamc: float,
    mu_jjm: int,
    mu_jjp: int,
    mu_wgts: float,
    mu_wgtn: float,
) -> int:
    lam_jjm = _liquid_lambda_jjm(g_lambda, nmu, nlambda, lamc, mu_jjm, mu_jjp, mu_wgts, mu_wgtn)
    first = _liquid_lambda_grid_value(g_lambda, nmu, 1, mu_jjm, mu_jjp, mu_wgts, mu_wgtn)
    second = _liquid_lambda_grid_value(g_lambda, nmu, 2, mu_jjm, mu_jjp, mu_wgts, mu_wgtn)
    increasing = True
    if first > second:
        increasing = False
    if lam_jjm == 1 and ((increasing and lamc <= first) or ((not increasing) and lamc > first)):
        return 1
    last = _liquid_lambda_grid_value(g_lambda, nmu, nlambda, mu_jjm, mu_jjp, mu_wgts, mu_wgtn)
    if lam_jjm == nlambda and ((increasing and lamc > last) or ((not increasing) and lamc <= last)):
        return nlambda
    return lam_jjm + 1


@inline
def _liquid_interp2(
    table: Ptr[float],
    nmu: int,
    nlambda: int,
    band: int,
    mu_jjm: int,
    mu_jjp: int,
    mu_wgts: float,
    mu_wgtn: float,
    lam_jjm: int,
    lam_jjp: int,
    lam_wgts: float,
    lam_wgtn: float,
) -> float:
    return (
        table[_idx3(mu_jjm, lam_jjm, band, nmu, nlambda)] * mu_wgts * lam_wgts
        + table[_idx3(mu_jjp, lam_jjm, band, nmu, nlambda)] * mu_wgtn * lam_wgts
        + table[_idx3(mu_jjm, lam_jjp, band, nmu, nlambda)] * mu_wgts * lam_wgtn
        + table[_idx3(mu_jjp, lam_jjp, band, nmu, nlambda)] * mu_wgtn * lam_wgtn
    )


@inline
def _rrtmg_src_level(k: int, pverp: int, num_rrtmg_levs: int) -> int:
    kk = k + (pverp - num_rrtmg_levs) - 1
    if kk < 1:
        return 1
    return kk


@inline
def _tint_at(
    i: int,
    k: int,
    pcols: int,
    pverp: int,
    stebol: float,
    t: Ptr[float],
    lnpint: Ptr[float],
    lnpmid: Ptr[float],
    lwup: Ptr[float],
) -> float:
    if k == 1:
        return t[_idx2(i, 1, pcols)]
    if k == pverp:
        return sqrt(sqrt(lwup[i - 1] / stebol))

    dy = (
        (lnpint[_idx2(i, k, pcols)] - lnpmid[_idx2(i, k, pcols)])
        / (lnpmid[_idx2(i, k - 1, pcols)] - lnpmid[_idx2(i, k, pcols)])
    )
    return t[_idx2(i, k, pcols)] - dy * (
        t[_idx2(i, k, pcols)] - t[_idx2(i, k - 1, pcols)]
    )


@export
def rrtmg_state_create_fields_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    num_rrtmg_levs: int,
    stebol: float,
    t_p: cobj,
    lnpint_p: cobj,
    lnpmid_p: cobj,
    pmid_p: cobj,
    pint_p: cobj,
    lwup_p: cobj,
    pmidmb_p: cobj,
    pintmb_p: cobj,
    tlay_p: cobj,
    tlev_p: cobj,
):
    t = Ptr[float](t_p)
    lnpint = Ptr[float](lnpint_p)
    lnpmid = Ptr[float](lnpmid_p)
    pmid = Ptr[float](pmid_p)
    pint = Ptr[float](pint_p)
    lwup = Ptr[float](lwup_p)
    pmidmb = Ptr[float](pmidmb_p)
    pintmb = Ptr[float](pintmb_p)
    tlay = Ptr[float](tlay_p)
    tlev = Ptr[float](tlev_p)

    for k in range(1, num_rrtmg_levs + 1):
        kk = _rrtmg_src_level(k, pverp, num_rrtmg_levs)
        for i in range(1, ncol + 1):
            pmidmb[_idx2(i, k, pcols)] = pmid[_idx2(i, kk, pcols)] * 1.0e-2
            pintmb[_idx2(i, k, pcols)] = pint[_idx2(i, kk, pcols)] * 1.0e-2
            tlay[_idx2(i, k, pcols)] = t[_idx2(i, kk, pcols)]
            tlev[_idx2(i, k, pcols)] = _tint_at(
                i, kk, pcols, pverp, stebol, t, lnpint, lnpmid, lwup
            )

    for i in range(1, ncol + 1):
        pintmb[_idx2(i, num_rrtmg_levs + 1, pcols)] = (
            pint[_idx2(i, pverp, pcols)] * 1.0e-2
        )
        tlev[_idx2(i, num_rrtmg_levs + 1, pcols)] = _tint_at(
            i, pverp, pcols, pverp, stebol, t, lnpint, lnpmid, lwup
        )

    if num_rrtmg_levs == pverp:
        for i in range(1, ncol + 1):
            pmidmb[_idx2(i, 1, pcols)] = 0.5 * pintmb[_idx2(i, 2, pcols)]
            pintmb[_idx2(i, 1, pcols)] = 1.0e-4


@export
def rrtmg_state_update_codon(
    ncol: int,
    pcols: int,
    pverp: int,
    num_rrtmg_levs: int,
    sp_hum_p: cobj,
    o2_p: cobj,
    o3_p: cobj,
    co2_p: cobj,
    n2o_p: cobj,
    ch4_p: cobj,
    cfc11_p: cobj,
    cfc12_p: cobj,
    ch4vmr_p: cobj,
    h2ovmr_p: cobj,
    o3vmr_p: cobj,
    co2vmr_p: cobj,
    o2vmr_p: cobj,
    n2ovmr_p: cobj,
    cfc11vmr_p: cobj,
    cfc12vmr_p: cobj,
    cfc22vmr_p: cobj,
    ccl4vmr_p: cobj,
):
    sp_hum = Ptr[float](sp_hum_p)
    o2 = Ptr[float](o2_p)
    o3 = Ptr[float](o3_p)
    co2 = Ptr[float](co2_p)
    n2o = Ptr[float](n2o_p)
    ch4 = Ptr[float](ch4_p)
    cfc11 = Ptr[float](cfc11_p)
    cfc12 = Ptr[float](cfc12_p)
    ch4vmr = Ptr[float](ch4vmr_p)
    h2ovmr = Ptr[float](h2ovmr_p)
    o3vmr = Ptr[float](o3vmr_p)
    co2vmr = Ptr[float](co2vmr_p)
    o2vmr = Ptr[float](o2vmr_p)
    n2ovmr = Ptr[float](n2ovmr_p)
    cfc11vmr = Ptr[float](cfc11vmr_p)
    cfc12vmr = Ptr[float](cfc12vmr_p)
    cfc22vmr = Ptr[float](cfc22vmr_p)
    ccl4vmr = Ptr[float](ccl4vmr_p)

    amdw = 1.607793
    amdc = 0.658114
    amdo = 0.603428
    amdm = 1.805423
    amdn = 0.658090
    amdo2 = 0.905140
    amdc1 = 0.210852
    amdc2 = 0.239546

    for k in range(1, num_rrtmg_levs + 1):
        kk = _rrtmg_src_level(k, pverp, num_rrtmg_levs)
        for i in range(1, ncol + 1):
            src = _idx2(i, kk, pcols)
            dst = _idx2(i, k, pcols)
            ch4vmr[dst] = ch4[src] * amdm
            h2ovmr[dst] = (sp_hum[src] / (1.0 - sp_hum[src])) * amdw
            o3vmr[dst] = o3[src] * amdo
            co2vmr[dst] = co2[src] * amdc
            ch4vmr[dst] = ch4[src] * amdm
            o2vmr[dst] = o2[src] * amdo2
            n2ovmr[dst] = n2o[src] * amdn
            cfc11vmr[dst] = cfc11[src] * amdc1
            cfc12vmr[dst] = cfc12[src] * amdc2
            cfc22vmr[dst] = 0.0
            ccl4vmr[dst] = 0.0


@export
def rrtmg_lw_zero_cloud_inputs_codon(
    ncol: int,
    pcols: int,
    nlay: int,
    cicewp_p: cobj,
    cliqwp_p: cobj,
    rei_p: cobj,
    rel_p: cobj,
):
    cicewp = Ptr[float](cicewp_p)
    cliqwp = Ptr[float](cliqwp_p)
    rei = Ptr[float](rei_p)
    rel = Ptr[float](rel_p)

    # Fortran declarations are (pcols, rrtmg_levs-1).
    for k in range(1, nlay + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            cicewp[idx] = 0.0
            cliqwp[idx] = 0.0
            rei[idx] = 0.0
            rel[idx] = 0.0


@export
def rrtmg_lw_pre_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    rrtmg_levs: int,
    nbndlw: int,
    aer_lw_abs_p: cobj,
    tlev_p: cobj,
    emis_p: cobj,
    tsfc_p: cobj,
    taua_lw_p: cobj,
):
    aer_lw_abs = Ptr[float](aer_lw_abs_p)
    tlev = Ptr[float](tlev_p)
    emis = Ptr[float](emis_p)
    tsfc = Ptr[float](tsfc_p)
    taua_lw = Ptr[float](taua_lw_p)

    for nbnd in range(1, nbndlw + 1):
        for i in range(1, ncol + 1):
            emis[_idx2(i, nbnd, pcols)] = 1.0

    for i in range(1, ncol + 1):
        tsfc[i - 1] = tlev[_idx2(i, rrtmg_levs + 1, pcols)]

    for nbnd in range(1, nbndlw + 1):
        for k in range(1, rrtmg_levs):
            src_k = pverp - rrtmg_levs + k
            for i in range(1, ncol + 1):
                taua_lw[_idx3(i, k, nbnd, pcols, rrtmg_levs - 1)] = aer_lw_abs[
                    _idx3(i, src_k, nbnd, pcols, pver)
                ]


@export
def rrtmg_lw_post_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    rrtmg_levs: int,
    ntoplw: int,
    cpair: float,
    uflx_p: cobj,
    dflx_p: cobj,
    hr_p: cobj,
    uflxc_p: cobj,
    dflxc_p: cobj,
    hrc_p: cobj,
    flwds_p: cobj,
    fldsc_p: cobj,
    flns_p: cobj,
    flnsc_p: cobj,
    flnt_p: cobj,
    flntc_p: cobj,
    flut_p: cobj,
    flutc_p: cobj,
    ful_p: cobj,
    fdl_p: cobj,
    fsul_p: cobj,
    fsdl_p: cobj,
    fnl_p: cobj,
    fcnl_p: cobj,
    qrl_p: cobj,
    qrlc_p: cobj,
):
    uflx = Ptr[float](uflx_p)
    dflx = Ptr[float](dflx_p)
    hr = Ptr[float](hr_p)
    uflxc = Ptr[float](uflxc_p)
    dflxc = Ptr[float](dflxc_p)
    hrc = Ptr[float](hrc_p)
    flwds = Ptr[float](flwds_p)
    fldsc = Ptr[float](fldsc_p)
    flns = Ptr[float](flns_p)
    flnsc = Ptr[float](flnsc_p)
    flnt = Ptr[float](flnt_p)
    flntc = Ptr[float](flntc_p)
    flut = Ptr[float](flut_p)
    flutc = Ptr[float](flutc_p)
    ful = Ptr[float](ful_p)
    fdl = Ptr[float](fdl_p)
    fsul = Ptr[float](fsul_p)
    fsdl = Ptr[float](fsdl_p)
    fnl = Ptr[float](fnl_p)
    fcnl = Ptr[float](fcnl_p)
    qrl = Ptr[float](qrl_p)
    qrlc = Ptr[float](qrlc_p)

    dps = 1.0 / 86400.0

    for i in range(1, ncol + 1):
        flwds[i - 1] = dflx[_idx2(i, 1, pcols)]
        fldsc[i - 1] = dflxc[_idx2(i, 1, pcols)]
        flns[i - 1] = uflx[_idx2(i, 1, pcols)] - dflx[_idx2(i, 1, pcols)]
        flnsc[i - 1] = uflxc[_idx2(i, 1, pcols)] - dflxc[_idx2(i, 1, pcols)]
        flnt[i - 1] = (
            uflx[_idx2(i, rrtmg_levs, pcols)]
            - dflx[_idx2(i, rrtmg_levs, pcols)]
        )
        flntc[i - 1] = (
            uflxc[_idx2(i, rrtmg_levs, pcols)]
            - dflxc[_idx2(i, rrtmg_levs, pcols)]
        )
        flut[i - 1] = uflx[_idx2(i, rrtmg_levs, pcols)]
        flutc[i - 1] = uflxc[_idx2(i, rrtmg_levs, pcols)]

    for k in range(1, pverp + 1):
        for i in range(1, pcols + 1):
            ful[_idx2(i, k, pcols)] = 0.0
            fdl[_idx2(i, k, pcols)] = 0.0
            fsul[_idx2(i, k, pcols)] = 0.0
            fsdl[_idx2(i, k, pcols)] = 0.0

    for j in range(1, rrtmg_levs + 1):
        cam_k = pverp - rrtmg_levs + j
        rrtmg_k = rrtmg_levs - j + 1
        for i in range(1, ncol + 1):
            ful[_idx2(i, cam_k, pcols)] = uflx[_idx2(i, rrtmg_k, pcols)]
            fdl[_idx2(i, cam_k, pcols)] = dflx[_idx2(i, rrtmg_k, pcols)]
            fsul[_idx2(i, cam_k, pcols)] = uflxc[_idx2(i, rrtmg_k, pcols)]
            fsdl[_idx2(i, cam_k, pcols)] = dflxc[_idx2(i, rrtmg_k, pcols)]

    for k in range(1, pverp + 1):
        for i in range(1, ncol + 1):
            fnl[_idx2(i, k, pcols)] = ful[_idx2(i, k, pcols)] - fdl[_idx2(i, k, pcols)]
            fcnl[_idx2(i, k, pcols)] = fsul[_idx2(i, k, pcols)] - fsdl[_idx2(i, k, pcols)]

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            qrl[_idx2(i, k, pcols)] = 0.0
            qrlc[_idx2(i, k, pcols)] = 0.0

    for j in range(1, rrtmg_levs):
        cam_k = pverp - rrtmg_levs + j
        rrtmg_k = rrtmg_levs - j
        for i in range(1, ncol + 1):
            qrl[_idx2(i, cam_k, pcols)] = hr[_idx2(i, rrtmg_k, pcols)] * cpair * dps
            qrlc[_idx2(i, cam_k, pcols)] = hrc[_idx2(i, rrtmg_k, pcols)] * cpair * dps

    if ntoplw > 1:
        for k in range(1, ntoplw):
            for i in range(1, ncol + 1):
                qrl[_idx2(i, k, pcols)] = 0.0
                qrlc[_idx2(i, k, pcols)] = 0.0


@inline
def _radsw_compact_1d(
    src: Ptr[float],
    dst: Ptr[float],
    nday: int,
    nnite: int,
    idxday: Ptr[int],
    idxnite: Ptr[int],
):
    for pos in range(1, nday + 1):
        dst[pos - 1] = src[idxday[pos - 1] - 1]
    for pos in range(1, nnite + 1):
        dst[nday + pos - 1] = src[idxnite[pos - 1] - 1]


@inline
def _radsw_compact_2d_same_levels(
    src: Ptr[float],
    dst: Ptr[float],
    nday: int,
    nnite: int,
    nlev: int,
    pcols: int,
    idxday: Ptr[int],
    idxnite: Ptr[int],
):
    for k in range(1, nlev + 1):
        for pos in range(1, nday + 1):
            dst[_idx2(pos, k, pcols)] = src[_idx2(idxday[pos - 1], k, pcols)]
        for pos in range(1, nnite + 1):
            dst[_idx2(nday + pos, k, pcols)] = src[_idx2(idxnite[pos - 1], k, pcols)]


@inline
def _radsw_compact_2d_level_offset(
    src: Ptr[float],
    dst: Ptr[float],
    nday: int,
    nnite: int,
    nlev: int,
    src_offset: int,
    pcols: int,
    idxday: Ptr[int],
    idxnite: Ptr[int],
):
    for k in range(1, nlev + 1):
        src_k = src_offset + k
        for pos in range(1, nday + 1):
            dst[_idx2(pos, k, pcols)] = src[_idx2(idxday[pos - 1], src_k, pcols)]
        for pos in range(1, nnite + 1):
            dst[_idx2(nday + pos, k, pcols)] = src[_idx2(idxnite[pos - 1], src_k, pcols)]


@inline
def _rrtmg_zero_1d_ncol(arr: Ptr[float], ncol: int):
    for i in range(1, ncol + 1):
        arr[i - 1] = 0.0


@inline
def _rrtmg_zero_2d_ncol(arr: Ptr[float], ncol: int, nlev: int, pcols: int):
    for k in range(1, nlev + 1):
        for i in range(1, ncol + 1):
            arr[_idx2(i, k, pcols)] = 0.0


@export
def rrtmg_sw_zero_outputs_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    solin_p: cobj,
    qrs_p: cobj,
    qrsc_p: cobj,
    fns_p: cobj,
    fcns_p: cobj,
    fsds_p: cobj,
    fsnirtoa_p: cobj,
    fsnrtoac_p: cobj,
    fsnrtoaq_p: cobj,
    fsns_p: cobj,
    fsnsc_p: cobj,
    fsdsc_p: cobj,
    fsnt_p: cobj,
    fsntc_p: cobj,
    fsntoa_p: cobj,
    fsutoa_p: cobj,
    fsntoac_p: cobj,
    sols_p: cobj,
    soll_p: cobj,
    solsd_p: cobj,
    solld_p: cobj,
):
    solin = Ptr[float](solin_p)
    qrs = Ptr[float](qrs_p)
    qrsc = Ptr[float](qrsc_p)
    fns = Ptr[float](fns_p)
    fcns = Ptr[float](fcns_p)
    fsds = Ptr[float](fsds_p)
    fsnirtoa = Ptr[float](fsnirtoa_p)
    fsnrtoac = Ptr[float](fsnrtoac_p)
    fsnrtoaq = Ptr[float](fsnrtoaq_p)
    fsns = Ptr[float](fsns_p)
    fsnsc = Ptr[float](fsnsc_p)
    fsdsc = Ptr[float](fsdsc_p)
    fsnt = Ptr[float](fsnt_p)
    fsntc = Ptr[float](fsntc_p)
    fsntoa = Ptr[float](fsntoa_p)
    fsutoa = Ptr[float](fsutoa_p)
    fsntoac = Ptr[float](fsntoac_p)
    sols = Ptr[float](sols_p)
    soll = Ptr[float](soll_p)
    solsd = Ptr[float](solsd_p)
    solld = Ptr[float](solld_p)

    _rrtmg_zero_1d_ncol(fsds, ncol)
    _rrtmg_zero_1d_ncol(fsnirtoa, ncol)
    _rrtmg_zero_1d_ncol(fsnrtoac, ncol)
    _rrtmg_zero_1d_ncol(fsnrtoaq, ncol)
    _rrtmg_zero_1d_ncol(fsns, ncol)
    _rrtmg_zero_1d_ncol(fsnsc, ncol)
    _rrtmg_zero_1d_ncol(fsdsc, ncol)
    _rrtmg_zero_1d_ncol(fsnt, ncol)
    _rrtmg_zero_1d_ncol(fsntc, ncol)
    _rrtmg_zero_1d_ncol(fsntoa, ncol)
    _rrtmg_zero_1d_ncol(fsutoa, ncol)
    _rrtmg_zero_1d_ncol(fsntoac, ncol)
    _rrtmg_zero_1d_ncol(solin, ncol)
    _rrtmg_zero_1d_ncol(sols, ncol)
    _rrtmg_zero_1d_ncol(soll, ncol)
    _rrtmg_zero_1d_ncol(solsd, ncol)
    _rrtmg_zero_1d_ncol(solld, ncol)
    _rrtmg_zero_2d_ncol(qrs, ncol, pver, pcols)
    _rrtmg_zero_2d_ncol(qrsc, ncol, pver, pcols)
    _rrtmg_zero_2d_ncol(fns, ncol, pverp, pcols)
    _rrtmg_zero_2d_ncol(fcns, ncol, pverp, pcols)


@export
def rrtmg_sw_zero_cloud_inputs_codon(
    ncol: int,
    pcols: int,
    nlay: int,
    cicewp_p: cobj,
    cliqwp_p: cobj,
    rel_p: cobj,
    rei_p: cobj,
):
    cicewp = Ptr[float](cicewp_p)
    cliqwp = Ptr[float](cliqwp_p)
    rel = Ptr[float](rel_p)
    rei = Ptr[float](rei_p)

    # Fortran declarations are (pcols, rrtmg_levs-1).
    for k in range(1, nlay + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            cicewp[idx] = 0.0
            cliqwp[idx] = 0.0
            rel[idx] = 0.0
            rei[idx] = 0.0


@export
def rrtmg_sw_compact_inputs_codon(
    nday: int,
    nnite: int,
    pcols: int,
    pverp: int,
    rrtmg_levs: int,
    idxday_p: cobj,
    idxnite_p: cobj,
    e_pmid_p: cobj,
    e_cld_p: cobj,
    state_pintmb_p: cobj,
    state_pmidmb_p: cobj,
    state_h2ovmr_p: cobj,
    state_o3vmr_p: cobj,
    state_co2vmr_p: cobj,
    e_coszrs_p: cobj,
    e_asdir_p: cobj,
    e_aldir_p: cobj,
    e_asdif_p: cobj,
    e_aldif_p: cobj,
    state_tlay_p: cobj,
    state_tlev_p: cobj,
    state_ch4vmr_p: cobj,
    state_o2vmr_p: cobj,
    state_n2ovmr_p: cobj,
    pmid_p: cobj,
    cld_p: cobj,
    pintmb_p: cobj,
    pmidmb_p: cobj,
    h2ovmr_p: cobj,
    o3vmr_p: cobj,
    co2vmr_p: cobj,
    coszrs_p: cobj,
    asdir_p: cobj,
    aldir_p: cobj,
    asdif_p: cobj,
    aldif_p: cobj,
    tlay_p: cobj,
    tlev_p: cobj,
    ch4vmr_p: cobj,
    o2vmr_p: cobj,
    n2ovmr_p: cobj,
):
    idxday = Ptr[int](idxday_p)
    idxnite = Ptr[int](idxnite_p)
    e_pmid = Ptr[float](e_pmid_p)
    e_cld = Ptr[float](e_cld_p)
    state_pintmb = Ptr[float](state_pintmb_p)
    state_pmidmb = Ptr[float](state_pmidmb_p)
    state_h2ovmr = Ptr[float](state_h2ovmr_p)
    state_o3vmr = Ptr[float](state_o3vmr_p)
    state_co2vmr = Ptr[float](state_co2vmr_p)
    e_coszrs = Ptr[float](e_coszrs_p)
    e_asdir = Ptr[float](e_asdir_p)
    e_aldir = Ptr[float](e_aldir_p)
    e_asdif = Ptr[float](e_asdif_p)
    e_aldif = Ptr[float](e_aldif_p)
    state_tlay = Ptr[float](state_tlay_p)
    state_tlev = Ptr[float](state_tlev_p)
    state_ch4vmr = Ptr[float](state_ch4vmr_p)
    state_o2vmr = Ptr[float](state_o2vmr_p)
    state_n2ovmr = Ptr[float](state_n2ovmr_p)
    pmid = Ptr[float](pmid_p)
    cld = Ptr[float](cld_p)
    pintmb = Ptr[float](pintmb_p)
    pmidmb = Ptr[float](pmidmb_p)
    h2ovmr = Ptr[float](h2ovmr_p)
    o3vmr = Ptr[float](o3vmr_p)
    co2vmr = Ptr[float](co2vmr_p)
    coszrs = Ptr[float](coszrs_p)
    asdir = Ptr[float](asdir_p)
    aldir = Ptr[float](aldir_p)
    asdif = Ptr[float](asdif_p)
    aldif = Ptr[float](aldif_p)
    tlay = Ptr[float](tlay_p)
    tlev = Ptr[float](tlev_p)
    ch4vmr = Ptr[float](ch4vmr_p)
    o2vmr = Ptr[float](o2vmr_p)
    n2ovmr = Ptr[float](n2ovmr_p)

    src_offset = pverp - rrtmg_levs
    _radsw_compact_2d_level_offset(e_pmid, pmid, nday, nnite, rrtmg_levs - 1, src_offset, pcols, idxday, idxnite)
    _radsw_compact_2d_level_offset(e_cld, cld, nday, nnite, rrtmg_levs - 1, src_offset, pcols, idxday, idxnite)
    _radsw_compact_2d_same_levels(state_pintmb, pintmb, nday, nnite, rrtmg_levs + 1, pcols, idxday, idxnite)
    _radsw_compact_2d_same_levels(state_pmidmb, pmidmb, nday, nnite, rrtmg_levs, pcols, idxday, idxnite)
    _radsw_compact_2d_same_levels(state_h2ovmr, h2ovmr, nday, nnite, rrtmg_levs, pcols, idxday, idxnite)
    _radsw_compact_2d_same_levels(state_o3vmr, o3vmr, nday, nnite, rrtmg_levs, pcols, idxday, idxnite)
    _radsw_compact_2d_same_levels(state_co2vmr, co2vmr, nday, nnite, rrtmg_levs, pcols, idxday, idxnite)
    _radsw_compact_1d(e_coszrs, coszrs, nday, nnite, idxday, idxnite)
    _radsw_compact_1d(e_asdir, asdir, nday, nnite, idxday, idxnite)
    _radsw_compact_1d(e_aldir, aldir, nday, nnite, idxday, idxnite)
    _radsw_compact_1d(e_asdif, asdif, nday, nnite, idxday, idxnite)
    _radsw_compact_1d(e_aldif, aldif, nday, nnite, idxday, idxnite)
    _radsw_compact_2d_same_levels(state_tlay, tlay, nday, nnite, rrtmg_levs, pcols, idxday, idxnite)
    _radsw_compact_2d_same_levels(state_tlev, tlev, nday, nnite, rrtmg_levs + 1, pcols, idxday, idxnite)
    _radsw_compact_2d_same_levels(state_ch4vmr, ch4vmr, nday, nnite, rrtmg_levs, pcols, idxday, idxnite)
    _radsw_compact_2d_same_levels(state_o2vmr, o2vmr, nday, nnite, rrtmg_levs, pcols, idxday, idxnite)
    _radsw_compact_2d_same_levels(state_n2ovmr, n2ovmr, nday, nnite, rrtmg_levs, pcols, idxday, idxnite)


@export
def rrtmg_sw_pre_codon(
    nday: int,
    pcols: int,
    pver: int,
    pverp: int,
    rrtmg_levs: int,
    nbndsw: int,
    eccf: float,
    e_aer_tau_p: cobj,
    e_aer_tau_w_p: cobj,
    e_aer_tau_w_g_p: cobj,
    idxday_p: cobj,
    tau_aer_sw_p: cobj,
    ssa_aer_sw_p: cobj,
    asm_aer_sw_p: cobj,
    tlev_p: cobj,
    sfac_p: cobj,
    solar_band_irrad_p: cobj,
    coszrs_p: cobj,
    tsfc_p: cobj,
    solvar_p: cobj,
    solin_p: cobj,
):
    e_aer_tau = Ptr[float](e_aer_tau_p)
    e_aer_tau_w = Ptr[float](e_aer_tau_w_p)
    e_aer_tau_w_g = Ptr[float](e_aer_tau_w_g_p)
    idxday = Ptr[int](idxday_p)
    tau_aer_sw = Ptr[float](tau_aer_sw_p)
    ssa_aer_sw = Ptr[float](ssa_aer_sw_p)
    asm_aer_sw = Ptr[float](asm_aer_sw_p)
    tlev = Ptr[float](tlev_p)
    sfac = Ptr[float](sfac_p)
    solar_band_irrad = Ptr[float](solar_band_irrad_p)
    coszrs = Ptr[float](coszrs_p)
    tsfc = Ptr[float](tsfc_p)
    solvar = Ptr[float](solvar_p)
    solin = Ptr[float](solin_p)

    for ns in range(1, nbndsw + 1):
        for k in range(1, rrtmg_levs):
            kk = (pverp - rrtmg_levs) + k
            for i in range(1, nday + 1):
                col = idxday[i - 1]
                src_tau_w = e_aer_tau_w[_idx3_lb0_dim2(col, kk, ns, pcols, pver)]
                src_tau = e_aer_tau[_idx3_lb0_dim2(col, kk, ns, pcols, pver)]
                dst = _idx3(i, k, ns, pcols, rrtmg_levs - 1)

                if src_tau_w > 1.0e-80:
                    asm_aer_sw[dst] = (
                        e_aer_tau_w_g[_idx3_lb0_dim2(col, kk, ns, pcols, pver)]
                        / src_tau_w
                    )
                else:
                    asm_aer_sw[dst] = 0.0

                if src_tau > 0.0:
                    ssa_aer_sw[dst] = src_tau_w / src_tau
                    tau_aer_sw[dst] = src_tau
                else:
                    ssa_aer_sw[dst] = 1.0
                    tau_aer_sw[dst] = 0.0

    for i in range(1, nday + 1):
        tsfc[i - 1] = tlev[_idx2(i, rrtmg_levs + 1, pcols)]

    for ns in range(1, nbndsw + 1):
        solvar[ns - 1] = sfac[ns - 1]

    for i in range(1, nday + 1):
        band_sum = 0.0
        for ns in range(1, nbndsw + 1):
            band_sum += sfac[ns - 1] * solar_band_irrad[ns - 1]
        solin[i - 1] = band_sum * eccf * coszrs[i - 1]


@export
def rrtmg_sw_cloud_optics_codon(
    nday: int,
    pcols: int,
    pver: int,
    pverp: int,
    rrtmg_levs: int,
    nbndsw: int,
    old_convert: int,
    e_cld_tau_p: cobj,
    e_cld_tau_w_p: cobj,
    e_cld_tau_w_g_p: cobj,
    e_cld_tau_w_f_p: cobj,
    idxday_p: cobj,
    tauc_sw_p: cobj,
    ssac_sw_p: cobj,
    asmc_sw_p: cobj,
    fsfc_sw_p: cobj,
):
    e_cld_tau = Ptr[float](e_cld_tau_p)
    e_cld_tau_w = Ptr[float](e_cld_tau_w_p)
    e_cld_tau_w_g = Ptr[float](e_cld_tau_w_g_p)
    e_cld_tau_w_f = Ptr[float](e_cld_tau_w_f_p)
    idxday = Ptr[int](idxday_p)
    tauc_sw = Ptr[float](tauc_sw_p)
    ssac_sw = Ptr[float](ssac_sw_p)
    asmc_sw = Ptr[float](asmc_sw_p)
    fsfc_sw = Ptr[float](fsfc_sw_p)

    if old_convert == 1:
        for i in range(1, nday + 1):
            col = idxday[i - 1]
            for k in range(1, rrtmg_levs):
                kk = (pverp - rrtmg_levs) + k
                for ns in range(1, nbndsw + 1):
                    src_idx = _idx3(ns, col, kk, nbndsw, pcols)
                    dst_idx = _idx3(ns, i, k, nbndsw, pcols)
                    tau_w = e_cld_tau_w[src_idx]
                    if tau_w > 0.0:
                        fsfc_sw[dst_idx] = e_cld_tau_w_f[src_idx] / tau_w
                        asmc_sw[dst_idx] = e_cld_tau_w_g[src_idx] / tau_w
                    else:
                        fsfc_sw[dst_idx] = 0.0
                        asmc_sw[dst_idx] = 0.0

                    tauc_sw[dst_idx] = e_cld_tau[src_idx]
                    if tauc_sw[dst_idx] > 0.0:
                        ssac_sw[dst_idx] = tau_w / tauc_sw[dst_idx]
                    else:
                        tauc_sw[dst_idx] = 0.0
                        fsfc_sw[dst_idx] = 0.0
                        asmc_sw[dst_idx] = 0.0
                        ssac_sw[dst_idx] = 1.0
    else:
        for i in range(1, nday + 1):
            col = idxday[i - 1]
            for k in range(1, rrtmg_levs):
                kk = (pverp - rrtmg_levs) + k
                for ns in range(1, nbndsw + 1):
                    src_idx = _idx3(ns, col, kk, nbndsw, pcols)
                    dst_idx = _idx3(ns, i, k, nbndsw, pcols)
                    tau_w = e_cld_tau_w[src_idx]
                    if tau_w > 0.0:
                        tau_w_floor = _max_1em80(tau_w)
                        fsfc_sw[dst_idx] = e_cld_tau_w_f[src_idx] / tau_w_floor
                        asmc_sw[dst_idx] = e_cld_tau_w_g[src_idx] / tau_w_floor
                    else:
                        fsfc_sw[dst_idx] = 0.0
                        asmc_sw[dst_idx] = 0.0

                    tauc_sw[dst_idx] = e_cld_tau[src_idx]
                    if tauc_sw[dst_idx] > 0.0:
                        ssac_sw[dst_idx] = _max_1em80(tau_w) / _max_1em80(tauc_sw[dst_idx])
                    else:
                        tauc_sw[dst_idx] = 0.0
                        fsfc_sw[dst_idx] = 0.0
                        asmc_sw[dst_idx] = 0.0
                        ssac_sw[dst_idx] = 1.0


@export
def rrtmg_cloud_ice_optics_sw_codon(
    ncol: int,
    pcols: int,
    pver: int,
    nswbands: int,
    ngd: int,
    iciwpth_p: cobj,
    dei_p: cobj,
    gd_p: cobj,
    ext_p: cobj,
    ssa_p: cobj,
    asm_p: cobj,
    tau_p: cobj,
    tau_w_p: cobj,
    tau_w_g_p: cobj,
    tau_w_f_p: cobj,
):
    iciwpth = Ptr[float](iciwpth_p)
    dei = Ptr[float](dei_p)
    gd = Ptr[float](gd_p)
    ext_table = Ptr[float](ext_p)
    ssa_table = Ptr[float](ssa_p)
    asm_table = Ptr[float](asm_p)
    tau = Ptr[float](tau_p)
    tau_w = Ptr[float](tau_w_p)
    tau_w_g = Ptr[float](tau_w_g_p)
    tau_w_f = Ptr[float](tau_w_f_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            cloud_ice = iciwpth[_idx2(i, k, pcols)]
            deff = dei[_idx2(i, k, pcols)]
            if cloud_ice < 1.0e-80 or deff == 0.0:
                for swband in range(1, nswbands + 1):
                    out_idx = _idx3(swband, i, k, nswbands, pcols)
                    tau[out_idx] = 0.0
                    tau_w[out_idx] = 0.0
                    tau_w_g[out_idx] = 0.0
                    tau_w_f[out_idx] = 0.0
            else:
                for swband in range(1, nswbands + 1):
                    ext = _interp1_bndry_table(gd, ext_table, ngd, swband, deff)
                    ssa = _interp1_bndry_table(gd, ssa_table, ngd, swband, deff)
                    asm = _interp1_bndry_table(gd, asm_table, ngd, swband, deff)
                    out_idx = _idx3(swband, i, k, nswbands, pcols)
                    tau[out_idx] = cloud_ice * ext
                    tau_w[out_idx] = tau[out_idx] * ssa
                    tau_w_g[out_idx] = tau_w[out_idx] * asm
                    tau_w_f[out_idx] = tau_w_g[out_idx] * asm


@export
def rrtmg_cloud_ice_optics_lw_codon(
    ncol: int,
    pcols: int,
    pver: int,
    nlwbands: int,
    ngd: int,
    iciwpth_p: cobj,
    dei_p: cobj,
    gd_p: cobj,
    absor_p: cobj,
    abs_od_p: cobj,
):
    iciwpth = Ptr[float](iciwpth_p)
    dei = Ptr[float](dei_p)
    gd = Ptr[float](gd_p)
    absor_table = Ptr[float](absor_p)
    abs_od = Ptr[float](abs_od_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            cloud_ice = iciwpth[_idx2(i, k, pcols)]
            deff = dei[_idx2(i, k, pcols)]
            if cloud_ice < 1.0e-80 or deff == 0.0:
                for lwband in range(1, nlwbands + 1):
                    abs_od[_idx3(lwband, i, k, nlwbands, pcols)] = 0.0
            else:
                for lwband in range(1, nlwbands + 1):
                    absor = _interp1_bndry_table(gd, absor_table, ngd, lwband, deff)
                    abs_od[_idx3(lwband, i, k, nlwbands, pcols)] = cloud_ice * absor


@export
def rrtmg_cloud_liquid_optics_sw_codon(
    ncol: int,
    pcols: int,
    pver: int,
    nswbands: int,
    nmu: int,
    nlambda: int,
    iclwpth_p: cobj,
    lamc_p: cobj,
    pgam_p: cobj,
    g_mu_p: cobj,
    g_lambda_p: cobj,
    ext_p: cobj,
    ssa_p: cobj,
    asm_p: cobj,
    tau_p: cobj,
    tau_w_p: cobj,
    tau_w_g_p: cobj,
    tau_w_f_p: cobj,
):
    iclwpth = Ptr[float](iclwpth_p)
    lamc = Ptr[float](lamc_p)
    pgam = Ptr[float](pgam_p)
    g_mu = Ptr[float](g_mu_p)
    g_lambda = Ptr[float](g_lambda_p)
    ext_table = Ptr[float](ext_p)
    ssa_table = Ptr[float](ssa_p)
    asm_table = Ptr[float](asm_p)
    tau = Ptr[float](tau_p)
    tau_w = Ptr[float](tau_w_p)
    tau_w_g = Ptr[float](tau_w_g_p)
    tau_w_f = Ptr[float](tau_w_f_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            cell_idx = _idx2(i, k, pcols)
            lambda_cell = lamc[cell_idx]
            cloud_liq = iclwpth[cell_idx]
            if lambda_cell > 0.0 and cloud_liq >= 1.0e-80:
                mu_cell = pgam[cell_idx]
                mu_jjm = _interp1_bndry_jjm(g_mu, nmu, mu_cell)
                mu_jjp = _interp1_bndry_jjp(g_mu, nmu, mu_cell)
                mu_wgts = _interp1_wgts(g_mu, mu_jjm, mu_jjp, mu_cell)
                mu_wgtn = _interp1_wgtn(g_mu, mu_jjm, mu_jjp, mu_cell)
                lam_jjm = _liquid_lambda_jjm(g_lambda, nmu, nlambda, lambda_cell, mu_jjm, mu_jjp, mu_wgts, mu_wgtn)
                lam_jjp = _liquid_lambda_jjp(g_lambda, nmu, nlambda, lambda_cell, mu_jjm, mu_jjp, mu_wgts, mu_wgtn)
                lam_lo = _liquid_lambda_grid_value(g_lambda, nmu, lam_jjm, mu_jjm, mu_jjp, mu_wgts, mu_wgtn)
                lam_hi = _liquid_lambda_grid_value(g_lambda, nmu, lam_jjp, mu_jjm, mu_jjp, mu_wgts, mu_wgtn)
                if lam_jjm == lam_jjp:
                    lam_wgts = 1.0
                    lam_wgtn = 0.0
                else:
                    lam_wgts = (lam_hi - lambda_cell) / (lam_hi - lam_lo)
                    lam_wgtn = (lambda_cell - lam_lo) / (lam_hi - lam_lo)

                for swband in range(1, nswbands + 1):
                    ext = _liquid_interp2(ext_table, nmu, nlambda, swband, mu_jjm, mu_jjp, mu_wgts, mu_wgtn, lam_jjm, lam_jjp, lam_wgts, lam_wgtn)
                    ssa = _liquid_interp2(ssa_table, nmu, nlambda, swband, mu_jjm, mu_jjp, mu_wgts, mu_wgtn, lam_jjm, lam_jjp, lam_wgts, lam_wgtn)
                    asym = _liquid_interp2(asm_table, nmu, nlambda, swband, mu_jjm, mu_jjp, mu_wgts, mu_wgtn, lam_jjm, lam_jjp, lam_wgts, lam_wgtn)
                    out_idx = _idx3(swband, i, k, nswbands, pcols)
                    tau[out_idx] = cloud_liq * ext
                    tau_w[out_idx] = tau[out_idx] * ssa
                    tau_w_g[out_idx] = tau_w[out_idx] * asym
                    tau_w_f[out_idx] = tau_w_g[out_idx] * asym
            else:
                for swband in range(1, nswbands + 1):
                    out_idx = _idx3(swband, i, k, nswbands, pcols)
                    tau[out_idx] = 0.0
                    tau_w[out_idx] = 0.0
                    tau_w_g[out_idx] = 0.0
                    tau_w_f[out_idx] = 0.0


@export
def rrtmg_cloud_liquid_optics_lw_codon(
    ncol: int,
    pcols: int,
    pver: int,
    nlwbands: int,
    nmu: int,
    nlambda: int,
    iclwpth_p: cobj,
    lamc_p: cobj,
    pgam_p: cobj,
    g_mu_p: cobj,
    g_lambda_p: cobj,
    abs_liq_p: cobj,
    abs_od_p: cobj,
):
    iclwpth = Ptr[float](iclwpth_p)
    lamc = Ptr[float](lamc_p)
    pgam = Ptr[float](pgam_p)
    g_mu = Ptr[float](g_mu_p)
    g_lambda = Ptr[float](g_lambda_p)
    abs_liq = Ptr[float](abs_liq_p)
    abs_od = Ptr[float](abs_od_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            cell_idx = _idx2(i, k, pcols)
            lambda_cell = lamc[cell_idx]
            cloud_liq = iclwpth[cell_idx]
            if lambda_cell > 0.0 and cloud_liq >= 1.0e-80:
                mu_cell = pgam[cell_idx]
                mu_jjm = _interp1_bndry_jjm(g_mu, nmu, mu_cell)
                mu_jjp = _interp1_bndry_jjp(g_mu, nmu, mu_cell)
                mu_wgts = _interp1_wgts(g_mu, mu_jjm, mu_jjp, mu_cell)
                mu_wgtn = _interp1_wgtn(g_mu, mu_jjm, mu_jjp, mu_cell)
                lam_jjm = _liquid_lambda_jjm(g_lambda, nmu, nlambda, lambda_cell, mu_jjm, mu_jjp, mu_wgts, mu_wgtn)
                lam_jjp = _liquid_lambda_jjp(g_lambda, nmu, nlambda, lambda_cell, mu_jjm, mu_jjp, mu_wgts, mu_wgtn)
                lam_lo = _liquid_lambda_grid_value(g_lambda, nmu, lam_jjm, mu_jjm, mu_jjp, mu_wgts, mu_wgtn)
                lam_hi = _liquid_lambda_grid_value(g_lambda, nmu, lam_jjp, mu_jjm, mu_jjp, mu_wgts, mu_wgtn)
                if lam_jjm == lam_jjp:
                    lam_wgts = 1.0
                    lam_wgtn = 0.0
                else:
                    lam_wgts = (lam_hi - lambda_cell) / (lam_hi - lam_lo)
                    lam_wgtn = (lambda_cell - lam_lo) / (lam_hi - lam_lo)

                for lwband in range(1, nlwbands + 1):
                    absor = _liquid_interp2(abs_liq, nmu, nlambda, lwband, mu_jjm, mu_jjp, mu_wgts, mu_wgtn, lam_jjm, lam_jjp, lam_wgts, lam_wgtn)
                    abs_od[_idx3(lwband, i, k, nlwbands, pcols)] = cloud_liq * absor
            else:
                for lwband in range(1, nlwbands + 1):
                    abs_od[_idx3(lwband, i, k, nlwbands, pcols)] = 0.0


@export
def rrtmg_sw_post_codon(
    nday: int,
    pcols: int,
    pver: int,
    pverp: int,
    rrtmg_levs: int,
    cpair: float,
    swuflx_p: cobj,
    swdflx_p: cobj,
    swhr_p: cobj,
    swuflxc_p: cobj,
    swdflxc_p: cobj,
    swhrc_p: cobj,
    dirdnuv_p: cobj,
    dirdnir_p: cobj,
    difdnuv_p: cobj,
    difdnir_p: cobj,
    ninflx_p: cobj,
    ninflxc_p: cobj,
    fsntoa_p: cobj,
    fsutoa_p: cobj,
    fsntoac_p: cobj,
    fsnirtoa_p: cobj,
    fsnrtoaq_p: cobj,
    fsnrtoac_p: cobj,
    fsnt_p: cobj,
    fsntc_p: cobj,
    fsds_p: cobj,
    fsdsc_p: cobj,
    fsns_p: cobj,
    fsnsc_p: cobj,
    sols_p: cobj,
    soll_p: cobj,
    solsd_p: cobj,
    solld_p: cobj,
    fns_p: cobj,
    fcns_p: cobj,
    fus_p: cobj,
    fds_p: cobj,
    fusc_p: cobj,
    fdsc_p: cobj,
    qrs_p: cobj,
    qrsc_p: cobj,
):
    swuflx = Ptr[float](swuflx_p)
    swdflx = Ptr[float](swdflx_p)
    swhr = Ptr[float](swhr_p)
    swuflxc = Ptr[float](swuflxc_p)
    swdflxc = Ptr[float](swdflxc_p)
    swhrc = Ptr[float](swhrc_p)
    dirdnuv = Ptr[float](dirdnuv_p)
    dirdnir = Ptr[float](dirdnir_p)
    difdnuv = Ptr[float](difdnuv_p)
    difdnir = Ptr[float](difdnir_p)
    ninflx = Ptr[float](ninflx_p)
    ninflxc = Ptr[float](ninflxc_p)
    fsntoa = Ptr[float](fsntoa_p)
    fsutoa = Ptr[float](fsutoa_p)
    fsntoac = Ptr[float](fsntoac_p)
    fsnirtoa = Ptr[float](fsnirtoa_p)
    fsnrtoaq = Ptr[float](fsnrtoaq_p)
    fsnrtoac = Ptr[float](fsnrtoac_p)
    fsnt = Ptr[float](fsnt_p)
    fsntc = Ptr[float](fsntc_p)
    fsds = Ptr[float](fsds_p)
    fsdsc = Ptr[float](fsdsc_p)
    fsns = Ptr[float](fsns_p)
    fsnsc = Ptr[float](fsnsc_p)
    sols = Ptr[float](sols_p)
    soll = Ptr[float](soll_p)
    solsd = Ptr[float](solsd_p)
    solld = Ptr[float](solld_p)
    fns = Ptr[float](fns_p)
    fcns = Ptr[float](fcns_p)
    fus = Ptr[float](fus_p)
    fds = Ptr[float](fds_p)
    fusc = Ptr[float](fusc_p)
    fdsc = Ptr[float](fdsc_p)
    qrs = Ptr[float](qrs_p)
    qrsc = Ptr[float](qrsc_p)

    dps = 1.0 / 86400.0

    for i in range(1, nday + 1):
        fsntoa[i - 1] = (
            swdflx[_idx2(i, rrtmg_levs + 1, pcols)]
            - swuflx[_idx2(i, rrtmg_levs + 1, pcols)]
        )
        fsutoa[i - 1] = swuflx[_idx2(i, rrtmg_levs + 1, pcols)]
        fsntoac[i - 1] = (
            swdflxc[_idx2(i, rrtmg_levs + 1, pcols)]
            - swuflxc[_idx2(i, rrtmg_levs + 1, pcols)]
        )

        fsnirtoa[i - 1] = ninflx[_idx2(i, rrtmg_levs, pcols)]
        fsnrtoaq[i - 1] = ninflx[_idx2(i, rrtmg_levs, pcols)]
        fsnrtoac[i - 1] = ninflxc[_idx2(i, rrtmg_levs, pcols)]

        fsnt[i - 1] = (
            swdflx[_idx2(i, rrtmg_levs, pcols)]
            - swuflx[_idx2(i, rrtmg_levs, pcols)]
        )
        fsntc[i - 1] = (
            swdflxc[_idx2(i, rrtmg_levs, pcols)]
            - swuflxc[_idx2(i, rrtmg_levs, pcols)]
        )

        fsds[i - 1] = swdflx[_idx2(i, 1, pcols)]
        fsdsc[i - 1] = swdflxc[_idx2(i, 1, pcols)]
        fsns[i - 1] = swdflx[_idx2(i, 1, pcols)] - swuflx[_idx2(i, 1, pcols)]
        fsnsc[i - 1] = swdflxc[_idx2(i, 1, pcols)] - swuflxc[_idx2(i, 1, pcols)]

        sols[i - 1] = dirdnuv[_idx2(i, 1, pcols)]
        soll[i - 1] = dirdnir[_idx2(i, 1, pcols)]
        solsd[i - 1] = difdnuv[_idx2(i, 1, pcols)]
        solld[i - 1] = difdnir[_idx2(i, 1, pcols)]

    for j in range(1, rrtmg_levs + 1):
        cam_k = pverp - rrtmg_levs + j
        rrtmg_k = rrtmg_levs - j + 1
        for i in range(1, nday + 1):
            fns[_idx2(i, cam_k, pcols)] = (
                swdflx[_idx2(i, rrtmg_k, pcols)]
                - swuflx[_idx2(i, rrtmg_k, pcols)]
            )
            fcns[_idx2(i, cam_k, pcols)] = (
                swdflxc[_idx2(i, rrtmg_k, pcols)]
                - swuflxc[_idx2(i, rrtmg_k, pcols)]
            )
            fus[_idx2(i, cam_k, pcols)] = swuflx[_idx2(i, rrtmg_k, pcols)]
            fusc[_idx2(i, cam_k, pcols)] = swuflxc[_idx2(i, rrtmg_k, pcols)]
            fds[_idx2(i, cam_k, pcols)] = swdflx[_idx2(i, rrtmg_k, pcols)]
            fdsc[_idx2(i, cam_k, pcols)] = swdflxc[_idx2(i, rrtmg_k, pcols)]

    for j in range(1, rrtmg_levs):
        cam_k = pverp - rrtmg_levs + j
        rrtmg_k = rrtmg_levs - j
        for i in range(1, nday + 1):
            qrs[_idx2(i, cam_k, pcols)] = (
                swhr[_idx2(i, rrtmg_k, pcols)] * cpair * dps
            )
            qrsc[_idx2(i, cam_k, pcols)] = (
                swhrc[_idx2(i, rrtmg_k, pcols)] * cpair * dps
            )


@inline
def _radsw_expand_1d_zero_night(
    arr: Ptr[float],
    nday: int,
    nnite: int,
    idxday: Ptr[int],
    idxnite: Ptr[int],
):
    # IdxDay/IdxNite are generated in ascending column order; scatter daylight
    # compact values backwards so no Fortran-owned scratch array is needed.
    for pos in range(nday, 0, -1):
        arr[idxday[pos - 1] - 1] = arr[pos - 1]
    for pos in range(1, nnite + 1):
        arr[idxnite[pos - 1] - 1] = 0.0


@inline
def _radsw_expand_2d_zero_night(
    arr: Ptr[float],
    nday: int,
    nnite: int,
    nlev: int,
    pcols: int,
    idxday: Ptr[int],
    idxnite: Ptr[int],
):
    for k in range(1, nlev + 1):
        for pos in range(nday, 0, -1):
            arr[_idx2(idxday[pos - 1], k, pcols)] = arr[_idx2(pos, k, pcols)]
        for pos in range(1, nnite + 1):
            arr[_idx2(idxnite[pos - 1], k, pcols)] = 0.0


@export
def rrtmg_sw_expand_outputs_codon(
    nday: int,
    nnite: int,
    pcols: int,
    pver: int,
    pverp: int,
    idxday_p: cobj,
    idxnite_p: cobj,
    solin_p: cobj,
    qrs_p: cobj,
    qrsc_p: cobj,
    fns_p: cobj,
    fcns_p: cobj,
    fsns_p: cobj,
    fsnt_p: cobj,
    fsntoa_p: cobj,
    fsutoa_p: cobj,
    fsds_p: cobj,
    fsnsc_p: cobj,
    fsdsc_p: cobj,
    fsntc_p: cobj,
    fsntoac_p: cobj,
    sols_p: cobj,
    soll_p: cobj,
    solsd_p: cobj,
    solld_p: cobj,
    fsnirtoa_p: cobj,
    fsnrtoac_p: cobj,
    fsnrtoaq_p: cobj,
):
    idxday = Ptr[int](idxday_p)
    idxnite = Ptr[int](idxnite_p)
    solin = Ptr[float](solin_p)
    qrs = Ptr[float](qrs_p)
    qrsc = Ptr[float](qrsc_p)
    fns = Ptr[float](fns_p)
    fcns = Ptr[float](fcns_p)
    fsns = Ptr[float](fsns_p)
    fsnt = Ptr[float](fsnt_p)
    fsntoa = Ptr[float](fsntoa_p)
    fsutoa = Ptr[float](fsutoa_p)
    fsds = Ptr[float](fsds_p)
    fsnsc = Ptr[float](fsnsc_p)
    fsdsc = Ptr[float](fsdsc_p)
    fsntc = Ptr[float](fsntc_p)
    fsntoac = Ptr[float](fsntoac_p)
    sols = Ptr[float](sols_p)
    soll = Ptr[float](soll_p)
    solsd = Ptr[float](solsd_p)
    solld = Ptr[float](solld_p)
    fsnirtoa = Ptr[float](fsnirtoa_p)
    fsnrtoac = Ptr[float](fsnrtoac_p)
    fsnrtoaq = Ptr[float](fsnrtoaq_p)

    _radsw_expand_1d_zero_night(solin, nday, nnite, idxday, idxnite)
    _radsw_expand_2d_zero_night(qrs, nday, nnite, pver, pcols, idxday, idxnite)
    _radsw_expand_2d_zero_night(qrsc, nday, nnite, pver, pcols, idxday, idxnite)
    _radsw_expand_2d_zero_night(fns, nday, nnite, pverp, pcols, idxday, idxnite)
    _radsw_expand_2d_zero_night(fcns, nday, nnite, pverp, pcols, idxday, idxnite)
    _radsw_expand_1d_zero_night(fsns, nday, nnite, idxday, idxnite)
    _radsw_expand_1d_zero_night(fsnt, nday, nnite, idxday, idxnite)
    _radsw_expand_1d_zero_night(fsntoa, nday, nnite, idxday, idxnite)
    _radsw_expand_1d_zero_night(fsutoa, nday, nnite, idxday, idxnite)
    _radsw_expand_1d_zero_night(fsds, nday, nnite, idxday, idxnite)
    _radsw_expand_1d_zero_night(fsnsc, nday, nnite, idxday, idxnite)
    _radsw_expand_1d_zero_night(fsdsc, nday, nnite, idxday, idxnite)
    _radsw_expand_1d_zero_night(fsntc, nday, nnite, idxday, idxnite)
    _radsw_expand_1d_zero_night(fsntoac, nday, nnite, idxday, idxnite)
    _radsw_expand_1d_zero_night(sols, nday, nnite, idxday, idxnite)
    _radsw_expand_1d_zero_night(soll, nday, nnite, idxday, idxnite)
    _radsw_expand_1d_zero_night(solsd, nday, nnite, idxday, idxnite)
    _radsw_expand_1d_zero_night(solld, nday, nnite, idxday, idxnite)
    _radsw_expand_1d_zero_night(fsnirtoa, nday, nnite, idxday, idxnite)
    _radsw_expand_1d_zero_night(fsnrtoac, nday, nnite, idxday, idxnite)
    _radsw_expand_1d_zero_night(fsnrtoaq, nday, nnite, idxday, idxnite)


@export
def rrtmg_solar_variability_codon(
    nsrc: int,
    ntrg: int,
    src_x_p: cobj,
    min_trg_p: cobj,
    max_trg_p: cobj,
    src_p: cobj,
    trg_p: cobj,
    ref_irrad_p: cobj,
    sfac_p: cobj,
):
    src_x = Ptr[float](src_x_p)
    min_trg = Ptr[float](min_trg_p)
    max_trg = Ptr[float](max_trg_p)
    src = Ptr[float](src_p)
    trg = Ptr[float](trg_p)
    ref_irrad = Ptr[float](ref_irrad_p)
    sfac = Ptr[float](sfac_p)

    for i in range(1, ntrg + 1):
        tl = min_trg[i - 1]
        tu = max_trg[i - 1]
        if tl < src_x[nsrc]:
            sil = 1
            for l in range(1, nsrc + 2):
                if tl <= src_x[l - 1]:
                    sil = l
                    break

            siu = 1
            for l in range(1, nsrc + 2):
                if tu <= src_x[l - 1]:
                    siu = l
                    break

            y = 0.0
            if sil < 2:
                sil = 2
            if siu > nsrc + 1:
                siu = nsrc + 1

            for si in range(sil, siu + 1):
                si1 = si - 1
                src_l = src_x[si1 - 1]
                if tl > src_l:
                    sl = tl
                else:
                    sl = src_l
                src_u = src_x[si - 1]
                if tu < src_u:
                    su = tu
                else:
                    su = src_u
                y = y + (su - sl) * src[si1 - 1]

            targ = y / (tu - tl)
        else:
            targ = 0.0

        trg[i - 1] = targ * (tu - tl)

    for i in range(1, ntrg + 1):
        sfac[i - 1] = trg[i - 1] / ref_irrad[i - 1]


@export
def rrtmg_sw_setcoef_codon(
    nlayers: int,
    mxmol: int,
    pavel_p: cobj,
    tavel_p: cobj,
    coldry_p: cobj,
    wkl_p: cobj,
    laytrop_p: cobj,
    layswtch_p: cobj,
    laylow_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    co2mult_p: cobj,
    colch4_p: cobj,
    colco2_p: cobj,
    colh2o_p: cobj,
    colmol_p: cobj,
    coln2o_p: cobj,
    colo2_p: cobj,
    colo3_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    selffac_p: cobj,
    selffrac_p: cobj,
    indself_p: cobj,
    forfac_p: cobj,
    forfrac_p: cobj,
    indfor_p: cobj,
    preflog_p: cobj,
    tref_p: cobj,
):
    pavel = Ptr[float](pavel_p)
    tavel = Ptr[float](tavel_p)
    coldry = Ptr[float](coldry_p)
    wkl = Ptr[float](wkl_p)
    laytrop = Ptr[int](laytrop_p)
    layswtch = Ptr[int](layswtch_p)
    laylow = Ptr[int](laylow_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    co2mult = Ptr[float](co2mult_p)
    colch4 = Ptr[float](colch4_p)
    colco2 = Ptr[float](colco2_p)
    colh2o = Ptr[float](colh2o_p)
    colmol = Ptr[float](colmol_p)
    coln2o = Ptr[float](coln2o_p)
    colo2 = Ptr[float](colo2_p)
    colo3 = Ptr[float](colo3_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    selffac = Ptr[float](selffac_p)
    selffrac = Ptr[float](selffrac_p)
    indself = Ptr[int](indself_p)
    forfac = Ptr[float](forfac_p)
    forfrac = Ptr[float](forfrac_p)
    indfor = Ptr[int](indfor_p)
    preflog = Ptr[float](preflog_p)
    tref = Ptr[float](tref_p)

    stpfac = 296.0 / 1013.0

    laytrop[0] = 0
    layswtch[0] = 0
    laylow[0] = 0

    for lay in range(1, nlayers + 1):
        idx = lay - 1
        plog = log(pavel[idx])
        jp[idx] = int(36.0 - 5.0 * (plog + 0.04))
        if jp[idx] < 1:
            jp[idx] = 1
        elif jp[idx] > 58:
            jp[idx] = 58

        jp1 = jp[idx] + 1
        fp = 5.0 * (preflog[jp[idx] - 1] - plog)

        jt[idx] = int(3.0 + (tavel[idx] - tref[jp[idx] - 1]) / 15.0)
        if jt[idx] < 1:
            jt[idx] = 1
        elif jt[idx] > 4:
            jt[idx] = 4
        ft = ((tavel[idx] - tref[jp[idx] - 1]) / 15.0) - float(jt[idx] - 3)

        jt1[idx] = int(3.0 + (tavel[idx] - tref[jp1 - 1]) / 15.0)
        if jt1[idx] < 1:
            jt1[idx] = 1
        elif jt1[idx] > 4:
            jt1[idx] = 4
        ft1 = ((tavel[idx] - tref[jp1 - 1]) / 15.0) - float(jt1[idx] - 3)

        water = wkl[_idx2(1, lay, mxmol)] / coldry[idx]
        scalefac = pavel[idx] * stpfac / tavel[idx]

        if plog <= 4.56:
            forfac[idx] = scalefac / (1.0 + water)
            factor = (tavel[idx] - 188.0) / 36.0
            indfor[idx] = 3
            forfrac[idx] = factor - 1.0

            colh2o[idx] = 1.0e-20 * wkl[_idx2(1, lay, mxmol)]
            colco2[idx] = 1.0e-20 * wkl[_idx2(2, lay, mxmol)]
            colo3[idx] = 1.0e-20 * wkl[_idx2(3, lay, mxmol)]
            coln2o[idx] = 1.0e-20 * wkl[_idx2(4, lay, mxmol)]
            colch4[idx] = 1.0e-20 * wkl[_idx2(6, lay, mxmol)]
            colo2[idx] = 1.0e-20 * wkl[_idx2(7, lay, mxmol)]
            colmol[idx] = 1.0e-20 * coldry[idx] + colh2o[idx]
            if colco2[idx] == 0.0:
                colco2[idx] = 1.0e-32 * coldry[idx]
            if coln2o[idx] == 0.0:
                coln2o[idx] = 1.0e-32 * coldry[idx]
            if colch4[idx] == 0.0:
                colch4[idx] = 1.0e-32 * coldry[idx]
            if colo2[idx] == 0.0:
                colo2[idx] = 1.0e-32 * coldry[idx]
            co2reg = 3.55e-24 * coldry[idx]
            co2mult[idx] = (
                (colco2[idx] - co2reg)
                * 272.63
                * exp(-1919.4 / tavel[idx])
                / (8.7604e-4 * tavel[idx])
            )

            selffac[idx] = 0.0
            selffrac[idx] = 0.0
            indself[idx] = 0
        else:
            laytrop[0] = laytrop[0] + 1
            if plog >= 6.62:
                laylow[0] = laylow[0] + 1

            forfac[idx] = scalefac / (1.0 + water)
            factor = (332.0 - tavel[idx]) / 36.0
            ind = int(factor)
            if ind < 1:
                ind = 1
            elif ind > 2:
                ind = 2
            indfor[idx] = ind
            forfrac[idx] = factor - float(indfor[idx])

            selffac[idx] = water * forfac[idx]
            factor = (tavel[idx] - 188.0) / 7.2
            ind = int(factor) - 7
            if ind < 1:
                ind = 1
            elif ind > 9:
                ind = 9
            indself[idx] = ind
            selffrac[idx] = factor - float(indself[idx] + 7)

            colh2o[idx] = 1.0e-20 * wkl[_idx2(1, lay, mxmol)]
            colco2[idx] = 1.0e-20 * wkl[_idx2(2, lay, mxmol)]
            colo3[idx] = 1.0e-20 * wkl[_idx2(3, lay, mxmol)]
            coln2o[idx] = 1.0e-20 * wkl[_idx2(4, lay, mxmol)]
            colch4[idx] = 1.0e-20 * wkl[_idx2(6, lay, mxmol)]
            colo2[idx] = 1.0e-20 * wkl[_idx2(7, lay, mxmol)]
            colmol[idx] = 1.0e-20 * coldry[idx] + colh2o[idx]
            if colco2[idx] == 0.0:
                colco2[idx] = 1.0e-32 * coldry[idx]
            if coln2o[idx] == 0.0:
                coln2o[idx] = 1.0e-32 * coldry[idx]
            if colch4[idx] == 0.0:
                colch4[idx] = 1.0e-32 * coldry[idx]
            if colo2[idx] == 0.0:
                colo2[idx] = 1.0e-32 * coldry[idx]
            co2reg = 3.55e-24 * coldry[idx]
            co2mult[idx] = (
                (colco2[idx] - co2reg)
                * 272.63
                * exp(-1919.4 / tavel[idx])
                / (8.7604e-4 * tavel[idx])
            )

        compfp = 1.0 - fp
        fac10[idx] = compfp * ft
        fac00[idx] = compfp * (1.0 - ft)
        fac11[idx] = fp * ft1
        fac01[idx] = fp * (1.0 - ft1)


@export
def rrtmg_sw_taumol26_codon(
    nlayers: int,
    laytrop: int,
    ng26: int,
    ngs25: int,
    colmol_p: cobj,
    sfluxref_p: cobj,
    rayl_p: cobj,
    sfluxzen_p: cobj,
    taug_p: cobj,
    taur_p: cobj,
):
    colmol = Ptr[float](colmol_p)
    sfluxref = Ptr[float](sfluxref_p)
    rayl = Ptr[float](rayl_p)
    sfluxzen = Ptr[float](sfluxzen_p)
    taug = Ptr[float](taug_p)
    taur = Ptr[float](taur_p)

    laysolfr = laytrop

    for lay in range(1, laytrop + 1):
        for ig in range(1, ng26 + 1):
            if lay == laysolfr:
                sfluxzen[ngs25 + ig - 1] = sfluxref[ig - 1]
            taug[_idx2(lay, ngs25 + ig, nlayers)] = 0.0
            taur[_idx2(lay, ngs25 + ig, nlayers)] = colmol[lay - 1] * rayl[ig - 1]

    for lay in range(laytrop + 1, nlayers + 1):
        for ig in range(1, ng26 + 1):
            taug[_idx2(lay, ngs25 + ig, nlayers)] = 0.0
            taur[_idx2(lay, ngs25 + ig, nlayers)] = colmol[lay - 1] * rayl[ig - 1]


@export
def rrtmg_sw_taumol23_codon(
    nlayers: int,
    laytrop: int,
    ng23: int,
    ngs22: int,
    nspa23: int,
    layreffr: int,
    givfac: float,
    colh2o_p: cobj,
    colmol_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    indself_p: cobj,
    indfor_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    selffac_p: cobj,
    selffrac_p: cobj,
    forfac_p: cobj,
    forfrac_p: cobj,
    absa_p: cobj,
    selfref_p: cobj,
    forref_p: cobj,
    sfluxref_p: cobj,
    rayl_p: cobj,
    sfluxzen_p: cobj,
    taug_p: cobj,
    taur_p: cobj,
):
    colh2o = Ptr[float](colh2o_p)
    colmol = Ptr[float](colmol_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    indself = Ptr[int](indself_p)
    indfor = Ptr[int](indfor_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    selffac = Ptr[float](selffac_p)
    selffrac = Ptr[float](selffrac_p)
    forfac = Ptr[float](forfac_p)
    forfrac = Ptr[float](forfrac_p)
    absa = Ptr[float](absa_p)
    selfref = Ptr[float](selfref_p)
    forref = Ptr[float](forref_p)
    sfluxref = Ptr[float](sfluxref_p)
    rayl = Ptr[float](rayl_p)
    sfluxzen = Ptr[float](sfluxzen_p)
    taug = Ptr[float](taug_p)
    taur = Ptr[float](taur_p)

    laysolfr = laytrop

    for lay in range(1, laytrop + 1):
        if jp[lay - 1] < layreffr and jp[lay] >= layreffr:
            if lay + 1 < laytrop:
                laysolfr = lay + 1
            else:
                laysolfr = laytrop
        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa23 + 1
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa23 + 1
        inds = indself[lay - 1]
        indf = indfor[lay - 1]

        for ig in range(1, ng23 + 1):
            tauray = colmol[lay - 1] * rayl[ig - 1]
            taug[_idx2(lay, ngs22 + ig, nlayers)] = colh2o[lay - 1] * (
                givfac
                * (
                    fac00[lay - 1] * absa[_idx2(ind0, ig, 65)]
                    + fac10[lay - 1] * absa[_idx2(ind0 + 1, ig, 65)]
                    + fac01[lay - 1] * absa[_idx2(ind1, ig, 65)]
                    + fac11[lay - 1] * absa[_idx2(ind1 + 1, ig, 65)]
                )
                + selffac[lay - 1]
                * (
                    selfref[_idx2(inds, ig, 10)]
                    + selffrac[lay - 1]
                    * (
                        selfref[_idx2(inds + 1, ig, 10)]
                        - selfref[_idx2(inds, ig, 10)]
                    )
                )
                + forfac[lay - 1]
                * (
                    forref[_idx2(indf, ig, 3)]
                    + forfrac[lay - 1]
                    * (
                        forref[_idx2(indf + 1, ig, 3)]
                        - forref[_idx2(indf, ig, 3)]
                    )
                )
            )
            if lay == laysolfr:
                sfluxzen[ngs22 + ig - 1] = sfluxref[ig - 1]
            taur[_idx2(lay, ngs22 + ig, nlayers)] = tauray

    for lay in range(laytrop + 1, nlayers + 1):
        for ig in range(1, ng23 + 1):
            taug[_idx2(lay, ngs22 + ig, nlayers)] = 0.0
            taur[_idx2(lay, ngs22 + ig, nlayers)] = colmol[lay - 1] * rayl[ig - 1]


@export
def rrtmg_sw_taumol24_codon(
    nlayers: int,
    laytrop: int,
    ng24: int,
    ngs23: int,
    nspa24: int,
    nspb24: int,
    layreffr: int,
    strrat: float,
    oneminus: float,
    colh2o_p: cobj,
    colo2_p: cobj,
    colo3_p: cobj,
    colmol_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    indself_p: cobj,
    indfor_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    selffac_p: cobj,
    selffrac_p: cobj,
    forfac_p: cobj,
    forfrac_p: cobj,
    absa_p: cobj,
    absb_p: cobj,
    selfref_p: cobj,
    forref_p: cobj,
    sfluxref_p: cobj,
    abso3a_p: cobj,
    abso3b_p: cobj,
    rayla_p: cobj,
    raylb_p: cobj,
    sfluxzen_p: cobj,
    taug_p: cobj,
    taur_p: cobj,
):
    colh2o = Ptr[float](colh2o_p)
    colo2 = Ptr[float](colo2_p)
    colo3 = Ptr[float](colo3_p)
    colmol = Ptr[float](colmol_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    indself = Ptr[int](indself_p)
    indfor = Ptr[int](indfor_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    selffac = Ptr[float](selffac_p)
    selffrac = Ptr[float](selffrac_p)
    forfac = Ptr[float](forfac_p)
    forfrac = Ptr[float](forfrac_p)
    absa = Ptr[float](absa_p)
    absb = Ptr[float](absb_p)
    selfref = Ptr[float](selfref_p)
    forref = Ptr[float](forref_p)
    sfluxref = Ptr[float](sfluxref_p)
    abso3a = Ptr[float](abso3a_p)
    abso3b = Ptr[float](abso3b_p)
    rayla = Ptr[float](rayla_p)
    raylb = Ptr[float](raylb_p)
    sfluxzen = Ptr[float](sfluxzen_p)
    taug = Ptr[float](taug_p)
    taur = Ptr[float](taur_p)

    laysolfr = laytrop

    for lay in range(1, laytrop + 1):
        if jp[lay - 1] < layreffr and jp[lay] >= layreffr:
            if lay + 1 < laytrop:
                laysolfr = lay + 1
            else:
                laysolfr = laytrop
        speccomb = colh2o[lay - 1] + strrat * colo2[lay - 1]
        specparm = colh2o[lay - 1] / speccomb
        if specparm >= oneminus:
            specparm = oneminus
        specmult = 8.0 * specparm
        js = 1 + int(specmult)
        fs = specmult - float(int(specmult))
        fac000 = (1.0 - fs) * fac00[lay - 1]
        fac010 = (1.0 - fs) * fac10[lay - 1]
        fac100 = fs * fac00[lay - 1]
        fac110 = fs * fac10[lay - 1]
        fac001 = (1.0 - fs) * fac01[lay - 1]
        fac011 = (1.0 - fs) * fac11[lay - 1]
        fac101 = fs * fac01[lay - 1]
        fac111 = fs * fac11[lay - 1]
        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa24 + js
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa24 + js
        inds = indself[lay - 1]
        indf = indfor[lay - 1]

        for ig in range(1, ng24 + 1):
            tauray = colmol[lay - 1] * (
                rayla[_idx2(ig, js, ng24)]
                + fs * (rayla[_idx2(ig, js + 1, ng24)] - rayla[_idx2(ig, js, ng24)])
            )
            taug[_idx2(lay, ngs23 + ig, nlayers)] = (
                speccomb
                * (
                    fac000 * absa[_idx2(ind0, ig, 585)]
                    + fac100 * absa[_idx2(ind0 + 1, ig, 585)]
                    + fac010 * absa[_idx2(ind0 + 9, ig, 585)]
                    + fac110 * absa[_idx2(ind0 + 10, ig, 585)]
                    + fac001 * absa[_idx2(ind1, ig, 585)]
                    + fac101 * absa[_idx2(ind1 + 1, ig, 585)]
                    + fac011 * absa[_idx2(ind1 + 9, ig, 585)]
                    + fac111 * absa[_idx2(ind1 + 10, ig, 585)]
                )
                + colo3[lay - 1] * abso3a[ig - 1]
                + colh2o[lay - 1]
                * (
                    selffac[lay - 1]
                    * (
                        selfref[_idx2(inds, ig, 10)]
                        + selffrac[lay - 1]
                        * (
                            selfref[_idx2(inds + 1, ig, 10)]
                            - selfref[_idx2(inds, ig, 10)]
                        )
                    )
                    + forfac[lay - 1]
                    * (
                        forref[_idx2(indf, ig, 3)]
                        + forfrac[lay - 1]
                        * (
                            forref[_idx2(indf + 1, ig, 3)]
                            - forref[_idx2(indf, ig, 3)]
                        )
                    )
                )
            )
            if lay == laysolfr:
                sfluxzen[ngs23 + ig - 1] = sfluxref[_idx2(ig, js, ng24)] + fs * (
                    sfluxref[_idx2(ig, js + 1, ng24)] - sfluxref[_idx2(ig, js, ng24)]
                )
            taur[_idx2(lay, ngs23 + ig, nlayers)] = tauray

    for lay in range(laytrop + 1, nlayers + 1):
        ind0 = ((jp[lay - 1] - 13) * 5 + (jt[lay - 1] - 1)) * nspb24 + 1
        ind1 = ((jp[lay - 1] - 12) * 5 + (jt1[lay - 1] - 1)) * nspb24 + 1

        for ig in range(1, ng24 + 1):
            tauray = colmol[lay - 1] * raylb[ig - 1]
            taug[_idx2(lay, ngs23 + ig, nlayers)] = (
                colo2[lay - 1]
                * (
                    fac00[lay - 1] * absb[_idx2(ind0, ig, 235)]
                    + fac10[lay - 1] * absb[_idx2(ind0 + 1, ig, 235)]
                    + fac01[lay - 1] * absb[_idx2(ind1, ig, 235)]
                    + fac11[lay - 1] * absb[_idx2(ind1 + 1, ig, 235)]
                )
                + colo3[lay - 1] * abso3b[ig - 1]
            )
            taur[_idx2(lay, ngs23 + ig, nlayers)] = tauray


@export
def rrtmg_sw_taumol25_codon(
    nlayers: int,
    laytrop: int,
    ng25: int,
    ngs24: int,
    nspa25: int,
    layreffr: int,
    colh2o_p: cobj,
    colo3_p: cobj,
    colmol_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    absa_p: cobj,
    sfluxref_p: cobj,
    abso3a_p: cobj,
    abso3b_p: cobj,
    rayl_p: cobj,
    sfluxzen_p: cobj,
    taug_p: cobj,
    taur_p: cobj,
):
    colh2o = Ptr[float](colh2o_p)
    colo3 = Ptr[float](colo3_p)
    colmol = Ptr[float](colmol_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    absa = Ptr[float](absa_p)
    sfluxref = Ptr[float](sfluxref_p)
    abso3a = Ptr[float](abso3a_p)
    abso3b = Ptr[float](abso3b_p)
    rayl = Ptr[float](rayl_p)
    sfluxzen = Ptr[float](sfluxzen_p)
    taug = Ptr[float](taug_p)
    taur = Ptr[float](taur_p)

    laysolfr = laytrop

    for lay in range(1, laytrop + 1):
        if jp[lay - 1] < layreffr and jp[lay] >= layreffr:
            if lay + 1 < laytrop:
                laysolfr = lay + 1
            else:
                laysolfr = laytrop
        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa25 + 1
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa25 + 1

        for ig in range(1, ng25 + 1):
            tauray = colmol[lay - 1] * rayl[ig - 1]
            taug[_idx2(lay, ngs24 + ig, nlayers)] = colh2o[lay - 1] * (
                fac00[lay - 1] * absa[_idx2(ind0, ig, 65)]
                + fac10[lay - 1] * absa[_idx2(ind0 + 1, ig, 65)]
                + fac01[lay - 1] * absa[_idx2(ind1, ig, 65)]
                + fac11[lay - 1] * absa[_idx2(ind1 + 1, ig, 65)]
            ) + colo3[lay - 1] * abso3a[ig - 1]
            if lay == laysolfr:
                sfluxzen[ngs24 + ig - 1] = sfluxref[ig - 1]
            taur[_idx2(lay, ngs24 + ig, nlayers)] = tauray

    for lay in range(laytrop + 1, nlayers + 1):
        for ig in range(1, ng25 + 1):
            tauray = colmol[lay - 1] * rayl[ig - 1]
            taug[_idx2(lay, ngs24 + ig, nlayers)] = colo3[lay - 1] * abso3b[ig - 1]
            taur[_idx2(lay, ngs24 + ig, nlayers)] = tauray


@export
def rrtmg_sw_taumol27_codon(
    nlayers: int,
    laytrop: int,
    ng27: int,
    ngs26: int,
    nspa27: int,
    nspb27: int,
    layreffr: int,
    scalekur: float,
    colo3_p: cobj,
    colmol_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    absa_p: cobj,
    absb_p: cobj,
    sfluxref_p: cobj,
    rayl_p: cobj,
    sfluxzen_p: cobj,
    taug_p: cobj,
    taur_p: cobj,
):
    colo3 = Ptr[float](colo3_p)
    colmol = Ptr[float](colmol_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    absa = Ptr[float](absa_p)
    absb = Ptr[float](absb_p)
    sfluxref = Ptr[float](sfluxref_p)
    rayl = Ptr[float](rayl_p)
    sfluxzen = Ptr[float](sfluxzen_p)
    taug = Ptr[float](taug_p)
    taur = Ptr[float](taur_p)

    for lay in range(1, laytrop + 1):
        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa27 + 1
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa27 + 1

        for ig in range(1, ng27 + 1):
            tauray = colmol[lay - 1] * rayl[ig - 1]
            taug[_idx2(lay, ngs26 + ig, nlayers)] = colo3[lay - 1] * (
                fac00[lay - 1] * absa[_idx2(ind0, ig, 65)]
                + fac10[lay - 1] * absa[_idx2(ind0 + 1, ig, 65)]
                + fac01[lay - 1] * absa[_idx2(ind1, ig, 65)]
                + fac11[lay - 1] * absa[_idx2(ind1 + 1, ig, 65)]
            )
            taur[_idx2(lay, ngs26 + ig, nlayers)] = tauray

    laysolfr = nlayers

    for lay in range(laytrop + 1, nlayers + 1):
        if jp[lay - 2] < layreffr and jp[lay - 1] >= layreffr:
            laysolfr = lay
        ind0 = ((jp[lay - 1] - 13) * 5 + (jt[lay - 1] - 1)) * nspb27 + 1
        ind1 = ((jp[lay - 1] - 12) * 5 + (jt1[lay - 1] - 1)) * nspb27 + 1

        for ig in range(1, ng27 + 1):
            tauray = colmol[lay - 1] * rayl[ig - 1]
            taug[_idx2(lay, ngs26 + ig, nlayers)] = colo3[lay - 1] * (
                fac00[lay - 1] * absb[_idx2(ind0, ig, 235)]
                + fac10[lay - 1] * absb[_idx2(ind0 + 1, ig, 235)]
                + fac01[lay - 1] * absb[_idx2(ind1, ig, 235)]
                + fac11[lay - 1] * absb[_idx2(ind1 + 1, ig, 235)]
            )
            if lay == laysolfr:
                sfluxzen[ngs26 + ig - 1] = scalekur * sfluxref[ig - 1]
            taur[_idx2(lay, ngs26 + ig, nlayers)] = tauray


@export
def rrtmg_sw_taumol28_codon(
    nlayers: int,
    laytrop: int,
    ng28: int,
    ngs27: int,
    nspa28: int,
    nspb28: int,
    layreffr: int,
    strrat: float,
    rayl: float,
    oneminus: float,
    colo2_p: cobj,
    colo3_p: cobj,
    colmol_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    absa_p: cobj,
    absb_p: cobj,
    sfluxref_p: cobj,
    sfluxzen_p: cobj,
    taug_p: cobj,
    taur_p: cobj,
):
    colo2 = Ptr[float](colo2_p)
    colo3 = Ptr[float](colo3_p)
    colmol = Ptr[float](colmol_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    absa = Ptr[float](absa_p)
    absb = Ptr[float](absb_p)
    sfluxref = Ptr[float](sfluxref_p)
    sfluxzen = Ptr[float](sfluxzen_p)
    taug = Ptr[float](taug_p)
    taur = Ptr[float](taur_p)

    for lay in range(1, laytrop + 1):
        speccomb = colo3[lay - 1] + strrat * colo2[lay - 1]
        specparm = colo3[lay - 1] / speccomb
        if specparm >= oneminus:
            specparm = oneminus
        specmult = 8.0 * specparm
        js = 1 + int(specmult)
        fs = specmult - float(int(specmult))
        fac000 = (1.0 - fs) * fac00[lay - 1]
        fac010 = (1.0 - fs) * fac10[lay - 1]
        fac100 = fs * fac00[lay - 1]
        fac110 = fs * fac10[lay - 1]
        fac001 = (1.0 - fs) * fac01[lay - 1]
        fac011 = (1.0 - fs) * fac11[lay - 1]
        fac101 = fs * fac01[lay - 1]
        fac111 = fs * fac11[lay - 1]
        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa28 + js
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa28 + js
        tauray = colmol[lay - 1] * rayl

        for ig in range(1, ng28 + 1):
            taug[_idx2(lay, ngs27 + ig, nlayers)] = speccomb * (
                fac000 * absa[_idx2(ind0, ig, 585)]
                + fac100 * absa[_idx2(ind0 + 1, ig, 585)]
                + fac010 * absa[_idx2(ind0 + 9, ig, 585)]
                + fac110 * absa[_idx2(ind0 + 10, ig, 585)]
                + fac001 * absa[_idx2(ind1, ig, 585)]
                + fac101 * absa[_idx2(ind1 + 1, ig, 585)]
                + fac011 * absa[_idx2(ind1 + 9, ig, 585)]
                + fac111 * absa[_idx2(ind1 + 10, ig, 585)]
            )
            taur[_idx2(lay, ngs27 + ig, nlayers)] = tauray

    laysolfr = nlayers

    for lay in range(laytrop + 1, nlayers + 1):
        if jp[lay - 2] < layreffr and jp[lay - 1] >= layreffr:
            laysolfr = lay
        speccomb = colo3[lay - 1] + strrat * colo2[lay - 1]
        specparm = colo3[lay - 1] / speccomb
        if specparm >= oneminus:
            specparm = oneminus
        specmult = 4.0 * specparm
        js = 1 + int(specmult)
        fs = specmult - float(int(specmult))
        fac000 = (1.0 - fs) * fac00[lay - 1]
        fac010 = (1.0 - fs) * fac10[lay - 1]
        fac100 = fs * fac00[lay - 1]
        fac110 = fs * fac10[lay - 1]
        fac001 = (1.0 - fs) * fac01[lay - 1]
        fac011 = (1.0 - fs) * fac11[lay - 1]
        fac101 = fs * fac01[lay - 1]
        fac111 = fs * fac11[lay - 1]
        ind0 = ((jp[lay - 1] - 13) * 5 + (jt[lay - 1] - 1)) * nspb28 + js
        ind1 = ((jp[lay - 1] - 12) * 5 + (jt1[lay - 1] - 1)) * nspb28 + js
        tauray = colmol[lay - 1] * rayl

        for ig in range(1, ng28 + 1):
            taug[_idx2(lay, ngs27 + ig, nlayers)] = speccomb * (
                fac000 * absb[_idx2(ind0, ig, 1175)]
                + fac100 * absb[_idx2(ind0 + 1, ig, 1175)]
                + fac010 * absb[_idx2(ind0 + 5, ig, 1175)]
                + fac110 * absb[_idx2(ind0 + 6, ig, 1175)]
                + fac001 * absb[_idx2(ind1, ig, 1175)]
                + fac101 * absb[_idx2(ind1 + 1, ig, 1175)]
                + fac011 * absb[_idx2(ind1 + 5, ig, 1175)]
                + fac111 * absb[_idx2(ind1 + 6, ig, 1175)]
            )
            if lay == laysolfr:
                sfluxzen[ngs27 + ig - 1] = sfluxref[_idx2(ig, js, ng28)] + fs * (
                    sfluxref[_idx2(ig, js + 1, ng28)] - sfluxref[_idx2(ig, js, ng28)]
                )
            taur[_idx2(lay, ngs27 + ig, nlayers)] = tauray


@export
def rrtmg_sw_taumol29_codon(
    nlayers: int,
    laytrop: int,
    ng29: int,
    ngs28: int,
    nspa29: int,
    nspb29: int,
    layreffr: int,
    rayl: float,
    colh2o_p: cobj,
    colco2_p: cobj,
    colmol_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    indself_p: cobj,
    indfor_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    selffac_p: cobj,
    selffrac_p: cobj,
    forfac_p: cobj,
    forfrac_p: cobj,
    absa_p: cobj,
    absb_p: cobj,
    selfref_p: cobj,
    forref_p: cobj,
    sfluxref_p: cobj,
    absh2o_p: cobj,
    absco2_p: cobj,
    sfluxzen_p: cobj,
    taug_p: cobj,
    taur_p: cobj,
):
    colh2o = Ptr[float](colh2o_p)
    colco2 = Ptr[float](colco2_p)
    colmol = Ptr[float](colmol_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    indself = Ptr[int](indself_p)
    indfor = Ptr[int](indfor_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    selffac = Ptr[float](selffac_p)
    selffrac = Ptr[float](selffrac_p)
    forfac = Ptr[float](forfac_p)
    forfrac = Ptr[float](forfrac_p)
    absa = Ptr[float](absa_p)
    absb = Ptr[float](absb_p)
    selfref = Ptr[float](selfref_p)
    forref = Ptr[float](forref_p)
    sfluxref = Ptr[float](sfluxref_p)
    absh2o = Ptr[float](absh2o_p)
    absco2 = Ptr[float](absco2_p)
    sfluxzen = Ptr[float](sfluxzen_p)
    taug = Ptr[float](taug_p)
    taur = Ptr[float](taur_p)

    for lay in range(1, laytrop + 1):
        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa29 + 1
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa29 + 1
        inds = indself[lay - 1]
        indf = indfor[lay - 1]
        tauray = colmol[lay - 1] * rayl

        for ig in range(1, ng29 + 1):
            taug[_idx2(lay, ngs28 + ig, nlayers)] = (
                colh2o[lay - 1]
                * (
                    (
                        fac00[lay - 1] * absa[_idx2(ind0, ig, 65)]
                        + fac10[lay - 1] * absa[_idx2(ind0 + 1, ig, 65)]
                        + fac01[lay - 1] * absa[_idx2(ind1, ig, 65)]
                        + fac11[lay - 1] * absa[_idx2(ind1 + 1, ig, 65)]
                    )
                    + selffac[lay - 1]
                    * (
                        selfref[_idx2(inds, ig, 10)]
                        + selffrac[lay - 1]
                        * (
                            selfref[_idx2(inds + 1, ig, 10)]
                            - selfref[_idx2(inds, ig, 10)]
                        )
                    )
                    + forfac[lay - 1]
                    * (
                        forref[_idx2(indf, ig, 4)]
                        + forfrac[lay - 1]
                        * (
                            forref[_idx2(indf + 1, ig, 4)]
                            - forref[_idx2(indf, ig, 4)]
                        )
                    )
                )
                + colco2[lay - 1] * absco2[ig - 1]
            )
            taur[_idx2(lay, ngs28 + ig, nlayers)] = tauray

    laysolfr = nlayers

    for lay in range(laytrop + 1, nlayers + 1):
        if jp[lay - 2] < layreffr and jp[lay - 1] >= layreffr:
            laysolfr = lay
        ind0 = ((jp[lay - 1] - 13) * 5 + (jt[lay - 1] - 1)) * nspb29 + 1
        ind1 = ((jp[lay - 1] - 12) * 5 + (jt1[lay - 1] - 1)) * nspb29 + 1
        tauray = colmol[lay - 1] * rayl

        for ig in range(1, ng29 + 1):
            taug[_idx2(lay, ngs28 + ig, nlayers)] = colco2[lay - 1] * (
                fac00[lay - 1] * absb[_idx2(ind0, ig, 235)]
                + fac10[lay - 1] * absb[_idx2(ind0 + 1, ig, 235)]
                + fac01[lay - 1] * absb[_idx2(ind1, ig, 235)]
                + fac11[lay - 1] * absb[_idx2(ind1 + 1, ig, 235)]
            ) + colh2o[lay - 1] * absh2o[ig - 1]
            if lay == laysolfr:
                sfluxzen[ngs28 + ig - 1] = sfluxref[ig - 1]
            taur[_idx2(lay, ngs28 + ig, nlayers)] = tauray


@export
def rrtmg_sw_taumol16_codon(
    nlayers: int,
    laytrop: int,
    ng16: int,
    nspa16: int,
    nspb16: int,
    layreffr: int,
    strrat1: float,
    rayl: float,
    oneminus: float,
    colh2o_p: cobj,
    colch4_p: cobj,
    colmol_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    indself_p: cobj,
    indfor_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    selffac_p: cobj,
    selffrac_p: cobj,
    forfac_p: cobj,
    forfrac_p: cobj,
    absa_p: cobj,
    absb_p: cobj,
    selfref_p: cobj,
    forref_p: cobj,
    sfluxref_p: cobj,
    sfluxzen_p: cobj,
    taug_p: cobj,
    taur_p: cobj,
):
    colh2o = Ptr[float](colh2o_p)
    colch4 = Ptr[float](colch4_p)
    colmol = Ptr[float](colmol_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    indself = Ptr[int](indself_p)
    indfor = Ptr[int](indfor_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    selffac = Ptr[float](selffac_p)
    selffrac = Ptr[float](selffrac_p)
    forfac = Ptr[float](forfac_p)
    forfrac = Ptr[float](forfrac_p)
    absa = Ptr[float](absa_p)
    absb = Ptr[float](absb_p)
    selfref = Ptr[float](selfref_p)
    forref = Ptr[float](forref_p)
    sfluxref = Ptr[float](sfluxref_p)
    sfluxzen = Ptr[float](sfluxzen_p)
    taug = Ptr[float](taug_p)
    taur = Ptr[float](taur_p)

    for lay in range(1, laytrop + 1):
        speccomb = colh2o[lay - 1] + strrat1 * colch4[lay - 1]
        specparm = colh2o[lay - 1] / speccomb
        if specparm >= oneminus:
            specparm = oneminus
        specmult = 8.0 * specparm
        js = 1 + int(specmult)
        fs = specmult - float(int(specmult))
        fac000 = (1.0 - fs) * fac00[lay - 1]
        fac010 = (1.0 - fs) * fac10[lay - 1]
        fac100 = fs * fac00[lay - 1]
        fac110 = fs * fac10[lay - 1]
        fac001 = (1.0 - fs) * fac01[lay - 1]
        fac011 = (1.0 - fs) * fac11[lay - 1]
        fac101 = fs * fac01[lay - 1]
        fac111 = fs * fac11[lay - 1]
        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa16 + js
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa16 + js
        inds = indself[lay - 1]
        indf = indfor[lay - 1]
        tauray = colmol[lay - 1] * rayl

        for ig in range(1, ng16 + 1):
            taug[_idx2(lay, ig, nlayers)] = (
                speccomb
                * (
                    fac000 * absa[_idx2(ind0, ig, 585)]
                    + fac100 * absa[_idx2(ind0 + 1, ig, 585)]
                    + fac010 * absa[_idx2(ind0 + 9, ig, 585)]
                    + fac110 * absa[_idx2(ind0 + 10, ig, 585)]
                    + fac001 * absa[_idx2(ind1, ig, 585)]
                    + fac101 * absa[_idx2(ind1 + 1, ig, 585)]
                    + fac011 * absa[_idx2(ind1 + 9, ig, 585)]
                    + fac111 * absa[_idx2(ind1 + 10, ig, 585)]
                )
                + colh2o[lay - 1]
                * (
                    selffac[lay - 1]
                    * (
                        selfref[_idx2(inds, ig, 10)]
                        + selffrac[lay - 1]
                        * (
                            selfref[_idx2(inds + 1, ig, 10)]
                            - selfref[_idx2(inds, ig, 10)]
                        )
                    )
                    + forfac[lay - 1]
                    * (
                        forref[_idx2(indf, ig, 3)]
                        + forfrac[lay - 1]
                        * (
                            forref[_idx2(indf + 1, ig, 3)]
                            - forref[_idx2(indf, ig, 3)]
                        )
                    )
                )
            )
            taur[_idx2(lay, ig, nlayers)] = tauray

    laysolfr = nlayers

    for lay in range(laytrop + 1, nlayers + 1):
        if jp[lay - 2] < layreffr and jp[lay - 1] >= layreffr:
            laysolfr = lay
        ind0 = ((jp[lay - 1] - 13) * 5 + (jt[lay - 1] - 1)) * nspb16 + 1
        ind1 = ((jp[lay - 1] - 12) * 5 + (jt1[lay - 1] - 1)) * nspb16 + 1
        tauray = colmol[lay - 1] * rayl

        for ig in range(1, ng16 + 1):
            taug[_idx2(lay, ig, nlayers)] = colch4[lay - 1] * (
                fac00[lay - 1] * absb[_idx2(ind0, ig, 235)]
                + fac10[lay - 1] * absb[_idx2(ind0 + 1, ig, 235)]
                + fac01[lay - 1] * absb[_idx2(ind1, ig, 235)]
                + fac11[lay - 1] * absb[_idx2(ind1 + 1, ig, 235)]
            )
            if lay == laysolfr:
                sfluxzen[ig - 1] = sfluxref[ig - 1]
            taur[_idx2(lay, ig, nlayers)] = tauray


@export
def rrtmg_sw_taumol17_codon(
    nlayers: int,
    laytrop: int,
    ng17: int,
    ngs16: int,
    nspa17: int,
    nspb17: int,
    layreffr: int,
    strrat: float,
    rayl: float,
    oneminus: float,
    colh2o_p: cobj,
    colco2_p: cobj,
    colmol_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    indself_p: cobj,
    indfor_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    selffac_p: cobj,
    selffrac_p: cobj,
    forfac_p: cobj,
    forfrac_p: cobj,
    absa_p: cobj,
    absb_p: cobj,
    selfref_p: cobj,
    forref_p: cobj,
    sfluxref_p: cobj,
    sfluxzen_p: cobj,
    taug_p: cobj,
    taur_p: cobj,
):
    colh2o = Ptr[float](colh2o_p)
    colco2 = Ptr[float](colco2_p)
    colmol = Ptr[float](colmol_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    indself = Ptr[int](indself_p)
    indfor = Ptr[int](indfor_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    selffac = Ptr[float](selffac_p)
    selffrac = Ptr[float](selffrac_p)
    forfac = Ptr[float](forfac_p)
    forfrac = Ptr[float](forfrac_p)
    absa = Ptr[float](absa_p)
    absb = Ptr[float](absb_p)
    selfref = Ptr[float](selfref_p)
    forref = Ptr[float](forref_p)
    sfluxref = Ptr[float](sfluxref_p)
    sfluxzen = Ptr[float](sfluxzen_p)
    taug = Ptr[float](taug_p)
    taur = Ptr[float](taur_p)

    for lay in range(1, laytrop + 1):
        speccomb = colh2o[lay - 1] + strrat * colco2[lay - 1]
        specparm = colh2o[lay - 1] / speccomb
        if specparm >= oneminus:
            specparm = oneminus
        specmult = 8.0 * specparm
        js = 1 + int(specmult)
        fs = specmult - float(int(specmult))
        fac000 = (1.0 - fs) * fac00[lay - 1]
        fac010 = (1.0 - fs) * fac10[lay - 1]
        fac100 = fs * fac00[lay - 1]
        fac110 = fs * fac10[lay - 1]
        fac001 = (1.0 - fs) * fac01[lay - 1]
        fac011 = (1.0 - fs) * fac11[lay - 1]
        fac101 = fs * fac01[lay - 1]
        fac111 = fs * fac11[lay - 1]
        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa17 + js
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa17 + js
        inds = indself[lay - 1]
        indf = indfor[lay - 1]
        tauray = colmol[lay - 1] * rayl

        for ig in range(1, ng17 + 1):
            taug[_idx2(lay, ngs16 + ig, nlayers)] = (
                speccomb
                * (
                    fac000 * absa[_idx2(ind0, ig, 585)]
                    + fac100 * absa[_idx2(ind0 + 1, ig, 585)]
                    + fac010 * absa[_idx2(ind0 + 9, ig, 585)]
                    + fac110 * absa[_idx2(ind0 + 10, ig, 585)]
                    + fac001 * absa[_idx2(ind1, ig, 585)]
                    + fac101 * absa[_idx2(ind1 + 1, ig, 585)]
                    + fac011 * absa[_idx2(ind1 + 9, ig, 585)]
                    + fac111 * absa[_idx2(ind1 + 10, ig, 585)]
                )
                + colh2o[lay - 1]
                * (
                    selffac[lay - 1]
                    * (
                        selfref[_idx2(inds, ig, 10)]
                        + selffrac[lay - 1]
                        * (
                            selfref[_idx2(inds + 1, ig, 10)]
                            - selfref[_idx2(inds, ig, 10)]
                        )
                    )
                    + forfac[lay - 1]
                    * (
                        forref[_idx2(indf, ig, 4)]
                        + forfrac[lay - 1]
                        * (
                            forref[_idx2(indf + 1, ig, 4)]
                            - forref[_idx2(indf, ig, 4)]
                        )
                    )
                )
            )
            taur[_idx2(lay, ngs16 + ig, nlayers)] = tauray

    laysolfr = nlayers

    for lay in range(laytrop + 1, nlayers + 1):
        if jp[lay - 2] < layreffr and jp[lay - 1] >= layreffr:
            laysolfr = lay
        speccomb = colh2o[lay - 1] + strrat * colco2[lay - 1]
        specparm = colh2o[lay - 1] / speccomb
        if specparm >= oneminus:
            specparm = oneminus
        specmult = 4.0 * specparm
        js = 1 + int(specmult)
        fs = specmult - float(int(specmult))
        fac000 = (1.0 - fs) * fac00[lay - 1]
        fac010 = (1.0 - fs) * fac10[lay - 1]
        fac100 = fs * fac00[lay - 1]
        fac110 = fs * fac10[lay - 1]
        fac001 = (1.0 - fs) * fac01[lay - 1]
        fac011 = (1.0 - fs) * fac11[lay - 1]
        fac101 = fs * fac01[lay - 1]
        fac111 = fs * fac11[lay - 1]
        ind0 = ((jp[lay - 1] - 13) * 5 + (jt[lay - 1] - 1)) * nspb17 + js
        ind1 = ((jp[lay - 1] - 12) * 5 + (jt1[lay - 1] - 1)) * nspb17 + js
        indf = indfor[lay - 1]
        tauray = colmol[lay - 1] * rayl

        for ig in range(1, ng17 + 1):
            taug[_idx2(lay, ngs16 + ig, nlayers)] = (
                speccomb
                * (
                    fac000 * absb[_idx2(ind0, ig, 1175)]
                    + fac100 * absb[_idx2(ind0 + 1, ig, 1175)]
                    + fac010 * absb[_idx2(ind0 + 5, ig, 1175)]
                    + fac110 * absb[_idx2(ind0 + 6, ig, 1175)]
                    + fac001 * absb[_idx2(ind1, ig, 1175)]
                    + fac101 * absb[_idx2(ind1 + 1, ig, 1175)]
                    + fac011 * absb[_idx2(ind1 + 5, ig, 1175)]
                    + fac111 * absb[_idx2(ind1 + 6, ig, 1175)]
                )
                + colh2o[lay - 1]
                * forfac[lay - 1]
                * (
                    forref[_idx2(indf, ig, 4)]
                    + forfrac[lay - 1]
                    * (
                        forref[_idx2(indf + 1, ig, 4)]
                        - forref[_idx2(indf, ig, 4)]
                    )
                )
            )
            if lay == laysolfr:
                sfluxzen[ngs16 + ig - 1] = sfluxref[_idx2(ig, js, ng17)] + fs * (
                    sfluxref[_idx2(ig, js + 1, ng17)] - sfluxref[_idx2(ig, js, ng17)]
                )
            taur[_idx2(lay, ngs16 + ig, nlayers)] = tauray


@export
def rrtmg_sw_taumol18_codon(
    nlayers: int,
    laytrop: int,
    ng18: int,
    ngs17: int,
    nspa18: int,
    nspb18: int,
    layreffr: int,
    strrat: float,
    rayl: float,
    oneminus: float,
    colh2o_p: cobj,
    colch4_p: cobj,
    colmol_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    indself_p: cobj,
    indfor_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    selffac_p: cobj,
    selffrac_p: cobj,
    forfac_p: cobj,
    forfrac_p: cobj,
    absa_p: cobj,
    absb_p: cobj,
    selfref_p: cobj,
    forref_p: cobj,
    sfluxref_p: cobj,
    sfluxzen_p: cobj,
    taug_p: cobj,
    taur_p: cobj,
):
    colh2o = Ptr[float](colh2o_p)
    colch4 = Ptr[float](colch4_p)
    colmol = Ptr[float](colmol_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    indself = Ptr[int](indself_p)
    indfor = Ptr[int](indfor_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    selffac = Ptr[float](selffac_p)
    selffrac = Ptr[float](selffrac_p)
    forfac = Ptr[float](forfac_p)
    forfrac = Ptr[float](forfrac_p)
    absa = Ptr[float](absa_p)
    absb = Ptr[float](absb_p)
    selfref = Ptr[float](selfref_p)
    forref = Ptr[float](forref_p)
    sfluxref = Ptr[float](sfluxref_p)
    sfluxzen = Ptr[float](sfluxzen_p)
    taug = Ptr[float](taug_p)
    taur = Ptr[float](taur_p)

    laysolfr = laytrop

    for lay in range(1, laytrop + 1):
        if jp[lay - 1] < layreffr and jp[lay] >= layreffr:
            if lay + 1 < laytrop:
                laysolfr = lay + 1
            else:
                laysolfr = laytrop
        speccomb = colh2o[lay - 1] + strrat * colch4[lay - 1]
        specparm = colh2o[lay - 1] / speccomb
        if specparm >= oneminus:
            specparm = oneminus
        specmult = 8.0 * specparm
        js = 1 + int(specmult)
        fs = specmult - float(int(specmult))
        fac000 = (1.0 - fs) * fac00[lay - 1]
        fac010 = (1.0 - fs) * fac10[lay - 1]
        fac100 = fs * fac00[lay - 1]
        fac110 = fs * fac10[lay - 1]
        fac001 = (1.0 - fs) * fac01[lay - 1]
        fac011 = (1.0 - fs) * fac11[lay - 1]
        fac101 = fs * fac01[lay - 1]
        fac111 = fs * fac11[lay - 1]
        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa18 + js
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa18 + js
        inds = indself[lay - 1]
        indf = indfor[lay - 1]
        tauray = colmol[lay - 1] * rayl

        for ig in range(1, ng18 + 1):
            taug[_idx2(lay, ngs17 + ig, nlayers)] = (
                speccomb
                * (
                    fac000 * absa[_idx2(ind0, ig, 585)]
                    + fac100 * absa[_idx2(ind0 + 1, ig, 585)]
                    + fac010 * absa[_idx2(ind0 + 9, ig, 585)]
                    + fac110 * absa[_idx2(ind0 + 10, ig, 585)]
                    + fac001 * absa[_idx2(ind1, ig, 585)]
                    + fac101 * absa[_idx2(ind1 + 1, ig, 585)]
                    + fac011 * absa[_idx2(ind1 + 9, ig, 585)]
                    + fac111 * absa[_idx2(ind1 + 10, ig, 585)]
                )
                + colh2o[lay - 1]
                * (
                    selffac[lay - 1]
                    * (
                        selfref[_idx2(inds, ig, 10)]
                        + selffrac[lay - 1]
                        * (
                            selfref[_idx2(inds + 1, ig, 10)]
                            - selfref[_idx2(inds, ig, 10)]
                        )
                    )
                    + forfac[lay - 1]
                    * (
                        forref[_idx2(indf, ig, 3)]
                        + forfrac[lay - 1]
                        * (
                            forref[_idx2(indf + 1, ig, 3)]
                            - forref[_idx2(indf, ig, 3)]
                        )
                    )
                )
            )
            if lay == laysolfr:
                sfluxzen[ngs17 + ig - 1] = sfluxref[_idx2(ig, js, ng18)] + fs * (
                    sfluxref[_idx2(ig, js + 1, ng18)] - sfluxref[_idx2(ig, js, ng18)]
                )
            taur[_idx2(lay, ngs17 + ig, nlayers)] = tauray

    for lay in range(laytrop + 1, nlayers + 1):
        ind0 = ((jp[lay - 1] - 13) * 5 + (jt[lay - 1] - 1)) * nspb18 + 1
        ind1 = ((jp[lay - 1] - 12) * 5 + (jt1[lay - 1] - 1)) * nspb18 + 1
        tauray = colmol[lay - 1] * rayl

        for ig in range(1, ng18 + 1):
            taug[_idx2(lay, ngs17 + ig, nlayers)] = colch4[lay - 1] * (
                fac00[lay - 1] * absb[_idx2(ind0, ig, 235)]
                + fac10[lay - 1] * absb[_idx2(ind0 + 1, ig, 235)]
                + fac01[lay - 1] * absb[_idx2(ind1, ig, 235)]
                + fac11[lay - 1] * absb[_idx2(ind1 + 1, ig, 235)]
            )
            taur[_idx2(lay, ngs17 + ig, nlayers)] = tauray


@export
def rrtmg_sw_taumol19_codon(
    nlayers: int,
    laytrop: int,
    ng19: int,
    ngs18: int,
    nspa19: int,
    nspb19: int,
    layreffr: int,
    strrat: float,
    rayl: float,
    oneminus: float,
    colh2o_p: cobj,
    colco2_p: cobj,
    colmol_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    indself_p: cobj,
    indfor_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    selffac_p: cobj,
    selffrac_p: cobj,
    forfac_p: cobj,
    forfrac_p: cobj,
    absa_p: cobj,
    absb_p: cobj,
    selfref_p: cobj,
    forref_p: cobj,
    sfluxref_p: cobj,
    sfluxzen_p: cobj,
    taug_p: cobj,
    taur_p: cobj,
):
    colh2o = Ptr[float](colh2o_p)
    colco2 = Ptr[float](colco2_p)
    colmol = Ptr[float](colmol_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    indself = Ptr[int](indself_p)
    indfor = Ptr[int](indfor_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    selffac = Ptr[float](selffac_p)
    selffrac = Ptr[float](selffrac_p)
    forfac = Ptr[float](forfac_p)
    forfrac = Ptr[float](forfrac_p)
    absa = Ptr[float](absa_p)
    absb = Ptr[float](absb_p)
    selfref = Ptr[float](selfref_p)
    forref = Ptr[float](forref_p)
    sfluxref = Ptr[float](sfluxref_p)
    sfluxzen = Ptr[float](sfluxzen_p)
    taug = Ptr[float](taug_p)
    taur = Ptr[float](taur_p)

    laysolfr = laytrop

    for lay in range(1, laytrop + 1):
        if jp[lay - 1] < layreffr and jp[lay] >= layreffr:
            if lay + 1 < laytrop:
                laysolfr = lay + 1
            else:
                laysolfr = laytrop
        speccomb = colh2o[lay - 1] + strrat * colco2[lay - 1]
        specparm = colh2o[lay - 1] / speccomb
        if specparm >= oneminus:
            specparm = oneminus
        specmult = 8.0 * specparm
        js = 1 + int(specmult)
        fs = specmult - float(int(specmult))
        fac000 = (1.0 - fs) * fac00[lay - 1]
        fac010 = (1.0 - fs) * fac10[lay - 1]
        fac100 = fs * fac00[lay - 1]
        fac110 = fs * fac10[lay - 1]
        fac001 = (1.0 - fs) * fac01[lay - 1]
        fac011 = (1.0 - fs) * fac11[lay - 1]
        fac101 = fs * fac01[lay - 1]
        fac111 = fs * fac11[lay - 1]
        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa19 + js
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa19 + js
        inds = indself[lay - 1]
        indf = indfor[lay - 1]
        tauray = colmol[lay - 1] * rayl

        for ig in range(1, ng19 + 1):
            taug[_idx2(lay, ngs18 + ig, nlayers)] = (
                speccomb
                * (
                    fac000 * absa[_idx2(ind0, ig, 585)]
                    + fac100 * absa[_idx2(ind0 + 1, ig, 585)]
                    + fac010 * absa[_idx2(ind0 + 9, ig, 585)]
                    + fac110 * absa[_idx2(ind0 + 10, ig, 585)]
                    + fac001 * absa[_idx2(ind1, ig, 585)]
                    + fac101 * absa[_idx2(ind1 + 1, ig, 585)]
                    + fac011 * absa[_idx2(ind1 + 9, ig, 585)]
                    + fac111 * absa[_idx2(ind1 + 10, ig, 585)]
                )
                + colh2o[lay - 1]
                * (
                    selffac[lay - 1]
                    * (
                        selfref[_idx2(inds, ig, 10)]
                        + selffrac[lay - 1]
                        * (
                            selfref[_idx2(inds + 1, ig, 10)]
                            - selfref[_idx2(inds, ig, 10)]
                        )
                    )
                    + forfac[lay - 1]
                    * (
                        forref[_idx2(indf, ig, 3)]
                        + forfrac[lay - 1]
                        * (
                            forref[_idx2(indf + 1, ig, 3)]
                            - forref[_idx2(indf, ig, 3)]
                        )
                    )
                )
            )
            if lay == laysolfr:
                sfluxzen[ngs18 + ig - 1] = sfluxref[_idx2(ig, js, ng19)] + fs * (
                    sfluxref[_idx2(ig, js + 1, ng19)] - sfluxref[_idx2(ig, js, ng19)]
                )
            taur[_idx2(lay, ngs18 + ig, nlayers)] = tauray

    for lay in range(laytrop + 1, nlayers + 1):
        ind0 = ((jp[lay - 1] - 13) * 5 + (jt[lay - 1] - 1)) * nspb19 + 1
        ind1 = ((jp[lay - 1] - 12) * 5 + (jt1[lay - 1] - 1)) * nspb19 + 1
        tauray = colmol[lay - 1] * rayl

        for ig in range(1, ng19 + 1):
            taug[_idx2(lay, ngs18 + ig, nlayers)] = colco2[lay - 1] * (
                fac00[lay - 1] * absb[_idx2(ind0, ig, 235)]
                + fac10[lay - 1] * absb[_idx2(ind0 + 1, ig, 235)]
                + fac01[lay - 1] * absb[_idx2(ind1, ig, 235)]
                + fac11[lay - 1] * absb[_idx2(ind1 + 1, ig, 235)]
            )
            taur[_idx2(lay, ngs18 + ig, nlayers)] = tauray


@export
def rrtmg_sw_taumol20_codon(
    nlayers: int,
    laytrop: int,
    ng20: int,
    ngs19: int,
    nspa20: int,
    nspb20: int,
    layreffr: int,
    rayl: float,
    colh2o_p: cobj,
    colch4_p: cobj,
    colmol_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    indself_p: cobj,
    indfor_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    selffac_p: cobj,
    selffrac_p: cobj,
    forfac_p: cobj,
    forfrac_p: cobj,
    absa_p: cobj,
    absb_p: cobj,
    selfref_p: cobj,
    forref_p: cobj,
    sfluxref_p: cobj,
    absch4_p: cobj,
    sfluxzen_p: cobj,
    taug_p: cobj,
    taur_p: cobj,
):
    colh2o = Ptr[float](colh2o_p)
    colch4 = Ptr[float](colch4_p)
    colmol = Ptr[float](colmol_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    indself = Ptr[int](indself_p)
    indfor = Ptr[int](indfor_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    selffac = Ptr[float](selffac_p)
    selffrac = Ptr[float](selffrac_p)
    forfac = Ptr[float](forfac_p)
    forfrac = Ptr[float](forfrac_p)
    absa = Ptr[float](absa_p)
    absb = Ptr[float](absb_p)
    selfref = Ptr[float](selfref_p)
    forref = Ptr[float](forref_p)
    sfluxref = Ptr[float](sfluxref_p)
    absch4 = Ptr[float](absch4_p)
    sfluxzen = Ptr[float](sfluxzen_p)
    taug = Ptr[float](taug_p)
    taur = Ptr[float](taur_p)

    laysolfr = laytrop

    for lay in range(1, laytrop + 1):
        if jp[lay - 1] < layreffr and jp[lay] >= layreffr:
            if lay + 1 < laytrop:
                laysolfr = lay + 1
            else:
                laysolfr = laytrop
        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa20 + 1
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa20 + 1
        inds = indself[lay - 1]
        indf = indfor[lay - 1]
        tauray = colmol[lay - 1] * rayl

        for ig in range(1, ng20 + 1):
            taug[_idx2(lay, ngs19 + ig, nlayers)] = colh2o[lay - 1] * (
                (
                    fac00[lay - 1] * absa[_idx2(ind0, ig, 65)]
                    + fac10[lay - 1] * absa[_idx2(ind0 + 1, ig, 65)]
                    + fac01[lay - 1] * absa[_idx2(ind1, ig, 65)]
                    + fac11[lay - 1] * absa[_idx2(ind1 + 1, ig, 65)]
                )
                + selffac[lay - 1]
                * (
                    selfref[_idx2(inds, ig, 10)]
                    + selffrac[lay - 1]
                    * (
                        selfref[_idx2(inds + 1, ig, 10)]
                        - selfref[_idx2(inds, ig, 10)]
                    )
                )
                + forfac[lay - 1]
                * (
                    forref[_idx2(indf, ig, 4)]
                    + forfrac[lay - 1]
                    * (
                        forref[_idx2(indf + 1, ig, 4)]
                        - forref[_idx2(indf, ig, 4)]
                    )
                )
            ) + colch4[lay - 1] * absch4[ig - 1]
            taur[_idx2(lay, ngs19 + ig, nlayers)] = tauray
            if lay == laysolfr:
                sfluxzen[ngs19 + ig - 1] = sfluxref[ig - 1]

    for lay in range(laytrop + 1, nlayers + 1):
        ind0 = ((jp[lay - 1] - 13) * 5 + (jt[lay - 1] - 1)) * nspb20 + 1
        ind1 = ((jp[lay - 1] - 12) * 5 + (jt1[lay - 1] - 1)) * nspb20 + 1
        indf = indfor[lay - 1]
        tauray = colmol[lay - 1] * rayl

        for ig in range(1, ng20 + 1):
            taug[_idx2(lay, ngs19 + ig, nlayers)] = colh2o[lay - 1] * (
                fac00[lay - 1] * absb[_idx2(ind0, ig, 235)]
                + fac10[lay - 1] * absb[_idx2(ind0 + 1, ig, 235)]
                + fac01[lay - 1] * absb[_idx2(ind1, ig, 235)]
                + fac11[lay - 1] * absb[_idx2(ind1 + 1, ig, 235)]
                + forfac[lay - 1]
                * (
                    forref[_idx2(indf, ig, 4)]
                    + forfrac[lay - 1]
                    * (
                        forref[_idx2(indf + 1, ig, 4)]
                        - forref[_idx2(indf, ig, 4)]
                    )
                )
            ) + colch4[lay - 1] * absch4[ig - 1]
            taur[_idx2(lay, ngs19 + ig, nlayers)] = tauray


@export
def rrtmg_sw_taumol21_codon(
    nlayers: int,
    laytrop: int,
    ng21: int,
    ngs20: int,
    nspa21: int,
    nspb21: int,
    layreffr: int,
    strrat: float,
    rayl: float,
    oneminus: float,
    colh2o_p: cobj,
    colco2_p: cobj,
    colmol_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    indself_p: cobj,
    indfor_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    selffac_p: cobj,
    selffrac_p: cobj,
    forfac_p: cobj,
    forfrac_p: cobj,
    absa_p: cobj,
    absb_p: cobj,
    selfref_p: cobj,
    forref_p: cobj,
    sfluxref_p: cobj,
    sfluxzen_p: cobj,
    taug_p: cobj,
    taur_p: cobj,
):
    colh2o = Ptr[float](colh2o_p)
    colco2 = Ptr[float](colco2_p)
    colmol = Ptr[float](colmol_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    indself = Ptr[int](indself_p)
    indfor = Ptr[int](indfor_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    selffac = Ptr[float](selffac_p)
    selffrac = Ptr[float](selffrac_p)
    forfac = Ptr[float](forfac_p)
    forfrac = Ptr[float](forfrac_p)
    absa = Ptr[float](absa_p)
    absb = Ptr[float](absb_p)
    selfref = Ptr[float](selfref_p)
    forref = Ptr[float](forref_p)
    sfluxref = Ptr[float](sfluxref_p)
    sfluxzen = Ptr[float](sfluxzen_p)
    taug = Ptr[float](taug_p)
    taur = Ptr[float](taur_p)

    laysolfr = laytrop

    for lay in range(1, laytrop + 1):
        if jp[lay - 1] < layreffr and jp[lay] >= layreffr:
            if lay + 1 < laytrop:
                laysolfr = lay + 1
            else:
                laysolfr = laytrop
        speccomb = colh2o[lay - 1] + strrat * colco2[lay - 1]
        specparm = colh2o[lay - 1] / speccomb
        if specparm >= oneminus:
            specparm = oneminus
        specmult = 8.0 * specparm
        js = 1 + int(specmult)
        fs = specmult - float(int(specmult))
        fac000 = (1.0 - fs) * fac00[lay - 1]
        fac010 = (1.0 - fs) * fac10[lay - 1]
        fac100 = fs * fac00[lay - 1]
        fac110 = fs * fac10[lay - 1]
        fac001 = (1.0 - fs) * fac01[lay - 1]
        fac011 = (1.0 - fs) * fac11[lay - 1]
        fac101 = fs * fac01[lay - 1]
        fac111 = fs * fac11[lay - 1]
        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa21 + js
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa21 + js
        inds = indself[lay - 1]
        indf = indfor[lay - 1]
        tauray = colmol[lay - 1] * rayl

        for ig in range(1, ng21 + 1):
            taug[_idx2(lay, ngs20 + ig, nlayers)] = (
                speccomb
                * (
                    fac000 * absa[_idx2(ind0, ig, 585)]
                    + fac100 * absa[_idx2(ind0 + 1, ig, 585)]
                    + fac010 * absa[_idx2(ind0 + 9, ig, 585)]
                    + fac110 * absa[_idx2(ind0 + 10, ig, 585)]
                    + fac001 * absa[_idx2(ind1, ig, 585)]
                    + fac101 * absa[_idx2(ind1 + 1, ig, 585)]
                    + fac011 * absa[_idx2(ind1 + 9, ig, 585)]
                    + fac111 * absa[_idx2(ind1 + 10, ig, 585)]
                )
                + colh2o[lay - 1]
                * (
                    selffac[lay - 1]
                    * (
                        selfref[_idx2(inds, ig, 10)]
                        + selffrac[lay - 1]
                        * (
                            selfref[_idx2(inds + 1, ig, 10)]
                            - selfref[_idx2(inds, ig, 10)]
                        )
                    )
                    + forfac[lay - 1]
                    * (
                        forref[_idx2(indf, ig, 4)]
                        + forfrac[lay - 1]
                        * (
                            forref[_idx2(indf + 1, ig, 4)]
                            - forref[_idx2(indf, ig, 4)]
                        )
                    )
                )
            )
            if lay == laysolfr:
                sfluxzen[ngs20 + ig - 1] = sfluxref[_idx2(ig, js, ng21)] + fs * (
                    sfluxref[_idx2(ig, js + 1, ng21)] - sfluxref[_idx2(ig, js, ng21)]
                )
            taur[_idx2(lay, ngs20 + ig, nlayers)] = tauray

    for lay in range(laytrop + 1, nlayers + 1):
        speccomb = colh2o[lay - 1] + strrat * colco2[lay - 1]
        specparm = colh2o[lay - 1] / speccomb
        if specparm >= oneminus:
            specparm = oneminus
        specmult = 4.0 * specparm
        js = 1 + int(specmult)
        fs = specmult - float(int(specmult))
        fac000 = (1.0 - fs) * fac00[lay - 1]
        fac010 = (1.0 - fs) * fac10[lay - 1]
        fac100 = fs * fac00[lay - 1]
        fac110 = fs * fac10[lay - 1]
        fac001 = (1.0 - fs) * fac01[lay - 1]
        fac011 = (1.0 - fs) * fac11[lay - 1]
        fac101 = fs * fac01[lay - 1]
        fac111 = fs * fac11[lay - 1]
        ind0 = ((jp[lay - 1] - 13) * 5 + (jt[lay - 1] - 1)) * nspb21 + js
        ind1 = ((jp[lay - 1] - 12) * 5 + (jt1[lay - 1] - 1)) * nspb21 + js
        indf = indfor[lay - 1]
        tauray = colmol[lay - 1] * rayl

        for ig in range(1, ng21 + 1):
            taug[_idx2(lay, ngs20 + ig, nlayers)] = (
                speccomb
                * (
                    fac000 * absb[_idx2(ind0, ig, 1175)]
                    + fac100 * absb[_idx2(ind0 + 1, ig, 1175)]
                    + fac010 * absb[_idx2(ind0 + 5, ig, 1175)]
                    + fac110 * absb[_idx2(ind0 + 6, ig, 1175)]
                    + fac001 * absb[_idx2(ind1, ig, 1175)]
                    + fac101 * absb[_idx2(ind1 + 1, ig, 1175)]
                    + fac011 * absb[_idx2(ind1 + 5, ig, 1175)]
                    + fac111 * absb[_idx2(ind1 + 6, ig, 1175)]
                )
                + colh2o[lay - 1]
                * forfac[lay - 1]
                * (
                    forref[_idx2(indf, ig, 4)]
                    + forfrac[lay - 1]
                    * (
                        forref[_idx2(indf + 1, ig, 4)]
                        - forref[_idx2(indf, ig, 4)]
                    )
                )
            )
            taur[_idx2(lay, ngs20 + ig, nlayers)] = tauray


@export
def rrtmg_sw_taumol22_codon(
    nlayers: int,
    laytrop: int,
    ng22: int,
    ngs21: int,
    nspa22: int,
    nspb22: int,
    layreffr: int,
    strrat: float,
    rayl: float,
    oneminus: float,
    colh2o_p: cobj,
    colo2_p: cobj,
    colmol_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    indself_p: cobj,
    indfor_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    selffac_p: cobj,
    selffrac_p: cobj,
    forfac_p: cobj,
    forfrac_p: cobj,
    absa_p: cobj,
    absb_p: cobj,
    selfref_p: cobj,
    forref_p: cobj,
    sfluxref_p: cobj,
    sfluxzen_p: cobj,
    taug_p: cobj,
    taur_p: cobj,
):
    colh2o = Ptr[float](colh2o_p)
    colo2 = Ptr[float](colo2_p)
    colmol = Ptr[float](colmol_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    indself = Ptr[int](indself_p)
    indfor = Ptr[int](indfor_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    selffac = Ptr[float](selffac_p)
    selffrac = Ptr[float](selffrac_p)
    forfac = Ptr[float](forfac_p)
    forfrac = Ptr[float](forfrac_p)
    absa = Ptr[float](absa_p)
    absb = Ptr[float](absb_p)
    selfref = Ptr[float](selfref_p)
    forref = Ptr[float](forref_p)
    sfluxref = Ptr[float](sfluxref_p)
    sfluxzen = Ptr[float](sfluxzen_p)
    taug = Ptr[float](taug_p)
    taur = Ptr[float](taur_p)
    o2adj = 1.6

    laysolfr = laytrop

    for lay in range(1, laytrop + 1):
        if jp[lay - 1] < layreffr and jp[lay] >= layreffr:
            if lay + 1 < laytrop:
                laysolfr = lay + 1
            else:
                laysolfr = laytrop
        o2cont = 4.35e-4 * colo2[lay - 1] / (350.0 * 2.0)
        speccomb = colh2o[lay - 1] + o2adj * strrat * colo2[lay - 1]
        specparm = colh2o[lay - 1] / speccomb
        if specparm >= oneminus:
            specparm = oneminus
        specmult = 8.0 * specparm
        js = 1 + int(specmult)
        fs = specmult - float(int(specmult))
        fac000 = (1.0 - fs) * fac00[lay - 1]
        fac010 = (1.0 - fs) * fac10[lay - 1]
        fac100 = fs * fac00[lay - 1]
        fac110 = fs * fac10[lay - 1]
        fac001 = (1.0 - fs) * fac01[lay - 1]
        fac011 = (1.0 - fs) * fac11[lay - 1]
        fac101 = fs * fac01[lay - 1]
        fac111 = fs * fac11[lay - 1]
        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa22 + js
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa22 + js
        inds = indself[lay - 1]
        indf = indfor[lay - 1]
        tauray = colmol[lay - 1] * rayl

        for ig in range(1, ng22 + 1):
            taug[_idx2(lay, ngs21 + ig, nlayers)] = (
                speccomb
                * (
                    fac000 * absa[_idx2(ind0, ig, 585)]
                    + fac100 * absa[_idx2(ind0 + 1, ig, 585)]
                    + fac010 * absa[_idx2(ind0 + 9, ig, 585)]
                    + fac110 * absa[_idx2(ind0 + 10, ig, 585)]
                    + fac001 * absa[_idx2(ind1, ig, 585)]
                    + fac101 * absa[_idx2(ind1 + 1, ig, 585)]
                    + fac011 * absa[_idx2(ind1 + 9, ig, 585)]
                    + fac111 * absa[_idx2(ind1 + 10, ig, 585)]
                )
                + colh2o[lay - 1]
                * (
                    selffac[lay - 1]
                    * (
                        selfref[_idx2(inds, ig, 10)]
                        + selffrac[lay - 1]
                        * (
                            selfref[_idx2(inds + 1, ig, 10)]
                            - selfref[_idx2(inds, ig, 10)]
                        )
                    )
                    + forfac[lay - 1]
                    * (
                        forref[_idx2(indf, ig, 3)]
                        + forfrac[lay - 1]
                        * (
                            forref[_idx2(indf + 1, ig, 3)]
                            - forref[_idx2(indf, ig, 3)]
                        )
                    )
                )
                + o2cont
            )
            if lay == laysolfr:
                sfluxzen[ngs21 + ig - 1] = sfluxref[_idx2(ig, js, ng22)] + fs * (
                    sfluxref[_idx2(ig, js + 1, ng22)] - sfluxref[_idx2(ig, js, ng22)]
                )
            taur[_idx2(lay, ngs21 + ig, nlayers)] = tauray

    for lay in range(laytrop + 1, nlayers + 1):
        o2cont = 4.35e-4 * colo2[lay - 1] / (350.0 * 2.0)
        ind0 = ((jp[lay - 1] - 13) * 5 + (jt[lay - 1] - 1)) * nspb22 + 1
        ind1 = ((jp[lay - 1] - 12) * 5 + (jt1[lay - 1] - 1)) * nspb22 + 1
        tauray = colmol[lay - 1] * rayl

        for ig in range(1, ng22 + 1):
            taug[_idx2(lay, ngs21 + ig, nlayers)] = (
                colo2[lay - 1]
                * o2adj
                * (
                    fac00[lay - 1] * absb[_idx2(ind0, ig, 235)]
                    + fac10[lay - 1] * absb[_idx2(ind0 + 1, ig, 235)]
                    + fac01[lay - 1] * absb[_idx2(ind1, ig, 235)]
                    + fac11[lay - 1] * absb[_idx2(ind1 + 1, ig, 235)]
                )
                + o2cont
            )
            taur[_idx2(lay, ngs21 + ig, nlayers)] = tauray


@export
def rrtmg_sw_subcol_prep_codon(
    ncol: int,
    nlay: int,
    ld_play: int,
    ld_size: int,
    ld_pmid: int,
    ld_out: int,
    play_p: cobj,
    rei_p: cobj,
    rel_p: cobj,
    pmid_p: cobj,
    reicmcl_p: cobj,
    relqmcl_p: cobj,
):
    play = Ptr[float](play_p)
    rei = Ptr[float](rei_p)
    rel = Ptr[float](rel_p)
    pmid = Ptr[float](pmid_p)
    reicmcl = Ptr[float](reicmcl_p)
    relqmcl = Ptr[float](relqmcl_p)

    for k in range(1, nlay + 1):
        for i in range(1, ncol + 1):
            pmid[_idx2(i, k, ld_pmid)] = play[_idx2(i, k, ld_play)] * 1.0e2
            reicmcl[_idx2(i, k, ld_out)] = rei[_idx2(i, k, ld_size)]
            relqmcl[_idx2(i, k, ld_out)] = rel[_idx2(i, k, ld_size)]


@export
def rrtmg_sw_subcol_cldf_prep_codon(
    ncol: int,
    nlay: int,
    ld_cld: int,
    ld_cldf: int,
    cldmin: float,
    cld_p: cobj,
    cldf_p: cobj,
):
    cld = Ptr[float](cld_p)
    cldf = Ptr[float](cldf_p)

    for k in range(1, nlay + 1):
        for i in range(1, ncol + 1):
            v = cld[_idx2(i, k, ld_cld)]
            if v < cldmin:
                v = 0.0
            cldf[_idx2(i, k, ld_cldf)] = v


@export
def rrtmg_sw_subcol_seed_init_codon(
    ncol: int,
    nlay: int,
    ld_pmid: int,
    pmid_p: cobj,
    seed1_p: cobj,
    seed2_p: cobj,
    seed3_p: cobj,
    seed4_p: cobj,
):
    pmid = Ptr[float](pmid_p)
    seed1 = Ptr[i32](seed1_p)
    seed2 = Ptr[i32](seed2_p)
    seed3 = Ptr[i32](seed3_p)
    seed4 = Ptr[i32](seed4_p)

    for i in range(1, ncol + 1):
        p1 = pmid[_idx2(i, nlay, ld_pmid)]
        p2 = pmid[_idx2(i, nlay - 1, ld_pmid)]
        p3 = pmid[_idx2(i, nlay - 2, ld_pmid)]
        p4 = pmid[_idx2(i, nlay - 3, ld_pmid)]
        seed1[i - 1] = i32(int((p1 - float(int(p1))) * 1000000000.0))
        seed2[i - 1] = i32(int((p2 - float(int(p2))) * 1000000000.0))
        seed3[i - 1] = i32(int((p3 - float(int(p3))) * 1000000000.0))
        seed4[i - 1] = i32(int((p4 - float(int(p4))) * 1000000000.0))


@inline
def _rrtmg_ishft_i32(k: i32, n: int) -> i32:
    if n >= 0:
        return i32(u32(k) << u32(n))
    return i32(u32(k) >> u32(-n))


@inline
def _rrtmg_ieor_ishft_i32(k: i32, n: int) -> i32:
    return i32(u32(k) ^ u32(_rrtmg_ishft_i32(k, n)))


@inline
def _rrtmg_low_i32(i: i64) -> i32:
    return i32(i)


@inline
def _rrtmg_kiss_next_i32(
    i0: int,
    seed1: Ptr[i32],
    seed2: Ptr[i32],
    seed3: Ptr[i32],
    seed4: Ptr[i32],
) -> float:
    kiss = i64(69069) * i64(seed1[i0]) + i64(1327217885)
    seed1[i0] = _rrtmg_low_i32(kiss)
    seed2[i0] = _rrtmg_ieor_ishft_i32(
        _rrtmg_ieor_ishft_i32(_rrtmg_ieor_ishft_i32(seed2[i0], 13), -17), 5
    )
    seed3[i0] = i32(
        18000 * int(u32(seed3[i0]) & u32(65535)) + int(u32(seed3[i0]) >> u32(16))
    )
    seed4[i0] = i32(
        30903 * int(u32(seed4[i0]) & u32(65535)) + int(u32(seed4[i0]) >> u32(16))
    )
    kiss = (
        i64(seed1[i0])
        + i64(seed2[i0])
        + i64(_rrtmg_ishft_i32(seed3[i0], 16))
        + i64(seed4[i0])
    )
    return float(_rrtmg_low_i32(kiss)) * 2.328306e-10 + 0.5


@export
def rrtmg_sw_subcol_kiss_random_case2_codon(
    ncol: int,
    nlay: int,
    nsubcol: int,
    seed1_p: cobj,
    seed2_p: cobj,
    seed3_p: cobj,
    seed4_p: cobj,
    cdf_p: cobj,
):
    seed1 = Ptr[i32](seed1_p)
    seed2 = Ptr[i32](seed2_p)
    seed3 = Ptr[i32](seed3_p)
    seed4 = Ptr[i32](seed4_p)
    cdf = Ptr[float](cdf_p)

    for isubcol in range(1, nsubcol + 1):
        for ilev in range(1, nlay + 1):
            for i in range(1, ncol + 1):
                cdf[_idx3(isubcol, i, ilev, nsubcol, ncol)] = _rrtmg_kiss_next_i32(
                    i - 1, seed1, seed2, seed3, seed4
                )


@export
def rrtmg_sw_subcol_kiss_advance_codon(
    ncol: int,
    nstep: int,
    seed1_p: cobj,
    seed2_p: cobj,
    seed3_p: cobj,
    seed4_p: cobj,
):
    seed1 = Ptr[i32](seed1_p)
    seed2 = Ptr[i32](seed2_p)
    seed3 = Ptr[i32](seed3_p)
    seed4 = Ptr[i32](seed4_p)

    for _ in range(1, nstep + 1):
        for i in range(1, ncol + 1):
            _rrtmg_kiss_next_i32(i - 1, seed1, seed2, seed3, seed4)


@export
def rrtmg_sw_subcol_overlap_codon(
    ncol: int,
    nlay: int,
    nsubcol: int,
    cdf_p: cobj,
    cldf_p: cobj,
):
    cdf = Ptr[float](cdf_p)
    cldf = Ptr[float](cldf_p)

    for ilev in range(2, nlay + 1):
        for i in range(1, ncol + 1):
            cldf_prev_idx = _idx2(i, ilev - 1, ncol)
            for isubcol in range(1, nsubcol + 1):
                prev_idx = _idx3(isubcol, i, ilev - 1, nsubcol, ncol)
                cur_idx = _idx3(isubcol, i, ilev, nsubcol, ncol)
                if cdf[prev_idx] > 1.0 - cldf[cldf_prev_idx]:
                    cdf[cur_idx] = cdf[prev_idx]
                else:
                    cdf[cur_idx] = cdf[cur_idx] * (1.0 - cldf[cldf_prev_idx])


@export
def rrtmg_sw_subcol_fill_codon(
    ncol: int,
    nlay: int,
    nsubcol: int,
    nbndsw: int,
    ld_cloud: int,
    ld_tauc_col: int,
    ld_out_col: int,
    cdf_p: cobj,
    cldf_p: cobj,
    clwp_p: cobj,
    ciwp_p: cobj,
    tauc_p: cobj,
    ssac_p: cobj,
    asmc_p: cobj,
    fsfc_p: cobj,
    ngb_p: cobj,
    cld_stoch_p: cobj,
    clwp_stoch_p: cobj,
    ciwp_stoch_p: cobj,
    tauc_stoch_p: cobj,
    ssac_stoch_p: cobj,
    asmc_stoch_p: cobj,
    fsfc_stoch_p: cobj,
):
    cdf = Ptr[float](cdf_p)
    cldf = Ptr[float](cldf_p)
    clwp = Ptr[float](clwp_p)
    ciwp = Ptr[float](ciwp_p)
    tauc = Ptr[float](tauc_p)
    ssac = Ptr[float](ssac_p)
    asmc = Ptr[float](asmc_p)
    fsfc = Ptr[float](fsfc_p)
    ngb = Ptr[int](ngb_p)
    cld_stoch = Ptr[float](cld_stoch_p)
    clwp_stoch = Ptr[float](clwp_stoch_p)
    ciwp_stoch = Ptr[float](ciwp_stoch_p)
    tauc_stoch = Ptr[float](tauc_stoch_p)
    ssac_stoch = Ptr[float](ssac_stoch_p)
    asmc_stoch = Ptr[float](asmc_stoch_p)
    fsfc_stoch = Ptr[float](fsfc_stoch_p)

    for ilev in range(1, nlay + 1):
        for i in range(1, ncol + 1):
            threshold = 1.0 - cldf[_idx2(i, ilev, ncol)]
            for isubcol in range(1, nsubcol + 1):
                cdf_idx = _idx3(isubcol, i, ilev, nsubcol, ncol)
                out_idx = _idx3(isubcol, i, ilev, nsubcol, ld_out_col)
                if cdf[cdf_idx] >= threshold:
                    cld_stoch[out_idx] = 1.0
                else:
                    cld_stoch[out_idx] = 0.0

    for ilev in range(1, nlay + 1):
        for i in range(1, ncol + 1):
            cfrac = cldf[_idx2(i, ilev, ncol)]
            threshold = 1.0 - cfrac
            for isubcol in range(1, nsubcol + 1):
                cdf_idx = _idx3(isubcol, i, ilev, nsubcol, ncol)
                out_idx = _idx3(isubcol, i, ilev, nsubcol, ld_out_col)
                if cdf[cdf_idx] >= threshold and cfrac > 0.0:
                    clwp_stoch[out_idx] = clwp[_idx2(i, ilev, ld_cloud)]
                    ciwp_stoch[out_idx] = ciwp[_idx2(i, ilev, ld_cloud)]
                else:
                    clwp_stoch[out_idx] = 0.0
                    ciwp_stoch[out_idx] = 0.0

    ngbm = ngb[0] - 1
    for ilev in range(1, nlay + 1):
        for i in range(1, ncol + 1):
            cfrac = cldf[_idx2(i, ilev, ncol)]
            threshold = 1.0 - cfrac
            for isubcol in range(1, nsubcol + 1):
                cdf_idx = _idx3(isubcol, i, ilev, nsubcol, ncol)
                out_idx = _idx3(isubcol, i, ilev, nsubcol, ld_out_col)
                if cdf[cdf_idx] >= threshold and cfrac > 0.0:
                    n = ngb[isubcol - 1] - ngbm
                    src_idx = _idx3(n, i, ilev, nbndsw, ld_tauc_col)
                    tauc_stoch[out_idx] = tauc[src_idx]
                    ssac_stoch[out_idx] = ssac[src_idx]
                    asmc_stoch[out_idx] = asmc[src_idx]
                    fsfc_stoch[out_idx] = fsfc[src_idx]
                else:
                    tauc_stoch[out_idx] = 0.0
                    ssac_stoch[out_idx] = 1.0
                    asmc_stoch[out_idx] = 0.0
                    fsfc_stoch[out_idx] = 0.0


@export
def rrtmg_lw_subcol_prep_codon(
    ncol: int,
    nlay: int,
    ld_play: int,
    ld_size: int,
    ld_pmid: int,
    ld_out: int,
    play_p: cobj,
    rei_p: cobj,
    rel_p: cobj,
    pmid_p: cobj,
    reicmcl_p: cobj,
    relqmcl_p: cobj,
):
    play = Ptr[float](play_p)
    rei = Ptr[float](rei_p)
    rel = Ptr[float](rel_p)
    pmid = Ptr[float](pmid_p)
    reicmcl = Ptr[float](reicmcl_p)
    relqmcl = Ptr[float](relqmcl_p)

    for k in range(1, nlay + 1):
        for i in range(1, ncol + 1):
            pmid[_idx2(i, k, ld_pmid)] = play[_idx2(i, k, ld_play)] * 1.0e2
            reicmcl[_idx2(i, k, ld_out)] = rei[_idx2(i, k, ld_size)]
            relqmcl[_idx2(i, k, ld_out)] = rel[_idx2(i, k, ld_size)]


@export
def rrtmg_lw_subcol_cldf_prep_codon(
    ncol: int,
    nlay: int,
    ld_cld: int,
    ld_cldf: int,
    cldmin: float,
    cld_p: cobj,
    cldf_p: cobj,
):
    cld = Ptr[float](cld_p)
    cldf = Ptr[float](cldf_p)

    for k in range(1, nlay + 1):
        for i in range(1, ncol + 1):
            v = cld[_idx2(i, k, ld_cld)]
            if v < cldmin:
                v = 0.0
            cldf[_idx2(i, k, ld_cldf)] = v


@export
def rrtmg_lw_subcol_seed_init_codon(
    ncol: int,
    nlay: int,
    ld_pmid: int,
    pmid_p: cobj,
    seed1_p: cobj,
    seed2_p: cobj,
    seed3_p: cobj,
    seed4_p: cobj,
):
    pmid = Ptr[float](pmid_p)
    seed1 = Ptr[i32](seed1_p)
    seed2 = Ptr[i32](seed2_p)
    seed3 = Ptr[i32](seed3_p)
    seed4 = Ptr[i32](seed4_p)

    for i in range(1, ncol + 1):
        p1 = pmid[_idx2(i, nlay, ld_pmid)]
        p2 = pmid[_idx2(i, nlay - 1, ld_pmid)]
        p3 = pmid[_idx2(i, nlay - 2, ld_pmid)]
        p4 = pmid[_idx2(i, nlay - 3, ld_pmid)]
        seed1[i - 1] = i32(int((p1 - float(int(p1))) * 1000000000.0))
        seed2[i - 1] = i32(int((p2 - float(int(p2))) * 1000000000.0))
        seed3[i - 1] = i32(int((p3 - float(int(p3))) * 1000000000.0))
        seed4[i - 1] = i32(int((p4 - float(int(p4))) * 1000000000.0))


@export
def rrtmg_lw_subcol_kiss_random_case2_codon(
    ncol: int,
    nlay: int,
    nsubcol: int,
    seed1_p: cobj,
    seed2_p: cobj,
    seed3_p: cobj,
    seed4_p: cobj,
    cdf_p: cobj,
):
    seed1 = Ptr[i32](seed1_p)
    seed2 = Ptr[i32](seed2_p)
    seed3 = Ptr[i32](seed3_p)
    seed4 = Ptr[i32](seed4_p)
    cdf = Ptr[float](cdf_p)

    for isubcol in range(1, nsubcol + 1):
        for ilev in range(1, nlay + 1):
            for i in range(1, ncol + 1):
                cdf[_idx3(isubcol, i, ilev, nsubcol, ncol)] = _rrtmg_kiss_next_i32(
                    i - 1, seed1, seed2, seed3, seed4
                )


@export
def rrtmg_lw_subcol_kiss_advance_codon(
    ncol: int,
    nstep: int,
    seed1_p: cobj,
    seed2_p: cobj,
    seed3_p: cobj,
    seed4_p: cobj,
):
    seed1 = Ptr[i32](seed1_p)
    seed2 = Ptr[i32](seed2_p)
    seed3 = Ptr[i32](seed3_p)
    seed4 = Ptr[i32](seed4_p)

    for _ in range(1, nstep + 1):
        for i in range(1, ncol + 1):
            _rrtmg_kiss_next_i32(i - 1, seed1, seed2, seed3, seed4)


@export
def rrtmg_lw_subcol_overlap_codon(
    ncol: int,
    nlay: int,
    nsubcol: int,
    cdf_p: cobj,
    cldf_p: cobj,
):
    cdf = Ptr[float](cdf_p)
    cldf = Ptr[float](cldf_p)

    for ilev in range(2, nlay + 1):
        for i in range(1, ncol + 1):
            cldf_prev_idx = _idx2(i, ilev - 1, ncol)
            for isubcol in range(1, nsubcol + 1):
                prev_idx = _idx3(isubcol, i, ilev - 1, nsubcol, ncol)
                cur_idx = _idx3(isubcol, i, ilev, nsubcol, ncol)
                if cdf[prev_idx] > 1.0 - cldf[cldf_prev_idx]:
                    cdf[cur_idx] = cdf[prev_idx]
                else:
                    cdf[cur_idx] = cdf[cur_idx] * (1.0 - cldf[cldf_prev_idx])


@export
def rrtmg_lw_subcol_fill_codon(
    ncol: int,
    nlay: int,
    nsubcol: int,
    nbndlw: int,
    ld_cloud: int,
    ld_tauc_col: int,
    ld_out_col: int,
    cdf_p: cobj,
    cldf_p: cobj,
    clwp_p: cobj,
    ciwp_p: cobj,
    tauc_p: cobj,
    ngb_p: cobj,
    cld_stoch_p: cobj,
    clwp_stoch_p: cobj,
    ciwp_stoch_p: cobj,
    tauc_stoch_p: cobj,
):
    cdf = Ptr[float](cdf_p)
    cldf = Ptr[float](cldf_p)
    clwp = Ptr[float](clwp_p)
    ciwp = Ptr[float](ciwp_p)
    tauc = Ptr[float](tauc_p)
    ngb = Ptr[int](ngb_p)
    cld_stoch = Ptr[float](cld_stoch_p)
    clwp_stoch = Ptr[float](clwp_stoch_p)
    ciwp_stoch = Ptr[float](ciwp_stoch_p)
    tauc_stoch = Ptr[float](tauc_stoch_p)

    for ilev in range(1, nlay + 1):
        for i in range(1, ncol + 1):
            threshold = 1.0 - cldf[_idx2(i, ilev, ncol)]
            for isubcol in range(1, nsubcol + 1):
                cdf_idx = _idx3(isubcol, i, ilev, nsubcol, ncol)
                out_idx = _idx3(isubcol, i, ilev, nsubcol, ld_out_col)
                if cdf[cdf_idx] >= threshold:
                    cld_stoch[out_idx] = 1.0
                else:
                    cld_stoch[out_idx] = 0.0

    for ilev in range(1, nlay + 1):
        for i in range(1, ncol + 1):
            cfrac = cldf[_idx2(i, ilev, ncol)]
            threshold = 1.0 - cfrac
            for isubcol in range(1, nsubcol + 1):
                cdf_idx = _idx3(isubcol, i, ilev, nsubcol, ncol)
                out_idx = _idx3(isubcol, i, ilev, nsubcol, ld_out_col)
                if cdf[cdf_idx] >= threshold and cfrac > 0.0:
                    clwp_stoch[out_idx] = clwp[_idx2(i, ilev, ld_cloud)]
                    ciwp_stoch[out_idx] = ciwp[_idx2(i, ilev, ld_cloud)]
                else:
                    clwp_stoch[out_idx] = 0.0
                    ciwp_stoch[out_idx] = 0.0

    for ilev in range(1, nlay + 1):
        for i in range(1, ncol + 1):
            cfrac = cldf[_idx2(i, ilev, ncol)]
            threshold = 1.0 - cfrac
            for isubcol in range(1, nsubcol + 1):
                cdf_idx = _idx3(isubcol, i, ilev, nsubcol, ncol)
                out_idx = _idx3(isubcol, i, ilev, nsubcol, ld_out_col)
                if cdf[cdf_idx] >= threshold and cfrac > 0.0:
                    n = ngb[isubcol - 1]
                    tauc_stoch[out_idx] = tauc[_idx3(n, i, ilev, nbndlw, ld_tauc_col)]
                else:
                    tauc_stoch[out_idx] = 0.0


@export
def rrtmg_lw_taugb10_codon(
    nlayers: int,
    laytrop: int,
    ng10: int,
    ngs9: int,
    nspa10: int,
    nspb10: int,
    colh2o_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    indself_p: cobj,
    indfor_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    selffac_p: cobj,
    selffrac_p: cobj,
    forfac_p: cobj,
    forfrac_p: cobj,
    fracrefa_p: cobj,
    fracrefb_p: cobj,
    absa_p: cobj,
    absb_p: cobj,
    selfref_p: cobj,
    forref_p: cobj,
    fracs_p: cobj,
    taug_p: cobj,
):
    colh2o = Ptr[float](colh2o_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    indself = Ptr[int](indself_p)
    indfor = Ptr[int](indfor_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    selffac = Ptr[float](selffac_p)
    selffrac = Ptr[float](selffrac_p)
    forfac = Ptr[float](forfac_p)
    forfrac = Ptr[float](forfrac_p)
    fracrefa = Ptr[float](fracrefa_p)
    fracrefb = Ptr[float](fracrefb_p)
    absa = Ptr[float](absa_p)
    absb = Ptr[float](absb_p)
    selfref = Ptr[float](selfref_p)
    forref = Ptr[float](forref_p)
    fracs = Ptr[float](fracs_p)
    taug = Ptr[float](taug_p)

    for lay in range(1, laytrop + 1):
        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa10 + 1
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa10 + 1
        inds = indself[lay - 1]
        indf = indfor[lay - 1]

        for ig in range(1, ng10 + 1):
            tauself = selffac[lay - 1] * (
                selfref[_idx2(inds, ig, 10)]
                + selffrac[lay - 1]
                * (selfref[_idx2(inds + 1, ig, 10)] - selfref[_idx2(inds, ig, 10)])
            )
            taufor = forfac[lay - 1] * (
                forref[_idx2(indf, ig, 4)]
                + forfrac[lay - 1]
                * (forref[_idx2(indf + 1, ig, 4)] - forref[_idx2(indf, ig, 4)])
            )
            taug[_idx2(lay, ngs9 + ig, nlayers)] = colh2o[lay - 1] * (
                fac00[lay - 1] * absa[_idx2(ind0, ig, 65)]
                + fac10[lay - 1] * absa[_idx2(ind0 + 1, ig, 65)]
                + fac01[lay - 1] * absa[_idx2(ind1, ig, 65)]
                + fac11[lay - 1] * absa[_idx2(ind1 + 1, ig, 65)]
            ) + tauself + taufor
            fracs[_idx2(lay, ngs9 + ig, nlayers)] = fracrefa[ig - 1]

    for lay in range(laytrop + 1, nlayers + 1):
        ind0 = ((jp[lay - 1] - 13) * 5 + (jt[lay - 1] - 1)) * nspb10 + 1
        ind1 = ((jp[lay - 1] - 12) * 5 + (jt1[lay - 1] - 1)) * nspb10 + 1
        indf = indfor[lay - 1]

        for ig in range(1, ng10 + 1):
            taufor = forfac[lay - 1] * (
                forref[_idx2(indf, ig, 4)]
                + forfrac[lay - 1]
                * (forref[_idx2(indf + 1, ig, 4)] - forref[_idx2(indf, ig, 4)])
            )
            taug[_idx2(lay, ngs9 + ig, nlayers)] = colh2o[lay - 1] * (
                fac00[lay - 1] * absb[_idx2(ind0, ig, 235)]
                + fac10[lay - 1] * absb[_idx2(ind0 + 1, ig, 235)]
                + fac01[lay - 1] * absb[_idx2(ind1, ig, 235)]
                + fac11[lay - 1] * absb[_idx2(ind1 + 1, ig, 235)]
            ) + taufor
            fracs[_idx2(lay, ngs9 + ig, nlayers)] = fracrefb[ig - 1]


@export
def rrtmg_lw_taugb11_codon(
    nlayers: int,
    laytrop: int,
    ng11: int,
    ngs10: int,
    nspa11: int,
    nspb11: int,
    colh2o_p: cobj,
    colo2_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    indself_p: cobj,
    indfor_p: cobj,
    indminor_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    selffac_p: cobj,
    selffrac_p: cobj,
    forfac_p: cobj,
    forfrac_p: cobj,
    scaleminor_p: cobj,
    minorfrac_p: cobj,
    fracrefa_p: cobj,
    fracrefb_p: cobj,
    absa_p: cobj,
    absb_p: cobj,
    ka_mo2_p: cobj,
    kb_mo2_p: cobj,
    selfref_p: cobj,
    forref_p: cobj,
    fracs_p: cobj,
    taug_p: cobj,
):
    colh2o = Ptr[float](colh2o_p)
    colo2 = Ptr[float](colo2_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    indself = Ptr[int](indself_p)
    indfor = Ptr[int](indfor_p)
    indminor = Ptr[int](indminor_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    selffac = Ptr[float](selffac_p)
    selffrac = Ptr[float](selffrac_p)
    forfac = Ptr[float](forfac_p)
    forfrac = Ptr[float](forfrac_p)
    scaleminor = Ptr[float](scaleminor_p)
    minorfrac = Ptr[float](minorfrac_p)
    fracrefa = Ptr[float](fracrefa_p)
    fracrefb = Ptr[float](fracrefb_p)
    absa = Ptr[float](absa_p)
    absb = Ptr[float](absb_p)
    ka_mo2 = Ptr[float](ka_mo2_p)
    kb_mo2 = Ptr[float](kb_mo2_p)
    selfref = Ptr[float](selfref_p)
    forref = Ptr[float](forref_p)
    fracs = Ptr[float](fracs_p)
    taug = Ptr[float](taug_p)

    for lay in range(1, laytrop + 1):
        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa11 + 1
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa11 + 1
        inds = indself[lay - 1]
        indf = indfor[lay - 1]
        indm = indminor[lay - 1]
        scaleo2 = colo2[lay - 1] * scaleminor[lay - 1]
        for ig in range(1, ng11 + 1):
            tauself = selffac[lay - 1] * (
                selfref[_idx2(inds, ig, 10)]
                + selffrac[lay - 1]
                * (selfref[_idx2(inds + 1, ig, 10)] - selfref[_idx2(inds, ig, 10)])
            )
            taufor = forfac[lay - 1] * (
                forref[_idx2(indf, ig, 4)]
                + forfrac[lay - 1]
                * (forref[_idx2(indf + 1, ig, 4)] - forref[_idx2(indf, ig, 4)])
            )
            tauo2 = scaleo2 * (
                ka_mo2[_idx2(indm, ig, 19)]
                + minorfrac[lay - 1]
                * (ka_mo2[_idx2(indm + 1, ig, 19)] - ka_mo2[_idx2(indm, ig, 19)])
            )
            taug[_idx2(lay, ngs10 + ig, nlayers)] = colh2o[lay - 1] * (
                fac00[lay - 1] * absa[_idx2(ind0, ig, 65)]
                + fac10[lay - 1] * absa[_idx2(ind0 + 1, ig, 65)]
                + fac01[lay - 1] * absa[_idx2(ind1, ig, 65)]
                + fac11[lay - 1] * absa[_idx2(ind1 + 1, ig, 65)]
            ) + tauself + taufor + tauo2
            fracs[_idx2(lay, ngs10 + ig, nlayers)] = fracrefa[ig - 1]

    for lay in range(laytrop + 1, nlayers + 1):
        ind0 = ((jp[lay - 1] - 13) * 5 + (jt[lay - 1] - 1)) * nspb11 + 1
        ind1 = ((jp[lay - 1] - 12) * 5 + (jt1[lay - 1] - 1)) * nspb11 + 1
        indf = indfor[lay - 1]
        indm = indminor[lay - 1]
        scaleo2 = colo2[lay - 1] * scaleminor[lay - 1]
        for ig in range(1, ng11 + 1):
            taufor = forfac[lay - 1] * (
                forref[_idx2(indf, ig, 4)]
                + forfrac[lay - 1]
                * (forref[_idx2(indf + 1, ig, 4)] - forref[_idx2(indf, ig, 4)])
            )
            tauo2 = scaleo2 * (
                kb_mo2[_idx2(indm, ig, 19)]
                + minorfrac[lay - 1]
                * (kb_mo2[_idx2(indm + 1, ig, 19)] - kb_mo2[_idx2(indm, ig, 19)])
            )
            taug[_idx2(lay, ngs10 + ig, nlayers)] = colh2o[lay - 1] * (
                fac00[lay - 1] * absb[_idx2(ind0, ig, 235)]
                + fac10[lay - 1] * absb[_idx2(ind0 + 1, ig, 235)]
                + fac01[lay - 1] * absb[_idx2(ind1, ig, 235)]
                + fac11[lay - 1] * absb[_idx2(ind1 + 1, ig, 235)]
            ) + taufor + tauo2
            fracs[_idx2(lay, ngs10 + ig, nlayers)] = fracrefb[ig - 1]


@export
def rrtmg_lw_taugb14_codon(
    nlayers: int,
    laytrop: int,
    ng14: int,
    ngs13: int,
    nspa14: int,
    nspb14: int,
    colco2_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    indself_p: cobj,
    indfor_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    selffac_p: cobj,
    selffrac_p: cobj,
    forfac_p: cobj,
    forfrac_p: cobj,
    fracrefa_p: cobj,
    fracrefb_p: cobj,
    absa_p: cobj,
    absb_p: cobj,
    selfref_p: cobj,
    forref_p: cobj,
    fracs_p: cobj,
    taug_p: cobj,
):
    colco2 = Ptr[float](colco2_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    indself = Ptr[int](indself_p)
    indfor = Ptr[int](indfor_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    selffac = Ptr[float](selffac_p)
    selffrac = Ptr[float](selffrac_p)
    forfac = Ptr[float](forfac_p)
    forfrac = Ptr[float](forfrac_p)
    fracrefa = Ptr[float](fracrefa_p)
    fracrefb = Ptr[float](fracrefb_p)
    absa = Ptr[float](absa_p)
    absb = Ptr[float](absb_p)
    selfref = Ptr[float](selfref_p)
    forref = Ptr[float](forref_p)
    fracs = Ptr[float](fracs_p)
    taug = Ptr[float](taug_p)

    for lay in range(1, laytrop + 1):
        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa14 + 1
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa14 + 1
        inds = indself[lay - 1]
        indf = indfor[lay - 1]
        for ig in range(1, ng14 + 1):
            tauself = selffac[lay - 1] * (
                selfref[_idx2(inds, ig, 10)]
                + selffrac[lay - 1]
                * (selfref[_idx2(inds + 1, ig, 10)] - selfref[_idx2(inds, ig, 10)])
            )
            taufor = forfac[lay - 1] * (
                forref[_idx2(indf, ig, 4)]
                + forfrac[lay - 1]
                * (forref[_idx2(indf + 1, ig, 4)] - forref[_idx2(indf, ig, 4)])
            )
            taug[_idx2(lay, ngs13 + ig, nlayers)] = colco2[lay - 1] * (
                fac00[lay - 1] * absa[_idx2(ind0, ig, 65)]
                + fac10[lay - 1] * absa[_idx2(ind0 + 1, ig, 65)]
                + fac01[lay - 1] * absa[_idx2(ind1, ig, 65)]
                + fac11[lay - 1] * absa[_idx2(ind1 + 1, ig, 65)]
            ) + tauself + taufor
            fracs[_idx2(lay, ngs13 + ig, nlayers)] = fracrefa[ig - 1]

    for lay in range(laytrop + 1, nlayers + 1):
        ind0 = ((jp[lay - 1] - 13) * 5 + (jt[lay - 1] - 1)) * nspb14 + 1
        ind1 = ((jp[lay - 1] - 12) * 5 + (jt1[lay - 1] - 1)) * nspb14 + 1
        for ig in range(1, ng14 + 1):
            taug[_idx2(lay, ngs13 + ig, nlayers)] = colco2[lay - 1] * (
                fac00[lay - 1] * absb[_idx2(ind0, ig, 235)]
                + fac10[lay - 1] * absb[_idx2(ind0 + 1, ig, 235)]
                + fac01[lay - 1] * absb[_idx2(ind1, ig, 235)]
                + fac11[lay - 1] * absb[_idx2(ind1 + 1, ig, 235)]
            )
            fracs[_idx2(lay, ngs13 + ig, nlayers)] = fracrefb[ig - 1]


@export
def rrtmg_lw_taugb1_codon(
    nlayers: int,
    laytrop: int,
    ng1: int,
    nspa1: int,
    nspb1: int,
    pavel_p: cobj,
    colh2o_p: cobj,
    colbrd_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    indself_p: cobj,
    indfor_p: cobj,
    indminor_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    selffac_p: cobj,
    selffrac_p: cobj,
    forfac_p: cobj,
    forfrac_p: cobj,
    scaleminorn2_p: cobj,
    minorfrac_p: cobj,
    fracrefa_p: cobj,
    fracrefb_p: cobj,
    absa_p: cobj,
    absb_p: cobj,
    ka_mn2_p: cobj,
    kb_mn2_p: cobj,
    selfref_p: cobj,
    forref_p: cobj,
    fracs_p: cobj,
    taug_p: cobj,
):
    pavel = Ptr[float](pavel_p)
    colh2o = Ptr[float](colh2o_p)
    colbrd = Ptr[float](colbrd_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    indself = Ptr[int](indself_p)
    indfor = Ptr[int](indfor_p)
    indminor = Ptr[int](indminor_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    selffac = Ptr[float](selffac_p)
    selffrac = Ptr[float](selffrac_p)
    forfac = Ptr[float](forfac_p)
    forfrac = Ptr[float](forfrac_p)
    scaleminorn2 = Ptr[float](scaleminorn2_p)
    minorfrac = Ptr[float](minorfrac_p)
    fracrefa = Ptr[float](fracrefa_p)
    fracrefb = Ptr[float](fracrefb_p)
    absa = Ptr[float](absa_p)
    absb = Ptr[float](absb_p)
    ka_mn2 = Ptr[float](ka_mn2_p)
    kb_mn2 = Ptr[float](kb_mn2_p)
    selfref = Ptr[float](selfref_p)
    forref = Ptr[float](forref_p)
    fracs = Ptr[float](fracs_p)
    taug = Ptr[float](taug_p)

    for lay in range(1, laytrop + 1):
        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa1 + 1
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa1 + 1
        inds = indself[lay - 1]
        indf = indfor[lay - 1]
        indm = indminor[lay - 1]
        pp = pavel[lay - 1]
        corradj = 1.0
        if pp < 250.0:
            corradj = 1.0 - 0.15 * (250.0 - pp) / 154.4

        scalen2 = colbrd[lay - 1] * scaleminorn2[lay - 1]
        for ig in range(1, ng1 + 1):
            tauself = selffac[lay - 1] * (
                selfref[_idx2(inds, ig, 10)]
                + selffrac[lay - 1]
                * (selfref[_idx2(inds + 1, ig, 10)] - selfref[_idx2(inds, ig, 10)])
            )
            taufor = forfac[lay - 1] * (
                forref[_idx2(indf, ig, 4)]
                + forfrac[lay - 1]
                * (forref[_idx2(indf + 1, ig, 4)] - forref[_idx2(indf, ig, 4)])
            )
            taun2 = scalen2 * (
                ka_mn2[_idx2(indm, ig, 19)]
                + minorfrac[lay - 1]
                * (ka_mn2[_idx2(indm + 1, ig, 19)] - ka_mn2[_idx2(indm, ig, 19)])
            )
            taug[_idx2(lay, ig, nlayers)] = corradj * (
                colh2o[lay - 1]
                * (
                    fac00[lay - 1] * absa[_idx2(ind0, ig, 65)]
                    + fac10[lay - 1] * absa[_idx2(ind0 + 1, ig, 65)]
                    + fac01[lay - 1] * absa[_idx2(ind1, ig, 65)]
                    + fac11[lay - 1] * absa[_idx2(ind1 + 1, ig, 65)]
                )
                + tauself
                + taufor
                + taun2
            )
            fracs[_idx2(lay, ig, nlayers)] = fracrefa[ig - 1]

    for lay in range(laytrop + 1, nlayers + 1):
        ind0 = ((jp[lay - 1] - 13) * 5 + (jt[lay - 1] - 1)) * nspb1 + 1
        ind1 = ((jp[lay - 1] - 12) * 5 + (jt1[lay - 1] - 1)) * nspb1 + 1
        indf = indfor[lay - 1]
        indm = indminor[lay - 1]
        pp = pavel[lay - 1]
        corradj = 1.0 - 0.15 * (pp / 95.6)

        scalen2 = colbrd[lay - 1] * scaleminorn2[lay - 1]
        for ig in range(1, ng1 + 1):
            taufor = forfac[lay - 1] * (
                forref[_idx2(indf, ig, 4)]
                + forfrac[lay - 1]
                * (forref[_idx2(indf + 1, ig, 4)] - forref[_idx2(indf, ig, 4)])
            )
            taun2 = scalen2 * (
                kb_mn2[_idx2(indm, ig, 19)]
                + minorfrac[lay - 1]
                * (kb_mn2[_idx2(indm + 1, ig, 19)] - kb_mn2[_idx2(indm, ig, 19)])
            )
            taug[_idx2(lay, ig, nlayers)] = corradj * (
                colh2o[lay - 1]
                * (
                    fac00[lay - 1] * absb[_idx2(ind0, ig, 235)]
                    + fac10[lay - 1] * absb[_idx2(ind0 + 1, ig, 235)]
                    + fac01[lay - 1] * absb[_idx2(ind1, ig, 235)]
                    + fac11[lay - 1] * absb[_idx2(ind1 + 1, ig, 235)]
                )
                + taufor
                + taun2
            )
            fracs[_idx2(lay, ig, nlayers)] = fracrefb[ig - 1]


@export
def rrtmg_lw_taugb2_codon(
    nlayers: int,
    laytrop: int,
    ng2: int,
    ngs1: int,
    nspa2: int,
    nspb2: int,
    pavel_p: cobj,
    colh2o_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    indself_p: cobj,
    indfor_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    selffac_p: cobj,
    selffrac_p: cobj,
    forfac_p: cobj,
    forfrac_p: cobj,
    fracrefa_p: cobj,
    fracrefb_p: cobj,
    absa_p: cobj,
    absb_p: cobj,
    selfref_p: cobj,
    forref_p: cobj,
    fracs_p: cobj,
    taug_p: cobj,
):
    pavel = Ptr[float](pavel_p)
    colh2o = Ptr[float](colh2o_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    indself = Ptr[int](indself_p)
    indfor = Ptr[int](indfor_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    selffac = Ptr[float](selffac_p)
    selffrac = Ptr[float](selffrac_p)
    forfac = Ptr[float](forfac_p)
    forfrac = Ptr[float](forfrac_p)
    fracrefa = Ptr[float](fracrefa_p)
    fracrefb = Ptr[float](fracrefb_p)
    absa = Ptr[float](absa_p)
    absb = Ptr[float](absb_p)
    selfref = Ptr[float](selfref_p)
    forref = Ptr[float](forref_p)
    fracs = Ptr[float](fracs_p)
    taug = Ptr[float](taug_p)

    for lay in range(1, laytrop + 1):
        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa2 + 1
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa2 + 1
        inds = indself[lay - 1]
        indf = indfor[lay - 1]
        pp = pavel[lay - 1]
        corradj = 1.0 - 0.05 * (pp - 100.0) / 900.0
        for ig in range(1, ng2 + 1):
            tauself = selffac[lay - 1] * (
                selfref[_idx2(inds, ig, 10)]
                + selffrac[lay - 1]
                * (selfref[_idx2(inds + 1, ig, 10)] - selfref[_idx2(inds, ig, 10)])
            )
            taufor = forfac[lay - 1] * (
                forref[_idx2(indf, ig, 4)]
                + forfrac[lay - 1]
                * (forref[_idx2(indf + 1, ig, 4)] - forref[_idx2(indf, ig, 4)])
            )
            taug[_idx2(lay, ngs1 + ig, nlayers)] = corradj * (
                colh2o[lay - 1]
                * (
                    fac00[lay - 1] * absa[_idx2(ind0, ig, 65)]
                    + fac10[lay - 1] * absa[_idx2(ind0 + 1, ig, 65)]
                    + fac01[lay - 1] * absa[_idx2(ind1, ig, 65)]
                    + fac11[lay - 1] * absa[_idx2(ind1 + 1, ig, 65)]
                )
                + tauself
                + taufor
            )
            fracs[_idx2(lay, ngs1 + ig, nlayers)] = fracrefa[ig - 1]

    for lay in range(laytrop + 1, nlayers + 1):
        ind0 = ((jp[lay - 1] - 13) * 5 + (jt[lay - 1] - 1)) * nspb2 + 1
        ind1 = ((jp[lay - 1] - 12) * 5 + (jt1[lay - 1] - 1)) * nspb2 + 1
        indf = indfor[lay - 1]
        for ig in range(1, ng2 + 1):
            taufor = forfac[lay - 1] * (
                forref[_idx2(indf, ig, 4)]
                + forfrac[lay - 1]
                * (forref[_idx2(indf + 1, ig, 4)] - forref[_idx2(indf, ig, 4)])
            )
            taug[_idx2(lay, ngs1 + ig, nlayers)] = colh2o[lay - 1] * (
                fac00[lay - 1] * absb[_idx2(ind0, ig, 235)]
                + fac10[lay - 1] * absb[_idx2(ind0 + 1, ig, 235)]
                + fac01[lay - 1] * absb[_idx2(ind1, ig, 235)]
                + fac11[lay - 1] * absb[_idx2(ind1 + 1, ig, 235)]
            ) + taufor
            fracs[_idx2(lay, ngs1 + ig, nlayers)] = fracrefb[ig - 1]


@export
def rrtmg_lw_taugb6_codon(
    nlayers: int,
    laytrop: int,
    ng6: int,
    ngs5: int,
    nspa6: int,
    maxxsec: int,
    colh2o_p: cobj,
    colco2_p: cobj,
    coldry_p: cobj,
    wx_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    indself_p: cobj,
    indfor_p: cobj,
    indminor_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    selffac_p: cobj,
    selffrac_p: cobj,
    forfac_p: cobj,
    forfrac_p: cobj,
    minorfrac_p: cobj,
    chi_mls_p: cobj,
    fracrefa_p: cobj,
    absa_p: cobj,
    ka_mco2_p: cobj,
    selfref_p: cobj,
    forref_p: cobj,
    cfc11adj_p: cobj,
    cfc12_p: cobj,
    fracs_p: cobj,
    taug_p: cobj,
):
    colh2o = Ptr[float](colh2o_p)
    colco2 = Ptr[float](colco2_p)
    coldry = Ptr[float](coldry_p)
    wx = Ptr[float](wx_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    indself = Ptr[int](indself_p)
    indfor = Ptr[int](indfor_p)
    indminor = Ptr[int](indminor_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    selffac = Ptr[float](selffac_p)
    selffrac = Ptr[float](selffrac_p)
    forfac = Ptr[float](forfac_p)
    forfrac = Ptr[float](forfrac_p)
    minorfrac = Ptr[float](minorfrac_p)
    chi_mls = Ptr[float](chi_mls_p)
    fracrefa = Ptr[float](fracrefa_p)
    absa = Ptr[float](absa_p)
    ka_mco2 = Ptr[float](ka_mco2_p)
    selfref = Ptr[float](selfref_p)
    forref = Ptr[float](forref_p)
    cfc11adj = Ptr[float](cfc11adj_p)
    cfc12 = Ptr[float](cfc12_p)
    fracs = Ptr[float](fracs_p)
    taug = Ptr[float](taug_p)

    for lay in range(1, laytrop + 1):
        chi_co2 = colco2[lay - 1] / coldry[lay - 1]
        ratco2 = 1.0e20 * chi_co2 / chi_mls[_idx2(2, jp[lay - 1] + 1, 7)]
        if ratco2 > 3.0:
            adjfac = 2.0 + (ratco2 - 2.0) ** 0.77
            adjcolco2 = (
                adjfac
                * chi_mls[_idx2(2, jp[lay - 1] + 1, 7)]
                * coldry[lay - 1]
                * 1.0e-20
            )
        else:
            adjcolco2 = colco2[lay - 1]

        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa6 + 1
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa6 + 1
        inds = indself[lay - 1]
        indf = indfor[lay - 1]
        indm = indminor[lay - 1]

        for ig in range(1, ng6 + 1):
            tauself = selffac[lay - 1] * (
                selfref[_idx2(inds, ig, 10)]
                + selffrac[lay - 1]
                * (selfref[_idx2(inds + 1, ig, 10)] - selfref[_idx2(inds, ig, 10)])
            )
            taufor = forfac[lay - 1] * (
                forref[_idx2(indf, ig, 4)]
                + forfrac[lay - 1]
                * (forref[_idx2(indf + 1, ig, 4)] - forref[_idx2(indf, ig, 4)])
            )
            absco2 = ka_mco2[_idx2(indm, ig, 19)] + minorfrac[lay - 1] * (
                ka_mco2[_idx2(indm + 1, ig, 19)] - ka_mco2[_idx2(indm, ig, 19)]
            )
            taug[_idx2(lay, ngs5 + ig, nlayers)] = (
                colh2o[lay - 1]
                * (
                    fac00[lay - 1] * absa[_idx2(ind0, ig, 65)]
                    + fac10[lay - 1] * absa[_idx2(ind0 + 1, ig, 65)]
                    + fac01[lay - 1] * absa[_idx2(ind1, ig, 65)]
                    + fac11[lay - 1] * absa[_idx2(ind1 + 1, ig, 65)]
                )
                + tauself
                + taufor
                + adjcolco2 * absco2
                + wx[_idx2(2, lay, maxxsec)] * cfc11adj[ig - 1]
                + wx[_idx2(3, lay, maxxsec)] * cfc12[ig - 1]
            )
            fracs[_idx2(lay, ngs5 + ig, nlayers)] = fracrefa[ig - 1]

    for lay in range(laytrop + 1, nlayers + 1):
        for ig in range(1, ng6 + 1):
            taug[_idx2(lay, ngs5 + ig, nlayers)] = (
                0.0
                + wx[_idx2(2, lay, maxxsec)] * cfc11adj[ig - 1]
                + wx[_idx2(3, lay, maxxsec)] * cfc12[ig - 1]
            )
            fracs[_idx2(lay, ngs5 + ig, nlayers)] = fracrefa[ig - 1]


@export
def rrtmg_lw_taugb8_codon(
    nlayers: int,
    laytrop: int,
    ng8: int,
    ngs7: int,
    nspa8: int,
    nspb8: int,
    maxxsec: int,
    colh2o_p: cobj,
    colco2_p: cobj,
    colo3_p: cobj,
    coln2o_p: cobj,
    coldry_p: cobj,
    wx_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    indself_p: cobj,
    indfor_p: cobj,
    indminor_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    selffac_p: cobj,
    selffrac_p: cobj,
    forfac_p: cobj,
    forfrac_p: cobj,
    minorfrac_p: cobj,
    chi_mls_p: cobj,
    fracrefa_p: cobj,
    fracrefb_p: cobj,
    absa_p: cobj,
    absb_p: cobj,
    ka_mco2_p: cobj,
    ka_mn2o_p: cobj,
    ka_mo3_p: cobj,
    kb_mco2_p: cobj,
    kb_mn2o_p: cobj,
    selfref_p: cobj,
    forref_p: cobj,
    cfc12_p: cobj,
    cfc22adj_p: cobj,
    fracs_p: cobj,
    taug_p: cobj,
):
    colh2o = Ptr[float](colh2o_p)
    colco2 = Ptr[float](colco2_p)
    colo3 = Ptr[float](colo3_p)
    coln2o = Ptr[float](coln2o_p)
    coldry = Ptr[float](coldry_p)
    wx = Ptr[float](wx_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    indself = Ptr[int](indself_p)
    indfor = Ptr[int](indfor_p)
    indminor = Ptr[int](indminor_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    selffac = Ptr[float](selffac_p)
    selffrac = Ptr[float](selffrac_p)
    forfac = Ptr[float](forfac_p)
    forfrac = Ptr[float](forfrac_p)
    minorfrac = Ptr[float](minorfrac_p)
    chi_mls = Ptr[float](chi_mls_p)
    fracrefa = Ptr[float](fracrefa_p)
    fracrefb = Ptr[float](fracrefb_p)
    absa = Ptr[float](absa_p)
    absb = Ptr[float](absb_p)
    ka_mco2 = Ptr[float](ka_mco2_p)
    ka_mn2o = Ptr[float](ka_mn2o_p)
    ka_mo3 = Ptr[float](ka_mo3_p)
    kb_mco2 = Ptr[float](kb_mco2_p)
    kb_mn2o = Ptr[float](kb_mn2o_p)
    selfref = Ptr[float](selfref_p)
    forref = Ptr[float](forref_p)
    cfc12 = Ptr[float](cfc12_p)
    cfc22adj = Ptr[float](cfc22adj_p)
    fracs = Ptr[float](fracs_p)
    taug = Ptr[float](taug_p)

    for lay in range(1, laytrop + 1):
        chi_co2 = colco2[lay - 1] / coldry[lay - 1]
        ratco2 = 1.0e20 * chi_co2 / chi_mls[_idx2(2, jp[lay - 1] + 1, 7)]
        if ratco2 > 3.0:
            adjfac = 2.0 + (ratco2 - 2.0) ** 0.65
            adjcolco2 = (
                adjfac
                * chi_mls[_idx2(2, jp[lay - 1] + 1, 7)]
                * coldry[lay - 1]
                * 1.0e-20
            )
        else:
            adjcolco2 = colco2[lay - 1]

        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa8 + 1
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa8 + 1
        inds = indself[lay - 1]
        indf = indfor[lay - 1]
        indm = indminor[lay - 1]

        for ig in range(1, ng8 + 1):
            tauself = selffac[lay - 1] * (
                selfref[_idx2(inds, ig, 10)]
                + selffrac[lay - 1]
                * (selfref[_idx2(inds + 1, ig, 10)] - selfref[_idx2(inds, ig, 10)])
            )
            taufor = forfac[lay - 1] * (
                forref[_idx2(indf, ig, 4)]
                + forfrac[lay - 1]
                * (forref[_idx2(indf + 1, ig, 4)] - forref[_idx2(indf, ig, 4)])
            )
            absco2 = ka_mco2[_idx2(indm, ig, 19)] + minorfrac[lay - 1] * (
                ka_mco2[_idx2(indm + 1, ig, 19)] - ka_mco2[_idx2(indm, ig, 19)]
            )
            abso3 = ka_mo3[_idx2(indm, ig, 19)] + minorfrac[lay - 1] * (
                ka_mo3[_idx2(indm + 1, ig, 19)] - ka_mo3[_idx2(indm, ig, 19)]
            )
            absn2o = ka_mn2o[_idx2(indm, ig, 19)] + minorfrac[lay - 1] * (
                ka_mn2o[_idx2(indm + 1, ig, 19)] - ka_mn2o[_idx2(indm, ig, 19)]
            )
            taug[_idx2(lay, ngs7 + ig, nlayers)] = (
                colh2o[lay - 1]
                * (
                    fac00[lay - 1] * absa[_idx2(ind0, ig, 65)]
                    + fac10[lay - 1] * absa[_idx2(ind0 + 1, ig, 65)]
                    + fac01[lay - 1] * absa[_idx2(ind1, ig, 65)]
                    + fac11[lay - 1] * absa[_idx2(ind1 + 1, ig, 65)]
                )
                + tauself
                + taufor
                + adjcolco2 * absco2
                + colo3[lay - 1] * abso3
                + coln2o[lay - 1] * absn2o
                + wx[_idx2(3, lay, maxxsec)] * cfc12[ig - 1]
                + wx[_idx2(4, lay, maxxsec)] * cfc22adj[ig - 1]
            )
            fracs[_idx2(lay, ngs7 + ig, nlayers)] = fracrefa[ig - 1]

    for lay in range(laytrop + 1, nlayers + 1):
        chi_co2 = colco2[lay - 1] / coldry[lay - 1]
        ratco2 = 1.0e20 * chi_co2 / chi_mls[_idx2(2, jp[lay - 1] + 1, 7)]
        if ratco2 > 3.0:
            adjfac = 2.0 + (ratco2 - 2.0) ** 0.65
            adjcolco2 = (
                adjfac
                * chi_mls[_idx2(2, jp[lay - 1] + 1, 7)]
                * coldry[lay - 1]
                * 1.0e-20
            )
        else:
            adjcolco2 = colco2[lay - 1]

        ind0 = ((jp[lay - 1] - 13) * 5 + (jt[lay - 1] - 1)) * nspb8 + 1
        ind1 = ((jp[lay - 1] - 12) * 5 + (jt1[lay - 1] - 1)) * nspb8 + 1
        indm = indminor[lay - 1]

        for ig in range(1, ng8 + 1):
            absco2 = kb_mco2[_idx2(indm, ig, 19)] + minorfrac[lay - 1] * (
                kb_mco2[_idx2(indm + 1, ig, 19)] - kb_mco2[_idx2(indm, ig, 19)]
            )
            absn2o = kb_mn2o[_idx2(indm, ig, 19)] + minorfrac[lay - 1] * (
                kb_mn2o[_idx2(indm + 1, ig, 19)] - kb_mn2o[_idx2(indm, ig, 19)]
            )
            taug[_idx2(lay, ngs7 + ig, nlayers)] = (
                colo3[lay - 1]
                * (
                    fac00[lay - 1] * absb[_idx2(ind0, ig, 235)]
                    + fac10[lay - 1] * absb[_idx2(ind0 + 1, ig, 235)]
                    + fac01[lay - 1] * absb[_idx2(ind1, ig, 235)]
                    + fac11[lay - 1] * absb[_idx2(ind1 + 1, ig, 235)]
                )
                + adjcolco2 * absco2
                + coln2o[lay - 1] * absn2o
                + wx[_idx2(3, lay, maxxsec)] * cfc12[ig - 1]
                + wx[_idx2(4, lay, maxxsec)] * cfc22adj[ig - 1]
            )
            fracs[_idx2(lay, ngs7 + ig, nlayers)] = fracrefb[ig - 1]


@export
def rrtmg_lw_taugb12_codon(
    nlayers: int,
    laytrop: int,
    ng12: int,
    ngs11: int,
    nspa12: int,
    oneminus: float,
    colh2o_p: cobj,
    colco2_p: cobj,
    rat_h2oco2_p: cobj,
    rat_h2oco2_1_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    indself_p: cobj,
    indfor_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    selffac_p: cobj,
    selffrac_p: cobj,
    forfac_p: cobj,
    forfrac_p: cobj,
    chi_mls_p: cobj,
    fracrefa_p: cobj,
    absa_p: cobj,
    selfref_p: cobj,
    forref_p: cobj,
    fracs_p: cobj,
    taug_p: cobj,
):
    colh2o = Ptr[float](colh2o_p)
    colco2 = Ptr[float](colco2_p)
    rat_h2oco2 = Ptr[float](rat_h2oco2_p)
    rat_h2oco2_1 = Ptr[float](rat_h2oco2_1_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    indself = Ptr[int](indself_p)
    indfor = Ptr[int](indfor_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    selffac = Ptr[float](selffac_p)
    selffrac = Ptr[float](selffrac_p)
    forfac = Ptr[float](forfac_p)
    forfrac = Ptr[float](forfrac_p)
    chi_mls = Ptr[float](chi_mls_p)
    fracrefa = Ptr[float](fracrefa_p)
    absa = Ptr[float](absa_p)
    selfref = Ptr[float](selfref_p)
    forref = Ptr[float](forref_p)
    fracs = Ptr[float](fracs_p)
    taug = Ptr[float](taug_p)

    refrat_planck_a = chi_mls[_idx2(1, 10, 7)] / chi_mls[_idx2(2, 10, 7)]

    for lay in range(1, laytrop + 1):
        speccomb = colh2o[lay - 1] + rat_h2oco2[lay - 1] * colco2[lay - 1]
        specparm = colh2o[lay - 1] / speccomb
        if specparm >= oneminus:
            specparm = oneminus
        specmult = 8.0 * specparm
        js = 1 + int(specmult)
        fs = specmult - float(int(specmult))

        speccomb1 = colh2o[lay - 1] + rat_h2oco2_1[lay - 1] * colco2[lay - 1]
        specparm1 = colh2o[lay - 1] / speccomb1
        if specparm1 >= oneminus:
            specparm1 = oneminus
        specmult1 = 8.0 * specparm1
        js1 = 1 + int(specmult1)
        fs1 = specmult1 - float(int(specmult1))

        speccomb_planck = colh2o[lay - 1] + refrat_planck_a * colco2[lay - 1]
        specparm_planck = colh2o[lay - 1] / speccomb_planck
        if specparm_planck >= oneminus:
            specparm_planck = oneminus
        specmult_planck = 8.0 * specparm_planck
        jpl = 1 + int(specmult_planck)
        fpl = specmult_planck - float(int(specmult_planck))

        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa12 + js
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa12 + js1
        inds = indself[lay - 1]
        indf = indfor[lay - 1]

        if specparm < 0.125:
            p = fs - 1.0
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac000 = fk0 * fac00[lay - 1]
            fac100 = fk1 * fac00[lay - 1]
            fac200 = fk2 * fac00[lay - 1]
            fac010 = fk0 * fac10[lay - 1]
            fac110 = fk1 * fac10[lay - 1]
            fac210 = fk2 * fac10[lay - 1]
        elif specparm > 0.875:
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac000 = fk0 * fac00[lay - 1]
            fac100 = fk1 * fac00[lay - 1]
            fac200 = fk2 * fac00[lay - 1]
            fac010 = fk0 * fac10[lay - 1]
            fac110 = fk1 * fac10[lay - 1]
            fac210 = fk2 * fac10[lay - 1]
        else:
            fac000 = (1.0 - fs) * fac00[lay - 1]
            fac010 = (1.0 - fs) * fac10[lay - 1]
            fac100 = fs * fac00[lay - 1]
            fac110 = fs * fac10[lay - 1]

        if specparm1 < 0.125:
            p = fs1 - 1.0
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac001 = fk0 * fac01[lay - 1]
            fac101 = fk1 * fac01[lay - 1]
            fac201 = fk2 * fac01[lay - 1]
            fac011 = fk0 * fac11[lay - 1]
            fac111 = fk1 * fac11[lay - 1]
            fac211 = fk2 * fac11[lay - 1]
        elif specparm1 > 0.875:
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac001 = fk0 * fac01[lay - 1]
            fac101 = fk1 * fac01[lay - 1]
            fac201 = fk2 * fac01[lay - 1]
            fac011 = fk0 * fac11[lay - 1]
            fac111 = fk1 * fac11[lay - 1]
            fac211 = fk2 * fac11[lay - 1]
        else:
            fac001 = (1.0 - fs1) * fac01[lay - 1]
            fac011 = (1.0 - fs1) * fac11[lay - 1]
            fac101 = fs1 * fac01[lay - 1]
            fac111 = fs1 * fac11[lay - 1]

        for ig in range(1, ng12 + 1):
            tauself = selffac[lay - 1] * (
                selfref[_idx2(inds, ig, 10)]
                + selffrac[lay - 1]
                * (selfref[_idx2(inds + 1, ig, 10)] - selfref[_idx2(inds, ig, 10)])
            )
            taufor = forfac[lay - 1] * (
                forref[_idx2(indf, ig, 4)]
                + forfrac[lay - 1]
                * (forref[_idx2(indf + 1, ig, 4)] - forref[_idx2(indf, ig, 4)])
            )

            if specparm < 0.125:
                tau_major = speccomb * (
                    fac000 * absa[_idx2(ind0, ig, 585)]
                    + fac100 * absa[_idx2(ind0 + 1, ig, 585)]
                    + fac200 * absa[_idx2(ind0 + 2, ig, 585)]
                    + fac010 * absa[_idx2(ind0 + 9, ig, 585)]
                    + fac110 * absa[_idx2(ind0 + 10, ig, 585)]
                    + fac210 * absa[_idx2(ind0 + 11, ig, 585)]
                )
            elif specparm > 0.875:
                tau_major = speccomb * (
                    fac200 * absa[_idx2(ind0 - 1, ig, 585)]
                    + fac100 * absa[_idx2(ind0, ig, 585)]
                    + fac000 * absa[_idx2(ind0 + 1, ig, 585)]
                    + fac210 * absa[_idx2(ind0 + 8, ig, 585)]
                    + fac110 * absa[_idx2(ind0 + 9, ig, 585)]
                    + fac010 * absa[_idx2(ind0 + 10, ig, 585)]
                )
            else:
                tau_major = speccomb * (
                    fac000 * absa[_idx2(ind0, ig, 585)]
                    + fac100 * absa[_idx2(ind0 + 1, ig, 585)]
                    + fac010 * absa[_idx2(ind0 + 9, ig, 585)]
                    + fac110 * absa[_idx2(ind0 + 10, ig, 585)]
                )

            if specparm1 < 0.125:
                tau_major1 = speccomb1 * (
                    fac001 * absa[_idx2(ind1, ig, 585)]
                    + fac101 * absa[_idx2(ind1 + 1, ig, 585)]
                    + fac201 * absa[_idx2(ind1 + 2, ig, 585)]
                    + fac011 * absa[_idx2(ind1 + 9, ig, 585)]
                    + fac111 * absa[_idx2(ind1 + 10, ig, 585)]
                    + fac211 * absa[_idx2(ind1 + 11, ig, 585)]
                )
            elif specparm1 > 0.875:
                tau_major1 = speccomb1 * (
                    fac201 * absa[_idx2(ind1 - 1, ig, 585)]
                    + fac101 * absa[_idx2(ind1, ig, 585)]
                    + fac001 * absa[_idx2(ind1 + 1, ig, 585)]
                    + fac211 * absa[_idx2(ind1 + 8, ig, 585)]
                    + fac111 * absa[_idx2(ind1 + 9, ig, 585)]
                    + fac011 * absa[_idx2(ind1 + 10, ig, 585)]
                )
            else:
                tau_major1 = speccomb1 * (
                    fac001 * absa[_idx2(ind1, ig, 585)]
                    + fac101 * absa[_idx2(ind1 + 1, ig, 585)]
                    + fac011 * absa[_idx2(ind1 + 9, ig, 585)]
                    + fac111 * absa[_idx2(ind1 + 10, ig, 585)]
                )

            taug[_idx2(lay, ngs11 + ig, nlayers)] = tau_major + tau_major1 + tauself + taufor
            fracs[_idx2(lay, ngs11 + ig, nlayers)] = fracrefa[
                _idx2(ig, jpl, ng12)
            ] + fpl * (
                fracrefa[_idx2(ig, jpl + 1, ng12)] - fracrefa[_idx2(ig, jpl, ng12)]
            )

    for lay in range(laytrop + 1, nlayers + 1):
        for ig in range(1, ng12 + 1):
            taug[_idx2(lay, ngs11 + ig, nlayers)] = 0.0
            fracs[_idx2(lay, ngs11 + ig, nlayers)] = 0.0


@export
def rrtmg_lw_taugb16_codon(
    nlayers: int,
    laytrop: int,
    ng16: int,
    ngs15: int,
    nspa16: int,
    nspb16: int,
    oneminus: float,
    colh2o_p: cobj,
    colch4_p: cobj,
    rat_h2och4_p: cobj,
    rat_h2och4_1_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    indself_p: cobj,
    indfor_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    selffac_p: cobj,
    selffrac_p: cobj,
    forfac_p: cobj,
    forfrac_p: cobj,
    chi_mls_p: cobj,
    fracrefa_p: cobj,
    fracrefb_p: cobj,
    absa_p: cobj,
    absb_p: cobj,
    selfref_p: cobj,
    forref_p: cobj,
    fracs_p: cobj,
    taug_p: cobj,
):
    colh2o = Ptr[float](colh2o_p)
    colch4 = Ptr[float](colch4_p)
    rat_h2och4 = Ptr[float](rat_h2och4_p)
    rat_h2och4_1 = Ptr[float](rat_h2och4_1_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    indself = Ptr[int](indself_p)
    indfor = Ptr[int](indfor_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    selffac = Ptr[float](selffac_p)
    selffrac = Ptr[float](selffrac_p)
    forfac = Ptr[float](forfac_p)
    forfrac = Ptr[float](forfrac_p)
    chi_mls = Ptr[float](chi_mls_p)
    fracrefa = Ptr[float](fracrefa_p)
    fracrefb = Ptr[float](fracrefb_p)
    absa = Ptr[float](absa_p)
    absb = Ptr[float](absb_p)
    selfref = Ptr[float](selfref_p)
    forref = Ptr[float](forref_p)
    fracs = Ptr[float](fracs_p)
    taug = Ptr[float](taug_p)

    refrat_planck_a = chi_mls[_idx2(1, 6, 7)] / chi_mls[_idx2(6, 6, 7)]

    for lay in range(1, laytrop + 1):
        speccomb = colh2o[lay - 1] + rat_h2och4[lay - 1] * colch4[lay - 1]
        specparm = colh2o[lay - 1] / speccomb
        if specparm >= oneminus:
            specparm = oneminus
        specmult = 8.0 * specparm
        js = 1 + int(specmult)
        fs = specmult - float(int(specmult))

        speccomb1 = colh2o[lay - 1] + rat_h2och4_1[lay - 1] * colch4[lay - 1]
        specparm1 = colh2o[lay - 1] / speccomb1
        if specparm1 >= oneminus:
            specparm1 = oneminus
        specmult1 = 8.0 * specparm1
        js1 = 1 + int(specmult1)
        fs1 = specmult1 - float(int(specmult1))

        speccomb_planck = colh2o[lay - 1] + refrat_planck_a * colch4[lay - 1]
        specparm_planck = colh2o[lay - 1] / speccomb_planck
        if specparm_planck >= oneminus:
            specparm_planck = oneminus
        specmult_planck = 8.0 * specparm_planck
        jpl = 1 + int(specmult_planck)
        fpl = specmult_planck - float(int(specmult_planck))

        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa16 + js
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa16 + js1
        inds = indself[lay - 1]
        indf = indfor[lay - 1]

        if specparm < 0.125:
            p = fs - 1.0
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac000 = fk0 * fac00[lay - 1]
            fac100 = fk1 * fac00[lay - 1]
            fac200 = fk2 * fac00[lay - 1]
            fac010 = fk0 * fac10[lay - 1]
            fac110 = fk1 * fac10[lay - 1]
            fac210 = fk2 * fac10[lay - 1]
        elif specparm > 0.875:
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac000 = fk0 * fac00[lay - 1]
            fac100 = fk1 * fac00[lay - 1]
            fac200 = fk2 * fac00[lay - 1]
            fac010 = fk0 * fac10[lay - 1]
            fac110 = fk1 * fac10[lay - 1]
            fac210 = fk2 * fac10[lay - 1]
        else:
            fac000 = (1.0 - fs) * fac00[lay - 1]
            fac010 = (1.0 - fs) * fac10[lay - 1]
            fac100 = fs * fac00[lay - 1]
            fac110 = fs * fac10[lay - 1]

        if specparm1 < 0.125:
            p = fs1 - 1.0
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac001 = fk0 * fac01[lay - 1]
            fac101 = fk1 * fac01[lay - 1]
            fac201 = fk2 * fac01[lay - 1]
            fac011 = fk0 * fac11[lay - 1]
            fac111 = fk1 * fac11[lay - 1]
            fac211 = fk2 * fac11[lay - 1]
        elif specparm1 > 0.875:
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac001 = fk0 * fac01[lay - 1]
            fac101 = fk1 * fac01[lay - 1]
            fac201 = fk2 * fac01[lay - 1]
            fac011 = fk0 * fac11[lay - 1]
            fac111 = fk1 * fac11[lay - 1]
            fac211 = fk2 * fac11[lay - 1]
        else:
            fac001 = (1.0 - fs1) * fac01[lay - 1]
            fac011 = (1.0 - fs1) * fac11[lay - 1]
            fac101 = fs1 * fac01[lay - 1]
            fac111 = fs1 * fac11[lay - 1]

        for ig in range(1, ng16 + 1):
            tauself = selffac[lay - 1] * (
                selfref[_idx2(inds, ig, 10)]
                + selffrac[lay - 1]
                * (selfref[_idx2(inds + 1, ig, 10)] - selfref[_idx2(inds, ig, 10)])
            )
            taufor = forfac[lay - 1] * (
                forref[_idx2(indf, ig, 4)]
                + forfrac[lay - 1]
                * (forref[_idx2(indf + 1, ig, 4)] - forref[_idx2(indf, ig, 4)])
            )

            if specparm < 0.125:
                tau_major = speccomb * (
                    fac000 * absa[_idx2(ind0, ig, 585)]
                    + fac100 * absa[_idx2(ind0 + 1, ig, 585)]
                    + fac200 * absa[_idx2(ind0 + 2, ig, 585)]
                    + fac010 * absa[_idx2(ind0 + 9, ig, 585)]
                    + fac110 * absa[_idx2(ind0 + 10, ig, 585)]
                    + fac210 * absa[_idx2(ind0 + 11, ig, 585)]
                )
            elif specparm > 0.875:
                tau_major = speccomb * (
                    fac200 * absa[_idx2(ind0 - 1, ig, 585)]
                    + fac100 * absa[_idx2(ind0, ig, 585)]
                    + fac000 * absa[_idx2(ind0 + 1, ig, 585)]
                    + fac210 * absa[_idx2(ind0 + 8, ig, 585)]
                    + fac110 * absa[_idx2(ind0 + 9, ig, 585)]
                    + fac010 * absa[_idx2(ind0 + 10, ig, 585)]
                )
            else:
                tau_major = speccomb * (
                    fac000 * absa[_idx2(ind0, ig, 585)]
                    + fac100 * absa[_idx2(ind0 + 1, ig, 585)]
                    + fac010 * absa[_idx2(ind0 + 9, ig, 585)]
                    + fac110 * absa[_idx2(ind0 + 10, ig, 585)]
                )

            if specparm1 < 0.125:
                tau_major1 = speccomb1 * (
                    fac001 * absa[_idx2(ind1, ig, 585)]
                    + fac101 * absa[_idx2(ind1 + 1, ig, 585)]
                    + fac201 * absa[_idx2(ind1 + 2, ig, 585)]
                    + fac011 * absa[_idx2(ind1 + 9, ig, 585)]
                    + fac111 * absa[_idx2(ind1 + 10, ig, 585)]
                    + fac211 * absa[_idx2(ind1 + 11, ig, 585)]
                )
            elif specparm1 > 0.875:
                tau_major1 = speccomb1 * (
                    fac201 * absa[_idx2(ind1 - 1, ig, 585)]
                    + fac101 * absa[_idx2(ind1, ig, 585)]
                    + fac001 * absa[_idx2(ind1 + 1, ig, 585)]
                    + fac211 * absa[_idx2(ind1 + 8, ig, 585)]
                    + fac111 * absa[_idx2(ind1 + 9, ig, 585)]
                    + fac011 * absa[_idx2(ind1 + 10, ig, 585)]
                )
            else:
                tau_major1 = speccomb1 * (
                    fac001 * absa[_idx2(ind1, ig, 585)]
                    + fac101 * absa[_idx2(ind1 + 1, ig, 585)]
                    + fac011 * absa[_idx2(ind1 + 9, ig, 585)]
                    + fac111 * absa[_idx2(ind1 + 10, ig, 585)]
                )

            taug[_idx2(lay, ngs15 + ig, nlayers)] = tau_major + tau_major1 + tauself + taufor
            fracs[_idx2(lay, ngs15 + ig, nlayers)] = fracrefa[
                _idx2(ig, jpl, ng16)
            ] + fpl * (
                fracrefa[_idx2(ig, jpl + 1, ng16)] - fracrefa[_idx2(ig, jpl, ng16)]
            )

    for lay in range(laytrop + 1, nlayers + 1):
        ind0 = ((jp[lay - 1] - 13) * 5 + (jt[lay - 1] - 1)) * nspb16 + 1
        ind1 = ((jp[lay - 1] - 12) * 5 + (jt1[lay - 1] - 1)) * nspb16 + 1
        for ig in range(1, ng16 + 1):
            taug[_idx2(lay, ngs15 + ig, nlayers)] = colch4[lay - 1] * (
                fac00[lay - 1] * absb[_idx2(ind0, ig, 235)]
                + fac10[lay - 1] * absb[_idx2(ind0 + 1, ig, 235)]
                + fac01[lay - 1] * absb[_idx2(ind1, ig, 235)]
                + fac11[lay - 1] * absb[_idx2(ind1 + 1, ig, 235)]
            )
            fracs[_idx2(lay, ngs15 + ig, nlayers)] = fracrefb[ig - 1]


@export
def rrtmg_lw_taugb4_codon(
    nlayers: int,
    laytrop: int,
    ng4: int,
    ngs3: int,
    nspa4: int,
    nspb4: int,
    oneminus: float,
    colh2o_p: cobj,
    colco2_p: cobj,
    colo3_p: cobj,
    rat_h2oco2_p: cobj,
    rat_h2oco2_1_p: cobj,
    rat_o3co2_p: cobj,
    rat_o3co2_1_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    indself_p: cobj,
    indfor_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    selffac_p: cobj,
    selffrac_p: cobj,
    forfac_p: cobj,
    forfrac_p: cobj,
    chi_mls_p: cobj,
    fracrefa_p: cobj,
    fracrefb_p: cobj,
    absa_p: cobj,
    absb_p: cobj,
    selfref_p: cobj,
    forref_p: cobj,
    fracs_p: cobj,
    taug_p: cobj,
):
    colh2o = Ptr[float](colh2o_p)
    colco2 = Ptr[float](colco2_p)
    colo3 = Ptr[float](colo3_p)
    rat_h2oco2 = Ptr[float](rat_h2oco2_p)
    rat_h2oco2_1 = Ptr[float](rat_h2oco2_1_p)
    rat_o3co2 = Ptr[float](rat_o3co2_p)
    rat_o3co2_1 = Ptr[float](rat_o3co2_1_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    indself = Ptr[int](indself_p)
    indfor = Ptr[int](indfor_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    selffac = Ptr[float](selffac_p)
    selffrac = Ptr[float](selffrac_p)
    forfac = Ptr[float](forfac_p)
    forfrac = Ptr[float](forfrac_p)
    chi_mls = Ptr[float](chi_mls_p)
    fracrefa = Ptr[float](fracrefa_p)
    fracrefb = Ptr[float](fracrefb_p)
    absa = Ptr[float](absa_p)
    absb = Ptr[float](absb_p)
    selfref = Ptr[float](selfref_p)
    forref = Ptr[float](forref_p)
    fracs = Ptr[float](fracs_p)
    taug = Ptr[float](taug_p)

    refrat_planck_a = chi_mls[_idx2(1, 11, 7)] / chi_mls[_idx2(2, 11, 7)]
    refrat_planck_b = chi_mls[_idx2(3, 13, 7)] / chi_mls[_idx2(2, 13, 7)]

    for lay in range(1, laytrop + 1):
        speccomb = colh2o[lay - 1] + rat_h2oco2[lay - 1] * colco2[lay - 1]
        specparm = colh2o[lay - 1] / speccomb
        if specparm >= oneminus:
            specparm = oneminus
        specmult = 8.0 * specparm
        js = 1 + int(specmult)
        fs = specmult - float(int(specmult))

        speccomb1 = colh2o[lay - 1] + rat_h2oco2_1[lay - 1] * colco2[lay - 1]
        specparm1 = colh2o[lay - 1] / speccomb1
        if specparm1 >= oneminus:
            specparm1 = oneminus
        specmult1 = 8.0 * specparm1
        js1 = 1 + int(specmult1)
        fs1 = specmult1 - float(int(specmult1))

        speccomb_planck = colh2o[lay - 1] + refrat_planck_a * colco2[lay - 1]
        specparm_planck = colh2o[lay - 1] / speccomb_planck
        if specparm_planck >= oneminus:
            specparm_planck = oneminus
        specmult_planck = 8.0 * specparm_planck
        jpl = 1 + int(specmult_planck)
        fpl = specmult_planck - float(int(specmult_planck))

        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa4 + js
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa4 + js1
        inds = indself[lay - 1]
        indf = indfor[lay - 1]
        fac200 = 0.0
        fac210 = 0.0
        fac201 = 0.0
        fac211 = 0.0

        if specparm < 0.125:
            p = fs - 1.0
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac000 = fk0 * fac00[lay - 1]
            fac100 = fk1 * fac00[lay - 1]
            fac200 = fk2 * fac00[lay - 1]
            fac010 = fk0 * fac10[lay - 1]
            fac110 = fk1 * fac10[lay - 1]
            fac210 = fk2 * fac10[lay - 1]
        elif specparm > 0.875:
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac000 = fk0 * fac00[lay - 1]
            fac100 = fk1 * fac00[lay - 1]
            fac200 = fk2 * fac00[lay - 1]
            fac010 = fk0 * fac10[lay - 1]
            fac110 = fk1 * fac10[lay - 1]
            fac210 = fk2 * fac10[lay - 1]
        else:
            fac000 = (1.0 - fs) * fac00[lay - 1]
            fac010 = (1.0 - fs) * fac10[lay - 1]
            fac100 = fs * fac00[lay - 1]
            fac110 = fs * fac10[lay - 1]

        if specparm1 < 0.125:
            p = fs1 - 1.0
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac001 = fk0 * fac01[lay - 1]
            fac101 = fk1 * fac01[lay - 1]
            fac201 = fk2 * fac01[lay - 1]
            fac011 = fk0 * fac11[lay - 1]
            fac111 = fk1 * fac11[lay - 1]
            fac211 = fk2 * fac11[lay - 1]
        elif specparm1 > 0.875:
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac001 = fk0 * fac01[lay - 1]
            fac101 = fk1 * fac01[lay - 1]
            fac201 = fk2 * fac01[lay - 1]
            fac011 = fk0 * fac11[lay - 1]
            fac111 = fk1 * fac11[lay - 1]
            fac211 = fk2 * fac11[lay - 1]
        else:
            fac001 = (1.0 - fs1) * fac01[lay - 1]
            fac011 = (1.0 - fs1) * fac11[lay - 1]
            fac101 = fs1 * fac01[lay - 1]
            fac111 = fs1 * fac11[lay - 1]

        for ig in range(1, ng4 + 1):
            tauself = selffac[lay - 1] * (
                selfref[_idx2(inds, ig, 10)]
                + selffrac[lay - 1]
                * (selfref[_idx2(inds + 1, ig, 10)] - selfref[_idx2(inds, ig, 10)])
            )
            taufor = forfac[lay - 1] * (
                forref[_idx2(indf, ig, 4)]
                + forfrac[lay - 1]
                * (forref[_idx2(indf + 1, ig, 4)] - forref[_idx2(indf, ig, 4)])
            )
            tau_major = _lw_tau_major_3pt(
                specparm,
                speccomb,
                ind0,
                ig,
                fac000,
                fac100,
                fac200,
                fac010,
                fac110,
                fac210,
                absa,
                585,
            )
            tau_major1 = _lw_tau_major_3pt(
                specparm1,
                speccomb1,
                ind1,
                ig,
                fac001,
                fac101,
                fac201,
                fac011,
                fac111,
                fac211,
                absa,
                585,
            )
            taug[_idx2(lay, ngs3 + ig, nlayers)] = tau_major + tau_major1 + tauself + taufor
            fracs[_idx2(lay, ngs3 + ig, nlayers)] = fracrefa[
                _idx2(ig, jpl, ng4)
            ] + fpl * (
                fracrefa[_idx2(ig, jpl + 1, ng4)] - fracrefa[_idx2(ig, jpl, ng4)]
            )

    for lay in range(laytrop + 1, nlayers + 1):
        speccomb = colo3[lay - 1] + rat_o3co2[lay - 1] * colco2[lay - 1]
        specparm = colo3[lay - 1] / speccomb
        if specparm >= oneminus:
            specparm = oneminus
        specmult = 4.0 * specparm
        js = 1 + int(specmult)
        fs = specmult - float(int(specmult))

        speccomb1 = colo3[lay - 1] + rat_o3co2_1[lay - 1] * colco2[lay - 1]
        specparm1 = colo3[lay - 1] / speccomb1
        if specparm1 >= oneminus:
            specparm1 = oneminus
        specmult1 = 4.0 * specparm1
        js1 = 1 + int(specmult1)
        fs1 = specmult1 - float(int(specmult1))

        fac000 = (1.0 - fs) * fac00[lay - 1]
        fac010 = (1.0 - fs) * fac10[lay - 1]
        fac100 = fs * fac00[lay - 1]
        fac110 = fs * fac10[lay - 1]
        fac001 = (1.0 - fs1) * fac01[lay - 1]
        fac011 = (1.0 - fs1) * fac11[lay - 1]
        fac101 = fs1 * fac01[lay - 1]
        fac111 = fs1 * fac11[lay - 1]

        speccomb_planck = colo3[lay - 1] + refrat_planck_b * colco2[lay - 1]
        specparm_planck = colo3[lay - 1] / speccomb_planck
        if specparm_planck >= oneminus:
            specparm_planck = oneminus
        specmult_planck = 4.0 * specparm_planck
        jpl = 1 + int(specmult_planck)
        fpl = specmult_planck - float(int(specmult_planck))

        ind0 = ((jp[lay - 1] - 13) * 5 + (jt[lay - 1] - 1)) * nspb4 + js
        ind1 = ((jp[lay - 1] - 12) * 5 + (jt1[lay - 1] - 1)) * nspb4 + js1

        for ig in range(1, ng4 + 1):
            taug[_idx2(lay, ngs3 + ig, nlayers)] = speccomb * (
                fac000 * absb[_idx2(ind0, ig, 1175)]
                + fac100 * absb[_idx2(ind0 + 1, ig, 1175)]
                + fac010 * absb[_idx2(ind0 + 5, ig, 1175)]
                + fac110 * absb[_idx2(ind0 + 6, ig, 1175)]
            ) + speccomb1 * (
                fac001 * absb[_idx2(ind1, ig, 1175)]
                + fac101 * absb[_idx2(ind1 + 1, ig, 1175)]
                + fac011 * absb[_idx2(ind1 + 5, ig, 1175)]
                + fac111 * absb[_idx2(ind1 + 6, ig, 1175)]
            )
            fracs[_idx2(lay, ngs3 + ig, nlayers)] = fracrefb[
                _idx2(ig, jpl, ng4)
            ] + fpl * (
                fracrefb[_idx2(ig, jpl + 1, ng4)] - fracrefb[_idx2(ig, jpl, ng4)]
            )

        taug[_idx2(lay, ngs3 + 8, nlayers)] = (
            taug[_idx2(lay, ngs3 + 8, nlayers)] * 0.9200000166893005
        )
        taug[_idx2(lay, ngs3 + 9, nlayers)] = (
            taug[_idx2(lay, ngs3 + 9, nlayers)] * 0.8799999952316284
        )
        taug[_idx2(lay, ngs3 + 10, nlayers)] = (
            taug[_idx2(lay, ngs3 + 10, nlayers)] * 1.0700000524520874
        )
        taug[_idx2(lay, ngs3 + 11, nlayers)] = (
            taug[_idx2(lay, ngs3 + 11, nlayers)] * 1.100000023841858
        )
        taug[_idx2(lay, ngs3 + 12, nlayers)] = (
            taug[_idx2(lay, ngs3 + 12, nlayers)] * 0.9900000095367432
        )
        taug[_idx2(lay, ngs3 + 13, nlayers)] = (
            taug[_idx2(lay, ngs3 + 13, nlayers)] * 0.8799999952316284
        )
        taug[_idx2(lay, ngs3 + 14, nlayers)] = (
            taug[_idx2(lay, ngs3 + 14, nlayers)] * 0.9430000185966492
        )


@export
def rrtmg_lw_taugb5_codon(
    nlayers: int,
    laytrop: int,
    ng5: int,
    ngs4: int,
    nspa5: int,
    nspb5: int,
    maxxsec: int,
    oneminus: float,
    colh2o_p: cobj,
    colco2_p: cobj,
    colo3_p: cobj,
    wx_p: cobj,
    rat_h2oco2_p: cobj,
    rat_h2oco2_1_p: cobj,
    rat_o3co2_p: cobj,
    rat_o3co2_1_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    indself_p: cobj,
    indfor_p: cobj,
    indminor_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    selffac_p: cobj,
    selffrac_p: cobj,
    forfac_p: cobj,
    forfrac_p: cobj,
    minorfrac_p: cobj,
    chi_mls_p: cobj,
    fracrefa_p: cobj,
    fracrefb_p: cobj,
    absa_p: cobj,
    absb_p: cobj,
    ka_mo3_p: cobj,
    selfref_p: cobj,
    forref_p: cobj,
    ccl4_p: cobj,
    fracs_p: cobj,
    taug_p: cobj,
):
    colh2o = Ptr[float](colh2o_p)
    colco2 = Ptr[float](colco2_p)
    colo3 = Ptr[float](colo3_p)
    wx = Ptr[float](wx_p)
    rat_h2oco2 = Ptr[float](rat_h2oco2_p)
    rat_h2oco2_1 = Ptr[float](rat_h2oco2_1_p)
    rat_o3co2 = Ptr[float](rat_o3co2_p)
    rat_o3co2_1 = Ptr[float](rat_o3co2_1_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    indself = Ptr[int](indself_p)
    indfor = Ptr[int](indfor_p)
    indminor = Ptr[int](indminor_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    selffac = Ptr[float](selffac_p)
    selffrac = Ptr[float](selffrac_p)
    forfac = Ptr[float](forfac_p)
    forfrac = Ptr[float](forfrac_p)
    minorfrac = Ptr[float](minorfrac_p)
    chi_mls = Ptr[float](chi_mls_p)
    fracrefa = Ptr[float](fracrefa_p)
    fracrefb = Ptr[float](fracrefb_p)
    absa = Ptr[float](absa_p)
    absb = Ptr[float](absb_p)
    ka_mo3 = Ptr[float](ka_mo3_p)
    selfref = Ptr[float](selfref_p)
    forref = Ptr[float](forref_p)
    ccl4 = Ptr[float](ccl4_p)
    fracs = Ptr[float](fracs_p)
    taug = Ptr[float](taug_p)

    refrat_planck_a = chi_mls[_idx2(1, 5, 7)] / chi_mls[_idx2(2, 5, 7)]
    refrat_planck_b = chi_mls[_idx2(3, 43, 7)] / chi_mls[_idx2(2, 43, 7)]
    refrat_m_a = chi_mls[_idx2(1, 7, 7)] / chi_mls[_idx2(2, 7, 7)]

    for lay in range(1, laytrop + 1):
        speccomb = colh2o[lay - 1] + rat_h2oco2[lay - 1] * colco2[lay - 1]
        specparm = colh2o[lay - 1] / speccomb
        if specparm >= oneminus:
            specparm = oneminus
        specmult = 8.0 * specparm
        js = 1 + int(specmult)
        fs = specmult - float(int(specmult))

        speccomb1 = colh2o[lay - 1] + rat_h2oco2_1[lay - 1] * colco2[lay - 1]
        specparm1 = colh2o[lay - 1] / speccomb1
        if specparm1 >= oneminus:
            specparm1 = oneminus
        specmult1 = 8.0 * specparm1
        js1 = 1 + int(specmult1)
        fs1 = specmult1 - float(int(specmult1))

        speccomb_mo3 = colh2o[lay - 1] + refrat_m_a * colco2[lay - 1]
        specparm_mo3 = colh2o[lay - 1] / speccomb_mo3
        if specparm_mo3 >= oneminus:
            specparm_mo3 = oneminus
        specmult_mo3 = 8.0 * specparm_mo3
        jmo3 = 1 + int(specmult_mo3)
        fmo3 = specmult_mo3 - float(int(specmult_mo3))

        speccomb_planck = colh2o[lay - 1] + refrat_planck_a * colco2[lay - 1]
        specparm_planck = colh2o[lay - 1] / speccomb_planck
        if specparm_planck >= oneminus:
            specparm_planck = oneminus
        specmult_planck = 8.0 * specparm_planck
        jpl = 1 + int(specmult_planck)
        fpl = specmult_planck - float(int(specmult_planck))

        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa5 + js
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa5 + js1
        inds = indself[lay - 1]
        indf = indfor[lay - 1]
        indm = indminor[lay - 1]
        fac200 = 0.0
        fac210 = 0.0
        fac201 = 0.0
        fac211 = 0.0

        if specparm < 0.125:
            p = fs - 1.0
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac000 = fk0 * fac00[lay - 1]
            fac100 = fk1 * fac00[lay - 1]
            fac200 = fk2 * fac00[lay - 1]
            fac010 = fk0 * fac10[lay - 1]
            fac110 = fk1 * fac10[lay - 1]
            fac210 = fk2 * fac10[lay - 1]
        elif specparm > 0.875:
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac000 = fk0 * fac00[lay - 1]
            fac100 = fk1 * fac00[lay - 1]
            fac200 = fk2 * fac00[lay - 1]
            fac010 = fk0 * fac10[lay - 1]
            fac110 = fk1 * fac10[lay - 1]
            fac210 = fk2 * fac10[lay - 1]
        else:
            fac000 = (1.0 - fs) * fac00[lay - 1]
            fac010 = (1.0 - fs) * fac10[lay - 1]
            fac100 = fs * fac00[lay - 1]
            fac110 = fs * fac10[lay - 1]

        if specparm1 < 0.125:
            p = fs1 - 1.0
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac001 = fk0 * fac01[lay - 1]
            fac101 = fk1 * fac01[lay - 1]
            fac201 = fk2 * fac01[lay - 1]
            fac011 = fk0 * fac11[lay - 1]
            fac111 = fk1 * fac11[lay - 1]
            fac211 = fk2 * fac11[lay - 1]
        elif specparm1 > 0.875:
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac001 = fk0 * fac01[lay - 1]
            fac101 = fk1 * fac01[lay - 1]
            fac201 = fk2 * fac01[lay - 1]
            fac011 = fk0 * fac11[lay - 1]
            fac111 = fk1 * fac11[lay - 1]
            fac211 = fk2 * fac11[lay - 1]
        else:
            fac001 = (1.0 - fs1) * fac01[lay - 1]
            fac011 = (1.0 - fs1) * fac11[lay - 1]
            fac101 = fs1 * fac01[lay - 1]
            fac111 = fs1 * fac11[lay - 1]

        for ig in range(1, ng5 + 1):
            tauself = selffac[lay - 1] * (
                selfref[_idx2(inds, ig, 10)]
                + selffrac[lay - 1]
                * (selfref[_idx2(inds + 1, ig, 10)] - selfref[_idx2(inds, ig, 10)])
            )
            taufor = forfac[lay - 1] * (
                forref[_idx2(indf, ig, 4)]
                + forfrac[lay - 1]
                * (forref[_idx2(indf + 1, ig, 4)] - forref[_idx2(indf, ig, 4)])
            )
            o3m1 = ka_mo3[_idx3(jmo3, indm, ig, 9, 19)] + fmo3 * (
                ka_mo3[_idx3(jmo3 + 1, indm, ig, 9, 19)]
                - ka_mo3[_idx3(jmo3, indm, ig, 9, 19)]
            )
            o3m2 = ka_mo3[_idx3(jmo3, indm + 1, ig, 9, 19)] + fmo3 * (
                ka_mo3[_idx3(jmo3 + 1, indm + 1, ig, 9, 19)]
                - ka_mo3[_idx3(jmo3, indm + 1, ig, 9, 19)]
            )
            abso3 = o3m1 + minorfrac[lay - 1] * (o3m2 - o3m1)
            tau_major = _lw_tau_major_3pt(
                specparm, speccomb, ind0, ig, fac000, fac100, fac200,
                fac010, fac110, fac210, absa, 585
            )
            tau_major1 = _lw_tau_major_3pt(
                specparm1, speccomb1, ind1, ig, fac001, fac101, fac201,
                fac011, fac111, fac211, absa, 585
            )
            taug[_idx2(lay, ngs4 + ig, nlayers)] = (
                tau_major
                + tau_major1
                + tauself
                + taufor
                + abso3 * colo3[lay - 1]
                + wx[_idx2(1, lay, maxxsec)] * ccl4[ig - 1]
            )
            fracs[_idx2(lay, ngs4 + ig, nlayers)] = fracrefa[
                _idx2(ig, jpl, ng5)
            ] + fpl * (
                fracrefa[_idx2(ig, jpl + 1, ng5)] - fracrefa[_idx2(ig, jpl, ng5)]
            )

    for lay in range(laytrop + 1, nlayers + 1):
        speccomb = colo3[lay - 1] + rat_o3co2[lay - 1] * colco2[lay - 1]
        specparm = colo3[lay - 1] / speccomb
        if specparm >= oneminus:
            specparm = oneminus
        specmult = 4.0 * specparm
        js = 1 + int(specmult)
        fs = specmult - float(int(specmult))

        speccomb1 = colo3[lay - 1] + rat_o3co2_1[lay - 1] * colco2[lay - 1]
        specparm1 = colo3[lay - 1] / speccomb1
        if specparm1 >= oneminus:
            specparm1 = oneminus
        specmult1 = 4.0 * specparm1
        js1 = 1 + int(specmult1)
        fs1 = specmult1 - float(int(specmult1))

        fac000 = (1.0 - fs) * fac00[lay - 1]
        fac010 = (1.0 - fs) * fac10[lay - 1]
        fac100 = fs * fac00[lay - 1]
        fac110 = fs * fac10[lay - 1]
        fac001 = (1.0 - fs1) * fac01[lay - 1]
        fac011 = (1.0 - fs1) * fac11[lay - 1]
        fac101 = fs1 * fac01[lay - 1]
        fac111 = fs1 * fac11[lay - 1]

        speccomb_planck = colo3[lay - 1] + refrat_planck_b * colco2[lay - 1]
        specparm_planck = colo3[lay - 1] / speccomb_planck
        if specparm_planck >= oneminus:
            specparm_planck = oneminus
        specmult_planck = 4.0 * specparm_planck
        jpl = 1 + int(specmult_planck)
        fpl = specmult_planck - float(int(specmult_planck))

        ind0 = ((jp[lay - 1] - 13) * 5 + (jt[lay - 1] - 1)) * nspb5 + js
        ind1 = ((jp[lay - 1] - 12) * 5 + (jt1[lay - 1] - 1)) * nspb5 + js1

        for ig in range(1, ng5 + 1):
            taug[_idx2(lay, ngs4 + ig, nlayers)] = (
                speccomb
                * (
                    fac000 * absb[_idx2(ind0, ig, 1175)]
                    + fac100 * absb[_idx2(ind0 + 1, ig, 1175)]
                    + fac010 * absb[_idx2(ind0 + 5, ig, 1175)]
                    + fac110 * absb[_idx2(ind0 + 6, ig, 1175)]
                )
                + speccomb1
                * (
                    fac001 * absb[_idx2(ind1, ig, 1175)]
                    + fac101 * absb[_idx2(ind1 + 1, ig, 1175)]
                    + fac011 * absb[_idx2(ind1 + 5, ig, 1175)]
                    + fac111 * absb[_idx2(ind1 + 6, ig, 1175)]
                )
                + wx[_idx2(1, lay, maxxsec)] * ccl4[ig - 1]
            )
            fracs[_idx2(lay, ngs4 + ig, nlayers)] = fracrefb[
                _idx2(ig, jpl, ng5)
            ] + fpl * (
                fracrefb[_idx2(ig, jpl + 1, ng5)] - fracrefb[_idx2(ig, jpl, ng5)]
            )


@export
def rrtmg_lw_taugb3_codon(
    nlayers: int,
    laytrop: int,
    ng3: int,
    ngs2: int,
    nspa3: int,
    nspb3: int,
    oneminus: float,
    colh2o_p: cobj,
    colco2_p: cobj,
    coln2o_p: cobj,
    coldry_p: cobj,
    rat_h2oco2_p: cobj,
    rat_h2oco2_1_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    indself_p: cobj,
    indfor_p: cobj,
    indminor_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    selffac_p: cobj,
    selffrac_p: cobj,
    forfac_p: cobj,
    forfrac_p: cobj,
    minorfrac_p: cobj,
    chi_mls_p: cobj,
    fracrefa_p: cobj,
    fracrefb_p: cobj,
    absa_p: cobj,
    absb_p: cobj,
    ka_mn2o_p: cobj,
    kb_mn2o_p: cobj,
    selfref_p: cobj,
    forref_p: cobj,
    fracs_p: cobj,
    taug_p: cobj,
):
    colh2o = Ptr[float](colh2o_p)
    colco2 = Ptr[float](colco2_p)
    coln2o = Ptr[float](coln2o_p)
    coldry = Ptr[float](coldry_p)
    rat_h2oco2 = Ptr[float](rat_h2oco2_p)
    rat_h2oco2_1 = Ptr[float](rat_h2oco2_1_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    indself = Ptr[int](indself_p)
    indfor = Ptr[int](indfor_p)
    indminor = Ptr[int](indminor_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    selffac = Ptr[float](selffac_p)
    selffrac = Ptr[float](selffrac_p)
    forfac = Ptr[float](forfac_p)
    forfrac = Ptr[float](forfrac_p)
    minorfrac = Ptr[float](minorfrac_p)
    chi_mls = Ptr[float](chi_mls_p)
    fracrefa = Ptr[float](fracrefa_p)
    fracrefb = Ptr[float](fracrefb_p)
    absa = Ptr[float](absa_p)
    absb = Ptr[float](absb_p)
    ka_mn2o = Ptr[float](ka_mn2o_p)
    kb_mn2o = Ptr[float](kb_mn2o_p)
    selfref = Ptr[float](selfref_p)
    forref = Ptr[float](forref_p)
    fracs = Ptr[float](fracs_p)
    taug = Ptr[float](taug_p)

    refrat_planck_a = chi_mls[_idx2(1, 9, 7)] / chi_mls[_idx2(2, 9, 7)]
    refrat_planck_b = chi_mls[_idx2(1, 13, 7)] / chi_mls[_idx2(2, 13, 7)]
    refrat_m_a = chi_mls[_idx2(1, 3, 7)] / chi_mls[_idx2(2, 3, 7)]
    refrat_m_b = chi_mls[_idx2(1, 13, 7)] / chi_mls[_idx2(2, 13, 7)]

    for lay in range(1, laytrop + 1):
        speccomb = colh2o[lay - 1] + rat_h2oco2[lay - 1] * colco2[lay - 1]
        specparm = colh2o[lay - 1] / speccomb
        if specparm >= oneminus:
            specparm = oneminus
        specmult = 8.0 * specparm
        js = 1 + int(specmult)
        fs = specmult - float(int(specmult))

        speccomb1 = colh2o[lay - 1] + rat_h2oco2_1[lay - 1] * colco2[lay - 1]
        specparm1 = colh2o[lay - 1] / speccomb1
        if specparm1 >= oneminus:
            specparm1 = oneminus
        specmult1 = 8.0 * specparm1
        js1 = 1 + int(specmult1)
        fs1 = specmult1 - float(int(specmult1))

        speccomb_mn2o = colh2o[lay - 1] + refrat_m_a * colco2[lay - 1]
        specparm_mn2o = colh2o[lay - 1] / speccomb_mn2o
        if specparm_mn2o >= oneminus:
            specparm_mn2o = oneminus
        specmult_mn2o = 8.0 * specparm_mn2o
        jmn2o = 1 + int(specmult_mn2o)
        fmn2o = specmult_mn2o - float(int(specmult_mn2o))

        chi_n2o = coln2o[lay - 1] / coldry[lay - 1]
        ratn2o = 1.0e20 * chi_n2o / chi_mls[_idx2(4, jp[lay - 1] + 1, 7)]
        if ratn2o > 1.5:
            adjfac = 0.5 + (ratn2o - 0.5) ** 0.65
            adjcoln2o = (
                adjfac
                * chi_mls[_idx2(4, jp[lay - 1] + 1, 7)]
                * coldry[lay - 1]
                * 1.0e-20
            )
        else:
            adjcoln2o = coln2o[lay - 1]

        speccomb_planck = colh2o[lay - 1] + refrat_planck_a * colco2[lay - 1]
        specparm_planck = colh2o[lay - 1] / speccomb_planck
        if specparm_planck >= oneminus:
            specparm_planck = oneminus
        specmult_planck = 8.0 * specparm_planck
        jpl = 1 + int(specmult_planck)
        fpl = specmult_planck - float(int(specmult_planck))

        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa3 + js
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa3 + js1
        inds = indself[lay - 1]
        indf = indfor[lay - 1]
        indm = indminor[lay - 1]
        fac200 = 0.0
        fac210 = 0.0
        fac201 = 0.0
        fac211 = 0.0

        if specparm < 0.125:
            p = fs - 1.0
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac000 = fk0 * fac00[lay - 1]
            fac100 = fk1 * fac00[lay - 1]
            fac200 = fk2 * fac00[lay - 1]
            fac010 = fk0 * fac10[lay - 1]
            fac110 = fk1 * fac10[lay - 1]
            fac210 = fk2 * fac10[lay - 1]
        elif specparm > 0.875:
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac000 = fk0 * fac00[lay - 1]
            fac100 = fk1 * fac00[lay - 1]
            fac200 = fk2 * fac00[lay - 1]
            fac010 = fk0 * fac10[lay - 1]
            fac110 = fk1 * fac10[lay - 1]
            fac210 = fk2 * fac10[lay - 1]
        else:
            fac000 = (1.0 - fs) * fac00[lay - 1]
            fac010 = (1.0 - fs) * fac10[lay - 1]
            fac100 = fs * fac00[lay - 1]
            fac110 = fs * fac10[lay - 1]

        if specparm1 < 0.125:
            p = fs1 - 1.0
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac001 = fk0 * fac01[lay - 1]
            fac101 = fk1 * fac01[lay - 1]
            fac201 = fk2 * fac01[lay - 1]
            fac011 = fk0 * fac11[lay - 1]
            fac111 = fk1 * fac11[lay - 1]
            fac211 = fk2 * fac11[lay - 1]
        elif specparm1 > 0.875:
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac001 = fk0 * fac01[lay - 1]
            fac101 = fk1 * fac01[lay - 1]
            fac201 = fk2 * fac01[lay - 1]
            fac011 = fk0 * fac11[lay - 1]
            fac111 = fk1 * fac11[lay - 1]
            fac211 = fk2 * fac11[lay - 1]
        else:
            fac001 = (1.0 - fs1) * fac01[lay - 1]
            fac011 = (1.0 - fs1) * fac11[lay - 1]
            fac101 = fs1 * fac01[lay - 1]
            fac111 = fs1 * fac11[lay - 1]

        for ig in range(1, ng3 + 1):
            tauself = selffac[lay - 1] * (
                selfref[_idx2(inds, ig, 10)]
                + selffrac[lay - 1]
                * (selfref[_idx2(inds + 1, ig, 10)] - selfref[_idx2(inds, ig, 10)])
            )
            taufor = forfac[lay - 1] * (
                forref[_idx2(indf, ig, 4)]
                + forfrac[lay - 1]
                * (forref[_idx2(indf + 1, ig, 4)] - forref[_idx2(indf, ig, 4)])
            )
            n2om1 = ka_mn2o[_idx3(jmn2o, indm, ig, 9, 19)] + fmn2o * (
                ka_mn2o[_idx3(jmn2o + 1, indm, ig, 9, 19)]
                - ka_mn2o[_idx3(jmn2o, indm, ig, 9, 19)]
            )
            n2om2 = ka_mn2o[_idx3(jmn2o, indm + 1, ig, 9, 19)] + fmn2o * (
                ka_mn2o[_idx3(jmn2o + 1, indm + 1, ig, 9, 19)]
                - ka_mn2o[_idx3(jmn2o, indm + 1, ig, 9, 19)]
            )
            absn2o = n2om1 + minorfrac[lay - 1] * (n2om2 - n2om1)
            tau_major = _lw_tau_major_3pt(
                specparm, speccomb, ind0, ig, fac000, fac100, fac200,
                fac010, fac110, fac210, absa, 585
            )
            tau_major1 = _lw_tau_major_3pt(
                specparm1, speccomb1, ind1, ig, fac001, fac101, fac201,
                fac011, fac111, fac211, absa, 585
            )
            taug[_idx2(lay, ngs2 + ig, nlayers)] = (
                tau_major + tau_major1 + tauself + taufor + adjcoln2o * absn2o
            )
            fracs[_idx2(lay, ngs2 + ig, nlayers)] = fracrefa[
                _idx2(ig, jpl, ng3)
            ] + fpl * (
                fracrefa[_idx2(ig, jpl + 1, ng3)] - fracrefa[_idx2(ig, jpl, ng3)]
            )

    for lay in range(laytrop + 1, nlayers + 1):
        speccomb = colh2o[lay - 1] + rat_h2oco2[lay - 1] * colco2[lay - 1]
        specparm = colh2o[lay - 1] / speccomb
        if specparm >= oneminus:
            specparm = oneminus
        specmult = 4.0 * specparm
        js = 1 + int(specmult)
        fs = specmult - float(int(specmult))

        speccomb1 = colh2o[lay - 1] + rat_h2oco2_1[lay - 1] * colco2[lay - 1]
        specparm1 = colh2o[lay - 1] / speccomb1
        if specparm1 >= oneminus:
            specparm1 = oneminus
        specmult1 = 4.0 * specparm1
        js1 = 1 + int(specmult1)
        fs1 = specmult1 - float(int(specmult1))

        fac000 = (1.0 - fs) * fac00[lay - 1]
        fac010 = (1.0 - fs) * fac10[lay - 1]
        fac100 = fs * fac00[lay - 1]
        fac110 = fs * fac10[lay - 1]
        fac001 = (1.0 - fs1) * fac01[lay - 1]
        fac011 = (1.0 - fs1) * fac11[lay - 1]
        fac101 = fs1 * fac01[lay - 1]
        fac111 = fs1 * fac11[lay - 1]

        speccomb_mn2o = colh2o[lay - 1] + refrat_m_b * colco2[lay - 1]
        specparm_mn2o = colh2o[lay - 1] / speccomb_mn2o
        if specparm_mn2o >= oneminus:
            specparm_mn2o = oneminus
        specmult_mn2o = 4.0 * specparm_mn2o
        jmn2o = 1 + int(specmult_mn2o)
        fmn2o = specmult_mn2o - float(int(specmult_mn2o))

        chi_n2o = coln2o[lay - 1] / coldry[lay - 1]
        ratn2o = 1.0000000200408773e20 * chi_n2o / chi_mls[_idx2(4, jp[lay - 1] + 1, 7)]
        if ratn2o > 1.5:
            adjfac = 0.5 + (ratn2o - 0.5) ** 0.65
            adjcoln2o = (
                adjfac
                * chi_mls[_idx2(4, jp[lay - 1] + 1, 7)]
                * coldry[lay - 1]
                * 1.0e-20
            )
        else:
            adjcoln2o = coln2o[lay - 1]

        speccomb_planck = colh2o[lay - 1] + refrat_planck_b * colco2[lay - 1]
        specparm_planck = colh2o[lay - 1] / speccomb_planck
        if specparm_planck >= oneminus:
            specparm_planck = oneminus
        specmult_planck = 4.0 * specparm_planck
        jpl = 1 + int(specmult_planck)
        fpl = specmult_planck - float(int(specmult_planck))

        ind0 = ((jp[lay - 1] - 13) * 5 + (jt[lay - 1] - 1)) * nspb3 + js
        ind1 = ((jp[lay - 1] - 12) * 5 + (jt1[lay - 1] - 1)) * nspb3 + js1
        indf = indfor[lay - 1]
        indm = indminor[lay - 1]

        for ig in range(1, ng3 + 1):
            taufor = forfac[lay - 1] * (
                forref[_idx2(indf, ig, 4)]
                + forfrac[lay - 1]
                * (forref[_idx2(indf + 1, ig, 4)] - forref[_idx2(indf, ig, 4)])
            )
            n2om1 = kb_mn2o[_idx3(jmn2o, indm, ig, 5, 19)] + fmn2o * (
                kb_mn2o[_idx3(jmn2o + 1, indm, ig, 5, 19)]
                - kb_mn2o[_idx3(jmn2o, indm, ig, 5, 19)]
            )
            n2om2 = kb_mn2o[_idx3(jmn2o, indm + 1, ig, 5, 19)] + fmn2o * (
                kb_mn2o[_idx3(jmn2o + 1, indm + 1, ig, 5, 19)]
                - kb_mn2o[_idx3(jmn2o, indm + 1, ig, 5, 19)]
            )
            absn2o = n2om1 + minorfrac[lay - 1] * (n2om2 - n2om1)
            taug[_idx2(lay, ngs2 + ig, nlayers)] = (
                speccomb
                * (
                    fac000 * absb[_idx2(ind0, ig, 1175)]
                    + fac100 * absb[_idx2(ind0 + 1, ig, 1175)]
                    + fac010 * absb[_idx2(ind0 + 5, ig, 1175)]
                    + fac110 * absb[_idx2(ind0 + 6, ig, 1175)]
                )
                + speccomb1
                * (
                    fac001 * absb[_idx2(ind1, ig, 1175)]
                    + fac101 * absb[_idx2(ind1 + 1, ig, 1175)]
                    + fac011 * absb[_idx2(ind1 + 5, ig, 1175)]
                    + fac111 * absb[_idx2(ind1 + 6, ig, 1175)]
                )
                + taufor
                + adjcoln2o * absn2o
            )
            fracs[_idx2(lay, ngs2 + ig, nlayers)] = fracrefb[
                _idx2(ig, jpl, ng3)
            ] + fpl * (
                fracrefb[_idx2(ig, jpl + 1, ng3)] - fracrefb[_idx2(ig, jpl, ng3)]
            )


@export
def rrtmg_lw_taugb7_codon(
    nlayers: int,
    laytrop: int,
    ng7: int,
    ngs6: int,
    nspa7: int,
    nspb7: int,
    oneminus: float,
    colh2o_p: cobj,
    colco2_p: cobj,
    colo3_p: cobj,
    coldry_p: cobj,
    rat_h2oo3_p: cobj,
    rat_h2oo3_1_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    indself_p: cobj,
    indfor_p: cobj,
    indminor_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    selffac_p: cobj,
    selffrac_p: cobj,
    forfac_p: cobj,
    forfrac_p: cobj,
    minorfrac_p: cobj,
    chi_mls_p: cobj,
    fracrefa_p: cobj,
    fracrefb_p: cobj,
    absa_p: cobj,
    absb_p: cobj,
    ka_mco2_p: cobj,
    kb_mco2_p: cobj,
    selfref_p: cobj,
    forref_p: cobj,
    fracs_p: cobj,
    taug_p: cobj,
):
    colh2o = Ptr[float](colh2o_p)
    colco2 = Ptr[float](colco2_p)
    colo3 = Ptr[float](colo3_p)
    coldry = Ptr[float](coldry_p)
    rat_h2oo3 = Ptr[float](rat_h2oo3_p)
    rat_h2oo3_1 = Ptr[float](rat_h2oo3_1_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    indself = Ptr[int](indself_p)
    indfor = Ptr[int](indfor_p)
    indminor = Ptr[int](indminor_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    selffac = Ptr[float](selffac_p)
    selffrac = Ptr[float](selffrac_p)
    forfac = Ptr[float](forfac_p)
    forfrac = Ptr[float](forfrac_p)
    minorfrac = Ptr[float](minorfrac_p)
    chi_mls = Ptr[float](chi_mls_p)
    fracrefa = Ptr[float](fracrefa_p)
    fracrefb = Ptr[float](fracrefb_p)
    absa = Ptr[float](absa_p)
    absb = Ptr[float](absb_p)
    ka_mco2 = Ptr[float](ka_mco2_p)
    kb_mco2 = Ptr[float](kb_mco2_p)
    selfref = Ptr[float](selfref_p)
    forref = Ptr[float](forref_p)
    fracs = Ptr[float](fracs_p)
    taug = Ptr[float](taug_p)

    refrat_planck_a = chi_mls[_idx2(1, 3, 7)] / chi_mls[_idx2(3, 3, 7)]
    refrat_m_a = chi_mls[_idx2(1, 3, 7)] / chi_mls[_idx2(3, 3, 7)]

    for lay in range(1, laytrop + 1):
        speccomb = colh2o[lay - 1] + rat_h2oo3[lay - 1] * colo3[lay - 1]
        specparm = colh2o[lay - 1] / speccomb
        if specparm >= oneminus:
            specparm = oneminus
        specmult = 8.0 * specparm
        js = 1 + int(specmult)
        fs = specmult - float(int(specmult))

        speccomb1 = colh2o[lay - 1] + rat_h2oo3_1[lay - 1] * colo3[lay - 1]
        specparm1 = colh2o[lay - 1] / speccomb1
        if specparm1 >= oneminus:
            specparm1 = oneminus
        specmult1 = 8.0 * specparm1
        js1 = 1 + int(specmult1)
        fs1 = specmult1 - float(int(specmult1))

        speccomb_mco2 = colh2o[lay - 1] + refrat_m_a * colo3[lay - 1]
        specparm_mco2 = colh2o[lay - 1] / speccomb_mco2
        if specparm_mco2 >= oneminus:
            specparm_mco2 = oneminus
        specmult_mco2 = 8.0 * specparm_mco2
        jmco2 = 1 + int(specmult_mco2)
        fmco2 = specmult_mco2 - float(int(specmult_mco2))

        chi_co2 = colco2[lay - 1] / coldry[lay - 1]
        ratco2 = 1.0000000200408773e20 * chi_co2 / chi_mls[_idx2(2, jp[lay - 1] + 1, 7)]
        if ratco2 > 3.0:
            adjfac = 3.0 + (ratco2 - 3.0) ** 0.79
            adjcolco2 = (
                adjfac
                * chi_mls[_idx2(2, jp[lay - 1] + 1, 7)]
                * coldry[lay - 1]
                * 1.0e-20
            )
        else:
            adjcolco2 = colco2[lay - 1]

        speccomb_planck = colh2o[lay - 1] + refrat_planck_a * colo3[lay - 1]
        specparm_planck = colh2o[lay - 1] / speccomb_planck
        if specparm_planck >= oneminus:
            specparm_planck = oneminus
        specmult_planck = 8.0 * specparm_planck
        jpl = 1 + int(specmult_planck)
        fpl = specmult_planck - float(int(specmult_planck))

        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa7 + js
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa7 + js1
        inds = indself[lay - 1]
        indf = indfor[lay - 1]
        indm = indminor[lay - 1]
        fac200 = 0.0
        fac210 = 0.0
        fac201 = 0.0
        fac211 = 0.0

        if specparm < 0.125:
            p = fs - 1.0
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac000 = fk0 * fac00[lay - 1]
            fac100 = fk1 * fac00[lay - 1]
            fac200 = fk2 * fac00[lay - 1]
            fac010 = fk0 * fac10[lay - 1]
            fac110 = fk1 * fac10[lay - 1]
            fac210 = fk2 * fac10[lay - 1]
        elif specparm > 0.875:
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac000 = fk0 * fac00[lay - 1]
            fac100 = fk1 * fac00[lay - 1]
            fac200 = fk2 * fac00[lay - 1]
            fac010 = fk0 * fac10[lay - 1]
            fac110 = fk1 * fac10[lay - 1]
            fac210 = fk2 * fac10[lay - 1]
        else:
            fac000 = (1.0 - fs) * fac00[lay - 1]
            fac010 = (1.0 - fs) * fac10[lay - 1]
            fac100 = fs * fac00[lay - 1]
            fac110 = fs * fac10[lay - 1]

        if specparm < 0.125:
            p = fs1 - 1.0
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac001 = fk0 * fac01[lay - 1]
            fac101 = fk1 * fac01[lay - 1]
            fac201 = fk2 * fac01[lay - 1]
            fac011 = fk0 * fac11[lay - 1]
            fac111 = fk1 * fac11[lay - 1]
            fac211 = fk2 * fac11[lay - 1]
        elif specparm1 > 0.875:
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac001 = fk0 * fac01[lay - 1]
            fac101 = fk1 * fac01[lay - 1]
            fac201 = fk2 * fac01[lay - 1]
            fac011 = fk0 * fac11[lay - 1]
            fac111 = fk1 * fac11[lay - 1]
            fac211 = fk2 * fac11[lay - 1]
        else:
            fac001 = (1.0 - fs1) * fac01[lay - 1]
            fac011 = (1.0 - fs1) * fac11[lay - 1]
            fac101 = fs1 * fac01[lay - 1]
            fac111 = fs1 * fac11[lay - 1]

        for ig in range(1, ng7 + 1):
            tauself = selffac[lay - 1] * (
                selfref[_idx2(inds, ig, 10)]
                + selffrac[lay - 1]
                * (selfref[_idx2(inds + 1, ig, 10)] - selfref[_idx2(inds, ig, 10)])
            )
            taufor = forfac[lay - 1] * (
                forref[_idx2(indf, ig, 4)]
                + forfrac[lay - 1]
                * (forref[_idx2(indf + 1, ig, 4)] - forref[_idx2(indf, ig, 4)])
            )
            co2m1 = ka_mco2[_idx3(jmco2, indm, ig, 9, 19)] + fmco2 * (
                ka_mco2[_idx3(jmco2 + 1, indm, ig, 9, 19)]
                - ka_mco2[_idx3(jmco2, indm, ig, 9, 19)]
            )
            co2m2 = ka_mco2[_idx3(jmco2, indm + 1, ig, 9, 19)] + fmco2 * (
                ka_mco2[_idx3(jmco2 + 1, indm + 1, ig, 9, 19)]
                - ka_mco2[_idx3(jmco2, indm + 1, ig, 9, 19)]
            )
            absco2 = co2m1 + minorfrac[lay - 1] * (co2m2 - co2m1)
            tau_major = _lw_tau_major_3pt(
                specparm, speccomb, ind0, ig, fac000, fac100, fac200,
                fac010, fac110, fac210, absa, 585
            )
            tau_major1 = _lw_tau_major_3pt(
                specparm1, speccomb1, ind1, ig, fac001, fac101, fac201,
                fac011, fac111, fac211, absa, 585
            )
            taug[_idx2(lay, ngs6 + ig, nlayers)] = (
                tau_major + tau_major1 + tauself + taufor + adjcolco2 * absco2
            )
            fracs[_idx2(lay, ngs6 + ig, nlayers)] = fracrefa[
                _idx2(ig, jpl, ng7)
            ] + fpl * (
                fracrefa[_idx2(ig, jpl + 1, ng7)] - fracrefa[_idx2(ig, jpl, ng7)]
            )

    for lay in range(laytrop + 1, nlayers + 1):
        chi_co2 = colco2[lay - 1] / coldry[lay - 1]
        ratco2 = 1.0000000200408773e20 * chi_co2 / chi_mls[_idx2(2, jp[lay - 1] + 1, 7)]
        if ratco2 > 3.0:
            adjfac = 2.0 + (ratco2 - 2.0) ** 0.79
            adjcolco2 = (
                adjfac
                * chi_mls[_idx2(2, jp[lay - 1] + 1, 7)]
                * coldry[lay - 1]
                * 1.0e-20
            )
        else:
            adjcolco2 = colco2[lay - 1]

        ind0 = ((jp[lay - 1] - 13) * 5 + (jt[lay - 1] - 1)) * nspb7 + 1
        ind1 = ((jp[lay - 1] - 12) * 5 + (jt1[lay - 1] - 1)) * nspb7 + 1
        indm = indminor[lay - 1]

        for ig in range(1, ng7 + 1):
            absco2 = kb_mco2[_idx2(indm, ig, 19)] + minorfrac[lay - 1] * (
                kb_mco2[_idx2(indm + 1, ig, 19)] - kb_mco2[_idx2(indm, ig, 19)]
            )
            taug[_idx2(lay, ngs6 + ig, nlayers)] = (
                colo3[lay - 1]
                * (
                    fac00[lay - 1] * absb[_idx2(ind0, ig, 235)]
                    + fac10[lay - 1] * absb[_idx2(ind0 + 1, ig, 235)]
                    + fac01[lay - 1] * absb[_idx2(ind1, ig, 235)]
                    + fac11[lay - 1] * absb[_idx2(ind1 + 1, ig, 235)]
                )
                + adjcolco2 * absco2
            )
            fracs[_idx2(lay, ngs6 + ig, nlayers)] = fracrefb[ig - 1]

        taug[_idx2(lay, ngs6 + 6, nlayers)] = taug[_idx2(lay, ngs6 + 6, nlayers)] * 0.92
        taug[_idx2(lay, ngs6 + 7, nlayers)] = taug[_idx2(lay, ngs6 + 7, nlayers)] * 0.88
        taug[_idx2(lay, ngs6 + 8, nlayers)] = taug[_idx2(lay, ngs6 + 8, nlayers)] * 1.07
        taug[_idx2(lay, ngs6 + 9, nlayers)] = taug[_idx2(lay, ngs6 + 9, nlayers)] * 1.1
        taug[_idx2(lay, ngs6 + 10, nlayers)] = taug[_idx2(lay, ngs6 + 10, nlayers)] * 0.99
        taug[_idx2(lay, ngs6 + 11, nlayers)] = taug[_idx2(lay, ngs6 + 11, nlayers)] * 0.855


@export
def rrtmg_lw_taugb9_codon(
    nlayers: int,
    laytrop: int,
    ng9: int,
    ngs8: int,
    nspa9: int,
    nspb9: int,
    oneminus: float,
    colh2o_p: cobj,
    colch4_p: cobj,
    coln2o_p: cobj,
    coldry_p: cobj,
    rat_h2och4_p: cobj,
    rat_h2och4_1_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    indself_p: cobj,
    indfor_p: cobj,
    indminor_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    selffac_p: cobj,
    selffrac_p: cobj,
    forfac_p: cobj,
    forfrac_p: cobj,
    minorfrac_p: cobj,
    chi_mls_p: cobj,
    fracrefa_p: cobj,
    fracrefb_p: cobj,
    absa_p: cobj,
    absb_p: cobj,
    ka_mn2o_p: cobj,
    kb_mn2o_p: cobj,
    selfref_p: cobj,
    forref_p: cobj,
    fracs_p: cobj,
    taug_p: cobj,
):
    colh2o = Ptr[float](colh2o_p)
    colch4 = Ptr[float](colch4_p)
    coln2o = Ptr[float](coln2o_p)
    coldry = Ptr[float](coldry_p)
    rat_h2och4 = Ptr[float](rat_h2och4_p)
    rat_h2och4_1 = Ptr[float](rat_h2och4_1_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    indself = Ptr[int](indself_p)
    indfor = Ptr[int](indfor_p)
    indminor = Ptr[int](indminor_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    selffac = Ptr[float](selffac_p)
    selffrac = Ptr[float](selffrac_p)
    forfac = Ptr[float](forfac_p)
    forfrac = Ptr[float](forfrac_p)
    minorfrac = Ptr[float](minorfrac_p)
    chi_mls = Ptr[float](chi_mls_p)
    fracrefa = Ptr[float](fracrefa_p)
    fracrefb = Ptr[float](fracrefb_p)
    absa = Ptr[float](absa_p)
    absb = Ptr[float](absb_p)
    ka_mn2o = Ptr[float](ka_mn2o_p)
    kb_mn2o = Ptr[float](kb_mn2o_p)
    selfref = Ptr[float](selfref_p)
    forref = Ptr[float](forref_p)
    fracs = Ptr[float](fracs_p)
    taug = Ptr[float](taug_p)

    refrat_planck_a = chi_mls[_idx2(1, 9, 7)] / chi_mls[_idx2(6, 9, 7)]
    refrat_m_a = chi_mls[_idx2(1, 3, 7)] / chi_mls[_idx2(6, 3, 7)]

    for lay in range(1, laytrop + 1):
        speccomb = colh2o[lay - 1] + rat_h2och4[lay - 1] * colch4[lay - 1]
        specparm = colh2o[lay - 1] / speccomb
        if specparm >= oneminus:
            specparm = oneminus
        specmult = 8.0 * specparm
        js = 1 + int(specmult)
        fs = specmult - float(int(specmult))

        speccomb1 = colh2o[lay - 1] + rat_h2och4_1[lay - 1] * colch4[lay - 1]
        specparm1 = colh2o[lay - 1] / speccomb1
        if specparm1 >= oneminus:
            specparm1 = oneminus
        specmult1 = 8.0 * specparm1
        js1 = 1 + int(specmult1)
        fs1 = specmult1 - float(int(specmult1))

        speccomb_mn2o = colh2o[lay - 1] + refrat_m_a * colch4[lay - 1]
        specparm_mn2o = colh2o[lay - 1] / speccomb_mn2o
        if specparm_mn2o >= oneminus:
            specparm_mn2o = oneminus
        specmult_mn2o = 8.0 * specparm_mn2o
        jmn2o = 1 + int(specmult_mn2o)
        fmn2o = specmult_mn2o - float(int(specmult_mn2o))

        chi_n2o = coln2o[lay - 1] / coldry[lay - 1]
        ratn2o = 1.0e20 * chi_n2o / chi_mls[_idx2(4, jp[lay - 1] + 1, 7)]
        if ratn2o > 1.5:
            adjfac = 0.5 + (ratn2o - 0.5) ** 0.65
            adjcoln2o = (
                adjfac
                * chi_mls[_idx2(4, jp[lay - 1] + 1, 7)]
                * coldry[lay - 1]
                * 1.0e-20
            )
        else:
            adjcoln2o = coln2o[lay - 1]

        speccomb_planck = colh2o[lay - 1] + refrat_planck_a * colch4[lay - 1]
        specparm_planck = colh2o[lay - 1] / speccomb_planck
        if specparm_planck >= oneminus:
            specparm_planck = oneminus
        specmult_planck = 8.0 * specparm_planck
        jpl = 1 + int(specmult_planck)
        fpl = specmult_planck - float(int(specmult_planck))

        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa9 + js
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa9 + js1
        inds = indself[lay - 1]
        indf = indfor[lay - 1]
        indm = indminor[lay - 1]
        fac200 = 0.0
        fac210 = 0.0
        fac201 = 0.0
        fac211 = 0.0

        if specparm < 0.125:
            p = fs - 1.0
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac000 = fk0 * fac00[lay - 1]
            fac100 = fk1 * fac00[lay - 1]
            fac200 = fk2 * fac00[lay - 1]
            fac010 = fk0 * fac10[lay - 1]
            fac110 = fk1 * fac10[lay - 1]
            fac210 = fk2 * fac10[lay - 1]
        elif specparm > 0.875:
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac000 = fk0 * fac00[lay - 1]
            fac100 = fk1 * fac00[lay - 1]
            fac200 = fk2 * fac00[lay - 1]
            fac010 = fk0 * fac10[lay - 1]
            fac110 = fk1 * fac10[lay - 1]
            fac210 = fk2 * fac10[lay - 1]
        else:
            fac000 = (1.0 - fs) * fac00[lay - 1]
            fac010 = (1.0 - fs) * fac10[lay - 1]
            fac100 = fs * fac00[lay - 1]
            fac110 = fs * fac10[lay - 1]

        if specparm1 < 0.125:
            p = fs1 - 1.0
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac001 = fk0 * fac01[lay - 1]
            fac101 = fk1 * fac01[lay - 1]
            fac201 = fk2 * fac01[lay - 1]
            fac011 = fk0 * fac11[lay - 1]
            fac111 = fk1 * fac11[lay - 1]
            fac211 = fk2 * fac11[lay - 1]
        elif specparm1 > 0.875:
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac001 = fk0 * fac01[lay - 1]
            fac101 = fk1 * fac01[lay - 1]
            fac201 = fk2 * fac01[lay - 1]
            fac011 = fk0 * fac11[lay - 1]
            fac111 = fk1 * fac11[lay - 1]
            fac211 = fk2 * fac11[lay - 1]
        else:
            fac001 = (1.0 - fs1) * fac01[lay - 1]
            fac011 = (1.0 - fs1) * fac11[lay - 1]
            fac101 = fs1 * fac01[lay - 1]
            fac111 = fs1 * fac11[lay - 1]

        for ig in range(1, ng9 + 1):
            tauself = selffac[lay - 1] * (
                selfref[_idx2(inds, ig, 10)]
                + selffrac[lay - 1]
                * (selfref[_idx2(inds + 1, ig, 10)] - selfref[_idx2(inds, ig, 10)])
            )
            taufor = forfac[lay - 1] * (
                forref[_idx2(indf, ig, 4)]
                + forfrac[lay - 1]
                * (forref[_idx2(indf + 1, ig, 4)] - forref[_idx2(indf, ig, 4)])
            )
            n2om1 = ka_mn2o[_idx3(jmn2o, indm, ig, 9, 19)] + fmn2o * (
                ka_mn2o[_idx3(jmn2o + 1, indm, ig, 9, 19)]
                - ka_mn2o[_idx3(jmn2o, indm, ig, 9, 19)]
            )
            n2om2 = ka_mn2o[_idx3(jmn2o, indm + 1, ig, 9, 19)] + fmn2o * (
                ka_mn2o[_idx3(jmn2o + 1, indm + 1, ig, 9, 19)]
                - ka_mn2o[_idx3(jmn2o, indm + 1, ig, 9, 19)]
            )
            absn2o = n2om1 + minorfrac[lay - 1] * (n2om2 - n2om1)
            tau_major = _lw_tau_major_3pt(
                specparm, speccomb, ind0, ig, fac000, fac100, fac200,
                fac010, fac110, fac210, absa, 585
            )
            tau_major1 = _lw_tau_major_3pt(
                specparm1, speccomb1, ind1, ig, fac001, fac101, fac201,
                fac011, fac111, fac211, absa, 585
            )
            taug[_idx2(lay, ngs8 + ig, nlayers)] = (
                tau_major + tau_major1 + tauself + taufor + adjcoln2o * absn2o
            )
            fracs[_idx2(lay, ngs8 + ig, nlayers)] = fracrefa[
                _idx2(ig, jpl, ng9)
            ] + fpl * (
                fracrefa[_idx2(ig, jpl + 1, ng9)] - fracrefa[_idx2(ig, jpl, ng9)]
            )

    for lay in range(laytrop + 1, nlayers + 1):
        chi_n2o = coln2o[lay - 1] / coldry[lay - 1]
        ratn2o = 1.0e20 * chi_n2o / chi_mls[_idx2(4, jp[lay - 1] + 1, 7)]
        if ratn2o > 1.5:
            adjfac = 0.5 + (ratn2o - 0.5) ** 0.65
            adjcoln2o = (
                adjfac
                * chi_mls[_idx2(4, jp[lay - 1] + 1, 7)]
                * coldry[lay - 1]
                * 1.0e-20
            )
        else:
            adjcoln2o = coln2o[lay - 1]

        ind0 = ((jp[lay - 1] - 13) * 5 + (jt[lay - 1] - 1)) * nspb9 + 1
        ind1 = ((jp[lay - 1] - 12) * 5 + (jt1[lay - 1] - 1)) * nspb9 + 1
        indm = indminor[lay - 1]

        for ig in range(1, ng9 + 1):
            absn2o = kb_mn2o[_idx2(indm, ig, 19)] + minorfrac[lay - 1] * (
                kb_mn2o[_idx2(indm + 1, ig, 19)] - kb_mn2o[_idx2(indm, ig, 19)]
            )
            taug[_idx2(lay, ngs8 + ig, nlayers)] = (
                colch4[lay - 1]
                * (
                    fac00[lay - 1] * absb[_idx2(ind0, ig, 235)]
                    + fac10[lay - 1] * absb[_idx2(ind0 + 1, ig, 235)]
                    + fac01[lay - 1] * absb[_idx2(ind1, ig, 235)]
                    + fac11[lay - 1] * absb[_idx2(ind1 + 1, ig, 235)]
                )
                + adjcoln2o * absn2o
            )
            fracs[_idx2(lay, ngs8 + ig, nlayers)] = fracrefb[ig - 1]


@export
def rrtmg_lw_taugb13_codon(
    nlayers: int,
    laytrop: int,
    ng13: int,
    ngs12: int,
    nspa13: int,
    oneminus: float,
    colh2o_p: cobj,
    coln2o_p: cobj,
    colco2_p: cobj,
    colco_p: cobj,
    colo3_p: cobj,
    coldry_p: cobj,
    rat_h2on2o_p: cobj,
    rat_h2on2o_1_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    indself_p: cobj,
    indfor_p: cobj,
    indminor_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    selffac_p: cobj,
    selffrac_p: cobj,
    forfac_p: cobj,
    forfrac_p: cobj,
    minorfrac_p: cobj,
    chi_mls_p: cobj,
    fracrefa_p: cobj,
    fracrefb_p: cobj,
    absa_p: cobj,
    ka_mco2_p: cobj,
    ka_mco_p: cobj,
    kb_mo3_p: cobj,
    selfref_p: cobj,
    forref_p: cobj,
    fracs_p: cobj,
    taug_p: cobj,
):
    colh2o = Ptr[float](colh2o_p)
    coln2o = Ptr[float](coln2o_p)
    colco2 = Ptr[float](colco2_p)
    colco = Ptr[float](colco_p)
    colo3 = Ptr[float](colo3_p)
    coldry = Ptr[float](coldry_p)
    rat_h2on2o = Ptr[float](rat_h2on2o_p)
    rat_h2on2o_1 = Ptr[float](rat_h2on2o_1_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    indself = Ptr[int](indself_p)
    indfor = Ptr[int](indfor_p)
    indminor = Ptr[int](indminor_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    selffac = Ptr[float](selffac_p)
    selffrac = Ptr[float](selffrac_p)
    forfac = Ptr[float](forfac_p)
    forfrac = Ptr[float](forfrac_p)
    minorfrac = Ptr[float](minorfrac_p)
    chi_mls = Ptr[float](chi_mls_p)
    fracrefa = Ptr[float](fracrefa_p)
    fracrefb = Ptr[float](fracrefb_p)
    absa = Ptr[float](absa_p)
    ka_mco2 = Ptr[float](ka_mco2_p)
    ka_mco = Ptr[float](ka_mco_p)
    kb_mo3 = Ptr[float](kb_mo3_p)
    selfref = Ptr[float](selfref_p)
    forref = Ptr[float](forref_p)
    fracs = Ptr[float](fracs_p)
    taug = Ptr[float](taug_p)

    refrat_planck_a = chi_mls[_idx2(1, 5, 7)] / chi_mls[_idx2(4, 5, 7)]
    refrat_m_a = chi_mls[_idx2(1, 1, 7)] / chi_mls[_idx2(4, 1, 7)]
    refrat_m_a3 = chi_mls[_idx2(1, 3, 7)] / chi_mls[_idx2(4, 3, 7)]

    for lay in range(1, laytrop + 1):
        speccomb = colh2o[lay - 1] + rat_h2on2o[lay - 1] * coln2o[lay - 1]
        specparm = colh2o[lay - 1] / speccomb
        if specparm >= oneminus:
            specparm = oneminus
        specmult = 8.0 * specparm
        js = 1 + int(specmult)
        fs = specmult - float(int(specmult))

        speccomb1 = colh2o[lay - 1] + rat_h2on2o_1[lay - 1] * coln2o[lay - 1]
        specparm1 = colh2o[lay - 1] / speccomb1
        if specparm1 >= oneminus:
            specparm1 = oneminus
        specmult1 = 8.0 * specparm1
        js1 = 1 + int(specmult1)
        fs1 = specmult1 - float(int(specmult1))

        speccomb_mco2 = colh2o[lay - 1] + refrat_m_a * coln2o[lay - 1]
        specparm_mco2 = colh2o[lay - 1] / speccomb_mco2
        if specparm_mco2 >= oneminus:
            specparm_mco2 = oneminus
        specmult_mco2 = 8.0 * specparm_mco2
        jmco2 = 1 + int(specmult_mco2)
        fmco2 = specmult_mco2 - float(int(specmult_mco2))

        chi_co2 = colco2[lay - 1] / coldry[lay - 1]
        ratco2 = 1.0e20 * chi_co2 / 3.55e-4
        if ratco2 > 3.0:
            adjfac = 2.0 + (ratco2 - 2.0) ** 0.68
            adjcolco2 = adjfac * 0.0003549999964889139 * coldry[lay - 1] * 1.0e-20
        else:
            adjcolco2 = colco2[lay - 1]

        speccomb_mco = colh2o[lay - 1] + refrat_m_a3 * coln2o[lay - 1]
        specparm_mco = colh2o[lay - 1] / speccomb_mco
        if specparm_mco >= oneminus:
            specparm_mco = oneminus
        specmult_mco = 8.0 * specparm_mco
        jmco = 1 + int(specmult_mco)
        fmco = specmult_mco - float(int(specmult_mco))

        speccomb_planck = colh2o[lay - 1] + refrat_planck_a * coln2o[lay - 1]
        specparm_planck = colh2o[lay - 1] / speccomb_planck
        if specparm_planck >= oneminus:
            specparm_planck = oneminus
        specmult_planck = 8.0 * specparm_planck
        jpl = 1 + int(specmult_planck)
        fpl = specmult_planck - float(int(specmult_planck))

        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa13 + js
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa13 + js1
        inds = indself[lay - 1]
        indf = indfor[lay - 1]
        indm = indminor[lay - 1]
        fac200 = 0.0
        fac210 = 0.0
        fac201 = 0.0
        fac211 = 0.0

        if specparm < 0.125:
            p = fs - 1.0
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac000 = fk0 * fac00[lay - 1]
            fac100 = fk1 * fac00[lay - 1]
            fac200 = fk2 * fac00[lay - 1]
            fac010 = fk0 * fac10[lay - 1]
            fac110 = fk1 * fac10[lay - 1]
            fac210 = fk2 * fac10[lay - 1]
        elif specparm > 0.875:
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac000 = fk0 * fac00[lay - 1]
            fac100 = fk1 * fac00[lay - 1]
            fac200 = fk2 * fac00[lay - 1]
            fac010 = fk0 * fac10[lay - 1]
            fac110 = fk1 * fac10[lay - 1]
            fac210 = fk2 * fac10[lay - 1]
        else:
            fac000 = (1.0 - fs) * fac00[lay - 1]
            fac010 = (1.0 - fs) * fac10[lay - 1]
            fac100 = fs * fac00[lay - 1]
            fac110 = fs * fac10[lay - 1]

        if specparm1 < 0.125:
            p = fs1 - 1.0
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac001 = fk0 * fac01[lay - 1]
            fac101 = fk1 * fac01[lay - 1]
            fac201 = fk2 * fac01[lay - 1]
            fac011 = fk0 * fac11[lay - 1]
            fac111 = fk1 * fac11[lay - 1]
            fac211 = fk2 * fac11[lay - 1]
        elif specparm1 > 0.875:
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac001 = fk0 * fac01[lay - 1]
            fac101 = fk1 * fac01[lay - 1]
            fac201 = fk2 * fac01[lay - 1]
            fac011 = fk0 * fac11[lay - 1]
            fac111 = fk1 * fac11[lay - 1]
            fac211 = fk2 * fac11[lay - 1]
        else:
            fac001 = (1.0 - fs1) * fac01[lay - 1]
            fac011 = (1.0 - fs1) * fac11[lay - 1]
            fac101 = fs1 * fac01[lay - 1]
            fac111 = fs1 * fac11[lay - 1]

        for ig in range(1, ng13 + 1):
            tauself = selffac[lay - 1] * (
                selfref[_idx2(inds, ig, 10)]
                + selffrac[lay - 1]
                * (selfref[_idx2(inds + 1, ig, 10)] - selfref[_idx2(inds, ig, 10)])
            )
            taufor = forfac[lay - 1] * (
                forref[_idx2(indf, ig, 4)]
                + forfrac[lay - 1]
                * (forref[_idx2(indf + 1, ig, 4)] - forref[_idx2(indf, ig, 4)])
            )
            co2m1 = ka_mco2[_idx3(jmco2, indm, ig, 9, 19)] + fmco2 * (
                ka_mco2[_idx3(jmco2 + 1, indm, ig, 9, 19)]
                - ka_mco2[_idx3(jmco2, indm, ig, 9, 19)]
            )
            co2m2 = ka_mco2[_idx3(jmco2, indm + 1, ig, 9, 19)] + fmco2 * (
                ka_mco2[_idx3(jmco2 + 1, indm + 1, ig, 9, 19)]
                - ka_mco2[_idx3(jmco2, indm + 1, ig, 9, 19)]
            )
            absco2 = co2m1 + minorfrac[lay - 1] * (co2m2 - co2m1)
            com1 = ka_mco[_idx3(jmco, indm, ig, 9, 19)] + fmco * (
                ka_mco[_idx3(jmco + 1, indm, ig, 9, 19)]
                - ka_mco[_idx3(jmco, indm, ig, 9, 19)]
            )
            com2 = ka_mco[_idx3(jmco, indm + 1, ig, 9, 19)] + fmco * (
                ka_mco[_idx3(jmco + 1, indm + 1, ig, 9, 19)]
                - ka_mco[_idx3(jmco, indm + 1, ig, 9, 19)]
            )
            absco = com1 + minorfrac[lay - 1] * (com2 - com1)
            tau_major = _lw_tau_major_3pt(
                specparm, speccomb, ind0, ig, fac000, fac100, fac200,
                fac010, fac110, fac210, absa, 585
            )
            tau_major1 = _lw_tau_major_3pt(
                specparm1, speccomb1, ind1, ig, fac001, fac101, fac201,
                fac011, fac111, fac211, absa, 585
            )
            taug[_idx2(lay, ngs12 + ig, nlayers)] = (
                tau_major
                + tau_major1
                + tauself
                + taufor
                + adjcolco2 * absco2
                + colco[lay - 1] * absco
            )
            fracs[_idx2(lay, ngs12 + ig, nlayers)] = fracrefa[
                _idx2(ig, jpl, ng13)
            ] + fpl * (
                fracrefa[_idx2(ig, jpl + 1, ng13)] - fracrefa[_idx2(ig, jpl, ng13)]
            )

    for lay in range(laytrop + 1, nlayers + 1):
        indm = indminor[lay - 1]
        for ig in range(1, ng13 + 1):
            abso3 = kb_mo3[_idx2(indm, ig, 19)] + minorfrac[lay - 1] * (
                kb_mo3[_idx2(indm + 1, ig, 19)] - kb_mo3[_idx2(indm, ig, 19)]
            )
            taug[_idx2(lay, ngs12 + ig, nlayers)] = colo3[lay - 1] * abso3
            fracs[_idx2(lay, ngs12 + ig, nlayers)] = fracrefb[ig - 1]


@export
def rrtmg_lw_taugb15_codon(
    nlayers: int,
    laytrop: int,
    ng15: int,
    ngs14: int,
    nspa15: int,
    oneminus: float,
    coln2o_p: cobj,
    colco2_p: cobj,
    colbrd_p: cobj,
    rat_n2oco2_p: cobj,
    rat_n2oco2_1_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    indself_p: cobj,
    indfor_p: cobj,
    indminor_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    selffac_p: cobj,
    selffrac_p: cobj,
    forfac_p: cobj,
    forfrac_p: cobj,
    scaleminor_p: cobj,
    minorfrac_p: cobj,
    chi_mls_p: cobj,
    fracrefa_p: cobj,
    absa_p: cobj,
    ka_mn2_p: cobj,
    selfref_p: cobj,
    forref_p: cobj,
    fracs_p: cobj,
    taug_p: cobj,
):
    coln2o = Ptr[float](coln2o_p)
    colco2 = Ptr[float](colco2_p)
    colbrd = Ptr[float](colbrd_p)
    rat_n2oco2 = Ptr[float](rat_n2oco2_p)
    rat_n2oco2_1 = Ptr[float](rat_n2oco2_1_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    indself = Ptr[int](indself_p)
    indfor = Ptr[int](indfor_p)
    indminor = Ptr[int](indminor_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    selffac = Ptr[float](selffac_p)
    selffrac = Ptr[float](selffrac_p)
    forfac = Ptr[float](forfac_p)
    forfrac = Ptr[float](forfrac_p)
    scaleminor = Ptr[float](scaleminor_p)
    minorfrac = Ptr[float](minorfrac_p)
    chi_mls = Ptr[float](chi_mls_p)
    fracrefa = Ptr[float](fracrefa_p)
    absa = Ptr[float](absa_p)
    ka_mn2 = Ptr[float](ka_mn2_p)
    selfref = Ptr[float](selfref_p)
    forref = Ptr[float](forref_p)
    fracs = Ptr[float](fracs_p)
    taug = Ptr[float](taug_p)

    refrat_planck_a = chi_mls[_idx2(4, 1, 7)] / chi_mls[_idx2(2, 1, 7)]
    refrat_m_a = chi_mls[_idx2(4, 1, 7)] / chi_mls[_idx2(2, 1, 7)]

    for lay in range(1, laytrop + 1):
        speccomb = coln2o[lay - 1] + rat_n2oco2[lay - 1] * colco2[lay - 1]
        specparm = coln2o[lay - 1] / speccomb
        if specparm >= oneminus:
            specparm = oneminus
        specmult = 8.0 * specparm
        js = 1 + int(specmult)
        fs = specmult - float(int(specmult))

        speccomb1 = coln2o[lay - 1] + rat_n2oco2_1[lay - 1] * colco2[lay - 1]
        specparm1 = coln2o[lay - 1] / speccomb1
        if specparm1 >= oneminus:
            specparm1 = oneminus
        specmult1 = 8.0 * specparm1
        js1 = 1 + int(specmult1)
        fs1 = specmult1 - float(int(specmult1))

        speccomb_mn2 = coln2o[lay - 1] + refrat_m_a * colco2[lay - 1]
        specparm_mn2 = coln2o[lay - 1] / speccomb_mn2
        if specparm_mn2 >= oneminus:
            specparm_mn2 = oneminus
        specmult_mn2 = 8.0 * specparm_mn2
        jmn2 = 1 + int(specmult_mn2)
        fmn2 = specmult_mn2 - float(int(specmult_mn2))

        speccomb_planck = coln2o[lay - 1] + refrat_planck_a * colco2[lay - 1]
        specparm_planck = coln2o[lay - 1] / speccomb_planck
        if specparm_planck >= oneminus:
            specparm_planck = oneminus
        specmult_planck = 8.0 * specparm_planck
        jpl = 1 + int(specmult_planck)
        fpl = specmult_planck - float(int(specmult_planck))

        ind0 = ((jp[lay - 1] - 1) * 5 + (jt[lay - 1] - 1)) * nspa15 + js
        ind1 = (jp[lay - 1] * 5 + (jt1[lay - 1] - 1)) * nspa15 + js1
        inds = indself[lay - 1]
        indf = indfor[lay - 1]
        indm = indminor[lay - 1]
        scalen2 = colbrd[lay - 1] * scaleminor[lay - 1]
        fac200 = 0.0
        fac210 = 0.0
        fac201 = 0.0
        fac211 = 0.0

        if specparm < 0.125:
            p = fs - 1.0
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac000 = fk0 * fac00[lay - 1]
            fac100 = fk1 * fac00[lay - 1]
            fac200 = fk2 * fac00[lay - 1]
            fac010 = fk0 * fac10[lay - 1]
            fac110 = fk1 * fac10[lay - 1]
            fac210 = fk2 * fac10[lay - 1]
        elif specparm > 0.875:
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac000 = fk0 * fac00[lay - 1]
            fac100 = fk1 * fac00[lay - 1]
            fac200 = fk2 * fac00[lay - 1]
            fac010 = fk0 * fac10[lay - 1]
            fac110 = fk1 * fac10[lay - 1]
            fac210 = fk2 * fac10[lay - 1]
        else:
            fac000 = (1.0 - fs) * fac00[lay - 1]
            fac010 = (1.0 - fs) * fac10[lay - 1]
            fac100 = fs * fac00[lay - 1]
            fac110 = fs * fac10[lay - 1]

        if specparm1 < 0.125:
            p = fs1 - 1.0
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac001 = fk0 * fac01[lay - 1]
            fac101 = fk1 * fac01[lay - 1]
            fac201 = fk2 * fac01[lay - 1]
            fac011 = fk0 * fac11[lay - 1]
            fac111 = fk1 * fac11[lay - 1]
            fac211 = fk2 * fac11[lay - 1]
        elif specparm1 > 0.875:
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1.0 - p - 2.0 * p4
            fk2 = p + p4
            fac001 = fk0 * fac01[lay - 1]
            fac101 = fk1 * fac01[lay - 1]
            fac201 = fk2 * fac01[lay - 1]
            fac011 = fk0 * fac11[lay - 1]
            fac111 = fk1 * fac11[lay - 1]
            fac211 = fk2 * fac11[lay - 1]
        else:
            fac001 = (1.0 - fs1) * fac01[lay - 1]
            fac011 = (1.0 - fs1) * fac11[lay - 1]
            fac101 = fs1 * fac01[lay - 1]
            fac111 = fs1 * fac11[lay - 1]

        for ig in range(1, ng15 + 1):
            tauself = selffac[lay - 1] * (
                selfref[_idx2(inds, ig, 10)]
                + selffrac[lay - 1]
                * (selfref[_idx2(inds + 1, ig, 10)] - selfref[_idx2(inds, ig, 10)])
            )
            taufor = forfac[lay - 1] * (
                forref[_idx2(indf, ig, 4)]
                + forfrac[lay - 1]
                * (forref[_idx2(indf + 1, ig, 4)] - forref[_idx2(indf, ig, 4)])
            )
            n2m1 = ka_mn2[_idx3(jmn2, indm, ig, 9, 19)] + fmn2 * (
                ka_mn2[_idx3(jmn2 + 1, indm, ig, 9, 19)]
                - ka_mn2[_idx3(jmn2, indm, ig, 9, 19)]
            )
            n2m2 = ka_mn2[_idx3(jmn2, indm + 1, ig, 9, 19)] + fmn2 * (
                ka_mn2[_idx3(jmn2 + 1, indm + 1, ig, 9, 19)]
                - ka_mn2[_idx3(jmn2, indm + 1, ig, 9, 19)]
            )
            taun2 = scalen2 * (n2m1 + minorfrac[lay - 1] * (n2m2 - n2m1))
            tau_major = _lw_tau_major_3pt(
                specparm, speccomb, ind0, ig, fac000, fac100, fac200,
                fac010, fac110, fac210, absa, 585
            )
            tau_major1 = _lw_tau_major_3pt(
                specparm1, speccomb1, ind1, ig, fac001, fac101, fac201,
                fac011, fac111, fac211, absa, 585
            )
            taug[_idx2(lay, ngs14 + ig, nlayers)] = (
                tau_major + tau_major1 + tauself + taufor + taun2
            )
            fracs[_idx2(lay, ngs14 + ig, nlayers)] = fracrefa[
                _idx2(ig, jpl, ng15)
            ] + fpl * (
                fracrefa[_idx2(ig, jpl + 1, ng15)] - fracrefa[_idx2(ig, jpl, ng15)]
            )

    for lay in range(laytrop + 1, nlayers + 1):
        for ig in range(1, ng15 + 1):
            taug[_idx2(lay, ngs14 + ig, nlayers)] = 0.0
            fracs[_idx2(lay, ngs14 + ig, nlayers)] = 0.0


@export
def rrtmg_lw_setcoef_codon(
    nlayers: int,
    istart: int,
    mxmol: int,
    nbndlw: int,
    pavel_p: cobj,
    tavel_p: cobj,
    tz_p: cobj,
    tbound: float,
    semiss_p: cobj,
    coldry_p: cobj,
    wkl_p: cobj,
    wbroad_p: cobj,
    laytrop_p: cobj,
    jp_p: cobj,
    jt_p: cobj,
    jt1_p: cobj,
    planklay_p: cobj,
    planklev_p: cobj,
    plankbnd_p: cobj,
    colh2o_p: cobj,
    colco2_p: cobj,
    colo3_p: cobj,
    coln2o_p: cobj,
    colco_p: cobj,
    colch4_p: cobj,
    colo2_p: cobj,
    colbrd_p: cobj,
    fac00_p: cobj,
    fac01_p: cobj,
    fac10_p: cobj,
    fac11_p: cobj,
    rat_h2oco2_p: cobj,
    rat_h2oco2_1_p: cobj,
    rat_h2oo3_p: cobj,
    rat_h2oo3_1_p: cobj,
    rat_h2on2o_p: cobj,
    rat_h2on2o_1_p: cobj,
    rat_h2och4_p: cobj,
    rat_h2och4_1_p: cobj,
    rat_n2oco2_p: cobj,
    rat_n2oco2_1_p: cobj,
    rat_o3co2_p: cobj,
    rat_o3co2_1_p: cobj,
    selffac_p: cobj,
    selffrac_p: cobj,
    indself_p: cobj,
    forfac_p: cobj,
    forfrac_p: cobj,
    indfor_p: cobj,
    minorfrac_p: cobj,
    scaleminor_p: cobj,
    scaleminorn2_p: cobj,
    indminor_p: cobj,
    totplnk_p: cobj,
    totplk16_p: cobj,
    preflog_p: cobj,
    tref_p: cobj,
    chi_mls_p: cobj,
):
    pavel = Ptr[float](pavel_p)
    tavel = Ptr[float](tavel_p)
    tz = Ptr[float](tz_p)
    semiss = Ptr[float](semiss_p)
    coldry = Ptr[float](coldry_p)
    wkl = Ptr[float](wkl_p)
    wbroad = Ptr[float](wbroad_p)
    laytrop = Ptr[int](laytrop_p)
    jp = Ptr[int](jp_p)
    jt = Ptr[int](jt_p)
    jt1 = Ptr[int](jt1_p)
    planklay = Ptr[float](planklay_p)
    planklev = Ptr[float](planklev_p)
    plankbnd = Ptr[float](plankbnd_p)
    colh2o = Ptr[float](colh2o_p)
    colco2 = Ptr[float](colco2_p)
    colo3 = Ptr[float](colo3_p)
    coln2o = Ptr[float](coln2o_p)
    colco = Ptr[float](colco_p)
    colch4 = Ptr[float](colch4_p)
    colo2 = Ptr[float](colo2_p)
    colbrd = Ptr[float](colbrd_p)
    fac00 = Ptr[float](fac00_p)
    fac01 = Ptr[float](fac01_p)
    fac10 = Ptr[float](fac10_p)
    fac11 = Ptr[float](fac11_p)
    rat_h2oco2 = Ptr[float](rat_h2oco2_p)
    rat_h2oco2_1 = Ptr[float](rat_h2oco2_1_p)
    rat_h2oo3 = Ptr[float](rat_h2oo3_p)
    rat_h2oo3_1 = Ptr[float](rat_h2oo3_1_p)
    rat_h2on2o = Ptr[float](rat_h2on2o_p)
    rat_h2on2o_1 = Ptr[float](rat_h2on2o_1_p)
    rat_h2och4 = Ptr[float](rat_h2och4_p)
    rat_h2och4_1 = Ptr[float](rat_h2och4_1_p)
    rat_n2oco2 = Ptr[float](rat_n2oco2_p)
    rat_n2oco2_1 = Ptr[float](rat_n2oco2_1_p)
    rat_o3co2 = Ptr[float](rat_o3co2_p)
    rat_o3co2_1 = Ptr[float](rat_o3co2_1_p)
    selffac = Ptr[float](selffac_p)
    selffrac = Ptr[float](selffrac_p)
    indself = Ptr[int](indself_p)
    forfac = Ptr[float](forfac_p)
    forfrac = Ptr[float](forfrac_p)
    indfor = Ptr[int](indfor_p)
    minorfrac = Ptr[float](minorfrac_p)
    scaleminor = Ptr[float](scaleminor_p)
    scaleminorn2 = Ptr[float](scaleminorn2_p)
    indminor = Ptr[int](indminor_p)
    totplnk = Ptr[float](totplnk_p)
    totplk16 = Ptr[float](totplk16_p)
    preflog = Ptr[float](preflog_p)
    tref = Ptr[float](tref_p)
    chi_mls = Ptr[float](chi_mls_p)

    stpfac = 296.0 / 1013.0

    indbound = int(tbound - 159.0)
    if indbound < 1:
        indbound = 1
    elif indbound > 180:
        indbound = 180
    tbndfrac = tbound - 159.0 - float(indbound)

    indlev0 = int(tz[0] - 159.0)
    if indlev0 < 1:
        indlev0 = 1
    elif indlev0 > 180:
        indlev0 = 180
    t0frac = tz[0] - 159.0 - float(indlev0)

    laytrop[0] = 0

    for lay in range(1, nlayers + 1):
        idx = lay - 1
        indlay = int(tavel[idx] - 159.0)
        if indlay < 1:
            indlay = 1
        elif indlay > 180:
            indlay = 180
        tlayfrac = tavel[idx] - 159.0 - float(indlay)

        indlev = int(tz[lay] - 159.0)
        if indlev < 1:
            indlev = 1
        elif indlev > 180:
            indlev = 180
        tlevfrac = tz[lay] - 159.0 - float(indlev)

        for iband in range(1, nbndlw):
            if lay == 1:
                dbdtlev = (
                    totplnk[_idx2(indbound + 1, iband, 181)]
                    - totplnk[_idx2(indbound, iband, 181)]
                )
                plankbnd[iband - 1] = semiss[iband - 1] * (
                    totplnk[_idx2(indbound, iband, 181)] + tbndfrac * dbdtlev
                )
                dbdtlev = (
                    totplnk[_idx2(indlev0 + 1, iband, 181)]
                    - totplnk[_idx2(indlev0, iband, 181)]
                )
                planklev[_idx2_dim1_lb0(0, iband, nlayers)] = (
                    totplnk[_idx2(indlev0, iband, 181)] + t0frac * dbdtlev
                )
            dbdtlev = (
                totplnk[_idx2(indlev + 1, iband, 181)]
                - totplnk[_idx2(indlev, iband, 181)]
            )
            dbdtlay = (
                totplnk[_idx2(indlay + 1, iband, 181)]
                - totplnk[_idx2(indlay, iband, 181)]
            )
            planklay[_idx2(lay, iband, nlayers)] = (
                totplnk[_idx2(indlay, iband, 181)] + tlayfrac * dbdtlay
            )
            planklev[_idx2_dim1_lb0(lay, iband, nlayers)] = (
                totplnk[_idx2(indlev, iband, 181)] + tlevfrac * dbdtlev
            )

        iband = nbndlw
        if istart == nbndlw:
            if lay == 1:
                dbdtlev = totplk16[indbound] - totplk16[indbound - 1]
                plankbnd[iband - 1] = semiss[iband - 1] * (
                    totplk16[indbound - 1] + tbndfrac * dbdtlev
                )
                dbdtlev = (
                    totplnk[_idx2(indlev0 + 1, iband, 181)]
                    - totplnk[_idx2(indlev0, iband, 181)]
                )
                planklev[_idx2_dim1_lb0(0, iband, nlayers)] = (
                    totplk16[indlev0 - 1] + t0frac * dbdtlev
                )
            dbdtlev = totplk16[indlev] - totplk16[indlev - 1]
            dbdtlay = totplk16[indlay] - totplk16[indlay - 1]
            planklay[_idx2(lay, iband, nlayers)] = (
                totplk16[indlay - 1] + tlayfrac * dbdtlay
            )
            planklev[_idx2_dim1_lb0(lay, iband, nlayers)] = (
                totplk16[indlev - 1] + tlevfrac * dbdtlev
            )
        else:
            if lay == 1:
                dbdtlev = (
                    totplnk[_idx2(indbound + 1, iband, 181)]
                    - totplnk[_idx2(indbound, iband, 181)]
                )
                plankbnd[iband - 1] = semiss[iband - 1] * (
                    totplnk[_idx2(indbound, iband, 181)] + tbndfrac * dbdtlev
                )
                dbdtlev = (
                    totplnk[_idx2(indlev0 + 1, iband, 181)]
                    - totplnk[_idx2(indlev0, iband, 181)]
                )
                planklev[_idx2_dim1_lb0(0, iband, nlayers)] = (
                    totplnk[_idx2(indlev0, iband, 181)] + t0frac * dbdtlev
                )
            dbdtlev = (
                totplnk[_idx2(indlev + 1, iband, 181)]
                - totplnk[_idx2(indlev, iband, 181)]
            )
            dbdtlay = (
                totplnk[_idx2(indlay + 1, iband, 181)]
                - totplnk[_idx2(indlay, iband, 181)]
            )
            planklay[_idx2(lay, iband, nlayers)] = (
                totplnk[_idx2(indlay, iband, 181)] + tlayfrac * dbdtlay
            )
            planklev[_idx2_dim1_lb0(lay, iband, nlayers)] = (
                totplnk[_idx2(indlev, iband, 181)] + tlevfrac * dbdtlev
            )

        plog = log(pavel[idx])
        jp[idx] = int(36.0 - 5.0 * (plog + 0.04))
        if jp[idx] < 1:
            jp[idx] = 1
        elif jp[idx] > 58:
            jp[idx] = 58
        jp1 = jp[idx] + 1
        fp = 5.0 * (preflog[jp[idx] - 1] - plog)

        jt[idx] = int(3.0 + (tavel[idx] - tref[jp[idx] - 1]) / 15.0)
        if jt[idx] < 1:
            jt[idx] = 1
        elif jt[idx] > 4:
            jt[idx] = 4
        ft = ((tavel[idx] - tref[jp[idx] - 1]) / 15.0) - float(jt[idx] - 3)

        jt1[idx] = int(3.0 + (tavel[idx] - tref[jp1 - 1]) / 15.0)
        if jt1[idx] < 1:
            jt1[idx] = 1
        elif jt1[idx] > 4:
            jt1[idx] = 4
        ft1 = ((tavel[idx] - tref[jp1 - 1]) / 15.0) - float(jt1[idx] - 3)

        water = wkl[_idx2(1, lay, mxmol)] / coldry[idx]
        scalefac = pavel[idx] * stpfac / tavel[idx]

        if plog > 4.56:
            laytrop[0] = laytrop[0] + 1

            forfac[idx] = scalefac / (1.0 + water)
            factor = (332.0 - tavel[idx]) / 36.0
            ind = int(factor)
            if ind < 1:
                ind = 1
            elif ind > 2:
                ind = 2
            indfor[idx] = ind
            forfrac[idx] = factor - float(indfor[idx])

            selffac[idx] = water * forfac[idx]
            factor = (tavel[idx] - 188.0) / 7.2
            ind = int(factor) - 7
            if ind < 1:
                ind = 1
            elif ind > 9:
                ind = 9
            indself[idx] = ind
            selffrac[idx] = factor - float(indself[idx] + 7)

            scaleminor[idx] = pavel[idx] / tavel[idx]
            scaleminorn2[idx] = (pavel[idx] / tavel[idx]) * (
                wbroad[idx] / (coldry[idx] + wkl[_idx2(1, lay, mxmol)])
            )
            factor = (tavel[idx] - 180.8) / 7.2
            ind = int(factor)
            if ind < 1:
                ind = 1
            elif ind > 18:
                ind = 18
            indminor[idx] = ind
            minorfrac[idx] = factor - float(indminor[idx])

            rat_h2oco2[idx] = (
                chi_mls[_idx2(1, jp[idx], 7)] / chi_mls[_idx2(2, jp[idx], 7)]
            )
            rat_h2oco2_1[idx] = (
                chi_mls[_idx2(1, jp[idx] + 1, 7)]
                / chi_mls[_idx2(2, jp[idx] + 1, 7)]
            )
            rat_h2oo3[idx] = (
                chi_mls[_idx2(1, jp[idx], 7)] / chi_mls[_idx2(3, jp[idx], 7)]
            )
            rat_h2oo3_1[idx] = (
                chi_mls[_idx2(1, jp[idx] + 1, 7)]
                / chi_mls[_idx2(3, jp[idx] + 1, 7)]
            )
            rat_h2on2o[idx] = (
                chi_mls[_idx2(1, jp[idx], 7)] / chi_mls[_idx2(4, jp[idx], 7)]
            )
            rat_h2on2o_1[idx] = (
                chi_mls[_idx2(1, jp[idx] + 1, 7)]
                / chi_mls[_idx2(4, jp[idx] + 1, 7)]
            )
            rat_h2och4[idx] = (
                chi_mls[_idx2(1, jp[idx], 7)] / chi_mls[_idx2(6, jp[idx], 7)]
            )
            rat_h2och4_1[idx] = (
                chi_mls[_idx2(1, jp[idx] + 1, 7)]
                / chi_mls[_idx2(6, jp[idx] + 1, 7)]
            )
            rat_n2oco2[idx] = (
                chi_mls[_idx2(4, jp[idx], 7)] / chi_mls[_idx2(2, jp[idx], 7)]
            )
            rat_n2oco2_1[idx] = (
                chi_mls[_idx2(4, jp[idx] + 1, 7)]
                / chi_mls[_idx2(2, jp[idx] + 1, 7)]
            )
        else:
            forfac[idx] = scalefac / (1.0 + water)
            factor = (tavel[idx] - 188.0) / 36.0
            indfor[idx] = 3
            forfrac[idx] = factor - 1.0

            selffac[idx] = water * forfac[idx]

            scaleminor[idx] = pavel[idx] / tavel[idx]
            scaleminorn2[idx] = (pavel[idx] / tavel[idx]) * (
                wbroad[idx] / (coldry[idx] + wkl[_idx2(1, lay, mxmol)])
            )
            factor = (tavel[idx] - 180.8) / 7.2
            ind = int(factor)
            if ind < 1:
                ind = 1
            elif ind > 18:
                ind = 18
            indminor[idx] = ind
            minorfrac[idx] = factor - float(indminor[idx])

            rat_h2oco2[idx] = (
                chi_mls[_idx2(1, jp[idx], 7)] / chi_mls[_idx2(2, jp[idx], 7)]
            )
            rat_h2oco2_1[idx] = (
                chi_mls[_idx2(1, jp[idx] + 1, 7)]
                / chi_mls[_idx2(2, jp[idx] + 1, 7)]
            )
            rat_o3co2[idx] = (
                chi_mls[_idx2(3, jp[idx], 7)] / chi_mls[_idx2(2, jp[idx], 7)]
            )
            rat_o3co2_1[idx] = (
                chi_mls[_idx2(3, jp[idx] + 1, 7)]
                / chi_mls[_idx2(2, jp[idx] + 1, 7)]
            )

        colh2o[idx] = 1.0e-20 * wkl[_idx2(1, lay, mxmol)]
        colco2[idx] = 1.0e-20 * wkl[_idx2(2, lay, mxmol)]
        colo3[idx] = 1.0e-20 * wkl[_idx2(3, lay, mxmol)]
        coln2o[idx] = 1.0e-20 * wkl[_idx2(4, lay, mxmol)]
        colco[idx] = 1.0e-20 * wkl[_idx2(5, lay, mxmol)]
        colch4[idx] = 1.0e-20 * wkl[_idx2(6, lay, mxmol)]
        colo2[idx] = 1.0e-20 * wkl[_idx2(7, lay, mxmol)]
        if colco2[idx] == 0.0:
            colco2[idx] = 1.0e-32 * coldry[idx]
        if colo3[idx] == 0.0:
            colo3[idx] = 1.0e-32 * coldry[idx]
        if coln2o[idx] == 0.0:
            coln2o[idx] = 1.0e-32 * coldry[idx]
        if colco[idx] == 0.0:
            colco[idx] = 1.0e-32 * coldry[idx]
        if colch4[idx] == 0.0:
            colch4[idx] = 1.0e-32 * coldry[idx]
        colbrd[idx] = 1.0e-20 * wbroad[idx]

        compfp = 1.0 - fp
        fac10[idx] = compfp * ft
        fac00[idx] = compfp * (1.0 - ft)
        fac11[idx] = fp * ft1
        fac01[idx] = fp * (1.0 - ft1)

        selffac[idx] = colh2o[idx] * selffac[idx]
        forfac[idx] = colh2o[idx] * forfac[idx]


@export
def rrtmg_lw_cldprmc_codon(
    nlayers: int,
    inflag: int,
    iceflag: int,
    liqflag: int,
    ngptlw: int,
    absliq0: float,
    cldfmc_p: cobj,
    ciwpmc_p: cobj,
    clwpmc_p: cobj,
    reicmc_p: cobj,
    dgesmc_p: cobj,
    relqmc_p: cobj,
    ncbands_p: cobj,
    taucmc_p: cobj,
    absice0_p: cobj,
    absice1_p: cobj,
    absice2_p: cobj,
    absice3_p: cobj,
    absliq1_p: cobj,
    ngb_p: cobj,
):
    cldfmc = Ptr[float](cldfmc_p)
    ciwpmc = Ptr[float](ciwpmc_p)
    clwpmc = Ptr[float](clwpmc_p)
    reicmc = Ptr[float](reicmc_p)
    dgesmc = Ptr[float](dgesmc_p)
    relqmc = Ptr[float](relqmc_p)
    ncbands = Ptr[int](ncbands_p)
    taucmc = Ptr[float](taucmc_p)
    absice0 = Ptr[float](absice0_p)
    absice1 = Ptr[float](absice1_p)
    absice2 = Ptr[float](absice2_p)
    absice3 = Ptr[float](absice3_p)
    absliq1 = Ptr[float](absliq1_p)
    ngb = Ptr[int](ngb_p)

    cldmin = 1.0e-80
    ncbands[0] = 1

    for lay in range(1, nlayers + 1):
        for ig in range(1, ngptlw + 1):
            mc_idx = _idx2(ig, lay, ngptlw)
            cwp = ciwpmc[mc_idx] + clwpmc[mc_idx]
            if cldfmc[mc_idx] >= cldmin and (
                cwp >= cldmin or taucmc[mc_idx] >= cldmin
            ):
                if inflag == 0:
                    return
                elif inflag == 2:
                    radice = reicmc[lay - 1]

                    if ciwpmc[mc_idx] == 0.0:
                        abscoice = 0.0
                    elif iceflag == 0:
                        abscoice = absice0[0] + absice0[1] / radice
                    elif iceflag == 1:
                        ncbands[0] = 5
                        ib = ngb[ig - 1]
                        abscoice = (
                            absice1[_idx2(1, ib, 2)]
                            + absice1[_idx2(2, ib, 2)] / radice
                        )
                    elif iceflag == 2:
                        if radice >= 5.0 and radice <= 131.0:
                            ncbands[0] = 16
                            factor = (radice - 2.0) / 3.0
                            index = int(factor)
                            if index == 43:
                                index = 42
                            fint = factor - float(index)
                            ib = ngb[ig - 1]
                            abscoice = absice2[_idx2(index, ib, 43)] + fint * (
                                absice2[_idx2(index + 1, ib, 43)]
                                - (absice2[_idx2(index, ib, 43)])
                            )
                        elif radice > 131.0:
                            abscoice = absice0[0] + absice0[1] / radice
                        else:
                            abscoice = 0.0
                    elif iceflag == 3:
                        dgeice = dgesmc[lay - 1]
                        if dgeice >= 5.0 and dgeice <= 140.0:
                            ncbands[0] = 16
                            factor = (dgeice - 2.0) / 3.0
                            index = int(factor)
                            if index == 46:
                                index = 45
                            fint = factor - float(index)
                            ib = ngb[ig - 1]
                            abscoice = absice3[_idx2(index, ib, 46)] + fint * (
                                absice3[_idx2(index + 1, ib, 46)]
                                - (absice3[_idx2(index, ib, 46)])
                            )
                        elif dgeice > 140.0:
                            abscoice = absice0[0] + absice0[1] / radice
                        else:
                            abscoice = 0.0
                    else:
                        abscoice = 0.0

                    if clwpmc[mc_idx] == 0.0:
                        abscoliq = 0.0
                    elif liqflag == 0:
                        abscoliq = absliq0
                    elif liqflag == 1:
                        radliq = relqmc[lay - 1]
                        index = int(radliq - 1.5)
                        if index == 58:
                            index = 57
                        if index == 0:
                            index = 1
                        fint = radliq - 1.5 - float(index)
                        ib = ngb[ig - 1]
                        abscoliq = absliq1[_idx2(index, ib, 58)] + fint * (
                            absliq1[_idx2(index + 1, ib, 58)]
                            - (absliq1[_idx2(index, ib, 58)])
                        )
                    else:
                        abscoliq = 0.0

                    taucmc[mc_idx] = (
                        ciwpmc[mc_idx] * abscoice + clwpmc[mc_idx] * abscoliq
                    )


@inline
def _rrtmg_sw_ec_icx(wavenum2: float) -> int:
    if wavenum2 > 1.43e04:
        return 1
    if wavenum2 > 7.7e03:
        return 2
    if wavenum2 > 5.3e03:
        return 3
    if wavenum2 > 4.0e03:
        return 4
    return 5


@export
def rrtmg_sw_cldprmc_codon(
    nlayers: int,
    inflag: int,
    iceflag: int,
    liqflag: int,
    ngptsw: int,
    jpb1: int,
    cldfmc_p: cobj,
    ciwpmc_p: cobj,
    clwpmc_p: cobj,
    reicmc_p: cobj,
    dgesmc_p: cobj,
    relqmc_p: cobj,
    taormc_p: cobj,
    taucmc_p: cobj,
    ssacmc_p: cobj,
    asmcmc_p: cobj,
    fsfcmc_p: cobj,
    extliq1_p: cobj,
    ssaliq1_p: cobj,
    asyliq1_p: cobj,
    extice2_p: cobj,
    ssaice2_p: cobj,
    asyice2_p: cobj,
    extice3_p: cobj,
    ssaice3_p: cobj,
    asyice3_p: cobj,
    fdlice3_p: cobj,
    abari_p: cobj,
    bbari_p: cobj,
    cbari_p: cobj,
    dbari_p: cobj,
    ebari_p: cobj,
    fbari_p: cobj,
    wavenum2_p: cobj,
    ngb_p: cobj,
):
    cldfmc = Ptr[float](cldfmc_p)
    ciwpmc = Ptr[float](ciwpmc_p)
    clwpmc = Ptr[float](clwpmc_p)
    reicmc = Ptr[float](reicmc_p)
    dgesmc = Ptr[float](dgesmc_p)
    relqmc = Ptr[float](relqmc_p)
    taormc = Ptr[float](taormc_p)
    taucmc = Ptr[float](taucmc_p)
    ssacmc = Ptr[float](ssacmc_p)
    asmcmc = Ptr[float](asmcmc_p)
    fsfcmc = Ptr[float](fsfcmc_p)
    extliq1 = Ptr[float](extliq1_p)
    ssaliq1 = Ptr[float](ssaliq1_p)
    asyliq1 = Ptr[float](asyliq1_p)
    extice2 = Ptr[float](extice2_p)
    ssaice2 = Ptr[float](ssaice2_p)
    asyice2 = Ptr[float](asyice2_p)
    extice3 = Ptr[float](extice3_p)
    ssaice3 = Ptr[float](ssaice3_p)
    asyice3 = Ptr[float](asyice3_p)
    fdlice3 = Ptr[float](fdlice3_p)
    abari = Ptr[float](abari_p)
    bbari = Ptr[float](bbari_p)
    cbari = Ptr[float](cbari_p)
    dbari = Ptr[float](dbari_p)
    ebari = Ptr[float](ebari_p)
    fbari = Ptr[float](fbari_p)
    wavenum2 = Ptr[float](wavenum2_p)
    ngb = Ptr[int](ngb_p)

    eps = 1.0e-06
    cldmin = 1.0e-80

    for lay in range(1, nlayers + 1):
        for ig in range(1, ngptsw + 1):
            mc_idx = _idx2(ig, lay, ngptsw)
            taormc[mc_idx] = taucmc[mc_idx]

    for lay in range(1, nlayers + 1):
        for ig in range(1, ngptsw + 1):
            mc_idx = _idx2(ig, lay, ngptsw)
            cwp = ciwpmc[mc_idx] + clwpmc[mc_idx]
            if cldfmc[mc_idx] >= cldmin and (
                cwp >= cldmin or taucmc[mc_idx] >= cldmin
            ):
                if inflag == 0:
                    taucldorig_a = taucmc[mc_idx]
                    ffp = fsfcmc[mc_idx]
                    ffp1 = 1.0 - ffp
                    ffpssa = 1.0 - ffp * ssacmc[mc_idx]
                    ssacloud_a = ffp1 * ssacmc[mc_idx] / ffpssa
                    taucloud_a = ffpssa * taucldorig_a

                    taormc[mc_idx] = taucldorig_a
                    ssacmc[mc_idx] = ssacloud_a
                    taucmc[mc_idx] = taucloud_a
                    asmcmc[mc_idx] = (asmcmc[mc_idx] - ffp) / (ffp1)
                elif inflag == 2:
                    radice = reicmc[lay - 1]
                    extcoice = 0.0
                    ssacoice = 0.0
                    gice = 0.0
                    forwice = 0.0

                    if ciwpmc[mc_idx] == 0.0:
                        extcoice = 0.0
                        ssacoice = 0.0
                        gice = 0.0
                        forwice = 0.0
                    elif iceflag == 1:
                        ib = ngb[ig - 1]
                        icx = _rrtmg_sw_ec_icx(wavenum2[ib - jpb1])
                        extcoice = abari[icx - 1] + bbari[icx - 1] / radice
                        ssacoice = 1.0 - cbari[icx - 1] - dbari[icx - 1] * radice
                        gice = ebari[icx - 1] + fbari[icx - 1] * radice
                        if gice >= 1.0:
                            gice = 1.0 - eps
                        forwice = gice * gice
                    elif iceflag == 2:
                        if radice >= 5.0 and radice <= 131.0:
                            factor = (radice - 2.0) / 3.0
                            index = int(factor)
                            if index == 43:
                                index = 42
                            fint = factor - float(index)
                            ib = ngb[ig - 1]
                            extcoice = extice2[
                                _idx2_dim2_lb(index, ib, 43, jpb1)
                            ] + fint * (
                                extice2[_idx2_dim2_lb(index + 1, ib, 43, jpb1)]
                                - extice2[_idx2_dim2_lb(index, ib, 43, jpb1)]
                            )
                            ssacoice = ssaice2[
                                _idx2_dim2_lb(index, ib, 43, jpb1)
                            ] + fint * (
                                ssaice2[_idx2_dim2_lb(index + 1, ib, 43, jpb1)]
                                - ssaice2[_idx2_dim2_lb(index, ib, 43, jpb1)]
                            )
                            gice = asyice2[
                                _idx2_dim2_lb(index, ib, 43, jpb1)
                            ] + fint * (
                                asyice2[_idx2_dim2_lb(index + 1, ib, 43, jpb1)]
                                - asyice2[_idx2_dim2_lb(index, ib, 43, jpb1)]
                            )
                            forwice = gice * gice
                        elif radice > 131.0:
                            ib = ngb[ig - 1]
                            icx = _rrtmg_sw_ec_icx(wavenum2[ib - jpb1])
                            extcoice = abari[icx - 1] + bbari[icx - 1] / radice
                            ssacoice = 1.0 - cbari[icx - 1] - dbari[icx - 1] * radice
                            gice = ebari[icx - 1] + fbari[icx - 1] * radice
                            if gice >= 1.0:
                                gice = 1.0 - eps
                            forwice = gice * gice
                    elif iceflag == 3:
                        dgeice = dgesmc[lay - 1]
                        if dgeice >= 5.0 and dgeice <= 140.0:
                            factor = (dgeice - 2.0) / 3.0
                            index = int(factor)
                            if index == 46:
                                index = 45
                            fint = factor - float(index)
                            ib = ngb[ig - 1]
                            extcoice = extice3[
                                _idx2_dim2_lb(index, ib, 46, jpb1)
                            ] + fint * (
                                extice3[_idx2_dim2_lb(index + 1, ib, 46, jpb1)]
                                - extice3[_idx2_dim2_lb(index, ib, 46, jpb1)]
                            )
                            ssacoice = ssaice3[
                                _idx2_dim2_lb(index, ib, 46, jpb1)
                            ] + fint * (
                                ssaice3[_idx2_dim2_lb(index + 1, ib, 46, jpb1)]
                                - ssaice3[_idx2_dim2_lb(index, ib, 46, jpb1)]
                            )
                            gice = asyice3[
                                _idx2_dim2_lb(index, ib, 46, jpb1)
                            ] + fint * (
                                asyice3[_idx2_dim2_lb(index + 1, ib, 46, jpb1)]
                                - asyice3[_idx2_dim2_lb(index, ib, 46, jpb1)]
                            )
                            fdelta = fdlice3[
                                _idx2_dim2_lb(index, ib, 46, jpb1)
                            ] + fint * (
                                fdlice3[_idx2_dim2_lb(index + 1, ib, 46, jpb1)]
                                - fdlice3[_idx2_dim2_lb(index, ib, 46, jpb1)]
                            )
                            forwice = fdelta + 0.5 / ssacoice
                            if forwice > gice:
                                forwice = gice
                        elif dgeice > 140.0:
                            ib = ngb[ig - 1]
                            icx = _rrtmg_sw_ec_icx(wavenum2[ib - jpb1])
                            extcoice = abari[icx - 1] + bbari[icx - 1] / radice
                            ssacoice = 1.0 - cbari[icx - 1] - dbari[icx - 1] * radice
                            gice = ebari[icx - 1] + fbari[icx - 1] * radice
                            if gice >= 1.0:
                                gice = 1.0 - eps
                            forwice = gice * gice

                    extcoliq = 0.0
                    ssacoliq = 0.0
                    gliq = 0.0
                    forwliq = 0.0

                    if clwpmc[mc_idx] == 0.0:
                        extcoliq = 0.0
                        ssacoliq = 0.0
                        gliq = 0.0
                        forwliq = 0.0
                    elif liqflag == 1:
                        radliq = relqmc[lay - 1]
                        index = int(radliq - 1.5)
                        if index == 0:
                            index = 1
                        if index == 58:
                            index = 57
                        fint = radliq - 1.5 - float(index)
                        ib = ngb[ig - 1]
                        extcoliq = extliq1[
                            _idx2_dim2_lb(index, ib, 58, jpb1)
                        ] + fint * (
                            extliq1[_idx2_dim2_lb(index + 1, ib, 58, jpb1)]
                            - extliq1[_idx2_dim2_lb(index, ib, 58, jpb1)]
                        )
                        ssacoliq = ssaliq1[
                            _idx2_dim2_lb(index, ib, 58, jpb1)
                        ] + fint * (
                            ssaliq1[_idx2_dim2_lb(index + 1, ib, 58, jpb1)]
                            - ssaliq1[_idx2_dim2_lb(index, ib, 58, jpb1)]
                        )
                        if fint < 0.0 and ssacoliq > 1.0:
                            ssacoliq = ssaliq1[
                                _idx2_dim2_lb(index, ib, 58, jpb1)
                            ]
                        gliq = asyliq1[
                            _idx2_dim2_lb(index, ib, 58, jpb1)
                        ] + fint * (
                            asyliq1[_idx2_dim2_lb(index + 1, ib, 58, jpb1)]
                            - asyliq1[_idx2_dim2_lb(index, ib, 58, jpb1)]
                        )
                        forwliq = gliq * gliq

                    tauliqorig = clwpmc[mc_idx] * extcoliq
                    tauiceorig = ciwpmc[mc_idx] * extcoice
                    taormc[mc_idx] = tauliqorig + tauiceorig

                    ssaliq = ssacoliq * (1.0 - forwliq) / (
                        1.0 - forwliq * ssacoliq
                    )
                    tauliq = (1.0 - forwliq * ssacoliq) * tauliqorig
                    ssaice = ssacoice * (1.0 - forwice) / (
                        1.0 - forwice * ssacoice
                    )
                    tauice = (1.0 - forwice * ssacoice) * tauiceorig

                    scatliq = ssaliq * tauliq
                    scatice = ssaice * tauice
                    taucmc[mc_idx] = tauliq + tauice

                    if taucmc[mc_idx] == 0.0:
                        taucmc[mc_idx] = cldmin
                    if scatice == 0.0:
                        scatice = cldmin

                    ssacmc[mc_idx] = (scatliq + scatice) / taucmc[mc_idx]

                    if iceflag == 3:
                        asmcmc[mc_idx] = (1.0 / (scatliq + scatice)) * (
                            scatliq * (gliq - forwliq) / (1.0 - forwliq)
                            + scatice * ((gice - forwice) / (1.0 - forwice))
                        )
                    else:
                        asmcmc[mc_idx] = (
                            scatliq * (gliq - forwliq) / (1.0 - forwliq)
                            + scatice * (gice - forwice) / (1.0 - forwice)
                        ) / (scatliq + scatice)


@inline
def _rrtmg_sw_exp_lookup(ze1: float, tblint: float, bpade: float, od_lo: float, exp_tbl: Ptr[float]) -> float:
    if ze1 <= od_lo:
        return 1.0 - ze1 + 0.5 * ze1 * ze1
    tblind = ze1 / (bpade + ze1)
    itind = int(tblint * tblind + 0.5)
    return exp_tbl[itind]


@export
def rrtmg_sw_spcvmc_pre_reftra_codon(
    klev: int,
    ngptsw: int,
    nbndsw: int,
    iw: int,
    ibm: int,
    icpr: int,
    idelm: int,
    prmu0: float,
    repclc: float,
    tblint: float,
    bpade: float,
    od_lo: float,
    pcldfmc_p: cobj,
    ptaucmc_p: cobj,
    pasycmc_p: cobj,
    pomgcmc_p: cobj,
    ptaormc_p: cobj,
    ptaua_p: cobj,
    pasya_p: cobj,
    pomga_p: cobj,
    ztaug_p: cobj,
    ztaur_p: cobj,
    lrtchkclr_p: cobj,
    lrtchkcld_p: cobj,
    ztauc_p: cobj,
    zomcc_p: cobj,
    zgcc_p: cobj,
    ztauo_p: cobj,
    zomco_p: cobj,
    zgco_p: cobj,
    zdbtc_nodel_p: cobj,
    ztdbtc_nodel_p: cobj,
    zdbt_nodel_p: cobj,
    ztdbt_nodel_p: cobj,
    exp_tbl_p: cobj,
):
    pcldfmc = Ptr[float](pcldfmc_p)
    ptaucmc = Ptr[float](ptaucmc_p)
    pasycmc = Ptr[float](pasycmc_p)
    pomgcmc = Ptr[float](pomgcmc_p)
    ptaormc = Ptr[float](ptaormc_p)
    ptaua = Ptr[float](ptaua_p)
    pasya = Ptr[float](pasya_p)
    pomga = Ptr[float](pomga_p)
    ztaug = Ptr[float](ztaug_p)
    ztaur = Ptr[float](ztaur_p)
    lrtchkclr = Ptr[int](lrtchkclr_p)
    lrtchkcld = Ptr[int](lrtchkcld_p)
    ztauc = Ptr[float](ztauc_p)
    zomcc = Ptr[float](zomcc_p)
    zgcc = Ptr[float](zgcc_p)
    ztauo = Ptr[float](ztauo_p)
    zomco = Ptr[float](zomco_p)
    zgco = Ptr[float](zgco_p)
    zdbtc_nodel = Ptr[float](zdbtc_nodel_p)
    ztdbtc_nodel = Ptr[float](ztdbtc_nodel_p)
    zdbt_nodel = Ptr[float](zdbt_nodel_p)
    ztdbt_nodel = Ptr[float](ztdbt_nodel_p)
    exp_tbl = Ptr[float](exp_tbl_p)

    for jk in range(1, klev + 1):
        ikl = klev + 1 - jk
        jidx = jk - 1
        mc_idx = _idx2(ikl, iw, klev)
        aer_idx = _idx2(ikl, ibm, klev)

        lrtchkclr[jidx] = 1
        lrtchkcld[jidx] = 0
        if pcldfmc[mc_idx] > repclc:
            lrtchkcld[jidx] = 1

        ztauc[jidx] = ztaur[mc_idx] + ztaug[mc_idx] + ptaua[aer_idx]
        zomcc[jidx] = ztaur[mc_idx] * 1.0 + ptaua[aer_idx] * pomga[aer_idx]
        zgcc[jidx] = pasya[aer_idx] * pomga[aer_idx] * ptaua[aer_idx] / zomcc[jidx]
        zomcc[jidx] = zomcc[jidx] / ztauc[jidx]

        if idelm == 0:
            zclear = 1.0 - pcldfmc[mc_idx]
            zcloud = pcldfmc[mc_idx]

            ze1 = ztauc[jidx] / prmu0
            zdbtmc = _rrtmg_sw_exp_lookup(ze1, tblint, bpade, od_lo, exp_tbl)
            zdbtc_nodel[jidx] = zdbtmc
            ztdbtc_nodel[jk] = zdbtc_nodel[jidx] * ztdbtc_nodel[jidx]

            tauorig = ztauc[jidx] + ptaormc[mc_idx]
            ze1 = tauorig / prmu0
            zdbtmo = _rrtmg_sw_exp_lookup(ze1, tblint, bpade, od_lo, exp_tbl)
            zdbt_nodel[jidx] = zclear * zdbtmc + zcloud * zdbtmo
            ztdbt_nodel[jk] = zdbt_nodel[jidx] * ztdbt_nodel[jidx]

        zf = zgcc[jidx] * zgcc[jidx]
        zwf = zomcc[jidx] * zf
        ztauc[jidx] = (1.0 - zwf) * ztauc[jidx]
        zomcc[jidx] = (zomcc[jidx] - zwf) / (1.0 - zwf)
        zgcc[jidx] = (zgcc[jidx] - zf) / (1.0 - zf)

        if icpr >= 1:
            ztauo[jidx] = ztauc[jidx] + ptaucmc[mc_idx]
            zomco[jidx] = ztauc[jidx] * zomcc[jidx] + ptaucmc[mc_idx] * pomgcmc[mc_idx]
            zgco[jidx] = (
                ptaucmc[mc_idx] * pomgcmc[mc_idx] * pasycmc[mc_idx]
                + ztauc[jidx] * zomcc[jidx] * zgcc[jidx]
            ) / zomco[jidx]
            zomco[jidx] = zomco[jidx] / ztauo[jidx]
        elif icpr == 0:
            ztauo[jidx] = ztaur[mc_idx] + ztaug[mc_idx] + ptaua[aer_idx] + ptaucmc[mc_idx]
            zomco[jidx] = (
                ptaua[aer_idx] * pomga[aer_idx]
                + ptaucmc[mc_idx] * pomgcmc[mc_idx]
                + ztaur[mc_idx] * 1.0
            )
            zgco[jidx] = (
                ptaucmc[mc_idx] * pomgcmc[mc_idx] * pasycmc[mc_idx]
                + ptaua[aer_idx] * pomga[aer_idx] * pasya[aer_idx]
            ) / zomco[jidx]
            zomco[jidx] = zomco[jidx] / ztauo[jidx]

            zf = zgco[jidx] * zgco[jidx]
            zwf = zomco[jidx] * zf
            ztauo[jidx] = (1.0 - zwf) * ztauo[jidx]
            zomco[jidx] = (zomco[jidx] - zwf) / (1.0 - zwf)
            zgco[jidx] = (zgco[jidx] - zf) / (1.0 - zf)


@export
def rrtmg_sw_spcvmc_post_reftra_codon(
    klev: int,
    ngptsw: int,
    iw: int,
    prmu0: float,
    tblint: float,
    bpade: float,
    od_lo: float,
    pcldfmc_p: cobj,
    ztauc_p: cobj,
    ztauo_p: cobj,
    zrefc_p: cobj,
    zrefdc_p: cobj,
    ztrac_p: cobj,
    ztradc_p: cobj,
    zrefo_p: cobj,
    zrefdo_p: cobj,
    ztrao_p: cobj,
    ztrado_p: cobj,
    zref_p: cobj,
    zrefd_p: cobj,
    ztra_p: cobj,
    ztrad_p: cobj,
    zdbtc_p: cobj,
    ztdbtc_p: cobj,
    zdbt_p: cobj,
    ztdbt_p: cobj,
    exp_tbl_p: cobj,
):
    pcldfmc = Ptr[float](pcldfmc_p)
    ztauc = Ptr[float](ztauc_p)
    ztauo = Ptr[float](ztauo_p)
    zrefc = Ptr[float](zrefc_p)
    zrefdc = Ptr[float](zrefdc_p)
    ztrac = Ptr[float](ztrac_p)
    ztradc = Ptr[float](ztradc_p)
    zrefo = Ptr[float](zrefo_p)
    zrefdo = Ptr[float](zrefdo_p)
    ztrao = Ptr[float](ztrao_p)
    ztrado = Ptr[float](ztrado_p)
    zref = Ptr[float](zref_p)
    zrefd = Ptr[float](zrefd_p)
    ztra = Ptr[float](ztra_p)
    ztrad = Ptr[float](ztrad_p)
    zdbtc = Ptr[float](zdbtc_p)
    ztdbtc = Ptr[float](ztdbtc_p)
    zdbt = Ptr[float](zdbt_p)
    ztdbt = Ptr[float](ztdbt_p)
    exp_tbl = Ptr[float](exp_tbl_p)

    for jk in range(1, klev + 1):
        ikl = klev + 1 - jk
        jidx = jk - 1
        mc_idx = _idx2(ikl, iw, klev)

        zclear = 1.0 - pcldfmc[mc_idx]
        zcloud = pcldfmc[mc_idx]

        zref[jidx] = zclear * zrefc[jidx] + zcloud * zrefo[jidx]
        zrefd[jidx] = zclear * zrefdc[jidx] + zcloud * zrefdo[jidx]
        ztra[jidx] = zclear * ztrac[jidx] + zcloud * ztrao[jidx]
        ztrad[jidx] = zclear * ztradc[jidx] + zcloud * ztrado[jidx]

        ze1 = ztauc[jidx] / prmu0
        zdbtmc = _rrtmg_sw_exp_lookup(ze1, tblint, bpade, od_lo, exp_tbl)
        zdbtc[jidx] = zdbtmc
        ztdbtc[jk] = zdbtc[jidx] * ztdbtc[jidx]

        ze1 = ztauo[jidx] / prmu0
        zdbtmo = _rrtmg_sw_exp_lookup(ze1, tblint, bpade, od_lo, exp_tbl)
        zdbt[jidx] = zclear * zdbtmc + zcloud * zdbtmo
        ztdbt[jk] = zdbt[jidx] * ztdbt[jidx]


@export
def rrtmg_sw_spcvmc_flux_codon(
    klev: int,
    ngptsw: int,
    nbndsw: int,
    iw: int,
    ibm: int,
    idelm: int,
    zincflx_p: cobj,
    zfu_p: cobj,
    zfd_p: cobj,
    zcu_p: cobj,
    zcd_p: cobj,
    ztdbt_nodel_p: cobj,
    ztdbtc_nodel_p: cobj,
    ztdbt_p: cobj,
    ztdbtc_p: cobj,
    pbbfsu_p: cobj,
    pbbfsd_p: cobj,
    pbbfu_p: cobj,
    pbbfd_p: cobj,
    pbbcu_p: cobj,
    pbbcd_p: cobj,
    pbbfddir_p: cobj,
    pbbcddir_p: cobj,
    puvcd_p: cobj,
    puvfd_p: cobj,
    puvcddir_p: cobj,
    puvfddir_p: cobj,
    pnicd_p: cobj,
    pnifd_p: cobj,
    pnicddir_p: cobj,
    pnifddir_p: cobj,
    pnicu_p: cobj,
    pnifu_p: cobj,
):
    zincflx = Ptr[float](zincflx_p)
    zfu = Ptr[float](zfu_p)
    zfd = Ptr[float](zfd_p)
    zcu = Ptr[float](zcu_p)
    zcd = Ptr[float](zcd_p)
    ztdbt_nodel = Ptr[float](ztdbt_nodel_p)
    ztdbtc_nodel = Ptr[float](ztdbtc_nodel_p)
    ztdbt = Ptr[float](ztdbt_p)
    ztdbtc = Ptr[float](ztdbtc_p)
    pbbfsu = Ptr[float](pbbfsu_p)
    pbbfsd = Ptr[float](pbbfsd_p)
    pbbfu = Ptr[float](pbbfu_p)
    pbbfd = Ptr[float](pbbfd_p)
    pbbcu = Ptr[float](pbbcu_p)
    pbbcd = Ptr[float](pbbcd_p)
    pbbfddir = Ptr[float](pbbfddir_p)
    pbbcddir = Ptr[float](pbbcddir_p)
    puvcd = Ptr[float](puvcd_p)
    puvfd = Ptr[float](puvfd_p)
    puvcddir = Ptr[float](puvcddir_p)
    puvfddir = Ptr[float](puvfddir_p)
    pnicd = Ptr[float](pnicd_p)
    pnifd = Ptr[float](pnifd_p)
    pnicddir = Ptr[float](pnicddir_p)
    pnifddir = Ptr[float](pnifddir_p)
    pnicu = Ptr[float](pnicu_p)
    pnifu = Ptr[float](pnifu_p)

    ldlev = klev + 1
    for jk in range(1, klev + 2):
        ikl = klev + 2 - jk
        ikl_idx = ikl - 1

        pbbfsu[_idx2(ibm, ikl, nbndsw)] = (
            pbbfsu[_idx2(ibm, ikl, nbndsw)]
            + zincflx[iw - 1] * zfu[_idx2(jk, iw, ldlev)]
        )
        pbbfsd[_idx2(ibm, ikl, nbndsw)] = (
            pbbfsd[_idx2(ibm, ikl, nbndsw)]
            + zincflx[iw - 1] * zfd[_idx2(jk, iw, ldlev)]
        )

        pbbfu[ikl_idx] = pbbfu[ikl_idx] + zincflx[iw - 1] * zfu[_idx2(jk, iw, ldlev)]
        pbbfd[ikl_idx] = pbbfd[ikl_idx] + zincflx[iw - 1] * zfd[_idx2(jk, iw, ldlev)]
        pbbcu[ikl_idx] = pbbcu[ikl_idx] + zincflx[iw - 1] * zcu[_idx2(jk, iw, ldlev)]
        pbbcd[ikl_idx] = pbbcd[ikl_idx] + zincflx[iw - 1] * zcd[_idx2(jk, iw, ldlev)]
        if idelm == 0:
            pbbfddir[ikl_idx] = pbbfddir[ikl_idx] + zincflx[iw - 1] * ztdbt_nodel[jk - 1]
            pbbcddir[ikl_idx] = pbbcddir[ikl_idx] + zincflx[iw - 1] * ztdbtc_nodel[jk - 1]
        elif idelm == 1:
            pbbfddir[ikl_idx] = pbbfddir[ikl_idx] + zincflx[iw - 1] * ztdbt[jk - 1]
            pbbcddir[ikl_idx] = pbbcddir[ikl_idx] + zincflx[iw - 1] * ztdbtc[jk - 1]

        if ibm >= 10 and ibm <= 13:
            puvcd[ikl_idx] = puvcd[ikl_idx] + zincflx[iw - 1] * zcd[_idx2(jk, iw, ldlev)]
            puvfd[ikl_idx] = puvfd[ikl_idx] + zincflx[iw - 1] * zfd[_idx2(jk, iw, ldlev)]
            if idelm == 0:
                puvfddir[ikl_idx] = (
                    puvfddir[ikl_idx] + zincflx[iw - 1] * ztdbt_nodel[jk - 1]
                )
                puvcddir[ikl_idx] = (
                    puvcddir[ikl_idx] + zincflx[iw - 1] * ztdbtc_nodel[jk - 1]
                )
            elif idelm == 1:
                puvfddir[ikl_idx] = (
                    puvfddir[ikl_idx] + zincflx[iw - 1] * ztdbt[jk - 1]
                )
                puvcddir[ikl_idx] = (
                    puvcddir[ikl_idx] + zincflx[iw - 1] * ztdbtc[jk - 1]
                )
        elif ibm == 9:
            puvcd[ikl_idx] = (
                puvcd[ikl_idx] + 0.5 * zincflx[iw - 1] * zcd[_idx2(jk, iw, ldlev)]
            )
            puvfd[ikl_idx] = (
                puvfd[ikl_idx] + 0.5 * zincflx[iw - 1] * zfd[_idx2(jk, iw, ldlev)]
            )
            pnicd[ikl_idx] = (
                pnicd[ikl_idx] + 0.5 * zincflx[iw - 1] * zcd[_idx2(jk, iw, ldlev)]
            )
            pnifd[ikl_idx] = (
                pnifd[ikl_idx] + 0.5 * zincflx[iw - 1] * zfd[_idx2(jk, iw, ldlev)]
            )
            if idelm == 0:
                puvfddir[ikl_idx] = (
                    puvfddir[ikl_idx] + 0.5 * zincflx[iw - 1] * ztdbt_nodel[jk - 1]
                )
                puvcddir[ikl_idx] = (
                    puvcddir[ikl_idx] + 0.5 * zincflx[iw - 1] * ztdbtc_nodel[jk - 1]
                )
                pnifddir[ikl_idx] = (
                    pnifddir[ikl_idx] + 0.5 * zincflx[iw - 1] * ztdbt_nodel[jk - 1]
                )
                pnicddir[ikl_idx] = (
                    pnicddir[ikl_idx] + 0.5 * zincflx[iw - 1] * ztdbtc_nodel[jk - 1]
                )
            elif idelm == 1:
                puvfddir[ikl_idx] = (
                    puvfddir[ikl_idx] + 0.5 * zincflx[iw - 1] * ztdbt[jk - 1]
                )
                puvcddir[ikl_idx] = (
                    puvcddir[ikl_idx] + 0.5 * zincflx[iw - 1] * ztdbtc[jk - 1]
                )
                pnifddir[ikl_idx] = (
                    pnifddir[ikl_idx] + 0.5 * zincflx[iw - 1] * ztdbt[jk - 1]
                )
                pnicddir[ikl_idx] = (
                    pnicddir[ikl_idx] + 0.5 * zincflx[iw - 1] * ztdbtc[jk - 1]
                )
            pnicu[ikl_idx] = (
                pnicu[ikl_idx] + 0.5 * zincflx[iw - 1] * zcu[_idx2(jk, iw, ldlev)]
            )
            pnifu[ikl_idx] = (
                pnifu[ikl_idx] + 0.5 * zincflx[iw - 1] * zfu[_idx2(jk, iw, ldlev)]
            )
        elif ibm == 14 or ibm <= 8:
            pnicd[ikl_idx] = pnicd[ikl_idx] + zincflx[iw - 1] * zcd[_idx2(jk, iw, ldlev)]
            pnifd[ikl_idx] = pnifd[ikl_idx] + zincflx[iw - 1] * zfd[_idx2(jk, iw, ldlev)]
            if idelm == 0:
                pnifddir[ikl_idx] = (
                    pnifddir[ikl_idx] + zincflx[iw - 1] * ztdbt_nodel[jk - 1]
                )
                pnicddir[ikl_idx] = (
                    pnicddir[ikl_idx] + zincflx[iw - 1] * ztdbtc_nodel[jk - 1]
                )
            elif idelm == 1:
                pnifddir[ikl_idx] = (
                    pnifddir[ikl_idx] + zincflx[iw - 1] * ztdbt[jk - 1]
                )
                pnicddir[ikl_idx] = (
                    pnicddir[ikl_idx] + zincflx[iw - 1] * ztdbtc[jk - 1]
                )
            pnicu[ikl_idx] = pnicu[ikl_idx] + zincflx[iw - 1] * zcu[_idx2(jk, iw, ldlev)]
            pnifu[ikl_idx] = pnifu[ikl_idx] + zincflx[iw - 1] * zfu[_idx2(jk, iw, ldlev)]


@export
def rrtmg_sw_spcvmc_zero_outputs_codon(
    klev: int,
    pbbcd_p: cobj,
    pbbcu_p: cobj,
    pbbfd_p: cobj,
    pbbfu_p: cobj,
    pbbcddir_p: cobj,
    pbbfddir_p: cobj,
    puvcd_p: cobj,
    puvfd_p: cobj,
    puvcddir_p: cobj,
    puvfddir_p: cobj,
    pnicd_p: cobj,
    pnifd_p: cobj,
    pnicddir_p: cobj,
    pnifddir_p: cobj,
    pnicu_p: cobj,
    pnifu_p: cobj,
):
    pbbcd = Ptr[float](pbbcd_p)
    pbbcu = Ptr[float](pbbcu_p)
    pbbfd = Ptr[float](pbbfd_p)
    pbbfu = Ptr[float](pbbfu_p)
    pbbcddir = Ptr[float](pbbcddir_p)
    pbbfddir = Ptr[float](pbbfddir_p)
    puvcd = Ptr[float](puvcd_p)
    puvfd = Ptr[float](puvfd_p)
    puvcddir = Ptr[float](puvcddir_p)
    puvfddir = Ptr[float](puvfddir_p)
    pnicd = Ptr[float](pnicd_p)
    pnifd = Ptr[float](pnifd_p)
    pnicddir = Ptr[float](pnicddir_p)
    pnifddir = Ptr[float](pnifddir_p)
    pnicu = Ptr[float](pnicu_p)
    pnifu = Ptr[float](pnifu_p)

    for jk in range(1, klev + 2):
        idx = jk - 1
        pbbcd[idx] = 0.0
        pbbcu[idx] = 0.0
        pbbfd[idx] = 0.0
        pbbfu[idx] = 0.0
        pbbcddir[idx] = 0.0
        pbbfddir[idx] = 0.0
        puvcd[idx] = 0.0
        puvfd[idx] = 0.0
        puvcddir[idx] = 0.0
        puvfddir[idx] = 0.0
        pnicd[idx] = 0.0
        pnifd[idx] = 0.0
        pnicddir[idx] = 0.0
        pnifddir[idx] = 0.0
        pnicu[idx] = 0.0
        pnifu[idx] = 0.0


@inline
def _rrtmg_lw_a0(iband: int) -> float:
    if iband == 2:
        return 1.55
    if iband == 3:
        return 1.58
    if iband == 5:
        return 1.54
    if iband == 6:
        return 1.454
    if iband == 7:
        return 1.89
    if iband == 8:
        return 1.33
    if iband == 9:
        return 1.668
    return 1.66


@inline
def _rrtmg_lw_a1(iband: int) -> float:
    if iband == 2:
        return 0.25
    if iband == 3:
        return 0.22
    if iband == 5:
        return 0.13
    if iband == 6:
        return 0.446
    if iband == 7:
        return -0.10
    if iband == 8:
        return 0.40
    if iband == 9:
        return -0.006
    return 0.00


@inline
def _rrtmg_lw_a2(iband: int) -> float:
    if iband == 2:
        return -12.0
    if iband == 3:
        return -11.7
    if iband == 5:
        return -0.72
    if iband == 6:
        return -0.243
    if iband == 7:
        return 0.19
    if iband == 8:
        return -0.062
    if iband == 9:
        return 0.414
    return 0.00


@export
def rrtmg_lw_rtrnmc_codon(
    nlayers: int,
    istart: int,
    iend: int,
    iout: int,
    ncbands: int,
    nbndlw: int,
    ngptlw: int,
    pwvcm: float,
    fluxfac: float,
    heatfac: float,
    tblint: float,
    bpade: float,
    pz_p: cobj,
    semiss_p: cobj,
    cldfmc_p: cobj,
    taucmc_p: cobj,
    planklay_p: cobj,
    planklev_p: cobj,
    plankbnd_p: cobj,
    fracs_p: cobj,
    taut_p: cobj,
    totuflux_p: cobj,
    totdflux_p: cobj,
    fnet_p: cobj,
    htr_p: cobj,
    totuclfl_p: cobj,
    totdclfl_p: cobj,
    fnetc_p: cobj,
    htrc_p: cobj,
    totufluxs_p: cobj,
    totdfluxs_p: cobj,
    tau_tbl_p: cobj,
    exp_tbl_p: cobj,
    tfn_tbl_p: cobj,
    delwave_p: cobj,
    ngs_p: cobj,
    ngb_p: cobj,
    abscld_p: cobj,
    atot_p: cobj,
    atrans_p: cobj,
    bbugas_p: cobj,
    bbutot_p: cobj,
    clrurad_p: cobj,
    clrdrad_p: cobj,
    efclfrac_p: cobj,
    uflux_p: cobj,
    dflux_p: cobj,
    urad_p: cobj,
    drad_p: cobj,
    uclfl_p: cobj,
    dclfl_p: cobj,
    odcld_p: cobj,
    secdiff_p: cobj,
    icldlyr_p: cobj,
):
    pz = Ptr[float](pz_p)
    semiss = Ptr[float](semiss_p)
    cldfmc = Ptr[float](cldfmc_p)
    taucmc = Ptr[float](taucmc_p)
    planklay = Ptr[float](planklay_p)
    planklev = Ptr[float](planklev_p)
    plankbnd = Ptr[float](plankbnd_p)
    fracs = Ptr[float](fracs_p)
    taut = Ptr[float](taut_p)
    totuflux = Ptr[float](totuflux_p)
    totdflux = Ptr[float](totdflux_p)
    fnet = Ptr[float](fnet_p)
    htr = Ptr[float](htr_p)
    totuclfl = Ptr[float](totuclfl_p)
    totdclfl = Ptr[float](totdclfl_p)
    fnetc = Ptr[float](fnetc_p)
    htrc = Ptr[float](htrc_p)
    totufluxs = Ptr[float](totufluxs_p)
    totdfluxs = Ptr[float](totdfluxs_p)
    tau_tbl = Ptr[float](tau_tbl_p)
    exp_tbl = Ptr[float](exp_tbl_p)
    tfn_tbl = Ptr[float](tfn_tbl_p)
    delwave = Ptr[float](delwave_p)
    ngs = Ptr[int](ngs_p)
    ngb = Ptr[int](ngb_p)
    abscld = Ptr[float](abscld_p)
    atot = Ptr[float](atot_p)
    atrans = Ptr[float](atrans_p)
    bbugas = Ptr[float](bbugas_p)
    bbutot = Ptr[float](bbutot_p)
    clrurad = Ptr[float](clrurad_p)
    clrdrad = Ptr[float](clrdrad_p)
    efclfrac = Ptr[float](efclfrac_p)
    uflux = Ptr[float](uflux_p)
    dflux = Ptr[float](dflux_p)
    urad = Ptr[float](urad_p)
    drad = Ptr[float](drad_p)
    uclfl = Ptr[float](uclfl_p)
    dclfl = Ptr[float](dclfl_p)
    odcld = Ptr[float](odcld_p)
    secdiff = Ptr[float](secdiff_p)
    icldlyr = Ptr[int](icldlyr_p)

    wtdiff = 0.5
    rec_6 = 0.166667

    for ibnd in range(1, nbndlw + 1):
        if ibnd == 1 or ibnd == 4 or ibnd >= 10:
            secdiff[ibnd - 1] = 1.66
        else:
            secdiff[ibnd - 1] = _rrtmg_lw_a0(ibnd) + _rrtmg_lw_a1(ibnd) * exp(
                _rrtmg_lw_a2(ibnd) * pwvcm
            )
            if secdiff[ibnd - 1] > 1.80:
                secdiff[ibnd - 1] = 1.80
            if secdiff[ibnd - 1] < 1.50:
                secdiff[ibnd - 1] = 1.50

    urad[0] = 0.0
    drad[0] = 0.0
    totuflux[0] = 0.0
    totdflux[0] = 0.0
    clrurad[0] = 0.0
    clrdrad[0] = 0.0
    totuclfl[0] = 0.0
    totdclfl[0] = 0.0

    for lay in range(1, nlayers + 1):
        urad[lay] = 0.0
        drad[lay] = 0.0
        totuflux[lay] = 0.0
        totdflux[lay] = 0.0
        clrurad[lay] = 0.0
        clrdrad[lay] = 0.0
        totuclfl[lay] = 0.0
        totdclfl[lay] = 0.0
        icldlyr[lay - 1] = 0

        for ig in range(1, ngptlw + 1):
            rt_idx = _idx2(lay, ig, nlayers)
            mc_idx = _idx2(ig, lay, ngptlw)
            if cldfmc[mc_idx] == 1.0:
                ib = ngb[ig - 1]
                odcld[rt_idx] = secdiff[ib - 1] * taucmc[mc_idx]
                transcld = exp(-odcld[rt_idx])
                abscld[rt_idx] = 1.0 - transcld
                efclfrac[rt_idx] = abscld[rt_idx] * cldfmc[mc_idx]
                icldlyr[lay - 1] = 1
            else:
                odcld[rt_idx] = 0.0
                abscld[rt_idx] = 0.0
                efclfrac[rt_idx] = 0.0

    igc = 1
    for iband in range(istart, iend + 1):
        if iout > 0 and iband >= 2:
            igc = ngs[iband - 2] + 1

        while True:
            radld = 0.0
            radclrd = 0.0
            iclddn = 0

            for lev in range(nlayers, 0, -1):
                plfrac = fracs[_idx2(lev, igc, nlayers)]
                blay = planklay[_idx2(lev, iband, nlayers)]
                dplankup = planklev[_idx2_dim1_lb0(lev, iband, nlayers)] - blay
                dplankdn = planklev[_idx2_dim1_lb0(lev - 1, iband, nlayers)] - blay
                odepth = secdiff[iband - 1] * taut[_idx2(lev, igc, nlayers)]
                if odepth < 0.0:
                    odepth = 0.0

                if icldlyr[lev - 1] == 1:
                    iclddn = 1
                    odtot = odepth + odcld[_idx2(lev, igc, nlayers)]
                    if odtot < 0.06:
                        atrans[lev - 1] = odepth - 0.5 * odepth * odepth
                        odepth_rec = rec_6 * odepth
                        gassrc = plfrac * (blay + dplankdn * odepth_rec) * atrans[lev - 1]

                        atot[lev - 1] = odtot - 0.5 * odtot * odtot
                        odtot_rec = rec_6 * odtot
                        bbdtot = plfrac * (blay + dplankdn * odtot_rec)
                        bbd = plfrac * (blay + dplankdn * odepth_rec)
                        radld = (
                            radld
                            - radld
                            * (
                                atrans[lev - 1]
                                + efclfrac[_idx2(lev, igc, nlayers)]
                                * (1.0 - atrans[lev - 1])
                            )
                            + gassrc
                            + cldfmc[_idx2(igc, lev, ngptlw)]
                            * (bbdtot * atot[lev - 1] - gassrc)
                        )
                        drad[lev - 1] = drad[lev - 1] + radld

                        bbugas[lev - 1] = plfrac * (blay + dplankup * odepth_rec)
                        bbutot[lev - 1] = plfrac * (blay + dplankup * odtot_rec)
                    elif odepth <= 0.06:
                        atrans[lev - 1] = odepth - 0.5 * odepth * odepth
                        odepth_rec = rec_6 * odepth
                        gassrc = plfrac * (blay + dplankdn * odepth_rec) * atrans[lev - 1]

                        odtot = odepth + odcld[_idx2(lev, igc, nlayers)]
                        tblind = odtot / (bpade + odtot)
                        ittot = int(tblint * tblind + 0.5)
                        tfactot = tfn_tbl[ittot]
                        bbdtot = plfrac * (blay + tfactot * dplankdn)
                        bbd = plfrac * (blay + dplankdn * odepth_rec)
                        atot[lev - 1] = 1.0 - exp_tbl[ittot]

                        radld = (
                            radld
                            - radld
                            * (
                                atrans[lev - 1]
                                + efclfrac[_idx2(lev, igc, nlayers)]
                                * (1.0 - atrans[lev - 1])
                            )
                            + gassrc
                            + cldfmc[_idx2(igc, lev, ngptlw)]
                            * (bbdtot * atot[lev - 1] - gassrc)
                        )
                        drad[lev - 1] = drad[lev - 1] + radld

                        bbugas[lev - 1] = plfrac * (blay + dplankup * odepth_rec)
                        bbutot[lev - 1] = plfrac * (blay + tfactot * dplankup)
                    else:
                        tblind = odepth / (bpade + odepth)
                        itgas = int(tblint * tblind + 0.5)
                        odepth = tau_tbl[itgas]
                        atrans[lev - 1] = 1.0 - exp_tbl[itgas]
                        tfacgas = tfn_tbl[itgas]
                        gassrc = (
                            atrans[lev - 1] * plfrac * (blay + tfacgas * dplankdn)
                        )

                        odtot = odepth + odcld[_idx2(lev, igc, nlayers)]
                        tblind = odtot / (bpade + odtot)
                        ittot = int(tblint * tblind + 0.5)
                        tfactot = tfn_tbl[ittot]
                        bbdtot = plfrac * (blay + tfactot * dplankdn)
                        bbd = plfrac * (blay + tfacgas * dplankdn)
                        atot[lev - 1] = 1.0 - exp_tbl[ittot]

                        radld = (
                            radld
                            - radld
                            * (
                                atrans[lev - 1]
                                + efclfrac[_idx2(lev, igc, nlayers)]
                                * (1.0 - atrans[lev - 1])
                            )
                            + gassrc
                            + cldfmc[_idx2(igc, lev, ngptlw)]
                            * (bbdtot * atot[lev - 1] - gassrc)
                        )
                        drad[lev - 1] = drad[lev - 1] + radld
                        bbugas[lev - 1] = plfrac * (blay + tfacgas * dplankup)
                        bbutot[lev - 1] = plfrac * (blay + tfactot * dplankup)
                else:
                    if odepth <= 0.06:
                        atrans[lev - 1] = odepth - 0.5 * odepth * odepth
                        odepth = rec_6 * odepth
                        bbd = plfrac * (blay + dplankdn * odepth)
                        bbugas[lev - 1] = plfrac * (blay + dplankup * odepth)
                    else:
                        tblind = odepth / (bpade + odepth)
                        itr = int(tblint * tblind + 0.5)
                        transc = exp_tbl[itr]
                        atrans[lev - 1] = 1.0 - transc
                        tausfac = tfn_tbl[itr]
                        bbd = plfrac * (blay + tausfac * dplankdn)
                        bbugas[lev - 1] = plfrac * (blay + tausfac * dplankup)
                    radld = radld + (bbd - radld) * atrans[lev - 1]
                    drad[lev - 1] = drad[lev - 1] + radld

                if iclddn == 1:
                    radclrd = radclrd + (bbd - radclrd) * atrans[lev - 1]
                    clrdrad[lev - 1] = clrdrad[lev - 1] + radclrd
                else:
                    radclrd = radld
                    clrdrad[lev - 1] = drad[lev - 1]

            rad0 = fracs[_idx2(1, igc, nlayers)] * plankbnd[iband - 1]
            reflect = 1.0 - semiss[iband - 1]
            radlu = rad0 + reflect * radld
            radclru = rad0 + reflect * radclrd

            urad[0] = urad[0] + radlu
            clrurad[0] = clrurad[0] + radclru

            for lev in range(1, nlayers + 1):
                if icldlyr[lev - 1] == 1:
                    gassrc = bbugas[lev - 1] * atrans[lev - 1]
                    radlu = (
                        radlu
                        - radlu
                        * (
                            atrans[lev - 1]
                            + efclfrac[_idx2(lev, igc, nlayers)]
                            * (1.0 - atrans[lev - 1])
                        )
                        + gassrc
                        + cldfmc[_idx2(igc, lev, ngptlw)]
                        * (bbutot[lev - 1] * atot[lev - 1] - gassrc)
                    )
                    urad[lev] = urad[lev] + radlu
                else:
                    radlu = radlu + (bbugas[lev - 1] - radlu) * atrans[lev - 1]
                    urad[lev] = urad[lev] + radlu

                if iclddn == 1:
                    radclru = radclru + (bbugas[lev - 1] - radclru) * atrans[lev - 1]
                    clrurad[lev] = clrurad[lev] + radclru
                else:
                    radclru = radlu
                    clrurad[lev] = urad[lev]

            igc = igc + 1
            if igc > ngs[iband - 1]:
                break

        for lev in range(nlayers, -1, -1):
            uflux[lev] = urad[lev] * wtdiff
            dflux[lev] = drad[lev] * wtdiff
            urad[lev] = 0.0
            drad[lev] = 0.0
            totuflux[lev] = totuflux[lev] + uflux[lev] * delwave[iband - 1]
            totdflux[lev] = totdflux[lev] + dflux[lev] * delwave[iband - 1]
            uclfl[lev] = clrurad[lev] * wtdiff
            dclfl[lev] = clrdrad[lev] * wtdiff
            clrurad[lev] = 0.0
            clrdrad[lev] = 0.0
            totuclfl[lev] = totuclfl[lev] + uclfl[lev] * delwave[iband - 1]
            totdclfl[lev] = totdclfl[lev] + dclfl[lev] * delwave[iband - 1]
            totufluxs[_idx2_lb0(iband, lev, nbndlw)] = (
                uflux[lev] * delwave[iband - 1]
            )
            totdfluxs[_idx2_lb0(iband, lev, nbndlw)] = (
                dflux[lev] * delwave[iband - 1]
            )

    totuflux[0] = totuflux[0] * fluxfac
    totdflux[0] = totdflux[0] * fluxfac
    for iband in range(1, nbndlw + 1):
        totufluxs[_idx2_lb0(iband, 0, nbndlw)] = (
            totufluxs[_idx2_lb0(iband, 0, nbndlw)] * fluxfac
        )
        totdfluxs[_idx2_lb0(iband, 0, nbndlw)] = (
            totdfluxs[_idx2_lb0(iband, 0, nbndlw)] * fluxfac
        )
    fnet[0] = totuflux[0] - totdflux[0]
    totuclfl[0] = totuclfl[0] * fluxfac
    totdclfl[0] = totdclfl[0] * fluxfac
    fnetc[0] = totuclfl[0] - totdclfl[0]

    for lev in range(1, nlayers + 1):
        totuflux[lev] = totuflux[lev] * fluxfac
        totdflux[lev] = totdflux[lev] * fluxfac
        for iband in range(1, nbndlw + 1):
            totufluxs[_idx2_lb0(iband, lev, nbndlw)] = (
                totufluxs[_idx2_lb0(iband, lev, nbndlw)] * fluxfac
            )
            totdfluxs[_idx2_lb0(iband, lev, nbndlw)] = (
                totdfluxs[_idx2_lb0(iband, lev, nbndlw)] * fluxfac
            )
        fnet[lev] = totuflux[lev] - totdflux[lev]
        totuclfl[lev] = totuclfl[lev] * fluxfac
        totdclfl[lev] = totdclfl[lev] * fluxfac
        fnetc[lev] = totuclfl[lev] - totdclfl[lev]
        l = lev - 1
        htr[l] = heatfac * (fnet[l] - fnet[lev]) / (pz[l] - pz[lev])
        htrc[l] = heatfac * (fnetc[l] - fnetc[lev]) / (pz[l] - pz[lev])

    htr[nlayers] = 0.0
    htrc[nlayers] = 0.0


@inline
def _rrtmg_sw_reftra_lookup(
    ze: float, tblint: float, bpade: float, od_lo: float, exp_tbl: Ptr[float]
) -> float:
    if ze <= od_lo:
        return 1.0 - ze + 0.5 * ze * ze
    tblind = ze / (bpade + ze)
    itind = int(tblint * tblind + 0.5)
    return exp_tbl[itind]


@export
def rrtmg_sw_reftra_codon(
    nlayers: int,
    prmuz: float,
    tblint: float,
    bpade: float,
    od_lo: float,
    lrtchk_mask_p: cobj,
    pgg_p: cobj,
    ptau_p: cobj,
    pw_p: cobj,
    pref_p: cobj,
    prefd_p: cobj,
    ptra_p: cobj,
    ptrad_p: cobj,
    exp_tbl_p: cobj,
):
    lrtchk_mask = Ptr[int](lrtchk_mask_p)
    pgg = Ptr[float](pgg_p)
    ptau = Ptr[float](ptau_p)
    pw = Ptr[float](pw_p)
    pref = Ptr[float](pref_p)
    prefd = Ptr[float](prefd_p)
    ptra = Ptr[float](ptra_p)
    ptrad = Ptr[float](ptrad_p)
    exp_tbl = Ptr[float](exp_tbl_p)

    zsr3 = sqrt(3.0)
    zwcrit = 0.9999995
    kmodts = 2
    eps = 1.0e-08

    for jk in range(1, nlayers + 1):
        idx = jk - 1
        if lrtchk_mask[idx] == 0:
            pref[idx] = 0.0
            ptra[idx] = 1.0
            prefd[idx] = 0.0
            ptrad[idx] = 1.0
        else:
            zto1 = ptau[idx]
            zw = pw[idx]
            zg = pgg[idx]

            zg3 = 3.0 * zg
            if kmodts == 1:
                zgamma1 = (7.0 - zw * (4.0 + zg3)) * 0.25
                zgamma2 = -(1.0 - zw * (4.0 - zg3)) * 0.25
                zgamma3 = (2.0 - zg3 * prmuz) * 0.25
            elif kmodts == 2:
                zgamma1 = (8.0 - zw * (5.0 + zg3)) * 0.25
                zgamma2 = 3.0 * (zw * (1.0 - zg)) * 0.25
                zgamma3 = (2.0 - zg3 * prmuz) * 0.25
            elif kmodts == 3:
                zgamma1 = zsr3 * (2.0 - zw * (1.0 + zg)) * 0.5
                zgamma2 = zsr3 * zw * (1.0 - zg) * 0.5
                zgamma3 = (1.0 - zsr3 * zg * prmuz) * 0.5
            else:
                zgamma1 = 0.0
                zgamma2 = 0.0
                zgamma3 = 0.0
            zgamma4 = 1.0 - zgamma3

            zwo = zw / (1.0 - (1.0 - zw) * (zg / (1.0 - zg)) ** 2)

            if zwo >= zwcrit:
                za = zgamma1 * prmuz
                za1 = za - zgamma3
                zgt = zgamma1 * zto1

                ze1 = min(zto1 / prmuz, 500.0)
                ze2 = _rrtmg_sw_reftra_lookup(ze1, tblint, bpade, od_lo, exp_tbl)

                pref[idx] = (zgt - za1 * (1.0 - ze2)) / (1.0 + zgt)
                ptra[idx] = 1.0 - pref[idx]

                prefd[idx] = zgt / (1.0 + zgt)
                ptrad[idx] = 1.0 - prefd[idx]

                if ze2 == 1.0:
                    pref[idx] = 0.0
                    ptra[idx] = 1.0
                    prefd[idx] = 0.0
                    ptrad[idx] = 1.0
            else:
                za1 = zgamma1 * zgamma4 + zgamma2 * zgamma3
                za2 = zgamma1 * zgamma3 + zgamma2 * zgamma4
                zrk = sqrt(zgamma1**2 - zgamma2**2)
                zrp = zrk * prmuz
                zrp1 = 1.0 + zrp
                zrm1 = 1.0 - zrp
                zrk2 = 2.0 * zrk
                zrpp = 1.0 - zrp * zrp
                zrkg = zrk + zgamma1
                zr1 = zrm1 * (za2 + zrk * zgamma3)
                zr2 = zrp1 * (za2 - zrk * zgamma3)
                zr3 = zrk2 * (zgamma3 - za2 * prmuz)
                zr4 = zrpp * zrkg
                zr5 = zrpp * (zrk - zgamma1)
                zt1 = zrp1 * (za1 + zrk * zgamma4)
                zt2 = zrm1 * (za1 - zrk * zgamma4)
                zt3 = zrk2 * (zgamma4 + za1 * prmuz)
                zt4 = zr4
                zt5 = zr5
                zbeta = (zgamma1 - zrk) / zrkg

                ze1 = min(zrk * zto1, 500.0)
                ze2 = min(zto1 / prmuz, 500.0)

                zem1 = _rrtmg_sw_reftra_lookup(ze1, tblint, bpade, od_lo, exp_tbl)
                zep1 = 1.0 / zem1

                zem2 = _rrtmg_sw_reftra_lookup(ze2, tblint, bpade, od_lo, exp_tbl)
                zep2 = 1.0 / zem2

                zdenr = zr4 * zep1 + zr5 * zem1
                zdent = zt4 * zep1 + zt5 * zem1
                if zdenr >= -eps and zdenr <= eps:
                    pref[idx] = eps
                    ptra[idx] = zem2
                else:
                    pref[idx] = zw * (zr1 * zep1 - zr2 * zem1 - zr3 * zem2) / zdenr
                    ptra[idx] = zem2 - zem2 * zw * (
                        zt1 * zep1 - zt2 * zem1 - zt3 * zep2
                    ) / zdent

                zemm = zem1 * zem1
                zdend = 1.0 / ((1.0 - zbeta * zemm) * zrkg)
                prefd[idx] = zgamma2 * (1.0 - zemm) * zdend
                ptrad[idx] = zrk2 * zem1 * zdend


@export
def rrtmg_sw_vrtqdr_codon(
    klev: int,
    kw: int,
    pref_p: cobj,
    prefd_p: cobj,
    ptra_p: cobj,
    ptrad_p: cobj,
    pdbt_p: cobj,
    prdnd_p: cobj,
    prup_p: cobj,
    prupd_p: cobj,
    ptdbt_p: cobj,
    pfd_p: cobj,
    pfu_p: cobj,
    ztdn_p: cobj,
):
    pref = Ptr[float](pref_p)
    prefd = Ptr[float](prefd_p)
    ptra = Ptr[float](ptra_p)
    ptrad = Ptr[float](ptrad_p)
    pdbt = Ptr[float](pdbt_p)
    prdnd = Ptr[float](prdnd_p)
    prup = Ptr[float](prup_p)
    prupd = Ptr[float](prupd_p)
    ptdbt = Ptr[float](ptdbt_p)
    pfd = Ptr[float](pfd_p)
    pfu = Ptr[float](pfu_p)
    ztdn = Ptr[float](ztdn_p)

    zreflect = 1.0 / (1.0 - prefd[klev] * prefd[klev - 1])
    prup[klev - 1] = pref[klev - 1] + (
        ptrad[klev - 1]
        * (
            (ptra[klev - 1] - pdbt[klev - 1]) * prefd[klev]
            + pdbt[klev - 1] * pref[klev]
        )
    ) * zreflect
    prupd[klev - 1] = (
        prefd[klev - 1]
        + ptrad[klev - 1] * ptrad[klev - 1] * prefd[klev] * zreflect
    )

    for jk in range(1, klev):
        ikp = klev + 1 - jk
        ikx = ikp - 1
        zreflect = 1.0 / (1.0 - prupd[ikp - 1] * prefd[ikx - 1])
        prup[ikx - 1] = pref[ikx - 1] + (
            ptrad[ikx - 1]
            * (
                (ptra[ikx - 1] - pdbt[ikx - 1]) * prupd[ikp - 1]
                + pdbt[ikx - 1] * prup[ikp - 1]
            )
        ) * zreflect
        prupd[ikx - 1] = (
            prefd[ikx - 1]
            + ptrad[ikx - 1] * ptrad[ikx - 1] * prupd[ikp - 1] * zreflect
        )

    ztdn[0] = 1.0
    prdnd[0] = 0.0
    ztdn[1] = ptra[0]
    prdnd[1] = prefd[0]

    for jk in range(2, klev + 1):
        ikp = jk + 1
        zreflect = 1.0 / (1.0 - prefd[jk - 1] * prdnd[jk - 1])
        ztdn[ikp - 1] = ptdbt[jk - 1] * ptra[jk - 1] + (
            ptrad[jk - 1]
            * (
                (ztdn[jk - 1] - ptdbt[jk - 1])
                + ptdbt[jk - 1] * pref[jk - 1] * prdnd[jk - 1]
            )
        ) * zreflect
        prdnd[ikp - 1] = (
            prefd[jk - 1]
            + ptrad[jk - 1] * ptrad[jk - 1] * prdnd[jk - 1] * zreflect
        )

    for jk in range(1, klev + 2):
        zreflect = 1.0 / (1.0 - prdnd[jk - 1] * prupd[jk - 1])
        pfu[_idx2(jk, kw, klev + 1)] = (
            ptdbt[jk - 1] * prup[jk - 1]
            + (ztdn[jk - 1] - ptdbt[jk - 1]) * prupd[jk - 1]
        ) * zreflect
        pfd[_idx2(jk, kw, klev + 1)] = (
            ptdbt[jk - 1]
            + (
                ztdn[jk - 1]
                - ptdbt[jk - 1]
                + ptdbt[jk - 1] * prup[jk - 1] * prdnd[jk - 1]
            )
            * zreflect
        )


@export
def rrtmg_sw_rad_setup_codon(
    nlay: int,
    ngptsw: int,
    nbndsw: int,
    icld: int,
    iaer: int,
    aldir_i: float,
    aldif_i: float,
    asdir_i: float,
    asdif_i: float,
    albdir_p: cobj,
    albdif_p: cobj,
    cldfmc_p: cobj,
    taucmc_p: cobj,
    taormc_p: cobj,
    asmcmc_p: cobj,
    ssacmc_p: cobj,
    zcldfmc_p: cobj,
    ztaucmc_p: cobj,
    ztaormc_p: cobj,
    zasycmc_p: cobj,
    zomgcmc_p: cobj,
    taua_p: cobj,
    ssaa_p: cobj,
    asma_p: cobj,
    ztaua_p: cobj,
    zasya_p: cobj,
    zomga_p: cobj,
):
    albdir = Ptr[float](albdir_p)
    albdif = Ptr[float](albdif_p)
    cldfmc = Ptr[float](cldfmc_p)
    taucmc = Ptr[float](taucmc_p)
    taormc = Ptr[float](taormc_p)
    asmcmc = Ptr[float](asmcmc_p)
    ssacmc = Ptr[float](ssacmc_p)
    zcldfmc = Ptr[float](zcldfmc_p)
    ztaucmc = Ptr[float](ztaucmc_p)
    ztaormc = Ptr[float](ztaormc_p)
    zasycmc = Ptr[float](zasycmc_p)
    zomgcmc = Ptr[float](zomgcmc_p)
    taua = Ptr[float](taua_p)
    ssaa = Ptr[float](ssaa_p)
    asma = Ptr[float](asma_p)
    ztaua = Ptr[float](ztaua_p)
    zasya = Ptr[float](zasya_p)
    zomga = Ptr[float](zomga_p)

    for ib in range(1, 9):
        albdir[ib - 1] = aldir_i
        albdif[ib - 1] = aldif_i
    albdir[nbndsw - 1] = aldir_i
    albdif[nbndsw - 1] = aldif_i
    albdir[8] = 0.5 * (aldir_i + asdir_i)
    albdif[8] = 0.5 * (aldif_i + asdif_i)
    for ib in range(10, 14):
        albdir[ib - 1] = asdir_i
        albdif[ib - 1] = asdif_i

    if icld == 0:
        for i in range(1, nlay + 1):
            for ig in range(1, ngptsw + 1):
                dst = _idx2(i, ig, nlay)
                zcldfmc[dst] = 0.0
                ztaucmc[dst] = 0.0
                ztaormc[dst] = 0.0
                zasycmc[dst] = 0.0
                zomgcmc[dst] = 1.0
    elif icld >= 1:
        for i in range(1, nlay + 1):
            for ig in range(1, ngptsw + 1):
                src = _idx2(ig, i, ngptsw)
                dst = _idx2(i, ig, nlay)
                zcldfmc[dst] = cldfmc[src]
                ztaucmc[dst] = taucmc[src]
                ztaormc[dst] = taormc[src]
                zasycmc[dst] = asmcmc[src]
                zomgcmc[dst] = ssacmc[src]

    if iaer == 0:
        for i in range(1, nlay + 1):
            for ib in range(1, nbndsw + 1):
                dst = _idx2(i, ib, nlay)
                ztaua[dst] = 0.0
                zasya[dst] = 0.0
                zomga[dst] = 1.0
    elif iaer == 10:
        for i in range(1, nlay + 1):
            for ib in range(1, nbndsw + 1):
                dst = _idx2(i, ib, nlay)
                ztaua[dst] = taua[dst]
                zasya[dst] = asma[dst]
                zomga[dst] = ssaa[dst]


@export
def rrtmg_sw_rad_zero_flux_codon(
    nlay: int,
    nbndsw: int,
    zbbcu_p: cobj,
    zbbcd_p: cobj,
    zbbfu_p: cobj,
    zbbfd_p: cobj,
    zbbcddir_p: cobj,
    zbbfddir_p: cobj,
    zuvcd_p: cobj,
    zuvfd_p: cobj,
    zuvcddir_p: cobj,
    zuvfddir_p: cobj,
    znicd_p: cobj,
    znifd_p: cobj,
    znicddir_p: cobj,
    znifddir_p: cobj,
    znicu_p: cobj,
    znifu_p: cobj,
    zbbfsu_p: cobj,
    zbbfsd_p: cobj,
):
    zbbcu = Ptr[float](zbbcu_p)
    zbbcd = Ptr[float](zbbcd_p)
    zbbfu = Ptr[float](zbbfu_p)
    zbbfd = Ptr[float](zbbfd_p)
    zbbcddir = Ptr[float](zbbcddir_p)
    zbbfddir = Ptr[float](zbbfddir_p)
    zuvcd = Ptr[float](zuvcd_p)
    zuvfd = Ptr[float](zuvfd_p)
    zuvcddir = Ptr[float](zuvcddir_p)
    zuvfddir = Ptr[float](zuvfddir_p)
    znicd = Ptr[float](znicd_p)
    znifd = Ptr[float](znifd_p)
    znicddir = Ptr[float](znicddir_p)
    znifddir = Ptr[float](znifddir_p)
    znicu = Ptr[float](znicu_p)
    znifu = Ptr[float](znifu_p)
    zbbfsu = Ptr[float](zbbfsu_p)
    zbbfsd = Ptr[float](zbbfsd_p)

    for i in range(1, nlay + 2):
        idx = i - 1
        zbbcu[idx] = 0.0
        zbbcd[idx] = 0.0
        zbbfu[idx] = 0.0
        zbbfd[idx] = 0.0
        zbbcddir[idx] = 0.0
        zbbfddir[idx] = 0.0
        zuvcd[idx] = 0.0
        zuvfd[idx] = 0.0
        zuvcddir[idx] = 0.0
        zuvfddir[idx] = 0.0
        znicd[idx] = 0.0
        znifd[idx] = 0.0
        znicddir[idx] = 0.0
        znifddir[idx] = 0.0
        znicu[idx] = 0.0
        znifu[idx] = 0.0
        for ib in range(1, nbndsw + 1):
            zbbfsu[_idx2(ib, i, nbndsw)] = 0.0
            zbbfsd[_idx2(ib, i, nbndsw)] = 0.0


@export
def rrtmg_sw_rad_store_flux_codon(
    nlay: int,
    nbndsw: int,
    ncol: int,
    iplon: int,
    heatfac: float,
    zbbcu_p: cobj,
    zbbcd_p: cobj,
    zbbfu_p: cobj,
    zbbfd_p: cobj,
    zbbcddir_p: cobj,
    zbbfddir_p: cobj,
    zuvfd_p: cobj,
    zuvfddir_p: cobj,
    znicd_p: cobj,
    znifd_p: cobj,
    znifddir_p: cobj,
    znicu_p: cobj,
    znifu_p: cobj,
    zbbfsu_p: cobj,
    zbbfsd_p: cobj,
    pdp_p: cobj,
    swuflxc_p: cobj,
    swdflxc_p: cobj,
    swuflx_p: cobj,
    swdflx_p: cobj,
    swuflxs_p: cobj,
    swdflxs_p: cobj,
    uvdflx_p: cobj,
    nidflx_p: cobj,
    dirdflux_p: cobj,
    difdflux_p: cobj,
    dirdnuv_p: cobj,
    difdnuv_p: cobj,
    dirdnir_p: cobj,
    difdnir_p: cobj,
    ninflx_p: cobj,
    ninflxc_p: cobj,
    swnflxc_p: cobj,
    swnflx_p: cobj,
    swhrc_p: cobj,
    swhr_p: cobj,
):
    zbbcu = Ptr[float](zbbcu_p)
    zbbcd = Ptr[float](zbbcd_p)
    zbbfu = Ptr[float](zbbfu_p)
    zbbfd = Ptr[float](zbbfd_p)
    zbbcddir = Ptr[float](zbbcddir_p)
    zbbfddir = Ptr[float](zbbfddir_p)
    zuvfd = Ptr[float](zuvfd_p)
    zuvfddir = Ptr[float](zuvfddir_p)
    znicd = Ptr[float](znicd_p)
    znifd = Ptr[float](znifd_p)
    znifddir = Ptr[float](znifddir_p)
    znicu = Ptr[float](znicu_p)
    znifu = Ptr[float](znifu_p)
    zbbfsu = Ptr[float](zbbfsu_p)
    zbbfsd = Ptr[float](zbbfsd_p)
    pdp = Ptr[float](pdp_p)
    swuflxc = Ptr[float](swuflxc_p)
    swdflxc = Ptr[float](swdflxc_p)
    swuflx = Ptr[float](swuflx_p)
    swdflx = Ptr[float](swdflx_p)
    swuflxs = Ptr[float](swuflxs_p)
    swdflxs = Ptr[float](swdflxs_p)
    uvdflx = Ptr[float](uvdflx_p)
    nidflx = Ptr[float](nidflx_p)
    dirdflux = Ptr[float](dirdflux_p)
    difdflux = Ptr[float](difdflux_p)
    dirdnuv = Ptr[float](dirdnuv_p)
    difdnuv = Ptr[float](difdnuv_p)
    dirdnir = Ptr[float](dirdnir_p)
    difdnir = Ptr[float](difdnir_p)
    ninflx = Ptr[float](ninflx_p)
    ninflxc = Ptr[float](ninflxc_p)
    swnflxc = Ptr[float](swnflxc_p)
    swnflx = Ptr[float](swnflx_p)
    swhrc = Ptr[float](swhrc_p)
    swhr = Ptr[float](swhr_p)

    for i in range(1, nlay + 2):
        jidx = i - 1
        dst2 = _idx2(iplon, i, ncol)
        swuflxc[dst2] = zbbcu[jidx]
        swdflxc[dst2] = zbbcd[jidx]
        swuflx[dst2] = zbbfu[jidx]
        swdflx[dst2] = zbbfd[jidx]
        for ib in range(1, nbndsw + 1):
            swuflxs[_idx3(ib, iplon, i, nbndsw, ncol)] = zbbfsu[_idx2(ib, i, nbndsw)]
            swdflxs[_idx3(ib, iplon, i, nbndsw, ncol)] = zbbfsd[_idx2(ib, i, nbndsw)]
        uvdflx[jidx] = zuvfd[jidx]
        nidflx[jidx] = znifd[jidx]
        dirdflux[jidx] = zbbfddir[jidx]
        difdflux[jidx] = swdflx[dst2] - dirdflux[jidx]
        dirdnuv[dst2] = zuvfddir[jidx]
        difdnuv[dst2] = zuvfd[jidx] - dirdnuv[dst2]
        dirdnir[dst2] = znifddir[jidx]
        difdnir[dst2] = znifd[jidx] - dirdnir[dst2]
        ninflx[dst2] = znifd[jidx] - znifu[jidx]
        ninflxc[dst2] = znicd[jidx] - znicu[jidx]

    for i in range(1, nlay + 2):
        dst2 = _idx2(iplon, i, ncol)
        swnflxc[i - 1] = swdflxc[dst2] - swuflxc[dst2]
        swnflx[i - 1] = swdflx[dst2] - swuflx[dst2]

    for i in range(1, nlay + 1):
        zdpgcp = heatfac / pdp[i - 1]
        dst2 = _idx2(iplon, i, ncol)
        swhrc[dst2] = (swnflxc[i] - swnflxc[i - 1]) * zdpgcp
        swhr[dst2] = (swnflx[i] - swnflx[i - 1]) * zdpgcp
    swhrc[_idx2(iplon, nlay, ncol)] = 0.0
    swhr[_idx2(iplon, nlay, ncol)] = 0.0


@export
def rrtmg_lw_rad_taut_codon(
    nlay: int,
    ngptlw: int,
    iaer: int,
    taug_p: cobj,
    taua_p: cobj,
    ngb_p: cobj,
    taut_p: cobj,
):
    taug = Ptr[float](taug_p)
    taua = Ptr[float](taua_p)
    ngb = Ptr[int](ngb_p)
    taut = Ptr[float](taut_p)

    if iaer == 0:
        for k in range(1, nlay + 1):
            for ig in range(1, ngptlw + 1):
                idx = _idx2(k, ig, nlay)
                taut[idx] = taug[idx]
    elif iaer == 10:
        for k in range(1, nlay + 1):
            for ig in range(1, ngptlw + 1):
                idx = _idx2(k, ig, nlay)
                taut[idx] = taug[idx] + taua[_idx2(k, ngb[ig - 1], nlay)]


@export
def rrtmg_lw_rad_store_flux_codon(
    nlay: int,
    nbndlw: int,
    ncol: int,
    iplon: int,
    totuflux_p: cobj,
    totdflux_p: cobj,
    totuclfl_p: cobj,
    totdclfl_p: cobj,
    totufluxs_p: cobj,
    totdfluxs_p: cobj,
    htr_p: cobj,
    htrc_p: cobj,
    uflx_p: cobj,
    dflx_p: cobj,
    uflxc_p: cobj,
    dflxc_p: cobj,
    uflxs_p: cobj,
    dflxs_p: cobj,
    hr_p: cobj,
    hrc_p: cobj,
):
    totuflux = Ptr[float](totuflux_p)
    totdflux = Ptr[float](totdflux_p)
    totuclfl = Ptr[float](totuclfl_p)
    totdclfl = Ptr[float](totdclfl_p)
    totufluxs = Ptr[float](totufluxs_p)
    totdfluxs = Ptr[float](totdfluxs_p)
    htr = Ptr[float](htr_p)
    htrc = Ptr[float](htrc_p)
    uflx = Ptr[float](uflx_p)
    dflx = Ptr[float](dflx_p)
    uflxc = Ptr[float](uflxc_p)
    dflxc = Ptr[float](dflxc_p)
    uflxs = Ptr[float](uflxs_p)
    dflxs = Ptr[float](dflxs_p)
    hr = Ptr[float](hr_p)
    hrc = Ptr[float](hrc_p)

    for k0 in range(0, nlay + 1):
        dst2 = _idx2(iplon, k0 + 1, ncol)
        uflx[dst2] = totuflux[k0]
        dflx[dst2] = totdflux[k0]
        uflxc[dst2] = totuclfl[k0]
        dflxc[dst2] = totdclfl[k0]
        for ib in range(1, nbndlw + 1):
            uflxs[_idx3(ib, iplon, k0 + 1, nbndlw, ncol)] = totufluxs[_idx2_dim2_lb(ib, k0, nbndlw, 0)]
            dflxs[_idx3(ib, iplon, k0 + 1, nbndlw, ncol)] = totdfluxs[_idx2_dim2_lb(ib, k0, nbndlw, 0)]

    for k0 in range(0, nlay):
        dst2 = _idx2(iplon, k0 + 1, ncol)
        hr[dst2] = htr[k0]
        hrc[dst2] = htrc[k0]


@export
def rrtmg_sw_inatm_codon(
    iplon: int,
    nlay: int,
    ldcol: int,
    icld: int,
    iaer: int,
    nbndsw: int,
    ngptsw: int,
    nmol: int,
    mxmol: int,
    jpband: int,
    jpb1: int,
    jpb2: int,
    grav: float,
    avogad: float,
    adjflx: float,
    play_p: cobj,
    plev_p: cobj,
    tlay_p: cobj,
    tlev_p: cobj,
    tsfc_p: cobj,
    h2ovmr_p: cobj,
    o3vmr_p: cobj,
    co2vmr_p: cobj,
    ch4vmr_p: cobj,
    o2vmr_p: cobj,
    n2ovmr_p: cobj,
    solvar_p: cobj,
    inflgsw: int,
    iceflgsw: int,
    liqflgsw: int,
    cldfmcl_p: cobj,
    taucmcl_p: cobj,
    ssacmcl_p: cobj,
    asmcmcl_p: cobj,
    fsfcmcl_p: cobj,
    ciwpmcl_p: cobj,
    clwpmcl_p: cobj,
    reicmcl_p: cobj,
    relqmcl_p: cobj,
    tauaer_p: cobj,
    ssaaer_p: cobj,
    asmaer_p: cobj,
    pavel_p: cobj,
    pz_p: cobj,
    pdp_p: cobj,
    tavel_p: cobj,
    tz_p: cobj,
    tbound_p: cobj,
    coldry_p: cobj,
    wkl_p: cobj,
    adjflux_p: cobj,
    inflag_p: cobj,
    iceflag_p: cobj,
    liqflag_p: cobj,
    cldfmc_p: cobj,
    taucmc_p: cobj,
    ssacmc_p: cobj,
    asmcmc_p: cobj,
    fsfcmc_p: cobj,
    ciwpmc_p: cobj,
    clwpmc_p: cobj,
    reicmc_p: cobj,
    dgesmc_p: cobj,
    relqmc_p: cobj,
    taua_p: cobj,
    ssaa_p: cobj,
    asma_p: cobj,
):
    play = Ptr[float](play_p)
    plev = Ptr[float](plev_p)
    tlay = Ptr[float](tlay_p)
    tlev = Ptr[float](tlev_p)
    tsfc = Ptr[float](tsfc_p)
    h2ovmr = Ptr[float](h2ovmr_p)
    o3vmr = Ptr[float](o3vmr_p)
    co2vmr = Ptr[float](co2vmr_p)
    ch4vmr = Ptr[float](ch4vmr_p)
    o2vmr = Ptr[float](o2vmr_p)
    n2ovmr = Ptr[float](n2ovmr_p)
    solvar = Ptr[float](solvar_p)
    cldfmcl = Ptr[float](cldfmcl_p)
    taucmcl = Ptr[float](taucmcl_p)
    ssacmcl = Ptr[float](ssacmcl_p)
    asmcmcl = Ptr[float](asmcmcl_p)
    fsfcmcl = Ptr[float](fsfcmcl_p)
    ciwpmcl = Ptr[float](ciwpmcl_p)
    clwpmcl = Ptr[float](clwpmcl_p)
    reicmcl = Ptr[float](reicmcl_p)
    relqmcl = Ptr[float](relqmcl_p)
    tauaer = Ptr[float](tauaer_p)
    ssaaer = Ptr[float](ssaaer_p)
    asmaer = Ptr[float](asmaer_p)
    pavel = Ptr[float](pavel_p)
    pz = Ptr[float](pz_p)
    pdp = Ptr[float](pdp_p)
    tavel = Ptr[float](tavel_p)
    tz = Ptr[float](tz_p)
    tbound = Ptr[float](tbound_p)
    coldry = Ptr[float](coldry_p)
    wkl = Ptr[float](wkl_p)
    adjflux = Ptr[float](adjflux_p)
    inflag = Ptr[int](inflag_p)
    iceflag = Ptr[int](iceflag_p)
    liqflag = Ptr[int](liqflag_p)
    cldfmc = Ptr[float](cldfmc_p)
    taucmc = Ptr[float](taucmc_p)
    ssacmc = Ptr[float](ssacmc_p)
    asmcmc = Ptr[float](asmcmc_p)
    fsfcmc = Ptr[float](fsfcmc_p)
    ciwpmc = Ptr[float](ciwpmc_p)
    clwpmc = Ptr[float](clwpmc_p)
    reicmc = Ptr[float](reicmc_p)
    dgesmc = Ptr[float](dgesmc_p)
    relqmc = Ptr[float](relqmc_p)
    taua = Ptr[float](taua_p)
    ssaa = Ptr[float](ssaa_p)
    asma = Ptr[float](asma_p)

    amd = 28.9660
    amw = 18.0160

    for l in range(1, nlay + 1):
        for imol in range(1, mxmol + 1):
            wkl[_idx2(imol, l, mxmol)] = 0.0
        reicmc[l - 1] = 0.0
        dgesmc[l - 1] = 0.0
        relqmc[l - 1] = 0.0
        for ig in range(1, ngptsw + 1):
            cldfmc[_idx2(ig, l, ngptsw)] = 0.0
            taucmc[_idx2(ig, l, ngptsw)] = 0.0
            ssacmc[_idx2(ig, l, ngptsw)] = 1.0
            asmcmc[_idx2(ig, l, ngptsw)] = 0.0
            fsfcmc[_idx2(ig, l, ngptsw)] = 0.0
            ciwpmc[_idx2(ig, l, ngptsw)] = 0.0
            clwpmc[_idx2(ig, l, ngptsw)] = 0.0
        for ib in range(1, nbndsw + 1):
            taua[_idx2(l, ib, nlay)] = 0.0
            ssaa[_idx2(l, ib, nlay)] = 1.0
            asma[_idx2(l, ib, nlay)] = 0.0

    for ib in range(1, jpband + 1):
        adjflux[ib - 1] = 0.0
    for ib in range(jpb1, jpb2 + 1):
        adjflux[ib - 1] = adjflx * solvar[ib - jpb1]

    tbound[0] = tsfc[iplon - 1]

    pz[0] = plev[_idx2(iplon, nlay + 1, ldcol)]
    tz[0] = tlev[_idx2(iplon, nlay + 1, ldcol)]
    amm = 0.0
    for l in range(1, nlay + 1):
        src_l = nlay - l + 1
        pavel[l - 1] = play[_idx2(iplon, src_l, ldcol)]
        tavel[l - 1] = tlay[_idx2(iplon, src_l, ldcol)]
        pz[l] = plev[_idx2(iplon, src_l, ldcol)]
        tz[l] = tlev[_idx2(iplon, src_l, ldcol)]
        pdp[l - 1] = pz[l - 1] - pz[l]
        wkl[_idx2(1, l, mxmol)] = h2ovmr[_idx2(iplon, src_l, ldcol)]
        wkl[_idx2(2, l, mxmol)] = co2vmr[_idx2(iplon, src_l, ldcol)]
        wkl[_idx2(3, l, mxmol)] = o3vmr[_idx2(iplon, src_l, ldcol)]
        wkl[_idx2(4, l, mxmol)] = n2ovmr[_idx2(iplon, src_l, ldcol)]
        wkl[_idx2(6, l, mxmol)] = ch4vmr[_idx2(iplon, src_l, ldcol)]
        wkl[_idx2(7, l, mxmol)] = o2vmr[_idx2(iplon, src_l, ldcol)]
        amm = (1.0 - wkl[_idx2(1, l, mxmol)]) * amd + wkl[_idx2(1, l, mxmol)] * amw
        coldry[l - 1] = (
            (pz[l - 1] - pz[l])
            * 1.0e3
            * avogad
            / (1.0e2 * grav * amm * (1.0 + wkl[_idx2(1, l, mxmol)]))
        )

    coldry[nlay - 1] = (
        pz[nlay - 1]
        * 1.0e3
        * avogad
        / (1.0e2 * grav * amm * (1.0 + wkl[_idx2(1, nlay - 1, mxmol)]))
    )

    for l in range(1, nlay + 1):
        for imol in range(1, nmol + 1):
            wkl[_idx2(imol, l, mxmol)] = coldry[l - 1] * wkl[_idx2(imol, l, mxmol)]

    if iaer >= 1:
        for l in range(1, nlay):
            src_l = nlay - l
            for ib in range(1, nbndsw + 1):
                taua[_idx2(l, ib, nlay)] = tauaer[_idx3(iplon, src_l, ib, ldcol, nlay - 1)]
                ssaa[_idx2(l, ib, nlay)] = ssaaer[_idx3(iplon, src_l, ib, ldcol, nlay - 1)]
                asma[_idx2(l, ib, nlay)] = asmaer[_idx3(iplon, src_l, ib, ldcol, nlay - 1)]

    if icld >= 1:
        inflag[0] = inflgsw
        iceflag[0] = iceflgsw
        liqflag[0] = liqflgsw
        for l in range(1, nlay):
            src_l = nlay - l
            for ig in range(1, ngptsw + 1):
                cldfmc[_idx2(ig, l, ngptsw)] = cldfmcl[_idx3(ig, iplon, src_l, ngptsw, ldcol)]
                taucmc[_idx2(ig, l, ngptsw)] = taucmcl[_idx3(ig, iplon, src_l, ngptsw, ldcol)]
                ssacmc[_idx2(ig, l, ngptsw)] = ssacmcl[_idx3(ig, iplon, src_l, ngptsw, ldcol)]
                asmcmc[_idx2(ig, l, ngptsw)] = asmcmcl[_idx3(ig, iplon, src_l, ngptsw, ldcol)]
                fsfcmc[_idx2(ig, l, ngptsw)] = fsfcmcl[_idx3(ig, iplon, src_l, ngptsw, ldcol)]
                ciwpmc[_idx2(ig, l, ngptsw)] = ciwpmcl[_idx3(ig, iplon, src_l, ngptsw, ldcol)]
                clwpmc[_idx2(ig, l, ngptsw)] = clwpmcl[_idx3(ig, iplon, src_l, ngptsw, ldcol)]
            reicmc[l - 1] = reicmcl[_idx2(iplon, src_l, ldcol)]
            if iceflgsw == 3:
                dgesmc[l - 1] = 1.5396 * reicmcl[_idx2(iplon, src_l, ldcol)]
            relqmc[l - 1] = relqmcl[_idx2(iplon, src_l, ldcol)]

        for ig in range(1, ngptsw + 1):
            cldfmc[_idx2(ig, nlay, ngptsw)] = 0.0
            taucmc[_idx2(ig, nlay, ngptsw)] = 0.0
            ssacmc[_idx2(ig, nlay, ngptsw)] = 1.0
            asmcmc[_idx2(ig, nlay, ngptsw)] = 0.0
            fsfcmc[_idx2(ig, nlay, ngptsw)] = 0.0
            ciwpmc[_idx2(ig, nlay, ngptsw)] = 0.0
            clwpmc[_idx2(ig, nlay, ngptsw)] = 0.0
        reicmc[nlay - 1] = 0.0
        dgesmc[nlay - 1] = 0.0
        relqmc[nlay - 1] = 0.0
        for ib in range(1, nbndsw + 1):
            taua[_idx2(nlay, ib, nlay)] = 0.0
            ssaa[_idx2(nlay, ib, nlay)] = 1.0
            asma[_idx2(nlay, ib, nlay)] = 0.0


@export
def rrtmg_lw_inatm_codon(
    iplon: int,
    nlay: int,
    ldcol: int,
    icld: int,
    iaer: int,
    nbndlw: int,
    ngptlw: int,
    nmol: int,
    maxxsec: int,
    mxmol: int,
    grav: float,
    avogad: float,
    play_p: cobj,
    plev_p: cobj,
    tlay_p: cobj,
    tlev_p: cobj,
    tsfc_p: cobj,
    h2ovmr_p: cobj,
    o3vmr_p: cobj,
    co2vmr_p: cobj,
    ch4vmr_p: cobj,
    o2vmr_p: cobj,
    n2ovmr_p: cobj,
    cfc11vmr_p: cobj,
    cfc12vmr_p: cobj,
    cfc22vmr_p: cobj,
    ccl4vmr_p: cobj,
    emis_p: cobj,
    inflglw: int,
    iceflglw: int,
    liqflglw: int,
    cldfmcl_p: cobj,
    taucmcl_p: cobj,
    ciwpmcl_p: cobj,
    clwpmcl_p: cobj,
    reicmcl_p: cobj,
    relqmcl_p: cobj,
    tauaer_p: cobj,
    pavel_p: cobj,
    pz_p: cobj,
    tavel_p: cobj,
    tz_p: cobj,
    tbound_p: cobj,
    semiss_p: cobj,
    coldry_p: cobj,
    wbrodl_p: cobj,
    wkl_p: cobj,
    wx_p: cobj,
    pwvcm_p: cobj,
    inflag_p: cobj,
    iceflag_p: cobj,
    liqflag_p: cobj,
    cldfmc_p: cobj,
    taucmc_p: cobj,
    ciwpmc_p: cobj,
    clwpmc_p: cobj,
    reicmc_p: cobj,
    dgesmc_p: cobj,
    relqmc_p: cobj,
    taua_p: cobj,
    ixindx_p: cobj,
):
    play = Ptr[float](play_p)
    plev = Ptr[float](plev_p)
    tlay = Ptr[float](tlay_p)
    tlev = Ptr[float](tlev_p)
    tsfc = Ptr[float](tsfc_p)
    h2ovmr = Ptr[float](h2ovmr_p)
    o3vmr = Ptr[float](o3vmr_p)
    co2vmr = Ptr[float](co2vmr_p)
    ch4vmr = Ptr[float](ch4vmr_p)
    o2vmr = Ptr[float](o2vmr_p)
    n2ovmr = Ptr[float](n2ovmr_p)
    cfc11vmr = Ptr[float](cfc11vmr_p)
    cfc12vmr = Ptr[float](cfc12vmr_p)
    cfc22vmr = Ptr[float](cfc22vmr_p)
    ccl4vmr = Ptr[float](ccl4vmr_p)
    emis = Ptr[float](emis_p)
    cldfmcl = Ptr[float](cldfmcl_p)
    taucmcl = Ptr[float](taucmcl_p)
    ciwpmcl = Ptr[float](ciwpmcl_p)
    clwpmcl = Ptr[float](clwpmcl_p)
    reicmcl = Ptr[float](reicmcl_p)
    relqmcl = Ptr[float](relqmcl_p)
    tauaer = Ptr[float](tauaer_p)
    pavel = Ptr[float](pavel_p)
    pz = Ptr[float](pz_p)
    tavel = Ptr[float](tavel_p)
    tz = Ptr[float](tz_p)
    tbound = Ptr[float](tbound_p)
    semiss = Ptr[float](semiss_p)
    coldry = Ptr[float](coldry_p)
    wbrodl = Ptr[float](wbrodl_p)
    wkl = Ptr[float](wkl_p)
    wx = Ptr[float](wx_p)
    pwvcm = Ptr[float](pwvcm_p)
    inflag = Ptr[int](inflag_p)
    iceflag = Ptr[int](iceflag_p)
    liqflag = Ptr[int](liqflag_p)
    cldfmc = Ptr[float](cldfmc_p)
    taucmc = Ptr[float](taucmc_p)
    ciwpmc = Ptr[float](ciwpmc_p)
    clwpmc = Ptr[float](clwpmc_p)
    reicmc = Ptr[float](reicmc_p)
    dgesmc = Ptr[float](dgesmc_p)
    relqmc = Ptr[float](relqmc_p)
    taua = Ptr[float](taua_p)
    ixindx = Ptr[int](ixindx_p)

    amd = 28.9660
    amw = 18.0160

    for l in range(1, nlay + 1):
        for imol in range(1, mxmol + 1):
            wkl[_idx2(imol, l, mxmol)] = 0.0
        for ix in range(1, maxxsec + 1):
            wx[_idx2(ix, l, maxxsec)] = 0.0
        reicmc[l - 1] = 0.0
        dgesmc[l - 1] = 0.0
        relqmc[l - 1] = 0.0
        for ig in range(1, ngptlw + 1):
            cldfmc[_idx2(ig, l, ngptlw)] = 0.0
            taucmc[_idx2(ig, l, ngptlw)] = 0.0
            ciwpmc[_idx2(ig, l, ngptlw)] = 0.0
            clwpmc[_idx2(ig, l, ngptlw)] = 0.0
        for ib in range(1, nbndlw + 1):
            taua[_idx2(l, ib, nlay)] = 0.0

    amttl = 0.0
    wvttl = 0.0
    tbound[0] = tsfc[iplon - 1]

    pz[0] = plev[_idx2(iplon, nlay + 1, ldcol)]
    tz[0] = tlev[_idx2(iplon, nlay + 1, ldcol)]
    amm = 0.0
    for l in range(1, nlay + 1):
        src_l = nlay - l + 1
        pavel[l - 1] = play[_idx2(iplon, src_l, ldcol)]
        tavel[l - 1] = tlay[_idx2(iplon, src_l, ldcol)]
        pz[l] = plev[_idx2(iplon, src_l, ldcol)]
        tz[l] = tlev[_idx2(iplon, src_l, ldcol)]
        wkl[_idx2(1, l, mxmol)] = h2ovmr[_idx2(iplon, src_l, ldcol)]
        wkl[_idx2(2, l, mxmol)] = co2vmr[_idx2(iplon, src_l, ldcol)]
        wkl[_idx2(3, l, mxmol)] = o3vmr[_idx2(iplon, src_l, ldcol)]
        wkl[_idx2(4, l, mxmol)] = n2ovmr[_idx2(iplon, src_l, ldcol)]
        wkl[_idx2(6, l, mxmol)] = ch4vmr[_idx2(iplon, src_l, ldcol)]
        wkl[_idx2(7, l, mxmol)] = o2vmr[_idx2(iplon, src_l, ldcol)]
        amm = (1.0 - wkl[_idx2(1, l, mxmol)]) * amd + wkl[_idx2(1, l, mxmol)] * amw
        coldry[l - 1] = (
            (pz[l - 1] - pz[l])
            * 1.0e3
            * avogad
            / (1.0e2 * grav * amm * (1.0 + wkl[_idx2(1, l, mxmol)]))
        )
        wx[_idx2(1, l, maxxsec)] = ccl4vmr[_idx2(iplon, src_l, ldcol)]
        wx[_idx2(2, l, maxxsec)] = cfc11vmr[_idx2(iplon, src_l, ldcol)]
        wx[_idx2(3, l, maxxsec)] = cfc12vmr[_idx2(iplon, src_l, ldcol)]
        wx[_idx2(4, l, maxxsec)] = cfc22vmr[_idx2(iplon, src_l, ldcol)]

    coldry[nlay - 1] = (
        pz[nlay - 1]
        * 1.0e3
        * avogad
        / (1.0e2 * grav * amm * (1.0 + wkl[_idx2(1, nlay - 1, mxmol)]))
    )

    for l in range(1, nlay + 1):
        summol = 0.0
        for imol in range(2, nmol + 1):
            summol = summol + wkl[_idx2(imol, l, mxmol)]
        wbrodl[l - 1] = coldry[l - 1] * (1.0 - summol)
        for imol in range(1, nmol + 1):
            wkl[_idx2(imol, l, mxmol)] = coldry[l - 1] * wkl[_idx2(imol, l, mxmol)]
        amttl = amttl + coldry[l - 1] + wkl[_idx2(1, l, mxmol)]
        wvttl = wvttl + wkl[_idx2(1, l, mxmol)]
        for ix in range(1, maxxsec + 1):
            ixout = ixindx[ix - 1]
            if ixout != 0:
                wx[_idx2(ixout, l, maxxsec)] = coldry[l - 1] * wx[_idx2(ix, l, maxxsec)] * 1.0e-20

    wvsh = (amw * wvttl) / (amd * amttl)
    pwvcm[0] = wvsh * (1.0e3 * pz[0]) / (1.0e2 * grav)

    for n in range(1, nbndlw + 1):
        semiss[n - 1] = emis[_idx2(iplon, n, ldcol)]

    if iaer >= 1:
        for l in range(1, nlay):
            src_l = nlay - l
            for ib in range(1, nbndlw + 1):
                taua[_idx2(l, ib, nlay)] = tauaer[_idx3(iplon, src_l, ib, ldcol, nlay - 1)]

    if icld >= 1:
        inflag[0] = inflglw
        iceflag[0] = iceflglw
        liqflag[0] = liqflglw
        for l in range(1, nlay):
            src_l = nlay - l
            for ig in range(1, ngptlw + 1):
                cldfmc[_idx2(ig, l, ngptlw)] = cldfmcl[_idx3(ig, iplon, src_l, ngptlw, ldcol)]
                taucmc[_idx2(ig, l, ngptlw)] = taucmcl[_idx3(ig, iplon, src_l, ngptlw, ldcol)]
                ciwpmc[_idx2(ig, l, ngptlw)] = ciwpmcl[_idx3(ig, iplon, src_l, ngptlw, ldcol)]
                clwpmc[_idx2(ig, l, ngptlw)] = clwpmcl[_idx3(ig, iplon, src_l, ngptlw, ldcol)]
            reicmc[l - 1] = reicmcl[_idx2(iplon, src_l, ldcol)]
            if iceflglw == 3:
                dgesmc[l - 1] = 1.5396 * reicmcl[_idx2(iplon, src_l, ldcol)]
            relqmc[l - 1] = relqmcl[_idx2(iplon, src_l, ldcol)]

        for ig in range(1, ngptlw + 1):
            cldfmc[_idx2(ig, nlay, ngptlw)] = 0.0
            taucmc[_idx2(ig, nlay, ngptlw)] = 0.0
            ciwpmc[_idx2(ig, nlay, ngptlw)] = 0.0
            clwpmc[_idx2(ig, nlay, ngptlw)] = 0.0
        reicmc[nlay - 1] = 0.0
        dgesmc[nlay - 1] = 0.0
        relqmc[nlay - 1] = 0.0
        for ib in range(1, nbndlw + 1):
            taua[_idx2(nlay, ib, nlay)] = 0.0
