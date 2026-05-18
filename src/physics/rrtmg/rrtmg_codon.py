from math import exp, log, sqrt


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
def _idx3_lb0_dim2(a: int, b0: int, c: int, ld1: int, ub2: int) -> int:
    """Fortran array declared as (ld1, 0:ub2, n3)."""
    return (a - 1) + b0 * ld1 + (c - 1) * ld1 * (ub2 + 1)


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


@export
def rrtmg_sw_pre_codon(
    nday: int,
    pcols: int,
    pver: int,
    pverp: int,
    rrtmg_levs: int,
    nbndsw: int,
    e_aer_tau_p: cobj,
    e_aer_tau_w_p: cobj,
    e_aer_tau_w_g_p: cobj,
    idxday_p: cobj,
    tau_aer_sw_p: cobj,
    ssa_aer_sw_p: cobj,
    asm_aer_sw_p: cobj,
    tlev_p: cobj,
    sfac_p: cobj,
    tsfc_p: cobj,
    solvar_p: cobj,
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
    tsfc = Ptr[float](tsfc_p)
    solvar = Ptr[float](solvar_p)

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
