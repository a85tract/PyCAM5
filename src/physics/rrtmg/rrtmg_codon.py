from math import sqrt


@inline
def _idx2(i: int, k: int, ld1: int) -> int:
    """Fortran array declared as (ld1, n2), both bounds starting at 1."""
    return (i - 1) + (k - 1) * ld1


@inline
def _idx2_lb0(i: int, k0: int, ld1: int) -> int:
    """Fortran array declared as (ld1, 0:n2)."""
    return (i - 1) + k0 * ld1


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
