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
