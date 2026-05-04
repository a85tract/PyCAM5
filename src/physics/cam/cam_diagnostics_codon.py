from math import sqrt


@inline
def _idx(i: int) -> int:
    return i - 1


@inline
def _idx2(i: int, k: int, ld1: int) -> int:
    return (i - 1) + (k - 1) * ld1


@inline
def _idx3(i: int, k: int, m: int, ld1: int, ld2: int) -> int:
    return (i - 1) + (k - 1) * ld1 + (m - 1) * ld1 * ld2


@export
def diag_conv_tend_ini_copy_s_codon(
    ncol: int,
    pcols: int,
    pver: int,
    state_s_p: cobj,
    dtcond_p: cobj,
):
    state_s = Ptr[float](state_s_p)
    dtcond = Ptr[float](dtcond_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            dtcond[_idx2(i, k, pcols)] = state_s[_idx2(i, k, pcols)]


@export
def diag_conv_tend_ini_copy_q_m_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    m: int,
    state_q_p: cobj,
    dqcond_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    dqcond = Ptr[float](dqcond_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            dqcond[_idx2(i, k, pcols)] = state_q[_idx3(i, k, m, pcols, pver)]


@export
def diag_conv_tend_ini_copy_2d_codon(
    ncol: int,
    pcols: int,
    pver: int,
    src_p: cobj,
    dst_p: cobj,
):
    src = Ptr[float](src_p)
    dst = Ptr[float](dst_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            dst[_idx2(i, k, pcols)] = src[_idx2(i, k, pcols)]


@export
def diag_conv_precip_codon(
    ncol: int,
    pcols: int,
    prec_dp_p: cobj,
    snow_dp_p: cobj,
    prec_sh_p: cobj,
    snow_sh_p: cobj,
    prec_sed_p: cobj,
    snow_sed_p: cobj,
    prec_pcw_p: cobj,
    snow_pcw_p: cobj,
    precc_p: cobj,
    precl_p: cobj,
    snowc_p: cobj,
    snowl_p: cobj,
    prect_p: cobj,
):
    prec_dp = Ptr[float](prec_dp_p)
    snow_dp = Ptr[float](snow_dp_p)
    prec_sh = Ptr[float](prec_sh_p)
    snow_sh = Ptr[float](snow_sh_p)
    prec_sed = Ptr[float](prec_sed_p)
    snow_sed = Ptr[float](snow_sed_p)
    prec_pcw = Ptr[float](prec_pcw_p)
    snow_pcw = Ptr[float](snow_pcw_p)
    precc = Ptr[float](precc_p)
    precl = Ptr[float](precl_p)
    snowc = Ptr[float](snowc_p)
    snowl = Ptr[float](snowl_p)
    prect = Ptr[float](prect_p)

    for i in range(1, ncol + 1):
        idx = _idx(i)
        precc[idx] = prec_dp[idx] + prec_sh[idx]
        precl[idx] = prec_sed[idx] + prec_pcw[idx]
        snowc[idx] = snow_dp[idx] + snow_sh[idx]
        snowl[idx] = snow_sed[idx] + snow_pcw[idx]
        prect[idx] = precc[idx] + precl[idx]


@export
def diag_conv_wtprect_codon(
    ncol: int,
    pcols: int,
    wtprec1_p: cobj,
    wtprec2_p: cobj,
    wtprec3_p: cobj,
    wtprec4_p: cobj,
    wtprect_p: cobj,
):
    wtprec1 = Ptr[float](wtprec1_p)
    wtprec2 = Ptr[float](wtprec2_p)
    wtprec3 = Ptr[float](wtprec3_p)
    wtprec4 = Ptr[float](wtprec4_p)
    wtprect = Ptr[float](wtprect_p)

    for i in range(1, ncol + 1):
        idx = _idx(i)
        wtprect[idx] = wtprec1[idx]
        wtprect[idx] = wtprect[idx] + wtprec2[idx]
        wtprect[idx] = wtprect[idx] + wtprec3[idx]
        wtprect[idx] = wtprect[idx] + wtprec4[idx]


@export
def diag_conv_dtcond_codon(
    ncol: int,
    pcols: int,
    pver: int,
    rtdt: float,
    cpair: float,
    state_s_p: cobj,
    dtcond_p: cobj,
):
    state_s = Ptr[float](state_s_p)
    dtcond = Ptr[float](dtcond_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            dtcond[idx] = (state_s[idx] - dtcond[idx]) * rtdt / cpair


@export
def diag_conv_dqcond_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    m: int,
    rtdt: float,
    state_q_p: cobj,
    dqcond_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    dqcond = Ptr[float](dqcond_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            dqcond[idx] = (state_q[_idx3(i, k, m, pcols, pver)] - dqcond[idx]) * rtdt


@export
def diag_surf_codon(
    ncol: int,
    pcols: int,
    end_day: int,
    qref_p: cobj,
    rhref_p: cobj,
    tref_p: cobj,
    trefmxav_p: cobj,
    trefmnav_p: cobj,
    trefmx_day_p: cobj,
    trefmn_day_p: cobj,
):
    qref = Ptr[float](qref_p)
    rhref = Ptr[float](rhref_p)
    tref = Ptr[float](tref_p)
    trefmxav = Ptr[float](trefmxav_p)
    trefmnav = Ptr[float](trefmnav_p)
    trefmx_day = Ptr[float](trefmx_day_p)
    trefmn_day = Ptr[float](trefmn_day_p)

    for i in range(1, ncol + 1):
        rhref[_idx(i)] = qref[_idx(i)] / rhref[_idx(i)] * 100.0
        trefmxav[_idx(i)] = max(tref[_idx(i)], trefmxav[_idx(i)])
        trefmnav[_idx(i)] = min(tref[_idx(i)], trefmnav[_idx(i)])

    if end_day != 0:
        for i in range(1, ncol + 1):
            trefmx_day[_idx(i)] = trefmxav[_idx(i)]
            trefmn_day[_idx(i)] = trefmnav[_idx(i)]
            trefmxav[_idx(i)] = -1.0e36
            trefmnav[_idx(i)] = 1.0e36


@export
def diag_physvar_ic_codon():
    return


@export
def diag_phys_writeout_z3_codon(
    ncol: int,
    pcols: int,
    pver: int,
    rga: float,
    zm_p: cobj,
    phis_p: cobj,
    z3_p: cobj,
):
    zm = Ptr[float](zm_p)
    phis = Ptr[float](phis_p)
    z3 = Ptr[float](z3_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            z3[_idx2(i, k, pcols)] = zm[_idx2(i, k, pcols)] + phis[_idx(i)] * rga


@export
def diag_phys_writeout_basic_2d_codon(
    mode: int,
    ncol: int,
    pcols: int,
    pver: int,
    scale: float,
    a_p: cobj,
    b_p: cobj,
    out_p: cobj,
):
    a = Ptr[float](a_p)
    b = Ptr[float](b_p)
    out = Ptr[float](out_p)

    if mode == 1:
        # z3 = zm + phis*rga; b is phis(pcols) for this mode.
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                out[_idx2(i, k, pcols)] = a[_idx2(i, k, pcols)] + b[_idx(i)] * scale
    elif mode == 2:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                out[_idx2(i, k, pcols)] = a[_idx2(i, k, pcols)] * b[_idx2(i, k, pcols)]
    elif mode == 3:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                out[_idx2(i, k, pcols)] = (
                    a[_idx2(i, k, pcols)] * b[_idx2(i, k, pcols)] * scale
                )
    elif mode == 4:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                val = a[_idx2(i, k, pcols)]
                out[_idx2(i, k, pcols)] = val * val
    elif mode == 5:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                uval = a[_idx2(i, k, pcols)]
                vval = b[_idx2(i, k, pcols)]
                out[_idx2(i, k, pcols)] = sqrt(uval * uval + vval * vval)


@export
def diag_phys_writeout_mul_codon(
    ncol: int,
    pcols: int,
    pver: int,
    a_p: cobj,
    b_p: cobj,
    out_p: cobj,
):
    a = Ptr[float](a_p)
    b = Ptr[float](b_p)
    out = Ptr[float](out_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            out[_idx2(i, k, pcols)] = a[_idx2(i, k, pcols)] * b[_idx2(i, k, pcols)]


@export
def diag_phys_writeout_mul_scalar_codon(
    ncol: int,
    pcols: int,
    pver: int,
    scale: float,
    a_p: cobj,
    b_p: cobj,
    out_p: cobj,
):
    a = Ptr[float](a_p)
    b = Ptr[float](b_p)
    out = Ptr[float](out_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            out[_idx2(i, k, pcols)] = (
                a[_idx2(i, k, pcols)] * b[_idx2(i, k, pcols)] * scale
            )


@export
def diag_phys_writeout_square_codon(
    ncol: int,
    pcols: int,
    pver: int,
    a_p: cobj,
    out_p: cobj,
):
    a = Ptr[float](a_p)
    out = Ptr[float](out_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            val = a[_idx2(i, k, pcols)]
            out[_idx2(i, k, pcols)] = val * val


@export
def diag_phys_writeout_wspeed_codon(
    ncol: int,
    pcols: int,
    pver: int,
    u_p: cobj,
    v_p: cobj,
    out_p: cobj,
):
    u = Ptr[float](u_p)
    v = Ptr[float](v_p)
    out = Ptr[float](out_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            uval = u[_idx2(i, k, pcols)]
            vval = v[_idx2(i, k, pcols)]
            out[_idx2(i, k, pcols)] = sqrt(uval * uval + vval * vval)


@export
def diag_phys_writeout_mass_and_tmq_codon(
    ncol: int,
    pcols: int,
    pver: int,
    rga: float,
    q_p: cobj,
    pdel_p: cobj,
    mq_p: cobj,
    tmq_p: cobj,
):
    q = Ptr[float](q_p)
    pdel = Ptr[float](pdel_p)
    mq = Ptr[float](mq_p)
    tmq = Ptr[float](tmq_p)

    for i in range(1, ncol + 1):
        total = 0.0
        for k in range(1, pver + 1):
            val = q[_idx2(i, k, pcols)] * pdel[_idx2(i, k, pcols)] * rga
            mq[_idx2(i, k, pcols)] = val
            total += val
        tmq[_idx(i)] = total


@export
def diag_phys_writeout_atmeint_codon(
    ncol: int,
    pcols: int,
    pver: int,
    cpair: float,
    latvap: float,
    gravit: float,
    t_p: cobj,
    q_p: cobj,
    u_p: cobj,
    v_p: cobj,
    pdel_p: cobj,
    phis_p: cobj,
    atmeint_p: cobj,
):
    t = Ptr[float](t_p)
    q = Ptr[float](q_p)
    u = Ptr[float](u_p)
    v = Ptr[float](v_p)
    pdel = Ptr[float](pdel_p)
    phis = Ptr[float](phis_p)
    atmeint = Ptr[float](atmeint_p)

    for i in range(1, ncol + 1):
        total = 0.0
        for k in range(1, pver + 1):
            uval = u[_idx2(i, k, pcols)]
            vval = v[_idx2(i, k, pcols)]
            total += (
                cpair * t[_idx2(i, k, pcols)]
                + phis[_idx(i)]
                + latvap * q[_idx2(i, k, pcols)]
                + 0.5 * (uval * uval + vval * vval)
            ) * (pdel[_idx2(i, k, pcols)] / gravit)
        atmeint[_idx(i)] = total


@export
def diag_phys_writeout_wtrc_column_codon(
    ncol: int,
    pcols: int,
    pver: int,
    mode: int,
    rga: float,
    qtr_p: cobj,
    wind_p: cobj,
    pdel_p: cobj,
    out_p: cobj,
):
    qtr = Ptr[float](qtr_p)
    wind = Ptr[float](wind_p)
    pdel = Ptr[float](pdel_p)
    out = Ptr[float](out_p)

    if mode == 1:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, pcols)
                out[idx] = qtr[idx] * pdel[idx] * rga
    else:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, pcols)
                out[idx] = wind[idx] * qtr[idx] * pdel[idx] * rga

    for k in range(2, pver + 1):
        for i in range(1, ncol + 1):
            out[_idx2(i, 1, pcols)] = out[_idx2(i, 1, pcols)] + out[_idx2(i, k, pcols)]


@export
def diag_phys_writeout_ivt_codon(
    ncol: int,
    pcols: int,
    pver: int,
    rga: float,
    q_p: cobj,
    u_p: cobj,
    v_p: cobj,
    pdel_p: cobj,
    uqdp_p: cobj,
    vqdp_p: cobj,
    ivt_p: cobj,
):
    q = Ptr[float](q_p)
    u = Ptr[float](u_p)
    v = Ptr[float](v_p)
    pdel = Ptr[float](pdel_p)
    uqdp = Ptr[float](uqdp_p)
    vqdp = Ptr[float](vqdp_p)
    ivt = Ptr[float](ivt_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            uqdp[idx] = q[idx] * u[idx] * pdel[idx] * rga
            vqdp[idx] = q[idx] * v[idx] * pdel[idx] * rga

    for k in range(2, pver + 1):
        for i in range(1, ncol + 1):
            idx1 = _idx2(i, 1, pcols)
            uqdp[idx1] = uqdp[idx1] + uqdp[_idx2(i, k, pcols)]
            vqdp[idx1] = vqdp[idx1] + vqdp[_idx2(i, k, pcols)]

    for i in range(1, ncol + 1):
        idx1 = _idx2(i, 1, pcols)
        ivt[idx1] = sqrt(uqdp[idx1] ** 2 + vqdp[idx1] ** 2)


@export
def diag_phys_writeout_copy_col1_codon(
    ncol: int,
    pcols: int,
    src_p: cobj,
    dst_p: cobj,
):
    src = Ptr[float](src_p)
    dst = Ptr[float](dst_p)

    for i in range(1, ncol + 1):
        dst[_idx2(i, 1, pcols)] = src[_idx2(i, 1, pcols)]


@export
def diag_phys_writeout_scale_relhum_codon(
    ncol: int,
    pcols: int,
    pver: int,
    q_p: cobj,
    rh_p: cobj,
):
    q = Ptr[float](q_p)
    rh = Ptr[float](rh_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            rh[idx] = q[idx] / rh[idx] * 100.0


@export
def diag_phys_writeout_rhi_rhcfmip_codon(
    ncol: int,
    pcols: int,
    pver: int,
    t_p: cobj,
    esl_p: cobj,
    esi_p: cobj,
    rhw_p: cobj,
    rhi_p: cobj,
    rhcfmip_p: cobj,
):
    t = Ptr[float](t_p)
    esl = Ptr[float](esl_p)
    esi = Ptr[float](esi_p)
    rhw = Ptr[float](rhw_p)
    rhi = Ptr[float](rhi_p)
    rhcfmip = Ptr[float](rhcfmip_p)

    for i in range(1, ncol + 1):
        for k in range(1, pver + 1):
            idx = _idx2(i, k, pcols)
            rhi[idx] = rhw[idx] * esl[idx] / esi[idx]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            rhcfmip[idx] = rhw[idx]

    for i in range(1, ncol + 1):
        for k in range(1, pver + 1):
            idx = _idx2(i, k, pcols)
            if t[idx] > 273.0:
                rhcfmip[idx] = rhw[idx]
            else:
                rhcfmip[idx] = rhi[idx]
