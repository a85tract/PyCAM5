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
def diag_conv_tend_ini_copy_batch_codon(
    mode: int,
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    m: int,
    src_p: cobj,
    dst_p: cobj,
):
    src = Ptr[float](src_p)
    dst = Ptr[float](dst_p)

    if mode == 1 or mode == 3:
        # mode 1: state%s -> dtcond; mode 3: state%t -> t_ttend.
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                dst[_idx2(i, k, pcols)] = src[_idx2(i, k, pcols)]
    elif mode == 2:
        # state%q(pcols,pver,pcnst) tracer m -> dqcond_work(pcols,pver).
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                dst[_idx2(i, k, pcols)] = src[_idx3(i, k, m, pcols, pver)]


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
def diag_conv_update_batch_codon(
    mode: int,
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    m: int,
    scalar1: float,
    scalar2: float,
    a_p: cobj,
    b_p: cobj,
    c_p: cobj,
    d_p: cobj,
    e_p: cobj,
    f_p: cobj,
    g_p: cobj,
    h_p: cobj,
    out1_p: cobj,
    out2_p: cobj,
    out3_p: cobj,
    out4_p: cobj,
    out5_p: cobj,
):
    a = Ptr[float](a_p)
    b = Ptr[float](b_p)
    c = Ptr[float](c_p)
    d = Ptr[float](d_p)
    e = Ptr[float](e_p)
    f = Ptr[float](f_p)
    g = Ptr[float](g_p)
    h = Ptr[float](h_p)
    out1 = Ptr[float](out1_p)
    out2 = Ptr[float](out2_p)
    out3 = Ptr[float](out3_p)
    out4 = Ptr[float](out4_p)
    out5 = Ptr[float](out5_p)

    if mode == 1:
        # precip: a=prec_dp, b=snow_dp, c=prec_sh, d=snow_sh,
        # e=prec_sed, f=snow_sed, g=prec_pcw, h=snow_pcw.
        for i in range(1, ncol + 1):
            idx = _idx(i)
            out1[idx] = a[idx] + c[idx]
            out2[idx] = e[idx] + g[idx]
            out3[idx] = b[idx] + d[idx]
            out4[idx] = f[idx] + h[idx]
            out5[idx] = out1[idx] + out2[idx]
    elif mode == 2:
        # water tracer precip total: a/b/c/d are the four precip components.
        for i in range(1, ncol + 1):
            idx = _idx(i)
            out1[idx] = a[idx]
            out1[idx] = out1[idx] + b[idx]
            out1[idx] = out1[idx] + c[idx]
            out1[idx] = out1[idx] + d[idx]
    elif mode == 3:
        # dtcond: a=state%s, out1=dtcond, scalar1=rtdt, scalar2=cpair.
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, pcols)
                out1[idx] = (a[idx] - out1[idx]) * scalar1 / scalar2
    elif mode == 4:
        # dqcond: a=state%q(pcols,pver,pcnst), out1=dqcond_work, scalar1=rtdt.
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, pcols)
                out1[idx] = (a[_idx3(i, k, m, pcols, pver)] - out1[idx]) * scalar1


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
def hub2atm_alloc_init_codon(
    pcols: int,
    pcnst: int,
    n_drydep: int,
    n_megan: int,
    posinf: float,
    init_bucket: int,
    has_ram1: int,
    has_fv: int,
    has_soilw: int,
    has_dstflx: int,
    has_meganflx: int,
    has_depvel: int,
    asdir_p: cobj,
    asdif_p: cobj,
    aldir_p: cobj,
    aldif_p: cobj,
    lwup_p: cobj,
    lhf_p: cobj,
    shf_p: cobj,
    wsx_p: cobj,
    wsy_p: cobj,
    tref_p: cobj,
    qref_p: cobj,
    u10_p: cobj,
    ts_p: cobj,
    sst_p: cobj,
    snowhland_p: cobj,
    snowhice_p: cobj,
    fco2_lnd_p: cobj,
    fco2_ocn_p: cobj,
    fdms_p: cobj,
    landfrac_p: cobj,
    icefrac_p: cobj,
    ocnfrac_p: cobj,
    ram1_p: cobj,
    fv_p: cobj,
    soilw_p: cobj,
    cflx_p: cobj,
    ustar_p: cobj,
    re_p: cobj,
    ssq_p: cobj,
    dstflx_p: cobj,
    meganflx_p: cobj,
    depvel_p: cobj,
    buckH_p: cobj,
    buck16_p: cobj,
    buckD_p: cobj,
    buck18_p: cobj,
):
    asdir = Ptr[float](asdir_p)
    asdif = Ptr[float](asdif_p)
    aldir = Ptr[float](aldir_p)
    aldif = Ptr[float](aldif_p)
    lwup = Ptr[float](lwup_p)
    lhf = Ptr[float](lhf_p)
    shf = Ptr[float](shf_p)
    wsx = Ptr[float](wsx_p)
    wsy = Ptr[float](wsy_p)
    tref = Ptr[float](tref_p)
    qref = Ptr[float](qref_p)
    u10 = Ptr[float](u10_p)
    ts = Ptr[float](ts_p)
    sst = Ptr[float](sst_p)
    snowhland = Ptr[float](snowhland_p)
    snowhice = Ptr[float](snowhice_p)
    fco2_lnd = Ptr[float](fco2_lnd_p)
    fco2_ocn = Ptr[float](fco2_ocn_p)
    fdms = Ptr[float](fdms_p)
    landfrac = Ptr[float](landfrac_p)
    icefrac = Ptr[float](icefrac_p)
    ocnfrac = Ptr[float](ocnfrac_p)
    ram1 = Ptr[float](ram1_p)
    fv = Ptr[float](fv_p)
    soilw = Ptr[float](soilw_p)
    cflx = Ptr[float](cflx_p)
    ustar = Ptr[float](ustar_p)
    re = Ptr[float](re_p)
    ssq = Ptr[float](ssq_p)
    dstflx = Ptr[float](dstflx_p)
    meganflx = Ptr[float](meganflx_p)
    depvel = Ptr[float](depvel_p)
    buckH = Ptr[float](buckH_p)
    buck16 = Ptr[float](buck16_p)
    buckD = Ptr[float](buckD_p)
    buck18 = Ptr[float](buck18_p)

    for i in range(1, pcols + 1):
        idx = _idx(i)
        asdir[idx] = 0.0
        asdif[idx] = 0.0
        aldir[idx] = 0.0
        aldif[idx] = 0.0
        lwup[idx] = 0.0
        lhf[idx] = 0.0
        shf[idx] = 0.0
        wsx[idx] = 0.0
        wsy[idx] = 0.0
        tref[idx] = 0.0
        qref[idx] = 0.0
        u10[idx] = 0.0
        ts[idx] = 0.0
        sst[idx] = 0.0
        snowhland[idx] = 0.0
        snowhice[idx] = 0.0
        fco2_lnd[idx] = 0.0
        fco2_ocn[idx] = 0.0
        fdms[idx] = 0.0
        landfrac[idx] = posinf
        icefrac[idx] = posinf
        ocnfrac[idx] = posinf
        ustar[idx] = 0.0
        re[idx] = 0.0
        ssq[idx] = 0.0
        if has_ram1 != 0:
            ram1[idx] = 0.1
        if has_fv != 0:
            fv[idx] = 0.1
        if has_soilw != 0:
            soilw[idx] = 0.0
        if init_bucket != 0:
            buckH[idx] = 0.021
            buck16[idx] = 0.021
            buckD[idx] = 0.021
            buck18[idx] = 0.021

    for m in range(1, pcnst + 1):
        for i in range(1, pcols + 1):
            cflx[_idx2(i, m, pcols)] = 0.0

    if has_dstflx != 0:
        for m in range(1, 4 + 1):
            for i in range(1, pcols + 1):
                dstflx[_idx2(i, m, pcols)] = 0.0

    if has_meganflx != 0:
        for m in range(1, n_megan + 1):
            for i in range(1, pcols + 1):
                meganflx[_idx2(i, m, pcols)] = 0.0

    if has_depvel != 0:
        for m in range(1, n_drydep + 1):
            for i in range(1, pcols + 1):
                depvel[_idx2(i, m, pcols)] = 0.0


@export
def atm2hub_alloc_init_codon(
    pcols: int,
    pcnst: int,
    tbot_p: cobj,
    zbot_p: cobj,
    ubot_p: cobj,
    vbot_p: cobj,
    qbot_p: cobj,
    pbot_p: cobj,
    rho_p: cobj,
    netsw_p: cobj,
    flwds_p: cobj,
    precsc_p: cobj,
    precsl_p: cobj,
    precc_p: cobj,
    precl_p: cobj,
    soll_p: cobj,
    sols_p: cobj,
    solld_p: cobj,
    solsd_p: cobj,
    thbot_p: cobj,
    co2prog_p: cobj,
    co2diag_p: cobj,
    psl_p: cobj,
    bcphidry_p: cobj,
    bcphodry_p: cobj,
    bcphiwet_p: cobj,
    ocphidry_p: cobj,
    ocphodry_p: cobj,
    ocphiwet_p: cobj,
    dstdry1_p: cobj,
    dstwet1_p: cobj,
    dstdry2_p: cobj,
    dstwet2_p: cobj,
    dstdry3_p: cobj,
    dstwet3_p: cobj,
    dstdry4_p: cobj,
    dstwet4_p: cobj,
    precrl_16O_p: cobj,
    precrl_HDO_p: cobj,
    precrl_18O_p: cobj,
    precsl_16O_p: cobj,
    precsl_HDO_p: cobj,
    precsl_18O_p: cobj,
    precrc_16O_p: cobj,
    precrc_HDO_p: cobj,
    precrc_18O_p: cobj,
    precsc_16O_p: cobj,
    precsc_HDO_p: cobj,
    precsc_18O_p: cobj,
):
    tbot = Ptr[float](tbot_p)
    zbot = Ptr[float](zbot_p)
    ubot = Ptr[float](ubot_p)
    vbot = Ptr[float](vbot_p)
    qbot = Ptr[float](qbot_p)
    pbot = Ptr[float](pbot_p)
    rho = Ptr[float](rho_p)
    netsw = Ptr[float](netsw_p)
    flwds = Ptr[float](flwds_p)
    precsc = Ptr[float](precsc_p)
    precsl = Ptr[float](precsl_p)
    precc = Ptr[float](precc_p)
    precl = Ptr[float](precl_p)
    soll = Ptr[float](soll_p)
    sols = Ptr[float](sols_p)
    solld = Ptr[float](solld_p)
    solsd = Ptr[float](solsd_p)
    thbot = Ptr[float](thbot_p)
    co2prog = Ptr[float](co2prog_p)
    co2diag = Ptr[float](co2diag_p)
    psl = Ptr[float](psl_p)
    bcphidry = Ptr[float](bcphidry_p)
    bcphodry = Ptr[float](bcphodry_p)
    bcphiwet = Ptr[float](bcphiwet_p)
    ocphidry = Ptr[float](ocphidry_p)
    ocphodry = Ptr[float](ocphodry_p)
    ocphiwet = Ptr[float](ocphiwet_p)
    dstdry1 = Ptr[float](dstdry1_p)
    dstwet1 = Ptr[float](dstwet1_p)
    dstdry2 = Ptr[float](dstdry2_p)
    dstwet2 = Ptr[float](dstwet2_p)
    dstdry3 = Ptr[float](dstdry3_p)
    dstwet3 = Ptr[float](dstwet3_p)
    dstdry4 = Ptr[float](dstdry4_p)
    dstwet4 = Ptr[float](dstwet4_p)
    precrl_16O = Ptr[float](precrl_16O_p)
    precrl_HDO = Ptr[float](precrl_HDO_p)
    precrl_18O = Ptr[float](precrl_18O_p)
    precsl_16O = Ptr[float](precsl_16O_p)
    precsl_HDO = Ptr[float](precsl_HDO_p)
    precsl_18O = Ptr[float](precsl_18O_p)
    precrc_16O = Ptr[float](precrc_16O_p)
    precrc_HDO = Ptr[float](precrc_HDO_p)
    precrc_18O = Ptr[float](precrc_18O_p)
    precsc_16O = Ptr[float](precsc_16O_p)
    precsc_HDO = Ptr[float](precsc_HDO_p)
    precsc_18O = Ptr[float](precsc_18O_p)

    for i in range(1, pcols + 1):
        idx = _idx(i)
        tbot[idx] = 0.0
        zbot[idx] = 0.0
        ubot[idx] = 0.0
        vbot[idx] = 0.0
        pbot[idx] = 0.0
        rho[idx] = 0.0
        netsw[idx] = 0.0
        flwds[idx] = 0.0
        precsc[idx] = 0.0
        precsl[idx] = 0.0
        precc[idx] = 0.0
        precl[idx] = 0.0
        soll[idx] = 0.0
        sols[idx] = 0.0
        solld[idx] = 0.0
        solsd[idx] = 0.0
        thbot[idx] = 0.0
        co2prog[idx] = 0.0
        co2diag[idx] = 0.0
        psl[idx] = 0.0
        bcphidry[idx] = 0.0
        bcphodry[idx] = 0.0
        bcphiwet[idx] = 0.0
        ocphidry[idx] = 0.0
        ocphodry[idx] = 0.0
        ocphiwet[idx] = 0.0
        dstdry1[idx] = 0.0
        dstwet1[idx] = 0.0
        dstdry2[idx] = 0.0
        dstwet2[idx] = 0.0
        dstdry3[idx] = 0.0
        dstwet3[idx] = 0.0
        dstdry4[idx] = 0.0
        dstwet4[idx] = 0.0
        precrl_16O[idx] = 0.0
        precrl_HDO[idx] = 0.0
        precrl_18O[idx] = 0.0
        precsl_16O[idx] = 0.0
        precsl_HDO[idx] = 0.0
        precsl_18O[idx] = 0.0
        precrc_16O[idx] = 0.0
        precrc_HDO[idx] = 0.0
        precrc_18O[idx] = 0.0
        precsc_16O[idx] = 0.0
        precsc_HDO[idx] = 0.0
        precsc_18O[idx] = 0.0

    for m in range(1, pcnst + 1):
        for i in range(1, pcols + 1):
            qbot[_idx2(i, m, pcols)] = 0.0


@export
def cam_export_core_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    rair: float,
    mwdry: float,
    mwco2: float,
    co2diag_val: float,
    co2_transport: int,
    co2_idx: int,
    trace_water: int,
    exist16: int,
    existD: int,
    exist18: int,
    state_t_p: cobj,
    state_exner_p: cobj,
    state_zm_p: cobj,
    state_u_p: cobj,
    state_v_p: cobj,
    state_pmid_p: cobj,
    state_q_p: cobj,
    state_ps_p: cobj,
    state_rpdel_p: cobj,
    psm1_p: cobj,
    srfrpdel_p: cobj,
    co2diag_p: cobj,
    co2prog_p: cobj,
    prcsnw_p: cobj,
    prec_dp_p: cobj,
    snow_dp_p: cobj,
    prec_sh_p: cobj,
    snow_sh_p: cobj,
    prec_sed_p: cobj,
    snow_sed_p: cobj,
    prec_pcw_p: cobj,
    snow_pcw_p: cobj,
    tbot_p: cobj,
    thbot_p: cobj,
    zbot_p: cobj,
    ubot_p: cobj,
    vbot_p: cobj,
    pbot_p: cobj,
    rho_p: cobj,
    qbot_p: cobj,
    precc_p: cobj,
    precl_p: cobj,
    precsc_p: cobj,
    precsl_p: cobj,
    precrl_16O_in_p: cobj,
    precsl_16O_in_p: cobj,
    precrc_16O_in_p: cobj,
    precsc_16O_in_p: cobj,
    precrl_HDO_in_p: cobj,
    precsl_HDO_in_p: cobj,
    precrc_HDO_in_p: cobj,
    precsc_HDO_in_p: cobj,
    precrl_18O_in_p: cobj,
    precsl_18O_in_p: cobj,
    precrc_18O_in_p: cobj,
    precsc_18O_in_p: cobj,
    precrl_16O_out_p: cobj,
    precsl_16O_out_p: cobj,
    precrc_16O_out_p: cobj,
    precsc_16O_out_p: cobj,
    precrl_HDO_out_p: cobj,
    precsl_HDO_out_p: cobj,
    precrc_HDO_out_p: cobj,
    precsc_HDO_out_p: cobj,
    precrl_18O_out_p: cobj,
    precsl_18O_out_p: cobj,
    precrc_18O_out_p: cobj,
    precsc_18O_out_p: cobj,
):
    state_t = Ptr[float](state_t_p)
    state_exner = Ptr[float](state_exner_p)
    state_zm = Ptr[float](state_zm_p)
    state_u = Ptr[float](state_u_p)
    state_v = Ptr[float](state_v_p)
    state_pmid = Ptr[float](state_pmid_p)
    state_q = Ptr[float](state_q_p)
    state_ps = Ptr[float](state_ps_p)
    state_rpdel = Ptr[float](state_rpdel_p)
    psm1 = Ptr[float](psm1_p)
    srfrpdel = Ptr[float](srfrpdel_p)
    co2diag = Ptr[float](co2diag_p)
    co2prog = Ptr[float](co2prog_p)
    prcsnw = Ptr[float](prcsnw_p)
    prec_dp = Ptr[float](prec_dp_p)
    snow_dp = Ptr[float](snow_dp_p)
    prec_sh = Ptr[float](prec_sh_p)
    snow_sh = Ptr[float](snow_sh_p)
    prec_sed = Ptr[float](prec_sed_p)
    snow_sed = Ptr[float](snow_sed_p)
    prec_pcw = Ptr[float](prec_pcw_p)
    snow_pcw = Ptr[float](snow_pcw_p)
    tbot = Ptr[float](tbot_p)
    thbot = Ptr[float](thbot_p)
    zbot = Ptr[float](zbot_p)
    ubot = Ptr[float](ubot_p)
    vbot = Ptr[float](vbot_p)
    pbot = Ptr[float](pbot_p)
    rho = Ptr[float](rho_p)
    qbot = Ptr[float](qbot_p)
    precc = Ptr[float](precc_p)
    precl = Ptr[float](precl_p)
    precsc = Ptr[float](precsc_p)
    precsl = Ptr[float](precsl_p)
    precrl_16O_in = Ptr[float](precrl_16O_in_p)
    precsl_16O_in = Ptr[float](precsl_16O_in_p)
    precrc_16O_in = Ptr[float](precrc_16O_in_p)
    precsc_16O_in = Ptr[float](precsc_16O_in_p)
    precrl_HDO_in = Ptr[float](precrl_HDO_in_p)
    precsl_HDO_in = Ptr[float](precsl_HDO_in_p)
    precrc_HDO_in = Ptr[float](precrc_HDO_in_p)
    precsc_HDO_in = Ptr[float](precsc_HDO_in_p)
    precrl_18O_in = Ptr[float](precrl_18O_in_p)
    precsl_18O_in = Ptr[float](precsl_18O_in_p)
    precrc_18O_in = Ptr[float](precrc_18O_in_p)
    precsc_18O_in = Ptr[float](precsc_18O_in_p)
    precrl_16O_out = Ptr[float](precrl_16O_out_p)
    precsl_16O_out = Ptr[float](precsl_16O_out_p)
    precrc_16O_out = Ptr[float](precrc_16O_out_p)
    precsc_16O_out = Ptr[float](precsc_16O_out_p)
    precrl_HDO_out = Ptr[float](precrl_HDO_out_p)
    precsl_HDO_out = Ptr[float](precsl_HDO_out_p)
    precrc_HDO_out = Ptr[float](precrc_HDO_out_p)
    precsc_HDO_out = Ptr[float](precsc_HDO_out_p)
    precrl_18O_out = Ptr[float](precrl_18O_out_p)
    precsl_18O_out = Ptr[float](precsl_18O_out_p)
    precrc_18O_out = Ptr[float](precrc_18O_out_p)
    precsc_18O_out = Ptr[float](precsc_18O_out_p)

    for i in range(1, ncol + 1):
        src = _idx2(i, pver, pcols)
        dst = _idx(i)
        tbot[dst] = state_t[src]
        thbot[dst] = state_t[src] * state_exner[src]
        zbot[dst] = state_zm[src]
        ubot[dst] = state_u[src]
        vbot[dst] = state_v[src]
        pbot[dst] = state_pmid[src]
        rho[dst] = pbot[dst] / (rair * tbot[dst])
        psm1[dst] = state_ps[dst]
        srfrpdel[dst] = state_rpdel[src]

    for m in range(1, pcnst + 1):
        for i in range(1, ncol + 1):
            qbot[_idx2(i, m, pcols)] = state_q[_idx3(i, pver, m, pcols, pver)]

    for i in range(1, ncol + 1):
        idx = _idx(i)
        precc[idx] = prec_dp[idx] + prec_sh[idx]
        precl[idx] = prec_sed[idx] + prec_pcw[idx]
        precsc[idx] = snow_dp[idx] + snow_sh[idx]
        precsl[idx] = snow_sed[idx] + snow_pcw[idx]
        if precc[idx] < 0.0:
            precc[idx] = 0.0
        if precl[idx] < 0.0:
            precl[idx] = 0.0
        if precsc[idx] < 0.0:
            precsc[idx] = 0.0
        if precsl[idx] < 0.0:
            precsl[idx] = 0.0
        if precsc[idx] > precc[idx]:
            precsc[idx] = precc[idx]
        if precsl[idx] > precl[idx]:
            precsl[idx] = precl[idx]
        co2diag[idx] = co2diag_val
        if co2_transport != 0:
            co2prog[idx] = state_q[_idx3(i, pver, co2_idx, pcols, pver)] * 1.0e6 * mwdry / mwco2
        if trace_water != 0:
            if exist16 != 0:
                precrl_16O_out[idx] = precrl_16O_in[idx]
                precsl_16O_out[idx] = precsl_16O_in[idx]
                precrc_16O_out[idx] = precrc_16O_in[idx]
                precsc_16O_out[idx] = precsc_16O_in[idx]
            if existD != 0:
                precrl_HDO_out[idx] = precrl_HDO_in[idx]
                precsl_HDO_out[idx] = precsl_HDO_in[idx]
                precrc_HDO_out[idx] = precrc_HDO_in[idx]
                precsc_HDO_out[idx] = precsc_HDO_in[idx]
            if exist18 != 0:
                precrl_18O_out[idx] = precrl_18O_in[idx]
                precsl_18O_out[idx] = precsl_18O_in[idx]
                precrc_18O_out[idx] = precrc_18O_in[idx]
                precsc_18O_out[idx] = precsc_18O_in[idx]
            if precrl_16O_out[idx] < 0.0:
                precrl_16O_out[idx] = 0.0
            if precrl_HDO_out[idx] < 0.0:
                precrl_HDO_out[idx] = 0.0
            if precrl_18O_out[idx] < 0.0:
                precrl_18O_out[idx] = 0.0
            if precsl_16O_out[idx] < 0.0:
                precsl_16O_out[idx] = 0.0
            if precsl_HDO_out[idx] < 0.0:
                precsl_HDO_out[idx] = 0.0
            if precsl_18O_out[idx] < 0.0:
                precsl_18O_out[idx] = 0.0
            if precrc_16O_out[idx] < 0.0:
                precrc_16O_out[idx] = 0.0
            if precrc_HDO_out[idx] < 0.0:
                precrc_HDO_out[idx] = 0.0
            if precrc_18O_out[idx] < 0.0:
                precrc_18O_out[idx] = 0.0
            if precsc_16O_out[idx] < 0.0:
                precsc_16O_out[idx] = 0.0
            if precsc_HDO_out[idx] < 0.0:
                precsc_HDO_out[idx] = 0.0
            if precsc_18O_out[idx] < 0.0:
                precsc_18O_out[idx] = 0.0
        prcsnw[idx] = precsc[idx] + precsl[idx]


@export
def diag_phys_tend_update_codon(
    mode: int,
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    m: int,
    scalar1: float,
    scalar2: float,
    a_p: cobj,
    b_p: cobj,
    out_p: cobj,
):
    a = Ptr[float](a_p)
    b = Ptr[float](b_p)
    out = Ptr[float](out_p)

    if mode == 1:
        # tmp_t = (tmp_t - state%t) / ztodt; a=state%t, out=tmp_t.
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, pcols)
                out[idx] = (out[idx] - a[idx]) / scalar1
    elif mode == 2:
        # ftem2(:ncol) = heat_glob / cpair.
        val = scalar1 / scalar2
        for i in range(1, ncol + 1):
            out[_idx(i)] = val
    elif mode == 3:
        # ftem3 = tend%dtdt or tend%dtdt - heat_glob/cpair when m==1.
        offset = 0.0
        if m == 1:
            offset = scalar1 / scalar2
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, pcols)
                out[idx] = a[idx] - offset
    elif mode == 4:
        # dry-mass q tendency: out=(state%q(:,:,m)-out)*rtdt.
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, pcols)
                out[idx] = (a[_idx3(i, k, m, pcols, pver)] - out[idx]) * scalar1
    elif mode == 5:
        # physics q tendency: out=(state%q(:,:,m)-initial_field)*rtdt.
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, pcols)
                out[idx] = (a[_idx3(i, k, m, pcols, pver)] - b[idx]) * scalar1
    elif mode == 6:
        # total temperature tendency: out=(state%t-t_ttend)/ztodt.
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, pcols)
                out[idx] = (a[idx] - b[idx]) / scalar1
    elif mode == 7:
        # copy 2D field, used for t_ttend update.
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, pcols)
                out[idx] = a[idx]


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
    elif mode == 6:
        for i in range(1, ncol + 1):
            out[_idx2(i, 1, pcols)] = a[_idx2(i, 1, pcols)]
    elif mode == 7:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, pcols)
                out[idx] = a[idx] / out[idx] * scale


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
def diag_phys_writeout_column_reduce_codon(
    mode: int,
    ncol: int,
    pcols: int,
    pver: int,
    scalar1: float,
    scalar2: float,
    scalar3: float,
    a_p: cobj,
    b_p: cobj,
    c_p: cobj,
    d_p: cobj,
    e_p: cobj,
    f_p: cobj,
    out2d_p: cobj,
    out1d_p: cobj,
):
    a = Ptr[float](a_p)
    b = Ptr[float](b_p)
    c = Ptr[float](c_p)
    d = Ptr[float](d_p)
    e = Ptr[float](e_p)
    f = Ptr[float](f_p)
    out2d = Ptr[float](out2d_p)
    out1d = Ptr[float](out1d_p)

    if mode == 1:
        # mass_and_tmq: a=q, b=pdel, scalar1=rga.
        for i in range(1, ncol + 1):
            total = 0.0
            for k in range(1, pver + 1):
                val = a[_idx2(i, k, pcols)] * b[_idx2(i, k, pcols)] * scalar1
                out2d[_idx2(i, k, pcols)] = val
                total += val
            out1d[_idx(i)] = total
    elif mode == 2:
        # atmeint: a=t, b=q, c=u, d=v, e=pdel, f=phis.
        cpair = scalar1
        latvap = scalar2
        gravit = scalar3
        for i in range(1, ncol + 1):
            total = 0.0
            for k in range(1, pver + 1):
                uval = c[_idx2(i, k, pcols)]
                vval = d[_idx2(i, k, pcols)]
                total += (
                    cpair * a[_idx2(i, k, pcols)]
                    + f[_idx(i)]
                    + latvap * b[_idx2(i, k, pcols)]
                    + 0.5 * (uval * uval + vval * vval)
                ) * (e[_idx2(i, k, pcols)] / gravit)
            out1d[_idx(i)] = total


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
def diag_phys_writeout_transport_moisture_codon(
    mode: int,
    submode: int,
    ncol: int,
    pcols: int,
    pver: int,
    scalar: float,
    a_p: cobj,
    b_p: cobj,
    c_p: cobj,
    d_p: cobj,
    out1_p: cobj,
    out2_p: cobj,
    out3_p: cobj,
):
    a = Ptr[float](a_p)
    b = Ptr[float](b_p)
    c = Ptr[float](c_p)
    d = Ptr[float](d_p)
    out1 = Ptr[float](out1_p)
    out2 = Ptr[float](out2_p)
    out3 = Ptr[float](out3_p)

    if mode == 1:
        # water tracer column: a=qtr, b=wind, c=pdel, out1=out.
        if submode == 1:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    idx = _idx2(i, k, pcols)
                    out1[idx] = a[idx] * c[idx] * scalar
        else:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    idx = _idx2(i, k, pcols)
                    out1[idx] = b[idx] * a[idx] * c[idx] * scalar

        for k in range(2, pver + 1):
            for i in range(1, ncol + 1):
                out1[_idx2(i, 1, pcols)] = out1[_idx2(i, 1, pcols)] + out1[
                    _idx2(i, k, pcols)
                ]
    elif mode == 2:
        # IVT: a=q, b=u, c=v, d=pdel, out1=uqdp, out2=vqdp, out3=ivt.
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, pcols)
                out1[idx] = a[idx] * b[idx] * d[idx] * scalar
                out2[idx] = a[idx] * c[idx] * d[idx] * scalar

        for k in range(2, pver + 1):
            for i in range(1, ncol + 1):
                idx1 = _idx2(i, 1, pcols)
                out1[idx1] = out1[idx1] + out1[_idx2(i, k, pcols)]
                out2[idx1] = out2[idx1] + out2[_idx2(i, k, pcols)]

        for i in range(1, ncol + 1):
            idx1 = _idx2(i, 1, pcols)
            out3[idx1] = sqrt(out1[idx1] ** 2 + out2[idx1] ** 2)
    elif mode == 3:
        # RHI/RHCFMIP: a=t, b=esl, c=esi, d=rhw, out1=rhi, out2=rhcfmip.
        for i in range(1, ncol + 1):
            for k in range(1, pver + 1):
                idx = _idx2(i, k, pcols)
                out1[idx] = d[idx] * b[idx] / c[idx]

        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, pcols)
                out2[idx] = d[idx]

        for i in range(1, ncol + 1):
            for k in range(1, pver + 1):
                idx = _idx2(i, k, pcols)
                if a[idx] > 273.0:
                    out2[idx] = d[idx]
                else:
                    out2[idx] = out1[idx]


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
