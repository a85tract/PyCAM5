from math import gamma, pi


def _idx2(i: int, k: int, ld1: int):
    return (k - 1) * ld1 + (i - 1)


def _idx3(i: int, k: int, m: int, ld1: int, ld2: int):
    return ((m - 1) * ld2 + (k - 1)) * ld1 + (i - 1)


@inline
def _process_rates_idx(
    i: int,
    k: int,
    idsttype: int,
    isrctype: int,
    rtype: int,
    pcols: int,
    pver: int,
    pwtype: int,
):
    return (
        (i - 1)
        + (k - 1) * pcols
        + (idsttype - 1) * pcols * pver
        + (isrctype - 1) * pcols * pver * pwtype
        + (rtype - 1) * pcols * pver * pwtype * pwtype
    )


@inline
def _add_process_rate(
    process_rates,
    i: int,
    k: int,
    isrctype: int,
    idsttype: int,
    rtype: int,
    rate_val: float,
    pcols: int,
    pver: int,
    pwtype: int,
):
    dst_idx = _process_rates_idx(i, k, idsttype, isrctype, rtype, pcols, pver, pwtype)
    process_rates[dst_idx] = process_rates[dst_idx] + rate_val
    if isrctype != idsttype:
        src_idx = _process_rates_idx(i, k, isrctype, idsttype, rtype, pcols, pver, pwtype)
        process_rates[src_idx] = process_rates[src_idx] - rate_val


@export
def micro_mg_cam_premg_diag_codon(
    ncol: int,
    psetcols: int,
    pver: int,
    top_lev: int,
    ixcldliq: int,
    ixcldice: int,
    mincld: float,
    gravit: float,
    state_q_p: cobj,
    state_pdel_p: cobj,
    ast_p: cobj,
    cldo_p: cobj,
    iclwpi_p: cobj,
    iciwpi_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    state_pdel = Ptr[float](state_pdel_p)
    ast = Ptr[float](ast_p)
    cldo = Ptr[float](cldo_p)
    iclwpi = Ptr[float](iclwpi_p)
    iciwpi = Ptr[float](iciwpi_p)

    for i in range(1, psetcols + 1):
        iclwpi[i - 1] = 0.0
        iciwpi[i - 1] = 0.0

    for i in range(1, ncol + 1):
        for k in range(top_lev, pver + 1):
            idx2 = _idx2(i, k, psetcols)
            idx_q_liq = _idx3(i, k, ixcldliq, psetcols, pver)
            idx_q_ice = _idx3(i, k, ixcldice, psetcols, pver)
            iclwpi[i - 1] = iclwpi[i - 1] + min(
                state_q[idx_q_liq] / max(mincld, ast[idx2]), 0.005
            ) * state_pdel[idx2] / gravit
            iciwpi[i - 1] = iciwpi[i - 1] + min(
                state_q[idx_q_ice] / max(mincld, ast[idx2]), 0.005
            ) * state_pdel[idx2] / gravit

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, psetcols)
            cldo[idx2] = ast[idx2]


@export
def micro_mg_cam_wtrc_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    pwtype: int,
    iwtvap: int,
    iwtliq: int,
    iwtice: int,
    iwtstrain: int,
    iwtstsnow: int,
    preo_grid_p: cobj,
    prdso_grid_p: cobj,
    cmeiout_grid_p: cobj,
    meltso_grid_p: cobj,
    qcsedten_grid_p: cobj,
    qisedten_grid_p: cobj,
    mnuccco_grid_p: cobj,
    mnuccto_grid_p: cobj,
    msacwio_grid_p: cobj,
    prao_grid_p: cobj,
    prco_grid_p: cobj,
    psacwso_grid_p: cobj,
    bergo_grid_p: cobj,
    bergso_grid_p: cobj,
    praio_grid_p: cobj,
    prcio_grid_p: cobj,
    pracso_grid_p: cobj,
    mnuccro_grid_p: cobj,
    qcreso_grid_p: cobj,
    qireso_grid_p: cobj,
    homoo_grid_p: cobj,
    melto_grid_p: cobj,
    pre_rates_grid_p: cobj,
    sed_rates_grid_p: cobj,
    post_rates_grid_p: cobj,
    pcmei_grid_p: cobj,
    ncmei_grid_p: cobj,
    pmelts_grid_p: cobj,
    nmelts_grid_p: cobj,
):
    preo_grid = Ptr[float](preo_grid_p)
    prdso_grid = Ptr[float](prdso_grid_p)
    cmeiout_grid = Ptr[float](cmeiout_grid_p)
    meltso_grid = Ptr[float](meltso_grid_p)
    qcsedten_grid = Ptr[float](qcsedten_grid_p)
    qisedten_grid = Ptr[float](qisedten_grid_p)
    mnuccco_grid = Ptr[float](mnuccco_grid_p)
    mnuccto_grid = Ptr[float](mnuccto_grid_p)
    msacwio_grid = Ptr[float](msacwio_grid_p)
    prao_grid = Ptr[float](prao_grid_p)
    prco_grid = Ptr[float](prco_grid_p)
    psacwso_grid = Ptr[float](psacwso_grid_p)
    bergo_grid = Ptr[float](bergo_grid_p)
    bergso_grid = Ptr[float](bergso_grid_p)
    praio_grid = Ptr[float](praio_grid_p)
    prcio_grid = Ptr[float](prcio_grid_p)
    pracso_grid = Ptr[float](pracso_grid_p)
    mnuccro_grid = Ptr[float](mnuccro_grid_p)
    qcreso_grid = Ptr[float](qcreso_grid_p)
    qireso_grid = Ptr[float](qireso_grid_p)
    homoo_grid = Ptr[float](homoo_grid_p)
    melto_grid = Ptr[float](melto_grid_p)
    pre_rates_grid = Ptr[float](pre_rates_grid_p)
    sed_rates_grid = Ptr[float](sed_rates_grid_p)
    post_rates_grid = Ptr[float](post_rates_grid_p)
    pcmei_grid = Ptr[float](pcmei_grid_p)
    ncmei_grid = Ptr[float](ncmei_grid_p)
    pmelts_grid = Ptr[float](pmelts_grid_p)
    nmelts_grid = Ptr[float](nmelts_grid_p)

    for rtype in range(1, pwtype + 1):
        for isrctype in range(1, pwtype + 1):
            for idsttype in range(1, pwtype + 1):
                for k in range(top_lev, pver + 1):
                    for i in range(1, pcols + 1):
                        pre_rates_grid[
                            _process_rates_idx(
                                i, k, idsttype, isrctype, rtype, pcols, pver, pwtype
                            )
                        ] = 0.0
                        post_rates_grid[
                            _process_rates_idx(
                                i, k, idsttype, isrctype, rtype, pcols, pver, pwtype
                            )
                        ] = 0.0

    for m in range(1, pwtype + 1):
        for k in range(top_lev, pver + 1):
            for i in range(1, pcols + 1):
                sed_rates_grid[_idx3(i, k, m, pcols, pver)] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, pcols + 1):
            idx2 = _idx2(i, k, pcols)
            pcmei_grid[idx2] = 0.0
            ncmei_grid[idx2] = 0.0
            pmelts_grid[idx2] = 0.0
            nmelts_grid[idx2] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            if cmeiout_grid[idx2] < 0.0:
                ncmei_grid[idx2] = cmeiout_grid[idx2]
            else:
                pcmei_grid[idx2] = cmeiout_grid[idx2]
            if meltso_grid[idx2] < 0.0:
                nmelts_grid[idx2] = meltso_grid[idx2]
            else:
                pmelts_grid[idx2] = meltso_grid[idx2]
            sed_rates_grid[_idx3(i, k, iwtliq, pcols, pver)] = qcsedten_grid[idx2]
            sed_rates_grid[_idx3(i, k, iwtice, pcols, pver)] = qisedten_grid[idx2]

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            _add_process_rate(
                pre_rates_grid, i, k, iwtvap, iwtice, iwtvap, pcmei_grid[idx2], pcols, pver, pwtype
            )
            _add_process_rate(
                pre_rates_grid, i, k, iwtvap, iwtice, iwtice, ncmei_grid[idx2], pcols, pver, pwtype
            )
            _add_process_rate(
                pre_rates_grid, i, k, iwtvap, iwtstrain, iwtstrain, preo_grid[idx2], pcols, pver, pwtype
            )
            _add_process_rate(
                pre_rates_grid, i, k, iwtvap, iwtstsnow, iwtstsnow, prdso_grid[idx2], pcols, pver, pwtype
            )

            rate_val = mnuccco_grid[idx2] + mnuccto_grid[idx2]
            rate_val = rate_val + msacwio_grid[idx2]
            _add_process_rate(
                pre_rates_grid, i, k, iwtliq, iwtice, iwtliq, rate_val, pcols, pver, pwtype
            )

            rate_val = prao_grid[idx2] + prco_grid[idx2]
            _add_process_rate(
                pre_rates_grid, i, k, iwtliq, iwtstrain, iwtliq, rate_val, pcols, pver, pwtype
            )
            _add_process_rate(
                pre_rates_grid, i, k, iwtliq, iwtstsnow, iwtliq, psacwso_grid[idx2], pcols, pver, pwtype
            )
            _add_process_rate(
                pre_rates_grid, i, k, iwtliq, iwtliq, iwtliq, bergo_grid[idx2], pcols, pver, pwtype
            )
            _add_process_rate(
                pre_rates_grid, i, k, iwtice, iwtice, iwtice, bergso_grid[idx2], pcols, pver, pwtype
            )

            rate_val = praio_grid[idx2] + prcio_grid[idx2]
            _add_process_rate(
                pre_rates_grid, i, k, iwtice, iwtstsnow, iwtice, rate_val, pcols, pver, pwtype
            )

            rate_val = pracso_grid[idx2] + mnuccro_grid[idx2]
            _add_process_rate(
                pre_rates_grid, i, k, iwtstrain, iwtstsnow, iwtstrain, rate_val, pcols, pver, pwtype
            )

            _add_process_rate(
                post_rates_grid, i, k, iwtvap, iwtliq, iwtvap, qcreso_grid[idx2], pcols, pver, pwtype
            )
            _add_process_rate(
                post_rates_grid, i, k, iwtvap, iwtice, iwtvap, qireso_grid[idx2], pcols, pver, pwtype
            )
            _add_process_rate(
                post_rates_grid, i, k, iwtliq, iwtice, iwtliq, homoo_grid[idx2], pcols, pver, pwtype
            )
            _add_process_rate(
                post_rates_grid, i, k, iwtice, iwtliq, iwtice, melto_grid[idx2], pcols, pver, pwtype
            )


def _size_dist_param_basic_codon(
    qsmall: float,
    qic: float,
    nic: float,
    eff_dim: float,
    shape_coef: float,
    lambda_lo: float,
    lambda_hi: float,
    min_mean_mass: float,
):
    lam = 0.0
    nic_out = nic

    if qic > qsmall:
        nic_out = min(nic_out, qic / min_mean_mass)
        lam = (shape_coef * nic_out / qic) ** (1.0 / eff_dim)

        if lam < lambda_lo:
            lam = lambda_lo
            nic_out = lam**eff_dim * qic / shape_coef
        elif lam > lambda_hi:
            lam = lambda_hi
            nic_out = lam**eff_dim * qic / shape_coef

    return lam, nic_out


def _avg_diameter_codon(q: float, n: float, rho_air: float, rho_sub: float):
    return (pi * rho_sub * n / (q * rho_air)) ** (-1.0 / 3.0)


def _size_dist_param_liq_codon(
    qsmall: float,
    qcic: float,
    ncic: float,
    rho_air: float,
    liq_rho: float,
    liq_eff_dim: float,
    liq_min_mean_mass: float,
):
    pgam = -100.0
    lamc = 0.0
    ncic_out = ncic

    if qcic > qsmall:
        pgam = 0.0005714 * (ncic / 1.0e6 * rho_air) + 0.2714
        pgam = 1.0 / (pgam ** 2) - 1.0
        pgam = max(pgam, 2.0)
        pgam = min(pgam, 15.0)

        shape_coef = pi * liq_rho / 6.0 * (gamma(pgam + 1.0 + liq_eff_dim) / gamma(pgam + 1.0))
        lambda_lo = (pgam + 1.0) * 1.0 / 50.0e-6
        lambda_hi = (pgam + 1.0) * 1.0 / 2.0e-6
        lamc, ncic_out = _size_dist_param_basic_codon(
            qsmall,
            qcic,
            ncic_out,
            liq_eff_dim,
            shape_coef,
            lambda_lo,
            lambda_hi,
            liq_min_mean_mass,
        )

    return pgam, lamc, ncic_out


@export
def micro_mg_cam_postmg_diag_codon(
    ncol: int,
    psetcols: int,
    pver: int,
    pverp: int,
    top_lev: int,
    micro_mg_version: int,
    rate1_cw2pr_st_idx: int,
    ixcldliq: int,
    ixcldice: int,
    ixnumliq: int,
    ixnumice: int,
    ixrain: int,
    ixsnow: int,
    cpair: float,
    gravit: float,
    mincld: float,
    qsmall: float,
    state_q_p: cobj,
    state_t_p: cobj,
    state_pmid_p: cobj,
    state_pdel_p: cobj,
    naai_p: cobj,
    naai_hom_p: cobj,
    mnuccdo_p: cobj,
    rflx_p: cobj,
    sflx_p: cobj,
    qrout_p: cobj,
    qsout_p: cobj,
    prect_p: cobj,
    preci_p: cobj,
    rate1cld_p: cobj,
    vtrmc_p: cobj,
    tlat_p: cobj,
    qvlat_p: cobj,
    qcten_p: cobj,
    qiten_p: cobj,
    ncten_p: cobj,
    niten_p: cobj,
    alst_mic_p: cobj,
    cmeliq_p: cobj,
    cmeiout_p: cobj,
    ast_p: cobj,
    cld_p: cobj,
    concld_p: cobj,
    mnuccdohet_p: cobj,
    mgflxprc_p: cobj,
    mgflxsnw_p: cobj,
    mgmrprc_p: cobj,
    mgmrsnw_p: cobj,
    cvreffliq_p: cobj,
    cvreffice_p: cobj,
    rate1ord_cw2pr_st_p: cobj,
    wsedl_p: cobj,
    cc_t_p: cobj,
    cc_qv_p: cobj,
    cc_ql_p: cobj,
    cc_qi_p: cobj,
    cc_nl_p: cobj,
    cc_ni_p: cobj,
    cc_qlst_p: cobj,
    qme_p: cobj,
    prec_pcw_p: cobj,
    snow_pcw_p: cobj,
    prec_sed_p: cobj,
    snow_sed_p: cobj,
    prec_str_p: cobj,
    snow_str_p: cobj,
    icecldf_p: cobj,
    liqcldf_p: cobj,
    icinc_p: cobj,
    icwnc_p: cobj,
    iciwpst_p: cobj,
    iclwpst_p: cobj,
    icswp_p: cobj,
    cldfsnow_p: cobj,
    icimrst_p: cobj,
    icwmrst_p: cobj,
    cldmax_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    state_t = Ptr[float](state_t_p)
    state_pmid = Ptr[float](state_pmid_p)
    state_pdel = Ptr[float](state_pdel_p)
    naai = Ptr[float](naai_p)
    naai_hom = Ptr[float](naai_hom_p)
    mnuccdo = Ptr[float](mnuccdo_p)
    rflx = Ptr[float](rflx_p)
    sflx = Ptr[float](sflx_p)
    qrout = Ptr[float](qrout_p)
    qsout = Ptr[float](qsout_p)
    prect = Ptr[float](prect_p)
    preci = Ptr[float](preci_p)
    rate1cld = Ptr[float](rate1cld_p)
    vtrmc = Ptr[float](vtrmc_p)
    tlat = Ptr[float](tlat_p)
    qvlat = Ptr[float](qvlat_p)
    qcten = Ptr[float](qcten_p)
    qiten = Ptr[float](qiten_p)
    ncten = Ptr[float](ncten_p)
    niten = Ptr[float](niten_p)
    alst_mic = Ptr[float](alst_mic_p)
    cmeliq = Ptr[float](cmeliq_p)
    cmeiout = Ptr[float](cmeiout_p)
    ast = Ptr[float](ast_p)
    cld = Ptr[float](cld_p)
    concld = Ptr[float](concld_p)
    mnuccdohet = Ptr[float](mnuccdohet_p)
    mgflxprc = Ptr[float](mgflxprc_p)
    mgflxsnw = Ptr[float](mgflxsnw_p)
    mgmrprc = Ptr[float](mgmrprc_p)
    mgmrsnw = Ptr[float](mgmrsnw_p)
    cvreffliq = Ptr[float](cvreffliq_p)
    cvreffice = Ptr[float](cvreffice_p)
    rate1ord_cw2pr_st = Ptr[float](rate1ord_cw2pr_st_p)
    wsedl = Ptr[float](wsedl_p)
    cc_t = Ptr[float](cc_t_p)
    cc_qv = Ptr[float](cc_qv_p)
    cc_ql = Ptr[float](cc_ql_p)
    cc_qi = Ptr[float](cc_qi_p)
    cc_nl = Ptr[float](cc_nl_p)
    cc_ni = Ptr[float](cc_ni_p)
    cc_qlst = Ptr[float](cc_qlst_p)
    qme = Ptr[float](qme_p)
    prec_pcw = Ptr[float](prec_pcw_p)
    snow_pcw = Ptr[float](snow_pcw_p)
    prec_sed = Ptr[float](prec_sed_p)
    snow_sed = Ptr[float](snow_sed_p)
    prec_str = Ptr[float](prec_str_p)
    snow_str = Ptr[float](snow_str_p)
    icecldf = Ptr[float](icecldf_p)
    liqcldf = Ptr[float](liqcldf_p)
    icinc = Ptr[float](icinc_p)
    icwnc = Ptr[float](icwnc_p)
    iciwpst = Ptr[float](iciwpst_p)
    iclwpst = Ptr[float](iclwpst_p)
    icswp = Ptr[float](icswp_p)
    cldfsnow = Ptr[float](cldfsnow_p)
    icimrst = Ptr[float](icimrst_p)
    icwmrst = Ptr[float](icwmrst_p)
    cldmax = Ptr[float](cldmax_p)

    for k in range(1, pver + 1):
        for i in range(1, psetcols + 1):
            idx2 = _idx2(i, k, psetcols)
            mnuccdohet[idx2] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, psetcols)
            if naai[idx2] > 0.0:
                mnuccdohet[idx2] = mnuccdo[idx2] - (naai_hom[idx2] / naai[idx2]) * mnuccdo[idx2]

    for k in range(top_lev, pverp + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, psetcols)
            mgflxprc[idx2] = rflx[idx2] + sflx[idx2]
            mgflxsnw[idx2] = sflx[idx2]

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, psetcols)
            mgmrprc[idx2] = qrout[idx2] + qsout[idx2]
            mgmrsnw[idx2] = qsout[idx2]

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, psetcols)
            cvreffliq[idx2] = 9.0
            cvreffice[idx2] = 37.0

    if rate1_cw2pr_st_idx > 0:
        for k in range(top_lev, pver + 1):
            for i in range(1, ncol + 1):
                idx2 = _idx2(i, k, psetcols)
                rate1ord_cw2pr_st[idx2] = rate1cld[idx2]

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, psetcols)
            wsedl[idx2] = vtrmc[idx2]
            cc_t[idx2] = tlat[idx2] / cpair
            cc_qv[idx2] = qvlat[idx2]
            cc_ql[idx2] = qcten[idx2]
            cc_qi[idx2] = qiten[idx2]
            cc_nl[idx2] = ncten[idx2]
            cc_ni[idx2] = niten[idx2]
            cc_qlst[idx2] = qcten[idx2] / max(0.01, alst_mic[idx2])
            qme[idx2] = cmeliq[idx2] + cmeiout[idx2]
            icecldf[idx2] = ast[idx2]
            liqcldf[idx2] = ast[idx2]

    for i in range(1, psetcols + 1):
        idx1 = i - 1
        prec_pcw[idx1] = prect[idx1]
        snow_pcw[idx1] = preci[idx1]
        prec_sed[idx1] = 0.0
        snow_sed[idx1] = 0.0
        prec_str[idx1] = prec_pcw[idx1] + prec_sed[idx1]
        snow_str[idx1] = snow_pcw[idx1] + snow_sed[idx1]

    for k in range(1, pver + 1):
        for i in range(1, psetcols + 1):
            idx2 = _idx2(i, k, psetcols)
            icinc[idx2] = 0.0
            icwnc[idx2] = 0.0
            iciwpst[idx2] = 0.0
            iclwpst[idx2] = 0.0
            icswp[idx2] = 0.0
            cldfsnow[idx2] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, psetcols)
            idx_q_liq = _idx3(i, k, ixcldliq, psetcols, pver)
            idx_q_ice = _idx3(i, k, ixcldice, psetcols, pver)
            idx_n_liq = _idx3(i, k, ixnumliq, psetcols, pver)
            idx_n_ice = _idx3(i, k, ixnumice, psetcols, pver)

            icimrst[idx2] = min(state_q[idx_q_ice] / max(mincld, icecldf[idx2]), 0.005)
            icwmrst[idx2] = min(state_q[idx_q_liq] / max(mincld, liqcldf[idx2]), 0.005)
            icinc[idx2] = state_q[idx_n_ice] / max(mincld, icecldf[idx2]) * state_pmid[idx2] / (287.15 * state_t[idx2])
            icwnc[idx2] = state_q[idx_n_liq] / max(mincld, liqcldf[idx2]) * state_pmid[idx2] / (287.15 * state_t[idx2])
            iciwpst[idx2] = min(state_q[idx_q_ice] / max(mincld, ast[idx2]), 0.005) * state_pdel[idx2] / gravit
            iclwpst[idx2] = min(state_q[idx_q_liq] / max(mincld, ast[idx2]), 0.005) * state_pdel[idx2] / gravit

            cldfsnow[idx2] = cld[idx2]
            if cldfsnow[idx2] > 1.0e-4 and concld[idx2] < 1.0e-4 and state_q[idx_q_liq] < 1.0e-10:
                cldfsnow[idx2] = 0.0
            if cldfsnow[idx2] <= 1.0e-4 and qsout[idx2] > 1.0e-6:
                cldfsnow[idx2] = 0.25

            icswp[idx2] = qsout[idx2] / max(mincld, cldfsnow[idx2]) * state_pdel[idx2] / gravit

    if micro_mg_version > 1:
        for k in range(1, pver + 1):
            for i in range(1, psetcols + 1):
                idx2 = _idx2(i, k, psetcols)
                cldmax[idx2] = max(mincld, ast[idx2])

        for k in range(top_lev + 1, pver + 1):
            for i in range(1, ncol + 1):
                idx2 = _idx2(i, k, psetcols)
                idx2_prev = _idx2(i, k - 1, psetcols)
                idx_rain_prev = _idx3(i, k - 1, ixrain, psetcols, pver)
                idx_snow_prev = _idx3(i, k - 1, ixsnow, psetcols, pver)
                if state_q[idx_rain_prev] >= qsmall or state_q[idx_snow_prev] >= qsmall:
                    cldmax[idx2] = max(cldmax[idx2_prev], cldmax[idx2])


@export
def micro_mg_cam_grid_diag_codon(
    ngrdcol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    minlwp: float,
    gravit: float,
    rhoh2o: float,
    iclwpst_grid_p: cobj,
    cld_grid_p: cobj,
    cmeliq_grid_p: cobj,
    pdel_grid_p: cobj,
    prec_str_grid_p: cobj,
    acgcme_grid_p: cobj,
    acprecl_grid_p: cobj,
    acnum_grid_p: cobj,
    prao_grid_p: cobj,
    prco_grid_p: cobj,
    nc_grid_p: cobj,
    liqcldf_grid_p: cobj,
    icwmrst_grid_p: cobj,
    rel_grid_p: cobj,
    icwnc_grid_p: cobj,
    icecldf_grid_p: cobj,
    icimrst_grid_p: cobj,
    rei_grid_p: cobj,
    icinc_grid_p: cobj,
    nevapr_grid_p: cobj,
    evpsnow_st_grid_p: cobj,
    tgliqwp_grid_p: cobj,
    tgcmeliq_grid_p: cobj,
    pe_grid_p: cobj,
    tpr_grid_p: cobj,
    pefrac_grid_p: cobj,
    vprao_grid_p: cobj,
    vprco_grid_p: cobj,
    racau_grid_p: cobj,
    cnt_grid_p: cobj,
    cdnumc_grid_p: cobj,
    efcout_grid_p: cobj,
    efiout_grid_p: cobj,
    ncout_grid_p: cobj,
    niout_grid_p: cobj,
    freql_grid_p: cobj,
    freqi_grid_p: cobj,
    icwmrst_grid_out_p: cobj,
    icimrst_grid_out_p: cobj,
    fcti_grid_p: cobj,
    fctl_grid_p: cobj,
    ctrel_grid_p: cobj,
    ctrei_grid_p: cobj,
    ctnl_grid_p: cobj,
    ctni_grid_p: cobj,
    evprain_st_grid_p: cobj,
):
    iclwpst_grid = Ptr[float](iclwpst_grid_p)
    cld_grid = Ptr[float](cld_grid_p)
    cmeliq_grid = Ptr[float](cmeliq_grid_p)
    pdel_grid = Ptr[float](pdel_grid_p)
    prec_str_grid = Ptr[float](prec_str_grid_p)
    acgcme_grid = Ptr[float](acgcme_grid_p)
    acprecl_grid = Ptr[float](acprecl_grid_p)
    acnum_grid = Ptr[i32](acnum_grid_p)
    prao_grid = Ptr[float](prao_grid_p)
    prco_grid = Ptr[float](prco_grid_p)
    nc_grid = Ptr[float](nc_grid_p)
    liqcldf_grid = Ptr[float](liqcldf_grid_p)
    icwmrst_grid = Ptr[float](icwmrst_grid_p)
    rel_grid = Ptr[float](rel_grid_p)
    icwnc_grid = Ptr[float](icwnc_grid_p)
    icecldf_grid = Ptr[float](icecldf_grid_p)
    icimrst_grid = Ptr[float](icimrst_grid_p)
    rei_grid = Ptr[float](rei_grid_p)
    icinc_grid = Ptr[float](icinc_grid_p)
    nevapr_grid = Ptr[float](nevapr_grid_p)
    evpsnow_st_grid = Ptr[float](evpsnow_st_grid_p)
    tgliqwp_grid = Ptr[float](tgliqwp_grid_p)
    tgcmeliq_grid = Ptr[float](tgcmeliq_grid_p)
    pe_grid = Ptr[float](pe_grid_p)
    tpr_grid = Ptr[float](tpr_grid_p)
    pefrac_grid = Ptr[float](pefrac_grid_p)
    vprao_grid = Ptr[float](vprao_grid_p)
    vprco_grid = Ptr[float](vprco_grid_p)
    racau_grid = Ptr[float](racau_grid_p)
    cnt_grid = Ptr[i32](cnt_grid_p)
    cdnumc_grid = Ptr[float](cdnumc_grid_p)
    efcout_grid = Ptr[float](efcout_grid_p)
    efiout_grid = Ptr[float](efiout_grid_p)
    ncout_grid = Ptr[float](ncout_grid_p)
    niout_grid = Ptr[float](niout_grid_p)
    freql_grid = Ptr[float](freql_grid_p)
    freqi_grid = Ptr[float](freqi_grid_p)
    icwmrst_grid_out = Ptr[float](icwmrst_grid_out_p)
    icimrst_grid_out = Ptr[float](icimrst_grid_out_p)
    fcti_grid = Ptr[float](fcti_grid_p)
    fctl_grid = Ptr[float](fctl_grid_p)
    ctrel_grid = Ptr[float](ctrel_grid_p)
    ctrei_grid = Ptr[float](ctrei_grid_p)
    ctnl_grid = Ptr[float](ctnl_grid_p)
    ctni_grid = Ptr[float](ctni_grid_p)
    evprain_st_grid = Ptr[float](evprain_st_grid_p)

    for i in range(1, ngrdcol + 1):
        idx1 = i - 1
        tgliqwp_grid[idx1] = 0.0
        tgcmeliq_grid[idx1] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ngrdcol + 1):
            idx1 = i - 1
            idx2 = _idx2(i, k, pcols)
            tgliqwp_grid[idx1] = tgliqwp_grid[idx1] + iclwpst_grid[idx2] * cld_grid[idx2]
            if cmeliq_grid[idx2] > 1.0e-12:
                tgcmeliq_grid[idx1] = tgcmeliq_grid[idx1] + cmeliq_grid[idx2] * (pdel_grid[idx2] / gravit) / rhoh2o

    for i in range(1, ngrdcol + 1):
        idx1 = i - 1
        pe_grid[idx1] = 0.0
        tpr_grid[idx1] = 0.0
        pefrac_grid[idx1] = 0.0

    for i in range(1, ngrdcol + 1):
        idx1 = i - 1
        acgcme_grid[idx1] = acgcme_grid[idx1] + tgcmeliq_grid[idx1]
        acprecl_grid[idx1] = acprecl_grid[idx1] + prec_str_grid[idx1]
        acnum_grid[idx1] = i32(int(acnum_grid[idx1]) + 1)

        if tgliqwp_grid[idx1] < minlwp:
            if acprecl_grid[idx1] > 5.0e-8:
                tpr_grid[idx1] = max(acprecl_grid[idx1] / int(acnum_grid[idx1]), 1.0e-15)
                if acgcme_grid[idx1] > 1.0e-10:
                    pe_grid[idx1] = min(max(acprecl_grid[idx1] / acgcme_grid[idx1], 1.0e-15), 1.0e5)
                    pefrac_grid[idx1] = 1.0

            acprecl_grid[idx1] = 0.0
            acgcme_grid[idx1] = 0.0
            acnum_grid[idx1] = i32(0)

        if int(acnum_grid[idx1]) > 1000:
            acnum_grid[idx1] = i32(0)
            acprecl_grid[idx1] = 0.0
            acgcme_grid[idx1] = 0.0

    for i in range(1, pcols + 1):
        idx1 = i - 1
        vprao_grid[idx1] = 0.0
        cnt_grid[idx1] = i32(0)

    for k in range(top_lev, pver + 1):
        for i in range(1, ngrdcol + 1):
            idx1 = i - 1
            idx2 = _idx2(i, k, pcols)
            vprao_grid[idx1] = vprao_grid[idx1] + prao_grid[idx2]
            if prao_grid[idx2] != 0.0:
                cnt_grid[idx1] = i32(int(cnt_grid[idx1]) + 1)

    for i in range(1, pcols + 1):
        idx1 = i - 1
        if int(cnt_grid[idx1]) > 0:
            vprao_grid[idx1] = vprao_grid[idx1] / int(cnt_grid[idx1])

    for i in range(1, pcols + 1):
        idx1 = i - 1
        vprco_grid[idx1] = 0.0
        cnt_grid[idx1] = i32(0)

    for k in range(top_lev, pver + 1):
        for i in range(1, ngrdcol + 1):
            idx1 = i - 1
            idx2 = _idx2(i, k, pcols)
            vprco_grid[idx1] = vprco_grid[idx1] + prco_grid[idx2]
            if prco_grid[idx2] != 0.0:
                cnt_grid[idx1] = i32(int(cnt_grid[idx1]) + 1)

    for i in range(1, pcols + 1):
        idx1 = i - 1
        if int(cnt_grid[idx1]) > 0:
            vprco_grid[idx1] = vprco_grid[idx1] / int(cnt_grid[idx1])
            racau_grid[idx1] = vprao_grid[idx1] / vprco_grid[idx1]
        else:
            racau_grid[idx1] = 0.0

    for i in range(1, pcols + 1):
        idx1 = i - 1
        racau_grid[idx1] = min(racau_grid[idx1], 1.0e10)

    for i in range(1, ngrdcol + 1):
        idx1 = i - 1
        cdnumc_grid[idx1] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ngrdcol + 1):
            idx1 = i - 1
            idx2 = _idx2(i, k, pcols)
            cdnumc_grid[idx1] = cdnumc_grid[idx1] + nc_grid[idx2] * pdel_grid[idx2] / gravit

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            idx2 = _idx2(i, k, pcols)
            efcout_grid[idx2] = 0.0
            efiout_grid[idx2] = 0.0
            ncout_grid[idx2] = 0.0
            niout_grid[idx2] = 0.0
            freql_grid[idx2] = 0.0
            freqi_grid[idx2] = 0.0
            icwmrst_grid_out[idx2] = 0.0
            icimrst_grid_out[idx2] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ngrdcol + 1):
            idx2 = _idx2(i, k, pcols)
            if liqcldf_grid[idx2] > 0.01 and icwmrst_grid[idx2] > 5.0e-5:
                efcout_grid[idx2] = rel_grid[idx2] * liqcldf_grid[idx2]
                ncout_grid[idx2] = icwnc_grid[idx2] * liqcldf_grid[idx2]
                freql_grid[idx2] = liqcldf_grid[idx2]
                icwmrst_grid_out[idx2] = icwmrst_grid[idx2]
            if icecldf_grid[idx2] > 0.01 and icimrst_grid[idx2] > 1.0e-6:
                efiout_grid[idx2] = rei_grid[idx2] * icecldf_grid[idx2]
                niout_grid[idx2] = icinc_grid[idx2] * icecldf_grid[idx2]
                freqi_grid[idx2] = icecldf_grid[idx2]
                icimrst_grid_out[idx2] = icimrst_grid[idx2]

    for i in range(1, pcols + 1):
        idx1 = i - 1
        fcti_grid[idx1] = 0.0
        fctl_grid[idx1] = 0.0
        ctrel_grid[idx1] = 0.0
        ctrei_grid[idx1] = 0.0
        ctnl_grid[idx1] = 0.0
        ctni_grid[idx1] = 0.0

    for i in range(1, ngrdcol + 1):
        idx1 = i - 1
        for k in range(top_lev, pver + 1):
            idx2 = _idx2(i, k, pcols)
            if liqcldf_grid[idx2] > 0.01 and icwmrst_grid[idx2] > 1.0e-7:
                ctrel_grid[idx1] = rel_grid[idx2] * liqcldf_grid[idx2]
                ctnl_grid[idx1] = icwnc_grid[idx2] * liqcldf_grid[idx2]
                fctl_grid[idx1] = liqcldf_grid[idx2]
                break
            if icecldf_grid[idx2] > 0.01 and icimrst_grid[idx2] > 1.0e-7:
                ctrei_grid[idx1] = rei_grid[idx2] * icecldf_grid[idx2]
                ctni_grid[idx1] = icinc_grid[idx2] * icecldf_grid[idx2]
                fcti_grid[idx1] = icecldf_grid[idx2]
                break

    for k in range(1, pver + 1):
        for i in range(1, ngrdcol + 1):
            idx2 = _idx2(i, k, pcols)
            evprain_st_grid[idx2] = nevapr_grid[idx2] - evpsnow_st_grid[idx2]

    for k in range(top_lev, pver + 1):
        for i in range(1, ngrdcol + 1):
            idx2 = _idx2(i, k, pcols)
            evprain_st_grid[idx2] = max(evprain_st_grid[idx2], 0.0)
            evpsnow_st_grid[idx2] = max(evpsnow_st_grid[idx2], 0.0)


@export
def micro_mg_cam_wtrc_prep_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    pwtype: int,
    iwtliq: int,
    iwtice: int,
    cmeiout_grid_p: cobj,
    meltso_grid_p: cobj,
    qcsedten_grid_p: cobj,
    qisedten_grid_p: cobj,
    pcmei_grid_p: cobj,
    ncmei_grid_p: cobj,
    pmelts_grid_p: cobj,
    nmelts_grid_p: cobj,
    sed_rates_grid_p: cobj,
):
    cmeiout_grid = Ptr[float](cmeiout_grid_p)
    meltso_grid = Ptr[float](meltso_grid_p)
    qcsedten_grid = Ptr[float](qcsedten_grid_p)
    qisedten_grid = Ptr[float](qisedten_grid_p)
    pcmei_grid = Ptr[float](pcmei_grid_p)
    ncmei_grid = Ptr[float](ncmei_grid_p)
    pmelts_grid = Ptr[float](pmelts_grid_p)
    nmelts_grid = Ptr[float](nmelts_grid_p)
    sed_rates_grid = Ptr[float](sed_rates_grid_p)

    for m in range(1, pwtype + 1):
        for k in range(top_lev, pver + 1):
            for i in range(1, pcols + 1):
                sed_rates_grid[_idx3(i, k, m, pcols, pver)] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, pcols + 1):
            idx2 = _idx2(i, k, pcols)
            pcmei_grid[idx2] = 0.0
            ncmei_grid[idx2] = 0.0
            pmelts_grid[idx2] = 0.0
            nmelts_grid[idx2] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            if cmeiout_grid[idx2] < 0.0:
                ncmei_grid[idx2] = cmeiout_grid[idx2]
            else:
                pcmei_grid[idx2] = cmeiout_grid[idx2]

            if meltso_grid[idx2] < 0.0:
                nmelts_grid[idx2] = meltso_grid[idx2]
            else:
                pmelts_grid[idx2] = meltso_grid[idx2]

            sed_rates_grid[_idx3(i, k, iwtliq, pcols, pver)] = qcsedten_grid[idx2]
            sed_rates_grid[_idx3(i, k, iwtice, pcols, pver)] = qisedten_grid[idx2]


@export
def micro_mg_cam_budget_diag_codon(
    mode: int,
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    qcreso_grid_p: cobj,
    melto_grid_p: cobj,
    mnuccco_grid_p: cobj,
    mnuccto_grid_p: cobj,
    bergo_grid_p: cobj,
    homoo_grid_p: cobj,
    msacwio_grid_p: cobj,
    prao_grid_p: cobj,
    prco_grid_p: cobj,
    psacwso_grid_p: cobj,
    bergso_grid_p: cobj,
    cmeiout_grid_p: cobj,
    qireso_grid_p: cobj,
    prcio_grid_p: cobj,
    praio_grid_p: cobj,
    ftem_grid_p: cobj,
):
    qcreso_grid = Ptr[float](qcreso_grid_p)
    melto_grid = Ptr[float](melto_grid_p)
    mnuccco_grid = Ptr[float](mnuccco_grid_p)
    mnuccto_grid = Ptr[float](mnuccto_grid_p)
    bergo_grid = Ptr[float](bergo_grid_p)
    homoo_grid = Ptr[float](homoo_grid_p)
    msacwio_grid = Ptr[float](msacwio_grid_p)
    prao_grid = Ptr[float](prao_grid_p)
    prco_grid = Ptr[float](prco_grid_p)
    psacwso_grid = Ptr[float](psacwso_grid_p)
    bergso_grid = Ptr[float](bergso_grid_p)
    cmeiout_grid = Ptr[float](cmeiout_grid_p)
    qireso_grid = Ptr[float](qireso_grid_p)
    prcio_grid = Ptr[float](prcio_grid_p)
    praio_grid = Ptr[float](praio_grid_p)
    ftem_grid = Ptr[float](ftem_grid_p)

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            ftem_grid[_idx2(i, k, pcols)] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)

            if mode == 1:
                ftem_grid[idx2] = qcreso_grid[idx2]
            elif mode == 2:
                tmp = melto_grid[idx2] - mnuccco_grid[idx2]
                tmp = tmp - mnuccto_grid[idx2]
                tmp = tmp - bergo_grid[idx2]
                tmp = tmp - homoo_grid[idx2]
                tmp = tmp - msacwio_grid[idx2]
                ftem_grid[idx2] = tmp
            elif mode == 3:
                tmp = -prao_grid[idx2]
                tmp = tmp - prco_grid[idx2]
                tmp = tmp - psacwso_grid[idx2]
                tmp = tmp - bergso_grid[idx2]
                ftem_grid[idx2] = tmp
            elif mode == 4:
                ftem_grid[idx2] = cmeiout_grid[idx2] + qireso_grid[idx2]
            elif mode == 5:
                tmp = -melto_grid[idx2] + mnuccco_grid[idx2]
                tmp = tmp + mnuccto_grid[idx2]
                tmp = tmp + bergo_grid[idx2]
                tmp = tmp + homoo_grid[idx2]
                tmp = tmp + msacwio_grid[idx2]
                ftem_grid[idx2] = tmp
            elif mode == 6:
                ftem_grid[idx2] = -prcio_grid[idx2] - praio_grid[idx2]


def _micro_mg_cam_pbuf_copy_fields(
    ncol: int,
    pcols: int,
    pver: int,
    copy_qrain: int,
    copy_qsnow: int,
    copy_nrain: int,
    copy_nsnow: int,
    qrout_grid_p: cobj,
    qsout_grid_p: cobj,
    nrout_grid_p: cobj,
    nsout_grid_p: cobj,
    qrout_grid_ptr_p: cobj,
    qsout_grid_ptr_p: cobj,
    nrout_grid_ptr_p: cobj,
    nsout_grid_ptr_p: cobj,
):
    qrout_grid = Ptr[float](qrout_grid_p)
    qsout_grid = Ptr[float](qsout_grid_p)
    nrout_grid = Ptr[float](nrout_grid_p)
    nsout_grid = Ptr[float](nsout_grid_p)
    qrout_grid_ptr = Ptr[float](qrout_grid_ptr_p)
    qsout_grid_ptr = Ptr[float](qsout_grid_ptr_p)
    nrout_grid_ptr = Ptr[float](nrout_grid_ptr_p)
    nsout_grid_ptr = Ptr[float](nsout_grid_ptr_p)

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            idx2 = _idx2(i, k, pcols)
            if copy_qrain == 1:
                qrout_grid_ptr[idx2] = qrout_grid[idx2]
            if copy_qsnow == 1:
                qsout_grid_ptr[idx2] = qsout_grid[idx2]
            if copy_nrain == 1:
                nrout_grid_ptr[idx2] = nrout_grid[idx2]
            if copy_nsnow == 1:
                nsout_grid_ptr[idx2] = nsout_grid[idx2]


@export
def micro_mg_cam_pbuf_copy_codon(
    ncol: int,
    pcols: int,
    pver: int,
    copy_qrain: int,
    copy_qsnow: int,
    copy_nrain: int,
    copy_nsnow: int,
    qrout_grid_p: cobj,
    qsout_grid_p: cobj,
    nrout_grid_p: cobj,
    nsout_grid_p: cobj,
    qrout_grid_ptr_p: cobj,
    qsout_grid_ptr_p: cobj,
    nrout_grid_ptr_p: cobj,
    nsout_grid_ptr_p: cobj,
):
    _micro_mg_cam_pbuf_copy_fields(
        ncol,
        pcols,
        pver,
        copy_qrain,
        copy_qsnow,
        copy_nrain,
        copy_nsnow,
        qrout_grid_p,
        qsout_grid_p,
        nrout_grid_p,
        nsout_grid_p,
        qrout_grid_ptr_p,
        qsout_grid_ptr_p,
        nrout_grid_ptr_p,
        nsout_grid_ptr_p,
    )


@export
def micro_mg_cam_tail_pbuf_copy_stage_dispatch_codon(
    ncol: int,
    pcols: int,
    pver: int,
    copy_qrain: int,
    copy_qsnow: int,
    copy_nrain: int,
    copy_nsnow: int,
    qrout_grid_p: cobj,
    qsout_grid_p: cobj,
    nrout_grid_p: cobj,
    nsout_grid_p: cobj,
    qrout_grid_ptr_p: cobj,
    qsout_grid_ptr_p: cobj,
    nrout_grid_ptr_p: cobj,
    nsout_grid_ptr_p: cobj,
):
    _micro_mg_cam_pbuf_copy_fields(
        ncol,
        pcols,
        pver,
        copy_qrain,
        copy_qsnow,
        copy_nrain,
        copy_nsnow,
        qrout_grid_p,
        qsout_grid_p,
        nrout_grid_p,
        nsout_grid_p,
        qrout_grid_ptr_p,
        qsout_grid_ptr_p,
        nrout_grid_ptr_p,
        nsout_grid_ptr_p,
    )

@export
def micro_mg_cam_tail_pbuf_copy_codon(
    ncol: int,
    pcols: int,
    pver: int,
    copy_qrain: int,
    copy_qsnow: int,
    copy_nrain: int,
    copy_nsnow: int,
    qrout_grid_p: cobj,
    qsout_grid_p: cobj,
    nrout_grid_p: cobj,
    nsout_grid_p: cobj,
    qrout_grid_ptr_p: cobj,
    qsout_grid_ptr_p: cobj,
    nrout_grid_ptr_p: cobj,
    nsout_grid_ptr_p: cobj,
):
    micro_mg_cam_tail_pbuf_copy_stage_dispatch_codon(
        ncol,
        pcols,
        pver,
        copy_qrain,
        copy_qsnow,
        copy_nrain,
        copy_nsnow,
        qrout_grid_p,
        qsout_grid_p,
        nrout_grid_p,
        nsout_grid_p,
        qrout_grid_ptr_p,
        qsout_grid_ptr_p,
        nrout_grid_ptr_p,
        nsout_grid_ptr_p,
    )


@export
def micro_mg_cam_reff_calc_codon(
    ngrdcol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    micro_mg_version: int,
    qsmall: float,
    mincld: float,
    liq_rho: float,
    liq_eff_dim: float,
    liq_min_mean_mass: float,
    ice_eff_dim: float,
    ice_shape_coef: float,
    ice_lambda_lo: float,
    ice_lambda_hi: float,
    ice_min_mean_mass: float,
    rhosn: float,
    rhoi: float,
    rhow: float,
    rhows: float,
    mucon: float,
    dcon: float,
    deicon: float,
    rho_grid_p: cobj,
    icwmrst_grid_p: cobj,
    liqcldf_grid_p: cobj,
    nc_grid_p: cobj,
    qr_grid_p: cobj,
    nr_grid_p: cobj,
    qs_grid_p: cobj,
    ns_grid_p: cobj,
    qrout_grid_p: cobj,
    nrout_grid_p: cobj,
    qsout_grid_p: cobj,
    nsout_grid_p: cobj,
    ni_grid_p: cobj,
    icecldf_grid_p: cobj,
    icimrst_grid_p: cobj,
    ast_grid_p: cobj,
    mu_grid_p: cobj,
    lambdac_grid_p: cobj,
    rel_fn_grid_p: cobj,
    ncic_grid_p: cobj,
    rel_grid_p: cobj,
    drout2_grid_p: cobj,
    reff_rain_grid_p: cobj,
    des_grid_p: cobj,
    dsout2_grid_p: cobj,
    reff_snow_grid_p: cobj,
    rei_grid_p: cobj,
    niic_grid_p: cobj,
    dei_grid_p: cobj,
    mgreffrain_grid_p: cobj,
    mgreffsnow_grid_p: cobj,
):
    rho_grid = Ptr[float](rho_grid_p)
    icwmrst_grid = Ptr[float](icwmrst_grid_p)
    liqcldf_grid = Ptr[float](liqcldf_grid_p)
    nc_grid = Ptr[float](nc_grid_p)
    qr_grid = Ptr[float](qr_grid_p)
    nr_grid = Ptr[float](nr_grid_p)
    qs_grid = Ptr[float](qs_grid_p)
    ns_grid = Ptr[float](ns_grid_p)
    qrout_grid = Ptr[float](qrout_grid_p)
    nrout_grid = Ptr[float](nrout_grid_p)
    qsout_grid = Ptr[float](qsout_grid_p)
    nsout_grid = Ptr[float](nsout_grid_p)
    ni_grid = Ptr[float](ni_grid_p)
    icecldf_grid = Ptr[float](icecldf_grid_p)
    icimrst_grid = Ptr[float](icimrst_grid_p)
    ast_grid = Ptr[float](ast_grid_p)
    mu_grid = Ptr[float](mu_grid_p)
    lambdac_grid = Ptr[float](lambdac_grid_p)
    rel_fn_grid = Ptr[float](rel_fn_grid_p)
    ncic_grid = Ptr[float](ncic_grid_p)
    rel_grid = Ptr[float](rel_grid_p)
    drout2_grid = Ptr[float](drout2_grid_p)
    reff_rain_grid = Ptr[float](reff_rain_grid_p)
    des_grid = Ptr[float](des_grid_p)
    dsout2_grid = Ptr[float](dsout2_grid_p)
    reff_snow_grid = Ptr[float](reff_snow_grid_p)
    rei_grid = Ptr[float](rei_grid_p)
    niic_grid = Ptr[float](niic_grid_p)
    dei_grid = Ptr[float](dei_grid_p)
    mgreffrain_grid = Ptr[float](mgreffrain_grid_p)
    mgreffsnow_grid = Ptr[float](mgreffsnow_grid_p)

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            idx2 = _idx2(i, k, pcols)
            mu_grid[idx2] = 0.0
            lambdac_grid[idx2] = 0.0
            rel_fn_grid[idx2] = 10.0
            ncic_grid[idx2] = 1.0e8

    for k in range(top_lev, pver + 1):
        for i in range(1, ngrdcol + 1):
            idx2 = _idx2(i, k, pcols)
            mu_grid[idx2], lambdac_grid[idx2], ncic_grid[idx2] = _size_dist_param_liq_codon(
                qsmall,
                icwmrst_grid[idx2],
                ncic_grid[idx2],
                rho_grid[idx2],
                liq_rho,
                liq_eff_dim,
                liq_min_mean_mass,
            )
            if icwmrst_grid[idx2] > qsmall:
                rel_fn_grid[idx2] = (mu_grid[idx2] + 3.0) / lambdac_grid[idx2] / 2.0 * 1.0e6

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            idx2 = _idx2(i, k, pcols)
            mu_grid[idx2] = 0.0
            lambdac_grid[idx2] = 0.0
            rel_grid[idx2] = 10.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ngrdcol + 1):
            idx2 = _idx2(i, k, pcols)
            ncic_grid[idx2] = nc_grid[idx2] / max(mincld, liqcldf_grid[idx2])
            mu_grid[idx2], lambdac_grid[idx2], ncic_grid[idx2] = _size_dist_param_liq_codon(
                qsmall,
                icwmrst_grid[idx2],
                ncic_grid[idx2],
                rho_grid[idx2],
                liq_rho,
                liq_eff_dim,
                liq_min_mean_mass,
            )
            if icwmrst_grid[idx2] >= qsmall:
                rel_grid[idx2] = (mu_grid[idx2] + 3.0) / lambdac_grid[idx2] / 2.0 * 1.0e6
            else:
                mu_grid[idx2] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            idx2 = _idx2(i, k, pcols)
            drout2_grid[idx2] = 0.0
            reff_rain_grid[idx2] = 0.0
            des_grid[idx2] = 0.0
            dsout2_grid[idx2] = 0.0
            reff_snow_grid[idx2] = 0.0

    if micro_mg_version > 1:
        for k in range(top_lev, pver + 1):
            for i in range(1, ngrdcol + 1):
                idx2 = _idx2(i, k, pcols)
                if qr_grid[idx2] >= 1.0e-7:
                    drout2_grid[idx2] = _avg_diameter_codon(
                        qr_grid[idx2],
                        nr_grid[idx2] * rho_grid[idx2],
                        rho_grid[idx2],
                        rhow,
                    )
                    reff_rain_grid[idx2] = drout2_grid[idx2] * 1.5 * 1.0e6

                if qs_grid[idx2] >= 1.0e-7:
                    dsout2_grid[idx2] = _avg_diameter_codon(
                        qs_grid[idx2],
                        ns_grid[idx2] * rho_grid[idx2],
                        rho_grid[idx2],
                        rhosn,
                    )
                    des_grid[idx2] = dsout2_grid[idx2] * 3.0 * rhosn / rhows
                    reff_snow_grid[idx2] = dsout2_grid[idx2] * 1.5 * 1.0e6
    else:
        for k in range(top_lev, pver + 1):
            for i in range(1, ngrdcol + 1):
                idx2 = _idx2(i, k, pcols)
                if qrout_grid[idx2] >= 1.0e-7:
                    drout2_grid[idx2] = _avg_diameter_codon(
                        qrout_grid[idx2],
                        nrout_grid[idx2] * rho_grid[idx2],
                        rho_grid[idx2],
                        rhow,
                    )
                    reff_rain_grid[idx2] = drout2_grid[idx2] * 1.5 * 1.0e6

                if qsout_grid[idx2] >= 1.0e-7:
                    dsout2_grid[idx2] = _avg_diameter_codon(
                        qsout_grid[idx2],
                        nsout_grid[idx2] * rho_grid[idx2],
                        rho_grid[idx2],
                        rhosn,
                    )
                    des_grid[idx2] = dsout2_grid[idx2] * 3.0 * rhosn / rhows
                    reff_snow_grid[idx2] = dsout2_grid[idx2] * 1.5 * 1.0e6

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            idx2 = _idx2(i, k, pcols)
            rei_grid[idx2] = 25.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ngrdcol + 1):
            idx2 = _idx2(i, k, pcols)
            niic_grid[idx2] = ni_grid[idx2] / max(mincld, icecldf_grid[idx2])
            rei_grid[idx2], niic_grid[idx2] = _size_dist_param_basic_codon(
                qsmall,
                icimrst_grid[idx2],
                niic_grid[idx2],
                ice_eff_dim,
                ice_shape_coef,
                ice_lambda_lo,
                ice_lambda_hi,
                ice_min_mean_mass,
            )
            if icimrst_grid[idx2] >= qsmall:
                rei_grid[idx2] = 1.5 / rei_grid[idx2] * 1.0e6
            else:
                rei_grid[idx2] = 25.0

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            idx2 = _idx2(i, k, pcols)
            dei_grid[idx2] = rei_grid[idx2] * rhoi / rhows * 2.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ngrdcol + 1):
            idx2 = _idx2(i, k, pcols)
            des_grid[idx2] = des_grid[idx2] * 1.0e6
            if ast_grid[idx2] < 1.0e-4:
                mu_grid[idx2] = mucon
                lambdac_grid[idx2] = (mucon + 1.0) / dcon
                dei_grid[idx2] = deicon

    for k in range(top_lev, pver + 1):
        for i in range(1, ngrdcol + 1):
            idx2 = _idx2(i, k, pcols)
            mgreffrain_grid[idx2] = reff_rain_grid[idx2]
            mgreffsnow_grid[idx2] = reff_snow_grid[idx2]


@export
def micro_mg_cam_rho_grid_codon(
    ncol: int,
    pcols: int,
    pver: int,
    psetcols: int,
    top_lev: int,
    rair: float,
    rho_p: cobj,
    pmid_p: cobj,
    t_p: cobj,
    rho_grid_p: cobj,
):
    rho = Ptr[float](rho_p)
    pmid = Ptr[float](pmid_p)
    t = Ptr[float](t_p)
    rho_grid = Ptr[float](rho_grid_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            rho[_idx2(i, k, psetcols)] = pmid[_idx2(i, k, psetcols)] / (
                rair * t[_idx2(i, k, psetcols)]
            )

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            rho_grid[_idx2(i, k, pcols)] = rho[_idx2(i, k, psetcols)]


@export
def micro_mg_cam_diag_stage_dispatch_codon(
    ngrdcol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    micro_mg_version: int,
    qsmall: float,
    mincld: float,
    liq_rho: float,
    liq_eff_dim: float,
    liq_min_mean_mass: float,
    ice_eff_dim: float,
    ice_shape_coef: float,
    ice_lambda_lo: float,
    ice_lambda_hi: float,
    ice_min_mean_mass: float,
    rhosn: float,
    rhoi: float,
    rhow: float,
    rhows: float,
    mucon: float,
    dcon: float,
    deicon: float,
    minlwp: float,
    gravit: float,
    rhoh2o: float,
    rho_grid_p: cobj,
    icwmrst_grid_p: cobj,
    liqcldf_grid_p: cobj,
    nc_grid_p: cobj,
    qr_grid_p: cobj,
    nr_grid_p: cobj,
    qs_grid_p: cobj,
    ns_grid_p: cobj,
    qrout_grid_p: cobj,
    nrout_grid_p: cobj,
    qsout_grid_p: cobj,
    nsout_grid_p: cobj,
    ni_grid_p: cobj,
    icecldf_grid_p: cobj,
    icimrst_grid_p: cobj,
    ast_grid_p: cobj,
    mu_grid_p: cobj,
    lambdac_grid_p: cobj,
    rel_fn_grid_p: cobj,
    ncic_grid_p: cobj,
    rel_grid_p: cobj,
    drout2_grid_p: cobj,
    reff_rain_grid_p: cobj,
    des_grid_p: cobj,
    dsout2_grid_p: cobj,
    reff_snow_grid_p: cobj,
    rei_grid_p: cobj,
    niic_grid_p: cobj,
    dei_grid_p: cobj,
    mgreffrain_grid_p: cobj,
    mgreffsnow_grid_p: cobj,
    iclwpst_grid_p: cobj,
    cld_grid_p: cobj,
    cmeliq_grid_p: cobj,
    pdel_grid_p: cobj,
    prec_str_grid_p: cobj,
    acgcme_grid_p: cobj,
    acprecl_grid_p: cobj,
    acnum_grid_p: cobj,
    prao_grid_p: cobj,
    prco_grid_p: cobj,
    icwnc_grid_p: cobj,
    icinc_grid_p: cobj,
    nevapr_grid_p: cobj,
    evpsnow_st_grid_p: cobj,
    tgliqwp_grid_p: cobj,
    tgcmeliq_grid_p: cobj,
    pe_grid_p: cobj,
    tpr_grid_p: cobj,
    pefrac_grid_p: cobj,
    vprao_grid_p: cobj,
    vprco_grid_p: cobj,
    racau_grid_p: cobj,
    cnt_grid_p: cobj,
    cdnumc_grid_p: cobj,
    efcout_grid_p: cobj,
    efiout_grid_p: cobj,
    ncout_grid_p: cobj,
    niout_grid_p: cobj,
    freql_grid_p: cobj,
    freqi_grid_p: cobj,
    icwmrst_grid_out_p: cobj,
    icimrst_grid_out_p: cobj,
    fcti_grid_p: cobj,
    fctl_grid_p: cobj,
    ctrel_grid_p: cobj,
    ctrei_grid_p: cobj,
    ctnl_grid_p: cobj,
    ctni_grid_p: cobj,
    evprain_st_grid_p: cobj,
    qcreso_grid_p: cobj,
    melto_grid_p: cobj,
    mnuccco_grid_p: cobj,
    mnuccto_grid_p: cobj,
    bergo_grid_p: cobj,
    homoo_grid_p: cobj,
    msacwio_grid_p: cobj,
    psacwso_grid_p: cobj,
    bergso_grid_p: cobj,
    cmeiout_grid_p: cobj,
    qireso_grid_p: cobj,
    prcio_grid_p: cobj,
    praio_grid_p: cobj,
    budget_ftem_grid_p: cobj,
    copy_qrain: int,
    copy_qsnow: int,
    copy_nrain: int,
    copy_nsnow: int,
    qrout_grid_ptr_p: cobj,
    qsout_grid_ptr_p: cobj,
    nrout_grid_ptr_p: cobj,
    nsout_grid_ptr_p: cobj,
):
    micro_mg_cam_reff_calc_codon(
        ngrdcol,
        pcols,
        pver,
        top_lev,
        micro_mg_version,
        qsmall,
        mincld,
        liq_rho,
        liq_eff_dim,
        liq_min_mean_mass,
        ice_eff_dim,
        ice_shape_coef,
        ice_lambda_lo,
        ice_lambda_hi,
        ice_min_mean_mass,
        rhosn,
        rhoi,
        rhow,
        rhows,
        mucon,
        dcon,
        deicon,
        rho_grid_p,
        icwmrst_grid_p,
        liqcldf_grid_p,
        nc_grid_p,
        qr_grid_p,
        nr_grid_p,
        qs_grid_p,
        ns_grid_p,
        qrout_grid_p,
        nrout_grid_p,
        qsout_grid_p,
        nsout_grid_p,
        ni_grid_p,
        icecldf_grid_p,
        icimrst_grid_p,
        ast_grid_p,
        mu_grid_p,
        lambdac_grid_p,
        rel_fn_grid_p,
        ncic_grid_p,
        rel_grid_p,
        drout2_grid_p,
        reff_rain_grid_p,
        des_grid_p,
        dsout2_grid_p,
        reff_snow_grid_p,
        rei_grid_p,
        niic_grid_p,
        dei_grid_p,
        mgreffrain_grid_p,
        mgreffsnow_grid_p,
    )

    micro_mg_cam_grid_diag_codon(
        ngrdcol,
        pcols,
        pver,
        top_lev,
        minlwp,
        gravit,
        rhoh2o,
        iclwpst_grid_p,
        cld_grid_p,
        cmeliq_grid_p,
        pdel_grid_p,
        prec_str_grid_p,
        acgcme_grid_p,
        acprecl_grid_p,
        acnum_grid_p,
        prao_grid_p,
        prco_grid_p,
        nc_grid_p,
        liqcldf_grid_p,
        icwmrst_grid_p,
        rel_grid_p,
        icwnc_grid_p,
        icecldf_grid_p,
        icimrst_grid_p,
        rei_grid_p,
        icinc_grid_p,
        nevapr_grid_p,
        evpsnow_st_grid_p,
        tgliqwp_grid_p,
        tgcmeliq_grid_p,
        pe_grid_p,
        tpr_grid_p,
        pefrac_grid_p,
        vprao_grid_p,
        vprco_grid_p,
        racau_grid_p,
        cnt_grid_p,
        cdnumc_grid_p,
        efcout_grid_p,
        efiout_grid_p,
        ncout_grid_p,
        niout_grid_p,
        freql_grid_p,
        freqi_grid_p,
        icwmrst_grid_out_p,
        icimrst_grid_out_p,
        fcti_grid_p,
        fctl_grid_p,
        ctrel_grid_p,
        ctrei_grid_p,
        ctnl_grid_p,
        ctni_grid_p,
        evprain_st_grid_p,
    )

    budget_ftem_grid = Ptr[float](budget_ftem_grid_p)
    qcreso_grid = Ptr[float](qcreso_grid_p)
    melto_grid = Ptr[float](melto_grid_p)
    mnuccco_grid = Ptr[float](mnuccco_grid_p)
    mnuccto_grid = Ptr[float](mnuccto_grid_p)
    bergo_grid = Ptr[float](bergo_grid_p)
    homoo_grid = Ptr[float](homoo_grid_p)
    msacwio_grid = Ptr[float](msacwio_grid_p)
    prao_grid = Ptr[float](prao_grid_p)
    prco_grid = Ptr[float](prco_grid_p)
    psacwso_grid = Ptr[float](psacwso_grid_p)
    bergso_grid = Ptr[float](bergso_grid_p)
    cmeiout_grid = Ptr[float](cmeiout_grid_p)
    qireso_grid = Ptr[float](qireso_grid_p)
    prcio_grid = Ptr[float](prcio_grid_p)
    praio_grid = Ptr[float](praio_grid_p)

    for mode in range(1, 7):
        for k in range(1, pver + 1):
            for i in range(1, pcols + 1):
                budget_ftem_grid[_idx3(i, k, mode, pcols, pver)] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ngrdcol + 1):
            idx2 = _idx2(i, k, pcols)

            budget_ftem_grid[_idx3(i, k, 1, pcols, pver)] = qcreso_grid[idx2]

            tmp = melto_grid[idx2] - mnuccco_grid[idx2]
            tmp = tmp - mnuccto_grid[idx2]
            tmp = tmp - bergo_grid[idx2]
            tmp = tmp - homoo_grid[idx2]
            tmp = tmp - msacwio_grid[idx2]
            budget_ftem_grid[_idx3(i, k, 2, pcols, pver)] = tmp

            tmp = -prao_grid[idx2]
            tmp = tmp - prco_grid[idx2]
            tmp = tmp - psacwso_grid[idx2]
            tmp = tmp - bergso_grid[idx2]
            budget_ftem_grid[_idx3(i, k, 3, pcols, pver)] = tmp

            budget_ftem_grid[_idx3(i, k, 4, pcols, pver)] = (
                cmeiout_grid[idx2] + qireso_grid[idx2]
            )

            tmp = -melto_grid[idx2] + mnuccco_grid[idx2]
            tmp = tmp + mnuccto_grid[idx2]
            tmp = tmp + bergo_grid[idx2]
            tmp = tmp + homoo_grid[idx2]
            tmp = tmp + msacwio_grid[idx2]
            budget_ftem_grid[_idx3(i, k, 5, pcols, pver)] = tmp

            budget_ftem_grid[_idx3(i, k, 6, pcols, pver)] = (
                -prcio_grid[idx2] - praio_grid[idx2]
            )

    if copy_qrain == 1 or copy_qsnow == 1 or copy_nrain == 1 or copy_nsnow == 1:
        _micro_mg_cam_pbuf_copy_fields(
            ngrdcol,
            pcols,
            pver,
            copy_qrain,
            copy_qsnow,
            copy_nrain,
            copy_nsnow,
            qrout_grid_p,
            qsout_grid_p,
            nrout_grid_p,
            nsout_grid_p,
            qrout_grid_ptr_p,
            qsout_grid_ptr_p,
            nrout_grid_ptr_p,
            nsout_grid_ptr_p,
        )

@export
def micro_mg_cam_diag_shell_codon(
    ngrdcol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    micro_mg_version: int,
    qsmall: float,
    mincld: float,
    liq_rho: float,
    liq_eff_dim: float,
    liq_min_mean_mass: float,
    ice_eff_dim: float,
    ice_shape_coef: float,
    ice_lambda_lo: float,
    ice_lambda_hi: float,
    ice_min_mean_mass: float,
    rhosn: float,
    rhoi: float,
    rhow: float,
    rhows: float,
    mucon: float,
    dcon: float,
    deicon: float,
    minlwp: float,
    gravit: float,
    rhoh2o: float,
    rho_grid_p: cobj,
    icwmrst_grid_p: cobj,
    liqcldf_grid_p: cobj,
    nc_grid_p: cobj,
    qr_grid_p: cobj,
    nr_grid_p: cobj,
    qs_grid_p: cobj,
    ns_grid_p: cobj,
    qrout_grid_p: cobj,
    nrout_grid_p: cobj,
    qsout_grid_p: cobj,
    nsout_grid_p: cobj,
    ni_grid_p: cobj,
    icecldf_grid_p: cobj,
    icimrst_grid_p: cobj,
    ast_grid_p: cobj,
    mu_grid_p: cobj,
    lambdac_grid_p: cobj,
    rel_fn_grid_p: cobj,
    ncic_grid_p: cobj,
    rel_grid_p: cobj,
    drout2_grid_p: cobj,
    reff_rain_grid_p: cobj,
    des_grid_p: cobj,
    dsout2_grid_p: cobj,
    reff_snow_grid_p: cobj,
    rei_grid_p: cobj,
    niic_grid_p: cobj,
    dei_grid_p: cobj,
    mgreffrain_grid_p: cobj,
    mgreffsnow_grid_p: cobj,
    iclwpst_grid_p: cobj,
    cld_grid_p: cobj,
    cmeliq_grid_p: cobj,
    pdel_grid_p: cobj,
    prec_str_grid_p: cobj,
    acgcme_grid_p: cobj,
    acprecl_grid_p: cobj,
    acnum_grid_p: cobj,
    prao_grid_p: cobj,
    prco_grid_p: cobj,
    icwnc_grid_p: cobj,
    icinc_grid_p: cobj,
    nevapr_grid_p: cobj,
    evpsnow_st_grid_p: cobj,
    tgliqwp_grid_p: cobj,
    tgcmeliq_grid_p: cobj,
    pe_grid_p: cobj,
    tpr_grid_p: cobj,
    pefrac_grid_p: cobj,
    vprao_grid_p: cobj,
    vprco_grid_p: cobj,
    racau_grid_p: cobj,
    cnt_grid_p: cobj,
    cdnumc_grid_p: cobj,
    efcout_grid_p: cobj,
    efiout_grid_p: cobj,
    ncout_grid_p: cobj,
    niout_grid_p: cobj,
    freql_grid_p: cobj,
    freqi_grid_p: cobj,
    icwmrst_grid_out_p: cobj,
    icimrst_grid_out_p: cobj,
    fcti_grid_p: cobj,
    fctl_grid_p: cobj,
    ctrel_grid_p: cobj,
    ctrei_grid_p: cobj,
    ctnl_grid_p: cobj,
    ctni_grid_p: cobj,
    evprain_st_grid_p: cobj,
    qcreso_grid_p: cobj,
    melto_grid_p: cobj,
    mnuccco_grid_p: cobj,
    mnuccto_grid_p: cobj,
    bergo_grid_p: cobj,
    homoo_grid_p: cobj,
    msacwio_grid_p: cobj,
    psacwso_grid_p: cobj,
    bergso_grid_p: cobj,
    cmeiout_grid_p: cobj,
    qireso_grid_p: cobj,
    prcio_grid_p: cobj,
    praio_grid_p: cobj,
    budget_ftem_grid_p: cobj,
    copy_qrain: int,
    copy_qsnow: int,
    copy_nrain: int,
    copy_nsnow: int,
    qrout_grid_ptr_p: cobj,
    qsout_grid_ptr_p: cobj,
    nrout_grid_ptr_p: cobj,
    nsout_grid_ptr_p: cobj,
):
    micro_mg_cam_diag_stage_dispatch_codon(
        ngrdcol,
        pcols,
        pver,
        top_lev,
        micro_mg_version,
        qsmall,
        mincld,
        liq_rho,
        liq_eff_dim,
        liq_min_mean_mass,
        ice_eff_dim,
        ice_shape_coef,
        ice_lambda_lo,
        ice_lambda_hi,
        ice_min_mean_mass,
        rhosn,
        rhoi,
        rhow,
        rhows,
        mucon,
        dcon,
        deicon,
        minlwp,
        gravit,
        rhoh2o,
        rho_grid_p,
        icwmrst_grid_p,
        liqcldf_grid_p,
        nc_grid_p,
        qr_grid_p,
        nr_grid_p,
        qs_grid_p,
        ns_grid_p,
        qrout_grid_p,
        nrout_grid_p,
        qsout_grid_p,
        nsout_grid_p,
        ni_grid_p,
        icecldf_grid_p,
        icimrst_grid_p,
        ast_grid_p,
        mu_grid_p,
        lambdac_grid_p,
        rel_fn_grid_p,
        ncic_grid_p,
        rel_grid_p,
        drout2_grid_p,
        reff_rain_grid_p,
        des_grid_p,
        dsout2_grid_p,
        reff_snow_grid_p,
        rei_grid_p,
        niic_grid_p,
        dei_grid_p,
        mgreffrain_grid_p,
        mgreffsnow_grid_p,
        iclwpst_grid_p,
        cld_grid_p,
        cmeliq_grid_p,
        pdel_grid_p,
        prec_str_grid_p,
        acgcme_grid_p,
        acprecl_grid_p,
        acnum_grid_p,
        prao_grid_p,
        prco_grid_p,
        icwnc_grid_p,
        icinc_grid_p,
        nevapr_grid_p,
        evpsnow_st_grid_p,
        tgliqwp_grid_p,
        tgcmeliq_grid_p,
        pe_grid_p,
        tpr_grid_p,
        pefrac_grid_p,
        vprao_grid_p,
        vprco_grid_p,
        racau_grid_p,
        cnt_grid_p,
        cdnumc_grid_p,
        efcout_grid_p,
        efiout_grid_p,
        ncout_grid_p,
        niout_grid_p,
        freql_grid_p,
        freqi_grid_p,
        icwmrst_grid_out_p,
        icimrst_grid_out_p,
        fcti_grid_p,
        fctl_grid_p,
        ctrel_grid_p,
        ctrei_grid_p,
        ctnl_grid_p,
        ctni_grid_p,
        evprain_st_grid_p,
        qcreso_grid_p,
        melto_grid_p,
        mnuccco_grid_p,
        mnuccto_grid_p,
        bergo_grid_p,
        homoo_grid_p,
        msacwio_grid_p,
        psacwso_grid_p,
        bergso_grid_p,
        cmeiout_grid_p,
        qireso_grid_p,
        prcio_grid_p,
        praio_grid_p,
        budget_ftem_grid_p,
        copy_qrain,
        copy_qsnow,
        copy_nrain,
        copy_nsnow,
        qrout_grid_ptr_p,
        qsout_grid_ptr_p,
        nrout_grid_ptr_p,
        nsout_grid_ptr_p,
    )


@export
def micro_mg_cam_stage_dispatch_codon(
    stage: int,
    ncol: int,
    psetcols: int,
    pcols: int,
    pver: int,
    pverp: int,
    top_lev: int,
    micro_mg_version: int,
    rate1_cw2pr_st_idx: int,
    ixcldliq: int,
    ixcldice: int,
    ixnumliq: int,
    ixnumice: int,
    ixrain: int,
    ixsnow: int,
    pwtype: int,
    iwtvap: int,
    iwtliq: int,
    iwtice: int,
    iwtstrain: int,
    iwtstsnow: int,
    mincld: float,
    gravit: float,
    cpair: float,
    qsmall: float,
    rair: float,
    p1: cobj,
    p2: cobj,
    p3: cobj,
    p4: cobj,
    p5: cobj,
    p6: cobj,
    p7: cobj,
    p8: cobj,
    p9: cobj,
    p10: cobj,
    p11: cobj,
    p12: cobj,
    p13: cobj,
    p14: cobj,
    p15: cobj,
    p16: cobj,
    p17: cobj,
    p18: cobj,
    p19: cobj,
    p20: cobj,
    p21: cobj,
    p22: cobj,
    p23: cobj,
    p24: cobj,
    p25: cobj,
    p26: cobj,
    p27: cobj,
    p28: cobj,
    p29: cobj,
    p30: cobj,
    p31: cobj,
    p32: cobj,
    p33: cobj,
    p34: cobj,
    p35: cobj,
    p36: cobj,
    p37: cobj,
    p38: cobj,
    p39: cobj,
    p40: cobj,
    p41: cobj,
    p42: cobj,
    p43: cobj,
    p44: cobj,
    p45: cobj,
    p46: cobj,
    p47: cobj,
    p48: cobj,
    p49: cobj,
    p50: cobj,
    p51: cobj,
    p52: cobj,
    p53: cobj,
    p54: cobj,
    p55: cobj,
    p56: cobj,
    p57: cobj,
    p58: cobj,
    p59: cobj,
    p60: cobj,
    p61: cobj,
):
    if stage == 1:
        micro_mg_cam_premg_diag_codon(
            ncol,
            psetcols,
            pver,
            top_lev,
            ixcldliq,
            ixcldice,
            mincld,
            gravit,
            p1,
            p2,
            p3,
            p4,
            p5,
            p6,
        )
    elif stage == 2:
        micro_mg_cam_postmg_diag_codon(
            ncol,
            psetcols,
            pver,
            pverp,
            top_lev,
            micro_mg_version,
            rate1_cw2pr_st_idx,
            ixcldliq,
            ixcldice,
            ixnumliq,
            ixnumice,
            ixrain,
            ixsnow,
            cpair,
            gravit,
            mincld,
            qsmall,
            p1,
            p2,
            p3,
            p4,
            p5,
            p6,
            p7,
            p8,
            p9,
            p10,
            p11,
            p12,
            p13,
            p14,
            p15,
            p16,
            p17,
            p18,
            p19,
            p20,
            p21,
            p22,
            p23,
            p24,
            p25,
            p26,
            p27,
            p28,
            p29,
            p30,
            p31,
            p32,
            p33,
            p34,
            p35,
            p36,
            p37,
            p38,
            p39,
            p40,
            p41,
            p42,
            p43,
            p44,
            p45,
            p46,
            p47,
            p48,
            p49,
            p50,
            p51,
            p52,
            p53,
            p54,
            p55,
            p56,
            p57,
            p58,
            p59,
            p60,
            p61,
        )
    elif stage == 3:
        micro_mg_cam_wtrc_shell_codon(
            ncol,
            pcols,
            pver,
            top_lev,
            pwtype,
            iwtvap,
            iwtliq,
            iwtice,
            iwtstrain,
            iwtstsnow,
            p1,
            p2,
            p3,
            p4,
            p5,
            p6,
            p7,
            p8,
            p9,
            p10,
            p11,
            p12,
            p13,
            p14,
            p15,
            p16,
            p17,
            p18,
            p19,
            p20,
            p21,
            p22,
            p23,
            p24,
            p25,
            p26,
            p27,
            p28,
            p29,
        )
    elif stage == 4:
        micro_mg_cam_rho_grid_codon(
            ncol,
            pcols,
            pver,
            psetcols,
            top_lev,
            rair,
            p1,
            p2,
            p3,
            p4,
        )
