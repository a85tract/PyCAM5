def _idx2(i: int, k: int, ld1: int):
    return (k - 1) * ld1 + (i - 1)


def _idx3(i: int, k: int, m: int, ld1: int, ld2: int):
    return ((m - 1) * ld2 + (k - 1)) * ld1 + (i - 1)


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
