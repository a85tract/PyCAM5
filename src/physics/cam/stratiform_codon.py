@export
def stratiform_implements_cnst_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def stratiform_select_branches_codon(
    use_shfrc: int,
    cam3: int,
    branch_mask_p: cobj,
):
    branch_mask = Ptr[int](branch_mask_p)

    mask = 0
    if use_shfrc != 0:
        mask |= 1
    if cam3 != 0:
        mask |= 2

    branch_mask[0] = mask


@inline
def _idx2(i: int, k: int, ld1: int) -> int:
    return (i - 1) + (k - 1) * ld1


@inline
def _idx3(i: int, k: int, m: int, ld1: int, ld2: int) -> int:
    return (i - 1) + (k - 1) * ld1 + (m - 1) * ld1 * ld2


@export
def stratiform_detrain_assign_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    ixcldliq: int,
    dlf_p: cobj,
    ptend_q_p: cobj,
):
    dlf = Ptr[float](dlf_p)
    ptend_q = Ptr[float](ptend_q_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            ptend_q[_idx3(i, k, ixcldliq, pcols, pver)] = dlf[_idx2(i, k, pcols)]


@export
def stratiform_sedimentation_diag_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gravit: float,
    state_pmid_p: cobj,
    state_t_p: cobj,
    pvliq_p: cobj,
    wsedl_p: cobj,
    rain_p: cobj,
    snow_sed_p: cobj,
    prec_sed_p: cobj,
):
    state_pmid = Ptr[float](state_pmid_p)
    state_t = Ptr[float](state_t_p)
    pvliq = Ptr[float](pvliq_p)
    wsedl = Ptr[float](wsedl_p)
    rain = Ptr[float](rain_p)
    snow_sed = Ptr[float](snow_sed_p)
    prec_sed = Ptr[float](prec_sed_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            wsedl[idx] = pvliq[idx] / gravit / (state_pmid[idx] / (287.15 * state_t[idx]))

    for i in range(1, ncol + 1):
        snow_sed[i - 1] = snow_sed[i - 1] / 1000.0
        rain[i - 1] = rain[i - 1] / 1000.0
        prec_sed[i - 1] = rain[i - 1] + snow_sed[i - 1]


@export
def stratiform_rhdfda_codon(
    ncol: int,
    pcols: int,
    pver: int,
    relhum_p: cobj,
    rhu00_p: cobj,
    cld_p: cobj,
    cld2_p: cobj,
    rhdfda_p: cobj,
):
    relhum = Ptr[float](relhum_p)
    rhu00 = Ptr[float](rhu00_p)
    cld = Ptr[float](cld_p)
    cld2 = Ptr[float](cld2_p)
    rhdfda = Ptr[float](rhdfda_p)

    for i in range(1, ncol + 1):
        rhu00[_idx2(i, 1, pcols)] = 2.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            if relhum[idx] < rhu00[idx]:
                rhdfda[idx] = 0.0
            elif relhum[idx] >= 1.0:
                rhdfda[idx] = 0.0
            else:
                if (cld2[idx] - cld[idx]) < 1.0e-4:
                    rhdfda[idx] = 0.01 * relhum[idx] * 1.0e4
                else:
                    rhdfda[idx] = 0.01 * relhum[idx] / (cld2[idx] - cld[idx])


@export
def stratiform_repartition_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    ixcldliq: int,
    ixcldice: int,
    rdtime: float,
    fice_p: cobj,
    state_q_p: cobj,
    totcw_p: cobj,
    repartht_p: cobj,
    ptend_q_p: cobj,
):
    fice = Ptr[float](fice_p)
    state_q = Ptr[float](state_q_p)
    totcw = Ptr[float](totcw_p)
    repartht = Ptr[float](repartht_p)
    ptend_q = Ptr[float](ptend_q_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            idx_ice = _idx3(i, k, ixcldice, pcols, pver)
            idx_liq = _idx3(i, k, ixcldliq, pcols, pver)
            totcw[idx2] = state_q[idx_ice] + state_q[idx_liq]
            repartht[idx2] = state_q[idx_ice]
            ptend_q[idx_ice] = rdtime * (totcw[idx2] * fice[idx2] - state_q[idx_ice])
            ptend_q[idx_liq] = rdtime * (totcw[idx2] * (1.0 - fice[idx2]) - state_q[idx_liq])


@export
def stratiform_forcing_prep_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    ixcldliq: int,
    ixcldice: int,
    dtime: float,
    latice: float,
    state_q_p: cobj,
    state_t_p: cobj,
    qcwat_p: cobj,
    tcwat_p: cobj,
    lcwat_p: cobj,
    totcw_p: cobj,
    repartht_p: cobj,
    qtend_p: cobj,
    ttend_p: cobj,
    ltend_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    state_t = Ptr[float](state_t_p)
    qcwat = Ptr[float](qcwat_p)
    tcwat = Ptr[float](tcwat_p)
    lcwat = Ptr[float](lcwat_p)
    totcw = Ptr[float](totcw_p)
    repartht = Ptr[float](repartht_p)
    qtend = Ptr[float](qtend_p)
    ttend = Ptr[float](ttend_p)
    ltend = Ptr[float](ltend_p)

    rdtime = 1.0 / dtime

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            idx_qv = _idx3(i, k, 1, pcols, pver)
            idx_ice = _idx3(i, k, ixcldice, pcols, pver)
            repartht[idx2] = (latice / dtime) * (state_q[idx_ice] - repartht[idx2])
            qtend[idx2] = (state_q[idx_qv] - qcwat[idx2]) * rdtime
            ttend[idx2] = (state_t[idx2] - tcwat[idx2]) * rdtime
            ltend[idx2] = (totcw[idx2] - lcwat[idx2]) * rdtime


@export
def stratiform_microphys_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    ixcldliq: int,
    ixcldice: int,
    dtime: float,
    latvap: float,
    latice: float,
    cld_p: cobj,
    concld_p: cobj,
    fice_p: cobj,
    qme_p: cobj,
    nevapr_p: cobj,
    evapheat_p: cobj,
    prfzheat_p: cobj,
    meltheat_p: cobj,
    ice2pr_p: cobj,
    liq2pr_p: cobj,
    state_q_p: cobj,
    prec_pcw_p: cobj,
    snow_pcw_p: cobj,
    repartht_p: cobj,
    ptend_q_p: cobj,
    ptend_s_p: cobj,
    ast_p: cobj,
    icimr_p: cobj,
    icwmr_p: cobj,
    cmeheat_p: cobj,
    cmeice_p: cobj,
    cmeliq_p: cobj,
    mr_ccliq_p: cobj,
    mr_ccice_p: cobj,
    mr_lsliq_p: cobj,
    mr_lsice_p: cobj,
):
    cld = Ptr[float](cld_p)
    concld = Ptr[float](concld_p)
    fice = Ptr[float](fice_p)
    qme = Ptr[float](qme_p)
    nevapr = Ptr[float](nevapr_p)
    evapheat = Ptr[float](evapheat_p)
    prfzheat = Ptr[float](prfzheat_p)
    meltheat = Ptr[float](meltheat_p)
    ice2pr = Ptr[float](ice2pr_p)
    liq2pr = Ptr[float](liq2pr_p)
    state_q = Ptr[float](state_q_p)
    prec_pcw = Ptr[float](prec_pcw_p)
    snow_pcw = Ptr[float](snow_pcw_p)
    repartht = Ptr[float](repartht_p)
    ptend_q = Ptr[float](ptend_q_p)
    ptend_s = Ptr[float](ptend_s_p)
    ast = Ptr[float](ast_p)
    icimr = Ptr[float](icimr_p)
    icwmr = Ptr[float](icwmr_p)
    cmeheat = Ptr[float](cmeheat_p)
    cmeice = Ptr[float](cmeice_p)
    cmeliq = Ptr[float](cmeliq_p)
    mr_ccliq = Ptr[float](mr_ccliq_p)
    mr_ccice = Ptr[float](mr_ccice_p)
    mr_lsliq = Ptr[float](mr_lsliq_p)
    mr_lsice = Ptr[float](mr_lsice_p)

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            idx2 = _idx2(i, k, pcols)
            mr_ccliq[idx2] = 0.0
            mr_ccice[idx2] = 0.0
            mr_lsliq[idx2] = 0.0
            mr_lsice[idx2] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            idx_qv = _idx3(i, k, 1, pcols, pver)
            idx_ice = _idx3(i, k, ixcldice, pcols, pver)
            idx_liq = _idx3(i, k, ixcldliq, pcols, pver)
            ptend_s[idx2] = (
                qme[idx2] * (latvap + latice * fice[idx2])
                + evapheat[idx2]
                + prfzheat[idx2]
                + meltheat[idx2]
                + repartht[idx2]
            )
            ptend_q[idx_qv] = -qme[idx2] + nevapr[idx2]
            ptend_q[idx_ice] = qme[idx2] * fice[idx2] - ice2pr[idx2]
            ptend_q[idx_liq] = qme[idx2] * (1.0 - fice[idx2]) - liq2pr[idx2]
            ast[idx2] = cld[idx2]
            denom = max(0.01, ast[idx2])
            icimr[idx2] = (state_q[idx_ice] + dtime * ptend_q[idx_ice]) / denom
            icwmr[idx2] = (state_q[idx_liq] + dtime * ptend_q[idx_liq]) / denom
            cmeheat[idx2] = qme[idx2] * (latvap + latice * fice[idx2])
            cmeice[idx2] = qme[idx2] * fice[idx2]
            cmeliq[idx2] = qme[idx2] * (1.0 - fice[idx2])

            if cld[idx2] > 0.0:
                mr_ccliq[idx2] = (state_q[idx_liq] / cld[idx2]) * concld[idx2]
                mr_ccice[idx2] = (state_q[idx_ice] / cld[idx2]) * concld[idx2]
                mr_lsliq[idx2] = (state_q[idx_liq] / cld[idx2]) * (cld[idx2] - concld[idx2])
                mr_lsice[idx2] = (state_q[idx_ice] / cld[idx2]) * (cld[idx2] - concld[idx2])

    for i in range(1, ncol + 1):
        snow_pcw[i - 1] = snow_pcw[i - 1] / 1000.0
        prec_pcw[i - 1] = prec_pcw[i - 1] / 1000.0


@export
def stratiform_microphys_tend_diag_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    ixcldliq: int,
    ixcldice: int,
    dtime: float,
    latvap: float,
    latice: float,
    cld_p: cobj,
    fice_p: cobj,
    qme_p: cobj,
    nevapr_p: cobj,
    evapheat_p: cobj,
    prfzheat_p: cobj,
    meltheat_p: cobj,
    ice2pr_p: cobj,
    liq2pr_p: cobj,
    state_q_p: cobj,
    prec_pcw_p: cobj,
    snow_pcw_p: cobj,
    repartht_p: cobj,
    ptend_q_p: cobj,
    ptend_s_p: cobj,
    ast_p: cobj,
    icimr_p: cobj,
    icwmr_p: cobj,
    cmeheat_p: cobj,
    cmeice_p: cobj,
    cmeliq_p: cobj,
):
    cld = Ptr[float](cld_p)
    fice = Ptr[float](fice_p)
    qme = Ptr[float](qme_p)
    nevapr = Ptr[float](nevapr_p)
    evapheat = Ptr[float](evapheat_p)
    prfzheat = Ptr[float](prfzheat_p)
    meltheat = Ptr[float](meltheat_p)
    ice2pr = Ptr[float](ice2pr_p)
    liq2pr = Ptr[float](liq2pr_p)
    state_q = Ptr[float](state_q_p)
    prec_pcw = Ptr[float](prec_pcw_p)
    snow_pcw = Ptr[float](snow_pcw_p)
    repartht = Ptr[float](repartht_p)
    ptend_q = Ptr[float](ptend_q_p)
    ptend_s = Ptr[float](ptend_s_p)
    ast = Ptr[float](ast_p)
    icimr = Ptr[float](icimr_p)
    icwmr = Ptr[float](icwmr_p)
    cmeheat = Ptr[float](cmeheat_p)
    cmeice = Ptr[float](cmeice_p)
    cmeliq = Ptr[float](cmeliq_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            idx_qv = _idx3(i, k, 1, pcols, pver)
            idx_ice = _idx3(i, k, ixcldice, pcols, pver)
            idx_liq = _idx3(i, k, ixcldliq, pcols, pver)
            ptend_s[idx2] = (
                qme[idx2] * (latvap + latice * fice[idx2])
                + evapheat[idx2]
                + prfzheat[idx2]
                + meltheat[idx2]
                + repartht[idx2]
            )
            ptend_q[idx_qv] = -qme[idx2] + nevapr[idx2]
            ptend_q[idx_ice] = qme[idx2] * fice[idx2] - ice2pr[idx2]
            ptend_q[idx_liq] = qme[idx2] * (1.0 - fice[idx2]) - liq2pr[idx2]
            ast[idx2] = cld[idx2]
            denom = max(0.01, ast[idx2])
            icimr[idx2] = (state_q[idx_ice] + dtime * ptend_q[idx_ice]) / denom
            icwmr[idx2] = (state_q[idx_liq] + dtime * ptend_q[idx_liq]) / denom
            cmeheat[idx2] = qme[idx2] * (latvap + latice * fice[idx2])
            cmeice[idx2] = qme[idx2] * fice[idx2]
            cmeliq[idx2] = qme[idx2] * (1.0 - fice[idx2])

    for i in range(1, ncol + 1):
        snow_pcw[i - 1] = snow_pcw[i - 1] / 1000.0
        prec_pcw[i - 1] = prec_pcw[i - 1] / 1000.0


@export
def stratiform_cloud_mixing_diag_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    ixcldliq: int,
    ixcldice: int,
    state_q_p: cobj,
    cld_p: cobj,
    concld_p: cobj,
    mr_ccliq_p: cobj,
    mr_ccice_p: cobj,
    mr_lsliq_p: cobj,
    mr_lsice_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    cld = Ptr[float](cld_p)
    concld = Ptr[float](concld_p)
    mr_ccliq = Ptr[float](mr_ccliq_p)
    mr_ccice = Ptr[float](mr_ccice_p)
    mr_lsliq = Ptr[float](mr_lsliq_p)
    mr_lsice = Ptr[float](mr_lsice_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            idx_ice = _idx3(i, k, ixcldice, pcols, pver)
            idx_liq = _idx3(i, k, ixcldliq, pcols, pver)
            if cld[idx2] > 0.0:
                mr_ccliq[idx2] = (state_q[idx_liq] / cld[idx2]) * concld[idx2]
                mr_ccice[idx2] = (state_q[idx_ice] / cld[idx2]) * concld[idx2]
                mr_lsliq[idx2] = (state_q[idx_liq] / cld[idx2]) * (cld[idx2] - concld[idx2])
                mr_lsice[idx2] = (state_q[idx_ice] / cld[idx2]) * (cld[idx2] - concld[idx2])
            else:
                mr_ccliq[idx2] = 0.0
                mr_ccice[idx2] = 0.0
                mr_lsliq[idx2] = 0.0
                mr_lsice[idx2] = 0.0


@export
def stratiform_postcloud_diag_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    ixcldliq: int,
    ixcldice: int,
    state_q_p: cobj,
    state_pmid_p: cobj,
    state_t_p: cobj,
    rhcloud_p: cobj,
    iwc_p: cobj,
    lwc_p: cobj,
    icimr_p: cobj,
    icwmr_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    state_pmid = Ptr[float](state_pmid_p)
    state_t = Ptr[float](state_t_p)
    rhcloud = Ptr[float](rhcloud_p)
    iwc = Ptr[float](iwc_p)
    lwc = Ptr[float](lwc_p)
    icimr = Ptr[float](icimr_p)
    icwmr = Ptr[float](icwmr_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            idx_ice = _idx3(i, k, ixcldice, pcols, pver)
            idx_liq = _idx3(i, k, ixcldliq, pcols, pver)
            iwc[idx2] = state_q[idx_ice] * state_pmid[idx2] / (287.15 * state_t[idx2])
            lwc[idx2] = state_q[idx_liq] * state_pmid[idx2] / (287.15 * state_t[idx2])
            icimr[idx2] = state_q[idx_ice] / max(0.01, rhcloud[idx2])
            icwmr[idx2] = state_q[idx_liq] / max(0.01, rhcloud[idx2])


@export
def stratiform_store_oldcloud_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    ixcldliq: int,
    ixcldice: int,
    state_q_p: cobj,
    state_t_p: cobj,
    qcwat_p: cobj,
    tcwat_p: cobj,
    lcwat_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    state_t = Ptr[float](state_t_p)
    qcwat = Ptr[float](qcwat_p)
    tcwat = Ptr[float](tcwat_p)
    lcwat = Ptr[float](lcwat_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            qcwat[idx2] = state_q[_idx3(i, k, 1, pcols, pver)]
            tcwat[idx2] = state_t[idx2]
            lcwat[idx2] = (
                state_q[_idx3(i, k, ixcldice, pcols, pver)]
                + state_q[_idx3(i, k, ixcldliq, pcols, pver)]
            )
