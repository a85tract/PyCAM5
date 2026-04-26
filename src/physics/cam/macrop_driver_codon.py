def _idx2(i: int, k: int, pcols: int):
    return (k - 1) * pcols + (i - 1)


@export
def macrop_driver_select_branches_codon(
    micro_do_icesupersat: int,
    trace_water: int,
    wtrc_detrain_in_macrop: int,
    cu_det_st: int,
    use_shfrc: int,
    do_cldice: int,
    do_cldliq: int,
    do_detrain: int,
    branch_mask_p: cobj,
):
    branch_mask = Ptr[int](branch_mask_p)

    mask = 0
    if micro_do_icesupersat != 0:
        mask |= 1
    if trace_water != 0:
        mask |= 2
    if wtrc_detrain_in_macrop != 0:
        mask |= 4
    if cu_det_st != 0:
        mask |= 8
    if use_shfrc != 0:
        mask |= 16
    if do_cldice != 0:
        mask |= 32
    if do_cldliq != 0:
        mask |= 64
    if do_detrain != 0:
        mask |= 128

    branch_mask[0] = mask

@export
def macrop_driver_clr_old_diag_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    concld_p: cobj,
    alst_p: cobj,
    ast_p: cobj,
    clrw_old_p: cobj,
    clri_old_p: cobj,
):
    concld = Ptr[float](concld_p)
    alst = Ptr[float](alst_p)
    ast = Ptr[float](ast_p)
    clrw_old = Ptr[float](clrw_old_p)
    clri_old = Ptr[float](clri_old_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            clrw_old[idx2] = max(0.0, min(1.0, 1.0 - concld[idx2] - alst[idx2]))
            clri_old[idx2] = max(0.0, min(1.0, 1.0 - concld[idx2] - ast[idx2]))


@export
def macrop_driver_forcing_prep_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    nstep: int,
    rdtime: float,
    state_t_p: cobj,
    state_qv_p: cobj,
    qc_p: cobj,
    qi_p: cobj,
    nc_p: cobj,
    ni_p: cobj,
    tcwat_p: cobj,
    qcwat_p: cobj,
    lcwat_p: cobj,
    iccwat_p: cobj,
    nlwat_p: cobj,
    niwat_p: cobj,
    cc_t_p: cobj,
    cc_qv_p: cobj,
    cc_ql_p: cobj,
    cc_qi_p: cobj,
    cc_nl_p: cobj,
    cc_ni_p: cobj,
    cc_qlst_p: cobj,
    ttend_p: cobj,
    qtend_p: cobj,
    ltend_p: cobj,
    itend_p: cobj,
    nltend_p: cobj,
    nitend_p: cobj,
    lmitend_p: cobj,
    t_inout_p: cobj,
    qv_inout_p: cobj,
    ql_inout_p: cobj,
    qi_inout_p: cobj,
    nl_inout_p: cobj,
    ni_inout_p: cobj,
):
    state_t = Ptr[float](state_t_p)
    state_qv = Ptr[float](state_qv_p)
    qc = Ptr[float](qc_p)
    qi = Ptr[float](qi_p)
    nc = Ptr[float](nc_p)
    ni = Ptr[float](ni_p)
    tcwat = Ptr[float](tcwat_p)
    qcwat = Ptr[float](qcwat_p)
    lcwat = Ptr[float](lcwat_p)
    iccwat = Ptr[float](iccwat_p)
    nlwat = Ptr[float](nlwat_p)
    niwat = Ptr[float](niwat_p)
    cc_t = Ptr[float](cc_t_p)
    cc_qv = Ptr[float](cc_qv_p)
    cc_ql = Ptr[float](cc_ql_p)
    cc_qi = Ptr[float](cc_qi_p)
    cc_nl = Ptr[float](cc_nl_p)
    cc_ni = Ptr[float](cc_ni_p)
    cc_qlst = Ptr[float](cc_qlst_p)
    ttend = Ptr[float](ttend_p)
    qtend = Ptr[float](qtend_p)
    ltend = Ptr[float](ltend_p)
    itend = Ptr[float](itend_p)
    nltend = Ptr[float](nltend_p)
    nitend = Ptr[float](nitend_p)
    lmitend = Ptr[float](lmitend_p)
    t_inout = Ptr[float](t_inout_p)
    qv_inout = Ptr[float](qv_inout_p)
    ql_inout = Ptr[float](ql_inout_p)
    qi_inout = Ptr[float](qi_inout_p)
    nl_inout = Ptr[float](nl_inout_p)
    ni_inout = Ptr[float](ni_inout_p)

    if nstep <= 1:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx2 = _idx2(i, k, pcols)
                tcwat[idx2] = state_t[idx2]
                qcwat[idx2] = state_qv[idx2]
                lcwat[idx2] = qc[idx2] + qi[idx2]
                iccwat[idx2] = qi[idx2]
                nlwat[idx2] = nc[idx2]
                niwat[idx2] = ni[idx2]
                ttend[idx2] = 0.0
                qtend[idx2] = 0.0
                ltend[idx2] = 0.0
                itend[idx2] = 0.0
                nltend[idx2] = 0.0
                nitend[idx2] = 0.0
                cc_t[idx2] = 0.0
                cc_qv[idx2] = 0.0
                cc_ql[idx2] = 0.0
                cc_qi[idx2] = 0.0
                cc_nl[idx2] = 0.0
                cc_ni[idx2] = 0.0
                cc_qlst[idx2] = 0.0
    else:
        for k in range(top_lev, pver + 1):
            for i in range(1, ncol + 1):
                idx2 = _idx2(i, k, pcols)
                ttend[idx2] = (state_t[idx2] - tcwat[idx2]) * rdtime - cc_t[idx2]
                qtend[idx2] = (state_qv[idx2] - qcwat[idx2]) * rdtime - cc_qv[idx2]
                ltend[idx2] = (qc[idx2] + qi[idx2] - lcwat[idx2]) * rdtime - (cc_ql[idx2] + cc_qi[idx2])
                itend[idx2] = (qi[idx2] - iccwat[idx2]) * rdtime - cc_qi[idx2]
                nltend[idx2] = (nc[idx2] - nlwat[idx2]) * rdtime - cc_nl[idx2]
                nitend[idx2] = (ni[idx2] - niwat[idx2]) * rdtime - cc_ni[idx2]

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            lmitend[idx2] = ltend[idx2] - itend[idx2]
            t_inout[idx2] = tcwat[idx2]
            qv_inout[idx2] = qcwat[idx2]
            ql_inout[idx2] = lcwat[idx2] - iccwat[idx2]
            qi_inout[idx2] = iccwat[idx2]
            nl_inout[idx2] = nlwat[idx2]
            ni_inout[idx2] = niwat[idx2]


@export
def macrop_driver_ptend_assign_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    tlat_p: cobj,
    qvlat_p: cobj,
    qcten_p: cobj,
    qiten_p: cobj,
    ncten_p: cobj,
    niten_p: cobj,
    ptend_s_p: cobj,
    ptend_qv_p: cobj,
    ptend_ql_p: cobj,
    ptend_qi_p: cobj,
    ptend_nl_p: cobj,
    ptend_ni_p: cobj,
):
    tlat = Ptr[float](tlat_p)
    qvlat = Ptr[float](qvlat_p)
    qcten = Ptr[float](qcten_p)
    qiten = Ptr[float](qiten_p)
    ncten = Ptr[float](ncten_p)
    niten = Ptr[float](niten_p)
    ptend_s = Ptr[float](ptend_s_p)
    ptend_qv = Ptr[float](ptend_qv_p)
    ptend_ql = Ptr[float](ptend_ql_p)
    ptend_qi = Ptr[float](ptend_qi_p)
    ptend_nl = Ptr[float](ptend_nl_p)
    ptend_ni = Ptr[float](ptend_ni_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            ptend_s[idx2] = tlat[idx2]
            ptend_qv[idx2] = qvlat[idx2]
            ptend_ql[idx2] = qcten[idx2]
            ptend_qi[idx2] = qiten[idx2]
            ptend_nl[idx2] = ncten[idx2]
            ptend_ni[idx2] = niten[idx2]


@export
def macrop_driver_wtrc_split_tend_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    qcten_p: cobj,
    qiten_p: cobj,
    pqctn_p: cobj,
    nqctn_p: cobj,
    pqitn_p: cobj,
    nqitn_p: cobj,
):
    qcten = Ptr[float](qcten_p)
    qiten = Ptr[float](qiten_p)
    pqctn = Ptr[float](pqctn_p)
    nqctn = Ptr[float](nqctn_p)
    pqitn = Ptr[float](pqitn_p)
    nqitn = Ptr[float](nqitn_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            if qcten[idx2] < 0.0:
                nqctn[idx2] = qcten[idx2]
            else:
                pqctn[idx2] = qcten[idx2]
            if qiten[idx2] < 0.0:
                nqitn[idx2] = qiten[idx2]
            else:
                pqitn[idx2] = qiten[idx2]


@export
def macrop_driver_cloud_mixing_diag_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    cld_p: cobj,
    state_ql_p: cobj,
    state_qi_p: cobj,
    mr_lsliq_p: cobj,
    mr_lsice_p: cobj,
):
    cld = Ptr[float](cld_p)
    state_ql = Ptr[float](state_ql_p)
    state_qi = Ptr[float](state_qi_p)
    mr_lsliq = Ptr[float](mr_lsliq_p)
    mr_lsice = Ptr[float](mr_lsice_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            if cld[idx2] > 0.0:
                mr_lsliq[idx2] = state_ql[idx2]
                mr_lsice[idx2] = state_qi[idx2]
            else:
                mr_lsliq[idx2] = 0.0
                mr_lsice[idx2] = 0.0


@export
def macrop_driver_store_state_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    tmelt: float,
    state_t_p: cobj,
    state_qv_p: cobj,
    state_ql_p: cobj,
    state_qi_p: cobj,
    state_nl_p: cobj,
    state_ni_p: cobj,
    tcwat_p: cobj,
    qcwat_p: cobj,
    lcwat_p: cobj,
    iccwat_p: cobj,
    nlwat_p: cobj,
    niwat_p: cobj,
    cldsice_p: cobj,
):
    state_t = Ptr[float](state_t_p)
    state_qv = Ptr[float](state_qv_p)
    state_ql = Ptr[float](state_ql_p)
    state_qi = Ptr[float](state_qi_p)
    state_nl = Ptr[float](state_nl_p)
    state_ni = Ptr[float](state_ni_p)
    tcwat = Ptr[float](tcwat_p)
    qcwat = Ptr[float](qcwat_p)
    lcwat = Ptr[float](lcwat_p)
    iccwat = Ptr[float](iccwat_p)
    nlwat = Ptr[float](nlwat_p)
    niwat = Ptr[float](niwat_p)
    cldsice = Ptr[float](cldsice_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            tcwat[idx2] = state_t[idx2]
            qcwat[idx2] = state_qv[idx2]
            lcwat[idx2] = state_ql[idx2] + state_qi[idx2]
            iccwat[idx2] = state_qi[idx2]
            nlwat[idx2] = state_nl[idx2]
            niwat[idx2] = state_ni[idx2]
            cldsice[idx2] = lcwat[idx2] * min(1.0, max(0.0, (tmelt - tcwat[idx2]) / 20.0))
