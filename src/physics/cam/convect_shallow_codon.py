@inline
def _idx2(i: int, k: int, ld1: int) -> int:
    return (i - 1) + (k - 1) * ld1


@inline
def _idx3(i: int, k: int, m: int, ld1: int, ld2: int) -> int:
    return (i - 1) + (k - 1) * ld1 + (m - 1) * ld1 * ld2


@export
def convect_shallow_select_scheme_codon(
    scheme_len: int,
    scheme_ascii_p: cobj,
    scheme_code_p: cobj,
    status_p: cobj,
):
    scheme_ascii = Ptr[int](scheme_ascii_p)
    scheme_code = Ptr[int](scheme_code_p)
    status = Ptr[int](status_p)

    status[0] = 0
    scheme_code[0] = 0

    n = scheme_len
    while n > 0 and scheme_ascii[n - 1] == 32:
        n -= 1

    if n == 2:
        c1 = scheme_ascii[0]
        c2 = scheme_ascii[1]

        if c1 >= 65 and c1 <= 90:
            c1 += 32
        if c2 >= 65 and c2 <= 90:
            c2 += 32

        if c1 == 117 and c2 == 119:
            scheme_code[0] = 3
            return

    if n == 3:
        c1 = scheme_ascii[0]
        c2 = scheme_ascii[1]
        c3 = scheme_ascii[2]

        if c1 >= 65 and c1 <= 90:
            c1 += 32
        if c2 >= 65 and c2 <= 90:
            c2 += 32
        if c3 >= 65 and c3 <= 90:
            c3 += 32

        if c1 == 111 and c2 == 102 and c3 == 102:
            scheme_code[0] = 1
            return

    if n == 4:
        c1 = scheme_ascii[0]
        c2 = scheme_ascii[1]
        c3 = scheme_ascii[2]
        c4 = scheme_ascii[3]

        if c1 >= 65 and c1 <= 90:
            c1 += 32
        if c2 >= 65 and c2 <= 90:
            c2 += 32
        if c3 >= 65 and c3 <= 90:
            c3 += 32
        if c4 >= 65 and c4 <= 90:
            c4 += 32

        if c1 == 104 and c2 == 97 and c3 == 99 and c4 == 107:
            scheme_code[0] = 2
            return

    if n == 6:
        c1 = scheme_ascii[0]
        c2 = scheme_ascii[1]
        c3 = scheme_ascii[2]
        c4 = scheme_ascii[3]
        c5 = scheme_ascii[4]
        c6 = scheme_ascii[5]

        if c1 >= 65 and c1 <= 90:
            c1 += 32
        if c2 >= 65 and c2 <= 90:
            c2 += 32
        if c3 >= 65 and c3 <= 90:
            c3 += 32
        if c4 >= 65 and c4 <= 90:
            c4 += 32
        if c5 >= 65 and c5 <= 90:
            c5 += 32
        if c6 >= 65 and c6 <= 90:
            c6 += 32

        if c1 == 117 and c2 == 110 and c3 == 105 and c4 == 99 and c5 == 111 and c6 == 110:
            scheme_code[0] = 4
            return

    if n == 9:
        c1 = scheme_ascii[0]
        c2 = scheme_ascii[1]
        c3 = scheme_ascii[2]
        c4 = scheme_ascii[3]
        c5 = scheme_ascii[4]
        c6 = scheme_ascii[5]
        c7 = scheme_ascii[6]
        c8 = scheme_ascii[7]
        c9 = scheme_ascii[8]

        if c1 >= 65 and c1 <= 90:
            c1 += 32
        if c2 >= 65 and c2 <= 90:
            c2 += 32
        if c3 >= 65 and c3 <= 90:
            c3 += 32
        if c4 >= 65 and c4 <= 90:
            c4 += 32
        if c5 >= 65 and c5 <= 90:
            c5 += 32
        if c6 >= 65 and c6 <= 90:
            c6 += 32
        if c7 >= 65 and c7 <= 90:
            c7 += 32
        if c8 >= 65 and c8 <= 90:
            c8 += 32
        if c9 >= 65 and c9 <= 90:
            c9 += 32

        if (
            c1 == 99 and c2 == 108 and c3 == 117 and c4 == 98 and c5 == 98
            and c6 == 95 and c7 == 115 and c8 == 103 and c9 == 115
        ):
            scheme_code[0] = 1
            return

    status[0] = 1


@export
def convect_shallow_init_shell_codon(
    ncol: int,
    pcols: int,
    tpert_p: cobj,
    landfracdum_p: cobj,
):
    tpert = Ptr[float](tpert_p)
    landfracdum = Ptr[float](landfracdum_p)

    for i in range(1, ncol + 1):
        idx1 = i - 1
        tpert[idx1] = 0.0
        landfracdum[idx1] = 0.0


@export
def convect_shallow_diag_shell_codon(
    mode: int,
    ncol: int,
    pcols: int,
    pver: int,
    ixcldliq: int,
    ixcldice: int,
    ztodt: float,
    latvap: float,
    latice: float,
    zvir: float,
    state_s_p: cobj,
    state_t_p: cobj,
    state_q_p: cobj,
    sat_rh_p: cobj,
    sl_p: cobj,
    qt_p: cobj,
    slv_p: cobj,
    t_precu_p: cobj,
    rh_precu_p: cobj,
    tten_p: cobj,
    rhten_p: cobj,
):
    state_s = Ptr[float](state_s_p)
    state_t = Ptr[float](state_t_p)
    state_q = Ptr[float](state_q_p)
    sat_rh = Ptr[float](sat_rh_p)
    sl = Ptr[float](sl_p)
    qt = Ptr[float](qt_p)
    slv = Ptr[float](slv_p)
    t_precu = Ptr[float](t_precu_p)
    rh_precu = Ptr[float](rh_precu_p)
    tten = Ptr[float](tten_p)
    rhten = Ptr[float](rhten_p)

    if mode == 1 or mode == 3:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, pcols)
                sl[idx] = (
                    state_s[idx]
                    - latvap * state_q[_idx3(i, k, ixcldliq, pcols, pver)]
                    - (latvap + latice) * state_q[_idx3(i, k, ixcldice, pcols, pver)]
                )
                qt[idx] = (
                    state_q[_idx3(i, k, 1, pcols, pver)]
                    + state_q[_idx3(i, k, ixcldliq, pcols, pver)]
                    + state_q[_idx3(i, k, ixcldice, pcols, pver)]
                )
                slv[idx] = sl[idx] * (1.0 + zvir * qt[idx])
                if mode == 1:
                    t_precu[idx] = state_t[idx]
    elif mode == 2:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, pcols)
                rh_precu[idx] = state_q[_idx3(i, k, 1, pcols, pver)] / sat_rh[idx] * 100.0
    elif mode == 4:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, pcols)
                sat_rh[idx] = state_q[_idx3(i, k, 1, pcols, pver)] / sat_rh[idx] * 100.0
                tten[idx] = (state_t[idx] - t_precu[idx]) / ztodt
                rhten[idx] = (sat_rh[idx] - rh_precu[idx]) / ztodt


@export
def convect_shallow_uw_post_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    pcnst: int,
    wtrc_nwset: int,
    latvap: float,
    cpair: float,
    state_pmid_p: cobj,
    cmfmc_p: cobj,
    cmfmc2_p: cobj,
    cnt_p: cobj,
    cnt2_p: cobj,
    cnb_p: cobj,
    cnb2_p: cobj,
    pcnt_p: cobj,
    pcnb_p: cobj,
    qc_p: cobj,
    qc2_p: cobj,
    rliq_p: cobj,
    rliq2_p: cobj,
    wtqc_p: cobj,
    wtdlf_p: cobj,
    freqsh_p: cobj,
    icwmr_p: cobj,
    iccmr_uw_p: cobj,
    rprdsh_p: cobj,
    cmfdqs_p: cobj,
    ptend_q_p: cobj,
    ptend_tracer_p: cobj,
    cmfsl_p: cobj,
    cmflq_p: cobj,
    slflx_p: cobj,
    qtflx_p: cobj,
    rprddp_p: cobj,
    rprdtot_p: cobj,
    ptend_s_p: cobj,
    ftem_p: cobj,
    liq_type_p: cobj,
    ice_type_p: cobj,
):
    state_pmid = Ptr[float](state_pmid_p)
    cmfmc = Ptr[float](cmfmc_p)
    cmfmc2 = Ptr[float](cmfmc2_p)
    cnt = Ptr[float](cnt_p)
    cnt2 = Ptr[float](cnt2_p)
    cnb = Ptr[float](cnb_p)
    cnb2 = Ptr[float](cnb2_p)
    pcnt = Ptr[float](pcnt_p)
    pcnb = Ptr[float](pcnb_p)
    qc = Ptr[float](qc_p)
    qc2 = Ptr[float](qc2_p)
    rliq = Ptr[float](rliq_p)
    rliq2 = Ptr[float](rliq2_p)
    wtqc = Ptr[float](wtqc_p)
    wtdlf = Ptr[float](wtdlf_p)
    freqsh = Ptr[float](freqsh_p)
    icwmr = Ptr[float](icwmr_p)
    iccmr_uw = Ptr[float](iccmr_uw_p)
    rprdsh = Ptr[float](rprdsh_p)
    cmfdqs = Ptr[float](cmfdqs_p)
    ptend_q = Ptr[float](ptend_q_p)
    ptend_tracer = Ptr[float](ptend_tracer_p)
    cmfsl = Ptr[float](cmfsl_p)
    cmflq = Ptr[float](cmflq_p)
    slflx = Ptr[float](slflx_p)
    qtflx = Ptr[float](qtflx_p)
    rprddp = Ptr[float](rprddp_p)
    rprdtot = Ptr[float](rprdtot_p)
    ptend_s = Ptr[float](ptend_s_p)
    ftem = Ptr[float](ftem_p)
    liq_type = Ptr[int](liq_type_p)
    ice_type = Ptr[int](ice_type_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            icwmr[idx] = iccmr_uw[idx]
            rprdsh[idx] = rprdsh[idx] + cmfdqs[idx]

    for m in range(4, pcnst + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx3(i, k, m, pcols, pver)
                ptend_q[idx] = ptend_tracer[idx]

    for k in range(1, pverp + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            cmfsl[idx] = slflx[idx]
            cmflq[idx] = qtflx[idx] * latvap

    for i in range(1, pcols + 1):
        freqsh[i - 1] = 0.0
    for i in range(1, ncol + 1):
        active = 0
        for k in range(1, pver + 1):
            if cmfmc2[_idx2(i, k, pcols)] > 0.0:
                active = 1
                break
        if active == 0:
            freqsh[i - 1] = 1.0

    for k in range(1, pverp + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            cmfmc[idx] = cmfmc[idx] + cmfmc2[idx]

    for i in range(1, ncol + 1):
        if cnt2[i - 1] < cnt[i - 1]:
            cnt[i - 1] = cnt2[i - 1]
        if cnb2[i - 1] > cnb[i - 1]:
            cnb[i - 1] = cnb2[i - 1]
        pcnt[i - 1] = state_pmid[_idx2(i, int(cnt[i - 1]), pcols)]
        pcnb[i - 1] = state_pmid[_idx2(i, int(cnb[i - 1]), pcols)]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            qc[idx] = qc[idx] + qc2[idx]

    for i in range(1, ncol + 1):
        rliq[i - 1] = rliq[i - 1] + rliq2[i - 1]

    for m in range(1, wtrc_nwset + 1):
        liq = liq_type[m - 1]
        ice = ice_type[m - 1]
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx_out = _idx3(i, k, m, pcols, pver)
                tmp_sum = (
                    wtqc[_idx3(i, k, liq, pcols, pver)]
                    + wtqc[_idx3(i, k, ice, pcols, pver)]
                )
                wtdlf[idx_out] = wtdlf[idx_out] + tmp_sum

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            rprdtot[idx] = rprdsh[idx] + rprddp[idx]
            ftem[idx] = ptend_s[idx] / cpair


@export
def convect_shallow_postmerge_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    pcnst: int,
    wtrc_nwset: int,
    state_pmid_p: cobj,
    cmfmc_p: cobj,
    cmfmc2_p: cobj,
    cnt_p: cobj,
    cnt2_p: cobj,
    cnb_p: cobj,
    cnb2_p: cobj,
    pcnt_p: cobj,
    pcnb_p: cobj,
    qc_p: cobj,
    qc2_p: cobj,
    rliq_p: cobj,
    rliq2_p: cobj,
    wtqc_p: cobj,
    wtdlf_p: cobj,
    freqsh_p: cobj,
    liq_type_p: cobj,
    ice_type_p: cobj,
):
    state_pmid = Ptr[float](state_pmid_p)
    cmfmc = Ptr[float](cmfmc_p)
    cmfmc2 = Ptr[float](cmfmc2_p)
    cnt = Ptr[float](cnt_p)
    cnt2 = Ptr[float](cnt2_p)
    cnb = Ptr[float](cnb_p)
    cnb2 = Ptr[float](cnb2_p)
    pcnt = Ptr[float](pcnt_p)
    pcnb = Ptr[float](pcnb_p)
    qc = Ptr[float](qc_p)
    qc2 = Ptr[float](qc2_p)
    rliq = Ptr[float](rliq_p)
    rliq2 = Ptr[float](rliq2_p)
    wtqc = Ptr[float](wtqc_p)
    wtdlf = Ptr[float](wtdlf_p)
    freqsh = Ptr[float](freqsh_p)
    liq_type = Ptr[int](liq_type_p)
    ice_type = Ptr[int](ice_type_p)

    for i in range(1, pcols + 1):
        freqsh[i - 1] = 0.0
    for i in range(1, ncol + 1):
        active = 0
        for k in range(1, pver + 1):
            if cmfmc2[_idx2(i, k, pcols)] > 0.0:
                active = 1
                break
        if active == 0:
            freqsh[i - 1] = 1.0

    for k in range(1, pverp + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            cmfmc[idx] = cmfmc[idx] + cmfmc2[idx]

    for i in range(1, ncol + 1):
        if cnt2[i - 1] < cnt[i - 1]:
            cnt[i - 1] = cnt2[i - 1]
        if cnb2[i - 1] > cnb[i - 1]:
            cnb[i - 1] = cnb2[i - 1]
        pcnt[i - 1] = state_pmid[_idx2(i, int(cnt[i - 1]), pcols)]
        pcnb[i - 1] = state_pmid[_idx2(i, int(cnb[i - 1]), pcols)]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            qc[idx] = qc[idx] + qc2[idx]

    for i in range(1, ncol + 1):
        rliq[i - 1] = rliq[i - 1] + rliq2[i - 1]

    for m in range(1, wtrc_nwset + 1):
        liq = liq_type[m - 1]
        ice = ice_type[m - 1]
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx_out = _idx3(i, k, m, pcols, pver)
                tmp_sum = (
                    wtqc[_idx3(i, k, liq, pcols, pver)]
                    + wtqc[_idx3(i, k, ice, pcols, pver)]
                )
                wtdlf[idx_out] = wtdlf[idx_out] + tmp_sum
