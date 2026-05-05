from math import log


@inline
def _state_q_idx(i: int, k: int, m: int, pcols: int, pver: int) -> int:
    """state%q(pcols, pver, pcnst)"""
    return (i - 1) + (k - 1) * pcols + (m - 1) * pcols * pver


@inline
def _field1_idx(i: int) -> int:
    """state%lat/state%phis(pcols)"""
    return i - 1


@inline
def _field2_idx(i: int, k: int, pcols: int) -> int:
    """state%pdel(pcols, pver)"""
    return (i - 1) + (k - 1) * pcols


@inline
def _wtrc_q1q2_3d_idx(i: int, k: int, m: int, pcols: int, pver: int) -> int:
    """work(pcols, pver, wtrc_nwset)"""
    return i + k * pcols + m * pcols * pver


@inline
def _wtrc_q1q2_q_idx(i: int, k: int, m: int, pcols: int, pver: int) -> int:
    """q/dqdt/wtrprd(pcols, pver, pcnst)"""
    return i + k * pcols + m * pcols * pver


@export
def wtrc_q1q2_init_qhat_shell_codon(
    lengath: int,
    pcols: int,
    pver: int,
    pcnst: int,
    pwtype: int,
    wtrc_nwset: int,
    msg: int,
    iwtvap: int,
    ideep_p: cobj,
    wtrc_iatype_p: cobj,
    q_p: cobj,
    dqdt_p: cobj,
    wtrprd_p: cobj,
    wtdlf_p: cobj,
    qhat_p: cobj,
    qu_p: cobj,
    qd_p: cobj,
    wtcu_p: cobj,
    wtevp_p: cobj,
    wthmn_p: cobj,
    hsat_p: cobj,
    qst_p: cobj,
    gamma_p: cobj,
    hsthat_p: cobj,
    qsthat_p: cobj,
    gamhat_p: cobj,
    hu_p: cobj,
    hd_p: cobj,
    totpcp_p: cobj,
    totevp_p: cobj,
    ru_p: cobj,
    rd_p: cobj,
    ql_p: cobj,
):
    ideep = Ptr[int](ideep_p)
    wtrc_iatype = Ptr[int](wtrc_iatype_p)
    q = Ptr[float](q_p)
    dqdt = Ptr[float](dqdt_p)
    wtrprd = Ptr[float](wtrprd_p)
    wtdlf = Ptr[float](wtdlf_p)
    qhat = Ptr[float](qhat_p)
    qu = Ptr[float](qu_p)
    qd = Ptr[float](qd_p)
    wtcu = Ptr[float](wtcu_p)
    wtevp = Ptr[float](wtevp_p)
    wthmn = Ptr[float](wthmn_p)
    hsat = Ptr[float](hsat_p)
    qst = Ptr[float](qst_p)
    gamma = Ptr[float](gamma_p)
    hsthat = Ptr[float](hsthat_p)
    qsthat = Ptr[float](qsthat_p)
    gamhat = Ptr[float](gamhat_p)
    hu = Ptr[float](hu_p)
    hd = Ptr[float](hd_p)
    totpcp = Ptr[float](totpcp_p)
    totevp = Ptr[float](totevp_p)
    ru = Ptr[float](ru_p)
    rd = Ptr[float](rd_p)
    ql = Ptr[float](ql_p)

    m = 0
    while m < wtrc_nwset:
        vap_idx = wtrc_iatype[m + (iwtvap - 1) * wtrc_nwset] - 1
        wt = 0
        while wt < pwtype:
            q_idx = wtrc_iatype[m + wt * wtrc_nwset] - 1
            k = 0
            while k < pver:
                i = 0
                while i < pcols:
                    dqdt[_wtrc_q1q2_q_idx(i, k, q_idx, pcols, pver)] = 0.0
                    i += 1
                k += 1
            wt += 1

        k = 0
        while k < pver:
            i = 0
            while i < pcols:
                qidx = _wtrc_q1q2_q_idx(i, k, vap_idx, pcols, pver)
                widx = _wtrc_q1q2_3d_idx(i, k, m, pcols, pver)
                wtrprd[qidx] = 0.0
                wtdlf[widx] = 0.0
                qhat[widx] = q[qidx]
                qu[widx] = 0.0
                qd[widx] = 0.0
                ql[widx] = 0.0
                wtcu[widx] = 0.0
                wtevp[widx] = 0.0
                wthmn[widx] = 0.0
                hsat[widx] = 0.0
                qst[widx] = 0.0
                gamma[widx] = 0.0
                hsthat[widx] = 0.0
                qsthat[widx] = 0.0
                gamhat[widx] = 0.0
                hu[widx] = 0.0
                hd[widx] = 0.0
                ru[widx] = 0.0
                rd[widx] = 0.0
                i += 1
            i = 0
            while i < lengath:
                ii = ideep[i] - 1
                widx = _wtrc_q1q2_3d_idx(i, k, m, pcols, pver)
                qidx = _wtrc_q1q2_q_idx(ii, k, vap_idx, pcols, pver)
                qu[widx] = q[qidx]
                qd[widx] = q[qidx]
                i += 1
            k += 1

        i = 0
        while i < pcols:
            totpcp[i + m * pcols] = 0.0
            totevp[i + m * pcols] = 0.0
            i += 1
        m += 1

    m = 0
    while m < wtrc_nwset:
        vap_idx = wtrc_iatype[m + (iwtvap - 1) * wtrc_nwset] - 1
        kk = msg + 2
        while kk <= pver:
            k = kk - 1
            km1 = kk - 2
            i = 0
            while i < lengath:
                ii = ideep[i] - 1
                idx = _wtrc_q1q2_3d_idx(i, k, m, pcols, pver)
                qidx = _wtrc_q1q2_q_idx(ii, k, vap_idx, pcols, pver)
                qm1idx = _wtrc_q1q2_q_idx(ii, km1, vap_idx, pcols, pver)
                qdifr = 0.0
                if q[qidx] > 0.0 or q[qm1idx] > 0.0:
                    qdifr = abs((q[qidx] - q[qm1idx]) / max(q[qm1idx], q[qidx]))
                if qdifr > 1.0e-6 and qdifr != 1.0:
                    qhat[idx] = log(q[qm1idx] / q[qidx]) * q[qm1idx] * q[qidx] / (q[qm1idx] - q[qidx])
                else:
                    qhat[idx] = 0.5 * (q[qidx] + q[qm1idx])
                i += 1
            kk += 1
        m += 1


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
) -> int:
    """process_rates(pcols, pver, pwtype, pwtype, pwtype)"""
    return (
        (i - 1)
        + (k - 1) * pcols
        + (idsttype - 1) * pcols * pver
        + (isrctype - 1) * pcols * pver * pwtype
        + (rtype - 1) * pcols * pver * pwtype * pwtype
    )


@inline
def _iatype_idx(m: int, p: int, wtrc_nwset: int) -> int:
    """wtrc_iatype64(wtrc_nwset, pwtype)"""
    return (m - 1) + (p - 1) * wtrc_nwset


@inline
def _iawset_idx(itype: int, iwset: int, pwtype: int) -> int:
    """wtrc_iawset64(pwtype, wtrc_nwset)"""
    return (itype - 1) + (iwset - 1) * pwtype


@inline
def _bulk_idx(p: int) -> int:
    """wtrc_bulk_indices64(pwtype)"""
    return p - 1


@inline
def _spec_idx(m: int) -> int:
    """rstd(pwtspec) and iwspec64(pcnst)"""
    return m - 1


@inline
def _wtrc_ratio(ispec: int, qtrc: float, qtot: float, wtrc_qmin: float, rstd) -> float:
    if abs(qtot) < wtrc_qmin:
        return rstd[_spec_idx(ispec)]
    return qtrc / qtot


@inline
def _field2p_idx(i: int, k: int, pcols: int) -> int:
    """field(pcols, pverp)"""
    return (i - 1) + (k - 1) * pcols


@inline
def _field3_pverp_idx(i: int, k: int, m: int, pcols: int, pverp: int) -> int:
    """field(pcols, pverp, nwvap)"""
    return (i - 1) + (k - 1) * pcols + (m - 1) * pcols * pverp


@export
def wtrc_precip_evap_init_shell_codon(
    pcols: int,
    pver: int,
    pverp: int,
    nwvap: int,
    evp_p: cobj,
    mlt_p: cobj,
    frz_p: cobj,
    sub_p: cobj,
    fice_p: cobj,
    fsnow_p: cobj,
    rnbulk_p: cobj,
    flxpr_p: cobj,
    flxsn_p: cobj,
    totrnfx_p: cobj,
    dz_p: cobj,
):
    evp = Ptr[float](evp_p)
    mlt = Ptr[float](mlt_p)
    frz = Ptr[float](frz_p)
    sub = Ptr[float](sub_p)
    fice = Ptr[float](fice_p)
    fsnow = Ptr[float](fsnow_p)
    rnbulk = Ptr[float](rnbulk_p)
    flxpr = Ptr[float](flxpr_p)
    flxsn = Ptr[float](flxsn_p)
    totrnfx = Ptr[float](totrnfx_p)
    dz = Ptr[float](dz_p)
    plane = pcols * pver
    flux_plane = pcols * pverp

    m = 0
    while m < nwvap:
        moff = m * plane
        foff = m * flux_plane
        k = 0
        while k < pver:
            i = 0
            while i < pcols:
                idx = i + k * pcols
                evp[idx + moff] = 0.0
                mlt[idx + moff] = 0.0
                frz[idx + moff] = 0.0
                sub[idx + moff] = 0.0
                i += 1
            k += 1
        k = 0
        while k < pverp:
            i = 0
            while i < pcols:
                idx = i + k * pcols
                flxpr[idx + foff] = 0.0
                flxsn[idx + foff] = 0.0
                i += 1
            k += 1
        m += 1

    k = 0
    while k < pver:
        i = 0
        while i < pcols:
            idx = i + k * pcols
            fice[idx] = 0.0
            fsnow[idx] = 0.0
            rnbulk[idx] = 0.0
            dz[idx] = 0.0
            i += 1
        k += 1

    k = 0
    while k < pverp:
        i = 0
        while i < pcols:
            totrnfx[i + k * pcols] = 0.0
            i += 1
        k += 1


@export
def wtrc_precip_evap_prep_shell_codon(
    pcols: int,
    pver: int,
    evpbulk_p: cobj,
    subbulk_p: cobj,
    rnbulk_p: cobj,
):
    evpbulk = Ptr[float](evpbulk_p)
    subbulk = Ptr[float](subbulk_p)
    rnbulk = Ptr[float](rnbulk_p)

    k = 0
    while k < pver:
        i = 0
        while i < pcols:
            idx = i + k * pcols
            rnbulk[idx] = evpbulk[idx] - subbulk[idx]
            i += 1
        k += 1


@export
def wtrc_precip_evap_tail_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    pcnst: int,
    nwvap: int,
    wtrc_nwset: int,
    iwtvap: int,
    wtrc_qmin: float,
    prec_p: cobj,
    snow_p: cobj,
    flxpr_p: cobj,
    flxsn_p: cobj,
    wtrc_iatype_p: cobj,
    iwspec_p: cobj,
    rstd_p: cobj,
):
    prec = Ptr[float](prec_p)
    snow = Ptr[float](snow_p)
    flxpr = Ptr[float](flxpr_p)
    flxsn = Ptr[float](flxsn_p)
    wtrc_iatype = Ptr[int](wtrc_iatype_p)
    iwspec = Ptr[int](iwspec_p)
    rstd = Ptr[float](rstd_p)

    m = 1
    while m <= wtrc_nwset:
        trc_m = wtrc_iatype[_iatype_idx(m, iwtvap, wtrc_nwset)]
        i = 1
        while i <= ncol:
            prec[_field2_idx(i, trc_m, pcols)] = flxpr[_field3_pverp_idx(i, pverp, m, pcols, pverp)] / 1000.0
            snow[_field2_idx(i, trc_m, pcols)] = flxsn[_field3_pverp_idx(i, pverp, m, pcols, pverp)] / 1000.0
            i += 1
        m += 1

    base1 = wtrc_iatype[_iatype_idx(1, iwtvap, wtrc_nwset)]
    base2 = wtrc_iatype[_iatype_idx(2, iwtvap, wtrc_nwset)]
    i = 1
    while i <= ncol:
        pmass0 = prec[_field2_idx(i, base2, pcols)]
        smass0 = snow[_field2_idx(i, base2, pcols)]
        pdiff = pmass0 - prec[_field2_idx(i, base1, pcols)]
        sdiff = smass0 - snow[_field2_idx(i, base1, pcols)]
        m = 2
        while m <= wtrc_nwset:
            trc_m = wtrc_iatype[_iatype_idx(m, iwtvap, wtrc_nwset)]
            prec_idx = _field2_idx(i, trc_m, pcols)
            snow_idx = _field2_idx(i, trc_m, pcols)
            rd = _wtrc_ratio(
                iwspec[_spec_idx(trc_m)],
                prec[prec_idx] - snow[snow_idx],
                pmass0,
                wtrc_qmin,
                rstd,
            )
            prec[prec_idx] = max(0.0, prec[prec_idx] - rd * pdiff)
            rd = _wtrc_ratio(
                iwspec[_spec_idx(trc_m)],
                snow[snow_idx],
                smass0,
                wtrc_qmin,
                rstd,
            )
            snow[snow_idx] = max(0.0, snow[snow_idx] - rd * sdiff)
            if prec[prec_idx] > 10.0 * prec[_field2_idx(i, base1, pcols)]:
                prec[prec_idx] = prec[_field2_idx(i, base1, pcols)]
                snow[snow_idx] = snow[_field2_idx(i, base1, pcols)]
            m += 1
        i += 1


@export
def wtrc_mass_fixer_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    pwtype: int,
    wtrc_nwset: int,
    isphdo: int,
    wisotope: int,
    wtrc_qmin: float,
    wtrc_limiter_18O_hgh: float,
    wtrc_limiter_18O_low: float,
    wtrc_limiter_HDO_hgh: float,
    wtrc_limiter_HDO_low: float,
    wtrc_limiter_phis_crit: float,
    radtodeg: float,
    state_q_p: cobj,
    state_lat_p: cobj,
    state_phis_p: cobj,
    wtrc_iatype_p: cobj,
    wtrc_bulk_indices_p: cobj,
    iwspec_p: cobj,
    rstd_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    state_lat = Ptr[float](state_lat_p)
    state_phis = Ptr[float](state_phis_p)
    wtrc_iatype = Ptr[int](wtrc_iatype_p)
    wtrc_bulk_indices = Ptr[int](wtrc_bulk_indices_p)
    iwspec = Ptr[int](iwspec_p)
    rstd = Ptr[float](rstd_p)

    for i in range(1, ncol + 1):
        for k in range(1, pver + 1):
            for p in range(1, pwtype + 1):
                base_m = wtrc_iatype[_iatype_idx(1, p, wtrc_nwset)]
                bulk_m = wtrc_bulk_indices[_bulk_idx(p)]
                base_idx = _state_q_idx(i, k, base_m, pcols, pver)
                bulk_idx = _state_q_idx(i, k, bulk_m, pcols, pver)
                if state_q[base_idx] != state_q[bulk_idx]:
                    diff = state_q[base_idx] - state_q[bulk_idx]
                    oval = state_q[base_idx]
                    for m in range(1, wtrc_nwset + 1):
                        trc_m = wtrc_iatype[_iatype_idx(m, p, wtrc_nwset)]
                        trc_idx = _state_q_idx(i, k, trc_m, pcols, pver)
                        ratio = _wtrc_ratio(
                            iwspec[_spec_idx(trc_m)],
                            state_q[trc_idx],
                            oval,
                            wtrc_qmin,
                            rstd,
                        )
                        state_q[trc_idx] = state_q[trc_idx] - ratio * diff

    if wisotope != 0:
        for i in range(1, ncol + 1):
            for k in range(1, pver + 1):
                for p in range(1, pwtype + 1):
                    base_m = wtrc_iatype[_iatype_idx(2, p, wtrc_nwset)]
                    bulk_m = wtrc_bulk_indices[_bulk_idx(p)]
                    base_idx = _state_q_idx(i, k, base_m, pcols, pver)
                    bulk_idx = _state_q_idx(i, k, bulk_m, pcols, pver)
                    if state_q[base_idx] != state_q[bulk_idx]:
                        diff = state_q[base_idx] - state_q[bulk_idx]
                        oval = state_q[base_idx]
                        for m in range(2, wtrc_nwset + 1):
                            trc_m = wtrc_iatype[_iatype_idx(m, p, wtrc_nwset)]
                            trc_idx = _state_q_idx(i, k, trc_m, pcols, pver)
                            ratio = _wtrc_ratio(
                                iwspec[_spec_idx(trc_m)],
                                state_q[trc_idx],
                                oval,
                                wtrc_qmin,
                                rstd,
                            )
                            state_q[trc_idx] = state_q[trc_idx] - ratio * diff

    for i in range(1, ncol + 1):
        wtlat = state_lat[_field1_idx(i)] * radtodeg
        wtphis = state_phis[_field1_idx(i)]
        for k in range(1, pver + 1):
            for p in range(1, pwtype + 1):
                base_m = wtrc_iatype[_iatype_idx(1, p, wtrc_nwset)]
                base_idx = _state_q_idx(i, k, base_m, pcols, pver)
                for m in range(2, wtrc_nwset + 1):
                    trc_m = wtrc_iatype[_iatype_idx(m, p, wtrc_nwset)]
                    trc_idx = _state_q_idx(i, k, trc_m, pcols, pver)
                    if m == isphdo:
                        diff_limit = wtrc_limiter_HDO_hgh
                    else:
                        diff_limit = wtrc_limiter_18O_hgh

                    if state_q[trc_idx] > diff_limit * state_q[base_idx]:
                        state_q[trc_idx] = state_q[base_idx]
                    else:
                        if abs(wtlat) < 60.0 and wtphis > wtrc_limiter_phis_crit:
                            if m == isphdo:
                                diff_limit = wtrc_limiter_HDO_low + 0.005 * (k - pver)
                            else:
                                diff_limit = wtrc_limiter_18O_low + 0.005 * (k - pver)
                            if state_q[trc_idx] < diff_limit * state_q[base_idx]:
                                state_q[trc_idx] = state_q[base_idx]


@export
def wtrc_check_h2o_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    pwtype: int,
    wtrc_nwset: int,
    wisotope: int,
    wtrc_check_total_h2o: int,
    wtrc_qchkmin: float,
    dtime: float,
    pstate_q_p: cobj,
    pstate_pdel_p: cobj,
    qloc_p: cobj,
    wtrc_bulk_indices_p: cobj,
    wtrc_iawset_p: cobj,
    iwspec_p: cobj,
    rstd_p: cobj,
    ptend_present: int,
    ptend_q_p: cobj,
    result_p: cobj,
    issue_p: cobj,
):
    pstate_q = Ptr[float](pstate_q_p)
    pstate_pdel = Ptr[float](pstate_pdel_p)
    qloc = Ptr[float](qloc_p)
    wtrc_bulk_indices = Ptr[int](wtrc_bulk_indices_p)
    wtrc_iawset = Ptr[int](wtrc_iawset_p)
    iwspec = Ptr[int](iwspec_p)
    rstd = Ptr[float](rstd_p)
    ptend_q = Ptr[float](ptend_q_p)
    result = Ptr[int](result_p)
    issue = Ptr[int](issue_p)

    ok = 1
    issue_found = 0

    if wtrc_check_total_h2o != 0:
        nwset = 1 if wisotope != 0 else wtrc_nwset

        for icol in range(1, ncol + 1):
            sbmass = 0.0

            for itype in range(1, pwtype + 1):
                bulk_m = wtrc_bulk_indices[_bulk_idx(itype)]
                bulk_mass = 0.0

                for k in range(1, pver + 1):
                    pdel_idx = _field2_idx(icol, k, pcols)
                    bulk_mass += (
                        pstate_q[_state_q_idx(icol, k, bulk_m, pcols, pver)]
                        * pstate_pdel[pdel_idx]
                    )

                if ptend_present != 0:
                    ptend_mass = 0.0
                    for k in range(1, pver + 1):
                        pdel_idx = _field2_idx(icol, k, pcols)
                        ptend_mass += (
                            ptend_q[_state_q_idx(icol, k, bulk_m, pcols, pver)]
                            * pstate_pdel[pdel_idx]
                        )
                    bulk_mass += ptend_mass * dtime

                sbmass += bulk_mass

            for iwset in range(1, nwset + 1):
                base_m = wtrc_iawset[_iawset_idx(1, iwset, pwtype)]
                scale = 1.0 / rstd[_spec_idx(iwspec[_spec_idx(base_m)])]
                stmass = 0.0

                for itype in range(1, pwtype + 1):
                    trc_m = wtrc_iawset[_iawset_idx(itype, iwset, pwtype)]
                    for k in range(1, pver + 1):
                        pdel_idx = _field2_idx(icol, k, pcols)
                        stmass += (
                            qloc[_state_q_idx(icol, k, trc_m, pcols, pver)]
                            * pstate_pdel[pdel_idx]
                        )

                if (stmass < 0.0) or (sbmass < 0.0):
                    ok = 0
                    issue_found = 1
                elif sbmass > 0.0:
                    ok = 0
                    if abs(((stmass * scale) - sbmass) / sbmass) > wtrc_qchkmin:
                        issue_found = 1

    result[0] = ok
    issue[0] = issue_found


@export
def wtrc_clear_precip_codon(ncol: int, srfpcp_p: cobj):
    srfpcp = Ptr[float](srfpcp_p)

    for i in range(1, ncol + 1):
        srfpcp[_field1_idx(i)] = 0.0


@export
def wtrc_diagnose_bulk_precip_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    bulk_idx: int,
    ptend_q_p: cobj,
):
    ptend_q = Ptr[float](ptend_q_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            ptend_q[_state_q_idx(i, k, bulk_idx, pcols, pver)] = 0.0


@export
def wtrc_add_rates_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pwtype: int,
    top_lev: int,
    isrctype: int,
    idsttype: int,
    rtype: int,
    do_reverse_present: int,
    do_reverse: int,
    process_rates_p: cobj,
    rate_p: cobj,
):
    process_rates = Ptr[float](process_rates_p)
    rate = Ptr[float](rate_p)

    ldo_reverse = True
    if do_reverse_present != 0:
        ldo_reverse = do_reverse != 0

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            rate_val = rate[_field2_idx(i, k, pcols)]
            dst_idx = _process_rates_idx(
                i, k, idsttype, isrctype, rtype, pcols, pver, pwtype
            )
            process_rates[dst_idx] = process_rates[dst_idx] + rate_val

            if isrctype != idsttype and ldo_reverse:
                src_idx = _process_rates_idx(
                    i, k, isrctype, idsttype, rtype, pcols, pver, pwtype
                )
                process_rates[src_idx] = process_rates[src_idx] - rate_val


@export
def wtrc_batch_mass_fixer_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    pwtype: int,
    wtrc_nwset: int,
    isphdo: int,
    wisotope: int,
    wtrc_qmin: float,
    wtrc_limiter_18O_hgh: float,
    wtrc_limiter_18O_low: float,
    wtrc_limiter_HDO_hgh: float,
    wtrc_limiter_HDO_low: float,
    wtrc_limiter_phis_crit: float,
    radtodeg: float,
    state_q_p: cobj,
    state_lat_p: cobj,
    state_phis_p: cobj,
    wtrc_iatype_p: cobj,
    wtrc_bulk_indices_p: cobj,
    iwspec_p: cobj,
    rstd_p: cobj,
):
    wtrc_mass_fixer_codon(
        ncol,
        pcols,
        pver,
        pcnst,
        pwtype,
        wtrc_nwset,
        isphdo,
        wisotope,
        wtrc_qmin,
        wtrc_limiter_18O_hgh,
        wtrc_limiter_18O_low,
        wtrc_limiter_HDO_hgh,
        wtrc_limiter_HDO_low,
        wtrc_limiter_phis_crit,
        radtodeg,
        state_q_p,
        state_lat_p,
        state_phis_p,
        wtrc_iatype_p,
        wtrc_bulk_indices_p,
        iwspec_p,
        rstd_p,
    )


@export
def wtrc_batch_check_h2o_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    pwtype: int,
    wtrc_nwset: int,
    wisotope: int,
    wtrc_check_total_h2o: int,
    wtrc_qchkmin: float,
    dtime: float,
    pstate_q_p: cobj,
    pstate_pdel_p: cobj,
    qloc_p: cobj,
    wtrc_bulk_indices_p: cobj,
    wtrc_iawset_p: cobj,
    iwspec_p: cobj,
    rstd_p: cobj,
    ptend_present: int,
    ptend_q_p: cobj,
    result_p: cobj,
    issue_p: cobj,
):
    wtrc_check_h2o_codon(
        ncol,
        pcols,
        pver,
        pcnst,
        pwtype,
        wtrc_nwset,
        wisotope,
        wtrc_check_total_h2o,
        wtrc_qchkmin,
        dtime,
        pstate_q_p,
        pstate_pdel_p,
        qloc_p,
        wtrc_bulk_indices_p,
        wtrc_iawset_p,
        iwspec_p,
        rstd_p,
        ptend_present,
        ptend_q_p,
        result_p,
        issue_p,
    )


@export
def wtrc_batch_clear_precip_codon(ncol: int, srfpcp_p: cobj):
    wtrc_clear_precip_codon(ncol, srfpcp_p)


@export
def wtrc_batch_diagnose_bulk_precip_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    bulk_idx: int,
    ptend_q_p: cobj,
):
    wtrc_diagnose_bulk_precip_codon(
        ncol,
        pcols,
        pver,
        top_lev,
        bulk_idx,
        ptend_q_p,
    )
