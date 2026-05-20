from math import log
import water_tracers_apply_rates_codon as _apply_rates


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

@export
def wtrc_int_eq_codon(value: int, expected: int) -> int:
    if value == expected:
        return 1
    return 0


@export
def wtrc_int_ne_codon(value: int, expected: int) -> int:
    if value != expected:
        return 1
    return 0


@export
def wtrc_bool_id_codon(value: int) -> int:
    if value != 0:
        return 1
    return 0


@export
def wtrc_select_real_codon(use_first: int, first: float, second: float) -> float:
    if use_first != 0:
        return first
    return second


@export
def wtrc_ratio_scalar_codon(qtrc: float, qtot: float, qmin: float, rstd: float) -> float:
    if abs(qtot) < qmin:
        return rstd
    return qtrc / qtot


@export
def wtrc_efac_scalar_codon(alpha: float, vapnew: float, liqnew: float, qmin: float, rstd_h2o: float) -> float:
    alov = wtrc_ratio_scalar_codon(vapnew, vapnew + liqnew, qmin, rstd_h2o)
    alov = alpha * (1.0 / alov - 1.0)
    efac = 1.0 / (alov + 1.0)
    efac = max(efac, 0.0)
    efac = min(efac, 1.0)
    return efac


@export
def wtrc_dqequil_scalar_codon(
    alpha: float,
    feq0: float,
    vtotnew: float,
    ltotnew: float,
    visoold: float,
    lisoold: float,
    qmin: float,
    rstd_h2o: float,
) -> float:
    qiso = visoold + lisoold
    vieql = qiso * wtrc_efac_scalar_codon(alpha, vtotnew, ltotnew, qmin, rstd_h2o)
    vinof = qiso * wtrc_efac_scalar_codon(1.0, vtotnew, ltotnew, qmin, rstd_h2o)
    visonew = feq0 * vieql + (1.0 - feq0) * vinof
    dviso = visonew - visoold
    if dviso < 0.0:
        dviso = max(dviso, -visoold)
    else:
        dviso = min(dviso, lisoold)
    return dviso


@export
def wtrc_liqvap_equil_scalar_codon(
    alpha: float,
    feq0: float,
    vaptot: float,
    liqtot: float,
    qmin: float,
    rstd_h2o: float,
    vapiso_p: cobj,
    liqiso_p: cobj,
    dliqiso_p: cobj,
) -> None:
    vapiso = Ptr[float](vapiso_p)
    liqiso = Ptr[float](liqiso_p)
    dliqiso = Ptr[float](dliqiso_p)

    qtiny = 1.0e-36
    dliqiso[0] = 0.0
    qtot = vaptot + liqtot
    qiso = vapiso[0] + liqiso[0]

    if qtot < qtiny:
        return
    if qiso < qtiny:
        return
    if liqtot < qtiny:
        dliqiso[0] = -liqiso[0]
        vapiso[0] = vapiso[0] - dliqiso[0]
        liqiso[0] = 0.0
        return
    if vaptot < qtiny:
        dliqiso[0] = vapiso[0]
        vapiso[0] = 0.0
        liqiso[0] = liqiso[0] + dliqiso[0]
        return

    dviso = wtrc_dqequil_scalar_codon(alpha, feq0, vaptot, liqtot, vapiso[0], liqiso[0], qmin, rstd_h2o)
    dliqiso[0] = -dviso
    liqiso[0] = liqiso[0] + dliqiso[0]
    vapiso[0] = vapiso[0] - dliqiso[0]


@export
def wtrc_init_rates_codon(pcols: int, pver: int, pwtype: int, top_lev: int, process_rates_p: cobj) -> None:
    process_rates = Ptr[float](process_rates_p)
    rtype = 1
    while rtype <= pwtype:
        isrc = 1
        while isrc <= pwtype:
            idst = 1
            while idst <= pwtype:
                k = top_lev
                while k <= pver:
                    i = 1
                    while i <= pcols:
                        idx = (
                            (i - 1)
                            + (k - 1) * pcols
                            + (idst - 1) * pcols * pver
                            + (isrc - 1) * pcols * pver * pwtype
                            + (rtype - 1) * pcols * pver * pwtype * pwtype
                        )
                        process_rates[idx] = 0.0
                        i += 1
                    k += 1
                idst += 1
            isrc += 1
        rtype += 1


@export
def wtrc_apply_rates_copy_state_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    top_lev: int,
    pstate_q_p: cobj,
    pstate_t_p: cobj,
    qloc_p: cobj,
    qloc0_p: cobj,
    tloc_p: cobj,
):
    _apply_rates.wtrc_apply_rates_copy_state_codon(
        ncol, pcols, pver, pcnst, top_lev, pstate_q_p, pstate_t_p, qloc_p, qloc0_p, tloc_p
    )


@export
def wtrc_apply_rates_zero_precip_codon(
    pcols: int,
    wtrc_nwset: int,
    rmass_p: cobj,
    smass_p: cobj,
    rmass0_p: cobj,
    smass0_p: cobj,
):
    _apply_rates.wtrc_apply_rates_zero_precip_codon(pcols, wtrc_nwset, rmass_p, smass_p, rmass0_p, smass0_p)


@export
def wtrc_apply_rates_copy_qloc0_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    top_lev: int,
    qloc_p: cobj,
    qloc0_p: cobj,
):
    _apply_rates.wtrc_apply_rates_copy_qloc0_codon(ncol, pcols, pver, pcnst, top_lev, qloc_p, qloc0_p)


@export
def wtrc_apply_rates_bulk_update_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pwtype: int,
    top_lev: int,
    dtime: float,
    bulk_indices_p: cobj,
    ptend_q_p: cobj,
    qloc_p: cobj,
):
    _apply_rates.wtrc_apply_rates_bulk_update_codon(
        ncol, pcols, pver, pwtype, top_lev, dtime, bulk_indices_p, ptend_q_p, qloc_p
    )


@export
def wtrc_apply_rates_net_tend_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    pwtype: int,
    wtrc_ncnst: int,
    top_lev: int,
    dtime: float,
    wtrc_indices_p: cobj,
    bulk_indices_p: cobj,
    pstate_q_p: cobj,
    ptend_q_p: cobj,
    qloc_p: cobj,
    diff_p: cobj,
):
    _apply_rates.wtrc_apply_rates_net_tend_codon(
        ncol,
        pcols,
        pver,
        pcnst,
        pwtype,
        wtrc_ncnst,
        top_lev,
        dtime,
        wtrc_indices_p,
        bulk_indices_p,
        pstate_q_p,
        ptend_q_p,
        qloc_p,
        diff_p,
    )


@export
def wtrc_apply_rates_first_correction_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pwtype: int,
    wtrc_nwset: int,
    top_lev: int,
    qmin: float,
    wtrc_iatype_p: cobj,
    bulk_indices_p: cobj,
    iwspec_p: cobj,
    rstd_p: cobj,
    ptend_q_p: cobj,
    diff_p: cobj,
):
    _apply_rates.wtrc_apply_rates_first_correction_codon(
        ncol,
        pcols,
        pver,
        pwtype,
        wtrc_nwset,
        top_lev,
        qmin,
        wtrc_iatype_p,
        bulk_indices_p,
        iwspec_p,
        rstd_p,
        ptend_q_p,
        diff_p,
    )


@export
def wtrc_apply_rates_second_correction_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pwtype: int,
    wtrc_nwset: int,
    top_lev: int,
    qmin: float,
    wtrc_iatype_p: cobj,
    bulk_indices_p: cobj,
    iwspec_p: cobj,
    rstd_p: cobj,
    ptend_q_p: cobj,
    diff_p: cobj,
):
    _apply_rates.wtrc_apply_rates_second_correction_codon(
        ncol,
        pcols,
        pver,
        pwtype,
        wtrc_nwset,
        top_lev,
        qmin,
        wtrc_iatype_p,
        bulk_indices_p,
        iwspec_p,
        rstd_p,
        ptend_q_p,
        diff_p,
    )


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


@export
def wtrc_q1q2_tail_shell_codon(
    lengath: int,
    pcols: int,
    pver: int,
    wtrc_nwset: int,
    msg: int,
    iwtvap: int,
    wtrc_qmin: float,
    ideep_p: cobj,
    wtrc_iatype_p: cobj,
    iwspec_p: cobj,
    rstd_p: cobj,
    jd_p: cobj,
    mx_p: cobj,
    jt_p: cobj,
    dp_p: cobj,
    dsubcld_p: cobj,
    dz_p: cobj,
    mdpc_p: cobj,
    qd_p: cobj,
    totpcp_p: cobj,
    totevp_p: cobj,
    wtcu_p: cobj,
    dupc_p: cobj,
    ql_p: cobj,
    c0mask_p: cobj,
    mupc_p: cobj,
    wtevp_p: cobj,
    rpdpc_p: cobj,
    wtrprd_p: cobj,
    cu_p: cobj,
    evp_p: cobj,
    mu_p: cobj,
    md_p: cobj,
    qu_p: cobj,
    qhat_p: cobj,
    dqdt_p: cobj,
    wtdlf_p: cobj,
    du_p: cobj,
    pevp_p: cobj,
    rprd_p: cobj,
    eps0_p: cobj,
):
    ideep = Ptr[int](ideep_p)
    wtrc_iatype = Ptr[int](wtrc_iatype_p)
    iwspec = Ptr[int](iwspec_p)
    rstd = Ptr[float](rstd_p)
    jd = Ptr[int](jd_p)
    mx = Ptr[int](mx_p)
    jt = Ptr[int](jt_p)
    dp = Ptr[float](dp_p)
    dsubcld = Ptr[float](dsubcld_p)
    dz = Ptr[float](dz_p)
    mdpc = Ptr[float](mdpc_p)
    qd = Ptr[float](qd_p)
    totpcp = Ptr[float](totpcp_p)
    totevp = Ptr[float](totevp_p)
    wtcu = Ptr[float](wtcu_p)
    dupc = Ptr[float](dupc_p)
    ql = Ptr[float](ql_p)
    c0mask = Ptr[float](c0mask_p)
    mupc = Ptr[float](mupc_p)
    wtevp = Ptr[float](wtevp_p)
    rpdpc = Ptr[float](rpdpc_p)
    wtrprd = Ptr[float](wtrprd_p)
    cu = Ptr[float](cu_p)
    evp = Ptr[float](evp_p)
    mu = Ptr[float](mu_p)
    md = Ptr[float](md_p)
    qu = Ptr[float](qu_p)
    qhat = Ptr[float](qhat_p)
    dqdt = Ptr[float](dqdt_p)
    wtdlf = Ptr[float](wtdlf_p)
    du = Ptr[float](du_p)
    pevp = Ptr[float](pevp_p)
    rprd = Ptr[float](rprd_p)
    eps0 = Ptr[float](eps0_p)

    m = 0
    while m < wtrc_nwset:
        i = 0
        while i < lengath:
            jdi = jd[i] - 1
            mxi = mx[i] - 1
            totevp[i + m * pcols] = (
                totevp[i + m * pcols]
                + mdpc[_field2_idx(i + 1, jdi + 1, pcols)] * qd[_wtrc_q1q2_3d_idx(i, jdi, m, pcols, pver)]
                - mdpc[_field2_idx(i + 1, mxi + 1, pcols)] * qd[_wtrc_q1q2_3d_idx(i, mxi, m, pcols, pver)]
            )
            i += 1
        m += 1

    m = 0
    while m < wtrc_nwset:
        vap_idx = wtrc_iatype[m + (iwtvap - 1) * wtrc_nwset] - 1
        ispec = iwspec[vap_idx]

        k = 0
        while k < pver:
            i = 0
            while i < pcols:
                pevp[_field2_idx(i + 1, k + 1, pcols)] = 0.0
                rprd[_field2_idx(i + 1, k + 1, pcols)] = 0.0
                i += 1
            k += 1

        kk = pver
        while kk >= msg + 2:
            k = kk - 1
            i = 0
            while i < lengath:
                if kk >= jt[i] and kk < mx[i] and eps0[i] > 0.0 and mupc[_field2_idx(i + 1, kk, pcols)] >= 0.0:
                    totpcp[i + m * pcols] = (
                        totpcp[i + m * pcols]
                        + dz[_field2_idx(i + 1, kk, pcols)]
                        * (
                            wtcu[_wtrc_q1q2_3d_idx(i, k, m, pcols, pver)]
                            - dupc[_field2_idx(i + 1, kk, pcols)] * ql[_wtrc_q1q2_3d_idx(i, k + 1, m, pcols, pver)]
                        )
                    )
                    rprd[_field2_idx(i + 1, kk, pcols)] = (
                        c0mask[i] * mupc[_field2_idx(i + 1, kk, pcols)] * ql[_wtrc_q1q2_3d_idx(i, k, m, pcols, pver)]
                    )
                i += 1
            kk -= 1

        kk = msg + 2
        while kk <= pver:
            k = kk - 1
            i = 0
            while i < lengath:
                ii = ideep[i] - 1
                totpcp[i + m * pcols] = max(totpcp[i + m * pcols], 0.0)
                totevp[i + m * pcols] = max(totevp[i + m * pcols], 0.0)
                if totevp[i + m * pcols] > 0.0 and totpcp[i + m * pcols] > 0.0:
                    pevp[_field2_idx(i + 1, kk, pcols)] = wtevp[_wtrc_q1q2_3d_idx(i, k, m, pcols, pver)] * min(
                        1.0, totpcp[i + m * pcols] / (totevp[i + m * pcols] + totpcp[i + m * pcols])
                    )
                else:
                    pevp[_field2_idx(i + 1, kk, pcols)] = 0.0
                rprd[_field2_idx(i + 1, kk, pcols)] = rprd[_field2_idx(i + 1, kk, pcols)] - pevp[
                    _field2_idx(i + 1, kk, pcols)
                ]
                rr = _wtrc_ratio(
                    ispec,
                    rprd[_field2_idx(i + 1, kk, pcols)],
                    rpdpc[_field2_idx(i + 1, kk, pcols)],
                    wtrc_qmin,
                    rstd,
                )
                wtrprd[_wtrc_q1q2_q_idx(ii, k, vap_idx, pcols, pver)] = (
                    rr * wtrprd[_wtrc_q1q2_q_idx(ii, k, 0, pcols, pver)]
                )
                i += 1
            kk += 1
        m += 1

    ktm = pver
    kbm = pver
    i = 0
    while i < lengath:
        ktm = min(ktm, jt[i])
        kbm = min(kbm, mx[i])
        i += 1

    m = 0
    while m < wtrc_nwset:
        vap_idx = wtrc_iatype[m + (iwtvap - 1) * wtrc_nwset] - 1
        ispec = iwspec[vap_idx]
        kk = ktm
        while kk <= pver - 1:
            k = kk - 1
            kp1 = k + 1
            i = 0
            while i < lengath:
                ii = ideep[i] - 1
                rc = _wtrc_ratio(
                    ispec,
                    wtcu[_wtrc_q1q2_3d_idx(i, k, m, pcols, pver)],
                    wtcu[_wtrc_q1q2_3d_idx(i, k, 0, pcols, pver)],
                    wtrc_qmin,
                    rstd,
                )
                re = _wtrc_ratio(
                    ispec,
                    wtevp[_wtrc_q1q2_3d_idx(i, k, m, pcols, pver)],
                    wtevp[_wtrc_q1q2_3d_idx(i, k, 0, pcols, pver)],
                    wtrc_qmin,
                    rstd,
                )
                emc = -rc * cu[_field2_idx(ii + 1, kk, pcols)] + re * evp[_field2_idx(ii + 1, kk, pcols)]
                if dp[_field2_idx(i + 1, kk, pcols)] > 0.0:
                    dqdt[_wtrc_q1q2_q_idx(ii, k, vap_idx, pcols, pver)] = emc + (
                        +mu[_field2_idx(i + 1, kk + 1, pcols)]
                        * (
                            qu[_wtrc_q1q2_3d_idx(i, kp1, m, pcols, pver)]
                            - qhat[_wtrc_q1q2_3d_idx(i, kp1, m, pcols, pver)]
                        )
                        - mu[_field2_idx(i + 1, kk, pcols)]
                        * (
                            qu[_wtrc_q1q2_3d_idx(i, k, m, pcols, pver)]
                            - qhat[_wtrc_q1q2_3d_idx(i, k, m, pcols, pver)]
                        )
                        + md[_field2_idx(i + 1, kk + 1, pcols)]
                        * (
                            qd[_wtrc_q1q2_3d_idx(i, kp1, m, pcols, pver)]
                            - qhat[_wtrc_q1q2_3d_idx(i, kp1, m, pcols, pver)]
                        )
                        - md[_field2_idx(i + 1, kk, pcols)]
                        * (
                            qd[_wtrc_q1q2_3d_idx(i, k, m, pcols, pver)]
                            - qhat[_wtrc_q1q2_3d_idx(i, k, m, pcols, pver)]
                        )
                    ) / dp[_field2_idx(i + 1, kk, pcols)]
                else:
                    dqdt[_wtrc_q1q2_q_idx(ii, k, vap_idx, pcols, pver)] = emc
                wtdlf[_wtrc_q1q2_3d_idx(ii, k, m, pcols, pver)] = du[_field2_idx(i + 1, kk, pcols)] * ql[
                    _wtrc_q1q2_3d_idx(i, kp1, m, pcols, pver)
                ]
                i += 1
            kk += 1

        kk = kbm
        while kk <= pver:
            k = kk - 1
            i = 0
            while i < lengath:
                ii = ideep[i] - 1
                if kk == mx[i]:
                    dqdt[_wtrc_q1q2_q_idx(ii, k, vap_idx, pcols, pver)] = (1.0 / dsubcld[i]) * (
                        -mu[_field2_idx(i + 1, kk, pcols)]
                        * (
                            qu[_wtrc_q1q2_3d_idx(i, k, m, pcols, pver)]
                            - qhat[_wtrc_q1q2_3d_idx(i, k, m, pcols, pver)]
                        )
                        - md[_field2_idx(i + 1, kk, pcols)]
                        * (
                            qd[_wtrc_q1q2_3d_idx(i, k, m, pcols, pver)]
                            - qhat[_wtrc_q1q2_3d_idx(i, k, m, pcols, pver)]
                        )
                    )
                elif kk > mx[i]:
                    dqdt[_wtrc_q1q2_q_idx(ii, k, vap_idx, pcols, pver)] = dqdt[
                        _wtrc_q1q2_q_idx(ii, k - 1, vap_idx, pcols, pver)
                    ]
                i += 1
            kk += 1
        m += 1


@export
def wtrc_q1q2_downdraft_shell_codon(
    lengath: int,
    pcols: int,
    pver: int,
    wtrc_nwset: int,
    msg: int,
    iwtvap: int,
    wtrc_qmin: float,
    ideep_p: cobj,
    wtrc_iatype_p: cobj,
    iwspec_p: cobj,
    rstd_p: cobj,
    jd_p: cobj,
    mx_p: cobj,
    eps0_p: cobj,
    qu_p: cobj,
    qds_p: cobj,
    rd_p: cobj,
    qd_p: cobj,
    dz_p: cobj,
    wtevp_p: cobj,
    evpc_p: cobj,
    ed_p: cobj,
    q_p: cobj,
    mdpc_p: cobj,
    qdb_p: cobj,
    totevp_p: cobj,
):
    ideep = Ptr[int](ideep_p)
    wtrc_iatype = Ptr[int](wtrc_iatype_p)
    iwspec = Ptr[int](iwspec_p)
    rstd = Ptr[float](rstd_p)
    jd = Ptr[int](jd_p)
    mx = Ptr[int](mx_p)
    eps0 = Ptr[float](eps0_p)
    qu = Ptr[float](qu_p)
    qds = Ptr[float](qds_p)
    rd = Ptr[float](rd_p)
    qd = Ptr[float](qd_p)
    dz = Ptr[float](dz_p)
    wtevp = Ptr[float](wtevp_p)
    evpc = Ptr[float](evpc_p)
    ed = Ptr[float](ed_p)
    q = Ptr[float](q_p)
    mdpc = Ptr[float](mdpc_p)
    qdb = Ptr[float](qdb_p)
    totevp = Ptr[float](totevp_p)

    i = 0
    while i < lengath:
        jdi = jd[i] - 1
        m = 0
        while m < wtrc_nwset:
            vap_idx = wtrc_iatype[m + (iwtvap - 1) * wtrc_nwset] - 1
            ispec = iwspec[vap_idx]
            rd[_wtrc_q1q2_3d_idx(i, jdi, m, pcols, pver)] = _wtrc_ratio(
                ispec,
                qu[_wtrc_q1q2_3d_idx(i, jdi, m, pcols, pver)],
                qu[_wtrc_q1q2_3d_idx(i, jdi, 0, pcols, pver)],
                wtrc_qmin,
                rstd,
            )
            qd[_wtrc_q1q2_3d_idx(i, jdi, m, pcols, pver)] = (
                rd[_wtrc_q1q2_3d_idx(i, jdi, m, pcols, pver)] * qds[_field2_idx(i + 1, jdi + 1, pcols)]
            )
            m += 1
        i += 1

    kk = msg + 2
    while kk <= pver:
        k = kk - 1
        kp1 = k + 1
        i = 0
        while i < lengath:
            if kk >= jd[i] and kk < mx[i] and eps0[i] > 0.0:
                oval = 0.0
                dqdiff = 0.0
                ii = ideep[i] - 1
                m = 0
                while m < wtrc_nwset:
                    vap_idx = wtrc_iatype[m + (iwtvap - 1) * wtrc_nwset] - 1
                    ispec = iwspec[vap_idx]
                    rd[_wtrc_q1q2_3d_idx(i, k, m, pcols, pver)] = _wtrc_ratio(
                        ispec,
                        qd[_wtrc_q1q2_3d_idx(i, k, m, pcols, pver)],
                        qd[_wtrc_q1q2_3d_idx(i, k, 0, pcols, pver)],
                        wtrc_qmin,
                        rstd,
                    )
                    wtevp[_wtrc_q1q2_3d_idx(i, k, m, pcols, pver)] = (
                        rd[_wtrc_q1q2_3d_idx(i, k, m, pcols, pver)] * evpc[_field2_idx(i + 1, kk, pcols)]
                    )
                    qd[_wtrc_q1q2_3d_idx(i, kp1, m, pcols, pver)] = (
                        dz[_field2_idx(i + 1, kk, pcols)]
                        * (
                            wtevp[_wtrc_q1q2_3d_idx(i, k, m, pcols, pver)]
                            + ed[_field2_idx(i + 1, kk, pcols)] * q[_wtrc_q1q2_q_idx(ii, k, vap_idx, pcols, pver)]
                        )
                        - mdpc[_field2_idx(i + 1, kk, pcols)] * qd[_wtrc_q1q2_3d_idx(i, k, m, pcols, pver)]
                    ) / (-mdpc[_field2_idx(i + 1, kk + 1, pcols)])
                    if m == 0:
                        oval = qd[_wtrc_q1q2_3d_idx(i, kp1, m, pcols, pver)]
                        dqdiff = qdb[_field2_idx(i + 1, kk + 1, pcols)] - qd[
                            _wtrc_q1q2_3d_idx(i, kp1, m, pcols, pver)
                        ]
                    rfix = _wtrc_ratio(
                        ispec,
                        qd[_wtrc_q1q2_3d_idx(i, kp1, m, pcols, pver)],
                        oval,
                        wtrc_qmin,
                        rstd,
                    )
                    qd[_wtrc_q1q2_3d_idx(i, kp1, m, pcols, pver)] = (
                        qd[_wtrc_q1q2_3d_idx(i, kp1, m, pcols, pver)] + rfix * dqdiff
                    )
                    totevp[i + m * pcols] = (
                        totevp[i + m * pcols]
                        - dz[_field2_idx(i + 1, kk, pcols)]
                        * ed[_field2_idx(i + 1, kk, pcols)]
                        * q[_wtrc_q1q2_q_idx(ii, k, vap_idx, pcols, pver)]
                    )
                    m += 1
            i += 1
        kk += 1


@export
def wtrc_q1q2_downdraft_tail_shell_codon(
    lengath: int,
    pcols: int,
    pver: int,
    wtrc_nwset: int,
    msg: int,
    iwtvap: int,
    wtrc_qmin: float,
    ideep_p: cobj,
    wtrc_iatype_p: cobj,
    iwspec_p: cobj,
    rstd_p: cobj,
    jd_p: cobj,
    mx_p: cobj,
    jt_p: cobj,
    eps0_p: cobj,
    qu_p: cobj,
    qds_p: cobj,
    rd_p: cobj,
    qd_p: cobj,
    dz_p: cobj,
    wtevp_p: cobj,
    evpc_p: cobj,
    ed_p: cobj,
    q_p: cobj,
    mdpc_p: cobj,
    qdb_p: cobj,
    totevp_p: cobj,
    dp_p: cobj,
    dsubcld_p: cobj,
    totpcp_p: cobj,
    wtcu_p: cobj,
    dupc_p: cobj,
    ql_p: cobj,
    c0mask_p: cobj,
    mupc_p: cobj,
    rpdpc_p: cobj,
    wtrprd_p: cobj,
    cu_p: cobj,
    evp_p: cobj,
    mu_p: cobj,
    md_p: cobj,
    qhat_p: cobj,
    dqdt_p: cobj,
    wtdlf_p: cobj,
    du_p: cobj,
    pevp_p: cobj,
    rprd_p: cobj,
):
    wtrc_q1q2_downdraft_shell_codon(
        lengath,
        pcols,
        pver,
        wtrc_nwset,
        msg,
        iwtvap,
        wtrc_qmin,
        ideep_p,
        wtrc_iatype_p,
        iwspec_p,
        rstd_p,
        jd_p,
        mx_p,
        eps0_p,
        qu_p,
        qds_p,
        rd_p,
        qd_p,
        dz_p,
        wtevp_p,
        evpc_p,
        ed_p,
        q_p,
        mdpc_p,
        qdb_p,
        totevp_p,
    )
    wtrc_q1q2_tail_shell_codon(
        lengath,
        pcols,
        pver,
        wtrc_nwset,
        msg,
        iwtvap,
        wtrc_qmin,
        ideep_p,
        wtrc_iatype_p,
        iwspec_p,
        rstd_p,
        jd_p,
        mx_p,
        jt_p,
        dp_p,
        dsubcld_p,
        dz_p,
        mdpc_p,
        qd_p,
        totpcp_p,
        totevp_p,
        wtcu_p,
        dupc_p,
        ql_p,
        c0mask_p,
        mupc_p,
        wtevp_p,
        rpdpc_p,
        wtrprd_p,
        cu_p,
        evp_p,
        mu_p,
        md_p,
        qu_p,
        qhat_p,
        dqdt_p,
        wtdlf_p,
        du_p,
        pevp_p,
        rprd_p,
        eps0_p,
    )


@export
def wtrc_q1q2_updraft_h2o_shell_codon(
    lengath: int,
    pcols: int,
    pver: int,
    wtrc_nwset: int,
    msg: int,
    iwtvap: int,
    wtrc_qmin: float,
    ideep_p: cobj,
    wtrc_iatype_p: cobj,
    iwspec_p: cobj,
    rstd_p: cobj,
    mx_p: cobj,
    jt_p: cobj,
    eps0_p: cobj,
    tu_p: cobj,
    mupc_p: cobj,
    qu_p: cobj,
    ru_p: cobj,
    wtcu_p: cobj,
    cupc_p: cobj,
    dz_p: cobj,
    eu_p: cobj,
    q_p: cobj,
    dupc_p: cobj,
    qstb_p: cobj,
    qub_p: cobj,
    ql_p: cobj,
    c0mask_p: cobj,
    oval_work_p: cobj,
    uqdiff_work_p: cobj,
):
    ideep = Ptr[int](ideep_p)
    wtrc_iatype = Ptr[int](wtrc_iatype_p)
    iwspec = Ptr[int](iwspec_p)
    rstd = Ptr[float](rstd_p)
    mx = Ptr[int](mx_p)
    jt = Ptr[int](jt_p)
    eps0 = Ptr[float](eps0_p)
    tu = Ptr[float](tu_p)
    mupc = Ptr[float](mupc_p)
    qu = Ptr[float](qu_p)
    ru = Ptr[float](ru_p)
    wtcu = Ptr[float](wtcu_p)
    cupc = Ptr[float](cupc_p)
    dz = Ptr[float](dz_p)
    eu = Ptr[float](eu_p)
    q = Ptr[float](q_p)
    dupc = Ptr[float](dupc_p)
    qstb = Ptr[float](qstb_p)
    qub = Ptr[float](qub_p)
    ql = Ptr[float](ql_p)
    c0mask = Ptr[float](c0mask_p)
    oval_work = Ptr[float](oval_work_p)
    uqdiff_work = Ptr[float](uqdiff_work_p)

    k = 0
    while k < pver:
        i = 0
        while i < pcols:
            oval_work[_field2_idx(i + 1, k + 1, pcols)] = 0.0
            uqdiff_work[_field2_idx(i + 1, k + 1, pcols)] = 0.0
            i += 1
        k += 1

    kk = pver
    while kk >= msg + 2:
        k = kk - 1
        kp1 = k + 1
        i = 0
        while i < lengath:
            ii = ideep[i] - 1
            if kk > 1 and kk < mx[i] and eps0[i] > 0.0:
                if mupc[_field2_idx(i + 1, kk, pcols)] > 0.0:
                    vap_idx = wtrc_iatype[iwtvap * wtrc_nwset - wtrc_nwset] - 1
                    ispec = iwspec[vap_idx]
                    ru[_wtrc_q1q2_3d_idx(i, k, 0, pcols, pver)] = _wtrc_ratio(
                        ispec,
                        qu[_wtrc_q1q2_3d_idx(i, kp1, 0, pcols, pver)],
                        qu[_wtrc_q1q2_3d_idx(i, kp1, 0, pcols, pver)],
                        wtrc_qmin,
                        rstd,
                    )
                    wtcu[_wtrc_q1q2_3d_idx(i, k, 0, pcols, pver)] = (
                        ru[_wtrc_q1q2_3d_idx(i, k, 0, pcols, pver)] * cupc[_field2_idx(i + 1, kk, pcols)]
                    )
                    qu[_wtrc_q1q2_3d_idx(i, k, 0, pcols, pver)] = (
                        mupc[_field2_idx(i + 1, kk + 1, pcols)]
                        / mupc[_field2_idx(i + 1, kk, pcols)]
                        * qu[_wtrc_q1q2_3d_idx(i, kp1, 0, pcols, pver)]
                        + dz[_field2_idx(i + 1, kk, pcols)]
                        / mupc[_field2_idx(i + 1, kk, pcols)]
                        * (
                            eu[_field2_idx(i + 1, kk, pcols)] * q[_wtrc_q1q2_q_idx(ii, k, vap_idx, pcols, pver)]
                            - dupc[_field2_idx(i + 1, kk, pcols)]
                            * ru[_wtrc_q1q2_3d_idx(i, k, 0, pcols, pver)]
                            * qstb[_field2_idx(i + 1, kk, pcols)]
                            - wtcu[_wtrc_q1q2_3d_idx(i, k, 0, pcols, pver)]
                        )
                    )
                    oval_work[_field2_idx(i + 1, kk, pcols)] = qu[_wtrc_q1q2_3d_idx(i, k, 0, pcols, pver)]
                    uqdiff_work[_field2_idx(i + 1, kk, pcols)] = qub[_field2_idx(i + 1, kk, pcols)] - qu[
                        _wtrc_q1q2_3d_idx(i, k, 0, pcols, pver)
                    ]
                    rfix = _wtrc_ratio(
                        ispec,
                        qu[_wtrc_q1q2_3d_idx(i, k, 0, pcols, pver)],
                        oval_work[_field2_idx(i + 1, kk, pcols)],
                        wtrc_qmin,
                        rstd,
                    )
                    qu[_wtrc_q1q2_3d_idx(i, k, 0, pcols, pver)] = (
                        qu[_wtrc_q1q2_3d_idx(i, k, 0, pcols, pver)]
                        + rfix * uqdiff_work[_field2_idx(i + 1, kk, pcols)]
                    )
                    if kk >= jt[i]:
                        ql1 = (
                            1.0
                            / mupc[_field2_idx(i + 1, kk, pcols)]
                            * (
                                mupc[_field2_idx(i + 1, kk + 1, pcols)]
                                * ql[_wtrc_q1q2_3d_idx(i, kp1, 0, pcols, pver)]
                                - dz[_field2_idx(i + 1, kk, pcols)]
                                * dupc[_field2_idx(i + 1, kk, pcols)]
                                * ql[_wtrc_q1q2_3d_idx(i, kp1, 0, pcols, pver)]
                                + dz[_field2_idx(i + 1, kk, pcols)]
                                * wtcu[_wtrc_q1q2_3d_idx(i, k, 0, pcols, pver)]
                            )
                        )
                        ql[_wtrc_q1q2_3d_idx(i, k, 0, pcols, pver)] = ql1 / (
                            1.0 + dz[_field2_idx(i + 1, kk, pcols)] * c0mask[i]
                        )
                    else:
                        ql[_wtrc_q1q2_3d_idx(i, k, 0, pcols, pver)] = 0.0
                else:
                    m = 0
                    while m < wtrc_nwset:
                        vap_idx = wtrc_iatype[m + (iwtvap - 1) * wtrc_nwset] - 1
                        qu[_wtrc_q1q2_3d_idx(i, k, m, pcols, pver)] = q[
                            _wtrc_q1q2_q_idx(ii, k, vap_idx, pcols, pver)
                        ]
                        m += 1
            i += 1
        kk -= 1


@export
def wtrc_q1q2_stage_dispatch_codon(
    stage: int,
    lengath: int,
    pcols: int,
    pver: int,
    pcnst: int,
    pwtype: int,
    wtrc_nwset: int,
    msg: int,
    iwtvap: int,
    wtrc_qmin: float,
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
):
    if stage == 1:
        wtrc_q1q2_init_qhat_shell_codon(
            lengath,
            pcols,
            pver,
            pcnst,
            pwtype,
            wtrc_nwset,
            msg,
            iwtvap,
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
        )
    elif stage == 2:
        wtrc_q1q2_updraft_h2o_shell_codon(
            lengath,
            pcols,
            pver,
            wtrc_nwset,
            msg,
            iwtvap,
            wtrc_qmin,
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
        )
    elif stage == 3:
        wtrc_q1q2_downdraft_shell_codon(
            lengath,
            pcols,
            pver,
            wtrc_nwset,
            msg,
            iwtvap,
            wtrc_qmin,
            p1,
            p2,
            p3,
            p4,
            p5,
            p6,
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
        )
        wtrc_q1q2_tail_shell_codon(
            lengath,
            pcols,
            pver,
            wtrc_nwset,
            msg,
            iwtvap,
            wtrc_qmin,
            p1,
            p2,
            p3,
            p4,
            p5,
            p6,
            p7,
            p21,
            p22,
            p13,
            p18,
            p12,
            p23,
            p20,
            p24,
            p25,
            p26,
            p27,
            p28,
            p14,
            p29,
            p30,
            p31,
            p32,
            p33,
            p34,
            p9,
            p35,
            p36,
            p37,
            p38,
            p39,
            p40,
            p8,
        )
    elif stage == 4:
        wtrc_q1q2_downdraft_shell_codon(
            lengath,
            pcols,
            pver,
            wtrc_nwset,
            msg,
            iwtvap,
            wtrc_qmin,
            p1,
            p2,
            p3,
            p4,
            p5,
            p6,
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
        )
    elif stage == 5:
        wtrc_q1q2_tail_shell_codon(
            lengath,
            pcols,
            pver,
            wtrc_nwset,
            msg,
            iwtvap,
            wtrc_qmin,
            p1,
            p2,
            p3,
            p4,
            p5,
            p6,
            p7,
            p21,
            p22,
            p13,
            p18,
            p12,
            p23,
            p20,
            p24,
            p25,
            p26,
            p27,
            p28,
            p14,
            p29,
            p30,
            p31,
            p32,
            p33,
            p34,
            p9,
            p35,
            p36,
            p37,
            p38,
            p39,
            p40,
            p8,
        )


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
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    state_q_p: cobj,
    qst_p: cobj,
    rh_p: cobj,
    state_zi_p: cobj,
    dz_p: cobj,
    evpbulk_p: cobj,
    subbulk_p: cobj,
    rnbulk_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    qst = Ptr[float](qst_p)
    rh = Ptr[float](rh_p)
    state_zi = Ptr[float](state_zi_p)
    dz = Ptr[float](dz_p)
    evpbulk = Ptr[float](evpbulk_p)
    subbulk = Ptr[float](subbulk_p)
    rnbulk = Ptr[float](rnbulk_p)

    k = 0
    while k < pver:
        i = 0
        while i < ncol:
            idx = i + k * pcols
            rh[idx] = state_q[idx] / qst[idx]
            dz[idx] = state_zi[i + k * pcols] - state_zi[i + (k + 1) * pcols]
            rnbulk[idx] = evpbulk[idx] - subbulk[idx]
            i += 1
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


@export
def wtrc_stage_dispatch_codon(
    stage: int,
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    pcnst: int,
    pwtype: int,
    wtrc_nwset: int,
    i1: int,
    i2: int,
    i3: int,
    i4: int,
    i5: int,
    r1: float,
    r2: float,
    r3: float,
    r4: float,
    r5: float,
    r6: float,
    r7: float,
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
):
    if stage == 1:
        wtrc_batch_mass_fixer_codon(
            ncol,
            pcols,
            pver,
            pcnst,
            pwtype,
            wtrc_nwset,
            i1,
            i2,
            r1,
            r2,
            r3,
            r4,
            r5,
            r6,
            r7,
            p1,
            p2,
            p3,
            p4,
            p5,
            p6,
            p7,
        )
    elif stage == 2:
        wtrc_batch_check_h2o_codon(
            ncol,
            pcols,
            pver,
            pcnst,
            pwtype,
            wtrc_nwset,
            i1,
            i2,
            r1,
            r2,
            p1,
            p2,
            p3,
            p4,
            p5,
            p6,
            p7,
            i3,
            p8,
            p9,
            p10,
        )
    elif stage == 3:
        wtrc_batch_clear_precip_codon(ncol, p1)
    elif stage == 4:
        wtrc_batch_diagnose_bulk_precip_codon(ncol, pcols, pver, i1, i2, p1)
    elif stage == 5:
        wtrc_precip_evap_init_shell_codon(
            pcols,
            pver,
            pverp,
            i1,
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
        )
    elif stage == 6:
        wtrc_precip_evap_prep_shell_codon(
            ncol,
            pcols,
            pver,
            pverp,
            p1,
            p2,
            p3,
            p4,
            p5,
            p6,
            p7,
            p8,
        )
    elif stage == 7:
        wtrc_precip_evap_tail_shell_codon(
            ncol,
            pcols,
            pver,
            pverp,
            pcnst,
            i1,
            wtrc_nwset,
            i2,
            r1,
            p1,
            p2,
            p3,
            p4,
            p5,
            p6,
            p7,
        )
