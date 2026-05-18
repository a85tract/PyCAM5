from math import exp, log10


def _idx2(i: int, k: int, pcols: int):
    return (k - 1) * pcols + (i - 1)


def _idx3(i: int, k: int, m: int, pcols: int, pver: int):
    return (m - 1) * pcols * pver + (k - 1) * pcols + (i - 1)


@export
def cldwat2m_rhcrit_const_codon(
    pcols: int,
    pver: int,
    rhmini_const: float,
    rhminl_const: float,
    rhminl_adj_land_const: float,
    rhminh_const: float,
    rhmini_p: cobj,
    rhminl_p: cobj,
    rhminl_adj_land_p: cobj,
    rhminh_p: cobj,
):
    rhmini = Ptr[float](rhmini_p)
    rhminl = Ptr[float](rhminl_p)
    rhminl_adj_land = Ptr[float](rhminl_adj_land_p)
    rhminh = Ptr[float](rhminh_p)

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            idx = _idx2(i, k, pcols)
            rhmini[idx] = rhmini_const
            rhminl[idx] = rhminl_const
            rhminl_adj_land[idx] = rhminl_adj_land_const
            rhminh[idx] = rhminh_const


@export
def cldwat2m_positive_moisture_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    do_cldice: int,
    dt: float,
    latvap: float,
    latice: float,
    cpair: float,
    dp_p: cobj,
    qvmin_p: cobj,
    qlmin_p: cobj,
    qimin_p: cobj,
    qv_p: cobj,
    ql_p: cobj,
    qi_p: cobj,
    t_p: cobj,
    qvten_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    tten_p: cobj,
):
    dp = Ptr[float](dp_p)
    qvmin = Ptr[float](qvmin_p)
    qlmin = Ptr[float](qlmin_p)
    qimin = Ptr[float](qimin_p)
    qv = Ptr[float](qv_p)
    ql = Ptr[float](ql_p)
    qi = Ptr[float](qi_p)
    t = Ptr[float](t_p)
    qvten = Ptr[float](qvten_p)
    qlten = Ptr[float](qlten_p)
    qiten = Ptr[float](qiten_p)
    tten = Ptr[float](tten_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            tten[idx] = 0.0
            qvten[idx] = 0.0
            qlten[idx] = 0.0
            qiten[idx] = 0.0

    for i in range(1, ncol + 1):
        needs_fix = False
        for k in range(top_lev, pver + 1):
            idx = _idx2(i, k, pcols)
            if qv[idx] < qvmin[idx] or ql[idx] < qlmin[idx] or qi[idx] < qimin[idx]:
                needs_fix = True
                break
        if not needs_fix:
            continue

        dqv = 0.0
        for k in range(top_lev, pver + 1):
            idx = _idx2(i, k, pcols)

            dql = qlmin[idx] - ql[idx]
            if dql < 0.0:
                dql = 0.0

            if do_cldice != 0:
                dqi = qimin[idx] - qi[idx]
                if dqi < 0.0:
                    dqi = 0.0
            else:
                dqi = 0.0

            qlten[idx] = qlten[idx] + dql / dt
            qiten[idx] = qiten[idx] + dqi / dt
            qvten[idx] = qvten[idx] - (dql + dqi) / dt
            tten[idx] = tten[idx] + (latvap / cpair) * (dql / dt) + ((latvap + latice) / cpair) * (dqi / dt)
            ql[idx] = ql[idx] + dql
            qi[idx] = qi[idx] + dqi
            qv[idx] = qv[idx] - dql - dqi
            t[idx] = t[idx] + (latvap * dql + (latvap + latice) * dqi) / cpair

            dqv = qvmin[idx] - qv[idx]
            if dqv < 0.0:
                dqv = 0.0
            qvten[idx] = qvten[idx] + dqv / dt
            qv[idx] = qv[idx] + dqv
            if k != pver:
                idx_next = _idx2(i, k + 1, pcols)
                transfer = dqv * dp[idx] / dp[idx_next]
                qv[idx_next] = qv[idx_next] - transfer
                qvten[idx_next] = qvten[idx_next] - transfer / dt

            if qv[idx] < qvmin[idx]:
                qv[idx] = qvmin[idx]
            if ql[idx] < qlmin[idx]:
                ql[idx] = qlmin[idx]
            if qi[idx] < qimin[idx]:
                qi[idx] = qimin[idx]

        if dqv > 1.0e-20:
            sum_val = 0.0
            for k in range(top_lev, pver + 1):
                idx = _idx2(i, k, pcols)
                if qv[idx] > 2.0 * qvmin[idx]:
                    sum_val = sum_val + qv[idx] * dp[idx]

            denom = sum_val
            if denom < 1.0e-20:
                denom = 1.0e-20
            aa = dqv * dp[_idx2(i, pver, pcols)] / denom
            if aa < 0.5:
                for k in range(top_lev, pver + 1):
                    idx = _idx2(i, k, pcols)
                    if qv[idx] > 2.0 * qvmin[idx]:
                        dum = aa * qv[idx]
                        qv[idx] = qv[idx] - dum
                        qvten[idx] = qvten[idx] - dum / dt


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
def _iawset_idx(itype: int, iwset: int, pwtype: int):
    return (itype - 1) + (iwset - 1) * pwtype


@inline
def _iatype_idx(iwset: int, itype: int, wtrc_nwset: int):
    return (iwset - 1) + (itype - 1) * wtrc_nwset


@inline
def _spec_idx(ispec: int):
    return ispec - 1


@inline
def _bulk_idx(itype: int):
    return itype - 1


@inline
def _wtrc_ratio(ispec: int, qtrc: float, qtot: float, wtrc_qmin: float, rstd) -> float:
    if abs(qtot) < wtrc_qmin:
        return rstd[_spec_idx(ispec)]
    return qtrc / qtot


@inline
def _wiso_alpl(ispec: int, tk: float) -> float:
    if ispec <= 2:
        return 1.0
    if ispec == 3:
        return exp(
            1158.8e-12 * tk**3
            + (-1620.1e-9) * tk**2
            + 794.84e-6 * tk
            + (-161.04e-3)
            + 2.9992e6 / tk**3
        )
    return exp(0.35041e6 / tk**3 + (-1.6664e3) / tk**2 + 6.7123 / tk + (-7.685e-3))


@inline
def _wiso_alpi(ispec: int, tk: float) -> float:
    if ispec <= 2:
        return 1.0
    if ispec == 3:
        return exp(16289.0 / tk**2 + (-9.45e-2))
    return exp(11.839 / tk + (-28.224e-3))


@inline
def _wiso_ssatf(tk: float) -> float:
    ssat = 1.0 + (-0.002) * (tk - 273.16)
    if ssat < 1.0:
        ssat = 1.0
    if ssat > 2.0:
        ssat = 2.0
    return ssat


@inline
def _wiso_akci(ispec: int, tk: float, alpeq: float) -> float:
    if tk >= 253.15:
        return alpeq
    sat1 = _wiso_ssatf(tk)
    difrmj = 1.0
    if ispec == 3:
        difrmj = 0.9757
    elif ispec == 4:
        difrmj = 0.9727
    dondi = 1.0 / difrmj
    return alpeq * sat1 / (alpeq * dondi * (sat1 - 1.0) + 1.0)


@inline
def _qsat_water(t: float, p: float, epsilo: float) -> float:
    tboil = 373.16
    es = 10.0 ** (
        -7.90298 * (tboil / t - 1.0)
        + 5.02808 * log10(tboil / t)
        - 1.3816e-7 * (10.0 ** (11.344 * (1.0 - t / tboil)) - 1.0)
        + 8.1328e-3 * (10.0 ** (-3.49149 * (tboil / t - 1.0)) - 1.0)
        + log10(1013.246)
    ) * 100.0
    if (p - es) <= 0.0:
        return 1.0
    return epsilo * es / (p - (1.0 - epsilo) * es)


@inline
def _wtrc_get_alpha(
    q: float,
    tk: float,
    ispec: int,
    isrctype: int,
    idsttype: int,
    rhclc: int,
    porqh: float,
    kin: int,
    wisotope: int,
    iwtvap: int,
    iwtliq: int,
    epsilo: float,
) -> float:
    if wisotope == 0:
        return 1.0

    rh = porqh
    if rhclc != 0:
        rh = q / _qsat_water(tk, porqh, epsilo)

    alpha = 1.0
    if isrctype != idsttype:
        if isrctype == iwtvap:
            if idsttype == iwtliq:
                alpha = _wiso_alpl(ispec, tk)
            else:
                alpha = _wiso_alpi(ispec, tk)
                if kin != 0:
                    alpha = _wiso_akci(ispec, tk, alpha)
        elif idsttype == iwtvap:
            if isrctype == iwtliq:
                alpha = _wiso_alpl(ispec, tk)
                alpha = 1.0 / alpha
            else:
                alpha = 1.0
    return alpha


@inline
def _wtrc_efac(alpha: float, vapnew: float, liqnew: float, wtrc_qmin: float, rstd) -> float:
    alov = _wtrc_ratio(1, vapnew, vapnew + liqnew, wtrc_qmin, rstd)
    alov = alpha * (1.0 / alov - 1.0)
    efac = 1.0 / (alov + 1.0)
    if efac < 0.0:
        efac = 0.0
    if efac > 1.0:
        efac = 1.0
    return efac


@inline
def _wtrc_dqequil(
    alpha: float,
    feq0: float,
    vtotnew: float,
    ltotnew: float,
    visoold: float,
    lisoold: float,
    wtrc_qmin: float,
    rstd,
) -> float:
    qiso = visoold + lisoold
    vieql = qiso * _wtrc_efac(alpha, vtotnew, ltotnew, wtrc_qmin, rstd)
    vinof = qiso * _wtrc_efac(1.0, vtotnew, ltotnew, wtrc_qmin, rstd)
    visonew = feq0 * vieql + (1.0 - feq0) * vinof
    dviso = visonew - visoold
    if dviso < 0.0:
        if dviso < (-visoold):
            dviso = -visoold
    else:
        if dviso > lisoold:
            dviso = lisoold
    return dviso


@inline
def _wtrc_liqvap_equil(
    alpha: float,
    feq0: float,
    vaptot: float,
    liqtot: float,
    vapiso: float,
    liqiso: float,
    wtrc_qmin: float,
    rstd,
):
    qtiny = 1.0e-36
    qtot = vaptot + liqtot
    qiso = vapiso + liqiso

    if qtot < qtiny or qiso < qtiny:
        return vapiso, liqiso

    if liqtot < qtiny:
        dliqiso = -liqiso
        vapiso = vapiso - dliqiso
        liqiso = 0.0
        return vapiso, liqiso

    if vaptot < qtiny:
        dliqiso = vapiso
        vapiso = 0.0
        liqiso = liqiso + dliqiso
        return vapiso, liqiso

    dviso = _wtrc_dqequil(alpha, feq0, vaptot, liqtot, vapiso, liqiso, wtrc_qmin, rstd)
    dliqiso = -dviso
    liqiso = liqiso + dliqiso
    vapiso = vapiso - dliqiso
    return vapiso, liqiso


@export
def macrop_driver_wtrc_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    pwtype: int,
    top_lev: int,
    wtrc_niter: int,
    wtrc_ncnst: int,
    wtrc_nwset: int,
    wisotope: int,
    iwtvap: int,
    iwtliq: int,
    iwtice: int,
    cpair: float,
    dtime: float,
    wtrc_qmin: float,
    epsilo: float,
    state_q_p: cobj,
    state_t_p: cobj,
    state_pmid_p: cobj,
    ptend_q_p: cobj,
    qvlat_p: cobj,
    qcten_p: cobj,
    qiten_p: cobj,
    prelat_p: cobj,
    process_rates_p: cobj,
    qloc_p: cobj,
    qloc0_p: cobj,
    tloc_p: cobj,
    diff_p: cobj,
    wtrc_iawset_p: cobj,
    wtrc_iatype_p: cobj,
    wtrc_bulk_indices_p: cobj,
    wtrc_indices_p: cobj,
    iwspec_p: cobj,
    rstd_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    state_t = Ptr[float](state_t_p)
    state_pmid = Ptr[float](state_pmid_p)
    ptend_q = Ptr[float](ptend_q_p)
    qvlat = Ptr[float](qvlat_p)
    qcten = Ptr[float](qcten_p)
    qiten = Ptr[float](qiten_p)
    prelat = Ptr[float](prelat_p)
    process_rates = Ptr[float](process_rates_p)
    qloc = Ptr[float](qloc_p)
    qloc0 = Ptr[float](qloc0_p)
    tloc = Ptr[float](tloc_p)
    diff = Ptr[float](diff_p)
    wtrc_iawset = Ptr[int](wtrc_iawset_p)
    wtrc_iatype = Ptr[int](wtrc_iatype_p)
    wtrc_bulk_indices = Ptr[int](wtrc_bulk_indices_p)
    wtrc_indices = Ptr[int](wtrc_indices_p)
    iwspec = Ptr[int](iwspec_p)
    rstd = Ptr[float](rstd_p)

    for rtype in range(1, pwtype + 1):
        for isrctype in range(1, pwtype + 1):
            for idsttype in range(1, pwtype + 1):
                for k in range(top_lev, pver + 1):
                    for i in range(1, pcols + 1):
                        process_rates[
                            _process_rates_idx(
                                i,
                                k,
                                idsttype,
                                isrctype,
                                rtype,
                                pcols,
                                pver,
                                pwtype,
                            )
                        ] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)

            rate_val = qvlat[idx2] + qcten[idx2]
            rate_val = rate_val + qiten[idx2]
            process_rates[
                _process_rates_idx(
                    i, k, iwtvap, iwtvap, iwtvap, pcols, pver, pwtype
                )
            ] = process_rates[
                _process_rates_idx(
                    i, k, iwtvap, iwtvap, iwtvap, pcols, pver, pwtype
                )
            ] + rate_val

            rate_val = qcten[idx2]
            if rate_val < 0.0:
                dst_idx = _process_rates_idx(
                    i, k, iwtliq, iwtvap, iwtliq, pcols, pver, pwtype
                )
                src_idx = _process_rates_idx(
                    i, k, iwtvap, iwtliq, iwtliq, pcols, pver, pwtype
                )
            else:
                dst_idx = _process_rates_idx(
                    i, k, iwtliq, iwtvap, iwtvap, pcols, pver, pwtype
                )
                src_idx = _process_rates_idx(
                    i, k, iwtvap, iwtliq, iwtvap, pcols, pver, pwtype
                )
            process_rates[dst_idx] = process_rates[dst_idx] + rate_val
            process_rates[src_idx] = process_rates[src_idx] - rate_val

            rate_val = qiten[idx2]
            if rate_val < 0.0:
                dst_idx = _process_rates_idx(
                    i, k, iwtice, iwtvap, iwtice, pcols, pver, pwtype
                )
                src_idx = _process_rates_idx(
                    i, k, iwtvap, iwtice, iwtice, pcols, pver, pwtype
                )
            else:
                dst_idx = _process_rates_idx(
                    i, k, iwtice, iwtvap, iwtvap, pcols, pver, pwtype
                )
                src_idx = _process_rates_idx(
                    i, k, iwtvap, iwtice, iwtvap, pcols, pver, pwtype
                )
            process_rates[dst_idx] = process_rates[dst_idx] + rate_val
            process_rates[src_idx] = process_rates[src_idx] - rate_val

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            tloc[idx2] = state_t[idx2]
            for icnst in range(1, wtrc_ncnst + 1):
                trc_idx = wtrc_indices[icnst - 1]
                idx3 = _idx3(i, k, trc_idx, pcols, pver)
                qloc[idx3] = state_q[idx3]
                qloc0[idx3] = qloc[idx3]

    for iter_idx in range(1, wtrc_niter + 1):
        for k in range(top_lev, pver + 1):
            for i in range(1, ncol + 1):
                idx2 = _idx2(i, k, pcols)
                tloc[idx2] = (tloc[idx2] + (tloc[idx2] + prelat[idx2] / cpair * dtime)) / 2.0

                for isrctype in range(1, pwtype + 1):
                    for idsttype in range(1, pwtype + 1):
                        rtype = isrctype
                        rate_val = process_rates[
                            _process_rates_idx(
                                i,
                                k,
                                idsttype,
                                isrctype,
                                rtype,
                                pcols,
                                pver,
                                pwtype,
                            )
                        ]
                        if rate_val > 0.0:
                            for iwset in range(1, wtrc_nwset + 1):
                                msrc = wtrc_iawset[_iawset_idx(isrctype, iwset, pwtype)]
                                mbase = wtrc_iawset[_iawset_idx(isrctype, 1, pwtype)]
                                mdst = wtrc_iawset[_iawset_idx(idsttype, iwset, pwtype)]

                                idx_msrc = _idx3(i, k, msrc, pcols, pver)
                                idx_mbase = _idx3(i, k, mbase, pcols, pver)
                                idx_mdst = _idx3(i, k, mdst, pcols, pver)

                                R = _wtrc_ratio(
                                    iwspec[msrc - 1],
                                    qloc0[idx_msrc],
                                    qloc0[idx_mbase],
                                    wtrc_qmin,
                                    rstd,
                                )

                                if (
                                    wisotope != 0
                                    and iwset != 1
                                    and isrctype == iwtvap
                                    and idsttype == iwtice
                                ):
                                    std_vap_idx = _idx3(
                                        i,
                                        k,
                                        wtrc_iawset[_iawset_idx(iwtvap, 1, pwtype)],
                                        pcols,
                                        pver,
                                    )
                                    ispec = iwspec[mdst - 1]
                                    alpha = _wtrc_get_alpha(
                                        qloc0[std_vap_idx],
                                        tloc[idx2],
                                        ispec,
                                        isrctype,
                                        idsttype,
                                        1,
                                        state_pmid[idx2],
                                        1,
                                        wisotope,
                                        iwtvap,
                                        iwtliq,
                                        epsilo,
                                    )
                                    fr = qloc[idx_mbase] / qloc0[idx_mbase]
                                    if fr < 0.0:
                                        fr = 0.0
                                    if fr > 1.0:
                                        fr = 1.0
                                    qloc[idx_msrc] = qloc0[idx_msrc] * (fr**alpha)
                                    qloc[idx_mdst] = qloc[idx_mdst] + (qloc0[idx_msrc] - qloc[idx_msrc])
                                else:
                                    qloc[idx_mdst] = (
                                        qloc[idx_mdst]
                                        + R * rate_val * dtime / wtrc_niter
                                    )
                                    if isrctype != idsttype:
                                        qloc[idx_msrc] = (
                                            qloc[idx_msrc]
                                            - R * rate_val * dtime / wtrc_niter
                                        )

                            for icnst in range(1, wtrc_ncnst + 1):
                                trc_idx = wtrc_indices[icnst - 1]
                                idx3 = _idx3(i, k, trc_idx, pcols, pver)
                                qloc0[idx3] = qloc[idx3]

                    if wisotope != 0:
                        for iwset in range(2, wtrc_nwset + 1):
                            std_vap = wtrc_iawset[_iawset_idx(iwtvap, 1, pwtype)]
                            std_liq = wtrc_iawset[_iawset_idx(iwtliq, 1, pwtype)]
                            iso_vap = wtrc_iawset[_iawset_idx(iwtvap, iwset, pwtype)]
                            iso_liq = wtrc_iawset[_iawset_idx(iwtliq, iwset, pwtype)]

                            idx_std_vap = _idx3(i, k, std_vap, pcols, pver)
                            idx_std_liq = _idx3(i, k, std_liq, pcols, pver)
                            idx_iso_vap = _idx3(i, k, iso_vap, pcols, pver)
                            idx_iso_liq = _idx3(i, k, iso_liq, pcols, pver)

                            alpha = _wtrc_get_alpha(
                                qloc0[idx_std_vap],
                                tloc[idx2],
                                iwspec[iso_vap - 1],
                                iwtvap,
                                iwtliq,
                                0,
                                1.0,
                                0,
                                wisotope,
                                iwtvap,
                                iwtliq,
                                epsilo,
                            )
                            vapiso, liqiso = _wtrc_liqvap_equil(
                                alpha,
                                1.0,
                                qloc[idx_std_vap],
                                qloc[idx_std_liq],
                                qloc[idx_iso_vap],
                                qloc[idx_iso_liq],
                                wtrc_qmin,
                                rstd,
                            )
                            qloc[idx_iso_vap] = vapiso
                            qloc[idx_iso_liq] = liqiso

                        for icnst in range(1, wtrc_ncnst + 1):
                            trc_idx = wtrc_indices[icnst - 1]
                            idx3 = _idx3(i, k, trc_idx, pcols, pver)
                            qloc0[idx3] = qloc[idx3]

                tloc[idx2] = state_t[idx2] + prelat[idx2] / cpair * dtime / wtrc_niter

    for itype in range(1, pwtype + 1):
        for k in range(top_lev, pver + 1):
            for i in range(1, ncol + 1):
                diff[_idx3(i, k, itype, pcols, pver)] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            for icnst in range(1, wtrc_ncnst + 1):
                trc_idx = wtrc_indices[icnst - 1]
                idx3 = _idx3(i, k, trc_idx, pcols, pver)
                ptend_q[idx3] = (qloc[idx3] - state_q[idx3]) / dtime

                if icnst <= pwtype:
                    bulk_idx = wtrc_bulk_indices[_bulk_idx(icnst)]
                    diff[_idx3(i, k, icnst, pcols, pver)] = ptend_q[idx3] - ptend_q[
                        _idx3(i, k, bulk_idx, pcols, pver)
                    ]

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            for itype in range(1, pwtype + 1):
                qtmp = 0.0
                diff_idx = _idx3(i, k, itype, pcols, pver)
                for iwset in range(1, wtrc_nwset + 1):
                    trc_idx = wtrc_iatype[_iatype_idx(iwset, itype, wtrc_nwset)]
                    idx3 = _idx3(i, k, trc_idx, pcols, pver)
                    if iwset == 1:
                        qtmp = ptend_q[idx3]
                    R = _wtrc_ratio(
                        iwspec[trc_idx - 1],
                        ptend_q[idx3],
                        qtmp,
                        wtrc_qmin,
                        rstd,
                    )
                    ptend_q[idx3] = ptend_q[idx3] - R * diff[diff_idx]

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            for itype in range(1, pwtype + 1):
                qtmp = 0.0
                diff_idx = _idx3(i, k, itype, pcols, pver)
                for iwset in range(1, wtrc_nwset + 1):
                    trc_idx = wtrc_iatype[_iatype_idx(iwset, itype, wtrc_nwset)]
                    idx3 = _idx3(i, k, trc_idx, pcols, pver)
                    if iwset == 1:
                        bulk_idx = wtrc_bulk_indices[_bulk_idx(itype)]
                        diff[diff_idx] = ptend_q[idx3] - ptend_q[
                            _idx3(i, k, bulk_idx, pcols, pver)
                        ]
                        qtmp = ptend_q[idx3]
                    R = _wtrc_ratio(
                        iwspec[trc_idx - 1],
                        ptend_q[idx3],
                        qtmp,
                        wtrc_qmin,
                        rstd,
                    )
                    ptend_q[idx3] = ptend_q[idx3] - R * diff[diff_idx]


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
def macrop_driver_ptend_lq_mask_shell_codon(
    mode: int,
    pcnst: int,
    wtrc_nwset: int,
    wtrc_ncnst: int,
    use_water_tracers: int,
    ixcldliq: int,
    ixcldice: int,
    ixnumliq: int,
    ixnumice: int,
    lq_mask_p: cobj,
    liq_type_p: cobj,
    ice_type_p: cobj,
    wtrc_indices_p: cobj,
):
    lq_mask = Ptr[int](lq_mask_p)
    liq_type = Ptr[int](liq_type_p)
    ice_type = Ptr[int](ice_type_p)
    wtrc_indices = Ptr[int](wtrc_indices_p)

    for m in range(1, pcnst + 1):
        lq_mask[m - 1] = 0

    if mode == 1:
        lq_mask[ixcldliq - 1] = 1
        lq_mask[ixcldice - 1] = 1
        lq_mask[ixnumliq - 1] = 1
        lq_mask[ixnumice - 1] = 1
        if use_water_tracers != 0:
            for m in range(1, wtrc_nwset + 1):
                lq_mask[liq_type[m - 1] - 1] = 1
                lq_mask[ice_type[m - 1] - 1] = 1
    elif mode == 2:
        lq_mask[0] = 1
        lq_mask[ixcldice - 1] = 1
        lq_mask[ixcldliq - 1] = 1
        lq_mask[ixnumliq - 1] = 1
        lq_mask[ixnumice - 1] = 1
        for m in range(1, wtrc_ncnst + 1):
            lq_mask[wtrc_indices[m - 1] - 1] = 1


@export
def macrop_driver_detrain_init_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    dlf_T_p: cobj,
    dlf_qv_p: cobj,
    dlf_ql_p: cobj,
    dlf_qi_p: cobj,
    dlf_nl_p: cobj,
    dlf_ni_p: cobj,
    det_s_p: cobj,
    det_ice_p: cobj,
    dpdlfliq_p: cobj,
    dpdlfice_p: cobj,
    shdlfliq_p: cobj,
    shdlfice_p: cobj,
    dpdlft_p: cobj,
    shdlft_p: cobj,
):
    dlf_T = Ptr[float](dlf_T_p)
    dlf_qv = Ptr[float](dlf_qv_p)
    dlf_ql = Ptr[float](dlf_ql_p)
    dlf_qi = Ptr[float](dlf_qi_p)
    dlf_nl = Ptr[float](dlf_nl_p)
    dlf_ni = Ptr[float](dlf_ni_p)
    det_s = Ptr[float](det_s_p)
    det_ice = Ptr[float](det_ice_p)
    dpdlfliq = Ptr[float](dpdlfliq_p)
    dpdlfice = Ptr[float](dpdlfice_p)
    shdlfliq = Ptr[float](shdlfliq_p)
    shdlfice = Ptr[float](shdlfice_p)
    dpdlft = Ptr[float](dpdlft_p)
    shdlft = Ptr[float](shdlft_p)

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            idx2 = _idx2(i, k, pcols)
            dlf_T[idx2] = 0.0
            dlf_qv[idx2] = 0.0
            dlf_ql[idx2] = 0.0
            dlf_qi[idx2] = 0.0
            dlf_nl[idx2] = 0.0
            dlf_ni[idx2] = 0.0
            dpdlfliq[idx2] = 0.0
            dpdlfice[idx2] = 0.0
            shdlfliq[idx2] = 0.0
            shdlfice[idx2] = 0.0
            dpdlft[idx2] = 0.0
            shdlft[idx2] = 0.0

    for i in range(1, pcols + 1):
        idx1 = i - 1
        det_s[idx1] = 0.0
        det_ice[idx1] = 0.0


@export
def macrop_driver_detrain_init_lq_mask_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    wtrc_nwset: int,
    wtrc_ncnst: int,
    use_water_tracers: int,
    ixcldliq: int,
    ixcldice: int,
    ixnumliq: int,
    ixnumice: int,
    dlf_T_p: cobj,
    dlf_qv_p: cobj,
    dlf_ql_p: cobj,
    dlf_qi_p: cobj,
    dlf_nl_p: cobj,
    dlf_ni_p: cobj,
    det_s_p: cobj,
    det_ice_p: cobj,
    dpdlfliq_p: cobj,
    dpdlfice_p: cobj,
    shdlfliq_p: cobj,
    shdlfice_p: cobj,
    dpdlft_p: cobj,
    shdlft_p: cobj,
    lq_mask_p: cobj,
    liq_type_p: cobj,
    ice_type_p: cobj,
    wtrc_indices_p: cobj,
):
    macrop_driver_detrain_init_shell_codon(
        ncol,
        pcols,
        pver,
        dlf_T_p,
        dlf_qv_p,
        dlf_ql_p,
        dlf_qi_p,
        dlf_nl_p,
        dlf_ni_p,
        det_s_p,
        det_ice_p,
        dpdlfliq_p,
        dpdlfice_p,
        shdlfliq_p,
        shdlfice_p,
        dpdlft_p,
        shdlft_p,
    )
    macrop_driver_ptend_lq_mask_shell_codon(
        1,
        pcnst,
        wtrc_nwset,
        wtrc_ncnst,
        use_water_tracers,
        ixcldliq,
        ixcldice,
        ixnumliq,
        ixnumice,
        lq_mask_p,
        liq_type_p,
        ice_type_p,
        wtrc_indices_p,
    )


@export
def macrop_driver_detrain_post_shell_codon(
    ncol: int,
    pcols: int,
    det_ice_p: cobj,
):
    det_ice = Ptr[float](det_ice_p)

    for i in range(1, ncol + 1):
        idx1 = i - 1
        det_ice[idx1] = det_ice[idx1] / 1000.0


@export
def macrop_driver_mmacro_input_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    top_lev: int,
    ixcldliq: int,
    ixcldice: int,
    ixnumliq: int,
    ixnumice: int,
    state_q_p: cobj,
    zeros_p: cobj,
    qc_p: cobj,
    qi_p: cobj,
    nc_p: cobj,
    ni_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    zeros = Ptr[float](zeros_p)
    qc = Ptr[float](qc_p)
    qi = Ptr[float](qi_p)
    nc = Ptr[float](nc_p)
    ni = Ptr[float](ni_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            zeros[idx2] = 0.0
            qc[idx2] = state_q[_idx3(i, k, ixcldliq, pcols, pver)]
            qi[idx2] = state_q[_idx3(i, k, ixcldice, pcols, pver)]
            nc[idx2] = state_q[_idx3(i, k, ixnumliq, pcols, pver)]
            ni[idx2] = state_q[_idx3(i, k, ixnumice, pcols, pver)]


@export
def macrop_driver_mmacro_post_fields_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    fice_p: cobj,
    alst_p: cobj,
    aist_p: cobj,
    fice_ql_p: cobj,
    ast_p: cobj,
):
    fice = Ptr[float](fice_p)
    alst = Ptr[float](alst_p)
    aist = Ptr[float](aist_p)
    fice_ql = Ptr[float](fice_ql_p)
    ast = Ptr[float](ast_p)

    for k in range(1, top_lev):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            fice_ql[idx2] = 0.0
            ast[idx2] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            fice_ql[idx2] = fice[idx2]
            ast[idx2] = max(alst[idx2], aist[idx2])


@export
def macrop_driver_mmacro_config_check_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    do_cldice: int,
    do_cldliq: int,
    qiten_p: cobj,
    niten_p: cobj,
    qcten_p: cobj,
    ncten_p: cobj,
    mask_p: cobj,
):
    qiten = Ptr[float](qiten_p)
    niten = Ptr[float](niten_p)
    qcten = Ptr[float](qcten_p)
    ncten = Ptr[float](ncten_p)
    mask_out = Ptr[int](mask_p)

    mask = 0
    if do_cldice == 0:
        qiten_nonzero = False
        niten_nonzero = False
        for k in range(top_lev, pver + 1):
            for i in range(1, ncol + 1):
                idx2 = _idx2(i, k, pcols)
                if qiten[idx2] != 0.0:
                    qiten_nonzero = True
                if niten[idx2] != 0.0:
                    niten_nonzero = True
        if qiten_nonzero:
            mask |= 1
        if niten_nonzero:
            mask |= 2

    if do_cldliq == 0:
        qcten_nonzero = False
        ncten_nonzero = False
        for k in range(top_lev, pver + 1):
            for i in range(1, ncol + 1):
                idx2 = _idx2(i, k, pcols)
                if qcten[idx2] != 0.0:
                    qcten_nonzero = True
                if ncten[idx2] != 0.0:
                    ncten_nonzero = True
        if qcten_nonzero:
            mask |= 4
        if ncten_nonzero:
            mask |= 8

    mask_out[0] = mask


@export
def macrop_driver_cfmip_diag_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    cld_p: cobj,
    state_ql_p: cobj,
    state_qi_p: cobj,
    mr_ccliq_p: cobj,
    mr_ccice_p: cobj,
    mr_lsliq_p: cobj,
    mr_lsice_p: cobj,
):
    cld = Ptr[float](cld_p)
    state_ql = Ptr[float](state_ql_p)
    state_qi = Ptr[float](state_qi_p)
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
def macrop_driver_wtrc_detrain_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    wtrc_nwset: int,
    state_t_p: cobj,
    wtdlf_p: cobj,
    liq_type_p: cobj,
    ice_type_p: cobj,
    ptend_q_p: cobj,
):
    state_t = Ptr[float](state_t_p)
    wtdlf = Ptr[float](wtdlf_p)
    liq_type = Ptr[int](liq_type_p)
    ice_type = Ptr[int](ice_type_p)
    ptend_q = Ptr[float](ptend_q_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            if state_t[idx2] > 268.15:
                dum1 = 0.0
            elif state_t[idx2] < 238.15:
                dum1 = 1.0
            else:
                dum1 = (268.15 - state_t[idx2]) / 30.0
            for m in range(1, wtrc_nwset + 1):
                idx_wtdlf = _idx3(i, k, m, pcols, pver)
                ptend_q[_idx3(i, k, liq_type[m - 1], pcols, pver)] = wtdlf[idx_wtdlf] * (1.0 - dum1)
                ptend_q[_idx3(i, k, ice_type[m - 1], pcols, pver)] = wtdlf[idx_wtdlf] * dum1


@export
def macrop_driver_clr_old_diag_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    concld_p: cobj,
    alst_p: cobj,
    ast_p: cobj,
    concld_old_p: cobj,
    clrw_old_p: cobj,
    clri_old_p: cobj,
):
    concld = Ptr[float](concld_p)
    alst = Ptr[float](alst_p)
    ast = Ptr[float](ast_p)
    concld_old = Ptr[float](concld_old_p)
    clrw_old = Ptr[float](clrw_old_p)
    clri_old = Ptr[float](clri_old_p)

    for k in range(1, top_lev):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            clrw_old[idx2] = 0.0
            clri_old[idx2] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            concld_old[idx2] = concld[idx2]
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
def macrop_driver_mmacro_prepare_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    top_lev: int,
    ixcldliq: int,
    ixcldice: int,
    ixnumliq: int,
    ixnumice: int,
    nstep: int,
    rdtime: float,
    state_q_p: cobj,
    state_t_p: cobj,
    state_qv_p: cobj,
    zeros_p: cobj,
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
    macrop_driver_mmacro_input_shell_codon(
        ncol,
        pcols,
        pver,
        pcnst,
        top_lev,
        ixcldliq,
        ixcldice,
        ixnumliq,
        ixnumice,
        state_q_p,
        zeros_p,
        qc_p,
        qi_p,
        nc_p,
        ni_p,
    )
    macrop_driver_forcing_prep_codon(
        ncol,
        pcols,
        pver,
        top_lev,
        nstep,
        rdtime,
        state_t_p,
        state_qv_p,
        qc_p,
        qi_p,
        nc_p,
        ni_p,
        tcwat_p,
        qcwat_p,
        lcwat_p,
        iccwat_p,
        nlwat_p,
        niwat_p,
        cc_t_p,
        cc_qv_p,
        cc_ql_p,
        cc_qi_p,
        cc_nl_p,
        cc_ni_p,
        cc_qlst_p,
        ttend_p,
        qtend_p,
        ltend_p,
        itend_p,
        nltend_p,
        nitend_p,
        lmitend_p,
        t_inout_p,
        qv_inout_p,
        ql_inout_p,
        qi_inout_p,
        nl_inout_p,
        ni_inout_p,
    )


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
def macrop_driver_ptend_config_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    do_cldice: int,
    do_cldliq: int,
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
    mask_p: cobj,
):
    macrop_driver_ptend_assign_codon(
        ncol,
        pcols,
        pver,
        top_lev,
        tlat_p,
        qvlat_p,
        qcten_p,
        qiten_p,
        ncten_p,
        niten_p,
        ptend_s_p,
        ptend_qv_p,
        ptend_ql_p,
        ptend_qi_p,
        ptend_nl_p,
        ptend_ni_p,
    )
    macrop_driver_mmacro_config_check_codon(
        ncol,
        pcols,
        pver,
        top_lev,
        do_cldice,
        do_cldliq,
        qiten_p,
        niten_p,
        qcten_p,
        ncten_p,
        mask_p,
    )


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
def macrop_driver_wtrc_process_rates_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pwtype: int,
    top_lev: int,
    iwtvap: int,
    iwtliq: int,
    iwtice: int,
    qvlat_p: cobj,
    qcten_p: cobj,
    qiten_p: cobj,
    process_rates_p: cobj,
):
    qvlat = Ptr[float](qvlat_p)
    qcten = Ptr[float](qcten_p)
    qiten = Ptr[float](qiten_p)
    process_rates = Ptr[float](process_rates_p)

    for rtype in range(1, pwtype + 1):
        for isrctype in range(1, pwtype + 1):
            for idsttype in range(1, pwtype + 1):
                for k in range(top_lev, pver + 1):
                    for i in range(1, pcols + 1):
                        process_rates[
                            _process_rates_idx(
                                i,
                                k,
                                idsttype,
                                isrctype,
                                rtype,
                                pcols,
                                pver,
                                pwtype,
                            )
                        ] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)

            rate_val = qvlat[idx2] + qcten[idx2]
            rate_val = rate_val + qiten[idx2]
            process_rates[
                _process_rates_idx(
                    i, k, iwtvap, iwtvap, iwtvap, pcols, pver, pwtype
                )
            ] = process_rates[
                _process_rates_idx(
                    i, k, iwtvap, iwtvap, iwtvap, pcols, pver, pwtype
                )
            ] + rate_val

            rate_val = qcten[idx2]
            if rate_val < 0.0:
                dst_idx = _process_rates_idx(
                    i, k, iwtliq, iwtvap, iwtliq, pcols, pver, pwtype
                )
                src_idx = _process_rates_idx(
                    i, k, iwtvap, iwtliq, iwtliq, pcols, pver, pwtype
                )
            else:
                dst_idx = _process_rates_idx(
                    i, k, iwtliq, iwtvap, iwtvap, pcols, pver, pwtype
                )
                src_idx = _process_rates_idx(
                    i, k, iwtvap, iwtliq, iwtvap, pcols, pver, pwtype
                )
            process_rates[dst_idx] = process_rates[dst_idx] + rate_val
            process_rates[src_idx] = process_rates[src_idx] - rate_val

            rate_val = qiten[idx2]
            if rate_val < 0.0:
                dst_idx = _process_rates_idx(
                    i, k, iwtice, iwtvap, iwtice, pcols, pver, pwtype
                )
                src_idx = _process_rates_idx(
                    i, k, iwtvap, iwtice, iwtice, pcols, pver, pwtype
                )
            else:
                dst_idx = _process_rates_idx(
                    i, k, iwtice, iwtvap, iwtvap, pcols, pver, pwtype
                )
                src_idx = _process_rates_idx(
                    i, k, iwtvap, iwtice, iwtvap, pcols, pver, pwtype
                )
            process_rates[dst_idx] = process_rates[dst_idx] + rate_val
            process_rates[src_idx] = process_rates[src_idx] - rate_val


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

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            cldsice[_idx2(i, k, pcols)] = 0.0

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


@export
def macrop_driver_detrain_core_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    do_detrain: int,
    cu_det_st: int,
    cpair: float,
    gravit: float,
    latice: float,
    nl_denom_a: float,
    nl_denom_b: float,
    ni_denom_a: float,
    ni_denom_b: float,
    state_t_p: cobj,
    state_pdel_p: cobj,
    dlf_p: cobj,
    dlf2_p: cobj,
    ptend_ql_p: cobj,
    ptend_qi_p: cobj,
    ptend_nl_p: cobj,
    ptend_ni_p: cobj,
    ptend_s_p: cobj,
    det_s_p: cobj,
    det_ice_p: cobj,
    dlf_t_p: cobj,
    dlf_qv_p: cobj,
    dlf_ql_p: cobj,
    dlf_qi_p: cobj,
    dlf_nl_p: cobj,
    dlf_ni_p: cobj,
    dpdlfliq_p: cobj,
    dpdlfice_p: cobj,
    shdlfliq_p: cobj,
    shdlfice_p: cobj,
    dpdlft_p: cobj,
    shdlft_p: cobj,
):
    # Fortran mappings: state_t/state_pdel/dlf/... are real(r8) arrays with shape (pcols,pver);
    # det_s/det_ice are real(r8) arrays with shape (pcols).
    state_t = Ptr[float](state_t_p)
    state_pdel = Ptr[float](state_pdel_p)
    dlf = Ptr[float](dlf_p)
    dlf2 = Ptr[float](dlf2_p)
    ptend_ql = Ptr[float](ptend_ql_p)
    ptend_qi = Ptr[float](ptend_qi_p)
    ptend_nl = Ptr[float](ptend_nl_p)
    ptend_ni = Ptr[float](ptend_ni_p)
    ptend_s = Ptr[float](ptend_s_p)
    det_s = Ptr[float](det_s_p)
    det_ice = Ptr[float](det_ice_p)
    dlf_t = Ptr[float](dlf_t_p)
    dlf_qv = Ptr[float](dlf_qv_p)
    dlf_ql = Ptr[float](dlf_ql_p)
    dlf_qi = Ptr[float](dlf_qi_p)
    dlf_nl = Ptr[float](dlf_nl_p)
    dlf_ni = Ptr[float](dlf_ni_p)
    dpdlfliq = Ptr[float](dpdlfliq_p)
    dpdlfice = Ptr[float](dpdlfice_p)
    shdlfliq = Ptr[float](shdlfliq_p)
    shdlfice = Ptr[float](shdlfice_p)
    dpdlft = Ptr[float](dpdlft_p)
    shdlft = Ptr[float](shdlft_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            idx1 = i - 1

            if state_t[idx2] > 268.15:
                dum1_local = 0.0
            elif state_t[idx2] < 238.15:
                dum1_local = 1.0
            else:
                dum1_local = (268.15 - state_t[idx2]) / 30.0

            if do_detrain != 0:
                ptend_ql[idx2] = dlf[idx2] * (1.0 - dum1_local)
                ptend_qi[idx2] = dlf[idx2] * dum1_local
                ptend_nl[idx2] = (
                    3.0
                    * (max(0.0, (dlf[idx2] - dlf2[idx2])) * (1.0 - dum1_local))
                    / nl_denom_a
                    + 3.0
                    * (dlf2[idx2] * (1.0 - dum1_local))
                    / nl_denom_b
                )
                ptend_ni[idx2] = (
                    3.0
                    * (max(0.0, (dlf[idx2] - dlf2[idx2])) * dum1_local)
                    / ni_denom_a
                    + 3.0
                    * (dlf2[idx2] * dum1_local)
                    / ni_denom_b
                )
                ptend_s[idx2] = dlf[idx2] * dum1_local * latice
            else:
                ptend_ql[idx2] = 0.0
                ptend_qi[idx2] = 0.0
                ptend_nl[idx2] = 0.0
                ptend_ni[idx2] = 0.0
                ptend_s[idx2] = 0.0

            det_s[idx1] = det_s[idx1] + ptend_s[idx2] * state_pdel[idx2] / gravit
            det_ice[idx1] = det_ice[idx1] - ptend_qi[idx2] * state_pdel[idx2] / gravit

            if cu_det_st != 0:
                dlf_t[idx2] = ptend_s[idx2] / cpair
                dlf_qv[idx2] = 0.0
                dlf_ql[idx2] = ptend_ql[idx2]
                dlf_qi[idx2] = ptend_qi[idx2]
                dlf_nl[idx2] = ptend_nl[idx2]
                dlf_ni[idx2] = ptend_ni[idx2]
                ptend_ql[idx2] = 0.0
                ptend_qi[idx2] = 0.0
                ptend_nl[idx2] = 0.0
                ptend_ni[idx2] = 0.0
                ptend_s[idx2] = 0.0
                dpdlfliq[idx2] = 0.0
                dpdlfice[idx2] = 0.0
                shdlfliq[idx2] = 0.0
                shdlfice[idx2] = 0.0
                dpdlft[idx2] = 0.0
                shdlft[idx2] = 0.0
            else:
                dpdlfliq[idx2] = (dlf[idx2] - dlf2[idx2]) * (1.0 - dum1_local)
                dpdlfice[idx2] = (dlf[idx2] - dlf2[idx2]) * dum1_local
                shdlfliq[idx2] = dlf2[idx2] * (1.0 - dum1_local)
                shdlfice[idx2] = dlf2[idx2] * dum1_local
                dpdlft[idx2] = (dlf[idx2] - dlf2[idx2]) * dum1_local * latice / cpair
                shdlft[idx2] = dlf2[idx2] * dum1_local * latice / cpair


@export
def macrop_driver_stage_dispatch_codon(
    stage: int,
    mode: int,
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    top_lev: int,
    wtrc_nwset: int,
    wtrc_ncnst: int,
    use_water_tracers: int,
    ixcldliq: int,
    ixcldice: int,
    ixnumliq: int,
    ixnumice: int,
    do_detrain: int,
    cu_det_st: int,
    do_cldice: int,
    do_cldliq: int,
    nstep: int,
    cpair: float,
    gravit: float,
    latice: float,
    nl_denom_a: float,
    nl_denom_b: float,
    ni_denom_a: float,
    ni_denom_b: float,
    rdtime: float,
    tmelt: float,
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
):
    if stage == 1:
        macrop_driver_ptend_lq_mask_shell_codon(
            mode,
            pcnst,
            wtrc_nwset,
            wtrc_ncnst,
            use_water_tracers,
            ixcldliq,
            ixcldice,
            ixnumliq,
            ixnumice,
            p1,
            p2,
            p3,
            p4,
        )
    elif stage == 2:
        macrop_driver_detrain_init_lq_mask_shell_codon(
            ncol,
            pcols,
            pver,
            pcnst,
            wtrc_nwset,
            wtrc_ncnst,
            use_water_tracers,
            ixcldliq,
            ixcldice,
            ixnumliq,
            ixnumice,
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
        )
    elif stage == 3:
        macrop_driver_detrain_core_codon(
            ncol,
            pcols,
            pver,
            top_lev,
            do_detrain,
            cu_det_st,
            cpair,
            gravit,
            latice,
            nl_denom_a,
            nl_denom_b,
            ni_denom_a,
            ni_denom_b,
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
    elif stage == 4:
        macrop_driver_detrain_post_shell_codon(ncol, pcols, p1)
    elif stage == 5:
        macrop_driver_mmacro_prepare_shell_codon(
            ncol,
            pcols,
            pver,
            pcnst,
            top_lev,
            ixcldliq,
            ixcldice,
            ixnumliq,
            ixnumice,
            nstep,
            rdtime,
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
        )
    elif stage == 6:
        macrop_driver_mmacro_post_fields_shell_codon(ncol, pcols, pver, top_lev, p1, p2, p3, p4, p5)
    elif stage == 7:
        macrop_driver_ptend_config_shell_codon(
            ncol,
            pcols,
            pver,
            top_lev,
            do_cldice,
            do_cldliq,
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
        )
    elif stage == 8:
        macrop_driver_cfmip_diag_shell_codon(ncol, pcols, pver, top_lev, p1, p2, p3, p4, p5, p6, p7)
    elif stage == 9:
        macrop_driver_clr_old_diag_codon(ncol, pcols, pver, top_lev, p1, p2, p3, p4, p5, p6)
    elif stage == 10:
        macrop_driver_store_state_codon(
            ncol,
            pcols,
            pver,
            top_lev,
            tmelt,
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
        )
