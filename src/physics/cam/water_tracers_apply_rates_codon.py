@inline
def _idx2(i: int, k: int, pcols: int) -> int:
    """water_tracers arrays declared as (pcols,pver)."""
    return (i - 1) + (k - 1) * pcols


@inline
def _idx3(i: int, k: int, m: int, pcols: int, pver: int) -> int:
    """water_tracers arrays declared as (pcols,pver,pcnst/pwtype)."""
    return (i - 1) + (k - 1) * pcols + (m - 1) * pcols * pver


@inline
def _idx_rmass(i: int, m: int, pcols: int) -> int:
    """rmass/smass arrays declared as (pcols,wtrc_nwset)."""
    return (i - 1) + (m - 1) * pcols


@inline
def _idx_iatype(m: int, icnst: int, wtrc_nwset: int) -> int:
    """wtrc_iatype64 declared as (wtrc_nwset,pwtype)."""
    return (m - 1) + (icnst - 1) * wtrc_nwset


@inline
def _idx_iawset(itype: int, iwset: int, pwtype: int) -> int:
    """wtrc_iawset64 declared as (pwtype,wtrc_nwset)."""
    return (itype - 1) + (iwset - 1) * pwtype


@inline
def _ratio_from_table(ispec: int, qtrc: float, qtot: float, qmin: float, rstd: Ptr[float]) -> float:
    if abs(qtot) < qmin:
        return rstd[ispec - 1]
    return qtrc / qtot


@inline
def _ratio_from_value(qtrc: float, qtot: float, qmin: float, rstd_value: float) -> float:
    if abs(qtot) < qmin:
        return rstd_value
    return qtrc / qtot


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
    pstate_q = Ptr[float](pstate_q_p)
    pstate_t = Ptr[float](pstate_t_p)
    qloc = Ptr[float](qloc_p)
    qloc0 = Ptr[float](qloc0_p)
    tloc = Ptr[float](tloc_p)

    for m in range(1, pcnst + 1):
        for k in range(top_lev, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx3(i, k, m, pcols, pver)
                qloc[idx] = pstate_q[idx]
                qloc0[idx] = qloc[idx]

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            tloc[idx2] = pstate_t[idx2]


def wtrc_apply_rates_zero_precip_codon(
    pcols: int,
    wtrc_nwset: int,
    rmass_p: cobj,
    smass_p: cobj,
    rmass0_p: cobj,
    smass0_p: cobj,
):
    rmass = Ptr[float](rmass_p)
    smass = Ptr[float](smass_p)
    rmass0 = Ptr[float](rmass0_p)
    smass0 = Ptr[float](smass0_p)

    for m in range(1, wtrc_nwset + 1):
        for i in range(1, pcols + 1):
            idx = _idx_rmass(i, m, pcols)
            rmass[idx] = 0.0
            smass[idx] = 0.0
            rmass0[idx] = 0.0
            smass0[idx] = 0.0


def wtrc_apply_rates_copy_qloc0_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    top_lev: int,
    qloc_p: cobj,
    qloc0_p: cobj,
):
    qloc = Ptr[float](qloc_p)
    qloc0 = Ptr[float](qloc0_p)

    for m in range(1, pcnst + 1):
        for k in range(top_lev, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx3(i, k, m, pcols, pver)
                qloc0[idx] = qloc[idx]


def wtrc_apply_rates_sync_level_state_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    pcnst: int,
    qloc_p: cobj,
    qloc0_p: cobj,
):
    qloc = Ptr[float](qloc_p)
    qloc0 = Ptr[float](qloc0_p)

    for m in range(1, pcnst + 1):
        idx = _idx3(i, k, m, pcols, pver)
        qloc0[idx] = qloc[idx]


def wtrc_apply_rates_sync_precip_column_codon(
    i: int,
    pcols: int,
    wtrc_nwset: int,
    rmass_p: cobj,
    smass_p: cobj,
    rmass0_p: cobj,
    smass0_p: cobj,
):
    rmass = Ptr[float](rmass_p)
    smass = Ptr[float](smass_p)
    rmass0 = Ptr[float](rmass0_p)
    smass0 = Ptr[float](smass0_p)

    for iwset in range(1, wtrc_nwset + 1):
        idx = _idx_rmass(i, iwset, pcols)
        rmass0[idx] = rmass[idx]
        smass0[idx] = smass[idx]


def wtrc_apply_rates_local_source_ratio_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    isrctype: int,
    iwset: int,
    iwtice: int,
    iwtstrain: int,
    msrc: int,
    mbase: int,
    qmin: float,
    rstd_value: float,
    qloc0_p: cobj,
    rmass0_p: cobj,
    smass0_p: cobj,
) -> float:
    qloc0 = Ptr[float](qloc0_p)
    rmass0 = Ptr[float](rmass0_p)
    smass0 = Ptr[float](smass0_p)

    if isrctype > iwtice:
        pidx = _idx_rmass(i, iwset, pcols)
        base_idx = _idx_rmass(i, 1, pcols)
        if isrctype == iwtstrain:
            return _ratio_from_value(rmass0[pidx], rmass0[base_idx], qmin, rstd_value)
        return _ratio_from_value(smass0[pidx], smass0[base_idx], qmin, rstd_value)

    qidx = _idx3(i, k, msrc, pcols, pver)
    qbase_idx = _idx3(i, k, mbase, pcols, pver)
    return _ratio_from_value(qloc0[qidx], qloc0[qbase_idx], qmin, rstd_value)


def wtrc_apply_rates_pre_temperature_begin_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    dtime: float,
    cpair: float,
    prelat_p: cobj,
    tloc_p: cobj,
):
    prelat = Ptr[float](prelat_p)
    tloc = Ptr[float](tloc_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            old_t = tloc[idx]
            tloc[idx] = (old_t + (old_t + prelat[idx] / cpair * dtime)) / 2.0


def wtrc_apply_rates_pre_temperature_end_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    dtime: float,
    cpair: float,
    niter: float,
    pstate_t_p: cobj,
    prelat_p: cobj,
    tloc_p: cobj,
):
    pstate_t = Ptr[float](pstate_t_p)
    prelat = Ptr[float](prelat_p)
    tloc = Ptr[float](tloc_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            tloc[idx] = pstate_t[idx] + prelat[idx] / cpair * dtime / niter


def wtrc_apply_rates_post_temperature_begin_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    dtime: float,
    cpair: float,
    postlat_p: cobj,
    tloc_p: cobj,
):
    postlat = Ptr[float](postlat_p)
    tloc = Ptr[float](tloc_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            old_t = tloc[idx]
            tloc[idx] = (old_t + (old_t + postlat[idx] / cpair * dtime)) / 2.0


def wtrc_apply_rates_post_temperature_end_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    dtime: float,
    cpair: float,
    niter: float,
    pstate_t_p: cobj,
    prelat_p: cobj,
    postlat_p: cobj,
    tloc_p: cobj,
):
    pstate_t = Ptr[float](pstate_t_p)
    prelat = Ptr[float](prelat_p)
    postlat = Ptr[float](postlat_p)
    tloc = Ptr[float](tloc_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            tloc[idx] = pstate_t[idx] + (prelat[idx] + postlat[idx]) / cpair * dtime / niter


def wtrc_apply_rates_precip_phase_codon(
    i: int,
    k: int,
    pcols: int,
    pwtype: int,
    wtrc_nwset: int,
    iwtstrain: int,
    dtime: float,
    qmin: float,
    meltso_ik: float,
    frzro_ik: float,
    wtrc_iawset_p: cobj,
    iwspec_p: cobj,
    rstd_p: cobj,
    rmass_p: cobj,
    smass_p: cobj,
    rmass0_p: cobj,
    smass0_p: cobj,
):
    wtrc_iawset = Ptr[int](wtrc_iawset_p)
    iwspec = Ptr[int](iwspec_p)
    rstd = Ptr[float](rstd_p)
    rmass = Ptr[float](rmass_p)
    smass = Ptr[float](smass_p)
    rmass0 = Ptr[float](rmass0_p)
    smass0 = Ptr[float](smass0_p)

    for iwset in range(1, wtrc_nwset + 1):
        mstrain = wtrc_iawset[_idx_iawset(iwtstrain, iwset, pwtype)]
        ispec = iwspec[mstrain - 1]

        idx = _idx_rmass(i, iwset, pcols)
        base_idx = _idx_rmass(i, 1, pcols)
        r_melt = _ratio_from_table(ispec, smass0[idx], smass0[base_idx], qmin, rstd)
        melt_delta = r_melt * meltso_ik * dtime
        rmass[idx] = rmass[idx] + melt_delta
        smass[idx] = smass[idx] - melt_delta

        r_freeze = _ratio_from_table(ispec, rmass0[idx], rmass0[base_idx], qmin, rstd)
        freeze_delta = r_freeze * frzro_ik * dtime
        smass[idx] = smass[idx] + freeze_delta
        rmass[idx] = rmass[idx] - freeze_delta

    for iwset in range(1, wtrc_nwset + 1):
        idx = _idx_rmass(i, iwset, pcols)
        rmass0[idx] = rmass[idx]
        smass0[idx] = smass[idx]


def wtrc_apply_rates_pre_normal_tendency_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    isrctype: int,
    idsttype: int,
    iwset: int,
    iwtstrain: int,
    iwtstsnow: int,
    msrc: int,
    mdst: int,
    ratio: float,
    rate: float,
    dtime: float,
    niter: float,
    pdel_ik: float,
    qloc_p: cobj,
    rmass_p: cobj,
    smass_p: cobj,
):
    qloc = Ptr[float](qloc_p)
    rmass = Ptr[float](rmass_p)
    smass = Ptr[float](smass_p)
    alpha = 1.0
    pidx = _idx_rmass(i, iwset, pcols)

    if idsttype == iwtstrain:
        rmass[pidx] = rmass[pidx] + (ratio * rate * dtime * pdel_ik) / niter
    elif idsttype == iwtstsnow:
        smass[pidx] = smass[pidx] + (ratio * rate * dtime * pdel_ik) / niter
    else:
        qdst = _idx3(i, k, mdst, pcols, pver)
        qloc[qdst] = qloc[qdst] + alpha * ratio * rate * dtime / niter

    if isrctype != idsttype:
        if isrctype == iwtstrain:
            rmass[pidx] = rmass[pidx] - (ratio * rate * dtime * pdel_ik) / niter
        elif isrctype == iwtstsnow:
            smass[pidx] = smass[pidx] - (ratio * rate * dtime * pdel_ik) / niter
        else:
            qsrc = _idx3(i, k, msrc, pcols, pver)
            qloc[qsrc] = qloc[qsrc] - alpha * ratio * rate * dtime / niter


def wtrc_apply_rates_pre_bergeron_direct_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    iwset: int,
    mdst: int,
    msrc: int,
    snow_mdst: int,
    ratio: float,
    rate: float,
    dtime: float,
    niter: float,
    pdel_ik: float,
    qloc_p: cobj,
    smass_p: cobj,
):
    qloc = Ptr[float](qloc_p)
    smass = Ptr[float](smass_p)

    if mdst == snow_mdst:
        pidx = _idx_rmass(i, iwset, pcols)
        smass[pidx] = smass[pidx] + (ratio * rate * dtime * pdel_ik) / niter
    else:
        qdst = _idx3(i, k, mdst, pcols, pver)
        qloc[qdst] = qloc[qdst] + ratio * rate * dtime / niter

    qsrc = _idx3(i, k, msrc, pcols, pver)
    qloc[qsrc] = qloc[qsrc] - ratio * rate * dtime / niter


def wtrc_apply_rates_post_normal_tendency_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    isrctype: int,
    idsttype: int,
    msrc: int,
    mdst: int,
    ratio: float,
    rate: float,
    dtime: float,
    niter: float,
    qloc_p: cobj,
):
    qloc = Ptr[float](qloc_p)
    alpha = 1.0

    qdst = _idx3(i, k, mdst, pcols, pver)
    qloc[qdst] = qloc[qdst] + alpha * ratio * rate * dtime / niter
    if isrctype != idsttype:
        qsrc = _idx3(i, k, msrc, pcols, pver)
        qloc[qsrc] = qloc[qsrc] - alpha * ratio * rate * dtime / niter


def wtrc_apply_rates_precip_error_correction_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    pcnst: int,
    pwtype: int,
    wtrc_nwset: int,
    iwtstrain: int,
    iwtvap: int,
    qmin: float,
    pdel_ik: float,
    wtrc_iawset_p: cobj,
    iwspec_p: cobj,
    rstd_p: cobj,
    qloc_p: cobj,
    qloc0_p: cobj,
    rmass_p: cobj,
    smass_p: cobj,
    rmass0_p: cobj,
    smass0_p: cobj,
):
    wtrc_iawset = Ptr[int](wtrc_iawset_p)
    iwspec = Ptr[int](iwspec_p)
    rstd = Ptr[float](rstd_p)
    qloc = Ptr[float](qloc_p)
    qloc0 = Ptr[float](qloc0_p)
    rmass = Ptr[float](rmass_p)
    smass = Ptr[float](smass_p)
    rmass0 = Ptr[float](rmass0_p)
    smass0 = Ptr[float](smass0_p)

    rain_pdiff = rmass[_idx_rmass(i, 2, pcols)] - rmass[_idx_rmass(i, 1, pcols)]
    for iwset in range(2, wtrc_nwset + 1):
        mstrain = wtrc_iawset[_idx_iawset(iwtstrain, iwset, pwtype)]
        ispec = iwspec[mstrain - 1]
        idx = _idx_rmass(i, iwset, pcols)
        r_corr = _ratio_from_table(ispec, rmass0[idx], rmass0[_idx_rmass(i, 2, pcols)], qmin, rstd)
        rain_delta = r_corr * rain_pdiff
        rmass[idx] = rmass[idx] - rain_delta
        mdst = wtrc_iawset[_idx_iawset(iwtvap, iwset, pwtype)]
        qidx = _idx3(i, k, mdst, pcols, pver)
        qloc[qidx] = qloc[qidx] + rain_delta / pdel_ik

    snow_pdiff = smass[_idx_rmass(i, 2, pcols)] - smass[_idx_rmass(i, 1, pcols)]
    for iwset in range(2, wtrc_nwset + 1):
        mstrain = wtrc_iawset[_idx_iawset(iwtstrain, iwset, pwtype)]
        ispec = iwspec[mstrain - 1]
        idx = _idx_rmass(i, iwset, pcols)
        r_corr = _ratio_from_table(ispec, smass0[idx], smass0[_idx_rmass(i, 2, pcols)], qmin, rstd)
        snow_delta = r_corr * snow_pdiff
        smass[idx] = smass[idx] - snow_delta
        mdst = wtrc_iawset[_idx_iawset(iwtvap, iwset, pwtype)]
        qidx = _idx3(i, k, mdst, pcols, pver)
        qloc[qidx] = qloc[qidx] + snow_delta / pdel_ik

    for m in range(1, pcnst + 1):
        qidx = _idx3(i, k, m, pcols, pver)
        qloc0[qidx] = qloc[qidx]

    for iwset in range(1, wtrc_nwset + 1):
        idx = _idx_rmass(i, iwset, pcols)
        rmass0[idx] = rmass[idx]
        smass0[idx] = smass[idx]


def wtrc_apply_rates_prepare_bulk_indices_codon(
    pwtype: int,
    bulk_indices_p: cobj,
    bulk_indices64_p: cobj,
):
    bulk_indices = Ptr[i32](bulk_indices_p)
    bulk_indices64 = Ptr[int](bulk_indices64_p)

    for idsttype in range(1, pwtype + 1):
        bulk_indices64[idsttype - 1] = int(bulk_indices[idsttype - 1])


def wtrc_apply_rates_prepare_net_indices_codon(
    pwtype: int,
    wtrc_ncnst: int,
    wtrc_indices_p: cobj,
    bulk_indices_p: cobj,
    wtrc_indices64_p: cobj,
    bulk_indices64_p: cobj,
):
    wtrc_indices = Ptr[i32](wtrc_indices_p)
    bulk_indices = Ptr[i32](bulk_indices_p)
    wtrc_indices64 = Ptr[int](wtrc_indices64_p)
    bulk_indices64 = Ptr[int](bulk_indices64_p)

    for icnst in range(1, wtrc_ncnst + 1):
        wtrc_indices64[icnst - 1] = int(wtrc_indices[icnst - 1])

    for icnst in range(1, pwtype + 1):
        bulk_indices64[icnst - 1] = int(bulk_indices[icnst - 1])


def wtrc_apply_rates_prepare_correction_indices_codon(
    pcnst: int,
    pwtype: int,
    wtrc_nwset: int,
    wtrc_max_cnst: int,
    wtrc_iatype_p: cobj,
    iwspec_p: cobj,
    wtrc_iatype64_p: cobj,
    iwspec64_p: cobj,
):
    wtrc_iatype = Ptr[i32](wtrc_iatype_p)
    iwspec = Ptr[i32](iwspec_p)
    wtrc_iatype64 = Ptr[int](wtrc_iatype64_p)
    iwspec64 = Ptr[int](iwspec64_p)

    for m in range(1, wtrc_nwset + 1):
        for icnst in range(1, pwtype + 1):
            src_idx = (m - 1) + (icnst - 1) * wtrc_max_cnst
            dst_idx = _idx_iatype(m, icnst, wtrc_nwset)
            wtrc_iatype64[dst_idx] = int(wtrc_iatype[src_idx])

    for ispec in range(1, pcnst + 1):
        iwspec64[ispec - 1] = int(iwspec[ispec - 1])


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
    bulk_indices = Ptr[int](bulk_indices_p)
    ptend_q = Ptr[float](ptend_q_p)
    qloc = Ptr[float](qloc_p)

    for idsttype in range(1, pwtype + 1):
        m = bulk_indices[idsttype - 1]
        for k in range(top_lev, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx3(i, k, m, pcols, pver)
                qloc[idx] = qloc[idx] + ptend_q[idx] * dtime


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
    wtrc_indices = Ptr[int](wtrc_indices_p)
    bulk_indices = Ptr[int](bulk_indices_p)
    pstate_q = Ptr[float](pstate_q_p)
    ptend_q = Ptr[float](ptend_q_p)
    qloc = Ptr[float](qloc_p)
    diff = Ptr[float](diff_p)

    for icnst in range(1, pwtype + 1):
        for k in range(1, pver + 1):
            for i in range(1, pcols + 1):
                diff[_idx3(i, k, icnst, pcols, pver)] = 0.0

    for icnst in range(1, wtrc_ncnst + 1):
        m = wtrc_indices[icnst - 1]
        for k in range(top_lev, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx3(i, k, m, pcols, pver)
                ptend_q[idx] = (qloc[idx] - pstate_q[idx]) / dtime
                if icnst <= pwtype:
                    bulk = bulk_indices[icnst - 1]
                    diff[_idx3(i, k, icnst, pcols, pver)] = ptend_q[idx] - ptend_q[
                        _idx3(i, k, bulk, pcols, pver)
                    ]


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
    wtrc_iatype = Ptr[int](wtrc_iatype_p)
    iwspec = Ptr[int](iwspec_p)
    rstd = Ptr[float](rstd_p)
    ptend_q = Ptr[float](ptend_q_p)
    diff = Ptr[float](diff_p)

    for i in range(1, ncol + 1):
        for k in range(top_lev, pver + 1):
            for icnst in range(1, pwtype + 1):
                qtmp = 0.0
                for m in range(1, wtrc_nwset + 1):
                    midx = wtrc_iatype[_idx_iatype(m, icnst, wtrc_nwset)]
                    qidx = _idx3(i, k, midx, pcols, pver)
                    if m == 1:
                        qtmp = ptend_q[qidx]

                    ispec = iwspec[midx - 1]
                    ratio = _ratio_from_table(ispec, ptend_q[qidx], qtmp, qmin, rstd)
                    ptend_q[qidx] = ptend_q[qidx] - ratio * diff[_idx3(i, k, icnst, pcols, pver)]


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
    wtrc_iatype = Ptr[int](wtrc_iatype_p)
    bulk_indices = Ptr[int](bulk_indices_p)
    iwspec = Ptr[int](iwspec_p)
    rstd = Ptr[float](rstd_p)
    ptend_q = Ptr[float](ptend_q_p)
    diff = Ptr[float](diff_p)

    for i in range(1, ncol + 1):
        for k in range(top_lev, pver + 1):
            for icnst in range(1, pwtype + 1):
                qtmp = 0.0
                for m in range(1, wtrc_nwset + 1):
                    midx = wtrc_iatype[_idx_iatype(m, icnst, wtrc_nwset)]
                    qidx = _idx3(i, k, midx, pcols, pver)
                    if m == 1:
                        bidx = bulk_indices[icnst - 1]
                        diff[_idx3(i, k, icnst, pcols, pver)] = ptend_q[qidx] - ptend_q[
                            _idx3(i, k, bidx, pcols, pver)
                        ]
                        qtmp = ptend_q[qidx]

                    ispec = iwspec[midx - 1]
                    ratio = _ratio_from_table(ispec, ptend_q[qidx], qtmp, qmin, rstd)
                    ptend_q[qidx] = ptend_q[qidx] - ratio * diff[_idx3(i, k, icnst, pcols, pver)]
