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
def _ratio_from_table(ispec: int, qtrc: float, qtot: float, qmin: float, rstd: Ptr[float]) -> float:
    if abs(qtot) < qmin:
        return rstd[ispec - 1]
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
