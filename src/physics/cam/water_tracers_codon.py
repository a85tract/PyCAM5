@inline
def _state_q_idx(i: int, k: int, m: int, pcols: int, pver: int) -> int:
    """state%q(pcols, pver, pcnst)"""
    return (i - 1) + (k - 1) * pcols + (m - 1) * pcols * pver


@inline
def _field1_idx(i: int) -> int:
    """state%lat/state%phis(pcols)"""
    return i - 1


@inline
def _iatype_idx(m: int, p: int, wtrc_nwset: int) -> int:
    """wtrc_iatype64(wtrc_nwset, pwtype)"""
    return (m - 1) + (p - 1) * wtrc_nwset


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
