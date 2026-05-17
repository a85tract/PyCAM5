@inline
def _col_idx(icol: int) -> int:
    """shflx/lhflx/srfrpdel/excess/indxexc declared as (pcols)"""
    return icol - 1


@inline
def _q_idx(icol: int, mconst: int, pcols: int) -> int:
    """qbot/qflx/qfxo declared as (pcols, pcnst)"""
    return (icol - 1) + (mconst - 1) * pcols


@inline
def _spec_idx(ispec: int) -> int:
    """rstd declared as (pwtspec)"""
    return ispec - 1


@inline
def _wtrc_ratio(ispec: int, qtrc: float, qtot: float, wtrc_qmin: float, rstd) -> float:
    if abs(qtot) < wtrc_qmin:
        return rstd[_spec_idx(ispec)]
    return qtrc / qtot


@export
def qneg4_codon(
    ncol: int,
    pcols: int,
    pcnst: int,
    trace_water_on: int,
    iwtvap_n: int,
    qmin1: float,
    ztodt: float,
    gravit: float,
    latvap: float,
    wtrc_qmin: float,
    qbot_p: cobj,
    srfrpdel_p: cobj,
    shflx_p: cobj,
    lhflx_p: cobj,
    qflx_p: cobj,
    indxexc_p: cobj,
    excess_p: cobj,
    qfxo_p: cobj,
    nptsexc_p: cobj,
    worst_p: cobj,
    iw_p: cobj,
    wtrc_iatype_vap_p: cobj,
    iwspec_p: cobj,
    rstd_p: cobj,
):
    qbot = Ptr[float](qbot_p)
    srfrpdel = Ptr[float](srfrpdel_p)
    shflx = Ptr[float](shflx_p)
    lhflx = Ptr[float](lhflx_p)
    qflx = Ptr[float](qflx_p)
    indxexc = Ptr[int](indxexc_p)
    excess = Ptr[float](excess_p)
    qfxo = Ptr[float](qfxo_p)
    nptsexc = Ptr[int](nptsexc_p)
    worst = Ptr[float](worst_p)
    iw = Ptr[int](iw_p)
    wtrc_iatype_vap = Ptr[int](wtrc_iatype_vap_p)
    iwspec = Ptr[int](iwspec_p)
    rstd = Ptr[float](rstd_p)

    for m in range(1, pcnst + 1):
        for i in range(1, ncol + 1):
            qfxo[_q_idx(i, m, pcols)] = qflx[_q_idx(i, m, pcols)]

    nptsexc[0] = 0

    for i in range(1, ncol + 1):
        excess[_col_idx(i)] = qflx[_q_idx(i, 1, pcols)] - (
            qmin1 - qbot[_q_idx(i, 1, pcols)]
        ) / (ztodt * gravit * srfrpdel[_col_idx(i)])

        if excess[_col_idx(i)] < 0.0:
            nptsexc[0] += 1
            indxexc[_col_idx(nptsexc[0])] = i
            qflx[_q_idx(i, 1, pcols)] = qflx[_q_idx(i, 1, pcols)] - excess[_col_idx(i)]
            lhflx[_col_idx(i)] = lhflx[_col_idx(i)] - excess[_col_idx(i)] * latvap
            shflx[_col_idx(i)] = shflx[_col_idx(i)] + excess[_col_idx(i)] * latvap

    worst[0] = 0.0
    iw[0] = 0
    if nptsexc[0] > 10:
        for ii in range(1, nptsexc[0] + 1):
            i = indxexc[_col_idx(ii)]
            if excess[_col_idx(i)] < worst[0]:
                worst[0] = excess[_col_idx(i)]
                iw[0] = i

    if trace_water_on != 0:
        base_m = wtrc_iatype_vap[_col_idx(1)]
        for ivap in range(1, iwtvap_n + 1):
            m = wtrc_iatype_vap[_col_idx(ivap)]
            for ii in range(1, nptsexc[0] + 1):
                i = indxexc[_col_idx(ii)]
                rat = _wtrc_ratio(
                    iwspec[_col_idx(m)],
                    qfxo[_q_idx(i, m, pcols)],
                    qfxo[_q_idx(i, base_m, pcols)],
                    wtrc_qmin,
                    rstd,
                )
                qflx[_q_idx(i, m, pcols)] = qflx[_q_idx(i, m, pcols)] - rat * excess[_col_idx(i)]


@export
def qneg_batch_4_stage_dispatch_codon(
    ncol: int,
    pcols: int,
    pcnst: int,
    trace_water_on: int,
    iwtvap_n: int,
    qmin1: float,
    ztodt: float,
    gravit: float,
    latvap: float,
    wtrc_qmin: float,
    qbot_p: cobj,
    srfrpdel_p: cobj,
    shflx_p: cobj,
    lhflx_p: cobj,
    qflx_p: cobj,
    indxexc_p: cobj,
    excess_p: cobj,
    qfxo_p: cobj,
    nptsexc_p: cobj,
    worst_p: cobj,
    iw_p: cobj,
    wtrc_iatype_vap_p: cobj,
    iwspec_p: cobj,
    rstd_p: cobj,
):
    qneg4_codon(
        ncol,
        pcols,
        pcnst,
        trace_water_on,
        iwtvap_n,
        qmin1,
        ztodt,
        gravit,
        latvap,
        wtrc_qmin,
        qbot_p,
        srfrpdel_p,
        shflx_p,
        lhflx_p,
        qflx_p,
        indxexc_p,
        excess_p,
        qfxo_p,
        nptsexc_p,
        worst_p,
        iw_p,
        wtrc_iatype_vap_p,
        iwspec_p,
        rstd_p,
    )

@export
def qneg_batch_4_codon(
    ncol: int,
    pcols: int,
    pcnst: int,
    trace_water_on: int,
    iwtvap_n: int,
    qmin1: float,
    ztodt: float,
    gravit: float,
    latvap: float,
    wtrc_qmin: float,
    qbot_p: cobj,
    srfrpdel_p: cobj,
    shflx_p: cobj,
    lhflx_p: cobj,
    qflx_p: cobj,
    indxexc_p: cobj,
    excess_p: cobj,
    qfxo_p: cobj,
    nptsexc_p: cobj,
    worst_p: cobj,
    iw_p: cobj,
    wtrc_iatype_vap_p: cobj,
    iwspec_p: cobj,
    rstd_p: cobj,
):
    qneg_batch_4_stage_dispatch_codon(
        ncol,
        pcols,
        pcnst,
        trace_water_on,
        iwtvap_n,
        qmin1,
        ztodt,
        gravit,
        latvap,
        wtrc_qmin,
        qbot_p,
        srfrpdel_p,
        shflx_p,
        lhflx_p,
        qflx_p,
        indxexc_p,
        excess_p,
        qfxo_p,
        nptsexc_p,
        worst_p,
        iw_p,
        wtrc_iatype_vap_p,
        iwspec_p,
        rstd_p,
    )
