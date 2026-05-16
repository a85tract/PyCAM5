@inline
def _ptend_idx(i: int, k: int, psetcols: int) -> int:
    """ptend%s(psetcols, pver)"""
    return (i - 1) + (k - 1) * psetcols


@inline
def _pint_idx(i: int, k: int, psetcols: int) -> int:
    """state%pint(psetcols, pver+1)"""
    return (i - 1) + (k - 1) * psetcols


@inline
def _state2_idx(i: int, k: int, psetcols: int) -> int:
    """state%u/v/s/pdel(psetcols, pver)"""
    return (i - 1) + (k - 1) * psetcols


@inline
def _stateq_idx(i: int, k: int, m: int, psetcols: int, pver: int) -> int:
    """state%q(psetcols, pver, pcnst)"""
    return (i - 1) + (k - 1) * psetcols + (m - 1) * psetcols * pver


@export
def check_energy_fix_codon(
    ncol: int,
    pcols: int,
    pver: int,
    psetcols: int,
    heat_glob: float,
    gravit: float,
    state_pint_p: cobj,
    ptend_s_p: cobj,
    eshflx_p: cobj,
):
    state_pint = Ptr[float](state_pint_p)
    ptend_s = Ptr[float](ptend_s_p)
    eshflx = Ptr[float](eshflx_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            ptend_s[_ptend_idx(i, k, psetcols)] = heat_glob

    for i in range(1, ncol + 1):
        eshflx[i - 1] = heat_glob * (
            state_pint[_pint_idx(i, pver + 1, psetcols)] - state_pint[_pint_idx(i, 1, psetcols)]
        ) / gravit


@export
def check_energy_timestep_init_codon(
    ncol: int,
    pver: int,
    psetcols: int,
    pcnst: int,
    latvap: float,
    latice: float,
    gravit: float,
    ixcldliq: int,
    ixcldice: int,
    ixrain: int,
    ixsnow: int,
    state_u_p: cobj,
    state_v_p: cobj,
    state_s_p: cobj,
    state_q_p: cobj,
    state_pdel_p: cobj,
    ke_p: cobj,
    se_p: cobj,
    wv_p: cobj,
    wl_p: cobj,
    wi_p: cobj,
    state_te_ini_p: cobj,
    state_tw_ini_p: cobj,
):
    state_u = Ptr[float](state_u_p)
    state_v = Ptr[float](state_v_p)
    state_s = Ptr[float](state_s_p)
    state_q = Ptr[float](state_q_p)
    state_pdel = Ptr[float](state_pdel_p)
    ke = Ptr[float](ke_p)
    se = Ptr[float](se_p)
    wv = Ptr[float](wv_p)
    wl = Ptr[float](wl_p)
    wi = Ptr[float](wi_p)
    state_te_ini = Ptr[float](state_te_ini_p)
    state_tw_ini = Ptr[float](state_tw_ini_p)

    for i in range(1, ncol + 1):
        ke[i - 1] = 0.0
        se[i - 1] = 0.0
        wv[i - 1] = 0.0
        wl[i - 1] = 0.0
        wi[i - 1] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            pdel = state_pdel[_state2_idx(i, k, psetcols)]
            u = state_u[_state2_idx(i, k, psetcols)]
            v = state_v[_state2_idx(i, k, psetcols)]
            ke[i - 1] = ke[i - 1] + 0.5 * (u * u + v * v) * pdel / gravit
            se[i - 1] = se[i - 1] + state_s[_state2_idx(i, k, psetcols)] * pdel / gravit
            wv[i - 1] = wv[i - 1] + state_q[_stateq_idx(i, k, 1, psetcols, pver)] * pdel / gravit

    if ixcldliq > 1 and ixcldice > 1:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                pdel = state_pdel[_state2_idx(i, k, psetcols)]
                wl[i - 1] = wl[i - 1] + state_q[_stateq_idx(i, k, ixcldliq, psetcols, pver)] * pdel / gravit
                wi[i - 1] = wi[i - 1] + state_q[_stateq_idx(i, k, ixcldice, psetcols, pver)] * pdel / gravit

    if ixrain > 1 and ixsnow > 1:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                pdel = state_pdel[_state2_idx(i, k, psetcols)]
                wl[i - 1] = wl[i - 1] + state_q[_stateq_idx(i, k, ixrain, psetcols, pver)] * pdel / gravit
                wi[i - 1] = wi[i - 1] + state_q[_stateq_idx(i, k, ixsnow, psetcols, pver)] * pdel / gravit

    for i in range(1, ncol + 1):
        state_te_ini[i - 1] = se[i - 1] + ke[i - 1] + (latvap + latice) * wv[i - 1] + latice * wl[i - 1]
        state_tw_ini[i - 1] = wv[i - 1] + wl[i - 1] + wi[i - 1]


@export
def check_energy_chng_codon(
    ncol: int,
    pver: int,
    psetcols: int,
    latvap: float,
    latice: float,
    gravit: float,
    ixcldliq: int,
    ixcldice: int,
    ixrain: int,
    ixsnow: int,
    state_u_p: cobj,
    state_v_p: cobj,
    state_s_p: cobj,
    state_q_p: cobj,
    state_pdel_p: cobj,
    flx_vap_p: cobj,
    flx_cnd_p: cobj,
    flx_ice_p: cobj,
    flx_sen_p: cobj,
    ke_p: cobj,
    se_p: cobj,
    wv_p: cobj,
    wl_p: cobj,
    wi_p: cobj,
    tend_te_tnd_p: cobj,
    tend_tw_tnd_p: cobj,
    state_te_cur_p: cobj,
    state_tw_cur_p: cobj,
):
    state_u = Ptr[float](state_u_p)
    state_v = Ptr[float](state_v_p)
    state_s = Ptr[float](state_s_p)
    state_q = Ptr[float](state_q_p)
    state_pdel = Ptr[float](state_pdel_p)
    flx_vap = Ptr[float](flx_vap_p)
    flx_cnd = Ptr[float](flx_cnd_p)
    flx_ice = Ptr[float](flx_ice_p)
    flx_sen = Ptr[float](flx_sen_p)
    ke = Ptr[float](ke_p)
    se = Ptr[float](se_p)
    wv = Ptr[float](wv_p)
    wl = Ptr[float](wl_p)
    wi = Ptr[float](wi_p)
    tend_te_tnd = Ptr[float](tend_te_tnd_p)
    tend_tw_tnd = Ptr[float](tend_tw_tnd_p)
    state_te_cur = Ptr[float](state_te_cur_p)
    state_tw_cur = Ptr[float](state_tw_cur_p)

    for i in range(1, ncol + 1):
        ke[i - 1] = 0.0
        se[i - 1] = 0.0
        wv[i - 1] = 0.0
        wl[i - 1] = 0.0
        wi[i - 1] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            pdel = state_pdel[_state2_idx(i, k, psetcols)]
            u = state_u[_state2_idx(i, k, psetcols)]
            v = state_v[_state2_idx(i, k, psetcols)]
            ke[i - 1] = ke[i - 1] + 0.5 * (u * u + v * v) * pdel / gravit
            se[i - 1] = se[i - 1] + state_s[_state2_idx(i, k, psetcols)] * pdel / gravit
            wv[i - 1] = wv[i - 1] + state_q[_stateq_idx(i, k, 1, psetcols, pver)] * pdel / gravit

    if ixcldliq > 1 and ixcldice > 1:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                pdel = state_pdel[_state2_idx(i, k, psetcols)]
                wl[i - 1] = wl[i - 1] + state_q[_stateq_idx(i, k, ixcldliq, psetcols, pver)] * pdel / gravit
                wi[i - 1] = wi[i - 1] + state_q[_stateq_idx(i, k, ixcldice, psetcols, pver)] * pdel / gravit

    if ixrain > 1 and ixsnow > 1:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                pdel = state_pdel[_state2_idx(i, k, psetcols)]
                wl[i - 1] = wl[i - 1] + state_q[_stateq_idx(i, k, ixrain, psetcols, pver)] * pdel / gravit
                wi[i - 1] = wi[i - 1] + state_q[_stateq_idx(i, k, ixsnow, psetcols, pver)] * pdel / gravit

    for i in range(1, ncol + 1):
        te = se[i - 1] + ke[i - 1] + (latvap + latice) * wv[i - 1] + latice * wl[i - 1]
        tw = wv[i - 1] + wl[i - 1] + wi[i - 1]
        tend_te_tnd[i - 1] = (
            tend_te_tnd[i - 1]
            + flx_vap[i - 1] * (latvap + latice)
            - (flx_cnd[i - 1] - flx_ice[i - 1]) * 1000.0 * latice
            + flx_sen[i - 1]
        )
        tend_tw_tnd[i - 1] = tend_tw_tnd[i - 1] + flx_vap[i - 1] - flx_cnd[i - 1] * 1000.0
        state_te_cur[i - 1] = te
        state_tw_cur[i - 1] = tw


@export
def check_tracers_init_codon():
    return


@export
def check_tracers_chng_codon():
    return


@export
def check_tracers_batch_init_codon():
    check_tracers_init_codon()


@export
def check_tracers_batch_chng_codon():
    check_tracers_chng_codon()


@export
def check_energy_gmean_fill_codon(
    ncol: int,
    state_te_ini_p: cobj,
    teout_p: cobj,
    pint_surf_p: cobj,
    te1_p: cobj,
    te2_p: cobj,
    te3_p: cobj,
):
    state_te_ini = Ptr[float](state_te_ini_p)
    teout = Ptr[float](teout_p)
    pint_surf = Ptr[float](pint_surf_p)
    te1 = Ptr[float](te1_p)
    te2 = Ptr[float](te2_p)
    te3 = Ptr[float](te3_p)

    for i in range(1, ncol + 1):
        te1[i - 1] = state_te_ini[i - 1]
        te2[i - 1] = teout[i - 1]
        te3[i - 1] = pint_surf[i - 1]


@export
def check_energy_batch_timestep_init_codon(
    ncol: int,
    pver: int,
    psetcols: int,
    pcnst: int,
    latvap: float,
    latice: float,
    gravit: float,
    ixcldliq: int,
    ixcldice: int,
    ixrain: int,
    ixsnow: int,
    state_u_p: cobj,
    state_v_p: cobj,
    state_s_p: cobj,
    state_q_p: cobj,
    state_pdel_p: cobj,
    ke_p: cobj,
    se_p: cobj,
    wv_p: cobj,
    wl_p: cobj,
    wi_p: cobj,
    state_te_ini_p: cobj,
    state_tw_ini_p: cobj,
):
    check_energy_timestep_init_codon(
        ncol,
        pver,
        psetcols,
        pcnst,
        latvap,
        latice,
        gravit,
        ixcldliq,
        ixcldice,
        ixrain,
        ixsnow,
        state_u_p,
        state_v_p,
        state_s_p,
        state_q_p,
        state_pdel_p,
        ke_p,
        se_p,
        wv_p,
        wl_p,
        wi_p,
        state_te_ini_p,
        state_tw_ini_p,
    )


@export
def check_energy_batch_chng_codon(
    ncol: int,
    pver: int,
    psetcols: int,
    latvap: float,
    latice: float,
    gravit: float,
    ixcldliq: int,
    ixcldice: int,
    ixrain: int,
    ixsnow: int,
    state_u_p: cobj,
    state_v_p: cobj,
    state_s_p: cobj,
    state_q_p: cobj,
    state_pdel_p: cobj,
    flx_vap_p: cobj,
    flx_cnd_p: cobj,
    flx_ice_p: cobj,
    flx_sen_p: cobj,
    ke_p: cobj,
    se_p: cobj,
    wv_p: cobj,
    wl_p: cobj,
    wi_p: cobj,
    tend_te_tnd_p: cobj,
    tend_tw_tnd_p: cobj,
    state_te_cur_p: cobj,
    state_tw_cur_p: cobj,
):
    check_energy_chng_codon(
        ncol,
        pver,
        psetcols,
        latvap,
        latice,
        gravit,
        ixcldliq,
        ixcldice,
        ixrain,
        ixsnow,
        state_u_p,
        state_v_p,
        state_s_p,
        state_q_p,
        state_pdel_p,
        flx_vap_p,
        flx_cnd_p,
        flx_ice_p,
        flx_sen_p,
        ke_p,
        se_p,
        wv_p,
        wl_p,
        wi_p,
        tend_te_tnd_p,
        tend_tw_tnd_p,
        state_te_cur_p,
        state_tw_cur_p,
    )


@export
def check_energy_batch_gmean_fill_codon(
    ncol: int,
    state_te_ini_p: cobj,
    teout_p: cobj,
    pint_surf_p: cobj,
    te1_p: cobj,
    te2_p: cobj,
    te3_p: cobj,
):
    check_energy_gmean_fill_codon(ncol, state_te_ini_p, teout_p, pint_surf_p, te1_p, te2_p, te3_p)


@export
def check_energy_batch_fix_codon(
    ncol: int,
    pcols: int,
    pver: int,
    psetcols: int,
    heat_glob: float,
    gravit: float,
    state_pint_p: cobj,
    ptend_s_p: cobj,
    eshflx_p: cobj,
):
    check_energy_fix_codon(ncol, pcols, pver, psetcols, heat_glob, gravit, state_pint_p, ptend_s_p, eshflx_p)


@export
def check_energy_batch_dispatch_codon(
    stage: int,
    ncol: int,
    pver: int,
    pcols: int,
    psetcols: int,
    pcnst: int,
    ixcldliq: int,
    ixcldice: int,
    ixrain: int,
    ixsnow: int,
    scalar1: float,
    scalar2: float,
    scalar3: float,
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
):
    if stage == 1:
        check_energy_batch_timestep_init_codon(
            ncol,
            pver,
            psetcols,
            pcnst,
            scalar1,
            scalar2,
            scalar3,
            ixcldliq,
            ixcldice,
            ixrain,
            ixsnow,
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
        )
    elif stage == 2:
        check_energy_batch_chng_codon(
            ncol,
            pver,
            psetcols,
            scalar1,
            scalar2,
            scalar3,
            ixcldliq,
            ixcldice,
            ixrain,
            ixsnow,
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
        check_energy_batch_gmean_fill_codon(ncol, p1, p2, p3, p4, p5, p6)
    elif stage == 4:
        check_energy_batch_fix_codon(ncol, pcols, pver, psetcols, scalar1, scalar2, p1, p2, p3)
    elif stage == 5:
        check_tracers_batch_init_codon()
    elif stage == 6:
        check_tracers_batch_chng_codon()
