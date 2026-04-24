from math import sqrt


@export
def vertical_diffusion_ts_init_codon():
    return


@export
def vertical_diffusion_tend_select_branches_codon(
    do_tms: int,
    do_molec_diff: int,
    use_diag_tke: int,
    use_hb_family: int,
    shallow_unicon: int,
    prog_modal_aero: int,
    do_pseudocon_diff: int,
    diff_cnsrv_mass_check: int,
    waccmx_special: int,
    branch_mask_p: cobj,
):
    branch_mask = Ptr[int](branch_mask_p)

    mask = 0
    if do_tms != 0:
        mask |= 1
    if do_molec_diff != 0:
        mask |= 2
    if use_diag_tke != 0:
        mask |= 4
    if use_hb_family != 0:
        mask |= 8
    if shallow_unicon != 0:
        mask |= 16
    if prog_modal_aero != 0:
        mask |= 32
    if do_pseudocon_diff != 0:
        mask |= 64
    if diff_cnsrv_mass_check != 0:
        mask |= 128
    if waccmx_special != 0:
        mask |= 256

    branch_mask[0] = mask


@inline
def _idx2(i: int, k: int, ld1: int) -> int:
    """Fortran arrays declared as (ld1, nlev)."""
    return (i - 1) + (k - 1) * ld1


@inline
def _idx3(i: int, k: int, m: int, ld1: int, ld2: int) -> int:
    """Fortran arrays declared as (ld1, ld2, nconst)."""
    return (i - 1) + (k - 1) * ld1 + (m - 1) * ld1 * ld2


@export
def vertical_diffusion_flux_diag_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    ixcldliq: int,
    ixcldice: int,
    latvap: float,
    latice: float,
    zvir: float,
    rair: float,
    gravit: float,
    cpair: float,
    q_tmp_p: cobj,
    s_tmp_p: cobj,
    u_tmp_p: cobj,
    v_tmp_p: cobj,
    pint_p: cobj,
    zi_p: cobj,
    zm_p: cobj,
    cflx_p: cobj,
    kvh_p: cobj,
    kvm_p: cobj,
    cgs_p: cobj,
    cgh_p: cobj,
    shflx_p: cobj,
    tautotx_p: cobj,
    tautoty_p: cobj,
    sl_p: cobj,
    qt_p: cobj,
    slv_p: cobj,
    slflx_p: cobj,
    qtflx_p: cobj,
    uflx_p: cobj,
    vflx_p: cobj,
    slflx_cg_p: cobj,
    qtflx_cg_p: cobj,
    uflx_cg_p: cobj,
    vflx_cg_p: cobj,
):
    q_tmp = Ptr[float](q_tmp_p)
    s_tmp = Ptr[float](s_tmp_p)
    u_tmp = Ptr[float](u_tmp_p)
    v_tmp = Ptr[float](v_tmp_p)
    pint = Ptr[float](pint_p)
    zi = Ptr[float](zi_p)
    zm = Ptr[float](zm_p)
    cflx = Ptr[float](cflx_p)
    kvh = Ptr[float](kvh_p)
    kvm = Ptr[float](kvm_p)
    cgs = Ptr[float](cgs_p)
    cgh = Ptr[float](cgh_p)
    shflx = Ptr[float](shflx_p)
    tautotx = Ptr[float](tautotx_p)
    tautoty = Ptr[float](tautoty_p)
    sl = Ptr[float](sl_p)
    qt = Ptr[float](qt_p)
    slv = Ptr[float](slv_p)
    slflx = Ptr[float](slflx_p)
    qtflx = Ptr[float](qtflx_p)
    uflx = Ptr[float](uflx_p)
    vflx = Ptr[float](vflx_p)
    slflx_cg = Ptr[float](slflx_cg_p)
    qtflx_cg = Ptr[float](qtflx_cg_p)
    uflx_cg = Ptr[float](uflx_cg_p)
    vflx_cg = Ptr[float](vflx_cg_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            sl[_idx2(i, k, pcols)] = (
                s_tmp[_idx2(i, k, pcols)]
                - latvap * q_tmp[_idx3(i, k, ixcldliq, pcols, pver)]
                - (latvap + latice) * q_tmp[_idx3(i, k, ixcldice, pcols, pver)]
            )

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            qt[_idx2(i, k, pcols)] = (
                q_tmp[_idx3(i, k, 1, pcols, pver)]
                + q_tmp[_idx3(i, k, ixcldliq, pcols, pver)]
                + q_tmp[_idx3(i, k, ixcldice, pcols, pver)]
            )

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            slv[_idx2(i, k, pcols)] = sl[_idx2(i, k, pcols)] * (
                1.0 + zvir * qt[_idx2(i, k, pcols)]
            )

    for i in range(1, ncol + 1):
        slflx[_idx2(i, 1, pcols)] = 0.0
        qtflx[_idx2(i, 1, pcols)] = 0.0
        uflx[_idx2(i, 1, pcols)] = 0.0
        vflx[_idx2(i, 1, pcols)] = 0.0
        slflx_cg[_idx2(i, 1, pcols)] = 0.0
        qtflx_cg[_idx2(i, 1, pcols)] = 0.0
        uflx_cg[_idx2(i, 1, pcols)] = 0.0
        vflx_cg[_idx2(i, 1, pcols)] = 0.0

    for k in range(2, pver + 1):
        for i in range(1, ncol + 1):
            rhoair = pint[_idx2(i, k, pcols)] / (
                rair
                * (
                    (
                        0.5
                        * (
                            slv[_idx2(i, k, pcols)]
                            + slv[_idx2(i, k - 1, pcols)]
                        )
                        - gravit * zi[_idx2(i, k, pcols)]
                    )
                    / cpair
                )
            )
            slflx[_idx2(i, k, pcols)] = kvh[_idx2(i, k, pcols)] * (
                -rhoair
                * (
                    sl[_idx2(i, k - 1, pcols)]
                    - sl[_idx2(i, k, pcols)]
                )
                / (
                    zm[_idx2(i, k - 1, pcols)]
                    - zm[_idx2(i, k, pcols)]
                )
                + cgh[_idx2(i, k, pcols)]
            )
            qtflx[_idx2(i, k, pcols)] = kvh[_idx2(i, k, pcols)] * (
                -rhoair
                * (
                    qt[_idx2(i, k - 1, pcols)]
                    - qt[_idx2(i, k, pcols)]
                )
                / (
                    zm[_idx2(i, k - 1, pcols)]
                    - zm[_idx2(i, k, pcols)]
                )
                + rhoair
                * (
                    cflx[_idx2(i, 1, pcols)]
                    + cflx[_idx2(i, ixcldliq, pcols)]
                    + cflx[_idx2(i, ixcldice, pcols)]
                )
                * cgs[_idx2(i, k, pcols)]
            )
            uflx[_idx2(i, k, pcols)] = kvm[_idx2(i, k, pcols)] * (
                -rhoair
                * (
                    u_tmp[_idx2(i, k - 1, pcols)]
                    - u_tmp[_idx2(i, k, pcols)]
                )
                / (
                    zm[_idx2(i, k - 1, pcols)]
                    - zm[_idx2(i, k, pcols)]
                )
            )
            vflx[_idx2(i, k, pcols)] = kvm[_idx2(i, k, pcols)] * (
                -rhoair
                * (
                    v_tmp[_idx2(i, k - 1, pcols)]
                    - v_tmp[_idx2(i, k, pcols)]
                )
                / (
                    zm[_idx2(i, k - 1, pcols)]
                    - zm[_idx2(i, k, pcols)]
                )
            )
            slflx_cg[_idx2(i, k, pcols)] = (
                kvh[_idx2(i, k, pcols)] * cgh[_idx2(i, k, pcols)]
            )
            qtflx_cg[_idx2(i, k, pcols)] = (
                kvh[_idx2(i, k, pcols)]
                * rhoair
                * (
                    cflx[_idx2(i, 1, pcols)]
                    + cflx[_idx2(i, ixcldliq, pcols)]
                    + cflx[_idx2(i, ixcldice, pcols)]
                )
                * cgs[_idx2(i, k, pcols)]
            )
            uflx_cg[_idx2(i, k, pcols)] = 0.0
            vflx_cg[_idx2(i, k, pcols)] = 0.0

    for i in range(1, ncol + 1):
        slflx[_idx2(i, pverp, pcols)] = shflx[i - 1]
        qtflx[_idx2(i, pverp, pcols)] = cflx[_idx2(i, 1, pcols)]
        uflx[_idx2(i, pverp, pcols)] = tautotx[i - 1]
        vflx[_idx2(i, pverp, pcols)] = tautoty[i - 1]
        slflx_cg[_idx2(i, pverp, pcols)] = 0.0
        qtflx_cg[_idx2(i, pverp, pcols)] = 0.0
        uflx_cg[_idx2(i, pverp, pcols)] = 0.0
        vflx_cg[_idx2(i, pverp, pcols)] = 0.0


@export
def vertical_diffusion_ptend_core_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    psetcols: int,
    rztodt: float,
    q_tmp_p: cobj,
    s_tmp_p: cobj,
    u_tmp_p: cobj,
    v_tmp_p: cobj,
    state_q_p: cobj,
    state_s_p: cobj,
    state_u_p: cobj,
    state_v_p: cobj,
    sl_p: cobj,
    qt_p: cobj,
    sl_prePBL_p: cobj,
    qt_prePBL_p: cobj,
    ptend_q_p: cobj,
    ptend_s_p: cobj,
    ptend_u_p: cobj,
    ptend_v_p: cobj,
    slten_p: cobj,
    qtten_p: cobj,
):
    q_tmp = Ptr[float](q_tmp_p)
    s_tmp = Ptr[float](s_tmp_p)
    u_tmp = Ptr[float](u_tmp_p)
    v_tmp = Ptr[float](v_tmp_p)
    state_q = Ptr[float](state_q_p)
    state_s = Ptr[float](state_s_p)
    state_u = Ptr[float](state_u_p)
    state_v = Ptr[float](state_v_p)
    sl = Ptr[float](sl_p)
    qt = Ptr[float](qt_p)
    sl_prePBL = Ptr[float](sl_prePBL_p)
    qt_prePBL = Ptr[float](qt_prePBL_p)
    ptend_q = Ptr[float](ptend_q_p)
    ptend_s = Ptr[float](ptend_s_p)
    ptend_u = Ptr[float](ptend_u_p)
    ptend_v = Ptr[float](ptend_v_p)
    slten = Ptr[float](slten_p)
    qtten = Ptr[float](qtten_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            ptend_s[_idx2(i, k, psetcols)] = (
                s_tmp[_idx2(i, k, pcols)] - state_s[_idx2(i, k, pcols)]
            ) * rztodt
            ptend_u[_idx2(i, k, psetcols)] = (
                u_tmp[_idx2(i, k, pcols)] - state_u[_idx2(i, k, pcols)]
            ) * rztodt
            ptend_v[_idx2(i, k, psetcols)] = (
                v_tmp[_idx2(i, k, pcols)] - state_v[_idx2(i, k, pcols)]
            ) * rztodt
            slten[_idx2(i, k, pcols)] = (
                sl[_idx2(i, k, pcols)] - sl_prePBL[_idx2(i, k, pcols)]
            ) * rztodt
            qtten[_idx2(i, k, pcols)] = (
                qt[_idx2(i, k, pcols)] - qt_prePBL[_idx2(i, k, pcols)]
            ) * rztodt
            for m in range(1, pcnst + 1):
                ptend_q[_idx3(i, k, m, psetcols, pver)] = (
                    q_tmp[_idx3(i, k, m, pcols, pver)]
                    - state_q[_idx3(i, k, m, pcols, pver)]
                ) * rztodt


@export
def vertical_diffusion_pre_pbl_diag_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    ixcldliq: int,
    ixcldice: int,
    latvap: float,
    latice: float,
    zvir: float,
    state_q_p: cobj,
    state_s_p: cobj,
    state_u_p: cobj,
    state_v_p: cobj,
    q_tmp_p: cobj,
    s_tmp_p: cobj,
    u_tmp_p: cobj,
    v_tmp_p: cobj,
    sl_prePBL_p: cobj,
    qt_prePBL_p: cobj,
    slv_prePBL_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    state_s = Ptr[float](state_s_p)
    state_u = Ptr[float](state_u_p)
    state_v = Ptr[float](state_v_p)
    q_tmp = Ptr[float](q_tmp_p)
    s_tmp = Ptr[float](s_tmp_p)
    u_tmp = Ptr[float](u_tmp_p)
    v_tmp = Ptr[float](v_tmp_p)
    sl_prePBL = Ptr[float](sl_prePBL_p)
    qt_prePBL = Ptr[float](qt_prePBL_p)
    slv_prePBL = Ptr[float](slv_prePBL_p)

    for m in range(1, pcnst + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                q_tmp[_idx3(i, k, m, pcols, pver)] = state_q[
                    _idx3(i, k, m, pcols, pver)
                ]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            s_tmp[_idx2(i, k, pcols)] = state_s[_idx2(i, k, pcols)]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            u_tmp[_idx2(i, k, pcols)] = state_u[_idx2(i, k, pcols)]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            v_tmp[_idx2(i, k, pcols)] = state_v[_idx2(i, k, pcols)]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            sl_prePBL[_idx2(i, k, pcols)] = (
                s_tmp[_idx2(i, k, pcols)]
                - latvap * q_tmp[_idx3(i, k, ixcldliq, pcols, pver)]
                - (latvap + latice) * q_tmp[_idx3(i, k, ixcldice, pcols, pver)]
            )

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            qt_prePBL[_idx2(i, k, pcols)] = (
                q_tmp[_idx3(i, k, 1, pcols, pver)]
                + q_tmp[_idx3(i, k, ixcldliq, pcols, pver)]
                + q_tmp[_idx3(i, k, ixcldice, pcols, pver)]
            )

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            slv_prePBL[_idx2(i, k, pcols)] = sl_prePBL[_idx2(i, k, pcols)] * (
                1.0 + zvir * qt_prePBL[_idx2(i, k, pcols)]
            )


@export
def vertical_diffusion_post_pbl_state_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    psetcols: int,
    ixcldliq: int,
    ixcldice: int,
    ztodt: float,
    gravit: float,
    cpair: float,
    state_q_p: cobj,
    state_s_p: cobj,
    state_u_p: cobj,
    state_v_p: cobj,
    state_zm_p: cobj,
    ptend_q_p: cobj,
    ptend_s_p: cobj,
    ptend_u_p: cobj,
    ptend_v_p: cobj,
    qv_aft_PBL_p: cobj,
    ql_aft_PBL_p: cobj,
    qi_aft_PBL_p: cobj,
    s_aft_PBL_p: cobj,
    t_aftPBL_p: cobj,
    u_aft_PBL_p: cobj,
    v_aft_PBL_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    state_s = Ptr[float](state_s_p)
    state_u = Ptr[float](state_u_p)
    state_v = Ptr[float](state_v_p)
    state_zm = Ptr[float](state_zm_p)
    ptend_q = Ptr[float](ptend_q_p)
    ptend_s = Ptr[float](ptend_s_p)
    ptend_u = Ptr[float](ptend_u_p)
    ptend_v = Ptr[float](ptend_v_p)
    qv_aft_PBL = Ptr[float](qv_aft_PBL_p)
    ql_aft_PBL = Ptr[float](ql_aft_PBL_p)
    qi_aft_PBL = Ptr[float](qi_aft_PBL_p)
    s_aft_PBL = Ptr[float](s_aft_PBL_p)
    t_aftPBL = Ptr[float](t_aftPBL_p)
    u_aft_PBL = Ptr[float](u_aft_PBL_p)
    v_aft_PBL = Ptr[float](v_aft_PBL_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            qv_aft_PBL[_idx2(i, k, pcols)] = (
                state_q[_idx3(i, k, 1, pcols, pver)]
                + ptend_q[_idx3(i, k, 1, psetcols, pver)] * ztodt
            )

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            ql_aft_PBL[_idx2(i, k, pcols)] = (
                state_q[_idx3(i, k, ixcldliq, pcols, pver)]
                + ptend_q[_idx3(i, k, ixcldliq, psetcols, pver)] * ztodt
            )

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            qi_aft_PBL[_idx2(i, k, pcols)] = (
                state_q[_idx3(i, k, ixcldice, pcols, pver)]
                + ptend_q[_idx3(i, k, ixcldice, psetcols, pver)] * ztodt
            )

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            s_aft_PBL[_idx2(i, k, pcols)] = (
                state_s[_idx2(i, k, pcols)]
                + ptend_s[_idx2(i, k, psetcols)] * ztodt
            )

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            t_aftPBL[_idx2(i, k, pcols)] = (
                s_aft_PBL[_idx2(i, k, pcols)]
                - gravit * state_zm[_idx2(i, k, pcols)]
            ) / cpair

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            u_aft_PBL[_idx2(i, k, pcols)] = (
                state_u[_idx2(i, k, pcols)]
                + ptend_u[_idx2(i, k, psetcols)] * ztodt
            )

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            v_aft_PBL[_idx2(i, k, pcols)] = (
                state_v[_idx2(i, k, pcols)]
                + ptend_v[_idx2(i, k, psetcols)] * ztodt
            )


@export
def vertical_diffusion_modal_aero_flux_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    pmam_ncnst: int,
    ztodt: float,
    gravit: float,
    state_rpdel_p: cobj,
    pmam_cnst_idx_p: cobj,
    cflx_p: cobj,
    q_tmp_p: cobj,
):
    state_rpdel = Ptr[float](state_rpdel_p)
    pmam_cnst_idx = Ptr[int](pmam_cnst_idx_p)
    cflx = Ptr[float](cflx_p)
    q_tmp = Ptr[float](q_tmp_p)

    for i in range(1, ncol + 1):
        tmp1 = ztodt * gravit * state_rpdel[_idx2(i, pver, pcols)]
        for m in range(1, pmam_ncnst + 1):
            l = pmam_cnst_idx[m - 1]
            q_tmp[_idx3(i, pver, l, pcols, pver)] = (
                q_tmp[_idx3(i, pver, l, pcols, pver)]
                + tmp1 * cflx[_idx2(i, l, pcols)]
            )


@export
def vertical_diffusion_pre_qsat_rh_codon(
    ncol: int,
    pcols: int,
    pver: int,
    state_q_p: cobj,
    ftem_p: cobj,
    ftem_prePBL_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    ftem = Ptr[float](ftem_p)
    ftem_prePBL = Ptr[float](ftem_prePBL_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            ftem_prePBL[_idx2(i, k, pcols)] = (
                state_q[_idx3(i, k, 1, pcols, pver)] / ftem[_idx2(i, k, pcols)] * 100.0
            )


@export
def austausch_atm_codon(
    ncol: int,
    pcols: int,
    pver: int,
    ntop_turb: int,
    nbot_turb: int,
    zkmin: float,
    ri_p: cobj,
    s2_p: cobj,
    ml2_p: cobj,
    kvf_p: cobj,
):
    ri = Ptr[float](ri_p)
    s2 = Ptr[float](s2_p)
    ml2 = Ptr[float](ml2_p)
    kvf = Ptr[float](kvf_p)

    for k in range(1, pver + 2):
        for i in range(1, ncol + 1):
            kvf[_idx2(i, k, pcols)] = 0.0

    for k in range(ntop_turb + 1, nbot_turb + 1):
        for i in range(1, ncol + 1):
            if ri[_idx2(i, k, pcols)] < 0.0:
                fofri = sqrt(max(1.0 - 18.0 * ri[_idx2(i, k, pcols)], 0.0))
            else:
                fofri = 1.0 / (
                    1.0
                    + 10.0
                    * ri[_idx2(i, k, pcols)]
                    * (1.0 + 8.0 * ri[_idx2(i, k, pcols)])
                )
            kvn = ml2[k - 1] * sqrt(s2[_idx2(i, k, pcols)])
            kvf[_idx2(i, k, pcols)] = max(zkmin, kvn * fofri)


@export
def eddy_diff_surface_stress_diag_codon(
    ncol: int,
    pcols: int,
    pver: int,
    rair: float,
    ustar_min: float,
    tfd_p: cobj,
    pmid_p: cobj,
    taux_p: cobj,
    tauy_p: cobj,
    ksrftms_p: cobj,
    ufd_p: cobj,
    vfd_p: cobj,
    rrho_p: cobj,
    ustar_p: cobj,
    minpblh_p: cobj,
):
    tfd = Ptr[float](tfd_p)
    pmid = Ptr[float](pmid_p)
    taux = Ptr[float](taux_p)
    tauy = Ptr[float](tauy_p)
    ksrftms = Ptr[float](ksrftms_p)
    ufd = Ptr[float](ufd_p)
    vfd = Ptr[float](vfd_p)
    rrho = Ptr[float](rrho_p)
    ustar = Ptr[float](ustar_p)
    minpblh = Ptr[float](minpblh_p)

    for i in range(1, ncol + 1):
        taux_eff = taux[i - 1] - ksrftms[i - 1] * ufd[_idx2(i, pver, pcols)]
        tauy_eff = tauy[i - 1] - ksrftms[i - 1] * vfd[_idx2(i, pver, pcols)]
        rrho[i - 1] = rair * tfd[_idx2(i, pver, pcols)] / pmid[_idx2(i, pver, pcols)]
        ustar[i - 1] = max(
            sqrt(sqrt(taux_eff * taux_eff + tauy_eff * tauy_eff) * rrho[i - 1]),
            ustar_min,
        )
        minpblh[i - 1] = 100.0 * ustar[i - 1]


@export
def eddy_diff_kv_init_codon(
    ncol: int,
    pcols: int,
    pver: int,
    iturb: int,
    kvinit: int,
    use_kvf: int,
    kvf_p: cobj,
    kvh_in_p: cobj,
    kvm_in_p: cobj,
    kvh_out_p: cobj,
    kvm_out_p: cobj,
    kvh_p: cobj,
    kvm_p: cobj,
):
    kvf = Ptr[float](kvf_p)
    kvh_in = Ptr[float](kvh_in_p)
    kvm_in = Ptr[float](kvm_in_p)
    kvh_out = Ptr[float](kvh_out_p)
    kvm_out = Ptr[float](kvm_out_p)
    kvh = Ptr[float](kvh_p)
    kvm = Ptr[float](kvm_p)

    if iturb == 1:
        if kvinit != 0:
            if use_kvf != 0:
                for k in range(1, pver + 2):
                    for i in range(1, ncol + 1):
                        kvh[_idx2(i, k, pcols)] = kvf[_idx2(i, k, pcols)]
                        kvm[_idx2(i, k, pcols)] = kvf[_idx2(i, k, pcols)]
            else:
                for k in range(1, pver + 2):
                    for i in range(1, ncol + 1):
                        kvh[_idx2(i, k, pcols)] = 0.0
                        kvm[_idx2(i, k, pcols)] = 0.0
        else:
            for k in range(1, pver + 2):
                for i in range(1, ncol + 1):
                    kvh[_idx2(i, k, pcols)] = kvh_in[_idx2(i, k, pcols)]
                    kvm[_idx2(i, k, pcols)] = kvm_in[_idx2(i, k, pcols)]
    else:
            for k in range(1, pver + 2):
                for i in range(1, ncol + 1):
                    kvh[_idx2(i, k, pcols)] = kvh_out[_idx2(i, k, pcols)]
                    kvm[_idx2(i, k, pcols)] = kvm_out[_idx2(i, k, pcols)]


@export
def eddy_diff_error_pbl_codon(
    ncol: int,
    pcols: int,
    pver: int,
    kvh_p: cobj,
    kvh_out_p: cobj,
    errorPBL_p: cobj,
):
    kvh = Ptr[float](kvh_p)
    kvh_out = Ptr[float](kvh_out_p)
    errorPBL = Ptr[float](errorPBL_p)

    for i in range(1, ncol + 1):
        errorPBL[i - 1] = 0.0
        for k in range(1, pver + 1):
            errorPBL[i - 1] = errorPBL[i - 1] + (
                kvh[_idx2(i, k, pcols)] - kvh_out[_idx2(i, k, pcols)]
            ) * (
                kvh[_idx2(i, k, pcols)] - kvh_out[_idx2(i, k, pcols)]
            )
        errorPBL[i - 1] = sqrt(errorPBL[i - 1] / pver)


@export
def eddy_diff_kv_relax_codon(
    ncol: int,
    pcols: int,
    pver: int,
    lambda_v: float,
    kvm_p: cobj,
    kvh_p: cobj,
    kvm_out_p: cobj,
    kvh_out_p: cobj,
):
    kvm = Ptr[float](kvm_p)
    kvh = Ptr[float](kvh_p)
    kvm_out = Ptr[float](kvm_out_p)
    kvh_out = Ptr[float](kvh_out_p)

    for k in range(1, pver + 2):
        for i in range(1, ncol + 1):
            kvm_out[_idx2(i, k, pcols)] = (
                lambda_v * kvm_out[_idx2(i, k, pcols)]
                + (1.0 - lambda_v) * kvm[_idx2(i, k, pcols)]
            )

    for k in range(1, pver + 2):
        for i in range(1, ncol + 1):
            kvh_out[_idx2(i, k, pcols)] = (
                lambda_v * kvh_out[_idx2(i, k, pcols)]
                + (1.0 - lambda_v) * kvh[_idx2(i, k, pcols)]
            )


@export
def vertical_diffusion_obklen_diag_codon(
    ncol: int,
    pcols: int,
    pver: int,
    cpair: float,
    zvir: float,
    gravit: float,
    karman: float,
    state_t_p: cobj,
    state_exner_p: cobj,
    state_q_p: cobj,
    cflx_p: cobj,
    shflx_p: cobj,
    rrho_p: cobj,
    ustar_p: cobj,
    th_p: cobj,
    thvs_p: cobj,
    khfs_p: cobj,
    kqfs_p: cobj,
    kbfs_p: cobj,
    obklen_p: cobj,
):
    state_t = Ptr[float](state_t_p)
    state_exner = Ptr[float](state_exner_p)
    state_q = Ptr[float](state_q_p)
    cflx = Ptr[float](cflx_p)
    shflx = Ptr[float](shflx_p)
    rrho = Ptr[float](rrho_p)
    ustar = Ptr[float](ustar_p)
    th = Ptr[float](th_p)
    thvs = Ptr[float](thvs_p)
    khfs = Ptr[float](khfs_p)
    kqfs = Ptr[float](kqfs_p)
    kbfs = Ptr[float](kbfs_p)
    obklen = Ptr[float](obklen_p)

    for i in range(1, ncol + 1):
        th[_idx2(i, pver, pcols)] = (
            state_t[_idx2(i, pver, pcols)] * state_exner[_idx2(i, pver, pcols)]
        )
        thvs[i - 1] = th[_idx2(i, pver, pcols)] * (
            1.0 + zvir * state_q[_idx3(i, pver, 1, pcols, pver)]
        )
        khfs[i - 1] = shflx[i - 1] * rrho[i - 1] / cpair
        kqfs[i - 1] = cflx[_idx2(i, 1, pcols)] * rrho[i - 1]
        kbfs[i - 1] = khfs[i - 1] + zvir * th[_idx2(i, pver, pcols)] * kqfs[i - 1]
        if kbfs[i - 1] < 0.0:
            signed_eps = -1.0e-10
        else:
            signed_eps = 1.0e-10
        obklen[i - 1] = (
            -thvs[i - 1]
            * (ustar[i - 1] * ustar[i - 1] * ustar[i - 1])
            / (gravit * karman * (kbfs[i - 1] + signed_eps))
        )


@export
def vertical_diffusion_post_qsat_diag_codon(
    ncol: int,
    pcols: int,
    pver: int,
    rztodt: float,
    state_t_p: cobj,
    qv_aft_PBL_p: cobj,
    ftem_prePBL_p: cobj,
    t_aftPBL_p: cobj,
    ftem_p: cobj,
    ftem_aftPBL_p: cobj,
    tten_p: cobj,
    rhten_p: cobj,
):
    state_t = Ptr[float](state_t_p)
    qv_aft_PBL = Ptr[float](qv_aft_PBL_p)
    ftem_prePBL = Ptr[float](ftem_prePBL_p)
    t_aftPBL = Ptr[float](t_aftPBL_p)
    ftem = Ptr[float](ftem_p)
    ftem_aftPBL = Ptr[float](ftem_aftPBL_p)
    tten = Ptr[float](tten_p)
    rhten = Ptr[float](rhten_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            ftem_aftPBL[_idx2(i, k, pcols)] = (
                qv_aft_PBL[_idx2(i, k, pcols)] / ftem[_idx2(i, k, pcols)] * 100.0
            )

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            tten[_idx2(i, k, pcols)] = (
                t_aftPBL[_idx2(i, k, pcols)] - state_t[_idx2(i, k, pcols)]
            ) * rztodt

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            rhten[_idx2(i, k, pcols)] = (
                ftem_aftPBL[_idx2(i, k, pcols)] - ftem_prePBL[_idx2(i, k, pcols)]
            ) * rztodt
