from C import eddy_diff_estblf_cb(float) -> float
from C import eddy_diff_svp_to_qsat_cb(float, float) -> float
from math import acos, cos, exp, log, sqrt


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
def eddy_diff_wstar_pbl_codon(
    ncol: int,
    pcols: int,
    ncvmax: int,
    ipbl_p: cobj,
    wstar_p: cobj,
    wstarPBL_p: cobj,
):
    ipbl = Ptr[int](ipbl_p)
    wstar = Ptr[float](wstar_p)
    wstarPBL = Ptr[float](wstarPBL_p)

    for i in range(1, ncol + 1):
        if ipbl[i - 1] == 1:
            wstarPBL[i - 1] = max(0.0, wstar[_idx2(i, 1, pcols)])
        else:
            wstarPBL[i - 1] = 0.0


@export
def eddy_diff_caleddy_init_codon(
    ncol: int,
    pcols: int,
    pver: int,
    qrlzero_mode: int,
    cldeff_mode: int,
    tkes_mode: int,
    use_kvf_mode: int,
    qmin: float,
    vk: float,
    ql_p: cobj,
    qrlin_p: cobj,
    cld_p: cobj,
    kvf_p: cobj,
    kvh_in_p: cobj,
    kvm_in_p: cobj,
    n2_p: cobj,
    s2_p: cobj,
    shflx_p: cobj,
    qflx_p: cobj,
    rrho_p: cobj,
    ustar_p: cobj,
    z_p: cobj,
    chu_p: cobj,
    chs_p: cobj,
    cmu_p: cobj,
    cms_p: cobj,
    sflh_p: cobj,
    qrlw_p: cobj,
    cldeff_p: cobj,
    kvh_p: cobj,
    kvm_p: cobj,
    bflxs_p: cobj,
    bprod_p: cobj,
    sprod_p: cobj,
    wcap_p: cobj,
    leng_p: cobj,
    tke_p: cobj,
    turbtype_p: cobj,
):
    ql = Ptr[float](ql_p)
    qrlin = Ptr[float](qrlin_p)
    cld = Ptr[float](cld_p)
    kvf = Ptr[float](kvf_p)
    kvh_in = Ptr[float](kvh_in_p)
    kvm_in = Ptr[float](kvm_in_p)
    n2 = Ptr[float](n2_p)
    s2 = Ptr[float](s2_p)
    shflx = Ptr[float](shflx_p)
    qflx = Ptr[float](qflx_p)
    rrho = Ptr[float](rrho_p)
    ustar = Ptr[float](ustar_p)
    z = Ptr[float](z_p)
    chu = Ptr[float](chu_p)
    chs = Ptr[float](chs_p)
    cmu = Ptr[float](cmu_p)
    cms = Ptr[float](cms_p)
    sflh = Ptr[float](sflh_p)
    qrlw = Ptr[float](qrlw_p)
    cldeff = Ptr[float](cldeff_p)
    kvh = Ptr[float](kvh_p)
    kvm = Ptr[float](kvm_p)
    bflxs = Ptr[float](bflxs_p)
    bprod = Ptr[float](bprod_p)
    sprod = Ptr[float](sprod_p)
    wcap = Ptr[float](wcap_p)
    leng = Ptr[float](leng_p)
    tke = Ptr[float](tke_p)
    turbtype = Ptr[i32](turbtype_p)

    if qrlzero_mode != 0:
        for k in range(1, pver + 1):
            for i in range(1, pcols + 1):
                qrlw[_idx2(i, k, pcols)] = 0.0
    else:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, pcols)
                qrlw[idx] = qrlin[idx]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            if cldeff_mode != 0:
                cldeff[idx] = cld[idx] * min(ql[idx] / qmin, 1.0)
            else:
                cldeff[idx] = cld[idx]

    for k in range(1, pver + 2):
        for i in range(1, pcols + 1):
            idx = _idx2(i, k, pcols)
            if use_kvf_mode != 0:
                kvh[idx] = kvf[idx]
                kvm[idx] = kvf[idx]
            else:
                kvh[idx] = 0.0
                kvm[idx] = 0.0

    for k in range(1, pver + 2):
        for i in range(1, pcols + 1):
            idx = _idx2(i, k, pcols)
            wcap[idx] = 0.0
            leng[idx] = 0.0
            tke[idx] = 0.0
            turbtype[idx] = i32(0)

    for k in range(2, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            bprod[idx] = -kvh_in[idx] * n2[idx]
            sprod[idx] = kvm_in[idx] * s2[idx]

    for i in range(1, ncol + 1):
        top_idx = _idx2(i, 1, pcols)
        surf_idx = _idx2(i, pver + 1, pcols)
        surf_layer_idx = _idx2(i, pver, pcols)

        bprod[top_idx] = 0.0
        sprod[top_idx] = 0.0

        ch = chu[surf_idx] * (1.0 - sflh[surf_layer_idx]) + chs[surf_idx] * sflh[surf_layer_idx]
        cm = cmu[surf_idx] * (1.0 - sflh[surf_layer_idx]) + cms[surf_idx] * sflh[surf_layer_idx]
        bflxs[i - 1] = ch * shflx[i - 1] * rrho[i - 1] + cm * qflx[i - 1] * rrho[i - 1]

        if tkes_mode != 0:
            bprod[surf_idx] = bflxs[i - 1]
        else:
            bprod[surf_idx] = 0.0

        sprod[surf_idx] = (ustar[i - 1] ** 3) / (vk * z[surf_layer_idx])


@export
def eddy_diff_caleddy_diaginit_codon(
    ncol: int,
    pcols: int,
    pver: int,
    ncvmax: int,
    went_p: cobj,
    wet_CL_p: cobj,
    web_CL_p: cobj,
    jtbu_CL_p: cobj,
    jbbu_CL_p: cobj,
    evhc_CL_p: cobj,
    jt2slv_CL_p: cobj,
    n2ht_CL_p: cobj,
    n2hb_CL_p: cobj,
    lwp_CL_p: cobj,
    opt_depth_CL_p: cobj,
    radinvfrac_CL_p: cobj,
    radf_CL_p: cobj,
    wstar_CL_p: cobj,
    wstar3fact_CL_p: cobj,
    ricl_p: cobj,
    ghcl_p: cobj,
    shcl_p: cobj,
    smcl_p: cobj,
    ebrk_p: cobj,
    wbrk_p: cobj,
    lbrk_p: cobj,
    gh_a_p: cobj,
    sh_a_p: cobj,
    sm_a_p: cobj,
    ri_a_p: cobj,
    sm_aw_p: cobj,
    ipbl_p: cobj,
    kpblh_p: cobj,
    wsed_CL_p: cobj,
):
    went = Ptr[float](went_p)
    wet_CL = Ptr[float](wet_CL_p)
    web_CL = Ptr[float](web_CL_p)
    jtbu_CL = Ptr[float](jtbu_CL_p)
    jbbu_CL = Ptr[float](jbbu_CL_p)
    evhc_CL = Ptr[float](evhc_CL_p)
    jt2slv_CL = Ptr[float](jt2slv_CL_p)
    n2ht_CL = Ptr[float](n2ht_CL_p)
    n2hb_CL = Ptr[float](n2hb_CL_p)
    lwp_CL = Ptr[float](lwp_CL_p)
    opt_depth_CL = Ptr[float](opt_depth_CL_p)
    radinvfrac_CL = Ptr[float](radinvfrac_CL_p)
    radf_CL = Ptr[float](radf_CL_p)
    wstar_CL = Ptr[float](wstar_CL_p)
    wstar3fact_CL = Ptr[float](wstar3fact_CL_p)
    ricl = Ptr[float](ricl_p)
    ghcl = Ptr[float](ghcl_p)
    shcl = Ptr[float](shcl_p)
    smcl = Ptr[float](smcl_p)
    ebrk = Ptr[float](ebrk_p)
    wbrk = Ptr[float](wbrk_p)
    lbrk = Ptr[float](lbrk_p)
    gh_a = Ptr[float](gh_a_p)
    sh_a = Ptr[float](sh_a_p)
    sm_a = Ptr[float](sm_a_p)
    ri_a = Ptr[float](ri_a_p)
    sm_aw = Ptr[float](sm_aw_p)
    ipbl = Ptr[i32](ipbl_p)
    kpblh = Ptr[i32](kpblh_p)
    wsed_CL = Ptr[float](wsed_CL_p)

    for i in range(1, ncol + 1):
        went[i - 1] = 0.0

        for ncv in range(1, ncvmax + 1):
            idx = _idx2(i, ncv, pcols)
            wet_CL[idx] = 0.0
            web_CL[idx] = 0.0
            jtbu_CL[idx] = 0.0
            jbbu_CL[idx] = 0.0
            evhc_CL[idx] = 0.0
            jt2slv_CL[idx] = 0.0
            n2ht_CL[idx] = 0.0
            n2hb_CL[idx] = 0.0
            lwp_CL[idx] = 0.0
            opt_depth_CL[idx] = 0.0
            radinvfrac_CL[idx] = 0.0
            radf_CL[idx] = 0.0
            wstar_CL[idx] = 0.0
            wstar3fact_CL[idx] = 0.0
            ricl[idx] = 0.0
            ghcl[idx] = 0.0
            shcl[idx] = 0.0
            smcl[idx] = 0.0
            ebrk[idx] = 0.0
            wbrk[idx] = 0.0
            lbrk[idx] = 0.0
            wsed_CL[idx] = 0.0

        for k in range(1, pver + 2):
            idx = _idx2(i, k, pcols)
            gh_a[idx] = 0.0
            sh_a[idx] = 0.0
            sm_a[idx] = 0.0
            ri_a[idx] = 0.0
            sm_aw[idx] = 0.0

        ipbl[i - 1] = i32(0)
        kpblh[i - 1] = i32(pver)


@export
def eddy_diff_caleddy_regime_diag_codon(
    i_col: int,
    pcols: int,
    ncvmax: int,
    kbase_p: cobj,
    ktop_p: cobj,
    ncvfin_p: cobj,
    kbase_diag_p: cobj,
    ktop_diag_p: cobj,
    ncvfin_diag_p: cobj,
):
    kbase = Ptr[i32](kbase_p)
    ktop = Ptr[i32](ktop_p)
    ncvfin = Ptr[i32](ncvfin_p)
    kbase_diag = Ptr[float](kbase_diag_p)
    ktop_diag = Ptr[float](ktop_diag_p)
    ncvfin_diag = Ptr[float](ncvfin_diag_p)

    i = i_col

    for ncv in range(1, ncvmax + 1):
        idx = _idx2(i, ncv, pcols)
        kbase_diag[idx] = float(kbase[idx])
        ktop_diag[idx] = float(ktop[idx])

    ncvfin_diag[i - 1] = float(ncvfin[i - 1])


@export
def eddy_diff_caleddy_stable_config_codon(
    ricrit: float,
    b1: float,
    alph2: float,
    alph3: float,
    alph4: float,
    alph5: float,
    alph4exs_p: cobj,
    ghmin_p: cobj,
    status_p: cobj,
):
    alph4exs = Ptr[float](alph4exs_p)
    ghmin = Ptr[float](ghmin_p)
    status = Ptr[i32](status_p)

    status[0] = i32(0)

    if ricrit == 0.19:
        alph4exs[0] = alph4
        ghmin[0] = -3.5334
    elif ricrit > 0.19:
        alph4exs[0] = -2.0 * b1 * alph2 / (alph3 - 2.0 * b1 * alph5) / ricrit
        ghmin[0] = -1.0e10
    else:
        status[0] = i32(1)


@export
def eddy_diff_caleddy_surface_tke_codon(
    i_col: int,
    pcols: int,
    pver: int,
    b1: float,
    vk: float,
    tkemax: float,
    z_p: cobj,
    bprod_p: cobj,
    sprod_p: cobj,
    tkes_p: cobj,
    tke_p: cobj,
    wcap_p: cobj,
):
    z = Ptr[float](z_p)
    bprod = Ptr[float](bprod_p)
    sprod = Ptr[float](sprod_p)
    tkes = Ptr[float](tkes_p)
    tke = Ptr[float](tke_p)
    wcap = Ptr[float](wcap_p)

    i = i_col
    surf_idx = _idx2(i, pver + 1, pcols)

    tkes[i - 1] = max(b1 * vk * z[_idx2(i, pver, pcols)] * (bprod[surf_idx] + sprod[surf_idx]), 1.0e-7) ** (2.0 / 3.0)
    tkes[i - 1] = min(tkes[i - 1], tkemax)
    tke[surf_idx] = tkes[i - 1]
    wcap[surf_idx] = tkes[i - 1] / b1


@inline
def _eddy_diff_zisocl_surface_energy_values(
    alph1: float,
    alph2: float,
    alph3: float,
    alph4: float,
    alph5: float,
    b1: float,
    vk: float,
    z_surf: float,
    bprod_surf: float,
    sprod_surf: float,
    tkes_surf: float,
):
    gg = 0.5 * vk * z_surf * bprod_surf / (tkes_surf ** (3.0 / 2.0))
    gh = gg / (alph5 - gg * alph3)
    gh = min(max(gh, -3.5334), 0.0233)
    sh = alph5 / (1.0 + alph3 * gh)
    sm = (alph1 + alph2 * gh) / (1.0 + alph3 * gh) / (1.0 + alph4 * gh)
    dlint_surf = z_surf
    dl2n2_surf = -vk * (z_surf ** 2.0) * bprod_surf / (sh * sqrt(tkes_surf))
    dl2s2_surf = vk * (z_surf ** 2.0) * sprod_surf / (sm * sqrt(tkes_surf))
    dw_surf = (tkes_surf / b1) * z_surf
    return gh, sh, sm, dlint_surf, dl2n2_surf, dl2s2_surf, dw_surf


@inline
def _eddy_diff_zisocl_stability_values(
    alph1: float,
    alph2: float,
    alph3: float,
    alph4: float,
    alph5: float,
    b1: float,
    ntzero: float,
    ricrit: float,
    l2n2: float,
    l2s2: float,
):
    ricll = min(l2n2 / max(l2s2, ntzero), ricrit)
    trma = alph3 * alph4 * ricll + 2.0 * b1 * (alph2 - alph4 * alph5 * ricll)
    trmb = ricll * (alph3 + alph4) + 2.0 * b1 * (-alph5 * ricll + alph1)
    trmc = ricll
    det = max(trmb * trmb - 4.0 * trma * trmc, 0.0)
    gh = (-trmb + sqrt(det)) / 2.0 / trma
    gh = min(max(gh, -3.5334), 0.0233)
    sh = alph5 / (1.0 + alph3 * gh)
    sm = (alph1 + alph2 * gh) / (1.0 + alph3 * gh) / (1.0 + alph4 * gh)
    return ricll, gh, sh, sm


@inline
def _eddy_diff_zisocl_layer_energy_values(
    tunl_mode: int,
    leng_mode: int,
    tunl: float,
    ctunl: float,
    cleng: float,
    lbulk: float,
    vk: float,
    i_col: int,
    k_ifc: int,
    pcols: int,
    z: Ptr[float],
    zi: Ptr[float],
    n2: Ptr[float],
    s2: Ptr[float],
    leng_max: Ptr[float],
):
    if tunl_mode == 1:
        tunlramp = 0.5 * (1.0 + ctunl) * tunl
    elif tunl_mode == 2:
        tunlramp = ctunl * tunl
    else:
        tunlramp = tunl

    if leng_mode == 0:
        lz = ((vk * zi[_idx2(i_col, k_ifc, pcols)]) ** (-cleng) + (tunlramp * lbulk) ** (-cleng)) ** (-1.0 / cleng)
    else:
        lz = min(vk * zi[_idx2(i_col, k_ifc, pcols)], tunlramp * lbulk)
    lz = min(leng_max[k_ifc - 1], lz)

    dzinc = z[_idx2(i_col, k_ifc - 1, pcols)] - z[_idx2(i_col, k_ifc, pcols)]
    dl2n2 = lz * lz * n2[_idx2(i_col, k_ifc, pcols)] * dzinc
    dl2s2 = lz * lz * s2[_idx2(i_col, k_ifc, pcols)] * dzinc
    return dzinc, dl2n2, dl2s2


@inline
def _eddy_diff_zisocl_interface_energy_values(
    tunl_mode: int,
    leng_mode: int,
    tunl: float,
    ctunl: float,
    cleng: float,
    lbulk: float,
    sh: float,
    sm: float,
    vk: float,
    i_col: int,
    k_ifc: int,
    pcols: int,
    z: Ptr[float],
    zi: Ptr[float],
    n2: Ptr[float],
    s2: Ptr[float],
    leng_max: Ptr[float],
):
    dzinc, dl2n2, dl2s2 = _eddy_diff_zisocl_layer_energy_values(
        tunl_mode, leng_mode, tunl, ctunl, cleng, lbulk, vk, i_col, k_ifc, pcols, z, zi, n2, s2, leng_max
    )
    dwinc = -sh * dl2n2 + sm * dl2s2
    return dzinc, dl2n2, dl2s2, dwinc


@export
def eddy_diff_zisocl_surface_energy_codon(
    alph1: float,
    alph2: float,
    alph3: float,
    alph4: float,
    alph5: float,
    b1: float,
    vk: float,
    z_surf: float,
    bprod_surf: float,
    sprod_surf: float,
    tkes_surf: float,
    gh_p: cobj,
    sh_p: cobj,
    sm_p: cobj,
    dlint_surf_p: cobj,
    dl2n2_surf_p: cobj,
    dl2s2_surf_p: cobj,
    dw_surf_p: cobj,
):
    gh = Ptr[float](gh_p)
    sh = Ptr[float](sh_p)
    sm = Ptr[float](sm_p)
    dlint_surf = Ptr[float](dlint_surf_p)
    dl2n2_surf = Ptr[float](dl2n2_surf_p)
    dl2s2_surf = Ptr[float](dl2s2_surf_p)
    dw_surf = Ptr[float](dw_surf_p)

    gh_val, sh_val, sm_val, dlint_surf_val, dl2n2_surf_val, dl2s2_surf_val, dw_surf_val = _eddy_diff_zisocl_surface_energy_values(
        alph1, alph2, alph3, alph4, alph5, b1, vk, z_surf, bprod_surf, sprod_surf, tkes_surf
    )
    gh[0] = gh_val
    sh[0] = sh_val
    sm[0] = sm_val
    dlint_surf[0] = dlint_surf_val
    dl2n2_surf[0] = dl2n2_surf_val
    dl2s2_surf[0] = dl2s2_surf_val
    dw_surf[0] = dw_surf_val


@export
def eddy_diff_zisocl_surface_state_codon(
    alph1: float,
    alph2: float,
    alph3: float,
    alph4: float,
    alph5: float,
    b1: float,
    vk: float,
    lbulk_max: float,
    kb_is_surface: int,
    use_dw_surf: int,
    zi_top: float,
    zi_base: float,
    z_surf: float,
    bflxs_surf: float,
    bprod_surf: float,
    sprod_surf: float,
    tkes_surf: float,
    lbulk_p: cobj,
    gh_p: cobj,
    sh_p: cobj,
    sm_p: cobj,
    dlint_surf_p: cobj,
    dl2n2_surf_p: cobj,
    dl2s2_surf_p: cobj,
    dw_surf_p: cobj,
    lint_p: cobj,
    l2n2_p: cobj,
    l2s2_p: cobj,
    wint_p: cobj,
):
    lbulk = Ptr[float](lbulk_p)
    gh = Ptr[float](gh_p)
    sh = Ptr[float](sh_p)
    sm = Ptr[float](sm_p)
    dlint_surf = Ptr[float](dlint_surf_p)
    dl2n2_surf = Ptr[float](dl2n2_surf_p)
    dl2s2_surf = Ptr[float](dl2s2_surf_p)
    dw_surf = Ptr[float](dw_surf_p)
    lint = Ptr[float](lint_p)
    l2n2 = Ptr[float](l2n2_p)
    l2s2 = Ptr[float](l2s2_p)
    wint = Ptr[float](wint_p)

    lbulk[0] = zi_top - zi_base
    lbulk[0] = min(lbulk[0], lbulk_max)
    dlint_surf[0] = 0.0
    dl2n2_surf[0] = 0.0
    dl2s2_surf[0] = 0.0
    dw_surf[0] = 0.0

    if kb_is_surface != 0:
        if bflxs_surf > 0.0:
            gh_val, sh_val, sm_val, dlint_surf_val, dl2n2_surf_val, dl2s2_surf_val, dw_surf_val = _eddy_diff_zisocl_surface_energy_values(
                alph1, alph2, alph3, alph4, alph5, b1, vk, z_surf, bprod_surf, sprod_surf, tkes_surf
            )
            gh[0] = gh_val
            sh[0] = sh_val
            sm[0] = sm_val
            dlint_surf[0] = dlint_surf_val
            dl2n2_surf[0] = dl2n2_surf_val
            dl2s2_surf[0] = dl2s2_surf_val
            dw_surf[0] = dw_surf_val
        else:
            lbulk[0] = zi_top - z_surf
            lbulk[0] = min(lbulk[0], lbulk_max)

    lint[0] = dlint_surf[0]
    l2n2[0] = dl2n2_surf[0]
    l2s2[0] = dl2s2_surf[0]
    wint[0] = dw_surf[0]

    if use_dw_surf != 0:
        l2n2[0] = 0.0
        l2s2[0] = 0.0
    else:
        wint[0] = 0.0


@export
def eddy_diff_zisocl_surface_extend_codon(
    alph1: float,
    alph2: float,
    alph3: float,
    alph4: float,
    alph5: float,
    b1: float,
    vk: float,
    tkemax: float,
    bflxs_surf: float,
    z_surf: float,
    bprod_surf: float,
    sprod_surf: float,
    tkes_surf: float,
    sh: float,
    dlint_surf_p: cobj,
    dl2n2_surf_p: cobj,
    dl2s2_surf_p: cobj,
    dw_surf_p: cobj,
    lint_p: cobj,
    l2n2_p: cobj,
    l2s2_p: cobj,
    wint_p: cobj,
):
    dlint_surf = Ptr[float](dlint_surf_p)
    dl2n2_surf = Ptr[float](dl2n2_surf_p)
    dl2s2_surf = Ptr[float](dl2s2_surf_p)
    dw_surf = Ptr[float](dw_surf_p)
    lint = Ptr[float](lint_p)
    l2n2 = Ptr[float](l2n2_p)
    l2s2 = Ptr[float](l2s2_p)
    wint = Ptr[float](wint_p)

    if bflxs_surf > 0.0:
        gh_surf, sh_surf, sm_surf, dlint_surf_val, dl2n2_surf_val, dl2s2_surf_val, dw_surf_val = _eddy_diff_zisocl_surface_energy_values(
            alph1, alph2, alph3, alph4, alph5, b1, vk, z_surf, bprod_surf, sprod_surf, tkes_surf
        )
        dlint_surf[0] = dlint_surf_val
        dl2n2_surf[0] = dl2n2_surf_val
        dl2s2_surf[0] = dl2s2_surf_val
        dw_surf[0] = dw_surf_val
    else:
        dlint_surf[0] = 0.0
        dl2n2_surf[0] = 0.0
        dl2s2_surf[0] = 0.0
        dw_surf[0] = 0.0

    lint[0] = lint[0] + dlint_surf[0]
    l2n2[0] = l2n2[0] + dl2n2_surf[0]
    l2n2[0] = -min(-l2n2[0], tkemax * lint[0] / (b1 * sh))
    l2s2[0] = l2s2[0] + dl2s2_surf[0]
    wint[0] = wint[0] + dw_surf[0]


@export
def eddy_diff_zisocl_sbcl_state_codon(
    choice_tkes_ebprod: int,
    sh: float,
    dlint_surf: float,
    dl2n2_surf: float,
    dl2s2_surf: float,
    dw_surf: float,
    lint_p: cobj,
    l2n2_p: cobj,
    l2s2_p: cobj,
    wint_p: cobj,
):
    lint = Ptr[float](lint_p)
    l2n2 = Ptr[float](l2n2_p)
    l2s2 = Ptr[float](l2s2_p)
    wint = Ptr[float](wint_p)

    lint[0] = dlint_surf
    l2n2[0] = dl2n2_surf
    l2s2[0] = dl2s2_surf
    wint[0] = dw_surf

    if choice_tkes_ebprod != 0:
        l2n2[0] = -wint[0] / sh


@export
def eddy_diff_zisocl_initial_state_codon(
    i_col: int,
    kt_ifc: int,
    kb_ifc: int,
    pcols: int,
    pver: int,
    kb_is_surface: int,
    use_dw_surf: int,
    choice_tkes_ebprod: int,
    tunl_mode: int,
    leng_mode: int,
    alph1: float,
    alph2: float,
    alph3: float,
    alph4: float,
    alph5: float,
    b1: float,
    vk: float,
    ntzero: float,
    ricrit: float,
    lbulk_max: float,
    tunl: float,
    ctunl: float,
    cleng: float,
    tkemax: float,
    z_surf: float,
    bflxs_surf: float,
    bprod_surf: float,
    sprod_surf: float,
    tkes_surf: float,
    z_p: cobj,
    zi_p: cobj,
    n2_p: cobj,
    s2_p: cobj,
    leng_max_p: cobj,
    lbulk_p: cobj,
    gh_p: cobj,
    sh_p: cobj,
    sm_p: cobj,
    dlint_surf_p: cobj,
    dl2n2_surf_p: cobj,
    dl2s2_surf_p: cobj,
    dw_surf_p: cobj,
    lint_p: cobj,
    l2n2_p: cobj,
    l2s2_p: cobj,
    wint_p: cobj,
    ricll_p: cobj,
):
    z = Ptr[float](z_p)
    zi = Ptr[float](zi_p)
    n2 = Ptr[float](n2_p)
    s2 = Ptr[float](s2_p)
    leng_max = Ptr[float](leng_max_p)
    lbulk = Ptr[float](lbulk_p)
    gh = Ptr[float](gh_p)
    sh = Ptr[float](sh_p)
    sm = Ptr[float](sm_p)
    dlint_surf = Ptr[float](dlint_surf_p)
    dl2n2_surf = Ptr[float](dl2n2_surf_p)
    dl2s2_surf = Ptr[float](dl2s2_surf_p)
    dw_surf = Ptr[float](dw_surf_p)
    lint = Ptr[float](lint_p)
    l2n2 = Ptr[float](l2n2_p)
    l2s2 = Ptr[float](l2s2_p)
    wint = Ptr[float](wint_p)
    ricll = Ptr[float](ricll_p)

    lbulk[0] = zi[_idx2(i_col, kt_ifc, pcols)] - zi[_idx2(i_col, kb_ifc, pcols)]
    lbulk[0] = min(lbulk[0], lbulk_max)
    dlint_surf[0] = 0.0
    dl2n2_surf[0] = 0.0
    dl2s2_surf[0] = 0.0
    dw_surf[0] = 0.0

    if kb_is_surface != 0:
        if bflxs_surf > 0.0:
            gh_val, sh_val, sm_val, dlint_surf_val, dl2n2_surf_val, dl2s2_surf_val, dw_surf_val = _eddy_diff_zisocl_surface_energy_values(
                alph1, alph2, alph3, alph4, alph5, b1, vk, z_surf, bprod_surf, sprod_surf, tkes_surf
            )
            gh[0] = gh_val
            sh[0] = sh_val
            sm[0] = sm_val
            dlint_surf[0] = dlint_surf_val
            dl2n2_surf[0] = dl2n2_surf_val
            dl2s2_surf[0] = dl2s2_surf_val
            dw_surf[0] = dw_surf_val
        else:
            lbulk[0] = zi[_idx2(i_col, kt_ifc, pcols)] - z_surf
            lbulk[0] = min(lbulk[0], lbulk_max)

    lint[0] = dlint_surf[0]
    l2n2[0] = dl2n2_surf[0]
    l2s2[0] = dl2s2_surf[0]
    wint[0] = dw_surf[0]

    if use_dw_surf != 0:
        l2n2[0] = 0.0
        l2s2[0] = 0.0
    else:
        wint[0] = 0.0

    if kb_ifc == pver + 1 and bflxs_surf > 0.0:
        ricll[0] = min(-(sm[0] / sh[0]) * (bprod_surf / sprod_surf), ricrit)

    if kt_ifc < kb_ifc - 1:
        k_ifc = kb_ifc - 1
        while k_ifc >= kt_ifc + 1:
            dzinc, dl2n2, dl2s2 = _eddy_diff_zisocl_layer_energy_values(
                tunl_mode, leng_mode, tunl, ctunl, cleng, lbulk[0], vk, i_col, k_ifc, pcols, z, zi, n2, s2, leng_max
            )

            l2n2[0] = l2n2[0] + dl2n2
            l2s2[0] = l2s2[0] + dl2s2
            lint[0] = lint[0] + dzinc

            k_ifc -= 1

        ricll_val, gh_val, sh_val, sm_val = _eddy_diff_zisocl_stability_values(
            alph1, alph2, alph3, alph4, alph5, b1, ntzero, ricrit, l2n2[0], l2s2[0]
        )
        ricll[0] = ricll_val
        gh[0] = gh_val
        sh[0] = sh_val
        sm[0] = sm_val
        wint[0] = wint[0] - sh[0] * l2n2[0] + sm[0] * l2s2[0]
    else:
        lint[0] = dlint_surf[0]
        l2n2[0] = dl2n2_surf[0]
        l2s2[0] = dl2s2_surf[0]
        wint[0] = dw_surf[0]
        if choice_tkes_ebprod != 0:
            l2n2[0] = -wint[0] / sh[0]

    l2n2[0] = -min(-l2n2[0], tkemax * lint[0] / (b1 * sh[0]))
    l2s2[0] = min(l2s2[0], tkemax * lint[0] / (b1 * sm[0]))


@export
def eddy_diff_zisocl_extended_state_codon(
    i_col: int,
    kt_ifc: int,
    kb_ifc: int,
    pcols: int,
    pver: int,
    use_dw_surf: int,
    tunl_mode: int,
    leng_mode: int,
    alph1: float,
    alph2: float,
    alph3: float,
    alph4: float,
    alph5: float,
    b1: float,
    vk: float,
    ntzero: float,
    ricrit: float,
    lbulk_max: float,
    tunl: float,
    ctunl: float,
    cleng: float,
    tkemax: float,
    zi_top: float,
    zi_base: float,
    z_surf: float,
    bflxs_surf: float,
    bprod_surf: float,
    sprod_surf: float,
    tkes_surf: float,
    z_p: cobj,
    zi_p: cobj,
    n2_p: cobj,
    s2_p: cobj,
    leng_max_p: cobj,
    gh_p: cobj,
    sh_p: cobj,
    sm_p: cobj,
    lint_p: cobj,
    wint_p: cobj,
    ricll_p: cobj,
):
    z = Ptr[float](z_p)
    zi = Ptr[float](zi_p)
    n2 = Ptr[float](n2_p)
    s2 = Ptr[float](s2_p)
    leng_max = Ptr[float](leng_max_p)
    gh = Ptr[float](gh_p)
    sh = Ptr[float](sh_p)
    sm = Ptr[float](sm_p)
    lint = Ptr[float](lint_p)
    wint = Ptr[float](wint_p)
    ricll = Ptr[float](ricll_p)

    lbulk = zi_top - zi_base
    lbulk = min(lbulk, lbulk_max)
    dlint_surf = 0.0
    dl2n2_surf = 0.0
    dl2s2_surf = 0.0
    dw_surf = 0.0

    if kb_ifc == pver + 1:
        if bflxs_surf > 0.0:
            gh_surf, sh_surf, sm_surf, dlint_surf, dl2n2_surf, dl2s2_surf, dw_surf = _eddy_diff_zisocl_surface_energy_values(
                alph1, alph2, alph3, alph4, alph5, b1, vk, z_surf, bprod_surf, sprod_surf, tkes_surf
            )
        else:
            lbulk = zi_top - z_surf
            lbulk = min(lbulk, lbulk_max)

    lint[0] = dlint_surf
    l2n2 = dl2n2_surf
    l2s2 = dl2s2_surf
    wint[0] = dw_surf

    if use_dw_surf != 0:
        l2n2 = 0.0
        l2s2 = 0.0
    else:
        wint[0] = 0.0

    for k_ifc in range(kt_ifc + 1, kb_ifc):
        dzinc, dl2n2, dl2s2 = _eddy_diff_zisocl_layer_energy_values(
            tunl_mode, leng_mode, tunl, ctunl, cleng, lbulk, vk, i_col, k_ifc, pcols, z, zi, n2, s2, leng_max
        )

        lint[0] = lint[0] + dzinc
        l2n2 = l2n2 + dl2n2
        l2s2 = l2s2 + dl2s2

    ricll_val, gh_val, sh_val, sm_val = _eddy_diff_zisocl_stability_values(
        alph1, alph2, alph3, alph4, alph5, b1, ntzero, ricrit, l2n2, l2s2
    )
    ricll[0] = ricll_val
    gh[0] = gh_val
    sh[0] = sh_val
    sm[0] = sm_val
    wint[0] = max(wint[0] - sh[0] * l2n2 + sm[0] * l2s2, 0.01)


@export
def eddy_diff_zisocl_non_sbcl_state_codon(
    i_col: int,
    kt_ifc: int,
    kb_ifc: int,
    pcols: int,
    pver: int,
    tunl_mode: int,
    leng_mode: int,
    alph1: float,
    alph2: float,
    alph3: float,
    alph4: float,
    alph5: float,
    b1: float,
    vk: float,
    ntzero: float,
    ricrit: float,
    tunl: float,
    ctunl: float,
    cleng: float,
    lbulk: float,
    z_p: cobj,
    zi_p: cobj,
    n2_p: cobj,
    s2_p: cobj,
    leng_max_p: cobj,
    lint_p: cobj,
    l2n2_p: cobj,
    l2s2_p: cobj,
    wint_p: cobj,
    ricll_p: cobj,
    gh_p: cobj,
    sh_p: cobj,
    sm_p: cobj,
):
    z = Ptr[float](z_p)
    zi = Ptr[float](zi_p)
    n2 = Ptr[float](n2_p)
    s2 = Ptr[float](s2_p)
    leng_max = Ptr[float](leng_max_p)
    lint = Ptr[float](lint_p)
    l2n2 = Ptr[float](l2n2_p)
    l2s2 = Ptr[float](l2s2_p)
    wint = Ptr[float](wint_p)
    ricll = Ptr[float](ricll_p)
    gh = Ptr[float](gh_p)
    sh = Ptr[float](sh_p)
    sm = Ptr[float](sm_p)

    k_ifc = kb_ifc - 1
    while k_ifc >= kt_ifc + 1:
        dzinc, dl2n2, dl2s2 = _eddy_diff_zisocl_layer_energy_values(
            tunl_mode, leng_mode, tunl, ctunl, cleng, lbulk, vk, i_col, k_ifc, pcols, z, zi, n2, s2, leng_max
        )

        l2n2[0] = l2n2[0] + dl2n2
        l2s2[0] = l2s2[0] + dl2s2
        lint[0] = lint[0] + dzinc

        k_ifc -= 1

    ricll_val, gh_val, sh_val, sm_val = _eddy_diff_zisocl_stability_values(
        alph1, alph2, alph3, alph4, alph5, b1, ntzero, ricrit, l2n2[0], l2s2[0]
    )
    ricll[0] = ricll_val
    gh[0] = gh_val
    sh[0] = sh_val
    sm[0] = sm_val
    wint[0] = wint[0] - sh[0] * l2n2[0] + sm[0] * l2s2[0]


@export
def eddy_diff_zisocl_upward_state_codon(
    i_col: int,
    pcols: int,
    pver: int,
    ntop_turb: int,
    ncv_col: int,
    tunl_mode: int,
    leng_mode: int,
    tunl: float,
    ctunl: float,
    cleng: float,
    vk: float,
    rinc: float,
    tkemax: float,
    b1: float,
    lbulk: float,
    sh: float,
    sm: float,
    z_p: cobj,
    zi_p: cobj,
    n2_p: cobj,
    s2_p: cobj,
    leng_max_p: cobj,
    ncvfin_p: cobj,
    kbase_p: cobj,
    ktop_p: cobj,
    kt_p: cobj,
    lint_p: cobj,
    l2n2_p: cobj,
    l2s2_p: cobj,
    wint_p: cobj,
    status_p: cobj,
):
    z = Ptr[float](z_p)
    zi = Ptr[float](zi_p)
    n2 = Ptr[float](n2_p)
    s2 = Ptr[float](s2_p)
    leng_max = Ptr[float](leng_max_p)
    ncvfin = Ptr[i32](ncvfin_p)
    kbase = Ptr[i32](kbase_p)
    ktop = Ptr[i32](ktop_p)
    kt_state = Ptr[int](kt_p)
    lint = Ptr[float](lint_p)
    l2n2 = Ptr[float](l2n2_p)
    l2s2 = Ptr[float](l2s2_p)
    wint = Ptr[float](wint_p)
    status = Ptr[int](status_p)

    kt = int(kt_state[0])
    cntu = 0
    status[0] = 0

    dzinc, dl2n2, dl2s2, dwinc = _eddy_diff_zisocl_interface_energy_values(
        tunl_mode, leng_mode, tunl, ctunl, cleng, lbulk, sh, sm, vk, i_col, kt, pcols, z, zi, n2, s2, leng_max
    )

    while (-dl2n2 > (-rinc * l2n2[0] / (1.0 - rinc))) and (
        kt > ntop_turb + 2 or z[_idx2(i_col, kt, pcols)] < 50000.0
    ):
        lint[0] = lint[0] + dzinc
        l2n2[0] = l2n2[0] + dl2n2
        l2n2[0] = -min(-l2n2[0], tkemax * lint[0] / (b1 * sh))
        l2s2[0] = l2s2[0] + dl2s2
        wint[0] = wint[0] + dwinc

        kt = kt - 1
        if kt == ntop_turb:
            kt_state[0] = kt
            status[0] = 1
            return

        ktinc = int(kbase[_idx2(i_col, ncv_col + cntu + 1, pcols)]) - 1

        if kt == ktinc:
            k_ifc = int(kbase[_idx2(i_col, ncv_col + cntu + 1, pcols)]) - 1
            while k_ifc >= int(ktop[_idx2(i_col, ncv_col + cntu + 1, pcols)]) + 1:
                dzinc, dl2n2, dl2s2, dwinc = _eddy_diff_zisocl_interface_energy_values(
                    tunl_mode, leng_mode, tunl, ctunl, cleng, lbulk, sh, sm, vk, i_col, k_ifc, pcols, z, zi, n2, s2, leng_max
                )

                lint[0] = lint[0] + dzinc
                l2n2[0] = l2n2[0] + dl2n2
                l2n2[0] = -min(-l2n2[0], tkemax * lint[0] / (b1 * sh))
                l2s2[0] = l2s2[0] + dl2s2
                wint[0] = wint[0] + dwinc

                k_ifc -= 1

            kt = int(ktop[_idx2(i_col, ncv_col + cntu + 1, pcols)])
            ncvfin[i_col - 1] = i32(int(ncvfin[i_col - 1]) - 1)
            cntu += 1

        dzinc, dl2n2, dl2s2, dwinc = _eddy_diff_zisocl_interface_energy_values(
            tunl_mode, leng_mode, tunl, ctunl, cleng, lbulk, sh, sm, vk, i_col, kt, pcols, z, zi, n2, s2, leng_max
        )

    if cntu > 0:
        incv = 1
        while incv <= int(ncvfin[i_col - 1]) - ncv_col:
            kbase[_idx2(i_col, ncv_col + incv, pcols)] = kbase[_idx2(i_col, ncv_col + cntu + incv, pcols)]
            ktop[_idx2(i_col, ncv_col + incv, pcols)] = ktop[_idx2(i_col, ncv_col + cntu + incv, pcols)]
            incv += 1

    kt_state[0] = kt


@export
def eddy_diff_zisocl_downward_state_codon(
    i_col: int,
    pcols: int,
    pver: int,
    ncvinit: int,
    tunl_mode: int,
    leng_mode: int,
    tunl: float,
    ctunl: float,
    cleng: float,
    vk: float,
    rinc: float,
    tkemax: float,
    b1: float,
    alph1: float,
    alph2: float,
    alph3: float,
    alph4: float,
    alph5: float,
    lbulk: float,
    sh: float,
    sm: float,
    z_surf: float,
    bflxs_surf: float,
    bprod_surf: float,
    sprod_surf: float,
    tkes_surf: float,
    z_p: cobj,
    zi_p: cobj,
    n2_p: cobj,
    s2_p: cobj,
    leng_max_p: cobj,
    dlint_surf_p: cobj,
    dl2n2_surf_p: cobj,
    dl2s2_surf_p: cobj,
    dw_surf_p: cobj,
    ncvfin_p: cobj,
    kbase_p: cobj,
    ktop_p: cobj,
    kb_p: cobj,
    ncv_p: cobj,
    lint_p: cobj,
    l2n2_p: cobj,
    l2s2_p: cobj,
    wint_p: cobj,
    status_p: cobj,
):
    z = Ptr[float](z_p)
    zi = Ptr[float](zi_p)
    n2 = Ptr[float](n2_p)
    s2 = Ptr[float](s2_p)
    leng_max = Ptr[float](leng_max_p)
    dlint_surf = Ptr[float](dlint_surf_p)
    dl2n2_surf = Ptr[float](dl2n2_surf_p)
    dl2s2_surf = Ptr[float](dl2s2_surf_p)
    dw_surf = Ptr[float](dw_surf_p)
    ncvfin = Ptr[i32](ncvfin_p)
    kbase = Ptr[i32](kbase_p)
    ktop = Ptr[i32](ktop_p)
    kb_state = Ptr[int](kb_p)
    ncv_state = Ptr[int](ncv_p)
    lint = Ptr[float](lint_p)
    l2n2 = Ptr[float](l2n2_p)
    l2s2 = Ptr[float](l2s2_p)
    wint = Ptr[float](wint_p)
    status = Ptr[int](status_p)

    kb = int(kb_state[0])
    ncv = int(ncv_state[0])
    cntd = 0
    status[0] = 0

    if kb == pver + 1:
        return

    dzinc, dl2n2, dl2s2, dwinc = _eddy_diff_zisocl_interface_energy_values(
        tunl_mode, leng_mode, tunl, ctunl, cleng, lbulk, sh, sm, vk, i_col, kb, pcols, z, zi, n2, s2, leng_max
    )

    while (-dl2n2 > (-rinc * l2n2[0] / (1.0 - rinc))) and (kb != pver + 1):
        lint[0] = lint[0] + dzinc
        l2n2[0] = l2n2[0] + dl2n2
        l2n2[0] = -min(-l2n2[0], tkemax * lint[0] / (b1 * sh))
        l2s2[0] = l2s2[0] + dl2s2
        wint[0] = wint[0] + dwinc

        kb = kb + 1

        kbinc = 0
        if ncv > 1:
            kbinc = int(ktop[_idx2(i_col, ncv - 1, pcols)]) + 1
        if kb == kbinc:
            k_ifc = int(ktop[_idx2(i_col, ncv - 1, pcols)]) + 1
            while k_ifc <= int(kbase[_idx2(i_col, ncv - 1, pcols)]) - 1:
                dzinc, dl2n2, dl2s2, dwinc = _eddy_diff_zisocl_interface_energy_values(
                    tunl_mode, leng_mode, tunl, ctunl, cleng, lbulk, sh, sm, vk, i_col, k_ifc, pcols, z, zi, n2, s2, leng_max
                )

                lint[0] = lint[0] + dzinc
                l2n2[0] = l2n2[0] + dl2n2
                l2n2[0] = -min(-l2n2[0], tkemax * lint[0] / (b1 * sh))
                l2s2[0] = l2s2[0] + dl2s2
                wint[0] = wint[0] + dwinc

                k_ifc += 1

            kb = int(kbase[_idx2(i_col, ncv - 1, pcols)])
            ncv = ncv - 1
            ncvfin[i_col - 1] = i32(int(ncvfin[i_col - 1]) - 1)
            cntd += 1

        if kb == pver + 1:
            if bflxs_surf > 0.0:
                gh_surf, sh_surf, sm_surf, dlint_surf_val, dl2n2_surf_val, dl2s2_surf_val, dw_surf_val = _eddy_diff_zisocl_surface_energy_values(
                    alph1, alph2, alph3, alph4, alph5, b1, vk, z_surf, bprod_surf, sprod_surf, tkes_surf
                )
                dlint_surf[0] = dlint_surf_val
                dl2n2_surf[0] = dl2n2_surf_val
                dl2s2_surf[0] = dl2s2_surf_val
                dw_surf[0] = dw_surf_val
            else:
                dlint_surf[0] = 0.0
                dl2n2_surf[0] = 0.0
                dl2s2_surf[0] = 0.0
                dw_surf[0] = 0.0

            lint[0] = lint[0] + dlint_surf[0]
            l2n2[0] = l2n2[0] + dl2n2_surf[0]
            l2n2[0] = -min(-l2n2[0], tkemax * lint[0] / (b1 * sh))
            l2s2[0] = l2s2[0] + dl2s2_surf[0]
            wint[0] = wint[0] + dw_surf[0]
        else:
            dzinc, dl2n2, dl2s2, dwinc = _eddy_diff_zisocl_interface_energy_values(
                tunl_mode, leng_mode, tunl, ctunl, cleng, lbulk, sh, sm, vk, i_col, kb, pcols, z, zi, n2, s2, leng_max
            )

    if kb == pver + 1 and ncv != 1:
        kb_state[0] = kb
        ncv_state[0] = ncv
        status[0] = 1
        return

    if cntd > 0:
        incv = 1
        while incv <= int(ncvfin[i_col - 1]) - ncv:
            kbase[_idx2(i_col, ncv + incv, pcols)] = kbase[_idx2(i_col, ncvinit + incv, pcols)]
            ktop[_idx2(i_col, ncv + incv, pcols)] = ktop[_idx2(i_col, ncvinit + incv, pcols)]
            incv += 1

    kb_state[0] = kb
    ncv_state[0] = ncv


@export
def eddy_diff_zisocl_stability_codon(
    alph1: float,
    alph2: float,
    alph3: float,
    alph4: float,
    alph5: float,
    b1: float,
    ntzero: float,
    ricrit: float,
    l2n2: float,
    l2s2: float,
    ricll_p: cobj,
    gh_p: cobj,
    sh_p: cobj,
    sm_p: cobj,
):
    ricll = Ptr[float](ricll_p)
    gh = Ptr[float](gh_p)
    sh = Ptr[float](sh_p)
    sm = Ptr[float](sm_p)

    ricll_val, gh_val, sh_val, sm_val = _eddy_diff_zisocl_stability_values(
        alph1, alph2, alph3, alph4, alph5, b1, ntzero, ricrit, l2n2, l2s2
    )
    ricll[0] = ricll_val
    gh[0] = gh_val
    sh[0] = sh_val
    sm[0] = sm_val


@export
def eddy_diff_zisocl_layer_energy_codon(
    i_col: int,
    k_ifc: int,
    pcols: int,
    pver: int,
    tunl_mode: int,
    leng_mode: int,
    tunl: float,
    ctunl: float,
    cleng: float,
    lbulk: float,
    vk: float,
    z_p: cobj,
    zi_p: cobj,
    n2_p: cobj,
    s2_p: cobj,
    leng_max_p: cobj,
    dzinc_p: cobj,
    dl2n2_p: cobj,
    dl2s2_p: cobj,
):
    z = Ptr[float](z_p)
    zi = Ptr[float](zi_p)
    n2 = Ptr[float](n2_p)
    s2 = Ptr[float](s2_p)
    leng_max = Ptr[float](leng_max_p)
    dzinc = Ptr[float](dzinc_p)
    dl2n2 = Ptr[float](dl2n2_p)
    dl2s2 = Ptr[float](dl2s2_p)

    dzinc_val, dl2n2_val, dl2s2_val = _eddy_diff_zisocl_layer_energy_values(
        tunl_mode, leng_mode, tunl, ctunl, cleng, lbulk, vk, i_col, k_ifc, pcols, z, zi, n2, s2, leng_max
    )
    dzinc[0] = dzinc_val
    dl2n2[0] = dl2n2_val
    dl2s2[0] = dl2s2_val


@export
def eddy_diff_zisocl_interface_energy_codon(
    i_col: int,
    k_ifc: int,
    pcols: int,
    pver: int,
    tunl_mode: int,
    leng_mode: int,
    tunl: float,
    ctunl: float,
    cleng: float,
    lbulk: float,
    sh: float,
    sm: float,
    vk: float,
    z_p: cobj,
    zi_p: cobj,
    n2_p: cobj,
    s2_p: cobj,
    leng_max_p: cobj,
    dzinc_p: cobj,
    dl2n2_p: cobj,
    dl2s2_p: cobj,
    dwinc_p: cobj,
):
    z = Ptr[float](z_p)
    zi = Ptr[float](zi_p)
    n2 = Ptr[float](n2_p)
    s2 = Ptr[float](s2_p)
    leng_max = Ptr[float](leng_max_p)
    dzinc = Ptr[float](dzinc_p)
    dl2n2 = Ptr[float](dl2n2_p)
    dl2s2 = Ptr[float](dl2s2_p)
    dwinc = Ptr[float](dwinc_p)

    dzinc_val, dl2n2_val, dl2s2_val, dwinc_val = _eddy_diff_zisocl_interface_energy_values(
        tunl_mode, leng_mode, tunl, ctunl, cleng, lbulk, sh, sm, vk, i_col, k_ifc, pcols, z, zi, n2, s2, leng_max
    )
    dzinc[0] = dzinc_val
    dl2n2[0] = dl2n2_val
    dl2s2[0] = dl2s2_val
    dwinc[0] = dwinc_val


@export
def eddy_diff_exacol_codon(
    ncol: int,
    pcols: int,
    pver: int,
    ncvmax: int,
    ntop_turb: int,
    ri_p: cobj,
    bflxs_p: cobj,
    ktop_p: cobj,
    kbase_p: cobj,
    ncvfin_p: cobj,
):
    ri = Ptr[float](ri_p)
    bflxs = Ptr[float](bflxs_p)
    ktop = Ptr[i32](ktop_p)
    kbase = Ptr[i32](kbase_p)
    ncvfin = Ptr[i32](ncvfin_p)

    rimaxentr = 0.0

    for i in range(1, ncol + 1):
        ncvfin[i - 1] = i32(0)
        for ncv in range(1, ncvmax + 1):
            ktop[_idx2(i, ncv, pcols)] = i32(0)
            kbase[_idx2(i, ncv, pcols)] = i32(0)

    for i in range(1, ncol + 1):
        ncv = 0
        k = pver + 1

        while k > ntop_turb + 1:
            if k == pver + 1:
                riex_k = rimaxentr - bflxs[i - 1]
            else:
                riex_k = ri[_idx2(i, k, pcols)]

            if riex_k < rimaxentr:
                ncv += 1
                kbase[_idx2(i, ncv, pcols)] = i32(min(k + 1, pver + 1))

                while k > ntop_turb + 1:
                    if k == pver + 1:
                        riex_k = rimaxentr - bflxs[i - 1]
                    else:
                        riex_k = ri[_idx2(i, k, pcols)]
                    if not (riex_k < rimaxentr):
                        break
                    k -= 1

                ktop[_idx2(i, ncv, pcols)] = i32(k)
            else:
                k -= 1

        ncvfin[i - 1] = i32(ncv)


@export
def eddy_diff_compute_radf_codon(
    i_col: int,
    pcols: int,
    pver: int,
    ncvmax: int,
    radf_mode: int,
    qmin: float,
    gravit: float,
    ncvfin_p: cobj,
    ktop_p: cobj,
    ql_p: cobj,
    pi_p: cobj,
    qrlw_p: cobj,
    cldeff_p: cobj,
    zi_p: cobj,
    chs_p: cobj,
    lwp_CL_p: cobj,
    opt_depth_CL_p: cobj,
    radinvfrac_CL_p: cobj,
    radf_CL_p: cobj,
):
    ncvfin = Ptr[i32](ncvfin_p)
    ktop = Ptr[i32](ktop_p)
    ql = Ptr[float](ql_p)
    pi = Ptr[float](pi_p)
    qrlw = Ptr[float](qrlw_p)
    cldeff = Ptr[float](cldeff_p)
    zi = Ptr[float](zi_p)
    chs = Ptr[float](chs_p)
    lwp_CL = Ptr[float](lwp_CL_p)
    opt_depth_CL = Ptr[float](opt_depth_CL_p)
    radinvfrac_CL = Ptr[float](radinvfrac_CL_p)
    radf_CL = Ptr[float](radf_CL_p)

    i = i_col
    ncvfin_i = int(ncvfin[i - 1])

    for ncv in range(1, ncvfin_i + 1):
        kt = int(ktop[_idx2(i, ncv, pcols)])

        lwp = 0.0
        opt_depth = 0.0
        radinvfrac = 0.0
        radf = 0.0

        if radf_mode == 0:
            if ql[_idx2(i, kt, pcols)] > qmin and ql[_idx2(i, kt - 1, pcols)] < qmin:
                lwp = ql[_idx2(i, kt, pcols)] * (pi[_idx2(i, kt + 1, pcols)] - pi[_idx2(i, kt, pcols)]) / gravit
                opt_depth = 156.0 * lwp
                radinvfrac = opt_depth * (4.0 + opt_depth) / (6.0 * (4.0 + opt_depth) + opt_depth**2)
                radf = qrlw[_idx2(i, kt, pcols)] / (pi[_idx2(i, kt, pcols)] - pi[_idx2(i, kt + 1, pcols)])
                radf = max(
                    radinvfrac * radf * (zi[_idx2(i, kt, pcols)] - zi[_idx2(i, kt + 1, pcols)]),
                    0.0,
                ) * chs[_idx2(i, kt, pcols)]

        elif radf_mode == 1:
            lwp = ql[_idx2(i, kt, pcols)] * (pi[_idx2(i, kt + 1, pcols)] - pi[_idx2(i, kt, pcols)]) / gravit
            opt_depth = 156.0 * lwp
            radinvfrac = opt_depth * (4.0 + opt_depth) / (6.0 * (4.0 + opt_depth) + opt_depth**2)
            radinvfrac = max(cldeff[_idx2(i, kt, pcols)] - cldeff[_idx2(i, kt - 1, pcols)], 0.0) * radinvfrac
            radf = qrlw[_idx2(i, kt, pcols)] / (pi[_idx2(i, kt, pcols)] - pi[_idx2(i, kt + 1, pcols)])
            radf = max(
                radinvfrac * radf * (zi[_idx2(i, kt, pcols)] - zi[_idx2(i, kt + 1, pcols)]),
                0.0,
            ) * chs[_idx2(i, kt, pcols)]

        else:
            lwp = ql[_idx2(i, kt, pcols)] * (pi[_idx2(i, kt + 1, pcols)] - pi[_idx2(i, kt, pcols)]) / gravit
            opt_depth = 156.0 * lwp
            radinvfrac = opt_depth * (4.0 + opt_depth) / (6.0 * (4.0 + opt_depth) + opt_depth**2)
            radf = max(
                radinvfrac
                * qrlw[_idx2(i, kt, pcols)]
                / (pi[_idx2(i, kt, pcols)] - pi[_idx2(i, kt + 1, pcols)])
                * (zi[_idx2(i, kt, pcols)] - zi[_idx2(i, kt + 1, pcols)]),
                0.0,
            )

            lwp = ql[_idx2(i, kt - 1, pcols)] * (pi[_idx2(i, kt, pcols)] - pi[_idx2(i, kt - 1, pcols)]) / gravit
            opt_depth = 156.0 * lwp
            radinvfrac = opt_depth * (4.0 + opt_depth) / (6.0 * (4.0 + opt_depth) + opt_depth**2)
            radf = radf + max(
                radinvfrac
                * qrlw[_idx2(i, kt - 1, pcols)]
                / (pi[_idx2(i, kt - 1, pcols)] - pi[_idx2(i, kt, pcols)])
                * (zi[_idx2(i, kt - 1, pcols)] - zi[_idx2(i, kt, pcols)]),
                0.0,
            )
            radf = max(radf, 0.0) * chs[_idx2(i, kt, pcols)]

        lwp_CL[_idx2(i, ncv, pcols)] = lwp
        opt_depth_CL[_idx2(i, ncv, pcols)] = opt_depth
        radinvfrac_CL[_idx2(i, ncv, pcols)] = radinvfrac
        radf_CL[_idx2(i, ncv, pcols)] = radf


@export
def eddy_diff_caleddy_clprep_codon(
    i_col: int,
    ncv_col: int,
    pcols: int,
    pver: int,
    ncvmax: int,
    tunl_mode: int,
    leng_mode: int,
    evhc_mode: int,
    tunl: float,
    ctunl: float,
    cleng: float,
    lbulk_max: float,
    qmin: float,
    gravit: float,
    vk: float,
    latvap: float,
    a2l: float,
    a3l: float,
    jbumin: float,
    evhcmax: float,
    ql_p: cobj,
    slv_p: cobj,
    sl_p: cobj,
    qt_p: cobj,
    u_p: cobj,
    v_p: cobj,
    zi_p: cobj,
    z_p: cobj,
    n2_p: cobj,
    s2_p: cobj,
    sfuh_p: cobj,
    sflh_p: cobj,
    chu_p: cobj,
    chs_p: cobj,
    cmu_p: cobj,
    cms_p: cobj,
    cldeff_p: cobj,
    bflxs_p: cobj,
    bprod_p: cobj,
    kbase_p: cobj,
    ktop_p: cobj,
    ricl_p: cobj,
    shcl_p: cobj,
    smcl_p: cobj,
    radf: float,
    leng_max_p: cobj,
    leng_p: cobj,
    wcap_p: cobj,
    clprep_state_p: cobj,
):
    ql = Ptr[float](ql_p)
    slv = Ptr[float](slv_p)
    sl = Ptr[float](sl_p)
    qt = Ptr[float](qt_p)
    u = Ptr[float](u_p)
    v = Ptr[float](v_p)
    zi = Ptr[float](zi_p)
    z = Ptr[float](z_p)
    n2 = Ptr[float](n2_p)
    s2 = Ptr[float](s2_p)
    sfuh = Ptr[float](sfuh_p)
    sflh = Ptr[float](sflh_p)
    chu = Ptr[float](chu_p)
    chs = Ptr[float](chs_p)
    cmu = Ptr[float](cmu_p)
    cms = Ptr[float](cms_p)
    cldeff = Ptr[float](cldeff_p)
    bflxs = Ptr[float](bflxs_p)
    bprod = Ptr[float](bprod_p)
    kbase = Ptr[i32](kbase_p)
    ktop = Ptr[i32](ktop_p)
    ricl = Ptr[float](ricl_p)
    shcl = Ptr[float](shcl_p)
    smcl = Ptr[float](smcl_p)
    leng_max = Ptr[float](leng_max_p)
    leng = Ptr[float](leng_p)
    wcap = Ptr[float](wcap_p)
    clprep_state = Ptr[float](clprep_state_p)

    i = i_col
    ncv = ncv_col
    kt = int(ktop[_idx2(i, ncv, pcols)])
    kb = int(kbase[_idx2(i, ncv, pcols)])

    if kb == pver + 1 and bflxs[i - 1] <= 0.0:
        lbulk = zi[_idx2(i, kt, pcols)] - z[_idx2(i, pver, pcols)]
    else:
        lbulk = zi[_idx2(i, kt, pcols)] - zi[_idx2(i, kb, pcols)]
    lbulk = min(lbulk, lbulk_max)

    for k in range(min(kb, pver), kt - 1, -1):
        if tunl_mode == 1:
            tunlramp = ctunl * tunl * (
                1.0
                - (1.0 - 1.0 / ctunl)
                * exp(min(0.0, ricl[_idx2(i, ncv, pcols)]))
            )
            tunlramp = min(max(tunlramp, tunl), ctunl * tunl)
        elif tunl_mode == 2:
            tunlramp = ctunl * tunl
        else:
            tunlramp = tunl

        if leng_mode == 0:
            leng[_idx2(i, k, pcols)] = (
                (vk * zi[_idx2(i, k, pcols)]) ** (-cleng)
                + (tunlramp * lbulk) ** (-cleng)
            ) ** (-1.0 / cleng)
        else:
            leng[_idx2(i, k, pcols)] = min(
                vk * zi[_idx2(i, k, pcols)], tunlramp * lbulk
            )
        leng[_idx2(i, k, pcols)] = min(leng_max[k - 1], leng[_idx2(i, k, pcols)])
        wcap[_idx2(i, k, pcols)] = (leng[_idx2(i, k, pcols)] ** 2) * (
            -shcl[_idx2(i, ncv, pcols)] * n2[_idx2(i, k, pcols)]
            + smcl[_idx2(i, ncv, pcols)] * s2[_idx2(i, k, pcols)]
        )

    if kb < pver + 1:
        jbzm = z[_idx2(i, kb - 1, pcols)] - z[_idx2(i, kb, pcols)]
        jbsl = sl[_idx2(i, kb - 1, pcols)] - sl[_idx2(i, kb, pcols)]
        jbqt = qt[_idx2(i, kb - 1, pcols)] - qt[_idx2(i, kb, pcols)]
        jbbu = n2[_idx2(i, kb, pcols)] * jbzm
        jbbu = max(jbbu, jbumin)
        jbu = u[_idx2(i, kb - 1, pcols)] - u[_idx2(i, kb, pcols)]
        jbv = v[_idx2(i, kb - 1, pcols)] - v[_idx2(i, kb, pcols)]
        ch = (1.0 - sflh[_idx2(i, kb - 1, pcols)]) * chu[_idx2(i, kb, pcols)] + sflh[
            _idx2(i, kb - 1, pcols)
        ] * chs[_idx2(i, kb, pcols)]
        cm = (1.0 - sflh[_idx2(i, kb - 1, pcols)]) * cmu[_idx2(i, kb, pcols)] + sflh[
            _idx2(i, kb - 1, pcols)
        ] * cms[_idx2(i, kb, pcols)]
        n2hb = (ch * jbsl + cm * jbqt) / jbzm
        vyb = n2hb * jbzm / jbbu
        vub = min(1.0, (jbu**2 + jbv**2) / (jbbu * jbzm))
    else:
        jbzm = 0.0
        jbbu = 0.0
        n2hb = 0.0
        vyb = 0.0
        vub = 0.0

    jtzm = z[_idx2(i, kt - 1, pcols)] - z[_idx2(i, kt, pcols)]
    jtsl = sl[_idx2(i, kt - 1, pcols)] - sl[_idx2(i, kt, pcols)]
    jtqt = qt[_idx2(i, kt - 1, pcols)] - qt[_idx2(i, kt, pcols)]
    jtbu = n2[_idx2(i, kt, pcols)] * jtzm
    jtbu = max(jtbu, jbumin)
    jtu = u[_idx2(i, kt - 1, pcols)] - u[_idx2(i, kt, pcols)]
    jtv = v[_idx2(i, kt - 1, pcols)] - v[_idx2(i, kt, pcols)]
    ch = (1.0 - sfuh[_idx2(i, kt, pcols)]) * chu[_idx2(i, kt, pcols)] + sfuh[
        _idx2(i, kt, pcols)
    ] * chs[_idx2(i, kt, pcols)]
    cm = (1.0 - sfuh[_idx2(i, kt, pcols)]) * cmu[_idx2(i, kt, pcols)] + sfuh[
        _idx2(i, kt, pcols)
    ] * cms[_idx2(i, kt, pcols)]
    n2ht = (ch * jtsl + cm * jtqt) / jtzm
    vyt = n2ht * jtzm / jtbu
    vut = min(1.0, (jtu**2 + jtv**2) / (jtbu * jtzm))

    evhc = 1.0
    jt2slv = 0.0

    if evhc_mode == 0:
        if ql[_idx2(i, kt, pcols)] > qmin and ql[_idx2(i, kt - 1, pcols)] < qmin:
            jt2slv = slv[_idx2(i, max(kt - 2, 1), pcols)] - slv[_idx2(i, kt, pcols)]
            jt2slv = max(jt2slv, jbumin * slv[_idx2(i, kt - 1, pcols)] / gravit)
            evhc = 1.0 + a2l * a3l * latvap * ql[_idx2(i, kt, pcols)] / jt2slv
            evhc = min(evhc, evhcmax)
    elif evhc_mode == 1:
        jt2slv = slv[_idx2(i, max(kt - 2, 1), pcols)] - slv[_idx2(i, kt, pcols)]
        jt2slv = max(jt2slv, jbumin * slv[_idx2(i, kt - 1, pcols)] / gravit)
        evhc = 1.0 + max(cldeff[_idx2(i, kt, pcols)] - cldeff[_idx2(i, kt - 1, pcols)], 0.0) * a2l * a3l * latvap * ql[
            _idx2(i, kt, pcols)
        ] / jt2slv
        evhc = min(evhc, evhcmax)
    else:
        qleff = max(ql[_idx2(i, kt - 1, pcols)], ql[_idx2(i, kt, pcols)])
        jt2slv = slv[_idx2(i, max(kt - 2, 1), pcols)] - slv[_idx2(i, kt, pcols)]
        jt2slv = max(jt2slv, jbumin * slv[_idx2(i, kt - 1, pcols)] / gravit)
        evhc = 1.0 + a2l * a3l * latvap * qleff / jt2slv
        evhc = min(evhc, evhcmax)

    dzht = zi[_idx2(i, kt, pcols)] - z[_idx2(i, kt, pcols)]
    dzhb = z[_idx2(i, kb - 1, pcols)] - zi[_idx2(i, kb, pcols)]
    wstar3 = radf * dzht
    for k in range(kt + 1, kb):
        wstar3 = wstar3 + bprod[_idx2(i, k, pcols)] * (
            z[_idx2(i, k - 1, pcols)] - z[_idx2(i, k, pcols)]
        )
    if kb == pver + 1 and bflxs[i - 1] > 0.0:
        wstar3 = wstar3 + bflxs[i - 1] * dzhb
    wstar3 = max(2.5 * wstar3, 0.0)

    # state: lbulk, jbzm, jbbu, n2hb, vyb, vub, jtzm, jtbu, jt2slv, n2ht,
    #        vyt, vut, evhc, dzht, dzhb, wstar3
    clprep_state[0] = lbulk
    clprep_state[1] = jbzm
    clprep_state[2] = jbbu
    clprep_state[3] = n2hb
    clprep_state[4] = vyb
    clprep_state[5] = vub
    clprep_state[6] = jtzm
    clprep_state[7] = jtbu
    clprep_state[8] = jt2slv
    clprep_state[9] = n2ht
    clprep_state[10] = vyt
    clprep_state[11] = vut
    clprep_state[12] = evhc
    clprep_state[13] = dzht
    clprep_state[14] = dzhb
    clprep_state[15] = wstar3


@inline
def _eddy_diff_compute_cubic(a: float, b: float, c: float) -> float:
    xmin = 1.0e-2
    qq = (a**2 - 3.0 * b) / 9.0
    rr = (2.0 * a**3 - 9.0 * a * b + 27.0 * c) / 54.0

    dd = rr**2 - qq**3
    if dd <= 0.0:
        theta = acos(rr / qq ** (3.0 / 2.0))
        x1 = -2.0 * sqrt(qq) * cos(theta / 3.0) - a / 3.0
        x2 = -2.0 * sqrt(qq) * cos((theta + 2.0 * 3.141592) / 3.0) - a / 3.0
        x3 = -2.0 * sqrt(qq) * cos((theta - 2.0 * 3.141592) / 3.0) - a / 3.0
        return max(max(max(x1, x2), x3), xmin)

    if rr >= 0.0:
        aa = -(sqrt(rr**2 - qq**3) + rr) ** (1.0 / 3.0)
    else:
        aa = (sqrt(rr**2 - qq**3) - rr) ** (1.0 / 3.0)

    if aa == 0.0:
        bb = 0.0
    else:
        bb = qq / aa

    return max((aa + bb) - a / 3.0, xmin)


@export
def eddy_diff_caleddy_closure_codon(
    i_col: int,
    pcols: int,
    pver: int,
    ncvmax: int,
    tunl_mode: int,
    leng_mode: int,
    evhc_mode: int,
    wstarent_mode: int,
    sedfact_mode: int,
    ncvsurf: int,
    tunl: float,
    ctunl: float,
    cleng: float,
    lbulk_max: float,
    tkemax: float,
    b1: float,
    ae: float,
    alph1: float,
    a1l: float,
    a1i: float,
    ccrit: float,
    wstar3factcrit: float,
    ntzero: float,
    onet: float,
    rcapmin: float,
    rcapmax: float,
    wfac: float,
    wpertmin: float,
    tfac: float,
    qmin: float,
    gravit: float,
    vk: float,
    cpair: float,
    latvap: float,
    a2l: float,
    a3l: float,
    jbumin: float,
    evhcmax: float,
    ased: float,
    ql_p: cobj,
    slv_p: cobj,
    sl_p: cobj,
    qt_p: cobj,
    u_p: cobj,
    v_p: cobj,
    pi_p: cobj,
    zi_p: cobj,
    z_p: cobj,
    n2_p: cobj,
    s2_p: cobj,
    shflx_p: cobj,
    qflx_p: cobj,
    rrho_p: cobj,
    sfuh_p: cobj,
    sflh_p: cobj,
    chu_p: cobj,
    chs_p: cobj,
    cmu_p: cobj,
    cms_p: cobj,
    cldeff_p: cobj,
    bflxs_p: cobj,
    bprod_p: cobj,
    sprod_p: cobj,
    wsedl_p: cobj,
    ncvfin_p: cobj,
    kbase_p: cobj,
    ktop_p: cobj,
    lbrk_p: cobj,
    ebrk_p: cobj,
    wbrk_p: cobj,
    ricl_p: cobj,
    shcl_p: cobj,
    smcl_p: cobj,
    radf_CL_p: cobj,
    wsed_CL_p: cobj,
    leng_max_p: cobj,
    wet_CL_p: cobj,
    web_CL_p: cobj,
    jtbu_CL_p: cobj,
    jbbu_CL_p: cobj,
    evhc_CL_p: cobj,
    jt2slv_CL_p: cobj,
    n2ht_CL_p: cobj,
    n2hb_CL_p: cobj,
    wstar_CL_p: cobj,
    wstar3fact_CL_p: cobj,
    leng_p: cobj,
    wcap_p: cobj,
    tke_p: cobj,
    kvh_p: cobj,
    kvm_p: cobj,
    turbtype_p: cobj,
    sm_aw_p: cobj,
    pblh_p: cobj,
    pblhp_p: cobj,
    wpert_p: cobj,
    tpert_p: cobj,
    qpert_p: cobj,
    ipbl_p: cobj,
    kpblh_p: cobj,
    went_p: cobj,
    zero_tke_mask_p: cobj,
    closure_status_p: cobj,
):
    ql = Ptr[float](ql_p)
    slv = Ptr[float](slv_p)
    sl = Ptr[float](sl_p)
    qt = Ptr[float](qt_p)
    u = Ptr[float](u_p)
    v = Ptr[float](v_p)
    pi = Ptr[float](pi_p)
    zi = Ptr[float](zi_p)
    z = Ptr[float](z_p)
    n2 = Ptr[float](n2_p)
    s2 = Ptr[float](s2_p)
    shflx = Ptr[float](shflx_p)
    qflx = Ptr[float](qflx_p)
    rrho = Ptr[float](rrho_p)
    sfuh = Ptr[float](sfuh_p)
    sflh = Ptr[float](sflh_p)
    chu = Ptr[float](chu_p)
    chs = Ptr[float](chs_p)
    cmu = Ptr[float](cmu_p)
    cms = Ptr[float](cms_p)
    cldeff = Ptr[float](cldeff_p)
    bflxs = Ptr[float](bflxs_p)
    bprod = Ptr[float](bprod_p)
    sprod = Ptr[float](sprod_p)
    wsedl = Ptr[float](wsedl_p)
    ncvfin = Ptr[i32](ncvfin_p)
    kbase = Ptr[i32](kbase_p)
    ktop = Ptr[i32](ktop_p)
    lbrk = Ptr[float](lbrk_p)
    ebrk = Ptr[float](ebrk_p)
    wbrk = Ptr[float](wbrk_p)
    ricl = Ptr[float](ricl_p)
    shcl = Ptr[float](shcl_p)
    smcl = Ptr[float](smcl_p)
    radf_CL = Ptr[float](radf_CL_p)
    wsed_CL = Ptr[float](wsed_CL_p)
    leng_max = Ptr[float](leng_max_p)
    wet_CL = Ptr[float](wet_CL_p)
    web_CL = Ptr[float](web_CL_p)
    jtbu_CL = Ptr[float](jtbu_CL_p)
    jbbu_CL = Ptr[float](jbbu_CL_p)
    evhc_CL = Ptr[float](evhc_CL_p)
    jt2slv_CL = Ptr[float](jt2slv_CL_p)
    n2ht_CL = Ptr[float](n2ht_CL_p)
    n2hb_CL = Ptr[float](n2hb_CL_p)
    wstar_CL = Ptr[float](wstar_CL_p)
    wstar3fact_CL = Ptr[float](wstar3fact_CL_p)
    leng = Ptr[float](leng_p)
    wcap = Ptr[float](wcap_p)
    tke = Ptr[float](tke_p)
    kvh = Ptr[float](kvh_p)
    kvm = Ptr[float](kvm_p)
    turbtype = Ptr[i32](turbtype_p)
    sm_aw = Ptr[float](sm_aw_p)
    pblh = Ptr[float](pblh_p)
    pblhp = Ptr[float](pblhp_p)
    wpert = Ptr[float](wpert_p)
    tpert = Ptr[float](tpert_p)
    qpert = Ptr[float](qpert_p)
    ipbl = Ptr[i32](ipbl_p)
    kpblh = Ptr[i32](kpblh_p)
    went = Ptr[float](went_p)
    zero_tke_mask = Ptr[i32](zero_tke_mask_p)
    closure_status = Ptr[i32](closure_status_p)

    i = i_col
    ktblw = 0

    for k in range(1, pver + 2):
        zero_tke_mask[k - 1] = i32(0)
    closure_status[0] = i32(0)
    closure_status[1] = i32(0)
    closure_status[2] = i32(0)

    ncvfin_i = int(ncvfin[i - 1])
    for ncv in range(1, ncvfin_i + 1):
        kt = int(ktop[_idx2(i, ncv, pcols)])
        kb = int(kbase[_idx2(i, ncv, pcols)])
        radf = radf_CL[_idx2(i, ncv, pcols)]

        if kb == pver + 1 and bflxs[i - 1] <= 0.0:
            lbulk = zi[_idx2(i, kt, pcols)] - z[_idx2(i, pver, pcols)]
        else:
            lbulk = zi[_idx2(i, kt, pcols)] - zi[_idx2(i, kb, pcols)]
        lbulk = min(lbulk, lbulk_max)

        for k in range(min(kb, pver), kt - 1, -1):
            if tunl_mode == 1:
                tunlramp = ctunl * tunl * (
                    1.0 - (1.0 - 1.0 / ctunl) * exp(min(0.0, ricl[_idx2(i, ncv, pcols)]))
                )
                tunlramp = min(max(tunlramp, tunl), ctunl * tunl)
            elif tunl_mode == 2:
                tunlramp = ctunl * tunl
            else:
                tunlramp = tunl

            if leng_mode == 0:
                leng[_idx2(i, k, pcols)] = (
                    (vk * zi[_idx2(i, k, pcols)]) ** (-cleng)
                    + (tunlramp * lbulk) ** (-cleng)
                ) ** (-1.0 / cleng)
            else:
                leng[_idx2(i, k, pcols)] = min(vk * zi[_idx2(i, k, pcols)], tunlramp * lbulk)
            leng[_idx2(i, k, pcols)] = min(leng_max[k - 1], leng[_idx2(i, k, pcols)])
            wcap[_idx2(i, k, pcols)] = (leng[_idx2(i, k, pcols)] ** 2) * (
                -shcl[_idx2(i, ncv, pcols)] * n2[_idx2(i, k, pcols)]
                + smcl[_idx2(i, ncv, pcols)] * s2[_idx2(i, k, pcols)]
            )

        if kb < pver + 1:
            jbzm = z[_idx2(i, kb - 1, pcols)] - z[_idx2(i, kb, pcols)]
            jbsl = sl[_idx2(i, kb - 1, pcols)] - sl[_idx2(i, kb, pcols)]
            jbqt = qt[_idx2(i, kb - 1, pcols)] - qt[_idx2(i, kb, pcols)]
            jbbu = n2[_idx2(i, kb, pcols)] * jbzm
            jbbu = max(jbbu, jbumin)
            jbu = u[_idx2(i, kb - 1, pcols)] - u[_idx2(i, kb, pcols)]
            jbv = v[_idx2(i, kb - 1, pcols)] - v[_idx2(i, kb, pcols)]
            ch = (1.0 - sflh[_idx2(i, kb - 1, pcols)]) * chu[_idx2(i, kb, pcols)] + sflh[
                _idx2(i, kb - 1, pcols)
            ] * chs[_idx2(i, kb, pcols)]
            cm = (1.0 - sflh[_idx2(i, kb - 1, pcols)]) * cmu[_idx2(i, kb, pcols)] + sflh[
                _idx2(i, kb - 1, pcols)
            ] * cms[_idx2(i, kb, pcols)]
            n2hb = (ch * jbsl + cm * jbqt) / jbzm
            vyb = n2hb * jbzm / jbbu
            vub = min(1.0, (jbu**2 + jbv**2) / (jbbu * jbzm))
        else:
            jbzm = 0.0
            jbbu = 0.0
            n2hb = 0.0
            vyb = 0.0
            vub = 0.0

        jtzm = z[_idx2(i, kt - 1, pcols)] - z[_idx2(i, kt, pcols)]
        jtsl = sl[_idx2(i, kt - 1, pcols)] - sl[_idx2(i, kt, pcols)]
        jtqt = qt[_idx2(i, kt - 1, pcols)] - qt[_idx2(i, kt, pcols)]
        jtbu = n2[_idx2(i, kt, pcols)] * jtzm
        jtbu = max(jtbu, jbumin)
        jtu = u[_idx2(i, kt - 1, pcols)] - u[_idx2(i, kt, pcols)]
        jtv = v[_idx2(i, kt - 1, pcols)] - v[_idx2(i, kt, pcols)]
        ch = (1.0 - sfuh[_idx2(i, kt, pcols)]) * chu[_idx2(i, kt, pcols)] + sfuh[
            _idx2(i, kt, pcols)
        ] * chs[_idx2(i, kt, pcols)]
        cm = (1.0 - sfuh[_idx2(i, kt, pcols)]) * cmu[_idx2(i, kt, pcols)] + sfuh[
            _idx2(i, kt, pcols)
        ] * cms[_idx2(i, kt, pcols)]
        n2ht = (ch * jtsl + cm * jtqt) / jtzm
        vyt = n2ht * jtzm / jtbu
        vut = min(1.0, (jtu**2 + jtv**2) / (jtbu * jtzm))

        evhc = 1.0
        jt2slv = 0.0
        if evhc_mode == 0:
            if ql[_idx2(i, kt, pcols)] > qmin and ql[_idx2(i, kt - 1, pcols)] < qmin:
                jt2slv = slv[_idx2(i, max(kt - 2, 1), pcols)] - slv[_idx2(i, kt, pcols)]
                jt2slv = max(jt2slv, jbumin * slv[_idx2(i, kt - 1, pcols)] / gravit)
                evhc = 1.0 + a2l * a3l * latvap * ql[_idx2(i, kt, pcols)] / jt2slv
                evhc = min(evhc, evhcmax)
        elif evhc_mode == 1:
            jt2slv = slv[_idx2(i, max(kt - 2, 1), pcols)] - slv[_idx2(i, kt, pcols)]
            jt2slv = max(jt2slv, jbumin * slv[_idx2(i, kt - 1, pcols)] / gravit)
            evhc = (
                1.0
                + max(cldeff[_idx2(i, kt, pcols)] - cldeff[_idx2(i, kt - 1, pcols)], 0.0)
                * a2l
                * a3l
                * latvap
                * ql[_idx2(i, kt, pcols)]
                / jt2slv
            )
            evhc = min(evhc, evhcmax)
        else:
            qleff = max(ql[_idx2(i, kt - 1, pcols)], ql[_idx2(i, kt, pcols)])
            jt2slv = slv[_idx2(i, max(kt - 2, 1), pcols)] - slv[_idx2(i, kt, pcols)]
            jt2slv = max(jt2slv, jbumin * slv[_idx2(i, kt - 1, pcols)] / gravit)
            evhc = 1.0 + a2l * a3l * latvap * qleff / jt2slv
            evhc = min(evhc, evhcmax)

        dzht = zi[_idx2(i, kt, pcols)] - z[_idx2(i, kt, pcols)]
        dzhb = z[_idx2(i, kb - 1, pcols)] - zi[_idx2(i, kb, pcols)]
        wstar3 = radf * dzht
        for k in range(kt + 1, kb):
            wstar3 = wstar3 + bprod[_idx2(i, k, pcols)] * (
                z[_idx2(i, k - 1, pcols)] - z[_idx2(i, k, pcols)]
            )
        if kb == pver + 1 and bflxs[i - 1] > 0.0:
            wstar3 = wstar3 + bflxs[i - 1] * dzhb
        wstar3 = max(2.5 * wstar3, 0.0)

        web = 0.0
        wstar = 0.0

        if sedfact_mode != 0:
            sedfact = exp(-ased * wsedl[_idx2(i, kt, pcols)] / (wstar3 ** (1.0 / 3.0) + 1.0e-6))
            wsed_CL[_idx2(i, ncv, pcols)] = wsedl[_idx2(i, kt, pcols)]
            if evhc_mode == 0:
                if ql[_idx2(i, kt, pcols)] > qmin and ql[_idx2(i, kt - 1, pcols)] < qmin:
                    jt2slv = slv[_idx2(i, max(kt - 2, 1), pcols)] - slv[_idx2(i, kt, pcols)]
                    jt2slv = max(jt2slv, jbumin * slv[_idx2(i, kt - 1, pcols)] / gravit)
                    evhc = 1.0 + sedfact * a2l * a3l * latvap * ql[_idx2(i, kt, pcols)] / jt2slv
                    evhc = min(evhc, evhcmax)
            elif evhc_mode == 1:
                jt2slv = slv[_idx2(i, max(kt - 2, 1), pcols)] - slv[_idx2(i, kt, pcols)]
                jt2slv = max(jt2slv, jbumin * slv[_idx2(i, kt - 1, pcols)] / gravit)
                evhc = (
                    1.0
                    + max(cldeff[_idx2(i, kt, pcols)] - cldeff[_idx2(i, kt - 1, pcols)], 0.0)
                    * sedfact
                    * a2l
                    * a3l
                    * latvap
                    * ql[_idx2(i, kt, pcols)]
                    / jt2slv
                )
                evhc = min(evhc, evhcmax)
            else:
                qleff = max(ql[_idx2(i, kt - 1, pcols)], ql[_idx2(i, kt, pcols)])
                jt2slv = slv[_idx2(i, max(kt - 2, 1), pcols)] - slv[_idx2(i, kt, pcols)]
                jt2slv = max(jt2slv, jbumin * slv[_idx2(i, kt - 1, pcols)] / gravit)
                evhc = 1.0 + sedfact * a2l * a3l * latvap * qleff / jt2slv
                evhc = min(evhc, evhcmax)

        if wstar3 > 0.0:
            cet = a1i * evhc / (jtbu * lbulk)
            if kb == pver + 1:
                wstar3fact = max(1.0 + 2.5 * cet * n2ht * jtzm * dzht, wstar3factcrit)
            else:
                ceb = a1i / (jbbu * lbulk)
                wstar3fact = max(
                    1.0 + 2.5 * cet * n2ht * jtzm * dzht + 2.5 * ceb * n2hb * jbzm * dzhb,
                    wstar3factcrit,
                )
            wstar3 = wstar3 / wstar3fact
        else:
            wstar3fact = 0.0
            cet = 0.0
            ceb = 0.0

        fact = (evhc * (-vyt + vut) * dzht + (-vyb + vub) * dzhb * leng[_idx2(i, kb, pcols)] / leng[_idx2(i, kt, pcols)]) / lbulk

        if wstarent_mode != 0:
            trma = 1.0
            trmp = ebrk[_idx2(i, ncv, pcols)] * (lbrk[_idx2(i, ncv, pcols)] / lbulk) / 3.0 + ntzero
            trmq = 0.5 * b1 * (leng[_idx2(i, kt, pcols)] / lbulk) * (radf * dzht + a1i * fact * wstar3)

            rmin = sqrt(trmp)
            fmin = rmin * (rmin * rmin - 3.0 * trmp) - 2.0 * trmq
            wstar = wstar3**onet
            rcrit = ccrit * wstar
            fcrit = rcrit * (rcrit * rcrit - 3.0 * trmp) - 2.0 * trmq
            noroot = ((rmin < rcrit) and (fcrit > 0.0)) or ((rmin >= rcrit) and (fmin > 0.0))
            if noroot:
                trma = 1.0 - b1 * (leng[_idx2(i, kt, pcols)] / lbulk) * a1i * fact / ccrit**3
                trma = max(trma, 0.5)
                trmp = trmp / trma
                trmq = 0.5 * b1 * (leng[_idx2(i, kt, pcols)] / lbulk) * radf * dzht / trma

            qq = trmq**2 - trmp**3
            if qq >= 0.0:
                rootp = (trmq + sqrt(qq)) ** (1.0 / 3.0) + (max(trmq - sqrt(qq), 0.0)) ** (1.0 / 3.0)
            else:
                rootp = 2.0 * sqrt(trmp) * cos(acos(trmq / sqrt(trmp**3)) / 3.0)

            if noroot:
                wstar3 = (rootp / ccrit) ** 3
            wet = cet * wstar3
            if kb < pver + 1:
                web = ceb * wstar3
        else:
            trma = 1.0 - b1 * a1l * fact
            trma = max(trma, 0.5)
            trmp = ebrk[_idx2(i, ncv, pcols)] * (lbrk[_idx2(i, ncv, pcols)] / lbulk) / (3.0 * trma)
            trmq = 0.5 * b1 * (leng[_idx2(i, kt, pcols)] / lbulk) * radf * dzht / trma

            qq = trmq**2 - trmp**3
            if qq >= 0.0:
                rootp = (trmq + sqrt(qq)) ** (1.0 / 3.0) + (max(trmq - sqrt(qq), 0.0)) ** (1.0 / 3.0)
            else:
                rootp = 2.0 * sqrt(trmp) * cos(acos(trmq / sqrt(trmp**3)) / 3.0)

            wet = a1l * rootp * min(evhc * rootp**2 / (leng[_idx2(i, kt, pcols)] * jtbu), 1.0)
            if kb < pver + 1:
                web = a1l * rootp * min(evhc * rootp**2 / (leng[_idx2(i, kb, pcols)] * jbbu), 1.0)

        ebrk[_idx2(i, ncv, pcols)] = rootp**2
        ebrk[_idx2(i, ncv, pcols)] = min(ebrk[_idx2(i, ncv, pcols)], tkemax)
        wbrk[_idx2(i, ncv, pcols)] = ebrk[_idx2(i, ncv, pcols)] / b1

        if ebrk[_idx2(i, ncv, pcols)] <= 0.0:
            closure_status[0] = i32(1)
            closure_status[1] = i32(kt)
            closure_status[2] = i32(kb)
            zero_tke_mask[kt - 1] = i32(1)
            zero_tke_mask[kb - 1] = i32(1)

        for k in range(kb - 1, kt, -1):
            rcap = (b1 * ae + wcap[_idx2(i, k, pcols)] / wbrk[_idx2(i, ncv, pcols)]) / (b1 * ae + 1.0)
            rcap = min(max(rcap, rcapmin), rcapmax)
            tke[_idx2(i, k, pcols)] = ebrk[_idx2(i, ncv, pcols)] * rcap
            tke[_idx2(i, k, pcols)] = min(tke[_idx2(i, k, pcols)], tkemax)
            kvh[_idx2(i, k, pcols)] = leng[_idx2(i, k, pcols)] * sqrt(tke[_idx2(i, k, pcols)]) * shcl[_idx2(i, ncv, pcols)]
            kvm[_idx2(i, k, pcols)] = leng[_idx2(i, k, pcols)] * sqrt(tke[_idx2(i, k, pcols)]) * smcl[_idx2(i, ncv, pcols)]
            bprod[_idx2(i, k, pcols)] = -kvh[_idx2(i, k, pcols)] * n2[_idx2(i, k, pcols)]
            sprod[_idx2(i, k, pcols)] = kvm[_idx2(i, k, pcols)] * s2[_idx2(i, k, pcols)]
            turbtype[_idx2(i, k, pcols)] = i32(2)
            sm_aw[_idx2(i, k, pcols)] = smcl[_idx2(i, ncv, pcols)] / alph1

        kentr = wet * jtzm
        kvh[_idx2(i, kt, pcols)] = kentr
        kvm[_idx2(i, kt, pcols)] = kentr
        bprod[_idx2(i, kt, pcols)] = -kentr * n2ht + radf
        sprod[_idx2(i, kt, pcols)] = kentr * s2[_idx2(i, kt, pcols)]
        turbtype[_idx2(i, kt, pcols)] = i32(4)
        trmp = -b1 * ae / (1.0 + b1 * ae)
        trmq = -(bprod[_idx2(i, kt, pcols)] + sprod[_idx2(i, kt, pcols)]) * b1 * leng[_idx2(i, kt, pcols)] / (
            1.0 + b1 * ae
        ) / (ebrk[_idx2(i, ncv, pcols)] ** (3.0 / 2.0))
        rcap = _eddy_diff_compute_cubic(0.0, trmp, trmq) ** 2
        rcap = min(max(rcap, rcapmin), rcapmax)
        tke[_idx2(i, kt, pcols)] = ebrk[_idx2(i, ncv, pcols)] * rcap
        tke[_idx2(i, kt, pcols)] = min(tke[_idx2(i, kt, pcols)], tkemax)
        sm_aw[_idx2(i, kt, pcols)] = smcl[_idx2(i, ncv, pcols)] / alph1

        if kb < pver + 1:
            kentr = web * jbzm
            if kb != ktblw:
                kvh[_idx2(i, kb, pcols)] = kentr
                kvm[_idx2(i, kb, pcols)] = kentr
                bprod[_idx2(i, kb, pcols)] = -kvh[_idx2(i, kb, pcols)] * n2hb
                sprod[_idx2(i, kb, pcols)] = kvm[_idx2(i, kb, pcols)] * s2[_idx2(i, kb, pcols)]
                turbtype[_idx2(i, kb, pcols)] = i32(3)
                trmp = -b1 * ae / (1.0 + b1 * ae)
                trmq = -(bprod[_idx2(i, kb, pcols)] + sprod[_idx2(i, kb, pcols)]) * b1 * leng[_idx2(i, kb, pcols)] / (
                    1.0 + b1 * ae
                ) / (ebrk[_idx2(i, ncv, pcols)] ** (3.0 / 2.0))
                rcap = _eddy_diff_compute_cubic(0.0, trmp, trmq) ** 2
                rcap = min(max(rcap, rcapmin), rcapmax)
                tke[_idx2(i, kb, pcols)] = ebrk[_idx2(i, ncv, pcols)] * rcap
                tke[_idx2(i, kb, pcols)] = min(tke[_idx2(i, kb, pcols)], tkemax)
            else:
                kvh[_idx2(i, kb, pcols)] = kvh[_idx2(i, kb, pcols)] + kentr
                kvm[_idx2(i, kb, pcols)] = kvm[_idx2(i, kb, pcols)] + kentr
                dzhb5 = z[_idx2(i, kb - 1, pcols)] - zi[_idx2(i, kb, pcols)]
                dzht5 = zi[_idx2(i, kb, pcols)] - z[_idx2(i, kb, pcols)]
                bprod[_idx2(i, kb, pcols)] = (
                    dzht5 * bprod[_idx2(i, kb, pcols)] - dzhb5 * kentr * n2hb
                ) / (dzhb5 + dzht5)
                sprod[_idx2(i, kb, pcols)] = (
                    dzht5 * sprod[_idx2(i, kb, pcols)] + dzhb5 * kentr * s2[_idx2(i, kb, pcols)]
                ) / (dzhb5 + dzht5)
                trmp = -b1 * ae / (1.0 + b1 * ae)
                trmq = -kentr * (s2[_idx2(i, kb, pcols)] - n2hb) * b1 * leng[_idx2(i, kb, pcols)] / (
                    1.0 + b1 * ae
                ) / (ebrk[_idx2(i, ncv, pcols)] ** (3.0 / 2.0))
                rcap = _eddy_diff_compute_cubic(0.0, trmp, trmq) ** 2
                rcap = min(max(rcap, rcapmin), rcapmax)
                tke_imsi = ebrk[_idx2(i, ncv, pcols)] * rcap
                tke_imsi = min(tke_imsi, tkemax)
                tke[_idx2(i, kb, pcols)] = (
                    dzht5 * tke[_idx2(i, kb, pcols)] + dzhb5 * tke_imsi
                ) / (dzhb5 + dzht5)
                tke[_idx2(i, kb, pcols)] = min(tke[_idx2(i, kb, pcols)], tkemax)
                turbtype[_idx2(i, kb, pcols)] = i32(5)
        else:
            rcap = (b1 * ae + wcap[_idx2(i, kb, pcols)] / wbrk[_idx2(i, ncv, pcols)]) / (b1 * ae + 1.0)
            rcap = min(max(rcap, rcapmin), rcapmax)
            tke[_idx2(i, kb, pcols)] = ebrk[_idx2(i, ncv, pcols)] * rcap
            tke[_idx2(i, kb, pcols)] = min(tke[_idx2(i, kb, pcols)], tkemax)

        sm_aw[_idx2(i, kb, pcols)] = smcl[_idx2(i, ncv, pcols)] / alph1
        wcap[_idx2(i, kt, pcols)] = (bprod[_idx2(i, kt, pcols)] + sprod[_idx2(i, kt, pcols)]) * leng[
            _idx2(i, kt, pcols)
        ] / sqrt(max(tke[_idx2(i, kt, pcols)], 1.0e-6))
        if kb < pver + 1:
            wcap[_idx2(i, kb, pcols)] = (bprod[_idx2(i, kb, pcols)] + sprod[_idx2(i, kb, pcols)]) * leng[
                _idx2(i, kb, pcols)
            ] / sqrt(max(tke[_idx2(i, kb, pcols)], 1.0e-6))

        ktblw = kt
        wet_CL[_idx2(i, ncv, pcols)] = wet
        web_CL[_idx2(i, ncv, pcols)] = web
        jtbu_CL[_idx2(i, ncv, pcols)] = jtbu
        jbbu_CL[_idx2(i, ncv, pcols)] = jbbu
        evhc_CL[_idx2(i, ncv, pcols)] = evhc
        jt2slv_CL[_idx2(i, ncv, pcols)] = jt2slv
        n2ht_CL[_idx2(i, ncv, pcols)] = n2ht
        n2hb_CL[_idx2(i, ncv, pcols)] = n2hb
        wstar_CL[_idx2(i, ncv, pcols)] = wstar
        wstar3fact_CL[_idx2(i, ncv, pcols)] = wstar3fact

    if ncvsurf > 0:
        ktopbl_local = int(ktop[_idx2(i, ncvsurf, pcols)])
        pblh[i - 1] = zi[_idx2(i, ktopbl_local, pcols)]
        pblhp[i - 1] = pi[_idx2(i, ktopbl_local, pcols)]
        wpert[i - 1] = max(wfac * sqrt(ebrk[_idx2(i, ncvsurf, pcols)]), wpertmin)
        tpert[i - 1] = max(abs(shflx[i - 1] * rrho[i - 1] / cpair) * tfac / wpert[i - 1], 0.0)
        qpert[i - 1] = max(abs(qflx[i - 1] * rrho[i - 1]) * tfac / wpert[i - 1], 0.0)
        if bflxs[i - 1] > 0.0:
            turbtype[_idx2(i, pver + 1, pcols)] = i32(2)
        else:
            turbtype[_idx2(i, pver + 1, pcols)] = i32(3)
        ipbl[i - 1] = i32(1)
        kpblh[i - 1] = i32(max(ktopbl_local - 1, 1))
        went[i - 1] = wet_CL[_idx2(i, ncvsurf, pcols)]


@export
def eddy_diff_caleddy_srcl_codon(
    i_col: int,
    pcols: int,
    pver: int,
    ncvmax: int,
    ntop_turb: int,
    nbot_turb: int,
    srcl_mode: int,
    qmin: float,
    ricrit: float,
    b1: float,
    vk: float,
    alph1: float,
    alph2: float,
    alph3: float,
    alph4exs: float,
    alph5: float,
    ghmin: float,
    ql_p: cobj,
    qrlw_p: cobj,
    ri_p: cobj,
    sfuh_p: cobj,
    chu_p: cobj,
    chs_p: cobj,
    cmu_p: cobj,
    cms_p: cobj,
    slslope_p: cobj,
    qtslope_p: cobj,
    z_p: cobj,
    bflxs_p: cobj,
    tkes_p: cobj,
    bprod_p: cobj,
    sprod_p: cobj,
    ncvfin_p: cobj,
    kbase_p: cobj,
    ktop_p: cobj,
    ricl_p: cobj,
    ghcl_p: cobj,
    shcl_p: cobj,
    smcl_p: cobj,
    lbrk_p: cobj,
    wbrk_p: cobj,
    ebrk_p: cobj,
    belong_mask_p: cobj,
    ncvsurf_p: cobj,
    srcl_status_p: cobj,
):
    ql = Ptr[float](ql_p)
    qrlw = Ptr[float](qrlw_p)
    ri = Ptr[float](ri_p)
    sfuh = Ptr[float](sfuh_p)
    chu = Ptr[float](chu_p)
    chs = Ptr[float](chs_p)
    cmu = Ptr[float](cmu_p)
    cms = Ptr[float](cms_p)
    slslope = Ptr[float](slslope_p)
    qtslope = Ptr[float](qtslope_p)
    z = Ptr[float](z_p)
    bflxs = Ptr[float](bflxs_p)
    tkes = Ptr[float](tkes_p)
    bprod = Ptr[float](bprod_p)
    sprod = Ptr[float](sprod_p)
    ncvfin = Ptr[i32](ncvfin_p)
    kbase = Ptr[i32](kbase_p)
    ktop = Ptr[i32](ktop_p)
    ricl = Ptr[float](ricl_p)
    ghcl = Ptr[float](ghcl_p)
    shcl = Ptr[float](shcl_p)
    smcl = Ptr[float](smcl_p)
    lbrk = Ptr[float](lbrk_p)
    wbrk = Ptr[float](wbrk_p)
    ebrk = Ptr[float](ebrk_p)
    belong_mask = Ptr[i32](belong_mask_p)
    ncvsurf = Ptr[i32](ncvsurf_p)
    srcl_status = Ptr[i32](srcl_status_p)

    i = i_col
    ncv = 1
    ncvf = int(ncvfin[i - 1])
    srcl_status[0] = i32(0)

    for k in range(1, pver + 2):
        belong_mask[k - 1] = i32(0)

    for ncv_mask in range(1, ncvf + 1):
        kt_mask = int(ktop[_idx2(i, ncv_mask, pcols)])
        kb_mask = int(kbase[_idx2(i, ncv_mask, pcols)])
        for k in range(kt_mask, kb_mask + 1):
            belong_mask[k - 1] = i32(1)

    if srcl_mode == 0:
        return

    for k in range(nbot_turb, ntop_turb, -1):
        if (
            ql[_idx2(i, k, pcols)] > qmin
            and ql[_idx2(i, k - 1, pcols)] < qmin
            and qrlw[_idx2(i, k, pcols)] < 0.0
            and ri[_idx2(i, k, pcols)] >= ricrit
        ):
            if srcl_mode == 2 and belong_mask[k] != i32(0):
                continue

            ch = (1.0 - sfuh[_idx2(i, k, pcols)]) * chu[_idx2(i, k, pcols)] + sfuh[
                _idx2(i, k, pcols)
            ] * chs[_idx2(i, k, pcols)]
            cm = (1.0 - sfuh[_idx2(i, k, pcols)]) * cmu[_idx2(i, k, pcols)] + sfuh[
                _idx2(i, k, pcols)
            ] * cms[_idx2(i, k, pcols)]

            n2htSRCL = ch * slslope[_idx2(i, k, pcols)] + cm * qtslope[_idx2(i, k, pcols)]

            if n2htSRCL <= 0.0:
                in_CL = False

                while ncv <= ncvf:
                    if int(ktop[_idx2(i, ncv, pcols)]) <= k:
                        if int(kbase[_idx2(i, ncv, pcols)]) > k:
                            in_CL = True
                        break
                    ncv += 1

                if not in_CL:
                    ncvnew = int(ncvfin[i - 1]) + 1
                    ncvfin[i - 1] = i32(ncvnew)
                    ktop[_idx2(i, ncvnew, pcols)] = i32(k)
                    kbase[_idx2(i, ncvnew, pcols)] = i32(k + 1)
                    belong_mask[k - 1] = i32(1)
                    belong_mask[k] = i32(1)

                    if k < pver:
                        wbrk[_idx2(i, ncvnew, pcols)] = 0.0
                        ebrk[_idx2(i, ncvnew, pcols)] = 0.0
                        lbrk[_idx2(i, ncvnew, pcols)] = 0.0
                        ghcl[_idx2(i, ncvnew, pcols)] = 0.0
                        shcl[_idx2(i, ncvnew, pcols)] = 0.0
                        smcl[_idx2(i, ncvnew, pcols)] = 0.0
                        ricl[_idx2(i, ncvnew, pcols)] = 0.0
                    else:
                        if bflxs[i - 1] > 0.0:
                            ebrk[_idx2(i, ncvnew, pcols)] = tkes[i - 1]
                            lbrk[_idx2(i, ncvnew, pcols)] = z[_idx2(i, pver, pcols)]
                            wbrk[_idx2(i, ncvnew, pcols)] = tkes[i - 1] / b1
                            srcl_status[0] = i32(1)
                            return
                        else:
                            ebrk[_idx2(i, ncvnew, pcols)] = 0.0
                            lbrk[_idx2(i, ncvnew, pcols)] = 0.0
                            wbrk[_idx2(i, ncvnew, pcols)] = 0.0

                        gg = (
                            0.5
                            * vk
                            * z[_idx2(i, pver, pcols)]
                            * bprod[_idx2(i, pver + 1, pcols)]
                            / (tkes[i - 1] ** (3.0 / 2.0))
                        )
                        if abs(alph5 - gg * alph3) <= 1.0e-7:
                            gh = ghmin
                        else:
                            gh = gg / (alph5 - gg * alph3)
                        gh = min(max(gh, ghmin), 0.0233)
                        ghcl[_idx2(i, ncvnew, pcols)] = gh
                        shcl[_idx2(i, ncvnew, pcols)] = max(0.0, alph5 / (1.0 + alph3 * gh))
                        smcl[_idx2(i, ncvnew, pcols)] = max(
                            0.0,
                            (alph1 + alph2 * gh)
                            / (1.0 + alph3 * gh)
                            / (1.0 + alph4exs * gh),
                        )
                        ricl[_idx2(i, ncvnew, pcols)] = -(
                            smcl[_idx2(i, ncvnew, pcols)] / shcl[_idx2(i, ncvnew, pcols)]
                        ) * (
                            bprod[_idx2(i, pver + 1, pcols)] / sprod[_idx2(i, pver + 1, pcols)]
                        )
                        ncvsurf[0] = i32(ncvnew)


@inline
def _eddy_diff_tunlramp_stl(tunl_mode: int, ctunl: float, tunl: float, ri_val: float, ricrit: float) -> float:
    if tunl_mode == 2:
        return max(1.0e-3, ctunl * tunl * exp(-log(ctunl) * ri_val / ricrit))
    return tunl


@inline
def _eddy_diff_leng_stl(
    leng_mode: int,
    vk: float,
    zi_val: float,
    tunlramp: float,
    lbulk: float,
    cleng: float,
    leng_max_val: float,
) -> float:
    if leng_mode == 0:
        leng_val = ((vk * zi_val) ** (-cleng) + (tunlramp * lbulk) ** (-cleng)) ** (-1.0 / cleng)
    else:
        leng_val = min(vk * zi_val, tunlramp * lbulk)
    return min(leng_max_val, leng_val)


@export
def eddy_diff_caleddy_stl_codon(
    i_col: int,
    pcols: int,
    pver: int,
    ncvmax: int,
    tunl_mode: int,
    leng_mode: int,
    ricrit: float,
    tunl: float,
    ctunl: float,
    cleng: float,
    lbulk_max: float,
    tkemax: float,
    b1: float,
    ae: float,
    alph1: float,
    alph2: float,
    alph3: float,
    alph4exs: float,
    alph5: float,
    ghmin: float,
    vk: float,
    fak: float,
    cpair: float,
    ri_p: cobj,
    z_p: cobj,
    zi_p: cobj,
    pi_p: cobj,
    n2_p: cobj,
    s2_p: cobj,
    shflx_p: cobj,
    qflx_p: cobj,
    rrho_p: cobj,
    ustar_p: cobj,
    leng_max_p: cobj,
    ncvfin_p: cobj,
    ktop_p: cobj,
    kbase_p: cobj,
    kvh_p: cobj,
    kvm_p: cobj,
    leng_p: cobj,
    tke_p: cobj,
    wcap_p: cobj,
    bprod_p: cobj,
    sprod_p: cobj,
    turbtype_p: cobj,
    sm_aw_p: cobj,
    pblh_p: cobj,
    pblhp_p: cobj,
    wpert_p: cobj,
    tpert_p: cobj,
    qpert_p: cobj,
    ipbl_p: cobj,
    kpblh_p: cobj,
    clmask_p: cobj,
    stlmask_p: cobj,
):
    ri = Ptr[float](ri_p)
    z = Ptr[float](z_p)
    zi = Ptr[float](zi_p)
    pi = Ptr[float](pi_p)
    n2 = Ptr[float](n2_p)
    s2 = Ptr[float](s2_p)
    shflx = Ptr[float](shflx_p)
    qflx = Ptr[float](qflx_p)
    rrho = Ptr[float](rrho_p)
    ustar = Ptr[float](ustar_p)
    leng_max = Ptr[float](leng_max_p)
    ncvfin = Ptr[i32](ncvfin_p)
    ktop = Ptr[i32](ktop_p)
    kbase = Ptr[i32](kbase_p)
    kvh = Ptr[float](kvh_p)
    kvm = Ptr[float](kvm_p)
    leng = Ptr[float](leng_p)
    tke = Ptr[float](tke_p)
    wcap = Ptr[float](wcap_p)
    bprod = Ptr[float](bprod_p)
    sprod = Ptr[float](sprod_p)
    turbtype = Ptr[i32](turbtype_p)
    sm_aw = Ptr[float](sm_aw_p)
    pblh = Ptr[float](pblh_p)
    pblhp = Ptr[float](pblhp_p)
    wpert = Ptr[float](wpert_p)
    tpert = Ptr[float](tpert_p)
    qpert = Ptr[float](qpert_p)
    ipbl = Ptr[i32](ipbl_p)
    kpblh = Ptr[i32](kpblh_p)
    clmask = Ptr[i32](clmask_p)
    stlmask = Ptr[i32](stlmask_p)

    i = i_col
    ncvfin_i = int(ncvfin[i - 1])

    for k in range(1, pver + 2):
        clmask[k - 1] = i32(0)
        stlmask[k - 1] = i32(0)

    for ncv in range(1, ncvfin_i + 1):
        kt = int(ktop[_idx2(i, ncv, pcols)])
        kb = int(kbase[_idx2(i, ncv, pcols)])
        for k in range(kt, kb + 1):
            clmask[k - 1] = i32(1)

    stlmask[0] = i32(0)
    for k in range(2, pver + 1):
        idx = _idx2(i, k, pcols)
        if ri[idx] < ricrit and clmask[k - 1] == i32(0):
            stlmask[k - 1] = i32(1)
        else:
            stlmask[k - 1] = i32(0)

        if stlmask[k - 1] != i32(0) and stlmask[k - 2] == i32(0):
            kt = k
        elif stlmask[k - 1] == i32(0) and stlmask[k - 2] != i32(0):
            kb = k - 1
            lbulk = z[_idx2(i, kt - 1, pcols)] - z[_idx2(i, kb, pcols)]
            lbulk = min(lbulk, lbulk_max)
            for ks in range(kt, kb + 1):
                tunlramp = _eddy_diff_tunlramp_stl(
                    tunl_mode,
                    ctunl,
                    tunl,
                    ri[_idx2(i, ks, pcols)],
                    ricrit,
                )
                leng[_idx2(i, ks, pcols)] = _eddy_diff_leng_stl(
                    leng_mode,
                    vk,
                    zi[_idx2(i, ks, pcols)],
                    tunlramp,
                    lbulk,
                    cleng,
                    leng_max[ks - 1],
                )

    if clmask[pver] == i32(0):
        stlmask[pver] = i32(1)
    else:
        stlmask[pver] = i32(0)

    if stlmask[pver] != i32(0):
        turbtype[_idx2(i, pver + 1, pcols)] = i32(1)

        if stlmask[pver - 1] != i32(0):
            lbulk = z[_idx2(i, kt - 1, pcols)]
        else:
            kt = pver + 1
            lbulk = z[_idx2(i, kt - 1, pcols)]
        lbulk = min(lbulk, lbulk_max)

        ktopbl = kt - 1
        pblh[i - 1] = z[_idx2(i, ktopbl, pcols)]
        pblhp[i - 1] = 0.5 * (pi[_idx2(i, ktopbl, pcols)] + pi[_idx2(i, ktopbl + 1, pcols)])

        for ks in range(kt, pver + 1):
            tunlramp = _eddy_diff_tunlramp_stl(
                tunl_mode,
                ctunl,
                tunl,
                ri[_idx2(i, ks, pcols)],
                ricrit,
            )
            leng[_idx2(i, ks, pcols)] = _eddy_diff_leng_stl(
                leng_mode,
                vk,
                zi[_idx2(i, ks, pcols)],
                tunlramp,
                lbulk,
                cleng,
                leng_max[ks - 1],
            )

        wpert[i - 1] = 0.0
        tpert[i - 1] = max(shflx[i - 1] * rrho[i - 1] / cpair * fak / ustar[i - 1], 0.0)
        qpert[i - 1] = max(qflx[i - 1] * rrho[i - 1] * fak / ustar[i - 1], 0.0)

        ipbl[i - 1] = i32(0)
        kpblh[i - 1] = i32(ktopbl)

    for k in range(2, pver + 1):
        if stlmask[k - 1] != i32(0):
            idx = _idx2(i, k, pcols)
            turbtype[idx] = i32(1)
            trma = alph3 * alph4exs * ri[idx] + 2.0 * b1 * (alph2 - alph4exs * alph5 * ri[idx])
            trmb = (alph3 + alph4exs) * ri[idx] + 2.0 * b1 * (-alph5 * ri[idx] + alph1)
            trmc = ri[idx]
            det = max(trmb * trmb - 4.0 * trma * trmc, 0.0)
            gh = (-trmb + sqrt(det)) / (2.0 * trma)
            gh = min(max(gh, ghmin), 0.0233)
            sh = max(0.0, alph5 / (1.0 + alph3 * gh))
            sm = max(0.0, (alph1 + alph2 * gh) / (1.0 + alph3 * gh) / (1.0 + alph4exs * gh))

            tke[idx] = b1 * (leng[idx] ** 2) * (-sh * n2[idx] + sm * s2[idx])
            tke[idx] = min(tke[idx], tkemax)
            wcap[idx] = tke[idx] / b1
            kvh[idx] = leng[idx] * sqrt(tke[idx]) * sh
            kvm[idx] = leng[idx] * sqrt(tke[idx]) * sm
            bprod[idx] = -kvh[idx] * n2[idx]
            sprod[idx] = kvm[idx] * s2[idx]
            sm_aw[idx] = sm / alph1

    for k in range(2, pver + 1):
        idx = _idx2(i, k, pcols)
        if turbtype[idx] == i32(3) or turbtype[idx] == i32(4) or turbtype[idx] == i32(5):
            trma = alph3 * alph4exs * ri[idx] + 2.0 * b1 * (alph2 - alph4exs * alph5 * ri[idx])
            trmb = (alph3 + alph4exs) * ri[idx] + 2.0 * b1 * (-alph5 * ri[idx] + alph1)
            trmc = ri[idx]
            det = max(trmb * trmb - 4.0 * trma * trmc, 0.0)
            gh = (-trmb + sqrt(det)) / (2.0 * trma)
            gh = min(max(gh, ghmin), 0.0233)
            sh = max(0.0, alph5 / (1.0 + alph3 * gh))
            sm = max(0.0, (alph1 + alph2 * gh) / (1.0 + alph3 * gh) / (1.0 + alph4exs * gh))

            lbulk = z[_idx2(i, k - 1, pcols)] - z[_idx2(i, k, pcols)]
            lbulk = min(lbulk, lbulk_max)
            tunlramp = _eddy_diff_tunlramp_stl(tunl_mode, ctunl, tunl, ri[idx], ricrit)
            leng_imsi = _eddy_diff_leng_stl(
                leng_mode,
                vk,
                zi[_idx2(i, k, pcols)],
                tunlramp,
                lbulk,
                cleng,
                leng_max[k - 1],
            )

            tke_imsi = b1 * (leng_imsi ** 2) * (-sh * n2[idx] + sm * s2[idx])
            tke_imsi = min(max(tke_imsi, 0.0), tkemax)
            kvh_imsi = leng_imsi * sqrt(tke_imsi) * sh
            kvm_imsi = leng_imsi * sqrt(tke_imsi) * sm

            if kvh[idx] < kvh_imsi:
                kvh[idx] = kvh_imsi
                kvm[idx] = kvm_imsi
                leng[idx] = leng_imsi
                tke[idx] = tke_imsi
                wcap[idx] = tke_imsi / b1
                bprod[idx] = -kvh_imsi * n2[idx]
                sprod[idx] = kvm_imsi * s2[idx]
                sm_aw[idx] = sm / alph1
                turbtype[idx] = i32(1)


@export
def eddy_diff_caleddy_diag_codon(
    i_col: int,
    pcols: int,
    pver: int,
    ricrit: float,
    b1: float,
    alph1: float,
    alph2: float,
    alph3: float,
    alph4: float,
    alph4exs: float,
    alph5: float,
    ghmin: float,
    vk: float,
    tkes_p: cobj,
    z_p: cobj,
    ri_p: cobj,
    bflxs_p: cobj,
    bprod_p: cobj,
    sprod_p: cobj,
    gh_a_p: cobj,
    sh_a_p: cobj,
    sm_a_p: cobj,
    ri_a_p: cobj,
    sm_aw_p: cobj,
):
    tkes = Ptr[float](tkes_p)
    z = Ptr[float](z_p)
    ri = Ptr[float](ri_p)
    bflxs = Ptr[float](bflxs_p)
    bprod = Ptr[float](bprod_p)
    sprod = Ptr[float](sprod_p)
    gh_a = Ptr[float](gh_a_p)
    sh_a = Ptr[float](sh_a_p)
    sm_a = Ptr[float](sm_a_p)
    ri_a = Ptr[float](ri_a_p)
    sm_aw = Ptr[float](sm_aw_p)

    i = i_col
    surf_idx = _idx2(i, pver + 1, pcols)

    bprod[surf_idx] = bflxs[i - 1]

    gg = 0.5 * vk * z[_idx2(i, pver, pcols)] * bprod[surf_idx] / (tkes[i - 1] ** (3.0 / 2.0))
    if abs(alph5 - gg * alph3) <= 1.0e-7:
        if bprod[surf_idx] > 0.0:
            gh = -3.5334
        else:
            gh = ghmin
    else:
        gh = gg / (alph5 - gg * alph3)

    if bprod[surf_idx] > 0.0:
        gh = min(max(gh, -3.5334), 0.0233)
    else:
        gh = min(max(gh, ghmin), 0.0233)

    gh_a[surf_idx] = gh
    sh_a[surf_idx] = max(0.0, alph5 / (1.0 + alph3 * gh))
    if bprod[surf_idx] > 0.0:
        sm_a[surf_idx] = max(0.0, (alph1 + alph2 * gh) / (1.0 + alph3 * gh) / (1.0 + alph4 * gh))
    else:
        sm_a[surf_idx] = max(0.0, (alph1 + alph2 * gh) / (1.0 + alph3 * gh) / (1.0 + alph4exs * gh))
    sm_aw[surf_idx] = sm_a[surf_idx] / alph1
    ri_a[surf_idx] = -(sm_a[surf_idx] / sh_a[surf_idx]) * (bprod[surf_idx] / sprod[surf_idx])

    for k in range(1, pver + 1):
        idx = _idx2(i, k, pcols)
        if ri[idx] < 0.0:
            trma = alph3 * alph4 * ri[idx] + 2.0 * b1 * (alph2 - alph4 * alph5 * ri[idx])
            trmb = (alph3 + alph4) * ri[idx] + 2.0 * b1 * (-alph5 * ri[idx] + alph1)
            trmc = ri[idx]
            det = max(trmb * trmb - 4.0 * trma * trmc, 0.0)
            gh = (-trmb + sqrt(det)) / (2.0 * trma)
            gh = min(max(gh, -3.5334), 0.0233)
            gh_a[idx] = gh
            sh_a[idx] = max(0.0, alph5 / (1.0 + alph3 * gh))
            sm_a[idx] = max(0.0, (alph1 + alph2 * gh) / (1.0 + alph3 * gh) / (1.0 + alph4 * gh))
            ri_a[idx] = ri[idx]
        else:
            if ri[idx] > ricrit:
                gh_a[idx] = ghmin
                sh_a[idx] = 0.0
                sm_a[idx] = 0.0
                ri_a[idx] = ri[idx]
            else:
                trma = alph3 * alph4exs * ri[idx] + 2.0 * b1 * (alph2 - alph4exs * alph5 * ri[idx])
                trmb = (alph3 + alph4exs) * ri[idx] + 2.0 * b1 * (-alph5 * ri[idx] + alph1)
                trmc = ri[idx]
                det = max(trmb * trmb - 4.0 * trma * trmc, 0.0)
                gh = (-trmb + sqrt(det)) / (2.0 * trma)
                gh = min(max(gh, ghmin), 0.0233)
                gh_a[idx] = gh
                sh_a[idx] = max(0.0, alph5 / (1.0 + alph3 * gh))
                sm_a[idx] = max(0.0, (alph1 + alph2 * gh) / (1.0 + alph3 * gh) / (1.0 + alph4exs * gh))
                ri_a[idx] = ri[idx]


@export
def eddy_diff_restore_fields_codon(
    ncol: int,
    pcols: int,
    pver: int,
    sl_p: cobj,
    qt_p: cobj,
    u_p: cobj,
    v_p: cobj,
    slfd_p: cobj,
    qtfd_p: cobj,
    ufd_p: cobj,
    vfd_p: cobj,
):
    sl = Ptr[float](sl_p)
    qt = Ptr[float](qt_p)
    u = Ptr[float](u_p)
    v = Ptr[float](v_p)
    slfd = Ptr[float](slfd_p)
    qtfd = Ptr[float](qtfd_p)
    ufd = Ptr[float](ufd_p)
    vfd = Ptr[float](vfd_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            slfd[_idx2(i, k, pcols)] = sl[_idx2(i, k, pcols)]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            qtfd[_idx2(i, k, pcols)] = qt[_idx2(i, k, pcols)]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            ufd[_idx2(i, k, pcols)] = u[_idx2(i, k, pcols)]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            vfd[_idx2(i, k, pcols)] = v[_idx2(i, k, pcols)]


@export
def eddy_diff_zero_nonlocal_codon(
    ncol: int,
    pcols: int,
    pver: int,
    cgh_p: cobj,
    cgs_p: cobj,
):
    cgh = Ptr[float](cgh_p)
    cgs = Ptr[float](cgs_p)

    for k in range(1, pver + 2):
        for i in range(1, ncol + 1):
            cgh[_idx2(i, k, pcols)] = 0.0

    for k in range(1, pver + 2):
        for i in range(1, ncol + 1):
            cgs[_idx2(i, k, pcols)] = 0.0


@export
def eddy_diff_init_fields_codon(
    ncol: int,
    pcols: int,
    pver: int,
    u_p: cobj,
    v_p: cobj,
    t_p: cobj,
    qv_p: cobj,
    ql_p: cobj,
    zero_p: cobj,
    zero2d_p: cobj,
    ufd_p: cobj,
    vfd_p: cobj,
    tfd_p: cobj,
    qvfd_p: cobj,
    qlfd_p: cobj,
):
    u = Ptr[float](u_p)
    v = Ptr[float](v_p)
    t = Ptr[float](t_p)
    qv = Ptr[float](qv_p)
    ql = Ptr[float](ql_p)
    zero = Ptr[float](zero_p)
    zero2d = Ptr[float](zero2d_p)
    ufd = Ptr[float](ufd_p)
    vfd = Ptr[float](vfd_p)
    tfd = Ptr[float](tfd_p)
    qvfd = Ptr[float](qvfd_p)
    qlfd = Ptr[float](qlfd_p)

    for i in range(1, pcols + 1):
        zero[i - 1] = 0.0

    for k in range(1, pver + 2):
        for i in range(1, pcols + 1):
            zero2d[_idx2(i, k, pcols)] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            ufd[_idx2(i, k, pcols)] = u[_idx2(i, k, pcols)]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            vfd[_idx2(i, k, pcols)] = v[_idx2(i, k, pcols)]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            tfd[_idx2(i, k, pcols)] = t[_idx2(i, k, pcols)]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            qvfd[_idx2(i, k, pcols)] = qv[_idx2(i, k, pcols)]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            qlfd[_idx2(i, k, pcols)] = ql[_idx2(i, k, pcols)]


@export
def eddy_diff_rebuild_thermo_codon(
    ncol: int,
    pcols: int,
    pver: int,
    cpair: float,
    latvap: float,
    latsub: float,
    gravit: float,
    rair: float,
    slfd_p: cobj,
    qtfd_p: cobj,
    qi_p: cobj,
    z_p: cobj,
    pmid_p: cobj,
    qlfd_p: cobj,
    qvfd_p: cobj,
    tfd_p: cobj,
):
    slfd = Ptr[float](slfd_p)
    qtfd = Ptr[float](qtfd_p)
    qi = Ptr[float](qi_p)
    z = Ptr[float](z_p)
    pmid = Ptr[float](pmid_p)
    qlfd = Ptr[float](qlfd_p)
    qvfd = Ptr[float](qvfd_p)
    tfd = Ptr[float](tfd_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            zval = z[idx]
            qival = qi[idx]
            qtval = qtfd[idx]
            slval = slfd[idx]
            pval = pmid[idx]

            templ = (slval - gravit * zval) / cpair
            es = eddy_diff_estblf_cb(templ)
            qs = eddy_diff_svp_to_qsat_cb(es, pval)
            templ_sq = templ * templ

            temps = templ + (qtval - qs) / (
                cpair / latvap + latvap * qs / (rair * templ_sq)
            )
            es = eddy_diff_estblf_cb(temps)
            qs = eddy_diff_svp_to_qsat_cb(es, pval)

            qlval = max(qtval - qival - qs, 0.0)
            qlfd[idx] = qlval
            qvfd[idx] = max(0.0, qtval - qival - qlval)
            tfd[idx] = (slval + latvap * qlval + latsub * qival - gravit * zval) / cpair


@export
def eddy_diff_trbintd_midpoint_codon(
    ncol: int,
    pcols: int,
    pver: int,
    cpair: float,
    latvap: float,
    latsub: float,
    gravit: float,
    zvir: float,
    t_p: cobj,
    z_p: cobj,
    qv_p: cobj,
    ql_p: cobj,
    qi_p: cobj,
    gam_p: cobj,
    qt_p: cobj,
    sl_p: cobj,
    slv_p: cobj,
    chu_p: cobj,
    chs_p: cobj,
    cmu_p: cobj,
    cms_p: cobj,
):
    t = Ptr[float](t_p)
    z = Ptr[float](z_p)
    qv = Ptr[float](qv_p)
    ql = Ptr[float](ql_p)
    qi = Ptr[float](qi_p)
    gam = Ptr[float](gam_p)
    qt = Ptr[float](qt_p)
    sl = Ptr[float](sl_p)
    slv = Ptr[float](slv_p)
    chu = Ptr[float](chu_p)
    chs = Ptr[float](chs_p)
    cmu = Ptr[float](cmu_p)
    cms = Ptr[float](cms_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            tval = t[idx]
            qvval = qv[idx]
            qlval = ql[idx]
            qival = qi[idx]

            qtval = qvval + qlval + qival
            qt[idx] = qtval

            slval = cpair * tval + gravit * z[idx] - latvap * qlval - latsub * qival
            sl[idx] = slval
            slv[idx] = slval * (1.0 + zvir * qtval)

            bfact = gravit / (tval * (1.0 + zvir * qvval - qlval - qival))
            chu[idx] = (1.0 + zvir * qtval) * bfact / cpair
            chs[idx] = (
                (1.0 + (1.0 + zvir) * gam[idx] * cpair * tval / latvap)
                / (1.0 + gam[idx])
            ) * bfact / cpair
            cmu[idx] = zvir * bfact * tval
            cms[idx] = latvap * chs[idx] - bfact * tval

    for i in range(1, ncol + 1):
        idx_last = _idx2(i, pver, pcols)
        idx_sfc = _idx2(i, pver + 1, pcols)
        chu[idx_sfc] = chu[idx_last]
        chs[idx_sfc] = chs[idx_last]
        cmu[idx_sfc] = cmu[idx_last]
        cms[idx_sfc] = cms[idx_last]


@export
def eddy_diff_trbintd_core_codon(
    ncol: int,
    pcols: int,
    pver: int,
    cpair: float,
    latvap: float,
    latsub: float,
    gravit: float,
    zvir: float,
    ntzero: float,
    t_p: cobj,
    z_p: cobj,
    u_p: cobj,
    v_p: cobj,
    qv_p: cobj,
    ql_p: cobj,
    qi_p: cobj,
    gam_p: cobj,
    pmid_p: cobj,
    cld_p: cobj,
    qt_p: cobj,
    sl_p: cobj,
    slv_p: cobj,
    slslope_p: cobj,
    qtslope_p: cobj,
    dsldp_b_p: cobj,
    dqtdp_b_p: cobj,
    chu_p: cobj,
    chs_p: cobj,
    cmu_p: cobj,
    cms_p: cobj,
    sfi_p: cobj,
    sfuh_p: cobj,
    sflh_p: cobj,
    n2_p: cobj,
    s2_p: cobj,
    ri_p: cobj,
):
    t = Ptr[float](t_p)
    z = Ptr[float](z_p)
    u = Ptr[float](u_p)
    v = Ptr[float](v_p)
    qv = Ptr[float](qv_p)
    ql = Ptr[float](ql_p)
    qi = Ptr[float](qi_p)
    gam = Ptr[float](gam_p)
    pmid = Ptr[float](pmid_p)
    cld = Ptr[float](cld_p)
    qt = Ptr[float](qt_p)
    sl = Ptr[float](sl_p)
    slv = Ptr[float](slv_p)
    slslope = Ptr[float](slslope_p)
    qtslope = Ptr[float](qtslope_p)
    dsldp_b = Ptr[float](dsldp_b_p)
    dqtdp_b = Ptr[float](dqtdp_b_p)
    chu = Ptr[float](chu_p)
    chs = Ptr[float](chs_p)
    cmu = Ptr[float](cmu_p)
    cms = Ptr[float](cms_p)
    sfi = Ptr[float](sfi_p)
    sfuh = Ptr[float](sfuh_p)
    sflh = Ptr[float](sflh_p)
    n2 = Ptr[float](n2_p)
    s2 = Ptr[float](s2_p)
    ri = Ptr[float](ri_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            tval = t[idx]
            qvval = qv[idx]
            qlval = ql[idx]
            qival = qi[idx]

            qtval = qvval + qlval + qival
            qt[idx] = qtval

            slval = cpair * tval + gravit * z[idx] - latvap * qlval - latsub * qival
            sl[idx] = slval
            slv[idx] = slval * (1.0 + zvir * qtval)

            bfact = gravit / (tval * (1.0 + zvir * qvval - qlval - qival))
            chu[idx] = (1.0 + zvir * qtval) * bfact / cpair
            chs[idx] = (
                (1.0 + (1.0 + zvir) * gam[idx] * cpair * tval / latvap)
                / (1.0 + gam[idx])
            ) * bfact / cpair
            cmu[idx] = zvir * bfact * tval
            cms[idx] = latvap * chs[idx] - bfact * tval

    for i in range(1, ncol + 1):
        idx_last = _idx2(i, pver, pcols)
        idx_sfc = _idx2(i, pver + 1, pcols)
        chu[idx_sfc] = chu[idx_last]
        chs[idx_sfc] = chs[idx_last]
        cmu[idx_sfc] = cmu[idx_last]
        cms[idx_sfc] = cms[idx_last]

    for i in range(1, ncol + 1):
        idx_pver = _idx2(i, pver, pcols)
        idx_pverm1 = _idx2(i, pver - 1, pcols)
        idx_1 = _idx2(i, 1, pcols)
        idx_2 = _idx2(i, 2, pcols)
        idx_i = i - 1

        slslope[idx_pver] = (sl[idx_pver] - sl[idx_pverm1]) / (pmid[idx_pver] - pmid[idx_pverm1])
        qtslope[idx_pver] = (qt[idx_pver] - qt[idx_pverm1]) / (pmid[idx_pver] - pmid[idx_pverm1])
        slslope[idx_1] = (sl[idx_2] - sl[idx_1]) / (pmid[idx_2] - pmid[idx_1])
        qtslope[idx_1] = (qt[idx_2] - qt[idx_1]) / (pmid[idx_2] - pmid[idx_1])
        dsldp_b[idx_i] = slslope[idx_1]
        dqtdp_b[idx_i] = qtslope[idx_1]

    for k in range(2, pver):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            idx_kp1 = _idx2(i, k + 1, pcols)
            idx_i = i - 1

            dsldp_a = dsldp_b[idx_i]
            dqtdp_a = dqtdp_b[idx_i]

            dsldp_b[idx_i] = (sl[idx_kp1] - sl[idx]) / (pmid[idx_kp1] - pmid[idx])
            dqtdp_b[idx_i] = (qt[idx_kp1] - qt[idx]) / (pmid[idx_kp1] - pmid[idx])

            product = dsldp_a * dsldp_b[idx_i]
            if product <= 0.0:
                slslope[idx] = 0.0
            elif dsldp_a < 0.0:
                slslope[idx] = max(dsldp_a, dsldp_b[idx_i])
            else:
                slslope[idx] = min(dsldp_a, dsldp_b[idx_i])

            product = dqtdp_a * dqtdp_b[idx_i]
            if product <= 0.0:
                qtslope[idx] = 0.0
            elif dqtdp_a < 0.0:
                qtslope[idx] = max(dqtdp_a, dqtdp_b[idx_i])
            else:
                qtslope[idx] = min(dqtdp_a, dqtdp_b[idx_i])

    for k in range(1, pver + 2):
        for i in range(1, ncol + 1):
            sfi[_idx2(i, k, pcols)] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            sfuh[_idx2(i, k, pcols)] = 0.0
            sflh[_idx2(i, k, pcols)] = 0.0

    for k in range(2, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            idx_km1 = _idx2(i, k - 1, pcols)
            cldval = cld[idx]
            sfuh[idx] = cldval
            sflh[idx] = cldval
            sfi[idx] = 0.5 * (sflh[idx_km1] + min(sfuh[idx], sflh[idx_km1]))

    for i in range(1, ncol + 1):
        sfi[_idx2(i, pver + 1, pcols)] = sflh[_idx2(i, pver, pcols)]

    for k in range(pver, 1, -1):
        km1 = k - 1
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            idx_km1 = _idx2(i, km1, pcols)

            rdz = 1.0 / (z[idx_km1] - z[idx])
            dsldz = (sl[idx_km1] - sl[idx]) * rdz
            dqtdz = (qt[idx_km1] - qt[idx]) * rdz

            chu[idx] = (chu[idx_km1] + chu[idx]) * 0.5
            chs[idx] = (chs[idx_km1] + chs[idx]) * 0.5
            cmu[idx] = (cmu[idx_km1] + cmu[idx]) * 0.5
            cms[idx] = (cms[idx_km1] + cms[idx]) * 0.5

            sfival = sfi[idx]
            ch = chu[idx] * (1.0 - sfival) + chs[idx] * sfival
            cm = cmu[idx] * (1.0 - sfival) + cms[idx] * sfival

            n2[idx] = ch * dsldz + cm * dqtdz

            du = u[idx_km1] - u[idx]
            dv = v[idx_km1] - v[idx]
            rdz_sq = rdz * rdz
            s2val = (du * du + dv * dv) * rdz_sq
            s2[idx] = max(ntzero, s2val)
            ri[idx] = n2[idx] / s2[idx]

    for i in range(1, ncol + 1):
        idx_1 = _idx2(i, 1, pcols)
        idx_2 = _idx2(i, 2, pcols)
        n2[idx_1] = n2[idx_2]
        s2[idx_1] = s2[idx_2]
        ri[idx_1] = ri[idx_2]


@export
def eddy_diff_trbintd_slopes_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pmid_p: cobj,
    sl_p: cobj,
    qt_p: cobj,
    slslope_p: cobj,
    qtslope_p: cobj,
    dsldp_b_p: cobj,
    dqtdp_b_p: cobj,
):
    pmid = Ptr[float](pmid_p)
    sl = Ptr[float](sl_p)
    qt = Ptr[float](qt_p)
    slslope = Ptr[float](slslope_p)
    qtslope = Ptr[float](qtslope_p)
    dsldp_b = Ptr[float](dsldp_b_p)
    dqtdp_b = Ptr[float](dqtdp_b_p)

    for i in range(1, ncol + 1):
        idx_pver = _idx2(i, pver, pcols)
        idx_pverm1 = _idx2(i, pver - 1, pcols)
        idx_1 = _idx2(i, 1, pcols)
        idx_2 = _idx2(i, 2, pcols)
        idx_i = i - 1

        slslope[idx_pver] = (sl[idx_pver] - sl[idx_pverm1]) / (pmid[idx_pver] - pmid[idx_pverm1])
        qtslope[idx_pver] = (qt[idx_pver] - qt[idx_pverm1]) / (pmid[idx_pver] - pmid[idx_pverm1])
        slslope[idx_1] = (sl[idx_2] - sl[idx_1]) / (pmid[idx_2] - pmid[idx_1])
        qtslope[idx_1] = (qt[idx_2] - qt[idx_1]) / (pmid[idx_2] - pmid[idx_1])
        dsldp_b[idx_i] = slslope[idx_1]
        dqtdp_b[idx_i] = qtslope[idx_1]

    for k in range(2, pver):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            idx_kp1 = _idx2(i, k + 1, pcols)
            idx_i = i - 1

            dsldp_a = dsldp_b[idx_i]
            dqtdp_a = dqtdp_b[idx_i]

            dsldp_b[idx_i] = (sl[idx_kp1] - sl[idx]) / (pmid[idx_kp1] - pmid[idx])
            dqtdp_b[idx_i] = (qt[idx_kp1] - qt[idx]) / (pmid[idx_kp1] - pmid[idx])

            product = dsldp_a * dsldp_b[idx_i]
            if product <= 0.0:
                slslope[idx] = 0.0
            elif dsldp_a < 0.0:
                slslope[idx] = max(dsldp_a, dsldp_b[idx_i])
            else:
                slslope[idx] = min(dsldp_a, dsldp_b[idx_i])

            product = dqtdp_a * dqtdp_b[idx_i]
            if product <= 0.0:
                qtslope[idx] = 0.0
            elif dqtdp_a < 0.0:
                qtslope[idx] = max(dqtdp_a, dqtdp_b[idx_i])
            else:
                qtslope[idx] = min(dqtdp_a, dqtdp_b[idx_i])


@export
def eddy_diff_trbintd_sfdiag_interface_codon(
    ncol: int,
    pcols: int,
    pver: int,
    ntzero: float,
    cld_p: cobj,
    u_p: cobj,
    v_p: cobj,
    z_p: cobj,
    sl_p: cobj,
    qt_p: cobj,
    chu_p: cobj,
    chs_p: cobj,
    cmu_p: cobj,
    cms_p: cobj,
    sfi_p: cobj,
    sfuh_p: cobj,
    sflh_p: cobj,
    n2_p: cobj,
    s2_p: cobj,
    ri_p: cobj,
):
    cld = Ptr[float](cld_p)
    u = Ptr[float](u_p)
    v = Ptr[float](v_p)
    z = Ptr[float](z_p)
    sl = Ptr[float](sl_p)
    qt = Ptr[float](qt_p)
    chu = Ptr[float](chu_p)
    chs = Ptr[float](chs_p)
    cmu = Ptr[float](cmu_p)
    cms = Ptr[float](cms_p)
    sfi = Ptr[float](sfi_p)
    sfuh = Ptr[float](sfuh_p)
    sflh = Ptr[float](sflh_p)
    n2 = Ptr[float](n2_p)
    s2 = Ptr[float](s2_p)
    ri = Ptr[float](ri_p)

    for k in range(1, pver + 2):
        for i in range(1, ncol + 1):
            sfi[_idx2(i, k, pcols)] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            sfuh[_idx2(i, k, pcols)] = 0.0
            sflh[_idx2(i, k, pcols)] = 0.0

    for k in range(2, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            idx_km1 = _idx2(i, k - 1, pcols)
            cldval = cld[idx]
            sfuh[idx] = cldval
            sflh[idx] = cldval
            sfi[idx] = 0.5 * (sflh[idx_km1] + min(sfuh[idx], sflh[idx_km1]))

    for i in range(1, ncol + 1):
        sfi[_idx2(i, pver + 1, pcols)] = sflh[_idx2(i, pver, pcols)]

    for k in range(pver, 1, -1):
        km1 = k - 1
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            idx_km1 = _idx2(i, km1, pcols)

            rdz = 1.0 / (z[idx_km1] - z[idx])
            dsldz = (sl[idx_km1] - sl[idx]) * rdz
            dqtdz = (qt[idx_km1] - qt[idx]) * rdz

            chu[idx] = (chu[idx_km1] + chu[idx]) * 0.5
            chs[idx] = (chs[idx_km1] + chs[idx]) * 0.5
            cmu[idx] = (cmu[idx_km1] + cmu[idx]) * 0.5
            cms[idx] = (cms[idx_km1] + cms[idx]) * 0.5

            sfival = sfi[idx]
            ch = chu[idx] * (1.0 - sfival) + chs[idx] * sfival
            cm = cmu[idx] * (1.0 - sfival) + cms[idx] * sfival

            n2[idx] = ch * dsldz + cm * dqtdz

            du = u[idx_km1] - u[idx]
            dv = v[idx_km1] - v[idx]
            rdz_sq = rdz * rdz
            s2val = (du * du + dv * dv) * rdz_sq
            s2[idx] = max(ntzero, s2val)
            ri[idx] = n2[idx] / s2[idx]

    for i in range(1, ncol + 1):
        idx_1 = _idx2(i, 1, pcols)
        idx_2 = _idx2(i, 2, pcols)
        n2[idx_1] = n2[idx_2]
        s2[idx_1] = s2[idx_2]
        ri[idx_1] = ri[idx_2]


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
