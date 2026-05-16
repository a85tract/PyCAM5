from math import exp, log

@export
def phys_timestep_init_select_branches_codon(
    cam3_aero_on: int,
    cam3_ozone_on: int,
    do_waccm_ions: int,
    branch_mask_p: cobj,
):
    branch_mask = Ptr[int](branch_mask_p)

    mask = 0
    if cam3_aero_on != 0:
        mask |= 1
    if cam3_ozone_on != 0:
        mask |= 2
    if do_waccm_ions != 0:
        mask |= 4

    branch_mask[0] = mask


@inline
def _idx(i: int) -> int:
    return i - 1


@inline
def _field2_idx(i: int, k: int, ld1: int) -> int:
    """state_t/tini/tend_dtdt/dtcore/tmp_t declared as (ld1, pver)"""
    return (i - 1) + (k - 1) * ld1


@inline
def _field3_idx(i: int, k: int, m: int, ld1: int, ld2: int) -> int:
    """state_q declared as (ld1, ld2, pcnst)"""
    return (i - 1) + (k - 1) * ld1 + (m - 1) * ld1 * ld2


@export
def tropopause_output_prep_codon(
    ncol: int,
    pcols: int,
    pver: int,
    notfound: int,
    fillvalue: float,
    trop_lev_p: cobj,
    trop_z_p: cobj,
    state_zm_p: cobj,
    trop_pdf_p: cobj,
    trop_found_p: cobj,
    trop_dz_p: cobj,
):
    trop_lev = Ptr[int](trop_lev_p)
    trop_z = Ptr[float](trop_z_p)
    state_zm = Ptr[float](state_zm_p)
    trop_pdf = Ptr[float](trop_pdf_p)
    trop_found = Ptr[float](trop_found_p)
    trop_dz = Ptr[float](trop_dz_p)

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            trop_pdf[_field2_idx(i, k, pcols)] = 0.0
            trop_dz[_field2_idx(i, k, pcols)] = fillvalue

    for i in range(1, pcols + 1):
        trop_found[_idx(i)] = 0.0

    for i in range(1, ncol + 1):
        lev = trop_lev[_idx(i)]
        if lev != notfound:
            trop_pdf[_field2_idx(i, lev, pcols)] = 1.0
            trop_found[_idx(i)] = 1.0
            for k in range(1, pver + 1):
                trop_dz[_field2_idx(i, k, pcols)] = (
                    state_zm[_field2_idx(i, k, pcols)] - trop_z[_idx(i)]
                )


@inline
def _tropopause_interp_t(
    i: int,
    lev: int,
    trop_p: float,
    pcols: int,
    pver: int,
    state_t: Ptr[float],
    state_pmid: Ptr[float],
) -> float:
    trop_t = 0.0
    if trop_p == state_pmid[_field2_idx(i, lev, pcols)]:
        trop_t = state_t[_field2_idx(i, lev, pcols)]
    elif trop_p < state_pmid[_field2_idx(i, lev, pcols)]:
        if lev > 1:
            dtdlogp = (
                state_t[_field2_idx(i, lev, pcols)]
                - state_t[_field2_idx(i, lev - 1, pcols)]
            ) / (
                log(state_pmid[_field2_idx(i, lev, pcols)])
                - log(state_pmid[_field2_idx(i, lev - 1, pcols)])
            )
            trop_t = state_t[_field2_idx(i, lev, pcols)] + (
                log(trop_p) - log(state_pmid[_field2_idx(i, lev, pcols)])
            ) * dtdlogp
    else:
        if lev < pver:
            dtdlogp = (
                state_t[_field2_idx(i, lev + 1, pcols)]
                - state_t[_field2_idx(i, lev, pcols)]
            ) / (
                log(state_pmid[_field2_idx(i, lev + 1, pcols)])
                - log(state_pmid[_field2_idx(i, lev, pcols)])
            )
            trop_t = state_t[_field2_idx(i, lev, pcols)] + (
                log(trop_p) - log(state_pmid[_field2_idx(i, lev, pcols)])
            ) * dtdlogp
    return trop_t


@inline
def _tropopause_interp_z(
    i: int,
    lev: int,
    trop_p: float,
    pcols: int,
    state_zm: Ptr[float],
    state_zi: Ptr[float],
    state_pmid: Ptr[float],
    state_pint: Ptr[float],
) -> float:
    trop_z = 0.0
    if trop_p == state_pmid[_field2_idx(i, lev, pcols)]:
        trop_z = state_zm[_field2_idx(i, lev, pcols)]
    elif trop_p < state_pmid[_field2_idx(i, lev, pcols)]:
        dzdlogp = (
            state_zm[_field2_idx(i, lev, pcols)] - state_zi[_field2_idx(i, lev, pcols)]
        ) / (
            log(state_pmid[_field2_idx(i, lev, pcols)])
            - log(state_pint[_field2_idx(i, lev, pcols)])
        )
        trop_z = state_zm[_field2_idx(i, lev, pcols)] + (
            log(trop_p) - log(state_pmid[_field2_idx(i, lev, pcols)])
        ) * dzdlogp
    else:
        dzdlogp = (
            state_zm[_field2_idx(i, lev, pcols)] - state_zi[_field2_idx(i, lev + 1, pcols)]
        ) / (
            log(state_pmid[_field2_idx(i, lev, pcols)])
            - log(state_pint[_field2_idx(i, lev + 1, pcols)])
        )
        trop_z = state_zm[_field2_idx(i, lev, pcols)] + (
            log(trop_p) - log(state_pmid[_field2_idx(i, lev, pcols)])
        ) * dzdlogp
    return trop_z


@inline
def _tropopause_twmo_pressure(
    i: int,
    pcols: int,
    pver: int,
    cnst_kap: float,
    cnst_ka1: float,
    cnst_faktor: float,
    state_t: Ptr[float],
    state_pmid: Ptr[float],
) -> float:
    gam = -0.002
    plimu = 45000.0
    pliml = 7500.0
    deltaz = 2000.0

    trp = -99.0
    level = pver
    pmk = 0.5 * (
        state_pmid[_field2_idx(i, level - 1, pcols)] ** cnst_kap
        + state_pmid[_field2_idx(i, level, pcols)] ** cnst_kap
    )
    pm = pmk ** (1 / cnst_kap)
    a = (
        state_t[_field2_idx(i, level - 1, pcols)]
        - state_t[_field2_idx(i, level, pcols)]
    ) / (
        state_pmid[_field2_idx(i, level - 1, pcols)] ** cnst_kap
        - state_pmid[_field2_idx(i, level, pcols)] ** cnst_kap
    )
    b = state_t[_field2_idx(i, level, pcols)] - (
        a * state_pmid[_field2_idx(i, level, pcols)] ** cnst_kap
    )
    tm = a * pmk + b
    dtdp = a * cnst_kap * (pm ** cnst_ka1)
    dtdz = cnst_faktor * dtdp * pm / tm

    for j in range(level - 1, 1, -1):
        pm0 = pm
        pmk0 = pmk
        dtdz0 = dtdz

        pmk = 0.5 * (
            state_pmid[_field2_idx(i, j - 1, pcols)] ** cnst_kap
            + state_pmid[_field2_idx(i, j, pcols)] ** cnst_kap
        )
        pm = pmk ** (1 / cnst_kap)
        a = (
            state_t[_field2_idx(i, j - 1, pcols)]
            - state_t[_field2_idx(i, j, pcols)]
        ) / (
            state_pmid[_field2_idx(i, j - 1, pcols)] ** cnst_kap
            - state_pmid[_field2_idx(i, j, pcols)] ** cnst_kap
        )
        b = state_t[_field2_idx(i, j, pcols)] - (
            a * state_pmid[_field2_idx(i, j, pcols)] ** cnst_kap
        )
        tm = a * pmk + b
        dtdp = a * cnst_kap * (pm ** cnst_ka1)
        dtdz = cnst_faktor * dtdp * pm / tm

        if dtdz <= gam:
            continue
        if pm > plimu:
            continue

        if dtdz0 < gam:
            ag = (dtdz - dtdz0) / (pmk - pmk0)
            bg = dtdz0 - (ag * pmk0)
            ptph = exp(log((gam - bg) / ag) / cnst_kap)
        else:
            ptph = pm

        if ptph < pliml:
            continue
        if ptph > plimu:
            continue

        p2km = ptph + deltaz * (pm / tm) * cnst_faktor
        asum = 0.0
        icount = 0
        valid = True

        for jj in range(j, 1, -1):
            pmk2 = 0.5 * (
                state_pmid[_field2_idx(i, jj - 1, pcols)] ** cnst_kap
                + state_pmid[_field2_idx(i, jj, pcols)] ** cnst_kap
            )
            pm2 = pmk2 ** (1 / cnst_kap)
            if pm2 > ptph:
                continue
            if pm2 < p2km:
                break

            a2 = state_t[_field2_idx(i, jj - 1, pcols)] - state_t[_field2_idx(i, jj, pcols)]
            a2 = a2 / (
                state_pmid[_field2_idx(i, jj - 1, pcols)] ** cnst_kap
                - state_pmid[_field2_idx(i, jj, pcols)] ** cnst_kap
            )
            b2 = state_t[_field2_idx(i, jj, pcols)] - (
                a2 * state_pmid[_field2_idx(i, jj, pcols)] ** cnst_kap
            )
            tm2 = a2 * pmk2 + b2
            dtdp2 = a2 * cnst_kap * (pm2 ** (cnst_kap - 1))
            dtdz2 = cnst_faktor * dtdp2 * pm2 / tm2
            asum = asum + dtdz2
            icount = icount + 1
            aquer = asum / float(icount)
            if aquer <= gam:
                valid = False
                break

        if not valid:
            continue

        trp = ptph
        break

    return trp


@export
def tropopause_twmo_codon(
    ncol: int,
    pcols: int,
    pver: int,
    notfound: int,
    write_tropp: int,
    write_tropt: int,
    write_tropz: int,
    cnst_kap: float,
    cnst_ka1: float,
    cnst_faktor: float,
    state_t_p: cobj,
    state_pmid_p: cobj,
    state_pint_p: cobj,
    state_zm_p: cobj,
    state_zi_p: cobj,
    trop_lev_p: cobj,
    trop_p_p: cobj,
    trop_t_p: cobj,
    trop_z_p: cobj,
):
    state_t = Ptr[float](state_t_p)
    state_pmid = Ptr[float](state_pmid_p)
    state_pint = Ptr[float](state_pint_p)
    state_zm = Ptr[float](state_zm_p)
    state_zi = Ptr[float](state_zi_p)
    trop_lev = Ptr[int](trop_lev_p)
    trop_p = Ptr[float](trop_p_p)
    trop_t = Ptr[float](trop_t_p)
    trop_z = Ptr[float](trop_z_p)

    for i in range(1, ncol + 1):
        if trop_lev[_idx(i)] == notfound:
            tp = _tropopause_twmo_pressure(
                i, pcols, pver, cnst_kap, cnst_ka1, cnst_faktor, state_t, state_pmid
            )
            if tp > 0.0:
                for k in range(pver, 1, -1):
                    if tp >= state_pint[_field2_idx(i, k, pcols)]:
                        trop_lev[_idx(i)] = k
                        break

                lev = trop_lev[_idx(i)]
                if write_tropp != 0:
                    trop_p[_idx(i)] = tp
                if write_tropt != 0:
                    trop_t[_idx(i)] = _tropopause_interp_t(i, lev, tp, pcols, pver, state_t, state_pmid)
                if write_tropz != 0:
                    trop_z[_idx(i)] = _tropopause_interp_z(
                        i, lev, tp, pcols, state_zm, state_zi, state_pmid, state_pint
                    )


@export
def tphysbc_precip_ops_codon(
    mode: int,
    ncol: int,
    pcols: int,
    cld_macmic_num_steps: int,
    prec_sed_macmic_p: cobj,
    snow_sed_macmic_p: cobj,
    prec_pcw_macmic_p: cobj,
    snow_pcw_macmic_p: cobj,
    prec_sed_p: cobj,
    snow_sed_p: cobj,
    prec_pcw_p: cobj,
    snow_pcw_p: cobj,
    prec_str_p: cobj,
    snow_str_p: cobj,
    prec_sed_carma_p: cobj,
    snow_sed_carma_p: cobj,
):
    prec_sed_macmic = Ptr[float](prec_sed_macmic_p)
    snow_sed_macmic = Ptr[float](snow_sed_macmic_p)
    prec_pcw_macmic = Ptr[float](prec_pcw_macmic_p)
    snow_pcw_macmic = Ptr[float](snow_pcw_macmic_p)
    prec_sed = Ptr[float](prec_sed_p)
    snow_sed = Ptr[float](snow_sed_p)
    prec_pcw = Ptr[float](prec_pcw_p)
    snow_pcw = Ptr[float](snow_pcw_p)
    prec_str = Ptr[float](prec_str_p)
    snow_str = Ptr[float](snow_str_p)
    prec_sed_carma = Ptr[float](prec_sed_carma_p)
    snow_sed_carma = Ptr[float](snow_sed_carma_p)

    if mode == 0:
        for i in range(1, pcols + 1):
            prec_sed_macmic[_idx(i)] = 0.0
            snow_sed_macmic[_idx(i)] = 0.0
            prec_pcw_macmic[_idx(i)] = 0.0
            snow_pcw_macmic[_idx(i)] = 0.0
    elif mode == 1:
        for i in range(1, ncol + 1):
            prec_sed_macmic[_idx(i)] = prec_sed_macmic[_idx(i)] + prec_sed[_idx(i)]
            snow_sed_macmic[_idx(i)] = snow_sed_macmic[_idx(i)] + snow_sed[_idx(i)]
            prec_pcw_macmic[_idx(i)] = prec_pcw_macmic[_idx(i)] + prec_pcw[_idx(i)]
            snow_pcw_macmic[_idx(i)] = snow_pcw_macmic[_idx(i)] + snow_pcw[_idx(i)]
    elif mode == 2:
        for i in range(1, ncol + 1):
            prec_sed[_idx(i)] = prec_sed_macmic[_idx(i)] / cld_macmic_num_steps
            snow_sed[_idx(i)] = snow_sed_macmic[_idx(i)] / cld_macmic_num_steps
            prec_pcw[_idx(i)] = prec_pcw_macmic[_idx(i)] / cld_macmic_num_steps
            snow_pcw[_idx(i)] = snow_pcw_macmic[_idx(i)] / cld_macmic_num_steps
            prec_str[_idx(i)] = prec_pcw[_idx(i)] + prec_sed[_idx(i)]
            snow_str[_idx(i)] = snow_pcw[_idx(i)] + snow_sed[_idx(i)]
    elif mode == 3:
        for i in range(1, ncol + 1):
            prec_sed[_idx(i)] = prec_sed[_idx(i)] + prec_sed_carma[_idx(i)]
            snow_sed[_idx(i)] = snow_sed[_idx(i)] + snow_sed_carma[_idx(i)]


@export
def tphysac_flx_net_update_codon(
    ncol: int,
    pcols: int,
    tend_flx_net_p: cobj,
    cam_in_shf_p: cobj,
    cam_out_precc_p: cobj,
    cam_out_precl_p: cobj,
    cam_out_precsc_p: cobj,
    cam_out_precsl_p: cobj,
    latvap_local: float,
    latice_local: float,
    rhoh2o_local: float,
):
    tend_flx_net = Ptr[float](tend_flx_net_p)
    cam_in_shf = Ptr[float](cam_in_shf_p)
    cam_out_precc = Ptr[float](cam_out_precc_p)
    cam_out_precl = Ptr[float](cam_out_precl_p)
    cam_out_precsc = Ptr[float](cam_out_precsc_p)
    cam_out_precsl = Ptr[float](cam_out_precsl_p)

    for i in range(1, ncol + 1):
        tend_flx_net[_idx(i)] = (
            tend_flx_net[_idx(i)]
            + cam_in_shf[_idx(i)]
            + (cam_out_precc[_idx(i)] + cam_out_precl[_idx(i)]) * latvap_local * rhoh2o_local
            + (cam_out_precsc[_idx(i)] + cam_out_precsl[_idx(i)]) * latice_local * rhoh2o_local
        )


@export
def tphysac_t_update_codon(
    ncol: int,
    pcols: int,
    pver: int,
    ztodt: float,
    state_t_p: cobj,
    tini_p: cobj,
    tend_dtdt_p: cobj,
    dtcore_p: cobj,
    tmp_t_p: cobj,
):
    state_t = Ptr[float](state_t_p)
    tini = Ptr[float](tini_p)
    tend_dtdt = Ptr[float](tend_dtdt_p)
    dtcore = Ptr[float](dtcore_p)
    tmp_t = Ptr[float](tmp_t_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            tmp_t[_field2_idx(i, k, pcols)] = state_t[_field2_idx(i, k, pcols)]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            state_t[_field2_idx(i, k, pcols)] = tini[_field2_idx(i, k, pcols)] + ztodt * tend_dtdt[
                _field2_idx(i, k, pcols)
            ]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            dtcore[_field2_idx(i, k, pcols)] = state_t[_field2_idx(i, k, pcols)]


@export
def tphysac_q_snapshot_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    ixcldliq: int,
    ixcldice: int,
    state_q_p: cobj,
    tmp_q_p: cobj,
    tmp_cldliq_p: cobj,
    tmp_cldice_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    tmp_q = Ptr[float](tmp_q_p)
    tmp_cldliq = Ptr[float](tmp_cldliq_p)
    tmp_cldice = Ptr[float](tmp_cldice_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            tmp_q[_field2_idx(i, k, pcols)] = state_q[_field3_idx(i, k, 1, pcols, pver)]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            tmp_cldliq[_field2_idx(i, k, pcols)] = state_q[_field3_idx(i, k, ixcldliq, pcols, pver)]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            tmp_cldice[_field2_idx(i, k, pcols)] = state_q[_field3_idx(i, k, ixcldice, pcols, pver)]


@export
def tphysbc_qini_snapshot_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    ixcldliq: int,
    ixcldice: int,
    state_q_p: cobj,
    qini_p: cobj,
    cldliqini_p: cobj,
    cldiceini_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    qini = Ptr[float](qini_p)
    cldliqini = Ptr[float](cldliqini_p)
    cldiceini = Ptr[float](cldiceini_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            qini[_field2_idx(i, k, pcols)] = state_q[_field3_idx(i, k, 1, pcols, pver)]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            cldliqini[_field2_idx(i, k, pcols)] = state_q[_field3_idx(i, k, ixcldliq, pcols, pver)]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            cldiceini[_field2_idx(i, k, pcols)] = state_q[_field3_idx(i, k, ixcldice, pcols, pver)]


@export
def tphysbc_dadadj_input_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    state_t_p: cobj,
    state_q_p: cobj,
    ptend_s_p: cobj,
    ptend_q_p: cobj,
):
    state_t = Ptr[float](state_t_p)
    state_q = Ptr[float](state_q_p)
    ptend_s = Ptr[float](ptend_s_p)
    ptend_q = Ptr[float](ptend_q_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            ptend_s[_field2_idx(i, k, pcols)] = state_t[_field2_idx(i, k, pcols)]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            ptend_q[_field3_idx(i, k, 1, pcols, pver)] = state_q[_field3_idx(i, k, 1, pcols, pver)]


@export
def tphysbc_dadadj_output_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    ztodt: float,
    cpair_local: float,
    state_t_p: cobj,
    state_q_p: cobj,
    ptend_s_p: cobj,
    ptend_q_p: cobj,
):
    state_t = Ptr[float](state_t_p)
    state_q = Ptr[float](state_q_p)
    ptend_s = Ptr[float](ptend_s_p)
    ptend_q = Ptr[float](ptend_q_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            ptend_s[_field2_idx(i, k, pcols)] = (
                (ptend_s[_field2_idx(i, k, pcols)] - state_t[_field2_idx(i, k, pcols)]) / ztodt
            ) * cpair_local

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            ptend_q[_field3_idx(i, k, 1, pcols, pver)] = (
                ptend_q[_field3_idx(i, k, 1, pcols, pver)] - state_q[_field3_idx(i, k, 1, pcols, pver)]
            ) / ztodt


@export
def tphysbc_dtcore_update_codon(
    ncol: int,
    pcols: int,
    pver: int,
    ztodt: float,
    tini_p: cobj,
    dtcore_p: cobj,
    tend_dtdt_p: cobj,
):
    tini = Ptr[float](tini_p)
    dtcore = Ptr[float](dtcore_p)
    tend_dtdt = Ptr[float](tend_dtdt_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            dtcore[_field2_idx(i, k, pcols)] = (
                (tini[_field2_idx(i, k, pcols)] - dtcore[_field2_idx(i, k, pcols)]) / ztodt
            ) + tend_dtdt[_field2_idx(i, k, pcols)]


@export
def tphysbc_init_fields_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    fracis_p: cobj,
    tend_dtdt_p: cobj,
    tend_dudt_p: cobj,
    tend_dvdt_p: cobj,
):
    fracis = Ptr[float](fracis_p)
    tend_dtdt = Ptr[float](tend_dtdt_p)
    tend_dudt = Ptr[float](tend_dudt_p)
    tend_dvdt = Ptr[float](tend_dvdt_p)

    for m in range(1, pcnst + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                fracis[_field3_idx(i, k, m, pcols, pver)] = 1.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            tend_dtdt[_field2_idx(i, k, pcols)] = 0.0
            tend_dudt[_field2_idx(i, k, pcols)] = 0.0
            tend_dvdt[_field2_idx(i, k, pcols)] = 0.0


@export
def tphysbc_tini_copy_codon(
    ncol: int,
    pcols: int,
    pver: int,
    state_t_p: cobj,
    tini_p: cobj,
):
    state_t = Ptr[float](state_t_p)
    tini = Ptr[float](tini_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            tini[_field2_idx(i, k, pcols)] = state_t[_field2_idx(i, k, pcols)]


@export
def tphysbc_flx_cnd_sum_codon(
    ncol: int,
    pcols: int,
    a_p: cobj,
    b_p: cobj,
    out_p: cobj,
):
    a = Ptr[float](a_p)
    b = Ptr[float](b_p)
    out = Ptr[float](out_p)

    for i in range(1, ncol + 1):
        out[_idx(i)] = a[_idx(i)] + b[_idx(i)]


@export
def tphysbc_macrop_fluxes_codon(
    mode: int,
    ncol: int,
    pcols: int,
    rliq_p: cobj,
    det_s_p: cobj,
    flx_cnd_p: cobj,
    flx_heat_p: cobj,
    shf_p: cobj,
):
    rliq = Ptr[float](rliq_p)
    det_s = Ptr[float](det_s_p)
    flx_cnd = Ptr[float](flx_cnd_p)
    flx_heat = Ptr[float](flx_heat_p)

    for i in range(1, ncol + 1):
        flx_cnd[_idx(i)] = -1.0 * rliq[_idx(i)]

    if mode == 1:
        for i in range(1, ncol + 1):
            flx_heat[_idx(i)] = det_s[_idx(i)]
    else:
        shf = Ptr[float](shf_p)
        for i in range(1, ncol + 1):
            flx_heat[_idx(i)] = shf[_idx(i)] + det_s[_idx(i)]


@export
def tphysbc_radheat_flx_net_codon(
    ncol: int,
    pcols: int,
    tend_flx_net_p: cobj,
    net_flx_p: cobj,
):
    tend_flx_net = Ptr[float](tend_flx_net_p)
    net_flx = Ptr[float](net_flx_p)

    for i in range(1, ncol + 1):
        tend_flx_net[_idx(i)] = net_flx[_idx(i)]


@export
def tphyspkg_flux_batch_precip_ops_codon(
    mode: int,
    ncol: int,
    pcols: int,
    cld_macmic_num_steps: int,
    prec_sed_macmic_p: cobj,
    snow_sed_macmic_p: cobj,
    prec_pcw_macmic_p: cobj,
    snow_pcw_macmic_p: cobj,
    prec_sed_p: cobj,
    snow_sed_p: cobj,
    prec_pcw_p: cobj,
    snow_pcw_p: cobj,
    prec_str_p: cobj,
    snow_str_p: cobj,
    prec_sed_carma_p: cobj,
    snow_sed_carma_p: cobj,
):
    tphysbc_precip_ops_codon(
        mode,
        ncol,
        pcols,
        cld_macmic_num_steps,
        prec_sed_macmic_p,
        snow_sed_macmic_p,
        prec_pcw_macmic_p,
        snow_pcw_macmic_p,
        prec_sed_p,
        snow_sed_p,
        prec_pcw_p,
        snow_pcw_p,
        prec_str_p,
        snow_str_p,
        prec_sed_carma_p,
        snow_sed_carma_p,
    )


@export
def tphyspkg_flux_batch_dispatch_codon(
    stage: int,
    mode: int,
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    cld_macmic_num_steps: int,
    ixcldliq: int,
    ixcldice: int,
    ztodt: float,
    latvap_local: float,
    latice_local: float,
    rhoh2o_local: float,
    prec_sed_macmic_p: cobj,
    snow_sed_macmic_p: cobj,
    prec_pcw_macmic_p: cobj,
    snow_pcw_macmic_p: cobj,
    prec_sed_p: cobj,
    snow_sed_p: cobj,
    prec_pcw_p: cobj,
    snow_pcw_p: cobj,
    prec_str_p: cobj,
    snow_str_p: cobj,
    prec_sed_carma_p: cobj,
    snow_sed_carma_p: cobj,
    tend_flx_net_p: cobj,
    cam_in_shf_p: cobj,
    cam_out_precc_p: cobj,
    cam_out_precl_p: cobj,
    cam_out_precsc_p: cobj,
    cam_out_precsl_p: cobj,
    state_t_p: cobj,
    tini_p: cobj,
    tend_dtdt_p: cobj,
    dtcore_p: cobj,
    tmp_t_p: cobj,
    state_q_p: cobj,
    tmp_q_p: cobj,
    tmp_cldliq_p: cobj,
    tmp_cldice_p: cobj,
    a_p: cobj,
    b_p: cobj,
    out_p: cobj,
    rliq_p: cobj,
    det_s_p: cobj,
    flx_cnd_p: cobj,
    flx_heat_p: cobj,
    shf_p: cobj,
    net_flx_p: cobj,
):
    if stage == 1:
        tphysbc_precip_ops_codon(
            mode,
            ncol,
            pcols,
            cld_macmic_num_steps,
            prec_sed_macmic_p,
            snow_sed_macmic_p,
            prec_pcw_macmic_p,
            snow_pcw_macmic_p,
            prec_sed_p,
            snow_sed_p,
            prec_pcw_p,
            snow_pcw_p,
            prec_str_p,
            snow_str_p,
            prec_sed_carma_p,
            snow_sed_carma_p,
        )
    elif stage == 2:
        tphysac_flx_net_update_codon(
            ncol,
            pcols,
            tend_flx_net_p,
            cam_in_shf_p,
            cam_out_precc_p,
            cam_out_precl_p,
            cam_out_precsc_p,
            cam_out_precsl_p,
            latvap_local,
            latice_local,
            rhoh2o_local,
        )
    elif stage == 3:
        tphysac_t_update_codon(ncol, pcols, pver, ztodt, state_t_p, tini_p, tend_dtdt_p, dtcore_p, tmp_t_p)
    elif stage == 4:
        tphysac_q_snapshot_codon(
            ncol, pcols, pver, pcnst, ixcldliq, ixcldice, state_q_p, tmp_q_p, tmp_cldliq_p, tmp_cldice_p
        )
    elif stage == 5:
        tphysbc_flx_cnd_sum_codon(ncol, pcols, a_p, b_p, out_p)
    elif stage == 6:
        tphysbc_macrop_fluxes_codon(mode, ncol, pcols, rliq_p, det_s_p, flx_cnd_p, flx_heat_p, shf_p)
    elif stage == 7:
        tphysbc_radheat_flx_net_codon(ncol, pcols, tend_flx_net_p, net_flx_p)


@export
def tphyspkg_flux_batch_flx_net_update_codon(
    ncol: int,
    pcols: int,
    tend_flx_net_p: cobj,
    cam_in_shf_p: cobj,
    cam_out_precc_p: cobj,
    cam_out_precl_p: cobj,
    cam_out_precsc_p: cobj,
    cam_out_precsl_p: cobj,
    latvap_local: float,
    latice_local: float,
    rhoh2o_local: float,
):
    tphysac_flx_net_update_codon(
        ncol,
        pcols,
        tend_flx_net_p,
        cam_in_shf_p,
        cam_out_precc_p,
        cam_out_precl_p,
        cam_out_precsc_p,
        cam_out_precsl_p,
        latvap_local,
        latice_local,
        rhoh2o_local,
    )


@export
def tphyspkg_flux_batch_t_update_codon(
    ncol: int,
    pcols: int,
    pver: int,
    ztodt: float,
    state_t_p: cobj,
    tini_p: cobj,
    tend_dtdt_p: cobj,
    dtcore_p: cobj,
    tmp_t_p: cobj,
):
    tphysac_t_update_codon(ncol, pcols, pver, ztodt, state_t_p, tini_p, tend_dtdt_p, dtcore_p, tmp_t_p)


@export
def tphyspkg_flux_batch_q_snapshot_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    ixcldliq: int,
    ixcldice: int,
    state_q_p: cobj,
    tmp_q_p: cobj,
    tmp_cldliq_p: cobj,
    tmp_cldice_p: cobj,
):
    tphysac_q_snapshot_codon(
        ncol, pcols, pver, pcnst, ixcldliq, ixcldice, state_q_p, tmp_q_p, tmp_cldliq_p, tmp_cldice_p
    )


@export
def tphyspkg_flux_batch_flx_cnd_sum_codon(
    ncol: int,
    pcols: int,
    a_p: cobj,
    b_p: cobj,
    out_p: cobj,
):
    tphysbc_flx_cnd_sum_codon(ncol, pcols, a_p, b_p, out_p)


@export
def tphyspkg_flux_batch_macrop_fluxes_codon(
    mode: int,
    ncol: int,
    pcols: int,
    rliq_p: cobj,
    det_s_p: cobj,
    flx_cnd_p: cobj,
    flx_heat_p: cobj,
    shf_p: cobj,
):
    tphysbc_macrop_fluxes_codon(mode, ncol, pcols, rliq_p, det_s_p, flx_cnd_p, flx_heat_p, shf_p)


@export
def tphyspkg_flux_batch_radheat_flx_net_codon(
    ncol: int,
    pcols: int,
    tend_flx_net_p: cobj,
    net_flx_p: cobj,
):
    tphysbc_radheat_flx_net_codon(ncol, pcols, tend_flx_net_p, net_flx_p)


@export
def tphysbc_zero_buffers_codon(
    pcols: int,
    pcnst: int,
    zero_sc_len: int,
    zero_tracers_p: cobj,
    zero_sc_p: cobj,
):
    zero_tracers = Ptr[float](zero_tracers_p)
    zero_sc = Ptr[float](zero_sc_p)

    for m in range(1, pcnst + 1):
        for i in range(1, pcols + 1):
            zero_tracers[_field2_idx(i, m, pcols)] = 0.0

    for i in range(1, zero_sc_len + 1):
        zero_sc[_idx(i)] = 0.0


@export
def tphysbc_trace_water_clip_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    pwtype: int,
    wtrc_nwset: int,
    wisotope_on: int,
    state_q_p: cobj,
    wtrc_iatype_p: cobj,
    tagged_p: cobj,
    rstd_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    wtrc_iatype = Ptr[int](wtrc_iatype_p)
    tagged = Ptr[int](tagged_p)
    rstd = Ptr[float](rstd_p)

    for i in range(1, ncol + 1):
        for k in range(1, pver + 1):
            for p in range(1, pwtype + 1):
                bulk_idx = wtrc_iatype[_field2_idx(1, p, wtrc_nwset)]
                bulk = state_q[_field3_idx(i, k, bulk_idx, pcols, pver)]
                for m in range(2, wtrc_nwset + 1):
                    tracer_idx = wtrc_iatype[_field2_idx(m, p, wtrc_nwset)]
                    state_idx = _field3_idx(i, k, tracer_idx, pcols, pver)
                    if wisotope_on != 0:
                        if state_q[state_idx] > 1.5 * bulk:
                            state_q[state_idx] = bulk
                    else:
                        if state_q[state_idx] > bulk:
                            if tagged[_idx(tracer_idx)] != 0:
                                state_q[state_idx] = bulk
                            else:
                                state_q[state_idx] = rstd[_idx(tracer_idx)] * bulk


@export
def tphysbc_dadadj_lq_init_codon(
    pcnst: int,
    lq_mask_p: cobj,
):
    lq_mask = Ptr[int](lq_mask_p)

    for m in range(1, pcnst + 1):
        lq_mask[_idx(m)] = 0

    if pcnst >= 1:
        lq_mask[0] = 1


@export
def tphysbc_state_batch_qini_snapshot_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    ixcldliq: int,
    ixcldice: int,
    state_q_p: cobj,
    qini_p: cobj,
    cldliqini_p: cobj,
    cldiceini_p: cobj,
):
    tphysbc_qini_snapshot_codon(
        ncol, pcols, pver, pcnst, ixcldliq, ixcldice, state_q_p, qini_p, cldliqini_p, cldiceini_p
    )


@export
def tphysbc_state_batch_dispatch_codon(
    stage: int,
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    ixcldliq: int,
    ixcldice: int,
    zero_sc_len: int,
    pwtype: int,
    wtrc_nwset: int,
    wisotope_on: int,
    ztodt: float,
    cpair_local: float,
    state_t_p: cobj,
    state_q_p: cobj,
    qini_p: cobj,
    cldliqini_p: cobj,
    cldiceini_p: cobj,
    tini_p: cobj,
    dtcore_p: cobj,
    tend_dtdt_p: cobj,
    tend_dudt_p: cobj,
    tend_dvdt_p: cobj,
    ptend_s_p: cobj,
    ptend_q_p: cobj,
    fracis_p: cobj,
    zero_tracers_p: cobj,
    zero_sc_p: cobj,
    wtrc_iatype_p: cobj,
    tagged_p: cobj,
    rstd_p: cobj,
    lq_mask_p: cobj,
):
    if stage == 1:
        tphysbc_zero_buffers_codon(pcols, pcnst, zero_sc_len, zero_tracers_p, zero_sc_p)
    elif stage == 2:
        tphysbc_trace_water_clip_codon(
            ncol, pcols, pver, pcnst, pwtype, wtrc_nwset, wisotope_on, state_q_p, wtrc_iatype_p, tagged_p, rstd_p
        )
    elif stage == 3:
        tphysbc_init_fields_codon(ncol, pcols, pver, pcnst, fracis_p, tend_dtdt_p, tend_dudt_p, tend_dvdt_p)
    elif stage == 4:
        tphysbc_tini_copy_codon(ncol, pcols, pver, state_t_p, tini_p)
    elif stage == 5:
        tphysbc_qini_snapshot_codon(
            ncol, pcols, pver, pcnst, ixcldliq, ixcldice, state_q_p, qini_p, cldliqini_p, cldiceini_p
        )
    elif stage == 6:
        tphysbc_dtcore_update_codon(ncol, pcols, pver, ztodt, tini_p, dtcore_p, tend_dtdt_p)
    elif stage == 7:
        tphysbc_dadadj_lq_init_codon(pcnst, lq_mask_p)
    elif stage == 8:
        tphysbc_dadadj_input_codon(ncol, pcols, pver, pcnst, state_t_p, state_q_p, ptend_s_p, ptend_q_p)
    elif stage == 9:
        tphysbc_dadadj_output_codon(
            ncol, pcols, pver, pcnst, ztodt, cpair_local, state_t_p, state_q_p, ptend_s_p, ptend_q_p
        )


@export
def tphysbc_state_batch_dadadj_input_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    state_t_p: cobj,
    state_q_p: cobj,
    ptend_s_p: cobj,
    ptend_q_p: cobj,
):
    tphysbc_dadadj_input_codon(ncol, pcols, pver, pcnst, state_t_p, state_q_p, ptend_s_p, ptend_q_p)


@export
def tphysbc_state_batch_dadadj_output_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    ztodt: float,
    cpair_local: float,
    state_t_p: cobj,
    state_q_p: cobj,
    ptend_s_p: cobj,
    ptend_q_p: cobj,
):
    tphysbc_dadadj_output_codon(
        ncol, pcols, pver, pcnst, ztodt, cpair_local, state_t_p, state_q_p, ptend_s_p, ptend_q_p
    )


@export
def tphysbc_state_batch_dtcore_update_codon(
    ncol: int,
    pcols: int,
    pver: int,
    ztodt: float,
    tini_p: cobj,
    dtcore_p: cobj,
    tend_dtdt_p: cobj,
):
    tphysbc_dtcore_update_codon(ncol, pcols, pver, ztodt, tini_p, dtcore_p, tend_dtdt_p)


@export
def tphysbc_state_batch_tini_copy_codon(
    ncol: int,
    pcols: int,
    pver: int,
    state_t_p: cobj,
    tini_p: cobj,
):
    tphysbc_tini_copy_codon(ncol, pcols, pver, state_t_p, tini_p)


@export
def tphysbc_state_batch_zero_buffers_codon(
    pcols: int,
    pcnst: int,
    zero_sc_len: int,
    zero_tracers_p: cobj,
    zero_sc_p: cobj,
):
    tphysbc_zero_buffers_codon(pcols, pcnst, zero_sc_len, zero_tracers_p, zero_sc_p)


@export
def tphysbc_state_batch_trace_water_clip_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    pwtype: int,
    wtrc_nwset: int,
    wisotope_on: int,
    state_q_p: cobj,
    wtrc_iatype_p: cobj,
    tagged_p: cobj,
    rstd_p: cobj,
):
    tphysbc_trace_water_clip_codon(
        ncol, pcols, pver, pcnst, pwtype, wtrc_nwset, wisotope_on, state_q_p, wtrc_iatype_p, tagged_p, rstd_p
    )


@export
def tphysbc_state_batch_dadadj_lq_init_codon(
    pcnst: int,
    lq_mask_p: cobj,
):
    tphysbc_dadadj_lq_init_codon(pcnst, lq_mask_p)


@export
def tphysbc_state_batch_init_fields_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    fracis_p: cobj,
    tend_dtdt_p: cobj,
    tend_dudt_p: cobj,
    tend_dvdt_p: cobj,
):
    tphysbc_init_fields_codon(ncol, pcols, pver, pcnst, fracis_p, tend_dtdt_p, tend_dudt_p, tend_dvdt_p)


@export
def phys_inidat_qpert_default_codon(
    pcols: int,
    chunk_count: int,
    found_i: int,
    tptr_p: cobj,
):
    if found_i != 0:
        return

    tptr = Ptr[float](tptr_p)

    for c in range(1, chunk_count + 1):
        for i in range(1, pcols + 1):
            tptr[_field2_idx(i, c, pcols)] = 0.0


@export
def phys_inidat_qpert_expand_codon(
    pcols: int,
    pcnst: int,
    chunk_count: int,
    tptr_p: cobj,
    tptr3d_2_p: cobj,
):
    tptr = Ptr[float](tptr_p)
    tptr3d_2 = Ptr[float](tptr3d_2_p)

    for c in range(1, chunk_count + 1):
        for m in range(1, pcnst + 1):
            for i in range(1, pcols + 1):
                tptr3d_2[_field3_idx(i, m, c, pcols, pcnst)] = 0.0

    for c in range(1, chunk_count + 1):
        for i in range(1, pcols + 1):
            tptr3d_2[_field3_idx(i, 1, c, pcols, pcnst)] = tptr[_field2_idx(i, c, pcols)]


@export
def phys_inidat_pblh_default_codon(
    pcols: int,
    chunk_count: int,
    found_i: int,
    tptr_p: cobj,
):
    if found_i != 0:
        return

    tptr = Ptr[float](tptr_p)

    for c in range(1, chunk_count + 1):
        for i in range(1, pcols + 1):
            tptr[_field2_idx(i, c, pcols)] = 0.0


@export
def phys_inidat_tpert_default_codon(
    pcols: int,
    chunk_count: int,
    found_i: int,
    tptr_p: cobj,
):
    if found_i != 0:
        return

    tptr = Ptr[float](tptr_p)

    for c in range(1, chunk_count + 1):
        for i in range(1, pcols + 1):
            tptr[_field2_idx(i, c, pcols)] = 0.0


@export
def phys_inidat_cush_default_codon(
    pcols: int,
    chunk_count: int,
    found_i: int,
    tptr_p: cobj,
    default_value: float,
):
    if found_i != 0:
        return

    tptr = Ptr[float](tptr_p)

    for c in range(1, chunk_count + 1):
        for i in range(1, pcols + 1):
            tptr[_field2_idx(i, c, pcols)] = default_value


@export
def phys_inidat_tke_default_codon(
    pcols: int,
    pverp: int,
    chunk_count: int,
    found_i: int,
    tptr3d_p: cobj,
    default_value: float,
):
    if found_i != 0:
        return

    tptr3d = Ptr[float](tptr3d_p)

    for c in range(1, chunk_count + 1):
        for k in range(1, pverp + 1):
            for i in range(1, pcols + 1):
                tptr3d[_field3_idx(i, k, c, pcols, pverp)] = default_value


@export
def phys_inidat_kvm_default_codon(
    pcols: int,
    pverp: int,
    chunk_count: int,
    found_i: int,
    tptr3d_p: cobj,
    default_value: float,
):
    if found_i != 0:
        return

    tptr3d = Ptr[float](tptr3d_p)

    for c in range(1, chunk_count + 1):
        for k in range(1, pverp + 1):
            for i in range(1, pcols + 1):
                tptr3d[_field3_idx(i, k, c, pcols, pverp)] = default_value


@export
def phys_inidat_kvh_default_codon(
    pcols: int,
    pverp: int,
    chunk_count: int,
    found_i: int,
    tptr3d_p: cobj,
    default_value: float,
):
    if found_i != 0:
        return

    tptr3d = Ptr[float](tptr3d_p)

    for c in range(1, chunk_count + 1):
        for k in range(1, pverp + 1):
            for i in range(1, pcols + 1):
                tptr3d[_field3_idx(i, k, c, pcols, pverp)] = default_value


@export
def phys_inidat_qcwat_default_codon(
    primary_found_i: int,
    fallback_found_i: int,
    init_source_p: cobj,
):
    init_source = Ptr[int](init_source_p)

    if primary_found_i != 0:
        init_source[0] = 1
    elif fallback_found_i != 0:
        init_source[0] = 2
    else:
        init_source[0] = 0


@export
def phys_inidat_iccwat_default_codon(
    primary_found_i: int,
    fallback_found_i: int,
    init_source_p: cobj,
):
    init_source = Ptr[int](init_source_p)

    if primary_found_i != 0:
        init_source[0] = 1
    elif fallback_found_i != 0:
        init_source[0] = 2
    else:
        init_source[0] = 0


@export
def phys_inidat_lcwat_default_codon(
    primary_found_i: int,
    cldice_found_i: int,
    cldliq_found_i: int,
    init_source_p: cobj,
):
    init_source = Ptr[int](init_source_p)

    if primary_found_i != 0:
        init_source[0] = 1
    elif cldice_found_i != 0 and cldliq_found_i != 0:
        init_source[0] = 2
    elif cldice_found_i != 0:
        init_source[0] = 3
    elif cldliq_found_i != 0:
        init_source[0] = 4
    else:
        init_source[0] = 0


@export
def phys_inidat_tcwat_default_codon(
    primary_found_i: int,
    init_source_p: cobj,
):
    init_source = Ptr[int](init_source_p)

    if primary_found_i != 0:
        init_source[0] = 1
    else:
        init_source[0] = 2


@export
def phys_inidat_cloud_default_codon(
    pcols: int,
    pver: int,
    chunk_count: int,
    found_i: int,
    tptr3d_p: cobj,
    default_value: float,
):
    if found_i != 0:
        return

    tptr3d = Ptr[float](tptr3d_p)

    for c in range(1, chunk_count + 1):
        for k in range(1, pver + 1):
            for i in range(1, pcols + 1):
                tptr3d[_field3_idx(i, k, c, pcols, pver)] = default_value


@export
def phys_inidat_concld_default_codon(
    pcols: int,
    pver: int,
    chunk_count: int,
    found_i: int,
    tptr3d_p: cobj,
    default_value: float,
):
    if found_i != 0:
        return

    tptr3d = Ptr[float](tptr3d_p)

    for c in range(1, chunk_count + 1):
        for k in range(1, pver + 1):
            for i in range(1, pcols + 1):
                tptr3d[_field3_idx(i, k, c, pcols, pver)] = default_value


@export
def phys_inidat_tbot_init_codon(
    pcols: int,
    tbot_p: cobj,
    posinf_local: float,
):
    tbot = Ptr[float](tbot_p)

    for i in range(1, pcols + 1):
        tbot[_idx(i)] = posinf_local


@export
def phys_inidat_batch_qpert_default_codon(
    pcols: int,
    chunk_count: int,
    found_i: int,
    tptr_p: cobj,
):
    phys_inidat_qpert_default_codon(pcols, chunk_count, found_i, tptr_p)


@export
def phys_inidat_batch_qpert_expand_codon(
    pcols: int,
    pcnst: int,
    chunk_count: int,
    tptr_p: cobj,
    tptr3d_2_p: cobj,
):
    phys_inidat_qpert_expand_codon(pcols, pcnst, chunk_count, tptr_p, tptr3d_2_p)


@export
def phys_inidat_batch_pblh_default_codon(
    pcols: int,
    chunk_count: int,
    found_i: int,
    tptr_p: cobj,
):
    phys_inidat_pblh_default_codon(pcols, chunk_count, found_i, tptr_p)


@export
def phys_inidat_batch_tpert_default_codon(
    pcols: int,
    chunk_count: int,
    found_i: int,
    tptr_p: cobj,
):
    phys_inidat_tpert_default_codon(pcols, chunk_count, found_i, tptr_p)


@export
def phys_inidat_batch_cush_default_codon(
    pcols: int,
    chunk_count: int,
    found_i: int,
    tptr_p: cobj,
    default_value: float,
):
    phys_inidat_cush_default_codon(pcols, chunk_count, found_i, tptr_p, default_value)


@export
def phys_inidat_batch_tke_default_codon(
    pcols: int,
    pverp: int,
    chunk_count: int,
    found_i: int,
    tptr3d_p: cobj,
    default_value: float,
):
    phys_inidat_tke_default_codon(pcols, pverp, chunk_count, found_i, tptr3d_p, default_value)


@export
def phys_inidat_batch_kvm_default_codon(
    pcols: int,
    pverp: int,
    chunk_count: int,
    found_i: int,
    tptr3d_p: cobj,
    default_value: float,
):
    phys_inidat_kvm_default_codon(pcols, pverp, chunk_count, found_i, tptr3d_p, default_value)


@export
def phys_inidat_batch_kvh_default_codon(
    pcols: int,
    pverp: int,
    chunk_count: int,
    found_i: int,
    tptr3d_p: cobj,
    default_value: float,
):
    phys_inidat_kvh_default_codon(pcols, pverp, chunk_count, found_i, tptr3d_p, default_value)


@export
def phys_inidat_batch_qcwat_default_codon(
    primary_found_i: int,
    fallback_found_i: int,
    init_source_p: cobj,
):
    phys_inidat_qcwat_default_codon(primary_found_i, fallback_found_i, init_source_p)


@export
def phys_inidat_batch_iccwat_default_codon(
    primary_found_i: int,
    fallback_found_i: int,
    init_source_p: cobj,
):
    phys_inidat_iccwat_default_codon(primary_found_i, fallback_found_i, init_source_p)


@export
def phys_inidat_batch_lcwat_default_codon(
    primary_found_i: int,
    cldice_found_i: int,
    cldliq_found_i: int,
    init_source_p: cobj,
):
    phys_inidat_lcwat_default_codon(primary_found_i, cldice_found_i, cldliq_found_i, init_source_p)


@export
def phys_inidat_batch_tcwat_default_codon(
    primary_found_i: int,
    init_source_p: cobj,
):
    phys_inidat_tcwat_default_codon(primary_found_i, init_source_p)


@export
def phys_inidat_batch_cloud_default_codon(
    pcols: int,
    pver: int,
    chunk_count: int,
    found_i: int,
    tptr3d_p: cobj,
    default_value: float,
):
    phys_inidat_cloud_default_codon(pcols, pver, chunk_count, found_i, tptr3d_p, default_value)


@export
def phys_inidat_batch_concld_default_codon(
    pcols: int,
    pver: int,
    chunk_count: int,
    found_i: int,
    tptr3d_p: cobj,
    default_value: float,
):
    phys_inidat_concld_default_codon(pcols, pver, chunk_count, found_i, tptr3d_p, default_value)


@export
def phys_inidat_batch_tbot_init_codon(
    pcols: int,
    tbot_p: cobj,
    posinf_local: float,
):
    phys_inidat_tbot_init_codon(pcols, tbot_p, posinf_local)
