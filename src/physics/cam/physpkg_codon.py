from math import exp, floor, log, sqrt

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


@export
def tropopause_output_prep_stage_dispatch_codon(
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
    tropopause_output_prep_codon(
        ncol,
        pcols,
        pver,
        notfound,
        fillvalue,
        trop_lev_p,
        trop_z_p,
        state_zm_p,
        trop_pdf_p,
        trop_found_p,
        trop_dz_p,
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
def tropopause_twmo_stage_dispatch_codon(
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
    tropopause_twmo_codon(
        ncol,
        pcols,
        pver,
        notfound,
        write_tropp,
        write_tropt,
        write_tropz,
        cnst_kap,
        cnst_ka1,
        cnst_faktor,
        state_t_p,
        state_pmid_p,
        state_pint_p,
        state_zm_p,
        state_zi_p,
        trop_lev_p,
        trop_p_p,
        trop_t_p,
        trop_z_p,
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
    _physpkg_tphys_batch_dispatch(
        1,
        stage,
        mode,
        ncol,
        pcols,
        pver,
        pcnst,
        cld_macmic_num_steps,
        ixcldliq,
        ixcldice,
        0,
        0,
        0,
        0,
        ztodt,
        latvap_local,
        latice_local,
        rhoh2o_local,
        0.0,
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
        tend_flx_net_p,
        cam_in_shf_p,
        cam_out_precc_p,
        cam_out_precl_p,
        cam_out_precsc_p,
        cam_out_precsl_p,
        state_t_p,
        tini_p,
        tend_dtdt_p,
        dtcore_p,
        tmp_t_p,
        state_q_p,
        tmp_q_p,
        tmp_cldliq_p,
        tmp_cldice_p,
        a_p,
        b_p,
        out_p,
        rliq_p,
        det_s_p,
        flx_cnd_p,
        flx_heat_p,
        shf_p,
        net_flx_p,
        state_t_p,
        state_q_p,
        tmp_q_p,
        tmp_cldliq_p,
        tmp_cldice_p,
        tini_p,
        dtcore_p,
        tend_dtdt_p,
        tend_dtdt_p,
        tend_dtdt_p,
        state_t_p,
        state_q_p,
        a_p,
        state_q_p,
        state_t_p,
        state_q_p,
        state_q_p,
        state_q_p,
        state_q_p,
    )


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


def _physpkg_tphys_batch_dispatch(
    group: int,
    stage: int,
    mode: int,
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    cld_macmic_num_steps: int,
    ixcldliq: int,
    ixcldice: int,
    zero_sc_len: int,
    pwtype: int,
    wtrc_nwset: int,
    wisotope_on: int,
    ztodt: float,
    latvap_local: float,
    latice_local: float,
    rhoh2o_local: float,
    cpair_local: float,
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
    flux_state_t_p: cobj,
    flux_tini_p: cobj,
    flux_tend_dtdt_p: cobj,
    flux_dtcore_p: cobj,
    tmp_t_p: cobj,
    flux_state_q_p: cobj,
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
    bc_state_t_p: cobj,
    bc_state_q_p: cobj,
    qini_p: cobj,
    cldliqini_p: cobj,
    cldiceini_p: cobj,
    bc_tini_p: cobj,
    bc_dtcore_p: cobj,
    bc_tend_dtdt_p: cobj,
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
    if group == 1:
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
            tphysac_t_update_codon(
                ncol,
                pcols,
                pver,
                ztodt,
                flux_state_t_p,
                flux_tini_p,
                flux_tend_dtdt_p,
                flux_dtcore_p,
                tmp_t_p,
            )
        elif stage == 4:
            tphysac_q_snapshot_codon(
                ncol,
                pcols,
                pver,
                pcnst,
                ixcldliq,
                ixcldice,
                flux_state_q_p,
                tmp_q_p,
                tmp_cldliq_p,
                tmp_cldice_p,
            )
        elif stage == 5:
            tphysbc_flx_cnd_sum_codon(ncol, pcols, a_p, b_p, out_p)
        elif stage == 6:
            tphysbc_macrop_fluxes_codon(mode, ncol, pcols, rliq_p, det_s_p, flx_cnd_p, flx_heat_p, shf_p)
        elif stage == 7:
            tphysbc_radheat_flx_net_codon(ncol, pcols, tend_flx_net_p, net_flx_p)
    elif group == 2:
        if stage == 1:
            tphysbc_zero_buffers_codon(pcols, pcnst, zero_sc_len, zero_tracers_p, zero_sc_p)
        elif stage == 2:
            tphysbc_trace_water_clip_codon(
                ncol,
                pcols,
                pver,
                pcnst,
                pwtype,
                wtrc_nwset,
                wisotope_on,
                bc_state_q_p,
                wtrc_iatype_p,
                tagged_p,
                rstd_p,
            )
        elif stage == 3:
            tphysbc_init_fields_codon(ncol, pcols, pver, pcnst, fracis_p, bc_tend_dtdt_p, tend_dudt_p, tend_dvdt_p)
        elif stage == 4:
            tphysbc_tini_copy_codon(ncol, pcols, pver, bc_state_t_p, bc_tini_p)
        elif stage == 5:
            tphysbc_qini_snapshot_codon(
                ncol,
                pcols,
                pver,
                pcnst,
                ixcldliq,
                ixcldice,
                bc_state_q_p,
                qini_p,
                cldliqini_p,
                cldiceini_p,
            )
        elif stage == 6:
            tphysbc_dtcore_update_codon(ncol, pcols, pver, ztodt, bc_tini_p, bc_dtcore_p, bc_tend_dtdt_p)
        elif stage == 7:
            tphysbc_dadadj_lq_init_codon(pcnst, lq_mask_p)
        elif stage == 8:
            tphysbc_dadadj_input_codon(ncol, pcols, pver, pcnst, bc_state_t_p, bc_state_q_p, ptend_s_p, ptend_q_p)
        elif stage == 9:
            tphysbc_dadadj_output_codon(
                ncol,
                pcols,
                pver,
                pcnst,
                ztodt,
                cpair_local,
                bc_state_t_p,
                bc_state_q_p,
                ptend_s_p,
                ptend_q_p,
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
    _physpkg_tphys_batch_dispatch(
        2,
        stage,
        0,
        ncol,
        pcols,
        pver,
        pcnst,
        0,
        ixcldliq,
        ixcldice,
        zero_sc_len,
        pwtype,
        wtrc_nwset,
        wisotope_on,
        ztodt,
        0.0,
        0.0,
        0.0,
        cpair_local,
        state_t_p,
        state_q_p,
        qini_p,
        cldliqini_p,
        cldiceini_p,
        tini_p,
        dtcore_p,
        tend_dtdt_p,
        tend_dudt_p,
        tend_dvdt_p,
        ptend_s_p,
        ptend_q_p,
        fracis_p,
        state_t_p,
        state_q_p,
        qini_p,
        cldliqini_p,
        cldiceini_p,
        state_t_p,
        tini_p,
        tend_dtdt_p,
        dtcore_p,
        tend_dudt_p,
        state_q_p,
        ptend_q_p,
        cldliqini_p,
        cldiceini_p,
        fracis_p,
        zero_tracers_p,
        zero_sc_p,
        qini_p,
        cldliqini_p,
        cldiceini_p,
        state_t_p,
        state_q_p,
        qini_p,
        state_t_p,
        state_q_p,
        qini_p,
        cldliqini_p,
        cldiceini_p,
        tini_p,
        dtcore_p,
        tend_dtdt_p,
        tend_dudt_p,
        tend_dvdt_p,
        ptend_s_p,
        ptend_q_p,
        fracis_p,
        zero_tracers_p,
        zero_sc_p,
        wtrc_iatype_p,
        tagged_p,
        rstd_p,
        lq_mask_p,
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
def phys_inidat_batch_dispatch_codon(
    stage: int,
    pcols: int,
    pver: int,
    pverp: int,
    pcnst: int,
    chunk_count: int,
    found1_i: int,
    found2_i: int,
    found3_i: int,
    default_value: float,
    tptr_p: cobj,
    tptr3d_p: cobj,
    tptr3d_2_p: cobj,
    init_source_p: cobj,
    tbot_p: cobj,
):
    if stage == 1:
        phys_inidat_qpert_default_codon(pcols, chunk_count, found1_i, tptr_p)
    elif stage == 2:
        phys_inidat_qpert_expand_codon(pcols, pcnst, chunk_count, tptr_p, tptr3d_2_p)
    elif stage == 3:
        phys_inidat_pblh_default_codon(pcols, chunk_count, found1_i, tptr_p)
    elif stage == 4:
        phys_inidat_tpert_default_codon(pcols, chunk_count, found1_i, tptr_p)
    elif stage == 5:
        phys_inidat_cush_default_codon(pcols, chunk_count, found1_i, tptr_p, default_value)
    elif stage == 6:
        phys_inidat_tke_default_codon(pcols, pverp, chunk_count, found1_i, tptr3d_p, default_value)
    elif stage == 7:
        phys_inidat_kvm_default_codon(pcols, pverp, chunk_count, found1_i, tptr3d_p, default_value)
    elif stage == 8:
        phys_inidat_kvh_default_codon(pcols, pverp, chunk_count, found1_i, tptr3d_p, default_value)
    elif stage == 9:
        phys_inidat_qcwat_default_codon(found1_i, found2_i, init_source_p)
    elif stage == 10:
        phys_inidat_iccwat_default_codon(found1_i, found2_i, init_source_p)
    elif stage == 11:
        phys_inidat_lcwat_default_codon(found1_i, found2_i, found3_i, init_source_p)
    elif stage == 12:
        phys_inidat_tcwat_default_codon(found1_i, init_source_p)
    elif stage == 13:
        phys_inidat_cloud_default_codon(pcols, pver, chunk_count, found1_i, tptr3d_p, default_value)
    elif stage == 14:
        phys_inidat_concld_default_codon(pcols, pver, chunk_count, found1_i, tptr3d_p, default_value)
    elif stage == 15:
        phys_inidat_tbot_init_codon(pcols, tbot_p, default_value)


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


@inline
def _physics_types_zero_1d(n: int, arr: Ptr[float]):
    for i in range(0, n):
        arr[i] = 0.0


@inline
def _physics_types_zero_2d(psetcols: int, pver: int, arr: Ptr[float]):
    for k in range(1, pver + 1):
        for i in range(1, psetcols + 1):
            arr[_field2_idx(i, k, psetcols)] = 0.0


@inline
def _physics_types_zero_3d(psetcols: int, pver: int, pcnst: int, arr: Ptr[float]):
    for m in range(1, pcnst + 1):
        for k in range(1, pver + 1):
            for i in range(1, psetcols + 1):
                arr[_field3_idx(i, k, m, psetcols, pver)] = 0.0


@inline
def _physics_types_copy_1d(n: int, src: Ptr[float], dst: Ptr[float]):
    for i in range(0, n):
        dst[i] = src[i]


@inline
def _physics_types_copy_2d(psetcols: int, pver: int, src: Ptr[float], dst: Ptr[float]):
    for k in range(1, pver + 1):
        for i in range(1, psetcols + 1):
            dst[_field2_idx(i, k, psetcols)] = src[_field2_idx(i, k, psetcols)]


@inline
def _physics_types_copy_3d(psetcols: int, pver: int, pcnst: int, src: Ptr[float], dst: Ptr[float]):
    for m in range(1, pcnst + 1):
        for k in range(1, pver + 1):
            for i in range(1, psetcols + 1):
                dst[_field3_idx(i, k, m, psetcols, pver)] = src[_field3_idx(i, k, m, psetcols, pver)]


@export
def physics_tend_init_codon(
    psetcols: int,
    pver: int,
    dtdt_p: cobj,
    dudt_p: cobj,
    dvdt_p: cobj,
    flx_net_p: cobj,
    te_tnd_p: cobj,
    tw_tnd_p: cobj,
):
    dtdt = Ptr[float](dtdt_p)
    dudt = Ptr[float](dudt_p)
    dvdt = Ptr[float](dvdt_p)
    flx_net = Ptr[float](flx_net_p)
    te_tnd = Ptr[float](te_tnd_p)
    tw_tnd = Ptr[float](tw_tnd_p)

    _physics_types_zero_2d(psetcols, pver, dtdt)
    _physics_types_zero_2d(psetcols, pver, dudt)
    _physics_types_zero_2d(psetcols, pver, dvdt)
    _physics_types_zero_1d(psetcols, flx_net)
    _physics_types_zero_1d(psetcols, te_tnd)
    _physics_types_zero_1d(psetcols, tw_tnd)


@export
def physics_ptend_reset_s_codon(psetcols: int, pver: int, s_p: cobj, hflux_srf_p: cobj, hflux_top_p: cobj):
    _physics_types_zero_2d(psetcols, pver, Ptr[float](s_p))
    _physics_types_zero_1d(psetcols, Ptr[float](hflux_srf_p))
    _physics_types_zero_1d(psetcols, Ptr[float](hflux_top_p))


@export
def physics_ptend_reset_u_codon(psetcols: int, pver: int, u_p: cobj, taux_srf_p: cobj, taux_top_p: cobj):
    _physics_types_zero_2d(psetcols, pver, Ptr[float](u_p))
    _physics_types_zero_1d(psetcols, Ptr[float](taux_srf_p))
    _physics_types_zero_1d(psetcols, Ptr[float](taux_top_p))


@export
def physics_ptend_reset_v_codon(psetcols: int, pver: int, v_p: cobj, tauy_srf_p: cobj, tauy_top_p: cobj):
    _physics_types_zero_2d(psetcols, pver, Ptr[float](v_p))
    _physics_types_zero_1d(psetcols, Ptr[float](tauy_srf_p))
    _physics_types_zero_1d(psetcols, Ptr[float](tauy_top_p))


@export
def physics_ptend_reset_q_codon(psetcols: int, pver: int, pcnst: int, q_p: cobj, cflx_srf_p: cobj, cflx_top_p: cobj):
    _physics_types_zero_3d(psetcols, pver, pcnst, Ptr[float](q_p))
    _physics_types_zero_2d(psetcols, pcnst, Ptr[float](cflx_srf_p))
    _physics_types_zero_2d(psetcols, pcnst, Ptr[float](cflx_top_p))


@export
def physics_ptend_copy_s_codon(
    psetcols: int,
    pver: int,
    src_s_p: cobj,
    src_hflux_srf_p: cobj,
    src_hflux_top_p: cobj,
    dst_s_p: cobj,
    dst_hflux_srf_p: cobj,
    dst_hflux_top_p: cobj,
):
    _physics_types_copy_2d(psetcols, pver, Ptr[float](src_s_p), Ptr[float](dst_s_p))
    _physics_types_copy_1d(psetcols, Ptr[float](src_hflux_srf_p), Ptr[float](dst_hflux_srf_p))
    _physics_types_copy_1d(psetcols, Ptr[float](src_hflux_top_p), Ptr[float](dst_hflux_top_p))


@export
def physics_ptend_copy_u_codon(
    psetcols: int,
    pver: int,
    src_u_p: cobj,
    src_taux_srf_p: cobj,
    src_taux_top_p: cobj,
    dst_u_p: cobj,
    dst_taux_srf_p: cobj,
    dst_taux_top_p: cobj,
):
    _physics_types_copy_2d(psetcols, pver, Ptr[float](src_u_p), Ptr[float](dst_u_p))
    _physics_types_copy_1d(psetcols, Ptr[float](src_taux_srf_p), Ptr[float](dst_taux_srf_p))
    _physics_types_copy_1d(psetcols, Ptr[float](src_taux_top_p), Ptr[float](dst_taux_top_p))


@export
def physics_ptend_copy_v_codon(
    psetcols: int,
    pver: int,
    src_v_p: cobj,
    src_tauy_srf_p: cobj,
    src_tauy_top_p: cobj,
    dst_v_p: cobj,
    dst_tauy_srf_p: cobj,
    dst_tauy_top_p: cobj,
):
    _physics_types_copy_2d(psetcols, pver, Ptr[float](src_v_p), Ptr[float](dst_v_p))
    _physics_types_copy_1d(psetcols, Ptr[float](src_tauy_srf_p), Ptr[float](dst_tauy_srf_p))
    _physics_types_copy_1d(psetcols, Ptr[float](src_tauy_top_p), Ptr[float](dst_tauy_top_p))


@export
def physics_ptend_copy_q_codon(
    psetcols: int,
    pver: int,
    pcnst: int,
    src_q_p: cobj,
    src_cflx_srf_p: cobj,
    src_cflx_top_p: cobj,
    dst_q_p: cobj,
    dst_cflx_srf_p: cobj,
    dst_cflx_top_p: cobj,
):
    _physics_types_copy_3d(psetcols, pver, pcnst, Ptr[float](src_q_p), Ptr[float](dst_q_p))
    _physics_types_copy_2d(psetcols, pcnst, Ptr[float](src_cflx_srf_p), Ptr[float](dst_cflx_srf_p))
    _physics_types_copy_2d(psetcols, pcnst, Ptr[float](src_cflx_top_p), Ptr[float](dst_cflx_top_p))


@export
def physics_ptend_scale_field_codon(
    ncol: int,
    psetcols: int,
    top_level: int,
    bot_level: int,
    fac: float,
    field_p: cobj,
    flx_srf_p: cobj,
    flx_top_p: cobj,
):
    field = Ptr[float](field_p)
    flx_srf = Ptr[float](flx_srf_p)
    flx_top = Ptr[float](flx_top_p)

    for k in range(top_level, bot_level + 1):
        for i in range(1, ncol + 1):
            field[_field2_idx(i, k, psetcols)] = field[_field2_idx(i, k, psetcols)] * fac

    for i in range(1, ncol + 1):
        flx_srf[_idx(i)] = flx_srf[_idx(i)] * fac
        flx_top[_idx(i)] = flx_top[_idx(i)] * fac


@export
def phys_grid_get_gcol_all_codon(ncols: int, out_dim: int, src_p: cobj, dst_p: cobj):
    src = Ptr[i32](src_p)
    dst = Ptr[i32](dst_p)

    for i in range(out_dim):
        dst[i] = i32(-1)

    for i in range(ncols):
        dst[i] = src[i]


@export
def phys_grid_get_gcol_vec_codon(lth: int, cols_p: cobj, src_p: cobj, dst_p: cobj):
    cols = Ptr[i32](cols_p)
    src = Ptr[i32](src_p)
    dst = Ptr[i32](dst_p)

    for i in range(lth):
        col = int(cols[i])
        dst[i] = src[col - 1]


@export
def phys_grid_get_int_all_codon(ncols: int, src_p: cobj, dst_p: cobj):
    src = Ptr[i32](src_p)
    dst = Ptr[i32](dst_p)

    for i in range(ncols):
        dst[i] = src[i]


@export
def phys_grid_get_int_vec_codon(lth: int, cols_p: cobj, src_p: cobj, dst_p: cobj):
    cols = Ptr[i32](cols_p)
    src = Ptr[i32](src_p)
    dst = Ptr[i32](dst_p)

    for i in range(lth):
        col = int(cols[i])
        dst[i] = src[col - 1]


@export
def phys_grid_get_lon_all_codon(
    ncols: int,
    lat_p: cobj,
    gcol_p: cobj,
    dyn_to_latlon_gcol_map_p: cobj,
    clat_p_idx_p: cobj,
    dst_p: cobj,
):
    lat_src = Ptr[i32](lat_p)
    gcol_src = Ptr[i32](gcol_p)
    dyn_to_latlon_gcol_map = Ptr[i32](dyn_to_latlon_gcol_map_p)
    clat_p_idx = Ptr[i32](clat_p_idx_p)
    dst = Ptr[i32](dst_p)

    for i in range(ncols):
        lat = int(lat_src[i])
        gcol = int(dyn_to_latlon_gcol_map[int(gcol_src[i]) - 1])
        dst[i] = i32((gcol - int(clat_p_idx[lat - 1])) + 1)


@export
def phys_grid_get_lon_vec_codon(
    lth: int,
    cols_p: cobj,
    lat_p: cobj,
    gcol_p: cobj,
    dyn_to_latlon_gcol_map_p: cobj,
    clat_p_idx_p: cobj,
    dst_p: cobj,
):
    cols = Ptr[i32](cols_p)
    lat_src = Ptr[i32](lat_p)
    gcol_src = Ptr[i32](gcol_p)
    dyn_to_latlon_gcol_map = Ptr[i32](dyn_to_latlon_gcol_map_p)
    clat_p_idx = Ptr[i32](clat_p_idx_p)
    dst = Ptr[i32](dst_p)

    for i in range(lth):
        col = int(cols[i])
        lat = int(lat_src[col - 1])
        gcol = int(dyn_to_latlon_gcol_map[int(gcol_src[i]) - 1])
        dst[i] = i32((gcol - int(clat_p_idx[lat - 1])) + 1)


@export
def phys_grid_get_real_all_codon(ncols: int, src_p: cobj, dst_p: cobj):
    src = Ptr[float](src_p)
    dst = Ptr[float](dst_p)

    for i in range(ncols):
        dst[i] = src[i]


@export
def phys_grid_get_lookup_real_all_codon(ncols: int, idx_p: cobj, lookup_p: cobj, dst_p: cobj):
    idx = Ptr[i32](idx_p)
    lookup = Ptr[float](lookup_p)
    dst = Ptr[float](dst_p)

    for i in range(ncols):
        dst[i] = lookup[int(idx[i]) - 1]


@export
def phys_grid_get_lookup_real_vec_codon(
    lth: int,
    cols_p: cobj,
    idx_p: cobj,
    lookup_p: cobj,
    dst_p: cobj,
):
    cols = Ptr[i32](cols_p)
    idx = Ptr[i32](idx_p)
    lookup = Ptr[float](lookup_p)
    dst = Ptr[float](dst_p)

    for i in range(lth):
        col = int(cols[i])
        dst[i] = lookup[int(idx[col - 1]) - 1]


@export
def restart_physics_fill_tail_codon(
    ncol: int,
    pcols: int,
    fillvalue: float,
    fsnt_p: cobj,
    fsns_p: cobj,
    fsds_p: cobj,
    flnt_p: cobj,
    flns_p: cobj,
    landm_p: cobj,
    sgh_p: cobj,
    sgh30_p: cobj,
    trefmxav_p: cobj,
    trefmnav_p: cobj,
):
    fsnt = Ptr[float](fsnt_p)
    fsns = Ptr[float](fsns_p)
    fsds = Ptr[float](fsds_p)
    flnt = Ptr[float](flnt_p)
    flns = Ptr[float](flns_p)
    landm = Ptr[float](landm_p)
    sgh = Ptr[float](sgh_p)
    sgh30 = Ptr[float](sgh30_p)
    trefmxav = Ptr[float](trefmxav_p)
    trefmnav = Ptr[float](trefmnav_p)

    for i in range(ncol + 1, pcols + 1):
        idx = i - 1
        fsnt[idx] = fillvalue
        fsns[idx] = fillvalue
        fsds[idx] = fillvalue
        flnt[idx] = fillvalue
        flns[idx] = fillvalue
        landm[idx] = fillvalue
        sgh[idx] = fillvalue
        sgh30[idx] = fillvalue
        trefmxav[idx] = fillvalue
        trefmnav[idx] = fillvalue


@export
def restart_physics_tmpfield_fill_codon(total_len: int, fillvalue: float, tmpfield_p: cobj):
    tmpfield = Ptr[float](tmpfield_p)

    for i in range(total_len):
        tmpfield[i] = fillvalue


@export
def restart_physics_pack_chunk_field_codon(
    ncol: int,
    pcols: int,
    chunk_pos: int,
    field_p: cobj,
    tmpfield_p: cobj,
):
    field = Ptr[float](field_p)
    tmpfield = Ptr[float](tmpfield_p)
    offset = (chunk_pos - 1) * pcols

    for j in range(1, pcols + 1):
        if j <= ncol:
            tmpfield[offset + j - 1] = field[j - 1]


@inline
def _geopotential_idx(i: int, k: int, ld: int) -> int:
    """geopotential arrays declared as (ld, vertical_level)"""
    return (i - 1) + (k - 1) * ld


@export
def geopotential_dse_codon(
    ncol: int,
    ld: int,
    pver: int,
    pverp: int,
    fvdyn: int,
    gravit: float,
    piln_p: cobj,
    pint_p: cobj,
    pmid_p: cobj,
    pdel_p: cobj,
    rpdel_p: cobj,
    dse_p: cobj,
    q_p: cobj,
    phis_p: cobj,
    rair_p: cobj,
    cpair_p: cobj,
    zvir_p: cobj,
    t_p: cobj,
    zi_p: cobj,
    zm_p: cobj,
):
    piln = Ptr[float](piln_p)
    pint = Ptr[float](pint_p)
    pmid = Ptr[float](pmid_p)
    pdel = Ptr[float](pdel_p)
    rpdel = Ptr[float](rpdel_p)
    dse = Ptr[float](dse_p)
    q = Ptr[float](q_p)
    phis = Ptr[float](phis_p)
    rair = Ptr[float](rair_p)
    cpair = Ptr[float](cpair_p)
    zvir = Ptr[float](zvir_p)
    t = Ptr[float](t_p)
    zi = Ptr[float](zi_p)
    zm = Ptr[float](zm_p)

    for i in range(1, ncol + 1):
        zi[_geopotential_idx(i, pverp, ld)] = 0.0

    for k in range(pver, 0, -1):
        for i in range(1, ncol + 1):
            idx = _geopotential_idx(i, k, ld)
            if fvdyn != 0:
                hkl = piln[_geopotential_idx(i, k + 1, ld)] - piln[idx]
                hkk = 1.0 - pint[idx] * hkl * rpdel[idx]
            else:
                hkl = pdel[idx] / pmid[idx]
                hkk = 0.5 * hkl

            tvfac = 1.0 + zvir[idx] * q[idx]
            rog = rair[idx] / gravit
            tv = (dse[idx] - phis[i - 1] - gravit * zi[_geopotential_idx(i, k + 1, ld)]) / (
                (cpair[idx] / tvfac) + rair[idx] * hkk
            )

            t[idx] = tv / tvfac
            zm[idx] = zi[_geopotential_idx(i, k + 1, ld)] + rog * tv * hkk
            zi[idx] = zi[_geopotential_idx(i, k + 1, ld)] + rog * tv * hkl


@export
def geopotential_t_codon(
    ncol: int,
    ld: int,
    pver: int,
    pverp: int,
    fvdyn: int,
    gravit: float,
    piln_p: cobj,
    pint_p: cobj,
    pmid_p: cobj,
    pdel_p: cobj,
    rpdel_p: cobj,
    t_p: cobj,
    q_p: cobj,
    rair_p: cobj,
    zvir_p: cobj,
    zi_p: cobj,
    zm_p: cobj,
):
    piln = Ptr[float](piln_p)
    pint = Ptr[float](pint_p)
    pmid = Ptr[float](pmid_p)
    pdel = Ptr[float](pdel_p)
    rpdel = Ptr[float](rpdel_p)
    t = Ptr[float](t_p)
    q = Ptr[float](q_p)
    rair = Ptr[float](rair_p)
    zvir = Ptr[float](zvir_p)
    zi = Ptr[float](zi_p)
    zm = Ptr[float](zm_p)

    for i in range(1, ncol + 1):
        zi[_geopotential_idx(i, pverp, ld)] = 0.0

    for k in range(pver, 0, -1):
        for i in range(1, ncol + 1):
            idx = _geopotential_idx(i, k, ld)
            if fvdyn != 0:
                hkl = piln[_geopotential_idx(i, k + 1, ld)] - piln[idx]
                hkk = 1.0 - pint[idx] * hkl * rpdel[idx]
            else:
                hkl = pdel[idx] / pmid[idx]
                hkk = 0.5 * hkl

            tvfac = 1.0 + zvir[idx] * q[idx]
            tv = t[idx] * tvfac
            rog = rair[idx] / gravit

            zm[idx] = zi[_geopotential_idx(i, k + 1, ld)] + rog * tv * hkk
            zi[idx] = zi[_geopotential_idx(i, k + 1, ld)] + rog * tv * hkl


@export
def comsrf_initialize_fields_codon(
    total_len: int,
    nan_value: float,
    landm_p: cobj,
    sgh_p: cobj,
    sgh30_p: cobj,
    fsns_p: cobj,
    fsds_p: cobj,
    fsnt_p: cobj,
    flns_p: cobj,
    flnt_p: cobj,
    srfrpdel_p: cobj,
    psm1_p: cobj,
    prcsnw_p: cobj,
    trefmxav_p: cobj,
    trefmnav_p: cobj,
):
    landm = Ptr[float](landm_p)
    sgh = Ptr[float](sgh_p)
    sgh30 = Ptr[float](sgh30_p)
    fsns = Ptr[float](fsns_p)
    fsds = Ptr[float](fsds_p)
    fsnt = Ptr[float](fsnt_p)
    flns = Ptr[float](flns_p)
    flnt = Ptr[float](flnt_p)
    srfrpdel = Ptr[float](srfrpdel_p)
    psm1 = Ptr[float](psm1_p)
    prcsnw = Ptr[float](prcsnw_p)
    trefmxav = Ptr[float](trefmxav_p)
    trefmnav = Ptr[float](trefmnav_p)

    for idx in range(total_len):
        landm[idx] = nan_value
        sgh[idx] = nan_value
        sgh30[idx] = nan_value
        fsns[idx] = nan_value
        fsds[idx] = nan_value
        fsnt[idx] = nan_value
        flns[idx] = nan_value
        flnt[idx] = nan_value
        srfrpdel[idx] = nan_value
        psm1[idx] = nan_value
        prcsnw[idx] = nan_value
        trefmxav[idx] = -1.0e36
        trefmnav[idx] = 1.0e36


@inline
def _ref_pres_press_lim_idx_top(pver: int, p: float, pref_mid: Ptr[float]) -> int:
    k_lim = pver + 1
    for k in range(1, pver + 1):
        if pref_mid[k - 1] > p:
            k_lim = k
            break
    return k_lim


@inline
def _ref_pres_press_lim_idx_bottom(pver: int, p: float, pref_mid: Ptr[float]) -> int:
    k_lim = 0
    for k in range(pver, 0, -1):
        if pref_mid[k - 1] < p:
            k_lim = k
            break
    return k_lim


@export
def ref_pres_init_finalize_codon(
    pver: int,
    pverp: int,
    trop_cloud_top_press: float,
    clim_modal_aero_top_press: float,
    do_molec_press: float,
    molec_diff_bot_press: float,
    pref_edge_p: cobj,
    pref_mid_p: cobj,
    pref_mid_norm_p: cobj,
    scalar_out_p: cobj,
    int_out_p: cobj,
    flag_out_p: cobj,
):
    pref_edge = Ptr[float](pref_edge_p)
    pref_mid = Ptr[float](pref_mid_p)
    pref_mid_norm = Ptr[float](pref_mid_norm_p)
    scalar_out = Ptr[float](scalar_out_p)
    int_out = Ptr[int](int_out_p)
    flag_out = Ptr[int](flag_out_p)

    ptop_ref = pref_edge[0]
    psurf_ref = pref_edge[pverp - 1]

    scalar_out[0] = ptop_ref
    scalar_out[1] = psurf_ref

    for k in range(1, pver + 1):
        pref_mid_norm[k - 1] = pref_mid[k - 1] / psurf_ref

    int_out[0] = _ref_pres_press_lim_idx_top(pver, trop_cloud_top_press, pref_mid)
    int_out[1] = _ref_pres_press_lim_idx_top(pver, clim_modal_aero_top_press, pref_mid)

    if ptop_ref < do_molec_press:
        flag_out[0] = 1
        int_out[2] = _ref_pres_press_lim_idx_bottom(pver, molec_diff_bot_press, pref_mid)
    else:
        flag_out[0] = 0
        int_out[2] = 0


@inline
def _trb_mtn_idx(i: int, k: int, pcols: int) -> int:
    """trb_mtn_stress arrays declared as (pcols,pver)."""
    return (i - 1) + (k - 1) * pcols


@export
def trb_mtn_stress_compute_codon(
    pcols: int,
    pver: int,
    ncol: int,
    orocnst: float,
    z0fac: float,
    karman: float,
    gravit: float,
    rair: float,
    u_p: cobj,
    v_p: cobj,
    t_p: cobj,
    pmid_p: cobj,
    exner_p: cobj,
    zm_p: cobj,
    sgh_p: cobj,
    landfrac_p: cobj,
    ksrf_p: cobj,
    taux_p: cobj,
    tauy_p: cobj,
):
    u = Ptr[float](u_p)
    v = Ptr[float](v_p)
    t = Ptr[float](t_p)
    pmid = Ptr[float](pmid_p)
    exner = Ptr[float](exner_p)
    zm = Ptr[float](zm_p)
    sgh = Ptr[float](sgh_p)
    landfrac = Ptr[float](landfrac_p)
    ksrf = Ptr[float](ksrf_p)
    taux = Ptr[float](taux_p)
    tauy = Ptr[float](tauy_p)

    horomin = 1.0
    z0max = 100.0
    dv2min = 0.01
    kt = pver - 1
    kb = pver

    for i in range(1, ncol + 1):
        horo = orocnst * sgh[i - 1]

        if horo < horomin:
            ksrf[i - 1] = 0.0
            taux[i - 1] = 0.0
            tauy[i - 1] = 0.0
        else:
            z0oro = min(z0fac * horo, z0max)
            cd = (karman / log((zm[_trb_mtn_idx(i, pver, pcols)] + z0oro) / z0oro)) ** 2

            dv2 = max(
                (u[_trb_mtn_idx(i, kt, pcols)] - u[_trb_mtn_idx(i, kb, pcols)]) ** 2
                + (v[_trb_mtn_idx(i, kt, pcols)] - v[_trb_mtn_idx(i, kb, pcols)]) ** 2,
                dv2min,
            )

            ri = (
                2.0
                * gravit
                * (
                    t[_trb_mtn_idx(i, kt, pcols)] * exner[_trb_mtn_idx(i, kt, pcols)]
                    - t[_trb_mtn_idx(i, kb, pcols)] * exner[_trb_mtn_idx(i, kb, pcols)]
                )
                * (zm[_trb_mtn_idx(i, kt, pcols)] - zm[_trb_mtn_idx(i, kb, pcols)])
                / (
                    (
                        t[_trb_mtn_idx(i, kt, pcols)] * exner[_trb_mtn_idx(i, kt, pcols)]
                        + t[_trb_mtn_idx(i, kb, pcols)] * exner[_trb_mtn_idx(i, kb, pcols)]
                    )
                    * dv2
                )
            )

            stabfri = max(0.0, min(1.0, 1.0 - ri))
            cd = cd * stabfri

            rho = pmid[_trb_mtn_idx(i, pver, pcols)] / (rair * t[_trb_mtn_idx(i, pver, pcols)])
            vmag = sqrt(u[_trb_mtn_idx(i, pver, pcols)] ** 2 + v[_trb_mtn_idx(i, pver, pcols)] ** 2)
            ksrf[i - 1] = rho * cd * vmag * landfrac[i - 1]
            taux[i - 1] = -ksrf[i - 1] * u[_trb_mtn_idx(i, pver, pcols)]
            tauy[i - 1] = -ksrf[i - 1] * v[_trb_mtn_idx(i, pver, pcols)]


@inline
def _diff_solver_pcols_idx(i: int, k: int, pcols: int) -> int:
    """diffusion_solver work arrays declared as (pcols,pver+1)."""
    return (i - 1) + (k - 1) * pcols


@inline
def _diff_solver_ncol_idx(i: int, k: int, ncol: int) -> int:
    """diffusion_solver Coords1D copies and dpidz_sq declared as (ncol,*)."""
    return (i - 1) + (k - 1) * ncol


@export
def diffusion_solver_setup_codon(
    pcols: int,
    pver: int,
    ncol: int,
    ztodt: float,
    gravit: float,
    rair: float,
    t_p: cobj,
    rairi_p: cobj,
    p_ifc_p: cobj,
    p_mid_p: cobj,
    p_rdel_p: cobj,
    p_rdst_p: cobj,
    tint_p: cobj,
    rhoi_p: cobj,
    dpidz_sq_p: cobj,
    tmpi2_p: cobj,
    rrho_p: cobj,
    tmp1_p: cobj,
):
    t = Ptr[float](t_p)
    rairi = Ptr[float](rairi_p)
    p_ifc = Ptr[float](p_ifc_p)
    p_mid = Ptr[float](p_mid_p)
    p_rdel = Ptr[float](p_rdel_p)
    p_rdst = Ptr[float](p_rdst_p)
    tint = Ptr[float](tint_p)
    rhoi = Ptr[float](rhoi_p)
    dpidz_sq = Ptr[float](dpidz_sq_p)
    tmpi2 = Ptr[float](tmpi2_p)
    rrho = Ptr[float](rrho_p)
    tmp1 = Ptr[float](tmp1_p)

    for i in range(1, ncol + 1):
        tint[_diff_solver_pcols_idx(i, 1, pcols)] = t[_diff_solver_pcols_idx(i, 1, pcols)]
        rhoi[_diff_solver_pcols_idx(i, 1, pcols)] = (
            p_ifc[_diff_solver_ncol_idx(i, 1, ncol)]
            / (
                rairi[_diff_solver_pcols_idx(i, 1, pcols)]
                * tint[_diff_solver_pcols_idx(i, 1, pcols)]
            )
        )

    for k in range(2, pver + 1):
        for i in range(1, ncol + 1):
            tint[_diff_solver_pcols_idx(i, k, pcols)] = 0.5 * (
                t[_diff_solver_pcols_idx(i, k, pcols)]
                + t[_diff_solver_pcols_idx(i, k - 1, pcols)]
            )
            rhoi[_diff_solver_pcols_idx(i, k, pcols)] = (
                p_ifc[_diff_solver_ncol_idx(i, k, ncol)]
                / (
                    rairi[_diff_solver_pcols_idx(i, k, pcols)]
                    * tint[_diff_solver_pcols_idx(i, k, pcols)]
                )
            )

    for i in range(1, ncol + 1):
        tint[_diff_solver_pcols_idx(i, pver + 1, pcols)] = t[_diff_solver_pcols_idx(i, pver, pcols)]
        rhoi[_diff_solver_pcols_idx(i, pver + 1, pcols)] = (
            p_ifc[_diff_solver_ncol_idx(i, pver + 1, ncol)]
            / (rair * tint[_diff_solver_pcols_idx(i, pver + 1, pcols)])
        )

    for k in range(1, pver + 2):
        for i in range(1, ncol + 1):
            idx_n = _diff_solver_ncol_idx(i, k, ncol)
            dpidz_sq[idx_n] = gravit * rhoi[_diff_solver_pcols_idx(i, k, pcols)]

    for k in range(1, pver + 2):
        for i in range(1, ncol + 1):
            idx_n = _diff_solver_ncol_idx(i, k, ncol)
            dpidz_sq[idx_n] = dpidz_sq[idx_n] * dpidz_sq[idx_n]

    for k in range(2, pver + 1):
        for i in range(1, ncol + 1):
            tmpi2[_diff_solver_pcols_idx(i, k, pcols)] = (
                ztodt
                * dpidz_sq[_diff_solver_ncol_idx(i, k, ncol)]
                * p_rdst[_diff_solver_ncol_idx(i, k - 1, ncol)]
            )

    for i in range(1, ncol + 1):
        rrho[i - 1] = rair * t[_diff_solver_pcols_idx(i, pver, pcols)] / p_mid[
            _diff_solver_ncol_idx(i, pver, ncol)
        ]
        tmp1[i - 1] = ztodt * gravit * p_rdel[_diff_solver_ncol_idx(i, pver, ncol)]


@export
def cldfrc2m_aist_vector_codon(
    pcols: int,
    ncol: int,
    rhmaxi: float,
    qv_p: cobj,
    qi_p: cobj,
    qsat_p: cobj,
    esl_p: cobj,
    esi_p: cobj,
    rhmini_p: cobj,
    aist_p: cobj,
):
    qv = Ptr[float](qv_p)
    qi = Ptr[float](qi_p)
    qsat = Ptr[float](qsat_p)
    esl = Ptr[float](esl_p)
    esi = Ptr[float](esi_p)
    rhmini = Ptr[float](rhmini_p)
    aist_out = Ptr[float](aist_p)

    qist_min = 1.0e-7
    qist_max = 5.0e-3
    minice = 1.0e-12
    mincld = 1.0e-4

    for i0 in range(pcols):
        aist_out[i0] = 0.0

    for i in range(1, ncol + 1):
        i0 = i - 1
        rhi = (qv[i0] + qi[i0]) / qsat[i0] * (esl[i0] / esi[i0])
        rhdif = (rhi - rhmini[i0]) / (rhmaxi - rhmini[i0])
        aist = min(1.0, max(rhdif, 0.0) ** 2)

        if qi[i0] < minice:
            aist = 0.0
        else:
            aist = max(mincld, aist)

        if qi[i0] >= minice:
            icimr = qi[i0] / aist

            if icimr < qist_min:
                aist = max(0.0, min(1.0, qi[i0] / qist_min))

            if icimr > qist_max:
                aist = max(0.0, min(1.0, qi[i0] / qist_max))

        aist_out[i0] = max(0.0, min(aist, 0.999))


@inline
def _phys_prop_interp_index(n: int, x: Ptr[float], y: float) -> int:
    k = 1
    if y <= x[0]:
        k = 1
    elif y >= x[n - 1]:
        k = n - 1
    else:
        k = 1
        while y > x[k] and k < n:
            k += 1
    return k


@export
def phys_prop_exp_interpol_codon(n: int, x_p: cobj, f_p: cobj, y: float) -> float:
    x = Ptr[float](x_p)
    f = Ptr[float](f_p)

    k = _phys_prop_interp_index(n, x, y)
    k0 = k - 1
    k1 = k

    a = (log(f[k1] / f[k0])) / (x[k1] - x[k0])
    return f[k0] * exp(a * (y - x[k0]))


@export
def phys_prop_lin_interpol_codon(n: int, x_p: cobj, f_p: cobj, y: float) -> float:
    x = Ptr[float](x_p)
    f = Ptr[float](f_p)

    k = _phys_prop_interp_index(n, x, y)
    k0 = k - 1
    k1 = k

    a = (f[k1] - f[k0]) / (x[k1] - x[k0])
    return f[k0] + a * (y - x[k0])


@inline
def _modal_aer_opt_2d_idx(i: int, k: int, pcols: int) -> int:
    """modal_aer_opt arrays declared as (pcols,pver)."""
    return (i - 1) + (k - 1) * pcols


@inline
def _modal_aer_opt_cheb_idx(nc: int, i: int, k: int, ncoef: int, pcols: int) -> int:
    """modal_aer_opt Chebyshev arrays declared as (ncoef,pcols,pver)."""
    return (nc - 1) + (i - 1) * ncoef + (k - 1) * ncoef * pcols


@inline
def _modal_aer_opt_table_idx(k: int, i: int, j: int, km: int, im: int) -> int:
    """binterp table declared as (km,im,jm)."""
    return (k - 1) + (i - 1) * km + (j - 1) * km * im


@inline
def _modal_aer_opt_out_idx(i: int, k: int, pcols: int) -> int:
    """binterp out declared as (pcols,km)."""
    return (i - 1) + (k - 1) * pcols


@export
def modal_aer_opt_size_parameters_codon(
    pcols: int,
    pver: int,
    top_lev: int,
    ncol: int,
    ncoef: int,
    sigma_logr_aer: float,
    xrmin: float,
    xrmax: float,
    dgnumwet_p: cobj,
    radsurf_p: cobj,
    logradsurf_p: cobj,
    cheb_p: cobj,
):
    dgnumwet = Ptr[float](dgnumwet_p)
    radsurf = Ptr[float](radsurf_p)
    logradsurf = Ptr[float](logradsurf_p)
    cheb = Ptr[float](cheb_p)

    alnsg_amode = log(sigma_logr_aer)
    explnsigma = exp(2.0 * alnsg_amode * alnsg_amode)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _modal_aer_opt_2d_idx(i, k, pcols)
            radsurf[idx2] = 0.5 * dgnumwet[idx2] * explnsigma
            logradsurf[idx2] = log(radsurf[idx2])
            xrad = max(logradsurf[idx2], xrmin)
            xrad = min(xrad, xrmax)
            xrad = (2.0 * xrad - xrmax - xrmin) / (xrmax - xrmin)
            cheb[_modal_aer_opt_cheb_idx(1, i, k, ncoef, pcols)] = 1.0
            cheb[_modal_aer_opt_cheb_idx(2, i, k, ncoef, pcols)] = xrad
            for nc in range(3, ncoef + 1):
                cheb[_modal_aer_opt_cheb_idx(nc, i, k, ncoef, pcols)] = (
                    2.0
                    * xrad
                    * cheb[_modal_aer_opt_cheb_idx(nc - 1, i, k, ncoef, pcols)]
                    - cheb[_modal_aer_opt_cheb_idx(nc - 2, i, k, ncoef, pcols)]
                )


@export
def modal_aer_opt_lw_size_parameters_codon(
    pcols: int,
    pver: int,
    top_lev: int,
    ncol: int,
    ncoef: int,
    sigma_logr_aer: float,
    xrmin: float,
    xrmax: float,
    dgnumwet_p: cobj,
    cheby_p: cobj,
):
    dgnumwet = Ptr[float](dgnumwet_p)
    cheby = Ptr[float](cheby_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            alnsg_amode = log(sigma_logr_aer)
            xrad = (
                log(0.5 * dgnumwet[_modal_aer_opt_2d_idx(i, k, pcols)])
                + 2.0 * alnsg_amode * alnsg_amode
            )
            xrad = max(xrad, xrmin)
            xrad = min(xrad, xrmax)
            xrad = (2.0 * xrad - xrmax - xrmin) / (xrmax - xrmin)
            cheby[_modal_aer_opt_cheb_idx(1, i, k, ncoef, pcols)] = 1.0
            cheby[_modal_aer_opt_cheb_idx(2, i, k, ncoef, pcols)] = xrad
            for nc in range(3, ncoef + 1):
                cheby[_modal_aer_opt_cheb_idx(nc, i, k, ncoef, pcols)] = (
                    2.0
                    * xrad
                    * cheby[_modal_aer_opt_cheb_idx(nc - 1, i, k, ncoef, pcols)]
                    - cheby[_modal_aer_opt_cheb_idx(nc - 2, i, k, ncoef, pcols)]
                )


@export
def modal_aer_opt_binterp_codon(
    pcols: int,
    ncol: int,
    km: int,
    im: int,
    jm: int,
    table_p: cobj,
    x_p: cobj,
    y_p: cobj,
    xtab_p: cobj,
    ytab_p: cobj,
    ix_p: cobj,
    jy_p: cobj,
    t_p: cobj,
    u_p: cobj,
    out_p: cobj,
):
    table = Ptr[float](table_p)
    x = Ptr[float](x_p)
    y = Ptr[float](y_p)
    xtab = Ptr[float](xtab_p)
    ytab = Ptr[float](ytab_p)
    ix = Ptr[i32](ix_p)
    jy = Ptr[i32](jy_p)
    t = Ptr[float](t_p)
    u = Ptr[float](u_p)
    out = Ptr[float](out_p)

    if int(ix[0]) <= 0:
        if im > 1:
            for ic in range(1, ncol + 1):
                found_i = im + 1
                for ii in range(1, im + 1):
                    if x[ic - 1] < xtab[ii - 1]:
                        found_i = ii
                        break
                ix[ic - 1] = i32(max(found_i - 1, 1))
                ip1 = min(int(ix[ic - 1]) + 1, im)
                dx = xtab[ip1 - 1] - xtab[int(ix[ic - 1]) - 1]
                if abs(dx) > 1.0e-20:
                    t[ic - 1] = (x[ic - 1] - xtab[int(ix[ic - 1]) - 1]) / dx
                else:
                    t[ic - 1] = 0.0
        else:
            for ic in range(1, ncol + 1):
                ix[ic - 1] = i32(1)
                t[ic - 1] = 0.0

        if jm > 1:
            for ic in range(1, ncol + 1):
                found_j = jm + 1
                for jj in range(1, jm + 1):
                    if y[ic - 1] < ytab[jj - 1]:
                        found_j = jj
                        break
                jy[ic - 1] = i32(max(found_j - 1, 1))
                jp1 = min(int(jy[ic - 1]) + 1, jm)
                dy = ytab[jp1 - 1] - ytab[int(jy[ic - 1]) - 1]
                if abs(dy) > 1.0e-20:
                    u[ic - 1] = (y[ic - 1] - ytab[int(jy[ic - 1]) - 1]) / dy
                else:
                    u[ic - 1] = 0.0
        else:
            for ic in range(1, ncol + 1):
                jy[ic - 1] = i32(1)
                u[ic - 1] = 0.0

    for ic in range(1, ncol + 1):
        ic0 = ic - 1
        tu = t[ic0] * u[ic0]
        tuc = t[ic0] - tu
        tcuc = 1.0 - tuc - u[ic0]
        tcu = u[ic0] - tu
        jp1 = min(int(jy[ic0]) + 1, jm)
        ip1 = min(int(ix[ic0]) + 1, im)
        for k in range(1, km + 1):
            value = (
                tcuc
                * table[_modal_aer_opt_table_idx(k, int(ix[ic0]), int(jy[ic0]), km, im)]
                + tuc * table[_modal_aer_opt_table_idx(k, ip1, int(jy[ic0]), km, im)]
            )
            value = value + tu * table[_modal_aer_opt_table_idx(k, ip1, jp1, km, im)]
            value = value + tcu * table[_modal_aer_opt_table_idx(k, int(ix[ic0]), jp1, km, im)]
            out[_modal_aer_opt_out_idx(ic, k, pcols)] = value


@export
def ndrop_mode_props_finalize_codon(
    nmode: int,
    pi: float,
    sigmag_p: cobj,
    dgnumlo_p: cobj,
    dgnumhi_p: cobj,
    alogsig_p: cobj,
    exp45logsig_p: cobj,
    f1_p: cobj,
    f2_p: cobj,
    voltonumblo_p: cobj,
    voltonumbhi_p: cobj,
):
    sigmag = Ptr[float](sigmag_p)
    dgnumlo = Ptr[float](dgnumlo_p)
    dgnumhi = Ptr[float](dgnumhi_p)
    alogsig = Ptr[float](alogsig_p)
    exp45logsig = Ptr[float](exp45logsig_p)
    f1 = Ptr[float](f1_p)
    f2 = Ptr[float](f2_p)
    voltonumblo = Ptr[float](voltonumblo_p)
    voltonumbhi = Ptr[float](voltonumbhi_p)

    for m in range(1, nmode + 1):
        m0 = m - 1
        alogsig[m0] = log(sigmag[m0])
        exp45logsig[m0] = exp(4.5 * alogsig[m0] * alogsig[m0])
        f1[m0] = 0.5 * exp(2.5 * alogsig[m0] * alogsig[m0])
        f2[m0] = 1.0 + 0.25 * alogsig[m0]

        voltonumblo[m0] = 1.0 / (
            (pi / 6.0) * (dgnumlo[m0] ** 3.0) * exp(4.5 * alogsig[m0] ** 2.0)
        )
        voltonumbhi[m0] = 1.0 / (
            (pi / 6.0) * (dgnumhi[m0] ** 3.0) * exp(4.5 * alogsig[m0] ** 2.0)
        )


@export
def ghg_data_mw_ratios_codon(
    mwdry: float,
    mwn2o: float,
    mwch4: float,
    mwf11: float,
    mwf12: float,
    mwco2: float,
    rmwn2o_p: cobj,
    rmwch4_p: cobj,
    rmwf11_p: cobj,
    rmwf12_p: cobj,
    rmwco2_p: cobj,
):
    rmwn2o = Ptr[float](rmwn2o_p)
    rmwch4 = Ptr[float](rmwch4_p)
    rmwf11 = Ptr[float](rmwf11_p)
    rmwf12 = Ptr[float](rmwf12_p)
    rmwco2 = Ptr[float](rmwco2_p)

    rmwn2o[0] = mwn2o / mwdry
    rmwch4[0] = mwch4 / mwdry
    rmwf11[0] = mwf11 / mwdry
    rmwf12[0] = mwf12 / mwdry
    rmwco2[0] = mwco2 / mwdry


@export
def rad_cnst_out_mass_cb_codon(
    ncol: int,
    pcols: int,
    pver: int,
    rga: float,
    mmr_p: cobj,
    pdeldry_p: cobj,
    mass_p: cobj,
    cb_p: cobj,
):
    mmr = Ptr[float](mmr_p)
    pdeldry = Ptr[float](pdeldry_p)
    mass = Ptr[float](mass_p)
    cb = Ptr[float](cb_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _field2_idx(i, k, pcols)
            mass[idx] = (mmr[idx] * pdeldry[idx]) * rga

    for i in range(1, ncol + 1):
        total = 0.0
        for k in range(1, pver + 1):
            total = total + mass[_field2_idx(i, k, pcols)]
        cb[i - 1] = total


@export
def phys_control_deepconv_pbl_codon(eddy_diag_tke: int, shallow_uw: int) -> int:
    if eddy_diag_tke != 0 or shallow_uw != 0:
        return 1
    return 0


@export
def phys_control_do_flux_avg_codon(srf_flux_avg: int) -> int:
    if srf_flux_avg == 1:
        return 1
    return 0


@export
def constituents_rgas_codon(r_universal: float, mwc: float) -> float:
    return r_universal * mwc


@export
def constituents_cv_codon(cpc: float, rgas: float) -> float:
    return cpc - rgas


@inline
def _aer_rad_sw_idx(i: int, k0: int, band: int, pcols: int, pverp: int) -> int:
    """tau/tau_w/tau_w_g/tau_w_f declared as (pcols,0:pver,nswbands)."""
    return (i - 1) + k0 * pcols + (band - 1) * pcols * pverp


@export
def aer_rad_props_sw_setup_codon(
    ncol: int,
    pcols: int,
    pver: int,
    nswbands: int,
    nrh: int,
    rga: float,
    pdeldry_p: cobj,
    qv_p: cobj,
    qs_p: cobj,
    mmr_to_mass_p: cobj,
    krh_p: cobj,
    wrh_p: cobj,
    tau_p: cobj,
    tau_w_p: cobj,
    tau_w_g_p: cobj,
    tau_w_f_p: cobj,
):
    pdeldry = Ptr[float](pdeldry_p)
    qv = Ptr[float](qv_p)
    qs = Ptr[float](qs_p)
    mmr_to_mass = Ptr[float](mmr_to_mass_p)
    krh = Ptr[i32](krh_p)
    wrh = Ptr[float](wrh_p)
    tau = Ptr[float](tau_p)
    tau_w = Ptr[float](tau_w_p)
    tau_w_g = Ptr[float](tau_w_g_p)
    tau_w_f = Ptr[float](tau_w_f_p)

    pverp = pver + 1
    for band in range(1, nswbands + 1):
        for k0 in range(0, pver + 1):
            for i in range(1, pcols + 1):
                idx3 = _aer_rad_sw_idx(i, k0, band, pcols, pverp)
                tau[idx3] = -100.0
                tau_w[idx3] = -100.0
                tau_w_g[idx3] = -100.0
                tau_w_f[idx3] = -100.0
            for i in range(1, ncol + 1):
                idx3 = _aer_rad_sw_idx(i, k0, band, pcols, pverp)
                tau[idx3] = 0.0
                tau_w[idx3] = 0.0
                tau_w_g[idx3] = 0.0
                tau_w_f[idx3] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _field2_idx(i, k, pcols)
            mmr_to_mass[idx2] = rga * pdeldry[idx2]
            rh = qv[idx2] / qs[idx2]
            rhtrunc = min(rh, 1.0)
            krh_value = min(int(floor(rhtrunc * float(nrh))) + 1, nrh - 1)
            krh[idx2] = i32(krh_value)
            wrh[idx2] = rhtrunc * float(nrh) - float(krh_value)


@export
def aer_rad_props_lw_setup_codon(
    ncol: int,
    pcols: int,
    pver: int,
    nrh: int,
    rga: float,
    pdeldry_p: cobj,
    qv_p: cobj,
    qs_p: cobj,
    mmr_to_mass_p: cobj,
    krh_p: cobj,
    wrh_p: cobj,
):
    pdeldry = Ptr[float](pdeldry_p)
    qv = Ptr[float](qv_p)
    qs = Ptr[float](qs_p)
    mmr_to_mass = Ptr[float](mmr_to_mass_p)
    krh = Ptr[i32](krh_p)
    wrh = Ptr[float](wrh_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _field2_idx(i, k, pcols)
            mmr_to_mass[idx2] = rga * pdeldry[idx2]
            rh = qv[idx2] / qs[idx2]
            rhtrunc = min(rh, 1.0)
            krh_value = min(int(floor(rhtrunc * float(nrh))) + 1, nrh - 1)
            krh[idx2] = i32(krh_value)
            wrh[idx2] = rhtrunc * float(nrh) - float(krh_value)


@export
def micro_mg_utils_init_scalars_codon(
    rh2o: float,
    cpair: float,
    tmelt_in: float,
    latvap: float,
    latice: float,
    rv_p: cobj,
    cpp_p: cobj,
    tmelt_p: cobj,
    xxlv_p: cobj,
    xlf_p: cobj,
    xxls_p: cobj,
):
    rv = Ptr[float](rv_p)
    cpp = Ptr[float](cpp_p)
    tmelt = Ptr[float](tmelt_p)
    xxlv = Ptr[float](xxlv_p)
    xlf = Ptr[float](xlf_p)
    xxls = Ptr[float](xxls_p)

    rv[0] = rh2o
    cpp[0] = cpair
    tmelt[0] = tmelt_in
    xxlv[0] = latvap
    xlf[0] = latice
    xxls[0] = xxlv[0] + xlf[0]


@export
def radiation_data_flag_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def cosp_set_values_basic_codon(
    nlr: int,
    use_vgrid: int,
    csat_vgrid: int,
    ncolumns: int,
    nradsteps: int,
    nht_current: int,
    nht_p: cobj,
    nscol_p: cobj,
    nradsteps_p: cobj,
    zstep_p: cobj,
):
    nht = Ptr[int](nht_p)
    nscol = Ptr[int](nscol_p)
    nradsteps_out = Ptr[int](nradsteps_p)
    zstep = Ptr[float](zstep_p)

    nht_value = nht_current
    zstep_value = 0.0
    if use_vgrid != 0:
        if csat_vgrid != 0:
            nht_value = 40
            zstep_value = 480.0
        else:
            nht_value = nlr
            zstep_value = 20000.0 / float(nlr)

    nht[0] = nht_value
    nscol[0] = ncolumns
    nradsteps_out[0] = nradsteps
    zstep[0] = zstep_value


@export
def carma_flags_bool_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def carma_flags_touch_codon() -> int:
    return 0


@export
def tidal_diag_int_codon(value: int, force_one: int) -> int:
    if force_one != 0:
        return 1
    return value


@export
def co2_cycle_flag_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def carma_intr_false_codon() -> int:
    return 0


@export
def carma_intr_touch_codon() -> int:
    return 0


@export
def subcol_touch_codon() -> int:
    return 0


@export
def cldwat_param_codon(value: float) -> float:
    return value


@export
def hkconv_param_codon(value: float) -> float:
    return value


@export
def cam3_aero_data_flag_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def cam3_ozone_data_flag_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def phys_debug_value_codon(value: float) -> float:
    return value


@export
def phys_debug_has_location_codon(lat_set: int, lon_set: int) -> int:
    if lat_set != 0 and lon_set != 0:
        return 1
    return 0


@export
def unicon_cam_flag_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def unicon_cam_int_codon(value: int) -> int:
    return value


@export
def iondrag_touch_codon() -> int:
    return 0


@export
def cld_sediment_param_codon(value: float) -> float:
    return value


@export
def tsinti_param_codon(value: float) -> float:
    return value


@export
def clubb_intr_flag_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def clubb_intr_touch_codon() -> int:
    return 0


@export
def radae_ntoplw_codon(pref_mid_p: cobj, nlev: int) -> int:
    pref_mid = Ptr[float](pref_mid_p)
    ntoplw = 1
    if pref_mid[0] < 0.1:
        for k in range(nlev):
            if pref_mid[k] < 1.0:
                ntoplw = k + 1
    else:
        ntoplw = 1
    return ntoplw


@export
def hirsbt_flag_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def hirsbt_freq_codon(freq: int, dtime: int) -> int:
    if freq < 0:
        value = (-freq * 3600.0) / float(dtime)
        return int(value + 0.5)
    return freq


@export
def hetfrz_classnuc_cam_flag_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def cpslec_codon(
    ncol: int,
    pmid_p: cobj,
    phis_p: cobj,
    ps_p: cobj,
    t_p: cobj,
    psl_p: cobj,
    gravit: float,
    rair: float,
    pcols: int,
    pver: int,
):
    pmid = Ptr[float](pmid_p)
    phis = Ptr[float](phis_p)
    ps = Ptr[float](ps_p)
    temp = Ptr[float](t_p)
    psl = Ptr[float](psl_p)

    xlapse = 6.5e-3
    alpha = rair * xlapse / gravit
    kbot_offset = (pver - 1) * pcols

    for i in range(ncol):
        if abs(phis[i] / gravit) < 1.0e-4:
            psl[i] = ps[i]
        else:
            tstar = temp[i + kbot_offset] * (1.0 + alpha * (ps[i] / pmid[i + kbot_offset] - 1.0))
            tt0 = tstar + xlapse * phis[i] / gravit

            alph = 0.0
            if tstar <= 290.5 and tt0 > 290.5:
                alph = rair / phis[i] * (290.5 - tstar)
            elif tstar > 290.5 and tt0 > 290.5:
                alph = 0.0
                tstar = 0.5 * (290.5 + tstar)
            else:
                alph = alpha
                if tstar < 255.0:
                    tstar = 0.5 * (255.0 + tstar)

            beta = phis[i] / (rair * tstar)
            ab = alph * beta
            psl[i] = ps[i] * exp(beta * (1.0 - alph * beta / 2.0 + (ab**2.0) / 3.0))


@export
def sslt_rebin_has_four_codon(i1: int, i2: int, i3: int, i4: int) -> int:
    if i1 > 0 and i2 > 0 and i3 > 0 and i4 > 0:
        return 1
    return 0


@export
def constituent_burden_flag_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def phys_gmean_normalize_codon(arr_p: cobj, nflds: int, pi_value: float):
    arr = Ptr[float](arr_p)
    denom = 4.0 * pi_value
    for i in range(nflds):
        arr[i] = arr[i] / denom


@export
def pbl_utils_value_codon(value: float) -> float:
    return value


@export
def wv_sat_methods_value_codon(value: float) -> float:
    return value


@export
def wv_sat_methods_omeps_codon(epsilo: float) -> float:
    return 1.0 - epsilo


@export
def wv_saturation_value_codon(value: float) -> float:
    return value
