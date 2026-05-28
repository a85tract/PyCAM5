import modal_aer_opt_codon as _modal_aer_opt
import ndrop_codon as _ndrop
import phys_grid_codon as _phys_grid
from C import micro_mg_utils_gamma_native_cb(float) -> float
from C import micro_mg_utils_shape_coef_native_cb(float, float) -> float
from C import micro_mg_utils_rising_factorial_native_cb(float, float) -> float
from C import micro_mg_utils_liq_pgam_native_cb(float, float) -> float
from C import micro_mg_utils_liq_shape_coef_native_cb(float, float, float) -> float
from C import micro_mg_utils_basic_lam_native_cb(float, float, float, float) -> float
from C import micro_mg_utils_basic_nic_native_cb(float, float, float, float) -> float
from C import micro_mg_utils_avg_diameter_native_cb(float, float, float, float) -> float
from C import goffgratch_svp_ice_native_cb(float) -> float
from C import goffgratch_svp_water_native_cb(float) -> float
from C import bolton_svp_water_native_cb(float) -> float
from C import murphykoop_svp_ice_native_cb(float) -> float
from C import murphykoop_svp_water_native_cb(float) -> float
from C import oldgoffgratch_svp_ice_native_cb(float) -> float
from C import oldgoffgratch_svp_water_native_cb(float) -> float
from C import zm_entropy_expr_native_cb(float, float, float, float, float, float, float, float, float, float, float, float) -> float
from math import acos, cos, exp, floor, log, sin, sqrt

@export
def physpkg_orch_stage_codon(stage: int, flag1: int, flag2: int, flag3: int) -> int:
    mask = stage
    if flag1 != 0:
        mask |= 16
    if flag2 != 0:
        mask |= 32
    if flag3 != 0:
        mask |= 64
    return mask


@export
def phys_register_codon(stage: int, flag1: int, flag2: int, flag3: int) -> int:
    return physpkg_orch_stage_codon(stage, flag1, flag2, flag3)


@export
def phys_init_codon(stage: int, flag1: int, flag2: int, flag3: int) -> int:
    return physpkg_orch_stage_codon(stage, flag1, flag2, flag3)


@export
def phys_run1_codon(stage: int, flag1: int, flag2: int, flag3: int) -> int:
    return physpkg_orch_stage_codon(stage, flag1, flag2, flag3)


@export
def phys_run2_codon(stage: int, flag1: int, flag2: int, flag3: int) -> int:
    return physpkg_orch_stage_codon(stage, flag1, flag2, flag3)


@export
def phys_inidat_shell_mask_codon(aqua_planet: int, unstructured: int, chunk_count: int, dyn_time_lvls: int) -> int:
    mask = 1
    if aqua_planet != 0:
        mask |= 2
    if unstructured != 0:
        mask |= 4
    if chunk_count > 0:
        mask |= 8
    if dyn_time_lvls > 1:
        mask |= 16
    return mask


@export
def phys_inidat_codon(aqua_planet: int, unstructured: int, chunk_count: int, dyn_time_lvls: int) -> int:
    return phys_inidat_shell_mask_codon(aqua_planet, unstructured, chunk_count, dyn_time_lvls)


@export
def tphys_shell_mask_codon(stage: int, ncol: int, flag1: int, flag2: int, flag3: int, flag4: int, flag5: int) -> int:
    mask = stage
    if ncol > 0:
        mask |= 16
    if flag1 != 0:
        mask |= 32
    if flag2 != 0:
        mask |= 64
    if flag3 != 0:
        mask |= 128
    if flag4 != 0:
        mask |= 256
    if flag5 != 0:
        mask |= 512
    return mask


@export
def tphysbc_codon(stage: int, ncol: int, flag1: int, flag2: int, flag3: int, flag4: int, flag5: int) -> int:
    return tphys_shell_mask_codon(stage, ncol, flag1, flag2, flag3, flag4, flag5)


@export
def tphysac_codon(stage: int, ncol: int, flag1: int, flag2: int, flag3: int, flag4: int, flag5: int) -> int:
    return tphys_shell_mask_codon(stage, ncol, flag1, flag2, flag3, flag4, flag5)


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


@export
def phys_timestep_init_codon(
    cam3_aero_on: int,
    cam3_ozone_on: int,
    do_waccm_ions: int,
    branch_mask_p: cobj,
):
    phys_timestep_init_select_branches_codon(cam3_aero_on, cam3_ozone_on, do_waccm_ions, branch_mask_p)


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


@export
def final_cam_cleanup_touch_codon(stage: int) -> int:
    return stage


@export
def init_restart_physics_codon(stage: int) -> int:
    return final_cam_cleanup_touch_codon(stage)


@export
def write_restart_physics_codon(stage: int) -> int:
    return final_cam_cleanup_touch_codon(stage)


@export
def tropopause_init_codon(stage: int) -> int:
    return final_cam_cleanup_touch_codon(stage)


@export
def tropopause_read_file_codon(stage: int) -> int:
    return final_cam_cleanup_touch_codon(stage)


@export
def tropopause_climate_codon(stage: int) -> int:
    return final_cam_cleanup_touch_codon(stage)


@export
def tropopause_find_codon(stage: int) -> int:
    return final_cam_cleanup_touch_codon(stage)


@export
def tropopause_findusing_codon(stage: int) -> int:
    return final_cam_cleanup_touch_codon(stage)


@export
def tropopause_output_codon(stage: int) -> int:
    return final_cam_cleanup_touch_codon(stage)


@export
def tropopause_readnl_codon(stage: int) -> int:
    return final_cam_cleanup_touch_codon(stage)


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


@export
def tropopause_interpolateT_codon(
    i: int,
    lev: int,
    trop_p: float,
    pcols: int,
    pver: int,
    state_t_p: cobj,
    state_pmid_p: cobj,
) -> float:
    return _tropopause_interp_t(
        i,
        lev,
        trop_p,
        pcols,
        pver,
        Ptr[float](state_t_p),
        Ptr[float](state_pmid_p),
    )


@export
def tropopause_interpolatet_codon(
    i: int,
    lev: int,
    trop_p: float,
    pcols: int,
    pver: int,
    state_t_p: cobj,
    state_pmid_p: cobj,
) -> float:
    return tropopause_interpolateT_codon(i, lev, trop_p, pcols, pver, state_t_p, state_pmid_p)


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


@export
def tropopause_interpolateZ_codon(
    i: int,
    lev: int,
    trop_p: float,
    pcols: int,
    state_zm_p: cobj,
    state_zi_p: cobj,
    state_pmid_p: cobj,
    state_pint_p: cobj,
) -> float:
    return _tropopause_interp_z(
        i,
        lev,
        trop_p,
        pcols,
        Ptr[float](state_zm_p),
        Ptr[float](state_zi_p),
        Ptr[float](state_pmid_p),
        Ptr[float](state_pint_p),
    )


@export
def tropopause_interpolatez_codon(
    i: int,
    lev: int,
    trop_p: float,
    pcols: int,
    state_zm_p: cobj,
    state_zi_p: cobj,
    state_pmid_p: cobj,
    state_pint_p: cobj,
) -> float:
    return tropopause_interpolateZ_codon(
        i,
        lev,
        trop_p,
        pcols,
        state_zm_p,
        state_zi_p,
        state_pmid_p,
        state_pint_p,
    )


@export
def tropopause_climate_find_codon(
    ncol: int,
    pcols: int,
    pver: int,
    chunk_pos: int,
    chunk_count: int,
    last_month: int,
    next_month: int,
    notfound: int,
    dels: float,
    pint_p: cobj,
    tropp_p_loc_p: cobj,
    trop_lev_p: cobj,
    trop_p_p: cobj,
    updated_p: cobj,
):
    pint = Ptr[float](pint_p)
    tropp_p_loc = Ptr[float](tropp_p_loc_p)
    trop_lev = Ptr[int](trop_lev_p)
    trop_p = Ptr[float](trop_p_p)
    updated = Ptr[int](updated_p)

    for i in range(1, pcols + 1):
        updated[i - 1] = 0

    for i in range(1, ncol + 1):
        if trop_lev[i - 1] == notfound:
            last_idx = (
                (i - 1)
                + (chunk_pos - 1) * pcols
                + (last_month - 1) * pcols * chunk_count
            )
            next_idx = (
                (i - 1)
                + (chunk_pos - 1) * pcols
                + (next_month - 1) * pcols * chunk_count
            )
            tp = tropp_p_loc[last_idx] + dels * (
                tropp_p_loc[next_idx] - tropp_p_loc[last_idx]
            )
            trop_p[i - 1] = tp
            for k in range(pver, 1, -1):
                if tp >= pint[_field2_idx(i, k, pcols)]:
                    trop_lev[i - 1] = k
                    updated[i - 1] = 1
                    break


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
def tropopause_twmo_pressure_profile_codon(
    level: int,
    cnst_kap: float,
    cnst_ka1: float,
    cnst_faktor: float,
    plimu: float,
    pliml: float,
    gam: float,
    t_p: cobj,
    p_p: cobj,
) -> float:
    t = Ptr[float](t_p)
    p = Ptr[float](p_p)
    deltaz = 2000.0

    trp = -99.0
    pmk = 0.5 * (p[level - 2] ** cnst_kap + p[level - 1] ** cnst_kap)
    pm = pmk ** (1.0 / cnst_kap)
    a = (t[level - 2] - t[level - 1]) / (p[level - 2] ** cnst_kap - p[level - 1] ** cnst_kap)
    b = t[level - 1] - (a * p[level - 1] ** cnst_kap)
    tm = a * pmk + b
    dtdp = a * cnst_kap * (pm ** cnst_ka1)
    dtdz = cnst_faktor * dtdp * pm / tm

    for j in range(level - 1, 1, -1):
        pm0 = pm
        pmk0 = pmk
        dtdz0 = dtdz

        pmk = 0.5 * (p[j - 2] ** cnst_kap + p[j - 1] ** cnst_kap)
        pm = pmk ** (1.0 / cnst_kap)
        a = (t[j - 2] - t[j - 1]) / (p[j - 2] ** cnst_kap - p[j - 1] ** cnst_kap)
        b = t[j - 1] - (a * p[j - 1] ** cnst_kap)
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
            pmk2 = 0.5 * (p[jj - 2] ** cnst_kap + p[jj - 1] ** cnst_kap)
            pm2 = pmk2 ** (1.0 / cnst_kap)
            if pm2 > ptph:
                continue
            if pm2 < p2km:
                break

            a2 = t[jj - 2] - t[jj - 1]
            a2 = a2 / (p[jj - 2] ** cnst_kap - p[jj - 1] ** cnst_kap)
            b2 = t[jj - 1] - (a2 * p[jj - 1] ** cnst_kap)
            tm2 = a2 * pmk2 + b2
            dtdp2 = a2 * cnst_kap * (pm2 ** (cnst_kap - 1.0))
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
def _physics_types_fill_1d(n: int, value: float, arr: Ptr[float]):
    for i in range(0, n):
        arr[i] = value


@inline
def _physics_types_fill_2d(psetcols: int, pver: int, value: float, arr: Ptr[float]):
    for k in range(1, pver + 1):
        for i in range(1, psetcols + 1):
            arr[_field2_idx(i, k, psetcols)] = value


@inline
def _physics_types_fill_3d(psetcols: int, pver: int, pcnst: int, value: float, arr: Ptr[float]):
    for m in range(1, pcnst + 1):
        for k in range(1, pver + 1):
            for i in range(1, psetcols + 1):
                arr[_field3_idx(i, k, m, psetcols, pver)] = value


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
def _physics_types_copy_2d_ncol(ncol: int, ld1: int, nlev: int, src: Ptr[float], dst: Ptr[float]):
    for k in range(1, nlev + 1):
        for i in range(1, ncol + 1):
            dst[_field2_idx(i, k, ld1)] = src[_field2_idx(i, k, ld1)]


@inline
def _physics_types_copy_3d(psetcols: int, pver: int, pcnst: int, src: Ptr[float], dst: Ptr[float]):
    for m in range(1, pcnst + 1):
        for k in range(1, pver + 1):
            for i in range(1, psetcols + 1):
                dst[_field3_idx(i, k, m, psetcols, pver)] = src[_field3_idx(i, k, m, psetcols, pver)]


@inline
def _physics_types_copy_3d_ncol(ncol: int, ld1: int, nlev: int, n3: int, src: Ptr[float], dst: Ptr[float]):
    for m in range(1, n3 + 1):
        for k in range(1, nlev + 1):
            for i in range(1, ncol + 1):
                dst[_field3_idx(i, k, m, ld1, nlev)] = src[_field3_idx(i, k, m, ld1, nlev)]


@export
def physics_dme_adjust_active_codon(is_lr: int) -> int:
    if is_lr != 0:
        return 1
    return 0


@export
def physics_dme_adjust_codon(is_lr: int) -> int:
    return physics_dme_adjust_active_codon(is_lr)


@export
def physics_types_touch_codon(stage: int) -> int:
    return stage


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
def physics_state_alloc_codon(
    psetcols: int,
    pver: int,
    pcnst: int,
    value: float,
    lat_p: cobj,
    lon_p: cobj,
    ulat_p: cobj,
    ulon_p: cobj,
    ps_p: cobj,
    psdry_p: cobj,
    phis_p: cobj,
    t_p: cobj,
    u_p: cobj,
    v_p: cobj,
    s_p: cobj,
    omega_p: cobj,
    pmid_p: cobj,
    pmiddry_p: cobj,
    pdel_p: cobj,
    pdeldry_p: cobj,
    rpdel_p: cobj,
    rpdeldry_p: cobj,
    lnpmid_p: cobj,
    lnpmiddry_p: cobj,
    exner_p: cobj,
    zm_p: cobj,
    q_p: cobj,
    pint_p: cobj,
    pintdry_p: cobj,
    lnpint_p: cobj,
    lnpintdry_p: cobj,
    zi_p: cobj,
    te_ini_p: cobj,
    te_cur_p: cobj,
    tw_ini_p: cobj,
    tw_cur_p: cobj,
):
    _physics_types_fill_1d(psetcols, value, Ptr[float](lat_p))
    _physics_types_fill_1d(psetcols, value, Ptr[float](lon_p))
    _physics_types_fill_1d(psetcols, value, Ptr[float](ulat_p))
    _physics_types_fill_1d(psetcols, value, Ptr[float](ulon_p))
    _physics_types_fill_1d(psetcols, value, Ptr[float](ps_p))
    _physics_types_fill_1d(psetcols, value, Ptr[float](psdry_p))
    _physics_types_fill_1d(psetcols, value, Ptr[float](phis_p))
    _physics_types_fill_2d(psetcols, pver, value, Ptr[float](t_p))
    _physics_types_fill_2d(psetcols, pver, value, Ptr[float](u_p))
    _physics_types_fill_2d(psetcols, pver, value, Ptr[float](v_p))
    _physics_types_fill_2d(psetcols, pver, value, Ptr[float](s_p))
    _physics_types_fill_2d(psetcols, pver, value, Ptr[float](omega_p))
    _physics_types_fill_2d(psetcols, pver, value, Ptr[float](pmid_p))
    _physics_types_fill_2d(psetcols, pver, value, Ptr[float](pmiddry_p))
    _physics_types_fill_2d(psetcols, pver, value, Ptr[float](pdel_p))
    _physics_types_fill_2d(psetcols, pver, value, Ptr[float](pdeldry_p))
    _physics_types_fill_2d(psetcols, pver, value, Ptr[float](rpdel_p))
    _physics_types_fill_2d(psetcols, pver, value, Ptr[float](rpdeldry_p))
    _physics_types_fill_2d(psetcols, pver, value, Ptr[float](lnpmid_p))
    _physics_types_fill_2d(psetcols, pver, value, Ptr[float](lnpmiddry_p))
    _physics_types_fill_2d(psetcols, pver, value, Ptr[float](exner_p))
    _physics_types_fill_2d(psetcols, pver, value, Ptr[float](zm_p))
    _physics_types_fill_3d(psetcols, pver, pcnst, value, Ptr[float](q_p))

    _physics_types_fill_2d(psetcols, pver + 1, value, Ptr[float](pint_p))
    _physics_types_fill_2d(psetcols, pver + 1, value, Ptr[float](pintdry_p))
    _physics_types_fill_2d(psetcols, pver + 1, value, Ptr[float](lnpint_p))
    _physics_types_fill_2d(psetcols, pver + 1, value, Ptr[float](lnpintdry_p))
    _physics_types_fill_2d(psetcols, pver + 1, value, Ptr[float](zi_p))

    _physics_types_fill_1d(psetcols, value, Ptr[float](te_ini_p))
    _physics_types_fill_1d(psetcols, value, Ptr[float](te_cur_p))
    _physics_types_fill_1d(psetcols, value, Ptr[float](tw_ini_p))
    _physics_types_fill_1d(psetcols, value, Ptr[float](tw_cur_p))


@export
def physics_tend_alloc_codon(
    psetcols: int,
    pver: int,
    value: float,
    dtdt_p: cobj,
    dudt_p: cobj,
    dvdt_p: cobj,
    flx_net_p: cobj,
    te_tnd_p: cobj,
    tw_tnd_p: cobj,
):
    _physics_types_fill_2d(psetcols, pver, value, Ptr[float](dtdt_p))
    _physics_types_fill_2d(psetcols, pver, value, Ptr[float](dudt_p))
    _physics_types_fill_2d(psetcols, pver, value, Ptr[float](dvdt_p))
    _physics_types_fill_1d(psetcols, value, Ptr[float](flx_net_p))
    _physics_types_fill_1d(psetcols, value, Ptr[float](te_tnd_p))
    _physics_types_fill_1d(psetcols, value, Ptr[float](tw_tnd_p))


@export
def physics_state_copy_codon(
    ncol: int,
    psetcols: int,
    pver: int,
    pcnst: int,
    src_lat_p: cobj,
    dst_lat_p: cobj,
    src_lon_p: cobj,
    dst_lon_p: cobj,
    src_ps_p: cobj,
    dst_ps_p: cobj,
    src_phis_p: cobj,
    dst_phis_p: cobj,
    src_te_ini_p: cobj,
    dst_te_ini_p: cobj,
    src_te_cur_p: cobj,
    dst_te_cur_p: cobj,
    src_tw_ini_p: cobj,
    dst_tw_ini_p: cobj,
    src_tw_cur_p: cobj,
    dst_tw_cur_p: cobj,
    src_psdry_p: cobj,
    dst_psdry_p: cobj,
    src_t_p: cobj,
    dst_t_p: cobj,
    src_u_p: cobj,
    dst_u_p: cobj,
    src_v_p: cobj,
    dst_v_p: cobj,
    src_s_p: cobj,
    dst_s_p: cobj,
    src_omega_p: cobj,
    dst_omega_p: cobj,
    src_pmid_p: cobj,
    dst_pmid_p: cobj,
    src_pdel_p: cobj,
    dst_pdel_p: cobj,
    src_rpdel_p: cobj,
    dst_rpdel_p: cobj,
    src_lnpmid_p: cobj,
    dst_lnpmid_p: cobj,
    src_exner_p: cobj,
    dst_exner_p: cobj,
    src_zm_p: cobj,
    dst_zm_p: cobj,
    src_lnpmiddry_p: cobj,
    dst_lnpmiddry_p: cobj,
    src_pmiddry_p: cobj,
    dst_pmiddry_p: cobj,
    src_pdeldry_p: cobj,
    dst_pdeldry_p: cobj,
    src_rpdeldry_p: cobj,
    dst_rpdeldry_p: cobj,
    src_pint_p: cobj,
    dst_pint_p: cobj,
    src_lnpint_p: cobj,
    dst_lnpint_p: cobj,
    src_zi_p: cobj,
    dst_zi_p: cobj,
    src_pintdry_p: cobj,
    dst_pintdry_p: cobj,
    src_lnpintdry_p: cobj,
    dst_lnpintdry_p: cobj,
    src_q_p: cobj,
    dst_q_p: cobj,
):
    _physics_types_copy_1d(ncol, Ptr[float](src_lat_p), Ptr[float](dst_lat_p))
    _physics_types_copy_1d(ncol, Ptr[float](src_lon_p), Ptr[float](dst_lon_p))
    _physics_types_copy_1d(ncol, Ptr[float](src_ps_p), Ptr[float](dst_ps_p))
    _physics_types_copy_1d(ncol, Ptr[float](src_phis_p), Ptr[float](dst_phis_p))
    _physics_types_copy_1d(ncol, Ptr[float](src_te_ini_p), Ptr[float](dst_te_ini_p))
    _physics_types_copy_1d(ncol, Ptr[float](src_te_cur_p), Ptr[float](dst_te_cur_p))
    _physics_types_copy_1d(ncol, Ptr[float](src_tw_ini_p), Ptr[float](dst_tw_ini_p))
    _physics_types_copy_1d(ncol, Ptr[float](src_tw_cur_p), Ptr[float](dst_tw_cur_p))
    _physics_types_copy_1d(ncol, Ptr[float](src_psdry_p), Ptr[float](dst_psdry_p))

    _physics_types_copy_2d_ncol(ncol, psetcols, pver, Ptr[float](src_t_p), Ptr[float](dst_t_p))
    _physics_types_copy_2d_ncol(ncol, psetcols, pver, Ptr[float](src_u_p), Ptr[float](dst_u_p))
    _physics_types_copy_2d_ncol(ncol, psetcols, pver, Ptr[float](src_v_p), Ptr[float](dst_v_p))
    _physics_types_copy_2d_ncol(ncol, psetcols, pver, Ptr[float](src_s_p), Ptr[float](dst_s_p))
    _physics_types_copy_2d_ncol(ncol, psetcols, pver, Ptr[float](src_omega_p), Ptr[float](dst_omega_p))
    _physics_types_copy_2d_ncol(ncol, psetcols, pver, Ptr[float](src_pmid_p), Ptr[float](dst_pmid_p))
    _physics_types_copy_2d_ncol(ncol, psetcols, pver, Ptr[float](src_pdel_p), Ptr[float](dst_pdel_p))
    _physics_types_copy_2d_ncol(ncol, psetcols, pver, Ptr[float](src_rpdel_p), Ptr[float](dst_rpdel_p))
    _physics_types_copy_2d_ncol(ncol, psetcols, pver, Ptr[float](src_lnpmid_p), Ptr[float](dst_lnpmid_p))
    _physics_types_copy_2d_ncol(ncol, psetcols, pver, Ptr[float](src_exner_p), Ptr[float](dst_exner_p))
    _physics_types_copy_2d_ncol(ncol, psetcols, pver, Ptr[float](src_zm_p), Ptr[float](dst_zm_p))
    _physics_types_copy_2d_ncol(ncol, psetcols, pver, Ptr[float](src_lnpmiddry_p), Ptr[float](dst_lnpmiddry_p))
    _physics_types_copy_2d_ncol(ncol, psetcols, pver, Ptr[float](src_pmiddry_p), Ptr[float](dst_pmiddry_p))
    _physics_types_copy_2d_ncol(ncol, psetcols, pver, Ptr[float](src_pdeldry_p), Ptr[float](dst_pdeldry_p))
    _physics_types_copy_2d_ncol(ncol, psetcols, pver, Ptr[float](src_rpdeldry_p), Ptr[float](dst_rpdeldry_p))

    _physics_types_copy_2d_ncol(ncol, psetcols, pver + 1, Ptr[float](src_pint_p), Ptr[float](dst_pint_p))
    _physics_types_copy_2d_ncol(ncol, psetcols, pver + 1, Ptr[float](src_lnpint_p), Ptr[float](dst_lnpint_p))
    _physics_types_copy_2d_ncol(ncol, psetcols, pver + 1, Ptr[float](src_zi_p), Ptr[float](dst_zi_p))
    _physics_types_copy_2d_ncol(ncol, psetcols, pver + 1, Ptr[float](src_pintdry_p), Ptr[float](dst_pintdry_p))
    _physics_types_copy_2d_ncol(ncol, psetcols, pver + 1, Ptr[float](src_lnpintdry_p), Ptr[float](dst_lnpintdry_p))

    _physics_types_copy_3d_ncol(ncol, psetcols, pver, pcnst, Ptr[float](src_q_p), Ptr[float](dst_q_p))


@export
def physics_ptend_init_codon(
    ls_present: int,
    ls_value: int,
    lu_present: int,
    lu_value: int,
    lv_present: int,
    lv_value: int,
    lq_present: int,
) -> int:
    if ls_present == 0 and lu_present == 0 and lv_present == 0 and lq_present == 0:
        return 8

    policy = 0
    if ls_present != 0 and ls_value != 0:
        policy = policy | 1
    if lu_present != 0 and lu_value != 0:
        policy = policy | 2
    if lv_present != 0 and lv_value != 0:
        policy = policy | 4
    return policy


@export
def physics_fill_real_1d_codon(n: int, value: float, arr_p: cobj):
    _physics_types_fill_1d(n, value, Ptr[float](arr_p))


@export
def physics_fill_real_2d_codon(psetcols: int, pver: int, value: float, arr_p: cobj):
    _physics_types_fill_2d(psetcols, pver, value, Ptr[float](arr_p))


@export
def physics_fill_real_3d_codon(psetcols: int, pver: int, pcnst: int, value: float, arr_p: cobj):
    _physics_types_fill_3d(psetcols, pver, pcnst, value, Ptr[float](arr_p))


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
def state_cnst_min_nz_codon(
    ncol: int,
    psetcols: int,
    pver: int,
    pcnst: int,
    q_p: cobj,
    lim: float,
    qix: int,
    numix: int,
):
    q = Ptr[float](q_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx_q = _field3_idx(i, k, qix, psetcols, pver)
            if q[idx_q] < lim:
                q[idx_q] = 0.0
                if numix > 0:
                    q[_field3_idx(i, k, numix, psetcols, pver)] = 0.0


@export
def physics_update_field_codon(
    ncol: int,
    psetcols: int,
    pver: int,
    top_level: int,
    bot_level: int,
    dt: float,
    update_tend: int,
    state_field_p: cobj,
    ptend_field_p: cobj,
    tend_field_p: cobj,
):
    state_field = Ptr[float](state_field_p)
    ptend_field = Ptr[float](ptend_field_p)
    tend_field = Ptr[float](tend_field_p)

    for k in range(top_level, bot_level + 1):
        for i in range(1, ncol + 1):
            idx = _field2_idx(i, k, psetcols)
            state_field[idx] = state_field[idx] + ptend_field[idx] * dt
            if update_tend != 0:
                tend_field[idx] = tend_field[idx] + ptend_field[idx]


@export
def physics_update_q_codon(
    ncol: int,
    psetcols: int,
    pver: int,
    pcnst: int,
    top_level: int,
    bot_level: int,
    dt: float,
    m: int,
    is_number: int,
    state_q_p: cobj,
    ptend_q_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    ptend_q = Ptr[float](ptend_q_p)

    for k in range(top_level, bot_level + 1):
        for i in range(1, ncol + 1):
            idx = _field3_idx(i, k, m, psetcols, pver)
            state_q[idx] = state_q[idx] + ptend_q[idx] * dt
            if is_number != 0:
                state_q[idx] = max(1.0e-12, state_q[idx])
                state_q[idx] = min(1.0e10, state_q[idx])


@export
def physics_update_s_codon(
    ncol: int,
    psetcols: int,
    pver: int,
    top_level: int,
    bot_level: int,
    dt: float,
    update_tend: int,
    state_s_p: cobj,
    ptend_s_p: cobj,
    tend_dtdt_p: cobj,
    cpairv_loc_p: cobj,
):
    state_s = Ptr[float](state_s_p)
    ptend_s = Ptr[float](ptend_s_p)
    tend_dtdt = Ptr[float](tend_dtdt_p)
    cpairv_loc = Ptr[float](cpairv_loc_p)

    for k in range(top_level, bot_level + 1):
        for i in range(1, ncol + 1):
            idx = _field2_idx(i, k, psetcols)
            state_s[idx] = state_s[idx] + ptend_s[idx] * dt
            if update_tend != 0:
                tend_dtdt[idx] = tend_dtdt[idx] + ptend_s[idx] / cpairv_loc[idx]


@export
def physics_ptend_sum_field_codon(
    ncol: int,
    psetcols: int,
    top_level: int,
    bot_level: int,
    src_field_p: cobj,
    dst_field_p: cobj,
    src_flx_srf_p: cobj,
    dst_flx_srf_p: cobj,
    src_flx_top_p: cobj,
    dst_flx_top_p: cobj,
):
    src_field = Ptr[float](src_field_p)
    dst_field = Ptr[float](dst_field_p)
    src_flx_srf = Ptr[float](src_flx_srf_p)
    dst_flx_srf = Ptr[float](dst_flx_srf_p)
    src_flx_top = Ptr[float](src_flx_top_p)
    dst_flx_top = Ptr[float](dst_flx_top_p)

    for k in range(top_level, bot_level + 1):
        for i in range(1, ncol + 1):
            idx = _field2_idx(i, k, psetcols)
            dst_field[idx] = dst_field[idx] + src_field[idx]

    for i in range(1, ncol + 1):
        idx = _idx(i)
        dst_flx_srf[idx] = dst_flx_srf[idx] + src_flx_srf[idx]
        dst_flx_top[idx] = dst_flx_top[idx] + src_flx_top[idx]


@inline
def _physics_ptend_zero_field(
    psetcols: int,
    pver: int,
    field_p: cobj,
    flx_srf_p: cobj,
    flx_top_p: cobj,
):
    _physics_types_zero_2d(psetcols, pver, Ptr[float](field_p))
    _physics_types_zero_1d(psetcols, Ptr[float](flx_srf_p))
    _physics_types_zero_1d(psetcols, Ptr[float](flx_top_p))


@inline
def _physics_ptend_scale_field(
    ncol: int,
    psetcols: int,
    pver: int,
    top_level: int,
    bot_level: int,
    fac: float,
    field_kind: int,
    qidx: int,
    field_p: cobj,
    flx_srf_p: cobj,
    flx_top_p: cobj,
):
    field = Ptr[float](field_p)
    flx_srf = Ptr[float](flx_srf_p)
    flx_top = Ptr[float](flx_top_p)

    if field_kind == 3:
        for k in range(top_level, bot_level + 1):
            for i in range(1, ncol + 1):
                idx = _field3_idx(i, k, qidx, psetcols, pver)
                field[idx] = field[idx] * fac
        for i in range(1, ncol + 1):
            flx_srf[_field2_idx(i, qidx, psetcols)] = flx_srf[_field2_idx(i, qidx, psetcols)] * fac
            flx_top[_field2_idx(i, qidx, psetcols)] = flx_top[_field2_idx(i, qidx, psetcols)] * fac
    else:
        for k in range(top_level, bot_level + 1):
            for i in range(1, ncol + 1):
                idx = _field2_idx(i, k, psetcols)
                field[idx] = field[idx] * fac
        for i in range(1, ncol + 1):
            flx_srf[_idx(i)] = flx_srf[_idx(i)] * fac
            flx_top[_idx(i)] = flx_top[_idx(i)] * fac


@inline
def _physics_ptend_sum_field(
    ncol: int,
    psetcols: int,
    pver: int,
    top_level: int,
    bot_level: int,
    field_kind: int,
    qidx: int,
    src_field_p: cobj,
    dst_field_p: cobj,
    src_flx_srf_p: cobj,
    dst_flx_srf_p: cobj,
    src_flx_top_p: cobj,
    dst_flx_top_p: cobj,
):
    src_field = Ptr[float](src_field_p)
    dst_field = Ptr[float](dst_field_p)
    src_flx_srf = Ptr[float](src_flx_srf_p)
    dst_flx_srf = Ptr[float](dst_flx_srf_p)
    src_flx_top = Ptr[float](src_flx_top_p)
    dst_flx_top = Ptr[float](dst_flx_top_p)

    if field_kind == 3:
        for k in range(top_level, bot_level + 1):
            for i in range(1, ncol + 1):
                idx = _field3_idx(i, k, qidx, psetcols, pver)
                dst_field[idx] = dst_field[idx] + src_field[idx]
        for i in range(1, ncol + 1):
            idx = _field2_idx(i, qidx, psetcols)
            dst_flx_srf[idx] = dst_flx_srf[idx] + src_flx_srf[idx]
            dst_flx_top[idx] = dst_flx_top[idx] + src_flx_top[idx]
    else:
        for k in range(top_level, bot_level + 1):
            for i in range(1, ncol + 1):
                idx = _field2_idx(i, k, psetcols)
                dst_field[idx] = dst_field[idx] + src_field[idx]
        for i in range(1, ncol + 1):
            idx = _idx(i)
            dst_flx_srf[idx] = dst_flx_srf[idx] + src_flx_srf[idx]
            dst_flx_top[idx] = dst_flx_top[idx] + src_flx_top[idx]


@export
def physics_ptend_reset_shell_codon(
    psetcols: int,
    pver: int,
    pcnst: int,
    ls: int,
    lu: int,
    lv: int,
    lq_p: cobj,
    s_p: cobj,
    u_p: cobj,
    v_p: cobj,
    q_p: cobj,
    hflux_srf_p: cobj,
    hflux_top_p: cobj,
    taux_srf_p: cobj,
    taux_top_p: cobj,
    tauy_srf_p: cobj,
    tauy_top_p: cobj,
    cflx_srf_p: cobj,
    cflx_top_p: cobj,
):
    lq = Ptr[int](lq_p)

    if ls != 0:
        _physics_ptend_zero_field(psetcols, pver, s_p, hflux_srf_p, hflux_top_p)
    if lu != 0:
        _physics_ptend_zero_field(psetcols, pver, u_p, taux_srf_p, taux_top_p)
    if lv != 0:
        _physics_ptend_zero_field(psetcols, pver, v_p, tauy_srf_p, tauy_top_p)

    any_lq = False
    for m in range(1, pcnst + 1):
        if lq[_idx(m)] != 0:
            any_lq = True
    if any_lq:
        _physics_types_zero_3d(psetcols, pver, pcnst, Ptr[float](q_p))
        _physics_types_zero_2d(psetcols, pcnst, Ptr[float](cflx_srf_p))
        _physics_types_zero_2d(psetcols, pcnst, Ptr[float](cflx_top_p))


@export
def physics_ptend_scale_shell_codon(
    ncol: int,
    psetcols: int,
    pver: int,
    pcnst: int,
    top_level: int,
    bot_level: int,
    fac: float,
    ls: int,
    lu: int,
    lv: int,
    lq_p: cobj,
    s_p: cobj,
    u_p: cobj,
    v_p: cobj,
    q_p: cobj,
    hflux_srf_p: cobj,
    hflux_top_p: cobj,
    taux_srf_p: cobj,
    taux_top_p: cobj,
    tauy_srf_p: cobj,
    tauy_top_p: cobj,
    cflx_srf_p: cobj,
    cflx_top_p: cobj,
):
    lq = Ptr[int](lq_p)

    if lu != 0:
        _physics_ptend_scale_field(ncol, psetcols, pver, top_level, bot_level, fac, 2, 0, u_p, taux_srf_p, taux_top_p)
    if lv != 0:
        _physics_ptend_scale_field(ncol, psetcols, pver, top_level, bot_level, fac, 2, 0, v_p, tauy_srf_p, tauy_top_p)
    if ls != 0:
        _physics_ptend_scale_field(ncol, psetcols, pver, top_level, bot_level, fac, 2, 0, s_p, hflux_srf_p, hflux_top_p)

    for m in range(1, pcnst + 1):
        if lq[_idx(m)] != 0:
            _physics_ptend_scale_field(ncol, psetcols, pver, top_level, bot_level, fac, 3, m, q_p, cflx_srf_p, cflx_top_p)


@export
def physics_ptend_sum_shell_codon(
    ncol: int,
    psetcols: int,
    pver: int,
    pcnst: int,
    top_level: int,
    bot_level: int,
    ls: int,
    lu: int,
    lv: int,
    lq_p: cobj,
    src_s_p: cobj,
    dst_s_p: cobj,
    src_u_p: cobj,
    dst_u_p: cobj,
    src_v_p: cobj,
    dst_v_p: cobj,
    src_q_p: cobj,
    dst_q_p: cobj,
    src_hflux_srf_p: cobj,
    dst_hflux_srf_p: cobj,
    src_hflux_top_p: cobj,
    dst_hflux_top_p: cobj,
    src_taux_srf_p: cobj,
    dst_taux_srf_p: cobj,
    src_taux_top_p: cobj,
    dst_taux_top_p: cobj,
    src_tauy_srf_p: cobj,
    dst_tauy_srf_p: cobj,
    src_tauy_top_p: cobj,
    dst_tauy_top_p: cobj,
    src_cflx_srf_p: cobj,
    dst_cflx_srf_p: cobj,
    src_cflx_top_p: cobj,
    dst_cflx_top_p: cobj,
):
    lq = Ptr[int](lq_p)

    if lu != 0:
        _physics_ptend_sum_field(ncol, psetcols, pver, top_level, bot_level, 2, 0,
                                 src_u_p, dst_u_p, src_taux_srf_p, dst_taux_srf_p,
                                 src_taux_top_p, dst_taux_top_p)
    if lv != 0:
        _physics_ptend_sum_field(ncol, psetcols, pver, top_level, bot_level, 2, 0,
                                 src_v_p, dst_v_p, src_tauy_srf_p, dst_tauy_srf_p,
                                 src_tauy_top_p, dst_tauy_top_p)
    if ls != 0:
        _physics_ptend_sum_field(ncol, psetcols, pver, top_level, bot_level, 2, 0,
                                 src_s_p, dst_s_p, src_hflux_srf_p, dst_hflux_srf_p,
                                 src_hflux_top_p, dst_hflux_top_p)

    for m in range(1, pcnst + 1):
        if lq[_idx(m)] != 0:
            _physics_ptend_sum_field(ncol, psetcols, pver, top_level, bot_level, 3, m,
                                     src_q_p, dst_q_p, src_cflx_srf_p, dst_cflx_srf_p,
                                     src_cflx_top_p, dst_cflx_top_p)


@export
def physics_set_state_pdry_codon(
    ncol: int,
    psetcols: int,
    pver: int,
    do_pdeld_calc: int,
    psdry_p: cobj,
    pint_p: cobj,
    pdel_p: cobj,
    q_p: cobj,
    pdeldry_p: cobj,
    pintdry_p: cobj,
    pmiddry_p: cobj,
    rpdeldry_p: cobj,
    lnpmiddry_p: cobj,
    lnpintdry_p: cobj,
):
    psdry = Ptr[float](psdry_p)
    pint = Ptr[float](pint_p)
    pdel = Ptr[float](pdel_p)
    q = Ptr[float](q_p)
    pdeldry = Ptr[float](pdeldry_p)
    pintdry = Ptr[float](pintdry_p)
    pmiddry = Ptr[float](pmiddry_p)
    rpdeldry = Ptr[float](rpdeldry_p)
    lnpmiddry = Ptr[float](lnpmiddry_p)
    lnpintdry = Ptr[float](lnpintdry_p)

    for i in range(1, ncol + 1):
        psdry[_idx(i)] = pint[_field2_idx(i, 1, psetcols)]
        pintdry[_field2_idx(i, 1, psetcols)] = pint[_field2_idx(i, 1, psetcols)]

    if do_pdeld_calc != 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx2 = _field2_idx(i, k, psetcols)
                pdeldry[idx2] = pdel[idx2] * (1.0 - q[_field3_idx(i, k, 1, psetcols, pver)])

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _field2_idx(i, k, psetcols)
            idx2p1 = _field2_idx(i, k + 1, psetcols)
            pintdry[idx2p1] = pintdry[idx2] + pdeldry[idx2]
            pmiddry[idx2] = (pintdry[idx2p1] + pintdry[idx2]) / 2.0
            psdry[_idx(i)] = psdry[_idx(i)] + pdeldry[idx2]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _field2_idx(i, k, psetcols)
            rpdeldry[idx2] = 1.0 / pdeldry[idx2]
            lnpmiddry[idx2] = log(pmiddry[idx2])

    for k in range(1, pver + 2):
        for i in range(1, ncol + 1):
            idx2 = _field2_idx(i, k, psetcols)
            lnpintdry[idx2] = log(pintdry[idx2])


@export
def physics_set_wet_to_dry_constituent_codon(
    ncol: int,
    psetcols: int,
    pver: int,
    q_p: cobj,
    pdel_p: cobj,
    pdeldry_p: cobj,
):
    q = Ptr[float](q_p)
    pdel = Ptr[float](pdel_p)
    pdeldry = Ptr[float](pdeldry_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _field2_idx(i, k, psetcols)
            q[idx2] = q[idx2] * pdel[idx2] / pdeldry[idx2]


@export
def physics_set_dry_to_wet_constituent_codon(
    ncol: int,
    psetcols: int,
    pver: int,
    q_p: cobj,
    pdeldry_p: cobj,
    pdel_p: cobj,
):
    q = Ptr[float](q_p)
    pdeldry = Ptr[float](pdeldry_p)
    pdel = Ptr[float](pdel_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _field2_idx(i, k, psetcols)
            q[idx2] = q[idx2] * pdeldry[idx2] / pdel[idx2]


@export
def physics_init_geo_unique_maps_codon(
    ncol: int,
    psetcols: int,
    lat_p: cobj,
    lon_p: cobj,
    ulat_p: cobj,
    ulon_p: cobj,
    latmapback_p: cobj,
    lonmapback_p: cobj,
    ulatcnt_p: cobj,
    uloncnt_p: cobj,
):
    lat = Ptr[float](lat_p)
    lon = Ptr[float](lon_p)
    ulat = Ptr[float](ulat_p)
    ulon = Ptr[float](ulon_p)
    latmapback = Ptr[i32](latmapback_p)
    lonmapback = Ptr[i32](lonmapback_p)
    ulatcnt_out = Ptr[int](ulatcnt_p)
    uloncnt_out = Ptr[int](uloncnt_p)

    for i in range(0, psetcols):
        ulat[i] = -999.0
        ulon[i] = -999.0
        latmapback[i] = i32(0)
        lonmapback[i] = i32(0)

    ulatcnt = 0
    uloncnt = 0
    match = False

    for i in range(1, ncol + 1):
        for j in range(1, ulatcnt + 1):
            if lat[_idx(i)] == ulat[_idx(j)]:
                match = True
                latmapback[_idx(i)] = i32(j)
        if not match:
            ulatcnt += 1
            ulat[_idx(ulatcnt)] = lat[_idx(i)]
            latmapback[_idx(i)] = i32(ulatcnt)

        match = False
        for j in range(1, uloncnt + 1):
            if lon[_idx(i)] == ulon[_idx(j)]:
                match = True
                lonmapback[_idx(i)] = i32(j)
        if not match:
            uloncnt += 1
            ulon[_idx(uloncnt)] = lon[_idx(i)]
            lonmapback[_idx(i)] = i32(uloncnt)
        match = False

    ulatcnt_out[0] = ulatcnt
    uloncnt_out[0] = uloncnt


@export
def physics_copy_real_1d_codon(n: int, src_p: cobj, dst_p: cobj):
    src = Ptr[float](src_p)
    dst = Ptr[float](dst_p)

    for i in range(0, n):
        dst[i] = src[i]


@export
def physics_copy_real_2d_codon(ncol: int, ld1: int, nlev: int, src_p: cobj, dst_p: cobj):
    src = Ptr[float](src_p)
    dst = Ptr[float](dst_p)

    for k in range(1, nlev + 1):
        for i in range(1, ncol + 1):
            idx = _field2_idx(i, k, ld1)
            dst[idx] = src[idx]


@export
def physics_copy_real_3d_codon(ncol: int, ld1: int, nlev: int, n3: int, src_p: cobj, dst_p: cobj):
    src = Ptr[float](src_p)
    dst = Ptr[float](dst_p)

    for m in range(1, n3 + 1):
        for k in range(1, nlev + 1):
            for i in range(1, ncol + 1):
                idx = _field3_idx(i, k, m, ld1, nlev)
                dst[idx] = src[idx]


@export
def phys_grid_get_gcol_all_codon(ncols: int, out_dim: int, src_p: cobj, dst_p: cobj):
    src = Ptr[i32](src_p)
    dst = Ptr[i32](dst_p)

    for i in range(out_dim):
        dst[i] = i32(-1)

    for i in range(ncols):
        dst[i] = src[i]


@export
def get_gcol_all_p_codon(ncols: int, out_dim: int, src_p: cobj, dst_p: cobj):
    phys_grid_get_gcol_all_codon(ncols, out_dim, src_p, dst_p)


@export
def phys_grid_int_scalar_codon(value: int) -> int:
    return value


@export
def get_nlcols_p_codon(value: int) -> int:
    return phys_grid_int_scalar_codon(value)


@export
def get_gcol_p_codon(value: int) -> int:
    return phys_grid_int_scalar_codon(value)


@export
def get_ncols_p_codon(value: int) -> int:
    return phys_grid_int_scalar_codon(value)


@export
def phys_grid_bool_scalar_codon(value: int) -> int:
    if value != 0:
        return 1
    return 0


@export
def phys_grid_initialized_codon(value: int) -> int:
    return phys_grid_bool_scalar_codon(value)


@export
def phys_grid_init_codon(stage: int) -> int:
    return stage


@export
def create_chunks_codon(stage: int) -> int:
    return stage


@export
def assign_chunks_codon(stage: int) -> int:
    return stage


@export
def phys_grid_defaultopts_codon(
    has_lbal: int,
    has_twin: int,
    has_alltoall: int,
    has_chunks: int,
    is_unstructured: int,
    def_lbal: int,
    def_twin_unstructured: int,
    def_twin_lonlat: int,
    def_alltoall: int,
    def_chunks: int,
    out_p: cobj,
):
    out = Ptr[int](out_p)

    if has_lbal != 0:
        out[0] = def_lbal
    if has_twin != 0:
        if is_unstructured != 0:
            out[1] = def_twin_unstructured
        else:
            out[1] = def_twin_lonlat
    if has_alltoall != 0:
        out[2] = def_alltoall
    if has_chunks != 0:
        out[3] = def_chunks


@export
def phys_grid_setopts_codon(
    has_lbal: int,
    lbal_in: int,
    has_twin: int,
    twin_in: int,
    has_alltoall: int,
    alltoall_in: int,
    has_chunks: int,
    chunks_in: int,
    current_lbal: int,
    current_twin: int,
    current_alltoall: int,
    current_chunks: int,
    min_lbal: int,
    max_lbal: int,
    min_twin: int,
    max_twin: int,
    min_alltoall: int,
    max_alltoall: int,
    allow_mod_alltoall: int,
    modmin_alltoall: int,
    modmax_alltoall: int,
    min_chunks: int,
    state_p: cobj,
    error_p: cobj,
):
    state = Ptr[int](state_p)
    error = Ptr[int](error_p)

    lbal = current_lbal
    twin = current_twin
    alltoall = current_alltoall
    chunks = current_chunks
    error_code = 0

    if has_lbal != 0:
        lbal = lbal_in
        if lbal < min_lbal or lbal > max_lbal:
            error_code = 1

    if error_code == 0 and has_twin != 0:
        twin = twin_in
        if twin < min_twin or twin > max_twin:
            error_code = 2

    if error_code == 0 and has_alltoall != 0:
        alltoall = alltoall_in
        standard_invalid = alltoall < min_alltoall or alltoall > max_alltoall
        mod_invalid = True
        if allow_mod_alltoall != 0:
            mod_invalid = alltoall < modmin_alltoall or alltoall > modmax_alltoall
        if standard_invalid and mod_invalid:
            error_code = 3

    if error_code == 0 and has_chunks != 0:
        chunks = chunks_in
        if chunks < min_chunks:
            error_code = 4

    state[0] = lbal
    state[1] = twin
    state[2] = alltoall
    state[3] = chunks
    if lbal == 3:
        state[4] = 1
    else:
        state[4] = 0
    error[0] = error_code


@export
def phys_grid_count_valid_cols_codon(ngcols: int, clon_d_p: cobj) -> int:
    return _phys_grid.phys_grid_count_valid_cols_codon(ngcols, clon_d_p)


@export
def phys_grid_count_unique_sorted_real_codon(ncols: int, cdex_p: cobj, coord_p: cobj) -> int:
    return _phys_grid.phys_grid_count_unique_sorted_real_codon(ncols, cdex_p, coord_p)


@export
def phys_grid_fill_unique_sorted_real_codon(
    ncols: int,
    cdex_p: cobj,
    coord_p: cobj,
    unique_p: cobj,
    counts_p: cobj,
):
    _phys_grid.phys_grid_fill_unique_sorted_real_codon(ncols, cdex_p, coord_p, unique_p, counts_p)


@export
def phys_grid_prefix_counts_codon(n: int, counts_p: cobj, idx_p: cobj):
    _phys_grid.phys_grid_prefix_counts_codon(n, counts_p, idx_p)


@export
def phys_grid_fill_real_pair_codon(
    n: int,
    first_value: float,
    second_value: float,
    first_p: cobj,
    second_p: cobj,
):
    _phys_grid.phys_grid_fill_real_pair_codon(n, first_value, second_value, first_p, second_p)


@export
def phys_grid_init_lat_map_codon(
    ngcols: int,
    ncols_p: int,
    clat_tot: int,
    has_latlon_map: int,
    cdex_p: cobj,
    clat_d_p: cobj,
    clat_p_p: cobj,
    lat_p_p: cobj,
    dyn_to_latlon_gcol_map_p: cobj,
    latlon_to_dyn_gcol_map_p: cobj,
):
    _phys_grid.phys_grid_init_lat_map_codon(
        ngcols,
        ncols_p,
        clat_tot,
        has_latlon_map,
        cdex_p,
        clat_d_p,
        clat_p_p,
        lat_p_p,
        dyn_to_latlon_gcol_map_p,
        latlon_to_dyn_gcol_map_p,
    )


@export
def phys_grid_init_lon_map_codon(
    ngcols: int,
    ncols_p: int,
    clon_tot: int,
    has_lonlat_map: int,
    cdex_p: cobj,
    clon_d_p: cobj,
    clon_p_p: cobj,
    lon_p_p: cobj,
    lonlat_to_dyn_gcol_map_p: cobj,
):
    _phys_grid.phys_grid_init_lon_map_codon(
        ngcols,
        ncols_p,
        clon_tot,
        has_lonlat_map,
        cdex_p,
        clon_d_p,
        clon_p_p,
        lon_p_p,
        lonlat_to_dyn_gcol_map_p,
    )


@export
def phys_grid_zero_proc_counts_codon(npes: int, chunk_counts_p: cobj, col_counts_p: cobj):
    _phys_grid.phys_grid_zero_proc_counts_codon(npes, chunk_counts_p, col_counts_p)


@export
def phys_grid_proc_prefix_offsets_codon(
    npes: int,
    start_value: int,
    set_final: int,
    chunk_counts_p: cobj,
    col_counts_p: cobj,
    pchunkid_p: cobj,
    gs_col_offset_p: cobj,
):
    _phys_grid.phys_grid_proc_prefix_offsets_codon(
        npes,
        start_value,
        set_final,
        chunk_counts_p,
        col_counts_p,
        pchunkid_p,
        gs_col_offset_p,
    )


@export
def phys_grid_process_bin_sort_codon(
    nchunks: int,
    lastblock: int,
    chunk_owner_p: cobj,
    chunk_ncols_p: cobj,
    pchunkid_p: cobj,
    gs_col_offset_p: cobj,
    chunk_lcid_p: cobj,
    pgcol_chunk_p: cobj,
    pgcol_ccol_p: cobj,
):
    _phys_grid.phys_grid_process_bin_sort_codon(
        nchunks,
        lastblock,
        chunk_owner_p,
        chunk_ncols_p,
        pchunkid_p,
        gs_col_offset_p,
        chunk_lcid_p,
        pgcol_chunk_p,
        pgcol_ccol_p,
    )


@export
def phys_grid_lchunk_gcol_copy_codon(ncols: int, src_gcol_p: cobj, dst_gcol_p: cobj):
    _phys_grid.phys_grid_lchunk_gcol_copy_codon(ncols, src_gcol_p, dst_gcol_p)


@export
def phys_grid_lchunk_area_wght_codon(
    ncols: int,
    gcol_p: cobj,
    area_d_p: cobj,
    wght_d_p: cobj,
    area_p: cobj,
    wght_p: cobj,
):
    _phys_grid.phys_grid_lchunk_area_wght_codon(ncols, gcol_p, area_d_p, wght_d_p, area_p, wght_p)


@export
def phys_grid_count_smp_procs_codon(
    npes: int,
    nsmpx: int,
    proc_smp_mapx_p: cobj,
    nsmpprocs_p: cobj,
) -> int:
    return _phys_grid.phys_grid_count_smp_procs_codon(npes, nsmpx, proc_smp_mapx_p, nsmpprocs_p)


@export
def phys_grid_create_chunks_thread_counts_codon(
    npes: int,
    nsmpx: int,
    proc_smp_mapx_p: cobj,
    npthreads_p: cobj,
    nsmpthreads_p: cobj,
):
    _phys_grid.phys_grid_create_chunks_thread_counts_codon(
        npes, nsmpx, proc_smp_mapx_p, npthreads_p, nsmpthreads_p
    )


@export
def phys_grid_create_chunks_shape_codon(
    nsmpx: int,
    pcols: int,
    chunks_per_thread: int,
    nsmpcolumns_p: cobj,
    nsmpthreads_p: cobj,
    nsmpchunks_p: cobj,
    maxcol_chk_p: cobj,
    maxcol_chks_p: cobj,
) -> int:
    return _phys_grid.phys_grid_create_chunks_shape_codon(
        nsmpx,
        pcols,
        chunks_per_thread,
        nsmpcolumns_p,
        nsmpthreads_p,
        nsmpchunks_p,
        maxcol_chk_p,
        maxcol_chks_p,
    )


@export
def phys_grid_create_chunks_prefix_codon(
    nsmpx: int,
    nsmpchunks_p: cobj,
    cid_offset_p: cobj,
    local_cid_p: cobj,
):
    _phys_grid.phys_grid_create_chunks_prefix_codon(nsmpx, nsmpchunks_p, cid_offset_p, local_cid_p)


@export
def phys_grid_count_smp_columns_codon(
    nsmpx: int,
    ngcols_p: int,
    latlon_to_dyn_gcol_map_p: cobj,
    col_smp_mapx_p: cobj,
    nsmpcolumns_p: cobj,
):
    _phys_grid.phys_grid_count_smp_columns_codon(
        nsmpx, ngcols_p, latlon_to_dyn_gcol_map_p, col_smp_mapx_p, nsmpcolumns_p
    )


@export
def phys_grid_zero_int_array_codon(n: int, values_p: cobj):
    _phys_grid.phys_grid_zero_int_array_codon(n, values_p)


@export
def phys_grid_assign_chunks_zero_column_count_codon(
    smp: int,
    nsmpx: int,
    max_nproc_smpx: int,
    ntsks_smpx_p: cobj,
    smp_proc_mapx_p: cobj,
    column_count_p: cobj,
):
    _phys_grid.phys_grid_assign_chunks_zero_column_count_codon(
        smp, nsmpx, max_nproc_smpx, ntsks_smpx_p, smp_proc_mapx_p, column_count_p
    )


@export
def phys_grid_assign_chunks_select_owner_codon(
    smp: int,
    nsmpx: int,
    max_nproc_smpx: int,
    ntsks_smpx_p: cobj,
    smp_proc_mapx_p: cobj,
    cur_npchunks_p: cobj,
    npchunks_p: cobj,
    column_count_p: cobj,
) -> int:
    return _phys_grid.phys_grid_assign_chunks_select_owner_codon(
        smp,
        nsmpx,
        max_nproc_smpx,
        ntsks_smpx_p,
        smp_proc_mapx_p,
        cur_npchunks_p,
        npchunks_p,
        column_count_p,
    )


@export
def phys_grid_assign_chunks_commit_owner_codon(
    owner: int,
    ncols: int,
    cur_npchunks_p: cobj,
    gs_col_num_p: cobj,
):
    _phys_grid.phys_grid_assign_chunks_commit_owner_codon(owner, ncols, cur_npchunks_p, gs_col_num_p)


@export
def phys_grid_assign_chunks_smp_setup_codon(
    npes: int,
    nsmpx: int,
    max_nproc_smpx: int,
    proc_smp_mapx_p: cobj,
    npthreads_p: cobj,
    nsmpthreads_p: cobj,
    nsmpchunks_p: cobj,
    ntsks_smpx_p: cobj,
    smp_proc_mapx_p: cobj,
    cid_offset_p: cobj,
    ntmp1_smp_p: cobj,
    ntmp2_smp_p: cobj,
    ntmp3_smp_p: cobj,
    ntmp4_smp_p: cobj,
    npchunks_p: cobj,
):
    _phys_grid.phys_grid_assign_chunks_smp_setup_codon(
        npes,
        nsmpx,
        max_nproc_smpx,
        proc_smp_mapx_p,
        npthreads_p,
        nsmpthreads_p,
        nsmpchunks_p,
        ntsks_smpx_p,
        smp_proc_mapx_p,
        cid_offset_p,
        ntmp1_smp_p,
        ntmp2_smp_p,
        ntmp3_smp_p,
        ntmp4_smp_p,
        npchunks_p,
    )


@export
def phys_grid_assign_block_no_twin_codon(
    blksiz: int,
    pcols: int,
    smp: int,
    cols_p: cobj,
    cid_offset_p: cobj,
    local_cid_p: cobj,
    nsmpchunks_p: cobj,
    maxcol_chk_p: cobj,
    maxcol_chks_p: cobj,
    dyn_to_latlon_gcol_map_p: cobj,
    lon_p_p: cobj,
    lat_p_p: cobj,
    chunk_ncols_p: cobj,
    chunk_gcol_p: cobj,
    chunk_lon_p: cobj,
    chunk_lat_p: cobj,
    knuhcs_chunkid_p: cobj,
    knuhcs_col_p: cobj,
):
    _phys_grid.phys_grid_assign_block_no_twin_codon(
        blksiz,
        pcols,
        smp,
        cols_p,
        cid_offset_p,
        local_cid_p,
        nsmpchunks_p,
        maxcol_chk_p,
        maxcol_chks_p,
        dyn_to_latlon_gcol_map_p,
        lon_p_p,
        lat_p_p,
        chunk_ncols_p,
        chunk_gcol_p,
        chunk_lon_p,
        chunk_lat_p,
        knuhcs_chunkid_p,
        knuhcs_col_p,
    )


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
def get_lat_all_p_codon(ncols: int, src_p: cobj, dst_p: cobj):
    phys_grid_get_int_all_codon(ncols, src_p, dst_p)


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
def get_lon_all_p_codon(
    ncols: int,
    lat_p: cobj,
    gcol_p: cobj,
    dyn_to_latlon_gcol_map_p: cobj,
    clat_p_idx_p: cobj,
    dst_p: cobj,
):
    phys_grid_get_lon_all_codon(ncols, lat_p, gcol_p, dyn_to_latlon_gcol_map_p, clat_p_idx_p, dst_p)


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
def get_area_all_p_codon(ncols: int, src_p: cobj, dst_p: cobj):
    phys_grid_get_real_all_codon(ncols, src_p, dst_p)


@export
def get_wght_all_p_codon(ncols: int, src_p: cobj, dst_p: cobj):
    phys_grid_get_real_all_codon(ncols, src_p, dst_p)


@export
def phys_grid_get_lookup_real_all_codon(ncols: int, idx_p: cobj, lookup_p: cobj, dst_p: cobj):
    idx = Ptr[i32](idx_p)
    lookup = Ptr[float](lookup_p)
    dst = Ptr[float](dst_p)

    for i in range(ncols):
        dst[i] = lookup[int(idx[i]) - 1]


@export
def get_rlat_all_p_codon(ncols: int, idx_p: cobj, lookup_p: cobj, dst_p: cobj):
    phys_grid_get_lookup_real_all_codon(ncols, idx_p, lookup_p, dst_p)


@export
def get_rlon_all_p_codon(ncols: int, idx_p: cobj, lookup_p: cobj, dst_p: cobj):
    phys_grid_get_lookup_real_all_codon(ncols, idx_p, lookup_p, dst_p)


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
def phys_grid_pter_offsets_codon(
    ncols: int,
    nlvls: int,
    fdim: int,
    ldim: int,
    record_size: int,
    src_p: cobj,
    dst_p: cobj,
):
    src = Ptr[i32](src_p)
    dst = Ptr[i32](dst_p)

    for k in range(1, nlvls + 1):
        for i in range(1, ncols + 1):
            dst[(i - 1) + (k - 1) * fdim] = i32(1 + record_size * int(src[(i - 1) + (k - 1) * ncols]))
        for i in range(ncols + 1, fdim + 1):
            dst[(i - 1) + (k - 1) * fdim] = i32(-1)

    for k in range(nlvls + 1, ldim + 1):
        for i in range(1, fdim + 1):
            dst[(i - 1) + (k - 1) * fdim] = i32(-1)


@export
def block_to_chunk_send_pters_codon(
    ncols: int,
    nlvls: int,
    fdim: int,
    ldim: int,
    record_size: int,
    src_p: cobj,
    dst_p: cobj,
):
    phys_grid_pter_offsets_codon(ncols, nlvls, fdim, ldim, record_size, src_p, dst_p)


@export
def block_to_chunk_recv_pters_codon(
    ncols: int,
    nlvls: int,
    fdim: int,
    ldim: int,
    record_size: int,
    src_p: cobj,
    dst_p: cobj,
):
    phys_grid_pter_offsets_codon(ncols, nlvls, fdim, ldim, record_size, src_p, dst_p)


@export
def chunk_to_block_send_pters_codon(
    ncols: int,
    nlvls: int,
    fdim: int,
    ldim: int,
    record_size: int,
    src_p: cobj,
    dst_p: cobj,
):
    phys_grid_pter_offsets_codon(ncols, nlvls, fdim, ldim, record_size, src_p, dst_p)


@export
def chunk_to_block_recv_pters_codon(
    ncols: int,
    nlvls: int,
    fdim: int,
    ldim: int,
    record_size: int,
    src_p: cobj,
    dst_p: cobj,
):
    phys_grid_pter_offsets_codon(ncols, nlvls, fdim, ldim, record_size, src_p, dst_p)


@export
def phys_grid_transpose_counts_codon(
    npes: int,
    record_size: int,
    direction: int,
    block_num_p: cobj,
    chunk_num_p: cobj,
    sndcnts_p: cobj,
    sdispls_p: cobj,
    rcvcnts_p: cobj,
    rdispls_p: cobj,
):
    _phys_grid.phys_grid_transpose_counts_codon(
        npes,
        record_size,
        direction,
        block_num_p,
        chunk_num_p,
        sndcnts_p,
        sdispls_p,
        rcvcnts_p,
        rdispls_p,
    )


@export
def transpose_block_to_chunk_codon(stage: int) -> int:
    return stage


@export
def transpose_chunk_to_block_codon(stage: int) -> int:
    return stage


@export
def phys_grid_transpose_lopt_codon(
    phys_alltoall: int,
    max_nproc_smpx: int,
    nproc_busy_d: int,
    npes: int,
    has_window: int,
) -> int:
    return _phys_grid.phys_grid_transpose_lopt_codon(
        phys_alltoall,
        max_nproc_smpx,
        nproc_busy_d,
        npes,
        has_window,
    )


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
def trb_mtn_stress_init_codon(
    oro_in: float,
    z0fac_in: float,
    karman_in: float,
    gravit_in: float,
    rair_in: float,
    orocnst_p: cobj,
    z0fac_p: cobj,
    karman_p: cobj,
    gravit_p: cobj,
    rair_p: cobj,
):
    orocnst = Ptr[float](orocnst_p)
    z0fac = Ptr[float](z0fac_p)
    karman = Ptr[float](karman_p)
    gravit = Ptr[float](gravit_p)
    rair = Ptr[float](rair_p)

    orocnst[0] = oro_in
    z0fac[0] = z0fac_in
    karman[0] = karman_in
    gravit[0] = gravit_in
    rair[0] = rair_in


@export
def init_tms_codon(
    oro_in: float,
    z0fac_in: float,
    karman_in: float,
    gravit_in: float,
    rair_in: float,
    orocnst_p: cobj,
    z0fac_p: cobj,
    karman_p: cobj,
    gravit_p: cobj,
    rair_p: cobj,
):
    trb_mtn_stress_init_codon(
        oro_in,
        z0fac_in,
        karman_in,
        gravit_in,
        rair_in,
        orocnst_p,
        z0fac_p,
        karman_p,
        gravit_p,
        rair_p,
    )


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


@inline
def _ptr_ascii_to_str(n: int, ptr_p: cobj) -> str:
    ptr = Ptr[int](ptr_p)
    out = ""
    for i in range(n):
        out += chr(ptr[i])
    return out.strip()


@inline
def _strip_fortran_comment(line: str) -> str:
    pos = line.find("!")
    if pos >= 0:
        return line[:pos]
    return line


@inline
def _parse_fortran_real(value: str) -> float:
    return float(value.strip().replace("D", "E").replace("d", "e"))


@export
def cldfrc2m_readnl_codon(path_len: int, path_ascii_p: cobj, reals_p: cobj) -> int:
    path = _ptr_ascii_to_str(path_len, path_ascii_p)
    reals = Ptr[float](reals_p)

    f = open(path, "r")
    text = f.read()
    f.close()

    in_group = False
    found_group = False
    assignments = ""

    for raw_line in text.split("\n"):
        line = _strip_fortran_comment(raw_line).strip()
        lowered = line.lower()
        if not in_group:
            if lowered.startswith("&cldfrc2m_nl"):
                in_group = True
                found_group = True
                rest = line[len("&cldfrc2m_nl") :]
                if rest:
                    assignments += rest + ","
            continue

        slash = line.find("/")
        if slash >= 0:
            assignments += line[:slash] + ","
            break
        assignments += line + ","

    if not found_group:
        return 0

    for item in assignments.split(","):
        if "=" not in item:
            continue
        parts = item.split("=", 1)
        key = parts[0].strip().lower()
        value = parts[1].strip()
        if key == "cldfrc2m_rhmini":
            reals[0] = _parse_fortran_real(value)
        elif key == "cldfrc2m_rhmaxi":
            reals[1] = _parse_fortran_real(value)

    return 0


@export
def cldfrc2m_init_codon(
    rhmini_in: float,
    rhmaxi_in: float,
    rhminl_in: float,
    rhminl_adj_land_in: float,
    rhminh_in: float,
    premit_in: float,
    premib_in: float,
    icecrit_in: float,
    iceopt_in: int,
    rhmini_p: cobj,
    rhmaxi_p: cobj,
    rhminl_p: cobj,
    rhminl_adj_land_p: cobj,
    rhminh_p: cobj,
    premit_p: cobj,
    premib_p: cobj,
    icecrit_p: cobj,
    iceopt_p: cobj,
):
    Ptr[float](rhmini_p)[0] = rhmini_in
    Ptr[float](rhmaxi_p)[0] = rhmaxi_in
    Ptr[float](rhminl_p)[0] = rhminl_in
    Ptr[float](rhminl_adj_land_p)[0] = rhminl_adj_land_in
    Ptr[float](rhminh_p)[0] = rhminh_in
    Ptr[float](premit_p)[0] = premit_in
    Ptr[float](premib_p)[0] = premib_in
    Ptr[float](icecrit_p)[0] = icecrit_in
    Ptr[i32](iceopt_p)[0] = i32(iceopt_in)


@export
def cldfrc2m_pressure_regime_codon(p: float, premib: float, premit: float) -> int:
    if p >= premib:
        return 0
    if p < premit:
        return 1
    return 2


@export
def cldfrc2m_astg_pdf_zero_codon(pcols: int, a_p: cobj, ga_p: cobj):
    a = Ptr[float](a_p)
    ga = Ptr[float](ga_p)
    for i in range(pcols):
        a[i] = 0.0
        ga[i] = 0.0


@inline
def _cldfrc2m_astg_pdf_core(
    u: float,
    p: float,
    qv: float,
    landfrac: float,
    snowh: float,
    rhminl: float,
    rhminl_adj_land: float,
    rhminh: float,
    premib: float,
    premit: float,
):
    cldrh = 1.0
    pressure_regime = cldfrc2m_pressure_regime_codon(p, premib, premit)

    if pressure_regime == 0:
        if int(landfrac + 0.5) == 1 and snowh <= 0.000001:
            rhmin = rhminl - rhminl_adj_land
        else:
            rhmin = rhminl
    elif pressure_regime == 1:
        rhmin = rhminh
    else:
        rhwght = (premib - max(p, premit)) / (premib - premit)
        rhmin = rhminh * rhwght + rhminl * (1.0 - rhwght)

    dv = cldrh - rhmin
    if u >= 1.0:
        a = 1.0
        ga = 1.0e10
    elif u > (cldrh - dv / 6.0) and u < 1.0:
        a = 1.0 - (-3.0 / sqrt(2.0) * (u - cldrh) / dv) ** (2.0 / 3.0)
        ga = dv / sqrt(2.0) * sqrt(1.0 - a)
    elif u > (cldrh - dv) and u <= (cldrh - dv / 6.0):
        a = 4.0 * (
            cos(
                (1.0 / 3.0)
                * (
                    acos((3.0 / 2.0 / sqrt(2.0)) * (1.0 + (u - cldrh) / dv))
                    - 2.0 * 3.141592
                )
            )
        ) ** 2.0
        ga = dv / sqrt(2.0) * (1.0 / sqrt(a) - sqrt(a))
    else:
        a = 0.0
        ga = 1.0e10

    return (a, ga, rhmin)


@export
def cldfrc2m_astg_pdf_single_codon(
    u: float,
    p: float,
    qv: float,
    landfrac: float,
    snowh: float,
    rhminl: float,
    rhminl_adj_land: float,
    rhminh: float,
    premib: float,
    premit: float,
    a_p: cobj,
    ga_p: cobj,
    rhmin_p: cobj,
):
    a_out = Ptr[float](a_p)
    ga_out = Ptr[float](ga_p)
    rhmin_out = Ptr[float](rhmin_p)

    a, ga, rhmin = _cldfrc2m_astg_pdf_core(
        u,
        p,
        qv,
        landfrac,
        snowh,
        rhminl,
        rhminl_adj_land,
        rhminh,
        premib,
        premit,
    )
    a_out[0] = a
    ga_out[0] = ga
    rhmin_out[0] = rhmin


@export
def cldfrc2m_astg_pdf_codon(
    pcols: int,
    ncol: int,
    premib: float,
    premit: float,
    u_p: cobj,
    p_p: cobj,
    qv_p: cobj,
    landfrac_p: cobj,
    snowh_p: cobj,
    rhminl_p: cobj,
    rhminl_adj_land_p: cobj,
    rhminh_p: cobj,
    a_p: cobj,
    ga_p: cobj,
):
    u = Ptr[float](u_p)
    pressure = Ptr[float](p_p)
    qv = Ptr[float](qv_p)
    landfrac = Ptr[float](landfrac_p)
    snowh = Ptr[float](snowh_p)
    rhminl = Ptr[float](rhminl_p)
    rhminl_adj_land = Ptr[float](rhminl_adj_land_p)
    rhminh = Ptr[float](rhminh_p)
    a_out = Ptr[float](a_p)
    ga_out = Ptr[float](ga_p)

    for i0 in range(pcols):
        a_out[i0] = 0.0
        ga_out[i0] = 0.0

    for i0 in range(ncol):
        a, ga, _ = _cldfrc2m_astg_pdf_core(
            u[i0],
            pressure[i0],
            qv[i0],
            landfrac[i0],
            snowh[i0],
            rhminl[i0],
            rhminl_adj_land[i0],
            rhminh[i0],
            premib,
            premit,
        )
        a_out[i0] = a
        ga_out[i0] = ga


@inline
def _cldfrc2m_aist_core(
    iceopt: int,
    qv: float,
    t: float,
    p: float,
    qi: float,
    qsat: float,
    esl: float,
    esi: float,
    rhmaxi: float,
    rhmini: float,
    rhminl: float,
    rhminl_adj_land: float,
    rhminh: float,
    landfrac: float,
    snowh: float,
    ni: float,
    premib: float,
    premit: float,
    icecrit: float,
    rair: float,
) -> float:
    a = 26.87
    b = 0.569
    c = 0.002892
    as_fit = -68.4202
    bs = 0.983917
    cs = 2.81795
    kc = 75.0
    minice = 1.0e-12
    mincld = 1.0e-4
    qist_min = 1.0e-7
    qist_max = 5.0e-3

    aist = 0.0
    if iceopt < 3:
        if iceopt == 1:
            ttmp = max(195.0, min(t, 253.0)) - 273.16
            icicval = a + b * ttmp + c * ttmp**2.0
            rho = p / (rair * t)
            icicval = icicval * 1.0e-6 / rho
        else:
            ttmp = max(190.0, min(t, 273.16))
            icicval = 10.0 ** (as_fit * bs**ttmp + cs)
            icicval = icicval * 1.0e-6 * 18.0 / 28.97
        aist = max(0.0, min(qi / icicval, 1.0))
    elif iceopt == 3:
        aist = 1.0 - exp(-kc * qi / (qsat * (esi / esl)))
        aist = max(0.0, min(aist, 1.0))
    elif iceopt == 4:
        if p >= premib:
            if int(landfrac + 0.5) == 1 and snowh <= 0.000001:
                rhmin = rhminl - rhminl_adj_land
            else:
                rhmin = rhminl
        elif p < premit:
            rhmin = rhminh
        else:
            rhwght = (premib - max(p, premit)) / (premib - premit)
            rhmin = rhminh * rhwght + rhminl * (1.0 - rhwght)

        ncf = qi / ((1.0 - icecrit) * qsat)
        if ncf <= 0.0:
            aist = 0.0
        elif ncf > 0.0 and ncf <= 1.0 / 6.0:
            aist = 0.5 * (6.0 * ncf) ** (2.0 / 3.0)
        elif ncf > 1.0 / 6.0 and ncf < 1.0:
            phi = (acos(3.0 * (1.0 - ncf) / 2.0 ** (3.0 / 2.0)) + 4.0 * 3.1415927) / 3.0
            aist = 1.0 - 4.0 * cos(phi) * cos(phi)
        else:
            aist = 1.0
        aist = max(0.0, min(aist, 1.0))
    elif iceopt == 5:
        rhi = (qv + qi) / qsat * (esl / esi)
        rhdif = (rhi - rhmini) / (rhmaxi - rhmini)
        aist = min(1.0, max(rhdif, 0.0) ** 2)
    elif iceopt == 6:
        ah = 6.73834e-08
        bh = 0.0533110
        ch = 0.3493813
        rho = p / (rair * t)
        nil = ni * rho / 1000.0
        icicval = ah * exp(bh * t) * nil**ch
        icicval = icicval / rho / 1000.0
        aist = max(0.0, min(qi / icicval, 1.0))
        aist = min(aist, 1.0)

    if iceopt == 5 or iceopt == 6:
        if qi < minice:
            aist = 0.0
        else:
            aist = max(mincld, aist)

        if qi >= minice:
            icimr = qi / aist
            if icimr < qist_min:
                aist = max(0.0, min(1.0, qi / qist_min))
            if icimr > qist_max:
                aist = max(0.0, min(1.0, qi / qist_max))

    return max(0.0, min(aist, 0.999))


@export
def cldfrc2m_aist_single_codon(
    iceopt: int,
    qv: float,
    t: float,
    p: float,
    qi: float,
    landfrac: float,
    snowh: float,
    qsat: float,
    esl: float,
    esi: float,
    rhmaxi: float,
    rhmini: float,
    rhminl: float,
    rhminl_adj_land: float,
    rhminh: float,
    premib: float,
    premit: float,
    icecrit: float,
    rair: float,
) -> float:
    return _cldfrc2m_aist_core(
        iceopt,
        qv,
        t,
        p,
        qi,
        qsat,
        esl,
        esi,
        rhmaxi,
        rhmini,
        rhminl,
        rhminl_adj_land,
        rhminh,
        landfrac,
        snowh,
        0.0,
        premib,
        premit,
        icecrit,
        rair,
    )


@export
def cldfrc2m_aist_single_option5_codon(
    qv: float,
    qi: float,
    qsat: float,
    esl: float,
    esi: float,
    rhmini: float,
    rhmaxi: float,
) -> float:
    rhi = (qv + qi) / qsat * (esl / esi)
    rhdif = (rhi - rhmini) / (rhmaxi - rhmini)
    aist = min(1.0, max(rhdif, 0.0) ** 2)

    if qi < 1.0e-12:
        aist = 0.0
    else:
        aist = max(1.0e-4, aist)

    if qi >= 1.0e-12:
        icimr = qi / aist
        if icimr < 1.0e-7:
            aist = max(0.0, min(1.0, qi / 1.0e-7))
        if icimr > 5.0e-3:
            aist = max(0.0, min(1.0, qi / 5.0e-3))

    return max(0.0, min(aist, 0.999))


@export
def cldfrc2m_aist_vector_codon(
    pcols: int,
    ncol: int,
    iceopt: int,
    rhmaxi: float,
    premib: float,
    premit: float,
    icecrit: float,
    rair: float,
    qv_p: cobj,
    t_p: cobj,
    p_p: cobj,
    qi_p: cobj,
    ni_p: cobj,
    landfrac_p: cobj,
    snowh_p: cobj,
    qsat_p: cobj,
    esl_p: cobj,
    esi_p: cobj,
    rhmini_p: cobj,
    rhminl_p: cobj,
    rhminl_adj_land_p: cobj,
    rhminh_p: cobj,
    aist_p: cobj,
):
    qv = Ptr[float](qv_p)
    t = Ptr[float](t_p)
    p = Ptr[float](p_p)
    qi = Ptr[float](qi_p)
    ni = Ptr[float](ni_p)
    landfrac = Ptr[float](landfrac_p)
    snowh = Ptr[float](snowh_p)
    qsat = Ptr[float](qsat_p)
    esl = Ptr[float](esl_p)
    esi = Ptr[float](esi_p)
    rhmini = Ptr[float](rhmini_p)
    rhminl = Ptr[float](rhminl_p)
    rhminl_adj_land = Ptr[float](rhminl_adj_land_p)
    rhminh = Ptr[float](rhminh_p)
    aist_out = Ptr[float](aist_p)

    for i0 in range(pcols):
        aist_out[i0] = 0.0

    for i in range(1, ncol + 1):
        i0 = i - 1
        aist_out[i0] = _cldfrc2m_aist_core(
            iceopt,
            qv[i0],
            t[i0],
            p[i0],
            qi[i0],
            qsat[i0],
            esl[i0],
            esi[i0],
            rhmaxi,
            rhmini[i0],
            rhminl[i0],
            rhminl_adj_land[i0],
            rhminh[i0],
            landfrac[i0],
            snowh[i0],
            ni[i0],
            premib,
            premit,
            icecrit,
            rair,
        )


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


@inline
def _phys_prop_trimmed_len(text: Ptr[int], n: int) -> int:
    out = n
    while out > 0 and text[out - 1] == 32:
        out -= 1
    return out


@inline
def _phys_prop_trimmed_eq(a: Ptr[int], a_len: int, b: Ptr[int], b_len: int) -> bool:
    a_trim = _phys_prop_trimmed_len(a, a_len)
    b_trim = _phys_prop_trimmed_len(b, b_len)
    if a_trim != b_trim:
        return False

    i = 0
    while i < a_trim:
        if a[i] != b[i]:
            return False
        i += 1
    return True


@inline
def _physprop_name_ptr(base: Ptr[int], name_len: int, idx1: int) -> Ptr[int]:
    return base + (idx1 - 1) * name_len


@inline
def _physprop_is_known(name: Ptr[int], name_len: int, names: Ptr[int], known_count: int) -> bool:
    i = 1
    while i <= known_count:
        if _phys_prop_trimmed_eq(name, name_len, _physprop_name_ptr(names, name_len, i), name_len):
            return True
        i += 1
    return False


@inline
def _phys_prop_trimmed_eq_lit(text: Ptr[int], text_len: int, literal: str) -> bool:
    trim_len = _phys_prop_trimmed_len(text, text_len)
    if trim_len != len(literal):
        return False

    i = 0
    while i < trim_len:
        if text[i] != ord(literal[i]):
            return False
        i += 1
    return True


@export
def aerosol_optics_init_dispatch_codon(optics_len: int, optics_ascii_p: cobj) -> int:
    optics_ascii = Ptr[int](optics_ascii_p)
    if _phys_prop_trimmed_eq_lit(optics_ascii, optics_len, "zero"):
        return 1
    if _phys_prop_trimmed_eq_lit(optics_ascii, optics_len, "hygro"):
        return 2
    if _phys_prop_trimmed_eq_lit(optics_ascii, optics_len, "hygroscopic"):
        return 3
    if _phys_prop_trimmed_eq_lit(optics_ascii, optics_len, "nonhygro"):
        return 4
    if _phys_prop_trimmed_eq_lit(optics_ascii, optics_len, "insoluble"):
        return 5
    if _phys_prop_trimmed_eq_lit(optics_ascii, optics_len, "volcanic_radius"):
        return 6
    if _phys_prop_trimmed_eq_lit(optics_ascii, optics_len, "volcanic"):
        return 7
    if _phys_prop_trimmed_eq_lit(optics_ascii, optics_len, "modal"):
        return 8
    return -1


@export
def aerosol_optics_init_codon(optics_len: int, optics_ascii_p: cobj) -> int:
    return aerosol_optics_init_dispatch_codon(optics_len, optics_ascii_p)


@export
def bulk_props_init_is_sulfate_codon(name_len: int, name_ascii_p: cobj) -> int:
    name_ascii = Ptr[int](name_ascii_p)
    if _phys_prop_trimmed_eq_lit(name_ascii, name_len, "SULFATE"):
        return 1
    return 0


@export
def bulk_props_init_codon(name_len: int, name_ascii_p: cobj) -> int:
    return bulk_props_init_is_sulfate_codon(name_len, name_ascii_p)


@export
def refindex_aer_init_have_pair_codon(istat1: int, istat2: int, noerr: int) -> int:
    if istat1 == noerr and istat2 == noerr:
        return 1
    return 0


@export
def refindex_aer_init_codon(istat1: int, istat2: int, noerr: int) -> int:
    return refindex_aer_init_have_pair_codon(istat1, istat2, noerr)


@export
def insoluble_optics_init_dim_mask_codon(nbnd: int, nlwbands: int, swbands: int, nswbands: int) -> int:
    mask = 0
    if nbnd != nlwbands:
        mask |= 1
    if swbands != nswbands:
        mask |= 2
    return mask


@export
def insoluble_optics_init_codon(nbnd: int, nlwbands: int, swbands: int, nswbands: int) -> int:
    return insoluble_optics_init_dim_mask_codon(nbnd, nlwbands, swbands, nswbands)


@export
def hygroscopic_optics_init_dim_mask_codon(nbnd: int, nlwbands: int, swbands: int, nswbands: int) -> int:
    mask = 0
    if nbnd != nlwbands:
        mask |= 1
    if swbands != nswbands:
        mask |= 2
    return mask


@export
def hygroscopic_optics_init_codon(nbnd: int, nlwbands: int, swbands: int, nswbands: int) -> int:
    return hygroscopic_optics_init_dim_mask_codon(nbnd, nlwbands, swbands, nswbands)


@export
def modal_optics_init_dim_mask_codon(lw_val: int, nlwbands: int, sw_val: int, nswbands: int) -> int:
    mask = 0
    if lw_val != nlwbands:
        mask |= 1
    if sw_val != nswbands:
        mask |= 2
    return mask


@export
def modal_optics_init_codon(lw_val: int, nlwbands: int, sw_val: int, nswbands: int) -> int:
    return modal_optics_init_dim_mask_codon(lw_val, nlwbands, sw_val, nswbands)


@export
def physprop_init_file_order_codon(numphysprops: int, file_order_p: cobj):
    file_order = Ptr[int](file_order_p)
    i = 0
    while i < numphysprops:
        file_order[i] = i + 1
        i += 1


@export
def physprop_init_codon(numphysprops: int, file_order_p: cobj):
    physprop_init_file_order_codon(numphysprops, file_order_p)


@export
def modal_optics_init_copy_mode1_codon(
    ncoef: int,
    prefr: int,
    prefi: int,
    nbands: int,
    src_p: cobj,
    dst_p: cobj,
):
    src = Ptr[float](src_p)
    dst = Ptr[float](dst_p)

    b = 0
    while b < nbands:
        i3 = 0
        while i3 < prefi:
            i2 = 0
            while i2 < prefr:
                i1 = 0
                while i1 < ncoef:
                    idx4 = i1 + i2 * ncoef + i3 * ncoef * prefr + b * ncoef * prefr * prefi
                    idx5 = idx4
                    dst[idx4] = src[idx5]
                    i1 += 1
                i2 += 1
            i3 += 1
        b += 1


@export
def refindex_aer_init_fill_complex_codon(n: int, ref_real_p: cobj, ref_im_p: cobj, refindex_p: cobj):
    ref_real = Ptr[float](ref_real_p)
    ref_im = Ptr[float](ref_im_p)
    refindex = Ptr[float](refindex_p)

    i = 0
    while i < n:
        refindex[2 * i] = ref_real[i]
        refindex[2 * i + 1] = abs(ref_im[i])
        i += 1


@export
def physprop_get_check_id_codon(id_value: int, numphysprops: int) -> int:
    if id_value <= 0 or id_value > numphysprops:
        return 1
    return 0


@export
def physprop_get_codon(id_value: int, numphysprops: int) -> int:
    return physprop_get_check_id_codon(id_value, numphysprops)


@export
def physprop_accum_unique_files_codon(
    ncnst: int,
    name_len: int,
    numphysprops: int,
    radname_ascii_p: cobj,
    type_ascii_p: cobj,
    names_ascii_p: cobj,
    append_flags_p: cobj,
):
    radname_ascii = Ptr[int](radname_ascii_p)
    type_ascii = Ptr[int](type_ascii_p)
    names_ascii = Ptr[int](names_ascii_p)
    append_flags = Ptr[int](append_flags_p)

    known_count = numphysprops
    i = 1
    while i <= ncnst:
        append_flags[i - 1] = 0
        t = type_ascii[i - 1]
        if t == 65 or t == 77:
            name = _physprop_name_ptr(radname_ascii, name_len, i)
            if not _physprop_is_known(name, name_len, names_ascii, known_count):
                append_flags[i - 1] = 1
                known_count += 1
                dst = _physprop_name_ptr(names_ascii, name_len, known_count)
                j = 0
                while j < name_len:
                    dst[j] = name[j]
                    j += 1
        i += 1


@export
def physprop_get_id_codon(filename_len: int, filename_ascii_p: cobj, names_len: int, names_ascii_p: cobj,
                          numphysprops: int) -> int:
    filename_ascii = Ptr[int](filename_ascii_p)
    names_ascii = Ptr[int](names_ascii_p)

    iphysprop = 1
    while iphysprop <= numphysprops:
        name_offset = (iphysprop - 1) * names_len
        if _phys_prop_trimmed_eq(filename_ascii, filename_len, names_ascii + name_offset, names_len):
            return iphysprop
        iphysprop += 1
    return -1


@export
def exp_interpol_codon(n: int, x_p: cobj, f_p: cobj, y: float) -> float:
    x = Ptr[float](x_p)
    f = Ptr[float](f_p)

    k = _phys_prop_interp_index(n, x, y)
    k0 = k - 1
    k1 = k

    a = (log(f[k1] / f[k0])) / (x[k1] - x[k0])
    return f[k0] * exp(a * (y - x[k0]))


@export
def lin_interpol_codon(n: int, x_p: cobj, f_p: cobj, y: float) -> float:
    x = Ptr[float](x_p)
    f = Ptr[float](f_p)

    k = _phys_prop_interp_index(n, x, y)
    k0 = k - 1
    k1 = k

    a = (f[k1] - f[k0]) / (x[k1] - x[k0])
    return f[k0] + a * (y - x[k0])


@export
def phys_prop_exp_interpol_codon(n: int, x_p: cobj, f_p: cobj, y: float) -> float:
    return exp_interpol_codon(n, x_p, f_p, y)


@export
def phys_prop_lin_interpol_codon(n: int, x_p: cobj, f_p: cobj, y: float) -> float:
    return lin_interpol_codon(n, x_p, f_p, y)


@export
def aer_optics_log_rh_codon(
    nrh: int,
    nrh_test: int,
    ext_p: cobj,
    ssa_p: cobj,
    asm_p: cobj,
    rh_test_p: cobj,
    exti_p: cobj,
    ssai_p: cobj,
    asmi_p: cobj,
):
    ext = Ptr[float](ext_p)
    ssa = Ptr[float](ssa_p)
    asm = Ptr[float](asm_p)
    rh_test = Ptr[float](rh_test_p)
    exti = Ptr[float](exti_p)
    ssai = Ptr[float](ssai_p)
    asmi = Ptr[float](asmi_p)

    for krh_test in range(1, nrh_test + 1):
        value = (float(krh_test) - 1.0) / float(nrh_test - 1)
        rh_test[krh_test - 1] = sqrt(sqrt(sqrt(sqrt(value))))

    for krh_test in range(1, nrh_test + 1):
        rh = rh_test[krh_test - 1]
        krh = min(int(floor(rh * float(nrh))) + 1, nrh - 1)
        wrh = rh * float(nrh) - float(krh)
        exti[krh_test - 1] = ext[krh] * (wrh + 1.0) - ext[krh - 1] * wrh
        ssai[krh_test - 1] = ssa[krh] * (wrh + 1.0) - ssa[krh - 1] * wrh
        asmi[krh_test - 1] = asm[krh] * (wrh + 1.0) - asm[krh - 1] * wrh


@export
def modal_aer_opt_init_mode_dims_mismatch_codon(
    m_ncoef: int,
    m_prefr: int,
    m_prefi: int,
    ncoef: int,
    prefr: int,
    prefi: int,
) -> int:
    return _modal_aer_opt.modal_aer_opt_init_mode_dims_mismatch_codon(
        m_ncoef, m_prefr, m_prefi, ncoef, prefr, prefi
    )


@export
def modal_aer_opt_sw_init_state_codon(
    ncol: int,
    pcols: int,
    pver: int,
    nswbands: int,
    rga: float,
    rair: float,
    pdeldry_p: cobj,
    pmid_p: cobj,
    state_t_p: cobj,
    tauxar_p: cobj,
    wa_p: cobj,
    ga_p: cobj,
    fa_p: cobj,
    mass_p: cobj,
    air_density_p: cobj,
):
    _modal_aer_opt.modal_aer_opt_sw_init_state_codon(
        ncol,
        pcols,
        pver,
        nswbands,
        rga,
        rair,
        pdeldry_p,
        pmid_p,
        state_t_p,
        tauxar_p,
        wa_p,
        ga_p,
        fa_p,
        mass_p,
        air_density_p,
    )


@export
def modal_aer_opt_sw_zero_diagnostics_codon(
    ncol: int,
    pcols: int,
    pver: int,
    extinct_p: cobj,
    absorb_p: cobj,
    extinctuv_p: cobj,
    extinctnir_p: cobj,
    aodvis_p: cobj,
    aodvisst_p: cobj,
    aodabs_p: cobj,
    aodabsbc_p: cobj,
    ssavis_p: cobj,
    burdendust_p: cobj,
    burdenso4_p: cobj,
    burdenpom_p: cobj,
    burdensoa_p: cobj,
    burdenbc_p: cobj,
    burdenseasalt_p: cobj,
    dustaod_p: cobj,
    so4aod_p: cobj,
    pomaod_p: cobj,
    soaaod_p: cobj,
    bcaod_p: cobj,
    seasaltaod_p: cobj,
    aoduv_p: cobj,
    aodnir_p: cobj,
    aoduvst_p: cobj,
    aodnirst_p: cobj,
):
    _modal_aer_opt.modal_aer_opt_sw_zero_diagnostics_codon(
        ncol,
        pcols,
        pver,
        extinct_p,
        absorb_p,
        extinctuv_p,
        extinctnir_p,
        aodvis_p,
        aodvisst_p,
        aodabs_p,
        aodabsbc_p,
        ssavis_p,
        burdendust_p,
        burdenso4_p,
        burdenpom_p,
        burdensoa_p,
        burdenbc_p,
        burdenseasalt_p,
        dustaod_p,
        so4aod_p,
        pomaod_p,
        soaaod_p,
        bcaod_p,
        seasaltaod_p,
        aoduv_p,
        aodnir_p,
        aoduvst_p,
        aodnirst_p,
    )


@export
def modal_aer_opt_sw_mode_diag_init_codon(
    ncol: int,
    burden_p: cobj,
    aodmode_p: cobj,
    dustaodmode_p: cobj,
):
    _modal_aer_opt.modal_aer_opt_sw_mode_diag_init_codon(
        ncol,
        burden_p,
        aodmode_p,
        dustaodmode_p,
    )


@export
def modal_aer_opt_sw_reset_layer_codon(
    ncol: int,
    dryvol_p: cobj,
    dustvol_p: cobj,
    scatdust_p: cobj,
    absdust_p: cobj,
    hygrodust_p: cobj,
    scatso4_p: cobj,
    absso4_p: cobj,
    hygroso4_p: cobj,
    scatbc_p: cobj,
    absbc_p: cobj,
    hygrobc_p: cobj,
    scatpom_p: cobj,
    abspom_p: cobj,
    hygropom_p: cobj,
    scatsoa_p: cobj,
    abssoa_p: cobj,
    hygrosoa_p: cobj,
    scatseasalt_p: cobj,
    absseasalt_p: cobj,
    hygroseasalt_p: cobj,
    crefin_re_p: cobj,
    crefin_im_p: cobj,
):
    _modal_aer_opt.modal_aer_opt_sw_reset_layer_codon(
        ncol,
        dryvol_p,
        dustvol_p,
        scatdust_p,
        absdust_p,
        hygrodust_p,
        scatso4_p,
        absso4_p,
        hygroso4_p,
        scatbc_p,
        absbc_p,
        hygrobc_p,
        scatpom_p,
        abspom_p,
        hygropom_p,
        scatsoa_p,
        abssoa_p,
        hygrosoa_p,
        scatseasalt_p,
        absseasalt_p,
        hygroseasalt_p,
        crefin_re_p,
        crefin_im_p,
    )


@export
def modal_aer_opt_sw_species_volume_codon(
    ncol: int,
    pcols: int,
    k: int,
    specdens: float,
    specrefr: float,
    specrefi: float,
    specmmr_p: cobj,
    vol_p: cobj,
    dryvol_p: cobj,
    crefin_re_p: cobj,
    crefin_im_p: cobj,
):
    _modal_aer_opt.modal_aer_opt_sw_species_volume_codon(
        ncol,
        pcols,
        k,
        specdens,
        specrefr,
        specrefi,
        specmmr_p,
        vol_p,
        dryvol_p,
        crefin_re_p,
        crefin_im_p,
    )


@export
def modal_aer_opt_sw_water_volume_codon(
    ncol: int,
    pcols: int,
    k: int,
    rhoh2o: float,
    qaerwat_p: cobj,
    dryvol_p: cobj,
    watervol_p: cobj,
    wetvol_p: cobj,
):
    _modal_aer_opt.modal_aer_opt_sw_water_volume_codon(
        ncol,
        pcols,
        k,
        rhoh2o,
        qaerwat_p,
        dryvol_p,
        watervol_p,
        wetvol_p,
    )


@export
def modal_aer_opt_sw_has_negative_water_codon(ncol: int, watervol_p: cobj) -> int:
    return _modal_aer_opt.modal_aer_opt_sw_has_negative_water_codon(ncol, watervol_p)


@export
def modal_aer_opt_sw_water_refr_fastpath_codon(
    ncol: int,
    pcols: int,
    k: int,
    rhoh2o: float,
    crefwsw_re: float,
    crefwsw_im: float,
    qaerwat_p: cobj,
    dryvol_p: cobj,
    watervol_p: cobj,
    wetvol_p: cobj,
    crefin_re_p: cobj,
    crefin_im_p: cobj,
    refr_p: cobj,
    refi_p: cobj,
) -> int:
    return _modal_aer_opt.modal_aer_opt_sw_water_refr_fastpath_codon(
        ncol,
        pcols,
        k,
        rhoh2o,
        crefwsw_re,
        crefwsw_im,
        qaerwat_p,
        dryvol_p,
        watervol_p,
        wetvol_p,
        crefin_re_p,
        crefin_im_p,
        refr_p,
        refi_p,
    )


@export
def modal_aer_opt_sw_finalize_refr_codon(
    ncol: int,
    crefwsw_re: float,
    crefwsw_im: float,
    watervol_p: cobj,
    wetvol_p: cobj,
    crefin_re_p: cobj,
    crefin_im_p: cobj,
    refr_p: cobj,
    refi_p: cobj,
):
    _modal_aer_opt.modal_aer_opt_sw_finalize_refr_codon(
        ncol,
        crefwsw_re,
        crefwsw_im,
        watervol_p,
        wetvol_p,
        crefin_re_p,
        crefin_im_p,
        refr_p,
        refi_p,
    )


@export
def modal_aer_opt_sw_species_vis_diag_codon(
    ncol: int,
    pcols: int,
    k: int,
    spectype_code: int,
    specrefr: float,
    specrefi: float,
    hygro_aer: float,
    specmmr_p: cobj,
    mass_p: cobj,
    vol_p: cobj,
    burden_p: cobj,
    burdendust_p: cobj,
    burdenso4_p: cobj,
    burdenbc_p: cobj,
    burdenpom_p: cobj,
    burdensoa_p: cobj,
    burdenseasalt_p: cobj,
    dustvol_p: cobj,
    scatdust_p: cobj,
    absdust_p: cobj,
    hygrodust_p: cobj,
    scatso4_p: cobj,
    absso4_p: cobj,
    hygroso4_p: cobj,
    scatbc_p: cobj,
    absbc_p: cobj,
    hygrobc_p: cobj,
    scatpom_p: cobj,
    abspom_p: cobj,
    hygropom_p: cobj,
    scatsoa_p: cobj,
    abssoa_p: cobj,
    hygrosoa_p: cobj,
    scatseasalt_p: cobj,
    absseasalt_p: cobj,
    hygroseasalt_p: cobj,
):
    _modal_aer_opt.modal_aer_opt_sw_species_vis_diag_codon(
        ncol,
        pcols,
        k,
        spectype_code,
        specrefr,
        specrefi,
        hygro_aer,
        specmmr_p,
        mass_p,
        vol_p,
        burden_p,
        burdendust_p,
        burdenso4_p,
        burdenbc_p,
        burdenpom_p,
        burdensoa_p,
        burdenseasalt_p,
        dustvol_p,
        scatdust_p,
        absdust_p,
        hygrodust_p,
        scatso4_p,
        absso4_p,
        hygroso4_p,
        scatbc_p,
        absbc_p,
        hygrobc_p,
        scatpom_p,
        abspom_p,
        hygropom_p,
        scatsoa_p,
        abssoa_p,
        hygrosoa_p,
        scatseasalt_p,
        absseasalt_p,
        hygroseasalt_p,
    )


@export
def modal_aer_opt_sw_species_layer_batch_codon(
    ncol: int,
    pcols: int,
    k: int,
    spectype_code: int,
    do_vis: int,
    specdens: float,
    specrefr: float,
    specrefi: float,
    hygro_aer: float,
    specmmr_p: cobj,
    mass_p: cobj,
    vol_p: cobj,
    dryvol_p: cobj,
    crefin_re_p: cobj,
    crefin_im_p: cobj,
    burden_p: cobj,
    burdendust_p: cobj,
    burdenso4_p: cobj,
    burdenbc_p: cobj,
    burdenpom_p: cobj,
    burdensoa_p: cobj,
    burdenseasalt_p: cobj,
    dustvol_p: cobj,
    scatdust_p: cobj,
    absdust_p: cobj,
    hygrodust_p: cobj,
    scatso4_p: cobj,
    absso4_p: cobj,
    hygroso4_p: cobj,
    scatbc_p: cobj,
    absbc_p: cobj,
    hygrobc_p: cobj,
    scatpom_p: cobj,
    abspom_p: cobj,
    hygropom_p: cobj,
    scatsoa_p: cobj,
    abssoa_p: cobj,
    hygrosoa_p: cobj,
    scatseasalt_p: cobj,
    absseasalt_p: cobj,
    hygroseasalt_p: cobj,
):
    _modal_aer_opt.modal_aer_opt_sw_species_layer_batch_codon(
        ncol,
        pcols,
        k,
        spectype_code,
        do_vis,
        specdens,
        specrefr,
        specrefi,
        hygro_aer,
        specmmr_p,
        mass_p,
        vol_p,
        dryvol_p,
        crefin_re_p,
        crefin_im_p,
        burden_p,
        burdendust_p,
        burdenso4_p,
        burdenbc_p,
        burdenpom_p,
        burdensoa_p,
        burdenseasalt_p,
        dustvol_p,
        scatdust_p,
        absdust_p,
        hygrodust_p,
        scatso4_p,
        absso4_p,
        hygroso4_p,
        scatbc_p,
        absbc_p,
        hygrobc_p,
        scatpom_p,
        abspom_p,
        hygropom_p,
        scatsoa_p,
        abssoa_p,
        hygrosoa_p,
        scatseasalt_p,
        absseasalt_p,
        hygroseasalt_p,
    )


@export
def modal_aer_opt_sw_optics_props_codon(
    ncol: int,
    pcols: int,
    k: int,
    ncoef: int,
    xrmax: float,
    rhoh2o: float,
    radsurf_p: cobj,
    logradsurf_p: cobj,
    cheb_p: cobj,
    cext_p: cobj,
    cabs_p: cobj,
    casm_p: cobj,
    wetvol_p: cobj,
    mass_p: cobj,
    pext_p: cobj,
    specpext_p: cobj,
    pabs_p: cobj,
    pasm_p: cobj,
    palb_p: cobj,
    dopaer_p: cobj,
):
    _modal_aer_opt.modal_aer_opt_sw_optics_props_codon(
        ncol,
        pcols,
        k,
        ncoef,
        xrmax,
        rhoh2o,
        radsurf_p,
        logradsurf_p,
        cheb_p,
        cext_p,
        cabs_p,
        casm_p,
        wetvol_p,
        mass_p,
        pext_p,
        specpext_p,
        pabs_p,
        pasm_p,
        palb_p,
        dopaer_p,
    )


@export
def modal_aer_opt_sw_mode_diag_night_codon(
    nnite: int,
    fillvalue: float,
    idxnite_p: cobj,
    burden_p: cobj,
    aodmode_p: cobj,
    dustaodmode_p: cobj,
):
    _modal_aer_opt.modal_aer_opt_sw_mode_diag_night_codon(
        nnite,
        fillvalue,
        idxnite_p,
        burden_p,
        aodmode_p,
        dustaodmode_p,
    )


@export
def modal_aer_opt_sw_sum_diag_night_codon(
    nnite: int,
    pcols: int,
    pver: int,
    fillvalue: float,
    idxnite_p: cobj,
    extinct_p: cobj,
    absorb_p: cobj,
    aodvis_p: cobj,
    aodabs_p: cobj,
    aodvisst_p: cobj,
):
    _modal_aer_opt.modal_aer_opt_sw_sum_diag_night_codon(
        nnite,
        pcols,
        pver,
        fillvalue,
        idxnite_p,
        extinct_p,
        absorb_p,
        aodvis_p,
        aodabs_p,
        aodvisst_p,
    )


@export
def modal_aer_opt_sw_finalize_ssavis_codon(
    ncol: int,
    aodvis_p: cobj,
    ssavis_p: cobj,
):
    _modal_aer_opt.modal_aer_opt_sw_finalize_ssavis_codon(ncol, aodvis_p, ssavis_p)


@export
def modal_aer_opt_sw_climate_diag_night_codon(
    nnite: int,
    pcols: int,
    pver: int,
    fillvalue: float,
    idxnite_p: cobj,
    ssavis_p: cobj,
    aoduv_p: cobj,
    aodnir_p: cobj,
    aoduvst_p: cobj,
    aodnirst_p: cobj,
    extinctuv_p: cobj,
    extinctnir_p: cobj,
    burdendust_p: cobj,
    burdenso4_p: cobj,
    burdenpom_p: cobj,
    burdensoa_p: cobj,
    burdenbc_p: cobj,
    burdenseasalt_p: cobj,
    aodabsbc_p: cobj,
    dustaod_p: cobj,
    so4aod_p: cobj,
    pomaod_p: cobj,
    soaaod_p: cobj,
    bcaod_p: cobj,
    seasaltaod_p: cobj,
):
    _modal_aer_opt.modal_aer_opt_sw_climate_diag_night_codon(
        nnite,
        pcols,
        pver,
        fillvalue,
        idxnite_p,
        ssavis_p,
        aoduv_p,
        aodnir_p,
        aoduvst_p,
        aodnirst_p,
        extinctuv_p,
        extinctnir_p,
        burdendust_p,
        burdenso4_p,
        burdenpom_p,
        burdensoa_p,
        burdenbc_p,
        burdenseasalt_p,
        aodabsbc_p,
        dustaod_p,
        so4aod_p,
        pomaod_p,
        soaaod_p,
        bcaod_p,
        seasaltaod_p,
    )


@export
def modal_aer_opt_sw_accumulate_diagnostics_codon(
    ncol: int,
    pcols: int,
    k: int,
    do_uv: int,
    do_nir: int,
    do_vis: int,
    crefwsw_re: float,
    crefwsw_im: float,
    troplev_p: cobj,
    mass_p: cobj,
    air_density_p: cobj,
    dopaer_p: cobj,
    pabs_p: cobj,
    palb_p: cobj,
    wetvol_p: cobj,
    watervol_p: cobj,
    dustvol_p: cobj,
    scatdust_p: cobj,
    scatso4_p: cobj,
    scatbc_p: cobj,
    scatpom_p: cobj,
    scatsoa_p: cobj,
    scatseasalt_p: cobj,
    absdust_p: cobj,
    absso4_p: cobj,
    absbc_p: cobj,
    abspom_p: cobj,
    abssoa_p: cobj,
    absseasalt_p: cobj,
    hygrodust_p: cobj,
    hygroso4_p: cobj,
    hygrobc_p: cobj,
    hygropom_p: cobj,
    hygrosoa_p: cobj,
    hygroseasalt_p: cobj,
    extinctuv_p: cobj,
    aoduv_p: cobj,
    aoduvst_p: cobj,
    extinctnir_p: cobj,
    aodnir_p: cobj,
    aodnirst_p: cobj,
    extinct_p: cobj,
    absorb_p: cobj,
    aodvis_p: cobj,
    aodabs_p: cobj,
    aodmode_p: cobj,
    ssavis_p: cobj,
    aodvisst_p: cobj,
    dustaodmode_p: cobj,
    aodabsbc_p: cobj,
    dustaod_p: cobj,
    so4aod_p: cobj,
    pomaod_p: cobj,
    soaaod_p: cobj,
    bcaod_p: cobj,
    seasaltaod_p: cobj,
):
    _modal_aer_opt.modal_aer_opt_sw_accumulate_diagnostics_codon(
        ncol,
        pcols,
        k,
        do_uv,
        do_nir,
        do_vis,
        crefwsw_re,
        crefwsw_im,
        troplev_p,
        mass_p,
        air_density_p,
        dopaer_p,
        pabs_p,
        palb_p,
        wetvol_p,
        watervol_p,
        dustvol_p,
        scatdust_p,
        scatso4_p,
        scatbc_p,
        scatpom_p,
        scatsoa_p,
        scatseasalt_p,
        absdust_p,
        absso4_p,
        absbc_p,
        abspom_p,
        abssoa_p,
        absseasalt_p,
        hygrodust_p,
        hygroso4_p,
        hygrobc_p,
        hygropom_p,
        hygrosoa_p,
        hygroseasalt_p,
        extinctuv_p,
        aoduv_p,
        aoduvst_p,
        extinctnir_p,
        aodnir_p,
        aodnirst_p,
        extinct_p,
        absorb_p,
        aodvis_p,
        aodabs_p,
        aodmode_p,
        ssavis_p,
        aodvisst_p,
        dustaodmode_p,
        aodabsbc_p,
        dustaod_p,
        so4aod_p,
        pomaod_p,
        soaaod_p,
        bcaod_p,
        seasaltaod_p,
    )


@export
def modal_aer_opt_sw_accumulate_tau_codon(
    ncol: int,
    pcols: int,
    pver: int,
    k: int,
    isw: int,
    dopaer_p: cobj,
    palb_p: cobj,
    pasm_p: cobj,
    tauxar_p: cobj,
    wa_p: cobj,
    ga_p: cobj,
    fa_p: cobj,
):
    _modal_aer_opt.modal_aer_opt_sw_accumulate_tau_codon(
        ncol,
        pcols,
        pver,
        k,
        isw,
        dopaer_p,
        palb_p,
        pasm_p,
        tauxar_p,
        wa_p,
        ga_p,
        fa_p,
    )


@export
def modal_aer_opt_sw_has_bad_dopaer_codon(ncol: int, dopaer_p: cobj) -> int:
    return _modal_aer_opt.modal_aer_opt_sw_has_bad_dopaer_codon(ncol, dopaer_p)


@export
def modal_aer_opt_sw_optics_diag_tau_batch_codon(
    ncol: int,
    pcols: int,
    pver: int,
    k: int,
    isw: int,
    ncoef: int,
    xrmax: float,
    rhoh2o: float,
    do_uv: int,
    do_nir: int,
    do_vis: int,
    crefwsw_re: float,
    crefwsw_im: float,
    radsurf_p: cobj,
    logradsurf_p: cobj,
    cheb_p: cobj,
    cext_p: cobj,
    cabs_p: cobj,
    casm_p: cobj,
    wetvol_p: cobj,
    mass_p: cobj,
    pext_p: cobj,
    specpext_p: cobj,
    pabs_p: cobj,
    pasm_p: cobj,
    palb_p: cobj,
    dopaer_p: cobj,
    troplev_p: cobj,
    air_density_p: cobj,
    watervol_p: cobj,
    dustvol_p: cobj,
    scatdust_p: cobj,
    scatso4_p: cobj,
    scatbc_p: cobj,
    scatpom_p: cobj,
    scatsoa_p: cobj,
    scatseasalt_p: cobj,
    absdust_p: cobj,
    absso4_p: cobj,
    absbc_p: cobj,
    abspom_p: cobj,
    abssoa_p: cobj,
    absseasalt_p: cobj,
    hygrodust_p: cobj,
    hygroso4_p: cobj,
    hygrobc_p: cobj,
    hygropom_p: cobj,
    hygrosoa_p: cobj,
    hygroseasalt_p: cobj,
    extinctuv_p: cobj,
    aoduv_p: cobj,
    aoduvst_p: cobj,
    extinctnir_p: cobj,
    aodnir_p: cobj,
    aodnirst_p: cobj,
    extinct_p: cobj,
    absorb_p: cobj,
    aodvis_p: cobj,
    aodabs_p: cobj,
    aodmode_p: cobj,
    ssavis_p: cobj,
    aodvisst_p: cobj,
    dustaodmode_p: cobj,
    aodabsbc_p: cobj,
    dustaod_p: cobj,
    so4aod_p: cobj,
    pomaod_p: cobj,
    soaaod_p: cobj,
    bcaod_p: cobj,
    seasaltaod_p: cobj,
    tauxar_p: cobj,
    wa_p: cobj,
    ga_p: cobj,
    fa_p: cobj,
) -> int:
    return _modal_aer_opt.modal_aer_opt_sw_optics_diag_tau_batch_codon(
        ncol,
        pcols,
        pver,
        k,
        isw,
        ncoef,
        xrmax,
        rhoh2o,
        do_uv,
        do_nir,
        do_vis,
        crefwsw_re,
        crefwsw_im,
        radsurf_p,
        logradsurf_p,
        cheb_p,
        cext_p,
        cabs_p,
        casm_p,
        wetvol_p,
        mass_p,
        pext_p,
        specpext_p,
        pabs_p,
        pasm_p,
        palb_p,
        dopaer_p,
        troplev_p,
        air_density_p,
        watervol_p,
        dustvol_p,
        scatdust_p,
        scatso4_p,
        scatbc_p,
        scatpom_p,
        scatsoa_p,
        scatseasalt_p,
        absdust_p,
        absso4_p,
        absbc_p,
        abspom_p,
        abssoa_p,
        absseasalt_p,
        hygrodust_p,
        hygroso4_p,
        hygrobc_p,
        hygropom_p,
        hygrosoa_p,
        hygroseasalt_p,
        extinctuv_p,
        aoduv_p,
        aoduvst_p,
        extinctnir_p,
        aodnir_p,
        aodnirst_p,
        extinct_p,
        absorb_p,
        aodvis_p,
        aodabs_p,
        aodmode_p,
        ssavis_p,
        aodvisst_p,
        dustaodmode_p,
        aodabsbc_p,
        dustaod_p,
        so4aod_p,
        pomaod_p,
        soaaod_p,
        bcaod_p,
        seasaltaod_p,
        tauxar_p,
        wa_p,
        ga_p,
        fa_p,
    )


@export
def modal_aer_opt_lw_init_state_codon(
    ncol: int,
    pcols: int,
    pver: int,
    nlwbands: int,
    rga: float,
    pdeldry_p: cobj,
    tauxar_p: cobj,
    mass_p: cobj,
):
    _modal_aer_opt.modal_aer_opt_lw_init_state_codon(
        ncol,
        pcols,
        pver,
        nlwbands,
        rga,
        pdeldry_p,
        tauxar_p,
        mass_p,
    )


@export
def modal_aer_opt_lw_optics_props_codon(
    ncol: int,
    pcols: int,
    pver: int,
    k: int,
    ncoef: int,
    rhoh2o: float,
    cheby_p: cobj,
    cabs_p: cobj,
    wetvol_p: cobj,
    mass_p: cobj,
    pabs_p: cobj,
    dopaer_p: cobj,
):
    _modal_aer_opt.modal_aer_opt_lw_optics_props_codon(
        ncol,
        pcols,
        pver,
        k,
        ncoef,
        rhoh2o,
        cheby_p,
        cabs_p,
        wetvol_p,
        mass_p,
        pabs_p,
        dopaer_p,
    )


@export
def modal_aer_opt_lw_accumulate_tau_codon(
    ncol: int,
    pcols: int,
    pver: int,
    k: int,
    ilw: int,
    dopaer_p: cobj,
    tauxar_p: cobj,
):
    _modal_aer_opt.modal_aer_opt_lw_accumulate_tau_codon(
        ncol,
        pcols,
        pver,
        k,
        ilw,
        dopaer_p,
        tauxar_p,
    )


@export
def modal_aero_lw_codon(
    stage: int,
    ncol: int,
    pcols: int,
    pver: int,
    k: int,
    ilw: int,
    ncoef: int,
    rhoh2o: float,
    crefwlw_re: float,
    crefwlw_im: float,
    specdens: float,
    specref_re: float,
    specref_im: float,
    specmmr_p: cobj,
    qaerwat_p: cobj,
    cheby_p: cobj,
    cabs_p: cobj,
    mass_p: cobj,
    vol_p: cobj,
    dryvol_p: cobj,
    watervol_p: cobj,
    wetvol_p: cobj,
    crefin_re_p: cobj,
    crefin_im_p: cobj,
    refr_p: cobj,
    refi_p: cobj,
    pabs_p: cobj,
    dopaer_p: cobj,
    tauxar_p: cobj,
):
    _modal_aer_opt.modal_aero_lw_codon(
        stage,
        ncol,
        pcols,
        pver,
        k,
        ilw,
        ncoef,
        rhoh2o,
        crefwlw_re,
        crefwlw_im,
        specdens,
        specref_re,
        specref_im,
        specmmr_p,
        qaerwat_p,
        cheby_p,
        cabs_p,
        mass_p,
        vol_p,
        dryvol_p,
        watervol_p,
        wetvol_p,
        crefin_re_p,
        crefin_im_p,
        refr_p,
        refi_p,
        pabs_p,
        dopaer_p,
        tauxar_p,
    )


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
def modal_aer_opt_sw_binterp3_codon(
    pcols: int,
    ncol: int,
    ncoef: int,
    prefr: int,
    prefi: int,
    extpsw_p: cobj,
    abspsw_p: cobj,
    asmpsw_p: cobj,
    refr_p: cobj,
    refi_p: cobj,
    refrtabsw_p: cobj,
    refitabsw_p: cobj,
    itab_p: cobj,
    jtab_p: cobj,
    ttab_p: cobj,
    utab_p: cobj,
    cext_p: cobj,
    cabs_p: cobj,
    casm_p: cobj,
):
    _modal_aer_opt.modal_aer_opt_sw_binterp3_codon(
        pcols,
        ncol,
        ncoef,
        prefr,
        prefi,
        extpsw_p,
        abspsw_p,
        asmpsw_p,
        refr_p,
        refi_p,
        refrtabsw_p,
        refitabsw_p,
        itab_p,
        jtab_p,
        ttab_p,
        utab_p,
        cext_p,
        cabs_p,
        casm_p,
    )


@export
def ndrop_init_scalars_codon(
    mwh2o: float,
    r_universal: float,
    rhoh2o: float,
    pi: float,
    scalars_p: cobj,
):
    _ndrop.ndrop_init_scalars_codon(mwh2o, r_universal, rhoh2o, pi, scalars_p)


@export
def ndrop_init_counts_codon(
    nmode: int,
    nspec_amode_p: cobj,
    nspec_max_p: cobj,
    ncnst_tot_p: cobj,
):
    _ndrop.ndrop_init_counts_codon(nmode, nspec_amode_p, nspec_max_p, ncnst_tot_p)


@export
def ndrop_init_mam_idx_codon(
    nmode: int,
    nspec_max: int,
    nspec_amode_p: cobj,
    mam_idx_p: cobj,
):
    _ndrop.ndrop_init_mam_idx_codon(nmode, nspec_max, nspec_amode_p, mam_idx_p)


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
def ndrop_dropmixnuc_zero_fields_codon(
    pcols: int,
    pver: int,
    ntot_amode: int,
    factnum_p: cobj,
    wtke_p: cobj,
):
    _ndrop.ndrop_dropmixnuc_zero_fields_codon(pcols, pver, ntot_amode, factnum_p, wtke_p)


@export
def ndrop_dropmixnuc_column_init_codon(
    i: int,
    pcols: int,
    pver: int,
    top_lev: int,
    ntot_amode: int,
    gravit: float,
    rair: float,
    zkmin: float,
    zkmax: float,
    wmixmin: float,
    ncldwtr_p: cobj,
    temp_p: cobj,
    pmid_p: cobj,
    pint_p: cobj,
    rpdel_p: cobj,
    zm_p: cobj,
    kvh_p: cobj,
    wsub_p: cobj,
    qcld_p: cobj,
    qncld_p: cobj,
    srcn_p: cobj,
    cs_p: cobj,
    dz_p: cobj,
    nact_p: cobj,
    mact_p: cobj,
    zn_p: cobj,
    ekd_p: cobj,
    csbot_p: cobj,
    csbot_cscen_p: cobj,
    wtke_cen_p: cobj,
    wtke_p: cobj,
    nsource_p: cobj,
    zs_p: cobj,
):
    _ndrop.ndrop_dropmixnuc_column_init_codon(
        i,
        pcols,
        pver,
        top_lev,
        ntot_amode,
        gravit,
        rair,
        zkmin,
        zkmax,
        wmixmin,
        ncldwtr_p,
        temp_p,
        pmid_p,
        pint_p,
        rpdel_p,
        zm_p,
        kvh_p,
        wsub_p,
        qcld_p,
        qncld_p,
        srcn_p,
        cs_p,
        dz_p,
        nact_p,
        mact_p,
        zn_p,
        ekd_p,
        csbot_p,
        csbot_cscen_p,
        wtke_cen_p,
        wtke_p,
        nsource_p,
        zs_p,
    )


@export
def ndrop_dropmixnuc_mix_setup_codon(
    i: int,
    pcols: int,
    pver: int,
    top_lev: int,
    ntot_amode: int,
    dtmicro: float,
    taumix_internal_pver_inv: float,
    cldn_p: cobj,
    zs_p: cobj,
    zn_p: cobj,
    csbot_p: cobj,
    ekd_p: cobj,
    nact_p: cobj,
    mact_p: cobj,
    ekk0_p: cobj,
    ekkp_p: cobj,
    ekkm_p: cobj,
    overlapp_p: cobj,
    overlapm_p: cobj,
    count_submix_p: cobj,
    nsubmix_p: cobj,
    dtmix_p: cobj,
):
    _ndrop.ndrop_dropmixnuc_mix_setup_codon(
        i,
        pcols,
        pver,
        top_lev,
        ntot_amode,
        dtmicro,
        taumix_internal_pver_inv,
        cldn_p,
        zs_p,
        zn_p,
        csbot_p,
        ekd_p,
        nact_p,
        mact_p,
        ekk0_p,
        ekkp_p,
        ekkm_p,
        overlapp_p,
        overlapm_p,
        count_submix_p,
        nsubmix_p,
        dtmix_p,
    )


@export
def ndrop_dropmixnuc_aero_column_copy_codon(
    i: int,
    pcols: int,
    pver: int,
    top_lev: int,
    ncnst_tot: int,
    mm: int,
    slot: int,
    zero_all: int,
    raer_fld_p: cobj,
    qqcw_fld_p: cobj,
    raercol_p: cobj,
    raercol_cw_p: cobj,
):
    _ndrop.ndrop_dropmixnuc_aero_column_copy_codon(
        i,
        pcols,
        pver,
        top_lev,
        ncnst_tot,
        mm,
        slot,
        zero_all,
        raer_fld_p,
        qqcw_fld_p,
        raercol_p,
        raercol_cw_p,
    )


@export
def ndrop_dropmixnuc_aero_column_copy_all_codon(
    i: int,
    pcols: int,
    pver: int,
    top_lev: int,
    ntot_amode: int,
    ncnst_tot: int,
    slot: int,
    raer_ptrs_p: cobj,
    qqcw_ptrs_p: cobj,
    nspec_amode_p: cobj,
    mam_idx_p: cobj,
    raercol_p: cobj,
    raercol_cw_p: cobj,
):
    _ndrop.ndrop_dropmixnuc_aero_column_copy_all_codon(
        i,
        pcols,
        pver,
        top_lev,
        ntot_amode,
        ncnst_tot,
        slot,
        raer_ptrs_p,
        qqcw_ptrs_p,
        nspec_amode_p,
        mam_idx_p,
        raercol_p,
        raercol_cw_p,
    )


@export
def ndrop_dropmixnuc_aero_tend_prepare_codon(
    i: int,
    pcols: int,
    pver: int,
    top_lev: int,
    ncnst_tot: int,
    mm: int,
    slot: int,
    dtinv: float,
    raer_fld_p: cobj,
    qqcw_fld_p: cobj,
    raercol_p: cobj,
    raercol_cw_p: cobj,
    raertend_p: cobj,
    qqcwtend_p: cobj,
):
    _ndrop.ndrop_dropmixnuc_aero_tend_prepare_codon(
        i,
        pcols,
        pver,
        top_lev,
        ncnst_tot,
        mm,
        slot,
        dtinv,
        raer_fld_p,
        qqcw_fld_p,
        raercol_p,
        raercol_cw_p,
        raertend_p,
        qqcwtend_p,
    )


@export
def ndrop_dropmixnuc_aero_tend_commit_qqcw_codon(
    i: int,
    pcols: int,
    pver: int,
    top_lev: int,
    ncnst_tot: int,
    mm: int,
    slot: int,
    qqcw_fld_p: cobj,
    raercol_cw_p: cobj,
):
    _ndrop.ndrop_dropmixnuc_aero_tend_commit_qqcw_codon(
        i,
        pcols,
        pver,
        top_lev,
        ncnst_tot,
        mm,
        slot,
        qqcw_fld_p,
        raercol_cw_p,
    )


@export
def ndrop_dropmixnuc_aero_tend_commit_ptend_codon(
    i: int,
    psetcols: int,
    pver: int,
    top_lev: int,
    pcnst: int,
    lptr: int,
    raertend_p: cobj,
    ptend_q_p: cobj,
):
    _ndrop.ndrop_dropmixnuc_aero_tend_commit_ptend_codon(
        i,
        psetcols,
        pver,
        top_lev,
        pcnst,
        lptr,
        raertend_p,
        ptend_q_p,
    )


@export
def ndrop_dropmixnuc_aero_coltend_codon(
    i: int,
    pcols: int,
    pver: int,
    mm: int,
    gravit: float,
    pdel_p: cobj,
    raertend_p: cobj,
    qqcwtend_p: cobj,
    coltend_out_p: cobj,
    coltend_cw_out_p: cobj,
):
    _ndrop.ndrop_dropmixnuc_aero_coltend_codon(
        i,
        pcols,
        pver,
        mm,
        gravit,
        pdel_p,
        raertend_p,
        qqcwtend_p,
        coltend_out_p,
        coltend_cw_out_p,
    )


@export
def ndrop_dropmixnuc_aero_tend_all_codon(
    i: int,
    pcols: int,
    psetcols: int,
    pver: int,
    top_lev: int,
    ntot_amode: int,
    ncnst_tot: int,
    slot: int,
    dtinv: float,
    gravit: float,
    raer_ptrs_p: cobj,
    qqcw_ptrs_p: cobj,
    nspec_amode_p: cobj,
    mam_idx_p: cobj,
    mam_cnst_idx_p: cobj,
    pdel_p: cobj,
    raercol_p: cobj,
    raercol_cw_p: cobj,
    coltend_p: cobj,
    coltend_cw_p: cobj,
    ptend_q_p: cobj,
):
    _ndrop.ndrop_dropmixnuc_aero_tend_all_codon(
        i,
        pcols,
        psetcols,
        pver,
        top_lev,
        ntot_amode,
        ncnst_tot,
        slot,
        dtinv,
        gravit,
        raer_ptrs_p,
        qqcw_ptrs_p,
        nspec_amode_p,
        mam_idx_p,
        mam_cnst_idx_p,
        pdel_p,
        raercol_p,
        raercol_cw_p,
        coltend_p,
        coltend_cw_p,
        ptend_q_p,
    )


@export
def ndrop_dropmixnuc_finalize_column_codon(
    i: int,
    pcols: int,
    pver: int,
    top_lev: int,
    dtinv: float,
    gravit: float,
    qcld_p: cobj,
    ncldwtr_p: cobj,
    pdel_p: cobj,
    nsource_p: cobj,
    ndropmix_p: cobj,
    tendnd_p: cobj,
    ndropcol_p: cobj,
):
    _ndrop.ndrop_dropmixnuc_finalize_column_codon(
        i,
        pcols,
        pver,
        top_lev,
        dtinv,
        gravit,
        qcld_p,
        ncldwtr_p,
        pdel_p,
        nsource_p,
        ndropmix_p,
        tendnd_p,
        ndropcol_p,
    )


@export
def ndrop_dropmixnuc_clear_old_cloud_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    ntot_amode: int,
    ncnst_tot: int,
    nsav: int,
    dtinv: float,
    qcld_p: cobj,
    nsource_p: cobj,
    nspec_amode_p: cobj,
    mam_idx_p: cobj,
    raercol_p: cobj,
    raercol_cw_p: cobj,
):
    _ndrop.ndrop_dropmixnuc_clear_old_cloud_codon(
        i,
        k,
        pcols,
        pver,
        ntot_amode,
        ncnst_tot,
        nsav,
        dtinv,
        qcld_p,
        nsource_p,
        nspec_amode_p,
        mam_idx_p,
        raercol_p,
        raercol_cw_p,
    )


@export
def ndrop_dropmixnuc_factnum_store_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    ntot_amode: int,
    fn_p: cobj,
    factnum_p: cobj,
):
    _ndrop.ndrop_dropmixnuc_factnum_store_codon(i, k, pcols, pver, ntot_amode, fn_p, factnum_p)


@export
def ndrop_dropmixnuc_shrink_cloud_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    ntot_amode: int,
    ncnst_tot: int,
    nsav: int,
    dtinv: float,
    cldn_tmp: float,
    cldo_tmp: float,
    qcld_p: cobj,
    nsource_p: cobj,
    nspec_amode_p: cobj,
    mam_idx_p: cobj,
    raercol_p: cobj,
    raercol_cw_p: cobj,
):
    _ndrop.ndrop_dropmixnuc_shrink_cloud_codon(
        i,
        k,
        pcols,
        pver,
        ntot_amode,
        ncnst_tot,
        nsav,
        dtinv,
        cldn_tmp,
        cldo_tmp,
        qcld_p,
        nsource_p,
        nspec_amode_p,
        mam_idx_p,
        raercol_p,
        raercol_cw_p,
    )


@export
def ndrop_dropmixnuc_grow_cloud_number_update_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    ncnst_tot: int,
    nsav: int,
    mm: int,
    dtinv: float,
    dumc: float,
    fn_m: float,
    raer_fld_p: cobj,
    qcld_p: cobj,
    nsource_p: cobj,
    raercol_p: cobj,
    raercol_cw_p: cobj,
):
    _ndrop.ndrop_dropmixnuc_grow_cloud_number_update_codon(
        i,
        k,
        pcols,
        pver,
        ncnst_tot,
        nsav,
        mm,
        dtinv,
        dumc,
        fn_m,
        raer_fld_p,
        qcld_p,
        nsource_p,
        raercol_p,
        raercol_cw_p,
    )


@export
def ndrop_dropmixnuc_grow_cloud_update_all_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    ntot_amode: int,
    ncnst_tot: int,
    nsav: int,
    dtinv: float,
    dumc: float,
    raer_ptrs_p: cobj,
    nspec_amode_p: cobj,
    mam_idx_p: cobj,
    fn_p: cobj,
    fm_p: cobj,
    qcld_p: cobj,
    nsource_p: cobj,
    raercol_p: cobj,
    raercol_cw_p: cobj,
    factnum_p: cobj,
):
    _ndrop.ndrop_dropmixnuc_grow_cloud_update_all_codon(
        i,
        k,
        pcols,
        pver,
        ntot_amode,
        ncnst_tot,
        nsav,
        dtinv,
        dumc,
        raer_ptrs_p,
        nspec_amode_p,
        mam_idx_p,
        fn_p,
        fm_p,
        qcld_p,
        nsource_p,
        raercol_p,
        raercol_cw_p,
        factnum_p,
    )


@export
def ndrop_dropmixnuc_grow_cloud_species_update_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    ncnst_tot: int,
    nsav: int,
    mm: int,
    dum: float,
    raer_fld_p: cobj,
    raercol_p: cobj,
    raercol_cw_p: cobj,
):
    _ndrop.ndrop_dropmixnuc_grow_cloud_species_update_codon(
        i,
        k,
        pcols,
        pver,
        ncnst_tot,
        nsav,
        mm,
        dum,
        raer_fld_p,
        raercol_p,
        raercol_cw_p,
    )


@export
def ndrop_dropmixnuc_old_cloud_activate_update_codon(
    i: int,
    k: int,
    kp1: int,
    pcols: int,
    pver: int,
    ntot_amode: int,
    ncnst_tot: int,
    nsav: int,
    dumc: float,
    dum: float,
    cs_ik: float,
    dz_ik: float,
    taumix_internal_pver_inv: float,
    fluxn_p: cobj,
    fluxm_p: cobj,
    nact_p: cobj,
    mact_p: cobj,
    mam_idx_p: cobj,
    raercol_p: cobj,
    raercol_cw_p: cobj,
    srcn_p: cobj,
    nsource_p: cobj,
):
    _ndrop.ndrop_dropmixnuc_old_cloud_activate_update_codon(
        i,
        k,
        kp1,
        pcols,
        pver,
        ntot_amode,
        ncnst_tot,
        nsav,
        dumc,
        dum,
        cs_ik,
        dz_ik,
        taumix_internal_pver_inv,
        fluxn_p,
        fluxm_p,
        nact_p,
        mact_p,
        mam_idx_p,
        raercol_p,
        raercol_cw_p,
        srcn_p,
        nsource_p,
    )


@export
def ndrop_dropmixnuc_srcn_from_nact_codon(
    pver: int,
    top_lev: int,
    ntot_amode: int,
    ncnst_tot: int,
    nsav: int,
    taumix_internal_pver_inv: float,
    nact_p: cobj,
    mam_idx_p: cobj,
    raercol_p: cobj,
    raercol_cw_p: cobj,
    srcn_p: cobj,
):
    _ndrop.ndrop_dropmixnuc_srcn_from_nact_codon(
        pver,
        top_lev,
        ntot_amode,
        ncnst_tot,
        nsav,
        taumix_internal_pver_inv,
        nact_p,
        mam_idx_p,
        raercol_p,
        raercol_cw_p,
        srcn_p,
    )


@export
def ndrop_dropmixnuc_source_from_act_codon(
    pver: int,
    top_lev: int,
    ncnst_tot: int,
    m: int,
    mm: int,
    nsav: int,
    taumix_internal_pver_inv: float,
    act_p: cobj,
    raercol_p: cobj,
    raercol_cw_p: cobj,
    source_p: cobj,
):
    _ndrop.ndrop_dropmixnuc_source_from_act_codon(
        pver,
        top_lev,
        ncnst_tot,
        m,
        mm,
        nsav,
        taumix_internal_pver_inv,
        act_p,
        raercol_p,
        raercol_cw_p,
        source_p,
    )


@export
def ndrop_dropmixnuc_evaporate_clear_layers_codon(
    i: int,
    pcols: int,
    pver: int,
    top_lev: int,
    ntot_amode: int,
    ncnst_tot: int,
    nnew: int,
    cldn_p: cobj,
    qcld_p: cobj,
    nspec_amode_p: cobj,
    mam_idx_p: cobj,
    raercol_p: cobj,
    raercol_cw_p: cobj,
):
    _ndrop.ndrop_dropmixnuc_evaporate_clear_layers_codon(
        i,
        pcols,
        pver,
        top_lev,
        ntot_amode,
        ncnst_tot,
        nnew,
        cldn_p,
        qcld_p,
        nspec_amode_p,
        mam_idx_p,
        raercol_p,
        raercol_cw_p,
    )


@export
def ndrop_dropmixnuc_swap_slots_codon(nsav_p: cobj, nnew_p: cobj):
    _ndrop.ndrop_dropmixnuc_swap_slots_codon(nsav_p, nnew_p)


@export
def ndrop_dropmixnuc_submix_iter_init_codon(
    pver: int,
    qcld_p: cobj,
    qncld_p: cobj,
    srcn_p: cobj,
    nsav_p: cobj,
    nnew_p: cobj,
):
    _ndrop.ndrop_dropmixnuc_submix_iter_init_codon(
        pver,
        qcld_p,
        qncld_p,
        srcn_p,
        nsav_p,
        nnew_p,
    )


@export
def ndrop_dropmixnuc_zero_tendencies_codon(
    pver: int,
    raertend_p: cobj,
    qqcwtend_p: cobj,
):
    _ndrop.ndrop_dropmixnuc_zero_tendencies_codon(pver, raertend_p, qqcwtend_p)


@export
def ndrop_loadaer_zero_codon(
    istart: int,
    istop: int,
    vaerosol_p: cobj,
    hygro_p: cobj,
):
    _ndrop.ndrop_loadaer_zero_codon(istart, istop, vaerosol_p, hygro_p)


@export
def ndrop_loadaer_species_accum_codon(
    istart: int,
    istop: int,
    k: int,
    pcols: int,
    phase: int,
    specdens: float,
    spechygro: float,
    raer_p: cobj,
    qqcw_p: cobj,
    vaerosol_p: cobj,
    hygro_p: cobj,
):
    _ndrop.ndrop_loadaer_species_accum_codon(
        istart,
        istop,
        k,
        pcols,
        phase,
        specdens,
        spechygro,
        raer_p,
        qqcw_p,
        vaerosol_p,
        hygro_p,
    )


@export
def ndrop_loadaer_species_batch_codon(
    istart: int,
    istop: int,
    k: int,
    pcols: int,
    nspec: int,
    phase: int,
    raer_ptrs_p: cobj,
    qqcw_ptrs_p: cobj,
    specdens_p: cobj,
    spechygro_p: cobj,
    vaerosol_p: cobj,
    hygro_p: cobj,
):
    _ndrop.ndrop_loadaer_species_batch_codon(
        istart,
        istop,
        k,
        pcols,
        nspec,
        phase,
        raer_ptrs_p,
        qqcw_ptrs_p,
        specdens_p,
        spechygro_p,
        vaerosol_p,
        hygro_p,
    )


@export
def ndrop_loadaer_finalize_volume_codon(
    istart: int,
    istop: int,
    k: int,
    pcols: int,
    cs_p: cobj,
    vaerosol_p: cobj,
    hygro_p: cobj,
):
    _ndrop.ndrop_loadaer_finalize_volume_codon(istart, istop, k, pcols, cs_p, vaerosol_p, hygro_p)


@export
def ndrop_loadaer_number_codon(
    istart: int,
    istop: int,
    k: int,
    pcols: int,
    phase: int,
    voltonumblo: float,
    voltonumbhi: float,
    raer_p: cobj,
    qqcw_p: cobj,
    cs_p: cobj,
    vaerosol_p: cobj,
    naerosol_p: cobj,
):
    _ndrop.ndrop_loadaer_number_codon(
        istart,
        istop,
        k,
        pcols,
        phase,
        voltonumblo,
        voltonumbhi,
        raer_p,
        qqcw_p,
        cs_p,
        vaerosol_p,
        naerosol_p,
    )


@export
def ndrop_loadaer_direct_codon(
    istart: int,
    istop: int,
    k: int,
    pcols: int,
    nspec: int,
    phase: int,
    voltonumblo: float,
    voltonumbhi: float,
    species_raer_ptrs_p: cobj,
    species_qqcw_ptrs_p: cobj,
    specdens_p: cobj,
    spechygro_p: cobj,
    num_raer_p: cobj,
    num_qqcw_p: cobj,
    cs_p: cobj,
    vaerosol_p: cobj,
    hygro_p: cobj,
    naerosol_p: cobj,
):
    _ndrop.ndrop_loadaer_direct_codon(
        istart,
        istop,
        k,
        pcols,
        nspec,
        phase,
        voltonumblo,
        voltonumbhi,
        species_raer_ptrs_p,
        species_qqcw_ptrs_p,
        specdens_p,
        spechygro_p,
        num_raer_p,
        num_qqcw_p,
        cs_p,
        vaerosol_p,
        hygro_p,
        naerosol_p,
    )


@export
def ndrop_ccncalc_zero_codon(
    pcols: int,
    pver: int,
    psat: int,
    ccn_p: cobj,
):
    _ndrop.ndrop_ccncalc_zero_codon(pcols, pver, psat, ccn_p)


@export
def ndrop_ccncalc_level_coeffs_codon(
    ncol: int,
    k: int,
    pcols: int,
    surften_coef: float,
    smcoefcoef: float,
    tair_p: cobj,
    smcoef_p: cobj,
):
    _ndrop.ndrop_ccncalc_level_coeffs_codon(ncol, k, pcols, surften_coef, smcoefcoef, tair_p, smcoef_p)


@export
def ndrop_ccncalc_mode_accum_codon(
    ncol: int,
    k: int,
    pcols: int,
    pver: int,
    psat: int,
    amcubecoef_m: float,
    argfactor_m: float,
    naerosol_p: cobj,
    vaerosol_p: cobj,
    hygro_p: cobj,
    smcoef_p: cobj,
    super_p: cobj,
    amcube_p: cobj,
    sm_p: cobj,
    arg_p: cobj,
    ccn_p: cobj,
):
    _ndrop.ndrop_ccncalc_mode_accum_codon(
        ncol,
        k,
        pcols,
        pver,
        psat,
        amcubecoef_m,
        argfactor_m,
        naerosol_p,
        vaerosol_p,
        hygro_p,
        smcoef_p,
        super_p,
        amcube_p,
        sm_p,
        arg_p,
        ccn_p,
    )


@export
def ndrop_ccncalc_scale_codon(
    ncol: int,
    pcols: int,
    pver: int,
    psat: int,
    ccn_p: cobj,
):
    _ndrop.ndrop_ccncalc_scale_codon(ncol, pcols, pver, psat, ccn_p)


@export
def ndrop_ccncalc_direct_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    psat: int,
    ntot_amode: int,
    ncnst_tot: int,
    pi_value: float,
    surften_coef: float,
    smcoefcoef: float,
    nspec_amode_p: cobj,
    species_raer_ptrs_p: cobj,
    species_qqcw_ptrs_p: cobj,
    specdens_p: cobj,
    spechygro_p: cobj,
    num_raer_ptrs_p: cobj,
    num_qqcw_ptrs_p: cobj,
    voltonumblo_p: cobj,
    voltonumbhi_p: cobj,
    alogsig_p: cobj,
    exp45logsig_p: cobj,
    tair_p: cobj,
    cs_p: cobj,
    super_p: cobj,
    naerosol_p: cobj,
    vaerosol_p: cobj,
    hygro_p: cobj,
    amcube_p: cobj,
    smcoef_p: cobj,
    sm_p: cobj,
    arg_p: cobj,
    ccn_p: cobj,
):
    _ndrop.ndrop_ccncalc_direct_codon(
        ncol,
        pcols,
        pver,
        top_lev,
        psat,
        ntot_amode,
        ncnst_tot,
        pi_value,
        surften_coef,
        smcoefcoef,
        nspec_amode_p,
        species_raer_ptrs_p,
        species_qqcw_ptrs_p,
        specdens_p,
        spechygro_p,
        num_raer_ptrs_p,
        num_qqcw_ptrs_p,
        voltonumblo_p,
        voltonumbhi_p,
        alogsig_p,
        exp45logsig_p,
        tair_p,
        cs_p,
        super_p,
        naerosol_p,
        vaerosol_p,
        hygro_p,
        amcube_p,
        smcoef_p,
        sm_p,
        arg_p,
        ccn_p,
    )


@export
def ndrop_activate_modal_core_codon(
    wbar: float,
    sigw: float,
    wdiab: float,
    wminf: float,
    wmaxf: float,
    tair: float,
    rhoair: float,
    qs: float,
    nmode: int,
    rair: float,
    p0: float,
    t0: float,
    rhoh2o: float,
    latvap: float,
    cpair: float,
    rh2o: float,
    gravit: float,
    pi_value: float,
    aten: float,
    twothird: float,
    sq2: float,
    sqpi: float,
    sixth: float,
    zero_value: float,
    na_p: cobj,
    volume_p: cobj,
    hygro_p: cobj,
    alogsig_p: cobj,
    exp45logsig_p: cobj,
    f1_p: cobj,
    f2_p: cobj,
    fn_p: cobj,
    fm_p: cobj,
    fluxn_p: cobj,
    fluxm_p: cobj,
    flux_fullact_p: cobj,
    zeta_p: cobj,
    eta_p: cobj,
    etafactor2_p: cobj,
    sqrtg_p: cobj,
    amcube_p: cobj,
    smc_p: cobj,
    lnsm_p: cobj,
    sumflxn_p: cobj,
    sumflxm_p: cobj,
    sumfn_p: cobj,
    sumfm_p: cobj,
    fnold_p: cobj,
    fmold_p: cobj,
) -> int:
    return _ndrop.ndrop_activate_modal_core_codon(
        wbar,
        sigw,
        wdiab,
        wminf,
        wmaxf,
        tair,
        rhoair,
        qs,
        nmode,
        rair,
        p0,
        t0,
        rhoh2o,
        latvap,
        cpair,
        rh2o,
        gravit,
        pi_value,
        aten,
        twothird,
        sq2,
        sqpi,
        sixth,
        zero_value,
        na_p,
        volume_p,
        hygro_p,
        alogsig_p,
        exp45logsig_p,
        f1_p,
        f2_p,
        fn_p,
        fm_p,
        fluxn_p,
        fluxm_p,
        flux_fullact_p,
        zeta_p,
        eta_p,
        etafactor2_p,
        sqrtg_p,
        amcube_p,
        smc_p,
        lnsm_p,
        sumflxn_p,
        sumflxm_p,
        sumfn_p,
        sumfm_p,
        fnold_p,
        fmold_p,
    )


@export
def ndrop_maxsat_codon(
    nmode: int,
    zeta_p: cobj,
    eta_p: cobj,
    smc_p: cobj,
    f1_p: cobj,
    f2_p: cobj,
) -> float:
    return _ndrop.ndrop_maxsat_codon(nmode, zeta_p, eta_p, smc_p, f1_p, f2_p)


@export
def ndrop_explmix_codon(
    pver: int,
    top_lev: int,
    surfrate: float,
    flxconv: float,
    dt: float,
    is_unact: int,
    q_p: cobj,
    src_p: cobj,
    ekkp_p: cobj,
    ekkm_p: cobj,
    overlapp_p: cobj,
    overlapm_p: cobj,
    qold_p: cobj,
    qactold_p: cobj,
):
    _ndrop.ndrop_explmix_codon(
        pver,
        top_lev,
        surfrate,
        flxconv,
        dt,
        is_unact,
        q_p,
        src_p,
        ekkp_p,
        ekkm_p,
        overlapp_p,
        overlapm_p,
        qold_p,
        qactold_p,
    )


@export
def ndrop_dropmixnuc_submix_all_codon(
    pver: int,
    top_lev: int,
    ntot_amode: int,
    ncnst_tot: int,
    dtmix: float,
    taumix_internal_pver_inv: float,
    nact_p: cobj,
    mact_p: cobj,
    mam_idx_p: cobj,
    nspec_amode_p: cobj,
    ekkp_p: cobj,
    ekkm_p: cobj,
    overlapp_p: cobj,
    overlapm_p: cobj,
    qcld_p: cobj,
    qncld_p: cobj,
    srcn_p: cobj,
    source_p: cobj,
    raercol_p: cobj,
    raercol_cw_p: cobj,
    nsav_p: cobj,
    nnew_p: cobj,
):
    _ndrop.ndrop_dropmixnuc_submix_all_codon(
        pver,
        top_lev,
        ntot_amode,
        ncnst_tot,
        dtmix,
        taumix_internal_pver_inv,
        nact_p,
        mact_p,
        mam_idx_p,
        nspec_amode_p,
        ekkp_p,
        ekkm_p,
        overlapp_p,
        overlapm_p,
        qcld_p,
        qncld_p,
        srcn_p,
        source_p,
        raercol_p,
        raercol_cw_p,
        nsav_p,
        nnew_p,
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
def ghg_data_trcmix_scale_codon(gas_id: int, dlat: float) -> float:
    if gas_id == 1:
        if dlat <= 45.0:
            return 0.2353
        return 0.2353 + 0.0225489 * (dlat - 45.0)

    if gas_id == 2:
        if dlat <= 45.0:
            return 0.3478 + 0.00116 * dlat
        return 0.4000 + 0.013333 * (dlat - 45.0)

    if gas_id == 3:
        if dlat <= 45.0:
            return 0.7273 + 0.00606 * dlat
        return 1.00 + 0.013333 * (dlat - 45.0)

    if dlat <= 45.0:
        return 0.4000 + 0.00222 * dlat
    return 0.50 + 0.024444 * (dlat - 45.0)


@export
def ghg_data_trcmix_codon(
    gas_id: int,
    ncol: int,
    pcols: int,
    pver: int,
    trop_mmr: float,
    constant_mmr: float,
    clat_p: cobj,
    pmid_p: cobj,
    q_p: cobj,
):
    clat = Ptr[float](clat_p)
    pmid = Ptr[float](pmid_p)
    q = Ptr[float](q_p)

    if gas_id >= 5:
        for k in range(1, pver + 1):
            for i in range(1, pcols + 1):
                q[_field2_idx(i, k, pcols)] = constant_mmr
        return

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            coslat = cos(clat[_idx(i)])
            dlat = abs(57.2958 * clat[_idx(i)])
            scale = ghg_data_trcmix_scale_codon(gas_id, dlat)
            ptrop = 250.0e2 - 150.0e2 * coslat ** 2.0
            idx = _field2_idx(i, k, pcols)
            if pmid[idx] >= ptrop:
                q[idx] = trop_mmr
            else:
                pratio = pmid[idx] / ptrop
                q[idx] = trop_mmr * (pratio) ** scale


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
def rad_constituents_touch_codon(stage: int) -> int:
    return stage


@export
def rad_cnst_readnl_codon(stage: int) -> int:
    return rad_constituents_touch_codon(stage)


@export
def rad_cnst_init_codon(stage: int) -> int:
    return rad_constituents_touch_codon(stage)


@export
def rad_cnst_get_gas_codon(stage: int) -> int:
    return rad_constituents_touch_codon(stage)


@export
def rad_cnst_get_info_codon(stage: int) -> int:
    return rad_constituents_touch_codon(stage)


@export
def rad_cnst_get_info_by_mode_codon(stage: int) -> int:
    return rad_constituents_touch_codon(stage)


@export
def rad_cnst_get_info_by_mode_spec_codon(stage: int) -> int:
    return rad_constituents_touch_codon(stage)


@export
def rad_cnst_get_info_by_spectype_codon(stage: int) -> int:
    return rad_constituents_touch_codon(stage)


@export
def rad_cnst_out_codon(stage: int) -> int:
    return rad_constituents_touch_codon(stage)


@export
def init_mode_comps_codon(stage: int) -> int:
    return rad_constituents_touch_codon(stage)


@export
def get_cam_idx_codon(stage: int) -> int:
    return rad_constituents_touch_codon(stage)


@export
def list_init1_codon(stage: int) -> int:
    return rad_constituents_touch_codon(stage)


@export
def list_init2_codon(stage: int) -> int:
    return rad_constituents_touch_codon(stage)


@export
def rad_gas_diag_init_codon(stage: int) -> int:
    return rad_constituents_touch_codon(stage)


@export
def parse_rad_specifier_codon(stage: int) -> int:
    return rad_constituents_touch_codon(stage)


@export
def rad_cnst_get_mam_mmr_by_idx_codon(stage: int) -> int:
    return rad_constituents_touch_codon(stage)


@export
def rad_cnst_get_mode_num_codon(stage: int) -> int:
    return rad_constituents_touch_codon(stage)


@export
def rad_cnst_get_mam_props_by_idx_codon(stage: int) -> int:
    return rad_constituents_touch_codon(stage)


@export
def rad_cnst_get_mode_props_codon(stage: int) -> int:
    return rad_constituents_touch_codon(stage)


@export
def print_modes_codon(stage: int) -> int:
    return rad_constituents_touch_codon(stage)


@export
def print_lists_codon(stage: int) -> int:
    return rad_constituents_touch_codon(stage)


@export
def rad_cnst_get_call_list_codon(n: int, active_p: cobj, call_list_p: cobj):
    active = Ptr[int](active_p)
    call_list = Ptr[int](call_list_p)

    for i in range(n):
        if active[i] != 0:
            call_list[i] = 1
        else:
            call_list[i] = 0


@inline
def _ascii_eq(text: Ptr[int], n: int, a0: int, a1: int, a2: int, a3: int, a4: int, a5: int, a6: int,
              a7: int, a8: int, a9: int, a10: int, a11: int, a12: int, a13: int) -> bool:
    if n >= 1 and text[0] != a0:
        return False
    if n >= 2 and text[1] != a1:
        return False
    if n >= 3 and text[2] != a2:
        return False
    if n >= 4 and text[3] != a3:
        return False
    if n >= 5 and text[4] != a4:
        return False
    if n >= 6 and text[5] != a5:
        return False
    if n >= 7 and text[6] != a6:
        return False
    if n >= 8 and text[7] != a7:
        return False
    if n >= 9 and text[8] != a8:
        return False
    if n >= 10 and text[9] != a9:
        return False
    if n >= 11 and text[10] != a10:
        return False
    if n >= 12 and text[11] != a11:
        return False
    if n >= 13 and text[12] != a12:
        return False
    if n >= 14 and text[13] != a13:
        return False
    return True


@export
def rad_cnst_check_specie_type_codon(n: int, text_p: cobj) -> int:
    text = Ptr[int](text_p)
    if n == 4:
        if _ascii_eq(text, n, 100, 117, 115, 116, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0):
            return 1
    if n == 7:
        if _ascii_eq(text, n, 110, 105, 116, 114, 97, 116, 101, 0, 0, 0, 0, 0, 0, 0):
            return 1
        if _ascii_eq(text, n, 115, 117, 108, 102, 97, 116, 101, 0, 0, 0, 0, 0, 0, 0):
            return 1
        if _ascii_eq(text, n, 115, 101, 97, 115, 97, 108, 116, 0, 0, 0, 0, 0, 0, 0):
            return 1
        if _ascii_eq(text, n, 98, 108, 97, 99, 107, 45, 99, 0, 0, 0, 0, 0, 0, 0):
            return 1
    if n == 8:
        if _ascii_eq(text, n, 97, 109, 109, 111, 110, 105, 117, 109, 0, 0, 0, 0, 0, 0):
            return 1
    if n == 9:
        if _ascii_eq(text, n, 112, 45, 111, 114, 103, 97, 110, 105, 99, 0, 0, 0, 0, 0):
            return 1
        if _ascii_eq(text, n, 115, 45, 111, 114, 103, 97, 110, 105, 99, 0, 0, 0, 0, 0):
            return 1
    return 0


@export
def check_specie_type_codon(n: int, text_p: cobj) -> int:
    return rad_cnst_check_specie_type_codon(n, text_p)


@export
def rad_cnst_check_mode_type_codon(n: int, text_p: cobj) -> int:
    text = Ptr[int](text_p)
    if n == 5:
        if _ascii_eq(text, n, 97, 99, 99, 117, 109, 0, 0, 0, 0, 0, 0, 0, 0, 0):
            return 1
    if n == 6:
        if _ascii_eq(text, n, 99, 111, 97, 114, 115, 101, 0, 0, 0, 0, 0, 0, 0, 0):
            return 1
        if _ascii_eq(text, n, 97, 105, 116, 107, 101, 110, 0, 0, 0, 0, 0, 0, 0, 0):
            return 1
    if n == 9:
        if _ascii_eq(text, n, 102, 105, 110, 101, 95, 100, 117, 115, 116, 0, 0, 0, 0, 0):
            return 1
    if n == 11:
        if _ascii_eq(text, n, 99, 111, 97, 114, 115, 101, 95, 100, 117, 115, 116, 0, 0, 0):
            return 1
    if n == 12:
        if _ascii_eq(text, n, 102, 105, 110, 101, 95, 115, 101, 97, 115, 97, 108, 116, 0, 0):
            return 1
    if n == 14:
        if _ascii_eq(text, n, 112, 114, 105, 109, 97, 114, 121, 95, 99, 97, 114, 98, 111, 110):
            return 1
        if _ascii_eq(text, n, 99, 111, 97, 114, 115, 101, 95, 115, 101, 97, 115, 97, 108, 116):
            return 1
    return 0


@export
def check_mode_type_codon(n: int, text_p: cobj) -> int:
    return rad_cnst_check_mode_type_codon(n, text_p)


@inline
def _parse_mode_defs_idx(row1: int, col1: int, ncols: int) -> int:
    return (row1 - 1) + (col1 - 1) * ncols


@inline
def _parse_mode_defs_field_idx(mode1: int, spec1: int, max_spec: int) -> int:
    return (mode1 - 1) * max_spec + (spec1 - 1)


@inline
def _parse_mode_defs_token_start(field: str, max_len: int, target: Ptr[int], offset: int):
    for i in range(max_len):
        target[offset + i] = 32
    n = min(len(field), max_len)
    for i in range(n):
        target[offset + i] = ord(field[i])


@inline
def _parse_mode_defs_clean_input(src: Ptr[int], row1: int, cs1: int, n_mode_str: int,
                                 clean: Ptr[int], clean_len: Ptr[int]) -> str:
    started = False
    out = ""
    for j in range(cs1):
        ch = src[_parse_mode_defs_idx(row1, j + 1, n_mode_str)]
        if ch == 0:
            ch = 32
        if not started:
            if ch == 32:
                clean[_parse_mode_defs_idx(row1, j + 1, n_mode_str)] = 32
                continue
            started = True
        if ch != 32:
            out += chr(ch)

    for j in range(cs1):
        clean[_parse_mode_defs_idx(row1, j + 1, n_mode_str)] = 32

    n = min(len(out), cs1)
    for j in range(n):
        clean[_parse_mode_defs_idx(row1, j + 1, n_mode_str)] = ord(out[j])
    clean_len[row1 - 1] = n
    return out


@inline
def _parse_mode_defs_write_error(msg: str, err: Ptr[int]):
    for i in range(256):
        err[i] = 32
    n = min(len(msg), 256)
    for i in range(n):
        err[i] = ord(msg[i])


@inline
def _parse_mode_defs_valid_source(field: str) -> bool:
    return field == "A" or field == "N" or field == "Z"


@inline
def _parse_mode_defs_valid_spec(field: str) -> bool:
    if len(field) == 4:
        return field == "dust"
    if len(field) == 7:
        return field == "nitrate" or field == "sulfate" or field == "seasalt" or field == "black-c"
    if len(field) == 8:
        return field == "ammonium"
    if len(field) == 9:
        return field == "p-organic" or field == "s-organic"
    return False


@inline
def _parse_mode_defs_valid_mode(field: str) -> bool:
    if len(field) == 5:
        return field == "accum"
    if len(field) == 6:
        return field == "coarse" or field == "aitken"
    if len(field) == 9:
        return field == "fine_dust"
    if len(field) == 11:
        return field == "coarse_dust"
    if len(field) == 12:
        return field == "fine_seasalt"
    if len(field) == 14:
        return field == "primary_carbon" or field == "coarse_seasalt"
    return False


@inline
def _parse_mode_defs_has_continue(fields: List[str]) -> bool:
    return len(fields) > 0 and fields[len(fields) - 1] == "+"


@export
def parse_mode_defs_plan_codon(
    n_mode_str: int,
    cs1: int,
    max_spec: int,
    nl_ascii_p: cobj,
    nl_clean_ascii_p: cobj,
    clean_len_p: cobj,
    nmodes_p: cobj,
    nstr_p: cobj,
    mode_nspec_p: cobj,
    mode_names_p: cobj,
    mode_types_p: cobj,
    source_num_a_p: cobj,
    camname_num_a_p: cobj,
    source_num_c_p: cobj,
    camname_num_c_p: cobj,
    source_mmr_a_p: cobj,
    camname_mmr_a_p: cobj,
    source_mmr_c_p: cobj,
    camname_mmr_c_p: cobj,
    spec_type_p: cobj,
    spec_props_p: cobj,
    status_p: cobj,
    err_string_p: cobj,
):
    nl_ascii = Ptr[int](nl_ascii_p)
    nl_clean_ascii = Ptr[int](nl_clean_ascii_p)
    clean_len = Ptr[int](clean_len_p)
    nmodes_out = Ptr[int](nmodes_p)
    nstr_out = Ptr[int](nstr_p)
    mode_nspec = Ptr[int](mode_nspec_p)
    mode_names = Ptr[int](mode_names_p)
    mode_types = Ptr[int](mode_types_p)
    source_num_a = Ptr[int](source_num_a_p)
    camname_num_a = Ptr[int](camname_num_a_p)
    source_num_c = Ptr[int](source_num_c_p)
    camname_num_c = Ptr[int](camname_num_c_p)
    source_mmr_a = Ptr[int](source_mmr_a_p)
    camname_mmr_a = Ptr[int](camname_mmr_a_p)
    source_mmr_c = Ptr[int](source_mmr_c_p)
    camname_mmr_c = Ptr[int](camname_mmr_c_p)
    spec_type = Ptr[int](spec_type_p)
    spec_props = Ptr[int](spec_props_p)
    status = Ptr[int](status_p)
    err_string = Ptr[int](err_string_p)

    status[0] = 0
    nmodes_out[0] = 0
    nstr_out[0] = 0
    _parse_mode_defs_write_error("", err_string)

    clean_strings = List[str]()
    for m in range(1, n_mode_str + 1):
        text = _parse_mode_defs_clean_input(nl_ascii, m, cs1, n_mode_str, nl_clean_ascii, clean_len)
        if len(text) == 0:
            break
        clean_strings.append(text)
        nstr_out[0] += 1
        if len(text) >= 2 and text[len(text) - 2:] == ":=":
            nmodes_out[0] += 1

    for mode in range(1, n_mode_str + 1):
        mode_nspec[mode - 1] = 0
        _parse_mode_defs_token_start("", 32, mode_names, (mode - 1) * 32)
        _parse_mode_defs_token_start("", 32, mode_types, (mode - 1) * 32)
        _parse_mode_defs_token_start("", 1, source_num_a, mode - 1)
        _parse_mode_defs_token_start("", 32, camname_num_a, (mode - 1) * 32)
        _parse_mode_defs_token_start("", 1, source_num_c, mode - 1)
        _parse_mode_defs_token_start("", 32, camname_num_c, (mode - 1) * 32)
        for spec in range(1, max_spec + 1):
            field_idx = _parse_mode_defs_field_idx(mode, spec, max_spec)
            _parse_mode_defs_token_start("", 1, source_mmr_a, field_idx)
            _parse_mode_defs_token_start("", 32, camname_mmr_a, field_idx * 32)
            _parse_mode_defs_token_start("", 1, source_mmr_c, field_idx)
            _parse_mode_defs_token_start("", 32, camname_mmr_c, field_idx * 32)
            _parse_mode_defs_token_start("", 32, spec_type, field_idx * 32)
            _parse_mode_defs_token_start("", cs1, spec_props, field_idx * cs1)

    if nmodes_out[0] == 0:
        return

    mcur = 0
    for mode in range(1, nmodes_out[0] + 1):
        if mcur >= len(clean_strings):
            status[0] = -1
            _parse_mode_defs_write_error("mode definition missing", err_string)
            return

        first = clean_strings[mcur]
        if len(first) < 2 or first[len(first) - 2:] != ":=":
            status[0] = -2
            _parse_mode_defs_write_error("= not found", err_string)
            return

        head = first[: len(first) - 2]
        parts = head.split(":")
        if len(parts) != 2 or len(parts[0]) == 0:
            status[0] = -3
            _parse_mode_defs_write_error("mode name/type not found", err_string)
            return
        if not _parse_mode_defs_valid_mode(parts[1]):
            status[0] = -4
            _parse_mode_defs_write_error("mode type not valid", err_string)
            return

        _parse_mode_defs_token_start(parts[0], 32, mode_names, (mode - 1) * 32)
        _parse_mode_defs_token_start(parts[1], 32, mode_types, (mode - 1) * 32)

        mcur += 1
        nspec = 0
        scan_idx = mcur
        while scan_idx < len(clean_strings):
            text = clean_strings[scan_idx]
            fields = text.split(":")
            if len(fields) < 6:
                status[0] = -5
                _parse_mode_defs_write_error("component field missing", err_string)
                return
            if _parse_mode_defs_has_continue(fields):
                nspec += 1
                scan_idx += 1
            else:
                break
        if nspec == 0:
            status[0] = -6
            _parse_mode_defs_write_error("mode must have at least one specie", err_string)
            return
        if nspec > max_spec:
            status[0] = -7
            _parse_mode_defs_write_error("too many species in mode", err_string)
            return
        mode_nspec[mode - 1] = nspec

        num_mr_found = False
        ispec = 0
        while mcur < len(clean_strings):
            text = clean_strings[mcur]
            fields = text.split(":")
            if len(fields) < 6:
                status[0] = -8
                _parse_mode_defs_write_error("component field missing", err_string)
                return
            src_a = fields[0]
            name_a = fields[1]
            src_c = fields[2]
            name_c = fields[3]
            comp_type = fields[4]
            if not _parse_mode_defs_valid_source(src_a) or not _parse_mode_defs_valid_source(src_c):
                status[0] = -9
                _parse_mode_defs_write_error("source must be A, N or Z", err_string)
                return

            if comp_type == "num_mr":
                if num_mr_found:
                    status[0] = -10
                    _parse_mode_defs_write_error("more than 1 number component", err_string)
                    return
                num_mr_found = True
                _parse_mode_defs_token_start(src_a, 1, source_num_a, mode - 1)
                _parse_mode_defs_token_start(name_a, 32, camname_num_a, (mode - 1) * 32)
                _parse_mode_defs_token_start(src_c, 1, source_num_c, mode - 1)
                _parse_mode_defs_token_start(name_c, 32, camname_num_c, (mode - 1) * 32)
            else:
                if not _parse_mode_defs_valid_spec(comp_type):
                    status[0] = -11
                    _parse_mode_defs_write_error("specie type not valid", err_string)
                    return
                if len(fields) < 6 or len(fields[5]) < 3 or fields[5][len(fields[5]) - 3:] != ".nc":
                    status[0] = -12
                    _parse_mode_defs_write_error("filename not valid", err_string)
                    return
                ispec += 1
                if ispec > nspec:
                    status[0] = -13
                    _parse_mode_defs_write_error("component parsing got wrong number of species", err_string)
                    return
                field_idx = _parse_mode_defs_field_idx(mode, ispec, max_spec)
                _parse_mode_defs_token_start(src_a, 1, source_mmr_a, field_idx)
                _parse_mode_defs_token_start(name_a, 32, camname_mmr_a, field_idx * 32)
                _parse_mode_defs_token_start(src_c, 1, source_mmr_c, field_idx)
                _parse_mode_defs_token_start(name_c, 32, camname_mmr_c, field_idx * 32)
                _parse_mode_defs_token_start(comp_type, 32, spec_type, field_idx * 32)
                _parse_mode_defs_token_start(fields[5], cs1, spec_props, field_idx * cs1)

            mcur += 1
            if _parse_mode_defs_has_continue(fields):
                continue
            break

        if not num_mr_found:
            status[0] = -14
            _parse_mode_defs_write_error("number component not found", err_string)
            return
        if ispec != nspec:
            status[0] = -15
            _parse_mode_defs_write_error("component parsing got wrong number of species", err_string)
            return


@export
def rad_cnst_mam_mmr_idx_codon(
    mode_idx: int,
    spec_idx: int,
    nmodes: int,
    nspec: int,
    idx_mmr_a: int,
    idx_p: cobj,
    status_p: cobj,
):
    idx = Ptr[int](idx_p)
    status = Ptr[int](status_p)

    idx[0] = 0
    status[0] = 0

    if mode_idx < 1 or mode_idx > nmodes:
        status[0] = -1
        return

    if nspec >= 0:
        if spec_idx < 1 or spec_idx > nspec:
            status[0] = -2
            return

        idx[0] = idx_mmr_a


@export
def rad_cnst_get_mam_mmr_idx_codon(
    mode_idx: int,
    spec_idx: int,
    nmodes: int,
    nspec: int,
    idx_mmr_a: int,
    idx_p: cobj,
    status_p: cobj,
):
    rad_cnst_mam_mmr_idx_codon(
        mode_idx,
        spec_idx,
        nmodes,
        nspec,
        idx_mmr_a,
        idx_p,
        status_p,
    )


@export
def rad_cnst_mode_num_idx_codon(
    mode_idx: int,
    nmodes: int,
    source_ascii: int,
    idx_num_a: int,
    idx_p: cobj,
    status_p: cobj,
):
    idx = Ptr[int](idx_p)
    status = Ptr[int](status_p)

    idx[0] = 0
    status[0] = 0

    if mode_idx < 1 or mode_idx > nmodes:
        status[0] = -1
        return

    if source_ascii >= 0:
        if source_ascii != 65:
            status[0] = -2
            return

        idx[0] = idx_num_a


@export
def rad_cnst_get_mode_num_idx_codon(
    mode_idx: int,
    nmodes: int,
    source_ascii: int,
    idx_num_a: int,
    idx_p: cobj,
    status_p: cobj,
):
    rad_cnst_mode_num_idx_codon(
        mode_idx,
        nmodes,
        source_ascii,
        idx_num_a,
        idx_p,
        status_p,
    )


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
def phys_control_bool_flag_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def phys_control_index_positive_codon(index_value: int) -> int:
    if index_value > 0:
        return 1
    return 0


@export
def phys_control_int_value_codon(value: int) -> int:
    return value


@inline
def _trimmed_ascii_len(text: Ptr[int], n: int) -> int:
    last = n
    while last > 0 and text[last - 1] == 32:
        last -= 1
    return last


@inline
def _trimmed_ascii_eq(name_len: int, name_p: cobj, value_len: int, value_p: cobj) -> int:
    name = Ptr[int](name_p)
    value = Ptr[int](value_p)

    n_name = _trimmed_ascii_len(name, name_len)
    n_value = _trimmed_ascii_len(value, value_len)
    if n_name != n_value:
        return 0

    for i in range(n_name):
        if name[i] != value[i]:
            return 0
    return 1


@export
def cam_physpkg_is_codon(name_len: int, name_p: cobj, pkg_len: int, pkg_p: cobj) -> int:
    return _trimmed_ascii_eq(name_len, name_p, pkg_len, pkg_p)


@export
def cam_chempkg_is_codon(name_len: int, name_p: cobj, pkg_len: int, pkg_p: cobj) -> int:
    return _trimmed_ascii_eq(name_len, name_p, pkg_len, pkg_p)


@export
def waccmx_is_codon(name_len: int, name_p: cobj, opt_len: int, opt_p: cobj) -> int:
    return _trimmed_ascii_eq(name_len, name_p, opt_len, opt_p)


@inline
def _ascii_first(name_len: int, name_ascii_p: cobj) -> int:
    if name_len <= 0:
        return 0
    name_ascii = Ptr[int](name_ascii_p)
    return name_ascii[0]


@inline
def _vdiff_field_idx(name_len: int, name_ascii_p: cobj, has_qindex: int, qindex: int) -> int:
    if name_len != 1:
        return 0
    first = _ascii_first(name_len, name_ascii_p)
    if first == 117 or first == 85:
        return 1
    if first == 118 or first == 86:
        return 2
    if first == 115 or first == 83:
        return 3
    if first == 113 or first == 81:
        if has_qindex != 0:
            return 3 + qindex
        return 4
    return 0


@export
def vdiff_select_codon(name_len: int, name_ascii_p: cobj, has_qindex: int, qindex: int) -> int:
    return _vdiff_field_idx(name_len, name_ascii_p, has_qindex, qindex)


@export
def diffuse_codon(name_len: int, name_ascii_p: cobj, has_qindex: int, qindex: int) -> int:
    return _vdiff_field_idx(name_len, name_ascii_p, has_qindex, qindex)


@export
def new_fieldlist_vdiff_codon(ncnst: int) -> int:
    return 3 + ncnst


@export
def my_any_codon(n: int, values_p: cobj) -> int:
    values = Ptr[int](values_p)
    for i in range(n):
        if values[i] != 0:
            return 1
    return 0


@export
def constituents_rgas_codon(r_universal: float, mwc: float) -> float:
    return r_universal * mwc


@export
def constituents_cv_codon(cpc: float, rgas: float) -> float:
    return cpc - rgas


@export
def cnst_read_iv_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def cnst_cam_outfld_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def cnst_add_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def cnst_chk_dim_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def cnst_get_ind_codon(name_len: int, name_ascii_p: cobj, cnst_name_len: int, cnst_names_ascii_p: cobj,
                       pcnst: int) -> int:
    name_ascii = Ptr[int](name_ascii_p)
    cnst_names_ascii = Ptr[int](cnst_names_ascii_p)

    m = 1
    while m <= pcnst:
        name_offset = (m - 1) * cnst_name_len
        if _phys_prop_trimmed_eq(name_ascii, name_len, cnst_names_ascii + name_offset, cnst_name_len):
            return m
        m += 1
    return -1


@inline
def _copy_ascii_fixed(n: int, src: Ptr[int], dst: Ptr[int]):
    i = 0
    while i < n:
        dst[i] = src[i]
        i += 1


@export
def cnst_get_type_byind_codon(ind: int, pcnst: int, cnst_type_ascii_p: cobj, out_ascii_p: cobj,
                              status_p: cobj):
    cnst_type_ascii = Ptr[int](cnst_type_ascii_p)
    out_ascii = Ptr[int](out_ascii_p)
    status = Ptr[int](status_p)

    if ind <= pcnst:
        _copy_ascii_fixed(3, cnst_type_ascii + (ind - 1) * 3, out_ascii)
        status[0] = 0
    else:
        status[0] = 1


@export
def cnst_get_molec_byind_codon(ind: int, pcnst: int, cnst_molec_ascii_p: cobj, out_ascii_p: cobj,
                               status_p: cobj):
    cnst_molec_ascii = Ptr[int](cnst_molec_ascii_p)
    out_ascii = Ptr[int](out_ascii_p)
    status = Ptr[int](status_p)

    if ind <= pcnst:
        _copy_ascii_fixed(5, cnst_molec_ascii + (ind - 1) * 5, out_ascii)
        status[0] = 0
    else:
        status[0] = 1


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
def aer_rad_props_init_plan_codon(
    history_amwg: int,
    history_aero_optics: int,
    prog_modal_aero: int,
    clim_modal_aero_top_lev: int,
    top_lev_p: cobj,
    add_amwg_default_p: cobj,
    add_aero_optics_default_p: cobj,
):
    top_lev = Ptr[int](top_lev_p)
    add_amwg_default = Ptr[int](add_amwg_default_p)
    add_aero_optics_default = Ptr[int](add_aero_optics_default_p)

    if prog_modal_aero != 0:
        top_lev[0] = clim_modal_aero_top_lev
    add_amwg_default[0] = 1 if history_amwg != 0 else 0
    add_aero_optics_default[0] = 1 if history_aero_optics != 0 else 0


@export
def aer_vis_diag_prepare_codon(
    ncol: int,
    pcols: int,
    tau_nlev: int,
    nnite: int,
    fillvalue: float,
    idxnite_p: cobj,
    tau_p: cobj,
    troplev_p: cobj,
    tmp_p: cobj,
    tmp2_p: cobj,
):
    idxnite = Ptr[int](idxnite_p)
    tau = Ptr[float](tau_p)
    troplev = Ptr[int](troplev_p)
    tmp = Ptr[float](tmp_p)
    tmp2 = Ptr[float](tmp2_p)

    for i in range(1, pcols + 1):
        tmp[i - 1] = 0.0
        tmp2[i - 1] = 0.0

    for i in range(1, ncol + 1):
        total = 0.0
        for k in range(1, tau_nlev + 1):
            total = total + tau[_field2_idx(i, k, pcols)]
        tmp[i - 1] = total

    for i in range(1, nnite + 1):
        idx = idxnite[i - 1]
        if idx >= 1 and idx <= pcols:
            tmp[idx - 1] = fillvalue

    for i in range(1, ncol + 1):
        trop = troplev[i - 1]
        if trop < 1:
            trop = 1
        if trop > tau_nlev:
            trop = tau_nlev
        strat_total = 0.0
        for k in range(1, trop + 1):
            strat_total = strat_total + tau[_field2_idx(i, k, pcols)]
        tmp2[i - 1] = strat_total


@export
def aer_vis_diag_out_codon(
    ncol: int,
    pcols: int,
    tau_nlev: int,
    nnite: int,
    fillvalue: float,
    idxnite_p: cobj,
    tau_p: cobj,
    troplev_p: cobj,
    tmp_p: cobj,
    tmp2_p: cobj,
):
    aer_vis_diag_prepare_codon(
        ncol,
        pcols,
        tau_nlev,
        nnite,
        fillvalue,
        idxnite_p,
        tau_p,
        troplev_p,
        tmp_p,
        tmp2_p,
    )


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
def micro_mg_utils_init_codon(
    kind: int,
    expected_kind: int,
    rh2o: float,
    cpair: float,
    tmelt_in: float,
    latvap: float,
    latice: float,
    dcs: float,
    pi: float,
    dsph: float,
    bs: float,
    br: float,
    rhow: float,
    rhoi: float,
    rhosn: float,
    min_mean_mass_liq: float,
    min_mean_mass_ice: float,
    no_limiter_bits: int,
    lam_bnd_rain1: float,
    lam_bnd_rain2: float,
    lam_bnd_snow1: float,
    lam_bnd_snow2: float,
    rv_p: cobj,
    cpp_p: cobj,
    tmelt_p: cobj,
    xxlv_p: cobj,
    xlf_p: cobj,
    xxls_p: cobj,
    gamma_bs_plus3_p: cobj,
    gamma_half_br_plus5_p: cobj,
    gamma_half_bs_plus5_p: cobj,
    liq_rho_p: cobj,
    liq_eff_dim_p: cobj,
    liq_shape_coef_p: cobj,
    liq_lambda_bounds_p: cobj,
    liq_min_mean_mass_p: cobj,
    ice_rho_p: cobj,
    ice_eff_dim_p: cobj,
    ice_shape_coef_p: cobj,
    ice_lambda_bounds_p: cobj,
    ice_min_mean_mass_p: cobj,
    rain_rho_p: cobj,
    rain_eff_dim_p: cobj,
    rain_shape_coef_p: cobj,
    rain_lambda_bounds_p: cobj,
    rain_min_mean_mass_p: cobj,
    snow_rho_p: cobj,
    snow_eff_dim_p: cobj,
    snow_shape_coef_p: cobj,
    snow_lambda_bounds_p: cobj,
    snow_min_mean_mass_p: cobj,
    status_p: cobj,
):
    rv = Ptr[float](rv_p)
    cpp = Ptr[float](cpp_p)
    tmelt = Ptr[float](tmelt_p)
    xxlv = Ptr[float](xxlv_p)
    xlf = Ptr[float](xlf_p)
    xxls = Ptr[float](xxls_p)
    gamma_bs_plus3 = Ptr[float](gamma_bs_plus3_p)
    gamma_half_br_plus5 = Ptr[float](gamma_half_br_plus5_p)
    gamma_half_bs_plus5 = Ptr[float](gamma_half_bs_plus5_p)
    liq_rho = Ptr[float](liq_rho_p)
    liq_eff_dim = Ptr[float](liq_eff_dim_p)
    liq_shape_coef = Ptr[float](liq_shape_coef_p)
    liq_lambda_bounds = Ptr[float](liq_lambda_bounds_p)
    liq_min_mean_mass = Ptr[float](liq_min_mean_mass_p)
    ice_rho = Ptr[float](ice_rho_p)
    ice_eff_dim = Ptr[float](ice_eff_dim_p)
    ice_shape_coef = Ptr[float](ice_shape_coef_p)
    ice_lambda_bounds = Ptr[float](ice_lambda_bounds_p)
    ice_min_mean_mass = Ptr[float](ice_min_mean_mass_p)
    rain_rho = Ptr[float](rain_rho_p)
    rain_eff_dim = Ptr[float](rain_eff_dim_p)
    rain_shape_coef = Ptr[float](rain_shape_coef_p)
    rain_lambda_bounds = Ptr[float](rain_lambda_bounds_p)
    rain_min_mean_mass = Ptr[float](rain_min_mean_mass_p)
    snow_rho = Ptr[float](snow_rho_p)
    snow_eff_dim = Ptr[float](snow_eff_dim_p)
    snow_shape_coef = Ptr[float](snow_shape_coef_p)
    snow_lambda_bounds = Ptr[float](snow_lambda_bounds_p)
    snow_min_mean_mass = Ptr[float](snow_min_mean_mass_p)
    status = Ptr[int](status_p)

    if kind != expected_kind:
        status[0] = 1
        return

    status[0] = 0
    rv[0] = rh2o
    cpp[0] = cpair
    tmelt[0] = tmelt_in
    xxlv[0] = latvap
    xlf[0] = latice
    xxls[0] = xxlv[0] + xlf[0]

    gamma_bs_plus3[0] = micro_mg_utils_gamma_native_cb(3.0 + bs)
    gamma_half_br_plus5[0] = micro_mg_utils_gamma_native_cb(5.0 / 2.0 + br / 2.0)
    gamma_half_bs_plus5[0] = micro_mg_utils_gamma_native_cb(5.0 / 2.0 + bs / 2.0)

    liq_rho[0] = rhow
    liq_eff_dim[0] = dsph
    liq_lambda_bounds[0] = float(no_limiter_bits)
    liq_lambda_bounds[1] = float(no_limiter_bits)
    Ptr[int](liq_lambda_bounds_p)[0] = no_limiter_bits
    Ptr[int](liq_lambda_bounds_p)[1] = no_limiter_bits
    liq_min_mean_mass[0] = min_mean_mass_liq
    liq_shape_coef[0] = rhow * pi * micro_mg_utils_gamma_native_cb(dsph + 1.0) / 6.0

    ice_rho[0] = rhoi
    ice_eff_dim[0] = dsph
    ice_lambda_bounds[0] = 1.0 / (2.0 * dcs)
    ice_lambda_bounds[1] = 1.0 / 10.0e-6
    ice_min_mean_mass[0] = min_mean_mass_ice
    ice_shape_coef[0] = rhoi * pi * micro_mg_utils_gamma_native_cb(dsph + 1.0) / 6.0

    rain_rho[0] = rhow
    rain_eff_dim[0] = dsph
    rain_lambda_bounds[0] = lam_bnd_rain1
    rain_lambda_bounds[1] = lam_bnd_rain2
    Ptr[int](rain_min_mean_mass_p)[0] = no_limiter_bits
    rain_shape_coef[0] = rhow * pi * micro_mg_utils_gamma_native_cb(dsph + 1.0) / 6.0

    snow_rho[0] = rhosn
    snow_eff_dim[0] = dsph
    snow_lambda_bounds[0] = lam_bnd_snow1
    snow_lambda_bounds[1] = lam_bnd_snow2
    Ptr[int](snow_min_mean_mass_p)[0] = no_limiter_bits
    snow_shape_coef[0] = rhosn * pi * micro_mg_utils_gamma_native_cb(dsph + 1.0) / 6.0


@export
def no_limiter_codon() -> int:
    return 0x7FF1111111111111


@export
def limiter_is_on_codon(bits: int, off_bits: int) -> int:
    if bits != off_bits:
        return 1
    return 0


@export
def newmghydrometeorprops_codon(
    rho_value: float,
    eff_dim_value: float,
    has_lambda_bounds: int,
    lambda_bounds1: float,
    lambda_bounds2: float,
    has_min_mean_mass: int,
    min_mean_mass_value: float,
    pi_value: float,
    no_limiter_bits: int,
    rho_p: cobj,
    eff_dim_p: cobj,
    shape_coef_p: cobj,
    lambda_bounds_p: cobj,
    min_mean_mass_p: cobj,
):
    rho = Ptr[float](rho_p)
    eff_dim = Ptr[float](eff_dim_p)
    shape_coef = Ptr[float](shape_coef_p)
    lambda_bounds = Ptr[float](lambda_bounds_p)
    min_mean_mass = Ptr[float](min_mean_mass_p)

    rho[0] = rho_value
    eff_dim[0] = eff_dim_value
    if has_lambda_bounds != 0:
        lambda_bounds[0] = lambda_bounds1
        lambda_bounds[1] = lambda_bounds2
    else:
        Ptr[int](lambda_bounds_p)[0] = no_limiter_bits
        Ptr[int](lambda_bounds_p)[1] = no_limiter_bits

    if has_min_mean_mass != 0:
        min_mean_mass[0] = min_mean_mass_value
    else:
        Ptr[int](min_mean_mass_p)[0] = no_limiter_bits

    shape_coef[0] = micro_mg_utils_shape_coef_native_cb(rho_value, eff_dim_value)


@export
def rising_factorial_codon(x: float, n: float) -> float:
    return micro_mg_utils_rising_factorial_native_cb(x, n)


@export
def size_dist_param_basic_codon(
    qsmall_value: float,
    no_limiter_bits: int,
    prop_eff_dim: float,
    prop_shape_coef: float,
    lambda_bounds1: float,
    lambda_bounds2: float,
    min_mean_mass: float,
    min_mean_mass_bits: int,
    qic: float,
    nic_in: float,
    want_n0: int,
    nic_out_c: cobj,
    lam_c: cobj,
    n0_c: cobj,
):
    nic_out = Ptr[float](nic_out_c)
    lam_out = Ptr[float](lam_c)
    n0_out = Ptr[float](n0_c)

    nic = nic_in
    lam = 0.0
    if qic > qsmall_value:
        if min_mean_mass_bits != no_limiter_bits:
            nic = min(nic, qic / min_mean_mass)

        lam = micro_mg_utils_basic_lam_native_cb(prop_shape_coef, nic, qic, prop_eff_dim)

        if lam < lambda_bounds1:
            lam = lambda_bounds1
            nic = micro_mg_utils_basic_nic_native_cb(lam, prop_eff_dim, qic, prop_shape_coef)
        elif lam > lambda_bounds2:
            lam = lambda_bounds2
            nic = micro_mg_utils_basic_nic_native_cb(lam, prop_eff_dim, qic, prop_shape_coef)

    nic_out[0] = nic
    lam_out[0] = lam
    if want_n0 != 0:
        n0_out[0] = nic * lam


@export
def size_dist_param_liq_codon(
    qsmall_value: float,
    pi_value: float,
    no_limiter_bits: int,
    prop_rho: float,
    prop_eff_dim: float,
    prop_min_mean_mass: float,
    prop_min_mean_mass_bits: int,
    qcic: float,
    ncic_in: float,
    rho: float,
    ncic_out_c: cobj,
    pgam_c: cobj,
    lamc_c: cobj,
):
    ncic_out = Ptr[float](ncic_out_c)
    pgam_out = Ptr[float](pgam_c)
    lamc_out = Ptr[float](lamc_c)

    ncic = ncic_in
    if qcic > qsmall_value:
        pgam = micro_mg_utils_liq_pgam_native_cb(ncic, rho)
        shape_coef = micro_mg_utils_liq_shape_coef_native_cb(prop_rho, pgam, prop_eff_dim)
        lambda_bounds1 = (pgam + 1.0) * 1.0 / 50.0e-6
        lambda_bounds2 = (pgam + 1.0) * 1.0 / 2.0e-6

        if prop_min_mean_mass_bits != no_limiter_bits:
            ncic = min(ncic, qcic / prop_min_mean_mass)

        lam = micro_mg_utils_basic_lam_native_cb(shape_coef, ncic, qcic, prop_eff_dim)

        if lam < lambda_bounds1:
            lam = lambda_bounds1
            ncic = micro_mg_utils_basic_nic_native_cb(lam, prop_eff_dim, qcic, shape_coef)
        elif lam > lambda_bounds2:
            lam = lambda_bounds2
            ncic = micro_mg_utils_basic_nic_native_cb(lam, prop_eff_dim, qcic, shape_coef)

        pgam_out[0] = pgam
        lamc_out[0] = lam
    else:
        pgam_out[0] = -100.0
        lamc_out[0] = 0.0

    ncic_out[0] = ncic


@export
def avg_diameter_codon(q: float, n: float, rho_air: float, rho_sub: float) -> float:
    return micro_mg_utils_avg_diameter_native_cb(q, n, rho_air, rho_sub)


@export
def radiation_data_flag_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def rad_data_init_codon(flag: int) -> int:
    return radiation_data_flag_codon(flag)


@export
def rad_data_readnl_codon(output: int, fdh: int) -> int:
    out = 0
    if output != 0:
        out += 1
    if fdh != 0:
        out += 2
    return out


@export
def rad_data_register_codon(flag: int) -> int:
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


@inline
def _idx2(r: int, c: int, ld1: int) -> int:
    return (r - 1) + (c - 1) * ld1


@export
def cosp_set_values_tables_codon(
    nht: int,
    nhtml: int,
    nscol: int,
    nprs: int,
    ntau: int,
    ntau_modis: int,
    ndbze: int,
    nsr: int,
    nhtmisr: int,
    zstep: float,
    use_vgrid: int,
    prslim_1d_p: cobj,
    taulim_1d_p: cobj,
    taulim_modis_1d_p: cobj,
    dbzelim_1d_p: cobj,
    srlim_1d_p: cobj,
    htmisrlim_1d_p: cobj,
    htlim_1d_p: cobj,
    prsmid_p: cobj,
    prslim_p: cobj,
    taumid_p: cobj,
    taulim_p: cobj,
    taumid_modis_p: cobj,
    taulim_modis_p: cobj,
    dbzemid_p: cobj,
    dbzelim_p: cobj,
    srmid_p: cobj,
    srlim_p: cobj,
    htmisrmid_p: cobj,
    htmisrlim_p: cobj,
    htmid_p: cobj,
    htlim_p: cobj,
    scol_p: cobj,
    htmlmid_p: cobj,
    prstau_p: cobj,
    prstau_modis_p: cobj,
    htdbze_p: cobj,
    htsr_p: cobj,
    htmlscol_p: cobj,
    htmisrtau_p: cobj,
    prstau_prsmid_p: cobj,
    prstau_taumid_p: cobj,
    prstau_prsmid_modis_p: cobj,
    prstau_taumid_modis_p: cobj,
    htdbze_htmid_p: cobj,
    htdbze_dbzemid_p: cobj,
    htsr_htmid_p: cobj,
    htsr_srmid_p: cobj,
    htmlscol_htmlmid_p: cobj,
    htmlscol_scol_p: cobj,
    htmisrtau_htmisrmid_p: cobj,
    htmisrtau_taumid_p: cobj,
):
    prslim_1d = Ptr[float](prslim_1d_p)
    taulim_1d = Ptr[float](taulim_1d_p)
    taulim_modis_1d = Ptr[float](taulim_modis_1d_p)
    dbzelim_1d = Ptr[float](dbzelim_1d_p)
    srlim_1d = Ptr[float](srlim_1d_p)
    htmisrlim_1d = Ptr[float](htmisrlim_1d_p)
    htlim_1d = Ptr[float](htlim_1d_p)
    prsmid = Ptr[float](prsmid_p)
    prslim = Ptr[float](prslim_p)
    taumid = Ptr[float](taumid_p)
    taulim = Ptr[float](taulim_p)
    taumid_modis = Ptr[float](taumid_modis_p)
    taulim_modis = Ptr[float](taulim_modis_p)
    dbzemid = Ptr[float](dbzemid_p)
    dbzelim = Ptr[float](dbzelim_p)
    srmid = Ptr[float](srmid_p)
    srlim = Ptr[float](srlim_p)
    htmisrmid = Ptr[float](htmisrmid_p)
    htmisrlim = Ptr[float](htmisrlim_p)
    htmid = Ptr[float](htmid_p)
    htlim = Ptr[float](htlim_p)
    scol = Ptr[i32](scol_p)
    htmlmid = Ptr[float](htmlmid_p)
    prstau = Ptr[i32](prstau_p)
    prstau_modis = Ptr[i32](prstau_modis_p)
    htdbze = Ptr[i32](htdbze_p)
    htsr = Ptr[i32](htsr_p)
    htmlscol = Ptr[i32](htmlscol_p)
    htmisrtau = Ptr[i32](htmisrtau_p)
    prstau_prsmid = Ptr[float](prstau_prsmid_p)
    prstau_taumid = Ptr[float](prstau_taumid_p)
    prstau_prsmid_modis = Ptr[float](prstau_prsmid_modis_p)
    prstau_taumid_modis = Ptr[float](prstau_taumid_modis_p)
    htdbze_htmid = Ptr[float](htdbze_htmid_p)
    htdbze_dbzemid = Ptr[float](htdbze_dbzemid_p)
    htsr_htmid = Ptr[float](htsr_htmid_p)
    htsr_srmid = Ptr[float](htsr_srmid_p)
    htmlscol_htmlmid = Ptr[float](htmlscol_htmlmid_p)
    htmlscol_scol = Ptr[float](htmlscol_scol_p)
    htmisrtau_htmisrmid = Ptr[float](htmisrtau_htmisrmid_p)
    htmisrtau_taumid = Ptr[float](htmisrtau_taumid_p)

    if use_vgrid != 0:
        htlim_1d[0] = 0.0
        for i in range(2, nht + 2):
            htlim_1d[i - 1] = float(i - 1) * zstep

    for k in range(1, nprs + 1):
        prsmid[k - 1] = 0.5 * (prslim_1d[k - 1] + prslim_1d[k])
        prslim[_idx2(1, k, 2)] = prslim_1d[k - 1]
        prslim[_idx2(2, k, 2)] = prslim_1d[k]

    for k in range(1, ntau + 1):
        taumid[k - 1] = 0.5 * (taulim_1d[k - 1] + taulim_1d[k])
        taulim[_idx2(1, k, 2)] = taulim_1d[k - 1]
        taulim[_idx2(2, k, 2)] = taulim_1d[k]

    for k in range(1, ntau_modis + 1):
        taumid_modis[k - 1] = 0.5 * (taulim_modis_1d[k - 1] + taulim_modis_1d[k])
        taulim_modis[_idx2(1, k, 2)] = taulim_modis_1d[k - 1]
        taulim_modis[_idx2(2, k, 2)] = taulim_modis_1d[k]

    for k in range(1, ndbze + 1):
        dbzemid[k - 1] = 0.5 * (dbzelim_1d[k - 1] + dbzelim_1d[k])
        dbzelim[_idx2(1, k, 2)] = dbzelim_1d[k - 1]
        dbzelim[_idx2(2, k, 2)] = dbzelim_1d[k]

    for k in range(1, nsr + 1):
        srmid[k - 1] = 0.5 * (srlim_1d[k - 1] + srlim_1d[k])
        srlim[_idx2(1, k, 2)] = srlim_1d[k - 1]
        srlim[_idx2(2, k, 2)] = srlim_1d[k]

    htmisrmid[0] = -99.0
    htmisrlim[_idx2(1, 1, 2)] = htmisrlim_1d[0]
    htmisrlim[_idx2(2, 1, 2)] = htmisrlim_1d[1]
    for k in range(2, nhtmisr + 1):
        htmisrmid[k - 1] = 0.5 * (htmisrlim_1d[k - 1] + htmisrlim_1d[k])
        htmisrlim[_idx2(1, k, 2)] = htmisrlim_1d[k - 1]
        htmisrlim[_idx2(2, k, 2)] = htmisrlim_1d[k]

    for k in range(1, nht + 1):
        htmid[k - 1] = 0.5 * (htlim_1d[k - 1] + htlim_1d[k])
        htlim[_idx2(1, k, 2)] = htlim_1d[k - 1]
        htlim[_idx2(2, k, 2)] = htlim_1d[k]

    for k in range(1, nscol + 1):
        scol[k - 1] = i32(k)

    for k in range(1, nhtml + 1):
        htmlmid[k - 1] = float(k)

    for k in range(1, nprs * ntau + 1):
        prstau[k - 1] = i32(k)
    for k in range(1, nprs * ntau_modis + 1):
        prstau_modis[k - 1] = i32(k)
    for k in range(1, nht * ndbze + 1):
        htdbze[k - 1] = i32(k)
    for k in range(1, nht * nsr + 1):
        htsr[k - 1] = i32(k)
    for k in range(1, nhtml * nscol + 1):
        htmlscol[k - 1] = i32(k)
    for k in range(1, nhtmisr * ntau + 1):
        htmisrtau[k - 1] = i32(k)

    for k in range(1, nprs + 1):
        for j in range(1, ntau + 1):
            idx = j + ntau * (k - 1) - 1
            prstau_taumid[idx] = taumid[j - 1]
            prstau_prsmid[idx] = prsmid[k - 1]
        for j in range(1, ntau_modis + 1):
            idx = j + ntau_modis * (k - 1) - 1
            prstau_taumid_modis[idx] = taumid_modis[j - 1]
            prstau_prsmid_modis[idx] = prsmid[k - 1]

    for k in range(1, nht + 1):
        for j in range(1, ndbze + 1):
            idx = j + ndbze * (k - 1) - 1
            htdbze_dbzemid[idx] = dbzemid[j - 1]
            htdbze_htmid[idx] = htmid[k - 1]

    for k in range(1, nht + 1):
        for j in range(1, nsr + 1):
            idx = j + nsr * (k - 1) - 1
            htsr_srmid[idx] = srmid[j - 1]
            htsr_htmid[idx] = htmid[k - 1]

    for k in range(1, nhtml + 1):
        for j in range(1, nscol + 1):
            idx = j + nscol * (k - 1) - 1
            htmlscol_scol[idx] = float(scol[j - 1])
            htmlscol_htmlmid[idx] = htmlmid[k - 1]

    for k in range(1, nhtmisr + 1):
        for j in range(1, ntau + 1):
            idx = j + ntau * (k - 1) - 1
            htmisrtau_taumid[idx] = taumid[j - 1]
            htmisrtau_htmisrmid[idx] = htmisrmid[k - 1]


@export
def carma_flags_bool_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def carma_flags_touch_codon() -> int:
    return 0


@export
def carma_model_readnl_codon() -> int:
    return 0


@export
def tidal_diag_int_codon(value: int, force_one: int) -> int:
    if force_one != 0:
        return 1
    return value


@export
def get_tidal_coeffs_codon(tod: int, pi: float, cday: float, dcoef_p: cobj):
    dcoef = Ptr[float](dcoef_p)
    gmtfrac = tod / cday
    pi_x_2 = 2.0 * pi
    pi_x_4 = 4.0 * pi

    dcoef[0] = 2.0 * sin(pi_x_2 * gmtfrac)
    dcoef[1] = 2.0 * cos(pi_x_2 * gmtfrac)
    dcoef[2] = 2.0 * sin(pi_x_4 * gmtfrac)
    dcoef[3] = 2.0 * cos(pi_x_4 * gmtfrac)


@export
def tidal_diag_scale_2d_codon(ncol: int, pcols: int, pver: int, coef: float, src_p: cobj, dst_p: cobj):
    src = Ptr[float](src_p)
    dst = Ptr[float](dst_p)
    for k in range(pver):
        base = k * pcols
        for i in range(ncol):
            idx = base + i
            dst[idx] = src[idx] * coef


@export
def tidal_diag_scale_1d_codon(ncol: int, coef: float, src_p: cobj, dst_p: cobj):
    src = Ptr[float](src_p)
    dst = Ptr[float](dst_p)
    for i in range(ncol):
        dst[i] = src[i] * coef


@export
def co2_cycle_flag_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def co2_cycle_readnl_codon(flag: int, read_ocn: int, read_fuel: int) -> int:
    out = 0
    if flag != 0:
        out += 1
    if read_ocn != 0:
        out += 2
    if read_fuel != 0:
        out += 4
    return out


@export
def co2_transport_codon(flag: int) -> int:
    return co2_cycle_flag_codon(flag)


@export
def co2_register_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@inline
def _co2_name_eq(name_len: int, name_ascii: Ptr[int], c0: int, c1: int, c2: int, c3: int, c4: int, c5: int, c6: int) -> int:
    if name_len > 7:
        i = 7
        while i < name_len:
            if name_ascii[i] != 32:
                return 0
            i += 1
    values = (c0, c1, c2, c3, c4, c5, c6)
    i = 0
    while i < 7:
        left = 32
        if i < name_len:
            left = name_ascii[i]
        if left != values[i]:
            return 0
        i += 1
    return 1


@export
def co2_implements_cnst_codon(flag: int, name_len: int, name_ascii_p: cobj) -> int:
    if flag == 0:
        return 0
    name_ascii = Ptr[int](name_ascii_p)
    if _co2_name_eq(name_len, name_ascii, 67, 79, 50, 95, 79, 67, 78) != 0:
        return 1
    if _co2_name_eq(name_len, name_ascii, 67, 79, 50, 95, 70, 70, 70) != 0:
        return 1
    if _co2_name_eq(name_len, name_ascii, 67, 79, 50, 95, 76, 78, 68) != 0:
        return 1
    if _co2_name_eq(name_len, name_ascii, 67, 79, 50, 32, 32, 32, 32) != 0:
        return 1
    return 0


@export
def co2_init_codon(flag: int) -> int:
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
def carma_register_codon() -> int:
    return 0


@export
def carma_implements_cnst_codon() -> int:
    return 0


@export
def carma_init_codon() -> int:
    return 0


@export
def carma_final_codon() -> int:
    return 0


@export
def carma_timestep_init_codon() -> int:
    return 0


@export
def carma_timestep_tend_codon(
    pcols: int,
    prec_str_present: int,
    snow_str_present: int,
    prec_sed_present: int,
    snow_sed_present: int,
    prec_str_p: cobj,
    snow_str_p: cobj,
    prec_sed_p: cobj,
    snow_sed_p: cobj,
):
    if prec_str_present != 0:
        prec_str = Ptr[float](prec_str_p)
        for i in range(pcols):
            prec_str[i] = 0.0
    if snow_str_present != 0:
        snow_str = Ptr[float](snow_str_p)
        for i in range(pcols):
            snow_str[i] = 0.0
    if prec_sed_present != 0:
        prec_sed = Ptr[float](prec_sed_p)
        for i in range(pcols):
            prec_sed[i] = 0.0
    if snow_sed_present != 0:
        snow_sed = Ptr[float](snow_sed_p)
        for i in range(pcols):
            snow_sed[i] = 0.0


@export
def carma_accumulate_stats_codon() -> int:
    return 0


@export
def subcol_touch_codon() -> int:
    return 0


@export
def subcol_readnl_codon() -> int:
    return 0


@export
def subcol_register_codon() -> int:
    return 0


@export
def subcol_init_codon() -> int:
    return 0


@export
def cldwat_param_codon(value: float) -> float:
    return value


@inline
def _cldwat_gaussj_idx(row: int, col: int, ld1: int) -> int:
    """gaussj a(np,np) and b(np,mp) use Fortran column-major order."""
    return (row - 1) + (col - 1) * ld1


@export
def gaussj_codon(
    n: int,
    np: int,
    m: int,
    mp: int,
    a_p: cobj,
    b_p: cobj,
    indxc_p: cobj,
    indxr_p: cobj,
    ipiv_p: cobj,
) -> int:
    a = Ptr[float](a_p)
    b = Ptr[float](b_p)
    indxc = Ptr[int](indxc_p)
    indxr = Ptr[int](indxr_p)
    ipiv = Ptr[int](ipiv_p)

    for j in range(1, n + 1):
        ipiv[j - 1] = 0

    for i in range(1, n + 1):
        big = 0.0
        irow = 1
        icol = 1
        for j in range(1, n + 1):
            if ipiv[j - 1] != 1:
                for k in range(1, n + 1):
                    if ipiv[k - 1] == 0:
                        value = abs(a[_cldwat_gaussj_idx(j, k, np)])
                        if value >= big:
                            big = value
                            irow = j
                            icol = k
                    elif ipiv[k - 1] > 1:
                        return 1

        ipiv[icol - 1] = ipiv[icol - 1] + 1
        if irow != icol:
            for l in range(1, n + 1):
                idx_row = _cldwat_gaussj_idx(irow, l, np)
                idx_col = _cldwat_gaussj_idx(icol, l, np)
                dum = a[idx_row]
                a[idx_row] = a[idx_col]
                a[idx_col] = dum
            for l in range(1, m + 1):
                idx_row = _cldwat_gaussj_idx(irow, l, np)
                idx_col = _cldwat_gaussj_idx(icol, l, np)
                dum = b[idx_row]
                b[idx_row] = b[idx_col]
                b[idx_col] = dum

        indxr[i - 1] = irow
        indxc[i - 1] = icol
        pivot_idx = _cldwat_gaussj_idx(icol, icol, np)
        if a[pivot_idx] == 0.0:
            return 2

        pivinv = 1.0 / a[pivot_idx]
        a[pivot_idx] = 1.0
        for l in range(1, n + 1):
            idx = _cldwat_gaussj_idx(icol, l, np)
            a[idx] = a[idx] * pivinv
        for l in range(1, m + 1):
            idx = _cldwat_gaussj_idx(icol, l, np)
            b[idx] = b[idx] * pivinv

        for ll in range(1, n + 1):
            if ll != icol:
                ll_icol_idx = _cldwat_gaussj_idx(ll, icol, np)
                dum = a[ll_icol_idx]
                a[ll_icol_idx] = 0.0
                for l in range(1, n + 1):
                    ll_l_idx = _cldwat_gaussj_idx(ll, l, np)
                    a[ll_l_idx] = a[ll_l_idx] - a[_cldwat_gaussj_idx(icol, l, np)] * dum
                for l in range(1, m + 1):
                    ll_l_idx = _cldwat_gaussj_idx(ll, l, np)
                    b[ll_l_idx] = b[ll_l_idx] - b[_cldwat_gaussj_idx(icol, l, np)] * dum

    for l in range(n, 0, -1):
        row = indxr[l - 1]
        col = indxc[l - 1]
        if row != col:
            for k in range(1, n + 1):
                row_idx = _cldwat_gaussj_idx(k, row, np)
                col_idx = _cldwat_gaussj_idx(k, col, np)
                dum = a[row_idx]
                a[row_idx] = a[col_idx]
                a[col_idx] = dum

    return 0


@export
def cldwat2m_gaussj_codon(
    n: int,
    np: int,
    m: int,
    mp: int,
    a_p: cobj,
    b_p: cobj,
    indxc_p: cobj,
    indxr_p: cobj,
    ipiv_p: cobj,
) -> int:
    return gaussj_codon(n, np, m, mp, a_p, b_p, indxc_p, indxr_p, ipiv_p)


@export
def hkconv_param_codon(value: float) -> float:
    return value


@export
def cam3_aero_data_flag_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def cam3_aero_data_readnl_codon(flag: int) -> int:
    return cam3_aero_data_flag_codon(flag)


@export
def cam3_ozone_data_flag_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def cam3_ozone_data_readnl_codon(flag: int, cyc: int) -> int:
    out = 0
    if flag != 0:
        out += 1
    if cyc != 0:
        out += 2
    return out


@export
def phys_debug_value_codon(value: float) -> float:
    return value


@export
def phys_debug_readnl_codon(lat_set: int, lon_set: int) -> int:
    if lat_set != 0 and lon_set != 0:
        return 1
    return 0


@export
def phys_debug_has_location_codon(lat_set: int, lon_set: int) -> int:
    if lat_set != 0 and lon_set != 0:
        return 1
    return 0


@export
def phys_debug_init_codon(lat_set: int, lon_set: int) -> int:
    if lat_set != 0 and lon_set != 0:
        return 1
    return 0


@export
def phys_debug_col_codon(chunk: int, debchunk: int, debcol: int) -> int:
    if chunk == debchunk:
        return debcol
    return 0


@export
def unicon_cam_flag_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def unicon_cam_readnl_codon(flag: int, hfile: int) -> int:
    out = hfile * 2
    if flag != 0:
        out += 1
    return out


@export
def unicon_cam_int_codon(value: int) -> int:
    return value


@export
def iondrag_touch_codon() -> int:
    return 0


@export
def iondrag_readnl_codon() -> int:
    return 0


@export
def iondrag_register_codon() -> int:
    return 0


@export
def iondrag_init_codon() -> int:
    return 0


@export
def iondrag_calc_ghg_codon() -> int:
    return 0


@export
def cld_sediment_param_codon(value: float) -> float:
    return value


@export
def tsinti_param_codon(value: float) -> float:
    return value


@export
def tsinti_codon(
    tmelt: float,
    latvap: float,
    rair: float,
    stebol: float,
    latice: float,
    values_p: cobj,
):
    values = Ptr[float](values_p)
    values[0] = latice
    values[1] = tmelt
    values[2] = latvap
    values[3] = rair
    values[4] = stebol
    values[5] = 10.0


@export
def clubb_intr_flag_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def clubb_intr_touch_codon() -> int:
    return 0


@export
def clubb_implements_cnst_codon(do_cnst: int, name_len: int, name_ascii_p: cobj) -> int:
    if do_cnst == 0:
        return 0
    name_ascii = Ptr[int](name_ascii_p)
    if _clubb_name_eq(name_len, name_ascii, 84, 72, 76, 80, 50, 32, 32, 32) != 0:
        return 1
    if _clubb_name_eq(name_len, name_ascii, 82, 84, 80, 50, 32, 32, 32, 32) != 0:
        return 1
    if _clubb_name_eq(name_len, name_ascii, 82, 84, 80, 84, 72, 76, 80, 32) != 0:
        return 1
    if _clubb_name_eq(name_len, name_ascii, 87, 80, 84, 72, 76, 80, 32, 32) != 0:
        return 1
    if _clubb_name_eq(name_len, name_ascii, 87, 80, 82, 84, 80, 32, 32, 32) != 0:
        return 1
    if _clubb_name_eq(name_len, name_ascii, 87, 80, 50, 32, 32, 32, 32, 32) != 0:
        return 1
    if _clubb_name_eq(name_len, name_ascii, 87, 80, 51, 32, 32, 32, 32, 32) != 0:
        return 1
    if _clubb_name_eq(name_len, name_ascii, 85, 80, 50, 32, 32, 32, 32, 32) != 0:
        return 1
    if _clubb_name_eq(name_len, name_ascii, 86, 80, 50, 32, 32, 32, 32, 32) != 0:
        return 1
    return 0


@inline
def _clubb_name_eq(
    name_len: int,
    name_ascii: Ptr[int],
    c0: int,
    c1: int,
    c2: int,
    c3: int,
    c4: int,
    c5: int,
    c6: int,
    c7: int,
) -> int:
    if name_len > 8:
        i = 8
        while i < name_len:
            if name_ascii[i] != 32:
                return 0
            i += 1
    values = (c0, c1, c2, c3, c4, c5, c6, c7)
    i = 0
    while i < 8:
        left = 32
        if i < name_len:
            left = name_ascii[i]
        if left != values[i]:
            return 0
        i += 1
    return 1


@export
def clubb_readnl_codon() -> int:
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
def hetfrz_classnuc_cam_readnl_codon(flag: int) -> int:
    return hetfrz_classnuc_cam_flag_codon(flag)


@export
def hetfrz_classnuc_cam_register_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def hetfrz_classnuc_cam_init_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def press_lim_idx_codon(p: float, top: int, pver: int, pref_mid_p: cobj) -> int:
    pref_mid = Ptr[float](pref_mid_p)
    if top != 0:
        k_lim = pver + 1
        for k in range(1, pver + 1):
            if pref_mid[k - 1] > p:
                k_lim = k
                break
        return k_lim

    k_lim = 0
    for k in range(pver, 0, -1):
        if pref_mid[k - 1] < p:
            k_lim = k
            break
    return k_lim


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
def sslt_rebin_active_codon(has_sslt: int) -> int:
    if has_sslt != 0:
        return 1
    return 0


@export
def sslt_rebin_register_codon(pcols: int, pver: int) -> int:
    return pcols * 1000 + pver


@export
def sslt_rebin_adv_codon(
    ncol: int,
    pver: int,
    pcols: int,
    wgt_sscm: float,
    sslt1_p: cobj,
    sslt2_p: cobj,
    sslt3_p: cobj,
    sslt4_p: cobj,
    sslta_p: cobj,
    ssltc_p: cobj,
):
    sslt1 = Ptr[float](sslt1_p)
    sslt2 = Ptr[float](sslt2_p)
    sslt3 = Ptr[float](sslt3_p)
    sslt4 = Ptr[float](sslt4_p)
    sslta = Ptr[float](sslta_p)
    ssltc = Ptr[float](ssltc_p)
    accum_wgt = 1.0 - wgt_sscm
    for k in range(pver):
        base = k * pcols
        for i in range(ncol):
            idx = base + i
            sslt_sum = sslt1[idx] + sslt2[idx] + sslt3[idx] + sslt4[idx]
            sslta[idx] = accum_wgt * sslt_sum
            ssltc[idx] = wgt_sscm * sslt_sum


@export
def constituent_burden_flag_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def constituent_burden_init_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def constituent_burden_comp_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def constituent_burden_comp_integral_codon(
    ncol: int,
    psetcols: int,
    pver: int,
    pcnst: int,
    m: int,
    dry: int,
    rga: float,
    q_p: cobj,
    pdel_p: cobj,
    pdeldry_p: cobj,
    ftem_p: cobj,
):
    q = Ptr[float](q_p)
    pdel = Ptr[float](pdel_p)
    pdeldry = Ptr[float](pdeldry_p)
    ftem = Ptr[float](ftem_p)
    m0 = m - 1
    q_m_base = m0 * psetcols * pver

    if dry != 0:
        for i in range(ncol):
            acc = 0.0
            for k in range(pver):
                idx2 = k * psetcols + i
                acc += q[q_m_base + idx2] * pdeldry[idx2]
            ftem[i] = acc * rga
    else:
        for i in range(ncol):
            acc = 0.0
            for k in range(pver):
                idx2 = k * psetcols + i
                acc += q[q_m_base + idx2] * pdel[idx2]
            ftem[i] = acc * rga


@export
def ghg_data_register_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def ghg_data_timestep_init_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def gmean_fixed_repro_codon(arr_p: cobj, nflds: int, pi_value: float):
    arr = Ptr[float](arr_p)
    denom = 4.0 * pi_value
    for i in range(nflds):
        arr[i] = arr[i] / denom


@export
def gmean_fixed_repro_preweight_chunk_codon(
    pcols: int,
    begchunk: int,
    endchunk: int,
    lchnk: int,
    ncols: int,
    ifld: int,
    nlcols: int,
    count: int,
    arr_p: cobj,
    wght_p: cobj,
    xfld_p: cobj,
) -> int:
    arr = Ptr[float](arr_p)
    wght = Ptr[float](wght_p)
    xfld = Ptr[float](xfld_p)

    field0 = ifld - 1
    chunk0 = lchnk - begchunk
    chunk_count = endchunk - begchunk + 1
    arr_base = field0 * pcols * chunk_count + chunk0 * pcols
    xfld_base = field0 * nlcols

    for i0 in range(ncols):
        xfld[xfld_base + count] = arr[arr_base + i0] * wght[i0]
        count += 1

    return count


@export
def gmean_arr_recompute_flags_codon(
    rel_diff_p: cobj,
    nflds: int,
    reldiffmax: float,
    recompute: int,
    recompute_flags_p: cobj,
):
    rel_diff = Ptr[float](rel_diff_p)
    recompute_flags = Ptr[int](recompute_flags_p)

    for ifld in range(nflds):
        flag = 0
        if recompute != 0 and rel_diff[ifld * 2] > reldiffmax:
            flag = 1
        recompute_flags[ifld] = flag


@export
def phys_gmean_normalize_codon(arr_p: cobj, nflds: int, pi_value: float):
    gmean_fixed_repro_codon(arr_p, nflds, pi_value)


@export
def pbl_utils_value_codon(value: float) -> float:
    return value


@export
def pbl_utils_init_codon(value: float) -> float:
    return value


@export
def calc_ustar_rrho_codon(rair: float, t: float, pmid: float) -> float:
    return rair * t / pmid


@export
def calc_ustar_codon(taux: float, tauy: float, rrho: float, ustar_min: float) -> float:
    return max(sqrt(sqrt(taux**2 + tauy**2) * rrho), ustar_min)


@export
def calc_obklen_khfs_codon(shflx: float, rrho: float, cpair: float) -> float:
    return shflx * rrho / cpair


@export
def calc_obklen_kqfs_codon(qflx: float, rrho: float) -> float:
    return qflx * rrho


@export
def calc_obklen_kbfs_codon(khfs: float, zvir: float, ths: float, kqfs: float) -> float:
    return khfs + zvir * ths * kqfs


@export
def calc_obklen_codon(thvs: float, ustar: float, g: float, vk: float, kbfs: float) -> float:
    sign_eps = 1.0e-10
    if kbfs < 0.0:
        sign_eps = -1.0e-10
    return -thvs * ustar**3 / (g * vk * (kbfs + sign_eps))


@export
def virtem_codon(t: float, q: float, zvir: float) -> float:
    return t * (1.0 + zvir * q)


@export
def compute_radf_codon(
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

    for ncv0 in range(0, ncvmax):
        lwp_CL[ncv0] = 0.0
        opt_depth_CL[ncv0] = 0.0
        radinvfrac_CL[ncv0] = 0.0
        radf_CL[ncv0] = 0.0

    i = i_col
    ncvfin_i = int(ncvfin[i - 1])

    for ncv in range(1, ncvfin_i + 1):
        kt = int(ktop[(i - 1) + (ncv - 1) * pcols])

        lwp = 0.0
        opt_depth = 0.0
        radinvfrac = 0.0
        radf = 0.0

        if radf_mode == 0:
            if ql[_field2_idx(i, kt, pcols)] > qmin and ql[_field2_idx(i, kt - 1, pcols)] < qmin:
                lwp = ql[_field2_idx(i, kt, pcols)] * (pi[_field2_idx(i, kt + 1, pcols)] - pi[_field2_idx(i, kt, pcols)]) / gravit
                opt_depth = 156.0 * lwp
                radinvfrac = opt_depth * (4.0 + opt_depth) / (6.0 * (4.0 + opt_depth) + opt_depth**2)
                radf = qrlw[_field2_idx(i, kt, pcols)] / (pi[_field2_idx(i, kt, pcols)] - pi[_field2_idx(i, kt + 1, pcols)])
                radf = max(
                    radinvfrac * radf * (zi[_field2_idx(i, kt, pcols)] - zi[_field2_idx(i, kt + 1, pcols)]),
                    0.0,
                ) * chs[_field2_idx(i, kt, pcols)]

        elif radf_mode == 1:
            lwp = ql[_field2_idx(i, kt, pcols)] * (pi[_field2_idx(i, kt + 1, pcols)] - pi[_field2_idx(i, kt, pcols)]) / gravit
            opt_depth = 156.0 * lwp
            radinvfrac = opt_depth * (4.0 + opt_depth) / (6.0 * (4.0 + opt_depth) + opt_depth**2)
            radinvfrac = max(cldeff[_field2_idx(i, kt, pcols)] - cldeff[_field2_idx(i, kt - 1, pcols)], 0.0) * radinvfrac
            radf = qrlw[_field2_idx(i, kt, pcols)] / (pi[_field2_idx(i, kt, pcols)] - pi[_field2_idx(i, kt + 1, pcols)])
            radf = max(
                radinvfrac * radf * (zi[_field2_idx(i, kt, pcols)] - zi[_field2_idx(i, kt + 1, pcols)]),
                0.0,
            ) * chs[_field2_idx(i, kt, pcols)]

        else:
            lwp = ql[_field2_idx(i, kt, pcols)] * (pi[_field2_idx(i, kt + 1, pcols)] - pi[_field2_idx(i, kt, pcols)]) / gravit
            opt_depth = 156.0 * lwp
            radinvfrac = opt_depth * (4.0 + opt_depth) / (6.0 * (4.0 + opt_depth) + opt_depth**2)
            radf = max(
                radinvfrac
                * qrlw[_field2_idx(i, kt, pcols)]
                / (pi[_field2_idx(i, kt, pcols)] - pi[_field2_idx(i, kt + 1, pcols)])
                * (zi[_field2_idx(i, kt, pcols)] - zi[_field2_idx(i, kt + 1, pcols)]),
                0.0,
            )

            lwp = ql[_field2_idx(i, kt - 1, pcols)] * (pi[_field2_idx(i, kt, pcols)] - pi[_field2_idx(i, kt - 1, pcols)]) / gravit
            opt_depth = 156.0 * lwp
            radinvfrac = opt_depth * (4.0 + opt_depth) / (6.0 * (4.0 + opt_depth) + opt_depth**2)
            radf = radf + max(
                radinvfrac
                * qrlw[_field2_idx(i, kt - 1, pcols)]
                / (pi[_field2_idx(i, kt - 1, pcols)] - pi[_field2_idx(i, kt, pcols)])
                * (zi[_field2_idx(i, kt - 1, pcols)] - zi[_field2_idx(i, kt, pcols)]),
                0.0,
            )
            radf = max(radf, 0.0) * chs[_field2_idx(i, kt, pcols)]

        lwp_CL[ncv - 1] = lwp
        opt_depth_CL[ncv - 1] = opt_depth
        radinvfrac_CL[ncv - 1] = radinvfrac
        radf_CL[ncv - 1] = radf


@export
def pbl_utils_compute_radf_codon(
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
    compute_radf_codon(
        i_col,
        pcols,
        pver,
        ncvmax,
        radf_mode,
        qmin,
        gravit,
        ncvfin_p,
        ktop_p,
        ql_p,
        pi_p,
        qrlw_p,
        cldeff_p,
        zi_p,
        chs_p,
        lwp_CL_p,
        opt_depth_CL_p,
        radinvfrac_CL_p,
        radf_CL_p,
    )


@export
def wv_sat_methods_value_codon(value: float) -> float:
    return value


@export
def wv_sat_methods_omeps_codon(epsilo: float) -> float:
    return 1.0 - epsilo


@export
def wv_sat_methods_init_codon(
    tmelt_in: float,
    h2otrip_in: float,
    tboil_in: float,
    ttrice_in: float,
    epsilo_in: float,
    tmelt_p: cobj,
    h2otrip_p: cobj,
    tboil_p: cobj,
    ttrice_p: cobj,
    epsilo_p: cobj,
    omeps_p: cobj,
):
    tmelt = Ptr[float](tmelt_p)
    h2otrip = Ptr[float](h2otrip_p)
    tboil = Ptr[float](tboil_p)
    ttrice = Ptr[float](ttrice_p)
    epsilo = Ptr[float](epsilo_p)
    omeps = Ptr[float](omeps_p)

    tmelt[0] = tmelt_in
    h2otrip[0] = h2otrip_in
    tboil[0] = tboil_in
    ttrice[0] = ttrice_in
    epsilo[0] = epsilo_in
    omeps[0] = 1.0 - epsilo_in


@export
def wv_sat_valid_idx_codon(idx: int) -> int:
    if idx != -1:
        return 1
    return 0


@export
def wv_sat_svp_to_qsat_codon(es: float, p: float, epsilo: float, omeps: float) -> float:
    if (p - es) <= 0.0:
        return 1.0
    return epsilo * es / (p - omeps * es)


WV_SAT_INVALID_IDX = -1
WV_SAT_OLD_GOFF_GRATCH_IDX = 0
WV_SAT_GOFF_GRATCH_IDX = 1
WV_SAT_MURPHY_KOOP_IDX = 2
WV_SAT_BOLTON_IDX = 3


@export
def wv_sat_get_scheme_idx_codon(name_len: int, name_ascii_p: cobj) -> int:
    name_ascii = Ptr[int](name_ascii_p)

    n = name_len
    while n > 0 and name_ascii[n - 1] == 32:
        n -= 1

    if n == 10:
        if (
            name_ascii[0] == 71
            and name_ascii[1] == 111
            and name_ascii[2] == 102
            and name_ascii[3] == 102
            and name_ascii[4] == 71
            and name_ascii[5] == 114
            and name_ascii[6] == 97
            and name_ascii[7] == 116
            and name_ascii[8] == 99
            and name_ascii[9] == 104
        ):
            return WV_SAT_GOFF_GRATCH_IDX
        if (
            name_ascii[0] == 77
            and name_ascii[1] == 117
            and name_ascii[2] == 114
            and name_ascii[3] == 112
            and name_ascii[4] == 104
            and name_ascii[5] == 121
            and name_ascii[6] == 75
            and name_ascii[7] == 111
            and name_ascii[8] == 111
            and name_ascii[9] == 112
        ):
            return WV_SAT_MURPHY_KOOP_IDX
    elif n == 13:
        if (
            name_ascii[0] == 79
            and name_ascii[1] == 108
            and name_ascii[2] == 100
            and name_ascii[3] == 71
            and name_ascii[4] == 111
            and name_ascii[5] == 102
            and name_ascii[6] == 102
            and name_ascii[7] == 71
            and name_ascii[8] == 114
            and name_ascii[9] == 97
            and name_ascii[10] == 116
            and name_ascii[11] == 99
            and name_ascii[12] == 104
        ):
            return WV_SAT_OLD_GOFF_GRATCH_IDX
    elif n == 6:
        if (
            name_ascii[0] == 66
            and name_ascii[1] == 111
            and name_ascii[2] == 108
            and name_ascii[3] == 116
            and name_ascii[4] == 111
            and name_ascii[5] == 110
        ):
            return WV_SAT_BOLTON_IDX

    return WV_SAT_INVALID_IDX


@export
def wv_sat_set_default_codon(tmp_idx: int, default_idx_p: cobj) -> int:
    default_idx = Ptr[i32](default_idx_p)

    if wv_sat_valid_idx_codon(tmp_idx) == 0:
        return 0

    default_idx[0] = i32(tmp_idx)
    return 1


@export
def goffgratch_svp_water_codon(t: float) -> float:
    return goffgratch_svp_water_native_cb(t)


@export
def goffgratch_svp_ice_codon(t: float) -> float:
    return goffgratch_svp_ice_native_cb(t)


@export
def murphykoop_svp_water_codon(t: float) -> float:
    return murphykoop_svp_water_native_cb(t)


@export
def murphykoop_svp_ice_codon(t: float) -> float:
    return murphykoop_svp_ice_native_cb(t)


@export
def oldgoffgratch_svp_water_codon(t: float) -> float:
    return oldgoffgratch_svp_water_native_cb(t)


@export
def oldgoffgratch_svp_ice_codon(t: float) -> float:
    return oldgoffgratch_svp_ice_native_cb(t)


@export
def bolton_svp_water_codon(t: float) -> float:
    return bolton_svp_water_native_cb(t)


@export
def wv_sat_svp_water_codon(t: float, idx: int) -> float:
    if idx == WV_SAT_GOFF_GRATCH_IDX:
        return goffgratch_svp_water_codon(t)
    if idx == WV_SAT_MURPHY_KOOP_IDX:
        return murphykoop_svp_water_codon(t)
    if idx == WV_SAT_OLD_GOFF_GRATCH_IDX:
        return oldgoffgratch_svp_water_codon(t)
    if idx == WV_SAT_BOLTON_IDX:
        return bolton_svp_water_codon(t)
    return 0.0


@export
def wv_sat_svp_ice_codon(t: float, idx: int) -> float:
    if idx == WV_SAT_GOFF_GRATCH_IDX:
        return goffgratch_svp_ice_codon(t)
    if idx == WV_SAT_MURPHY_KOOP_IDX:
        return murphykoop_svp_ice_codon(t)
    if idx == WV_SAT_OLD_GOFF_GRATCH_IDX:
        return oldgoffgratch_svp_ice_codon(t)
    if idx == WV_SAT_BOLTON_IDX:
        return bolton_svp_water_codon(t)
    return 0.0


@export
def wv_sat_svp_trans_codon(t: float, idx: int, tmelt: float, ttrice: float) -> float:
    es = 0.0
    if t >= (tmelt - ttrice):
        es = wv_sat_svp_water_codon(t, idx)

    if t < tmelt:
        esice = wv_sat_svp_ice_codon(t, idx)
        weight = 0.0
        if (tmelt - t) > ttrice:
            weight = 1.0
        else:
            weight = (tmelt - t) / ttrice
        es = weight * esice + (1.0 - weight) * es

    return es


@export
def wv_sat_init_table_codon(plenest: int, tmin: float, tmelt: float, ttrice: float, idx: int, estbl_p: cobj):
    estbl = Ptr[float](estbl_p)
    for i in range(plenest):
        estbl[i] = wv_sat_svp_trans_codon(tmin + float(i), idx, tmelt, ttrice)


@export
def svp_water_codon(t: float, idx: int) -> float:
    return wv_sat_svp_water_codon(t, idx)


@export
def svp_ice_codon(t: float, idx: int) -> float:
    return wv_sat_svp_ice_codon(t, idx)


@export
def svp_trans_codon(t: float, idx: int, tmelt: float, ttrice: float) -> float:
    return wv_sat_svp_trans_codon(t, idx, tmelt, ttrice)


@export
def estblf_codon(t: float, tmin: float, tmax: float, estbl_p: cobj, plenest: int) -> float:
    estbl = Ptr[float](estbl_p)

    t_limited = t
    if t_limited > tmax:
        t_limited = tmax

    t_tmp = t_limited - tmin
    if t_tmp < 0.0:
        t_tmp = 0.0

    i0 = int(t_tmp)
    if i0 < 0:
        i0 = 0
    if i0 > plenest - 2:
        i0 = plenest - 2

    weight = t_tmp - float(i0)
    return (1.0 - weight) * estbl[i0] + weight * estbl[i0 + 1]


@export
def wv_sat_qsat_water_codon(
    t: float,
    p: float,
    idx: int,
    epsilo: float,
    omeps: float,
    es_p: cobj,
    qs_p: cobj,
):
    es_out = Ptr[float](es_p)
    qs_out = Ptr[float](qs_p)

    es = wv_sat_svp_water_codon(t, idx)
    qs = wv_sat_svp_to_qsat_codon(es, p, epsilo, omeps)
    if p < es:
        es = p

    es_out[0] = es
    qs_out[0] = qs


@inline
def _zm_qsat_hpa_ptr_codon(
    t: float,
    p_hpa: float,
    idx: int,
    epsilo: float,
    omeps: float,
    es_out: Ptr[float],
    qm_out: Ptr[float],
):
    p_pa = p_hpa * 100.0
    es = wv_sat_svp_water_codon(t, idx)
    qm = wv_sat_svp_to_qsat_codon(es, p_pa, epsilo, omeps)
    if p_pa < es:
        es = p_pa

    es_out[0] = es * 0.01
    qm_out[0] = qm


@export
def zm_qsat_hpa_codon(
    t: float,
    p_hpa: float,
    idx: int,
    epsilo: float,
    omeps: float,
    es_p: cobj,
    qm_p: cobj,
):
    _zm_qsat_hpa_ptr_codon(t, p_hpa, idx, epsilo, omeps, Ptr[float](es_p), Ptr[float](qm_p))


@export
def qsat_hpa_codon(
    t: float,
    p_hpa: float,
    idx: int,
    epsilo: float,
    omeps: float,
    es_p: cobj,
    qm_p: cobj,
):
    _zm_qsat_hpa_ptr_codon(t, p_hpa, idx, epsilo, omeps, Ptr[float](es_p), Ptr[float](qm_p))


@export
def zm_entropy_codon(
    tk: float,
    p_hpa: float,
    qtot: float,
    rl: float,
    cpliq: float,
    cpwv: float,
    tfreez: float,
    cpres: float,
    rgas: float,
    eps1: float,
    rh2o: float,
    idx: int,
    epsilo: float,
    omeps: float,
) -> float:
    est = 0.0
    qst = 0.0
    _zm_qsat_hpa_ptr_codon(tk, p_hpa, idx, epsilo, omeps, __ptr__(est), __ptr__(qst))
    return zm_entropy_expr_native_cb(
        tk,
        p_hpa,
        qtot,
        qst,
        rl,
        cpliq,
        cpwv,
        tfreez,
        cpres,
        rgas,
        eps1,
        rh2o,
    )


@export
def entropy_codon(
    tk: float,
    p_hpa: float,
    qtot: float,
    rl: float,
    cpliq: float,
    cpwv: float,
    tfreez: float,
    cpres: float,
    rgas: float,
    eps1: float,
    rh2o: float,
    idx: int,
    epsilo: float,
    omeps: float,
) -> float:
    return zm_entropy_codon(
        tk,
        p_hpa,
        qtot,
        rl,
        cpliq,
        cpwv,
        tfreez,
        cpres,
        rgas,
        eps1,
        rh2o,
        idx,
        epsilo,
        omeps,
    )


@export
def ientropy_codon(
    s: float,
    p_hpa: float,
    qt: float,
    tfg: float,
    rl: float,
    cpliq: float,
    cpwv: float,
    tfreez: float,
    cpres: float,
    rgas: float,
    eps1: float,
    rh2o: float,
    idx: int,
    epsilo: float,
    omeps: float,
    t_p: cobj,
    qst_p: cobj,
    converged_p: cobj,
):
    t_out = Ptr[float](t_p)
    qst_out = Ptr[float](qst_p)
    converged_out = Ptr[int](converged_p)

    loopmax = 100
    eps = 3.0e-8
    tol = 0.001

    a = tfg - 10.0
    b = tfg + 10.0
    fa = zm_entropy_codon(a, p_hpa, qt, rl, cpliq, cpwv, tfreez, cpres, rgas, eps1, rh2o, idx, epsilo, omeps) - s
    fb = zm_entropy_codon(b, p_hpa, qt, rl, cpliq, cpwv, tfreez, cpres, rgas, eps1, rh2o, idx, epsilo, omeps) - s
    c = b
    fc = fb
    d = 0.0
    ebr = 0.0
    converged = False

    for _ in range(0, loopmax + 1):
        if (fb > 0.0 and fc > 0.0) or (fb < 0.0 and fc < 0.0):
            c = a
            fc = fa
            d = b - a
            ebr = d

        if abs(fc) < abs(fb):
            a = b
            b = c
            c = a
            fa = fb
            fb = fc
            fc = fa

        tol1 = 2.0 * eps * abs(b) + 0.5 * tol
        xm = 0.5 * (c - b)
        converged = abs(xm) <= tol1 or fb == 0.0
        if converged:
            break

        if abs(ebr) >= tol1 and abs(fa) > abs(fb):
            sbr = fb / fa
            if a == c:
                pbr = 2.0 * xm * sbr
                qbr = 1.0 - sbr
            else:
                qbr = fa / fc
                rbr = fb / fc
                pbr = sbr * (2.0 * xm * qbr * (qbr - rbr) - (b - a) * (rbr - 1.0))
                qbr = (qbr - 1.0) * (rbr - 1.0) * (sbr - 1.0)

            if pbr > 0.0:
                qbr = -qbr
            pbr = abs(pbr)
            if 2.0 * pbr < min(3.0 * xm * qbr - abs(tol1 * qbr), abs(ebr * qbr)):
                ebr = d
                d = pbr / qbr
            else:
                d = xm
                ebr = d
        else:
            d = xm
            ebr = d

        a = b
        fa = fb
        if abs(d) > tol1:
            b = b + d
        else:
            if xm >= 0.0:
                b = b + abs(tol1)
            else:
                b = b - abs(tol1)

        fb = zm_entropy_codon(b, p_hpa, qt, rl, cpliq, cpwv, tfreez, cpres, rgas, eps1, rh2o, idx, epsilo, omeps) - s

    est = 0.0
    qst = 0.0
    _zm_qsat_hpa_ptr_codon(b, p_hpa, idx, epsilo, omeps, __ptr__(est), __ptr__(qst))
    t_out[0] = b
    qst_out[0] = qst
    if converged:
        converged_out[0] = 1
    else:
        converged_out[0] = 0


@export
def zm_convi_codon(
    limcnv_in: int,
    zmconv_c0_lnd: float,
    zmconv_c0_ocn: float,
    zmconv_ke: float,
    zmconv_ke_lnd: float,
    zmconv_org: int,
    no_deep_present: int,
    no_deep_value: int,
    tmelt: float,
    epsilo: float,
    latvap: float,
    cpair: float,
    gravit: float,
    rair: float,
    limcnv_p: cobj,
    tfreez_p: cobj,
    eps1_p: cobj,
    rl_p: cobj,
    cpres_p: cobj,
    rgrav_p: cobj,
    rgas_p: cobj,
    grav_p: cobj,
    cp_p: cobj,
    c0_lnd_p: cobj,
    c0_ocn_p: cobj,
    ke_p: cobj,
    ke_lnd_p: cobj,
    zm_org_p: cobj,
    no_deep_pbl_p: cobj,
    tau_p: cobj,
):
    limcnv = Ptr[int](limcnv_p)
    tfreez = Ptr[float](tfreez_p)
    eps1 = Ptr[float](eps1_p)
    rl = Ptr[float](rl_p)
    cpres = Ptr[float](cpres_p)
    rgrav = Ptr[float](rgrav_p)
    rgas = Ptr[float](rgas_p)
    grav = Ptr[float](grav_p)
    cp = Ptr[float](cp_p)
    c0_lnd = Ptr[float](c0_lnd_p)
    c0_ocn = Ptr[float](c0_ocn_p)
    ke = Ptr[float](ke_p)
    ke_lnd = Ptr[float](ke_lnd_p)
    zm_org_out = Ptr[int](zm_org_p)
    no_deep_pbl = Ptr[int](no_deep_pbl_p)
    tau = Ptr[float](tau_p)

    limcnv[0] = limcnv_in
    tfreez[0] = tmelt
    eps1[0] = epsilo
    rl[0] = latvap
    cpres[0] = cpair
    rgrav[0] = 1.0 / gravit
    rgas[0] = rair
    grav[0] = gravit
    cp[0] = cpres[0]
    c0_lnd[0] = zmconv_c0_lnd
    c0_ocn[0] = zmconv_c0_ocn
    ke[0] = zmconv_ke
    ke_lnd[0] = zmconv_ke_lnd
    zm_org_out[0] = zmconv_org
    if no_deep_present != 0:
        no_deep_pbl[0] = no_deep_value
    else:
        no_deep_pbl[0] = 0
    tau[0] = 3600.0


@export
def wv_saturation_value_codon(value: float) -> float:
    return value


@export
def wv_sat_readnl_codon() -> int:
    return 1


@export
def wv_saturation_touch_codon(stage: int) -> int:
    return stage


@export
def wv_saturation_limit_es_codon(es: float, p: float) -> float:
    if p < es:
        return p
    return es


@export
def wv_saturation_tq_enthalpy_codon(cpair: float, t: float, q: float, hltalt: float) -> float:
    return cpair * t + hltalt * q


@export
def wv_saturation_no_ip_hltalt_codon(t: float, tmelt: float, latvap: float) -> float:
    hltalt = latvap
    if t >= tmelt:
        hltalt = hltalt - 2369.0 * (t - tmelt)
    return hltalt


@export
def wv_saturation_calc_hltalt_codon(
    t: float,
    tmelt: float,
    ttrice: float,
    latvap: float,
    latice: float,
    pcf1: float,
    pcf2: float,
    pcf3: float,
    pcf4: float,
    pcf5: float,
    hltalt_p: cobj,
    tterm_p: cobj,
):
    hltalt_out = Ptr[float](hltalt_p)
    tterm_out = Ptr[float](tterm_p)

    hltalt = latvap
    if t >= tmelt:
        hltalt = hltalt - 2369.0 * (t - tmelt)

    tterm = 0.0
    if t < tmelt:
        tc = t - tmelt
        if tc >= -ttrice:
            weight = -tc / ttrice
            tterm = pcf5 + tc * tterm
            tterm = pcf4 + tc * tterm
            tterm = pcf3 + tc * tterm
            tterm = pcf2 + tc * tterm
            tterm = pcf1 + tc * tterm
            tterm = tterm / ttrice
        else:
            weight = 1.0
        hltalt = hltalt + weight * latice

    hltalt_out[0] = hltalt
    tterm_out[0] = tterm


@export
def wv_saturation_deriv_dqsdt_codon(
    t: float,
    p: float,
    es: float,
    qs: float,
    hltalt: float,
    tterm: float,
    rh2o: float,
    omeps: float,
) -> float:
    if qs == 1.0:
        return 0.0
    den = rh2o * t
    den = den * t
    desdt = hltalt * es
    desdt = desdt / den
    desdt = desdt + tterm
    out = qs * p
    out = out * desdt
    den2 = p - omeps * es
    den2 = es * den2
    return out / den2


@inline
def _wv_saturation_qsat_ptr_codon(
    t: float,
    p: float,
    epsilo: float,
    omeps: float,
    tmin: float,
    tmax: float,
    estbl_p: cobj,
    plenest: int,
    es_p: Ptr[float],
    qs_p: Ptr[float],
):
    es = estblf_codon(t, tmin, tmax, estbl_p, plenest)
    qs = wv_sat_svp_to_qsat_codon(es, p, epsilo, omeps)
    if p < es:
        es = p
    es_p[0] = es
    qs_p[0] = qs


@inline
def _wv_saturation_qsat_water_ptr_codon(
    t: float,
    p: float,
    idx: int,
    epsilo: float,
    omeps: float,
    es_p: Ptr[float],
    qs_p: Ptr[float],
):
    es = wv_sat_svp_water_codon(t, idx)
    qs = wv_sat_svp_to_qsat_codon(es, p, epsilo, omeps)
    if p < es:
        es = p
    es_p[0] = es
    qs_p[0] = qs


@inline
def _wv_saturation_calc_hltalt_vals_codon(
    t: float,
    tmelt: float,
    ttrice: float,
    latvap: float,
    latice: float,
    pcf1: float,
    pcf2: float,
    pcf3: float,
    pcf4: float,
    pcf5: float,
    hltalt_p: Ptr[float],
    tterm_p: Ptr[float],
):
    hltalt = latvap
    if t >= tmelt:
        hltalt = hltalt - 2369.0 * (t - tmelt)

    tterm = 0.0
    if t < tmelt:
        tc = t - tmelt
        if tc >= -ttrice:
            weight = -tc / ttrice
            tterm = pcf5 + tc * tterm
            tterm = pcf4 + tc * tterm
            tterm = pcf3 + tc * tterm
            tterm = pcf2 + tc * tterm
            tterm = pcf1 + tc * tterm
            tterm = tterm / ttrice
        else:
            weight = 1.0
        hltalt = hltalt + weight * latice

    hltalt_p[0] = hltalt
    tterm_p[0] = tterm


@inline
def _wv_saturation_no_ip_hltalt_val_codon(t: float, tmelt: float, latvap: float) -> float:
    hltalt = latvap
    if t >= tmelt:
        hltalt = hltalt - 2369.0 * (t - tmelt)
    return hltalt


@inline
def _wv_saturation_qsat_gam_enthalpy_codon(
    use_ice: int,
    t: float,
    p: float,
    idx: int,
    epsilo: float,
    omeps: float,
    cpair: float,
    tmelt: float,
    ttrice: float,
    latvap: float,
    latice: float,
    rh2o: float,
    pcf1: float,
    pcf2: float,
    pcf3: float,
    pcf4: float,
    pcf5: float,
    tmin: float,
    tmax: float,
    estbl_p: cobj,
    plenest: int,
    es_p: Ptr[float],
    qs_p: Ptr[float],
    gam_p: Ptr[float],
    enthalpy_p: Ptr[float],
):
    if use_ice != 0:
        _wv_saturation_qsat_ptr_codon(t, p, epsilo, omeps, tmin, tmax, estbl_p, plenest, es_p, qs_p)
        hltalt = 0.0
        tterm = 0.0
        _wv_saturation_calc_hltalt_vals_codon(
            t,
            tmelt,
            ttrice,
            latvap,
            latice,
            pcf1,
            pcf2,
            pcf3,
            pcf4,
            pcf5,
            __ptr__(hltalt),
            __ptr__(tterm),
        )
    else:
        _wv_saturation_qsat_water_ptr_codon(t, p, idx, epsilo, omeps, es_p, qs_p)
        hltalt = _wv_saturation_no_ip_hltalt_val_codon(t, tmelt, latvap)
        tterm = 0.0

    enthalpy_p[0] = cpair * t + hltalt * qs_p[0]
    dqsdt = wv_saturation_deriv_dqsdt_codon(t, p, es_p[0], qs_p[0], hltalt, tterm, rh2o, omeps)
    gam_p[0] = dqsdt * (hltalt / cpair)


@inline
def _wv_saturation_findsp_impl(
    q: float,
    t: float,
    p: float,
    use_ice: int,
    idx: int,
    epsilo: float,
    omeps: float,
    cpair: float,
    c3: float,
    tmelt: float,
    ttrice: float,
    latvap: float,
    latice: float,
    rh2o: float,
    pcf1: float,
    pcf2: float,
    pcf3: float,
    pcf4: float,
    pcf5: float,
    tmin: float,
    tmax: float,
    estbl_p: cobj,
    plenest: int,
    tsp_out: Ptr[float],
    qsp_out: Ptr[float],
    status_out: Ptr[i32],
):
    es = 0.0
    qs = 0.0
    if use_ice != 0:
        _wv_saturation_qsat_ptr_codon(t, p, epsilo, omeps, tmin, tmax, estbl_p, plenest, __ptr__(es), __ptr__(qs))
    else:
        _wv_saturation_qsat_water_ptr_codon(t, p, idx, epsilo, omeps, __ptr__(es), __ptr__(qs))

    if p <= 5.0 * es or qs <= 0.0 or qs >= 0.5 or t < tmin or t > tmax:
        status_out[0] = i32(1)
        tsp_out[0] = t
        qsp_out[0] = q
        return

    status = 2

    if use_ice != 0:
        hltalt = 0.0
        tterm = 0.0
        _wv_saturation_calc_hltalt_vals_codon(
            t,
            tmelt,
            ttrice,
            latvap,
            latice,
            pcf1,
            pcf2,
            pcf3,
            pcf4,
            pcf5,
            __ptr__(hltalt),
            __ptr__(tterm),
        )
    else:
        hltalt = _wv_saturation_no_ip_hltalt_val_codon(t, tmelt, latvap)

    enin = cpair * t + hltalt * q

    c1 = hltalt * c3
    c2 = (t + 36.0) ** 2
    r1b = c2 / (c2 + c1 * qs)
    qvd = r1b * (q - qs)
    tsp = t + ((hltalt / cpair) * qvd)

    gam = 0.0
    enout = 0.0
    _wv_saturation_qsat_gam_enthalpy_codon(
        use_ice,
        tsp,
        p,
        idx,
        epsilo,
        omeps,
        cpair,
        tmelt,
        ttrice,
        latvap,
        latice,
        rh2o,
        pcf1,
        pcf2,
        pcf3,
        pcf4,
        pcf5,
        tmin,
        tmax,
        estbl_p,
        plenest,
        __ptr__(es),
        qsp_out,
        __ptr__(gam),
        __ptr__(enout),
    )

    for _ in range(1, 9):
        g = enin - enout
        dgdt = -cpair * (1.0 + gam)

        t1 = tsp - g / dgdt
        dt = abs(t1 - tsp) / t1
        tsp = t1

        if tsp < tmin:
            tsp = tmin
            if use_ice != 0:
                hltalt = 0.0
                tterm = 0.0
                _wv_saturation_calc_hltalt_vals_codon(
                    tsp,
                    tmelt,
                    ttrice,
                    latvap,
                    latice,
                    pcf1,
                    pcf2,
                    pcf3,
                    pcf4,
                    pcf5,
                    __ptr__(hltalt),
                    __ptr__(tterm),
                )
            else:
                hltalt = _wv_saturation_no_ip_hltalt_val_codon(tsp, tmelt, latvap)
            qsp_out[0] = (enin - cpair * tsp) / hltalt
            enout = cpair * tsp + hltalt * qsp_out[0]
            status = 4
            break

        q1 = 0.0
        _wv_saturation_qsat_gam_enthalpy_codon(
            use_ice,
            tsp,
            p,
            idx,
            epsilo,
            omeps,
            cpair,
            tmelt,
            ttrice,
            latvap,
            latice,
            rh2o,
            pcf1,
            pcf2,
            pcf3,
            pcf4,
            pcf5,
            tmin,
            tmax,
            estbl_p,
            plenest,
            __ptr__(es),
            __ptr__(q1),
            __ptr__(gam),
            __ptr__(enout),
        )
        dq = abs(q1 - qsp_out[0]) / max(q1, 1.0e-12)
        qsp_out[0] = q1

        if dt < 1.0e-4 and dq < 1.0e-4:
            status = 0
            break

    if abs((enin - enout) / (enin + enout)) > 1.0e-4:
        status = 8

    tsp_out[0] = tsp
    status_out[0] = i32(status)


@export
def wv_saturation_findsp_codon(
    q: float,
    t: float,
    p: float,
    use_ice: int,
    idx: int,
    epsilo: float,
    omeps: float,
    cpair: float,
    c3: float,
    tmelt: float,
    ttrice: float,
    latvap: float,
    latice: float,
    rh2o: float,
    pcf1: float,
    pcf2: float,
    pcf3: float,
    pcf4: float,
    pcf5: float,
    tmin: float,
    tmax: float,
    estbl_p: cobj,
    plenest: int,
    tsp_p: cobj,
    qsp_p: cobj,
    status_p: cobj,
):
    _wv_saturation_findsp_impl(
        q,
        t,
        p,
        use_ice,
        idx,
        epsilo,
        omeps,
        cpair,
        c3,
        tmelt,
        ttrice,
        latvap,
        latice,
        rh2o,
        pcf1,
        pcf2,
        pcf3,
        pcf4,
        pcf5,
        tmin,
        tmax,
        estbl_p,
        plenest,
        Ptr[float](tsp_p),
        Ptr[float](qsp_p),
        Ptr[i32](status_p),
    )


@export
def findsp_vc_codon(
    n: int,
    use_ice: int,
    idx: int,
    epsilo: float,
    omeps: float,
    cpair: float,
    c3: float,
    tmelt: float,
    ttrice: float,
    latvap: float,
    latice: float,
    rh2o: float,
    pcf1: float,
    pcf2: float,
    pcf3: float,
    pcf4: float,
    pcf5: float,
    tmin: float,
    tmax: float,
    estbl_p: cobj,
    plenest: int,
    q_p: cobj,
    t_p: cobj,
    p_p: cobj,
    tsp_p: cobj,
    qsp_p: cobj,
    status_p: cobj,
):
    q = Ptr[float](q_p)
    t = Ptr[float](t_p)
    p = Ptr[float](p_p)
    tsp = Ptr[float](tsp_p)
    qsp = Ptr[float](qsp_p)
    status = Ptr[i32](status_p)

    for i in range(n):
        _wv_saturation_findsp_impl(
            q[i],
            t[i],
            p[i],
            use_ice,
            idx,
            epsilo,
            omeps,
            cpair,
            c3,
            tmelt,
            ttrice,
            latvap,
            latice,
            rh2o,
            pcf1,
            pcf2,
            pcf3,
            pcf4,
            pcf5,
            tmin,
            tmax,
            estbl_p,
            plenest,
            tsp + i,
            qsp + i,
            status + i,
        )


@export
def init_vdiff_codon(
    kind: int,
    expected_kind: int,
    do_iss_in: int,
    rair_in: float,
    gravit_in: float,
    rair_p: cobj,
    gravit_p: cobj,
    do_iss_p: cobj,
    status_p: cobj,
):
    rair = Ptr[float](rair_p)
    gravit = Ptr[float](gravit_p)
    do_iss = Ptr[int](do_iss_p)
    status = Ptr[int](status_p)

    if kind != expected_kind:
        status[0] = 1
        return

    rair[0] = rair_in
    gravit[0] = gravit_in
    if do_iss_in != 0:
        do_iss[0] = 1
    else:
        do_iss[0] = 0
    status[0] = 0


@export
def fin_vol_lu_decomp_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def vdiff_lu_solver_flag_codon(flag: int) -> int:
    return fin_vol_lu_decomp_codon(flag)


@export
def gw_utils_get_unit_vector_codon(u_p: cobj, v_p: cobj, u_n_p: cobj, v_n_p: cobj, mag_p: cobj, n: int):
    u = Ptr[float](u_p)
    v = Ptr[float](v_p)
    u_n = Ptr[float](u_n_p)
    v_n = Ptr[float](v_n_p)
    mag = Ptr[float](mag_p)

    for i in range(n):
        mag[i] = sqrt(u[i] * u[i] + v[i] * v[i])

    for i in range(n):
        if mag[i] > 0.0:
            u_n[i] = u[i] / mag[i]
            v_n[i] = v[i] / mag[i]
        else:
            u_n[i] = 0.0
            v_n[i] = 0.0


@export
def gw_utils_dot_2d_codon(u1_p: cobj, v1_p: cobj, u2_p: cobj, v2_p: cobj, out_p: cobj, n: int):
    u1 = Ptr[float](u1_p)
    v1 = Ptr[float](v1_p)
    u2 = Ptr[float](u2_p)
    v2 = Ptr[float](v2_p)
    out = Ptr[float](out_p)

    for i in range(n):
        out[i] = u1[i] * u2[i] + v1[i] * v2[i]


@export
def gw_utils_midpoint_interp_codon(arr_p: cobj, interp_p: cobj, n1: int, n2: int):
    arr = Ptr[float](arr_p)
    interp = Ptr[float](interp_p)

    for k in range(n2 - 1):
        for i in range(n1):
            interp[i + k * n1] = 0.5 * (arr[i + k * n1] + arr[i + (k + 1) * n1])


@inline
def _cmparray_idx3(i: int, j: int, k: int, il1: int, il2: int, il3: int, ld1: int, ld2: int) -> int:
    """cmparray arrays declared as (il1:iu1, il2:iu2, il3:iu3)."""
    return (i - il1) + (j - il2) * ld1 + (k - il3) * ld1 * ld2


@export
def cmparray_daynite_copy_real_codon(
    in_array_p: cobj,
    out_array_p: cobj,
    idxday_p: cobj,
    idxnite_p: cobj,
    nday: int,
    nnite: int,
    il1: int,
    iu1: int,
    il2: int,
    iu2: int,
    il3: int,
    iu3: int,
):
    in_array = Ptr[float](in_array_p)
    out_array = Ptr[float](out_array_p)
    idxday = Ptr[i32](idxday_p)
    idxnite = Ptr[i32](idxnite_p)
    ld1 = iu1 - il1 + 1
    ld2 = iu2 - il2 + 1

    for k in range(il3, iu3 + 1):
        for j in range(il2, iu2 + 1):
            for pos in range(nday):
                i_out = il1 + pos
                i_src = int(idxday[pos])
                out_array[_cmparray_idx3(i_out, j, k, il1, il2, il3, ld1, ld2)] = in_array[
                    _cmparray_idx3(i_src, j, k, il1, il2, il3, ld1, ld2)
                ]
            for pos in range(nnite):
                i_out = il1 + nday + pos
                i_src = int(idxnite[pos])
                out_array[_cmparray_idx3(i_out, j, k, il1, il2, il3, ld1, ld2)] = in_array[
                    _cmparray_idx3(i_src, j, k, il1, il2, il3, ld1, ld2)
                ]


@export
def cmpdaynite_1d_r_copy_codon(
    in_array_p: cobj,
    out_array_p: cobj,
    idxday_p: cobj,
    idxnite_p: cobj,
    nday: int,
    nnite: int,
    il1: int,
    iu1: int,
    il2: int,
    iu2: int,
    il3: int,
    iu3: int,
):
    cmparray_daynite_copy_real_codon(
        in_array_p, out_array_p, idxday_p, idxnite_p, nday, nnite, il1, iu1, il2, iu2, il3, iu3
    )


@export
def cmpdaynite_2d_r_copy_codon(
    in_array_p: cobj,
    out_array_p: cobj,
    idxday_p: cobj,
    idxnite_p: cobj,
    nday: int,
    nnite: int,
    il1: int,
    iu1: int,
    il2: int,
    iu2: int,
    il3: int,
    iu3: int,
):
    cmparray_daynite_copy_real_codon(
        in_array_p, out_array_p, idxday_p, idxnite_p, nday, nnite, il1, iu1, il2, iu2, il3, iu3
    )


@export
def cmparray_exp_daynite_real_codon(
    array_p: cobj,
    tmp_p: cobj,
    idxday_p: cobj,
    idxnite_p: cobj,
    nday: int,
    nnite: int,
    il1: int,
    iu1: int,
    il2: int,
    iu2: int,
    il3: int,
    iu3: int,
):
    array = Ptr[float](array_p)
    tmp = Ptr[float](tmp_p)
    idxday = Ptr[i32](idxday_p)
    idxnite = Ptr[i32](idxnite_p)
    ld1 = iu1 - il1 + 1
    ld2 = iu2 - il2 + 1

    for k in range(il3, iu3 + 1):
        for j in range(il2, iu2 + 1):
            for pos in range(nday):
                i_src = il1 + pos
                tmp[pos] = array[_cmparray_idx3(i_src, j, k, il1, il2, il3, ld1, ld2)]
            for pos in range(nnite):
                i_dst = int(idxnite[pos])
                i_src = il1 + nday + pos
                array[_cmparray_idx3(i_dst, j, k, il1, il2, il3, ld1, ld2)] = array[
                    _cmparray_idx3(i_src, j, k, il1, il2, il3, ld1, ld2)
                ]
            for pos in range(nday):
                i_dst = int(idxday[pos])
                array[_cmparray_idx3(i_dst, j, k, il1, il2, il3, ld1, ld2)] = tmp[pos]


@export
def expdaynite_1d_r_codon(
    array_p: cobj,
    tmp_p: cobj,
    idxday_p: cobj,
    idxnite_p: cobj,
    nday: int,
    nnite: int,
    il1: int,
    iu1: int,
    il2: int,
    iu2: int,
    il3: int,
    iu3: int,
):
    cmparray_exp_daynite_real_codon(
        array_p, tmp_p, idxday_p, idxnite_p, nday, nnite, il1, iu1, il2, iu2, il3, iu3
    )


@export
def expdaynite_2d_r_codon(
    array_p: cobj,
    tmp_p: cobj,
    idxday_p: cobj,
    idxnite_p: cobj,
    nday: int,
    nnite: int,
    il1: int,
    iu1: int,
    il2: int,
    iu2: int,
    il3: int,
    iu3: int,
):
    cmparray_exp_daynite_real_codon(
        array_p, tmp_p, idxday_p, idxnite_p, nday, nnite, il1, iu1, il2, iu2, il3, iu3
    )
