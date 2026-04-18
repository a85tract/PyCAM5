from math import copysign


@inline
def _hy_idx(klev: int) -> int:
    """hvcoord%hyai / hvcoord%hybi declared as (nlev+1)."""
    return klev - 1


@inline
def _plane_idx(iidx: int, jidx: int, np: int) -> int:
    """ps_v and dp_np1 declared as (np,np)."""
    return (iidx - 1) + (jidx - 1) * np


@inline
def _q_idx(iidx: int, jidx: int, klev: int, qidx: int, np: int, nlev: int) -> int:
    """q and qdp slices declared as (np,np,nlev,qsize)."""
    return (
        (iidx - 1)
        + (jidx - 1) * np
        + (klev - 1) * np * np
        + (qidx - 1) * np * np * nlev
    )


@inline
def _q_tl_idx(iidx: int, jidx: int, klev: int, qidx: int, tlidx: int, np: int, nlev: int, qsize: int) -> int:
    """state%Qdp declared as (np,np,nlev,qsize,2)."""
    return _q_idx(iidx, jidx, klev, qidx, np, nlev) + (tlidx - 1) * np * np * nlev * qsize


@inline
def _vol_idx(iidx: int, jidx: int, klev: int, np: int) -> int:
    """dp3d, dp, dp_star declared as (np,np,nlev)."""
    return (iidx - 1) + (jidx - 1) * np + (klev - 1) * np * np


@inline
def _field_vol_idx(iidx: int, jidx: int, klev: int, fidx: int, np: int, nlev: int) -> int:
    """ttmp declared as (np,np,nlev,2)."""
    return _vol_idx(iidx, jidx, klev, np) + (fidx - 1) * np * np * nlev


@inline
def _v_idx(iidx: int, jidx: int, comp: int, klev: int, np: int) -> int:
    """state%v slice declared as (np,np,2,nlev)."""
    return (iidx - 1) + (jidx - 1) * np + (comp - 1) * np * np + (klev - 1) * np * np * 2


@inline
def _ghost_col_idx(klev: int) -> int:
    """dpo declared as (-1:nlev+2)."""
    return klev + 1


@inline
def _ppm_grid_idx(row: int, jidx: int) -> int:
    """ppmdx declared as (10,0:nlev+1)."""
    return (row - 1) + jidx * 10


@inline
def _ppm_scratch_idx(jidx: int) -> int:
    """ppm_ai declared as (0:nlev), ppm_dma declared as (0:nlev+1)."""
    return jidx


@inline
def _ppm_coef_idx(comp: int, jidx: int) -> int:
    """coefs declared as (0:2,nlev)."""
    return comp + (jidx - 1) * 3


@inline
def _col_idx(klev: int) -> int:
    """pio declared as (nlev+2), pin declared as (nlev+1), z1/z2/kid declared as (nlev)."""
    return klev - 1


@export
def prim_subcycle_dp3d_init_codon(
    np: int,
    nlev: int,
    ps0: float,
    hyai_p: cobj,
    hybi_p: cobj,
    ps_v_p: cobj,
    dp3d_p: cobj,
):
    hyai = Ptr[float](hyai_p)
    hybi = Ptr[float](hybi_p)
    ps_v = Ptr[float](ps_v_p)
    dp3d = Ptr[float](dp3d_p)

    for k in range(1, nlev + 1):
        hyai_delta = hyai[_hy_idx(k + 1)] - hyai[_hy_idx(k)]
        hybi_delta = hybi[_hy_idx(k + 1)] - hybi[_hy_idx(k)]

        for j in range(1, np + 1):
            for i in range(1, np + 1):
                plane_idx = _plane_idx(i, j, np)
                vol_idx = _vol_idx(i, j, k, np)
                dp3d[vol_idx] = hyai_delta * ps0 + hybi_delta * ps_v[plane_idx]


@export
def prim_subcycle_q_update_codon(
    np: int,
    nlev: int,
    qsize: int,
    ps0: float,
    hyai_p: cobj,
    hybi_p: cobj,
    ps_v_p: cobj,
    dp_np1_p: cobj,
    qdp_p: cobj,
    q_p: cobj,
):
    hyai = Ptr[float](hyai_p)
    hybi = Ptr[float](hybi_p)
    ps_v = Ptr[float](ps_v_p)
    dp_np1 = Ptr[float](dp_np1_p)
    qdp = Ptr[float](qdp_p)
    q_out = Ptr[float](q_p)

    for k in range(1, nlev + 1):
        hyai_delta = hyai[_hy_idx(k + 1)] - hyai[_hy_idx(k)]
        hybi_delta = hybi[_hy_idx(k + 1)] - hybi[_hy_idx(k)]

        for j in range(1, np + 1):
            for i in range(1, np + 1):
                plane_idx = _plane_idx(i, j, np)
                dp_np1[plane_idx] = hyai_delta * ps0 + hybi_delta * ps_v[plane_idx]

        for qidx in range(1, qsize + 1):
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    plane_idx = _plane_idx(i, j, np)
                    q_idx = _q_idx(i, j, k, qidx, np, nlev)
                    q_out[q_idx] = qdp[q_idx] / dp_np1[plane_idx]


@export
def qdp_time_avg_codon(
    np: int,
    nlev: int,
    qsize: int,
    rkstage: int,
    n0_qdp: int,
    np1_qdp: int,
    qdp_p: cobj,
):
    qdp = Ptr[float](qdp_p)
    rkstage_minus_one = rkstage - 1

    for qidx in range(1, qsize + 1):
        for k in range(1, nlev + 1):
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    n0_idx = _q_tl_idx(i, j, k, qidx, n0_qdp, np, nlev, qsize)
                    np1_idx = _q_tl_idx(i, j, k, qidx, np1_qdp, np, nlev, qsize)
                    qdp[np1_idx] = (qdp[n0_idx] + rkstage_minus_one * qdp[np1_idx]) / rkstage


@export
def euler_step_vstar_prepare_codon(
    np: int,
    nlev: int,
    dt: float,
    rhs_multiplier: int,
    dp_in_p: cobj,
    divdp_proj_p: cobj,
    vn0_p: cobj,
    dp_out_p: cobj,
    vstar_p: cobj,
):
    dp_in = Ptr[float](dp_in_p)
    divdp_proj = Ptr[float](divdp_proj_p)
    vn0 = Ptr[float](vn0_p)
    dp_out = Ptr[float](dp_out_p)
    vstar = Ptr[float](vstar_p)

    for k in range(1, nlev + 1):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                vol_idx = _vol_idx(i, j, k, np)
                dp_val = dp_in[vol_idx] - rhs_multiplier * dt * divdp_proj[vol_idx]
                dp_out[vol_idx] = dp_val
                vstar[_v_idx(i, j, 1, k, np)] = vn0[_v_idx(i, j, 1, k, np)] / dp_val
                vstar[_v_idx(i, j, 2, k, np)] = vn0[_v_idx(i, j, 2, k, np)] / dp_val


@export
def euler_step_limiter_dpstar_codon(
    np: int,
    nlev: int,
    dt: float,
    rhs_viss: int,
    nu_q: float,
    nu_p: float,
    dp_in_p: cobj,
    divdp_p: cobj,
    dpdiss_biharmonic_p: cobj,
    spheremp_p: cobj,
    dp_star_p: cobj,
):
    dp_in = Ptr[float](dp_in_p)
    divdp = Ptr[float](divdp_p)
    dpdiss_biharmonic = Ptr[float](dpdiss_biharmonic_p)
    spheremp = Ptr[float](spheremp_p)
    dp_star = Ptr[float](dp_star_p)

    use_dpdiss = nu_p > 0.0 and rhs_viss != 0

    for k in range(1, nlev + 1):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                vol_idx = _vol_idx(i, j, k, np)
                plane_idx = _plane_idx(i, j, np)
                dp_val = dp_in[vol_idx] - dt * divdp[vol_idx]
                if use_dpdiss:
                    dp_val = dp_val - rhs_viss * dt * nu_q * dpdiss_biharmonic[vol_idx] / spheremp[plane_idx]
                dp_star[vol_idx] = dp_val


@export
def euler_step_qdp_writeback_codon(
    np: int,
    nlev: int,
    qsize: int,
    qidx: int,
    np1_qdp: int,
    qdp_p: cobj,
    spheremp_p: cobj,
    qtens_p: cobj,
):
    qdp = Ptr[float](qdp_p)
    spheremp = Ptr[float](spheremp_p)
    qtens = Ptr[float](qtens_p)

    for k in range(1, nlev + 1):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                plane_idx = _plane_idx(i, j, np)
                vol_idx = _vol_idx(i, j, k, np)
                qdp[_q_tl_idx(i, j, k, qidx, np1_qdp, np, nlev, qsize)] = spheremp[plane_idx] * qtens[vol_idx]


@export
def euler_step_qdp_restore_codon(
    np: int,
    nlev: int,
    qdp_p: cobj,
    rspheremp_p: cobj,
):
    qdp = Ptr[float](qdp_p)
    rspheremp = Ptr[float](rspheremp_p)

    for k in range(1, nlev + 1):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                plane_idx = _plane_idx(i, j, np)
                vol_idx = _vol_idx(i, j, k, np)
                qdp[vol_idx] = rspheremp[plane_idx] * qdp[vol_idx]


@export
def advance_hypervis_qtens_prepare_codon(
    np: int,
    nlev: int,
    qsize: int,
    ps0: float,
    dt2: float,
    nu_p: float,
    hyai_p: cobj,
    hybi_p: cobj,
    dp_in_p: cobj,
    divdp_proj_p: cobj,
    dpdiss_ave_p: cobj,
    qdp_p: cobj,
    dp_out_p: cobj,
    qtens_p: cobj,
):
    hyai = Ptr[float](hyai_p)
    hybi = Ptr[float](hybi_p)
    dp_in = Ptr[float](dp_in_p)
    divdp_proj = Ptr[float](divdp_proj_p)
    dpdiss_ave = Ptr[float](dpdiss_ave_p)
    qdp = Ptr[float](qdp_p)
    dp_out = Ptr[float](dp_out_p)
    qtens = Ptr[float](qtens_p)

    use_dpdiss = nu_p > 0.0

    for k in range(1, nlev + 1):
        dp0 = (hyai[_hy_idx(k + 1)] - hyai[_hy_idx(k)]) * ps0 + (hybi[_hy_idx(k + 1)] - hybi[_hy_idx(k)]) * ps0
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                vol_idx = _vol_idx(i, j, k, np)
                dp_val = dp_in[vol_idx] - dt2 * divdp_proj[vol_idx]
                dp_out[vol_idx] = dp_val
                if use_dpdiss:
                    for q in range(1, qsize + 1):
                        q_idx = _q_idx(i, j, k, q, np, nlev)
                        qtens[q_idx] = dpdiss_ave[vol_idx] * qdp[q_idx] / dp_val
                else:
                    for q in range(1, qsize + 1):
                        q_idx = _q_idx(i, j, k, q, np, nlev)
                        qtens[q_idx] = dp0 * qdp[q_idx] / dp_val


@export
def advance_hypervis_qdp_update_codon(
    np: int,
    nlev: int,
    dt: float,
    nu_q: float,
    qdp_p: cobj,
    spheremp_p: cobj,
    qtens_p: cobj,
):
    qdp = Ptr[float](qdp_p)
    spheremp = Ptr[float](spheremp_p)
    qtens = Ptr[float](qtens_p)

    for k in range(1, nlev + 1):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                plane_idx = _plane_idx(i, j, np)
                vol_idx = _vol_idx(i, j, k, np)
                qdp[vol_idx] = qdp[vol_idx] * spheremp[plane_idx] - dt * nu_q * qtens[vol_idx]


@export
def euler_step_dssvar_restore_codon(
    np: int,
    nlev: int,
    dssvar_p: cobj,
    rspheremp_p: cobj,
):
    dssvar = Ptr[float](dssvar_p)
    rspheremp = Ptr[float](rspheremp_p)

    for k in range(1, nlev + 1):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                plane_idx = _plane_idx(i, j, np)
                vol_idx = _vol_idx(i, j, k, np)
                dssvar[vol_idx] = dssvar[vol_idx] * rspheremp[plane_idx]


@export
def euler_step_dssvar_pack_codon(
    np: int,
    nlev: int,
    dssvar_p: cobj,
    spheremp_p: cobj,
):
    dssvar = Ptr[float](dssvar_p)
    spheremp = Ptr[float](spheremp_p)

    for k in range(1, nlev + 1):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                plane_idx = _plane_idx(i, j, np)
                vol_idx = _vol_idx(i, j, k, np)
                dssvar[vol_idx] = spheremp[plane_idx] * dssvar[vol_idx]


@export
def euler_step_qtens_base_codon(
    np: int,
    dt: float,
    qdp_p: cobj,
    dp_star_p: cobj,
    qtens_p: cobj,
):
    qdp = Ptr[float](qdp_p)
    dp_star = Ptr[float](dp_star_p)
    qtens = Ptr[float](qtens_p)

    for j in range(1, np + 1):
        for i in range(1, np + 1):
            plane_idx = _plane_idx(i, j, np)
            qtens[plane_idx] = qdp[plane_idx] - dt * dp_star[plane_idx]


@export
def euler_step_gradq_prepare_codon(
    np: int,
    vstar1_p: cobj,
    vstar2_p: cobj,
    qdp_p: cobj,
    gradq1_p: cobj,
    gradq2_p: cobj,
):
    vstar1 = Ptr[float](vstar1_p)
    vstar2 = Ptr[float](vstar2_p)
    qdp = Ptr[float](qdp_p)
    gradq1 = Ptr[float](gradq1_p)
    gradq2 = Ptr[float](gradq2_p)

    for j in range(1, np + 1):
        for i in range(1, np + 1):
            plane_idx = _plane_idx(i, j, np)
            gradq1[plane_idx] = vstar1[plane_idx] * qdp[plane_idx]
            gradq2[plane_idx] = vstar2[plane_idx] * qdp[plane_idx]


@export
def euler_step_qtens_biharmonic_add_codon(
    np: int,
    qtens_p: cobj,
    qtens_biharmonic_p: cobj,
):
    qtens = Ptr[float](qtens_p)
    qtens_biharmonic = Ptr[float](qtens_biharmonic_p)

    for j in range(1, np + 1):
        for i in range(1, np + 1):
            plane_idx = _plane_idx(i, j, np)
            qtens[plane_idx] = qtens[plane_idx] + qtens_biharmonic[plane_idx]


@export
def euler_step_qtens_biharmonic_init_codon(
    np: int,
    nlev: int,
    qsize: int,
    dt: float,
    rhs_multiplier: int,
    dp_in_p: cobj,
    divdp_proj_p: cobj,
    qdp_p: cobj,
    dp_out_p: cobj,
    qtens_biharmonic_p: cobj,
):
    dp_in = Ptr[float](dp_in_p)
    divdp_proj = Ptr[float](divdp_proj_p)
    qdp = Ptr[float](qdp_p)
    dp_out = Ptr[float](dp_out_p)
    qtens_biharmonic = Ptr[float](qtens_biharmonic_p)

    for k in range(1, nlev + 1):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                vol_idx = _vol_idx(i, j, k, np)
                dp_val = dp_in[vol_idx] - rhs_multiplier * dt * divdp_proj[vol_idx]
                dp_out[vol_idx] = dp_val
                for q in range(1, qsize + 1):
                    q_idx = _q_idx(i, j, k, q, np, nlev)
                    qtens_biharmonic[q_idx] = qdp[q_idx] / dp_val


@export
def euler_step_qtens_biharmonic_scale_codon(
    np: int,
    nlev: int,
    qsize: int,
    qtens_biharmonic_p: cobj,
    dpdiss_ave_p: cobj,
    dp0_p: cobj,
):
    qtens_biharmonic = Ptr[float](qtens_biharmonic_p)
    dpdiss_ave = Ptr[float](dpdiss_ave_p)
    dp0 = Ptr[float](dp0_p)

    for q in range(1, qsize + 1):
        for k in range(1, nlev + 1):
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    q_idx = _q_idx(i, j, k, q, np, nlev)
                    vol_idx = _vol_idx(i, j, k, np)
                    qtens_biharmonic[q_idx] = qtens_biharmonic[q_idx] * dpdiss_ave[vol_idx] / dp0[_col_idx(k)]


@export
def euler_step_qtens_biharmonic_unapply_codon(
    np: int,
    nlev: int,
    rhs_viss: int,
    dt: float,
    nu_q: float,
    qtens_biharmonic_p: cobj,
    spheremp_p: cobj,
    dp0_p: cobj,
):
    qtens_biharmonic = Ptr[float](qtens_biharmonic_p)
    spheremp = Ptr[float](spheremp_p)
    dp0 = Ptr[float](dp0_p)

    for k in range(1, nlev + 1):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                vol_idx = _vol_idx(i, j, k, np)
                plane_idx = _plane_idx(i, j, np)
                qtens_biharmonic[vol_idx] = -rhs_viss * dt * nu_q * dp0[_col_idx(k)] * qtens_biharmonic[vol_idx] / spheremp[plane_idx]


@export
def vertical_remap_rsplit_prepare_codon(
    np: int,
    nlev: int,
    ps0: float,
    hyai_p: cobj,
    hybi_p: cobj,
    ps_v_p: cobj,
    dp3d_p: cobj,
    dp_p: cobj,
    dp_star_p: cobj,
):
    hyai = Ptr[float](hyai_p)
    hybi = Ptr[float](hybi_p)
    ps_v = Ptr[float](ps_v_p)
    dp3d = Ptr[float](dp3d_p)
    dp = Ptr[float](dp_p)
    dp_star = Ptr[float](dp_star_p)

    for k in range(1, nlev + 1):
        hyai_delta = hyai[_hy_idx(k + 1)] - hyai[_hy_idx(k)]
        hybi_delta = hybi[_hy_idx(k + 1)] - hybi[_hy_idx(k)]

        for j in range(1, np + 1):
            for i in range(1, np + 1):
                plane_idx = _plane_idx(i, j, np)
                vol_idx = _vol_idx(i, j, k, np)
                dp[vol_idx] = hyai_delta * ps0 + hybi_delta * ps_v[plane_idx]
                dp_star[vol_idx] = dp3d[vol_idx]


@export
def vertical_remap_t_scale_codon(
    np: int,
    nlev: int,
    t_p: cobj,
    dp_star_p: cobj,
    ttmp_p: cobj,
):
    t = Ptr[float](t_p)
    dp_star = Ptr[float](dp_star_p)
    ttmp = Ptr[float](ttmp_p)

    for k in range(1, nlev + 1):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                vol_idx = _vol_idx(i, j, k, np)
                ttmp[vol_idx] = t[vol_idx] * dp_star[vol_idx]


@export
def vertical_remap_t_unscale_codon(
    np: int,
    nlev: int,
    ttmp_p: cobj,
    dp_p: cobj,
    t_p: cobj,
):
    ttmp = Ptr[float](ttmp_p)
    dp = Ptr[float](dp_p)
    t = Ptr[float](t_p)

    for k in range(1, nlev + 1):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                vol_idx = _vol_idx(i, j, k, np)
                t[vol_idx] = ttmp[vol_idx] / dp[vol_idx]


@export
def vertical_remap_v_scale_codon(
    np: int,
    nlev: int,
    v_p: cobj,
    dp_star_p: cobj,
    ttmp_p: cobj,
):
    v = Ptr[float](v_p)
    dp_star = Ptr[float](dp_star_p)
    ttmp = Ptr[float](ttmp_p)

    for k in range(1, nlev + 1):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                vol_idx = _vol_idx(i, j, k, np)
                ttmp[_field_vol_idx(i, j, k, 1, np, nlev)] = v[_v_idx(i, j, 1, k, np)] * dp_star[vol_idx]
                ttmp[_field_vol_idx(i, j, k, 2, np, nlev)] = v[_v_idx(i, j, 2, k, np)] * dp_star[vol_idx]


@export
def vertical_remap_v_unscale_codon(
    np: int,
    nlev: int,
    ttmp_p: cobj,
    dp_p: cobj,
    v_p: cobj,
):
    ttmp = Ptr[float](ttmp_p)
    dp = Ptr[float](dp_p)
    v = Ptr[float](v_p)

    for k in range(1, nlev + 1):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                vol_idx = _vol_idx(i, j, k, np)
                v[_v_idx(i, j, 1, k, np)] = ttmp[_field_vol_idx(i, j, k, 1, np, nlev)] / dp[vol_idx]
                v[_v_idx(i, j, 2, k, np)] = ttmp[_field_vol_idx(i, j, k, 2, np, nlev)] / dp[vol_idx]


@export
def vertical_remap_ps_v_update_codon(
    np: int,
    nlev: int,
    hyai1: float,
    ps0: float,
    dp3d_p: cobj,
    ps_v_p: cobj,
):
    dp3d = Ptr[float](dp3d_p)
    ps_v = Ptr[float](ps_v_p)

    for j in range(1, np + 1):
        for i in range(1, np + 1):
            total = 0.0
            for k in range(1, nlev + 1):
                total += dp3d[_vol_idx(i, j, k, np)]
            ps_v[_plane_idx(i, j, np)] = hyai1 * ps0 + total


@export
def remap_q_ppm_interval_setup_codon(
    nlev: int,
    pio_p: cobj,
    pin_p: cobj,
    dpo_p: cobj,
    kid_p: cobj,
    z1_p: cobj,
    z2_p: cobj,
):
    pio = Ptr[float](pio_p)
    pin = Ptr[float](pin_p)
    dpo = Ptr[float](dpo_p)
    kid = Ptr[int](kid_p)
    z1 = Ptr[float](z1_p)
    z2 = Ptr[float](z2_p)

    pio[_col_idx(nlev + 2)] = pio[_col_idx(nlev + 1)] + 1.0
    pin[_col_idx(nlev + 1)] = pio[_col_idx(nlev + 1)]

    for k in range(1, 3):
        dpo[_ghost_col_idx(1 - k)] = dpo[_ghost_col_idx(k)]
        dpo[_ghost_col_idx(nlev + k)] = dpo[_ghost_col_idx(nlev + 1 - k)]

    for k in range(1, nlev + 1):
        kk = k
        while pio[_col_idx(kk)] <= pin[_col_idx(k + 1)]:
            kk += 1
        kk -= 1
        if kk == nlev + 1:
            kk = nlev
        kid[_col_idx(k)] = kk
        z1[_col_idx(k)] = -0.5
        z2[_col_idx(k)] = (pin[_col_idx(k + 1)] - (pio[_col_idx(kk)] + pio[_col_idx(kk + 1)]) * 0.5) / dpo[
            _ghost_col_idx(kk)
        ]


@export
def remap_q_ppm_mass_prep_codon(
    nx: int,
    nlev: int,
    qsize: int,
    iidx: int,
    jidx: int,
    qidx: int,
    qdp_p: cobj,
    dpo_p: cobj,
    masso_p: cobj,
    ao_p: cobj,
):
    qdp = Ptr[float](qdp_p)
    dpo = Ptr[float](dpo_p)
    masso = Ptr[float](masso_p)
    ao = Ptr[float](ao_p)

    masso[_col_idx(1)] = 0.0
    for k in range(1, nlev + 1):
        ao[_ghost_col_idx(k)] = qdp[_q_idx(iidx, jidx, k, qidx, nx, nlev)]
        masso[_col_idx(k + 1)] = masso[_col_idx(k)] + ao[_ghost_col_idx(k)]
        ao[_ghost_col_idx(k)] = ao[_ghost_col_idx(k)] / dpo[_ghost_col_idx(k)]

    for k in range(1, 3):
        ao[_ghost_col_idx(1 - k)] = ao[_ghost_col_idx(k)]
        ao[_ghost_col_idx(nlev + k)] = ao[_ghost_col_idx(nlev + 1 - k)]


@export
def remap_q_ppm_compute_ppm_grids_codon(
    nlev: int,
    vert_remap_q_alg: int,
    dx_p: cobj,
    rslt_p: cobj,
):
    dx = Ptr[float](dx_p)
    rslt = Ptr[float](rslt_p)

    if vert_remap_q_alg == 2:
        indB = 2
        indE = nlev - 1
    else:
        indB = 0
        indE = nlev + 1

    for j in range(indB, indE + 1):
        dxjm1 = dx[_ghost_col_idx(j - 1)]
        dxj = dx[_ghost_col_idx(j)]
        dxjp1 = dx[_ghost_col_idx(j + 1)]
        rslt[_ppm_grid_idx(1, j)] = dxj / ((dxjm1 + dxj) + dxjp1)
        rslt[_ppm_grid_idx(2, j)] = ((2.0 * dxjm1) + dxj) / (dxjp1 + dxj)
        rslt[_ppm_grid_idx(3, j)] = (dxj + (2.0 * dxjp1)) / (dxjm1 + dxj)

    if vert_remap_q_alg == 2:
        indB = 2
        indE = nlev - 2
    else:
        indB = 0
        indE = nlev

    for j in range(indB, indE + 1):
        dxjm1 = dx[_ghost_col_idx(j - 1)]
        dxj = dx[_ghost_col_idx(j)]
        dxjp1 = dx[_ghost_col_idx(j + 1)]
        dxjp2 = dx[_ghost_col_idx(j + 2)]
        rslt[_ppm_grid_idx(4, j)] = dxj / (dxj + dxjp1)
        rslt[_ppm_grid_idx(5, j)] = 1.0 / (((dxjm1 + dxj) + dxjp1) + dxjp2)
        rslt[_ppm_grid_idx(6, j)] = ((2.0 * dxjp1) * dxj) / (dxj + dxjp1)
        rslt[_ppm_grid_idx(7, j)] = (dxjm1 + dxj) / ((2.0 * dxj) + dxjp1)
        rslt[_ppm_grid_idx(8, j)] = (dxjp2 + dxjp1) / ((2.0 * dxjp1) + dxj)
        rslt[_ppm_grid_idx(9, j)] = dxj * (dxjm1 + dxj) / ((2.0 * dxj) + dxjp1)
        rslt[_ppm_grid_idx(10, j)] = dxjp1 * (dxjp1 + dxjp2) / (dxj + (2.0 * dxjp1))


@export
def remap_q_ppm_compute_ppm_codon(
    nlev: int,
    vert_remap_q_alg: int,
    a_p: cobj,
    dx_p: cobj,
    ai_p: cobj,
    dma_p: cobj,
    coefs_p: cobj,
):
    a = Ptr[float](a_p)
    dx = Ptr[float](dx_p)
    ai = Ptr[float](ai_p)
    dma = Ptr[float](dma_p)
    coefs = Ptr[float](coefs_p)

    if vert_remap_q_alg == 2:
        indB = 2
        indE = nlev - 1
    else:
        indB = 0
        indE = nlev + 1

    for j in range(indB, indE + 1):
        ajm1 = a[_ghost_col_idx(j - 1)]
        aj = a[_ghost_col_idx(j)]
        ajp1 = a[_ghost_col_idx(j + 1)]
        da = dx[_ppm_grid_idx(1, j)] * (
            dx[_ppm_grid_idx(2, j)] * (ajp1 - aj) + dx[_ppm_grid_idx(3, j)] * (aj - ajm1)
        )
        dma_val = min(
            abs(da),
            min(2.0 * abs(aj - ajm1), 2.0 * abs(ajp1 - aj)),
        ) * copysign(1.0, da)
        if (ajp1 - aj) * (aj - ajm1) <= 0.0:
            dma_val = 0.0
        dma[_ppm_scratch_idx(j)] = dma_val

    if vert_remap_q_alg == 2:
        indB = 2
        indE = nlev - 2
    else:
        indB = 0
        indE = nlev

    for j in range(indB, indE + 1):
        aj = a[_ghost_col_idx(j)]
        ajp1 = a[_ghost_col_idx(j + 1)]
        ai[_ppm_scratch_idx(j)] = aj + dx[_ppm_grid_idx(4, j)] * (ajp1 - aj) + dx[_ppm_grid_idx(5, j)] * (
            dx[_ppm_grid_idx(6, j)] * (dx[_ppm_grid_idx(7, j)] - dx[_ppm_grid_idx(8, j)]) * (ajp1 - aj)
            - dx[_ppm_grid_idx(9, j)] * dma[_ppm_scratch_idx(j + 1)]
            + dx[_ppm_grid_idx(10, j)] * dma[_ppm_scratch_idx(j)]
        )

    if vert_remap_q_alg == 2:
        indB = 3
        indE = nlev - 2
    else:
        indB = 1
        indE = nlev

    for j in range(indB, indE + 1):
        aj = a[_ghost_col_idx(j)]
        al = ai[_ppm_scratch_idx(j - 1)]
        ar = ai[_ppm_scratch_idx(j)]
        if (ar - aj) * (aj - al) <= 0.0:
            al = aj
            ar = aj
        if (ar - al) * (aj - (al + ar) / 2.0) > ((ar - al) * (ar - al)) / 6.0:
            al = 3.0 * aj - 2.0 * ar
        if (ar - al) * (aj - (al + ar) / 2.0) < -(((ar - al) * (ar - al)) / 6.0):
            ar = 3.0 * aj - 2.0 * al
        coefs[_ppm_coef_idx(0, j)] = 1.5 * aj - (al + ar) / 4.0
        coefs[_ppm_coef_idx(1, j)] = ar - al
        coefs[_ppm_coef_idx(2, j)] = -6.0 * aj + 3.0 * (al + ar)

    if vert_remap_q_alg == 2:
        for j in range(1, 3):
            aj = a[_ghost_col_idx(j)]
            coefs[_ppm_coef_idx(0, j)] = aj
            coefs[_ppm_coef_idx(1, j)] = 0.0
            coefs[_ppm_coef_idx(2, j)] = 0.0
        for j in range(nlev - 1, nlev + 1):
            aj = a[_ghost_col_idx(j)]
            coefs[_ppm_coef_idx(0, j)] = aj
            coefs[_ppm_coef_idx(1, j)] = 0.0
            coefs[_ppm_coef_idx(2, j)] = 0.0


@export
def remap_q_ppm_mass_apply_codon(
    nx: int,
    nlev: int,
    iidx: int,
    jidx: int,
    qidx: int,
    kid_p: cobj,
    masso_p: cobj,
    coefs_p: cobj,
    z1_p: cobj,
    z2_p: cobj,
    dpo_p: cobj,
    qdp_p: cobj,
):
    kid = Ptr[int](kid_p)
    masso = Ptr[float](masso_p)
    coefs = Ptr[float](coefs_p)
    z1 = Ptr[float](z1_p)
    z2 = Ptr[float](z2_p)
    dpo = Ptr[float](dpo_p)
    qdp = Ptr[float](qdp_p)

    massn1 = 0.0
    for k in range(1, nlev + 1):
        kk = kid[_col_idx(k)]
        x1 = z1[_col_idx(k)]
        x2 = z2[_col_idx(k)]
        x1_sq = x1 * x1
        x2_sq = x2 * x2
        mass = (
            coefs[_ppm_coef_idx(0, kk)] * (x2 - x1)
            + coefs[_ppm_coef_idx(1, kk)] * (x2_sq - x1_sq) / 2.0
            + coefs[_ppm_coef_idx(2, kk)] * ((x2_sq * x2) - (x1_sq * x1)) / 3.0
        )
        massn2 = masso[_col_idx(kk)] + mass * dpo[_ghost_col_idx(kk)]
        qdp[_q_idx(iidx, jidx, k, qidx, nx, nlev)] = massn2 - massn1
        massn1 = massn2
