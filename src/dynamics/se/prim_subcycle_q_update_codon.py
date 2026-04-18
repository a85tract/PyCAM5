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
def _vec2_idx(iidx: int, jidx: int, comp: int, np: int) -> int:
    """v and gv declared as (np,np,2)."""
    return (iidx - 1) + (jidx - 1) * np + (comp - 1) * np * np


@inline
def _vec3_idx(iidx: int, jidx: int, comp: int, np: int) -> int:
    """dum_cart declared as (np,np,3)."""
    return (iidx - 1) + (jidx - 1) * np + (comp - 1) * np * np


@inline
def _mat22_idx(iidx: int, jidx: int, row: int, col: int, np: int) -> int:
    """Dinv and D declared as (np,np,2,2)."""
    return (iidx - 1) + (jidx - 1) * np + (row - 1) * np * np + (col - 1) * np * np * 2


@inline
def _mat32_idx(iidx: int, jidx: int, row: int, col: int, np: int) -> int:
    """vec_sphere2cart declared as (np,np,3,2)."""
    return (iidx - 1) + (jidx - 1) * np + (row - 1) * np * np + (col - 1) * np * np * 3


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


@inline
def _lev_q_idx(klev: int, qidx: int, nlev: int) -> int:
    """qmin and qmax declared as (nlev,qsize)."""
    return (klev - 1) + (qidx - 1) * nlev


@inline
def _cell_lev_idx(cell: int, klev: int, np: int) -> int:
    """ptens, dpmass, and workspaces declared as (np*np,nlev)."""
    ncols = np * np
    return (cell - 1) + (klev - 1) * ncols


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
def euler_step_qminmax_update_codon(
    np: int,
    nlev: int,
    qsize: int,
    rhs_multiplier: int,
    qtens_biharmonic_p: cobj,
    qmin_p: cobj,
    qmax_p: cobj,
):
    qtens_biharmonic = Ptr[float](qtens_biharmonic_p)
    qmin = Ptr[float](qmin_p)
    qmax = Ptr[float](qmax_p)

    for q in range(1, qsize + 1):
        for k in range(1, nlev + 1):
            qmin_val = 1.0e24
            qmax_val = -1.0e24
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    val = qtens_biharmonic[_q_idx(i, j, k, q, np, nlev)]
                    qmin_val = min(qmin_val, val)
                    qmax_val = max(qmax_val, val)

            lev_q_idx = _lev_q_idx(k, q, nlev)
            if rhs_multiplier == 1:
                qmin[lev_q_idx] = min(qmin[lev_q_idx], qmin_val)
                qmin[lev_q_idx] = max(qmin[lev_q_idx], 0.0)
                qmax[lev_q_idx] = max(qmax[lev_q_idx], qmax_val)
            else:
                qmin[lev_q_idx] = max(qmin_val, 0.0)
                qmax[lev_q_idx] = qmax_val


@export
def limiter2d_zero_codon(
    np: int,
    nlev: int,
    q_p: cobj,
):
    q = Ptr[float](q_p)

    for k in range(nlev, 0, -1):
        mass = 0.0
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                plane_idx = _vol_idx(i, j, k, np)
                mass = mass + q[plane_idx]

        if mass < 0.0:
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    plane_idx = _vol_idx(i, j, k, np)
                    q[plane_idx] = -q[plane_idx]

        mass_new = 0.0
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                plane_idx = _vol_idx(i, j, k, np)
                if q[plane_idx] < 0.0:
                    q[plane_idx] = 0.0
                else:
                    mass_new = mass_new + q[plane_idx]

        if mass_new > 0.0:
            scale = abs(mass) / mass_new
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    plane_idx = _vol_idx(i, j, k, np)
                    q[plane_idx] = q[plane_idx] * scale

        if mass < 0.0:
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    plane_idx = _vol_idx(i, j, k, np)
                    q[plane_idx] = -q[plane_idx]


@export
def limiter_optim_iter_full_codon(
    np: int,
    nlev: int,
    ptens_p: cobj,
    sphweights_p: cobj,
    minp_p: cobj,
    maxp_p: cobj,
    dpmass_p: cobj,
    weights_p: cobj,
    whois_neg_p: cobj,
    whois_pos_p: cobj,
    x_p: cobj,
    c_p: cobj,
    al_neg_p: cobj,
    al_pos_p: cobj,
):
    ptens = Ptr[float](ptens_p)
    sphweights = Ptr[float](sphweights_p)
    minp = Ptr[float](minp_p)
    maxp = Ptr[float](maxp_p)
    dpmass = Ptr[float](dpmass_p)
    weights = Ptr[float](weights_p)
    whois_neg = Ptr[int](whois_neg_p)
    whois_pos = Ptr[int](whois_pos_p)
    x = Ptr[float](x_p)
    c = Ptr[float](c_p)
    al_neg = Ptr[float](al_neg_p)
    al_pos = Ptr[float](al_pos_p)

    ncols = np * np
    tol_limiter = 1.0e-15
    maxiter = 5

    for k in range(1, nlev + 1):
        for k1 in range(1, ncols + 1):
            wk_idx = _cell_lev_idx(k1, k, np)
            weights[wk_idx] = sphweights[k1 - 1] * dpmass[wk_idx]
            ptens[wk_idx] = ptens[wk_idx] / dpmass[wk_idx]

    for k in range(1, nlev + 1):
        mass = 0.0
        sumc = 0.0
        for k1 in range(1, ncols + 1):
            wk_idx = _cell_lev_idx(k1, k, np)
            c[k1 - 1] = weights[wk_idx]
            x[k1 - 1] = ptens[wk_idx]
            mass = mass + c[k1 - 1] * x[k1 - 1]
            sumc = sumc + c[k1 - 1]

        if (mass / sumc) < minp[_col_idx(k)]:
            minp[_col_idx(k)] = mass / sumc
        if (mass / sumc) > maxp[_col_idx(k)]:
            maxp[_col_idx(k)] = mass / sumc

        addmass = 0.0
        pos_counter = 0
        neg_counter = 0

        for k1 in range(1, ncols + 1):
            if x[k1 - 1] >= maxp[_col_idx(k)]:
                addmass = addmass + (x[k1 - 1] - maxp[_col_idx(k)]) * c[k1 - 1]
                x[k1 - 1] = maxp[_col_idx(k)]
                whois_pos[k1 - 1] = -1
            else:
                pos_counter = pos_counter + 1
                whois_pos[pos_counter - 1] = k1

            if x[k1 - 1] <= minp[_col_idx(k)]:
                addmass = addmass - (minp[_col_idx(k)] - x[k1 - 1]) * c[k1 - 1]
                x[k1 - 1] = minp[_col_idx(k)]
                whois_neg[k1 - 1] = -1
            else:
                neg_counter = neg_counter + 1
                whois_neg[neg_counter - 1] = k1

        if addmass > 0.0:
            for _iter in range(1, maxiter + 1):
                weightssum = 0.0
                for k1 in range(1, pos_counter + 1):
                    i1 = whois_pos[k1 - 1]
                    weightssum = weightssum + c[i1 - 1]
                    al_pos[i1 - 1] = maxp[_col_idx(k)] - x[i1 - 1]

                if pos_counter > 0 and addmass > tol_limiter * abs(mass):
                    for k1 in range(1, pos_counter + 1):
                        i1 = whois_pos[k1 - 1]
                        howmuch = addmass / weightssum
                        if howmuch > al_pos[i1 - 1]:
                            howmuch = al_pos[i1 - 1]
                            whois_pos[k1 - 1] = -1
                        addmass = addmass - howmuch * c[i1 - 1]
                        weightssum = weightssum - c[i1 - 1]
                        x[i1 - 1] = x[i1 - 1] + howmuch

                    neg_counter = pos_counter
                    for k1 in range(1, ncols + 1):
                        whois_neg[k1 - 1] = whois_pos[k1 - 1]
                        whois_pos[k1 - 1] = -1
                    pos_counter = 0
                    for k1 in range(1, neg_counter + 1):
                        if whois_neg[k1 - 1] != -1:
                            pos_counter = pos_counter + 1
                            whois_pos[pos_counter - 1] = whois_neg[k1 - 1]
                else:
                    break
        else:
            for _iter in range(1, maxiter + 1):
                weightssum = 0.0
                for k1 in range(1, neg_counter + 1):
                    i1 = whois_neg[k1 - 1]
                    weightssum = weightssum + c[i1 - 1]
                    al_neg[i1 - 1] = x[i1 - 1] - minp[_col_idx(k)]

                if neg_counter > 0 and (-addmass) > tol_limiter * abs(mass):
                    for k1 in range(1, neg_counter + 1):
                        i1 = whois_neg[k1 - 1]
                        howmuch = -addmass / weightssum
                        if howmuch > al_neg[i1 - 1]:
                            howmuch = al_neg[i1 - 1]
                            whois_neg[k1 - 1] = -1
                        addmass = addmass + howmuch * c[i1 - 1]
                        weightssum = weightssum - c[i1 - 1]
                        x[i1 - 1] = x[i1 - 1] - howmuch

                    pos_counter = neg_counter
                    for k1 in range(1, ncols + 1):
                        whois_pos[k1 - 1] = whois_neg[k1 - 1]
                        whois_neg[k1 - 1] = -1
                    neg_counter = 0
                    for k1 in range(1, pos_counter + 1):
                        if whois_pos[k1 - 1] != -1:
                            neg_counter = neg_counter + 1
                            whois_neg[neg_counter - 1] = whois_pos[k1 - 1]
                else:
                    break

        for k1 in range(1, ncols + 1):
            ptens[_cell_lev_idx(k1, k, np)] = x[k1 - 1]

    for k in range(1, nlev + 1):
        for k1 in range(1, ncols + 1):
            wk_idx = _cell_lev_idx(k1, k, np)
            ptens[wk_idx] = ptens[wk_idx] * dpmass[wk_idx]


@export
def divergence_sphere_codon(
    np: int,
    rrearth: float,
    v_p: cobj,
    dvv_p: cobj,
    metdet_p: cobj,
    dinv_p: cobj,
    rmetdet_p: cobj,
    gv_p: cobj,
    vvtemp_p: cobj,
    div_p: cobj,
):
    v = Ptr[float](v_p)
    dvv = Ptr[float](dvv_p)
    metdet = Ptr[float](metdet_p)
    dinv = Ptr[float](dinv_p)
    rmetdet = Ptr[float](rmetdet_p)
    gv = Ptr[float](gv_p)
    vvtemp = Ptr[float](vvtemp_p)
    div = Ptr[float](div_p)

    for j in range(1, np + 1):
        for i in range(1, np + 1):
            plane_idx = _plane_idx(i, j, np)
            v1 = v[_vec2_idx(i, j, 1, np)]
            v2 = v[_vec2_idx(i, j, 2, np)]
            gv[_vec2_idx(i, j, 1, np)] = metdet[plane_idx] * (
                dinv[_mat22_idx(i, j, 1, 1, np)] * v1 + dinv[_mat22_idx(i, j, 1, 2, np)] * v2
            )
            gv[_vec2_idx(i, j, 2, np)] = metdet[plane_idx] * (
                dinv[_mat22_idx(i, j, 2, 1, np)] * v1 + dinv[_mat22_idx(i, j, 2, 2, np)] * v2
            )

    for j in range(1, np + 1):
        for l in range(1, np + 1):
            dudx00 = 0.0
            dvdy00 = 0.0
            for i in range(1, np + 1):
                dudx00 = dudx00 + dvv[_plane_idx(i, l, np)] * gv[_vec2_idx(i, j, 1, np)]
                dvdy00 = dvdy00 + dvv[_plane_idx(i, l, np)] * gv[_vec2_idx(j, i, 2, np)]
            div[_plane_idx(l, j, np)] = dudx00
            vvtemp[_plane_idx(j, l, np)] = dvdy00

    for j in range(1, np + 1):
        for i in range(1, np + 1):
            plane_idx = _plane_idx(i, j, np)
            div[plane_idx] = (div[plane_idx] + vvtemp[plane_idx]) * (rmetdet[plane_idx] * rrearth)


def _divergence_sphere_wk_vec(
    np: int,
    rrearth: float,
    v: Ptr[float],
    dvv: Ptr[float],
    spheremp: Ptr[float],
    dinv: Ptr[float],
    vtemp: Ptr[float],
    div: Ptr[float],
):
    for j in range(1, np + 1):
        for i in range(1, np + 1):
            v1 = v[_vec2_idx(i, j, 1, np)]
            v2 = v[_vec2_idx(i, j, 2, np)]
            vtemp[_vec2_idx(i, j, 1, np)] = (
                dinv[_mat22_idx(i, j, 1, 1, np)] * v1 + dinv[_mat22_idx(i, j, 1, 2, np)] * v2
            )
            vtemp[_vec2_idx(i, j, 2, np)] = (
                dinv[_mat22_idx(i, j, 2, 1, np)] * v1 + dinv[_mat22_idx(i, j, 2, 2, np)] * v2
            )

    for n in range(1, np + 1):
        for m in range(1, np + 1):
            div_idx = _plane_idx(m, n, np)
            div[div_idx] = 0.0
            for j in range(1, np + 1):
                div[div_idx] = div[div_idx] - (
                    spheremp[_plane_idx(j, n, np)] * vtemp[_vec2_idx(j, n, 1, np)] * dvv[_plane_idx(m, j, np)]
                    + spheremp[_plane_idx(m, j, np)] * vtemp[_vec2_idx(m, j, 2, np)] * dvv[_plane_idx(n, j, np)]
                ) * rrearth


def _gradient_sphere_field(
    np: int,
    rrearth: float,
    s: Ptr[float],
    s_offset: int,
    dvv: Ptr[float],
    dinv: Ptr[float],
    v1: Ptr[float],
    v2: Ptr[float],
    ds: Ptr[float],
):
    for j in range(1, np + 1):
        for l in range(1, np + 1):
            dsdx00 = 0.0
            dsdy00 = 0.0
            for i in range(1, np + 1):
                dsdx00 = dsdx00 + dvv[_plane_idx(i, l, np)] * s[s_offset + _plane_idx(i, j, np)]
                dsdy00 = dsdy00 + dvv[_plane_idx(i, l, np)] * s[s_offset + _plane_idx(j, i, np)]
            v1[_plane_idx(l, j, np)] = dsdx00 * rrearth
            v2[_plane_idx(j, l, np)] = dsdy00 * rrearth

    for j in range(1, np + 1):
        for i in range(1, np + 1):
            plane_idx = _plane_idx(i, j, np)
            v1_val = v1[plane_idx]
            v2_val = v2[plane_idx]
            ds[_vec2_idx(i, j, 1, np)] = (
                dinv[_mat22_idx(i, j, 1, 1, np)] * v1_val + dinv[_mat22_idx(i, j, 2, 1, np)] * v2_val
            )
            ds[_vec2_idx(i, j, 2, np)] = (
                dinv[_mat22_idx(i, j, 1, 2, np)] * v1_val + dinv[_mat22_idx(i, j, 2, 2, np)] * v2_val
            )


def _laplace_sphere_wk_field(
    np: int,
    rrearth: float,
    hypervis_power: int,
    hypervis_scaling: int,
    var_coef: int,
    s: Ptr[float],
    s_offset: int,
    dvv: Ptr[float],
    spheremp: Ptr[float],
    dinv: Ptr[float],
    variable_hyperviscosity: Ptr[float],
    tensorvisc: Ptr[float],
    grads: Ptr[float],
    oldgrads: Ptr[float],
    v1: Ptr[float],
    v2: Ptr[float],
    laplace: Ptr[float],
):
    _gradient_sphere_field(np, rrearth, s, s_offset, dvv, dinv, v1, v2, grads)

    if var_coef != 0:
        if hypervis_power != 0:
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    plane_idx = _plane_idx(i, j, np)
                    scale = variable_hyperviscosity[plane_idx]
                    grads[_vec2_idx(i, j, 1, np)] = grads[_vec2_idx(i, j, 1, np)] * scale
                    grads[_vec2_idx(i, j, 2, np)] = grads[_vec2_idx(i, j, 2, np)] * scale
        elif hypervis_scaling != 0:
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    oldgrads[_vec2_idx(i, j, 1, np)] = grads[_vec2_idx(i, j, 1, np)]
                    oldgrads[_vec2_idx(i, j, 2, np)] = grads[_vec2_idx(i, j, 2, np)]

            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    oldgrad1 = oldgrads[_vec2_idx(i, j, 1, np)]
                    oldgrad2 = oldgrads[_vec2_idx(i, j, 2, np)]
                    grads[_vec2_idx(i, j, 1, np)] = (
                        oldgrad1 * tensorvisc[_mat22_idx(i, j, 1, 1, np)]
                        + oldgrad2 * tensorvisc[_mat22_idx(i, j, 1, 2, np)]
                    )
                    grads[_vec2_idx(i, j, 2, np)] = (
                        oldgrad1 * tensorvisc[_mat22_idx(i, j, 2, 1, np)]
                        + oldgrad2 * tensorvisc[_mat22_idx(i, j, 2, 2, np)]
                    )

    _divergence_sphere_wk_vec(np, rrearth, grads, dvv, spheremp, dinv, oldgrads, laplace)


def _curl_sphere_wk_testcov(
    np: int,
    rrearth: float,
    s: Ptr[float],
    dvv: Ptr[float],
    mp: Ptr[float],
    d: Ptr[float],
    dscontra: Ptr[float],
    ds: Ptr[float],
):
    for n in range(1, np + 1):
        for m in range(1, np + 1):
            dscontra[_vec2_idx(m, n, 1, np)] = 0.0
            dscontra[_vec2_idx(m, n, 2, np)] = 0.0
            for j in range(1, np + 1):
                dscontra[_vec2_idx(m, n, 1, np)] = dscontra[_vec2_idx(m, n, 1, np)] - (
                    mp[_plane_idx(m, j, np)] * s[_plane_idx(m, j, np)] * dvv[_plane_idx(n, j, np)]
                ) * rrearth
                dscontra[_vec2_idx(m, n, 2, np)] = dscontra[_vec2_idx(m, n, 2, np)] + (
                    mp[_plane_idx(j, n, np)] * s[_plane_idx(j, n, np)] * dvv[_plane_idx(m, j, np)]
                ) * rrearth

    for j in range(1, np + 1):
        for i in range(1, np + 1):
            ds[_vec2_idx(i, j, 1, np)] = (
                d[_mat22_idx(i, j, 1, 1, np)] * dscontra[_vec2_idx(i, j, 1, np)]
                + d[_mat22_idx(i, j, 1, 2, np)] * dscontra[_vec2_idx(i, j, 2, np)]
            )
            ds[_vec2_idx(i, j, 2, np)] = (
                d[_mat22_idx(i, j, 2, 1, np)] * dscontra[_vec2_idx(i, j, 1, np)]
                + d[_mat22_idx(i, j, 2, 2, np)] * dscontra[_vec2_idx(i, j, 2, np)]
            )


def _gradient_sphere_wk_testcov(
    np: int,
    rrearth: float,
    s: Ptr[float],
    dvv: Ptr[float],
    mp: Ptr[float],
    metinv: Ptr[float],
    metdet: Ptr[float],
    d: Ptr[float],
    dscontra: Ptr[float],
    ds: Ptr[float],
):
    for n in range(1, np + 1):
        for m in range(1, np + 1):
            dscontra[_vec2_idx(m, n, 1, np)] = 0.0
            dscontra[_vec2_idx(m, n, 2, np)] = 0.0
            for j in range(1, np + 1):
                plane_idx = _plane_idx(m, n, np)
                dscontra[_vec2_idx(m, n, 1, np)] = dscontra[_vec2_idx(m, n, 1, np)] - (
                    (
                        mp[_plane_idx(j, n, np)]
                        * metinv[_mat22_idx(m, n, 1, 1, np)]
                        * metdet[plane_idx]
                        * s[_plane_idx(j, n, np)]
                        * dvv[_plane_idx(m, j, np)]
                    )
                    + (
                        mp[_plane_idx(m, j, np)]
                        * metinv[_mat22_idx(m, n, 2, 1, np)]
                        * metdet[plane_idx]
                        * s[_plane_idx(m, j, np)]
                        * dvv[_plane_idx(n, j, np)]
                    )
                ) * rrearth
                dscontra[_vec2_idx(m, n, 2, np)] = dscontra[_vec2_idx(m, n, 2, np)] - (
                    (
                        mp[_plane_idx(j, n, np)]
                        * metinv[_mat22_idx(m, n, 1, 2, np)]
                        * metdet[plane_idx]
                        * s[_plane_idx(j, n, np)]
                        * dvv[_plane_idx(m, j, np)]
                    )
                    + (
                        mp[_plane_idx(m, j, np)]
                        * metinv[_mat22_idx(m, n, 2, 2, np)]
                        * metdet[plane_idx]
                        * s[_plane_idx(m, j, np)]
                        * dvv[_plane_idx(n, j, np)]
                    )
                ) * rrearth

    for j in range(1, np + 1):
        for i in range(1, np + 1):
            ds[_vec2_idx(i, j, 1, np)] = (
                d[_mat22_idx(i, j, 1, 1, np)] * dscontra[_vec2_idx(i, j, 1, np)]
                + d[_mat22_idx(i, j, 1, 2, np)] * dscontra[_vec2_idx(i, j, 2, np)]
            )
            ds[_vec2_idx(i, j, 2, np)] = (
                d[_mat22_idx(i, j, 2, 1, np)] * dscontra[_vec2_idx(i, j, 1, np)]
                + d[_mat22_idx(i, j, 2, 2, np)] * dscontra[_vec2_idx(i, j, 2, np)]
            )


@export
def divergence_sphere_wk_codon(
    np: int,
    rrearth: float,
    v_p: cobj,
    dvv_p: cobj,
    spheremp_p: cobj,
    dinv_p: cobj,
    vtemp_p: cobj,
    div_p: cobj,
):
    _divergence_sphere_wk_vec(
        np,
        rrearth,
        Ptr[float](v_p),
        Ptr[float](dvv_p),
        Ptr[float](spheremp_p),
        Ptr[float](dinv_p),
        Ptr[float](vtemp_p),
        Ptr[float](div_p),
    )


@export
def laplace_sphere_wk_codon(
    np: int,
    rrearth: float,
    hypervis_power: int,
    hypervis_scaling: int,
    var_coef: int,
    s_p: cobj,
    dvv_p: cobj,
    spheremp_p: cobj,
    dinv_p: cobj,
    variable_hyperviscosity_p: cobj,
    tensorvisc_p: cobj,
    grads_p: cobj,
    oldgrads_p: cobj,
    v1_p: cobj,
    v2_p: cobj,
    laplace_p: cobj,
):
    _laplace_sphere_wk_field(
        np,
        rrearth,
        hypervis_power,
        hypervis_scaling,
        var_coef,
        Ptr[float](s_p),
        0,
        Ptr[float](dvv_p),
        Ptr[float](spheremp_p),
        Ptr[float](dinv_p),
        Ptr[float](variable_hyperviscosity_p),
        Ptr[float](tensorvisc_p),
        Ptr[float](grads_p),
        Ptr[float](oldgrads_p),
        Ptr[float](v1_p),
        Ptr[float](v2_p),
        Ptr[float](laplace_p),
    )


@export
def vlaplace_sphere_wk_codon(
    np: int,
    rrearth: float,
    hypervis_power: int,
    hypervis_scaling: int,
    var_coef: int,
    has_nu_ratio: int,
    nu_ratio: float,
    v_p: cobj,
    dvv_p: cobj,
    mp_p: cobj,
    spheremp_p: cobj,
    metinv_p: cobj,
    metdet_p: cobj,
    rmetdet_p: cobj,
    d_p: cobj,
    dinv_p: cobj,
    variable_hyperviscosity_p: cobj,
    tensorvisc_p: cobj,
    vec_sphere2cart_p: cobj,
    dum_cart_p: cobj,
    dum_tmp_p: cobj,
    div_p: cobj,
    vor_p: cobj,
    lap_tmp_p: cobj,
    lap_tmp2_p: cobj,
    work1_p: cobj,
    work2_p: cobj,
    v1_p: cobj,
    v2_p: cobj,
    laplace_p: cobj,
):
    v = Ptr[float](v_p)
    mp = Ptr[float](mp_p)
    spheremp = Ptr[float](spheremp_p)
    metinv = Ptr[float](metinv_p)
    metdet = Ptr[float](metdet_p)
    d = Ptr[float](d_p)
    variable_hyperviscosity = Ptr[float](variable_hyperviscosity_p)
    vec_sphere2cart = Ptr[float](vec_sphere2cart_p)
    dum_cart = Ptr[float](dum_cart_p)
    dum_tmp = Ptr[float](dum_tmp_p)
    div = Ptr[float](div_p)
    vor = Ptr[float](vor_p)
    lap_tmp = Ptr[float](lap_tmp_p)
    lap_tmp2 = Ptr[float](lap_tmp2_p)
    work1 = Ptr[float](work1_p)
    work2 = Ptr[float](work2_p)
    laplace = Ptr[float](laplace_p)

    if hypervis_scaling != 0 and var_coef != 0:
        for component in range(1, 4):
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    dum_cart[_vec3_idx(i, j, component, np)] = (
                        vec_sphere2cart[_mat32_idx(i, j, component, 1, np)] * v[_vec2_idx(i, j, 1, np)]
                        + vec_sphere2cart[_mat32_idx(i, j, component, 2, np)] * v[_vec2_idx(i, j, 2, np)]
                    )

        for component in range(1, 4):
            _laplace_sphere_wk_field(
                np,
                rrearth,
                hypervis_power,
                hypervis_scaling,
                var_coef,
                dum_cart,
                (component - 1) * np * np,
                Ptr[float](dvv_p),
                spheremp,
                Ptr[float](dinv_p),
                variable_hyperviscosity,
                Ptr[float](tensorvisc_p),
                work1,
                work2,
                Ptr[float](v1_p),
                Ptr[float](v2_p),
                dum_tmp,
            )
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    dum_cart[_vec3_idx(i, j, component, np)] = dum_tmp[_plane_idx(i, j, np)]

        for component in range(1, 3):
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    laplace[_vec2_idx(i, j, component, np)] = (
                        dum_cart[_vec3_idx(i, j, 1, np)] * vec_sphere2cart[_mat32_idx(i, j, 1, component, np)]
                        + dum_cart[_vec3_idx(i, j, 2, np)] * vec_sphere2cart[_mat32_idx(i, j, 2, component, np)]
                        + dum_cart[_vec3_idx(i, j, 3, np)] * vec_sphere2cart[_mat32_idx(i, j, 3, component, np)]
                    )
    else:
        divergence_sphere_codon(np, rrearth, v_p, dvv_p, metdet_p, dinv_p, rmetdet_p, lap_tmp_p, dum_tmp_p, div_p)
        vorticity_sphere_codon(np, rrearth, v_p, dvv_p, d_p, rmetdet_p, lap_tmp2_p, dum_tmp_p, vor_p)

        if var_coef != 0 and hypervis_power != 0:
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    plane_idx = _plane_idx(i, j, np)
                    scale = variable_hyperviscosity[plane_idx]
                    div[plane_idx] = div[plane_idx] * scale
                    vor[plane_idx] = vor[plane_idx] * scale

        if has_nu_ratio != 0:
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    plane_idx = _plane_idx(i, j, np)
                    div[plane_idx] = nu_ratio * div[plane_idx]

        _gradient_sphere_wk_testcov(np, rrearth, div, Ptr[float](dvv_p), mp, metinv, metdet, d, work1, lap_tmp)
        _curl_sphere_wk_testcov(np, rrearth, vor, Ptr[float](dvv_p), mp, d, work2, lap_tmp2)

        rrearth_sq = rrearth * rrearth
        for n in range(1, np + 1):
            for m in range(1, np + 1):
                plane_idx = _plane_idx(m, n, np)
                laplace[_vec2_idx(m, n, 1, np)] = (
                    lap_tmp[_vec2_idx(m, n, 1, np)]
                    - lap_tmp2[_vec2_idx(m, n, 1, np)]
                    + 2.0 * spheremp[plane_idx] * v[_vec2_idx(m, n, 1, np)] * rrearth_sq
                )
                laplace[_vec2_idx(m, n, 2, np)] = (
                    lap_tmp[_vec2_idx(m, n, 2, np)]
                    - lap_tmp2[_vec2_idx(m, n, 2, np)]
                    + 2.0 * spheremp[plane_idx] * v[_vec2_idx(m, n, 2, np)] * rrearth_sq
                )


@export
def vorticity_sphere_codon(
    np: int,
    rrearth: float,
    v_p: cobj,
    dvv_p: cobj,
    d_p: cobj,
    rmetdet_p: cobj,
    vco_p: cobj,
    vtemp_p: cobj,
    vort_p: cobj,
):
    v = Ptr[float](v_p)
    dvv = Ptr[float](dvv_p)
    d = Ptr[float](d_p)
    rmetdet = Ptr[float](rmetdet_p)
    vco = Ptr[float](vco_p)
    vtemp = Ptr[float](vtemp_p)
    vort = Ptr[float](vort_p)

    for j in range(1, np + 1):
        for i in range(1, np + 1):
            v1 = v[_vec2_idx(i, j, 1, np)]
            v2 = v[_vec2_idx(i, j, 2, np)]
            vco[_vec2_idx(i, j, 1, np)] = (
                d[_mat22_idx(i, j, 1, 1, np)] * v1 + d[_mat22_idx(i, j, 2, 1, np)] * v2
            )
            vco[_vec2_idx(i, j, 2, np)] = (
                d[_mat22_idx(i, j, 1, 2, np)] * v1 + d[_mat22_idx(i, j, 2, 2, np)] * v2
            )

    for j in range(1, np + 1):
        for l in range(1, np + 1):
            dudy00 = 0.0
            dvdx00 = 0.0
            for i in range(1, np + 1):
                dvdx00 = dvdx00 + dvv[_plane_idx(i, l, np)] * vco[_vec2_idx(i, j, 2, np)]
                dudy00 = dudy00 + dvv[_plane_idx(i, l, np)] * vco[_vec2_idx(j, i, 1, np)]
            vort[_plane_idx(l, j, np)] = dvdx00
            vtemp[_plane_idx(j, l, np)] = dudy00

    for j in range(1, np + 1):
        for i in range(1, np + 1):
            plane_idx = _plane_idx(i, j, np)
            vort[plane_idx] = (vort[plane_idx] - vtemp[plane_idx]) * (rmetdet[plane_idx] * rrearth)


@export
def gradient_sphere_codon(
    np: int,
    rrearth: float,
    s_p: cobj,
    dvv_p: cobj,
    dinv_p: cobj,
    v1_p: cobj,
    v2_p: cobj,
    ds_p: cobj,
):
    s = Ptr[float](s_p)
    dvv = Ptr[float](dvv_p)
    dinv = Ptr[float](dinv_p)
    v1 = Ptr[float](v1_p)
    v2 = Ptr[float](v2_p)
    ds = Ptr[float](ds_p)

    for j in range(1, np + 1):
        for l in range(1, np + 1):
            dsdx00 = 0.0
            dsdy00 = 0.0
            for i in range(1, np + 1):
                dsdx00 = dsdx00 + dvv[_plane_idx(i, l, np)] * s[_plane_idx(i, j, np)]
                dsdy00 = dsdy00 + dvv[_plane_idx(i, l, np)] * s[_plane_idx(j, i, np)]
            v1[_plane_idx(l, j, np)] = dsdx00 * rrearth
            v2[_plane_idx(j, l, np)] = dsdy00 * rrearth

    for j in range(1, np + 1):
        for i in range(1, np + 1):
            plane_idx = _plane_idx(i, j, np)
            v1_val = v1[plane_idx]
            v2_val = v2[plane_idx]
            ds[_vec2_idx(i, j, 1, np)] = (
                dinv[_mat22_idx(i, j, 1, 1, np)] * v1_val + dinv[_mat22_idx(i, j, 2, 1, np)] * v2_val
            )
            ds[_vec2_idx(i, j, 2, np)] = (
                dinv[_mat22_idx(i, j, 1, 2, np)] * v1_val + dinv[_mat22_idx(i, j, 2, 2, np)] * v2_val
            )


@export
def curl_sphere_codon(
    np: int,
    rrearth: float,
    s_p: cobj,
    dvv_p: cobj,
    d_p: cobj,
    metdet_p: cobj,
    v1_p: cobj,
    v2_p: cobj,
    ds_p: cobj,
):
    s = Ptr[float](s_p)
    dvv = Ptr[float](dvv_p)
    d = Ptr[float](d_p)
    metdet = Ptr[float](metdet_p)
    v1 = Ptr[float](v1_p)
    v2 = Ptr[float](v2_p)
    ds = Ptr[float](ds_p)

    for j in range(1, np + 1):
        for l in range(1, np + 1):
            dsdx00 = 0.0
            dsdy00 = 0.0
            for i in range(1, np + 1):
                dsdx00 = dsdx00 + dvv[_plane_idx(i, l, np)] * s[_plane_idx(i, j, np)]
                dsdy00 = dsdy00 + dvv[_plane_idx(i, l, np)] * s[_plane_idx(j, i, np)]
            v2[_plane_idx(l, j, np)] = -dsdx00 * rrearth
            v1[_plane_idx(j, l, np)] = dsdy00 * rrearth

    for j in range(1, np + 1):
        for i in range(1, np + 1):
            plane_idx = _plane_idx(i, j, np)
            v1_val = v1[plane_idx]
            v2_val = v2[plane_idx]
            ds[_vec2_idx(i, j, 1, np)] = (
                d[_mat22_idx(i, j, 1, 1, np)] * v1_val + d[_mat22_idx(i, j, 1, 2, np)] * v2_val
            ) / metdet[plane_idx]
            ds[_vec2_idx(i, j, 2, np)] = (
                d[_mat22_idx(i, j, 2, 1, np)] * v1_val + d[_mat22_idx(i, j, 2, 2, np)] * v2_val
            ) / metdet[plane_idx]


@export
def ugradv_sphere_codon(
    np: int,
    rrearth: float,
    u_p: cobj,
    v_p: cobj,
    dvv_p: cobj,
    dinv_p: cobj,
    vec_sphere2cart_p: cobj,
    dum_cart_p: cobj,
    tmp_p: cobj,
    v1_p: cobj,
    v2_p: cobj,
    ugradv_p: cobj,
):
    u = Ptr[float](u_p)
    v = Ptr[float](v_p)
    dvv = Ptr[float](dvv_p)
    dinv = Ptr[float](dinv_p)
    vec_sphere2cart = Ptr[float](vec_sphere2cart_p)
    dum_cart = Ptr[float](dum_cart_p)
    tmp = Ptr[float](tmp_p)
    v1 = Ptr[float](v1_p)
    v2 = Ptr[float](v2_p)
    ugradv = Ptr[float](ugradv_p)

    for component in range(1, 4):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                dum_cart[_vec3_idx(i, j, component, np)] = (
                    vec_sphere2cart[_mat32_idx(i, j, component, 1, np)] * v[_vec2_idx(i, j, 1, np)]
                    + vec_sphere2cart[_mat32_idx(i, j, component, 2, np)] * v[_vec2_idx(i, j, 2, np)]
                )

    for component in range(1, 4):
        for j in range(1, np + 1):
            for l in range(1, np + 1):
                dsdx00 = 0.0
                dsdy00 = 0.0
                for i in range(1, np + 1):
                    dsdx00 = dsdx00 + dvv[_plane_idx(i, l, np)] * dum_cart[_vec3_idx(i, j, component, np)]
                    dsdy00 = dsdy00 + dvv[_plane_idx(i, l, np)] * dum_cart[_vec3_idx(j, i, component, np)]
                v1[_plane_idx(l, j, np)] = dsdx00 * rrearth
                v2[_plane_idx(j, l, np)] = dsdy00 * rrearth

        for j in range(1, np + 1):
            for i in range(1, np + 1):
                plane_idx = _plane_idx(i, j, np)
                v1_val = v1[plane_idx]
                v2_val = v2[plane_idx]
                tmp[_vec2_idx(i, j, 1, np)] = (
                    dinv[_mat22_idx(i, j, 1, 1, np)] * v1_val + dinv[_mat22_idx(i, j, 2, 1, np)] * v2_val
                )
                tmp[_vec2_idx(i, j, 2, np)] = (
                    dinv[_mat22_idx(i, j, 1, 2, np)] * v1_val + dinv[_mat22_idx(i, j, 2, 2, np)] * v2_val
                )
                dum_cart[_vec3_idx(i, j, component, np)] = (
                    u[_vec2_idx(i, j, 1, np)] * tmp[_vec2_idx(i, j, 1, np)]
                    + u[_vec2_idx(i, j, 2, np)] * tmp[_vec2_idx(i, j, 2, np)]
                )

    for component in range(1, 3):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                ugradv[_vec2_idx(i, j, component, np)] = (
                    dum_cart[_vec3_idx(i, j, 1, np)] * vec_sphere2cart[_mat32_idx(i, j, 1, component, np)]
                    + dum_cart[_vec3_idx(i, j, 2, np)] * vec_sphere2cart[_mat32_idx(i, j, 2, component, np)]
                    + dum_cart[_vec3_idx(i, j, 3, np)] * vec_sphere2cart[_mat32_idx(i, j, 3, component, np)]
                )


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
