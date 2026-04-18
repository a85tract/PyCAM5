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
