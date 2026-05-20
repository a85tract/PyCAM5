from math import copysign
import se_dynamics_common_codon as _common

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

    pio[_common._col_idx(nlev + 2)] = pio[_common._col_idx(nlev + 1)] + 1.0
    pin[_common._col_idx(nlev + 1)] = pio[_common._col_idx(nlev + 1)]

    for k in range(1, 3):
        dpo[_common._ghost_col_idx(1 - k)] = dpo[_common._ghost_col_idx(k)]
        dpo[_common._ghost_col_idx(nlev + k)] = dpo[_common._ghost_col_idx(nlev + 1 - k)]

    for k in range(1, nlev + 1):
        kk = k
        while pio[_common._col_idx(kk)] <= pin[_common._col_idx(k + 1)]:
            kk += 1
        kk -= 1
        if kk == nlev + 1:
            kk = nlev
        kid[_common._col_idx(k)] = kk
        z1[_common._col_idx(k)] = -0.5
        z2[_common._col_idx(k)] = (pin[_common._col_idx(k + 1)] - (pio[_common._col_idx(kk)] + pio[_common._col_idx(kk + 1)]) * 0.5) / dpo[
            _common._ghost_col_idx(kk)
        ]


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

    masso[_common._col_idx(1)] = 0.0
    for k in range(1, nlev + 1):
        ao[_common._ghost_col_idx(k)] = qdp[_common._q_idx(iidx, jidx, k, qidx, nx, nlev)]
        masso[_common._col_idx(k + 1)] = masso[_common._col_idx(k)] + ao[_common._ghost_col_idx(k)]
        ao[_common._ghost_col_idx(k)] = ao[_common._ghost_col_idx(k)] / dpo[_common._ghost_col_idx(k)]

    for k in range(1, 3):
        ao[_common._ghost_col_idx(1 - k)] = ao[_common._ghost_col_idx(k)]
        ao[_common._ghost_col_idx(nlev + k)] = ao[_common._ghost_col_idx(nlev + 1 - k)]


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
        dxjm1 = dx[_common._ghost_col_idx(j - 1)]
        dxj = dx[_common._ghost_col_idx(j)]
        dxjp1 = dx[_common._ghost_col_idx(j + 1)]
        rslt[_common._ppm_grid_idx(1, j)] = dxj / ((dxjm1 + dxj) + dxjp1)
        rslt[_common._ppm_grid_idx(2, j)] = ((2.0 * dxjm1) + dxj) / (dxjp1 + dxj)
        rslt[_common._ppm_grid_idx(3, j)] = (dxj + (2.0 * dxjp1)) / (dxjm1 + dxj)

    if vert_remap_q_alg == 2:
        indB = 2
        indE = nlev - 2
    else:
        indB = 0
        indE = nlev

    for j in range(indB, indE + 1):
        dxjm1 = dx[_common._ghost_col_idx(j - 1)]
        dxj = dx[_common._ghost_col_idx(j)]
        dxjp1 = dx[_common._ghost_col_idx(j + 1)]
        dxjp2 = dx[_common._ghost_col_idx(j + 2)]
        rslt[_common._ppm_grid_idx(4, j)] = dxj / (dxj + dxjp1)
        rslt[_common._ppm_grid_idx(5, j)] = 1.0 / (((dxjm1 + dxj) + dxjp1) + dxjp2)
        rslt[_common._ppm_grid_idx(6, j)] = ((2.0 * dxjp1) * dxj) / (dxj + dxjp1)
        rslt[_common._ppm_grid_idx(7, j)] = (dxjm1 + dxj) / ((2.0 * dxj) + dxjp1)
        rslt[_common._ppm_grid_idx(8, j)] = (dxjp2 + dxjp1) / ((2.0 * dxjp1) + dxj)
        rslt[_common._ppm_grid_idx(9, j)] = dxj * (dxjm1 + dxj) / ((2.0 * dxj) + dxjp1)
        rslt[_common._ppm_grid_idx(10, j)] = dxjp1 * (dxjp1 + dxjp2) / (dxj + (2.0 * dxjp1))


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
        ajm1 = a[_common._ghost_col_idx(j - 1)]
        aj = a[_common._ghost_col_idx(j)]
        ajp1 = a[_common._ghost_col_idx(j + 1)]
        da = dx[_common._ppm_grid_idx(1, j)] * (
            dx[_common._ppm_grid_idx(2, j)] * (ajp1 - aj) + dx[_common._ppm_grid_idx(3, j)] * (aj - ajm1)
        )
        dma_val = min(
            abs(da),
            min(2.0 * abs(aj - ajm1), 2.0 * abs(ajp1 - aj)),
        ) * copysign(1.0, da)
        if (ajp1 - aj) * (aj - ajm1) <= 0.0:
            dma_val = 0.0
        dma[_common._ppm_scratch_idx(j)] = dma_val

    if vert_remap_q_alg == 2:
        indB = 2
        indE = nlev - 2
    else:
        indB = 0
        indE = nlev

    for j in range(indB, indE + 1):
        aj = a[_common._ghost_col_idx(j)]
        ajp1 = a[_common._ghost_col_idx(j + 1)]
        ai[_common._ppm_scratch_idx(j)] = aj + dx[_common._ppm_grid_idx(4, j)] * (ajp1 - aj) + dx[_common._ppm_grid_idx(5, j)] * (
            dx[_common._ppm_grid_idx(6, j)] * (dx[_common._ppm_grid_idx(7, j)] - dx[_common._ppm_grid_idx(8, j)]) * (ajp1 - aj)
            - dx[_common._ppm_grid_idx(9, j)] * dma[_common._ppm_scratch_idx(j + 1)]
            + dx[_common._ppm_grid_idx(10, j)] * dma[_common._ppm_scratch_idx(j)]
        )

    if vert_remap_q_alg == 2:
        indB = 3
        indE = nlev - 2
    else:
        indB = 1
        indE = nlev

    for j in range(indB, indE + 1):
        aj = a[_common._ghost_col_idx(j)]
        al = ai[_common._ppm_scratch_idx(j - 1)]
        ar = ai[_common._ppm_scratch_idx(j)]
        if (ar - aj) * (aj - al) <= 0.0:
            al = aj
            ar = aj
        if (ar - al) * (aj - (al + ar) / 2.0) > ((ar - al) * (ar - al)) / 6.0:
            al = 3.0 * aj - 2.0 * ar
        if (ar - al) * (aj - (al + ar) / 2.0) < -(((ar - al) * (ar - al)) / 6.0):
            ar = 3.0 * aj - 2.0 * al
        coefs[_common._ppm_coef_idx(0, j)] = 1.5 * aj - (al + ar) / 4.0
        coefs[_common._ppm_coef_idx(1, j)] = ar - al
        coefs[_common._ppm_coef_idx(2, j)] = -6.0 * aj + 3.0 * (al + ar)

    if vert_remap_q_alg == 2:
        for j in range(1, 3):
            aj = a[_common._ghost_col_idx(j)]
            coefs[_common._ppm_coef_idx(0, j)] = aj
            coefs[_common._ppm_coef_idx(1, j)] = 0.0
            coefs[_common._ppm_coef_idx(2, j)] = 0.0
        for j in range(nlev - 1, nlev + 1):
            aj = a[_common._ghost_col_idx(j)]
            coefs[_common._ppm_coef_idx(0, j)] = aj
            coefs[_common._ppm_coef_idx(1, j)] = 0.0
            coefs[_common._ppm_coef_idx(2, j)] = 0.0


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
        kk = kid[_common._col_idx(k)]
        x1 = z1[_common._col_idx(k)]
        x2 = z2[_common._col_idx(k)]
        x1_sq = x1 * x1
        x2_sq = x2 * x2
        mass = (
            coefs[_common._ppm_coef_idx(0, kk)] * (x2 - x1)
            + coefs[_common._ppm_coef_idx(1, kk)] * (x2_sq - x1_sq) / 2.0
            + coefs[_common._ppm_coef_idx(2, kk)] * ((x2_sq * x2) - (x1_sq * x1)) / 3.0
        )
        massn2 = masso[_common._col_idx(kk)] + mass * dpo[_common._ghost_col_idx(kk)]
        qdp[_common._q_idx(iidx, jidx, k, qidx, nx, nlev)] = massn2 - massn1
        massn1 = massn2
