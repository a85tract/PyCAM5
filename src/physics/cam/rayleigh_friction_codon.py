from math import tanh


@inline
def _state_idx(icol: int, klev: int, psetcols: int) -> int:
    """state/ptend arrays declared as (psetcols, pver)"""
    return (icol - 1) + (klev - 1) * psetcols


@inline
def _otau_idx(klev: int) -> int:
    """otau(pver)"""
    return klev - 1


@export
def rayleigh_friction_init_codon(
    rayk0: int,
    raykrange: float,
    raytau0: float,
    pver: int,
    krange_p: cobj,
    tau0_p: cobj,
    otau0_p: cobj,
    otau_p: cobj,
):
    krange = Ptr[float](krange_p)
    tau0 = Ptr[float](tau0_p)
    otau0 = Ptr[float](otau0_p)
    otau = Ptr[float](otau_p)

    krange[0] = raykrange
    if raykrange == 0.0:
        krange[0] = float(rayk0 - 1) / 2.0

    tau0[0] = 86400.0 * raytau0
    otau0[0] = 0.0
    if tau0[0] != 0.0:
        otau0[0] = 1.0 / tau0[0]

    for k in range(1, pver + 1):
        x = float(rayk0 - k) / krange[0]
        otau[_otau_idx(k)] = otau0[0] * (1.0 + tanh(x)) / 2.0


@export
def rayleigh_friction_tend_codon(
    ztodt: float,
    pver: int,
    psetcols: int,
    ncol: int,
    otau_p: cobj,
    state_u_p: cobj,
    state_v_p: cobj,
    ptend_u_p: cobj,
    ptend_v_p: cobj,
    ptend_s_p: cobj,
):
    otau = Ptr[float](otau_p)
    state_u = Ptr[float](state_u_p)
    state_v = Ptr[float](state_v_p)
    ptend_u = Ptr[float](ptend_u_p)
    ptend_v = Ptr[float](ptend_v_p)
    ptend_s = Ptr[float](ptend_s_p)

    rztodt = 1.0 / ztodt

    for k in range(1, pver + 1):
        otau_k = otau[_otau_idx(k)]
        c2 = 1.0 / (1.0 + otau_k * ztodt)
        c1 = -otau_k * c2
        c3 = 0.5 * (1.0 - c2 * c2) * rztodt

        for i in range(1, ncol + 1):
            idx = _state_idx(i, k, psetcols)
            u_val = state_u[idx]
            v_val = state_v[idx]
            ptend_u[idx] = c1 * u_val
            ptend_v[idx] = c1 * v_val
            ptend_s[idx] = c3 * (u_val * u_val + v_val * v_val)
