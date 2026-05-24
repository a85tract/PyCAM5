from math import sin


@export
def aoa_tracers_readnl_codon(flag: int, read_from_ic: int) -> int:
    out = 0
    if flag != 0:
        out += 1
    if read_from_ic != 0:
        out += 2
    return out


@export
def aoa_tracers_flag_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def aoa_tracers_implements_cnst_codon(flag: int, name_len: int, name_ascii_p: cobj) -> int:
    if flag == 0:
        return 0
    name_ascii = Ptr[int](name_ascii_p)
    if _name_eq8(name_len, name_ascii, 65, 79, 65, 49, 32, 32, 32, 32) != 0:
        return 1
    if _name_eq8(name_len, name_ascii, 65, 79, 65, 50, 32, 32, 32, 32) != 0:
        return 1
    if _name_eq8(name_len, name_ascii, 72, 79, 82, 90, 32, 32, 32, 32) != 0:
        return 1
    if _name_eq8(name_len, name_ascii, 86, 69, 82, 84, 32, 32, 32, 32) != 0:
        return 1
    return 0


@inline
def _name_eq8(
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
def aoa_tracers_register_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def aoa_tracers_init_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@inline
def _col_idx(icol: int) -> int:
    """state%lat declared as (pcols)"""
    return icol - 1


@inline
def _vec_idx(idx1: int) -> int:
    """qrel_vert declared as (pver)"""
    return idx1 - 1


@inline
def _field3_idx(icol: int, klev: int, mconst: int, ld1: int, ld2: int) -> int:
    """state%q/ptend%q declared as (ld1, pver, pcnst)"""
    return (icol - 1) + (klev - 1) * ld1 + (mconst - 1) * ld1 * ld2


@inline
def _flux_idx(icol: int, mconst: int, pcols: int) -> int:
    """cflx declared as (pcols, pcnst)"""
    return (icol - 1) + (mconst - 1) * pcols


@export
def aoa_tracers_tstep_init_codon(
    ncol: int,
    pcols: int,
    pver: int,
    ixht: int,
    ixvt: int,
    qrel_vert_p: cobj,
    state_lat_p: cobj,
    state_q_p: cobj,
):
    qrel_vert = Ptr[float](qrel_vert_p)
    state_lat = Ptr[float](state_lat_p)
    state_q = Ptr[float](state_q_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            state_q[_field3_idx(i, k, ixht, pcols, pver)] = 2.0 + sin(
                state_lat[_col_idx(i)]
            )
            state_q[_field3_idx(i, k, ixvt, pcols, pver)] = qrel_vert[
                _vec_idx(k)
            ]


@export
def aoa_tracers_timestep_tend_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    psetcols: int,
    ixaoa1: int,
    ixaoa2: int,
    ixht: int,
    ixvt: int,
    nstep: int,
    dt: float,
    qrel_vert_p: cobj,
    state_lat_p: cobj,
    state_q_p: cobj,
    ptend_q_p: cobj,
    cflx_p: cobj,
    landfrac_p: cobj,
):
    qrel_vert = Ptr[float](qrel_vert_p)
    state_lat = Ptr[float](state_lat_p)
    state_q = Ptr[float](state_q_p)
    ptend_q = Ptr[float](ptend_q_p)
    cflx = Ptr[float](cflx_p)
    landfrac = Ptr[float](landfrac_p)

    teul = 0.5 * dt / (86400.0 * 15.0)
    wimp = 1.0 / (1.0 + teul)
    wsrc = teul * wimp

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            ptend_q[_field3_idx(i, k, ixaoa1, psetcols, pver)] = 0.0
            ptend_q[_field3_idx(i, k, ixaoa2, psetcols, pver)] = 0.0

            qrel = 2.0 + sin(state_lat[_col_idx(i)])
            xhorz = (
                state_q[_field3_idx(i, k, ixht, pcols, pver)] * wimp
                + wsrc * qrel
            )
            ptend_q[_field3_idx(i, k, ixht, psetcols, pver)] = (
                xhorz - state_q[_field3_idx(i, k, ixht, pcols, pver)]
            ) / dt

            qrel = qrel_vert[_vec_idx(k)]
            xvert = (
                wimp * state_q[_field3_idx(i, k, ixvt, pcols, pver)]
                + wsrc * qrel
            )
            ptend_q[_field3_idx(i, k, ixvt, psetcols, pver)] = (
                xvert - state_q[_field3_idx(i, k, ixvt, pcols, pver)]
            ) / dt

    for i in range(1, ncol + 1):
        cflx[_flux_idx(i, ixaoa1, pcols)] = 1.0e-6

        if (
            landfrac[_col_idx(i)] == 1.0
            and state_lat[_col_idx(i)] > 0.35
        ):
            cflx[_flux_idx(i, ixaoa2, pcols)] = (
                1.0e-6 + 1.0e-6 * 0.0434 * float(nstep) * dt / (86400.0 * 365.0)
            )
        else:
            cflx[_flux_idx(i, ixaoa2, pcols)] = 0.0

        cflx[_flux_idx(i, ixht, pcols)] = 0.0
        cflx[_flux_idx(i, ixvt, pcols)] = 0.0
