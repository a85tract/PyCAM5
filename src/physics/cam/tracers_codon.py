@inline
def _field3_idx(icol: int, klev: int, mconst: int, ld1: int, ld2: int) -> int:
    """ptend%q declared as (ld1, pver, pcnst)"""
    return (icol - 1) + (klev - 1) * ld1 + (mconst - 1) * ld1 * ld2


@export
def tracers_flag_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def tracers_implements_cnst_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def tracers_register_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def tracers_init_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@inline
def _flux_idx(icol: int, mconst: int, pcols: int) -> int:
    """cflx declared as (pcols, pcnst)"""
    return (icol - 1) + (mconst - 1) * pcols


@export
def tracers_timestep_init_codon():
    return


@export
def tracers_timestep_tend_codon(
    ncol: int,
    pcols: int,
    pver: int,
    psetcols: int,
    ixtrct: int,
    trac_ncnst: int,
    ptend_q_p: cobj,
    cflx_p: cobj,
):
    ptend_q = Ptr[float](ptend_q_p)
    cflx = Ptr[float](cflx_p)

    for m in range(1, trac_ncnst + 1):
        mconst = ixtrct + m - 1

        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                ptend_q[_field3_idx(i, k, mconst, psetcols, pver)] = 0.0

        for i in range(1, ncol + 1):
            cflx[_flux_idx(i, mconst, pcols)] = 0.0
