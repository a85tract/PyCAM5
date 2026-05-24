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
def tracers_implements_cnst_codon(flag: int, name_len: int, name_ascii_p: cobj, ncnst: int) -> int:
    if flag == 0:
        return 0
    name_ascii = Ptr[int](name_ascii_p)
    m = 1
    while m <= ncnst:
        if _tracers_name_match(name_len, name_ascii, m) != 0:
            return 1
        m += 1
    return 0


@inline
def _tracers_name_match(name_len: int, name_ascii: Ptr[int], m: int) -> int:
    nbase = (m - 1) % 5 + 1
    ncopy = (m - 1) // 5
    if ncopy == 0:
        if nbase == 1:
            return _name_eq5(name_len, name_ascii, 84, 84, 95, 76, 87)
        if nbase == 2:
            return _name_eq5(name_len, name_ascii, 84, 84, 95, 77, 68)
        if nbase == 3:
            return _name_eq5(name_len, name_ascii, 84, 84, 95, 72, 73)
        if nbase == 4:
            return _name_eq5(name_len, name_ascii, 84, 84, 82, 77, 68)
        if nbase == 5:
            return _name_eq5(name_len, name_ascii, 84, 84, 95, 85, 78)
    if ncopy >= 1 and ncopy <= 9:
        digit = 48 + ncopy
        if nbase == 1:
            return _name_eq6(name_len, name_ascii, 84, 84, 95, 76, 87, digit)
        if nbase == 2:
            return _name_eq6(name_len, name_ascii, 84, 84, 95, 77, 68, digit)
        if nbase == 3:
            return _name_eq6(name_len, name_ascii, 84, 84, 95, 72, 73, digit)
        if nbase == 4:
            return _name_eq6(name_len, name_ascii, 84, 84, 82, 77, 68, digit)
        if nbase == 5:
            return _name_eq6(name_len, name_ascii, 84, 84, 95, 85, 78, digit)
    if ncopy >= 10 and ncopy <= 99:
        d1 = 48 + ncopy // 10
        d2 = 48 + ncopy % 10
        if nbase == 1:
            return _name_eq7(name_len, name_ascii, 84, 84, 95, 76, 87, d1, d2)
        if nbase == 2:
            return _name_eq7(name_len, name_ascii, 84, 84, 95, 77, 68, d1, d2)
        if nbase == 3:
            return _name_eq7(name_len, name_ascii, 84, 84, 95, 72, 73, d1, d2)
        if nbase == 4:
            return _name_eq7(name_len, name_ascii, 84, 84, 82, 77, 68, d1, d2)
        if nbase == 5:
            return _name_eq7(name_len, name_ascii, 84, 84, 95, 85, 78, d1, d2)
    if ncopy >= 100 and ncopy <= 999:
        d1 = 48 + ncopy // 100
        d2 = 48 + (ncopy // 10) % 10
        d3 = 48 + ncopy % 10
        if nbase == 1:
            return _name_eq8(name_len, name_ascii, 84, 84, 95, 76, 87, d1, d2, d3)
        if nbase == 2:
            return _name_eq8(name_len, name_ascii, 84, 84, 95, 77, 68, d1, d2, d3)
        if nbase == 3:
            return _name_eq8(name_len, name_ascii, 84, 84, 95, 72, 73, d1, d2, d3)
        if nbase == 4:
            return _name_eq8(name_len, name_ascii, 84, 84, 82, 77, 68, d1, d2, d3)
        if nbase == 5:
            return _name_eq8(name_len, name_ascii, 84, 84, 95, 85, 78, d1, d2, d3)
    return 0


@inline
def _name_eq5(name_len: int, name_ascii: Ptr[int], c0: int, c1: int, c2: int, c3: int, c4: int) -> int:
    if name_len > 8:
        return _name_eq8(name_len, name_ascii, c0, c1, c2, c3, c4, 32, 32, 32)
    return _name_eq8(name_len, name_ascii, c0, c1, c2, c3, c4, 32, 32, 32)


@inline
def _name_eq6(name_len: int, name_ascii: Ptr[int], c0: int, c1: int, c2: int, c3: int, c4: int, c5: int) -> int:
    return _name_eq8(name_len, name_ascii, c0, c1, c2, c3, c4, c5, 32, 32)


@inline
def _name_eq7(name_len: int, name_ascii: Ptr[int], c0: int, c1: int, c2: int, c3: int, c4: int, c5: int, c6: int) -> int:
    return _name_eq8(name_len, name_ascii, c0, c1, c2, c3, c4, c5, c6, 32)


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
