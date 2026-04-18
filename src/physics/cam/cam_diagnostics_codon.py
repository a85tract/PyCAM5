from math import sqrt


@inline
def _idx(i: int) -> int:
    return i - 1


@inline
def _idx2(i: int, k: int, ld1: int) -> int:
    return (i - 1) + (k - 1) * ld1


@export
def diag_surf_codon(
    ncol: int,
    pcols: int,
    end_day: int,
    qref_p: cobj,
    rhref_p: cobj,
    tref_p: cobj,
    trefmxav_p: cobj,
    trefmnav_p: cobj,
    trefmx_day_p: cobj,
    trefmn_day_p: cobj,
):
    qref = Ptr[float](qref_p)
    rhref = Ptr[float](rhref_p)
    tref = Ptr[float](tref_p)
    trefmxav = Ptr[float](trefmxav_p)
    trefmnav = Ptr[float](trefmnav_p)
    trefmx_day = Ptr[float](trefmx_day_p)
    trefmn_day = Ptr[float](trefmn_day_p)

    for i in range(1, ncol + 1):
        rhref[_idx(i)] = qref[_idx(i)] / rhref[_idx(i)] * 100.0
        trefmxav[_idx(i)] = max(tref[_idx(i)], trefmxav[_idx(i)])
        trefmnav[_idx(i)] = min(tref[_idx(i)], trefmnav[_idx(i)])

    if end_day != 0:
        for i in range(1, ncol + 1):
            trefmx_day[_idx(i)] = trefmxav[_idx(i)]
            trefmn_day[_idx(i)] = trefmnav[_idx(i)]
            trefmxav[_idx(i)] = -1.0e36
            trefmnav[_idx(i)] = 1.0e36


@export
def diag_physvar_ic_codon():
    return


@export
def diag_phys_writeout_z3_codon(
    ncol: int,
    pcols: int,
    pver: int,
    rga: float,
    zm_p: cobj,
    phis_p: cobj,
    z3_p: cobj,
):
    zm = Ptr[float](zm_p)
    phis = Ptr[float](phis_p)
    z3 = Ptr[float](z3_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            z3[_idx2(i, k, pcols)] = zm[_idx2(i, k, pcols)] + phis[_idx(i)] * rga


@export
def diag_phys_writeout_mul_codon(
    ncol: int,
    pcols: int,
    pver: int,
    a_p: cobj,
    b_p: cobj,
    out_p: cobj,
):
    a = Ptr[float](a_p)
    b = Ptr[float](b_p)
    out = Ptr[float](out_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            out[_idx2(i, k, pcols)] = a[_idx2(i, k, pcols)] * b[_idx2(i, k, pcols)]


@export
def diag_phys_writeout_mul_scalar_codon(
    ncol: int,
    pcols: int,
    pver: int,
    scale: float,
    a_p: cobj,
    b_p: cobj,
    out_p: cobj,
):
    a = Ptr[float](a_p)
    b = Ptr[float](b_p)
    out = Ptr[float](out_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            out[_idx2(i, k, pcols)] = (
                a[_idx2(i, k, pcols)] * b[_idx2(i, k, pcols)] * scale
            )


@export
def diag_phys_writeout_square_codon(
    ncol: int,
    pcols: int,
    pver: int,
    a_p: cobj,
    out_p: cobj,
):
    a = Ptr[float](a_p)
    out = Ptr[float](out_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            val = a[_idx2(i, k, pcols)]
            out[_idx2(i, k, pcols)] = val * val


@export
def diag_phys_writeout_wspeed_codon(
    ncol: int,
    pcols: int,
    pver: int,
    u_p: cobj,
    v_p: cobj,
    out_p: cobj,
):
    u = Ptr[float](u_p)
    v = Ptr[float](v_p)
    out = Ptr[float](out_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            uval = u[_idx2(i, k, pcols)]
            vval = v[_idx2(i, k, pcols)]
            out[_idx2(i, k, pcols)] = sqrt(uval * uval + vval * vval)


@export
def diag_phys_writeout_mass_and_tmq_codon(
    ncol: int,
    pcols: int,
    pver: int,
    rga: float,
    q_p: cobj,
    pdel_p: cobj,
    mq_p: cobj,
    tmq_p: cobj,
):
    q = Ptr[float](q_p)
    pdel = Ptr[float](pdel_p)
    mq = Ptr[float](mq_p)
    tmq = Ptr[float](tmq_p)

    for i in range(1, ncol + 1):
        total = 0.0
        for k in range(1, pver + 1):
            val = q[_idx2(i, k, pcols)] * pdel[_idx2(i, k, pcols)] * rga
            mq[_idx2(i, k, pcols)] = val
            total += val
        tmq[_idx(i)] = total


@export
def diag_phys_writeout_atmeint_codon(
    ncol: int,
    pcols: int,
    pver: int,
    cpair: float,
    latvap: float,
    gravit: float,
    t_p: cobj,
    q_p: cobj,
    u_p: cobj,
    v_p: cobj,
    pdel_p: cobj,
    phis_p: cobj,
    atmeint_p: cobj,
):
    t = Ptr[float](t_p)
    q = Ptr[float](q_p)
    u = Ptr[float](u_p)
    v = Ptr[float](v_p)
    pdel = Ptr[float](pdel_p)
    phis = Ptr[float](phis_p)
    atmeint = Ptr[float](atmeint_p)

    for i in range(1, ncol + 1):
        total = 0.0
        for k in range(1, pver + 1):
            uval = u[_idx2(i, k, pcols)]
            vval = v[_idx2(i, k, pcols)]
            total += (
                cpair * t[_idx2(i, k, pcols)]
                + phis[_idx(i)]
                + latvap * q[_idx2(i, k, pcols)]
                + 0.5 * (uval * uval + vval * vval)
            ) * (pdel[_idx2(i, k, pcols)] / gravit)
        atmeint[_idx(i)] = total
