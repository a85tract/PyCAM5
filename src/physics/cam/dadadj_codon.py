NITER = 15


@inline
def _col_lev_idx(icol: int, klev: int, pcols: int) -> int:
    """pmid/pdel/t/q declared as (pcols, pver)"""
    return (icol - 1) + (klev - 1) * pcols


@inline
def _col_int_idx(icol: int, kint: int, pcols: int) -> int:
    """pint declared as (pcols, pverp)"""
    return (icol - 1) + (kint - 1) * pcols


@inline
def _lev_idx(klev: int) -> int:
    """c1dad/c2dad/c3dad/c4dad declared as (pver)"""
    return klev - 1


@inline
def _col_idx(icol: int) -> int:
    """dodad declared as (pcols)"""
    return icol - 1


@export
def dadadj_codon(
    ncol: int,
    pcols: int,
    nlvdry: int,
    cappa: float,
    pmid_p: cobj,
    pint_p: cobj,
    pdel_p: cobj,
    t_p: cobj,
    q_p: cobj,
    c1dad_p: cobj,
    c2dad_p: cobj,
    c3dad_p: cobj,
    c4dad_p: cobj,
    dodad_p: cobj,
    status_p: cobj,
    zeps_fail_p: cobj,
    fail_i_p: cobj,
):
    pmid = Ptr[float](pmid_p)
    pint = Ptr[float](pint_p)
    pdel = Ptr[float](pdel_p)
    t = Ptr[float](t_p)
    q = Ptr[float](q_p)
    c1dad = Ptr[float](c1dad_p)
    c2dad = Ptr[float](c2dad_p)
    c3dad = Ptr[float](c3dad_p)
    c4dad = Ptr[float](c4dad_p)
    dodad = Ptr[int](dodad_p)
    status = Ptr[int](status_p)
    zeps_fail = Ptr[float](zeps_fail_p)
    fail_i = Ptr[int](fail_i_p)

    status[0] = 0
    zeps_fail[0] = 0.0
    fail_i[0] = 0

    zeps = 2.0e-5

    for i in range(1, ncol + 1):
        gammad = cappa * 0.5 * (
            t[_col_lev_idx(i, 2, pcols)] + t[_col_lev_idx(i, 1, pcols)]
        ) / pint[_col_int_idx(i, 2, pcols)]
        dtdp = (
            t[_col_lev_idx(i, 2, pcols)] - t[_col_lev_idx(i, 1, pcols)]
        ) / (
            pmid[_col_lev_idx(i, 2, pcols)] - pmid[_col_lev_idx(i, 1, pcols)]
        )
        if (dtdp + zeps) > gammad:
            dodad[_col_idx(i)] = 1
        else:
            dodad[_col_idx(i)] = 0

    for k in range(2, nlvdry + 1):
        for i in range(1, ncol + 1):
            gammad = cappa * 0.5 * (
                t[_col_lev_idx(i, k + 1, pcols)] + t[_col_lev_idx(i, k, pcols)]
            ) / pint[_col_int_idx(i, k + 1, pcols)]
            dtdp = (
                t[_col_lev_idx(i, k + 1, pcols)] - t[_col_lev_idx(i, k, pcols)]
            ) / (
                pmid[_col_lev_idx(i, k + 1, pcols)] - pmid[_col_lev_idx(i, k, pcols)]
            )
            if (dtdp + zeps) > gammad:
                dodad[_col_idx(i)] = 1

    for i in range(1, ncol + 1):
        if dodad[_col_idx(i)] != 0:
            zeps = 2.0e-5

            for k in range(1, nlvdry + 1):
                c1dad[_lev_idx(k)] = cappa * 0.5 * (
                    pmid[_col_lev_idx(i, k + 1, pcols)] - pmid[_col_lev_idx(i, k, pcols)]
                ) / pint[_col_int_idx(i, k + 1, pcols)]
                c2dad[_lev_idx(k)] = (1.0 - c1dad[_lev_idx(k)]) / (1.0 + c1dad[_lev_idx(k)])
                rdenom = 1.0 / (
                    pdel[_col_lev_idx(i, k, pcols)] * c2dad[_lev_idx(k)]
                    + pdel[_col_lev_idx(i, k + 1, pcols)]
                )
                c3dad[_lev_idx(k)] = rdenom * pdel[_col_lev_idx(i, k, pcols)]
                c4dad[_lev_idx(k)] = rdenom * pdel[_col_lev_idx(i, k + 1, pcols)]

            while True:
                ilconv = 1

                for jiter in range(1, NITER + 1):
                    ilconv = 1
                    for k in range(1, nlvdry + 1):
                        zepsdp = zeps * (
                            pmid[_col_lev_idx(i, k + 1, pcols)] - pmid[_col_lev_idx(i, k, pcols)]
                        )
                        zgamma = c1dad[_lev_idx(k)] * (
                            t[_col_lev_idx(i, k, pcols)] + t[_col_lev_idx(i, k + 1, pcols)]
                        )
                        if (
                            t[_col_lev_idx(i, k + 1, pcols)] - t[_col_lev_idx(i, k, pcols)]
                        ) >= (zgamma + zepsdp):
                            ilconv = 0
                            t[_col_lev_idx(i, k + 1, pcols)] = (
                                t[_col_lev_idx(i, k, pcols)] * c3dad[_lev_idx(k)]
                                + t[_col_lev_idx(i, k + 1, pcols)] * c4dad[_lev_idx(k)]
                            )
                            t[_col_lev_idx(i, k, pcols)] = (
                                c2dad[_lev_idx(k)] * t[_col_lev_idx(i, k + 1, pcols)]
                            )
                            qave = (
                                pdel[_col_lev_idx(i, k + 1, pcols)] * q[_col_lev_idx(i, k + 1, pcols)]
                                + pdel[_col_lev_idx(i, k, pcols)] * q[_col_lev_idx(i, k, pcols)]
                            ) / (
                                pdel[_col_lev_idx(i, k + 1, pcols)]
                                + pdel[_col_lev_idx(i, k, pcols)]
                            )
                            q[_col_lev_idx(i, k + 1, pcols)] = qave
                            q[_col_lev_idx(i, k, pcols)] = qave
                    if ilconv != 0:
                        break

                if ilconv != 0:
                    break

                zeps = zeps + zeps
                if zeps > 1.0e-4:
                    status[0] = 1
                    zeps_fail[0] = zeps
                    fail_i[0] = i
                    return
