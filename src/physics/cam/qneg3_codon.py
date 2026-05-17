@inline
def _q_idx(icol: int, klev: int, mloc: int, ncold: int, lver: int) -> int:
    """q declared as (ncold, lver, nconst) with mloc = m - lconst_beg + 1"""
    return (icol - 1) + (klev - 1) * ncold + (mloc - 1) * ncold * lver


@inline
def _qmin_idx(mloc: int) -> int:
    """qmin declared as (nconst)"""
    return mloc - 1


@inline
def _indx_idx(nn: int, klev: int, ncol: int) -> int:
    """indx declared as (ncol, lver)"""
    return (nn - 1) + (klev - 1) * ncol


@inline
def _nval_idx(klev: int) -> int:
    """nval declared as (lver)"""
    return klev - 1


@inline
def _stat_idx(mloc: int) -> int:
    """per-constituent stats declared as (nconst)"""
    return mloc - 1


@export
def qneg3_codon(
    ncol: int,
    ncold: int,
    lver: int,
    nconst: int,
    qmin_p: cobj,
    q_p: cobj,
    indx_p: cobj,
    nval_p: cobj,
    nvals_p: cobj,
    worst_p: cobj,
    iw_p: cobj,
    kw_p: cobj,
):
    qmin = Ptr[float](qmin_p)
    q = Ptr[float](q_p)
    indx = Ptr[int](indx_p)
    nval = Ptr[int](nval_p)
    nvals = Ptr[int](nvals_p)
    worst = Ptr[float](worst_p)
    iw = Ptr[int](iw_p)
    kw = Ptr[int](kw_p)

    for mloc in range(1, nconst + 1):
        nvals[_stat_idx(mloc)] = 0
        worst[_stat_idx(mloc)] = 1.0e35
        iw[_stat_idx(mloc)] = -1
        kw[_stat_idx(mloc)] = -1

        for k in range(1, lver + 1):
            nn = 0
            for i in range(1, ncol + 1):
                if q[_q_idx(i, k, mloc, ncold, lver)] < qmin[_qmin_idx(mloc)]:
                    nn += 1
                    indx[_indx_idx(nn, k, ncol)] = i
            nval[_nval_idx(k)] = nn

        for k in range(1, lver + 1):
            if nval[_nval_idx(k)] > 0:
                nvals[_stat_idx(mloc)] += nval[_nval_idx(k)]
                iwtmp = -1

                for ii in range(1, nval[_nval_idx(k)] + 1):
                    i = indx[_indx_idx(ii, k, ncol)]
                    q_val = q[_q_idx(i, k, mloc, ncold, lver)]
                    if q_val < worst[_stat_idx(mloc)]:
                        worst[_stat_idx(mloc)] = q_val
                        iwtmp = ii

                if iwtmp != -1:
                    kw[_stat_idx(mloc)] = k
                    iw[_stat_idx(mloc)] = indx[_indx_idx(iwtmp, k, ncol)]

                for ii in range(1, nval[_nval_idx(k)] + 1):
                    i = indx[_indx_idx(ii, k, ncol)]
                    q[_q_idx(i, k, mloc, ncold, lver)] = qmin[_qmin_idx(mloc)]


@export
def qneg_batch_3_stage_dispatch_codon(
    ncol: int,
    ncold: int,
    lver: int,
    nconst: int,
    qmin_p: cobj,
    q_p: cobj,
    indx_p: cobj,
    nval_p: cobj,
    nvals_p: cobj,
    worst_p: cobj,
    iw_p: cobj,
    kw_p: cobj,
):
    qneg3_codon(
        ncol,
        ncold,
        lver,
        nconst,
        qmin_p,
        q_p,
        indx_p,
        nval_p,
        nvals_p,
        worst_p,
        iw_p,
        kw_p,
    )

@export
def qneg_batch_3_codon(
    ncol: int,
    ncold: int,
    lver: int,
    nconst: int,
    qmin_p: cobj,
    q_p: cobj,
    indx_p: cobj,
    nval_p: cobj,
    nvals_p: cobj,
    worst_p: cobj,
    iw_p: cobj,
    kw_p: cobj,
):
    qneg_batch_3_stage_dispatch_codon(
        ncol,
        ncold,
        lver,
        nconst,
        qmin_p,
        q_p,
        indx_p,
        nval_p,
        nvals_p,
        worst_p,
        iw_p,
        kw_p,
    )
