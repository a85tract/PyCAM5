@inline
def _field_idx(icol: int, klev: int, ld1: int) -> int:
    """qrl/qrs declared as (pcols, pver); ptend%s declared as (psetcols, pver)"""
    return (icol - 1) + (klev - 1) * ld1


@inline
def _col_idx(icol: int) -> int:
    """fsns/fsnt/flns/flnt/net_flx declared as (pcols)"""
    return icol - 1


@export
def radheat_timestep_init_codon():
    return


@export
def radheat_tend_codon(
    ncol: int,
    pcols: int,
    pver: int,
    psetcols: int,
    qrl_p: cobj,
    qrs_p: cobj,
    ptend_s_p: cobj,
    fsns_p: cobj,
    fsnt_p: cobj,
    flns_p: cobj,
    flnt_p: cobj,
    net_flx_p: cobj,
):
    qrl = Ptr[float](qrl_p)
    qrs = Ptr[float](qrs_p)
    ptend_s = Ptr[float](ptend_s_p)
    fsns = Ptr[float](fsns_p)
    fsnt = Ptr[float](fsnt_p)
    flns = Ptr[float](flns_p)
    flnt = Ptr[float](flnt_p)
    net_flx = Ptr[float](net_flx_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            ptend_s[_field_idx(i, k, psetcols)] = (
                qrs[_field_idx(i, k, pcols)] + qrl[_field_idx(i, k, pcols)]
            )

    for i in range(1, ncol + 1):
        net_flx[_col_idx(i)] = (
            fsnt[_col_idx(i)]
            - fsns[_col_idx(i)]
            - flnt[_col_idx(i)]
            + flns[_col_idx(i)]
        )
