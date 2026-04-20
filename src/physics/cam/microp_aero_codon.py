@inline
def _idx2(i: int, k: int, ld1: int) -> int:
    return (i - 1) + (k - 1) * ld1


@inline
def _idx3(i: int, k: int, m: int, ld1: int, ld2: int) -> int:
    return (i - 1) + (k - 1) * ld1 + (m - 1) * ld1 * ld2


@export
def microp_aero_init_fields_codon(
    ncol: int,
    pcols: int,
    pver: int,
    rn_dst1: float,
    rn_dst2: float,
    rn_dst3: float,
    rn_dst4: float,
    npccn_p: cobj,
    nacon_p: cobj,
    rndst_p: cobj,
):
    npccn = Ptr[float](npccn_p)
    nacon = Ptr[float](nacon_p)
    rndst = Ptr[float](rndst_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            npccn[_idx2(i, k, pcols)] = 0.0
            for m in range(1, 5):
                nacon[_idx3(i, k, m, pcols, pver)] = 0.0
            rndst[_idx3(i, k, 1, pcols, pver)] = rn_dst1
            rndst[_idx3(i, k, 2, pcols, pver)] = rn_dst2
            rndst[_idx3(i, k, 3, pcols, pver)] = rn_dst3
            rndst[_idx3(i, k, 4, pcols, pver)] = rn_dst4
