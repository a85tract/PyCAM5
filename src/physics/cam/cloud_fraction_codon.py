@inline
def _field2_idx(i: int, k: int, ld1: int) -> int:
    """t/fice/fsnow declared as (ld1, pver)"""
    return (i - 1) + (k - 1) * ld1


@export
def cldfrc_fice_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    t_p: cobj,
    fice_p: cobj,
    fsnow_p: cobj,
    tmax_fice: float,
    tmin_fice: float,
    tmax_fsnow: float,
    tmin_fsnow: float,
):
    t = Ptr[float](t_p)
    fice = Ptr[float](fice_p)
    fsnow = Ptr[float](fsnow_p)

    for k in range(1, top_lev):
        for i in range(1, pcols + 1):
            fice[_field2_idx(i, k, pcols)] = 0.0
            fsnow[_field2_idx(i, k, pcols)] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            t_ik = t[_field2_idx(i, k, pcols)]

            if t_ik > tmax_fice:
                fice[_field2_idx(i, k, pcols)] = 0.0
            elif t_ik < tmin_fice:
                fice[_field2_idx(i, k, pcols)] = 1.0
            else:
                fice[_field2_idx(i, k, pcols)] = (tmax_fice - t_ik) / (tmax_fice - tmin_fice)

            if t_ik > tmax_fsnow:
                fsnow[_field2_idx(i, k, pcols)] = 0.0
            elif t_ik < tmin_fsnow:
                fsnow[_field2_idx(i, k, pcols)] = 1.0
            else:
                fsnow[_field2_idx(i, k, pcols)] = (tmax_fsnow - t_ik) / (tmax_fsnow - tmin_fsnow)
