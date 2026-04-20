from math import log


@inline
def _field2_idx(i: int, k: int, ld1: int) -> int:
    """t/fice/fsnow declared as (ld1, pver)"""
    return (i - 1) + (k - 1) * ld1


@export
def cldfrc_convective_cover_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    use_shfrc_i: int,
    sh1: float,
    sh2: float,
    dp1: float,
    dp2: float,
    shfrc_p: cobj,
    cmfmc_p: cobj,
    cmfmc2_p: cobj,
    shallowcu_p: cobj,
    deepcu_p: cobj,
    concld_p: cobj,
    rh_p: cobj,
):
    shfrc = Ptr[float](shfrc_p)
    cmfmc = Ptr[float](cmfmc_p)
    cmfmc2 = Ptr[float](cmfmc2_p)
    shallowcu = Ptr[float](shallowcu_p)
    deepcu = Ptr[float](deepcu_p)
    concld = Ptr[float](concld_p)
    rh = Ptr[float](rh_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx = _field2_idx(i, k, pcols)
            idx_kp1 = _field2_idx(i, k + 1, pcols)

            if use_shfrc_i == 0:
                shallow_tmp = sh1 * log(1.0 + sh2 * cmfmc2[idx_kp1])
                if shallow_tmp < 0.0:
                    shallowcu[idx] = 0.0
                elif shallow_tmp > 0.30:
                    shallowcu[idx] = 0.30
                else:
                    shallowcu[idx] = shallow_tmp
            else:
                shallowcu[idx] = shfrc[idx]

            deep_tmp = dp1 * log(1.0 + dp2 * (cmfmc[idx_kp1] - cmfmc2[idx_kp1]))
            if deep_tmp < 0.0:
                deepcu[idx] = 0.0
            elif deep_tmp > 0.60:
                deepcu[idx] = 0.60
            else:
                deepcu[idx] = deep_tmp

            concld_tmp = shallowcu[idx] + deepcu[idx]
            if concld_tmp > 0.80:
                concld[idx] = 0.80
            else:
                concld[idx] = concld_tmp

            rh[idx] = (rh[idx] - concld[idx]) / (1.0 - concld[idx])


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
