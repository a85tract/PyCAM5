from math import log


@inline
def _field2_idx(i: int, k: int, ld1: int) -> int:
    """t/fice/fsnow declared as (ld1, pver)"""
    return (i - 1) + (k - 1) * ld1


@export
def cldfrc_state_init_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    dindex: int,
    rhpert: float,
    q_p: cobj,
    qs_p: cobj,
    relhum_p: cobj,
    rh_p: cobj,
    cloud_p: cobj,
    icecldf_p: cobj,
    liqcldf_p: cobj,
    rhcloud_p: cobj,
    cldst_p: cobj,
    concld_p: cobj,
):
    q = Ptr[float](q_p)
    qs = Ptr[float](qs_p)
    relhum = Ptr[float](relhum_p)
    rh = Ptr[float](rh_p)
    cloud = Ptr[float](cloud_p)
    icecldf = Ptr[float](icecldf_p)
    liqcldf = Ptr[float](liqcldf_p)
    rhcloud = Ptr[float](rhcloud_p)
    cldst = Ptr[float](cldst_p)
    concld = Ptr[float](concld_p)

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            idx = _field2_idx(i, k, pcols)
            cloud[idx] = 0.0
            icecldf[idx] = 0.0
            liqcldf[idx] = 0.0
            rhcloud[idx] = 0.0
            cldst[idx] = 0.0
            concld[idx] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx = _field2_idx(i, k, pcols)
            rh[idx] = q[idx] / qs[idx] * (1.0 + dindex * rhpert)
            relhum[idx] = rh[idx]


@export
def cldfrc_total_cloud_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    rhcloud_p: cobj,
    cldst_p: cobj,
    concld_p: cobj,
    cloud_p: cobj,
):
    rhcloud = Ptr[float](rhcloud_p)
    cldst = Ptr[float](cldst_p)
    concld = Ptr[float](concld_p)
    cloud = Ptr[float](cloud_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx = _field2_idx(i, k, pcols)
            if rhcloud[idx] > cldst[idx]:
                cloud[idx] = rhcloud[idx]
            else:
                cloud[idx] = cldst[idx]

            cloud_sum = cloud[idx] + concld[idx]
            if cloud_sum < 1.0:
                cloud[idx] = cloud_sum
            else:
                cloud[idx] = 1.0


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
