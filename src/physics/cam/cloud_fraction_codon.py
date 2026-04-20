from math import log


@inline
def _field2_idx(i: int, k: int, ld1: int) -> int:
    """t/fice/fsnow declared as (ld1, pver)"""
    return (i - 1) + (k - 1) * ld1


@inline
def _field1_idx(i: int) -> int:
    return i - 1


@inline
def _relhum_min(
    press: float,
    lat: float,
    rhminh: float,
    rhminp: float,
    cldfrc_rhminp_botmb: float,
    unset_r8: float,
    pi: float,
) -> float:
    rh = rhminh
    if rhminp == unset_r8:
        return rh

    if press < cldfrc_rhminp_botmb * 1.0e2:
        if abs(lat * 180.0 / pi) > 60.0:
            rh = rhminp

    return rh


@export
def cldfrc_layer_rh_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    cldfrc_freeze_dry_i: int,
    premib: float,
    premit: float,
    rhminl: float,
    rhminl_adj_land: float,
    rhminh: float,
    rhminp: float,
    cldfrc_rhminp_botmb: float,
    unset_r8: float,
    pi: float,
    landfrac_p: cobj,
    snowh_p: cobj,
    clat_p: cobj,
    pmid_p: cobj,
    pref_mid_p: cobj,
    q_p: cobj,
    rh_p: cobj,
    rhcloud_p: cobj,
    rhu00_p: cobj,
):
    landfrac = Ptr[float](landfrac_p)
    snowh = Ptr[float](snowh_p)
    clat = Ptr[float](clat_p)
    pmid = Ptr[float](pmid_p)
    pref_mid = Ptr[float](pref_mid_p)
    q = Ptr[float](q_p)
    rh = Ptr[float](rh_p)
    rhcloud = Ptr[float](rhcloud_p)
    rhu00 = Ptr[float](rhu00_p)

    for k in range(top_lev + 1, pver + 1):
        pref_mid_k = pref_mid[_field1_idx(k)]
        for i in range(1, ncol + 1):
            idx = _field2_idx(i, k, pcols)
            i1 = _field1_idx(i)
            land_i = int(landfrac[i1] + 0.5) == 1
            pmid_ik = pmid[idx]

            if pmid_ik >= premib:
                if land_i and (snowh[i1] <= 0.000001):
                    rhlim = rhminl - rhminl_adj_land
                else:
                    rhlim = rhminl

                rhdif = (rh[idx] - rhlim) / (1.0 - rhlim)
                if rhdif > 0.0:
                    rhcloud_val = rhdif * rhdif
                else:
                    rhcloud_val = 0.0

                if rhcloud_val > 0.999:
                    rhcloud[idx] = 0.999
                else:
                    rhcloud[idx] = rhcloud_val

                if cldfrc_freeze_dry_i != 0:
                    qscale = q[idx] / 0.0030
                    if qscale > 1.0:
                        qscale = 1.0
                    elif qscale < 0.15:
                        qscale = 0.15
                    rhcloud[idx] = rhcloud[idx] * qscale

            elif pmid_ik < premit:
                rhlim = _relhum_min(pref_mid_k, clat[i1], rhminh, rhminp, cldfrc_rhminp_botmb, unset_r8, pi)
                rhdif = (rh[idx] - rhlim) / (1.0 - rhlim)
                if rhdif > 0.0:
                    rhcloud_val = rhdif * rhdif
                else:
                    rhcloud_val = 0.0

                if rhcloud_val > 0.999:
                    rhcloud[idx] = 0.999
                else:
                    rhcloud[idx] = rhcloud_val

            else:
                if pmid_ik > premit:
                    pmid_or_premit = pmid_ik
                else:
                    pmid_or_premit = premit
                rhwght = (premib - pmid_or_premit) / (premib - premit)

                rhlim_high = _relhum_min(pref_mid_k, clat[i1], rhminh, rhminp, cldfrc_rhminp_botmb, unset_r8, pi)
                if land_i and (snowh[i1] <= 0.000001):
                    rhlim = rhlim_high * rhwght + (rhminl - rhminl_adj_land) * (1.0 - rhwght)
                else:
                    rhlim = rhlim_high * rhwght + rhminl * (1.0 - rhwght)

                rhdif = (rh[idx] - rhlim) / (1.0 - rhlim)
                if rhdif > 0.0:
                    rhcloud_val = rhdif * rhdif
                else:
                    rhcloud_val = 0.0

                if rhcloud_val > 0.999:
                    rhcloud[idx] = 0.999
                else:
                    rhcloud[idx] = rhcloud_val

            rhu00[idx] = rhlim


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
