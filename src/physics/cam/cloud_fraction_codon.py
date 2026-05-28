from math import acos, cos, log


@export
def cldfrc_register_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def cldfrc_codon(stage: int) -> int:
    return stage


@inline
def _ascii_ptr_to_str(n: int, ptr_p: cobj) -> str:
    ptr = Ptr[int](ptr_p)
    out = ""
    for i in range(n):
        out += chr(ptr[i])
    return out.strip()


@inline
def _strip_comment(line: str) -> str:
    pos = line.find("!")
    if pos >= 0:
        return line[:pos]
    return line


@inline
def _parse_fortran_float(value: str) -> float:
    return float(value.strip().replace("D", "E").replace("d", "e"))


@inline
def _parse_fortran_int(value: str) -> int:
    return int(value.strip())


@inline
def _parse_fortran_bool(value: str) -> int:
    text = value.strip().lower()
    if text == ".true." or text == "t" or text == "true":
        return 1
    return 0


@export
def cldfrc_readnl_codon(
    path_len: int,
    path_ascii_p: cobj,
    reals_p: cobj,
    ints_p: cobj,
    logicals_p: cobj,
) -> int:
    path = _ascii_ptr_to_str(path_len, path_ascii_p)
    reals = Ptr[float](reals_p)
    ints = Ptr[i32](ints_p)
    logicals = Ptr[i32](logicals_p)

    f = open(path, "r")
    text = f.read()
    f.close()

    in_group = False
    found_group = False
    assignments = ""

    for raw_line in text.split("\n"):
        line = _strip_comment(raw_line).strip()
        lowered = line.lower()
        if not in_group:
            if lowered.startswith("&cldfrc_nl"):
                in_group = True
                found_group = True
                rest = line[len("&cldfrc_nl") :]
                if rest:
                    assignments += rest + ","
            continue

        slash = line.find("/")
        if slash >= 0:
            assignments += line[:slash] + ","
            break
        assignments += line + ","

    if not found_group:
        return 0

    for item in assignments.split(","):
        if "=" not in item:
            continue
        parts = item.split("=", 1)
        key = parts[0].strip().lower()
        value = parts[1].strip()
        if key == "cldfrc_rhminl":
            reals[0] = _parse_fortran_float(value)
        elif key == "cldfrc_rhminl_adj_land":
            reals[1] = _parse_fortran_float(value)
        elif key == "cldfrc_rhminh":
            reals[2] = _parse_fortran_float(value)
        elif key == "cldfrc_rhminp":
            reals[3] = _parse_fortran_float(value)
        elif key == "cldfrc_rhminp_botmb":
            reals[4] = _parse_fortran_float(value)
        elif key == "cldfrc_sh1":
            reals[5] = _parse_fortran_float(value)
        elif key == "cldfrc_sh2":
            reals[6] = _parse_fortran_float(value)
        elif key == "cldfrc_dp1":
            reals[7] = _parse_fortran_float(value)
        elif key == "cldfrc_dp2":
            reals[8] = _parse_fortran_float(value)
        elif key == "cldfrc_premit":
            reals[9] = _parse_fortran_float(value)
        elif key == "cldfrc_premib":
            reals[10] = _parse_fortran_float(value)
        elif key == "cldfrc_icecrit":
            reals[11] = _parse_fortran_float(value)
        elif key == "cldfrc_iceopt":
            ints[0] = i32(_parse_fortran_int(value))
        elif key == "cldfrc_freeze_dry":
            logicals[0] = i32(_parse_fortran_bool(value))
        elif key == "cldfrc_ice":
            logicals[1] = i32(_parse_fortran_bool(value))

    return 0


@export
def cldfrc_init_codon(
    pver: int,
    macrop_rk_i: int,
    eddy_diag_tke_i: int,
    shallow_uw_i: int,
    trop_cloud_top_lev: int,
    pref_mid_p: cobj,
    top_lev_p: cobj,
    inversion_cld_off_p: cobj,
    k700_p: cobj,
) -> int:
    pref_mid = Ptr[float](pref_mid_p)
    top_lev_out = Ptr[i32](top_lev_p)
    inversion_out = Ptr[i32](inversion_cld_off_p)
    k700_out = Ptr[i32](k700_p)

    top = 1
    if macrop_rk_i == 0:
        top = trop_cloud_top_lev
    top_lev_out[0] = i32(top)

    if eddy_diag_tke_i != 0 or shallow_uw_i != 0:
        inversion_out[0] = i32(1)
    else:
        inversion_out[0] = i32(0)

    if pref_mid[top - 1] > 7.0e4:
        return 1

    best_section_index = 1
    best_abs = abs(pref_mid[top - 1] - 7.0e4)
    section_index = 1
    for k in range(top + 1, pver + 1):
        section_index += 1
        diff = abs(pref_mid[k - 1] - 7.0e4)
        if diff < best_abs:
            best_abs = diff
            best_section_index = section_index

    k700_out[0] = i32(best_section_index)
    return 0


@inline
def _field2_idx(i: int, k: int, ld1: int) -> int:
    """t/fice/fsnow declared as (ld1, pver)"""
    return (i - 1) + (k - 1) * ld1


@inline
def _field1_idx(i: int) -> int:
    return i - 1


@export
def cldfrc_getparams_codon(
    flags: int,
    rhminl_in: float,
    rhminl_adj_land_in: float,
    rhminh_in: float,
    rhminp_in: float,
    premit_in: float,
    premib_in: float,
    iceopt_in: int,
    icecrit_in: float,
    rhminl_p: cobj,
    rhminl_adj_land_p: cobj,
    rhminh_p: cobj,
    rhminp_p: cobj,
    premit_p: cobj,
    premib_p: cobj,
    iceopt_p: cobj,
    icecrit_p: cobj,
):
    if flags & 1:
        Ptr[float](rhminl_p)[0] = rhminl_in
    if flags & 2:
        Ptr[float](rhminl_adj_land_p)[0] = rhminl_adj_land_in
    if flags & 4:
        Ptr[float](rhminh_p)[0] = rhminh_in
    if flags & 8:
        Ptr[float](rhminp_p)[0] = rhminp_in
    if flags & 16:
        Ptr[float](premit_p)[0] = premit_in
    if flags & 32:
        Ptr[float](premib_p)[0] = premib_in
    if flags & 64:
        Ptr[i32](iceopt_p)[0] = i32(iceopt_in)
    if flags & 128:
        Ptr[float](icecrit_p)[0] = icecrit_in


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
def relhum_min_codon(
    press: float,
    lat: float,
    rhminh: float,
    rhminp: float,
    cldfrc_rhminp_botmb: float,
    unset_r8: float,
    pi: float,
) -> float:
    return _relhum_min(press, lat, rhminh, rhminp, cldfrc_rhminp_botmb, unset_r8, pi)


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
def cldfrc_ice_wilson_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    icecrit: float,
    one_sixth: float,
    two_thirds: float,
    two_pow_three_halves: float,
    phi_offset: float,
    cldice_p: cobj,
    qs_p: cobj,
    rhcloud_p: cobj,
    icecldf_p: cobj,
    liqcldf_p: cobj,
    cloud_p: cobj,
):
    cldice = Ptr[float](cldice_p)
    qs = Ptr[float](qs_p)
    rhcloud = Ptr[float](rhcloud_p)
    icecldf = Ptr[float](icecldf_p)
    liqcldf = Ptr[float](liqcldf_p)
    cloud = Ptr[float](cloud_p)

    for k in range(top_lev + 1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _field2_idx(i, k, pcols)
            ncf = cldice[idx] / ((1.0 - icecrit) * qs[idx])

            if ncf <= 0.0:
                icecldf[idx] = 0.0
            elif ncf <= one_sixth:
                icecldf[idx] = 0.5 * (6.0 * ncf) ** two_thirds
            elif ncf < 1.0:
                phi = (acos(3.0 * (1.0 - ncf) / two_pow_three_halves) + phi_offset) / 3.0
                icecldf[idx] = 1.0 - 4.0 * cos(phi) * cos(phi)
            else:
                icecldf[idx] = 1.0

            if icecldf[idx] < 0.0:
                icecldf[idx] = 0.0
            elif icecldf[idx] > 1.0:
                icecldf[idx] = 1.0

            liqcldf[idx] = (1.0 - icecldf[idx]) * rhcloud[idx]
            cloud[idx] = liqcldf[idx] + icecldf[idx]


@export
def cldfrc_marine_stratus_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    k700: int,
    inversion_cld_off: int,
    premib: float,
    theta_p: cobj,
    thetas_p: cobj,
    pmid_p: cobj,
    ps_p: cobj,
    ocnfrac_p: cobj,
    rh_p: cobj,
    rpdeli_p: cobj,
    cldst_p: cobj,
    dthdpmn_p: cobj,
    kdthdp_p: cobj,
):
    theta = Ptr[float](theta_p)
    thetas = Ptr[float](thetas_p)
    pmid = Ptr[float](pmid_p)
    ps = Ptr[float](ps_p)
    ocnfrac = Ptr[float](ocnfrac_p)
    rh = Ptr[float](rh_p)
    rpdeli = Ptr[float](rpdeli_p)
    cldst = Ptr[float](cldst_p)
    dthdpmn = Ptr[float](dthdpmn_p)
    kdthdp = Ptr[int](kdthdp_p)

    if inversion_cld_off != 0:
        return

    for i in range(1, ncol + 1):
        dthdpmn[i - 1] = -0.125
        kdthdp[i - 1] = 0

    for k in range(top_lev + 1, pver + 1):
        for i in range(1, ncol + 1):
            i0 = i - 1
            idx = _field2_idx(i, k, pcols)
            if pmid[idx] >= premib and ocnfrac[i0] > 0.01:
                dthdp = 100.0 * (theta[idx] - theta[_field2_idx(i, k - 1, pcols)]) * rpdeli[_field2_idx(i, k - 1, pcols)]
                if dthdp < dthdpmn[i0]:
                    dthdpmn[i0] = dthdp
                    kdthdp[i0] = k

    for i in range(1, ncol + 1):
        i0 = i - 1
        if kdthdp[i0] == 0 and ocnfrac[i0] > 0.01:
            dthdp = 100.0 * (thetas[i0] - theta[_field2_idx(i, pver, pcols)]) / (ps[i0] - pmid[_field2_idx(i, pver, pcols)])
            if dthdp < dthdpmn[i0]:
                dthdpmn[i0] = dthdp
                kdthdp[i0] = pver

    for i in range(1, ncol + 1):
        i0 = i - 1
        k = kdthdp[i0]
        if k != 0:
            kp1 = k + 1
            if kp1 > pver:
                kp1 = pver
            strat = ocnfrac[i0] * ((theta[_field2_idx(i, k700, pcols)] - thetas[i0]) * 0.057 - 0.5573)
            if strat < 0.0:
                strat = 0.0
            elif strat > 1.0:
                strat = 1.0
            max_rh = rh[_field2_idx(i, k, pcols)]
            rh_kp1 = rh[_field2_idx(i, kp1, pcols)]
            if rh_kp1 > max_rh:
                max_rh = rh_kp1
            if max_rh < strat:
                cldst[_field2_idx(i, k, pcols)] = max_rh
            else:
                cldst[_field2_idx(i, k, pcols)] = strat


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


@export
def cldfrc_batch_layer_rh_codon(
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
    cldfrc_layer_rh_codon(
        ncol,
        pcols,
        pver,
        top_lev,
        cldfrc_freeze_dry_i,
        premib,
        premit,
        rhminl,
        rhminl_adj_land,
        rhminh,
        rhminp,
        cldfrc_rhminp_botmb,
        unset_r8,
        pi,
        landfrac_p,
        snowh_p,
        clat_p,
        pmid_p,
        pref_mid_p,
        q_p,
        rh_p,
        rhcloud_p,
        rhu00_p,
    )


@export
def cldfrc_batch_ice_wilson_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    icecrit: float,
    one_sixth: float,
    two_thirds: float,
    two_pow_three_halves: float,
    phi_offset: float,
    cldice_p: cobj,
    qs_p: cobj,
    rhcloud_p: cobj,
    icecldf_p: cobj,
    liqcldf_p: cobj,
    cloud_p: cobj,
):
    cldfrc_ice_wilson_codon(
        ncol,
        pcols,
        pver,
        top_lev,
        icecrit,
        one_sixth,
        two_thirds,
        two_pow_three_halves,
        phi_offset,
        cldice_p,
        qs_p,
        rhcloud_p,
        icecldf_p,
        liqcldf_p,
        cloud_p,
    )


@export
def cldfrc_batch_state_init_codon(
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
    cldfrc_state_init_codon(
        ncol,
        pcols,
        pver,
        top_lev,
        dindex,
        rhpert,
        q_p,
        qs_p,
        relhum_p,
        rh_p,
        cloud_p,
        icecldf_p,
        liqcldf_p,
        rhcloud_p,
        cldst_p,
        concld_p,
    )


@export
def cldfrc_batch_total_cloud_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    rhcloud_p: cobj,
    cldst_p: cobj,
    concld_p: cobj,
    cloud_p: cobj,
):
    cldfrc_total_cloud_codon(ncol, pcols, pver, top_lev, rhcloud_p, cldst_p, concld_p, cloud_p)


@export
def cldfrc_batch_convective_cover_codon(
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
    cldfrc_convective_cover_codon(
        ncol,
        pcols,
        pver,
        top_lev,
        use_shfrc_i,
        sh1,
        sh2,
        dp1,
        dp2,
        shfrc_p,
        cmfmc_p,
        cmfmc2_p,
        shallowcu_p,
        deepcu_p,
        concld_p,
        rh_p,
    )


@export
def cldfrc_batch_fice_codon(
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
    cldfrc_fice_codon(
        ncol,
        pcols,
        pver,
        top_lev,
        t_p,
        fice_p,
        fsnow_p,
        tmax_fice,
        tmin_fice,
        tmax_fsnow,
        tmin_fsnow,
    )


@export
def cldfrc_batch_dispatch_codon(
    stage: int,
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    flag: int,
    scalar1: float,
    scalar2: float,
    scalar3: float,
    scalar4: float,
    scalar5: float,
    scalar6: float,
    scalar7: float,
    scalar8: float,
    scalar9: float,
    p1: cobj,
    p2: cobj,
    p3: cobj,
    p4: cobj,
    p5: cobj,
    p6: cobj,
    p7: cobj,
    p8: cobj,
    p9: cobj,
    p10: cobj,
):
    if stage == 1:
        cldfrc_batch_state_init_codon(
            ncol,
            pcols,
            pver,
            top_lev,
            flag,
            scalar1,
            p1,
            p2,
            p3,
            p4,
            p5,
            p6,
            p7,
            p8,
            p9,
            p10,
        )
    elif stage == 2:
        cldfrc_batch_layer_rh_codon(
            ncol,
            pcols,
            pver,
            top_lev,
            flag,
            scalar1,
            scalar2,
            scalar3,
            scalar4,
            scalar5,
            scalar6,
            scalar7,
            scalar8,
            scalar9,
            p1,
            p2,
            p3,
            p4,
            p5,
            p6,
            p7,
            p8,
            p9,
        )
    elif stage == 3:
        cldfrc_batch_ice_wilson_codon(
            ncol,
            pcols,
            pver,
            top_lev,
            scalar1,
            scalar2,
            scalar3,
            scalar4,
            scalar5,
            p1,
            p2,
            p3,
            p4,
            p5,
            p6,
        )
    elif stage == 4:
        cldfrc_batch_convective_cover_codon(
            ncol,
            pcols,
            pver,
            top_lev,
            flag,
            scalar1,
            scalar2,
            scalar3,
            scalar4,
            p1,
            p2,
            p3,
            p4,
            p5,
            p6,
            p7,
        )
    elif stage == 5:
        cldfrc_batch_total_cloud_codon(ncol, pcols, pver, top_lev, p1, p2, p3, p4)
    elif stage == 6:
        cldfrc_batch_fice_codon(
            ncol,
            pcols,
            pver,
            top_lev,
            p1,
            p2,
            p3,
            scalar1,
            scalar2,
            scalar3,
            scalar4,
        )
