from chemistry_common_codon import _idx2, _idx3

def chm_diags_zero_codon(
    ncol: int,
    pver: int,
    vmr_nox_p: cobj,
    vmr_noy_p: cobj,
    vmr_clox_p: cobj,
    vmr_cloy_p: cobj,
    vmr_tcly_p: cobj,
    vmr_brox_p: cobj,
    vmr_broy_p: cobj,
    vmr_toth_p: cobj,
    vmr_tbry_p: cobj,
    vmr_foy_p: cobj,
    vmr_tfy_p: cobj,
    mmr_noy_p: cobj,
    mmr_sox_p: cobj,
    mmr_nhx_p: cobj,
    df_noy_p: cobj,
    df_sox_p: cobj,
    df_nhx_p: cobj,
):
    arrays2 = (
        Ptr[float](vmr_nox_p),
        Ptr[float](vmr_noy_p),
        Ptr[float](vmr_clox_p),
        Ptr[float](vmr_cloy_p),
        Ptr[float](vmr_tcly_p),
        Ptr[float](vmr_brox_p),
        Ptr[float](vmr_broy_p),
        Ptr[float](vmr_toth_p),
        Ptr[float](vmr_tbry_p),
        Ptr[float](vmr_foy_p),
        Ptr[float](vmr_tfy_p),
        Ptr[float](mmr_noy_p),
        Ptr[float](mmr_sox_p),
        Ptr[float](mmr_nhx_p),
    )
    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, ncol)
            for arr in arrays2:
                arr[idx] = 0.0

    arrays1 = (Ptr[float](df_noy_p), Ptr[float](df_sox_p), Ptr[float](df_nhx_p))
    for i in range(1, ncol + 1):
        idx = i - 1
        for arr in arrays1:
            arr[idx] = 0.0

def chm_diags_mass_codon(
    ncol: int,
    pver: int,
    rgrav: float,
    rearth: float,
    area_p: cobj,
    pdel_p: cobj,
    mass_p: cobj,
):
    area = Ptr[float](area_p)
    pdel = Ptr[float](pdel_p)
    mass = Ptr[float](mass_p)
    rearth2 = rearth**2

    for i in range(1, ncol + 1):
        area[i - 1] = area[i - 1] * rearth2

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, ncol)
            mass[idx] = pdel[idx] * area[i - 1] * rgrav

@inline
def _first_weight(
    m: int,
    id_cfc12: int,
    id_hcfc22: int,
    id_cf2clbr: int,
    id_h1202: int,
    id_hcfc142b: int,
    id_cof2: int,
    id_cfc113: int,
    id_cf3br: int,
    id_cfc114: int,
    id_h2402: int,
    id_cfc115: int,
) -> float:
    if (
        m == id_cfc12
        or m == id_hcfc22
        or m == id_cf2clbr
        or m == id_h1202
        or m == id_hcfc142b
        or m == id_cof2
    ):
        return 2.0
    if m == id_cfc113 or m == id_cf3br:
        return 3.0
    if m == id_cfc114 or m == id_h2402:
        return 4.0
    if m == id_cfc115:
        return 5.0
    return 1.0

@inline
def _second_weight(
    m: int,
    id_ch4: int,
    id_n2o5: int,
    id_cfc12: int,
    id_cl2: int,
    id_cl2o2: int,
    id_cfc114: int,
    id_hcfc141b: int,
    id_h1202: int,
    id_h2402: int,
    id_ch2br2: int,
    id_cfc11: int,
    id_cfc113: int,
    id_ch3ccl3: int,
    id_chbr3: int,
    id_ccl4: int,
) -> float:
    if m == id_ch4 or m == id_n2o5 or m == id_cfc12 or m == id_cl2 or m == id_cl2o2:
        return 2.0
    if m == id_cfc114 or m == id_hcfc141b or m == id_h1202 or m == id_h2402 or m == id_ch2br2:
        return 2.0
    if m == id_cfc11 or m == id_cfc113 or m == id_ch3ccl3 or m == id_chbr3:
        return 3.0
    if m == id_ccl4:
        return 4.0
    return 1.0

def chm_diags_species_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    m: int,
    do3chm_flag: int,
    fillvalue: float,
    n_molwgt: float,
    s_molwgt: float,
    id_ch4: int,
    id_n2o5: int,
    id_cfc12: int,
    id_cl2: int,
    id_cl2o2: int,
    id_cfc114: int,
    id_hcfc141b: int,
    id_h1202: int,
    id_h2402: int,
    id_ch2br2: int,
    id_cfc11: int,
    id_cfc113: int,
    id_ch3ccl3: int,
    id_chbr3: int,
    id_ccl4: int,
    id_hcfc22: int,
    id_cf2clbr: int,
    id_hcfc142b: int,
    id_cof2: int,
    id_cf3br: int,
    id_cfc115: int,
    in_foy: int,
    in_tfy: int,
    in_nox: int,
    in_noy: int,
    in_sox: int,
    in_nhx: int,
    in_clox: int,
    in_cloy: int,
    in_tcly: int,
    in_brox: int,
    in_broy: int,
    in_tbry: int,
    in_toth: int,
    vmr_p: cobj,
    mmr_p: cobj,
    depflx_p: cobj,
    mmr_tend_p: cobj,
    mass_p: cobj,
    pmid_p: cobj,
    adv_mass_p: cobj,
    ltrop_p: cobj,
    vmr_nox_p: cobj,
    vmr_noy_p: cobj,
    vmr_clox_p: cobj,
    vmr_cloy_p: cobj,
    vmr_tcly_p: cobj,
    vmr_brox_p: cobj,
    vmr_broy_p: cobj,
    vmr_toth_p: cobj,
    vmr_tbry_p: cobj,
    vmr_foy_p: cobj,
    vmr_tfy_p: cobj,
    mmr_noy_p: cobj,
    mmr_sox_p: cobj,
    mmr_nhx_p: cobj,
    df_noy_p: cobj,
    df_sox_p: cobj,
    df_nhx_p: cobj,
    net_chem_p: cobj,
    do3chm_trp_p: cobj,
    do3chm_lms_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    mmr = Ptr[float](mmr_p)
    depflx = Ptr[float](depflx_p)
    mmr_tend = Ptr[float](mmr_tend_p)
    mass = Ptr[float](mass_p)
    pmid = Ptr[float](pmid_p)
    adv_mass = Ptr[float](adv_mass_p)
    ltrop = Ptr[int](ltrop_p)

    vmr_nox = Ptr[float](vmr_nox_p)
    vmr_noy = Ptr[float](vmr_noy_p)
    vmr_clox = Ptr[float](vmr_clox_p)
    vmr_cloy = Ptr[float](vmr_cloy_p)
    vmr_tcly = Ptr[float](vmr_tcly_p)
    vmr_brox = Ptr[float](vmr_brox_p)
    vmr_broy = Ptr[float](vmr_broy_p)
    vmr_toth = Ptr[float](vmr_toth_p)
    vmr_tbry = Ptr[float](vmr_tbry_p)
    vmr_foy = Ptr[float](vmr_foy_p)
    vmr_tfy = Ptr[float](vmr_tfy_p)
    mmr_noy = Ptr[float](mmr_noy_p)
    mmr_sox = Ptr[float](mmr_sox_p)
    mmr_nhx = Ptr[float](mmr_nhx_p)
    df_noy = Ptr[float](df_noy_p)
    df_sox = Ptr[float](df_sox_p)
    df_nhx = Ptr[float](df_nhx_p)
    net_chem = Ptr[float](net_chem_p)
    do3chm_trp = Ptr[float](do3chm_trp_p)
    do3chm_lms = Ptr[float](do3chm_lms_p)

    wgt = _first_weight(
        m, id_cfc12, id_hcfc22, id_cf2clbr, id_h1202, id_hcfc142b, id_cof2, id_cfc113,
        id_cf3br, id_cfc114, id_h2402, id_cfc115
    )
    if in_foy != 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, ncol)
                vmr_foy[idx] = vmr_foy[idx] + wgt * vmr[_idx3(i, k, m, ncol, pver)]
    if in_tfy != 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, ncol)
                vmr_tfy[idx] = vmr_tfy[idx] + wgt * vmr[_idx3(i, k, m, ncol, pver)]

    wgt = _second_weight(
        m, id_ch4, id_n2o5, id_cfc12, id_cl2, id_cl2o2, id_cfc114, id_hcfc141b,
        id_h1202, id_h2402, id_ch2br2, id_cfc11, id_cfc113, id_ch3ccl3, id_chbr3,
        id_ccl4
    )
    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, ncol)
            vmr_val = vmr[_idx3(i, k, m, ncol, pver)]
            mmr_val = mmr[_idx3(i, k, m, ncol, pver)]
            if in_nox != 0:
                vmr_nox[idx] = vmr_nox[idx] + wgt * vmr_val
            if in_noy != 0:
                vmr_noy[idx] = vmr_noy[idx] + wgt * vmr_val
                mmr_noy[idx] = mmr_noy[idx] + wgt * mmr_val
            if in_sox != 0:
                mmr_sox[idx] = mmr_sox[idx] + wgt * mmr_val
            if in_nhx != 0:
                mmr_nhx[idx] = mmr_nhx[idx] + wgt * mmr_val
            if in_clox != 0:
                vmr_clox[idx] = vmr_clox[idx] + wgt * vmr_val
            if in_cloy != 0:
                vmr_cloy[idx] = vmr_cloy[idx] + wgt * vmr_val
            if in_tcly != 0:
                vmr_tcly[idx] = vmr_tcly[idx] + wgt * vmr_val
            if in_brox != 0:
                vmr_brox[idx] = vmr_brox[idx] + wgt * vmr_val
            if in_broy != 0:
                vmr_broy[idx] = vmr_broy[idx] + wgt * vmr_val
            if in_tbry != 0:
                vmr_tbry[idx] = vmr_tbry[idx] + wgt * vmr_val
            if in_toth != 0:
                vmr_toth[idx] = vmr_toth[idx] + wgt * vmr_val
            net_chem[idx] = mmr_tend[_idx3(i, k, m, ncol, pver)] * mass[idx]

    adv = adv_mass[m - 1]
    if in_noy != 0:
        for i in range(1, ncol + 1):
            idx = i - 1
            df_noy[idx] = df_noy[idx] + wgt * depflx[_idx2(i, m, ncol)] * n_molwgt / adv
    if in_sox != 0:
        for i in range(1, ncol + 1):
            idx = i - 1
            df_sox[idx] = df_sox[idx] + wgt * depflx[_idx2(i, m, ncol)] * s_molwgt / adv
    if in_nhx != 0:
        for i in range(1, ncol + 1):
            idx = i - 1
            df_nhx[idx] = df_nhx[idx] + wgt * depflx[_idx2(i, m, ncol)] * n_molwgt / adv

    if do3chm_flag != 0:
        for i in range(1, ncol + 1):
            acc = 0.0
            for k in range(ltrop[i - 1], pver + 1):
                acc = acc + net_chem[_idx2(i, k, ncol)]
            if acc == 0.0:
                acc = fillvalue
            do3chm_trp[i - 1] = acc

        for i in range(1, ncol + 1):
            acc = 0.0
            for k in range(1, pver + 1):
                if pmid[_idx2(i, k, ncol)] > 100.0e2 and k < ltrop[i - 1]:
                    acc = acc + net_chem[_idx2(i, k, ncol)]
            if acc == 0.0:
                acc = fillvalue
            do3chm_lms[i - 1] = acc

def chm_diags_species_packed_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    m: int,
    do3chm_flag: int,
    fillvalue: float,
    n_molwgt: float,
    s_molwgt: float,
    ids_p: cobj,
    flags_p: cobj,
    vmr_p: cobj,
    mmr_p: cobj,
    depflx_p: cobj,
    mmr_tend_p: cobj,
    mass_p: cobj,
    pmid_p: cobj,
    adv_mass_p: cobj,
    ltrop_p: cobj,
    vmr_nox_p: cobj,
    vmr_noy_p: cobj,
    vmr_clox_p: cobj,
    vmr_cloy_p: cobj,
    vmr_tcly_p: cobj,
    vmr_brox_p: cobj,
    vmr_broy_p: cobj,
    vmr_toth_p: cobj,
    vmr_tbry_p: cobj,
    vmr_foy_p: cobj,
    vmr_tfy_p: cobj,
    mmr_noy_p: cobj,
    mmr_sox_p: cobj,
    mmr_nhx_p: cobj,
    df_noy_p: cobj,
    df_sox_p: cobj,
    df_nhx_p: cobj,
    net_chem_p: cobj,
    do3chm_trp_p: cobj,
    do3chm_lms_p: cobj,
):
    ids = Ptr[int](ids_p)
    flags = Ptr[int](flags_p)
    return chm_diags_species_codon(
        ncol,
        pver,
        gas_pcnst,
        m,
        do3chm_flag,
        fillvalue,
        n_molwgt,
        s_molwgt,
        ids[0],
        ids[1],
        ids[2],
        ids[3],
        ids[4],
        ids[5],
        ids[6],
        ids[7],
        ids[8],
        ids[9],
        ids[10],
        ids[11],
        ids[12],
        ids[13],
        ids[14],
        ids[15],
        ids[16],
        ids[17],
        ids[18],
        ids[19],
        ids[20],
        flags[0],
        flags[1],
        flags[2],
        flags[3],
        flags[4],
        flags[5],
        flags[6],
        flags[7],
        flags[8],
        flags[9],
        flags[10],
        flags[11],
        flags[12],
        vmr_p,
        mmr_p,
        depflx_p,
        mmr_tend_p,
        mass_p,
        pmid_p,
        adv_mass_p,
        ltrop_p,
        vmr_nox_p,
        vmr_noy_p,
        vmr_clox_p,
        vmr_cloy_p,
        vmr_tcly_p,
        vmr_brox_p,
        vmr_broy_p,
        vmr_toth_p,
        vmr_tbry_p,
        vmr_foy_p,
        vmr_tfy_p,
        mmr_noy_p,
        mmr_sox_p,
        mmr_nhx_p,
        df_noy_p,
        df_sox_p,
        df_nhx_p,
        net_chem_p,
        do3chm_trp_p,
        do3chm_lms_p,
    )

@inline
def _rxt(rxt_rates: Ptr[float], i: int, k: int, rid: int, ncol: int, pver: int) -> float:
    return rxt_rates[_idx3(i, k, rid, ncol, pver)]

def chm_diags_euv_codon(
    stage: int,
    ncol: int,
    pver: int,
    indexm: int,
    id_o: int,
    id_o2: int,
    id_h: int,
    id_n: int,
    id_no: int,
    rid_scalar: int,
    rid_jeuv_p: cobj,
    vmr_p: cobj,
    rxt_rates_p: cobj,
    invariants_p: cobj,
    wrk_p: cobj,
):
    rid_jeuv = Ptr[int](rid_jeuv_p)
    vmr = Ptr[float](vmr_p)
    rxt_rates = Ptr[float](rxt_rates_p)
    invariants = Ptr[float](invariants_p)
    wrk = Ptr[float](wrk_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            inv = invariants[_idx3(i, k, indexm, ncol, pver)]
            o = vmr[_idx3(i, k, id_o, ncol, pver)]
            o2 = vmr[_idx3(i, k, id_o2, ncol, pver)]
            h = vmr[_idx3(i, k, id_h, ncol, pver)]
            n = vmr[_idx3(i, k, id_n, ncol, pver)]
            un2 = 1.0 - (o + o2 + h)
            val = 0.0
            if stage == 1:
                val = (
                    o * (
                        _rxt(rxt_rates, i, k, rid_jeuv[0], ncol, pver)
                        + _rxt(rxt_rates, i, k, rid_jeuv[1], ncol, pver)
                        + _rxt(rxt_rates, i, k, rid_jeuv[2], ncol, pver)
                        + _rxt(rxt_rates, i, k, rid_jeuv[13], ncol, pver)
                        + _rxt(rxt_rates, i, k, rid_jeuv[14], ncol, pver)
                        + _rxt(rxt_rates, i, k, rid_jeuv[15], ncol, pver)
                    )
                    + n * _rxt(rxt_rates, i, k, rid_jeuv[3], ncol, pver)
                    + o2 * (
                        _rxt(rxt_rates, i, k, rid_jeuv[4], ncol, pver)
                        + _rxt(rxt_rates, i, k, rid_jeuv[6], ncol, pver)
                        + _rxt(rxt_rates, i, k, rid_jeuv[7], ncol, pver)
                        + _rxt(rxt_rates, i, k, rid_jeuv[8], ncol, pver)
                        + _rxt(rxt_rates, i, k, rid_jeuv[16], ncol, pver)
                        + _rxt(rxt_rates, i, k, rid_jeuv[18], ncol, pver)
                        + _rxt(rxt_rates, i, k, rid_jeuv[19], ncol, pver)
                        + _rxt(rxt_rates, i, k, rid_jeuv[20], ncol, pver)
                    )
                    + un2 * (
                        _rxt(rxt_rates, i, k, rid_jeuv[5], ncol, pver)
                        + _rxt(rxt_rates, i, k, rid_jeuv[9], ncol, pver)
                        + _rxt(rxt_rates, i, k, rid_jeuv[10], ncol, pver)
                        + _rxt(rxt_rates, i, k, rid_jeuv[17], ncol, pver)
                        + _rxt(rxt_rates, i, k, rid_jeuv[21], ncol, pver)
                        + _rxt(rxt_rates, i, k, rid_jeuv[22], ncol, pver)
                    )
                )
            elif stage == 2:
                val = o * (
                    _rxt(rxt_rates, i, k, rid_jeuv[0], ncol, pver)
                    + _rxt(rxt_rates, i, k, rid_jeuv[1], ncol, pver)
                    + _rxt(rxt_rates, i, k, rid_jeuv[2], ncol, pver)
                )
            elif stage == 3:
                val = o * (
                    _rxt(rxt_rates, i, k, rid_jeuv[13], ncol, pver)
                    + _rxt(rxt_rates, i, k, rid_jeuv[14], ncol, pver)
                    + _rxt(rxt_rates, i, k, rid_jeuv[15], ncol, pver)
                )
            elif stage == 4:
                val = n * _rxt(rxt_rates, i, k, rid_jeuv[3], ncol, pver)
            elif stage == 5:
                val = o2 * (
                    _rxt(rxt_rates, i, k, rid_jeuv[4], ncol, pver)
                    + _rxt(rxt_rates, i, k, rid_jeuv[6], ncol, pver)
                    + _rxt(rxt_rates, i, k, rid_jeuv[7], ncol, pver)
                    + _rxt(rxt_rates, i, k, rid_jeuv[8], ncol, pver)
                )
            elif stage == 6:
                val = o2 * (
                    _rxt(rxt_rates, i, k, rid_jeuv[16], ncol, pver)
                    + _rxt(rxt_rates, i, k, rid_jeuv[18], ncol, pver)
                    + _rxt(rxt_rates, i, k, rid_jeuv[19], ncol, pver)
                    + _rxt(rxt_rates, i, k, rid_jeuv[20], ncol, pver)
                )
            elif stage == 7:
                val = un2 * (
                    _rxt(rxt_rates, i, k, rid_jeuv[5], ncol, pver)
                    + _rxt(rxt_rates, i, k, rid_jeuv[9], ncol, pver)
                    + _rxt(rxt_rates, i, k, rid_jeuv[10], ncol, pver)
                )
            elif stage == 8:
                val = un2 * (
                    _rxt(rxt_rates, i, k, rid_jeuv[17], ncol, pver)
                    + _rxt(rxt_rates, i, k, rid_jeuv[21], ncol, pver)
                    + _rxt(rxt_rates, i, k, rid_jeuv[22], ncol, pver)
                )
            elif stage == 9:
                val = un2 * (
                    _rxt(rxt_rates, i, k, rid_jeuv[10], ncol, pver)
                    + _rxt(rxt_rates, i, k, rid_jeuv[12], ncol, pver)
                )
            elif stage == 10:
                val = un2 * (
                    _rxt(rxt_rates, i, k, rid_jeuv[22], ncol, pver)
                    + _rxt(rxt_rates, i, k, rid_jeuv[24], ncol, pver)
                )
            elif stage == 11:
                val = vmr[_idx3(i, k, id_no, ncol, pver)] * _rxt(rxt_rates, i, k, rid_scalar, ncol, pver)
            wrk[_idx2(i, k, ncol)] = val * inv
