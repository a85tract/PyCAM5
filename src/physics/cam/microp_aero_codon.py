from math import erf, exp, log, sqrt


@export
def microp_aero_register_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def microp_aero_init_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def microp_aero_readnl_codon(value: float) -> float:
    return value


@export
def nucleate_ice_cam_readnl_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def nucleate_ice_cam_register_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def nucleate_ice_cam_init_mincld_codon(value: float) -> float:
    return value


@export
def nucleate_ice_cam_init_bulk_scale_codon(value: float) -> float:
    return value


@export
def nucleati_init_scalars_codon(
    use_preexisting_ice_in: int,
    use_hetfrz_classnuc_in: int,
    iulog_in: int,
    pi_in: float,
    mincld_in: float,
    subgrid_in: float,
    rhoice: float,
    use_preexisting_ice_p: cobj,
    use_hetfrz_classnuc_p: cobj,
    iulog_p: cobj,
    pi_p: cobj,
    mincld_p: cobj,
    subgrid_p: cobj,
    ci_p: cobj,
):
    use_preexisting_ice = Ptr[i32](use_preexisting_ice_p)
    use_hetfrz_classnuc = Ptr[i32](use_hetfrz_classnuc_p)
    iulog = Ptr[i32](iulog_p)
    pi = Ptr[float](pi_p)
    mincld = Ptr[float](mincld_p)
    subgrid = Ptr[float](subgrid_p)
    ci = Ptr[float](ci_p)

    use_preexisting_ice[0] = i32(use_preexisting_ice_in)
    use_hetfrz_classnuc[0] = i32(use_hetfrz_classnuc_in)
    iulog[0] = i32(iulog_in)
    pi[0] = pi_in
    mincld[0] = mincld_in
    subgrid[0] = subgrid_in
    ci[0] = rhoice * pi[0] / 6.0


@export
def nucleati_init_codon(
    use_preexisting_ice_in: int,
    use_hetfrz_classnuc_in: int,
    iulog_in: int,
    pi_in: float,
    mincld_in: float,
    subgrid_in: float,
    rhoice: float,
    use_preexisting_ice_p: cobj,
    use_hetfrz_classnuc_p: cobj,
    iulog_p: cobj,
    pi_p: cobj,
    mincld_p: cobj,
    subgrid_p: cobj,
    ci_p: cobj,
):
    nucleati_init_scalars_codon(
        use_preexisting_ice_in,
        use_hetfrz_classnuc_in,
        iulog_in,
        pi_in,
        mincld_in,
        subgrid_in,
        rhoice,
        use_preexisting_ice_p,
        use_hetfrz_classnuc_p,
        iulog_p,
        pi_p,
        mincld_p,
        subgrid_p,
        ci_p,
    )


@inline
def _nucleate_ice_vpreice(p_in: float, t_in: float, r_in: float, c_in: float, s_in: float) -> float:
    alphac = 0.5
    fa1c = 0.601272523
    fa2c = 0.000342181855
    fa3c = 1.49236645e-12
    wvp1c = 3.6e10
    wvp2c = 6145.0
    fvthc = 11713803.0
    thoubkc = 7.24637701e18
    svolc = 3.23e-23
    fdc = 249.239822
    fpivolc = 3.89051704e23

    t = t_in
    p = p_in * 1.0e-2
    if s_in < 1.0:
        s = 2.349 - (t / 259.0)
    else:
        s = s_in
    r = r_in * 1.0e2
    c = c_in * 1.0e-6
    t_1 = 1.0 / t
    pice = wvp1c * exp(-(wvp2c * t_1))
    alp4 = 0.25 * alphac
    flux = alp4 * sqrt(fvthc * t)
    cisat = thoubkc * pice * t_1
    a1 = (fa1c * t_1 - fa2c) * t_1
    a2 = 1.0 / cisat
    a3 = fa3c * t_1 / p
    b1 = flux * svolc * cisat * (s - 1.0)
    b2 = flux * fdc * p * (t_1 ** 1.94)
    dloss = fpivolc * c * b1 * (r ** 2.0) / (1.0 + b2 * r)
    vice = (a2 + a3 * s) * dloss / (a1 * s)
    return vice * 1.0e-2


@inline
def _nucleate_ice_frachom(tmean: float, rhimean: float, detat: float, pi_value: float) -> float:
    seta = 6132.9
    nbin = 200
    sihom = 2.349 - tmean / 259.0
    fhom = 0.0

    for i in range(nbin, 0, -1):
        deta = (float(i) - 0.5 - float(nbin) / 2.0) * 6.0 / float(nbin)
        sbin = rhimean * exp(deta * detat * seta / (tmean ** 2.0))
        pdf_t = exp(-(deta ** 2.0) / 2.0) * 6.0 / (sqrt(2.0 * pi_value) * float(nbin))
        if sbin >= sihom:
            fhom = fhom + pdf_t
        else:
            break

    return fhom / 0.997


@inline
def _nucleate_ice_hetero_nis(t: float, ww: float, ns: float) -> float:
    a11 = 0.0263
    a12 = -0.0185
    a21 = 2.758
    a22 = 1.3221
    b11 = -0.008
    b12 = -0.0468
    b21 = -0.2667
    b22 = -1.4588

    b = (a11 + b11 * log(ns)) * log(ww) + (a12 + b12 * log(ns))
    c = a21 + b21 * log(ns)

    nis = exp(a22) * (ns ** b22) * exp(b * t) * (ww ** c)
    return min(nis, ns)


@inline
def _nucleate_ice_hf_value(t: float, ww: float, rh: float, na: float, subgrid: float) -> float:
    a1_fast = 0.0231
    a21_fast = -1.6387
    a22_fast = -6.045
    b1_fast = -0.008
    b21_fast = -0.042
    b22_fast = -0.112
    c1_fast = 0.0739
    c2_fast = 1.2372

    a1_slow = -0.3949
    a2_slow = 1.282
    b1_slow = -0.0156
    b2_slow = 0.0111
    b3_slow = 0.0217
    c1_slow = 0.120
    c2_slow = 2.312

    ni = 0.0
    a = 6.0e-4 * log(ww) + 6.6e-3
    b = 6.0e-2 * log(ww) + 1.052
    c = 1.68 * log(ww) + 129.35
    rhw = (a * t * t + b * t + c) * 0.01

    if (t <= -37.0) and ((rh * subgrid) >= rhw):
        regm = 6.07 * log(ww) - 55.0
        if t >= regm:
            a2_fast = a21_fast
            b2_fast = b21_fast
            if t <= -64.0:
                a2_fast = a22_fast
                b2_fast = b22_fast

            k1_fast = exp(a2_fast + b2_fast * t + c2_fast * log(ww))
            k2_fast = a1_fast + b1_fast * t + c1_fast * log(ww)

            ni = k1_fast * (na ** k2_fast)
            ni = min(ni, na)
        else:
            k1_slow = exp(a2_slow + (b2_slow + b3_slow * log(ww)) * t + c2_slow * log(ww))
            k2_slow = a1_slow + b1_slow * t + c1_slow * log(ww)

            ni = k1_slow * (na ** k2_slow)
            ni = min(ni, na)

    return ni


@export
def nucleati_codon(
    wbar: float,
    tair: float,
    pmid: float,
    relhum: float,
    cldn: float,
    qc: float,
    qi: float,
    ni_in: float,
    rhoair: float,
    so4_num: float,
    dst_num: float,
    soot_num: float,
    svp_water_tair: float,
    svp_ice_tair: float,
    use_preexisting_ice: int,
    use_hetfrz_classnuc: int,
    mincld: float,
    subgrid: float,
    ci: float,
    shet: float,
    minweff: float,
    gamma4: float,
    pi_value: float,
    nuci_p: cobj,
    onihf_p: cobj,
    oniimm_p: cobj,
    onidep_p: cobj,
    onimey_p: cobj,
    wpice_p: cobj,
    weff_p: cobj,
    fhom_p: cobj,
    warn_p: cobj,
    warn_ni_p: cobj,
    warn_nihf_p: cobj,
    warn_niimm_p: cobj,
    warn_nidep_p: cobj,
    warn_deles_p: cobj,
    warn_esi_p: cobj,
):
    nuci_out = Ptr[float](nuci_p)
    onihf_out = Ptr[float](onihf_p)
    oniimm_out = Ptr[float](oniimm_p)
    onidep_out = Ptr[float](onidep_p)
    onimey_out = Ptr[float](onimey_p)
    wpice_out = Ptr[float](wpice_p)
    weff_out = Ptr[float](weff_p)
    fhom_out = Ptr[float](fhom_p)
    warn = Ptr[int](warn_p)
    warn_ni = Ptr[float](warn_ni_p)
    warn_nihf = Ptr[float](warn_nihf_p)
    warn_niimm = Ptr[float](warn_niimm_p)
    warn_nidep = Ptr[float](warn_nidep_p)
    warn_deles = Ptr[float](warn_deles_p)
    warn_esi = Ptr[float](warn_esi_p)

    wbar1 = wbar
    wbar2 = wbar
    wpice = 0.0
    weff = 0.0
    fhom = 0.0
    wpicehet = 0.0
    ni_preice = 0.0

    if use_preexisting_ice != 0:
        ni_preice = ni_in * rhoair
        ni_preice = ni_preice / max(mincld, cldn)

        if ni_preice > 10.0:
            shom = -1.5
            lami = (gamma4 * ci * ni_in / qi) ** (1.0 / 3.0)
            ri_preice = 0.5 / lami
            ri_preice = max(ri_preice, 1.0e-8)
            wpice = _nucleate_ice_vpreice(pmid, tair, ri_preice, ni_preice, shom)
            wpicehet = _nucleate_ice_vpreice(pmid, tair, ri_preice, ni_preice, shet)
        else:
            wpice = 0.0
            wpicehet = 0.0

        weff = max(wbar - wpice, minweff)
        wpice = min(wpice, wbar)
        weffhet = max(wbar - wpicehet, minweff)
        wpicehet = min(wpicehet, wbar)

        wbar1 = weff
        wbar2 = weffhet

        detat = wbar / 0.23
        rhimean = 1.0
        fhom = _nucleate_ice_frachom(tair, rhimean, detat, pi_value)

    ni = 0.0
    tc = tair - 273.15
    niimm = 0.0
    nidep = 0.0
    nihf = 0.0
    deles = 0.0
    esi = 0.0

    if so4_num >= 1.0e-10 and (soot_num + dst_num) >= 1.0e-10 and cldn > 0.0:
        rhi_gate = relhum * svp_water_tair / svp_ice_tair * subgrid
        if (tc <= -35.0) and (rhi_gate >= 1.2):
            a = -1.4938 * log(soot_num + dst_num) + 12.884
            b = -10.41 * log(soot_num + dst_num) - 67.69
            regm = a * log(wbar1) + b

            if tc > regm:
                if tc < -40.0 and wbar1 > 1.0:
                    nihf = _nucleate_ice_hf_value(tc, wbar1, relhum, so4_num, subgrid)
                    niimm = 0.0
                    nidep = 0.0
                    if use_preexisting_ice != 0:
                        if nihf > 1.0e-3:
                            niimm = min(dst_num, ni_preice * 1.0e-6)
                            nihf = nihf + ni_preice * 1.0e-6 - niimm
                        nihf = nihf * fhom
                        n1 = nihf + niimm
                    else:
                        n1 = nihf
                else:
                    niimm = _nucleate_ice_hetero_nis(tc, wbar2, soot_num + dst_num)
                    nidep = 0.0
                    if use_preexisting_ice != 0:
                        if niimm > 1.0e-6:
                            niimm = niimm + ni_preice * 1.0e-6
                            niimm = min(dst_num, niimm)
                    nihf = 0.0
                    n1 = niimm + nidep
            elif tc < regm - 5.0:
                nihf = _nucleate_ice_hf_value(tc, wbar1, relhum, so4_num, subgrid)
                niimm = 0.0
                nidep = 0.0
                if use_preexisting_ice != 0:
                    if nihf > 1.0e-3:
                        niimm = min(dst_num, ni_preice * 1.0e-6)
                        nihf = nihf + ni_preice * 1.0e-6 - niimm
                    nihf = nihf * fhom
                    n1 = nihf + niimm
                else:
                    n1 = nihf
            else:
                if tc < -40.0 and wbar1 > 1.0:
                    nihf = _nucleate_ice_hf_value(tc, wbar1, relhum, so4_num, subgrid)
                    niimm = 0.0
                    nidep = 0.0
                    if use_preexisting_ice != 0:
                        if nihf > 1.0e-3:
                            niimm = min(dst_num, ni_preice * 1.0e-6)
                            nihf = nihf + ni_preice * 1.0e-6 - niimm
                        nihf = nihf * fhom
                        n1 = nihf + niimm
                    else:
                        n1 = nihf
                else:
                    nihf = _nucleate_ice_hf_value(regm - 5.0, wbar1, relhum, so4_num, subgrid)
                    niimm = _nucleate_ice_hetero_nis(regm, wbar2, soot_num + dst_num)
                    nidep = 0.0

                    if use_preexisting_ice != 0:
                        nihf = nihf * fhom

                    if nihf <= (niimm + nidep):
                        n1 = nihf
                    else:
                        n1 = (niimm + nidep) * (((niimm + nidep) / nihf) ** ((tc - regm) / 5.0))

                    if use_preexisting_ice != 0:
                        if n1 > 1.0e-3:
                            n1 = n1 + ni_preice * 1.0e-6
                            niimm = min(dst_num, n1)
                            nihf = n1 - niimm
                        else:
                            n1 = 0.0
                            niimm = 0.0
                            nihf = 0.0
            ni = n1

    if tc < 0.0 and tc > -37.0 and qc > 1.0e-12:
        esl = svp_water_tair
        esi = svp_ice_tair
        deles = esl - esi
        nimey = 1.0e-3 * exp(12.96 * deles / esi - 0.639)
    else:
        nimey = 0.0

    if use_hetfrz_classnuc != 0:
        nimey = 0.0

    nuci = ni + nimey
    warn[0] = 0
    warn_ni[0] = ni
    warn_nihf[0] = nihf
    warn_niimm[0] = niimm
    warn_nidep[0] = nidep
    warn_deles[0] = deles
    warn_esi[0] = esi
    if nuci > 9999.0 or nuci < 0.0:
        warn[0] = 1
        nuci = 0.0

    nuci_out[0] = nuci * 1.0e6 / rhoair
    onimey_out[0] = nimey * 1.0e6 / rhoair
    onidep_out[0] = nidep * 1.0e6 / rhoair
    oniimm_out[0] = niimm * 1.0e6 / rhoair
    onihf_out[0] = nihf * 1.0e6 / rhoair
    wpice_out[0] = wpice
    weff_out[0] = weff
    fhom_out[0] = fhom


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


@export
def microp_aero_rho_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    rair: float,
    pmid_p: cobj,
    t_p: cobj,
    rho_p: cobj,
):
    pmid = Ptr[float](pmid_p)
    t = Ptr[float](t_p)
    rho = Ptr[float](rho_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            rho[_idx2(i, k, pcols)] = pmid[_idx2(i, k, pcols)] / (
                rair * t[_idx2(i, k, pcols)]
            )


@export
def microp_aero_diag_tke_wsub_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    use_preexisting_ice_flag: int,
    tke_p: cobj,
    wsub_p: cobj,
    wsubi_p: cobj,
):
    tke = Ptr[float](tke_p)
    wsub = Ptr[float](wsub_p)
    wsubi = Ptr[float](wsubi_p)

    for k in range(1, top_lev):
        for i in range(1, ncol + 1):
            wsub[_idx2(i, k, pcols)] = 0.20
            wsubi[_idx2(i, k, pcols)] = 0.001

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            val = sqrt(
                0.5
                * (tke[_idx2(i, k, pcols)] + tke[_idx2(i, k + 1, pcols)])
                * (2.0 / 3.0)
            )
            if val > 10.0:
                val = 10.0

            ice_val = val
            if ice_val < 0.001:
                ice_val = 0.001
            if use_preexisting_ice_flag == 0:
                if ice_val > 0.2:
                    ice_val = 0.2

            if val < 0.20:
                val = 0.20

            wsub[_idx2(i, k, pcols)] = val
            wsubi[_idx2(i, k, pcols)] = ice_val


@export
def microp_aero_lcldm_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    mincld: float,
    ast_p: cobj,
    lcldm_p: cobj,
):
    ast = Ptr[float](ast_p)
    lcldm = Ptr[float](lcldm_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            val = ast[_idx2(i, k, pcols)]
            if val < mincld:
                val = mincld
            lcldm[_idx2(i, k, pcols)] = val


@export
def microp_aero_modal_lcloud_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    qsmall: float,
    qc_p: cobj,
    qi_p: cobj,
    cldn_p: cobj,
    cldo_p: cobj,
    lcldn_p: cobj,
    lcldo_p: cobj,
):
    qc = Ptr[float](qc_p)
    qi = Ptr[float](qi_p)
    cldn = Ptr[float](cldn_p)
    cldo = Ptr[float](cldo_p)
    lcldn = Ptr[float](lcldn_p)
    lcldo = Ptr[float](lcldo_p)

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            lcldn[_idx2(i, k, pcols)] = 0.0
            lcldo[_idx2(i, k, pcols)] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            qcld = qc[_idx2(i, k, pcols)] + qi[_idx2(i, k, pcols)]
            if qcld > qsmall:
                lcldn[_idx2(i, k, pcols)] = (
                    cldn[_idx2(i, k, pcols)] * qc[_idx2(i, k, pcols)] / qcld
                )
                lcldo[_idx2(i, k, pcols)] = (
                    cldo[_idx2(i, k, pcols)] * qc[_idx2(i, k, pcols)] / qcld
                )


@export
def microp_aero_modal_contact_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    nmodes: int,
    mode_coarse_dst_idx: int,
    separate_dust_flag: int,
    rn_dst3: float,
    t_p: cobj,
    rho_p: cobj,
    coarse_dust_p: cobj,
    coarse_nacl_p: cobj,
    num_coarse_p: cobj,
    dgnumwet_p: cobj,
    nacon_p: cobj,
    rndst_p: cobj,
):
    t = Ptr[float](t_p)
    rho = Ptr[float](rho_p)
    coarse_dust = Ptr[float](coarse_dust_p)
    coarse_nacl = Ptr[float](coarse_nacl_p)
    num_coarse = Ptr[float](num_coarse_p)
    dgnumwet = Ptr[float](dgnumwet_p)
    nacon = Ptr[float](nacon_p)
    rndst = Ptr[float](rndst_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            if t[_idx2(i, k, pcols)] < 269.15:
                dmc = coarse_dust[_idx2(i, k, pcols)]
                ssmc = coarse_nacl[_idx2(i, k, pcols)]

                if separate_dust_flag != 0:
                    wght = 1.0
                else:
                    wght = dmc / (ssmc + dmc)

                if dmc > 0.0:
                    nacon[_idx3(i, k, 3, pcols, pver)] = (
                        wght * num_coarse[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)]
                    )
                else:
                    nacon[_idx3(i, k, 3, pcols, pver)] = 0.0

                radius = 0.5 * dgnumwet[_idx3(i, k, mode_coarse_dst_idx, pcols, pver)]
                if radius <= 0.0:
                    radius = rn_dst3
                rndst[_idx3(i, k, 3, pcols, pver)] = radius


@export
def microp_aero_npccn_copy_codon(
    ncol: int,
    pcols: int,
    pver: int,
    nctend_mixnuc_p: cobj,
    npccn_p: cobj,
):
    nctend_mixnuc = Ptr[float](nctend_mixnuc_p)
    npccn = Ptr[float](npccn_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            npccn[_idx2(i, k, pcols)] = nctend_mixnuc[_idx2(i, k, pcols)]


@export
def nucleate_ice_cam_rho_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    rair: float,
    pmid_p: cobj,
    t_p: cobj,
    rho_p: cobj,
):
    pmid = Ptr[float](pmid_p)
    t = Ptr[float](t_p)
    rho = Ptr[float](rho_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            rho[_idx2(i, k, pcols)] = pmid[_idx2(i, k, pcols)] / (
                rair * t[_idx2(i, k, pcols)]
            )


@export
def nucleate_ice_cam_icecldf_codon(
    ncol: int,
    pcols: int,
    pver: int,
    ast_p: cobj,
    icecldf_p: cobj,
):
    ast = Ptr[float](ast_p)
    icecldf = Ptr[float](icecldf_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            icecldf[_idx2(i, k, pcols)] = ast[_idx2(i, k, pcols)]


@export
def nucleate_ice_cam_zero_outputs_codon(
    ncol: int,
    pcols: int,
    pver: int,
    use_preexisting_ice_flag: int,
    naai_p: cobj,
    naai_hom_p: cobj,
    nihf_p: cobj,
    niimm_p: cobj,
    nidep_p: cobj,
    nimey_p: cobj,
    fhom_p: cobj,
    wice_p: cobj,
    weff_p: cobj,
    innso4_p: cobj,
    innbc_p: cobj,
    inndust_p: cobj,
    inhet_p: cobj,
    inhom_p: cobj,
    infrehom_p: cobj,
    infrein_p: cobj,
):
    naai = Ptr[float](naai_p)
    naai_hom = Ptr[float](naai_hom_p)
    nihf = Ptr[float](nihf_p)
    niimm = Ptr[float](niimm_p)
    nidep = Ptr[float](nidep_p)
    nimey = Ptr[float](nimey_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            naai[_idx2(i, k, pcols)] = 0.0
            naai_hom[_idx2(i, k, pcols)] = 0.0
            nihf[_idx2(i, k, pcols)] = 0.0
            niimm[_idx2(i, k, pcols)] = 0.0
            nidep[_idx2(i, k, pcols)] = 0.0
            nimey[_idx2(i, k, pcols)] = 0.0

    if use_preexisting_ice_flag != 0:
        fhom = Ptr[float](fhom_p)
        wice = Ptr[float](wice_p)
        weff = Ptr[float](weff_p)
        innso4 = Ptr[float](innso4_p)
        innbc = Ptr[float](innbc_p)
        inndust = Ptr[float](inndust_p)
        inhet = Ptr[float](inhet_p)
        inhom = Ptr[float](inhom_p)
        infrehom = Ptr[float](infrehom_p)
        infrein = Ptr[float](infrein_p)

        for k in range(1, pver + 1):
            for i in range(1, pcols + 1):
                fhom[_idx2(i, k, pcols)] = 0.0
                wice[_idx2(i, k, pcols)] = 0.0
                weff[_idx2(i, k, pcols)] = 0.0
                innso4[_idx2(i, k, pcols)] = 0.0
                innbc[_idx2(i, k, pcols)] = 0.0
                inndust[_idx2(i, k, pcols)] = 0.0
                inhet[_idx2(i, k, pcols)] = 0.0
                inhom[_idx2(i, k, pcols)] = 0.0
                infrehom[_idx2(i, k, pcols)] = 0.0
                infrein[_idx2(i, k, pcols)] = 0.0


@export
def nucleate_ice_cam_relhum_codon(
    ncol: int,
    pcols: int,
    k: int,
    mincld: float,
    qn_p: cobj,
    qs_p: cobj,
    icecldf_p: cobj,
    relhum_p: cobj,
    icldm_p: cobj,
):
    qn = Ptr[float](qn_p)
    qs = Ptr[float](qs_p)
    icecldf = Ptr[float](icecldf_p)
    relhum = Ptr[float](relhum_p)
    icldm = Ptr[float](icldm_p)

    for i in range(1, ncol + 1):
        relhum[_idx2(i, k, pcols)] = qn[_idx2(i, k, pcols)] / qs[i - 1]
        val = icecldf[_idx2(i, k, pcols)]
        if val < mincld:
            val = mincld
        icldm[_idx2(i, k, pcols)] = val


@export
def nucleate_ice_cam_post_nucleati_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    tmelt: float,
    t_p: cobj,
    rho_p: cobj,
    naai_hom_p: cobj,
    nihf_p: cobj,
    niimm_p: cobj,
    nidep_p: cobj,
    nimey_p: cobj,
):
    t = Ptr[float](t_p)
    rho = Ptr[float](rho_p)
    naai_hom = Ptr[float](naai_hom_p)
    nihf = Ptr[float](nihf_p)
    niimm = Ptr[float](niimm_p)
    nidep = Ptr[float](nidep_p)
    nimey = Ptr[float](nimey_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            if t[idx] < tmelt - 5.0:
                rho_val = rho[idx]
                nihf_val = nihf[idx]
                naai_hom[idx] = nihf_val
                nihf[idx] = nihf_val * rho_val
                niimm[idx] = niimm[idx] * rho_val
                nidep[idx] = nidep[idx] * rho_val
                nimey[idx] = nimey[idx] * rho_val


@export
def nucleate_ice_cam_modal_dst_num_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    separate_dust_flag: int,
    rho_p: cobj,
    coarse_dust_p: cobj,
    coarse_nacl_p: cobj,
    num_coarse_p: cobj,
    dst_num_p: cobj,
):
    rho = Ptr[float](rho_p)
    coarse_dust = Ptr[float](coarse_dust_p)
    coarse_nacl = Ptr[float](coarse_nacl_p)
    num_coarse = Ptr[float](num_coarse_p)
    dst_num = Ptr[float](dst_num_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            dmc = coarse_dust[idx] * rho[idx]
            ssmc = coarse_nacl[idx] * rho[idx]

            if dmc > 0.0:
                if separate_dust_flag != 0:
                    wght = 1.0
                else:
                    wght = dmc / (ssmc + dmc)
                dst_num[idx] = wght * num_coarse[idx] * rho[idx] * 1.0e-6
            else:
                dst_num[idx] = 0.0


@export
def nucleate_ice_cam_modal_so4_num_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    mode_aitken_idx: int,
    tmelt: float,
    sigmag_aitken: float,
    t_p: cobj,
    rho_p: cobj,
    num_aitken_p: cobj,
    dgnum_p: cobj,
    so4_num_p: cobj,
):
    t = Ptr[float](t_p)
    rho = Ptr[float](rho_p)
    num_aitken = Ptr[float](num_aitken_p)
    dgnum = Ptr[float](dgnum_p)
    so4_num = Ptr[float](so4_num_p)

    log_sigmag = log(sigmag_aitken)
    sqrt_two = 2.0 ** 0.5

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            val = 0.0
            if t[idx] < tmelt - 5.0:
                dg = dgnum[_idx3(i, k, mode_aitken_idx, pcols, pver)]
                if dg > 0.0:
                    val = (
                        num_aitken[idx]
                        * rho[idx]
                        * 1.0e-6
                        * (
                            0.5
                            - 0.5
                            * erf(log(0.1e-6 / dg) / (sqrt_two * log_sigmag))
                        )
                    )
                    if val < 0.0:
                        val = 0.0
            so4_num[idx] = val


@export
def nucleate_ice_hetero_codon(T: float, ww: float, Ns: float, Nis_p: cobj, Nid_p: cobj):
    Nis = Ptr[float](Nis_p)
    Nid = Ptr[float](Nid_p)

    A11 = 0.0263
    A12 = -0.0185
    A21 = 2.758
    A22 = 1.3221
    B11 = -0.008
    B12 = -0.0468
    B21 = -0.2667
    B22 = -1.4588

    B = (A11 + B11 * log(Ns)) * log(ww) + (A12 + B12 * log(Ns))
    C = A21 + B21 * log(Ns)

    nis_val = exp(A22) * (Ns ** B22) * exp(B * T) * (ww ** C)
    if nis_val > Ns:
        nis_val = Ns
    Nis[0] = nis_val
    Nid[0] = 0.0


@export
def nucleate_ice_hf_codon(T: float, ww: float, RH: float, Na: float, subgrid: float, Ni_p: cobj):
    Ni = Ptr[float](Ni_p)

    A1_fast = 0.0231
    A21_fast = -1.6387
    A22_fast = -6.045
    B1_fast = -0.008
    B21_fast = -0.042
    B22_fast = -0.112
    C1_fast = 0.0739
    C2_fast = 1.2372

    A1_slow = -0.3949
    A2_slow = 1.282
    B1_slow = -0.0156
    B2_slow = 0.0111
    B3_slow = 0.0217
    C1_slow = 0.120
    C2_slow = 2.312

    ni_val = 0.0

    A = 6.0e-4 * log(ww) + 6.6e-3
    B = 6.0e-2 * log(ww) + 1.052
    C = 1.68 * log(ww) + 129.35
    RHw = (A * T * T + B * T + C) * 0.01

    if (T <= -37.0) and ((RH * subgrid) >= RHw):
        regm = 6.07 * log(ww) - 55.0

        if T >= regm:
            if T > -64.0:
                A2_fast = A21_fast
                B2_fast = B21_fast
            else:
                A2_fast = A22_fast
                B2_fast = B22_fast

            k1_fast = exp(A2_fast + B2_fast * T + C2_fast * log(ww))
            k2_fast = A1_fast + B1_fast * T + C1_fast * log(ww)

            ni_val = k1_fast * (Na ** k2_fast)
            if ni_val > Na:
                ni_val = Na
        else:
            k1_slow = exp(A2_slow + (B2_slow + B3_slow * log(ww)) * T + C2_slow * log(ww))
            k2_slow = A1_slow + B1_slow * T + C1_slow * log(ww)

            ni_val = k1_slow * (Na ** k2_slow)
            if ni_val > Na:
                ni_val = Na

    Ni[0] = ni_val
