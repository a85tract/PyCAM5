from math import sqrt


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
