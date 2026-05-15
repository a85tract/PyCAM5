@inline
def _idx2(i: int, k: int, ld1: int) -> int:
    return (i - 1) + (k - 1) * ld1


@inline
def _idx3(i: int, k: int, m: int, ld1: int, ld2: int) -> int:
    return (i - 1) + (k - 1) * ld1 + (m - 1) * ld1 * ld2


@export
def convect_shallow_select_scheme_codon(
    scheme_len: int,
    scheme_ascii_p: cobj,
    scheme_code_p: cobj,
    status_p: cobj,
):
    scheme_ascii = Ptr[int](scheme_ascii_p)
    scheme_code = Ptr[int](scheme_code_p)
    status = Ptr[int](status_p)

    status[0] = 0
    scheme_code[0] = 0

    n = scheme_len
    while n > 0 and scheme_ascii[n - 1] == 32:
        n -= 1

    if n == 2:
        c1 = scheme_ascii[0]
        c2 = scheme_ascii[1]

        if c1 >= 65 and c1 <= 90:
            c1 += 32
        if c2 >= 65 and c2 <= 90:
            c2 += 32

        if c1 == 117 and c2 == 119:
            scheme_code[0] = 3
            return

    if n == 3:
        c1 = scheme_ascii[0]
        c2 = scheme_ascii[1]
        c3 = scheme_ascii[2]

        if c1 >= 65 and c1 <= 90:
            c1 += 32
        if c2 >= 65 and c2 <= 90:
            c2 += 32
        if c3 >= 65 and c3 <= 90:
            c3 += 32

        if c1 == 111 and c2 == 102 and c3 == 102:
            scheme_code[0] = 1
            return

    if n == 4:
        c1 = scheme_ascii[0]
        c2 = scheme_ascii[1]
        c3 = scheme_ascii[2]
        c4 = scheme_ascii[3]

        if c1 >= 65 and c1 <= 90:
            c1 += 32
        if c2 >= 65 and c2 <= 90:
            c2 += 32
        if c3 >= 65 and c3 <= 90:
            c3 += 32
        if c4 >= 65 and c4 <= 90:
            c4 += 32

        if c1 == 104 and c2 == 97 and c3 == 99 and c4 == 107:
            scheme_code[0] = 2
            return

    if n == 6:
        c1 = scheme_ascii[0]
        c2 = scheme_ascii[1]
        c3 = scheme_ascii[2]
        c4 = scheme_ascii[3]
        c5 = scheme_ascii[4]
        c6 = scheme_ascii[5]

        if c1 >= 65 and c1 <= 90:
            c1 += 32
        if c2 >= 65 and c2 <= 90:
            c2 += 32
        if c3 >= 65 and c3 <= 90:
            c3 += 32
        if c4 >= 65 and c4 <= 90:
            c4 += 32
        if c5 >= 65 and c5 <= 90:
            c5 += 32
        if c6 >= 65 and c6 <= 90:
            c6 += 32

        if c1 == 117 and c2 == 110 and c3 == 105 and c4 == 99 and c5 == 111 and c6 == 110:
            scheme_code[0] = 4
            return

    if n == 9:
        c1 = scheme_ascii[0]
        c2 = scheme_ascii[1]
        c3 = scheme_ascii[2]
        c4 = scheme_ascii[3]
        c5 = scheme_ascii[4]
        c6 = scheme_ascii[5]
        c7 = scheme_ascii[6]
        c8 = scheme_ascii[7]
        c9 = scheme_ascii[8]

        if c1 >= 65 and c1 <= 90:
            c1 += 32
        if c2 >= 65 and c2 <= 90:
            c2 += 32
        if c3 >= 65 and c3 <= 90:
            c3 += 32
        if c4 >= 65 and c4 <= 90:
            c4 += 32
        if c5 >= 65 and c5 <= 90:
            c5 += 32
        if c6 >= 65 and c6 <= 90:
            c6 += 32
        if c7 >= 65 and c7 <= 90:
            c7 += 32
        if c8 >= 65 and c8 <= 90:
            c8 += 32
        if c9 >= 65 and c9 <= 90:
            c9 += 32

        if (
            c1 == 99 and c2 == 108 and c3 == 117 and c4 == 98 and c5 == 98
            and c6 == 95 and c7 == 115 and c8 == 103 and c9 == 115
        ):
            scheme_code[0] = 1
            return

    status[0] = 1


@export
def convect_shallow_init_shell_codon(
    ncol: int,
    pcols: int,
    tpert_p: cobj,
    landfracdum_p: cobj,
):
    tpert = Ptr[float](tpert_p)
    landfracdum = Ptr[float](landfracdum_p)

    for i in range(1, ncol + 1):
        idx1 = i - 1
        tpert[idx1] = 0.0
        landfracdum[idx1] = 0.0


@export
def convect_shallow_ptend_lq_mask_shell_codon(pcnst: int, lq_mask_p: cobj):
    lq_mask = Ptr[int](lq_mask_p)

    for m in range(0, pcnst):
        lq_mask[m] = 1


@export
def uwshcu_output_init_shell_codon(
    mix: int,
    mkx: int,
    iend: int,
    ncnst: int,
    umf_p: cobj,
    slflx_p: cobj,
    qtflx_p: cobj,
    flxprc1_p: cobj,
    flxsnow1_p: cobj,
    qvten_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    sten_p: cobj,
    uten_p: cobj,
    vten_p: cobj,
    qrten_p: cobj,
    qsten_p: cobj,
    evapc_p: cobj,
    cufrc_p: cobj,
    qcu_p: cobj,
    qlu_p: cobj,
    qiu_p: cobj,
    fer_p: cobj,
    fdr_p: cobj,
    qc_p: cobj,
    qtten_p: cobj,
    slten_p: cobj,
    ufrc_p: cobj,
    uflx_p: cobj,
    vflx_p: cobj,
    trten_p: cobj,
    trflx_p: cobj,
    wtqc_p: cobj,
    wtprec_p: cobj,
    wtsnow_p: cobj,
    precip_p: cobj,
    snow_p: cobj,
    cinh_p: cobj,
    cinlclh_p: cobj,
    cbmf_p: cobj,
    rliq_p: cobj,
    cnt_p: cobj,
    cnb_p: cobj,
):
    umf = Ptr[float](umf_p)
    slflx = Ptr[float](slflx_p)
    qtflx = Ptr[float](qtflx_p)
    flxprc1 = Ptr[float](flxprc1_p)
    flxsnow1 = Ptr[float](flxsnow1_p)
    qvten = Ptr[float](qvten_p)
    qlten = Ptr[float](qlten_p)
    qiten = Ptr[float](qiten_p)
    sten = Ptr[float](sten_p)
    uten = Ptr[float](uten_p)
    vten = Ptr[float](vten_p)
    qrten = Ptr[float](qrten_p)
    qsten = Ptr[float](qsten_p)
    evapc = Ptr[float](evapc_p)
    cufrc = Ptr[float](cufrc_p)
    qcu = Ptr[float](qcu_p)
    qlu = Ptr[float](qlu_p)
    qiu = Ptr[float](qiu_p)
    fer = Ptr[float](fer_p)
    fdr = Ptr[float](fdr_p)
    qc = Ptr[float](qc_p)
    qtten = Ptr[float](qtten_p)
    slten = Ptr[float](slten_p)
    ufrc = Ptr[float](ufrc_p)
    uflx = Ptr[float](uflx_p)
    vflx = Ptr[float](vflx_p)
    trten = Ptr[float](trten_p)
    trflx = Ptr[float](trflx_p)
    wtqc = Ptr[float](wtqc_p)
    wtprec = Ptr[float](wtprec_p)
    wtsnow = Ptr[float](wtsnow_p)
    precip = Ptr[float](precip_p)
    snow = Ptr[float](snow_p)
    cinh = Ptr[float](cinh_p)
    cinlclh = Ptr[float](cinlclh_p)
    cbmf = Ptr[float](cbmf_p)
    rliq = Ptr[float](rliq_p)
    cnt = Ptr[float](cnt_p)
    cnb = Ptr[float](cnb_p)

    k = 0
    while k <= mkx:
        i = 0
        while i < iend:
            idx = i + k * mix
            umf[idx] = 0.0
            slflx[idx] = 0.0
            qtflx[idx] = 0.0
            flxprc1[idx] = 0.0
            flxsnow1[idx] = 0.0
            ufrc[idx] = 0.0
            uflx[idx] = 0.0
            vflx[idx] = 0.0
            i += 1
        k += 1

    k = 0
    while k < mkx:
        i = 0
        while i < iend:
            idx = i + k * mix
            qvten[idx] = 0.0
            qlten[idx] = 0.0
            qiten[idx] = 0.0
            sten[idx] = 0.0
            uten[idx] = 0.0
            vten[idx] = 0.0
            qrten[idx] = 0.0
            qsten[idx] = 0.0
            evapc[idx] = 0.0
            cufrc[idx] = 0.0
            qcu[idx] = 0.0
            qlu[idx] = 0.0
            qiu[idx] = 0.0
            fer[idx] = 0.0
            fdr[idx] = 0.0
            qc[idx] = 0.0
            qtten[idx] = 0.0
            slten[idx] = 0.0
            i += 1
        k += 1

    i = 0
    while i < iend:
        precip[i] = 0.0
        snow[i] = 0.0
        cinh[i] = -1.0
        cinlclh[i] = -1.0
        cbmf[i] = 0.0
        rliq[i] = 0.0
        cnt[i] = float(mkx)
        cnb[i] = 0.0
        i += 1

    m = 0
    while m < ncnst:
        offset_2d = m * mix
        offset_3d = m * mix * mkx
        offset_3dz = m * mix * (mkx + 1)
        i = 0
        while i < iend:
            wtprec[i + offset_2d] = 0.0
            wtsnow[i + offset_2d] = 0.0
            i += 1
        k = 0
        while k < mkx:
            i = 0
            while i < iend:
                idx = i + k * mix
                trten[idx + offset_3d] = 0.0
                wtqc[idx + offset_3d] = 0.0
                i += 1
            k += 1
        k = 0
        while k <= mkx:
            i = 0
            while i < iend:
                trflx[i + k * mix + offset_3dz] = 0.0
                i += 1
            k += 1
        m += 1


@export
def uwshcu_diag_init_shell_codon(
    mix: int,
    mkx: int,
    iend: int,
    ncnst: int,
    ufrcinvbase_p: cobj,
    ufrclcl_p: cobj,
    winvbase_p: cobj,
    wlcl_p: cobj,
    plcl_p: cobj,
    pinv_p: cobj,
    plfc_p: cobj,
    pbup_p: cobj,
    ppen_p: cobj,
    qtsrc_p: cobj,
    thlsrc_p: cobj,
    thvlsrc_p: cobj,
    emfkbup_p: cobj,
    cbmflimit_p: cobj,
    tkeavg_p: cobj,
    zinv_p: cobj,
    rcwp_p: cobj,
    rlwp_p: cobj,
    riwp_p: cobj,
    wu_p: cobj,
    qtu_p: cobj,
    thlu_p: cobj,
    thvu_p: cobj,
    uu_p: cobj,
    vu_p: cobj,
    qtu_emf_p: cobj,
    thlu_emf_p: cobj,
    uu_emf_p: cobj,
    vu_emf_p: cobj,
    uemf_p: cobj,
    tru_p: cobj,
    tru_emf_p: cobj,
    dwten_p: cobj,
    diten_p: cobj,
    flxrain_p: cobj,
    flxsnow_p: cobj,
    ntraprd_p: cobj,
    ntsnprd_p: cobj,
    excessu_p: cobj,
    excess0_p: cobj,
    xc_p: cobj,
    aquad_p: cobj,
    bquad_p: cobj,
    cquad_p: cobj,
    bogbot_p: cobj,
    bogtop_p: cobj,
    exit_uwcu_p: cobj,
    exit_conden_p: cobj,
    exit_klclmkx_p: cobj,
    exit_klfcmkx_p: cobj,
    exit_ufrc_p: cobj,
    exit_wtw_p: cobj,
    exit_drycore_p: cobj,
    exit_wu_p: cobj,
    exit_cufilter_p: cobj,
    exit_kinv1_p: cobj,
    exit_rei_p: cobj,
    limit_shcu_p: cobj,
    limit_negcon_p: cobj,
    limit_ufrc_p: cobj,
    limit_ppen_p: cobj,
    limit_emf_p: cobj,
    limit_cinlcl_p: cobj,
    limit_cin_p: cobj,
    limit_cbmf_p: cobj,
    limit_rei_p: cobj,
    ind_delcin_p: cobj,
):
    ufrcinvbase = Ptr[float](ufrcinvbase_p)
    ufrclcl = Ptr[float](ufrclcl_p)
    winvbase = Ptr[float](winvbase_p)
    wlcl = Ptr[float](wlcl_p)
    plcl = Ptr[float](plcl_p)
    pinv = Ptr[float](pinv_p)
    plfc = Ptr[float](plfc_p)
    pbup = Ptr[float](pbup_p)
    ppen = Ptr[float](ppen_p)
    qtsrc = Ptr[float](qtsrc_p)
    thlsrc = Ptr[float](thlsrc_p)
    thvlsrc = Ptr[float](thvlsrc_p)
    emfkbup = Ptr[float](emfkbup_p)
    cbmflimit = Ptr[float](cbmflimit_p)
    tkeavg = Ptr[float](tkeavg_p)
    zinv = Ptr[float](zinv_p)
    rcwp = Ptr[float](rcwp_p)
    rlwp = Ptr[float](rlwp_p)
    riwp = Ptr[float](riwp_p)
    wu = Ptr[float](wu_p)
    qtu = Ptr[float](qtu_p)
    thlu = Ptr[float](thlu_p)
    thvu = Ptr[float](thvu_p)
    uu = Ptr[float](uu_p)
    vu = Ptr[float](vu_p)
    qtu_emf = Ptr[float](qtu_emf_p)
    thlu_emf = Ptr[float](thlu_emf_p)
    uu_emf = Ptr[float](uu_emf_p)
    vu_emf = Ptr[float](vu_emf_p)
    uemf = Ptr[float](uemf_p)
    tru = Ptr[float](tru_p)
    tru_emf = Ptr[float](tru_emf_p)
    dwten = Ptr[float](dwten_p)
    diten = Ptr[float](diten_p)
    flxrain = Ptr[float](flxrain_p)
    flxsnow = Ptr[float](flxsnow_p)
    ntraprd = Ptr[float](ntraprd_p)
    ntsnprd = Ptr[float](ntsnprd_p)
    excessu = Ptr[float](excessu_p)
    excess0 = Ptr[float](excess0_p)
    xc = Ptr[float](xc_p)
    aquad = Ptr[float](aquad_p)
    bquad = Ptr[float](bquad_p)
    cquad = Ptr[float](cquad_p)
    bogbot = Ptr[float](bogbot_p)
    bogtop = Ptr[float](bogtop_p)
    exit_uwcu = Ptr[float](exit_uwcu_p)
    exit_conden = Ptr[float](exit_conden_p)
    exit_klclmkx = Ptr[float](exit_klclmkx_p)
    exit_klfcmkx = Ptr[float](exit_klfcmkx_p)
    exit_ufrc = Ptr[float](exit_ufrc_p)
    exit_wtw = Ptr[float](exit_wtw_p)
    exit_drycore = Ptr[float](exit_drycore_p)
    exit_wu = Ptr[float](exit_wu_p)
    exit_cufilter = Ptr[float](exit_cufilter_p)
    exit_kinv1 = Ptr[float](exit_kinv1_p)
    exit_rei = Ptr[float](exit_rei_p)
    limit_shcu = Ptr[float](limit_shcu_p)
    limit_negcon = Ptr[float](limit_negcon_p)
    limit_ufrc = Ptr[float](limit_ufrc_p)
    limit_ppen = Ptr[float](limit_ppen_p)
    limit_emf = Ptr[float](limit_emf_p)
    limit_cinlcl = Ptr[float](limit_cinlcl_p)
    limit_cin = Ptr[float](limit_cin_p)
    limit_cbmf = Ptr[float](limit_cbmf_p)
    limit_rei = Ptr[float](limit_rei_p)
    ind_delcin = Ptr[float](ind_delcin_p)

    i = 0
    while i < iend:
        ufrcinvbase[i] = 0.0
        ufrclcl[i] = 0.0
        winvbase[i] = 0.0
        wlcl[i] = 0.0
        plcl[i] = 0.0
        pinv[i] = 0.0
        plfc[i] = 0.0
        pbup[i] = 0.0
        ppen[i] = 0.0
        qtsrc[i] = 0.0
        thlsrc[i] = 0.0
        thvlsrc[i] = 0.0
        emfkbup[i] = 0.0
        cbmflimit[i] = 0.0
        tkeavg[i] = 0.0
        zinv[i] = 0.0
        rcwp[i] = 0.0
        rlwp[i] = 0.0
        riwp[i] = 0.0
        exit_uwcu[i] = 0.0
        exit_conden[i] = 0.0
        exit_klclmkx[i] = 0.0
        exit_klfcmkx[i] = 0.0
        exit_ufrc[i] = 0.0
        exit_wtw[i] = 0.0
        exit_drycore[i] = 0.0
        exit_wu[i] = 0.0
        exit_cufilter[i] = 0.0
        exit_kinv1[i] = 0.0
        exit_rei[i] = 0.0
        limit_shcu[i] = 0.0
        limit_negcon[i] = 0.0
        limit_ufrc[i] = 0.0
        limit_ppen[i] = 0.0
        limit_emf[i] = 0.0
        limit_cinlcl[i] = 0.0
        limit_cin[i] = 0.0
        limit_cbmf[i] = 0.0
        limit_rei[i] = 0.0
        ind_delcin[i] = 0.0
        i += 1

    k = 0
    while k <= mkx:
        i = 0
        while i < iend:
            idx = i + k * mix
            wu[idx] = 0.0
            qtu[idx] = 0.0
            thlu[idx] = 0.0
            thvu[idx] = 0.0
            uu[idx] = 0.0
            vu[idx] = 0.0
            qtu_emf[idx] = 0.0
            thlu_emf[idx] = 0.0
            uu_emf[idx] = 0.0
            vu_emf[idx] = 0.0
            uemf[idx] = 0.0
            flxrain[idx] = 0.0
            flxsnow[idx] = 0.0
            i += 1
        k += 1

    k = 0
    while k < mkx:
        i = 0
        while i < iend:
            idx = i + k * mix
            dwten[idx] = 0.0
            diten[idx] = 0.0
            excessu[idx] = 0.0
            excess0[idx] = 0.0
            xc[idx] = 0.0
            aquad[idx] = 0.0
            bquad[idx] = 0.0
            cquad[idx] = 0.0
            bogbot[idx] = 0.0
            bogtop[idx] = 0.0
            i += 1
        k += 1

    last_k_idx = mkx - 1
    i = 0
    while i < iend:
        ntraprd[i + last_k_idx * mix] = 0.0
        ntsnprd[i + last_k_idx * mix] = 0.0
        i += 1

    m = 0
    while m < ncnst:
        offset = m * mix * (mkx + 1)
        k = 0
        while k <= mkx:
            i = 0
            while i < iend:
                idx = i + k * mix + offset
                tru[idx] = 0.0
                tru_emf[idx] = 0.0
                i += 1
            k += 1
        m += 1


@export
def uwshcu_output_diag_init_shell_codon(
    mix: int,
    mkx: int,
    iend: int,
    ncnst: int,
    umf_p: cobj,
    slflx_p: cobj,
    qtflx_p: cobj,
    flxprc1_p: cobj,
    flxsnow1_p: cobj,
    qvten_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    sten_p: cobj,
    uten_p: cobj,
    vten_p: cobj,
    qrten_p: cobj,
    qsten_p: cobj,
    evapc_p: cobj,
    cufrc_p: cobj,
    qcu_p: cobj,
    qlu_p: cobj,
    qiu_p: cobj,
    fer_p: cobj,
    fdr_p: cobj,
    qc_p: cobj,
    qtten_p: cobj,
    slten_p: cobj,
    ufrc_p: cobj,
    uflx_p: cobj,
    vflx_p: cobj,
    trten_p: cobj,
    trflx_p: cobj,
    wtqc_p: cobj,
    wtprec_p: cobj,
    wtsnow_p: cobj,
    precip_p: cobj,
    snow_p: cobj,
    cinh_p: cobj,
    cinlclh_p: cobj,
    cbmf_p: cobj,
    rliq_p: cobj,
    cnt_p: cobj,
    cnb_p: cobj,
    ufrcinvbase_p: cobj,
    ufrclcl_p: cobj,
    winvbase_p: cobj,
    wlcl_p: cobj,
    plcl_p: cobj,
    pinv_p: cobj,
    plfc_p: cobj,
    pbup_p: cobj,
    ppen_p: cobj,
    qtsrc_p: cobj,
    thlsrc_p: cobj,
    thvlsrc_p: cobj,
    emfkbup_p: cobj,
    cbmflimit_p: cobj,
    tkeavg_p: cobj,
    zinv_p: cobj,
    rcwp_p: cobj,
    rlwp_p: cobj,
    riwp_p: cobj,
    wu_p: cobj,
    qtu_p: cobj,
    thlu_p: cobj,
    thvu_p: cobj,
    uu_p: cobj,
    vu_p: cobj,
    qtu_emf_p: cobj,
    thlu_emf_p: cobj,
    uu_emf_p: cobj,
    vu_emf_p: cobj,
    uemf_p: cobj,
    tru_p: cobj,
    tru_emf_p: cobj,
    dwten_p: cobj,
    diten_p: cobj,
    flxrain_p: cobj,
    flxsnow_p: cobj,
    ntraprd_p: cobj,
    ntsnprd_p: cobj,
    excessu_p: cobj,
    excess0_p: cobj,
    xc_p: cobj,
    aquad_p: cobj,
    bquad_p: cobj,
    cquad_p: cobj,
    bogbot_p: cobj,
    bogtop_p: cobj,
    exit_uwcu_p: cobj,
    exit_conden_p: cobj,
    exit_klclmkx_p: cobj,
    exit_klfcmkx_p: cobj,
    exit_ufrc_p: cobj,
    exit_wtw_p: cobj,
    exit_drycore_p: cobj,
    exit_wu_p: cobj,
    exit_cufilter_p: cobj,
    exit_kinv1_p: cobj,
    exit_rei_p: cobj,
    limit_shcu_p: cobj,
    limit_negcon_p: cobj,
    limit_ufrc_p: cobj,
    limit_ppen_p: cobj,
    limit_emf_p: cobj,
    limit_cinlcl_p: cobj,
    limit_cin_p: cobj,
    limit_cbmf_p: cobj,
    limit_rei_p: cobj,
    ind_delcin_p: cobj,
):
    uwshcu_output_init_shell_codon(
        mix,
        mkx,
        iend,
        ncnst,
        umf_p,
        slflx_p,
        qtflx_p,
        flxprc1_p,
        flxsnow1_p,
        qvten_p,
        qlten_p,
        qiten_p,
        sten_p,
        uten_p,
        vten_p,
        qrten_p,
        qsten_p,
        evapc_p,
        cufrc_p,
        qcu_p,
        qlu_p,
        qiu_p,
        fer_p,
        fdr_p,
        qc_p,
        qtten_p,
        slten_p,
        ufrc_p,
        uflx_p,
        vflx_p,
        trten_p,
        trflx_p,
        wtqc_p,
        wtprec_p,
        wtsnow_p,
        precip_p,
        snow_p,
        cinh_p,
        cinlclh_p,
        cbmf_p,
        rliq_p,
        cnt_p,
        cnb_p,
    )
    uwshcu_diag_init_shell_codon(
        mix,
        mkx,
        iend,
        ncnst,
        ufrcinvbase_p,
        ufrclcl_p,
        winvbase_p,
        wlcl_p,
        plcl_p,
        pinv_p,
        plfc_p,
        pbup_p,
        ppen_p,
        qtsrc_p,
        thlsrc_p,
        thvlsrc_p,
        emfkbup_p,
        cbmflimit_p,
        tkeavg_p,
        zinv_p,
        rcwp_p,
        rlwp_p,
        riwp_p,
        wu_p,
        qtu_p,
        thlu_p,
        thvu_p,
        uu_p,
        vu_p,
        qtu_emf_p,
        thlu_emf_p,
        uu_emf_p,
        vu_emf_p,
        uemf_p,
        tru_p,
        tru_emf_p,
        dwten_p,
        diten_p,
        flxrain_p,
        flxsnow_p,
        ntraprd_p,
        ntsnprd_p,
        excessu_p,
        excess0_p,
        xc_p,
        aquad_p,
        bquad_p,
        cquad_p,
        bogbot_p,
        bogtop_p,
        exit_uwcu_p,
        exit_conden_p,
        exit_klclmkx_p,
        exit_klfcmkx_p,
        exit_ufrc_p,
        exit_wtw_p,
        exit_drycore_p,
        exit_wu_p,
        exit_cufilter_p,
        exit_kinv1_p,
        exit_rei_p,
        limit_shcu_p,
        limit_negcon_p,
        limit_ufrc_p,
        limit_ppen_p,
        limit_emf_p,
        limit_cinlcl_p,
        limit_cin_p,
        limit_cbmf_p,
        limit_rei_p,
        ind_delcin_p,
    )


@export
def uwshcu_diag_post_shell_codon(
    mix: int,
    mkx: int,
    i_col: int,
    ncnst: int,
    cin: float,
    cinlcl: float,
    ufrcinvbase_v: float,
    ufrclcl_v: float,
    winvbase_v: float,
    wlcl_v: float,
    plcl_v: float,
    pinv_v: float,
    plfc_v: float,
    pbup_v: float,
    ppen_v: float,
    qtsrc_v: float,
    thlsrc_v: float,
    thvlsrc_v: float,
    emfkbup_v: float,
    cbmflimit_v: float,
    tkeavg_v: float,
    zinv_v: float,
    rcwp_v: float,
    rlwp_v: float,
    riwp_v: float,
    fer_p: cobj,
    fdr_p: cobj,
    qtten_p: cobj,
    slten_p: cobj,
    ufrc_p: cobj,
    uflx_p: cobj,
    vflx_p: cobj,
    wu_p: cobj,
    qtu_p: cobj,
    thlu_p: cobj,
    thvu_p: cobj,
    uu_p: cobj,
    vu_p: cobj,
    qtu_emf_p: cobj,
    thlu_emf_p: cobj,
    uu_emf_p: cobj,
    vu_emf_p: cobj,
    uemf_p: cobj,
    dwten_p: cobj,
    diten_p: cobj,
    flxrain_p: cobj,
    flxsnow_p: cobj,
    ntraprd_p: cobj,
    ntsnprd_p: cobj,
    excessu_p: cobj,
    excess0_p: cobj,
    xc_p: cobj,
    aquad_p: cobj,
    bquad_p: cobj,
    cquad_p: cobj,
    bogbot_p: cobj,
    bogtop_p: cobj,
    trflx_p: cobj,
    tru_p: cobj,
    tru_emf_p: cobj,
    fer_out_p: cobj,
    fdr_out_p: cobj,
    cinh_out_p: cobj,
    cinlclh_out_p: cobj,
    qtten_out_p: cobj,
    slten_out_p: cobj,
    ufrc_out_p: cobj,
    uflx_out_p: cobj,
    vflx_out_p: cobj,
    ufrcinvbase_out_p: cobj,
    ufrclcl_out_p: cobj,
    winvbase_out_p: cobj,
    wlcl_out_p: cobj,
    plcl_out_p: cobj,
    pinv_out_p: cobj,
    plfc_out_p: cobj,
    pbup_out_p: cobj,
    ppen_out_p: cobj,
    qtsrc_out_p: cobj,
    thlsrc_out_p: cobj,
    thvlsrc_out_p: cobj,
    emfkbup_out_p: cobj,
    cbmflimit_out_p: cobj,
    tkeavg_out_p: cobj,
    zinv_out_p: cobj,
    rcwp_out_p: cobj,
    rlwp_out_p: cobj,
    riwp_out_p: cobj,
    wu_out_p: cobj,
    qtu_out_p: cobj,
    thlu_out_p: cobj,
    thvu_out_p: cobj,
    uu_out_p: cobj,
    vu_out_p: cobj,
    qtu_emf_out_p: cobj,
    thlu_emf_out_p: cobj,
    uu_emf_out_p: cobj,
    vu_emf_out_p: cobj,
    uemf_out_p: cobj,
    dwten_out_p: cobj,
    diten_out_p: cobj,
    flxrain_out_p: cobj,
    flxsnow_out_p: cobj,
    ntraprd_out_p: cobj,
    ntsnprd_out_p: cobj,
    excessu_out_p: cobj,
    excess0_out_p: cobj,
    xc_out_p: cobj,
    aquad_out_p: cobj,
    bquad_out_p: cobj,
    cquad_out_p: cobj,
    bogbot_out_p: cobj,
    bogtop_out_p: cobj,
    trflx_out_p: cobj,
    tru_out_p: cobj,
    tru_emf_out_p: cobj,
):
    fer = Ptr[float](fer_p)
    fdr = Ptr[float](fdr_p)
    qtten = Ptr[float](qtten_p)
    slten = Ptr[float](slten_p)
    ufrc = Ptr[float](ufrc_p)
    uflx = Ptr[float](uflx_p)
    vflx = Ptr[float](vflx_p)
    wu = Ptr[float](wu_p)
    qtu = Ptr[float](qtu_p)
    thlu = Ptr[float](thlu_p)
    thvu = Ptr[float](thvu_p)
    uu = Ptr[float](uu_p)
    vu = Ptr[float](vu_p)
    qtu_emf = Ptr[float](qtu_emf_p)
    thlu_emf = Ptr[float](thlu_emf_p)
    uu_emf = Ptr[float](uu_emf_p)
    vu_emf = Ptr[float](vu_emf_p)
    uemf = Ptr[float](uemf_p)
    dwten = Ptr[float](dwten_p)
    diten = Ptr[float](diten_p)
    flxrain = Ptr[float](flxrain_p)
    flxsnow = Ptr[float](flxsnow_p)
    ntraprd = Ptr[float](ntraprd_p)
    ntsnprd = Ptr[float](ntsnprd_p)
    excessu = Ptr[float](excessu_p)
    excess0 = Ptr[float](excess0_p)
    xc = Ptr[float](xc_p)
    aquad = Ptr[float](aquad_p)
    bquad = Ptr[float](bquad_p)
    cquad = Ptr[float](cquad_p)
    bogbot = Ptr[float](bogbot_p)
    bogtop = Ptr[float](bogtop_p)
    trflx = Ptr[float](trflx_p)
    tru = Ptr[float](tru_p)
    tru_emf = Ptr[float](tru_emf_p)
    fer_out = Ptr[float](fer_out_p)
    fdr_out = Ptr[float](fdr_out_p)
    cinh_out = Ptr[float](cinh_out_p)
    cinlclh_out = Ptr[float](cinlclh_out_p)
    qtten_out = Ptr[float](qtten_out_p)
    slten_out = Ptr[float](slten_out_p)
    ufrc_out = Ptr[float](ufrc_out_p)
    uflx_out = Ptr[float](uflx_out_p)
    vflx_out = Ptr[float](vflx_out_p)
    ufrcinvbase_out = Ptr[float](ufrcinvbase_out_p)
    ufrclcl_out = Ptr[float](ufrclcl_out_p)
    winvbase_out = Ptr[float](winvbase_out_p)
    wlcl_out = Ptr[float](wlcl_out_p)
    plcl_out = Ptr[float](plcl_out_p)
    pinv_out = Ptr[float](pinv_out_p)
    plfc_out = Ptr[float](plfc_out_p)
    pbup_out = Ptr[float](pbup_out_p)
    ppen_out = Ptr[float](ppen_out_p)
    qtsrc_out = Ptr[float](qtsrc_out_p)
    thlsrc_out = Ptr[float](thlsrc_out_p)
    thvlsrc_out = Ptr[float](thvlsrc_out_p)
    emfkbup_out = Ptr[float](emfkbup_out_p)
    cbmflimit_out = Ptr[float](cbmflimit_out_p)
    tkeavg_out = Ptr[float](tkeavg_out_p)
    zinv_out = Ptr[float](zinv_out_p)
    rcwp_out = Ptr[float](rcwp_out_p)
    rlwp_out = Ptr[float](rlwp_out_p)
    riwp_out = Ptr[float](riwp_out_p)
    wu_out = Ptr[float](wu_out_p)
    qtu_out = Ptr[float](qtu_out_p)
    thlu_out = Ptr[float](thlu_out_p)
    thvu_out = Ptr[float](thvu_out_p)
    uu_out = Ptr[float](uu_out_p)
    vu_out = Ptr[float](vu_out_p)
    qtu_emf_out = Ptr[float](qtu_emf_out_p)
    thlu_emf_out = Ptr[float](thlu_emf_out_p)
    uu_emf_out = Ptr[float](uu_emf_out_p)
    vu_emf_out = Ptr[float](vu_emf_out_p)
    uemf_out = Ptr[float](uemf_out_p)
    dwten_out = Ptr[float](dwten_out_p)
    diten_out = Ptr[float](diten_out_p)
    flxrain_out = Ptr[float](flxrain_out_p)
    flxsnow_out = Ptr[float](flxsnow_out_p)
    ntraprd_out = Ptr[float](ntraprd_out_p)
    ntsnprd_out = Ptr[float](ntsnprd_out_p)
    excessu_out = Ptr[float](excessu_out_p)
    excess0_out = Ptr[float](excess0_out_p)
    xc_out = Ptr[float](xc_out_p)
    aquad_out = Ptr[float](aquad_out_p)
    bquad_out = Ptr[float](bquad_out_p)
    cquad_out = Ptr[float](cquad_out_p)
    bogbot_out = Ptr[float](bogbot_out_p)
    bogtop_out = Ptr[float](bogtop_out_p)
    trflx_out = Ptr[float](trflx_out_p)
    tru_out = Ptr[float](tru_out_p)
    tru_emf_out = Ptr[float](tru_emf_out_p)

    col = i_col - 1
    cinh_out[col] = cin
    cinlclh_out[col] = cinlcl
    ufrcinvbase_out[col] = ufrcinvbase_v
    ufrclcl_out[col] = ufrclcl_v
    winvbase_out[col] = winvbase_v
    wlcl_out[col] = wlcl_v
    plcl_out[col] = plcl_v
    pinv_out[col] = pinv_v
    plfc_out[col] = plfc_v
    pbup_out[col] = pbup_v
    ppen_out[col] = ppen_v
    qtsrc_out[col] = qtsrc_v
    thlsrc_out[col] = thlsrc_v
    thvlsrc_out[col] = thvlsrc_v
    emfkbup_out[col] = emfkbup_v
    cbmflimit_out[col] = cbmflimit_v
    tkeavg_out[col] = tkeavg_v
    zinv_out[col] = zinv_v
    rcwp_out[col] = rcwp_v
    rlwp_out[col] = rlwp_v
    riwp_out[col] = riwp_v

    k = 0
    while k < mkx:
        dst_k = mkx - k - 1
        dst = col + dst_k * mix
        fer_out[dst] = fer[k]
        fdr_out[dst] = fdr[k]
        qtten_out[dst] = qtten[k]
        slten_out[dst] = slten[k]
        dwten_out[dst] = dwten[k]
        diten_out[dst] = diten[k]
        ntraprd_out[dst] = ntraprd[k]
        ntsnprd_out[dst] = ntsnprd[k]
        excessu_out[dst] = excessu[k]
        excess0_out[dst] = excess0[k]
        xc_out[dst] = xc[k]
        aquad_out[dst] = aquad[k]
        bquad_out[dst] = bquad[k]
        cquad_out[dst] = cquad[k]
        bogbot_out[dst] = bogbot[k]
        bogtop_out[dst] = bogtop[k]
        k += 1

    k = 0
    while k <= mkx:
        dst_k = mkx - k
        dst = col + dst_k * mix
        ufrc_out[dst] = ufrc[k]
        uflx_out[dst] = uflx[k]
        vflx_out[dst] = vflx[k]
        wu_out[dst] = wu[k]
        qtu_out[dst] = qtu[k]
        thlu_out[dst] = thlu[k]
        thvu_out[dst] = thvu[k]
        uu_out[dst] = uu[k]
        vu_out[dst] = vu[k]
        qtu_emf_out[dst] = qtu_emf[k]
        thlu_emf_out[dst] = thlu_emf[k]
        uu_emf_out[dst] = uu_emf[k]
        vu_emf_out[dst] = vu_emf[k]
        uemf_out[dst] = uemf[k]
        flxrain_out[dst] = flxrain[k]
        flxsnow_out[dst] = flxsnow[k]
        k += 1

    m = 0
    while m < ncnst:
        offset = m * mix * (mkx + 1)
        k = 0
        while k <= mkx:
            dst_k = mkx - k
            dst = col + dst_k * mix + offset
            src = k + m * (mkx + 1)
            trflx_out[dst] = trflx[src]
            tru_out[dst] = tru[src]
            tru_emf_out[dst] = tru_emf[src]
            k += 1
        m += 1


@export
def uwshcu_main_post_shell_codon(
    mix: int,
    mkx: int,
    i_col: int,
    ncnst: int,
    precip_v: float,
    snow_v: float,
    cush_v: float,
    cbmf_v: float,
    rliq_v: float,
    cnt_v: float,
    cnb_v: float,
    umf_p: cobj,
    slflx_p: cobj,
    qtflx_p: cobj,
    flxrain_p: cobj,
    flxsnow_p: cobj,
    qvten_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    sten_p: cobj,
    uten_p: cobj,
    vten_p: cobj,
    qrten_p: cobj,
    qsten_p: cobj,
    evapc_p: cobj,
    cufrc_p: cobj,
    qcu_p: cobj,
    qlu_p: cobj,
    qiu_p: cobj,
    qc_p: cobj,
    trten_p: cobj,
    umf_out_p: cobj,
    slflx_out_p: cobj,
    qtflx_out_p: cobj,
    flxprc1_out_p: cobj,
    flxsnow1_out_p: cobj,
    qvten_out_p: cobj,
    qlten_out_p: cobj,
    qiten_out_p: cobj,
    sten_out_p: cobj,
    uten_out_p: cobj,
    vten_out_p: cobj,
    qrten_out_p: cobj,
    qsten_out_p: cobj,
    precip_out_p: cobj,
    snow_out_p: cobj,
    evapc_out_p: cobj,
    cufrc_out_p: cobj,
    qcu_out_p: cobj,
    qlu_out_p: cobj,
    qiu_out_p: cobj,
    cush_out_p: cobj,
    cbmf_out_p: cobj,
    rliq_out_p: cobj,
    qc_out_p: cobj,
    cnt_out_p: cobj,
    cnb_out_p: cobj,
    trten_out_p: cobj,
):
    umf = Ptr[float](umf_p)
    slflx = Ptr[float](slflx_p)
    qtflx = Ptr[float](qtflx_p)
    flxrain = Ptr[float](flxrain_p)
    flxsnow = Ptr[float](flxsnow_p)
    qvten = Ptr[float](qvten_p)
    qlten = Ptr[float](qlten_p)
    qiten = Ptr[float](qiten_p)
    sten = Ptr[float](sten_p)
    uten = Ptr[float](uten_p)
    vten = Ptr[float](vten_p)
    qrten = Ptr[float](qrten_p)
    qsten = Ptr[float](qsten_p)
    evapc = Ptr[float](evapc_p)
    cufrc = Ptr[float](cufrc_p)
    qcu = Ptr[float](qcu_p)
    qlu = Ptr[float](qlu_p)
    qiu = Ptr[float](qiu_p)
    qc = Ptr[float](qc_p)
    trten = Ptr[float](trten_p)
    umf_out = Ptr[float](umf_out_p)
    slflx_out = Ptr[float](slflx_out_p)
    qtflx_out = Ptr[float](qtflx_out_p)
    flxprc1_out = Ptr[float](flxprc1_out_p)
    flxsnow1_out = Ptr[float](flxsnow1_out_p)
    qvten_out = Ptr[float](qvten_out_p)
    qlten_out = Ptr[float](qlten_out_p)
    qiten_out = Ptr[float](qiten_out_p)
    sten_out = Ptr[float](sten_out_p)
    uten_out = Ptr[float](uten_out_p)
    vten_out = Ptr[float](vten_out_p)
    qrten_out = Ptr[float](qrten_out_p)
    qsten_out = Ptr[float](qsten_out_p)
    precip_out = Ptr[float](precip_out_p)
    snow_out = Ptr[float](snow_out_p)
    evapc_out = Ptr[float](evapc_out_p)
    cufrc_out = Ptr[float](cufrc_out_p)
    qcu_out = Ptr[float](qcu_out_p)
    qlu_out = Ptr[float](qlu_out_p)
    qiu_out = Ptr[float](qiu_out_p)
    cush_out = Ptr[float](cush_out_p)
    cbmf_out = Ptr[float](cbmf_out_p)
    rliq_out = Ptr[float](rliq_out_p)
    qc_out = Ptr[float](qc_out_p)
    cnt_out = Ptr[float](cnt_out_p)
    cnb_out = Ptr[float](cnb_out_p)
    trten_out = Ptr[float](trten_out_p)

    col = i_col - 1
    precip_out[col] = precip_v
    snow_out[col] = snow_v
    cush_out[col] = cush_v
    cbmf_out[col] = cbmf_v
    rliq_out[col] = rliq_v
    cnt_out[col] = cnt_v
    cnb_out[col] = cnb_v

    k = 0
    while k <= mkx:
        dst = col + k * mix
        umf_out[dst] = umf[k]
        slflx_out[dst] = slflx[k]
        qtflx_out[dst] = qtflx[k]
        flxprc1_out[dst] = flxrain[k] + flxsnow[k]
        flxsnow1_out[dst] = flxsnow[k]
        k += 1

    k = 0
    while k < mkx:
        dst = col + k * mix
        qvten_out[dst] = qvten[k]
        qlten_out[dst] = qlten[k]
        qiten_out[dst] = qiten[k]
        sten_out[dst] = sten[k]
        uten_out[dst] = uten[k]
        vten_out[dst] = vten[k]
        qrten_out[dst] = qrten[k]
        qsten_out[dst] = qsten[k]
        evapc_out[dst] = evapc[k]
        cufrc_out[dst] = cufrc[k]
        qcu_out[dst] = qcu[k]
        qlu_out[dst] = qlu[k]
        qiu_out[dst] = qiu[k]
        qc_out[dst] = qc[k]
        k += 1

    m = 0
    while m < ncnst:
        offset = m * mix * mkx
        src_offset = m * mkx
        k = 0
        while k < mkx:
            trten_out[col + k * mix + offset] = trten[k + src_offset]
            k += 1
        m += 1


@export
def uwshcu_wtrc_post_shell_codon(
    mix: int,
    mkx: int,
    i_col: int,
    ncnst: int,
    wtrc_nwset: int,
    wtqc_liq_p: cobj,
    wtqc_ice_p: cobj,
    wtprec_p: cobj,
    wtsnow_p: cobj,
    wtrc_iatype_p: cobj,
    wtqc_out_p: cobj,
    wtprec_out_p: cobj,
    wtsnow_out_p: cobj,
):
    wtqc_liq = Ptr[float](wtqc_liq_p)
    wtqc_ice = Ptr[float](wtqc_ice_p)
    wtprec = Ptr[float](wtprec_p)
    wtsnow = Ptr[float](wtsnow_p)
    wtrc_iatype = Ptr[int](wtrc_iatype_p)
    wtqc_out = Ptr[float](wtqc_out_p)
    wtprec_out = Ptr[float](wtprec_out_p)
    wtsnow_out = Ptr[float](wtsnow_out_p)

    col = i_col - 1
    m = 0
    while m < wtrc_nwset:
        vap = wtrc_iatype[m] - 1
        liq = wtrc_iatype[m + wtrc_nwset] - 1
        ice = wtrc_iatype[m + 2 * wtrc_nwset] - 1
        src_offset = m * mkx
        liq_offset = liq * mix * mkx
        ice_offset = ice * mix * mkx
        k = 0
        while k < mkx:
            wtqc_out[col + k * mix + liq_offset] = wtqc_liq[k + src_offset]
            wtqc_out[col + k * mix + ice_offset] = wtqc_ice[k + src_offset]
            k += 1
        wtprec_out[col + vap * mix] = wtprec[m]
        wtsnow_out[col + vap * mix] = wtsnow[m]
        m += 1


@export
def uwshcu_main_wtrc_post_shell_codon(
    mix: int,
    mkx: int,
    i_col: int,
    ncnst: int,
    wtrc_nwset: int,
    precip_v: float,
    snow_v: float,
    cush_v: float,
    cbmf_v: float,
    rliq_v: float,
    cnt_v: float,
    cnb_v: float,
    umf_p: cobj,
    slflx_p: cobj,
    qtflx_p: cobj,
    flxrain_p: cobj,
    flxsnow_p: cobj,
    qvten_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    sten_p: cobj,
    uten_p: cobj,
    vten_p: cobj,
    qrten_p: cobj,
    qsten_p: cobj,
    evapc_p: cobj,
    cufrc_p: cobj,
    qcu_p: cobj,
    qlu_p: cobj,
    qiu_p: cobj,
    qc_p: cobj,
    trten_p: cobj,
    wtqc_liq_p: cobj,
    wtqc_ice_p: cobj,
    wtprec_p: cobj,
    wtsnow_p: cobj,
    wtrc_iatype_p: cobj,
    umf_out_p: cobj,
    slflx_out_p: cobj,
    qtflx_out_p: cobj,
    flxprc1_out_p: cobj,
    flxsnow1_out_p: cobj,
    qvten_out_p: cobj,
    qlten_out_p: cobj,
    qiten_out_p: cobj,
    sten_out_p: cobj,
    uten_out_p: cobj,
    vten_out_p: cobj,
    qrten_out_p: cobj,
    qsten_out_p: cobj,
    precip_out_p: cobj,
    snow_out_p: cobj,
    evapc_out_p: cobj,
    cufrc_out_p: cobj,
    qcu_out_p: cobj,
    qlu_out_p: cobj,
    qiu_out_p: cobj,
    cush_out_p: cobj,
    cbmf_out_p: cobj,
    rliq_out_p: cobj,
    qc_out_p: cobj,
    cnt_out_p: cobj,
    cnb_out_p: cobj,
    trten_out_p: cobj,
    wtqc_out_p: cobj,
    wtprec_out_p: cobj,
    wtsnow_out_p: cobj,
):
    uwshcu_main_post_shell_codon(
        mix,
        mkx,
        i_col,
        ncnst,
        precip_v,
        snow_v,
        cush_v,
        cbmf_v,
        rliq_v,
        cnt_v,
        cnb_v,
        umf_p,
        slflx_p,
        qtflx_p,
        flxrain_p,
        flxsnow_p,
        qvten_p,
        qlten_p,
        qiten_p,
        sten_p,
        uten_p,
        vten_p,
        qrten_p,
        qsten_p,
        evapc_p,
        cufrc_p,
        qcu_p,
        qlu_p,
        qiu_p,
        qc_p,
        trten_p,
        umf_out_p,
        slflx_out_p,
        qtflx_out_p,
        flxprc1_out_p,
        flxsnow1_out_p,
        qvten_out_p,
        qlten_out_p,
        qiten_out_p,
        sten_out_p,
        uten_out_p,
        vten_out_p,
        qrten_out_p,
        qsten_out_p,
        precip_out_p,
        snow_out_p,
        evapc_out_p,
        cufrc_out_p,
        qcu_out_p,
        qlu_out_p,
        qiu_out_p,
        cush_out_p,
        cbmf_out_p,
        rliq_out_p,
        qc_out_p,
        cnt_out_p,
        cnb_out_p,
        trten_out_p,
    )
    if wtrc_nwset > 0:
        uwshcu_wtrc_post_shell_codon(
            mix,
            mkx,
            i_col,
            ncnst,
            wtrc_nwset,
            wtqc_liq_p,
            wtqc_ice_p,
            wtprec_p,
            wtsnow_p,
            wtrc_iatype_p,
            wtqc_out_p,
            wtprec_out_p,
            wtsnow_out_p,
        )


@export
def uwshcu_main_diag_post_all_shell_codon(
    mix: int,
    mkx: int,
    i_col: int,
    ncnst: int,
    wtrc_nwset: int,
    precip_v: float,
    snow_v: float,
    cush_v: float,
    cbmf_v: float,
    rliq_v: float,
    cnt_v: float,
    cnb_v: float,
    umf_p: cobj,
    slflx_p: cobj,
    qtflx_p: cobj,
    flxrain_p: cobj,
    flxsnow_p: cobj,
    qvten_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    sten_p: cobj,
    uten_p: cobj,
    vten_p: cobj,
    qrten_p: cobj,
    qsten_p: cobj,
    evapc_p: cobj,
    cufrc_p: cobj,
    qcu_p: cobj,
    qlu_p: cobj,
    qiu_p: cobj,
    qc_p: cobj,
    trten_p: cobj,
    wtqc_liq_p: cobj,
    wtqc_ice_p: cobj,
    wtprec_p: cobj,
    wtsnow_p: cobj,
    wtrc_iatype_p: cobj,
    umf_out_p: cobj,
    slflx_out_p: cobj,
    qtflx_out_p: cobj,
    flxprc1_out_p: cobj,
    flxsnow1_out_p: cobj,
    qvten_out_p: cobj,
    qlten_out_p: cobj,
    qiten_out_p: cobj,
    sten_out_p: cobj,
    uten_out_p: cobj,
    vten_out_p: cobj,
    qrten_out_p: cobj,
    qsten_out_p: cobj,
    precip_out_p: cobj,
    snow_out_p: cobj,
    evapc_out_p: cobj,
    cufrc_out_p: cobj,
    qcu_out_p: cobj,
    qlu_out_p: cobj,
    qiu_out_p: cobj,
    cush_out_p: cobj,
    cbmf_out_p: cobj,
    rliq_out_p: cobj,
    qc_out_p: cobj,
    cnt_out_p: cobj,
    cnb_out_p: cobj,
    trten_out_p: cobj,
    wtqc_out_p: cobj,
    wtprec_out_p: cobj,
    wtsnow_out_p: cobj,
    cin: float,
    cinlcl: float,
    ufrcinvbase_v: float,
    ufrclcl_v: float,
    winvbase_v: float,
    wlcl_v: float,
    plcl_v: float,
    pinv_v: float,
    plfc_v: float,
    pbup_v: float,
    ppen_v: float,
    qtsrc_v: float,
    thlsrc_v: float,
    thvlsrc_v: float,
    emfkbup_v: float,
    cbmflimit_v: float,
    tkeavg_v: float,
    zinv_v: float,
    rcwp_v: float,
    rlwp_v: float,
    riwp_v: float,
    fer_p: cobj,
    fdr_p: cobj,
    qtten_p: cobj,
    slten_p: cobj,
    ufrc_p: cobj,
    uflx_p: cobj,
    vflx_p: cobj,
    wu_p: cobj,
    qtu_p: cobj,
    thlu_p: cobj,
    thvu_p: cobj,
    uu_p: cobj,
    vu_p: cobj,
    qtu_emf_p: cobj,
    thlu_emf_p: cobj,
    uu_emf_p: cobj,
    vu_emf_p: cobj,
    uemf_p: cobj,
    dwten_p: cobj,
    diten_p: cobj,
    flxrain_diag_p: cobj,
    flxsnow_diag_p: cobj,
    ntraprd_p: cobj,
    ntsnprd_p: cobj,
    excessu_p: cobj,
    excess0_p: cobj,
    xc_p: cobj,
    aquad_p: cobj,
    bquad_p: cobj,
    cquad_p: cobj,
    bogbot_p: cobj,
    bogtop_p: cobj,
    trflx_p: cobj,
    tru_p: cobj,
    tru_emf_p: cobj,
    fer_out_p: cobj,
    fdr_out_p: cobj,
    cinh_out_p: cobj,
    cinlclh_out_p: cobj,
    qtten_out_p: cobj,
    slten_out_p: cobj,
    ufrc_out_p: cobj,
    uflx_out_p: cobj,
    vflx_out_p: cobj,
    ufrcinvbase_out_p: cobj,
    ufrclcl_out_p: cobj,
    winvbase_out_p: cobj,
    wlcl_out_p: cobj,
    plcl_out_p: cobj,
    pinv_out_p: cobj,
    plfc_out_p: cobj,
    pbup_out_p: cobj,
    ppen_out_p: cobj,
    qtsrc_out_p: cobj,
    thlsrc_out_p: cobj,
    thvlsrc_out_p: cobj,
    emfkbup_out_p: cobj,
    cbmflimit_out_p: cobj,
    tkeavg_out_p: cobj,
    zinv_out_p: cobj,
    rcwp_out_p: cobj,
    rlwp_out_p: cobj,
    riwp_out_p: cobj,
    wu_out_p: cobj,
    qtu_out_p: cobj,
    thlu_out_p: cobj,
    thvu_out_p: cobj,
    uu_out_p: cobj,
    vu_out_p: cobj,
    qtu_emf_out_p: cobj,
    thlu_emf_out_p: cobj,
    uu_emf_out_p: cobj,
    vu_emf_out_p: cobj,
    uemf_out_p: cobj,
    dwten_out_p: cobj,
    diten_out_p: cobj,
    flxrain_out_p: cobj,
    flxsnow_out_p: cobj,
    ntraprd_out_p: cobj,
    ntsnprd_out_p: cobj,
    excessu_out_p: cobj,
    excess0_out_p: cobj,
    xc_out_p: cobj,
    aquad_out_p: cobj,
    bquad_out_p: cobj,
    cquad_out_p: cobj,
    bogbot_out_p: cobj,
    bogtop_out_p: cobj,
    trflx_out_p: cobj,
    tru_out_p: cobj,
    tru_emf_out_p: cobj,
):
    uwshcu_main_wtrc_post_shell_codon(
        mix, mkx, i_col, ncnst, wtrc_nwset,
        precip_v, snow_v, cush_v, cbmf_v, rliq_v, cnt_v, cnb_v,
        umf_p, slflx_p, qtflx_p, flxrain_p, flxsnow_p, qvten_p, qlten_p,
        qiten_p, sten_p, uten_p, vten_p, qrten_p, qsten_p, evapc_p,
        cufrc_p, qcu_p, qlu_p, qiu_p, qc_p, trten_p, wtqc_liq_p,
        wtqc_ice_p, wtprec_p, wtsnow_p, wtrc_iatype_p, umf_out_p,
        slflx_out_p, qtflx_out_p, flxprc1_out_p, flxsnow1_out_p,
        qvten_out_p, qlten_out_p, qiten_out_p, sten_out_p, uten_out_p,
        vten_out_p, qrten_out_p, qsten_out_p, precip_out_p, snow_out_p,
        evapc_out_p, cufrc_out_p, qcu_out_p, qlu_out_p, qiu_out_p,
        cush_out_p, cbmf_out_p, rliq_out_p, qc_out_p, cnt_out_p, cnb_out_p,
        trten_out_p, wtqc_out_p, wtprec_out_p, wtsnow_out_p,
    )
    uwshcu_diag_post_shell_codon(
        mix, mkx, i_col, ncnst,
        cin, cinlcl, ufrcinvbase_v, ufrclcl_v, winvbase_v, wlcl_v,
        plcl_v, pinv_v, plfc_v, pbup_v, ppen_v, qtsrc_v, thlsrc_v,
        thvlsrc_v, emfkbup_v, cbmflimit_v, tkeavg_v, zinv_v, rcwp_v,
        rlwp_v, riwp_v, fer_p, fdr_p, qtten_p, slten_p, ufrc_p, uflx_p,
        vflx_p, wu_p, qtu_p, thlu_p, thvu_p, uu_p, vu_p, qtu_emf_p,
        thlu_emf_p, uu_emf_p, vu_emf_p, uemf_p, dwten_p, diten_p,
        flxrain_diag_p, flxsnow_diag_p, ntraprd_p, ntsnprd_p, excessu_p,
        excess0_p, xc_p, aquad_p, bquad_p, cquad_p, bogbot_p, bogtop_p,
        trflx_p, tru_p, tru_emf_p, fer_out_p, fdr_out_p, cinh_out_p,
        cinlclh_out_p, qtten_out_p, slten_out_p, ufrc_out_p, uflx_out_p,
        vflx_out_p, ufrcinvbase_out_p, ufrclcl_out_p, winvbase_out_p,
        wlcl_out_p, plcl_out_p, pinv_out_p, plfc_out_p, pbup_out_p,
        ppen_out_p, qtsrc_out_p, thlsrc_out_p, thvlsrc_out_p,
        emfkbup_out_p, cbmflimit_out_p, tkeavg_out_p, zinv_out_p,
        rcwp_out_p, rlwp_out_p, riwp_out_p, wu_out_p, qtu_out_p,
        thlu_out_p, thvu_out_p, uu_out_p, vu_out_p, qtu_emf_out_p,
        thlu_emf_out_p, uu_emf_out_p, vu_emf_out_p, uemf_out_p,
        dwten_out_p, diten_out_p, flxrain_out_p, flxsnow_out_p,
        ntraprd_out_p, ntsnprd_out_p, excessu_out_p, excess0_out_p,
        xc_out_p, aquad_out_p, bquad_out_p, cquad_out_p, bogbot_out_p,
        bogtop_out_p, trflx_out_p, tru_out_p, tru_emf_out_p,
    )


@export
def uwshcu_exit_main_zero_shell_codon(
    mix: int,
    mkx: int,
    i_col: int,
    ncnst: int,
    umf_out_p: cobj,
    slflx_out_p: cobj,
    qtflx_out_p: cobj,
    qvten_out_p: cobj,
    qlten_out_p: cobj,
    qiten_out_p: cobj,
    sten_out_p: cobj,
    uten_out_p: cobj,
    vten_out_p: cobj,
    qrten_out_p: cobj,
    qsten_out_p: cobj,
    precip_out_p: cobj,
    snow_out_p: cobj,
    evapc_out_p: cobj,
    cufrc_out_p: cobj,
    qcu_out_p: cobj,
    qlu_out_p: cobj,
    qiu_out_p: cobj,
    cush_out_p: cobj,
    cbmf_out_p: cobj,
    rliq_out_p: cobj,
    qc_out_p: cobj,
    cnt_out_p: cobj,
    cnb_out_p: cobj,
    trten_out_p: cobj,
    wtqc_out_p: cobj,
    wtprec_out_p: cobj,
    wtsnow_out_p: cobj,
):
    umf_out = Ptr[float](umf_out_p)
    slflx_out = Ptr[float](slflx_out_p)
    qtflx_out = Ptr[float](qtflx_out_p)
    qvten_out = Ptr[float](qvten_out_p)
    qlten_out = Ptr[float](qlten_out_p)
    qiten_out = Ptr[float](qiten_out_p)
    sten_out = Ptr[float](sten_out_p)
    uten_out = Ptr[float](uten_out_p)
    vten_out = Ptr[float](vten_out_p)
    qrten_out = Ptr[float](qrten_out_p)
    qsten_out = Ptr[float](qsten_out_p)
    precip_out = Ptr[float](precip_out_p)
    snow_out = Ptr[float](snow_out_p)
    evapc_out = Ptr[float](evapc_out_p)
    cufrc_out = Ptr[float](cufrc_out_p)
    qcu_out = Ptr[float](qcu_out_p)
    qlu_out = Ptr[float](qlu_out_p)
    qiu_out = Ptr[float](qiu_out_p)
    cush_out = Ptr[float](cush_out_p)
    cbmf_out = Ptr[float](cbmf_out_p)
    rliq_out = Ptr[float](rliq_out_p)
    qc_out = Ptr[float](qc_out_p)
    cnt_out = Ptr[float](cnt_out_p)
    cnb_out = Ptr[float](cnb_out_p)
    trten_out = Ptr[float](trten_out_p)
    wtqc_out = Ptr[float](wtqc_out_p)
    wtprec_out = Ptr[float](wtprec_out_p)
    wtsnow_out = Ptr[float](wtsnow_out_p)

    col = i_col - 1
    k = 0
    while k <= mkx:
        idx = col + k * mix
        umf_out[idx] = 0.0
        slflx_out[idx] = 0.0
        qtflx_out[idx] = 0.0
        k += 1

    k = 0
    while k < mkx:
        idx = col + k * mix
        qvten_out[idx] = 0.0
        qlten_out[idx] = 0.0
        qiten_out[idx] = 0.0
        sten_out[idx] = 0.0
        uten_out[idx] = 0.0
        vten_out[idx] = 0.0
        qrten_out[idx] = 0.0
        qsten_out[idx] = 0.0
        evapc_out[idx] = 0.0
        cufrc_out[idx] = 0.0
        qcu_out[idx] = 0.0
        qlu_out[idx] = 0.0
        qiu_out[idx] = 0.0
        qc_out[idx] = 0.0
        k += 1

    precip_out[col] = 0.0
    snow_out[col] = 0.0
    cush_out[col] = -1.0
    cbmf_out[col] = 0.0
    rliq_out[col] = 0.0
    cnt_out[col] = 1.0
    cnb_out[col] = float(mkx)

    m = 0
    while m < ncnst:
        layer_offset = m * mix * mkx
        column_offset = m * mix
        k = 0
        while k < mkx:
            idx = col + k * mix + layer_offset
            trten_out[idx] = 0.0
            wtqc_out[idx] = 0.0
            k += 1
        wtprec_out[col + column_offset] = 0.0
        wtsnow_out[col + column_offset] = 0.0
        m += 1


@export
def uwshcu_exit_diag_zero_shell_codon(
    mix: int,
    mkx: int,
    i_col: int,
    ncnst: int,
    exit_uwcu_p: cobj,
    fer_out_p: cobj,
    fdr_out_p: cobj,
    cinh_out_p: cobj,
    cinlclh_out_p: cobj,
    qtten_out_p: cobj,
    slten_out_p: cobj,
    ufrc_out_p: cobj,
    uflx_out_p: cobj,
    vflx_out_p: cobj,
    ufrcinvbase_out_p: cobj,
    ufrclcl_out_p: cobj,
    winvbase_out_p: cobj,
    wlcl_out_p: cobj,
    plcl_out_p: cobj,
    pinv_out_p: cobj,
    plfc_out_p: cobj,
    pbup_out_p: cobj,
    ppen_out_p: cobj,
    qtsrc_out_p: cobj,
    thlsrc_out_p: cobj,
    thvlsrc_out_p: cobj,
    emfkbup_out_p: cobj,
    cbmflimit_out_p: cobj,
    tkeavg_out_p: cobj,
    zinv_out_p: cobj,
    rcwp_out_p: cobj,
    rlwp_out_p: cobj,
    riwp_out_p: cobj,
    wu_out_p: cobj,
    qtu_out_p: cobj,
    thlu_out_p: cobj,
    thvu_out_p: cobj,
    uu_out_p: cobj,
    vu_out_p: cobj,
    qtu_emf_out_p: cobj,
    thlu_emf_out_p: cobj,
    uu_emf_out_p: cobj,
    vu_emf_out_p: cobj,
    uemf_out_p: cobj,
    dwten_out_p: cobj,
    diten_out_p: cobj,
    flxrain_out_p: cobj,
    flxsnow_out_p: cobj,
    ntraprd_out_p: cobj,
    ntsnprd_out_p: cobj,
    excessu_out_p: cobj,
    excess0_out_p: cobj,
    xc_out_p: cobj,
    aquad_out_p: cobj,
    bquad_out_p: cobj,
    cquad_out_p: cobj,
    bogbot_out_p: cobj,
    bogtop_out_p: cobj,
    trflx_out_p: cobj,
    tru_out_p: cobj,
    tru_emf_out_p: cobj,
):
    exit_uwcu = Ptr[float](exit_uwcu_p)
    fer_out = Ptr[float](fer_out_p)
    fdr_out = Ptr[float](fdr_out_p)
    cinh_out = Ptr[float](cinh_out_p)
    cinlclh_out = Ptr[float](cinlclh_out_p)
    qtten_out = Ptr[float](qtten_out_p)
    slten_out = Ptr[float](slten_out_p)
    ufrc_out = Ptr[float](ufrc_out_p)
    uflx_out = Ptr[float](uflx_out_p)
    vflx_out = Ptr[float](vflx_out_p)
    ufrcinvbase_out = Ptr[float](ufrcinvbase_out_p)
    ufrclcl_out = Ptr[float](ufrclcl_out_p)
    winvbase_out = Ptr[float](winvbase_out_p)
    wlcl_out = Ptr[float](wlcl_out_p)
    plcl_out = Ptr[float](plcl_out_p)
    pinv_out = Ptr[float](pinv_out_p)
    plfc_out = Ptr[float](plfc_out_p)
    pbup_out = Ptr[float](pbup_out_p)
    ppen_out = Ptr[float](ppen_out_p)
    qtsrc_out = Ptr[float](qtsrc_out_p)
    thlsrc_out = Ptr[float](thlsrc_out_p)
    thvlsrc_out = Ptr[float](thvlsrc_out_p)
    emfkbup_out = Ptr[float](emfkbup_out_p)
    cbmflimit_out = Ptr[float](cbmflimit_out_p)
    tkeavg_out = Ptr[float](tkeavg_out_p)
    zinv_out = Ptr[float](zinv_out_p)
    rcwp_out = Ptr[float](rcwp_out_p)
    rlwp_out = Ptr[float](rlwp_out_p)
    riwp_out = Ptr[float](riwp_out_p)
    wu_out = Ptr[float](wu_out_p)
    qtu_out = Ptr[float](qtu_out_p)
    thlu_out = Ptr[float](thlu_out_p)
    thvu_out = Ptr[float](thvu_out_p)
    uu_out = Ptr[float](uu_out_p)
    vu_out = Ptr[float](vu_out_p)
    qtu_emf_out = Ptr[float](qtu_emf_out_p)
    thlu_emf_out = Ptr[float](thlu_emf_out_p)
    uu_emf_out = Ptr[float](uu_emf_out_p)
    vu_emf_out = Ptr[float](vu_emf_out_p)
    uemf_out = Ptr[float](uemf_out_p)
    dwten_out = Ptr[float](dwten_out_p)
    diten_out = Ptr[float](diten_out_p)
    flxrain_out = Ptr[float](flxrain_out_p)
    flxsnow_out = Ptr[float](flxsnow_out_p)
    ntraprd_out = Ptr[float](ntraprd_out_p)
    ntsnprd_out = Ptr[float](ntsnprd_out_p)
    excessu_out = Ptr[float](excessu_out_p)
    excess0_out = Ptr[float](excess0_out_p)
    xc_out = Ptr[float](xc_out_p)
    aquad_out = Ptr[float](aquad_out_p)
    bquad_out = Ptr[float](bquad_out_p)
    cquad_out = Ptr[float](cquad_out_p)
    bogbot_out = Ptr[float](bogbot_out_p)
    bogtop_out = Ptr[float](bogtop_out_p)
    trflx_out = Ptr[float](trflx_out_p)
    tru_out = Ptr[float](tru_out_p)
    tru_emf_out = Ptr[float](tru_emf_out_p)

    col = i_col - 1
    exit_uwcu[col] = 1.0
    cinh_out[col] = -1.0
    cinlclh_out[col] = -1.0
    ufrcinvbase_out[col] = 0.0
    ufrclcl_out[col] = 0.0
    winvbase_out[col] = 0.0
    wlcl_out[col] = 0.0
    plcl_out[col] = 0.0
    pinv_out[col] = 0.0
    plfc_out[col] = 0.0
    pbup_out[col] = 0.0
    ppen_out[col] = 0.0
    qtsrc_out[col] = 0.0
    thlsrc_out[col] = 0.0
    thvlsrc_out[col] = 0.0
    emfkbup_out[col] = 0.0
    cbmflimit_out[col] = 0.0
    tkeavg_out[col] = 0.0
    zinv_out[col] = 0.0
    rcwp_out[col] = 0.0
    rlwp_out[col] = 0.0
    riwp_out[col] = 0.0

    k = 0
    while k < mkx:
        idx = col + k * mix
        fer_out[idx] = 0.0
        fdr_out[idx] = 0.0
        qtten_out[idx] = 0.0
        slten_out[idx] = 0.0
        dwten_out[idx] = 0.0
        diten_out[idx] = 0.0
        ntraprd_out[idx] = 0.0
        ntsnprd_out[idx] = 0.0
        excessu_out[idx] = 0.0
        excess0_out[idx] = 0.0
        xc_out[idx] = 0.0
        aquad_out[idx] = 0.0
        bquad_out[idx] = 0.0
        cquad_out[idx] = 0.0
        bogbot_out[idx] = 0.0
        bogtop_out[idx] = 0.0
        k += 1

    k = 0
    while k <= mkx:
        idx = col + k * mix
        ufrc_out[idx] = 0.0
        uflx_out[idx] = 0.0
        vflx_out[idx] = 0.0
        wu_out[idx] = 0.0
        qtu_out[idx] = 0.0
        thlu_out[idx] = 0.0
        thvu_out[idx] = 0.0
        uu_out[idx] = 0.0
        vu_out[idx] = 0.0
        qtu_emf_out[idx] = 0.0
        thlu_emf_out[idx] = 0.0
        uu_emf_out[idx] = 0.0
        vu_emf_out[idx] = 0.0
        uemf_out[idx] = 0.0
        flxrain_out[idx] = 0.0
        flxsnow_out[idx] = 0.0
        k += 1

    m = 0
    while m < ncnst:
        offset = m * mix * (mkx + 1)
        k = 0
        while k <= mkx:
            idx = col + k * mix + offset
            trflx_out[idx] = 0.0
            tru_out[idx] = 0.0
            tru_emf_out[idx] = 0.0
            k += 1
        m += 1


@export
def uwshcu_exit_zero_all_shell_codon(
    mix: int,
    mkx: int,
    i_col: int,
    ncnst: int,
    umf_out_p: cobj,
    slflx_out_p: cobj,
    qtflx_out_p: cobj,
    qvten_out_p: cobj,
    qlten_out_p: cobj,
    qiten_out_p: cobj,
    sten_out_p: cobj,
    uten_out_p: cobj,
    vten_out_p: cobj,
    qrten_out_p: cobj,
    qsten_out_p: cobj,
    precip_out_p: cobj,
    snow_out_p: cobj,
    evapc_out_p: cobj,
    cufrc_out_p: cobj,
    qcu_out_p: cobj,
    qlu_out_p: cobj,
    qiu_out_p: cobj,
    cush_out_p: cobj,
    cbmf_out_p: cobj,
    rliq_out_p: cobj,
    qc_out_p: cobj,
    cnt_out_p: cobj,
    cnb_out_p: cobj,
    trten_out_p: cobj,
    wtqc_out_p: cobj,
    wtprec_out_p: cobj,
    wtsnow_out_p: cobj,
    exit_uwcu_p: cobj,
    fer_out_p: cobj,
    fdr_out_p: cobj,
    cinh_out_p: cobj,
    cinlclh_out_p: cobj,
    qtten_out_p: cobj,
    slten_out_p: cobj,
    ufrc_out_p: cobj,
    uflx_out_p: cobj,
    vflx_out_p: cobj,
    ufrcinvbase_out_p: cobj,
    ufrclcl_out_p: cobj,
    winvbase_out_p: cobj,
    wlcl_out_p: cobj,
    plcl_out_p: cobj,
    pinv_out_p: cobj,
    plfc_out_p: cobj,
    pbup_out_p: cobj,
    ppen_out_p: cobj,
    qtsrc_out_p: cobj,
    thlsrc_out_p: cobj,
    thvlsrc_out_p: cobj,
    emfkbup_out_p: cobj,
    cbmflimit_out_p: cobj,
    tkeavg_out_p: cobj,
    zinv_out_p: cobj,
    rcwp_out_p: cobj,
    rlwp_out_p: cobj,
    riwp_out_p: cobj,
    wu_out_p: cobj,
    qtu_out_p: cobj,
    thlu_out_p: cobj,
    thvu_out_p: cobj,
    uu_out_p: cobj,
    vu_out_p: cobj,
    qtu_emf_out_p: cobj,
    thlu_emf_out_p: cobj,
    uu_emf_out_p: cobj,
    vu_emf_out_p: cobj,
    uemf_out_p: cobj,
    dwten_out_p: cobj,
    diten_out_p: cobj,
    flxrain_out_p: cobj,
    flxsnow_out_p: cobj,
    ntraprd_out_p: cobj,
    ntsnprd_out_p: cobj,
    excessu_out_p: cobj,
    excess0_out_p: cobj,
    xc_out_p: cobj,
    aquad_out_p: cobj,
    bquad_out_p: cobj,
    cquad_out_p: cobj,
    bogbot_out_p: cobj,
    bogtop_out_p: cobj,
    trflx_out_p: cobj,
    tru_out_p: cobj,
    tru_emf_out_p: cobj,
):
    uwshcu_exit_main_zero_shell_codon(
        mix,
        mkx,
        i_col,
        ncnst,
        umf_out_p,
        slflx_out_p,
        qtflx_out_p,
        qvten_out_p,
        qlten_out_p,
        qiten_out_p,
        sten_out_p,
        uten_out_p,
        vten_out_p,
        qrten_out_p,
        qsten_out_p,
        precip_out_p,
        snow_out_p,
        evapc_out_p,
        cufrc_out_p,
        qcu_out_p,
        qlu_out_p,
        qiu_out_p,
        cush_out_p,
        cbmf_out_p,
        rliq_out_p,
        qc_out_p,
        cnt_out_p,
        cnb_out_p,
        trten_out_p,
        wtqc_out_p,
        wtprec_out_p,
        wtsnow_out_p,
    )
    uwshcu_exit_diag_zero_shell_codon(
        mix,
        mkx,
        i_col,
        ncnst,
        exit_uwcu_p,
        fer_out_p,
        fdr_out_p,
        cinh_out_p,
        cinlclh_out_p,
        qtten_out_p,
        slten_out_p,
        ufrc_out_p,
        uflx_out_p,
        vflx_out_p,
        ufrcinvbase_out_p,
        ufrclcl_out_p,
        winvbase_out_p,
        wlcl_out_p,
        plcl_out_p,
        pinv_out_p,
        plfc_out_p,
        pbup_out_p,
        ppen_out_p,
        qtsrc_out_p,
        thlsrc_out_p,
        thvlsrc_out_p,
        emfkbup_out_p,
        cbmflimit_out_p,
        tkeavg_out_p,
        zinv_out_p,
        rcwp_out_p,
        rlwp_out_p,
        riwp_out_p,
        wu_out_p,
        qtu_out_p,
        thlu_out_p,
        thvu_out_p,
        uu_out_p,
        vu_out_p,
        qtu_emf_out_p,
        thlu_emf_out_p,
        uu_emf_out_p,
        vu_emf_out_p,
        uemf_out_p,
        dwten_out_p,
        diten_out_p,
        flxrain_out_p,
        flxsnow_out_p,
        ntraprd_out_p,
        ntsnprd_out_p,
        excessu_out_p,
        excess0_out_p,
        xc_out_p,
        aquad_out_p,
        bquad_out_p,
        cquad_out_p,
        bogbot_out_p,
        bogtop_out_p,
        trflx_out_p,
        tru_out_p,
        tru_emf_out_p,
    )


@export
def uwshcu_iter_restore_main_shell_codon(
    mix: int,
    mkx: int,
    i_col: int,
    ncnst: int,
    wtrc_nwset: int,
    precip_v: float,
    snow_v: float,
    cush_v: float,
    cbmf_v: float,
    rliq_v: float,
    cnt_v: float,
    cnb_v: float,
    umf_p: cobj,
    qvten_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    sten_p: cobj,
    uten_p: cobj,
    vten_p: cobj,
    qrten_p: cobj,
    qsten_p: cobj,
    evapc_p: cobj,
    cufrc_p: cobj,
    slflx_p: cobj,
    qtflx_p: cobj,
    qcu_p: cobj,
    qlu_p: cobj,
    qiu_p: cobj,
    qc_p: cobj,
    trten_p: cobj,
    wtqc_liq_p: cobj,
    wtqc_ice_p: cobj,
    wtprec_p: cobj,
    wtsnow_p: cobj,
    wtrc_iatype_p: cobj,
    umf_out_p: cobj,
    qvten_out_p: cobj,
    qlten_out_p: cobj,
    qiten_out_p: cobj,
    sten_out_p: cobj,
    uten_out_p: cobj,
    vten_out_p: cobj,
    qrten_out_p: cobj,
    qsten_out_p: cobj,
    precip_out_p: cobj,
    snow_out_p: cobj,
    evapc_out_p: cobj,
    cush_out_p: cobj,
    cufrc_out_p: cobj,
    slflx_out_p: cobj,
    qtflx_out_p: cobj,
    qcu_out_p: cobj,
    qlu_out_p: cobj,
    qiu_out_p: cobj,
    cbmf_out_p: cobj,
    qc_out_p: cobj,
    rliq_out_p: cobj,
    cnt_out_p: cobj,
    cnb_out_p: cobj,
    trten_out_p: cobj,
    wtqc_out_p: cobj,
    wtprec_out_p: cobj,
    wtsnow_out_p: cobj,
):
    umf = Ptr[float](umf_p)
    qvten = Ptr[float](qvten_p)
    qlten = Ptr[float](qlten_p)
    qiten = Ptr[float](qiten_p)
    sten = Ptr[float](sten_p)
    uten = Ptr[float](uten_p)
    vten = Ptr[float](vten_p)
    qrten = Ptr[float](qrten_p)
    qsten = Ptr[float](qsten_p)
    evapc = Ptr[float](evapc_p)
    cufrc = Ptr[float](cufrc_p)
    slflx = Ptr[float](slflx_p)
    qtflx = Ptr[float](qtflx_p)
    qcu = Ptr[float](qcu_p)
    qlu = Ptr[float](qlu_p)
    qiu = Ptr[float](qiu_p)
    qc = Ptr[float](qc_p)
    trten = Ptr[float](trten_p)
    wtqc_liq = Ptr[float](wtqc_liq_p)
    wtqc_ice = Ptr[float](wtqc_ice_p)
    wtprec = Ptr[float](wtprec_p)
    wtsnow = Ptr[float](wtsnow_p)
    wtrc_iatype = Ptr[int](wtrc_iatype_p)
    umf_out = Ptr[float](umf_out_p)
    qvten_out = Ptr[float](qvten_out_p)
    qlten_out = Ptr[float](qlten_out_p)
    qiten_out = Ptr[float](qiten_out_p)
    sten_out = Ptr[float](sten_out_p)
    uten_out = Ptr[float](uten_out_p)
    vten_out = Ptr[float](vten_out_p)
    qrten_out = Ptr[float](qrten_out_p)
    qsten_out = Ptr[float](qsten_out_p)
    precip_out = Ptr[float](precip_out_p)
    snow_out = Ptr[float](snow_out_p)
    evapc_out = Ptr[float](evapc_out_p)
    cush_out = Ptr[float](cush_out_p)
    cufrc_out = Ptr[float](cufrc_out_p)
    slflx_out = Ptr[float](slflx_out_p)
    qtflx_out = Ptr[float](qtflx_out_p)
    qcu_out = Ptr[float](qcu_out_p)
    qlu_out = Ptr[float](qlu_out_p)
    qiu_out = Ptr[float](qiu_out_p)
    cbmf_out = Ptr[float](cbmf_out_p)
    qc_out = Ptr[float](qc_out_p)
    rliq_out = Ptr[float](rliq_out_p)
    cnt_out = Ptr[float](cnt_out_p)
    cnb_out = Ptr[float](cnb_out_p)
    trten_out = Ptr[float](trten_out_p)
    wtqc_out = Ptr[float](wtqc_out_p)
    wtprec_out = Ptr[float](wtprec_out_p)
    wtsnow_out = Ptr[float](wtsnow_out_p)

    col = i_col - 1
    precip_out[col] = precip_v
    snow_out[col] = snow_v
    cush_out[col] = cush_v
    cbmf_out[col] = cbmf_v
    rliq_out[col] = rliq_v
    cnt_out[col] = cnt_v
    cnb_out[col] = cnb_v

    k = 0
    while k <= mkx:
        dst = col + k * mix
        umf_out[dst] = umf[k]
        slflx_out[dst] = slflx[k]
        qtflx_out[dst] = qtflx[k]
        k += 1

    k = 0
    while k < mkx:
        dst = col + k * mix
        qvten_out[dst] = qvten[k]
        qlten_out[dst] = qlten[k]
        qiten_out[dst] = qiten[k]
        sten_out[dst] = sten[k]
        uten_out[dst] = uten[k]
        vten_out[dst] = vten[k]
        qrten_out[dst] = qrten[k]
        qsten_out[dst] = qsten[k]
        evapc_out[dst] = evapc[k]
        cufrc_out[dst] = cufrc[k]
        qcu_out[dst] = qcu[k]
        qlu_out[dst] = qlu[k]
        qiu_out[dst] = qiu[k]
        qc_out[dst] = qc[k]
        k += 1

    m = 0
    while m < ncnst:
        offset = m * mix * mkx
        src_offset = m * mkx
        k = 0
        while k < mkx:
            trten_out[col + k * mix + offset] = trten[k + src_offset]
            k += 1
        m += 1

    m = 0
    while m < wtrc_nwset:
        vap = wtrc_iatype[m] - 1
        liq = wtrc_iatype[m + wtrc_nwset] - 1
        ice = wtrc_iatype[m + 2 * wtrc_nwset] - 1
        src_offset = m * mkx
        liq_offset = liq * mix * mkx
        ice_offset = ice * mix * mkx
        k = 0
        while k < mkx:
            wtqc_out[col + k * mix + liq_offset] = wtqc_liq[k + src_offset]
            wtqc_out[col + k * mix + ice_offset] = wtqc_ice[k + src_offset]
            k += 1
        wtprec_out[col + vap * mix] = wtprec[m]
        wtsnow_out[col + vap * mix] = wtsnow[m]
        m += 1


@export
def uwshcu_iter_restore_diag_shell_codon(
    mix: int,
    mkx: int,
    i_col: int,
    ncnst: int,
    cin_v: float,
    cinlcl_v: float,
    ufrcinvbase_v: float,
    ufrclcl_v: float,
    winvbase_v: float,
    wlcl_v: float,
    plcl_v: float,
    pinv_v: float,
    plfc_v: float,
    pbup_v: float,
    ppen_v: float,
    qtsrc_v: float,
    thlsrc_v: float,
    thvlsrc_v: float,
    emfkbup_v: float,
    cbmflimit_v: float,
    tkeavg_v: float,
    zinv_v: float,
    rcwp_v: float,
    rlwp_v: float,
    riwp_v: float,
    fer_p: cobj,
    fdr_p: cobj,
    qtten_p: cobj,
    slten_p: cobj,
    ufrc_p: cobj,
    uflx_p: cobj,
    vflx_p: cobj,
    wu_p: cobj,
    qtu_p: cobj,
    thlu_p: cobj,
    thvu_p: cobj,
    uu_p: cobj,
    vu_p: cobj,
    qtu_emf_p: cobj,
    thlu_emf_p: cobj,
    uu_emf_p: cobj,
    vu_emf_p: cobj,
    uemf_p: cobj,
    dwten_p: cobj,
    diten_p: cobj,
    flxrain_p: cobj,
    flxsnow_p: cobj,
    ntraprd_p: cobj,
    ntsnprd_p: cobj,
    excessu_p: cobj,
    excess0_p: cobj,
    xc_p: cobj,
    aquad_p: cobj,
    bquad_p: cobj,
    cquad_p: cobj,
    bogbot_p: cobj,
    bogtop_p: cobj,
    trflx_p: cobj,
    tru_p: cobj,
    tru_emf_p: cobj,
    fer_out_p: cobj,
    fdr_out_p: cobj,
    cinh_out_p: cobj,
    cinlclh_out_p: cobj,
    qtten_out_p: cobj,
    slten_out_p: cobj,
    ufrc_out_p: cobj,
    uflx_out_p: cobj,
    vflx_out_p: cobj,
    ufrcinvbase_out_p: cobj,
    ufrclcl_out_p: cobj,
    winvbase_out_p: cobj,
    wlcl_out_p: cobj,
    plcl_out_p: cobj,
    pinv_out_p: cobj,
    plfc_out_p: cobj,
    pbup_out_p: cobj,
    ppen_out_p: cobj,
    qtsrc_out_p: cobj,
    thlsrc_out_p: cobj,
    thvlsrc_out_p: cobj,
    emfkbup_out_p: cobj,
    cbmflimit_out_p: cobj,
    tkeavg_out_p: cobj,
    zinv_out_p: cobj,
    rcwp_out_p: cobj,
    rlwp_out_p: cobj,
    riwp_out_p: cobj,
    wu_out_p: cobj,
    qtu_out_p: cobj,
    thlu_out_p: cobj,
    thvu_out_p: cobj,
    uu_out_p: cobj,
    vu_out_p: cobj,
    qtu_emf_out_p: cobj,
    thlu_emf_out_p: cobj,
    uu_emf_out_p: cobj,
    vu_emf_out_p: cobj,
    uemf_out_p: cobj,
    dwten_out_p: cobj,
    diten_out_p: cobj,
    flxrain_out_p: cobj,
    flxsnow_out_p: cobj,
    ntraprd_out_p: cobj,
    ntsnprd_out_p: cobj,
    excessu_out_p: cobj,
    excess0_out_p: cobj,
    xc_out_p: cobj,
    aquad_out_p: cobj,
    bquad_out_p: cobj,
    cquad_out_p: cobj,
    bogbot_out_p: cobj,
    bogtop_out_p: cobj,
    trflx_out_p: cobj,
    tru_out_p: cobj,
    tru_emf_out_p: cobj,
):
    fer = Ptr[float](fer_p)
    fdr = Ptr[float](fdr_p)
    qtten = Ptr[float](qtten_p)
    slten = Ptr[float](slten_p)
    ufrc = Ptr[float](ufrc_p)
    uflx = Ptr[float](uflx_p)
    vflx = Ptr[float](vflx_p)
    wu = Ptr[float](wu_p)
    qtu = Ptr[float](qtu_p)
    thlu = Ptr[float](thlu_p)
    thvu = Ptr[float](thvu_p)
    uu = Ptr[float](uu_p)
    vu = Ptr[float](vu_p)
    qtu_emf = Ptr[float](qtu_emf_p)
    thlu_emf = Ptr[float](thlu_emf_p)
    uu_emf = Ptr[float](uu_emf_p)
    vu_emf = Ptr[float](vu_emf_p)
    uemf = Ptr[float](uemf_p)
    dwten = Ptr[float](dwten_p)
    diten = Ptr[float](diten_p)
    flxrain = Ptr[float](flxrain_p)
    flxsnow = Ptr[float](flxsnow_p)
    ntraprd = Ptr[float](ntraprd_p)
    ntsnprd = Ptr[float](ntsnprd_p)
    excessu = Ptr[float](excessu_p)
    excess0 = Ptr[float](excess0_p)
    xc = Ptr[float](xc_p)
    aquad = Ptr[float](aquad_p)
    bquad = Ptr[float](bquad_p)
    cquad = Ptr[float](cquad_p)
    bogbot = Ptr[float](bogbot_p)
    bogtop = Ptr[float](bogtop_p)
    trflx = Ptr[float](trflx_p)
    tru = Ptr[float](tru_p)
    tru_emf = Ptr[float](tru_emf_p)
    fer_out = Ptr[float](fer_out_p)
    fdr_out = Ptr[float](fdr_out_p)
    cinh_out = Ptr[float](cinh_out_p)
    cinlclh_out = Ptr[float](cinlclh_out_p)
    qtten_out = Ptr[float](qtten_out_p)
    slten_out = Ptr[float](slten_out_p)
    ufrc_out = Ptr[float](ufrc_out_p)
    uflx_out = Ptr[float](uflx_out_p)
    vflx_out = Ptr[float](vflx_out_p)
    ufrcinvbase_out = Ptr[float](ufrcinvbase_out_p)
    ufrclcl_out = Ptr[float](ufrclcl_out_p)
    winvbase_out = Ptr[float](winvbase_out_p)
    wlcl_out = Ptr[float](wlcl_out_p)
    plcl_out = Ptr[float](plcl_out_p)
    pinv_out = Ptr[float](pinv_out_p)
    plfc_out = Ptr[float](plfc_out_p)
    pbup_out = Ptr[float](pbup_out_p)
    ppen_out = Ptr[float](ppen_out_p)
    qtsrc_out = Ptr[float](qtsrc_out_p)
    thlsrc_out = Ptr[float](thlsrc_out_p)
    thvlsrc_out = Ptr[float](thvlsrc_out_p)
    emfkbup_out = Ptr[float](emfkbup_out_p)
    cbmflimit_out = Ptr[float](cbmflimit_out_p)
    tkeavg_out = Ptr[float](tkeavg_out_p)
    zinv_out = Ptr[float](zinv_out_p)
    rcwp_out = Ptr[float](rcwp_out_p)
    rlwp_out = Ptr[float](rlwp_out_p)
    riwp_out = Ptr[float](riwp_out_p)
    wu_out = Ptr[float](wu_out_p)
    qtu_out = Ptr[float](qtu_out_p)
    thlu_out = Ptr[float](thlu_out_p)
    thvu_out = Ptr[float](thvu_out_p)
    uu_out = Ptr[float](uu_out_p)
    vu_out = Ptr[float](vu_out_p)
    qtu_emf_out = Ptr[float](qtu_emf_out_p)
    thlu_emf_out = Ptr[float](thlu_emf_out_p)
    uu_emf_out = Ptr[float](uu_emf_out_p)
    vu_emf_out = Ptr[float](vu_emf_out_p)
    uemf_out = Ptr[float](uemf_out_p)
    dwten_out = Ptr[float](dwten_out_p)
    diten_out = Ptr[float](diten_out_p)
    flxrain_out = Ptr[float](flxrain_out_p)
    flxsnow_out = Ptr[float](flxsnow_out_p)
    ntraprd_out = Ptr[float](ntraprd_out_p)
    ntsnprd_out = Ptr[float](ntsnprd_out_p)
    excessu_out = Ptr[float](excessu_out_p)
    excess0_out = Ptr[float](excess0_out_p)
    xc_out = Ptr[float](xc_out_p)
    aquad_out = Ptr[float](aquad_out_p)
    bquad_out = Ptr[float](bquad_out_p)
    cquad_out = Ptr[float](cquad_out_p)
    bogbot_out = Ptr[float](bogbot_out_p)
    bogtop_out = Ptr[float](bogtop_out_p)
    trflx_out = Ptr[float](trflx_out_p)
    tru_out = Ptr[float](tru_out_p)
    tru_emf_out = Ptr[float](tru_emf_out_p)

    col = i_col - 1
    cinh_out[col] = cin_v
    cinlclh_out[col] = cinlcl_v
    ufrcinvbase_out[col] = ufrcinvbase_v
    ufrclcl_out[col] = ufrclcl_v
    winvbase_out[col] = winvbase_v
    wlcl_out[col] = wlcl_v
    plcl_out[col] = plcl_v
    pinv_out[col] = pinv_v
    plfc_out[col] = plfc_v
    pbup_out[col] = pbup_v
    ppen_out[col] = ppen_v
    qtsrc_out[col] = qtsrc_v
    thlsrc_out[col] = thlsrc_v
    thvlsrc_out[col] = thvlsrc_v
    emfkbup_out[col] = emfkbup_v
    cbmflimit_out[col] = cbmflimit_v
    tkeavg_out[col] = tkeavg_v
    zinv_out[col] = zinv_v
    rcwp_out[col] = rcwp_v
    rlwp_out[col] = rlwp_v
    riwp_out[col] = riwp_v

    k = 0
    while k < mkx:
        dst_k = mkx - k - 1
        dst = col + dst_k * mix
        fer_out[dst] = fer[k]
        fdr_out[dst] = fdr[k]
        qtten_out[dst] = qtten[k]
        slten_out[dst] = slten[k]
        dwten_out[dst] = dwten[k]
        diten_out[dst] = diten[k]
        ntraprd_out[dst] = ntraprd[k]
        ntsnprd_out[dst] = ntsnprd[k]
        excessu_out[dst] = excessu[k]
        excess0_out[dst] = excess0[k]
        xc_out[dst] = xc[k]
        aquad_out[dst] = aquad[k]
        bquad_out[dst] = bquad[k]
        cquad_out[dst] = cquad[k]
        bogbot_out[dst] = bogbot[k]
        bogtop_out[dst] = bogtop[k]
        k += 1

    k = 0
    while k <= mkx:
        dst_k = mkx - k
        dst = col + dst_k * mix
        ufrc_out[dst] = ufrc[k]
        uflx_out[dst] = uflx[k]
        vflx_out[dst] = vflx[k]
        wu_out[dst] = wu[k]
        qtu_out[dst] = qtu[k]
        thlu_out[dst] = thlu[k]
        thvu_out[dst] = thvu[k]
        uu_out[dst] = uu[k]
        vu_out[dst] = vu[k]
        qtu_emf_out[dst] = qtu_emf[k]
        thlu_emf_out[dst] = thlu_emf[k]
        uu_emf_out[dst] = uu_emf[k]
        vu_emf_out[dst] = vu_emf[k]
        uemf_out[dst] = uemf[k]
        flxrain_out[dst] = flxrain[k]
        flxsnow_out[dst] = flxsnow[k]
        k += 1

    m = 0
    while m < ncnst:
        offset = m * mix * (mkx + 1)
        src_offset = m * (mkx + 1)
        k = 0
        while k <= mkx:
            dst_k = mkx - k
            dst = col + dst_k * mix + offset
            src = k + src_offset
            trflx_out[dst] = trflx[src]
            tru_out[dst] = tru[src]
            tru_emf_out[dst] = tru_emf[src]
            k += 1
        m += 1


@export
def uwshcu_iter_restore_all_shell_codon(
    mix: int,
    mkx: int,
    i_col: int,
    ncnst: int,
    wtrc_nwset: int,
    precip_v: float,
    snow_v: float,
    cush_v: float,
    cbmf_v: float,
    rliq_v: float,
    cnt_v: float,
    cnb_v: float,
    umf_p: cobj,
    qvten_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    sten_p: cobj,
    uten_p: cobj,
    vten_p: cobj,
    qrten_p: cobj,
    qsten_p: cobj,
    evapc_p: cobj,
    cufrc_p: cobj,
    slflx_p: cobj,
    qtflx_p: cobj,
    qcu_p: cobj,
    qlu_p: cobj,
    qiu_p: cobj,
    qc_p: cobj,
    trten_p: cobj,
    wtqc_liq_p: cobj,
    wtqc_ice_p: cobj,
    wtprec_p: cobj,
    wtsnow_p: cobj,
    wtrc_iatype_p: cobj,
    umf_out_p: cobj,
    qvten_out_p: cobj,
    qlten_out_p: cobj,
    qiten_out_p: cobj,
    sten_out_p: cobj,
    uten_out_p: cobj,
    vten_out_p: cobj,
    qrten_out_p: cobj,
    qsten_out_p: cobj,
    precip_out_p: cobj,
    snow_out_p: cobj,
    evapc_out_p: cobj,
    cush_out_p: cobj,
    cufrc_out_p: cobj,
    slflx_out_p: cobj,
    qtflx_out_p: cobj,
    qcu_out_p: cobj,
    qlu_out_p: cobj,
    qiu_out_p: cobj,
    cbmf_out_p: cobj,
    qc_out_p: cobj,
    rliq_out_p: cobj,
    cnt_out_p: cobj,
    cnb_out_p: cobj,
    trten_out_p: cobj,
    wtqc_out_p: cobj,
    wtprec_out_p: cobj,
    wtsnow_out_p: cobj,
    cin_v: float,
    cinlcl_v: float,
    ufrcinvbase_v: float,
    ufrclcl_v: float,
    winvbase_v: float,
    wlcl_v: float,
    plcl_v: float,
    pinv_v: float,
    plfc_v: float,
    pbup_v: float,
    ppen_v: float,
    qtsrc_v: float,
    thlsrc_v: float,
    thvlsrc_v: float,
    emfkbup_v: float,
    cbmflimit_v: float,
    tkeavg_v: float,
    zinv_v: float,
    rcwp_v: float,
    rlwp_v: float,
    riwp_v: float,
    fer_p: cobj,
    fdr_p: cobj,
    qtten_p: cobj,
    slten_p: cobj,
    ufrc_p: cobj,
    uflx_p: cobj,
    vflx_p: cobj,
    wu_p: cobj,
    qtu_p: cobj,
    thlu_p: cobj,
    thvu_p: cobj,
    uu_p: cobj,
    vu_p: cobj,
    qtu_emf_p: cobj,
    thlu_emf_p: cobj,
    uu_emf_p: cobj,
    vu_emf_p: cobj,
    uemf_p: cobj,
    dwten_p: cobj,
    diten_p: cobj,
    flxrain_p: cobj,
    flxsnow_p: cobj,
    ntraprd_p: cobj,
    ntsnprd_p: cobj,
    excessu_p: cobj,
    excess0_p: cobj,
    xc_p: cobj,
    aquad_p: cobj,
    bquad_p: cobj,
    cquad_p: cobj,
    bogbot_p: cobj,
    bogtop_p: cobj,
    trflx_p: cobj,
    tru_p: cobj,
    tru_emf_p: cobj,
    fer_out_p: cobj,
    fdr_out_p: cobj,
    cinh_out_p: cobj,
    cinlclh_out_p: cobj,
    qtten_out_p: cobj,
    slten_out_p: cobj,
    ufrc_out_p: cobj,
    uflx_out_p: cobj,
    vflx_out_p: cobj,
    ufrcinvbase_out_p: cobj,
    ufrclcl_out_p: cobj,
    winvbase_out_p: cobj,
    wlcl_out_p: cobj,
    plcl_out_p: cobj,
    pinv_out_p: cobj,
    plfc_out_p: cobj,
    pbup_out_p: cobj,
    ppen_out_p: cobj,
    qtsrc_out_p: cobj,
    thlsrc_out_p: cobj,
    thvlsrc_out_p: cobj,
    emfkbup_out_p: cobj,
    cbmflimit_out_p: cobj,
    tkeavg_out_p: cobj,
    zinv_out_p: cobj,
    rcwp_out_p: cobj,
    rlwp_out_p: cobj,
    riwp_out_p: cobj,
    wu_out_p: cobj,
    qtu_out_p: cobj,
    thlu_out_p: cobj,
    thvu_out_p: cobj,
    uu_out_p: cobj,
    vu_out_p: cobj,
    qtu_emf_out_p: cobj,
    thlu_emf_out_p: cobj,
    uu_emf_out_p: cobj,
    vu_emf_out_p: cobj,
    uemf_out_p: cobj,
    dwten_out_p: cobj,
    diten_out_p: cobj,
    flxrain_out_p: cobj,
    flxsnow_out_p: cobj,
    ntraprd_out_p: cobj,
    ntsnprd_out_p: cobj,
    excessu_out_p: cobj,
    excess0_out_p: cobj,
    xc_out_p: cobj,
    aquad_out_p: cobj,
    bquad_out_p: cobj,
    cquad_out_p: cobj,
    bogbot_out_p: cobj,
    bogtop_out_p: cobj,
    trflx_out_p: cobj,
    tru_out_p: cobj,
    tru_emf_out_p: cobj,
):
    uwshcu_iter_restore_main_shell_codon(
        mix, mkx, i_col, ncnst, wtrc_nwset, precip_v, snow_v, cush_v, cbmf_v,
        rliq_v, cnt_v, cnb_v, umf_p, qvten_p, qlten_p, qiten_p, sten_p, uten_p,
        vten_p, qrten_p, qsten_p, evapc_p, cufrc_p, slflx_p, qtflx_p, qcu_p,
        qlu_p, qiu_p, qc_p, trten_p, wtqc_liq_p, wtqc_ice_p, wtprec_p,
        wtsnow_p, wtrc_iatype_p, umf_out_p, qvten_out_p, qlten_out_p,
        qiten_out_p, sten_out_p, uten_out_p, vten_out_p, qrten_out_p,
        qsten_out_p, precip_out_p, snow_out_p, evapc_out_p, cush_out_p,
        cufrc_out_p, slflx_out_p, qtflx_out_p, qcu_out_p, qlu_out_p,
        qiu_out_p, cbmf_out_p, qc_out_p, rliq_out_p, cnt_out_p, cnb_out_p,
        trten_out_p, wtqc_out_p, wtprec_out_p, wtsnow_out_p
    )
    uwshcu_iter_restore_diag_shell_codon(
        mix, mkx, i_col, ncnst, cin_v, cinlcl_v, ufrcinvbase_v, ufrclcl_v,
        winvbase_v, wlcl_v, plcl_v, pinv_v, plfc_v, pbup_v, ppen_v, qtsrc_v,
        thlsrc_v, thvlsrc_v, emfkbup_v, cbmflimit_v, tkeavg_v, zinv_v, rcwp_v,
        rlwp_v, riwp_v, fer_p, fdr_p, qtten_p, slten_p, ufrc_p, uflx_p, vflx_p,
        wu_p, qtu_p, thlu_p, thvu_p, uu_p, vu_p, qtu_emf_p, thlu_emf_p,
        uu_emf_p, vu_emf_p, uemf_p, dwten_p, diten_p, flxrain_p, flxsnow_p,
        ntraprd_p, ntsnprd_p, excessu_p, excess0_p, xc_p, aquad_p, bquad_p,
        cquad_p, bogbot_p, bogtop_p, trflx_p, tru_p, tru_emf_p, fer_out_p,
        fdr_out_p, cinh_out_p, cinlclh_out_p, qtten_out_p, slten_out_p,
        ufrc_out_p, uflx_out_p, vflx_out_p, ufrcinvbase_out_p, ufrclcl_out_p,
        winvbase_out_p, wlcl_out_p, plcl_out_p, pinv_out_p, plfc_out_p,
        pbup_out_p, ppen_out_p, qtsrc_out_p, thlsrc_out_p, thvlsrc_out_p,
        emfkbup_out_p, cbmflimit_out_p, tkeavg_out_p, zinv_out_p, rcwp_out_p,
        rlwp_out_p, riwp_out_p, wu_out_p, qtu_out_p, thlu_out_p, thvu_out_p,
        uu_out_p, vu_out_p, qtu_emf_out_p, thlu_emf_out_p, uu_emf_out_p,
        vu_emf_out_p, uemf_out_p, dwten_out_p, diten_out_p, flxrain_out_p,
        flxsnow_out_p, ntraprd_out_p, ntsnprd_out_p, excessu_out_p,
        excess0_out_p, xc_out_p, aquad_out_p, bquad_out_p, cquad_out_p,
        bogbot_out_p, bogtop_out_p, trflx_out_p, tru_out_p, tru_emf_out_p
    )


@export
def uwshcu_column_thermo_state_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    xlv: float,
    xls: float,
    cp: float,
    zvir: float,
    qv0_p: cobj,
    ql0_p: cobj,
    qi0_p: cobj,
    t0_p: cobj,
    exn0_p: cobj,
    tr0_p: cobj,
    wtrc_iatype_p: cobj,
    qt0_p: cobj,
    thl0_p: cobj,
    thvl0_p: cobj,
    wt0_p: cobj,
):
    qv0 = Ptr[float](qv0_p)
    ql0 = Ptr[float](ql0_p)
    qi0 = Ptr[float](qi0_p)
    t0 = Ptr[float](t0_p)
    exn0 = Ptr[float](exn0_p)
    tr0 = Ptr[float](tr0_p)
    wtrc_iatype = Ptr[int](wtrc_iatype_p)
    qt0 = Ptr[float](qt0_p)
    thl0 = Ptr[float](thl0_p)
    thvl0 = Ptr[float](thvl0_p)
    wt0 = Ptr[float](wt0_p)

    k = 0
    while k < mkx:
        qt0[k] = qv0[k] + ql0[k] + qi0[k]
        thl0[k] = (t0[k] - xlv * ql0[k] / cp - xls * qi0[k] / cp) / exn0[k]
        thvl0[k] = (1.0 + zvir * qt0[k]) * thl0[k]
        k += 1

    m = 0
    while m < wtrc_nwset:
        vap = wtrc_iatype[m] - 1
        liq = wtrc_iatype[m + wtrc_nwset] - 1
        ice = wtrc_iatype[m + 2 * wtrc_nwset] - 1
        dst_offset = m * mkx
        vap_offset = vap * mkx
        liq_offset = liq * mkx
        ice_offset = ice * mkx
        k = 0
        while k < mkx:
            wt0[k + dst_offset] = tr0[k + vap_offset] + tr0[k + liq_offset] + tr0[k + ice_offset]
            k += 1
        m += 1


@export
def uwshcu_pbl_precheck_shell_codon(
    mkx: int,
    pblh: float,
    zs0_p: cobj,
    cush_p: cobj,
    tscaleh_p: cobj,
    kinv_out_p: cobj,
    exit_code_p: cobj,
):
    zs0 = Ptr[float](zs0_p)
    cush = Ptr[float](cush_p)
    tscaleh = Ptr[float](tscaleh_p)
    kinv_out = Ptr[int](kinv_out_p)
    exit_code = Ptr[int](exit_code_p)

    tscaleh[0] = cush[0]
    cush[0] = -1.0

    kinv = 1
    k = mkx - 1
    while k >= 1:
        if (pblh + 5.0 - zs0[k]) * (pblh + 5.0 - zs0[k + 1]) < 0.0:
            kinv = k + 1
            break
        k -= 1

    kinv_out[0] = kinv
    exit_code[0] = 0
    if kinv <= 1:
        exit_code[0] = 1


@export
def uwshcu_pbl_source_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    kinv: int,
    zvir: float,
    ps0_p: cobj,
    p0_p: cobj,
    tke_p: cobj,
    thvl0bot_p: cobj,
    thvl0top_p: cobj,
    qt0_p: cobj,
    u0_p: cobj,
    v0_p: cobj,
    ssu0_p: cobj,
    ssv0_p: cobj,
    tr0_p: cobj,
    wt0_p: cobj,
    tkeavg_p: cobj,
    thvlmin_p: cobj,
    qtsrc_p: cobj,
    thvlsrc_p: cobj,
    thlsrc_p: cobj,
    usrc_p: cobj,
    vsrc_p: cobj,
    trsrc_p: cobj,
    wtsrc_p: cobj,
):
    ps0 = Ptr[float](ps0_p)
    p0 = Ptr[float](p0_p)
    tke = Ptr[float](tke_p)
    thvl0bot = Ptr[float](thvl0bot_p)
    thvl0top = Ptr[float](thvl0top_p)
    qt0 = Ptr[float](qt0_p)
    u0 = Ptr[float](u0_p)
    v0 = Ptr[float](v0_p)
    ssu0 = Ptr[float](ssu0_p)
    ssv0 = Ptr[float](ssv0_p)
    tr0 = Ptr[float](tr0_p)
    wt0 = Ptr[float](wt0_p)
    tkeavg = Ptr[float](tkeavg_p)
    thvlmin = Ptr[float](thvlmin_p)
    qtsrc = Ptr[float](qtsrc_p)
    thvlsrc = Ptr[float](thvlsrc_p)
    thlsrc = Ptr[float](thlsrc_p)
    usrc = Ptr[float](usrc_p)
    vsrc = Ptr[float](vsrc_p)
    trsrc = Ptr[float](trsrc_p)
    wtsrc = Ptr[float](wtsrc_p)

    dpsum = 0.0
    tkeavg_val = 0.0
    thvlmin_val = 1000.0
    k = 0
    while k < kinv:
        if k == 0:
            dpi = ps0[0] - p0[0]
        elif k == kinv - 1:
            dpi = p0[kinv - 2] - ps0[kinv - 1]
        else:
            dpi = p0[k - 1] - p0[k]
        dpsum = dpsum + dpi
        tkeavg_val = tkeavg_val + dpi * tke[k]
        if k != 0:
            thvl_layer_min = thvl0bot[k - 1]
            if thvl0top[k - 1] < thvl_layer_min:
                thvl_layer_min = thvl0top[k - 1]
            if thvl_layer_min < thvlmin_val:
                thvlmin_val = thvl_layer_min
        k += 1

    tkeavg[0] = tkeavg_val / dpsum
    thvlmin[0] = thvlmin_val
    qtsrc[0] = qt0[0]
    thvlsrc[0] = thvlmin[0]
    thlsrc[0] = thvlsrc[0] / (1.0 + zvir * qtsrc[0])
    usrc[0] = u0[kinv - 2] + ssu0[kinv - 2] * (ps0[kinv - 1] - p0[kinv - 2])
    vsrc[0] = v0[kinv - 2] + ssv0[kinv - 2] * (ps0[kinv - 1] - p0[kinv - 2])

    m = 0
    while m < ncnst:
        trsrc[m] = tr0[m * mkx]
        m += 1

    m = 0
    while m < wtrc_nwset:
        wtsrc[m] = wt0[m * mkx]
        m += 1


@export
def uwshcu_pbl_precheck_source_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    pblh: float,
    zvir: float,
    zs0_p: cobj,
    cush_p: cobj,
    tscaleh_p: cobj,
    kinv_out_p: cobj,
    exit_code_p: cobj,
    ps0_p: cobj,
    p0_p: cobj,
    tke_p: cobj,
    thvl0bot_p: cobj,
    thvl0top_p: cobj,
    qt0_p: cobj,
    u0_p: cobj,
    v0_p: cobj,
    ssu0_p: cobj,
    ssv0_p: cobj,
    tr0_p: cobj,
    wt0_p: cobj,
    tkeavg_p: cobj,
    thvlmin_p: cobj,
    qtsrc_p: cobj,
    thvlsrc_p: cobj,
    thlsrc_p: cobj,
    usrc_p: cobj,
    vsrc_p: cobj,
    trsrc_p: cobj,
    wtsrc_p: cobj,
):
    uwshcu_pbl_precheck_shell_codon(
        mkx,
        pblh,
        zs0_p,
        cush_p,
        tscaleh_p,
        kinv_out_p,
        exit_code_p,
    )

    kinv_out = Ptr[int](kinv_out_p)
    exit_code = Ptr[int](exit_code_p)
    if exit_code[0] != 0:
        return

    uwshcu_pbl_source_shell_codon(
        mkx,
        ncnst,
        wtrc_nwset,
        kinv_out[0],
        zvir,
        ps0_p,
        p0_p,
        tke_p,
        thvl0bot_p,
        thvl0top_p,
        qt0_p,
        u0_p,
        v0_p,
        ssu0_p,
        ssv0_p,
        tr0_p,
        wt0_p,
        tkeavg_p,
        thvlmin_p,
        qtsrc_p,
        thvlsrc_p,
        thlsrc_p,
        usrc_p,
        vsrc_p,
        trsrc_p,
        wtsrc_p,
    )


@export
def uwshcu_lcl_prep_shell_codon(
    mkx: int,
    plcl: float,
    ps0_p: cobj,
    p0_p: cobj,
    thl0_p: cobj,
    ssthl0_p: cobj,
    qt0_p: cobj,
    ssqt0_p: cobj,
    klcl_out_p: cobj,
    exit_code_p: cobj,
    thl0lcl_p: cobj,
    qt0lcl_p: cobj,
):
    ps0 = Ptr[float](ps0_p)
    p0 = Ptr[float](p0_p)
    thl0 = Ptr[float](thl0_p)
    ssthl0 = Ptr[float](ssthl0_p)
    qt0 = Ptr[float](qt0_p)
    ssqt0 = Ptr[float](ssqt0_p)
    klcl_out = Ptr[int](klcl_out_p)
    exit_code = Ptr[int](exit_code_p)
    thl0lcl = Ptr[float](thl0lcl_p)
    qt0lcl = Ptr[float](qt0lcl_p)

    klcl = mkx
    k = 0
    while k <= mkx:
        if ps0[k] < plcl:
            klcl = k
            break
        k += 1
    if klcl < 1:
        klcl = 1

    klcl_out[0] = klcl
    exit_code[0] = 0
    if plcl < 30000.0:
        exit_code[0] = 1
        return

    idx = klcl - 1
    thl0lcl[0] = thl0[idx] + ssthl0[idx] * (plcl - p0[idx])
    qt0lcl[0] = qt0[idx] + ssqt0[idx] * (plcl - p0[idx])


@export
def uwshcu_interface_thv_shell_codon(
    k_fortran: int,
    zvir: float,
    thj: float,
    qvj: float,
    qlj: float,
    qij: float,
    thl0edge: float,
    qt0edge: float,
    thv0_p: cobj,
    thvl0_p: cobj,
):
    thv0 = Ptr[float](thv0_p)
    thvl0 = Ptr[float](thvl0_p)

    idx = k_fortran - 1
    thv0[idx] = thj * (1.0 + zvir * qvj - qlj - qij)
    thvl0[idx] = thl0edge * (1.0 + zvir * qt0edge)


@export
def uwshcu_interface_conden_exit_shell_codon(
    id_check: int,
    exit_conden_p: cobj,
    exit_code_p: cobj,
):
    exit_conden = Ptr[float](exit_conden_p)
    exit_code = Ptr[int](exit_code_p)

    exit_code[0] = 0
    if id_check == 1:
        exit_conden[0] = 1.0
        exit_code[0] = 1


@export
def uwshcu_cin_conden_exit_shell_codon(
    id_check: int,
    exit_conden_p: cobj,
    exit_code_p: cobj,
):
    exit_conden = Ptr[float](exit_conden_p)
    exit_code = Ptr[int](exit_code_p)

    exit_code[0] = 0
    if id_check == 1:
        exit_conden[0] = 1.0
        exit_code[0] = 1


@export
def uwshcu_thv_scalar_shell_codon(
    zvir: float,
    thj: float,
    qvj: float,
    qlj: float,
    qij: float,
    thv_p: cobj,
):
    thv = Ptr[float](thv_p)

    thv[0] = thj * (1.0 + zvir * qvj - qlj - qij)


@export
def uwshcu_conden_exit_thv_batch_shell_codon(
    kind: int,
    k_fortran: int,
    id_check: int,
    zvir: float,
    thj: float,
    qvj: float,
    qlj: float,
    qij: float,
    thl0edge: float,
    qt0edge: float,
    exit_conden_p: cobj,
    exit_code_p: cobj,
    thv_p: cobj,
    thvl_p: cobj,
):
    exit_conden = Ptr[float](exit_conden_p)
    exit_code = Ptr[int](exit_code_p)

    exit_code[0] = 0
    if id_check == 1:
        exit_conden[0] = 1.0
        exit_code[0] = 1
        return

    if kind == 1:
        thv = Ptr[float](thv_p)
        thv[0] = thj * (1.0 + zvir * qvj - qlj - qij)
    elif kind == 2:
        thv = Ptr[float](thv_p)
        thvl = Ptr[float](thvl_p)
        idx = k_fortran - 1
        thv[idx] = thj * (1.0 + zvir * qvj - qlj - qij)
        thvl[idx] = thl0edge * (1.0 + zvir * qt0edge)


@export
def uwshcu_buoy_env_pre_qsat_shell_codon(
    zvir: float,
    r_v: float,
    pe_v: float,
    thj: float,
    qvj: float,
    qlj: float,
    qij: float,
    exne: float,
    thle: float,
    thv0j_p: cobj,
    rho0j_p: cobj,
    qsat_arg_p: cobj,
):
    thv0j = Ptr[float](thv0j_p)
    rho0j = Ptr[float](rho0j_p)
    qsat_arg = Ptr[float](qsat_arg_p)

    thv0j[0] = thj * (1.0 + zvir * qvj - qlj - qij)
    rho0j[0] = pe_v / (r_v * thv0j[0] * exne)
    qsat_arg[0] = thle * exne


@export
def uwshcu_buoy_excess_shell_codon(
    qt_v: float,
    qs_v: float,
    excess_p: cobj,
):
    excess = Ptr[float](excess_p)

    excess[0] = qt_v - qs_v


@export
def uwshcu_buoy_detrain_excess_shell_codon(
    criqc: float,
    xlv: float,
    xls: float,
    cp: float,
    exne: float,
    qlj: float,
    qij: float,
    thlue_p: cobj,
    qtue_p: cobj,
):
    thlue = Ptr[float](thlue_p)
    qtue = Ptr[float](qtue_p)

    if (qlj + qij) > criqc:
        exql = ((qlj + qij) - criqc) * qlj / (qlj + qij)
        exqi = ((qlj + qij) - criqc) * qij / (qlj + qij)
        qtue[0] = qtue[0] - exql - exqi
        thlue[0] = thlue[0] + (xlv / cp / exne) * exql + (xls / cp / exne) * exqi


@export
def uwshcu_buoy_up_pre_qsat_shell_codon(
    zvir: float,
    exne: float,
    thj: float,
    qvj: float,
    qlj: float,
    qij: float,
    thlue: float,
    thvj_p: cobj,
    tj_p: cobj,
    qsat_arg_p: cobj,
):
    thvj = Ptr[float](thvj_p)
    tj = Ptr[float](tj_p)
    qsat_arg = Ptr[float](qsat_arg_p)

    thvj[0] = thj * (1.0 + zvir * qvj - qlj - qij)
    tj[0] = thj * exne
    qsat_arg[0] = thlue * exne


@export
def uwshcu_buoy_conden_scalar_batch_shell_codon(
    kind: int,
    id_check: int,
    zvir: float,
    r_v: float,
    pe_v: float,
    thj: float,
    qvj: float,
    qlj: float,
    qij: float,
    exne: float,
    thle: float,
    criqc: float,
    xlv: float,
    xls: float,
    cp_v: float,
    v1: float,
    v2: float,
    exit_conden_p: cobj,
    exit_code_p: cobj,
    thlue_p: cobj,
    qtue_p: cobj,
    thv0j_p: cobj,
    rho0j_p: cobj,
    thvj_p: cobj,
    tj_p: cobj,
    thvxsat_p: cobj,
    qsat_arg_p: cobj,
    excess_p: cobj,
):
    if kind == 5:
        excess = Ptr[float](excess_p)
        excess[0] = v1 - v2
        return

    exit_conden = Ptr[float](exit_conden_p)
    exit_code = Ptr[int](exit_code_p)
    exit_code[0] = 0
    if id_check == 1:
        exit_conden[0] = 1.0
        exit_code[0] = 1
        return

    if kind == 1:
        thv0j = Ptr[float](thv0j_p)
        rho0j = Ptr[float](rho0j_p)
        qsat_arg = Ptr[float](qsat_arg_p)
        thv0j[0] = thj * (1.0 + zvir * qvj - qlj - qij)
        rho0j[0] = pe_v / (r_v * thv0j[0] * exne)
        qsat_arg[0] = thle * exne
    elif kind == 2:
        thlue = Ptr[float](thlue_p)
        qtue = Ptr[float](qtue_p)
        if (qlj + qij) > criqc:
            exql = ((qlj + qij) - criqc) * qlj / (qlj + qij)
            exqi = ((qlj + qij) - criqc) * qij / (qlj + qij)
            qtue[0] = qtue[0] - exql - exqi
            thlue[0] = thlue[0] + (xlv / cp_v / exne) * exql + (xls / cp_v / exne) * exqi
    elif kind == 3:
        thvj = Ptr[float](thvj_p)
        tj = Ptr[float](tj_p)
        qsat_arg = Ptr[float](qsat_arg_p)
        thvj[0] = thj * (1.0 + zvir * qvj - qlj - qij)
        tj[0] = thj * exne
        qsat_arg[0] = v1 * exne
    elif kind == 4:
        thvxsat = Ptr[float](thvxsat_p)
        thvxsat[0] = thj * (1.0 + zvir * qvj - qlj - qij)


@export
def uwshcu_buoy_top_expel_shell_codon(
    k_fortran: int,
    criqc: float,
    xlv: float,
    xls: float,
    cp: float,
    exn: float,
    qlj: float,
    qij: float,
    qtu_p: cobj,
    thlu_p: cobj,
    dwten_p: cobj,
    diten_p: cobj,
):
    # qtu/thlu: Fortran real(r8) qtu(0:mkx), thlu(0:mkx).
    # dwten/diten: Fortran real(r8) dwten(mkx), diten(mkx).
    qtu = Ptr[float](qtu_p)
    thlu = Ptr[float](thlu_p)
    dwten = Ptr[float](dwten_p)
    diten = Ptr[float](diten_p)

    iface_idx = k_fortran
    layer_idx = k_fortran - 1
    condensate = qlj + qij
    if condensate > criqc:
        exql = (condensate - criqc) * qlj / condensate
        exqi = (condensate - criqc) * qij / condensate
        qtu[iface_idx] = qtu[iface_idx] - exql - exqi
        thlu[iface_idx] = thlu[iface_idx] + (xlv / cp / exn) * exql + (xls / cp / exn) * exqi
        dwten[layer_idx] = exql
        diten[layer_idx] = exqi
    else:
        dwten[layer_idx] = 0.0
        diten[layer_idx] = 0.0


@export
def uwshcu_buoy_top_state_shell_codon(
    mkx: int,
    wtrc_nwset: int,
    kpen: int,
    trace_water: int,
    linear_branch: int,
    ppen: float,
    top_expfac: float,
    fer_kpen: float,
    thl0_kpen: float,
    ssthl0_kpen: float,
    qt0_kpen: float,
    ssqt0_kpen: float,
    thlu_p: cobj,
    qtu_p: cobj,
    wt0_p: cobj,
    sswt0_p: cobj,
    wtu_p: cobj,
    thlu_top_p: cobj,
    qtu_top_p: cobj,
    wtu_top_p: cobj,
):
    # thlu/qtu/wtu use Fortran interface lower bound 0; wt0/sswt0/wtu_top lower bound 1.
    thlu = Ptr[float](thlu_p)
    qtu = Ptr[float](qtu_p)
    wt0 = Ptr[float](wt0_p)
    sswt0 = Ptr[float](sswt0_p)
    wtu = Ptr[float](wtu_p)
    thlu_top = Ptr[float](thlu_top_p)
    qtu_top = Ptr[float](qtu_top_p)
    wtu_top = Ptr[float](wtu_top_p)

    km1 = kpen - 1
    neg_ppen = -ppen
    iface_stride = mkx + 1

    if linear_branch != 0:
        thlu_top[0] = thlu[km1] + (thl0_kpen + ssthl0_kpen * neg_ppen / 2.0 - thlu[km1]) * fer_kpen * neg_ppen
        qtu_top[0] = qtu[km1] + (qt0_kpen + ssqt0_kpen * neg_ppen / 2.0 - qtu[km1]) * fer_kpen * neg_ppen

        if trace_water != 0:
            m = 0
            while m < wtrc_nwset:
                layer_idx = km1 + m * mkx
                iface_idx = km1 + m * iface_stride
                wtu_top[m] = wtu[iface_idx] + (wt0[layer_idx] + sswt0[layer_idx] * neg_ppen / 2.0 - wtu[iface_idx]) * fer_kpen * neg_ppen
                m += 1
    else:
        thlu_top[0] = (
            thl0_kpen + ssthl0_kpen / fer_kpen - ssthl0_kpen * neg_ppen / 2.0
        ) - (
            thl0_kpen + ssthl0_kpen * neg_ppen / 2.0 - thlu[km1] + ssthl0_kpen / fer_kpen
        ) * top_expfac
        qtu_top[0] = (
            qt0_kpen + ssqt0_kpen / fer_kpen - ssqt0_kpen * neg_ppen / 2.0
        ) - (
            qt0_kpen + ssqt0_kpen * neg_ppen / 2.0 - qtu[km1] + ssqt0_kpen / fer_kpen
        ) * top_expfac

        if trace_water != 0:
            m = 0
            while m < wtrc_nwset:
                layer_idx = km1 + m * mkx
                iface_idx = km1 + m * iface_stride
                wtu_top[m] = (
                    wt0[layer_idx] + sswt0[layer_idx] / fer_kpen - sswt0[layer_idx] * neg_ppen / 2.0
                ) - (
                    wt0[layer_idx] + sswt0[layer_idx] * neg_ppen / 2.0 - wtu[iface_idx] + sswt0[layer_idx] / fer_kpen
                ) * top_expfac
                m += 1


@export
def uwshcu_buoy_updraft_state_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    k_fortran: int,
    trace_water: int,
    linear_branch: int,
    dpe: float,
    expfac: float,
    fer_k: float,
    PGFc: float,
    thle: float,
    qte: float,
    ue: float,
    ve: float,
    thlu_p: cobj,
    qtu_p: cobj,
    uu_p: cobj,
    vu_p: cobj,
    tru_p: cobj,
    wtu_p: cobj,
    ssthl0_p: cobj,
    ssqt0_p: cobj,
    ssu0_p: cobj,
    ssv0_p: cobj,
    tre_p: cobj,
    sstr0_p: cobj,
    wte_p: cobj,
    sswt0_p: cobj,
):
    # thlu/qtu/uu/vu/wtu/tru have first-dimension lower bound 0; slope arrays use layer lower bound 1.
    thlu = Ptr[float](thlu_p)
    qtu = Ptr[float](qtu_p)
    uu = Ptr[float](uu_p)
    vu = Ptr[float](vu_p)
    tru = Ptr[float](tru_p)
    wtu = Ptr[float](wtu_p)
    ssthl0 = Ptr[float](ssthl0_p)
    ssqt0 = Ptr[float](ssqt0_p)
    ssu0 = Ptr[float](ssu0_p)
    ssv0 = Ptr[float](ssv0_p)
    tre = Ptr[float](tre_p)
    sstr0 = Ptr[float](sstr0_p)
    wte = Ptr[float](wte_p)
    sswt0 = Ptr[float](sswt0_p)

    layer_idx = k_fortran - 1
    iface_idx = k_fortran
    km1_idx = k_fortran - 1
    iface_stride = mkx + 1

    if linear_branch != 0:
        thlu[iface_idx] = thlu[km1_idx] + (thle + ssthl0[layer_idx] * dpe / 2.0 - thlu[km1_idx]) * fer_k * dpe
        qtu[iface_idx] = qtu[km1_idx] + (qte + ssqt0[layer_idx] * dpe / 2.0 - qtu[km1_idx]) * fer_k * dpe
        uu[iface_idx] = (
            uu[km1_idx] + (ue + ssu0[layer_idx] * dpe / 2.0 - uu[km1_idx]) * fer_k * dpe
            - PGFc * ssu0[layer_idx] * dpe
        )
        vu[iface_idx] = (
            vu[km1_idx] + (ve + ssv0[layer_idx] * dpe / 2.0 - vu[km1_idx]) * fer_k * dpe
            - PGFc * ssv0[layer_idx] * dpe
        )

        m = 0
        while m < ncnst:
            layer_m_idx = layer_idx + m * mkx
            iface_m_idx = iface_idx + m * iface_stride
            km1_m_idx = km1_idx + m * iface_stride
            tru[iface_m_idx] = (
                tru[km1_m_idx] + (tre[m] + sstr0[layer_m_idx] * dpe / 2.0 - tru[km1_m_idx]) * fer_k * dpe
            )
            m += 1

        if trace_water != 0:
            m = 0
            while m < wtrc_nwset:
                layer_m_idx = layer_idx + m * mkx
                iface_m_idx = iface_idx + m * iface_stride
                km1_m_idx = km1_idx + m * iface_stride
                wtu[iface_m_idx] = (
                    wtu[km1_m_idx] + (wte[m] + sswt0[layer_m_idx] * dpe / 2.0 - wtu[km1_m_idx]) * fer_k * dpe
                )
                m += 1
    else:
        thlu[iface_idx] = (
            thle + ssthl0[layer_idx] / fer_k - ssthl0[layer_idx] * dpe / 2.0
        ) - (
            thle + ssthl0[layer_idx] * dpe / 2.0 - thlu[km1_idx] + ssthl0[layer_idx] / fer_k
        ) * expfac
        qtu[iface_idx] = (
            qte + ssqt0[layer_idx] / fer_k - ssqt0[layer_idx] * dpe / 2.0
        ) - (
            qte + ssqt0[layer_idx] * dpe / 2.0 - qtu[km1_idx] + ssqt0[layer_idx] / fer_k
        ) * expfac
        uu[iface_idx] = (
            ue + (1.0 - PGFc) * ssu0[layer_idx] / fer_k - ssu0[layer_idx] * dpe / 2.0
        ) - (
            ue + ssu0[layer_idx] * dpe / 2.0 - uu[km1_idx] + (1.0 - PGFc) * ssu0[layer_idx] / fer_k
        ) * expfac
        vu[iface_idx] = (
            ve + (1.0 - PGFc) * ssv0[layer_idx] / fer_k - ssv0[layer_idx] * dpe / 2.0
        ) - (
            ve + ssv0[layer_idx] * dpe / 2.0 - vu[km1_idx] + (1.0 - PGFc) * ssv0[layer_idx] / fer_k
        ) * expfac

        m = 0
        while m < ncnst:
            layer_m_idx = layer_idx + m * mkx
            iface_m_idx = iface_idx + m * iface_stride
            km1_m_idx = km1_idx + m * iface_stride
            tru[iface_m_idx] = (
                tre[m] + sstr0[layer_m_idx] / fer_k - sstr0[layer_m_idx] * dpe / 2.0
            ) - (
                tre[m] + sstr0[layer_m_idx] * dpe / 2.0 - tru[km1_m_idx] + sstr0[layer_m_idx] / fer_k
            ) * expfac
            m += 1

        if trace_water != 0:
            m = 0
            while m < wtrc_nwset:
                layer_m_idx = layer_idx + m * mkx
                iface_m_idx = iface_idx + m * iface_stride
                km1_m_idx = km1_idx + m * iface_stride
                wtu[iface_m_idx] = (
                    wte[m] + sswt0[layer_m_idx] / fer_k - sswt0[layer_m_idx] * dpe / 2.0
                ) - (
                    wte[m] + sswt0[layer_m_idx] * dpe / 2.0 - wtu[km1_m_idx] + sswt0[layer_m_idx] / fer_k
                ) * expfac
                m += 1


@export
def uwshcu_buoy_velocity_shell_codon(
    rbuoy: float,
    thvu_km1: float,
    thvebot: float,
    thvu_k: float,
    thv0top_k: float,
    drage: float,
    dpe: float,
    expfac: float,
    rho0j: float,
    bogbot_p: cobj,
    bogtop_p: cobj,
    delbog_p: cobj,
    wtwb_p: cobj,
    wtw_p: cobj,
):
    bogbot = Ptr[float](bogbot_p)
    bogtop = Ptr[float](bogtop_p)
    delbog = Ptr[float](delbog_p)
    wtwb = Ptr[float](wtwb_p)
    wtw = Ptr[float](wtw_p)

    bogbot[0] = rbuoy * (thvu_km1 / thvebot - 1.0)
    bogtop[0] = rbuoy * (thvu_k / thv0top_k - 1.0)
    delbog[0] = bogtop[0] - bogbot[0]
    wtwb[0] = wtw[0]

    if drage * dpe > 1.0e-3:
        wtw[0] = wtw[0] * expfac + (
            delbog[0] + (1.0 - expfac) * (bogbot[0] + delbog[0] / (-2.0 * drage * dpe))
        ) / (rho0j * drage)
    else:
        wtw[0] = wtw[0] + dpe * (bogbot[0] + bogtop[0]) / rho0j


@export
def uwshcu_buoy_midstate_shell_codon(
    k_fortran: int,
    thlu_p: cobj,
    qtu_p: cobj,
    thlue_p: cobj,
    qtue_p: cobj,
):
    # thlu/qtu: Fortran real(r8) thlu(0:mkx), qtu(0:mkx).
    thlu = Ptr[float](thlu_p)
    qtu = Ptr[float](qtu_p)
    thlue = Ptr[float](thlue_p)
    qtue = Ptr[float](qtue_p)

    thlue[0] = 0.5 * (thlu[k_fortran - 1] + thlu[k_fortran])
    qtue[0] = 0.5 * (qtu[k_fortran - 1] + qtu[k_fortran])


@export
def uwshcu_buoy_self_detrain_shell_codon(
    k_fortran: int,
    use_self_detrain: int,
    expfac: float,
    umf_p: cobj,
    wtw_p: cobj,
):
    # umf: Fortran real(r8) umf(0:mkx); wtw is a scalar.
    umf = Ptr[float](umf_p)
    wtw = Ptr[float](wtw_p)

    if use_self_detrain != 0:
        umf[k_fortran] = umf[k_fortran] * expfac

    if umf[k_fortran] == 0.0:
        wtw[0] = -1.0


@export
def uwshcu_buoy_ufrc_init_shell_codon(
    k_fortran: int,
    r_v: float,
    ps0_p: cobj,
    thv0bot_p: cobj,
    thv0top_p: cobj,
    exns0_p: cobj,
    umf_p: cobj,
    wu_p: cobj,
    ufrc_p: cobj,
    rhos0j_p: cobj,
):
    # ps0/exns0/umf/wu/ufrc: Fortran real(r8) arrays with lower bound 0.
    # thv0bot/thv0top: Fortran real(r8) arrays with lower bound 1.
    ps0 = Ptr[float](ps0_p)
    thv0bot = Ptr[float](thv0bot_p)
    thv0top = Ptr[float](thv0top_p)
    exns0 = Ptr[float](exns0_p)
    umf = Ptr[float](umf_p)
    wu = Ptr[float](wu_p)
    ufrc = Ptr[float](ufrc_p)
    rhos0j = Ptr[float](rhos0j_p)

    rhos0j[0] = ps0[k_fortran] / (r_v * 0.5 * (thv0bot[k_fortran] + thv0top[k_fortran - 1]) * exns0[k_fortran])
    ufrc[k_fortran] = umf[k_fortran] / (rhos0j[0] * wu[k_fortran])


@export
def uwshcu_buoy_ufrc_limit_shell_codon(
    k_fortran: int,
    rmaxfrac_v: float,
    rhos0j_v: float,
    ufrc_p: cobj,
    umf_p: cobj,
    wu_p: cobj,
    limit_ufrc_p: cobj,
    limit_code_p: cobj,
):
    # ufrc/umf/wu: Fortran real(r8) arrays with lower bound 0.
    ufrc = Ptr[float](ufrc_p)
    umf = Ptr[float](umf_p)
    wu = Ptr[float](wu_p)
    limit_ufrc = Ptr[float](limit_ufrc_p)
    limit_code = Ptr[int](limit_code_p)

    limit_code[0] = 0
    if ufrc[k_fortran] > rmaxfrac_v:
        limit_ufrc[0] = 1.0
        ufrc[k_fortran] = rmaxfrac_v
        umf[k_fortran] = rmaxfrac_v * rhos0j_v * wu[k_fortran]
        limit_code[0] = 1


@export
def uwshcu_buoy_ppen_limit_shell_codon(
    ppen: float,
    dp0_kpen: float,
    limit_ppen_p: cobj,
):
    limit_ppen = Ptr[float](limit_ppen_p)

    if ppen == -dp0_kpen or ppen == 0.0:
        limit_ppen[0] = 1.0


@export
def uwshcu_buoy_conden_exit_shell_codon(
    id_check: int,
    exit_conden_p: cobj,
    exit_code_p: cobj,
):
    exit_conden = Ptr[float](exit_conden_p)
    exit_code = Ptr[int](exit_code_p)

    exit_code[0] = 0
    if id_check == 1:
        exit_conden[0] = 1.0
        exit_code[0] = 1


@export
def uwshcu_buoy_top_conden_exit_shell_codon(
    id_check: int,
    exit_conden_p: cobj,
    exit_code_p: cobj,
):
    exit_conden = Ptr[float](exit_conden_p)
    exit_code = Ptr[int](exit_code_p)

    exit_code[0] = 0
    if id_check == 1:
        exit_conden[0] = 1.0
        exit_code[0] = 1


@export
def uwshcu_buoy_top_expel_final_shell_codon(
    kpen: int,
    criqc: float,
    xlv: float,
    xls: float,
    cp: float,
    exntop: float,
    qlj: float,
    qij: float,
    thlu_top_p: cobj,
    qtu_top_p: cobj,
    dwten_p: cobj,
    diten_p: cobj,
):
    # dwten/diten: Fortran real(r8) dwten(mkx), diten(mkx).
    thlu_top = Ptr[float](thlu_top_p)
    qtu_top = Ptr[float](qtu_top_p)
    dwten = Ptr[float](dwten_p)
    diten = Ptr[float](diten_p)

    layer_idx = kpen - 1
    condensate = qlj + qij
    if condensate > criqc:
        dwten[layer_idx] = (condensate - criqc) * qlj / condensate
        diten[layer_idx] = (condensate - criqc) * qij / condensate
        qtu_top[0] = qtu_top[0] - dwten[layer_idx] - diten[layer_idx]
        thlu_top[0] = thlu_top[0] + (xlv / cp / exntop) * dwten[layer_idx] + (xls / cp / exntop) * diten[layer_idx]
    else:
        dwten[layer_idx] = 0.0
        diten[layer_idx] = 0.0


@export
def uwshcu_buoy_scaleh_shell_codon(
    kpen: int,
    r_v: float,
    g_v: float,
    ppen: float,
    ps0_p: cobj,
    zs0_p: cobj,
    thv0bot_p: cobj,
    thv0top_p: cobj,
    exns0_p: cobj,
    cush_p: cobj,
    scaleh_p: cobj,
):
    # ps0/zs0/exns0: Fortran real(r8) ps0(0:mkx), zs0(0:mkx), exns0(0:mkx).
    # thv0bot/thv0top: Fortran real(r8) thv0bot(mkx), thv0top(mkx).
    ps0 = Ptr[float](ps0_p)
    zs0 = Ptr[float](zs0_p)
    thv0bot = Ptr[float](thv0bot_p)
    thv0top = Ptr[float](thv0top_p)
    exns0 = Ptr[float](exns0_p)
    cush = Ptr[float](cush_p)
    scaleh = Ptr[float](scaleh_p)

    iface_idx = kpen - 1
    rho = ps0[iface_idx] / (r_v * 0.5 * (thv0bot[kpen - 1] + thv0top[kpen - 2]) * exns0[iface_idx])
    cush[0] = zs0[iface_idx] - ppen / rho / g_v
    scaleh[0] = cush[0]


@export
def uwshcu_buoy_diag_update_shell_codon(
    k_fortran: int,
    excessu_v: float,
    excess0_v: float,
    xc_v: float,
    aquad_v: float,
    bquad_v: float,
    cquad_v: float,
    bogbot_v: float,
    bogtop_v: float,
    excessu_arr_p: cobj,
    excess0_arr_p: cobj,
    xc_arr_p: cobj,
    aquad_arr_p: cobj,
    bquad_arr_p: cobj,
    cquad_arr_p: cobj,
    bogbot_arr_p: cobj,
    bogtop_arr_p: cobj,
):
    # Diagnostic arrays: Fortran real(r8) *_arr(mkx).
    excessu_arr = Ptr[float](excessu_arr_p)
    excess0_arr = Ptr[float](excess0_arr_p)
    xc_arr = Ptr[float](xc_arr_p)
    aquad_arr = Ptr[float](aquad_arr_p)
    bquad_arr = Ptr[float](bquad_arr_p)
    cquad_arr = Ptr[float](cquad_arr_p)
    bogbot_arr = Ptr[float](bogbot_arr_p)
    bogtop_arr = Ptr[float](bogtop_arr_p)

    idx = k_fortran - 1
    excessu_arr[idx] = excessu_v
    excess0_arr[idx] = excess0_v
    xc_arr[idx] = xc_v
    aquad_arr[idx] = aquad_v
    bquad_arr[idx] = bquad_v
    cquad_arr[idx] = cquad_v
    bogbot_arr[idx] = bogbot_v
    bogtop_arr[idx] = bogtop_v


@export
def uwshcu_buoy_reach_update_shell_codon(
    k_fortran: int,
    bogtop: float,
    wtw: float,
    kbup_v: int,
    kpen_v: int,
    kbup_p: cobj,
    kpen_p: cobj,
    exit_code_p: cobj,
):
    kbup = Ptr[int](kbup_p)
    kpen = Ptr[int](kpen_p)
    exit_code = Ptr[int](exit_code_p)

    kbup[0] = kbup_v
    kpen[0] = kpen_v
    exit_code[0] = 0

    if bogtop > 0.0 and wtw > 0.0:
        kbup[0] = k_fortran

    if wtw <= 0.0:
        kpen[0] = k_fortran
        exit_code[0] = 1


@export
def uwshcu_buoy_state_batch_shell_codon(
    kind: int,
    k_fortran: int,
    mkx: int,
    wtrc_nwset: int,
    i1: int,
    i2: int,
    flag1: int,
    flag2: int,
    v1: float,
    v2: float,
    v3: float,
    v4: float,
    v5: float,
    v6: float,
    v7: float,
    v8: float,
    v9: float,
    v10: float,
    v11: float,
    v12: float,
    p1: cobj,
    p2: cobj,
    p3: cobj,
    p4: cobj,
    p5: cobj,
    p6: cobj,
    p7: cobj,
    p8: cobj,
):
    # Batch dispatcher for already validated buoy helpers. Pointer slots are per-kind:
    # 1 top_expel: qtu, thlu, dwten, diten.
    # 2 velocity: bogbot, bogtop, delbog, wtwb, wtw.
    # 3 midstate: thlu, qtu, thlue, qtue.
    # 4 self_detrain: umf, wtw.
    # 5 diag_update: excessu/excess0/xc/aquad/bquad/cquad/bogbot/bogtop arrays.
    # 6 reach_update: kbup, kpen, exit_code.
    # 7 ufrc_init: ps0, thv0bot, thv0top, exns0, umf, wu, ufrc, rhos0j.
    # 8 top_state: thlu, qtu, wt0, sswt0, wtu, thlu_top, qtu_top, wtu_top.
    # 9 top_expel_final: thlu_top, qtu_top, dwten, diten.
    # 10 scaleh: ps0, zs0, thv0bot, thv0top, exns0, cush, scaleh.
    if kind == 1:
        uwshcu_buoy_top_expel_shell_codon(k_fortran, v1, v2, v3, v4, v5, v6, v7, p1, p2, p3, p4)
    elif kind == 2:
        uwshcu_buoy_velocity_shell_codon(v1, v2, v3, v4, v5, v6, v7, v8, v9, p1, p2, p3, p4, p5)
    elif kind == 3:
        uwshcu_buoy_midstate_shell_codon(k_fortran, p1, p2, p3, p4)
    elif kind == 4:
        uwshcu_buoy_self_detrain_shell_codon(k_fortran, flag1, v1, p1, p2)
    elif kind == 5:
        uwshcu_buoy_diag_update_shell_codon(k_fortran, v1, v2, v3, v4, v5, v6, v7, v8, p1, p2, p3, p4, p5, p6, p7, p8)
    elif kind == 6:
        uwshcu_buoy_reach_update_shell_codon(k_fortran, v1, v2, i1, i2, p1, p2, p3)
    elif kind == 7:
        uwshcu_buoy_ufrc_init_shell_codon(k_fortran, v1, p1, p2, p3, p4, p5, p6, p7, p8)
    elif kind == 8:
        uwshcu_buoy_top_state_shell_codon(
            mkx, wtrc_nwset, k_fortran, flag1, flag2, v1, v2, v3, v4, v5, v6, v7, p1, p2, p3, p4, p5, p6, p7, p8
        )
    elif kind == 9:
        uwshcu_buoy_top_expel_final_shell_codon(k_fortran, v1, v2, v3, v4, v5, v6, v7, p1, p2, p3, p4)
    elif kind == 10:
        uwshcu_buoy_scaleh_shell_codon(k_fortran, v1, v2, v3, p1, p2, p3, p4, p5, p6, p7)


@export
def uwshcu_scalar_exit_limit_batch_shell_codon(
    kind: int,
    idx: int,
    code_in: int,
    v1: float,
    v2: float,
    v3: float,
    field1_p: cobj,
    field2_p: cobj,
    field3_p: cobj,
    flag1_p: cobj,
    flag2_p: cobj,
    code_out_p: cobj,
):
    if kind == 1:
        flag1 = Ptr[float](flag1_p)
        code_out = Ptr[int](code_out_p)
        code_out[0] = 0
        if code_in == 1:
            flag1[0] = 1.0
            code_out[0] = 1
    elif kind == 2:
        code_out = Ptr[int](code_out_p)
        code_out[0] = 0
        if code_in == 1:
            code_out[0] = 1
    elif kind == 3:
        code_out = Ptr[int](code_out_p)
        code_out[0] = 0
        if v1 >= 3.0:
            code_out[0] = 1
    elif kind == 4:
        flag1 = Ptr[float](flag1_p)
        if v1 == v2:
            flag1[0] = 1.0
    elif kind == 5:
        flag1 = Ptr[float](flag1_p)
        flag2 = Ptr[float](flag2_p)
        if v1 == v2:
            flag1[0] = 1.0
        if v1 == v3:
            flag2[0] = 1.0
    elif kind == 6:
        flag1 = Ptr[float](flag1_p)
        code_out = Ptr[int](code_out_p)
        code_out[0] = 0
        if v1 <= 0.0:
            flag1[0] = 1.0
            code_out[0] = 1
    elif kind == 7:
        flag1 = Ptr[float](flag1_p)
        code_out = Ptr[int](code_out_p)
        code_out[0] = 0
        if v1 <= 0.0001:
            flag1[0] = 1.0
            code_out[0] = 1
    elif kind == 8:
        flag1 = Ptr[float](flag1_p)
        code_out = Ptr[int](code_out_p)
        code_out[0] = 0
        if v1 > 100.0:
            flag1[0] = 1.0
            code_out[0] = 1
    elif kind == 9:
        ufrc = Ptr[float](field1_p)
        umf = Ptr[float](field2_p)
        wu = Ptr[float](field3_p)
        flag1 = Ptr[float](flag1_p)
        code_out = Ptr[int](code_out_p)
        code_out[0] = 0
        if ufrc[idx] > v1:
            flag1[0] = 1.0
            ufrc[idx] = v1
            umf[idx] = v1 * v2 * wu[idx]
            code_out[0] = 1
    elif kind == 10:
        flag1 = Ptr[float](flag1_p)
        if v1 == -v2 or v1 == 0.0:
            flag1[0] = 1.0
    elif kind == 11:
        flag1 = Ptr[float](flag1_p)
        code_out = Ptr[int](code_out_p)
        code_out[0] = 0
        if code_in != 0:
            flag1[0] = 1.0
            code_out[0] = 1


@export
def uwshcu_buoy_wu_exit_shell_codon(
    wu_v: float,
    exit_wu_p: cobj,
    exit_code_p: cobj,
):
    exit_wu = Ptr[float](exit_wu_p)
    exit_code = Ptr[int](exit_code_p)

    exit_code[0] = 0
    if wu_v > 100.0:
        exit_wu[0] = 1.0
        exit_code[0] = 1


@export
def uwshcu_release_mu_exit_shell_codon(
    mu_v: float,
    exit_code_p: cobj,
):
    exit_code = Ptr[int](exit_code_p)

    exit_code[0] = 0
    if mu_v >= 3.0:
        exit_code[0] = 1


@export
def uwshcu_release_mumin2_limit_shell_codon(
    mu_v: float,
    mumin2_v: float,
    limit_ufrc_p: cobj,
):
    limit_ufrc = Ptr[float](limit_ufrc_p)

    if mu_v == mumin2_v:
        limit_ufrc[0] = 1.0


@export
def uwshcu_release_mu_limit_flags_shell_codon(
    mu_v: float,
    mumin0_v: float,
    mumin1_v: float,
    limit_cbmf_p: cobj,
    limit_ufrc_p: cobj,
):
    limit_cbmf = Ptr[float](limit_cbmf_p)
    limit_ufrc = Ptr[float](limit_ufrc_p)

    if mu_v == mumin0_v:
        limit_cbmf[0] = 1.0
    if mu_v == mumin1_v:
        limit_ufrc[0] = 1.0


@export
def uwshcu_release_wtw_exit_shell_codon(
    wtw_v: float,
    exit_wtw_p: cobj,
    exit_code_p: cobj,
):
    exit_wtw = Ptr[float](exit_wtw_p)
    exit_code = Ptr[int](exit_code_p)

    exit_code[0] = 0
    if wtw_v <= 0.0:
        exit_wtw[0] = 1.0
        exit_code[0] = 1


@export
def uwshcu_release_ufrc_exit_shell_codon(
    ufrclcl_v: float,
    exit_ufrc_p: cobj,
    exit_code_p: cobj,
):
    exit_ufrc = Ptr[float](exit_ufrc_p)
    exit_code = Ptr[int](exit_code_p)

    exit_code[0] = 0
    if ufrclcl_v <= 0.0001:
        exit_ufrc[0] = 1.0
        exit_code[0] = 1


@export
def uwshcu_buoy_next_env_load_shell_codon(
    k_fortran: int,
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    p0_p: cobj,
    dp0_p: cobj,
    exn0_p: cobj,
    thv0bot_p: cobj,
    thl0_p: cobj,
    qt0_p: cobj,
    u0_p: cobj,
    v0_p: cobj,
    tr0_p: cobj,
    wt0_p: cobj,
    pe_p: cobj,
    dpe_p: cobj,
    exne_p: cobj,
    thvebot_p: cobj,
    thle_p: cobj,
    qte_p: cobj,
    ue_p: cobj,
    ve_p: cobj,
    tre_p: cobj,
    wte_p: cobj,
):
    # p0/dp0/exn0/thv0bot/thl0/qt0/u0/v0: Fortran real(r8) arrays with extent mkx.
    # tr0/wt0: Fortran real(r8) tr0(mkx,ncnst), wt0(mkx,wtrc_nwset).
    p0 = Ptr[float](p0_p)
    dp0 = Ptr[float](dp0_p)
    exn0 = Ptr[float](exn0_p)
    thv0bot = Ptr[float](thv0bot_p)
    thl0 = Ptr[float](thl0_p)
    qt0 = Ptr[float](qt0_p)
    u0 = Ptr[float](u0_p)
    v0 = Ptr[float](v0_p)
    tr0 = Ptr[float](tr0_p)
    wt0 = Ptr[float](wt0_p)
    pe = Ptr[float](pe_p)
    dpe = Ptr[float](dpe_p)
    exne = Ptr[float](exne_p)
    thvebot = Ptr[float](thvebot_p)
    thle = Ptr[float](thle_p)
    qte = Ptr[float](qte_p)
    ue = Ptr[float](ue_p)
    ve = Ptr[float](ve_p)
    tre = Ptr[float](tre_p)
    wte = Ptr[float](wte_p)

    layer_idx = k_fortran
    pe[0] = p0[layer_idx]
    dpe[0] = dp0[layer_idx]
    exne[0] = exn0[layer_idx]
    thvebot[0] = thv0bot[layer_idx]
    thle[0] = thl0[layer_idx]
    qte[0] = qt0[layer_idx]
    ue[0] = u0[layer_idx]
    ve[0] = v0[layer_idx]

    m = 0
    while m < ncnst:
        tre[m] = tr0[layer_idx + m * mkx]
        m += 1

    m = 0
    while m < wtrc_nwset:
        wte[m] = wt0[layer_idx + m * mkx]
        m += 1


@export
def uwshcu_buoy_loop_batch_shell_codon(
    kind: int,
    k_fortran: int,
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    flag1: int,
    flag2: int,
    id_check: int,
    v1: float,
    v2: float,
    v3: float,
    v4: float,
    v5: float,
    v6: float,
    v7: float,
    v8: float,
    v9: float,
    v10: float,
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
    p11: cobj,
    p12: cobj,
    p13: cobj,
    p14: cobj,
    p15: cobj,
    p16: cobj,
    p17: cobj,
    p18: cobj,
    p19: cobj,
    p20: cobj,
):
    # Batch dispatcher for validated buoy updraft-loop shell helpers.
    # kind 0: updraft state; 1: post-conden exit; 2: post-conden thv;
    # kind 3: wu exit; 4: ufrc limiter; 5: next environment load.
    if kind == 0:
        uwshcu_buoy_updraft_state_shell_codon(
            mkx,
            ncnst,
            wtrc_nwset,
            k_fortran,
            flag1,
            flag2,
            v1,
            v2,
            v3,
            v4,
            v5,
            v6,
            v7,
            v8,
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
            p11,
            p12,
            p13,
            p14,
        )
    elif kind == 1:
        uwshcu_scalar_exit_limit_batch_shell_codon(
            1,
            0,
            id_check,
            0.0,
            0.0,
            0.0,
            p1,
            p1,
            p1,
            p1,
            p1,
            p2,
        )
    elif kind == 2:
        uwshcu_conden_exit_thv_batch_shell_codon(
            1,
            0,
            id_check,
            v1,
            v2,
            v3,
            v4,
            v5,
            0.0,
            0.0,
            p1,
            p2,
            p3,
            p3,
        )
    elif kind == 3:
        uwshcu_scalar_exit_limit_batch_shell_codon(
            8,
            0,
            0,
            v1,
            0.0,
            0.0,
            p1,
            p1,
            p1,
            p1,
            p1,
            p2,
        )
    elif kind == 4:
        uwshcu_scalar_exit_limit_batch_shell_codon(
            9,
            k_fortran,
            0,
            v1,
            v2,
            0.0,
            p1,
            p2,
            p3,
            p4,
            p4,
            p5,
        )
    elif kind == 5:
        uwshcu_buoy_next_env_load_shell_codon(
            k_fortran,
            mkx,
            ncnst,
            wtrc_nwset,
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
            p11,
            p12,
            p13,
            p14,
            p15,
            p16,
            p17,
            p18,
            p19,
            p20,
        )


@export
def uwshcu_cin_lcl_init_shell_codon(
    mkx: int,
    zvir: float,
    thj: float,
    qvj: float,
    qlj: float,
    qij: float,
    thv0lcl_p: cobj,
    cin_p: cobj,
    cinlcl_p: cobj,
    plfc_p: cobj,
    klfc_p: cobj,
):
    thv0lcl = Ptr[float](thv0lcl_p)
    cin = Ptr[float](cin_p)
    cinlcl = Ptr[float](cinlcl_p)
    plfc = Ptr[float](plfc_p)
    klfc = Ptr[int](klfc_p)

    thv0lcl[0] = thj * (1.0 + zvir * qvj - qlj - qij)
    cin[0] = 0.0
    cinlcl[0] = 0.0
    plfc[0] = 0.0
    klfc[0] = mkx


@export
def uwshcu_cin_prep_batch_shell_codon(
    kind: int,
    k_fortran: int,
    mkx: int,
    id_check: int,
    plcl: float,
    zvir: float,
    thj: float,
    qvj: float,
    qlj: float,
    qij: float,
    thl0edge: float,
    qt0edge: float,
    ps0_p: cobj,
    p0_p: cobj,
    thl0_p: cobj,
    ssthl0_p: cobj,
    qt0_p: cobj,
    ssqt0_p: cobj,
    klcl_out_p: cobj,
    lcl_exit_code_p: cobj,
    thl0lcl_p: cobj,
    qt0lcl_p: cobj,
    exit_conden_p: cobj,
    exit_code_p: cobj,
    thv_p: cobj,
    thvl_p: cobj,
    cin_p: cobj,
    cinlcl_p: cobj,
    plfc_p: cobj,
    klfc_p: cobj,
):
    if kind == 0:
        uwshcu_conden_exit_thv_batch_shell_codon(
            2,
            k_fortran,
            id_check,
            zvir,
            thj,
            qvj,
            qlj,
            qij,
            thl0edge,
            qt0edge,
            exit_conden_p,
            exit_code_p,
            thv_p,
            thvl_p,
        )
    elif kind == 1:
        uwshcu_lcl_prep_shell_codon(
            mkx,
            plcl,
            ps0_p,
            p0_p,
            thl0_p,
            ssthl0_p,
            qt0_p,
            ssqt0_p,
            klcl_out_p,
            lcl_exit_code_p,
            thl0lcl_p,
            qt0lcl_p,
        )
    elif kind == 2:
        exit_conden = Ptr[float](exit_conden_p)
        exit_code = Ptr[int](exit_code_p)
        exit_code[0] = 0
        if id_check == 1:
            exit_conden[0] = 1.0
            exit_code[0] = 1
    elif kind == 3:
        uwshcu_cin_lcl_init_shell_codon(
            mkx,
            zvir,
            thj,
            qvj,
            qlj,
            qij,
            thv_p,
            cin_p,
            cinlcl_p,
            plfc_p,
            klfc_p,
        )
    elif kind == 4:
        uwshcu_conden_exit_thv_batch_shell_codon(
            1,
            0,
            id_check,
            zvir,
            thj,
            qvj,
            qlj,
            qij,
            0.0,
            0.0,
            exit_conden_p,
            exit_code_p,
            thv_p,
            thvl_p,
        )


@export
def uwshcu_cin_state_save_shell_codon(
    ncnst: int,
    cin_v: float,
    cinlcl_v: float,
    rbuoy: float,
    rkfre: float,
    tkeavg_v: float,
    epsvarw: float,
    kinv_v: int,
    klcl_v: int,
    klfc_v: int,
    plcl_v: float,
    plfc_v: float,
    thvlmin_v: float,
    qtsrc_v: float,
    thvlsrc_v: float,
    thlsrc_v: float,
    usrc_v: float,
    vsrc_v: float,
    thv0lcl_v: float,
    trsrc_p: cobj,
    cin_i_p: cobj,
    cinlcl_i_p: cobj,
    ke_p: cobj,
    kinv_o_p: cobj,
    klcl_o_p: cobj,
    klfc_o_p: cobj,
    plcl_o_p: cobj,
    plfc_o_p: cobj,
    tkeavg_o_p: cobj,
    thvlmin_o_p: cobj,
    qtsrc_o_p: cobj,
    thvlsrc_o_p: cobj,
    thlsrc_o_p: cobj,
    usrc_o_p: cobj,
    vsrc_o_p: cobj,
    thv0lcl_o_p: cobj,
    trsrc_o_p: cobj,
):
    trsrc = Ptr[float](trsrc_p)
    cin_i = Ptr[float](cin_i_p)
    cinlcl_i = Ptr[float](cinlcl_i_p)
    ke = Ptr[float](ke_p)
    kinv_o = Ptr[int](kinv_o_p)
    klcl_o = Ptr[int](klcl_o_p)
    klfc_o = Ptr[int](klfc_o_p)
    plcl_o = Ptr[float](plcl_o_p)
    plfc_o = Ptr[float](plfc_o_p)
    tkeavg_o = Ptr[float](tkeavg_o_p)
    thvlmin_o = Ptr[float](thvlmin_o_p)
    qtsrc_o = Ptr[float](qtsrc_o_p)
    thvlsrc_o = Ptr[float](thvlsrc_o_p)
    thlsrc_o = Ptr[float](thlsrc_o_p)
    usrc_o = Ptr[float](usrc_o_p)
    vsrc_o = Ptr[float](vsrc_o_p)
    thv0lcl_o = Ptr[float](thv0lcl_o_p)
    trsrc_o = Ptr[float](trsrc_o_p)

    cin_i[0] = cin_v
    cinlcl_i[0] = cinlcl_v
    ke[0] = rbuoy / (rkfre * tkeavg_v + epsvarw)
    kinv_o[0] = kinv_v
    klcl_o[0] = klcl_v
    klfc_o[0] = klfc_v
    plcl_o[0] = plcl_v
    plfc_o[0] = plfc_v
    tkeavg_o[0] = tkeavg_v
    thvlmin_o[0] = thvlmin_v
    qtsrc_o[0] = qtsrc_v
    thvlsrc_o[0] = thvlsrc_v
    thlsrc_o[0] = thlsrc_v
    usrc_o[0] = usrc_v
    vsrc_o[0] = vsrc_v
    thv0lcl_o[0] = thv0lcl_v

    m = 0
    while m < ncnst:
        trsrc_o[m] = trsrc[m]
        m += 1


@export
def uwshcu_cin_postcheck_save_shell_codon(
    iter_v: int,
    ncnst: int,
    mkx: int,
    cinlcl_v: float,
    rbuoy: float,
    rkfre: float,
    tkeavg_v: float,
    epsvarw: float,
    kinv_v: int,
    klcl_v: int,
    klfc_v: int,
    plcl_v: float,
    plfc_v: float,
    thvlmin_v: float,
    qtsrc_v: float,
    thvlsrc_v: float,
    thlsrc_v: float,
    usrc_v: float,
    vsrc_v: float,
    thv0lcl_v: float,
    trsrc_p: cobj,
    limit_cin_p: cobj,
    cin_p: cobj,
    klfc_p: cobj,
    exit_code_p: cobj,
    cin_i_p: cobj,
    cinlcl_i_p: cobj,
    ke_p: cobj,
    kinv_o_p: cobj,
    klcl_o_p: cobj,
    klfc_o_p: cobj,
    plcl_o_p: cobj,
    plfc_o_p: cobj,
    tkeavg_o_p: cobj,
    thvlmin_o_p: cobj,
    qtsrc_o_p: cobj,
    thvlsrc_o_p: cobj,
    thlsrc_o_p: cobj,
    usrc_o_p: cobj,
    vsrc_o_p: cobj,
    thv0lcl_o_p: cobj,
    trsrc_o_p: cobj,
):
    limit_cin = Ptr[float](limit_cin_p)
    cin = Ptr[float](cin_p)
    klfc = Ptr[int](klfc_p)
    exit_code = Ptr[int](exit_code_p)

    exit_code[0] = 0
    if cin[0] < 0.0:
        limit_cin[0] = 1.0
        cin[0] = 0.0

    if klfc_v >= mkx:
        klfc[0] = mkx
        exit_code[0] = 1
        return

    klfc[0] = klfc_v

    if iter_v == 1:
        uwshcu_cin_state_save_shell_codon(
            ncnst,
            cin[0],
            cinlcl_v,
            rbuoy,
            rkfre,
            tkeavg_v,
            epsvarw,
            kinv_v,
            klcl_v,
            klfc_v,
            plcl_v,
            plfc_v,
            thvlmin_v,
            qtsrc_v,
            thvlsrc_v,
            thlsrc_v,
            usrc_v,
            vsrc_v,
            thv0lcl_v,
            trsrc_p,
            cin_i_p,
            cinlcl_i_p,
            ke_p,
            kinv_o_p,
            klcl_o_p,
            klfc_o_p,
            plcl_o_p,
            plfc_o_p,
            tkeavg_o_p,
            thvlmin_o_p,
            qtsrc_o_p,
            thvlsrc_o_p,
            thlsrc_o_p,
            usrc_o_p,
            vsrc_o_p,
            thv0lcl_o_p,
            trsrc_o_p,
        )


@export
def uwshcu_cin_state_restore_shell_codon(
    ncnst: int,
    use_cincin: int,
    cin_i_v: float,
    cinlcl_i_v: float,
    alpha: float,
    del_cin: float,
    del_cinlcl: float,
    kinv_o_v: int,
    klcl_o_v: int,
    klfc_o_v: int,
    plcl_o_v: float,
    plfc_o_v: float,
    tkeavg_o_v: float,
    thvlmin_o_v: float,
    qtsrc_o_v: float,
    thvlsrc_o_v: float,
    thlsrc_o_v: float,
    usrc_o_v: float,
    vsrc_o_v: float,
    thv0lcl_o_v: float,
    trsrc_o_p: cobj,
    cin_p: cobj,
    cinlcl_p: cobj,
    kinv_p: cobj,
    klcl_p: cobj,
    klfc_p: cobj,
    plcl_p: cobj,
    plfc_p: cobj,
    tkeavg_p: cobj,
    thvlmin_p: cobj,
    qtsrc_p: cobj,
    thvlsrc_p: cobj,
    thlsrc_p: cobj,
    usrc_p: cobj,
    vsrc_p: cobj,
    thv0lcl_p: cobj,
    trsrc_p: cobj,
):
    trsrc_o = Ptr[float](trsrc_o_p)
    cin = Ptr[float](cin_p)
    cinlcl = Ptr[float](cinlcl_p)
    kinv = Ptr[int](kinv_p)
    klcl = Ptr[int](klcl_p)
    klfc = Ptr[int](klfc_p)
    plcl = Ptr[float](plcl_p)
    plfc = Ptr[float](plfc_p)
    tkeavg = Ptr[float](tkeavg_p)
    thvlmin = Ptr[float](thvlmin_p)
    qtsrc = Ptr[float](qtsrc_p)
    thvlsrc = Ptr[float](thvlsrc_p)
    thlsrc = Ptr[float](thlsrc_p)
    usrc = Ptr[float](usrc_p)
    vsrc = Ptr[float](vsrc_p)
    thv0lcl = Ptr[float](thv0lcl_p)
    trsrc = Ptr[float](trsrc_p)

    cin[0] = cin_i_v + alpha * del_cin
    if use_cincin != 0:
        cinlcl[0] = cinlcl_i_v
    else:
        cinlcl[0] = cinlcl_i_v + alpha * del_cinlcl

    kinv[0] = kinv_o_v
    klcl[0] = klcl_o_v
    klfc[0] = klfc_o_v
    plcl[0] = plcl_o_v
    plfc[0] = plfc_o_v
    tkeavg[0] = tkeavg_o_v
    thvlmin[0] = thvlmin_o_v
    qtsrc[0] = qtsrc_o_v
    thvlsrc[0] = thvlsrc_o_v
    thlsrc[0] = thlsrc_o_v
    usrc[0] = usrc_o_v
    vsrc[0] = vsrc_o_v
    thv0lcl[0] = thv0lcl_o_v

    m = 0
    while m < ncnst:
        trsrc[m] = trsrc_o[m]
        m += 1


@export
def uwshcu_release_level_shell_codon(
    kinv: int,
    klcl: int,
    plcl: float,
    thv0lcl_v: float,
    ps0_p: cobj,
    thv0bot_p: cobj,
    krel_p: cobj,
    prel_p: cobj,
    thv0rel_p: cobj,
):
    ps0 = Ptr[float](ps0_p)
    thv0bot = Ptr[float](thv0bot_p)
    krel_out = Ptr[int](krel_p)
    prel = Ptr[float](prel_p)
    thv0rel = Ptr[float](thv0rel_p)

    if klcl < kinv:
        krel = kinv
        prel[0] = ps0[krel - 1]
        thv0rel[0] = thv0bot[krel - 1]
    else:
        krel = klcl
        prel[0] = plcl
        thv0rel[0] = thv0lcl_v
    krel_out[0] = krel


@export
def uwshcu_release_base_shell_codon(
    mkx: int,
    kinv: int,
    krel: int,
    cbmf: float,
    wrel: float,
    winv: float,
    ufrcinv: float,
    ufrclcl: float,
    thlsrc: float,
    qtsrc: float,
    prel: float,
    ps0_p: cobj,
    ufrc_p: cobj,
    umf_p: cobj,
    wu_p: cobj,
    emf_p: cobj,
    thlu_p: cobj,
    qtu_p: cobj,
    ufrcinvbase_p: cobj,
    winvbase_p: cobj,
    pe_p: cobj,
    dpe_p: cobj,
):
    ps0 = Ptr[float](ps0_p)
    ufrc = Ptr[float](ufrc_p)
    umf = Ptr[float](umf_p)
    wu = Ptr[float](wu_p)
    emf = Ptr[float](emf_p)
    thlu = Ptr[float](thlu_p)
    qtu = Ptr[float](qtu_p)
    ufrcinvbase = Ptr[float](ufrcinvbase_p)
    winvbase = Ptr[float](winvbase_p)
    pe = Ptr[float](pe_p)
    dpe = Ptr[float](dpe_p)

    km1_rel = krel - 1
    ufrc[km1_rel] = ufrclcl
    ufrcinvbase[0] = ufrcinv
    winvbase[0] = winv

    k = kinv - 1
    while k <= km1_rel:
        umf[k] = cbmf
        wu[k] = winv
        k += 1

    emf[km1_rel] = 0.0
    umf[km1_rel] = cbmf
    wu[km1_rel] = wrel
    thlu[km1_rel] = thlsrc
    qtu[km1_rel] = qtsrc
    pe[0] = 0.5 * (prel + ps0[krel])
    dpe[0] = prel - ps0[krel]


@export
def uwshcu_release_scaleh_batch_shell_codon(
    kind: int,
    mkx: int,
    kinv: int,
    klcl: int,
    krel: int,
    code_in: int,
    plcl: float,
    thv0lcl_v: float,
    mu_v: float,
    mumin2_v: float,
    mumin0_v: float,
    mumin1_v: float,
    wtw_v: float,
    ufrclcl_v: float,
    cbmf: float,
    wrel: float,
    winv: float,
    ufrcinv: float,
    thlsrc: float,
    qtsrc: float,
    prel_v: float,
    ppen_v: float,
    dp0_kpen: float,
    ps0_p: cobj,
    thv0bot_p: cobj,
    krel_p: cobj,
    prel_p: cobj,
    thv0rel_p: cobj,
    limit_ufrc_p: cobj,
    limit_cbmf_p: cobj,
    limit_ppen_p: cobj,
    exit_wtw_p: cobj,
    exit_ufrc_p: cobj,
    exit_conden_p: cobj,
    exit_cufilter_p: cobj,
    exit_code_p: cobj,
    ufrc_p: cobj,
    umf_p: cobj,
    wu_p: cobj,
    emf_p: cobj,
    thlu_p: cobj,
    qtu_p: cobj,
    ufrcinvbase_p: cobj,
    winvbase_p: cobj,
    pe_p: cobj,
    dpe_p: cobj,
):
    if kind == 0:
        uwshcu_release_level_shell_codon(
            kinv,
            klcl,
            plcl,
            thv0lcl_v,
            ps0_p,
            thv0bot_p,
            krel_p,
            prel_p,
            thv0rel_p,
        )
    elif kind == 1:
        exit_code = Ptr[int](exit_code_p)
        exit_code[0] = 0
        if mu_v >= 3.0:
            exit_code[0] = 1
    elif kind == 2:
        limit_ufrc = Ptr[float](limit_ufrc_p)
        if mu_v == mumin2_v:
            limit_ufrc[0] = 1.0
    elif kind == 3:
        limit_cbmf = Ptr[float](limit_cbmf_p)
        limit_ufrc = Ptr[float](limit_ufrc_p)
        if mu_v == mumin0_v:
            limit_cbmf[0] = 1.0
        if mu_v == mumin1_v:
            limit_ufrc[0] = 1.0
    elif kind == 4:
        exit_wtw = Ptr[float](exit_wtw_p)
        exit_code = Ptr[int](exit_code_p)
        exit_code[0] = 0
        if wtw_v <= 0.0:
            exit_wtw[0] = 1.0
            exit_code[0] = 1
    elif kind == 5:
        exit_ufrc = Ptr[float](exit_ufrc_p)
        exit_code = Ptr[int](exit_code_p)
        exit_code[0] = 0
        if ufrclcl_v <= 0.0001:
            exit_ufrc[0] = 1.0
            exit_code[0] = 1
    elif kind == 6:
        uwshcu_release_base_shell_codon(
            mkx,
            kinv,
            krel,
            cbmf,
            wrel,
            winv,
            ufrcinv,
            ufrclcl_v,
            thlsrc,
            qtsrc,
            prel_v,
            ps0_p,
            ufrc_p,
            umf_p,
            wu_p,
            emf_p,
            thlu_p,
            qtu_p,
            ufrcinvbase_p,
            winvbase_p,
            pe_p,
            dpe_p,
        )
    elif kind == 7:
        exit_conden = Ptr[float](exit_conden_p)
        exit_code = Ptr[int](exit_code_p)
        exit_code[0] = 0
        if code_in == 1:
            exit_conden[0] = 1.0
            exit_code[0] = 1
    elif kind == 8:
        limit_ppen = Ptr[float](limit_ppen_p)
        if ppen_v == -dp0_kpen or ppen_v == 0.0:
            limit_ppen[0] = 1.0
    elif kind == 9:
        exit_cufilter = Ptr[float](exit_cufilter_p)
        exit_code = Ptr[int](exit_code_p)
        exit_code[0] = 0
        if code_in != 0:
            exit_cufilter[0] = 1.0
            exit_code[0] = 1


@export
def uwshcu_release_conden_exit_shell_codon(
    id_check: int,
    exit_conden_p: cobj,
    exit_code_p: cobj,
):
    exit_conden = Ptr[float](exit_conden_p)
    exit_code = Ptr[int](exit_code_p)

    exit_code[0] = 0
    if id_check == 1:
        exit_conden[0] = 1.0
        exit_code[0] = 1


@export
def uwshcu_release_env_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    kinv: int,
    krel: int,
    zvir: float,
    pgfc: float,
    usrc: float,
    vsrc: float,
    prel: float,
    pe_v: float,
    thv0rel_v: float,
    thj: float,
    qvj: float,
    qlj: float,
    qij: float,
    ps0_p: cobj,
    p0_p: cobj,
    thl0_p: cobj,
    ssthl0_p: cobj,
    qt0_p: cobj,
    ssqt0_p: cobj,
    u0_p: cobj,
    ssu0_p: cobj,
    v0_p: cobj,
    ssv0_p: cobj,
    tr0_p: cobj,
    sstr0_p: cobj,
    wt0_p: cobj,
    sswt0_p: cobj,
    trsrc_p: cobj,
    wtsrc_p: cobj,
    thvu_p: cobj,
    uu_p: cobj,
    vu_p: cobj,
    tru_p: cobj,
    wtu_p: cobj,
    thvebot_p: cobj,
    thle_p: cobj,
    qte_p: cobj,
    ue_p: cobj,
    ve_p: cobj,
    tre_p: cobj,
    wte_p: cobj,
):
    ps0 = Ptr[float](ps0_p)
    p0 = Ptr[float](p0_p)
    thl0 = Ptr[float](thl0_p)
    ssthl0 = Ptr[float](ssthl0_p)
    qt0 = Ptr[float](qt0_p)
    ssqt0 = Ptr[float](ssqt0_p)
    u0 = Ptr[float](u0_p)
    ssu0 = Ptr[float](ssu0_p)
    v0 = Ptr[float](v0_p)
    ssv0 = Ptr[float](ssv0_p)
    tr0 = Ptr[float](tr0_p)
    sstr0 = Ptr[float](sstr0_p)
    wt0 = Ptr[float](wt0_p)
    sswt0 = Ptr[float](sswt0_p)
    trsrc = Ptr[float](trsrc_p)
    wtsrc = Ptr[float](wtsrc_p)
    thvu = Ptr[float](thvu_p)
    uu = Ptr[float](uu_p)
    vu = Ptr[float](vu_p)
    tru = Ptr[float](tru_p)
    wtu = Ptr[float](wtu_p)
    thvebot = Ptr[float](thvebot_p)
    thle = Ptr[float](thle_p)
    qte = Ptr[float](qte_p)
    ue = Ptr[float](ue_p)
    ve = Ptr[float](ve_p)
    tre = Ptr[float](tre_p)
    wte = Ptr[float](wte_p)

    km1_rel = krel - 1
    thvu[km1_rel] = thj * (1.0 + zvir * qvj - qlj - qij)

    uplus = 0.0
    vplus = 0.0
    if krel == kinv:
        uplus = pgfc * ssu0[kinv - 1] * (prel - ps0[kinv - 1])
        vplus = pgfc * ssv0[kinv - 1] * (prel - ps0[kinv - 1])
    else:
        k = kinv
        while k <= krel - 1:
            uplus = uplus + pgfc * ssu0[k - 1] * (ps0[k] - ps0[k - 1])
            vplus = vplus + pgfc * ssv0[k - 1] * (ps0[k] - ps0[k - 1])
            k += 1
        uplus = uplus + pgfc * ssu0[krel - 1] * (prel - ps0[krel - 1])
        vplus = vplus + pgfc * ssv0[krel - 1] * (prel - ps0[krel - 1])

    uu[km1_rel] = usrc + uplus
    vu[km1_rel] = vsrc + vplus

    m = 0
    iface_stride = mkx + 1
    while m < ncnst:
        tru[km1_rel + m * iface_stride] = trsrc[m]
        m += 1

    m = 0
    while m < wtrc_nwset:
        wtu[km1_rel + m * iface_stride] = wtsrc[m]
        m += 1

    layer_idx = krel - 1
    thvebot[0] = thv0rel_v
    thle[0] = thl0[layer_idx] + ssthl0[layer_idx] * (pe_v - p0[layer_idx])
    qte[0] = qt0[layer_idx] + ssqt0[layer_idx] * (pe_v - p0[layer_idx])
    ue[0] = u0[layer_idx] + ssu0[layer_idx] * (pe_v - p0[layer_idx])
    ve[0] = v0[layer_idx] + ssv0[layer_idx] * (pe_v - p0[layer_idx])

    m = 0
    while m < ncnst:
        idx = layer_idx + m * mkx
        tre[m] = tr0[idx] + sstr0[idx] * (pe_v - p0[layer_idx])
        m += 1

    m = 0
    while m < wtrc_nwset:
        idx = layer_idx + m * mkx
        wte[m] = wt0[idx] + sswt0[idx] * (pe_v - p0[layer_idx])
        m += 1


@export
def uwshcu_scaleh_set_codon(
    tscaleh: float,
    scaleh_p: cobj,
):
    scaleh = Ptr[float](scaleh_p)

    scaleh[0] = tscaleh
    if tscaleh < 0.0:
        scaleh[0] = 1000.0


@export
def uwshcu_scaleh_iter_init_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    iter_scaleh: int,
    krel: int,
    tscaleh: float,
    wlcl: float,
    prel: float,
    thv0rel_v: float,
    scaleh_p: cobj,
    ps0_p: cobj,
    p0_p: cobj,
    thl0_p: cobj,
    ssthl0_p: cobj,
    qt0_p: cobj,
    ssqt0_p: cobj,
    u0_p: cobj,
    ssu0_p: cobj,
    v0_p: cobj,
    ssv0_p: cobj,
    tr0_p: cobj,
    sstr0_p: cobj,
    wt0_p: cobj,
    sswt0_p: cobj,
    kbup_p: cobj,
    kpen_p: cobj,
    wtw_p: cobj,
    pe_p: cobj,
    dpe_p: cobj,
    thvebot_p: cobj,
    thle_p: cobj,
    qte_p: cobj,
    ue_p: cobj,
    ve_p: cobj,
    tre_p: cobj,
    wte_p: cobj,
):
    scaleh = Ptr[float](scaleh_p)
    ps0 = Ptr[float](ps0_p)
    p0 = Ptr[float](p0_p)
    thl0 = Ptr[float](thl0_p)
    ssthl0 = Ptr[float](ssthl0_p)
    qt0 = Ptr[float](qt0_p)
    ssqt0 = Ptr[float](ssqt0_p)
    u0 = Ptr[float](u0_p)
    ssu0 = Ptr[float](ssu0_p)
    v0 = Ptr[float](v0_p)
    ssv0 = Ptr[float](ssv0_p)
    tr0 = Ptr[float](tr0_p)
    sstr0 = Ptr[float](sstr0_p)
    wt0 = Ptr[float](wt0_p)
    sswt0 = Ptr[float](sswt0_p)
    kbup = Ptr[int](kbup_p)
    kpen = Ptr[int](kpen_p)
    wtw = Ptr[float](wtw_p)
    pe = Ptr[float](pe_p)
    dpe = Ptr[float](dpe_p)
    thvebot = Ptr[float](thvebot_p)
    thle = Ptr[float](thle_p)
    qte = Ptr[float](qte_p)
    ue = Ptr[float](ue_p)
    ve = Ptr[float](ve_p)
    tre = Ptr[float](tre_p)
    wte = Ptr[float](wte_p)

    if iter_scaleh == 1:
        scaleh[0] = tscaleh
        if tscaleh < 0.0:
            scaleh[0] = 1000.0

    kbup[0] = krel
    kpen[0] = krel
    wtw[0] = wlcl * wlcl
    pe[0] = 0.5 * (prel + ps0[krel])
    dpe[0] = prel - ps0[krel]
    thvebot[0] = thv0rel_v

    layer_idx = krel - 1
    thle[0] = thl0[layer_idx] + ssthl0[layer_idx] * (pe[0] - p0[layer_idx])
    qte[0] = qt0[layer_idx] + ssqt0[layer_idx] * (pe[0] - p0[layer_idx])
    ue[0] = u0[layer_idx] + ssu0[layer_idx] * (pe[0] - p0[layer_idx])
    ve[0] = v0[layer_idx] + ssv0[layer_idx] * (pe[0] - p0[layer_idx])

    m = 0
    while m < ncnst:
        idx = layer_idx + m * mkx
        tre[m] = tr0[idx] + sstr0[idx] * (pe[0] - p0[layer_idx])
        m += 1

    m = 0
    while m < wtrc_nwset:
        idx = layer_idx + m * mkx
        wte[m] = wt0[idx] + sswt0[idx] * (pe[0] - p0[layer_idx])
        m += 1


@export
def uwshcu_release_env_scaleh_iter_init_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    iter_scaleh: int,
    kinv: int,
    krel: int,
    zvir: float,
    pgfc: float,
    usrc: float,
    vsrc: float,
    prel: float,
    pe_v: float,
    thv0rel_v: float,
    thj: float,
    qvj: float,
    qlj: float,
    qij: float,
    tscaleh: float,
    wlcl: float,
    ps0_p: cobj,
    p0_p: cobj,
    thl0_p: cobj,
    ssthl0_p: cobj,
    qt0_p: cobj,
    ssqt0_p: cobj,
    u0_p: cobj,
    ssu0_p: cobj,
    v0_p: cobj,
    ssv0_p: cobj,
    tr0_p: cobj,
    sstr0_p: cobj,
    wt0_p: cobj,
    sswt0_p: cobj,
    trsrc_p: cobj,
    wtsrc_p: cobj,
    scaleh_p: cobj,
    thvu_p: cobj,
    uu_p: cobj,
    vu_p: cobj,
    tru_p: cobj,
    wtu_p: cobj,
    kbup_p: cobj,
    kpen_p: cobj,
    wtw_p: cobj,
    pe_p: cobj,
    dpe_p: cobj,
    thvebot_p: cobj,
    thle_p: cobj,
    qte_p: cobj,
    ue_p: cobj,
    ve_p: cobj,
    tre_p: cobj,
    wte_p: cobj,
):
    if iter_scaleh == 1:
        uwshcu_release_env_shell_codon(
            mkx,
            ncnst,
            wtrc_nwset,
            kinv,
            krel,
            zvir,
            pgfc,
            usrc,
            vsrc,
            prel,
            pe_v,
            thv0rel_v,
            thj,
            qvj,
            qlj,
            qij,
            ps0_p,
            p0_p,
            thl0_p,
            ssthl0_p,
            qt0_p,
            ssqt0_p,
            u0_p,
            ssu0_p,
            v0_p,
            ssv0_p,
            tr0_p,
            sstr0_p,
            wt0_p,
            sswt0_p,
            trsrc_p,
            wtsrc_p,
            thvu_p,
            uu_p,
            vu_p,
            tru_p,
            wtu_p,
            thvebot_p,
            thle_p,
            qte_p,
            ue_p,
            ve_p,
            tre_p,
            wte_p,
        )

    uwshcu_scaleh_iter_init_shell_codon(
        mkx,
        ncnst,
        wtrc_nwset,
        iter_scaleh,
        krel,
        tscaleh,
        wlcl,
        prel,
        thv0rel_v,
        scaleh_p,
        ps0_p,
        p0_p,
        thl0_p,
        ssthl0_p,
        qt0_p,
        ssqt0_p,
        u0_p,
        ssu0_p,
        v0_p,
        ssv0_p,
        tr0_p,
        sstr0_p,
        wt0_p,
        sswt0_p,
        kbup_p,
        kpen_p,
        wtw_p,
        pe_p,
        dpe_p,
        thvebot_p,
        thle_p,
        qte_p,
        ue_p,
        ve_p,
        tre_p,
        wte_p,
    )


@export
def uwshcu_penent_prep_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    kbup: int,
    kpen: int,
    r: float,
    g: float,
    dt: float,
    rpen: float,
    ppen: float,
    ps0_p: cobj,
    p0_p: cobj,
    dp0_p: cobj,
    thv0bot_p: cobj,
    thv0top_p: cobj,
    exns0_p: cobj,
    thl0_p: cobj,
    ssthl0_p: cobj,
    qt0_p: cobj,
    ssqt0_p: cobj,
    u0_p: cobj,
    ssu0_p: cobj,
    v0_p: cobj,
    ssv0_p: cobj,
    tr0_p: cobj,
    sstr0_p: cobj,
    wt0_p: cobj,
    sswt0_p: cobj,
    umf_p: cobj,
    emf_p: cobj,
    ufrc_p: cobj,
    dwten_p: cobj,
    diten_p: cobj,
    fer_p: cobj,
    fdr_p: cobj,
    rei_p: cobj,
    thlu_p: cobj,
    qtu_p: cobj,
    uu_p: cobj,
    vu_p: cobj,
    tru_p: cobj,
    wtu_p: cobj,
    thlu_emf_p: cobj,
    qtu_emf_p: cobj,
    uu_emf_p: cobj,
    vu_emf_p: cobj,
    tru_emf_p: cobj,
    wtu_emf_p: cobj,
    wtdwten_p: cobj,
    wtditen_p: cobj,
    limit_emf_p: cobj,
):
    ps0 = Ptr[float](ps0_p)
    p0 = Ptr[float](p0_p)
    dp0 = Ptr[float](dp0_p)
    thv0bot = Ptr[float](thv0bot_p)
    thv0top = Ptr[float](thv0top_p)
    exns0 = Ptr[float](exns0_p)
    thl0 = Ptr[float](thl0_p)
    ssthl0 = Ptr[float](ssthl0_p)
    qt0 = Ptr[float](qt0_p)
    ssqt0 = Ptr[float](ssqt0_p)
    u0 = Ptr[float](u0_p)
    ssu0 = Ptr[float](ssu0_p)
    v0 = Ptr[float](v0_p)
    ssv0 = Ptr[float](ssv0_p)
    tr0 = Ptr[float](tr0_p)
    sstr0 = Ptr[float](sstr0_p)
    wt0 = Ptr[float](wt0_p)
    sswt0 = Ptr[float](sswt0_p)
    umf = Ptr[float](umf_p)
    emf = Ptr[float](emf_p)
    ufrc = Ptr[float](ufrc_p)
    dwten = Ptr[float](dwten_p)
    diten = Ptr[float](diten_p)
    fer = Ptr[float](fer_p)
    fdr = Ptr[float](fdr_p)
    rei = Ptr[float](rei_p)
    thlu = Ptr[float](thlu_p)
    qtu = Ptr[float](qtu_p)
    uu = Ptr[float](uu_p)
    vu = Ptr[float](vu_p)
    tru = Ptr[float](tru_p)
    wtu = Ptr[float](wtu_p)
    thlu_emf = Ptr[float](thlu_emf_p)
    qtu_emf = Ptr[float](qtu_emf_p)
    uu_emf = Ptr[float](uu_emf_p)
    vu_emf = Ptr[float](vu_emf_p)
    tru_emf = Ptr[float](tru_emf_p)
    wtu_emf = Ptr[float](wtu_emf_p)
    wtdwten = Ptr[float](wtdwten_p)
    wtditen = Ptr[float](wtditen_p)
    limit_emf = Ptr[float](limit_emf_p)

    iface_stride = mkx + 1

    k = kpen
    while k <= mkx:
        umf[k] = 0.0
        emf[k] = 0.0
        ufrc[k] = 0.0
        k += 1

    k = kpen + 1
    while k <= mkx:
        layer_idx = k - 1
        dwten[layer_idx] = 0.0
        diten[layer_idx] = 0.0
        fer[layer_idx] = 0.0
        fdr[layer_idx] = 0.0
        m = 0
        while m < wtrc_nwset:
            wt_idx = layer_idx + m * mkx
            wtdwten[wt_idx] = 0.0
            wtditen[wt_idx] = 0.0
            m += 1
        k += 1

    k = 0
    while k <= mkx:
        thlu_emf[k] = thlu[k]
        qtu_emf[k] = qtu[k]
        uu_emf[k] = uu[k]
        vu_emf[k] = vu[k]
        m = 0
        while m < ncnst:
            idx = k + m * iface_stride
            tru_emf[idx] = tru[idx]
            m += 1
        m = 0
        while m < wtrc_nwset:
            idx = k + m * iface_stride
            wtu_emf[idx] = wtu[idx]
            m += 1
        k += 1

    k = kpen - 1
    while k >= kbup:
        rhos0j = ps0[k] / (r * 0.5 * (thv0bot[k] + thv0top[k - 1]) * exns0[k])

        if k == kpen - 1:
            emf_trial = umf[k] * ppen * rei[kpen - 1] * rpen
            if emf_trial < -0.1 * rhos0j:
                limit_emf[0] = 1.0
            if emf_trial < -0.9 * dp0[kpen - 1] / g / dt:
                limit_emf[0] = 1.0

            emf[k] = max(max(emf_trial, -0.1 * rhos0j), -0.9 * dp0[kpen - 1] / g / dt)
            layer_idx = kpen - 1
            p_delta = ps0[k] - p0[layer_idx]
            thlu_emf[k] = thl0[layer_idx] + ssthl0[layer_idx] * p_delta
            qtu_emf[k] = qt0[layer_idx] + ssqt0[layer_idx] * p_delta
            uu_emf[k] = u0[layer_idx] + ssu0[layer_idx] * p_delta
            vu_emf[k] = v0[layer_idx] + ssv0[layer_idx] * p_delta

            m = 0
            while m < ncnst:
                idx_layer = layer_idx + m * mkx
                idx_iface = k + m * iface_stride
                tru_emf[idx_iface] = tr0[idx_layer] + sstr0[idx_layer] * p_delta
                m += 1

            m = 0
            while m < wtrc_nwset:
                idx_layer = layer_idx + m * mkx
                idx_iface = k + m * iface_stride
                wtu_emf[idx_iface] = wt0[idx_layer] + sswt0[idx_layer] * p_delta
                m += 1
        else:
            layer_idx = k
            emf_trial = emf[k + 1] - umf[k] * dp0[layer_idx] * rei[layer_idx] * rpen
            if emf_trial < -0.1 * rhos0j:
                limit_emf[0] = 1.0
            if emf_trial < -0.9 * dp0[layer_idx] / g / dt:
                limit_emf[0] = 1.0
            emf[k] = max(max(emf_trial, -0.1 * rhos0j), -0.9 * dp0[layer_idx] / g / dt)

            if abs(emf[k]) > abs(emf[k + 1]):
                emf_delta = emf[k] - emf[k + 1]
                thlu_emf[k] = (thlu_emf[k + 1] * emf[k + 1] + thl0[layer_idx] * emf_delta) / emf[k]
                qtu_emf[k] = (qtu_emf[k + 1] * emf[k + 1] + qt0[layer_idx] * emf_delta) / emf[k]
                uu_emf[k] = (uu_emf[k + 1] * emf[k + 1] + u0[layer_idx] * emf_delta) / emf[k]
                vu_emf[k] = (vu_emf[k + 1] * emf[k + 1] + v0[layer_idx] * emf_delta) / emf[k]

                m = 0
                while m < ncnst:
                    idx_layer = layer_idx + m * mkx
                    idx_iface = k + m * iface_stride
                    idx_next = k + 1 + m * iface_stride
                    tru_emf[idx_iface] = (tru_emf[idx_next] * emf[k + 1] + tr0[idx_layer] * emf_delta) / emf[k]
                    m += 1

                m = 0
                while m < wtrc_nwset:
                    idx_layer = layer_idx + m * mkx
                    idx_iface = k + m * iface_stride
                    idx_next = k + 1 + m * iface_stride
                    wtu_emf[idx_iface] = (wtu_emf[idx_next] * emf[k + 1] + wt0[idx_layer] * emf_delta) / emf[k]
                    m += 1
            else:
                thlu_emf[k] = thl0[layer_idx]
                qtu_emf[k] = qt0[layer_idx]
                uu_emf[k] = u0[layer_idx]
                vu_emf[k] = v0[layer_idx]

                m = 0
                while m < ncnst:
                    tru_emf[k + m * iface_stride] = tr0[layer_idx + m * mkx]
                    m += 1

                m = 0
                while m < wtrc_nwset:
                    wtu_emf[k + m * iface_stride] = wt0[layer_idx + m * mkx]
                    m += 1

        k -= 1


def _uwshcu_fluxbelowinv_value(
    cbmf: float,
    g_v: float,
    dt_v: float,
    ps0: Ptr[float],
    kinv: int,
    xsrc: float,
    xmean: float,
    xtop_in: float,
    xbot_in: float,
    k: int,
) -> float:
    dp = ps0[kinv - 1] - ps0[kinv]
    xbot = xbot_in
    xtop = xtop_in

    xtop_ori = xtop
    xbot_ori = xbot
    rcbmf = (cbmf * g_v * dt_v) / dp

    if xbot >= xtop:
        rpeff = (xmean - xtop) / max(1.0e-20, xbot - xtop)
    else:
        rpeff = (xmean - xtop) / min(-1.0e-20, xbot - xtop)

    rpeff = min(max(0.0, rpeff), 1.0)
    if rpeff == 0.0 or rpeff == 1.0:
        xbot = xmean
        xtop = xmean

    rr = rpeff / rcbmf
    pinv = ps0[kinv - 1] - rpeff * dp
    xflx = cbmf * (xsrc - xbot) * (ps0[0] - ps0[k]) / (ps0[0] - pinv)
    if k == kinv - 1 and rr <= 1.0:
        xflx = xflx - (1.0 - rr) * cbmf * (xtop_ori - xbot_ori)
    return xflx


@export
def uwshcu_turbulent_flux_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    kinv: int,
    krel: int,
    kbup: int,
    kpen: int,
    use_momenflx: int,
    cbmf: float,
    g_v: float,
    dt_v: float,
    cp: float,
    pgfc: float,
    qtsrc: float,
    thlsrc: float,
    usrc: float,
    vsrc: float,
    ps0_p: cobj,
    p0_p: cobj,
    exns0_p: cobj,
    qt0_p: cobj,
    ssqt0_p: cobj,
    thl0_p: cobj,
    ssthl0_p: cobj,
    u0_p: cobj,
    ssu0_p: cobj,
    v0_p: cobj,
    ssv0_p: cobj,
    tr0_p: cobj,
    sstr0_p: cobj,
    wt0_p: cobj,
    sswt0_p: cobj,
    trsrc_p: cobj,
    wtsrc_p: cobj,
    umf_p: cobj,
    emf_p: cobj,
    thlu_p: cobj,
    qtu_p: cobj,
    uu_p: cobj,
    vu_p: cobj,
    tru_p: cobj,
    wtu_p: cobj,
    thlu_emf_p: cobj,
    qtu_emf_p: cobj,
    uu_emf_p: cobj,
    vu_emf_p: cobj,
    tru_emf_p: cobj,
    wtu_emf_p: cobj,
    slflx_p: cobj,
    qtflx_p: cobj,
    uflx_p: cobj,
    vflx_p: cobj,
    trflx_p: cobj,
    wtflx_p: cobj,
):
    ps0 = Ptr[float](ps0_p)
    p0 = Ptr[float](p0_p)
    exns0 = Ptr[float](exns0_p)
    qt0 = Ptr[float](qt0_p)
    ssqt0 = Ptr[float](ssqt0_p)
    thl0 = Ptr[float](thl0_p)
    ssthl0 = Ptr[float](ssthl0_p)
    u0 = Ptr[float](u0_p)
    ssu0 = Ptr[float](ssu0_p)
    v0 = Ptr[float](v0_p)
    ssv0 = Ptr[float](ssv0_p)
    tr0 = Ptr[float](tr0_p)
    sstr0 = Ptr[float](sstr0_p)
    wt0 = Ptr[float](wt0_p)
    sswt0 = Ptr[float](sswt0_p)
    trsrc = Ptr[float](trsrc_p)
    wtsrc = Ptr[float](wtsrc_p)
    umf = Ptr[float](umf_p)
    emf = Ptr[float](emf_p)
    thlu = Ptr[float](thlu_p)
    qtu = Ptr[float](qtu_p)
    uu = Ptr[float](uu_p)
    vu = Ptr[float](vu_p)
    tru = Ptr[float](tru_p)
    wtu = Ptr[float](wtu_p)
    thlu_emf = Ptr[float](thlu_emf_p)
    qtu_emf = Ptr[float](qtu_emf_p)
    uu_emf = Ptr[float](uu_emf_p)
    vu_emf = Ptr[float](vu_emf_p)
    tru_emf = Ptr[float](tru_emf_p)
    wtu_emf = Ptr[float](wtu_emf_p)
    slflx = Ptr[float](slflx_p)
    qtflx = Ptr[float](qtflx_p)
    uflx = Ptr[float](uflx_p)
    vflx = Ptr[float](vflx_p)
    trflx = Ptr[float](trflx_p)
    wtflx = Ptr[float](wtflx_p)

    iface_stride = mkx + 1

    xmean = qt0[kinv - 1]
    xtop = qt0[kinv] + ssqt0[kinv] * (ps0[kinv] - p0[kinv])
    xbot = qt0[kinv - 2] + ssqt0[kinv - 2] * (ps0[kinv - 1] - p0[kinv - 2])
    k = 0
    while k <= kinv - 1:
        qtflx[k] = _uwshcu_fluxbelowinv_value(cbmf, g_v, dt_v, ps0, kinv, qtsrc, xmean, xtop, xbot, k)
        k += 1

    xmean = thl0[kinv - 1]
    xtop = thl0[kinv] + ssthl0[kinv] * (ps0[kinv] - p0[kinv])
    xbot = thl0[kinv - 2] + ssthl0[kinv - 2] * (ps0[kinv - 1] - p0[kinv - 2])
    k = 0
    while k <= kinv - 1:
        xflx = _uwshcu_fluxbelowinv_value(cbmf, g_v, dt_v, ps0, kinv, thlsrc, xmean, xtop, xbot, k)
        slflx[k] = cp * exns0[k] * xflx
        k += 1

    xmean = u0[kinv - 1]
    xtop = u0[kinv] + ssu0[kinv] * (ps0[kinv] - p0[kinv])
    xbot = u0[kinv - 2] + ssu0[kinv - 2] * (ps0[kinv - 1] - p0[kinv - 2])
    k = 0
    while k <= kinv - 1:
        uflx[k] = _uwshcu_fluxbelowinv_value(cbmf, g_v, dt_v, ps0, kinv, usrc, xmean, xtop, xbot, k)
        k += 1

    xmean = v0[kinv - 1]
    xtop = v0[kinv] + ssv0[kinv] * (ps0[kinv] - p0[kinv])
    xbot = v0[kinv - 2] + ssv0[kinv - 2] * (ps0[kinv - 1] - p0[kinv - 2])
    k = 0
    while k <= kinv - 1:
        vflx[k] = _uwshcu_fluxbelowinv_value(cbmf, g_v, dt_v, ps0, kinv, vsrc, xmean, xtop, xbot, k)
        k += 1

    m = 0
    while m < ncnst:
        xmean = tr0[kinv - 1 + m * mkx]
        xtop = tr0[kinv + m * mkx] + sstr0[kinv + m * mkx] * (ps0[kinv] - p0[kinv])
        xbot = tr0[kinv - 2 + m * mkx] + sstr0[kinv - 2 + m * mkx] * (ps0[kinv - 1] - p0[kinv - 2])
        k = 0
        while k <= kinv - 1:
            trflx[k + m * iface_stride] = _uwshcu_fluxbelowinv_value(
                cbmf, g_v, dt_v, ps0, kinv, trsrc[m], xmean, xtop, xbot, k
            )
            k += 1
        m += 1

    m = 0
    while m < wtrc_nwset:
        xmean = wt0[kinv - 1 + m * mkx]
        xtop = wt0[kinv + m * mkx] + sswt0[kinv + m * mkx] * (ps0[kinv] - p0[kinv])
        xbot = wt0[kinv - 2 + m * mkx] + sswt0[kinv - 2 + m * mkx] * (ps0[kinv - 1] - p0[kinv - 2])
        k = 0
        while k <= kinv - 1:
            wtflx[k + m * iface_stride] = _uwshcu_fluxbelowinv_value(
                cbmf, g_v, dt_v, ps0, kinv, wtsrc[m], xmean, xtop, xbot, k
            )
            k += 1
        m += 1

    uplus = 0.0
    vplus = 0.0
    k = kinv
    while k <= krel - 1:
        kp1_idx = k
        p_delta = ps0[k] - p0[kp1_idx]
        qtflx[k] = cbmf * (qtsrc - (qt0[kp1_idx] + ssqt0[kp1_idx] * p_delta))
        slflx[k] = cbmf * (thlsrc - (thl0[kp1_idx] + ssthl0[kp1_idx] * p_delta)) * cp * exns0[k]
        uplus = uplus + pgfc * ssu0[k - 1] * (ps0[k] - ps0[k - 1])
        vplus = vplus + pgfc * ssv0[k - 1] * (ps0[k] - ps0[k - 1])
        uflx[k] = cbmf * (usrc + uplus - (u0[kp1_idx] + ssu0[kp1_idx] * p_delta))
        vflx[k] = cbmf * (vsrc + vplus - (v0[kp1_idx] + ssv0[kp1_idx] * p_delta))

        m = 0
        while m < ncnst:
            trflx[k + m * iface_stride] = cbmf * (
                trsrc[m] - (tr0[kp1_idx + m * mkx] + sstr0[kp1_idx + m * mkx] * p_delta)
            )
            m += 1

        m = 0
        while m < wtrc_nwset:
            wtflx[k + m * iface_stride] = cbmf * (
                wtsrc[m] - (wt0[kp1_idx + m * mkx] + sswt0[kp1_idx + m * mkx] * p_delta)
            )
            m += 1
        k += 1

    k = krel
    while k <= kbup - 1:
        kp1_idx = k
        p_delta = ps0[k] - p0[kp1_idx]
        slflx[k] = cp * exns0[k] * umf[k] * (thlu[k] - (thl0[kp1_idx] + ssthl0[kp1_idx] * p_delta))
        qtflx[k] = umf[k] * (qtu[k] - (qt0[kp1_idx] + ssqt0[kp1_idx] * p_delta))
        uflx[k] = umf[k] * (uu[k] - (u0[kp1_idx] + ssu0[kp1_idx] * p_delta))
        vflx[k] = umf[k] * (vu[k] - (v0[kp1_idx] + ssv0[kp1_idx] * p_delta))

        m = 0
        while m < ncnst:
            trflx[k + m * iface_stride] = umf[k] * (
                tru[k + m * iface_stride] - (tr0[kp1_idx + m * mkx] + sstr0[kp1_idx + m * mkx] * p_delta)
            )
            m += 1

        m = 0
        while m < wtrc_nwset:
            wtflx[k + m * iface_stride] = umf[k] * (
                wtu[k + m * iface_stride] - (wt0[kp1_idx + m * mkx] + sswt0[kp1_idx + m * mkx] * p_delta)
            )
            m += 1
        k += 1

    k = kbup
    while k <= kpen - 1:
        layer_idx = k - 1
        p_delta = ps0[k] - p0[layer_idx]
        slflx[k] = cp * exns0[k] * emf[k] * (thlu_emf[k] - (thl0[layer_idx] + ssthl0[layer_idx] * p_delta))
        qtflx[k] = emf[k] * (qtu_emf[k] - (qt0[layer_idx] + ssqt0[layer_idx] * p_delta))
        uflx[k] = emf[k] * (uu_emf[k] - (u0[layer_idx] + ssu0[layer_idx] * p_delta))
        vflx[k] = emf[k] * (vu_emf[k] - (v0[layer_idx] + ssv0[layer_idx] * p_delta))

        m = 0
        while m < ncnst:
            trflx[k + m * iface_stride] = emf[k] * (
                tru_emf[k + m * iface_stride] - (tr0[layer_idx + m * mkx] + sstr0[layer_idx + m * mkx] * p_delta)
            )
            m += 1

        m = 0
        while m < wtrc_nwset:
            wtflx[k + m * iface_stride] = emf[k] * (
                wtu_emf[k + m * iface_stride] - (wt0[layer_idx + m * mkx] + sswt0[layer_idx + m * mkx] * p_delta)
            )
            m += 1
        k += 1

    if use_momenflx == 0:
        k = 0
        while k <= mkx:
            uflx[k] = 0.0
            vflx[k] = 0.0
            k += 1


@export
def uwshcu_massflux_comsub_shell_codon(
    mkx: int,
    kinv: int,
    krel: int,
    kbup: int,
    kpen: int,
    cbmf: float,
    ps0_p: cobj,
    umf_p: cobj,
    emf_p: cobj,
    uemf_p: cobj,
    comsub_p: cobj,
):
    ps0 = Ptr[float](ps0_p)
    umf = Ptr[float](umf_p)
    emf = Ptr[float](emf_p)
    uemf = Ptr[float](uemf_p)
    comsub = Ptr[float](comsub_p)

    k = 0
    while k <= mkx:
        uemf[k] = 0.0
        k += 1

    k = 0
    while k <= kinv - 2:
        uemf[k] = cbmf * (ps0[0] - ps0[k]) / (ps0[0] - ps0[kinv - 1])
        k += 1

    k = kinv - 1
    while k <= krel - 1:
        uemf[k] = cbmf
        k += 1

    k = krel
    while k <= kbup - 1:
        uemf[k] = umf[k]
        k += 1

    k = kbup
    while k <= kpen - 1:
        uemf[k] = emf[k]
        k += 1

    k = 0
    while k < mkx:
        comsub[k] = 0.0
        k += 1

    k = 1
    while k <= kpen:
        comsub[k - 1] = 0.5 * (uemf[k] + uemf[k - 1])
        k += 1


@export
def uwshcu_momentum_detrainment_shell_codon(
    mkx: int,
    wtrc_nwset: int,
    kpen: int,
    g_v: float,
    dt_v: float,
    dp0_p: cobj,
    u0_p: cobj,
    v0_p: cobj,
    uflx_p: cobj,
    vflx_p: cobj,
    umf_p: cobj,
    uten_p: cobj,
    vten_p: cobj,
    uf_p: cobj,
    vf_p: cobj,
    dwten_p: cobj,
    diten_p: cobj,
    wtdwten_p: cobj,
    wtditen_p: cobj,
):
    dp0 = Ptr[float](dp0_p)
    u0 = Ptr[float](u0_p)
    v0 = Ptr[float](v0_p)
    uflx = Ptr[float](uflx_p)
    vflx = Ptr[float](vflx_p)
    uten = Ptr[float](uten_p)
    vten = Ptr[float](vten_p)
    uf = Ptr[float](uf_p)
    vf = Ptr[float](vf_p)

    k = 1
    while k <= kpen:
        km1 = k - 1
        layer_idx = k - 1
        uten[layer_idx] = (uflx[km1] - uflx[k]) * g_v / dp0[layer_idx]
        vten[layer_idx] = (vflx[km1] - vflx[k]) * g_v / dp0[layer_idx]
        uf[layer_idx] = u0[layer_idx] + uten[layer_idx] * dt_v
        vf[layer_idx] = v0[layer_idx] + vten[layer_idx] * dt_v
        k += 1


@export
def uwshcu_comp_sub_tendency_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    kpen: int,
    ixnumliq: int,
    ixnumice: int,
    g_v: float,
    p0_p: cobj,
    thl0_p: cobj,
    qt0_p: cobj,
    ql0_p: cobj,
    qi0_p: cobj,
    tr0_p: cobj,
    wtrc_iatype_p: cobj,
    comsub_p: cobj,
    thlten_sub_p: cobj,
    qtten_sub_p: cobj,
    qlten_sub_p: cobj,
    qiten_sub_p: cobj,
    nlten_sub_p: cobj,
    niten_sub_p: cobj,
    wtlten_sub_p: cobj,
    wtiten_sub_p: cobj,
):
    p0 = Ptr[float](p0_p)
    thl0 = Ptr[float](thl0_p)
    qt0 = Ptr[float](qt0_p)
    ql0 = Ptr[float](ql0_p)
    qi0 = Ptr[float](qi0_p)
    tr0 = Ptr[float](tr0_p)
    wtrc_iatype = Ptr[int](wtrc_iatype_p)
    comsub = Ptr[float](comsub_p)
    thlten_sub = Ptr[float](thlten_sub_p)
    qtten_sub = Ptr[float](qtten_sub_p)
    qlten_sub = Ptr[float](qlten_sub_p)
    qiten_sub = Ptr[float](qiten_sub_p)
    nlten_sub = Ptr[float](nlten_sub_p)
    niten_sub = Ptr[float](niten_sub_p)
    wtlten_sub = Ptr[float](wtlten_sub_p)
    wtiten_sub = Ptr[float](wtiten_sub_p)

    liq_idx = ixnumliq - 1
    ice_idx = ixnumice - 1
    k_fortran = 1
    while k_fortran <= kpen:
        k = k_fortran - 1
        if comsub[k] >= 0.0:
            if k_fortran == mkx:
                thlten_sub[k] = 0.0
                qtten_sub[k] = 0.0
                qlten_sub[k] = 0.0
                qiten_sub[k] = 0.0
                nlten_sub[k] = 0.0
                niten_sub[k] = 0.0
                m = 0
                while m < wtrc_nwset:
                    wt_idx = k + m * mkx
                    wtlten_sub[wt_idx] = 0.0
                    wtiten_sub[wt_idx] = 0.0
                    m += 1
            else:
                kp1 = k + 1
                denom = p0[k] - p0[kp1]
                thlten_sub[k] = g_v * comsub[k] * (thl0[kp1] - thl0[k]) / denom
                qtten_sub[k] = g_v * comsub[k] * (qt0[kp1] - qt0[k]) / denom
                qlten_sub[k] = g_v * comsub[k] * (ql0[kp1] - ql0[k]) / denom
                qiten_sub[k] = g_v * comsub[k] * (qi0[kp1] - qi0[k]) / denom
                nlten_sub[k] = g_v * comsub[k] * (tr0[kp1 + liq_idx * mkx] - tr0[k + liq_idx * mkx]) / denom
                niten_sub[k] = g_v * comsub[k] * (tr0[kp1 + ice_idx * mkx] - tr0[k + ice_idx * mkx]) / denom
                m = 0
                while m < wtrc_nwset:
                    liq = wtrc_iatype[m + wtrc_nwset] - 1
                    ice = wtrc_iatype[m + 2 * wtrc_nwset] - 1
                    wt_idx = k + m * mkx
                    wtlten_sub[wt_idx] = g_v * comsub[k] * (tr0[kp1 + liq * mkx] - tr0[k + liq * mkx]) / denom
                    wtiten_sub[wt_idx] = g_v * comsub[k] * (tr0[kp1 + ice * mkx] - tr0[k + ice * mkx]) / denom
                    m += 1
        else:
            if k_fortran == 1:
                thlten_sub[k] = 0.0
                qtten_sub[k] = 0.0
                qlten_sub[k] = 0.0
                qiten_sub[k] = 0.0
                nlten_sub[k] = 0.0
                niten_sub[k] = 0.0
                m = 0
                while m < wtrc_nwset:
                    wt_idx = k + m * mkx
                    wtlten_sub[wt_idx] = 0.0
                    wtiten_sub[wt_idx] = 0.0
                    m += 1
            else:
                km1 = k - 1
                denom = p0[km1] - p0[k]
                thlten_sub[k] = g_v * comsub[k] * (thl0[k] - thl0[km1]) / denom
                qtten_sub[k] = g_v * comsub[k] * (qt0[k] - qt0[km1]) / denom
                qlten_sub[k] = g_v * comsub[k] * (ql0[k] - ql0[km1]) / denom
                qiten_sub[k] = g_v * comsub[k] * (qi0[k] - qi0[km1]) / denom
                nlten_sub[k] = g_v * comsub[k] * (tr0[k + liq_idx * mkx] - tr0[km1 + liq_idx * mkx]) / denom
                niten_sub[k] = g_v * comsub[k] * (tr0[k + ice_idx * mkx] - tr0[km1 + ice_idx * mkx]) / denom
                m = 0
                while m < wtrc_nwset:
                    liq = wtrc_iatype[m + wtrc_nwset] - 1
                    ice = wtrc_iatype[m + 2 * wtrc_nwset] - 1
                    wt_idx = k + m * mkx
                    wtlten_sub[wt_idx] = g_v * comsub[k] * (tr0[k + liq * mkx] - tr0[km1 + liq * mkx]) / denom
                    wtiten_sub[wt_idx] = g_v * comsub[k] * (tr0[k + ice * mkx] - tr0[km1 + ice * mkx]) / denom
                    m += 1
        k_fortran += 1


@export
def uwshcu_comp_sub_prepare_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    kinv: int,
    krel: int,
    kbup: int,
    kpen: int,
    ixnumliq: int,
    ixnumice: int,
    cbmf: float,
    g_v: float,
    ps0_p: cobj,
    umf_p: cobj,
    emf_p: cobj,
    uemf_p: cobj,
    comsub_p: cobj,
    p0_p: cobj,
    thl0_p: cobj,
    qt0_p: cobj,
    ql0_p: cobj,
    qi0_p: cobj,
    tr0_p: cobj,
    wtrc_iatype_p: cobj,
    thlten_sub_p: cobj,
    qtten_sub_p: cobj,
    qlten_sub_p: cobj,
    qiten_sub_p: cobj,
    nlten_sub_p: cobj,
    niten_sub_p: cobj,
    wtlten_sub_p: cobj,
    wtiten_sub_p: cobj,
):
    uwshcu_massflux_comsub_shell_codon(
        mkx,
        kinv,
        krel,
        kbup,
        kpen,
        cbmf,
        ps0_p,
        umf_p,
        emf_p,
        uemf_p,
        comsub_p,
    )
    uwshcu_comp_sub_tendency_shell_codon(
        mkx,
        ncnst,
        wtrc_nwset,
        kpen,
        ixnumliq,
        ixnumice,
        g_v,
        p0_p,
        thl0_p,
        qt0_p,
        ql0_p,
        qi0_p,
        tr0_p,
        wtrc_iatype_p,
        comsub_p,
        thlten_sub_p,
        qtten_sub_p,
        qlten_sub_p,
        qiten_sub_p,
        nlten_sub_p,
        niten_sub_p,
        wtlten_sub_p,
        wtiten_sub_p,
    )


@export
def uwshcu_penent_flux_comp_sub_prep_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    kinv: int,
    krel: int,
    kbup: int,
    kpen: int,
    ixnumliq: int,
    ixnumice: int,
    use_momenflx: int,
    r: float,
    g_v: float,
    dt_v: float,
    rpen: float,
    ppen: float,
    cbmf: float,
    cp: float,
    pgfc: float,
    qtsrc: float,
    thlsrc: float,
    usrc: float,
    vsrc: float,
    ps0_p: cobj,
    p0_p: cobj,
    dp0_p: cobj,
    thv0bot_p: cobj,
    thv0top_p: cobj,
    exns0_p: cobj,
    thl0_p: cobj,
    ssthl0_p: cobj,
    qt0_p: cobj,
    ssqt0_p: cobj,
    u0_p: cobj,
    ssu0_p: cobj,
    v0_p: cobj,
    ssv0_p: cobj,
    tr0_p: cobj,
    sstr0_p: cobj,
    wt0_p: cobj,
    sswt0_p: cobj,
    umf_p: cobj,
    emf_p: cobj,
    ufrc_p: cobj,
    dwten_p: cobj,
    diten_p: cobj,
    fer_p: cobj,
    fdr_p: cobj,
    rei_p: cobj,
    thlu_p: cobj,
    qtu_p: cobj,
    uu_p: cobj,
    vu_p: cobj,
    tru_p: cobj,
    wtu_p: cobj,
    thlu_emf_p: cobj,
    qtu_emf_p: cobj,
    uu_emf_p: cobj,
    vu_emf_p: cobj,
    tru_emf_p: cobj,
    wtu_emf_p: cobj,
    wtdwten_p: cobj,
    wtditen_p: cobj,
    limit_emf_p: cobj,
    trsrc_p: cobj,
    wtsrc_p: cobj,
    slflx_p: cobj,
    qtflx_p: cobj,
    uflx_p: cobj,
    vflx_p: cobj,
    trflx_p: cobj,
    wtflx_p: cobj,
    uemf_p: cobj,
    comsub_p: cobj,
    ql0_p: cobj,
    qi0_p: cobj,
    wtrc_iatype_p: cobj,
    thlten_sub_p: cobj,
    qtten_sub_p: cobj,
    qlten_sub_p: cobj,
    qiten_sub_p: cobj,
    nlten_sub_p: cobj,
    niten_sub_p: cobj,
    wtlten_sub_p: cobj,
    wtiten_sub_p: cobj,
):
    uwshcu_penent_prep_shell_codon(
        mkx,
        ncnst,
        wtrc_nwset,
        kbup,
        kpen,
        r,
        g_v,
        dt_v,
        rpen,
        ppen,
        ps0_p,
        p0_p,
        dp0_p,
        thv0bot_p,
        thv0top_p,
        exns0_p,
        thl0_p,
        ssthl0_p,
        qt0_p,
        ssqt0_p,
        u0_p,
        ssu0_p,
        v0_p,
        ssv0_p,
        tr0_p,
        sstr0_p,
        wt0_p,
        sswt0_p,
        umf_p,
        emf_p,
        ufrc_p,
        dwten_p,
        diten_p,
        fer_p,
        fdr_p,
        rei_p,
        thlu_p,
        qtu_p,
        uu_p,
        vu_p,
        tru_p,
        wtu_p,
        thlu_emf_p,
        qtu_emf_p,
        uu_emf_p,
        vu_emf_p,
        tru_emf_p,
        wtu_emf_p,
        wtdwten_p,
        wtditen_p,
        limit_emf_p,
    )
    uwshcu_turbulent_flux_shell_codon(
        mkx,
        ncnst,
        wtrc_nwset,
        kinv,
        krel,
        kbup,
        kpen,
        use_momenflx,
        cbmf,
        g_v,
        dt_v,
        cp,
        pgfc,
        qtsrc,
        thlsrc,
        usrc,
        vsrc,
        ps0_p,
        p0_p,
        exns0_p,
        qt0_p,
        ssqt0_p,
        thl0_p,
        ssthl0_p,
        u0_p,
        ssu0_p,
        v0_p,
        ssv0_p,
        tr0_p,
        sstr0_p,
        wt0_p,
        sswt0_p,
        trsrc_p,
        wtsrc_p,
        umf_p,
        emf_p,
        thlu_p,
        qtu_p,
        uu_p,
        vu_p,
        tru_p,
        wtu_p,
        thlu_emf_p,
        qtu_emf_p,
        uu_emf_p,
        vu_emf_p,
        tru_emf_p,
        wtu_emf_p,
        slflx_p,
        qtflx_p,
        uflx_p,
        vflx_p,
        trflx_p,
        wtflx_p,
    )
    uwshcu_comp_sub_prepare_shell_codon(
        mkx,
        ncnst,
        wtrc_nwset,
        kinv,
        krel,
        kbup,
        kpen,
        ixnumliq,
        ixnumice,
        cbmf,
        g_v,
        ps0_p,
        umf_p,
        emf_p,
        uemf_p,
        comsub_p,
        p0_p,
        thl0_p,
        qt0_p,
        ql0_p,
        qi0_p,
        tr0_p,
        wtrc_iatype_p,
        thlten_sub_p,
        qtten_sub_p,
        qlten_sub_p,
        qiten_sub_p,
        nlten_sub_p,
        niten_sub_p,
        wtlten_sub_p,
        wtiten_sub_p,
    )


@export
def uwshcu_scaleh_filter_penent_flux_comp_sub_prep_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    kinv: int,
    krel: int,
    kbup: int,
    kpen: int,
    ixnumliq: int,
    ixnumice: int,
    use_momenflx: int,
    r: float,
    g_v: float,
    dt_v: float,
    rpen: float,
    ppen: float,
    cbmf: float,
    cp: float,
    pgfc: float,
    qtsrc: float,
    thlsrc: float,
    usrc: float,
    vsrc: float,
    ps0_p: cobj,
    p0_p: cobj,
    dp0_p: cobj,
    thv0bot_p: cobj,
    thv0top_p: cobj,
    exns0_p: cobj,
    thl0_p: cobj,
    ssthl0_p: cobj,
    qt0_p: cobj,
    ssqt0_p: cobj,
    u0_p: cobj,
    ssu0_p: cobj,
    v0_p: cobj,
    ssv0_p: cobj,
    tr0_p: cobj,
    sstr0_p: cobj,
    wt0_p: cobj,
    sswt0_p: cobj,
    umf_p: cobj,
    emf_p: cobj,
    ufrc_p: cobj,
    dwten_p: cobj,
    diten_p: cobj,
    fer_p: cobj,
    fdr_p: cobj,
    rei_p: cobj,
    thlu_p: cobj,
    qtu_p: cobj,
    uu_p: cobj,
    vu_p: cobj,
    tru_p: cobj,
    wtu_p: cobj,
    thlu_emf_p: cobj,
    qtu_emf_p: cobj,
    uu_emf_p: cobj,
    vu_emf_p: cobj,
    tru_emf_p: cobj,
    wtu_emf_p: cobj,
    wtdwten_p: cobj,
    wtditen_p: cobj,
    limit_emf_p: cobj,
    limit_shcu_p: cobj,
    exit_code_p: cobj,
    trsrc_p: cobj,
    wtsrc_p: cobj,
    slflx_p: cobj,
    qtflx_p: cobj,
    uflx_p: cobj,
    vflx_p: cobj,
    trflx_p: cobj,
    wtflx_p: cobj,
    uemf_p: cobj,
    comsub_p: cobj,
    ql0_p: cobj,
    qi0_p: cobj,
    wtrc_iatype_p: cobj,
    thlten_sub_p: cobj,
    qtten_sub_p: cobj,
    qlten_sub_p: cobj,
    qiten_sub_p: cobj,
    nlten_sub_p: cobj,
    niten_sub_p: cobj,
    wtlten_sub_p: cobj,
    wtiten_sub_p: cobj,
):
    limit_shcu = Ptr[float](limit_shcu_p)
    exit_code = Ptr[int](exit_code_p)

    exit_code[0] = 0
    if kbup == krel:
        limit_shcu[0] = 1.0
        exit_code[0] = 1
        return
    limit_shcu[0] = 0.0

    uwshcu_penent_flux_comp_sub_prep_shell_codon(
        mkx,
        ncnst,
        wtrc_nwset,
        kinv,
        krel,
        kbup,
        kpen,
        ixnumliq,
        ixnumice,
        use_momenflx,
        r,
        g_v,
        dt_v,
        rpen,
        ppen,
        cbmf,
        cp,
        pgfc,
        qtsrc,
        thlsrc,
        usrc,
        vsrc,
        ps0_p,
        p0_p,
        dp0_p,
        thv0bot_p,
        thv0top_p,
        exns0_p,
        thl0_p,
        ssthl0_p,
        qt0_p,
        ssqt0_p,
        u0_p,
        ssu0_p,
        v0_p,
        ssv0_p,
        tr0_p,
        sstr0_p,
        wt0_p,
        sswt0_p,
        umf_p,
        emf_p,
        ufrc_p,
        dwten_p,
        diten_p,
        fer_p,
        fdr_p,
        rei_p,
        thlu_p,
        qtu_p,
        uu_p,
        vu_p,
        tru_p,
        wtu_p,
        thlu_emf_p,
        qtu_emf_p,
        uu_emf_p,
        vu_emf_p,
        tru_emf_p,
        wtu_emf_p,
        wtdwten_p,
        wtditen_p,
        limit_emf_p,
        trsrc_p,
        wtsrc_p,
        slflx_p,
        qtflx_p,
        uflx_p,
        vflx_p,
        trflx_p,
        wtflx_p,
        uemf_p,
        comsub_p,
        ql0_p,
        qi0_p,
        wtrc_iatype_p,
        thlten_sub_p,
        qtten_sub_p,
        qlten_sub_p,
        qiten_sub_p,
        nlten_sub_p,
        niten_sub_p,
        wtlten_sub_p,
        wtiten_sub_p,
    )


@export
def uwshcu_scaleh_cufilter_exit_shell_codon(
    post_exit_code: int,
    exit_cufilter_p: cobj,
    exit_code_p: cobj,
):
    exit_cufilter = Ptr[float](exit_cufilter_p)
    exit_code = Ptr[int](exit_code_p)

    exit_code[0] = 0
    if post_exit_code != 0:
        exit_cufilter[0] = 1.0
        exit_code[0] = 1


@export
def uwshcu_comp_sub_sink_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    kpen: int,
    ixnumliq: int,
    ixnumice: int,
    dt_v: float,
    ql0_p: cobj,
    qi0_p: cobj,
    tr0_p: cobj,
    wtrc_iatype_p: cobj,
    qlten_sub_p: cobj,
    qiten_sub_p: cobj,
    nlten_sub_p: cobj,
    niten_sub_p: cobj,
    wtlten_sub_p: cobj,
    wtiten_sub_p: cobj,
    qlten_sink_p: cobj,
    qiten_sink_p: cobj,
    nlten_sink_p: cobj,
    niten_sink_p: cobj,
    wtten_sink_liq_p: cobj,
    wtten_sink_ice_p: cobj,
):
    ql0 = Ptr[float](ql0_p)
    qi0 = Ptr[float](qi0_p)
    tr0 = Ptr[float](tr0_p)
    wtrc_iatype = Ptr[int](wtrc_iatype_p)
    qlten_sub = Ptr[float](qlten_sub_p)
    qiten_sub = Ptr[float](qiten_sub_p)
    nlten_sub = Ptr[float](nlten_sub_p)
    niten_sub = Ptr[float](niten_sub_p)
    wtlten_sub = Ptr[float](wtlten_sub_p)
    wtiten_sub = Ptr[float](wtiten_sub_p)
    qlten_sink = Ptr[float](qlten_sink_p)
    qiten_sink = Ptr[float](qiten_sink_p)
    nlten_sink = Ptr[float](nlten_sink_p)
    niten_sink = Ptr[float](niten_sink_p)
    wtten_sink_liq = Ptr[float](wtten_sink_liq_p)
    wtten_sink_ice = Ptr[float](wtten_sink_ice_p)

    liq_idx = ixnumliq - 1
    ice_idx = ixnumice - 1
    k = 0
    while k < kpen:
        ql_floor = -ql0[k] / dt_v
        qi_floor = -qi0[k] / dt_v
        nl_floor = -tr0[k + liq_idx * mkx] / dt_v
        ni_floor = -tr0[k + ice_idx * mkx] / dt_v

        if qlten_sub[k] > ql_floor:
            qlten_sink[k] = qlten_sub[k]
        else:
            qlten_sink[k] = ql_floor

        if qiten_sub[k] > qi_floor:
            qiten_sink[k] = qiten_sub[k]
        else:
            qiten_sink[k] = qi_floor

        if nlten_sub[k] > nl_floor:
            nlten_sink[k] = nlten_sub[k]
        else:
            nlten_sink[k] = nl_floor

        if niten_sub[k] > ni_floor:
            niten_sink[k] = niten_sub[k]
        else:
            niten_sink[k] = ni_floor

        m = 0
        while m < wtrc_nwset:
            liq = wtrc_iatype[m + wtrc_nwset] - 1
            ice = wtrc_iatype[m + 2 * wtrc_nwset] - 1
            wt_idx = k + m * mkx
            wt_liq_floor = -tr0[k + liq * mkx] / dt_v
            wt_ice_floor = -tr0[k + ice * mkx] / dt_v
            if wtlten_sub[wt_idx] > wt_liq_floor:
                wtten_sink_liq[wt_idx] = wtlten_sub[wt_idx]
            else:
                wtten_sink_liq[wt_idx] = wt_liq_floor
            if wtiten_sub[wt_idx] > wt_ice_floor:
                wtten_sink_ice[wt_idx] = wtiten_sub[wt_idx]
            else:
                wtten_sink_ice[wt_idx] = wt_ice_floor
            m += 1

        k += 1


@export
def uwshcu_comp_sub_conden_exit_shell_codon(
    id_check: int,
    exit_code_p: cobj,
):
    exit_code = Ptr[int](exit_code_p)

    exit_code[0] = 0
    if id_check == 1:
        exit_code[0] = 1


@export
def uwshcu_thermo_conden_exit_shell_codon(
    id_check: int,
    exit_conden_p: cobj,
    exit_code_p: cobj,
):
    exit_conden = Ptr[float](exit_conden_p)
    exit_code = Ptr[int](exit_code_p)

    exit_code[0] = 0
    if id_check == 1:
        exit_conden[0] = 1.0
        exit_code[0] = 1


@export
def uwshcu_thermo_emf_conden_exit_shell_codon(
    id_check: int,
    exit_code_p: cobj,
):
    exit_code = Ptr[int](exit_code_p)

    exit_code[0] = 0
    if id_check == 1:
        exit_code[0] = 1


@export
def uwshcu_thermo_emf_kbup_state_shell_codon(
    k: int,
    mkx: int,
    ixnumliq: int,
    ixnumice: int,
    ql_emf_kbup: float,
    qi_emf_kbup: float,
    tru_emf_p: cobj,
    nl_emf_kbup_p: cobj,
    ni_emf_kbup_p: cobj,
):
    tru_emf = Ptr[float](tru_emf_p)
    nl_emf_kbup = Ptr[float](nl_emf_kbup_p)
    ni_emf_kbup = Ptr[float](ni_emf_kbup_p)

    iface_stride = mkx + 1
    liq_idx = ixnumliq - 1
    ice_idx = ixnumice - 1

    nl_emf_kbup[0] = 0.0
    if ql_emf_kbup > 0.0:
        nl_emf_kbup[0] = tru_emf[k + liq_idx * iface_stride]

    ni_emf_kbup[0] = 0.0
    if qi_emf_kbup > 0.0:
        ni_emf_kbup[0] = tru_emf[k + ice_idx * iface_stride]


@export
def uwshcu_thermo_emf_kbup_tendency_shell_codon(
    g_v: float,
    emf_k: float,
    ql_emf_kbup: float,
    qi_emf_kbup: float,
    nl_emf_kbup_v: float,
    ni_emf_kbup_v: float,
    ql0_k: float,
    qi0_k: float,
    tr0_liq_k: float,
    tr0_ice_k: float,
    ps0_km1: float,
    ps0_k: float,
    qc_lm_p: cobj,
    qc_im_p: cobj,
    nc_lm_p: cobj,
    nc_im_p: cobj,
):
    qc_lm = Ptr[float](qc_lm_p)
    qc_im = Ptr[float](qc_im_p)
    nc_lm = Ptr[float](nc_lm_p)
    nc_im = Ptr[float](nc_im_p)

    denom = ps0_km1 - ps0_k
    qc_lm[0] = qc_lm[0] - g_v * emf_k * (ql_emf_kbup - ql0_k) / denom
    qc_im[0] = qc_im[0] - g_v * emf_k * (qi_emf_kbup - qi0_k) / denom
    nc_lm[0] = nc_lm[0] - g_v * emf_k * (nl_emf_kbup_v - tr0_liq_k) / denom
    nc_im[0] = nc_im[0] - g_v * emf_k * (ni_emf_kbup_v - tr0_ice_k) / denom


@export
def uwshcu_thermo_sustain_shell_codon(
    frc_rasn: float,
    dwten_k: float,
    diten_k: float,
    qc_l_k_p: cobj,
    qc_i_k_p: cobj,
):
    qc_l_k = Ptr[float](qc_l_k_p)
    qc_i_k = Ptr[float](qc_i_k_p)

    qc_l_k[0] = (1.0 - frc_rasn) * dwten_k
    qc_i_k[0] = (1.0 - frc_rasn) * diten_k


@export
def uwshcu_thermo_detrain_shell_codon(
    k_le_kbup: int,
    g_v: float,
    umf_km1: float,
    umf_k: float,
    fdr_k: float,
    qlu_mid: float,
    qiu_mid: float,
    ql0_k: float,
    qi0_k: float,
    tr0_liq_k: float,
    tr0_ice_k: float,
    qc_l_k_p: cobj,
    qc_i_k_p: cobj,
    qc_lm_p: cobj,
    qc_im_p: cobj,
    nc_lm_p: cobj,
    nc_im_p: cobj,
):
    qc_l_k = Ptr[float](qc_l_k_p)
    qc_i_k = Ptr[float](qc_i_k_p)
    qc_lm = Ptr[float](qc_lm_p)
    qc_im = Ptr[float](qc_im_p)
    nc_lm = Ptr[float](nc_lm_p)
    nc_im = Ptr[float](nc_im_p)

    if k_le_kbup != 0:
        qc_l_k[0] = qc_l_k[0] + g_v * 0.5 * (umf_km1 + umf_k) * fdr_k * qlu_mid
        qc_i_k[0] = qc_i_k[0] + g_v * 0.5 * (umf_km1 + umf_k) * fdr_k * qiu_mid
        qc_lm[0] = -g_v * 0.5 * (umf_km1 + umf_k) * fdr_k * ql0_k
        qc_im[0] = -g_v * 0.5 * (umf_km1 + umf_k) * fdr_k * qi0_k
        nc_lm[0] = -g_v * 0.5 * (umf_km1 + umf_k) * fdr_k * tr0_liq_k
        nc_im[0] = -g_v * 0.5 * (umf_km1 + umf_k) * fdr_k * tr0_ice_k
    else:
        qc_lm[0] = 0.0
        qc_im[0] = 0.0
        nc_lm[0] = 0.0
        nc_im[0] = 0.0


@export
def uwshcu_thermo_detached_shell_codon(
    g_v: float,
    umf_k: float,
    qlj: float,
    qij: float,
    ql0_k: float,
    qi0_k: float,
    tr0_liq_k: float,
    tr0_ice_k: float,
    ps0_km1: float,
    ps0_k: float,
    qc_l_k_p: cobj,
    qc_i_k_p: cobj,
    qc_lm_p: cobj,
    qc_im_p: cobj,
    nc_lm_p: cobj,
    nc_im_p: cobj,
):
    qc_l_k = Ptr[float](qc_l_k_p)
    qc_i_k = Ptr[float](qc_i_k_p)
    qc_lm = Ptr[float](qc_lm_p)
    qc_im = Ptr[float](qc_im_p)
    nc_lm = Ptr[float](nc_lm_p)
    nc_im = Ptr[float](nc_im_p)

    qc_l_k[0] = qc_l_k[0] + g_v * umf_k * qlj / (ps0_km1 - ps0_k)
    qc_i_k[0] = qc_i_k[0] + g_v * umf_k * qij / (ps0_km1 - ps0_k)
    qc_lm[0] = qc_lm[0] - g_v * umf_k * ql0_k / (ps0_km1 - ps0_k)
    qc_im[0] = qc_im[0] - g_v * umf_k * qi0_k / (ps0_km1 - ps0_k)
    nc_lm[0] = nc_lm[0] - g_v * umf_k * tr0_liq_k / (ps0_km1 - ps0_k)
    nc_im[0] = nc_im[0] - g_v * umf_k * tr0_ice_k / (ps0_km1 - ps0_k)


@export
def uwshcu_thermo_condensate_batch_shell_codon(
    kind: int,
    k_fortran: int,
    mkx: int,
    ixnumliq: int,
    ixnumice: int,
    flag1: int,
    flag2: int,
    frc_rasn: float,
    g_v: float,
    dwten_k: float,
    diten_k: float,
    umf_km1: float,
    umf_k: float,
    fdr_k: float,
    qlu_mid: float,
    qiu_mid: float,
    qlj: float,
    qij: float,
    ql0_k: float,
    qi0_k: float,
    tr0_liq_k: float,
    tr0_ice_k: float,
    ps0_km1: float,
    ps0_k: float,
    emf_k: float,
    ql_emf_kbup: float,
    qi_emf_kbup: float,
    tru_emf_p: cobj,
    qc_l_k_p: cobj,
    qc_i_k_p: cobj,
    qc_lm_p: cobj,
    qc_im_p: cobj,
    nc_lm_p: cobj,
    nc_im_p: cobj,
    nl_emf_kbup_p: cobj,
    ni_emf_kbup_p: cobj,
):
    qc_l_k = Ptr[float](qc_l_k_p)
    qc_i_k = Ptr[float](qc_i_k_p)
    qc_lm = Ptr[float](qc_lm_p)
    qc_im = Ptr[float](qc_im_p)
    nc_lm = Ptr[float](nc_lm_p)
    nc_im = Ptr[float](nc_im_p)

    if kind == 1:
        qc_l_k[0] = (1.0 - frc_rasn) * dwten_k
        qc_i_k[0] = (1.0 - frc_rasn) * diten_k

        if flag1 != 0:
            qc_l_k[0] = qc_l_k[0] + g_v * 0.5 * (umf_km1 + umf_k) * fdr_k * qlu_mid
            qc_i_k[0] = qc_i_k[0] + g_v * 0.5 * (umf_km1 + umf_k) * fdr_k * qiu_mid
            qc_lm[0] = -g_v * 0.5 * (umf_km1 + umf_k) * fdr_k * ql0_k
            qc_im[0] = -g_v * 0.5 * (umf_km1 + umf_k) * fdr_k * qi0_k
            nc_lm[0] = -g_v * 0.5 * (umf_km1 + umf_k) * fdr_k * tr0_liq_k
            nc_im[0] = -g_v * 0.5 * (umf_km1 + umf_k) * fdr_k * tr0_ice_k
        else:
            qc_lm[0] = 0.0
            qc_im[0] = 0.0
            nc_lm[0] = 0.0
            nc_im[0] = 0.0

        if flag2 != 0:
            denom = ps0_km1 - ps0_k
            qc_l_k[0] = qc_l_k[0] + g_v * umf_k * qlj / denom
            qc_i_k[0] = qc_i_k[0] + g_v * umf_k * qij / denom
            qc_lm[0] = qc_lm[0] - g_v * umf_k * ql0_k / denom
            qc_im[0] = qc_im[0] - g_v * umf_k * qi0_k / denom
            nc_lm[0] = nc_lm[0] - g_v * umf_k * tr0_liq_k / denom
            nc_im[0] = nc_im[0] - g_v * umf_k * tr0_ice_k / denom
    elif kind == 2:
        tru_emf = Ptr[float](tru_emf_p)
        nl_emf_kbup = Ptr[float](nl_emf_kbup_p)
        ni_emf_kbup = Ptr[float](ni_emf_kbup_p)

        iface_stride = mkx + 1
        liq_idx = ixnumliq - 1
        ice_idx = ixnumice - 1

        nl_emf_kbup[0] = 0.0
        if ql_emf_kbup > 0.0:
            nl_emf_kbup[0] = tru_emf[k_fortran + liq_idx * iface_stride]

        ni_emf_kbup[0] = 0.0
        if qi_emf_kbup > 0.0:
            ni_emf_kbup[0] = tru_emf[k_fortran + ice_idx * iface_stride]

        denom = ps0_km1 - ps0_k
        qc_lm[0] = qc_lm[0] - g_v * emf_k * (ql_emf_kbup - ql0_k) / denom
        qc_im[0] = qc_im[0] - g_v * emf_k * (qi_emf_kbup - qi0_k) / denom
        nc_lm[0] = nc_lm[0] - g_v * emf_k * (nl_emf_kbup[0] - tr0_liq_k) / denom
        nc_im[0] = nc_im[0] - g_v * emf_k * (ni_emf_kbup[0] - tr0_ice_k) / denom


@export
def uwshcu_thermo_conden_condensate_batch_shell_codon(
    kind: int,
    k_fortran: int,
    mkx: int,
    ixnumliq: int,
    ixnumice: int,
    id_check: int,
    flag1: int,
    flag2: int,
    frc_rasn: float,
    g_v: float,
    dwten_k: float,
    diten_k: float,
    umf_km1: float,
    umf_k: float,
    fdr_k: float,
    qlj: float,
    qij: float,
    ql0_k: float,
    qi0_k: float,
    tr0_liq_k: float,
    tr0_ice_k: float,
    ps0_km1: float,
    ps0_k: float,
    prel_v: float,
    ppen_v: float,
    emf_k: float,
    ql_emf_kbup: float,
    qi_emf_kbup: float,
    qlubelow_p: cobj,
    qiubelow_p: cobj,
    qlu_mid_p: cobj,
    qiu_mid_p: cobj,
    qlu_top_p: cobj,
    qiu_top_p: cobj,
    exit_conden_p: cobj,
    exit_code_p: cobj,
    tru_emf_p: cobj,
    qc_l_k_p: cobj,
    qc_i_k_p: cobj,
    qc_lm_p: cobj,
    qc_im_p: cobj,
    nc_lm_p: cobj,
    nc_im_p: cobj,
    nl_emf_kbup_p: cobj,
    ni_emf_kbup_p: cobj,
):
    exit_code = Ptr[int](exit_code_p)

    exit_code[0] = 0
    if id_check == 1:
        if kind != 5:
            exit_conden = Ptr[float](exit_conden_p)
            exit_conden[0] = 1.0
        exit_code[0] = 1
        return

    if kind == 5:
        uwshcu_thermo_condensate_batch_shell_codon(
            2,
            k_fortran,
            mkx,
            ixnumliq,
            ixnumice,
            0,
            0,
            0.0,
            g_v,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            ql0_k,
            qi0_k,
            tr0_liq_k,
            tr0_ice_k,
            ps0_km1,
            ps0_k,
            emf_k,
            ql_emf_kbup,
            qi_emf_kbup,
            tru_emf_p,
            qc_l_k_p,
            qc_i_k_p,
            qc_lm_p,
            qc_im_p,
            nc_lm_p,
            nc_im_p,
            nl_emf_kbup_p,
            ni_emf_kbup_p,
        )
        return

    qlubelow = Ptr[float](qlubelow_p)
    qiubelow = Ptr[float](qiubelow_p)
    qlu_mid = Ptr[float](qlu_mid_p)
    qiu_mid = Ptr[float](qiu_mid_p)

    if kind == 0:
        qlu_mid[0] = 0.0
        qiu_mid[0] = 0.0
    elif kind == 1:
        qlubelow[0] = qlj
        qiubelow[0] = qij
        return
    elif kind == 2:
        qlu_mid[0] = (0.5 * (qlubelow[0] + qlj) * (prel_v - ps0_k)) / (ps0_km1 - ps0_k)
        qiu_mid[0] = (0.5 * (qiubelow[0] + qij) * (prel_v - ps0_k)) / (ps0_km1 - ps0_k)
    elif kind == 3:
        qlu_top = Ptr[float](qlu_top_p)
        qiu_top = Ptr[float](qiu_top_p)
        qlu_mid[0] = (0.5 * (qlubelow[0] + qlj) * (-ppen_v)) / (ps0_km1 - ps0_k)
        qiu_mid[0] = (0.5 * (qiubelow[0] + qij) * (-ppen_v)) / (ps0_km1 - ps0_k)
        qlu_top[0] = qlj
        qiu_top[0] = qij
    elif kind == 4:
        qlu_mid[0] = 0.5 * (qlubelow[0] + qlj)
        qiu_mid[0] = 0.5 * (qiubelow[0] + qij)

    qlubelow[0] = qlj
    qiubelow[0] = qij

    uwshcu_thermo_condensate_batch_shell_codon(
        1,
        k_fortran,
        mkx,
        ixnumliq,
        ixnumice,
        flag1,
        flag2,
        frc_rasn,
        g_v,
        dwten_k,
        diten_k,
        umf_km1,
        umf_k,
        fdr_k,
        qlu_mid[0],
        qiu_mid[0],
        qlj,
        qij,
        ql0_k,
        qi0_k,
        tr0_liq_k,
        tr0_ice_k,
        ps0_km1,
        ps0_k,
        0.0,
        0.0,
        0.0,
        tru_emf_p,
        qc_l_k_p,
        qc_i_k_p,
        qc_lm_p,
        qc_im_p,
        nc_lm_p,
        nc_im_p,
        nl_emf_kbup_p,
        ni_emf_kbup_p,
    )


@export
def uwshcu_thermo_prelim_shell_codon(
    mkx: int,
    wtrc_nwset: int,
    kpen: int,
    frc_rasn: float,
    g_v: float,
    dp0_p: cobj,
    umf_p: cobj,
    dwten_p: cobj,
    diten_p: cobj,
    wtdwten_p: cobj,
    wtditen_p: cobj,
    qrten_p: cobj,
    qsten_p: cobj,
    wtrpten_p: cobj,
    wtspten_p: cobj,
    slflx_p: cobj,
    qtflx_p: cobj,
    uflx_p: cobj,
    vflx_p: cobj,
    uf_p: cobj,
    vf_p: cobj,
    u0_p: cobj,
    v0_p: cobj,
    wtflx_p: cobj,
    slten_p: cobj,
    qtten_p: cobj,
    wttotten_p: cobj,
    rainflx_p: cobj,
    snowflx_p: cobj,
):
    dp0 = Ptr[float](dp0_p)
    umf = Ptr[float](umf_p)
    dwten = Ptr[float](dwten_p)
    diten = Ptr[float](diten_p)
    wtdwten = Ptr[float](wtdwten_p)
    wtditen = Ptr[float](wtditen_p)
    qrten = Ptr[float](qrten_p)
    qsten = Ptr[float](qsten_p)
    wtrpten = Ptr[float](wtrpten_p)
    wtspten = Ptr[float](wtspten_p)
    slflx = Ptr[float](slflx_p)
    qtflx = Ptr[float](qtflx_p)
    uflx = Ptr[float](uflx_p)
    vflx = Ptr[float](vflx_p)
    uf = Ptr[float](uf_p)
    vf = Ptr[float](vf_p)
    u0 = Ptr[float](u0_p)
    v0 = Ptr[float](v0_p)
    wtflx = Ptr[float](wtflx_p)
    slten = Ptr[float](slten_p)
    qtten = Ptr[float](qtten_p)
    wttotten = Ptr[float](wttotten_p)
    rainflx = Ptr[float](rainflx_p)
    snowflx = Ptr[float](snowflx_p)

    rainflx[0] = 0.0
    snowflx[0] = 0.0

    k_fortran = 1
    while k_fortran <= kpen:
        k = k_fortran - 1
        km1 = k_fortran - 1

        dwten[k] = dwten[k] * 0.5 * (umf[k_fortran - 1] + umf[k_fortran]) * g_v / dp0[k]
        diten[k] = diten[k] * 0.5 * (umf[k_fortran - 1] + umf[k_fortran]) * g_v / dp0[k]

        m = 0
        while m < wtrc_nwset:
            wt_idx = k + m * mkx
            wtdwten[wt_idx] = wtdwten[wt_idx] * 0.5 * (umf[k_fortran - 1] + umf[k_fortran]) * g_v / dp0[k]
            wtditen[wt_idx] = wtditen[wt_idx] * 0.5 * (umf[k_fortran - 1] + umf[k_fortran]) * g_v / dp0[k]
            m += 1

        qrten[k] = frc_rasn * dwten[k]
        qsten[k] = frc_rasn * diten[k]

        m = 0
        while m < wtrc_nwset:
            wt_idx = k + m * mkx
            wtrpten[wt_idx] = frc_rasn * wtdwten[wt_idx]
            wtspten[wt_idx] = frc_rasn * wtditen[wt_idx]
            m += 1

        rainflx[0] = rainflx[0] + qrten[k] * dp0[k] / g_v
        snowflx[0] = snowflx[0] + qsten[k] * dp0[k] / g_v

        slten[k] = (slflx[km1] - slflx[k_fortran]) * g_v / dp0[k]
        if k_fortran == 1:
            slten[k] = slten[k] - g_v / 4.0 / dp0[k] * (
                uflx[k_fortran] * (uf[k + 1] - uf[k] + u0[k + 1] - u0[k])
                + vflx[k_fortran] * (vf[k + 1] - vf[k] + v0[k + 1] - v0[k])
            )
        elif k_fortran >= 2 and k_fortran <= kpen - 1:
            slten[k] = slten[k] - g_v / 4.0 / dp0[k] * (
                uflx[k_fortran] * (uf[k + 1] - uf[k] + u0[k + 1] - u0[k])
                + uflx[k_fortran - 1] * (uf[k] - uf[k - 1] + u0[k] - u0[k - 1])
                + vflx[k_fortran] * (vf[k + 1] - vf[k] + v0[k + 1] - v0[k])
                + vflx[k_fortran - 1] * (vf[k] - vf[k - 1] + v0[k] - v0[k - 1])
            )
        elif k_fortran == kpen:
            slten[k] = slten[k] - g_v / 4.0 / dp0[k] * (
                uflx[k_fortran - 1] * (uf[k] - uf[k - 1] + u0[k] - u0[k - 1])
                + vflx[k_fortran - 1] * (vf[k] - vf[k - 1] + v0[k] - v0[k - 1])
            )

        qtten[k] = (qtflx[km1] - qtflx[k_fortran]) * g_v / dp0[k]

        m = 0
        while m < wtrc_nwset:
            wt_idx = k + m * mkx
            wtflx_offset = m * (mkx + 1)
            wttotten[wt_idx] = (wtflx[km1 + wtflx_offset] - wtflx[k_fortran + wtflx_offset]) * g_v / dp0[k]
            m += 1

        k_fortran += 1


@export
def uwshcu_tendency_prep_shell_codon(
    mkx: int,
    wtrc_nwset: int,
    kpen: int,
    frc_rasn: float,
    g_v: float,
    dt_v: float,
    dp0_p: cobj,
    u0_p: cobj,
    v0_p: cobj,
    uflx_p: cobj,
    vflx_p: cobj,
    umf_p: cobj,
    uten_p: cobj,
    vten_p: cobj,
    uf_p: cobj,
    vf_p: cobj,
    dwten_p: cobj,
    diten_p: cobj,
    wtdwten_p: cobj,
    wtditen_p: cobj,
    qrten_p: cobj,
    qsten_p: cobj,
    wtrpten_p: cobj,
    wtspten_p: cobj,
    slflx_p: cobj,
    qtflx_p: cobj,
    wtflx_p: cobj,
    slten_p: cobj,
    qtten_p: cobj,
    wttotten_p: cobj,
    rainflx_p: cobj,
    snowflx_p: cobj,
):
    uwshcu_momentum_detrainment_shell_codon(
        mkx,
        wtrc_nwset,
        kpen,
        g_v,
        dt_v,
        dp0_p,
        u0_p,
        v0_p,
        uflx_p,
        vflx_p,
        umf_p,
        uten_p,
        vten_p,
        uf_p,
        vf_p,
        dwten_p,
        diten_p,
        wtdwten_p,
        wtditen_p,
    )
    uwshcu_thermo_prelim_shell_codon(
        mkx,
        wtrc_nwset,
        kpen,
        frc_rasn,
        g_v,
        dp0_p,
        umf_p,
        dwten_p,
        diten_p,
        wtdwten_p,
        wtditen_p,
        qrten_p,
        qsten_p,
        wtrpten_p,
        wtspten_p,
        slflx_p,
        qtflx_p,
        uflx_p,
        vflx_p,
        uf_p,
        vf_p,
        u0_p,
        v0_p,
        wtflx_p,
        slten_p,
        qtten_p,
        wttotten_p,
        rainflx_p,
        snowflx_p,
    )


@export
def uwshcu_comp_sub_sink_thermo_prelim_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    kpen: int,
    ixnumliq: int,
    ixnumice: int,
    frc_rasn: float,
    g_v: float,
    dt_v: float,
    ql0_p: cobj,
    qi0_p: cobj,
    tr0_p: cobj,
    wtrc_iatype_p: cobj,
    qlten_sub_p: cobj,
    qiten_sub_p: cobj,
    nlten_sub_p: cobj,
    niten_sub_p: cobj,
    wtlten_sub_p: cobj,
    wtiten_sub_p: cobj,
    qlten_sink_p: cobj,
    qiten_sink_p: cobj,
    nlten_sink_p: cobj,
    niten_sink_p: cobj,
    wtten_sink_liq_p: cobj,
    wtten_sink_ice_p: cobj,
    dp0_p: cobj,
    u0_p: cobj,
    v0_p: cobj,
    uflx_p: cobj,
    vflx_p: cobj,
    umf_p: cobj,
    uten_p: cobj,
    vten_p: cobj,
    uf_p: cobj,
    vf_p: cobj,
    dwten_p: cobj,
    diten_p: cobj,
    wtdwten_p: cobj,
    wtditen_p: cobj,
    qrten_p: cobj,
    qsten_p: cobj,
    wtrpten_p: cobj,
    wtspten_p: cobj,
    slflx_p: cobj,
    qtflx_p: cobj,
    wtflx_p: cobj,
    slten_p: cobj,
    qtten_p: cobj,
    wttotten_p: cobj,
    rliq_p: cobj,
    rainflx_p: cobj,
    snowflx_p: cobj,
):
    rliq = Ptr[float](rliq_p)
    rainflx = Ptr[float](rainflx_p)
    snowflx = Ptr[float](snowflx_p)

    rliq[0] = 0.0
    rainflx[0] = 0.0
    snowflx[0] = 0.0

    uwshcu_comp_sub_sink_shell_codon(
        mkx,
        ncnst,
        wtrc_nwset,
        kpen,
        ixnumliq,
        ixnumice,
        dt_v,
        ql0_p,
        qi0_p,
        tr0_p,
        wtrc_iatype_p,
        qlten_sub_p,
        qiten_sub_p,
        nlten_sub_p,
        niten_sub_p,
        wtlten_sub_p,
        wtiten_sub_p,
        qlten_sink_p,
        qiten_sink_p,
        nlten_sink_p,
        niten_sink_p,
        wtten_sink_liq_p,
        wtten_sink_ice_p,
    )
    uwshcu_tendency_prep_shell_codon(
        mkx,
        wtrc_nwset,
        kpen,
        frc_rasn,
        g_v,
        dt_v,
        dp0_p,
        u0_p,
        v0_p,
        uflx_p,
        vflx_p,
        umf_p,
        uten_p,
        vten_p,
        uf_p,
        vf_p,
        dwten_p,
        diten_p,
        wtdwten_p,
        wtditen_p,
        qrten_p,
        qsten_p,
        wtrpten_p,
        wtspten_p,
        slflx_p,
        qtflx_p,
        wtflx_p,
        slten_p,
        qtten_p,
        wttotten_p,
        rainflx_p,
        snowflx_p,
    )


@export
def uwshcu_thermo_final_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    k_fortran: int,
    use_expconten: int,
    use_unicondet: int,
    ixnumliq: int,
    ixnumice: int,
    frc_rasn: float,
    dt_v: float,
    xlv_v: float,
    xls_v: float,
    g_v: float,
    qc_lm_v: float,
    qc_im_v: float,
    nc_lm_v: float,
    nc_im_v: float,
    dp0_p: cobj,
    qt0_p: cobj,
    ql0_p: cobj,
    qi0_p: cobj,
    dwten_p: cobj,
    diten_p: cobj,
    qtten_p: cobj,
    slten_p: cobj,
    qlten_sink_p: cobj,
    qiten_sink_p: cobj,
    nlten_sink_p: cobj,
    niten_sink_p: cobj,
    qc_l_p: cobj,
    qc_i_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    qvten_p: cobj,
    sten_p: cobj,
    tr0_p: cobj,
    trten_p: cobj,
    wtrc_iatype_p: cobj,
    wt0_p: cobj,
    wtdwten_p: cobj,
    wtditen_p: cobj,
    wttotten_p: cobj,
    wtten_sink_liq_p: cobj,
    wtten_sink_ice_p: cobj,
    wtqc_liq_p: cobj,
    wtqc_ice_p: cobj,
    wtqcm_liq_p: cobj,
    wtqcm_ice_p: cobj,
    wtlten_det_p: cobj,
    wtiten_det_p: cobj,
    qc_p: cobj,
    rliq_p: cobj,
):
    dp0 = Ptr[float](dp0_p)
    qt0 = Ptr[float](qt0_p)
    ql0 = Ptr[float](ql0_p)
    qi0 = Ptr[float](qi0_p)
    dwten = Ptr[float](dwten_p)
    diten = Ptr[float](diten_p)
    qtten = Ptr[float](qtten_p)
    slten = Ptr[float](slten_p)
    qlten_sink = Ptr[float](qlten_sink_p)
    qiten_sink = Ptr[float](qiten_sink_p)
    nlten_sink = Ptr[float](nlten_sink_p)
    niten_sink = Ptr[float](niten_sink_p)
    qc_l = Ptr[float](qc_l_p)
    qc_i = Ptr[float](qc_i_p)
    qlten = Ptr[float](qlten_p)
    qiten = Ptr[float](qiten_p)
    qvten = Ptr[float](qvten_p)
    sten = Ptr[float](sten_p)
    tr0 = Ptr[float](tr0_p)
    trten = Ptr[float](trten_p)
    wtrc_iatype = Ptr[int](wtrc_iatype_p)
    wt0 = Ptr[float](wt0_p)
    wtdwten = Ptr[float](wtdwten_p)
    wtditen = Ptr[float](wtditen_p)
    wttotten = Ptr[float](wttotten_p)
    wtten_sink_liq = Ptr[float](wtten_sink_liq_p)
    wtten_sink_ice = Ptr[float](wtten_sink_ice_p)
    wtqc_liq = Ptr[float](wtqc_liq_p)
    wtqc_ice = Ptr[float](wtqc_ice_p)
    wtqcm_liq = Ptr[float](wtqcm_liq_p)
    wtqcm_ice = Ptr[float](wtqcm_ice_p)
    wtlten_det = Ptr[float](wtlten_det_p)
    wtiten_det = Ptr[float](wtiten_det_p)
    qc = Ptr[float](qc_p)
    rliq = Ptr[float](rliq_p)

    k = k_fortran - 1
    liq_idx = ixnumliq - 1
    ice_idx = ixnumice - 1
    qlten_det = qc_l[k] + qc_lm_v
    qiten_det = qc_i[k] + qc_im_v

    m = 0
    while m < wtrc_nwset:
        wt_idx = k + m * mkx
        wtlten_det[wt_idx] = wtqc_liq[wt_idx] + wtqcm_liq[m]
        wtiten_det[wt_idx] = wtqc_ice[wt_idx] + wtqcm_ice[m]
        m += 1

    if use_expconten != 0:
        if use_unicondet != 0:
            qc_l[k] = 0.0
            qc_i[k] = 0.0
            qlten[k] = frc_rasn * dwten[k] + qlten_sink[k] + qlten_det
            qiten[k] = frc_rasn * diten[k] + qiten_sink[k] + qiten_det
            m = 0
            while m < wtrc_nwset:
                wt_idx = k + m * mkx
                liq = wtrc_iatype[m + wtrc_nwset] - 1
                ice = wtrc_iatype[m + 2 * wtrc_nwset] - 1
                wtqc_liq[wt_idx] = 0.0
                wtqc_ice[wt_idx] = 0.0
                trten[k + liq * mkx] = frc_rasn * wtdwten[wt_idx] + wtten_sink_liq[wt_idx] + wtlten_det[wt_idx]
                trten[k + ice * mkx] = frc_rasn * wtditen[wt_idx] + wtten_sink_ice[wt_idx] + wtiten_det[wt_idx]
                m += 1
        else:
            ql_tmp = ql0[k] + (qc_lm_v + qlten_sink[k]) * dt_v
            if ql_tmp > 0.0:
                ql_pos = ql_tmp
            else:
                ql_pos = 0.0
            qi_tmp = qi0[k] + (qc_im_v + qiten_sink[k]) * dt_v
            if qi_tmp > 0.0:
                qi_pos = qi_tmp
            else:
                qi_pos = 0.0
            qlten[k] = qc_l[k] + frc_rasn * dwten[k] + (ql_pos - ql0[k]) / dt_v
            qiten[k] = qc_i[k] + frc_rasn * diten[k] + (qi_pos - qi0[k]) / dt_v

            nl_val = nc_lm_v + nlten_sink[k]
            nl_floor = -tr0[k + liq_idx * mkx] / dt_v
            if nl_val > nl_floor:
                trten[k + liq_idx * mkx] = nl_val
            else:
                trten[k + liq_idx * mkx] = nl_floor
            ni_val = nc_im_v + niten_sink[k]
            ni_floor = -tr0[k + ice_idx * mkx] / dt_v
            if ni_val > ni_floor:
                trten[k + ice_idx * mkx] = ni_val
            else:
                trten[k + ice_idx * mkx] = ni_floor

            m = 0
            while m < wtrc_nwset:
                wt_idx = k + m * mkx
                liq = wtrc_iatype[m + wtrc_nwset] - 1
                ice = wtrc_iatype[m + 2 * wtrc_nwset] - 1
                wt_liq_tmp = tr0[k + liq * mkx] + (wtqcm_liq[m] + wtten_sink_liq[wt_idx]) * dt_v
                if wt_liq_tmp > 0.0:
                    wt_liq_pos = wt_liq_tmp
                else:
                    wt_liq_pos = 0.0
                wt_ice_tmp = tr0[k + ice * mkx] + (wtqcm_ice[m] + wtten_sink_ice[wt_idx]) * dt_v
                if wt_ice_tmp > 0.0:
                    wt_ice_pos = wt_ice_tmp
                else:
                    wt_ice_pos = 0.0
                trten[k + liq * mkx] = wtqc_liq[wt_idx] + frc_rasn * wtdwten[wt_idx] + (
                    wt_liq_pos - tr0[k + liq * mkx]
                ) / dt_v
                trten[k + ice * mkx] = wtqc_ice[wt_idx] + frc_rasn * wtditen[wt_idx] + (
                    wt_ice_pos - tr0[k + ice * mkx]
                ) / dt_v
                m += 1
    else:
        if use_unicondet != 0:
            qc_l[k] = 0.0
            qc_i[k] = 0.0
            m = 0
            while m < wtrc_nwset:
                wt_idx = k + m * mkx
                wtqc_liq[wt_idx] = 0.0
                wtqc_ice[wt_idx] = 0.0
                m += 1
        qlten[k] = dwten[k] + (qtten[k] - dwten[k] - diten[k]) * (ql0[k] / qt0[k])
        qiten[k] = diten[k] + (qtten[k] - dwten[k] - diten[k]) * (qi0[k] / qt0[k])
        m = 0
        while m < wtrc_nwset:
            wt_idx = k + m * mkx
            liq = wtrc_iatype[m + wtrc_nwset] - 1
            ice = wtrc_iatype[m + 2 * wtrc_nwset] - 1
            trten[k + liq * mkx] = wtdwten[wt_idx] + (wttotten[wt_idx] - wtdwten[wt_idx] - wtditen[wt_idx]) * (
                tr0[k + liq * mkx] / wt0[wt_idx]
            )
            trten[k + ice * mkx] = wtditen[wt_idx] + (wttotten[wt_idx] - wtdwten[wt_idx] - wtditen[wt_idx]) * (
                tr0[k + ice * mkx] / wt0[wt_idx]
            )
            m += 1

    qvten[k] = qtten[k] - qlten[k] - qiten[k]
    sten[k] = slten[k] + xlv_v * qlten[k] + xls_v * qiten[k]

    m = 0
    while m < wtrc_nwset:
        wt_idx = k + m * mkx
        vap = wtrc_iatype[m] - 1
        liq = wtrc_iatype[m + wtrc_nwset] - 1
        ice = wtrc_iatype[m + 2 * wtrc_nwset] - 1
        trten[k + vap * mkx] = wttotten[wt_idx] - trten[k + liq * mkx] - trten[k + ice * mkx]
        m += 1

    qc[k] = qc_l[k] + qc_i[k]
    rliq[0] = rliq[0] + qc[k] * dp0[k] / g_v / 1000.0


@export
def uwshcu_thermo_final_precip_bulk_init_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    k_fortran: int,
    kpen: int,
    use_expconten: int,
    use_unicondet: int,
    ixnumliq: int,
    ixnumice: int,
    trace_water: int,
    frc_rasn: float,
    dt_v: float,
    xlv_v: float,
    xls_v: float,
    g_v: float,
    qc_lm_v: float,
    qc_im_v: float,
    nc_lm_v: float,
    nc_im_v: float,
    rainflx: float,
    snowflx: float,
    dp0_p: cobj,
    qt0_p: cobj,
    ql0_p: cobj,
    qi0_p: cobj,
    dwten_p: cobj,
    diten_p: cobj,
    qtten_p: cobj,
    slten_p: cobj,
    qlten_sink_p: cobj,
    qiten_sink_p: cobj,
    nlten_sink_p: cobj,
    niten_sink_p: cobj,
    qc_l_p: cobj,
    qc_i_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    qvten_p: cobj,
    sten_p: cobj,
    tr0_p: cobj,
    trten_p: cobj,
    wtrc_iatype_p: cobj,
    wt0_p: cobj,
    wtdwten_p: cobj,
    wtditen_p: cobj,
    wttotten_p: cobj,
    wtten_sink_liq_p: cobj,
    wtten_sink_ice_p: cobj,
    wtqc_liq_p: cobj,
    wtqc_ice_p: cobj,
    wtqcm_liq_p: cobj,
    wtqcm_ice_p: cobj,
    wtlten_det_p: cobj,
    wtiten_det_p: cobj,
    qc_p: cobj,
    rliq_p: cobj,
    precip_p: cobj,
    snow_p: cobj,
    evpint_rain_p: cobj,
    evpint_snow_p: cobj,
    flxrain_p: cobj,
    flxsnow_p: cobj,
    ntraprd_p: cobj,
    ntsnprd_p: cobj,
    wtflxrn_p: cobj,
    wtflxsn_p: cobj,
):
    uwshcu_thermo_final_shell_codon(
        mkx, ncnst, wtrc_nwset, k_fortran, use_expconten,
        use_unicondet, ixnumliq, ixnumice, frc_rasn, dt_v, xlv_v,
        xls_v, g_v, qc_lm_v, qc_im_v, nc_lm_v, nc_im_v, dp0_p,
        qt0_p, ql0_p, qi0_p, dwten_p, diten_p, qtten_p, slten_p,
        qlten_sink_p, qiten_sink_p, nlten_sink_p, niten_sink_p,
        qc_l_p, qc_i_p, qlten_p, qiten_p, qvten_p, sten_p, tr0_p,
        trten_p, wtrc_iatype_p, wt0_p, wtdwten_p, wtditen_p,
        wttotten_p, wtten_sink_liq_p, wtten_sink_ice_p, wtqc_liq_p,
        wtqc_ice_p, wtqcm_liq_p, wtqcm_ice_p, wtlten_det_p,
        wtiten_det_p, qc_p, rliq_p,
    )
    if k_fortran == kpen:
        uwshcu_precip_bulk_init_shell_codon(
            mkx, wtrc_nwset, trace_water, rainflx, snowflx, precip_p,
            snow_p, evpint_rain_p, evpint_snow_p, flxrain_p,
            flxsnow_p, ntraprd_p, ntsnprd_p, wtflxrn_p, wtflxsn_p,
        )


@export
def uwshcu_reserved_condensate_adjust_shell_codon(
    mkx: int,
    wtrc_nwset: int,
    kpen: int,
    xlv_v: float,
    xls_v: float,
    qtten_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    slten_p: cobj,
    sten_p: cobj,
    qc_p: cobj,
    qc_l_p: cobj,
    qc_i_p: cobj,
    trten_p: cobj,
    wtrc_iatype_p: cobj,
    wtqc_liq_p: cobj,
    wtqc_ice_p: cobj,
):
    qtten = Ptr[float](qtten_p)
    qlten = Ptr[float](qlten_p)
    qiten = Ptr[float](qiten_p)
    slten = Ptr[float](slten_p)
    sten = Ptr[float](sten_p)
    qc = Ptr[float](qc_p)
    qc_l = Ptr[float](qc_l_p)
    qc_i = Ptr[float](qc_i_p)
    trten = Ptr[float](trten_p)
    wtrc_iatype = Ptr[int](wtrc_iatype_p)
    wtqc_liq = Ptr[float](wtqc_liq_p)
    wtqc_ice = Ptr[float](wtqc_ice_p)

    k = 0
    while k < kpen:
        qtten[k] = qtten[k] - qc[k]
        qlten[k] = qlten[k] - qc_l[k]
        qiten[k] = qiten[k] - qc_i[k]
        slten[k] = slten[k] + (xlv_v * qc_l[k] + xls_v * qc_i[k])
        sten[k] = sten[k] - (xls_v - xlv_v) * qc_i[k]

        m = 0
        while m < wtrc_nwset:
            wt_idx = k + m * mkx
            liq = wtrc_iatype[m + wtrc_nwset] - 1
            ice = wtrc_iatype[m + 2 * wtrc_nwset] - 1
            trten[k + liq * mkx] = trten[k + liq * mkx] - wtqc_liq[wt_idx]
            trten[k + ice * mkx] = trten[k + ice * mkx] - wtqc_ice[wt_idx]
            m += 1

        k += 1


@export
def uwshcu_post_positive_thermo_shell_codon(
    mkx: int,
    xlv_v: float,
    xls_v: float,
    qvten_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    sten_p: cobj,
    qtten_p: cobj,
    slten_p: cobj,
):
    qvten = Ptr[float](qvten_p)
    qlten = Ptr[float](qlten_p)
    qiten = Ptr[float](qiten_p)
    sten = Ptr[float](sten_p)
    qtten = Ptr[float](qtten_p)
    slten = Ptr[float](slten_p)

    k = 0
    while k < mkx:
        qtten[k] = qvten[k] + qlten[k] + qiten[k]
        slten[k] = sten[k] - xlv_v * qlten[k] - xls_v * qiten[k]
        k += 1


@export
def uwshcu_tracer_limiter_shell_codon(
    mkx: int,
    ncnst: int,
    g_v: float,
    dt_v: float,
    ixnumliq: int,
    ixnumice: int,
    dp0_p: cobj,
    dpdry0_p: cobj,
    tr0_p: cobj,
    trflx_p: cobj,
    trten_p: cobj,
    trflx_d_p: cobj,
    trflx_u_p: cobj,
    qmin_p: cobj,
    is_water_p: cobj,
    wet_p: cobj,
):
    dp0 = Ptr[float](dp0_p)
    dpdry0 = Ptr[float](dpdry0_p)
    tr0 = Ptr[float](tr0_p)
    trflx = Ptr[float](trflx_p)
    trten = Ptr[float](trten_p)
    trflx_d = Ptr[float](trflx_d_p)
    trflx_u = Ptr[float](trflx_u_p)
    qmin = Ptr[float](qmin_p)
    is_water = Ptr[int](is_water_p)
    wet = Ptr[int](wet_p)

    sp = 3
    while sp < ncnst:
        m_fortran = sp + 1
        if m_fortran != ixnumliq and m_fortran != ixnumice and is_water[sp] == 0:
            trmin = qmin[sp]

            iface = 0
            while iface <= mkx:
                trflx_d[iface] = 0.0
                trflx_u[iface] = 0.0
                iface += 1

            k = 0
            while k < mkx - 1:
                if wet[sp] != 0:
                    pdelx = dp0[k]
                else:
                    pdelx = dpdry0[k]
                km1 = k
                col = sp * (mkx + 1)
                dum = (tr0[k + sp * mkx] - trmin) * pdelx / g_v / dt_v + trflx[km1 + col] - trflx[k + 1 + col] + trflx_d[km1]
                if dum < 0.0:
                    trflx_d[k + 1] = dum
                else:
                    trflx_d[k + 1] = 0.0
                k += 1

            k = mkx - 1
            while k >= 1:
                if wet[sp] != 0:
                    pdelx = dp0[k]
                else:
                    pdelx = dpdry0[k]
                km1 = k
                col = sp * (mkx + 1)
                dum = (tr0[k + sp * mkx] - trmin) * pdelx / g_v / dt_v + trflx[km1 + col] - trflx[k + 1 + col] + trflx_d[km1] - trflx_d[k + 1] - trflx_u[k + 1]
                if -dum > 0.0:
                    trflx_u[km1] = -dum
                else:
                    trflx_u[km1] = 0.0
                k -= 1

            k = 0
            while k < mkx:
                if wet[sp] != 0:
                    pdelx = dp0[k]
                else:
                    pdelx = dpdry0[k]
                km1 = k
                col = sp * (mkx + 1)
                trten[k + sp * mkx] = (
                    trflx[km1 + col] - trflx[k + 1 + col]
                    + trflx_d[km1] - trflx_d[k + 1]
                    + trflx_u[km1] - trflx_u[k + 1]
                ) * g_v / pdelx
                k += 1
        sp += 1


@export
def uwshcu_cloud_diag_init_shell_codon(
    qlj_v: float,
    qij_v: float,
    qcubelow_p: cobj,
    qlubelow_p: cobj,
    qiubelow_p: cobj,
    rcwp_p: cobj,
    rlwp_p: cobj,
    riwp_p: cobj,
):
    qcubelow = Ptr[float](qcubelow_p)
    qlubelow = Ptr[float](qlubelow_p)
    qiubelow = Ptr[float](qiubelow_p)
    rcwp = Ptr[float](rcwp_p)
    rlwp = Ptr[float](rlwp_p)
    riwp = Ptr[float](riwp_p)

    qcubelow[0] = qlj_v + qij_v
    qlubelow[0] = qlj_v
    qiubelow[0] = qij_v
    rcwp[0] = 0.0
    rlwp[0] = 0.0
    riwp[0] = 0.0


@export
def uwshcu_cloud_diag_layer_shell_codon(
    mkx: int,
    k_fortran: int,
    krel: int,
    kpen: int,
    qlj_v: float,
    qij_v: float,
    criqc_v: float,
    prel_v: float,
    ppen_v: float,
    ufrclcl_v: float,
    g_v: float,
    ps0_p: cobj,
    ufrc_p: cobj,
    qcu_p: cobj,
    qlu_p: cobj,
    qiu_p: cobj,
    cufrc_p: cobj,
    qcubelow_p: cobj,
    qlubelow_p: cobj,
    qiubelow_p: cobj,
    rcwp_p: cobj,
    rlwp_p: cobj,
    riwp_p: cobj,
    cnt_p: cobj,
    cnb_p: cobj,
):
    ps0 = Ptr[float](ps0_p)
    ufrc = Ptr[float](ufrc_p)
    qcu = Ptr[float](qcu_p)
    qlu = Ptr[float](qlu_p)
    qiu = Ptr[float](qiu_p)
    cufrc = Ptr[float](cufrc_p)
    qcubelow = Ptr[float](qcubelow_p)
    qlubelow = Ptr[float](qlubelow_p)
    qiubelow = Ptr[float](qiubelow_p)
    rcwp = Ptr[float](rcwp_p)
    rlwp = Ptr[float](rlwp_p)
    riwp = Ptr[float](riwp_p)
    cnt = Ptr[float](cnt_p)
    cnb = Ptr[float](cnb_p)

    k = k_fortran
    idx = k - 1
    qcu[idx] = 0.5 * (qcubelow[0] + qlj_v + qij_v)
    qlu[idx] = 0.5 * (qlubelow[0] + qlj_v)
    qiu[idx] = 0.5 * (qiubelow[0] + qij_v)
    cufrc[idx] = ufrc[k - 1] + ufrc[k]
    if k == krel:
        cufrc[idx] = (ufrclcl_v + ufrc[k]) * (prel_v - ps0[k]) / (ps0[k - 1] - ps0[k])
    elif k == kpen:
        cufrc[idx] = (ufrc[k - 1] + 0.0) * (-ppen_v) / (ps0[k - 1] - ps0[k])
        if qlj_v + qij_v > criqc_v:
            qcu[idx] = 0.5 * (qcubelow[0] + criqc_v)
            qlu[idx] = 0.5 * (qlubelow[0] + criqc_v * qlj_v / (qlj_v + qij_v))
            qiu[idx] = 0.5 * (qiubelow[0] + criqc_v * qij_v / (qlj_v + qij_v))

    rcwp[0] = rcwp[0] + (qlu[idx] + qiu[idx]) * (ps0[k - 1] - ps0[k]) / g_v * cufrc[idx]
    rlwp[0] = rlwp[0] + qlu[idx] * (ps0[k - 1] - ps0[k]) / g_v * cufrc[idx]
    riwp[0] = riwp[0] + qiu[idx] * (ps0[k - 1] - ps0[k]) / g_v * cufrc[idx]
    qcubelow[0] = qlj_v + qij_v
    qlubelow[0] = qlj_v
    qiubelow[0] = qij_v
    if k == kpen:
        cnt[0] = float(kpen)
        cnb[0] = float(krel - 1)


@export
def uwshcu_cloud_diag_all_shell_codon(
    mkx: int,
    krel: int,
    kpen: int,
    qlj_base: float,
    qij_base: float,
    criqc_v: float,
    prel_v: float,
    ppen_v: float,
    ufrclcl_v: float,
    g_v: float,
    cloud_qlj_p: cobj,
    cloud_qij_p: cobj,
    ps0_p: cobj,
    ufrc_p: cobj,
    qcu_p: cobj,
    qlu_p: cobj,
    qiu_p: cobj,
    cufrc_p: cobj,
    qcubelow_p: cobj,
    qlubelow_p: cobj,
    qiubelow_p: cobj,
    rcwp_p: cobj,
    rlwp_p: cobj,
    riwp_p: cobj,
    cnt_p: cobj,
    cnb_p: cobj,
):
    cloud_qlj = Ptr[float](cloud_qlj_p)
    cloud_qij = Ptr[float](cloud_qij_p)

    uwshcu_cloud_diag_init_shell_codon(
        qlj_base,
        qij_base,
        qcubelow_p,
        qlubelow_p,
        qiubelow_p,
        rcwp_p,
        rlwp_p,
        riwp_p,
    )

    k = krel
    while k <= kpen:
        idx = k - 1
        uwshcu_cloud_diag_layer_shell_codon(
            mkx,
            k,
            krel,
            kpen,
            cloud_qlj[idx],
            cloud_qij[idx],
            criqc_v,
            prel_v,
            ppen_v,
            ufrclcl_v,
            g_v,
            ps0_p,
            ufrc_p,
            qcu_p,
            qlu_p,
            qiu_p,
            cufrc_p,
            qcubelow_p,
            qlubelow_p,
            qiubelow_p,
            rcwp_p,
            rlwp_p,
            riwp_p,
            cnt_p,
            cnb_p,
        )
        k += 1


@export
def uwshcu_cloud_diag_index_shell_codon(
    kpen: int,
    krel: int,
    cnt_p: cobj,
    cnb_p: cobj,
):
    cnt = Ptr[float](cnt_p)
    cnb = Ptr[float](cnb_p)

    cnt[0] = float(kpen)
    cnb[0] = float(krel - 1)


@export
def uwshcu_cloud_diag_conden_exit_shell_codon(
    id_check: int,
    exit_conden_p: cobj,
    exit_code_p: cobj,
):
    exit_conden = Ptr[float](exit_conden_p)
    exit_code = Ptr[int](exit_code_p)

    exit_code[0] = 0
    if id_check == 1:
        exit_conden[0] = 1.0
        exit_code[0] = 1


@export
def uwshcu_positive_moisture_prep_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    dt_v: float,
    qv0_p: cobj,
    ql0_p: cobj,
    qi0_p: cobj,
    s0_p: cobj,
    qvten_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    sten_p: cobj,
    tr0_p: cobj,
    trten_p: cobj,
    wtrc_iatype_p: cobj,
    qv0_star_p: cobj,
    ql0_star_p: cobj,
    qi0_star_p: cobj,
    s0_star_p: cobj,
    wt0_star_p: cobj,
):
    qv0 = Ptr[float](qv0_p)
    ql0 = Ptr[float](ql0_p)
    qi0 = Ptr[float](qi0_p)
    s0 = Ptr[float](s0_p)
    qvten = Ptr[float](qvten_p)
    qlten = Ptr[float](qlten_p)
    qiten = Ptr[float](qiten_p)
    sten = Ptr[float](sten_p)
    tr0 = Ptr[float](tr0_p)
    trten = Ptr[float](trten_p)
    wtrc_iatype = Ptr[int](wtrc_iatype_p)
    qv0_star = Ptr[float](qv0_star_p)
    ql0_star = Ptr[float](ql0_star_p)
    qi0_star = Ptr[float](qi0_star_p)
    s0_star = Ptr[float](s0_star_p)
    wt0_star = Ptr[float](wt0_star_p)

    k = 0
    while k < mkx:
        qv0_star[k] = qv0[k] + qvten[k] * dt_v
        ql0_star[k] = ql0[k] + qlten[k] * dt_v
        qi0_star[k] = qi0[k] + qiten[k] * dt_v
        s0_star[k] = s0[k] + sten[k] * dt_v
        k += 1

    if wtrc_nwset > 0:
        comp = 0
        while comp < 3:
            m = 0
            while m < wtrc_nwset:
                k = 0
                while k < mkx:
                    wt0_star[k + m * mkx + comp * mkx * wtrc_nwset] = 0.0
                    k += 1
                m += 1
            comp += 1

        m = 0
        while m < wtrc_nwset:
            vap = wtrc_iatype[m] - 1
            liq = wtrc_iatype[m + wtrc_nwset] - 1
            ice = wtrc_iatype[m + 2 * wtrc_nwset] - 1
            k = 0
            while k < mkx:
                wt0_star[k + m * mkx] = tr0[k + vap * mkx] + trten[k + vap * mkx] * dt_v
                wt0_star[k + m * mkx + mkx * wtrc_nwset] = tr0[k + liq * mkx] + trten[k + liq * mkx] * dt_v
                wt0_star[k + m * mkx + 2 * mkx * wtrc_nwset] = tr0[k + ice * mkx] + trten[k + ice * mkx] * dt_v
                k += 1
            m += 1


@export
def uwshcu_post_precip_positive_prep_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    kpen: int,
    xlv_v: float,
    xls_v: float,
    dt_v: float,
    qtten_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    slten_p: cobj,
    sten_p: cobj,
    qc_p: cobj,
    qc_l_p: cobj,
    qc_i_p: cobj,
    trten_p: cobj,
    wtrc_iatype_p: cobj,
    wtqc_liq_p: cobj,
    wtqc_ice_p: cobj,
    qv0_p: cobj,
    ql0_p: cobj,
    qi0_p: cobj,
    s0_p: cobj,
    qvten_p: cobj,
    tr0_p: cobj,
    qv0_star_p: cobj,
    ql0_star_p: cobj,
    qi0_star_p: cobj,
    s0_star_p: cobj,
    wt0_star_p: cobj,
):
    uwshcu_reserved_condensate_adjust_shell_codon(
        mkx,
        wtrc_nwset,
        kpen,
        xlv_v,
        xls_v,
        qtten_p,
        qlten_p,
        qiten_p,
        slten_p,
        sten_p,
        qc_p,
        qc_l_p,
        qc_i_p,
        trten_p,
        wtrc_iatype_p,
        wtqc_liq_p,
        wtqc_ice_p,
    )
    uwshcu_positive_moisture_prep_shell_codon(
        mkx,
        ncnst,
        wtrc_nwset,
        dt_v,
        qv0_p,
        ql0_p,
        qi0_p,
        s0_p,
        qvten_p,
        qlten_p,
        qiten_p,
        sten_p,
        tr0_p,
        trten_p,
        wtrc_iatype_p,
        qv0_star_p,
        ql0_star_p,
        qi0_star_p,
        s0_star_p,
        wt0_star_p,
    )


@export
def uwshcu_post_positive_tracer_limiter_shell_codon(
    mkx: int,
    ncnst: int,
    g_v: float,
    dt_v: float,
    ixnumliq: int,
    ixnumice: int,
    xlv_v: float,
    xls_v: float,
    qvten_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    sten_p: cobj,
    qtten_p: cobj,
    slten_p: cobj,
    dp0_p: cobj,
    dpdry0_p: cobj,
    tr0_p: cobj,
    trflx_p: cobj,
    trten_p: cobj,
    trflx_d_p: cobj,
    trflx_u_p: cobj,
    qmin_p: cobj,
    is_water_p: cobj,
    wet_p: cobj,
):
    uwshcu_post_positive_thermo_shell_codon(
        mkx,
        xlv_v,
        xls_v,
        qvten_p,
        qlten_p,
        qiten_p,
        sten_p,
        qtten_p,
        slten_p,
    )
    uwshcu_tracer_limiter_shell_codon(
        mkx,
        ncnst,
        g_v,
        dt_v,
        ixnumliq,
        ixnumice,
        dp0_p,
        dpdry0_p,
        tr0_p,
        trflx_p,
        trten_p,
        trflx_d_p,
        trflx_u_p,
        qmin_p,
        is_water_p,
        wet_p,
    )


@export
def uwshcu_precip_surface_finalize_shell_codon(
    mkx: int,
    wtrc_nwset: int,
    flxrain_p: cobj,
    flxsnow_p: cobj,
    wtflxrn_p: cobj,
    wtflxsn_p: cobj,
    precip_p: cobj,
    snow_p: cobj,
    wtprec_p: cobj,
    wtsnow_p: cobj,
):
    flxrain = Ptr[float](flxrain_p)
    flxsnow = Ptr[float](flxsnow_p)
    wtflxrn = Ptr[float](wtflxrn_p)
    wtflxsn = Ptr[float](wtflxsn_p)
    precip = Ptr[float](precip_p)
    snow = Ptr[float](snow_p)
    wtprec = Ptr[float](wtprec_p)
    wtsnow = Ptr[float](wtsnow_p)

    precip[0] = (flxrain[0] + flxsnow[0]) / 1000.0
    snow[0] = flxsnow[0] / 1000.0

    m = 0
    while m < wtrc_nwset:
        wtprec[m] = (wtflxrn[m * (mkx + 1)] + wtflxsn[m * (mkx + 1)]) / 1000.0
        wtsnow[m] = wtflxsn[m * (mkx + 1)] / 1000.0
        m += 1


@export
def uwshcu_precip_surface_positive_prep_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    kpen: int,
    xlv_v: float,
    xls_v: float,
    dt_v: float,
    flxrain_p: cobj,
    flxsnow_p: cobj,
    wtflxrn_p: cobj,
    wtflxsn_p: cobj,
    precip_p: cobj,
    snow_p: cobj,
    wtprec_p: cobj,
    wtsnow_p: cobj,
    qtten_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    slten_p: cobj,
    sten_p: cobj,
    qc_p: cobj,
    qc_l_p: cobj,
    qc_i_p: cobj,
    trten_p: cobj,
    wtrc_iatype_p: cobj,
    wtqc_liq_p: cobj,
    wtqc_ice_p: cobj,
    qv0_p: cobj,
    ql0_p: cobj,
    qi0_p: cobj,
    s0_p: cobj,
    qvten_p: cobj,
    tr0_p: cobj,
    qv0_star_p: cobj,
    ql0_star_p: cobj,
    qi0_star_p: cobj,
    s0_star_p: cobj,
    wt0_star_p: cobj,
):
    uwshcu_precip_surface_finalize_shell_codon(
        mkx, wtrc_nwset, flxrain_p, flxsnow_p, wtflxrn_p, wtflxsn_p,
        precip_p, snow_p, wtprec_p, wtsnow_p,
    )
    uwshcu_post_precip_positive_prep_shell_codon(
        mkx, ncnst, wtrc_nwset, kpen, xlv_v, xls_v, dt_v, qtten_p,
        qlten_p, qiten_p, slten_p, sten_p, qc_p, qc_l_p, qc_i_p,
        trten_p, wtrc_iatype_p, wtqc_liq_p, wtqc_ice_p, qv0_p, ql0_p,
        qi0_p, s0_p, qvten_p, tr0_p, qv0_star_p, ql0_star_p,
        qi0_star_p, s0_star_p, wt0_star_p,
    )


@export
def uwshcu_precip_bulk_init_shell_codon(
    mkx: int,
    wtrc_nwset: int,
    trace_water: int,
    rainflx: float,
    snowflx: float,
    precip_p: cobj,
    snow_p: cobj,
    evpint_rain_p: cobj,
    evpint_snow_p: cobj,
    flxrain_p: cobj,
    flxsnow_p: cobj,
    ntraprd_p: cobj,
    ntsnprd_p: cobj,
    wtflxrn_p: cobj,
    wtflxsn_p: cobj,
):
    precip = Ptr[float](precip_p)
    snow = Ptr[float](snow_p)
    evpint_rain = Ptr[float](evpint_rain_p)
    evpint_snow = Ptr[float](evpint_snow_p)
    flxrain = Ptr[float](flxrain_p)
    flxsnow = Ptr[float](flxsnow_p)
    ntraprd = Ptr[float](ntraprd_p)
    ntsnprd = Ptr[float](ntsnprd_p)
    wtflxrn = Ptr[float](wtflxrn_p)
    wtflxsn = Ptr[float](wtflxsn_p)

    precip[0] = rainflx + snowflx
    snow[0] = snowflx
    evpint_rain[0] = 0.0
    evpint_snow[0] = 0.0

    iface = 0
    while iface <= mkx:
        flxrain[iface] = 0.0
        flxsnow[iface] = 0.0
        iface += 1

    k = 0
    while k < mkx:
        ntraprd[k] = 0.0
        ntsnprd[k] = 0.0
        k += 1

    if trace_water != 0:
        m = 0
        while m < wtrc_nwset:
            iface = 0
            while iface <= mkx:
                idx = iface + m * (mkx + 1)
                wtflxrn[idx] = 0.0
                wtflxsn[idx] = 0.0
                iface += 1
            m += 1


@export
def uwshcu_precip_bulk_layer_shell_codon(
    mkx: int,
    mix: int,
    i: int,
    k: int,
    wtrc_nwset: int,
    trace_water: int,
    t0_v: float,
    rainflx: float,
    snowflx: float,
    snowmlt: float,
    evprain: float,
    evpsnow: float,
    g: float,
    dt: float,
    xlv: float,
    xls: float,
    qmin_vap: float,
    qmin_liq: float,
    qmin_ice: float,
    dp0_p: cobj,
    qv0_p: cobj,
    ql0_p: cobj,
    qi0_p: cobj,
    qrten_p: cobj,
    qsten_p: cobj,
    evapc_p: cobj,
    evpint_rain_p: cobj,
    evpint_snow_p: cobj,
    ntraprd_p: cobj,
    ntsnprd_p: cobj,
    flxrain_p: cobj,
    flxsnow_p: cobj,
    qvten_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    qtten_p: cobj,
    sten_p: cobj,
    slten_p: cobj,
    limit_negcon_p: cobj,
    wtrc_iatype_p: cobj,
    wtrpten_p: cobj,
    wtspten_p: cobj,
    wtevp_p: cobj,
    wtsub_p: cobj,
    wtflxrn_p: cobj,
    wtflxsn_p: cobj,
    trten_p: cobj,
):
    dp0 = Ptr[float](dp0_p)
    qv0 = Ptr[float](qv0_p)
    ql0 = Ptr[float](ql0_p)
    qi0 = Ptr[float](qi0_p)
    qrten = Ptr[float](qrten_p)
    qsten = Ptr[float](qsten_p)
    evapc = Ptr[float](evapc_p)
    evpint_rain = Ptr[float](evpint_rain_p)
    evpint_snow = Ptr[float](evpint_snow_p)
    ntraprd = Ptr[float](ntraprd_p)
    ntsnprd = Ptr[float](ntsnprd_p)
    flxrain = Ptr[float](flxrain_p)
    flxsnow = Ptr[float](flxsnow_p)
    qvten = Ptr[float](qvten_p)
    qlten = Ptr[float](qlten_p)
    qiten = Ptr[float](qiten_p)
    qtten = Ptr[float](qtten_p)
    sten = Ptr[float](sten_p)
    slten = Ptr[float](slten_p)
    limit_negcon = Ptr[float](limit_negcon_p)
    wtrc_iatype = Ptr[int](wtrc_iatype_p)
    wtrpten = Ptr[float](wtrpten_p)
    wtspten = Ptr[float](wtspten_p)
    wtevp = Ptr[float](wtevp_p)
    wtsub = Ptr[float](wtsub_p)
    wtflxrn = Ptr[float](wtflxrn_p)
    wtflxsn = Ptr[float](wtflxsn_p)
    trten = Ptr[float](trten_p)

    kk = k - 1
    km1 = k - 1

    evapc[kk] = evprain + evpsnow

    evpint_rain[0] = evpint_rain[0] + evprain * dp0[kk] / g
    evpint_snow[0] = evpint_snow[0] + evpsnow * dp0[kk] / g

    ntraprd[kk] = qrten[kk] - evprain + snowmlt
    ntsnprd[kk] = qsten[kk] - evpsnow - snowmlt

    flxrain[km1] = flxrain[k] + ntraprd[kk] * dp0[kk] / g
    flxsnow[km1] = flxsnow[k] + ntsnprd[kk] * dp0[kk] / g
    if flxrain[km1] < 0.0:
        flxrain[km1] = 0.0
    if flxrain[km1] == 0.0:
        ntraprd[kk] = -flxrain[k] * g / dp0[kk]
    if flxsnow[km1] < 0.0:
        flxsnow[km1] = 0.0
    if flxsnow[km1] == 0.0:
        ntsnprd[kk] = -flxsnow[k] * g / dp0[kk]

    qlten[kk] = qlten[kk] - qrten[kk]
    qiten[kk] = qiten[kk] - qsten[kk]
    qvten[kk] = qvten[kk] + evprain + evpsnow
    qtten[kk] = qlten[kk] + qiten[kk] + qvten[kk]
    if (
        (qv0[kk] + qvten[kk] * dt) < qmin_vap
        or (ql0[kk] + qlten[kk] * dt) < qmin_liq
        or (qi0[kk] + qiten[kk] * dt) < qmin_ice
    ):
        limit_negcon[i - 1] = 1.0
    sten[kk] = sten[kk] - xlv * evprain - xls * evpsnow - (xls - xlv) * snowmlt
    slten[kk] = sten[kk] - xlv * qlten[kk] - xls * qiten[kk]

    if trace_water != 0:
        iface_stride = mkx + 1
        m = 0
        while m < wtrc_nwset:
            wt_layer = kk + m * mkx
            wt_top = k + m * iface_stride
            wt_bottom = kk + m * iface_stride

            if t0_v > 273.16:
                wtsnwmlt = wtflxsn[wt_top] * g / dp0[kk]
                if wtsnwmlt < 0.0:
                    wtsnwmlt = 0.0
            else:
                wtsnwmlt = 0.0

            wtflxrn[wt_bottom] = wtflxrn[wt_top] + (wtrpten[wt_layer] - wtevp[wt_layer] + wtsnwmlt) * dp0[kk] / g
            wtflxsn[wt_bottom] = wtflxsn[wt_top] + (wtspten[wt_layer] - wtsub[wt_layer] - wtsnwmlt) * dp0[kk] / g
            if wtflxrn[wt_bottom] < 0.0:
                wtflxrn[wt_bottom] = 0.0
            if wtflxsn[wt_bottom] < 0.0:
                wtflxsn[wt_bottom] = 0.0

            vap = wtrc_iatype[m] - 1
            liq = wtrc_iatype[m + wtrc_nwset] - 1
            ice = wtrc_iatype[m + 2 * wtrc_nwset] - 1
            trten[kk + liq * mkx] = trten[kk + liq * mkx] - wtrpten[wt_layer]
            trten[kk + ice * mkx] = trten[kk + ice * mkx] - wtspten[wt_layer]
            trten[kk + vap * mkx] = trten[kk + vap * mkx] + wtevp[wt_layer] + wtsub[wt_layer]
            m += 1


@export
def uwshcu_thermo_post_batch_shell_codon(
    kind: int,
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    k_fortran: int,
    kpen: int,
    mix: int,
    i_col: int,
    use_expconten: int,
    use_unicondet: int,
    ixnumliq: int,
    ixnumice: int,
    trace_water: int,
    frc_rasn: float,
    dt_v: float,
    xlv_v: float,
    xls_v: float,
    g_v: float,
    qc_lm_v: float,
    qc_im_v: float,
    nc_lm_v: float,
    nc_im_v: float,
    rainflx_v: float,
    snowflx_v: float,
    t0_v: float,
    snowmlt_v: float,
    evprain_v: float,
    evpsnow_v: float,
    qmin_vap: float,
    qmin_liq: float,
    qmin_ice: float,
    dp0_p: cobj,
    qv0_p: cobj,
    qt0_p: cobj,
    ql0_p: cobj,
    qi0_p: cobj,
    s0_p: cobj,
    qrten_p: cobj,
    qsten_p: cobj,
    dwten_p: cobj,
    diten_p: cobj,
    qtten_p: cobj,
    slten_p: cobj,
    qlten_sink_p: cobj,
    qiten_sink_p: cobj,
    nlten_sink_p: cobj,
    niten_sink_p: cobj,
    qc_l_p: cobj,
    qc_i_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    qvten_p: cobj,
    sten_p: cobj,
    tr0_p: cobj,
    trten_p: cobj,
    wtrc_iatype_p: cobj,
    wt0_p: cobj,
    wtdwten_p: cobj,
    wtditen_p: cobj,
    wttotten_p: cobj,
    wtten_sink_liq_p: cobj,
    wtten_sink_ice_p: cobj,
    wtqc_liq_p: cobj,
    wtqc_ice_p: cobj,
    wtqcm_liq_p: cobj,
    wtqcm_ice_p: cobj,
    wtlten_det_p: cobj,
    wtiten_det_p: cobj,
    qc_p: cobj,
    rliq_p: cobj,
    precip_p: cobj,
    snow_p: cobj,
    evapc_p: cobj,
    evpint_rain_p: cobj,
    evpint_snow_p: cobj,
    flxrain_p: cobj,
    flxsnow_p: cobj,
    ntraprd_p: cobj,
    ntsnprd_p: cobj,
    limit_negcon_p: cobj,
    wtrpten_p: cobj,
    wtspten_p: cobj,
    wtevp_p: cobj,
    wtsub_p: cobj,
    wtflxrn_p: cobj,
    wtflxsn_p: cobj,
    wtprec_p: cobj,
    wtsnow_p: cobj,
    qv0_star_p: cobj,
    ql0_star_p: cobj,
    qi0_star_p: cobj,
    s0_star_p: cobj,
    wt0_star_p: cobj,
    dpdry0_p: cobj,
    trflx_p: cobj,
    trflx_d_p: cobj,
    trflx_u_p: cobj,
    qmin_p: cobj,
    is_water_p: cobj,
    wet_p: cobj,
):
    if kind == 0:
        uwshcu_thermo_final_precip_bulk_init_shell_codon(
            mkx,
            ncnst,
            wtrc_nwset,
            k_fortran,
            kpen,
            use_expconten,
            use_unicondet,
            ixnumliq,
            ixnumice,
            trace_water,
            frc_rasn,
            dt_v,
            xlv_v,
            xls_v,
            g_v,
            qc_lm_v,
            qc_im_v,
            nc_lm_v,
            nc_im_v,
            rainflx_v,
            snowflx_v,
            dp0_p,
            qt0_p,
            ql0_p,
            qi0_p,
            dwten_p,
            diten_p,
            qtten_p,
            slten_p,
            qlten_sink_p,
            qiten_sink_p,
            nlten_sink_p,
            niten_sink_p,
            qc_l_p,
            qc_i_p,
            qlten_p,
            qiten_p,
            qvten_p,
            sten_p,
            tr0_p,
            trten_p,
            wtrc_iatype_p,
            wt0_p,
            wtdwten_p,
            wtditen_p,
            wttotten_p,
            wtten_sink_liq_p,
            wtten_sink_ice_p,
            wtqc_liq_p,
            wtqc_ice_p,
            wtqcm_liq_p,
            wtqcm_ice_p,
            wtlten_det_p,
            wtiten_det_p,
            qc_p,
            rliq_p,
            precip_p,
            snow_p,
            evpint_rain_p,
            evpint_snow_p,
            flxrain_p,
            flxsnow_p,
            ntraprd_p,
            ntsnprd_p,
            wtflxrn_p,
            wtflxsn_p,
        )
    elif kind == 1:
        uwshcu_precip_bulk_layer_shell_codon(
            mkx,
            mix,
            i_col,
            k_fortran,
            wtrc_nwset,
            trace_water,
            t0_v,
            rainflx_v,
            snowflx_v,
            snowmlt_v,
            evprain_v,
            evpsnow_v,
            g_v,
            dt_v,
            xlv_v,
            xls_v,
            qmin_vap,
            qmin_liq,
            qmin_ice,
            dp0_p,
            qv0_p,
            ql0_p,
            qi0_p,
            qrten_p,
            qsten_p,
            evapc_p,
            evpint_rain_p,
            evpint_snow_p,
            ntraprd_p,
            ntsnprd_p,
            flxrain_p,
            flxsnow_p,
            qvten_p,
            qlten_p,
            qiten_p,
            qtten_p,
            sten_p,
            slten_p,
            limit_negcon_p,
            wtrc_iatype_p,
            wtrpten_p,
            wtspten_p,
            wtevp_p,
            wtsub_p,
            wtflxrn_p,
            wtflxsn_p,
            trten_p,
        )
    elif kind == 2:
        uwshcu_precip_surface_positive_prep_shell_codon(
            mkx,
            ncnst,
            wtrc_nwset,
            kpen,
            xlv_v,
            xls_v,
            dt_v,
            flxrain_p,
            flxsnow_p,
            wtflxrn_p,
            wtflxsn_p,
            precip_p,
            snow_p,
            wtprec_p,
            wtsnow_p,
            qtten_p,
            qlten_p,
            qiten_p,
            slten_p,
            sten_p,
            qc_p,
            qc_l_p,
            qc_i_p,
            trten_p,
            wtrc_iatype_p,
            wtqc_liq_p,
            wtqc_ice_p,
            qv0_p,
            ql0_p,
            qi0_p,
            s0_p,
            qvten_p,
            tr0_p,
            qv0_star_p,
            ql0_star_p,
            qi0_star_p,
            s0_star_p,
            wt0_star_p,
        )
    elif kind == 3:
        uwshcu_post_positive_tracer_limiter_shell_codon(
            mkx,
            ncnst,
            g_v,
            dt_v,
            ixnumliq,
            ixnumice,
            xlv_v,
            xls_v,
            qvten_p,
            qlten_p,
            qiten_p,
            sten_p,
            qtten_p,
            slten_p,
            dp0_p,
            dpdry0_p,
            tr0_p,
            trflx_p,
            trten_p,
            trflx_d_p,
            trflx_u_p,
            qmin_p,
            is_water_p,
            wet_p,
        )


def _uwshcu_slope_column(
    mkx: int,
    p0: Ptr[float],
    field: Ptr[float],
    field_offset: int,
    field_stride: int,
    out: Ptr[float],
    out_offset: int,
    out_stride: int,
):
    below = (field[field_offset + field_stride] - field[field_offset]) / (p0[1] - p0[0])

    k = 2
    while k <= mkx:
        cur_idx = field_offset + (k - 1) * field_stride
        prev_idx = field_offset + (k - 2) * field_stride
        above = (field[cur_idx] - field[prev_idx]) / (p0[k - 1] - p0[k - 2])
        out_idx = out_offset + (k - 2) * out_stride
        if above > 0.0:
            out[out_idx] = max(0.0, min(above, below))
        else:
            out[out_idx] = min(0.0, max(above, below))
        below = above
        k += 1

    out[out_offset + (mkx - 1) * out_stride] = out[out_offset + (mkx - 2) * out_stride]


@export
def uwshcu_slope_reconstruction_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    p0_p: cobj,
    thl0_p: cobj,
    qt0_p: cobj,
    u0_p: cobj,
    v0_p: cobj,
    tr0_p: cobj,
    wt0_p: cobj,
    ssthl0_p: cobj,
    ssqt0_p: cobj,
    ssu0_p: cobj,
    ssv0_p: cobj,
    sstr0_p: cobj,
    sswt0_p: cobj,
):
    p0 = Ptr[float](p0_p)
    thl0 = Ptr[float](thl0_p)
    qt0 = Ptr[float](qt0_p)
    u0 = Ptr[float](u0_p)
    v0 = Ptr[float](v0_p)
    tr0 = Ptr[float](tr0_p)
    wt0 = Ptr[float](wt0_p)
    ssthl0 = Ptr[float](ssthl0_p)
    ssqt0 = Ptr[float](ssqt0_p)
    ssu0 = Ptr[float](ssu0_p)
    ssv0 = Ptr[float](ssv0_p)
    sstr0 = Ptr[float](sstr0_p)
    sswt0 = Ptr[float](sswt0_p)

    _uwshcu_slope_column(mkx, p0, thl0, 0, 1, ssthl0, 0, 1)
    _uwshcu_slope_column(mkx, p0, qt0, 0, 1, ssqt0, 0, 1)
    _uwshcu_slope_column(mkx, p0, u0, 0, 1, ssu0, 0, 1)
    _uwshcu_slope_column(mkx, p0, v0, 0, 1, ssv0, 0, 1)

    m = 0
    while m < ncnst:
        _uwshcu_slope_column(mkx, p0, tr0, m * mkx, 1, sstr0, m * mkx, 1)
        m += 1

    m = 0
    while m < wtrc_nwset:
        _uwshcu_slope_column(mkx, p0, wt0, m * mkx, 1, sswt0, m * mkx, 1)
        m += 1


@export
def uwshcu_column_thermo_slope_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    xlv: float,
    xls: float,
    cp: float,
    zvir: float,
    qv0_p: cobj,
    ql0_p: cobj,
    qi0_p: cobj,
    t0_p: cobj,
    exn0_p: cobj,
    tr0_p: cobj,
    wtrc_iatype_p: cobj,
    p0_p: cobj,
    u0_p: cobj,
    v0_p: cobj,
    qt0_p: cobj,
    thl0_p: cobj,
    thvl0_p: cobj,
    wt0_p: cobj,
    ssthl0_p: cobj,
    ssqt0_p: cobj,
    ssu0_p: cobj,
    ssv0_p: cobj,
    sstr0_p: cobj,
    sswt0_p: cobj,
):
    uwshcu_column_thermo_state_shell_codon(
        mkx, ncnst, wtrc_nwset, xlv, xls, cp, zvir, qv0_p, ql0_p, qi0_p,
        t0_p, exn0_p, tr0_p, wtrc_iatype_p, qt0_p, thl0_p, thvl0_p, wt0_p
    )
    uwshcu_slope_reconstruction_shell_codon(
        mkx, ncnst, wtrc_nwset, p0_p, thl0_p, qt0_p, u0_p, v0_p, tr0_p,
        wt0_p, ssthl0_p, ssqt0_p, ssu0_p, ssv0_p, sstr0_p, sswt0_p
    )


@export
def uwshcu_column_input_load_shell_codon(
    mix: int,
    mkx: int,
    i_col: int,
    ncnst: int,
    pblh_v: float,
    cush_v: float,
    ps0_in_p: cobj,
    zs0_in_p: cobj,
    p0_in_p: cobj,
    z0_in_p: cobj,
    dp0_in_p: cobj,
    dpdry0_in_p: cobj,
    u0_in_p: cobj,
    v0_in_p: cobj,
    qv0_in_p: cobj,
    ql0_in_p: cobj,
    qi0_in_p: cobj,
    t0_in_p: cobj,
    s0_in_p: cobj,
    tke_in_p: cobj,
    cldfrct_in_p: cobj,
    concldfrct_in_p: cobj,
    tr0_in_p: cobj,
    ps0_p: cobj,
    zs0_p: cobj,
    p0_p: cobj,
    z0_p: cobj,
    dp0_p: cobj,
    dpdry0_p: cobj,
    u0_p: cobj,
    v0_p: cobj,
    qv0_p: cobj,
    ql0_p: cobj,
    qi0_p: cobj,
    t0_p: cobj,
    s0_p: cobj,
    tke_p: cobj,
    cldfrct_p: cobj,
    concldfrct_p: cobj,
    tr0_p: cobj,
    pblh_p: cobj,
    cush_p: cobj,
):
    ps0_in = Ptr[float](ps0_in_p)
    zs0_in = Ptr[float](zs0_in_p)
    p0_in = Ptr[float](p0_in_p)
    z0_in = Ptr[float](z0_in_p)
    dp0_in = Ptr[float](dp0_in_p)
    dpdry0_in = Ptr[float](dpdry0_in_p)
    u0_in = Ptr[float](u0_in_p)
    v0_in = Ptr[float](v0_in_p)
    qv0_in = Ptr[float](qv0_in_p)
    ql0_in = Ptr[float](ql0_in_p)
    qi0_in = Ptr[float](qi0_in_p)
    t0_in = Ptr[float](t0_in_p)
    s0_in = Ptr[float](s0_in_p)
    tke_in = Ptr[float](tke_in_p)
    cldfrct_in = Ptr[float](cldfrct_in_p)
    concldfrct_in = Ptr[float](concldfrct_in_p)
    tr0_in = Ptr[float](tr0_in_p)
    ps0 = Ptr[float](ps0_p)
    zs0 = Ptr[float](zs0_p)
    p0 = Ptr[float](p0_p)
    z0 = Ptr[float](z0_p)
    dp0 = Ptr[float](dp0_p)
    dpdry0 = Ptr[float](dpdry0_p)
    u0 = Ptr[float](u0_p)
    v0 = Ptr[float](v0_p)
    qv0 = Ptr[float](qv0_p)
    ql0 = Ptr[float](ql0_p)
    qi0 = Ptr[float](qi0_p)
    t0 = Ptr[float](t0_p)
    s0 = Ptr[float](s0_p)
    tke = Ptr[float](tke_p)
    cldfrct = Ptr[float](cldfrct_p)
    concldfrct = Ptr[float](concldfrct_p)
    tr0 = Ptr[float](tr0_p)
    pblh = Ptr[float](pblh_p)
    cush = Ptr[float](cush_p)

    col = i_col - 1

    pblh[0] = pblh_v
    cush[0] = cush_v

    k = 0
    while k <= mkx:
        src = col + k * mix
        ps0[k] = ps0_in[src]
        zs0[k] = zs0_in[src]
        tke[k] = tke_in[src]
        k += 1

    k = 0
    while k < mkx:
        src = col + k * mix
        p0[k] = p0_in[src]
        z0[k] = z0_in[src]
        dp0[k] = dp0_in[src]
        dpdry0[k] = dpdry0_in[src]
        u0[k] = u0_in[src]
        v0[k] = v0_in[src]
        qv0[k] = qv0_in[src]
        ql0[k] = ql0_in[src]
        qi0[k] = qi0_in[src]
        t0[k] = t0_in[src]
        s0[k] = s0_in[src]
        cldfrct[k] = cldfrct_in[src]
        concldfrct[k] = concldfrct_in[src]
        k += 1

    m = 0
    while m < ncnst:
        src_offset = m * mix * mkx
        dst_offset = m * mkx
        k = 0
        while k < mkx:
            tr0[k + dst_offset] = tr0_in[col + k * mix + src_offset]
            k += 1
        m += 1


@export
def uwshcu_column_env_save_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    qv0_p: cobj,
    ql0_p: cobj,
    qi0_p: cobj,
    t0_p: cobj,
    s0_p: cobj,
    u0_p: cobj,
    v0_p: cobj,
    qt0_p: cobj,
    thl0_p: cobj,
    thvl0_p: cobj,
    ssthl0_p: cobj,
    ssqt0_p: cobj,
    thv0bot_p: cobj,
    thv0top_p: cobj,
    thvl0bot_p: cobj,
    thvl0top_p: cobj,
    ssu0_p: cobj,
    ssv0_p: cobj,
    tr0_p: cobj,
    sstr0_p: cobj,
    sswt0_p: cobj,
    qv0_o_p: cobj,
    ql0_o_p: cobj,
    qi0_o_p: cobj,
    t0_o_p: cobj,
    s0_o_p: cobj,
    u0_o_p: cobj,
    v0_o_p: cobj,
    qt0_o_p: cobj,
    thl0_o_p: cobj,
    thvl0_o_p: cobj,
    ssthl0_o_p: cobj,
    ssqt0_o_p: cobj,
    thv0bot_o_p: cobj,
    thv0top_o_p: cobj,
    thvl0bot_o_p: cobj,
    thvl0top_o_p: cobj,
    ssu0_o_p: cobj,
    ssv0_o_p: cobj,
    tr0_o_p: cobj,
    sstr0_o_p: cobj,
    sswt0_o_p: cobj,
):
    qv0 = Ptr[float](qv0_p)
    ql0 = Ptr[float](ql0_p)
    qi0 = Ptr[float](qi0_p)
    t0 = Ptr[float](t0_p)
    s0 = Ptr[float](s0_p)
    u0 = Ptr[float](u0_p)
    v0 = Ptr[float](v0_p)
    qt0 = Ptr[float](qt0_p)
    thl0 = Ptr[float](thl0_p)
    thvl0 = Ptr[float](thvl0_p)
    ssthl0 = Ptr[float](ssthl0_p)
    ssqt0 = Ptr[float](ssqt0_p)
    thv0bot = Ptr[float](thv0bot_p)
    thv0top = Ptr[float](thv0top_p)
    thvl0bot = Ptr[float](thvl0bot_p)
    thvl0top = Ptr[float](thvl0top_p)
    ssu0 = Ptr[float](ssu0_p)
    ssv0 = Ptr[float](ssv0_p)
    tr0 = Ptr[float](tr0_p)
    sstr0 = Ptr[float](sstr0_p)
    sswt0 = Ptr[float](sswt0_p)
    qv0_o = Ptr[float](qv0_o_p)
    ql0_o = Ptr[float](ql0_o_p)
    qi0_o = Ptr[float](qi0_o_p)
    t0_o = Ptr[float](t0_o_p)
    s0_o = Ptr[float](s0_o_p)
    u0_o = Ptr[float](u0_o_p)
    v0_o = Ptr[float](v0_o_p)
    qt0_o = Ptr[float](qt0_o_p)
    thl0_o = Ptr[float](thl0_o_p)
    thvl0_o = Ptr[float](thvl0_o_p)
    ssthl0_o = Ptr[float](ssthl0_o_p)
    ssqt0_o = Ptr[float](ssqt0_o_p)
    thv0bot_o = Ptr[float](thv0bot_o_p)
    thv0top_o = Ptr[float](thv0top_o_p)
    thvl0bot_o = Ptr[float](thvl0bot_o_p)
    thvl0top_o = Ptr[float](thvl0top_o_p)
    ssu0_o = Ptr[float](ssu0_o_p)
    ssv0_o = Ptr[float](ssv0_o_p)
    tr0_o = Ptr[float](tr0_o_p)
    sstr0_o = Ptr[float](sstr0_o_p)
    sswt0_o = Ptr[float](sswt0_o_p)

    k = 0
    while k < mkx:
        qv0_o[k] = qv0[k]
        ql0_o[k] = ql0[k]
        qi0_o[k] = qi0[k]
        t0_o[k] = t0[k]
        s0_o[k] = s0[k]
        u0_o[k] = u0[k]
        v0_o[k] = v0[k]
        qt0_o[k] = qt0[k]
        thl0_o[k] = thl0[k]
        thvl0_o[k] = thvl0[k]
        ssthl0_o[k] = ssthl0[k]
        ssqt0_o[k] = ssqt0[k]
        thv0bot_o[k] = thv0bot[k]
        thv0top_o[k] = thv0top[k]
        thvl0bot_o[k] = thvl0bot[k]
        thvl0top_o[k] = thvl0top[k]
        ssu0_o[k] = ssu0[k]
        ssv0_o[k] = ssv0[k]
        k += 1

    m = 0
    while m < ncnst:
        offset = m * mkx
        k = 0
        while k < mkx:
            idx = k + offset
            tr0_o[idx] = tr0[idx]
            sstr0_o[idx] = sstr0[idx]
            k += 1
        m += 1

    m = 0
    while m < wtrc_nwset:
        offset = m * mkx
        k = 0
        while k < mkx:
            idx = k + offset
            sswt0_o[idx] = sswt0[idx]
            k += 1
        m += 1


@export
def uwshcu_column_init_all_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    qv0_p: cobj,
    ql0_p: cobj,
    qi0_p: cobj,
    t0_p: cobj,
    s0_p: cobj,
    u0_p: cobj,
    v0_p: cobj,
    qt0_p: cobj,
    thl0_p: cobj,
    thvl0_p: cobj,
    ssthl0_p: cobj,
    ssqt0_p: cobj,
    thv0bot_p: cobj,
    thv0top_p: cobj,
    thvl0bot_p: cobj,
    thvl0top_p: cobj,
    ssu0_p: cobj,
    ssv0_p: cobj,
    tr0_p: cobj,
    sstr0_p: cobj,
    sswt0_p: cobj,
    qv0_o_p: cobj,
    ql0_o_p: cobj,
    qi0_o_p: cobj,
    t0_o_p: cobj,
    s0_o_p: cobj,
    u0_o_p: cobj,
    v0_o_p: cobj,
    qt0_o_p: cobj,
    thl0_o_p: cobj,
    thvl0_o_p: cobj,
    ssthl0_o_p: cobj,
    ssqt0_o_p: cobj,
    thv0bot_o_p: cobj,
    thv0top_o_p: cobj,
    thvl0bot_o_p: cobj,
    thvl0top_o_p: cobj,
    ssu0_o_p: cobj,
    ssv0_o_p: cobj,
    tr0_o_p: cobj,
    sstr0_o_p: cobj,
    sswt0_o_p: cobj,
    umf_p: cobj,
    emf_p: cobj,
    slflx_p: cobj,
    qtflx_p: cobj,
    uflx_p: cobj,
    vflx_p: cobj,
    qvten_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    sten_p: cobj,
    uten_p: cobj,
    vten_p: cobj,
    qrten_p: cobj,
    qsten_p: cobj,
    dwten_p: cobj,
    diten_p: cobj,
    evapc_p: cobj,
    cufrc_p: cobj,
    qcu_p: cobj,
    qlu_p: cobj,
    qiu_p: cobj,
    fer_p: cobj,
    fdr_p: cobj,
    qc_p: cobj,
    qc_l_p: cobj,
    qc_i_p: cobj,
    qtten_p: cobj,
    slten_p: cobj,
    ufrc_p: cobj,
    thlu_p: cobj,
    qtu_p: cobj,
    uu_p: cobj,
    vu_p: cobj,
    wu_p: cobj,
    thvu_p: cobj,
    thlu_emf_p: cobj,
    qtu_emf_p: cobj,
    uu_emf_p: cobj,
    vu_emf_p: cobj,
    trflx_p: cobj,
    trten_p: cobj,
    tru_p: cobj,
    tru_emf_p: cobj,
    wtdwten_p: cobj,
    wtditen_p: cobj,
    wtrpten_p: cobj,
    wtspten_p: cobj,
    wtqc_liq_p: cobj,
    wtqc_ice_p: cobj,
    wtu_p: cobj,
    wtu_emf_p: cobj,
    wtflx_p: cobj,
    wttotten_p: cobj,
    wtprec_p: cobj,
    wtsnow_p: cobj,
    excessu_p: cobj,
    excess0_p: cobj,
    xc_p: cobj,
    aquad_p: cobj,
    bquad_p: cobj,
    cquad_p: cobj,
    bogbot_p: cobj,
    bogtop_p: cobj,
    comsub_p: cobj,
    qlten_sink_p: cobj,
    qiten_sink_p: cobj,
    nlten_sink_p: cobj,
    niten_sink_p: cobj,
    wtten_sink_liq_p: cobj,
    wtten_sink_ice_p: cobj,
    wtevp_p: cobj,
    wtsub_p: cobj,
    dz_p: cobj,
    uemf_p: cobj,
    precip_p: cobj,
    snow_p: cobj,
    cin_p: cobj,
    cbmf_p: cobj,
    rliq_p: cobj,
    cnt_p: cobj,
    cnb_p: cobj,
    ufrcinvbase_p: cobj,
    ufrclcl_p: cobj,
    winvbase_p: cobj,
    wlcl_p: cobj,
    emfkbup_p: cobj,
    cbmflimit_p: cobj,
):
    uwshcu_column_env_save_shell_codon(
        mkx, ncnst, wtrc_nwset, qv0_p, ql0_p, qi0_p, t0_p, s0_p, u0_p, v0_p,
        qt0_p, thl0_p, thvl0_p, ssthl0_p, ssqt0_p, thv0bot_p, thv0top_p,
        thvl0bot_p, thvl0top_p, ssu0_p, ssv0_p, tr0_p, sstr0_p, sswt0_p,
        qv0_o_p, ql0_o_p, qi0_o_p, t0_o_p, s0_o_p, u0_o_p, v0_o_p, qt0_o_p,
        thl0_o_p, thvl0_o_p, ssthl0_o_p, ssqt0_o_p, thv0bot_o_p, thv0top_o_p,
        thvl0bot_o_p, thvl0top_o_p, ssu0_o_p, ssv0_o_p, tr0_o_p, sstr0_o_p,
        sswt0_o_p,
    )
    uwshcu_initial_workspace_reset_shell_codon(
        mkx, ncnst, wtrc_nwset, umf_p, emf_p, slflx_p, qtflx_p, uflx_p, vflx_p,
        qvten_p, qlten_p, qiten_p, sten_p, uten_p, vten_p, qrten_p, qsten_p,
        dwten_p, diten_p, evapc_p, cufrc_p, qcu_p, qlu_p, qiu_p, fer_p, fdr_p,
        qc_p, qc_l_p, qc_i_p, qtten_p, slten_p, ufrc_p, thlu_p, qtu_p, uu_p,
        vu_p, wu_p, thvu_p, thlu_emf_p, qtu_emf_p, uu_emf_p, vu_emf_p,
        trflx_p, trten_p, tru_p, tru_emf_p, wtdwten_p, wtditen_p, wtrpten_p,
        wtspten_p, wtqc_liq_p, wtqc_ice_p, wtu_p, wtu_emf_p, wtflx_p,
        wttotten_p, wtprec_p, wtsnow_p, excessu_p, excess0_p, xc_p, aquad_p,
        bquad_p, cquad_p, bogbot_p, bogtop_p, comsub_p, qlten_sink_p,
        qiten_sink_p, nlten_sink_p, niten_sink_p, wtten_sink_liq_p,
        wtten_sink_ice_p, wtevp_p, wtsub_p, dz_p, uemf_p, precip_p, snow_p,
        cin_p, cbmf_p, rliq_p, cnt_p, cnb_p, ufrcinvbase_p, ufrclcl_p,
        winvbase_p, wlcl_p, emfkbup_p, cbmflimit_p,
    )


@export
def uwshcu_column_extra_workspace_reset_shell_codon(
    mkx: int,
    wtrc_nwset: int,
    comsub_p: cobj,
    qlten_sink_p: cobj,
    qiten_sink_p: cobj,
    nlten_sink_p: cobj,
    niten_sink_p: cobj,
    wtten_sink_liq_p: cobj,
    wtten_sink_ice_p: cobj,
    wtevp_p: cobj,
    wtsub_p: cobj,
    dz_p: cobj,
    uemf_p: cobj,
):
    comsub = Ptr[float](comsub_p)
    qlten_sink = Ptr[float](qlten_sink_p)
    qiten_sink = Ptr[float](qiten_sink_p)
    nlten_sink = Ptr[float](nlten_sink_p)
    niten_sink = Ptr[float](niten_sink_p)
    wtten_sink_liq = Ptr[float](wtten_sink_liq_p)
    wtten_sink_ice = Ptr[float](wtten_sink_ice_p)
    wtevp = Ptr[float](wtevp_p)
    wtsub = Ptr[float](wtsub_p)
    dz = Ptr[float](dz_p)
    uemf = Ptr[float](uemf_p)

    k = 0
    while k < mkx:
        comsub[k] = 0.0
        qlten_sink[k] = 0.0
        qiten_sink[k] = 0.0
        nlten_sink[k] = 0.0
        niten_sink[k] = 0.0
        k += 1

    k = 0
    while k <= mkx:
        dz[k] = 0.0
        uemf[k] = 0.0
        k += 1

    m = 0
    while m < wtrc_nwset:
        offset = m * mkx
        k = 0
        while k < mkx:
            idx = k + offset
            wtten_sink_liq[idx] = 0.0
            wtten_sink_ice[idx] = 0.0
            wtevp[idx] = 0.0
            wtsub[idx] = 0.0
            k += 1
        m += 1


@export
def uwshcu_initial_workspace_reset_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    umf_p: cobj,
    emf_p: cobj,
    slflx_p: cobj,
    qtflx_p: cobj,
    uflx_p: cobj,
    vflx_p: cobj,
    qvten_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    sten_p: cobj,
    uten_p: cobj,
    vten_p: cobj,
    qrten_p: cobj,
    qsten_p: cobj,
    dwten_p: cobj,
    diten_p: cobj,
    evapc_p: cobj,
    cufrc_p: cobj,
    qcu_p: cobj,
    qlu_p: cobj,
    qiu_p: cobj,
    fer_p: cobj,
    fdr_p: cobj,
    qc_p: cobj,
    qc_l_p: cobj,
    qc_i_p: cobj,
    qtten_p: cobj,
    slten_p: cobj,
    ufrc_p: cobj,
    thlu_p: cobj,
    qtu_p: cobj,
    uu_p: cobj,
    vu_p: cobj,
    wu_p: cobj,
    thvu_p: cobj,
    thlu_emf_p: cobj,
    qtu_emf_p: cobj,
    uu_emf_p: cobj,
    vu_emf_p: cobj,
    trflx_p: cobj,
    trten_p: cobj,
    tru_p: cobj,
    tru_emf_p: cobj,
    wtdwten_p: cobj,
    wtditen_p: cobj,
    wtrpten_p: cobj,
    wtspten_p: cobj,
    wtqc_liq_p: cobj,
    wtqc_ice_p: cobj,
    wtu_p: cobj,
    wtu_emf_p: cobj,
    wtflx_p: cobj,
    wttotten_p: cobj,
    wtprec_p: cobj,
    wtsnow_p: cobj,
    excessu_p: cobj,
    excess0_p: cobj,
    xc_p: cobj,
    aquad_p: cobj,
    bquad_p: cobj,
    cquad_p: cobj,
    bogbot_p: cobj,
    bogtop_p: cobj,
    comsub_p: cobj,
    qlten_sink_p: cobj,
    qiten_sink_p: cobj,
    nlten_sink_p: cobj,
    niten_sink_p: cobj,
    wtten_sink_liq_p: cobj,
    wtten_sink_ice_p: cobj,
    wtevp_p: cobj,
    wtsub_p: cobj,
    dz_p: cobj,
    uemf_p: cobj,
    precip_p: cobj,
    snow_p: cobj,
    cin_p: cobj,
    cbmf_p: cobj,
    rliq_p: cobj,
    cnt_p: cobj,
    cnb_p: cobj,
    ufrcinvbase_p: cobj,
    ufrclcl_p: cobj,
    winvbase_p: cobj,
    wlcl_p: cobj,
    emfkbup_p: cobj,
    cbmflimit_p: cobj,
):
    uwshcu_delcin_workspace_reset_shell_codon(
        mkx,
        ncnst,
        wtrc_nwset,
        umf_p,
        emf_p,
        slflx_p,
        qtflx_p,
        uflx_p,
        vflx_p,
        qvten_p,
        qlten_p,
        qiten_p,
        sten_p,
        uten_p,
        vten_p,
        qrten_p,
        qsten_p,
        dwten_p,
        diten_p,
        evapc_p,
        cufrc_p,
        qcu_p,
        qlu_p,
        qiu_p,
        fer_p,
        fdr_p,
        qc_p,
        qc_l_p,
        qc_i_p,
        qtten_p,
        slten_p,
        ufrc_p,
        thlu_p,
        qtu_p,
        uu_p,
        vu_p,
        wu_p,
        thvu_p,
        thlu_emf_p,
        qtu_emf_p,
        uu_emf_p,
        vu_emf_p,
        trflx_p,
        trten_p,
        tru_p,
        tru_emf_p,
        wtdwten_p,
        wtditen_p,
        wtrpten_p,
        wtspten_p,
        wtqc_liq_p,
        wtqc_ice_p,
        wtu_p,
        wtu_emf_p,
        wtflx_p,
        wttotten_p,
        wtprec_p,
        wtsnow_p,
        excessu_p,
        excess0_p,
        xc_p,
        aquad_p,
        bquad_p,
        cquad_p,
        bogbot_p,
        bogtop_p,
    )
    uwshcu_column_extra_workspace_reset_shell_codon(
        mkx,
        wtrc_nwset,
        comsub_p,
        qlten_sink_p,
        qiten_sink_p,
        nlten_sink_p,
        niten_sink_p,
        wtten_sink_liq_p,
        wtten_sink_ice_p,
        wtevp_p,
        wtsub_p,
        dz_p,
        uemf_p,
    )

    Ptr[float](precip_p)[0] = 0.0
    Ptr[float](snow_p)[0] = 0.0
    Ptr[float](cin_p)[0] = 0.0
    Ptr[float](cbmf_p)[0] = 0.0
    Ptr[float](rliq_p)[0] = 0.0
    Ptr[float](cnt_p)[0] = float(mkx)
    Ptr[float](cnb_p)[0] = 0.0
    Ptr[float](ufrcinvbase_p)[0] = 0.0
    Ptr[float](ufrclcl_p)[0] = 0.0
    Ptr[float](winvbase_p)[0] = 0.0
    Ptr[float](wlcl_p)[0] = 0.0
    Ptr[float](emfkbup_p)[0] = 0.0
    Ptr[float](cbmflimit_p)[0] = 0.0


@export
def uwshcu_iter_save_env_shell_codon(
    mkx: int,
    ncnst: int,
    dt: float,
    cp: float,
    qv0_p: cobj,
    ql0_p: cobj,
    qi0_p: cobj,
    s0_p: cobj,
    u0_p: cobj,
    v0_p: cobj,
    t0_p: cobj,
    qvten_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    sten_p: cobj,
    uten_p: cobj,
    vten_p: cobj,
    tr0_p: cobj,
    trten_p: cobj,
    qv0_s_p: cobj,
    ql0_s_p: cobj,
    qi0_s_p: cobj,
    s0_s_p: cobj,
    u0_s_p: cobj,
    v0_s_p: cobj,
    qt0_s_p: cobj,
    t0_s_p: cobj,
    tr0_s_p: cobj,
):
    qv0 = Ptr[float](qv0_p)
    ql0 = Ptr[float](ql0_p)
    qi0 = Ptr[float](qi0_p)
    s0 = Ptr[float](s0_p)
    u0 = Ptr[float](u0_p)
    v0 = Ptr[float](v0_p)
    t0 = Ptr[float](t0_p)
    qvten = Ptr[float](qvten_p)
    qlten = Ptr[float](qlten_p)
    qiten = Ptr[float](qiten_p)
    sten = Ptr[float](sten_p)
    uten = Ptr[float](uten_p)
    vten = Ptr[float](vten_p)
    tr0 = Ptr[float](tr0_p)
    trten = Ptr[float](trten_p)
    qv0_s = Ptr[float](qv0_s_p)
    ql0_s = Ptr[float](ql0_s_p)
    qi0_s = Ptr[float](qi0_s_p)
    s0_s = Ptr[float](s0_s_p)
    u0_s = Ptr[float](u0_s_p)
    v0_s = Ptr[float](v0_s_p)
    qt0_s = Ptr[float](qt0_s_p)
    t0_s = Ptr[float](t0_s_p)
    tr0_s = Ptr[float](tr0_s_p)

    k = 0
    while k < mkx:
        qv0_s[k] = qv0[k] + qvten[k] * dt
        ql0_s[k] = ql0[k] + qlten[k] * dt
        qi0_s[k] = qi0[k] + qiten[k] * dt
        s0_s[k] = s0[k] + sten[k] * dt
        u0_s[k] = u0[k] + uten[k] * dt
        v0_s[k] = v0[k] + vten[k] * dt
        qt0_s[k] = qv0_s[k] + ql0_s[k] + qi0_s[k]
        t0_s[k] = t0[k] + sten[k] * dt / cp
        k += 1

    m = 0
    while m < ncnst:
        offset = m * mkx
        k = 0
        while k < mkx:
            idx = k + offset
            tr0_s[idx] = tr0[idx] + trten[idx] * dt
            k += 1
        m += 1


@export
def uwshcu_iter_env_restore_state_shell_codon(
    mkx: int,
    qv0_s_p: cobj,
    ql0_s_p: cobj,
    qi0_s_p: cobj,
    s0_s_p: cobj,
    t0_s_p: cobj,
    qv0_p: cobj,
    ql0_p: cobj,
    qi0_p: cobj,
    s0_p: cobj,
    t0_p: cobj,
):
    qv0_s = Ptr[float](qv0_s_p)
    ql0_s = Ptr[float](ql0_s_p)
    qi0_s = Ptr[float](qi0_s_p)
    s0_s = Ptr[float](s0_s_p)
    t0_s = Ptr[float](t0_s_p)
    qv0 = Ptr[float](qv0_p)
    ql0 = Ptr[float](ql0_p)
    qi0 = Ptr[float](qi0_p)
    s0 = Ptr[float](s0_p)
    t0 = Ptr[float](t0_p)

    k = 0
    while k < mkx:
        qv0[k] = qv0_s[k]
        ql0[k] = ql0_s[k]
        qi0[k] = qi0_s[k]
        s0[k] = s0_s[k]
        t0[k] = t0_s[k]
        k += 1


@export
def uwshcu_iter_env_restore_thermo_slope_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    xlv: float,
    xls: float,
    cp: float,
    zvir: float,
    qv0_s_p: cobj,
    ql0_s_p: cobj,
    qi0_s_p: cobj,
    s0_s_p: cobj,
    t0_s_p: cobj,
    qv0_p: cobj,
    ql0_p: cobj,
    qi0_p: cobj,
    s0_p: cobj,
    t0_p: cobj,
    exn0_p: cobj,
    tr0_p: cobj,
    wtrc_iatype_p: cobj,
    p0_p: cobj,
    u0_p: cobj,
    v0_p: cobj,
    qt0_p: cobj,
    thl0_p: cobj,
    thvl0_p: cobj,
    wt0_p: cobj,
    ssthl0_p: cobj,
    ssqt0_p: cobj,
    ssu0_p: cobj,
    ssv0_p: cobj,
    sstr0_p: cobj,
    sswt0_p: cobj,
):
    uwshcu_iter_env_restore_state_shell_codon(
        mkx, qv0_s_p, ql0_s_p, qi0_s_p, s0_s_p, t0_s_p,
        qv0_p, ql0_p, qi0_p, s0_p, t0_p
    )
    uwshcu_column_thermo_slope_shell_codon(
        mkx, ncnst, wtrc_nwset, xlv, xls, cp, zvir, qv0_p, ql0_p,
        qi0_p, t0_p, exn0_p, tr0_p, wtrc_iatype_p, p0_p, u0_p, v0_p,
        qt0_p, thl0_p, thvl0_p, wt0_p, ssthl0_p, ssqt0_p, ssu0_p,
        ssv0_p, sstr0_p, sswt0_p
    )


@export
def uwshcu_iter_save_main_arrays_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    umf_p: cobj,
    qvten_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    sten_p: cobj,
    uten_p: cobj,
    vten_p: cobj,
    qrten_p: cobj,
    qsten_p: cobj,
    evapc_p: cobj,
    cufrc_p: cobj,
    slflx_p: cobj,
    qtflx_p: cobj,
    qcu_p: cobj,
    qlu_p: cobj,
    qiu_p: cobj,
    qc_p: cobj,
    trten_p: cobj,
    wtqc_liq_p: cobj,
    wtqc_ice_p: cobj,
    wtprec_p: cobj,
    wtsnow_p: cobj,
    umf_s_p: cobj,
    qvten_s_p: cobj,
    qlten_s_p: cobj,
    qiten_s_p: cobj,
    sten_s_p: cobj,
    uten_s_p: cobj,
    vten_s_p: cobj,
    qrten_s_p: cobj,
    qsten_s_p: cobj,
    evapc_s_p: cobj,
    cufrc_s_p: cobj,
    slflx_s_p: cobj,
    qtflx_s_p: cobj,
    qcu_s_p: cobj,
    qlu_s_p: cobj,
    qiu_s_p: cobj,
    qc_s_p: cobj,
    trten_s_p: cobj,
    wtqc_liq_s_p: cobj,
    wtqc_ice_s_p: cobj,
    wtprec_s_p: cobj,
    wtsnow_s_p: cobj,
):
    umf = Ptr[float](umf_p)
    qvten = Ptr[float](qvten_p)
    qlten = Ptr[float](qlten_p)
    qiten = Ptr[float](qiten_p)
    sten = Ptr[float](sten_p)
    uten = Ptr[float](uten_p)
    vten = Ptr[float](vten_p)
    qrten = Ptr[float](qrten_p)
    qsten = Ptr[float](qsten_p)
    evapc = Ptr[float](evapc_p)
    cufrc = Ptr[float](cufrc_p)
    slflx = Ptr[float](slflx_p)
    qtflx = Ptr[float](qtflx_p)
    qcu = Ptr[float](qcu_p)
    qlu = Ptr[float](qlu_p)
    qiu = Ptr[float](qiu_p)
    qc = Ptr[float](qc_p)
    trten = Ptr[float](trten_p)
    wtqc_liq = Ptr[float](wtqc_liq_p)
    wtqc_ice = Ptr[float](wtqc_ice_p)
    wtprec = Ptr[float](wtprec_p)
    wtsnow = Ptr[float](wtsnow_p)
    umf_s = Ptr[float](umf_s_p)
    qvten_s = Ptr[float](qvten_s_p)
    qlten_s = Ptr[float](qlten_s_p)
    qiten_s = Ptr[float](qiten_s_p)
    sten_s = Ptr[float](sten_s_p)
    uten_s = Ptr[float](uten_s_p)
    vten_s = Ptr[float](vten_s_p)
    qrten_s = Ptr[float](qrten_s_p)
    qsten_s = Ptr[float](qsten_s_p)
    evapc_s = Ptr[float](evapc_s_p)
    cufrc_s = Ptr[float](cufrc_s_p)
    slflx_s = Ptr[float](slflx_s_p)
    qtflx_s = Ptr[float](qtflx_s_p)
    qcu_s = Ptr[float](qcu_s_p)
    qlu_s = Ptr[float](qlu_s_p)
    qiu_s = Ptr[float](qiu_s_p)
    qc_s = Ptr[float](qc_s_p)
    trten_s = Ptr[float](trten_s_p)
    wtqc_liq_s = Ptr[float](wtqc_liq_s_p)
    wtqc_ice_s = Ptr[float](wtqc_ice_s_p)
    wtprec_s = Ptr[float](wtprec_s_p)
    wtsnow_s = Ptr[float](wtsnow_s_p)

    k = 0
    while k <= mkx:
        umf_s[k] = umf[k]
        slflx_s[k] = slflx[k]
        qtflx_s[k] = qtflx[k]
        k += 1

    k = 0
    while k < mkx:
        qvten_s[k] = qvten[k]
        qlten_s[k] = qlten[k]
        qiten_s[k] = qiten[k]
        sten_s[k] = sten[k]
        uten_s[k] = uten[k]
        vten_s[k] = vten[k]
        qrten_s[k] = qrten[k]
        qsten_s[k] = qsten[k]
        evapc_s[k] = evapc[k]
        cufrc_s[k] = cufrc[k]
        qcu_s[k] = qcu[k]
        qlu_s[k] = qlu[k]
        qiu_s[k] = qiu[k]
        qc_s[k] = qc[k]
        k += 1

    m = 0
    while m < ncnst:
        offset = m * mkx
        k = 0
        while k < mkx:
            idx = k + offset
            trten_s[idx] = trten[idx]
            k += 1
        m += 1

    m = 0
    while m < wtrc_nwset:
        offset = m * mkx
        k = 0
        while k < mkx:
            idx = k + offset
            wtqc_liq_s[idx] = wtqc_liq[idx]
            wtqc_ice_s[idx] = wtqc_ice[idx]
            k += 1
        wtprec_s[m] = wtprec[m]
        wtsnow_s[m] = wtsnow[m]
        m += 1


@export
def uwshcu_iter_save_diag_arrays_shell_codon(
    mkx: int,
    ncnst: int,
    fer_p: cobj,
    fdr_p: cobj,
    qtten_p: cobj,
    slten_p: cobj,
    ufrc_p: cobj,
    uflx_p: cobj,
    vflx_p: cobj,
    wu_p: cobj,
    qtu_p: cobj,
    thlu_p: cobj,
    thvu_p: cobj,
    uu_p: cobj,
    vu_p: cobj,
    qtu_emf_p: cobj,
    thlu_emf_p: cobj,
    uu_emf_p: cobj,
    vu_emf_p: cobj,
    uemf_p: cobj,
    dwten_p: cobj,
    diten_p: cobj,
    flxrain_p: cobj,
    flxsnow_p: cobj,
    ntraprd_p: cobj,
    ntsnprd_p: cobj,
    excessu_p: cobj,
    excess0_p: cobj,
    xc_p: cobj,
    aquad_p: cobj,
    bquad_p: cobj,
    cquad_p: cobj,
    bogbot_p: cobj,
    bogtop_p: cobj,
    trflx_p: cobj,
    tru_p: cobj,
    tru_emf_p: cobj,
    fer_s_p: cobj,
    fdr_s_p: cobj,
    qtten_s_p: cobj,
    slten_s_p: cobj,
    ufrc_s_p: cobj,
    uflx_s_p: cobj,
    vflx_s_p: cobj,
    wu_s_p: cobj,
    qtu_s_p: cobj,
    thlu_s_p: cobj,
    thvu_s_p: cobj,
    uu_s_p: cobj,
    vu_s_p: cobj,
    qtu_emf_s_p: cobj,
    thlu_emf_s_p: cobj,
    uu_emf_s_p: cobj,
    vu_emf_s_p: cobj,
    uemf_s_p: cobj,
    dwten_s_p: cobj,
    diten_s_p: cobj,
    flxrain_s_p: cobj,
    flxsnow_s_p: cobj,
    ntraprd_s_p: cobj,
    ntsnprd_s_p: cobj,
    excessu_s_p: cobj,
    excess0_s_p: cobj,
    xc_s_p: cobj,
    aquad_s_p: cobj,
    bquad_s_p: cobj,
    cquad_s_p: cobj,
    bogbot_s_p: cobj,
    bogtop_s_p: cobj,
    trflx_s_p: cobj,
    tru_s_p: cobj,
    tru_emf_s_p: cobj,
):
    fer = Ptr[float](fer_p)
    fdr = Ptr[float](fdr_p)
    qtten = Ptr[float](qtten_p)
    slten = Ptr[float](slten_p)
    ufrc = Ptr[float](ufrc_p)
    uflx = Ptr[float](uflx_p)
    vflx = Ptr[float](vflx_p)
    wu = Ptr[float](wu_p)
    qtu = Ptr[float](qtu_p)
    thlu = Ptr[float](thlu_p)
    thvu = Ptr[float](thvu_p)
    uu = Ptr[float](uu_p)
    vu = Ptr[float](vu_p)
    qtu_emf = Ptr[float](qtu_emf_p)
    thlu_emf = Ptr[float](thlu_emf_p)
    uu_emf = Ptr[float](uu_emf_p)
    vu_emf = Ptr[float](vu_emf_p)
    uemf = Ptr[float](uemf_p)
    dwten = Ptr[float](dwten_p)
    diten = Ptr[float](diten_p)
    flxrain = Ptr[float](flxrain_p)
    flxsnow = Ptr[float](flxsnow_p)
    ntraprd = Ptr[float](ntraprd_p)
    ntsnprd = Ptr[float](ntsnprd_p)
    excessu = Ptr[float](excessu_p)
    excess0 = Ptr[float](excess0_p)
    xc = Ptr[float](xc_p)
    aquad = Ptr[float](aquad_p)
    bquad = Ptr[float](bquad_p)
    cquad = Ptr[float](cquad_p)
    bogbot = Ptr[float](bogbot_p)
    bogtop = Ptr[float](bogtop_p)
    trflx = Ptr[float](trflx_p)
    tru = Ptr[float](tru_p)
    tru_emf = Ptr[float](tru_emf_p)
    fer_s = Ptr[float](fer_s_p)
    fdr_s = Ptr[float](fdr_s_p)
    qtten_s = Ptr[float](qtten_s_p)
    slten_s = Ptr[float](slten_s_p)
    ufrc_s = Ptr[float](ufrc_s_p)
    uflx_s = Ptr[float](uflx_s_p)
    vflx_s = Ptr[float](vflx_s_p)
    wu_s = Ptr[float](wu_s_p)
    qtu_s = Ptr[float](qtu_s_p)
    thlu_s = Ptr[float](thlu_s_p)
    thvu_s = Ptr[float](thvu_s_p)
    uu_s = Ptr[float](uu_s_p)
    vu_s = Ptr[float](vu_s_p)
    qtu_emf_s = Ptr[float](qtu_emf_s_p)
    thlu_emf_s = Ptr[float](thlu_emf_s_p)
    uu_emf_s = Ptr[float](uu_emf_s_p)
    vu_emf_s = Ptr[float](vu_emf_s_p)
    uemf_s = Ptr[float](uemf_s_p)
    dwten_s = Ptr[float](dwten_s_p)
    diten_s = Ptr[float](diten_s_p)
    flxrain_s = Ptr[float](flxrain_s_p)
    flxsnow_s = Ptr[float](flxsnow_s_p)
    ntraprd_s = Ptr[float](ntraprd_s_p)
    ntsnprd_s = Ptr[float](ntsnprd_s_p)
    excessu_s = Ptr[float](excessu_s_p)
    excess0_s = Ptr[float](excess0_s_p)
    xc_s = Ptr[float](xc_s_p)
    aquad_s = Ptr[float](aquad_s_p)
    bquad_s = Ptr[float](bquad_s_p)
    cquad_s = Ptr[float](cquad_s_p)
    bogbot_s = Ptr[float](bogbot_s_p)
    bogtop_s = Ptr[float](bogtop_s_p)
    trflx_s = Ptr[float](trflx_s_p)
    tru_s = Ptr[float](tru_s_p)
    tru_emf_s = Ptr[float](tru_emf_s_p)

    k = 0
    while k < mkx:
        fer_s[k] = fer[k]
        fdr_s[k] = fdr[k]
        qtten_s[k] = qtten[k]
        slten_s[k] = slten[k]
        dwten_s[k] = dwten[k]
        diten_s[k] = diten[k]
        ntraprd_s[k] = ntraprd[k]
        ntsnprd_s[k] = ntsnprd[k]
        excessu_s[k] = excessu[k]
        excess0_s[k] = excess0[k]
        xc_s[k] = xc[k]
        aquad_s[k] = aquad[k]
        bquad_s[k] = bquad[k]
        cquad_s[k] = cquad[k]
        bogbot_s[k] = bogbot[k]
        bogtop_s[k] = bogtop[k]
        k += 1

    k = 0
    while k <= mkx:
        ufrc_s[k] = ufrc[k]
        uflx_s[k] = uflx[k]
        vflx_s[k] = vflx[k]
        wu_s[k] = wu[k]
        qtu_s[k] = qtu[k]
        thlu_s[k] = thlu[k]
        thvu_s[k] = thvu[k]
        uu_s[k] = uu[k]
        vu_s[k] = vu[k]
        qtu_emf_s[k] = qtu_emf[k]
        thlu_emf_s[k] = thlu_emf[k]
        uu_emf_s[k] = uu_emf[k]
        vu_emf_s[k] = vu_emf[k]
        uemf_s[k] = uemf[k]
        flxrain_s[k] = flxrain[k]
        flxsnow_s[k] = flxsnow[k]
        k += 1

    m = 0
    while m < ncnst:
        offset = m * (mkx + 1)
        k = 0
        while k <= mkx:
            idx = k + offset
            trflx_s[idx] = trflx[idx]
            tru_s[idx] = tru[idx]
            tru_emf_s[idx] = tru_emf[idx]
            k += 1
        m += 1


@export
def uwshcu_iter_save_all_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    kinv: int,
    kbup: int,
    kpen: int,
    dt: float,
    cp: float,
    ppen: float,
    precip_v: float,
    snow_v: float,
    cush_v: float,
    cin_v: float,
    cinlcl_v: float,
    cbmf_v: float,
    rliq_v: float,
    cnt_v: float,
    cnb_v: float,
    ufrcinvbase_v: float,
    ufrclcl_v: float,
    winvbase_v: float,
    wlcl_v: float,
    plcl_v: float,
    plfc_v: float,
    qtsrc_v: float,
    thlsrc_v: float,
    thvlsrc_v: float,
    cbmflimit_v: float,
    tkeavg_v: float,
    rcwp_v: float,
    rlwp_v: float,
    riwp_v: float,
    qv0_p: cobj,
    ql0_p: cobj,
    qi0_p: cobj,
    s0_p: cobj,
    u0_p: cobj,
    v0_p: cobj,
    t0_p: cobj,
    qvten_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    sten_p: cobj,
    uten_p: cobj,
    vten_p: cobj,
    tr0_p: cobj,
    trten_p: cobj,
    qv0_s_p: cobj,
    ql0_s_p: cobj,
    qi0_s_p: cobj,
    s0_s_p: cobj,
    u0_s_p: cobj,
    v0_s_p: cobj,
    qt0_s_p: cobj,
    t0_s_p: cobj,
    tr0_s_p: cobj,
    umf_p: cobj,
    qrten_p: cobj,
    qsten_p: cobj,
    evapc_p: cobj,
    cufrc_p: cobj,
    slflx_p: cobj,
    qtflx_p: cobj,
    qcu_p: cobj,
    qlu_p: cobj,
    qiu_p: cobj,
    qc_p: cobj,
    wtqc_liq_p: cobj,
    wtqc_ice_p: cobj,
    wtprec_p: cobj,
    wtsnow_p: cobj,
    umf_s_p: cobj,
    qvten_s_p: cobj,
    qlten_s_p: cobj,
    qiten_s_p: cobj,
    sten_s_p: cobj,
    uten_s_p: cobj,
    vten_s_p: cobj,
    qrten_s_p: cobj,
    qsten_s_p: cobj,
    evapc_s_p: cobj,
    cufrc_s_p: cobj,
    slflx_s_p: cobj,
    qtflx_s_p: cobj,
    qcu_s_p: cobj,
    qlu_s_p: cobj,
    qiu_s_p: cobj,
    qc_s_p: cobj,
    trten_s_p: cobj,
    wtqc_liq_s_p: cobj,
    wtqc_ice_s_p: cobj,
    wtprec_s_p: cobj,
    wtsnow_s_p: cobj,
    fer_p: cobj,
    fdr_p: cobj,
    qtten_p: cobj,
    slten_p: cobj,
    ufrc_p: cobj,
    uflx_p: cobj,
    vflx_p: cobj,
    wu_p: cobj,
    qtu_p: cobj,
    thlu_p: cobj,
    thvu_p: cobj,
    uu_p: cobj,
    vu_p: cobj,
    qtu_emf_p: cobj,
    thlu_emf_p: cobj,
    uu_emf_p: cobj,
    vu_emf_p: cobj,
    uemf_p: cobj,
    dwten_p: cobj,
    diten_p: cobj,
    flxrain_p: cobj,
    flxsnow_p: cobj,
    ntraprd_p: cobj,
    ntsnprd_p: cobj,
    excessu_p: cobj,
    excess0_p: cobj,
    xc_p: cobj,
    aquad_p: cobj,
    bquad_p: cobj,
    cquad_p: cobj,
    bogbot_p: cobj,
    bogtop_p: cobj,
    trflx_p: cobj,
    tru_p: cobj,
    tru_emf_p: cobj,
    fer_s_p: cobj,
    fdr_s_p: cobj,
    qtten_s_p: cobj,
    slten_s_p: cobj,
    ufrc_s_p: cobj,
    uflx_s_p: cobj,
    vflx_s_p: cobj,
    wu_s_p: cobj,
    qtu_s_p: cobj,
    thlu_s_p: cobj,
    thvu_s_p: cobj,
    uu_s_p: cobj,
    vu_s_p: cobj,
    qtu_emf_s_p: cobj,
    thlu_emf_s_p: cobj,
    uu_emf_s_p: cobj,
    vu_emf_s_p: cobj,
    uemf_s_p: cobj,
    dwten_s_p: cobj,
    diten_s_p: cobj,
    flxrain_s_p: cobj,
    flxsnow_s_p: cobj,
    ntraprd_s_p: cobj,
    ntsnprd_s_p: cobj,
    excessu_s_p: cobj,
    excess0_s_p: cobj,
    xc_s_p: cobj,
    aquad_s_p: cobj,
    bquad_s_p: cobj,
    cquad_s_p: cobj,
    bogbot_s_p: cobj,
    bogtop_s_p: cobj,
    trflx_s_p: cobj,
    tru_s_p: cobj,
    tru_emf_s_p: cobj,
    ps0_p: cobj,
    zs0_p: cobj,
    emf_p: cobj,
    precip_s_p: cobj,
    snow_s_p: cobj,
    cush_s_p: cobj,
    cin_s_p: cobj,
    cinlcl_s_p: cobj,
    cbmf_s_p: cobj,
    rliq_s_p: cobj,
    cnt_s_p: cobj,
    cnb_s_p: cobj,
    ufrcinvbase_s_p: cobj,
    ufrclcl_s_p: cobj,
    winvbase_s_p: cobj,
    wlcl_s_p: cobj,
    plcl_s_p: cobj,
    pinv_s_p: cobj,
    plfc_s_p: cobj,
    pbup_s_p: cobj,
    ppen_s_p: cobj,
    qtsrc_s_p: cobj,
    thlsrc_s_p: cobj,
    thvlsrc_s_p: cobj,
    emfkbup_s_p: cobj,
    cbmflimit_s_p: cobj,
    tkeavg_s_p: cobj,
    zinv_s_p: cobj,
    rcwp_s_p: cobj,
    rlwp_s_p: cobj,
    riwp_s_p: cobj,
    xlv: float,
    xls: float,
    zvir: float,
    exn0_p: cobj,
    wtrc_iatype_p: cobj,
    p0_p: cobj,
    qt0_p: cobj,
    thl0_p: cobj,
    thvl0_p: cobj,
    wt0_p: cobj,
    ssthl0_p: cobj,
    ssqt0_p: cobj,
    ssu0_p: cobj,
    ssv0_p: cobj,
    sstr0_p: cobj,
    sswt0_p: cobj,
):
    uwshcu_iter_save_env_shell_codon(
        mkx, ncnst, dt, cp, qv0_p, ql0_p, qi0_p, s0_p, u0_p, v0_p, t0_p,
        qvten_p, qlten_p, qiten_p, sten_p, uten_p, vten_p, tr0_p, trten_p,
        qv0_s_p, ql0_s_p, qi0_s_p, s0_s_p, u0_s_p, v0_s_p, qt0_s_p,
        t0_s_p, tr0_s_p,
    )
    uwshcu_iter_save_main_arrays_shell_codon(
        mkx, ncnst, wtrc_nwset, umf_p, qvten_p, qlten_p, qiten_p, sten_p,
        uten_p, vten_p, qrten_p, qsten_p, evapc_p, cufrc_p, slflx_p,
        qtflx_p, qcu_p, qlu_p, qiu_p, qc_p, trten_p, wtqc_liq_p,
        wtqc_ice_p, wtprec_p, wtsnow_p, umf_s_p, qvten_s_p, qlten_s_p,
        qiten_s_p, sten_s_p, uten_s_p, vten_s_p, qrten_s_p, qsten_s_p,
        evapc_s_p, cufrc_s_p, slflx_s_p, qtflx_s_p, qcu_s_p, qlu_s_p,
        qiu_s_p, qc_s_p, trten_s_p, wtqc_liq_s_p, wtqc_ice_s_p,
        wtprec_s_p, wtsnow_s_p,
    )

    ps0 = Ptr[float](ps0_p)
    zs0 = Ptr[float](zs0_p)
    emf = Ptr[float](emf_p)
    Ptr[float](precip_s_p)[0] = precip_v
    Ptr[float](snow_s_p)[0] = snow_v
    Ptr[float](cush_s_p)[0] = cush_v
    Ptr[float](cin_s_p)[0] = cin_v
    Ptr[float](cinlcl_s_p)[0] = cinlcl_v
    Ptr[float](cbmf_s_p)[0] = cbmf_v
    Ptr[float](rliq_s_p)[0] = rliq_v
    Ptr[float](cnt_s_p)[0] = cnt_v
    Ptr[float](cnb_s_p)[0] = cnb_v
    Ptr[float](ufrcinvbase_s_p)[0] = ufrcinvbase_v
    Ptr[float](ufrclcl_s_p)[0] = ufrclcl_v
    Ptr[float](winvbase_s_p)[0] = winvbase_v
    Ptr[float](wlcl_s_p)[0] = wlcl_v
    Ptr[float](plcl_s_p)[0] = plcl_v
    Ptr[float](pinv_s_p)[0] = ps0[kinv - 1]
    Ptr[float](plfc_s_p)[0] = plfc_v
    Ptr[float](pbup_s_p)[0] = ps0[kbup]
    Ptr[float](ppen_s_p)[0] = ps0[kpen - 1] + ppen
    Ptr[float](qtsrc_s_p)[0] = qtsrc_v
    Ptr[float](thlsrc_s_p)[0] = thlsrc_v
    Ptr[float](thvlsrc_s_p)[0] = thvlsrc_v
    Ptr[float](emfkbup_s_p)[0] = emf[kbup]
    Ptr[float](cbmflimit_s_p)[0] = cbmflimit_v
    Ptr[float](tkeavg_s_p)[0] = tkeavg_v
    Ptr[float](zinv_s_p)[0] = zs0[kinv - 1]
    Ptr[float](rcwp_s_p)[0] = rcwp_v
    Ptr[float](rlwp_s_p)[0] = rlwp_v
    Ptr[float](riwp_s_p)[0] = riwp_v

    uwshcu_iter_save_diag_arrays_shell_codon(
        mkx, ncnst, fer_p, fdr_p, qtten_p, slten_p, ufrc_p, uflx_p, vflx_p,
        wu_p, qtu_p, thlu_p, thvu_p, uu_p, vu_p, qtu_emf_p, thlu_emf_p,
        uu_emf_p, vu_emf_p, uemf_p, dwten_p, diten_p, flxrain_p, flxsnow_p,
        ntraprd_p, ntsnprd_p, excessu_p, excess0_p, xc_p, aquad_p, bquad_p,
        cquad_p, bogbot_p, bogtop_p, trflx_p, tru_p, tru_emf_p, fer_s_p,
        fdr_s_p, qtten_s_p, slten_s_p, ufrc_s_p, uflx_s_p, vflx_s_p, wu_s_p,
        qtu_s_p, thlu_s_p, thvu_s_p, uu_s_p, vu_s_p, qtu_emf_s_p,
        thlu_emf_s_p, uu_emf_s_p, vu_emf_s_p, uemf_s_p, dwten_s_p,
        diten_s_p, flxrain_s_p, flxsnow_s_p, ntraprd_s_p, ntsnprd_s_p,
        excessu_s_p, excess0_s_p, xc_s_p, aquad_s_p, bquad_s_p, cquad_s_p,
        bogbot_s_p, bogtop_s_p, trflx_s_p, tru_s_p, tru_emf_s_p,
    )

    uwshcu_iter_env_restore_thermo_slope_shell_codon(
        mkx, ncnst, 0, xlv, xls, cp, zvir, qv0_s_p, ql0_s_p, qi0_s_p,
        s0_s_p, t0_s_p, qv0_p, ql0_p, qi0_p, s0_p, t0_p, exn0_p, tr0_p,
        wtrc_iatype_p, p0_p, u0_p, v0_p, qt0_p, thl0_p, thvl0_p, wt0_p,
        ssthl0_p, ssqt0_p, ssu0_p, ssv0_p, sstr0_p, sswt0_p,
    )


@export
def uwshcu_delcin_env_restore_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    qv0_o_p: cobj,
    ql0_o_p: cobj,
    qi0_o_p: cobj,
    t0_o_p: cobj,
    s0_o_p: cobj,
    u0_o_p: cobj,
    v0_o_p: cobj,
    qt0_o_p: cobj,
    thl0_o_p: cobj,
    thvl0_o_p: cobj,
    ssthl0_o_p: cobj,
    ssqt0_o_p: cobj,
    thv0bot_o_p: cobj,
    thv0top_o_p: cobj,
    thvl0bot_o_p: cobj,
    thvl0top_o_p: cobj,
    ssu0_o_p: cobj,
    ssv0_o_p: cobj,
    tr0_o_p: cobj,
    sstr0_o_p: cobj,
    sswt0_o_p: cobj,
    qv0_p: cobj,
    ql0_p: cobj,
    qi0_p: cobj,
    t0_p: cobj,
    s0_p: cobj,
    u0_p: cobj,
    v0_p: cobj,
    qt0_p: cobj,
    thl0_p: cobj,
    thvl0_p: cobj,
    ssthl0_p: cobj,
    ssqt0_p: cobj,
    thv0bot_p: cobj,
    thv0top_p: cobj,
    thvl0bot_p: cobj,
    thvl0top_p: cobj,
    ssu0_p: cobj,
    ssv0_p: cobj,
    tr0_p: cobj,
    sstr0_p: cobj,
    sswt0_p: cobj,
):
    qv0_o = Ptr[float](qv0_o_p)
    ql0_o = Ptr[float](ql0_o_p)
    qi0_o = Ptr[float](qi0_o_p)
    t0_o = Ptr[float](t0_o_p)
    s0_o = Ptr[float](s0_o_p)
    u0_o = Ptr[float](u0_o_p)
    v0_o = Ptr[float](v0_o_p)
    qt0_o = Ptr[float](qt0_o_p)
    thl0_o = Ptr[float](thl0_o_p)
    thvl0_o = Ptr[float](thvl0_o_p)
    ssthl0_o = Ptr[float](ssthl0_o_p)
    ssqt0_o = Ptr[float](ssqt0_o_p)
    thv0bot_o = Ptr[float](thv0bot_o_p)
    thv0top_o = Ptr[float](thv0top_o_p)
    thvl0bot_o = Ptr[float](thvl0bot_o_p)
    thvl0top_o = Ptr[float](thvl0top_o_p)
    ssu0_o = Ptr[float](ssu0_o_p)
    ssv0_o = Ptr[float](ssv0_o_p)
    tr0_o = Ptr[float](tr0_o_p)
    sstr0_o = Ptr[float](sstr0_o_p)
    sswt0_o = Ptr[float](sswt0_o_p)
    qv0 = Ptr[float](qv0_p)
    ql0 = Ptr[float](ql0_p)
    qi0 = Ptr[float](qi0_p)
    t0 = Ptr[float](t0_p)
    s0 = Ptr[float](s0_p)
    u0 = Ptr[float](u0_p)
    v0 = Ptr[float](v0_p)
    qt0 = Ptr[float](qt0_p)
    thl0 = Ptr[float](thl0_p)
    thvl0 = Ptr[float](thvl0_p)
    ssthl0 = Ptr[float](ssthl0_p)
    ssqt0 = Ptr[float](ssqt0_p)
    thv0bot = Ptr[float](thv0bot_p)
    thv0top = Ptr[float](thv0top_p)
    thvl0bot = Ptr[float](thvl0bot_p)
    thvl0top = Ptr[float](thvl0top_p)
    ssu0 = Ptr[float](ssu0_p)
    ssv0 = Ptr[float](ssv0_p)
    tr0 = Ptr[float](tr0_p)
    sstr0 = Ptr[float](sstr0_p)
    sswt0 = Ptr[float](sswt0_p)

    k = 0
    while k < mkx:
        qv0[k] = qv0_o[k]
        ql0[k] = ql0_o[k]
        qi0[k] = qi0_o[k]
        t0[k] = t0_o[k]
        s0[k] = s0_o[k]
        u0[k] = u0_o[k]
        v0[k] = v0_o[k]
        qt0[k] = qt0_o[k]
        thl0[k] = thl0_o[k]
        thvl0[k] = thvl0_o[k]
        ssthl0[k] = ssthl0_o[k]
        ssqt0[k] = ssqt0_o[k]
        thv0bot[k] = thv0bot_o[k]
        thv0top[k] = thv0top_o[k]
        thvl0bot[k] = thvl0bot_o[k]
        thvl0top[k] = thvl0top_o[k]
        ssu0[k] = ssu0_o[k]
        ssv0[k] = ssv0_o[k]
        k += 1

    m = 0
    while m < ncnst:
        offset = m * mkx
        k = 0
        while k < mkx:
            idx = k + offset
            tr0[idx] = tr0_o[idx]
            sstr0[idx] = sstr0_o[idx]
            k += 1
        m += 1

    m = 0
    while m < wtrc_nwset:
        offset = m * mkx
        k = 0
        while k < mkx:
            idx = k + offset
            sswt0[idx] = sswt0_o[idx]
            k += 1
        m += 1


@export
def uwshcu_delcin_workspace_reset_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    umf_p: cobj,
    emf_p: cobj,
    slflx_p: cobj,
    qtflx_p: cobj,
    uflx_p: cobj,
    vflx_p: cobj,
    qvten_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    sten_p: cobj,
    uten_p: cobj,
    vten_p: cobj,
    qrten_p: cobj,
    qsten_p: cobj,
    dwten_p: cobj,
    diten_p: cobj,
    evapc_p: cobj,
    cufrc_p: cobj,
    qcu_p: cobj,
    qlu_p: cobj,
    qiu_p: cobj,
    fer_p: cobj,
    fdr_p: cobj,
    qc_p: cobj,
    qc_l_p: cobj,
    qc_i_p: cobj,
    qtten_p: cobj,
    slten_p: cobj,
    ufrc_p: cobj,
    thlu_p: cobj,
    qtu_p: cobj,
    uu_p: cobj,
    vu_p: cobj,
    wu_p: cobj,
    thvu_p: cobj,
    thlu_emf_p: cobj,
    qtu_emf_p: cobj,
    uu_emf_p: cobj,
    vu_emf_p: cobj,
    trflx_p: cobj,
    trten_p: cobj,
    tru_p: cobj,
    tru_emf_p: cobj,
    wtdwten_p: cobj,
    wtditen_p: cobj,
    wtrpten_p: cobj,
    wtspten_p: cobj,
    wtqc_liq_p: cobj,
    wtqc_ice_p: cobj,
    wtu_p: cobj,
    wtu_emf_p: cobj,
    wtflx_p: cobj,
    wttotten_p: cobj,
    wtprec_p: cobj,
    wtsnow_p: cobj,
    excessu_p: cobj,
    excess0_p: cobj,
    xc_p: cobj,
    aquad_p: cobj,
    bquad_p: cobj,
    cquad_p: cobj,
    bogbot_p: cobj,
    bogtop_p: cobj,
):
    umf = Ptr[float](umf_p)
    emf = Ptr[float](emf_p)
    slflx = Ptr[float](slflx_p)
    qtflx = Ptr[float](qtflx_p)
    uflx = Ptr[float](uflx_p)
    vflx = Ptr[float](vflx_p)
    qvten = Ptr[float](qvten_p)
    qlten = Ptr[float](qlten_p)
    qiten = Ptr[float](qiten_p)
    sten = Ptr[float](sten_p)
    uten = Ptr[float](uten_p)
    vten = Ptr[float](vten_p)
    qrten = Ptr[float](qrten_p)
    qsten = Ptr[float](qsten_p)
    dwten = Ptr[float](dwten_p)
    diten = Ptr[float](diten_p)
    evapc = Ptr[float](evapc_p)
    cufrc = Ptr[float](cufrc_p)
    qcu = Ptr[float](qcu_p)
    qlu = Ptr[float](qlu_p)
    qiu = Ptr[float](qiu_p)
    fer = Ptr[float](fer_p)
    fdr = Ptr[float](fdr_p)
    qc = Ptr[float](qc_p)
    qc_l = Ptr[float](qc_l_p)
    qc_i = Ptr[float](qc_i_p)
    qtten = Ptr[float](qtten_p)
    slten = Ptr[float](slten_p)
    ufrc = Ptr[float](ufrc_p)
    thlu = Ptr[float](thlu_p)
    qtu = Ptr[float](qtu_p)
    uu = Ptr[float](uu_p)
    vu = Ptr[float](vu_p)
    wu = Ptr[float](wu_p)
    thvu = Ptr[float](thvu_p)
    thlu_emf = Ptr[float](thlu_emf_p)
    qtu_emf = Ptr[float](qtu_emf_p)
    uu_emf = Ptr[float](uu_emf_p)
    vu_emf = Ptr[float](vu_emf_p)
    trflx = Ptr[float](trflx_p)
    trten = Ptr[float](trten_p)
    tru = Ptr[float](tru_p)
    tru_emf = Ptr[float](tru_emf_p)
    wtdwten = Ptr[float](wtdwten_p)
    wtditen = Ptr[float](wtditen_p)
    wtrpten = Ptr[float](wtrpten_p)
    wtspten = Ptr[float](wtspten_p)
    wtqc_liq = Ptr[float](wtqc_liq_p)
    wtqc_ice = Ptr[float](wtqc_ice_p)
    wtu = Ptr[float](wtu_p)
    wtu_emf = Ptr[float](wtu_emf_p)
    wtflx = Ptr[float](wtflx_p)
    wttotten = Ptr[float](wttotten_p)
    wtprec = Ptr[float](wtprec_p)
    wtsnow = Ptr[float](wtsnow_p)
    excessu = Ptr[float](excessu_p)
    excess0 = Ptr[float](excess0_p)
    xc = Ptr[float](xc_p)
    aquad = Ptr[float](aquad_p)
    bquad = Ptr[float](bquad_p)
    cquad = Ptr[float](cquad_p)
    bogbot = Ptr[float](bogbot_p)
    bogtop = Ptr[float](bogtop_p)

    k = 0
    while k <= mkx:
        umf[k] = 0.0
        emf[k] = 0.0
        slflx[k] = 0.0
        qtflx[k] = 0.0
        uflx[k] = 0.0
        vflx[k] = 0.0
        ufrc[k] = 0.0
        thlu[k] = 0.0
        qtu[k] = 0.0
        uu[k] = 0.0
        vu[k] = 0.0
        wu[k] = 0.0
        thvu[k] = 0.0
        thlu_emf[k] = 0.0
        qtu_emf[k] = 0.0
        uu_emf[k] = 0.0
        vu_emf[k] = 0.0
        k += 1

    k = 0
    while k < mkx:
        qvten[k] = 0.0
        qlten[k] = 0.0
        qiten[k] = 0.0
        sten[k] = 0.0
        uten[k] = 0.0
        vten[k] = 0.0
        qrten[k] = 0.0
        qsten[k] = 0.0
        dwten[k] = 0.0
        diten[k] = 0.0
        evapc[k] = 0.0
        cufrc[k] = 0.0
        qcu[k] = 0.0
        qlu[k] = 0.0
        qiu[k] = 0.0
        fer[k] = 0.0
        fdr[k] = 0.0
        qc[k] = 0.0
        qc_l[k] = 0.0
        qc_i[k] = 0.0
        qtten[k] = 0.0
        slten[k] = 0.0
        excessu[k] = 0.0
        excess0[k] = 0.0
        xc[k] = 0.0
        aquad[k] = 0.0
        bquad[k] = 0.0
        cquad[k] = 0.0
        bogbot[k] = 0.0
        bogtop[k] = 0.0
        k += 1

    m = 0
    while m < ncnst:
        layer_offset = m * mkx
        iface_offset = m * (mkx + 1)
        k = 0
        while k < mkx:
            trten[k + layer_offset] = 0.0
            k += 1
        k = 0
        while k <= mkx:
            idx = k + iface_offset
            trflx[idx] = 0.0
            tru[idx] = 0.0
            tru_emf[idx] = 0.0
            k += 1
        m += 1

    m = 0
    while m < wtrc_nwset:
        layer_offset = m * mkx
        iface_offset = m * (mkx + 1)
        k = 0
        while k < mkx:
            idx = k + layer_offset
            wtdwten[idx] = 0.0
            wtditen[idx] = 0.0
            wtrpten[idx] = 0.0
            wtspten[idx] = 0.0
            wtqc_liq[idx] = 0.0
            wtqc_ice[idx] = 0.0
            wttotten[idx] = 0.0
            k += 1
        k = 0
        while k <= mkx:
            idx = k + iface_offset
            wtu[idx] = 0.0
            wtu_emf[idx] = 0.0
            wtflx[idx] = 0.0
            k += 1
        wtprec[m] = 0.0
        wtsnow[m] = 0.0
        m += 1


@export
def uwshcu_delcin_restore_reset_all_shell_codon(
    mkx: int,
    ncnst: int,
    wtrc_nwset: int,
    use_cincin: int,
    cin_i_v: float,
    cinlcl_i_v: float,
    alpha: float,
    del_cin: float,
    del_cinlcl: float,
    kinv_o_v: int,
    klcl_o_v: int,
    klfc_o_v: int,
    plcl_o_v: float,
    plfc_o_v: float,
    tkeavg_o_v: float,
    thvlmin_o_v: float,
    qtsrc_o_v: float,
    thvlsrc_o_v: float,
    thlsrc_o_v: float,
    usrc_o_v: float,
    vsrc_o_v: float,
    thv0lcl_o_v: float,
    cnt_reset_v: float,
    trsrc_o_p: cobj,
    cin_p: cobj,
    cinlcl_p: cobj,
    kinv_p: cobj,
    klcl_p: cobj,
    klfc_p: cobj,
    plcl_p: cobj,
    plfc_p: cobj,
    tkeavg_p: cobj,
    thvlmin_p: cobj,
    qtsrc_p: cobj,
    thvlsrc_p: cobj,
    thlsrc_p: cobj,
    usrc_p: cobj,
    vsrc_p: cobj,
    thv0lcl_p: cobj,
    trsrc_p: cobj,
    qv0_o_p: cobj,
    ql0_o_p: cobj,
    qi0_o_p: cobj,
    t0_o_p: cobj,
    s0_o_p: cobj,
    u0_o_p: cobj,
    v0_o_p: cobj,
    qt0_o_p: cobj,
    thl0_o_p: cobj,
    thvl0_o_p: cobj,
    ssthl0_o_p: cobj,
    ssqt0_o_p: cobj,
    thv0bot_o_p: cobj,
    thv0top_o_p: cobj,
    thvl0bot_o_p: cobj,
    thvl0top_o_p: cobj,
    ssu0_o_p: cobj,
    ssv0_o_p: cobj,
    tr0_o_p: cobj,
    sstr0_o_p: cobj,
    sswt0_o_p: cobj,
    qv0_p: cobj,
    ql0_p: cobj,
    qi0_p: cobj,
    t0_p: cobj,
    s0_p: cobj,
    u0_p: cobj,
    v0_p: cobj,
    qt0_p: cobj,
    thl0_p: cobj,
    thvl0_p: cobj,
    ssthl0_p: cobj,
    ssqt0_p: cobj,
    thv0bot_p: cobj,
    thv0top_p: cobj,
    thvl0bot_p: cobj,
    thvl0top_p: cobj,
    ssu0_p: cobj,
    ssv0_p: cobj,
    tr0_p: cobj,
    sstr0_p: cobj,
    sswt0_p: cobj,
    umf_p: cobj,
    emf_p: cobj,
    slflx_p: cobj,
    qtflx_p: cobj,
    uflx_p: cobj,
    vflx_p: cobj,
    qvten_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    sten_p: cobj,
    uten_p: cobj,
    vten_p: cobj,
    qrten_p: cobj,
    qsten_p: cobj,
    dwten_p: cobj,
    diten_p: cobj,
    evapc_p: cobj,
    cufrc_p: cobj,
    qcu_p: cobj,
    qlu_p: cobj,
    qiu_p: cobj,
    fer_p: cobj,
    fdr_p: cobj,
    qc_p: cobj,
    qc_l_p: cobj,
    qc_i_p: cobj,
    qtten_p: cobj,
    slten_p: cobj,
    ufrc_p: cobj,
    thlu_p: cobj,
    qtu_p: cobj,
    uu_p: cobj,
    vu_p: cobj,
    wu_p: cobj,
    thvu_p: cobj,
    thlu_emf_p: cobj,
    qtu_emf_p: cobj,
    uu_emf_p: cobj,
    vu_emf_p: cobj,
    trflx_p: cobj,
    trten_p: cobj,
    tru_p: cobj,
    tru_emf_p: cobj,
    wtdwten_p: cobj,
    wtditen_p: cobj,
    wtrpten_p: cobj,
    wtspten_p: cobj,
    wtqc_liq_p: cobj,
    wtqc_ice_p: cobj,
    wtu_p: cobj,
    wtu_emf_p: cobj,
    wtflx_p: cobj,
    wttotten_p: cobj,
    wtprec_p: cobj,
    wtsnow_p: cobj,
    excessu_p: cobj,
    excess0_p: cobj,
    xc_p: cobj,
    aquad_p: cobj,
    bquad_p: cobj,
    cquad_p: cobj,
    bogbot_p: cobj,
    bogtop_p: cobj,
    precip_p: cobj,
    snow_p: cobj,
    rliq_p: cobj,
    cbmf_p: cobj,
    cnt_p: cobj,
    cnb_p: cobj,
    ufrcinvbase_p: cobj,
    ufrclcl_p: cobj,
    winvbase_p: cobj,
    wlcl_p: cobj,
    emfkbup_p: cobj,
    cbmflimit_p: cobj,
):
    uwshcu_cin_state_restore_shell_codon(
        ncnst, use_cincin, cin_i_v, cinlcl_i_v, alpha, del_cin, del_cinlcl,
        kinv_o_v, klcl_o_v, klfc_o_v, plcl_o_v, plfc_o_v, tkeavg_o_v,
        thvlmin_o_v, qtsrc_o_v, thvlsrc_o_v, thlsrc_o_v, usrc_o_v, vsrc_o_v,
        thv0lcl_o_v, trsrc_o_p, cin_p, cinlcl_p, kinv_p, klcl_p, klfc_p,
        plcl_p, plfc_p, tkeavg_p, thvlmin_p, qtsrc_p, thvlsrc_p, thlsrc_p,
        usrc_p, vsrc_p, thv0lcl_p, trsrc_p
    )
    uwshcu_delcin_env_restore_shell_codon(
        mkx, ncnst, wtrc_nwset, qv0_o_p, ql0_o_p, qi0_o_p, t0_o_p, s0_o_p,
        u0_o_p, v0_o_p, qt0_o_p, thl0_o_p, thvl0_o_p, ssthl0_o_p, ssqt0_o_p,
        thv0bot_o_p, thv0top_o_p, thvl0bot_o_p, thvl0top_o_p, ssu0_o_p,
        ssv0_o_p, tr0_o_p, sstr0_o_p, sswt0_o_p, qv0_p, ql0_p, qi0_p, t0_p,
        s0_p, u0_p, v0_p, qt0_p, thl0_p, thvl0_p, ssthl0_p, ssqt0_p,
        thv0bot_p, thv0top_p, thvl0bot_p, thvl0top_p, ssu0_p, ssv0_p, tr0_p,
        sstr0_p, sswt0_p
    )
    uwshcu_delcin_workspace_reset_shell_codon(
        mkx, ncnst, wtrc_nwset, umf_p, emf_p, slflx_p, qtflx_p, uflx_p, vflx_p,
        qvten_p, qlten_p, qiten_p, sten_p, uten_p, vten_p, qrten_p, qsten_p,
        dwten_p, diten_p, evapc_p, cufrc_p, qcu_p, qlu_p, qiu_p, fer_p, fdr_p,
        qc_p, qc_l_p, qc_i_p, qtten_p, slten_p, ufrc_p, thlu_p, qtu_p, uu_p,
        vu_p, wu_p, thvu_p, thlu_emf_p, qtu_emf_p, uu_emf_p, vu_emf_p,
        trflx_p, trten_p, tru_p, tru_emf_p, wtdwten_p, wtditen_p, wtrpten_p,
        wtspten_p, wtqc_liq_p, wtqc_ice_p, wtu_p, wtu_emf_p, wtflx_p,
        wttotten_p, wtprec_p, wtsnow_p, excessu_p, excess0_p, xc_p, aquad_p,
        bquad_p, cquad_p, bogbot_p, bogtop_p
    )

    precip = Ptr[float](precip_p)
    snow = Ptr[float](snow_p)
    rliq = Ptr[float](rliq_p)
    cbmf = Ptr[float](cbmf_p)
    cnt = Ptr[float](cnt_p)
    cnb = Ptr[float](cnb_p)
    ufrcinvbase = Ptr[float](ufrcinvbase_p)
    ufrclcl = Ptr[float](ufrclcl_p)
    winvbase = Ptr[float](winvbase_p)
    wlcl = Ptr[float](wlcl_p)
    emfkbup = Ptr[float](emfkbup_p)
    cbmflimit = Ptr[float](cbmflimit_p)

    precip[0] = 0.0
    snow[0] = 0.0
    rliq[0] = 0.0
    cbmf[0] = 0.0
    cnt[0] = cnt_reset_v
    cnb[0] = 0.0
    ufrcinvbase[0] = 0.0
    ufrclcl[0] = 0.0
    winvbase[0] = 0.0
    wlcl[0] = 0.0
    emfkbup[0] = 0.0
    cbmflimit[0] = 0.0


@export
def uwshcu_inv_prep_shell_codon(
    mix: int,
    mkx: int,
    iend: int,
    ncnst: int,
    ps0_inv_p: cobj,
    zs0_inv_p: cobj,
    p0_inv_p: cobj,
    z0_inv_p: cobj,
    dp0_inv_p: cobj,
    dpdry0_inv_p: cobj,
    u0_inv_p: cobj,
    v0_inv_p: cobj,
    qv0_inv_p: cobj,
    ql0_inv_p: cobj,
    qi0_inv_p: cobj,
    t0_inv_p: cobj,
    s0_inv_p: cobj,
    tr0_inv_p: cobj,
    tke_inv_p: cobj,
    cldfrct_inv_p: cobj,
    concldfrct_inv_p: cobj,
    ps0_p: cobj,
    zs0_p: cobj,
    p0_p: cobj,
    z0_p: cobj,
    dp0_p: cobj,
    dpdry0_p: cobj,
    u0_p: cobj,
    v0_p: cobj,
    qv0_p: cobj,
    ql0_p: cobj,
    qi0_p: cobj,
    t0_p: cobj,
    s0_p: cobj,
    tr0_p: cobj,
    tke_p: cobj,
    cldfrct_p: cobj,
    concldfrct_p: cobj,
):
    ps0_inv = Ptr[float](ps0_inv_p)
    zs0_inv = Ptr[float](zs0_inv_p)
    p0_inv = Ptr[float](p0_inv_p)
    z0_inv = Ptr[float](z0_inv_p)
    dp0_inv = Ptr[float](dp0_inv_p)
    dpdry0_inv = Ptr[float](dpdry0_inv_p)
    u0_inv = Ptr[float](u0_inv_p)
    v0_inv = Ptr[float](v0_inv_p)
    qv0_inv = Ptr[float](qv0_inv_p)
    ql0_inv = Ptr[float](ql0_inv_p)
    qi0_inv = Ptr[float](qi0_inv_p)
    t0_inv = Ptr[float](t0_inv_p)
    s0_inv = Ptr[float](s0_inv_p)
    tr0_inv = Ptr[float](tr0_inv_p)
    tke_inv = Ptr[float](tke_inv_p)
    cldfrct_inv = Ptr[float](cldfrct_inv_p)
    concldfrct_inv = Ptr[float](concldfrct_inv_p)
    ps0 = Ptr[float](ps0_p)
    zs0 = Ptr[float](zs0_p)
    p0 = Ptr[float](p0_p)
    z0 = Ptr[float](z0_p)
    dp0 = Ptr[float](dp0_p)
    dpdry0 = Ptr[float](dpdry0_p)
    u0 = Ptr[float](u0_p)
    v0 = Ptr[float](v0_p)
    qv0 = Ptr[float](qv0_p)
    ql0 = Ptr[float](ql0_p)
    qi0 = Ptr[float](qi0_p)
    t0 = Ptr[float](t0_p)
    s0 = Ptr[float](s0_p)
    tr0 = Ptr[float](tr0_p)
    tke = Ptr[float](tke_p)
    cldfrct = Ptr[float](cldfrct_p)
    concldfrct = Ptr[float](concldfrct_p)

    k = 0
    while k < mkx:
        src_k = mkx - k - 1
        i = 0
        while i < iend:
            dst = i + k * mix
            src = i + src_k * mix
            p0[dst] = p0_inv[src]
            u0[dst] = u0_inv[src]
            v0[dst] = v0_inv[src]
            z0[dst] = z0_inv[src]
            dp0[dst] = dp0_inv[src]
            dpdry0[dst] = dpdry0_inv[src]
            qv0[dst] = qv0_inv[src]
            ql0[dst] = ql0_inv[src]
            qi0[dst] = qi0_inv[src]
            t0[dst] = t0_inv[src]
            s0[dst] = s0_inv[src]
            cldfrct[dst] = cldfrct_inv[src]
            concldfrct[dst] = concldfrct_inv[src]
            i += 1
        m = 0
        while m < ncnst:
            offset = m * mix * mkx
            i = 0
            while i < iend:
                tr0[i + k * mix + offset] = tr0_inv[i + src_k * mix + offset]
                i += 1
            m += 1
        k += 1

    k = 0
    while k <= mkx:
        src_k = mkx - k
        i = 0
        while i < iend:
            dst = i + k * mix
            src = i + src_k * mix
            ps0[dst] = ps0_inv[src]
            zs0[dst] = zs0_inv[src]
            tke[dst] = tke_inv[src]
            i += 1
        k += 1


@export
def uwshcu_inv_post_shell_codon(
    mix: int,
    mkx: int,
    iend: int,
    ncnst: int,
    umf_p: cobj,
    qvten_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    sten_p: cobj,
    uten_p: cobj,
    vten_p: cobj,
    trten_p: cobj,
    qrten_p: cobj,
    qsten_p: cobj,
    evapc_p: cobj,
    slflx_p: cobj,
    qtflx_p: cobj,
    flxprc1_p: cobj,
    flxsnow1_p: cobj,
    cufrc_p: cobj,
    qcu_p: cobj,
    qlu_p: cobj,
    qiu_p: cobj,
    qc_p: cobj,
    cnt_p: cobj,
    cnb_p: cobj,
    wtqc_p: cobj,
    umf_inv_p: cobj,
    qvten_inv_p: cobj,
    qlten_inv_p: cobj,
    qiten_inv_p: cobj,
    sten_inv_p: cobj,
    uten_inv_p: cobj,
    vten_inv_p: cobj,
    trten_inv_p: cobj,
    qrten_inv_p: cobj,
    qsten_inv_p: cobj,
    evapc_inv_p: cobj,
    slflx_inv_p: cobj,
    qtflx_inv_p: cobj,
    flxprc1_inv_p: cobj,
    flxsnow1_inv_p: cobj,
    cufrc_inv_p: cobj,
    qcu_inv_p: cobj,
    qlu_inv_p: cobj,
    qiu_inv_p: cobj,
    qc_inv_p: cobj,
    cnt_inv_p: cobj,
    cnb_inv_p: cobj,
    wtqc_inv_p: cobj,
):
    umf = Ptr[float](umf_p)
    qvten = Ptr[float](qvten_p)
    qlten = Ptr[float](qlten_p)
    qiten = Ptr[float](qiten_p)
    sten = Ptr[float](sten_p)
    uten = Ptr[float](uten_p)
    vten = Ptr[float](vten_p)
    trten = Ptr[float](trten_p)
    qrten = Ptr[float](qrten_p)
    qsten = Ptr[float](qsten_p)
    evapc = Ptr[float](evapc_p)
    slflx = Ptr[float](slflx_p)
    qtflx = Ptr[float](qtflx_p)
    flxprc1 = Ptr[float](flxprc1_p)
    flxsnow1 = Ptr[float](flxsnow1_p)
    cufrc = Ptr[float](cufrc_p)
    qcu = Ptr[float](qcu_p)
    qlu = Ptr[float](qlu_p)
    qiu = Ptr[float](qiu_p)
    qc = Ptr[float](qc_p)
    cnt = Ptr[float](cnt_p)
    cnb = Ptr[float](cnb_p)
    wtqc = Ptr[float](wtqc_p)
    umf_inv = Ptr[float](umf_inv_p)
    qvten_inv = Ptr[float](qvten_inv_p)
    qlten_inv = Ptr[float](qlten_inv_p)
    qiten_inv = Ptr[float](qiten_inv_p)
    sten_inv = Ptr[float](sten_inv_p)
    uten_inv = Ptr[float](uten_inv_p)
    vten_inv = Ptr[float](vten_inv_p)
    trten_inv = Ptr[float](trten_inv_p)
    qrten_inv = Ptr[float](qrten_inv_p)
    qsten_inv = Ptr[float](qsten_inv_p)
    evapc_inv = Ptr[float](evapc_inv_p)
    slflx_inv = Ptr[float](slflx_inv_p)
    qtflx_inv = Ptr[float](qtflx_inv_p)
    flxprc1_inv = Ptr[float](flxprc1_inv_p)
    flxsnow1_inv = Ptr[float](flxsnow1_inv_p)
    cufrc_inv = Ptr[float](cufrc_inv_p)
    qcu_inv = Ptr[float](qcu_inv_p)
    qlu_inv = Ptr[float](qlu_inv_p)
    qiu_inv = Ptr[float](qiu_inv_p)
    qc_inv = Ptr[float](qc_inv_p)
    cnt_inv = Ptr[float](cnt_inv_p)
    cnb_inv = Ptr[float](cnb_inv_p)
    wtqc_inv = Ptr[float](wtqc_inv_p)

    i = 0
    while i < iend:
        cnt_inv[i] = float(mkx + 1) - cnt[i]
        cnb_inv[i] = float(mkx + 1) - cnb[i]
        i += 1

    k = 0
    while k <= mkx:
        dst_k = mkx - k
        i = 0
        while i < iend:
            src = i + k * mix
            dst = i + dst_k * mix
            umf_inv[dst] = umf[src]
            slflx_inv[dst] = slflx[src]
            qtflx_inv[dst] = qtflx[src]
            flxprc1_inv[dst] = flxprc1[src]
            flxsnow1_inv[dst] = flxsnow1[src]
            i += 1
        k += 1

    k = 0
    while k < mkx:
        dst_k = mkx - k - 1
        i = 0
        while i < iend:
            src = i + k * mix
            dst = i + dst_k * mix
            qvten_inv[dst] = qvten[src]
            qlten_inv[dst] = qlten[src]
            qiten_inv[dst] = qiten[src]
            sten_inv[dst] = sten[src]
            uten_inv[dst] = uten[src]
            vten_inv[dst] = vten[src]
            qrten_inv[dst] = qrten[src]
            qsten_inv[dst] = qsten[src]
            evapc_inv[dst] = evapc[src]
            cufrc_inv[dst] = cufrc[src]
            qcu_inv[dst] = qcu[src]
            qlu_inv[dst] = qlu[src]
            qiu_inv[dst] = qiu[src]
            qc_inv[dst] = qc[src]
            i += 1
        m = 0
        while m < ncnst:
            offset = m * mix * mkx
            i = 0
            while i < iend:
                trten_inv[i + dst_k * mix + offset] = trten[i + k * mix + offset]
                wtqc_inv[i + dst_k * mix + offset] = wtqc[i + k * mix + offset]
                i += 1
            m += 1
        k += 1


@export
def convect_shallow_diag_shell_codon(
    mode: int,
    ncol: int,
    pcols: int,
    pver: int,
    ixcldliq: int,
    ixcldice: int,
    ztodt: float,
    latvap: float,
    latice: float,
    zvir: float,
    state_s_p: cobj,
    state_t_p: cobj,
    state_q_p: cobj,
    sat_rh_p: cobj,
    sl_p: cobj,
    qt_p: cobj,
    slv_p: cobj,
    t_precu_p: cobj,
    rh_precu_p: cobj,
    tten_p: cobj,
    rhten_p: cobj,
):
    state_s = Ptr[float](state_s_p)
    state_t = Ptr[float](state_t_p)
    state_q = Ptr[float](state_q_p)
    sat_rh = Ptr[float](sat_rh_p)
    sl = Ptr[float](sl_p)
    qt = Ptr[float](qt_p)
    slv = Ptr[float](slv_p)
    t_precu = Ptr[float](t_precu_p)
    rh_precu = Ptr[float](rh_precu_p)
    tten = Ptr[float](tten_p)
    rhten = Ptr[float](rhten_p)

    if mode == 1 or mode == 3:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, pcols)
                sl[idx] = (
                    state_s[idx]
                    - latvap * state_q[_idx3(i, k, ixcldliq, pcols, pver)]
                    - (latvap + latice) * state_q[_idx3(i, k, ixcldice, pcols, pver)]
                )
                qt[idx] = (
                    state_q[_idx3(i, k, 1, pcols, pver)]
                    + state_q[_idx3(i, k, ixcldliq, pcols, pver)]
                    + state_q[_idx3(i, k, ixcldice, pcols, pver)]
                )
                slv[idx] = sl[idx] * (1.0 + zvir * qt[idx])
                if mode == 1:
                    t_precu[idx] = state_t[idx]
    elif mode == 2:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, pcols)
                rh_precu[idx] = state_q[_idx3(i, k, 1, pcols, pver)] / sat_rh[idx] * 100.0
    elif mode == 4:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, pcols)
                sat_rh[idx] = state_q[_idx3(i, k, 1, pcols, pver)] / sat_rh[idx] * 100.0
                tten[idx] = (state_t[idx] - t_precu[idx]) / ztodt
                rhten[idx] = (sat_rh[idx] - rh_precu[idx]) / ztodt


@export
def convect_shallow_wtrc_precip_shell_codon(
    pcols: int,
    vap_idx: int,
    wtprect_p: cobj,
    wtsnowt_p: cobj,
    wtprec_p: cobj,
    wtsnow_p: cobj,
):
    wtprect = Ptr[float](wtprect_p)
    wtsnowt = Ptr[float](wtsnowt_p)
    wtprec = Ptr[float](wtprec_p)
    wtsnow = Ptr[float](wtsnow_p)

    i = 1
    while i <= pcols:
        src_idx = _idx2(i, vap_idx, pcols)
        dst_idx = i - 1
        wtprec[dst_idx] = wtprec[dst_idx] + (wtprect[src_idx] - wtsnowt[src_idx])
        wtsnow[dst_idx] = wtsnow[dst_idx] + wtsnowt[src_idx]
        i += 1


@export
def convect_shallow_uw_post_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    pcnst: int,
    wtrc_nwset: int,
    latvap: float,
    cpair: float,
    state_pmid_p: cobj,
    cmfmc_p: cobj,
    cmfmc2_p: cobj,
    cnt_p: cobj,
    cnt2_p: cobj,
    cnb_p: cobj,
    cnb2_p: cobj,
    pcnt_p: cobj,
    pcnb_p: cobj,
    qc_p: cobj,
    qc2_p: cobj,
    rliq_p: cobj,
    rliq2_p: cobj,
    wtqc_p: cobj,
    wtdlf_p: cobj,
    freqsh_p: cobj,
    icwmr_p: cobj,
    iccmr_uw_p: cobj,
    rprdsh_p: cobj,
    cmfdqs_p: cobj,
    ptend_q_p: cobj,
    ptend_tracer_p: cobj,
    cmfsl_p: cobj,
    cmflq_p: cobj,
    slflx_p: cobj,
    qtflx_p: cobj,
    rprddp_p: cobj,
    rprdtot_p: cobj,
    ptend_s_p: cobj,
    ftem_p: cobj,
    liq_type_p: cobj,
    ice_type_p: cobj,
):
    state_pmid = Ptr[float](state_pmid_p)
    cmfmc = Ptr[float](cmfmc_p)
    cmfmc2 = Ptr[float](cmfmc2_p)
    cnt = Ptr[float](cnt_p)
    cnt2 = Ptr[float](cnt2_p)
    cnb = Ptr[float](cnb_p)
    cnb2 = Ptr[float](cnb2_p)
    pcnt = Ptr[float](pcnt_p)
    pcnb = Ptr[float](pcnb_p)
    qc = Ptr[float](qc_p)
    qc2 = Ptr[float](qc2_p)
    rliq = Ptr[float](rliq_p)
    rliq2 = Ptr[float](rliq2_p)
    wtqc = Ptr[float](wtqc_p)
    wtdlf = Ptr[float](wtdlf_p)
    freqsh = Ptr[float](freqsh_p)
    icwmr = Ptr[float](icwmr_p)
    iccmr_uw = Ptr[float](iccmr_uw_p)
    rprdsh = Ptr[float](rprdsh_p)
    cmfdqs = Ptr[float](cmfdqs_p)
    ptend_q = Ptr[float](ptend_q_p)
    ptend_tracer = Ptr[float](ptend_tracer_p)
    cmfsl = Ptr[float](cmfsl_p)
    cmflq = Ptr[float](cmflq_p)
    slflx = Ptr[float](slflx_p)
    qtflx = Ptr[float](qtflx_p)
    rprddp = Ptr[float](rprddp_p)
    rprdtot = Ptr[float](rprdtot_p)
    ptend_s = Ptr[float](ptend_s_p)
    ftem = Ptr[float](ftem_p)
    liq_type = Ptr[int](liq_type_p)
    ice_type = Ptr[int](ice_type_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            icwmr[idx] = iccmr_uw[idx]
            rprdsh[idx] = rprdsh[idx] + cmfdqs[idx]

    for m in range(4, pcnst + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx3(i, k, m, pcols, pver)
                ptend_q[idx] = ptend_tracer[idx]

    for k in range(1, pverp + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            cmfsl[idx] = slflx[idx]
            cmflq[idx] = qtflx[idx] * latvap

    for i in range(1, pcols + 1):
        freqsh[i - 1] = 0.0
    for i in range(1, ncol + 1):
        active = 0
        for k in range(1, pver + 1):
            if cmfmc2[_idx2(i, k, pcols)] > 0.0:
                active = 1
                break
        if active == 0:
            freqsh[i - 1] = 1.0

    for k in range(1, pverp + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            cmfmc[idx] = cmfmc[idx] + cmfmc2[idx]

    for i in range(1, ncol + 1):
        if cnt2[i - 1] < cnt[i - 1]:
            cnt[i - 1] = cnt2[i - 1]
        if cnb2[i - 1] > cnb[i - 1]:
            cnb[i - 1] = cnb2[i - 1]
        pcnt[i - 1] = state_pmid[_idx2(i, int(cnt[i - 1]), pcols)]
        pcnb[i - 1] = state_pmid[_idx2(i, int(cnb[i - 1]), pcols)]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            qc[idx] = qc[idx] + qc2[idx]

    for i in range(1, ncol + 1):
        rliq[i - 1] = rliq[i - 1] + rliq2[i - 1]

    for m in range(1, wtrc_nwset + 1):
        liq = liq_type[m - 1]
        ice = ice_type[m - 1]
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx_out = _idx3(i, k, m, pcols, pver)
                tmp_sum = (
                    wtqc[_idx3(i, k, liq, pcols, pver)]
                    + wtqc[_idx3(i, k, ice, pcols, pver)]
                )
                wtdlf[idx_out] = wtdlf[idx_out] + tmp_sum

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            rprdtot[idx] = rprdsh[idx] + rprddp[idx]
            ftem[idx] = ptend_s[idx] / cpair


@export
def convect_shallow_postmerge_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    pcnst: int,
    wtrc_nwset: int,
    state_pmid_p: cobj,
    cmfmc_p: cobj,
    cmfmc2_p: cobj,
    cnt_p: cobj,
    cnt2_p: cobj,
    cnb_p: cobj,
    cnb2_p: cobj,
    pcnt_p: cobj,
    pcnb_p: cobj,
    qc_p: cobj,
    qc2_p: cobj,
    rliq_p: cobj,
    rliq2_p: cobj,
    wtqc_p: cobj,
    wtdlf_p: cobj,
    freqsh_p: cobj,
    liq_type_p: cobj,
    ice_type_p: cobj,
):
    state_pmid = Ptr[float](state_pmid_p)
    cmfmc = Ptr[float](cmfmc_p)
    cmfmc2 = Ptr[float](cmfmc2_p)
    cnt = Ptr[float](cnt_p)
    cnt2 = Ptr[float](cnt2_p)
    cnb = Ptr[float](cnb_p)
    cnb2 = Ptr[float](cnb2_p)
    pcnt = Ptr[float](pcnt_p)
    pcnb = Ptr[float](pcnb_p)
    qc = Ptr[float](qc_p)
    qc2 = Ptr[float](qc2_p)
    rliq = Ptr[float](rliq_p)
    rliq2 = Ptr[float](rliq2_p)
    wtqc = Ptr[float](wtqc_p)
    wtdlf = Ptr[float](wtdlf_p)
    freqsh = Ptr[float](freqsh_p)
    liq_type = Ptr[int](liq_type_p)
    ice_type = Ptr[int](ice_type_p)

    for i in range(1, pcols + 1):
        freqsh[i - 1] = 0.0
    for i in range(1, ncol + 1):
        active = 0
        for k in range(1, pver + 1):
            if cmfmc2[_idx2(i, k, pcols)] > 0.0:
                active = 1
                break
        if active == 0:
            freqsh[i - 1] = 1.0

    for k in range(1, pverp + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            cmfmc[idx] = cmfmc[idx] + cmfmc2[idx]

    for i in range(1, ncol + 1):
        if cnt2[i - 1] < cnt[i - 1]:
            cnt[i - 1] = cnt2[i - 1]
        if cnb2[i - 1] > cnb[i - 1]:
            cnb[i - 1] = cnb2[i - 1]
        pcnt[i - 1] = state_pmid[_idx2(i, int(cnt[i - 1]), pcols)]
        pcnb[i - 1] = state_pmid[_idx2(i, int(cnb[i - 1]), pcols)]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            qc[idx] = qc[idx] + qc2[idx]

    for i in range(1, ncol + 1):
        rliq[i - 1] = rliq[i - 1] + rliq2[i - 1]

    for m in range(1, wtrc_nwset + 1):
        liq = liq_type[m - 1]
        ice = ice_type[m - 1]
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx_out = _idx3(i, k, m, pcols, pver)
                tmp_sum = (
                    wtqc[_idx3(i, k, liq, pcols, pver)]
                    + wtqc[_idx3(i, k, ice, pcols, pver)]
                )
                wtdlf[idx_out] = wtdlf[idx_out] + tmp_sum
