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
