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
