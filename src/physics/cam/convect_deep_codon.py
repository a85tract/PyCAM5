from math import log, sqrt


@export
def convect_deep_select_scheme_codon(
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

        if c1 == 122 and c2 == 109:
            scheme_code[0] = 1
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
            scheme_code[0] = 3
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
            scheme_code[0] = 4
            return

    status[0] = 1


@export
def zm_convr_post_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    lengath: int,
    gravit: float,
    cpair: float,
    mcon_p: cobj,
    mu_p: cobj,
    md_p: cobj,
    ideep_p: cobj,
    jt_p: cobj,
    maxg_p: cobj,
    ptend_s_p: cobj,
    state_ps_p: cobj,
    state_pmid_p: cobj,
    freqzm_p: cobj,
    mu_out_p: cobj,
    md_out_p: cobj,
    ftem_p: cobj,
    pcont_p: cobj,
    pconb_p: cobj,
):
    mcon = Ptr[float](mcon_p)
    mu = Ptr[float](mu_p)
    md = Ptr[float](md_p)
    ideep = Ptr[int](ideep_p)
    jt = Ptr[int](jt_p)
    maxg = Ptr[int](maxg_p)
    ptend_s = Ptr[float](ptend_s_p)
    state_ps = Ptr[float](state_ps_p)
    state_pmid = Ptr[float](state_pmid_p)
    freqzm = Ptr[float](freqzm_p)
    mu_out = Ptr[float](mu_out_p)
    md_out = Ptr[float](md_out_p)
    ftem = Ptr[float](ftem_p)
    pcont = Ptr[float](pcont_p)
    pconb = Ptr[float](pconb_p)

    i = 0
    while i < pcols:
        freqzm[i] = 0.0
        i += 1

    i = 0
    while i < lengath:
        freqzm[ideep[i] - 1] = 1.0
        i += 1

    k = 0
    while k < pver:
        i = 0
        while i < ncol:
            idx = i + k * pcols
            mcon[idx] = mcon[idx] * 100.0 / gravit
            ftem[idx] = ptend_s[idx] / cpair
            i += 1
        k += 1

    k = 0
    while k < pver:
        i = 0
        while i < pcols:
            idx = i + k * pcols
            mu_out[idx] = 0.0
            md_out[idx] = 0.0
            i += 1
        k += 1

    i = 0
    while i < lengath:
        ii = ideep[i] - 1
        k = 0
        while k < pver:
            gathered_idx = i + k * pcols
            out_idx = ii + k * pcols
            mu_out[out_idx] = mu[gathered_idx] * 100.0 / gravit
            md_out[out_idx] = md[gathered_idx] * 100.0 / gravit
            k += 1
        i += 1

    i = 0
    while i < ncol:
        pcont[i] = state_ps[i]
        pconb[i] = state_ps[i]
        i += 1

    i = 0
    while i < lengath:
        if maxg[i] > jt[i]:
            ii = ideep[i] - 1
            pcont[ii] = state_pmid[ii + (jt[i] - 1) * pcols]
            pconb[ii] = state_pmid[ii + (maxg[i] - 1) * pcols]
        i += 1


@export
def zm_conv_workspace_init_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    ftem_p: cobj,
    mu_out_p: cobj,
    md_out_p: cobj,
    wind_tends_p: cobj,
):
    ftem = Ptr[float](ftem_p)
    mu_out = Ptr[float](mu_out_p)
    md_out = Ptr[float](md_out_p)
    wind_tends = Ptr[float](wind_tends_p)
    plane = pcols * pver

    k = 0
    while k < pver:
        i = 0
        while i < pcols:
            idx = i + k * pcols
            ftem[idx] = 0.0
            mu_out[idx] = 0.0
            md_out[idx] = 0.0
            i += 1
        k += 1

    k = 0
    while k < pver:
        i = 0
        while i < ncol:
            idx = i + k * pcols
            wind_tends[idx] = 0.0
            wind_tends[idx + plane] = 0.0
            i += 1
        k += 1


@export
def zm_wtrc_convr_prep_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    wtrc_nwset: int,
    wtdlf_p: cobj,
    wtrprd_p: cobj,
    wtprect_p: cobj,
    rprd_p: cobj,
    prec_p: cobj,
):
    wtdlf = Ptr[float](wtdlf_p)
    wtrprd = Ptr[float](wtrprd_p)
    wtprect = Ptr[float](wtprect_p)
    rprd = Ptr[float](rprd_p)
    prec = Ptr[float](prec_p)
    plane = pcols * pver

    m = 0
    while m < wtrc_nwset:
        k = 0
        while k < pver:
            i = 0
            while i < ncol:
                wtdlf[i + k * pcols + m * plane] = 0.0
                i += 1
            k += 1
        m += 1

    m = 0
    while m < pcnst:
        k = 0
        while k < pver:
            i = 0
            while i < ncol:
                wtrprd[i + k * pcols + m * plane] = 0.0
                i += 1
            k += 1
        i = 0
        while i < ncol:
            wtprect[i + m * pcols] = 0.0
            i += 1
        m += 1

    k = 0
    while k < pver:
        i = 0
        while i < ncol:
            idx = i + k * pcols
            wtrprd[idx] = rprd[idx]
            i += 1
        k += 1

    i = 0
    while i < ncol:
        wtprect[i] = prec[i]
        i += 1


@export
def zm_wtrc_precip_assign_shell_codon(
    pcols: int,
    vap_type: int,
    wtprect_p: cobj,
    wtsnowt_p: cobj,
    wtprec_p: cobj,
    wtsnow_p: cobj,
):
    wtprect = Ptr[float](wtprect_p)
    wtsnowt = Ptr[float](wtsnowt_p)
    wtprec = Ptr[float](wtprec_p)
    wtsnow = Ptr[float](wtsnow_p)
    offset = (vap_type - 1) * pcols

    i = 0
    while i < pcols:
        snow_val = wtsnowt[i + offset]
        wtprec[i] = wtprect[i + offset] - snow_val
        wtsnow[i] = snow_val
        i += 1


@export
def zm_conv_evap_prep_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    dp_cldliq_p: cobj,
    dp_cldice_p: cobj,
):
    dp_cldliq = Ptr[float](dp_cldliq_p)
    dp_cldice = Ptr[float](dp_cldice_p)

    k = 0
    while k < pver:
        i = 0
        while i < ncol:
            idx = i + k * pcols
            dp_cldliq[idx] = 0.0
            dp_cldice[idx] = 0.0
            i += 1
        k += 1


@export
def zm_conv_evap_post_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    ixorg: int,
    do_org: int,
    ztodt: float,
    evapcdp_p: cobj,
    ptend_q_p: cobj,
    state_q_p: cobj,
):
    evapcdp = Ptr[float](evapcdp_p)
    ptend_q = Ptr[float](ptend_q_p)
    state_q = Ptr[float](state_q_p)
    plane = pcols * pver

    k = 0
    while k < pver:
        i = 0
        while i < ncol:
            idx = i + k * pcols
            evapcdp[idx] = ptend_q[idx]
            i += 1
        k += 1

    if do_org != 0:
        org_offset = (ixorg - 1) * plane
        k = 0
        while k < pver:
            i = 0
            while i < ncol:
                idx = i + k * pcols
                org_idx = idx + org_offset
                val = (50.0 * 1000.0 * 1000.0 * abs(evapcdp[idx])) - (state_q[org_idx] / 10800.0)
                if val < 0.0:
                    val = 0.0
                if val > 1.0:
                    val = 1.0
                ptend_q[org_idx] = (val - state_q[org_idx]) / ztodt
                i += 1
            k += 1


@export
def zm_conv_evap_hist_shell_codon(
    mode: int,
    ncol: int,
    pcols: int,
    pver: int,
    cpair: float,
    ptend_s_p: cobj,
    tend_s_snwprd_p: cobj,
    tend_s_snwevmlt_p: cobj,
    ftem_p: cobj,
):
    ptend_s = Ptr[float](ptend_s_p)
    tend_s_snwprd = Ptr[float](tend_s_snwprd_p)
    tend_s_snwevmlt = Ptr[float](tend_s_snwevmlt_p)
    ftem = Ptr[float](ftem_p)

    k = 0
    while k < pver:
        i = 0
        while i < ncol:
            idx = i + k * pcols
            if mode == 1:
                ftem[idx] = ptend_s[idx] / cpair
            elif mode == 2:
                ftem[idx] = tend_s_snwprd[idx] / cpair
            else:
                ftem[idx] = tend_s_snwevmlt[idx] / cpair
            i += 1
        k += 1


@export
def zm_conv_evap_main_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    pergro: int,
    do_org: int,
    gravit: float,
    latvap: float,
    latice: float,
    tmelt: float,
    ke: float,
    ke_lnd: float,
    deltat: float,
    t_p: cobj,
    q_p: cobj,
    pdel_p: cobj,
    landfrac_p: cobj,
    prdprec_p: cobj,
    cldfrc_p: cobj,
    qs_p: cobj,
    fsnow_conv_p: cobj,
    prec_p: cobj,
    snow_p: cobj,
    tend_s_p: cobj,
    tend_q_p: cobj,
    tend_s_snwprd_p: cobj,
    tend_s_snwevmlt_p: cobj,
    evpstore_p: cobj,
    substore_p: cobj,
    ntprprd_p: cobj,
    ntsnprd_p: cobj,
    flxprec_p: cobj,
    flxsnow_p: cobj,
    evpvint_p: cobj,
):
    t = Ptr[float](t_p)
    q = Ptr[float](q_p)
    pdel = Ptr[float](pdel_p)
    landfrac = Ptr[float](landfrac_p)
    prdprec = Ptr[float](prdprec_p)
    cldfrc = Ptr[float](cldfrc_p)
    qs = Ptr[float](qs_p)
    fsnow_conv = Ptr[float](fsnow_conv_p)
    prec = Ptr[float](prec_p)
    snow = Ptr[float](snow_p)
    tend_s = Ptr[float](tend_s_p)
    tend_q = Ptr[float](tend_q_p)
    tend_s_snwprd = Ptr[float](tend_s_snwprd_p)
    tend_s_snwevmlt = Ptr[float](tend_s_snwevmlt_p)
    evpstore = Ptr[float](evpstore_p)
    substore = Ptr[float](substore_p)
    ntprprd = Ptr[float](ntprprd_p)
    ntsnprd = Ptr[float](ntsnprd_p)
    flxprec = Ptr[float](flxprec_p)
    flxsnow = Ptr[float](flxsnow_p)
    evpvint = Ptr[float](evpvint_p)

    i = 0
    while i < ncol:
        flxprec[i] = 0.0
        flxsnow[i] = 0.0
        evpvint[i] = 0.0
        i += 1

    k = 0
    while k < pver:
        i = 0
        while i < pcols:
            idx = i + k * pcols
            evpstore[idx] = 0.0
            substore[idx] = 0.0
            i += 1
        k += 1

    k = 0
    while k < pver:
        i = 0
        while i < ncol:
            idx = i + k * pcols
            flx_idx = i + k * pcols
            flx_next_idx = i + (k + 1) * pcols

            if t[idx] > tmelt:
                flxsntm = 0.0
                snowmlt = flxsnow[flx_idx] * gravit / pdel[idx]
            else:
                flxsntm = flxsnow[flx_idx]
                snowmlt = 0.0

            evplimit = 1.0 - q[idx] / qs[idx]
            if evplimit < 0.0:
                evplimit = 0.0

            if do_org != 0:
                kemask = ke * (1.0 - landfrac[i]) + ke_lnd * landfrac[i]
            else:
                kemask = ke

            evpprec = kemask * (1.0 - cldfrc[idx]) * evplimit * sqrt(flxprec[flx_idx])

            evplimit = (qs[idx] - q[idx]) / deltat
            if evplimit < 0.0:
                evplimit = 0.0

            limit2 = flxprec[flx_idx] * gravit / pdel[idx]
            if limit2 < evplimit:
                evplimit = limit2

            limit2 = (prec[i] - evpvint[i]) * gravit / pdel[idx]
            if limit2 < evplimit:
                evplimit = limit2

            if evplimit < evpprec:
                evpprec = evplimit

            if flxprec[flx_idx] > 0.0:
                work1 = flxsntm / flxprec[flx_idx]
                if work1 < 0.0:
                    work1 = 0.0
                if work1 > 1.0:
                    work1 = 1.0
                evpsnow = evpprec * work1
            else:
                evpsnow = 0.0

            evpstore[idx] = evpprec
            substore[idx] = evpsnow

            evpvint[i] = evpvint[i] + evpprec * pdel[idx] / gravit
            ntprprd[idx] = prdprec[idx] - evpprec

            if pergro != 0:
                work1 = flxsnow[flx_idx] / (flxprec[flx_idx] + 8.64e-11)
                if work1 < 0.0:
                    work1 = 0.0
                if work1 > 1.0:
                    work1 = 1.0
            else:
                if flxprec[flx_idx] > 0.0:
                    work1 = flxsnow[flx_idx] / flxprec[flx_idx]
                    if work1 < 0.0:
                        work1 = 0.0
                    if work1 > 1.0:
                        work1 = 1.0
                else:
                    work1 = 0.0

            if fsnow_conv[idx] > work1:
                work2 = fsnow_conv[idx]
            else:
                work2 = work1
            if snowmlt > 0.0:
                work2 = 0.0

            ntsnprd[idx] = prdprec[idx] * work2 - evpsnow - snowmlt
            tend_s_snwprd[idx] = prdprec[idx] * work2 * latice
            tend_s_snwevmlt[idx] = -(evpsnow + snowmlt) * latice

            flxprec[flx_next_idx] = flxprec[flx_idx] + ntprprd[idx] * pdel[idx] / gravit
            flxsnow[flx_next_idx] = flxsnow[flx_idx] + ntsnprd[idx] * pdel[idx] / gravit

            if flxprec[flx_next_idx] < 0.0:
                flxprec[flx_next_idx] = 0.0
            if flxsnow[flx_next_idx] < 0.0:
                flxsnow[flx_next_idx] = 0.0

            tend_s[idx] = -evpprec * latvap + ntsnprd[idx] * latice
            tend_q[idx] = evpprec
            i += 1
        k += 1

    i = 0
    while i < ncol:
        bottom_idx = i + pver * pcols
        prec[i] = flxprec[bottom_idx] / 1000.0
        snow[i] = flxsnow[bottom_idx] / 1000.0
        i += 1


@export
def zm_momtran_prep_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    state_u_p: cobj,
    state_v_p: cobj,
    winds_p: cobj,
):
    state_u = Ptr[float](state_u_p)
    state_v = Ptr[float](state_v_p)
    winds = Ptr[float](winds_p)
    plane = pcols * pver

    k = 0
    while k < pver:
        i = 0
        while i < ncol:
            idx = i + k * pcols
            winds[idx] = state_u[idx]
            winds[idx + plane] = state_v[idx]
            i += 1
        k += 1


@export
def zm_momtran_post_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    cpair: float,
    wind_tends_p: cobj,
    seten_p: cobj,
    ptend_u_p: cobj,
    ptend_v_p: cobj,
    ptend_s_p: cobj,
    ftem_p: cobj,
):
    wind_tends = Ptr[float](wind_tends_p)
    seten = Ptr[float](seten_p)
    ptend_u = Ptr[float](ptend_u_p)
    ptend_v = Ptr[float](ptend_v_p)
    ptend_s = Ptr[float](ptend_s_p)
    ftem = Ptr[float](ftem_p)
    plane = pcols * pver

    k = 0
    while k < pver:
        i = 0
        while i < ncol:
            idx = i + k * pcols
            ptend_u[idx] = wind_tends[idx]
            ptend_v[idx] = wind_tends[idx + plane]
            ptend_s[idx] = seten[idx]
            ftem[idx] = seten[idx] / cpair
            i += 1
        k += 1


@export
def zm_momtran_main_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    ncnst: int,
    il1g: int,
    il2g: int,
    dt: float,
    domomtran_p: cobj,
    q_p: cobj,
    mu_p: cobj,
    md_p: cobj,
    du_p: cobj,
    eu_p: cobj,
    ed_p: cobj,
    dp_p: cobj,
    dsubcld_p: cobj,
    jt_p: cobj,
    mx_p: cobj,
    ideep_p: cobj,
    dqdt_p: cobj,
    pguall_p: cobj,
    pgdall_p: cobj,
    icwu_p: cobj,
    icwd_p: cobj,
    seten_p: cobj,
    chat_p: cobj,
    cond_p: cobj,
    const_p: cobj,
    conu_p: cobj,
    dcondt_p: cobj,
    mududp_p: cobj,
    mddudp_p: cobj,
    pgu_p: cobj,
    pgd_p: cobj,
    gseten_p: cobj,
    mflux_p: cobj,
    wind0_p: cobj,
    windf_p: cobj,
):
    domomtran = Ptr[int](domomtran_p)
    q = Ptr[float](q_p)
    mu = Ptr[float](mu_p)
    md = Ptr[float](md_p)
    du = Ptr[float](du_p)
    eu = Ptr[float](eu_p)
    ed = Ptr[float](ed_p)
    dp = Ptr[float](dp_p)
    dsubcld = Ptr[float](dsubcld_p)
    jt = Ptr[int](jt_p)
    mx = Ptr[int](mx_p)
    ideep = Ptr[int](ideep_p)
    dqdt = Ptr[float](dqdt_p)
    pguall = Ptr[float](pguall_p)
    pgdall = Ptr[float](pgdall_p)
    icwu = Ptr[float](icwu_p)
    icwd = Ptr[float](icwd_p)
    seten = Ptr[float](seten_p)
    chat = Ptr[float](chat_p)
    cond = Ptr[float](cond_p)
    const = Ptr[float](const_p)
    conu = Ptr[float](conu_p)
    dcondt = Ptr[float](dcondt_p)
    mududp = Ptr[float](mududp_p)
    mddudp = Ptr[float](mddudp_p)
    pgu = Ptr[float](pgu_p)
    pgd = Ptr[float](pgd_p)
    gseten = Ptr[float](gseten_p)
    mflux = Ptr[float](mflux_p)
    wind0 = Ptr[float](wind0_p)
    windf = Ptr[float](windf_p)
    plane = pcols * pver
    flux_plane = pcols * pverp
    ilo = il1g - 1
    ihi = il2g - 1

    m = 0
    while m < ncnst:
        moff = m * plane
        foff = m * flux_plane
        k = 0
        while k < pver:
            i = 0
            while i < pcols:
                idx = i + k * pcols
                pguall[idx + moff] = 0.0
                pgdall[idx + moff] = 0.0
                wind0[idx + moff] = 0.0
                windf[idx + moff] = 0.0
                if i < ncol:
                    icwu[idx + moff] = q[idx + moff]
                    icwd[idx + moff] = q[idx + moff]
                i += 1
            k += 1
        k = 0
        while k < pverp:
            i = 0
            while i < pcols:
                mflux[i + k * pcols + foff] = 0.0
                i += 1
            k += 1
        m += 1

    k = 0
    while k < pver:
        i = 0
        while i < pcols:
            idx = i + k * pcols
            seten[idx] = 0.0
            gseten[idx] = 0.0
            i += 1
        k += 1

    momcu = 0.4
    momcd = 0.4
    mbsth = 1.0e-15

    ktm = pver
    kbm = pver
    i = ilo
    while i <= ihi:
        if jt[i] < ktm:
            ktm = jt[i]
        if mx[i] < kbm:
            kbm = mx[i]
        i += 1

    m = 0
    while m < ncnst:
        moff = m * plane
        foff = m * flux_plane
        if domomtran[m] != 0:
            k = 0
            while k < pver:
                i = ilo
                while i <= ihi:
                    ii = ideep[i] - 1
                    idx = i + k * pcols
                    const[idx] = q[ii + k * pcols + moff]
                    wind0[idx + moff] = const[idx]
                    i += 1
                k += 1

            k = 0
            while k < pver:
                km1 = k - 1
                if km1 < 0:
                    km1 = 0
                i = ilo
                while i <= ihi:
                    idx = i + k * pcols
                    chat[idx] = 0.5 * (const[idx] + const[i + km1 * pcols])
                    conu[idx] = chat[idx]
                    cond[idx] = chat[idx]
                    dcondt[idx] = 0.0
                    i += 1
                k += 1

            i = 0
            while i < il2g:
                pgu[i] = 0.0
                pgd[i] = 0.0
                i += 1

            k = 1
            while k < pver - 1:
                km1 = k - 1
                kp1 = k + 1
                i = ilo
                while i <= ihi:
                    idx = i + k * pcols
                    mududp[idx] = (
                        mu[idx] * (const[idx] - const[i + km1 * pcols]) / dp[i + km1 * pcols]
                        + mu[i + kp1 * pcols] * (const[i + kp1 * pcols] - const[idx]) / dp[idx]
                    )
                    pgu[idx] = -momcu * 0.5 * mududp[idx]
                    mddudp[idx] = (
                        md[idx] * (const[idx] - const[i + km1 * pcols]) / dp[i + km1 * pcols]
                        + md[i + kp1 * pcols] * (const[i + kp1 * pcols] - const[idx]) / dp[idx]
                    )
                    pgd[idx] = -momcd * 0.5 * mddudp[idx]
                    i += 1
                k += 1

            k = pver - 1
            km1 = k - 1
            i = ilo
            while i <= ihi:
                idx = i + k * pcols
                mududp[idx] = mu[idx] * (const[idx] - const[i + km1 * pcols]) / dp[i + km1 * pcols]
                pgu[idx] = -momcu * mududp[idx]
                mddudp[idx] = md[idx] * (const[idx] - const[i + km1 * pcols]) / dp[i + km1 * pcols]
                pgd[idx] = -momcd * mddudp[idx]
                i += 1

            k = 1
            km1 = 0
            kk = pver - 1
            i = ilo
            while i <= ihi:
                kkidx = i + kk * pcols
                mupdudp = mu[kkidx] + du[kkidx] * dp[kkidx]
                if mupdudp > mbsth:
                    conu[kkidx] = (eu[kkidx] * const[kkidx] * dp[kkidx] + pgu[kkidx] * dp[kkidx]) / mupdudp
                idx = i + k * pcols
                if md[idx] < -mbsth:
                    cond[idx] = (-ed[i + km1 * pcols] * const[i + km1 * pcols] * dp[i + km1 * pcols]) - pgd[
                        i + km1 * pcols
                    ] * dp[i + km1 * pcols] / md[idx]
                i += 1

            kk = pver - 2
            while kk >= 0:
                kkp1 = kk + 1
                i = ilo
                while i <= ihi:
                    idx = i + kk * pcols
                    mupdudp = mu[idx] + du[idx] * dp[idx]
                    if mupdudp > mbsth:
                        conu[idx] = (
                            mu[i + kkp1 * pcols] * conu[i + kkp1 * pcols]
                            + eu[idx] * const[idx] * dp[idx]
                            + pgu[idx] * dp[idx]
                        ) / mupdudp
                    i += 1
                kk -= 1

            k = 2
            while k < pver:
                km1 = k - 1
                i = ilo
                while i <= ihi:
                    idx = i + k * pcols
                    if md[idx] < -mbsth:
                        cond[idx] = (
                            md[i + km1 * pcols] * cond[i + km1 * pcols]
                            - ed[i + km1 * pcols] * const[i + km1 * pcols] * dp[i + km1 * pcols]
                            - pgd[i + km1 * pcols] * dp[i + km1 * pcols]
                        ) / md[idx]
                    i += 1
                k += 1

            k = ktm - 1
            while k < pver:
                kp1 = k + 1
                if kp1 >= pver:
                    kp1 = pver - 1
                i = ilo
                while i <= ihi:
                    idx = i + k * pcols
                    dcondt[idx] = (
                        +(
                            mu[i + kp1 * pcols] * (conu[i + kp1 * pcols] - chat[i + kp1 * pcols])
                            - mu[idx] * (conu[idx] - chat[idx])
                            + md[i + kp1 * pcols] * (cond[i + kp1 * pcols] - chat[i + kp1 * pcols])
                            - md[idx] * (cond[idx] - chat[idx])
                        )
                        / dp[idx]
                    )
                    i += 1
                k += 1

            k = kbm - 1
            while k < pver:
                i = ilo
                while i <= ihi:
                    if k + 1 == mx[i]:
                        idx = i + k * pcols
                        dcondt[idx] = (1.0 / dp[idx]) * (
                            -mu[idx] * (conu[idx] - chat[idx]) - md[idx] * (cond[idx] - chat[idx])
                        )
                    i += 1
                k += 1

            k = 0
            while k < pver:
                i = 0
                while i < pcols:
                    dqdt[i + k * pcols + moff] = 0.0
                    i += 1
                k += 1

            k = 0
            while k < pver:
                i = ilo
                while i <= ihi:
                    ii = ideep[i] - 1
                    src_idx = i + k * pcols
                    dst_idx = ii + k * pcols + moff
                    dqdt[dst_idx] = dcondt[src_idx]
                    pguall[dst_idx] = -pgu[src_idx]
                    pgdall[dst_idx] = -pgd[src_idx]
                    icwu[dst_idx] = conu[src_idx]
                    icwd[dst_idx] = cond[src_idx]
                    i += 1
                k += 1

            k = ktm - 1
            while k < pver:
                i = ilo
                while i <= ihi:
                    idx = i + k * pcols
                    mflux[i + k * pcols + foff] = -mu[idx] * (conu[idx] - chat[idx]) - md[idx] * (
                        cond[idx] - chat[idx]
                    )
                    i += 1
                k += 1

            k = ktm - 1
            while k < pver:
                i = ilo
                while i <= ihi:
                    idx = i + k * pcols
                    windf[idx + moff] = const[idx] - (
                        mflux[i + (k + 1) * pcols + foff] - mflux[i + k * pcols + foff]
                    ) * dt / dp[idx]
                    i += 1
                k += 1
        m += 1

    k = ktm - 1
    while k < pver:
        km1 = k - 1
        if km1 < 0:
            km1 = 0
        kp1 = k + 1
        if kp1 >= pver:
            kp1 = pver - 1
        i = ilo
        while i <= ihi:
            idx = i + k * pcols
            utop = (wind0[idx] + wind0[i + km1 * pcols]) / 2.0
            vtop = (wind0[idx + plane] + wind0[i + km1 * pcols + plane]) / 2.0
            ubot = (wind0[i + kp1 * pcols] + wind0[idx]) / 2.0
            vbot = (wind0[i + kp1 * pcols + plane] + wind0[idx + plane]) / 2.0
            fket = utop * mflux[i + k * pcols] + vtop * mflux[i + k * pcols + flux_plane]
            fkeb = ubot * mflux[i + (k + 1) * pcols] + vbot * mflux[i + (k + 1) * pcols + flux_plane]
            ketend_cons = (fket - fkeb) / dp[idx]
            ketend = ((windf[idx] ** 2 + windf[idx + plane] ** 2) - (wind0[idx] ** 2 + wind0[idx + plane] ** 2)) * 0.5 / dt
            gseten[idx] = ketend_cons - ketend
            i += 1
        k += 1

    k = 0
    while k < pver:
        i = ilo
        while i <= ihi:
            ii = ideep[i] - 1
            seten[ii + k * pcols] = gseten[i + k * pcols]
            i += 1
        k += 1


@export
def zm_convtran_main_codon(
    pcols: int,
    pver: int,
    ncnst: int,
    il1g: int,
    il2g: int,
    doconvtran_p: cobj,
    is_dry_p: cobj,
    q_p: cobj,
    mu_p: cobj,
    md_p: cobj,
    du_p: cobj,
    eu_p: cobj,
    ed_p: cobj,
    dp_p: cobj,
    fracis_p: cobj,
    dpdry_p: cobj,
    jt_p: cobj,
    mx_p: cobj,
    ideep_p: cobj,
    dqdt_p: cobj,
    chat_p: cobj,
    cond_p: cobj,
    const_p: cobj,
    fisg_p: cobj,
    conu_p: cobj,
    dcondt_p: cobj,
    dutmp_p: cobj,
    eutmp_p: cobj,
    edtmp_p: cobj,
    dptmp_p: cobj,
):
    doconvtran = Ptr[int](doconvtran_p)
    is_dry = Ptr[int](is_dry_p)
    q = Ptr[float](q_p)
    mu = Ptr[float](mu_p)
    md = Ptr[float](md_p)
    du = Ptr[float](du_p)
    eu = Ptr[float](eu_p)
    ed = Ptr[float](ed_p)
    dp = Ptr[float](dp_p)
    fracis = Ptr[float](fracis_p)
    dpdry = Ptr[float](dpdry_p)
    jt = Ptr[int](jt_p)
    mx = Ptr[int](mx_p)
    ideep = Ptr[int](ideep_p)
    dqdt = Ptr[float](dqdt_p)
    chat = Ptr[float](chat_p)
    cond = Ptr[float](cond_p)
    const = Ptr[float](const_p)
    fisg = Ptr[float](fisg_p)
    conu = Ptr[float](conu_p)
    dcondt = Ptr[float](dcondt_p)
    dutmp = Ptr[float](dutmp_p)
    eutmp = Ptr[float](eutmp_p)
    edtmp = Ptr[float](edtmp_p)
    dptmp = Ptr[float](dptmp_p)
    plane = pcols * pver
    ilo = il1g - 1
    ihi = il2g - 1
    small = 1.0e-36
    mbsth = 1.0e-15

    ktm = pver
    kbm = pver
    i = ilo
    while i <= ihi:
        if jt[i] < ktm:
            ktm = jt[i]
        if mx[i] < kbm:
            kbm = mx[i]
        i += 1

    m = 1
    while m < ncnst:
        moff = m * plane
        if doconvtran[m] != 0:
            if is_dry[m] != 0:
                k = 0
                while k < pver:
                    i = ilo
                    while i <= ihi:
                        idx = i + k * pcols
                        dptmp[idx] = dpdry[idx]
                        dutmp[idx] = du[idx] * dp[idx] / dpdry[idx]
                        eutmp[idx] = eu[idx] * dp[idx] / dpdry[idx]
                        edtmp[idx] = ed[idx] * dp[idx] / dpdry[idx]
                        i += 1
                    k += 1
            else:
                k = 0
                while k < pver:
                    i = ilo
                    while i <= ihi:
                        idx = i + k * pcols
                        dptmp[idx] = dp[idx]
                        dutmp[idx] = du[idx]
                        eutmp[idx] = eu[idx]
                        edtmp[idx] = ed[idx]
                        i += 1
                    k += 1

            k = 0
            while k < pver:
                i = ilo
                while i <= ihi:
                    ii = ideep[i] - 1
                    idx = i + k * pcols
                    const[idx] = q[ii + k * pcols + moff]
                    fisg[idx] = fracis[ii + k * pcols + moff]
                    i += 1
                k += 1

            k = 0
            while k < pver:
                km1 = k - 1
                if km1 < 0:
                    km1 = 0
                i = ilo
                while i <= ihi:
                    idx = i + k * pcols
                    km1idx = i + km1 * pcols
                    minc = min(const[km1idx], const[idx])
                    maxc = max(const[km1idx], const[idx])
                    if minc < 0.0:
                        cdifr = 0.0
                    else:
                        cdifr = abs(const[idx] - const[km1idx]) / max(maxc, small)
                    if cdifr > 1.0e-6:
                        cabv = max(const[km1idx], maxc * 1.0e-12)
                        cbel = max(const[idx], maxc * 1.0e-12)
                        chat[idx] = log(cabv / cbel) / (cabv - cbel) * cabv * cbel
                    else:
                        chat[idx] = 0.5 * (const[idx] + const[km1idx])
                    conu[idx] = chat[idx]
                    cond[idx] = chat[idx]
                    dcondt[idx] = 0.0
                    i += 1
                k += 1

            k = 1
            km1 = 0
            kk = pver - 1
            i = ilo
            while i <= ihi:
                kkidx = i + kk * pcols
                mupdudp = mu[kkidx] + dutmp[kkidx] * dptmp[kkidx]
                if mupdudp > mbsth:
                    conu[kkidx] = (+eutmp[kkidx] * fisg[kkidx] * const[kkidx] * dptmp[kkidx]) / mupdudp
                idx = i + k * pcols
                km1idx = i + km1 * pcols
                if md[idx] < -mbsth:
                    cond[idx] = (-edtmp[km1idx] * fisg[km1idx] * const[km1idx] * dptmp[km1idx]) / md[idx]
                i += 1

            kk = pver - 2
            while kk >= 0:
                kkp1 = kk + 1
                i = ilo
                while i <= ihi:
                    idx = i + kk * pcols
                    mupdudp = mu[idx] + dutmp[idx] * dptmp[idx]
                    if mupdudp > mbsth:
                        conu[idx] = (
                            mu[i + kkp1 * pcols] * conu[i + kkp1 * pcols]
                            + eutmp[idx] * fisg[idx] * const[idx] * dptmp[idx]
                        ) / mupdudp
                    i += 1
                kk -= 1

            k = 2
            while k < pver:
                km1 = k - 1
                i = ilo
                while i <= ihi:
                    idx = i + k * pcols
                    km1idx = i + km1 * pcols
                    if md[idx] < -mbsth:
                        cond[idx] = (
                            md[km1idx] * cond[km1idx]
                            - edtmp[km1idx] * fisg[km1idx] * const[km1idx] * dptmp[km1idx]
                        ) / md[idx]
                    i += 1
                k += 1

            k = ktm - 1
            while k < pver:
                km1 = k - 1
                if km1 < 0:
                    km1 = 0
                kp1 = k + 1
                if kp1 >= pver:
                    kp1 = pver - 1
                i = ilo
                while i <= ihi:
                    idx = i + k * pcols
                    km1idx = i + km1 * pcols
                    kp1idx = i + kp1 * pcols
                    fluxin = mu[kp1idx] * conu[kp1idx] + mu[idx] * min(chat[idx], const[km1idx]) - (
                        md[idx] * cond[idx] + md[kp1idx] * min(chat[kp1idx], const[kp1idx])
                    )
                    fluxout = mu[idx] * conu[idx] + mu[kp1idx] * min(chat[kp1idx], const[idx]) - (
                        md[kp1idx] * cond[kp1idx] + md[idx] * min(chat[idx], const[idx])
                    )
                    netflux = fluxin - fluxout
                    if abs(netflux) < max(fluxin, fluxout) * 1.0e-12:
                        netflux = 0.0
                    dcondt[idx] = netflux / dptmp[idx]
                    i += 1
                k += 1

            k = kbm - 1
            while k < pver:
                km1 = k - 1
                if km1 < 0:
                    km1 = 0
                i = ilo
                while i <= ihi:
                    idx = i + k * pcols
                    km1idx = i + km1 * pcols
                    if k + 1 == mx[i]:
                        fluxin = mu[idx] * min(chat[idx], const[km1idx]) - md[idx] * cond[idx]
                        fluxout = mu[idx] * conu[idx] - md[idx] * min(chat[idx], const[idx])
                        netflux = fluxin - fluxout
                        if abs(netflux) < max(fluxin, fluxout) * 1.0e-12:
                            netflux = 0.0
                        dcondt[idx] = netflux / dptmp[idx]
                    elif k + 1 > mx[i]:
                        dcondt[idx] = 0.0
                    i += 1
                k += 1

            k = 0
            while k < pver:
                i = 0
                while i < pcols:
                    dqdt[i + k * pcols + moff] = 0.0
                    i += 1
                k += 1

            k = 0
            while k < pver:
                i = ilo
                while i <= ihi:
                    ii = ideep[i] - 1
                    dqdt[ii + k * pcols + moff] = dcondt[i + k * pcols]
                    i += 1
                k += 1
        m += 1


@export
def zm_convtran1_prep_shell_codon(
    pcols: int,
    pver: int,
    fake_dpdry_p: cobj,
):
    fake_dpdry = Ptr[float](fake_dpdry_p)

    k = 0
    while k < pver:
        i = 0
        while i < pcols:
            fake_dpdry[i + k * pcols] = 0.0
            i += 1
        k += 1


@export
def zm_convtran1_ratio_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    wtrc_nwset: int,
    Rwt_p: cobj,
    ptend_q_p: cobj,
    liq_type_p: cobj,
    ice_type_p: cobj,
):
    Rwt = Ptr[float](Rwt_p)
    ptend_q = Ptr[float](ptend_q_p)
    liq_type = Ptr[int](liq_type_p)
    ice_type = Ptr[int](ice_type_p)
    plane = pcols * pver

    m = 1
    while m < wtrc_nwset:
        liq_idx = (liq_type[m] - 1) * plane
        ice_idx = (ice_type[m] - 1) * plane
        base_liq_idx = (liq_type[0] - 1) * plane
        base_ice_idx = (ice_type[0] - 1) * plane
        ratio_m_offset = m * plane
        ratio_ice_offset = wtrc_nwset * plane + m * plane
        k = 0
        while k < pver:
            i = 0
            while i < ncol:
                idx = i + k * pcols
                ptend_q[idx + liq_idx] = Rwt[idx + ratio_m_offset] * ptend_q[idx + base_liq_idx]
                ptend_q[idx + ice_idx] = Rwt[idx + ratio_ice_offset] * ptend_q[idx + base_ice_idx]
                i += 1
            k += 1
        m += 1


@export
def zm_convtran2_dpdry_shell_codon(
    pcols: int,
    pver: int,
    lengath: int,
    state_pdeldry_p: cobj,
    ideep_p: cobj,
    dpdry_p: cobj,
):
    state_pdeldry = Ptr[float](state_pdeldry_p)
    ideep = Ptr[int](ideep_p)
    dpdry = Ptr[float](dpdry_p)

    k = 0
    while k < pver:
        i = 0
        while i < pcols:
            dpdry[i + k * pcols] = 0.0
            i += 1
        k += 1

    i = 0
    while i < lengath:
        src_i = ideep[i] - 1
        k = 0
        while k < pver:
            dpdry[i + k * pcols] = state_pdeldry[src_i + k * pcols] / 100.0
            k += 1
        i += 1
