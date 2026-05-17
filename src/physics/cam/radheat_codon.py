from math import exp, sqrt


@inline
def _field_idx(icol: int, klev: int, ld1: int) -> int:
    """qrl/qrs declared as (pcols, pver); ptend%s declared as (psetcols, pver)"""
    return (icol - 1) + (klev - 1) * ld1


@inline
def _band3_idx(ibnd: int, icol: int, klev: int, nbnd: int, pcols: int) -> int:
    """cloud optics arrays declared as (nbnd, pcols, pver)"""
    return (ibnd - 1) + (icol - 1) * nbnd + (klev - 1) * nbnd * pcols


@inline
def _col_idx(icol: int) -> int:
    """fsns/fsnt/flns/flnt/net_flx declared as (pcols)"""
    return icol - 1


@export
def radheat_timestep_init_codon():
    return


@export
def radheat_batch_timestep_init_codon():
    radheat_timestep_init_codon()


@export
def radheat_batch_timestep_init_stage_dispatch_codon():
    radheat_batch_timestep_init_codon()


@export
def radheat_tend_codon(
    ncol: int,
    pcols: int,
    pver: int,
    psetcols: int,
    qrl_p: cobj,
    qrs_p: cobj,
    ptend_s_p: cobj,
    fsns_p: cobj,
    fsnt_p: cobj,
    flns_p: cobj,
    flnt_p: cobj,
    net_flx_p: cobj,
):
    qrl = Ptr[float](qrl_p)
    qrs = Ptr[float](qrs_p)
    ptend_s = Ptr[float](ptend_s_p)
    fsns = Ptr[float](fsns_p)
    fsnt = Ptr[float](fsnt_p)
    flns = Ptr[float](flns_p)
    flnt = Ptr[float](flnt_p)
    net_flx = Ptr[float](net_flx_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            ptend_s[_field_idx(i, k, psetcols)] = (
                qrs[_field_idx(i, k, pcols)] + qrl[_field_idx(i, k, pcols)]
            )

    for i in range(1, ncol + 1):
        net_flx[_col_idx(i)] = (
            fsnt[_col_idx(i)]
            - fsns[_col_idx(i)]
            - flnt[_col_idx(i)]
            + flns[_col_idx(i)]
        )


@export
def radheat_batch_tend_codon(
    ncol: int,
    pcols: int,
    pver: int,
    psetcols: int,
    qrl_p: cobj,
    qrs_p: cobj,
    ptend_s_p: cobj,
    fsns_p: cobj,
    fsnt_p: cobj,
    flns_p: cobj,
    flnt_p: cobj,
    net_flx_p: cobj,
):
    radheat_tend_codon(
        ncol,
        pcols,
        pver,
        psetcols,
        qrl_p,
        qrs_p,
        ptend_s_p,
        fsns_p,
        fsnt_p,
        flns_p,
        flnt_p,
        net_flx_p,
    )


@export
def radheat_batch_tend_stage_dispatch_codon(
    ncol: int,
    pcols: int,
    pver: int,
    psetcols: int,
    qrl_p: cobj,
    qrs_p: cobj,
    ptend_s_p: cobj,
    fsns_p: cobj,
    fsnt_p: cobj,
    flns_p: cobj,
    flnt_p: cobj,
    net_flx_p: cobj,
):
    radheat_batch_tend_codon(
        ncol,
        pcols,
        pver,
        psetcols,
        qrl_p,
        qrs_p,
        ptend_s_p,
        fsns_p,
        fsnt_p,
        flns_p,
        flnt_p,
        net_flx_p,
    )


@export
def radiation_diag_prep_codon(
    stage: int,
    mode: int,
    ncol: int,
    pcols: int,
    pver: int,
    nday: int,
    nnite: int,
    cpair: float,
    cgs2mks: float,
    fillvalue: float,
    scalar_p: cobj,
    a_p: cobj,
    b_p: cobj,
    c_p: cobj,
    d_p: cobj,
    e_p: cobj,
    f_p: cobj,
    g_p: cobj,
    h_p: cobj,
    i_p: cobj,
    j_p: cobj,
    k_p: cobj,
    l_p: cobj,
    m_p: cobj,
    n_p: cobj,
    o_p: cobj,
    p_p: cobj,
    q_p: cobj,
    r_p: cobj,
    idxday_p: cobj,
    idxnite_p: cobj,
):
    scalar = Ptr[float](scalar_p)
    a = Ptr[float](a_p)
    b = Ptr[float](b_p)
    c = Ptr[float](c_p)
    d = Ptr[float](d_p)
    e = Ptr[float](e_p)
    f = Ptr[float](f_p)
    g = Ptr[float](g_p)
    h = Ptr[float](h_p)
    iarr = Ptr[float](i_p)
    jarr = Ptr[float](j_p)
    karr = Ptr[float](k_p)
    larr = Ptr[float](l_p)
    marr = Ptr[float](m_p)
    narr = Ptr[float](n_p)
    oarr = Ptr[float](o_p)
    parr = Ptr[float](p_p)
    qarr = Ptr[float](q_p)
    rarr = Ptr[float](r_p)
    idxday = Ptr[int](idxday_p)
    idxnite = Ptr[int](idxnite_p)

    if stage == 1:
        # SW CGS-to-MKS conversion plus shortwave cloud forcing.
        for i in range(1, ncol + 1):
            a[_col_idx(i)] = a[_col_idx(i)] * cgs2mks
            b[_col_idx(i)] = b[_col_idx(i)] * cgs2mks
            c[_col_idx(i)] = c[_col_idx(i)] * cgs2mks
            d[_col_idx(i)] = d[_col_idx(i)] * cgs2mks
            e[_col_idx(i)] = e[_col_idx(i)] * cgs2mks
            f[_col_idx(i)] = f[_col_idx(i)] * cgs2mks
            g[_col_idx(i)] = g[_col_idx(i)] * cgs2mks
            h[_col_idx(i)] = h[_col_idx(i)] * cgs2mks
            iarr[_col_idx(i)] = iarr[_col_idx(i)] * cgs2mks
            jarr[_col_idx(i)] = jarr[_col_idx(i)] * cgs2mks
            karr[_col_idx(i)] = karr[_col_idx(i)] * cgs2mks
            larr[_col_idx(i)] = larr[_col_idx(i)] * cgs2mks
            marr[_col_idx(i)] = marr[_col_idx(i)] * cgs2mks
            narr[_col_idx(i)] = narr[_col_idx(i)] * cgs2mks
            oarr[_col_idx(i)] = oarr[_col_idx(i)] * cgs2mks
            parr[_col_idx(i)] = parr[_col_idx(i)] * cgs2mks
            qarr[_col_idx(i)] = qarr[_col_idx(i)] * cgs2mks
            rarr[_col_idx(i)] = larr[_col_idx(i)] - marr[_col_idx(i)]
    elif stage == 2:
        # ftem = qrs/cpair or qrsc/cpair, depending on caller-provided input.
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                b[_field_idx(i, k, pcols)] = a[_field_idx(i, k, pcols)] / cpair
    elif stage == 3:
        # Visible cloud optical-depth diagnostics.
        for k in range(1, pver + 1):
            for i in range(1, pcols + 1):
                a[_field_idx(i, k, pcols)] = fillvalue
                b[_field_idx(i, k, pcols)] = fillvalue

        for ii in range(1, nday + 1):
            col = idxday[_col_idx(ii)]
            for k in range(1, pver + 1):
                total = c[_field_idx(col, k, pcols)] + d[_field_idx(col, k, pcols)]
                b[_field_idx(col, k, pcols)] = total
                a[_field_idx(col, k, pcols)] = total * e[_field_idx(col, k, pcols)]

        for ii in range(1, nnite + 1):
            col = idxnite[_col_idx(ii)]
            for k in range(1, pver + 1):
                c[_field_idx(col, k, pcols)] = fillvalue
                d[_field_idx(col, k, pcols)] = fillvalue
    elif stage == 4:
        # LW CGS-to-MKS conversion plus longwave cloud forcing.
        for i in range(1, ncol + 1):
            a[_col_idx(i)] = a[_col_idx(i)] * cgs2mks
            b[_col_idx(i)] = b[_col_idx(i)] * cgs2mks
            c[_col_idx(i)] = c[_col_idx(i)] * cgs2mks
            d[_col_idx(i)] = d[_col_idx(i)] * cgs2mks
            e[_col_idx(i)] = e[_col_idx(i)] * cgs2mks
            f[_col_idx(i)] = g[_col_idx(i)] * cgs2mks
            h[_col_idx(i)] = h[_col_idx(i)] * cgs2mks
            iarr[_col_idx(i)] = iarr[_col_idx(i)] * cgs2mks
            jarr[_col_idx(i)] = jarr[_col_idx(i)] * cgs2mks
            g[_col_idx(i)] = g[_col_idx(i)] * cgs2mks
            karr[_col_idx(i)] = c[_col_idx(i)] - b[_col_idx(i)]
            larr[_col_idx(i)] = larr[_col_idx(i)] * cgs2mks
    elif stage == 5:
        # Convert qrs/qrl between Q and Q*dp for energy conservation.
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                if nday == 1:
                    a[_field_idx(i, k, pcols)] = a[_field_idx(i, k, pcols)] / b[_field_idx(i, k, pcols)]
                    c[_field_idx(i, k, pcols)] = c[_field_idx(i, k, pcols)] / b[_field_idx(i, k, pcols)]
                else:
                    a[_field_idx(i, k, pcols)] = a[_field_idx(i, k, pcols)] * b[_field_idx(i, k, pcols)]
                    c[_field_idx(i, k, pcols)] = c[_field_idx(i, k, pcols)] * b[_field_idx(i, k, pcols)]
    elif stage == 6:
        # Cloud forcing diagnostics: out = all-sky clear/up flux difference.
        for i in range(1, ncol + 1):
            c[_col_idx(i)] = a[_col_idx(i)] - b[_col_idx(i)]
    elif stage == 7:
        # radinp pressure conversion: pmid/pint Pa -> dynes/cm2.
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                c[_field_idx(i, k, pcols)] = a[_field_idx(i, k, pcols)] * cpair
                d[_field_idx(i, k, pcols)] = b[_field_idx(i, k, pcols)] * cpair

        for i in range(1, ncol + 1):
            d[_field_idx(i, pver + 1, pcols)] = b[_field_idx(i, pver + 1, pcols)] * cpair
    elif stage == 8:
        # calc_col_mean: preserve Fortran k-outer accumulation order.
        for i in range(1, pcols + 1):
            c[_col_idx(i)] = 0.0
            d[_col_idx(i)] = 0.0

        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                c[_col_idx(i)] += a[_field_idx(i, k, pcols)] * b[_field_idx(i, k, pcols)]
                d[_col_idx(i)] += b[_field_idx(i, k, pcols)]

        for i in range(1, ncol + 1):
            c[_col_idx(i)] = c[_col_idx(i)] / d[_col_idx(i)]
    elif stage == 9:
        # SW cloud/snow optics combine. a=cld, b=cldfsnow, c-f=cld_tau*, g-j=snow_tau*,
        # k-n=c_cld_tau*, o=cldfprime. nday carries nbndsw, nnite is has_snow.
        nbnd = nday
        has_snow = nnite != 0
        if has_snow:
            for i in range(1, ncol + 1):
                for kk in range(1, pver + 1):
                    idx2 = _field_idx(i, kk, pcols)
                    cldval = a[idx2]
                    snowval = b[idx2]
                    prime = cldval
                    if snowval > prime:
                        prime = snowval
                    oarr[idx2] = prime
                    if prime > 0.0:
                        for ib in range(1, nbnd + 1):
                            idx3 = _band3_idx(ib, i, kk, nbnd, pcols)
                            karr[idx3] = (snowval * g[idx3] + cldval * c[idx3]) / prime
                            larr[idx3] = (snowval * h[idx3] + cldval * d[idx3]) / prime
                            marr[idx3] = (snowval * iarr[idx3] + cldval * e[idx3]) / prime
                            narr[idx3] = (snowval * jarr[idx3] + cldval * f[idx3]) / prime
                    else:
                        for ib in range(1, nbnd + 1):
                            idx3 = _band3_idx(ib, i, kk, nbnd, pcols)
                            karr[idx3] = 0.0
                            larr[idx3] = 0.0
                            marr[idx3] = 0.0
                            narr[idx3] = 0.0
        else:
            for i in range(1, ncol + 1):
                for kk in range(1, pver + 1):
                    idx2 = _field_idx(i, kk, pcols)
                    oarr[idx2] = a[idx2]
                    for ib in range(1, nbnd + 1):
                        idx3 = _band3_idx(ib, i, kk, nbnd, pcols)
                        karr[idx3] = c[idx3]
                        larr[idx3] = d[idx3]
                        marr[idx3] = e[idx3]
                        narr[idx3] = f[idx3]
    elif stage == 10:
        # LW cloud/snow absorption combine. a=cld, b=cldfsnow, c=cld_lw_abs,
        # d=snow_lw_abs, e=c_cld_lw_abs, o=cldfprime. nday carries nbndlw.
        nbnd = nday
        has_snow = nnite != 0
        if has_snow:
            for i in range(1, ncol + 1):
                for kk in range(1, pver + 1):
                    idx2 = _field_idx(i, kk, pcols)
                    cldval = a[idx2]
                    snowval = b[idx2]
                    prime = cldval
                    if snowval > prime:
                        prime = snowval
                    oarr[idx2] = prime
                    if prime > 0.0:
                        for ib in range(1, nbnd + 1):
                            idx3 = _band3_idx(ib, i, kk, nbnd, pcols)
                            e[idx3] = (snowval * d[idx3] + cldval * c[idx3]) / prime
                    else:
                        for ib in range(1, nbnd + 1):
                            e[_band3_idx(ib, i, kk, nbnd, pcols)] = 0.0
        else:
            for i in range(1, ncol + 1):
                for kk in range(1, pver + 1):
                    idx2 = _field_idx(i, kk, pcols)
                    oarr[idx2] = a[idx2]
                    for ib in range(1, nbnd + 1):
                        idx3 = _band3_idx(ib, i, kk, nbnd, pcols)
                        e[idx3] = c[idx3]
    elif stage == 11:
        # COSP snow diagnostic band extract. a=cldfsnow, b=band3 input, c=output2D.
        # nday carries the band3 leading dimension; nnite carries selected band or 0.
        nbnd = nday
        selected_band = nnite
        for kk in range(1, pver + 1):
            for i in range(1, pcols + 1):
                c[_field_idx(i, kk, pcols)] = 0.0

        if selected_band != 0:
            for i in range(1, ncol + 1):
                for kk in range(1, pver + 1):
                    idx2 = _field_idx(i, kk, pcols)
                    snowfrac = a[idx2]
                    if snowfrac > 0.0:
                        c[idx2] = b[_band3_idx(selected_band, i, kk, nbnd, pcols)] * snowfrac
    elif stage == 12:
        # Visible optical-depth diagnostics. a=c_cld_tau, b=liq_tau, c=ice_tau,
        # d=snow_tau, e=cldfprime, f=tot_cld, g=tot_icld, h=liq_icld,
        # i=ice_icld, j=snow_icld. abs(nday)=nbndsw, nday>0 when snow is active,
        # nnite=selected band, mode=Nnite, idxnite is an int64 copy of IdxNite.
        nbnd = nday
        if nbnd < 0:
            nbnd = -nbnd
        selected_band = nnite
        has_snow = nday > 0
        for kk in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx2 = _field_idx(i, kk, pcols)
                band_idx = _band3_idx(selected_band, i, kk, nbnd, pcols)
                g[idx2] = a[band_idx]
                h[idx2] = b[band_idx]
                iarr[idx2] = c[band_idx]
                if has_snow:
                    jarr[idx2] = d[band_idx]
                f[idx2] = a[band_idx] * e[idx2]

        for ii in range(1, mode + 1):
            col = idxnite[_col_idx(ii)]
            for kk in range(1, pver + 1):
                idx2 = _field_idx(col, kk, pcols)
                f[idx2] = fillvalue
                g[idx2] = fillvalue
                h[idx2] = fillvalue
                iarr[idx2] = fillvalue
                if has_snow:
                    jarr[idx2] = fillvalue
    elif stage == 13:
        # SW cloud optics sum. a-d=liq_tau*, e-h=ice_tau*, i-l=cld_tau*.
        nbnd = nday
        for kk in range(1, pver + 1):
            for i in range(1, ncol + 1):
                for ib in range(1, nbnd + 1):
                    idx3 = _band3_idx(ib, i, kk, nbnd, pcols)
                    iarr[idx3] = a[idx3] + e[idx3]
                    jarr[idx3] = b[idx3] + f[idx3]
                    karr[idx3] = c[idx3] + g[idx3]
                    larr[idx3] = d[idx3] + h[idx3]
    elif stage == 14:
        # LW cloud absorption sum. a=liq_lw_abs, b=ice_lw_abs, c=cld_lw_abs.
        nbnd = nday
        for kk in range(1, pver + 1):
            for i in range(1, ncol + 1):
                for ib in range(1, nbnd + 1):
                    idx3 = _band3_idx(ib, i, kk, nbnd, pcols)
                    c[idx3] = a[idx3] + b[idx3]
    elif stage == 15:
        # Compact longwave history workspace: input has pcols leading dimension,
        # output matches the original compiler temporary shape (ncol, pver).
        out_ld = nday
        for kk in range(1, pver + 1):
            for i in range(1, ncol + 1):
                b[_field_idx(i, kk, out_ld)] = a[_field_idx(i, kk, pcols)] / cpair
    elif stage == 16:
        # EMIS diagnostic: a=cld_lw_abs(nbndlw,pcols,pver), b=emis(pcols,pver).
        nbnd = nday
        selected_band = nnite
        for kk in range(1, pver + 1):
            for i in range(1, pcols + 1):
                b[_field_idx(i, kk, pcols)] = 0.0

        for kk in range(1, pver + 1):
            for i in range(1, ncol + 1):
                b[_field_idx(i, kk, pcols)] = 1.0 - exp(
                    -a[_band3_idx(selected_band, i, kk, nbnd, pcols)]
                )
    elif stage == 17:
        # Interface temperature prep. a=state%t, b=state%lnpint,
        # c=state%lnpmid, d=cam_in%lwup, e=tint.
        for i in range(1, ncol + 1):
            e[_field_idx(i, 1, pcols)] = a[_field_idx(i, 1, pcols)]
            e[_field_idx(i, pver + 1, pcols)] = sqrt(sqrt(d[_col_idx(i)] / cpair))
            for kk in range(2, pver + 1):
                dy = (b[_field_idx(i, kk, pcols)] - c[_field_idx(i, kk, pcols)]) / (
                    c[_field_idx(i, kk - 1, pcols)] - c[_field_idx(i, kk, pcols)]
                )
                e[_field_idx(i, kk, pcols)] = a[_field_idx(i, kk, pcols)] - dy * (
                    a[_field_idx(i, kk, pcols)] - a[_field_idx(i, kk - 1, pcols)]
                )
    elif stage == 18:
        # HIRS prep. a=cam_in%lwup, b=landfrac, c=state%pint,
        # d=ts, e=oro, f=pintmb.
        for i in range(1, ncol + 1):
            d[_col_idx(i)] = sqrt(sqrt(a[_col_idx(i)] / cpair))
            if b[_col_idx(i)] >= 0.001:
                e[_col_idx(i)] = 1.0
            else:
                e[_col_idx(i)] = 0.0
            for kk in range(1, pver + 1):
                f[_field_idx(i, kk, pcols)] = c[_field_idx(i, kk, pcols)] * 1.0e-2
            f[_field_idx(i, pver + 1, pcols)] = c[_field_idx(i, pver + 1, pcols)] * 1.0e-2
    elif stage == 19:
        # Longwave upward flux CGS conversion. a=cam_in%lwup, b=lwupcgs.
        for i in range(1, ncol + 1):
            b[_col_idx(i)] = a[_col_idx(i)] * 1000.0
    elif stage == 20:
        # Final net shortwave export copy. a=fsns, b=cam_out%netsw.
        for i in range(1, ncol + 1):
            b[_col_idx(i)] = a[_col_idx(i)]
    elif stage == 21:
        # HR diagnostic workspace. a=qrs, b=qrl, c=native-computed theta factor, d=ftem.
        for kk in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx2 = _field_idx(i, kk, pcols)
                d[idx2] = (a[idx2] + b[idx2]) / cpair * c[idx2]


@export
def radiation_diag_prep_stage_dispatch_codon(
    stage: int,
    mode: int,
    ncol: int,
    pcols: int,
    pver: int,
    nday: int,
    nnite: int,
    cpair: float,
    cgs2mks: float,
    fillvalue: float,
    scalar_p: cobj,
    a_p: cobj,
    b_p: cobj,
    c_p: cobj,
    d_p: cobj,
    e_p: cobj,
    f_p: cobj,
    g_p: cobj,
    h_p: cobj,
    i_p: cobj,
    j_p: cobj,
    k_p: cobj,
    l_p: cobj,
    m_p: cobj,
    n_p: cobj,
    o_p: cobj,
    p_p: cobj,
    q_p: cobj,
    r_p: cobj,
    idxday_p: cobj,
    idxnite_p: cobj,
):
    radiation_diag_prep_codon(
        stage,
        mode,
        ncol,
        pcols,
        pver,
        nday,
        nnite,
        cpair,
        cgs2mks,
        fillvalue,
        scalar_p,
        a_p,
        b_p,
        c_p,
        d_p,
        e_p,
        f_p,
        g_p,
        h_p,
        i_p,
        j_p,
        k_p,
        l_p,
        m_p,
        n_p,
        o_p,
        p_p,
        q_p,
        r_p,
        idxday_p,
        idxnite_p,
    )
