@inline
def _field_idx(icol: int, klev: int, ld1: int) -> int:
    """qrl/qrs declared as (pcols, pver); ptend%s declared as (psetcols, pver)"""
    return (icol - 1) + (klev - 1) * ld1


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
def radiation_diag_prep_codon(
    stage: int,
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
