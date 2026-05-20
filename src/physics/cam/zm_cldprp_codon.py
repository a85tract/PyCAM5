@inline
def _idx1(i: int) -> int:
    return i - 1


@inline
def _idx2(i: int, k: int, pcols: int) -> int:
    """cldprp arrays declared as (pcols,pver) or (pcols,pverp)."""
    return (i - 1) + (k - 1) * pcols


def zm_cldprp_init_arrays_codon(
    il2g: int,
    pcols: int,
    pver: int,
    c0_ocn: float,
    c0_lnd: float,
    landfrac_p: cobj,
    zf_p: cobj,
    q_p: cobj,
    s_p: cobj,
    ftemp_p: cobj,
    expnum_p: cobj,
    expdif_p: cobj,
    c0mask_p: cobj,
    dz_p: cobj,
    pflx_p: cobj,
    k1_p: cobj,
    i2_p: cobj,
    i3_p: cobj,
    i4_p: cobj,
    mu_p: cobj,
    f_p: cobj,
    eps_p: cobj,
    eu_p: cobj,
    du_p: cobj,
    ql_p: cobj,
    cu_p: cobj,
    evp_p: cobj,
    wtevp_p: cobj,
    cmeg_p: cobj,
    qds_p: cobj,
    md_p: cobj,
    ed_p: cobj,
    sd_p: cobj,
    qd_p: cobj,
    mc_p: cobj,
    qu_p: cobj,
    su_p: cobj,
    rprd_p: cobj,
    totpcp_p: cobj,
    totevp_p: cobj,
):
    landfrac = Ptr[float](landfrac_p)
    zf = Ptr[float](zf_p)
    q = Ptr[float](q_p)
    s = Ptr[float](s_p)
    ftemp = Ptr[float](ftemp_p)
    expnum = Ptr[float](expnum_p)
    expdif = Ptr[float](expdif_p)
    c0mask = Ptr[float](c0mask_p)
    dz = Ptr[float](dz_p)
    pflx = Ptr[float](pflx_p)
    k1 = Ptr[float](k1_p)
    i2 = Ptr[float](i2_p)
    i3 = Ptr[float](i3_p)
    i4 = Ptr[float](i4_p)
    mu = Ptr[float](mu_p)
    f = Ptr[float](f_p)
    eps = Ptr[float](eps_p)
    eu = Ptr[float](eu_p)
    du = Ptr[float](du_p)
    ql = Ptr[float](ql_p)
    cu = Ptr[float](cu_p)
    evp = Ptr[float](evp_p)
    wtevp = Ptr[float](wtevp_p)
    cmeg = Ptr[float](cmeg_p)
    qds = Ptr[float](qds_p)
    md = Ptr[float](md_p)
    ed = Ptr[float](ed_p)
    sd = Ptr[float](sd_p)
    qd = Ptr[float](qd_p)
    mc = Ptr[float](mc_p)
    qu = Ptr[float](qu_p)
    su = Ptr[float](su_p)
    rprd = Ptr[float](rprd_p)
    totpcp = Ptr[float](totpcp_p)
    totevp = Ptr[float](totevp_p)

    for i in range(1, il2g + 1):
        idx = _idx1(i)
        ftemp[idx] = 0.0
        expnum[idx] = 0.0
        expdif[idx] = 0.0
        c0mask[idx] = c0_ocn * (1.0 - landfrac[idx]) + c0_lnd * landfrac[idx]
        totpcp[idx] = 0.0
        totevp[idx] = 0.0
        pflx[_idx2(i, 1, pcols)] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, il2g + 1):
            idx = _idx2(i, k, pcols)
            dz[idx] = zf[_idx2(i, k, pcols)] - zf[_idx2(i, k + 1, pcols)]
            k1[idx] = 0.0
            i2[idx] = 0.0
            i3[idx] = 0.0
            i4[idx] = 0.0
            mu[idx] = 0.0
            f[idx] = 0.0
            eps[idx] = 0.0
            eu[idx] = 0.0
            du[idx] = 0.0
            ql[idx] = 0.0
            cu[idx] = 0.0
            evp[idx] = 0.0
            wtevp[idx] = 0.0
            cmeg[idx] = 0.0
            qds[idx] = q[idx]
            md[idx] = 0.0
            ed[idx] = 0.0
            sd[idx] = s[idx]
            qd[idx] = q[idx]
            mc[idx] = 0.0
            qu[idx] = q[idx]
            su[idx] = s[idx]
            rprd[idx] = 0.0


def zm_cldprp_index_setup_codon(
    il2g: int,
    pcols: int,
    pver: int,
    msg: int,
    limcnv: int,
    cp_tiedke_add: float,
    tiedke_add: float,
    jb_p: cobj,
    lel_p: cobj,
    mx_p: cobj,
    hsat_p: cobj,
    hmn_p: cobj,
    s_p: cobj,
    jt_p: cobj,
    jd_p: cobj,
    jlcl_p: cobj,
    hmin_p: cobj,
    j0_p: cobj,
    hu_p: cobj,
    su_p: cobj,
):
    jb = Ptr[i32](jb_p)
    lel = Ptr[i32](lel_p)
    mx = Ptr[i32](mx_p)
    hsat = Ptr[float](hsat_p)
    hmn = Ptr[float](hmn_p)
    s = Ptr[float](s_p)
    jt = Ptr[i32](jt_p)
    jd = Ptr[i32](jd_p)
    jlcl = Ptr[i32](jlcl_p)
    hmin = Ptr[float](hmin_p)
    j0 = Ptr[i32](j0_p)
    hu = Ptr[float](hu_p)
    su = Ptr[float](su_p)

    for i in range(1, pcols + 1):
        jt[_idx1(i)] = i32(pver)

    for i in range(1, il2g + 1):
        idx = _idx1(i)
        top = int(lel[idx])
        if limcnv + 1 > top:
            top = limcnv + 1
        if top > pver:
            top = pver
        jt[idx] = i32(top)
        jd[idx] = i32(pver)
        jlcl[idx] = lel[idx]
        hmin[idx] = 1.0e6

    for k in range(msg + 1, pver + 1):
        for i in range(1, il2g + 1):
            i0 = _idx1(i)
            idx = _idx2(i, k, pcols)
            if hsat[idx] <= hmin[i0] and k >= int(jt[i0]) and k <= int(jb[i0]):
                hmin[i0] = hsat[idx]
                j0[i0] = i32(k)

    for i in range(1, il2g + 1):
        i0 = _idx1(i)
        if int(j0[i0]) > int(jb[i0]) - 2:
            j0[i0] = i32(int(jb[i0]) - 2)
        if int(j0[i0]) < int(jt[i0]) + 2:
            j0[i0] = i32(int(jt[i0]) + 2)
        if int(j0[i0]) > pver:
            j0[i0] = i32(pver)

    for k in range(msg + 1, pver + 1):
        for i in range(1, il2g + 1):
            i0 = _idx1(i)
            if k >= int(jt[i0]) and k <= int(jb[i0]):
                idx = _idx2(i, k, pcols)
                mx_idx = _idx2(i, int(mx[i0]), pcols)
                hu[idx] = hmn[mx_idx] + cp_tiedke_add
                su[idx] = s[mx_idx] + tiedke_add


def zm_cldprp_copy_mass_fields_codon(
    pcols: int,
    pver: int,
    ed_p: cobj,
    md_p: cobj,
    mu_p: cobj,
    du_p: cobj,
    eu_p: cobj,
    wted_p: cobj,
    wtmd_p: cobj,
    wtmu_p: cobj,
    wtdu_p: cobj,
    wteu_p: cobj,
):
    ed = Ptr[float](ed_p)
    md = Ptr[float](md_p)
    mu = Ptr[float](mu_p)
    du = Ptr[float](du_p)
    eu = Ptr[float](eu_p)
    wted = Ptr[float](wted_p)
    wtmd = Ptr[float](wtmd_p)
    wtmu = Ptr[float](wtmu_p)
    wtdu = Ptr[float](wtdu_p)
    wteu = Ptr[float](wteu_p)

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            idx = _idx2(i, k, pcols)
            wted[idx] = ed[idx]
            wtmd[idx] = md[idx]
            wtmu[idx] = mu[idx]
            wtdu[idx] = du[idx]
            wteu[idx] = eu[idx]


def zm_cldprp_copy_2d_codon(
    pcols: int,
    pver: int,
    src_p: cobj,
    dst_p: cobj,
):
    src = Ptr[float](src_p)
    dst = Ptr[float](dst_p)

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            idx = _idx2(i, k, pcols)
            dst[idx] = src[idx]


def zm_cldprp_evap_finalize_codon(
    il2g: int,
    pcols: int,
    pver: int,
    pverp: int,
    msg: int,
    jd_p: cobj,
    jb_p: cobj,
    md_p: cobj,
    qd_p: cobj,
    totevp_p: cobj,
    totpcp_p: cobj,
    ed_p: cobj,
    evp_p: cobj,
    cu_p: cobj,
    cmeg_p: cobj,
    rprd_p: cobj,
    dz_p: cobj,
    pflx_p: cobj,
    mc_p: cobj,
    mu_p: cobj,
    wtevp_p: cobj,
):
    jd = Ptr[i32](jd_p)
    jb = Ptr[i32](jb_p)
    md = Ptr[float](md_p)
    qd = Ptr[float](qd_p)
    totevp = Ptr[float](totevp_p)
    totpcp = Ptr[float](totpcp_p)
    ed = Ptr[float](ed_p)
    evp = Ptr[float](evp_p)
    cu = Ptr[float](cu_p)
    cmeg = Ptr[float](cmeg_p)
    rprd = Ptr[float](rprd_p)
    dz = Ptr[float](dz_p)
    pflx = Ptr[float](pflx_p)
    mc = Ptr[float](mc_p)
    mu = Ptr[float](mu_p)
    wtevp = Ptr[float](wtevp_p)

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            idx = _idx2(i, k, pcols)
            wtevp[idx] = evp[idx]

    for i in range(1, il2g + 1):
        i0 = _idx1(i)
        jd_i = int(jd[i0])
        jb_i = int(jb[i0])
        totevp[i0] = (
            totevp[i0]
            + md[_idx2(i, jd_i, pcols)] * qd[_idx2(i, jd_i, pcols)]
            - md[_idx2(i, jb_i, pcols)] * qd[_idx2(i, jb_i, pcols)]
        )

    for i in range(1, il2g + 1):
        i0 = _idx1(i)
        totpcp[i0] = max(totpcp[i0], 0.0)
        totevp[i0] = max(totevp[i0], 0.0)

    for k in range(msg + 2, pver + 1):
        for i in range(1, il2g + 1):
            i0 = _idx1(i)
            idx = _idx2(i, k, pcols)
            if totevp[i0] > 0.0 and totpcp[i0] > 0.0:
                md[idx] = md[idx] * min(1.0, totpcp[i0] / (totevp[i0] + totpcp[i0]))
                ed[idx] = ed[idx] * min(1.0, totpcp[i0] / (totevp[i0] + totpcp[i0]))
                evp[idx] = evp[idx] * min(1.0, totpcp[i0] / (totevp[i0] + totpcp[i0]))
            else:
                md[idx] = 0.0
                ed[idx] = 0.0
                evp[idx] = 0.0
            cmeg[idx] = cu[idx] - evp[idx]
            rprd[idx] = rprd[idx] - evp[idx]

    for i in range(1, il2g + 1):
        pflx[_idx2(i, 1, pcols)] = 0.0

    for k in range(2, pverp + 1):
        for i in range(1, il2g + 1):
            pflx[_idx2(i, k, pcols)] = (
                pflx[_idx2(i, k - 1, pcols)]
                + rprd[_idx2(i, k - 1, pcols)] * dz[_idx2(i, k - 1, pcols)]
            )

    for k in range(msg + 1, pver + 1):
        for i in range(1, il2g + 1):
            idx = _idx2(i, k, pcols)
            mc[idx] = mu[idx] + md[idx]
