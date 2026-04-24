from math import acos, cos, exp, log, sin, sqrt


@inline
def _idx2(i: int, k: int, ld1: int) -> int:
    return (i - 1) + (k - 1) * ld1


@inline
def _idx3(i: int, k: int, m: int, ld1: int, ld2: int) -> int:
    return (i - 1) + (k - 1) * ld1 + (m - 1) * ld1 * ld2


@inline
def _idx3_k0(i: int, k: int, m: int, ld1: int, nk: int) -> int:
    return (i - 1) + k * ld1 + (m - 1) * ld1 * nk


@inline
def _idx4(i1: int, i2: int, i3: int, i4: int, ld1: int, ld2: int, ld3: int) -> int:
    return (
        (i1 - 1)
        + (i2 - 1) * ld1
        + (i3 - 1) * ld1 * ld2
        + (i4 - 1) * ld1 * ld2 * ld3
    )


@inline
def _idx5(i1: int, i2: int, i3: int, i4: int, i5: int, ld1: int, ld2: int, ld3: int, ld4: int) -> int:
    return (
        (i1 - 1)
        + (i2 - 1) * ld1
        + (i3 - 1) * ld1 * ld2
        + (i4 - 1) * ld1 * ld2 * ld3
        + (i5 - 1) * ld1 * ld2 * ld3 * ld4
    )


@inline
def _flux_idx(i: int, m: int, pcols: int) -> int:
    return (i - 1) + (m - 1) * pcols


@inline
def _rebin_core(
    nsrc: int,
    ntrg: int,
    src_x: Ptr[float],
    trg_x: Ptr[float],
    src: Ptr[float],
    trg: Ptr[float],
):
    for i in range(1, ntrg + 1):
        tl = trg_x[i - 1]
        if tl < src_x[nsrc]:
            sil = nsrc + 2
            for idx in range(1, nsrc + 2):
                if tl <= src_x[idx - 1]:
                    sil = idx
                    break

            tu = trg_x[i]
            siu = nsrc + 2
            for idx in range(1, nsrc + 2):
                if tu <= src_x[idx - 1]:
                    siu = idx
                    break

            y = 0.0
            sil = max(sil, 2)
            siu = min(siu, nsrc + 1)
            for si in range(sil, siu + 1):
                si1 = si - 1
                sl = max(tl, src_x[si1 - 1])
                su = min(tu, src_x[si - 1])
                y = y + (su - sl) * src[si1 - 1]

            trg[i - 1] = y / (trg_x[i] - trg_x[i - 1])
        else:
            trg[i - 1] = 0.0


@export
def rebin_codon(
    nsrc: int,
    ntrg: int,
    src_x_p: cobj,
    trg_x_p: cobj,
    src_p: cobj,
    trg_p: cobj,
):
    src_x = Ptr[float](src_x_p)
    trg_x = Ptr[float](trg_x_p)
    src = Ptr[float](src_p)
    trg = Ptr[float](trg_p)

    _rebin_core(nsrc, ntrg, src_x, trg_x, src, trg)


@export
def jlong_timestep_init_codon(
    jlong_used_flag: int,
    nsrc: int,
    ntrg: int,
    src_x_p: cobj,
    trg_x_p: cobj,
    src_p: cobj,
    trg_p: cobj,
):
    if jlong_used_flag == 0:
        return

    src_x = Ptr[float](src_x_p)
    trg_x = Ptr[float](trg_x_p)
    src = Ptr[float](src_p)
    trg = Ptr[float](trg_p)

    _rebin_core(nsrc, ntrg, src_x, trg_x, src, trg)


@export
def jlong_init_set_we_codon(
    nw: int,
    wc_p: cobj,
    wlintv_p: cobj,
    we_p: cobj,
):
    wc = Ptr[float](wc_p)
    wlintv = Ptr[float](wlintv_p)
    we = Ptr[float](we_p)

    for w in range(1, nw + 1):
        we[w - 1] = wc[w - 1] - 0.5 * wlintv[w - 1]
    we[nw] = wc[nw - 1] + 0.5 * wlintv[nw - 1]


@export
def jlong_get_xsqy_numj_codon(
    phtcnt: int,
    lng_indexer_p: cobj,
    numj_p: cobj,
):
    lng_indexer = Ptr[int](lng_indexer_p)
    numj_out = Ptr[int](numj_p)

    count = 0
    for m in range(1, phtcnt + 1):
        value = lng_indexer[m - 1]
        if value > 0:
            seen = 0
            for i in range(1, m):
                if lng_indexer[i - 1] == value:
                    seen = 1
                    break
            if seen == 0:
                count += 1

    numj_out[0] = count


@export
def jlong_get_xsqy_read_order_codon(
    phtcnt: int,
    numj: int,
    lng_indexer_p: cobj,
    read_varids_p: cobj,
):
    lng_indexer = Ptr[int](lng_indexer_p)
    read_varids = Ptr[int](read_varids_p)

    ndx = 0
    for m in range(1, phtcnt + 1):
        value = lng_indexer[m - 1]
        if value > 0:
            seen = 0
            for i in range(1, m):
                if lng_indexer[i - 1] == value:
                    seen = 1
                    break
            if seen == 0:
                if ndx < numj:
                    read_varids[ndx] = value
                ndx += 1


@export
def jlong_get_xsqy_index_map_codon(
    phtcnt: int,
    lng_indexer_p: cobj,
    wrk_ndx_p: cobj,
):
    lng_indexer = Ptr[int](lng_indexer_p)
    wrk_ndx = Ptr[int](wrk_ndx_p)

    ndx = 0
    for m in range(1, phtcnt + 1):
        if wrk_ndx[m - 1] > 0:
            ndx += 1
            value = wrk_ndx[m - 1]
            for i in range(1, phtcnt + 1):
                if wrk_ndx[i - 1] == value:
                    lng_indexer[i - 1] = ndx
                    wrk_ndx[i - 1] = -100000


@export
def jlong_get_xsqy_dprs_codon(
    np_xs: int,
    prs_p: cobj,
    dprs_p: cobj,
):
    prs = Ptr[float](prs_p)
    dprs = Ptr[float](dprs_p)

    for i in range(1, np_xs):
        dprs[i - 1] = 1.0 / (prs[i - 1] - prs[i])


@export
def jlong_get_rsf_scale_codon(
    nw: int,
    nump: int,
    numsza: int,
    numcolo3: int,
    numalb: int,
    wlintv_p: cobj,
    rsf_tab_p: cobj,
):
    wlintv = Ptr[float](wlintv_p)
    rsf_tab = Ptr[float32](rsf_tab_p)

    for w in range(1, nw + 1):
        wrk = wlintv[w - 1]
        for ial in range(1, numalb + 1):
            for iv in range(1, numcolo3 + 1):
                for is_idx in range(1, numsza + 1):
                    for iz in range(1, nump + 1):
                        idx = _idx5(w, iz, is_idx, iv, ial, nw, nump, numsza, numcolo3)
                        rsf_tab[idx] = float32(wrk * float(rsf_tab[idx]))


@export
def jlong_get_rsf_deltas_codon(
    nump: int,
    numsza: int,
    numalb: int,
    numcolo3: int,
    p_p: cobj,
    sza_p: cobj,
    alb_p: cobj,
    o3rat_p: cobj,
    del_p_p: cobj,
    del_sza_p: cobj,
    del_alb_p: cobj,
    del_o3rat_p: cobj,
):
    p = Ptr[float](p_p)
    sza = Ptr[float](sza_p)
    alb = Ptr[float](alb_p)
    o3rat = Ptr[float](o3rat_p)
    del_p = Ptr[float](del_p_p)
    del_sza = Ptr[float](del_sza_p)
    del_alb = Ptr[float](del_alb_p)
    del_o3rat = Ptr[float](del_o3rat_p)

    for i in range(1, nump):
        del_p[i - 1] = 1.0 / abs(p[i - 1] - p[i])
    for i in range(1, numsza):
        del_sza[i - 1] = 1.0 / (sza[i] - sza[i - 1])
    for i in range(1, numalb):
        del_alb[i - 1] = 1.0 / (alb[i] - alb[i - 1])
    for i in range(1, numcolo3):
        del_o3rat[i - 1] = 1.0 / (o3rat[i] - o3rat[i - 1])


@export
def jlong_get_rsf_bde_codon(
    nw: int,
    use_bde_flag: int,
    hc_val: float,
    wc_o2_b_val: float,
    wc_o3_a_val: float,
    wc_o3_b_val: float,
    wc_p: cobj,
    bde_o2_b_p: cobj,
    bde_o3_a_p: cobj,
    bde_o3_b_p: cobj,
):
    wc = Ptr[float](wc_p)
    bde_o2_b = Ptr[float](bde_o2_b_p)
    bde_o3_a = Ptr[float](bde_o3_a_p)
    bde_o3_b = Ptr[float](bde_o3_b_p)

    if use_bde_flag != 0:
        for i in range(1, nw + 1):
            wc_i = wc[i - 1]
            bde_o2_b[i - 1] = max(0.0, hc_val * (wc_o2_b_val - wc_i) / (wc_o2_b_val * wc_i))
            bde_o3_a[i - 1] = max(0.0, hc_val * (wc_o3_a_val - wc_i) / (wc_o3_a_val * wc_i))
            bde_o3_b[i - 1] = max(0.0, hc_val * (wc_o3_b_val - wc_i) / (wc_o3_b_val * wc_i))
    else:
        for i in range(1, nw + 1):
            wc_i = wc[i - 1]
            value = hc_val / wc_i
            bde_o2_b[i - 1] = value
            bde_o3_a[i - 1] = value
            bde_o3_b[i - 1] = value


@export
def zenith_codon(
    ncol: int,
    calday: float,
    pi_val: float,
    delta: float,
    clat_p: cobj,
    clon_p: cobj,
    coszrs_p: cobj,
):
    clat = Ptr[float](clat_p)
    clon = Ptr[float](clon_p)
    coszrs = Ptr[float](coszrs_p)

    for i in range(1, ncol + 1):
        coszrs[i - 1] = sin(clat[i - 1]) * sin(delta) - cos(clat[i - 1]) * cos(delta) * cos(
            calday * 2.0 * pi_val + clon[i - 1]
        )


@export
def jlong_interpolate_rsf_codon(
    nw: int,
    nump: int,
    numsza: int,
    numalb: int,
    numcolo3: int,
    kbot: int,
    sza_in: float,
    alb_in_p: cobj,
    p_in_p: cobj,
    colo3_in_p: cobj,
    p_grid_p: cobj,
    del_p_p: cobj,
    sza_grid_p: cobj,
    del_sza_p: cobj,
    alb_grid_p: cobj,
    del_alb_p: cobj,
    o3rat_p: cobj,
    del_o3rat_p: cobj,
    colo3_grid_p: cobj,
    rsf_tab_p: cobj,
    etfphot_p: cobj,
    psum_l_p: cobj,
    rsf_p: cobj,
):
    alb_in = Ptr[float](alb_in_p)
    p_in = Ptr[float](p_in_p)
    colo3_in = Ptr[float](colo3_in_p)
    p_grid = Ptr[float](p_grid_p)
    del_p = Ptr[float](del_p_p)
    sza_grid = Ptr[float](sza_grid_p)
    del_sza = Ptr[float](del_sza_p)
    alb_grid = Ptr[float](alb_grid_p)
    del_alb = Ptr[float](del_alb_p)
    o3rat = Ptr[float](o3rat_p)
    del_o3rat = Ptr[float](del_o3rat_p)
    colo3_grid = Ptr[float](colo3_grid_p)
    rsf_tab = Ptr[float32](rsf_tab_p)
    etfphot = Ptr[float](etfphot_p)
    psum_l = Ptr[float](psum_l_p)
    rsf = Ptr[float](rsf_p)

    is_idx = numsza + 1
    for idx in range(1, numsza + 1):
        if sza_grid[idx - 1] > sza_in:
            is_idx = idx
            break
    is_idx = max(min(is_idx, numsza) - 1, 1)
    isp1 = is_idx + 1
    dels1 = max(0.0, min(1.0, (sza_in - sza_grid[is_idx - 1]) * del_sza[is_idx - 1]))
    wrk0 = 1.0 - dels1

    izl = 2
    for k in range(kbot, 0, -1):
        ial = numalb + 1
        for idx in range(1, numalb + 1):
            if alb_grid[idx - 1] > alb_in[k - 1]:
                ial = idx
                break
        albind = max(min(ial, numalb) - 1, 1)

        if p_in[k - 1] > p_grid[0]:
            pind = 2
            wght1 = 1.0
        elif p_in[k - 1] <= p_grid[nump - 1]:
            pind = nump
            wght1 = 0.0
        else:
            iz = nump + 1
            for idx in range(izl, nump + 1):
                if p_grid[idx - 1] < p_in[k - 1]:
                    iz = idx
                    izl = idx
                    break
            pind = max(min(iz, nump), 2)
            wght1 = max(0.0, min(1.0, (p_in[k - 1] - p_grid[pind - 1]) * del_p[pind - 2]))

        v3ratu = colo3_in[k - 1] / colo3_grid[pind - 2]
        iv = numcolo3 + 1
        for idx in range(1, numcolo3 + 1):
            if o3rat[idx - 1] > v3ratu:
                iv = idx
                break
        ratindu = max(min(iv, numcolo3) - 1, 1)

        if colo3_grid[pind - 1] != 0.0:
            v3ratl = colo3_in[k - 1] / colo3_grid[pind - 1]
            iv = numcolo3 + 1
            for idx in range(1, numcolo3 + 1):
                if o3rat[idx - 1] > v3ratl:
                    iv = idx
                    break
            ratindl = max(min(iv, numcolo3) - 1, 1)
        else:
            ratindl = ratindu
            v3ratl = o3rat[ratindu - 1]

        ial = albind
        ialp1 = ial + 1
        iv = ratindl

        dels2 = max(0.0, min(1.0, (v3ratl - o3rat[iv - 1]) * del_o3rat[iv - 1]))
        dels3 = max(0.0, min(1.0, (alb_in[k - 1] - alb_grid[ial - 1]) * del_alb[ial - 1]))

        wrk1 = (1.0 - dels2) * (1.0 - dels3)
        wghtl000 = wrk0 * wrk1
        wghtl100 = dels1 * wrk1
        wrk1 = (1.0 - dels2) * dels3
        wghtl001 = wrk0 * wrk1
        wghtl101 = dels1 * wrk1
        wrk1 = dels2 * (1.0 - dels3)
        wghtl010 = wrk0 * wrk1
        wghtl110 = dels1 * wrk1
        wrk1 = dels2 * dels3
        wghtl011 = wrk0 * wrk1
        wghtl111 = dels1 * wrk1

        iv = ratindu
        dels2 = max(0.0, min(1.0, (v3ratu - o3rat[iv - 1]) * del_o3rat[iv - 1]))

        wrk1 = (1.0 - dels2) * (1.0 - dels3)
        wghtu000 = wrk0 * wrk1
        wghtu100 = dels1 * wrk1
        wrk1 = (1.0 - dels2) * dels3
        wghtu001 = wrk0 * wrk1
        wghtu101 = dels1 * wrk1
        wrk1 = dels2 * (1.0 - dels3)
        wghtu010 = wrk0 * wrk1
        wghtu110 = dels1 * wrk1
        wrk1 = dels2 * dels3
        wghtu011 = wrk0 * wrk1
        wghtu111 = dels1 * wrk1

        iz = pind
        iv = ratindl
        ivp1 = iv + 1
        for wn in range(1, nw + 1):
            psum_l[wn - 1] = (
                wghtl000 * float(rsf_tab[_idx5(wn, iz, is_idx, iv, ial, nw, nump, numsza, numcolo3)])
                + wghtl001 * float(rsf_tab[_idx5(wn, iz, is_idx, iv, ialp1, nw, nump, numsza, numcolo3)])
                + wghtl010 * float(rsf_tab[_idx5(wn, iz, is_idx, ivp1, ial, nw, nump, numsza, numcolo3)])
                + wghtl011 * float(rsf_tab[_idx5(wn, iz, is_idx, ivp1, ialp1, nw, nump, numsza, numcolo3)])
                + wghtl100 * float(rsf_tab[_idx5(wn, iz, isp1, iv, ial, nw, nump, numsza, numcolo3)])
                + wghtl101 * float(rsf_tab[_idx5(wn, iz, isp1, iv, ialp1, nw, nump, numsza, numcolo3)])
                + wghtl110 * float(rsf_tab[_idx5(wn, iz, isp1, ivp1, ial, nw, nump, numsza, numcolo3)])
                + wghtl111 * float(rsf_tab[_idx5(wn, iz, isp1, ivp1, ialp1, nw, nump, numsza, numcolo3)])
            )

        iz = iz - 1
        iv = ratindu
        ivp1 = iv + 1
        for wn in range(1, nw + 1):
            psum_u = (
                wghtu000 * float(rsf_tab[_idx5(wn, iz, is_idx, iv, ial, nw, nump, numsza, numcolo3)])
                + wghtu001 * float(rsf_tab[_idx5(wn, iz, is_idx, iv, ialp1, nw, nump, numsza, numcolo3)])
                + wghtu010 * float(rsf_tab[_idx5(wn, iz, is_idx, ivp1, ial, nw, nump, numsza, numcolo3)])
                + wghtu011 * float(rsf_tab[_idx5(wn, iz, is_idx, ivp1, ialp1, nw, nump, numsza, numcolo3)])
                + wghtu100 * float(rsf_tab[_idx5(wn, iz, isp1, iv, ial, nw, nump, numsza, numcolo3)])
                + wghtu101 * float(rsf_tab[_idx5(wn, iz, isp1, iv, ialp1, nw, nump, numsza, numcolo3)])
                + wghtu110 * float(rsf_tab[_idx5(wn, iz, isp1, ivp1, ial, nw, nump, numsza, numcolo3)])
                + wghtu111 * float(rsf_tab[_idx5(wn, iz, isp1, ivp1, ialp1, nw, nump, numsza, numcolo3)])
            )
            rsf[_idx2(wn, k, nw)] = psum_l[wn - 1] + wght1 * (psum_u - psum_l[wn - 1])

        for wn in range(1, nw + 1):
            rsf[_idx2(wn, k, nw)] = etfphot[wn - 1] * rsf[_idx2(wn, k, nw)]


@export
def jlong_photo_fill_xswk_codon(
    numj: int,
    nw: int,
    nt: int,
    np_xs: int,
    k: int,
    p_in_p: cobj,
    t_in_p: cobj,
    prs_p: cobj,
    dprs_p: cobj,
    xsqy_p: cobj,
    xswk_p: cobj,
):
    p_in = Ptr[float](p_in_p)
    t_in = Ptr[float](t_in_p)
    prs = Ptr[float](prs_p)
    dprs = Ptr[float](dprs_p)
    xsqy = Ptr[float32](xsqy_p)
    xswk = Ptr[float](xswk_p)

    t_index = int(t_in[k - 1] - 148.5)
    t_index = min(201, max(t_index, 1))
    ptarget = p_in[k - 1]

    if ptarget >= prs[0]:
        pndx = 1
        for wn in range(1, nw + 1):
            for m in range(1, numj + 1):
                xswk[_idx2(m, wn, numj)] = float(xsqy[_idx4(m, wn, t_index, pndx, numj, nw, nt)])
    elif ptarget <= prs[np_xs - 1]:
        pndx = np_xs
        for wn in range(1, nw + 1):
            for m in range(1, numj + 1):
                xswk[_idx2(m, wn, numj)] = float(xsqy[_idx4(m, wn, t_index, pndx, numj, nw, nt)])
    else:
        pndx = np_xs - 1
        delp = 0.0
        for km in range(2, np_xs + 1):
            if ptarget >= prs[km - 1]:
                pndx = km - 1
                delp = (prs[pndx - 1] - ptarget) * dprs[pndx - 1]
                break
        for wn in range(1, nw + 1):
            for m in range(1, numj + 1):
                lo = float(xsqy[_idx4(m, wn, t_index, pndx, numj, nw, nt)])
                hi = float(xsqy[_idx4(m, wn, t_index, pndx + 1, numj, nw, nt)])
                xswk[_idx2(m, wn, numj)] = lo + delp * (hi - lo)


@export
def jlong_photo_accum_codon(
    numj: int,
    nw: int,
    xswk_p: cobj,
    rsf_col_p: cobj,
    j_long_col_p: cobj,
):
    xswk = Ptr[float](xswk_p)
    rsf_col = Ptr[float](rsf_col_p)
    j_long_col = Ptr[float](j_long_col_p)

    for m in range(1, numj + 1):
        acc = 0.0
        for wn in range(1, nw + 1):
            acc = acc + xswk[_idx2(m, wn, numj)] * rsf_col[wn - 1]
        j_long_col[m - 1] = acc


@export
def jlong_photo_loop_codon(
    numj: int,
    nw: int,
    nt: int,
    np_xs: int,
    nlev: int,
    p_in_p: cobj,
    t_in_p: cobj,
    prs_p: cobj,
    dprs_p: cobj,
    xsqy_p: cobj,
    rsf_p: cobj,
    xswk_p: cobj,
    j_long_p: cobj,
):
    p_in = Ptr[float](p_in_p)
    t_in = Ptr[float](t_in_p)
    prs = Ptr[float](prs_p)
    dprs = Ptr[float](dprs_p)
    xsqy = Ptr[float32](xsqy_p)
    rsf = Ptr[float](rsf_p)
    xswk = Ptr[float](xswk_p)
    j_long = Ptr[float](j_long_p)

    for k in range(1, nlev + 1):
        t_index = int(t_in[k - 1] - 148.5)
        t_index = min(201, max(t_index, 1))
        ptarget = p_in[k - 1]

        if ptarget >= prs[0]:
            pndx = 1
            for wn in range(1, nw + 1):
                for m in range(1, numj + 1):
                    xswk[_idx2(m, wn, numj)] = float(xsqy[_idx4(m, wn, t_index, pndx, numj, nw, nt)])
        elif ptarget <= prs[np_xs - 1]:
            pndx = np_xs
            for wn in range(1, nw + 1):
                for m in range(1, numj + 1):
                    xswk[_idx2(m, wn, numj)] = float(xsqy[_idx4(m, wn, t_index, pndx, numj, nw, nt)])
        else:
            pndx = np_xs - 1
            delp = 0.0
            for km in range(2, np_xs + 1):
                if ptarget >= prs[km - 1]:
                    pndx = km - 1
                    delp = (prs[pndx - 1] - ptarget) * dprs[pndx - 1]
                    break
            for wn in range(1, nw + 1):
                for m in range(1, numj + 1):
                    lo = float(xsqy[_idx4(m, wn, t_index, pndx, numj, nw, nt)])
                    hi = float(xsqy[_idx4(m, wn, t_index, pndx + 1, numj, nw, nt)])
                    xswk[_idx2(m, wn, numj)] = lo + delp * (hi - lo)

        for m in range(1, numj + 1):
            acc = 0.0
            for wn in range(1, nw + 1):
                acc = acc + xswk[_idx2(m, wn, numj)] * rsf[_idx2(wn, k, nw)]
            j_long[_idx2(m, k, numj)] = acc


@export
def jlong_photo_codon(
    numj: int,
    nw: int,
    nt: int,
    np_xs: int,
    nump: int,
    numsza: int,
    numalb: int,
    numcolo3: int,
    nlev: int,
    sza_in: float,
    alb_in_p: cobj,
    p_in_p: cobj,
    t_in_p: cobj,
    colo3_in_p: cobj,
    p_grid_p: cobj,
    del_p_p: cobj,
    sza_grid_p: cobj,
    del_sza_p: cobj,
    alb_grid_p: cobj,
    del_alb_p: cobj,
    o3rat_p: cobj,
    del_o3rat_p: cobj,
    colo3_grid_p: cobj,
    rsf_tab_p: cobj,
    etfphot_p: cobj,
    prs_p: cobj,
    dprs_p: cobj,
    xsqy_p: cobj,
    rsf_p: cobj,
    psum_l_p: cobj,
    xswk_p: cobj,
    j_long_p: cobj,
):
    jlong_interpolate_rsf_codon(
        nw,
        nump,
        numsza,
        numalb,
        numcolo3,
        nlev,
        sza_in,
        alb_in_p,
        p_in_p,
        colo3_in_p,
        p_grid_p,
        del_p_p,
        sza_grid_p,
        del_sza_p,
        alb_grid_p,
        del_alb_p,
        o3rat_p,
        del_o3rat_p,
        colo3_grid_p,
        rsf_tab_p,
        etfphot_p,
        psum_l_p,
        rsf_p,
    )
    jlong_photo_loop_codon(
        numj,
        nw,
        nt,
        np_xs,
        nlev,
        p_in_p,
        t_in_p,
        prs_p,
        dprs_p,
        xsqy_p,
        rsf_p,
        xswk_p,
        j_long_p,
    )


@export
def chem_emissions_zero_cflx_codon(
    pcols: int,
    pcnst: int,
    map2chm_p: cobj,
    cflx_p: cobj,
):
    map2chm = Ptr[int](map2chm_p)
    cflx = Ptr[float](cflx_p)

    for m in range(2, pcnst + 1):
        if map2chm[m - 1] > 0:
            for i in range(1, pcols + 1):
                cflx[_flux_idx(i, m, pcols)] = 0.0


@export
def gas_phase_chemdr_prepare_sza_codon(
    ncol: int,
    rad2deg: float,
    zen_angle_p: cobj,
    sza_p: cobj,
):
    zen_angle = Ptr[float](zen_angle_p)
    sza = Ptr[float](sza_p)

    for i in range(1, ncol + 1):
        z = acos(zen_angle[i - 1])
        zen_angle[i - 1] = z
        sza[i - 1] = z * rad2deg


@export
def table_photo_zero_photos_codon(
    ncol: int,
    pver: int,
    phtcnt: int,
    photos_p: cobj,
):
    photos = Ptr[float](photos_p)

    for m in range(1, phtcnt + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                photos[_idx3(i, k, m, ncol, pver)] = 0.0


@export
def table_photo_daylight_setup_codon(
    ncol: int,
    pcols: int,
    pver: int,
    ncol_abs: int,
    i_col: int,
    p1: int,
    p2: int,
    pa2mb: float,
    pmid_p: cobj,
    pdel_p: cobj,
    col_dens_p: cobj,
    lwc_p: cobj,
    clouds_p: cobj,
    temper_p: cobj,
    zmid_p: cobj,
    parg_p: cobj,
    colo3_p: cobj,
    fac1_p: cobj,
    lwc_line_p: cobj,
    cld_line_p: cobj,
    tline_p: cobj,
    zarg_p: cobj,
):
    pmid = Ptr[float](pmid_p)
    pdel = Ptr[float](pdel_p)
    col_dens = Ptr[float](col_dens_p)
    lwc = Ptr[float](lwc_p)
    clouds = Ptr[float](clouds_p)
    temper = Ptr[float](temper_p)
    zmid = Ptr[float](zmid_p)
    parg = Ptr[float](parg_p)
    colo3 = Ptr[float](colo3_p)
    fac1 = Ptr[float](fac1_p)
    lwc_line = Ptr[float](lwc_line_p)
    cld_line = Ptr[float](cld_line_p)
    tline = Ptr[float](tline_p)
    zarg = Ptr[float](zarg_p)

    for k in range(1, pver + 1):
        parg[k - 1] = pa2mb * pmid[_idx2(i_col, k, pcols)]
        colo3[k - 1] = col_dens[_idx3(i_col, k, 1, ncol, pver)]
        fac1[k - 1] = pdel[_idx2(i_col, k, pcols)]
        lwc_line[k - 1] = lwc[_idx2(i_col, k, ncol)]
        cld_line[k - 1] = clouds[_idx2(i_col, k, ncol)]

    for k in range(p1, p2 + 1):
        src_k = k - p1 + 1
        tline[k - 1] = temper[_idx2(i_col, src_k, pcols)]
        zarg[k - 1] = zmid[_idx2(i_col, src_k, ncol)]


@export
def table_photo_scale_cld_mult_codon(
    pver: int,
    esfact: float,
    cld_mult_p: cobj,
):
    cld_mult = Ptr[float](cld_mult_p)

    for k in range(1, pver + 1):
        cld_mult[k - 1] = esfact * cld_mult[k - 1]


@export
def mmr2vmr_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    mbar_p: cobj,
    mmr_p: cobj,
    adv_mass_p: cobj,
    vmr_p: cobj,
):
    mbar = Ptr[float](mbar_p)
    mmr = Ptr[float](mmr_p)
    adv_mass = Ptr[float](adv_mass_p)
    vmr = Ptr[float](vmr_p)

    for m in range(1, gas_pcnst + 1):
        adv = adv_mass[m - 1]
        if adv != 0.0:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    vmr[_idx3(i, k, m, ncol, pver)] = mbar[_idx2(i, k, ncol)] * mmr[
                        _idx3(i, k, m, pcols, pver)
                    ] / adv


@export
def vmr2mmr_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    mbar_p: cobj,
    vmr_p: cobj,
    adv_mass_p: cobj,
    mmr_p: cobj,
):
    mbar = Ptr[float](mbar_p)
    vmr = Ptr[float](vmr_p)
    adv_mass = Ptr[float](adv_mass_p)
    mmr = Ptr[float](mmr_p)

    for m in range(1, gas_pcnst + 1):
        adv = adv_mass[m - 1]
        if adv != 0.0:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    mmr[_idx3(i, k, m, pcols, pver)] = adv * vmr[_idx3(i, k, m, ncol, pver)] / mbar[
                        _idx2(i, k, ncol)
                    ]


@export
def h2o_to_vmr_codon(
    ncol: int,
    pcols: int,
    pver: int,
    adv_mass_h2o: float,
    h2o_mmr_p: cobj,
    mbar_p: cobj,
    h2o_vmr_p: cobj,
):
    h2o_mmr = Ptr[float](h2o_mmr_p)
    mbar = Ptr[float](mbar_p)
    h2o_vmr = Ptr[float](h2o_vmr_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            h2o_vmr[_idx2(i, k, ncol)] = mbar[_idx2(i, k, ncol)] * h2o_mmr[_idx2(i, k, pcols)] / adv_mass_h2o


@export
def set_mean_mass_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    id_o2: int,
    id_o: int,
    id_h: int,
    id_n: int,
    fixed_mbar: int,
    mwdry: float,
    mmr_p: cobj,
    adv_mass_p: cobj,
    mbar_p: cobj,
):
    mmr = Ptr[float](mmr_p)
    adv_mass = Ptr[float](adv_mass_p)
    mbar = Ptr[float](mbar_p)

    if fixed_mbar != 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                mbar[_idx2(i, k, ncol)] = mwdry
        return

    adv_n = adv_mass[id_n - 1]
    adv_o2 = adv_mass[id_o2 - 1]
    adv_o = adv_mass[id_o - 1]
    adv_h = adv_mass[id_h - 1]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            xn2 = 1.0 - (
                mmr[_idx3(i, k, id_o2, pcols, pver)]
                + mmr[_idx3(i, k, id_o, pcols, pver)]
                + mmr[_idx3(i, k, id_h, pcols, pver)]
            )
            fn2 = 0.5 * xn2 / adv_n
            fo2 = mmr[_idx3(i, k, id_o2, pcols, pver)] / adv_o2
            fo = mmr[_idx3(i, k, id_o, pcols, pver)] / adv_o
            fh = mmr[_idx3(i, k, id_h, pcols, pver)] / adv_h
            mbar[_idx2(i, k, ncol)] = 1.0 / (fn2 + fo2 + fo + fh)


@export
def setinv_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    nfs: int,
    m_ndx: int,
    n2_ndx: int,
    o2_ndx: int,
    h2o_ndx: int,
    id_o: int,
    id_o2: int,
    id_h: int,
    has_n2: int,
    has_o2: int,
    has_h2o: int,
    has_var_o2: int,
    pa_xfac: float,
    boltz_cgs: float,
    tfld_p: cobj,
    h2ovmr_p: cobj,
    vmr_p: cobj,
    pmid_p: cobj,
    invariants_p: cobj,
):
    tfld = Ptr[float](tfld_p)
    h2ovmr = Ptr[float](h2ovmr_p)
    vmr = Ptr[float](vmr_p)
    pmid = Ptr[float](pmid_p)
    invariants = Ptr[float](invariants_p)

    for m in range(1, nfs + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                invariants[_idx3(i, k, m, ncol, pver)] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            invariants[_idx3(i, k, m_ndx, ncol, pver)] = (
                pa_xfac * pmid[_idx2(i, k, pcols)] / (boltz_cgs * tfld[_idx2(i, k, pcols)])
            )

    if has_n2 != 0:
        if has_var_o2 != 0:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    sum1 = (
                        vmr[_idx3(i, k, id_o, ncol, pver)]
                        + vmr[_idx3(i, k, id_o2, ncol, pver)]
                        + vmr[_idx3(i, k, id_h, ncol, pver)]
                    )
                    invariants[_idx3(i, k, n2_ndx, ncol, pver)] = (
                        (1.0 - sum1) * invariants[_idx3(i, k, m_ndx, ncol, pver)]
                    )
        else:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    invariants[_idx3(i, k, n2_ndx, ncol, pver)] = (
                        0.79 * invariants[_idx3(i, k, m_ndx, ncol, pver)]
                    )

    if has_o2 != 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                invariants[_idx3(i, k, o2_ndx, ncol, pver)] = (
                    0.21 * invariants[_idx3(i, k, m_ndx, ncol, pver)]
                )

    if has_h2o != 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                invariants[_idx3(i, k, h2o_ndx, ncol, pver)] = (
                    h2ovmr[_idx2(i, k, ncol)] * invariants[_idx3(i, k, m_ndx, ncol, pver)]
                )


@export
def charge_balance_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    np_ndx: int,
    n2p_ndx: int,
    op_ndx: int,
    o2p_ndx: int,
    nop_ndx: int,
    conc_p: cobj,
    wrk_p: cobj,
):
    conc = Ptr[float](conc_p)
    wrk = Ptr[float](wrk_p)

    if np_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                wrk[_idx2(i, k, ncol)] = wrk[_idx2(i, k, ncol)] + conc[_idx3(i, k, np_ndx, ncol, pver)]

    if n2p_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                wrk[_idx2(i, k, ncol)] = wrk[_idx2(i, k, ncol)] + conc[_idx3(i, k, n2p_ndx, ncol, pver)]

    if op_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                wrk[_idx2(i, k, ncol)] = wrk[_idx2(i, k, ncol)] + conc[_idx3(i, k, op_ndx, ncol, pver)]

    if o2p_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                wrk[_idx2(i, k, ncol)] = wrk[_idx2(i, k, ncol)] + conc[_idx3(i, k, o2p_ndx, ncol, pver)]

    if nop_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                wrk[_idx2(i, k, ncol)] = wrk[_idx2(i, k, ncol)] + conc[_idx3(i, k, nop_ndx, ncol, pver)]


@export
def setcol_codon(
    ncol: int,
    pver: int,
    ncol_abs: int,
    col_delta_p: cobj,
    col_dens_p: cobj,
):
    col_delta = Ptr[float](col_delta_p)
    col_dens = Ptr[float](col_dens_p)

    for m in range(1, ncol_abs + 1):
        for i in range(1, ncol + 1):
            col_dens[_idx3(i, 1, m, ncol, pver)] = col_delta[_idx3_k0(i, 0, m, ncol, pver + 1)] + 0.5 * col_delta[
                _idx3_k0(i, 1, m, ncol, pver + 1)
            ]

        for k in range(2, pver + 1):
            km1 = k - 1
            for i in range(1, ncol + 1):
                col_dens[_idx3(i, k, m, ncol, pver)] = col_dens[_idx3(i, km1, m, ncol, pver)] + 0.5 * (
                    col_delta[_idx3_k0(i, km1, m, ncol, pver + 1)] + col_delta[_idx3_k0(i, k, m, ncol, pver + 1)]
                )


@export
def set_ub_col_codon(
    ncol: int,
    pcols: int,
    pver: int,
    ncol_abs: int,
    indexm: int,
    o3rad_ndx: int,
    ox_ndx: int,
    o3_ndx: int,
    o3_inv_ndx: int,
    o2_ndx: int,
    o2_is_inv: int,
    xfactor: float,
    pdel_p: cobj,
    vmr_p: cobj,
    invariants_p: cobj,
    o2_exo_col_p: cobj,
    o3_exo_col_p: cobj,
    col_delta_p: cobj,
):
    pdel = Ptr[float](pdel_p)
    vmr = Ptr[float](vmr_p)
    invariants = Ptr[float](invariants_p)
    o2_exo_col = Ptr[float](o2_exo_col_p)
    o3_exo_col = Ptr[float](o3_exo_col_p)
    col_delta = Ptr[float](col_delta_p)

    spc_ndx = o3rad_ndx
    if spc_ndx <= 0:
        spc_ndx = ox_ndx
    if spc_ndx < 1:
        spc_ndx = o3_ndx

    if spc_ndx > 0:
        for i in range(1, ncol + 1):
            col_delta[_idx3_k0(i, 0, 1, ncol, pver + 1)] = o3_exo_col[i - 1]
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                col_delta[_idx3_k0(i, k, 1, ncol, pver + 1)] = (
                    xfactor * pdel[_idx2(i, k, pcols)] * vmr[_idx3(i, k, spc_ndx, ncol, pver)]
                )
    elif o3_inv_ndx > 0:
        for i in range(1, ncol + 1):
            col_delta[_idx3_k0(i, 0, 1, ncol, pver + 1)] = o3_exo_col[i - 1]
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                col_delta[_idx3_k0(i, k, 1, ncol, pver + 1)] = (
                    xfactor
                    * pdel[_idx2(i, k, pcols)]
                    * invariants[_idx3(i, k, o3_inv_ndx, ncol, pver)]
                    / invariants[_idx3(i, k, indexm, ncol, pver)]
                )
    else:
        for k in range(0, pver + 1):
            for i in range(1, ncol + 1):
                col_delta[_idx3_k0(i, k, 1, ncol, pver + 1)] = 0.0

    if ncol_abs > 1:
        if o2_ndx > 1:
            for i in range(1, ncol + 1):
                col_delta[_idx3_k0(i, 0, 2, ncol, pver + 1)] = o2_exo_col[i - 1]
            if o2_is_inv != 0:
                for k in range(1, pver + 1):
                    for i in range(1, ncol + 1):
                        col_delta[_idx3_k0(i, k, 2, ncol, pver + 1)] = (
                            xfactor
                            * pdel[_idx2(i, k, pcols)]
                            * invariants[_idx3(i, k, o2_ndx, ncol, pver)]
                            / invariants[_idx3(i, k, indexm, ncol, pver)]
                        )
            else:
                for k in range(1, pver + 1):
                    for i in range(1, ncol + 1):
                        col_delta[_idx3_k0(i, k, 2, ncol, pver + 1)] = (
                            xfactor * pdel[_idx2(i, k, pcols)] * vmr[_idx3(i, k, o2_ndx, ncol, pver)]
                        )
        else:
            for k in range(0, pver + 1):
                for i in range(1, ncol + 1):
                    col_delta[_idx3_k0(i, k, 2, ncol, pver + 1)] = 0.0


@export
def cloud_mod_codon(
    pver: int,
    zen_angle: float,
    srf_alb: float,
    rgrav: float,
    clouds_p: cobj,
    lwc_p: cobj,
    delp_p: cobj,
    eff_alb_p: cobj,
    cld_mult_p: cobj,
    del_lwp_p: cobj,
    del_tau_p: cobj,
    above_tau_p: cobj,
    below_tau_p: cobj,
    above_cld_p: cobj,
    below_cld_p: cobj,
    above_tra_p: cobj,
    below_tra_p: cobj,
    fac1_p: cobj,
    fac2_p: cobj,
):
    clouds = Ptr[float](clouds_p)
    lwc = Ptr[float](lwc_p)
    delp = Ptr[float](delp_p)
    eff_alb = Ptr[float](eff_alb_p)
    cld_mult = Ptr[float](cld_mult_p)
    del_lwp = Ptr[float](del_lwp_p)
    del_tau = Ptr[float](del_tau_p)
    above_tau = Ptr[float](above_tau_p)
    below_tau = Ptr[float](below_tau_p)
    above_cld = Ptr[float](above_cld_p)
    below_cld = Ptr[float](below_cld_p)
    above_tra = Ptr[float](above_tra_p)
    below_tra = Ptr[float](below_tra_p)
    fac1 = Ptr[float](fac1_p)
    fac2 = Ptr[float](fac2_p)

    for k in range(1, pver + 1):
        if clouds[k - 1] != 0.0:
            del_lwp[k - 1] = rgrav * lwc[k - 1] * delp[k - 1] * 1.0e3 / clouds[k - 1]
        else:
            del_lwp[k - 1] = 0.0

    for k in range(1, pver + 1):
        if clouds[k - 1] != 0.0:
            del_tau[k - 1] = del_lwp[k - 1] * 0.155 * (clouds[k - 1] ** 1.5)
        else:
            del_tau[k - 1] = 0.0

    above_tau[0] = 0.0
    for k in range(1, pver):
        above_tau[k] = del_tau[k - 1] + above_tau[k - 1]

    below_tau[pver - 1] = 0.0
    for k in range(pver - 1, 0, -1):
        below_tau[k - 1] = del_tau[k] + below_tau[k]

    above_cld[0] = 0.0
    for k in range(1, pver):
        above_cld[k] = clouds[k - 1] * del_tau[k - 1] + above_cld[k - 1]
    for k in range(1, pver):
        if above_tau[k] != 0.0:
            above_cld[k] = above_cld[k] / above_tau[k]
        else:
            above_cld[k] = above_cld[k - 1]

    below_cld[pver - 1] = 0.0
    for k in range(pver - 1, 0, -1):
        below_cld[k - 1] = clouds[k] * del_tau[k] + below_cld[k]
    for k in range(pver - 2, -1, -1):
        if below_tau[k] != 0.0:
            below_cld[k] = below_cld[k] / below_tau[k]
        else:
            below_cld[k] = below_cld[k + 1]

    for k in range(1, pver):
        if above_cld[k] != 0.0:
            above_tau[k] = above_tau[k] / above_cld[k]
    for k in range(0, pver - 1):
        if below_cld[k] != 0.0:
            below_tau[k] = below_tau[k] / below_cld[k]

    for k in range(1, pver):
        if above_tau[k] < 5.0:
            above_cld[k] = 0.0
    for k in range(0, pver - 1):
        if below_tau[k] < 5.0:
            below_cld[k] = 0.0

    for k in range(1, pver + 1):
        above_tra[k - 1] = 11.905 / (9.524 + above_tau[k - 1])
        below_tra[k - 1] = 11.905 / (9.524 + below_tau[k - 1])

    for k in range(1, pver + 1):
        if below_cld[k - 1] != 0.0:
            eff_alb[k - 1] = srf_alb + below_cld[k - 1] * (1.0 - below_tra[k - 1]) * (1.0 - srf_alb)
        else:
            eff_alb[k - 1] = srf_alb

    coschi = cos(zen_angle)
    if coschi < 0.5:
        coschi = 0.5

    for k in range(1, pver + 1):
        if del_lwp[k - 1] * 0.155 < 5.0:
            fac1[k - 1] = 0.0
        else:
            fac1[k - 1] = 1.4 * coschi - 1.0

    for k in range(1, pver + 1):
        fac2_val = 1.6 * coschi * above_tra[k - 1] - 1.0
        if fac2_val > 0.0:
            fac2[k - 1] = 0.0
        else:
            fac2[k - 1] = fac2_val

    for k in range(1, pver + 1):
        cld_mult[k - 1] = 1.0 + fac1[k - 1] * clouds[k - 1] + fac2[k - 1] * above_cld[k - 1]
        if cld_mult[k - 1] < 0.05:
            cld_mult[k - 1] = 0.05


@export
def photo_inti_fixed_press_setup_codon(
    pinterp: float,
    n_exo_levs: int,
    levs_p: cobj,
    ki_p: cobj,
    delp_p: cobj,
):
    levs = Ptr[float](levs_p)
    ki_out = Ptr[int](ki_p)
    delp_out = Ptr[float](delp_p)

    if pinterp <= levs[0]:
        ki_out[0] = 1
        delp_out[0] = 0.0
        return

    ki_val = 2
    for idx in range(2, n_exo_levs + 1):
        ki_val = idx
        if pinterp <= levs[idx - 1]:
            ki_out[0] = idx
            delp_out[0] = log(pinterp / levs[idx - 2]) / log(levs[idx - 1] / levs[idx - 2])
            return

    ki_out[0] = ki_val
    delp_out[0] = 0.0


@export
def photo_timestep_init_exo_time_codon(
    calday: float,
    days_p: cobj,
    next_p: cobj,
    last_p: cobj,
    dels_p: cobj,
):
    days = Ptr[float](days_p)
    next_v = Ptr[int](next_p)
    last_v = Ptr[int](last_p)
    dels_v = Ptr[float](dels_p)

    if calday < days[0]:
        next_v[0] = 1
        last_v[0] = 12
        dels_v[0] = (365.0 + calday - days[11]) / (365.0 + days[0] - days[11])
        return

    if calday >= days[11]:
        next_v[0] = 1
        last_v[0] = 12
        dels_v[0] = (calday - days[11]) / (365.0 + days[0] - days[11])
        return

    m = 0
    for idx in range(10, -1, -1):
        if calday >= days[idx]:
            m = idx
            break

    last_v[0] = m + 1
    next_v[0] = m + 2
    dels_v[0] = (calday - days[m]) / (days[m + 1] - days[m])


@export
def table_photo_jlong_apply_codon(
    ncol: int,
    pver: int,
    phtcnt: int,
    nlng: int,
    i_col: int,
    photos_p: cobj,
    lng_prates_p: cobj,
    cld_mult_p: cobj,
    lng_indexer_p: cobj,
    alias_mult2_p: cobj,
):
    photos = Ptr[float](photos_p)
    lng_prates = Ptr[float](lng_prates_p)
    cld_mult = Ptr[float](cld_mult_p)
    lng_indexer = Ptr[int](lng_indexer_p)
    alias_mult2 = Ptr[float](alias_mult2_p)

    for m in range(1, phtcnt + 1):
        if lng_indexer[m - 1] > 0:
            alias_factor = alias_mult2[m - 1]
            idx_lng = lng_indexer[m - 1]
            if alias_factor == 1.0:
                for k in range(1, pver + 1):
                    photos[_idx3(i_col, k, m, ncol, pver)] = (
                        photos[_idx3(i_col, k, m, ncol, pver)] + lng_prates[_idx2(idx_lng, k, nlng)]
                    ) * cld_mult[k - 1]
            else:
                for k in range(1, pver + 1):
                    photos[_idx3(i_col, k, m, ncol, pver)] = (
                        photos[_idx3(i_col, k, m, ncol, pver)] + alias_factor * lng_prates[_idx2(idx_lng, k, nlng)]
                    ) * cld_mult[k - 1]


@export
def table_photo_jno_ho2no2_codon(
    ncol: int,
    pver: int,
    phtcnt: int,
    i_col: int,
    jno_ndx: int,
    jho2no2_ndx: int,
    do_jshort: int,
    has_o2_col: int,
    has_o3_col: int,
    zen_angle: float,
    photos_p: cobj,
    col_dens_p: cobj,
    cld_mult_p: cobj,
):
    photos = Ptr[float](photos_p)
    col_dens = Ptr[float](col_dens_p)
    cld_mult = Ptr[float](cld_mult_p)

    if jno_ndx > 0 and do_jshort == 0:
        if has_o2_col != 0 and has_o3_col != 0:
            for k in range(1, pver + 1):
                fac1 = 1.0e-8 * (abs(col_dens[_idx3(i_col, k, 2, ncol, pver)] / cos(zen_angle))) ** 0.38
                fac2 = 5.0e-19 * abs(col_dens[_idx3(i_col, k, 1, ncol, pver)] / cos(zen_angle))
                photos[_idx3(i_col, k, jno_ndx, ncol, pver)] = (
                    photos[_idx3(i_col, k, jno_ndx, ncol, pver)] + 4.5e-6 * exp(-(fac1 + fac2))
                )

    if jho2no2_ndx > 0:
        for k in range(1, pver + 1):
            photos[_idx3(i_col, k, jho2no2_ndx, ncol, pver)] = (
                photos[_idx3(i_col, k, jho2no2_ndx, ncol, pver)] + 1.0e-5 * cld_mult[k - 1]
            )


@export
def gas_phase_chemdr_zero_sulfate_codon(
    ncol: int,
    pver: int,
    sulfate_p: cobj,
):
    sulfate = Ptr[float](sulfate_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            sulfate[_idx2(i, k, ncol)] = 0.0


@export
def gas_phase_chemdr_load_prognostic_sulfate_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    so4_ndx: int,
    vmr_p: cobj,
    sulfate_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    sulfate = Ptr[float](sulfate_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            sulfate[_idx2(i, k, ncol)] = vmr[_idx3(i, k, so4_ndx, ncol, pver)]


@export
def chem_emissions_megan_flux_codon(
    ncol: int,
    pcols: int,
    megan_index: int,
    megan_weight: float,
    meganflx_p: cobj,
    cflx_p: cobj,
    megflx_p: cobj,
):
    meganflx = Ptr[float](meganflx_p)
    cflx = Ptr[float](cflx_p)
    megflx = Ptr[float](megflx_p)

    for i in range(1, ncol + 1):
        flux = -meganflx[i - 1] * megan_weight
        megflx[i - 1] = flux
        cflx[_flux_idx(i, megan_index, pcols)] += flux


@export
def chem_emissions_add_sflx_codon(
    ncol: int,
    pcols: int,
    pcnst: int,
    h2o_ndx: int,
    map2chm_p: cobj,
    cflx_p: cobj,
    sflx_p: cobj,
):
    map2chm = Ptr[int](map2chm_p)
    cflx = Ptr[float](cflx_p)
    sflx = Ptr[float](sflx_p)

    for m in range(1, pcnst + 1):
        n = map2chm[m - 1]
        if n > 0 and n != h2o_ndx:
            for i in range(1, ncol + 1):
                cflx[_flux_idx(i, m, pcols)] += sflx[_flux_idx(i, n, pcols)]


@export
def aero_model_gasaerexch_column_flux_codon(
    ncol: int,
    pcols: int,
    pver: int,
    adv_mass: float,
    gravit: float,
    field_p: cobj,
    mbar_p: cobj,
    pdel_p: cobj,
    wrk_p: cobj,
):
    field = Ptr[float](field_p)
    mbar = Ptr[float](mbar_p)
    pdel = Ptr[float](pdel_p)
    wrk = Ptr[float](wrk_p)

    for i in range(1, ncol + 1):
        wrk[i - 1] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            wrk[i - 1] += field[_idx2(i, k, ncol)] * adv_mass / mbar[_idx2(i, k, pcols)] * pdel[
                _idx2(i, k, pcols)
            ] / gravit


@export
def aero_model_gasaerexch_gas_tend_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    delt: float,
    vmr0_p: cobj,
    vmr_p: cobj,
    dvmrdt_p: cobj,
):
    vmr0 = Ptr[float](vmr0_p)
    vmr = Ptr[float](vmr_p)
    dvmrdt = Ptr[float](dvmrdt_p)

    for m in range(1, gas_pcnst + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx3(i, k, m, ncol, pver)
                dvmrdt[idx] = (vmr[idx] - vmr0[idx]) / delt


@export
def aero_model_gasaerexch_aq_tend_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    delt: float,
    vmr_p: cobj,
    vmrcw_p: cobj,
    dvmrdt_p: cobj,
    dvmrcwdt_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    vmrcw = Ptr[float](vmrcw_p)
    dvmrdt = Ptr[float](dvmrdt_p)
    dvmrcwdt = Ptr[float](dvmrcwdt_p)

    for m in range(1, gas_pcnst + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx3(i, k, m, ncol, pver)
                dvmrdt[idx] = (vmr[idx] - dvmrdt[idx]) / delt

    for m in range(1, gas_pcnst + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx3(i, k, m, ncol, pver)
                dvmrcwdt[idx] = (vmrcw[idx] - dvmrcwdt[idx]) / delt


@export
def aero_model_emissions_accumulate_sflx_codon(
    ncol: int,
    pcols: int,
    nindices: int,
    indices_p: cobj,
    cflx_p: cobj,
    sflx_p: cobj,
):
    indices = Ptr[int](indices_p)
    cflx = Ptr[float](cflx_p)
    sflx = Ptr[float](sflx_p)

    for i in range(1, pcols + 1):
        sflx[i - 1] = 0.0

    for m in range(1, nindices + 1):
        idx = indices[m - 1]
        for i in range(1, ncol + 1):
            sflx[i - 1] += cflx[_flux_idx(i, idx, pcols)]


@export
def aero_model_emissions_seasalt_wind_codon(
    ncol: int,
    pcols: int,
    pver: int,
    z0: float,
    state_u_p: cobj,
    state_v_p: cobj,
    state_zm_p: cobj,
    u10cubed_p: cobj,
):
    state_u = Ptr[float](state_u_p)
    state_v = Ptr[float](state_v_p)
    state_zm = Ptr[float](state_zm_p)
    u10cubed = Ptr[float](u10cubed_p)

    for i in range(1, ncol + 1):
        wind = sqrt(state_u[_idx2(i, pver, pcols)] ** 2 + state_v[_idx2(i, pver, pcols)] ** 2)
        wind = wind * log(10.0 / z0) / log(state_zm[_idx2(i, pver, pcols)] / z0)
        u10cubed[i - 1] = wind ** 3.41


@export
def chem_timestep_init_should_run_codon(
    nstep: int,
    chem_freq: int,
    chem_step_flag_p: cobj,
):
    chem_step_flag = Ptr[int](chem_step_flag_p)
    chem_step_flag[0] = 1 if nstep % chem_freq == 0 else 0


@export
def chem_timestep_tend_fill_cloud_fields_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    ixcldliq: int,
    ixcldice: int,
    ixndrop: int,
    state_q_p: cobj,
    cldw_p: cobj,
    ncldwtr_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    cldw = Ptr[float](cldw_p)
    ncldwtr = Ptr[float](ncldwtr_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            cldw[_idx2(i, k, pcols)] = state_q[_idx3(i, k, ixcldliq, pcols, pver)] + state_q[
                _idx3(i, k, ixcldice, pcols, pver)
            ]
            if ixndrop > 0:
                ncldwtr[_idx2(i, k, pcols)] = state_q[
                    _idx3(i, k, ixndrop, pcols, pver)
                ]


@export
def chem_timestep_tend_init_lq_codon(
    pcnst: int,
    ghg_chem: int,
    map2chm_p: cobj,
    lq_mask_p: cobj,
):
    map2chm = Ptr[int](map2chm_p)
    lq_mask = Ptr[int](lq_mask_p)

    for n in range(1, pcnst + 1):
        lq_mask[n - 1] = 1 if map2chm[n - 1] > 0 else 0

    if ghg_chem != 0 and pcnst > 0:
        lq_mask[0] = 1


@export
def chem_timestep_tend_apply_depflux_codon(
    ncol: int,
    pcols: int,
    idx_cb1: int,
    idx_cb2: int,
    idx_oc1: int,
    idx_oc2: int,
    drydepflx_p: cobj,
    bcphodry_p: cobj,
    bcphidry_p: cobj,
    ocphodry_p: cobj,
    ocphidry_p: cobj,
):
    drydepflx = Ptr[float](drydepflx_p)
    bcphodry = Ptr[float](bcphodry_p)
    bcphidry = Ptr[float](bcphidry_p)
    ocphodry = Ptr[float](ocphodry_p)
    ocphidry = Ptr[float](ocphidry_p)

    if idx_cb1 > 0:
        for i in range(1, ncol + 1):
            bcphodry[i - 1] = max(drydepflx[_flux_idx(i, idx_cb1, pcols)], 0.0)

    if idx_cb2 > 0:
        for i in range(1, ncol + 1):
            bcphidry[i - 1] = max(drydepflx[_flux_idx(i, idx_cb2, pcols)], 0.0)

    if idx_oc1 > 0:
        for i in range(1, ncol + 1):
            ocphodry[i - 1] = max(drydepflx[_flux_idx(i, idx_oc1, pcols)], 0.0)

    if idx_oc2 > 0:
        for i in range(1, ncol + 1):
            ocphidry[i - 1] = max(drydepflx[_flux_idx(i, idx_oc2, pcols)], 0.0)


@export
def chem_timestep_tend_sum_fh2o_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gravit: float,
    ptend_q1_p: cobj,
    pdel_p: cobj,
    fh2o_p: cobj,
):
    ptend_q1 = Ptr[float](ptend_q1_p)
    pdel = Ptr[float](pdel_p)
    fh2o = Ptr[float](fh2o_p)

    for i in range(1, ncol + 1):
        total = 0.0
        for k in range(1, pver + 1):
            total += ptend_q1[_idx2(i, k, pcols)] * pdel[_idx2(i, k, pcols)] / gravit
        fh2o[i - 1] = total


@export
def gas_phase_chemdr_finalize_tendencies_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    pcnst: int,
    delt_inverse: float,
    map2chm_p: cobj,
    mmr_p: cobj,
    mmr_tend_p: cobj,
    mmr_new_p: cobj,
    qtend_p: cobj,
):
    map2chm = Ptr[int](map2chm_p)
    mmr = Ptr[float](mmr_p)
    mmr_tend = Ptr[float](mmr_tend_p)
    mmr_new = Ptr[float](mmr_new_p)
    qtend = Ptr[float](qtend_p)

    for m in range(1, gas_pcnst + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx3(i, k, m, pcols, pver)
                mmr_new[idx] = mmr_tend[idx]
                mmr_tend[idx] = (mmr_tend[idx] - mmr[idx]) * delt_inverse

    for m in range(1, pcnst + 1):
        n = map2chm[m - 1]
        if n > 0:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    qtend[_idx3(i, k, m, pcols, pver)] += mmr_tend[_idx3(i, k, n, pcols, pver)]


@export
def gas_phase_chemdr_prepare_state_codon(
    ncol: int,
    pcols: int,
    pver: int,
    rga: float,
    m2km: float,
    pa2mb: float,
    phis_p: cobj,
    zi_p: cobj,
    zm_p: cobj,
    pmid_p: cobj,
    zsurf_p: cobj,
    zintr_p: cobj,
    zmidr_p: cobj,
    zmid_p: cobj,
    zint_p: cobj,
    pmb_p: cobj,
):
    phis = Ptr[float](phis_p)
    zi = Ptr[float](zi_p)
    zm = Ptr[float](zm_p)
    pmid = Ptr[float](pmid_p)
    zsurf = Ptr[float](zsurf_p)
    zintr = Ptr[float](zintr_p)
    zmidr = Ptr[float](zmidr_p)
    zmid = Ptr[float](zmid_p)
    zint = Ptr[float](zint_p)
    pmb = Ptr[float](pmb_p)

    for i in range(1, ncol + 1):
        zsurf[i - 1] = rga * phis[i - 1]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            zi_in_idx = _idx2(i, k, pcols)
            zm_in_idx = _idx2(i, k, pcols)
            out_idx = _idx2(i, k, ncol)
            zsurf_val = zsurf[i - 1]
            zintr[out_idx] = m2km * zi[zi_in_idx]
            zmidr[out_idx] = m2km * zm[zm_in_idx]
            zmid[out_idx] = m2km * (zm[zm_in_idx] + zsurf_val)
            zint[out_idx] = m2km * (zi[zi_in_idx] + zsurf_val)
            pmb[out_idx] = pa2mb * pmid[zm_in_idx]

    for i in range(1, ncol + 1):
        zi_in_idx = _idx2(i, pver + 1, pcols)
        zi_out_idx = _idx2(i, pver + 1, ncol)
        zint[zi_out_idx] = m2km * (zi[zi_in_idx] + zsurf[i - 1])
        zintr[zi_out_idx] = m2km * zi[zi_in_idx]


@export
def gas_phase_chemdr_load_mmr_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    map2chm_p: cobj,
    q_p: cobj,
    mmr_p: cobj,
):
    map2chm = Ptr[int](map2chm_p)
    q = Ptr[float](q_p)
    mmr = Ptr[float](mmr_p)

    for m in range(1, pcnst + 1):
        n = map2chm[m - 1]
        if n > 0:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    idx_q = _idx3(i, k, m, pcols, pver)
                    idx_mmr = _idx3(i, k, n, pcols, pver)
                    mmr[idx_mmr] = q[idx_q]


@export
def gas_phase_chemdr_init_reaction_rates_codon(
    ncol: int,
    pver: int,
    rxntot: int,
    nan_value_p: cobj,
    reaction_rates_p: cobj,
):
    nan_value = Ptr[float](nan_value_p)[0]
    reaction_rates = Ptr[float](reaction_rates_p)

    for m in range(1, max(1, rxntot) + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                reaction_rates[_idx3(i, k, m, ncol, pver)] = nan_value


@export
def gas_phase_chemdr_clip_sulfate_codon(
    ncol: int,
    pcols: int,
    pver: int,
    troplev_p: cobj,
    sulfate_p: cobj,
):
    troplev = Ptr[int](troplev_p)
    sulfate = Ptr[float](sulfate_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            if k < troplev[i - 1]:
                sulfate[_idx2(i, k, ncol)] = 0.0


@export
def gas_phase_chemdr_zero_het_rates_codon(
    ncol: int,
    pver: int,
    gas_pcnst_dim: int,
    het_rates_p: cobj,
):
    het_rates = Ptr[float](het_rates_p)

    for m in range(1, gas_pcnst_dim + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                het_rates[_idx3(i, k, m, ncol, pver)] = 0.0


@export
def gas_phase_chemdr_load_oxygen_mmr_codon(
    ncol: int,
    pcols: int,
    pver: int,
    o2_ndx: int,
    o_ndx: int,
    mmr_p: cobj,
    o2mmr_p: cobj,
    ommr_p: cobj,
):
    mmr = Ptr[float](mmr_p)
    o2mmr = Ptr[float](o2mmr_p)
    ommr = Ptr[float](ommr_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            o2mmr[_idx2(i, k, ncol)] = mmr[_idx3(i, k, o2_ndx, pcols, pver)]
            ommr[_idx2(i, k, ncol)] = mmr[_idx3(i, k, o_ndx, pcols, pver)]


@export
def gas_phase_chemdr_set_ltrop_sol_codon(
    ncol: int,
    has_linoz_data_flag: int,
    troplev_p: cobj,
    ltrop_sol_p: cobj,
):
    troplev = Ptr[int](troplev_p)
    ltrop_sol = Ptr[int](ltrop_sol_p)

    if has_linoz_data_flag != 0:
        for i in range(1, ncol + 1):
            ltrop_sol[i - 1] = troplev[i - 1]
    else:
        for i in range(1, ncol + 1):
            ltrop_sol[i - 1] = 0


@export
def gas_phase_chemdr_zero_st80_tau_codon(
    ncol: int,
    pver: int,
    rxntot: int,
    st80_25_tau_ndx: int,
    troplev_p: cobj,
    reaction_rates_p: cobj,
):
    troplev = Ptr[int](troplev_p)
    reaction_rates = Ptr[float](reaction_rates_p)

    if st80_25_tau_ndx > 0:
        for i in range(1, ncol + 1):
            for k in range(1, troplev[i - 1] + 1):
                reaction_rates[_idx3(i, k, st80_25_tau_ndx, ncol, pver)] = 0.0


@export
def gas_phase_chemdr_compute_relhum_codon(
    ncol: int,
    pver: int,
    h2ovmr_p: cobj,
    satq_p: cobj,
    relhum_p: cobj,
):
    h2ovmr = Ptr[float](h2ovmr_p)
    satq = Ptr[float](satq_p)
    relhum = Ptr[float](relhum_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, ncol)
            value = 0.622 * h2ovmr[idx]
            value = value / satq[idx]
            if value < 0.0:
                value = 0.0
            elif value > 1.0:
                value = 1.0
            relhum[idx] = value


@export
def gas_phase_chemdr_restore_strat_gases_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    hno3_ndx: int,
    h2o_ndx: int,
    delt_inverse: float,
    vmr_p: cobj,
    hno3_gas_p: cobj,
    h2o_gas_p: cobj,
    h2ovmr_p: cobj,
    wrk_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    hno3_gas = Ptr[float](hno3_gas_p)
    h2o_gas = Ptr[float](h2o_gas_p)
    h2ovmr = Ptr[float](h2ovmr_p)
    wrk = Ptr[float](wrk_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, ncol)
            vmr[_idx3(i, k, hno3_ndx, ncol, pver)] = hno3_gas[idx2]
            h2ovmr[idx2] = h2o_gas[idx2]
            vmr[_idx3(i, k, h2o_ndx, ncol, pver)] = h2o_gas[idx2]
            wrk[idx2] = (h2ovmr[idx2] - wrk[idx2]) * delt_inverse


@export
def gas_phase_chemdr_restore_hcl_gas_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    hcl_ndx: int,
    vmr_p: cobj,
    hcl_gas_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    hcl_gas = Ptr[float](hcl_gas_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, ncol)
            vmr[_idx3(i, k, hcl_ndx, ncol, pver)] = hcl_gas[idx2]


@export
def gas_phase_chemdr_init_dust_vmr_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    ndust: int,
    dst_ndx: int,
    vmr_p: cobj,
    dust_vmr_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    dust_vmr = Ptr[float](dust_vmr_p)

    if dst_ndx > 0:
        for m in range(1, ndust + 1):
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    dust_vmr[_idx3(i, k, m, ncol, pver)] = vmr[_idx3(i, k, dst_ndx + m - 1, ncol, pver)]
    else:
        for m in range(1, ndust + 1):
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    dust_vmr[_idx3(i, k, m, ncol, pver)] = 0.0


@export
def gas_phase_chemdr_reset_ste_tracer_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    st80_25_ndx: int,
    pmid_threshold: float,
    st80_vmr: float,
    pmid_p: cobj,
    vmr_p: cobj,
):
    pmid = Ptr[float](pmid_p)
    vmr = Ptr[float](vmr_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            if pmid[_idx2(i, k, pcols)] < pmid_threshold:
                vmr[_idx3(i, k, st80_25_ndx, ncol, pver)] = st80_vmr


@export
def gas_phase_chemdr_zero_sflx_codon(
    pcols: int,
    gas_pcnst: int,
    sflx_p: cobj,
):
    sflx = Ptr[float](sflx_p)

    for m in range(1, gas_pcnst + 1):
        for i in range(1, pcols + 1):
            sflx[_flux_idx(i, m, pcols)] = 0.0


@export
def gas_phase_chemdr_compute_wind_speed_codon(
    ncol: int,
    pcols: int,
    pver: int,
    ufld_p: cobj,
    vfld_p: cobj,
    wind_speed_p: cobj,
):
    ufld = Ptr[float](ufld_p)
    vfld = Ptr[float](vfld_p)
    wind_speed = Ptr[float](wind_speed_p)

    for i in range(1, ncol + 1):
        uval = ufld[_idx2(i, pver, pcols)]
        vval = vfld[_idx2(i, pver, pcols)]
        wind_speed[i - 1] = sqrt(uval * uval + vval * vval)


@export
def gas_phase_chemdr_compute_prect_codon(
    ncol: int,
    pcols: int,
    precc_p: cobj,
    precl_p: cobj,
    prect_p: cobj,
):
    precc = Ptr[float](precc_p)
    precl = Ptr[float](precl_p)
    prect = Ptr[float](prect_p)

    for i in range(1, ncol + 1):
        prect[i - 1] = precc[i - 1] + precl[i - 1]


@export
def gas_phase_chemdr_compute_tvs_codon(
    ncol: int,
    pcols: int,
    pver: int,
    tfld_p: cobj,
    qh2o_p: cobj,
    tvs_p: cobj,
):
    tfld = Ptr[float](tfld_p)
    qh2o = Ptr[float](qh2o_p)
    tvs = Ptr[float](tvs_p)

    for i in range(1, ncol + 1):
        tvs[i - 1] = tfld[_idx2(i, pver, pcols)] * (1.0 + qh2o[_idx2(i, pver, pcols)])


@export
def gas_phase_chemdr_copy_cldw_to_cwat_codon(
    ncol: int,
    pcols: int,
    pver: int,
    cldw_p: cobj,
    cwat_p: cobj,
):
    cldw = Ptr[float](cldw_p)
    cwat = Ptr[float](cwat_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            cwat[_idx2(i, k, ncol)] = cldw[_idx2(i, k, pcols)]


@export
def gas_phase_chemdr_load_h2o_fields_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    h2o_ndx: int,
    mmr_p: cobj,
    vmr_p: cobj,
    qh2o_p: cobj,
    h2ovmr_p: cobj,
):
    mmr = Ptr[float](mmr_p)
    vmr = Ptr[float](vmr_p)
    qh2o = Ptr[float](qh2o_p)
    h2ovmr = Ptr[float](h2ovmr_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            qh2o[_idx2(i, k, pcols)] = mmr[_idx3(i, k, h2o_ndx, pcols, pver)]
            h2ovmr[_idx2(i, k, ncol)] = vmr[_idx3(i, k, h2o_ndx, ncol, pver)]


@export
def gas_phase_chemdr_copy_o3_to_o3s_trop_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    troplev_p: cobj,
    o3_ndx: int,
    o3s_ndx: int,
    vmr_p: cobj,
):
    troplev = Ptr[int](troplev_p)
    vmr = Ptr[float](vmr_p)

    for i in range(1, ncol + 1):
        for k in range(1, troplev[i - 1] + 1):
            vmr[_idx3(i, k, o3s_ndx, ncol, pver)] = vmr[_idx3(i, k, o3_ndx, ncol, pver)]


@export
def gas_phase_chemdr_copy_h2o_to_wrk_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    h2o_ndx: int,
    vmr_p: cobj,
    wrk_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    wrk = Ptr[float](wrk_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            wrk[_idx2(i, k, ncol)] = vmr[_idx3(i, k, h2o_ndx, ncol, pver)]


@export
def gas_phase_chemdr_update_qdsett_wrk_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    h2o_ndx: int,
    delt_inverse: float,
    vmr_p: cobj,
    wrk_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    wrk = Ptr[float](wrk_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, ncol)
            wrk[idx2] = (vmr[_idx3(i, k, h2o_ndx, ncol, pver)] - wrk[idx2]) * delt_inverse


@export
def gas_phase_chemdr_update_qdchem_wrk_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    h2o_ndx: int,
    delt_inverse: float,
    vmr_p: cobj,
    wrk_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    wrk = Ptr[float](wrk_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, ncol)
            wrk[idx2] = (vmr[_idx3(i, k, h2o_ndx, ncol, pver)] - wrk[idx2]) * delt_inverse


@export
def gas_phase_chemdr_init_stratchem_state_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    pcnst: int,
    hno3_ndx: int,
    hcl_ndx: int,
    cldice_ndx: int,
    vmr_p: cobj,
    h2ovmr_p: cobj,
    q_p: cobj,
    hcl_cond_p: cobj,
    hcl_gas_p: cobj,
    hno3_gas_p: cobj,
    h2o_gas_p: cobj,
    wrk_p: cobj,
    cldice_p: cobj,
    hno3_cond_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    h2ovmr = Ptr[float](h2ovmr_p)
    q = Ptr[float](q_p)
    hcl_cond = Ptr[float](hcl_cond_p)
    hcl_gas = Ptr[float](hcl_gas_p)
    hno3_gas = Ptr[float](hno3_gas_p)
    h2o_gas = Ptr[float](h2o_gas_p)
    wrk = Ptr[float](wrk_p)
    cldice = Ptr[float](cldice_p)
    hno3_cond = Ptr[float](hno3_cond_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, ncol)
            hcl_cond[idx2] = 0.0
            hno3_cond[_idx3(i, k, 1, ncol, pver)] = 0.0
            hno3_cond[_idx3(i, k, 2, ncol, pver)] = 0.0
            hno3_gas[idx2] = vmr[_idx3(i, k, hno3_ndx, ncol, pver)]
            h2o_gas[idx2] = h2ovmr[idx2]
            hcl_gas[idx2] = vmr[_idx3(i, k, hcl_ndx, ncol, pver)]
            wrk[idx2] = h2ovmr[idx2]
            cldice[_idx2(i, k, pcols)] = q[_idx3(i, k, cldice_ndx, pcols, pver)]


@export
def gas_phase_chemdr_init_h2so4_gasprod_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    ndx_h2so4: int,
    vmr_p: cobj,
    del_h2so4_gasprod_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    del_h2so4_gasprod = Ptr[float](del_h2so4_gasprod_p)

    if ndx_h2so4 > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                del_h2so4_gasprod[_idx2(i, k, ncol)] = vmr[_idx3(i, k, ndx_h2so4, ncol, pver)]
    else:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                del_h2so4_gasprod[_idx2(i, k, ncol)] = 0.0


@export
def gas_phase_chemdr_store_vmr0_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    vmr_p: cobj,
    vmr0_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    vmr0 = Ptr[float](vmr0_p)

    for m in range(1, gas_pcnst + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx3(i, k, m, ncol, pver)
                vmr0[idx] = vmr[idx]


@export
def gas_phase_chemdr_update_h2so4_gasprod_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    ndx_h2so4: int,
    vmr_p: cobj,
    del_h2so4_gasprod_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    del_h2so4_gasprod = Ptr[float](del_h2so4_gasprod_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, ncol)
            del_h2so4_gasprod[idx2] = vmr[_idx3(i, k, ndx_h2so4, ncol, pver)] - del_h2so4_gasprod[idx2]


@export
def gas_phase_chemdr_reform_hno3_hcl_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    hno3_ndx: int,
    hcl_ndx: int,
    vmr_p: cobj,
    hno3_cond_p: cobj,
    hcl_cond_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    hno3_cond = Ptr[float](hno3_cond_p)
    hcl_cond = Ptr[float](hcl_cond_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            hno3_idx = _idx3(i, k, hno3_ndx, ncol, pver)
            hno3_total = vmr[hno3_idx] + hno3_cond[_idx3(i, k, 1, ncol, pver)]
            vmr[hno3_idx] = hno3_total + hno3_cond[_idx3(i, k, 2, ncol, pver)]

            hcl_idx = _idx3(i, k, hcl_ndx, ncol, pver)
            vmr[hcl_idx] = vmr[hcl_idx] + hcl_cond[_idx2(i, k, ncol)]


@export
def gas_phase_chemdr_normalize_extfrc_codon(
    ncol: int,
    pver: int,
    extcnt: int,
    synoz_ndx: int,
    aoa_nh_ext_ndx: int,
    indexm: int,
    extfrc_p: cobj,
    invariants_p: cobj,
):
    extfrc = Ptr[float](extfrc_p)
    invariants = Ptr[float](invariants_p)

    for m in range(1, extcnt + 1):
        if m != synoz_ndx and m != aoa_nh_ext_ndx:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    extfrc[_idx3(i, k, m, ncol, pver)] = extfrc[_idx3(i, k, m, ncol, pver)] / invariants[
                        _idx3(i, k, indexm, ncol, pver)
                    ]


@export
def gas_phase_chemdr_store_drydep_codon(
    ncol: int,
    pcols: int,
    gas_pcnst: int,
    pcnst: int,
    map2chm_p: cobj,
    sflx_p: cobj,
    cflx_p: cobj,
    drydepflx_p: cobj,
):
    map2chm = Ptr[int](map2chm_p)
    sflx = Ptr[float](sflx_p)
    cflx = Ptr[float](cflx_p)
    drydepflx = Ptr[float](drydepflx_p)

    for m in range(1, pcnst + 1):
        for i in range(1, pcols + 1):
            drydepflx[_flux_idx(i, m, pcols)] = 0.0

    for m in range(1, pcnst + 1):
        n = map2chm[m - 1]
        if n > 0:
            for i in range(1, ncol + 1):
                src_idx = _flux_idx(i, n, pcols)
                dst_idx = _flux_idx(i, m, pcols)
                cflx[dst_idx] = cflx[dst_idx] - sflx[src_idx]
                drydepflx[dst_idx] = sflx[src_idx]


@export
def set_xnox_photo_codon(
    ncol: int,
    pver: int,
    photos_p: cobj,
    jno2a_ndx: int,
    jno2_ndx: int,
    jn2o5a_ndx: int,
    jn2o5_ndx: int,
    jn2o5b_ndx: int,
    jhno3a_ndx: int,
    jhno3_ndx: int,
    jno3a_ndx: int,
    jno3_ndx: int,
    jho2no2a_ndx: int,
    jho2no2_ndx: int,
    jmpana_ndx: int,
    jmpan_ndx: int,
    jpana_ndx: int,
    jpan_ndx: int,
    jonitra_ndx: int,
    jonitr_ndx: int,
    jo1da_ndx: int,
    jo1d_ndx: int,
    jo3pa_ndx: int,
    jo3p_ndx: int,
):
    photos = Ptr[float](photos_p)

    if jno2a_ndx > 0 and jno2_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                photos[_idx3(i, k, jno2a_ndx, ncol, pver)] = photos[_idx3(i, k, jno2_ndx, ncol, pver)]

    if jn2o5a_ndx > 0 and jn2o5_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                photos[_idx3(i, k, jn2o5a_ndx, ncol, pver)] = photos[_idx3(i, k, jn2o5_ndx, ncol, pver)]

    if jn2o5b_ndx > 0 and jn2o5_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                photos[_idx3(i, k, jn2o5b_ndx, ncol, pver)] = photos[_idx3(i, k, jn2o5_ndx, ncol, pver)]

    if jhno3a_ndx > 0 and jhno3_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                photos[_idx3(i, k, jhno3a_ndx, ncol, pver)] = photos[_idx3(i, k, jhno3_ndx, ncol, pver)]

    if jno3a_ndx > 0 and jno3_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                photos[_idx3(i, k, jno3a_ndx, ncol, pver)] = photos[_idx3(i, k, jno3_ndx, ncol, pver)]

    if jho2no2a_ndx > 0 and jho2no2_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                photos[_idx3(i, k, jho2no2a_ndx, ncol, pver)] = photos[_idx3(i, k, jho2no2_ndx, ncol, pver)]

    if jmpana_ndx > 0 and jmpan_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                photos[_idx3(i, k, jmpana_ndx, ncol, pver)] = photos[_idx3(i, k, jmpan_ndx, ncol, pver)]

    if jpana_ndx > 0 and jpan_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                photos[_idx3(i, k, jpana_ndx, ncol, pver)] = photos[_idx3(i, k, jpan_ndx, ncol, pver)]

    if jonitra_ndx > 0 and jonitr_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                photos[_idx3(i, k, jonitra_ndx, ncol, pver)] = photos[_idx3(i, k, jonitr_ndx, ncol, pver)]

    if jo1da_ndx > 0 and jo1d_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                photos[_idx3(i, k, jo1da_ndx, ncol, pver)] = photos[_idx3(i, k, jo1d_ndx, ncol, pver)]

    if jo3pa_ndx > 0 and jo3p_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                photos[_idx3(i, k, jo3pa_ndx, ncol, pver)] = photos[_idx3(i, k, jo3p_ndx, ncol, pver)]


@export
def adjrxt_codon(
    ncol: int,
    pver: int,
    rate_p: cobj,
    inv_p: cobj,
    m_p: cobj,
):
    rate = Ptr[float](rate_p)
    inv = Ptr[float](inv_p)
    m = Ptr[float](m_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx_r2 = _idx3(i, k, 2, ncol, pver)
            idx_r3 = _idx3(i, k, 3, ncol, pver)
            idx_r4 = _idx3(i, k, 4, ncol, pver)
            idx_r5 = _idx3(i, k, 5, ncol, pver)
            idx_r6 = _idx3(i, k, 6, ncol, pver)
            idx_r7 = _idx3(i, k, 7, ncol, pver)
            idx_i6 = _idx3(i, k, 6, ncol, pver)
            idx_i7 = _idx3(i, k, 7, ncol, pver)
            idx_i8 = _idx3(i, k, 8, ncol, pver)
            idx_m = _idx2(i, k, ncol)

            inv6 = inv[idx_i6]
            inv7 = inv[idx_i7]
            inv8 = inv[idx_i8]
            im = 1.0 / m[idx_m]

            rate[idx_r3] = rate[idx_r3] * inv6
            rate[idx_r4] = rate[idx_r4] * inv6
            rate[idx_r5] = rate[idx_r5] * inv6
            rate[idx_r6] = rate[idx_r6] * inv6
            rate[idx_r7] = rate[idx_r7] * inv7

            tmp = rate[idx_r2] * inv8
            tmp = tmp * inv8
            rate[idx_r2] = tmp * im


@export
def setrxt_codon(
    ncol: int,
    pcols: int,
    pver: int,
    temp_p: cobj,
    rate_p: cobj,
):
    temp = Ptr[float](temp_p)
    rate = Ptr[float](rate_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            itemp = 1.0 / temp[_idx2(i, k, pcols)]
            rate[_idx3(i, k, 3, ncol, pver)] = 2.9e-12 * exp(-160.0 * itemp)
            rate[_idx3(i, k, 5, ncol, pver)] = 9.6e-12 * exp(-234.0 * itemp)
            rate[_idx3(i, k, 7, ncol, pver)] = 1.9e-13 * exp(520.0 * itemp)


@export
def lu_slv_codon(
    lu_p: cobj,
    b_p: cobj,
):
    lu = Ptr[float](lu_p)
    b = Ptr[float](b_p)

    b[19] = b[19] * lu[21]
    b[18] = b[18] * lu[20]
    b[17] = b[17] * lu[19]
    b[16] = b[16] * lu[18]
    b[15] = b[15] * lu[17]
    b[14] = b[14] * lu[16]
    b[13] = b[13] * lu[15]
    b[12] = b[12] * lu[14]
    b[11] = b[11] * lu[13]
    b[10] = b[10] * lu[12]
    b[9] = b[9] * lu[11]
    b[8] = b[8] * lu[10]
    b[7] = b[7] * lu[9]
    b[6] = b[6] * lu[8]
    b[5] = b[5] * lu[7]
    b[4] = b[4] * lu[6]
    b[3] = b[3] * lu[5]
    b[2] = b[2] - lu[4] * b[3]
    b[2] = b[2] * lu[3]
    b[1] = b[1] - lu[2] * b[2]
    b[1] = b[1] * lu[1]
    b[0] = b[0] * lu[0]


@export
def lu_fac_codon(
    lu_p: cobj,
):
    lu = Ptr[float](lu_p)

    lu[0] = 1.0 / lu[0]
    lu[1] = 1.0 / lu[1]
    lu[3] = 1.0 / lu[3]
    lu[5] = 1.0 / lu[5]
    lu[6] = 1.0 / lu[6]
    lu[7] = 1.0 / lu[7]
    lu[8] = 1.0 / lu[8]
    lu[9] = 1.0 / lu[9]
    lu[10] = 1.0 / lu[10]
    lu[11] = 1.0 / lu[11]
    lu[12] = 1.0 / lu[12]
    lu[13] = 1.0 / lu[13]
    lu[14] = 1.0 / lu[14]
    lu[15] = 1.0 / lu[15]
    lu[16] = 1.0 / lu[16]
    lu[17] = 1.0 / lu[17]
    lu[18] = 1.0 / lu[18]
    lu[19] = 1.0 / lu[19]
    lu[20] = 1.0 / lu[20]
    lu[21] = 1.0 / lu[21]


@export
def linmat_codon(
    mat_p: cobj,
    rxt_p: cobj,
    het_rates_p: cobj,
):
    mat = Ptr[float](mat_p)
    rxt = Ptr[float](rxt_p)
    het_rates = Ptr[float](het_rates_p)

    mat[0] = -(rxt[0] + rxt[2] + het_rates[0])
    mat[1] = -(het_rates[1])
    mat[2] = rxt[3]
    mat[3] = -(rxt[3] + het_rates[2])
    mat[4] = rxt[4] + 0.5 * rxt[5] + rxt[6]
    mat[5] = -(rxt[4] + rxt[5] + rxt[6] + het_rates[3])
    mat[6] = -(het_rates[4])
    mat[7] = -(het_rates[5])
    mat[8] = -(het_rates[6])
    mat[9] = -(het_rates[7])
    mat[10] = -(het_rates[8])
    mat[11] = -(het_rates[9])
    mat[12] = -(het_rates[10])
    mat[13] = -(het_rates[11])
    mat[14] = -(het_rates[12])
    mat[15] = -(het_rates[13])
    mat[16] = -(het_rates[14])
    mat[17] = -(het_rates[15])
    mat[18] = -(het_rates[16])
    mat[19] = -(het_rates[17])
    mat[20] = -(het_rates[18])
    mat[21] = -(het_rates[19])


@export
def imp_prod_loss_codon(
    prod_p: cobj,
    loss_p: cobj,
    y_p: cobj,
    rxt_p: cobj,
    het_rates_p: cobj,
):
    prod = Ptr[float](prod_p)
    loss = Ptr[float](loss_p)
    y = Ptr[float](y_p)
    rxt = Ptr[float](rxt_p)
    het_rates = Ptr[float](het_rates_p)

    loss[0] = (rxt[0] + rxt[2] + het_rates[0]) * y[0]
    prod[0] = 0.0
    loss[1] = het_rates[1] * y[1]
    prod[1] = rxt[3] * y[2]
    loss[2] = (rxt[3] + het_rates[2]) * y[2]
    prod[2] = (rxt[4] + 0.5 * rxt[5] + rxt[6]) * y[3]
    loss[3] = (rxt[4] + rxt[5] + rxt[6] + het_rates[3]) * y[3]
    prod[3] = 0.0
    loss[4] = het_rates[4] * y[4]
    prod[4] = 0.0
    loss[5] = het_rates[5] * y[5]
    prod[5] = 0.0
    loss[6] = het_rates[6] * y[6]
    prod[6] = 0.0
    loss[7] = het_rates[7] * y[7]
    prod[7] = 0.0
    loss[8] = het_rates[8] * y[8]
    prod[8] = 0.0
    loss[9] = het_rates[9] * y[9]
    prod[9] = 0.0
    loss[10] = het_rates[10] * y[10]
    prod[10] = 0.0
    loss[11] = het_rates[11] * y[11]
    prod[11] = 0.0
    loss[12] = het_rates[12] * y[12]
    prod[12] = 0.0
    loss[13] = het_rates[13] * y[13]
    prod[13] = 0.0
    loss[14] = het_rates[14] * y[14]
    prod[14] = 0.0
    loss[15] = het_rates[15] * y[15]
    prod[15] = 0.0
    loss[16] = het_rates[16] * y[16]
    prod[16] = 0.0
    loss[17] = het_rates[17] * y[17]
    prod[17] = 0.0
    loss[18] = het_rates[18] * y[18]
    prod[18] = 0.0
    loss[19] = het_rates[19] * y[19]
    prod[19] = 0.0


@export
def nlnmat_codon(
    mat_p: cobj,
    lmat_p: cobj,
    dti: float,
):
    mat = Ptr[float](mat_p)
    lmat = Ptr[float](lmat_p)

    mat[0] = lmat[0]
    mat[1] = lmat[1]
    mat[2] = lmat[2]
    mat[3] = lmat[3]
    mat[4] = lmat[4]
    mat[5] = lmat[5]
    mat[6] = lmat[6]
    mat[7] = lmat[7]
    mat[8] = lmat[8]
    mat[9] = lmat[9]
    mat[10] = lmat[10]
    mat[11] = lmat[11]
    mat[12] = lmat[12]
    mat[13] = lmat[13]
    mat[14] = lmat[14]
    mat[15] = lmat[15]
    mat[16] = lmat[16]
    mat[17] = lmat[17]
    mat[18] = lmat[18]
    mat[19] = lmat[19]
    mat[20] = lmat[20]
    mat[21] = lmat[21]

    mat[0] = mat[0] - dti
    mat[1] = mat[1] - dti
    mat[3] = mat[3] - dti
    mat[5] = mat[5] - dti
    mat[6] = mat[6] - dti
    mat[7] = mat[7] - dti
    mat[8] = mat[8] - dti
    mat[9] = mat[9] - dti
    mat[10] = mat[10] - dti
    mat[11] = mat[11] - dti
    mat[12] = mat[12] - dti
    mat[13] = mat[13] - dti
    mat[14] = mat[14] - dti
    mat[15] = mat[15] - dti
    mat[16] = mat[16] - dti
    mat[17] = mat[17] - dti
    mat[18] = mat[18] - dti
    mat[19] = mat[19] - dti
    mat[20] = mat[20] - dti
    mat[21] = mat[21] - dti


@export
def indprd_codon(
    class_id: int,
    ncol: int,
    pver: int,
    nprod: int,
    rxt_p: cobj,
    extfrc_p: cobj,
    prod_p: cobj,
):
    rxt = Ptr[float](rxt_p)
    extfrc = Ptr[float](extfrc_p)
    prod = Ptr[float](prod_p)

    if class_id != 4:
        return

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            prod[_idx3(i, k, 1, ncol, pver)] = rxt[_idx3(i, k, 2, ncol, pver)]
            prod[_idx3(i, k, 2, ncol, pver)] = 0.0
            prod[_idx3(i, k, 3, ncol, pver)] = extfrc[_idx3(i, k, 1, ncol, pver)]
            prod[_idx3(i, k, 4, ncol, pver)] = 0.0
            prod[_idx3(i, k, 5, ncol, pver)] = 0.0
            prod[_idx3(i, k, 6, ncol, pver)] = extfrc[_idx3(i, k, 2, ncol, pver)]
            prod[_idx3(i, k, 7, ncol, pver)] = extfrc[_idx3(i, k, 4, ncol, pver)]
            prod[_idx3(i, k, 8, ncol, pver)] = 0.0
            prod[_idx3(i, k, 9, ncol, pver)] = extfrc[_idx3(i, k, 5, ncol, pver)]
            prod[_idx3(i, k, 10, ncol, pver)] = 0.0
            prod[_idx3(i, k, 11, ncol, pver)] = 0.0
            prod[_idx3(i, k, 12, ncol, pver)] = extfrc[_idx3(i, k, 6, ncol, pver)]
            prod[_idx3(i, k, 13, ncol, pver)] = extfrc[_idx3(i, k, 3, ncol, pver)]
            prod[_idx3(i, k, 14, ncol, pver)] = 0.0
            prod[_idx3(i, k, 15, ncol, pver)] = 0.0
            prod[_idx3(i, k, 16, ncol, pver)] = extfrc[_idx3(i, k, 7, ncol, pver)]
            prod[_idx3(i, k, 17, ncol, pver)] = 0.0
            prod[_idx3(i, k, 18, ncol, pver)] = 0.0
            prod[_idx3(i, k, 19, ncol, pver)] = 0.0
            prod[_idx3(i, k, 20, ncol, pver)] = 0.0


@export
def negtrc_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    fld_p: cobj,
):
    fld = Ptr[float](fld_p)

    for m in range(1, gas_pcnst + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx3(i, k, m, ncol, pver)
                if fld[idx] < 0.0:
                    fld[idx] = 0.0


@export
def O1D_to_2OH_adj_codon(
    ncol: int,
    pcols: int,
    pver: int,
    rxntot: int,
    nfs: int,
    jo1d_ndx: int,
    n2_ndx: int,
    o2_ndx: int,
    h2o_ndx: int,
    p_rate_p: cobj,
    inv_p: cobj,
    tfld_p: cobj,
):
    if jo1d_ndx < 1:
        return

    p_rate = Ptr[float](p_rate_p)
    inv = Ptr[float](inv_p)
    tfld = Ptr[float](tfld_p)

    x1 = 2.15e-11
    x2 = 3.30e-11
    x3 = 1.63e-10
    y1 = 110.0
    y2 = 55.0
    y3 = 60.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            temp = tfld[_idx2(i, k, pcols)]
            n2_rate = x1 * exp(y1 / temp) * inv[_idx3(i, k, n2_ndx, ncol, pver)]
            o2_rate = x2 * exp(y2 / temp) * inv[_idx3(i, k, o2_ndx, ncol, pver)]
            h2o_rate = x3 * exp(y3 / temp) * inv[_idx3(i, k, h2o_ndx, ncol, pver)]
            denom = h2o_rate + n2_rate + o2_rate
            p_rate[_idx3(i, k, jo1d_ndx, ncol, pver)] = (
                p_rate[_idx3(i, k, jo1d_ndx, ncol, pver)] * (h2o_rate / denom)
            )
