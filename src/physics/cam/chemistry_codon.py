from math import acos, cos, exp, log, log10, sin, sqrt
from C import neu_wetdep_dempirical_native_cb(float, float) -> float
from C import neu_wetdep_gamma_native_cb(float) -> float


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


@inline
def _aero_model_gasaerexch_column_flux(
    ncol: int,
    pcols: int,
    pver: int,
    adv_mass: float,
    gravit: float,
    field: Ptr[float],
    mbar: Ptr[float],
    pdel: Ptr[float],
    wrk: Ptr[float],
):
    for i in range(1, ncol + 1):
        wrk[i - 1] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            wrk[i - 1] += field[_idx2(i, k, ncol)] * adv_mass / mbar[_idx2(i, k, pcols)] * pdel[
                _idx2(i, k, pcols)
            ] / gravit


@inline
def _aero_model_gasaerexch_all_column_fluxes(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    gravit: float,
    field: Ptr[float],
    mbar: Ptr[float],
    pdel: Ptr[float],
    adv_mass: Ptr[float],
    wrk: Ptr[float],
):
    for m in range(1, gas_pcnst + 1):
        mass = adv_mass[m - 1]
        for i in range(1, ncol + 1):
            wrk[_idx2(i, m, ncol)] = 0.0

        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                wrk[_idx2(i, m, ncol)] += (
                    field[_idx3(i, k, m, ncol, pver)]
                    * mass
                    / mbar[_idx2(i, k, pcols)]
                    * pdel[_idx2(i, k, pcols)]
                    / gravit
                )


@inline
def _aero_model_gasaerexch_h2so4_save_or_delta(
    ncol: int,
    pver: int,
    ndx_h2so4: int,
    stage3_mode: int,
    vmr: Ptr[float],
    del_h2so4_aeruptk: Ptr[float],
):
    if stage3_mode == 0:
        if ndx_h2so4 > 0:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    del_h2so4_aeruptk[_idx2(i, k, ncol)] = vmr[_idx3(i, k, ndx_h2so4, ncol, pver)]
        else:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    del_h2so4_aeruptk[_idx2(i, k, ncol)] = 0.0
        return

    if ndx_h2so4 > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, ncol)
                del_h2so4_aeruptk[idx] = vmr[_idx3(i, k, ndx_h2so4, ncol, pver)] - del_h2so4_aeruptk[idx]


@inline
def _aero_model_gasaerexch_gas_tend(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    delt: float,
    vmr0: Ptr[float],
    vmr: Ptr[float],
    dvmrdt: Ptr[float],
):
    for m in range(1, gas_pcnst + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx3(i, k, m, ncol, pver)
                dvmrdt[idx] = (vmr[idx] - vmr0[idx]) / delt


@inline
def _aero_model_gasaerexch_aq_tend(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    delt: float,
    vmr: Ptr[float],
    vmrcw: Ptr[float],
    dvmrdt: Ptr[float],
    dvmrcwdt: Ptr[float],
):
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
def aero_model_gasaerexch_codon(
    stage: int,
    stage3_mode: int,
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    ndx_h2so4: int,
    delt: float,
    gravit: float,
    vmr0_p: cobj,
    vmr_p: cobj,
    vmrcw_p: cobj,
    dvmrdt_p: cobj,
    dvmrcwdt_p: cobj,
    mbar_p: cobj,
    pdel_p: cobj,
    adv_mass_p: cobj,
    wrk_p: cobj,
    del_h2so4_aeruptk_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    dvmrdt = Ptr[float](dvmrdt_p)

    if stage == 1:
        vmr0 = Ptr[float](vmr0_p)
        mbar = Ptr[float](mbar_p)
        pdel = Ptr[float](pdel_p)
        adv_mass = Ptr[float](adv_mass_p)
        wrk = Ptr[float](wrk_p)

        _aero_model_gasaerexch_gas_tend(ncol, pver, gas_pcnst, delt, vmr0, vmr, dvmrdt)
        _aero_model_gasaerexch_all_column_fluxes(ncol, pcols, pver, gas_pcnst, gravit, dvmrdt, mbar, pdel, adv_mass, wrk)
        return

    if stage == 2:
        vmrcw = Ptr[float](vmrcw_p)
        dvmrcwdt = Ptr[float](dvmrcwdt_p)
        mbar = Ptr[float](mbar_p)
        pdel = Ptr[float](pdel_p)
        adv_mass = Ptr[float](adv_mass_p)
        wrk = Ptr[float](wrk_p)

        _aero_model_gasaerexch_aq_tend(ncol, pver, gas_pcnst, delt, vmr, vmrcw, dvmrdt, dvmrcwdt)
        _aero_model_gasaerexch_all_column_fluxes(ncol, pcols, pver, gas_pcnst, gravit, dvmrdt, mbar, pdel, adv_mass, wrk)
        return

    if stage == 3:
        del_h2so4_aeruptk = Ptr[float](del_h2so4_aeruptk_p)
        _aero_model_gasaerexch_h2so4_save_or_delta(
            ncol, pver, ndx_h2so4, stage3_mode, vmr, del_h2so4_aeruptk
        )


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

    _aero_model_gasaerexch_column_flux(ncol, pcols, pver, adv_mass, gravit, field, mbar, pdel, wrk)


@export
def aero_model_gasaerexch_h2so4_save_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    ndx_h2so4: int,
    vmr_p: cobj,
    del_h2so4_aeruptk_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    del_h2so4_aeruptk = Ptr[float](del_h2so4_aeruptk_p)

    _aero_model_gasaerexch_h2so4_save_or_delta(ncol, pver, ndx_h2so4, 0, vmr, del_h2so4_aeruptk)


@export
def aero_model_gasaerexch_h2so4_delta_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    ndx_h2so4: int,
    vmr_p: cobj,
    del_h2so4_aeruptk_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    del_h2so4_aeruptk = Ptr[float](del_h2so4_aeruptk_p)

    _aero_model_gasaerexch_h2so4_save_or_delta(ncol, pver, ndx_h2so4, 1, vmr, del_h2so4_aeruptk)


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

    _aero_model_gasaerexch_gas_tend(ncol, pver, gas_pcnst, delt, vmr0, vmr, dvmrdt)


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

    _aero_model_gasaerexch_aq_tend(ncol, pver, gas_pcnst, delt, vmr, vmrcw, dvmrdt, dvmrcwdt)


@export
def neu_wetdep_aux_prepare_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    gas_cnt: int,
    index_cldice: int,
    index_cldliq: int,
    gravit: float,
    mapping_to_mmr_p: cobj,
    area_p: cobj,
    mmr_p: cobj,
    pmid_p: cobj,
    pdel_p: cobj,
    zint_p: cobj,
    tfld_p: cobj,
    prain_p: cobj,
    nevapr_p: cobj,
    cld_p: cobj,
    cmfdqr_p: cobj,
    mass_in_layer_p: cobj,
    cldice_p: cobj,
    cldliq_p: cobj,
    cldfrc_p: cobj,
    totprec_p: cobj,
    totevap_p: cobj,
    delz_p: cobj,
    delp_p: cobj,
    p_p: cobj,
    rls_p: cobj,
    evaprate_p: cobj,
    temp_p: cobj,
    trc_mass_p: cobj,
    dtwr_p: cobj,
):
    mapping_to_mmr = Ptr[int](mapping_to_mmr_p)
    area = Ptr[float](area_p)
    mmr = Ptr[float](mmr_p)
    pmid = Ptr[float](pmid_p)
    pdel = Ptr[float](pdel_p)
    zint = Ptr[float](zint_p)
    tfld = Ptr[float](tfld_p)
    prain = Ptr[float](prain_p)
    nevapr = Ptr[float](nevapr_p)
    cld = Ptr[float](cld_p)
    cmfdqr = Ptr[float](cmfdqr_p)
    mass_in_layer = Ptr[float](mass_in_layer_p)
    cldice = Ptr[float](cldice_p)
    cldliq = Ptr[float](cldliq_p)
    cldfrc = Ptr[float](cldfrc_p)
    totprec = Ptr[float](totprec_p)
    totevap = Ptr[float](totevap_p)
    delz = Ptr[float](delz_p)
    delp = Ptr[float](delp_p)
    press = Ptr[float](p_p)
    rls = Ptr[float](rls_p)
    evaprate = Ptr[float](evaprate_p)
    temp = Ptr[float](temp_p)
    trc_mass = Ptr[float](trc_mass_p)
    dtwr = Ptr[float](dtwr_p)

    for k in range(1, pver + 1):
        kk = pver - k + 1
        for i in range(1, ncol + 1):
            idx_rev_ncol = _idx2(i, k, ncol)
            idx_kk_pcols = _idx2(i, kk, pcols)
            layer_mass = area[i - 1] * pdel[idx_kk_pcols] / gravit
            mass_in_layer[idx_rev_ncol] = layer_mass

            cldice[idx_rev_ncol] = mmr[_idx3(i, kk, index_cldice, pcols, pver)]
            cldliq[idx_rev_ncol] = mmr[_idx3(i, kk, index_cldliq, pcols, pver)]
            cldfrc[idx_rev_ncol] = cld[_idx2(i, kk, ncol)]

            totprec[idx_rev_ncol] = (prain[_idx2(i, kk, ncol)] + cmfdqr[_idx2(i, kk, ncol)]) * layer_mass
            totevap[idx_rev_ncol] = nevapr[_idx2(i, kk, ncol)] * layer_mass

            delz[idx_rev_ncol] = zint[_idx2(i, kk, pcols)] - zint[_idx2(i, kk + 1, pcols)]
            temp[idx_rev_ncol] = tfld[idx_kk_pcols]

            for m in range(1, gas_cnt + 1):
                spc = mapping_to_mmr[m - 1]
                trc_mass[_idx3(i, k, m, ncol, pver)] = mmr[_idx3(i, kk, spc, pcols, pver)] * layer_mass

            delp[idx_rev_ncol] = pdel[idx_kk_pcols] * 0.01
            press[idx_rev_ncol] = pmid[idx_kk_pcols] * 0.01

    for m in range(1, gas_cnt + 1):
        spc = mapping_to_mmr[m - 1]
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                dtwr[_idx3(i, k, m, ncol, pver)] = mmr[_idx3(i, k, spc, pcols, pver)]

    for i in range(1, ncol + 1):
        rls[_idx2(i, pver, ncol)] = 0.0
        evaprate[_idx2(i, pver, ncol)] = 0.0

    for k in range(pver - 1, 0, -1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, ncol)
            next_idx = _idx2(i, k + 1, ncol)
            rls[idx] = max(0.0, totprec[idx] - totevap[idx] + rls[next_idx])
            evaprate[idx] = min(1.0, totevap[idx] / (rls[next_idx] + 1.0e-36))


@export
def neu_wetdep_aux_finish_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    gas_cnt: int,
    delt: float,
    pi: float,
    mapping_to_mmr_p: cobj,
    lats_p: cobj,
    pmid_p: cobj,
    mass_in_layer_p: cobj,
    trc_mass_p: cobj,
    dtwr_p: cobj,
    wd_mmr_p: cobj,
    wd_tend_p: cobj,
):
    mapping_to_mmr = Ptr[int](mapping_to_mmr_p)
    lats = Ptr[float](lats_p)
    pmid = Ptr[float](pmid_p)
    mass_in_layer = Ptr[float](mass_in_layer_p)
    trc_mass = Ptr[float](trc_mass_p)
    dtwr = Ptr[float](dtwr_p)
    wd_mmr = Ptr[float](wd_mmr_p)
    wd_tend = Ptr[float](wd_tend_p)

    for k in range(1, pver + 1):
        kk = pver - k + 1
        for i in range(1, ncol + 1):
            layer_mass = mass_in_layer[_idx2(i, k, ncol)]
            for m in range(1, gas_cnt + 1):
                wd_mmr[_idx3(i, kk, m, ncol, pver)] = trc_mass[_idx3(i, k, m, ncol, pver)] / layer_mass

    for m in range(1, gas_cnt + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx3(i, k, m, ncol, pver)
                dtwr[idx] = (wd_mmr[idx] - dtwr[idx]) / delt

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            if abs(lats[i - 1] * 180.0 / pi) > 60.0:
                if pmid[_idx2(i, k, pcols)] < 20000.0:
                    for m in range(1, gas_cnt + 1):
                        dtwr[_idx3(i, k, m, ncol, pver)] = 0.0

    for m in range(1, gas_cnt + 1):
        spc = mapping_to_mmr[m - 1]
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                wd_tend[_idx3(i, k, spc, pcols, pver)] += dtwr[_idx3(i, k, m, ncol, pver)]


@export
def neu_wetdep_henry_flags_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_cnt: int,
    nh3_ndx: int,
    co2_ndx: int,
    t0: float,
    ph: float,
    ph_inv: float,
    mapping_to_heff_p: cobj,
    dheff_p: cobj,
    tfld_p: cobj,
    heff_p: cobj,
    wrk_p: cobj,
    dk1s_p: cobj,
    dk2s_p: cobj,
    tckaqb_p: cobj,
):
    mapping_to_heff = Ptr[int](mapping_to_heff_p)
    dheff = Ptr[float](dheff_p)
    tfld = Ptr[float](tfld_p)
    heff = Ptr[float](heff_p)
    wrk = Ptr[float](wrk_p)
    dk1s = Ptr[float](dk1s_p)
    dk2s = Ptr[float](dk2s_p)
    tckaqb = Ptr[int](tckaqb_p)

    # Fortran declarations: tfld(pcols,pver), heff(ncol,pver,gas_wetdep_cnt).
    for m in range(1, gas_cnt + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                heff[_idx3(i, k, m, ncol, pver)] = 0.0

    for k in range(1, pver + 1):
        kk = pver - k + 1

        for i in range(1, ncol + 1):
            temp = tfld[_idx2(i, kk, pcols)]
            wrk[i - 1] = (t0 - temp) / (t0 * temp)

        for m in range(1, gas_cnt + 1):
            l = mapping_to_heff[m - 1]
            base = 6 * (l - 1)
            e298 = dheff[base]
            dhr = dheff[base + 1]

            for i in range(1, ncol + 1):
                heff[_idx3(i, k, m, ncol, pver)] = e298 * exp(dhr * wrk[i - 1])

            if dheff[base + 2] != 0.0 and dheff[base + 4] == 0.0:
                e298 = dheff[base + 2]
                dhr = dheff[base + 3]
                for i in range(1, ncol + 1):
                    dk1s[i - 1] = e298 * exp(dhr * wrk[i - 1])

                for i in range(1, ncol + 1):
                    idx = _idx3(i, k, m, ncol, pver)
                    if heff[idx] != 0.0:
                        heff[idx] = heff[idx] * (1.0 + dk1s[i - 1] * ph_inv)
                    else:
                        heff[idx] = dk1s[i - 1] * ph_inv

            if dheff[base + 4] != 0.0:
                if nh3_ndx > 0 or co2_ndx > 0:
                    e298 = dheff[base + 2]
                    dhr = dheff[base + 3]
                    for i in range(1, ncol + 1):
                        dk1s[i - 1] = e298 * exp(dhr * wrk[i - 1])

                    e298 = dheff[base + 4]
                    dhr = dheff[base + 5]
                    for i in range(1, ncol + 1):
                        dk2s[i - 1] = e298 * exp(dhr * wrk[i - 1])

                    if m == co2_ndx:
                        for i in range(1, ncol + 1):
                            idx = _idx3(i, k, m, ncol, pver)
                            heff[idx] = heff[idx] * (1.0 + dk1s[i - 1] * ph_inv) * (
                                1.0 + dk2s[i - 1] * ph_inv
                            )
                    elif m == nh3_ndx:
                        for i in range(1, ncol + 1):
                            idx = _idx3(i, k, m, ncol, pver)
                            heff[idx] = heff[idx] * (1.0 + dk1s[i - 1] * ph / dk2s[i - 1])

    for m in range(1, gas_cnt + 1):
        max_heff = heff[_idx3(1, 1, m, ncol, pver)]
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                val = heff[_idx3(i, k, m, ncol, pver)]
                if val > max_heff:
                    max_heff = val

        if max_heff > 1.0e4:
            tckaqb[m - 1] = 1
        else:
            tckaqb[m - 1] = 0


@inline
def _neu_wetdep_disgas_core(
    clwx: float,
    cfx: float,
    molmass: float,
    hstar: float,
    tm: float,
    pr: float,
    qm: float,
    qt: float,
) -> float:
    tmix = 258.0
    reteff = 0.5

    if tm >= 263.0:
        return (hstar * (qt / (qm * cfx)) * 0.029 * (pr / 1.0e3)) * (clwx * qm)
    elif tm <= tmix:
        muemp = exp(-14.2252 + (1.55704e-1 * tm) - (7.1929e-4 * (tm ** 2.0)))
        return muemp * (molmass / 18.0) * (clwx * qm)

    return reteff * ((hstar * (qt / (qm * cfx)) * 0.029 * (pr / 1.0e3)) * (clwx * qm))


@inline
def _neu_wetdep_raingas_core(
    rrain: float,
    dtscav: float,
    clwx: float,
    cfx: float,
    qm: float,
    qt: float,
    qtdis: float,
) -> float:
    qtdisstar = (qtdis * (qt * cfx)) / (qtdis + (qt * cfx))
    qtlf = (rrain * qtdisstar) / (clwx * qm * qt * cfx)
    return qt * cfx * (1.0 - exp(-dtscav * qtlf))


@inline
def _neu_wetdep_dempirical_core(cwater: float, rrate: float) -> float:
    rratex = rrate * 3600.0
    wx = cwater * 1.0e3

    if rratex > 0.04:
        theta = exp(-1.43 * log10(7.0 * rratex)) + 2.8
    else:
        theta = 5.0

    phi = rratex / (3600.0 * 10.0)
    eta = exp((3.01 * theta) - 10.5)
    beta = theta / (1.0 + 0.638)
    alpha = exp(4.0 * (beta - 3.5))
    bee = (0.638 * theta / (1.0 + 0.638)) - 1.0
    gamtheta = neu_wetdep_gamma_native_cb(theta)
    gambeta = neu_wetdep_gamma_native_cb(beta + 1.0)
    return (((wx * eta * gamtheta) / (1.0e6 * alpha * phi * gambeta)) ** (-1.0 / bee)) * 10.0


@export
def neu_wetdep_disgas_codon(
    clwx: float,
    cfx: float,
    molmass: float,
    hstar: float,
    tm: float,
    pr: float,
    qm: float,
    qt: float,
    qtdis_p: cobj,
):
    qtdis = Ptr[float](qtdis_p)
    qtdis[0] = _neu_wetdep_disgas_core(clwx, cfx, molmass, hstar, tm, pr, qm, qt)


@export
def neu_wetdep_raingas_codon(
    rrain: float,
    dtscav: float,
    clwx: float,
    cfx: float,
    qm: float,
    qt: float,
    qtdis: float,
    qtrain_p: cobj,
):
    qtrain = Ptr[float](qtrain_p)
    qtrain[0] = _neu_wetdep_raingas_core(rrain, dtscav, clwx, cfx, qm, qt, qtdis)


@export
def neu_wetdep_washo_codon(
    lpar: int,
    ntrace: int,
    hno3_ndx: int,
    do_diag: int,
    dempirical_impl: int,
    dtscav: float,
    garea: float,
    adj_factor: float,
    qttjfl_p: cobj,
    qm_p: cobj,
    pofl_p: cobj,
    delz_p: cobj,
    rls_p: cobj,
    clwc_p: cobj,
    ciwc_p: cobj,
    cfr_p: cobj,
    tem_p: cobj,
    evaprate_p: cobj,
    hstar_p: cobj,
    tcmass_p: cobj,
    tckaqb_p: cobj,
    tcnion_p: cobj,
    qt_rain_p: cobj,
    qt_rime_p: cobj,
    qt_wash_p: cobj,
    qt_evap_p: cobj,
    cfxx_p: cobj,
    qtt_p: cobj,
    qttnew_p: cobj,
):
    qttjfl = Ptr[float](qttjfl_p)
    qm = Ptr[float](qm_p)
    pofl = Ptr[float](pofl_p)
    delz = Ptr[float](delz_p)
    rls = Ptr[float](rls_p)
    clwc = Ptr[float](clwc_p)
    ciwc = Ptr[float](ciwc_p)
    cfr = Ptr[float](cfr_p)
    tem = Ptr[float](tem_p)
    evaprate = Ptr[float](evaprate_p)
    hstar = Ptr[float](hstar_p)
    tcmass = Ptr[float](tcmass_p)
    tckaqb = Ptr[int](tckaqb_p)
    tcnion = Ptr[int](tcnion_p)
    qt_rain = Ptr[float](qt_rain_p)
    qt_rime = Ptr[float](qt_rime_p)
    qt_wash = Ptr[float](qt_wash_p)
    qt_evap = Ptr[float](qt_evap_p)
    cfxx = Ptr[float](cfxx_p)
    qtt = Ptr[float](qtt_p)
    qttnew = Ptr[float](qttnew_p)

    zero = 0.0
    one = 1.0
    cfmin = 0.1
    cwmin = 1.0e-5
    dmin = 1.0e-1
    volpow = 1.0 / 3.0
    rhorain = 1.0e3
    rhosnowfix = 1.0e2
    coleffrain = 0.7
    tmix = 258.0
    tfroz = 240.0
    coleffaer = 0.05
    tice = 263.0
    four = 4.0

    le = lpar - 1
    n = 1
    while n <= ntrace:
        ll = 1
        while ll <= lpar:
            ln_idx = _idx2(ll, n, lpar)
            qtt[ll - 1] = qttjfl[ln_idx]
            qttnew[ll - 1] = qttjfl[ln_idx]
            ll += 1

        is_hno3 = 0
        if n == hno3_ndx:
            is_hno3 = 1
            ll = 1
            while ll <= lpar:
                qt_rain[ll - 1] = zero
                qt_rime[ll - 1] = zero
                qt_wash[ll - 1] = zero
                qt_evap[ll - 1] = zero
                ll += 1

        if tckaqb[n - 1] != 0:
            lwashtyp = 1
        else:
            lwashtyp = 2

        if tcnion[n - 1] != 0:
            licetyp = 1
        else:
            licetyp = 2

        qttopaa = zero
        qttopca = zero
        rca = zero
        fca = zero
        dca = zero
        rama = zero
        fama = zero
        dama = zero
        ampct = zero
        amclpct = zero
        clnewpct = zero
        clnewampct = zero
        cloldpct = zero
        cloldampct = zero

        if le >= 1:
            if rls[le - 1] > zero:
                cfxx[le - 1] = max(cfmin, cfr[le - 1])
            else:
                cfxx[le - 1] = cfr[le - 1]

        skip_species = 0
        l = le
        while l >= 1:
            lm1 = l - 1
            ln_idx = _idx2(l, n, lpar)
            hstar_ln = hstar[ln_idx]

            fax = zero
            rax = zero
            dax = zero
            fcxa = zero
            fcxb = zero
            dcxa = zero
            dcxb = zero
            rcxa = zero
            rcxb = zero
            qtdiscf = zero
            qtdisrime = zero
            qtdiscxa = zero
            qtevapaxp = zero
            qtevapaxw = zero
            qtevapax = zero
            qtwashax = zero
            qtevapcxap = zero
            qtevapcxaw = zero
            qtevapcxa = zero
            qtrimecxa = zero
            qtwashcxa = zero
            qtraincxa = zero
            qtraincxb = zero
            rampct = zero
            rprecip = zero
            deltarimemass = zero
            deltarime = zero
            dor = zero
            dnew = zero
            qttopaax = zero
            qttopcax = zero
            freezing_l = 0
            if tem[l - 1] < tice:
                freezing_l = 1

            if rls[l - 1] > zero:
                fax = max(zero, fama * (one - evaprate[l - 1]))
                rax = rama
                if fama > zero:
                    if freezing_l != 0:
                        dax = dama
                    else:
                        dax = four
                else:
                    dax = zero

                if rama > zero:
                    qtevapaxp = min(qttopaa, evaprate[l - 1] * qttopaa)
                else:
                    qtevapaxp = zero

                wrk = rax * fax + rca * fca
                if wrk > 0.0:
                    rnew_tst = rls[l - 1] / (garea * wrk)
                else:
                    rnew_tst = 10.0
                rnew = (rls[l - 1] / garea) - (rax * fax + rca * fca)

                if (rls[l - 1] / garea) > adj_factor * (rax * fax + rca * fca):
                    if cfxx[l - 1] == zero:
                        ll = 1
                        while ll <= lpar:
                            qttjfl[_idx2(ll, n, lpar)] = qtt[ll - 1]
                            ll += 1
                        skip_species = 1
                        break

                    clwx = max(clwc[l - 1] + ciwc[l - 1], cwmin * cfxx[l - 1])
                    fcxa = fca
                    fcxb = max(zero, cfxx[l - 1] - fcxa)

                    if freezing_l != 0:
                        coleffsnow = exp(2.5e-2 * (tem[l - 1] - tice))
                        if tem[l - 1] <= tfroz:
                            rhosnow = rhosnowfix
                        else:
                            rhosnow = 0.303 * (tem[l - 1] - tfroz) * rhosnowfix

                        if fcxa > zero:
                            if dca > zero:
                                deltarimemass = clwx * qm[l - 1] * (fcxa / cfxx[l - 1]) * (
                                    one
                                    - exp(
                                        (-coleffsnow / (dca * 1.0e-3))
                                        * ((rca) / (2.0 * rhosnow))
                                        * dtscav
                                    )
                                )
                            else:
                                deltarimemass = zero
                        else:
                            deltarimemass = zero

                        if fcxa > zero:
                            deltarime = min(rnew / fcxa, deltarimemass / (fcxa * garea * dtscav))
                        else:
                            deltarime = zero

                        if rca > zero:
                            dor = max(dmin, (((rca + deltarime) / rca) ** volpow) * dca)
                        else:
                            dor = zero

                        rprecip = (rnew - (deltarime * fcxa)) / cfxx[l - 1]
                        rcxa = rca + deltarime + rprecip
                        rcxb = rprecip

                        if rprecip > zero:
                            wemp = (clwx * qm[l - 1]) / (garea * cfxx[l - 1] * delz[l - 1])
                            remp = rprecip / (rhorain / 1.0e3)
                            if dempirical_impl == 0:
                                dnew = neu_wetdep_dempirical_native_cb(wemp, remp)
                            else:
                                dnew = _neu_wetdep_dempirical_core(wemp, remp)
                            dnew = max(dmin, dnew)
                            if fcxb > zero:
                                dcxb = dnew
                            else:
                                dcxb = zero
                        else:
                            dcxb = zero

                        if fcxa > zero:
                            wemp = (clwx * qm[l - 1] * (fcxa / cfxx[l - 1])) / (garea * fcxa * delz[l - 1])
                            remp = rcxa / (rhorain / 1.0e3)
                            if dempirical_impl == 0:
                                demp = neu_wetdep_dempirical_native_cb(wemp, remp)
                            else:
                                demp = _neu_wetdep_dempirical_core(wemp, remp)
                            dcxa = ((rca + deltarime) / rcxa) * dor + (rprecip / rcxa) * dnew
                            dcxa = max(demp, dcxa)
                            dcxa = max(dmin, dcxa)
                        else:
                            dcxa = zero

                        if qtt[l - 1] > zero:
                            if rprecip > zero:
                                if licetyp == 1:
                                    rrain = rprecip * garea
                                    qtdiscf = _neu_wetdep_disgas_core(
                                        clwx,
                                        cfxx[l - 1],
                                        tcmass[n - 1],
                                        hstar_ln,
                                        tem[l - 1],
                                        pofl[l - 1],
                                        qm[l - 1],
                                        qtt[l - 1] * cfxx[l - 1],
                                    )
                                    qtrain = _neu_wetdep_raingas_core(
                                        rrain,
                                        dtscav,
                                        clwx,
                                        cfxx[l - 1],
                                        qm[l - 1],
                                        qtt[l - 1],
                                        qtdiscf,
                                    )
                                    wrk = qtrain / cfxx[l - 1]
                                    qtraincxa = fcxa * wrk
                                    qtraincxb = fcxb * wrk
                                else:
                                    qtraincxa = zero
                                    qtraincxb = zero

                            if deltarime > zero:
                                if licetyp == 1:
                                    if tem[l - 1] <= tfroz:
                                        rhosnow = rhosnowfix
                                    else:
                                        rhosnow = 0.303 * (tem[l - 1] - tfroz) * rhosnowfix
                                    qtcxa = qtt[l - 1] * fcxa
                                    qtdisrime = _neu_wetdep_disgas_core(
                                        clwx * (fcxa / cfxx[l - 1]),
                                        fcxa,
                                        tcmass[n - 1],
                                        hstar_ln,
                                        tem[l - 1],
                                        pofl[l - 1],
                                        qm[l - 1],
                                        qtcxa,
                                    )
                                    qtdisstar = (qtdisrime * qtcxa) / (qtdisrime + qtcxa)
                                    qtrimecxa = qtcxa * (
                                        one
                                        - exp(
                                            (-coleffsnow / (dca * 1.0e-3))
                                            * (rca / (2.0 * rhosnow))
                                            * (qtdisstar / qtcxa)
                                            * dtscav
                                        )
                                    )
                                    qtrimecxa = min(
                                        qtrimecxa,
                                        ((rnew * garea * dtscav) / (clwx * qm[l - 1] * (fcxa / cfxx[l - 1])))
                                        * qtdisstar,
                                    )
                                else:
                                    qtrimecxa = zero
                        else:
                            qtraincxa = zero
                            qtraincxb = zero
                            qtrimecxa = zero

                        qtwashcxa = zero
                        qtevapcxa = zero
                    else:
                        if fcxa > zero:
                            deltarimemass = (clwx * qm[l - 1]) * (fcxa / cfxx[l - 1]) * (
                                one - exp(-0.24 * coleffrain * ((rca) ** 0.75) * dtscav)
                            )
                        else:
                            deltarimemass = zero

                        if fcxa > zero:
                            deltarime = min(rnew / fcxa, deltarimemass / (fcxa * garea * dtscav))
                        else:
                            deltarime = zero

                        rprecip = (rnew - (deltarime * fcxa)) / cfxx[l - 1]
                        rcxa = rca + deltarime + rprecip
                        rcxb = rprecip
                        dcxa = four
                        if fcxb > zero:
                            dcxb = four
                        else:
                            dcxb = zero

                        if qtt[l - 1] > zero:
                            if rprecip > zero:
                                rrain = rprecip * garea
                                qtdiscf = _neu_wetdep_disgas_core(
                                    clwx,
                                    cfxx[l - 1],
                                    tcmass[n - 1],
                                    hstar_ln,
                                    tem[l - 1],
                                    pofl[l - 1],
                                    qm[l - 1],
                                    qtt[l - 1] * cfxx[l - 1],
                                )
                                qtrain = _neu_wetdep_raingas_core(
                                    rrain,
                                    dtscav,
                                    clwx,
                                    cfxx[l - 1],
                                    qm[l - 1],
                                    qtt[l - 1],
                                    qtdiscf,
                                )
                                wrk = qtrain / cfxx[l - 1]
                                qtraincxa = fcxa * wrk
                                qtraincxb = fcxb * wrk

                            if deltarime > zero:
                                qtcxa = qtt[l - 1] * fcxa
                                qtdisrime = _neu_wetdep_disgas_core(
                                    clwx * (fcxa / cfxx[l - 1]),
                                    fcxa,
                                    tcmass[n - 1],
                                    hstar_ln,
                                    tem[l - 1],
                                    pofl[l - 1],
                                    qm[l - 1],
                                    qtcxa,
                                )
                                qtdisstar = (qtdisrime * qtcxa) / (qtdisrime + qtcxa)
                                qtrimecxa = qtcxa * (
                                    one
                                    - exp(
                                        -0.24 * coleffrain * ((rca) ** 0.75) * (qtdisstar / qtcxa) * dtscav
                                    )
                                )
                                qtrimecxa = min(
                                    qtrimecxa,
                                    ((rnew * garea * dtscav) / (clwx * qm[l - 1] * (fcxa / cfxx[l - 1])))
                                    * qtdisstar,
                                )
                            else:
                                qtrimecxa = zero
                        else:
                            qtraincxa = zero
                            qtraincxb = zero
                            qtrimecxa = zero

                        if rca > zero:
                            qtprecip = fcxa * qtt[l - 1] - qtdisrime
                            if lwashtyp == 1:
                                if qtprecip > zero:
                                    qtwashcxa = qtprecip * (
                                        one - exp(-0.24 * coleffaer * ((rca) ** 0.75) * dtscav)
                                    )
                                else:
                                    qtwashcxa = zero
                                qtevapcxa = zero
                            else:
                                rwash = rca * garea
                                if qtprecip > zero:
                                    if fca == zero:
                                        qtwashcxa = zero
                                        qtevapcxa = zero
                                    else:
                                        fwash = (rwash * hstar_ln * 29.0e-6 * pofl[l - 1]) / (qm[l - 1] * fca)
                                        qtmax = qtprecip * fwash * dtscav
                                        if qtmax > (qttopca + qtrimecxa):
                                            qtdif = min(qtprecip, qtmax - (qttopca + qtrimecxa))
                                            qtwashcxa = qtdif * (one - exp(-dtscav * fwash))
                                            qtevapcxa = zero
                                        else:
                                            qtwashcxa = zero
                                            qtevapcxa = (qttopca + qtrimecxa) - qtmax
                                else:
                                    qtwashcxa = zero
                                    qtevapcxa = zero
                else:
                    clwx = clwc[l - 1] + ciwc[l - 1]
                    fcxa = fca
                    fcxb = max(zero, cfxx[l - 1] - fcxa)
                    rcxb = zero
                    dcxb = zero
                    qtraincxa = zero
                    qtraincxb = zero
                    qtrimecxa = zero

                    if fcxa > zero:
                        rcxa = min(rca, rls[l - 1] / (garea * fcxa))
                        if fax > zero and ((rcxa + 1.0e-12) < rls[l - 1] / (garea * fcxa)):
                            raxadjf = rls[l - 1] / garea - rcxa * fcxa
                            rampct = raxadjf / (rax * fax)
                            faxadj = rampct * fax
                            if faxadj > zero:
                                raxadj = raxadjf / faxadj
                            else:
                                raxadj = zero
                        else:
                            raxadj = zero
                            rampct = zero
                            faxadj = zero
                    else:
                        rcxa = zero
                        if fax > zero:
                            raxadjf = rls[l - 1] / garea
                            rampct = raxadjf / (rax * fax)
                            faxadj = rampct * fax
                            if faxadj > zero:
                                raxadj = raxadjf / faxadj
                            else:
                                raxadj = zero
                        else:
                            raxadj = zero
                            rampct = zero
                            faxadj = zero

                    qtevapaxp = min(qttopaa, qttopaa - (rampct * (qttopaa - qtevapaxp)))
                    fax = faxadj
                    rax = raxadj

                    if rcxa <= zero:
                        qtevapcxa = qttopca
                        rcxa = zero
                        dcxa = zero
                    else:
                        if freezing_l != 0:
                            qtwashcxa = zero
                            dcxa = ((rcxa / rca) ** volpow) * dca
                            if licetyp == 1:
                                if tem[l - 1] <= tmix:
                                    massloss = (rca - rcxa) * fcxa * garea * dtscav
                                    qtevapcxa = _neu_wetdep_disgas_core(
                                        massloss / qm[l - 1],
                                        fcxa,
                                        tcmass[n - 1],
                                        hstar_ln,
                                        tem[l - 1],
                                        pofl[l - 1],
                                        qm[l - 1],
                                        qtt[l - 1],
                                    )
                                    qtevapcxa = min(qttopca, qtevapcxa)
                                else:
                                    qtevapcxa = zero
                            else:
                                qtevapcxa = zero
                        else:
                            qtevapcxap = (rca - rcxa) / rca * qttopca
                            dcxa = four
                            qtcxa = fcxa * qtt[l - 1]
                            if lwashtyp == 1:
                                if qtt[l - 1] > zero:
                                    qtdiscxa = _neu_wetdep_disgas_core(
                                        clwx * (fcxa / cfxx[l - 1]),
                                        fcxa,
                                        tcmass[n - 1],
                                        hstar_ln,
                                        tem[l - 1],
                                        pofl[l - 1],
                                        qm[l - 1],
                                        qtcxa,
                                    )
                                    if qtcxa > qtdiscxa:
                                        qtwashcxa = (qtcxa - qtdiscxa) * (
                                            one - exp(-0.24 * coleffaer * ((rcxa) ** 0.75) * dtscav)
                                        )
                                    else:
                                        qtwashcxa = zero
                                    qtevapcxaw = zero
                                else:
                                    qtwashcxa = zero
                                    qtevapcxaw = zero
                            else:
                                rwash = rcxa * garea
                                if fcxa == zero:
                                    qtwashcxa = zero
                                    qtevapcxaw = zero
                                else:
                                    fwash = (rwash * hstar_ln * 29.0e-6 * pofl[l - 1]) / (qm[l - 1] * fcxa)
                                    qtmax = (qtcxa - qtdiscxa) * fwash * dtscav
                                    if qtmax > qttopca:
                                        qtdif = min(qtcxa - qtdiscxa, qtmax - qttopca)
                                        qtwashcxa = qtdif * (one - exp(-dtscav * fwash))
                                        qtevapcxaw = zero
                                    else:
                                        qtwashcxa = zero
                                        qtevapcxaw = qttopca - qtmax
                            qtevapcxa = qtevapcxap + qtevapcxaw
            else:
                rnew = zero
                qtevapcxa = qttopca
                qtevapax = qttopaa
                if l > 1:
                    if rls[lm1 - 1] > zero:
                        cfxx[lm1 - 1] = max(cfmin, cfr[lm1 - 1])
                    else:
                        cfxx[lm1 - 1] = cfr[lm1 - 1]
                ampct = zero
                amclpct = zero
                clnewpct = zero
                clnewampct = zero
                cloldpct = zero
                cloldampct = zero
                rca = zero
                rama = zero
                fca = zero
                fama = zero
                dca = zero
                dama = zero

            if skip_species == 0:
                if rls[l - 1] > zero:
                    if rax > zero:
                        if freezing_l == 0:
                            qtax = fax * qtt[l - 1]
                            if lwashtyp == 1:
                                qtwashax = qtax * (one - exp(-0.24 * coleffaer * ((rax) ** 0.75) * dtscav))
                                qtevapaxw = zero
                            else:
                                rwash = rax * garea
                                if fax == zero:
                                    qtwashax = zero
                                    qtevapaxw = zero
                                else:
                                    fwash = (rwash * hstar_ln * 29.0e-6 * pofl[l - 1]) / (qm[l - 1] * fax)
                                    qtmax = qtax * fwash * dtscav
                                    if qtmax > qttopaa:
                                        qtdif = min(qtax, qtmax - qttopaa)
                                        qtwashax = qtdif * (one - exp(-dtscav * fwash))
                                        qtevapaxw = zero
                                    else:
                                        qtwashax = zero
                                        qtevapaxw = qttopaa - qtmax
                        else:
                            qtevapaxw = zero
                            qtwashax = zero
                    else:
                        qtevapaxw = zero
                        qtwashax = zero
                    qtevapax = qtevapaxp + qtevapaxw

                    if l > 1:
                        fama = max(fcxa + fcxb + fax - cfr[lm1 - 1], zero)

                        if cfr[lm1 - 1] >= cfmin:
                            cfxx[lm1 - 1] = cfr[lm1 - 1]
                        else:
                            if adj_factor * (rls[lm1 - 1] / garea) >= (
                                (rcxa * fcxa + rcxb * fcxb + rax * fax) * (one - evaprate[lm1 - 1])
                            ):
                                cfxx[lm1 - 1] = cfmin
                            else:
                                cfxx[lm1 - 1] = cfr[lm1 - 1]

                        if fax > zero:
                            ampct = max(zero, min(one, (cfxx[l - 1] + fax - cfxx[lm1 - 1]) / fax))
                            amclpct = one - ampct
                        else:
                            ampct = zero
                            amclpct = zero

                        if fcxb > zero:
                            clnewpct = max(zero, min((cfxx[lm1 - 1] - fcxa) / fcxb, one))
                            clnewampct = one - clnewpct
                        else:
                            clnewpct = zero
                            clnewampct = zero

                        if fcxa > zero:
                            cloldpct = max(zero, min(cfxx[lm1 - 1] / fcxa, one))
                            cloldampct = one - cloldpct
                        else:
                            cloldpct = zero
                            cloldampct = zero

                        fca = min(cfxx[lm1 - 1], fcxa * cloldpct + clnewpct * fcxb + amclpct * fax)
                        if fca > zero:
                            rca = (rcxa * fcxa * cloldpct + rcxb * fcxb * clnewpct + rax * fax * amclpct) / fca
                            if rca > zero:
                                dca = (rcxa * fcxa * cloldpct) / (rca * fca) * dcxa + (
                                    rcxb * fcxb * clnewpct
                                ) / (rca * fca) * dcxb + (rax * fax * amclpct) / (rca * fca) * dax
                            else:
                                dca = zero
                                fca = zero
                        else:
                            fca = zero
                            dca = zero
                            rca = zero

                        fama = fcxa + fcxb + fax - cfxx[lm1 - 1]
                        if fama > zero:
                            rama = (
                                rcxa * fcxa * cloldampct + rcxb * fcxb * clnewampct + rax * fax * ampct
                            ) / fama
                            if rama > zero:
                                dama = (rcxa * fcxa * cloldampct) / (rama * fama) * dcxa + (
                                    rcxb * fcxb * clnewampct
                                ) / (rama * fama) * dcxb + (rax * fax * ampct) / (rama * fama) * dax
                            else:
                                fama = zero
                                dama = zero
                        else:
                            fama = zero
                            dama = zero
                            rama = zero
                    else:
                        ampct = zero
                        amclpct = zero
                        clnewpct = zero
                        clnewampct = zero
                        cloldpct = zero
                        cloldampct = zero

                qtnetlcxa = qtraincxa + qtrimecxa + qtwashcxa - qtevapcxa
                qtnetlcxa = min(qtt[l - 1] * fcxa, qtnetlcxa)
                qtnetlcxb = qtraincxb
                qtnetlcxb = min(qtt[l - 1] * fcxb, qtnetlcxb)
                qtnetlax = qtwashax - qtevapax
                qtnetlax = min(qtt[l - 1] * fax, qtnetlax)
                qttnew[l - 1] = qtt[l - 1] - (qtnetlcxa + qtnetlcxb + qtnetlax)

                if do_diag != 0 and is_hno3 != 0:
                    qt_rain[l - 1] = qtraincxa + qtraincxb
                    qt_rime[l - 1] = qtrimecxa
                    qt_wash[l - 1] = qtwashcxa + qtwashax
                    qt_evap[l - 1] = qtevapcxa + qtevapax

                qttopcax = (qttopca + qtnetlcxa) * cloldpct + qtnetlcxb * clnewpct + (qttopaa + qtnetlax) * amclpct
                qttopaax = (qttopca + qtnetlcxa) * cloldampct + qtnetlcxb * clnewampct + (qttopaa + qtnetlax) * ampct
                qttopca = qttopcax
                qttopaa = qttopaax

            l -= 1

        if skip_species == 0:
            ll = 1
            while ll <= le:
                qttjfl[_idx2(ll, n, lpar)] = qttnew[ll - 1]
                ll += 1

        n += 1


@export
def setsox_init_fields_codon(
    stage: int,
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    nfs: int,
    cloud_borne_flag: int,
    inv_so2_flag: int,
    inv_h2o2_flag: int,
    inv_o3_flag: int,
    inv_ho2_flag: int,
    id_so2: int,
    id_hno3: int,
    id_h2o2: int,
    id_nh3: int,
    id_o3: int,
    id_ho2: int,
    id_h2so4: int,
    id_so4: int,
    id_msa: int,
    ph0: float,
    xhnm_p: cobj,
    invariants_p: cobj,
    qin_p: cobj,
    cfact_p: cobj,
    xph_p: cobj,
    xso2_p: cobj,
    xhno3_p: cobj,
    xh2o2_p: cobj,
    xnh3_p: cobj,
    xo3_p: cobj,
    xho2_p: cobj,
    xh2so4_p: cobj,
    xso4_p: cobj,
    xno3_p: cobj,
    xnh4_p: cobj,
    xmsa_p: cobj,
):
    xhnm = Ptr[float](xhnm_p)
    cfact = Ptr[float](cfact_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, ncol)
            cfact[idx] = xhnm[idx] * 1.0e6 * 1.38e-23 / 287.0 * 1.0e-3

    if stage == 1:
        return

    invariants = Ptr[float](invariants_p)
    qin = Ptr[float](qin_p)
    xph = Ptr[float](xph_p)
    xso2 = Ptr[float](xso2_p)
    xhno3 = Ptr[float](xhno3_p)
    xh2o2 = Ptr[float](xh2o2_p)
    xnh3 = Ptr[float](xnh3_p)
    xo3 = Ptr[float](xo3_p)
    xho2 = Ptr[float](xho2_p)
    xh2so4 = Ptr[float](xh2so4_p)
    xso4 = Ptr[float](xso4_p)
    xno3 = Ptr[float](xno3_p)
    xnh4 = Ptr[float](xnh4_p)
    xmsa = Ptr[float](xmsa_p)
    xph0 = 10.0 ** (-ph0)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, ncol)
            xso4[idx] = 0.0
            xno3[idx] = 0.0
            xnh4[idx] = 0.0
            xph[idx] = xph0

            if inv_so2_flag != 0:
                xso2[idx] = invariants[_idx3(i, k, id_so2, ncol, pver)] / xhnm[idx]
            else:
                xso2[idx] = qin[_idx3(i, k, id_so2, ncol, pver)]

            if id_hno3 > 0:
                xhno3[idx] = qin[_idx3(i, k, id_hno3, ncol, pver)]
            else:
                xhno3[idx] = 0.0

            if inv_h2o2_flag != 0:
                xh2o2[idx] = invariants[_idx3(i, k, id_h2o2, ncol, pver)] / xhnm[idx]
            else:
                xh2o2[idx] = qin[_idx3(i, k, id_h2o2, ncol, pver)]

            if id_nh3 > 0:
                xnh3[idx] = qin[_idx3(i, k, id_nh3, ncol, pver)]
            else:
                xnh3[idx] = 0.0

            if inv_o3_flag != 0:
                xo3[idx] = invariants[_idx3(i, k, id_o3, ncol, pver)] / xhnm[idx]
            else:
                xo3[idx] = qin[_idx3(i, k, id_o3, ncol, pver)]

            if inv_ho2_flag != 0:
                xho2[idx] = invariants[_idx3(i, k, id_ho2, ncol, pver)] / xhnm[idx]
            else:
                xho2[idx] = qin[_idx3(i, k, id_ho2, ncol, pver)]

            if cloud_borne_flag != 0:
                xh2so4[idx] = qin[_idx3(i, k, id_h2so4, ncol, pver)]
            else:
                xso4[idx] = qin[_idx3(i, k, id_so4, ncol, pver)]

            if id_msa > 0:
                xmsa[idx] = qin[_idx3(i, k, id_msa, ncol, pver)]


@export
def setsox_ph_solve_codon(
    ncol: int,
    pcols: int,
    pver: int,
    itermax: int,
    cloud_borne_flag: int,
    const0: float,
    ra: float,
    xkw: float,
    so4_fact: float,
    press_p: cobj,
    tfld_p: cobj,
    cldfrc_p: cobj,
    xhnm_p: cobj,
    xlwc_p: cobj,
    xso4c_p: cobj,
    xnh4c_p: cobj,
    xno3c_p: cobj,
    xso4_p: cobj,
    xnh4_p: cobj,
    xno3_p: cobj,
    xso2_p: cobj,
    xhno3_p: cobj,
    xnh3_p: cobj,
    xph_p: cobj,
):
    press = Ptr[float](press_p)
    tfld = Ptr[float](tfld_p)
    cldfrc = Ptr[float](cldfrc_p)
    xhnm = Ptr[float](xhnm_p)
    xlwc = Ptr[float](xlwc_p)
    xso4c = Ptr[float](xso4c_p)
    xnh4c = Ptr[float](xnh4c_p)
    xno3c = Ptr[float](xno3c_p)
    xso4 = Ptr[float](xso4_p)
    xnh4 = Ptr[float](xnh4_p)
    xno3 = Ptr[float](xno3_p)
    xso2 = Ptr[float](xso2_p)
    xhno3 = Ptr[float](xhno3_p)
    xnh3 = Ptr[float](xnh3_p)
    xph = Ptr[float](xph_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, ncol)
            idxp = _idx2(i, k, pcols)
            if cloud_borne_flag != 0 and cldfrc[idxp] > 0.0:
                xso4[idx] = xso4c[idxp] / cldfrc[idxp]
                xnh4[idx] = xnh4c[idxp] / cldfrc[idxp]
                xno3[idx] = xno3c[idxp] / cldfrc[idxp]

            xl = xlwc[idxp]
            if xl >= 1.0e-8:
                work1 = 1.0 / tfld[idxp] - 1.0 / 298.0
                pz = 0.01 * press[idxp]
                tz = tfld[idxp]
                patm = pz / 1013.0

                xk = 2.1e5 * exp(8700.0 * work1)
                xe = 15.4
                fact1_hno3 = xk * xe * patm * xhno3[idx]
                fact2_hno3 = xk * ra * tz * xl
                fact3_hno3 = xe

                xk = 1.23 * exp(3120.0 * work1)
                xe = 1.7e-2 * exp(2090.0 * work1)
                x2 = 6.0e-8 * exp(1120.0 * work1)
                fact1_so2 = xk * xe * patm * xso2[idx]
                fact2_so2 = xk * ra * tz * xl
                fact3_so2 = xe
                fact4_so2 = x2

                xk = 58.0 * exp(4085.0 * work1)
                xe = 1.7e-5 * exp(-4325.0 * work1)
                fact1_nh3 = (xk * xe * patm / xkw) * (xnh3[idx] + xnh4[idx])
                fact2_nh3 = xk * ra * tz * xl
                fact3_nh3 = xe / xkw

                eh2o = xkw
                co2g = 330.0e-6
                xk = 3.1e-2 * exp(2423.0 * work1)
                xe = 4.3e-7 * exp(-913.0 * work1)
                eco2 = xk * xe * co2g * patm
                eso4 = xso4[idx] * xhnm[idx] * const0 / xl

                converged = 0
                yph_lo = 0.0
                yph_hi = 0.0
                ynetpos_lo = 0.0
                ynetpos_hi = 0.0
                for iter in range(1, itermax + 1):
                    if iter == 1:
                        yph_lo = 2.0
                        yph_hi = yph_lo
                        yph = yph_lo
                    elif iter == 2:
                        yph_hi = 7.0
                        yph = yph_hi
                    else:
                        yph = 0.5 * (yph_lo + yph_hi)

                    xph[idx] = 10.0 ** (-yph)
                    ehno3 = fact1_hno3 / (1.0 + fact2_hno3 * (1.0 + fact3_hno3 / xph[idx]))
                    eso2 = fact1_so2 / (
                        1.0
                        + fact2_so2
                        * (1.0 + (fact3_so2 / xph[idx]) * (1.0 + fact4_so2 / xph[idx]))
                    )
                    enh3 = fact1_nh3 / (1.0 + fact2_nh3 * (1.0 + fact3_nh3 * xph[idx]))

                    tmp_nh4 = enh3 * xph[idx]
                    tmp_hso3 = eso2 / xph[idx]
                    tmp_so3 = tmp_hso3 * 2.0 * fact4_so2 / xph[idx]
                    tmp_hco3 = eco2 / xph[idx]
                    tmp_oh = eh2o / xph[idx]
                    tmp_no3 = ehno3 / xph[idx]
                    tmp_so4 = so4_fact * eso4
                    tmp_pos = xph[idx] + tmp_nh4
                    tmp_neg = tmp_oh + tmp_hco3 + tmp_no3 + tmp_hso3 + tmp_so3 + tmp_so4
                    ynetpos = tmp_pos - tmp_neg

                    if iter > 2:
                        if ynetpos == 0.0:
                            converged = 1
                            break
                        elif ynetpos >= 0.0:
                            yph_lo = yph
                            ynetpos_lo = ynetpos
                        else:
                            yph_hi = yph
                            ynetpos_hi = ynetpos

                        if abs(yph_hi - yph_lo) <= 0.005:
                            yph = 0.5 * (yph_hi + yph_lo)
                            xph[idx] = 10.0 ** (-yph)
                            converged = 1
                            break
                    elif iter == 1:
                        if ynetpos <= 0.0:
                            converged = 1
                            break
                        ynetpos_lo = ynetpos
                    else:
                        if ynetpos >= 0.0:
                            converged = 1
                            break
                        ynetpos_hi = ynetpos
            else:
                xph[idx] = 1.0e-7


@export
def setsox_aqchem_predict_codon(
    ncol: int,
    pcols: int,
    pver: int,
    cloud_borne_flag: int,
    modal_aerosols_flag: int,
    id_nh3: int,
    dtime: float,
    const0: float,
    kh0: float,
    kh1: float,
    kh2: float,
    kh3: float,
    ra: float,
    xkw: float,
    press_p: cobj,
    tfld_p: cobj,
    xhnm_p: cobj,
    xlwc_p: cobj,
    xph_p: cobj,
    xho2_p: cobj,
    xhno3_p: cobj,
    xno3_p: cobj,
    xh2o2_p: cobj,
    xso2_p: cobj,
    xo3_p: cobj,
    xnh3_p: cobj,
    xnh4_p: cobj,
    xso4_p: cobj,
    xso4_init_p: cobj,
    xdelso4hp_p: cobj,
    hno3g_p: cobj,
    nh3g_p: cobj,
    hehno3_p: cobj,
    heh2o2_p: cobj,
    heso2_p: cobj,
    henh3_p: cobj,
    heo3_p: cobj,
):
    press = Ptr[float](press_p)
    tfld = Ptr[float](tfld_p)
    xhnm = Ptr[float](xhnm_p)
    xlwc = Ptr[float](xlwc_p)
    xph = Ptr[float](xph_p)
    xho2 = Ptr[float](xho2_p)
    xhno3 = Ptr[float](xhno3_p)
    xno3 = Ptr[float](xno3_p)
    xh2o2 = Ptr[float](xh2o2_p)
    xso2 = Ptr[float](xso2_p)
    xo3 = Ptr[float](xo3_p)
    xnh3 = Ptr[float](xnh3_p)
    xnh4 = Ptr[float](xnh4_p)
    xso4 = Ptr[float](xso4_p)
    xso4_init = Ptr[float](xso4_init_p)
    xdelso4hp = Ptr[float](xdelso4hp_p)
    hno3g = Ptr[float](hno3g_p)
    nh3g = Ptr[float](nh3g_p)
    hehno3 = Ptr[float](hehno3_p)
    heh2o2 = Ptr[float](heh2o2_p)
    heso2 = Ptr[float](heso2_p)
    henh3 = Ptr[float](henh3_p)
    heo3 = Ptr[float](heo3_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, ncol)
            idxp = _idx2(i, k, pcols)
            work1 = 1.0 / tfld[idxp] - 1.0 / 298.0
            tz = tfld[idxp]
            xl = xlwc[idxp]
            patm = press[idxp] / 101300.0
            xam = press[idxp] / (1.38e-23 * tz)

            xk = 2.1e5 * exp(8700.0 * work1)
            xe = 15.4
            hehno3[idx] = xk * (1.0 + xe / xph[idx])

            xk = 7.4e4 * exp(6621.0 * work1)
            xe = 2.2e-12 * exp(-3730.0 * work1)
            heh2o2[idx] = xk * (1.0 + xe / xph[idx])

            xk = 1.23 * exp(3120.0 * work1)
            xe = 1.7e-2 * exp(2090.0 * work1)
            x2 = 6.0e-8 * exp(1120.0 * work1)
            wrk = xe / xph[idx]
            heso2[idx] = xk * (1.0 + wrk * (1.0 + x2 / xph[idx]))

            xk = 58.0 * exp(4085.0 * work1)
            xe = 1.7e-5 * exp(-4325.0 * work1)
            henh3[idx] = xk * (1.0 + xe * xph[idx] / xkw)

            xk = 1.15e-2 * exp(2560.0 * work1)
            heo3[idx] = xk

            kh4 = (kh2 + kh3 * kh1 / xph[idx]) / ((1.0 + kh1 / xph[idx]) ** 2)
            ho2s = kh0 * xho2[idx] * patm * (1.0 + kh1 / xph[idx])
            r1h2o2 = kh4 * ho2s * ho2s

            if cloud_borne_flag != 0:
                r2h2o2 = r1h2o2 * xl / const0 * 1.0e6 / xam
            else:
                r2h2o2 = r1h2o2 * xl * const0 / xam

            if modal_aerosols_flag == 0:
                xh2o2[idx] = xh2o2[idx] + r2h2o2 * dtime

            px = hehno3[idx] * ra * tz * xl
            hno3g[idx] = (xhno3[idx] + xno3[idx]) / (1.0 + px)

            px = heh2o2[idx] * ra * tz * xl
            h2o2g = xh2o2[idx] / (1.0 + px)

            px = heso2[idx] * ra * tz * xl
            so2g = xso2[idx] / (1.0 + px)

            px = heo3[idx] * ra * tz * xl
            o3g = xo3[idx] / (1.0 + px)

            px = henh3[idx] * ra * tz * xl
            if id_nh3 > 0:
                nh3g[idx] = (xnh3[idx] + xnh4[idx]) / (1.0 + px)
            else:
                nh3g[idx] = 0.0

            rah2o2 = 8.0e4 * exp(-3650.0 * work1) / (0.1 + xph[idx])
            rao3 = 4.39e11 * exp(-4131.0 / tz) + 2.56e3 * exp(-996.0 / tz) / xph[idx]

            if xl >= 1.0e-8:
                if cloud_borne_flag != 0:
                    patm_x = patm
                else:
                    patm_x = 1.0

                if modal_aerosols_flag != 0:
                    pso4 = (
                        rah2o2
                        * 7.4e4
                        * exp(6621.0 * work1)
                        * h2o2g
                        * patm_x
                        * 1.23
                        * exp(3120.0 * work1)
                        * so2g
                        * patm_x
                    )
                else:
                    pso4 = rah2o2 * heh2o2[idx] * h2o2g * patm_x * heso2[idx] * so2g * patm_x

                pso4 = pso4 * xl / const0 / xhnm[idx]
                ccc = pso4 * dtime
                ccc = max(ccc, 1.0e-30)
                xso4_init[idx] = xso4[idx]

                if xh2o2[idx] > xso2[idx]:
                    if ccc > xso2[idx]:
                        xso4[idx] = xso4[idx] + xso2[idx]
                        if cloud_borne_flag != 0:
                            xh2o2[idx] = xh2o2[idx] - xso2[idx]
                            xso2[idx] = 1.0e-20
                        else:
                            xso2[idx] = 1.0e-20
                            xh2o2[idx] = xh2o2[idx] - xso2[idx]
                    else:
                        xso4[idx] = xso4[idx] + ccc
                        xh2o2[idx] = xh2o2[idx] - ccc
                        xso2[idx] = xso2[idx] - ccc
                else:
                    if ccc > xh2o2[idx]:
                        xso4[idx] = xso4[idx] + xh2o2[idx]
                        xso2[idx] = xso2[idx] - xh2o2[idx]
                        xh2o2[idx] = 1.0e-20
                    else:
                        xso4[idx] = xso4[idx] + ccc
                        xh2o2[idx] = xh2o2[idx] - ccc
                        xso2[idx] = xso2[idx] - ccc

                if modal_aerosols_flag != 0:
                    xdelso4hp[idx] = xso4[idx] - xso4_init[idx]

                pso4 = rao3 * heo3[idx] * o3g * patm_x * heso2[idx] * so2g * patm_x
                pso4 = pso4 * xl / const0 / xhnm[idx]
                ccc = pso4 * dtime
                ccc = max(ccc, 1.0e-30)
                xso4_init[idx] = xso4[idx]

                if ccc > xso2[idx]:
                    xso4[idx] = xso4[idx] + xso2[idx]
                    xso2[idx] = 1.0e-20
                else:
                    xso4[idx] = xso4[idx] + ccc
                    xso2[idx] = xso2[idx] - ccc


@export
def setsox_xph_lwc_diag_codon(
    ncol: int,
    pcols: int,
    pver: int,
    cldfrc_p: cobj,
    lwc_p: cobj,
    xph_p: cobj,
    xphlwc_p: cobj,
):
    cldfrc = Ptr[float](cldfrc_p)
    lwc = Ptr[float](lwc_p)
    xph = Ptr[float](xph_p)
    xphlwc = Ptr[float](xphlwc_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, ncol)
            idxp = _idx2(i, k, pcols)
            xphlwc[idx] = 0.0
            if cldfrc[idxp] >= 1.0e-5 and lwc[idx] >= 1.0e-8:
                xphlwc[idx] = -1.0 * log10(xph[idx]) * lwc[idx]


@inline
def _sox_cldaero_uptakerate(
    xl: float,
    cldnum: float,
    cfact: float,
    cldfrc: float,
    tfld: float,
    press: float,
    pi_val: float,
) -> float:
    num_cd = 1.0e-3 * cldnum * cfact / cldfrc
    num_cd = max(num_cd, 0.0)
    volx34pi_cd = xl * 0.75 / pi_val
    radxnum_cd = (volx34pi_cd * num_cd * num_cd) ** 0.3333333
    if radxnum_cd <= volx34pi_cd * 4.0e4:
        radxnum_cd = volx34pi_cd * 4.0e4
        rad_cd = 50.0e-4
    elif radxnum_cd >= volx34pi_cd * 4.0e8:
        radxnum_cd = volx34pi_cd * 4.0e8
        rad_cd = 0.5e-4
    else:
        rad_cd = radxnum_cd / num_cd

    gasdiffus = 0.557 * (tfld ** 1.75) / press
    gasspeed = 1.455e4 * sqrt(tfld / 98.0)
    knudsen = 3.0 * gasdiffus / (gasspeed * rad_cd)
    fuchs_sutugin = (0.4875 * (1.0 + knudsen)) / (knudsen * (1.184 + knudsen) + 0.4875)
    return 12.56637 * radxnum_cd * gasdiffus * fuchs_sutugin


@export
def sox_cldaero_update_core_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    ntot_amode: int,
    loffset: int,
    id_msa: int,
    id_h2so4: int,
    id_so2: int,
    id_h2o2: int,
    id_nh3: int,
    modeptr_accum: int,
    dtime: float,
    pi_val: float,
    cldfrc_p: cobj,
    xlwc_p: cobj,
    cldnum_p: cobj,
    cfact_p: cobj,
    tfld_p: cobj,
    press_p: cobj,
    delso4_hprxn_p: cobj,
    xh2so4_p: cobj,
    xso4_p: cobj,
    xso4_init_p: cobj,
    nh3g_p: cobj,
    xnh3_p: cobj,
    xnh4c_p: cobj,
    xmsa_p: cobj,
    xso2_p: cobj,
    xh2o2_p: cobj,
    qcw_p: cobj,
    qin_p: cobj,
    dqdt_aqso4_p: cobj,
    dqdt_aqh2so4_p: cobj,
    dqdt_aqhprxn_p: cobj,
    dqdt_aqo3rxn_p: cobj,
    faqgain_msa_p: cobj,
    faqgain_so4_p: cobj,
    qnum_c_p: cobj,
    numptrcw_amode_p: cobj,
    lptr_so4_cw_amode_p: cobj,
    lptr_msa_cw_amode_p: cobj,
    lptr_nh4_cw_amode_p: cobj,
):
    cldfrc = Ptr[float](cldfrc_p)
    xlwc = Ptr[float](xlwc_p)
    cldnum = Ptr[float](cldnum_p)
    cfact = Ptr[float](cfact_p)
    tfld = Ptr[float](tfld_p)
    press = Ptr[float](press_p)
    delso4_hprxn = Ptr[float](delso4_hprxn_p)
    xh2so4 = Ptr[float](xh2so4_p)
    xso4 = Ptr[float](xso4_p)
    xso4_init = Ptr[float](xso4_init_p)
    nh3g = Ptr[float](nh3g_p)
    xnh3 = Ptr[float](xnh3_p)
    xnh4c = Ptr[float](xnh4c_p)
    xmsa = Ptr[float](xmsa_p)
    xso2 = Ptr[float](xso2_p)
    xh2o2 = Ptr[float](xh2o2_p)
    qcw = Ptr[float](qcw_p)
    qin = Ptr[float](qin_p)
    dqdt_aqso4 = Ptr[float](dqdt_aqso4_p)
    dqdt_aqh2so4 = Ptr[float](dqdt_aqh2so4_p)
    dqdt_aqhprxn = Ptr[float](dqdt_aqhprxn_p)
    dqdt_aqo3rxn = Ptr[float](dqdt_aqo3rxn_p)
    faqgain_msa = Ptr[float](faqgain_msa_p)
    faqgain_so4 = Ptr[float](faqgain_so4_p)
    qnum_c = Ptr[float](qnum_c_p)
    numptrcw_amode = Ptr[int](numptrcw_amode_p)
    lptr_so4_cw_amode = Ptr[int](lptr_so4_cw_amode_p)
    lptr_msa_cw_amode = Ptr[int](lptr_msa_cw_amode_p)
    lptr_nh4_cw_amode = Ptr[int](lptr_nh4_cw_amode_p)

    for m in range(1, gas_pcnst + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                dqdt_aqso4[_idx3(i, k, m, ncol, pver)] = 0.0
                dqdt_aqh2so4[_idx3(i, k, m, ncol, pver)] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            dqdt_aqhprxn[_idx2(i, k, ncol)] = 0.0
            dqdt_aqo3rxn[_idx2(i, k, ncol)] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, ncol)
            idxp = _idx2(i, k, pcols)
            if cldfrc[idxp] >= 1.0e-5:
                xl = xlwc[idxp]
                if xl >= 1.0e-8:
                    delso4_o3rxn = xso4[idx] - xso4_init[idx]
                    if id_nh3 > 0:
                        delnh3 = nh3g[idx] - xnh3[idx]
                        delnh4 = -delnh3
                    else:
                        delnh3 = 0.0
                        delnh4 = 0.0

                    for n in range(1, ntot_amode + 1):
                        qnum_c[n - 1] = 0.0
                        l = numptrcw_amode[n - 1] - loffset
                        if l > 0:
                            qnum_c[n - 1] = max(0.0, qcw[_idx3(i, k, l, ncol, pver)])

                    n_accum = modeptr_accum
                    if n_accum <= 0:
                        n_accum = 1
                    qnum_c[n_accum - 1] = max(1.0e-10, qnum_c[n_accum - 1])

                    sumf = 0.0
                    for n in range(1, ntot_amode + 1):
                        faqgain_so4[n - 1] = 0.0
                        if lptr_so4_cw_amode[n - 1] > 0:
                            faqgain_so4[n - 1] = qnum_c[n - 1]
                            sumf = sumf + faqgain_so4[n - 1]

                    if sumf > 0.0:
                        for n in range(1, ntot_amode + 1):
                            faqgain_so4[n - 1] = faqgain_so4[n - 1] / sumf

                    ntot_msa_c = 0
                    sumf = 0.0
                    for n in range(1, ntot_amode + 1):
                        faqgain_msa[n - 1] = 0.0
                        if lptr_msa_cw_amode[n - 1] > 0:
                            faqgain_msa[n - 1] = qnum_c[n - 1]
                            ntot_msa_c = ntot_msa_c + 1
                        sumf = sumf + faqgain_msa[n - 1]

                    if sumf > 0.0:
                        for n in range(1, ntot_amode + 1):
                            faqgain_msa[n - 1] = faqgain_msa[n - 1] / sumf

                    uptkrate = _sox_cldaero_uptakerate(
                        xl, cldnum[idxp], cfact[idx], cldfrc[idxp], tfld[idxp], press[idxp], pi_val
                    )
                    uptkrate = (1.0 - exp(-min(100.0, dtime * uptkrate))) / dtime

                    dso4dt_gasuptk = xh2so4[idx] * uptkrate
                    if id_msa > 0:
                        dmsadt_gasuptk = xmsa[idx] * uptkrate
                    else:
                        dmsadt_gasuptk = 0.0

                    dmsadt_gasuptk_toso4 = 0.0
                    dmsadt_gasuptk_tomsa = dmsadt_gasuptk
                    if ntot_msa_c == 0:
                        dmsadt_gasuptk_tomsa = 0.0
                        dmsadt_gasuptk_toso4 = dmsadt_gasuptk

                    dso4dt_aqrxn = (delso4_o3rxn + delso4_hprxn[idx]) / dtime
                    dso4dt_hprxn = delso4_hprxn[idx] / dtime
                    fwetrem = 0.0

                    for n in range(1, ntot_amode + 1):
                        l = lptr_so4_cw_amode[n - 1] - loffset
                        if l > 0:
                            qidx = _idx3(i, k, l, ncol, pver)
                            dqdt_aqso4[qidx] = faqgain_so4[n - 1] * dso4dt_aqrxn * cldfrc[idxp]
                            dqdt_aqh2so4[qidx] = (
                                faqgain_so4[n - 1] * (dso4dt_gasuptk + dmsadt_gasuptk_toso4) * cldfrc[idxp]
                            )
                            dqdt_aq = dqdt_aqso4[qidx] + dqdt_aqh2so4[qidx]
                            dqdt_wr = -fwetrem * dqdt_aq
                            dqdt = dqdt_aq + dqdt_wr
                            qcw[qidx] = qcw[qidx] + dqdt * dtime

                        l = lptr_msa_cw_amode[n - 1] - loffset
                        if l > 0:
                            qidx = _idx3(i, k, l, ncol, pver)
                            dqdt_aq = faqgain_msa[n - 1] * dmsadt_gasuptk_tomsa * cldfrc[idxp]
                            dqdt_wr = -fwetrem * dqdt_aq
                            dqdt = dqdt_aq + dqdt_wr
                            qcw[qidx] = qcw[qidx] + dqdt * dtime

                        l = lptr_nh4_cw_amode[n - 1] - loffset
                        if l > 0:
                            qidx = _idx3(i, k, l, ncol, pver)
                            if delnh4 > 0.0:
                                dqdt_aq = faqgain_so4[n - 1] * delnh4 / dtime * cldfrc[idxp]
                                dqdt = dqdt_aq
                                qcw[qidx] = qcw[qidx] + dqdt * dtime
                            else:
                                dqdt = (
                                    qcw[qidx]
                                    / max(xnh4c[idxp], 1.0e-35)
                                    * delnh4
                                    / dtime
                                    * cldfrc[idxp]
                                )
                                qcw[qidx] = qcw[qidx] + dqdt * dtime

                    qin[_idx3(i, k, id_h2so4, ncol, pver)] = (
                        qin[_idx3(i, k, id_h2so4, ncol, pver)] - dso4dt_gasuptk * dtime * cldfrc[idxp]
                    )
                    if id_msa > 0:
                        qin[_idx3(i, k, id_msa, ncol, pver)] = (
                            qin[_idx3(i, k, id_msa, ncol, pver)] - dmsadt_gasuptk * dtime * cldfrc[idxp]
                        )

                    fwetrem = 0.0
                    dqdt_wr = -fwetrem * xso2[idx] / dtime * cldfrc[idxp]
                    dqdt_aq = -dso4dt_aqrxn * cldfrc[idxp]
                    dqdt = dqdt_aq + dqdt_wr
                    qin[_idx3(i, k, id_so2, ncol, pver)] = qin[_idx3(i, k, id_so2, ncol, pver)] + dqdt * dtime

                    fwetrem = 0.0
                    dqdt_wr = -fwetrem * xh2o2[idx] / dtime * cldfrc[idxp]
                    dqdt_aq = -dso4dt_hprxn * cldfrc[idxp]
                    dqdt = dqdt_aq + dqdt_wr
                    qin[_idx3(i, k, id_h2o2, ncol, pver)] = (
                        qin[_idx3(i, k, id_h2o2, ncol, pver)] + dqdt * dtime
                    )

                    if id_nh3 > 0:
                        dqdt_aq = delnh3 / dtime * cldfrc[idxp]
                        dqdt = dqdt_aq
                        qin[_idx3(i, k, id_nh3, ncol, pver)] = (
                            qin[_idx3(i, k, id_nh3, ncol, pver)] + dqdt * dtime
                        )

                    dqdt_aqhprxn[idx] = dso4dt_hprxn * cldfrc[idxp]
                    dqdt_aqo3rxn[idx] = (dso4dt_aqrxn - dso4dt_hprxn) * cldfrc[idxp]


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


@inline
def _gas_phase_chemdr_shell_h2o_setup(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    h2o_ndx: int,
    st80_25_ndx: int,
    aoa_nh_ndx: int,
    nh_5_ndx: int,
    nh_50_ndx: int,
    nh_50w_ndx: int,
    rad2deg: float,
    pmid_p: cobj,
    vmr_p: cobj,
    rlats_p: cobj,
    mmr_p: cobj,
    qh2o_p: cobj,
    h2ovmr_p: cobj,
):
    if st80_25_ndx > 0:
        gas_phase_chemdr_reset_ste_tracer_codon(
            ncol, pcols, pver, gas_pcnst, st80_25_ndx, 80.0e2, 200.0e-9, pmid_p, vmr_p
        )

    if aoa_nh_ndx > 0 and nh_5_ndx > 0 and nh_50_ndx > 0 and nh_50w_ndx > 0:
        rlats = Ptr[float](rlats_p)
        vmr = Ptr[float](vmr_p)
        for j in range(1, ncol + 1):
            xlat = rlats[j - 1] * rad2deg
            if xlat >= 30.0 and xlat <= 50.0:
                vmr[_idx3(j, pver, nh_5_ndx, ncol, pver)] = 100.0e-9
                vmr[_idx3(j, pver, nh_50_ndx, ncol, pver)] = 100.0e-9
                vmr[_idx3(j, pver, nh_50w_ndx, ncol, pver)] = 100.0e-9
                vmr[_idx3(j, pver, aoa_nh_ndx, ncol, pver)] = 0.0

    if h2o_ndx > 0:
        gas_phase_chemdr_load_h2o_fields_codon(
            ncol, pcols, pver, gas_pcnst, h2o_ndx, mmr_p, vmr_p, qh2o_p, h2ovmr_p
        )


@inline
def _gas_phase_chemdr_shell_post_solver(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    ndx_h2so4: int,
    o3_ndx: int,
    o3s_ndx: int,
    delt: float,
    troplev_p: cobj,
    vmr_p: cobj,
    o3s_loss_p: cobj,
    del_h2so4_gasprod_p: cobj,
):
    if o3_ndx > 0 and o3s_ndx > 0:
        gas_phase_chemdr_copy_o3_to_o3s_trop_codon(
            ncol, pcols, pver, gas_pcnst, troplev_p, o3_ndx, o3s_ndx, vmr_p
        )
        troplev = Ptr[int](troplev_p)
        vmr = Ptr[float](vmr_p)
        o3s_loss = Ptr[float](o3s_loss_p)
        for i in range(1, ncol + 1):
            for k in range(troplev[i - 1] + 1, pver + 1):
                idx = _idx3(i, k, o3s_ndx, ncol, pver)
                vmr[idx] = vmr[idx] * exp(-delt * o3s_loss[_idx2(i, k, ncol)])

    if ndx_h2so4 > 0:
        gas_phase_chemdr_update_h2so4_gasprod_codon(
            ncol, pver, gas_pcnst, ndx_h2so4, vmr_p, del_h2so4_gasprod_p
        )


@export
def gas_phase_chemdr_shell_codon(
    stage: int,
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    pcnst: int,
    rxntot: int,
    extcnt: int,
    nfs: int,
    indexm: int,
    has_linoz_data_flag: int,
    h2o_ndx: int,
    st80_25_ndx: int,
    aoa_nh_ndx: int,
    nh_5_ndx: int,
    nh_50_ndx: int,
    nh_50w_ndx: int,
    so4_ndx: int,
    st80_25_tau_ndx: int,
    ndx_h2so4: int,
    o3_ndx: int,
    o3s_ndx: int,
    synoz_ndx: int,
    aoa_nh_ext_ndx: int,
    rad2deg: float,
    delt: float,
    delt_inverse: float,
    rga: float,
    m2km: float,
    pa2mb: float,
    map2chm_p: cobj,
    troplev_p: cobj,
    ltrop_sol_p: cobj,
    zen_angle_p: cobj,
    sza_p: cobj,
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
    q_p: cobj,
    mmr_p: cobj,
    vmr_p: cobj,
    qh2o_p: cobj,
    h2ovmr_p: cobj,
    rlats_p: cobj,
    sulfate_p: cobj,
    satq_p: cobj,
    relhum_p: cobj,
    cldw_p: cobj,
    cwat_p: cobj,
    extfrc_p: cobj,
    invariants_p: cobj,
    het_rates_p: cobj,
    reaction_rates_p: cobj,
    del_h2so4_gasprod_p: cobj,
    vmr0_p: cobj,
    o3s_loss_p: cobj,
    mmr_tend_p: cobj,
    mmr_new_p: cobj,
    qtend_p: cobj,
    tfld_p: cobj,
    tvs_p: cobj,
    sflx_p: cobj,
    ufld_p: cobj,
    vfld_p: cobj,
    wind_speed_p: cobj,
    precc_p: cobj,
    precl_p: cobj,
    prect_p: cobj,
    cflx_p: cobj,
    drydepflx_p: cobj,
):
    if stage == 1:
        gas_phase_chemdr_prepare_sza_codon(ncol, rad2deg, zen_angle_p, sza_p)
    elif stage == 2:
        gas_phase_chemdr_prepare_state_codon(
            ncol, pcols, pver, rga, m2km, pa2mb, phis_p, zi_p, zm_p, pmid_p, zsurf_p, zintr_p, zmidr_p,
            zmid_p, zint_p, pmb_p
        )
        gas_phase_chemdr_load_mmr_codon(ncol, pcols, pver, pcnst, map2chm_p, q_p, mmr_p)
    elif stage == 3:
        _gas_phase_chemdr_shell_h2o_setup(
            ncol, pcols, pver, gas_pcnst, h2o_ndx, st80_25_ndx, aoa_nh_ndx, nh_5_ndx, nh_50_ndx,
            nh_50w_ndx, rad2deg, pmid_p, vmr_p, rlats_p, mmr_p, qh2o_p, h2ovmr_p
        )
    elif stage == 4:
        gas_phase_chemdr_zero_sulfate_codon(ncol, pver, sulfate_p)
    elif stage == 5:
        gas_phase_chemdr_load_prognostic_sulfate_codon(ncol, pver, gas_pcnst, so4_ndx, vmr_p, sulfate_p)
    elif stage == 6:
        gas_phase_chemdr_clip_sulfate_codon(ncol, pcols, pver, troplev_p, sulfate_p)
    elif stage == 7:
        gas_phase_chemdr_compute_relhum_codon(ncol, pver, h2ovmr_p, satq_p, relhum_p)
        gas_phase_chemdr_copy_cldw_to_cwat_codon(ncol, pcols, pver, cldw_p, cwat_p)
    elif stage == 8:
        gas_phase_chemdr_normalize_extfrc_codon(
            ncol, pver, extcnt, synoz_ndx, aoa_nh_ext_ndx, indexm, extfrc_p, invariants_p
        )
    elif stage == 9:
        gas_phase_chemdr_zero_het_rates_codon(ncol, pver, gas_pcnst, het_rates_p)
    elif stage == 10:
        gas_phase_chemdr_zero_st80_tau_codon(ncol, pver, rxntot, st80_25_tau_ndx, troplev_p, reaction_rates_p)
    elif stage == 11:
        gas_phase_chemdr_set_ltrop_sol_codon(ncol, has_linoz_data_flag, troplev_p, ltrop_sol_p)
        gas_phase_chemdr_init_h2so4_gasprod_codon(
            ncol, pver, gas_pcnst, ndx_h2so4, vmr_p, del_h2so4_gasprod_p
        )
        gas_phase_chemdr_store_vmr0_codon(ncol, pver, gas_pcnst, vmr_p, vmr0_p)
    elif stage == 12:
        _gas_phase_chemdr_shell_post_solver(
            ncol, pcols, pver, gas_pcnst, ndx_h2so4, o3_ndx, o3s_ndx, delt, troplev_p, vmr_p, o3s_loss_p,
            del_h2so4_gasprod_p
        )
    elif stage == 13:
        gas_phase_chemdr_finalize_tendencies_codon(
            ncol, pcols, pver, gas_pcnst, pcnst, delt_inverse, map2chm_p, mmr_p, mmr_tend_p, mmr_new_p,
            qtend_p
        )
        gas_phase_chemdr_compute_tvs_codon(ncol, pcols, pver, tfld_p, qh2o_p, tvs_p)
        gas_phase_chemdr_zero_sflx_codon(pcols, gas_pcnst, sflx_p)
    elif stage == 14:
        gas_phase_chemdr_compute_wind_speed_codon(ncol, pcols, pver, ufld_p, vfld_p, wind_speed_p)
        gas_phase_chemdr_compute_prect_codon(ncol, pcols, precc_p, precl_p, prect_p)
    elif stage == 15:
        gas_phase_chemdr_store_drydep_codon(
            ncol, pcols, gas_pcnst, pcnst, map2chm_p, sflx_p, cflx_p, drydepflx_p
        )


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
