from math import cos, exp, log, sin
from chemistry_common_codon import _idx2, _idx3, _idx4, _idx5

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

def jlong_init_solar_batch_codon(
    data_nw: int,
    nw: int,
    data_we_p: cobj,
    wc_p: cobj,
    wlintv_p: cobj,
    we_p: cobj,
    data_etf_p: cobj,
    etfphot_p: cobj,
):
    jlong_init_set_we_codon(nw, wc_p, wlintv_p, we_p)
    rebin_codon(data_nw, nw, data_we_p, we_p, data_etf_p, etfphot_p)

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

def jlong_get_xsqy_meta_batch_codon(
    phtcnt: int,
    lng_indexer_p: cobj,
    numj_p: cobj,
    read_varids_p: cobj,
):
    lng_indexer = Ptr[int](lng_indexer_p)
    numj_out = Ptr[int](numj_p)
    read_varids = Ptr[int](read_varids_p)

    count = 0
    for m in range(1, phtcnt + 1):
        if lng_indexer[m - 1] > 0:
            seen = 0
            for i in range(1, m):
                if lng_indexer[i - 1] == lng_indexer[m - 1]:
                    seen = 1
                    break
            if seen == 0:
                count += 1
                read_varids[count - 1] = lng_indexer[m - 1]

    numj_out[0] = count

    for m in range(1, phtcnt + 1):
        value = lng_indexer[m - 1]
        if value > 0:
            for ndx in range(1, count + 1):
                if read_varids[ndx - 1] == value:
                    lng_indexer[m - 1] = ndx
                    break

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

def jlong_get_xsqy_dprs_codon(
    np_xs: int,
    prs_p: cobj,
    dprs_p: cobj,
):
    prs = Ptr[float](prs_p)
    dprs = Ptr[float](dprs_p)

    for i in range(1, np_xs):
        dprs[i - 1] = 1.0 / (prs[i - 1] - prs[i])

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

def jlong_get_rsf_postread_batch_codon(
    nw: int,
    nump: int,
    numsza: int,
    numalb: int,
    numcolo3: int,
    use_bde_flag: int,
    hc_val: float,
    wc_o2_b_val: float,
    wc_o3_a_val: float,
    wc_o3_b_val: float,
    wc_p: cobj,
    p_p: cobj,
    sza_p: cobj,
    alb_p: cobj,
    o3rat_p: cobj,
    bde_o2_b_p: cobj,
    bde_o3_a_p: cobj,
    bde_o3_b_p: cobj,
    del_p_p: cobj,
    del_sza_p: cobj,
    del_alb_p: cobj,
    del_o3rat_p: cobj,
):
    wc = Ptr[float](wc_p)
    p = Ptr[float](p_p)
    sza = Ptr[float](sza_p)
    alb = Ptr[float](alb_p)
    o3rat = Ptr[float](o3rat_p)
    bde_o2_b = Ptr[float](bde_o2_b_p)
    bde_o3_a = Ptr[float](bde_o3_a_p)
    bde_o3_b = Ptr[float](bde_o3_b_p)
    del_p = Ptr[float](del_p_p)
    del_sza = Ptr[float](del_sza_p)
    del_alb = Ptr[float](del_alb_p)
    del_o3rat = Ptr[float](del_o3rat_p)

    if use_bde_flag != 0:
        for i in range(1, nw + 1):
            wc_i = wc[i - 1]
            bde_o2_b[i - 1] = max(0.0, hc_val * (wc_o2_b_val - wc_i) / (wc_o2_b_val * wc_i))
            bde_o3_a[i - 1] = max(0.0, hc_val * (wc_o3_a_val - wc_i) / (wc_o3_a_val * wc_i))
            bde_o3_b[i - 1] = max(0.0, hc_val * (wc_o3_b_val - wc_i) / (wc_o3_b_val * wc_i))
    else:
        for i in range(1, nw + 1):
            wc_i = wc[i - 1]
            bde_o2_b[i - 1] = hc_val / wc_i
            bde_o3_a[i - 1] = hc_val / wc_i
            bde_o3_b[i - 1] = hc_val / wc_i

    for i in range(1, nump):
        del_p[i - 1] = 1.0 / abs(p[i - 1] - p[i])
    for i in range(1, numsza):
        del_sza[i - 1] = 1.0 / (sza[i] - sza[i - 1])
    for i in range(1, numalb):
        del_alb[i - 1] = 1.0 / (alb[i] - alb[i - 1])
    for i in range(1, numcolo3):
        del_o3rat[i - 1] = 1.0 / (o3rat[i] - o3rat[i - 1])


def jlong_prep_init_solar_batch_codon(
    data_nw: int,
    nw: int,
    data_we_p: cobj,
    wc_p: cobj,
    wlintv_p: cobj,
    we_p: cobj,
    data_etf_p: cobj,
    etfphot_p: cobj,
):
    jlong_init_solar_batch_codon(data_nw, nw, data_we_p, wc_p, wlintv_p, we_p, data_etf_p, etfphot_p)


def jlong_prep_get_xsqy_meta_batch_codon(
    phtcnt: int,
    lng_indexer_p: cobj,
    numj_p: cobj,
    read_varids_p: cobj,
):
    jlong_get_xsqy_meta_batch_codon(phtcnt, lng_indexer_p, numj_p, read_varids_p)


def jlong_prep_get_xsqy_dprs_codon(
    np_xs: int,
    prs_p: cobj,
    dprs_p: cobj,
):
    jlong_get_xsqy_dprs_codon(np_xs, prs_p, dprs_p)


def jlong_prep_get_rsf_scale_codon(
    nw: int,
    nump: int,
    numsza: int,
    numcolo3: int,
    numalb: int,
    wlintv_p: cobj,
    rsf_tab_p: cobj,
):
    jlong_get_rsf_scale_codon(nw, nump, numsza, numcolo3, numalb, wlintv_p, rsf_tab_p)


def jlong_prep_get_rsf_postread_batch_codon(
    nw: int,
    nump: int,
    numsza: int,
    numalb: int,
    numcolo3: int,
    use_bde_flag: int,
    hc_val: float,
    wc_o2_b_val: float,
    wc_o3_a_val: float,
    wc_o3_b_val: float,
    wc_p: cobj,
    p_p: cobj,
    sza_p: cobj,
    alb_p: cobj,
    o3rat_p: cobj,
    bde_o2_b_p: cobj,
    bde_o3_a_p: cobj,
    bde_o3_b_p: cobj,
    del_p_p: cobj,
    del_sza_p: cobj,
    del_alb_p: cobj,
    del_o3rat_p: cobj,
):
    jlong_get_rsf_postread_batch_codon(
        nw,
        nump,
        numsza,
        numalb,
        numcolo3,
        use_bde_flag,
        hc_val,
        wc_o2_b_val,
        wc_o3_a_val,
        wc_o3_b_val,
        wc_p,
        p_p,
        sza_p,
        alb_p,
        o3rat_p,
        bde_o2_b_p,
        bde_o3_a_p,
        bde_o3_b_p,
        del_p_p,
        del_sza_p,
        del_alb_p,
        del_o3rat_p,
    )


def jlong_prep_batch_codon(
    stage: int,
    data_nw: int,
    nw: int,
    phtcnt: int,
    np_xs: int,
    nump: int,
    numsza: int,
    numalb: int,
    numcolo3: int,
    use_bde_flag: int,
    hc_val: float,
    wc_o2_b_val: float,
    wc_o3_a_val: float,
    wc_o3_b_val: float,
    data_we_p: cobj,
    wc_p: cobj,
    wlintv_p: cobj,
    we_p: cobj,
    data_etf_p: cobj,
    etfphot_p: cobj,
    lng_indexer_p: cobj,
    numj_p: cobj,
    read_varids_p: cobj,
    prs_p: cobj,
    dprs_p: cobj,
    rsf_tab_p: cobj,
    p_p: cobj,
    sza_p: cobj,
    alb_p: cobj,
    o3rat_p: cobj,
    bde_o2_b_p: cobj,
    bde_o3_a_p: cobj,
    bde_o3_b_p: cobj,
    del_p_p: cobj,
    del_sza_p: cobj,
    del_alb_p: cobj,
    del_o3rat_p: cobj,
):
    if stage == 1:
        jlong_prep_init_solar_batch_codon(data_nw, nw, data_we_p, wc_p, wlintv_p, we_p, data_etf_p, etfphot_p)
    elif stage == 2:
        jlong_prep_get_xsqy_meta_batch_codon(phtcnt, lng_indexer_p, numj_p, read_varids_p)
    elif stage == 3:
        jlong_prep_get_xsqy_dprs_codon(np_xs, prs_p, dprs_p)
    elif stage == 4:
        jlong_prep_get_rsf_scale_codon(nw, nump, numsza, numcolo3, numalb, wlintv_p, rsf_tab_p)
    elif stage == 5:
        jlong_prep_get_rsf_postread_batch_codon(
            nw,
            nump,
            numsza,
            numalb,
            numcolo3,
            use_bde_flag,
            hc_val,
            wc_o2_b_val,
            wc_o3_a_val,
            wc_o3_b_val,
            wc_p,
            p_p,
            sza_p,
            alb_p,
            o3rat_p,
            bde_o2_b_p,
            bde_o3_a_p,
            bde_o3_b_p,
            del_p_p,
            del_sza_p,
            del_alb_p,
            del_o3rat_p,
        )


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

def table_photo_zero_finalize_batch_codon(
    stage: int,
    ncol: int,
    pver: int,
    phtcnt: int,
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

    if stage == 0:
        for m in range(1, phtcnt + 1):
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    photos[_idx3(i, k, m, ncol, pver)] = 0.0
        return

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

def table_photo_daylight_prepare_batch_codon(
    ncol: int,
    pcols: int,
    pver: int,
    ncol_abs: int,
    nfs: int,
    gas_pcnst: int,
    i_col: int,
    p1: int,
    p2: int,
    do_jshort_flag: int,
    ptop_gt_10_flag: int,
    o_is_inv_flag: int,
    o2_is_inv_flag: int,
    o3_is_inv_flag: int,
    n2_is_inv_flag: int,
    no_is_inv_flag: int,
    o_ndx: int,
    o2_ndx: int,
    o3_inv_ndx: int,
    o3_ndx: int,
    n2_ndx: int,
    no_ndx: int,
    indexm: int,
    pa2mb: float,
    pmid_p: cobj,
    pdel_p: cobj,
    col_dens_p: cobj,
    lwc_p: cobj,
    clouds_p: cobj,
    temper_p: cobj,
    zmid_p: cobj,
    zint_p: cobj,
    vmr_p: cobj,
    invariants_p: cobj,
    parg_p: cobj,
    colo3_p: cobj,
    fac1_p: cobj,
    lwc_line_p: cobj,
    cld_line_p: cobj,
    tline_p: cobj,
    zarg_p: cobj,
    o_den_p: cobj,
    o2_den_p: cobj,
    o3_den_p: cobj,
    no_den_p: cobj,
    n2_den_p: cobj,
):
    pmid = Ptr[float](pmid_p)
    pdel = Ptr[float](pdel_p)
    col_dens = Ptr[float](col_dens_p)
    lwc = Ptr[float](lwc_p)
    clouds = Ptr[float](clouds_p)
    temper = Ptr[float](temper_p)
    zmid = Ptr[float](zmid_p)
    zint = Ptr[float](zint_p)
    vmr = Ptr[float](vmr_p)
    invariants = Ptr[float](invariants_p)
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

    if do_jshort_flag == 0:
        return

    o_den = Ptr[float](o_den_p)
    o2_den = Ptr[float](o2_den_p)
    o3_den = Ptr[float](o3_den_p)
    no_den = Ptr[float](no_den_p)
    n2_den = Ptr[float](n2_den_p)

    for k in range(1, pver + 1):
        dst_k = p1 + k - 1
        if o_is_inv_flag != 0:
            o_den[dst_k - 1] = invariants[_idx3(i_col, k, o_ndx, ncol, pver)]
        else:
            o_den[dst_k - 1] = (
                vmr[_idx3(i_col, k, o_ndx, ncol, pver)] * invariants[_idx3(i_col, k, indexm, ncol, pver)]
            )

        if o2_is_inv_flag != 0:
            o2_den[dst_k - 1] = invariants[_idx3(i_col, k, o2_ndx, ncol, pver)]
        else:
            o2_den[dst_k - 1] = (
                vmr[_idx3(i_col, k, o2_ndx, ncol, pver)] * invariants[_idx3(i_col, k, indexm, ncol, pver)]
            )

        if o3_is_inv_flag != 0:
            o3_den[dst_k - 1] = invariants[_idx3(i_col, k, o3_inv_ndx, ncol, pver)]
        else:
            o3_den[dst_k - 1] = (
                vmr[_idx3(i_col, k, o3_ndx, ncol, pver)] * invariants[_idx3(i_col, k, indexm, ncol, pver)]
            )

        if n2_is_inv_flag != 0:
            n2_den[dst_k - 1] = invariants[_idx3(i_col, k, n2_ndx, ncol, pver)]
        else:
            n2_den[dst_k - 1] = (
                vmr[_idx3(i_col, k, n2_ndx, ncol, pver)] * invariants[_idx3(i_col, k, indexm, ncol, pver)]
            )

        if no_is_inv_flag != 0:
            no_den[dst_k - 1] = invariants[_idx3(i_col, k, no_ndx, ncol, pver)]
        else:
            no_den[dst_k - 1] = (
                vmr[_idx3(i_col, k, no_ndx, ncol, pver)] * invariants[_idx3(i_col, k, indexm, ncol, pver)]
            )

    if ptop_gt_10_flag != 0:
        ideltaZkm = 1.0 / (zint[_idx2(i_col, 1, ncol)] - zint[_idx2(i_col, 2, ncol)])
        o3_den[0] = o3_den[1] * 7.0 * ideltaZkm
        o2_den[0] = o2_den[1] * 7.0 * ideltaZkm
        no_den[0] = no_den[1] * 0.9
        n2_den[0] = n2_den[1] * 0.9
        tline[0] = tline[1] + 5.0
        zarg[0] = zarg[1] + (zint[_idx2(i_col, 1, ncol)] - zint[_idx2(i_col, 2, ncol)])

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

def table_photo_scale_cld_mult_codon(
    pver: int,
    esfact: float,
    cld_mult_p: cobj,
):
    cld_mult = Ptr[float](cld_mult_p)

    for k in range(1, pver + 1):
        cld_mult[k - 1] = esfact * cld_mult[k - 1]

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

def photo_prep_batch_codon(
    stage: int,
    pinterp: float,
    calday: float,
    n_exo_levs: int,
    levs_p: cobj,
    days_p: cobj,
    ki_p: cobj,
    next_p: cobj,
    last_p: cobj,
    delp_p: cobj,
    dels_p: cobj,
):
    if stage == 1:
        photo_inti_fixed_press_setup_codon(pinterp, n_exo_levs, levs_p, ki_p, delp_p)
    elif stage == 2:
        photo_timestep_init_exo_time_codon(calday, days_p, next_p, last_p, dels_p)


def photo_prep_fixed_press_setup_codon(
    pinterp: float,
    n_exo_levs: int,
    levs_p: cobj,
    ki_p: cobj,
    delp_p: cobj,
):
    photo_inti_fixed_press_setup_codon(pinterp, n_exo_levs, levs_p, ki_p, delp_p)


def photo_prep_timestep_init_exo_time_codon(
    calday: float,
    days_p: cobj,
    next_p: cobj,
    last_p: cobj,
    dels_p: cobj,
):
    photo_timestep_init_exo_time_codon(calday, days_p, next_p, last_p, dels_p)


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

def table_photo_prejlong_batch_codon(
    ncol: int,
    pcols: int,
    pver: int,
    ncol_abs: int,
    nfs: int,
    gas_pcnst: int,
    i_col: int,
    p1: int,
    p2: int,
    do_jshort_flag: int,
    ptop_gt_10_flag: int,
    o_is_inv_flag: int,
    o2_is_inv_flag: int,
    o3_is_inv_flag: int,
    n2_is_inv_flag: int,
    no_is_inv_flag: int,
    o_ndx: int,
    o2_ndx: int,
    o3_inv_ndx: int,
    o3_ndx: int,
    n2_ndx: int,
    no_ndx: int,
    indexm: int,
    pa2mb: float,
    zen_angle: float,
    srf_alb: float,
    rgrav: float,
    pmid_p: cobj,
    pdel_p: cobj,
    col_dens_p: cobj,
    lwc_p: cobj,
    clouds_p: cobj,
    temper_p: cobj,
    zmid_p: cobj,
    zint_p: cobj,
    vmr_p: cobj,
    invariants_p: cobj,
    parg_p: cobj,
    colo3_p: cobj,
    fac1_p: cobj,
    lwc_line_p: cobj,
    cld_line_p: cobj,
    tline_p: cobj,
    zarg_p: cobj,
    o_den_p: cobj,
    o2_den_p: cobj,
    o3_den_p: cobj,
    no_den_p: cobj,
    n2_den_p: cobj,
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
    cloud_fac1_p: cobj,
    cloud_fac2_p: cobj,
):
    table_photo_daylight_prepare_batch_codon(
        ncol,
        pcols,
        pver,
        ncol_abs,
        nfs,
        gas_pcnst,
        i_col,
        p1,
        p2,
        do_jshort_flag,
        ptop_gt_10_flag,
        o_is_inv_flag,
        o2_is_inv_flag,
        o3_is_inv_flag,
        n2_is_inv_flag,
        no_is_inv_flag,
        o_ndx,
        o2_ndx,
        o3_inv_ndx,
        o3_ndx,
        n2_ndx,
        no_ndx,
        indexm,
        pa2mb,
        pmid_p,
        pdel_p,
        col_dens_p,
        lwc_p,
        clouds_p,
        temper_p,
        zmid_p,
        zint_p,
        vmr_p,
        invariants_p,
        parg_p,
        colo3_p,
        fac1_p,
        lwc_line_p,
        cld_line_p,
        tline_p,
        zarg_p,
        o_den_p,
        o2_den_p,
        o3_den_p,
        no_den_p,
        n2_den_p,
    )
    table_photo_cloud_mod_batch_codon(
        pver,
        zen_angle,
        srf_alb,
        rgrav,
        cld_line_p,
        lwc_line_p,
        fac1_p,
        eff_alb_p,
        cld_mult_p,
        del_lwp_p,
        del_tau_p,
        above_tau_p,
        below_tau_p,
        above_cld_p,
        below_cld_p,
        above_tra_p,
        below_tra_p,
        cloud_fac1_p,
        cloud_fac2_p,
    )

def table_photo_cloud_mod_batch_codon(
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

def table_photo_postcloud_batch_codon(
    ncol: int,
    pver: int,
    phtcnt: int,
    nlng: int,
    i_col: int,
    jno_ndx: int,
    jho2no2_ndx: int,
    do_jshort_flag: int,
    has_o2_col_flag: int,
    has_o3_col_flag: int,
    zen_angle: float,
    esfact: float,
    photos_p: cobj,
    lng_prates_p: cobj,
    cld_mult_p: cobj,
    col_dens_p: cobj,
    lng_indexer_p: cobj,
    alias_mult2_p: cobj,
):
    photos = Ptr[float](photos_p)
    lng_prates = Ptr[float](lng_prates_p)
    cld_mult = Ptr[float](cld_mult_p)
    col_dens = Ptr[float](col_dens_p)
    lng_indexer = Ptr[int](lng_indexer_p)
    alias_mult2 = Ptr[float](alias_mult2_p)

    for k in range(1, pver + 1):
        cld_mult[k - 1] = esfact * cld_mult[k - 1]

    if nlng > 0:
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
                            photos[_idx3(i_col, k, m, ncol, pver)]
                            + alias_factor * lng_prates[_idx2(idx_lng, k, nlng)]
                        ) * cld_mult[k - 1]

    if jno_ndx > 0 and do_jshort_flag == 0:
        if has_o2_col_flag != 0 and has_o3_col_flag != 0:
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
