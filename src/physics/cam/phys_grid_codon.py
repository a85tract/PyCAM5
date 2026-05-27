def phys_grid_count_valid_cols_codon(ngcols: int, clon_d_p: cobj) -> int:
    clon_d = Ptr[float](clon_d_p)

    ngcols_p = 0
    for i in range(ngcols):
        if clon_d[i] < 100000.0:
            ngcols_p += 1
    return ngcols_p


def phys_grid_count_unique_sorted_real_codon(ncols: int, cdex_p: cobj, coord_p: cobj) -> int:
    cdex = Ptr[i32](cdex_p)
    coord = Ptr[float](coord_p)

    if ncols <= 0:
        return 0

    unique_count = 1
    prev = coord[int(cdex[0]) - 1]
    for i in range(1, ncols):
        value = coord[int(cdex[i]) - 1]
        if value > prev:
            unique_count += 1
            prev = value
    return unique_count


def phys_grid_fill_unique_sorted_real_codon(
    ncols: int,
    cdex_p: cobj,
    coord_p: cobj,
    unique_p: cobj,
    counts_p: cobj,
):
    cdex = Ptr[i32](cdex_p)
    coord = Ptr[float](coord_p)
    unique = Ptr[float](unique_p)
    counts = Ptr[i32](counts_p)

    if ncols <= 0:
        return

    pre_i = 1
    unique_tot = 1
    unique[0] = coord[int(cdex[0]) - 1]

    for i0 in range(1, ncols):
        i_fortran = i0 + 1
        value = coord[int(cdex[i0]) - 1]
        if value > unique[unique_tot - 1]:
            counts[unique_tot - 1] = i32(i_fortran - pre_i)
            pre_i = i_fortran
            unique_tot += 1
            unique[unique_tot - 1] = value

    counts[unique_tot - 1] = i32((ncols + 1) - pre_i)


def phys_grid_prefix_counts_codon(n: int, counts_p: cobj, idx_p: cobj):
    counts = Ptr[i32](counts_p)
    idx = Ptr[i32](idx_p)

    if n <= 0:
        return

    idx[0] = i32(1)
    for j in range(1, n):
        idx[j] = i32(int(idx[j - 1]) + int(counts[j - 1]))


def phys_grid_fill_real_pair_codon(
    n: int,
    first_value: float,
    second_value: float,
    first_p: cobj,
    second_p: cobj,
):
    first = Ptr[float](first_p)
    second = Ptr[float](second_p)

    for i in range(n):
        first[i] = first_value
        second[i] = second_value


def phys_grid_init_lat_map_codon(
    ngcols: int,
    ncols_p: int,
    clat_tot: int,
    has_latlon_map: int,
    cdex_p: cobj,
    clat_d_p: cobj,
    clat_p_p: cobj,
    lat_p_p: cobj,
    dyn_to_latlon_gcol_map_p: cobj,
    latlon_to_dyn_gcol_map_p: cobj,
):
    cdex = Ptr[i32](cdex_p)
    clat_d = Ptr[float](clat_d_p)
    clat_p = Ptr[float](clat_p_p)
    lat_p = Ptr[i32](lat_p_p)
    dyn_to_latlon_gcol_map = Ptr[i32](dyn_to_latlon_gcol_map_p)
    latlon_to_dyn_gcol_map = Ptr[i32](latlon_to_dyn_gcol_map_p)

    for i in range(ngcols):
        lat_p[i] = i32(-1)
        dyn_to_latlon_gcol_map[i] = i32(-1)

    clat_p_dex = 1
    for i in range(ncols_p):
        dyn_idx = int(cdex[i])
        if has_latlon_map != 0:
            latlon_to_dyn_gcol_map[i] = i32(dyn_idx)
        dyn_to_latlon_gcol_map[dyn_idx - 1] = i32(i + 1)

        while clat_p[clat_p_dex - 1] < clat_d[dyn_idx - 1] and clat_p_dex < clat_tot:
            clat_p_dex += 1
        lat_p[dyn_idx - 1] = i32(clat_p_dex)


def phys_grid_init_lon_map_codon(
    ngcols: int,
    ncols_p: int,
    clon_tot: int,
    has_lonlat_map: int,
    cdex_p: cobj,
    clon_d_p: cobj,
    clon_p_p: cobj,
    lon_p_p: cobj,
    lonlat_to_dyn_gcol_map_p: cobj,
):
    cdex = Ptr[i32](cdex_p)
    clon_d = Ptr[float](clon_d_p)
    clon_p = Ptr[float](clon_p_p)
    lon_p = Ptr[i32](lon_p_p)
    lonlat_to_dyn_gcol_map = Ptr[i32](lonlat_to_dyn_gcol_map_p)

    for i in range(ngcols):
        lon_p[i] = i32(-1)

    clon_p_dex = 1
    for i in range(ncols_p):
        dyn_idx = int(cdex[i])
        if has_lonlat_map != 0:
            lonlat_to_dyn_gcol_map[i] = i32(dyn_idx)
        while clon_p[clon_p_dex - 1] < clon_d[dyn_idx - 1] and clon_p_dex < clon_tot:
            clon_p_dex += 1
        lon_p[dyn_idx - 1] = i32(clon_p_dex)


def phys_grid_zero_proc_counts_codon(npes: int, chunk_counts_p: cobj, col_counts_p: cobj):
    chunk_counts = Ptr[i32](chunk_counts_p)
    col_counts = Ptr[i32](col_counts_p)

    for p in range(npes):
        chunk_counts[p] = i32(0)
        col_counts[p] = i32(0)


def phys_grid_proc_prefix_offsets_codon(
    npes: int,
    start_value: int,
    set_final: int,
    chunk_counts_p: cobj,
    col_counts_p: cobj,
    pchunkid_p: cobj,
    gs_col_offset_p: cobj,
):
    chunk_counts = Ptr[i32](chunk_counts_p)
    col_counts = Ptr[i32](col_counts_p)
    pchunkid = Ptr[i32](pchunkid_p)
    gs_col_offset = Ptr[i32](gs_col_offset_p)

    if npes <= 0:
        return

    pchunkid[0] = i32(start_value)
    gs_col_offset[0] = i32(start_value)
    for p in range(1, npes):
        pchunkid[p] = i32(int(pchunkid[p - 1]) + int(chunk_counts[p - 1]))
        gs_col_offset[p] = i32(int(gs_col_offset[p - 1]) + int(col_counts[p - 1]))

    if set_final != 0:
        pchunkid[npes] = i32(int(pchunkid[npes - 1]) + int(chunk_counts[npes - 1]))
        gs_col_offset[npes] = i32(int(gs_col_offset[npes - 1]) + int(col_counts[npes - 1]))


def phys_grid_process_bin_sort_codon(
    nchunks: int,
    lastblock: int,
    chunk_owner_p: cobj,
    chunk_ncols_p: cobj,
    pchunkid_p: cobj,
    gs_col_offset_p: cobj,
    chunk_lcid_p: cobj,
    pgcol_chunk_p: cobj,
    pgcol_ccol_p: cobj,
):
    chunk_owner = Ptr[i32](chunk_owner_p)
    chunk_ncols = Ptr[i32](chunk_ncols_p)
    pchunkid = Ptr[i32](pchunkid_p)
    gs_col_offset = Ptr[i32](gs_col_offset_p)
    chunk_lcid = Ptr[i32](chunk_lcid_p)
    pgcol_chunk = Ptr[i32](pgcol_chunk_p)
    pgcol_ccol = Ptr[i32](pgcol_ccol_p)

    for cid0 in range(nchunks):
        cid = cid0 + 1
        p = int(chunk_owner[cid0])
        pchunkid[p] = i32(int(pchunkid[p]) + 1)
        chunk_lcid[cid0] = i32(int(pchunkid[p]) + lastblock)

        curgcol = int(gs_col_offset[p])
        for i in range(1, int(chunk_ncols[cid0]) + 1):
            curgcol += 1
            pgcol_idx = curgcol - 1
            pgcol_chunk[pgcol_idx] = i32(cid)
            pgcol_ccol[pgcol_idx] = i32(i)
        gs_col_offset[p] = i32(curgcol)


def phys_grid_lchunk_gcol_copy_codon(ncols: int, src_gcol_p: cobj, dst_gcol_p: cobj):
    src_gcol = Ptr[i32](src_gcol_p)
    dst_gcol = Ptr[i32](dst_gcol_p)

    for i in range(ncols):
        dst_gcol[i] = src_gcol[i]


def phys_grid_lchunk_area_wght_codon(
    ncols: int,
    gcol_p: cobj,
    area_d_p: cobj,
    wght_d_p: cobj,
    area_p: cobj,
    wght_p: cobj,
):
    gcol = Ptr[i32](gcol_p)
    area_d = Ptr[float](area_d_p)
    wght_d = Ptr[float](wght_d_p)
    area = Ptr[float](area_p)
    wght = Ptr[float](wght_p)

    for i in range(ncols):
        dyn_idx = int(gcol[i]) - 1
        area[i] = area_d[dyn_idx]
        wght[i] = wght_d[dyn_idx]


def phys_grid_count_smp_procs_codon(
    npes: int,
    nsmpx: int,
    proc_smp_mapx_p: cobj,
    nsmpprocs_p: cobj,
) -> int:
    proc_smp_mapx = Ptr[i32](proc_smp_mapx_p)
    nsmpprocs = Ptr[i32](nsmpprocs_p)

    for smp in range(nsmpx):
        nsmpprocs[smp] = i32(0)

    for p in range(npes):
        smp = int(proc_smp_mapx[p])
        nsmpprocs[smp] = i32(int(nsmpprocs[smp]) + 1)

    max_count = 0
    for smp in range(nsmpx):
        if int(nsmpprocs[smp]) > max_count:
            max_count = int(nsmpprocs[smp])
    return max_count


def phys_grid_create_chunks_thread_counts_codon(
    npes: int,
    nsmpx: int,
    proc_smp_mapx_p: cobj,
    npthreads_p: cobj,
    nsmpthreads_p: cobj,
):
    proc_smp_mapx = Ptr[i32](proc_smp_mapx_p)
    npthreads = Ptr[i32](npthreads_p)
    nsmpthreads = Ptr[i32](nsmpthreads_p)

    for smp in range(nsmpx):
        nsmpthreads[smp] = i32(0)

    for p in range(npes):
        smp = int(proc_smp_mapx[p])
        nsmpthreads[smp] = i32(int(nsmpthreads[smp]) + int(npthreads[p]))


def phys_grid_create_chunks_shape_codon(
    nsmpx: int,
    pcols: int,
    chunks_per_thread: int,
    nsmpcolumns_p: cobj,
    nsmpthreads_p: cobj,
    nsmpchunks_p: cobj,
    maxcol_chk_p: cobj,
    maxcol_chks_p: cobj,
) -> int:
    nsmpcolumns = Ptr[i32](nsmpcolumns_p)
    nsmpthreads = Ptr[i32](nsmpthreads_p)
    nsmpchunks = Ptr[i32](nsmpchunks_p)
    maxcol_chk = Ptr[i32](maxcol_chk_p)
    maxcol_chks = Ptr[i32](maxcol_chks_p)

    nchunks = 0
    for smp in range(nsmpx):
        cols = int(nsmpcolumns[smp])
        threads = int(nsmpthreads[smp])
        chunks = cols // pcols
        if cols % pcols != 0:
            chunks += 1
        min_chunks = chunks_per_thread * threads
        if chunks < min_chunks:
            chunks = min_chunks
        while chunks % threads != 0:
            chunks += 1
        if chunks > cols:
            chunks = cols
        nsmpchunks[smp] = i32(chunks)
        nchunks += chunks

    for smp in range(nsmpx):
        chunks = int(nsmpchunks[smp])
        cols = int(nsmpcolumns[smp])
        if chunks != 0:
            ntmp1 = cols // chunks
            ntmp2 = cols % chunks
            if ntmp2 > 0:
                maxcol_chk[smp] = i32(ntmp1 + 1)
                maxcol_chks[smp] = i32(ntmp2)
            else:
                maxcol_chk[smp] = i32(ntmp1)
                maxcol_chks[smp] = i32(chunks)
        else:
            maxcol_chk[smp] = i32(0)
            maxcol_chks[smp] = i32(0)
    return nchunks


def phys_grid_create_chunks_prefix_codon(
    nsmpx: int,
    nsmpchunks_p: cobj,
    cid_offset_p: cobj,
    local_cid_p: cobj,
):
    nsmpchunks = Ptr[i32](nsmpchunks_p)
    cid_offset = Ptr[i32](cid_offset_p)
    local_cid = Ptr[i32](local_cid_p)

    if nsmpx <= 0:
        return
    cid_offset[0] = i32(1)
    local_cid[0] = i32(0)
    for smp in range(1, nsmpx):
        cid_offset[smp] = i32(int(cid_offset[smp - 1]) + int(nsmpchunks[smp - 1]))
        local_cid[smp] = i32(0)


def phys_grid_count_smp_columns_codon(
    nsmpx: int,
    ngcols_p: int,
    latlon_to_dyn_gcol_map_p: cobj,
    col_smp_mapx_p: cobj,
    nsmpcolumns_p: cobj,
):
    latlon_to_dyn_gcol_map = Ptr[i32](latlon_to_dyn_gcol_map_p)
    col_smp_mapx = Ptr[i32](col_smp_mapx_p)
    nsmpcolumns = Ptr[i32](nsmpcolumns_p)

    for smp in range(nsmpx):
        nsmpcolumns[smp] = i32(0)

    for i in range(ngcols_p):
        curgcol = int(latlon_to_dyn_gcol_map[i])
        smp = int(col_smp_mapx[curgcol - 1])
        nsmpcolumns[smp] = i32(int(nsmpcolumns[smp]) + 1)


def phys_grid_zero_int_array_codon(n: int, values_p: cobj):
    values = Ptr[i32](values_p)

    for i in range(n):
        values[i] = i32(0)


def get_gcol_all_p_codon(ncols: int, out_dim: int, src_p: cobj, dst_p: cobj):
    phys_grid_get_gcol_all_codon(ncols, out_dim, src_p, dst_p)


def get_lat_all_p_codon(ncols: int, src_p: cobj, dst_p: cobj):
    phys_grid_get_int_all_codon(ncols, src_p, dst_p)


def get_lon_all_p_codon(
    ncols: int,
    lat_p: cobj,
    gcol_p: cobj,
    map_p: cobj,
    clat_idx_p: cobj,
    dst_p: cobj,
):
    phys_grid_get_lon_all_codon(ncols, lat_p, gcol_p, map_p, clat_idx_p, dst_p)


def get_rlat_all_p_codon(ncols: int, idx_p: cobj, lookup_p: cobj, dst_p: cobj):
    phys_grid_get_lookup_real_all_codon(ncols, idx_p, lookup_p, dst_p)


def get_area_all_p_codon(ncols: int, src_p: cobj, dst_p: cobj):
    phys_grid_get_real_all_codon(ncols, src_p, dst_p)


def get_wght_all_p_codon(ncols: int, src_p: cobj, dst_p: cobj):
    phys_grid_get_real_all_codon(ncols, src_p, dst_p)


def get_rlon_all_p_codon(ncols: int, idx_p: cobj, lookup_p: cobj, dst_p: cobj):
    phys_grid_get_lookup_real_all_codon(ncols, idx_p, lookup_p, dst_p)


def phys_grid_assign_chunks_zero_column_count_codon(
    smp: int,
    nsmpx: int,
    max_nproc_smpx: int,
    ntsks_smpx_p: cobj,
    smp_proc_mapx_p: cobj,
    column_count_p: cobj,
):
    ntsks_smpx = Ptr[i32](ntsks_smpx_p)
    smp_proc_mapx = Ptr[i32](smp_proc_mapx_p)
    column_count = Ptr[i32](column_count_p)

    tasks = int(ntsks_smpx[smp])
    for i in range(1, tasks + 1):
        p = int(smp_proc_mapx[smp + (i - 1) * nsmpx])
        column_count[p] = i32(0)


def phys_grid_assign_chunks_select_owner_codon(
    smp: int,
    nsmpx: int,
    max_nproc_smpx: int,
    ntsks_smpx_p: cobj,
    smp_proc_mapx_p: cobj,
    cur_npchunks_p: cobj,
    npchunks_p: cobj,
    column_count_p: cobj,
) -> int:
    ntsks_smpx = Ptr[i32](ntsks_smpx_p)
    smp_proc_mapx = Ptr[i32](smp_proc_mapx_p)
    cur_npchunks = Ptr[i32](cur_npchunks_p)
    npchunks = Ptr[i32](npchunks_p)
    column_count = Ptr[i32](column_count_p)

    tasks = int(ntsks_smpx[smp])
    for i in range(1, tasks + 1):
        p = int(smp_proc_mapx[smp + (i - 1) * nsmpx])
        if int(cur_npchunks[p]) == int(npchunks[p]):
            column_count[p] = i32(-1)

    best_count = -1
    best_proc = -1
    for i in range(1, tasks + 1):
        p = int(smp_proc_mapx[smp + (i - 1) * nsmpx])
        count = int(column_count[p])
        if count > best_count:
            best_count = count
            best_proc = p

    return best_proc


def phys_grid_assign_chunks_commit_owner_codon(
    owner: int,
    ncols: int,
    cur_npchunks_p: cobj,
    gs_col_num_p: cobj,
):
    cur_npchunks = Ptr[i32](cur_npchunks_p)
    gs_col_num = Ptr[i32](gs_col_num_p)

    cur_npchunks[owner] = i32(int(cur_npchunks[owner]) + 1)
    gs_col_num[owner] = i32(int(gs_col_num[owner]) + ncols)


def phys_grid_assign_chunks_smp_setup_codon(
    npes: int,
    nsmpx: int,
    max_nproc_smpx: int,
    proc_smp_mapx_p: cobj,
    npthreads_p: cobj,
    nsmpthreads_p: cobj,
    nsmpchunks_p: cobj,
    ntsks_smpx_p: cobj,
    smp_proc_mapx_p: cobj,
    cid_offset_p: cobj,
    ntmp1_smp_p: cobj,
    ntmp2_smp_p: cobj,
    ntmp3_smp_p: cobj,
    ntmp4_smp_p: cobj,
    npchunks_p: cobj,
):
    proc_smp_mapx = Ptr[i32](proc_smp_mapx_p)
    npthreads = Ptr[i32](npthreads_p)
    nsmpthreads = Ptr[i32](nsmpthreads_p)
    nsmpchunks = Ptr[i32](nsmpchunks_p)
    ntsks_smpx = Ptr[i32](ntsks_smpx_p)
    smp_proc_mapx = Ptr[i32](smp_proc_mapx_p)
    cid_offset = Ptr[i32](cid_offset_p)
    ntmp1_smp = Ptr[i32](ntmp1_smp_p)
    ntmp2_smp = Ptr[i32](ntmp2_smp_p)
    ntmp3_smp = Ptr[i32](ntmp3_smp_p)
    ntmp4_smp = Ptr[i32](ntmp4_smp_p)
    npchunks = Ptr[i32](npchunks_p)

    for smp in range(nsmpx):
        ntsks_smpx[smp] = i32(0)
    for j in range(max_nproc_smpx):
        for smp in range(nsmpx):
            smp_proc_mapx[smp + j * nsmpx] = i32(-1)

    for p in range(npes):
        smp = int(proc_smp_mapx[p])
        task = int(ntsks_smpx[smp]) + 1
        ntsks_smpx[smp] = i32(task)
        smp_proc_mapx[smp + (task - 1) * nsmpx] = i32(p)

    cid_offset[0] = i32(1)
    for smp in range(1, nsmpx + 1):
        cid_offset[smp] = i32(int(cid_offset[smp - 1]) + int(nsmpchunks[smp - 1]))

    for smp in range(nsmpx):
        ntmp1_smp[smp] = i32(int(nsmpchunks[smp]) // int(nsmpthreads[smp]))
        ntmp2_smp[smp] = i32(int(nsmpchunks[smp]) % int(nsmpthreads[smp]))
        ntmp3_smp[smp] = i32(int(ntmp2_smp[smp]) % int(ntsks_smpx[smp]))
        ntmp4_smp[smp] = i32(int(ntmp2_smp[smp]) // int(ntsks_smpx[smp]))
        if int(ntmp3_smp[smp]) > 0:
            ntmp4_smp[smp] = i32(int(ntmp4_smp[smp]) + 1)

    for p in range(npes):
        smp = int(proc_smp_mapx[p])
        if int(ntmp2_smp[smp]) > int(ntmp4_smp[smp]):
            ntmp2_smp[smp] = i32(int(ntmp2_smp[smp]) - int(ntmp4_smp[smp]))
        else:
            ntmp4_smp[smp] = ntmp2_smp[smp]
            ntmp2_smp[smp] = i32(0)
            ntmp3_smp[smp] = i32(0)

        npchunks[p] = i32(int(ntmp1_smp[smp]) * int(npthreads[p]) + int(ntmp4_smp[smp]))

        if int(ntmp3_smp[smp]) > 0:
            ntmp3_smp[smp] = i32(int(ntmp3_smp[smp]) - 1)
            if int(ntmp3_smp[smp]) == 0:
                ntmp4_smp[smp] = i32(int(ntmp4_smp[smp]) - 1)


def phys_grid_assign_block_no_twin_codon(
    blksiz: int,
    pcols: int,
    smp: int,
    cols_p: cobj,
    cid_offset_p: cobj,
    local_cid_p: cobj,
    nsmpchunks_p: cobj,
    maxcol_chk_p: cobj,
    maxcol_chks_p: cobj,
    dyn_to_latlon_gcol_map_p: cobj,
    lon_p_p: cobj,
    lat_p_p: cobj,
    chunk_ncols_p: cobj,
    chunk_gcol_p: cobj,
    chunk_lon_p: cobj,
    chunk_lat_p: cobj,
    knuhcs_chunkid_p: cobj,
    knuhcs_col_p: cobj,
):
    cols = Ptr[i32](cols_p)
    cid_offset = Ptr[i32](cid_offset_p)
    local_cid = Ptr[i32](local_cid_p)
    nsmpchunks = Ptr[i32](nsmpchunks_p)
    maxcol_chk = Ptr[i32](maxcol_chk_p)
    maxcol_chks = Ptr[i32](maxcol_chks_p)
    dyn_to_latlon_gcol_map = Ptr[i32](dyn_to_latlon_gcol_map_p)
    lon_p = Ptr[i32](lon_p_p)
    lat_p = Ptr[i32](lat_p_p)
    chunk_ncols = Ptr[i32](chunk_ncols_p)
    chunk_gcol = Ptr[i32](chunk_gcol_p)
    chunk_lon = Ptr[i32](chunk_lon_p)
    chunk_lat = Ptr[i32](chunk_lat_p)
    knuhcs_chunkid = Ptr[i32](knuhcs_chunkid_p)
    knuhcs_col = Ptr[i32](knuhcs_col_p)

    for ib in range(blksiz):
        curgcol = int(cols[ib])
        curgcol0 = curgcol - 1
        if int(dyn_to_latlon_gcol_map[curgcol0]) != -1 and int(knuhcs_chunkid[curgcol0]) == -1:
            cid = int(cid_offset[smp]) + int(local_cid[smp])
            cid0 = cid - 1
            if int(maxcol_chks[smp]) > 0:
                while int(chunk_ncols[cid0]) >= int(maxcol_chk[smp]):
                    local_cid[smp] = i32((int(local_cid[smp]) + 1) % int(nsmpchunks[smp]))
                    cid = int(cid_offset[smp]) + int(local_cid[smp])
                    cid0 = cid - 1
            else:
                while int(chunk_ncols[cid0]) >= int(maxcol_chk[smp]) - 1:
                    local_cid[smp] = i32((int(local_cid[smp]) + 1) % int(nsmpchunks[smp]))
                    cid = int(cid_offset[smp]) + int(local_cid[smp])
                    cid0 = cid - 1

            ncols = int(chunk_ncols[cid0]) + 1
            chunk_ncols[cid0] = i32(ncols)
            if ncols == int(maxcol_chk[smp]):
                maxcol_chks[smp] = i32(int(maxcol_chks[smp]) - 1)

            slot = (ncols - 1) + cid0 * pcols
            chunk_gcol[slot] = i32(curgcol)
            chunk_lon[slot] = lon_p[curgcol0]
            chunk_lat[slot] = lat_p[curgcol0]
            knuhcs_chunkid[curgcol0] = i32(cid)
            knuhcs_col[curgcol0] = i32(ncols)

            local_cid[smp] = i32((int(local_cid[smp]) + 1) % int(nsmpchunks[smp]))


def phys_grid_transpose_counts_codon(
    npes: int,
    record_size: int,
    direction: int,
    block_num_p: cobj,
    chunk_num_p: cobj,
    sndcnts_p: cobj,
    sdispls_p: cobj,
    rcvcnts_p: cobj,
    rdispls_p: cobj,
):
    block_num = Ptr[i32](block_num_p)
    chunk_num = Ptr[i32](chunk_num_p)
    sndcnts = Ptr[i32](sndcnts_p)
    sdispls = Ptr[i32](sdispls_p)
    rcvcnts = Ptr[i32](rcvcnts_p)
    rdispls = Ptr[i32](rdispls_p)

    if npes <= 0:
        return

    send_num = block_num
    recv_num = chunk_num
    if direction != 1:
        send_num = chunk_num
        recv_num = block_num

    sdispls[0] = i32(0)
    sndcnts[0] = i32(record_size * int(send_num[0]))
    rdispls[0] = i32(0)
    rcvcnts[0] = i32(record_size * int(recv_num[0]))

    for p in range(1, npes):
        sdispls[p] = i32(int(sdispls[p - 1]) + int(sndcnts[p - 1]))
        sndcnts[p] = i32(record_size * int(send_num[p]))
        rdispls[p] = i32(int(rdispls[p - 1]) + int(rcvcnts[p - 1]))
        rcvcnts[p] = i32(record_size * int(recv_num[p]))


def phys_grid_transpose_lopt_codon(
    phys_alltoall: int,
    max_nproc_smpx: int,
    nproc_busy_d: int,
    npes: int,
    has_window: int,
) -> int:
    if phys_alltoall < 0:
        if max_nproc_smpx > npes // 2 and nproc_busy_d > npes // 2:
            return 0
        return 1

    lopt = phys_alltoall
    if lopt == 2 and has_window == 0:
        lopt = 1
    return lopt
