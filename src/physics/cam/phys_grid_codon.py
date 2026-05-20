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
