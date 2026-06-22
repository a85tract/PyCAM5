from math import exp

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

def chemistry_misc_touch_codon(tag: int) -> int:
    return tag

@inline
def _copy_i64_chars(n: int, src: Ptr[int], dst: Ptr[int]):
    for i in range(n):
        dst[i] = src[i]

@inline
def _trimmed_i64_len(n: int, text: Ptr[int]) -> int:
    trimmed = n
    while trimmed > 0 and (text[trimmed - 1] == 32 or text[trimmed - 1] == 0):
        trimmed -= 1
    return trimmed

def prescribed_field_readnl_codon(
    name_len: int,
    name_p: cobj,
    file_len: int,
    file_p: cobj,
    filelist_len: int,
    filelist_p: cobj,
    datapath_len: int,
    datapath_p: cobj,
    type_len: int,
    type_p: cobj,
    rmfile: int,
    cycle_yr: int,
    fixed_ymd: int,
    fixed_tod: int,
    name_out_p: cobj,
    file_out_p: cobj,
    filelist_out_p: cobj,
    datapath_out_p: cobj,
    type_out_p: cobj,
    scalar_out_p: cobj,
):
    name = Ptr[int](name_p)
    filename = Ptr[int](file_p)
    filelist = Ptr[int](filelist_p)
    datapath = Ptr[int](datapath_p)
    data_type = Ptr[int](type_p)
    name_out = Ptr[int](name_out_p)
    filename_out = Ptr[int](file_out_p)
    filelist_out = Ptr[int](filelist_out_p)
    datapath_out = Ptr[int](datapath_out_p)
    data_type_out = Ptr[int](type_out_p)
    scalar_out = Ptr[int](scalar_out_p)

    _copy_i64_chars(name_len, name, name_out)
    _copy_i64_chars(file_len, filename, filename_out)
    _copy_i64_chars(filelist_len, filelist, filelist_out)
    _copy_i64_chars(datapath_len, datapath, datapath_out)
    _copy_i64_chars(type_len, data_type, data_type_out)

    scalar_out[0] = 1 if rmfile != 0 else 0
    scalar_out[1] = cycle_yr
    scalar_out[2] = fixed_ymd
    scalar_out[3] = fixed_tod
    scalar_out[4] = 1 if _trimmed_i64_len(file_len, filename) > 0 else 0

def tracer_defaultopts_codon(
    file_len: int,
    file_p: cobj,
    filelist_len: int,
    filelist_p: cobj,
    datapath_len: int,
    datapath_p: cobj,
    type_len: int,
    type_p: cobj,
    specifier_len: int,
    specifier_count: int,
    specifier_p: cobj,
    rmfile: int,
    cycle_yr: int,
    fixed_ymd: int,
    fixed_tod: int,
    present_file: int,
    file_out_p: cobj,
    present_filelist: int,
    filelist_out_p: cobj,
    present_datapath: int,
    datapath_out_p: cobj,
    present_type: int,
    type_out_p: cobj,
    present_specifier: int,
    specifier_out_p: cobj,
    present_rmfile: int,
    present_cycle_yr: int,
    present_fixed_ymd: int,
    present_fixed_tod: int,
    scalar_out_p: cobj,
):
    if present_file != 0:
        _copy_i64_chars(file_len, Ptr[int](file_p), Ptr[int](file_out_p))

    if present_filelist != 0:
        _copy_i64_chars(filelist_len, Ptr[int](filelist_p), Ptr[int](filelist_out_p))

    if present_datapath != 0:
        _copy_i64_chars(datapath_len, Ptr[int](datapath_p), Ptr[int](datapath_out_p))

    if present_type != 0:
        _copy_i64_chars(type_len, Ptr[int](type_p), Ptr[int](type_out_p))

    if present_specifier != 0:
        specifier = Ptr[int](specifier_p)
        specifier_out = Ptr[int](specifier_out_p)
        total = specifier_len * specifier_count
        for i in range(total):
            specifier_out[i] = specifier[i]

    scalar_out = Ptr[int](scalar_out_p)
    if present_rmfile != 0:
        scalar_out[0] = 1 if rmfile != 0 else 0
    if present_cycle_yr != 0:
        scalar_out[1] = cycle_yr
    if present_fixed_ymd != 0:
        scalar_out[2] = fixed_ymd
    if present_fixed_tod != 0:
        scalar_out[3] = fixed_tod

def linoz_defaultopts_codon(
    file_len: int,
    file_p: cobj,
    filelist_len: int,
    filelist_p: cobj,
    datapath_len: int,
    datapath_p: cobj,
    type_len: int,
    type_p: cobj,
    rmfile: int,
    cycle_yr: int,
    fixed_ymd: int,
    fixed_tod: int,
    present_file: int,
    file_out_p: cobj,
    present_filelist: int,
    filelist_out_p: cobj,
    present_datapath: int,
    datapath_out_p: cobj,
    present_type: int,
    type_out_p: cobj,
    present_rmfile: int,
    present_cycle_yr: int,
    present_fixed_ymd: int,
    present_fixed_tod: int,
    scalar_out_p: cobj,
):
    if present_file != 0:
        _copy_i64_chars(file_len, Ptr[int](file_p), Ptr[int](file_out_p))

    if present_filelist != 0:
        _copy_i64_chars(filelist_len, Ptr[int](filelist_p), Ptr[int](filelist_out_p))

    if present_datapath != 0:
        _copy_i64_chars(datapath_len, Ptr[int](datapath_p), Ptr[int](datapath_out_p))

    if present_type != 0:
        _copy_i64_chars(type_len, Ptr[int](type_p), Ptr[int](type_out_p))

    scalar_out = Ptr[int](scalar_out_p)
    if present_rmfile != 0:
        scalar_out[0] = 1 if rmfile != 0 else 0
    if present_cycle_yr != 0:
        scalar_out[1] = cycle_yr
    if present_fixed_ymd != 0:
        scalar_out[2] = fixed_ymd
    if present_fixed_tod != 0:
        scalar_out[3] = fixed_tod

def spedata_defaultopts_codon(
    spe_data_file_len: int,
    spe_data_file_p: cobj,
    spe_remove_file: int,
    spe_filenames_list_len: int,
    spe_filenames_list_p: cobj,
    present_data_file: int,
    spe_data_file_out_p: cobj,
    present_remove_file: int,
    spe_remove_file_out_p: cobj,
    present_filenames_list: int,
    spe_filenames_list_out_p: cobj,
):
    spe_data_file = Ptr[int](spe_data_file_p)
    spe_filenames_list = Ptr[int](spe_filenames_list_p)
    spe_data_file_out = Ptr[int](spe_data_file_out_p)
    spe_remove_file_out = Ptr[int](spe_remove_file_out_p)
    spe_filenames_list_out = Ptr[int](spe_filenames_list_out_p)

    if present_data_file != 0:
        for i in range(spe_data_file_len):
            spe_data_file_out[i] = spe_data_file[i]

    if present_remove_file != 0:
        spe_remove_file_out[0] = 1 if spe_remove_file != 0 else 0

    if present_filenames_list != 0:
        for i in range(spe_filenames_list_len):
            spe_filenames_list_out[i] = spe_filenames_list[i]


@inline
def _trimmed_equal(name: Ptr[int], name_len: int, items: Ptr[int], item_len: int, item_index: int) -> bool:
    name_trim = name_len
    while name_trim > 0 and name[name_trim - 1] == 32:
        name_trim -= 1

    item_trim = item_len
    item_offset = item_index * item_len
    while item_trim > 0 and items[item_offset + item_trim - 1] == 32:
        item_trim -= 1

    if name_trim != item_trim:
        return False

    for i in range(name_trim):
        if name[i] != items[item_offset + i]:
            return False

    return True


def chem_lookup_name_codon(name_len: int, name_ascii_p: cobj, list_len: int, list_ascii_p: cobj, list_count: int) -> int:
    name_ascii = Ptr[int](name_ascii_p)
    list_ascii = Ptr[int](list_ascii_p)

    for i in range(list_count):
        if _trimmed_equal(name_ascii, name_len, list_ascii, list_len, i):
            return i + 1

    return -1


def has_drydep_codon(name_len: int, name_ascii_p: cobj, list_len: int, list_ascii_p: cobj, list_count: int) -> int:
    return 1 if chem_lookup_name_codon(name_len, name_ascii_p, list_len, list_ascii_p, list_count) > 0 else 0


def chem_lookup_mapped_name_codon(
    name_len: int,
    name_ascii_p: cobj,
    list_len: int,
    list_ascii_p: cobj,
    map_p: cobj,
    list_count: int,
) -> int:
    name_ascii = Ptr[int](name_ascii_p)
    list_ascii = Ptr[int](list_ascii_p)
    item_map = Ptr[int](map_p)

    for i in range(list_count):
        if _trimmed_equal(name_ascii, name_len, list_ascii, list_len, i):
            return item_map[i]

    return -1

def init_mean_mass_ids_codon(lookup_ids_p: cobj, species_ids_p: cobj):
    lookup_ids = Ptr[int](lookup_ids_p)
    species_ids = Ptr[int](species_ids_p)
    for i in range(4):
        species_ids[i] = lookup_ids[i]

def init_hrates_ids_codon(lookup_ids_p: cobj, ptop_ref: float, psurf_ref: float, ids_p: cobj, has_hrates_p: cobj):
    lookup_ids = Ptr[int](lookup_ids_p)
    ids = Ptr[int](ids_p)
    has_hrates = Ptr[int](has_hrates_p)
    all_present = True
    for i in range(9):
        ids[i] = lookup_ids[i]
        if ids[i] <= 0:
            all_present = False
    has_hrates[0] = 1 if all_present and ptop_ref < 0.0004 * psurf_ref else 0

def clybry_fam_init_ids_codon(lookup_ids_p: cobj, ids_p: cobj, has_clybry_p: cobj):
    lookup_ids = Ptr[int](lookup_ids_p)
    ids = Ptr[int](ids_p)
    has_clybry = Ptr[int](has_clybry_p)
    all_present = True
    for i in range(16):
        ids[i] = lookup_ids[i]
        if ids[i] <= 0:
            all_present = False
    has_clybry[0] = 1 if all_present else 0

def init_strato_rates_ids_codon(lookup_ids_p: cobj, ids_p: cobj, has_strato_chem_p: cobj):
    lookup_ids = Ptr[int](lookup_ids_p)
    ids = Ptr[int](ids_p)
    has_strato_chem = Ptr[int](has_strato_chem_p)
    all_present = True
    for i in range(23):
        ids[i] = lookup_ids[i]
        if ids[i] <= 0:
            all_present = False
    has_strato_chem[0] = 1 if all_present else 0

def setinv_inti_ids_codon(lookup_ids_p: cobj, ids_p: cobj, flags_p: cobj):
    lookup_ids = Ptr[int](lookup_ids_p)
    ids = Ptr[int](ids_p)
    flags = Ptr[int](flags_p)
    for i in range(8):
        ids[i] = lookup_ids[i]
    flags[0] = 1 if ids[6] > 0 and ids[5] > 0 and ids[7] > 0 else 0
    flags[1] = 1 if ids[1] > 0 else 0
    flags[2] = 1 if ids[2] > 0 else 0
    flags[3] = 1 if ids[3] > 0 else 0
    flags[4] = 1 if ids[4] > 0 else 0

def gas_wetdep_readnl_status_codon(pcnst: int, list_p: cobj, method_p: cobj, status_p: cobj):
    wetdep_list = Ptr[byte](list_p)
    method = Ptr[byte](method_p)
    status = Ptr[int](status_p)
    count = 0
    for i in range(pcnst):
        has_name = False
        for j in range(8):
            ch = int(wetdep_list[i * 8 + j])
            if ch != 32 and ch != 0:
                has_name = True
        if has_name:
            count += 1
    c0 = int(method[0])
    c1 = int(method[1])
    c2 = int(method[2])
    is_moz = c0 == 77 and c1 == 79 and c2 == 90
    is_neu = c0 == 78 and c1 == 69 and c2 == 85
    status[0] = count
    status[1] = 1 if count > 0 and not (is_moz or is_neu) else 0

def tracer_cnst_init_codon() -> int:
    return 189

@inline
def _has_four_char_prefix(chars: Ptr[int], offset: int, length: int, c0: int, c1: int, c2: int, c3: int) -> bool:
    for j in range(max(0, length - 3)):
        if (
            chars[offset + j] == c0
            and chars[offset + j + 1] == c1
            and chars[offset + j + 2] == c2
            and chars[offset + j + 3] == c3
        ):
            return True
    return False

def rate_diags_init_codon(
    tag_len: int,
    fieldname_len: int,
    tag_count: int,
    tag_ascii_p: cobj,
    rate_name_ascii_p: cobj,
) -> int:
    tag_ascii = Ptr[int](tag_ascii_p)
    rate_name_ascii = Ptr[int](rate_name_ascii_p)

    for item in range(tag_count):
        tag_offset = item * tag_len
        name_offset = item * fieldname_len

        for j in range(fieldname_len):
            rate_name_ascii[name_offset + j] = 32

        start = 0
        if (
            _has_four_char_prefix(tag_ascii, tag_offset, tag_len, 116, 97, 103, 95)
            or _has_four_char_prefix(tag_ascii, tag_offset, tag_len, 117, 115, 114, 95)
            or _has_four_char_prefix(tag_ascii, tag_offset, tag_len, 99, 112, 104, 95)
            or _has_four_char_prefix(tag_ascii, tag_offset, tag_len, 105, 111, 110, 95)
        ):
            start = 4

        tag_end = tag_len
        while tag_end > start and tag_ascii[tag_offset + tag_end - 1] == 32:
            tag_end -= 1

        out_pos = 0
        if fieldname_len > 0:
            rate_name_ascii[name_offset] = 114
            out_pos = 1
        if fieldname_len > 1:
            rate_name_ascii[name_offset + 1] = 95
            out_pos = 2

        tag_pos = start
        while tag_pos < tag_end and out_pos < fieldname_len:
            rate_name_ascii[name_offset + out_pos] = tag_ascii[tag_offset + tag_pos]
            tag_pos += 1
            out_pos += 1

    return tag_count

def noy_ubc_readnl_codon() -> int:
    return 173

def noy_ubc_active_codon(active: int) -> int:
    return active

def spedata_active_codon(active: int) -> int:
    return active

def lightning_inti_active_codon(no_ndx: int, xno_ndx: int) -> int:
    return 1 if no_ndx > 0 or xno_ndx > 0 else 0

def euvac_set_etf_active_codon(active: int) -> int:
    return active

def spe_prod_zero_codon(active: int, ncol: int, pver: int, noxprod_p: cobj, hoxprod_p: cobj) -> int:
    noxprod = Ptr[float](noxprod_p)
    hoxprod = Ptr[float](hoxprod_p)
    for k in range(pver):
        for i in range(ncol):
            idx = i + k * ncol
            noxprod[idx] = 0.0
            hoxprod[idx] = 0.0
    return active

def gcr_ionization_noxhox_zero_codon(active: int, ncol: int, pver: int, gcr_nox_p: cobj, gcr_hox_p: cobj) -> int:
    gcr_nox = Ptr[float](gcr_nox_p)
    gcr_hox = Ptr[float](gcr_hox_p)
    for k in range(pver):
        for i in range(ncol):
            idx = i + k * ncol
            gcr_nox[idx] = 0.0
            gcr_hox[idx] = 0.0
    return active

def airpl_src_active_codon(active: int) -> int:
    return active

def airpl_set_zero_codon(active: int, ncol: int, pver: int, no_air_p: cobj, co_air_p: cobj) -> int:
    no_air = Ptr[float](no_air_p)
    co_air = Ptr[float](co_air_p)
    for k in range(pver):
        for i in range(ncol):
            idx = i + k * ncol
            no_air[idx] = 0.0
            co_air[idx] = 0.0
    return active

def sulf_inti_active_codon(active: int) -> int:
    return active

def fstrat_inti_active_codon(active: int) -> int:
    return active

def set_fstrat_vals_active_codon(active: int) -> int:
    return active

def set_fstrat_h2o_active_codon(active: int) -> int:
    return active

def jeuv_init_active_codon(active: int) -> int:
    return active

def charge_fix_active_codon(active: int) -> int:
    return active

def o1d_to_2oh_adj_init_active_codon(active: int) -> int:
    return active

def init_airglow_active_codon(active: int) -> int:
    return active

def register_cfc11star_active_codon(active: int) -> int:
    return active

def update_cfc11star_active_codon(active: int) -> int:
    return active

def chlorine_loading_init_active_codon(active: int) -> int:
    return active

def parse_rate_sums_active_codon(active: int) -> int:
    return active

def chem_final_codon() -> int:
    return 0

def gcr_ionization_init_codon() -> int:
    return 0

def gcr_ionization_adv_codon() -> int:
    return 0

def aurora_timestep_init_codon() -> int:
    return 0

def set_sulf_time_codon() -> int:
    return 0

def neu_wetdep_init_active_codon(method_is_neu: int, wetdep_count: int) -> int:
    return 1 if method_is_neu != 0 and wetdep_count > 0 else 0

def usrrxt_inti_has_ion_codon(
    ion1: int,
    ion2: int,
    ion3: int,
    elec1: int,
    elec2: int,
    elec3: int,
) -> int:
    return 1 if ion1 > 0 and ion2 > 0 and ion3 > 0 and elec1 > 0 and elec2 > 0 and elec3 > 0 else 0

def comp_exp_codon(x_p: cobj, y_p: cobj, n: int):
    x = Ptr[float](x_p)
    y = Ptr[float](y_p)
    for i in range(n):
        x[i] = exp(y[i])

def heatnirco2_init_xspara_codon(ndpara: int, zppara_p: cobj, xspara_p: cobj):
    zppara = Ptr[float](zppara_p)
    xspara = Ptr[float](xspara_p)
    for k in range(ndpara):
        xspara[k] = 5.0e-7 * exp(-zppara[k])
