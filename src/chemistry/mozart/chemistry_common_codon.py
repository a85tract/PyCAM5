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
