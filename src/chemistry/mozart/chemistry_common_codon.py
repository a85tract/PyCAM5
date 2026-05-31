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
