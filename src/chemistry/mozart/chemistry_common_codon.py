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
