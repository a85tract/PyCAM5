def se_misc_touch_codon(tag: int) -> int:
    return tag


def get_block_gcol_d_codon(size: int, unique_pt_offset: int, cdex_p: cobj):
    cdex = Ptr[i32](cdex_p)
    for ic in range(size):
        cdex[ic] = i32(unique_pt_offset + ic)


def get_block_owner_d_codon(owner: int) -> int:
    return owner


def latlon_interpolation_codon(t: int, n: int, value: int) -> int:
    if t <= n:
        return value
    return 0


def dycore_is_codon(is_match: int) -> int:
    return is_match


def isfactorable_codon(n: int) -> int:
    tmp = n
    while (tmp // 2) * 2 == tmp:
        tmp = tmp // 2
    while (tmp // 3) * 3 == tmp:
        tmp = tmp // 3
    while (tmp // 5) * 5 == tmp:
        tmp = tmp // 5
    if tmp == 1:
        return 1
    return 0


def genlocaldof_codon(ig: int, npts: int, ldof_p: cobj):
    ldof = Ptr[i32](ldof_p)
    npts2 = npts * npts
    for j in range(1, npts + 1):
        for i in range(1, npts + 1):
            ldof[(i - 1) + (j - 1) * npts] = i32((ig - 1) * npts2 + (j - 1) * npts + i)


def uniquepoints2d_codon(num_unique_pts: int, ia_p: cobj, ja_p: cobj, ni: int, src_p: cobj, dest_p: cobj):
    ia = Ptr[i32](ia_p)
    ja = Ptr[i32](ja_p)
    src = Ptr[float](src_p)
    dest = Ptr[float](dest_p)
    for ii in range(1, num_unique_pts + 1):
        i = int(ia[ii - 1])
        j = int(ja[ii - 1])
        dest[ii - 1] = src[(i - 1) + (j - 1) * ni]


def convert_gbl_index_codon(number: int, ne: int, ie_p: cobj, je_p: cobj, face_no_p: cobj):
    ie = Ptr[i32](ie_p)
    je = Ptr[i32](je_p)
    face_no = Ptr[i32](face_no_p)
    face = ((number - 1) // (ne * ne)) + 1
    ie[0] = i32((number - 1) % ne)
    je[0] = i32((number - 1) // ne - (face - 1) * ne)
    face_no[0] = i32(face)


def gridedge_type_codon(head_processor: int, tail_processor: int, internal_edge: int, external_edge: int) -> int:
    if head_processor == tail_processor:
        return internal_edge
    return external_edge


from C import gbarrier_initialize(cobj, i32)
from C import gbarrier_free(cobj)
from C import gbarrier_synchronize(cobj, i32)


def gbarrier_init_codon(c_barrier_p: cobj, nthreads: int):
    gbarrier_initialize(c_barrier_p, i32(nthreads))


def gbarrier_delete_codon(c_barrier_p: cobj):
    gbarrier_free(c_barrier_p)


def gbarrier_synchronize_codon(c_barrier: cobj, thread: int):
    gbarrier_synchronize(c_barrier, i32(thread))
