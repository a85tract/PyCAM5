def se_misc_touch_codon(tag: int) -> int:
    return tag


def get_block_gcol_d_codon(size: int, unique_pt_offset: int, cdex_p: cobj):
    cdex = Ptr[i32](cdex_p)
    for ic in range(size):
        cdex[ic] = i32(unique_pt_offset + ic)


def get_block_owner_d_codon(owner: int) -> int:
    return owner


def get_dyn_grid_parm_real2d_codon(name_code: int) -> int:
    if name_code == 1:
        return 1
    if name_code == 2:
        return 2
    return 0


def get_dyn_grid_parm_real1d_codon(name_code: int) -> int:
    if name_code == 3:
        return 1
    if name_code == 4:
        return 2
    if name_code == 5:
        return 3
    return 0


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


def copy_buffer_codon(
    nthreads: int,
    ithr: int,
    len_move_ptr: int,
    buf_p: cobj,
    receive_p: cobj,
    move_ptr_p: cobj,
    move_length_p: cobj,
):
    buf = Ptr[float](buf_p)
    receive = Ptr[float](receive_p)
    move_ptr = Ptr[i32](move_ptr_p)
    move_length = Ptr[i32](move_length_p)
    if len_move_ptr == nthreads:
        iptr = int(move_ptr[ithr])
        length = int(move_length[ithr])
        for i in range(0, length):
            receive[iptr + i - 1] = buf[iptr + i - 1]
    elif ithr == 0:
        for j in range(0, len_move_ptr):
            iptr = int(move_ptr[j])
            length = int(move_length[j])
            for i in range(0, length):
                receive[iptr + i - 1] = buf[iptr + i - 1]


def var_is_vector_codon(
    name_len: int,
    name_ascii_p: cobj,
    entry_len: int,
    entries_ascii_p: cobj,
    nentries: int,
) -> int:
    name_ascii = Ptr[int](name_ascii_p)
    entries = Ptr[int](entries_ascii_p)
    trimmed_name_len = name_len
    while trimmed_name_len > 0 and name_ascii[trimmed_name_len - 1] == 32:
        trimmed_name_len -= 1

    for entry in range(1, nentries + 1):
        base = (entry - 1) * entry_len
        trimmed_entry_len = entry_len
        while trimmed_entry_len > 0 and entries[base + trimmed_entry_len - 1] == 32:
            trimmed_entry_len -= 1
        if trimmed_entry_len == 0:
            return 0
        if trimmed_entry_len == trimmed_name_len:
            match = True
            for i in range(0, trimmed_name_len):
                if entries[base + i] != name_ascii[i]:
                    match = False
                    break
            if match:
                return entry
    return 0


def reduction_max_r_local_codon(
    buf_p: cobj,
    ctr_p: cobj,
    redp_p: cobj,
    length: int,
    nthreads: int,
):
    buf = Ptr[float](buf_p)
    ctr = Ptr[i32](ctr_p)
    redp = Ptr[float](redp_p)
    if int(ctr[0]) == 0:
        for k in range(0, length):
            buf[k] = -9.11e30
    if int(ctr[0]) < nthreads:
        for k in range(0, length):
            if redp[k] > buf[k]:
                buf[k] = redp[k]
        ctr[0] += i32(1)
    if int(ctr[0]) == nthreads:
        ctr[0] = i32(0)


def reduction_min_r_local_codon(
    buf_p: cobj,
    ctr_p: cobj,
    redp_p: cobj,
    length: int,
    nthreads: int,
):
    buf = Ptr[float](buf_p)
    ctr = Ptr[i32](ctr_p)
    redp = Ptr[float](redp_p)
    if int(ctr[0]) == 0:
        for k in range(0, length):
            buf[k] = 9.11e30
    if int(ctr[0]) < nthreads:
        for k in range(0, length):
            if redp[k] < buf[k]:
                buf[k] = redp[k]
        ctr[0] += i32(1)
    if int(ctr[0]) == nthreads:
        ctr[0] = i32(0)


def copy_par_codon(
    rank2_p: cobj,
    root2_p: cobj,
    nprocs2_p: cobj,
    comm2_p: cobj,
    intercomm2_p: cobj,
    intracomm2_p: cobj,
    intracommsize2_p: cobj,
    intracommrank2_p: cobj,
    comm_graph_full2_p: cobj,
    comm_graph_inter2_p: cobj,
    comm_graph_intra2_p: cobj,
    group_graph_full2_p: cobj,
    masterproc2_p: cobj,
    rank1_p: cobj,
    root1_p: cobj,
    nprocs1_p: cobj,
    comm1_p: cobj,
    intercomm1_p: cobj,
    intracomm1_p: cobj,
    intracommsize1_p: cobj,
    intracommrank1_p: cobj,
    comm_graph_full1_p: cobj,
    comm_graph_inter1_p: cobj,
    comm_graph_intra1_p: cobj,
    group_graph_full1_p: cobj,
    masterproc1_p: cobj,
):
    Ptr[i32](rank2_p)[0] = Ptr[i32](rank1_p)[0]
    Ptr[i32](root2_p)[0] = Ptr[i32](root1_p)[0]
    Ptr[i32](nprocs2_p)[0] = Ptr[i32](nprocs1_p)[0]
    Ptr[i32](comm2_p)[0] = Ptr[i32](comm1_p)[0]
    Ptr[i32](intercomm2_p)[0] = Ptr[i32](intercomm1_p)[0]
    Ptr[i32](intracomm2_p)[0] = Ptr[i32](intracomm1_p)[0]
    Ptr[i32](intracommsize2_p)[0] = Ptr[i32](intracommsize1_p)[0]
    Ptr[i32](intracommrank2_p)[0] = Ptr[i32](intracommrank1_p)[0]
    Ptr[i32](comm_graph_full2_p)[0] = Ptr[i32](comm_graph_full1_p)[0]
    Ptr[i32](comm_graph_inter2_p)[0] = Ptr[i32](comm_graph_inter1_p)[0]
    Ptr[i32](comm_graph_intra2_p)[0] = Ptr[i32](comm_graph_intra1_p)[0]
    Ptr[i32](group_graph_full2_p)[0] = Ptr[i32](group_graph_full1_p)[0]
    Ptr[i32](masterproc2_p)[0] = Ptr[i32](masterproc1_p)[0]


def init_edge_buffer_i8_header_codon(
    np: int,
    max_corner_elem: int,
    nelemd: int,
    nlyr: int,
    nlyr_p: cobj,
    nbuf_p: cobj,
):
    Ptr[i32](nlyr_p)[0] = i32(nlyr)
    Ptr[i32](nbuf_p)[0] = i32(4 * (np + max_corner_elem) * nelemd)


def zero_i32_buffer_codon(n: int, buf_p: cobj):
    buf = Ptr[i32](buf_p)
    for i in range(0, n):
        buf[i] = i32(0)


from C import gbarrier_initialize(cobj, i32)
from C import gbarrier_free(cobj)
from C import gbarrier_synchronize(cobj, i32)


def gbarrier_init_codon(c_barrier_p: cobj, nthreads: int):
    gbarrier_initialize(c_barrier_p, i32(nthreads))


def gbarrier_delete_codon(c_barrier_p: cobj):
    gbarrier_free(c_barrier_p)


def gbarrier_synchronize_codon(c_barrier: cobj, thread: int):
    gbarrier_synchronize(c_barrier, i32(thread))


def se_log2_codon(n: int) -> int:
    ans = 1
    tmp = n
    while tmp // 2 != 1:
        tmp = tmp // 2
        ans += 1
    return ans


def se_calcsegmentlength_codon(lenp: int, lens: int, mpattern: int, nlyr: int, hme_mpattern_s: int, hme_mpattern_p: int) -> int:
    if mpattern == hme_mpattern_s:
        ans = nlyr * lens
    elif mpattern == hme_mpattern_p:
        ans = nlyr * lenp
    else:
        ans = nlyr * lenp
    return ans


def se_timelevel_init_default_codon(nm1_p: cobj, n0_p: cobj, np1_p: cobj, nstep_p: cobj, nstep0_p: cobj):
    Ptr[i32](nm1_p)[0] = i32(1)
    Ptr[i32](n0_p)[0] = i32(2)
    Ptr[i32](np1_p)[0] = i32(3)
    Ptr[i32](nstep_p)[0] = i32(0)
    Ptr[i32](nstep0_p)[0] = i32(2)


def se_timelevel_update_codon(nm1_p: cobj, n0_p: cobj, np1_p: cobj, nstep_p: cobj, uptype_code: int) -> int:
    nm1 = Ptr[i32](nm1_p)
    n0 = Ptr[i32](n0_p)
    np1 = Ptr[i32](np1_p)
    nstep = Ptr[i32](nstep_p)
    if uptype_code == 1:
        ntmp = int(np1[0])
        np1[0] = nm1[0]
        nm1[0] = n0[0]
        n0[0] = i32(ntmp)
    elif uptype_code == 2:
        ntmp = int(np1[0])
        np1[0] = n0[0]
        n0[0] = i32(ntmp)
    else:
        return 1
    nstep[0] = i32(int(nstep[0]) + 1)
    return 0
