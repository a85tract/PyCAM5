from C import MPI_Alltoall(Ptr[byte], i32, i32, Ptr[byte], i32, i32, i32) -> i32
from C import MPI_Allgather(Ptr[byte], i32, i32, Ptr[byte], i32, i32, i32) -> i32


@export
def cam_misc_touch_codon(tag: int) -> int:
    return tag


@export
def get_hist_restart_filepath_codon(tag: int) -> int:
    return tag


@export
def get_ptapes_codon(tag: int) -> int:
    return tag


@export
def init_masterlinkedlist_codon(tag: int) -> int:
    return tag


@export
def hist_fld_active_codon(tag: int) -> int:
    return tag


@export
def add_entry_to_master_codon(tag: int) -> int:
    return tag


@export
def get_entry_by_name_codon(tag: int) -> int:
    return tag


@export
def restartvar_getdesc_codon(tag: int) -> int:
    return tag


@export
def restart_dims_setnames_codon(tag: int) -> int:
    return tag


@export
def atm_final_mct_codon(tag: int) -> int:
    return tag


@export
def atm_setgsmap_mct_codon(tag: int) -> int:
    return tag


@export
def atm_write_srfrest_mct_codon(tag: int) -> int:
    return tag


@export
def getunit_codon(tag: int) -> int:
    return tag


@export
def freeunit_codon(tag: int) -> int:
    return tag


@export
def restart_defaultopts_codon(tag: int) -> int:
    return tag


@export
def restart_printopts_codon(tag: int) -> int:
    return tag


@export
def write_rest_pfile_codon(tag: int) -> int:
    return tag


@export
def init_pio_subsystem_codon(tag: int) -> int:
    return tag


@export
def cam_pio_openfile_codon(tag: int) -> int:
    return tag


@export
def cam_pio_createfile_codon(tag: int) -> int:
    return tag


@export
def get_phys_decomp_md1d_codon(tag: int) -> int:
    return tag


@export
def clean_iodesc_list_codon(tag: int) -> int:
    return tag


@export
def find_iodesc_codon(tag: int) -> int:
    return tag


@export
def get_decomp_codon(tag: int) -> int:
    return tag


@export
def write_restart_hycoef_codon(tag: int) -> int:
    return tag


@export
def init_restart_hycoef_codon(tag: int) -> int:
    return tag


@export
def hycoef_read_codon(tag: int) -> int:
    return tag


@export
def get_hist_coord_index_codon(tag: int) -> int:
    return tag


@export
def add_vert_coord_codon(tag: int) -> int:
    return tag


@export
def write_hist_coord_vars_codon(tag: int) -> int:
    return tag


@export
def write_hist_coord_attrs_codon(tag: int) -> int:
    return tag


@export
def initial_conds_codon(tag: int) -> int:
    return tag


@export
def allocate_coords_codon(tag: int) -> int:
    return tag


@export
def finalize_codon(tag: int) -> int:
    return tag


@export
def get_curr_time_codon(tag: int) -> int:
    return tag


@export
def timevars_set_names_codon(tag: int) -> int:
    return tag


@export
def timemgr_init_codon(tag: int) -> int:
    return tag


@export
def get_curr_date_codon(tag: int) -> int:
    return tag


@export
def set_time_float_from_date_codon(tag: int) -> int:
    return tag


@export
def get_curr_calday_codon(tag: int) -> int:
    return tag


@export
def initialize_clock_codon(tag: int) -> int:
    return tag


@export
def bldfld_codon(tag: int) -> int:
    return tag


@export
def opnfil_codon(tag: int) -> int:
    return tag


@export
def hub2atm_deallocate_codon(tag: int) -> int:
    return tag


@export
def spmd_utils_readnl_codon(tag: int) -> int:
    return tag


@inline
def _cam_char_is_none(chars: Ptr[int], n: int) -> int:
    if n < 4:
        return 0
    if chars[0] != 78 or chars[1] != 79 or chars[2] != 78 or chars[3] != 69:
        return 0
    for i in range(4, n):
        if chars[i] != 32:
            return 0
    return 1


@inline
def _cam_blank_fill(out: Ptr[int], n: int):
    for i in range(n):
        out[i] = 32


@export
def glc_codon(chars_p: cobj, n: int) -> int:
    chars = Ptr[int](chars_p)
    if n == 0:
        return 0
    for i in range(n - 1, -1, -1):
        if chars[i] != 32 and chars[i] != 0:
            return i + 1
    return 0


@export
def to_upper_codon(in_p: cobj, out_p: cobj, n: int):
    inp = Ptr[int](in_p)
    out = Ptr[int](out_p)
    for i in range(n):
        ch = inp[i]
        if ch >= 97 and ch <= 122:
            ch -= 32
        out[i] = ch


@export
def to_lower_codon(in_p: cobj, out_p: cobj, n: int):
    inp = Ptr[int](in_p)
    out = Ptr[int](out_p)
    for i in range(n):
        ch = inp[i]
        if ch >= 65 and ch <= 90:
            ch += 32
        out[i] = ch


@export
def strip_suffix_codon(name_p: cobj, name_len: int, out_p: cobj, out_len: int, field_len: int):
    name = Ptr[int](name_p)
    out = Ptr[int](out_p)
    _cam_blank_fill(out, out_len)
    limit = field_len
    if limit > out_len:
        limit = out_len

    for i in range(limit):
        out[i] = name[i] if i < name_len else 32
        next_i = i + 1
        if next_i >= name_len:
            return
        if name[next_i] == 32:
            return
        if next_i + 2 < name_len and name[next_i] == 38 and name[next_i + 1] == 73 and name[next_i + 2] == 67:
            return

    for i in range(field_len, out_len):
        out[i] = name[i] if i < name_len else 32


@export
def is_initfile_codon(inithist_p: cobj, inithist_len: int, has_file: int, file_index: int, ptapes: int) -> int:
    inithist = Ptr[int](inithist_p)
    if _cam_char_is_none(inithist, inithist_len) != 0:
        return 0
    if has_file != 0:
        return 1 if file_index == ptapes else 0
    return 1


@export
def gen_hash_key_codon(chars_p: cobj, n: int) -> int:
    chars = Ptr[int](chars_p)
    keys = (61, 59, 53, 47, 43, 41, 37, 31, 29, 23, 17, 13, 11, 7, 3, 1)
    hash_value = 0x000053DB
    if n != 19:
        for i in range(n):
            hash_value = hash_value ^ (chars[i] * keys[i & 15])
    else:
        hash_value = hash_value ^ (chars[0] * 61)
        hash_value = hash_value ^ (chars[1] * 59)
        hash_value = hash_value ^ (chars[2] * 53)
        hash_value = hash_value ^ (chars[3] * 47)
        hash_value = hash_value ^ (chars[4] * 43)
        hash_value = hash_value ^ (chars[5] * 41)
        hash_value = hash_value ^ (chars[6] * 37)
        hash_value = hash_value ^ (chars[7] * 31)
        hash_value = hash_value ^ (chars[8] * 29)
        hash_value = hash_value ^ (chars[9] * 23)
        hash_value = hash_value ^ (chars[10] * 17)
        hash_value = hash_value ^ (chars[11] * 13)
        hash_value = hash_value ^ (chars[12] * 11)
        hash_value = hash_value ^ (chars[13] * 7)
        hash_value = hash_value ^ (chars[14] * 3)
        hash_value = hash_value ^ chars[15]
        hash_value = hash_value ^ (chars[16] * 61)
        hash_value = hash_value ^ (chars[17] * 59)
        hash_value = hash_value ^ (chars[18] * 53)
    return hash_value & 65535


@inline
def _cam_list_name_is_blank(list_chars: Ptr[int], offset: int, field_len: int) -> int:
    for i in range(field_len):
        ch = list_chars[offset + i]
        if ch == 58:
            return 1
        if ch != 32:
            return 0
    return 1


@inline
def _cam_list_name_matches(list_chars: Ptr[int], offset: int, name: Ptr[int], name_len: int, field_len: int) -> int:
    for i in range(field_len):
        list_ch = list_chars[offset + i]
        if list_ch == 58:
            list_ch = 32
        name_ch = name[i] if i < name_len else 32
        if list_ch != name_ch:
            return 0
        if list_chars[offset + i] == 58:
            for j in range(i + 1, field_len):
                name_tail = name[j] if j < name_len else 32
                if name_tail != 32:
                    return 0
            return 1
    for i in range(field_len, name_len):
        if name[i] != 32:
            return 0
    return 1


@export
def list_index_codon(list_p: cobj, list_len: int, name_p: cobj, name_len: int, pflds: int, field_len: int) -> int:
    list_chars = Ptr[int](list_p)
    name = Ptr[int](name_p)
    for f in range(pflds):
        offset = f * list_len
        if _cam_list_name_is_blank(list_chars, offset, field_len) != 0:
            return 0
        if _cam_list_name_matches(list_chars, offset, name, name_len, field_len) != 0:
            return f + 1
    return 0


@export
def date2yyyymmdd_codon(date: int, out_p: cobj) -> int:
    if date < 0:
        return 0
    out = Ptr[int](out_p)
    year = date // 10000
    month = (date - year * 10000) // 100
    day = date - year * 10000 - month * 100
    out[0] = 48 + (year // 1000) % 10
    out[1] = 48 + (year // 100) % 10
    out[2] = 48 + (year // 10) % 10
    out[3] = 48 + year % 10
    out[4] = 45
    out[5] = 48 + (month // 10) % 10
    out[6] = 48 + month % 10
    out[7] = 45
    out[8] = 48 + (day // 10) % 10
    out[9] = 48 + day % 10
    return 1


@inline
def _write_two_digits(out: Ptr[int], off: int, value: int):
    out[off] = 48 + value // 10
    out[off + 1] = 48 + value % 10


@export
def datetime_codon(values_p: cobj, cdate_p: cobj, ctime_p: cobj):
    values = Ptr[int](values_p)
    cdate = Ptr[int](cdate_p)
    ctime = Ptr[int](ctime_p)

    year = values[0]
    month = values[1]
    day = values[2]
    hour = values[4]
    minute = values[5]
    second = values[6]

    _write_two_digits(cdate, 0, month)
    cdate[2] = 47
    _write_two_digits(cdate, 3, day)
    cdate[5] = 47
    _write_two_digits(cdate, 6, year % 100)

    _write_two_digits(ctime, 0, hour)
    ctime[2] = 58
    _write_two_digits(ctime, 3, minute)
    ctime[5] = 58
    _write_two_digits(ctime, 6, second)


@export
def mpialltoallint_codon(
    sendbuf_p: cobj,
    sendcnt: int,
    recvbuf_p: cobj,
    recvcnt: int,
    comm: int,
    mpiint_arg: int,
) -> int:
    return int(
        MPI_Alltoall(
            Ptr[byte](sendbuf_p),
            i32(sendcnt),
            i32(mpiint_arg),
            Ptr[byte](recvbuf_p),
            i32(recvcnt),
            i32(mpiint_arg),
            i32(comm),
        )
    )


@export
def mpiallgatherint_codon(
    sendbuf_p: cobj,
    scount: int,
    recvbuf_p: cobj,
    rcount: int,
    comm: int,
    mpiint_arg: int,
) -> int:
    return int(
        MPI_Allgather(
            Ptr[byte](sendbuf_p),
            i32(scount),
            i32(mpiint_arg),
            Ptr[byte](recvbuf_p),
            i32(rcount),
            i32(mpiint_arg),
            i32(comm),
        )
    )


@export
def atm2hub_deallocate_codon(tag: int) -> int:
    return tag


@export
def is_satfile_codon(file_index: int, sat_tape_num: int) -> int:
    return 1 if file_index == sat_tape_num else 0


@export
def findplb_codon(x: Ptr[float], nx: int, xval: float) -> int:
    if xval < x[0] or xval >= x[nx - 1]:
        return nx

    for i in range(1, nx):
        if xval < x[i]:
            return i

    return nx


@export
def lininterp_full1d_codon(
    arrin_p: cobj,
    yin_p: cobj,
    yout_p: cobj,
    arrout_p: cobj,
    nin: int,
    nout: int,
) -> int:
    arrin = Ptr[float](arrin_p)
    yin = Ptr[float](yin_p)
    yout = Ptr[float](yout_p)
    arrout = Ptr[float](arrout_p)

    if nin < 2:
        return 1

    icount = 0
    for j in range(nin - 1):
        if yin[j] > yin[j + 1]:
            icount += 1

    increasing = True
    if icount == nin - 1:
        increasing = False
        icount = 0

    if icount > 0:
        return 2

    for j in range(nout):
        jjm = 0
        jjp = 0
        wgts = 0.0
        wgtn = 0.0

        if increasing:
            if yout[j] <= yin[0]:
                jjm = 1
                jjp = 1
                wgts = 1.0
                wgtn = 0.0
            elif yout[j] > yin[nin - 1]:
                jjm = nin
                jjp = nin
                wgts = 1.0
                wgtn = 0.0
        else:
            if yout[j] > yin[0]:
                jjm = 1
                jjp = 1
                wgts = 1.0
                wgtn = 0.0
            elif yout[j] <= yin[nin - 1]:
                jjm = nin
                jjp = nin
                wgts = 1.0
                wgtn = 0.0

        if increasing:
            for jj in range(nin - 1):
                if yout[j] > yin[jj] and yout[j] <= yin[jj + 1]:
                    jjm = jj + 1
                    jjp = jj + 2
                    wgts = (yin[jj + 1] - yout[j]) / (yin[jj + 1] - yin[jj])
                    wgtn = (yout[j] - yin[jj]) / (yin[jj + 1] - yin[jj])
                    break
        else:
            for jj in range(nin - 1):
                if yout[j] <= yin[jj] and yout[j] > yin[jj + 1]:
                    jjm = jj + 1
                    jjp = jj + 2
                    wgts = (yin[jj + 1] - yout[j]) / (yin[jj + 1] - yin[jj])
                    wgtn = (yout[j] - yin[jj]) / (yin[jj + 1] - yin[jj])
                    break

        if jjm == 0 or jjp == 0:
            return 4

        ratio = wgts + wgtn
        if ratio < 0.9 or ratio > 1.1:
            return 3

        arrout[j] = arrin[jjm - 1] * wgts + arrin[jjp - 1] * wgtn

    return 0


@export
def lininterp1d_codon(
    arrin_p: cobj,
    arrout_p: cobj,
    wgts_p: cobj,
    wgtn_p: cobj,
    jjm_p: cobj,
    jjp_p: cobj,
    m1: int,
):
    arrin = Ptr[float](arrin_p)
    arrout = Ptr[float](arrout_p)
    wgts = Ptr[float](wgts_p)
    wgtn = Ptr[float](wgtn_p)
    jjm = Ptr[int](jjm_p)
    jjp = Ptr[int](jjp_p)

    for j in range(m1):
        arrout[j] = arrin[jjm[j] - 1] * wgts[j] + arrin[jjp[j] - 1] * wgtn[j]


@export
def vertinterp_codon(
    ncol: int,
    ncold: int,
    nlev: int,
    pmid_p: cobj,
    pout: float,
    arrin_p: cobj,
    arrout_p: cobj,
) -> int:
    pmid = Ptr[float](pmid_p)
    arrin = Ptr[float](arrin_p)
    arrout = Ptr[float](arrout_p)

    for i in range(ncol):
        found = False
        kupper = 1
        for k in range(nlev - 1):
            if (not found) and pmid[i + ncold * k] < pout and pout <= pmid[i + ncold * (k + 1)]:
                found = True
                kupper = k + 1

        if pout <= pmid[i]:
            arrout[i] = arrin[i]
        elif pout >= pmid[i + ncold * (nlev - 1)]:
            arrout[i] = arrin[i + ncold * (nlev - 1)]
        elif found:
            k0 = kupper - 1
            dpu = pout - pmid[i + ncold * k0]
            dpl = pmid[i + ncold * (k0 + 1)] - pout
            arrout[i] = (
                arrin[i + ncold * k0] * dpl + arrin[i + ncold * (k0 + 1)] * dpu
            ) / (dpl + dpu)
        else:
            return 1

    return 0


@export
def lininterp2d1d_codon(
    arrin_p: cobj,
    arrout_p: cobj,
    wgtw_p: cobj,
    wgte_p: cobj,
    wgts_p: cobj,
    wgtn_p: cobj,
    iim_p: cobj,
    iip_p: cobj,
    jjm_p: cobj,
    jjp_p: cobj,
    n1: int,
    m1: int,
):
    arrin = Ptr[float](arrin_p)
    arrout = Ptr[float](arrout_p)
    wgtw = Ptr[float](wgtw_p)
    wgte = Ptr[float](wgte_p)
    wgts = Ptr[float](wgts_p)
    wgtn = Ptr[float](wgtn_p)
    iim = Ptr[int](iim_p)
    iip = Ptr[int](iip_p)
    jjm = Ptr[int](jjm_p)
    jjp = Ptr[int](jjp_p)

    for i in range(m1):
        ii_m = iim[i] - 1
        ii_p = iip[i] - 1
        jj_m = jjm[i] - 1
        jj_p = jjp[i] - 1
        arrout[i] = (
            arrin[ii_m + n1 * jj_m] * wgtw[i] * wgts[i]
            + arrin[ii_p + n1 * jj_m] * wgte[i] * wgts[i]
            + arrin[ii_m + n1 * jj_p] * wgtw[i] * wgtn[i]
            + arrin[ii_p + n1 * jj_p] * wgte[i] * wgtn[i]
        )


@export
def lininterp3d2d_codon(
    arrin_p: cobj,
    arrout_p: cobj,
    wgtw_p: cobj,
    wgte_p: cobj,
    wgts_p: cobj,
    wgtn_p: cobj,
    iim_p: cobj,
    iip_p: cobj,
    jjm_p: cobj,
    jjp_p: cobj,
    n1: int,
    n2: int,
    len1: int,
    n3: int,
    m1: int,
):
    arrin = Ptr[float](arrin_p)
    arrout = Ptr[float](arrout_p)
    wgtw = Ptr[float](wgtw_p)
    wgte = Ptr[float](wgte_p)
    wgts = Ptr[float](wgts_p)
    wgtn = Ptr[float](wgtn_p)
    iim = Ptr[int](iim_p)
    iip = Ptr[int](iip_p)
    jjm = Ptr[int](jjm_p)
    jjp = Ptr[int](jjp_p)

    for k in range(n3):
        arrin_k = n1 * n2 * k
        arrout_k = len1 * k
        for i in range(m1):
            ii_m = iim[i] - 1
            ii_p = iip[i] - 1
            jj_m = jjm[i] - 1
            jj_p = jjp[i] - 1
            arrout[arrout_k + i] = (
                arrin[ii_m + n1 * jj_m + arrin_k] * wgtw[i] * wgts[i]
                + arrin[ii_p + n1 * jj_m + arrin_k] * wgte[i] * wgts[i]
                + arrin[ii_m + n1 * jj_p + arrin_k] * wgtw[i] * wgtn[i]
                + arrin[ii_p + n1 * jj_p + arrin_k] * wgte[i] * wgtn[i]
            )


@inline
def _hbuf_idx(i: int, k: int, ld: int) -> int:
    return i + k * ld


@inline
def _field_idx(i: int, k: int, idim: int) -> int:
    return i + k * idim


@export
def hbuf_accum_inst_codon(
    buf8_p: cobj,
    field_p: cobj,
    nacs_p: cobj,
    buf8_ld: int,
    idim: int,
    ieu: int,
    jeu: int,
    flag_xyfill: int,
    fillvalue: float,
):
    buf8 = Ptr[float](buf8_p)
    field = Ptr[float](field_p)
    nacs = Ptr[int](nacs_p)

    for k in range(jeu):
        for i in range(ieu):
            buf8[_hbuf_idx(i, k, buf8_ld)] = field[_field_idx(i, k, idim)]

    if flag_xyfill != 0:
        for i in range(ieu):
            if field[_field_idx(i, 0, idim)] == fillvalue:
                nacs[i] = 0
            else:
                nacs[i] = 1
    else:
        nacs[0] = 1


@export
def handle_pio_error_codon(ierr: int, pio_noerr: int) -> int:
    return 1 if ierr == pio_noerr else 0


@export
def cam_initfile_getter_touch_codon(tag: int) -> int:
    return tag


@export
def initial_file_get_id_codon(tag: int) -> int:
    return cam_initfile_getter_touch_codon(tag)


@export
def topo_file_get_id_codon(tag: int) -> int:
    return cam_initfile_getter_touch_codon(tag)


@export
def getname_codon(tag: int) -> int:
    return tag


@export
def sec2hms_codon(tag: int) -> int:
    return tag


@export
def write_hist_coord_att_codon(tag: int) -> int:
    return tag


@export
def read_initial_codon(tag: int) -> int:
    return tag


@export
def preset_codon(tag: int) -> int:
    return tag


@export
def cam_initfiles_close_codon(tag: int) -> int:
    return tag


@export
def cam_initfiles_open_codon(tag: int) -> int:
    return tag


@export
def restart_setopts_codon(tag: int) -> int:
    return tag


@export
def write_inithist_codon(tag: int) -> int:
    return tag


@export
def lookup_hist_coord_indices_codon(tag: int) -> int:
    return tag


@export
def check_var_codon(tag: int) -> int:
    return tag


@export
def scam_default_opts_codon(tag: int) -> int:
    return tag


@export
def formula_terms_copy_codon(tag: int) -> int:
    return tag


@export
def get_masterlist_indx_codon(tag: int) -> int:
    return tag


@export
def bld_htapefld_indices_codon(tag: int) -> int:
    return tag


@export
def inifld_codon(tag: int) -> int:
    return tag


@export
def setup_interpolation_and_define_vector_compliments_codon(tag: int) -> int:
    return tag


@export
def add_hist_coord_regonly_codon(tag: int) -> int:
    return tag


@export
def add_default_codon(tag: int) -> int:
    return tag


@export
def init_restart_history_codon(tag: int) -> int:
    return tag


@export
def sat_hist_init_codon(has_sat_hist: int) -> int:
    return 1 if has_sat_hist == 0 else 0


@export
def scalar_add_tridiag_codon(diag: Ptr[float], nsys: int, ncel: int, constant: float):
    for i in range(nsys * ncel):
        diag[i] += constant


@export
def new_boundaryfixedlayer_codon(
    bndry_type: Ptr[int],
    edge_width: Ptr[float],
    width: Ptr[float],
    n: int,
):
    bndry_type[0] = 3
    for i in range(n):
        edge_width[i] = width[i]


@export
def advance_timestep_codon() -> int:
    return 0


@export
def is_end_curr_day_codon(tod: int) -> int:
    return 1 if tod == 0 else 0


@export
def is_first_restart_step_codon(first_restart_step: int) -> int:
    return 1 if first_restart_step != 0 else 0


@export
def chkrc_codon(rc: int, success: int) -> int:
    return 1 if rc == success else 0


@export
def alloc_err_codon(istat: int) -> int:
    return 1 if istat != 0 else 0


@export
def handle_err_codon(istat: int) -> int:
    return 1 if istat != 0 else 0


@export
def handle_ncerr_codon(ret: int, noerr: int) -> int:
    return 1 if ret != noerr else 0


@export
def pair_codon(np: int, p: int, k: int) -> int:
    q = p ^ k
    return -1 if q > np - 1 else q


@export
def get_step_size_codon(step_seconds: int) -> int:
    return step_seconds


@export
def get_nstep_codon(step_no: int) -> int:
    return step_no


@export
def get_ref_date_codon(yr: int, mon: int, day: int, tod: int, out_p: Ptr[int]):
    out_p[0] = yr
    out_p[1] = mon
    out_p[2] = day
    out_p[3] = tod


@export
def is_first_step_codon(step_no: int) -> int:
    return 1 if step_no == 0 else 0


@export
def get_calday_codon(day_of_year: float, is_gregorian: int) -> float:
    if day_of_year > 366.0 and day_of_year <= 367.0 and is_gregorian != 0:
        return day_of_year - 1.0
    return day_of_year


@export
def init_calendar_codon(cal_code: int) -> int:
    return cal_code


@inline
def _cam_char_trim_len(chars: Ptr[byte], n: int) -> int:
    last = n
    while last > 0:
        ch = int(chars[last - 1])
        if ch != 32 and ch != 0:
            break
        last -= 1
    return last


@inline
def _cam_ascii_upper(ch: int) -> int:
    if ch >= 97 and ch <= 122:
        return ch - 32
    return ch


@inline
def _cam_trim_equals_noleap(chars: Ptr[byte], trim_len: int) -> int:
    if trim_len != 7:
        return 0
    if _cam_ascii_upper(int(chars[0])) != 78:
        return 0
    if _cam_ascii_upper(int(chars[1])) != 79:
        return 0
    if int(chars[2]) != 95:
        return 0
    if _cam_ascii_upper(int(chars[3])) != 76:
        return 0
    if _cam_ascii_upper(int(chars[4])) != 69:
        return 0
    if _cam_ascii_upper(int(chars[5])) != 65:
        return 0
    if _cam_ascii_upper(int(chars[6])) != 80:
        return 0
    return 1


@inline
def _cam_trim_equals_gregorian(chars: Ptr[byte], trim_len: int) -> int:
    if trim_len != 9:
        return 0
    if _cam_ascii_upper(int(chars[0])) != 71:
        return 0
    if _cam_ascii_upper(int(chars[1])) != 82:
        return 0
    if _cam_ascii_upper(int(chars[2])) != 69:
        return 0
    if _cam_ascii_upper(int(chars[3])) != 71:
        return 0
    if _cam_ascii_upper(int(chars[4])) != 79:
        return 0
    if _cam_ascii_upper(int(chars[5])) != 82:
        return 0
    if _cam_ascii_upper(int(chars[6])) != 73:
        return 0
    if _cam_ascii_upper(int(chars[7])) != 65:
        return 0
    if _cam_ascii_upper(int(chars[8])) != 78:
        return 0
    return 1


@export
def timemgr_get_calendar_cf_codon(calendar_p: cobj, calendar_len: int) -> int:
    calendar_chars = Ptr[byte](calendar_p)
    calendar_trim = _cam_char_trim_len(calendar_chars, calendar_len)

    if _cam_trim_equals_noleap(calendar_chars, calendar_trim) != 0:
        return 1
    if _cam_trim_equals_gregorian(calendar_chars, calendar_trim) != 0:
        return 2
    return 0


@export
def timemgr_is_caltype_codon(calendar_p: cobj, calendar_len: int, cal_in_p: cobj, cal_in_len: int) -> int:
    calendar_chars = Ptr[byte](calendar_p)
    cal_in_chars = Ptr[byte](cal_in_p)

    calendar_trim = _cam_char_trim_len(calendar_chars, calendar_len)
    cal_in_trim = _cam_char_trim_len(cal_in_chars, cal_in_len)
    if calendar_trim != cal_in_trim:
        return 0

    for i in range(calendar_trim):
        if _cam_ascii_upper(int(calendar_chars[i])) != _cam_ascii_upper(int(cal_in_chars[i])):
            return 0
    return 1


@export
def timemgr_init_restart_codon(varcnt: int) -> int:
    return varcnt


@export
def timesetymd_codon(ymd: int, tod: int, out_p: Ptr[int]) -> int:
    if ymd < 0 or tod < 0 or tod > 24 * 3600:
        return 0
    yr = ymd // 10000
    mon = (ymd - yr * 10000) // 100
    day = ymd - yr * 10000 - mon * 100
    out_p[0] = yr
    out_p[1] = mon
    out_p[2] = day
    return 1


@export
def timegetymd_codon(yr: int, mon: int, day: int) -> int:
    if yr < 0:
        return -1
    return yr * 10000 + mon * 100 + day


@export
def make_tridiag_deriv_diag_codon(
    spr: Ptr[float],
    sub: Ptr[float],
    diag: Ptr[float],
    left_bound: Ptr[float],
    right_bound: Ptr[float],
    nsys: int,
    ncel: int,
):
    for k in range(ncel - 1):
        for i in range(nsys):
            diag[i + k * nsys] = -spr[i + k * nsys]

    last = ncel - 1
    for i in range(nsys):
        diag[i + last * nsys] = -right_bound[i]
        diag[i] = diag[i] - left_bound[i]

    for k in range(1, ncel):
        for i in range(nsys):
            diag[i + k * nsys] = diag[i + k * nsys] - sub[i + (k - 1) * nsys]


@export
def add_in_place_tridiag_ops_codon(
    spr: Ptr[float],
    sub: Ptr[float],
    diag: Ptr[float],
    left_bound: Ptr[float],
    right_bound: Ptr[float],
    other_spr: Ptr[float],
    other_sub: Ptr[float],
    other_diag: Ptr[float],
    other_left_bound: Ptr[float],
    other_right_bound: Ptr[float],
    nsys: int,
    ncel: int,
):
    off_count = nsys * (ncel - 1)
    diag_count = nsys * ncel
    for i in range(off_count):
        spr[i] += other_spr[i]
        sub[i] += other_sub[i]
    for i in range(diag_count):
        diag[i] += other_diag[i]
    for i in range(nsys):
        left_bound[i] += other_left_bound[i]
        right_bound[i] += other_right_bound[i]


@export
def scalar_lmult_tridiag_codon(
    spr: Ptr[float],
    sub: Ptr[float],
    diag: Ptr[float],
    left_bound: Ptr[float],
    right_bound: Ptr[float],
    nsys: int,
    ncel: int,
    constant: float,
):
    off_count = nsys * (ncel - 1)
    diag_count = nsys * ncel
    for i in range(off_count):
        spr[i] = spr[i] * constant
        sub[i] = sub[i] * constant
    for i in range(diag_count):
        diag[i] = diag[i] * constant
    for i in range(nsys):
        left_bound[i] = left_bound[i] * constant
        right_bound[i] = right_bound[i] * constant


@export
def finalize_dims_codon(dims_p: Ptr[int]):
    dims_p[0] = 0
    dims_p[1] = 0


@export
def tridiag_finalize_codon(dims_p: Ptr[int]):
    finalize_dims_codon(dims_p)


@export
def decomp_finalize_codon(dims_p: Ptr[int]):
    finalize_dims_codon(dims_p)


@export
def diagonal_operator_codon(
    spr: Ptr[float],
    sub: Ptr[float],
    diag: Ptr[float],
    left_bound: Ptr[float],
    right_bound: Ptr[float],
    in_diag: Ptr[float],
    nsys: int,
    ncel: int,
):
    off_count = nsys * (ncel - 1)
    diag_count = nsys * ncel
    for i in range(off_count):
        spr[i] = 0.0
        sub[i] = 0.0
    for i in range(diag_count):
        diag[i] = in_diag[i]
    for i in range(nsys):
        left_bound[i] = 0.0
        right_bound[i] = 0.0


@export
def diffusion_operator_codon(
    spr: Ptr[float],
    sub: Ptr[float],
    diag: Ptr[float],
    left_bound: Ptr[float],
    right_bound: Ptr[float],
    d_coef: Ptr[float],
    rdst: Ptr[float],
    rdel: Ptr[float],
    delta: Ptr[float],
    l_edge_width: Ptr[float],
    r_edge_width: Ptr[float],
    nsys: int,
    ncel: int,
    l_bndry_type: int,
    r_bndry_type: int,
):
    fixed_layer_bndry = 3

    if l_bndry_type == fixed_layer_bndry:
        for i in range(nsys):
            left_bound[i] = (
                2.0 * d_coef[i] * rdel[i] / (l_edge_width[i] + delta[i])
            )
    else:
        for i in range(nsys):
            left_bound[i] = 0.0

    for k in range(ncel - 1):
        off = k * nsys
        dcoef_off = (k + 1) * nsys
        next_cell_off = (k + 1) * nsys
        for i in range(nsys):
            flux_term = d_coef[i + dcoef_off] * rdst[i + off]
            spr[i + off] = flux_term * rdel[i + off]
            sub[i + off] = flux_term * rdel[i + next_cell_off]

    if r_bndry_type == fixed_layer_bndry:
        dcoef_right_off = ncel * nsys
        cell_right_off = (ncel - 1) * nsys
        for i in range(nsys):
            right_bound[i] = (
                2.0
                * d_coef[i + dcoef_right_off]
                * rdel[i + cell_right_off]
                / (r_edge_width[i] + delta[i + cell_right_off])
            )
    else:
        for i in range(nsys):
            right_bound[i] = 0.0

    make_tridiag_deriv_diag_codon(
        spr, sub, diag, left_bound, right_bound, nsys, ncel
    )


@export
def new_tridiagdecomp_codon(
    ca: Ptr[float],
    dnom: Ptr[float],
    ze: Ptr[float],
    op_spr: Ptr[float],
    op_sub: Ptr[float],
    op_diag: Ptr[float],
    op_left_bound: Ptr[float],
    op_right_bound: Ptr[float],
    nsys: int,
    ncel: int,
):
    for k in range(ncel - 1):
        off = k * nsys
        for i in range(nsys):
            ca[i + off] = -op_spr[i + off]

    last = ncel - 1
    last_off = last * nsys
    for i in range(nsys):
        ca[i + last_off] = -op_right_bound[i]
        dnom[i + last_off] = 1.0 / op_diag[i + last_off]

    for k in range(ncel - 2, -1, -1):
        off = k * nsys
        next_off = (k + 1) * nsys
        for i in range(nsys):
            ze[i + next_off] = -op_sub[i + off] * dnom[i + next_off]
            dnom[i + off] = 1.0 / (
                op_diag[i + off] - ca[i + off] * ze[i + next_off]
            )

    for i in range(nsys):
        ze[i] = -op_left_bound[i]


@export
def new_tridiagdecomp_graft_codon(
    ca: Ptr[float],
    dnom: Ptr[float],
    ze: Ptr[float],
    op_spr: Ptr[float],
    op_sub: Ptr[float],
    op_diag: Ptr[float],
    op_left_bound: Ptr[float],
    op_right_bound: Ptr[float],
    graft_ca: Ptr[float],
    graft_dnom: Ptr[float],
    graft_ze: Ptr[float],
    nsys: int,
    op_ncel: int,
    decomp_ncel: int,
):
    for k in range(op_ncel - 1):
        off = k * nsys
        for i in range(nsys):
            ca[i + off] = -op_spr[i + off]

    edge = op_ncel - 1
    edge_off = edge * nsys
    for i in range(nsys):
        ca[i + edge_off] = -op_right_bound[i]

    for k in range(op_ncel, decomp_ncel):
        off = k * nsys
        for i in range(nsys):
            ca[i + off] = graft_ca[i + off]
            dnom[i + off] = graft_dnom[i + off]
            ze[i + off] = graft_ze[i + off]

    next_off = op_ncel * nsys
    for i in range(nsys):
        dnom[i + edge_off] = 1.0 / (
            op_diag[i + edge_off] - ca[i + edge_off] * ze[i + next_off]
        )

    for k in range(op_ncel - 2, -1, -1):
        off = k * nsys
        next_k_off = (k + 1) * nsys
        for i in range(nsys):
            ze[i + next_k_off] = -op_sub[i + off] * dnom[i + next_k_off]
            dnom[i + off] = 1.0 / (
                op_diag[i + off] - ca[i + off] * ze[i + next_k_off]
            )

    for i in range(nsys):
        ze[i] = -op_left_bound[i]


@export
def decomp_left_div_codon(
    ca: Ptr[float],
    dnom: Ptr[float],
    ze: Ptr[float],
    q: Ptr[float],
    zf: Ptr[float],
    l_edge_data: Ptr[float],
    r_edge_data: Ptr[float],
    nsys: int,
    ncel: int,
    has_l_cond: int,
    l_cond_type: int,
    has_r_cond: int,
    r_cond_type: int,
):
    no_data_cond = 0
    data_cond = 1
    flux_cond = 2

    if has_l_cond != 0:
        if l_cond_type == no_data_cond:
            third_off = 2 * nsys
            for i in range(nsys):
                q[i] = q[i] + ze[i] * q[i + third_off]
        elif l_cond_type == data_cond:
            for i in range(nsys):
                q[i] = q[i] + ze[i] * l_edge_data[i]
        elif l_cond_type == flux_cond:
            for i in range(nsys):
                q[i] = q[i] + l_edge_data[i]

    last = ncel - 1
    last_off = last * nsys
    if has_r_cond != 0:
        if r_cond_type == no_data_cond:
            edge_data_off = (ncel - 3) * nsys
            for i in range(nsys):
                q[i + last_off] = q[i + last_off] + ca[i + last_off] * q[i + edge_data_off]
        elif r_cond_type == data_cond:
            for i in range(nsys):
                q[i + last_off] = q[i + last_off] + ca[i + last_off] * r_edge_data[i]
        elif r_cond_type == flux_cond:
            for i in range(nsys):
                q[i + last_off] = q[i + last_off] + r_edge_data[i]

    for i in range(nsys):
        zf[i + last_off] = q[i + last_off] * dnom[i + last_off]

    for k in range(ncel - 2, -1, -1):
        off = k * nsys
        next_off = (k + 1) * nsys
        for i in range(nsys):
            zf[i + off] = (q[i + off] + ca[i + off] * zf[i + next_off]) * dnom[i + off]

    for i in range(nsys):
        q[i] = zf[i]

    for k in range(1, ncel):
        off = k * nsys
        prev_off = (k - 1) * nsys
        for i in range(nsys):
            q[i + off] = zf[i + off] + ze[i + off] * q[i + prev_off]


@export
def ceil2_codon(n: int) -> int:
    p = 1
    while p < n:
        p = p * 2
    return p


@export
def physconst_init_codon(
    count: int,
    cpair: float,
    rair: float,
    mwdry: float,
    cpairv: Ptr[float],
    rairv: Ptr[float],
    cappav: Ptr[float],
    mbarv: Ptr[float],
):
    ratio = rair / cpair
    for i in range(count):
        cpairv[i] = cpair
        rairv[i] = rair
        cappav[i] = ratio
        mbarv[i] = mwdry


@export
def lininterp_finish_codon(tag: int) -> int:
    return tag

@export
def radsw_init_codon(tag: int) -> int:
    return tag

@export
def radiation_printopts_codon(tag: int) -> int:
    return tag

@export
def radiation_register_codon(tag: int) -> int:
    return tag

@export
def radlw_init_codon(tag: int) -> int:
    return tag

@export
def distnl_codon(tag: int) -> int:
    return tag

@export
def bld_outfld_hash_tbls_codon(tag: int) -> int:
    return tag

@export
def add_hist_coord_r8_codon(tag: int) -> int:
    return tag

@export
def getfil_codon(tag: int) -> int:
    return tag

@export
def cam_write_restart_codon(tag: int) -> int:
    return tag

@export
def atm_domain_mct_codon(tag: int) -> int:
    return tag

@export
def timemgr_write_restart_codon(tag: int) -> int:
    return tag

@export
def timemgr_print_codon(tag: int) -> int:
    return tag

@export
def sat_hist_readnl_codon(tag: int) -> int:
    return tag

@export
def physconst_readnl_codon(tag: int) -> int:
    return tag

@export
def get_field_properties_codon(tag: int) -> int:
    return tag

@export
def get_phys_ldof_codon(tag: int) -> int:
    return tag

@export
def wrapup_codon(tag: int) -> int:
    return tag

@export
def dump_field_codon(tag: int) -> int:
    return tag

@export
def intht_codon(tag: int) -> int:
    return tag

@export
def get_phys_decomp_mdnd_codon(tag: int) -> int:
    return tag

@export
def infld_real_2d_codon(tag: int) -> int:
    return tag

@export
def hycoef_init_codon(tag: int) -> int:
    return tag

@export
def write_restart_history_codon(tag: int) -> int:
    return tag

@export
def restart_vars_setnames_codon(tag: int) -> int:
    return tag

@export
def atm_export_codon(tag: int) -> int:
    return tag

@export
def column_init_codon(tag: int) -> int:
    return tag

@export
def infld_real_3d_codon(tag: int) -> int:
    return tag

@export
def lininterp_init_codon(tag: int) -> int:
    return tag

@export
def wshist_codon(tag: int) -> int:
    return tag

@export
def addfld_codon(tag: int) -> int:
    return tag

@export
def fldlst_codon(tag: int) -> int:
    return tag

@export
def read_namelist_codon(tag: int) -> int:
    return tag

@export
def h_define_codon(tag: int) -> int:
    return tag

@export
def check_hist_coord_all_codon(tag: int) -> int:
    return tag

@export
def new_coords1d_from_int_codon(tag: int) -> int:
    return tag

@export
def cam_run1_codon(tag: int) -> int:
    return tag
