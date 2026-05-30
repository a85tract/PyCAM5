@export
def cam_misc_touch_codon(tag: int) -> int:
    return tag


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


@export
def handle_pio_error_ok_codon(ierr: int, pio_noerr: int) -> int:
    return 1 if ierr == pio_noerr else 0


@export
def cam_initfile_getter_touch_codon(tag: int) -> int:
    return tag


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
def sat_hist_init_noop_codon(has_sat_hist: int) -> int:
    return 1 if has_sat_hist == 0 else 0


@export
def scalar_add_tridiag_codon(diag: Ptr[float], nsys: int, ncel: int, constant: float):
    for i in range(nsys * ncel):
        diag[i] += constant


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
