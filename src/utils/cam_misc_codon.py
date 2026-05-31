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
