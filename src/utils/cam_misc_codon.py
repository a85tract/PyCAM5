@export
def cam_misc_touch_codon(tag: int) -> int:
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
def handle_pio_error_ok_codon(ierr: int, pio_noerr: int) -> int:
    return 1 if ierr == pio_noerr else 0


@export
def cam_initfile_getter_touch_codon(tag: int) -> int:
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
