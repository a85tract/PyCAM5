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
