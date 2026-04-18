@export
def vertical_diffusion_ts_init_codon():
    return


@export
def vertical_diffusion_tend_select_branches_codon(
    do_tms: int,
    do_molec_diff: int,
    use_diag_tke: int,
    use_hb_family: int,
    shallow_unicon: int,
    prog_modal_aero: int,
    do_pseudocon_diff: int,
    diff_cnsrv_mass_check: int,
    waccmx_special: int,
    branch_mask_p: cobj,
):
    branch_mask = Ptr[int](branch_mask_p)

    mask = 0
    if do_tms != 0:
        mask |= 1
    if do_molec_diff != 0:
        mask |= 2
    if use_diag_tke != 0:
        mask |= 4
    if use_hb_family != 0:
        mask |= 8
    if shallow_unicon != 0:
        mask |= 16
    if prog_modal_aero != 0:
        mask |= 32
    if do_pseudocon_diff != 0:
        mask |= 64
    if diff_cnsrv_mass_check != 0:
        mask |= 128
    if waccmx_special != 0:
        mask |= 256

    branch_mask[0] = mask
