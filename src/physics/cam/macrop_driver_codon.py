@export
def macrop_driver_select_branches_codon(
    micro_do_icesupersat: int,
    trace_water: int,
    wtrc_detrain_in_macrop: int,
    cu_det_st: int,
    use_shfrc: int,
    do_cldice: int,
    do_cldliq: int,
    do_detrain: int,
    branch_mask_p: cobj,
):
    branch_mask = Ptr[int](branch_mask_p)

    mask = 0
    if micro_do_icesupersat != 0:
        mask |= 1
    if trace_water != 0:
        mask |= 2
    if wtrc_detrain_in_macrop != 0:
        mask |= 4
    if cu_det_st != 0:
        mask |= 8
    if use_shfrc != 0:
        mask |= 16
    if do_cldice != 0:
        mask |= 32
    if do_cldliq != 0:
        mask |= 64
    if do_detrain != 0:
        mask |= 128

    branch_mask[0] = mask
