@export
def gw_tend_select_branches_codon(
    do_molec_diff: int,
    use_gw_convect_dp: int,
    use_gw_convect_sh: int,
    use_gw_front: int,
    use_gw_front_igw: int,
    use_gw_oro: int,
    branch_mask_p: cobj,
):
    branch_mask = Ptr[int](branch_mask_p)

    mask = 0
    if do_molec_diff != 0:
        mask |= 1
    if use_gw_convect_dp != 0:
        mask |= 2
    if use_gw_convect_sh != 0:
        mask |= 4
    if use_gw_front != 0:
        mask |= 8
    if use_gw_front_igw != 0:
        mask |= 16
    if use_gw_oro != 0:
        mask |= 32

    branch_mask[0] = mask
