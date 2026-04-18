@export
def stratiform_select_branches_codon(
    use_shfrc: int,
    cam3: int,
    branch_mask_p: cobj,
):
    branch_mask = Ptr[int](branch_mask_p)

    mask = 0
    if use_shfrc != 0:
        mask |= 1
    if cam3 != 0:
        mask |= 2

    branch_mask[0] = mask
