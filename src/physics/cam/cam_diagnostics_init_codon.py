def cam_diag_init_dqcond_num_codon(history_budget: int, conv_tend_code: int, pcnst: int) -> int:
    if history_budget != 0:
        return pcnst
    if conv_tend_code == 1:
        return 1
    if conv_tend_code == 2:
        return pcnst
    return 0
