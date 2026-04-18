@inline
def _idx2(i: int, k: int, ld1: int) -> int:
    return (i - 1) + (k - 1) * ld1


@inline
def _idx3(i: int, k: int, m: int, ld1: int, ld2: int) -> int:
    return (i - 1) + (k - 1) * ld1 + (m - 1) * ld1 * ld2


@export
def aero_model_drydep_select_branches_codon(
    apply_srf_drydep: int,
    branch_mask_p: cobj,
):
    branch_mask = Ptr[int](branch_mask_p)

    mask = 0
    if apply_srf_drydep != 0:
        mask |= 1

    branch_mask[0] = mask


@export
def aero_model_wetdep_f_act_conv_coarse_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    dt: float,
    lcoardust: int,
    lcoarnacl: int,
    state_q_p: cobj,
    ptend_q_p: cobj,
    f_act_conv_coarse_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    ptend_q = Ptr[float](ptend_q_p)
    f_act_conv_coarse = Ptr[float](f_act_conv_coarse_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            f_act_conv_coarse[_idx2(i, k, pcols)] = 0.60

    if lcoardust <= 0 or lcoarnacl <= 0:
        return

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            tmpdust = state_q[_idx3(i, k, lcoardust, pcols, pver)] + dt * ptend_q[
                _idx3(i, k, lcoardust, pcols, pver)
            ]
            if tmpdust < 0.0:
                tmpdust = 0.0
            tmpnacl = state_q[_idx3(i, k, lcoarnacl, pcols, pver)] + dt * ptend_q[
                _idx3(i, k, lcoarnacl, pcols, pver)
            ]
            if tmpnacl < 0.0:
                tmpnacl = 0.0
            if tmpdust + tmpnacl > 1.0e-30:
                f_act_conv_coarse[_idx2(i, k, pcols)] = (
                    0.40 * tmpdust + 0.80 * tmpnacl
                ) / (tmpdust + tmpnacl)


@export
def aero_model_wetdep_column_flux_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gravit: float,
    field_p: cobj,
    pdel_p: cobj,
    sflx_p: cobj,
):
    field = Ptr[float](field_p)
    pdel = Ptr[float](pdel_p)
    sflx = Ptr[float](sflx_p)

    for i in range(1, ncol + 1):
        total = 0.0
        for k in range(1, pver + 1):
            total += field[_idx2(i, k, pcols)] * pdel[_idx2(i, k, pcols)] / gravit
        sflx[i - 1] = total
