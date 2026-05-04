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


@inline
def _idx2(i: int, k: int, ld1: int) -> int:
    """Fortran array(i,k) with first dimension ld1."""
    return (i - 1) + (k - 1) * ld1


@inline
def _idx3(i: int, k: int, m: int, ld1: int, ld2: int) -> int:
    """Fortran array(i,k,m) with dimensions (ld1,ld2,*)."""
    return (i - 1) + (k - 1) * ld1 + (m - 1) * ld1 * ld2


@export
def gw_tend_oro_post_codon(
    stage: int,
    ncol: int,
    psetcols: int,
    pcols: int,
    pver: int,
    pcnst: int,
    pverp: int,
    egwdffi_tot_p: cobj,
    egwdffi_p: cobj,
    utgw_p: cobj,
    vtgw_p: cobj,
    ttgw_p: cobj,
    cpairv_p: cobj,
    ptend_u_p: cobj,
    ptend_v_p: cobj,
    ptend_s_p: cobj,
    qtgw_p: cobj,
    ptend_q_p: cobj,
    tau_p: cobj,
    xv_p: cobj,
    yv_p: cobj,
    tau0x_p: cobj,
    tau0y_p: cobj,
):
    egwdffi_tot = Ptr[float](egwdffi_tot_p)
    egwdffi = Ptr[float](egwdffi_p)
    utgw = Ptr[float](utgw_p)
    vtgw = Ptr[float](vtgw_p)
    ttgw = Ptr[float](ttgw_p)
    cpairv = Ptr[float](cpairv_p)
    ptend_u = Ptr[float](ptend_u_p)
    ptend_v = Ptr[float](ptend_v_p)
    ptend_s = Ptr[float](ptend_s_p)
    qtgw = Ptr[float](qtgw_p)
    ptend_q = Ptr[float](ptend_q_p)
    tau = Ptr[float](tau_p)
    xv = Ptr[float](xv_p)
    yv = Ptr[float](yv_p)
    tau0x = Ptr[float](tau0x_p)
    tau0y = Ptr[float](tau0y_p)

    if stage == 1:
        for k in range(1, pverp + 1):
            for i in range(1, ncol + 1):
                idx_local = _idx2(i, k, ncol)
                egwdffi_tot[idx_local] = egwdffi_tot[idx_local] + egwdffi[idx_local]

        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx_local = _idx2(i, k, ncol)
                idx_ptend = _idx2(i, k, psetcols)
                ptend_u[idx_ptend] = ptend_u[idx_ptend] + utgw[idx_local]
                ptend_v[idx_ptend] = ptend_v[idx_ptend] + vtgw[idx_local]
                ptend_s[idx_ptend] = ptend_s[idx_ptend] + ttgw[idx_local]
                ttgw[idx_local] = ttgw[idx_local] / cpairv[_idx2(i, k, pcols)]

    elif stage == 2:
        for m in range(1, pcnst + 1):
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    ptend_idx = _idx3(i, k, m, psetcols, pver)
                    local_idx = _idx3(i, k, m, ncol, pver)
                    ptend_q[ptend_idx] = ptend_q[ptend_idx] + qtgw[local_idx]

        tau_surface_offset = pver * ncol
        for i in range(1, ncol + 1):
            tau_surface = tau[(i - 1) + tau_surface_offset]
            tau0x[i - 1] = tau_surface * xv[i - 1]
            tau0y[i - 1] = tau_surface * yv[i - 1]
