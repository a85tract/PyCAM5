from math import sqrt


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


@inline
def _idx_tau(i: int, l: int, k: int, ncol: int, ngwv: int) -> int:
    """Fortran tau(i,l,k) with dimensions (ncol,-ngwv:ngwv,*)."""
    return (i - 1) + (l + ngwv) * ncol + (k - 1) * ncol * (2 * ngwv + 1)


@inline
def _idx_c(i: int, l: int, ncol: int, ngwv: int) -> int:
    """Fortran c(i,l) with dimensions (ncol,-ngwv:ngwv)."""
    return (i - 1) + (l + ngwv) * ncol


@export
def gw_oro_src_codon(
    ncol: int,
    pver: int,
    ngwv: int,
    fcrit2: float,
    kwv: float,
    rair: float,
    p_mid_p: cobj,
    p_del_p: cobj,
    p_ifc_p: cobj,
    u_p: cobj,
    v_p: cobj,
    t_p: cobj,
    sgh_p: cobj,
    zm_p: cobj,
    nm_p: cobj,
    src_level_p: cobj,
    tend_level_p: cobj,
    tau_p: cobj,
    ubm_p: cobj,
    ubi_p: cobj,
    xv_p: cobj,
    yv_p: cobj,
    c_p: cobj,
    hdsp_p: cobj,
    tauoro_p: cobj,
    nsrc_p: cobj,
    rsrc_p: cobj,
    usrc_p: cobj,
    vsrc_p: cobj,
    dpsrc_p: cobj,
):
    p_mid = Ptr[float](p_mid_p)
    p_del = Ptr[float](p_del_p)
    p_ifc = Ptr[float](p_ifc_p)
    u = Ptr[float](u_p)
    v = Ptr[float](v_p)
    t = Ptr[float](t_p)
    sgh = Ptr[float](sgh_p)
    zm = Ptr[float](zm_p)
    nm = Ptr[float](nm_p)
    src_level = Ptr[int](src_level_p)
    tend_level = Ptr[int](tend_level_p)
    tau = Ptr[float](tau_p)
    ubm = Ptr[float](ubm_p)
    ubi = Ptr[float](ubi_p)
    xv = Ptr[float](xv_p)
    yv = Ptr[float](yv_p)
    c = Ptr[float](c_p)
    hdsp = Ptr[float](hdsp_p)
    tauoro = Ptr[float](tauoro_p)
    nsrc = Ptr[float](nsrc_p)
    rsrc = Ptr[float](rsrc_p)
    usrc = Ptr[float](usrc_p)
    vsrc = Ptr[float](vsrc_p)
    dpsrc = Ptr[float](dpsrc_p)

    orohmin = 10.0
    orovmin = 2.0

    for i in range(1, ncol + 1):
        hdsp[i - 1] = 2.0 * sgh[i - 1]

    k = pver
    for i in range(1, ncol + 1):
        src_level[i - 1] = k - 1
        rsrc[i - 1] = (p_mid[_idx2(i, k, ncol)] / (rair * t[_idx2(i, k, ncol)])) * p_del[_idx2(i, k, ncol)]
        usrc[i - 1] = u[_idx2(i, k, ncol)] * p_del[_idx2(i, k, ncol)]
        vsrc[i - 1] = v[_idx2(i, k, ncol)] * p_del[_idx2(i, k, ncol)]
        nsrc[i - 1] = nm[_idx2(i, k, ncol)] * p_del[_idx2(i, k, ncol)]

    for k in range(pver - 1, 0, -1):
        for i in range(1, ncol + 1):
            if hdsp[i - 1] > sqrt(zm[_idx2(i, k, ncol)] * zm[_idx2(i, k + 1, ncol)]):
                src_level[i - 1] = k - 1
                rsrc[i - 1] = rsrc[i - 1] + (
                    p_mid[_idx2(i, k, ncol)] / (rair * t[_idx2(i, k, ncol)])
                ) * p_del[_idx2(i, k, ncol)]
                usrc[i - 1] = usrc[i - 1] + u[_idx2(i, k, ncol)] * p_del[_idx2(i, k, ncol)]
                vsrc[i - 1] = vsrc[i - 1] + v[_idx2(i, k, ncol)] * p_del[_idx2(i, k, ncol)]
                nsrc[i - 1] = nsrc[i - 1] + nm[_idx2(i, k, ncol)] * p_del[_idx2(i, k, ncol)]

        all_found = True
        for i in range(1, ncol + 1):
            if src_level[i - 1] < k:
                all_found = False
        if all_found:
            break

    for i in range(1, ncol + 1):
        dpsrc[i - 1] = p_ifc[_idx2(i, pver + 1, ncol)] - p_ifc[_idx2(i, src_level[i - 1] + 1, ncol)]

    for i in range(1, ncol + 1):
        rsrc[i - 1] = rsrc[i - 1] / dpsrc[i - 1]
        usrc[i - 1] = usrc[i - 1] / dpsrc[i - 1]
        vsrc[i - 1] = vsrc[i - 1] / dpsrc[i - 1]
        nsrc[i - 1] = nsrc[i - 1] / dpsrc[i - 1]

    for i in range(1, ncol + 1):
        ubi[_idx2(i, pver + 1, ncol)] = sqrt(usrc[i - 1] * usrc[i - 1] + vsrc[i - 1] * vsrc[i - 1])

    for i in range(1, ncol + 1):
        mag = ubi[_idx2(i, pver + 1, ncol)]
        if mag > 0.0:
            xv[i - 1] = usrc[i - 1] / mag
            yv[i - 1] = vsrc[i - 1] / mag
        else:
            xv[i - 1] = 0.0
            yv[i - 1] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            ubm[_idx2(i, k, ncol)] = u[_idx2(i, k, ncol)] * xv[i - 1] + v[_idx2(i, k, ncol)] * yv[i - 1]

    for i in range(1, ncol + 1):
        ubi[_idx2(i, 1, ncol)] = ubm[_idx2(i, 1, ncol)]

    for k in range(2, pver + 1):
        for i in range(1, ncol + 1):
            ubi[_idx2(i, k, ncol)] = 0.5 * (ubm[_idx2(i, k - 1, ncol)] + ubm[_idx2(i, k, ncol)])

    for i in range(1, ncol + 1):
        ubi_surface = ubi[_idx2(i, pver + 1, ncol)]
        if (ubi_surface > orovmin) and (hdsp[i - 1] > orohmin):
            sghmax = fcrit2 * ((ubi_surface / nsrc[i - 1]) ** 2)
            tauoro[i - 1] = 0.5 * kwv * min(hdsp[i - 1] ** 2, sghmax) * rsrc[i - 1] * nsrc[i - 1] * ubi_surface
        else:
            tauoro[i - 1] = 0.0
            src_level[i - 1] = pver

    for k in range(1, pver + 2):
        for l in range(-ngwv, ngwv + 1):
            for i in range(1, ncol + 1):
                tau[_idx_tau(i, l, k, ncol, ngwv)] = 0.0

    min_src_level = src_level[0]
    for i in range(2, ncol + 1):
        if src_level[i - 1] < min_src_level:
            min_src_level = src_level[i - 1]

    for k in range(pver, min_src_level - 1, -1):
        for i in range(1, ncol + 1):
            if src_level[i - 1] <= k:
                tau[_idx_tau(i, 0, k + 1, ncol, ngwv)] = tauoro[i - 1]

    for i in range(1, ncol + 1):
        tend_level[i - 1] = pver

    for l in range(-ngwv, ngwv + 1):
        for i in range(1, ncol + 1):
            c[_idx_c(i, l, ncol, ngwv)] = 0.0


@export
def gw_tend_prep_codon(
    stage: int,
    ncol: int,
    psetcols: int,
    pcols: int,
    pver: int,
    pverp: int,
    pcnst: int,
    effgw_oro: float,
    eps: float,
    state_s_p: cobj,
    state_t_p: cobj,
    state_u_p: cobj,
    state_v_p: cobj,
    state_q_p: cobj,
    state_lnpint_p: cobj,
    state_zm_p: cobj,
    dse_p: cobj,
    t_p: cobj,
    u_p: cobj,
    v_p: cobj,
    q_p: cobj,
    piln_p: cobj,
    zm_p: cobj,
    egwdffi_tot_p: cobj,
    flx_heat_p: cobj,
    landfrac_p: cobj,
    sgh_p: cobj,
    effgw_p: cobj,
    sgh_scaled_p: cobj,
):
    state_s = Ptr[float](state_s_p)
    state_t = Ptr[float](state_t_p)
    state_u = Ptr[float](state_u_p)
    state_v = Ptr[float](state_v_p)
    state_q = Ptr[float](state_q_p)
    state_lnpint = Ptr[float](state_lnpint_p)
    state_zm = Ptr[float](state_zm_p)
    dse = Ptr[float](dse_p)
    t = Ptr[float](t_p)
    u = Ptr[float](u_p)
    v = Ptr[float](v_p)
    q = Ptr[float](q_p)
    piln = Ptr[float](piln_p)
    zm = Ptr[float](zm_p)
    egwdffi_tot = Ptr[float](egwdffi_tot_p)
    flx_heat = Ptr[float](flx_heat_p)
    landfrac = Ptr[float](landfrac_p)
    sgh = Ptr[float](sgh_p)
    effgw = Ptr[float](effgw_p)
    sgh_scaled = Ptr[float](sgh_scaled_p)

    if stage == 1:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                src_idx = _idx2(i, k, psetcols)
                dst_idx = _idx2(i, k, ncol)
                dse[dst_idx] = state_s[src_idx]
                t[dst_idx] = state_t[src_idx]
                u[dst_idx] = state_u[src_idx]
                v[dst_idx] = state_v[src_idx]
                zm[dst_idx] = state_zm[src_idx]

        for k in range(1, pverp + 1):
            for i in range(1, ncol + 1):
                piln[_idx2(i, k, ncol)] = state_lnpint[_idx2(i, k, psetcols)]
                egwdffi_tot[_idx2(i, k, ncol)] = 0.0

        for m in range(1, pcnst + 1):
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    q[_idx3(i, k, m, ncol, pver)] = state_q[_idx3(i, k, m, psetcols, pver)]

        for i in range(1, pcols + 1):
            flx_heat[i - 1] = 0.0

    elif stage == 2:
        for i in range(1, ncol + 1):
            landfrac_i = landfrac[i - 1]
            if landfrac_i >= eps:
                effgw[i - 1] = effgw_oro * landfrac_i
                sgh_scaled[i - 1] = sgh[i - 1] / sqrt(landfrac_i)
            else:
                effgw[i - 1] = 0.0
                sgh_scaled[i - 1] = 0.0


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
