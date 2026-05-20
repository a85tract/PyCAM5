import se_dynamics_common_codon as _common

def euler_step_vstar_prepare_codon(
    np: int,
    nlev: int,
    dt: float,
    rhs_multiplier: int,
    dp_in_p: cobj,
    divdp_proj_p: cobj,
    vn0_p: cobj,
    dp_out_p: cobj,
    vstar_p: cobj,
):
    dp_in = Ptr[float](dp_in_p)
    divdp_proj = Ptr[float](divdp_proj_p)
    vn0 = Ptr[float](vn0_p)
    dp_out = Ptr[float](dp_out_p)
    vstar = Ptr[float](vstar_p)

    for k in range(1, nlev + 1):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                vol_idx = _common._vol_idx(i, j, k, np)
                dp_val = dp_in[vol_idx] - rhs_multiplier * dt * divdp_proj[vol_idx]
                dp_out[vol_idx] = dp_val
                vstar[_common._v_idx(i, j, 1, k, np)] = vn0[_common._v_idx(i, j, 1, k, np)] / dp_val
                vstar[_common._v_idx(i, j, 2, k, np)] = vn0[_common._v_idx(i, j, 2, k, np)] / dp_val


def euler_step_limiter_dpstar_codon(
    np: int,
    nlev: int,
    dt: float,
    rhs_viss: int,
    nu_q: float,
    nu_p: float,
    dp_in_p: cobj,
    divdp_p: cobj,
    dpdiss_biharmonic_p: cobj,
    spheremp_p: cobj,
    dp_star_p: cobj,
):
    dp_in = Ptr[float](dp_in_p)
    divdp = Ptr[float](divdp_p)
    dpdiss_biharmonic = Ptr[float](dpdiss_biharmonic_p)
    spheremp = Ptr[float](spheremp_p)
    dp_star = Ptr[float](dp_star_p)

    use_dpdiss = nu_p > 0.0 and rhs_viss != 0

    for k in range(1, nlev + 1):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                vol_idx = _common._vol_idx(i, j, k, np)
                plane_idx = _common._plane_idx(i, j, np)
                dp_val = dp_in[vol_idx] - dt * divdp[vol_idx]
                if use_dpdiss:
                    dp_val = dp_val - rhs_viss * dt * nu_q * dpdiss_biharmonic[vol_idx] / spheremp[plane_idx]
                dp_star[vol_idx] = dp_val


def euler_step_qdp_writeback_codon(
    np: int,
    nlev: int,
    qsize: int,
    qidx: int,
    np1_qdp: int,
    qdp_p: cobj,
    spheremp_p: cobj,
    qtens_p: cobj,
):
    qdp = Ptr[float](qdp_p)
    spheremp = Ptr[float](spheremp_p)
    qtens = Ptr[float](qtens_p)

    for k in range(1, nlev + 1):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                plane_idx = _common._plane_idx(i, j, np)
                vol_idx = _common._vol_idx(i, j, k, np)
                qdp[_common._q_tl_idx(i, j, k, qidx, np1_qdp, np, nlev, qsize)] = spheremp[plane_idx] * qtens[vol_idx]


def euler_step_qdp_restore_codon(
    np: int,
    nlev: int,
    qdp_p: cobj,
    rspheremp_p: cobj,
):
    qdp = Ptr[float](qdp_p)
    rspheremp = Ptr[float](rspheremp_p)

    for k in range(1, nlev + 1):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                plane_idx = _common._plane_idx(i, j, np)
                vol_idx = _common._vol_idx(i, j, k, np)
                qdp[vol_idx] = rspheremp[plane_idx] * qdp[vol_idx]


def advance_hypervis_qtens_prepare_codon(
    np: int,
    nlev: int,
    qsize: int,
    ps0: float,
    dt2: float,
    nu_p: float,
    hyai_p: cobj,
    hybi_p: cobj,
    dp_in_p: cobj,
    divdp_proj_p: cobj,
    dpdiss_ave_p: cobj,
    qdp_p: cobj,
    dp_out_p: cobj,
    qtens_p: cobj,
):
    hyai = Ptr[float](hyai_p)
    hybi = Ptr[float](hybi_p)
    dp_in = Ptr[float](dp_in_p)
    divdp_proj = Ptr[float](divdp_proj_p)
    dpdiss_ave = Ptr[float](dpdiss_ave_p)
    qdp = Ptr[float](qdp_p)
    dp_out = Ptr[float](dp_out_p)
    qtens = Ptr[float](qtens_p)

    use_dpdiss = nu_p > 0.0

    for k in range(1, nlev + 1):
        dp0 = (hyai[_common._hy_idx(k + 1)] - hyai[_common._hy_idx(k)]) * ps0 + (hybi[_common._hy_idx(k + 1)] - hybi[_common._hy_idx(k)]) * ps0
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                vol_idx = _common._vol_idx(i, j, k, np)
                dp_val = dp_in[vol_idx] - dt2 * divdp_proj[vol_idx]
                dp_out[vol_idx] = dp_val
                if use_dpdiss:
                    for q in range(1, qsize + 1):
                        q_idx = _common._q_idx(i, j, k, q, np, nlev)
                        qtens[q_idx] = dpdiss_ave[vol_idx] * qdp[q_idx] / dp_val
                else:
                    for q in range(1, qsize + 1):
                        q_idx = _common._q_idx(i, j, k, q, np, nlev)
                        qtens[q_idx] = dp0 * qdp[q_idx] / dp_val


def advance_hypervis_qdp_update_codon(
    np: int,
    nlev: int,
    dt: float,
    nu_q: float,
    qdp_p: cobj,
    spheremp_p: cobj,
    qtens_p: cobj,
):
    qdp = Ptr[float](qdp_p)
    spheremp = Ptr[float](spheremp_p)
    qtens = Ptr[float](qtens_p)

    for k in range(1, nlev + 1):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                plane_idx = _common._plane_idx(i, j, np)
                vol_idx = _common._vol_idx(i, j, k, np)
                qdp[vol_idx] = qdp[vol_idx] * spheremp[plane_idx] - dt * nu_q * qtens[vol_idx]


def euler_step_dssvar_restore_codon(
    np: int,
    nlev: int,
    dssvar_p: cobj,
    rspheremp_p: cobj,
):
    dssvar = Ptr[float](dssvar_p)
    rspheremp = Ptr[float](rspheremp_p)

    for k in range(1, nlev + 1):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                plane_idx = _common._plane_idx(i, j, np)
                vol_idx = _common._vol_idx(i, j, k, np)
                dssvar[vol_idx] = dssvar[vol_idx] * rspheremp[plane_idx]


def euler_step_dssvar_pack_codon(
    np: int,
    nlev: int,
    dssvar_p: cobj,
    spheremp_p: cobj,
):
    dssvar = Ptr[float](dssvar_p)
    spheremp = Ptr[float](spheremp_p)

    for k in range(1, nlev + 1):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                plane_idx = _common._plane_idx(i, j, np)
                vol_idx = _common._vol_idx(i, j, k, np)
                dssvar[vol_idx] = spheremp[plane_idx] * dssvar[vol_idx]


def euler_step_qtens_base_codon(
    np: int,
    dt: float,
    qdp_p: cobj,
    dp_star_p: cobj,
    qtens_p: cobj,
):
    qdp = Ptr[float](qdp_p)
    dp_star = Ptr[float](dp_star_p)
    qtens = Ptr[float](qtens_p)

    for j in range(1, np + 1):
        for i in range(1, np + 1):
            plane_idx = _common._plane_idx(i, j, np)
            qtens[plane_idx] = qdp[plane_idx] - dt * dp_star[plane_idx]


def euler_step_gradq_prepare_codon(
    np: int,
    vstar1_p: cobj,
    vstar2_p: cobj,
    qdp_p: cobj,
    gradq1_p: cobj,
    gradq2_p: cobj,
):
    vstar1 = Ptr[float](vstar1_p)
    vstar2 = Ptr[float](vstar2_p)
    qdp = Ptr[float](qdp_p)
    gradq1 = Ptr[float](gradq1_p)
    gradq2 = Ptr[float](gradq2_p)

    for j in range(1, np + 1):
        for i in range(1, np + 1):
            plane_idx = _common._plane_idx(i, j, np)
            gradq1[plane_idx] = vstar1[plane_idx] * qdp[plane_idx]
            gradq2[plane_idx] = vstar2[plane_idx] * qdp[plane_idx]


def euler_step_qtens_biharmonic_add_codon(
    np: int,
    qtens_p: cobj,
    qtens_biharmonic_p: cobj,
):
    qtens = Ptr[float](qtens_p)
    qtens_biharmonic = Ptr[float](qtens_biharmonic_p)

    for j in range(1, np + 1):
        for i in range(1, np + 1):
            plane_idx = _common._plane_idx(i, j, np)
            qtens[plane_idx] = qtens[plane_idx] + qtens_biharmonic[plane_idx]


def euler_step_qtens_biharmonic_init_codon(
    np: int,
    nlev: int,
    qsize: int,
    dt: float,
    rhs_multiplier: int,
    dp_in_p: cobj,
    divdp_proj_p: cobj,
    qdp_p: cobj,
    dp_out_p: cobj,
    qtens_biharmonic_p: cobj,
):
    dp_in = Ptr[float](dp_in_p)
    divdp_proj = Ptr[float](divdp_proj_p)
    qdp = Ptr[float](qdp_p)
    dp_out = Ptr[float](dp_out_p)
    qtens_biharmonic = Ptr[float](qtens_biharmonic_p)

    for k in range(1, nlev + 1):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                vol_idx = _common._vol_idx(i, j, k, np)
                dp_val = dp_in[vol_idx] - rhs_multiplier * dt * divdp_proj[vol_idx]
                dp_out[vol_idx] = dp_val
                for q in range(1, qsize + 1):
                    q_idx = _common._q_idx(i, j, k, q, np, nlev)
                    qtens_biharmonic[q_idx] = qdp[q_idx] / dp_val


def euler_step_qtens_biharmonic_scale_codon(
    np: int,
    nlev: int,
    qsize: int,
    qtens_biharmonic_p: cobj,
    dpdiss_ave_p: cobj,
    dp0_p: cobj,
):
    qtens_biharmonic = Ptr[float](qtens_biharmonic_p)
    dpdiss_ave = Ptr[float](dpdiss_ave_p)
    dp0 = Ptr[float](dp0_p)

    for q in range(1, qsize + 1):
        for k in range(1, nlev + 1):
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    q_idx = _common._q_idx(i, j, k, q, np, nlev)
                    vol_idx = _common._vol_idx(i, j, k, np)
                    qtens_biharmonic[q_idx] = qtens_biharmonic[q_idx] * dpdiss_ave[vol_idx] / dp0[_common._col_idx(k)]


def euler_step_qtens_biharmonic_unapply_codon(
    np: int,
    nlev: int,
    rhs_viss: int,
    dt: float,
    nu_q: float,
    qtens_biharmonic_p: cobj,
    spheremp_p: cobj,
    dp0_p: cobj,
):
    qtens_biharmonic = Ptr[float](qtens_biharmonic_p)
    spheremp = Ptr[float](spheremp_p)
    dp0 = Ptr[float](dp0_p)

    for k in range(1, nlev + 1):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                vol_idx = _common._vol_idx(i, j, k, np)
                plane_idx = _common._plane_idx(i, j, np)
                qtens_biharmonic[vol_idx] = -rhs_viss * dt * nu_q * dp0[_common._col_idx(k)] * qtens_biharmonic[vol_idx] / spheremp[plane_idx]


def euler_step_qminmax_update_codon(
    np: int,
    nlev: int,
    qsize: int,
    rhs_multiplier: int,
    qtens_biharmonic_p: cobj,
    qmin_p: cobj,
    qmax_p: cobj,
):
    qtens_biharmonic = Ptr[float](qtens_biharmonic_p)
    qmin = Ptr[float](qmin_p)
    qmax = Ptr[float](qmax_p)

    for q in range(1, qsize + 1):
        for k in range(1, nlev + 1):
            qmin_val = 1.0e24
            qmax_val = -1.0e24
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    val = qtens_biharmonic[_common._q_idx(i, j, k, q, np, nlev)]
                    qmin_val = min(qmin_val, val)
                    qmax_val = max(qmax_val, val)

            lev_q_idx = _common._lev_q_idx(k, q, nlev)
            if rhs_multiplier == 1:
                qmin[lev_q_idx] = min(qmin[lev_q_idx], qmin_val)
                qmin[lev_q_idx] = max(qmin[lev_q_idx], 0.0)
                qmax[lev_q_idx] = max(qmax[lev_q_idx], qmax_val)
            else:
                qmin[lev_q_idx] = max(qmin_val, 0.0)
                qmax[lev_q_idx] = qmax_val


def limiter2d_zero_codon(
    np: int,
    nlev: int,
    q_p: cobj,
):
    q = Ptr[float](q_p)

    for k in range(nlev, 0, -1):
        mass = 0.0
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                plane_idx = _common._vol_idx(i, j, k, np)
                mass = mass + q[plane_idx]

        if mass < 0.0:
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    plane_idx = _common._vol_idx(i, j, k, np)
                    q[plane_idx] = -q[plane_idx]

        mass_new = 0.0
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                plane_idx = _common._vol_idx(i, j, k, np)
                if q[plane_idx] < 0.0:
                    q[plane_idx] = 0.0
                else:
                    mass_new = mass_new + q[plane_idx]

        if mass_new > 0.0:
            scale = abs(mass) / mass_new
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    plane_idx = _common._vol_idx(i, j, k, np)
                    q[plane_idx] = q[plane_idx] * scale

        if mass < 0.0:
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    plane_idx = _common._vol_idx(i, j, k, np)
                    q[plane_idx] = -q[plane_idx]


def limiter_optim_iter_full_codon(
    np: int,
    nlev: int,
    ptens_p: cobj,
    sphweights_p: cobj,
    minp_p: cobj,
    maxp_p: cobj,
    dpmass_p: cobj,
    weights_p: cobj,
    whois_neg_p: cobj,
    whois_pos_p: cobj,
    x_p: cobj,
    c_p: cobj,
    al_neg_p: cobj,
    al_pos_p: cobj,
):
    ptens = Ptr[float](ptens_p)
    sphweights = Ptr[float](sphweights_p)
    minp = Ptr[float](minp_p)
    maxp = Ptr[float](maxp_p)
    dpmass = Ptr[float](dpmass_p)
    weights = Ptr[float](weights_p)
    whois_neg = Ptr[int](whois_neg_p)
    whois_pos = Ptr[int](whois_pos_p)
    x = Ptr[float](x_p)
    c = Ptr[float](c_p)
    al_neg = Ptr[float](al_neg_p)
    al_pos = Ptr[float](al_pos_p)

    ncols = np * np
    tol_limiter = 1.0e-15
    maxiter = 5

    for k in range(1, nlev + 1):
        for k1 in range(1, ncols + 1):
            wk_idx = _common._cell_lev_idx(k1, k, np)
            weights[wk_idx] = sphweights[k1 - 1] * dpmass[wk_idx]
            ptens[wk_idx] = ptens[wk_idx] / dpmass[wk_idx]

    for k in range(1, nlev + 1):
        mass = 0.0
        sumc = 0.0
        for k1 in range(1, ncols + 1):
            wk_idx = _common._cell_lev_idx(k1, k, np)
            c[k1 - 1] = weights[wk_idx]
            x[k1 - 1] = ptens[wk_idx]
            mass = mass + c[k1 - 1] * x[k1 - 1]
            sumc = sumc + c[k1 - 1]

        if (mass / sumc) < minp[_common._col_idx(k)]:
            minp[_common._col_idx(k)] = mass / sumc
        if (mass / sumc) > maxp[_common._col_idx(k)]:
            maxp[_common._col_idx(k)] = mass / sumc

        addmass = 0.0
        pos_counter = 0
        neg_counter = 0

        for k1 in range(1, ncols + 1):
            if x[k1 - 1] >= maxp[_common._col_idx(k)]:
                addmass = addmass + (x[k1 - 1] - maxp[_common._col_idx(k)]) * c[k1 - 1]
                x[k1 - 1] = maxp[_common._col_idx(k)]
                whois_pos[k1 - 1] = -1
            else:
                pos_counter = pos_counter + 1
                whois_pos[pos_counter - 1] = k1

            if x[k1 - 1] <= minp[_common._col_idx(k)]:
                addmass = addmass - (minp[_common._col_idx(k)] - x[k1 - 1]) * c[k1 - 1]
                x[k1 - 1] = minp[_common._col_idx(k)]
                whois_neg[k1 - 1] = -1
            else:
                neg_counter = neg_counter + 1
                whois_neg[neg_counter - 1] = k1

        if addmass > 0.0:
            for _iter in range(1, maxiter + 1):
                weightssum = 0.0
                for k1 in range(1, pos_counter + 1):
                    i1 = whois_pos[k1 - 1]
                    weightssum = weightssum + c[i1 - 1]
                    al_pos[i1 - 1] = maxp[_common._col_idx(k)] - x[i1 - 1]

                if pos_counter > 0 and addmass > tol_limiter * abs(mass):
                    for k1 in range(1, pos_counter + 1):
                        i1 = whois_pos[k1 - 1]
                        howmuch = addmass / weightssum
                        if howmuch > al_pos[i1 - 1]:
                            howmuch = al_pos[i1 - 1]
                            whois_pos[k1 - 1] = -1
                        addmass = addmass - howmuch * c[i1 - 1]
                        weightssum = weightssum - c[i1 - 1]
                        x[i1 - 1] = x[i1 - 1] + howmuch

                    neg_counter = pos_counter
                    for k1 in range(1, ncols + 1):
                        whois_neg[k1 - 1] = whois_pos[k1 - 1]
                        whois_pos[k1 - 1] = -1
                    pos_counter = 0
                    for k1 in range(1, neg_counter + 1):
                        if whois_neg[k1 - 1] != -1:
                            pos_counter = pos_counter + 1
                            whois_pos[pos_counter - 1] = whois_neg[k1 - 1]
                else:
                    break
        else:
            for _iter in range(1, maxiter + 1):
                weightssum = 0.0
                for k1 in range(1, neg_counter + 1):
                    i1 = whois_neg[k1 - 1]
                    weightssum = weightssum + c[i1 - 1]
                    al_neg[i1 - 1] = x[i1 - 1] - minp[_common._col_idx(k)]

                if neg_counter > 0 and (-addmass) > tol_limiter * abs(mass):
                    for k1 in range(1, neg_counter + 1):
                        i1 = whois_neg[k1 - 1]
                        howmuch = -addmass / weightssum
                        if howmuch > al_neg[i1 - 1]:
                            howmuch = al_neg[i1 - 1]
                            whois_neg[k1 - 1] = -1
                        addmass = addmass + howmuch * c[i1 - 1]
                        weightssum = weightssum - c[i1 - 1]
                        x[i1 - 1] = x[i1 - 1] - howmuch

                    pos_counter = neg_counter
                    for k1 in range(1, ncols + 1):
                        whois_pos[k1 - 1] = whois_neg[k1 - 1]
                        whois_neg[k1 - 1] = -1
                    neg_counter = 0
                    for k1 in range(1, pos_counter + 1):
                        if whois_pos[k1 - 1] != -1:
                            neg_counter = neg_counter + 1
                            whois_neg[neg_counter - 1] = whois_pos[k1 - 1]
                else:
                    break

        for k1 in range(1, ncols + 1):
            ptens[_common._cell_lev_idx(k1, k, np)] = x[k1 - 1]

    for k in range(1, nlev + 1):
        for k1 in range(1, ncols + 1):
            wk_idx = _common._cell_lev_idx(k1, k, np)
            ptens[wk_idx] = ptens[wk_idx] * dpmass[wk_idx]
