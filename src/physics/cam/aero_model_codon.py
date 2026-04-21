from math import atan, copysign, erfc, exp, log, sqrt

@inline
def _idx2(i: int, k: int, ld1: int) -> int:
    """Fortran array declared as (ld1, *)."""
    return (i - 1) + (k - 1) * ld1


@inline
def _idx3(i: int, k: int, m: int, ld1: int, ld2: int) -> int:
    """Fortran array declared as (ld1, ld2, *)."""
    return (i - 1) + (k - 1) * ld1 + (m - 1) * ld1 * ld2


@inline
def _max3(a: float, b: float, c: float) -> float:
    if a >= b:
        if a >= c:
            return a
        return c
    if b >= c:
        return b
    return c


@inline
def _min4(a: float, b: float, c: float, d: float) -> float:
    return min(min(a, b), min(c, d))


@inline
def _max4(a: float, b: float, c: float, d: float) -> float:
    return max(max(a, b), max(c, d))


@inline
def _sign_one(x: float) -> float:
    return copysign(1.0, x)


@inline
def _minmod(a: float, b: float) -> float:
    return 0.5 * (_sign_one(a) + _sign_one(b)) * min(abs(a), abs(b))


@inline
def _medan(a: float, b: float, c: float) -> float:
    return a + _minmod(b - a, c - a)


@export
def modal_aero_depvel_part_codon(
    ncol: int,
    pcols: int,
    pver: int,
    n_land_type: int,
    moment: int,
    pi: float,
    boltz: float,
    gravit: float,
    rair: float,
    t_p: cobj,
    pmid_p: cobj,
    ram1_p: cobj,
    fv_p: cobj,
    vlc_dry_p: cobj,
    vlc_trb_p: cobj,
    vlc_grv_p: cobj,
    radius_part_p: cobj,
    density_part_p: cobj,
    sig_part_p: cobj,
    fraction_landuse_p: cobj,
    vsc_dyn_atm_p: cobj,
    vsc_knm_atm_p: cobj,
    mfp_atm_p: cobj,
    slp_crc_p: cobj,
    radius_moment_p: cobj,
):
    t = Ptr[float](t_p)
    pmid = Ptr[float](pmid_p)
    ram1 = Ptr[float](ram1_p)
    fv = Ptr[float](fv_p)
    vlc_dry = Ptr[float](vlc_dry_p)
    vlc_trb = Ptr[float](vlc_trb_p)
    vlc_grv = Ptr[float](vlc_grv_p)
    radius_part = Ptr[float](radius_part_p)
    density_part = Ptr[float](density_part_p)
    sig_part = Ptr[float](sig_part_p)
    fraction_landuse = Ptr[float](fraction_landuse_p)
    vsc_dyn_atm = Ptr[float](vsc_dyn_atm_p)
    vsc_knm_atm = Ptr[float](vsc_knm_atm_p)
    mfp_atm = Ptr[float](mfp_atm_p)
    slp_crc = Ptr[float](slp_crc_p)
    radius_moment = Ptr[float](radius_moment_p)

    gamma = (0.56, 0.54, 0.54, 0.56, 0.56, 0.56, 0.50, 0.54, 0.54, 0.54, 0.54)
    alpha = (1.50, 1.20, 1.20, 0.80, 1.00, 0.80, 100.00, 50.00, 2.00, 1.20, 50.00)
    radius_collector = (
        10.00e-03,
        3.50e-03,
        3.50e-03,
        5.10e-03,
        2.00e-03,
        5.00e-03,
        -1.00,
        -1.00,
        10.00e-03,
        3.50e-03,
        -1.00,
    )
    iwet = (-1, -1, -1, -1, -1, -1, 1, -1, 1, -1, -1)

    moment_f = float(moment)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            lnsig = log(sig_part[_idx2(i, k, pcols)])
            radius_moment[_idx2(i, k, pcols)] = min(
                50.0e-6, radius_part[_idx2(i, k, pcols)]
            ) * exp((moment_f - 1.5) * lnsig * lnsig)
            dispersion = exp(2.0 * lnsig * lnsig)

            rho = pmid[_idx2(i, k, pcols)] / rair / t[_idx2(i, k, pcols)]

            vsc_dyn_atm[_idx2(i, k, pcols)] = (
                1.72e-5
                * ((t[_idx2(i, k, pcols)] / 273.0) ** 1.5)
                * 393.0
                / (t[_idx2(i, k, pcols)] + 120.0)
            )
            mfp_atm[_idx2(i, k, pcols)] = (
                2.0
                * vsc_dyn_atm[_idx2(i, k, pcols)]
                / (
                    pmid[_idx2(i, k, pcols)]
                    * sqrt(8.0 / (pi * rair * t[_idx2(i, k, pcols)]))
                )
            )
            vsc_knm_atm[_idx2(i, k, pcols)] = vsc_dyn_atm[_idx2(i, k, pcols)] / rho

            slp_crc[_idx2(i, k, pcols)] = (
                1.0
                + mfp_atm[_idx2(i, k, pcols)]
                * (
                    1.257
                    + 0.4
                    * exp(
                        -1.1
                        * radius_moment[_idx2(i, k, pcols)]
                        / mfp_atm[_idx2(i, k, pcols)]
                    )
                )
                / radius_moment[_idx2(i, k, pcols)]
            )
            vlc_grv[_idx2(i, k, pcols)] = (
                (4.0 / 18.0)
                * radius_moment[_idx2(i, k, pcols)]
                * radius_moment[_idx2(i, k, pcols)]
                * density_part[_idx2(i, k, pcols)]
                * gravit
                * slp_crc[_idx2(i, k, pcols)]
                / vsc_dyn_atm[_idx2(i, k, pcols)]
            )
            vlc_grv[_idx2(i, k, pcols)] = vlc_grv[_idx2(i, k, pcols)] * dispersion

            vlc_dry[_idx2(i, k, pcols)] = vlc_grv[_idx2(i, k, pcols)]

    k = pver
    for i in range(1, ncol + 1):
        dff_aer = (
            boltz
            * t[_idx2(i, k, pcols)]
            * slp_crc[_idx2(i, k, pcols)]
            / (6.0 * pi * vsc_dyn_atm[_idx2(i, k, pcols)] * radius_moment[_idx2(i, k, pcols)])
        )
        shm_nbr = vsc_knm_atm[_idx2(i, k, pcols)] / dff_aer

        wrk2 = 0.0
        wrk3 = 0.0
        for lt in range(1, n_land_type + 1):
            lnd_frc = fraction_landuse[_idx2(i, lt, pcols)]
            if lnd_frc != 0.0:
                brownian = shm_nbr ** (-gamma[lt - 1])
                if radius_collector[lt - 1] > 0.0:
                    stk_nbr = (
                        vlc_grv[_idx2(i, k, pcols)] * fv[i - 1] / (gravit * radius_collector[lt - 1])
                    )
                    interception = 2.0 * (
                        radius_moment[_idx2(i, k, pcols)] / radius_collector[lt - 1]
                    ) ** 2.0
                else:
                    stk_nbr = (
                        vlc_grv[_idx2(i, k, pcols)]
                        * fv[i - 1]
                        * fv[i - 1]
                        / (gravit * vsc_knm_atm[_idx2(i, k, pcols)])
                    )
                    interception = 0.0
                impaction = (stk_nbr / (alpha[lt - 1] + stk_nbr)) ** 2.0

                if iwet[lt - 1] > 0:
                    stickfrac = 1.0
                else:
                    stickfrac = exp(-sqrt(stk_nbr))
                    if stickfrac < 1.0e-10:
                        stickfrac = 1.0e-10
                rss_lmn = 1.0 / (
                    3.0 * fv[i - 1] * stickfrac * (brownian + interception + impaction)
                )
                rss_trb = (
                    ram1[i - 1]
                    + rss_lmn
                    + ram1[i - 1] * rss_lmn * vlc_grv[_idx2(i, k, pcols)]
                )

                wrk1 = 1.0 / rss_trb
                wrk2 = wrk2 + lnd_frc * wrk1
                wrk3 = wrk3 + lnd_frc * (wrk1 + vlc_grv[_idx2(i, k, pcols)])
        vlc_trb[i - 1] = wrk2
        vlc_dry[_idx2(i, k, pcols)] = wrk3


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
def modal_aero_bcscavcoef_get_codon(
    m: int,
    ncol: int,
    pcols: int,
    pver: int,
    ntot_amode: int,
    nimptblgrow_mind: int,
    nimptblgrow_maxd: int,
    dlndg_nimptblgrow: float,
    dgnum_mode: float,
    isprx_mask_p: cobj,
    dgn_awet_p: cobj,
    scavimptblnum_mode_p: cobj,
    scavimptblvol_mode_p: cobj,
    scavcoefnum_p: cobj,
    scavcoefvol_p: cobj,
):
    isprx_mask = Ptr[int](isprx_mask_p)
    dgn_awet = Ptr[float](dgn_awet_p)
    scavimptblnum_mode = Ptr[float](scavimptblnum_mode_p)
    scavimptblvol_mode = Ptr[float](scavimptblvol_mode_p)
    scavcoefnum = Ptr[float](scavcoefnum_p)
    scavcoefvol = Ptr[float](scavcoefvol_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            if isprx_mask[_idx2(i, k, pcols)] != 0:
                dumdgratio = dgn_awet[_idx3(i, k, m, pcols, pver)] / dgnum_mode

                if dumdgratio >= 0.99 and dumdgratio <= 1.01:
                    tbl_idx = 0 - nimptblgrow_mind
                    scavimpvol = scavimptblvol_mode[tbl_idx]
                    scavimpnum = scavimptblnum_mode[tbl_idx]
                else:
                    xgrow = log(dumdgratio) / dlndg_nimptblgrow
                    jgrow = int(xgrow)
                    if xgrow < 0.0:
                        jgrow = jgrow - 1
                    if jgrow < nimptblgrow_mind:
                        jgrow = nimptblgrow_mind
                        xgrow = float(jgrow)
                    else:
                        jgrow = min(jgrow, nimptblgrow_maxd - 1)

                    dumfhi = xgrow - jgrow
                    dumflo = 1.0 - dumfhi
                    tbl_idx = jgrow - nimptblgrow_mind

                    scavimpvol = (
                        dumflo * scavimptblvol_mode[tbl_idx]
                        + dumfhi * scavimptblvol_mode[tbl_idx + 1]
                    )
                    scavimpnum = (
                        dumflo * scavimptblnum_mode[tbl_idx]
                        + dumfhi * scavimptblnum_mode[tbl_idx + 1]
                    )

                scavcoefvol[_idx2(i, k, pcols)] = exp(scavimpvol)
                scavcoefnum[_idx2(i, k, pcols)] = exp(scavimpnum)
            else:
                scavcoefvol[_idx2(i, k, pcols)] = 0.0
                scavcoefnum[_idx2(i, k, pcols)] = 0.0


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


@export
def calcram_codon(
    ncol: int,
    pcols: int,
    rair: float,
    gravit: float,
    ram1in_p: cobj,
    fvin_p: cobj,
    ram1_p: cobj,
    fv_p: cobj,
    obklen_p: cobj,
    ustar_p: cobj,
    landfrac_p: cobj,
    icefrac_p: cobj,
    ocnfrac_p: cobj,
    t_p: cobj,
    pmid_p: cobj,
    pdel_p: cobj,
):
    ram1in = Ptr[float](ram1in_p)
    fvin = Ptr[float](fvin_p)
    ram1 = Ptr[float](ram1_p)
    fv = Ptr[float](fv_p)
    obklen = Ptr[float](obklen_p)
    ustar = Ptr[float](ustar_p)
    landfrac = Ptr[float](landfrac_p)
    icefrac = Ptr[float](icefrac_p)
    ocnfrac = Ptr[float](ocnfrac_p)
    t = Ptr[float](t_p)
    pmid = Ptr[float](pmid_p)
    pdel = Ptr[float](pdel_p)

    zzocen = 0.0001
    zzsice = 0.0400
    xkar = 0.4

    for i in range(1, ncol + 1):
        z = pdel[i - 1] * rair * t[i - 1] / pmid[i - 1] / gravit / 2.0
        if obklen[i - 1] == 0.0:
            psi = 0.0
            psi0 = 0.0
        else:
            psi = min(max(z / obklen[i - 1], -1.0), 1.0)
            psi0 = min(max(zzocen / obklen[i - 1], -1.0), 1.0)

        temp = z / zzocen
        if icefrac[i - 1] > 0.5:
            if obklen[i - 1] > 0.0:
                psi0 = min(max(zzsice / obklen[i - 1], -1.0), 1.0)
            else:
                psi0 = 0.0
            temp = z / zzsice

        if psi > 0.0:
            ram = 1.0 / xkar / ustar[i - 1] * (log(temp) + 4.7 * (psi - psi0))
        else:
            nu = (1.0 - 15.0 * psi) ** 0.25
            nu0 = (1.0 - 15.0 * psi0) ** 0.25
            if ustar[i - 1] != 0.0:
                ram = 1.0 / xkar / ustar[i - 1] * (
                    log(temp)
                    + log(
                        ((nu0**2.0 + 1.0) * (nu0 + 1.0) ** 2.0)
                        / ((nu**2.0 + 1.0) * (nu + 1.0) ** 2.0)
                    )
                    + 2.0 * (atan(nu) - atan(nu0))
                )
            else:
                ram = 0.0

        if landfrac[i - 1] < 0.000000001:
            fv[i - 1] = ustar[i - 1]
            ram1[i - 1] = ram
        else:
            fv[i - 1] = fvin[i - 1]
            ram1[i - 1] = ram1in[i - 1]

    for i in range(1, ncol + 1):
        if fv[i - 1] == 0.0:
            fv[i - 1] = 1.0e-12


def _dust_cfint2(
    ncol: int,
    pcols: int,
    pverp: int,
    xin_k: int,
    x_p: cobj,
    f_p: cobj,
    fdot_p: cobj,
    xxk_p: cobj,
    fxdot_p: cobj,
    fxdd_p: cobj,
    psistar_p: cobj,
    xins_p: cobj,
    intz_p: cobj,
    status_p: cobj,
    fail_i_p: cobj,
    fail_k_p: cobj,
):
    x = Ptr[float](x_p)
    f = Ptr[float](f_p)
    fdot = Ptr[float](fdot_p)
    xxk = Ptr[float](xxk_p)
    fxdot = Ptr[float](fxdot_p)
    fxdd = Ptr[float](fxdd_p)
    psistar = Ptr[float](psistar_p)
    xins = Ptr[float](xins_p)
    intz = Ptr[int](intz_p)
    status = Ptr[int](status_p)
    fail_i = Ptr[int](fail_i_p)
    fail_k = Ptr[int](fail_k_p)

    for i in range(1, ncol + 1):
        xins[i - 1] = _medan(
            x[_idx2(i, 1, pcols)],
            xxk[_idx2(i, xin_k, pcols)],
            x[_idx2(i, pverp, pcols)],
        )
        intz[i - 1] = 0

    for k in range(1, pverp):
        for i in range(1, ncol + 1):
            if (
                (xins[i - 1] - x[_idx2(i, k, pcols)])
                * (x[_idx2(i, k + 1, pcols)] - xins[i - 1])
            ) >= 0.0:
                intz[i - 1] = k

    for i in range(1, ncol + 1):
        if intz[i - 1] == 0:
            status[0] = 1
            fail_i[0] = i
            fail_k[0] = xin_k
            return

    for i in range(1, ncol + 1):
        k = int(intz[i - 1])
        dx = x[_idx2(i, k + 1, pcols)] - x[_idx2(i, k, pcols)]
        s = (f[_idx2(i, k + 1, pcols)] - f[_idx2(i, k, pcols)]) / dx
        c2 = (3.0 * s - 2.0 * fdot[_idx2(i, k, pcols)] - fdot[_idx2(i, k + 1, pcols)]) / dx
        c3 = (
            fdot[_idx2(i, k, pcols)] + fdot[_idx2(i, k + 1, pcols)] - 2.0 * s
        ) / (dx * dx)
        xx = xins[i - 1] - x[_idx2(i, k, pcols)]
        fxdot[i - 1] = (3.0 * c3 * xx + 2.0 * c2) * xx + fdot[_idx2(i, k, pcols)]
        fxdd[i - 1] = 6.0 * c3 * xx + 2.0 * c2
        cfint = ((c3 * xx + c2) * xx + fdot[_idx2(i, k, pcols)]) * xx + f[_idx2(i, k, pcols)]

        psi1 = f[_idx2(i, k, pcols)] + (
            (f[_idx2(i, k + 1, pcols)] - f[_idx2(i, k, pcols)]) * xx / dx
        )
        if k == 1:
            psi2 = f[_idx2(i, 1, pcols)]
        else:
            psi2 = f[_idx2(i, k, pcols)] + (
                (f[_idx2(i, k, pcols)] - f[_idx2(i, k - 1, pcols)])
                * xx
                / (x[_idx2(i, k, pcols)] - x[_idx2(i, k - 1, pcols)])
            )

        if (k + 1) == pverp:
            psi3 = f[_idx2(i, pverp, pcols)]
        else:
            psi3 = f[_idx2(i, k + 1, pcols)] - (
                (f[_idx2(i, k + 2, pcols)] - f[_idx2(i, k + 1, pcols)])
                * (dx - xx)
                / (x[_idx2(i, k + 2, pcols)] - x[_idx2(i, k + 1, pcols)])
            )

        psim = _medan(psi1, psi2, psi3)
        cfnew = _medan(cfint, psi1, psim)
        psistar[i - 1] = cfnew


def _dust_cfdotmc_pro(
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    x_p: cobj,
    f_p: cobj,
    fdot_p: cobj,
    s_p: cobj,
    sh_p: cobj,
    d_p: cobj,
    dh_p: cobj,
    e_p: cobj,
    eh_p: cobj,
    ppl_p: cobj,
    ppr_p: cobj,
    delxh_p: cobj,
):
    x = Ptr[float](x_p)
    f = Ptr[float](f_p)
    fdot = Ptr[float](fdot_p)
    s = Ptr[float](s_p)
    sh = Ptr[float](sh_p)
    d = Ptr[float](d_p)
    dh = Ptr[float](dh_p)
    e = Ptr[float](e_p)
    eh = Ptr[float](eh_p)
    ppl = Ptr[float](ppl_p)
    ppr = Ptr[float](ppr_p)
    delxh = Ptr[float](delxh_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            delxh[_idx2(i, k, pcols)] = x[_idx2(i, k + 1, pcols)] - x[_idx2(i, k, pcols)]
            sh[_idx2(i, k, pcols)] = (
                f[_idx2(i, k + 1, pcols)] - f[_idx2(i, k, pcols)]
            ) / delxh[_idx2(i, k, pcols)]

        if k >= 2:
            for i in range(1, ncol + 1):
                d[_idx2(i, k, pcols)] = (
                    sh[_idx2(i, k, pcols)] - sh[_idx2(i, k - 1, pcols)]
                ) / (x[_idx2(i, k + 1, pcols)] - x[_idx2(i, k - 1, pcols)])
                s[_idx2(i, k, pcols)] = _minmod(
                    sh[_idx2(i, k, pcols)], sh[_idx2(i, k - 1, pcols)]
                )

    for k in range(2, pver):
        for i in range(1, ncol + 1):
            eh[_idx2(i, k, pcols)] = (
                d[_idx2(i, k + 1, pcols)] - d[_idx2(i, k, pcols)]
            ) / (x[_idx2(i, k + 2, pcols)] - x[_idx2(i, k - 1, pcols)])
            dh[_idx2(i, k, pcols)] = _minmod(
                d[_idx2(i, k, pcols)], d[_idx2(i, k + 1, pcols)]
            )

    for i in range(1, ncol + 1):
        e[_idx2(i, 2, pcols)] = eh[_idx2(i, 2, pcols)]
        e[_idx2(i, pver, pcols)] = eh[_idx2(i, pver - 1, pcols)]

        fdot[_idx2(i, 1, pcols)] = (
            sh[_idx2(i, 1, pcols)]
            - d[_idx2(i, 2, pcols)] * delxh[_idx2(i, 1, pcols)]
            - eh[_idx2(i, 2, pcols)]
            * delxh[_idx2(i, 1, pcols)]
            * (x[_idx2(i, 1, pcols)] - x[_idx2(i, 3, pcols)])
        )
        fdot[_idx2(i, 1, pcols)] = _minmod(
            fdot[_idx2(i, 1, pcols)], 3.0 * sh[_idx2(i, 1, pcols)]
        )

        fdot[_idx2(i, pverp, pcols)] = (
            sh[_idx2(i, pver, pcols)]
            + d[_idx2(i, pver, pcols)] * delxh[_idx2(i, pver, pcols)]
            + eh[_idx2(i, pver - 1, pcols)]
            * delxh[_idx2(i, pver, pcols)]
            * (x[_idx2(i, pverp, pcols)] - x[_idx2(i, pver - 1, pcols)])
        )
        fdot[_idx2(i, pverp, pcols)] = _minmod(
            fdot[_idx2(i, pverp, pcols)], 3.0 * sh[_idx2(i, pver, pcols)]
        )

        fdot[_idx2(i, 2, pcols)] = (
            sh[_idx2(i, 1, pcols)]
            + d[_idx2(i, 2, pcols)] * delxh[_idx2(i, 1, pcols)]
            - eh[_idx2(i, 2, pcols)]
            * delxh[_idx2(i, 1, pcols)]
            * delxh[_idx2(i, 2, pcols)]
        )
        fdot[_idx2(i, 2, pcols)] = _minmod(
            fdot[_idx2(i, 2, pcols)], 3.0 * s[_idx2(i, 2, pcols)]
        )

        fdot[_idx2(i, pver, pcols)] = (
            sh[_idx2(i, pver, pcols)]
            - d[_idx2(i, pver, pcols)] * delxh[_idx2(i, pver, pcols)]
            - eh[_idx2(i, pver - 1, pcols)]
            * delxh[_idx2(i, pver, pcols)]
            * delxh[_idx2(i, pver - 1, pcols)]
        )
        fdot[_idx2(i, pver, pcols)] = _minmod(
            fdot[_idx2(i, pver, pcols)], 3.0 * s[_idx2(i, pver, pcols)]
        )

    for k in range(3, pver):
        for i in range(1, ncol + 1):
            e[_idx2(i, k, pcols)] = _minmod(
                eh[_idx2(i, k, pcols)], eh[_idx2(i, k - 1, pcols)]
            )

    for k in range(3, pver):
        for i in range(1, ncol + 1):
            ppl[_idx2(i, k, pcols)] = (
                sh[_idx2(i, k - 1, pcols)] + dh[_idx2(i, k - 1, pcols)] * delxh[_idx2(i, k - 1, pcols)]
            )
            ppr[_idx2(i, k, pcols)] = (
                sh[_idx2(i, k, pcols)] - dh[_idx2(i, k, pcols)] * delxh[_idx2(i, k, pcols)]
            )

            t = _minmod(ppl[_idx2(i, k, pcols)], ppr[_idx2(i, k, pcols)])

            pp = sh[_idx2(i, k - 1, pcols)] + d[_idx2(i, k, pcols)] * delxh[_idx2(i, k - 1, pcols)]

            fdot[_idx2(i, k, pcols)] = pp - (
                delxh[_idx2(i, k - 1, pcols)]
                * delxh[_idx2(i, k, pcols)]
                * (
                    eh[_idx2(i, k - 1, pcols)] * (x[_idx2(i, k + 2, pcols)] - x[_idx2(i, k, pcols)])
                    + eh[_idx2(i, k, pcols)] * (x[_idx2(i, k, pcols)] - x[_idx2(i, k - 2, pcols)])
                )
                / (x[_idx2(i, k + 2, pcols)] - x[_idx2(i, k - 2, pcols)])
            )

            qpl = sh[_idx2(i, k - 1, pcols)] + delxh[_idx2(i, k - 1, pcols)] * _minmod(
                d[_idx2(i, k - 1, pcols)]
                + e[_idx2(i, k - 1, pcols)] * (x[_idx2(i, k, pcols)] - x[_idx2(i, k - 2, pcols)]),
                d[_idx2(i, k, pcols)] - e[_idx2(i, k, pcols)] * delxh[_idx2(i, k, pcols)],
            )
            qpr = sh[_idx2(i, k, pcols)] + delxh[_idx2(i, k, pcols)] * _minmod(
                d[_idx2(i, k, pcols)] + e[_idx2(i, k, pcols)] * delxh[_idx2(i, k - 1, pcols)],
                d[_idx2(i, k + 1, pcols)]
                + e[_idx2(i, k + 1, pcols)] * (x[_idx2(i, k, pcols)] - x[_idx2(i, k + 2, pcols)]),
            )

            fdot[_idx2(i, k, pcols)] = _medan(fdot[_idx2(i, k, pcols)], qpl, qpr)

            ttt = _minmod(qpl, qpr)
            tmin = _min4(
                0.0,
                3.0 * s[_idx2(i, k, pcols)],
                1.5 * t,
                ttt,
            )
            tmax = _max4(
                0.0,
                3.0 * s[_idx2(i, k, pcols)],
                1.5 * t,
                ttt,
            )

            fdot[_idx2(i, k, pcols)] = fdot[_idx2(i, k, pcols)] + _minmod(
                tmin - fdot[_idx2(i, k, pcols)],
                tmax - fdot[_idx2(i, k, pcols)],
            )


def _dust_getflx(
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    deltat: float,
    xw_p: cobj,
    phi_p: cobj,
    vel_p: cobj,
    flux_p: cobj,
    psi_p: cobj,
    fdot_p: cobj,
    xxk_p: cobj,
    fxdot_p: cobj,
    fxdd_p: cobj,
    psistar_p: cobj,
    xins_p: cobj,
    intz_p: cobj,
    status_p: cobj,
    fail_i_p: cobj,
    fail_k_p: cobj,
    s_p: cobj,
    sh_p: cobj,
    d_p: cobj,
    dh_p: cobj,
    e_p: cobj,
    eh_p: cobj,
    ppl_p: cobj,
    ppr_p: cobj,
    delxh_p: cobj,
):
    xw = Ptr[float](xw_p)
    phi = Ptr[float](phi_p)
    vel = Ptr[float](vel_p)
    flux = Ptr[float](flux_p)
    psi = Ptr[float](psi_p)
    xxk = Ptr[float](xxk_p)
    psistar = Ptr[float](psistar_p)
    status = Ptr[int](status_p)

    for i in range(1, ncol + 1):
        psi[_idx2(i, 1, pcols)] = 0.0
        flux[_idx2(i, 1, pcols)] = 0.0
        flux[_idx2(i, pverp, pcols)] = 0.0

    for k in range(2, pverp + 1):
        for i in range(1, ncol + 1):
            psi[_idx2(i, k, pcols)] = (
                phi[_idx2(i, k - 1, pcols)]
                * (xw[_idx2(i, k, pcols)] - xw[_idx2(i, k - 1, pcols)])
                + psi[_idx2(i, k - 1, pcols)]
            )

    _dust_cfdotmc_pro(
        ncol, pcols, pver, pverp, xw_p, psi_p, fdot_p, s_p, sh_p, d_p, dh_p, e_p, eh_p, ppl_p, ppr_p, delxh_p
    )

    for k in range(2, pver + 1):
        for i in range(1, ncol + 1):
            xxk[_idx2(i, k, pcols)] = xw[_idx2(i, k, pcols)] - vel[_idx2(i, k, pcols)] * deltat

    for k in range(2, pver + 1):
        _dust_cfint2(
            ncol,
            pcols,
            pverp,
            k,
            xw_p,
            psi_p,
            fdot_p,
            xxk_p,
            fxdot_p,
            fxdd_p,
            psistar_p,
            xins_p,
            intz_p,
            status_p,
            fail_i_p,
            fail_k_p,
        )
        if status[0] != 0:
            return
        for i in range(1, ncol + 1):
            flux[_idx2(i, k, pcols)] = psi[_idx2(i, k, pcols)] - psistar[i - 1]


@export
def dust_sediment_tend_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    dtime: float,
    mxsedfac: float,
    gravit: float,
    pint_p: cobj,
    pdel_p: cobj,
    dustmr_p: cobj,
    pvdust_p: cobj,
    dusttend_p: cobj,
    sfdust_p: cobj,
    fxdust_p: cobj,
    psi_p: cobj,
    fdot_p: cobj,
    xxk_p: cobj,
    fxdot_p: cobj,
    fxdd_p: cobj,
    psistar_p: cobj,
    s_p: cobj,
    sh_p: cobj,
    d_p: cobj,
    dh_p: cobj,
    e_p: cobj,
    eh_p: cobj,
    ppl_p: cobj,
    ppr_p: cobj,
    delxh_p: cobj,
    xins_p: cobj,
    intz_p: cobj,
    status_p: cobj,
    fail_i_p: cobj,
    fail_k_p: cobj,
):
    pdel = Ptr[float](pdel_p)
    dustmr = Ptr[float](dustmr_p)
    pvdust = Ptr[float](pvdust_p)
    fxdust = Ptr[float](fxdust_p)
    dusttend = Ptr[float](dusttend_p)
    sfdust = Ptr[float](sfdust_p)
    status = Ptr[int](status_p)
    fail_i = Ptr[int](fail_i_p)
    fail_k = Ptr[int](fail_k_p)

    status[0] = 0
    fail_i[0] = 0
    fail_k[0] = 0

    for k in range(1, pverp + 1):
        for i in range(1, ncol + 1):
            fxdust[_idx2(i, k, pcols)] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            dusttend[_idx2(i, k, pcols)] = 0.0

    for i in range(1, ncol + 1):
        sfdust[i - 1] = 0.0

    _dust_getflx(
        ncol,
        pcols,
        pver,
        pverp,
        dtime,
        pint_p,
        dustmr_p,
        pvdust_p,
        fxdust_p,
        psi_p,
        fdot_p,
        xxk_p,
        fxdot_p,
        fxdd_p,
        psistar_p,
        xins_p,
        intz_p,
        status_p,
        fail_i_p,
        fail_k_p,
        s_p,
        sh_p,
        d_p,
        dh_p,
        e_p,
        eh_p,
        ppl_p,
        ppr_p,
        delxh_p,
    )
    if status[0] != 0:
        return

    for i in range(1, ncol + 1):
        fxdust[_idx2(i, 1, pcols)] = 0.0
        fxdust[_idx2(i, pverp, pcols)] = (
            dustmr[_idx2(i, pver, pcols)] * pvdust[_idx2(i, pverp, pcols)] * dtime
        )

    for k in range(2, pver + 1):
        for i in range(1, ncol + 1):
            fxdust[_idx2(i, k, pcols)] = max(0.0, fxdust[_idx2(i, k, pcols)])

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            fxdust[_idx2(i, k + 1, pcols)] = min(
                fxdust[_idx2(i, k + 1, pcols)],
                mxsedfac * dustmr[_idx2(i, k, pcols)] * pdel[_idx2(i, k, pcols)],
            )

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            dusttend[_idx2(i, k, pcols)] = (
                fxdust[_idx2(i, k, pcols)] - fxdust[_idx2(i, k + 1, pcols)]
            ) / (dtime * pdel[_idx2(i, k, pcols)])

    for i in range(1, ncol + 1):
        sfdust[i - 1] = fxdust[_idx2(i, pverp, pcols)] / (dtime * gravit)


@export
def qqcw2vmr_codon(
    ncol: int,
    pver: int,
    fldcw_ld1: int,
    mbar_p: cobj,
    fldcw_p: cobj,
    adv_mass: float,
    vmr_p: cobj,
):
    mbar = Ptr[float](mbar_p)
    fldcw = Ptr[float](fldcw_p)
    vmr = Ptr[float](vmr_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            vmr[_idx2(i, k, ncol)] = (
                mbar[_idx2(i, k, ncol)] * fldcw[_idx2(i, k, fldcw_ld1)] / adv_mass
            )


@export
def vmr2qqcw_codon(
    ncol: int,
    pver: int,
    fldcw_ld1: int,
    vmr_p: cobj,
    mbar_p: cobj,
    adv_mass: float,
    fldcw_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    mbar = Ptr[float](mbar_p)
    fldcw = Ptr[float](fldcw_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            fldcw[_idx2(i, k, fldcw_ld1)] = (
                adv_mass * vmr[_idx2(i, k, ncol)] / mbar[_idx2(i, k, ncol)]
            )


@export
def gas_aer_uptkrates_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    ntot_amode: int,
    q_p: cobj,
    t_p: cobj,
    pmid_p: cobj,
    dgncur_awet_p: cobj,
    numptr_p: cobj,
    sigmag_p: cobj,
    mwdry: float,
    rair: float,
    uptkrate_p: cobj,
):
    q = Ptr[float](q_p)
    t = Ptr[float](t_p)
    pmid = Ptr[float](pmid_p)
    dgncur_awet = Ptr[float](dgncur_awet_p)
    numptr = Ptr[int](numptr_p)
    sigmag = Ptr[float](sigmag_p)
    uptkrate = Ptr[float](uptkrate_p)

    tworootpi = 3.5449077
    root2 = 1.4142135
    beta = 2.0
    xghq0 = 0.70710678
    xghq1 = -0.70710678
    wghq0 = 0.88622693
    wghq1 = 0.88622693

    for n in range(1, ntot_amode + 1):
        lnsg = log(sigmag[n - 1])
        beta_lnsg_sq = beta * (lnsg**2.0)
        half_beta_lnsg_sq = 0.5 * ((beta * lnsg) ** 2.0)
        numptr_idx = numptr[n - 1]

        for k in range(top_lev, pver + 1):
            for i in range(1, ncol + 1):
                temp = t[_idx2(i, k, pcols)]
                pmid_ik = pmid[_idx2(i, k, pcols)]
                rhoair = pmid_ik / (rair * temp)
                aircon = rhoair / mwdry
                num_a = q[_idx3(i, k, numptr_idx, ncol, pver)] * aircon

                gasdiffus = 0.557e-4 * (temp**1.75) / pmid_ik
                gasspeed = 1.470e1 * sqrt(temp)
                freepathx2 = 6.0 * gasdiffus / gasspeed

                lndpgn = log(dgncur_awet[_idx3(i, k, n, pcols, pver)])
                const = tworootpi * num_a * exp(beta * lndpgn + half_beta_lnsg_sq)

                lndp = lndpgn + beta_lnsg_sq + root2 * lnsg * xghq0
                dp = exp(lndp)
                knudsen = freepathx2 / dp
                fuchs_sutugin = (0.4875 * (1.0 + knudsen)) / (
                    knudsen * (1.184 + knudsen) + 0.4875
                )
                sumghq = wghq0 * dp * fuchs_sutugin / (dp**beta)

                lndp = lndpgn + beta_lnsg_sq + root2 * lnsg * xghq1
                dp = exp(lndp)
                knudsen = freepathx2 / dp
                fuchs_sutugin = (0.4875 * (1.0 + knudsen)) / (
                    knudsen * (1.184 + knudsen) + 0.4875
                )
                sumghq += wghq1 * dp * fuchs_sutugin / (dp**beta)

                uptkrate[_idx3(n, i, k, ntot_amode, pcols)] = const * gasdiffus * sumghq


@export
def modal_aero_soaexch_codon(
    dtfull: float,
    temp: float,
    pres: float,
    niter_max: int,
    ntot_soamode: int,
    g_soa_in: float,
    a_soa_in_p: cobj,
    a_poa_in_p: cobj,
    xferrate_p: cobj,
    rgas: float,
    a_opoa_p: cobj,
    a_soa_p: cobj,
    beta_p: cobj,
    g_star_p: cobj,
    phi_p: cobj,
    sat_p: cobj,
    niter_p: cobj,
    g_soa_tend_p: cobj,
    a_soa_tend_p: cobj,
):
    a_soa_in = Ptr[float](a_soa_in_p)
    a_poa_in = Ptr[float](a_poa_in_p)
    xferrate = Ptr[float](xferrate_p)
    a_opoa = Ptr[float](a_opoa_p)
    a_soa = Ptr[float](a_soa_p)
    beta = Ptr[float](beta_p)
    g_star = Ptr[float](g_star_p)
    phi = Ptr[float](phi_p)
    sat = Ptr[float](sat_p)
    niter_out = Ptr[int](niter_p)
    g_soa_tend_out = Ptr[float](g_soa_tend_p)
    a_soa_tend = Ptr[float](a_soa_tend_p)

    alpha = 0.05
    g_min1 = 1.0e-20
    opoa_frac = 0.1
    delh_vap_soa = 156.0e3
    p0_soa_298 = 1.0e-10

    g_soa = g_soa_in
    if g_soa < 0.0:
        g_soa = 0.0
    tot_soa = g_soa

    for m in range(1, ntot_soamode + 1):
        a_soa_val = a_soa_in[m - 1]
        if a_soa_val < 0.0:
            a_soa_val = 0.0
        a_soa[m - 1] = a_soa_val
        tot_soa += a_soa_val

        a_opoa_val = opoa_frac * a_poa_in[m - 1]
        if a_opoa_val < 1.0e-20:
            a_opoa_val = 1.0e-20
        a_opoa[m - 1] = a_opoa_val

    p0_soa = p0_soa_298 * exp(
        -(delh_vap_soa / rgas) * ((1.0 / temp) - (1.0 / 298.0))
    )
    g0_soa = 1.01325e5 * p0_soa / pres
    g0_soa = g0_soa * (150.0 / 12.0)

    niter = 0
    tcur = 0.0
    dtcur = 0.0
    for m in range(1, ntot_soamode + 1):
        phi[m - 1] = 0.0
        g_star[m - 1] = 0.0

    while tcur < dtfull - 1.0e-3:
        niter += 1
        if niter > niter_max:
            break

        tmpa = 0.0
        for m in range(1, ntot_soamode + 1):
            sat[m - 1] = g0_soa / (a_soa[m - 1] + a_opoa[m - 1])
            g_star[m - 1] = sat[m - 1] * a_soa[m - 1]
            denom = _max3(g_soa, g_star[m - 1], g_min1)
            phi[m - 1] = (g_soa - g_star[m - 1]) / denom
            tmpa += xferrate[m - 1] * abs(phi[m - 1])

        dtmax = dtfull - tcur
        if dtmax * tmpa <= alpha:
            dtcur = dtmax
            tcur = dtfull
        else:
            dtcur = alpha / tmpa
            tcur += dtcur

        for m in range(1, ntot_soamode + 1):
            beta[m - 1] = dtcur * xferrate[m - 1]
            tmpa = g_soa - g_star[m - 1]
            if tmpa > 0.0:
                a_soa_tmp = a_soa[m - 1] + beta[m - 1] * tmpa
                sat[m - 1] = g0_soa / (a_soa_tmp + a_opoa[m - 1])
                g_star[m - 1] = sat[m - 1] * a_soa_tmp

        tmpa = 0.0
        tmpb = 0.0
        for m in range(1, ntot_soamode + 1):
            denom = 1.0 + beta[m - 1] * sat[m - 1]
            tmpa += a_soa[m - 1] / denom
            tmpb += beta[m - 1] / denom

        g_soa = (tot_soa - tmpa) / (1.0 + tmpb)
        if g_soa < 0.0:
            g_soa = 0.0
        for m in range(1, ntot_soamode + 1):
            a_soa[m - 1] = (a_soa[m - 1] + beta[m - 1] * g_soa) / (
                1.0 + beta[m - 1] * sat[m - 1]
            )

    g_soa_tend_out[0] = (g_soa - g_soa_in) / dtfull
    for m in range(1, ntot_soamode + 1):
        a_soa_tend[m - 1] = (a_soa[m - 1] - a_soa_in[m - 1]) / dtfull
    niter_out[0] = niter


@export
def modal_aero_rename_no_acc_crs_dryvols_codon(
    ncol: int,
    pver: int,
    pcnstxx: int,
    ntot_amode: int,
    maxspec_renamexf: int,
    loffset: int,
    deltat: float,
    idomode_p: cobj,
    q_p: cobj,
    qqcw_p: cobj,
    dqdt_p: cobj,
    dqdt_other_p: cobj,
    dqqcwdt_p: cobj,
    dqqcwdt_other_p: cobj,
    nspec_amode_p: cobj,
    lspectype_amode_p: cobj,
    specmw_amode_p: cobj,
    specdens_amode_p: cobj,
    lmassptr_amode_p: cobj,
    lmassptrcw_amode_p: cobj,
    dryvol_a_p: cobj,
    dryvol_c_p: cobj,
    deldryvol_a_p: cobj,
    deldryvol_c_p: cobj,
):
    idomode = Ptr[int](idomode_p)
    q = Ptr[float](q_p)
    qqcw = Ptr[float](qqcw_p)
    dqdt = Ptr[float](dqdt_p)
    dqdt_other = Ptr[float](dqdt_other_p)
    dqqcwdt = Ptr[float](dqqcwdt_p)
    dqqcwdt_other = Ptr[float](dqqcwdt_other_p)
    nspec_amode = Ptr[int](nspec_amode_p)
    lspectype_amode = Ptr[int](lspectype_amode_p)
    specmw_amode = Ptr[float](specmw_amode_p)
    specdens_amode = Ptr[float](specdens_amode_p)
    lmassptr_amode = Ptr[int](lmassptr_amode_p)
    lmassptrcw_amode = Ptr[int](lmassptrcw_amode_p)
    dryvol_a = Ptr[float](dryvol_a_p)
    dryvol_c = Ptr[float](dryvol_c_p)
    deldryvol_a = Ptr[float](deldryvol_a_p)
    deldryvol_c = Ptr[float](deldryvol_c_p)

    for n in range(1, ntot_amode + 1):
        if idomode[n - 1] > 0:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    dryvol_a[_idx3(i, k, n, ncol, pver)] = 0.0
                    dryvol_c[_idx3(i, k, n, ncol, pver)] = 0.0
                    deldryvol_a[_idx3(i, k, n, ncol, pver)] = 0.0
                    deldryvol_c[_idx3(i, k, n, ncol, pver)] = 0.0

            for l1 in range(1, nspec_amode[n - 1] + 1):
                l2 = lspectype_amode[_idx2(l1, n, maxspec_renamexf)]
                dum_m2v = specmw_amode[l2 - 1] / specdens_amode[l2 - 1]
                dum_m2vdt = dum_m2v * deltat

                la = lmassptr_amode[_idx2(l1, n, maxspec_renamexf)] - loffset
                if la > 0:
                    for k in range(1, pver + 1):
                        for i in range(1, ncol + 1):
                            qold = q[_idx3(i, k, la, ncol, pver)] - deltat * dqdt_other[
                                _idx3(i, k, la, ncol, pver)
                            ]
                            if qold < 0.0:
                                qold = 0.0
                            dryvol_a[_idx3(i, k, n, ncol, pver)] += dum_m2v * qold
                            deldryvol_a[_idx3(i, k, n, ncol, pver)] += (
                                dqdt_other[_idx3(i, k, la, ncol, pver)]
                                + dqdt[_idx3(i, k, la, ncol, pver)]
                            ) * dum_m2vdt

                lc = lmassptrcw_amode[_idx2(l1, n, maxspec_renamexf)] - loffset
                if lc > 0:
                    for k in range(1, pver + 1):
                        for i in range(1, ncol + 1):
                            qqcwold = qqcw[_idx3(i, k, lc, ncol, pver)] - deltat * dqqcwdt_other[
                                _idx3(i, k, lc, ncol, pver)
                            ]
                            if qqcwold < 0.0:
                                qqcwold = 0.0
                            dryvol_c[_idx3(i, k, n, ncol, pver)] += dum_m2v * qqcwold
                            deldryvol_c[_idx3(i, k, n, ncol, pver)] += (
                                dqqcwdt_other[_idx3(i, k, lc, ncol, pver)]
                                + dqqcwdt[_idx3(i, k, lc, ncol, pver)]
                            ) * dum_m2vdt


@export
def modal_aero_rename_no_acc_crs_xferfracs_codon(
    ncol: int,
    pver: int,
    pcnstxx: int,
    ntot_amode: int,
    maxpair_renamexf: int,
    loffset: int,
    npair_renamexf: int,
    q_p: cobj,
    qqcw_p: cobj,
    dryvol_a_p: cobj,
    dryvol_c_p: cobj,
    deldryvol_a_p: cobj,
    deldryvol_c_p: cobj,
    modefrm_renamexf_p: cobj,
    modetoo_renamexf_p: cobj,
    numptr_amode_p: cobj,
    numptrcw_amode_p: cobj,
    dgnum_amode_p: cobj,
    factoraa_p: cobj,
    factoryy_p: cobj,
    dryvol_smallest_p: cobj,
    v2nlorlx_p: cobj,
    v2nhirlx_p: cobj,
    dum3alnsg2_p: cobj,
    dp_cut_p: cobj,
    lndp_cut_p: cobj,
    dp_belowcut_p: cobj,
    onethird: float,
    xferfrac_max: float,
    xferfrac_vol_p: cobj,
    xferfrac_num_p: cobj,
):
    q = Ptr[float](q_p)
    qqcw = Ptr[float](qqcw_p)
    dryvol_a = Ptr[float](dryvol_a_p)
    dryvol_c = Ptr[float](dryvol_c_p)
    deldryvol_a = Ptr[float](deldryvol_a_p)
    deldryvol_c = Ptr[float](deldryvol_c_p)
    modefrm_renamexf = Ptr[int](modefrm_renamexf_p)
    modetoo_renamexf = Ptr[int](modetoo_renamexf_p)
    numptr_amode = Ptr[int](numptr_amode_p)
    numptrcw_amode = Ptr[int](numptrcw_amode_p)
    dgnum_amode = Ptr[float](dgnum_amode_p)
    factoraa = Ptr[float](factoraa_p)
    factoryy = Ptr[float](factoryy_p)
    dryvol_smallest = Ptr[float](dryvol_smallest_p)
    v2nlorlx = Ptr[float](v2nlorlx_p)
    v2nhirlx = Ptr[float](v2nhirlx_p)
    dum3alnsg2 = Ptr[float](dum3alnsg2_p)
    dp_cut = Ptr[float](dp_cut_p)
    lndp_cut = Ptr[float](lndp_cut_p)
    dp_belowcut = Ptr[float](dp_belowcut_p)
    xferfrac_vol = Ptr[float](xferfrac_vol_p)
    xferfrac_num = Ptr[float](xferfrac_num_p)

    for ipair in range(1, maxpair_renamexf + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                xferfrac_vol[_idx3(i, k, ipair, ncol, pver)] = 0.0
                xferfrac_num[_idx3(i, k, ipair, ncol, pver)] = 0.0

    for ipair in range(1, npair_renamexf + 1):
        mfrm = modefrm_renamexf[ipair - 1]
        mtoo = modetoo_renamexf[ipair - 1]

        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                dryvol_t_old = dryvol_a[_idx3(i, k, mfrm, ncol, pver)] + dryvol_c[
                    _idx3(i, k, mfrm, ncol, pver)
                ]
                dryvol_t_del = deldryvol_a[_idx3(i, k, mfrm, ncol, pver)] + deldryvol_c[
                    _idx3(i, k, mfrm, ncol, pver)
                ]
                dryvol_t_new = dryvol_t_old + dryvol_t_del
                dryvol_t_oldbnd = max(dryvol_t_old, dryvol_smallest[mfrm - 1])

                if dryvol_t_new <= dryvol_smallest[mfrm - 1]:
                    continue
                if dryvol_t_del <= 1.0e-6 * dryvol_t_oldbnd:
                    continue

                num_t_old = q[_idx3(i, k, numptr_amode[mfrm - 1] - loffset, ncol, pver)]
                num_t_old += qqcw[
                    _idx3(i, k, numptrcw_amode[mfrm - 1] - loffset, ncol, pver)
                ]
                if num_t_old < 0.0:
                    num_t_old = 0.0

                dryvol_t_oldbnd = max(dryvol_t_old, dryvol_smallest[mfrm - 1])
                num_t_oldbnd = min(dryvol_t_oldbnd * v2nlorlx[mfrm - 1], num_t_old)
                num_t_oldbnd = max(dryvol_t_oldbnd * v2nhirlx[mfrm - 1], num_t_oldbnd)

                dgn_t_new = (dryvol_t_new / (num_t_oldbnd * factoraa[mfrm - 1])) ** onethird
                if dgn_t_new <= dgnum_amode[mfrm - 1]:
                    continue

                lndgn_new = log(dgn_t_new)
                lndgv_new = lndgn_new + dum3alnsg2[ipair - 1]
                yn_tail = (lndp_cut[ipair - 1] - lndgn_new) * factoryy[mfrm - 1]
                yv_tail = (lndp_cut[ipair - 1] - lndgv_new) * factoryy[mfrm - 1]
                tailfr_numnew = 0.5 * erfc(yn_tail)
                tailfr_volnew = 0.5 * erfc(yv_tail)

                dgn_t_old = (dryvol_t_oldbnd / (num_t_oldbnd * factoraa[mfrm - 1])) ** onethird
                if dgn_t_new >= dp_cut[ipair - 1]:
                    dgn_t_old = min(dgn_t_old, dp_belowcut[ipair - 1])

                lndgn_old = log(dgn_t_old)
                lndgv_old = lndgn_old + dum3alnsg2[ipair - 1]
                yn_tail = (lndp_cut[ipair - 1] - lndgn_old) * factoryy[mfrm - 1]
                yv_tail = (lndp_cut[ipair - 1] - lndgv_old) * factoryy[mfrm - 1]
                tailfr_numold = 0.5 * erfc(yn_tail)
                tailfr_volold = 0.5 * erfc(yv_tail)

                dum = tailfr_volnew * dryvol_t_new - tailfr_volold * dryvol_t_old
                if dum <= 0.0:
                    continue

                xferfrac_vol_val = min(dum, dryvol_t_new) / dryvol_t_new
                xferfrac_vol_val = min(xferfrac_vol_val, xferfrac_max)
                xferfrac_num_val = tailfr_numnew - tailfr_numold
                xferfrac_num_val = max(0.0, min(xferfrac_num_val, xferfrac_vol_val))

                xferfrac_vol[_idx3(i, k, ipair, ncol, pver)] = xferfrac_vol_val
                xferfrac_num[_idx3(i, k, ipair, ncol, pver)] = xferfrac_num_val


@export
def modal_aero_rename_no_acc_crs_tendencies_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnstxx: int,
    maxpair_renamexf: int,
    maxspec_renamexf: int,
    loffset: int,
    npair_renamexf: int,
    jsrflx_rename: int,
    nsrflx: int,
    is_dorename_atik: int,
    deltat: float,
    deltatinv: float,
    gravit: float,
    pdel_p: cobj,
    dorename_atik_p: cobj,
    q_p: cobj,
    qqcw_p: cobj,
    dqdt_p: cobj,
    dqqcwdt_p: cobj,
    qsrflx_p: cobj,
    qqcwsrflx_p: cobj,
    xferfrac_vol_p: cobj,
    xferfrac_num_p: cobj,
    nspecfrm_renamexf_p: cobj,
    lspecfrma_renamexf_p: cobj,
    lspecfrmc_renamexf_p: cobj,
    lspectooa_renamexf_p: cobj,
    lspectooc_renamexf_p: cobj,
):
    pdel = Ptr[float](pdel_p)
    dorename_atik = Ptr[int](dorename_atik_p)
    q = Ptr[float](q_p)
    qqcw = Ptr[float](qqcw_p)
    dqdt = Ptr[float](dqdt_p)
    dqqcwdt = Ptr[float](dqqcwdt_p)
    qsrflx = Ptr[float](qsrflx_p)
    qqcwsrflx = Ptr[float](qqcwsrflx_p)
    xferfrac_vol = Ptr[float](xferfrac_vol_p)
    xferfrac_num = Ptr[float](xferfrac_num_p)
    nspecfrm_renamexf = Ptr[int](nspecfrm_renamexf_p)
    lspecfrma_renamexf = Ptr[int](lspecfrma_renamexf_p)
    lspecfrmc_renamexf = Ptr[int](lspecfrmc_renamexf_p)
    lspectooa_renamexf = Ptr[int](lspectooa_renamexf_p)
    lspectooc_renamexf = Ptr[int](lspectooc_renamexf_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            if is_dorename_atik != 0:
                if dorename_atik[_idx2(i, k, ncol)] == 0:
                    continue

            pdel_fac = pdel[_idx2(i, k, pcols)] / gravit

            for ipair in range(1, npair_renamexf + 1):
                xferfrac_vol_local = xferfrac_vol[_idx3(i, k, ipair, ncol, pver)]
                xferfrac_num_local = xferfrac_num[_idx3(i, k, ipair, ncol, pver)]
                if xferfrac_vol_local <= 0.0:
                    continue

                for iq in range(1, nspecfrm_renamexf[ipair - 1] + 1):
                    xfercoef = xferfrac_vol_local * deltatinv
                    if iq == 1:
                        xfercoef = xferfrac_num_local * deltatinv

                    lsfrma = lspecfrma_renamexf[_idx2(iq, ipair, maxspec_renamexf)] - loffset
                    lsfrmc = lspecfrmc_renamexf[_idx2(iq, ipair, maxspec_renamexf)] - loffset
                    lstooa = lspectooa_renamexf[_idx2(iq, ipair, maxspec_renamexf)] - loffset
                    lstooc = lspectooc_renamexf[_idx2(iq, ipair, maxspec_renamexf)] - loffset

                    if lsfrma > 0:
                        xfertend = xfercoef * max(
                            0.0,
                            q[_idx3(i, k, lsfrma, ncol, pver)]
                            + dqdt[_idx3(i, k, lsfrma, ncol, pver)] * deltat,
                        )
                        dqdt[_idx3(i, k, lsfrma, ncol, pver)] -= xfertend
                        qsrflx[_idx3(i, lsfrma, jsrflx_rename, pcols, pcnstxx)] -= (
                            xfertend * pdel_fac
                        )
                        if lstooa > 0:
                            dqdt[_idx3(i, k, lstooa, ncol, pver)] += xfertend
                            qsrflx[_idx3(i, lstooa, jsrflx_rename, pcols, pcnstxx)] += (
                                xfertend * pdel_fac
                            )

                    if lsfrmc > 0:
                        xfertend = xfercoef * max(
                            0.0,
                            qqcw[_idx3(i, k, lsfrmc, ncol, pver)]
                            + dqqcwdt[_idx3(i, k, lsfrmc, ncol, pver)] * deltat,
                        )
                        dqqcwdt[_idx3(i, k, lsfrmc, ncol, pver)] -= xfertend
                        qqcwsrflx[_idx3(i, lsfrmc, jsrflx_rename, pcols, pcnstxx)] -= (
                            xfertend * pdel_fac
                        )
                        if lstooc > 0:
                            dqqcwdt[_idx3(i, k, lstooc, ncol, pver)] += xfertend
                            qqcwsrflx[_idx3(i, lstooc, jsrflx_rename, pcols, pcnstxx)] += (
                                xfertend * pdel_fac
                            )


@export
def modal_aero_rename_no_acc_crs_sub_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnstxx: int,
    ntot_amode: int,
    maxpair_renamexf: int,
    maxspec_renamexf: int,
    loffset: int,
    npair_renamexf: int,
    is_dorename_atik: int,
    jsrflx_rename: int,
    nsrflx: int,
    deltat: float,
    deltatinv: float,
    onethird: float,
    xferfrac_max: float,
    pi: float,
    gravit: float,
    pdel_p: cobj,
    dorename_atik_p: cobj,
    q_p: cobj,
    qqcw_p: cobj,
    dqdt_p: cobj,
    dqdt_other_p: cobj,
    dqqcwdt_p: cobj,
    dqqcwdt_other_p: cobj,
    qsrflx_p: cobj,
    qqcwsrflx_p: cobj,
    modefrm_renamexf_p: cobj,
    modetoo_renamexf_p: cobj,
    nspec_amode_p: cobj,
    lspectype_amode_p: cobj,
    specmw_amode_p: cobj,
    specdens_amode_p: cobj,
    lmassptr_amode_p: cobj,
    lmassptrcw_amode_p: cobj,
    numptr_amode_p: cobj,
    numptrcw_amode_p: cobj,
    dgnum_amode_p: cobj,
    alnsg_amode_p: cobj,
    voltonumblo_amode_p: cobj,
    voltonumbhi_amode_p: cobj,
    nspecfrm_renamexf_p: cobj,
    lspecfrma_renamexf_p: cobj,
    lspecfrmc_renamexf_p: cobj,
    lspectooa_renamexf_p: cobj,
    lspectooc_renamexf_p: cobj,
    idomode_p: cobj,
    dryvol_a_p: cobj,
    dryvol_c_p: cobj,
    deldryvol_a_p: cobj,
    deldryvol_c_p: cobj,
    factoraa_p: cobj,
    factoryy_p: cobj,
    dryvol_smallest_p: cobj,
    v2nlorlx_p: cobj,
    v2nhirlx_p: cobj,
    dum3alnsg2_p: cobj,
    dp_cut_p: cobj,
    lndp_cut_p: cobj,
    dp_belowcut_p: cobj,
    xferfrac_vol_p: cobj,
    xferfrac_num_p: cobj,
    dotendrn_p: cobj,
    dotendqqcwrn_p: cobj,
):
    idomode = Ptr[int](idomode_p)
    modefrm_renamexf = Ptr[int](modefrm_renamexf_p)
    modetoo_renamexf = Ptr[int](modetoo_renamexf_p)
    nspec_amode = Ptr[int](nspec_amode_p)
    numptr_amode = Ptr[int](numptr_amode_p)
    numptrcw_amode = Ptr[int](numptrcw_amode_p)
    dgnum_amode = Ptr[float](dgnum_amode_p)
    alnsg_amode = Ptr[float](alnsg_amode_p)
    voltonumblo_amode = Ptr[float](voltonumblo_amode_p)
    voltonumbhi_amode = Ptr[float](voltonumbhi_amode_p)
    factoraa = Ptr[float](factoraa_p)
    factoryy = Ptr[float](factoryy_p)
    dryvol_smallest = Ptr[float](dryvol_smallest_p)
    v2nlorlx = Ptr[float](v2nlorlx_p)
    v2nhirlx = Ptr[float](v2nhirlx_p)
    dum3alnsg2 = Ptr[float](dum3alnsg2_p)
    dp_cut = Ptr[float](dp_cut_p)
    lndp_cut = Ptr[float](lndp_cut_p)
    dp_belowcut = Ptr[float](dp_belowcut_p)

    frelax = 27.0

    for n in range(1, ntot_amode + 1):
        idomode[n - 1] = 0

    for ipair in range(1, npair_renamexf + 1):
        mfrm = modefrm_renamexf[ipair - 1]
        mtoo = modetoo_renamexf[ipair - 1]
        idomode[mfrm - 1] = 1

        factoraa[mfrm - 1] = (pi / 6.0) * exp(4.5 * (alnsg_amode[mfrm - 1] ** 2))
        factoraa[mtoo - 1] = (pi / 6.0) * exp(4.5 * (alnsg_amode[mtoo - 1] ** 2))
        factoryy[mfrm - 1] = sqrt(0.5) / alnsg_amode[mfrm - 1]
        dryvol_smallest[mfrm - 1] = 1.0e-25
        v2nlorlx[mfrm - 1] = voltonumblo_amode[mfrm - 1] * frelax
        v2nhirlx[mfrm - 1] = voltonumbhi_amode[mfrm - 1] / frelax

        dum3alnsg2[ipair - 1] = 3.0 * (alnsg_amode[mfrm - 1] ** 2)
        dp_cut[ipair - 1] = sqrt(
            dgnum_amode[mfrm - 1] * exp(1.5 * (alnsg_amode[mfrm - 1] ** 2))
            * dgnum_amode[mtoo - 1]
            * exp(1.5 * (alnsg_amode[mtoo - 1] ** 2))
        )
        lndp_cut[ipair - 1] = log(dp_cut[ipair - 1])
        dp_belowcut[ipair - 1] = 0.99 * dp_cut[ipair - 1]

    modal_aero_rename_no_acc_crs_dryvols_codon(
        ncol,
        pver,
        pcnstxx,
        ntot_amode,
        maxspec_renamexf,
        loffset,
        deltat,
        idomode_p,
        q_p,
        qqcw_p,
        dqdt_p,
        dqdt_other_p,
        dqqcwdt_p,
        dqqcwdt_other_p,
        nspec_amode_p,
        lspectype_amode_p,
        specmw_amode_p,
        specdens_amode_p,
        lmassptr_amode_p,
        lmassptrcw_amode_p,
        dryvol_a_p,
        dryvol_c_p,
        deldryvol_a_p,
        deldryvol_c_p,
    )

    modal_aero_rename_no_acc_crs_xferfracs_codon(
        ncol,
        pver,
        pcnstxx,
        ntot_amode,
        maxpair_renamexf,
        loffset,
        npair_renamexf,
        q_p,
        qqcw_p,
        dryvol_a_p,
        dryvol_c_p,
        deldryvol_a_p,
        deldryvol_c_p,
        modefrm_renamexf_p,
        modetoo_renamexf_p,
        numptr_amode_p,
        numptrcw_amode_p,
        dgnum_amode_p,
        factoraa_p,
        factoryy_p,
        dryvol_smallest_p,
        v2nlorlx_p,
        v2nhirlx_p,
        dum3alnsg2_p,
        dp_cut_p,
        lndp_cut_p,
        dp_belowcut_p,
        onethird,
        xferfrac_max,
        xferfrac_vol_p,
        xferfrac_num_p,
    )

    modal_aero_rename_no_acc_crs_tendencies_codon(
        ncol,
        pcols,
        pver,
        pcnstxx,
        maxpair_renamexf,
        maxspec_renamexf,
        loffset,
        npair_renamexf,
        jsrflx_rename,
        nsrflx,
        is_dorename_atik,
        deltat,
        deltatinv,
        gravit,
        pdel_p,
        dorename_atik_p,
        q_p,
        qqcw_p,
        dqdt_p,
        dqqcwdt_p,
        qsrflx_p,
        qqcwsrflx_p,
        xferfrac_vol_p,
        xferfrac_num_p,
        nspecfrm_renamexf_p,
        lspecfrma_renamexf_p,
        lspecfrmc_renamexf_p,
        lspectooa_renamexf_p,
        lspectooc_renamexf_p,
    )

    modal_aero_rename_set_dotend_flags_codon(
        pcnstxx,
        maxpair_renamexf,
        maxspec_renamexf,
        loffset,
        npair_renamexf,
        nspecfrm_renamexf_p,
        lspecfrma_renamexf_p,
        lspecfrmc_renamexf_p,
        lspectooa_renamexf_p,
        lspectooc_renamexf_p,
        dotendrn_p,
        dotendqqcwrn_p,
    )


@export
def modal_aero_rename_acc_crs_dryvols_codon(
    ncol: int,
    pver: int,
    pcnstxx: int,
    maxspec_renamexf: int,
    loffset: int,
    ixferable_all: int,
    nspec_mfrm: int,
    deltat: float,
    q_p: cobj,
    qqcw_p: cobj,
    dqdt_p: cobj,
    dqdt_other_p: cobj,
    dqqcwdt_p: cobj,
    dqqcwdt_other_p: cobj,
    lspectype_mfrm_p: cobj,
    specmw_amode_p: cobj,
    specdens_amode_p: cobj,
    lmassptr_mfrm_p: cobj,
    lmassptrcw_mfrm_p: cobj,
    ixferable_a_p: cobj,
    ixferable_c_p: cobj,
    dryvol_a_p: cobj,
    dryvol_c_p: cobj,
    deldryvol_a_p: cobj,
    deldryvol_c_p: cobj,
    dryvol_a_xfab_p: cobj,
    dryvol_c_xfab_p: cobj,
):
    q = Ptr[float](q_p)
    qqcw = Ptr[float](qqcw_p)
    dqdt = Ptr[float](dqdt_p)
    dqdt_other = Ptr[float](dqdt_other_p)
    dqqcwdt = Ptr[float](dqqcwdt_p)
    dqqcwdt_other = Ptr[float](dqqcwdt_other_p)
    lspectype_mfrm = Ptr[int](lspectype_mfrm_p)
    specmw_amode = Ptr[float](specmw_amode_p)
    specdens_amode = Ptr[float](specdens_amode_p)
    lmassptr_mfrm = Ptr[int](lmassptr_mfrm_p)
    lmassptrcw_mfrm = Ptr[int](lmassptrcw_mfrm_p)
    ixferable_a = Ptr[int](ixferable_a_p)
    ixferable_c = Ptr[int](ixferable_c_p)
    dryvol_a = Ptr[float](dryvol_a_p)
    dryvol_c = Ptr[float](dryvol_c_p)
    deldryvol_a = Ptr[float](deldryvol_a_p)
    deldryvol_c = Ptr[float](deldryvol_c_p)
    dryvol_a_xfab = Ptr[float](dryvol_a_xfab_p)
    dryvol_c_xfab = Ptr[float](dryvol_c_xfab_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            dryvol_a[_idx2(i, k, ncol)] = 0.0
            dryvol_c[_idx2(i, k, ncol)] = 0.0
            deldryvol_a[_idx2(i, k, ncol)] = 0.0
            deldryvol_c[_idx2(i, k, ncol)] = 0.0
            dryvol_a_xfab[_idx2(i, k, ncol)] = 0.0
            dryvol_c_xfab[_idx2(i, k, ncol)] = 0.0

    for l1 in range(1, nspec_mfrm + 1):
        l2 = lspectype_mfrm[l1 - 1]
        tmp_m2v = specmw_amode[l2 - 1] / specdens_amode[l2 - 1]
        tmp_m2vdt = tmp_m2v * deltat

        la = lmassptr_mfrm[l1 - 1] - loffset
        if la > 0:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    dryvol_a[_idx2(i, k, ncol)] += tmp_m2v * max(
                        0.0,
                        q[_idx3(i, k, la, ncol, pver)]
                        - deltat * dqdt_other[_idx3(i, k, la, ncol, pver)],
                    )
                    deldryvol_a[_idx2(i, k, ncol)] += (
                        dqdt_other[_idx3(i, k, la, ncol, pver)]
                        + dqdt[_idx3(i, k, la, ncol, pver)]
                    ) * tmp_m2vdt
                    if ixferable_all <= 0 and ixferable_a[l1 - 1] > 0:
                        dryvol_a_xfab[_idx2(i, k, ncol)] += tmp_m2v * max(
                            0.0,
                            q[_idx3(i, k, la, ncol, pver)]
                            + deltat * dqdt[_idx3(i, k, la, ncol, pver)],
                        )

        lc = lmassptrcw_mfrm[l1 - 1] - loffset
        if lc > 0:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    dryvol_c[_idx2(i, k, ncol)] += tmp_m2v * max(
                        0.0,
                        qqcw[_idx3(i, k, lc, ncol, pver)]
                        - deltat * dqqcwdt_other[_idx3(i, k, lc, ncol, pver)],
                    )
                    deldryvol_c[_idx2(i, k, ncol)] += (
                        dqqcwdt_other[_idx3(i, k, lc, ncol, pver)]
                        + dqqcwdt[_idx3(i, k, lc, ncol, pver)]
                    ) * tmp_m2vdt
                    if ixferable_all <= 0 and ixferable_c[l1 - 1] > 0:
                        dryvol_c_xfab[_idx2(i, k, ncol)] += tmp_m2v * max(
                            0.0,
                            qqcw[_idx3(i, k, lc, ncol, pver)]
                            + deltat * dqqcwdt[_idx3(i, k, lc, ncol, pver)],
                        )


@export
def modal_aero_rename_acc_crs_tendencies_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnstxx: int,
    maxspec_renamexf: int,
    loffset: int,
    nspecfrm_renamexf: int,
    jsrflx_rename: int,
    nsrflx: int,
    is_dorename_atik: int,
    l_dqdt_rnpos: int,
    deltat: float,
    deltatinv: float,
    gravit: float,
    pdel_p: cobj,
    dorename_atik_p: cobj,
    q_p: cobj,
    qqcw_p: cobj,
    dqdt_p: cobj,
    dqqcwdt_p: cobj,
    qsrflx_p: cobj,
    qqcwsrflx_p: cobj,
    xferfrac_vol_p: cobj,
    xferfrac_num_p: cobj,
    lspecfrma_renamexf_p: cobj,
    lspecfrmc_renamexf_p: cobj,
    lspectooa_renamexf_p: cobj,
    lspectooc_renamexf_p: cobj,
    dqdt_rnpos_p: cobj,
):
    pdel = Ptr[float](pdel_p)
    dorename_atik = Ptr[int](dorename_atik_p)
    q = Ptr[float](q_p)
    qqcw = Ptr[float](qqcw_p)
    dqdt = Ptr[float](dqdt_p)
    dqqcwdt = Ptr[float](dqqcwdt_p)
    qsrflx = Ptr[float](qsrflx_p)
    qqcwsrflx = Ptr[float](qqcwsrflx_p)
    xferfrac_vol = Ptr[float](xferfrac_vol_p)
    xferfrac_num = Ptr[float](xferfrac_num_p)
    lspecfrma_renamexf = Ptr[int](lspecfrma_renamexf_p)
    lspecfrmc_renamexf = Ptr[int](lspecfrmc_renamexf_p)
    lspectooa_renamexf = Ptr[int](lspectooa_renamexf_p)
    lspectooc_renamexf = Ptr[int](lspectooc_renamexf_p)
    dqdt_rnpos = Ptr[float](dqdt_rnpos_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            if is_dorename_atik != 0:
                if dorename_atik[_idx2(i, k, ncol)] == 0:
                    continue

            xferfrac_vol_local = xferfrac_vol[_idx2(i, k, ncol)]
            xferfrac_num_local = xferfrac_num[_idx2(i, k, ncol)]
            if xferfrac_vol_local <= 0.0:
                continue

            pdel_fac = pdel[_idx2(i, k, pcols)] / gravit

            for iq in range(1, nspecfrm_renamexf + 1):
                xfercoef = xferfrac_vol_local * deltatinv
                if iq == 1:
                    xfercoef = xferfrac_num_local * deltatinv

                lsfrma = lspecfrma_renamexf[iq - 1] - loffset
                lsfrmc = lspecfrmc_renamexf[iq - 1] - loffset
                lstooa = lspectooa_renamexf[iq - 1] - loffset
                lstooc = lspectooc_renamexf[iq - 1] - loffset

                if lsfrma > 0:
                    xfertend = xfercoef * max(
                        0.0,
                        q[_idx3(i, k, lsfrma, ncol, pver)]
                        + dqdt[_idx3(i, k, lsfrma, ncol, pver)] * deltat,
                    )
                    dqdt[_idx3(i, k, lsfrma, ncol, pver)] -= xfertend
                    qsrflx[_idx3(i, lsfrma, jsrflx_rename, pcols, pcnstxx)] -= (
                        xfertend * pdel_fac
                    )
                    if lstooa > 0:
                        dqdt[_idx3(i, k, lstooa, ncol, pver)] += xfertend
                        qsrflx[_idx3(i, lstooa, jsrflx_rename, pcols, pcnstxx)] += (
                            xfertend * pdel_fac
                        )
                        if l_dqdt_rnpos != 0:
                            dqdt_rnpos[_idx3(i, k, lstooa, ncol, pver)] += xfertend

                if lsfrmc > 0:
                    xfertend = xfercoef * max(
                        0.0,
                        qqcw[_idx3(i, k, lsfrmc, ncol, pver)]
                        + dqqcwdt[_idx3(i, k, lsfrmc, ncol, pver)] * deltat,
                    )
                    dqqcwdt[_idx3(i, k, lsfrmc, ncol, pver)] -= xfertend
                    qqcwsrflx[_idx3(i, lsfrmc, jsrflx_rename, pcols, pcnstxx)] -= (
                        xfertend * pdel_fac
                    )
                    if lstooc > 0:
                        dqqcwdt[_idx3(i, k, lstooc, ncol, pver)] += xfertend
                        qqcwsrflx[_idx3(i, lstooc, jsrflx_rename, pcols, pcnstxx)] += (
                            xfertend * pdel_fac
                        )


@export
def modal_aero_rename_acc_crs_pair_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnstxx: int,
    maxspec_renamexf: int,
    loffset: int,
    is_dorename_atik: int,
    l_dqdt_rnpos: int,
    jsrflx_rename: int,
    nsrflx: int,
    ixferable_all: int,
    nspec_mfrm: int,
    mfrm: int,
    numptr_amode_mfrm: int,
    numptrcw_amode_mfrm: int,
    igrow_shrink: int,
    method_optbb: int,
    flagaa_shrink: int,
    nspecfrm_renamexf: int,
    deltat: float,
    deltatinv: float,
    gravit: float,
    dgnum_amode_mfrm: float,
    factoraa: float,
    factoryy: float,
    dryvol_smallest: float,
    v2nlorlx: float,
    v2nhirlx: float,
    factor_3alnsg2: float,
    dp_cut: float,
    lndp_cut: float,
    dp_belowcut: float,
    dp_xfernone_thresh: float,
    dp_xferall_thresh: float,
    onethird: float,
    xferfrac_max: float,
    troplev_p: cobj,
    pdel_p: cobj,
    dorename_atik_p: cobj,
    q_p: cobj,
    qqcw_p: cobj,
    dqdt_p: cobj,
    dqdt_other_p: cobj,
    dqqcwdt_p: cobj,
    dqqcwdt_other_p: cobj,
    qsrflx_p: cobj,
    qqcwsrflx_p: cobj,
    lspectype_mfrm_p: cobj,
    specmw_amode_p: cobj,
    specdens_amode_p: cobj,
    lmassptr_mfrm_p: cobj,
    lmassptrcw_mfrm_p: cobj,
    ixferable_a_p: cobj,
    ixferable_c_p: cobj,
    lspecfrma_renamexf_p: cobj,
    lspecfrmc_renamexf_p: cobj,
    lspectooa_renamexf_p: cobj,
    lspectooc_renamexf_p: cobj,
    dryvol_a_p: cobj,
    dryvol_c_p: cobj,
    deldryvol_a_p: cobj,
    deldryvol_c_p: cobj,
    dryvol_a_xfab_p: cobj,
    dryvol_c_xfab_p: cobj,
    xferfrac_vol_p: cobj,
    xferfrac_num_p: cobj,
    dqdt_rnpos_p: cobj,
):
    modal_aero_rename_acc_crs_dryvols_codon(
        ncol,
        pver,
        pcnstxx,
        maxspec_renamexf,
        loffset,
        ixferable_all,
        nspec_mfrm,
        deltat,
        q_p,
        qqcw_p,
        dqdt_p,
        dqdt_other_p,
        dqqcwdt_p,
        dqqcwdt_other_p,
        lspectype_mfrm_p,
        specmw_amode_p,
        specdens_amode_p,
        lmassptr_mfrm_p,
        lmassptrcw_mfrm_p,
        ixferable_a_p,
        ixferable_c_p,
        dryvol_a_p,
        dryvol_c_p,
        deldryvol_a_p,
        deldryvol_c_p,
        dryvol_a_xfab_p,
        dryvol_c_xfab_p,
    )

    modal_aero_rename_acc_crs_xferfracs_codon(
        ncol,
        pcols,
        pver,
        pcnstxx,
        loffset,
        mfrm,
        numptr_amode_mfrm,
        numptrcw_amode_mfrm,
        igrow_shrink,
        ixferable_all,
        method_optbb,
        flagaa_shrink,
        dgnum_amode_mfrm,
        factoraa,
        factoryy,
        dryvol_smallest,
        v2nlorlx,
        v2nhirlx,
        factor_3alnsg2,
        dp_cut,
        lndp_cut,
        dp_belowcut,
        dp_xfernone_thresh,
        dp_xferall_thresh,
        onethird,
        xferfrac_max,
        troplev_p,
        q_p,
        qqcw_p,
        dryvol_a_p,
        dryvol_c_p,
        deldryvol_a_p,
        deldryvol_c_p,
        dryvol_a_xfab_p,
        dryvol_c_xfab_p,
        xferfrac_vol_p,
        xferfrac_num_p,
    )

    modal_aero_rename_acc_crs_tendencies_codon(
        ncol,
        pcols,
        pver,
        pcnstxx,
        maxspec_renamexf,
        loffset,
        nspecfrm_renamexf,
        jsrflx_rename,
        nsrflx,
        is_dorename_atik,
        l_dqdt_rnpos,
        deltat,
        deltatinv,
        gravit,
        pdel_p,
        dorename_atik_p,
        q_p,
        qqcw_p,
        dqdt_p,
        dqqcwdt_p,
        qsrflx_p,
        qqcwsrflx_p,
        xferfrac_vol_p,
        xferfrac_num_p,
        lspecfrma_renamexf_p,
        lspecfrmc_renamexf_p,
        lspectooa_renamexf_p,
        lspectooc_renamexf_p,
        dqdt_rnpos_p,
    )


def modal_aero_rename_acc_crs_dryvols_full_codon(
    ncol: int,
    pver: int,
    maxspec_renamexf: int,
    loffset: int,
    ipair: int,
    mfrm: int,
    nspec_mfrm: int,
    ixferable_all: int,
    deltat: float,
    q_p: cobj,
    qqcw_p: cobj,
    dqdt_p: cobj,
    dqdt_other_p: cobj,
    dqqcwdt_p: cobj,
    dqqcwdt_other_p: cobj,
    lspectype_amode_p: cobj,
    specmw_amode_p: cobj,
    specdens_amode_p: cobj,
    lmassptr_amode_p: cobj,
    lmassptrcw_amode_p: cobj,
    ixferable_a_renamexf_p: cobj,
    ixferable_c_renamexf_p: cobj,
    dryvol_a_p: cobj,
    dryvol_c_p: cobj,
    deldryvol_a_p: cobj,
    deldryvol_c_p: cobj,
    dryvol_a_xfab_p: cobj,
    dryvol_c_xfab_p: cobj,
):
    q = Ptr[float](q_p)
    qqcw = Ptr[float](qqcw_p)
    dqdt = Ptr[float](dqdt_p)
    dqdt_other = Ptr[float](dqdt_other_p)
    dqqcwdt = Ptr[float](dqqcwdt_p)
    dqqcwdt_other = Ptr[float](dqqcwdt_other_p)
    lspectype_amode = Ptr[int](lspectype_amode_p)
    specmw_amode = Ptr[float](specmw_amode_p)
    specdens_amode = Ptr[float](specdens_amode_p)
    lmassptr_amode = Ptr[int](lmassptr_amode_p)
    lmassptrcw_amode = Ptr[int](lmassptrcw_amode_p)
    ixferable_a_renamexf = Ptr[int](ixferable_a_renamexf_p)
    ixferable_c_renamexf = Ptr[int](ixferable_c_renamexf_p)
    dryvol_a = Ptr[float](dryvol_a_p)
    dryvol_c = Ptr[float](dryvol_c_p)
    deldryvol_a = Ptr[float](deldryvol_a_p)
    deldryvol_c = Ptr[float](deldryvol_c_p)
    dryvol_a_xfab = Ptr[float](dryvol_a_xfab_p)
    dryvol_c_xfab = Ptr[float](dryvol_c_xfab_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            dryvol_a[_idx2(i, k, ncol)] = 0.0
            dryvol_c[_idx2(i, k, ncol)] = 0.0
            deldryvol_a[_idx2(i, k, ncol)] = 0.0
            deldryvol_c[_idx2(i, k, ncol)] = 0.0
            dryvol_a_xfab[_idx2(i, k, ncol)] = 0.0
            dryvol_c_xfab[_idx2(i, k, ncol)] = 0.0

    for l1 in range(1, nspec_mfrm + 1):
        l2 = lspectype_amode[_idx2(l1, mfrm, maxspec_renamexf)]
        tmp_m2v = specmw_amode[l2 - 1] / specdens_amode[l2 - 1]
        tmp_m2vdt = tmp_m2v * deltat

        la = lmassptr_amode[_idx2(l1, mfrm, maxspec_renamexf)] - loffset
        if la > 0:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    dryvol_a[_idx2(i, k, ncol)] += tmp_m2v * max(
                        0.0,
                        q[_idx3(i, k, la, ncol, pver)]
                        - deltat * dqdt_other[_idx3(i, k, la, ncol, pver)],
                    )
                    deldryvol_a[_idx2(i, k, ncol)] += (
                        dqdt_other[_idx3(i, k, la, ncol, pver)]
                        + dqdt[_idx3(i, k, la, ncol, pver)]
                    ) * tmp_m2vdt
                    if (
                        ixferable_all <= 0
                        and ixferable_a_renamexf[_idx2(l1, ipair, maxspec_renamexf)] > 0
                    ):
                        dryvol_a_xfab[_idx2(i, k, ncol)] += tmp_m2v * max(
                            0.0,
                            q[_idx3(i, k, la, ncol, pver)]
                            + deltat * dqdt[_idx3(i, k, la, ncol, pver)],
                        )

        lc = lmassptrcw_amode[_idx2(l1, mfrm, maxspec_renamexf)] - loffset
        if lc > 0:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    dryvol_c[_idx2(i, k, ncol)] += tmp_m2v * max(
                        0.0,
                        qqcw[_idx3(i, k, lc, ncol, pver)]
                        - deltat * dqqcwdt_other[_idx3(i, k, lc, ncol, pver)],
                    )
                    deldryvol_c[_idx2(i, k, ncol)] += (
                        dqqcwdt_other[_idx3(i, k, lc, ncol, pver)]
                        + dqqcwdt[_idx3(i, k, lc, ncol, pver)]
                    ) * tmp_m2vdt
                    if (
                        ixferable_all <= 0
                        and ixferable_c_renamexf[_idx2(l1, ipair, maxspec_renamexf)] > 0
                    ):
                        dryvol_c_xfab[_idx2(i, k, ncol)] += tmp_m2v * max(
                            0.0,
                            qqcw[_idx3(i, k, lc, ncol, pver)]
                            + deltat * dqqcwdt[_idx3(i, k, lc, ncol, pver)],
                        )


def modal_aero_rename_acc_crs_tendencies_full_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnstxx: int,
    maxspec_renamexf: int,
    loffset: int,
    ipair: int,
    nspecfrm_ipair: int,
    jsrflx_rename: int,
    is_dorename_atik: int,
    l_dqdt_rnpos: int,
    deltat: float,
    deltatinv: float,
    gravit: float,
    pdel_p: cobj,
    dorename_atik_p: cobj,
    q_p: cobj,
    qqcw_p: cobj,
    dqdt_p: cobj,
    dqqcwdt_p: cobj,
    qsrflx_p: cobj,
    qqcwsrflx_p: cobj,
    xferfrac_vol_p: cobj,
    xferfrac_num_p: cobj,
    lspecfrma_renamexf_p: cobj,
    lspecfrmc_renamexf_p: cobj,
    lspectooa_renamexf_p: cobj,
    lspectooc_renamexf_p: cobj,
    dqdt_rnpos_p: cobj,
):
    pdel = Ptr[float](pdel_p)
    dorename_atik = Ptr[int](dorename_atik_p)
    q = Ptr[float](q_p)
    qqcw = Ptr[float](qqcw_p)
    dqdt = Ptr[float](dqdt_p)
    dqqcwdt = Ptr[float](dqqcwdt_p)
    qsrflx = Ptr[float](qsrflx_p)
    qqcwsrflx = Ptr[float](qqcwsrflx_p)
    xferfrac_vol = Ptr[float](xferfrac_vol_p)
    xferfrac_num = Ptr[float](xferfrac_num_p)
    lspecfrma_renamexf = Ptr[int](lspecfrma_renamexf_p)
    lspecfrmc_renamexf = Ptr[int](lspecfrmc_renamexf_p)
    lspectooa_renamexf = Ptr[int](lspectooa_renamexf_p)
    lspectooc_renamexf = Ptr[int](lspectooc_renamexf_p)
    dqdt_rnpos = Ptr[float](dqdt_rnpos_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            if is_dorename_atik != 0:
                if dorename_atik[_idx2(i, k, ncol)] == 0:
                    continue

            xferfrac_vol_local = xferfrac_vol[_idx2(i, k, ncol)]
            xferfrac_num_local = xferfrac_num[_idx2(i, k, ncol)]
            if xferfrac_vol_local <= 0.0:
                continue

            pdel_fac = pdel[_idx2(i, k, pcols)] / gravit

            for iq in range(1, nspecfrm_ipair + 1):
                xfercoef = xferfrac_vol_local * deltatinv
                if iq == 1:
                    xfercoef = xferfrac_num_local * deltatinv

                lsfrma = lspecfrma_renamexf[_idx2(iq, ipair, maxspec_renamexf)] - loffset
                lsfrmc = lspecfrmc_renamexf[_idx2(iq, ipair, maxspec_renamexf)] - loffset
                lstooa = lspectooa_renamexf[_idx2(iq, ipair, maxspec_renamexf)] - loffset
                lstooc = lspectooc_renamexf[_idx2(iq, ipair, maxspec_renamexf)] - loffset

                if lsfrma > 0:
                    xfertend = xfercoef * max(
                        0.0,
                        q[_idx3(i, k, lsfrma, ncol, pver)]
                        + dqdt[_idx3(i, k, lsfrma, ncol, pver)] * deltat,
                    )
                    dqdt[_idx3(i, k, lsfrma, ncol, pver)] -= xfertend
                    qsrflx[_idx3(i, lsfrma, jsrflx_rename, pcols, pcnstxx)] -= (
                        xfertend * pdel_fac
                    )
                    if lstooa > 0:
                        dqdt[_idx3(i, k, lstooa, ncol, pver)] += xfertend
                        qsrflx[_idx3(i, lstooa, jsrflx_rename, pcols, pcnstxx)] += (
                            xfertend * pdel_fac
                        )
                        if l_dqdt_rnpos != 0:
                            dqdt_rnpos[_idx3(i, k, lstooa, ncol, pver)] += xfertend

                if lsfrmc > 0:
                    xfertend = xfercoef * max(
                        0.0,
                        qqcw[_idx3(i, k, lsfrmc, ncol, pver)]
                        + dqqcwdt[_idx3(i, k, lsfrmc, ncol, pver)] * deltat,
                    )
                    dqqcwdt[_idx3(i, k, lsfrmc, ncol, pver)] -= xfertend
                    qqcwsrflx[_idx3(i, lsfrmc, jsrflx_rename, pcols, pcnstxx)] -= (
                        xfertend * pdel_fac
                    )
                    if lstooc > 0:
                        dqqcwdt[_idx3(i, k, lstooc, ncol, pver)] += xfertend
                        qqcwsrflx[_idx3(i, lstooc, jsrflx_rename, pcols, pcnstxx)] += (
                            xfertend * pdel_fac
                        )


@export
def modal_aero_rename_acc_crs_sub_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnstxx: int,
    maxpair_renamexf: int,
    maxspec_renamexf: int,
    loffset: int,
    npair_renamexf: int,
    is_dorename_atik: int,
    l_dqdt_rnpos: int,
    jsrflx_rename: int,
    nsrflx: int,
    modeptr_coarse: int,
    modeptr_accum: int,
    method_optbb: int,
    deltat: float,
    deltatinv: float,
    onethird: float,
    xferfrac_max: float,
    gravit: float,
    troplev_p: cobj,
    pdel_p: cobj,
    dorename_atik_p: cobj,
    q_p: cobj,
    qqcw_p: cobj,
    dqdt_p: cobj,
    dqdt_other_p: cobj,
    dqqcwdt_p: cobj,
    dqqcwdt_other_p: cobj,
    qsrflx_p: cobj,
    qqcwsrflx_p: cobj,
    modefrm_renamexf_p: cobj,
    modetoo_renamexf_p: cobj,
    nspec_amode_p: cobj,
    lspectype_amode_p: cobj,
    specmw_amode_p: cobj,
    specdens_amode_p: cobj,
    lmassptr_amode_p: cobj,
    lmassptrcw_amode_p: cobj,
    numptr_amode_p: cobj,
    numptrcw_amode_p: cobj,
    dgnum_amode_p: cobj,
    factoraa_p: cobj,
    factoryy_p: cobj,
    dryvol_smallest_p: cobj,
    v2nlorlx_p: cobj,
    v2nhirlx_p: cobj,
    factor_3alnsg2_p: cobj,
    dp_cut_p: cobj,
    lndp_cut_p: cobj,
    dp_belowcut_p: cobj,
    dp_xfernone_threshaa_p: cobj,
    dp_xferall_thresh_p: cobj,
    igrow_shrink_renamexf_p: cobj,
    ixferable_all_renamexf_p: cobj,
    ixferable_a_renamexf_p: cobj,
    ixferable_c_renamexf_p: cobj,
    nspecfrm_renamexf_p: cobj,
    lspecfrma_renamexf_p: cobj,
    lspecfrmc_renamexf_p: cobj,
    lspectooa_renamexf_p: cobj,
    lspectooc_renamexf_p: cobj,
    dryvol_a_p: cobj,
    dryvol_c_p: cobj,
    deldryvol_a_p: cobj,
    deldryvol_c_p: cobj,
    dryvol_a_xfab_p: cobj,
    dryvol_c_xfab_p: cobj,
    xferfrac_vol_p: cobj,
    xferfrac_num_p: cobj,
    dotendrn_p: cobj,
    dotendqqcwrn_p: cobj,
    dqdt_rnpos_p: cobj,
):
    modefrm_renamexf = Ptr[int](modefrm_renamexf_p)
    modetoo_renamexf = Ptr[int](modetoo_renamexf_p)
    nspec_amode = Ptr[int](nspec_amode_p)
    numptr_amode = Ptr[int](numptr_amode_p)
    numptrcw_amode = Ptr[int](numptrcw_amode_p)
    dgnum_amode = Ptr[float](dgnum_amode_p)
    factoraa = Ptr[float](factoraa_p)
    factoryy = Ptr[float](factoryy_p)
    dryvol_smallest = Ptr[float](dryvol_smallest_p)
    v2nlorlx = Ptr[float](v2nlorlx_p)
    v2nhirlx = Ptr[float](v2nhirlx_p)
    factor_3alnsg2 = Ptr[float](factor_3alnsg2_p)
    dp_cut = Ptr[float](dp_cut_p)
    lndp_cut = Ptr[float](lndp_cut_p)
    dp_belowcut = Ptr[float](dp_belowcut_p)
    dp_xfernone_threshaa = Ptr[float](dp_xfernone_threshaa_p)
    dp_xferall_thresh = Ptr[float](dp_xferall_thresh_p)
    igrow_shrink_renamexf = Ptr[int](igrow_shrink_renamexf_p)
    ixferable_all_renamexf = Ptr[int](ixferable_all_renamexf_p)
    nspecfrm_renamexf = Ptr[int](nspecfrm_renamexf_p)

    for ipair in range(1, npair_renamexf + 1):
        mfrm = modefrm_renamexf[ipair - 1]
        mtoo = modetoo_renamexf[ipair - 1]
        nspec_mfrm = nspec_amode[mfrm - 1]
        ixferable_all = ixferable_all_renamexf[ipair - 1]
        igrow_shrink = igrow_shrink_renamexf[ipair - 1]

        flagaa_shrink = 0
        if mfrm == modeptr_coarse and mtoo == modeptr_accum:
            flagaa_shrink = 1

        modal_aero_rename_acc_crs_dryvols_full_codon(
            ncol,
            pver,
            maxspec_renamexf,
            loffset,
            ipair,
            mfrm,
            nspec_mfrm,
            ixferable_all,
            deltat,
            q_p,
            qqcw_p,
            dqdt_p,
            dqdt_other_p,
            dqqcwdt_p,
            dqqcwdt_other_p,
            lspectype_amode_p,
            specmw_amode_p,
            specdens_amode_p,
            lmassptr_amode_p,
            lmassptrcw_amode_p,
            ixferable_a_renamexf_p,
            ixferable_c_renamexf_p,
            dryvol_a_p,
            dryvol_c_p,
            deldryvol_a_p,
            deldryvol_c_p,
            dryvol_a_xfab_p,
            dryvol_c_xfab_p,
        )

        modal_aero_rename_acc_crs_xferfracs_codon(
            ncol,
            pcols,
            pver,
            pcnstxx,
            loffset,
            mfrm,
            numptr_amode[mfrm - 1],
            numptrcw_amode[mfrm - 1],
            igrow_shrink,
            ixferable_all,
            method_optbb,
            flagaa_shrink,
            dgnum_amode[mfrm - 1],
            factoraa[mfrm - 1],
            factoryy[mfrm - 1],
            dryvol_smallest[mfrm - 1],
            v2nlorlx[mfrm - 1],
            v2nhirlx[mfrm - 1],
            factor_3alnsg2[ipair - 1],
            dp_cut[ipair - 1],
            lndp_cut[ipair - 1],
            dp_belowcut[ipair - 1],
            dp_xfernone_threshaa[ipair - 1],
            dp_xferall_thresh[ipair - 1],
            onethird,
            xferfrac_max,
            troplev_p,
            q_p,
            qqcw_p,
            dryvol_a_p,
            dryvol_c_p,
            deldryvol_a_p,
            deldryvol_c_p,
            dryvol_a_xfab_p,
            dryvol_c_xfab_p,
            xferfrac_vol_p,
            xferfrac_num_p,
        )

        modal_aero_rename_acc_crs_tendencies_full_codon(
            ncol,
            pcols,
            pver,
            pcnstxx,
            maxspec_renamexf,
            loffset,
            ipair,
            nspecfrm_renamexf[ipair - 1],
            jsrflx_rename,
            is_dorename_atik,
            l_dqdt_rnpos,
            deltat,
            deltatinv,
            gravit,
            pdel_p,
            dorename_atik_p,
            q_p,
            qqcw_p,
            dqdt_p,
            dqqcwdt_p,
            qsrflx_p,
            qqcwsrflx_p,
            xferfrac_vol_p,
            xferfrac_num_p,
            lspecfrma_renamexf_p,
            lspecfrmc_renamexf_p,
            lspectooa_renamexf_p,
            lspectooc_renamexf_p,
            dqdt_rnpos_p,
        )

    modal_aero_rename_set_dotend_flags_codon(
        pcnstxx,
        maxpair_renamexf,
        maxspec_renamexf,
        loffset,
        npair_renamexf,
        nspecfrm_renamexf_p,
        lspecfrma_renamexf_p,
        lspecfrmc_renamexf_p,
        lspectooa_renamexf_p,
        lspectooc_renamexf_p,
        dotendrn_p,
        dotendqqcwrn_p,
    )


@export
def modal_aero_rename_set_dotend_flags_codon(
    pcnstxx: int,
    maxpair_renamexf: int,
    maxspec_renamexf: int,
    loffset: int,
    npair_renamexf: int,
    nspecfrm_renamexf_p: cobj,
    lspecfrma_renamexf_p: cobj,
    lspecfrmc_renamexf_p: cobj,
    lspectooa_renamexf_p: cobj,
    lspectooc_renamexf_p: cobj,
    dotendrn_p: cobj,
    dotendqqcwrn_p: cobj,
):
    nspecfrm_renamexf = Ptr[int](nspecfrm_renamexf_p)
    lspecfrma_renamexf = Ptr[int](lspecfrma_renamexf_p)
    lspecfrmc_renamexf = Ptr[int](lspecfrmc_renamexf_p)
    lspectooa_renamexf = Ptr[int](lspectooa_renamexf_p)
    lspectooc_renamexf = Ptr[int](lspectooc_renamexf_p)
    dotendrn = Ptr[int](dotendrn_p)
    dotendqqcwrn = Ptr[int](dotendqqcwrn_p)

    for l in range(1, pcnstxx + 1):
        dotendrn[l - 1] = 0
        dotendqqcwrn[l - 1] = 0

    for ipair in range(1, npair_renamexf + 1):
        for iq in range(1, nspecfrm_renamexf[ipair - 1] + 1):
            lsfrma = lspecfrma_renamexf[_idx2(iq, ipair, maxspec_renamexf)] - loffset
            lsfrmc = lspecfrmc_renamexf[_idx2(iq, ipair, maxspec_renamexf)] - loffset
            lstooa = lspectooa_renamexf[_idx2(iq, ipair, maxspec_renamexf)] - loffset
            lstooc = lspectooc_renamexf[_idx2(iq, ipair, maxspec_renamexf)] - loffset

            if lsfrma > 0:
                dotendrn[lsfrma - 1] = 1
                if lstooa > 0:
                    dotendrn[lstooa - 1] = 1

            if lsfrmc > 0:
                dotendqqcwrn[lsfrmc - 1] = 1
                if lstooc > 0:
                    dotendqqcwrn[lstooc - 1] = 1


@export
def modal_aero_rename_acc_crs_xferfracs_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnstxx: int,
    loffset: int,
    mfrm: int,
    numptr_amode_mfrm: int,
    numptrcw_amode_mfrm: int,
    igrow_shrink: int,
    ixferable_all: int,
    method_optbb: int,
    flagaa_shrink: int,
    dgnum_amode_mfrm: float,
    factoraa: float,
    factoryy: float,
    dryvol_smallest: float,
    v2nlorlx: float,
    v2nhirlx: float,
    factor_3alnsg2: float,
    dp_cut: float,
    lndp_cut: float,
    dp_belowcut: float,
    dp_xfernone_thresh: float,
    dp_xferall_thresh: float,
    onethird: float,
    xferfrac_max: float,
    troplev_p: cobj,
    q_p: cobj,
    qqcw_p: cobj,
    dryvol_a_p: cobj,
    dryvol_c_p: cobj,
    deldryvol_a_p: cobj,
    deldryvol_c_p: cobj,
    dryvol_a_xfab_p: cobj,
    dryvol_c_xfab_p: cobj,
    xferfrac_vol_p: cobj,
    xferfrac_num_p: cobj,
):
    troplev = Ptr[int](troplev_p)
    q = Ptr[float](q_p)
    qqcw = Ptr[float](qqcw_p)
    dryvol_a = Ptr[float](dryvol_a_p)
    dryvol_c = Ptr[float](dryvol_c_p)
    deldryvol_a = Ptr[float](deldryvol_a_p)
    deldryvol_c = Ptr[float](deldryvol_c_p)
    dryvol_a_xfab = Ptr[float](dryvol_a_xfab_p)
    dryvol_c_xfab = Ptr[float](dryvol_c_xfab_p)
    xferfrac_vol = Ptr[float](xferfrac_vol_p)
    xferfrac_num = Ptr[float](xferfrac_num_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            xferfrac_vol[_idx2(i, k, ncol)] = 0.0
            xferfrac_num[_idx2(i, k, ncol)] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            dryvol_t_old = dryvol_a[_idx2(i, k, ncol)] + dryvol_c[_idx2(i, k, ncol)]
            dryvol_t_del = deldryvol_a[_idx2(i, k, ncol)] + deldryvol_c[_idx2(i, k, ncol)]
            dryvol_t_new = dryvol_t_old + dryvol_t_del
            dryvol_t_oldbnd = max(dryvol_t_old, dryvol_smallest)

            if igrow_shrink > 0:
                if dryvol_t_new <= dryvol_smallest:
                    continue
                if method_optbb != 2:
                    if dryvol_t_del <= 1.0e-6 * dryvol_t_oldbnd:
                        continue

                num_t_old = q[_idx3(i, k, numptr_amode_mfrm - loffset, ncol, pver)]
                num_t_old += qqcw[
                    _idx3(i, k, numptrcw_amode_mfrm - loffset, ncol, pver)
                ]
                num_t_old = max(0.0, num_t_old)
                dryvol_t_oldbnd = max(dryvol_t_old, dryvol_smallest)
                num_t_oldbnd = min(dryvol_t_oldbnd * v2nlorlx, num_t_old)
                num_t_oldbnd = max(dryvol_t_oldbnd * v2nhirlx, num_t_oldbnd)

                dgn_t_new = (dryvol_t_new / (num_t_oldbnd * factoraa)) ** onethird
                if dgn_t_new <= dp_xfernone_thresh:
                    continue

                dgn_t_old = (dryvol_t_oldbnd / (num_t_oldbnd * factoraa)) ** onethird
                dgn_t_oldb = dgn_t_old
                dryvol_t_oldb = dryvol_t_old
                if method_optbb == 2:
                    if dgn_t_old >= dp_cut:
                        dryvol_t_oldb = dryvol_t_old * (dp_belowcut / dgn_t_old) ** 3
                        dgn_t_oldb = dp_belowcut
                    if dgn_t_new < dp_xferall_thresh:
                        if (dryvol_t_new - dryvol_t_oldb) <= 1.0e-6 * dryvol_t_oldbnd:
                            continue
                elif dgn_t_new >= dp_cut:
                    dgn_t_oldb = min(dgn_t_oldb, dp_belowcut)

                lndgn_new = log(dgn_t_new)
                lndgv_new = lndgn_new + factor_3alnsg2
                yn_tail = (lndp_cut - lndgn_new) * factoryy
                yv_tail = (lndp_cut - lndgv_new) * factoryy
                tailfr_numnew = 0.5 * erfc(yn_tail)
                tailfr_volnew = 0.5 * erfc(yv_tail)

                lndgn_old = log(dgn_t_oldb)
                lndgv_old = lndgn_old + factor_3alnsg2
                yn_tail = (lndp_cut - lndgn_old) * factoryy
                yv_tail = (lndp_cut - lndgv_old) * factoryy
                tailfr_numold = 0.5 * erfc(yn_tail)
                tailfr_volold = 0.5 * erfc(yv_tail)

                if method_optbb == 2 and dgn_t_new >= dp_xferall_thresh:
                    dryvol_xferamt = dryvol_t_new
                else:
                    dryvol_xferamt = (
                        tailfr_volnew * dryvol_t_new - tailfr_volold * dryvol_t_oldb
                    )
                if dryvol_xferamt <= 0.0:
                    continue

                xferfrac_vol_val = max(0.0, dryvol_xferamt / dryvol_t_new)
                if method_optbb == 2 and xferfrac_vol_val >= xferfrac_max:
                    xferfrac_vol_val = 1.0
                    xferfrac_num_val = 1.0
                else:
                    xferfrac_vol_val = min(xferfrac_vol_val, xferfrac_max)
                    xferfrac_num_val = tailfr_numnew - tailfr_numold
                    xferfrac_num_val = max(
                        0.0, min(xferfrac_num_val, xferfrac_vol_val)
                    )

                if ixferable_all <= 0:
                    dryvol_t_new_xfab = max(
                        0.0,
                        dryvol_a_xfab[_idx2(i, k, ncol)] + dryvol_c_xfab[_idx2(i, k, ncol)],
                    )
                    dryvol_xferamt = xferfrac_vol_val * dryvol_t_new
                    if dryvol_t_new_xfab >= 0.999999 * dryvol_xferamt:
                        xferfrac_vol_val = min(1.0, dryvol_xferamt / dryvol_t_new_xfab)
                    elif dryvol_t_new_xfab >= 1.0e-7 * dryvol_xferamt:
                        xferfrac_vol_val = 1.0
                        xferfrac_num_val = xferfrac_num_val * (
                            dryvol_t_new_xfab / dryvol_xferamt
                        )
                    else:
                        continue

            else:
                if dryvol_t_old <= dryvol_smallest:
                    continue

                if dryvol_t_del >= -1.0e-6 * dryvol_t_oldbnd:
                    if flagaa_shrink != 0 and k < troplev[i - 1]:
                        flagbb_shrink = 1
                    else:
                        continue
                else:
                    flagbb_shrink = 0

                num_t_old = q[_idx3(i, k, numptr_amode_mfrm - loffset, ncol, pver)]
                num_t_old += qqcw[
                    _idx3(i, k, numptrcw_amode_mfrm - loffset, ncol, pver)
                ]
                num_t_old = max(0.0, num_t_old)
                dryvol_t_oldbnd = max(dryvol_t_old, dryvol_smallest)
                num_t_oldbnd = min(dryvol_t_oldbnd * v2nlorlx, num_t_old)
                num_t_oldbnd = max(dryvol_t_oldbnd * v2nhirlx, num_t_oldbnd)

                dgn_t_new = (dryvol_t_new / (num_t_oldbnd * factoraa)) ** onethird
                if dgn_t_new >= dp_xfernone_thresh:
                    continue
                if flagbb_shrink != 0:
                    if dgn_t_new > dp_cut:
                        continue

                if dgn_t_new <= dp_xferall_thresh:
                    tailfr_numnew = 1.0
                    tailfr_volnew = 1.0
                else:
                    lndgn_new = log(dgn_t_new)
                    lndgv_new = lndgn_new + factor_3alnsg2
                    yn_tail = (lndp_cut - lndgn_new) * factoryy
                    yv_tail = (lndp_cut - lndgv_new) * factoryy
                    tailfr_numnew = 1.0 - 0.5 * erfc(yn_tail)
                    tailfr_volnew = 1.0 - 0.5 * erfc(yv_tail)

                dgn_t_old = (dryvol_t_oldbnd / (num_t_oldbnd * factoraa)) ** onethird
                dgn_t_oldb = dgn_t_old
                dryvol_t_oldb = dryvol_t_old
                tailfr_numold = 0.0
                tailfr_volold = 0.0

                xferfrac_vol_val = tailfr_volnew
                if xferfrac_vol_val <= 0.0:
                    continue
                xferfrac_num_val = tailfr_numnew

                if xferfrac_vol_val >= xferfrac_max:
                    xferfrac_vol_val = 1.0
                    xferfrac_num_val = 1.0
                else:
                    xferfrac_vol_val = min(xferfrac_vol_val, xferfrac_max)
                    xferfrac_num_val = max(xferfrac_num_val, xferfrac_vol_val)
                    xferfrac_num_val = min(xferfrac_max, xferfrac_num_val)

                if ixferable_all <= 0:
                    dryvol_t_new_xfab = max(
                        0.0,
                        dryvol_a_xfab[_idx2(i, k, ncol)] + dryvol_c_xfab[_idx2(i, k, ncol)],
                    )
                    dryvol_xferamt = xferfrac_vol_val * dryvol_t_new
                    if dryvol_t_new_xfab >= 0.999999 * dryvol_xferamt:
                        xferfrac_vol_val = min(1.0, dryvol_xferamt / dryvol_t_new_xfab)
                    elif dryvol_t_new_xfab >= 1.0e-7 * dryvol_xferamt:
                        xferfrac_vol_val = 1.0
                        xferfrac_num_val = xferfrac_num_val * (
                            dryvol_t_new_xfab / dryvol_xferamt
                        )
                    else:
                        continue

            xferfrac_vol[_idx2(i, k, ncol)] = xferfrac_vol_val
            xferfrac_num[_idx2(i, k, ncol)] = xferfrac_num_val
