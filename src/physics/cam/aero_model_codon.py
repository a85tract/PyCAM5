from math import exp, log, sqrt

@inline
def _idx2(i: int, k: int, ld1: int) -> int:
    return (i - 1) + (k - 1) * ld1


@inline
def _idx3(i: int, k: int, m: int, ld1: int, ld2: int) -> int:
    return (i - 1) + (k - 1) * ld1 + (m - 1) * ld1 * ld2


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
