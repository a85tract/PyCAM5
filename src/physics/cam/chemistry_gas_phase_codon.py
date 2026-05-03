from math import acos, cos, exp, sqrt
from chemistry_common_codon import _idx2, _idx3, _idx3_k0, _flux_idx

def gas_phase_chemdr_prepare_sza_codon(
    ncol: int,
    rad2deg: float,
    zen_angle_p: cobj,
    sza_p: cobj,
):
    zen_angle = Ptr[float](zen_angle_p)
    sza = Ptr[float](sza_p)

    for i in range(1, ncol + 1):
        z = acos(zen_angle[i - 1])
        zen_angle[i - 1] = z
        sza[i - 1] = z * rad2deg

def mmr2vmr_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    mbar_p: cobj,
    mmr_p: cobj,
    adv_mass_p: cobj,
    vmr_p: cobj,
):
    mbar = Ptr[float](mbar_p)
    mmr = Ptr[float](mmr_p)
    adv_mass = Ptr[float](adv_mass_p)
    vmr = Ptr[float](vmr_p)

    for m in range(1, gas_pcnst + 1):
        adv = adv_mass[m - 1]
        if adv != 0.0:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    vmr[_idx3(i, k, m, ncol, pver)] = mbar[_idx2(i, k, ncol)] * mmr[
                        _idx3(i, k, m, pcols, pver)
                    ] / adv

def vmr2mmr_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    mbar_p: cobj,
    vmr_p: cobj,
    adv_mass_p: cobj,
    mmr_p: cobj,
):
    mbar = Ptr[float](mbar_p)
    vmr = Ptr[float](vmr_p)
    adv_mass = Ptr[float](adv_mass_p)
    mmr = Ptr[float](mmr_p)

    for m in range(1, gas_pcnst + 1):
        adv = adv_mass[m - 1]
        if adv != 0.0:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    mmr[_idx3(i, k, m, pcols, pver)] = adv * vmr[_idx3(i, k, m, ncol, pver)] / mbar[
                        _idx2(i, k, ncol)
                    ]

def h2o_to_vmr_codon(
    ncol: int,
    pcols: int,
    pver: int,
    adv_mass_h2o: float,
    h2o_mmr_p: cobj,
    mbar_p: cobj,
    h2o_vmr_p: cobj,
):
    h2o_mmr = Ptr[float](h2o_mmr_p)
    mbar = Ptr[float](mbar_p)
    h2o_vmr = Ptr[float](h2o_vmr_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            h2o_vmr[_idx2(i, k, ncol)] = mbar[_idx2(i, k, ncol)] * h2o_mmr[_idx2(i, k, pcols)] / adv_mass_h2o

@inline
def _gas_phase_chemdr_shell_h2o_from_q(
    ncol: int,
    pcols: int,
    pver: int,
    adv_mass_h2o: float,
    q_p: cobj,
    mbar_p: cobj,
    qh2o_p: cobj,
    h2ovmr_p: cobj,
):
    q = Ptr[float](q_p)
    qh2o = Ptr[float](qh2o_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            qh2o[_idx2(i, k, pcols)] = q[_idx3(i, k, 1, pcols, pver)]

    h2o_to_vmr_codon(ncol, pcols, pver, adv_mass_h2o, q_p, mbar_p, h2ovmr_p)

def set_mean_mass_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    id_o2: int,
    id_o: int,
    id_h: int,
    id_n: int,
    fixed_mbar: int,
    mwdry: float,
    mmr_p: cobj,
    adv_mass_p: cobj,
    mbar_p: cobj,
):
    mmr = Ptr[float](mmr_p)
    adv_mass = Ptr[float](adv_mass_p)
    mbar = Ptr[float](mbar_p)

    if fixed_mbar != 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                mbar[_idx2(i, k, ncol)] = mwdry
        return

    adv_n = adv_mass[id_n - 1]
    adv_o2 = adv_mass[id_o2 - 1]
    adv_o = adv_mass[id_o - 1]
    adv_h = adv_mass[id_h - 1]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            xn2 = 1.0 - (
                mmr[_idx3(i, k, id_o2, pcols, pver)]
                + mmr[_idx3(i, k, id_o, pcols, pver)]
                + mmr[_idx3(i, k, id_h, pcols, pver)]
            )
            fn2 = 0.5 * xn2 / adv_n
            fo2 = mmr[_idx3(i, k, id_o2, pcols, pver)] / adv_o2
            fo = mmr[_idx3(i, k, id_o, pcols, pver)] / adv_o
            fh = mmr[_idx3(i, k, id_h, pcols, pver)] / adv_h
            mbar[_idx2(i, k, ncol)] = 1.0 / (fn2 + fo2 + fo + fh)

@inline
def _gas_phase_chemdr_shell_mass_vmr(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    o2_ndx: int,
    o_ndx: int,
    h_ndx: int,
    n_ndx: int,
    fixed_mbar: int,
    mwdry: float,
    mmr_p: cobj,
    adv_mass_p: cobj,
    mbar_p: cobj,
    vmr_p: cobj,
):
    set_mean_mass_codon(
        ncol, pcols, pver, gas_pcnst, o2_ndx, o_ndx, h_ndx, n_ndx, fixed_mbar, mwdry, mmr_p,
        adv_mass_p, mbar_p
    )
    mmr2vmr_codon(ncol, pcols, pver, gas_pcnst, mbar_p, mmr_p, adv_mass_p, vmr_p)

def setinv_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    nfs: int,
    m_ndx: int,
    n2_ndx: int,
    o2_ndx: int,
    h2o_ndx: int,
    id_o: int,
    id_o2: int,
    id_h: int,
    has_n2: int,
    has_o2: int,
    has_h2o: int,
    has_var_o2: int,
    pa_xfac: float,
    boltz_cgs: float,
    tfld_p: cobj,
    h2ovmr_p: cobj,
    vmr_p: cobj,
    pmid_p: cobj,
    invariants_p: cobj,
):
    tfld = Ptr[float](tfld_p)
    h2ovmr = Ptr[float](h2ovmr_p)
    vmr = Ptr[float](vmr_p)
    pmid = Ptr[float](pmid_p)
    invariants = Ptr[float](invariants_p)

    for m in range(1, nfs + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                invariants[_idx3(i, k, m, ncol, pver)] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            invariants[_idx3(i, k, m_ndx, ncol, pver)] = (
                pa_xfac * pmid[_idx2(i, k, pcols)] / (boltz_cgs * tfld[_idx2(i, k, pcols)])
            )

    if has_n2 != 0:
        if has_var_o2 != 0:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    sum1 = (
                        vmr[_idx3(i, k, id_o, ncol, pver)]
                        + vmr[_idx3(i, k, id_o2, ncol, pver)]
                        + vmr[_idx3(i, k, id_h, ncol, pver)]
                    )
                    invariants[_idx3(i, k, n2_ndx, ncol, pver)] = (
                        (1.0 - sum1) * invariants[_idx3(i, k, m_ndx, ncol, pver)]
                    )
        else:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    invariants[_idx3(i, k, n2_ndx, ncol, pver)] = (
                        0.79 * invariants[_idx3(i, k, m_ndx, ncol, pver)]
                    )

    if has_o2 != 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                invariants[_idx3(i, k, o2_ndx, ncol, pver)] = (
                    0.21 * invariants[_idx3(i, k, m_ndx, ncol, pver)]
                )

    if has_h2o != 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                invariants[_idx3(i, k, h2o_ndx, ncol, pver)] = (
                    h2ovmr[_idx2(i, k, ncol)] * invariants[_idx3(i, k, m_ndx, ncol, pver)]
                )

def charge_balance_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    np_ndx: int,
    n2p_ndx: int,
    op_ndx: int,
    o2p_ndx: int,
    nop_ndx: int,
    conc_p: cobj,
    wrk_p: cobj,
):
    conc = Ptr[float](conc_p)
    wrk = Ptr[float](wrk_p)

    if np_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                wrk[_idx2(i, k, ncol)] = wrk[_idx2(i, k, ncol)] + conc[_idx3(i, k, np_ndx, ncol, pver)]

    if n2p_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                wrk[_idx2(i, k, ncol)] = wrk[_idx2(i, k, ncol)] + conc[_idx3(i, k, n2p_ndx, ncol, pver)]

    if op_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                wrk[_idx2(i, k, ncol)] = wrk[_idx2(i, k, ncol)] + conc[_idx3(i, k, op_ndx, ncol, pver)]

    if o2p_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                wrk[_idx2(i, k, ncol)] = wrk[_idx2(i, k, ncol)] + conc[_idx3(i, k, o2p_ndx, ncol, pver)]

    if nop_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                wrk[_idx2(i, k, ncol)] = wrk[_idx2(i, k, ncol)] + conc[_idx3(i, k, nop_ndx, ncol, pver)]

@inline
def _gas_phase_chemdr_shell_charge_balance(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    elec_ndx: int,
    np_ndx: int,
    n2p_ndx: int,
    op_ndx: int,
    o2p_ndx: int,
    nop_ndx: int,
    vmr_p: cobj,
    wrk_p: cobj,
):
    if elec_ndx > 0:
        wrk = Ptr[float](wrk_p)
        vmr = Ptr[float](vmr_p)
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                wrk[_idx2(i, k, ncol)] = 0.0

        charge_balance_codon(ncol, pver, gas_pcnst, np_ndx, n2p_ndx, op_ndx, o2p_ndx, nop_ndx, vmr_p, wrk_p)

        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                vmr[_idx3(i, k, elec_ndx, ncol, pver)] = wrk[_idx2(i, k, ncol)]

def setcol_codon(
    ncol: int,
    pver: int,
    ncol_abs: int,
    col_delta_p: cobj,
    col_dens_p: cobj,
):
    col_delta = Ptr[float](col_delta_p)
    col_dens = Ptr[float](col_dens_p)

    for m in range(1, ncol_abs + 1):
        for i in range(1, ncol + 1):
            col_dens[_idx3(i, 1, m, ncol, pver)] = col_delta[_idx3_k0(i, 0, m, ncol, pver + 1)] + 0.5 * col_delta[
                _idx3_k0(i, 1, m, ncol, pver + 1)
            ]

        for k in range(2, pver + 1):
            km1 = k - 1
            for i in range(1, ncol + 1):
                col_dens[_idx3(i, k, m, ncol, pver)] = col_dens[_idx3(i, km1, m, ncol, pver)] + 0.5 * (
                    col_delta[_idx3_k0(i, km1, m, ncol, pver + 1)] + col_delta[_idx3_k0(i, k, m, ncol, pver + 1)]
                )

@inline
def _gas_phase_chemdr_shell_adjrxt_setcol(
    ncol: int,
    pver: int,
    ncol_abs: int,
    indexm: int,
    reaction_rates_p: cobj,
    invariants_p: cobj,
    col_delta_p: cobj,
    col_dens_p: cobj,
):
    rate = Ptr[float](reaction_rates_p)
    inv = Ptr[float](invariants_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx_r2 = _idx3(i, k, 2, ncol, pver)
            idx_r3 = _idx3(i, k, 3, ncol, pver)
            idx_r4 = _idx3(i, k, 4, ncol, pver)
            idx_r5 = _idx3(i, k, 5, ncol, pver)
            idx_r6 = _idx3(i, k, 6, ncol, pver)
            idx_r7 = _idx3(i, k, 7, ncol, pver)
            inv6 = inv[_idx3(i, k, 6, ncol, pver)]
            inv7 = inv[_idx3(i, k, 7, ncol, pver)]
            inv8 = inv[_idx3(i, k, 8, ncol, pver)]
            im = 1.0 / inv[_idx3(i, k, indexm, ncol, pver)]

            rate[idx_r3] = rate[idx_r3] * inv6
            rate[idx_r4] = rate[idx_r4] * inv6
            rate[idx_r5] = rate[idx_r5] * inv6
            rate[idx_r6] = rate[idx_r6] * inv6
            rate[idx_r7] = rate[idx_r7] * inv7
            rate[idx_r2] = rate[idx_r2] * inv8 * inv8 * im

    setcol_codon(ncol, pver, ncol_abs, col_delta_p, col_dens_p)

def set_ub_col_codon(
    ncol: int,
    pcols: int,
    pver: int,
    ncol_abs: int,
    indexm: int,
    o3rad_ndx: int,
    ox_ndx: int,
    o3_ndx: int,
    o3_inv_ndx: int,
    o2_ndx: int,
    o2_is_inv: int,
    xfactor: float,
    pdel_p: cobj,
    vmr_p: cobj,
    invariants_p: cobj,
    o2_exo_col_p: cobj,
    o3_exo_col_p: cobj,
    col_delta_p: cobj,
):
    pdel = Ptr[float](pdel_p)
    vmr = Ptr[float](vmr_p)
    invariants = Ptr[float](invariants_p)
    o2_exo_col = Ptr[float](o2_exo_col_p)
    o3_exo_col = Ptr[float](o3_exo_col_p)
    col_delta = Ptr[float](col_delta_p)

    spc_ndx = o3rad_ndx
    if spc_ndx <= 0:
        spc_ndx = ox_ndx
    if spc_ndx < 1:
        spc_ndx = o3_ndx

    if spc_ndx > 0:
        for i in range(1, ncol + 1):
            col_delta[_idx3_k0(i, 0, 1, ncol, pver + 1)] = o3_exo_col[i - 1]
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                col_delta[_idx3_k0(i, k, 1, ncol, pver + 1)] = (
                    xfactor * pdel[_idx2(i, k, pcols)] * vmr[_idx3(i, k, spc_ndx, ncol, pver)]
                )
    elif o3_inv_ndx > 0:
        for i in range(1, ncol + 1):
            col_delta[_idx3_k0(i, 0, 1, ncol, pver + 1)] = o3_exo_col[i - 1]
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                col_delta[_idx3_k0(i, k, 1, ncol, pver + 1)] = (
                    xfactor
                    * pdel[_idx2(i, k, pcols)]
                    * invariants[_idx3(i, k, o3_inv_ndx, ncol, pver)]
                    / invariants[_idx3(i, k, indexm, ncol, pver)]
                )
    else:
        for k in range(0, pver + 1):
            for i in range(1, ncol + 1):
                col_delta[_idx3_k0(i, k, 1, ncol, pver + 1)] = 0.0

    if ncol_abs > 1:
        if o2_ndx > 1:
            for i in range(1, ncol + 1):
                col_delta[_idx3_k0(i, 0, 2, ncol, pver + 1)] = o2_exo_col[i - 1]
            if o2_is_inv != 0:
                for k in range(1, pver + 1):
                    for i in range(1, ncol + 1):
                        col_delta[_idx3_k0(i, k, 2, ncol, pver + 1)] = (
                            xfactor
                            * pdel[_idx2(i, k, pcols)]
                            * invariants[_idx3(i, k, o2_ndx, ncol, pver)]
                            / invariants[_idx3(i, k, indexm, ncol, pver)]
                        )
            else:
                for k in range(1, pver + 1):
                    for i in range(1, ncol + 1):
                        col_delta[_idx3_k0(i, k, 2, ncol, pver + 1)] = (
                            xfactor * pdel[_idx2(i, k, pcols)] * vmr[_idx3(i, k, o2_ndx, ncol, pver)]
                        )
        else:
            for k in range(0, pver + 1):
                for i in range(1, ncol + 1):
                    col_delta[_idx3_k0(i, k, 2, ncol, pver + 1)] = 0.0

def cloud_mod_codon(
    pver: int,
    zen_angle: float,
    srf_alb: float,
    rgrav: float,
    clouds_p: cobj,
    lwc_p: cobj,
    delp_p: cobj,
    eff_alb_p: cobj,
    cld_mult_p: cobj,
    del_lwp_p: cobj,
    del_tau_p: cobj,
    above_tau_p: cobj,
    below_tau_p: cobj,
    above_cld_p: cobj,
    below_cld_p: cobj,
    above_tra_p: cobj,
    below_tra_p: cobj,
    fac1_p: cobj,
    fac2_p: cobj,
):
    clouds = Ptr[float](clouds_p)
    lwc = Ptr[float](lwc_p)
    delp = Ptr[float](delp_p)
    eff_alb = Ptr[float](eff_alb_p)
    cld_mult = Ptr[float](cld_mult_p)
    del_lwp = Ptr[float](del_lwp_p)
    del_tau = Ptr[float](del_tau_p)
    above_tau = Ptr[float](above_tau_p)
    below_tau = Ptr[float](below_tau_p)
    above_cld = Ptr[float](above_cld_p)
    below_cld = Ptr[float](below_cld_p)
    above_tra = Ptr[float](above_tra_p)
    below_tra = Ptr[float](below_tra_p)
    fac1 = Ptr[float](fac1_p)
    fac2 = Ptr[float](fac2_p)

    for k in range(1, pver + 1):
        if clouds[k - 1] != 0.0:
            del_lwp[k - 1] = rgrav * lwc[k - 1] * delp[k - 1] * 1.0e3 / clouds[k - 1]
        else:
            del_lwp[k - 1] = 0.0

    for k in range(1, pver + 1):
        if clouds[k - 1] != 0.0:
            del_tau[k - 1] = del_lwp[k - 1] * 0.155 * (clouds[k - 1] ** 1.5)
        else:
            del_tau[k - 1] = 0.0

    above_tau[0] = 0.0
    for k in range(1, pver):
        above_tau[k] = del_tau[k - 1] + above_tau[k - 1]

    below_tau[pver - 1] = 0.0
    for k in range(pver - 1, 0, -1):
        below_tau[k - 1] = del_tau[k] + below_tau[k]

    above_cld[0] = 0.0
    for k in range(1, pver):
        above_cld[k] = clouds[k - 1] * del_tau[k - 1] + above_cld[k - 1]
    for k in range(1, pver):
        if above_tau[k] != 0.0:
            above_cld[k] = above_cld[k] / above_tau[k]
        else:
            above_cld[k] = above_cld[k - 1]

    below_cld[pver - 1] = 0.0
    for k in range(pver - 1, 0, -1):
        below_cld[k - 1] = clouds[k] * del_tau[k] + below_cld[k]
    for k in range(pver - 2, -1, -1):
        if below_tau[k] != 0.0:
            below_cld[k] = below_cld[k] / below_tau[k]
        else:
            below_cld[k] = below_cld[k + 1]

    for k in range(1, pver):
        if above_cld[k] != 0.0:
            above_tau[k] = above_tau[k] / above_cld[k]
    for k in range(0, pver - 1):
        if below_cld[k] != 0.0:
            below_tau[k] = below_tau[k] / below_cld[k]

    for k in range(1, pver):
        if above_tau[k] < 5.0:
            above_cld[k] = 0.0
    for k in range(0, pver - 1):
        if below_tau[k] < 5.0:
            below_cld[k] = 0.0

    for k in range(1, pver + 1):
        above_tra[k - 1] = 11.905 / (9.524 + above_tau[k - 1])
        below_tra[k - 1] = 11.905 / (9.524 + below_tau[k - 1])

    for k in range(1, pver + 1):
        if below_cld[k - 1] != 0.0:
            eff_alb[k - 1] = srf_alb + below_cld[k - 1] * (1.0 - below_tra[k - 1]) * (1.0 - srf_alb)
        else:
            eff_alb[k - 1] = srf_alb

    coschi = cos(zen_angle)
    if coschi < 0.5:
        coschi = 0.5

    for k in range(1, pver + 1):
        if del_lwp[k - 1] * 0.155 < 5.0:
            fac1[k - 1] = 0.0
        else:
            fac1[k - 1] = 1.4 * coschi - 1.0

    for k in range(1, pver + 1):
        fac2_val = 1.6 * coschi * above_tra[k - 1] - 1.0
        if fac2_val > 0.0:
            fac2[k - 1] = 0.0
        else:
            fac2[k - 1] = fac2_val

    for k in range(1, pver + 1):
        cld_mult[k - 1] = 1.0 + fac1[k - 1] * clouds[k - 1] + fac2[k - 1] * above_cld[k - 1]
        if cld_mult[k - 1] < 0.05:
            cld_mult[k - 1] = 0.05

def gas_phase_chemdr_zero_sulfate_codon(
    ncol: int,
    pver: int,
    sulfate_p: cobj,
):
    sulfate = Ptr[float](sulfate_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            sulfate[_idx2(i, k, ncol)] = 0.0

def gas_phase_chemdr_load_prognostic_sulfate_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    so4_ndx: int,
    vmr_p: cobj,
    sulfate_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    sulfate = Ptr[float](sulfate_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            sulfate[_idx2(i, k, ncol)] = vmr[_idx3(i, k, so4_ndx, ncol, pver)]

def chem_timestep_init_should_run_codon(
    nstep: int,
    chem_freq: int,
    chem_step_flag_p: cobj,
):
    chem_step_flag = Ptr[int](chem_step_flag_p)
    chem_step_flag[0] = 1 if nstep % chem_freq == 0 else 0

def chem_timestep_tend_fill_cloud_fields_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    ixcldliq: int,
    ixcldice: int,
    ixndrop: int,
    state_q_p: cobj,
    cldw_p: cobj,
    ncldwtr_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    cldw = Ptr[float](cldw_p)
    ncldwtr = Ptr[float](ncldwtr_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            cldw[_idx2(i, k, pcols)] = state_q[_idx3(i, k, ixcldliq, pcols, pver)] + state_q[
                _idx3(i, k, ixcldice, pcols, pver)
            ]
            if ixndrop > 0:
                ncldwtr[_idx2(i, k, pcols)] = state_q[
                    _idx3(i, k, ixndrop, pcols, pver)
                ]

def chem_timestep_tend_init_lq_codon(
    pcnst: int,
    ghg_chem: int,
    map2chm_p: cobj,
    lq_mask_p: cobj,
):
    map2chm = Ptr[int](map2chm_p)
    lq_mask = Ptr[int](lq_mask_p)

    for n in range(1, pcnst + 1):
        lq_mask[n - 1] = 1 if map2chm[n - 1] > 0 else 0

    if ghg_chem != 0 and pcnst > 0:
        lq_mask[0] = 1

def chem_timestep_tend_apply_depflux_codon(
    ncol: int,
    pcols: int,
    idx_cb1: int,
    idx_cb2: int,
    idx_oc1: int,
    idx_oc2: int,
    drydepflx_p: cobj,
    bcphodry_p: cobj,
    bcphidry_p: cobj,
    ocphodry_p: cobj,
    ocphidry_p: cobj,
):
    drydepflx = Ptr[float](drydepflx_p)
    bcphodry = Ptr[float](bcphodry_p)
    bcphidry = Ptr[float](bcphidry_p)
    ocphodry = Ptr[float](ocphodry_p)
    ocphidry = Ptr[float](ocphidry_p)

    if idx_cb1 > 0:
        for i in range(1, ncol + 1):
            bcphodry[i - 1] = max(drydepflx[_flux_idx(i, idx_cb1, pcols)], 0.0)

    if idx_cb2 > 0:
        for i in range(1, ncol + 1):
            bcphidry[i - 1] = max(drydepflx[_flux_idx(i, idx_cb2, pcols)], 0.0)

    if idx_oc1 > 0:
        for i in range(1, ncol + 1):
            ocphodry[i - 1] = max(drydepflx[_flux_idx(i, idx_oc1, pcols)], 0.0)

    if idx_oc2 > 0:
        for i in range(1, ncol + 1):
            ocphidry[i - 1] = max(drydepflx[_flux_idx(i, idx_oc2, pcols)], 0.0)

def chem_timestep_tend_sum_fh2o_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gravit: float,
    ptend_q1_p: cobj,
    pdel_p: cobj,
    fh2o_p: cobj,
):
    ptend_q1 = Ptr[float](ptend_q1_p)
    pdel = Ptr[float](pdel_p)
    fh2o = Ptr[float](fh2o_p)

    for i in range(1, ncol + 1):
        total = 0.0
        for k in range(1, pver + 1):
            total += ptend_q1[_idx2(i, k, pcols)] * pdel[_idx2(i, k, pcols)] / gravit
        fh2o[i - 1] = total

def gas_phase_chemdr_finalize_tendencies_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    pcnst: int,
    delt_inverse: float,
    map2chm_p: cobj,
    mmr_p: cobj,
    mmr_tend_p: cobj,
    mmr_new_p: cobj,
    qtend_p: cobj,
):
    map2chm = Ptr[int](map2chm_p)
    mmr = Ptr[float](mmr_p)
    mmr_tend = Ptr[float](mmr_tend_p)
    mmr_new = Ptr[float](mmr_new_p)
    qtend = Ptr[float](qtend_p)

    for m in range(1, gas_pcnst + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx3(i, k, m, pcols, pver)
                mmr_new[idx] = mmr_tend[idx]
                mmr_tend[idx] = (mmr_tend[idx] - mmr[idx]) * delt_inverse

    for m in range(1, pcnst + 1):
        n = map2chm[m - 1]
        if n > 0:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    qtend[_idx3(i, k, m, pcols, pver)] += mmr_tend[_idx3(i, k, n, pcols, pver)]

def gas_phase_chemdr_prepare_state_codon(
    ncol: int,
    pcols: int,
    pver: int,
    rga: float,
    m2km: float,
    pa2mb: float,
    phis_p: cobj,
    zi_p: cobj,
    zm_p: cobj,
    pmid_p: cobj,
    zsurf_p: cobj,
    zintr_p: cobj,
    zmidr_p: cobj,
    zmid_p: cobj,
    zint_p: cobj,
    pmb_p: cobj,
):
    phis = Ptr[float](phis_p)
    zi = Ptr[float](zi_p)
    zm = Ptr[float](zm_p)
    pmid = Ptr[float](pmid_p)
    zsurf = Ptr[float](zsurf_p)
    zintr = Ptr[float](zintr_p)
    zmidr = Ptr[float](zmidr_p)
    zmid = Ptr[float](zmid_p)
    zint = Ptr[float](zint_p)
    pmb = Ptr[float](pmb_p)

    for i in range(1, ncol + 1):
        zsurf[i - 1] = rga * phis[i - 1]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            zi_in_idx = _idx2(i, k, pcols)
            zm_in_idx = _idx2(i, k, pcols)
            out_idx = _idx2(i, k, ncol)
            zsurf_val = zsurf[i - 1]
            zintr[out_idx] = m2km * zi[zi_in_idx]
            zmidr[out_idx] = m2km * zm[zm_in_idx]
            zmid[out_idx] = m2km * (zm[zm_in_idx] + zsurf_val)
            zint[out_idx] = m2km * (zi[zi_in_idx] + zsurf_val)
            pmb[out_idx] = pa2mb * pmid[zm_in_idx]

    for i in range(1, ncol + 1):
        zi_in_idx = _idx2(i, pver + 1, pcols)
        zi_out_idx = _idx2(i, pver + 1, ncol)
        zint[zi_out_idx] = m2km * (zi[zi_in_idx] + zsurf[i - 1])
        zintr[zi_out_idx] = m2km * zi[zi_in_idx]

def gas_phase_chemdr_load_mmr_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    map2chm_p: cobj,
    q_p: cobj,
    mmr_p: cobj,
):
    map2chm = Ptr[int](map2chm_p)
    q = Ptr[float](q_p)
    mmr = Ptr[float](mmr_p)

    for m in range(1, pcnst + 1):
        n = map2chm[m - 1]
        if n > 0:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    idx_q = _idx3(i, k, m, pcols, pver)
                    idx_mmr = _idx3(i, k, n, pcols, pver)
                    mmr[idx_mmr] = q[idx_q]

def gas_phase_chemdr_init_reaction_rates_codon(
    ncol: int,
    pver: int,
    rxntot: int,
    nan_value_p: cobj,
    reaction_rates_p: cobj,
):
    nan_value = Ptr[float](nan_value_p)[0]
    reaction_rates = Ptr[float](reaction_rates_p)

    for m in range(1, max(1, rxntot) + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                reaction_rates[_idx3(i, k, m, ncol, pver)] = nan_value

def gas_phase_chemdr_clip_sulfate_codon(
    ncol: int,
    pcols: int,
    pver: int,
    troplev_p: cobj,
    sulfate_p: cobj,
):
    troplev = Ptr[int](troplev_p)
    sulfate = Ptr[float](sulfate_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            if k < troplev[i - 1]:
                sulfate[_idx2(i, k, ncol)] = 0.0

def gas_phase_chemdr_zero_het_rates_codon(
    ncol: int,
    pver: int,
    gas_pcnst_dim: int,
    het_rates_p: cobj,
):
    het_rates = Ptr[float](het_rates_p)

    for m in range(1, gas_pcnst_dim + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                het_rates[_idx3(i, k, m, ncol, pver)] = 0.0

def gas_phase_chemdr_load_oxygen_mmr_codon(
    ncol: int,
    pcols: int,
    pver: int,
    o2_ndx: int,
    o_ndx: int,
    mmr_p: cobj,
    o2mmr_p: cobj,
    ommr_p: cobj,
):
    mmr = Ptr[float](mmr_p)
    o2mmr = Ptr[float](o2mmr_p)
    ommr = Ptr[float](ommr_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            o2mmr[_idx2(i, k, ncol)] = mmr[_idx3(i, k, o2_ndx, pcols, pver)]
            ommr[_idx2(i, k, ncol)] = mmr[_idx3(i, k, o_ndx, pcols, pver)]

def gas_phase_chemdr_set_ltrop_sol_codon(
    ncol: int,
    has_linoz_data_flag: int,
    troplev_p: cobj,
    ltrop_sol_p: cobj,
):
    troplev = Ptr[int](troplev_p)
    ltrop_sol = Ptr[int](ltrop_sol_p)

    if has_linoz_data_flag != 0:
        for i in range(1, ncol + 1):
            ltrop_sol[i - 1] = troplev[i - 1]
    else:
        for i in range(1, ncol + 1):
            ltrop_sol[i - 1] = 0

def gas_phase_chemdr_zero_st80_tau_codon(
    ncol: int,
    pver: int,
    rxntot: int,
    st80_25_tau_ndx: int,
    troplev_p: cobj,
    reaction_rates_p: cobj,
):
    troplev = Ptr[int](troplev_p)
    reaction_rates = Ptr[float](reaction_rates_p)

    if st80_25_tau_ndx > 0:
        for i in range(1, ncol + 1):
            for k in range(1, troplev[i - 1] + 1):
                reaction_rates[_idx3(i, k, st80_25_tau_ndx, ncol, pver)] = 0.0

def gas_phase_chemdr_compute_relhum_codon(
    ncol: int,
    pver: int,
    h2ovmr_p: cobj,
    satq_p: cobj,
    relhum_p: cobj,
):
    h2ovmr = Ptr[float](h2ovmr_p)
    satq = Ptr[float](satq_p)
    relhum = Ptr[float](relhum_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, ncol)
            value = 0.622 * h2ovmr[idx]
            value = value / satq[idx]
            if value < 0.0:
                value = 0.0
            elif value > 1.0:
                value = 1.0
            relhum[idx] = value

def gas_phase_chemdr_restore_strat_gases_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    hno3_ndx: int,
    h2o_ndx: int,
    delt_inverse: float,
    vmr_p: cobj,
    hno3_gas_p: cobj,
    h2o_gas_p: cobj,
    h2ovmr_p: cobj,
    wrk_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    hno3_gas = Ptr[float](hno3_gas_p)
    h2o_gas = Ptr[float](h2o_gas_p)
    h2ovmr = Ptr[float](h2ovmr_p)
    wrk = Ptr[float](wrk_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, ncol)
            vmr[_idx3(i, k, hno3_ndx, ncol, pver)] = hno3_gas[idx2]
            h2ovmr[idx2] = h2o_gas[idx2]
            vmr[_idx3(i, k, h2o_ndx, ncol, pver)] = h2o_gas[idx2]
            wrk[idx2] = (h2ovmr[idx2] - wrk[idx2]) * delt_inverse

def gas_phase_chemdr_restore_hcl_gas_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    hcl_ndx: int,
    vmr_p: cobj,
    hcl_gas_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    hcl_gas = Ptr[float](hcl_gas_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, ncol)
            vmr[_idx3(i, k, hcl_ndx, ncol, pver)] = hcl_gas[idx2]

def gas_phase_chemdr_init_dust_vmr_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    ndust: int,
    dst_ndx: int,
    vmr_p: cobj,
    dust_vmr_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    dust_vmr = Ptr[float](dust_vmr_p)

    if dst_ndx > 0:
        for m in range(1, ndust + 1):
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    dust_vmr[_idx3(i, k, m, ncol, pver)] = vmr[_idx3(i, k, dst_ndx + m - 1, ncol, pver)]
    else:
        for m in range(1, ndust + 1):
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    dust_vmr[_idx3(i, k, m, ncol, pver)] = 0.0

def gas_phase_chemdr_reset_ste_tracer_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    st80_25_ndx: int,
    pmid_threshold: float,
    st80_vmr: float,
    pmid_p: cobj,
    vmr_p: cobj,
):
    pmid = Ptr[float](pmid_p)
    vmr = Ptr[float](vmr_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            if pmid[_idx2(i, k, pcols)] < pmid_threshold:
                vmr[_idx3(i, k, st80_25_ndx, ncol, pver)] = st80_vmr

def gas_phase_chemdr_zero_sflx_codon(
    pcols: int,
    gas_pcnst: int,
    sflx_p: cobj,
):
    sflx = Ptr[float](sflx_p)

    for m in range(1, gas_pcnst + 1):
        for i in range(1, pcols + 1):
            sflx[_flux_idx(i, m, pcols)] = 0.0

def gas_phase_chemdr_compute_wind_speed_codon(
    ncol: int,
    pcols: int,
    pver: int,
    ufld_p: cobj,
    vfld_p: cobj,
    wind_speed_p: cobj,
):
    ufld = Ptr[float](ufld_p)
    vfld = Ptr[float](vfld_p)
    wind_speed = Ptr[float](wind_speed_p)

    for i in range(1, ncol + 1):
        uval = ufld[_idx2(i, pver, pcols)]
        vval = vfld[_idx2(i, pver, pcols)]
        wind_speed[i - 1] = sqrt(uval * uval + vval * vval)

def gas_phase_chemdr_compute_prect_codon(
    ncol: int,
    pcols: int,
    precc_p: cobj,
    precl_p: cobj,
    prect_p: cobj,
):
    precc = Ptr[float](precc_p)
    precl = Ptr[float](precl_p)
    prect = Ptr[float](prect_p)

    for i in range(1, ncol + 1):
        prect[i - 1] = precc[i - 1] + precl[i - 1]

def gas_phase_chemdr_compute_tvs_codon(
    ncol: int,
    pcols: int,
    pver: int,
    tfld_p: cobj,
    qh2o_p: cobj,
    tvs_p: cobj,
):
    tfld = Ptr[float](tfld_p)
    qh2o = Ptr[float](qh2o_p)
    tvs = Ptr[float](tvs_p)

    for i in range(1, ncol + 1):
        tvs[i - 1] = tfld[_idx2(i, pver, pcols)] * (1.0 + qh2o[_idx2(i, pver, pcols)])

def gas_phase_chemdr_copy_cldw_to_cwat_codon(
    ncol: int,
    pcols: int,
    pver: int,
    cldw_p: cobj,
    cwat_p: cobj,
):
    cldw = Ptr[float](cldw_p)
    cwat = Ptr[float](cwat_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            cwat[_idx2(i, k, ncol)] = cldw[_idx2(i, k, pcols)]

def gas_phase_chemdr_load_h2o_fields_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    h2o_ndx: int,
    mmr_p: cobj,
    vmr_p: cobj,
    qh2o_p: cobj,
    h2ovmr_p: cobj,
):
    mmr = Ptr[float](mmr_p)
    vmr = Ptr[float](vmr_p)
    qh2o = Ptr[float](qh2o_p)
    h2ovmr = Ptr[float](h2ovmr_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            qh2o[_idx2(i, k, pcols)] = mmr[_idx3(i, k, h2o_ndx, pcols, pver)]
            h2ovmr[_idx2(i, k, ncol)] = vmr[_idx3(i, k, h2o_ndx, ncol, pver)]

def gas_phase_chemdr_copy_o3_to_o3s_trop_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    troplev_p: cobj,
    o3_ndx: int,
    o3s_ndx: int,
    vmr_p: cobj,
):
    troplev = Ptr[int](troplev_p)
    vmr = Ptr[float](vmr_p)

    for i in range(1, ncol + 1):
        for k in range(1, troplev[i - 1] + 1):
            vmr[_idx3(i, k, o3s_ndx, ncol, pver)] = vmr[_idx3(i, k, o3_ndx, ncol, pver)]

def gas_phase_chemdr_copy_h2o_to_wrk_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    h2o_ndx: int,
    vmr_p: cobj,
    wrk_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    wrk = Ptr[float](wrk_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            wrk[_idx2(i, k, ncol)] = vmr[_idx3(i, k, h2o_ndx, ncol, pver)]

def gas_phase_chemdr_update_qdsett_wrk_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    h2o_ndx: int,
    delt_inverse: float,
    vmr_p: cobj,
    wrk_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    wrk = Ptr[float](wrk_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, ncol)
            wrk[idx2] = (vmr[_idx3(i, k, h2o_ndx, ncol, pver)] - wrk[idx2]) * delt_inverse

def gas_phase_chemdr_update_qdchem_wrk_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    h2o_ndx: int,
    delt_inverse: float,
    vmr_p: cobj,
    wrk_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    wrk = Ptr[float](wrk_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, ncol)
            wrk[idx2] = (vmr[_idx3(i, k, h2o_ndx, ncol, pver)] - wrk[idx2]) * delt_inverse

def gas_phase_chemdr_init_stratchem_state_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    pcnst: int,
    hno3_ndx: int,
    hcl_ndx: int,
    cldice_ndx: int,
    vmr_p: cobj,
    h2ovmr_p: cobj,
    q_p: cobj,
    hcl_cond_p: cobj,
    hcl_gas_p: cobj,
    hno3_gas_p: cobj,
    h2o_gas_p: cobj,
    wrk_p: cobj,
    cldice_p: cobj,
    hno3_cond_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    h2ovmr = Ptr[float](h2ovmr_p)
    q = Ptr[float](q_p)
    hcl_cond = Ptr[float](hcl_cond_p)
    hcl_gas = Ptr[float](hcl_gas_p)
    hno3_gas = Ptr[float](hno3_gas_p)
    h2o_gas = Ptr[float](h2o_gas_p)
    wrk = Ptr[float](wrk_p)
    cldice = Ptr[float](cldice_p)
    hno3_cond = Ptr[float](hno3_cond_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, ncol)
            hcl_cond[idx2] = 0.0
            hno3_cond[_idx3(i, k, 1, ncol, pver)] = 0.0
            hno3_cond[_idx3(i, k, 2, ncol, pver)] = 0.0
            hno3_gas[idx2] = vmr[_idx3(i, k, hno3_ndx, ncol, pver)]
            h2o_gas[idx2] = h2ovmr[idx2]
            hcl_gas[idx2] = vmr[_idx3(i, k, hcl_ndx, ncol, pver)]
            wrk[idx2] = h2ovmr[idx2]
            cldice[_idx2(i, k, pcols)] = q[_idx3(i, k, cldice_ndx, pcols, pver)]

def gas_phase_chemdr_init_h2so4_gasprod_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    ndx_h2so4: int,
    vmr_p: cobj,
    del_h2so4_gasprod_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    del_h2so4_gasprod = Ptr[float](del_h2so4_gasprod_p)

    if ndx_h2so4 > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                del_h2so4_gasprod[_idx2(i, k, ncol)] = vmr[_idx3(i, k, ndx_h2so4, ncol, pver)]
    else:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                del_h2so4_gasprod[_idx2(i, k, ncol)] = 0.0

def gas_phase_chemdr_store_vmr0_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    vmr_p: cobj,
    vmr0_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    vmr0 = Ptr[float](vmr0_p)

    for m in range(1, gas_pcnst + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx3(i, k, m, ncol, pver)
                vmr0[idx] = vmr[idx]

def gas_phase_chemdr_update_h2so4_gasprod_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    ndx_h2so4: int,
    vmr_p: cobj,
    del_h2so4_gasprod_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    del_h2so4_gasprod = Ptr[float](del_h2so4_gasprod_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, ncol)
            del_h2so4_gasprod[idx2] = vmr[_idx3(i, k, ndx_h2so4, ncol, pver)] - del_h2so4_gasprod[idx2]

def gas_phase_chemdr_reform_hno3_hcl_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    hno3_ndx: int,
    hcl_ndx: int,
    vmr_p: cobj,
    hno3_cond_p: cobj,
    hcl_cond_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    hno3_cond = Ptr[float](hno3_cond_p)
    hcl_cond = Ptr[float](hcl_cond_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            hno3_idx = _idx3(i, k, hno3_ndx, ncol, pver)
            hno3_total = vmr[hno3_idx] + hno3_cond[_idx3(i, k, 1, ncol, pver)]
            vmr[hno3_idx] = hno3_total + hno3_cond[_idx3(i, k, 2, ncol, pver)]

            hcl_idx = _idx3(i, k, hcl_ndx, ncol, pver)
            vmr[hcl_idx] = vmr[hcl_idx] + hcl_cond[_idx2(i, k, ncol)]

def gas_phase_chemdr_stratchem_finalize_batch_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    hno3_ndx: int,
    hcl_ndx: int,
    h2o_ndx: int,
    delt_inverse: float,
    vmr_p: cobj,
    hno3_cond_p: cobj,
    hcl_cond_p: cobj,
    wrk_p: cobj,
):
    gas_phase_chemdr_reform_hno3_hcl_codon(
        ncol, pver, gas_pcnst, hno3_ndx, hcl_ndx, vmr_p, hno3_cond_p, hcl_cond_p
    )
    gas_phase_chemdr_update_qdsett_wrk_codon(
        ncol, pver, gas_pcnst, h2o_ndx, delt_inverse, vmr_p, wrk_p
    )

def gas_phase_chemdr_normalize_extfrc_codon(
    ncol: int,
    pver: int,
    extcnt: int,
    synoz_ndx: int,
    aoa_nh_ext_ndx: int,
    indexm: int,
    extfrc_p: cobj,
    invariants_p: cobj,
):
    extfrc = Ptr[float](extfrc_p)
    invariants = Ptr[float](invariants_p)

    for m in range(1, extcnt + 1):
        if m != synoz_ndx and m != aoa_nh_ext_ndx:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    extfrc[_idx3(i, k, m, ncol, pver)] = extfrc[_idx3(i, k, m, ncol, pver)] / invariants[
                        _idx3(i, k, indexm, ncol, pver)
                    ]

def gas_phase_chemdr_store_drydep_codon(
    ncol: int,
    pcols: int,
    gas_pcnst: int,
    pcnst: int,
    map2chm_p: cobj,
    sflx_p: cobj,
    cflx_p: cobj,
    drydepflx_p: cobj,
):
    map2chm = Ptr[int](map2chm_p)
    sflx = Ptr[float](sflx_p)
    cflx = Ptr[float](cflx_p)
    drydepflx = Ptr[float](drydepflx_p)

    for m in range(1, pcnst + 1):
        for i in range(1, pcols + 1):
            drydepflx[_flux_idx(i, m, pcols)] = 0.0

    for m in range(1, pcnst + 1):
        n = map2chm[m - 1]
        if n > 0:
            for i in range(1, ncol + 1):
                src_idx = _flux_idx(i, n, pcols)
                dst_idx = _flux_idx(i, m, pcols)
                cflx[dst_idx] = cflx[dst_idx] - sflx[src_idx]
                drydepflx[dst_idx] = sflx[src_idx]

@inline
def _gas_phase_chemdr_shell_h2o_setup(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    h2o_ndx: int,
    st80_25_ndx: int,
    aoa_nh_ndx: int,
    nh_5_ndx: int,
    nh_50_ndx: int,
    nh_50w_ndx: int,
    rad2deg: float,
    pmid_p: cobj,
    vmr_p: cobj,
    rlats_p: cobj,
    mmr_p: cobj,
    qh2o_p: cobj,
    h2ovmr_p: cobj,
):
    if st80_25_ndx > 0:
        gas_phase_chemdr_reset_ste_tracer_codon(
            ncol, pcols, pver, gas_pcnst, st80_25_ndx, 80.0e2, 200.0e-9, pmid_p, vmr_p
        )

    if aoa_nh_ndx > 0 and nh_5_ndx > 0 and nh_50_ndx > 0 and nh_50w_ndx > 0:
        rlats = Ptr[float](rlats_p)
        vmr = Ptr[float](vmr_p)
        for j in range(1, ncol + 1):
            xlat = rlats[j - 1] * rad2deg
            if xlat >= 30.0 and xlat <= 50.0:
                vmr[_idx3(j, pver, nh_5_ndx, ncol, pver)] = 100.0e-9
                vmr[_idx3(j, pver, nh_50_ndx, ncol, pver)] = 100.0e-9
                vmr[_idx3(j, pver, nh_50w_ndx, ncol, pver)] = 100.0e-9
                vmr[_idx3(j, pver, aoa_nh_ndx, ncol, pver)] = 0.0

    if h2o_ndx > 0:
        gas_phase_chemdr_load_h2o_fields_codon(
            ncol, pcols, pver, gas_pcnst, h2o_ndx, mmr_p, vmr_p, qh2o_p, h2ovmr_p
        )

@inline
def _gas_phase_chemdr_shell_post_solver(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    ndx_h2so4: int,
    o3_ndx: int,
    o3s_ndx: int,
    delt: float,
    troplev_p: cobj,
    vmr_p: cobj,
    o3s_loss_p: cobj,
    del_h2so4_gasprod_p: cobj,
):
    if o3_ndx > 0 and o3s_ndx > 0:
        gas_phase_chemdr_copy_o3_to_o3s_trop_codon(
            ncol, pcols, pver, gas_pcnst, troplev_p, o3_ndx, o3s_ndx, vmr_p
        )
        troplev = Ptr[int](troplev_p)
        vmr = Ptr[float](vmr_p)
        o3s_loss = Ptr[float](o3s_loss_p)
        for i in range(1, ncol + 1):
            for k in range(troplev[i - 1] + 1, pver + 1):
                idx = _idx3(i, k, o3s_ndx, ncol, pver)
                vmr[idx] = vmr[idx] * exp(-delt * o3s_loss[_idx2(i, k, ncol)])

    if ndx_h2so4 > 0:
        gas_phase_chemdr_update_h2so4_gasprod_codon(
            ncol, pver, gas_pcnst, ndx_h2so4, vmr_p, del_h2so4_gasprod_p
        )

def sulf_interp_codon(
    ncol: int,
    pcols: int,
    pver: int,
    begchunk: int,
    lchnk: int,
    read_sulf_flag: int,
    fields_data_p: cobj,
    ccm_sulf_p: cobj,
):
    ccm_sulf = Ptr[float](ccm_sulf_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            ccm_sulf[(i - 1) + ncol * (k - 1)] = 0.0

    if read_sulf_flag == 0:
        return

    fields_data = Ptr[float](fields_data_p)
    chunk0 = lchnk - begchunk
    chunk_offset = pcols * pver * chunk0
    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            ccm_sulf[(i - 1) + ncol * (k - 1)] = fields_data[chunk_offset + (i - 1) + pcols * (k - 1)]

def gas_phase_chemdr_shell_codon(
    stage: int,
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    pcnst: int,
    rxntot: int,
    extcnt: int,
    nfs: int,
    indexm: int,
    ncol_abs: int,
    jo1d_adj_ndx: int,
    inv_n2_ndx: int,
    inv_o2_ndx: int,
    inv_h2o_ndx: int,
    has_linoz_data_flag: int,
    h2o_ndx: int,
    o2_ndx: int,
    o_ndx: int,
    hno3_ndx: int,
    hcl_ndx: int,
    cldice_ndx: int,
    st80_25_ndx: int,
    aoa_nh_ndx: int,
    nh_5_ndx: int,
    nh_50_ndx: int,
    nh_50w_ndx: int,
    so4_ndx: int,
    st80_25_tau_ndx: int,
    ndx_h2so4: int,
    o3_ndx: int,
    o3s_ndx: int,
    synoz_ndx: int,
    aoa_nh_ext_ndx: int,
    h_ndx: int,
    n_ndx: int,
    elec_ndx: int,
    np_ndx: int,
    n2p_ndx: int,
    op_ndx: int,
    o2p_ndx: int,
    nop_ndx: int,
    fixed_mbar: int,
    rad2deg: float,
    delt: float,
    delt_inverse: float,
    rga: float,
    m2km: float,
    pa2mb: float,
    mwdry: float,
    adv_mass_h2o: float,
    map2chm_p: cobj,
    adv_mass_p: cobj,
    troplev_p: cobj,
    ltrop_sol_p: cobj,
    zen_angle_p: cobj,
    sza_p: cobj,
    phis_p: cobj,
    zi_p: cobj,
    zm_p: cobj,
    pmid_p: cobj,
    zsurf_p: cobj,
    zintr_p: cobj,
    zmidr_p: cobj,
    zmid_p: cobj,
    zint_p: cobj,
    pmb_p: cobj,
    q_p: cobj,
    mmr_p: cobj,
    mbar_p: cobj,
    vmr_p: cobj,
    qh2o_p: cobj,
    h2ovmr_p: cobj,
    rlats_p: cobj,
    sulfate_p: cobj,
    satq_p: cobj,
    relhum_p: cobj,
    cldw_p: cobj,
    cwat_p: cobj,
    extfrc_p: cobj,
    invariants_p: cobj,
    het_rates_p: cobj,
    reaction_rates_p: cobj,
    col_delta_p: cobj,
    col_dens_p: cobj,
    del_h2so4_gasprod_p: cobj,
    vmr0_p: cobj,
    o3s_loss_p: cobj,
    mmr_tend_p: cobj,
    mmr_new_p: cobj,
    qtend_p: cobj,
    tfld_p: cobj,
    tvs_p: cobj,
    sflx_p: cobj,
    ufld_p: cobj,
    vfld_p: cobj,
    wind_speed_p: cobj,
    precc_p: cobj,
    precl_p: cobj,
    prect_p: cobj,
    cflx_p: cobj,
    drydepflx_p: cobj,
    o2mmr_p: cobj,
    ommr_p: cobj,
    hcl_cond_p: cobj,
    hcl_gas_p: cobj,
    hno3_gas_p: cobj,
    h2o_gas_p: cobj,
    wrk_p: cobj,
    cldice_p: cobj,
    hno3_cond_p: cobj,
):
    if stage == 1:
        gas_phase_chemdr_prepare_sza_codon(ncol, rad2deg, zen_angle_p, sza_p)
    elif stage == 2:
        gas_phase_chemdr_prepare_state_codon(
            ncol, pcols, pver, rga, m2km, pa2mb, phis_p, zi_p, zm_p, pmid_p, zsurf_p, zintr_p, zmidr_p,
            zmid_p, zint_p, pmb_p
        )
        gas_phase_chemdr_load_mmr_codon(ncol, pcols, pver, pcnst, map2chm_p, q_p, mmr_p)
    elif stage == 3:
        _gas_phase_chemdr_shell_h2o_setup(
            ncol, pcols, pver, gas_pcnst, h2o_ndx, st80_25_ndx, aoa_nh_ndx, nh_5_ndx, nh_50_ndx,
            nh_50w_ndx, rad2deg, pmid_p, vmr_p, rlats_p, mmr_p, qh2o_p, h2ovmr_p
        )
    elif stage == 4:
        gas_phase_chemdr_zero_sulfate_codon(ncol, pver, sulfate_p)
    elif stage == 5:
        gas_phase_chemdr_load_prognostic_sulfate_codon(ncol, pver, gas_pcnst, so4_ndx, vmr_p, sulfate_p)
    elif stage == 6:
        gas_phase_chemdr_clip_sulfate_codon(ncol, pcols, pver, troplev_p, sulfate_p)
    elif stage == 7:
        gas_phase_chemdr_compute_relhum_codon(ncol, pver, h2ovmr_p, satq_p, relhum_p)
        gas_phase_chemdr_copy_cldw_to_cwat_codon(ncol, pcols, pver, cldw_p, cwat_p)
    elif stage == 8:
        gas_phase_chemdr_normalize_extfrc_codon(
            ncol, pver, extcnt, synoz_ndx, aoa_nh_ext_ndx, indexm, extfrc_p, invariants_p
        )
    elif stage == 9:
        gas_phase_chemdr_zero_het_rates_codon(ncol, pver, gas_pcnst, het_rates_p)
    elif stage == 10:
        gas_phase_chemdr_zero_st80_tau_codon(ncol, pver, rxntot, st80_25_tau_ndx, troplev_p, reaction_rates_p)
    elif stage == 11:
        gas_phase_chemdr_set_ltrop_sol_codon(ncol, has_linoz_data_flag, troplev_p, ltrop_sol_p)
        gas_phase_chemdr_init_h2so4_gasprod_codon(
            ncol, pver, gas_pcnst, ndx_h2so4, vmr_p, del_h2so4_gasprod_p
        )
        gas_phase_chemdr_store_vmr0_codon(ncol, pver, gas_pcnst, vmr_p, vmr0_p)
    elif stage == 12:
        _gas_phase_chemdr_shell_post_solver(
            ncol, pcols, pver, gas_pcnst, ndx_h2so4, o3_ndx, o3s_ndx, delt, troplev_p, vmr_p, o3s_loss_p,
            del_h2so4_gasprod_p
        )
    elif stage == 13:
        gas_phase_chemdr_finalize_tendencies_codon(
            ncol, pcols, pver, gas_pcnst, pcnst, delt_inverse, map2chm_p, mmr_p, mmr_tend_p, mmr_new_p,
            qtend_p
        )
        gas_phase_chemdr_compute_tvs_codon(ncol, pcols, pver, tfld_p, qh2o_p, tvs_p)
        gas_phase_chemdr_zero_sflx_codon(pcols, gas_pcnst, sflx_p)
    elif stage == 14:
        gas_phase_chemdr_compute_wind_speed_codon(ncol, pcols, pver, ufld_p, vfld_p, wind_speed_p)
        gas_phase_chemdr_compute_prect_codon(ncol, pcols, precc_p, precl_p, prect_p)
    elif stage == 15:
        gas_phase_chemdr_store_drydep_codon(
            ncol, pcols, gas_pcnst, pcnst, map2chm_p, sflx_p, cflx_p, drydepflx_p
        )
    elif stage == 16:
        gas_phase_chemdr_init_stratchem_state_codon(
            ncol, pcols, pver, gas_pcnst, pcnst, hno3_ndx, hcl_ndx, cldice_ndx, vmr_p, h2ovmr_p, q_p,
            hcl_cond_p, hcl_gas_p, hno3_gas_p, h2o_gas_p, wrk_p, cldice_p, hno3_cond_p
        )
    elif stage == 17:
        gas_phase_chemdr_restore_strat_gases_codon(
            ncol, pver, gas_pcnst, hno3_ndx, h2o_ndx, delt_inverse, vmr_p, hno3_gas_p, h2o_gas_p, h2ovmr_p,
            wrk_p
        )
    elif stage == 18:
        gas_phase_chemdr_restore_hcl_gas_codon(ncol, pver, gas_pcnst, hcl_ndx, vmr_p, hcl_gas_p)
    elif stage == 19:
        gas_phase_chemdr_update_qdchem_wrk_codon(ncol, pver, gas_pcnst, h2o_ndx, delt_inverse, vmr_p, wrk_p)
    elif stage == 20:
        gas_phase_chemdr_copy_h2o_to_wrk_codon(ncol, pver, gas_pcnst, h2o_ndx, vmr_p, wrk_p)
    elif stage == 21:
        gas_phase_chemdr_stratchem_finalize_batch_codon(
            ncol, pver, gas_pcnst, hno3_ndx, hcl_ndx, h2o_ndx, delt_inverse, vmr_p, hno3_cond_p, hcl_cond_p,
            wrk_p
        )
    elif stage == 22:
        _gas_phase_chemdr_shell_mass_vmr(
            ncol, pcols, pver, gas_pcnst, o2_ndx, o_ndx, h_ndx, n_ndx, fixed_mbar, mwdry, mmr_p, adv_mass_p,
            mbar_p, vmr_p
        )
    elif stage == 23:
        _gas_phase_chemdr_shell_h2o_from_q(ncol, pcols, pver, adv_mass_h2o, q_p, mbar_p, qh2o_p, h2ovmr_p)
    elif stage == 24:
        _gas_phase_chemdr_shell_charge_balance(
            ncol, pver, gas_pcnst, elec_ndx, np_ndx, n2p_ndx, op_ndx, o2p_ndx, nop_ndx, vmr_p, wrk_p
        )
    elif stage == 25:
        gas_phase_chemdr_load_oxygen_mmr_codon(ncol, pcols, pver, o2_ndx, o_ndx, mmr_p, o2mmr_p, ommr_p)
    elif stage == 26:
        _gas_phase_chemdr_shell_adjrxt_setcol(
            ncol, pver, ncol_abs, indexm, reaction_rates_p, invariants_p, col_delta_p, col_dens_p
        )
    elif stage == 27:
        O1D_to_2OH_adj_codon(
            ncol, pcols, pver, rxntot, nfs, jo1d_adj_ndx, inv_n2_ndx, inv_o2_ndx, inv_h2o_ndx,
            reaction_rates_p, invariants_p, tfld_p
        )
    elif stage == 28:
        setrxt_codon(ncol, pcols, pver, tfld_p, reaction_rates_p)
    elif stage == 29:
        negtrc_codon(ncol, pver, gas_pcnst, vmr_p)
    elif stage == 30:
        vmr2mmr_codon(ncol, pcols, pver, gas_pcnst, mbar_p, vmr_p, adv_mass_p, mmr_tend_p)

def set_xnox_photo_codon(
    ncol: int,
    pver: int,
    photos_p: cobj,
    jno2a_ndx: int,
    jno2_ndx: int,
    jn2o5a_ndx: int,
    jn2o5_ndx: int,
    jn2o5b_ndx: int,
    jhno3a_ndx: int,
    jhno3_ndx: int,
    jno3a_ndx: int,
    jno3_ndx: int,
    jho2no2a_ndx: int,
    jho2no2_ndx: int,
    jmpana_ndx: int,
    jmpan_ndx: int,
    jpana_ndx: int,
    jpan_ndx: int,
    jonitra_ndx: int,
    jonitr_ndx: int,
    jo1da_ndx: int,
    jo1d_ndx: int,
    jo3pa_ndx: int,
    jo3p_ndx: int,
):
    photos = Ptr[float](photos_p)

    if jno2a_ndx > 0 and jno2_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                photos[_idx3(i, k, jno2a_ndx, ncol, pver)] = photos[_idx3(i, k, jno2_ndx, ncol, pver)]

    if jn2o5a_ndx > 0 and jn2o5_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                photos[_idx3(i, k, jn2o5a_ndx, ncol, pver)] = photos[_idx3(i, k, jn2o5_ndx, ncol, pver)]

    if jn2o5b_ndx > 0 and jn2o5_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                photos[_idx3(i, k, jn2o5b_ndx, ncol, pver)] = photos[_idx3(i, k, jn2o5_ndx, ncol, pver)]

    if jhno3a_ndx > 0 and jhno3_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                photos[_idx3(i, k, jhno3a_ndx, ncol, pver)] = photos[_idx3(i, k, jhno3_ndx, ncol, pver)]

    if jno3a_ndx > 0 and jno3_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                photos[_idx3(i, k, jno3a_ndx, ncol, pver)] = photos[_idx3(i, k, jno3_ndx, ncol, pver)]

    if jho2no2a_ndx > 0 and jho2no2_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                photos[_idx3(i, k, jho2no2a_ndx, ncol, pver)] = photos[_idx3(i, k, jho2no2_ndx, ncol, pver)]

    if jmpana_ndx > 0 and jmpan_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                photos[_idx3(i, k, jmpana_ndx, ncol, pver)] = photos[_idx3(i, k, jmpan_ndx, ncol, pver)]

    if jpana_ndx > 0 and jpan_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                photos[_idx3(i, k, jpana_ndx, ncol, pver)] = photos[_idx3(i, k, jpan_ndx, ncol, pver)]

    if jonitra_ndx > 0 and jonitr_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                photos[_idx3(i, k, jonitra_ndx, ncol, pver)] = photos[_idx3(i, k, jonitr_ndx, ncol, pver)]

    if jo1da_ndx > 0 and jo1d_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                photos[_idx3(i, k, jo1da_ndx, ncol, pver)] = photos[_idx3(i, k, jo1d_ndx, ncol, pver)]

    if jo3pa_ndx > 0 and jo3p_ndx > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                photos[_idx3(i, k, jo3pa_ndx, ncol, pver)] = photos[_idx3(i, k, jo3p_ndx, ncol, pver)]

def adjrxt_codon(
    ncol: int,
    pver: int,
    rate_p: cobj,
    inv_p: cobj,
    m_p: cobj,
):
    rate = Ptr[float](rate_p)
    inv = Ptr[float](inv_p)
    m = Ptr[float](m_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx_r2 = _idx3(i, k, 2, ncol, pver)
            idx_r3 = _idx3(i, k, 3, ncol, pver)
            idx_r4 = _idx3(i, k, 4, ncol, pver)
            idx_r5 = _idx3(i, k, 5, ncol, pver)
            idx_r6 = _idx3(i, k, 6, ncol, pver)
            idx_r7 = _idx3(i, k, 7, ncol, pver)
            idx_i6 = _idx3(i, k, 6, ncol, pver)
            idx_i7 = _idx3(i, k, 7, ncol, pver)
            idx_i8 = _idx3(i, k, 8, ncol, pver)
            idx_m = _idx2(i, k, ncol)

            inv6 = inv[idx_i6]
            inv7 = inv[idx_i7]
            inv8 = inv[idx_i8]
            im = 1.0 / m[idx_m]

            rate[idx_r3] = rate[idx_r3] * inv6
            rate[idx_r4] = rate[idx_r4] * inv6
            rate[idx_r5] = rate[idx_r5] * inv6
            rate[idx_r6] = rate[idx_r6] * inv6
            rate[idx_r7] = rate[idx_r7] * inv7

            tmp = rate[idx_r2] * inv8
            tmp = tmp * inv8
            rate[idx_r2] = tmp * im

def setrxt_codon(
    ncol: int,
    pcols: int,
    pver: int,
    temp_p: cobj,
    rate_p: cobj,
):
    temp = Ptr[float](temp_p)
    rate = Ptr[float](rate_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            itemp = 1.0 / temp[_idx2(i, k, pcols)]
            rate[_idx3(i, k, 3, ncol, pver)] = 2.9e-12 * exp(-160.0 * itemp)
            rate[_idx3(i, k, 5, ncol, pver)] = 9.6e-12 * exp(-234.0 * itemp)
            rate[_idx3(i, k, 7, ncol, pver)] = 1.9e-13 * exp(520.0 * itemp)

def lu_slv_codon(
    lu_p: cobj,
    b_p: cobj,
):
    lu = Ptr[float](lu_p)
    b = Ptr[float](b_p)

    b[19] = b[19] * lu[21]
    b[18] = b[18] * lu[20]
    b[17] = b[17] * lu[19]
    b[16] = b[16] * lu[18]
    b[15] = b[15] * lu[17]
    b[14] = b[14] * lu[16]
    b[13] = b[13] * lu[15]
    b[12] = b[12] * lu[14]
    b[11] = b[11] * lu[13]
    b[10] = b[10] * lu[12]
    b[9] = b[9] * lu[11]
    b[8] = b[8] * lu[10]
    b[7] = b[7] * lu[9]
    b[6] = b[6] * lu[8]
    b[5] = b[5] * lu[7]
    b[4] = b[4] * lu[6]
    b[3] = b[3] * lu[5]
    b[2] = b[2] - lu[4] * b[3]
    b[2] = b[2] * lu[3]
    b[1] = b[1] - lu[2] * b[2]
    b[1] = b[1] * lu[1]
    b[0] = b[0] * lu[0]

def lu_fac_codon(
    lu_p: cobj,
):
    lu = Ptr[float](lu_p)

    lu[0] = 1.0 / lu[0]
    lu[1] = 1.0 / lu[1]
    lu[3] = 1.0 / lu[3]
    lu[5] = 1.0 / lu[5]
    lu[6] = 1.0 / lu[6]
    lu[7] = 1.0 / lu[7]
    lu[8] = 1.0 / lu[8]
    lu[9] = 1.0 / lu[9]
    lu[10] = 1.0 / lu[10]
    lu[11] = 1.0 / lu[11]
    lu[12] = 1.0 / lu[12]
    lu[13] = 1.0 / lu[13]
    lu[14] = 1.0 / lu[14]
    lu[15] = 1.0 / lu[15]
    lu[16] = 1.0 / lu[16]
    lu[17] = 1.0 / lu[17]
    lu[18] = 1.0 / lu[18]
    lu[19] = 1.0 / lu[19]
    lu[20] = 1.0 / lu[20]
    lu[21] = 1.0 / lu[21]

def linmat_codon(
    mat_p: cobj,
    rxt_p: cobj,
    het_rates_p: cobj,
):
    mat = Ptr[float](mat_p)
    rxt = Ptr[float](rxt_p)
    het_rates = Ptr[float](het_rates_p)

    mat[0] = -(rxt[0] + rxt[2] + het_rates[0])
    mat[1] = -(het_rates[1])
    mat[2] = rxt[3]
    mat[3] = -(rxt[3] + het_rates[2])
    mat[4] = rxt[4] + 0.5 * rxt[5] + rxt[6]
    mat[5] = -(rxt[4] + rxt[5] + rxt[6] + het_rates[3])
    mat[6] = -(het_rates[4])
    mat[7] = -(het_rates[5])
    mat[8] = -(het_rates[6])
    mat[9] = -(het_rates[7])
    mat[10] = -(het_rates[8])
    mat[11] = -(het_rates[9])
    mat[12] = -(het_rates[10])
    mat[13] = -(het_rates[11])
    mat[14] = -(het_rates[12])
    mat[15] = -(het_rates[13])
    mat[16] = -(het_rates[14])
    mat[17] = -(het_rates[15])
    mat[18] = -(het_rates[16])
    mat[19] = -(het_rates[17])
    mat[20] = -(het_rates[18])
    mat[21] = -(het_rates[19])

def imp_prod_loss_codon(
    prod_p: cobj,
    loss_p: cobj,
    y_p: cobj,
    rxt_p: cobj,
    het_rates_p: cobj,
):
    prod = Ptr[float](prod_p)
    loss = Ptr[float](loss_p)
    y = Ptr[float](y_p)
    rxt = Ptr[float](rxt_p)
    het_rates = Ptr[float](het_rates_p)

    loss[0] = (rxt[0] + rxt[2] + het_rates[0]) * y[0]
    prod[0] = 0.0
    loss[1] = het_rates[1] * y[1]
    prod[1] = rxt[3] * y[2]
    loss[2] = (rxt[3] + het_rates[2]) * y[2]
    prod[2] = (rxt[4] + 0.5 * rxt[5] + rxt[6]) * y[3]
    loss[3] = (rxt[4] + rxt[5] + rxt[6] + het_rates[3]) * y[3]
    prod[3] = 0.0
    loss[4] = het_rates[4] * y[4]
    prod[4] = 0.0
    loss[5] = het_rates[5] * y[5]
    prod[5] = 0.0
    loss[6] = het_rates[6] * y[6]
    prod[6] = 0.0
    loss[7] = het_rates[7] * y[7]
    prod[7] = 0.0
    loss[8] = het_rates[8] * y[8]
    prod[8] = 0.0
    loss[9] = het_rates[9] * y[9]
    prod[9] = 0.0
    loss[10] = het_rates[10] * y[10]
    prod[10] = 0.0
    loss[11] = het_rates[11] * y[11]
    prod[11] = 0.0
    loss[12] = het_rates[12] * y[12]
    prod[12] = 0.0
    loss[13] = het_rates[13] * y[13]
    prod[13] = 0.0
    loss[14] = het_rates[14] * y[14]
    prod[14] = 0.0
    loss[15] = het_rates[15] * y[15]
    prod[15] = 0.0
    loss[16] = het_rates[16] * y[16]
    prod[16] = 0.0
    loss[17] = het_rates[17] * y[17]
    prod[17] = 0.0
    loss[18] = het_rates[18] * y[18]
    prod[18] = 0.0
    loss[19] = het_rates[19] * y[19]
    prod[19] = 0.0

def nlnmat_codon(
    mat_p: cobj,
    lmat_p: cobj,
    dti: float,
):
    mat = Ptr[float](mat_p)
    lmat = Ptr[float](lmat_p)

    mat[0] = lmat[0]
    mat[1] = lmat[1]
    mat[2] = lmat[2]
    mat[3] = lmat[3]
    mat[4] = lmat[4]
    mat[5] = lmat[5]
    mat[6] = lmat[6]
    mat[7] = lmat[7]
    mat[8] = lmat[8]
    mat[9] = lmat[9]
    mat[10] = lmat[10]
    mat[11] = lmat[11]
    mat[12] = lmat[12]
    mat[13] = lmat[13]
    mat[14] = lmat[14]
    mat[15] = lmat[15]
    mat[16] = lmat[16]
    mat[17] = lmat[17]
    mat[18] = lmat[18]
    mat[19] = lmat[19]
    mat[20] = lmat[20]
    mat[21] = lmat[21]

    mat[0] = mat[0] - dti
    mat[1] = mat[1] - dti
    mat[3] = mat[3] - dti
    mat[5] = mat[5] - dti
    mat[6] = mat[6] - dti
    mat[7] = mat[7] - dti
    mat[8] = mat[8] - dti
    mat[9] = mat[9] - dti
    mat[10] = mat[10] - dti
    mat[11] = mat[11] - dti
    mat[12] = mat[12] - dti
    mat[13] = mat[13] - dti
    mat[14] = mat[14] - dti
    mat[15] = mat[15] - dti
    mat[16] = mat[16] - dti
    mat[17] = mat[17] - dti
    mat[18] = mat[18] - dti
    mat[19] = mat[19] - dti
    mat[20] = mat[20] - dti
    mat[21] = mat[21] - dti

def indprd_codon(
    class_id: int,
    ncol: int,
    pver: int,
    nprod: int,
    rxt_p: cobj,
    extfrc_p: cobj,
    prod_p: cobj,
):
    rxt = Ptr[float](rxt_p)
    extfrc = Ptr[float](extfrc_p)
    prod = Ptr[float](prod_p)

    if class_id != 4:
        return

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            prod[_idx3(i, k, 1, ncol, pver)] = rxt[_idx3(i, k, 2, ncol, pver)]
            prod[_idx3(i, k, 2, ncol, pver)] = 0.0
            prod[_idx3(i, k, 3, ncol, pver)] = extfrc[_idx3(i, k, 1, ncol, pver)]
            prod[_idx3(i, k, 4, ncol, pver)] = 0.0
            prod[_idx3(i, k, 5, ncol, pver)] = 0.0
            prod[_idx3(i, k, 6, ncol, pver)] = extfrc[_idx3(i, k, 2, ncol, pver)]
            prod[_idx3(i, k, 7, ncol, pver)] = extfrc[_idx3(i, k, 4, ncol, pver)]
            prod[_idx3(i, k, 8, ncol, pver)] = 0.0
            prod[_idx3(i, k, 9, ncol, pver)] = extfrc[_idx3(i, k, 5, ncol, pver)]
            prod[_idx3(i, k, 10, ncol, pver)] = 0.0
            prod[_idx3(i, k, 11, ncol, pver)] = 0.0
            prod[_idx3(i, k, 12, ncol, pver)] = extfrc[_idx3(i, k, 6, ncol, pver)]
            prod[_idx3(i, k, 13, ncol, pver)] = extfrc[_idx3(i, k, 3, ncol, pver)]
            prod[_idx3(i, k, 14, ncol, pver)] = 0.0
            prod[_idx3(i, k, 15, ncol, pver)] = 0.0
            prod[_idx3(i, k, 16, ncol, pver)] = extfrc[_idx3(i, k, 7, ncol, pver)]
            prod[_idx3(i, k, 17, ncol, pver)] = 0.0
            prod[_idx3(i, k, 18, ncol, pver)] = 0.0
            prod[_idx3(i, k, 19, ncol, pver)] = 0.0
            prod[_idx3(i, k, 20, ncol, pver)] = 0.0

def negtrc_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    fld_p: cobj,
):
    fld = Ptr[float](fld_p)

    for m in range(1, gas_pcnst + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx3(i, k, m, ncol, pver)
                if fld[idx] < 0.0:
                    fld[idx] = 0.0

def O1D_to_2OH_adj_codon(
    ncol: int,
    pcols: int,
    pver: int,
    rxntot: int,
    nfs: int,
    jo1d_ndx: int,
    n2_ndx: int,
    o2_ndx: int,
    h2o_ndx: int,
    p_rate_p: cobj,
    inv_p: cobj,
    tfld_p: cobj,
):
    if jo1d_ndx < 1:
        return

    p_rate = Ptr[float](p_rate_p)
    inv = Ptr[float](inv_p)
    tfld = Ptr[float](tfld_p)

    x1 = 2.15e-11
    x2 = 3.30e-11
    x3 = 1.63e-10
    y1 = 110.0
    y2 = 55.0
    y3 = 60.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            temp = tfld[_idx2(i, k, pcols)]
            n2_rate = x1 * exp(y1 / temp) * inv[_idx3(i, k, n2_ndx, ncol, pver)]
            o2_rate = x2 * exp(y2 / temp) * inv[_idx3(i, k, o2_ndx, ncol, pver)]
            h2o_rate = x3 * exp(y3 / temp) * inv[_idx3(i, k, h2o_ndx, ncol, pver)]
            denom = h2o_rate + n2_rate + o2_rate
            p_rate[_idx3(i, k, jo1d_ndx, ncol, pver)] = (
                p_rate[_idx3(i, k, jo1d_ndx, ncol, pver)] * (h2o_rate / denom)
            )
