from math import exp, log


KABSL = 0.090361
IC_LIMIT = 1.0e-12


@export
def conv_water_readnl_codon(value: float) -> float:
    return value


@export
def conv_water_register_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def conv_water_init_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@inline
def _idx(i: int, k: int, pcols: int) -> int:
    """Fortran arrays declared as (pcols, pver)."""
    return (i - 1) + (k - 1) * pcols


@export
def conv_water_4rad_codon(
    ncol: int,
    pcols: int,
    pver: int,
    conv_water_mode: int,
    microp_is_rk: int,
    frac_limit: float,
    gravit: float,
    pdel_p: cobj,
    ls_liq_p: cobj,
    ls_ice_p: cobj,
    ast_p: cobj,
    sh_frac_p: cobj,
    dp_frac_p: cobj,
    rei_p: cobj,
    dp_icwmr_p: cobj,
    sh_icwmr_p: cobj,
    fice_p: cobj,
    totg_liq_p: cobj,
    totg_ice_p: cobj,
    conv_ice_p: cobj,
    conv_liq_p: cobj,
    tot_ice_p: cobj,
    tot_liq_p: cobj,
    totg_ice_sh_p: cobj,
    totg_liq_sh_p: cobj,
    totg_ice_dp_p: cobj,
    totg_liq_dp_p: cobj,
    fresh_p: cobj,
    fredp_p: cobj,
    frecu_p: cobj,
    fretot_p: cobj,
    sh_cldliq_p: cobj,
    sh_cldice_p: cobj,
):
    pdel = Ptr[float](pdel_p)
    ls_liq = Ptr[float](ls_liq_p)
    ls_ice = Ptr[float](ls_ice_p)
    ast = Ptr[float](ast_p)
    sh_frac = Ptr[float](sh_frac_p)
    dp_frac = Ptr[float](dp_frac_p)
    rei = Ptr[float](rei_p)
    dp_icwmr = Ptr[float](dp_icwmr_p)
    sh_icwmr = Ptr[float](sh_icwmr_p)
    fice = Ptr[float](fice_p)
    totg_liq = Ptr[float](totg_liq_p)
    totg_ice = Ptr[float](totg_ice_p)
    conv_ice = Ptr[float](conv_ice_p)
    conv_liq = Ptr[float](conv_liq_p)
    tot_ice = Ptr[float](tot_ice_p)
    tot_liq = Ptr[float](tot_liq_p)
    totg_ice_sh = Ptr[float](totg_ice_sh_p)
    totg_liq_sh = Ptr[float](totg_liq_sh_p)
    totg_ice_dp = Ptr[float](totg_ice_dp_p)
    totg_liq_dp = Ptr[float](totg_liq_dp_p)
    fresh = Ptr[float](fresh_p)
    fredp = Ptr[float](fredp_p)
    frecu = Ptr[float](frecu_p)
    fretot = Ptr[float](fretot_p)
    sh_cldliq = Ptr[float](sh_cldliq_p)
    sh_cldice = Ptr[float](sh_cldice_p)

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            idx = _idx(i, k, pcols)
            fresh[idx] = 0.0
            fredp[idx] = 0.0
            frecu[idx] = 0.0
            fretot[idx] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx(i, k, pcols)

            sh0_frac = 0.0
            if not (sh_frac[idx] <= frac_limit or sh_icwmr[idx] <= IC_LIMIT):
                sh0_frac = sh_frac[idx]

            dp0_frac = 0.0
            if not (dp_frac[idx] <= frac_limit or dp_icwmr[idx] <= IC_LIMIT):
                dp0_frac = dp_frac[idx]

            cu0_frac = sh0_frac + dp0_frac
            wrk1 = min(1.0, max(0.0, ls_ice[idx] / (ls_ice[idx] + ls_liq[idx] + 1.0e-36)))

            if cu0_frac < frac_limit or (sh_icwmr[idx] + dp_icwmr[idx]) < IC_LIMIT:
                cu0_frac = 0.0
                cu_icwmr = 0.0

                ls_frac = ast[idx]
                if ls_frac < frac_limit:
                    ls_frac = 0.0
                    ls_icwmr = 0.0
                else:
                    ls_icwmr = (ls_liq[idx] + ls_ice[idx]) / max(frac_limit, ls_frac)

                tot0_frac = ls_frac
                tot_icwmr = ls_icwmr
            else:
                if microp_is_rk == 1:
                    kabsi = 0.005 + 1.0 / rei[idx]
                else:
                    kabsi = 0.005 + 1.0 / min(max(13.0, rei[idx]), 130.0)

                kabs = KABSL * (1.0 - wrk1) + kabsi * wrk1
                alpha = -1.66 * kabs * pdel[idx] / gravit * 1000.0

                if conv_water_mode == 1:
                    cu_icwmr = (sh0_frac * sh_icwmr[idx] + dp0_frac * dp_icwmr[idx]) / max(frac_limit, cu0_frac)
                else:
                    sh0 = exp(alpha * sh_icwmr[idx])
                    dp0 = exp(alpha * dp_icwmr[idx])
                    cu_icwmr = log((sh0_frac * sh0 + dp0_frac * dp0) / max(frac_limit, cu0_frac))
                    cu_icwmr = cu_icwmr / alpha

                ls_frac = ast[idx]
                ls_icwmr = (ls_liq[idx] + ls_ice[idx]) / max(frac_limit, ls_frac)
                tot0_frac = ls_frac + cu0_frac

                if conv_water_mode == 1:
                    tot_icwmr = (ls_frac * ls_icwmr + cu0_frac * cu_icwmr) / max(frac_limit, tot0_frac)
                else:
                    tot_icwmr = log(
                        (ls_frac * exp(alpha * ls_icwmr) + cu0_frac * exp(alpha * cu_icwmr))
                        / max(frac_limit, tot0_frac)
                    )
                    tot_icwmr = tot_icwmr / alpha

            conv_ice[idx] = cu_icwmr * wrk1
            conv_liq[idx] = cu_icwmr * (1.0 - wrk1)
            tot_ice[idx] = tot_icwmr * wrk1
            tot_liq[idx] = tot_icwmr * (1.0 - wrk1)
            totg_ice[idx] = tot0_frac * tot_icwmr * wrk1
            totg_liq[idx] = tot0_frac * tot_icwmr * (1.0 - wrk1)

            totg_ice_sh[idx] = sh0_frac * sh_icwmr[idx] * wrk1
            totg_ice_dp[idx] = dp0_frac * dp_icwmr[idx] * wrk1
            totg_liq_sh[idx] = sh0_frac * sh_icwmr[idx] * (1.0 - wrk1)
            totg_liq_dp[idx] = dp0_frac * dp_icwmr[idx] * (1.0 - wrk1)

            if sh0_frac > frac_limit:
                fresh[idx] = 1.0
            if dp0_frac > frac_limit:
                fredp[idx] = 1.0
            if cu0_frac > frac_limit:
                frecu[idx] = 1.0
            if tot0_frac > frac_limit:
                fretot[idx] = 1.0

            sh_cldliq[idx] = sh_icwmr[idx] * (1.0 - fice[idx]) * sh_frac[idx]
            sh_cldice[idx] = sh_icwmr[idx] * fice[idx] * sh_frac[idx]
