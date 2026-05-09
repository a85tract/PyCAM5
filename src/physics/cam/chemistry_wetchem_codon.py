from math import exp, gamma, log10, sqrt
from chemistry_common_codon import _idx2, _idx3
from C import neu_wetdep_dempirical_native_cb(float, float) -> float

def neu_wetdep_aux_prepare_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    gas_cnt: int,
    index_cldice: int,
    index_cldliq: int,
    gravit: float,
    mapping_to_mmr_p: cobj,
    area_p: cobj,
    mmr_p: cobj,
    pmid_p: cobj,
    pdel_p: cobj,
    zint_p: cobj,
    tfld_p: cobj,
    prain_p: cobj,
    nevapr_p: cobj,
    cld_p: cobj,
    cmfdqr_p: cobj,
    mass_in_layer_p: cobj,
    cldice_p: cobj,
    cldliq_p: cobj,
    cldfrc_p: cobj,
    totprec_p: cobj,
    totevap_p: cobj,
    delz_p: cobj,
    delp_p: cobj,
    p_p: cobj,
    rls_p: cobj,
    evaprate_p: cobj,
    temp_p: cobj,
    trc_mass_p: cobj,
    dtwr_p: cobj,
):
    mapping_to_mmr = Ptr[int](mapping_to_mmr_p)
    area = Ptr[float](area_p)
    mmr = Ptr[float](mmr_p)
    pmid = Ptr[float](pmid_p)
    pdel = Ptr[float](pdel_p)
    zint = Ptr[float](zint_p)
    tfld = Ptr[float](tfld_p)
    prain = Ptr[float](prain_p)
    nevapr = Ptr[float](nevapr_p)
    cld = Ptr[float](cld_p)
    cmfdqr = Ptr[float](cmfdqr_p)
    mass_in_layer = Ptr[float](mass_in_layer_p)
    cldice = Ptr[float](cldice_p)
    cldliq = Ptr[float](cldliq_p)
    cldfrc = Ptr[float](cldfrc_p)
    totprec = Ptr[float](totprec_p)
    totevap = Ptr[float](totevap_p)
    delz = Ptr[float](delz_p)
    delp = Ptr[float](delp_p)
    press = Ptr[float](p_p)
    rls = Ptr[float](rls_p)
    evaprate = Ptr[float](evaprate_p)
    temp = Ptr[float](temp_p)
    trc_mass = Ptr[float](trc_mass_p)
    dtwr = Ptr[float](dtwr_p)

    for k in range(1, pver + 1):
        kk = pver - k + 1
        for i in range(1, ncol + 1):
            idx_rev_ncol = _idx2(i, k, ncol)
            idx_kk_pcols = _idx2(i, kk, pcols)
            layer_mass = area[i - 1] * pdel[idx_kk_pcols] / gravit
            mass_in_layer[idx_rev_ncol] = layer_mass

            cldice[idx_rev_ncol] = mmr[_idx3(i, kk, index_cldice, pcols, pver)]
            cldliq[idx_rev_ncol] = mmr[_idx3(i, kk, index_cldliq, pcols, pver)]
            cldfrc[idx_rev_ncol] = cld[_idx2(i, kk, ncol)]

            totprec[idx_rev_ncol] = (prain[_idx2(i, kk, ncol)] + cmfdqr[_idx2(i, kk, ncol)]) * layer_mass
            totevap[idx_rev_ncol] = nevapr[_idx2(i, kk, ncol)] * layer_mass

            delz[idx_rev_ncol] = zint[_idx2(i, kk, pcols)] - zint[_idx2(i, kk + 1, pcols)]
            temp[idx_rev_ncol] = tfld[idx_kk_pcols]

            for m in range(1, gas_cnt + 1):
                spc = mapping_to_mmr[m - 1]
                trc_mass[_idx3(i, k, m, ncol, pver)] = mmr[_idx3(i, kk, spc, pcols, pver)] * layer_mass

            delp[idx_rev_ncol] = pdel[idx_kk_pcols] * 0.01
            press[idx_rev_ncol] = pmid[idx_kk_pcols] * 0.01

    for m in range(1, gas_cnt + 1):
        spc = mapping_to_mmr[m - 1]
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                dtwr[_idx3(i, k, m, ncol, pver)] = mmr[_idx3(i, k, spc, pcols, pver)]

    for i in range(1, ncol + 1):
        rls[_idx2(i, pver, ncol)] = 0.0
        evaprate[_idx2(i, pver, ncol)] = 0.0

    for k in range(pver - 1, 0, -1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, ncol)
            next_idx = _idx2(i, k + 1, ncol)
            rls[idx] = max(0.0, totprec[idx] - totevap[idx] + rls[next_idx])
            evaprate[idx] = min(1.0, totevap[idx] / (rls[next_idx] + 1.0e-36))

def neu_wetdep_aux_finish_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    gas_cnt: int,
    delt: float,
    pi: float,
    mapping_to_mmr_p: cobj,
    lats_p: cobj,
    pmid_p: cobj,
    mass_in_layer_p: cobj,
    trc_mass_p: cobj,
    dtwr_p: cobj,
    wd_mmr_p: cobj,
    wd_tend_p: cobj,
):
    mapping_to_mmr = Ptr[int](mapping_to_mmr_p)
    lats = Ptr[float](lats_p)
    pmid = Ptr[float](pmid_p)
    mass_in_layer = Ptr[float](mass_in_layer_p)
    trc_mass = Ptr[float](trc_mass_p)
    dtwr = Ptr[float](dtwr_p)
    wd_mmr = Ptr[float](wd_mmr_p)
    wd_tend = Ptr[float](wd_tend_p)

    for k in range(1, pver + 1):
        kk = pver - k + 1
        for i in range(1, ncol + 1):
            layer_mass = mass_in_layer[_idx2(i, k, ncol)]
            for m in range(1, gas_cnt + 1):
                wd_mmr[_idx3(i, kk, m, ncol, pver)] = trc_mass[_idx3(i, k, m, ncol, pver)] / layer_mass

    for m in range(1, gas_cnt + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx3(i, k, m, ncol, pver)
                dtwr[idx] = (wd_mmr[idx] - dtwr[idx]) / delt

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            if abs(lats[i - 1] * 180.0 / pi) > 60.0:
                if pmid[_idx2(i, k, pcols)] < 20000.0:
                    for m in range(1, gas_cnt + 1):
                        dtwr[_idx3(i, k, m, ncol, pver)] = 0.0

    for m in range(1, gas_cnt + 1):
        spc = mapping_to_mmr[m - 1]
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                wd_tend[_idx3(i, k, spc, pcols, pver)] += dtwr[_idx3(i, k, m, ncol, pver)]

def neu_wetdep_henry_flags_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_cnt: int,
    nh3_ndx: int,
    co2_ndx: int,
    t0: float,
    ph: float,
    ph_inv: float,
    mapping_to_heff_p: cobj,
    dheff_p: cobj,
    tfld_p: cobj,
    heff_p: cobj,
    wrk_p: cobj,
    dk1s_p: cobj,
    dk2s_p: cobj,
    tckaqb_p: cobj,
):
    mapping_to_heff = Ptr[int](mapping_to_heff_p)
    dheff = Ptr[float](dheff_p)
    tfld = Ptr[float](tfld_p)
    heff = Ptr[float](heff_p)
    wrk = Ptr[float](wrk_p)
    dk1s = Ptr[float](dk1s_p)
    dk2s = Ptr[float](dk2s_p)
    tckaqb = Ptr[int](tckaqb_p)

    # Fortran declarations: tfld(pcols,pver), heff(ncol,pver,gas_wetdep_cnt).
    for m in range(1, gas_cnt + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                heff[_idx3(i, k, m, ncol, pver)] = 0.0

    for k in range(1, pver + 1):
        kk = pver - k + 1

        for i in range(1, ncol + 1):
            temp = tfld[_idx2(i, kk, pcols)]
            wrk[i - 1] = (t0 - temp) / (t0 * temp)

        for m in range(1, gas_cnt + 1):
            l = mapping_to_heff[m - 1]
            base = 6 * (l - 1)
            e298 = dheff[base]
            dhr = dheff[base + 1]

            for i in range(1, ncol + 1):
                heff[_idx3(i, k, m, ncol, pver)] = e298 * exp(dhr * wrk[i - 1])

            if dheff[base + 2] != 0.0 and dheff[base + 4] == 0.0:
                e298 = dheff[base + 2]
                dhr = dheff[base + 3]
                for i in range(1, ncol + 1):
                    dk1s[i - 1] = e298 * exp(dhr * wrk[i - 1])

                for i in range(1, ncol + 1):
                    idx = _idx3(i, k, m, ncol, pver)
                    if heff[idx] != 0.0:
                        heff[idx] = heff[idx] * (1.0 + dk1s[i - 1] * ph_inv)
                    else:
                        heff[idx] = dk1s[i - 1] * ph_inv

            if dheff[base + 4] != 0.0:
                if nh3_ndx > 0 or co2_ndx > 0:
                    e298 = dheff[base + 2]
                    dhr = dheff[base + 3]
                    for i in range(1, ncol + 1):
                        dk1s[i - 1] = e298 * exp(dhr * wrk[i - 1])

                    e298 = dheff[base + 4]
                    dhr = dheff[base + 5]
                    for i in range(1, ncol + 1):
                        dk2s[i - 1] = e298 * exp(dhr * wrk[i - 1])

                    if m == co2_ndx:
                        for i in range(1, ncol + 1):
                            idx = _idx3(i, k, m, ncol, pver)
                            heff[idx] = heff[idx] * (1.0 + dk1s[i - 1] * ph_inv) * (
                                1.0 + dk2s[i - 1] * ph_inv
                            )
                    elif m == nh3_ndx:
                        for i in range(1, ncol + 1):
                            idx = _idx3(i, k, m, ncol, pver)
                            heff[idx] = heff[idx] * (1.0 + dk1s[i - 1] * ph / dk2s[i - 1])

    for m in range(1, gas_cnt + 1):
        max_heff = heff[_idx3(1, 1, m, ncol, pver)]
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                val = heff[_idx3(i, k, m, ncol, pver)]
                if val > max_heff:
                    max_heff = val

        if max_heff > 1.0e4:
            tckaqb[m - 1] = 1
        else:
            tckaqb[m - 1] = 0

def neu_wetdep_prepare_henry_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    gas_cnt: int,
    index_cldice: int,
    index_cldliq: int,
    nh3_ndx: int,
    co2_ndx: int,
    gravit: float,
    t0: float,
    ph: float,
    ph_inv: float,
    mapping_to_mmr_p: cobj,
    mapping_to_heff_p: cobj,
    area_p: cobj,
    mmr_p: cobj,
    pmid_p: cobj,
    pdel_p: cobj,
    zint_p: cobj,
    tfld_p: cobj,
    prain_p: cobj,
    nevapr_p: cobj,
    cld_p: cobj,
    cmfdqr_p: cobj,
    dheff_p: cobj,
    mass_in_layer_p: cobj,
    cldice_p: cobj,
    cldliq_p: cobj,
    cldfrc_p: cobj,
    totprec_p: cobj,
    totevap_p: cobj,
    delz_p: cobj,
    delp_p: cobj,
    p_p: cobj,
    rls_p: cobj,
    evaprate_p: cobj,
    temp_p: cobj,
    trc_mass_p: cobj,
    dtwr_p: cobj,
    heff_p: cobj,
    wrk_p: cobj,
    dk1s_p: cobj,
    dk2s_p: cobj,
    tckaqb_p: cobj,
):
    neu_wetdep_aux_prepare_codon(
        ncol,
        pcols,
        pver,
        pcnst,
        gas_cnt,
        index_cldice,
        index_cldliq,
        gravit,
        mapping_to_mmr_p,
        area_p,
        mmr_p,
        pmid_p,
        pdel_p,
        zint_p,
        tfld_p,
        prain_p,
        nevapr_p,
        cld_p,
        cmfdqr_p,
        mass_in_layer_p,
        cldice_p,
        cldliq_p,
        cldfrc_p,
        totprec_p,
        totevap_p,
        delz_p,
        delp_p,
        p_p,
        rls_p,
        evaprate_p,
        temp_p,
        trc_mass_p,
        dtwr_p,
    )
    neu_wetdep_henry_flags_codon(
        ncol,
        pcols,
        pver,
        gas_cnt,
        nh3_ndx,
        co2_ndx,
        t0,
        ph,
        ph_inv,
        mapping_to_heff_p,
        dheff_p,
        tfld_p,
        heff_p,
        wrk_p,
        dk1s_p,
        dk2s_p,
        tckaqb_p,
    )

@inline
def _neu_wetdep_disgas_core(
    clwx: float,
    cfx: float,
    molmass: float,
    hstar: float,
    tm: float,
    pr: float,
    qm: float,
    qt: float,
) -> float:
    tmix = 258.0
    reteff = 0.5

    if tm >= 263.0:
        return (hstar * (qt / (qm * cfx)) * 0.029 * (pr / 1.0e3)) * (clwx * qm)
    elif tm <= tmix:
        muemp = exp(-14.2252 + (1.55704e-1 * tm) - (7.1929e-4 * (tm ** 2.0)))
        return muemp * (molmass / 18.0) * (clwx * qm)

    return reteff * ((hstar * (qt / (qm * cfx)) * 0.029 * (pr / 1.0e3)) * (clwx * qm))

@inline
def _neu_wetdep_raingas_core(
    rrain: float,
    dtscav: float,
    clwx: float,
    cfx: float,
    qm: float,
    qt: float,
    qtdis: float,
) -> float:
    qtdisstar = (qtdis * (qt * cfx)) / (qtdis + (qt * cfx))
    qtlf = (rrain * qtdisstar) / (clwx * qm * qt * cfx)
    return qt * cfx * (1.0 - exp(-dtscav * qtlf))

@inline
def _neu_wetdep_dempirical_core(cwater: float, rrate: float) -> float:
    rratex = rrate * 3600.0
    wx = cwater * 1.0e3

    if rratex > 0.04:
        theta = exp(-1.43 * log10(7.0 * rratex)) + 2.8
    else:
        theta = 5.0

    phi = rratex / (3600.0 * 10.0)
    eta = exp((3.01 * theta) - 10.5)
    beta = theta / (1.0 + 0.638)
    alpha = exp(4.0 * (beta - 3.5))
    bee = (0.638 * theta / (1.0 + 0.638)) - 1.0
    gamtheta = gamma(theta)
    gambeta = gamma(beta + 1.0)
    return (((wx * eta * gamtheta) / (1.0e6 * alpha * phi * gambeta)) ** (-1.0 / bee)) * 10.0

@inline
def _neu_wetdep_dempirical_eval(cwater: float, rrate: float, dempirical_impl: int) -> float:
    if dempirical_impl == 0:
        return neu_wetdep_dempirical_native_cb(cwater, rrate)

    return _neu_wetdep_dempirical_core(cwater, rrate)

def neu_wetdep_dempirical_codon(
    cwater: float,
    rrate: float,
    dempirical_p: cobj,
):
    dempirical = Ptr[float](dempirical_p)
    dempirical[0] = _neu_wetdep_dempirical_core(cwater, rrate)

@inline
def _neu_wetdep_washgas_core(
    rwash: float,
    boxf: float,
    dtscav: float,
    qtrtop: float,
    hstar: float,
    tm: float,
    pr: float,
    qm: float,
    qt: float,
    qtwash_p: Ptr[float],
    qtevap_p: Ptr[float],
):
    if boxf == 0.0:
        qtwash_p[0] = 0.0
        qtevap_p[0] = 0.0
        return

    fwash = (rwash * hstar * 29.0e-6 * pr) / (qm * boxf)
    qtmax = qt * fwash * dtscav

    if qtmax > qtrtop:
        qtdif = min(qt, qtmax - qtrtop)
        qtwash_p[0] = qtdif * (1.0 - exp(-dtscav * fwash))
        qtevap_p[0] = 0.0
    else:
        qtwash_p[0] = 0.0
        qtevap_p[0] = qtrtop - qtmax

@inline
def _neu_wetdep_new_precip_scavenging(
    scavenging_active: int,
    rprecip: float,
    garea: float,
    dtscav: float,
    clwx: float,
    cfxx_l: float,
    tcmass_n: float,
    hstar_ln: float,
    tem_l: float,
    pofl_l: float,
    qm_l: float,
    qtt_l: float,
    fcxa: float,
    fcxb: float,
    qtraincxa_p: Ptr[float],
    qtraincxb_p: Ptr[float],
):
    if rprecip > 0.0:
        if scavenging_active != 0:
            rrain = rprecip * garea
            qtdiscf = _neu_wetdep_disgas_core(
                clwx,
                cfxx_l,
                tcmass_n,
                hstar_ln,
                tem_l,
                pofl_l,
                qm_l,
                qtt_l * cfxx_l,
            )
            qtrain = _neu_wetdep_raingas_core(
                rrain,
                dtscav,
                clwx,
                cfxx_l,
                qm_l,
                qtt_l,
                qtdiscf,
            )
            wrk = qtrain / cfxx_l
            qtraincxa_p[0] = fcxa * wrk
            qtraincxb_p[0] = fcxb * wrk
        else:
            qtraincxa_p[0] = 0.0
            qtraincxb_p[0] = 0.0
    else:
        qtraincxa_p[0] = 0.0
        qtraincxb_p[0] = 0.0

@inline
def _neu_wetdep_ice_riming_scavenging(
    scavenging_active: int,
    tem_l: float,
    tfroz: float,
    rhosnowfix: float,
    coleffsnow: float,
    dca: float,
    rca: float,
    qtt_l: float,
    fcxa: float,
    clwx: float,
    cfxx_l: float,
    tcmass_n: float,
    hstar_ln: float,
    pofl_l: float,
    qm_l: float,
    rnew: float,
    garea: float,
    dtscav: float,
    qtrimecxa_p: Ptr[float],
):
    if scavenging_active != 0:
        if tem_l <= tfroz:
            rhosnow = rhosnowfix
        else:
            rhosnow = 0.303 * (tem_l - tfroz) * rhosnowfix

        qtcxa = qtt_l * fcxa
        qtdisrime = _neu_wetdep_disgas_core(
            clwx * (fcxa / cfxx_l),
            fcxa,
            tcmass_n,
            hstar_ln,
            tem_l,
            pofl_l,
            qm_l,
            qtcxa,
        )
        qtdisstar = (qtdisrime * qtcxa) / (qtdisrime + qtcxa)
        qtrimecxa_p[0] = qtcxa * (
            1.0
            - exp(
                (-coleffsnow / (dca * 1.0e-3))
                * (rca / (2.0 * rhosnow))
                * (qtdisstar / qtcxa)
                * dtscav
            )
        )
        qtrimecxa_p[0] = min(
            qtrimecxa_p[0],
            ((rnew * garea * dtscav) / (clwx * qm_l * (fcxa / cfxx_l))) * qtdisstar,
        )
    else:
        qtrimecxa_p[0] = 0.0

@inline
def _neu_wetdep_rain_riming_scavenging(
    scavenging_active: int,
    coleffrain: float,
    rca: float,
    qtt_l: float,
    fcxa: float,
    clwx: float,
    cfxx_l: float,
    tcmass_n: float,
    hstar_ln: float,
    tem_l: float,
    pofl_l: float,
    qm_l: float,
    rnew: float,
    garea: float,
    dtscav: float,
    qtdisrime_p: Ptr[float],
    qtrimecxa_p: Ptr[float],
):
    if scavenging_active != 0:
        qtcxa = qtt_l * fcxa
        qtdisrime_p[0] = _neu_wetdep_disgas_core(
            clwx * (fcxa / cfxx_l),
            fcxa,
            tcmass_n,
            hstar_ln,
            tem_l,
            pofl_l,
            qm_l,
            qtcxa,
        )
        qtdisstar = (qtdisrime_p[0] * qtcxa) / (qtdisrime_p[0] + qtcxa)
        qtrimecxa_p[0] = qtcxa * (
            1.0
            - exp(-0.24 * coleffrain * ((rca) ** 0.75) * (qtdisstar / qtcxa) * dtscav)
        )
        qtrimecxa_p[0] = min(
            qtrimecxa_p[0],
            ((rnew * garea * dtscav) / (clwx * qm_l * (fcxa / cfxx_l))) * qtdisstar,
        )
    else:
        qtdisrime_p[0] = 0.0
        qtrimecxa_p[0] = 0.0

@inline
def _neu_wetdep_impaction_washout(qt: float, rlocal: float, dtscav: float, coleffaer: float) -> float:
    if qt > 0.0:
        return qt * (1.0 - exp(-0.24 * coleffaer * ((rlocal) ** 0.75) * dtscav))

    return 0.0

@inline
def _neu_wetdep_gas_washout(
    rwash: float,
    boxf: float,
    dtscav: float,
    qtrtop: float,
    hstar: float,
    tm: float,
    pr: float,
    qm: float,
    qt: float,
    qtwash_p: Ptr[float],
    qtevap_p: Ptr[float],
):
    if qt > 0.0:
        _neu_wetdep_washgas_core(rwash, boxf, dtscav, qtrtop, hstar, tm, pr, qm, qt, qtwash_p, qtevap_p)
    else:
        qtwash_p[0] = 0.0
        qtevap_p[0] = 0.0

@inline
def _neu_wetdep_rnew_freezing_regime(
    licetyp: int,
    tem_l: float,
    tice: float,
    tfroz: float,
    rhosnowfix: float,
    dmin: float,
    volpow: float,
    rhorain: float,
    dempirical_impl: int,
    dca: float,
    rca: float,
    qtt_l: float,
    fcxa: float,
    fcxb: float,
    clwx: float,
    cfxx_l: float,
    tcmass_n: float,
    hstar_ln: float,
    pofl_l: float,
    qm_l: float,
    rnew: float,
    garea: float,
    dtscav: float,
    delz_l: float,
    rcxa_p: Ptr[float],
    rcxb_p: Ptr[float],
    dcxa_p: Ptr[float],
    dcxb_p: Ptr[float],
    qtraincxa_p: Ptr[float],
    qtraincxb_p: Ptr[float],
    qtrimecxa_p: Ptr[float],
    qtwashcxa_p: Ptr[float],
    qtevapcxa_p: Ptr[float],
):
    qtraincxa_p[0] = 0.0
    qtraincxb_p[0] = 0.0
    qtrimecxa_p[0] = 0.0
    qtwashcxa_p[0] = 0.0
    qtevapcxa_p[0] = 0.0
    dcxb_p[0] = 0.0
    dcxa_p[0] = 0.0
    rcxb_p[0] = 0.0
    rcxa_p[0] = 0.0

    deltarimemass = 0.0
    deltarime = 0.0
    dor = 0.0
    dnew = 0.0

    coleffsnow = exp(2.5e-2 * (tem_l - tice))
    if tem_l <= tfroz:
        rhosnow = rhosnowfix
    else:
        rhosnow = 0.303 * (tem_l - tfroz) * rhosnowfix

    if fcxa > 0.0:
        if dca > 0.0:
            deltarimemass = clwx * qm_l * (fcxa / cfxx_l) * (
                1.0
                - exp(
                    (-coleffsnow / (dca * 1.0e-3))
                    * ((rca) / (2.0 * rhosnow))
                    * dtscav
                )
            )
        else:
            deltarimemass = 0.0
    else:
        deltarimemass = 0.0

    if fcxa > 0.0:
        deltarime = min(rnew / fcxa, deltarimemass / (fcxa * garea * dtscav))
    else:
        deltarime = 0.0

    if rca > 0.0:
        dor = max(dmin, (((rca + deltarime) / rca) ** volpow) * dca)
    else:
        dor = 0.0

    rprecip = (rnew - (deltarime * fcxa)) / cfxx_l
    rcxa_p[0] = rca + deltarime + rprecip
    rcxb_p[0] = rprecip

    if rprecip > 0.0:
        wemp = (clwx * qm_l) / (garea * cfxx_l * delz_l)
        remp = rprecip / (rhorain / 1.0e3)
        dnew = _neu_wetdep_dempirical_eval(wemp, remp, dempirical_impl)
        dnew = max(dmin, dnew)
        if fcxb > 0.0:
            dcxb_p[0] = dnew
        else:
            dcxb_p[0] = 0.0
    else:
        dcxb_p[0] = 0.0

    if fcxa > 0.0:
        wemp = (clwx * qm_l * (fcxa / cfxx_l)) / (garea * fcxa * delz_l)
        remp = rcxa_p[0] / (rhorain / 1.0e3)
        demp = _neu_wetdep_dempirical_eval(wemp, remp, dempirical_impl)
        dcxa_p[0] = ((rca + deltarime) / rcxa_p[0]) * dor + (rprecip / rcxa_p[0]) * dnew
        dcxa_p[0] = max(demp, dcxa_p[0])
        dcxa_p[0] = max(dmin, dcxa_p[0])
    else:
        dcxa_p[0] = 0.0

    if qtt_l > 0.0:
        if rprecip > 0.0:
            _neu_wetdep_new_precip_scavenging(
                1 if licetyp == 1 else 0,
                rprecip,
                garea,
                dtscav,
                clwx,
                cfxx_l,
                tcmass_n,
                hstar_ln,
                tem_l,
                pofl_l,
                qm_l,
                qtt_l,
                fcxa,
                fcxb,
                qtraincxa_p,
                qtraincxb_p,
            )

        if deltarime > 0.0:
            _neu_wetdep_ice_riming_scavenging(
                1 if licetyp == 1 else 0,
                tem_l,
                tfroz,
                rhosnowfix,
                coleffsnow,
                dca,
                rca,
                qtt_l,
                fcxa,
                clwx,
                cfxx_l,
                tcmass_n,
                hstar_ln,
                pofl_l,
                qm_l,
                rnew,
                garea,
                dtscav,
                qtrimecxa_p,
            )
        else:
            qtrimecxa_p[0] = 0.0
    else:
        qtraincxa_p[0] = 0.0
        qtraincxb_p[0] = 0.0
        qtrimecxa_p[0] = 0.0

@inline
def _neu_wetdep_rnew_rain_regime(
    lwashtyp: int,
    coleffrain: float,
    coleffaer: float,
    four: float,
    rhorain: float,
    dca: float,
    rca: float,
    qtt_l: float,
    fca: float,
    fcxa: float,
    fcxb: float,
    qttopca: float,
    clwx: float,
    cfxx_l: float,
    tcmass_n: float,
    hstar_ln: float,
    tem_l: float,
    pofl_l: float,
    qm_l: float,
    rnew: float,
    garea: float,
    dtscav: float,
    delz_l: float,
    rcxa_p: Ptr[float],
    rcxb_p: Ptr[float],
    dcxa_p: Ptr[float],
    dcxb_p: Ptr[float],
    qtraincxa_p: Ptr[float],
    qtraincxb_p: Ptr[float],
    qtrimecxa_p: Ptr[float],
    qtwashcxa_p: Ptr[float],
    qtevapcxa_p: Ptr[float],
):
    qtraincxa_p[0] = 0.0
    qtraincxb_p[0] = 0.0
    qtrimecxa_p[0] = 0.0
    qtwashcxa_p[0] = 0.0
    qtevapcxa_p[0] = 0.0
    dcxb_p[0] = 0.0
    dcxa_p[0] = 0.0
    rcxb_p[0] = 0.0
    rcxa_p[0] = 0.0

    deltarimemass = 0.0
    deltarime = 0.0
    qtdisrime = 0.0

    if fcxa > 0.0:
        deltarimemass = (clwx * qm_l) * (fcxa / cfxx_l) * (
            1.0 - exp(-0.24 * coleffrain * ((rca) ** 0.75) * dtscav)
        )
    else:
        deltarimemass = 0.0

    if fcxa > 0.0:
        deltarime = min(rnew / fcxa, deltarimemass / (fcxa * garea * dtscav))
    else:
        deltarime = 0.0

    rprecip = (rnew - (deltarime * fcxa)) / cfxx_l
    rcxa_p[0] = rca + deltarime + rprecip
    rcxb_p[0] = rprecip
    dcxa_p[0] = four
    if fcxb > 0.0:
        dcxb_p[0] = four
    else:
        dcxb_p[0] = 0.0

    if qtt_l > 0.0:
        if rprecip > 0.0:
            _neu_wetdep_new_precip_scavenging(
                1,
                rprecip,
                garea,
                dtscav,
                clwx,
                cfxx_l,
                tcmass_n,
                hstar_ln,
                tem_l,
                pofl_l,
                qm_l,
                qtt_l,
                fcxa,
                fcxb,
                qtraincxa_p,
                qtraincxb_p,
            )

        if deltarime > 0.0:
            _neu_wetdep_rain_riming_scavenging(
                1,
                coleffrain,
                rca,
                qtt_l,
                fcxa,
                clwx,
                cfxx_l,
                tcmass_n,
                hstar_ln,
                tem_l,
                pofl_l,
                qm_l,
                rnew,
                garea,
                dtscav,
                __ptr__(qtdisrime),
                qtrimecxa_p,
            )
        else:
            qtrimecxa_p[0] = 0.0
    else:
        qtraincxa_p[0] = 0.0
        qtraincxb_p[0] = 0.0
        qtrimecxa_p[0] = 0.0

    if rca > 0.0:
        qtprecip = fcxa * qtt_l - qtdisrime
        if lwashtyp == 1:
            qtwashcxa_p[0] = _neu_wetdep_impaction_washout(qtprecip, rca, dtscav, coleffaer)
            qtevapcxa_p[0] = 0.0
        else:
            rwash = rca * garea
            wash_qtwash = 0.0
            wash_qtevap = 0.0
            _neu_wetdep_gas_washout(
                rwash,
                fca,
                dtscav,
                qttopca + qtrimecxa_p[0],
                hstar_ln,
                tem_l,
                pofl_l,
                qm_l,
                qtprecip,
                __ptr__(wash_qtwash),
                __ptr__(wash_qtevap),
            )
            qtwashcxa_p[0] = wash_qtwash
            qtevapcxa_p[0] = wash_qtevap

@inline
def _neu_wetdep_existing_precip_regime(
    clwc_l: float,
    ciwc_l: float,
    cfxx_l: float,
    fca: float,
    rls_l: float,
    garea: float,
    fax_in: float,
    rax_in: float,
    rca: float,
    dca: float,
    qttopaa: float,
    qttopca: float,
    qtevapaxp_in: float,
    freezing_l: int,
    licetyp: int,
    lwashtyp: int,
    tmix: float,
    volpow: float,
    four: float,
    coleffaer: float,
    dtscav: float,
    tcmass_n: float,
    hstar_ln: float,
    tem_l: float,
    pofl_l: float,
    qm_l: float,
    qtt_l: float,
    clwx_p: Ptr[float],
    fcxa_p: Ptr[float],
    fcxb_p: Ptr[float],
    rcxb_p: Ptr[float],
    dcxb_p: Ptr[float],
    qtraincxa_p: Ptr[float],
    qtraincxb_p: Ptr[float],
    qtrimecxa_p: Ptr[float],
    rcxa_p: Ptr[float],
    qtevapaxp_p: Ptr[float],
    fax_p: Ptr[float],
    rax_p: Ptr[float],
    qtevapcxa_p: Ptr[float],
    dcxa_p: Ptr[float],
    qtwashcxa_p: Ptr[float],
):
    clwx_p[0] = clwc_l + ciwc_l
    fcxa_p[0] = fca
    fcxb_p[0] = max(0.0, cfxx_l - fcxa_p[0])
    rcxb_p[0] = 0.0
    dcxb_p[0] = 0.0
    qtraincxa_p[0] = 0.0
    qtraincxb_p[0] = 0.0
    qtrimecxa_p[0] = 0.0
    qtwashcxa_p[0] = 0.0

    if fcxa_p[0] > 0.0:
        rcxa_p[0] = min(rca, rls_l / (garea * fcxa_p[0]))
        if fax_in > 0.0 and ((rcxa_p[0] + 1.0e-12) < rls_l / (garea * fcxa_p[0])):
            raxadjf = rls_l / garea - rcxa_p[0] * fcxa_p[0]
            rampct = raxadjf / (rax_in * fax_in)
            faxadj = rampct * fax_in
            if faxadj > 0.0:
                raxadj = raxadjf / faxadj
            else:
                raxadj = 0.0
        else:
            raxadj = 0.0
            rampct = 0.0
            faxadj = 0.0
    else:
        rcxa_p[0] = 0.0
        if fax_in > 0.0:
            raxadjf = rls_l / garea
            rampct = raxadjf / (rax_in * fax_in)
            faxadj = rampct * fax_in
            if faxadj > 0.0:
                raxadj = raxadjf / faxadj
            else:
                raxadj = 0.0
        else:
            raxadj = 0.0
            rampct = 0.0
            faxadj = 0.0

    qtevapaxp_p[0] = min(qttopaa, qttopaa - (rampct * (qttopaa - qtevapaxp_in)))
    fax_p[0] = faxadj
    rax_p[0] = raxadj

    if rcxa_p[0] <= 0.0:
        qtevapcxa_p[0] = qttopca
        rcxa_p[0] = 0.0
        dcxa_p[0] = 0.0
    else:
        if freezing_l != 0:
            dcxa_p[0] = ((rcxa_p[0] / rca) ** volpow) * dca
            if licetyp == 1:
                if tem_l <= tmix:
                    massloss = (rca - rcxa_p[0]) * fcxa_p[0] * garea * dtscav
                    qtevapcxa_p[0] = _neu_wetdep_disgas_core(
                        massloss / qm_l,
                        fcxa_p[0],
                        tcmass_n,
                        hstar_ln,
                        tem_l,
                        pofl_l,
                        qm_l,
                        qtt_l,
                    )
                    qtevapcxa_p[0] = min(qttopca, qtevapcxa_p[0])
                else:
                    qtevapcxa_p[0] = 0.0
            else:
                qtevapcxa_p[0] = 0.0
        else:
            qtevapcxap = (rca - rcxa_p[0]) / rca * qttopca
            dcxa_p[0] = four
            qtcxa = fcxa_p[0] * qtt_l
            qtdiscxa = 0.0
            if lwashtyp == 1:
                if qtt_l > 0.0:
                    qtdiscxa = _neu_wetdep_disgas_core(
                        clwx_p[0] * (fcxa_p[0] / cfxx_l),
                        fcxa_p[0],
                        tcmass_n,
                        hstar_ln,
                        tem_l,
                        pofl_l,
                        qm_l,
                        qtcxa,
                    )
                    qtwashcxa_p[0] = _neu_wetdep_impaction_washout(
                        qtcxa - qtdiscxa, rcxa_p[0], dtscav, coleffaer
                    )
                    qtevapcxaw = 0.0
                else:
                    qtwashcxa_p[0] = 0.0
                    qtevapcxaw = 0.0
            else:
                rwash = rcxa_p[0] * garea
                wash_qtwash = 0.0
                wash_qtevap = 0.0
                _neu_wetdep_gas_washout(
                    rwash,
                    fcxa_p[0],
                    dtscav,
                    qttopca,
                    hstar_ln,
                    tem_l,
                    pofl_l,
                    qm_l,
                    qtcxa - qtdiscxa,
                    __ptr__(wash_qtwash),
                    __ptr__(wash_qtevap),
                )
                qtwashcxa_p[0] = wash_qtwash
                qtevapcxaw = wash_qtevap
            qtevapcxa_p[0] = qtevapcxap + qtevapcxaw

@inline
def _neu_wetdep_ambient_washout_finalized(
    rax: float,
    freezing_l: int,
    fax: float,
    qtt_l: float,
    lwashtyp: int,
    dtscav: float,
    coleffaer: float,
    garea: float,
    qttopaa: float,
    hstar_ln: float,
    tem_l: float,
    pofl_l: float,
    qm_l: float,
    qtevapaxp: float,
    qtwashax_p: Ptr[float],
    qtevapaxw_p: Ptr[float],
    qtevapax_p: Ptr[float],
):
    if rax > 0.0:
        if freezing_l == 0:
            qtax = fax * qtt_l
            if lwashtyp == 1:
                qtwashax_p[0] = _neu_wetdep_impaction_washout(qtax, rax, dtscav, coleffaer)
                qtevapaxw_p[0] = 0.0
            else:
                rwash = rax * garea
                wash_qtwash = 0.0
                wash_qtevap = 0.0
                _neu_wetdep_gas_washout(
                    rwash,
                    fax,
                    dtscav,
                    qttopaa,
                    hstar_ln,
                    tem_l,
                    pofl_l,
                    qm_l,
                    qtax,
                    __ptr__(wash_qtwash),
                    __ptr__(wash_qtevap),
                )
                qtwashax_p[0] = wash_qtwash
                qtevapaxw_p[0] = wash_qtevap
        else:
            qtevapaxw_p[0] = 0.0
            qtwashax_p[0] = 0.0
    else:
        qtevapaxw_p[0] = 0.0
        qtwashax_p[0] = 0.0

    qtevapax_p[0] = qtevapaxp + qtevapaxw_p[0]

@inline
def _neu_wetdep_upper_level_redistribute(
    l: int,
    lm1: int,
    cfmin: float,
    adj_factor: float,
    garea: float,
    cfxx: Ptr[float],
    cfr: Ptr[float],
    rls: Ptr[float],
    evaprate: Ptr[float],
    fcxa: float,
    fcxb: float,
    fax: float,
    rcxa: float,
    rcxb: float,
    rax: float,
    dcxa: float,
    dcxb: float,
    dax: float,
    ampct_p: Ptr[float],
    amclpct_p: Ptr[float],
    clnewpct_p: Ptr[float],
    clnewampct_p: Ptr[float],
    cloldpct_p: Ptr[float],
    cloldampct_p: Ptr[float],
    fca_p: Ptr[float],
    rca_p: Ptr[float],
    dca_p: Ptr[float],
    fama_p: Ptr[float],
    rama_p: Ptr[float],
    dama_p: Ptr[float],
):
    fama_p[0] = max(fcxa + fcxb + fax - cfr[lm1 - 1], 0.0)

    if cfr[lm1 - 1] >= cfmin:
        cfxx[lm1 - 1] = cfr[lm1 - 1]
    else:
        if adj_factor * (rls[lm1 - 1] / garea) >= ((rcxa * fcxa + rcxb * fcxb + rax * fax) * (1.0 - evaprate[lm1 - 1])):
            cfxx[lm1 - 1] = cfmin
        else:
            cfxx[lm1 - 1] = cfr[lm1 - 1]

    if fax > 0.0:
        ampct_p[0] = max(0.0, min(1.0, (cfxx[l - 1] + fax - cfxx[lm1 - 1]) / fax))
        amclpct_p[0] = 1.0 - ampct_p[0]
    else:
        ampct_p[0] = 0.0
        amclpct_p[0] = 0.0

    if fcxb > 0.0:
        clnewpct_p[0] = max(0.0, min((cfxx[lm1 - 1] - fcxa) / fcxb, 1.0))
        clnewampct_p[0] = 1.0 - clnewpct_p[0]
    else:
        clnewpct_p[0] = 0.0
        clnewampct_p[0] = 0.0

    if fcxa > 0.0:
        cloldpct_p[0] = max(0.0, min(cfxx[lm1 - 1] / fcxa, 1.0))
        cloldampct_p[0] = 1.0 - cloldpct_p[0]
    else:
        cloldpct_p[0] = 0.0
        cloldampct_p[0] = 0.0

    fca_p[0] = min(cfxx[lm1 - 1], fcxa * cloldpct_p[0] + clnewpct_p[0] * fcxb + amclpct_p[0] * fax)
    if fca_p[0] > 0.0:
        rca_p[0] = (rcxa * fcxa * cloldpct_p[0] + rcxb * fcxb * clnewpct_p[0] + rax * fax * amclpct_p[0]) / fca_p[0]
        if rca_p[0] > 0.0:
            dca_p[0] = (rcxa * fcxa * cloldpct_p[0]) / (rca_p[0] * fca_p[0]) * dcxa + (
                rcxb * fcxb * clnewpct_p[0]
            ) / (rca_p[0] * fca_p[0]) * dcxb + (rax * fax * amclpct_p[0]) / (rca_p[0] * fca_p[0]) * dax
        else:
            dca_p[0] = 0.0
            fca_p[0] = 0.0
    else:
        fca_p[0] = 0.0
        dca_p[0] = 0.0
        rca_p[0] = 0.0

    fama_p[0] = fcxa + fcxb + fax - cfxx[lm1 - 1]
    if fama_p[0] > 0.0:
        rama_p[0] = (rcxa * fcxa * cloldampct_p[0] + rcxb * fcxb * clnewampct_p[0] + rax * fax * ampct_p[0]) / fama_p[0]
        if rama_p[0] > 0.0:
            dama_p[0] = (rcxa * fcxa * cloldampct_p[0]) / (rama_p[0] * fama_p[0]) * dcxa + (
                rcxb * fcxb * clnewampct_p[0]
            ) / (rama_p[0] * fama_p[0]) * dcxb + (rax * fax * ampct_p[0]) / (rama_p[0] * fama_p[0]) * dax
        else:
            fama_p[0] = 0.0
            dama_p[0] = 0.0
    else:
        fama_p[0] = 0.0
        dama_p[0] = 0.0
        rama_p[0] = 0.0

@inline
def _neu_wetdep_washo_level(
    l: int,
    lm1: int,
    do_diag: int,
    is_hno3: int,
    licetyp: int,
    lwashtyp: int,
    dempirical_impl: int,
    cfmin: float,
    cwmin: float,
    dmin: float,
    volpow: float,
    rhorain: float,
    rhosnowfix: float,
    coleffrain: float,
    tmix: float,
    tfroz: float,
    coleffaer: float,
    tice: float,
    four: float,
    dtscav: float,
    garea: float,
    adj_factor: float,
    qtt_l: float,
    qm_l: float,
    pofl_l: float,
    delz_l: float,
    rls_l: float,
    clwc_l: float,
    ciwc_l: float,
    tem_l: float,
    evaprate_l: float,
    hstar_ln: float,
    tcmass_n: float,
    cfxx: Ptr[float],
    cfr: Ptr[float],
    rls: Ptr[float],
    evaprate: Ptr[float],
    qt_rain: Ptr[float],
    qt_rime: Ptr[float],
    qt_wash: Ptr[float],
    qt_evap: Ptr[float],
    qttnew: Ptr[float],
    qttopaa_p: Ptr[float],
    qttopca_p: Ptr[float],
    rca_p: Ptr[float],
    fca_p: Ptr[float],
    dca_p: Ptr[float],
    rama_p: Ptr[float],
    fama_p: Ptr[float],
    dama_p: Ptr[float],
):
    fax = 0.0
    rax = 0.0
    dax = 0.0
    clwx = 0.0
    fcxa = 0.0
    fcxb = 0.0
    dcxa = 0.0
    dcxb = 0.0
    rcxa = 0.0
    rcxb = 0.0
    qtevapaxp = 0.0
    qtevapaxw = 0.0
    qtevapax = 0.0
    qtwashax = 0.0
    qtevapcxa = 0.0
    qtrimecxa = 0.0
    qtwashcxa = 0.0
    qtraincxa = 0.0
    qtraincxb = 0.0
    qttopaax = 0.0
    qttopcax = 0.0
    ampct = 0.0
    amclpct = 0.0
    clnewpct = 0.0
    clnewampct = 0.0
    cloldpct = 0.0
    cloldampct = 0.0

    freezing_l = 0
    if tem_l < tice:
        freezing_l = 1

    if rls_l > 0.0:
        fax = max(0.0, fama_p[0] * (1.0 - evaprate_l))
        rax = rama_p[0]
        if fama_p[0] > 0.0:
            if freezing_l != 0:
                dax = dama_p[0]
            else:
                dax = four
        else:
            dax = 0.0

        if rama_p[0] > 0.0:
            qtevapaxp = min(qttopaa_p[0], evaprate_l * qttopaa_p[0])
        else:
            qtevapaxp = 0.0

        wrk = rax * fax + rca_p[0] * fca_p[0]
        if wrk > 0.0:
            rnew_tst = rls_l / (garea * wrk)
        else:
            rnew_tst = 10.0
        rnew = (rls_l / garea) - (rax * fax + rca_p[0] * fca_p[0])

        if (rls_l / garea) > adj_factor * (rax * fax + rca_p[0] * fca_p[0]):
            if cfxx[l - 1] == 0.0:
                return 1

            clwx = max(clwc_l + ciwc_l, cwmin * cfxx[l - 1])
            fcxa = fca_p[0]
            fcxb = max(0.0, cfxx[l - 1] - fcxa)

            if freezing_l != 0:
                _neu_wetdep_rnew_freezing_regime(
                    licetyp,
                    tem_l,
                    tice,
                    tfroz,
                    rhosnowfix,
                    dmin,
                    volpow,
                    rhorain,
                    dempirical_impl,
                    dca_p[0],
                    rca_p[0],
                    qtt_l,
                    fcxa,
                    fcxb,
                    clwx,
                    cfxx[l - 1],
                    tcmass_n,
                    hstar_ln,
                    pofl_l,
                    qm_l,
                    rnew,
                    garea,
                    dtscav,
                    delz_l,
                    __ptr__(rcxa),
                    __ptr__(rcxb),
                    __ptr__(dcxa),
                    __ptr__(dcxb),
                    __ptr__(qtraincxa),
                    __ptr__(qtraincxb),
                    __ptr__(qtrimecxa),
                    __ptr__(qtwashcxa),
                    __ptr__(qtevapcxa),
                )
            else:
                _neu_wetdep_rnew_rain_regime(
                    lwashtyp,
                    coleffrain,
                    coleffaer,
                    four,
                    rhorain,
                    dca_p[0],
                    rca_p[0],
                    qtt_l,
                    fca_p[0],
                    fcxa,
                    fcxb,
                    qttopca_p[0],
                    clwx,
                    cfxx[l - 1],
                    tcmass_n,
                    hstar_ln,
                    tem_l,
                    pofl_l,
                    qm_l,
                    rnew,
                    garea,
                    dtscav,
                    delz_l,
                    __ptr__(rcxa),
                    __ptr__(rcxb),
                    __ptr__(dcxa),
                    __ptr__(dcxb),
                    __ptr__(qtraincxa),
                    __ptr__(qtraincxb),
                    __ptr__(qtrimecxa),
                    __ptr__(qtwashcxa),
                    __ptr__(qtevapcxa),
                )
        else:
            _neu_wetdep_existing_precip_regime(
                clwc_l,
                ciwc_l,
                cfxx[l - 1],
                fca_p[0],
                rls_l,
                garea,
                fax,
                rax,
                rca_p[0],
                dca_p[0],
                qttopaa_p[0],
                qttopca_p[0],
                qtevapaxp,
                freezing_l,
                licetyp,
                lwashtyp,
                tmix,
                volpow,
                four,
                coleffaer,
                dtscav,
                tcmass_n,
                hstar_ln,
                tem_l,
                pofl_l,
                qm_l,
                qtt_l,
                __ptr__(clwx),
                __ptr__(fcxa),
                __ptr__(fcxb),
                __ptr__(rcxb),
                __ptr__(dcxb),
                __ptr__(qtraincxa),
                __ptr__(qtraincxb),
                __ptr__(qtrimecxa),
                __ptr__(rcxa),
                __ptr__(qtevapaxp),
                __ptr__(fax),
                __ptr__(rax),
                __ptr__(qtevapcxa),
                __ptr__(dcxa),
                __ptr__(qtwashcxa),
            )
    else:
        qtevapcxa = qttopca_p[0]
        qtevapax = qttopaa_p[0]
        if l > 1:
            if rls[lm1 - 1] > 0.0:
                cfxx[lm1 - 1] = max(cfmin, cfr[lm1 - 1])
            else:
                cfxx[lm1 - 1] = cfr[lm1 - 1]
        rca_p[0] = 0.0
        rama_p[0] = 0.0
        fca_p[0] = 0.0
        fama_p[0] = 0.0
        dca_p[0] = 0.0
        dama_p[0] = 0.0

    if rls_l > 0.0:
        _neu_wetdep_ambient_washout_finalized(
            rax,
            freezing_l,
            fax,
            qtt_l,
            lwashtyp,
            dtscav,
            coleffaer,
            garea,
            qttopaa_p[0],
            hstar_ln,
            tem_l,
            pofl_l,
            qm_l,
            qtevapaxp,
            __ptr__(qtwashax),
            __ptr__(qtevapaxw),
            __ptr__(qtevapax),
        )

        if l > 1:
            _neu_wetdep_upper_level_redistribute(
                l,
                lm1,
                cfmin,
                adj_factor,
                garea,
                cfxx,
                cfr,
                rls,
                evaprate,
                fcxa,
                fcxb,
                fax,
                rcxa,
                rcxb,
                rax,
                dcxa,
                dcxb,
                dax,
                __ptr__(ampct),
                __ptr__(amclpct),
                __ptr__(clnewpct),
                __ptr__(clnewampct),
                __ptr__(cloldpct),
                __ptr__(cloldampct),
                fca_p,
                rca_p,
                dca_p,
                fama_p,
                rama_p,
                dama_p,
            )
        else:
            ampct = 0.0
            amclpct = 0.0
            clnewpct = 0.0
            clnewampct = 0.0
            cloldpct = 0.0
            cloldampct = 0.0

    qtnetlcxa = qtraincxa + qtrimecxa + qtwashcxa - qtevapcxa
    qtnetlcxa = min(qtt_l * fcxa, qtnetlcxa)
    qtnetlcxb = qtraincxb
    qtnetlcxb = min(qtt_l * fcxb, qtnetlcxb)
    qtnetlax = qtwashax - qtevapax
    qtnetlax = min(qtt_l * fax, qtnetlax)
    qttnew[l - 1] = qtt_l - (qtnetlcxa + qtnetlcxb + qtnetlax)

    if do_diag != 0 and is_hno3 != 0:
        qt_rain[l - 1] = qtraincxa + qtraincxb
        qt_rime[l - 1] = qtrimecxa
        qt_wash[l - 1] = qtwashcxa + qtwashax
        qt_evap[l - 1] = qtevapcxa + qtevapax

    qttopcax = (qttopca_p[0] + qtnetlcxa) * cloldpct + qtnetlcxb * clnewpct + (qttopaa_p[0] + qtnetlax) * amclpct
    qttopaax = (qttopca_p[0] + qtnetlcxa) * cloldampct + qtnetlcxb * clnewampct + (qttopaa_p[0] + qtnetlax) * ampct
    qttopca_p[0] = qttopcax
    qttopaa_p[0] = qttopaax
    return 0

@inline
def _neu_wetdep_washo_species(
    n: int,
    le: int,
    lpar: int,
    hno3_ndx: int,
    do_diag: int,
    dempirical_impl: int,
    dtscav: float,
    garea: float,
    adj_factor: float,
    cfmin: float,
    cwmin: float,
    dmin: float,
    volpow: float,
    rhorain: float,
    rhosnowfix: float,
    coleffrain: float,
    tmix: float,
    tfroz: float,
    coleffaer: float,
    tice: float,
    four: float,
    qttjfl: Ptr[float],
    qm: Ptr[float],
    pofl: Ptr[float],
    delz: Ptr[float],
    rls: Ptr[float],
    clwc: Ptr[float],
    ciwc: Ptr[float],
    cfr: Ptr[float],
    tem: Ptr[float],
    evaprate: Ptr[float],
    hstar: Ptr[float],
    tcmass: Ptr[float],
    tckaqb: Ptr[int],
    tcnion: Ptr[int],
    qt_rain: Ptr[float],
    qt_rime: Ptr[float],
    qt_wash: Ptr[float],
    qt_evap: Ptr[float],
    cfxx: Ptr[float],
    qtt: Ptr[float],
    qttnew: Ptr[float],
):
    ll = 1
    while ll <= lpar:
        ln_idx = _idx2(ll, n, lpar)
        qtt[ll - 1] = qttjfl[ln_idx]
        qttnew[ll - 1] = qttjfl[ln_idx]
        ll += 1

    is_hno3 = 0
    if n == hno3_ndx:
        is_hno3 = 1
        ll = 1
        while ll <= lpar:
            qt_rain[ll - 1] = 0.0
            qt_rime[ll - 1] = 0.0
            qt_wash[ll - 1] = 0.0
            qt_evap[ll - 1] = 0.0
            ll += 1

    if tckaqb[n - 1] != 0:
        lwashtyp = 1
    else:
        lwashtyp = 2

    if tcnion[n - 1] != 0:
        licetyp = 1
    else:
        licetyp = 2

    qttopaa = 0.0
    qttopca = 0.0
    rca = 0.0
    fca = 0.0
    dca = 0.0
    rama = 0.0
    fama = 0.0
    dama = 0.0

    if le >= 1:
        if rls[le - 1] > 0.0:
            cfxx[le - 1] = max(cfmin, cfr[le - 1])
        else:
            cfxx[le - 1] = cfr[le - 1]

    l = le
    while l >= 1:
        lm1 = l - 1
        ln_idx = _idx2(l, n, lpar)
        hstar_ln = hstar[ln_idx]

        if (
            _neu_wetdep_washo_level(
                l,
                lm1,
                do_diag,
                is_hno3,
                licetyp,
                lwashtyp,
                dempirical_impl,
                cfmin,
                cwmin,
                dmin,
                volpow,
                rhorain,
                rhosnowfix,
                coleffrain,
                tmix,
                tfroz,
                coleffaer,
                tice,
                four,
                dtscav,
                garea,
                adj_factor,
                qtt[l - 1],
                qm[l - 1],
                pofl[l - 1],
                delz[l - 1],
                rls[l - 1],
                clwc[l - 1],
                ciwc[l - 1],
                tem[l - 1],
                evaprate[l - 1],
                hstar_ln,
                tcmass[n - 1],
                cfxx,
                cfr,
                rls,
                evaprate,
                qt_rain,
                qt_rime,
                qt_wash,
                qt_evap,
                qttnew,
                __ptr__(qttopaa),
                __ptr__(qttopca),
                __ptr__(rca),
                __ptr__(fca),
                __ptr__(dca),
                __ptr__(rama),
                __ptr__(fama),
                __ptr__(dama),
            )
            != 0
        ):
            ll = 1
            while ll <= lpar:
                qttjfl[_idx2(ll, n, lpar)] = qtt[ll - 1]
                ll += 1
            return

        l -= 1

    ll = 1
    while ll <= le:
        qttjfl[_idx2(ll, n, lpar)] = qttnew[ll - 1]
        ll += 1

def neu_wetdep_disgas_codon(
    clwx: float,
    cfx: float,
    molmass: float,
    hstar: float,
    tm: float,
    pr: float,
    qm: float,
    qt: float,
    qtdis_p: cobj,
):
    qtdis = Ptr[float](qtdis_p)
    qtdis[0] = _neu_wetdep_disgas_core(clwx, cfx, molmass, hstar, tm, pr, qm, qt)

def neu_wetdep_raingas_codon(
    rrain: float,
    dtscav: float,
    clwx: float,
    cfx: float,
    qm: float,
    qt: float,
    qtdis: float,
    qtrain_p: cobj,
):
    qtrain = Ptr[float](qtrain_p)
    qtrain[0] = _neu_wetdep_raingas_core(rrain, dtscav, clwx, cfx, qm, qt, qtdis)

def neu_wetdep_washgas_codon(
    rwash: float,
    boxf: float,
    dtscav: float,
    qtrtop: float,
    hstar: float,
    tm: float,
    pr: float,
    qm: float,
    qt: float,
    qtwash_p: cobj,
    qtevap_p: cobj,
):
    qtwash = Ptr[float](qtwash_p)
    qtevap = Ptr[float](qtevap_p)
    _neu_wetdep_washgas_core(rwash, boxf, dtscav, qtrtop, hstar, tm, pr, qm, qt, qtwash, qtevap)

def neu_wetdep_washo_codon(
    lpar: int,
    ntrace: int,
    hno3_ndx: int,
    do_diag: int,
    dempirical_impl: int,
    dtscav: float,
    garea: float,
    adj_factor: float,
    qttjfl_p: cobj,
    qm_p: cobj,
    pofl_p: cobj,
    delz_p: cobj,
    rls_p: cobj,
    clwc_p: cobj,
    ciwc_p: cobj,
    cfr_p: cobj,
    tem_p: cobj,
    evaprate_p: cobj,
    hstar_p: cobj,
    tcmass_p: cobj,
    tckaqb_p: cobj,
    tcnion_p: cobj,
    qt_rain_p: cobj,
    qt_rime_p: cobj,
    qt_wash_p: cobj,
    qt_evap_p: cobj,
    cfxx_p: cobj,
    qtt_p: cobj,
    qttnew_p: cobj,
):
    qttjfl = Ptr[float](qttjfl_p)
    qm = Ptr[float](qm_p)
    pofl = Ptr[float](pofl_p)
    delz = Ptr[float](delz_p)
    rls = Ptr[float](rls_p)
    clwc = Ptr[float](clwc_p)
    ciwc = Ptr[float](ciwc_p)
    cfr = Ptr[float](cfr_p)
    tem = Ptr[float](tem_p)
    evaprate = Ptr[float](evaprate_p)
    hstar = Ptr[float](hstar_p)
    tcmass = Ptr[float](tcmass_p)
    tckaqb = Ptr[int](tckaqb_p)
    tcnion = Ptr[int](tcnion_p)
    qt_rain = Ptr[float](qt_rain_p)
    qt_rime = Ptr[float](qt_rime_p)
    qt_wash = Ptr[float](qt_wash_p)
    qt_evap = Ptr[float](qt_evap_p)
    cfxx = Ptr[float](cfxx_p)
    qtt = Ptr[float](qtt_p)
    qttnew = Ptr[float](qttnew_p)

    zero = 0.0
    one = 1.0
    cfmin = 0.1
    cwmin = 1.0e-5
    dmin = 1.0e-1
    volpow = 1.0 / 3.0
    rhorain = 1.0e3
    rhosnowfix = 1.0e2
    coleffrain = 0.7
    tmix = 258.0
    tfroz = 240.0
    coleffaer = 0.05
    tice = 263.0
    four = 4.0

    le = lpar - 1
    n = 1
    while n <= ntrace:
        _neu_wetdep_washo_species(
            n,
            le,
            lpar,
            hno3_ndx,
            do_diag,
            dempirical_impl,
            dtscav,
            garea,
            adj_factor,
            cfmin,
            cwmin,
            dmin,
            volpow,
            rhorain,
            rhosnowfix,
            coleffrain,
            tmix,
            tfroz,
            coleffaer,
            tice,
            four,
            qttjfl,
            qm,
            pofl,
            delz,
            rls,
            clwc,
            ciwc,
            cfr,
            tem,
            evaprate,
            hstar,
            tcmass,
            tckaqb,
            tcnion,
            qt_rain,
            qt_rime,
            qt_wash,
            qt_evap,
            cfxx,
            qtt,
            qttnew,
        )

        n += 1

def setsox_init_fields_codon(
    stage: int,
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    nfs: int,
    cloud_borne_flag: int,
    inv_so2_flag: int,
    inv_h2o2_flag: int,
    inv_o3_flag: int,
    inv_ho2_flag: int,
    id_so2: int,
    id_hno3: int,
    id_h2o2: int,
    id_nh3: int,
    id_o3: int,
    id_ho2: int,
    id_h2so4: int,
    id_so4: int,
    id_msa: int,
    ph0: float,
    xhnm_p: cobj,
    invariants_p: cobj,
    qin_p: cobj,
    cfact_p: cobj,
    xph_p: cobj,
    xso2_p: cobj,
    xhno3_p: cobj,
    xh2o2_p: cobj,
    xnh3_p: cobj,
    xo3_p: cobj,
    xho2_p: cobj,
    xh2so4_p: cobj,
    xso4_p: cobj,
    xno3_p: cobj,
    xnh4_p: cobj,
    xmsa_p: cobj,
):
    xhnm = Ptr[float](xhnm_p)
    cfact = Ptr[float](cfact_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, ncol)
            cfact[idx] = xhnm[idx] * 1.0e6 * 1.38e-23 / 287.0 * 1.0e-3

    if stage == 1:
        return

    invariants = Ptr[float](invariants_p)
    qin = Ptr[float](qin_p)
    xph = Ptr[float](xph_p)
    xso2 = Ptr[float](xso2_p)
    xhno3 = Ptr[float](xhno3_p)
    xh2o2 = Ptr[float](xh2o2_p)
    xnh3 = Ptr[float](xnh3_p)
    xo3 = Ptr[float](xo3_p)
    xho2 = Ptr[float](xho2_p)
    xh2so4 = Ptr[float](xh2so4_p)
    xso4 = Ptr[float](xso4_p)
    xno3 = Ptr[float](xno3_p)
    xnh4 = Ptr[float](xnh4_p)
    xmsa = Ptr[float](xmsa_p)
    xph0 = 10.0 ** (-ph0)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, ncol)
            xso4[idx] = 0.0
            xno3[idx] = 0.0
            xnh4[idx] = 0.0
            xph[idx] = xph0

            if inv_so2_flag != 0:
                xso2[idx] = invariants[_idx3(i, k, id_so2, ncol, pver)] / xhnm[idx]
            else:
                xso2[idx] = qin[_idx3(i, k, id_so2, ncol, pver)]

            if id_hno3 > 0:
                xhno3[idx] = qin[_idx3(i, k, id_hno3, ncol, pver)]
            else:
                xhno3[idx] = 0.0

            if inv_h2o2_flag != 0:
                xh2o2[idx] = invariants[_idx3(i, k, id_h2o2, ncol, pver)] / xhnm[idx]
            else:
                xh2o2[idx] = qin[_idx3(i, k, id_h2o2, ncol, pver)]

            if id_nh3 > 0:
                xnh3[idx] = qin[_idx3(i, k, id_nh3, ncol, pver)]
            else:
                xnh3[idx] = 0.0

            if inv_o3_flag != 0:
                xo3[idx] = invariants[_idx3(i, k, id_o3, ncol, pver)] / xhnm[idx]
            else:
                xo3[idx] = qin[_idx3(i, k, id_o3, ncol, pver)]

            if inv_ho2_flag != 0:
                xho2[idx] = invariants[_idx3(i, k, id_ho2, ncol, pver)] / xhnm[idx]
            else:
                xho2[idx] = qin[_idx3(i, k, id_ho2, ncol, pver)]

            if cloud_borne_flag != 0:
                xh2so4[idx] = qin[_idx3(i, k, id_h2so4, ncol, pver)]
            else:
                xso4[idx] = qin[_idx3(i, k, id_so4, ncol, pver)]

            if id_msa > 0:
                xmsa[idx] = qin[_idx3(i, k, id_msa, ncol, pver)]

def setsox_ph_solve_codon(
    ncol: int,
    pcols: int,
    pver: int,
    itermax: int,
    cloud_borne_flag: int,
    const0: float,
    ra: float,
    xkw: float,
    so4_fact: float,
    press_p: cobj,
    tfld_p: cobj,
    cldfrc_p: cobj,
    xhnm_p: cobj,
    xlwc_p: cobj,
    xso4c_p: cobj,
    xnh4c_p: cobj,
    xno3c_p: cobj,
    xso4_p: cobj,
    xnh4_p: cobj,
    xno3_p: cobj,
    xso2_p: cobj,
    xhno3_p: cobj,
    xnh3_p: cobj,
    xph_p: cobj,
):
    press = Ptr[float](press_p)
    tfld = Ptr[float](tfld_p)
    cldfrc = Ptr[float](cldfrc_p)
    xhnm = Ptr[float](xhnm_p)
    xlwc = Ptr[float](xlwc_p)
    xso4c = Ptr[float](xso4c_p)
    xnh4c = Ptr[float](xnh4c_p)
    xno3c = Ptr[float](xno3c_p)
    xso4 = Ptr[float](xso4_p)
    xnh4 = Ptr[float](xnh4_p)
    xno3 = Ptr[float](xno3_p)
    xso2 = Ptr[float](xso2_p)
    xhno3 = Ptr[float](xhno3_p)
    xnh3 = Ptr[float](xnh3_p)
    xph = Ptr[float](xph_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, ncol)
            idxp = _idx2(i, k, pcols)
            if cloud_borne_flag != 0 and cldfrc[idxp] > 0.0:
                xso4[idx] = xso4c[idxp] / cldfrc[idxp]
                xnh4[idx] = xnh4c[idxp] / cldfrc[idxp]
                xno3[idx] = xno3c[idxp] / cldfrc[idxp]

            xl = xlwc[idxp]
            if xl >= 1.0e-8:
                work1 = 1.0 / tfld[idxp] - 1.0 / 298.0
                pz = 0.01 * press[idxp]
                tz = tfld[idxp]
                patm = pz / 1013.0

                xk = 2.1e5 * exp(8700.0 * work1)
                xe = 15.4
                fact1_hno3 = xk * xe * patm * xhno3[idx]
                fact2_hno3 = xk * ra * tz * xl
                fact3_hno3 = xe

                xk = 1.23 * exp(3120.0 * work1)
                xe = 1.7e-2 * exp(2090.0 * work1)
                x2 = 6.0e-8 * exp(1120.0 * work1)
                fact1_so2 = xk * xe * patm * xso2[idx]
                fact2_so2 = xk * ra * tz * xl
                fact3_so2 = xe
                fact4_so2 = x2

                xk = 58.0 * exp(4085.0 * work1)
                xe = 1.7e-5 * exp(-4325.0 * work1)
                fact1_nh3 = (xk * xe * patm / xkw) * (xnh3[idx] + xnh4[idx])
                fact2_nh3 = xk * ra * tz * xl
                fact3_nh3 = xe / xkw

                eh2o = xkw
                co2g = 330.0e-6
                xk = 3.1e-2 * exp(2423.0 * work1)
                xe = 4.3e-7 * exp(-913.0 * work1)
                eco2 = xk * xe * co2g * patm
                eso4 = xso4[idx] * xhnm[idx] * const0 / xl

                converged = 0
                yph_lo = 0.0
                yph_hi = 0.0
                ynetpos_lo = 0.0
                ynetpos_hi = 0.0
                for iter in range(1, itermax + 1):
                    if iter == 1:
                        yph_lo = 2.0
                        yph_hi = yph_lo
                        yph = yph_lo
                    elif iter == 2:
                        yph_hi = 7.0
                        yph = yph_hi
                    else:
                        yph = 0.5 * (yph_lo + yph_hi)

                    xph[idx] = 10.0 ** (-yph)
                    ehno3 = fact1_hno3 / (1.0 + fact2_hno3 * (1.0 + fact3_hno3 / xph[idx]))
                    eso2 = fact1_so2 / (
                        1.0
                        + fact2_so2
                        * (1.0 + (fact3_so2 / xph[idx]) * (1.0 + fact4_so2 / xph[idx]))
                    )
                    enh3 = fact1_nh3 / (1.0 + fact2_nh3 * (1.0 + fact3_nh3 * xph[idx]))

                    tmp_nh4 = enh3 * xph[idx]
                    tmp_hso3 = eso2 / xph[idx]
                    tmp_so3 = tmp_hso3 * 2.0 * fact4_so2 / xph[idx]
                    tmp_hco3 = eco2 / xph[idx]
                    tmp_oh = eh2o / xph[idx]
                    tmp_no3 = ehno3 / xph[idx]
                    tmp_so4 = so4_fact * eso4
                    tmp_pos = xph[idx] + tmp_nh4
                    tmp_neg = tmp_oh + tmp_hco3 + tmp_no3 + tmp_hso3 + tmp_so3 + tmp_so4
                    ynetpos = tmp_pos - tmp_neg

                    if iter > 2:
                        if ynetpos == 0.0:
                            converged = 1
                            break
                        elif ynetpos >= 0.0:
                            yph_lo = yph
                            ynetpos_lo = ynetpos
                        else:
                            yph_hi = yph
                            ynetpos_hi = ynetpos

                        if abs(yph_hi - yph_lo) <= 0.005:
                            yph = 0.5 * (yph_hi + yph_lo)
                            xph[idx] = 10.0 ** (-yph)
                            converged = 1
                            break
                    elif iter == 1:
                        if ynetpos <= 0.0:
                            converged = 1
                            break
                        ynetpos_lo = ynetpos
                    else:
                        if ynetpos >= 0.0:
                            converged = 1
                            break
                        ynetpos_hi = ynetpos
            else:
                xph[idx] = 1.0e-7

def setsox_aqchem_predict_codon(
    ncol: int,
    pcols: int,
    pver: int,
    cloud_borne_flag: int,
    modal_aerosols_flag: int,
    id_nh3: int,
    dtime: float,
    const0: float,
    kh0: float,
    kh1: float,
    kh2: float,
    kh3: float,
    ra: float,
    xkw: float,
    press_p: cobj,
    tfld_p: cobj,
    xhnm_p: cobj,
    xlwc_p: cobj,
    xph_p: cobj,
    xho2_p: cobj,
    xhno3_p: cobj,
    xno3_p: cobj,
    xh2o2_p: cobj,
    xso2_p: cobj,
    xo3_p: cobj,
    xnh3_p: cobj,
    xnh4_p: cobj,
    xso4_p: cobj,
    xso4_init_p: cobj,
    xdelso4hp_p: cobj,
    hno3g_p: cobj,
    nh3g_p: cobj,
    hehno3_p: cobj,
    heh2o2_p: cobj,
    heso2_p: cobj,
    henh3_p: cobj,
    heo3_p: cobj,
):
    press = Ptr[float](press_p)
    tfld = Ptr[float](tfld_p)
    xhnm = Ptr[float](xhnm_p)
    xlwc = Ptr[float](xlwc_p)
    xph = Ptr[float](xph_p)
    xho2 = Ptr[float](xho2_p)
    xhno3 = Ptr[float](xhno3_p)
    xno3 = Ptr[float](xno3_p)
    xh2o2 = Ptr[float](xh2o2_p)
    xso2 = Ptr[float](xso2_p)
    xo3 = Ptr[float](xo3_p)
    xnh3 = Ptr[float](xnh3_p)
    xnh4 = Ptr[float](xnh4_p)
    xso4 = Ptr[float](xso4_p)
    xso4_init = Ptr[float](xso4_init_p)
    xdelso4hp = Ptr[float](xdelso4hp_p)
    hno3g = Ptr[float](hno3g_p)
    nh3g = Ptr[float](nh3g_p)
    hehno3 = Ptr[float](hehno3_p)
    heh2o2 = Ptr[float](heh2o2_p)
    heso2 = Ptr[float](heso2_p)
    henh3 = Ptr[float](henh3_p)
    heo3 = Ptr[float](heo3_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, ncol)
            idxp = _idx2(i, k, pcols)
            work1 = 1.0 / tfld[idxp] - 1.0 / 298.0
            tz = tfld[idxp]
            xl = xlwc[idxp]
            patm = press[idxp] / 101300.0
            xam = press[idxp] / (1.38e-23 * tz)

            xk = 2.1e5 * exp(8700.0 * work1)
            xe = 15.4
            hehno3[idx] = xk * (1.0 + xe / xph[idx])

            xk = 7.4e4 * exp(6621.0 * work1)
            xe = 2.2e-12 * exp(-3730.0 * work1)
            heh2o2[idx] = xk * (1.0 + xe / xph[idx])

            xk = 1.23 * exp(3120.0 * work1)
            xe = 1.7e-2 * exp(2090.0 * work1)
            x2 = 6.0e-8 * exp(1120.0 * work1)
            wrk = xe / xph[idx]
            heso2[idx] = xk * (1.0 + wrk * (1.0 + x2 / xph[idx]))

            xk = 58.0 * exp(4085.0 * work1)
            xe = 1.7e-5 * exp(-4325.0 * work1)
            henh3[idx] = xk * (1.0 + xe * xph[idx] / xkw)

            xk = 1.15e-2 * exp(2560.0 * work1)
            heo3[idx] = xk

            kh4 = (kh2 + kh3 * kh1 / xph[idx]) / ((1.0 + kh1 / xph[idx]) ** 2)
            ho2s = kh0 * xho2[idx] * patm * (1.0 + kh1 / xph[idx])
            r1h2o2 = kh4 * ho2s * ho2s

            if cloud_borne_flag != 0:
                r2h2o2 = r1h2o2 * xl / const0 * 1.0e6 / xam
            else:
                r2h2o2 = r1h2o2 * xl * const0 / xam

            if modal_aerosols_flag == 0:
                xh2o2[idx] = xh2o2[idx] + r2h2o2 * dtime

            px = hehno3[idx] * ra * tz * xl
            hno3g[idx] = (xhno3[idx] + xno3[idx]) / (1.0 + px)

            px = heh2o2[idx] * ra * tz * xl
            h2o2g = xh2o2[idx] / (1.0 + px)

            px = heso2[idx] * ra * tz * xl
            so2g = xso2[idx] / (1.0 + px)

            px = heo3[idx] * ra * tz * xl
            o3g = xo3[idx] / (1.0 + px)

            px = henh3[idx] * ra * tz * xl
            if id_nh3 > 0:
                nh3g[idx] = (xnh3[idx] + xnh4[idx]) / (1.0 + px)
            else:
                nh3g[idx] = 0.0

            rah2o2 = 8.0e4 * exp(-3650.0 * work1) / (0.1 + xph[idx])
            rao3 = 4.39e11 * exp(-4131.0 / tz) + 2.56e3 * exp(-996.0 / tz) / xph[idx]

            if xl >= 1.0e-8:
                if cloud_borne_flag != 0:
                    patm_x = patm
                else:
                    patm_x = 1.0

                if modal_aerosols_flag != 0:
                    pso4 = (
                        rah2o2
                        * 7.4e4
                        * exp(6621.0 * work1)
                        * h2o2g
                        * patm_x
                        * 1.23
                        * exp(3120.0 * work1)
                        * so2g
                        * patm_x
                    )
                else:
                    pso4 = rah2o2 * heh2o2[idx] * h2o2g * patm_x * heso2[idx] * so2g * patm_x

                pso4 = pso4 * xl / const0 / xhnm[idx]
                ccc = pso4 * dtime
                ccc = max(ccc, 1.0e-30)
                xso4_init[idx] = xso4[idx]

                if xh2o2[idx] > xso2[idx]:
                    if ccc > xso2[idx]:
                        xso4[idx] = xso4[idx] + xso2[idx]
                        if cloud_borne_flag != 0:
                            xh2o2[idx] = xh2o2[idx] - xso2[idx]
                            xso2[idx] = 1.0e-20
                        else:
                            xso2[idx] = 1.0e-20
                            xh2o2[idx] = xh2o2[idx] - xso2[idx]
                    else:
                        xso4[idx] = xso4[idx] + ccc
                        xh2o2[idx] = xh2o2[idx] - ccc
                        xso2[idx] = xso2[idx] - ccc
                else:
                    if ccc > xh2o2[idx]:
                        xso4[idx] = xso4[idx] + xh2o2[idx]
                        xso2[idx] = xso2[idx] - xh2o2[idx]
                        xh2o2[idx] = 1.0e-20
                    else:
                        xso4[idx] = xso4[idx] + ccc
                        xh2o2[idx] = xh2o2[idx] - ccc
                        xso2[idx] = xso2[idx] - ccc

                if modal_aerosols_flag != 0:
                    xdelso4hp[idx] = xso4[idx] - xso4_init[idx]

                pso4 = rao3 * heo3[idx] * o3g * patm_x * heso2[idx] * so2g * patm_x
                pso4 = pso4 * xl / const0 / xhnm[idx]
                ccc = pso4 * dtime
                ccc = max(ccc, 1.0e-30)
                xso4_init[idx] = xso4[idx]

                if ccc > xso2[idx]:
                    xso4[idx] = xso4[idx] + xso2[idx]
                    xso2[idx] = 1.0e-20
                else:
                    xso4[idx] = xso4[idx] + ccc
                    xso2[idx] = xso2[idx] - ccc

def setsox_xph_lwc_diag_codon(
    ncol: int,
    pcols: int,
    pver: int,
    cldfrc_p: cobj,
    lwc_p: cobj,
    xph_p: cobj,
    xphlwc_p: cobj,
):
    cldfrc = Ptr[float](cldfrc_p)
    lwc = Ptr[float](lwc_p)
    xph = Ptr[float](xph_p)
    xphlwc = Ptr[float](xphlwc_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, ncol)
            idxp = _idx2(i, k, pcols)
            xphlwc[idx] = 0.0
            if cldfrc[idxp] >= 1.0e-5 and lwc[idx] >= 1.0e-8:
                xphlwc[idx] = -1.0 * log10(xph[idx]) * lwc[idx]

@inline
def _sox_cldaero_uptakerate(
    xl: float,
    cldnum: float,
    cfact: float,
    cldfrc: float,
    tfld: float,
    press: float,
    pi_val: float,
) -> float:
    num_cd = 1.0e-3 * cldnum * cfact / cldfrc
    num_cd = max(num_cd, 0.0)
    volx34pi_cd = xl * 0.75 / pi_val
    radxnum_cd = (volx34pi_cd * num_cd * num_cd) ** 0.3333333
    if radxnum_cd <= volx34pi_cd * 4.0e4:
        radxnum_cd = volx34pi_cd * 4.0e4
        rad_cd = 50.0e-4
    elif radxnum_cd >= volx34pi_cd * 4.0e8:
        radxnum_cd = volx34pi_cd * 4.0e8
        rad_cd = 0.5e-4
    else:
        rad_cd = radxnum_cd / num_cd

    gasdiffus = 0.557 * (tfld ** 1.75) / press
    gasspeed = 1.455e4 * sqrt(tfld / 98.0)
    knudsen = 3.0 * gasdiffus / (gasspeed * rad_cd)
    fuchs_sutugin = (0.4875 * (1.0 + knudsen)) / (knudsen * (1.184 + knudsen) + 0.4875)
    return 12.56637 * radxnum_cd * gasdiffus * fuchs_sutugin

def sox_cldaero_update_core_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    ntot_amode: int,
    loffset: int,
    id_msa: int,
    id_h2so4: int,
    id_so2: int,
    id_h2o2: int,
    id_nh3: int,
    modeptr_accum: int,
    dtime: float,
    pi_val: float,
    cldfrc_p: cobj,
    xlwc_p: cobj,
    cldnum_p: cobj,
    cfact_p: cobj,
    tfld_p: cobj,
    press_p: cobj,
    delso4_hprxn_p: cobj,
    xh2so4_p: cobj,
    xso4_p: cobj,
    xso4_init_p: cobj,
    nh3g_p: cobj,
    xnh3_p: cobj,
    xnh4c_p: cobj,
    xmsa_p: cobj,
    xso2_p: cobj,
    xh2o2_p: cobj,
    qcw_p: cobj,
    qin_p: cobj,
    dqdt_aqso4_p: cobj,
    dqdt_aqh2so4_p: cobj,
    dqdt_aqhprxn_p: cobj,
    dqdt_aqo3rxn_p: cobj,
    faqgain_msa_p: cobj,
    faqgain_so4_p: cobj,
    qnum_c_p: cobj,
    numptrcw_amode_p: cobj,
    lptr_so4_cw_amode_p: cobj,
    lptr_msa_cw_amode_p: cobj,
    lptr_nh4_cw_amode_p: cobj,
):
    cldfrc = Ptr[float](cldfrc_p)
    xlwc = Ptr[float](xlwc_p)
    cldnum = Ptr[float](cldnum_p)
    cfact = Ptr[float](cfact_p)
    tfld = Ptr[float](tfld_p)
    press = Ptr[float](press_p)
    delso4_hprxn = Ptr[float](delso4_hprxn_p)
    xh2so4 = Ptr[float](xh2so4_p)
    xso4 = Ptr[float](xso4_p)
    xso4_init = Ptr[float](xso4_init_p)
    nh3g = Ptr[float](nh3g_p)
    xnh3 = Ptr[float](xnh3_p)
    xnh4c = Ptr[float](xnh4c_p)
    xmsa = Ptr[float](xmsa_p)
    xso2 = Ptr[float](xso2_p)
    xh2o2 = Ptr[float](xh2o2_p)
    qcw = Ptr[float](qcw_p)
    qin = Ptr[float](qin_p)
    dqdt_aqso4 = Ptr[float](dqdt_aqso4_p)
    dqdt_aqh2so4 = Ptr[float](dqdt_aqh2so4_p)
    dqdt_aqhprxn = Ptr[float](dqdt_aqhprxn_p)
    dqdt_aqo3rxn = Ptr[float](dqdt_aqo3rxn_p)
    faqgain_msa = Ptr[float](faqgain_msa_p)
    faqgain_so4 = Ptr[float](faqgain_so4_p)
    qnum_c = Ptr[float](qnum_c_p)
    numptrcw_amode = Ptr[int](numptrcw_amode_p)
    lptr_so4_cw_amode = Ptr[int](lptr_so4_cw_amode_p)
    lptr_msa_cw_amode = Ptr[int](lptr_msa_cw_amode_p)
    lptr_nh4_cw_amode = Ptr[int](lptr_nh4_cw_amode_p)

    for m in range(1, gas_pcnst + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                dqdt_aqso4[_idx3(i, k, m, ncol, pver)] = 0.0
                dqdt_aqh2so4[_idx3(i, k, m, ncol, pver)] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            dqdt_aqhprxn[_idx2(i, k, ncol)] = 0.0
            dqdt_aqo3rxn[_idx2(i, k, ncol)] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, ncol)
            idxp = _idx2(i, k, pcols)
            if cldfrc[idxp] >= 1.0e-5:
                xl = xlwc[idxp]
                if xl >= 1.0e-8:
                    delso4_o3rxn = xso4[idx] - xso4_init[idx]
                    if id_nh3 > 0:
                        delnh3 = nh3g[idx] - xnh3[idx]
                        delnh4 = -delnh3
                    else:
                        delnh3 = 0.0
                        delnh4 = 0.0

                    for n in range(1, ntot_amode + 1):
                        qnum_c[n - 1] = 0.0
                        l = numptrcw_amode[n - 1] - loffset
                        if l > 0:
                            qnum_c[n - 1] = max(0.0, qcw[_idx3(i, k, l, ncol, pver)])

                    n_accum = modeptr_accum
                    if n_accum <= 0:
                        n_accum = 1
                    qnum_c[n_accum - 1] = max(1.0e-10, qnum_c[n_accum - 1])

                    sumf = 0.0
                    for n in range(1, ntot_amode + 1):
                        faqgain_so4[n - 1] = 0.0
                        if lptr_so4_cw_amode[n - 1] > 0:
                            faqgain_so4[n - 1] = qnum_c[n - 1]
                            sumf = sumf + faqgain_so4[n - 1]

                    if sumf > 0.0:
                        for n in range(1, ntot_amode + 1):
                            faqgain_so4[n - 1] = faqgain_so4[n - 1] / sumf

                    ntot_msa_c = 0
                    sumf = 0.0
                    for n in range(1, ntot_amode + 1):
                        faqgain_msa[n - 1] = 0.0
                        if lptr_msa_cw_amode[n - 1] > 0:
                            faqgain_msa[n - 1] = qnum_c[n - 1]
                            ntot_msa_c = ntot_msa_c + 1
                        sumf = sumf + faqgain_msa[n - 1]

                    if sumf > 0.0:
                        for n in range(1, ntot_amode + 1):
                            faqgain_msa[n - 1] = faqgain_msa[n - 1] / sumf

                    uptkrate = _sox_cldaero_uptakerate(
                        xl, cldnum[idxp], cfact[idx], cldfrc[idxp], tfld[idxp], press[idxp], pi_val
                    )
                    uptkrate = (1.0 - exp(-min(100.0, dtime * uptkrate))) / dtime

                    dso4dt_gasuptk = xh2so4[idx] * uptkrate
                    if id_msa > 0:
                        dmsadt_gasuptk = xmsa[idx] * uptkrate
                    else:
                        dmsadt_gasuptk = 0.0

                    dmsadt_gasuptk_toso4 = 0.0
                    dmsadt_gasuptk_tomsa = dmsadt_gasuptk
                    if ntot_msa_c == 0:
                        dmsadt_gasuptk_tomsa = 0.0
                        dmsadt_gasuptk_toso4 = dmsadt_gasuptk

                    dso4dt_aqrxn = (delso4_o3rxn + delso4_hprxn[idx]) / dtime
                    dso4dt_hprxn = delso4_hprxn[idx] / dtime
                    fwetrem = 0.0

                    for n in range(1, ntot_amode + 1):
                        l = lptr_so4_cw_amode[n - 1] - loffset
                        if l > 0:
                            qidx = _idx3(i, k, l, ncol, pver)
                            dqdt_aqso4[qidx] = faqgain_so4[n - 1] * dso4dt_aqrxn * cldfrc[idxp]
                            dqdt_aqh2so4[qidx] = (
                                faqgain_so4[n - 1] * (dso4dt_gasuptk + dmsadt_gasuptk_toso4) * cldfrc[idxp]
                            )
                            dqdt_aq = dqdt_aqso4[qidx] + dqdt_aqh2so4[qidx]
                            dqdt_wr = -fwetrem * dqdt_aq
                            dqdt = dqdt_aq + dqdt_wr
                            qcw[qidx] = qcw[qidx] + dqdt * dtime

                        l = lptr_msa_cw_amode[n - 1] - loffset
                        if l > 0:
                            qidx = _idx3(i, k, l, ncol, pver)
                            dqdt_aq = faqgain_msa[n - 1] * dmsadt_gasuptk_tomsa * cldfrc[idxp]
                            dqdt_wr = -fwetrem * dqdt_aq
                            dqdt = dqdt_aq + dqdt_wr
                            qcw[qidx] = qcw[qidx] + dqdt * dtime

                        l = lptr_nh4_cw_amode[n - 1] - loffset
                        if l > 0:
                            qidx = _idx3(i, k, l, ncol, pver)
                            if delnh4 > 0.0:
                                dqdt_aq = faqgain_so4[n - 1] * delnh4 / dtime * cldfrc[idxp]
                                dqdt = dqdt_aq
                                qcw[qidx] = qcw[qidx] + dqdt * dtime
                            else:
                                dqdt = (
                                    qcw[qidx]
                                    / max(xnh4c[idxp], 1.0e-35)
                                    * delnh4
                                    / dtime
                                    * cldfrc[idxp]
                                )
                                qcw[qidx] = qcw[qidx] + dqdt * dtime

                    qin[_idx3(i, k, id_h2so4, ncol, pver)] = (
                        qin[_idx3(i, k, id_h2so4, ncol, pver)] - dso4dt_gasuptk * dtime * cldfrc[idxp]
                    )
                    if id_msa > 0:
                        qin[_idx3(i, k, id_msa, ncol, pver)] = (
                            qin[_idx3(i, k, id_msa, ncol, pver)] - dmsadt_gasuptk * dtime * cldfrc[idxp]
                        )

                    fwetrem = 0.0
                    dqdt_wr = -fwetrem * xso2[idx] / dtime * cldfrc[idxp]
                    dqdt_aq = -dso4dt_aqrxn * cldfrc[idxp]
                    dqdt = dqdt_aq + dqdt_wr
                    qin[_idx3(i, k, id_so2, ncol, pver)] = qin[_idx3(i, k, id_so2, ncol, pver)] + dqdt * dtime

                    fwetrem = 0.0
                    dqdt_wr = -fwetrem * xh2o2[idx] / dtime * cldfrc[idxp]
                    dqdt_aq = -dso4dt_hprxn * cldfrc[idxp]
                    dqdt = dqdt_aq + dqdt_wr
                    qin[_idx3(i, k, id_h2o2, ncol, pver)] = (
                        qin[_idx3(i, k, id_h2o2, ncol, pver)] + dqdt * dtime
                    )

                    if id_nh3 > 0:
                        dqdt_aq = delnh3 / dtime * cldfrc[idxp]
                        dqdt = dqdt_aq
                        qin[_idx3(i, k, id_nh3, ncol, pver)] = (
                            qin[_idx3(i, k, id_nh3, ncol, pver)] + dqdt * dtime
                        )

                    dqdt_aqhprxn[idx] = dso4dt_hprxn * cldfrc[idxp]
                    dqdt_aqo3rxn[idx] = (dso4dt_aqrxn - dso4dt_hprxn) * cldfrc[idxp]

def sox_cldaero_finalize_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    ntot_amode: int,
    loffset: int,
    id_so2: int,
    id_nh3: int,
    small_value: float,
    specmw_so4_amode: float,
    gravit: float,
    mbar_p: cobj,
    pdel_p: cobj,
    qcw_p: cobj,
    qin_p: cobj,
    dqdt_aqso4_p: cobj,
    dqdt_aqh2so4_p: cobj,
    dqdt_aqhprxn_p: cobj,
    dqdt_aqo3rxn_p: cobj,
    sflx_aqso4_p: cobj,
    sflx_aqh2so4_p: cobj,
    sflx_aqhprxn_p: cobj,
    sflx_aqo3rxn_p: cobj,
    adv_mass_p: cobj,
    lptr_so4_cw_amode_p: cobj,
    lptr_msa_cw_amode_p: cobj,
    lptr_nh4_cw_amode_p: cobj,
):
    mbar = Ptr[float](mbar_p)
    pdel = Ptr[float](pdel_p)
    qcw = Ptr[float](qcw_p)
    qin = Ptr[float](qin_p)
    dqdt_aqso4 = Ptr[float](dqdt_aqso4_p)
    dqdt_aqh2so4 = Ptr[float](dqdt_aqh2so4_p)
    dqdt_aqhprxn = Ptr[float](dqdt_aqhprxn_p)
    dqdt_aqo3rxn = Ptr[float](dqdt_aqo3rxn_p)
    sflx_aqso4 = Ptr[float](sflx_aqso4_p)
    sflx_aqh2so4 = Ptr[float](sflx_aqh2so4_p)
    sflx_aqhprxn = Ptr[float](sflx_aqhprxn_p)
    sflx_aqo3rxn = Ptr[float](sflx_aqo3rxn_p)
    adv_mass = Ptr[float](adv_mass_p)
    lptr_so4_cw_amode = Ptr[int](lptr_so4_cw_amode_p)
    lptr_msa_cw_amode = Ptr[int](lptr_msa_cw_amode_p)
    lptr_nh4_cw_amode = Ptr[int](lptr_nh4_cw_amode_p)

    for n in range(1, ntot_amode + 1):
        for i in range(1, ncol + 1):
            sflx_aqso4[_idx2(i, n, ncol)] = 0.0
            sflx_aqh2so4[_idx2(i, n, ncol)] = 0.0

    for i in range(1, ncol + 1):
        sflx_aqhprxn[i - 1] = 0.0
        sflx_aqo3rxn[i - 1] = 0.0

    for k in range(1, pver + 1):
        for n in range(1, ntot_amode + 1):
            l = lptr_so4_cw_amode[n - 1] - loffset
            if l > 0:
                for i in range(1, ncol + 1):
                    idx = _idx3(i, k, l, ncol, pver)
                    qcw[idx] = max(qcw[idx], small_value)

            l = lptr_msa_cw_amode[n - 1] - loffset
            if l > 0:
                for i in range(1, ncol + 1):
                    idx = _idx3(i, k, l, ncol, pver)
                    qcw[idx] = max(qcw[idx], small_value)

            l = lptr_nh4_cw_amode[n - 1] - loffset
            if l > 0:
                for i in range(1, ncol + 1):
                    idx = _idx3(i, k, l, ncol, pver)
                    qcw[idx] = max(qcw[idx], small_value)

        for i in range(1, ncol + 1):
            idx_so2 = _idx3(i, k, id_so2, ncol, pver)
            qin[idx_so2] = max(qin[idx_so2], small_value)
            if id_nh3 > 0:
                idx_nh3 = _idx3(i, k, id_nh3, ncol, pver)
                qin[idx_nh3] = max(qin[idx_nh3], small_value)

    for n in range(1, ntot_amode + 1):
        m = lptr_so4_cw_amode[n - 1]
        l = m - loffset
        if l > 0:
            adv = adv_mass[l - 1]
            for i in range(1, ncol + 1):
                sum_aqso4 = 0.0
                sum_aqh2so4 = 0.0
                for k in range(1, pver + 1):
                    midx = _idx2(i, k, ncol)
                    idx = _idx3(i, k, l, ncol, pver)
                    sum_aqso4 = sum_aqso4 + dqdt_aqso4[idx] * adv / mbar[midx] * pdel[midx] / gravit
                    sum_aqh2so4 = sum_aqh2so4 + dqdt_aqh2so4[idx] * adv / mbar[midx] * pdel[midx] / gravit
                sflx_aqso4[_idx2(i, n, ncol)] = sum_aqso4
                sflx_aqh2so4[_idx2(i, n, ncol)] = sum_aqh2so4

    for i in range(1, ncol + 1):
        sum_hprxn = 0.0
        sum_o3rxn = 0.0
        for k in range(1, pver + 1):
            midx = _idx2(i, k, ncol)
            idx = _idx2(i, k, ncol)
            sum_hprxn = sum_hprxn + dqdt_aqhprxn[idx] * specmw_so4_amode / mbar[midx] * pdel[midx] / gravit
            sum_o3rxn = sum_o3rxn + dqdt_aqo3rxn[idx] * specmw_so4_amode / mbar[midx] * pdel[midx] / gravit
        sflx_aqhprxn[i - 1] = sum_hprxn
        sflx_aqo3rxn[i - 1] = sum_o3rxn

def setsox_shell_codon(
    stage: int,
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    nfs: int,
    ntot_amode: int,
    loffset: int,
    itermax: int,
    cloud_borne_flag: int,
    modal_aerosols_flag: int,
    inv_so2_flag: int,
    inv_h2o2_flag: int,
    inv_o3_flag: int,
    inv_ho2_flag: int,
    id_so2: int,
    id_hno3: int,
    id_h2o2: int,
    id_nh3: int,
    id_o3: int,
    id_ho2: int,
    id_h2so4: int,
    id_so4: int,
    id_msa: int,
    modeptr_accum: int,
    dtime: float,
    ph0: float,
    const0: float,
    kh0: float,
    kh1: float,
    kh2: float,
    kh3: float,
    ra: float,
    xkw: float,
    so4_fact: float,
    pi_val: float,
    xhnm_p: cobj,
    invariants_p: cobj,
    qin_p: cobj,
    cfact_p: cobj,
    xph_p: cobj,
    xso2_p: cobj,
    xhno3_p: cobj,
    xh2o2_p: cobj,
    xnh3_p: cobj,
    xo3_p: cobj,
    xho2_p: cobj,
    xh2so4_p: cobj,
    xso4_p: cobj,
    xno3_p: cobj,
    xnh4_p: cobj,
    xmsa_p: cobj,
    press_p: cobj,
    tfld_p: cobj,
    cldfrc_p: cobj,
    cldnum_p: cobj,
    lwc_p: cobj,
    xlwc_p: cobj,
    xso4c_p: cobj,
    xnh4c_p: cobj,
    xno3c_p: cobj,
    xso4_init_p: cobj,
    xdelso4hp_p: cobj,
    hno3g_p: cobj,
    nh3g_p: cobj,
    hehno3_p: cobj,
    heh2o2_p: cobj,
    heso2_p: cobj,
    henh3_p: cobj,
    heo3_p: cobj,
    xphlwc_p: cobj,
    qcw_p: cobj,
    dqdt_aqso4_p: cobj,
    dqdt_aqh2so4_p: cobj,
    dqdt_aqhprxn_p: cobj,
    dqdt_aqo3rxn_p: cobj,
    faqgain_msa_p: cobj,
    faqgain_so4_p: cobj,
    qnum_c_p: cobj,
    numptrcw_amode_p: cobj,
    lptr_so4_cw_amode_p: cobj,
    lptr_msa_cw_amode_p: cobj,
    lptr_nh4_cw_amode_p: cobj,
):
    setsox_init_fields_codon(
        stage,
        ncol,
        pcols,
        pver,
        gas_pcnst,
        nfs,
        cloud_borne_flag,
        inv_so2_flag,
        inv_h2o2_flag,
        inv_o3_flag,
        inv_ho2_flag,
        id_so2,
        id_hno3,
        id_h2o2,
        id_nh3,
        id_o3,
        id_ho2,
        id_h2so4,
        id_so4,
        id_msa,
        ph0,
        xhnm_p,
        invariants_p,
        qin_p,
        cfact_p,
        xph_p,
        xso2_p,
        xhno3_p,
        xh2o2_p,
        xnh3_p,
        xo3_p,
        xho2_p,
        xh2so4_p,
        xso4_p,
        xno3_p,
        xnh4_p,
        xmsa_p,
    )
    if stage == 1:
        return

    setsox_ph_solve_codon(
        ncol,
        pcols,
        pver,
        itermax,
        cloud_borne_flag,
        const0,
        ra,
        xkw,
        so4_fact,
        press_p,
        tfld_p,
        cldfrc_p,
        xhnm_p,
        xlwc_p,
        xso4c_p,
        xnh4c_p,
        xno3c_p,
        xso4_p,
        xnh4_p,
        xno3_p,
        xso2_p,
        xhno3_p,
        xnh3_p,
        xph_p,
    )
    setsox_aqchem_predict_codon(
        ncol,
        pcols,
        pver,
        cloud_borne_flag,
        modal_aerosols_flag,
        id_nh3,
        dtime,
        const0,
        kh0,
        kh1,
        kh2,
        kh3,
        ra,
        xkw,
        press_p,
        tfld_p,
        xhnm_p,
        xlwc_p,
        xph_p,
        xho2_p,
        xhno3_p,
        xno3_p,
        xh2o2_p,
        xso2_p,
        xo3_p,
        xnh3_p,
        xnh4_p,
        xso4_p,
        xso4_init_p,
        xdelso4hp_p,
        hno3g_p,
        nh3g_p,
        hehno3_p,
        heh2o2_p,
        heso2_p,
        henh3_p,
        heo3_p,
    )
    sox_cldaero_update_core_codon(
        ncol,
        pcols,
        pver,
        gas_pcnst,
        ntot_amode,
        loffset,
        id_msa,
        id_h2so4,
        id_so2,
        id_h2o2,
        id_nh3,
        modeptr_accum,
        dtime,
        pi_val,
        cldfrc_p,
        xlwc_p,
        cldnum_p,
        cfact_p,
        tfld_p,
        press_p,
        xdelso4hp_p,
        xh2so4_p,
        xso4_p,
        xso4_init_p,
        nh3g_p,
        xnh3_p,
        xnh4c_p,
        xmsa_p,
        xso2_p,
        xh2o2_p,
        qcw_p,
        qin_p,
        dqdt_aqso4_p,
        dqdt_aqh2so4_p,
        dqdt_aqhprxn_p,
        dqdt_aqo3rxn_p,
        faqgain_msa_p,
        faqgain_so4_p,
        qnum_c_p,
        numptrcw_amode_p,
        lptr_so4_cw_amode_p,
        lptr_msa_cw_amode_p,
        lptr_nh4_cw_amode_p,
    )
    setsox_xph_lwc_diag_codon(ncol, pcols, pver, cldfrc_p, lwc_p, xph_p, xphlwc_p)
