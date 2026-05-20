from math import exp

@inline
def _idx1(i: int) -> int:
    return i - 1


@inline
def _idx2(i: int, k: int, pcols: int) -> int:
    """modal_aer_opt arrays declared as (pcols,pver)."""
    return (i - 1) + (k - 1) * pcols


@inline
def _sw_idx(i: int, k0: int, band: int, pcols: int, pverp: int) -> int:
    """modal_aero_sw output arrays declared as (pcols,0:pver,nswbands)."""
    return (i - 1) + k0 * pcols + (band - 1) * pcols * pverp


@inline
def _cheb_idx(nc: int, i: int, k: int, ncoef: int, pcols: int) -> int:
    """modal_aero_sw Chebyshev arrays declared as (ncoef,pcols,pver)."""
    return (nc - 1) + (i - 1) * ncoef + (k - 1) * ncoef * pcols


def modal_aer_opt_sw_init_state_codon(
    ncol: int,
    pcols: int,
    pver: int,
    nswbands: int,
    rga: float,
    rair: float,
    pdeldry_p: cobj,
    pmid_p: cobj,
    state_t_p: cobj,
    tauxar_p: cobj,
    wa_p: cobj,
    ga_p: cobj,
    fa_p: cobj,
    mass_p: cobj,
    air_density_p: cobj,
):
    pdeldry = Ptr[float](pdeldry_p)
    pmid = Ptr[float](pmid_p)
    state_t = Ptr[float](state_t_p)
    tauxar = Ptr[float](tauxar_p)
    wa = Ptr[float](wa_p)
    ga = Ptr[float](ga_p)
    fa = Ptr[float](fa_p)
    mass = Ptr[float](mass_p)
    air_density = Ptr[float](air_density_p)

    pverp = pver + 1
    for isw in range(1, nswbands + 1):
        for k0 in range(0, pver + 1):
            for i in range(1, ncol + 1):
                idx3 = _sw_idx(i, k0, isw, pcols, pverp)
                tauxar[idx3] = 0.0
                wa[idx3] = 0.0
                ga[idx3] = 0.0
                fa[idx3] = 0.0

        for i in range(1, ncol + 1):
            idx0 = _sw_idx(i, 0, isw, pcols, pverp)
            tauxar[idx0] = 0.0
            wa[idx0] = 0.925
            ga[idx0] = 0.850
            fa[idx0] = 0.7225

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            mass[idx] = pdeldry[idx] * rga
            air_density[idx] = pmid[idx] / (rair * state_t[idx])


def modal_aer_opt_sw_zero_diagnostics_codon(
    ncol: int,
    pcols: int,
    pver: int,
    extinct_p: cobj,
    absorb_p: cobj,
    extinctuv_p: cobj,
    extinctnir_p: cobj,
    aodvis_p: cobj,
    aodvisst_p: cobj,
    aodabs_p: cobj,
    aodabsbc_p: cobj,
    ssavis_p: cobj,
    burdendust_p: cobj,
    burdenso4_p: cobj,
    burdenpom_p: cobj,
    burdensoa_p: cobj,
    burdenbc_p: cobj,
    burdenseasalt_p: cobj,
    dustaod_p: cobj,
    so4aod_p: cobj,
    pomaod_p: cobj,
    soaaod_p: cobj,
    bcaod_p: cobj,
    seasaltaod_p: cobj,
    aoduv_p: cobj,
    aodnir_p: cobj,
    aoduvst_p: cobj,
    aodnirst_p: cobj,
):
    extinct = Ptr[float](extinct_p)
    absorb = Ptr[float](absorb_p)
    extinctuv = Ptr[float](extinctuv_p)
    extinctnir = Ptr[float](extinctnir_p)
    aodvis = Ptr[float](aodvis_p)
    aodvisst = Ptr[float](aodvisst_p)
    aodabs = Ptr[float](aodabs_p)
    aodabsbc = Ptr[float](aodabsbc_p)
    ssavis = Ptr[float](ssavis_p)
    burdendust = Ptr[float](burdendust_p)
    burdenso4 = Ptr[float](burdenso4_p)
    burdenpom = Ptr[float](burdenpom_p)
    burdensoa = Ptr[float](burdensoa_p)
    burdenbc = Ptr[float](burdenbc_p)
    burdenseasalt = Ptr[float](burdenseasalt_p)
    dustaod = Ptr[float](dustaod_p)
    so4aod = Ptr[float](so4aod_p)
    pomaod = Ptr[float](pomaod_p)
    soaaod = Ptr[float](soaaod_p)
    bcaod = Ptr[float](bcaod_p)
    seasaltaod = Ptr[float](seasaltaod_p)
    aoduv = Ptr[float](aoduv_p)
    aodnir = Ptr[float](aodnir_p)
    aoduvst = Ptr[float](aoduvst_p)
    aodnirst = Ptr[float](aodnirst_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            extinct[idx] = 0.0
            absorb[idx] = 0.0
            extinctuv[idx] = 0.0
            extinctnir[idx] = 0.0

    for i in range(1, ncol + 1):
        idx = _idx1(i)
        aodvis[idx] = 0.0
        aodvisst[idx] = 0.0
        aodabs[idx] = 0.0
        aodabsbc[idx] = 0.0
        ssavis[idx] = 0.0
        burdendust[idx] = 0.0
        burdenso4[idx] = 0.0
        burdenpom[idx] = 0.0
        burdensoa[idx] = 0.0
        burdenbc[idx] = 0.0
        burdenseasalt[idx] = 0.0
        dustaod[idx] = 0.0
        so4aod[idx] = 0.0
        pomaod[idx] = 0.0
        soaaod[idx] = 0.0
        bcaod[idx] = 0.0
        seasaltaod[idx] = 0.0
        aoduv[idx] = 0.0
        aodnir[idx] = 0.0
        aoduvst[idx] = 0.0
        aodnirst[idx] = 0.0


def modal_aer_opt_sw_mode_diag_init_codon(
    ncol: int,
    burden_p: cobj,
    aodmode_p: cobj,
    dustaodmode_p: cobj,
):
    burden = Ptr[float](burden_p)
    aodmode = Ptr[float](aodmode_p)
    dustaodmode = Ptr[float](dustaodmode_p)

    for i in range(1, ncol + 1):
        idx = _idx1(i)
        burden[idx] = 0.0
        aodmode[idx] = 0.0
        dustaodmode[idx] = 0.0


def modal_aer_opt_sw_reset_layer_codon(
    ncol: int,
    dryvol_p: cobj,
    dustvol_p: cobj,
    scatdust_p: cobj,
    absdust_p: cobj,
    hygrodust_p: cobj,
    scatso4_p: cobj,
    absso4_p: cobj,
    hygroso4_p: cobj,
    scatbc_p: cobj,
    absbc_p: cobj,
    hygrobc_p: cobj,
    scatpom_p: cobj,
    abspom_p: cobj,
    hygropom_p: cobj,
    scatsoa_p: cobj,
    abssoa_p: cobj,
    hygrosoa_p: cobj,
    scatseasalt_p: cobj,
    absseasalt_p: cobj,
    hygroseasalt_p: cobj,
):
    dryvol = Ptr[float](dryvol_p)
    dustvol = Ptr[float](dustvol_p)
    scatdust = Ptr[float](scatdust_p)
    absdust = Ptr[float](absdust_p)
    hygrodust = Ptr[float](hygrodust_p)
    scatso4 = Ptr[float](scatso4_p)
    absso4 = Ptr[float](absso4_p)
    hygroso4 = Ptr[float](hygroso4_p)
    scatbc = Ptr[float](scatbc_p)
    absbc = Ptr[float](absbc_p)
    hygrobc = Ptr[float](hygrobc_p)
    scatpom = Ptr[float](scatpom_p)
    abspom = Ptr[float](abspom_p)
    hygropom = Ptr[float](hygropom_p)
    scatsoa = Ptr[float](scatsoa_p)
    abssoa = Ptr[float](abssoa_p)
    hygrosoa = Ptr[float](hygrosoa_p)
    scatseasalt = Ptr[float](scatseasalt_p)
    absseasalt = Ptr[float](absseasalt_p)
    hygroseasalt = Ptr[float](hygroseasalt_p)

    for i in range(1, ncol + 1):
        idx = _idx1(i)
        dryvol[idx] = 0.0
        dustvol[idx] = 0.0
        scatdust[idx] = 0.0
        absdust[idx] = 0.0
        hygrodust[idx] = 0.0
        scatso4[idx] = 0.0
        absso4[idx] = 0.0
        hygroso4[idx] = 0.0
        scatbc[idx] = 0.0
        absbc[idx] = 0.0
        hygrobc[idx] = 0.0
        scatpom[idx] = 0.0
        abspom[idx] = 0.0
        hygropom[idx] = 0.0
        scatsoa[idx] = 0.0
        abssoa[idx] = 0.0
        hygrosoa[idx] = 0.0
        scatseasalt[idx] = 0.0
        absseasalt[idx] = 0.0
        hygroseasalt[idx] = 0.0


def modal_aer_opt_sw_species_volume_codon(
    ncol: int,
    pcols: int,
    k: int,
    specdens: float,
    specmmr_p: cobj,
    vol_p: cobj,
    dryvol_p: cobj,
):
    specmmr = Ptr[float](specmmr_p)
    vol = Ptr[float](vol_p)
    dryvol = Ptr[float](dryvol_p)

    for i in range(1, ncol + 1):
        idx1 = _idx1(i)
        idx2 = _idx2(i, k, pcols)
        vol[idx1] = specmmr[idx2] / specdens
        dryvol[idx1] = dryvol[idx1] + vol[idx1]


def modal_aer_opt_sw_species_vis_diag_codon(
    ncol: int,
    pcols: int,
    k: int,
    spectype_code: int,
    specrefr: float,
    specrefi: float,
    hygro_aer: float,
    specmmr_p: cobj,
    mass_p: cobj,
    vol_p: cobj,
    burden_p: cobj,
    burdendust_p: cobj,
    burdenso4_p: cobj,
    burdenbc_p: cobj,
    burdenpom_p: cobj,
    burdensoa_p: cobj,
    burdenseasalt_p: cobj,
    dustvol_p: cobj,
    scatdust_p: cobj,
    absdust_p: cobj,
    hygrodust_p: cobj,
    scatso4_p: cobj,
    absso4_p: cobj,
    hygroso4_p: cobj,
    scatbc_p: cobj,
    absbc_p: cobj,
    hygrobc_p: cobj,
    scatpom_p: cobj,
    abspom_p: cobj,
    hygropom_p: cobj,
    scatsoa_p: cobj,
    abssoa_p: cobj,
    hygrosoa_p: cobj,
    scatseasalt_p: cobj,
    absseasalt_p: cobj,
    hygroseasalt_p: cobj,
):
    specmmr = Ptr[float](specmmr_p)
    mass = Ptr[float](mass_p)
    vol = Ptr[float](vol_p)
    burden = Ptr[float](burden_p)
    burdendust = Ptr[float](burdendust_p)
    burdenso4 = Ptr[float](burdenso4_p)
    burdenbc = Ptr[float](burdenbc_p)
    burdenpom = Ptr[float](burdenpom_p)
    burdensoa = Ptr[float](burdensoa_p)
    burdenseasalt = Ptr[float](burdenseasalt_p)
    dustvol = Ptr[float](dustvol_p)
    scatdust = Ptr[float](scatdust_p)
    absdust = Ptr[float](absdust_p)
    hygrodust = Ptr[float](hygrodust_p)
    scatso4 = Ptr[float](scatso4_p)
    absso4 = Ptr[float](absso4_p)
    hygroso4 = Ptr[float](hygroso4_p)
    scatbc = Ptr[float](scatbc_p)
    absbc = Ptr[float](absbc_p)
    hygrobc = Ptr[float](hygrobc_p)
    scatpom = Ptr[float](scatpom_p)
    abspom = Ptr[float](abspom_p)
    hygropom = Ptr[float](hygropom_p)
    scatsoa = Ptr[float](scatsoa_p)
    abssoa = Ptr[float](abssoa_p)
    hygrosoa = Ptr[float](hygrosoa_p)
    scatseasalt = Ptr[float](scatseasalt_p)
    absseasalt = Ptr[float](absseasalt_p)
    hygroseasalt = Ptr[float](hygroseasalt_p)

    for i in range(1, ncol + 1):
        idx1 = _idx1(i)
        idx2 = _idx2(i, k, pcols)
        burden[idx1] = burden[idx1] + specmmr[idx2] * mass[idx2]

    if spectype_code == 1:
        for i in range(1, ncol + 1):
            idx1 = _idx1(i)
            idx2 = _idx2(i, k, pcols)
            burdendust[idx1] = burdendust[idx1] + specmmr[idx2] * mass[idx2]
            dustvol[idx1] = vol[idx1]
            scatdust[idx1] = vol[idx1] * specrefr
            absdust[idx1] = -vol[idx1] * specrefi
            hygrodust[idx1] = vol[idx1] * hygro_aer
    if spectype_code == 2:
        for i in range(1, ncol + 1):
            idx1 = _idx1(i)
            idx2 = _idx2(i, k, pcols)
            burdenso4[idx1] = burdenso4[idx1] + specmmr[idx2] * mass[idx2]
            scatso4[idx1] = vol[idx1] * specrefr
            absso4[idx1] = -vol[idx1] * specrefi
            hygroso4[idx1] = vol[idx1] * hygro_aer
    if spectype_code == 3:
        for i in range(1, ncol + 1):
            idx1 = _idx1(i)
            idx2 = _idx2(i, k, pcols)
            burdenbc[idx1] = burdenbc[idx1] + specmmr[idx2] * mass[idx2]
            scatbc[idx1] = vol[idx1] * specrefr
            absbc[idx1] = -vol[idx1] * specrefi
            hygrobc[idx1] = vol[idx1] * hygro_aer
    if spectype_code == 4:
        for i in range(1, ncol + 1):
            idx1 = _idx1(i)
            idx2 = _idx2(i, k, pcols)
            burdenpom[idx1] = burdenpom[idx1] + specmmr[idx2] * mass[idx2]
            scatpom[idx1] = vol[idx1] * specrefr
            abspom[idx1] = -vol[idx1] * specrefi
            hygropom[idx1] = vol[idx1] * hygro_aer
    if spectype_code == 5:
        for i in range(1, ncol + 1):
            idx1 = _idx1(i)
            idx2 = _idx2(i, k, pcols)
            burdensoa[idx1] = burdensoa[idx1] + specmmr[idx2] * mass[idx2]
            scatsoa[idx1] = vol[idx1] * specrefr
            abssoa[idx1] = -vol[idx1] * specrefi
            hygrosoa[idx1] = vol[idx1] * hygro_aer
    if spectype_code == 6:
        for i in range(1, ncol + 1):
            idx1 = _idx1(i)
            idx2 = _idx2(i, k, pcols)
            burdenseasalt[idx1] = burdenseasalt[idx1] + specmmr[idx2] * mass[idx2]
            scatseasalt[idx1] = vol[idx1] * specrefr
            absseasalt[idx1] = -vol[idx1] * specrefi
            hygroseasalt[idx1] = vol[idx1] * hygro_aer


def modal_aer_opt_sw_optics_props_codon(
    ncol: int,
    pcols: int,
    k: int,
    ncoef: int,
    xrmax: float,
    rhoh2o: float,
    radsurf_p: cobj,
    logradsurf_p: cobj,
    cheb_p: cobj,
    cext_p: cobj,
    cabs_p: cobj,
    casm_p: cobj,
    wetvol_p: cobj,
    mass_p: cobj,
    pext_p: cobj,
    specpext_p: cobj,
    pabs_p: cobj,
    pasm_p: cobj,
    palb_p: cobj,
    dopaer_p: cobj,
):
    radsurf = Ptr[float](radsurf_p)
    logradsurf = Ptr[float](logradsurf_p)
    cheb = Ptr[float](cheb_p)
    cext = Ptr[float](cext_p)
    cabs = Ptr[float](cabs_p)
    casm = Ptr[float](casm_p)
    wetvol = Ptr[float](wetvol_p)
    mass = Ptr[float](mass_p)
    pext = Ptr[float](pext_p)
    specpext = Ptr[float](specpext_p)
    pabs = Ptr[float](pabs_p)
    pasm = Ptr[float](pasm_p)
    palb = Ptr[float](palb_p)
    dopaer = Ptr[float](dopaer_p)

    for i in range(1, ncol + 1):
        idx1 = _idx1(i)
        idx2 = _idx2(i, k, pcols)

        if logradsurf[idx2] <= xrmax:
            pext[idx1] = 0.5 * cext[_idx2(i, 1, pcols)]
            for nc in range(2, ncoef + 1):
                pext[idx1] = pext[idx1] + cheb[_cheb_idx(nc, i, k, ncoef, pcols)] * cext[_idx2(i, nc, pcols)]
            pext[idx1] = exp(pext[idx1])
        else:
            pext[idx1] = 1.5 / (radsurf[idx2] * rhoh2o)

        specpext[idx1] = pext[idx1]
        pext[idx1] = pext[idx1] * wetvol[idx1] * rhoh2o
        pabs[idx1] = 0.5 * cabs[_idx2(i, 1, pcols)]
        pasm[idx1] = 0.5 * casm[_idx2(i, 1, pcols)]
        for nc in range(2, ncoef + 1):
            pabs[idx1] = pabs[idx1] + cheb[_cheb_idx(nc, i, k, ncoef, pcols)] * cabs[_idx2(i, nc, pcols)]
            pasm[idx1] = pasm[idx1] + cheb[_cheb_idx(nc, i, k, ncoef, pcols)] * casm[_idx2(i, nc, pcols)]
        pabs[idx1] = pabs[idx1] * wetvol[idx1] * rhoh2o
        pabs[idx1] = max(0.0, pabs[idx1])
        pabs[idx1] = min(pext[idx1], pabs[idx1])

        palb[idx1] = 1.0 - pabs[idx1] / max(pext[idx1], 1.0e-40)
        palb[idx1] = 1.0 - pabs[idx1] / max(pext[idx1], 1.0e-40)

        dopaer[idx1] = pext[idx1] * mass[idx2]


def modal_aer_opt_sw_mode_diag_night_codon(
    nnite: int,
    fillvalue: float,
    idxnite_p: cobj,
    burden_p: cobj,
    aodmode_p: cobj,
    dustaodmode_p: cobj,
):
    idxnite = Ptr[i32](idxnite_p)
    burden = Ptr[float](burden_p)
    aodmode = Ptr[float](aodmode_p)
    dustaodmode = Ptr[float](dustaodmode_p)

    for i in range(nnite):
        idx = int(idxnite[i]) - 1
        burden[idx] = fillvalue
        aodmode[idx] = fillvalue
        dustaodmode[idx] = fillvalue


def modal_aer_opt_sw_sum_diag_night_codon(
    nnite: int,
    pcols: int,
    pver: int,
    fillvalue: float,
    idxnite_p: cobj,
    extinct_p: cobj,
    absorb_p: cobj,
    aodvis_p: cobj,
    aodabs_p: cobj,
    aodvisst_p: cobj,
):
    idxnite = Ptr[i32](idxnite_p)
    extinct = Ptr[float](extinct_p)
    absorb = Ptr[float](absorb_p)
    aodvis = Ptr[float](aodvis_p)
    aodabs = Ptr[float](aodabs_p)
    aodvisst = Ptr[float](aodvisst_p)

    for i in range(nnite):
        col = int(idxnite[i])
        idx1 = col - 1
        aodvis[idx1] = fillvalue
        aodabs[idx1] = fillvalue
        aodvisst[idx1] = fillvalue
        for k in range(1, pver + 1):
            idx2 = _idx2(col, k, pcols)
            extinct[idx2] = fillvalue
            absorb[idx2] = fillvalue


def modal_aer_opt_sw_finalize_ssavis_codon(
    ncol: int,
    aodvis_p: cobj,
    ssavis_p: cobj,
):
    aodvis = Ptr[float](aodvis_p)
    ssavis = Ptr[float](ssavis_p)

    for i in range(1, ncol + 1):
        idx = _idx1(i)
        if aodvis[idx] > 1.0e-10:
            ssavis[idx] = ssavis[idx] / aodvis[idx]
        else:
            ssavis[idx] = 0.925


def modal_aer_opt_sw_climate_diag_night_codon(
    nnite: int,
    pcols: int,
    pver: int,
    fillvalue: float,
    idxnite_p: cobj,
    ssavis_p: cobj,
    aoduv_p: cobj,
    aodnir_p: cobj,
    aoduvst_p: cobj,
    aodnirst_p: cobj,
    extinctuv_p: cobj,
    extinctnir_p: cobj,
    burdendust_p: cobj,
    burdenso4_p: cobj,
    burdenpom_p: cobj,
    burdensoa_p: cobj,
    burdenbc_p: cobj,
    burdenseasalt_p: cobj,
    aodabsbc_p: cobj,
    dustaod_p: cobj,
    so4aod_p: cobj,
    pomaod_p: cobj,
    soaaod_p: cobj,
    bcaod_p: cobj,
    seasaltaod_p: cobj,
):
    idxnite = Ptr[i32](idxnite_p)
    ssavis = Ptr[float](ssavis_p)
    aoduv = Ptr[float](aoduv_p)
    aodnir = Ptr[float](aodnir_p)
    aoduvst = Ptr[float](aoduvst_p)
    aodnirst = Ptr[float](aodnirst_p)
    extinctuv = Ptr[float](extinctuv_p)
    extinctnir = Ptr[float](extinctnir_p)
    burdendust = Ptr[float](burdendust_p)
    burdenso4 = Ptr[float](burdenso4_p)
    burdenpom = Ptr[float](burdenpom_p)
    burdensoa = Ptr[float](burdensoa_p)
    burdenbc = Ptr[float](burdenbc_p)
    burdenseasalt = Ptr[float](burdenseasalt_p)
    aodabsbc = Ptr[float](aodabsbc_p)
    dustaod = Ptr[float](dustaod_p)
    so4aod = Ptr[float](so4aod_p)
    pomaod = Ptr[float](pomaod_p)
    soaaod = Ptr[float](soaaod_p)
    bcaod = Ptr[float](bcaod_p)
    seasaltaod = Ptr[float](seasaltaod_p)

    for i in range(nnite):
        col = int(idxnite[i])
        idx1 = col - 1
        ssavis[idx1] = fillvalue
        aoduv[idx1] = fillvalue
        aodnir[idx1] = fillvalue
        aoduvst[idx1] = fillvalue
        aodnirst[idx1] = fillvalue
        burdendust[idx1] = fillvalue
        burdenso4[idx1] = fillvalue
        burdenpom[idx1] = fillvalue
        burdensoa[idx1] = fillvalue
        burdenbc[idx1] = fillvalue
        burdenseasalt[idx1] = fillvalue
        aodabsbc[idx1] = fillvalue
        dustaod[idx1] = fillvalue
        so4aod[idx1] = fillvalue
        pomaod[idx1] = fillvalue
        soaaod[idx1] = fillvalue
        bcaod[idx1] = fillvalue
        seasaltaod[idx1] = fillvalue

        for k in range(1, pver + 1):
            idx2 = _idx2(col, k, pcols)
            extinctuv[idx2] = fillvalue
            extinctnir[idx2] = fillvalue


def modal_aer_opt_sw_accumulate_diagnostics_codon(
    ncol: int,
    pcols: int,
    k: int,
    do_uv: int,
    do_nir: int,
    do_vis: int,
    crefwsw_re: float,
    crefwsw_im: float,
    troplev_p: cobj,
    mass_p: cobj,
    air_density_p: cobj,
    dopaer_p: cobj,
    pabs_p: cobj,
    palb_p: cobj,
    wetvol_p: cobj,
    watervol_p: cobj,
    dustvol_p: cobj,
    scatdust_p: cobj,
    scatso4_p: cobj,
    scatbc_p: cobj,
    scatpom_p: cobj,
    scatsoa_p: cobj,
    scatseasalt_p: cobj,
    absdust_p: cobj,
    absso4_p: cobj,
    absbc_p: cobj,
    abspom_p: cobj,
    abssoa_p: cobj,
    absseasalt_p: cobj,
    hygrodust_p: cobj,
    hygroso4_p: cobj,
    hygrobc_p: cobj,
    hygropom_p: cobj,
    hygrosoa_p: cobj,
    hygroseasalt_p: cobj,
    extinctuv_p: cobj,
    aoduv_p: cobj,
    aoduvst_p: cobj,
    extinctnir_p: cobj,
    aodnir_p: cobj,
    aodnirst_p: cobj,
    extinct_p: cobj,
    absorb_p: cobj,
    aodvis_p: cobj,
    aodabs_p: cobj,
    aodmode_p: cobj,
    ssavis_p: cobj,
    aodvisst_p: cobj,
    dustaodmode_p: cobj,
    aodabsbc_p: cobj,
    dustaod_p: cobj,
    so4aod_p: cobj,
    pomaod_p: cobj,
    soaaod_p: cobj,
    bcaod_p: cobj,
    seasaltaod_p: cobj,
):
    troplev = Ptr[i32](troplev_p)
    mass = Ptr[float](mass_p)
    air_density = Ptr[float](air_density_p)
    dopaer = Ptr[float](dopaer_p)
    pabs = Ptr[float](pabs_p)
    palb = Ptr[float](palb_p)
    wetvol = Ptr[float](wetvol_p)
    watervol = Ptr[float](watervol_p)
    dustvol = Ptr[float](dustvol_p)
    scatdust = Ptr[float](scatdust_p)
    scatso4 = Ptr[float](scatso4_p)
    scatbc = Ptr[float](scatbc_p)
    scatpom = Ptr[float](scatpom_p)
    scatsoa = Ptr[float](scatsoa_p)
    scatseasalt = Ptr[float](scatseasalt_p)
    absdust = Ptr[float](absdust_p)
    absso4 = Ptr[float](absso4_p)
    absbc = Ptr[float](absbc_p)
    abspom = Ptr[float](abspom_p)
    abssoa = Ptr[float](abssoa_p)
    absseasalt = Ptr[float](absseasalt_p)
    hygrodust = Ptr[float](hygrodust_p)
    hygroso4 = Ptr[float](hygroso4_p)
    hygrobc = Ptr[float](hygrobc_p)
    hygropom = Ptr[float](hygropom_p)
    hygrosoa = Ptr[float](hygrosoa_p)
    hygroseasalt = Ptr[float](hygroseasalt_p)
    extinctuv = Ptr[float](extinctuv_p)
    aoduv = Ptr[float](aoduv_p)
    aoduvst = Ptr[float](aoduvst_p)
    extinctnir = Ptr[float](extinctnir_p)
    aodnir = Ptr[float](aodnir_p)
    aodnirst = Ptr[float](aodnirst_p)
    extinct = Ptr[float](extinct_p)
    absorb = Ptr[float](absorb_p)
    aodvis = Ptr[float](aodvis_p)
    aodabs = Ptr[float](aodabs_p)
    aodmode = Ptr[float](aodmode_p)
    ssavis = Ptr[float](ssavis_p)
    aodvisst = Ptr[float](aodvisst_p)
    dustaodmode = Ptr[float](dustaodmode_p)
    aodabsbc = Ptr[float](aodabsbc_p)
    dustaod = Ptr[float](dustaod_p)
    so4aod = Ptr[float](so4aod_p)
    pomaod = Ptr[float](pomaod_p)
    soaaod = Ptr[float](soaaod_p)
    bcaod = Ptr[float](bcaod_p)
    seasaltaod = Ptr[float](seasaltaod_p)

    for i in range(1, ncol + 1):
        idx1 = _idx1(i)
        idx2 = _idx2(i, k, pcols)

        if do_uv != 0:
            extinctuv[idx2] = extinctuv[idx2] + dopaer[idx1] * air_density[idx2] / mass[idx2]
            aoduv[idx1] = aoduv[idx1] + dopaer[idx1]
            if k <= int(troplev[idx1]):
                aoduvst[idx1] = aoduvst[idx1] + dopaer[idx1]

        if do_nir != 0:
            extinctnir[idx2] = extinctnir[idx2] + dopaer[idx1] * air_density[idx2] / mass[idx2]
            aodnir[idx1] = aodnir[idx1] + dopaer[idx1]
            if k <= int(troplev[idx1]):
                aodnirst[idx1] = aodnirst[idx1] + dopaer[idx1]

        if do_vis != 0:
            extinct[idx2] = extinct[idx2] + dopaer[idx1] * air_density[idx2] / mass[idx2]
            absorb[idx2] = absorb[idx2] + pabs[idx1] * air_density[idx2]
            aodvis[idx1] = aodvis[idx1] + dopaer[idx1]
            aodabs[idx1] = aodabs[idx1] + pabs[idx1] * mass[idx2]
            aodmode[idx1] = aodmode[idx1] + dopaer[idx1]
            ssavis[idx1] = ssavis[idx1] + dopaer[idx1] * palb[idx1]
            if k <= int(troplev[idx1]):
                aodvisst[idx1] = aodvisst[idx1] + dopaer[idx1]

            if wetvol[idx1] > 1.0e-40:
                dustaodmode[idx1] = dustaodmode[idx1] + dopaer[idx1] * dustvol[idx1] / wetvol[idx1]

                scath2o = watervol[idx1] * crefwsw_re
                absh2o = -watervol[idx1] * crefwsw_im
                sumscat = scatso4[idx1] + scatpom[idx1]
                sumscat = sumscat + scatsoa[idx1]
                sumscat = sumscat + scatbc[idx1]
                sumscat = sumscat + scatdust[idx1]
                sumscat = sumscat + scatseasalt[idx1]
                sumscat = sumscat + scath2o
                sumabs = absso4[idx1] + abspom[idx1]
                sumabs = sumabs + abssoa[idx1]
                sumabs = sumabs + absbc[idx1]
                sumabs = sumabs + absdust[idx1]
                sumabs = sumabs + absseasalt[idx1]
                sumabs = sumabs + absh2o
                sumhygro = hygroso4[idx1] + hygropom[idx1]
                sumhygro = sumhygro + hygrosoa[idx1]
                sumhygro = sumhygro + hygrobc[idx1]
                sumhygro = sumhygro + hygrodust[idx1]
                sumhygro = sumhygro + hygroseasalt[idx1]

                scatdust[idx1] = (scatdust[idx1] + scath2o * hygrodust[idx1] / sumhygro) / sumscat
                absdust[idx1] = (absdust[idx1] + absh2o * hygrodust[idx1] / sumhygro) / sumabs

                scatso4[idx1] = (scatso4[idx1] + scath2o * hygroso4[idx1] / sumhygro) / sumscat
                absso4[idx1] = (absso4[idx1] + absh2o * hygroso4[idx1] / sumhygro) / sumabs

                scatpom[idx1] = (scatpom[idx1] + scath2o * hygropom[idx1] / sumhygro) / sumscat
                abspom[idx1] = (abspom[idx1] + absh2o * hygropom[idx1] / sumhygro) / sumabs

                scatsoa[idx1] = (scatsoa[idx1] + scath2o * hygrosoa[idx1] / sumhygro) / sumscat
                abssoa[idx1] = (abssoa[idx1] + absh2o * hygrosoa[idx1] / sumhygro) / sumabs

                scatbc[idx1] = (scatbc[idx1] + scath2o * hygrobc[idx1] / sumhygro) / sumscat
                absbc[idx1] = (absbc[idx1] + absh2o * hygrobc[idx1] / sumhygro) / sumabs

                scatseasalt[idx1] = (scatseasalt[idx1] + scath2o * hygroseasalt[idx1] / sumhygro) / sumscat
                absseasalt[idx1] = (absseasalt[idx1] + absh2o * hygroseasalt[idx1] / sumhygro) / sumabs

                aodabsbc[idx1] = aodabsbc[idx1] + absbc[idx1] * dopaer[idx1] * (1.0 - palb[idx1])

                aodc = (absdust[idx1] * (1.0 - palb[idx1]) + palb[idx1] * scatdust[idx1]) * dopaer[idx1]
                dustaod[idx1] = dustaod[idx1] + aodc

                aodc = (absso4[idx1] * (1.0 - palb[idx1]) + palb[idx1] * scatso4[idx1]) * dopaer[idx1]
                so4aod[idx1] = so4aod[idx1] + aodc

                aodc = (abspom[idx1] * (1.0 - palb[idx1]) + palb[idx1] * scatpom[idx1]) * dopaer[idx1]
                pomaod[idx1] = pomaod[idx1] + aodc

                aodc = (abssoa[idx1] * (1.0 - palb[idx1]) + palb[idx1] * scatsoa[idx1]) * dopaer[idx1]
                soaaod[idx1] = soaaod[idx1] + aodc

                aodc = (absbc[idx1] * (1.0 - palb[idx1]) + palb[idx1] * scatbc[idx1]) * dopaer[idx1]
                bcaod[idx1] = bcaod[idx1] + aodc

                aodc = (absseasalt[idx1] * (1.0 - palb[idx1]) + palb[idx1] * scatseasalt[idx1]) * dopaer[idx1]
                seasaltaod[idx1] = seasaltaod[idx1] + aodc


def modal_aer_opt_sw_accumulate_tau_codon(
    ncol: int,
    pcols: int,
    pver: int,
    k: int,
    isw: int,
    dopaer_p: cobj,
    palb_p: cobj,
    pasm_p: cobj,
    tauxar_p: cobj,
    wa_p: cobj,
    ga_p: cobj,
    fa_p: cobj,
):
    dopaer = Ptr[float](dopaer_p)
    palb = Ptr[float](palb_p)
    pasm = Ptr[float](pasm_p)
    tauxar = Ptr[float](tauxar_p)
    wa = Ptr[float](wa_p)
    ga = Ptr[float](ga_p)
    fa = Ptr[float](fa_p)

    pverp = pver + 1
    for i in range(1, ncol + 1):
        idx1 = _idx1(i)
        idx3 = _sw_idx(i, k, isw, pcols, pverp)
        tauxar[idx3] = tauxar[idx3] + dopaer[idx1]
        wa[idx3] = wa[idx3] + dopaer[idx1] * palb[idx1]
        ga[idx3] = ga[idx3] + dopaer[idx1] * palb[idx1] * pasm[idx1]
        fa[idx3] = fa[idx3] + dopaer[idx1] * palb[idx1] * pasm[idx1] * pasm[idx1]
