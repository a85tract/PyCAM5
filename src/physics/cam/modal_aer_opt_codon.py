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
