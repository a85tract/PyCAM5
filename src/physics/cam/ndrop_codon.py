@inline
def _idx2(i: int, k: int, pcols: int) -> int:
    """dropmixnuc arrays declared as (pcols,pver)."""
    return (i - 1) + (k - 1) * pcols


@inline
def _idx3(i: int, k: int, m: int, pcols: int, pver: int) -> int:
    """factnum declared as (pcols,pver,ntot_amode)."""
    return (i - 1) + (k - 1) * pcols + (m - 1) * pcols * pver


@inline
def _mode_idx(k: int, m: int, pver: int) -> int:
    """nact/mact arrays declared as (pver,ntot_amode)."""
    return (k - 1) + (m - 1) * pver


@inline
def _aero_col_idx(k: int, mm: int, slot: int, pver: int, ncnst_tot: int) -> int:
    """raercol arrays declared as (pver,ncnst_tot,2)."""
    return (k - 1) + (mm - 1) * pver + (slot - 1) * pver * ncnst_tot


@inline
def _mam_idx(m: int, l: int, ntot_amode: int) -> int:
    """mam_idx is passed from mam_idx(1,0), declared as (ntot_amode,0:nspec_max)."""
    return (m - 1) + l * ntot_amode


def ndrop_dropmixnuc_zero_fields_codon(
    pcols: int,
    pver: int,
    ntot_amode: int,
    factnum_p: cobj,
    wtke_p: cobj,
):
    factnum = Ptr[float](factnum_p)
    wtke = Ptr[float](wtke_p)

    for m in range(1, ntot_amode + 1):
        for k in range(1, pver + 1):
            for i in range(1, pcols + 1):
                factnum[_idx3(i, k, m, pcols, pver)] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            wtke[_idx2(i, k, pcols)] = 0.0


def ndrop_dropmixnuc_column_init_codon(
    i: int,
    pcols: int,
    pver: int,
    top_lev: int,
    ntot_amode: int,
    gravit: float,
    rair: float,
    zkmin: float,
    zkmax: float,
    wmixmin: float,
    ncldwtr_p: cobj,
    temp_p: cobj,
    pmid_p: cobj,
    pint_p: cobj,
    rpdel_p: cobj,
    zm_p: cobj,
    kvh_p: cobj,
    wsub_p: cobj,
    qcld_p: cobj,
    qncld_p: cobj,
    srcn_p: cobj,
    cs_p: cobj,
    dz_p: cobj,
    nact_p: cobj,
    mact_p: cobj,
    zn_p: cobj,
    ekd_p: cobj,
    csbot_p: cobj,
    csbot_cscen_p: cobj,
    wtke_cen_p: cobj,
    wtke_p: cobj,
    nsource_p: cobj,
    zs_p: cobj,
):
    ncldwtr = Ptr[float](ncldwtr_p)
    temp = Ptr[float](temp_p)
    pmid = Ptr[float](pmid_p)
    pint = Ptr[float](pint_p)
    rpdel = Ptr[float](rpdel_p)
    zm = Ptr[float](zm_p)
    kvh = Ptr[float](kvh_p)
    wsub = Ptr[float](wsub_p)
    qcld = Ptr[float](qcld_p)
    qncld = Ptr[float](qncld_p)
    srcn = Ptr[float](srcn_p)
    cs = Ptr[float](cs_p)
    dz = Ptr[float](dz_p)
    nact = Ptr[float](nact_p)
    mact = Ptr[float](mact_p)
    zn = Ptr[float](zn_p)
    ekd = Ptr[float](ekd_p)
    csbot = Ptr[float](csbot_p)
    csbot_cscen = Ptr[float](csbot_cscen_p)
    wtke_cen = Ptr[float](wtke_cen_p)
    wtke = Ptr[float](wtke_p)
    nsource = Ptr[float](nsource_p)
    zs = Ptr[float](zs_p)

    for k in range(top_lev, pver):
        zs[k - 1] = 1.0 / (zm[_idx2(i, k, pcols)] - zm[_idx2(i, k + 1, pcols)])
    zs[pver - 1] = zs[pver - 2]

    for k in range(top_lev, pver + 1):
        idx = _idx2(i, k, pcols)
        qcld[k - 1] = ncldwtr[idx]
        qncld[k - 1] = 0.0
        srcn[k - 1] = 0.0
        cs[idx] = pmid[idx] / (rair * temp[idx])
        dz[idx] = 1.0 / (cs[idx] * gravit * rpdel[idx])

        for m in range(1, ntot_amode + 1):
            nact[_mode_idx(k, m, pver)] = 0.0
            mact[_mode_idx(k, m, pver)] = 0.0

        zn[k - 1] = gravit * rpdel[idx]

        if k < pver:
            ekd[k - 1] = kvh[_idx2(i, k + 1, pcols)]
            ekd[k - 1] = max(ekd[k - 1], zkmin)
            ekd[k - 1] = min(ekd[k - 1], zkmax)
            csbot[k - 1] = 2.0 * pint[_idx2(i, k + 1, pcols)] / (
                rair * (temp[idx] + temp[_idx2(i, k + 1, pcols)])
            )
            csbot_cscen[k - 1] = csbot[k - 1] / cs[idx]
        else:
            ekd[k - 1] = 0.0
            csbot[k - 1] = cs[idx]
            csbot_cscen[k - 1] = 1.0

        wtke_cen[idx] = wsub[idx]
        wtke[idx] = wsub[idx]
        wtke_cen[idx] = max(wtke_cen[idx], wmixmin)
        wtke[idx] = max(wtke[idx], wmixmin)

        nsource[idx] = 0.0


def ndrop_dropmixnuc_mix_setup_codon(
    i: int,
    pcols: int,
    pver: int,
    top_lev: int,
    ntot_amode: int,
    dtmicro: float,
    taumix_internal_pver_inv: float,
    cldn_p: cobj,
    zs_p: cobj,
    zn_p: cobj,
    csbot_p: cobj,
    ekd_p: cobj,
    nact_p: cobj,
    mact_p: cobj,
    ekk0_p: cobj,
    ekkp_p: cobj,
    ekkm_p: cobj,
    overlapp_p: cobj,
    overlapm_p: cobj,
    count_submix_p: cobj,
    nsubmix_p: cobj,
    dtmix_p: cobj,
):
    cldn = Ptr[float](cldn_p)
    zs = Ptr[float](zs_p)
    zn = Ptr[float](zn_p)
    csbot = Ptr[float](csbot_p)
    ekd = Ptr[float](ekd_p)
    nact = Ptr[float](nact_p)
    mact = Ptr[float](mact_p)
    ekk0 = Ptr[float](ekk0_p)
    ekkp = Ptr[float](ekkp_p)
    ekkm = Ptr[float](ekkm_p)
    overlapp = Ptr[float](overlapp_p)
    overlapm = Ptr[float](overlapm_p)
    count_submix = Ptr[i32](count_submix_p)
    nsubmix_out = Ptr[i32](nsubmix_p)
    dtmix_out = Ptr[float](dtmix_p)

    dtmin = dtmicro
    ekk0[top_lev - 1] = 0.0
    ekk0[pver] = 0.0
    for k in range(top_lev, pver):
        ekk0[k] = ekd[k - 1] * csbot[k - 1]

    for k in range(top_lev, pver + 1):
        km1 = max(k - 1, top_lev)
        ekkp[k - 1] = zn[k - 1] * ekk0[k] * zs[k - 1]
        ekkm[k - 1] = zn[k - 1] * ekk0[k - 1] * zs[km1 - 1]
        tinv = ekkp[k - 1] + ekkm[k - 1]
        if k == pver:
            tinv = tinv + taumix_internal_pver_inv

        if tinv > 1.0e-6:
            dtt = 1.0 / tinv
            dtmin = min(dtmin, dtt)

    dtmix = 0.9 * dtmin
    nsubmix = int(dtmicro / dtmix + 1.0)
    if nsubmix > 100:
        nsubmix_bnd = 100
    else:
        nsubmix_bnd = nsubmix
    count_submix[nsubmix_bnd - 1] = i32(int(count_submix[nsubmix_bnd - 1]) + 1)
    dtmix = dtmicro / float(nsubmix)
    nsubmix_out[0] = i32(nsubmix)
    dtmix_out[0] = dtmix

    for k in range(top_lev, pver + 1):
        kp1 = min(k + 1, pver)
        km1 = max(k - 1, top_lev)
        if cldn[_idx2(i, kp1, pcols)] > 1.0e-10:
            overlapp[k - 1] = min(cldn[_idx2(i, k, pcols)] / cldn[_idx2(i, kp1, pcols)], 1.0)
        else:
            overlapp[k - 1] = 1.0
        if cldn[_idx2(i, km1, pcols)] > 1.0e-10:
            overlapm[k - 1] = min(cldn[_idx2(i, k, pcols)] / cldn[_idx2(i, km1, pcols)], 1.0)
        else:
            overlapm[k - 1] = 1.0

    for k in range(top_lev, pver):
        for m in range(1, ntot_amode + 1):
            idx = _mode_idx(k, m, pver)
            nact[idx] = min(nact[idx], ekkp[k - 1])
            mact[idx] = min(mact[idx], ekkp[k - 1])


def ndrop_dropmixnuc_aero_column_copy_codon(
    i: int,
    pcols: int,
    pver: int,
    top_lev: int,
    ncnst_tot: int,
    mm: int,
    slot: int,
    zero_all: int,
    raer_fld_p: cobj,
    qqcw_fld_p: cobj,
    raercol_p: cobj,
    raercol_cw_p: cobj,
):
    raer_fld = Ptr[float](raer_fld_p)
    qqcw_fld = Ptr[float](qqcw_fld_p)
    raercol = Ptr[float](raercol_p)
    raercol_cw = Ptr[float](raercol_cw_p)

    if zero_all != 0:
        for k in range(1, pver + 1):
            col_idx = _aero_col_idx(k, mm, slot, pver, ncnst_tot)
            raercol_cw[col_idx] = 0.0
            raercol[col_idx] = 0.0

    for k in range(top_lev, pver + 1):
        col_idx = _aero_col_idx(k, mm, slot, pver, ncnst_tot)
        field_idx = _idx2(i, k, pcols)
        raercol_cw[col_idx] = qqcw_fld[field_idx]
        raercol[col_idx] = raer_fld[field_idx]


def ndrop_dropmixnuc_aero_tend_prepare_codon(
    i: int,
    pcols: int,
    pver: int,
    top_lev: int,
    ncnst_tot: int,
    mm: int,
    slot: int,
    dtinv: float,
    raer_fld_p: cobj,
    qqcw_fld_p: cobj,
    raercol_p: cobj,
    raercol_cw_p: cobj,
    raertend_p: cobj,
    qqcwtend_p: cobj,
):
    raer_fld = Ptr[float](raer_fld_p)
    qqcw_fld = Ptr[float](qqcw_fld_p)
    raercol = Ptr[float](raercol_p)
    raercol_cw = Ptr[float](raercol_cw_p)
    raertend = Ptr[float](raertend_p)
    qqcwtend = Ptr[float](qqcwtend_p)

    for k in range(top_lev, pver + 1):
        col_idx = _aero_col_idx(k, mm, slot, pver, ncnst_tot)
        field_idx = _idx2(i, k, pcols)
        tend_idx = k - 1
        raertend[tend_idx] = (raercol[col_idx] - raer_fld[field_idx]) * dtinv
        qqcwtend[tend_idx] = (raercol_cw[col_idx] - qqcw_fld[field_idx]) * dtinv


def ndrop_dropmixnuc_aero_tend_commit_qqcw_codon(
    i: int,
    pcols: int,
    pver: int,
    top_lev: int,
    ncnst_tot: int,
    mm: int,
    slot: int,
    qqcw_fld_p: cobj,
    raercol_cw_p: cobj,
):
    qqcw_fld = Ptr[float](qqcw_fld_p)
    raercol_cw = Ptr[float](raercol_cw_p)

    for k in range(1, pver + 1):
        field_idx = _idx2(i, k, pcols)
        qqcw_fld[field_idx] = 0.0

    for k in range(top_lev, pver + 1):
        col_idx = _aero_col_idx(k, mm, slot, pver, ncnst_tot)
        field_idx = _idx2(i, k, pcols)
        qqcw_fld[field_idx] = raercol_cw[col_idx]


def ndrop_dropmixnuc_aero_coltend_codon(
    i: int,
    pcols: int,
    pver: int,
    mm: int,
    gravit: float,
    pdel_p: cobj,
    raertend_p: cobj,
    qqcwtend_p: cobj,
    coltend_out_p: cobj,
    coltend_cw_out_p: cobj,
):
    pdel = Ptr[float](pdel_p)
    raertend = Ptr[float](raertend_p)
    qqcwtend = Ptr[float](qqcwtend_p)
    coltend_out = Ptr[float](coltend_out_p)
    coltend_cw_out = Ptr[float](coltend_cw_out_p)

    sum_raer = 0.0
    sum_cw = 0.0
    for k in range(1, pver + 1):
        pdel_ik = pdel[_idx2(i, k, pcols)]
        sum_raer = sum_raer + pdel_ik * raertend[k - 1]
        sum_cw = sum_cw + pdel_ik * qqcwtend[k - 1]

    coltend_out[0] = sum_raer / gravit
    coltend_cw_out[0] = sum_cw / gravit


def ndrop_dropmixnuc_finalize_column_codon(
    i: int,
    pcols: int,
    pver: int,
    top_lev: int,
    dtinv: float,
    gravit: float,
    qcld_p: cobj,
    ncldwtr_p: cobj,
    pdel_p: cobj,
    nsource_p: cobj,
    ndropmix_p: cobj,
    tendnd_p: cobj,
    ndropcol_p: cobj,
):
    qcld = Ptr[float](qcld_p)
    ncldwtr = Ptr[float](ncldwtr_p)
    pdel = Ptr[float](pdel_p)
    nsource = Ptr[float](nsource_p)
    ndropmix = Ptr[float](ndropmix_p)
    tendnd = Ptr[float](tendnd_p)
    ndropcol = Ptr[float](ndropcol_p)

    ndropcol[i - 1] = 0.0
    for k in range(top_lev, pver + 1):
        idx = _idx2(i, k, pcols)
        ndropmix[idx] = (qcld[k - 1] - ncldwtr[idx]) * dtinv - nsource[idx]
        tendnd[idx] = (max(qcld[k - 1], 1.0e-6) - ncldwtr[idx]) * dtinv
        ndropcol[i - 1] = ndropcol[i - 1] + ncldwtr[idx] * pdel[idx]
    ndropcol[i - 1] = ndropcol[i - 1] / gravit


def ndrop_dropmixnuc_clear_old_cloud_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    ntot_amode: int,
    ncnst_tot: int,
    nsav: int,
    dtinv: float,
    qcld_p: cobj,
    nsource_p: cobj,
    nspec_amode_p: cobj,
    mam_idx_p: cobj,
    raercol_p: cobj,
    raercol_cw_p: cobj,
):
    qcld = Ptr[float](qcld_p)
    nsource = Ptr[float](nsource_p)
    nspec_amode = Ptr[i32](nspec_amode_p)
    mam_idx = Ptr[i32](mam_idx_p)
    raercol = Ptr[float](raercol_p)
    raercol_cw = Ptr[float](raercol_cw_p)

    nsource[_idx2(i, k, pcols)] = nsource[_idx2(i, k, pcols)] - qcld[k - 1] * dtinv
    qcld[k - 1] = 0.0

    for m in range(1, ntot_amode + 1):
        mm = int(mam_idx[_mam_idx(m, 0, ntot_amode)])
        col_idx = _aero_col_idx(k, mm, nsav, pver, ncnst_tot)
        raercol[col_idx] = raercol[col_idx] + raercol_cw[col_idx]
        raercol_cw[col_idx] = 0.0

        for l in range(1, int(nspec_amode[m - 1]) + 1):
            mm = int(mam_idx[_mam_idx(m, l, ntot_amode)])
            col_idx = _aero_col_idx(k, mm, nsav, pver, ncnst_tot)
            raercol[col_idx] = raercol[col_idx] + raercol_cw[col_idx]
            raercol_cw[col_idx] = 0.0


def ndrop_dropmixnuc_factnum_store_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    ntot_amode: int,
    fn_p: cobj,
    factnum_p: cobj,
):
    fn = Ptr[float](fn_p)
    factnum = Ptr[float](factnum_p)

    for m in range(1, ntot_amode + 1):
        factnum[_idx3(i, k, m, pcols, pver)] = fn[m - 1]


def ndrop_dropmixnuc_shrink_cloud_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    ntot_amode: int,
    ncnst_tot: int,
    nsav: int,
    dtinv: float,
    cldn_tmp: float,
    cldo_tmp: float,
    qcld_p: cobj,
    nsource_p: cobj,
    nspec_amode_p: cobj,
    mam_idx_p: cobj,
    raercol_p: cobj,
    raercol_cw_p: cobj,
):
    qcld = Ptr[float](qcld_p)
    nsource = Ptr[float](nsource_p)
    nspec_amode = Ptr[i32](nspec_amode_p)
    mam_idx = Ptr[i32](mam_idx_p)
    raercol = Ptr[float](raercol_p)
    raercol_cw = Ptr[float](raercol_cw_p)

    nsource[_idx2(i, k, pcols)] = nsource[_idx2(i, k, pcols)] + (
        qcld[k - 1] * (cldn_tmp - cldo_tmp) / cldo_tmp * dtinv
    )
    qcld[k - 1] = qcld[k - 1] * (1.0 + (cldn_tmp - cldo_tmp) / cldo_tmp)

    dumc = (cldn_tmp - cldo_tmp) / cldo_tmp
    for m in range(1, ntot_amode + 1):
        mm = int(mam_idx[_mam_idx(m, 0, ntot_amode)])
        col_idx = _aero_col_idx(k, mm, nsav, pver, ncnst_tot)
        dact = raercol_cw[col_idx] * dumc
        raercol_cw[col_idx] = raercol_cw[col_idx] + dact
        raercol[col_idx] = raercol[col_idx] - dact

        for l in range(1, int(nspec_amode[m - 1]) + 1):
            mm = int(mam_idx[_mam_idx(m, l, ntot_amode)])
            col_idx = _aero_col_idx(k, mm, nsav, pver, ncnst_tot)
            dact = raercol_cw[col_idx] * dumc
            raercol_cw[col_idx] = raercol_cw[col_idx] + dact
            raercol[col_idx] = raercol[col_idx] - dact


def ndrop_dropmixnuc_grow_cloud_number_update_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    ncnst_tot: int,
    nsav: int,
    mm: int,
    dtinv: float,
    dumc: float,
    fn_m: float,
    raer_fld_p: cobj,
    qcld_p: cobj,
    nsource_p: cobj,
    raercol_p: cobj,
    raercol_cw_p: cobj,
):
    raer_fld = Ptr[float](raer_fld_p)
    qcld = Ptr[float](qcld_p)
    nsource = Ptr[float](nsource_p)
    raercol = Ptr[float](raercol_p)
    raercol_cw = Ptr[float](raercol_cw_p)

    dact = dumc * fn_m * raer_fld[_idx2(i, k, pcols)]
    qcld[k - 1] = qcld[k - 1] + dact
    nsource[_idx2(i, k, pcols)] = nsource[_idx2(i, k, pcols)] + dact * dtinv
    col_idx = _aero_col_idx(k, mm, nsav, pver, ncnst_tot)
    raercol_cw[col_idx] = raercol_cw[col_idx] + dact
    raercol[col_idx] = raercol[col_idx] - dact


def ndrop_dropmixnuc_grow_cloud_species_update_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    ncnst_tot: int,
    nsav: int,
    mm: int,
    dum: float,
    raer_fld_p: cobj,
    raercol_p: cobj,
    raercol_cw_p: cobj,
):
    raer_fld = Ptr[float](raer_fld_p)
    raercol = Ptr[float](raercol_p)
    raercol_cw = Ptr[float](raercol_cw_p)

    dact = dum * raer_fld[_idx2(i, k, pcols)]
    col_idx = _aero_col_idx(k, mm, nsav, pver, ncnst_tot)
    raercol_cw[col_idx] = raercol_cw[col_idx] + dact
    raercol[col_idx] = raercol[col_idx] - dact


def ndrop_dropmixnuc_old_cloud_activate_update_codon(
    i: int,
    k: int,
    kp1: int,
    pcols: int,
    pver: int,
    ntot_amode: int,
    ncnst_tot: int,
    nsav: int,
    dumc: float,
    dum: float,
    cs_ik: float,
    dz_ik: float,
    taumix_internal_pver_inv: float,
    fluxn_p: cobj,
    fluxm_p: cobj,
    nact_p: cobj,
    mact_p: cobj,
    mam_idx_p: cobj,
    raercol_p: cobj,
    raercol_cw_p: cobj,
    srcn_p: cobj,
    nsource_p: cobj,
):
    fluxn = Ptr[float](fluxn_p)
    fluxm = Ptr[float](fluxm_p)
    nact = Ptr[float](nact_p)
    mact = Ptr[float](mact_p)
    mam_idx = Ptr[i32](mam_idx_p)
    raercol = Ptr[float](raercol_p)
    raercol_cw = Ptr[float](raercol_cw_p)
    srcn = Ptr[float](srcn_p)
    nsource = Ptr[float](nsource_p)

    fluxntot = 0.0
    for m in range(1, ntot_amode + 1):
        mm = int(mam_idx[_mam_idx(m, 0, ntot_amode)])
        fluxn[m - 1] = fluxn[m - 1] * dumc
        fluxm[m - 1] = fluxm[m - 1] * dumc
        nact[_mode_idx(k, m, pver)] = nact[_mode_idx(k, m, pver)] + fluxn[m - 1] * dum
        mact[_mode_idx(k, m, pver)] = mact[_mode_idx(k, m, pver)] + fluxm[m - 1] * dum
        if k < pver:
            fluxntot = fluxntot + (
                fluxn[m - 1] * raercol[_aero_col_idx(kp1, mm, nsav, pver, ncnst_tot)] * cs_ik
            )
        else:
            tmpa = (
                raercol[_aero_col_idx(kp1, mm, nsav, pver, ncnst_tot)] * fluxn[m - 1]
                + raercol_cw[_aero_col_idx(kp1, mm, nsav, pver, ncnst_tot)]
                * (fluxn[m - 1] - taumix_internal_pver_inv * dz_ik)
            )
            fluxntot = fluxntot + max(0.0, tmpa) * cs_ik

    srcn[k - 1] = srcn[k - 1] + fluxntot / (cs_ik * dz_ik)
    nsource[_idx2(i, k, pcols)] = nsource[_idx2(i, k, pcols)] + fluxntot / (cs_ik * dz_ik)


def ndrop_dropmixnuc_srcn_from_nact_codon(
    pver: int,
    top_lev: int,
    ntot_amode: int,
    ncnst_tot: int,
    nsav: int,
    taumix_internal_pver_inv: float,
    nact_p: cobj,
    mam_idx_p: cobj,
    raercol_p: cobj,
    raercol_cw_p: cobj,
    srcn_p: cobj,
):
    nact = Ptr[float](nact_p)
    mam_idx = Ptr[i32](mam_idx_p)
    raercol = Ptr[float](raercol_p)
    raercol_cw = Ptr[float](raercol_cw_p)
    srcn = Ptr[float](srcn_p)

    for m in range(1, ntot_amode + 1):
        mm = int(mam_idx[_mam_idx(m, 0, ntot_amode)])
        for k in range(top_lev, pver):
            srcn[k - 1] = srcn[k - 1] + nact[_mode_idx(k, m, pver)] * (
                raercol[_aero_col_idx(k + 1, mm, nsav, pver, ncnst_tot)]
            )

        tmpa = (
            raercol[_aero_col_idx(pver, mm, nsav, pver, ncnst_tot)]
            * nact[_mode_idx(pver, m, pver)]
            + raercol_cw[_aero_col_idx(pver, mm, nsav, pver, ncnst_tot)]
            * (nact[_mode_idx(pver, m, pver)] - taumix_internal_pver_inv)
        )
        srcn[pver - 1] = srcn[pver - 1] + max(0.0, tmpa)


def ndrop_dropmixnuc_source_from_act_codon(
    pver: int,
    top_lev: int,
    ncnst_tot: int,
    m: int,
    mm: int,
    nsav: int,
    taumix_internal_pver_inv: float,
    act_p: cobj,
    raercol_p: cobj,
    raercol_cw_p: cobj,
    source_p: cobj,
):
    act = Ptr[float](act_p)
    raercol = Ptr[float](raercol_p)
    raercol_cw = Ptr[float](raercol_cw_p)
    source = Ptr[float](source_p)

    for k in range(top_lev, pver):
        source[k - 1] = act[_mode_idx(k, m, pver)] * (
            raercol[_aero_col_idx(k + 1, mm, nsav, pver, ncnst_tot)]
        )

    tmpa = (
        raercol[_aero_col_idx(pver, mm, nsav, pver, ncnst_tot)]
        * act[_mode_idx(pver, m, pver)]
        + raercol_cw[_aero_col_idx(pver, mm, nsav, pver, ncnst_tot)]
        * (act[_mode_idx(pver, m, pver)] - taumix_internal_pver_inv)
    )
    source[pver - 1] = max(0.0, tmpa)


def ndrop_dropmixnuc_evaporate_clear_layers_codon(
    i: int,
    pcols: int,
    pver: int,
    top_lev: int,
    ntot_amode: int,
    ncnst_tot: int,
    nnew: int,
    cldn_p: cobj,
    qcld_p: cobj,
    nspec_amode_p: cobj,
    mam_idx_p: cobj,
    raercol_p: cobj,
    raercol_cw_p: cobj,
):
    cldn = Ptr[float](cldn_p)
    qcld = Ptr[float](qcld_p)
    nspec_amode = Ptr[i32](nspec_amode_p)
    mam_idx = Ptr[i32](mam_idx_p)
    raercol = Ptr[float](raercol_p)
    raercol_cw = Ptr[float](raercol_cw_p)

    for k in range(top_lev, pver + 1):
        if cldn[_idx2(i, k, pcols)] == 0.0:
            qcld[k - 1] = 0.0
            for m in range(1, ntot_amode + 1):
                mm = int(mam_idx[_mam_idx(m, 0, ntot_amode)])
                col_idx = _aero_col_idx(k, mm, nnew, pver, ncnst_tot)
                raercol[col_idx] = raercol[col_idx] + raercol_cw[col_idx]
                raercol_cw[col_idx] = 0.0

                for l in range(1, int(nspec_amode[m - 1]) + 1):
                    mm = int(mam_idx[_mam_idx(m, l, ntot_amode)])
                    col_idx = _aero_col_idx(k, mm, nnew, pver, ncnst_tot)
                    raercol[col_idx] = raercol[col_idx] + raercol_cw[col_idx]
                    raercol_cw[col_idx] = 0.0


def ndrop_dropmixnuc_swap_slots_codon(nsav_p: cobj, nnew_p: cobj):
    nsav = Ptr[i32](nsav_p)
    nnew = Ptr[i32](nnew_p)

    ntemp = nsav[0]
    nsav[0] = nnew[0]
    nnew[0] = ntemp


def ndrop_dropmixnuc_submix_iter_init_codon(
    pver: int,
    qcld_p: cobj,
    qncld_p: cobj,
    srcn_p: cobj,
    nsav_p: cobj,
    nnew_p: cobj,
):
    qcld = Ptr[float](qcld_p)
    qncld = Ptr[float](qncld_p)
    srcn = Ptr[float](srcn_p)

    for k in range(1, pver + 1):
        qncld[k - 1] = qcld[k - 1]
        srcn[k - 1] = 0.0

    ndrop_dropmixnuc_swap_slots_codon(nsav_p, nnew_p)


def ndrop_dropmixnuc_zero_tendencies_codon(
    pver: int,
    raertend_p: cobj,
    qqcwtend_p: cobj,
):
    raertend = Ptr[float](raertend_p)
    qqcwtend = Ptr[float](qqcwtend_p)

    for k in range(1, pver + 1):
        raertend[k - 1] = 0.0
        qqcwtend[k - 1] = 0.0


def ndrop_explmix_codon(
    pver: int,
    top_lev: int,
    surfrate: float,
    flxconv: float,
    dt: float,
    is_unact: int,
    q_p: cobj,
    src_p: cobj,
    ekkp_p: cobj,
    ekkm_p: cobj,
    overlapp_p: cobj,
    overlapm_p: cobj,
    qold_p: cobj,
    qactold_p: cobj,
):
    q = Ptr[float](q_p)
    src = Ptr[float](src_p)
    ekkp = Ptr[float](ekkp_p)
    ekkm = Ptr[float](ekkm_p)
    overlapp = Ptr[float](overlapp_p)
    overlapm = Ptr[float](overlapm_p)
    qold = Ptr[float](qold_p)

    if is_unact != 0:
        qactold = Ptr[float](qactold_p)
        for k in range(top_lev, pver + 1):
            kp1 = min(k + 1, pver)
            km1 = max(k - 1, top_lev)
            k0 = k - 1
            kp10 = kp1 - 1
            km10 = km1 - 1
            q[k0] = qold[k0] + dt * (
                -src[k0]
                + ekkp[k0]
                * (qold[kp10] - qold[k0] + qactold[kp10] * (1.0 - overlapp[k0]))
                + ekkm[k0]
                * (qold[km10] - qold[k0] + qactold[km10] * (1.0 - overlapm[k0]))
            )
            q[k0] = max(q[k0], 0.0)

        q[pver - 1] = q[pver - 1] - surfrate * qold[pver - 1] * dt + flxconv * dt
        q[pver - 1] = max(q[pver - 1], 0.0)
    else:
        for k in range(top_lev, pver + 1):
            kp1 = min(k + 1, pver)
            km1 = max(k - 1, top_lev)
            k0 = k - 1
            kp10 = kp1 - 1
            km10 = km1 - 1
            q[k0] = qold[k0] + dt * (
                src[k0]
                + ekkp[k0] * (overlapp[k0] * qold[kp10] - qold[k0])
                + ekkm[k0] * (overlapm[k0] * qold[km10] - qold[k0])
            )
            q[k0] = max(q[k0], 0.0)

        q[pver - 1] = q[pver - 1] - surfrate * qold[pver - 1] * dt + flxconv * dt
        q[pver - 1] = max(q[pver - 1], 0.0)
