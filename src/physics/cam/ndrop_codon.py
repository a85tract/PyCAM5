from math import erf, exp, log, sqrt


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


def ndrop_init_scalars_codon(
    mwh2o: float,
    r_universal: float,
    rhoh2o: float,
    pi: float,
    scalars_p: cobj,
):
    scalars = Ptr[float](scalars_p)

    zero = 0.0
    third = 1.0 / 3.0
    twothird = 2.0 * third
    sixth = 1.0 / 6.0
    sq2 = sqrt(2.0)
    sqpi = sqrt(pi)

    t0 = 273.0
    surften = 0.076
    aten = 2.0 * mwh2o * surften / (r_universal * t0 * rhoh2o)
    alogaten = log(aten)
    alog2 = log(2.0)
    alog3 = log(3.0)

    scalars[0] = zero
    scalars[1] = third
    scalars[2] = twothird
    scalars[3] = sixth
    scalars[4] = sq2
    scalars[5] = sqpi
    scalars[6] = t0
    scalars[7] = surften
    scalars[8] = aten
    scalars[9] = alogaten
    scalars[10] = alog2
    scalars[11] = alog3


def ndrop_init_counts_codon(
    nmode: int,
    nspec_amode_p: cobj,
    nspec_max_p: cobj,
    ncnst_tot_p: cobj,
):
    nspec_amode = Ptr[i32](nspec_amode_p)
    nspec_max_out = Ptr[int](nspec_max_p)
    ncnst_tot_out = Ptr[int](ncnst_tot_p)

    nspec_max = int(nspec_amode[0])
    ncnst_tot = int(nspec_amode[0]) + 1
    for m in range(2, nmode + 1):
        nspec_m = int(nspec_amode[m - 1])
        nspec_max = max(nspec_max, nspec_m)
        ncnst_tot = ncnst_tot + nspec_m + 1

    nspec_max_out[0] = nspec_max
    ncnst_tot_out[0] = ncnst_tot


def ndrop_init_mam_idx_codon(
    nmode: int,
    nspec_max: int,
    nspec_amode_p: cobj,
    mam_idx_p: cobj,
):
    nspec_amode = Ptr[i32](nspec_amode_p)
    mam_idx = Ptr[i32](mam_idx_p)

    ii = 0
    for m in range(1, nmode + 1):
        for l in range(0, int(nspec_amode[m - 1]) + 1):
            ii += 1
            mam_idx[(m - 1) + l * nmode] = i32(ii)


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


def ndrop_dropmixnuc_aero_column_copy_all_codon(
    i: int,
    pcols: int,
    pver: int,
    top_lev: int,
    ntot_amode: int,
    ncnst_tot: int,
    slot: int,
    raer_ptrs_p: cobj,
    qqcw_ptrs_p: cobj,
    nspec_amode_p: cobj,
    mam_idx_p: cobj,
    raercol_p: cobj,
    raercol_cw_p: cobj,
):
    raer_ptrs = Ptr[cobj](raer_ptrs_p)
    qqcw_ptrs = Ptr[cobj](qqcw_ptrs_p)
    nspec_amode = Ptr[i32](nspec_amode_p)
    mam_idx = Ptr[i32](mam_idx_p)
    raercol = Ptr[float](raercol_p)
    raercol_cw = Ptr[float](raercol_cw_p)

    for m in range(1, ntot_amode + 1):
        mm = int(mam_idx[_mam_idx(m, 0, ntot_amode)])
        raer_fld = Ptr[float](raer_ptrs[mm - 1])
        qqcw_fld = Ptr[float](qqcw_ptrs[mm - 1])

        for k in range(1, pver + 1):
            col_idx = _aero_col_idx(k, mm, slot, pver, ncnst_tot)
            raercol_cw[col_idx] = 0.0
            raercol[col_idx] = 0.0

        for k in range(top_lev, pver + 1):
            col_idx = _aero_col_idx(k, mm, slot, pver, ncnst_tot)
            field_idx = _idx2(i, k, pcols)
            raercol_cw[col_idx] = qqcw_fld[field_idx]
            raercol[col_idx] = raer_fld[field_idx]

        for l in range(1, int(nspec_amode[m - 1]) + 1):
            mm = int(mam_idx[_mam_idx(m, l, ntot_amode)])
            raer_fld = Ptr[float](raer_ptrs[mm - 1])
            qqcw_fld = Ptr[float](qqcw_ptrs[mm - 1])
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


def ndrop_dropmixnuc_aero_tend_commit_ptend_codon(
    i: int,
    psetcols: int,
    pver: int,
    top_lev: int,
    pcnst: int,
    lptr: int,
    raertend_p: cobj,
    ptend_q_p: cobj,
):
    raertend = Ptr[float](raertend_p)
    ptend_q = Ptr[float](ptend_q_p)

    for k in range(1, pver + 1):
        ptend_q[_idx3(i, k, lptr, psetcols, pver)] = 0.0

    for k in range(top_lev, pver + 1):
        ptend_q[_idx3(i, k, lptr, psetcols, pver)] = raertend[k - 1]


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


def ndrop_dropmixnuc_aero_tend_all_codon(
    i: int,
    pcols: int,
    psetcols: int,
    pver: int,
    top_lev: int,
    ntot_amode: int,
    ncnst_tot: int,
    slot: int,
    dtinv: float,
    gravit: float,
    raer_ptrs_p: cobj,
    qqcw_ptrs_p: cobj,
    nspec_amode_p: cobj,
    mam_idx_p: cobj,
    mam_cnst_idx_p: cobj,
    pdel_p: cobj,
    raercol_p: cobj,
    raercol_cw_p: cobj,
    coltend_p: cobj,
    coltend_cw_p: cobj,
    ptend_q_p: cobj,
):
    raer_ptrs = Ptr[cobj](raer_ptrs_p)
    qqcw_ptrs = Ptr[cobj](qqcw_ptrs_p)
    nspec_amode = Ptr[i32](nspec_amode_p)
    mam_idx = Ptr[i32](mam_idx_p)
    mam_cnst_idx = Ptr[i32](mam_cnst_idx_p)
    pdel = Ptr[float](pdel_p)
    raercol = Ptr[float](raercol_p)
    raercol_cw = Ptr[float](raercol_cw_p)
    coltend = Ptr[float](coltend_p)
    coltend_cw = Ptr[float](coltend_cw_p)
    ptend_q = Ptr[float](ptend_q_p)

    for m in range(1, ntot_amode + 1):
        for l in range(0, int(nspec_amode[m - 1]) + 1):
            mm = int(mam_idx[_mam_idx(m, l, ntot_amode)])
            lptr = int(mam_cnst_idx[_mam_idx(m, l, ntot_amode)])
            raer_fld = Ptr[float](raer_ptrs[mm - 1])
            qqcw_fld = Ptr[float](qqcw_ptrs[mm - 1])

            for k in range(1, pver + 1):
                ptend_q[_idx3(i, k, lptr, psetcols, pver)] = 0.0

            sum_raer = 0.0
            sum_cw = 0.0
            for k in range(1, pver + 1):
                if k >= top_lev:
                    col_idx = _aero_col_idx(k, mm, slot, pver, ncnst_tot)
                    field_idx = _idx2(i, k, pcols)
                    raertend = (raercol[col_idx] - raer_fld[field_idx]) * dtinv
                    qqcwtend = (raercol_cw[col_idx] - qqcw_fld[field_idx]) * dtinv
                    ptend_q[_idx3(i, k, lptr, psetcols, pver)] = raertend
                else:
                    raertend = 0.0
                    qqcwtend = 0.0

                pdel_ik = pdel[_idx2(i, k, pcols)]
                sum_raer = sum_raer + pdel_ik * raertend
                sum_cw = sum_cw + pdel_ik * qqcwtend

            coltend[_idx2(i, mm, pcols)] = sum_raer / gravit
            coltend_cw[_idx2(i, mm, pcols)] = sum_cw / gravit

            for k in range(1, pver + 1):
                qqcw_fld[_idx2(i, k, pcols)] = 0.0

            for k in range(top_lev, pver + 1):
                qqcw_fld[_idx2(i, k, pcols)] = raercol_cw[_aero_col_idx(k, mm, slot, pver, ncnst_tot)]


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


def ndrop_dropmixnuc_grow_cloud_update_all_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    ntot_amode: int,
    ncnst_tot: int,
    nsav: int,
    dtinv: float,
    dumc: float,
    raer_ptrs_p: cobj,
    nspec_amode_p: cobj,
    mam_idx_p: cobj,
    fn_p: cobj,
    fm_p: cobj,
    qcld_p: cobj,
    nsource_p: cobj,
    raercol_p: cobj,
    raercol_cw_p: cobj,
    factnum_p: cobj,
):
    raer_ptrs = Ptr[cobj](raer_ptrs_p)
    nspec_amode = Ptr[i32](nspec_amode_p)
    mam_idx = Ptr[i32](mam_idx_p)
    fn = Ptr[float](fn_p)
    fm = Ptr[float](fm_p)
    qcld = Ptr[float](qcld_p)
    nsource = Ptr[float](nsource_p)
    raercol = Ptr[float](raercol_p)
    raercol_cw = Ptr[float](raercol_cw_p)
    factnum = Ptr[float](factnum_p)

    for m in range(1, ntot_amode + 1):
        factnum[_idx3(i, k, m, pcols, pver)] = fn[m - 1]

    for m in range(1, ntot_amode + 1):
        mm = int(mam_idx[_mam_idx(m, 0, ntot_amode)])
        raer_fld = Ptr[float](raer_ptrs[mm - 1])
        dact = dumc * fn[m - 1] * raer_fld[_idx2(i, k, pcols)]
        qcld[k - 1] = qcld[k - 1] + dact
        nsource[_idx2(i, k, pcols)] = nsource[_idx2(i, k, pcols)] + dact * dtinv
        col_idx = _aero_col_idx(k, mm, nsav, pver, ncnst_tot)
        raercol_cw[col_idx] = raercol_cw[col_idx] + dact
        raercol[col_idx] = raercol[col_idx] - dact

        dum = dumc * fm[m - 1]
        for l in range(1, int(nspec_amode[m - 1]) + 1):
            mm = int(mam_idx[_mam_idx(m, l, ntot_amode)])
            raer_fld = Ptr[float](raer_ptrs[mm - 1])
            dact = dum * raer_fld[_idx2(i, k, pcols)]
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


def ndrop_loadaer_zero_codon(
    istart: int,
    istop: int,
    vaerosol_p: cobj,
    hygro_p: cobj,
):
    vaerosol = Ptr[float](vaerosol_p)
    hygro = Ptr[float](hygro_p)

    for i in range(istart, istop + 1):
        vaerosol[i - 1] = 0.0
        hygro[i - 1] = 0.0


def ndrop_loadaer_species_accum_codon(
    istart: int,
    istop: int,
    k: int,
    pcols: int,
    phase: int,
    specdens: float,
    spechygro: float,
    raer_p: cobj,
    qqcw_p: cobj,
    vaerosol_p: cobj,
    hygro_p: cobj,
):
    raer = Ptr[float](raer_p)
    qqcw = Ptr[float](qqcw_p)
    vaerosol = Ptr[float](vaerosol_p)
    hygro = Ptr[float](hygro_p)

    for i in range(istart, istop + 1):
        idx1 = i - 1
        idx2 = _idx2(i, k, pcols)
        if phase == 3:
            vol = max(raer[idx2] + qqcw[idx2], 0.0) / specdens
        elif phase == 2:
            vol = max(qqcw[idx2], 0.0) / specdens
        else:
            vol = max(raer[idx2], 0.0) / specdens
        vaerosol[idx1] = vaerosol[idx1] + vol
        hygro[idx1] = hygro[idx1] + vol * spechygro


def ndrop_loadaer_species_batch_codon(
    istart: int,
    istop: int,
    k: int,
    pcols: int,
    nspec: int,
    phase: int,
    raer_ptrs_p: cobj,
    qqcw_ptrs_p: cobj,
    specdens_p: cobj,
    spechygro_p: cobj,
    vaerosol_p: cobj,
    hygro_p: cobj,
):
    raer_ptrs = Ptr[cobj](raer_ptrs_p)
    qqcw_ptrs = Ptr[cobj](qqcw_ptrs_p)
    specdens = Ptr[float](specdens_p)
    spechygro = Ptr[float](spechygro_p)
    vaerosol = Ptr[float](vaerosol_p)
    hygro = Ptr[float](hygro_p)

    for l in range(1, nspec + 1):
        raer = Ptr[float](raer_ptrs[l - 1])
        qqcw = Ptr[float](qqcw_ptrs[l - 1])
        specdens_l = specdens[l - 1]
        spechygro_l = spechygro[l - 1]
        for i in range(istart, istop + 1):
            idx1 = i - 1
            idx2 = _idx2(i, k, pcols)
            if phase == 3:
                vol = max(raer[idx2] + qqcw[idx2], 0.0) / specdens_l
            elif phase == 2:
                vol = max(qqcw[idx2], 0.0) / specdens_l
            else:
                vol = max(raer[idx2], 0.0) / specdens_l
            vaerosol[idx1] = vaerosol[idx1] + vol
            hygro[idx1] = hygro[idx1] + vol * spechygro_l


def ndrop_loadaer_finalize_volume_codon(
    istart: int,
    istop: int,
    k: int,
    pcols: int,
    cs_p: cobj,
    vaerosol_p: cobj,
    hygro_p: cobj,
):
    cs = Ptr[float](cs_p)
    vaerosol = Ptr[float](vaerosol_p)
    hygro = Ptr[float](hygro_p)

    for i in range(istart, istop + 1):
        idx1 = i - 1
        if vaerosol[idx1] > 1.0e-30:
            hygro[idx1] = hygro[idx1] / vaerosol[idx1]
            vaerosol[idx1] = vaerosol[idx1] * cs[_idx2(i, k, pcols)]
        else:
            hygro[idx1] = 0.0
            vaerosol[idx1] = 0.0


def ndrop_loadaer_number_codon(
    istart: int,
    istop: int,
    k: int,
    pcols: int,
    phase: int,
    voltonumblo: float,
    voltonumbhi: float,
    raer_p: cobj,
    qqcw_p: cobj,
    cs_p: cobj,
    vaerosol_p: cobj,
    naerosol_p: cobj,
):
    raer = Ptr[float](raer_p)
    qqcw = Ptr[float](qqcw_p)
    cs = Ptr[float](cs_p)
    vaerosol = Ptr[float](vaerosol_p)
    naerosol = Ptr[float](naerosol_p)

    for i in range(istart, istop + 1):
        idx1 = i - 1
        idx2 = _idx2(i, k, pcols)
        if phase == 3:
            naerosol[idx1] = (raer[idx2] + qqcw[idx2]) * cs[idx2]
        elif phase == 2:
            naerosol[idx1] = qqcw[idx2] * cs[idx2]
        else:
            naerosol[idx1] = raer[idx2] * cs[idx2]
        naerosol[idx1] = max(naerosol[idx1], vaerosol[idx1] * voltonumbhi)
        naerosol[idx1] = min(naerosol[idx1], vaerosol[idx1] * voltonumblo)


def ndrop_ccncalc_zero_codon(
    pcols: int,
    pver: int,
    psat: int,
    ccn_p: cobj,
):
    ccn = Ptr[float](ccn_p)

    for l in range(1, psat + 1):
        for k in range(1, pver + 1):
            for i in range(1, pcols + 1):
                ccn[_idx3(i, k, l, pcols, pver)] = 0.0


def ndrop_ccncalc_level_coeffs_codon(
    ncol: int,
    k: int,
    pcols: int,
    surften_coef: float,
    smcoefcoef: float,
    tair_p: cobj,
    smcoef_p: cobj,
):
    tair = Ptr[float](tair_p)
    smcoef = Ptr[float](smcoef_p)

    for i in range(1, ncol + 1):
        a = surften_coef / tair[_idx2(i, k, pcols)]
        smcoef[i - 1] = smcoefcoef * a * sqrt(a)


def ndrop_ccncalc_mode_accum_codon(
    ncol: int,
    k: int,
    pcols: int,
    pver: int,
    psat: int,
    amcubecoef_m: float,
    argfactor_m: float,
    naerosol_p: cobj,
    vaerosol_p: cobj,
    hygro_p: cobj,
    smcoef_p: cobj,
    super_p: cobj,
    amcube_p: cobj,
    sm_p: cobj,
    arg_p: cobj,
    ccn_p: cobj,
):
    naerosol = Ptr[float](naerosol_p)
    vaerosol = Ptr[float](vaerosol_p)
    hygro = Ptr[float](hygro_p)
    smcoef = Ptr[float](smcoef_p)
    super_sat = Ptr[float](super_p)
    amcube = Ptr[float](amcube_p)
    sm = Ptr[float](sm_p)
    arg = Ptr[float](arg_p)
    ccn = Ptr[float](ccn_p)

    for i in range(1, ncol + 1):
        idx1 = i - 1
        if naerosol[idx1] > 1.0e-3:
            amcube[idx1] = amcubecoef_m * vaerosol[idx1] / naerosol[idx1]
            sm[idx1] = smcoef[idx1] / sqrt(hygro[idx1] * amcube[idx1])
        else:
            sm[idx1] = 1.0

    for l in range(1, psat + 1):
        for i in range(1, ncol + 1):
            idx1 = i - 1
            arg[idx1] = argfactor_m * log(sm[idx1] / super_sat[l - 1])
            ccn[_idx3(i, k, l, pcols, pver)] = ccn[_idx3(i, k, l, pcols, pver)] + (
                naerosol[idx1] * 0.5 * (1.0 - erf(arg[idx1]))
            )


def ndrop_ccncalc_scale_codon(
    ncol: int,
    pcols: int,
    pver: int,
    psat: int,
    ccn_p: cobj,
):
    ccn = Ptr[float](ccn_p)

    for l in range(1, psat + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx3(i, k, l, pcols, pver)
                ccn[idx] = ccn[idx] * 1.0e-6


@inline
def _ndrop_maxsat_codon(
    nmode: int,
    zeta: Ptr[float],
    eta: Ptr[float],
    smc: Ptr[float],
    f1: Ptr[float],
    f2: Ptr[float],
) -> float:
    smax = 1.0e-20
    any_active = False
    for m in range(1, nmode + 1):
        m0 = m - 1
        if zeta[m0] > 1.0e5 * eta[m0] or smc[m0] * smc[m0] > 1.0e5 * eta[m0]:
            smax = 1.0e-20
        else:
            any_active = True
            break

    if not any_active:
        return smax

    total = 0.0
    for m in range(1, nmode + 1):
        m0 = m - 1
        if eta[m0] > 1.0e-20:
            g1 = zeta[m0] / eta[m0]
            g1sqrt = sqrt(g1)
            g1 = g1sqrt * g1
            g2 = smc[m0] / sqrt(eta[m0] + 3.0 * zeta[m0])
            g2sqrt = sqrt(g2)
            g2 = g2sqrt * g2
            total = total + (f1[m0] * g1 + f2[m0] * g2) / (smc[m0] * smc[m0])
        else:
            total = 1.0e20

    return 1.0 / sqrt(total)


def ndrop_maxsat_codon(
    nmode: int,
    zeta_p: cobj,
    eta_p: cobj,
    smc_p: cobj,
    f1_p: cobj,
    f2_p: cobj,
) -> float:
    zeta = Ptr[float](zeta_p)
    eta = Ptr[float](eta_p)
    smc = Ptr[float](smc_p)
    f1 = Ptr[float](f1_p)
    f2 = Ptr[float](f2_p)
    return _ndrop_maxsat_codon(nmode, zeta, eta, smc, f1, f2)


def ndrop_activate_modal_core_codon(
    wbar: float,
    sigw: float,
    wdiab: float,
    wminf: float,
    wmaxf: float,
    tair: float,
    rhoair: float,
    qs: float,
    nmode: int,
    rair: float,
    p0: float,
    t0: float,
    rhoh2o: float,
    latvap: float,
    cpair: float,
    rh2o: float,
    gravit: float,
    pi_value: float,
    aten: float,
    twothird: float,
    sq2: float,
    sqpi: float,
    sixth: float,
    zero_value: float,
    na_p: cobj,
    volume_p: cobj,
    hygro_p: cobj,
    alogsig_p: cobj,
    exp45logsig_p: cobj,
    f1_p: cobj,
    f2_p: cobj,
    fn_p: cobj,
    fm_p: cobj,
    fluxn_p: cobj,
    fluxm_p: cobj,
    flux_fullact_p: cobj,
    zeta_p: cobj,
    eta_p: cobj,
    etafactor2_p: cobj,
    sqrtg_p: cobj,
    amcube_p: cobj,
    smc_p: cobj,
    lnsm_p: cobj,
    sumflxn_p: cobj,
    sumflxm_p: cobj,
    sumfn_p: cobj,
    sumfm_p: cobj,
    fnold_p: cobj,
    fmold_p: cobj,
) -> int:
    na = Ptr[float](na_p)
    volume = Ptr[float](volume_p)
    hygro = Ptr[float](hygro_p)
    alogsig = Ptr[float](alogsig_p)
    exp45logsig = Ptr[float](exp45logsig_p)
    f1 = Ptr[float](f1_p)
    f2 = Ptr[float](f2_p)
    fn = Ptr[float](fn_p)
    fm = Ptr[float](fm_p)
    fluxn = Ptr[float](fluxn_p)
    fluxm = Ptr[float](fluxm_p)
    flux_fullact = Ptr[float](flux_fullact_p)
    zeta = Ptr[float](zeta_p)
    eta = Ptr[float](eta_p)
    etafactor2 = Ptr[float](etafactor2_p)
    sqrtg = Ptr[float](sqrtg_p)
    amcube = Ptr[float](amcube_p)
    smc = Ptr[float](smc_p)
    lnsm = Ptr[float](lnsm_p)
    sumflxn = Ptr[float](sumflxn_p)
    sumflxm = Ptr[float](sumflxm_p)
    sumfn = Ptr[float](sumfn_p)
    sumfm = Ptr[float](sumfm_p)
    fnold = Ptr[float](fnold_p)
    fmold = Ptr[float](fmold_p)

    for m in range(1, nmode + 1):
        m0 = m - 1
        fn[m0] = 0.0
        fm[m0] = 0.0
        fluxn[m0] = 0.0
        fluxm[m0] = 0.0
    flux_fullact[0] = 0.0

    if nmode == 1 and na[0] < 1.0e-20:
        return 0
    if sigw <= 1.0e-5 and wbar <= 0.0:
        return 0

    pres = rair * rhoair * tair
    p0_over_pres = p0 / pres
    tair_over_t0 = tair / t0
    diff0 = 0.211e-4 * p0_over_pres * (tair_over_t0 ** 1.94)
    conduct0 = (5.69 + 0.017 * (tair - t0)) * 4.186e2 * 1.0e-5
    dqsdt = latvap / (rh2o * tair * tair) * qs
    alpha = gravit * (latvap / (cpair * rh2o * tair * tair) - 1.0 / (rair * tair))
    gamma = (1.0 + latvap / cpair * dqsdt) / (rhoair * qs)
    etafactor2max = 1.0e10 / ((alpha * wmaxf) ** 1.5)

    for m in range(1, nmode + 1):
        m0 = m - 1
        if volume[m0] > 1.0e-39 and na[m0] > 1.0e-39:
            amcube[m0] = 3.0 * volume[m0] / (4.0 * pi_value * exp45logsig[m0] * na[m0])
            growth = 1.0 / (
                rhoh2o / (diff0 * rhoair * qs)
                + latvap * rhoh2o / (conduct0 * tair) * (latvap / (rh2o * tair) - 1.0)
            )
            sqrtg[m0] = sqrt(growth)
            beta = 2.0 * pi_value * rhoh2o * growth * gamma
            etafactor2[m0] = 1.0 / (na[m0] * beta * sqrtg[m0])
            if hygro[m0] > 1.0e-10:
                smc[m0] = 2.0 * aten * sqrt(aten / (27.0 * hygro[m0] * amcube[m0]))
            else:
                smc[m0] = 100.0
        else:
            growth = 1.0 / (
                rhoh2o / (diff0 * rhoair * qs)
                + latvap * rhoh2o / (conduct0 * tair) * (latvap / (rh2o * tair) - 1.0)
            )
            sqrtg[m0] = sqrt(growth)
            smc[m0] = 1.0
            etafactor2[m0] = etafactor2max
        lnsm[m0] = log(smc[m0])

    eps = 0.3
    fmax = 0.99
    sds = 3.0

    if sigw > 1.0e-5:
        wmax = min(wmaxf, wbar + sds * sigw)
        wmin = max(wminf, -wdiab)
        wmin = max(wmin, wbar - sds * sigw)
        w = wmin
        dwmax = eps * sigw
        dw = dwmax
        dfmax = 0.2
        dfmin = 0.1
        if wmax <= w:
            return 0

        for m in range(1, nmode + 1):
            m0 = m - 1
            sumflxn[m0] = 0.0
            sumfn[m0] = 0.0
            fnold[m0] = 0.0
            sumflxm[m0] = 0.0
            sumfm[m0] = 0.0
            fmold[m0] = 0.0
        sumflx_fullact = 0.0

        fold = 0.0
        wold = 0.0
        gold = 0.0
        dwmin = min(dwmax, 0.01)

        n = 1
        reached_exit = False
        while n <= 200:
            while True:
                wnuc = w + wdiab
                alw = alpha * wnuc
                sqrtalw = sqrt(alw)
                etafactor1 = alw * sqrtalw

                for m in range(1, nmode + 1):
                    m0 = m - 1
                    eta[m0] = etafactor1 * etafactor2[m0]
                    zeta[m0] = twothird * sqrtalw * aten / sqrtg[m0]

                smax = _ndrop_maxsat_codon(nmode, zeta, eta, smc, f1, f2)
                lnsmax = log(smax)
                x = twothird * (lnsm[nmode - 1] - lnsmax) / (sq2 * alogsig[nmode - 1])
                fnew = 0.5 * (1.0 - erf(x))

                dwnew = dw
                if fnew - fold > dfmax and n > 1:
                    if dw > 1.01 * dwmin:
                        dw = 0.7 * dw
                        dw = max(dw, dwmin)
                        w = wold + dw
                        continue
                    else:
                        dwnew = dwmin
                break

            fold = fnew
            z = (w - wbar) / (sigw * sq2)
            g = exp(-z * z)
            fnmin = 1.0
            for m in range(1, nmode + 1):
                m0 = m - 1
                x = twothird * (lnsm[m0] - lnsmax) / (sq2 * alogsig[m0])
                fn[m0] = 0.5 * (1.0 - erf(x))
                fnmin = min(fn[m0], fnmin)
                fnbar = fn[m0] * g + fnold[m0] * gold
                arg = x - 1.5 * sq2 * alogsig[m0]
                fm[m0] = 0.5 * (1.0 - erf(arg))
                fmbar = fm[m0] * g + fmold[m0] * gold
                wb = w + wold
                if w > 0.0:
                    sumflxn[m0] = sumflxn[m0] + sixth * (
                        wb * fnbar + (fn[m0] * g * w + fnold[m0] * gold * wold)
                    ) * dw
                    sumflxm[m0] = sumflxm[m0] + sixth * (
                        wb * fmbar + (fm[m0] * g * w + fmold[m0] * gold * wold)
                    ) * dw
                sumfn[m0] = sumfn[m0] + 0.5 * fnbar * dw
                fnold[m0] = fn[m0]
                sumfm[m0] = sumfm[m0] + 0.5 * fmbar * dw
                fmold[m0] = fm[m0]

            sumflx_fullact = sumflx_fullact + sixth * (
                wb * (g + gold) + (g * w + gold * wold)
            ) * dw
            gold = g
            wold = w
            dw = dwnew
            if n > 1 and (w > wmax or fnmin > fmax):
                reached_exit = True
                break
            w = w + dw
            n += 1

        if not reached_exit:
            return 1

        if w < wmaxf:
            z1 = (w - wbar) / (sigw * sq2)
            z2 = (wmaxf - wbar) / (sigw * sq2)
            g = exp(-z1 * z1)
            integ = sigw * 0.5 * sq2 * sqpi * (erf(z2) - erf(z1))
            wf1 = max(w, zero_value)
            zf1 = (wf1 - wbar) / (sigw * sq2)
            gf1 = exp(-zf1 * zf1)
            wf2 = max(wmaxf, zero_value)
            zf2 = (wf2 - wbar) / (sigw * sq2)
            gf2 = exp(-zf2 * zf2)
            gf = gf1 - gf2
            integf = wbar * sigw * 0.5 * sq2 * sqpi * (erf(zf2) - erf(zf1)) + sigw * sigw * gf

            for m in range(1, nmode + 1):
                m0 = m - 1
                sumflxn[m0] = sumflxn[m0] + integf * fn[m0]
                sumfn[m0] = sumfn[m0] + fn[m0] * integ
                sumflxm[m0] = sumflxm[m0] + integf * fm[m0]
                sumfm[m0] = sumfm[m0] + fm[m0] * integ
            sumflx_fullact = sumflx_fullact + integf

        norm = sq2 * sqpi * sigw
        for m in range(1, nmode + 1):
            m0 = m - 1
            fn[m0] = sumfn[m0] / norm
            if fn[m0] > 1.01:
                return 2
            fluxn[m0] = sumflxn[m0] / norm
            fm[m0] = sumfm[m0] / norm
            fluxm[m0] = sumflxm[m0] / norm
        flux_fullact[0] = sumflx_fullact / norm
    else:
        wnuc = wbar + wdiab
        if wnuc > 0.0:
            w = wbar
            alw = alpha * wnuc
            sqrtalw = sqrt(alw)
            etafactor1 = alw * sqrtalw

            for m in range(1, nmode + 1):
                m0 = m - 1
                eta[m0] = etafactor1 * etafactor2[m0]
                zeta[m0] = twothird * sqrtalw * aten / sqrtg[m0]

            smax = _ndrop_maxsat_codon(nmode, zeta, eta, smc, f1, f2)
            lnsmax = log(smax)

            for m in range(1, nmode + 1):
                m0 = m - 1
                x = twothird * (lnsm[m0] - lnsmax) / (sq2 * alogsig[m0])
                fn[m0] = 0.5 * (1.0 - erf(x))
                arg = x - 1.5 * sq2 * alogsig[m0]
                fm[m0] = 0.5 * (1.0 - erf(arg))
                if wbar > 0.0:
                    fluxn[m0] = fn[m0] * w
                    fluxm[m0] = fm[m0] * w
            flux_fullact[0] = w

    return 0


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


@inline
def _dropmixnuc_explmix_ptr(
    pver: int,
    top_lev: int,
    surfrate: float,
    flxconv: float,
    dt: float,
    is_unact: int,
    q: Ptr[float],
    q_base: int,
    src: Ptr[float],
    ekkp: Ptr[float],
    ekkm: Ptr[float],
    overlapp: Ptr[float],
    overlapm: Ptr[float],
    qold: Ptr[float],
    qold_base: int,
    qactold: Ptr[float],
    qactold_base: int,
):
    if is_unact != 0:
        for k in range(top_lev, pver + 1):
            kp1 = min(k + 1, pver)
            km1 = max(k - 1, top_lev)
            k0 = k - 1
            kp10 = kp1 - 1
            km10 = km1 - 1
            q[q_base + k0] = qold[qold_base + k0] + dt * (
                -src[k0]
                + ekkp[k0]
                * (
                    qold[qold_base + kp10]
                    - qold[qold_base + k0]
                    + qactold[qactold_base + kp10] * (1.0 - overlapp[k0])
                )
                + ekkm[k0]
                * (
                    qold[qold_base + km10]
                    - qold[qold_base + k0]
                    + qactold[qactold_base + km10] * (1.0 - overlapm[k0])
                )
            )
            q[q_base + k0] = max(q[q_base + k0], 0.0)

        q[q_base + pver - 1] = (
            q[q_base + pver - 1] - surfrate * qold[qold_base + pver - 1] * dt + flxconv * dt
        )
        q[q_base + pver - 1] = max(q[q_base + pver - 1], 0.0)
    else:
        for k in range(top_lev, pver + 1):
            kp1 = min(k + 1, pver)
            km1 = max(k - 1, top_lev)
            k0 = k - 1
            kp10 = kp1 - 1
            km10 = km1 - 1
            q[q_base + k0] = qold[qold_base + k0] + dt * (
                src[k0]
                + ekkp[k0] * (overlapp[k0] * qold[qold_base + kp10] - qold[qold_base + k0])
                + ekkm[k0] * (overlapm[k0] * qold[qold_base + km10] - qold[qold_base + k0])
            )
            q[q_base + k0] = max(q[q_base + k0], 0.0)

        q[q_base + pver - 1] = (
            q[q_base + pver - 1] - surfrate * qold[qold_base + pver - 1] * dt + flxconv * dt
        )
        q[q_base + pver - 1] = max(q[q_base + pver - 1], 0.0)


@inline
def _dropmixnuc_source_from_act_ptr(
    pver: int,
    top_lev: int,
    ncnst_tot: int,
    m: int,
    mm: int,
    nsav: int,
    taumix_internal_pver_inv: float,
    act: Ptr[float],
    raercol: Ptr[float],
    raercol_cw: Ptr[float],
    source: Ptr[float],
):
    for k in range(top_lev, pver):
        source[k - 1] = act[_mode_idx(k, m, pver)] * (
            raercol[_aero_col_idx(k + 1, mm, nsav, pver, ncnst_tot)]
        )

    tmpa = (
        raercol[_aero_col_idx(pver, mm, nsav, pver, ncnst_tot)] * act[_mode_idx(pver, m, pver)]
        + raercol_cw[_aero_col_idx(pver, mm, nsav, pver, ncnst_tot)]
        * (act[_mode_idx(pver, m, pver)] - taumix_internal_pver_inv)
    )
    source[pver - 1] = max(0.0, tmpa)


def ndrop_dropmixnuc_submix_all_codon(
    pver: int,
    top_lev: int,
    ntot_amode: int,
    ncnst_tot: int,
    dtmix: float,
    taumix_internal_pver_inv: float,
    nact_p: cobj,
    mact_p: cobj,
    mam_idx_p: cobj,
    nspec_amode_p: cobj,
    ekkp_p: cobj,
    ekkm_p: cobj,
    overlapp_p: cobj,
    overlapm_p: cobj,
    qcld_p: cobj,
    qncld_p: cobj,
    srcn_p: cobj,
    source_p: cobj,
    raercol_p: cobj,
    raercol_cw_p: cobj,
    nsav_p: cobj,
    nnew_p: cobj,
):
    nact = Ptr[float](nact_p)
    mact = Ptr[float](mact_p)
    mam_idx = Ptr[i32](mam_idx_p)
    nspec_amode = Ptr[i32](nspec_amode_p)
    ekkp = Ptr[float](ekkp_p)
    ekkm = Ptr[float](ekkm_p)
    overlapp = Ptr[float](overlapp_p)
    overlapm = Ptr[float](overlapm_p)
    qcld = Ptr[float](qcld_p)
    qncld = Ptr[float](qncld_p)
    srcn = Ptr[float](srcn_p)
    source = Ptr[float](source_p)
    raercol = Ptr[float](raercol_p)
    raercol_cw = Ptr[float](raercol_cw_p)
    nsav_ref = Ptr[i32](nsav_p)
    nnew_ref = Ptr[i32](nnew_p)

    for k in range(1, pver + 1):
        qncld[k - 1] = qcld[k - 1]
        srcn[k - 1] = 0.0

    ntemp = nsav_ref[0]
    nsav_ref[0] = nnew_ref[0]
    nnew_ref[0] = ntemp
    nsav = int(nsav_ref[0])
    nnew = int(nnew_ref[0])

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

    _dropmixnuc_explmix_ptr(
        pver,
        top_lev,
        0.0,
        0.0,
        dtmix,
        0,
        qcld,
        0,
        srcn,
        ekkp,
        ekkm,
        overlapp,
        overlapm,
        qncld,
        0,
        qncld,
        0,
    )

    for m in range(1, ntot_amode + 1):
        mm = int(mam_idx[_mam_idx(m, 0, ntot_amode)])
        _dropmixnuc_source_from_act_ptr(
            pver,
            top_lev,
            ncnst_tot,
            m,
            mm,
            nsav,
            taumix_internal_pver_inv,
            nact,
            raercol,
            raercol_cw,
            source,
        )

        cw_new = _aero_col_idx(1, mm, nnew, pver, ncnst_tot)
        cw_sav = _aero_col_idx(1, mm, nsav, pver, ncnst_tot)
        raer_new = _aero_col_idx(1, mm, nnew, pver, ncnst_tot)
        raer_sav = _aero_col_idx(1, mm, nsav, pver, ncnst_tot)
        _dropmixnuc_explmix_ptr(
            pver,
            top_lev,
            0.0,
            0.0,
            dtmix,
            0,
            raercol_cw,
            cw_new,
            source,
            ekkp,
            ekkm,
            overlapp,
            overlapm,
            raercol_cw,
            cw_sav,
            raercol_cw,
            cw_sav,
        )
        _dropmixnuc_explmix_ptr(
            pver,
            top_lev,
            0.0,
            0.0,
            dtmix,
            1,
            raercol,
            raer_new,
            source,
            ekkp,
            ekkm,
            overlapp,
            overlapm,
            raercol,
            raer_sav,
            raercol_cw,
            cw_sav,
        )

        for l in range(1, int(nspec_amode[m - 1]) + 1):
            mm = int(mam_idx[_mam_idx(m, l, ntot_amode)])
            _dropmixnuc_source_from_act_ptr(
                pver,
                top_lev,
                ncnst_tot,
                m,
                mm,
                nsav,
                taumix_internal_pver_inv,
                mact,
                raercol,
                raercol_cw,
                source,
            )

            cw_new = _aero_col_idx(1, mm, nnew, pver, ncnst_tot)
            cw_sav = _aero_col_idx(1, mm, nsav, pver, ncnst_tot)
            raer_new = _aero_col_idx(1, mm, nnew, pver, ncnst_tot)
            raer_sav = _aero_col_idx(1, mm, nsav, pver, ncnst_tot)
            _dropmixnuc_explmix_ptr(
                pver,
                top_lev,
                0.0,
                0.0,
                dtmix,
                0,
                raercol_cw,
                cw_new,
                source,
                ekkp,
                ekkm,
                overlapp,
                overlapm,
                raercol_cw,
                cw_sav,
                raercol_cw,
                cw_sav,
            )
            _dropmixnuc_explmix_ptr(
                pver,
                top_lev,
                0.0,
                0.0,
                dtmix,
                1,
                raercol,
                raer_new,
                source,
                ekkp,
                ekkm,
                overlapp,
                overlapm,
                raercol,
                raer_sav,
                raercol_cw,
                cw_sav,
            )
