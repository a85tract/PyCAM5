@export
def phys_timestep_init_select_branches_codon(
    cam3_aero_on: int,
    cam3_ozone_on: int,
    do_waccm_ions: int,
    branch_mask_p: cobj,
):
    branch_mask = Ptr[int](branch_mask_p)

    mask = 0
    if cam3_aero_on != 0:
        mask |= 1
    if cam3_ozone_on != 0:
        mask |= 2
    if do_waccm_ions != 0:
        mask |= 4

    branch_mask[0] = mask


@inline
def _idx(i: int) -> int:
    return i - 1


@inline
def _field2_idx(i: int, k: int, ld1: int) -> int:
    """state_t/tini/tend_dtdt/dtcore/tmp_t declared as (ld1, pver)"""
    return (i - 1) + (k - 1) * ld1


@inline
def _field3_idx(i: int, k: int, m: int, ld1: int, ld2: int) -> int:
    """state_q declared as (ld1, ld2, pcnst)"""
    return (i - 1) + (k - 1) * ld1 + (m - 1) * ld1 * ld2


@export
def tphysbc_precip_ops_codon(
    mode: int,
    ncol: int,
    pcols: int,
    cld_macmic_num_steps: int,
    prec_sed_macmic_p: cobj,
    snow_sed_macmic_p: cobj,
    prec_pcw_macmic_p: cobj,
    snow_pcw_macmic_p: cobj,
    prec_sed_p: cobj,
    snow_sed_p: cobj,
    prec_pcw_p: cobj,
    snow_pcw_p: cobj,
    prec_str_p: cobj,
    snow_str_p: cobj,
    prec_sed_carma_p: cobj,
    snow_sed_carma_p: cobj,
):
    prec_sed_macmic = Ptr[float](prec_sed_macmic_p)
    snow_sed_macmic = Ptr[float](snow_sed_macmic_p)
    prec_pcw_macmic = Ptr[float](prec_pcw_macmic_p)
    snow_pcw_macmic = Ptr[float](snow_pcw_macmic_p)
    prec_sed = Ptr[float](prec_sed_p)
    snow_sed = Ptr[float](snow_sed_p)
    prec_pcw = Ptr[float](prec_pcw_p)
    snow_pcw = Ptr[float](snow_pcw_p)
    prec_str = Ptr[float](prec_str_p)
    snow_str = Ptr[float](snow_str_p)
    prec_sed_carma = Ptr[float](prec_sed_carma_p)
    snow_sed_carma = Ptr[float](snow_sed_carma_p)

    if mode == 0:
        for i in range(1, pcols + 1):
            prec_sed_macmic[_idx(i)] = 0.0
            snow_sed_macmic[_idx(i)] = 0.0
            prec_pcw_macmic[_idx(i)] = 0.0
            snow_pcw_macmic[_idx(i)] = 0.0
    elif mode == 1:
        for i in range(1, ncol + 1):
            prec_sed_macmic[_idx(i)] = prec_sed_macmic[_idx(i)] + prec_sed[_idx(i)]
            snow_sed_macmic[_idx(i)] = snow_sed_macmic[_idx(i)] + snow_sed[_idx(i)]
            prec_pcw_macmic[_idx(i)] = prec_pcw_macmic[_idx(i)] + prec_pcw[_idx(i)]
            snow_pcw_macmic[_idx(i)] = snow_pcw_macmic[_idx(i)] + snow_pcw[_idx(i)]
    elif mode == 2:
        for i in range(1, ncol + 1):
            prec_sed[_idx(i)] = prec_sed_macmic[_idx(i)] / cld_macmic_num_steps
            snow_sed[_idx(i)] = snow_sed_macmic[_idx(i)] / cld_macmic_num_steps
            prec_pcw[_idx(i)] = prec_pcw_macmic[_idx(i)] / cld_macmic_num_steps
            snow_pcw[_idx(i)] = snow_pcw_macmic[_idx(i)] / cld_macmic_num_steps
            prec_str[_idx(i)] = prec_pcw[_idx(i)] + prec_sed[_idx(i)]
            snow_str[_idx(i)] = snow_pcw[_idx(i)] + snow_sed[_idx(i)]
    elif mode == 3:
        for i in range(1, ncol + 1):
            prec_sed[_idx(i)] = prec_sed[_idx(i)] + prec_sed_carma[_idx(i)]
            snow_sed[_idx(i)] = snow_sed[_idx(i)] + snow_sed_carma[_idx(i)]


@export
def tphysac_t_update_codon(
    ncol: int,
    pcols: int,
    pver: int,
    ztodt: float,
    state_t_p: cobj,
    tini_p: cobj,
    tend_dtdt_p: cobj,
    dtcore_p: cobj,
    tmp_t_p: cobj,
):
    state_t = Ptr[float](state_t_p)
    tini = Ptr[float](tini_p)
    tend_dtdt = Ptr[float](tend_dtdt_p)
    dtcore = Ptr[float](dtcore_p)
    tmp_t = Ptr[float](tmp_t_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            tmp_t[_field2_idx(i, k, pcols)] = state_t[_field2_idx(i, k, pcols)]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            state_t[_field2_idx(i, k, pcols)] = tini[_field2_idx(i, k, pcols)] + ztodt * tend_dtdt[
                _field2_idx(i, k, pcols)
            ]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            dtcore[_field2_idx(i, k, pcols)] = state_t[_field2_idx(i, k, pcols)]


@export
def tphysac_q_snapshot_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    ixcldliq: int,
    ixcldice: int,
    state_q_p: cobj,
    tmp_q_p: cobj,
    tmp_cldliq_p: cobj,
    tmp_cldice_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    tmp_q = Ptr[float](tmp_q_p)
    tmp_cldliq = Ptr[float](tmp_cldliq_p)
    tmp_cldice = Ptr[float](tmp_cldice_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            tmp_q[_field2_idx(i, k, pcols)] = state_q[_field3_idx(i, k, 1, pcols, pver)]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            tmp_cldliq[_field2_idx(i, k, pcols)] = state_q[_field3_idx(i, k, ixcldliq, pcols, pver)]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            tmp_cldice[_field2_idx(i, k, pcols)] = state_q[_field3_idx(i, k, ixcldice, pcols, pver)]


@export
def tphysbc_qini_snapshot_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    ixcldliq: int,
    ixcldice: int,
    state_q_p: cobj,
    qini_p: cobj,
    cldliqini_p: cobj,
    cldiceini_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    qini = Ptr[float](qini_p)
    cldliqini = Ptr[float](cldliqini_p)
    cldiceini = Ptr[float](cldiceini_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            qini[_field2_idx(i, k, pcols)] = state_q[_field3_idx(i, k, 1, pcols, pver)]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            cldliqini[_field2_idx(i, k, pcols)] = state_q[_field3_idx(i, k, ixcldliq, pcols, pver)]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            cldiceini[_field2_idx(i, k, pcols)] = state_q[_field3_idx(i, k, ixcldice, pcols, pver)]


@export
def tphysbc_dadadj_input_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    state_t_p: cobj,
    state_q_p: cobj,
    ptend_s_p: cobj,
    ptend_q_p: cobj,
):
    state_t = Ptr[float](state_t_p)
    state_q = Ptr[float](state_q_p)
    ptend_s = Ptr[float](ptend_s_p)
    ptend_q = Ptr[float](ptend_q_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            ptend_s[_field2_idx(i, k, pcols)] = state_t[_field2_idx(i, k, pcols)]

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            ptend_q[_field3_idx(i, k, 1, pcols, pver)] = state_q[_field3_idx(i, k, 1, pcols, pver)]


@export
def tphysbc_dadadj_output_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    ztodt: float,
    cpair_local: float,
    state_t_p: cobj,
    state_q_p: cobj,
    ptend_s_p: cobj,
    ptend_q_p: cobj,
):
    state_t = Ptr[float](state_t_p)
    state_q = Ptr[float](state_q_p)
    ptend_s = Ptr[float](ptend_s_p)
    ptend_q = Ptr[float](ptend_q_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            ptend_s[_field2_idx(i, k, pcols)] = (
                (ptend_s[_field2_idx(i, k, pcols)] - state_t[_field2_idx(i, k, pcols)]) / ztodt
            ) * cpair_local

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            ptend_q[_field3_idx(i, k, 1, pcols, pver)] = (
                ptend_q[_field3_idx(i, k, 1, pcols, pver)] - state_q[_field3_idx(i, k, 1, pcols, pver)]
            ) / ztodt


@export
def tphysbc_dtcore_update_codon(
    ncol: int,
    pcols: int,
    pver: int,
    ztodt: float,
    tini_p: cobj,
    dtcore_p: cobj,
    tend_dtdt_p: cobj,
):
    tini = Ptr[float](tini_p)
    dtcore = Ptr[float](dtcore_p)
    tend_dtdt = Ptr[float](tend_dtdt_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            dtcore[_field2_idx(i, k, pcols)] = (
                (tini[_field2_idx(i, k, pcols)] - dtcore[_field2_idx(i, k, pcols)]) / ztodt
            ) + tend_dtdt[_field2_idx(i, k, pcols)]


@export
def tphysbc_init_fields_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    fracis_p: cobj,
    tend_dtdt_p: cobj,
    tend_dudt_p: cobj,
    tend_dvdt_p: cobj,
):
    fracis = Ptr[float](fracis_p)
    tend_dtdt = Ptr[float](tend_dtdt_p)
    tend_dudt = Ptr[float](tend_dudt_p)
    tend_dvdt = Ptr[float](tend_dvdt_p)

    for m in range(1, pcnst + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                fracis[_field3_idx(i, k, m, pcols, pver)] = 1.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            tend_dtdt[_field2_idx(i, k, pcols)] = 0.0
            tend_dudt[_field2_idx(i, k, pcols)] = 0.0
            tend_dvdt[_field2_idx(i, k, pcols)] = 0.0


@export
def tphysbc_tini_copy_codon(
    ncol: int,
    pcols: int,
    pver: int,
    state_t_p: cobj,
    tini_p: cobj,
):
    state_t = Ptr[float](state_t_p)
    tini = Ptr[float](tini_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            tini[_field2_idx(i, k, pcols)] = state_t[_field2_idx(i, k, pcols)]


@export
def tphysbc_flx_cnd_sum_codon(
    ncol: int,
    pcols: int,
    a_p: cobj,
    b_p: cobj,
    out_p: cobj,
):
    a = Ptr[float](a_p)
    b = Ptr[float](b_p)
    out = Ptr[float](out_p)

    for i in range(1, ncol + 1):
        out[_idx(i)] = a[_idx(i)] + b[_idx(i)]


@export
def tphysbc_macrop_fluxes_codon(
    mode: int,
    ncol: int,
    pcols: int,
    rliq_p: cobj,
    det_s_p: cobj,
    flx_cnd_p: cobj,
    flx_heat_p: cobj,
    shf_p: cobj,
):
    rliq = Ptr[float](rliq_p)
    det_s = Ptr[float](det_s_p)
    flx_cnd = Ptr[float](flx_cnd_p)
    flx_heat = Ptr[float](flx_heat_p)

    for i in range(1, ncol + 1):
        flx_cnd[_idx(i)] = -1.0 * rliq[_idx(i)]

    if mode == 1:
        for i in range(1, ncol + 1):
            flx_heat[_idx(i)] = det_s[_idx(i)]
    else:
        shf = Ptr[float](shf_p)
        for i in range(1, ncol + 1):
            flx_heat[_idx(i)] = shf[_idx(i)] + det_s[_idx(i)]
