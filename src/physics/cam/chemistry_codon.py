@inline
def _idx2(i: int, k: int, ld1: int) -> int:
    return (i - 1) + (k - 1) * ld1


@inline
def _idx3(i: int, k: int, m: int, ld1: int, ld2: int) -> int:
    return (i - 1) + (k - 1) * ld1 + (m - 1) * ld1 * ld2


@inline
def _flux_idx(i: int, m: int, pcols: int) -> int:
    return (i - 1) + (m - 1) * pcols


@export
def chem_emissions_zero_cflx_codon(
    pcols: int,
    pcnst: int,
    map2chm_p: cobj,
    cflx_p: cobj,
):
    map2chm = Ptr[int](map2chm_p)
    cflx = Ptr[float](cflx_p)

    for m in range(2, pcnst + 1):
        if map2chm[m - 1] > 0:
            for i in range(1, pcols + 1):
                cflx[_flux_idx(i, m, pcols)] = 0.0


@export
def chem_emissions_megan_flux_codon(
    ncol: int,
    pcols: int,
    megan_index: int,
    megan_weight: float,
    meganflx_p: cobj,
    cflx_p: cobj,
    megflx_p: cobj,
):
    meganflx = Ptr[float](meganflx_p)
    cflx = Ptr[float](cflx_p)
    megflx = Ptr[float](megflx_p)

    for i in range(1, ncol + 1):
        flux = -meganflx[i - 1] * megan_weight
        megflx[i - 1] = flux
        cflx[_flux_idx(i, megan_index, pcols)] += flux


@export
def chem_emissions_add_sflx_codon(
    ncol: int,
    pcols: int,
    pcnst: int,
    h2o_ndx: int,
    map2chm_p: cobj,
    cflx_p: cobj,
    sflx_p: cobj,
):
    map2chm = Ptr[int](map2chm_p)
    cflx = Ptr[float](cflx_p)
    sflx = Ptr[float](sflx_p)

    for m in range(1, pcnst + 1):
        n = map2chm[m - 1]
        if n > 0 and n != h2o_ndx:
            for i in range(1, ncol + 1):
                cflx[_flux_idx(i, m, pcols)] += sflx[_flux_idx(i, n, pcols)]


@export
def chem_timestep_init_should_run_codon(
    nstep: int,
    chem_freq: int,
    chem_step_flag_p: cobj,
):
    chem_step_flag = Ptr[int](chem_step_flag_p)
    chem_step_flag[0] = 1 if nstep % chem_freq == 0 else 0


@export
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


@export
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


@export
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


@export
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


@export
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


@export
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


@export
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


@export
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


@export
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


@export
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


@export
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


@export
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


@export
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


@export
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


@export
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


@export
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


@export
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


@export
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


@export
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


@export
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


@export
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


@export
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


@export
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


@export
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


@export
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


@export
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


@export
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


@export
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


@export
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


@export
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


@export
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
