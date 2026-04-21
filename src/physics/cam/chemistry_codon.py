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
