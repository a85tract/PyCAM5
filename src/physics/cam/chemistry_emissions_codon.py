from math import log, sqrt

from chemistry_common_codon import _idx2, _flux_idx


def _chem_emissions_zero_cflx(
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


def _chem_emissions_megan_flux_single(
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


def _chem_emissions_megan_flux_shell(
    ncol: int,
    pcols: int,
    nmegan: int,
    megan_indices_map_p: cobj,
    megan_wght_factors_p: cobj,
    meganflx_p: cobj,
    cflx_p: cobj,
    megflx_p: cobj,
):
    megan_indices_map = Ptr[int](megan_indices_map_p)
    megan_wght_factors = Ptr[float](megan_wght_factors_p)
    meganflx = Ptr[float](meganflx_p)
    cflx = Ptr[float](cflx_p)
    megflx = Ptr[float](megflx_p)

    for n in range(1, nmegan + 1):
        megan_index = megan_indices_map[n - 1]
        megan_weight = megan_wght_factors[n - 1]
        for i in range(1, ncol + 1):
            flux = -meganflx[_idx2(i, n, pcols)] * megan_weight
            megflx[_idx2(i, n, pcols)] = flux
            cflx[_flux_idx(i, megan_index, pcols)] += flux


def _chem_emissions_add_sflx(
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


def chem_emissions_zero_cflx_codon(
    pcols: int,
    pcnst: int,
    map2chm_p: cobj,
    cflx_p: cobj,
):
    _chem_emissions_zero_cflx(
        pcols,
        pcnst,
        map2chm_p,
        cflx_p,
    )


def chem_emissions_megan_flux_codon(
    ncol: int,
    pcols: int,
    megan_index: int,
    megan_weight: float,
    meganflx_p: cobj,
    cflx_p: cobj,
    megflx_p: cobj,
):
    _chem_emissions_megan_flux_single(
        ncol,
        pcols,
        megan_index,
        megan_weight,
        meganflx_p,
        cflx_p,
        megflx_p,
    )


def chem_emissions_add_sflx_codon(
    ncol: int,
    pcols: int,
    pcnst: int,
    h2o_ndx: int,
    map2chm_p: cobj,
    cflx_p: cobj,
    sflx_p: cobj,
):
    _chem_emissions_add_sflx(
        ncol,
        pcols,
        pcnst,
        h2o_ndx,
        map2chm_p,
        cflx_p,
        sflx_p,
    )


def chem_emissions_shell_codon(
    stage: int,
    ncol: int,
    pcols: int,
    pcnst: int,
    nmegan: int,
    h2o_ndx: int,
    map2chm_p: cobj,
    megan_indices_map_p: cobj,
    megan_wght_factors_p: cobj,
    meganflx_p: cobj,
    cflx_p: cobj,
    megflx_p: cobj,
    sflx_p: cobj,
):
    if stage == 1:
        _chem_emissions_zero_cflx(
            pcols,
            pcnst,
            map2chm_p,
            cflx_p,
        )
    elif stage == 2:
        _chem_emissions_megan_flux_shell(
            ncol,
            pcols,
            nmegan,
            megan_indices_map_p,
            megan_wght_factors_p,
            meganflx_p,
            cflx_p,
            megflx_p,
        )
    elif stage == 3:
        _chem_emissions_add_sflx(
            ncol,
            pcols,
            pcnst,
            h2o_ndx,
            map2chm_p,
            cflx_p,
            sflx_p,
        )


def _aero_model_emissions_accumulate_sflx(
    ncol: int,
    pcols: int,
    nindices: int,
    indices_p: cobj,
    cflx_p: cobj,
    sflx_p: cobj,
):
    indices = Ptr[int](indices_p)
    cflx = Ptr[float](cflx_p)
    sflx = Ptr[float](sflx_p)

    for i in range(1, pcols + 1):
        sflx[i - 1] = 0.0

    for m in range(1, nindices + 1):
        idx = indices[m - 1]
        for i in range(1, ncol + 1):
            sflx[i - 1] += cflx[_flux_idx(i, idx, pcols)]


def aero_model_emissions_accumulate_sflx_codon(
    ncol: int,
    pcols: int,
    nindices: int,
    indices_p: cobj,
    cflx_p: cobj,
    sflx_p: cobj,
):
    _aero_model_emissions_accumulate_sflx(
        ncol,
        pcols,
        nindices,
        indices_p,
        cflx_p,
        sflx_p,
    )


def aero_model_emissions_seasalt_wind_codon(
    ncol: int,
    pcols: int,
    pver: int,
    z0: float,
    state_u_p: cobj,
    state_v_p: cobj,
    state_zm_p: cobj,
    u10cubed_p: cobj,
):
    state_u = Ptr[float](state_u_p)
    state_v = Ptr[float](state_v_p)
    state_zm = Ptr[float](state_zm_p)
    u10cubed = Ptr[float](u10cubed_p)

    for i in range(1, ncol + 1):
        wind = sqrt(state_u[_idx2(i, pver, pcols)] ** 2 + state_v[_idx2(i, pver, pcols)] ** 2)
        wind = wind * log(10.0 / z0) / log(state_zm[_idx2(i, pver, pcols)] / z0)
        u10cubed[i - 1] = wind ** 3.41


def _aero_model_emissions_dust_shell(
    ncol: int,
    pcols: int,
    dust_nbin: int,
    soil_erod_fact: float,
    soil_erod_p: cobj,
    dust_flux_sum_p: cobj,
    cflx_p: cobj,
    sflx_p: cobj,
    dust_indices_p: cobj,
    dust_emis_sclfctr_p: cobj,
    dust_x_mton_p: cobj,
):
    soil_erod = Ptr[float](soil_erod_p)
    dust_flux_sum = Ptr[float](dust_flux_sum_p)
    cflx = Ptr[float](cflx_p)
    sflx = Ptr[float](sflx_p)
    dust_indices = Ptr[int](dust_indices_p)
    dust_emis_sclfctr = Ptr[float](dust_emis_sclfctr_p)
    dust_x_mton = Ptr[float](dust_x_mton_p)

    for i in range(1, ncol + 1):
        for m in range(1, dust_nbin + 1):
            idst = dust_indices[m - 1]
            value = dust_flux_sum[i - 1] * dust_emis_sclfctr[m - 1]
            value = value * soil_erod[i - 1] / soil_erod_fact
            value = value * 1.15
            cflx[_flux_idx(i, idst, pcols)] = value

            inum = dust_indices[dust_nbin + m - 1]
            cflx[_flux_idx(i, inum, pcols)] = value * dust_x_mton[m - 1]

    _aero_model_emissions_accumulate_sflx(
        ncol,
        pcols,
        dust_nbin,
        dust_indices_p,
        cflx_p,
        sflx_p,
    )


def _aero_model_emissions_seasalt_shell(
    ncol: int,
    pcols: int,
    seasalt_nbin: int,
    nsections: int,
    seasalt_emis_scale: float,
    pi_val: float,
    seasalt_density: float,
    fi_p: cobj,
    ocnfrac_p: cobj,
    cflx_p: cobj,
    sflx_p: cobj,
    seasalt_indices_p: cobj,
    seasalt_sz_range_lo_p: cobj,
    seasalt_sz_range_hi_p: cobj,
    dg_p: cobj,
    rdry_p: cobj,
):
    fi = Ptr[float](fi_p)
    ocnfrac = Ptr[float](ocnfrac_p)
    cflx = Ptr[float](cflx_p)
    seasalt_indices = Ptr[int](seasalt_indices_p)
    seasalt_sz_range_lo = Ptr[float](seasalt_sz_range_lo_p)
    seasalt_sz_range_hi = Ptr[float](seasalt_sz_range_hi_p)
    dg = Ptr[float](dg_p)
    rdry = Ptr[float](rdry_p)

    for ibin in range(1, seasalt_nbin + 1):
        mm = seasalt_indices[ibin - 1]
        mn = seasalt_indices[seasalt_nbin + ibin - 1]

        if mn > 0:
            for isec in range(1, nsections + 1):
                if dg[isec - 1] >= seasalt_sz_range_lo[ibin - 1] and dg[isec - 1] < seasalt_sz_range_hi[ibin - 1]:
                    for i in range(1, ncol + 1):
                        term = fi[_idx2(i, isec, pcols)] * ocnfrac[i - 1]
                        term = term * seasalt_emis_scale
                        cflx[_flux_idx(i, mn, pcols)] += term

        for i in range(1, ncol + 1):
            cflx[_flux_idx(i, mm, pcols)] = 0.0

        for isec in range(1, nsections + 1):
            if dg[isec - 1] >= seasalt_sz_range_lo[ibin - 1] and dg[isec - 1] < seasalt_sz_range_hi[ibin - 1]:
                for i in range(1, ncol + 1):
                    term = fi[_idx2(i, isec, pcols)] * ocnfrac[i - 1]
                    term = term * seasalt_emis_scale
                    term = term * 4.0
                    term = term / 3.0
                    term = term * pi_val
                    term = term * rdry[isec - 1] ** 3
                    term = term * seasalt_density
                    cflx[_flux_idx(i, mm, pcols)] += term

    _aero_model_emissions_accumulate_sflx(
        ncol,
        pcols,
        seasalt_nbin,
        seasalt_indices_p,
        cflx_p,
        sflx_p,
    )


def aero_model_emissions_shell_codon(
    stage: int,
    ncol: int,
    pcols: int,
    ndstflx: int,
    dust_nbin: int,
    seasalt_nbin: int,
    nsections: int,
    soil_erod_fact: float,
    seasalt_emis_scale: float,
    pi_val: float,
    seasalt_density: float,
    soil_erod_p: cobj,
    dust_flux_sum_p: cobj,
    fi_p: cobj,
    ocnfrac_p: cobj,
    cflx_p: cobj,
    sflx_p: cobj,
    dust_indices_p: cobj,
    dust_emis_sclfctr_p: cobj,
    dust_x_mton_p: cobj,
    seasalt_indices_p: cobj,
    seasalt_sz_range_lo_p: cobj,
    seasalt_sz_range_hi_p: cobj,
    dg_p: cobj,
    rdry_p: cobj,
):
    if stage == 1:
        _aero_model_emissions_dust_shell(
            ncol,
            pcols,
            dust_nbin,
            soil_erod_fact,
            soil_erod_p,
            dust_flux_sum_p,
            cflx_p,
            sflx_p,
            dust_indices_p,
            dust_emis_sclfctr_p,
            dust_x_mton_p,
        )
    elif stage == 2:
        _aero_model_emissions_seasalt_shell(
            ncol,
            pcols,
            seasalt_nbin,
            nsections,
            seasalt_emis_scale,
            pi_val,
            seasalt_density,
            fi_p,
            ocnfrac_p,
            cflx_p,
            sflx_p,
            seasalt_indices_p,
            seasalt_sz_range_lo_p,
            seasalt_sz_range_hi_p,
            dg_p,
            rdry_p,
        )
