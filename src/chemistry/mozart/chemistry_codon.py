import chemistry_common_codon as _common
import chemistry_aero_bridge_codon as _aero_bridge
import chemistry_diags_codon as _diags
import chemistry_emissions_codon as _emissions
import chemistry_gas_phase_codon as _gas_phase
import chemistry_photolysis_codon as _photolysis
import chemistry_wetchem_codon as _wetchem

@export
def chemistry_misc_touch_codon(tag: int) -> int:
    return _common.chemistry_misc_touch_codon(tag)

@export
def chem_lookup_name_codon(name_len: int, name_ascii_p: cobj, list_len: int, list_ascii_p: cobj, list_count: int) -> int:
    return _common.chem_lookup_name_codon(name_len, name_ascii_p, list_len, list_ascii_p, list_count)

@export
def has_drydep_codon(name_len: int, name_ascii_p: cobj, list_len: int, list_ascii_p: cobj, list_count: int) -> int:
    return _common.has_drydep_codon(name_len, name_ascii_p, list_len, list_ascii_p, list_count)

@export
def chem_lookup_mapped_name_codon(
    name_len: int,
    name_ascii_p: cobj,
    list_len: int,
    list_ascii_p: cobj,
    map_p: cobj,
    list_count: int,
) -> int:
    return _common.chem_lookup_mapped_name_codon(name_len, name_ascii_p, list_len, list_ascii_p, map_p, list_count)

@export
def get_spc_ndx_codon(name_len: int, name_ascii_p: cobj, list_len: int, list_ascii_p: cobj, list_count: int) -> int:
    return _common.chem_lookup_name_codon(name_len, name_ascii_p, list_len, list_ascii_p, list_count)

@export
def get_inv_ndx_codon(name_len: int, name_ascii_p: cobj, list_len: int, list_ascii_p: cobj, list_count: int) -> int:
    return _common.chem_lookup_name_codon(name_len, name_ascii_p, list_len, list_ascii_p, list_count)

@export
def get_extfrc_ndx_codon(name_len: int, name_ascii_p: cobj, list_len: int, list_ascii_p: cobj, list_count: int) -> int:
    return _common.chem_lookup_name_codon(name_len, name_ascii_p, list_len, list_ascii_p, list_count)

@export
def get_rxt_ndx_codon(
    name_len: int,
    name_ascii_p: cobj,
    list_len: int,
    list_ascii_p: cobj,
    map_p: cobj,
    list_count: int,
) -> int:
    return _common.chem_lookup_mapped_name_codon(name_len, name_ascii_p, list_len, list_ascii_p, map_p, list_count)

@inline
def _chem_name_is_h2o(name_len: int, name_ascii_p: cobj) -> bool:
    name_ascii = Ptr[int](name_ascii_p)
    name_trim = name_len
    while name_trim > 0 and name_ascii[name_trim - 1] == 32:
        name_trim -= 1
    return (
        name_trim == 3
        and name_ascii[0] == 72
        and name_ascii[1] == 50
        and name_ascii[2] == 79
    )

@export
def chem_is_active_codon(active: int) -> int:
    return active

@export
def chem_is_codon(active: int) -> int:
    return chem_is_active_codon(active)

@export
def chem_implements_cnst_codon(
    name_len: int,
    name_ascii_p: cobj,
    solsym_len: int,
    solsym_ascii_p: cobj,
    gas_pcnst: int,
    inv_len: int,
    inv_ascii_p: cobj,
    nfs: int,
) -> int:
    if _chem_name_is_h2o(name_len, name_ascii_p):
        return 0
    if _common.chem_lookup_name_codon(name_len, name_ascii_p, solsym_len, solsym_ascii_p, gas_pcnst) > 0:
        return 1
    if _common.chem_lookup_name_codon(name_len, name_ascii_p, inv_len, inv_ascii_p, nfs) > 0:
        return 1
    return 0

@export
def slvd_index_codon(name_len: int, name_ascii_p: cobj, list_len: int, list_ascii_p: cobj, list_count: int) -> int:
    return _common.chem_lookup_name_codon(name_len, name_ascii_p, list_len, list_ascii_p, list_count)

@export
def init_cfc11star_codon(active: int) -> int:
    return active

@export
def tracer_cnst_adv_codon(active: int) -> int:
    return active

@export
def get_cnst_data_codon(active: int) -> int:
    return active

@export
def tracer_srcs_adv_codon(active: int) -> int:
    return active

@export
def linoz_data_adv_codon(active: int) -> int:
    return active

@export
def noy_ubc_init_codon(active: int) -> int:
    return active

@export
def noy_ubc_set_codon(active: int) -> int:
    return active

@export
def register_short_lived_species_codon(active: int) -> int:
    return active

@export
def get_short_lived_species_codon(active: int) -> int:
    return active

@export
def init_cph_codon(active: int) -> int:
    return active

@export
def spe_prod_codon(active: int, ncol: int, pver: int, noxprod_p: cobj, hoxprod_p: cobj) -> int:
    noxprod = Ptr[float](noxprod_p)
    hoxprod = Ptr[float](hoxprod_p)
    total = ncol * pver
    for i in range(total):
        noxprod[i] = 0.0
        hoxprod[i] = 0.0
    return active

@export
def tracer_cnst_defaultopts_codon(
    file_len: int,
    file_p: cobj,
    filelist_len: int,
    filelist_p: cobj,
    datapath_len: int,
    datapath_p: cobj,
    type_len: int,
    type_p: cobj,
    specifier_len: int,
    specifier_count: int,
    specifier_p: cobj,
    rmfile: int,
    cycle_yr: int,
    fixed_ymd: int,
    fixed_tod: int,
    present_file: int,
    file_out_p: cobj,
    present_filelist: int,
    filelist_out_p: cobj,
    present_datapath: int,
    datapath_out_p: cobj,
    present_type: int,
    type_out_p: cobj,
    present_specifier: int,
    specifier_out_p: cobj,
    present_rmfile: int,
    present_cycle_yr: int,
    present_fixed_ymd: int,
    present_fixed_tod: int,
    scalar_out_p: cobj,
):
    return _common.tracer_defaultopts_codon(
        file_len,
        file_p,
        filelist_len,
        filelist_p,
        datapath_len,
        datapath_p,
        type_len,
        type_p,
        specifier_len,
        specifier_count,
        specifier_p,
        rmfile,
        cycle_yr,
        fixed_ymd,
        fixed_tod,
        present_file,
        file_out_p,
        present_filelist,
        filelist_out_p,
        present_datapath,
        datapath_out_p,
        present_type,
        type_out_p,
        present_specifier,
        specifier_out_p,
        present_rmfile,
        present_cycle_yr,
        present_fixed_ymd,
        present_fixed_tod,
        scalar_out_p,
    )

@export
def tracer_srcs_defaultopts_codon(
    file_len: int,
    file_p: cobj,
    filelist_len: int,
    filelist_p: cobj,
    datapath_len: int,
    datapath_p: cobj,
    type_len: int,
    type_p: cobj,
    specifier_len: int,
    specifier_count: int,
    specifier_p: cobj,
    rmfile: int,
    cycle_yr: int,
    fixed_ymd: int,
    fixed_tod: int,
    present_file: int,
    file_out_p: cobj,
    present_filelist: int,
    filelist_out_p: cobj,
    present_datapath: int,
    datapath_out_p: cobj,
    present_type: int,
    type_out_p: cobj,
    present_specifier: int,
    specifier_out_p: cobj,
    present_rmfile: int,
    present_cycle_yr: int,
    present_fixed_ymd: int,
    present_fixed_tod: int,
    scalar_out_p: cobj,
):
    return _common.tracer_defaultopts_codon(
        file_len,
        file_p,
        filelist_len,
        filelist_p,
        datapath_len,
        datapath_p,
        type_len,
        type_p,
        specifier_len,
        specifier_count,
        specifier_p,
        rmfile,
        cycle_yr,
        fixed_ymd,
        fixed_tod,
        present_file,
        file_out_p,
        present_filelist,
        filelist_out_p,
        present_datapath,
        datapath_out_p,
        present_type,
        type_out_p,
        present_specifier,
        specifier_out_p,
        present_rmfile,
        present_cycle_yr,
        present_fixed_ymd,
        present_fixed_tod,
        scalar_out_p,
    )

@export
def linoz_data_defaultopts_codon(
    file_len: int,
    file_p: cobj,
    filelist_len: int,
    filelist_p: cobj,
    datapath_len: int,
    datapath_p: cobj,
    type_len: int,
    type_p: cobj,
    rmfile: int,
    cycle_yr: int,
    fixed_ymd: int,
    fixed_tod: int,
    present_file: int,
    file_out_p: cobj,
    present_filelist: int,
    filelist_out_p: cobj,
    present_datapath: int,
    datapath_out_p: cobj,
    present_type: int,
    type_out_p: cobj,
    present_rmfile: int,
    present_cycle_yr: int,
    present_fixed_ymd: int,
    present_fixed_tod: int,
    scalar_out_p: cobj,
):
    return _common.linoz_defaultopts_codon(
        file_len,
        file_p,
        filelist_len,
        filelist_p,
        datapath_len,
        datapath_p,
        type_len,
        type_p,
        rmfile,
        cycle_yr,
        fixed_ymd,
        fixed_tod,
        present_file,
        file_out_p,
        present_filelist,
        filelist_out_p,
        present_datapath,
        datapath_out_p,
        present_type,
        type_out_p,
        present_rmfile,
        present_cycle_yr,
        present_fixed_ymd,
        present_fixed_tod,
        scalar_out_p,
    )

@export
def drydep_update_codon(active: int) -> int:
    return active

@export
def aircraft_emit_init_codon(active: int) -> int:
    return active

@export
def prescribed_aero_init_codon(active: int) -> int:
    return active

@export
def prescribed_ghg_init_codon(active: int) -> int:
    return active

@export
def prescribed_ozone_init_codon(active: int) -> int:
    return active

@export
def prescribed_ozone_adv_codon(active: int) -> int:
    return active

@export
def prescribed_ozone_readnl_codon(
    name_len: int,
    name_p: cobj,
    file_len: int,
    file_p: cobj,
    filelist_len: int,
    filelist_p: cobj,
    datapath_len: int,
    datapath_p: cobj,
    type_len: int,
    type_p: cobj,
    rmfile: int,
    cycle_yr: int,
    fixed_ymd: int,
    fixed_tod: int,
    name_out_p: cobj,
    file_out_p: cobj,
    filelist_out_p: cobj,
    datapath_out_p: cobj,
    type_out_p: cobj,
    scalar_out_p: cobj,
):
    _common.prescribed_field_readnl_codon(
        name_len,
        name_p,
        file_len,
        file_p,
        filelist_len,
        filelist_p,
        datapath_len,
        datapath_p,
        type_len,
        type_p,
        rmfile,
        cycle_yr,
        fixed_ymd,
        fixed_tod,
        name_out_p,
        file_out_p,
        filelist_out_p,
        datapath_out_p,
        type_out_p,
        scalar_out_p,
    )

@export
def prescribed_strataero_init_codon(active: int) -> int:
    return active

@export
def prescribed_volcaero_init_codon(active: int) -> int:
    return active

@export
def prescribed_volcaero_readnl_codon(
    name_len: int,
    name_p: cobj,
    file_len: int,
    file_p: cobj,
    filelist_len: int,
    filelist_p: cobj,
    datapath_len: int,
    datapath_p: cobj,
    type_len: int,
    type_p: cobj,
    rmfile: int,
    cycle_yr: int,
    fixed_ymd: int,
    fixed_tod: int,
    name_out_p: cobj,
    file_out_p: cobj,
    filelist_out_p: cobj,
    datapath_out_p: cobj,
    type_out_p: cobj,
    scalar_out_p: cobj,
):
    _common.prescribed_field_readnl_codon(
        name_len,
        name_p,
        file_len,
        file_p,
        filelist_len,
        filelist_p,
        datapath_len,
        datapath_p,
        type_len,
        type_p,
        rmfile,
        cycle_yr,
        fixed_ymd,
        fixed_tod,
        name_out_p,
        file_out_p,
        filelist_out_p,
        datapath_out_p,
        type_out_p,
        scalar_out_p,
    )

@export
def clybry_fam_set_codon(active: int) -> int:
    return active

@export
def clybry_fam_adj_codon(active: int) -> int:
    return active

@export
def aurora_inti_codon(active: int) -> int:
    return active

@inline
def _set_rates_idx(i: int, k: int, m: int, d1: int, d2: int) -> int:
    return i + d1 * (k + d2 * m)

@export
def set_rates_codon(
    ncol: int,
    rxt_d1: int,
    rxt_d2: int,
    sol_d1: int,
    sol_d2: int,
    rxt_rates_p: cobj,
    sol_p: cobj,
):
    rxt_rates = Ptr[float](rxt_rates_p)
    sol = Ptr[float](sol_p)
    k = 0
    while k < rxt_d2:
        i = 0
        while i < ncol:
            rxt_rates[_set_rates_idx(i, k, 0, rxt_d1, rxt_d2)] = (
                rxt_rates[_set_rates_idx(i, k, 0, rxt_d1, rxt_d2)]
                * sol[_set_rates_idx(i, k, 0, sol_d1, sol_d2)]
            )
            rxt_rates[_set_rates_idx(i, k, 2, rxt_d1, rxt_d2)] = (
                rxt_rates[_set_rates_idx(i, k, 2, rxt_d1, rxt_d2)]
                * sol[_set_rates_idx(i, k, 0, sol_d1, sol_d2)]
            )
            rxt_rates[_set_rates_idx(i, k, 3, rxt_d1, rxt_d2)] = (
                rxt_rates[_set_rates_idx(i, k, 3, rxt_d1, rxt_d2)]
                * sol[_set_rates_idx(i, k, 2, sol_d1, sol_d2)]
            )
            rxt_rates[_set_rates_idx(i, k, 4, rxt_d1, rxt_d2)] = (
                rxt_rates[_set_rates_idx(i, k, 4, rxt_d1, rxt_d2)]
                * sol[_set_rates_idx(i, k, 3, sol_d1, sol_d2)]
            )
            rxt_rates[_set_rates_idx(i, k, 5, rxt_d1, rxt_d2)] = (
                rxt_rates[_set_rates_idx(i, k, 5, rxt_d1, rxt_d2)]
                * sol[_set_rates_idx(i, k, 3, sol_d1, sol_d2)]
            )
            rxt_rates[_set_rates_idx(i, k, 6, rxt_d1, rxt_d2)] = (
                rxt_rates[_set_rates_idx(i, k, 6, rxt_d1, rxt_d2)]
                * sol[_set_rates_idx(i, k, 3, sol_d1, sol_d2)]
            )
            i += 1
        k += 1

@export
def init_mean_mass_ids_codon(lookup_ids_p: cobj, species_ids_p: cobj):
    return _common.init_mean_mass_ids_codon(lookup_ids_p, species_ids_p)

@export
def init_mean_mass_codon(lookup_ids_p: cobj, species_ids_p: cobj):
    return init_mean_mass_ids_codon(lookup_ids_p, species_ids_p)

@export
def init_hrates_ids_codon(lookup_ids_p: cobj, ptop_ref: float, psurf_ref: float, ids_p: cobj, has_hrates_p: cobj):
    return _common.init_hrates_ids_codon(lookup_ids_p, ptop_ref, psurf_ref, ids_p, has_hrates_p)

@export
def init_hrates_codon(lookup_ids_p: cobj, ptop_ref: float, psurf_ref: float, ids_p: cobj, has_hrates_p: cobj):
    return init_hrates_ids_codon(lookup_ids_p, ptop_ref, psurf_ref, ids_p, has_hrates_p)

@export
def clybry_fam_init_ids_codon(lookup_ids_p: cobj, ids_p: cobj, has_clybry_p: cobj):
    return _common.clybry_fam_init_ids_codon(lookup_ids_p, ids_p, has_clybry_p)

@export
def clybry_fam_init_codon(lookup_ids_p: cobj, ids_p: cobj, has_clybry_p: cobj):
    return clybry_fam_init_ids_codon(lookup_ids_p, ids_p, has_clybry_p)

@export
def init_strato_rates_ids_codon(lookup_ids_p: cobj, ids_p: cobj, has_strato_chem_p: cobj):
    return _common.init_strato_rates_ids_codon(lookup_ids_p, ids_p, has_strato_chem_p)

@export
def init_strato_rates_codon(lookup_ids_p: cobj, ids_p: cobj, has_strato_chem_p: cobj):
    return init_strato_rates_ids_codon(lookup_ids_p, ids_p, has_strato_chem_p)

@export
def setinv_inti_ids_codon(lookup_ids_p: cobj, ids_p: cobj, flags_p: cobj):
    return _common.setinv_inti_ids_codon(lookup_ids_p, ids_p, flags_p)

@export
def setinv_inti_codon(lookup_ids_p: cobj, ids_p: cobj, flags_p: cobj):
    return setinv_inti_ids_codon(lookup_ids_p, ids_p, flags_p)

@export
def gas_wetdep_readnl_status_codon(pcnst: int, list_p: cobj, method_p: cobj, status_p: cobj):
    return _common.gas_wetdep_readnl_status_codon(pcnst, list_p, method_p, status_p)

@export
def gas_wetdep_readnl_codon(pcnst: int, list_p: cobj, method_p: cobj, status_p: cobj):
    return gas_wetdep_readnl_status_codon(pcnst, list_p, method_p, status_p)

@export
def tracer_cnst_init_codon() -> int:
    return _common.tracer_cnst_init_codon()

@export
def rate_diags_init_codon(
    tag_len: int,
    fieldname_len: int,
    tag_count: int,
    tag_ascii_p: cobj,
    rate_name_ascii_p: cobj,
) -> int:
    return _common.rate_diags_init_codon(tag_len, fieldname_len, tag_count, tag_ascii_p, rate_name_ascii_p)

@export
def noy_ubc_readnl_codon() -> int:
    return _common.noy_ubc_readnl_codon()

@export
def phtadj_codon() -> int:
    return 1

@export
def spedata_setopts_codon() -> int:
    return 1

@export
def spedata_defaultopts_codon(
    spe_data_file_len: int,
    spe_data_file_p: cobj,
    spe_remove_file: int,
    spe_filenames_list_len: int,
    spe_filenames_list_p: cobj,
    present_data_file: int,
    spe_data_file_out_p: cobj,
    present_remove_file: int,
    spe_remove_file_out_p: cobj,
    present_filenames_list: int,
    spe_filenames_list_out_p: cobj,
):
    return _common.spedata_defaultopts_codon(
        spe_data_file_len,
        spe_data_file_p,
        spe_remove_file,
        spe_filenames_list_len,
        spe_filenames_list_p,
        present_data_file,
        spe_data_file_out_p,
        present_remove_file,
        spe_remove_file_out_p,
        present_filenames_list,
        spe_filenames_list_out_p,
    )

@export
def aerodep_flx_prescribed_codon(active: int) -> int:
    return active

@export
def init_prescribed_ghg_restart_codon(stage: int) -> int:
    return stage

@export
def write_prescribed_ghg_restart_codon(stage: int) -> int:
    return stage

@export
def init_prescribed_ozone_restart_codon(stage: int) -> int:
    return stage

@export
def write_prescribed_ozone_restart_codon(stage: int) -> int:
    return stage

@export
def init_prescribed_volcaero_restart_codon(stage: int) -> int:
    return stage

@export
def write_prescribed_volcaero_restart_codon(stage: int) -> int:
    return stage

@export
def modal_aero_calcsize_reg_codon(stage: int) -> int:
    return stage

@export
def init_prescribed_aero_restart_codon(stage: int) -> int:
    return stage

@export
def write_prescribed_aero_restart_codon(stage: int) -> int:
    return stage

@export
def aircraft_emit_adv_codon(active: int) -> int:
    return active

@export
def aircraft_emit_readnl_codon(active: int) -> int:
    return active

@export
def flbc_chk_codon(active: int) -> int:
    return active

@export
def flbc_set_codon(flbc_count: int) -> int:
    if flbc_count < 1:
        return 1
    return 0

@export
def noy_ubc_active_codon(active: int) -> int:
    return _common.noy_ubc_active_codon(active)

@export
def noy_ubc_advance_codon(active: int) -> int:
    return noy_ubc_active_codon(active)

def spedata_active_codon(active: int) -> int:
    return _common.spedata_active_codon(active)

@export
def spedata_init_codon(active: int) -> int:
    return spedata_active_codon(active)

@export
def advance_spedata_codon(active: int) -> int:
    return spedata_active_codon(active)

@export
def lightning_inti_active_codon(no_ndx: int, xno_ndx: int) -> int:
    return _common.lightning_inti_active_codon(no_ndx, xno_ndx)

@export
def lightning_inti_codon(no_ndx: int, xno_ndx: int) -> int:
    return lightning_inti_active_codon(no_ndx, xno_ndx)

@export
def euvac_set_etf_active_codon(active: int) -> int:
    return _common.euvac_set_etf_active_codon(active)

@export
def euvac_set_etf_codon(active: int) -> int:
    return euvac_set_etf_active_codon(active)

@export
def spe_prod_zero_codon(active: int, ncol: int, pver: int, noxprod_p: cobj, hoxprod_p: cobj) -> int:
    return _common.spe_prod_zero_codon(active, ncol, pver, noxprod_p, hoxprod_p)

@export
def gcr_ionization_noxhox_zero_codon(active: int, ncol: int, pver: int, gcr_nox_p: cobj, gcr_hox_p: cobj) -> int:
    return _common.gcr_ionization_noxhox_zero_codon(active, ncol, pver, gcr_nox_p, gcr_hox_p)

@export
def gcr_ionization_noxhox_codon(active: int, ncol: int, pver: int, gcr_nox_p: cobj, gcr_hox_p: cobj) -> int:
    return gcr_ionization_noxhox_zero_codon(active, ncol, pver, gcr_nox_p, gcr_hox_p)

@export
def airpl_src_active_codon(active: int) -> int:
    return _common.airpl_src_active_codon(active)

@export
def airpl_src_codon(active: int) -> int:
    return airpl_src_active_codon(active)

@export
def airpl_set_zero_codon(active: int, ncol: int, pver: int, no_air_p: cobj, co_air_p: cobj) -> int:
    return _common.airpl_set_zero_codon(active, ncol, pver, no_air_p, co_air_p)

@export
def airpl_set_codon(active: int, ncol: int, pver: int, no_air_p: cobj, co_air_p: cobj) -> int:
    return _common.airpl_set_codon(active, ncol, pver, no_air_p, co_air_p)

@export
def sulf_inti_active_codon(active: int) -> int:
    return _common.sulf_inti_active_codon(active)

@export
def sulf_inti_codon(active: int) -> int:
    return sulf_inti_active_codon(active)

@export
def sox_inti_codon(active: int) -> int:
    return active

@export
def sox_cldaero_init_codon(active: int) -> int:
    return active

@export
def fstrat_inti_active_codon(active: int) -> int:
    return _common.fstrat_inti_active_codon(active)

@export
def fstrat_inti_codon(active: int) -> int:
    return fstrat_inti_active_codon(active)

def set_fstrat_vals_active_codon(active: int) -> int:
    return _common.set_fstrat_vals_active_codon(active)

@export
def set_fstrat_vals_codon(active: int) -> int:
    return set_fstrat_vals_active_codon(active)

def set_fstrat_h2o_active_codon(active: int) -> int:
    return _common.set_fstrat_h2o_active_codon(active)

@export
def set_fstrat_h2o_codon(active: int) -> int:
    return set_fstrat_h2o_active_codon(active)

@export
def jeuv_init_active_codon(active: int) -> int:
    return _common.jeuv_init_active_codon(active)

@export
def jeuv_init_codon(active: int) -> int:
    return jeuv_init_active_codon(active)

@export
def charge_fix_active_codon(active: int) -> int:
    return _common.charge_fix_active_codon(active)

@export
def charge_fix_codon(active: int) -> int:
    return charge_fix_active_codon(active)

@export
def o1d_to_2oh_adj_init_active_codon(active: int) -> int:
    return _common.o1d_to_2oh_adj_init_active_codon(active)

@export
def o1d_to_2oh_adj_init_codon(active: int) -> int:
    return o1d_to_2oh_adj_init_active_codon(active)

@export
def init_airglow_active_codon(active: int) -> int:
    return _common.init_airglow_active_codon(active)

@export
def init_airglow_codon(active: int) -> int:
    return init_airglow_active_codon(active)

@export
def register_cfc11star_active_codon(active: int) -> int:
    return _common.register_cfc11star_active_codon(active)

@export
def register_cfc11star_codon(active: int) -> int:
    return register_cfc11star_active_codon(active)

def update_cfc11star_active_codon(active: int) -> int:
    return _common.update_cfc11star_active_codon(active)

@export
def update_cfc11star_codon(active: int) -> int:
    return update_cfc11star_active_codon(active)

@export
def chlorine_loading_init_active_codon(active: int) -> int:
    return _common.chlorine_loading_init_active_codon(active)

@export
def chlorine_loading_init_codon(active: int) -> int:
    return chlorine_loading_init_active_codon(active)

@export
def parse_rate_sums_active_codon(active: int) -> int:
    return _common.parse_rate_sums_active_codon(active)

@export
def parse_rate_sums_codon(active: int) -> int:
    return parse_rate_sums_active_codon(active)

@export
def prescribed_aero_adv_codon(active: int) -> int:
    return active

@export
def prescribed_ghg_adv_codon(active: int) -> int:
    return active

@export
def aircraft_emit_register_codon(active: int) -> int:
    return active

@export
def prescribed_aero_register_codon(active: int) -> int:
    return active

@export
def prescribed_ghg_register_codon(active: int) -> int:
    return active

@export
def prescribed_strataero_register_codon(active: int) -> int:
    return active

@export
def prescribed_volcaero_adv_codon(active: int) -> int:
    return active

@export
def chem_final_codon() -> int:
    return _common.chem_final_codon()

@export
def gcr_ionization_init_codon() -> int:
    return _common.gcr_ionization_init_codon()

@export
def gcr_ionization_adv_codon() -> int:
    return _common.gcr_ionization_adv_codon()

@export
def aurora_timestep_init_codon() -> int:
    return _common.aurora_timestep_init_codon()

@export
def set_sulf_time_codon() -> int:
    return _common.set_sulf_time_codon()

@export
def neu_wetdep_init_active_codon(method_is_neu: int, wetdep_count: int) -> int:
    return _common.neu_wetdep_init_active_codon(method_is_neu, wetdep_count)

@export
def neu_wetdep_init_codon(method_is_neu: int, wetdep_count: int) -> int:
    return neu_wetdep_init_active_codon(method_is_neu, wetdep_count)

@export
def dvel_inti_fromlnd_codon() -> int:
    return 1

@export
def linoz_data_setopts_codon() -> int:
    return 1

@export
def tracer_cnst_setopts_codon() -> int:
    return 1

@export
def tracer_srcs_setopts_codon() -> int:
    return 1

@export
def chem_init_restart_codon() -> int:
    return 1

@export
def chem_write_restart_codon() -> int:
    return 1

@export
def chlorine_loading_advance_codon(active: int) -> int:
    return active

@export
def euvac_init_codon(active: int) -> int:
    return active

@export
def gcr_ionization_readnl_codon() -> int:
    return 1

@export
def lin_strat_chem_inti_codon(active: int) -> int:
    return active

@export
def init_linoz_data_restart_codon() -> int:
    return 1

@export
def write_linoz_data_restart_codon() -> int:
    return 1

@export
def chm_diags_zero_codon(
    ncol: int,
    pver: int,
    vmr_nox_p: cobj,
    vmr_noy_p: cobj,
    vmr_clox_p: cobj,
    vmr_cloy_p: cobj,
    vmr_tcly_p: cobj,
    vmr_brox_p: cobj,
    vmr_broy_p: cobj,
    vmr_toth_p: cobj,
    vmr_tbry_p: cobj,
    vmr_foy_p: cobj,
    vmr_tfy_p: cobj,
    mmr_noy_p: cobj,
    mmr_sox_p: cobj,
    mmr_nhx_p: cobj,
    df_noy_p: cobj,
    df_sox_p: cobj,
    df_nhx_p: cobj,
):
    return _diags.chm_diags_zero_codon(
        ncol,
        pver,
        vmr_nox_p,
        vmr_noy_p,
        vmr_clox_p,
        vmr_cloy_p,
        vmr_tcly_p,
        vmr_brox_p,
        vmr_broy_p,
        vmr_toth_p,
        vmr_tbry_p,
        vmr_foy_p,
        vmr_tfy_p,
        mmr_noy_p,
        mmr_sox_p,
        mmr_nhx_p,
        df_noy_p,
        df_sox_p,
        df_nhx_p,
    )

@export
def chm_diags_mass_codon(
    ncol: int,
    pver: int,
    rgrav: float,
    rearth: float,
    area_p: cobj,
    pdel_p: cobj,
    mass_p: cobj,
):
    return _diags.chm_diags_mass_codon(ncol, pver, rgrav, rearth, area_p, pdel_p, mass_p)

@export
def chm_diags_species_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    m: int,
    do3chm_flag: int,
    fillvalue: float,
    n_molwgt: float,
    s_molwgt: float,
    id_ch4: int,
    id_n2o5: int,
    id_cfc12: int,
    id_cl2: int,
    id_cl2o2: int,
    id_cfc114: int,
    id_hcfc141b: int,
    id_h1202: int,
    id_h2402: int,
    id_ch2br2: int,
    id_cfc11: int,
    id_cfc113: int,
    id_ch3ccl3: int,
    id_chbr3: int,
    id_ccl4: int,
    id_hcfc22: int,
    id_cf2clbr: int,
    id_hcfc142b: int,
    id_cof2: int,
    id_cf3br: int,
    id_cfc115: int,
    in_foy: int,
    in_tfy: int,
    in_nox: int,
    in_noy: int,
    in_sox: int,
    in_nhx: int,
    in_clox: int,
    in_cloy: int,
    in_tcly: int,
    in_brox: int,
    in_broy: int,
    in_tbry: int,
    in_toth: int,
    vmr_p: cobj,
    mmr_p: cobj,
    depflx_p: cobj,
    mmr_tend_p: cobj,
    mass_p: cobj,
    pmid_p: cobj,
    adv_mass_p: cobj,
    ltrop_p: cobj,
    vmr_nox_p: cobj,
    vmr_noy_p: cobj,
    vmr_clox_p: cobj,
    vmr_cloy_p: cobj,
    vmr_tcly_p: cobj,
    vmr_brox_p: cobj,
    vmr_broy_p: cobj,
    vmr_toth_p: cobj,
    vmr_tbry_p: cobj,
    vmr_foy_p: cobj,
    vmr_tfy_p: cobj,
    mmr_noy_p: cobj,
    mmr_sox_p: cobj,
    mmr_nhx_p: cobj,
    df_noy_p: cobj,
    df_sox_p: cobj,
    df_nhx_p: cobj,
    net_chem_p: cobj,
    do3chm_trp_p: cobj,
    do3chm_lms_p: cobj,
):
    return _diags.chm_diags_species_codon(
        ncol,
        pver,
        gas_pcnst,
        m,
        do3chm_flag,
        fillvalue,
        n_molwgt,
        s_molwgt,
        id_ch4,
        id_n2o5,
        id_cfc12,
        id_cl2,
        id_cl2o2,
        id_cfc114,
        id_hcfc141b,
        id_h1202,
        id_h2402,
        id_ch2br2,
        id_cfc11,
        id_cfc113,
        id_ch3ccl3,
        id_chbr3,
        id_ccl4,
        id_hcfc22,
        id_cf2clbr,
        id_hcfc142b,
        id_cof2,
        id_cf3br,
        id_cfc115,
        in_foy,
        in_tfy,
        in_nox,
        in_noy,
        in_sox,
        in_nhx,
        in_clox,
        in_cloy,
        in_tcly,
        in_brox,
        in_broy,
        in_tbry,
        in_toth,
        vmr_p,
        mmr_p,
        depflx_p,
        mmr_tend_p,
        mass_p,
        pmid_p,
        adv_mass_p,
        ltrop_p,
        vmr_nox_p,
        vmr_noy_p,
        vmr_clox_p,
        vmr_cloy_p,
        vmr_tcly_p,
        vmr_brox_p,
        vmr_broy_p,
        vmr_toth_p,
        vmr_tbry_p,
        vmr_foy_p,
        vmr_tfy_p,
        mmr_noy_p,
        mmr_sox_p,
        mmr_nhx_p,
        df_noy_p,
        df_sox_p,
        df_nhx_p,
        net_chem_p,
        do3chm_trp_p,
        do3chm_lms_p,
    )

@export
def chm_diags_species_packed_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    m: int,
    do3chm_flag: int,
    fillvalue: float,
    n_molwgt: float,
    s_molwgt: float,
    ids_p: cobj,
    flags_p: cobj,
    vmr_p: cobj,
    mmr_p: cobj,
    depflx_p: cobj,
    mmr_tend_p: cobj,
    mass_p: cobj,
    pmid_p: cobj,
    adv_mass_p: cobj,
    ltrop_p: cobj,
    vmr_nox_p: cobj,
    vmr_noy_p: cobj,
    vmr_clox_p: cobj,
    vmr_cloy_p: cobj,
    vmr_tcly_p: cobj,
    vmr_brox_p: cobj,
    vmr_broy_p: cobj,
    vmr_toth_p: cobj,
    vmr_tbry_p: cobj,
    vmr_foy_p: cobj,
    vmr_tfy_p: cobj,
    mmr_noy_p: cobj,
    mmr_sox_p: cobj,
    mmr_nhx_p: cobj,
    df_noy_p: cobj,
    df_sox_p: cobj,
    df_nhx_p: cobj,
    net_chem_p: cobj,
    do3chm_trp_p: cobj,
    do3chm_lms_p: cobj,
):
    return _diags.chm_diags_species_packed_codon(
        ncol,
        pver,
        gas_pcnst,
        m,
        do3chm_flag,
        fillvalue,
        n_molwgt,
        s_molwgt,
        ids_p,
        flags_p,
        vmr_p,
        mmr_p,
        depflx_p,
        mmr_tend_p,
        mass_p,
        pmid_p,
        adv_mass_p,
        ltrop_p,
        vmr_nox_p,
        vmr_noy_p,
        vmr_clox_p,
        vmr_cloy_p,
        vmr_tcly_p,
        vmr_brox_p,
        vmr_broy_p,
        vmr_toth_p,
        vmr_tbry_p,
        vmr_foy_p,
        vmr_tfy_p,
        mmr_noy_p,
        mmr_sox_p,
        mmr_nhx_p,
        df_noy_p,
        df_sox_p,
        df_nhx_p,
        net_chem_p,
        do3chm_trp_p,
        do3chm_lms_p,
    )

@export
def chm_diags_euv_codon(
    stage: int,
    ncol: int,
    pver: int,
    indexm: int,
    id_o: int,
    id_o2: int,
    id_h: int,
    id_n: int,
    id_no: int,
    rid_scalar: int,
    rid_jeuv_p: cobj,
    vmr_p: cobj,
    rxt_rates_p: cobj,
    invariants_p: cobj,
    wrk_p: cobj,
):
    return _diags.chm_diags_euv_codon(
        stage,
        ncol,
        pver,
        indexm,
        id_o,
        id_o2,
        id_h,
        id_n,
        id_no,
        rid_scalar,
        rid_jeuv_p,
        vmr_p,
        rxt_rates_p,
        invariants_p,
        wrk_p,
    )
@export
def usrrxt_inti_codon(stage: int) -> int:
    return stage

@export
def usrrxt_inti_has_ion_codon(
    ion1: int,
    ion2: int,
    ion3: int,
    elec1: int,
    elec2: int,
    elec3: int,
) -> int:
    return _common.usrrxt_inti_has_ion_codon(ion1, ion2, ion3, elec1, elec2, elec3)

@export
def comp_exp_codon(x_p: cobj, y_p: cobj, n: int):
    return _common.comp_exp_codon(x_p, y_p, n)

@export
def heatnirco2_init_codon(stage: int) -> int:
    return stage

@export
def heatnirco2_init_xspara_codon(ndpara: int, zppara_p: cobj, xspara_p: cobj):
    return _common.heatnirco2_init_xspara_codon(ndpara, zppara_p, xspara_p)

@export
def rebin_codon(
    nsrc: int,
    ntrg: int,
    src_x_p: cobj,
    trg_x_p: cobj,
    src_p: cobj,
    trg_p: cobj,
):
    return _photolysis.rebin_codon(
        nsrc,
        ntrg,
        src_x_p,
        trg_x_p,
        src_p,
        trg_p,
    )

@export
def sulf_interp_codon(
    ncol: int,
    pcols: int,
    pver: int,
    begchunk: int,
    lchnk: int,
    read_sulf_flag: int,
    fields_data_p: cobj,
    ccm_sulf_p: cobj,
):
    return _gas_phase.sulf_interp_codon(
        ncol,
        pcols,
        pver,
        begchunk,
        lchnk,
        read_sulf_flag,
        fields_data_p,
        ccm_sulf_p,
    )

@export
def sulf_interp_stage_dispatch_codon(
    ncol: int,
    pcols: int,
    pver: int,
    begchunk: int,
    lchnk: int,
    read_sulf_flag: int,
    fields_data_p: cobj,
    ccm_sulf_p: cobj,
):
    return sulf_interp_codon(
        ncol,
        pcols,
        pver,
        begchunk,
        lchnk,
        read_sulf_flag,
        fields_data_p,
        ccm_sulf_p,
    )

@export
def jlong_timestep_init_codon(
    jlong_used_flag: int,
    nsrc: int,
    ntrg: int,
    src_x_p: cobj,
    trg_x_p: cobj,
    src_p: cobj,
    trg_p: cobj,
):
    return _photolysis.jlong_timestep_init_codon(
        jlong_used_flag,
        nsrc,
        ntrg,
        src_x_p,
        trg_x_p,
        src_p,
        trg_p,
    )

@export
def jlong_init_set_we_codon(
    nw: int,
    wc_p: cobj,
    wlintv_p: cobj,
    we_p: cobj,
):
    return _photolysis.jlong_init_set_we_codon(
        nw,
        wc_p,
        wlintv_p,
        we_p,
    )

@export
def jlong_init_solar_batch_codon(
    data_nw: int,
    nw: int,
    data_we_p: cobj,
    wc_p: cobj,
    wlintv_p: cobj,
    we_p: cobj,
    data_etf_p: cobj,
    etfphot_p: cobj,
):
    return _photolysis.jlong_init_solar_batch_codon(
        data_nw,
        nw,
        data_we_p,
        wc_p,
        wlintv_p,
        we_p,
        data_etf_p,
        etfphot_p,
    )

@export
def jlong_get_xsqy_numj_codon(
    phtcnt: int,
    lng_indexer_p: cobj,
    numj_p: cobj,
):
    return _photolysis.jlong_get_xsqy_numj_codon(
        phtcnt,
        lng_indexer_p,
        numj_p,
    )

@export
def jlong_get_xsqy_read_order_codon(
    phtcnt: int,
    numj: int,
    lng_indexer_p: cobj,
    read_varids_p: cobj,
):
    return _photolysis.jlong_get_xsqy_read_order_codon(
        phtcnt,
        numj,
        lng_indexer_p,
        read_varids_p,
    )

@export
def jlong_get_xsqy_meta_batch_codon(
    phtcnt: int,
    lng_indexer_p: cobj,
    numj_p: cobj,
    read_varids_p: cobj,
):
    return _photolysis.jlong_get_xsqy_meta_batch_codon(
        phtcnt,
        lng_indexer_p,
        numj_p,
        read_varids_p,
    )

@export
def jlong_get_xsqy_index_map_codon(
    phtcnt: int,
    lng_indexer_p: cobj,
    wrk_ndx_p: cobj,
):
    return _photolysis.jlong_get_xsqy_index_map_codon(
        phtcnt,
        lng_indexer_p,
        wrk_ndx_p,
    )

@export
def jlong_get_xsqy_dprs_codon(
    np_xs: int,
    prs_p: cobj,
    dprs_p: cobj,
):
    return _photolysis.jlong_get_xsqy_dprs_codon(
        np_xs,
        prs_p,
        dprs_p,
    )

@export
def jlong_get_rsf_scale_codon(
    nw: int,
    nump: int,
    numsza: int,
    numcolo3: int,
    numalb: int,
    wlintv_p: cobj,
    rsf_tab_p: cobj,
):
    return _photolysis.jlong_get_rsf_scale_codon(
        nw,
        nump,
        numsza,
        numcolo3,
        numalb,
        wlintv_p,
        rsf_tab_p,
    )

@export
def jlong_get_rsf_deltas_codon(
    nump: int,
    numsza: int,
    numalb: int,
    numcolo3: int,
    p_p: cobj,
    sza_p: cobj,
    alb_p: cobj,
    o3rat_p: cobj,
    del_p_p: cobj,
    del_sza_p: cobj,
    del_alb_p: cobj,
    del_o3rat_p: cobj,
):
    return _photolysis.jlong_get_rsf_deltas_codon(
        nump,
        numsza,
        numalb,
        numcolo3,
        p_p,
        sza_p,
        alb_p,
        o3rat_p,
        del_p_p,
        del_sza_p,
        del_alb_p,
        del_o3rat_p,
    )

@export
def jlong_get_rsf_bde_codon(
    nw: int,
    use_bde_flag: int,
    hc_val: float,
    wc_o2_b_val: float,
    wc_o3_a_val: float,
    wc_o3_b_val: float,
    wc_p: cobj,
    bde_o2_b_p: cobj,
    bde_o3_a_p: cobj,
    bde_o3_b_p: cobj,
):
    return _photolysis.jlong_get_rsf_bde_codon(
        nw,
        use_bde_flag,
        hc_val,
        wc_o2_b_val,
        wc_o3_a_val,
        wc_o3_b_val,
        wc_p,
        bde_o2_b_p,
        bde_o3_a_p,
        bde_o3_b_p,
    )

@export
def jlong_get_rsf_postread_batch_codon(
    nw: int,
    nump: int,
    numsza: int,
    numalb: int,
    numcolo3: int,
    use_bde_flag: int,
    hc_val: float,
    wc_o2_b_val: float,
    wc_o3_a_val: float,
    wc_o3_b_val: float,
    wc_p: cobj,
    p_p: cobj,
    sza_p: cobj,
    alb_p: cobj,
    o3rat_p: cobj,
    bde_o2_b_p: cobj,
    bde_o3_a_p: cobj,
    bde_o3_b_p: cobj,
    del_p_p: cobj,
    del_sza_p: cobj,
    del_alb_p: cobj,
    del_o3rat_p: cobj,
):
    return _photolysis.jlong_get_rsf_postread_batch_codon(
        nw,
        nump,
        numsza,
        numalb,
        numcolo3,
        use_bde_flag,
        hc_val,
        wc_o2_b_val,
        wc_o3_a_val,
        wc_o3_b_val,
        wc_p,
        p_p,
        sza_p,
        alb_p,
        o3rat_p,
        bde_o2_b_p,
        bde_o3_a_p,
        bde_o3_b_p,
        del_p_p,
        del_sza_p,
        del_alb_p,
        del_o3rat_p,
    )

@export
def jlong_prep_init_solar_batch_codon(
    data_nw: int,
    nw: int,
    data_we_p: cobj,
    wc_p: cobj,
    wlintv_p: cobj,
    we_p: cobj,
    data_etf_p: cobj,
    etfphot_p: cobj,
):
    return _photolysis.jlong_prep_init_solar_batch_codon(
        data_nw,
        nw,
        data_we_p,
        wc_p,
        wlintv_p,
        we_p,
        data_etf_p,
        etfphot_p,
    )

@export
def jlong_prep_get_xsqy_meta_batch_codon(
    phtcnt: int,
    lng_indexer_p: cobj,
    numj_p: cobj,
    read_varids_p: cobj,
):
    return _photolysis.jlong_prep_get_xsqy_meta_batch_codon(
        phtcnt,
        lng_indexer_p,
        numj_p,
        read_varids_p,
    )

@export
def jlong_prep_get_xsqy_dprs_codon(
    np_xs: int,
    prs_p: cobj,
    dprs_p: cobj,
):
    return _photolysis.jlong_prep_get_xsqy_dprs_codon(
        np_xs,
        prs_p,
        dprs_p,
    )

@export
def jlong_prep_get_rsf_scale_codon(
    nw: int,
    nump: int,
    numsza: int,
    numcolo3: int,
    numalb: int,
    wlintv_p: cobj,
    rsf_tab_p: cobj,
):
    return _photolysis.jlong_prep_get_rsf_scale_codon(
        nw,
        nump,
        numsza,
        numcolo3,
        numalb,
        wlintv_p,
        rsf_tab_p,
    )

@export
def jlong_prep_get_rsf_postread_batch_codon(
    nw: int,
    nump: int,
    numsza: int,
    numalb: int,
    numcolo3: int,
    use_bde_flag: int,
    hc_val: float,
    wc_o2_b_val: float,
    wc_o3_a_val: float,
    wc_o3_b_val: float,
    wc_p: cobj,
    p_p: cobj,
    sza_p: cobj,
    alb_p: cobj,
    o3rat_p: cobj,
    bde_o2_b_p: cobj,
    bde_o3_a_p: cobj,
    bde_o3_b_p: cobj,
    del_p_p: cobj,
    del_sza_p: cobj,
    del_alb_p: cobj,
    del_o3rat_p: cobj,
):
    return _photolysis.jlong_prep_get_rsf_postread_batch_codon(
        nw,
        nump,
        numsza,
        numalb,
        numcolo3,
        use_bde_flag,
        hc_val,
        wc_o2_b_val,
        wc_o3_a_val,
        wc_o3_b_val,
        wc_p,
        p_p,
        sza_p,
        alb_p,
        o3rat_p,
        bde_o2_b_p,
        bde_o3_a_p,
        bde_o3_b_p,
        del_p_p,
        del_sza_p,
        del_alb_p,
        del_o3rat_p,
    )

@export
def jlong_prep_batch_codon(
    stage: int,
    data_nw: int,
    nw: int,
    phtcnt: int,
    np_xs: int,
    nump: int,
    numsza: int,
    numalb: int,
    numcolo3: int,
    use_bde_flag: int,
    hc_val: float,
    wc_o2_b_val: float,
    wc_o3_a_val: float,
    wc_o3_b_val: float,
    data_we_p: cobj,
    wc_p: cobj,
    wlintv_p: cobj,
    we_p: cobj,
    data_etf_p: cobj,
    etfphot_p: cobj,
    lng_indexer_p: cobj,
    numj_p: cobj,
    read_varids_p: cobj,
    prs_p: cobj,
    dprs_p: cobj,
    rsf_tab_p: cobj,
    p_p: cobj,
    sza_p: cobj,
    alb_p: cobj,
    o3rat_p: cobj,
    bde_o2_b_p: cobj,
    bde_o3_a_p: cobj,
    bde_o3_b_p: cobj,
    del_p_p: cobj,
    del_sza_p: cobj,
    del_alb_p: cobj,
    del_o3rat_p: cobj,
):
    return _photolysis.jlong_prep_batch_codon(
        stage,
        data_nw,
        nw,
        phtcnt,
        np_xs,
        nump,
        numsza,
        numalb,
        numcolo3,
        use_bde_flag,
        hc_val,
        wc_o2_b_val,
        wc_o3_a_val,
        wc_o3_b_val,
        data_we_p,
        wc_p,
        wlintv_p,
        we_p,
        data_etf_p,
        etfphot_p,
        lng_indexer_p,
        numj_p,
        read_varids_p,
        prs_p,
        dprs_p,
        rsf_tab_p,
        p_p,
        sza_p,
        alb_p,
        o3rat_p,
        bde_o2_b_p,
        bde_o3_a_p,
        bde_o3_b_p,
        del_p_p,
        del_sza_p,
        del_alb_p,
        del_o3rat_p,
    )

@export
def jlong_prep_stage_dispatch_codon(
    stage: int,
    data_nw: int,
    nw: int,
    phtcnt: int,
    np_xs: int,
    nump: int,
    numsza: int,
    numalb: int,
    numcolo3: int,
    use_bde_flag: int,
    hc_val: float,
    wc_o2_b_val: float,
    wc_o3_a_val: float,
    wc_o3_b_val: float,
    data_we_p: cobj,
    wc_p: cobj,
    wlintv_p: cobj,
    we_p: cobj,
    data_etf_p: cobj,
    etfphot_p: cobj,
    lng_indexer_p: cobj,
    numj_p: cobj,
    read_varids_p: cobj,
    prs_p: cobj,
    dprs_p: cobj,
    rsf_tab_p: cobj,
    p_p: cobj,
    sza_p: cobj,
    alb_p: cobj,
    o3rat_p: cobj,
    bde_o2_b_p: cobj,
    bde_o3_a_p: cobj,
    bde_o3_b_p: cobj,
    del_p_p: cobj,
    del_sza_p: cobj,
    del_alb_p: cobj,
    del_o3rat_p: cobj,
):
    return jlong_prep_batch_codon(
        stage,
        data_nw,
        nw,
        phtcnt,
        np_xs,
        nump,
        numsza,
        numalb,
        numcolo3,
        use_bde_flag,
        hc_val,
        wc_o2_b_val,
        wc_o3_a_val,
        wc_o3_b_val,
        data_we_p,
        wc_p,
        wlintv_p,
        we_p,
        data_etf_p,
        etfphot_p,
        lng_indexer_p,
        numj_p,
        read_varids_p,
        prs_p,
        dprs_p,
        rsf_tab_p,
        p_p,
        sza_p,
        alb_p,
        o3rat_p,
        bde_o2_b_p,
        bde_o3_a_p,
        bde_o3_b_p,
        del_p_p,
        del_sza_p,
        del_alb_p,
        del_o3rat_p,
    )

@export
def jlong_init_codon(
    stage: int,
    data_nw: int,
    nw: int,
    phtcnt: int,
    np_xs: int,
    nump: int,
    numsza: int,
    numalb: int,
    numcolo3: int,
    use_bde_flag: int,
    hc_val: float,
    wc_o2_b_val: float,
    wc_o3_a_val: float,
    wc_o3_b_val: float,
    data_we_p: cobj,
    wc_p: cobj,
    wlintv_p: cobj,
    we_p: cobj,
    data_etf_p: cobj,
    etfphot_p: cobj,
    lng_indexer_p: cobj,
    numj_p: cobj,
    read_varids_p: cobj,
    prs_p: cobj,
    dprs_p: cobj,
    rsf_tab_p: cobj,
    p_p: cobj,
    sza_p: cobj,
    alb_p: cobj,
    o3rat_p: cobj,
    bde_o2_b_p: cobj,
    bde_o3_a_p: cobj,
    bde_o3_b_p: cobj,
    del_p_p: cobj,
    del_sza_p: cobj,
    del_alb_p: cobj,
    del_o3rat_p: cobj,
):
    return jlong_prep_stage_dispatch_codon(
        stage,
        data_nw,
        nw,
        phtcnt,
        np_xs,
        nump,
        numsza,
        numalb,
        numcolo3,
        use_bde_flag,
        hc_val,
        wc_o2_b_val,
        wc_o3_a_val,
        wc_o3_b_val,
        data_we_p,
        wc_p,
        wlintv_p,
        we_p,
        data_etf_p,
        etfphot_p,
        lng_indexer_p,
        numj_p,
        read_varids_p,
        prs_p,
        dprs_p,
        rsf_tab_p,
        p_p,
        sza_p,
        alb_p,
        o3rat_p,
        bde_o2_b_p,
        bde_o3_a_p,
        bde_o3_b_p,
        del_p_p,
        del_sza_p,
        del_alb_p,
        del_o3rat_p,
    )

@export
def zenith_codon(
    ncol: int,
    calday: float,
    pi_val: float,
    delta: float,
    clat_p: cobj,
    clon_p: cobj,
    coszrs_p: cobj,
):
    return _photolysis.zenith_codon(
        ncol,
        calday,
        pi_val,
        delta,
        clat_p,
        clon_p,
        coszrs_p,
    )

@export
def jlong_interpolate_rsf_codon(
    nw: int,
    nump: int,
    numsza: int,
    numalb: int,
    numcolo3: int,
    kbot: int,
    sza_in: float,
    alb_in_p: cobj,
    p_in_p: cobj,
    colo3_in_p: cobj,
    p_grid_p: cobj,
    del_p_p: cobj,
    sza_grid_p: cobj,
    del_sza_p: cobj,
    alb_grid_p: cobj,
    del_alb_p: cobj,
    o3rat_p: cobj,
    del_o3rat_p: cobj,
    colo3_grid_p: cobj,
    rsf_tab_p: cobj,
    etfphot_p: cobj,
    psum_l_p: cobj,
    rsf_p: cobj,
):
    return _photolysis.jlong_interpolate_rsf_codon(
        nw,
        nump,
        numsza,
        numalb,
        numcolo3,
        kbot,
        sza_in,
        alb_in_p,
        p_in_p,
        colo3_in_p,
        p_grid_p,
        del_p_p,
        sza_grid_p,
        del_sza_p,
        alb_grid_p,
        del_alb_p,
        o3rat_p,
        del_o3rat_p,
        colo3_grid_p,
        rsf_tab_p,
        etfphot_p,
        psum_l_p,
        rsf_p,
    )

@export
def interpolate_rsf_codon(
    nw: int,
    nump: int,
    numsza: int,
    numalb: int,
    numcolo3: int,
    kbot: int,
    sza_in: float,
    alb_in_p: cobj,
    p_in_p: cobj,
    colo3_in_p: cobj,
    p_grid_p: cobj,
    del_p_p: cobj,
    sza_grid_p: cobj,
    del_sza_p: cobj,
    alb_grid_p: cobj,
    del_alb_p: cobj,
    o3rat_p: cobj,
    del_o3rat_p: cobj,
    colo3_grid_p: cobj,
    rsf_tab_p: cobj,
    etfphot_p: cobj,
    psum_l_p: cobj,
    rsf_p: cobj,
):
    return jlong_interpolate_rsf_codon(
        nw,
        nump,
        numsza,
        numalb,
        numcolo3,
        kbot,
        sza_in,
        alb_in_p,
        p_in_p,
        colo3_in_p,
        p_grid_p,
        del_p_p,
        sza_grid_p,
        del_sza_p,
        alb_grid_p,
        del_alb_p,
        o3rat_p,
        del_o3rat_p,
        colo3_grid_p,
        rsf_tab_p,
        etfphot_p,
        psum_l_p,
        rsf_p,
    )

@export
def jlong_photo_fill_xswk_codon(
    numj: int,
    nw: int,
    nt: int,
    np_xs: int,
    k: int,
    p_in_p: cobj,
    t_in_p: cobj,
    prs_p: cobj,
    dprs_p: cobj,
    xsqy_p: cobj,
    xswk_p: cobj,
):
    return _photolysis.jlong_photo_fill_xswk_codon(
        numj,
        nw,
        nt,
        np_xs,
        k,
        p_in_p,
        t_in_p,
        prs_p,
        dprs_p,
        xsqy_p,
        xswk_p,
    )

@export
def jlong_photo_accum_codon(
    numj: int,
    nw: int,
    xswk_p: cobj,
    rsf_col_p: cobj,
    j_long_col_p: cobj,
):
    return _photolysis.jlong_photo_accum_codon(
        numj,
        nw,
        xswk_p,
        rsf_col_p,
        j_long_col_p,
    )

@export
def jlong_photo_loop_codon(
    numj: int,
    nw: int,
    nt: int,
    np_xs: int,
    nlev: int,
    p_in_p: cobj,
    t_in_p: cobj,
    prs_p: cobj,
    dprs_p: cobj,
    xsqy_p: cobj,
    rsf_p: cobj,
    xswk_p: cobj,
    j_long_p: cobj,
):
    return _photolysis.jlong_photo_loop_codon(
        numj,
        nw,
        nt,
        np_xs,
        nlev,
        p_in_p,
        t_in_p,
        prs_p,
        dprs_p,
        xsqy_p,
        rsf_p,
        xswk_p,
        j_long_p,
    )

@export
def jlong_photo_codon(
    numj: int,
    nw: int,
    nt: int,
    np_xs: int,
    nump: int,
    numsza: int,
    numalb: int,
    numcolo3: int,
    nlev: int,
    sza_in: float,
    alb_in_p: cobj,
    p_in_p: cobj,
    t_in_p: cobj,
    colo3_in_p: cobj,
    p_grid_p: cobj,
    del_p_p: cobj,
    sza_grid_p: cobj,
    del_sza_p: cobj,
    alb_grid_p: cobj,
    del_alb_p: cobj,
    o3rat_p: cobj,
    del_o3rat_p: cobj,
    colo3_grid_p: cobj,
    rsf_tab_p: cobj,
    etfphot_p: cobj,
    prs_p: cobj,
    dprs_p: cobj,
    xsqy_p: cobj,
    rsf_p: cobj,
    psum_l_p: cobj,
    xswk_p: cobj,
    j_long_p: cobj,
):
    return _photolysis.jlong_photo_codon(
        numj,
        nw,
        nt,
        np_xs,
        nump,
        numsza,
        numalb,
        numcolo3,
        nlev,
        sza_in,
        alb_in_p,
        p_in_p,
        t_in_p,
        colo3_in_p,
        p_grid_p,
        del_p_p,
        sza_grid_p,
        del_sza_p,
        alb_grid_p,
        del_alb_p,
        o3rat_p,
        del_o3rat_p,
        colo3_grid_p,
        rsf_tab_p,
        etfphot_p,
        prs_p,
        dprs_p,
        xsqy_p,
        rsf_p,
        psum_l_p,
        xswk_p,
        j_long_p,
    )

@export
def chem_emissions_zero_cflx_codon(
    pcols: int,
    pcnst: int,
    map2chm_p: cobj,
    cflx_p: cobj,
):
    return _emissions.chem_emissions_zero_cflx_codon(
        pcols,
        pcnst,
        map2chm_p,
        cflx_p,
    )

@export
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
    return _emissions.chem_emissions_shell_codon(
        stage,
        ncol,
        pcols,
        pcnst,
        nmegan,
        h2o_ndx,
        map2chm_p,
        megan_indices_map_p,
        megan_wght_factors_p,
        meganflx_p,
        cflx_p,
        megflx_p,
        sflx_p,
    )

@export
def chem_emissions_codon(
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
    return chem_emissions_shell_codon(
        stage,
        ncol,
        pcols,
        pcnst,
        nmegan,
        h2o_ndx,
        map2chm_p,
        megan_indices_map_p,
        megan_wght_factors_p,
        meganflx_p,
        cflx_p,
        megflx_p,
        sflx_p,
    )

@export
def gas_phase_chemdr_prepare_sza_codon(
    ncol: int,
    rad2deg: float,
    zen_angle_p: cobj,
    sza_p: cobj,
):
    return _gas_phase.gas_phase_chemdr_prepare_sza_codon(
        ncol,
        rad2deg,
        zen_angle_p,
        sza_p,
    )

@export
def table_photo_zero_photos_codon(
    ncol: int,
    pver: int,
    phtcnt: int,
    photos_p: cobj,
):
    return _photolysis.table_photo_zero_photos_codon(
        ncol,
        pver,
        phtcnt,
        photos_p,
    )

@export
def table_photo_zero_finalize_batch_codon(
    stage: int,
    ncol: int,
    pver: int,
    phtcnt: int,
    photos_p: cobj,
    jno2a_ndx: int,
    jno2_ndx: int,
    jn2o5a_ndx: int,
    jn2o5_ndx: int,
    jn2o5b_ndx: int,
    jhno3a_ndx: int,
    jhno3_ndx: int,
    jno3a_ndx: int,
    jno3_ndx: int,
    jho2no2a_ndx: int,
    jho2no2_ndx: int,
    jmpana_ndx: int,
    jmpan_ndx: int,
    jpana_ndx: int,
    jpan_ndx: int,
    jonitra_ndx: int,
    jonitr_ndx: int,
    jo1da_ndx: int,
    jo1d_ndx: int,
    jo3pa_ndx: int,
    jo3p_ndx: int,
):
    return _photolysis.table_photo_zero_finalize_batch_codon(
        stage,
        ncol,
        pver,
        phtcnt,
        photos_p,
        jno2a_ndx,
        jno2_ndx,
        jn2o5a_ndx,
        jn2o5_ndx,
        jn2o5b_ndx,
        jhno3a_ndx,
        jhno3_ndx,
        jno3a_ndx,
        jno3_ndx,
        jho2no2a_ndx,
        jho2no2_ndx,
        jmpana_ndx,
        jmpan_ndx,
        jpana_ndx,
        jpan_ndx,
        jonitra_ndx,
        jonitr_ndx,
        jo1da_ndx,
        jo1d_ndx,
        jo3pa_ndx,
        jo3p_ndx,
    )

@export
def table_photo_daylight_prepare_batch_codon(
    ncol: int,
    pcols: int,
    pver: int,
    ncol_abs: int,
    nfs: int,
    gas_pcnst: int,
    i_col: int,
    p1: int,
    p2: int,
    do_jshort_flag: int,
    ptop_gt_10_flag: int,
    o_is_inv_flag: int,
    o2_is_inv_flag: int,
    o3_is_inv_flag: int,
    n2_is_inv_flag: int,
    no_is_inv_flag: int,
    o_ndx: int,
    o2_ndx: int,
    o3_inv_ndx: int,
    o3_ndx: int,
    n2_ndx: int,
    no_ndx: int,
    indexm: int,
    pa2mb: float,
    pmid_p: cobj,
    pdel_p: cobj,
    col_dens_p: cobj,
    lwc_p: cobj,
    clouds_p: cobj,
    temper_p: cobj,
    zmid_p: cobj,
    zint_p: cobj,
    vmr_p: cobj,
    invariants_p: cobj,
    parg_p: cobj,
    colo3_p: cobj,
    fac1_p: cobj,
    lwc_line_p: cobj,
    cld_line_p: cobj,
    tline_p: cobj,
    zarg_p: cobj,
    o_den_p: cobj,
    o2_den_p: cobj,
    o3_den_p: cobj,
    no_den_p: cobj,
    n2_den_p: cobj,
):
    return _photolysis.table_photo_daylight_prepare_batch_codon(
        ncol,
        pcols,
        pver,
        ncol_abs,
        nfs,
        gas_pcnst,
        i_col,
        p1,
        p2,
        do_jshort_flag,
        ptop_gt_10_flag,
        o_is_inv_flag,
        o2_is_inv_flag,
        o3_is_inv_flag,
        n2_is_inv_flag,
        no_is_inv_flag,
        o_ndx,
        o2_ndx,
        o3_inv_ndx,
        o3_ndx,
        n2_ndx,
        no_ndx,
        indexm,
        pa2mb,
        pmid_p,
        pdel_p,
        col_dens_p,
        lwc_p,
        clouds_p,
        temper_p,
        zmid_p,
        zint_p,
        vmr_p,
        invariants_p,
        parg_p,
        colo3_p,
        fac1_p,
        lwc_line_p,
        cld_line_p,
        tline_p,
        zarg_p,
        o_den_p,
        o2_den_p,
        o3_den_p,
        no_den_p,
        n2_den_p,
    )

@export
def table_photo_daylight_setup_codon(
    ncol: int,
    pcols: int,
    pver: int,
    ncol_abs: int,
    i_col: int,
    p1: int,
    p2: int,
    pa2mb: float,
    pmid_p: cobj,
    pdel_p: cobj,
    col_dens_p: cobj,
    lwc_p: cobj,
    clouds_p: cobj,
    temper_p: cobj,
    zmid_p: cobj,
    parg_p: cobj,
    colo3_p: cobj,
    fac1_p: cobj,
    lwc_line_p: cobj,
    cld_line_p: cobj,
    tline_p: cobj,
    zarg_p: cobj,
):
    return _photolysis.table_photo_daylight_setup_codon(
        ncol,
        pcols,
        pver,
        ncol_abs,
        i_col,
        p1,
        p2,
        pa2mb,
        pmid_p,
        pdel_p,
        col_dens_p,
        lwc_p,
        clouds_p,
        temper_p,
        zmid_p,
        parg_p,
        colo3_p,
        fac1_p,
        lwc_line_p,
        cld_line_p,
        tline_p,
        zarg_p,
    )

@export
def table_photo_scale_cld_mult_codon(
    pver: int,
    esfact: float,
    cld_mult_p: cobj,
):
    return _photolysis.table_photo_scale_cld_mult_codon(
        pver,
        esfact,
        cld_mult_p,
    )

@export
def mmr2vmr_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    mbar_p: cobj,
    mmr_p: cobj,
    adv_mass_p: cobj,
    vmr_p: cobj,
):
    return _gas_phase.mmr2vmr_codon(
        ncol,
        pcols,
        pver,
        gas_pcnst,
        mbar_p,
        mmr_p,
        adv_mass_p,
        vmr_p,
    )

@export
def vmr2mmr_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    mbar_p: cobj,
    vmr_p: cobj,
    adv_mass_p: cobj,
    mmr_p: cobj,
):
    return _gas_phase.vmr2mmr_codon(
        ncol,
        pcols,
        pver,
        gas_pcnst,
        mbar_p,
        vmr_p,
        adv_mass_p,
        mmr_p,
    )

@export
def h2o_to_vmr_codon(
    ncol: int,
    pcols: int,
    pver: int,
    adv_mass_h2o: float,
    h2o_mmr_p: cobj,
    mbar_p: cobj,
    h2o_vmr_p: cobj,
):
    return _gas_phase.h2o_to_vmr_codon(
        ncol,
        pcols,
        pver,
        adv_mass_h2o,
        h2o_mmr_p,
        mbar_p,
        h2o_vmr_p,
    )

@export
def init_mass_xforms_codon(list_len: int, list_ascii_p: cobj, list_count: int, adv_mass_p: cobj) -> float:
    list_ascii = Ptr[int](list_ascii_p)
    adv_mass = Ptr[float](adv_mass_p)

    for item in range(list_count):
        offset = item * list_len
        name_trim = list_len
        while name_trim > 0 and list_ascii[offset + name_trim - 1] == 32:
            name_trim -= 1
        if (
            name_trim == 3
            and list_ascii[offset] == 72
            and list_ascii[offset + 1] == 50
            and list_ascii[offset + 2] == 79
        ):
            return adv_mass[item]

    return 18.0

@export
def mass_xforms_batch_mmr2vmr_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    mbar_p: cobj,
    mmr_p: cobj,
    adv_mass_p: cobj,
    vmr_p: cobj,
):
    return _gas_phase.mass_xforms_batch_mmr2vmr_codon(
        ncol,
        pcols,
        pver,
        gas_pcnst,
        mbar_p,
        mmr_p,
        adv_mass_p,
        vmr_p,
    )

@export
def mass_xforms_batch_vmr2mmr_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    mbar_p: cobj,
    vmr_p: cobj,
    adv_mass_p: cobj,
    mmr_p: cobj,
):
    return _gas_phase.mass_xforms_batch_vmr2mmr_codon(
        ncol,
        pcols,
        pver,
        gas_pcnst,
        mbar_p,
        vmr_p,
        adv_mass_p,
        mmr_p,
    )

@export
def mass_xforms_batch_h2o_to_vmr_codon(
    ncol: int,
    pcols: int,
    pver: int,
    adv_mass_h2o: float,
    h2o_mmr_p: cobj,
    mbar_p: cobj,
    h2o_vmr_p: cobj,
):
    return _gas_phase.mass_xforms_batch_h2o_to_vmr_codon(
        ncol,
        pcols,
        pver,
        adv_mass_h2o,
        h2o_mmr_p,
        mbar_p,
        h2o_vmr_p,
    )

@export
def set_mean_mass_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    id_o2: int,
    id_o: int,
    id_h: int,
    id_n: int,
    fixed_mbar: int,
    mwdry: float,
    mmr_p: cobj,
    adv_mass_p: cobj,
    mbar_p: cobj,
):
    return _gas_phase.set_mean_mass_codon(
        ncol,
        pcols,
        pver,
        gas_pcnst,
        id_o2,
        id_o,
        id_h,
        id_n,
        fixed_mbar,
        mwdry,
        mmr_p,
        adv_mass_p,
        mbar_p,
    )

@export
def setinv_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    nfs: int,
    m_ndx: int,
    n2_ndx: int,
    o2_ndx: int,
    h2o_ndx: int,
    id_o: int,
    id_o2: int,
    id_h: int,
    has_n2: int,
    has_o2: int,
    has_h2o: int,
    has_var_o2: int,
    pa_xfac: float,
    boltz_cgs: float,
    tfld_p: cobj,
    h2ovmr_p: cobj,
    vmr_p: cobj,
    pmid_p: cobj,
    invariants_p: cobj,
):
    return _gas_phase.setinv_codon(
        ncol,
        pcols,
        pver,
        gas_pcnst,
        nfs,
        m_ndx,
        n2_ndx,
        o2_ndx,
        h2o_ndx,
        id_o,
        id_o2,
        id_h,
        has_n2,
        has_o2,
        has_h2o,
        has_var_o2,
        pa_xfac,
        boltz_cgs,
        tfld_p,
        h2ovmr_p,
        vmr_p,
        pmid_p,
        invariants_p,
    )

@export
def setinv_apply_tracer_cnst_codon(
    ncol: int,
    pver: int,
    nfs: int,
    ndx: int,
    m_ndx: int,
    cnst_offline_p: cobj,
    invariants_p: cobj,
):
    return _gas_phase.setinv_apply_tracer_cnst_codon(
        ncol,
        pver,
        nfs,
        ndx,
        m_ndx,
        cnst_offline_p,
        invariants_p,
    )

@export
def setinv_apply_tracer_cnst_stage_dispatch_codon(
    ncol: int,
    pver: int,
    nfs: int,
    ndx: int,
    m_ndx: int,
    cnst_offline_p: cobj,
    invariants_p: cobj,
):
    return setinv_apply_tracer_cnst_codon(
        ncol,
        pver,
        nfs,
        ndx,
        m_ndx,
        cnst_offline_p,
        invariants_p,
    )

@export
def setinv_copy_invariant_codon(
    ncol: int,
    pver: int,
    nfs: int,
    inv_ndx: int,
    invariants_p: cobj,
    tmp_out_p: cobj,
):
    return _gas_phase.setinv_copy_invariant_codon(
        ncol,
        pver,
        nfs,
        inv_ndx,
        invariants_p,
        tmp_out_p,
    )

@export
def setinv_vmr_output_codon(
    ncol: int,
    pver: int,
    nfs: int,
    inv_ndx: int,
    m_ndx: int,
    invariants_p: cobj,
    tmp_out_p: cobj,
):
    return _gas_phase.setinv_vmr_output_codon(
        ncol,
        pver,
        nfs,
        inv_ndx,
        m_ndx,
        invariants_p,
        tmp_out_p,
    )

@export
def setinv_output_pair_codon(
    ncol: int,
    pver: int,
    nfs: int,
    inv_ndx: int,
    m_ndx: int,
    invariants_p: cobj,
    tmp_dens_p: cobj,
    tmp_vmr_p: cobj,
):
    return _gas_phase.setinv_output_pair_codon(
        ncol,
        pver,
        nfs,
        inv_ndx,
        m_ndx,
        invariants_p,
        tmp_dens_p,
        tmp_vmr_p,
    )

@export
def setinv_output_pair_stage_dispatch_codon(
    ncol: int,
    pver: int,
    nfs: int,
    inv_ndx: int,
    m_ndx: int,
    invariants_p: cobj,
    tmp_dens_p: cobj,
    tmp_vmr_p: cobj,
):
    return setinv_output_pair_codon(
        ncol,
        pver,
        nfs,
        inv_ndx,
        m_ndx,
        invariants_p,
        tmp_dens_p,
        tmp_vmr_p,
    )

@export
def charge_balance_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    np_ndx: int,
    n2p_ndx: int,
    op_ndx: int,
    o2p_ndx: int,
    nop_ndx: int,
    conc_p: cobj,
    wrk_p: cobj,
):
    return _gas_phase.charge_balance_codon(
        ncol,
        pver,
        gas_pcnst,
        np_ndx,
        n2p_ndx,
        op_ndx,
        o2p_ndx,
        nop_ndx,
        conc_p,
        wrk_p,
    )

@export
def setcol_codon(
    ncol: int,
    pver: int,
    ncol_abs: int,
    col_delta_p: cobj,
    col_dens_p: cobj,
):
    return _gas_phase.setcol_codon(
        ncol,
        pver,
        ncol_abs,
        col_delta_p,
        col_dens_p,
    )

@export
def set_ub_col_codon(
    ncol: int,
    pcols: int,
    pver: int,
    ncol_abs: int,
    indexm: int,
    o3rad_ndx: int,
    ox_ndx: int,
    o3_ndx: int,
    o3_inv_ndx: int,
    o2_ndx: int,
    o2_is_inv: int,
    xfactor: float,
    pdel_p: cobj,
    vmr_p: cobj,
    invariants_p: cobj,
    o2_exo_col_p: cobj,
    o3_exo_col_p: cobj,
    col_delta_p: cobj,
):
    return _gas_phase.set_ub_col_codon(
        ncol,
        pcols,
        pver,
        ncol_abs,
        indexm,
        o3rad_ndx,
        ox_ndx,
        o3_ndx,
        o3_inv_ndx,
        o2_ndx,
        o2_is_inv,
        xfactor,
        pdel_p,
        vmr_p,
        invariants_p,
        o2_exo_col_p,
        o3_exo_col_p,
        col_delta_p,
    )

@export
def cloud_mod_codon(
    pver: int,
    zen_angle: float,
    srf_alb: float,
    rgrav: float,
    clouds_p: cobj,
    lwc_p: cobj,
    delp_p: cobj,
    eff_alb_p: cobj,
    cld_mult_p: cobj,
    del_lwp_p: cobj,
    del_tau_p: cobj,
    above_tau_p: cobj,
    below_tau_p: cobj,
    above_cld_p: cobj,
    below_cld_p: cobj,
    above_tra_p: cobj,
    below_tra_p: cobj,
    fac1_p: cobj,
    fac2_p: cobj,
):
    return _gas_phase.cloud_mod_codon(
        pver,
        zen_angle,
        srf_alb,
        rgrav,
        clouds_p,
        lwc_p,
        delp_p,
        eff_alb_p,
        cld_mult_p,
        del_lwp_p,
        del_tau_p,
        above_tau_p,
        below_tau_p,
        above_cld_p,
        below_cld_p,
        above_tra_p,
        below_tra_p,
        fac1_p,
        fac2_p,
    )

@export
def photo_inti_fixed_press_setup_codon(
    pinterp: float,
    n_exo_levs: int,
    levs_p: cobj,
    ki_p: cobj,
    delp_p: cobj,
):
    return _photolysis.photo_inti_fixed_press_setup_codon(
        pinterp,
        n_exo_levs,
        levs_p,
        ki_p,
        delp_p,
    )

@export
def photo_timestep_init_exo_time_codon(
    calday: float,
    days_p: cobj,
    next_p: cobj,
    last_p: cobj,
    dels_p: cobj,
):
    return _photolysis.photo_timestep_init_exo_time_codon(
        calday,
        days_p,
        next_p,
        last_p,
        dels_p,
    )

@export
def photo_prep_batch_codon(
    stage: int,
    pinterp: float,
    calday: float,
    n_exo_levs: int,
    levs_p: cobj,
    days_p: cobj,
    ki_p: cobj,
    next_p: cobj,
    last_p: cobj,
    delp_p: cobj,
    dels_p: cobj,
):
    return _photolysis.photo_prep_batch_codon(
        stage,
        pinterp,
        calday,
        n_exo_levs,
        levs_p,
        days_p,
        ki_p,
        next_p,
        last_p,
        delp_p,
        dels_p,
    )

@export
def photo_prep_stage_dispatch_codon(
    stage: int,
    pinterp: float,
    calday: float,
    n_exo_levs: int,
    levs_p: cobj,
    days_p: cobj,
    ki_p: cobj,
    next_p: cobj,
    last_p: cobj,
    delp_p: cobj,
    dels_p: cobj,
):
    return photo_prep_batch_codon(
        stage,
        pinterp,
        calday,
        n_exo_levs,
        levs_p,
        days_p,
        ki_p,
        next_p,
        last_p,
        delp_p,
        dels_p,
    )

@export
def photo_timestep_init_codon(
    stage: int,
    pinterp: float,
    calday: float,
    n_exo_levs: int,
    levs_p: cobj,
    days_p: cobj,
    ki_p: cobj,
    next_p: cobj,
    last_p: cobj,
    delp_p: cobj,
    dels_p: cobj,
):
    return photo_prep_stage_dispatch_codon(
        stage,
        pinterp,
        calday,
        n_exo_levs,
        levs_p,
        days_p,
        ki_p,
        next_p,
        last_p,
        delp_p,
        dels_p,
    )

@export
def photo_prep_fixed_press_setup_codon(
    pinterp: float,
    n_exo_levs: int,
    levs_p: cobj,
    ki_p: cobj,
    delp_p: cobj,
):
    return _photolysis.photo_prep_fixed_press_setup_codon(
        pinterp,
        n_exo_levs,
        levs_p,
        ki_p,
        delp_p,
    )

@export
def photo_prep_timestep_init_exo_time_codon(
    calday: float,
    days_p: cobj,
    next_p: cobj,
    last_p: cobj,
    dels_p: cobj,
):
    return _photolysis.photo_prep_timestep_init_exo_time_codon(
        calday,
        days_p,
        next_p,
        last_p,
        dels_p,
    )

@export
def table_photo_jlong_apply_codon(
    ncol: int,
    pver: int,
    phtcnt: int,
    nlng: int,
    i_col: int,
    photos_p: cobj,
    lng_prates_p: cobj,
    cld_mult_p: cobj,
    lng_indexer_p: cobj,
    alias_mult2_p: cobj,
):
    return _photolysis.table_photo_jlong_apply_codon(
        ncol,
        pver,
        phtcnt,
        nlng,
        i_col,
        photos_p,
        lng_prates_p,
        cld_mult_p,
        lng_indexer_p,
        alias_mult2_p,
    )

@export
def table_photo_jno_ho2no2_codon(
    ncol: int,
    pver: int,
    phtcnt: int,
    i_col: int,
    jno_ndx: int,
    jho2no2_ndx: int,
    do_jshort: int,
    has_o2_col: int,
    has_o3_col: int,
    zen_angle: float,
    photos_p: cobj,
    col_dens_p: cobj,
    cld_mult_p: cobj,
):
    return _photolysis.table_photo_jno_ho2no2_codon(
        ncol,
        pver,
        phtcnt,
        i_col,
        jno_ndx,
        jho2no2_ndx,
        do_jshort,
        has_o2_col,
        has_o3_col,
        zen_angle,
        photos_p,
        col_dens_p,
        cld_mult_p,
    )

@export
def table_photo_postcloud_batch_codon(
    ncol: int,
    pver: int,
    phtcnt: int,
    nlng: int,
    i_col: int,
    jno_ndx: int,
    jho2no2_ndx: int,
    do_jshort_flag: int,
    has_o2_col_flag: int,
    has_o3_col_flag: int,
    zen_angle: float,
    esfact: float,
    photos_p: cobj,
    lng_prates_p: cobj,
    cld_mult_p: cobj,
    col_dens_p: cobj,
    lng_indexer_p: cobj,
    alias_mult2_p: cobj,
):
    return _photolysis.table_photo_postcloud_batch_codon(
        ncol,
        pver,
        phtcnt,
        nlng,
        i_col,
        jno_ndx,
        jho2no2_ndx,
        do_jshort_flag,
        has_o2_col_flag,
        has_o3_col_flag,
        zen_angle,
        esfact,
        photos_p,
        lng_prates_p,
        cld_mult_p,
        col_dens_p,
        lng_indexer_p,
        alias_mult2_p,
    )

@export
def table_photo_direct_batch_codon(
    stage: int,
    ncol: int,
    pcols: int,
    pver: int,
    ncol_abs: int,
    nfs: int,
    gas_pcnst: int,
    phtcnt: int,
    nlng: int,
    i_col: int,
    p1: int,
    p2: int,
    do_jshort_flag: int,
    ptop_gt_10_flag: int,
    o_is_inv_flag: int,
    o2_is_inv_flag: int,
    o3_is_inv_flag: int,
    n2_is_inv_flag: int,
    no_is_inv_flag: int,
    o_ndx: int,
    o2_ndx: int,
    o3_inv_ndx: int,
    o3_ndx: int,
    n2_ndx: int,
    no_ndx: int,
    indexm: int,
    jno_ndx: int,
    jho2no2_ndx: int,
    has_o2_col_flag: int,
    has_o3_col_flag: int,
    pa2mb: float,
    zen_angle: float,
    srf_alb: float,
    rgrav: float,
    esfact: float,
    pmid_p: cobj,
    pdel_p: cobj,
    col_dens_p: cobj,
    lwc_p: cobj,
    clouds_p: cobj,
    temper_p: cobj,
    zmid_p: cobj,
    zint_p: cobj,
    vmr_p: cobj,
    invariants_p: cobj,
    parg_p: cobj,
    colo3_p: cobj,
    fac1_p: cobj,
    lwc_line_p: cobj,
    cld_line_p: cobj,
    tline_p: cobj,
    zarg_p: cobj,
    o_den_p: cobj,
    o2_den_p: cobj,
    o3_den_p: cobj,
    no_den_p: cobj,
    n2_den_p: cobj,
    eff_alb_p: cobj,
    cld_mult_p: cobj,
    del_lwp_p: cobj,
    del_tau_p: cobj,
    above_tau_p: cobj,
    below_tau_p: cobj,
    above_cld_p: cobj,
    below_cld_p: cobj,
    above_tra_p: cobj,
    below_tra_p: cobj,
    cloud_fac1_p: cobj,
    cloud_fac2_p: cobj,
    photos_p: cobj,
    lng_prates_p: cobj,
    lng_indexer_p: cobj,
    alias_mult2_p: cobj,
    jno2a_ndx: int,
    jno2_ndx: int,
    jn2o5a_ndx: int,
    jn2o5_ndx: int,
    jn2o5b_ndx: int,
    jhno3a_ndx: int,
    jhno3_ndx: int,
    jno3a_ndx: int,
    jno3_ndx: int,
    jho2no2a_ndx: int,
    jho2no2_finalize_ndx: int,
    jmpana_ndx: int,
    jmpan_ndx: int,
    jpana_ndx: int,
    jpan_ndx: int,
    jonitra_ndx: int,
    jonitr_ndx: int,
    jo1da_ndx: int,
    jo1d_ndx: int,
    jo3pa_ndx: int,
    jo3p_ndx: int,
):
    return _photolysis.table_photo_direct_batch_codon(
        stage,
        ncol,
        pcols,
        pver,
        ncol_abs,
        nfs,
        gas_pcnst,
        phtcnt,
        nlng,
        i_col,
        p1,
        p2,
        do_jshort_flag,
        ptop_gt_10_flag,
        o_is_inv_flag,
        o2_is_inv_flag,
        o3_is_inv_flag,
        n2_is_inv_flag,
        no_is_inv_flag,
        o_ndx,
        o2_ndx,
        o3_inv_ndx,
        o3_ndx,
        n2_ndx,
        no_ndx,
        indexm,
        jno_ndx,
        jho2no2_ndx,
        has_o2_col_flag,
        has_o3_col_flag,
        pa2mb,
        zen_angle,
        srf_alb,
        rgrav,
        esfact,
        pmid_p,
        pdel_p,
        col_dens_p,
        lwc_p,
        clouds_p,
        temper_p,
        zmid_p,
        zint_p,
        vmr_p,
        invariants_p,
        parg_p,
        colo3_p,
        fac1_p,
        lwc_line_p,
        cld_line_p,
        tline_p,
        zarg_p,
        o_den_p,
        o2_den_p,
        o3_den_p,
        no_den_p,
        n2_den_p,
        eff_alb_p,
        cld_mult_p,
        del_lwp_p,
        del_tau_p,
        above_tau_p,
        below_tau_p,
        above_cld_p,
        below_cld_p,
        above_tra_p,
        below_tra_p,
        cloud_fac1_p,
        cloud_fac2_p,
        photos_p,
        lng_prates_p,
        lng_indexer_p,
        alias_mult2_p,
        jno2a_ndx,
        jno2_ndx,
        jn2o5a_ndx,
        jn2o5_ndx,
        jn2o5b_ndx,
        jhno3a_ndx,
        jhno3_ndx,
        jno3a_ndx,
        jno3_ndx,
        jho2no2a_ndx,
        jho2no2_finalize_ndx,
        jmpana_ndx,
        jmpan_ndx,
        jpana_ndx,
        jpan_ndx,
        jonitra_ndx,
        jonitr_ndx,
        jo1da_ndx,
        jo1d_ndx,
        jo3pa_ndx,
        jo3p_ndx,
    )

@export
def table_photo_stage_dispatch_codon(
    stage: int,
    ncol: int,
    pcols: int,
    pver: int,
    ncol_abs: int,
    nfs: int,
    gas_pcnst: int,
    phtcnt: int,
    nlng: int,
    i_col: int,
    p1: int,
    p2: int,
    do_jshort_flag: int,
    ptop_gt_10_flag: int,
    o_is_inv_flag: int,
    o2_is_inv_flag: int,
    o3_is_inv_flag: int,
    n2_is_inv_flag: int,
    no_is_inv_flag: int,
    o_ndx: int,
    o2_ndx: int,
    o3_inv_ndx: int,
    o3_ndx: int,
    n2_ndx: int,
    no_ndx: int,
    indexm: int,
    jno_ndx: int,
    jho2no2_ndx: int,
    has_o2_col_flag: int,
    has_o3_col_flag: int,
    pa2mb: float,
    zen_angle: float,
    srf_alb: float,
    rgrav: float,
    esfact: float,
    pmid_p: cobj,
    pdel_p: cobj,
    col_dens_p: cobj,
    lwc_p: cobj,
    clouds_p: cobj,
    temper_p: cobj,
    zmid_p: cobj,
    zint_p: cobj,
    vmr_p: cobj,
    invariants_p: cobj,
    parg_p: cobj,
    colo3_p: cobj,
    fac1_p: cobj,
    lwc_line_p: cobj,
    cld_line_p: cobj,
    tline_p: cobj,
    zarg_p: cobj,
    o_den_p: cobj,
    o2_den_p: cobj,
    o3_den_p: cobj,
    no_den_p: cobj,
    n2_den_p: cobj,
    eff_alb_p: cobj,
    cld_mult_p: cobj,
    del_lwp_p: cobj,
    del_tau_p: cobj,
    above_tau_p: cobj,
    below_tau_p: cobj,
    above_cld_p: cobj,
    below_cld_p: cobj,
    above_tra_p: cobj,
    below_tra_p: cobj,
    cloud_fac1_p: cobj,
    cloud_fac2_p: cobj,
    photos_p: cobj,
    lng_prates_p: cobj,
    lng_indexer_p: cobj,
    alias_mult2_p: cobj,
    jno2a_ndx: int,
    jno2_ndx: int,
    jn2o5a_ndx: int,
    jn2o5_ndx: int,
    jn2o5b_ndx: int,
    jhno3a_ndx: int,
    jhno3_ndx: int,
    jno3a_ndx: int,
    jno3_ndx: int,
    jho2no2a_ndx: int,
    jho2no2_finalize_ndx: int,
    jmpana_ndx: int,
    jmpan_ndx: int,
    jpana_ndx: int,
    jpan_ndx: int,
    jonitra_ndx: int,
    jonitr_ndx: int,
    jo1da_ndx: int,
    jo1d_ndx: int,
    jo3pa_ndx: int,
    jo3p_ndx: int,
):
    return table_photo_direct_batch_codon(
        stage,
        ncol,
        pcols,
        pver,
        ncol_abs,
        nfs,
        gas_pcnst,
        phtcnt,
        nlng,
        i_col,
        p1,
        p2,
        do_jshort_flag,
        ptop_gt_10_flag,
        o_is_inv_flag,
        o2_is_inv_flag,
        o3_is_inv_flag,
        n2_is_inv_flag,
        no_is_inv_flag,
        o_ndx,
        o2_ndx,
        o3_inv_ndx,
        o3_ndx,
        n2_ndx,
        no_ndx,
        indexm,
        jno_ndx,
        jho2no2_ndx,
        has_o2_col_flag,
        has_o3_col_flag,
        pa2mb,
        zen_angle,
        srf_alb,
        rgrav,
        esfact,
        pmid_p,
        pdel_p,
        col_dens_p,
        lwc_p,
        clouds_p,
        temper_p,
        zmid_p,
        zint_p,
        vmr_p,
        invariants_p,
        parg_p,
        colo3_p,
        fac1_p,
        lwc_line_p,
        cld_line_p,
        tline_p,
        zarg_p,
        o_den_p,
        o2_den_p,
        o3_den_p,
        no_den_p,
        n2_den_p,
        eff_alb_p,
        cld_mult_p,
        del_lwp_p,
        del_tau_p,
        above_tau_p,
        below_tau_p,
        above_cld_p,
        below_cld_p,
        above_tra_p,
        below_tra_p,
        cloud_fac1_p,
        cloud_fac2_p,
        photos_p,
        lng_prates_p,
        lng_indexer_p,
        alias_mult2_p,
        jno2a_ndx,
        jno2_ndx,
        jn2o5a_ndx,
        jn2o5_ndx,
        jn2o5b_ndx,
        jhno3a_ndx,
        jhno3_ndx,
        jno3a_ndx,
        jno3_ndx,
        jho2no2a_ndx,
        jho2no2_finalize_ndx,
        jmpana_ndx,
        jmpan_ndx,
        jpana_ndx,
        jpan_ndx,
        jonitra_ndx,
        jonitr_ndx,
        jo1da_ndx,
        jo1d_ndx,
        jo3pa_ndx,
        jo3p_ndx,
    )

@export
def table_photo_prejlong_batch_codon(
    ncol: int,
    pcols: int,
    pver: int,
    ncol_abs: int,
    nfs: int,
    gas_pcnst: int,
    i_col: int,
    p1: int,
    p2: int,
    do_jshort_flag: int,
    ptop_gt_10_flag: int,
    o_is_inv_flag: int,
    o2_is_inv_flag: int,
    o3_is_inv_flag: int,
    n2_is_inv_flag: int,
    no_is_inv_flag: int,
    o_ndx: int,
    o2_ndx: int,
    o3_inv_ndx: int,
    o3_ndx: int,
    n2_ndx: int,
    no_ndx: int,
    indexm: int,
    pa2mb: float,
    zen_angle: float,
    srf_alb: float,
    rgrav: float,
    pmid_p: cobj,
    pdel_p: cobj,
    col_dens_p: cobj,
    lwc_p: cobj,
    clouds_p: cobj,
    temper_p: cobj,
    zmid_p: cobj,
    zint_p: cobj,
    vmr_p: cobj,
    invariants_p: cobj,
    parg_p: cobj,
    colo3_p: cobj,
    fac1_p: cobj,
    lwc_line_p: cobj,
    cld_line_p: cobj,
    tline_p: cobj,
    zarg_p: cobj,
    o_den_p: cobj,
    o2_den_p: cobj,
    o3_den_p: cobj,
    no_den_p: cobj,
    n2_den_p: cobj,
    eff_alb_p: cobj,
    cld_mult_p: cobj,
    del_lwp_p: cobj,
    del_tau_p: cobj,
    above_tau_p: cobj,
    below_tau_p: cobj,
    above_cld_p: cobj,
    below_cld_p: cobj,
    above_tra_p: cobj,
    below_tra_p: cobj,
    cloud_fac1_p: cobj,
    cloud_fac2_p: cobj,
):
    return _photolysis.table_photo_prejlong_batch_codon(
        ncol,
        pcols,
        pver,
        ncol_abs,
        nfs,
        gas_pcnst,
        i_col,
        p1,
        p2,
        do_jshort_flag,
        ptop_gt_10_flag,
        o_is_inv_flag,
        o2_is_inv_flag,
        o3_is_inv_flag,
        n2_is_inv_flag,
        no_is_inv_flag,
        o_ndx,
        o2_ndx,
        o3_inv_ndx,
        o3_ndx,
        n2_ndx,
        no_ndx,
        indexm,
        pa2mb,
        zen_angle,
        srf_alb,
        rgrav,
        pmid_p,
        pdel_p,
        col_dens_p,
        lwc_p,
        clouds_p,
        temper_p,
        zmid_p,
        zint_p,
        vmr_p,
        invariants_p,
        parg_p,
        colo3_p,
        fac1_p,
        lwc_line_p,
        cld_line_p,
        tline_p,
        zarg_p,
        o_den_p,
        o2_den_p,
        o3_den_p,
        no_den_p,
        n2_den_p,
        eff_alb_p,
        cld_mult_p,
        del_lwp_p,
        del_tau_p,
        above_tau_p,
        below_tau_p,
        above_cld_p,
        below_cld_p,
        above_tra_p,
        below_tra_p,
        cloud_fac1_p,
        cloud_fac2_p,
    )

@export
def table_photo_cloud_mod_batch_codon(
    pver: int,
    zen_angle: float,
    srf_alb: float,
    rgrav: float,
    clouds_p: cobj,
    lwc_p: cobj,
    delp_p: cobj,
    eff_alb_p: cobj,
    cld_mult_p: cobj,
    del_lwp_p: cobj,
    del_tau_p: cobj,
    above_tau_p: cobj,
    below_tau_p: cobj,
    above_cld_p: cobj,
    below_cld_p: cobj,
    above_tra_p: cobj,
    below_tra_p: cobj,
    fac1_p: cobj,
    fac2_p: cobj,
):
    return _photolysis.table_photo_cloud_mod_batch_codon(
        pver,
        zen_angle,
        srf_alb,
        rgrav,
        clouds_p,
        lwc_p,
        delp_p,
        eff_alb_p,
        cld_mult_p,
        del_lwp_p,
        del_tau_p,
        above_tau_p,
        below_tau_p,
        above_cld_p,
        below_cld_p,
        above_tra_p,
        below_tra_p,
        fac1_p,
        fac2_p,
    )

@export
def gas_phase_chemdr_zero_sulfate_codon(
    ncol: int,
    pver: int,
    sulfate_p: cobj,
):
    return _gas_phase.gas_phase_chemdr_zero_sulfate_codon(
        ncol,
        pver,
        sulfate_p,
    )

@export
def gas_phase_chemdr_load_prognostic_sulfate_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    so4_ndx: int,
    vmr_p: cobj,
    sulfate_p: cobj,
):
    return _gas_phase.gas_phase_chemdr_load_prognostic_sulfate_codon(
        ncol,
        pver,
        gas_pcnst,
        so4_ndx,
        vmr_p,
        sulfate_p,
    )

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
    return _emissions.chem_emissions_megan_flux_codon(
        ncol,
        pcols,
        megan_index,
        megan_weight,
        meganflx_p,
        cflx_p,
        megflx_p,
    )

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
    return _emissions.chem_emissions_add_sflx_codon(
        ncol,
        pcols,
        pcnst,
        h2o_ndx,
        map2chm_p,
        cflx_p,
        sflx_p,
    )

@export
def aero_model_gasaerexch_codon(
    stage: int,
    stage3_mode: int,
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    ndx_h2so4: int,
    delt: float,
    gravit: float,
    vmr0_p: cobj,
    vmr_p: cobj,
    vmrcw_p: cobj,
    dvmrdt_p: cobj,
    dvmrcwdt_p: cobj,
    mbar_p: cobj,
    pdel_p: cobj,
    adv_mass_p: cobj,
    wrk_p: cobj,
    del_h2so4_aeruptk_p: cobj,
):
    return _aero_bridge.aero_model_gasaerexch_codon(
        stage,
        stage3_mode,
        ncol,
        pcols,
        pver,
        gas_pcnst,
        ndx_h2so4,
        delt,
        gravit,
        vmr0_p,
        vmr_p,
        vmrcw_p,
        dvmrdt_p,
        dvmrcwdt_p,
        mbar_p,
        pdel_p,
        adv_mass_p,
        wrk_p,
        del_h2so4_aeruptk_p,
    )

@export
def aero_model_gasaerexch_stage_dispatch_codon(
    stage: int,
    stage3_mode: int,
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    ndx_h2so4: int,
    delt: float,
    gravit: float,
    vmr0_p: cobj,
    vmr_p: cobj,
    vmrcw_p: cobj,
    dvmrdt_p: cobj,
    dvmrcwdt_p: cobj,
    mbar_p: cobj,
    pdel_p: cobj,
    adv_mass_p: cobj,
    wrk_p: cobj,
    del_h2so4_aeruptk_p: cobj,
):
    return aero_model_gasaerexch_codon(
        stage,
        stage3_mode,
        ncol,
        pcols,
        pver,
        gas_pcnst,
        ndx_h2so4,
        delt,
        gravit,
        vmr0_p,
        vmr_p,
        vmrcw_p,
        dvmrdt_p,
        dvmrcwdt_p,
        mbar_p,
        pdel_p,
        adv_mass_p,
        wrk_p,
        del_h2so4_aeruptk_p,
    )

@export
def aero_model_gasaerexch_presetsox_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    delt: float,
    gravit: float,
    vmr0_p: cobj,
    vmr_p: cobj,
    dvmrdt_p: cobj,
    mbar_p: cobj,
    pdel_p: cobj,
    adv_mass_p: cobj,
    wrk_p: cobj,
):
    return _aero_bridge.aero_model_gasaerexch_presetsox_shell_codon(
        ncol,
        pcols,
        pver,
        gas_pcnst,
        delt,
        gravit,
        vmr0_p,
        vmr_p,
        dvmrdt_p,
        mbar_p,
        pdel_p,
        adv_mass_p,
        wrk_p,
    )

@export
def aero_model_gasaerexch_column_flux_codon(
    ncol: int,
    pcols: int,
    pver: int,
    adv_mass: float,
    gravit: float,
    field_p: cobj,
    mbar_p: cobj,
    pdel_p: cobj,
    wrk_p: cobj,
):
    return _aero_bridge.aero_model_gasaerexch_column_flux_codon(
        ncol,
        pcols,
        pver,
        adv_mass,
        gravit,
        field_p,
        mbar_p,
        pdel_p,
        wrk_p,
    )

@export
def aero_model_gasaerexch_h2so4_save_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    ndx_h2so4: int,
    vmr_p: cobj,
    del_h2so4_aeruptk_p: cobj,
):
    return _aero_bridge.aero_model_gasaerexch_h2so4_save_codon(
        ncol,
        pver,
        gas_pcnst,
        ndx_h2so4,
        vmr_p,
        del_h2so4_aeruptk_p,
    )

@export
def aero_model_gasaerexch_h2so4_delta_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    ndx_h2so4: int,
    vmr_p: cobj,
    del_h2so4_aeruptk_p: cobj,
):
    return _aero_bridge.aero_model_gasaerexch_h2so4_delta_codon(
        ncol,
        pver,
        gas_pcnst,
        ndx_h2so4,
        vmr_p,
        del_h2so4_aeruptk_p,
    )

@export
def aero_model_gasaerexch_gas_tend_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    delt: float,
    vmr0_p: cobj,
    vmr_p: cobj,
    dvmrdt_p: cobj,
):
    return _aero_bridge.aero_model_gasaerexch_gas_tend_codon(
        ncol,
        pver,
        gas_pcnst,
        delt,
        vmr0_p,
        vmr_p,
        dvmrdt_p,
    )

@export
def aero_model_gasaerexch_aq_tend_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    delt: float,
    vmr_p: cobj,
    vmrcw_p: cobj,
    dvmrdt_p: cobj,
    dvmrcwdt_p: cobj,
):
    return _aero_bridge.aero_model_gasaerexch_aq_tend_codon(
        ncol,
        pver,
        gas_pcnst,
        delt,
        vmr_p,
        vmrcw_p,
        dvmrdt_p,
        dvmrcwdt_p,
    )

@export
def aero_model_gasaerexch_vmrcw_batch_codon(
    mode: int,
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    qqcw_offset: int,
    mbar_ld1: int,
    qqcw_ptrs_p: cobj,
    qqcw_present_p: cobj,
    mbar_p: cobj,
    adv_mass_p: cobj,
    vmr_p: cobj,
):
    return _aero_bridge.aero_model_gasaerexch_vmrcw_batch_codon(
        mode,
        ncol,
        pcols,
        pver,
        gas_pcnst,
        qqcw_offset,
        mbar_ld1,
        qqcw_ptrs_p,
        qqcw_present_p,
        mbar_p,
        adv_mass_p,
        vmr_p,
    )

@export
def aero_model_gasaerexch_load_snapshot_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    qqcw_offset: int,
    mbar_ld1: int,
    qqcw_ptrs_p: cobj,
    qqcw_present_p: cobj,
    mbar_p: cobj,
    adv_mass_p: cobj,
    vmr_p: cobj,
    vmrcw_p: cobj,
    dvmrdt_p: cobj,
    dvmrcwdt_p: cobj,
):
    return _aero_bridge.aero_model_gasaerexch_load_snapshot_codon(
        ncol,
        pcols,
        pver,
        gas_pcnst,
        qqcw_offset,
        mbar_ld1,
        qqcw_ptrs_p,
        qqcw_present_p,
        mbar_p,
        adv_mass_p,
        vmr_p,
        vmrcw_p,
        dvmrdt_p,
        dvmrcwdt_p,
    )

@export
def aero_model_gasaerexch_preset_load_stage_codon(
    stage: int,
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    qqcw_offset: int,
    mbar_ld1: int,
    delt: float,
    gravit: float,
    qqcw_ptrs_p: cobj,
    qqcw_present_p: cobj,
    vmr0_p: cobj,
    vmr_p: cobj,
    vmrcw_p: cobj,
    dvmrdt_p: cobj,
    dvmrcwdt_p: cobj,
    mbar_p: cobj,
    pdel_p: cobj,
    adv_mass_p: cobj,
    wrk_p: cobj,
):
    return _aero_bridge.aero_model_gasaerexch_preset_load_stage_codon(
        stage,
        ncol,
        pcols,
        pver,
        gas_pcnst,
        qqcw_offset,
        mbar_ld1,
        delt,
        gravit,
        qqcw_ptrs_p,
        qqcw_present_p,
        vmr0_p,
        vmr_p,
        vmrcw_p,
        dvmrdt_p,
        dvmrcwdt_p,
        mbar_p,
        pdel_p,
        adv_mass_p,
        wrk_p,
    )

@export
def aero_model_gasaerexch_preset_load_stage_dispatch_codon(
    stage: int,
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    qqcw_offset: int,
    mbar_ld1: int,
    delt: float,
    gravit: float,
    qqcw_ptrs_p: cobj,
    qqcw_present_p: cobj,
    vmr0_p: cobj,
    vmr_p: cobj,
    vmrcw_p: cobj,
    dvmrdt_p: cobj,
    dvmrcwdt_p: cobj,
    mbar_p: cobj,
    pdel_p: cobj,
    adv_mass_p: cobj,
    wrk_p: cobj,
):
    return aero_model_gasaerexch_preset_load_stage_codon(
        stage,
        ncol,
        pcols,
        pver,
        gas_pcnst,
        qqcw_offset,
        mbar_ld1,
        delt,
        gravit,
        qqcw_ptrs_p,
        qqcw_present_p,
        vmr0_p,
        vmr_p,
        vmrcw_p,
        dvmrdt_p,
        dvmrcwdt_p,
        mbar_p,
        pdel_p,
        adv_mass_p,
        wrk_p,
    )

@export
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
    return _wetchem.neu_wetdep_aux_prepare_codon(
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

@export
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
    return _wetchem.neu_wetdep_aux_finish_codon(
        ncol,
        pcols,
        pver,
        pcnst,
        gas_cnt,
        delt,
        pi,
        mapping_to_mmr_p,
        lats_p,
        pmid_p,
        mass_in_layer_p,
        trc_mass_p,
        dtwr_p,
        wd_mmr_p,
        wd_tend_p,
    )

@export
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
    return _wetchem.neu_wetdep_henry_flags_codon(
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

@export
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
    return _wetchem.neu_wetdep_prepare_henry_shell_codon(
        ncol,
        pcols,
        pver,
        pcnst,
        gas_cnt,
        index_cldice,
        index_cldliq,
        nh3_ndx,
        co2_ndx,
        gravit,
        t0,
        ph,
        ph_inv,
        mapping_to_mmr_p,
        mapping_to_heff_p,
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
        dheff_p,
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
        heff_p,
        wrk_p,
        dk1s_p,
        dk2s_p,
        tckaqb_p,
    )

@export
def neu_wetdep_prepare_henry_shell_stage_dispatch_codon(
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
    return neu_wetdep_prepare_henry_shell_codon(
        ncol,
        pcols,
        pver,
        pcnst,
        gas_cnt,
        index_cldice,
        index_cldliq,
        nh3_ndx,
        co2_ndx,
        gravit,
        t0,
        ph,
        ph_inv,
        mapping_to_mmr_p,
        mapping_to_heff_p,
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
        dheff_p,
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
        heff_p,
        wrk_p,
        dk1s_p,
        dk2s_p,
        tckaqb_p,
    )

@export
def dempirical_codon(
    cwater: float,
    rrate: float,
    dempirical_p: cobj,
):
    return _wetchem.dempirical_codon(
        cwater,
        rrate,
        dempirical_p,
    )

@export
def disgas_codon(
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
    return _wetchem.disgas_codon(
        clwx,
        cfx,
        molmass,
        hstar,
        tm,
        pr,
        qm,
        qt,
        qtdis_p,
    )

@export
def raingas_codon(
    rrain: float,
    dtscav: float,
    clwx: float,
    cfx: float,
    qm: float,
    qt: float,
    qtdis: float,
    qtrain_p: cobj,
):
    return _wetchem.raingas_codon(
        rrain,
        dtscav,
        clwx,
        cfx,
        qm,
        qt,
        qtdis,
        qtrain_p,
    )

@export
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
    return _wetchem.neu_wetdep_washgas_codon(
        rwash,
        boxf,
        dtscav,
        qtrtop,
        hstar,
        tm,
        pr,
        qm,
        qt,
        qtwash_p,
        qtevap_p,
    )

@export
def washo_codon(
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
    return _wetchem.washo_codon(
        lpar,
        ntrace,
        hno3_ndx,
        do_diag,
        dempirical_impl,
        dtscav,
        garea,
        adj_factor,
        qttjfl_p,
        qm_p,
        pofl_p,
        delz_p,
        rls_p,
        clwc_p,
        ciwc_p,
        cfr_p,
        tem_p,
        evaprate_p,
        hstar_p,
        tcmass_p,
        tckaqb_p,
        tcnion_p,
        qt_rain_p,
        qt_rime_p,
        qt_wash_p,
        qt_evap_p,
        cfxx_p,
        qtt_p,
        qttnew_p,
    )

@export
def neu_wetdep_washo_columns_codon(
    ncol: int,
    lpar: int,
    ntrace: int,
    hno3_ndx: int,
    do_diag: int,
    dempirical_impl: int,
    dtscav: float,
    adj_factor: float,
    trc_mass_p: cobj,
    qm_p: cobj,
    pofl_p: cobj,
    delz_p: cobj,
    rls_p: cobj,
    clwc_p: cobj,
    ciwc_p: cobj,
    cfr_p: cobj,
    tem_p: cobj,
    evaprate_p: cobj,
    garea_p: cobj,
    hstar_p: cobj,
    tcmass_p: cobj,
    tckaqb_p: cobj,
    tcnion_p: cobj,
    qt_rain_p: cobj,
    qt_rime_p: cobj,
    qt_wash_p: cobj,
    qt_evap_p: cobj,
    qttjfl_work_p: cobj,
    hstar_work_p: cobj,
    qm_work_p: cobj,
    pofl_work_p: cobj,
    delz_work_p: cobj,
    rls_work_p: cobj,
    clwc_work_p: cobj,
    ciwc_work_p: cobj,
    cfr_work_p: cobj,
    tem_work_p: cobj,
    evaprate_work_p: cobj,
    qt_rain_work_p: cobj,
    qt_rime_work_p: cobj,
    qt_wash_work_p: cobj,
    qt_evap_work_p: cobj,
    cfxx_work_p: cobj,
    qtt_work_p: cobj,
    qttnew_work_p: cobj,
):
    return _wetchem.neu_wetdep_washo_columns_codon(
        ncol,
        lpar,
        ntrace,
        hno3_ndx,
        do_diag,
        dempirical_impl,
        dtscav,
        adj_factor,
        trc_mass_p,
        qm_p,
        pofl_p,
        delz_p,
        rls_p,
        clwc_p,
        ciwc_p,
        cfr_p,
        tem_p,
        evaprate_p,
        garea_p,
        hstar_p,
        tcmass_p,
        tckaqb_p,
        tcnion_p,
        qt_rain_p,
        qt_rime_p,
        qt_wash_p,
        qt_evap_p,
        qttjfl_work_p,
        hstar_work_p,
        qm_work_p,
        pofl_work_p,
        delz_work_p,
        rls_work_p,
        clwc_work_p,
        ciwc_work_p,
        cfr_work_p,
        tem_work_p,
        evaprate_work_p,
        qt_rain_work_p,
        qt_rime_work_p,
        qt_wash_work_p,
        qt_evap_work_p,
        cfxx_work_p,
        qtt_work_p,
        qttnew_work_p,
    )

@export
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
    return _wetchem.setsox_init_fields_codon(
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

@export
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
    return _wetchem.setsox_ph_solve_codon(
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

@export
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
    return _wetchem.setsox_aqchem_predict_codon(
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

@export
def setsox_xph_lwc_diag_codon(
    ncol: int,
    pcols: int,
    pver: int,
    cldfrc_p: cobj,
    lwc_p: cobj,
    xph_p: cobj,
    xphlwc_p: cobj,
):
    return _wetchem.setsox_xph_lwc_diag_codon(
        ncol,
        pcols,
        pver,
        cldfrc_p,
        lwc_p,
        xph_p,
        xphlwc_p,
    )

@export
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
    return _wetchem.sox_cldaero_update_core_codon(
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
        delso4_hprxn_p,
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

@export
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
    return _wetchem.sox_cldaero_finalize_codon(
        ncol,
        pver,
        gas_pcnst,
        ntot_amode,
        loffset,
        id_so2,
        id_nh3,
        small_value,
        specmw_so4_amode,
        gravit,
        mbar_p,
        pdel_p,
        qcw_p,
        qin_p,
        dqdt_aqso4_p,
        dqdt_aqh2so4_p,
        dqdt_aqhprxn_p,
        dqdt_aqo3rxn_p,
        sflx_aqso4_p,
        sflx_aqh2so4_p,
        sflx_aqhprxn_p,
        sflx_aqo3rxn_p,
        adv_mass_p,
        lptr_so4_cw_amode_p,
        lptr_msa_cw_amode_p,
        lptr_nh4_cw_amode_p,
    )

@export
def setsox_codon(stage: int) -> int:
    return stage


@export
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
    return _wetchem.setsox_shell_codon(
        stage,
        ncol,
        pcols,
        pver,
        gas_pcnst,
        nfs,
        ntot_amode,
        loffset,
        itermax,
        cloud_borne_flag,
        modal_aerosols_flag,
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
        modeptr_accum,
        dtime,
        ph0,
        const0,
        kh0,
        kh1,
        kh2,
        kh3,
        ra,
        xkw,
        so4_fact,
        pi_val,
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
        press_p,
        tfld_p,
        cldfrc_p,
        cldnum_p,
        lwc_p,
        xlwc_p,
        xso4c_p,
        xnh4c_p,
        xno3c_p,
        xso4_init_p,
        xdelso4hp_p,
        hno3g_p,
        nh3g_p,
        hehno3_p,
        heh2o2_p,
        heso2_p,
        henh3_p,
        heo3_p,
        xphlwc_p,
        qcw_p,
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


@export
def setsox_shell_stage_dispatch_codon(
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
    return setsox_shell_codon(
        stage,
        ncol,
        pcols,
        pver,
        gas_pcnst,
        nfs,
        ntot_amode,
        loffset,
        itermax,
        cloud_borne_flag,
        modal_aerosols_flag,
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
        modeptr_accum,
        dtime,
        ph0,
        const0,
        kh0,
        kh1,
        kh2,
        kh3,
        ra,
        xkw,
        so4_fact,
        pi_val,
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
        press_p,
        tfld_p,
        cldfrc_p,
        cldnum_p,
        lwc_p,
        xlwc_p,
        xso4c_p,
        xnh4c_p,
        xno3c_p,
        xso4_init_p,
        xdelso4hp_p,
        hno3g_p,
        nh3g_p,
        hehno3_p,
        heh2o2_p,
        heso2_p,
        henh3_p,
        heo3_p,
        xphlwc_p,
        qcw_p,
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


@export
def setsox_shell_finalize_codon(
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
    small_value: float,
    specmw_so4_amode: float,
    gravit: float,
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
    mbar_p: cobj,
    pdel_p: cobj,
    sflx_aqso4_p: cobj,
    sflx_aqh2so4_p: cobj,
    sflx_aqhprxn_p: cobj,
    sflx_aqo3rxn_p: cobj,
    adv_mass_p: cobj,
):
    return _wetchem.setsox_shell_finalize_codon(
        ncol,
        pcols,
        pver,
        gas_pcnst,
        nfs,
        ntot_amode,
        loffset,
        itermax,
        cloud_borne_flag,
        modal_aerosols_flag,
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
        modeptr_accum,
        dtime,
        ph0,
        const0,
        kh0,
        kh1,
        kh2,
        kh3,
        ra,
        xkw,
        so4_fact,
        pi_val,
        small_value,
        specmw_so4_amode,
        gravit,
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
        press_p,
        tfld_p,
        cldfrc_p,
        cldnum_p,
        lwc_p,
        xlwc_p,
        xso4c_p,
        xnh4c_p,
        xno3c_p,
        xso4_init_p,
        xdelso4hp_p,
        hno3g_p,
        nh3g_p,
        hehno3_p,
        heh2o2_p,
        heso2_p,
        henh3_p,
        heo3_p,
        xphlwc_p,
        qcw_p,
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
        mbar_p,
        pdel_p,
        sflx_aqso4_p,
        sflx_aqh2so4_p,
        sflx_aqhprxn_p,
        sflx_aqo3rxn_p,
        adv_mass_p,
    )

@export
def aero_model_emissions_codon(stage: int) -> int:
    return stage

@export
def aero_model_emissions_accumulate_sflx_codon(
    ncol: int,
    pcols: int,
    nindices: int,
    indices_p: cobj,
    cflx_p: cobj,
    sflx_p: cobj,
):
    return _emissions.aero_model_emissions_accumulate_sflx_codon(
        ncol,
        pcols,
        nindices,
        indices_p,
        cflx_p,
        sflx_p,
    )

@export
def dust_emis_codon(
    ncol: int,
    pcols: int,
    ndstflx: int,
    dust_nbin: int,
    soil_erod_fact: float,
    pi_val: float,
    dust_density: float,
    soil_erod_threshold: float,
    dust_flux_in_p: cobj,
    cflx_p: cobj,
    soil_erod_p: cobj,
    soil_erodibility_p: cobj,
    dust_indices_p: cobj,
    dust_emis_sclfctr_p: cobj,
    dust_dmt_vwr_p: cobj,
):
    return _emissions.dust_emis_codon(
        ncol,
        pcols,
        ndstflx,
        dust_nbin,
        soil_erod_fact,
        pi_val,
        dust_density,
        soil_erod_threshold,
        dust_flux_in_p,
        cflx_p,
        soil_erod_p,
        soil_erodibility_p,
        dust_indices_p,
        dust_emis_sclfctr_p,
        dust_dmt_vwr_p,
    )

@export
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
    return _emissions.aero_model_emissions_seasalt_wind_codon(
        ncol,
        pcols,
        pver,
        z0,
        state_u_p,
        state_v_p,
        state_zm_p,
        u10cubed_p,
    )

@export
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
    dust_density: float,
    seasalt_density: float,
    dstflx_p: cobj,
    soil_erod_p: cobj,
    dust_flux_sum_p: cobj,
    fi_p: cobj,
    ocnfrac_p: cobj,
    cflx_p: cobj,
    dust_sflx_p: cobj,
    seasalt_sflx_p: cobj,
    dust_indices_p: cobj,
    dust_emis_sclfctr_p: cobj,
    dust_dmt_vwr_p: cobj,
    seasalt_indices_p: cobj,
    seasalt_sz_range_lo_p: cobj,
    seasalt_sz_range_hi_p: cobj,
    dg_p: cobj,
    rdry_p: cobj,
    soil_erodibility_p: cobj,
    soil_erod_threshold: float,
    sst_p: cobj,
    u10cubed_p: cobj,
    whitecap_p: cobj,
    consta_p: cobj,
    constb_p: cobj,
):
    return _emissions.aero_model_emissions_shell_codon(
        stage,
        ncol,
        pcols,
        ndstflx,
        dust_nbin,
        seasalt_nbin,
        nsections,
        soil_erod_fact,
        seasalt_emis_scale,
        pi_val,
        dust_density,
        seasalt_density,
        dstflx_p,
        soil_erod_p,
        dust_flux_sum_p,
        fi_p,
        ocnfrac_p,
        cflx_p,
        dust_sflx_p,
        seasalt_sflx_p,
        dust_indices_p,
        dust_emis_sclfctr_p,
        dust_dmt_vwr_p,
        seasalt_indices_p,
        seasalt_sz_range_lo_p,
        seasalt_sz_range_hi_p,
        dg_p,
        rdry_p,
        soil_erodibility_p,
        soil_erod_threshold,
        sst_p,
        u10cubed_p,
        whitecap_p,
        consta_p,
        constb_p,
    )


@export
def aero_model_emissions_shell_wind_codon(
    stage: int,
    ncol: int,
    pcols: int,
    pver: int,
    ndstflx: int,
    dust_nbin: int,
    seasalt_nbin: int,
    nsections: int,
    soil_erod_fact: float,
    seasalt_emis_scale: float,
    pi_val: float,
    dust_density: float,
    seasalt_density: float,
    dstflx_p: cobj,
    soil_erod_p: cobj,
    dust_flux_sum_p: cobj,
    fi_p: cobj,
    ocnfrac_p: cobj,
    cflx_p: cobj,
    dust_sflx_p: cobj,
    seasalt_sflx_p: cobj,
    dust_indices_p: cobj,
    dust_emis_sclfctr_p: cobj,
    dust_dmt_vwr_p: cobj,
    seasalt_indices_p: cobj,
    seasalt_sz_range_lo_p: cobj,
    seasalt_sz_range_hi_p: cobj,
    dg_p: cobj,
    rdry_p: cobj,
    soil_erodibility_p: cobj,
    soil_erod_threshold: float,
    sst_p: cobj,
    u10cubed_p: cobj,
    whitecap_p: cobj,
    consta_p: cobj,
    constb_p: cobj,
    compute_wind: int,
    z0: float,
    state_u_p: cobj,
    state_v_p: cobj,
    state_zm_p: cobj,
):
    return _emissions.aero_model_emissions_shell_wind_codon(
        stage,
        ncol,
        pcols,
        pver,
        ndstflx,
        dust_nbin,
        seasalt_nbin,
        nsections,
        soil_erod_fact,
        seasalt_emis_scale,
        pi_val,
        dust_density,
        seasalt_density,
        dstflx_p,
        soil_erod_p,
        dust_flux_sum_p,
        fi_p,
        ocnfrac_p,
        cflx_p,
        dust_sflx_p,
        seasalt_sflx_p,
        dust_indices_p,
        dust_emis_sclfctr_p,
        dust_dmt_vwr_p,
        seasalt_indices_p,
        seasalt_sz_range_lo_p,
        seasalt_sz_range_hi_p,
        dg_p,
        rdry_p,
        soil_erodibility_p,
        soil_erod_threshold,
        sst_p,
        u10cubed_p,
        whitecap_p,
        consta_p,
        constb_p,
        compute_wind,
        z0,
        state_u_p,
        state_v_p,
        state_zm_p,
    )

@export
def aero_model_emissions_shell_wind_stage_dispatch_codon(
    stage: int,
    ncol: int,
    pcols: int,
    pver: int,
    ndstflx: int,
    dust_nbin: int,
    seasalt_nbin: int,
    nsections: int,
    soil_erod_fact: float,
    seasalt_emis_scale: float,
    pi_val: float,
    dust_density: float,
    seasalt_density: float,
    dstflx_p: cobj,
    soil_erod_p: cobj,
    dust_flux_sum_p: cobj,
    fi_p: cobj,
    ocnfrac_p: cobj,
    cflx_p: cobj,
    dust_sflx_p: cobj,
    seasalt_sflx_p: cobj,
    dust_indices_p: cobj,
    dust_emis_sclfctr_p: cobj,
    dust_dmt_vwr_p: cobj,
    seasalt_indices_p: cobj,
    seasalt_sz_range_lo_p: cobj,
    seasalt_sz_range_hi_p: cobj,
    dg_p: cobj,
    rdry_p: cobj,
    soil_erodibility_p: cobj,
    soil_erod_threshold: float,
    sst_p: cobj,
    u10cubed_p: cobj,
    whitecap_p: cobj,
    consta_p: cobj,
    constb_p: cobj,
    compute_wind: int,
    z0: float,
    state_u_p: cobj,
    state_v_p: cobj,
    state_zm_p: cobj,
):
    return aero_model_emissions_shell_wind_codon(
        stage,
        ncol,
        pcols,
        pver,
        ndstflx,
        dust_nbin,
        seasalt_nbin,
        nsections,
        soil_erod_fact,
        seasalt_emis_scale,
        pi_val,
        dust_density,
        seasalt_density,
        dstflx_p,
        soil_erod_p,
        dust_flux_sum_p,
        fi_p,
        ocnfrac_p,
        cflx_p,
        dust_sflx_p,
        seasalt_sflx_p,
        dust_indices_p,
        dust_emis_sclfctr_p,
        dust_dmt_vwr_p,
        seasalt_indices_p,
        seasalt_sz_range_lo_p,
        seasalt_sz_range_hi_p,
        dg_p,
        rdry_p,
        soil_erodibility_p,
        soil_erod_threshold,
        sst_p,
        u10cubed_p,
        whitecap_p,
        consta_p,
        constb_p,
        compute_wind,
        z0,
        state_u_p,
        state_v_p,
        state_zm_p,
    )

@export
def chem_timestep_init_should_run_codon(
    nstep: int,
    chem_freq: int,
    chem_step_flag_p: cobj,
):
    return _gas_phase.chem_timestep_init_should_run_codon(
        nstep,
        chem_freq,
        chem_step_flag_p,
    )

@export
def chem_timestep_init_codon(
    nstep: int,
    chem_freq: int,
    chem_step_flag_p: cobj,
):
    return chem_timestep_init_should_run_codon(
        nstep,
        chem_freq,
        chem_step_flag_p,
    )

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
    return _gas_phase.chem_timestep_tend_fill_cloud_fields_codon(
        ncol,
        pcols,
        pver,
        pcnst,
        ixcldliq,
        ixcldice,
        ixndrop,
        state_q_p,
        cldw_p,
        ncldwtr_p,
    )

@export
def chem_timestep_tend_init_lq_codon(
    pcnst: int,
    ghg_chem: int,
    map2chm_p: cobj,
    lq_mask_p: cobj,
):
    return _gas_phase.chem_timestep_tend_init_lq_codon(
        pcnst,
        ghg_chem,
        map2chm_p,
        lq_mask_p,
    )

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
    return _gas_phase.chem_timestep_tend_apply_depflux_codon(
        ncol,
        pcols,
        idx_cb1,
        idx_cb2,
        idx_oc1,
        idx_oc2,
        drydepflx_p,
        bcphodry_p,
        bcphidry_p,
        ocphodry_p,
        ocphidry_p,
    )

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
    return _gas_phase.chem_timestep_tend_sum_fh2o_codon(
        ncol,
        pcols,
        pver,
        gravit,
        ptend_q1_p,
        pdel_p,
        fh2o_p,
    )

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
    return _gas_phase.gas_phase_chemdr_finalize_tendencies_codon(
        ncol,
        pcols,
        pver,
        gas_pcnst,
        pcnst,
        delt_inverse,
        map2chm_p,
        mmr_p,
        mmr_tend_p,
        mmr_new_p,
        qtend_p,
    )

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
    return _gas_phase.gas_phase_chemdr_prepare_state_codon(
        ncol,
        pcols,
        pver,
        rga,
        m2km,
        pa2mb,
        phis_p,
        zi_p,
        zm_p,
        pmid_p,
        zsurf_p,
        zintr_p,
        zmidr_p,
        zmid_p,
        zint_p,
        pmb_p,
    )

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
    return _gas_phase.gas_phase_chemdr_load_mmr_codon(
        ncol,
        pcols,
        pver,
        pcnst,
        map2chm_p,
        q_p,
        mmr_p,
    )

@export
def rate_diags_batch_codon(
    ncol: int,
    pver: int,
    rxntot: int,
    rxt_tag_cnt: int,
    rxt_rates_p: cobj,
    vmr_p: cobj,
    m_p: cobj,
    rxt_tag_map_p: cobj,
):
    return _gas_phase.rate_diags_batch_codon(
        ncol,
        pver,
        rxntot,
        rxt_tag_cnt,
        rxt_rates_p,
        vmr_p,
        m_p,
        rxt_tag_map_p,
    )

@export
def rate_diags_batch_stage_dispatch_codon(
    ncol: int,
    pver: int,
    rxntot: int,
    rxt_tag_cnt: int,
    rxt_rates_p: cobj,
    vmr_p: cobj,
    m_p: cobj,
    rxt_tag_map_p: cobj,
):
    return rate_diags_batch_codon(
        ncol,
        pver,
        rxntot,
        rxt_tag_cnt,
        rxt_rates_p,
        vmr_p,
        m_p,
        rxt_tag_map_p,
    )

@export
def rate_diags_calc_codon(
    ncol: int,
    pver: int,
    rxntot: int,
    rxt_tag_cnt: int,
    rxt_rates_p: cobj,
    vmr_p: cobj,
    m_p: cobj,
    rxt_tag_map_p: cobj,
    ngrps: int,
    max_group_members: int,
    grp_nm_p: cobj,
    grp_map_p: cobj,
    grp_mult_p: cobj,
    group_rates_p: cobj,
):
    return _gas_phase.rate_diags_calc_codon(
        ncol,
        pver,
        rxntot,
        rxt_tag_cnt,
        rxt_rates_p,
        vmr_p,
        m_p,
        rxt_tag_map_p,
        ngrps,
        max_group_members,
        grp_nm_p,
        grp_map_p,
        grp_mult_p,
        group_rates_p,
    )

@export
def gas_phase_chemdr_init_reaction_rates_codon(
    ncol: int,
    pver: int,
    rxntot: int,
    nan_value_p: cobj,
    reaction_rates_p: cobj,
):
    return _gas_phase.gas_phase_chemdr_init_reaction_rates_codon(
        ncol,
        pver,
        rxntot,
        nan_value_p,
        reaction_rates_p,
    )

@export
def gas_phase_chemdr_clip_sulfate_codon(
    ncol: int,
    pcols: int,
    pver: int,
    troplev_p: cobj,
    sulfate_p: cobj,
):
    return _gas_phase.gas_phase_chemdr_clip_sulfate_codon(
        ncol,
        pcols,
        pver,
        troplev_p,
        sulfate_p,
    )

@export
def gas_phase_chemdr_zero_het_rates_codon(
    ncol: int,
    pver: int,
    gas_pcnst_dim: int,
    het_rates_p: cobj,
):
    return _gas_phase.gas_phase_chemdr_zero_het_rates_codon(
        ncol,
        pver,
        gas_pcnst_dim,
        het_rates_p,
    )

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
    return _gas_phase.gas_phase_chemdr_load_oxygen_mmr_codon(
        ncol,
        pcols,
        pver,
        o2_ndx,
        o_ndx,
        mmr_p,
        o2mmr_p,
        ommr_p,
    )

@export
def gas_phase_chemdr_set_ltrop_sol_codon(
    ncol: int,
    has_linoz_data_flag: int,
    troplev_p: cobj,
    ltrop_sol_p: cobj,
):
    return _gas_phase.gas_phase_chemdr_set_ltrop_sol_codon(
        ncol,
        has_linoz_data_flag,
        troplev_p,
        ltrop_sol_p,
    )

@export
def gas_phase_chemdr_zero_st80_tau_codon(
    ncol: int,
    pver: int,
    rxntot: int,
    st80_25_tau_ndx: int,
    troplev_p: cobj,
    reaction_rates_p: cobj,
):
    return _gas_phase.gas_phase_chemdr_zero_st80_tau_codon(
        ncol,
        pver,
        rxntot,
        st80_25_tau_ndx,
        troplev_p,
        reaction_rates_p,
    )

@export
def gas_phase_chemdr_compute_relhum_codon(
    ncol: int,
    pver: int,
    h2ovmr_p: cobj,
    satq_p: cobj,
    relhum_p: cobj,
):
    return _gas_phase.gas_phase_chemdr_compute_relhum_codon(
        ncol,
        pver,
        h2ovmr_p,
        satq_p,
        relhum_p,
    )

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
    return _gas_phase.gas_phase_chemdr_restore_strat_gases_codon(
        ncol,
        pver,
        gas_pcnst,
        hno3_ndx,
        h2o_ndx,
        delt_inverse,
        vmr_p,
        hno3_gas_p,
        h2o_gas_p,
        h2ovmr_p,
        wrk_p,
    )

@export
def gas_phase_chemdr_restore_hcl_gas_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    hcl_ndx: int,
    vmr_p: cobj,
    hcl_gas_p: cobj,
):
    return _gas_phase.gas_phase_chemdr_restore_hcl_gas_codon(
        ncol,
        pver,
        gas_pcnst,
        hcl_ndx,
        vmr_p,
        hcl_gas_p,
    )

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
    return _gas_phase.gas_phase_chemdr_init_dust_vmr_codon(
        ncol,
        pver,
        gas_pcnst,
        ndust,
        dst_ndx,
        vmr_p,
        dust_vmr_p,
    )

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
    return _gas_phase.gas_phase_chemdr_reset_ste_tracer_codon(
        ncol,
        pcols,
        pver,
        gas_pcnst,
        st80_25_ndx,
        pmid_threshold,
        st80_vmr,
        pmid_p,
        vmr_p,
    )

@export
def gas_phase_chemdr_zero_sflx_codon(
    pcols: int,
    gas_pcnst: int,
    sflx_p: cobj,
):
    return _gas_phase.gas_phase_chemdr_zero_sflx_codon(
        pcols,
        gas_pcnst,
        sflx_p,
    )

@export
def gas_phase_chemdr_compute_wind_speed_codon(
    ncol: int,
    pcols: int,
    pver: int,
    ufld_p: cobj,
    vfld_p: cobj,
    wind_speed_p: cobj,
):
    return _gas_phase.gas_phase_chemdr_compute_wind_speed_codon(
        ncol,
        pcols,
        pver,
        ufld_p,
        vfld_p,
        wind_speed_p,
    )

@export
def gas_phase_chemdr_compute_prect_codon(
    ncol: int,
    pcols: int,
    precc_p: cobj,
    precl_p: cobj,
    prect_p: cobj,
):
    return _gas_phase.gas_phase_chemdr_compute_prect_codon(
        ncol,
        pcols,
        precc_p,
        precl_p,
        prect_p,
    )

@export
def gas_phase_chemdr_compute_tvs_codon(
    ncol: int,
    pcols: int,
    pver: int,
    tfld_p: cobj,
    qh2o_p: cobj,
    tvs_p: cobj,
):
    return _gas_phase.gas_phase_chemdr_compute_tvs_codon(
        ncol,
        pcols,
        pver,
        tfld_p,
        qh2o_p,
        tvs_p,
    )

@export
def gas_phase_chemdr_surface_diag_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    pm25_flag: int,
    pm25_soa_flag: int,
    cb1_ndx: int,
    cb2_ndx: int,
    oc1_ndx: int,
    oc2_ndx: int,
    dst1_ndx: int,
    dst2_ndx: int,
    sslt1_ndx: int,
    sslt2_ndx: int,
    soa_ndx: int,
    soam_ndx: int,
    soai_ndx: int,
    soat_ndx: int,
    soab_ndx: int,
    soax_ndx: int,
    so4_ndx: int,
    mmr_new_p: cobj,
    qh2o_p: cobj,
    ufld_p: cobj,
    vfld_p: cobj,
    pm25_p: cobj,
    pm25_soa_p: cobj,
    q_srf_p: cobj,
    u_srf_p: cobj,
    v_srf_p: cobj,
):
    return _gas_phase.gas_phase_chemdr_surface_diag_codon(
        ncol,
        pcols,
        pver,
        gas_pcnst,
        pm25_flag,
        pm25_soa_flag,
        cb1_ndx,
        cb2_ndx,
        oc1_ndx,
        oc2_ndx,
        dst1_ndx,
        dst2_ndx,
        sslt1_ndx,
        sslt2_ndx,
        soa_ndx,
        soam_ndx,
        soai_ndx,
        soat_ndx,
        soab_ndx,
        soax_ndx,
        so4_ndx,
        mmr_new_p,
        qh2o_p,
        ufld_p,
        vfld_p,
        pm25_p,
        pm25_soa_p,
        q_srf_p,
        u_srf_p,
        v_srf_p,
    )

@export
def gas_phase_chemdr_copy_cldw_to_cwat_codon(
    ncol: int,
    pcols: int,
    pver: int,
    cldw_p: cobj,
    cwat_p: cobj,
):
    return _gas_phase.gas_phase_chemdr_copy_cldw_to_cwat_codon(
        ncol,
        pcols,
        pver,
        cldw_p,
        cwat_p,
    )

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
    return _gas_phase.gas_phase_chemdr_load_h2o_fields_codon(
        ncol,
        pcols,
        pver,
        gas_pcnst,
        h2o_ndx,
        mmr_p,
        vmr_p,
        qh2o_p,
        h2ovmr_p,
    )

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
    return _gas_phase.gas_phase_chemdr_copy_o3_to_o3s_trop_codon(
        ncol,
        pcols,
        pver,
        gas_pcnst,
        troplev_p,
        o3_ndx,
        o3s_ndx,
        vmr_p,
    )

@export
def gas_phase_chemdr_copy_h2o_to_wrk_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    h2o_ndx: int,
    vmr_p: cobj,
    wrk_p: cobj,
):
    return _gas_phase.gas_phase_chemdr_copy_h2o_to_wrk_codon(
        ncol,
        pver,
        gas_pcnst,
        h2o_ndx,
        vmr_p,
        wrk_p,
    )

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
    return _gas_phase.gas_phase_chemdr_update_qdsett_wrk_codon(
        ncol,
        pver,
        gas_pcnst,
        h2o_ndx,
        delt_inverse,
        vmr_p,
        wrk_p,
    )

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
    return _gas_phase.gas_phase_chemdr_update_qdchem_wrk_codon(
        ncol,
        pver,
        gas_pcnst,
        h2o_ndx,
        delt_inverse,
        vmr_p,
        wrk_p,
    )

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
    return _gas_phase.gas_phase_chemdr_init_stratchem_state_codon(
        ncol,
        pcols,
        pver,
        gas_pcnst,
        pcnst,
        hno3_ndx,
        hcl_ndx,
        cldice_ndx,
        vmr_p,
        h2ovmr_p,
        q_p,
        hcl_cond_p,
        hcl_gas_p,
        hno3_gas_p,
        h2o_gas_p,
        wrk_p,
        cldice_p,
        hno3_cond_p,
    )

@export
def gas_phase_chemdr_init_h2so4_gasprod_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    ndx_h2so4: int,
    vmr_p: cobj,
    del_h2so4_gasprod_p: cobj,
):
    return _gas_phase.gas_phase_chemdr_init_h2so4_gasprod_codon(
        ncol,
        pver,
        gas_pcnst,
        ndx_h2so4,
        vmr_p,
        del_h2so4_gasprod_p,
    )

@export
def gas_phase_chemdr_store_vmr0_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    vmr_p: cobj,
    vmr0_p: cobj,
):
    return _gas_phase.gas_phase_chemdr_store_vmr0_codon(
        ncol,
        pver,
        gas_pcnst,
        vmr_p,
        vmr0_p,
    )

@export
def gas_phase_chemdr_update_h2so4_gasprod_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    ndx_h2so4: int,
    vmr_p: cobj,
    del_h2so4_gasprod_p: cobj,
):
    return _gas_phase.gas_phase_chemdr_update_h2so4_gasprod_codon(
        ncol,
        pver,
        gas_pcnst,
        ndx_h2so4,
        vmr_p,
        del_h2so4_gasprod_p,
    )

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
    return _gas_phase.gas_phase_chemdr_reform_hno3_hcl_codon(
        ncol,
        pver,
        gas_pcnst,
        hno3_ndx,
        hcl_ndx,
        vmr_p,
        hno3_cond_p,
        hcl_cond_p,
    )

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
    return _gas_phase.gas_phase_chemdr_normalize_extfrc_codon(
        ncol,
        pver,
        extcnt,
        synoz_ndx,
        aoa_nh_ext_ndx,
        indexm,
        extfrc_p,
        invariants_p,
    )

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
    return _gas_phase.gas_phase_chemdr_store_drydep_codon(
        ncol,
        pcols,
        gas_pcnst,
        pcnst,
        map2chm_p,
        sflx_p,
        cflx_p,
        drydepflx_p,
    )

@export
def gas_phase_chemdr_stage_dispatch_codon(
    stage: int,
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    pcnst: int,
    rxntot: int,
    extcnt: int,
    nfs: int,
    indexm: int,
    ncol_abs: int,
    jo1d_adj_ndx: int,
    inv_n2_ndx: int,
    inv_o2_ndx: int,
    inv_h2o_ndx: int,
    has_linoz_data_flag: int,
    h2o_ndx: int,
    o2_ndx: int,
    o_ndx: int,
    hno3_ndx: int,
    hcl_ndx: int,
    cldice_ndx: int,
    st80_25_ndx: int,
    aoa_nh_ndx: int,
    nh_5_ndx: int,
    nh_50_ndx: int,
    nh_50w_ndx: int,
    so4_ndx: int,
    st80_25_tau_ndx: int,
    ndx_h2so4: int,
    o3_ndx: int,
    o3s_ndx: int,
    synoz_ndx: int,
    aoa_nh_ext_ndx: int,
    h_ndx: int,
    n_ndx: int,
    elec_ndx: int,
    np_ndx: int,
    n2p_ndx: int,
    op_ndx: int,
    o2p_ndx: int,
    nop_ndx: int,
    fixed_mbar: int,
    rad2deg: float,
    delt: float,
    delt_inverse: float,
    rga: float,
    m2km: float,
    pa2mb: float,
    mwdry: float,
    adv_mass_h2o: float,
    map2chm_p: cobj,
    adv_mass_p: cobj,
    troplev_p: cobj,
    ltrop_sol_p: cobj,
    zen_angle_p: cobj,
    sza_p: cobj,
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
    q_p: cobj,
    mmr_p: cobj,
    mbar_p: cobj,
    vmr_p: cobj,
    qh2o_p: cobj,
    h2ovmr_p: cobj,
    rlats_p: cobj,
    sulfate_p: cobj,
    satq_p: cobj,
    relhum_p: cobj,
    cldw_p: cobj,
    cwat_p: cobj,
    extfrc_p: cobj,
    invariants_p: cobj,
    het_rates_p: cobj,
    reaction_rates_p: cobj,
    col_delta_p: cobj,
    col_dens_p: cobj,
    del_h2so4_gasprod_p: cobj,
    vmr0_p: cobj,
    o3s_loss_p: cobj,
    mmr_tend_p: cobj,
    mmr_new_p: cobj,
    qtend_p: cobj,
    tfld_p: cobj,
    tvs_p: cobj,
    sflx_p: cobj,
    ufld_p: cobj,
    vfld_p: cobj,
    wind_speed_p: cobj,
    precc_p: cobj,
    precl_p: cobj,
    prect_p: cobj,
    cflx_p: cobj,
    drydepflx_p: cobj,
    o2mmr_p: cobj,
    ommr_p: cobj,
    hcl_cond_p: cobj,
    hcl_gas_p: cobj,
    hno3_gas_p: cobj,
    h2o_gas_p: cobj,
    wrk_p: cobj,
    cldice_p: cobj,
    hno3_cond_p: cobj,
):
    return _gas_phase.gas_phase_chemdr_shell_codon(
        stage,
        ncol,
        pcols,
        pver,
        gas_pcnst,
        pcnst,
        rxntot,
        extcnt,
        nfs,
        indexm,
        ncol_abs,
        jo1d_adj_ndx,
        inv_n2_ndx,
        inv_o2_ndx,
        inv_h2o_ndx,
        has_linoz_data_flag,
        h2o_ndx,
        o2_ndx,
        o_ndx,
        hno3_ndx,
        hcl_ndx,
        cldice_ndx,
        st80_25_ndx,
        aoa_nh_ndx,
        nh_5_ndx,
        nh_50_ndx,
        nh_50w_ndx,
        so4_ndx,
        st80_25_tau_ndx,
        ndx_h2so4,
        o3_ndx,
        o3s_ndx,
        synoz_ndx,
        aoa_nh_ext_ndx,
        h_ndx,
        n_ndx,
        elec_ndx,
        np_ndx,
        n2p_ndx,
        op_ndx,
        o2p_ndx,
        nop_ndx,
        fixed_mbar,
        rad2deg,
        delt,
        delt_inverse,
        rga,
        m2km,
        pa2mb,
        mwdry,
        adv_mass_h2o,
        map2chm_p,
        adv_mass_p,
        troplev_p,
        ltrop_sol_p,
        zen_angle_p,
        sza_p,
        phis_p,
        zi_p,
        zm_p,
        pmid_p,
        zsurf_p,
        zintr_p,
        zmidr_p,
        zmid_p,
        zint_p,
        pmb_p,
        q_p,
        mmr_p,
        mbar_p,
        vmr_p,
        qh2o_p,
        h2ovmr_p,
        rlats_p,
        sulfate_p,
        satq_p,
        relhum_p,
        cldw_p,
        cwat_p,
        extfrc_p,
        invariants_p,
        het_rates_p,
        reaction_rates_p,
        col_delta_p,
        col_dens_p,
        del_h2so4_gasprod_p,
        vmr0_p,
        o3s_loss_p,
        mmr_tend_p,
        mmr_new_p,
        qtend_p,
        tfld_p,
        tvs_p,
        sflx_p,
        ufld_p,
        vfld_p,
        wind_speed_p,
        precc_p,
        precl_p,
        prect_p,
        cflx_p,
        drydepflx_p,
        o2mmr_p,
        ommr_p,
        hcl_cond_p,
        hcl_gas_p,
        hno3_gas_p,
        h2o_gas_p,
        wrk_p,
        cldice_p,
        hno3_cond_p,
    )

@export
def gas_phase_chemdr_shell_codon(
    stage: int,
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    pcnst: int,
    rxntot: int,
    extcnt: int,
    nfs: int,
    indexm: int,
    ncol_abs: int,
    jo1d_adj_ndx: int,
    inv_n2_ndx: int,
    inv_o2_ndx: int,
    inv_h2o_ndx: int,
    has_linoz_data_flag: int,
    h2o_ndx: int,
    o2_ndx: int,
    o_ndx: int,
    hno3_ndx: int,
    hcl_ndx: int,
    cldice_ndx: int,
    st80_25_ndx: int,
    aoa_nh_ndx: int,
    nh_5_ndx: int,
    nh_50_ndx: int,
    nh_50w_ndx: int,
    so4_ndx: int,
    st80_25_tau_ndx: int,
    ndx_h2so4: int,
    o3_ndx: int,
    o3s_ndx: int,
    synoz_ndx: int,
    aoa_nh_ext_ndx: int,
    h_ndx: int,
    n_ndx: int,
    elec_ndx: int,
    np_ndx: int,
    n2p_ndx: int,
    op_ndx: int,
    o2p_ndx: int,
    nop_ndx: int,
    fixed_mbar: int,
    rad2deg: float,
    delt: float,
    delt_inverse: float,
    rga: float,
    m2km: float,
    pa2mb: float,
    mwdry: float,
    adv_mass_h2o: float,
    map2chm_p: cobj,
    adv_mass_p: cobj,
    troplev_p: cobj,
    ltrop_sol_p: cobj,
    zen_angle_p: cobj,
    sza_p: cobj,
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
    q_p: cobj,
    mmr_p: cobj,
    mbar_p: cobj,
    vmr_p: cobj,
    qh2o_p: cobj,
    h2ovmr_p: cobj,
    rlats_p: cobj,
    sulfate_p: cobj,
    satq_p: cobj,
    relhum_p: cobj,
    cldw_p: cobj,
    cwat_p: cobj,
    extfrc_p: cobj,
    invariants_p: cobj,
    het_rates_p: cobj,
    reaction_rates_p: cobj,
    col_delta_p: cobj,
    col_dens_p: cobj,
    del_h2so4_gasprod_p: cobj,
    vmr0_p: cobj,
    o3s_loss_p: cobj,
    mmr_tend_p: cobj,
    mmr_new_p: cobj,
    qtend_p: cobj,
    tfld_p: cobj,
    tvs_p: cobj,
    sflx_p: cobj,
    ufld_p: cobj,
    vfld_p: cobj,
    wind_speed_p: cobj,
    precc_p: cobj,
    precl_p: cobj,
    prect_p: cobj,
    cflx_p: cobj,
    drydepflx_p: cobj,
    o2mmr_p: cobj,
    ommr_p: cobj,
    hcl_cond_p: cobj,
    hcl_gas_p: cobj,
    hno3_gas_p: cobj,
    h2o_gas_p: cobj,
    wrk_p: cobj,
    cldice_p: cobj,
    hno3_cond_p: cobj,
):
    return gas_phase_chemdr_stage_dispatch_codon(
        stage,
        ncol,
        pcols,
        pver,
        gas_pcnst,
        pcnst,
        rxntot,
        extcnt,
        nfs,
        indexm,
        ncol_abs,
        jo1d_adj_ndx,
        inv_n2_ndx,
        inv_o2_ndx,
        inv_h2o_ndx,
        has_linoz_data_flag,
        h2o_ndx,
        o2_ndx,
        o_ndx,
        hno3_ndx,
        hcl_ndx,
        cldice_ndx,
        st80_25_ndx,
        aoa_nh_ndx,
        nh_5_ndx,
        nh_50_ndx,
        nh_50w_ndx,
        so4_ndx,
        st80_25_tau_ndx,
        ndx_h2so4,
        o3_ndx,
        o3s_ndx,
        synoz_ndx,
        aoa_nh_ext_ndx,
        h_ndx,
        n_ndx,
        elec_ndx,
        np_ndx,
        n2p_ndx,
        op_ndx,
        o2p_ndx,
        nop_ndx,
        fixed_mbar,
        rad2deg,
        delt,
        delt_inverse,
        rga,
        m2km,
        pa2mb,
        mwdry,
        adv_mass_h2o,
        map2chm_p,
        adv_mass_p,
        troplev_p,
        ltrop_sol_p,
        zen_angle_p,
        sza_p,
        phis_p,
        zi_p,
        zm_p,
        pmid_p,
        zsurf_p,
        zintr_p,
        zmidr_p,
        zmid_p,
        zint_p,
        pmb_p,
        q_p,
        mmr_p,
        mbar_p,
        vmr_p,
        qh2o_p,
        h2ovmr_p,
        rlats_p,
        sulfate_p,
        satq_p,
        relhum_p,
        cldw_p,
        cwat_p,
        extfrc_p,
        invariants_p,
        het_rates_p,
        reaction_rates_p,
        col_delta_p,
        col_dens_p,
        del_h2so4_gasprod_p,
        vmr0_p,
        o3s_loss_p,
        mmr_tend_p,
        mmr_new_p,
        qtend_p,
        tfld_p,
        tvs_p,
        sflx_p,
        ufld_p,
        vfld_p,
        wind_speed_p,
        precc_p,
        precl_p,
        prect_p,
        cflx_p,
        drydepflx_p,
        o2mmr_p,
        ommr_p,
        hcl_cond_p,
        hcl_gas_p,
        hno3_gas_p,
        h2o_gas_p,
        wrk_p,
        cldice_p,
        hno3_cond_p,
    )


@export
def set_xnox_photo_codon(
    ncol: int,
    pver: int,
    photos_p: cobj,
    jno2a_ndx: int,
    jno2_ndx: int,
    jn2o5a_ndx: int,
    jn2o5_ndx: int,
    jn2o5b_ndx: int,
    jhno3a_ndx: int,
    jhno3_ndx: int,
    jno3a_ndx: int,
    jno3_ndx: int,
    jho2no2a_ndx: int,
    jho2no2_ndx: int,
    jmpana_ndx: int,
    jmpan_ndx: int,
    jpana_ndx: int,
    jpan_ndx: int,
    jonitra_ndx: int,
    jonitr_ndx: int,
    jo1da_ndx: int,
    jo1d_ndx: int,
    jo3pa_ndx: int,
    jo3p_ndx: int,
):
    return _gas_phase.set_xnox_photo_codon(
        ncol,
        pver,
        photos_p,
        jno2a_ndx,
        jno2_ndx,
        jn2o5a_ndx,
        jn2o5_ndx,
        jn2o5b_ndx,
        jhno3a_ndx,
        jhno3_ndx,
        jno3a_ndx,
        jno3_ndx,
        jho2no2a_ndx,
        jho2no2_ndx,
        jmpana_ndx,
        jmpan_ndx,
        jpana_ndx,
        jpan_ndx,
        jonitra_ndx,
        jonitr_ndx,
        jo1da_ndx,
        jo1d_ndx,
        jo3pa_ndx,
        jo3p_ndx,
    )

@export
def adjrxt_codon(
    ncol: int,
    pver: int,
    rate_p: cobj,
    inv_p: cobj,
    m_p: cobj,
):
    return _gas_phase.adjrxt_codon(
        ncol,
        pver,
        rate_p,
        inv_p,
        m_p,
    )

@export
def setrxt_codon(
    ncol: int,
    pcols: int,
    pver: int,
    temp_p: cobj,
    rate_p: cobj,
):
    return _gas_phase.setrxt_codon(
        ncol,
        pcols,
        pver,
        temp_p,
        rate_p,
    )

@export
def lu_slv_codon(
    lu_p: cobj,
    b_p: cobj,
):
    return _gas_phase.lu_slv_codon(
        lu_p,
        b_p,
    )

@export
def lu_slv01_codon(
    lu_p: cobj,
    b_p: cobj,
):
    return _gas_phase.lu_slv_codon(
        lu_p,
        b_p,
    )

@export
def lu_fac_codon(
    lu_p: cobj,
):
    return _gas_phase.lu_fac_codon(
        lu_p,
    )

@export
def lu_fac01_codon(
    lu_p: cobj,
):
    return _gas_phase.lu_fac_codon(
        lu_p,
    )

@export
def linmat_codon(
    mat_p: cobj,
    rxt_p: cobj,
    het_rates_p: cobj,
):
    return _gas_phase.linmat_codon(
        mat_p,
        rxt_p,
        het_rates_p,
    )

@export
def linmat01_codon(
    mat_p: cobj,
    rxt_p: cobj,
    het_rates_p: cobj,
):
    return _gas_phase.linmat_codon(
        mat_p,
        rxt_p,
        het_rates_p,
    )

@export
def imp_prod_loss_codon(
    prod_p: cobj,
    loss_p: cobj,
    y_p: cobj,
    rxt_p: cobj,
    het_rates_p: cobj,
):
    return _gas_phase.imp_prod_loss_codon(
        prod_p,
        loss_p,
        y_p,
        rxt_p,
        het_rates_p,
    )

@export
def exp_prod_loss_codon():
    return _gas_phase.exp_prod_loss_codon()

@export
def nlnmat_codon(
    mat_p: cobj,
    lmat_p: cobj,
    dti: float,
):
    return _gas_phase.nlnmat_codon(
        mat_p,
        lmat_p,
        dti,
    )

@export
def nlnmat_finit_codon(
    mat_p: cobj,
    lmat_p: cobj,
    dti: float,
):
    return _gas_phase.nlnmat_codon(
        mat_p,
        lmat_p,
        dti,
    )

@export
def imp_sol_inner_batch_codon(
    mode: int,
    factor_flag: int,
    clscnt4: int,
    dti: float,
    lin_jac_p: cobj,
    sys_jac_p: cobj,
    prod_p: cobj,
    loss_p: cobj,
    lsol_p: cobj,
    lrxt_p: cobj,
    lhet_p: cobj,
    solution_p: cobj,
    iter_invariant_p: cobj,
    forcing_p: cobj,
):
    return _gas_phase.imp_sol_inner_batch_codon(
        mode,
        factor_flag,
        clscnt4,
        dti,
        lin_jac_p,
        sys_jac_p,
        prod_p,
        loss_p,
        lsol_p,
        lrxt_p,
        lhet_p,
        solution_p,
        iter_invariant_p,
        forcing_p,
    )

@export
def imp_sol_outer_batch_codon(
    mode: int,
    i: int,
    lev: int,
    nr_iter: int,
    has_independent: int,
    ncol: int,
    pver: int,
    gas_pcnst: int,
    rxntot: int,
    extcnt: int,
    clscnt4: int,
    dti: float,
    small: float,
    base_sol_p: cobj,
    reaction_rates_p: cobj,
    het_rates_p: cobj,
    extfrc_p: cobj,
    ind_prd_p: cobj,
    clsmap4_p: cobj,
    permute4_p: cobj,
    epsilon_p: cobj,
    max_delta_p: cobj,
    converged_code_p: cobj,
    convergence_code_p: cobj,
    lrxt_p: cobj,
    lhet_p: cobj,
    lsol_p: cobj,
    solution_p: cobj,
    iter_invariant_p: cobj,
    forcing_p: cobj,
):
    return _gas_phase.imp_sol_outer_batch_codon(
        mode,
        i,
        lev,
        nr_iter,
        has_independent,
        ncol,
        pver,
        gas_pcnst,
        rxntot,
        extcnt,
        clscnt4,
        dti,
        small,
        base_sol_p,
        reaction_rates_p,
        het_rates_p,
        extfrc_p,
        ind_prd_p,
        clsmap4_p,
        permute4_p,
        epsilon_p,
        max_delta_p,
        converged_code_p,
        convergence_code_p,
        lrxt_p,
        lhet_p,
        lsol_p,
        solution_p,
        iter_invariant_p,
        forcing_p,
    )

@export
def imp_sol_outer_batch_stage_dispatch_codon(
    mode: int,
    i: int,
    lev: int,
    nr_iter: int,
    has_independent: int,
    ncol: int,
    pver: int,
    gas_pcnst: int,
    rxntot: int,
    extcnt: int,
    clscnt4: int,
    dti: float,
    small: float,
    base_sol_p: cobj,
    reaction_rates_p: cobj,
    het_rates_p: cobj,
    extfrc_p: cobj,
    ind_prd_p: cobj,
    clsmap4_p: cobj,
    permute4_p: cobj,
    epsilon_p: cobj,
    max_delta_p: cobj,
    converged_code_p: cobj,
    convergence_code_p: cobj,
    lrxt_p: cobj,
    lhet_p: cobj,
    lsol_p: cobj,
    solution_p: cobj,
    iter_invariant_p: cobj,
    forcing_p: cobj,
):
    return imp_sol_outer_batch_codon(
        mode,
        i,
        lev,
        nr_iter,
        has_independent,
        ncol,
        pver,
        gas_pcnst,
        rxntot,
        extcnt,
        clscnt4,
        dti,
        small,
        base_sol_p,
        reaction_rates_p,
        het_rates_p,
        extfrc_p,
        ind_prd_p,
        clsmap4_p,
        permute4_p,
        epsilon_p,
        max_delta_p,
        converged_code_p,
        convergence_code_p,
        lrxt_p,
        lhet_p,
        lsol_p,
        solution_p,
        iter_invariant_p,
        forcing_p,
    )

@export
def indprd_codon(
    class_id: int,
    ncol: int,
    pver: int,
    nprod: int,
    rxt_p: cobj,
    extfrc_p: cobj,
    prod_p: cobj,
):
    return _gas_phase.indprd_codon(
        class_id,
        ncol,
        pver,
        nprod,
        rxt_p,
        extfrc_p,
        prod_p,
    )

@export
def negtrc_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    fld_p: cobj,
):
    return _gas_phase.negtrc_codon(
        ncol,
        pver,
        gas_pcnst,
        fld_p,
    )

@export
def O1D_to_2OH_adj_codon(
    ncol: int,
    pcols: int,
    pver: int,
    rxntot: int,
    nfs: int,
    jo1d_ndx: int,
    n2_ndx: int,
    o2_ndx: int,
    h2o_ndx: int,
    p_rate_p: cobj,
    inv_p: cobj,
    tfld_p: cobj,
):
    return _gas_phase.O1D_to_2OH_adj_codon(
        ncol,
        pcols,
        pver,
        rxntot,
        nfs,
        jo1d_ndx,
        n2_ndx,
        o2_ndx,
        h2o_ndx,
        p_rate_p,
        inv_p,
        tfld_p,
    )

@export
def spe_init_codon(tag: int) -> int:
    return tag

@export
def init_tracer_cnst_restart_codon(tag: int) -> int:
    return tag

@export
def write_tracer_cnst_restart_codon(tag: int) -> int:
    return tag

@export
def init_tracer_srcs_restart_codon(tag: int) -> int:
    return tag

@export
def write_tracer_srcs_restart_codon(tag: int) -> int:
    return tag

@export
def prescribed_ozone_register_codon(tag: int) -> int:
    return tag

@export
def prescribed_volcaero_register_codon(tag: int) -> int:
    return tag

@export
def get_model_time_codon(tag: int) -> int:
    return tag

@export
def solar_data_get_model_time_codon(
    time_value: float,
    yr: int,
    mn: int,
    dy: int,
    sc: int,
    has_year: int,
    has_month: int,
    has_day: int,
    has_seconds: int,
    time_p: cobj,
    year_p: cobj,
    month_p: cobj,
    day_p: cobj,
    seconds_p: cobj,
):
    time_out = Ptr[float](time_p)
    time_out[0] = time_value

    if has_year != 0:
        year_out = Ptr[int](year_p)
        year_out[0] = yr
    if has_month != 0:
        month_out = Ptr[int](month_p)
        month_out[0] = mn
    if has_day != 0:
        day_out = Ptr[int](day_p)
        day_out[0] = dy
    if has_seconds != 0:
        seconds_out = Ptr[int](seconds_p)
        seconds_out[0] = sc

@export
def convert_date_codon(tag: int) -> int:
    return tag

@export
def modal_aero_wateruptake_reg_codon(tag: int) -> int:
    return tag

@export
def exp_sol_inti_codon(tag: int) -> int:
    return tag

@export
def set_short_lived_species_codon(tag: int) -> int:
    return tag

@export
def convert_dates_codon(tag: int) -> int:
    return tag

@export
def aerodep_flx_adv_codon(tag: int) -> int:
    return tag

@export
def initialize_short_lived_species_codon(tag: int) -> int:
    return tag

@export
def solar_parms_init_codon(tag: int) -> int:
    return tag

@export
def solar_parms_get_codon(tag: int) -> int:
    return tag

@export
def solar_parms_timestep_init_codon(tag: int) -> int:
    return tag

@export
def aerodep_flx_init_codon(tag: int) -> int:
    return tag

@export
def sethet_inti_codon(tag: int) -> int:
    return tag

@export
def write_trc_restart_codon(tag: int) -> int:
    return tag

@export
def seasalt_init_codon(tag: int) -> int:
    return tag

@export
def get_fld_data_codon(tag: int) -> int:
    return tag

@export
def sad_inti_codon(tag: int) -> int:
    return tag

@export
def rate_diags_readnl_codon(tag: int) -> int:
    return tag

@export
def solar_parms_readnl_codon(tag: int) -> int:
    return tag

@export
def tracer_srcs_init_codon(tag: int) -> int:
    return tag

@export
def dust_readnl_codon(tag: int) -> int:
    return tag

@export
def modal_aero_bcscavcoef_init_codon(tag: int) -> int:
    return tag

@export
def advance_trcdata_codon(tag: int) -> int:
    return tag

@export
def init_trc_restart_codon(tag: int) -> int:
    return tag

@inline
def _tracer_idx2(i: int, k: int, ld1: int) -> int:
    return (i - 1) + (k - 1) * ld1

@export
def vert_interp_codon(
    ncol: int,
    pcols: int,
    pver: int,
    levsiz: int,
    pin_p: cobj,
    pmid_p: cobj,
    datain_p: cobj,
    dataout_p: cobj,
    kupper_p: cobj,
):
    pin = Ptr[float](pin_p)
    pmid = Ptr[float](pmid_p)
    datain = Ptr[float](datain_p)
    dataout = Ptr[float](dataout_p)
    kupper = Ptr[int](kupper_p)

    for i in range(1, ncol + 1):
        kupper[i - 1] = 1

    for k in range(1, pver + 1):
        kkstart = levsiz
        for i in range(1, ncol + 1):
            if kupper[i - 1] < kkstart:
                kkstart = kupper[i - 1]

        for kk in range(kkstart, levsiz):
            for i in range(1, ncol + 1):
                if (
                    pin[_tracer_idx2(i, kk, pcols)] < pmid[_tracer_idx2(i, k, pcols)]
                    and pmid[_tracer_idx2(i, k, pcols)] <= pin[_tracer_idx2(i, kk + 1, pcols)]
                ):
                    kupper[i - 1] = kk

        for i in range(1, ncol + 1):
            if pmid[_tracer_idx2(i, k, pcols)] < pin[_tracer_idx2(i, 1, pcols)]:
                dataout[_tracer_idx2(i, k, pcols)] = (
                    datain[_tracer_idx2(i, 1, pcols)]
                    * pmid[_tracer_idx2(i, k, pcols)]
                    / pin[_tracer_idx2(i, 1, pcols)]
                )
            elif pmid[_tracer_idx2(i, k, pcols)] > pin[_tracer_idx2(i, levsiz, pcols)]:
                dataout[_tracer_idx2(i, k, pcols)] = datain[_tracer_idx2(i, levsiz, pcols)]
            else:
                kup = int(kupper[i - 1])
                dpu = pmid[_tracer_idx2(i, k, pcols)] - pin[_tracer_idx2(i, kup, pcols)]
                dpl = pin[_tracer_idx2(i, kup + 1, pcols)] - pmid[_tracer_idx2(i, k, pcols)]
                dataout[_tracer_idx2(i, k, pcols)] = (
                    datain[_tracer_idx2(i, kup, pcols)] * dpl
                    + datain[_tracer_idx2(i, kup + 1, pcols)] * dpu
                ) / (dpl + dpu)

@export
def get_dimension_codon(tag: int) -> int:
    return tag

@export
def ubc_defaultopts_codon(tag: int) -> int:
    return tag

@export
def ubc_setopts_codon(tag: int) -> int:
    return tag

@export
def modal_aero_wateruptake_init_codon(tag: int) -> int:
    return tag

@export
def specify_fields_codon(tag: int) -> int:
    return tag

@export
def sslt_sections_init_codon(tag: int) -> int:
    return tag

@export
def cldaero_uptakerate_codon(
    xl: float,
    cldnum: float,
    cfact: float,
    cldfrc: float,
    tfld: float,
    press: float,
    pi_val: float,
) -> float:
    return _wetchem._sox_cldaero_uptakerate(xl, cldnum, cfact, cldfrc, tfld, press, pi_val)

@export
def fluxes_codon(
    ncol: int,
    nsections: int,
    sst_p: cobj,
    u10cubed_p: cobj,
    consta_p: cobj,
    constb_p: cobj,
    fi_p: cobj,
):
    sst = Ptr[float](sst_p)
    u10cubed = Ptr[float](u10cubed_p)
    consta = Ptr[float](consta_p)
    constb = Ptr[float](constb_p)
    fi = Ptr[float](fi_p)

    for m0 in range(nsections):
        row = 3
        if m0 <= 8:
            row = 0
        elif m0 <= 12:
            row = 1
        elif m0 <= 20:
            row = 2
        for i in range(ncol):
            if m0 <= 20:
                w = (3.84e-6 * u10cubed[i]) * 0.1
                fi[i + m0 * ncol] = w * ((sst[i] * consta[row + m0 * 4]) + constb[row + m0 * 4])
            else:
                fi[i + m0 * ncol] = consta[row + m0 * 4] * u10cubed[i]

@export
def modal_aero_newnuc_init_codon(tag: int) -> int:
    return tag

@export
def setext_inti_codon(tag: int) -> int:
    return tag

@export
def solar_data_readnl_codon(tag: int) -> int:
    return tag

@export
def sulf_readnl_codon(tag: int) -> int:
    return tag

@export
def aerodep_flx_readnl_codon(tag: int) -> int:
    return tag

@export
def prescribed_ghg_readnl_codon(tag: int) -> int:
    return tag

@export
def prescribed_strataero_readnl_codon(tag: int) -> int:
    return tag

@export
def dust_set_params_codon(tag: int) -> int:
    return tag

@export
def set_sim_dat_codon(tag: int) -> int:
    return tag

@export
def soil_erod_init_codon(tag: int) -> int:
    return tag
