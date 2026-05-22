from C import uwshcu_compute_native_from_c_cb(int, int, int, int, float, Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], int, Ptr[float], Ptr[float], Ptr[float], Ptr[float], int, Ptr[float], Ptr[float], int, Ptr[int], int, Ptr[int], int, Ptr[int], Ptr[int], int) -> None
from C import uwshcu_conden_scalar_from_c_cb(float, float, float, Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[int], int) -> None
from C import uwshcu_top_conden_from_c_cb(int, int, int, float, float, float, float, float, Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[int], Ptr[float], Ptr[float], Ptr[float]) -> None
from C import uwshcu_thermo_conden_from_c_cb(int, int, int, int, int, int, float, float, float, Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[int], Ptr[float], Ptr[float], Ptr[float]) -> None
from C import uwshcu_cnst_indices_from_c_cb(Ptr[int]) -> None
from C import uwshcu_findsp_layer_from_c_cb(int, Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float]) -> None
from C import uwshcu_qsinvert_from_c_cb(float, float, float) -> float
from C import uwshcu_qsat_from_c_cb(float, float, Ptr[float], Ptr[float]) -> None
from C import uwshcu_select_init_shell_from_c_cb(Ptr[int]) -> None
from C import uwshcu_wtrc_metadata_from_c_cb(Ptr[int], Ptr[int]) -> None
from C import uwshcu_wtrc_ratio_type_from_c_cb(int, float, float) -> float
from C import uwshcu_positive_moisture_single_from_c_cb(int, int, int, float, float, float, float, float, float, Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[float], Ptr[int]) -> None


def uwshcu_compute_native_from_c_dispatch(
    mix: int,
    mkx: int,
    iend: int,
    ncnst: int,
    dt: float,
    ps0_p: cobj,
    zs0_p: cobj,
    p0_p: cobj,
    z0_p: cobj,
    dp0_p: cobj,
    u0_p: cobj,
    v0_p: cobj,
    qv0_p: cobj,
    ql0_p: cobj,
    qi0_p: cobj,
    t0_p: cobj,
    s0_p: cobj,
    tr0_p: cobj,
    tke_p: cobj,
    cldfrct_p: cobj,
    concldfrct_p: cobj,
    pblh_p: cobj,
    cush_p: cobj,
    umf_p: cobj,
    slflx_p: cobj,
    qtflx_p: cobj,
    flxprc1_p: cobj,
    flxsnow1_p: cobj,
    qvten_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    sten_p: cobj,
    uten_p: cobj,
    vten_p: cobj,
    trten_p: cobj,
    qrten_p: cobj,
    qsten_p: cobj,
    precip_p: cobj,
    snow_p: cobj,
    evapc_p: cobj,
    cufrc_p: cobj,
    qcu_p: cobj,
    qlu_p: cobj,
    qiu_p: cobj,
    cbmf_p: cobj,
    qc_p: cobj,
    rliq_p: cobj,
    cnt_p: cobj,
    cnb_p: cobj,
    lchnk: int,
    dpdry0_p: cobj,
    wtprec_p: cobj,
    wtsnow_p: cobj,
    wtqc_p: cobj,
    wetbulb_precomputed: int,
    tw0_precomputed_p: cobj,
    qw0_precomputed_p: cobj,
    constituent_indices_precomputed: int,
    constituent_indices_p: cobj,
    init_shell_preselected: int,
    init_shell_flags_p: cobj,
    wtrc_metadata_precomputed: int,
    wtrc_metadata_flags_p: cobj,
    wtrc_iatype_p: cobj,
    public_outputs_preinitialized: int,
):
    uwshcu_compute_native_from_c_cb(
        mix,
        mkx,
        iend,
        ncnst,
        dt,
        Ptr[float](ps0_p),
        Ptr[float](zs0_p),
        Ptr[float](p0_p),
        Ptr[float](z0_p),
        Ptr[float](dp0_p),
        Ptr[float](u0_p),
        Ptr[float](v0_p),
        Ptr[float](qv0_p),
        Ptr[float](ql0_p),
        Ptr[float](qi0_p),
        Ptr[float](t0_p),
        Ptr[float](s0_p),
        Ptr[float](tr0_p),
        Ptr[float](tke_p),
        Ptr[float](cldfrct_p),
        Ptr[float](concldfrct_p),
        Ptr[float](pblh_p),
        Ptr[float](cush_p),
        Ptr[float](umf_p),
        Ptr[float](slflx_p),
        Ptr[float](qtflx_p),
        Ptr[float](flxprc1_p),
        Ptr[float](flxsnow1_p),
        Ptr[float](qvten_p),
        Ptr[float](qlten_p),
        Ptr[float](qiten_p),
        Ptr[float](sten_p),
        Ptr[float](uten_p),
        Ptr[float](vten_p),
        Ptr[float](trten_p),
        Ptr[float](qrten_p),
        Ptr[float](qsten_p),
        Ptr[float](precip_p),
        Ptr[float](snow_p),
        Ptr[float](evapc_p),
        Ptr[float](cufrc_p),
        Ptr[float](qcu_p),
        Ptr[float](qlu_p),
        Ptr[float](qiu_p),
        Ptr[float](cbmf_p),
        Ptr[float](qc_p),
        Ptr[float](rliq_p),
        Ptr[float](cnt_p),
        Ptr[float](cnb_p),
        lchnk,
        Ptr[float](dpdry0_p),
        Ptr[float](wtprec_p),
        Ptr[float](wtsnow_p),
        Ptr[float](wtqc_p),
        wetbulb_precomputed,
        Ptr[float](tw0_precomputed_p),
        Ptr[float](qw0_precomputed_p),
        constituent_indices_precomputed,
        Ptr[int](constituent_indices_p),
        init_shell_preselected,
        Ptr[int](init_shell_flags_p),
        wtrc_metadata_precomputed,
        Ptr[int](wtrc_metadata_flags_p),
        Ptr[int](wtrc_iatype_p),
        public_outputs_preinitialized,
    )


def uwshcu_cnst_indices_from_c_dispatch(indices_p: cobj):
    uwshcu_cnst_indices_from_c_cb(Ptr[int](indices_p))


def uwshcu_conden_scalar_from_c_dispatch(
    p: float,
    thl: float,
    qt: float,
    th_p: cobj,
    qv_p: cobj,
    ql_p: cobj,
    qi_p: cobj,
    qse_p: cobj,
    id_check_p: cobj,
    ncnst: int,
):
    uwshcu_conden_scalar_from_c_cb(
        p,
        thl,
        qt,
        Ptr[float](th_p),
        Ptr[float](qv_p),
        Ptr[float](ql_p),
        Ptr[float](qi_p),
        Ptr[float](qse_p),
        Ptr[int](id_check_p),
        ncnst,
    )


def uwshcu_top_conden_from_c_dispatch(
    trace_water: int,
    wtrc_nwset: int,
    ncnst: int,
    pressure: float,
    thl: float,
    qt: float,
    p00: float,
    rovcp: float,
    th_p: cobj,
    qv_p: cobj,
    ql_p: cobj,
    qi_p: cobj,
    qse_p: cobj,
    id_check_p: cobj,
    exntop_p: cobj,
    wtu_top_p: cobj,
    wtout_p: cobj,
):
    uwshcu_top_conden_from_c_cb(
        trace_water,
        wtrc_nwset,
        ncnst,
        pressure,
        thl,
        qt,
        p00,
        rovcp,
        Ptr[float](th_p),
        Ptr[float](qv_p),
        Ptr[float](ql_p),
        Ptr[float](qi_p),
        Ptr[float](qse_p),
        Ptr[int](id_check_p),
        Ptr[float](exntop_p),
        Ptr[float](wtu_top_p),
        Ptr[float](wtout_p),
    )


def uwshcu_thermo_conden_from_c_dispatch(
    trace_water: int,
    wtrc_nwset: int,
    ncnst: int,
    mkx: int,
    wtu_row: int,
    use_top: int,
    p: float,
    thl: float,
    qt: float,
    th_p: cobj,
    qv_p: cobj,
    ql_p: cobj,
    qi_p: cobj,
    qse_p: cobj,
    id_check_p: cobj,
    wtu_p: cobj,
    wtu_top_p: cobj,
    wtout_p: cobj,
):
    uwshcu_thermo_conden_from_c_cb(
        trace_water,
        wtrc_nwset,
        ncnst,
        mkx,
        wtu_row,
        use_top,
        p,
        thl,
        qt,
        Ptr[float](th_p),
        Ptr[float](qv_p),
        Ptr[float](ql_p),
        Ptr[float](qi_p),
        Ptr[float](qse_p),
        Ptr[int](id_check_p),
        Ptr[float](wtu_p),
        Ptr[float](wtu_top_p),
        Ptr[float](wtout_p),
    )


def uwshcu_select_init_shell_from_c_dispatch(flags_p: cobj):
    uwshcu_select_init_shell_from_c_cb(Ptr[int](flags_p))


def uwshcu_wtrc_metadata_from_c_dispatch(flags_p: cobj, iatype_p: cobj):
    uwshcu_wtrc_metadata_from_c_cb(Ptr[int](flags_p), Ptr[int](iatype_p))


def uwshcu_findsp_layer_from_c_dispatch(
    iend: int,
    qv0_p: Ptr[float],
    t0_p: Ptr[float],
    p0_p: Ptr[float],
    tw0_p: Ptr[float],
    qw0_p: Ptr[float],
):
    uwshcu_findsp_layer_from_c_cb(iend, qv0_p, t0_p, p0_p, tw0_p, qw0_p)


def uwshcu_qsinvert_from_c_dispatch(qt: float, thl: float, psfc: float) -> float:
    return uwshcu_qsinvert_from_c_cb(qt, thl, psfc)


def uwshcu_qsat_from_c_dispatch(t: float, p: float, es_p: cobj, qs_p: cobj):
    uwshcu_qsat_from_c_cb(t, p, Ptr[float](es_p), Ptr[float](qs_p))


def uwshcu_wtrc_ratio_type_from_c_dispatch(iatype: int, qtrc: float, qtot: float) -> float:
    return uwshcu_wtrc_ratio_type_from_c_cb(iatype, qtrc, qtot)


def uwshcu_positive_moisture_single_from_c_dispatch(
    mkx: int,
    ncnst: int,
    trace_water: int,
    xlv: float,
    xls: float,
    dt: float,
    qvmin: float,
    qlmin: float,
    qimin: float,
    dp_p: cobj,
    qv_p: cobj,
    ql_p: cobj,
    qi_p: cobj,
    s_p: cobj,
    qvten_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    sten_p: cobj,
    wtr_p: cobj,
    wtten_p: cobj,
    status_p: cobj,
):
    uwshcu_positive_moisture_single_from_c_cb(
        mkx,
        ncnst,
        trace_water,
        xlv,
        xls,
        dt,
        qvmin,
        qlmin,
        qimin,
        Ptr[float](dp_p),
        Ptr[float](qv_p),
        Ptr[float](ql_p),
        Ptr[float](qi_p),
        Ptr[float](s_p),
        Ptr[float](qvten_p),
        Ptr[float](qlten_p),
        Ptr[float](qiten_p),
        Ptr[float](sten_p),
        Ptr[float](wtr_p),
        Ptr[float](wtten_p),
        Ptr[int](status_p),
    )
