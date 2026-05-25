from math import exp, log10
from C import cldwat2m_qsat_water_native_cb(float, float, Ptr[float], Ptr[float], Ptr[float]) -> None
from C import cldwat2m_astg_single_native_cb(
    int,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    Ptr[float],
    Ptr[float],
) -> None
from C import cldwat2m_aist_single_native_cb(
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    Ptr[float],
) -> None


@export
def macrop_driver_readnl_codon(flag: int) -> int:
    return flag


@export
def macrop_driver_register_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def macrop_driver_init_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def cldwat2m_ini_macro_codon(
    rhminl_opt: int,
    rhmini_opt: int,
    rhminl_in: float,
    rhminl_adj_land_in: float,
    rhminh_in: float,
    premit_in: float,
    premib_in: float,
    i_rhminl_p: cobj,
    i_rhmini_p: cobj,
    rhminl_p: cobj,
    rhminl_adj_land_p: cobj,
    rhminh_p: cobj,
    premit_p: cobj,
    premib_p: cobj,
):
    i_rhminl = Ptr[int](i_rhminl_p)
    i_rhmini = Ptr[int](i_rhmini_p)
    rhminl = Ptr[float](rhminl_p)
    rhminl_adj_land = Ptr[float](rhminl_adj_land_p)
    rhminh = Ptr[float](rhminh_p)
    premit = Ptr[float](premit_p)
    premib = Ptr[float](premib_p)

    i_rhminl[0] = rhminl_opt
    i_rhmini[0] = rhmini_opt
    rhminl[0] = rhminl_in
    rhminl_adj_land[0] = rhminl_adj_land_in
    rhminh[0] = rhminh_in
    premit[0] = premit_in
    premib[0] = premib_in


def _idx2(i: int, k: int, pcols: int):
    return (k - 1) * pcols + (i - 1)


def _idx3(i: int, k: int, m: int, pcols: int, pver: int):
    return (m - 1) * pcols * pver + (k - 1) * pcols + (i - 1)


@inline
def _cldwat2m_qsat_water_native(t: float, p: float):
    es = 0.0
    qs = 0.0
    dqsdt = 0.0
    cldwat2m_qsat_water_native_cb(t, p, __ptr__(es), __ptr__(qs), __ptr__(dqsdt))
    return es, qs, dqsdt


@inline
def _cldwat2m_astg_single_native(
    use_rhu: int,
    u: float,
    p: float,
    qv: float,
    landfrac: float,
    snowh: float,
    rhminl: float,
    rhminl_adj_land: float,
    rhminh: float,
):
    a = 0.0
    ga = 0.0
    cldwat2m_astg_single_native_cb(
        use_rhu,
        u,
        p,
        qv,
        landfrac,
        snowh,
        rhminl,
        rhminl_adj_land,
        rhminh,
        __ptr__(a),
        __ptr__(ga),
    )
    return a, ga


@inline
def _cldwat2m_aist_single_native(
    qv: float,
    t: float,
    p: float,
    qi: float,
    landfrac: float,
    snowh: float,
    rhmaxi: float,
    rhmini: float,
    rhminl: float,
    rhminl_adj_land: float,
    rhminh: float,
):
    aist = 0.0
    cldwat2m_aist_single_native_cb(
        qv,
        t,
        p,
        qi,
        landfrac,
        snowh,
        rhmaxi,
        rhmini,
        rhminl,
        rhminl_adj_land,
        rhminh,
        __ptr__(aist),
    )
    return aist


@export
def cldwat2m_gridmean_rh_codon(
    p: float,
    latvap: float,
    cpair: float,
    t_p: cobj,
    qv_p: cobj,
    ql_p: cobj,
    qi_p: cobj,
    a_dc: float,
    ql_dc: float,
    qi_dc: float,
    a_sc: float,
    ql_sc: float,
    qi_sc: float,
):
    t = Ptr[float](t_p)
    qv = Ptr[float](qv_p)
    ql = Ptr[float](ql_p)
    qi = Ptr[float](qi_p)

    _cldwat2m_gridmean_rh_calc(
        p,
        latvap,
        cpair,
        t,
        qv,
        ql,
        qi,
        a_dc,
        ql_dc,
        qi_dc,
        a_sc,
        ql_sc,
        qi_sc,
    )


@inline
def _cldwat2m_gridmean_rh_calc(
    p: float,
    latvap: float,
    cpair: float,
    t: Ptr[float],
    qv: Ptr[float],
    ql: Ptr[float],
    qi: Ptr[float],
    a_dc: float,
    ql_dc: float,
    qi_dc: float,
    a_sc: float,
    ql_sc: float,
    qi_sc: float,
):
    ql_nc0 = max(0.0, ql[0] - a_dc * ql_dc - a_sc * ql_sc)
    qi_nc0 = max(0.0, qi[0] - a_dc * qi_dc - a_sc * qi_sc)
    qc_nc0 = max(0.0, ql[0] + qi[0] - a_dc * (ql_dc + qi_dc) - a_sc * (ql_sc + qi_sc))
    tc = t[0] - (latvap / cpair) * ql[0]
    qt = qv[0] + ql[0]

    for _ in range(20):
        _, qs, dqsdt = _cldwat2m_qsat_water_native(t[0], p)
        tscale = latvap / cpair
        qc = (t[0] - tc) / tscale
        dqcdt = 1.0 / tscale
        f = qs + qc - qt
        fg = dqsdt + dqcdt
        fg_abs = abs(fg)
        if fg_abs < 1.0e-10:
            fg_abs = 1.0e-10
        if fg < 0.0:
            fg = -fg_abs
        else:
            fg = fg_abs
        if qc >= 0.0 and (qt - qc) >= 0.999 * qs and (qt - qc) <= 1.0 * qs:
            break
        t[0] = t[0] - f / fg

    _, qs_final, _ = _cldwat2m_qsat_water_native(t[0], p)
    qv[0] = min(qt, qs_final)
    ql[0] = qt - qv[0]
    t[0] = tc + (latvap / cpair) * ql[0]


@export
def cldwat2m_funcd_instratus_codon(
    t: float,
    p: float,
    t0: float,
    qv0: float,
    ql0: float,
    qi0: float,
    fice0: float,
    muq0: float,
    qc_nc0: float,
    a_dc: float,
    ql_dc: float,
    qi_dc: float,
    a_sc: float,
    ql_sc: float,
    qi_sc: float,
    ai_st: float,
    qcst_crit: float,
    landfrac: float,
    snowh: float,
    rhminl: float,
    rhminl_adj_land: float,
    rhminh: float,
    cpair: float,
    latvap: float,
    camstfrac: int,
    f_p: cobj,
    fg_p: cobj,
    qc_nc_p: cobj,
    fice_p: cobj,
    al_st_p: cobj,
):
    f = Ptr[float](f_p)
    fg = Ptr[float](fg_p)
    qc_nc_out = Ptr[float](qc_nc_p)
    fice_out = Ptr[float](fice_p)
    al_st_out = Ptr[float](al_st_p)

    _cldwat2m_funcd_instratus_calc(
        t,
        p,
        t0,
        qv0,
        ql0,
        qi0,
        fice0,
        muq0,
        qc_nc0,
        a_dc,
        ql_dc,
        qi_dc,
        a_sc,
        ql_sc,
        qi_sc,
        ai_st,
        qcst_crit,
        landfrac,
        snowh,
        rhminl,
        rhminl_adj_land,
        rhminh,
        cpair,
        latvap,
        camstfrac,
        f,
        fg,
        qc_nc_out,
        fice_out,
        al_st_out,
    )


@inline
def _cldwat2m_funcd_instratus_calc(
    t: float,
    p: float,
    t0: float,
    qv0: float,
    ql0: float,
    qi0: float,
    fice0: float,
    muq0: float,
    qc_nc0: float,
    a_dc: float,
    ql_dc: float,
    qi_dc: float,
    a_sc: float,
    ql_sc: float,
    qi_sc: float,
    ai_st: float,
    qcst_crit: float,
    landfrac: float,
    snowh: float,
    rhminl: float,
    rhminl_adj_land: float,
    rhminh: float,
    cpair: float,
    latvap: float,
    camstfrac: int,
    f_out: Ptr[float],
    fg_out: Ptr[float],
    qc_nc_out: Ptr[float],
    fice_out: Ptr[float],
    al_st_out: Ptr[float],
):
    _, qs, dqsdt = _cldwat2m_qsat_water_native(t, p)

    fice = fice0
    qc_nc = (cpair / latvap) * (t - t0) + muq0 * qc_nc0
    dqcncdt = cpair / latvap
    qv = qv0 + ql0 + qi0 - (qc_nc + a_dc * (ql_dc + qi_dc) + a_sc * (ql_sc + qi_sc))
    alpha = 1.0 / qs
    beta = (qv / qs**2.0) * dqsdt

    u = qv / qs
    u_nc = u
    al_st_nc, g_nc = _cldwat2m_astg_single_native(
        camstfrac,
        u_nc,
        p,
        qv,
        landfrac,
        snowh,
        rhminl,
        rhminl_adj_land,
        rhminh,
    )
    al_st = (1.0 - a_dc - a_sc) * al_st_nc
    dudt = -(alpha * dqcncdt + beta)
    dalstdt = (1.0 / g_nc) * dudt
    if u_nc == 1.0:
        dalstdt = 0.0

    f_out[0] = qc_nc - qcst_crit * al_st
    fg_out[0] = dqcncdt - qcst_crit * dalstdt
    qc_nc_out[0] = qc_nc
    fice_out[0] = fice
    al_st_out[0] = al_st


@export
def cldwat2m_instratus_core_codon(
    p: float,
    t0: float,
    qv0: float,
    ql0: float,
    qi0: float,
    a_dc: float,
    ql_dc: float,
    qi_dc: float,
    a_sc: float,
    ql_sc: float,
    qi_sc: float,
    ai_st: float,
    qcst_crit: float,
    tmin: float,
    tmax: float,
    landfrac: float,
    snowh: float,
    rhminl: float,
    rhminl_adj_land: float,
    rhminh: float,
    cpair: float,
    latvap: float,
    qlst_min: float,
    qlst_max: float,
    camstfrac: int,
    t_p: cobj,
    qv_p: cobj,
    ql_p: cobj,
    qi_p: cobj,
):
    t_out = Ptr[float](t_p)
    qv_out = Ptr[float](qv_p)
    ql_out = Ptr[float](ql_p)
    qi_out = Ptr[float](qi_p)

    _cldwat2m_instratus_core_calc(
        p,
        t0,
        qv0,
        ql0,
        qi0,
        a_dc,
        ql_dc,
        qi_dc,
        a_sc,
        ql_sc,
        qi_sc,
        ai_st,
        qcst_crit,
        tmin,
        tmax,
        landfrac,
        snowh,
        rhminl,
        rhminl_adj_land,
        rhminh,
        cpair,
        latvap,
        qlst_min,
        qlst_max,
        camstfrac,
        t_out,
        qv_out,
        ql_out,
        qi_out,
    )


@inline
def _cldwat2m_instratus_core_calc(
    p: float,
    t0: float,
    qv0: float,
    ql0: float,
    qi0: float,
    a_dc: float,
    ql_dc: float,
    qi_dc: float,
    a_sc: float,
    ql_sc: float,
    qi_sc: float,
    ai_st: float,
    qcst_crit: float,
    tmin: float,
    tmax: float,
    landfrac: float,
    snowh: float,
    rhminl: float,
    rhminl_adj_land: float,
    rhminh: float,
    cpair: float,
    latvap: float,
    qlst_min: float,
    qlst_max: float,
    camstfrac: int,
    t_out: Ptr[float],
    qv_out: Ptr[float],
    ql_out: Ptr[float],
    qi_out: Ptr[float],
):
    ql_nc0 = max(0.0, ql0 - a_dc * ql_dc - a_sc * ql_sc)
    qi_nc0 = max(0.0, qi0 - a_dc * qi_dc - a_sc * qi_sc)
    qc_nc0 = max(0.0, ql0 + qi0 - a_dc * (ql_dc + qi_dc) - a_sc * (ql_sc + qi_sc))
    fice0 = 0.0
    ficeg0 = 0.0
    muq0 = 1.0

    df = 0.0
    f = 0.0
    fh = 0.0
    fl = 0.0
    qc_nc = 0.0
    fice = 0.0
    al_st = 0.0

    x1 = tmin
    x2 = tmax
    _cldwat2m_funcd_instratus_calc(
        x1,
        p,
        t0,
        qv0,
        ql0,
        qi0,
        fice0,
        muq0,
        qc_nc0,
        a_dc,
        ql_dc,
        qi_dc,
        a_sc,
        ql_sc,
        qi_sc,
        ai_st,
        qcst_crit,
        landfrac,
        snowh,
        rhminl,
        rhminl_adj_land,
        rhminh,
        cpair,
        latvap,
        camstfrac,
        __ptr__(fl),
        __ptr__(df),
        __ptr__(qc_nc),
        __ptr__(fice),
        __ptr__(al_st),
    )
    _cldwat2m_funcd_instratus_calc(
        x2,
        p,
        t0,
        qv0,
        ql0,
        qi0,
        fice0,
        muq0,
        qc_nc0,
        a_dc,
        ql_dc,
        qi_dc,
        a_sc,
        ql_sc,
        qi_sc,
        ai_st,
        qcst_crit,
        landfrac,
        snowh,
        rhminl,
        rhminl_adj_land,
        rhminh,
        cpair,
        latvap,
        camstfrac,
        __ptr__(fh),
        __ptr__(df),
        __ptr__(qc_nc),
        __ptr__(fice),
        __ptr__(al_st),
    )

    rtsafe = 0.0
    if (fl > 0.0 and fh > 0.0) or (fl < 0.0 and fh < 0.0):
        _cldwat2m_funcd_instratus_calc(
            t0,
            p,
            t0,
            qv0,
            ql0,
            qi0,
            fice0,
            muq0,
            qc_nc0,
            a_dc,
            ql_dc,
            qi_dc,
            a_sc,
            ql_sc,
            qi_sc,
            ai_st,
            qcst_crit,
            landfrac,
            snowh,
            rhminl,
            rhminl_adj_land,
            rhminh,
            cpair,
            latvap,
            camstfrac,
            __ptr__(fl),
            __ptr__(df),
            __ptr__(qc_nc),
            __ptr__(fice),
            __ptr__(al_st),
        )
        rtsafe = t0
    elif fl == 0.0:
        rtsafe = x1
    elif fh == 0.0:
        rtsafe = x2
    else:
        if fl < 0.0:
            xl = x1
            xh = x2
        else:
            xh = x1
            xl = x2
        rtsafe = 0.5 * (x1 + x2)
        dxold = abs(x2 - x1)
        dx = dxold
        _cldwat2m_funcd_instratus_calc(
            rtsafe,
            p,
            t0,
            qv0,
            ql0,
            qi0,
            fice0,
            muq0,
            qc_nc0,
            a_dc,
            ql_dc,
            qi_dc,
            a_sc,
            ql_sc,
            qi_sc,
            ai_st,
            qcst_crit,
            landfrac,
            snowh,
            rhminl,
            rhminl_adj_land,
            rhminh,
            cpair,
            latvap,
            camstfrac,
            __ptr__(f),
            __ptr__(df),
            __ptr__(qc_nc),
            __ptr__(fice),
            __ptr__(al_st),
        )
        for _ in range(20):
            if ((rtsafe - xh) * df - f) * ((rtsafe - xl) * df - f) > 0.0 or abs(2.0 * f) > abs(dxold * df):
                dxold = dx
                dx = 0.5 * (xh - xl)
                rtsafe = xl + dx
                if xl == rtsafe:
                    break
            else:
                dxold = dx
                dx = f / df
                temp = rtsafe
                rtsafe = rtsafe - dx
                if temp == rtsafe:
                    break
            _cldwat2m_funcd_instratus_calc(
                rtsafe,
                p,
                t0,
                qv0,
                ql0,
                qi0,
                fice0,
                muq0,
                qc_nc0,
                a_dc,
                ql_dc,
                qi_dc,
                a_sc,
                ql_sc,
                qi_sc,
                ai_st,
                qcst_crit,
                landfrac,
                snowh,
                rhminl,
                rhminl_adj_land,
                rhminh,
                cpair,
                latvap,
                camstfrac,
                __ptr__(f),
                __ptr__(df),
                __ptr__(qc_nc),
                __ptr__(fice),
                __ptr__(al_st),
            )
            if qcst_crit < 0.5 * (qlst_min + qlst_max):
                if qc_nc * (1.0 - fice) > qlst_min * al_st and qc_nc * (1.0 - fice) < 1.1 * qlst_min * al_st:
                    break
            else:
                if qc_nc * (1.0 - fice) > 0.9 * qlst_max * al_st and qc_nc * (1.0 - fice) < qlst_max * al_st:
                    break
            if f < 0.0:
                xl = rtsafe
            else:
                xh = rtsafe

    qc_nc = max(0.0, qc_nc)

    t_out[0] = rtsafe
    ql_out[0] = qc_nc * (1.0 - fice) + a_dc * ql_dc + a_sc * ql_sc
    qi_out[0] = qc_nc * fice + a_dc * qi_dc + a_sc * qi_sc
    qv_out[0] = qv0 + ql0 + qi0 - (qc_nc + a_dc * (ql_dc + qi_dc) + a_sc * (ql_sc + qi_sc))
    qv_out[0] = max(qv_out[0], 1.0e-12)


@export
def cldwat2m_instratus_condensate_codon(
    ncol: int,
    pcols: int,
    camstfrac: int,
    cpair: float,
    latvap: float,
    latice: float,
    qlst_min: float,
    qlst_max: float,
    rhmaxi: float,
    p_p: cobj,
    t0_p: cobj,
    qv0_p: cobj,
    ql0_p: cobj,
    qi0_p: cobj,
    ni0_p: cobj,
    a_dc_p: cobj,
    ql_dc_p: cobj,
    qi_dc_p: cobj,
    a_sc_p: cobj,
    ql_sc_p: cobj,
    qi_sc_p: cobj,
    landfrac_p: cobj,
    snowh_p: cobj,
    rhmini_p: cobj,
    rhminl_p: cobj,
    rhminl_adj_land_p: cobj,
    rhminh_p: cobj,
    t_out_p: cobj,
    qv_out_p: cobj,
    ql_out_p: cobj,
    qi_out_p: cobj,
    al_st_out_p: cobj,
    ai_st_out_p: cobj,
    ql_st_out_p: cobj,
    qi_st_out_p: cobj,
    status_p: cobj,
):
    p_in = Ptr[float](p_p)
    t0_in = Ptr[float](t0_p)
    qv0_in = Ptr[float](qv0_p)
    ql0_in = Ptr[float](ql0_p)
    qi0_in = Ptr[float](qi0_p)
    ni0_in = Ptr[float](ni0_p)
    a_dc_in = Ptr[float](a_dc_p)
    ql_dc_in = Ptr[float](ql_dc_p)
    qi_dc_in = Ptr[float](qi_dc_p)
    a_sc_in = Ptr[float](a_sc_p)
    ql_sc_in = Ptr[float](ql_sc_p)
    qi_sc_in = Ptr[float](qi_sc_p)
    landfrac = Ptr[float](landfrac_p)
    snowh = Ptr[float](snowh_p)
    rhmini_in = Ptr[float](rhmini_p)
    rhminl_in = Ptr[float](rhminl_p)
    rhminl_adj_land_in = Ptr[float](rhminl_adj_land_p)
    rhminh_in = Ptr[float](rhminh_p)
    t_out = Ptr[float](t_out_p)
    qv_out = Ptr[float](qv_out_p)
    ql_out = Ptr[float](ql_out_p)
    qi_out = Ptr[float](qi_out_p)
    al_st_out = Ptr[float](al_st_out_p)
    ai_st_out = Ptr[float](ai_st_out_p)
    ql_st_out = Ptr[float](ql_st_out_p)
    qi_st_out = Ptr[float](qi_st_out_p)
    status = Ptr[int](status_p)

    status[0] = 0
    for i in range(1, ncol + 1):
        idx = i - 1
        p = p_in[idx]

        t0 = t0_in[idx]
        qv0 = qv0_in[idx]
        ql0 = ql0_in[idx]
        qi0 = qi0_in[idx]

        a_dc = a_dc_in[idx]
        ql_dc = ql_dc_in[idx]
        qi_dc = qi_dc_in[idx]

        a_sc = a_sc_in[idx]
        ql_sc = ql_sc_in[idx]
        qi_sc = qi_sc_in[idx]

        ql_dc = 0.0
        qi_dc = 0.0
        ql_sc = 0.0
        qi_sc = 0.0

        _, qs, _ = _cldwat2m_qsat_water_native(t0, p)

        rhmini = rhmini_in[idx]
        rhminl = rhminl_in[idx]
        rhminl_adj_land = rhminl_adj_land_in[idx]
        rhminh = rhminh_in[idx]

        idxmod = 0

        u0 = qv0 / qs
        u0_nc = u0
        al0_st_nc, g0_nc = _cldwat2m_astg_single_native(
            camstfrac,
            u0_nc,
            p,
            qv0,
            landfrac[idx],
            snowh[idx],
            rhminl,
            rhminl_adj_land,
            rhminh,
        )
        ai0_st_nc = _cldwat2m_aist_single_native(
            qv0,
            t0,
            p,
            qi0,
            landfrac[idx],
            snowh[idx],
            rhmaxi,
            rhmini,
            rhminl,
            rhminl_adj_land,
            rhminh,
        )

        if qv0 > qs:
            _cldwat2m_gridmean_rh_calc(
                p,
                latvap,
                cpair,
                __ptr__(t0),
                __ptr__(qv0),
                __ptr__(ql0),
                __ptr__(qi0),
                a_dc,
                ql_dc,
                qi_dc,
                a_sc,
                ql_sc,
                qi_sc,
            )
            _, qsat0, _ = _cldwat2m_qsat_water_native(t0, p)
            u0 = qv0 / qsat0
            u0_nc = u0
            al0_st_nc, g0_nc = _cldwat2m_astg_single_native(
                camstfrac,
                u0_nc,
                p,
                qv0,
                landfrac[idx],
                snowh[idx],
                rhminl,
                rhminl_adj_land,
                rhminh,
            )
            ai0_st_nc = _cldwat2m_aist_single_native(
                qv0,
                t0,
                p,
                qi0,
                landfrac[idx],
                snowh[idx],
                rhmaxi,
                rhmini,
                rhminl,
                rhminl_adj_land,
                rhminh,
            )
            ai0_st = (1.0 - a_dc - a_sc) * ai0_st_nc
            al0_st = (1.0 - a_dc - a_sc) * al0_st_nc
            a0_st = max(ai0_st, al0_st)
            idxmod = 1
        else:
            ai0_st = (1.0 - a_dc - a_sc) * ai0_st_nc
            al0_st = (1.0 - a_dc - a_sc) * al0_st_nc
        a0_st = max(ai0_st, al0_st)

        ql0_nc = max(0.0, ql0 - a_dc * ql_dc - a_sc * ql_sc)
        qi0_nc = max(0.0, qi0 - a_dc * qi_dc - a_sc * qi_sc)
        qc0_nc = ql0_nc + qi0_nc

        tmin0 = t0 - (latvap / cpair) * ql0
        tmax0 = t0 + ((latvap + latice) / cpair) * qv0

        t = 0.0
        qv = 0.0
        ql = 0.0
        qi = 0.0
        al_st = 0.0
        ai_st = 0.0
        ql_st = 0.0
        qi_st = 0.0

        if ql0_nc >= qlst_min * al0_st and ql0_nc <= qlst_max * al0_st:
            t = t0
            qv = qv0
            ql = ql0
            qi = qi0
        else:
            if al0_st == 0.0 and ql0_nc > 0.0:
                t = tmin0
                qv = qv0 + ql0
                _, qs, _ = _cldwat2m_qsat_water_native(t, p)
                u = qv / qs
                u_nc = u
                al_st_nc, g_nc = _cldwat2m_astg_single_native(
                    camstfrac,
                    u_nc,
                    p,
                    qv,
                    landfrac[idx],
                    snowh[idx],
                    rhminl,
                    rhminl_adj_land,
                    rhminh,
                )
                al_st = (1.0 - a_dc - a_sc) * al_st_nc

                if al_st == 0.0:
                    ql = 0.0
                    qi = qi0
                    idxmod = 1
                else:
                    _cldwat2m_instratus_core_calc(
                        p,
                        t0,
                        qv0,
                        ql0,
                        0.0,
                        a_dc,
                        ql_dc,
                        qi_dc,
                        a_sc,
                        ql_sc,
                        qi_sc,
                        ai0_st,
                        qlst_max,
                        tmin0,
                        t0,
                        landfrac[idx],
                        snowh[idx],
                        rhminl,
                        rhminl_adj_land,
                        rhminh,
                        cpair,
                        latvap,
                        qlst_min,
                        qlst_max,
                        camstfrac,
                        __ptr__(t),
                        __ptr__(qv),
                        __ptr__(ql),
                        __ptr__(qi),
                    )
                    idxmod = 1
            elif al0_st > 0.0 and ql0_nc == 0.0:
                _cldwat2m_instratus_core_calc(
                    p,
                    t0,
                    qv0,
                    ql0,
                    0.0,
                    a_dc,
                    ql_dc,
                    qi_dc,
                    a_sc,
                    ql_sc,
                    qi_sc,
                    ai0_st,
                    qlst_min,
                    tmin0,
                    tmax0,
                    landfrac[idx],
                    snowh[idx],
                    rhminl,
                    rhminl_adj_land,
                    rhminh,
                    cpair,
                    latvap,
                    qlst_min,
                    qlst_max,
                    camstfrac,
                    __ptr__(t),
                    __ptr__(qv),
                    __ptr__(ql),
                    __ptr__(qi),
                )
                idxmod = 1
            elif al0_st > 0.0 and ql0_nc > 0.0:
                if ql0_nc > qlst_max * al0_st:
                    _cldwat2m_instratus_core_calc(
                        p,
                        t0,
                        qv0,
                        ql0,
                        0.0,
                        a_dc,
                        ql_dc,
                        qi_dc,
                        a_sc,
                        ql_sc,
                        qi_sc,
                        ai0_st,
                        qlst_max,
                        tmin0,
                        tmax0,
                        landfrac[idx],
                        snowh[idx],
                        rhminl,
                        rhminl_adj_land,
                        rhminh,
                        cpair,
                        latvap,
                        qlst_min,
                        qlst_max,
                        camstfrac,
                        __ptr__(t),
                        __ptr__(qv),
                        __ptr__(ql),
                        __ptr__(qi),
                    )
                    idxmod = 1
                elif ql0_nc < qlst_min * al0_st:
                    _cldwat2m_instratus_core_calc(
                        p,
                        t0,
                        qv0,
                        ql0,
                        0.0,
                        a_dc,
                        ql_dc,
                        qi_dc,
                        a_sc,
                        ql_sc,
                        qi_sc,
                        ai0_st,
                        qlst_min,
                        tmin0,
                        tmax0,
                        landfrac[idx],
                        snowh[idx],
                        rhminl,
                        rhminl_adj_land,
                        rhminh,
                        cpair,
                        latvap,
                        qlst_min,
                        qlst_max,
                        camstfrac,
                        __ptr__(t),
                        __ptr__(qv),
                        __ptr__(ql),
                        __ptr__(qi),
                    )
                    idxmod = 1
                else:
                    status[0] = 1
                    return
            else:
                status[0] = 2
                return

        qi = qi0

        if idxmod == 1:
            ai_st_nc = _cldwat2m_aist_single_native(
                qv,
                t,
                p,
                qi,
                landfrac[idx],
                snowh[idx],
                rhmaxi,
                rhmini,
                rhminl,
                rhminl_adj_land,
                rhminh,
            )
            ai_st = (1.0 - a_dc - a_sc) * ai_st_nc
            _, qs, _ = _cldwat2m_qsat_water_native(t, p)
            u = qv / qs
            u_nc = u
            al_st_nc, g_nc = _cldwat2m_astg_single_native(
                camstfrac,
                u_nc,
                p,
                qv,
                landfrac[idx],
                snowh[idx],
                rhminl,
                rhminl_adj_land,
                rhminh,
            )
            al_st = (1.0 - a_dc - a_sc) * al_st_nc
        else:
            ai_st = (1.0 - a_dc - a_sc) * ai0_st_nc
            al_st = (1.0 - a_dc - a_sc) * al0_st_nc

        a_st = max(ai_st, al_st)

        if al_st == 0.0:
            ql_st = 0.0
        else:
            ql_st = ql / al_st
            ql_st = min(qlst_max, max(qlst_min, ql_st))
        if ai_st == 0.0:
            qi_st = 0.0
        else:
            qi_st = qi / ai_st

        qi = ai_st * qi_st
        ql = al_st * ql_st

        t = t0 - (latvap / cpair) * (ql0 - ql) - ((latvap + latice) / cpair) * (qi0 - qi)
        qv = qv0 + ql0 - ql + qi0 - qi

        t_out[idx] = t
        qv_out[idx] = qv
        ql_out[idx] = ql
        qi_out[idx] = qi
        al_st_out[idx] = al_st
        ai_st_out[idx] = ai_st
        ql_st_out[idx] = ql_st
        qi_st_out[idx] = qi_st


@export
def cldwat2m_instratus_tendency_codon(
    stage: int,
    ncol: int,
    dt: float,
    qsmall: float,
    cone: float,
    ql_new_p: cobj,
    qi_new_p: cobj,
    ql_old_p: cobj,
    qi_old_p: cobj,
    nl_old_p: cobj,
    ni_old_p: cobj,
    al_st_new_p: cobj,
    ai_st_new_p: cobj,
    a_st_new_p: cobj,
    qqw_p: cobj,
    qqi_p: cobj,
    qqnl_p: cobj,
    qqni_p: cobj,
    nl_new_p: cobj,
    ni_new_p: cobj,
):
    ql_new = Ptr[float](ql_new_p)
    qi_new = Ptr[float](qi_new_p)
    ql_old = Ptr[float](ql_old_p)
    qi_old = Ptr[float](qi_old_p)
    nl_old = Ptr[float](nl_old_p)
    ni_old = Ptr[float](ni_old_p)
    al_st_new = Ptr[float](al_st_new_p)
    ai_st_new = Ptr[float](ai_st_new_p)
    a_st_new = Ptr[float](a_st_new_p)
    qqw = Ptr[float](qqw_p)
    qqi = Ptr[float](qqi_p)
    qqnl = Ptr[float](qqnl_p)
    qqni = Ptr[float](qqni_p)
    nl_new = Ptr[float](nl_new_p)
    ni_new = Ptr[float](ni_new_p)

    for i in range(1, ncol + 1):
        idx = i - 1

        a_st_new[idx] = max(al_st_new[idx], ai_st_new[idx])
        qqw[idx] = (ql_new[idx] - ql_old[idx]) / dt
        qqi[idx] = (qi_new[idx] - qi_old[idx]) / dt
        qqnl[idx] = 0.0
        qqni[idx] = 0.0

        if qqw[idx] <= 0.0:
            ql_available = False
            if stage == 2:
                ql_available = ql_old[idx] >= qsmall
            else:
                ql_available = ql_old[idx] > qsmall
            if ql_available:
                val = qqw[idx] * nl_old[idx] / ql_old[idx]
                lower = -nl_old[idx] / dt
                if val < lower:
                    val = lower
                if val < 0.0:
                    qqnl[idx] = cone * val
                else:
                    qqnl[idx] = 0.0
            else:
                qqnl[idx] = 0.0

        if qqi[idx] <= 0.0:
            if qi_old[idx] > qsmall:
                val_i = qqi[idx] * ni_old[idx] / qi_old[idx]
                lower_i = -ni_old[idx] / dt
                if val_i < lower_i:
                    val_i = lower_i
                if val_i < 0.0:
                    qqni[idx] = cone * val_i
                else:
                    qqni[idx] = 0.0
            else:
                qqni[idx] = 0.0

        nl_val = nl_old[idx] + qqnl[idx] * dt
        if nl_val < 0.0:
            nl_val = 0.0
        nl_new[idx] = nl_val

        ni_val = ni_old[idx] + qqni[idx] * dt
        if ni_val < 0.0:
            ni_val = 0.0
        ni_new[idx] = ni_val


@export
def cldwat2m_dropnum_limit_codon(
    stage: int,
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    dt: float,
    qsmall: float,
    ql_p: cobj,
    qi_p: cobj,
    nl_p: cobj,
    ni_p: cobj,
    nlten_p: cobj,
    niten_p: cobj,
):
    ql = Ptr[float](ql_p)
    qi = Ptr[float](qi_p)
    nl = Ptr[float](nl_p)
    ni = Ptr[float](ni_p)
    nlten = Ptr[float](nlten_p)
    niten = Ptr[float](niten_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            nlten[idx] = 0.0
            niten[idx] = 0.0
            if ql[idx] < qsmall:
                nlten[idx] = -nl[idx] / dt
                nl[idx] = 0.0
            if qi[idx] < qsmall:
                niten[idx] = -ni[idx] / dt
                ni[idx] = 0.0


@export
def cldwat2m_ref_state_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    t_src_p: cobj,
    qv_src_p: cobj,
    ql_src_p: cobj,
    qi_src_p: cobj,
    al_st_src_p: cobj,
    ai_st_src_p: cobj,
    a_st_src_p: cobj,
    ql_st_src_p: cobj,
    qi_st_src_p: cobj,
    nl_src_p: cobj,
    ni_src_p: cobj,
    t_dst_p: cobj,
    qv_dst_p: cobj,
    ql_dst_p: cobj,
    qi_dst_p: cobj,
    al_st_dst_p: cobj,
    ai_st_dst_p: cobj,
    a_st_dst_p: cobj,
    ql_st_dst_p: cobj,
    qi_st_dst_p: cobj,
    nl_dst_p: cobj,
    ni_dst_p: cobj,
):
    t_src = Ptr[float](t_src_p)
    qv_src = Ptr[float](qv_src_p)
    ql_src = Ptr[float](ql_src_p)
    qi_src = Ptr[float](qi_src_p)
    al_st_src = Ptr[float](al_st_src_p)
    ai_st_src = Ptr[float](ai_st_src_p)
    a_st_src = Ptr[float](a_st_src_p)
    ql_st_src = Ptr[float](ql_st_src_p)
    qi_st_src = Ptr[float](qi_st_src_p)
    nl_src = Ptr[float](nl_src_p)
    ni_src = Ptr[float](ni_src_p)
    t_dst = Ptr[float](t_dst_p)
    qv_dst = Ptr[float](qv_dst_p)
    ql_dst = Ptr[float](ql_dst_p)
    qi_dst = Ptr[float](qi_dst_p)
    al_st_dst = Ptr[float](al_st_dst_p)
    ai_st_dst = Ptr[float](ai_st_dst_p)
    a_st_dst = Ptr[float](a_st_dst_p)
    ql_st_dst = Ptr[float](ql_st_dst_p)
    qi_st_dst = Ptr[float](qi_st_dst_p)
    nl_dst = Ptr[float](nl_dst_p)
    ni_dst = Ptr[float](ni_dst_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            t_dst[idx] = t_src[idx]
            qv_dst[idx] = qv_src[idx]
            ql_dst[idx] = ql_src[idx]
            qi_dst[idx] = qi_src[idx]
            al_st_dst[idx] = al_st_src[idx]
            ai_st_dst[idx] = ai_st_src[idx]
            a_st_dst[idx] = a_st_src[idx]
            ql_st_dst[idx] = ql_st_src[idx]
            qi_st_dst[idx] = qi_st_src[idx]
            nl_dst[idx] = nl_src[idx]
            ni_dst[idx] = ni_src[idx]


@export
def cldwat2m_input_state_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    dt: float,
    a_cu0_p: cobj,
    a_cud_p: cobj,
    t0_p: cobj,
    qv0_p: cobj,
    ql0_p: cobj,
    qi0_p: cobj,
    nl0_p: cobj,
    ni0_p: cobj,
    dacudt_p: cobj,
    t1_p: cobj,
    qv1_p: cobj,
    ql1_p: cobj,
    qi1_p: cobj,
    nl1_p: cobj,
    ni1_p: cobj,
):
    a_cu0 = Ptr[float](a_cu0_p)
    a_cud = Ptr[float](a_cud_p)
    t0 = Ptr[float](t0_p)
    qv0 = Ptr[float](qv0_p)
    ql0 = Ptr[float](ql0_p)
    qi0 = Ptr[float](qi0_p)
    nl0 = Ptr[float](nl0_p)
    ni0 = Ptr[float](ni0_p)
    dacudt = Ptr[float](dacudt_p)
    t1 = Ptr[float](t1_p)
    qv1 = Ptr[float](qv1_p)
    ql1 = Ptr[float](ql1_p)
    qi1 = Ptr[float](qi1_p)
    nl1 = Ptr[float](nl1_p)
    ni1 = Ptr[float](ni1_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            dacudt[idx] = (a_cu0[idx] - a_cud[idx]) / dt

    for k in range(1, top_lev):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            ql0[idx] = 0.0
            qi0[idx] = 0.0
            nl0[idx] = 0.0
            ni0[idx] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            t1[idx] = t0[idx]
            qv1[idx] = qv0[idx]
            ql1[idx] = ql0[idx]
            qi1[idx] = qi0[idx]
            nl1[idx] = nl0[idx]
            ni1[idx] = ni0[idx]


@export
def cldwat2m_qmin_fill_codon(
    ncol: int,
    pcols: int,
    pver: int,
    qvmin: float,
    qlmin: float,
    qimin: float,
    qmin1_p: cobj,
    qmin2_p: cobj,
    qmin3_p: cobj,
):
    qmin1 = Ptr[float](qmin1_p)
    qmin2 = Ptr[float](qmin2_p)
    qmin3 = Ptr[float](qmin3_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            qmin1[idx] = qvmin
            qmin2[idx] = qlmin
            qmin3[idx] = qimin


@export
def cldwat2m_detrain_state_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    dt: float,
    t_prime0_p: cobj,
    qv_prime0_p: cobj,
    ql_prime0_p: cobj,
    qi_prime0_p: cobj,
    nl_prime0_p: cobj,
    ni_prime0_p: cobj,
    d_t_p: cobj,
    d_qv_p: cobj,
    d_ql_p: cobj,
    d_qi_p: cobj,
    d_nl_p: cobj,
    d_ni_p: cobj,
    t_dprime_p: cobj,
    qv_dprime_p: cobj,
    ql_dprime_p: cobj,
    qi_dprime_p: cobj,
    nl_dprime_p: cobj,
    ni_dprime_p: cobj,
):
    t_prime0 = Ptr[float](t_prime0_p)
    qv_prime0 = Ptr[float](qv_prime0_p)
    ql_prime0 = Ptr[float](ql_prime0_p)
    qi_prime0 = Ptr[float](qi_prime0_p)
    nl_prime0 = Ptr[float](nl_prime0_p)
    ni_prime0 = Ptr[float](ni_prime0_p)
    d_t = Ptr[float](d_t_p)
    d_qv = Ptr[float](d_qv_p)
    d_ql = Ptr[float](d_ql_p)
    d_qi = Ptr[float](d_qi_p)
    d_nl = Ptr[float](d_nl_p)
    d_ni = Ptr[float](d_ni_p)
    t_dprime = Ptr[float](t_dprime_p)
    qv_dprime = Ptr[float](qv_dprime_p)
    ql_dprime = Ptr[float](ql_dprime_p)
    qi_dprime = Ptr[float](qi_dprime_p)
    nl_dprime = Ptr[float](nl_dprime_p)
    ni_dprime = Ptr[float](ni_dprime_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            t_dprime[idx] = t_prime0[idx]
            qv_dprime[idx] = qv_prime0[idx]
            ql_dprime[idx] = ql_prime0[idx]
            qi_dprime[idx] = qi_prime0[idx]
            nl_dprime[idx] = nl_prime0[idx]
            ni_dprime[idx] = ni_prime0[idx]

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            t_dprime[idx] = t_dprime[idx] + d_t[idx] * dt
            qv_dprime[idx] = qv_dprime[idx] + d_qv[idx] * dt
            ql_dprime[idx] = ql_dprime[idx] + d_ql[idx] * dt
            qi_dprime[idx] = qi_dprime[idx] + d_qi[idx] * dt
            nl_dprime[idx] = nl_dprime[idx] + d_nl[idx] * dt
            ni_dprime[idx] = ni_dprime[idx] + d_ni[idx] * dt


@export
def cldwat2m_qq_limiter_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    dt: float,
    qsmall: float,
    cone: float,
    qvmin: float,
    qv_05_p: cobj,
    ql_05_p: cobj,
    qi_05_p: cobj,
    nl_05_p: cobj,
    ni_05_p: cobj,
    qsat_a_p: cobj,
    qvwb_aw_p: cobj,
    qq_p: cobj,
    qqw_p: cobj,
    qqi_p: cobj,
    qqnl_p: cobj,
    qqni_p: cobj,
):
    qv_05 = Ptr[float](qv_05_p)
    ql_05 = Ptr[float](ql_05_p)
    qi_05 = Ptr[float](qi_05_p)
    nl_05 = Ptr[float](nl_05_p)
    ni_05 = Ptr[float](ni_05_p)
    qsat_a = Ptr[float](qsat_a_p)
    qvwb_aw = Ptr[float](qvwb_aw_p)
    qq = Ptr[float](qq_p)
    qqw = Ptr[float](qqw_p)
    qqi = Ptr[float](qqi_p)
    qqnl = Ptr[float](qqnl_p)
    qqni = Ptr[float](qqni_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            qqnl[idx] = 0.0
            qqni[idx] = 0.0

            if qq[idx] >= 0.0:
                qqmax = (qv_05[idx] - qvmin) / dt
                if qqmax < 0.0:
                    qqmax = 0.0
                if qq[idx] > qqmax:
                    qq[idx] = qqmax
                qqw[idx] = qq[idx]
                qqi[idx] = 0.0
            else:
                qqmin = 0.0
                if qv_05[idx] < qsat_a[idx]:
                    qqmin = cone * (qv_05[idx] - qvwb_aw[idx]) / dt
                    if qqmin > 0.0:
                        qqmin = 0.0
                if qq[idx] < qqmin:
                    qq[idx] = qqmin
                qqw[idx] = qq[idx]
                qqi[idx] = 0.0

                qqwmin = -cone * ql_05[idx] / dt
                if qqwmin > 0.0:
                    qqwmin = 0.0
                qqimin = -cone * qi_05[idx] / dt
                if qqimin > 0.0:
                    qqimin = 0.0

                qqw_val = qqw[idx]
                if qqw_val < qqwmin:
                    qqw_val = qqwmin
                if qqw_val < 0.0:
                    qqw[idx] = qqw_val
                else:
                    qqw[idx] = 0.0

                qqi_val = qqi[idx]
                if qqi_val < qqimin:
                    qqi_val = qqimin
                if qqi_val < 0.0:
                    qqi[idx] = qqi_val
                else:
                    qqi[idx] = 0.0

            if qqw[idx] < 0.0:
                if ql_05[idx] > qsmall:
                    qqnl_val = qqw[idx] * nl_05[idx] / ql_05[idx]
                    lower = -nl_05[idx] / dt
                    if qqnl_val < lower:
                        qqnl_val = lower
                    qqnl_val = cone * qqnl_val
                    if qqnl_val < 0.0:
                        qqnl[idx] = qqnl_val
                    else:
                        qqnl[idx] = 0.0
                else:
                    qqnl[idx] = 0.0

            if qqi[idx] < 0.0:
                if qi_05[idx] > qsmall:
                    qqni_val = qqi[idx] * ni_05[idx] / qi_05[idx]
                    lower_i = -ni_05[idx] / dt
                    if qqni_val < lower_i:
                        qqni_val = lower_i
                    qqni_val = cone * qqni_val
                    if qqni_val < 0.0:
                        qqni[idx] = qqni_val
                    else:
                        qqni[idx] = 0.0
                else:
                    qqni[idx] = 0.0


@export
def cldwat2m_iter_zero_codon(
    pcols: int,
    pver: int,
    qq_p: cobj,
    qqw_p: cobj,
    qqi_p: cobj,
    qqnl_p: cobj,
    qqni_p: cobj,
    qqw2_p: cobj,
    qqi2_p: cobj,
    qqnl2_p: cobj,
    qqni2_p: cobj,
    nlten_pwi2_p: cobj,
    niten_pwi2_p: cobj,
    acnl_p: cobj,
    acni_p: cobj,
    aa_p: cobj,
    bb_p: cobj,
):
    qq = Ptr[float](qq_p)
    qqw = Ptr[float](qqw_p)
    qqi = Ptr[float](qqi_p)
    qqnl = Ptr[float](qqnl_p)
    qqni = Ptr[float](qqni_p)
    qqw2 = Ptr[float](qqw2_p)
    qqi2 = Ptr[float](qqi2_p)
    qqnl2 = Ptr[float](qqnl2_p)
    qqni2 = Ptr[float](qqni2_p)
    nlten_pwi2 = Ptr[float](nlten_pwi2_p)
    niten_pwi2 = Ptr[float](niten_pwi2_p)
    acnl = Ptr[float](acnl_p)
    acni = Ptr[float](acni_p)
    aa = Ptr[float](aa_p)
    bb = Ptr[float](bb_p)

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            idx = _idx2(i, k, pcols)
            qq[idx] = 0.0
            qqw[idx] = 0.0
            qqi[idx] = 0.0
            qqnl[idx] = 0.0
            qqni[idx] = 0.0
            qqw2[idx] = 0.0
            qqi2[idx] = 0.0
            qqnl2[idx] = 0.0
            qqni2[idx] = 0.0
            nlten_pwi2[idx] = 0.0
            niten_pwi2[idx] = 0.0
            acnl[idx] = 0.0
            acni[idx] = 0.0

    for idx_small in range(0, 4):
        aa[idx_small] = 0.0
    for idx_small in range(0, 2):
        bb[idx_small] = 0.0


@export
def cldwat2m_zero16_ncol_codon(
    ncol: int,
    pcols: int,
    pver: int,
    arr01_p: cobj,
    arr02_p: cobj,
    arr03_p: cobj,
    arr04_p: cobj,
    arr05_p: cobj,
    arr06_p: cobj,
    arr07_p: cobj,
    arr08_p: cobj,
    arr09_p: cobj,
    arr10_p: cobj,
    arr11_p: cobj,
    arr12_p: cobj,
    arr13_p: cobj,
    arr14_p: cobj,
    arr15_p: cobj,
    arr16_p: cobj,
):
    arr01 = Ptr[float](arr01_p)
    arr02 = Ptr[float](arr02_p)
    arr03 = Ptr[float](arr03_p)
    arr04 = Ptr[float](arr04_p)
    arr05 = Ptr[float](arr05_p)
    arr06 = Ptr[float](arr06_p)
    arr07 = Ptr[float](arr07_p)
    arr08 = Ptr[float](arr08_p)
    arr09 = Ptr[float](arr09_p)
    arr10 = Ptr[float](arr10_p)
    arr11 = Ptr[float](arr11_p)
    arr12 = Ptr[float](arr12_p)
    arr13 = Ptr[float](arr13_p)
    arr14 = Ptr[float](arr14_p)
    arr15 = Ptr[float](arr15_p)
    arr16 = Ptr[float](arr16_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            arr01[idx] = 0.0
            arr02[idx] = 0.0
            arr03[idx] = 0.0
            arr04[idx] = 0.0
            arr05[idx] = 0.0
            arr06[idx] = 0.0
            arr07[idx] = 0.0
            arr08[idx] = 0.0
            arr09[idx] = 0.0
            arr10[idx] = 0.0
            arr11[idx] = 0.0
            arr12[idx] = 0.0
            arr13[idx] = 0.0
            arr14[idx] = 0.0
            arr15[idx] = 0.0
            arr16[idx] = 0.0


@export
def cldwat2m_iter_column_state_codon(
    iter: int,
    ncol: int,
    qsat_b_p: cobj,
    qv_p: cobj,
    a_cud_p: cobj,
    a_cu0_p: cobj,
    a_cu_p: cobj,
    u_p: cobj,
    u_nc_p: cobj,
):
    qsat_b = Ptr[float](qsat_b_p)
    qv = Ptr[float](qv_p)
    a_cud = Ptr[float](a_cud_p)
    a_cu0 = Ptr[float](a_cu0_p)
    a_cu = Ptr[float](a_cu_p)
    u = Ptr[float](u_p)
    u_nc = Ptr[float](u_nc_p)

    if iter == 1:
        for idx in range(ncol):
            a_cu[idx] = a_cud[idx]
    else:
        for idx in range(ncol):
            a_cu[idx] = a_cu0[idx]

    for idx in range(ncol):
        u[idx] = qv[idx] / qsat_b[idx]
        u_nc[idx] = u[idx]


@export
def cldwat2m_iter_column_stratus_codon(
    ncol: int,
    a_cu_p: cobj,
    al_st_nc_p: cobj,
    ai_st_nc_p: cobj,
    al_st_p: cobj,
    ai_st_p: cobj,
    a_st_p: cobj,
):
    a_cu = Ptr[float](a_cu_p)
    al_st_nc = Ptr[float](al_st_nc_p)
    ai_st_nc = Ptr[float](ai_st_nc_p)
    al_st = Ptr[float](al_st_p)
    ai_st = Ptr[float](ai_st_p)
    a_st = Ptr[float](a_st_p)

    for idx in range(ncol):
        ai_st[idx] = (1.0 - a_cu[idx]) * ai_st_nc[idx]
    for idx in range(ncol):
        al_st[idx] = (1.0 - a_cu[idx]) * al_st_nc[idx]
    for idx in range(ncol):
        a_st[idx] = max(al_st[idx], ai_st[idx])


@export
def cldwat2m_qq_coeff_solve_codon(
    k: int,
    ncol: int,
    pcols: int,
    pver: int,
    latvap: float,
    latice: float,
    cpair: float,
    cc: float,
    qsat_b_p: cobj,
    dqsdT_b_p: cobj,
    qv_p: cobj,
    a_t_p: cobj,
    a_t_adj_p: cobj,
    a_ql_p: cobj,
    a_ql_adj_p: cobj,
    a_qi_p: cobj,
    a_qi_adj_p: cobj,
    a_qv_p: cobj,
    a_qv_adj_p: cobj,
    c_t_p: cobj,
    c_ql_p: cobj,
    c_qi_p: cobj,
    c_qv_p: cobj,
    c_qlst_p: cobj,
    a_cu_p: cobj,
    g_nc_p: cobj,
    al_st_p: cobj,
    ql_st_p: cobj,
    al_st_nc_p: cobj,
    dacudt_p: cobj,
    f_nc_p: cobj,
    qq_p: cobj,
) -> int:
    qsat_b = Ptr[float](qsat_b_p)
    dqsdT_b = Ptr[float](dqsdT_b_p)
    qv = Ptr[float](qv_p)
    a_t = Ptr[float](a_t_p)
    a_t_adj = Ptr[float](a_t_adj_p)
    a_ql = Ptr[float](a_ql_p)
    a_ql_adj = Ptr[float](a_ql_adj_p)
    a_qi = Ptr[float](a_qi_p)
    a_qi_adj = Ptr[float](a_qi_adj_p)
    a_qv = Ptr[float](a_qv_p)
    a_qv_adj = Ptr[float](a_qv_adj_p)
    c_t = Ptr[float](c_t_p)
    c_ql = Ptr[float](c_ql_p)
    c_qi = Ptr[float](c_qi_p)
    c_qv = Ptr[float](c_qv_p)
    c_qlst = Ptr[float](c_qlst_p)
    a_cu = Ptr[float](a_cu_p)
    g_nc = Ptr[float](g_nc_p)
    al_st = Ptr[float](al_st_p)
    ql_st = Ptr[float](ql_st_p)
    al_st_nc = Ptr[float](al_st_nc_p)
    dacudt = Ptr[float](dacudt_p)
    f_nc = Ptr[float](f_nc_p)
    qq = Ptr[float](qq_p)

    for i in range(1, ncol + 1):
        idx = _idx2(i, k, pcols)
        idx1 = i - 1

        alpha = 1.0 / qsat_b[idx1]
        beta = dqsdT_b[idx1] * (qv[idx] / (qsat_b[idx1] ** 2))
        betast = alpha * dqsdT_b[idx1]
        gammal = alpha + (latvap / cpair) * beta
        gammai = alpha + ((latvap + latice) / cpair) * beta
        gammaQ = alpha + (latvap / cpair) * beta
        deltal = 1.0 + al_st[idx] * (latvap / cpair) * (betast / alpha)
        deltai = 1.0 + al_st[idx] * ((latvap + latice) / cpair) * (betast / alpha)
        a_tc = a_t[idx] + a_t_adj[idx] - (latvap / cpair) * (a_ql[idx] + a_ql_adj[idx]) - (
            (latvap + latice) / cpair
        ) * (a_qi[idx] + a_qi_adj[idx])
        a_qt = a_qv[idx] + a_qv_adj[idx] + a_ql[idx] + a_ql_adj[idx] + a_qi[idx] + a_qi_adj[idx]
        c_tc = c_t[idx] - (latvap / cpair) * c_ql[idx] - ((latvap + latice) / cpair) * c_qi[idx]
        c_qt = c_qv[idx] + c_ql[idx] + c_qi[idx]
        dTcdt = a_tc + c_tc
        dqtdt = a_qt + c_qt
        dqtstldt = a_qt - a_qi[idx] - a_qi_adj[idx] + c_qlst[idx]
        dqidt = a_qi[idx] + a_qi_adj[idx] + c_qi[idx]

        anic = max(1.0e-8, (1.0 - a_cu[idx]))
        gg = g_nc[idx] / anic
        a11 = gammal * al_st[idx]
        a12 = gg + gammal * cc * ql_st[idx]
        a21 = alpha + (latvap / cpair) * betast * al_st[idx]
        a22 = (latvap / cpair) * betast * cc * ql_st[idx]
        b1 = alpha * dqtdt - beta * dTcdt - gammai * dqidt - gg * al_st_nc[idx] * dacudt[idx] + f_nc[idx]
        b2 = alpha * dqtstldt - betast * (dTcdt + ((latvap + latice) / cpair) * dqidt)

        ipiv1 = 0
        ipiv2 = 0
        for iter_idx in range(1, 3):
            big = 0.0
            irow = 1
            icol = 1

            if ipiv1 != 1:
                if ipiv1 == 0:
                    value = abs(a11)
                    if value >= big:
                        big = value
                        irow = 1
                        icol = 1
                elif ipiv1 > 1:
                    return 1
                if ipiv2 == 0:
                    value = abs(a12)
                    if value >= big:
                        big = value
                        irow = 1
                        icol = 2
                elif ipiv2 > 1:
                    return 1

            if ipiv2 != 1:
                if ipiv1 == 0:
                    value = abs(a21)
                    if value >= big:
                        big = value
                        irow = 2
                        icol = 1
                elif ipiv1 > 1:
                    return 1
                if ipiv2 == 0:
                    value = abs(a22)
                    if value >= big:
                        big = value
                        irow = 2
                        icol = 2
                elif ipiv2 > 1:
                    return 1

            if icol == 1:
                ipiv1 += 1
            else:
                ipiv2 += 1

            if irow != icol:
                dum = a11
                a11 = a21
                a21 = dum
                dum = a12
                a12 = a22
                a22 = dum
                dum = b1
                b1 = b2
                b2 = dum

            if icol == 1:
                if a11 == 0.0:
                    return 2
                pivinv = 1.0 / a11
                a11 = 1.0
                a11 = a11 * pivinv
                a12 = a12 * pivinv
                b1 = b1 * pivinv
                dum = a21
                a21 = 0.0
                a21 = a21 - a11 * dum
                a22 = a22 - a12 * dum
                b2 = b2 - b1 * dum
            else:
                if a22 == 0.0:
                    return 2
                pivinv = 1.0 / a22
                a22 = 1.0
                a21 = a21 * pivinv
                a22 = a22 * pivinv
                b2 = b2 * pivinv
                dum = a12
                a12 = 0.0
                a11 = a11 - a21 * dum
                a12 = a12 - a22 * dum
                b1 = b1 - b2 * dum

        dqlstdt = b1
        dalstdt = b2
        qq[idx] = al_st[idx] * dqlstdt + cc * ql_st[idx] * dalstdt - (a_ql[idx] + a_ql_adj[idx] + c_ql[idx])

        gammaQ = gammaQ
        deltal = deltal
        deltai = deltai

    return 0


@export
def cldwat2m_advective_state_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    dt: float,
    t_0_p: cobj,
    qv_0_p: cobj,
    ql_0_p: cobj,
    qi_0_p: cobj,
    nl_0_p: cobj,
    ni_0_p: cobj,
    a_t_p: cobj,
    c_t_p: cobj,
    a_qv_p: cobj,
    c_qv_p: cobj,
    a_ql_p: cobj,
    c_ql_p: cobj,
    a_qi_p: cobj,
    c_qi_p: cobj,
    a_nl_p: cobj,
    c_nl_p: cobj,
    a_ni_p: cobj,
    c_ni_p: cobj,
    t_05_p: cobj,
    qv_05_p: cobj,
    ql_05_p: cobj,
    qi_05_p: cobj,
    nl_05_p: cobj,
    ni_05_p: cobj,
):
    t_0 = Ptr[float](t_0_p)
    qv_0 = Ptr[float](qv_0_p)
    ql_0 = Ptr[float](ql_0_p)
    qi_0 = Ptr[float](qi_0_p)
    nl_0 = Ptr[float](nl_0_p)
    ni_0 = Ptr[float](ni_0_p)
    a_t = Ptr[float](a_t_p)
    c_t = Ptr[float](c_t_p)
    a_qv = Ptr[float](a_qv_p)
    c_qv = Ptr[float](c_qv_p)
    a_ql = Ptr[float](a_ql_p)
    c_ql = Ptr[float](c_ql_p)
    a_qi = Ptr[float](a_qi_p)
    c_qi = Ptr[float](c_qi_p)
    a_nl = Ptr[float](a_nl_p)
    c_nl = Ptr[float](c_nl_p)
    a_ni = Ptr[float](a_ni_p)
    c_ni = Ptr[float](c_ni_p)
    t_05 = Ptr[float](t_05_p)
    qv_05 = Ptr[float](qv_05_p)
    ql_05 = Ptr[float](ql_05_p)
    qi_05 = Ptr[float](qi_05_p)
    nl_05 = Ptr[float](nl_05_p)
    ni_05 = Ptr[float](ni_05_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            t_05[idx] = t_0[idx] + (a_t[idx] + c_t[idx]) * dt
            qv_05[idx] = qv_0[idx] + (a_qv[idx] + c_qv[idx]) * dt
            ql_05[idx] = ql_0[idx] + (a_ql[idx] + c_ql[idx]) * dt
            qi_05[idx] = qi_0[idx] + (a_qi[idx] + c_qi[idx]) * dt

            nl_val = nl_0[idx] + (a_nl[idx] + c_nl[idx]) * dt
            if nl_val < 0.0:
                nl_val = 0.0
            nl_05[idx] = nl_val

            ni_val = ni_0[idx] + (a_ni[idx] + c_ni[idx]) * dt
            if ni_val < 0.0:
                ni_val = 0.0
            ni_05[idx] = ni_val


@export
def cldwat2m_iter_state_codon(
    iter_num: int,
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    dt: float,
    ramda: float,
    qsmall: float,
    latvap: float,
    latice: float,
    cpair: float,
    qqw_p: cobj,
    qqi_p: cobj,
    qqnl_p: cobj,
    qqni_p: cobj,
    qqw_prev_p: cobj,
    qqi_prev_p: cobj,
    qqnl_prev_p: cobj,
    qqni_prev_p: cobj,
    qqw_prog_p: cobj,
    qqi_prog_p: cobj,
    qqnl_prog_p: cobj,
    qqni_prog_p: cobj,
    t_0_p: cobj,
    qv_0_p: cobj,
    ql_0_p: cobj,
    qi_0_p: cobj,
    nl_0_p: cobj,
    ni_0_p: cobj,
    a_t_p: cobj,
    a_t_adj_p: cobj,
    c_t_p: cobj,
    a_qv_p: cobj,
    a_qv_adj_p: cobj,
    c_qv_p: cobj,
    a_ql_p: cobj,
    a_ql_adj_p: cobj,
    c_ql_p: cobj,
    a_qi_p: cobj,
    a_qi_adj_p: cobj,
    c_qi_p: cobj,
    a_nl_p: cobj,
    c_nl_p: cobj,
    a_ni_p: cobj,
    c_ni_p: cobj,
    t_prime0_p: cobj,
    qv_prime0_p: cobj,
    ql_prime0_p: cobj,
    qi_prime0_p: cobj,
    nl_prime0_p: cobj,
    ni_prime0_p: cobj,
):
    qqw = Ptr[float](qqw_p)
    qqi = Ptr[float](qqi_p)
    qqnl = Ptr[float](qqnl_p)
    qqni = Ptr[float](qqni_p)
    qqw_prev = Ptr[float](qqw_prev_p)
    qqi_prev = Ptr[float](qqi_prev_p)
    qqnl_prev = Ptr[float](qqnl_prev_p)
    qqni_prev = Ptr[float](qqni_prev_p)
    qqw_prog = Ptr[float](qqw_prog_p)
    qqi_prog = Ptr[float](qqi_prog_p)
    qqnl_prog = Ptr[float](qqnl_prog_p)
    qqni_prog = Ptr[float](qqni_prog_p)
    t_0 = Ptr[float](t_0_p)
    qv_0 = Ptr[float](qv_0_p)
    ql_0 = Ptr[float](ql_0_p)
    qi_0 = Ptr[float](qi_0_p)
    nl_0 = Ptr[float](nl_0_p)
    ni_0 = Ptr[float](ni_0_p)
    a_t = Ptr[float](a_t_p)
    a_t_adj = Ptr[float](a_t_adj_p)
    c_t = Ptr[float](c_t_p)
    a_qv = Ptr[float](a_qv_p)
    a_qv_adj = Ptr[float](a_qv_adj_p)
    c_qv = Ptr[float](c_qv_p)
    a_ql = Ptr[float](a_ql_p)
    a_ql_adj = Ptr[float](a_ql_adj_p)
    c_ql = Ptr[float](c_ql_p)
    a_qi = Ptr[float](a_qi_p)
    a_qi_adj = Ptr[float](a_qi_adj_p)
    c_qi = Ptr[float](c_qi_p)
    a_nl = Ptr[float](a_nl_p)
    c_nl = Ptr[float](c_nl_p)
    a_ni = Ptr[float](a_ni_p)
    c_ni = Ptr[float](c_ni_p)
    t_prime0 = Ptr[float](t_prime0_p)
    qv_prime0 = Ptr[float](qv_prime0_p)
    ql_prime0 = Ptr[float](ql_prime0_p)
    qi_prime0 = Ptr[float](qi_prime0_p)
    nl_prime0 = Ptr[float](nl_prime0_p)
    ni_prime0 = Ptr[float](ni_prime0_p)

    if iter_num == 1:
        for k in range(top_lev, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, pcols)
                qqw_prev[idx] = qqw[idx]
                qqi_prev[idx] = qqi[idx]
                qqnl_prev[idx] = qqnl[idx]
                qqni_prev[idx] = qqni[idx]

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)

            qqw_prog[idx] = ramda * qqw[idx] + (1.0 - ramda) * qqw_prev[idx]
            qqi_prog[idx] = ramda * qqi[idx] + (1.0 - ramda) * qqi_prev[idx]
            qqnl_prog[idx] = ramda * qqnl[idx] + (1.0 - ramda) * qqnl_prev[idx]
            qqni_prog[idx] = ramda * qqni[idx] + (1.0 - ramda) * qqni_prev[idx]

            qqw_prev[idx] = qqw_prog[idx]
            qqi_prev[idx] = qqi_prog[idx]
            qqnl_prev[idx] = qqnl_prog[idx]
            qqni_prev[idx] = qqni_prog[idx]

            latent_term = (latvap * qqw_prog[idx] + (latvap + latice) * qqi_prog[idx]) / cpair
            t_prime0[idx] = t_0[idx] + dt * (a_t[idx] + a_t_adj[idx] + c_t[idx] + latent_term)
            qv_prime0[idx] = qv_0[idx] + dt * (
                a_qv[idx] + a_qv_adj[idx] + c_qv[idx] - qqw_prog[idx] - qqi_prog[idx]
            )
            ql_prime0[idx] = ql_0[idx] + dt * (a_ql[idx] + a_ql_adj[idx] + c_ql[idx] + qqw_prog[idx])
            qi_prime0[idx] = qi_0[idx] + dt * (a_qi[idx] + a_qi_adj[idx] + c_qi[idx] + qqi_prog[idx])

            nl_val = nl_0[idx] + dt * (a_nl[idx] + c_nl[idx] + qqnl_prog[idx])
            if nl_val < 0.0:
                nl_val = 0.0
            nl_prime0[idx] = nl_val

            ni_val = ni_0[idx] + dt * (a_ni[idx] + c_ni[idx] + qqni_prog[idx])
            if ni_val < 0.0:
                ni_val = 0.0
            ni_prime0[idx] = ni_val

            if ql_prime0[idx] < qsmall:
                nl_prime0[idx] = 0.0
            if qi_prime0[idx] < qsmall:
                ni_prime0[idx] = 0.0


@export
def cldwat2m_final_tendency_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    do_cldice: int,
    dt: float,
    cpair: float,
    qqw_prog_p: cobj,
    qqi_prog_p: cobj,
    qqnl_prog_p: cobj,
    qqni_prog_p: cobj,
    qqw1_p: cobj,
    qqi1_p: cobj,
    qqw2_p: cobj,
    qqi2_p: cobj,
    qlten_pwi1_p: cobj,
    qlten_pwi2_p: cobj,
    qiten_pwi1_p: cobj,
    qiten_pwi2_p: cobj,
    a_ql_adj_p: cobj,
    a_qi_adj_p: cobj,
    qqnl1_p: cobj,
    qqni1_p: cobj,
    qqnl2_p: cobj,
    qqni2_p: cobj,
    nlten_pwi1_p: cobj,
    nlten_pwi2_p: cobj,
    niten_pwi1_p: cobj,
    niten_pwi2_p: cobj,
    acnl_p: cobj,
    acni_p: cobj,
    a_nl_adj_p: cobj,
    a_ni_adj_p: cobj,
    qvten_pwi1_p: cobj,
    qvten_pwi2_p: cobj,
    a_qv_adj_p: cobj,
    t_star_p: cobj,
    qv_star_p: cobj,
    ql_star_p: cobj,
    qi_star_p: cobj,
    nl_star_p: cobj,
    ni_star_p: cobj,
    a_t_p: cobj,
    c_t_p: cobj,
    a_qv_p: cobj,
    c_qv_p: cobj,
    a_ql_p: cobj,
    c_ql_p: cobj,
    a_qi_p: cobj,
    c_qi_p: cobj,
    a_nl_p: cobj,
    c_nl_p: cobj,
    a_ni_p: cobj,
    c_ni_p: cobj,
    a_st_star_p: cobj,
    a_cu0_p: cobj,
    qqw_final_p: cobj,
    qqi_final_p: cobj,
    qq_final_p: cobj,
    qqw_all_p: cobj,
    qqi_all_p: cobj,
    qq_all_p: cobj,
    qqnl_final_p: cobj,
    qqni_final_p: cobj,
    qqn_final_p: cobj,
    qqnl_all_p: cobj,
    qqni_all_p: cobj,
    qqn_all_p: cobj,
    qme_p: cobj,
    qvadj_p: cobj,
    qladj_p: cobj,
    qiadj_p: cobj,
    qllim_p: cobj,
    qilim_p: cobj,
    s_tendout_p: cobj,
    qv_tendout_p: cobj,
    ql_tendout_p: cobj,
    qi_tendout_p: cobj,
    nl_tendout_p: cobj,
    ni_tendout_p: cobj,
    cld_p: cobj,
    t0_p: cobj,
    qv0_p: cobj,
    ql0_p: cobj,
    qi0_p: cobj,
    nl0_p: cobj,
    ni0_p: cobj,
):
    qqw_prog = Ptr[float](qqw_prog_p)
    qqi_prog = Ptr[float](qqi_prog_p)
    qqnl_prog = Ptr[float](qqnl_prog_p)
    qqni_prog = Ptr[float](qqni_prog_p)
    qqw1 = Ptr[float](qqw1_p)
    qqi1 = Ptr[float](qqi1_p)
    qqw2 = Ptr[float](qqw2_p)
    qqi2 = Ptr[float](qqi2_p)
    qlten_pwi1 = Ptr[float](qlten_pwi1_p)
    qlten_pwi2 = Ptr[float](qlten_pwi2_p)
    qiten_pwi1 = Ptr[float](qiten_pwi1_p)
    qiten_pwi2 = Ptr[float](qiten_pwi2_p)
    a_ql_adj = Ptr[float](a_ql_adj_p)
    a_qi_adj = Ptr[float](a_qi_adj_p)
    qqnl1 = Ptr[float](qqnl1_p)
    qqni1 = Ptr[float](qqni1_p)
    qqnl2 = Ptr[float](qqnl2_p)
    qqni2 = Ptr[float](qqni2_p)
    nlten_pwi1 = Ptr[float](nlten_pwi1_p)
    nlten_pwi2 = Ptr[float](nlten_pwi2_p)
    niten_pwi1 = Ptr[float](niten_pwi1_p)
    niten_pwi2 = Ptr[float](niten_pwi2_p)
    acnl = Ptr[float](acnl_p)
    acni = Ptr[float](acni_p)
    a_nl_adj = Ptr[float](a_nl_adj_p)
    a_ni_adj = Ptr[float](a_ni_adj_p)
    qvten_pwi1 = Ptr[float](qvten_pwi1_p)
    qvten_pwi2 = Ptr[float](qvten_pwi2_p)
    a_qv_adj = Ptr[float](a_qv_adj_p)
    t_star = Ptr[float](t_star_p)
    qv_star = Ptr[float](qv_star_p)
    ql_star = Ptr[float](ql_star_p)
    qi_star = Ptr[float](qi_star_p)
    nl_star = Ptr[float](nl_star_p)
    ni_star = Ptr[float](ni_star_p)
    a_t = Ptr[float](a_t_p)
    c_t = Ptr[float](c_t_p)
    a_qv = Ptr[float](a_qv_p)
    c_qv = Ptr[float](c_qv_p)
    a_ql = Ptr[float](a_ql_p)
    c_ql = Ptr[float](c_ql_p)
    a_qi = Ptr[float](a_qi_p)
    c_qi = Ptr[float](c_qi_p)
    a_nl = Ptr[float](a_nl_p)
    c_nl = Ptr[float](c_nl_p)
    a_ni = Ptr[float](a_ni_p)
    c_ni = Ptr[float](c_ni_p)
    a_st_star = Ptr[float](a_st_star_p)
    a_cu0 = Ptr[float](a_cu0_p)
    qqw_final = Ptr[float](qqw_final_p)
    qqi_final = Ptr[float](qqi_final_p)
    qq_final = Ptr[float](qq_final_p)
    qqw_all = Ptr[float](qqw_all_p)
    qqi_all = Ptr[float](qqi_all_p)
    qq_all = Ptr[float](qq_all_p)
    qqnl_final = Ptr[float](qqnl_final_p)
    qqni_final = Ptr[float](qqni_final_p)
    qqn_final = Ptr[float](qqn_final_p)
    qqnl_all = Ptr[float](qqnl_all_p)
    qqni_all = Ptr[float](qqni_all_p)
    qqn_all = Ptr[float](qqn_all_p)
    qme = Ptr[float](qme_p)
    qvadj = Ptr[float](qvadj_p)
    qladj = Ptr[float](qladj_p)
    qiadj = Ptr[float](qiadj_p)
    qllim = Ptr[float](qllim_p)
    qilim = Ptr[float](qilim_p)
    s_tendout = Ptr[float](s_tendout_p)
    qv_tendout = Ptr[float](qv_tendout_p)
    ql_tendout = Ptr[float](ql_tendout_p)
    qi_tendout = Ptr[float](qi_tendout_p)
    nl_tendout = Ptr[float](nl_tendout_p)
    ni_tendout = Ptr[float](ni_tendout_p)
    cld = Ptr[float](cld_p)
    t0 = Ptr[float](t0_p)
    qv0 = Ptr[float](qv0_p)
    ql0 = Ptr[float](ql0_p)
    qi0 = Ptr[float](qi0_p)
    nl0 = Ptr[float](nl0_p)
    ni0 = Ptr[float](ni0_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)

            qqw_final[idx] = qqw_prog[idx]
            qqi_final[idx] = qqi_prog[idx]
            qq_final[idx] = qqw_final[idx] + qqi_final[idx]

            qqw_all_val = qqw_prog[idx]
            qqw_all_val = qqw_all_val + qqw1[idx]
            qqw_all_val = qqw_all_val + qqw2[idx]
            qqw_all_val = qqw_all_val + qlten_pwi1[idx]
            qqw_all_val = qqw_all_val + qlten_pwi2[idx]
            qqw_all_val = qqw_all_val + a_ql_adj[idx]
            qqw_all[idx] = qqw_all_val

            qqi_all_val = qqi_prog[idx]
            qqi_all_val = qqi_all_val + qqi1[idx]
            qqi_all_val = qqi_all_val + qqi2[idx]
            qqi_all_val = qqi_all_val + qiten_pwi1[idx]
            qqi_all_val = qqi_all_val + qiten_pwi2[idx]
            qqi_all_val = qqi_all_val + a_qi_adj[idx]
            qqi_all[idx] = qqi_all_val

            qq_all[idx] = qqw_all[idx] + qqi_all[idx]

            qqnl_final[idx] = qqnl_prog[idx]
            qqni_final[idx] = qqni_prog[idx]
            qqn_final[idx] = qqnl_final[idx] + qqni_final[idx]

            qqnl_all_val = qqnl_prog[idx]
            qqnl_all_val = qqnl_all_val + qqnl1[idx]
            qqnl_all_val = qqnl_all_val + qqnl2[idx]
            qqnl_all_val = qqnl_all_val + nlten_pwi1[idx]
            qqnl_all_val = qqnl_all_val + nlten_pwi2[idx]
            qqnl_all_val = qqnl_all_val + acnl[idx]
            qqnl_all_val = qqnl_all_val + a_nl_adj[idx]
            qqnl_all[idx] = qqnl_all_val

            qqni_all_val = qqni_prog[idx]
            qqni_all_val = qqni_all_val + qqni1[idx]
            qqni_all_val = qqni_all_val + qqni2[idx]
            qqni_all_val = qqni_all_val + niten_pwi1[idx]
            qqni_all_val = qqni_all_val + niten_pwi2[idx]
            qqni_all_val = qqni_all_val + acni[idx]
            qqni_all_val = qqni_all_val + a_ni_adj[idx]
            qqni_all[idx] = qqni_all_val

            qqn_all[idx] = qqnl_all[idx] + qqni_all[idx]
            qme[idx] = qq_final[idx]

            qvadj_val = qvten_pwi1[idx]
            qvadj_val = qvadj_val + qvten_pwi2[idx]
            qvadj_val = qvadj_val + a_qv_adj[idx]
            qvadj[idx] = qvadj_val

            qladj_val = qlten_pwi1[idx]
            qladj_val = qladj_val + qlten_pwi2[idx]
            qladj_val = qladj_val + a_ql_adj[idx]
            qladj[idx] = qladj_val

            qiadj_val = qiten_pwi1[idx]
            qiadj_val = qiadj_val + qiten_pwi2[idx]
            qiadj_val = qiadj_val + a_qi_adj[idx]
            qiadj[idx] = qiadj_val

            qllim[idx] = qqw1[idx] + qqw2[idx]
            qilim[idx] = qqi1[idx] + qqi2[idx]

            s_tendout[idx] = cpair * (t_star[idx] - t0[idx]) / dt - cpair * (a_t[idx] + c_t[idx])
            qv_tendout[idx] = (qv_star[idx] - qv0[idx]) / dt - (a_qv[idx] + c_qv[idx])
            ql_tendout[idx] = (ql_star[idx] - ql0[idx]) / dt - (a_ql[idx] + c_ql[idx])
            qi_tendout[idx] = (qi_star[idx] - qi0[idx]) / dt - (a_qi[idx] + c_qi[idx])
            nl_tendout[idx] = (nl_star[idx] - nl0[idx]) / dt - (a_nl[idx] + c_nl[idx])
            ni_tendout[idx] = (ni_star[idx] - ni0[idx]) / dt - (a_ni[idx] + c_ni[idx])

            if do_cldice == 0:
                qi_tendout[idx] = 0.0
                ni_tendout[idx] = 0.0

            cld[idx] = a_st_star[idx] + a_cu0[idx]

            t0[idx] = t_star[idx]
            qv0[idx] = qv_star[idx]
            ql0[idx] = ql_star[idx]
            qi0[idx] = qi_star[idx]
            nl0[idx] = nl_star[idx]
            ni0[idx] = ni_star[idx]


@export
def cldwat2m_rhcrit_const_codon(
    pcols: int,
    pver: int,
    rhmini_const: float,
    rhminl_const: float,
    rhminl_adj_land_const: float,
    rhminh_const: float,
    rhmini_p: cobj,
    rhminl_p: cobj,
    rhminl_adj_land_p: cobj,
    rhminh_p: cobj,
):
    rhmini = Ptr[float](rhmini_p)
    rhminl = Ptr[float](rhminl_p)
    rhminl_adj_land = Ptr[float](rhminl_adj_land_p)
    rhminh = Ptr[float](rhminh_p)

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            idx = _idx2(i, k, pcols)
            rhmini[idx] = rhmini_const
            rhminl[idx] = rhminl_const
            rhminl_adj_land[idx] = rhminl_adj_land_const
            rhminh[idx] = rhminh_const


@export
def cldwat2m_positive_moisture_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    do_cldice: int,
    dt: float,
    latvap: float,
    latice: float,
    cpair: float,
    dp_p: cobj,
    qvmin_p: cobj,
    qlmin_p: cobj,
    qimin_p: cobj,
    qv_p: cobj,
    ql_p: cobj,
    qi_p: cobj,
    t_p: cobj,
    qvten_p: cobj,
    qlten_p: cobj,
    qiten_p: cobj,
    tten_p: cobj,
):
    dp = Ptr[float](dp_p)
    qvmin = Ptr[float](qvmin_p)
    qlmin = Ptr[float](qlmin_p)
    qimin = Ptr[float](qimin_p)
    qv = Ptr[float](qv_p)
    ql = Ptr[float](ql_p)
    qi = Ptr[float](qi_p)
    t = Ptr[float](t_p)
    qvten = Ptr[float](qvten_p)
    qlten = Ptr[float](qlten_p)
    qiten = Ptr[float](qiten_p)
    tten = Ptr[float](tten_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            tten[idx] = 0.0
            qvten[idx] = 0.0
            qlten[idx] = 0.0
            qiten[idx] = 0.0

    for i in range(1, ncol + 1):
        needs_fix = False
        for k in range(top_lev, pver + 1):
            idx = _idx2(i, k, pcols)
            if qv[idx] < qvmin[idx] or ql[idx] < qlmin[idx] or qi[idx] < qimin[idx]:
                needs_fix = True
                break
        if not needs_fix:
            continue

        dqv = 0.0
        for k in range(top_lev, pver + 1):
            idx = _idx2(i, k, pcols)

            dql = qlmin[idx] - ql[idx]
            if dql < 0.0:
                dql = 0.0

            if do_cldice != 0:
                dqi = qimin[idx] - qi[idx]
                if dqi < 0.0:
                    dqi = 0.0
            else:
                dqi = 0.0

            qlten[idx] = qlten[idx] + dql / dt
            qiten[idx] = qiten[idx] + dqi / dt
            qvten[idx] = qvten[idx] - (dql + dqi) / dt
            tten[idx] = tten[idx] + (latvap / cpair) * (dql / dt) + ((latvap + latice) / cpair) * (dqi / dt)
            ql[idx] = ql[idx] + dql
            qi[idx] = qi[idx] + dqi
            qv[idx] = qv[idx] - dql - dqi
            t[idx] = t[idx] + (latvap * dql + (latvap + latice) * dqi) / cpair

            dqv = qvmin[idx] - qv[idx]
            if dqv < 0.0:
                dqv = 0.0
            qvten[idx] = qvten[idx] + dqv / dt
            qv[idx] = qv[idx] + dqv
            if k != pver:
                idx_next = _idx2(i, k + 1, pcols)
                transfer = dqv * dp[idx] / dp[idx_next]
                qv[idx_next] = qv[idx_next] - transfer
                qvten[idx_next] = qvten[idx_next] - transfer / dt

            if qv[idx] < qvmin[idx]:
                qv[idx] = qvmin[idx]
            if ql[idx] < qlmin[idx]:
                ql[idx] = qlmin[idx]
            if qi[idx] < qimin[idx]:
                qi[idx] = qimin[idx]

        if dqv > 1.0e-20:
            sum_val = 0.0
            for k in range(top_lev, pver + 1):
                idx = _idx2(i, k, pcols)
                if qv[idx] > 2.0 * qvmin[idx]:
                    sum_val = sum_val + qv[idx] * dp[idx]

            denom = sum_val
            if denom < 1.0e-20:
                denom = 1.0e-20
            aa = dqv * dp[_idx2(i, pver, pcols)] / denom
            if aa < 0.5:
                for k in range(top_lev, pver + 1):
                    idx = _idx2(i, k, pcols)
                    if qv[idx] > 2.0 * qvmin[idx]:
                        dum = aa * qv[idx]
                        qv[idx] = qv[idx] - dum
                        qvten[idx] = qvten[idx] - dum / dt


def _process_rates_idx(
    i: int,
    k: int,
    idsttype: int,
    isrctype: int,
    rtype: int,
    pcols: int,
    pver: int,
    pwtype: int,
):
    return (
        (i - 1)
        + (k - 1) * pcols
        + (idsttype - 1) * pcols * pver
        + (isrctype - 1) * pcols * pver * pwtype
        + (rtype - 1) * pcols * pver * pwtype * pwtype
    )


@inline
def _iawset_idx(itype: int, iwset: int, pwtype: int):
    return (itype - 1) + (iwset - 1) * pwtype


@inline
def _iatype_idx(iwset: int, itype: int, wtrc_nwset: int):
    return (iwset - 1) + (itype - 1) * wtrc_nwset


@inline
def _spec_idx(ispec: int):
    return ispec - 1


@inline
def _bulk_idx(itype: int):
    return itype - 1


@inline
def _wtrc_ratio(ispec: int, qtrc: float, qtot: float, wtrc_qmin: float, rstd) -> float:
    if abs(qtot) < wtrc_qmin:
        return rstd[_spec_idx(ispec)]
    return qtrc / qtot


@inline
def _wiso_alpl(ispec: int, tk: float) -> float:
    if ispec <= 2:
        return 1.0
    if ispec == 3:
        return exp(
            1158.8e-12 * tk**3
            + (-1620.1e-9) * tk**2
            + 794.84e-6 * tk
            + (-161.04e-3)
            + 2.9992e6 / tk**3
        )
    return exp(0.35041e6 / tk**3 + (-1.6664e3) / tk**2 + 6.7123 / tk + (-7.685e-3))


@inline
def _wiso_alpi(ispec: int, tk: float) -> float:
    if ispec <= 2:
        return 1.0
    if ispec == 3:
        return exp(16289.0 / tk**2 + (-9.45e-2))
    return exp(11.839 / tk + (-28.224e-3))


@inline
def _wiso_ssatf(tk: float) -> float:
    ssat = 1.0 + (-0.002) * (tk - 273.16)
    if ssat < 1.0:
        ssat = 1.0
    if ssat > 2.0:
        ssat = 2.0
    return ssat


@inline
def _wiso_akci(ispec: int, tk: float, alpeq: float) -> float:
    if tk >= 253.15:
        return alpeq
    sat1 = _wiso_ssatf(tk)
    difrmj = 1.0
    if ispec == 3:
        difrmj = 0.9757
    elif ispec == 4:
        difrmj = 0.9727
    dondi = 1.0 / difrmj
    return alpeq * sat1 / (alpeq * dondi * (sat1 - 1.0) + 1.0)


@inline
def _qsat_water(t: float, p: float, epsilo: float) -> float:
    tboil = 373.16
    es = 10.0 ** (
        -7.90298 * (tboil / t - 1.0)
        + 5.02808 * log10(tboil / t)
        - 1.3816e-7 * (10.0 ** (11.344 * (1.0 - t / tboil)) - 1.0)
        + 8.1328e-3 * (10.0 ** (-3.49149 * (tboil / t - 1.0)) - 1.0)
        + log10(1013.246)
    ) * 100.0
    if (p - es) <= 0.0:
        return 1.0
    return epsilo * es / (p - (1.0 - epsilo) * es)


@inline
def _wtrc_get_alpha(
    q: float,
    tk: float,
    ispec: int,
    isrctype: int,
    idsttype: int,
    rhclc: int,
    porqh: float,
    kin: int,
    wisotope: int,
    iwtvap: int,
    iwtliq: int,
    epsilo: float,
) -> float:
    if wisotope == 0:
        return 1.0

    rh = porqh
    if rhclc != 0:
        rh = q / _qsat_water(tk, porqh, epsilo)

    alpha = 1.0
    if isrctype != idsttype:
        if isrctype == iwtvap:
            if idsttype == iwtliq:
                alpha = _wiso_alpl(ispec, tk)
            else:
                alpha = _wiso_alpi(ispec, tk)
                if kin != 0:
                    alpha = _wiso_akci(ispec, tk, alpha)
        elif idsttype == iwtvap:
            if isrctype == iwtliq:
                alpha = _wiso_alpl(ispec, tk)
                alpha = 1.0 / alpha
            else:
                alpha = 1.0
    return alpha


@inline
def _wtrc_efac(alpha: float, vapnew: float, liqnew: float, wtrc_qmin: float, rstd) -> float:
    alov = _wtrc_ratio(1, vapnew, vapnew + liqnew, wtrc_qmin, rstd)
    alov = alpha * (1.0 / alov - 1.0)
    efac = 1.0 / (alov + 1.0)
    if efac < 0.0:
        efac = 0.0
    if efac > 1.0:
        efac = 1.0
    return efac


@inline
def _wtrc_dqequil(
    alpha: float,
    feq0: float,
    vtotnew: float,
    ltotnew: float,
    visoold: float,
    lisoold: float,
    wtrc_qmin: float,
    rstd,
) -> float:
    qiso = visoold + lisoold
    vieql = qiso * _wtrc_efac(alpha, vtotnew, ltotnew, wtrc_qmin, rstd)
    vinof = qiso * _wtrc_efac(1.0, vtotnew, ltotnew, wtrc_qmin, rstd)
    visonew = feq0 * vieql + (1.0 - feq0) * vinof
    dviso = visonew - visoold
    if dviso < 0.0:
        if dviso < (-visoold):
            dviso = -visoold
    else:
        if dviso > lisoold:
            dviso = lisoold
    return dviso


@inline
def _wtrc_liqvap_equil(
    alpha: float,
    feq0: float,
    vaptot: float,
    liqtot: float,
    vapiso: float,
    liqiso: float,
    wtrc_qmin: float,
    rstd,
):
    qtiny = 1.0e-36
    qtot = vaptot + liqtot
    qiso = vapiso + liqiso

    if qtot < qtiny or qiso < qtiny:
        return vapiso, liqiso

    if liqtot < qtiny:
        dliqiso = -liqiso
        vapiso = vapiso - dliqiso
        liqiso = 0.0
        return vapiso, liqiso

    if vaptot < qtiny:
        dliqiso = vapiso
        vapiso = 0.0
        liqiso = liqiso + dliqiso
        return vapiso, liqiso

    dviso = _wtrc_dqequil(alpha, feq0, vaptot, liqtot, vapiso, liqiso, wtrc_qmin, rstd)
    dliqiso = -dviso
    liqiso = liqiso + dliqiso
    vapiso = vapiso - dliqiso
    return vapiso, liqiso


@export
def macrop_driver_wtrc_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    pwtype: int,
    top_lev: int,
    wtrc_niter: int,
    wtrc_ncnst: int,
    wtrc_nwset: int,
    wisotope: int,
    iwtvap: int,
    iwtliq: int,
    iwtice: int,
    cpair: float,
    dtime: float,
    wtrc_qmin: float,
    epsilo: float,
    state_q_p: cobj,
    state_t_p: cobj,
    state_pmid_p: cobj,
    ptend_q_p: cobj,
    qvlat_p: cobj,
    qcten_p: cobj,
    qiten_p: cobj,
    prelat_p: cobj,
    process_rates_p: cobj,
    qloc_p: cobj,
    qloc0_p: cobj,
    tloc_p: cobj,
    diff_p: cobj,
    wtrc_iawset_p: cobj,
    wtrc_iatype_p: cobj,
    wtrc_bulk_indices_p: cobj,
    wtrc_indices_p: cobj,
    iwspec_p: cobj,
    rstd_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    state_t = Ptr[float](state_t_p)
    state_pmid = Ptr[float](state_pmid_p)
    ptend_q = Ptr[float](ptend_q_p)
    qvlat = Ptr[float](qvlat_p)
    qcten = Ptr[float](qcten_p)
    qiten = Ptr[float](qiten_p)
    prelat = Ptr[float](prelat_p)
    process_rates = Ptr[float](process_rates_p)
    qloc = Ptr[float](qloc_p)
    qloc0 = Ptr[float](qloc0_p)
    tloc = Ptr[float](tloc_p)
    diff = Ptr[float](diff_p)
    wtrc_iawset = Ptr[int](wtrc_iawset_p)
    wtrc_iatype = Ptr[int](wtrc_iatype_p)
    wtrc_bulk_indices = Ptr[int](wtrc_bulk_indices_p)
    wtrc_indices = Ptr[int](wtrc_indices_p)
    iwspec = Ptr[int](iwspec_p)
    rstd = Ptr[float](rstd_p)

    for rtype in range(1, pwtype + 1):
        for isrctype in range(1, pwtype + 1):
            for idsttype in range(1, pwtype + 1):
                for k in range(top_lev, pver + 1):
                    for i in range(1, pcols + 1):
                        process_rates[
                            _process_rates_idx(
                                i,
                                k,
                                idsttype,
                                isrctype,
                                rtype,
                                pcols,
                                pver,
                                pwtype,
                            )
                        ] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)

            rate_val = qvlat[idx2] + qcten[idx2]
            rate_val = rate_val + qiten[idx2]
            process_rates[
                _process_rates_idx(
                    i, k, iwtvap, iwtvap, iwtvap, pcols, pver, pwtype
                )
            ] = process_rates[
                _process_rates_idx(
                    i, k, iwtvap, iwtvap, iwtvap, pcols, pver, pwtype
                )
            ] + rate_val

            rate_val = qcten[idx2]
            if rate_val < 0.0:
                dst_idx = _process_rates_idx(
                    i, k, iwtliq, iwtvap, iwtliq, pcols, pver, pwtype
                )
                src_idx = _process_rates_idx(
                    i, k, iwtvap, iwtliq, iwtliq, pcols, pver, pwtype
                )
            else:
                dst_idx = _process_rates_idx(
                    i, k, iwtliq, iwtvap, iwtvap, pcols, pver, pwtype
                )
                src_idx = _process_rates_idx(
                    i, k, iwtvap, iwtliq, iwtvap, pcols, pver, pwtype
                )
            process_rates[dst_idx] = process_rates[dst_idx] + rate_val
            process_rates[src_idx] = process_rates[src_idx] - rate_val

            rate_val = qiten[idx2]
            if rate_val < 0.0:
                dst_idx = _process_rates_idx(
                    i, k, iwtice, iwtvap, iwtice, pcols, pver, pwtype
                )
                src_idx = _process_rates_idx(
                    i, k, iwtvap, iwtice, iwtice, pcols, pver, pwtype
                )
            else:
                dst_idx = _process_rates_idx(
                    i, k, iwtice, iwtvap, iwtvap, pcols, pver, pwtype
                )
                src_idx = _process_rates_idx(
                    i, k, iwtvap, iwtice, iwtvap, pcols, pver, pwtype
                )
            process_rates[dst_idx] = process_rates[dst_idx] + rate_val
            process_rates[src_idx] = process_rates[src_idx] - rate_val

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            tloc[idx2] = state_t[idx2]
            for icnst in range(1, wtrc_ncnst + 1):
                trc_idx = wtrc_indices[icnst - 1]
                idx3 = _idx3(i, k, trc_idx, pcols, pver)
                qloc[idx3] = state_q[idx3]
                qloc0[idx3] = qloc[idx3]

    for iter_idx in range(1, wtrc_niter + 1):
        for k in range(top_lev, pver + 1):
            for i in range(1, ncol + 1):
                idx2 = _idx2(i, k, pcols)
                tloc[idx2] = (tloc[idx2] + (tloc[idx2] + prelat[idx2] / cpair * dtime)) / 2.0

                for isrctype in range(1, pwtype + 1):
                    for idsttype in range(1, pwtype + 1):
                        rtype = isrctype
                        rate_val = process_rates[
                            _process_rates_idx(
                                i,
                                k,
                                idsttype,
                                isrctype,
                                rtype,
                                pcols,
                                pver,
                                pwtype,
                            )
                        ]
                        if rate_val > 0.0:
                            for iwset in range(1, wtrc_nwset + 1):
                                msrc = wtrc_iawset[_iawset_idx(isrctype, iwset, pwtype)]
                                mbase = wtrc_iawset[_iawset_idx(isrctype, 1, pwtype)]
                                mdst = wtrc_iawset[_iawset_idx(idsttype, iwset, pwtype)]

                                idx_msrc = _idx3(i, k, msrc, pcols, pver)
                                idx_mbase = _idx3(i, k, mbase, pcols, pver)
                                idx_mdst = _idx3(i, k, mdst, pcols, pver)

                                R = _wtrc_ratio(
                                    iwspec[msrc - 1],
                                    qloc0[idx_msrc],
                                    qloc0[idx_mbase],
                                    wtrc_qmin,
                                    rstd,
                                )

                                if (
                                    wisotope != 0
                                    and iwset != 1
                                    and isrctype == iwtvap
                                    and idsttype == iwtice
                                ):
                                    std_vap_idx = _idx3(
                                        i,
                                        k,
                                        wtrc_iawset[_iawset_idx(iwtvap, 1, pwtype)],
                                        pcols,
                                        pver,
                                    )
                                    ispec = iwspec[mdst - 1]
                                    alpha = _wtrc_get_alpha(
                                        qloc0[std_vap_idx],
                                        tloc[idx2],
                                        ispec,
                                        isrctype,
                                        idsttype,
                                        1,
                                        state_pmid[idx2],
                                        1,
                                        wisotope,
                                        iwtvap,
                                        iwtliq,
                                        epsilo,
                                    )
                                    fr = qloc[idx_mbase] / qloc0[idx_mbase]
                                    if fr < 0.0:
                                        fr = 0.0
                                    if fr > 1.0:
                                        fr = 1.0
                                    qloc[idx_msrc] = qloc0[idx_msrc] * (fr**alpha)
                                    qloc[idx_mdst] = qloc[idx_mdst] + (qloc0[idx_msrc] - qloc[idx_msrc])
                                else:
                                    qloc[idx_mdst] = (
                                        qloc[idx_mdst]
                                        + R * rate_val * dtime / wtrc_niter
                                    )
                                    if isrctype != idsttype:
                                        qloc[idx_msrc] = (
                                            qloc[idx_msrc]
                                            - R * rate_val * dtime / wtrc_niter
                                        )

                            for icnst in range(1, wtrc_ncnst + 1):
                                trc_idx = wtrc_indices[icnst - 1]
                                idx3 = _idx3(i, k, trc_idx, pcols, pver)
                                qloc0[idx3] = qloc[idx3]

                    if wisotope != 0:
                        for iwset in range(2, wtrc_nwset + 1):
                            std_vap = wtrc_iawset[_iawset_idx(iwtvap, 1, pwtype)]
                            std_liq = wtrc_iawset[_iawset_idx(iwtliq, 1, pwtype)]
                            iso_vap = wtrc_iawset[_iawset_idx(iwtvap, iwset, pwtype)]
                            iso_liq = wtrc_iawset[_iawset_idx(iwtliq, iwset, pwtype)]

                            idx_std_vap = _idx3(i, k, std_vap, pcols, pver)
                            idx_std_liq = _idx3(i, k, std_liq, pcols, pver)
                            idx_iso_vap = _idx3(i, k, iso_vap, pcols, pver)
                            idx_iso_liq = _idx3(i, k, iso_liq, pcols, pver)

                            alpha = _wtrc_get_alpha(
                                qloc0[idx_std_vap],
                                tloc[idx2],
                                iwspec[iso_vap - 1],
                                iwtvap,
                                iwtliq,
                                0,
                                1.0,
                                0,
                                wisotope,
                                iwtvap,
                                iwtliq,
                                epsilo,
                            )
                            vapiso, liqiso = _wtrc_liqvap_equil(
                                alpha,
                                1.0,
                                qloc[idx_std_vap],
                                qloc[idx_std_liq],
                                qloc[idx_iso_vap],
                                qloc[idx_iso_liq],
                                wtrc_qmin,
                                rstd,
                            )
                            qloc[idx_iso_vap] = vapiso
                            qloc[idx_iso_liq] = liqiso

                        for icnst in range(1, wtrc_ncnst + 1):
                            trc_idx = wtrc_indices[icnst - 1]
                            idx3 = _idx3(i, k, trc_idx, pcols, pver)
                            qloc0[idx3] = qloc[idx3]

                tloc[idx2] = state_t[idx2] + prelat[idx2] / cpair * dtime / wtrc_niter

    for itype in range(1, pwtype + 1):
        for k in range(top_lev, pver + 1):
            for i in range(1, ncol + 1):
                diff[_idx3(i, k, itype, pcols, pver)] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            for icnst in range(1, wtrc_ncnst + 1):
                trc_idx = wtrc_indices[icnst - 1]
                idx3 = _idx3(i, k, trc_idx, pcols, pver)
                ptend_q[idx3] = (qloc[idx3] - state_q[idx3]) / dtime

                if icnst <= pwtype:
                    bulk_idx = wtrc_bulk_indices[_bulk_idx(icnst)]
                    diff[_idx3(i, k, icnst, pcols, pver)] = ptend_q[idx3] - ptend_q[
                        _idx3(i, k, bulk_idx, pcols, pver)
                    ]

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            for itype in range(1, pwtype + 1):
                qtmp = 0.0
                diff_idx = _idx3(i, k, itype, pcols, pver)
                for iwset in range(1, wtrc_nwset + 1):
                    trc_idx = wtrc_iatype[_iatype_idx(iwset, itype, wtrc_nwset)]
                    idx3 = _idx3(i, k, trc_idx, pcols, pver)
                    if iwset == 1:
                        qtmp = ptend_q[idx3]
                    R = _wtrc_ratio(
                        iwspec[trc_idx - 1],
                        ptend_q[idx3],
                        qtmp,
                        wtrc_qmin,
                        rstd,
                    )
                    ptend_q[idx3] = ptend_q[idx3] - R * diff[diff_idx]

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            for itype in range(1, pwtype + 1):
                qtmp = 0.0
                diff_idx = _idx3(i, k, itype, pcols, pver)
                for iwset in range(1, wtrc_nwset + 1):
                    trc_idx = wtrc_iatype[_iatype_idx(iwset, itype, wtrc_nwset)]
                    idx3 = _idx3(i, k, trc_idx, pcols, pver)
                    if iwset == 1:
                        bulk_idx = wtrc_bulk_indices[_bulk_idx(itype)]
                        diff[diff_idx] = ptend_q[idx3] - ptend_q[
                            _idx3(i, k, bulk_idx, pcols, pver)
                        ]
                        qtmp = ptend_q[idx3]
                    R = _wtrc_ratio(
                        iwspec[trc_idx - 1],
                        ptend_q[idx3],
                        qtmp,
                        wtrc_qmin,
                        rstd,
                    )
                    ptend_q[idx3] = ptend_q[idx3] - R * diff[diff_idx]


@export
def macrop_driver_select_branches_codon(
    micro_do_icesupersat: int,
    trace_water: int,
    wtrc_detrain_in_macrop: int,
    cu_det_st: int,
    use_shfrc: int,
    do_cldice: int,
    do_cldliq: int,
    do_detrain: int,
    branch_mask_p: cobj,
):
    branch_mask = Ptr[int](branch_mask_p)

    mask = 0
    if micro_do_icesupersat != 0:
        mask |= 1
    if trace_water != 0:
        mask |= 2
    if wtrc_detrain_in_macrop != 0:
        mask |= 4
    if cu_det_st != 0:
        mask |= 8
    if use_shfrc != 0:
        mask |= 16
    if do_cldice != 0:
        mask |= 32
    if do_cldliq != 0:
        mask |= 64
    if do_detrain != 0:
        mask |= 128

    branch_mask[0] = mask


@export
def macrop_driver_ptend_lq_mask_shell_codon(
    mode: int,
    pcnst: int,
    wtrc_nwset: int,
    wtrc_ncnst: int,
    use_water_tracers: int,
    ixcldliq: int,
    ixcldice: int,
    ixnumliq: int,
    ixnumice: int,
    lq_mask_p: cobj,
    liq_type_p: cobj,
    ice_type_p: cobj,
    wtrc_indices_p: cobj,
):
    lq_mask = Ptr[int](lq_mask_p)
    liq_type = Ptr[int](liq_type_p)
    ice_type = Ptr[int](ice_type_p)
    wtrc_indices = Ptr[int](wtrc_indices_p)

    for m in range(1, pcnst + 1):
        lq_mask[m - 1] = 0

    if mode == 1:
        lq_mask[ixcldliq - 1] = 1
        lq_mask[ixcldice - 1] = 1
        lq_mask[ixnumliq - 1] = 1
        lq_mask[ixnumice - 1] = 1
        if use_water_tracers != 0:
            for m in range(1, wtrc_nwset + 1):
                lq_mask[liq_type[m - 1] - 1] = 1
                lq_mask[ice_type[m - 1] - 1] = 1
    elif mode == 2:
        lq_mask[0] = 1
        lq_mask[ixcldice - 1] = 1
        lq_mask[ixcldliq - 1] = 1
        lq_mask[ixnumliq - 1] = 1
        lq_mask[ixnumice - 1] = 1
        for m in range(1, wtrc_ncnst + 1):
            lq_mask[wtrc_indices[m - 1] - 1] = 1


@export
def macrop_driver_detrain_init_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    dlf_T_p: cobj,
    dlf_qv_p: cobj,
    dlf_ql_p: cobj,
    dlf_qi_p: cobj,
    dlf_nl_p: cobj,
    dlf_ni_p: cobj,
    det_s_p: cobj,
    det_ice_p: cobj,
    dpdlfliq_p: cobj,
    dpdlfice_p: cobj,
    shdlfliq_p: cobj,
    shdlfice_p: cobj,
    dpdlft_p: cobj,
    shdlft_p: cobj,
):
    dlf_T = Ptr[float](dlf_T_p)
    dlf_qv = Ptr[float](dlf_qv_p)
    dlf_ql = Ptr[float](dlf_ql_p)
    dlf_qi = Ptr[float](dlf_qi_p)
    dlf_nl = Ptr[float](dlf_nl_p)
    dlf_ni = Ptr[float](dlf_ni_p)
    det_s = Ptr[float](det_s_p)
    det_ice = Ptr[float](det_ice_p)
    dpdlfliq = Ptr[float](dpdlfliq_p)
    dpdlfice = Ptr[float](dpdlfice_p)
    shdlfliq = Ptr[float](shdlfliq_p)
    shdlfice = Ptr[float](shdlfice_p)
    dpdlft = Ptr[float](dpdlft_p)
    shdlft = Ptr[float](shdlft_p)

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            idx2 = _idx2(i, k, pcols)
            dlf_T[idx2] = 0.0
            dlf_qv[idx2] = 0.0
            dlf_ql[idx2] = 0.0
            dlf_qi[idx2] = 0.0
            dlf_nl[idx2] = 0.0
            dlf_ni[idx2] = 0.0
            dpdlfliq[idx2] = 0.0
            dpdlfice[idx2] = 0.0
            shdlfliq[idx2] = 0.0
            shdlfice[idx2] = 0.0
            dpdlft[idx2] = 0.0
            shdlft[idx2] = 0.0

    for i in range(1, pcols + 1):
        idx1 = i - 1
        det_s[idx1] = 0.0
        det_ice[idx1] = 0.0


@export
def macrop_driver_detrain_init_lq_mask_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    wtrc_nwset: int,
    wtrc_ncnst: int,
    use_water_tracers: int,
    ixcldliq: int,
    ixcldice: int,
    ixnumliq: int,
    ixnumice: int,
    dlf_T_p: cobj,
    dlf_qv_p: cobj,
    dlf_ql_p: cobj,
    dlf_qi_p: cobj,
    dlf_nl_p: cobj,
    dlf_ni_p: cobj,
    det_s_p: cobj,
    det_ice_p: cobj,
    dpdlfliq_p: cobj,
    dpdlfice_p: cobj,
    shdlfliq_p: cobj,
    shdlfice_p: cobj,
    dpdlft_p: cobj,
    shdlft_p: cobj,
    lq_mask_p: cobj,
    liq_type_p: cobj,
    ice_type_p: cobj,
    wtrc_indices_p: cobj,
):
    macrop_driver_detrain_init_shell_codon(
        ncol,
        pcols,
        pver,
        dlf_T_p,
        dlf_qv_p,
        dlf_ql_p,
        dlf_qi_p,
        dlf_nl_p,
        dlf_ni_p,
        det_s_p,
        det_ice_p,
        dpdlfliq_p,
        dpdlfice_p,
        shdlfliq_p,
        shdlfice_p,
        dpdlft_p,
        shdlft_p,
    )
    macrop_driver_ptend_lq_mask_shell_codon(
        1,
        pcnst,
        wtrc_nwset,
        wtrc_ncnst,
        use_water_tracers,
        ixcldliq,
        ixcldice,
        ixnumliq,
        ixnumice,
        lq_mask_p,
        liq_type_p,
        ice_type_p,
        wtrc_indices_p,
    )


@export
def macrop_driver_detrain_post_shell_codon(
    ncol: int,
    pcols: int,
    det_ice_p: cobj,
):
    det_ice = Ptr[float](det_ice_p)

    for i in range(1, ncol + 1):
        idx1 = i - 1
        det_ice[idx1] = det_ice[idx1] / 1000.0


@export
def macrop_driver_mmacro_input_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    top_lev: int,
    ixcldliq: int,
    ixcldice: int,
    ixnumliq: int,
    ixnumice: int,
    state_q_p: cobj,
    zeros_p: cobj,
    qc_p: cobj,
    qi_p: cobj,
    nc_p: cobj,
    ni_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    zeros = Ptr[float](zeros_p)
    qc = Ptr[float](qc_p)
    qi = Ptr[float](qi_p)
    nc = Ptr[float](nc_p)
    ni = Ptr[float](ni_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            zeros[idx2] = 0.0
            qc[idx2] = state_q[_idx3(i, k, ixcldliq, pcols, pver)]
            qi[idx2] = state_q[_idx3(i, k, ixcldice, pcols, pver)]
            nc[idx2] = state_q[_idx3(i, k, ixnumliq, pcols, pver)]
            ni[idx2] = state_q[_idx3(i, k, ixnumice, pcols, pver)]


@export
def macrop_driver_mmacro_post_fields_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    fice_p: cobj,
    alst_p: cobj,
    aist_p: cobj,
    fice_ql_p: cobj,
    ast_p: cobj,
):
    fice = Ptr[float](fice_p)
    alst = Ptr[float](alst_p)
    aist = Ptr[float](aist_p)
    fice_ql = Ptr[float](fice_ql_p)
    ast = Ptr[float](ast_p)

    for k in range(1, top_lev):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            fice_ql[idx2] = 0.0
            ast[idx2] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            fice_ql[idx2] = fice[idx2]
            ast[idx2] = max(alst[idx2], aist[idx2])


@export
def macrop_driver_mmacro_config_check_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    do_cldice: int,
    do_cldliq: int,
    qiten_p: cobj,
    niten_p: cobj,
    qcten_p: cobj,
    ncten_p: cobj,
    mask_p: cobj,
):
    qiten = Ptr[float](qiten_p)
    niten = Ptr[float](niten_p)
    qcten = Ptr[float](qcten_p)
    ncten = Ptr[float](ncten_p)
    mask_out = Ptr[int](mask_p)

    mask = 0
    if do_cldice == 0:
        qiten_nonzero = False
        niten_nonzero = False
        for k in range(top_lev, pver + 1):
            for i in range(1, ncol + 1):
                idx2 = _idx2(i, k, pcols)
                if qiten[idx2] != 0.0:
                    qiten_nonzero = True
                if niten[idx2] != 0.0:
                    niten_nonzero = True
        if qiten_nonzero:
            mask |= 1
        if niten_nonzero:
            mask |= 2

    if do_cldliq == 0:
        qcten_nonzero = False
        ncten_nonzero = False
        for k in range(top_lev, pver + 1):
            for i in range(1, ncol + 1):
                idx2 = _idx2(i, k, pcols)
                if qcten[idx2] != 0.0:
                    qcten_nonzero = True
                if ncten[idx2] != 0.0:
                    ncten_nonzero = True
        if qcten_nonzero:
            mask |= 4
        if ncten_nonzero:
            mask |= 8

    mask_out[0] = mask


@export
def macrop_driver_cfmip_diag_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    cld_p: cobj,
    state_ql_p: cobj,
    state_qi_p: cobj,
    mr_ccliq_p: cobj,
    mr_ccice_p: cobj,
    mr_lsliq_p: cobj,
    mr_lsice_p: cobj,
):
    cld = Ptr[float](cld_p)
    state_ql = Ptr[float](state_ql_p)
    state_qi = Ptr[float](state_qi_p)
    mr_ccliq = Ptr[float](mr_ccliq_p)
    mr_ccice = Ptr[float](mr_ccice_p)
    mr_lsliq = Ptr[float](mr_lsliq_p)
    mr_lsice = Ptr[float](mr_lsice_p)

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            idx2 = _idx2(i, k, pcols)
            mr_ccliq[idx2] = 0.0
            mr_ccice[idx2] = 0.0
            mr_lsliq[idx2] = 0.0
            mr_lsice[idx2] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            if cld[idx2] > 0.0:
                mr_lsliq[idx2] = state_ql[idx2]
                mr_lsice[idx2] = state_qi[idx2]
            else:
                mr_lsliq[idx2] = 0.0
                mr_lsice[idx2] = 0.0


@export
def macrop_driver_wtrc_detrain_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    wtrc_nwset: int,
    state_t_p: cobj,
    wtdlf_p: cobj,
    liq_type_p: cobj,
    ice_type_p: cobj,
    ptend_q_p: cobj,
):
    state_t = Ptr[float](state_t_p)
    wtdlf = Ptr[float](wtdlf_p)
    liq_type = Ptr[int](liq_type_p)
    ice_type = Ptr[int](ice_type_p)
    ptend_q = Ptr[float](ptend_q_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            if state_t[idx2] > 268.15:
                dum1 = 0.0
            elif state_t[idx2] < 238.15:
                dum1 = 1.0
            else:
                dum1 = (268.15 - state_t[idx2]) / 30.0
            for m in range(1, wtrc_nwset + 1):
                idx_wtdlf = _idx3(i, k, m, pcols, pver)
                ptend_q[_idx3(i, k, liq_type[m - 1], pcols, pver)] = wtdlf[idx_wtdlf] * (1.0 - dum1)
                ptend_q[_idx3(i, k, ice_type[m - 1], pcols, pver)] = wtdlf[idx_wtdlf] * dum1


@export
def macrop_driver_clr_old_diag_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    concld_p: cobj,
    alst_p: cobj,
    ast_p: cobj,
    concld_old_p: cobj,
    clrw_old_p: cobj,
    clri_old_p: cobj,
):
    concld = Ptr[float](concld_p)
    alst = Ptr[float](alst_p)
    ast = Ptr[float](ast_p)
    concld_old = Ptr[float](concld_old_p)
    clrw_old = Ptr[float](clrw_old_p)
    clri_old = Ptr[float](clri_old_p)

    for k in range(1, top_lev):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            clrw_old[idx2] = 0.0
            clri_old[idx2] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            concld_old[idx2] = concld[idx2]
            clrw_old[idx2] = max(0.0, min(1.0, 1.0 - concld[idx2] - alst[idx2]))
            clri_old[idx2] = max(0.0, min(1.0, 1.0 - concld[idx2] - ast[idx2]))


@export
def macrop_driver_forcing_prep_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    nstep: int,
    rdtime: float,
    state_t_p: cobj,
    state_qv_p: cobj,
    qc_p: cobj,
    qi_p: cobj,
    nc_p: cobj,
    ni_p: cobj,
    tcwat_p: cobj,
    qcwat_p: cobj,
    lcwat_p: cobj,
    iccwat_p: cobj,
    nlwat_p: cobj,
    niwat_p: cobj,
    cc_t_p: cobj,
    cc_qv_p: cobj,
    cc_ql_p: cobj,
    cc_qi_p: cobj,
    cc_nl_p: cobj,
    cc_ni_p: cobj,
    cc_qlst_p: cobj,
    ttend_p: cobj,
    qtend_p: cobj,
    ltend_p: cobj,
    itend_p: cobj,
    nltend_p: cobj,
    nitend_p: cobj,
    lmitend_p: cobj,
    t_inout_p: cobj,
    qv_inout_p: cobj,
    ql_inout_p: cobj,
    qi_inout_p: cobj,
    nl_inout_p: cobj,
    ni_inout_p: cobj,
):
    state_t = Ptr[float](state_t_p)
    state_qv = Ptr[float](state_qv_p)
    qc = Ptr[float](qc_p)
    qi = Ptr[float](qi_p)
    nc = Ptr[float](nc_p)
    ni = Ptr[float](ni_p)
    tcwat = Ptr[float](tcwat_p)
    qcwat = Ptr[float](qcwat_p)
    lcwat = Ptr[float](lcwat_p)
    iccwat = Ptr[float](iccwat_p)
    nlwat = Ptr[float](nlwat_p)
    niwat = Ptr[float](niwat_p)
    cc_t = Ptr[float](cc_t_p)
    cc_qv = Ptr[float](cc_qv_p)
    cc_ql = Ptr[float](cc_ql_p)
    cc_qi = Ptr[float](cc_qi_p)
    cc_nl = Ptr[float](cc_nl_p)
    cc_ni = Ptr[float](cc_ni_p)
    cc_qlst = Ptr[float](cc_qlst_p)
    ttend = Ptr[float](ttend_p)
    qtend = Ptr[float](qtend_p)
    ltend = Ptr[float](ltend_p)
    itend = Ptr[float](itend_p)
    nltend = Ptr[float](nltend_p)
    nitend = Ptr[float](nitend_p)
    lmitend = Ptr[float](lmitend_p)
    t_inout = Ptr[float](t_inout_p)
    qv_inout = Ptr[float](qv_inout_p)
    ql_inout = Ptr[float](ql_inout_p)
    qi_inout = Ptr[float](qi_inout_p)
    nl_inout = Ptr[float](nl_inout_p)
    ni_inout = Ptr[float](ni_inout_p)

    if nstep <= 1:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx2 = _idx2(i, k, pcols)
                tcwat[idx2] = state_t[idx2]
                qcwat[idx2] = state_qv[idx2]
                lcwat[idx2] = qc[idx2] + qi[idx2]
                iccwat[idx2] = qi[idx2]
                nlwat[idx2] = nc[idx2]
                niwat[idx2] = ni[idx2]
                ttend[idx2] = 0.0
                qtend[idx2] = 0.0
                ltend[idx2] = 0.0
                itend[idx2] = 0.0
                nltend[idx2] = 0.0
                nitend[idx2] = 0.0
                cc_t[idx2] = 0.0
                cc_qv[idx2] = 0.0
                cc_ql[idx2] = 0.0
                cc_qi[idx2] = 0.0
                cc_nl[idx2] = 0.0
                cc_ni[idx2] = 0.0
                cc_qlst[idx2] = 0.0
    else:
        for k in range(top_lev, pver + 1):
            for i in range(1, ncol + 1):
                idx2 = _idx2(i, k, pcols)
                ttend[idx2] = (state_t[idx2] - tcwat[idx2]) * rdtime - cc_t[idx2]
                qtend[idx2] = (state_qv[idx2] - qcwat[idx2]) * rdtime - cc_qv[idx2]
                ltend[idx2] = (qc[idx2] + qi[idx2] - lcwat[idx2]) * rdtime - (cc_ql[idx2] + cc_qi[idx2])
                itend[idx2] = (qi[idx2] - iccwat[idx2]) * rdtime - cc_qi[idx2]
                nltend[idx2] = (nc[idx2] - nlwat[idx2]) * rdtime - cc_nl[idx2]
                nitend[idx2] = (ni[idx2] - niwat[idx2]) * rdtime - cc_ni[idx2]

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            lmitend[idx2] = ltend[idx2] - itend[idx2]
            t_inout[idx2] = tcwat[idx2]
            qv_inout[idx2] = qcwat[idx2]
            ql_inout[idx2] = lcwat[idx2] - iccwat[idx2]
            qi_inout[idx2] = iccwat[idx2]
            nl_inout[idx2] = nlwat[idx2]
            ni_inout[idx2] = niwat[idx2]


@export
def macrop_driver_mmacro_prepare_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    top_lev: int,
    ixcldliq: int,
    ixcldice: int,
    ixnumliq: int,
    ixnumice: int,
    nstep: int,
    rdtime: float,
    state_q_p: cobj,
    state_t_p: cobj,
    state_qv_p: cobj,
    zeros_p: cobj,
    qc_p: cobj,
    qi_p: cobj,
    nc_p: cobj,
    ni_p: cobj,
    tcwat_p: cobj,
    qcwat_p: cobj,
    lcwat_p: cobj,
    iccwat_p: cobj,
    nlwat_p: cobj,
    niwat_p: cobj,
    cc_t_p: cobj,
    cc_qv_p: cobj,
    cc_ql_p: cobj,
    cc_qi_p: cobj,
    cc_nl_p: cobj,
    cc_ni_p: cobj,
    cc_qlst_p: cobj,
    ttend_p: cobj,
    qtend_p: cobj,
    ltend_p: cobj,
    itend_p: cobj,
    nltend_p: cobj,
    nitend_p: cobj,
    lmitend_p: cobj,
    t_inout_p: cobj,
    qv_inout_p: cobj,
    ql_inout_p: cobj,
    qi_inout_p: cobj,
    nl_inout_p: cobj,
    ni_inout_p: cobj,
):
    macrop_driver_mmacro_input_shell_codon(
        ncol,
        pcols,
        pver,
        pcnst,
        top_lev,
        ixcldliq,
        ixcldice,
        ixnumliq,
        ixnumice,
        state_q_p,
        zeros_p,
        qc_p,
        qi_p,
        nc_p,
        ni_p,
    )
    macrop_driver_forcing_prep_codon(
        ncol,
        pcols,
        pver,
        top_lev,
        nstep,
        rdtime,
        state_t_p,
        state_qv_p,
        qc_p,
        qi_p,
        nc_p,
        ni_p,
        tcwat_p,
        qcwat_p,
        lcwat_p,
        iccwat_p,
        nlwat_p,
        niwat_p,
        cc_t_p,
        cc_qv_p,
        cc_ql_p,
        cc_qi_p,
        cc_nl_p,
        cc_ni_p,
        cc_qlst_p,
        ttend_p,
        qtend_p,
        ltend_p,
        itend_p,
        nltend_p,
        nitend_p,
        lmitend_p,
        t_inout_p,
        qv_inout_p,
        ql_inout_p,
        qi_inout_p,
        nl_inout_p,
        ni_inout_p,
    )


@export
def macrop_driver_ptend_assign_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    tlat_p: cobj,
    qvlat_p: cobj,
    qcten_p: cobj,
    qiten_p: cobj,
    ncten_p: cobj,
    niten_p: cobj,
    ptend_s_p: cobj,
    ptend_qv_p: cobj,
    ptend_ql_p: cobj,
    ptend_qi_p: cobj,
    ptend_nl_p: cobj,
    ptend_ni_p: cobj,
):
    tlat = Ptr[float](tlat_p)
    qvlat = Ptr[float](qvlat_p)
    qcten = Ptr[float](qcten_p)
    qiten = Ptr[float](qiten_p)
    ncten = Ptr[float](ncten_p)
    niten = Ptr[float](niten_p)
    ptend_s = Ptr[float](ptend_s_p)
    ptend_qv = Ptr[float](ptend_qv_p)
    ptend_ql = Ptr[float](ptend_ql_p)
    ptend_qi = Ptr[float](ptend_qi_p)
    ptend_nl = Ptr[float](ptend_nl_p)
    ptend_ni = Ptr[float](ptend_ni_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            ptend_s[idx2] = tlat[idx2]
            ptend_qv[idx2] = qvlat[idx2]
            ptend_ql[idx2] = qcten[idx2]
            ptend_qi[idx2] = qiten[idx2]
            ptend_nl[idx2] = ncten[idx2]
            ptend_ni[idx2] = niten[idx2]


@export
def macrop_driver_ptend_config_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    do_cldice: int,
    do_cldliq: int,
    tlat_p: cobj,
    qvlat_p: cobj,
    qcten_p: cobj,
    qiten_p: cobj,
    ncten_p: cobj,
    niten_p: cobj,
    ptend_s_p: cobj,
    ptend_qv_p: cobj,
    ptend_ql_p: cobj,
    ptend_qi_p: cobj,
    ptend_nl_p: cobj,
    ptend_ni_p: cobj,
    mask_p: cobj,
):
    macrop_driver_ptend_assign_codon(
        ncol,
        pcols,
        pver,
        top_lev,
        tlat_p,
        qvlat_p,
        qcten_p,
        qiten_p,
        ncten_p,
        niten_p,
        ptend_s_p,
        ptend_qv_p,
        ptend_ql_p,
        ptend_qi_p,
        ptend_nl_p,
        ptend_ni_p,
    )
    macrop_driver_mmacro_config_check_codon(
        ncol,
        pcols,
        pver,
        top_lev,
        do_cldice,
        do_cldliq,
        qiten_p,
        niten_p,
        qcten_p,
        ncten_p,
        mask_p,
    )


@export
def macrop_driver_wtrc_split_tend_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    qcten_p: cobj,
    qiten_p: cobj,
    pqctn_p: cobj,
    nqctn_p: cobj,
    pqitn_p: cobj,
    nqitn_p: cobj,
):
    qcten = Ptr[float](qcten_p)
    qiten = Ptr[float](qiten_p)
    pqctn = Ptr[float](pqctn_p)
    nqctn = Ptr[float](nqctn_p)
    pqitn = Ptr[float](pqitn_p)
    nqitn = Ptr[float](nqitn_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            if qcten[idx2] < 0.0:
                nqctn[idx2] = qcten[idx2]
            else:
                pqctn[idx2] = qcten[idx2]
            if qiten[idx2] < 0.0:
                nqitn[idx2] = qiten[idx2]
            else:
                pqitn[idx2] = qiten[idx2]


@export
def macrop_driver_wtrc_process_rates_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pwtype: int,
    top_lev: int,
    iwtvap: int,
    iwtliq: int,
    iwtice: int,
    qvlat_p: cobj,
    qcten_p: cobj,
    qiten_p: cobj,
    process_rates_p: cobj,
):
    qvlat = Ptr[float](qvlat_p)
    qcten = Ptr[float](qcten_p)
    qiten = Ptr[float](qiten_p)
    process_rates = Ptr[float](process_rates_p)

    for rtype in range(1, pwtype + 1):
        for isrctype in range(1, pwtype + 1):
            for idsttype in range(1, pwtype + 1):
                for k in range(top_lev, pver + 1):
                    for i in range(1, pcols + 1):
                        process_rates[
                            _process_rates_idx(
                                i,
                                k,
                                idsttype,
                                isrctype,
                                rtype,
                                pcols,
                                pver,
                                pwtype,
                            )
                        ] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)

            rate_val = qvlat[idx2] + qcten[idx2]
            rate_val = rate_val + qiten[idx2]
            process_rates[
                _process_rates_idx(
                    i, k, iwtvap, iwtvap, iwtvap, pcols, pver, pwtype
                )
            ] = process_rates[
                _process_rates_idx(
                    i, k, iwtvap, iwtvap, iwtvap, pcols, pver, pwtype
                )
            ] + rate_val

            rate_val = qcten[idx2]
            if rate_val < 0.0:
                dst_idx = _process_rates_idx(
                    i, k, iwtliq, iwtvap, iwtliq, pcols, pver, pwtype
                )
                src_idx = _process_rates_idx(
                    i, k, iwtvap, iwtliq, iwtliq, pcols, pver, pwtype
                )
            else:
                dst_idx = _process_rates_idx(
                    i, k, iwtliq, iwtvap, iwtvap, pcols, pver, pwtype
                )
                src_idx = _process_rates_idx(
                    i, k, iwtvap, iwtliq, iwtvap, pcols, pver, pwtype
                )
            process_rates[dst_idx] = process_rates[dst_idx] + rate_val
            process_rates[src_idx] = process_rates[src_idx] - rate_val

            rate_val = qiten[idx2]
            if rate_val < 0.0:
                dst_idx = _process_rates_idx(
                    i, k, iwtice, iwtvap, iwtice, pcols, pver, pwtype
                )
                src_idx = _process_rates_idx(
                    i, k, iwtvap, iwtice, iwtice, pcols, pver, pwtype
                )
            else:
                dst_idx = _process_rates_idx(
                    i, k, iwtice, iwtvap, iwtvap, pcols, pver, pwtype
                )
                src_idx = _process_rates_idx(
                    i, k, iwtvap, iwtice, iwtvap, pcols, pver, pwtype
                )
            process_rates[dst_idx] = process_rates[dst_idx] + rate_val
            process_rates[src_idx] = process_rates[src_idx] - rate_val


@export
def macrop_driver_cloud_mixing_diag_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    cld_p: cobj,
    state_ql_p: cobj,
    state_qi_p: cobj,
    mr_lsliq_p: cobj,
    mr_lsice_p: cobj,
):
    cld = Ptr[float](cld_p)
    state_ql = Ptr[float](state_ql_p)
    state_qi = Ptr[float](state_qi_p)
    mr_lsliq = Ptr[float](mr_lsliq_p)
    mr_lsice = Ptr[float](mr_lsice_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            if cld[idx2] > 0.0:
                mr_lsliq[idx2] = state_ql[idx2]
                mr_lsice[idx2] = state_qi[idx2]
            else:
                mr_lsliq[idx2] = 0.0
                mr_lsice[idx2] = 0.0


@export
def macrop_driver_store_state_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    tmelt: float,
    state_t_p: cobj,
    state_qv_p: cobj,
    state_ql_p: cobj,
    state_qi_p: cobj,
    state_nl_p: cobj,
    state_ni_p: cobj,
    tcwat_p: cobj,
    qcwat_p: cobj,
    lcwat_p: cobj,
    iccwat_p: cobj,
    nlwat_p: cobj,
    niwat_p: cobj,
    cldsice_p: cobj,
):
    state_t = Ptr[float](state_t_p)
    state_qv = Ptr[float](state_qv_p)
    state_ql = Ptr[float](state_ql_p)
    state_qi = Ptr[float](state_qi_p)
    state_nl = Ptr[float](state_nl_p)
    state_ni = Ptr[float](state_ni_p)
    tcwat = Ptr[float](tcwat_p)
    qcwat = Ptr[float](qcwat_p)
    lcwat = Ptr[float](lcwat_p)
    iccwat = Ptr[float](iccwat_p)
    nlwat = Ptr[float](nlwat_p)
    niwat = Ptr[float](niwat_p)
    cldsice = Ptr[float](cldsice_p)

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            cldsice[_idx2(i, k, pcols)] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            tcwat[idx2] = state_t[idx2]
            qcwat[idx2] = state_qv[idx2]
            lcwat[idx2] = state_ql[idx2] + state_qi[idx2]
            iccwat[idx2] = state_qi[idx2]
            nlwat[idx2] = state_nl[idx2]
            niwat[idx2] = state_ni[idx2]
            cldsice[idx2] = lcwat[idx2] * min(1.0, max(0.0, (tmelt - tcwat[idx2]) / 20.0))


@export
def macrop_driver_detrain_core_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    do_detrain: int,
    cu_det_st: int,
    cpair: float,
    gravit: float,
    latice: float,
    nl_denom_a: float,
    nl_denom_b: float,
    ni_denom_a: float,
    ni_denom_b: float,
    state_t_p: cobj,
    state_pdel_p: cobj,
    dlf_p: cobj,
    dlf2_p: cobj,
    ptend_ql_p: cobj,
    ptend_qi_p: cobj,
    ptend_nl_p: cobj,
    ptend_ni_p: cobj,
    ptend_s_p: cobj,
    det_s_p: cobj,
    det_ice_p: cobj,
    dlf_t_p: cobj,
    dlf_qv_p: cobj,
    dlf_ql_p: cobj,
    dlf_qi_p: cobj,
    dlf_nl_p: cobj,
    dlf_ni_p: cobj,
    dpdlfliq_p: cobj,
    dpdlfice_p: cobj,
    shdlfliq_p: cobj,
    shdlfice_p: cobj,
    dpdlft_p: cobj,
    shdlft_p: cobj,
):
    # Fortran mappings: state_t/state_pdel/dlf/... are real(r8) arrays with shape (pcols,pver);
    # det_s/det_ice are real(r8) arrays with shape (pcols).
    state_t = Ptr[float](state_t_p)
    state_pdel = Ptr[float](state_pdel_p)
    dlf = Ptr[float](dlf_p)
    dlf2 = Ptr[float](dlf2_p)
    ptend_ql = Ptr[float](ptend_ql_p)
    ptend_qi = Ptr[float](ptend_qi_p)
    ptend_nl = Ptr[float](ptend_nl_p)
    ptend_ni = Ptr[float](ptend_ni_p)
    ptend_s = Ptr[float](ptend_s_p)
    det_s = Ptr[float](det_s_p)
    det_ice = Ptr[float](det_ice_p)
    dlf_t = Ptr[float](dlf_t_p)
    dlf_qv = Ptr[float](dlf_qv_p)
    dlf_ql = Ptr[float](dlf_ql_p)
    dlf_qi = Ptr[float](dlf_qi_p)
    dlf_nl = Ptr[float](dlf_nl_p)
    dlf_ni = Ptr[float](dlf_ni_p)
    dpdlfliq = Ptr[float](dpdlfliq_p)
    dpdlfice = Ptr[float](dpdlfice_p)
    shdlfliq = Ptr[float](shdlfliq_p)
    shdlfice = Ptr[float](shdlfice_p)
    dpdlft = Ptr[float](dpdlft_p)
    shdlft = Ptr[float](shdlft_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            idx1 = i - 1

            if state_t[idx2] > 268.15:
                dum1_local = 0.0
            elif state_t[idx2] < 238.15:
                dum1_local = 1.0
            else:
                dum1_local = (268.15 - state_t[idx2]) / 30.0

            if do_detrain != 0:
                ptend_ql[idx2] = dlf[idx2] * (1.0 - dum1_local)
                ptend_qi[idx2] = dlf[idx2] * dum1_local
                ptend_nl[idx2] = (
                    3.0
                    * (max(0.0, (dlf[idx2] - dlf2[idx2])) * (1.0 - dum1_local))
                    / nl_denom_a
                    + 3.0
                    * (dlf2[idx2] * (1.0 - dum1_local))
                    / nl_denom_b
                )
                ptend_ni[idx2] = (
                    3.0
                    * (max(0.0, (dlf[idx2] - dlf2[idx2])) * dum1_local)
                    / ni_denom_a
                    + 3.0
                    * (dlf2[idx2] * dum1_local)
                    / ni_denom_b
                )
                ptend_s[idx2] = dlf[idx2] * dum1_local * latice
            else:
                ptend_ql[idx2] = 0.0
                ptend_qi[idx2] = 0.0
                ptend_nl[idx2] = 0.0
                ptend_ni[idx2] = 0.0
                ptend_s[idx2] = 0.0

            det_s[idx1] = det_s[idx1] + ptend_s[idx2] * state_pdel[idx2] / gravit
            det_ice[idx1] = det_ice[idx1] - ptend_qi[idx2] * state_pdel[idx2] / gravit

            if cu_det_st != 0:
                dlf_t[idx2] = ptend_s[idx2] / cpair
                dlf_qv[idx2] = 0.0
                dlf_ql[idx2] = ptend_ql[idx2]
                dlf_qi[idx2] = ptend_qi[idx2]
                dlf_nl[idx2] = ptend_nl[idx2]
                dlf_ni[idx2] = ptend_ni[idx2]
                ptend_ql[idx2] = 0.0
                ptend_qi[idx2] = 0.0
                ptend_nl[idx2] = 0.0
                ptend_ni[idx2] = 0.0
                ptend_s[idx2] = 0.0
                dpdlfliq[idx2] = 0.0
                dpdlfice[idx2] = 0.0
                shdlfliq[idx2] = 0.0
                shdlfice[idx2] = 0.0
                dpdlft[idx2] = 0.0
                shdlft[idx2] = 0.0
            else:
                dpdlfliq[idx2] = (dlf[idx2] - dlf2[idx2]) * (1.0 - dum1_local)
                dpdlfice[idx2] = (dlf[idx2] - dlf2[idx2]) * dum1_local
                shdlfliq[idx2] = dlf2[idx2] * (1.0 - dum1_local)
                shdlfice[idx2] = dlf2[idx2] * dum1_local
                dpdlft[idx2] = (dlf[idx2] - dlf2[idx2]) * dum1_local * latice / cpair
                shdlft[idx2] = dlf2[idx2] * dum1_local * latice / cpair


@export
def macrop_driver_stage_dispatch_codon(
    stage: int,
    mode: int,
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    top_lev: int,
    wtrc_nwset: int,
    wtrc_ncnst: int,
    use_water_tracers: int,
    ixcldliq: int,
    ixcldice: int,
    ixnumliq: int,
    ixnumice: int,
    do_detrain: int,
    cu_det_st: int,
    do_cldice: int,
    do_cldliq: int,
    nstep: int,
    cpair: float,
    gravit: float,
    latice: float,
    nl_denom_a: float,
    nl_denom_b: float,
    ni_denom_a: float,
    ni_denom_b: float,
    rdtime: float,
    tmelt: float,
    p1: cobj,
    p2: cobj,
    p3: cobj,
    p4: cobj,
    p5: cobj,
    p6: cobj,
    p7: cobj,
    p8: cobj,
    p9: cobj,
    p10: cobj,
    p11: cobj,
    p12: cobj,
    p13: cobj,
    p14: cobj,
    p15: cobj,
    p16: cobj,
    p17: cobj,
    p18: cobj,
    p19: cobj,
    p20: cobj,
    p21: cobj,
    p22: cobj,
    p23: cobj,
    p24: cobj,
    p25: cobj,
    p26: cobj,
    p27: cobj,
    p28: cobj,
    p29: cobj,
    p30: cobj,
    p31: cobj,
    p32: cobj,
    p33: cobj,
    p34: cobj,
):
    if stage == 1:
        macrop_driver_ptend_lq_mask_shell_codon(
            mode,
            pcnst,
            wtrc_nwset,
            wtrc_ncnst,
            use_water_tracers,
            ixcldliq,
            ixcldice,
            ixnumliq,
            ixnumice,
            p1,
            p2,
            p3,
            p4,
        )
    elif stage == 2:
        macrop_driver_detrain_init_lq_mask_shell_codon(
            ncol,
            pcols,
            pver,
            pcnst,
            wtrc_nwset,
            wtrc_ncnst,
            use_water_tracers,
            ixcldliq,
            ixcldice,
            ixnumliq,
            ixnumice,
            p1,
            p2,
            p3,
            p4,
            p5,
            p6,
            p7,
            p8,
            p9,
            p10,
            p11,
            p12,
            p13,
            p14,
            p15,
            p16,
            p17,
            p18,
        )
    elif stage == 3:
        macrop_driver_detrain_core_codon(
            ncol,
            pcols,
            pver,
            top_lev,
            do_detrain,
            cu_det_st,
            cpair,
            gravit,
            latice,
            nl_denom_a,
            nl_denom_b,
            ni_denom_a,
            ni_denom_b,
            p1,
            p2,
            p3,
            p4,
            p5,
            p6,
            p7,
            p8,
            p9,
            p10,
            p11,
            p12,
            p13,
            p14,
            p15,
            p16,
            p17,
            p18,
            p19,
            p20,
            p21,
            p22,
            p23,
        )
    elif stage == 4:
        macrop_driver_detrain_post_shell_codon(ncol, pcols, p1)
    elif stage == 5:
        macrop_driver_mmacro_prepare_shell_codon(
            ncol,
            pcols,
            pver,
            pcnst,
            top_lev,
            ixcldliq,
            ixcldice,
            ixnumliq,
            ixnumice,
            nstep,
            rdtime,
            p1,
            p2,
            p3,
            p4,
            p5,
            p6,
            p7,
            p8,
            p9,
            p10,
            p11,
            p12,
            p13,
            p14,
            p15,
            p16,
            p17,
            p18,
            p19,
            p20,
            p21,
            p22,
            p23,
            p24,
            p25,
            p26,
            p27,
            p28,
            p29,
            p30,
            p31,
            p32,
            p33,
            p34,
        )
    elif stage == 6:
        macrop_driver_mmacro_post_fields_shell_codon(ncol, pcols, pver, top_lev, p1, p2, p3, p4, p5)
    elif stage == 7:
        macrop_driver_ptend_config_shell_codon(
            ncol,
            pcols,
            pver,
            top_lev,
            do_cldice,
            do_cldliq,
            p1,
            p2,
            p3,
            p4,
            p5,
            p6,
            p7,
            p8,
            p9,
            p10,
            p11,
            p12,
            p13,
        )
    elif stage == 8:
        macrop_driver_cfmip_diag_shell_codon(ncol, pcols, pver, top_lev, p1, p2, p3, p4, p5, p6, p7)
    elif stage == 9:
        macrop_driver_clr_old_diag_codon(ncol, pcols, pver, top_lev, p1, p2, p3, p4, p5, p6)
    elif stage == 10:
        macrop_driver_store_state_codon(
            ncol,
            pcols,
            pver,
            top_lev,
            tmelt,
            p1,
            p2,
            p3,
            p4,
            p5,
            p6,
            p7,
            p8,
            p9,
            p10,
            p11,
            p12,
            p13,
        )
