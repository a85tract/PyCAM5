from C import bolton_svp_water_native_cb(float) -> float
from C import goffgratch_svp_water_native_cb(float) -> float
from C import murphykoop_svp_water_native_cb(float) -> float
from C import oldgoffgratch_svp_water_native_cb(float) -> float
from C import zm_entropy_expr_native_cb(float, float, float, float, float, float, float, float, float, float, float, float) -> float
from math import log, sqrt
import zm_cldprp_codon as _zm_cldprp


WV_SAT_OLD_GOFF_GRATCH_IDX = 0
WV_SAT_GOFF_GRATCH_IDX = 1
WV_SAT_MURPHY_KOOP_IDX = 2
WV_SAT_BOLTON_IDX = 3


@inline
def _zm_wv_sat_svp_to_qsat_codon(es: float, p: float, epsilo: float, omeps: float) -> float:
    if (p - es) <= 0.0:
        return 1.0
    return epsilo * es / (p - omeps * es)


@inline
def _zm_wv_sat_svp_water_codon(t: float, idx: int) -> float:
    if idx == WV_SAT_GOFF_GRATCH_IDX:
        return goffgratch_svp_water_native_cb(t)
    if idx == WV_SAT_MURPHY_KOOP_IDX:
        return murphykoop_svp_water_native_cb(t)
    if idx == WV_SAT_OLD_GOFF_GRATCH_IDX:
        return oldgoffgratch_svp_water_native_cb(t)
    if idx == WV_SAT_BOLTON_IDX:
        return bolton_svp_water_native_cb(t)
    return 0.0


@inline
def _zm_qsat_hpa_ptr_codon(
    t: float,
    p_hpa: float,
    idx: int,
    epsilo: float,
    omeps: float,
    es_out: Ptr[float],
    qm_out: Ptr[float],
):
    p_pa = p_hpa * 100.0
    es = _zm_wv_sat_svp_water_codon(t, idx)
    qm = _zm_wv_sat_svp_to_qsat_codon(es, p_pa, epsilo, omeps)
    if p_pa < es:
        es = p_pa

    es_out[0] = es * 0.01
    qm_out[0] = qm


@inline
def _zm_entropy_codon(
    tk: float,
    p_hpa: float,
    qtot: float,
    rl: float,
    cpliq: float,
    cpwv: float,
    tfreez: float,
    cpres: float,
    rgas: float,
    eps1: float,
    rh2o: float,
    idx: int,
    epsilo: float,
    omeps: float,
) -> float:
    est = 0.0
    qst = 0.0
    _zm_qsat_hpa_ptr_codon(tk, p_hpa, idx, epsilo, omeps, __ptr__(est), __ptr__(qst))
    return zm_entropy_expr_native_cb(
        tk,
        p_hpa,
        qtot,
        qst,
        rl,
        cpliq,
        cpwv,
        tfreez,
        cpres,
        rgas,
        eps1,
        rh2o,
    )


@inline
def _zm_ientropy_codon(
    s: float,
    p_hpa: float,
    qt: float,
    tfg: float,
    rl: float,
    cpliq: float,
    cpwv: float,
    tfreez: float,
    cpres: float,
    rgas: float,
    eps1: float,
    rh2o: float,
    idx: int,
    epsilo: float,
    omeps: float,
    t_out: Ptr[float],
    qst_out: Ptr[float],
    converged_out: Ptr[int],
):
    loopmax = 100
    eps = 3.0e-8
    tol = 0.001

    a = tfg - 10.0
    b = tfg + 10.0
    fa = _zm_entropy_codon(a, p_hpa, qt, rl, cpliq, cpwv, tfreez, cpres, rgas, eps1, rh2o, idx, epsilo, omeps) - s
    fb = _zm_entropy_codon(b, p_hpa, qt, rl, cpliq, cpwv, tfreez, cpres, rgas, eps1, rh2o, idx, epsilo, omeps) - s
    c = b
    fc = fb
    d = 0.0
    ebr = 0.0
    converged = False

    ii = 0
    while ii <= loopmax:
        if (fb > 0.0 and fc > 0.0) or (fb < 0.0 and fc < 0.0):
            c = a
            fc = fa
            d = b - a
            ebr = d

        if abs(fc) < abs(fb):
            a = b
            b = c
            c = a
            fa = fb
            fb = fc
            fc = fa

        tol1 = 2.0 * eps * abs(b) + 0.5 * tol
        xm = 0.5 * (c - b)
        converged = abs(xm) <= tol1 or fb == 0.0
        if converged:
            break

        if abs(ebr) >= tol1 and abs(fa) > abs(fb):
            sbr = fb / fa
            if a == c:
                pbr = 2.0 * xm * sbr
                qbr = 1.0 - sbr
            else:
                qbr = fa / fc
                rbr = fb / fc
                pbr = sbr * (2.0 * xm * qbr * (qbr - rbr) - (b - a) * (rbr - 1.0))
                qbr = (qbr - 1.0) * (rbr - 1.0) * (sbr - 1.0)

            if pbr > 0.0:
                qbr = -qbr
            pbr = abs(pbr)
            if 2.0 * pbr < min(3.0 * xm * qbr - abs(tol1 * qbr), abs(ebr * qbr)):
                ebr = d
                d = pbr / qbr
            else:
                d = xm
                ebr = d
        else:
            d = xm
            ebr = d

        a = b
        fa = fb
        if abs(d) > tol1:
            b = b + d
        else:
            if xm >= 0.0:
                b = b + abs(tol1)
            else:
                b = b - abs(tol1)

        fb = _zm_entropy_codon(b, p_hpa, qt, rl, cpliq, cpwv, tfreez, cpres, rgas, eps1, rh2o, idx, epsilo, omeps) - s
        ii += 1

    est = 0.0
    qst = 0.0
    _zm_qsat_hpa_ptr_codon(b, p_hpa, idx, epsilo, omeps, __ptr__(est), __ptr__(qst))
    t_out[0] = b
    qst_out[0] = qst
    if converged:
        converged_out[0] = 1
    else:
        converged_out[0] = 0


@export
def convect_deep_register_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def convect_deep_init_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def zm_conv_register_codon(zmconv_org: int) -> int:
    if zmconv_org != 0:
        return 1
    return 0


@export
def zm_conv_readnl_codon(
    c0_lnd_set: int,
    c0_ocn_set: int,
    ke_set: int,
    ke_lnd_set: int,
    org_set: int,
    flags_p: cobj,
):
    flags = Ptr[int](flags_p)
    mask = 0
    if c0_lnd_set != 0:
        mask |= 1
    if c0_ocn_set != 0:
        mask |= 2
    if ke_set != 0:
        mask |= 4
    if ke_lnd_set != 0:
        mask |= 8
    if org_set != 0:
        mask |= 16
    flags[0] = mask


@export
def zm_conv_init_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def convect_deep_tend_2_action_codon(scheme_code: int) -> int:
    if scheme_code == 1:
        return 1
    return 0


@export
def convect_deep_tend_action_codon(scheme_code: int) -> int:
    if scheme_code == 1:
        return 1
    if scheme_code == 2 or scheme_code == 3 or scheme_code == 4:
        return 2
    return 0


@export
def zm_conv_init_limcnv_codon(plev: int, pref_edge_p: cobj) -> int:
    pref_edge = Ptr[float](pref_edge_p)
    threshold = 4000.0

    if pref_edge[0] >= threshold:
        return 1

    for k in range(1, plev + 1):
        if pref_edge[k - 1] < threshold and pref_edge[k] >= threshold:
            return k

    return plev + 1


@export
def convect_deep_select_scheme_codon(
    scheme_len: int,
    scheme_ascii_p: cobj,
    scheme_code_p: cobj,
    status_p: cobj,
):
    scheme_ascii = Ptr[int](scheme_ascii_p)
    scheme_code = Ptr[int](scheme_code_p)
    status = Ptr[int](status_p)

    status[0] = 0
    scheme_code[0] = 0

    n = scheme_len
    while n > 0 and scheme_ascii[n - 1] == 32:
        n -= 1

    if n == 2:
        c1 = scheme_ascii[0]
        c2 = scheme_ascii[1]

        if c1 >= 65 and c1 <= 90:
            c1 += 32
        if c2 >= 65 and c2 <= 90:
            c2 += 32

        if c1 == 122 and c2 == 109:
            scheme_code[0] = 1
            return

        if c1 == 107 and c2 == 101:
            scheme_code[0] = 5
            return

    if n == 3:
        c1 = scheme_ascii[0]
        c2 = scheme_ascii[1]
        c3 = scheme_ascii[2]

        if c1 >= 65 and c1 <= 90:
            c1 += 32
        if c2 >= 65 and c2 <= 90:
            c2 += 32
        if c3 >= 65 and c3 <= 90:
            c3 += 32

        if c1 == 111 and c2 == 102 and c3 == 102:
            scheme_code[0] = 2
            return

    if n == 6:
        c1 = scheme_ascii[0]
        c2 = scheme_ascii[1]
        c3 = scheme_ascii[2]
        c4 = scheme_ascii[3]
        c5 = scheme_ascii[4]
        c6 = scheme_ascii[5]

        if c1 >= 65 and c1 <= 90:
            c1 += 32
        if c2 >= 65 and c2 <= 90:
            c2 += 32
        if c3 >= 65 and c3 <= 90:
            c3 += 32
        if c4 >= 65 and c4 <= 90:
            c4 += 32
        if c5 >= 65 and c5 <= 90:
            c5 += 32
        if c6 >= 65 and c6 <= 90:
            c6 += 32

        if c1 == 117 and c2 == 110 and c3 == 105 and c4 == 99 and c5 == 111 and c6 == 110:
            scheme_code[0] = 3
            return

    if n == 9:
        c1 = scheme_ascii[0]
        c2 = scheme_ascii[1]
        c3 = scheme_ascii[2]
        c4 = scheme_ascii[3]
        c5 = scheme_ascii[4]
        c6 = scheme_ascii[5]
        c7 = scheme_ascii[6]
        c8 = scheme_ascii[7]
        c9 = scheme_ascii[8]

        if c1 >= 65 and c1 <= 90:
            c1 += 32
        if c2 >= 65 and c2 <= 90:
            c2 += 32
        if c3 >= 65 and c3 <= 90:
            c3 += 32
        if c4 >= 65 and c4 <= 90:
            c4 += 32
        if c5 >= 65 and c5 <= 90:
            c5 += 32
        if c6 >= 65 and c6 <= 90:
            c6 += 32
        if c7 >= 65 and c7 <= 90:
            c7 += 32
        if c8 >= 65 and c8 <= 90:
            c8 += 32
        if c9 >= 65 and c9 <= 90:
            c9 += 32

        if (
            c1 == 99 and c2 == 108 and c3 == 117 and c4 == 98 and c5 == 98
            and c6 == 95 and c7 == 115 and c8 == 103 and c9 == 115
        ):
            scheme_code[0] = 4
            return

    status[0] = 1


@export
def deep_scheme_does_scav_trans_codon(scheme_code: int) -> int:
    if scheme_code == 5:
        return 1
    return 0


@export
def zm_cldprp_init_arrays_codon(
    il2g: int,
    pcols: int,
    pver: int,
    c0_ocn: float,
    c0_lnd: float,
    landfrac_p: cobj,
    zf_p: cobj,
    q_p: cobj,
    s_p: cobj,
    ftemp_p: cobj,
    expnum_p: cobj,
    expdif_p: cobj,
    c0mask_p: cobj,
    dz_p: cobj,
    pflx_p: cobj,
    k1_p: cobj,
    i2_p: cobj,
    i3_p: cobj,
    i4_p: cobj,
    mu_p: cobj,
    f_p: cobj,
    eps_p: cobj,
    eu_p: cobj,
    du_p: cobj,
    ql_p: cobj,
    cu_p: cobj,
    evp_p: cobj,
    wtevp_p: cobj,
    cmeg_p: cobj,
    qds_p: cobj,
    md_p: cobj,
    ed_p: cobj,
    sd_p: cobj,
    qd_p: cobj,
    mc_p: cobj,
    qu_p: cobj,
    su_p: cobj,
    rprd_p: cobj,
    totpcp_p: cobj,
    totevp_p: cobj,
):
    _zm_cldprp.zm_cldprp_init_arrays_codon(
        il2g, pcols, pver, c0_ocn, c0_lnd, landfrac_p, zf_p, q_p, s_p, ftemp_p, expnum_p,
        expdif_p, c0mask_p, dz_p, pflx_p, k1_p, i2_p, i3_p, i4_p, mu_p, f_p, eps_p, eu_p,
        du_p, ql_p, cu_p, evp_p, wtevp_p, cmeg_p, qds_p, md_p, ed_p, sd_p, qd_p, mc_p, qu_p,
        su_p, rprd_p, totpcp_p, totevp_p,
    )


@export
def zm_cldprp_thermo_level_codon(
    il2g: int,
    pcols: int,
    msg: int,
    k: int,
    eps1: float,
    rl: float,
    rd: float,
    cp: float,
    grav: float,
    t_p: cobj,
    p_p: cobj,
    z_p: cobj,
    zf_p: cobj,
    q_p: cobj,
    qst_p: cobj,
    est_p: cobj,
    gamma_p: cobj,
    hmn_p: cobj,
    hsat_p: cobj,
    hu_p: cobj,
    hd_p: cobj,
    sd_p: cobj,
    su_p: cobj,
    tdt_p: cobj,
    tut_p: cobj,
    rprd_p: cobj,
    hsthat_p: cobj,
    qsthat_p: cobj,
    gamhat_p: cobj,
):
    _zm_cldprp.zm_cldprp_thermo_level_codon(
        il2g,
        pcols,
        msg,
        k,
        eps1,
        rl,
        rd,
        cp,
        grav,
        t_p,
        p_p,
        z_p,
        zf_p,
        q_p,
        qst_p,
        est_p,
        gamma_p,
        hmn_p,
        hsat_p,
        hu_p,
        hd_p,
        sd_p,
        su_p,
        tdt_p,
        tut_p,
        rprd_p,
        hsthat_p,
        qsthat_p,
        gamhat_p,
    )


@export
def zm_cldprp_interface_interp_codon(
    il2g: int,
    pcols: int,
    pver: int,
    msg: int,
    cp: float,
    rl: float,
    qst_p: cobj,
    gamma_p: cobj,
    shat_p: cobj,
    qsthat_p: cobj,
    hsthat_p: cobj,
    gamhat_p: cobj,
    totpcp_p: cobj,
    totevp_p: cobj,
):
    _zm_cldprp.zm_cldprp_interface_interp_codon(
        il2g, pcols, pver, msg, cp, rl, qst_p, gamma_p, shat_p, qsthat_p, hsthat_p,
        gamhat_p, totpcp_p, totevp_p,
    )


@export
def zm_cldprp_index_setup_codon(
    il2g: int,
    pcols: int,
    pver: int,
    msg: int,
    limcnv: int,
    cp_tiedke_add: float,
    tiedke_add: float,
    jb_p: cobj,
    lel_p: cobj,
    mx_p: cobj,
    hsat_p: cobj,
    hmn_p: cobj,
    s_p: cobj,
    jt_p: cobj,
    jd_p: cobj,
    jlcl_p: cobj,
    hmin_p: cobj,
    j0_p: cobj,
    hu_p: cobj,
    su_p: cobj,
):
    _zm_cldprp.zm_cldprp_index_setup_codon(
        il2g, pcols, pver, msg, limcnv, cp_tiedke_add, tiedke_add, jb_p, lel_p, mx_p, hsat_p, hmn_p, s_p,
        jt_p, jd_p, jlcl_p, hmin_p, j0_p, hu_p, su_p,
    )


@export
def zm_cldprp_copy_mass_fields_codon(
    pcols: int,
    pver: int,
    ed_p: cobj,
    md_p: cobj,
    mu_p: cobj,
    du_p: cobj,
    eu_p: cobj,
    wted_p: cobj,
    wtmd_p: cobj,
    wtmu_p: cobj,
    wtdu_p: cobj,
    wteu_p: cobj,
):
    _zm_cldprp.zm_cldprp_copy_mass_fields_codon(
        pcols, pver, ed_p, md_p, mu_p, du_p, eu_p, wted_p, wtmd_p, wtmu_p, wtdu_p, wteu_p,
    )


@export
def zm_cldprp_copy_2d_codon(
    pcols: int,
    pver: int,
    src_p: cobj,
    dst_p: cobj,
):
    _zm_cldprp.zm_cldprp_copy_2d_codon(pcols, pver, src_p, dst_p)


@export
def zm_cldprp_taylor_hmin_codon(
    il2g: int,
    pcols: int,
    pver: int,
    msg: int,
    jt_p: cobj,
    jb_p: cobj,
    mx_p: cobj,
    j0_p: cobj,
    hmn_p: cobj,
    dz_p: cobj,
    k1_p: cobj,
    ihat_p: cobj,
    i2_p: cobj,
    idag_p: cobj,
    i3_p: cobj,
    iprm_p: cobj,
    i4_p: cobj,
    hmin_p: cobj,
    expdif_p: cobj,
):
    _zm_cldprp.zm_cldprp_taylor_hmin_codon(
        il2g, pcols, pver, msg, jt_p, jb_p, mx_p, j0_p, hmn_p, dz_p, k1_p, ihat_p,
        i2_p, idag_p, i3_p, iprm_p, i4_p, hmin_p, expdif_p,
    )


@export
def zm_cldprp_taylor_f_codon(
    il2g: int,
    pcols: int,
    pver: int,
    msg: int,
    jt_p: cobj,
    jb_p: cobj,
    mx_p: cobj,
    expdif_p: cobj,
    expnum_p: cobj,
    ftemp_p: cobj,
    hmn_p: cobj,
    hsat_p: cobj,
    zf_p: cobj,
    z_p: cobj,
    k1_p: cobj,
    i2_p: cobj,
    i3_p: cobj,
    i4_p: cobj,
    dz_p: cobj,
    f_p: cobj,
):
    _zm_cldprp.zm_cldprp_taylor_f_codon(
        il2g, pcols, pver, msg, jt_p, jb_p, mx_p, expdif_p, expnum_p, ftemp_p, hmn_p,
        hsat_p, zf_p, z_p, k1_p, i2_p, i3_p, i4_p, dz_p, f_p,
    )


@export
def zm_cldprp_eps_profile_codon(
    il2g: int,
    pcols: int,
    pver: int,
    msg: int,
    jt_p: cobj,
    jb_p: cobj,
    j0_p: cobj,
    f_p: cobj,
    eps_p: cobj,
    eps0_p: cobj,
):
    _zm_cldprp.zm_cldprp_eps_profile_codon(il2g, pcols, pver, msg, jt_p, jb_p, j0_p, f_p, eps_p, eps0_p)


@export
def zm_cldprp_updraft_mass_energy_codon(
    il2g: int,
    pcols: int,
    pver: int,
    pverp: int,
    msg: int,
    jb_p: cobj,
    jt_p: cobj,
    lel_p: cobj,
    eps0_p: cobj,
    eps_p: cobj,
    zf_p: cobj,
    dz_p: cobj,
    mu_p: cobj,
    eu_p: cobj,
    du_p: cobj,
    hmn_p: cobj,
    hsat_p: cobj,
    hu_p: cobj,
):
    _zm_cldprp.zm_cldprp_updraft_mass_energy_codon(
        il2g, pcols, pver, pverp, msg, jb_p, jt_p, lel_p, eps0_p, eps_p, zf_p, dz_p,
        mu_p, eu_p, du_p, hmn_p, hsat_p, hu_p,
    )


@export
def zm_cldprp_cloud_top_reset_codon(
    il2g: int,
    pcols: int,
    pver: int,
    pverp: int,
    msg: int,
    lel_p: cobj,
    jb_p: cobj,
    jt_p: cobj,
    eps0_p: cobj,
    mu_p: cobj,
    eu_p: cobj,
    du_p: cobj,
    hu_p: cobj,
    hmn_p: cobj,
    hsthat_p: cobj,
    dz_p: cobj,
):
    _zm_cldprp.zm_cldprp_cloud_top_reset_codon(
        il2g, pcols, pver, pverp, msg, lel_p, jb_p, jt_p, eps0_p, mu_p, eu_p, du_p,
        hu_p, hmn_p, hsthat_p, dz_p,
    )


@export
def zm_cldprp_downdraft_init_codon(
    il2g: int,
    pcols: int,
    jt_p: cobj,
    jb_p: cobj,
    j0_p: cobj,
    jd_p: cobj,
    hmn_p: cobj,
    hd_p: cobj,
    eps0_p: cobj,
    epsm_p: cobj,
    alfa_p: cobj,
    md_p: cobj,
):
    _zm_cldprp.zm_cldprp_downdraft_init_codon(
        il2g, pcols, jt_p, jb_p, j0_p, jd_p, hmn_p, hd_p, eps0_p, epsm_p, alfa_p, md_p,
    )


@export
def zm_cldprp_downdraft_mass_profile_codon(
    il2g: int,
    pcols: int,
    pver: int,
    msg: int,
    jd_p: cobj,
    jb_p: cobj,
    eps0_p: cobj,
    epsm_p: cobj,
    alfa_p: cobj,
    zf_p: cobj,
    md_p: cobj,
):
    _zm_cldprp.zm_cldprp_downdraft_mass_profile_codon(
        il2g, pcols, pver, msg, jd_p, jb_p, eps0_p, epsm_p, alfa_p, zf_p, md_p,
    )


@export
def zm_cldprp_downdraft_scale_energy_codon(
    il2g: int,
    pcols: int,
    pver: int,
    msg: int,
    small: float,
    jt_p: cobj,
    jb_p: cobj,
    jd_p: cobj,
    eps0_p: cobj,
    mu_p: cobj,
    md_p: cobj,
    dz_p: cobj,
    ed_p: cobj,
    hd_p: cobj,
    hmn_p: cobj,
):
    _zm_cldprp.zm_cldprp_downdraft_scale_energy_codon(
        il2g, pcols, pver, msg, small, jt_p, jb_p, jd_p, eps0_p, mu_p, md_p, dz_p, ed_p,
        hd_p, hmn_p,
    )


@export
def zm_cldprp_qds_codon(
    il2g: int,
    pcols: int,
    pver: int,
    msg: int,
    rl: float,
    jd_p: cobj,
    jb_p: cobj,
    eps0_p: cobj,
    qds_p: cobj,
    qsthat_p: cobj,
    gamhat_p: cobj,
    hd_p: cobj,
    hsthat_p: cobj,
):
    _zm_cldprp.zm_cldprp_qds_codon(
        il2g, pcols, pver, msg, rl, jd_p, jb_p, eps0_p, qds_p, qsthat_p, gamhat_p, hd_p, hsthat_p,
    )


@export
def zm_cldprp_updraft_seed_codon(
    il2g: int,
    pcols: int,
    rl: float,
    cp: float,
    jb_p: cobj,
    mx_p: cobj,
    eps0_p: cobj,
    q_p: cobj,
    hu_p: cobj,
    qu_p: cobj,
    su_p: cobj,
):
    _zm_cldprp.zm_cldprp_updraft_seed_codon(
        il2g, pcols, rl, cp, jb_p, mx_p, eps0_p, q_p, hu_p, qu_p, su_p,
    )


@export
def zm_cldprp_updraft_lcl_reset_codon(
    il2g: int,
    done_p: cobj,
    active_p: cobj,
    found_p: cobj,
    tu_p: cobj,
    qstu_p: cobj,
    kount_p: cobj,
):
    _zm_cldprp.zm_cldprp_updraft_lcl_reset_codon(
        il2g, done_p, active_p, found_p, tu_p, qstu_p, kount_p,
    )


@export
def zm_cldprp_updraft_lcl_prepare_codon(
    il2g: int,
    pcols: int,
    k: int,
    cp: float,
    grav: float,
    jt_p: cobj,
    jb_p: cobj,
    done_p: cobj,
    eps0_p: cobj,
    mu_p: cobj,
    dz_p: cobj,
    eu_p: cobj,
    du_p: cobj,
    s_p: cobj,
    q_p: cobj,
    qst_p: cobj,
    zf_p: cobj,
    su_p: cobj,
    qu_p: cobj,
    active_p: cobj,
    found_p: cobj,
    tu_p: cobj,
):
    _zm_cldprp.zm_cldprp_updraft_lcl_prepare_codon(
        il2g, pcols, k, cp, grav, jt_p, jb_p, done_p, eps0_p, mu_p, dz_p, eu_p,
        du_p, s_p, q_p, qst_p, zf_p, su_p, qu_p, active_p, found_p, tu_p,
    )


@export
def zm_cldprp_updraft_lcl_finalize_codon(
    il2g: int,
    pcols: int,
    k: int,
    active_p: cobj,
    qstu_p: cobj,
    tu_p: cobj,
    qu_p: cobj,
    tut_p: cobj,
    jlcl_p: cobj,
    done_p: cobj,
    found_p: cobj,
    kount_p: cobj,
):
    _zm_cldprp.zm_cldprp_updraft_lcl_finalize_codon(
        il2g, pcols, k, active_p, qstu_p, tu_p, qu_p, tut_p, jlcl_p, done_p,
        found_p, kount_p,
    )


@export
def zm_cldprp_updraft_saturation_adjust_codon(
    il2g: int,
    pcols: int,
    pver: int,
    msg: int,
    cp: float,
    grav: float,
    rl: float,
    jt_p: cobj,
    jlcl_p: cobj,
    eps0_p: cobj,
    shat_p: cobj,
    hu_p: cobj,
    hsthat_p: cobj,
    gamhat_p: cobj,
    zf_p: cobj,
    qsthat_p: cobj,
    su_p: cobj,
    tut_p: cobj,
    qu_p: cobj,
):
    _zm_cldprp.zm_cldprp_updraft_saturation_adjust_codon(
        il2g, pcols, pver, msg, cp, grav, rl, jt_p, jlcl_p, eps0_p, shat_p, hu_p,
        hsthat_p, gamhat_p, zf_p, qsthat_p, su_p, tut_p, qu_p,
    )


@export
def zm_cldprp_condensation_codon(
    il2g: int,
    pcols: int,
    pver: int,
    msg: int,
    rl: float,
    cp: float,
    jt_p: cobj,
    jb_p: cobj,
    eps0_p: cobj,
    mu_p: cobj,
    su_p: cobj,
    dz_p: cobj,
    eu_p: cobj,
    du_p: cobj,
    s_p: cobj,
    cu_p: cobj,
):
    _zm_cldprp.zm_cldprp_condensation_codon(
        il2g, pcols, pver, msg, rl, cp, jt_p, jb_p, eps0_p, mu_p, su_p, dz_p, eu_p, du_p, s_p, cu_p,
    )


@export
def zm_cldprp_rain_production_codon(
    il2g: int,
    pcols: int,
    pver: int,
    msg: int,
    jt_p: cobj,
    jb_p: cobj,
    eps0_p: cobj,
    mu_p: cobj,
    ql_p: cobj,
    dz_p: cobj,
    du_p: cobj,
    cu_p: cobj,
    c0mask_p: cobj,
    totpcp_p: cobj,
    rprd_p: cobj,
):
    _zm_cldprp.zm_cldprp_rain_production_codon(
        il2g, pcols, pver, msg, jt_p, jb_p, eps0_p, mu_p, ql_p, dz_p, du_p, cu_p,
        c0mask_p, totpcp_p, rprd_p,
    )


@export
def zm_cldprp_downdraft_seed_codon(
    il2g: int,
    pcols: int,
    cp: float,
    grav: float,
    rl: float,
    jd_p: cobj,
    qds_p: cobj,
    qd_p: cobj,
    hd_p: cobj,
    sd_p: cobj,
    tdt_p: cobj,
    zf_p: cobj,
):
    _zm_cldprp.zm_cldprp_downdraft_seed_codon(
        il2g, pcols, cp, grav, rl, jd_p, qds_p, qd_p, hd_p, sd_p, tdt_p, zf_p,
    )


@export
def zm_cldprp_downdraft_evap_codon(
    il2g: int,
    pcols: int,
    pver: int,
    msg: int,
    cp: float,
    grav: float,
    rl: float,
    small: float,
    jd_p: cobj,
    jb_p: cobj,
    eps0_p: cobj,
    q_p: cobj,
    s_p: cobj,
    zf_p: cobj,
    dz_p: cobj,
    ed_p: cobj,
    md_p: cobj,
    qd_p: cobj,
    qds_p: cobj,
    sd_p: cobj,
    tdt_p: cobj,
    evp_p: cobj,
    totevp_p: cobj,
):
    _zm_cldprp.zm_cldprp_downdraft_evap_codon(
        il2g, pcols, pver, msg, cp, grav, rl, small, jd_p, jb_p, eps0_p, q_p, s_p,
        zf_p, dz_p, ed_p, md_p, qd_p, qds_p, sd_p, tdt_p, evp_p, totevp_p,
    )


@export
def zm_cldprp_evap_finalize_codon(
    il2g: int,
    pcols: int,
    pver: int,
    pverp: int,
    msg: int,
    jd_p: cobj,
    jb_p: cobj,
    md_p: cobj,
    qd_p: cobj,
    totevp_p: cobj,
    totpcp_p: cobj,
    ed_p: cobj,
    evp_p: cobj,
    cu_p: cobj,
    cmeg_p: cobj,
    rprd_p: cobj,
    dz_p: cobj,
    pflx_p: cobj,
    mc_p: cobj,
    mu_p: cobj,
    wtevp_p: cobj,
):
    _zm_cldprp.zm_cldprp_evap_finalize_codon(
        il2g, pcols, pver, pverp, msg, jd_p, jb_p, md_p, qd_p, totevp_p, totpcp_p, ed_p,
        evp_p, cu_p, cmeg_p, rprd_p, dz_p, pflx_p, mc_p, mu_p, wtevp_p,
    )


@export
def zm_cldprp_pre_lcl_batch_codon(
    il2g: int,
    pcols: int,
    pver: int,
    pverp: int,
    msg: int,
    limcnv: int,
    cp_tiedke_add: float,
    tiedke_add: float,
    small: float,
    rl: float,
    cp: float,
    jb_p: cobj,
    lel_p: cobj,
    mx_p: cobj,
    hsat_p: cobj,
    hmn_p: cobj,
    s_p: cobj,
    jt_p: cobj,
    jd_p: cobj,
    jlcl_p: cobj,
    hmin_p: cobj,
    j0_p: cobj,
    hu_p: cobj,
    su_p: cobj,
    dz_p: cobj,
    k1_p: cobj,
    ihat_p: cobj,
    i2_p: cobj,
    idag_p: cobj,
    i3_p: cobj,
    iprm_p: cobj,
    i4_p: cobj,
    expdif_p: cobj,
    expnum_p: cobj,
    ftemp_p: cobj,
    zf_p: cobj,
    z_p: cobj,
    f_p: cobj,
    eps_p: cobj,
    eps0_p: cobj,
    eu_p: cobj,
    du_p: cobj,
    mu_p: cobj,
    hd_p: cobj,
    epsm_p: cobj,
    alfa_p: cobj,
    md_p: cobj,
    ed_p: cobj,
    wted_p: cobj,
    wtmd_p: cobj,
    wtmu_p: cobj,
    wtdu_p: cobj,
    wteu_p: cobj,
    qds_p: cobj,
    qsthat_p: cobj,
    gamhat_p: cobj,
    hsthat_p: cobj,
    q_p: cobj,
    qu_p: cobj,
    lcl_done_p: cobj,
    lcl_active_p: cobj,
    lcl_found_p: cobj,
    lcl_tu_p: cobj,
    lcl_qstu_p: cobj,
    lcl_kount_p: cobj,
):
    _zm_cldprp.zm_cldprp_pre_lcl_batch_codon(
        il2g, pcols, pver, pverp, msg, limcnv, cp_tiedke_add, tiedke_add, small,
        rl, cp, jb_p, lel_p, mx_p, hsat_p, hmn_p, s_p, jt_p, jd_p, jlcl_p,
        hmin_p, j0_p, hu_p, su_p, dz_p, k1_p, ihat_p, i2_p, idag_p, i3_p,
        iprm_p, i4_p, expdif_p, expnum_p, ftemp_p, zf_p, z_p, f_p, eps_p,
        eps0_p, eu_p, du_p, mu_p, hd_p, epsm_p, alfa_p, md_p, ed_p, wted_p,
        wtmd_p, wtmu_p, wtdu_p, wteu_p, qds_p, qsthat_p, gamhat_p, hsthat_p,
        q_p, qu_p, lcl_done_p, lcl_active_p, lcl_found_p, lcl_tu_p, lcl_qstu_p,
        lcl_kount_p,
    )


@export
def zm_cldprp_post_lcl_batch_codon(
    il2g: int,
    pcols: int,
    pver: int,
    pverp: int,
    msg: int,
    cp: float,
    grav: float,
    rl: float,
    small: float,
    jt_p: cobj,
    jlcl_p: cobj,
    eps0_p: cobj,
    shat_p: cobj,
    hu_p: cobj,
    hsthat_p: cobj,
    gamhat_p: cobj,
    zf_p: cobj,
    qsthat_p: cobj,
    su_p: cobj,
    tut_p: cobj,
    qu_p: cobj,
    jb_p: cobj,
    mu_p: cobj,
    dz_p: cobj,
    eu_p: cobj,
    du_p: cobj,
    s_p: cobj,
    cu_p: cobj,
    ql_p: cobj,
    c0mask_p: cobj,
    totpcp_p: cobj,
    rprd_p: cobj,
    rppe_p: cobj,
    jd_p: cobj,
    qds_p: cobj,
    qd_p: cobj,
    hd_p: cobj,
    sd_p: cobj,
    tdt_p: cobj,
    q_p: cobj,
    ed_p: cobj,
    md_p: cobj,
    evp_p: cobj,
    totevp_p: cobj,
    cmeg_p: cobj,
    pflx_p: cobj,
    mc_p: cobj,
    wtevp_p: cobj,
):
    _zm_cldprp.zm_cldprp_post_lcl_batch_codon(
        il2g, pcols, pver, pverp, msg, cp, grav, rl, small, jt_p, jlcl_p,
        eps0_p, shat_p, hu_p, hsthat_p, gamhat_p, zf_p, qsthat_p, su_p,
        tut_p, qu_p, jb_p, mu_p, dz_p, eu_p, du_p, s_p, cu_p, ql_p,
        c0mask_p, totpcp_p, rprd_p, rppe_p, jd_p, qds_p, qd_p, hd_p, sd_p,
        tdt_p, q_p, ed_p, md_p, evp_p, totevp_p, cmeg_p, pflx_p, mc_p,
        wtevp_p,
    )


@export
def zm_conv_ptend_lq_mask_shell_codon(
    pcnst: int,
    wtrc_nwset: int,
    org_enabled: int,
    ixorg: int,
    vap_type_p: cobj,
    lq_mask_p: cobj,
):
    vap_type = Ptr[int](vap_type_p)
    lq_mask = Ptr[int](lq_mask_p)

    for m in range(0, pcnst):
        lq_mask[m] = 0

    lq_mask[0] = 1

    if org_enabled != 0:
        lq_mask[ixorg - 1] = 1

    for m in range(0, wtrc_nwset):
        lq_mask[vap_type[m] - 1] = 1


@export
def zm_convr_post_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    lengath: int,
    gravit: float,
    cpair: float,
    mcon_p: cobj,
    mu_p: cobj,
    md_p: cobj,
    ideep_p: cobj,
    jt_p: cobj,
    maxg_p: cobj,
    ptend_s_p: cobj,
    state_ps_p: cobj,
    state_pmid_p: cobj,
    freqzm_p: cobj,
    mu_out_p: cobj,
    md_out_p: cobj,
    ftem_p: cobj,
    pcont_p: cobj,
    pconb_p: cobj,
):
    mcon = Ptr[float](mcon_p)
    mu = Ptr[float](mu_p)
    md = Ptr[float](md_p)
    ideep = Ptr[int](ideep_p)
    jt = Ptr[int](jt_p)
    maxg = Ptr[int](maxg_p)
    ptend_s = Ptr[float](ptend_s_p)
    state_ps = Ptr[float](state_ps_p)
    state_pmid = Ptr[float](state_pmid_p)
    freqzm = Ptr[float](freqzm_p)
    mu_out = Ptr[float](mu_out_p)
    md_out = Ptr[float](md_out_p)
    ftem = Ptr[float](ftem_p)
    pcont = Ptr[float](pcont_p)
    pconb = Ptr[float](pconb_p)

    i = 0
    while i < pcols:
        freqzm[i] = 0.0
        i += 1

    i = 0
    while i < lengath:
        freqzm[ideep[i] - 1] = 1.0
        i += 1

    k = 0
    while k < pver:
        i = 0
        while i < ncol:
            idx = i + k * pcols
            mcon[idx] = mcon[idx] * 100.0 / gravit
            ftem[idx] = ptend_s[idx] / cpair
            i += 1
        k += 1

    k = 0
    while k < pver:
        i = 0
        while i < pcols:
            idx = i + k * pcols
            mu_out[idx] = 0.0
            md_out[idx] = 0.0
            i += 1
        k += 1

    i = 0
    while i < lengath:
        ii = ideep[i] - 1
        k = 0
        while k < pver:
            gathered_idx = i + k * pcols
            out_idx = ii + k * pcols
            mu_out[out_idx] = mu[gathered_idx] * 100.0 / gravit
            md_out[out_idx] = md[gathered_idx] * 100.0 / gravit
            k += 1
        i += 1

    i = 0
    while i < ncol:
        pcont[i] = state_ps[i]
        pconb[i] = state_ps[i]
        i += 1

    i = 0
    while i < lengath:
        if maxg[i] > jt[i]:
            ii = ideep[i] - 1
            pcont[ii] = state_pmid[ii + (jt[i] - 1) * pcols]
            pconb[ii] = state_pmid[ii + (maxg[i] - 1) * pcols]
        i += 1


@export
def zm_convr_init_stage_dispatch_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    qtnd_p: cobj,
    heat_p: cobj,
    mcon_p: cobj,
    rliq_p: cobj,
    prec_p: cobj,
    dqdt_p: cobj,
    dsdt_p: cobj,
    dudt_p: cobj,
    dvdt_p: cobj,
    pflx_p: cobj,
    pflxg_p: cobj,
    cme_p: cobj,
    rprd_p: cobj,
    zdu_p: cobj,
    ql_p: cobj,
    qlg_p: cobj,
    dlf_p: cobj,
    dlg_p: cobj,
    tug_p: cobj,
    tdg_p: cobj,
    tu_p: cobj,
    td_p: cobj,
    cu_p: cobj,
    evp_p: cobj,
    wtcu_p: cobj,
    t_p: cobj,
    pblt_p: cobj,
    dsubcld_p: cobj,
    jctop_p: cobj,
    jcbot_p: cobj,
):
    qtnd = Ptr[float](qtnd_p)
    heat = Ptr[float](heat_p)
    mcon = Ptr[float](mcon_p)
    rliq = Ptr[float](rliq_p)
    prec = Ptr[float](prec_p)
    dqdt = Ptr[float](dqdt_p)
    dsdt = Ptr[float](dsdt_p)
    dudt = Ptr[float](dudt_p)
    dvdt = Ptr[float](dvdt_p)
    pflx = Ptr[float](pflx_p)
    pflxg = Ptr[float](pflxg_p)
    cme = Ptr[float](cme_p)
    rprd = Ptr[float](rprd_p)
    zdu = Ptr[float](zdu_p)
    ql = Ptr[float](ql_p)
    qlg = Ptr[float](qlg_p)
    dlf = Ptr[float](dlf_p)
    dlg = Ptr[float](dlg_p)
    tug = Ptr[float](tug_p)
    tdg = Ptr[float](tdg_p)
    tu = Ptr[float](tu_p)
    td = Ptr[float](td_p)
    cu = Ptr[float](cu_p)
    evp = Ptr[float](evp_p)
    wtcu = Ptr[float](wtcu_p)
    t = Ptr[float](t_p)
    pblt = Ptr[float](pblt_p)
    dsubcld = Ptr[float](dsubcld_p)
    jctop = Ptr[float](jctop_p)
    jcbot = Ptr[float](jcbot_p)

    k = 0
    while k < pver:
        i = 0
        while i < pcols:
            idx = i + k * pcols
            qtnd[idx] = 0.0
            heat[idx] = 0.0
            i += 1
        k += 1
    k = 0
    while k < pverp:
        i = 0
        while i < pcols:
            mcon[i + k * pcols] = 0.0
            i += 1
        k += 1
    i = 0
    while i < ncol:
        rliq[i] = 0.0
        prec[i] = 0.0
        i += 1

    k = 0
    while k < pver:
        i = 0
        while i < ncol:
            idx = i + k * pcols
            dqdt[idx] = 0.0
            dsdt[idx] = 0.0
            dudt[idx] = 0.0
            dvdt[idx] = 0.0
            pflx[idx] = 0.0
            pflxg[idx] = 0.0
            cme[idx] = 0.0
            rprd[idx] = 0.0
            zdu[idx] = 0.0
            ql[idx] = 0.0
            qlg[idx] = 0.0
            dlf[idx] = 0.0
            dlg[idx] = 0.0
            tug[idx] = 0.0
            tdg[idx] = 0.0
            tu[idx] = t[idx]
            td[idx] = t[idx]
            cu[idx] = 0.0
            evp[idx] = 0.0
            wtcu[idx] = 0.0
            i += 1
        k += 1

    i = 0
    while i < ncol:
        pflx[i + (pverp - 1) * pcols] = 0.0
        pflxg[i + (pverp - 1) * pcols] = 0.0
        pblt[i] = float(pver)
        dsubcld[i] = 0.0
        jctop[i] = float(pver)
        jcbot[i] = 1.0
        i += 1

@export
def zm_convr_init_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    qtnd_p: cobj,
    heat_p: cobj,
    mcon_p: cobj,
    rliq_p: cobj,
    prec_p: cobj,
    dqdt_p: cobj,
    dsdt_p: cobj,
    dudt_p: cobj,
    dvdt_p: cobj,
    pflx_p: cobj,
    pflxg_p: cobj,
    cme_p: cobj,
    rprd_p: cobj,
    zdu_p: cobj,
    ql_p: cobj,
    qlg_p: cobj,
    dlf_p: cobj,
    dlg_p: cobj,
    tug_p: cobj,
    tdg_p: cobj,
    tu_p: cobj,
    td_p: cobj,
    cu_p: cobj,
    evp_p: cobj,
    wtcu_p: cobj,
    t_p: cobj,
    pblt_p: cobj,
    dsubcld_p: cobj,
    jctop_p: cobj,
    jcbot_p: cobj,
):
    zm_convr_init_stage_dispatch_codon(
        ncol,
        pcols,
        pver,
        pverp,
        qtnd_p,
        heat_p,
        mcon_p,
        rliq_p,
        prec_p,
        dqdt_p,
        dsdt_p,
        dudt_p,
        dvdt_p,
        pflx_p,
        pflxg_p,
        cme_p,
        rprd_p,
        zdu_p,
        ql_p,
        qlg_p,
        dlf_p,
        dlg_p,
        tug_p,
        tdg_p,
        tu_p,
        td_p,
        cu_p,
        evp_p,
        wtcu_p,
        t_p,
        pblt_p,
        dsubcld_p,
        jctop_p,
        jcbot_p,
    )


@export
def zm_convr_pressure_state_stage_dispatch_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    msg: int,
    rgrav: float,
    grav: float,
    cpres: float,
    geos_p: cobj,
    zi_p: cobj,
    paph_p: cobj,
    pap_p: cobj,
    zm_p: cobj,
    pblh_p: cobj,
    t_p: cobj,
    qh_p: cobj,
    zs_p: cobj,
    pf_p: cobj,
    zf_p: cobj,
    p_p: cobj,
    z_p: cobj,
    q_p: cobj,
    s_p: cobj,
    tp_p: cobj,
    shat_p: cobj,
    qhat_p: cobj,
    pblt_p: cobj,
):
    geos = Ptr[float](geos_p)
    zi = Ptr[float](zi_p)
    paph = Ptr[float](paph_p)
    pap = Ptr[float](pap_p)
    zm = Ptr[float](zm_p)
    pblh = Ptr[float](pblh_p)
    t = Ptr[float](t_p)
    qh = Ptr[float](qh_p)
    zs = Ptr[float](zs_p)
    pf = Ptr[float](pf_p)
    zf = Ptr[float](zf_p)
    p = Ptr[float](p_p)
    z = Ptr[float](z_p)
    q = Ptr[float](q_p)
    s = Ptr[float](s_p)
    tp = Ptr[float](tp_p)
    shat = Ptr[float](shat_p)
    qhat = Ptr[float](qhat_p)
    pblt = Ptr[float](pblt_p)

    i = 0
    while i < ncol:
        zs[i] = geos[i] * rgrav
        pf[i + (pverp - 1) * pcols] = paph[i + (pverp - 1) * pcols] * 0.01
        zf[i + (pverp - 1) * pcols] = zi[i + (pverp - 1) * pcols] + zs[i]
        i += 1

    k = 0
    while k < pver:
        i = 0
        while i < ncol:
            idx = i + k * pcols
            p[idx] = pap[idx] * 0.01
            pf[idx] = paph[idx] * 0.01
            z[idx] = zm[idx] + zs[i]
            zf[idx] = zi[idx] + zs[i]
            i += 1
        k += 1

    kk = pver - 1
    while kk >= msg + 1:
        k = kk - 1
        i = 0
        while i < ncol:
            if abs(z[i + k * pcols] - zs[i] - pblh[i]) < (zf[i + k * pcols] - zf[i + (k + 1) * pcols]) * 0.5:
                pblt[i] = float(kk)
            i += 1
        kk -= 1

    k = 0
    while k < pver:
        i = 0
        while i < ncol:
            idx = i + k * pcols
            q[idx] = qh[idx]
            s[idx] = t[idx] + (grav / cpres) * z[idx]
            tp[idx] = 0.0
            shat[idx] = s[idx]
            qhat[idx] = q[idx]
            i += 1
        k += 1

@export
def zm_convr_pressure_state_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    msg: int,
    rgrav: float,
    grav: float,
    cpres: float,
    geos_p: cobj,
    zi_p: cobj,
    paph_p: cobj,
    pap_p: cobj,
    zm_p: cobj,
    pblh_p: cobj,
    t_p: cobj,
    qh_p: cobj,
    zs_p: cobj,
    pf_p: cobj,
    zf_p: cobj,
    p_p: cobj,
    z_p: cobj,
    q_p: cobj,
    s_p: cobj,
    tp_p: cobj,
    shat_p: cobj,
    qhat_p: cobj,
    pblt_p: cobj,
):
    zm_convr_pressure_state_stage_dispatch_codon(
        ncol,
        pcols,
        pver,
        pverp,
        msg,
        rgrav,
        grav,
        cpres,
        geos_p,
        zi_p,
        paph_p,
        pap_p,
        zm_p,
        pblh_p,
        t_p,
        qh_p,
        zs_p,
        pf_p,
        zf_p,
        p_p,
        z_p,
        q_p,
        s_p,
        tp_p,
        shat_p,
        qhat_p,
        pblt_p,
    )


@export
def zm_convr_gather_interface_stage_dispatch_codon(
    lengath: int,
    pcols: int,
    pver: int,
    pverp: int,
    msg: int,
    dpp_p: cobj,
    q_p: cobj,
    t_p: cobj,
    p_p: cobj,
    z_p: cobj,
    s_p: cobj,
    tp_p: cobj,
    zf_p: cobj,
    qstp_p: cobj,
    ideep_p: cobj,
    dp_p: cobj,
    qg_p: cobj,
    tg_p: cobj,
    pg_p: cobj,
    zg_p: cobj,
    sg_p: cobj,
    tpg_p: cobj,
    zfg_p: cobj,
    qstpg_p: cobj,
    ug_p: cobj,
    vg_p: cobj,
    shat_p: cobj,
    qhat_p: cobj,
):
    dpp = Ptr[float](dpp_p)
    q = Ptr[float](q_p)
    t = Ptr[float](t_p)
    p = Ptr[float](p_p)
    z = Ptr[float](z_p)
    s = Ptr[float](s_p)
    tp = Ptr[float](tp_p)
    zf = Ptr[float](zf_p)
    qstp = Ptr[float](qstp_p)
    ideep = Ptr[int](ideep_p)
    dp = Ptr[float](dp_p)
    qg = Ptr[float](qg_p)
    tg = Ptr[float](tg_p)
    pg = Ptr[float](pg_p)
    zg = Ptr[float](zg_p)
    sg = Ptr[float](sg_p)
    tpg = Ptr[float](tpg_p)
    zfg = Ptr[float](zfg_p)
    qstpg = Ptr[float](qstpg_p)
    ug = Ptr[float](ug_p)
    vg = Ptr[float](vg_p)
    shat = Ptr[float](shat_p)
    qhat = Ptr[float](qhat_p)

    k = 0
    while k < pver:
        i = 0
        while i < lengath:
            ii = ideep[i] - 1
            gidx = i + k * pcols
            sidx = ii + k * pcols
            dp[gidx] = 0.01 * dpp[sidx]
            qg[gidx] = q[sidx]
            tg[gidx] = t[sidx]
            pg[gidx] = p[sidx]
            zg[gidx] = z[sidx]
            sg[gidx] = s[sidx]
            tpg[gidx] = tp[sidx]
            zfg[gidx] = zf[sidx]
            qstpg[gidx] = qstp[sidx]
            ug[gidx] = 0.0
            vg[gidx] = 0.0
            i += 1
        k += 1

    i = 0
    while i < lengath:
        zfg[i + (pverp - 1) * pcols] = zf[(ideep[i] - 1) + (pverp - 1) * pcols]
        i += 1

    kk = msg + 2
    while kk <= pver:
        k = kk - 1
        km1 = kk - 2
        i = 0
        while i < lengath:
            idx = i + k * pcols
            km1idx = i + km1 * pcols
            sdifr = 0.0
            qdifr = 0.0
            if sg[idx] > 0.0 or sg[km1idx] > 0.0:
                sdifr = abs((sg[idx] - sg[km1idx]) / max(sg[km1idx], sg[idx]))
            if qg[idx] > 0.0 or qg[km1idx] > 0.0:
                qdifr = abs((qg[idx] - qg[km1idx]) / max(qg[km1idx], qg[idx]))
            if sdifr > 1.0e-6:
                shat[idx] = log(sg[km1idx] / sg[idx]) * sg[km1idx] * sg[idx] / (sg[km1idx] - sg[idx])
            else:
                shat[idx] = 0.5 * (sg[idx] + sg[km1idx])
            if qdifr > 1.0e-6:
                qhat[idx] = log(qg[km1idx] / qg[idx]) * qg[km1idx] * qg[idx] / (qg[km1idx] - qg[idx])
            else:
                qhat[idx] = 0.5 * (qg[idx] + qg[km1idx])
            i += 1
        kk += 1

@export
def zm_convr_gather_interface_shell_codon(
    lengath: int,
    pcols: int,
    pver: int,
    pverp: int,
    msg: int,
    dpp_p: cobj,
    q_p: cobj,
    t_p: cobj,
    p_p: cobj,
    z_p: cobj,
    s_p: cobj,
    tp_p: cobj,
    zf_p: cobj,
    qstp_p: cobj,
    ideep_p: cobj,
    dp_p: cobj,
    qg_p: cobj,
    tg_p: cobj,
    pg_p: cobj,
    zg_p: cobj,
    sg_p: cobj,
    tpg_p: cobj,
    zfg_p: cobj,
    qstpg_p: cobj,
    ug_p: cobj,
    vg_p: cobj,
    shat_p: cobj,
    qhat_p: cobj,
):
    zm_convr_gather_interface_stage_dispatch_codon(
        lengath,
        pcols,
        pver,
        pverp,
        msg,
        dpp_p,
        q_p,
        t_p,
        p_p,
        z_p,
        s_p,
        tp_p,
        zf_p,
        qstp_p,
        ideep_p,
        dp_p,
        qg_p,
        tg_p,
        pg_p,
        zg_p,
        sg_p,
        tpg_p,
        zfg_p,
        qstpg_p,
        ug_p,
        vg_p,
        shat_p,
        qhat_p,
    )


@export
def zm_convr_unit_stage_dispatch_codon(
    lengath: int,
    pcols: int,
    pver: int,
    msg: int,
    zfg_p: cobj,
    dp_p: cobj,
    du_p: cobj,
    eu_p: cobj,
    ed_p: cobj,
    cug_p: cobj,
    cmeg_p: cobj,
    rprdg_p: cobj,
    evpg_p: cobj,
    rppe_p: cobj,
):
    zfg = Ptr[float](zfg_p)
    dp = Ptr[float](dp_p)
    du = Ptr[float](du_p)
    eu = Ptr[float](eu_p)
    ed = Ptr[float](ed_p)
    cug = Ptr[float](cug_p)
    cmeg = Ptr[float](cmeg_p)
    rprdg = Ptr[float](rprdg_p)
    evpg = Ptr[float](evpg_p)
    rppe = Ptr[float](rppe_p)

    kk = msg + 1
    while kk <= pver:
        k = kk - 1
        i = 0
        while i < lengath:
            idx = i + k * pcols
            next_idx = i + (k + 1) * pcols
            du[idx] = du[idx] * (zfg[idx] - zfg[next_idx]) / dp[idx]
            eu[idx] = eu[idx] * (zfg[idx] - zfg[next_idx]) / dp[idx]
            ed[idx] = ed[idx] * (zfg[idx] - zfg[next_idx]) / dp[idx]
            cug[idx] = cug[idx] * (zfg[idx] - zfg[next_idx]) / dp[idx]
            cmeg[idx] = cmeg[idx] * (zfg[idx] - zfg[next_idx]) / dp[idx]
            rprdg[idx] = rprdg[idx] * (zfg[idx] - zfg[next_idx]) / dp[idx]
            evpg[idx] = evpg[idx] * (zfg[idx] - zfg[next_idx]) / dp[idx]
            rppe[idx] = rppe[idx] * (zfg[idx] - zfg[next_idx]) / dp[idx]
            i += 1
        kk += 1

@export
def zm_convr_unit_shell_codon(
    lengath: int,
    pcols: int,
    pver: int,
    msg: int,
    zfg_p: cobj,
    dp_p: cobj,
    du_p: cobj,
    eu_p: cobj,
    ed_p: cobj,
    cug_p: cobj,
    cmeg_p: cobj,
    rprdg_p: cobj,
    evpg_p: cobj,
    rppe_p: cobj,
):
    zm_convr_unit_stage_dispatch_codon(
        lengath,
        pcols,
        pver,
        msg,
        zfg_p,
        dp_p,
        du_p,
        eu_p,
        ed_p,
        cug_p,
        cmeg_p,
        rprdg_p,
        evpg_p,
        rppe_p,
    )


@export
def zm_convr_deep_state_stage_dispatch_codon(
    lengath: int,
    pcols: int,
    pver: int,
    msg: int,
    ideep_p: cobj,
    maxg_p: cobj,
    cape_p: cobj,
    tl_p: cobj,
    landfrac_p: cobj,
    dp_p: cobj,
    capeg_p: cobj,
    tlg_p: cobj,
    landfracg_p: cobj,
    dsubcld_p: cobj,
):
    ideep = Ptr[int](ideep_p)
    maxg = Ptr[int](maxg_p)
    cape = Ptr[float](cape_p)
    tl = Ptr[float](tl_p)
    landfrac = Ptr[float](landfrac_p)
    dp = Ptr[float](dp_p)
    capeg = Ptr[float](capeg_p)
    tlg = Ptr[float](tlg_p)
    landfracg = Ptr[float](landfracg_p)
    dsubcld = Ptr[float](dsubcld_p)

    i = 0
    while i < lengath:
        ii = ideep[i] - 1
        capeg[i] = cape[ii]
        tlg[i] = tl[ii]
        landfracg[i] = landfrac[ii]
        i += 1

    kk = msg + 1
    while kk <= pver:
        k = kk - 1
        i = 0
        while i < lengath:
            if kk >= maxg[i]:
                dsubcld[i] = dsubcld[i] + dp[i + k * pcols]
            i += 1
        kk += 1

@export
def zm_convr_deep_state_shell_codon(
    lengath: int,
    pcols: int,
    pver: int,
    msg: int,
    ideep_p: cobj,
    maxg_p: cobj,
    cape_p: cobj,
    tl_p: cobj,
    landfrac_p: cobj,
    dp_p: cobj,
    capeg_p: cobj,
    tlg_p: cobj,
    landfracg_p: cobj,
    dsubcld_p: cobj,
):
    zm_convr_deep_state_stage_dispatch_codon(
        lengath,
        pcols,
        pver,
        msg,
        ideep_p,
        maxg_p,
        cape_p,
        tl_p,
        landfrac_p,
        dp_p,
        capeg_p,
        tlg_p,
        landfracg_p,
        dsubcld_p,
    )


@export
def zm_convr_control_stage_dispatch_codon(
    stage: int,
    ncol: int,
    pcols: int,
    pver: int,
    capelmt: float,
    cape_p: cobj,
    lcl_p: cobj,
    lel_p: cobj,
    maxi_p: cobj,
    jt_p: cobj,
    capeg_p: cobj,
    lclg_p: cobj,
    lelg_p: cobj,
    maxg_p: cobj,
    tlg_p: cobj,
    dsubcld_p: cobj,
    lengath_p: cobj,
    index_p: cobj,
    ideep_p: cobj,
    ideep64_p: cobj,
    jt64_p: cobj,
    maxg64_p: cobj,
):
    if stage == 1:
        capeg = Ptr[float](capeg_p)
        lclg = Ptr[i32](lclg_p)
        lelg = Ptr[i32](lelg_p)
        maxg = Ptr[i32](maxg_p)
        tlg = Ptr[float](tlg_p)
        dsubcld = Ptr[float](dsubcld_p)

        i = 0
        while i < ncol:
            capeg[i] = 0.0
            lclg[i] = i32(1)
            lelg[i] = i32(pver)
            maxg[i] = i32(1)
            tlg[i] = 400.0
            dsubcld[i] = 0.0
            i += 1
    elif stage == 2:
        cape = Ptr[float](cape_p)
        lengath = Ptr[i32](lengath_p)
        index = Ptr[i32](index_p)
        ideep = Ptr[i32](ideep_p)
        ideep64 = Ptr[int](ideep64_p)

        count = 0
        i = 0
        while i < ncol:
            if cape[i] > capelmt:
                index[count] = i32(i + 1)
                count += 1
            i += 1

        lengath[0] = i32(count)
        ii = 0
        while ii < count:
            col = index[ii]
            ideep[ii] = col
            ideep64[ii] = int(col)
            ii += 1
    elif stage == 3:
        lcl = Ptr[i32](lcl_p)
        lel = Ptr[i32](lel_p)
        maxi = Ptr[i32](maxi_p)
        lclg = Ptr[i32](lclg_p)
        lelg = Ptr[i32](lelg_p)
        maxg = Ptr[i32](maxg_p)
        lengath = Ptr[i32](lengath_p)
        ideep = Ptr[i32](ideep_p)
        maxg64 = Ptr[int](maxg64_p)

        count = int(lengath[0])
        i = 0
        while i < count:
            col = int(ideep[i]) - 1
            lclg[i] = lcl[col]
            lelg[i] = lel[col]
            maxg[i] = maxi[col]
            maxg64[i] = int(maxg[i])
            i += 1
    elif stage == 4:
        jt = Ptr[i32](jt_p)
        maxg = Ptr[i32](maxg_p)
        lengath = Ptr[i32](lengath_p)
        jt64 = Ptr[int](jt64_p)
        maxg64 = Ptr[int](maxg64_p)

        count = int(lengath[0])
        i = 0
        while i < count:
            jt64[i] = int(jt[i])
            maxg64[i] = int(maxg[i])
            i += 1

@export
def zm_convr_control_shell_codon(
    stage: int,
    ncol: int,
    pcols: int,
    pver: int,
    capelmt: float,
    cape_p: cobj,
    lcl_p: cobj,
    lel_p: cobj,
    maxi_p: cobj,
    jt_p: cobj,
    capeg_p: cobj,
    lclg_p: cobj,
    lelg_p: cobj,
    maxg_p: cobj,
    tlg_p: cobj,
    dsubcld_p: cobj,
    lengath_p: cobj,
    index_p: cobj,
    ideep_p: cobj,
    ideep64_p: cobj,
    jt64_p: cobj,
    maxg64_p: cobj,
):
    zm_convr_control_stage_dispatch_codon(
        stage,
        ncol,
        pcols,
        pver,
        capelmt,
        cape_p,
        lcl_p,
        lel_p,
        maxi_p,
        jt_p,
        capeg_p,
        lclg_p,
        lelg_p,
        maxg_p,
        tlg_p,
        dsubcld_p,
        lengath_p,
        index_p,
        ideep_p,
        ideep64_p,
        jt64_p,
        maxg64_p,
    )


@export
def zm_convr_cloud_copy_stage_dispatch_codon(
    pcols: int,
    pver: int,
    cug_p: cobj,
    rprdg_p: cobj,
    wtcu_p: cobj,
    wtrpd_p: cobj,
):
    cug = Ptr[float](cug_p)
    rprdg = Ptr[float](rprdg_p)
    wtcu = Ptr[float](wtcu_p)
    wtrpd = Ptr[float](wtrpd_p)

    k = 0
    while k < pver:
        i = 0
        while i < pcols:
            idx = i + k * pcols
            wtcu[idx] = cug[idx]
            wtrpd[idx] = rprdg[idx]
            i += 1
        k += 1

@export
def zm_convr_cloud_copy_shell_codon(
    pcols: int,
    pver: int,
    cug_p: cobj,
    rprdg_p: cobj,
    wtcu_p: cobj,
    wtrpd_p: cobj,
):
    zm_convr_cloud_copy_stage_dispatch_codon(
        pcols,
        pver,
        cug_p,
        rprdg_p,
        wtcu_p,
        wtrpd_p,
    )


@export
def zm_convr_mflux_shell_codon(
    lengath: int,
    pcols: int,
    pver: int,
    pverp: int,
    msg: int,
    gravit: float,
    mb_p: cobj,
    mu_p: cobj,
    md_p: cobj,
    mc_p: cobj,
    du_p: cobj,
    eu_p: cobj,
    ed_p: cobj,
    cmeg_p: cobj,
    rprdg_p: cobj,
    cug_p: cobj,
    evpg_p: cobj,
    pflxg_p: cobj,
    rppe_p: cobj,
):
    mb = Ptr[float](mb_p)
    mu = Ptr[float](mu_p)
    md = Ptr[float](md_p)
    mc = Ptr[float](mc_p)
    du = Ptr[float](du_p)
    eu = Ptr[float](eu_p)
    ed = Ptr[float](ed_p)
    cmeg = Ptr[float](cmeg_p)
    rprdg = Ptr[float](rprdg_p)
    cug = Ptr[float](cug_p)
    evpg = Ptr[float](evpg_p)
    pflxg = Ptr[float](pflxg_p)
    rppe = Ptr[float](rppe_p)

    kk = msg + 1
    while kk <= pver:
        k = kk - 1
        i = 0
        while i < lengath:
            idx = i + k * pcols
            mu[idx] = mu[idx] * mb[i]
            md[idx] = md[idx] * mb[i]
            mc[idx] = mc[idx] * mb[i]
            du[idx] = du[idx] * mb[i]
            eu[idx] = eu[idx] * mb[i]
            ed[idx] = ed[idx] * mb[i]
            cmeg[idx] = cmeg[idx] * mb[i]
            rprdg[idx] = rprdg[idx] * mb[i]
            cug[idx] = cug[idx] * mb[i]
            evpg[idx] = evpg[idx] * mb[i]
            pflxg[i + (k + 1) * pcols] = pflxg[i + (k + 1) * pcols] * mb[i] * 100.0 / gravit
            rppe[idx] = rppe[idx] * mb[i]
            i += 1
        kk += 1


@export
def zm_convr_closure_limit_shell_codon(
    lengath: int,
    pcols: int,
    pver: int,
    msg: int,
    no_deep_pbl_flag: int,
    delt: float,
    ideep_p: cobj,
    jt_p: cobj,
    zm_p: cobj,
    pblh_p: cobj,
    mu_p: cobj,
    dp_p: cobj,
    mb_p: cobj,
    mumax_p: cobj,
):
    ideep = Ptr[int](ideep_p)
    jt = Ptr[int](jt_p)
    zm = Ptr[float](zm_p)
    pblh = Ptr[float](pblh_p)
    mu = Ptr[float](mu_p)
    dp = Ptr[float](dp_p)
    mb = Ptr[float](mb_p)
    mumax = Ptr[float](mumax_p)

    i = 0
    while i < lengath:
        mumax[i] = 0.0
        i += 1

    kk = msg + 2
    while kk <= pver:
        k = kk - 1
        i = 0
        while i < lengath:
            idx = i + k * pcols
            mumax[i] = max(mumax[i], mu[idx] / dp[idx])
            i += 1
        kk += 1

    i = 0
    while i < lengath:
        if mumax[i] > 0.0:
            mb[i] = min(mb[i], 0.5 / (delt * mumax[i]))
        else:
            mb[i] = 0.0
        i += 1

    if no_deep_pbl_flag != 0:
        i = 0
        while i < lengath:
            ii = ideep[i] - 1
            k = jt[i] - 1
            if zm[ii + k * pcols] < pblh[ii]:
                mb[i] = 0.0
            i += 1


@export
def zm_closure_codon(
    pcols: int,
    pver: int,
    il1g: int,
    il2g: int,
    msg: int,
    rd: float,
    grav: float,
    cp: float,
    rl: float,
    eps1: float,
    tau: float,
    capelmt: float,
    q_p: cobj,
    t_p: cobj,
    p_p: cobj,
    z_p: cobj,
    s_p: cobj,
    tp_p: cobj,
    qs_p: cobj,
    qu_p: cobj,
    su_p: cobj,
    mc_p: cobj,
    du_p: cobj,
    mu_p: cobj,
    md_p: cobj,
    qd_p: cobj,
    sd_p: cobj,
    qhat_p: cobj,
    shat_p: cobj,
    dp_p: cobj,
    qstp_p: cobj,
    zf_p: cobj,
    ql_p: cobj,
    dsubcld_p: cobj,
    mb_p: cobj,
    cape_p: cobj,
    tl_p: cobj,
    lcl_p: cobj,
    lel_p: cobj,
    jt_p: cobj,
    mx_p: cobj,
    dtpdt_p: cobj,
    dqsdtp_p: cobj,
    dtmdt_p: cobj,
    dqmdt_p: cobj,
    dboydt_p: cobj,
    thetavp_p: cobj,
    thetavm_p: cobj,
    dtbdt_p: cobj,
    dqbdt_p: cobj,
    dtldt_p: cobj,
    dadt_p: cobj,
):
    q = Ptr[float](q_p)
    t = Ptr[float](t_p)
    p = Ptr[float](p_p)
    z = Ptr[float](z_p)
    s = Ptr[float](s_p)
    tp = Ptr[float](tp_p)
    qs = Ptr[float](qs_p)
    qu = Ptr[float](qu_p)
    su = Ptr[float](su_p)
    mc = Ptr[float](mc_p)
    du = Ptr[float](du_p)
    mu = Ptr[float](mu_p)
    md = Ptr[float](md_p)
    qd = Ptr[float](qd_p)
    sd = Ptr[float](sd_p)
    qhat = Ptr[float](qhat_p)
    shat = Ptr[float](shat_p)
    dp = Ptr[float](dp_p)
    qstp = Ptr[float](qstp_p)
    zf = Ptr[float](zf_p)
    ql = Ptr[float](ql_p)
    dsubcld = Ptr[float](dsubcld_p)
    mb = Ptr[float](mb_p)
    cape = Ptr[float](cape_p)
    tl = Ptr[float](tl_p)
    lcl = Ptr[i32](lcl_p)
    lel = Ptr[i32](lel_p)
    jt = Ptr[i32](jt_p)
    mx = Ptr[i32](mx_p)
    dtpdt = Ptr[float](dtpdt_p)
    dqsdtp = Ptr[float](dqsdtp_p)
    dtmdt = Ptr[float](dtmdt_p)
    dqmdt = Ptr[float](dqmdt_p)
    dboydt = Ptr[float](dboydt_p)
    thetavp = Ptr[float](thetavp_p)
    thetavm = Ptr[float](thetavm_p)
    dtbdt = Ptr[float](dtbdt_p)
    dqbdt = Ptr[float](dqbdt_p)
    dtldt = Ptr[float](dtldt_p)
    dadt = Ptr[float](dadt_p)

    ii0 = il1g - 1
    ii1 = il2g - 1

    ii = ii0
    while ii <= ii1:
        mx_i = int(mx[ii])
        mxidx = ii + (mx_i - 1) * pcols
        mb[ii] = 0.0
        eb = p[mxidx] * q[mxidx] / (eps1 + q[mxidx])
        dtbdt[ii] = (1.0 / dsubcld[ii]) * (
            mu[mxidx] * (shat[mxidx] - su[mxidx])
            + md[mxidx] * (shat[mxidx] - sd[mxidx])
        )
        dqbdt[ii] = (1.0 / dsubcld[ii]) * (
            mu[mxidx] * (qhat[mxidx] - qu[mxidx])
            + md[mxidx] * (qhat[mxidx] - qd[mxidx])
        )
        debdt = eps1 * p[mxidx] / (eps1 + q[mxidx]) ** 2 * dqbdt[ii]
        dtldt[ii] = (
            -2840.0
            * (3.5 / t[mxidx] * dtbdt[ii] - debdt / eb)
            / (3.5 * log(t[mxidx]) - log(eb) - 4.805) ** 2
        )
        ii += 1

    kk = msg + 1
    while kk <= pver:
        k = kk - 1
        ii = ii0
        while ii <= ii1:
            idx = ii + k * pcols
            dtmdt[idx] = 0.0
            dqmdt[idx] = 0.0
            ii += 1
        kk += 1

    kk = msg + 1
    while kk <= pver - 1:
        k = kk - 1
        kp1 = kk
        ii = ii0
        while ii <= ii1:
            if kk == int(jt[ii]):
                idx = ii + k * pcols
                kp1idx = ii + kp1 * pcols
                dtmdt[idx] = (1.0 / dp[idx]) * (
                    mu[kp1idx] * (su[kp1idx] - shat[kp1idx] - rl / cp * ql[kp1idx])
                    + md[kp1idx] * (sd[kp1idx] - shat[kp1idx])
                )
                dqmdt[idx] = (1.0 / dp[idx]) * (
                    mu[kp1idx] * (qu[kp1idx] - qhat[kp1idx] + ql[kp1idx])
                    + md[kp1idx] * (qd[kp1idx] - qhat[kp1idx])
                )
            ii += 1
        kk += 1

    beta = 0.0
    kk = msg + 1
    while kk <= pver - 1:
        k = kk - 1
        kp1 = kk
        ii = ii0
        while ii <= ii1:
            if kk > int(jt[ii]) and kk < int(mx[ii]):
                idx = ii + k * pcols
                kp1idx = ii + kp1 * pcols
                dtmdt[idx] = (
                    mc[idx] * (shat[idx] - s[idx])
                    + mc[kp1idx] * (s[idx] - shat[kp1idx])
                ) / dp[idx] - rl / cp * du[idx] * (beta * ql[idx] + (1.0 - beta) * ql[kp1idx])
                dqmdt[idx] = (
                    mu[kp1idx] * (qu[kp1idx] - qhat[kp1idx] + cp / rl * (su[kp1idx] - s[idx]))
                    - mu[idx] * (qu[idx] - qhat[idx] + cp / rl * (su[idx] - s[idx]))
                    + md[kp1idx] * (qd[kp1idx] - qhat[kp1idx] + cp / rl * (sd[kp1idx] - s[idx]))
                    - md[idx] * (qd[idx] - qhat[idx] + cp / rl * (sd[idx] - s[idx]))
                ) / dp[idx] + du[idx] * (beta * ql[idx] + (1.0 - beta) * ql[kp1idx])
            ii += 1
        kk += 1

    kk = msg + 1
    while kk <= pver:
        k = kk - 1
        ii = ii0
        while ii <= ii1:
            if kk >= int(lel[ii]) and kk <= int(lcl[ii]):
                idx = ii + k * pcols
                mxidx = ii + (int(mx[ii]) - 1) * pcols
                thetavp[idx] = tp[idx] * (1000.0 / p[idx]) ** (rd / cp) * (1.0 + 1.608 * qstp[idx] - q[mxidx])
                thetavm[idx] = t[idx] * (1000.0 / p[idx]) ** (rd / cp) * (1.0 + 0.608 * q[idx])
                dqsdtp[idx] = qstp[idx] * (1.0 + qstp[idx] / eps1) * eps1 * rl / (rd * tp[idx] ** 2)
                dtpdt[idx] = tp[idx] / (1.0 + rl / cp * (dqsdtp[idx] - qstp[idx] / tp[idx])) * (
                    dtbdt[ii] / t[mxidx]
                    + rl / cp * (dqbdt[ii] / tl[ii] - q[mxidx] / tl[ii] ** 2 * dtldt[ii])
                )
                dboydt[idx] = (
                    (
                        dtpdt[idx] / tp[idx]
                        + 1.0 / (1.0 + 1.608 * qstp[idx] - q[mxidx])
                        * (1.608 * dqsdtp[idx] * dtpdt[idx] - dqbdt[ii])
                    )
                    - (
                        dtmdt[idx] / t[idx]
                        + 0.608 / (1.0 + 0.608 * q[idx]) * dqmdt[idx]
                    )
                ) * grav * thetavp[idx] / thetavm[idx]
            ii += 1
        kk += 1

    kk = msg + 1
    while kk <= pver:
        k = kk - 1
        ii = ii0
        while ii <= ii1:
            if kk > int(lcl[ii]) and kk < int(mx[ii]):
                idx = ii + k * pcols
                mxidx = ii + (int(mx[ii]) - 1) * pcols
                thetavp[idx] = tp[idx] * (1000.0 / p[idx]) ** (rd / cp) * (1.0 + 0.608 * q[mxidx])
                thetavm[idx] = t[idx] * (1000.0 / p[idx]) ** (rd / cp) * (1.0 + 0.608 * q[idx])
                dboydt[idx] = (
                    dtbdt[ii] / t[mxidx]
                    + 0.608 / (1.0 + 0.608 * q[mxidx]) * dqbdt[ii]
                    - dtmdt[idx] / t[idx]
                    - 0.608 / (1.0 + 0.608 * q[idx]) * dqmdt[idx]
                ) * grav * thetavp[idx] / thetavm[idx]
            ii += 1
        kk += 1

    ii = ii0
    while ii <= ii1:
        dadt[ii] = 0.0
        ii += 1

    kmin = int(lel[ii0])
    kmax = int(mx[ii0])
    ii = ii0 + 1
    while ii <= ii1:
        kmin = min(kmin, int(lel[ii]))
        kmax = max(kmax, int(mx[ii]))
        ii += 1
    kmax -= 1

    kk = kmin
    while kk <= kmax:
        k = kk - 1
        kp1 = kk
        ii = ii0
        while ii <= ii1:
            if kk >= int(lel[ii]) and kk <= int(mx[ii]) - 1:
                idx = ii + k * pcols
                dadt[ii] = dadt[ii] + dboydt[idx] * (zf[idx] - zf[ii + kp1 * pcols])
            ii += 1
        kk += 1

    ii = ii0
    while ii <= ii1:
        dltaa = -1.0 * (cape[ii] - capelmt)
        if dadt[ii] != 0.0:
            mb[ii] = max(dltaa / tau / dadt[ii], 0.0)
        ii += 1


@export
def zm_q1q2_pjr_codon(
    lengath: int,
    pcols: int,
    pver: int,
    msg: int,
    cp: float,
    rl: float,
    q_p: cobj,
    qu_p: cobj,
    su_p: cobj,
    du_p: cobj,
    qhat_p: cobj,
    shat_p: cobj,
    dp_p: cobj,
    mu_p: cobj,
    md_p: cobj,
    sd_p: cobj,
    qd_p: cobj,
    ql_p: cobj,
    dsubcld_p: cobj,
    jt_p: cobj,
    mx_p: cobj,
    dl_p: cobj,
    evp_p: cobj,
    cu_p: cobj,
    dqdt_p: cobj,
    dsdt_p: cobj,
):
    q = Ptr[float](q_p)
    qu = Ptr[float](qu_p)
    su = Ptr[float](su_p)
    du = Ptr[float](du_p)
    qhat = Ptr[float](qhat_p)
    shat = Ptr[float](shat_p)
    dp = Ptr[float](dp_p)
    mu = Ptr[float](mu_p)
    md = Ptr[float](md_p)
    sd = Ptr[float](sd_p)
    qd = Ptr[float](qd_p)
    ql = Ptr[float](ql_p)
    dsubcld = Ptr[float](dsubcld_p)
    jt = Ptr[int](jt_p)
    mx = Ptr[int](mx_p)
    dl = Ptr[float](dl_p)
    evp = Ptr[float](evp_p)
    cu = Ptr[float](cu_p)
    dqdt = Ptr[float](dqdt_p)
    dsdt = Ptr[float](dsdt_p)

    kk = msg + 1
    while kk <= pver:
        k = kk - 1
        i = 0
        while i < lengath:
            idx = i + k * pcols
            dsdt[idx] = 0.0
            dqdt[idx] = 0.0
            dl[idx] = 0.0
            i += 1
        kk += 1

    ktm = pver
    kbm = pver
    i = 0
    while i < lengath:
        ktm = min(ktm, jt[i])
        kbm = min(kbm, mx[i])
        i += 1

    kk = ktm
    while kk <= pver - 1:
        k = kk - 1
        kp1 = kk
        i = 0
        while i < lengath:
            idx = i + k * pcols
            kp1idx = i + kp1 * pcols
            emc = -cu[idx] + evp[idx]
            dsdt[idx] = -rl / cp * emc + (
                +mu[kp1idx] * (su[kp1idx] - shat[kp1idx])
                - mu[idx] * (su[idx] - shat[idx])
                + md[kp1idx] * (sd[kp1idx] - shat[kp1idx])
                - md[idx] * (sd[idx] - shat[idx])
            ) / dp[idx]
            dqdt[idx] = emc + (
                +mu[kp1idx] * (qu[kp1idx] - qhat[kp1idx])
                - mu[idx] * (qu[idx] - qhat[idx])
                + md[kp1idx] * (qd[kp1idx] - qhat[kp1idx])
                - md[idx] * (qd[idx] - qhat[idx])
            ) / dp[idx]
            dl[idx] = du[idx] * ql[kp1idx]
            i += 1
        kk += 1

    kk = kbm
    while kk <= pver:
        k = kk - 1
        i = 0
        while i < lengath:
            idx = i + k * pcols
            if kk == mx[i]:
                dsdt[idx] = (1.0 / dsubcld[i]) * (
                    -mu[idx] * (su[idx] - shat[idx])
                    - md[idx] * (sd[idx] - shat[idx])
                )
                dqdt[idx] = (1.0 / dsubcld[i]) * (
                    -mu[idx] * (qu[idx] - qhat[idx])
                    - md[idx] * (qd[idx] - qhat[idx])
                )
            elif kk > mx[i]:
                km1idx = i + (k - 1) * pcols
                dsdt[idx] = dsdt[km1idx]
                dqdt[idx] = dqdt[km1idx]
            i += 1
        kk += 1


@export
def zm_convr_tail_shell_codon(
    ncol: int,
    lengath: int,
    pcols: int,
    pver: int,
    pverp: int,
    msg: int,
    delt: float,
    cpres: float,
    rgrav: float,
    gravit: float,
    ideep_p: cobj,
    jt_p: cobj,
    maxg_p: cobj,
    qh_p: cobj,
    dpp_p: cobj,
    q_p: cobj,
    dqdt_p: cobj,
    dsdt_p: cobj,
    cmeg_p: cobj,
    rprdg_p: cobj,
    du_p: cobj,
    mc_p: cobj,
    dlg_p: cobj,
    cug_p: cobj,
    evpg_p: cobj,
    tug_p: cobj,
    tdg_p: cobj,
    pflxg_p: cobj,
    qlg_p: cobj,
    qtnd_p: cobj,
    cme_p: cobj,
    rprd_p: cobj,
    zdu_p: cobj,
    mcon_p: cobj,
    heat_p: cobj,
    dlf_p: cobj,
    cu_p: cobj,
    evp_p: cobj,
    tu_p: cobj,
    td_p: cobj,
    pflx_p: cobj,
    ql_p: cobj,
    jctop_p: cobj,
    jcbot_p: cobj,
    prec_p: cobj,
    rliq_p: cobj,
):
    ideep = Ptr[int](ideep_p)
    jt = Ptr[int](jt_p)
    maxg = Ptr[int](maxg_p)
    qh = Ptr[float](qh_p)
    dpp = Ptr[float](dpp_p)
    q = Ptr[float](q_p)
    dqdt = Ptr[float](dqdt_p)
    dsdt = Ptr[float](dsdt_p)
    cmeg = Ptr[float](cmeg_p)
    rprdg = Ptr[float](rprdg_p)
    du = Ptr[float](du_p)
    mc = Ptr[float](mc_p)
    dlg = Ptr[float](dlg_p)
    cug = Ptr[float](cug_p)
    evpg = Ptr[float](evpg_p)
    tug = Ptr[float](tug_p)
    tdg = Ptr[float](tdg_p)
    pflxg = Ptr[float](pflxg_p)
    qlg = Ptr[float](qlg_p)
    qtnd = Ptr[float](qtnd_p)
    cme = Ptr[float](cme_p)
    rprd = Ptr[float](rprd_p)
    zdu = Ptr[float](zdu_p)
    mcon = Ptr[float](mcon_p)
    heat = Ptr[float](heat_p)
    dlf = Ptr[float](dlf_p)
    cu = Ptr[float](cu_p)
    evp = Ptr[float](evp_p)
    tu = Ptr[float](tu_p)
    td = Ptr[float](td_p)
    pflx = Ptr[float](pflx_p)
    ql = Ptr[float](ql_p)
    jctop = Ptr[float](jctop_p)
    jcbot = Ptr[float](jcbot_p)
    prec = Ptr[float](prec_p)
    rliq = Ptr[float](rliq_p)

    kk = msg + 1
    while kk <= pver:
        k = kk - 1
        i = 0
        while i < lengath:
            ii = ideep[i] - 1
            gidx = i + k * pcols
            sidx = ii + k * pcols
            q[sidx] = qh[sidx] + 2.0 * delt * dqdt[gidx]
            qtnd[sidx] = dqdt[gidx]
            cme[sidx] = cmeg[gidx]
            rprd[sidx] = rprdg[gidx]
            zdu[sidx] = du[gidx]
            mcon[sidx] = mc[gidx]
            heat[sidx] = dsdt[gidx] * cpres
            dlf[sidx] = dlg[gidx]
            cu[sidx] = cug[gidx]
            evp[sidx] = evpg[gidx]
            tu[sidx] = tug[gidx]
            td[sidx] = tdg[gidx]
            pflx[sidx] = pflxg[gidx]
            ql[sidx] = qlg[gidx]
            i += 1
        kk += 1

    i = 0
    while i < lengath:
        ii = ideep[i] - 1
        jctop[ii] = float(jt[i])
        jcbot[ii] = float(maxg[i])
        pflx[ii + (pverp - 1) * pcols] = pflxg[i + (pverp - 1) * pcols]
        i += 1

    kk = pver
    while kk >= msg + 1:
        k = kk - 1
        i = 0
        while i < ncol:
            idx = i + k * pcols
            prec[i] = prec[i] - dpp[idx] * (q[idx] - qh[idx]) - dpp[idx] * dlf[idx] * 2 * delt
            i += 1
        kk -= 1

    i = 0
    while i < ncol:
        prec[i] = rgrav * max(prec[i], 0.0) / (2.0 * delt) / 1000.0
        i += 1

    k = 0
    while k < pver:
        i = 0
        while i < ncol:
            idx = i + k * pcols
            rliq[i] = rliq[i] + dlf[idx] * dpp[idx] / gravit
            i += 1
        k += 1

    i = 0
    while i < ncol:
        rliq[i] = rliq[i] / 1000.0
        i += 1


@export
def zm_convr_finish_stage_dispatch_codon(
    ncol: int,
    lengath: int,
    pcols: int,
    pver: int,
    pverp: int,
    msg: int,
    no_deep_pbl_flag: int,
    delt: float,
    cpres: float,
    rl: float,
    rgrav: float,
    grav: float,
    gravit: float,
    ideep_p: cobj,
    jt_p: cobj,
    maxg_p: cobj,
    zm_p: cobj,
    pblh_p: cobj,
    mb_p: cobj,
    mumax_p: cobj,
    qg_p: cobj,
    qu_p: cobj,
    su_p: cobj,
    qhat_p: cobj,
    shat_p: cobj,
    dp_p: cobj,
    mu_p: cobj,
    md_p: cobj,
    sd_p: cobj,
    qd_p: cobj,
    qlg_p: cobj,
    dsubcld_p: cobj,
    dlg_p: cobj,
    evpg_p: cobj,
    cug_p: cobj,
    dqdt_p: cobj,
    dsdt_p: cobj,
    mc_p: cobj,
    du_p: cobj,
    eu_p: cobj,
    ed_p: cobj,
    cmeg_p: cobj,
    rprdg_p: cobj,
    pflxg_p: cobj,
    rppe_p: cobj,
    qh_p: cobj,
    dpp_p: cobj,
    q_p: cobj,
    tug_p: cobj,
    tdg_p: cobj,
    qtnd_p: cobj,
    cme_p: cobj,
    rprd_p: cobj,
    zdu_p: cobj,
    mcon_p: cobj,
    heat_p: cobj,
    dlf_p: cobj,
    cu_p: cobj,
    evp_p: cobj,
    tu_p: cobj,
    td_p: cobj,
    pflx_p: cobj,
    ql_p: cobj,
    jctop_p: cobj,
    jcbot_p: cobj,
    prec_p: cobj,
    rliq_p: cobj,
):
    zm_convr_closure_limit_shell_codon(
        lengath,
        pcols,
        pver,
        msg,
        no_deep_pbl_flag,
        delt,
        ideep_p,
        jt_p,
        zm_p,
        pblh_p,
        mu_p,
        dp_p,
        mb_p,
        mumax_p,
    )
    zm_convr_mflux_shell_codon(
        lengath,
        pcols,
        pver,
        pverp,
        msg,
        grav,
        mb_p,
        mu_p,
        md_p,
        mc_p,
        du_p,
        eu_p,
        ed_p,
        cmeg_p,
        rprdg_p,
        cug_p,
        evpg_p,
        pflxg_p,
        rppe_p,
    )
    zm_q1q2_pjr_codon(
        lengath,
        pcols,
        pver,
        msg,
        cpres,
        rl,
        qg_p,
        qu_p,
        su_p,
        du_p,
        qhat_p,
        shat_p,
        dp_p,
        mu_p,
        md_p,
        sd_p,
        qd_p,
        qlg_p,
        dsubcld_p,
        jt_p,
        maxg_p,
        dlg_p,
        evpg_p,
        cug_p,
        dqdt_p,
        dsdt_p,
    )
    zm_convr_tail_shell_codon(
        ncol,
        lengath,
        pcols,
        pver,
        pverp,
        msg,
        delt,
        cpres,
        rgrav,
        gravit,
        ideep_p,
        jt_p,
        maxg_p,
        qh_p,
        dpp_p,
        q_p,
        dqdt_p,
        dsdt_p,
        cmeg_p,
        rprdg_p,
        du_p,
        mc_p,
        dlg_p,
        cug_p,
        evpg_p,
        tug_p,
        tdg_p,
        pflxg_p,
        qlg_p,
        qtnd_p,
        cme_p,
        rprd_p,
        zdu_p,
        mcon_p,
        heat_p,
        dlf_p,
        cu_p,
        evp_p,
        tu_p,
        td_p,
        pflx_p,
        ql_p,
        jctop_p,
        jcbot_p,
        prec_p,
        rliq_p,
    )

@export
def zm_convr_finish_shell_codon(
    ncol: int,
    lengath: int,
    pcols: int,
    pver: int,
    pverp: int,
    msg: int,
    no_deep_pbl_flag: int,
    delt: float,
    cpres: float,
    rl: float,
    rgrav: float,
    grav: float,
    gravit: float,
    ideep_p: cobj,
    jt_p: cobj,
    maxg_p: cobj,
    zm_p: cobj,
    pblh_p: cobj,
    mb_p: cobj,
    mumax_p: cobj,
    qg_p: cobj,
    qu_p: cobj,
    su_p: cobj,
    qhat_p: cobj,
    shat_p: cobj,
    dp_p: cobj,
    mu_p: cobj,
    md_p: cobj,
    sd_p: cobj,
    qd_p: cobj,
    qlg_p: cobj,
    dsubcld_p: cobj,
    dlg_p: cobj,
    evpg_p: cobj,
    cug_p: cobj,
    dqdt_p: cobj,
    dsdt_p: cobj,
    mc_p: cobj,
    du_p: cobj,
    eu_p: cobj,
    ed_p: cobj,
    cmeg_p: cobj,
    rprdg_p: cobj,
    pflxg_p: cobj,
    rppe_p: cobj,
    qh_p: cobj,
    dpp_p: cobj,
    q_p: cobj,
    tug_p: cobj,
    tdg_p: cobj,
    qtnd_p: cobj,
    cme_p: cobj,
    rprd_p: cobj,
    zdu_p: cobj,
    mcon_p: cobj,
    heat_p: cobj,
    dlf_p: cobj,
    cu_p: cobj,
    evp_p: cobj,
    tu_p: cobj,
    td_p: cobj,
    pflx_p: cobj,
    ql_p: cobj,
    jctop_p: cobj,
    jcbot_p: cobj,
    prec_p: cobj,
    rliq_p: cobj,
):
    zm_convr_finish_stage_dispatch_codon(
        ncol,
        lengath,
        pcols,
        pver,
        pverp,
        msg,
        no_deep_pbl_flag,
        delt,
        cpres,
        rl,
        rgrav,
        grav,
        gravit,
        ideep_p,
        jt_p,
        maxg_p,
        zm_p,
        pblh_p,
        mb_p,
        mumax_p,
        qg_p,
        qu_p,
        su_p,
        qhat_p,
        shat_p,
        dp_p,
        mu_p,
        md_p,
        sd_p,
        qd_p,
        qlg_p,
        dsubcld_p,
        dlg_p,
        evpg_p,
        cug_p,
        dqdt_p,
        dsdt_p,
        mc_p,
        du_p,
        eu_p,
        ed_p,
        cmeg_p,
        rprdg_p,
        pflxg_p,
        rppe_p,
        qh_p,
        dpp_p,
        q_p,
        tug_p,
        tdg_p,
        qtnd_p,
        cme_p,
        rprd_p,
        zdu_p,
        mcon_p,
        heat_p,
        dlf_p,
        cu_p,
        evp_p,
        tu_p,
        td_p,
        pflx_p,
        ql_p,
        jctop_p,
        jcbot_p,
        prec_p,
        rliq_p,
    )


@export
def zm_conv_workspace_init_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    ftem_p: cobj,
    mu_out_p: cobj,
    md_out_p: cobj,
    wind_tends_p: cobj,
):
    ftem = Ptr[float](ftem_p)
    mu_out = Ptr[float](mu_out_p)
    md_out = Ptr[float](md_out_p)
    wind_tends = Ptr[float](wind_tends_p)
    plane = pcols * pver

    k = 0
    while k < pver:
        i = 0
        while i < pcols:
            idx = i + k * pcols
            ftem[idx] = 0.0
            mu_out[idx] = 0.0
            md_out[idx] = 0.0
            i += 1
        k += 1

    k = 0
    while k < pver:
        i = 0
        while i < ncol:
            idx = i + k * pcols
            wind_tends[idx] = 0.0
            wind_tends[idx + plane] = 0.0
            i += 1
        k += 1


@export
def zm_wtrc_convr_prep_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    wtrc_nwset: int,
    wtdlf_p: cobj,
    wtrprd_p: cobj,
    wtprect_p: cobj,
    rprd_p: cobj,
    prec_p: cobj,
):
    wtdlf = Ptr[float](wtdlf_p)
    wtrprd = Ptr[float](wtrprd_p)
    wtprect = Ptr[float](wtprect_p)
    rprd = Ptr[float](rprd_p)
    prec = Ptr[float](prec_p)
    plane = pcols * pver

    m = 0
    while m < wtrc_nwset:
        k = 0
        while k < pver:
            i = 0
            while i < ncol:
                wtdlf[i + k * pcols + m * plane] = 0.0
                i += 1
            k += 1
        m += 1

    m = 0
    while m < pcnst:
        k = 0
        while k < pver:
            i = 0
            while i < ncol:
                wtrprd[i + k * pcols + m * plane] = 0.0
                i += 1
            k += 1
        i = 0
        while i < ncol:
            wtprect[i + m * pcols] = 0.0
            i += 1
        m += 1

    k = 0
    while k < pver:
        i = 0
        while i < ncol:
            idx = i + k * pcols
            wtrprd[idx] = rprd[idx]
            i += 1
        k += 1

    i = 0
    while i < ncol:
        wtprect[i] = prec[i]
        i += 1


@export
def zm_wtrc_precip_assign_shell_codon(
    pcols: int,
    vap_type: int,
    wtprect_p: cobj,
    wtsnowt_p: cobj,
    wtprec_p: cobj,
    wtsnow_p: cobj,
):
    wtprect = Ptr[float](wtprect_p)
    wtsnowt = Ptr[float](wtsnowt_p)
    wtprec = Ptr[float](wtprec_p)
    wtsnow = Ptr[float](wtsnow_p)
    offset = (vap_type - 1) * pcols

    i = 0
    while i < pcols:
        snow_val = wtsnowt[i + offset]
        wtprec[i] = wtprect[i + offset] - snow_val
        wtsnow[i] = snow_val
        i += 1


@export
def zm_conv_evap_prep_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    dp_cldliq_p: cobj,
    dp_cldice_p: cobj,
):
    dp_cldliq = Ptr[float](dp_cldliq_p)
    dp_cldice = Ptr[float](dp_cldice_p)

    k = 0
    while k < pver:
        i = 0
        while i < ncol:
            idx = i + k * pcols
            dp_cldliq[idx] = 0.0
            dp_cldice[idx] = 0.0
            i += 1
        k += 1


@export
def zm_conv_evap_post_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    ixorg: int,
    do_org: int,
    ztodt: float,
    evapcdp_p: cobj,
    ptend_q_p: cobj,
    state_q_p: cobj,
):
    evapcdp = Ptr[float](evapcdp_p)
    ptend_q = Ptr[float](ptend_q_p)
    state_q = Ptr[float](state_q_p)
    plane = pcols * pver

    k = 0
    while k < pver:
        i = 0
        while i < ncol:
            idx = i + k * pcols
            evapcdp[idx] = ptend_q[idx]
            i += 1
        k += 1

    if do_org != 0:
        org_offset = (ixorg - 1) * plane
        k = 0
        while k < pver:
            i = 0
            while i < ncol:
                idx = i + k * pcols
                org_idx = idx + org_offset
                val = (50.0 * 1000.0 * 1000.0 * abs(evapcdp[idx])) - (state_q[org_idx] / 10800.0)
                if val < 0.0:
                    val = 0.0
                if val > 1.0:
                    val = 1.0
                ptend_q[org_idx] = (val - state_q[org_idx]) / ztodt
                i += 1
            k += 1


@export
def zm_conv_evap_hist_shell_codon(
    mode: int,
    ncol: int,
    pcols: int,
    pver: int,
    cpair: float,
    ptend_s_p: cobj,
    tend_s_snwprd_p: cobj,
    tend_s_snwevmlt_p: cobj,
    ftem_p: cobj,
):
    ptend_s = Ptr[float](ptend_s_p)
    tend_s_snwprd = Ptr[float](tend_s_snwprd_p)
    tend_s_snwevmlt = Ptr[float](tend_s_snwevmlt_p)
    ftem = Ptr[float](ftem_p)

    k = 0
    while k < pver:
        i = 0
        while i < ncol:
            idx = i + k * pcols
            if mode == 1:
                ftem[idx] = ptend_s[idx] / cpair
            elif mode == 2:
                ftem[idx] = tend_s_snwprd[idx] / cpair
            else:
                ftem[idx] = tend_s_snwevmlt[idx] / cpair
            i += 1
        k += 1


@export
def zm_conv_evap_main_stage_dispatch_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    pergro: int,
    do_org: int,
    gravit: float,
    latvap: float,
    latice: float,
    tmelt: float,
    ke: float,
    ke_lnd: float,
    deltat: float,
    t_p: cobj,
    q_p: cobj,
    pdel_p: cobj,
    landfrac_p: cobj,
    prdprec_p: cobj,
    cldfrc_p: cobj,
    qs_p: cobj,
    fsnow_conv_p: cobj,
    prec_p: cobj,
    snow_p: cobj,
    tend_s_p: cobj,
    tend_q_p: cobj,
    tend_s_snwprd_p: cobj,
    tend_s_snwevmlt_p: cobj,
    evpstore_p: cobj,
    substore_p: cobj,
    ntprprd_p: cobj,
    ntsnprd_p: cobj,
    flxprec_p: cobj,
    flxsnow_p: cobj,
    evpvint_p: cobj,
):
    t = Ptr[float](t_p)
    q = Ptr[float](q_p)
    pdel = Ptr[float](pdel_p)
    landfrac = Ptr[float](landfrac_p)
    prdprec = Ptr[float](prdprec_p)
    cldfrc = Ptr[float](cldfrc_p)
    qs = Ptr[float](qs_p)
    fsnow_conv = Ptr[float](fsnow_conv_p)
    prec = Ptr[float](prec_p)
    snow = Ptr[float](snow_p)
    tend_s = Ptr[float](tend_s_p)
    tend_q = Ptr[float](tend_q_p)
    tend_s_snwprd = Ptr[float](tend_s_snwprd_p)
    tend_s_snwevmlt = Ptr[float](tend_s_snwevmlt_p)
    evpstore = Ptr[float](evpstore_p)
    substore = Ptr[float](substore_p)
    ntprprd = Ptr[float](ntprprd_p)
    ntsnprd = Ptr[float](ntsnprd_p)
    flxprec = Ptr[float](flxprec_p)
    flxsnow = Ptr[float](flxsnow_p)
    evpvint = Ptr[float](evpvint_p)

    i = 0
    while i < ncol:
        prec[i] = prec[i] * 1000.0
        flxprec[i] = 0.0
        flxsnow[i] = 0.0
        evpvint[i] = 0.0
        i += 1

    k = 0
    while k < pver:
        i = 0
        while i < pcols:
            idx = i + k * pcols
            evpstore[idx] = 0.0
            substore[idx] = 0.0
            i += 1
        k += 1

    k = 0
    while k < pver:
        i = 0
        while i < ncol:
            idx = i + k * pcols
            flx_idx = i + k * pcols
            flx_next_idx = i + (k + 1) * pcols

            if t[idx] > tmelt:
                flxsntm = 0.0
                snowmlt = flxsnow[flx_idx] * gravit / pdel[idx]
            else:
                flxsntm = flxsnow[flx_idx]
                snowmlt = 0.0

            evplimit = 1.0 - q[idx] / qs[idx]
            if evplimit < 0.0:
                evplimit = 0.0

            if do_org != 0:
                kemask = ke * (1.0 - landfrac[i]) + ke_lnd * landfrac[i]
            else:
                kemask = ke

            evpprec = kemask * (1.0 - cldfrc[idx]) * evplimit * sqrt(flxprec[flx_idx])

            evplimit = (qs[idx] - q[idx]) / deltat
            if evplimit < 0.0:
                evplimit = 0.0

            limit2 = flxprec[flx_idx] * gravit / pdel[idx]
            if limit2 < evplimit:
                evplimit = limit2

            limit2 = (prec[i] - evpvint[i]) * gravit / pdel[idx]
            if limit2 < evplimit:
                evplimit = limit2

            if evplimit < evpprec:
                evpprec = evplimit

            if flxprec[flx_idx] > 0.0:
                work1 = flxsntm / flxprec[flx_idx]
                if work1 < 0.0:
                    work1 = 0.0
                if work1 > 1.0:
                    work1 = 1.0
                evpsnow = evpprec * work1
            else:
                evpsnow = 0.0

            evpstore[idx] = evpprec
            substore[idx] = evpsnow

            evpvint[i] = evpvint[i] + evpprec * pdel[idx] / gravit
            ntprprd[idx] = prdprec[idx] - evpprec

            if pergro != 0:
                work1 = flxsnow[flx_idx] / (flxprec[flx_idx] + 8.64e-11)
                if work1 < 0.0:
                    work1 = 0.0
                if work1 > 1.0:
                    work1 = 1.0
            else:
                if flxprec[flx_idx] > 0.0:
                    work1 = flxsnow[flx_idx] / flxprec[flx_idx]
                    if work1 < 0.0:
                        work1 = 0.0
                    if work1 > 1.0:
                        work1 = 1.0
                else:
                    work1 = 0.0

            if fsnow_conv[idx] > work1:
                work2 = fsnow_conv[idx]
            else:
                work2 = work1
            if snowmlt > 0.0:
                work2 = 0.0

            ntsnprd[idx] = prdprec[idx] * work2 - evpsnow - snowmlt
            tend_s_snwprd[idx] = prdprec[idx] * work2 * latice
            tend_s_snwevmlt[idx] = -(evpsnow + snowmlt) * latice

            flxprec[flx_next_idx] = flxprec[flx_idx] + ntprprd[idx] * pdel[idx] / gravit
            flxsnow[flx_next_idx] = flxsnow[flx_idx] + ntsnprd[idx] * pdel[idx] / gravit

            if flxprec[flx_next_idx] < 0.0:
                flxprec[flx_next_idx] = 0.0
            if flxsnow[flx_next_idx] < 0.0:
                flxsnow[flx_next_idx] = 0.0

            tend_s[idx] = -evpprec * latvap + ntsnprd[idx] * latice
            tend_q[idx] = evpprec
            i += 1
        k += 1

    i = 0
    while i < ncol:
        bottom_idx = i + pver * pcols
        prec[i] = flxprec[bottom_idx] / 1000.0
        snow[i] = flxsnow[bottom_idx] / 1000.0
        i += 1

@export
def zm_conv_evap_main_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    pergro: int,
    do_org: int,
    gravit: float,
    latvap: float,
    latice: float,
    tmelt: float,
    ke: float,
    ke_lnd: float,
    deltat: float,
    t_p: cobj,
    q_p: cobj,
    pdel_p: cobj,
    landfrac_p: cobj,
    prdprec_p: cobj,
    cldfrc_p: cobj,
    qs_p: cobj,
    fsnow_conv_p: cobj,
    prec_p: cobj,
    snow_p: cobj,
    tend_s_p: cobj,
    tend_q_p: cobj,
    tend_s_snwprd_p: cobj,
    tend_s_snwevmlt_p: cobj,
    evpstore_p: cobj,
    substore_p: cobj,
    ntprprd_p: cobj,
    ntsnprd_p: cobj,
    flxprec_p: cobj,
    flxsnow_p: cobj,
    evpvint_p: cobj,
):
    zm_conv_evap_main_stage_dispatch_codon(
        ncol,
        pcols,
        pver,
        pverp,
        pergro,
        do_org,
        gravit,
        latvap,
        latice,
        tmelt,
        ke,
        ke_lnd,
        deltat,
        t_p,
        q_p,
        pdel_p,
        landfrac_p,
        prdprec_p,
        cldfrc_p,
        qs_p,
        fsnow_conv_p,
        prec_p,
        snow_p,
        tend_s_p,
        tend_q_p,
        tend_s_snwprd_p,
        tend_s_snwevmlt_p,
        evpstore_p,
        substore_p,
        ntprprd_p,
        ntsnprd_p,
        flxprec_p,
        flxsnow_p,
        evpvint_p,
    )


@export
def zm_momtran_prep_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    state_u_p: cobj,
    state_v_p: cobj,
    winds_p: cobj,
):
    state_u = Ptr[float](state_u_p)
    state_v = Ptr[float](state_v_p)
    winds = Ptr[float](winds_p)
    plane = pcols * pver

    k = 0
    while k < pver:
        i = 0
        while i < ncol:
            idx = i + k * pcols
            winds[idx] = state_u[idx]
            winds[idx + plane] = state_v[idx]
            i += 1
        k += 1


@export
def zm_momtran_post_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    cpair: float,
    wind_tends_p: cobj,
    seten_p: cobj,
    ptend_u_p: cobj,
    ptend_v_p: cobj,
    ptend_s_p: cobj,
    ftem_p: cobj,
):
    wind_tends = Ptr[float](wind_tends_p)
    seten = Ptr[float](seten_p)
    ptend_u = Ptr[float](ptend_u_p)
    ptend_v = Ptr[float](ptend_v_p)
    ptend_s = Ptr[float](ptend_s_p)
    ftem = Ptr[float](ftem_p)
    plane = pcols * pver

    k = 0
    while k < pver:
        i = 0
        while i < ncol:
            idx = i + k * pcols
            ptend_u[idx] = wind_tends[idx]
            ptend_v[idx] = wind_tends[idx + plane]
            ptend_s[idx] = seten[idx]
            ftem[idx] = seten[idx] / cpair
            i += 1
        k += 1


@export
def zm_momtran_main_stage_dispatch_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    ncnst: int,
    il1g: int,
    il2g: int,
    dt: float,
    domomtran_p: cobj,
    q_p: cobj,
    mu_p: cobj,
    md_p: cobj,
    du_p: cobj,
    eu_p: cobj,
    ed_p: cobj,
    dp_p: cobj,
    dsubcld_p: cobj,
    jt_p: cobj,
    mx_p: cobj,
    ideep_p: cobj,
    dqdt_p: cobj,
    pguall_p: cobj,
    pgdall_p: cobj,
    icwu_p: cobj,
    icwd_p: cobj,
    seten_p: cobj,
    chat_p: cobj,
    cond_p: cobj,
    const_p: cobj,
    conu_p: cobj,
    dcondt_p: cobj,
    mududp_p: cobj,
    mddudp_p: cobj,
    pgu_p: cobj,
    pgd_p: cobj,
    gseten_p: cobj,
    mflux_p: cobj,
    wind0_p: cobj,
    windf_p: cobj,
):
    domomtran = Ptr[int](domomtran_p)
    q = Ptr[float](q_p)
    mu = Ptr[float](mu_p)
    md = Ptr[float](md_p)
    du = Ptr[float](du_p)
    eu = Ptr[float](eu_p)
    ed = Ptr[float](ed_p)
    dp = Ptr[float](dp_p)
    dsubcld = Ptr[float](dsubcld_p)
    jt = Ptr[int](jt_p)
    mx = Ptr[int](mx_p)
    ideep = Ptr[int](ideep_p)
    dqdt = Ptr[float](dqdt_p)
    pguall = Ptr[float](pguall_p)
    pgdall = Ptr[float](pgdall_p)
    icwu = Ptr[float](icwu_p)
    icwd = Ptr[float](icwd_p)
    seten = Ptr[float](seten_p)
    chat = Ptr[float](chat_p)
    cond = Ptr[float](cond_p)
    const = Ptr[float](const_p)
    conu = Ptr[float](conu_p)
    dcondt = Ptr[float](dcondt_p)
    mududp = Ptr[float](mududp_p)
    mddudp = Ptr[float](mddudp_p)
    pgu = Ptr[float](pgu_p)
    pgd = Ptr[float](pgd_p)
    gseten = Ptr[float](gseten_p)
    mflux = Ptr[float](mflux_p)
    wind0 = Ptr[float](wind0_p)
    windf = Ptr[float](windf_p)
    plane = pcols * pver
    flux_plane = pcols * pverp
    ilo = il1g - 1
    ihi = il2g - 1

    m = 0
    while m < ncnst:
        moff = m * plane
        foff = m * flux_plane
        k = 0
        while k < pver:
            i = 0
            while i < pcols:
                idx = i + k * pcols
                pguall[idx + moff] = 0.0
                pgdall[idx + moff] = 0.0
                wind0[idx + moff] = 0.0
                windf[idx + moff] = 0.0
                if i < ncol:
                    icwu[idx + moff] = q[idx + moff]
                    icwd[idx + moff] = q[idx + moff]
                i += 1
            k += 1
        k = 0
        while k < pverp:
            i = 0
            while i < pcols:
                mflux[i + k * pcols + foff] = 0.0
                i += 1
            k += 1
        m += 1

    k = 0
    while k < pver:
        i = 0
        while i < pcols:
            idx = i + k * pcols
            seten[idx] = 0.0
            gseten[idx] = 0.0
            i += 1
        k += 1

    momcu = 0.4
    momcd = 0.4
    mbsth = 1.0e-15

    ktm = pver
    kbm = pver
    i = ilo
    while i <= ihi:
        if jt[i] < ktm:
            ktm = jt[i]
        if mx[i] < kbm:
            kbm = mx[i]
        i += 1

    m = 0
    while m < ncnst:
        moff = m * plane
        foff = m * flux_plane
        if domomtran[m] != 0:
            k = 0
            while k < pver:
                i = ilo
                while i <= ihi:
                    ii = ideep[i] - 1
                    idx = i + k * pcols
                    const[idx] = q[ii + k * pcols + moff]
                    wind0[idx + moff] = const[idx]
                    i += 1
                k += 1

            k = 0
            while k < pver:
                km1 = k - 1
                if km1 < 0:
                    km1 = 0
                i = ilo
                while i <= ihi:
                    idx = i + k * pcols
                    chat[idx] = 0.5 * (const[idx] + const[i + km1 * pcols])
                    conu[idx] = chat[idx]
                    cond[idx] = chat[idx]
                    dcondt[idx] = 0.0
                    i += 1
                k += 1

            i = 0
            while i < il2g:
                pgu[i] = 0.0
                pgd[i] = 0.0
                i += 1

            k = 1
            while k < pver - 1:
                km1 = k - 1
                kp1 = k + 1
                i = ilo
                while i <= ihi:
                    idx = i + k * pcols
                    mududp[idx] = (
                        mu[idx] * (const[idx] - const[i + km1 * pcols]) / dp[i + km1 * pcols]
                        + mu[i + kp1 * pcols] * (const[i + kp1 * pcols] - const[idx]) / dp[idx]
                    )
                    pgu[idx] = -momcu * 0.5 * mududp[idx]
                    mddudp[idx] = (
                        md[idx] * (const[idx] - const[i + km1 * pcols]) / dp[i + km1 * pcols]
                        + md[i + kp1 * pcols] * (const[i + kp1 * pcols] - const[idx]) / dp[idx]
                    )
                    pgd[idx] = -momcd * 0.5 * mddudp[idx]
                    i += 1
                k += 1

            k = pver - 1
            km1 = k - 1
            i = ilo
            while i <= ihi:
                idx = i + k * pcols
                mududp[idx] = mu[idx] * (const[idx] - const[i + km1 * pcols]) / dp[i + km1 * pcols]
                pgu[idx] = -momcu * mududp[idx]
                mddudp[idx] = md[idx] * (const[idx] - const[i + km1 * pcols]) / dp[i + km1 * pcols]
                pgd[idx] = -momcd * mddudp[idx]
                i += 1

            k = 1
            km1 = 0
            kk = pver - 1
            i = ilo
            while i <= ihi:
                kkidx = i + kk * pcols
                mupdudp = mu[kkidx] + du[kkidx] * dp[kkidx]
                if mupdudp > mbsth:
                    conu[kkidx] = (eu[kkidx] * const[kkidx] * dp[kkidx] + pgu[kkidx] * dp[kkidx]) / mupdudp
                idx = i + k * pcols
                if md[idx] < -mbsth:
                    cond[idx] = (-ed[i + km1 * pcols] * const[i + km1 * pcols] * dp[i + km1 * pcols]) - pgd[
                        i + km1 * pcols
                    ] * dp[i + km1 * pcols] / md[idx]
                i += 1

            kk = pver - 2
            while kk >= 0:
                kkp1 = kk + 1
                i = ilo
                while i <= ihi:
                    idx = i + kk * pcols
                    mupdudp = mu[idx] + du[idx] * dp[idx]
                    if mupdudp > mbsth:
                        conu[idx] = (
                            mu[i + kkp1 * pcols] * conu[i + kkp1 * pcols]
                            + eu[idx] * const[idx] * dp[idx]
                            + pgu[idx] * dp[idx]
                        ) / mupdudp
                    i += 1
                kk -= 1

            k = 2
            while k < pver:
                km1 = k - 1
                i = ilo
                while i <= ihi:
                    idx = i + k * pcols
                    if md[idx] < -mbsth:
                        cond[idx] = (
                            md[i + km1 * pcols] * cond[i + km1 * pcols]
                            - ed[i + km1 * pcols] * const[i + km1 * pcols] * dp[i + km1 * pcols]
                            - pgd[i + km1 * pcols] * dp[i + km1 * pcols]
                        ) / md[idx]
                    i += 1
                k += 1

            k = ktm - 1
            while k < pver:
                kp1 = k + 1
                if kp1 >= pver:
                    kp1 = pver - 1
                i = ilo
                while i <= ihi:
                    idx = i + k * pcols
                    dcondt[idx] = (
                        +(
                            mu[i + kp1 * pcols] * (conu[i + kp1 * pcols] - chat[i + kp1 * pcols])
                            - mu[idx] * (conu[idx] - chat[idx])
                            + md[i + kp1 * pcols] * (cond[i + kp1 * pcols] - chat[i + kp1 * pcols])
                            - md[idx] * (cond[idx] - chat[idx])
                        )
                        / dp[idx]
                    )
                    i += 1
                k += 1

            k = kbm - 1
            while k < pver:
                i = ilo
                while i <= ihi:
                    if k + 1 == mx[i]:
                        idx = i + k * pcols
                        dcondt[idx] = (1.0 / dp[idx]) * (
                            -mu[idx] * (conu[idx] - chat[idx]) - md[idx] * (cond[idx] - chat[idx])
                        )
                    i += 1
                k += 1

            k = 0
            while k < pver:
                i = 0
                while i < pcols:
                    dqdt[i + k * pcols + moff] = 0.0
                    i += 1
                k += 1

            k = 0
            while k < pver:
                i = ilo
                while i <= ihi:
                    ii = ideep[i] - 1
                    src_idx = i + k * pcols
                    dst_idx = ii + k * pcols + moff
                    dqdt[dst_idx] = dcondt[src_idx]
                    pguall[dst_idx] = -pgu[src_idx]
                    pgdall[dst_idx] = -pgd[src_idx]
                    icwu[dst_idx] = conu[src_idx]
                    icwd[dst_idx] = cond[src_idx]
                    i += 1
                k += 1

            k = ktm - 1
            while k < pver:
                i = ilo
                while i <= ihi:
                    idx = i + k * pcols
                    mflux[i + k * pcols + foff] = -mu[idx] * (conu[idx] - chat[idx]) - md[idx] * (
                        cond[idx] - chat[idx]
                    )
                    i += 1
                k += 1

            k = ktm - 1
            while k < pver:
                i = ilo
                while i <= ihi:
                    idx = i + k * pcols
                    windf[idx + moff] = const[idx] - (
                        mflux[i + (k + 1) * pcols + foff] - mflux[i + k * pcols + foff]
                    ) * dt / dp[idx]
                    i += 1
                k += 1
        m += 1

    k = ktm - 1
    while k < pver:
        km1 = k - 1
        if km1 < 0:
            km1 = 0
        kp1 = k + 1
        if kp1 >= pver:
            kp1 = pver - 1
        i = ilo
        while i <= ihi:
            idx = i + k * pcols
            utop = (wind0[idx] + wind0[i + km1 * pcols]) / 2.0
            vtop = (wind0[idx + plane] + wind0[i + km1 * pcols + plane]) / 2.0
            ubot = (wind0[i + kp1 * pcols] + wind0[idx]) / 2.0
            vbot = (wind0[i + kp1 * pcols + plane] + wind0[idx + plane]) / 2.0
            fket = utop * mflux[i + k * pcols] + vtop * mflux[i + k * pcols + flux_plane]
            fkeb = ubot * mflux[i + (k + 1) * pcols] + vbot * mflux[i + (k + 1) * pcols + flux_plane]
            ketend_cons = (fket - fkeb) / dp[idx]
            ketend = ((windf[idx] ** 2 + windf[idx + plane] ** 2) - (wind0[idx] ** 2 + wind0[idx + plane] ** 2)) * 0.5 / dt
            gseten[idx] = ketend_cons - ketend
            i += 1
        k += 1

    k = 0
    while k < pver:
        i = ilo
        while i <= ihi:
            ii = ideep[i] - 1
            seten[ii + k * pcols] = gseten[i + k * pcols]
            i += 1
        k += 1

@export
def zm_momtran_main_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    ncnst: int,
    il1g: int,
    il2g: int,
    dt: float,
    domomtran_p: cobj,
    q_p: cobj,
    mu_p: cobj,
    md_p: cobj,
    du_p: cobj,
    eu_p: cobj,
    ed_p: cobj,
    dp_p: cobj,
    dsubcld_p: cobj,
    jt_p: cobj,
    mx_p: cobj,
    ideep_p: cobj,
    dqdt_p: cobj,
    pguall_p: cobj,
    pgdall_p: cobj,
    icwu_p: cobj,
    icwd_p: cobj,
    seten_p: cobj,
    chat_p: cobj,
    cond_p: cobj,
    const_p: cobj,
    conu_p: cobj,
    dcondt_p: cobj,
    mududp_p: cobj,
    mddudp_p: cobj,
    pgu_p: cobj,
    pgd_p: cobj,
    gseten_p: cobj,
    mflux_p: cobj,
    wind0_p: cobj,
    windf_p: cobj,
):
    zm_momtran_main_stage_dispatch_codon(
        ncol,
        pcols,
        pver,
        pverp,
        ncnst,
        il1g,
        il2g,
        dt,
        domomtran_p,
        q_p,
        mu_p,
        md_p,
        du_p,
        eu_p,
        ed_p,
        dp_p,
        dsubcld_p,
        jt_p,
        mx_p,
        ideep_p,
        dqdt_p,
        pguall_p,
        pgdall_p,
        icwu_p,
        icwd_p,
        seten_p,
        chat_p,
        cond_p,
        const_p,
        conu_p,
        dcondt_p,
        mududp_p,
        mddudp_p,
        pgu_p,
        pgd_p,
        gseten_p,
        mflux_p,
        wind0_p,
        windf_p,
    )


@export
def zm_convtran_main_stage_dispatch_codon(
    pcols: int,
    pver: int,
    ncnst: int,
    il1g: int,
    il2g: int,
    doconvtran_p: cobj,
    is_dry_p: cobj,
    q_p: cobj,
    mu_p: cobj,
    md_p: cobj,
    du_p: cobj,
    eu_p: cobj,
    ed_p: cobj,
    dp_p: cobj,
    fracis_p: cobj,
    dpdry_p: cobj,
    jt_p: cobj,
    mx_p: cobj,
    ideep_p: cobj,
    dqdt_p: cobj,
    chat_p: cobj,
    cond_p: cobj,
    const_p: cobj,
    fisg_p: cobj,
    conu_p: cobj,
    dcondt_p: cobj,
    dutmp_p: cobj,
    eutmp_p: cobj,
    edtmp_p: cobj,
    dptmp_p: cobj,
    trace_water: int,
    nwt_liq: int,
    nwt_ice: int,
    wtrc_qmin: float,
    wtrc_liq_iatype_p: cobj,
    wtrc_ice_iatype_p: cobj,
    iwspec_p: cobj,
    rstd_p: cobj,
    Rwt_p: cobj,
):
    doconvtran = Ptr[int](doconvtran_p)
    is_dry = Ptr[int](is_dry_p)
    q = Ptr[float](q_p)
    mu = Ptr[float](mu_p)
    md = Ptr[float](md_p)
    du = Ptr[float](du_p)
    eu = Ptr[float](eu_p)
    ed = Ptr[float](ed_p)
    dp = Ptr[float](dp_p)
    fracis = Ptr[float](fracis_p)
    dpdry = Ptr[float](dpdry_p)
    jt = Ptr[int](jt_p)
    mx = Ptr[int](mx_p)
    ideep = Ptr[int](ideep_p)
    dqdt = Ptr[float](dqdt_p)
    chat = Ptr[float](chat_p)
    cond = Ptr[float](cond_p)
    const = Ptr[float](const_p)
    fisg = Ptr[float](fisg_p)
    conu = Ptr[float](conu_p)
    dcondt = Ptr[float](dcondt_p)
    dutmp = Ptr[float](dutmp_p)
    eutmp = Ptr[float](eutmp_p)
    edtmp = Ptr[float](edtmp_p)
    dptmp = Ptr[float](dptmp_p)
    wtrc_liq_iatype = Ptr[int](wtrc_liq_iatype_p)
    wtrc_ice_iatype = Ptr[int](wtrc_ice_iatype_p)
    iwspec = Ptr[int](iwspec_p)
    rstd = Ptr[float](rstd_p)
    Rwt = Ptr[float](Rwt_p)
    plane = pcols * pver
    ilo = il1g - 1
    ihi = il2g - 1
    small = 1.0e-36
    mbsth = 1.0e-15

    ktm = pver
    kbm = pver
    i = ilo
    while i <= ihi:
        if jt[i] < ktm:
            ktm = jt[i]
        if mx[i] < kbm:
            kbm = mx[i]
        i += 1

    m = 1
    while m < ncnst:
        moff = m * plane
        if doconvtran[m] != 0:
            if is_dry[m] != 0:
                k = 0
                while k < pver:
                    i = ilo
                    while i <= ihi:
                        idx = i + k * pcols
                        dptmp[idx] = dpdry[idx]
                        dutmp[idx] = du[idx] * dp[idx] / dpdry[idx]
                        eutmp[idx] = eu[idx] * dp[idx] / dpdry[idx]
                        edtmp[idx] = ed[idx] * dp[idx] / dpdry[idx]
                        i += 1
                    k += 1
            else:
                k = 0
                while k < pver:
                    i = ilo
                    while i <= ihi:
                        idx = i + k * pcols
                        dptmp[idx] = dp[idx]
                        dutmp[idx] = du[idx]
                        eutmp[idx] = eu[idx]
                        edtmp[idx] = ed[idx]
                        i += 1
                    k += 1

            k = 0
            while k < pver:
                i = ilo
                while i <= ihi:
                    ii = ideep[i] - 1
                    idx = i + k * pcols
                    const[idx] = q[ii + k * pcols + moff]
                    fisg[idx] = fracis[ii + k * pcols + moff]
                    i += 1
                k += 1

            k = 0
            while k < pver:
                km1 = k - 1
                if km1 < 0:
                    km1 = 0
                i = ilo
                while i <= ihi:
                    idx = i + k * pcols
                    km1idx = i + km1 * pcols
                    minc = min(const[km1idx], const[idx])
                    maxc = max(const[km1idx], const[idx])
                    if minc < 0.0:
                        cdifr = 0.0
                    else:
                        cdifr = abs(const[idx] - const[km1idx]) / max(maxc, small)
                    if cdifr > 1.0e-6:
                        cabv = max(const[km1idx], maxc * 1.0e-12)
                        cbel = max(const[idx], maxc * 1.0e-12)
                        chat[idx] = log(cabv / cbel) / (cabv - cbel) * cabv * cbel
                    else:
                        chat[idx] = 0.5 * (const[idx] + const[km1idx])
                    conu[idx] = chat[idx]
                    cond[idx] = chat[idx]
                    dcondt[idx] = 0.0
                    i += 1
                k += 1

            k = 1
            km1 = 0
            kk = pver - 1
            i = ilo
            while i <= ihi:
                kkidx = i + kk * pcols
                mupdudp = mu[kkidx] + dutmp[kkidx] * dptmp[kkidx]
                if mupdudp > mbsth:
                    conu[kkidx] = (+eutmp[kkidx] * fisg[kkidx] * const[kkidx] * dptmp[kkidx]) / mupdudp
                idx = i + k * pcols
                km1idx = i + km1 * pcols
                if md[idx] < -mbsth:
                    cond[idx] = (-edtmp[km1idx] * fisg[km1idx] * const[km1idx] * dptmp[km1idx]) / md[idx]
                i += 1

            kk = pver - 2
            while kk >= 0:
                kkp1 = kk + 1
                i = ilo
                while i <= ihi:
                    idx = i + kk * pcols
                    mupdudp = mu[idx] + dutmp[idx] * dptmp[idx]
                    if mupdudp > mbsth:
                        conu[idx] = (
                            mu[i + kkp1 * pcols] * conu[i + kkp1 * pcols]
                            + eutmp[idx] * fisg[idx] * const[idx] * dptmp[idx]
                        ) / mupdudp
                    i += 1
                kk -= 1

            k = 2
            while k < pver:
                km1 = k - 1
                i = ilo
                while i <= ihi:
                    idx = i + k * pcols
                    km1idx = i + km1 * pcols
                    if md[idx] < -mbsth:
                        cond[idx] = (
                            md[km1idx] * cond[km1idx]
                            - edtmp[km1idx] * fisg[km1idx] * const[km1idx] * dptmp[km1idx]
                        ) / md[idx]
                    i += 1
                k += 1

            k = ktm - 1
            while k < pver:
                km1 = k - 1
                if km1 < 0:
                    km1 = 0
                kp1 = k + 1
                if kp1 >= pver:
                    kp1 = pver - 1
                i = ilo
                while i <= ihi:
                    idx = i + k * pcols
                    km1idx = i + km1 * pcols
                    kp1idx = i + kp1 * pcols
                    fluxin = mu[kp1idx] * conu[kp1idx] + mu[idx] * min(chat[idx], const[km1idx]) - (
                        md[idx] * cond[idx] + md[kp1idx] * min(chat[kp1idx], const[kp1idx])
                    )
                    fluxout = mu[idx] * conu[idx] + mu[kp1idx] * min(chat[kp1idx], const[idx]) - (
                        md[kp1idx] * cond[kp1idx] + md[idx] * min(chat[idx], const[idx])
                    )
                    netflux = fluxin - fluxout
                    if abs(netflux) < max(fluxin, fluxout) * 1.0e-12:
                        netflux = 0.0
                    dcondt[idx] = netflux / dptmp[idx]
                    i += 1
                k += 1

            k = kbm - 1
            while k < pver:
                km1 = k - 1
                if km1 < 0:
                    km1 = 0
                i = ilo
                while i <= ihi:
                    idx = i + k * pcols
                    km1idx = i + km1 * pcols
                    if k + 1 == mx[i]:
                        fluxin = mu[idx] * min(chat[idx], const[km1idx]) - md[idx] * cond[idx]
                        fluxout = mu[idx] * conu[idx] - md[idx] * min(chat[idx], const[idx])
                        netflux = fluxin - fluxout
                        if abs(netflux) < max(fluxin, fluxout) * 1.0e-12:
                            netflux = 0.0
                        dcondt[idx] = netflux / dptmp[idx]
                    elif k + 1 > mx[i]:
                        dcondt[idx] = 0.0
                    i += 1
                k += 1

            k = 0
            while k < pver:
                i = 0
                while i < pcols:
                    dqdt[i + k * pcols + moff] = 0.0
                    i += 1
                k += 1

            k = 0
            while k < pver:
                i = ilo
                while i <= ihi:
                    ii = ideep[i] - 1
                    dqdt[ii + k * pcols + moff] = dcondt[i + k * pcols]
                    i += 1
                k += 1
        m += 1

    if trace_water != 0:
        n = 0
        while n < plane * nwt_ice * 2:
            Rwt[n] = 1.0
            n += 1

        base_liq = wtrc_liq_iatype[0] - 1
        if doconvtran[base_liq] != 0:
            base_ice = wtrc_ice_iatype[0] - 1
            m = 1
            while m < nwt_liq:
                liq_trc = wtrc_liq_iatype[m] - 1
                ice_trc = wtrc_ice_iatype[m] - 1
                liq_ispec = iwspec[liq_trc]
                ice_ispec = iwspec[ice_trc]
                k = 0
                while k < pver:
                    i = ilo
                    while i <= ihi:
                        ii = ideep[i] - 1
                        idx = ii + k * pcols
                        liq_base_val = dqdt[idx + base_liq * plane]
                        ice_base_val = dqdt[idx + base_ice * plane]
                        liq_ratio_idx = idx + m * plane
                        ice_ratio_idx = idx + m * plane + nwt_ice * plane
                        if abs(liq_base_val) < wtrc_qmin:
                            Rwt[liq_ratio_idx] = rstd[liq_ispec - 1]
                        else:
                            Rwt[liq_ratio_idx] = dqdt[idx + liq_trc * plane] / liq_base_val
                        if abs(ice_base_val) < wtrc_qmin:
                            Rwt[ice_ratio_idx] = rstd[ice_ispec - 1]
                        else:
                            Rwt[ice_ratio_idx] = dqdt[idx + ice_trc * plane] / ice_base_val
                        i += 1
                    k += 1
                m += 1

@export
def zm_convtran_main_codon(
    pcols: int,
    pver: int,
    ncnst: int,
    il1g: int,
    il2g: int,
    doconvtran_p: cobj,
    is_dry_p: cobj,
    q_p: cobj,
    mu_p: cobj,
    md_p: cobj,
    du_p: cobj,
    eu_p: cobj,
    ed_p: cobj,
    dp_p: cobj,
    fracis_p: cobj,
    dpdry_p: cobj,
    jt_p: cobj,
    mx_p: cobj,
    ideep_p: cobj,
    dqdt_p: cobj,
    chat_p: cobj,
    cond_p: cobj,
    const_p: cobj,
    fisg_p: cobj,
    conu_p: cobj,
    dcondt_p: cobj,
    dutmp_p: cobj,
    eutmp_p: cobj,
    edtmp_p: cobj,
    dptmp_p: cobj,
    trace_water: int,
    nwt_liq: int,
    nwt_ice: int,
    wtrc_qmin: float,
    wtrc_liq_iatype_p: cobj,
    wtrc_ice_iatype_p: cobj,
    iwspec_p: cobj,
    rstd_p: cobj,
    Rwt_p: cobj,
):
    zm_convtran_main_stage_dispatch_codon(
        pcols,
        pver,
        ncnst,
        il1g,
        il2g,
        doconvtran_p,
        is_dry_p,
        q_p,
        mu_p,
        md_p,
        du_p,
        eu_p,
        ed_p,
        dp_p,
        fracis_p,
        dpdry_p,
        jt_p,
        mx_p,
        ideep_p,
        dqdt_p,
        chat_p,
        cond_p,
        const_p,
        fisg_p,
        conu_p,
        dcondt_p,
        dutmp_p,
        eutmp_p,
        edtmp_p,
        dptmp_p,
        trace_water,
        nwt_liq,
        nwt_ice,
        wtrc_qmin,
        wtrc_liq_iatype_p,
        wtrc_ice_iatype_p,
        iwspec_p,
        rstd_p,
        Rwt_p,
    )


@export
def zm_convtran1_prep_shell_codon(
    pcols: int,
    pver: int,
    fake_dpdry_p: cobj,
):
    fake_dpdry = Ptr[float](fake_dpdry_p)

    k = 0
    while k < pver:
        i = 0
        while i < pcols:
            fake_dpdry[i + k * pcols] = 0.0
            i += 1
        k += 1


@export
def zm_convtran1_ratio_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    wtrc_nwset: int,
    Rwt_p: cobj,
    ptend_q_p: cobj,
    liq_type_p: cobj,
    ice_type_p: cobj,
):
    Rwt = Ptr[float](Rwt_p)
    ptend_q = Ptr[float](ptend_q_p)
    liq_type = Ptr[int](liq_type_p)
    ice_type = Ptr[int](ice_type_p)
    plane = pcols * pver

    m = 1
    while m < wtrc_nwset:
        liq_idx = (liq_type[m] - 1) * plane
        ice_idx = (ice_type[m] - 1) * plane
        base_liq_idx = (liq_type[0] - 1) * plane
        base_ice_idx = (ice_type[0] - 1) * plane
        ratio_m_offset = m * plane
        ratio_ice_offset = wtrc_nwset * plane + m * plane
        k = 0
        while k < pver:
            i = 0
            while i < ncol:
                idx = i + k * pcols
                ptend_q[idx + liq_idx] = Rwt[idx + ratio_m_offset] * ptend_q[idx + base_liq_idx]
                ptend_q[idx + ice_idx] = Rwt[idx + ratio_ice_offset] * ptend_q[idx + base_ice_idx]
                i += 1
            k += 1
        m += 1


@export
def zm_convtran2_dpdry_shell_codon(
    pcols: int,
    pver: int,
    lengath: int,
    state_pdeldry_p: cobj,
    ideep_p: cobj,
    dpdry_p: cobj,
):
    state_pdeldry = Ptr[float](state_pdeldry_p)
    ideep = Ptr[int](ideep_p)
    dpdry = Ptr[float](dpdry_p)

    k = 0
    while k < pver:
        i = 0
        while i < pcols:
            dpdry[i + k * pcols] = 0.0
            i += 1
        k += 1

    i = 0
    while i < lengath:
        src_i = ideep[i] - 1
        k = 0
        while k < pver:
            dpdry[i + k * pcols] = state_pdeldry[src_i + k * pcols] / 100.0
            k += 1
        i += 1


@export
def zm_parcel_dilute_codon(
    lchnk: int,
    ncol: int,
    msg: int,
    pcols: int,
    pver: int,
    zm_org: int,
    grav: float,
    rgas: float,
    cpliq: float,
    tfreez: float,
    latice: float,
    rl: float,
    cpwv: float,
    cpres: float,
    eps1: float,
    rh2o: float,
    epsilo: float,
    omeps: float,
    wv_idx: int,
    klaunch_p: cobj,
    p_p: cobj,
    t_p: cobj,
    q_p: cobj,
    tpert_p: cobj,
    tp_p: cobj,
    tpv_p: cobj,
    qstp_p: cobj,
    pl_p: cobj,
    tl_p: cobj,
    lcl_p: cobj,
    org_p: cobj,
    landfrac_p: cobj,
    tmix_p: cobj,
    qtmix_p: cobj,
    qsmix_p: cobj,
    smix_p: cobj,
    xsh2o_p: cobj,
    ds_xsh2o_p: cobj,
    ds_freeze_p: cobj,
    mp_p: cobj,
    qtp_p: cobj,
    sp_p: cobj,
    sp0_p: cobj,
    qtp0_p: cobj,
    mp0_p: cobj,
    status_p: cobj,
):
    klaunch = Ptr[i32](klaunch_p)
    p = Ptr[float](p_p)
    t = Ptr[float](t_p)
    q = Ptr[float](q_p)
    tpert = Ptr[float](tpert_p)
    tp = Ptr[float](tp_p)
    tpv = Ptr[float](tpv_p)
    qstp = Ptr[float](qstp_p)
    pl = Ptr[float](pl_p)
    tl = Ptr[float](tl_p)
    lcl = Ptr[i32](lcl_p)
    org = Ptr[float](org_p)
    landfrac = Ptr[float](landfrac_p)
    tmix = Ptr[float](tmix_p)
    qtmix = Ptr[float](qtmix_p)
    qsmix = Ptr[float](qsmix_p)
    smix = Ptr[float](smix_p)
    xsh2o = Ptr[float](xsh2o_p)
    ds_xsh2o = Ptr[float](ds_xsh2o_p)
    ds_freeze = Ptr[float](ds_freeze_p)
    mp = Ptr[float](mp_p)
    qtp = Ptr[float](qtp_p)
    sp = Ptr[float](sp_p)
    sp0 = Ptr[float](sp0_p)
    qtp0 = Ptr[float](qtp0_p)
    mp0 = Ptr[float](mp0_p)
    status = Ptr[int](status_p)

    status[0] = 1

    org2rkm = 0.0
    org2tpert = 0.0
    if zm_org != 0:
        org2rkm = 10.0
        org2tpert = 0.0

    nit_lheat = 2
    dmpdz = -1.0e-3
    dmpdz_lnd = -1.0e-3
    lwmax = 1.0e-3
    tscool = 0.0

    kk0 = 0
    while kk0 < pver:
        i0 = 0
        while i0 < pcols:
            idx0 = i0 + kk0 * pcols
            qtmix[idx0] = 0.0
            smix[idx0] = 0.0
            i0 += 1
        kk0 += 1

    i0 = 0
    while i0 < pcols:
        qtp0[i0] = 0.0
        sp0[i0] = 0.0
        mp0[i0] = 0.0
        qtp[i0] = 0.0
        sp[i0] = 0.0
        mp[i0] = 0.0
        i0 += 1

    kk = pver
    while kk >= msg + 1:
        k0 = kk - 1
        i = 0
        while i < ncol:
            idx = i + k0 * pcols
            launch_i = int(klaunch[i])

            if kk == launch_i:
                qtp0[i] = q[idx]
                sp0[i] = _zm_entropy_codon(
                    t[idx],
                    p[idx],
                    qtp0[i],
                    rl,
                    cpliq,
                    cpwv,
                    tfreez,
                    cpres,
                    rgas,
                    eps1,
                    rh2o,
                    wv_idx,
                    epsilo,
                    omeps,
                )
                mp0[i] = 1.0
                smix[idx] = sp0[i]
                qtmix[idx] = qtp0[i]
                tmix_val = 0.0
                qsmix_val = 0.0
                conv = 0
                _zm_ientropy_codon(
                    smix[idx],
                    p[idx],
                    qtmix[idx],
                    t[idx],
                    rl,
                    cpliq,
                    cpwv,
                    tfreez,
                    cpres,
                    rgas,
                    eps1,
                    rh2o,
                    wv_idx,
                    epsilo,
                    omeps,
                    __ptr__(tmix_val),
                    __ptr__(qsmix_val),
                    __ptr__(conv),
                )
                if conv == 0:
                    status[0] = 0
                    return
                tmix[idx] = tmix_val
                qsmix[idx] = qsmix_val

            if kk < launch_i:
                idxp1 = i + kk * pcols
                dp = p[idx] - p[idxp1]
                qtenv = 0.5 * (q[idx] + q[idxp1])
                tenv = 0.5 * (t[idx] + t[idxp1])
                penv = 0.5 * (p[idx] + p[idxp1])
                senv = _zm_entropy_codon(
                    tenv,
                    penv,
                    qtenv,
                    rl,
                    cpliq,
                    cpwv,
                    tfreez,
                    cpres,
                    rgas,
                    eps1,
                    rh2o,
                    wv_idx,
                    epsilo,
                    omeps,
                )
                dpdz = -(penv * grav) / (rgas * tenv)
                dzdp = 1.0 / dpdz
                if zm_org != 0:
                    dmpdz_mask = landfrac[i] * dmpdz_lnd + (1.0 - landfrac[i]) * dmpdz
                    dmpdp = (dmpdz_mask / (1.0 + org[idx] * org2rkm)) * dzdp
                else:
                    dmpdp = dmpdz * dzdp

                sp[i] = sp[i] - dmpdp * dp * senv
                qtp[i] = qtp[i] - dmpdp * dp * qtenv
                mp[i] = mp[i] - dmpdp * dp

                smix[idx] = (sp0[i] + sp[i]) / (mp0[i] + mp[i])
                qtmix[idx] = (qtp0[i] + qtp[i]) / (mp0[i] + mp[i])

                tmix_val = 0.0
                qsmix_val = 0.0
                conv = 0
                _zm_ientropy_codon(
                    smix[idx],
                    p[idx],
                    qtmix[idx],
                    tmix[idxp1],
                    rl,
                    cpliq,
                    cpwv,
                    tfreez,
                    cpres,
                    rgas,
                    eps1,
                    rh2o,
                    wv_idx,
                    epsilo,
                    omeps,
                    __ptr__(tmix_val),
                    __ptr__(qsmix_val),
                    __ptr__(conv),
                )
                if conv == 0:
                    status[0] = 0
                    return
                tmix[idx] = tmix_val
                qsmix[idx] = qsmix_val

                if qsmix[idx] <= qtmix[idx] and qsmix[idxp1] > qtmix[idxp1]:
                    lcl[i] = i32(kk)
                    qxsk = qtmix[idx] - qsmix[idx]
                    qxskp1 = qtmix[idxp1] - qsmix[idxp1]
                    dqxsdp = (qxsk - qxskp1) / dp
                    pl[i] = p[idxp1] - qxskp1 / dqxsdp
                    dsdp = (smix[idx] - smix[idxp1]) / dp
                    dqtdp = (qtmix[idx] - qtmix[idxp1]) / dp
                    slcl = smix[idxp1] + dsdp * (pl[i] - p[idxp1])
                    qtlcl = qtmix[idxp1] + dqtdp * (pl[i] - p[idxp1])

                    tl_val = 0.0
                    qslcl = 0.0
                    conv = 0
                    _zm_ientropy_codon(
                        slcl,
                        pl[i],
                        qtlcl,
                        tmix[idx],
                        rl,
                        cpliq,
                        cpwv,
                        tfreez,
                        cpres,
                        rgas,
                        eps1,
                        rh2o,
                        wv_idx,
                        epsilo,
                        omeps,
                        __ptr__(tl_val),
                        __ptr__(qslcl),
                        __ptr__(conv),
                    )
                    if conv == 0:
                        status[0] = 0
                        return
                    tl[i] = tl_val

            i += 1
        kk -= 1

    kk0 = 0
    while kk0 < pver:
        i0 = 0
        while i0 < pcols:
            idx0 = i0 + kk0 * pcols
            xsh2o[idx0] = 0.0
            ds_xsh2o[idx0] = 0.0
            ds_freeze[idx0] = 0.0
            i0 += 1
        kk0 += 1

    kk = pver
    while kk >= msg + 1:
        k0 = kk - 1
        i = 0
        while i < ncol:
            idx = i + k0 * pcols
            launch_i = int(klaunch[i])

            if kk == launch_i:
                tp[idx] = tmix[idx]
                qstp[idx] = q[idx]
                if zm_org != 0:
                    tpv[idx] = (
                        (tp[idx] + (org2tpert * org[idx] + tpert[i]))
                        * (1.0 + 1.608 * qstp[idx])
                        / (1.0 + qstp[idx])
                    )
                else:
                    tpv[idx] = (
                        (tp[idx] + tpert[i])
                        * (1.0 + 1.608 * qstp[idx])
                        / (1.0 + qstp[idx])
                    )

            if kk < launch_i:
                idxp1 = i + kk * pcols
                new_q = 0.0
                new_s = 0.0
                ii = 0
                while ii <= nit_lheat - 1:
                    xsh2o[idx] = max(0.0, qtmix[idx] - qsmix[idx] - lwmax)
                    ds_xsh2o[idx] = (
                        ds_xsh2o[idxp1]
                        - cpliq * log(tmix[idx] / tfreez) * max(0.0, (xsh2o[idx] - xsh2o[idxp1]))
                    )

                    if tmix[idx] <= tfreez + tscool and ds_freeze[idxp1] == 0.0:
                        ds_freeze[idx] = (
                            (latice / tmix[idx])
                            * max(0.0, qtmix[idx] - qsmix[idx] - xsh2o[idx])
                        )

                    if tmix[idx] <= tfreez + tscool and ds_freeze[idxp1] != 0.0:
                        ds_freeze[idx] = (
                            ds_freeze[idxp1]
                            + (latice / tmix[idx]) * max(0.0, (qsmix[idxp1] - qsmix[idx]))
                        )

                    new_s = smix[idx] + ds_xsh2o[idx] + ds_freeze[idx]
                    new_q = qtmix[idx] - xsh2o[idx]

                    tmix_val = 0.0
                    qsmix_val = 0.0
                    conv = 0
                    _zm_ientropy_codon(
                        new_s,
                        p[idx],
                        new_q,
                        tmix[idx],
                        rl,
                        cpliq,
                        cpwv,
                        tfreez,
                        cpres,
                        rgas,
                        eps1,
                        rh2o,
                        wv_idx,
                        epsilo,
                        omeps,
                        __ptr__(tmix_val),
                        __ptr__(qsmix_val),
                        __ptr__(conv),
                    )
                    if conv == 0:
                        status[0] = 0
                        return
                    tmix[idx] = tmix_val
                    qsmix[idx] = qsmix_val
                    ii += 1

                tp[idx] = tmix[idx]

                if new_q > qsmix[idx]:
                    qstp[idx] = qsmix[idx]
                else:
                    qstp[idx] = new_q

                if zm_org != 0:
                    tpv[idx] = (
                        (tp[idx] + (org2tpert * org[idx] + tpert[i]))
                        * (1.0 + 1.608 * qstp[idx])
                        / (1.0 + new_q)
                    )
                else:
                    tpv[idx] = (
                        (tp[idx] + tpert[i])
                        * (1.0 + 1.608 * qstp[idx])
                        / (1.0 + new_q)
                    )

            i += 1
        kk -= 1


@inline
def _zm_nint_codon(value: float) -> int:
    if value >= 0.0:
        return int(value + 0.5)
    return int(value - 0.5)


@export
def zm_buoyan_dilute_codon(
    lchnk: int,
    ncol: int,
    msg: int,
    pcols: int,
    pver: int,
    pverp: int,
    zm_org: int,
    tiedke_add: float,
    rl: float,
    rd: float,
    grav: float,
    cp: float,
    rgas: float,
    cpliq: float,
    tfreez: float,
    latice: float,
    cpwv: float,
    cpres: float,
    eps1: float,
    rh2o: float,
    epsilo: float,
    omeps: float,
    wv_idx: int,
    q_p: cobj,
    t_p: cobj,
    p_p: cobj,
    z_p: cobj,
    pf_p: cobj,
    pblt_p: cobj,
    tpert_p: cobj,
    tp_p: cobj,
    qstp_p: cobj,
    tl_p: cobj,
    cape_p: cobj,
    lcl_p: cobj,
    lel_p: cobj,
    lon_p: cobj,
    mx_p: cobj,
    org_p: cobj,
    landfrac_p: cobj,
    capeten_p: cobj,
    tv_p: cobj,
    tpv_p: cobj,
    buoy_p: cobj,
    pl_p: cobj,
    hmax_p: cobj,
    hmn_p: cobj,
    knt_p: cobj,
    lelten_p: cobj,
    tmix_p: cobj,
    qtmix_p: cobj,
    qsmix_p: cobj,
    smix_p: cobj,
    xsh2o_p: cobj,
    ds_xsh2o_p: cobj,
    ds_freeze_p: cobj,
    mp_p: cobj,
    qtp_p: cobj,
    sp_p: cobj,
    sp0_p: cobj,
    qtp0_p: cobj,
    mp0_p: cobj,
    status_p: cobj,
):
    q = Ptr[float](q_p)
    t = Ptr[float](t_p)
    p = Ptr[float](p_p)
    z = Ptr[float](z_p)
    pf = Ptr[float](pf_p)
    pblt = Ptr[float](pblt_p)
    tpert = Ptr[float](tpert_p)
    tp = Ptr[float](tp_p)
    qstp = Ptr[float](qstp_p)
    tl = Ptr[float](tl_p)
    cape = Ptr[float](cape_p)
    lcl = Ptr[i32](lcl_p)
    lel = Ptr[i32](lel_p)
    lon = Ptr[i32](lon_p)
    mx = Ptr[i32](mx_p)
    org = Ptr[float](org_p)
    landfrac = Ptr[float](landfrac_p)
    capeten = Ptr[float](capeten_p)
    tv = Ptr[float](tv_p)
    tpv = Ptr[float](tpv_p)
    buoy = Ptr[float](buoy_p)
    pl = Ptr[float](pl_p)
    hmax = Ptr[float](hmax_p)
    hmn = Ptr[float](hmn_p)
    knt = Ptr[i32](knt_p)
    lelten = Ptr[i32](lelten_p)
    status = Ptr[int](status_p)

    status[0] = 1

    n0 = 0
    while n0 < 5:
        i = 0
        while i < ncol:
            lelten[i + n0 * pcols] = i32(pver)
            capeten[i + n0 * pcols] = 0.0
            i += 1
        n0 += 1

    i = 0
    while i < ncol:
        lon[i] = i32(pver)
        knt[i] = i32(0)
        lel[i] = i32(pver)
        mx[i] = lon[i]
        cape[i] = 0.0
        hmax[i] = 0.0
        i += 1

    k0 = 0
    while k0 < pver:
        i = 0
        while i < ncol:
            idx = i + k0 * pcols
            tp[idx] = t[idx]
            qstp[idx] = q[idx]
            tv[idx] = t[idx] * (1.0 + 1.608 * q[idx]) / (1.0 + q[idx])
            tpv[idx] = tv[idx]
            buoy[idx] = 0.0
            i += 1
        k0 += 1

    kk = pver
    while kk >= msg + 1:
        k0 = kk - 1
        i = 0
        while i < ncol:
            idx = i + k0 * pcols
            hmn[i] = cp * t[idx] + grav * z[idx] + rl * q[idx]
            if kk >= _zm_nint_codon(pblt[i]) and kk <= int(lon[i]) and hmn[i] > hmax[i]:
                hmax[i] = hmn[i]
                mx[i] = i32(kk)
            i += 1
        kk -= 1

    i = 0
    while i < ncol:
        mx_i = int(mx[i])
        mx_idx = i + (mx_i - 1) * pcols
        lcl[i] = mx[i]
        tl[i] = t[mx_idx]
        pl[i] = p[mx_idx]
        i += 1

    zm_parcel_dilute_codon(
        lchnk,
        ncol,
        msg,
        pcols,
        pver,
        zm_org,
        grav,
        rgas,
        cpliq,
        tfreez,
        latice,
        rl,
        cpwv,
        cpres,
        eps1,
        rh2o,
        epsilo,
        omeps,
        wv_idx,
        mx_p,
        p_p,
        t_p,
        q_p,
        tpert_p,
        tp_p,
        tpv_p,
        qstp_p,
        pl_p,
        tl_p,
        lcl_p,
        org_p,
        landfrac_p,
        tmix_p,
        qtmix_p,
        qsmix_p,
        smix_p,
        xsh2o_p,
        ds_xsh2o_p,
        ds_freeze_p,
        mp_p,
        qtp_p,
        sp_p,
        sp0_p,
        qtp0_p,
        mp0_p,
        status_p,
    )
    if status[0] == 0:
        return

    kk = pver
    while kk >= msg + 1:
        k0 = kk - 1
        i = 0
        while i < ncol:
            idx = i + k0 * pcols
            if kk <= int(mx[i]) and pl[i] >= 600.0:
                tv[idx] = t[idx] * (1.0 + 1.608 * q[idx]) / (1.0 + q[idx])
                buoy[idx] = tpv[idx] - tv[idx] + tiedke_add
            else:
                qstp[idx] = q[idx]
                tp[idx] = t[idx]
                tpv[idx] = tv[idx]
            i += 1
        kk -= 1

    kk = msg + 2
    while kk <= pver:
        k0 = kk - 1
        i = 0
        while i < ncol:
            if kk < int(lcl[i]) and pl[i] >= 600.0:
                idx = i + k0 * pcols
                idxp1 = i + kk * pcols
                if buoy[idxp1] > 0.0 and buoy[idx] <= 0.0:
                    knt_i = min(5, int(knt[i]) + 1)
                    knt[i] = i32(knt_i)
                    lelten[i + (knt_i - 1) * pcols] = i32(kk)
            i += 1
        kk += 1

    n0 = 0
    while n0 < 5:
        kk = msg + 1
        while kk <= pver:
            k0 = kk - 1
            i = 0
            while i < ncol:
                if pl[i] >= 600.0 and kk <= int(mx[i]) and kk > int(lelten[i + n0 * pcols]):
                    idx = i + k0 * pcols
                    capeten[i + n0 * pcols] = (
                        capeten[i + n0 * pcols] + rd * buoy[idx] * log(pf[i + kk * pcols] / pf[idx])
                    )
                i += 1
            kk += 1
        n0 += 1

    n0 = 0
    while n0 < 5:
        i = 0
        while i < ncol:
            capeten_idx = i + n0 * pcols
            if capeten[capeten_idx] > cape[i]:
                cape[i] = capeten[capeten_idx]
                lel[i] = lelten[capeten_idx]
            i += 1
        n0 += 1

    i = 0
    while i < ncol:
        cape[i] = max(cape[i], 0.0)
        i += 1


@export
def zm_conv_tend_2_lq_mask_codon(
    pcnst: int,
    wtrc_nwset: int,
    iwtvap: int,
    iwtcvsnow: int,
    trace_water: int,
    convtran1_mask_p: cobj,
    wtrc_iatype_p: cobj,
    lq_mask_p: cobj,
):
    convtran1_mask = Ptr[int](convtran1_mask_p)
    wtrc_iatype = Ptr[int](wtrc_iatype_p)
    lq_mask = Ptr[int](lq_mask_p)

    i = 0
    while i < pcnst:
        if convtran1_mask[i] != 0:
            lq_mask[i] = 0
        else:
            lq_mask[i] = 1
        i += 1

    if trace_water != 0:
        m = iwtvap
        while m <= iwtcvsnow:
            moff = (m - iwtvap) * wtrc_nwset
            i = 0
            while i < wtrc_nwset:
                idx = wtrc_iatype[i + moff] - 1
                lq_mask[idx] = 0
                i += 1
            m += 1


@export
def zm_conv_post_stage_dispatch_codon(
    stage: int,
    mode: int,
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    pcnst: int,
    wtrc_nwset: int,
    lengath: int,
    ixorg: int,
    do_org: int,
    vap_type: int,
    gravit: float,
    cpair: float,
    ztodt: float,
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
):
    if stage == 1:
        zm_conv_workspace_init_shell_codon(ncol, pcols, pver, p1, p2, p3, p4)
    elif stage == 2:
        zm_conv_ptend_lq_mask_shell_codon(pcnst, wtrc_nwset, do_org, ixorg, p1, p2)
    elif stage == 3:
        zm_wtrc_convr_prep_shell_codon(ncol, pcols, pver, pcnst, wtrc_nwset, p1, p2, p3, p4, p5)
    elif stage == 4:
        zm_wtrc_precip_assign_shell_codon(pcols, vap_type, p1, p2, p3, p4)
    elif stage == 5:
        zm_convr_post_shell_codon(
            ncol,
            pcols,
            pver,
            pverp,
            lengath,
            gravit,
            cpair,
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
        )
    elif stage == 6:
        zm_conv_evap_prep_shell_codon(ncol, pcols, pver, p1, p2)
    elif stage == 7:
        zm_conv_evap_post_shell_codon(ncol, pcols, pver, ixorg, do_org, ztodt, p1, p2, p3)
    elif stage == 8:
        zm_conv_evap_hist_shell_codon(mode, ncol, pcols, pver, cpair, p1, p2, p3, p4)
    elif stage == 9:
        zm_momtran_prep_shell_codon(ncol, pcols, pver, p1, p2, p3)
    elif stage == 10:
        zm_momtran_post_shell_codon(ncol, pcols, pver, cpair, p1, p2, p3, p4, p5, p6)
    elif stage == 11:
        zm_convtran1_prep_shell_codon(pcols, pver, p1)
    elif stage == 12:
        zm_convtran1_ratio_shell_codon(ncol, pcols, pver, wtrc_nwset, p1, p2, p3, p4)
    elif stage == 13:
        zm_convtran2_dpdry_shell_codon(pcols, pver, lengath, p1, p2, p3)
