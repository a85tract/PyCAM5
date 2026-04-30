from math import atan, copysign, erfc, exp, log, pi, sqrt
from C import modal_aero_kohler_native_cb(float, float, float) -> float
from C import modal_aero_kohler_cubic_real_root_native_cb(float, float, float, float, float) -> float
from C import modal_aero_kohler_quartic_real_root_native_cb(float, float, float, float, float, float) -> float
from C import modal_aero_complex_sqrt_native_cb(float, float, Ptr[float], Ptr[float]) -> None
from C import modal_aero_complex_pow_third_native_cb(float, float, Ptr[float], Ptr[float]) -> None
from C import modal_aero_vol_from_radius_native_cb(float) -> float

@inline
def _idx2(i: int, k: int, ld1: int) -> int:
    """Fortran array declared as (ld1, *)."""
    return (i - 1) + (k - 1) * ld1


@inline
def _idx3(i: int, k: int, m: int, ld1: int, ld2: int) -> int:
    """Fortran array declared as (ld1, ld2, *)."""
    return (i - 1) + (k - 1) * ld1 + (m - 1) * ld1 * ld2


@inline
def _idx4(i: int, j: int, k: int, l: int, ld1: int, ld2: int, ld3: int) -> int:
    """Fortran array declared as (ld1, ld2, ld3, *)."""
    return (i - 1) + (j - 1) * ld1 + (k - 1) * ld1 * ld2 + (l - 1) * ld1 * ld2 * ld3


@inline
def _modal_aero_v2ncur(dgncur_a: float, pi_const: float, alnsg: float) -> float:
    return 1.0 / ((pi_const / 6.0) * (dgncur_a**3.0) * exp(4.5 * (alnsg**2.0)))


@inline
def _modal_aero_radius_from_vol(vol: float, pi43_const: float) -> float:
    return (vol / pi43_const) ** (1.0 / 3.0)


def _complex_sqrt(z):
    return z ** 0.5


def _complex_sqrt_native(z):
    root_re = 0.0
    root_im = 0.0
    modal_aero_complex_sqrt_native_cb(
        z.real, z.imag, __ptr__(root_re), __ptr__(root_im)
    )
    return complex(root_re, root_im)


def _complex_pow_third_native(z):
    root_re = 0.0
    root_im = 0.0
    modal_aero_complex_pow_third_native_cb(
        z.real, z.imag, __ptr__(root_re), __ptr__(root_im)
    )
    return complex(root_re, root_im)


def _modal_aero_kohler_cubic_real_root(p2: float, p1: float, p0: float, rdry: float, eps: float) -> float:
    third = 1.0 / 3.0
    ci = complex(0.0, 1.0)
    sqrt3 = sqrt(3.0)
    cw = 0.5 * (-1.0 + ci * sqrt3)
    cwsq = 0.5 * (-1.0 - ci * sqrt3)

    if p1 == 0.0:
        root = (-p0) ** third
        return root

    q = p1 / 3.0
    r = p0 / 2.0
    crad = complex(r * r + q * q * q, 0.0)
    crad = _complex_sqrt(crad)

    cy = complex(r, 0.0) - crad
    if abs(cy) > eps:
        cy = cy ** third
    cq = complex(q, 0.0)
    cz = -cq / cy

    cx1 = -cy - cz
    cx2 = -cw * cy - cwsq * cz
    cx3 = -cwsq * cy - cw * cz

    root = 1000.0 * rdry
    nsol = 0

    xr = cx1.real
    xi = cx1.imag
    if abs(xi) <= abs(xr) * eps:
        if xr <= root:
            if xr >= rdry * (1.0 - eps):
                if xr == xr:
                    root = xr
                    nsol = 1

    xr = cx2.real
    xi = cx2.imag
    if abs(xi) <= abs(xr) * eps:
        if xr <= root:
            if xr >= rdry * (1.0 - eps):
                if xr == xr:
                    root = xr
                    nsol = 2

    xr = cx3.real
    xi = cx3.imag
    if abs(xi) <= abs(xr) * eps:
        if xr <= root:
            if xr >= rdry * (1.0 - eps):
                if xr == xr:
                    root = xr
                    nsol = 3

    if nsol == 0:
        root = rdry

    return root


def _modal_aero_kohler_quartic_real_root(
    p3: float, p2: float, p1: float, p0: float, rdry: float, eps: float
) -> float:
    third = 1.0 / 3.0
    czero = complex(0.0, 0.0)

    q = -(p2 * p2) / 36.0 + (p3 * p1 - 4 * p0) / 12.0
    r = -((p2 / 6) ** 3) + p2 * (p3 * p1 - 4 * p0) / 48.0 + (
        4 * p0 * p2 - p0 * p3 * p3 - p1 * p1
    ) / 16.0

    crad = complex(r * r + q * q * q, 0.0)
    crad = _complex_sqrt(crad)

    cb = complex(r, 0.0) - crad
    if cb == czero:
        cx1 = complex((-p1) ** third, 0.0)
        cx2 = cx1
        cx3 = cx1
        cx4 = cx1
    else:
        cb = cb ** third

        cy = -cb + q / cb + p2 / 6

        cb0 = _complex_sqrt(cy * cy - p0)
        cb1 = (p3 * cy - p1) / (2 * cb0)

        cb = p3 / 2 + cb1
        crad = cb * cb - 4 * (cy + cb0)
        crad = _complex_sqrt(crad)
        cx1 = (-cb + crad) / 2.0
        cx2 = (-cb - crad) / 2.0

        cb = p3 / 2 - cb1
        crad = cb * cb - 4 * (cy - cb0)
        crad = _complex_sqrt(crad)
        cx3 = (-cb + crad) / 2.0
        cx4 = (-cb - crad) / 2.0

    root = 1000.0 * rdry
    nsol = 0

    xr = cx1.real
    xi = cx1.imag
    if abs(xi) <= abs(xr) * eps:
        if xr <= root:
            if xr >= rdry * (1.0 - eps):
                if xr == xr:
                    root = xr
                    nsol = 1

    xr = cx2.real
    xi = cx2.imag
    if abs(xi) <= abs(xr) * eps:
        if xr <= root:
            if xr >= rdry * (1.0 - eps):
                if xr == xr:
                    root = xr
                    nsol = 2

    xr = cx3.real
    xi = cx3.imag
    if abs(xi) <= abs(xr) * eps:
        if xr <= root:
            if xr >= rdry * (1.0 - eps):
                if xr == xr:
                    root = xr
                    nsol = 3

    xr = cx4.real
    xi = cx4.imag
    if abs(xi) <= abs(xr) * eps:
        if xr <= root:
            if xr >= rdry * (1.0 - eps):
                if xr == xr:
                    root = xr
                    nsol = 4

    if nsol == 0:
        root = rdry

    return root


def _modal_aero_kohler_quartic_real_root_sqrt_native(
    p3: float, p2: float, p1: float, p0: float, rdry: float, eps: float
) -> float:
    third = 1.0 / 3.0
    czero = complex(0.0, 0.0)

    q = -(p2 * p2) / 36.0 + (p3 * p1 - 4 * p0) / 12.0
    r = -((p2 / 6) ** 3) + p2 * (p3 * p1 - 4 * p0) / 48.0 + (
        4 * p0 * p2 - p0 * p3 * p3 - p1 * p1
    ) / 16.0

    crad = complex(r * r + q * q * q, 0.0)
    crad = _complex_sqrt_native(crad)

    cb = complex(r, 0.0) - crad
    if cb == czero:
        cx1 = complex((-p1) ** third, 0.0)
        cx2 = cx1
        cx3 = cx1
        cx4 = cx1
    else:
        cb = cb ** third

        cy = -cb + q / cb + p2 / 6

        cb0 = _complex_sqrt_native(cy * cy - p0)
        cb1 = (p3 * cy - p1) / (2 * cb0)

        cb = p3 / 2 + cb1
        crad = cb * cb - 4 * (cy + cb0)
        crad = _complex_sqrt_native(crad)
        cx1 = (-cb + crad) / 2.0
        cx2 = (-cb - crad) / 2.0

        cb = p3 / 2 - cb1
        crad = cb * cb - 4 * (cy - cb0)
        crad = _complex_sqrt_native(crad)
        cx3 = (-cb + crad) / 2.0
        cx4 = (-cb - crad) / 2.0

    root = 1000.0 * rdry
    nsol = 0

    xr = cx1.real
    xi = cx1.imag
    if abs(xi) <= abs(xr) * eps:
        if xr <= root:
            if xr >= rdry * (1.0 - eps):
                if xr == xr:
                    root = xr
                    nsol = 1

    xr = cx2.real
    xi = cx2.imag
    if abs(xi) <= abs(xr) * eps:
        if xr <= root:
            if xr >= rdry * (1.0 - eps):
                if xr == xr:
                    root = xr
                    nsol = 2

    xr = cx3.real
    xi = cx3.imag
    if abs(xi) <= abs(xr) * eps:
        if xr <= root:
            if xr >= rdry * (1.0 - eps):
                if xr == xr:
                    root = xr
                    nsol = 3

    xr = cx4.real
    xi = cx4.imag
    if abs(xi) <= abs(xr) * eps:
        if xr <= root:
            if xr >= rdry * (1.0 - eps):
                if xr == xr:
                    root = xr
                    nsol = 4

    if nsol == 0:
        root = rdry

    return root


def _modal_aero_kohler_quartic_real_root_pow_native(
    p3: float, p2: float, p1: float, p0: float, rdry: float, eps: float
) -> float:
    third = 1.0 / 3.0
    czero = complex(0.0, 0.0)

    q = -(p2 * p2) / 36.0 + (p3 * p1 - 4 * p0) / 12.0
    r = -((p2 / 6) ** 3) + p2 * (p3 * p1 - 4 * p0) / 48.0 + (
        4 * p0 * p2 - p0 * p3 * p3 - p1 * p1
    ) / 16.0

    crad = complex(r * r + q * q * q, 0.0)
    crad = _complex_sqrt(crad)

    cb = complex(r, 0.0) - crad
    if cb == czero:
        cx1 = complex((-p1) ** third, 0.0)
        cx2 = cx1
        cx3 = cx1
        cx4 = cx1
    else:
        cb = _complex_pow_third_native(cb)

        cy = -cb + q / cb + p2 / 6

        cb0 = _complex_sqrt(cy * cy - p0)
        cb1 = (p3 * cy - p1) / (2 * cb0)

        cb = p3 / 2 + cb1
        crad = cb * cb - 4 * (cy + cb0)
        crad = _complex_sqrt(crad)
        cx1 = (-cb + crad) / 2.0
        cx2 = (-cb - crad) / 2.0

        cb = p3 / 2 - cb1
        crad = cb * cb - 4 * (cy - cb0)
        crad = _complex_sqrt(crad)
        cx3 = (-cb + crad) / 2.0
        cx4 = (-cb - crad) / 2.0

    root = 1000.0 * rdry
    nsol = 0

    xr = cx1.real
    xi = cx1.imag
    if abs(xi) <= abs(xr) * eps:
        if xr <= root:
            if xr >= rdry * (1.0 - eps):
                if xr == xr:
                    root = xr
                    nsol = 1

    xr = cx2.real
    xi = cx2.imag
    if abs(xi) <= abs(xr) * eps:
        if xr <= root:
            if xr >= rdry * (1.0 - eps):
                if xr == xr:
                    root = xr
                    nsol = 2

    xr = cx3.real
    xi = cx3.imag
    if abs(xi) <= abs(xr) * eps:
        if xr <= root:
            if xr >= rdry * (1.0 - eps):
                if xr == xr:
                    root = xr
                    nsol = 3

    xr = cx4.real
    xi = cx4.imag
    if abs(xi) <= abs(xr) * eps:
        if xr <= root:
            if xr >= rdry * (1.0 - eps):
                if xr == xr:
                    root = xr
                    nsol = 4

    if nsol == 0:
        root = rdry

    return root


def _modal_aero_kohler_quartic_real_root_sqrt_pow_native(
    p3: float, p2: float, p1: float, p0: float, rdry: float, eps: float
) -> float:
    third = 1.0 / 3.0
    czero = complex(0.0, 0.0)

    q = -(p2 * p2) / 36.0 + (p3 * p1 - 4 * p0) / 12.0
    r = -((p2 / 6) ** 3) + p2 * (p3 * p1 - 4 * p0) / 48.0 + (
        4 * p0 * p2 - p0 * p3 * p3 - p1 * p1
    ) / 16.0

    crad = complex(r * r + q * q * q, 0.0)
    crad = _complex_sqrt_native(crad)

    cb = complex(r, 0.0) - crad
    if cb == czero:
        cx1 = complex((-p1) ** third, 0.0)
        cx2 = cx1
        cx3 = cx1
        cx4 = cx1
    else:
        cb = _complex_pow_third_native(cb)

        cy = -cb + q / cb + p2 / 6

        cb0 = _complex_sqrt_native(cy * cy - p0)
        cb1 = (p3 * cy - p1) / (2 * cb0)

        cb = p3 / 2 + cb1
        crad = cb * cb - 4 * (cy + cb0)
        crad = _complex_sqrt_native(crad)
        cx1 = (-cb + crad) / 2.0
        cx2 = (-cb - crad) / 2.0

        cb = p3 / 2 - cb1
        crad = cb * cb - 4 * (cy - cb0)
        crad = _complex_sqrt_native(crad)
        cx3 = (-cb + crad) / 2.0
        cx4 = (-cb - crad) / 2.0

    root = 1000.0 * rdry
    nsol = 0

    xr = cx1.real
    xi = cx1.imag
    if abs(xi) <= abs(xr) * eps:
        if xr <= root:
            if xr >= rdry * (1.0 - eps):
                if xr == xr:
                    root = xr
                    nsol = 1

    xr = cx2.real
    xi = cx2.imag
    if abs(xi) <= abs(xr) * eps:
        if xr <= root:
            if xr >= rdry * (1.0 - eps):
                if xr == xr:
                    root = xr
                    nsol = 2

    xr = cx3.real
    xi = cx3.imag
    if abs(xi) <= abs(xr) * eps:
        if xr <= root:
            if xr >= rdry * (1.0 - eps):
                if xr == xr:
                    root = xr
                    nsol = 3

    xr = cx4.real
    xi = cx4.imag
    if abs(xi) <= abs(xr) * eps:
        if xr <= root:
            if xr >= rdry * (1.0 - eps):
                if xr == xr:
                    root = xr
                    nsol = 4

    if nsol == 0:
        root = rdry

    return root


def _modal_aero_kohler_scalar_selective(
    dryrad_in: float, hygro: float, s: float, quartic_mode: int, cubic_mode: int
) -> float:
    eps = 1.0e-4
    mw = 18.0
    rhow = 1.0
    surften = 76.0
    tair = 273.0
    third = 1.0 / 3.0
    ugascon = 8.3e7

    a = 2.0e4 * mw * surften / (ugascon * tair * rhow)

    rdry = dryrad_in * 1.0e6
    vol = rdry**3
    b = vol * hygro

    ss = min(s, 1.0 - eps)
    ss = max(ss, 1.0e-10)
    slog = log(ss)
    p43 = -a / slog
    p42 = 0.0
    p41 = b / slog - vol
    p40 = a * vol / slog
    p32 = 0.0
    p31 = -b / a
    p30 = -vol
    r = rdry
    r3 = rdry
    r4 = rdry

    if vol <= 1.0e-12:
        r = rdry
    else:
        p = abs(p31) / (rdry * rdry)
        if p < eps:
            r = rdry * (1.0 + p * third / (1.0 - slog * rdry / a))
        else:
            if quartic_mode == 1:
                r = modal_aero_kohler_quartic_real_root_native_cb(
                    p43, p42, p41, p40, rdry, eps
                )
            elif quartic_mode == 2:
                r = _modal_aero_kohler_quartic_real_root_sqrt_native(
                    p43, p42, p41, p40, rdry, eps
                )
            elif quartic_mode == 3:
                r = _modal_aero_kohler_quartic_real_root_pow_native(
                    p43, p42, p41, p40, rdry, eps
                )
            elif quartic_mode == 4:
                r = _modal_aero_kohler_quartic_real_root_sqrt_pow_native(
                    p43, p42, p41, p40, rdry, eps
                )
            else:
                r = _modal_aero_kohler_quartic_real_root(
                    p43, p42, p41, p40, rdry, eps
                )

    if s > 1.0 - eps:
        r4 = r
        p = abs(p31) / (rdry * rdry)
        if p < eps:
            r = rdry * (1.0 + p * third)
        else:
            if cubic_mode == 1:
                r = modal_aero_kohler_cubic_real_root_native_cb(
                    p32, p31, p30, rdry, eps
                )
            else:
                r = _modal_aero_kohler_cubic_real_root(p32, p31, p30, rdry, eps)
        r3 = r
        r = (r4 * (1.0 - s) + r3 * (s - 1.0 + eps)) / eps

    r = min(r, 30.0)
    return r * 1.0e-6


def _modal_aero_kohler_scalar_all_codon(dryrad_in: float, hygro: float, s: float) -> float:
    return _modal_aero_kohler_scalar_selective(dryrad_in, hygro, s, 0, 0)


def _modal_aero_kohler_scalar_native_roots(dryrad_in: float, hygro: float, s: float) -> float:
    return _modal_aero_kohler_scalar_selective(dryrad_in, hygro, s, 1, 1)


def _modal_aero_kohler_scalar_quartic_native(dryrad_in: float, hygro: float, s: float) -> float:
    return _modal_aero_kohler_scalar_selective(dryrad_in, hygro, s, 1, 0)


def _modal_aero_kohler_scalar_cubic_native(dryrad_in: float, hygro: float, s: float) -> float:
    return _modal_aero_kohler_scalar_selective(dryrad_in, hygro, s, 0, 1)


def _modal_aero_kohler_scalar_quartic_sqrt_native(
    dryrad_in: float, hygro: float, s: float
) -> float:
    return _modal_aero_kohler_scalar_selective(dryrad_in, hygro, s, 2, 0)


def _modal_aero_kohler_scalar_quartic_pow_native(
    dryrad_in: float, hygro: float, s: float
) -> float:
    return _modal_aero_kohler_scalar_selective(dryrad_in, hygro, s, 3, 0)


def _modal_aero_kohler_scalar_quartic_sqrt_pow_native(
    dryrad_in: float, hygro: float, s: float
) -> float:
    return _modal_aero_kohler_scalar_selective(dryrad_in, hygro, s, 4, 0)


def _modal_aero_kohler_scalar_sat_native(dryrad_in: float, hygro: float, s: float) -> float:
    eps = 1.0e-4

    if s > 1.0 - eps:
        return modal_aero_kohler_native_cb(dryrad_in, hygro, s)

    return _modal_aero_kohler_scalar_native_roots(dryrad_in, hygro, s)


def _modal_aero_kohler_scalar_subsat_native(dryrad_in: float, hygro: float, s: float) -> float:
    eps = 1.0e-4

    if s > 1.0 - eps:
        return _modal_aero_kohler_scalar_native_roots(dryrad_in, hygro, s)

    return modal_aero_kohler_native_cb(dryrad_in, hygro, s)


@inline
def _max3(a: float, b: float, c: float) -> float:
    if a >= b:
        if a >= c:
            return a
        return c
    if b >= c:
        return b
    return c


@inline
def _min4(a: float, b: float, c: float, d: float) -> float:
    return min(min(a, b), min(c, d))


@inline
def _max4(a: float, b: float, c: float, d: float) -> float:
    return max(max(a, b), max(c, d))


@inline
def _sign_one(x: float) -> float:
    return copysign(1.0, x)


@inline
def _minmod(a: float, b: float) -> float:
    return 0.5 * (_sign_one(a) + _sign_one(b)) * min(abs(a), abs(b))


@inline
def _medan(a: float, b: float, c: float) -> float:
    return a + _minmod(b - a, c - a)


@export
def modal_aero_depvel_part_codon(
    ncol: int,
    pcols: int,
    pver: int,
    n_land_type: int,
    moment: int,
    pi: float,
    boltz: float,
    gravit: float,
    rair: float,
    t_p: cobj,
    pmid_p: cobj,
    ram1_p: cobj,
    fv_p: cobj,
    vlc_dry_p: cobj,
    vlc_trb_p: cobj,
    vlc_grv_p: cobj,
    radius_part_p: cobj,
    density_part_p: cobj,
    sig_part_p: cobj,
    fraction_landuse_p: cobj,
    vsc_dyn_atm_p: cobj,
    vsc_knm_atm_p: cobj,
    mfp_atm_p: cobj,
    slp_crc_p: cobj,
    radius_moment_p: cobj,
):
    t = Ptr[float](t_p)
    pmid = Ptr[float](pmid_p)
    ram1 = Ptr[float](ram1_p)
    fv = Ptr[float](fv_p)
    vlc_dry = Ptr[float](vlc_dry_p)
    vlc_trb = Ptr[float](vlc_trb_p)
    vlc_grv = Ptr[float](vlc_grv_p)
    radius_part = Ptr[float](radius_part_p)
    density_part = Ptr[float](density_part_p)
    sig_part = Ptr[float](sig_part_p)
    fraction_landuse = Ptr[float](fraction_landuse_p)
    vsc_dyn_atm = Ptr[float](vsc_dyn_atm_p)
    vsc_knm_atm = Ptr[float](vsc_knm_atm_p)
    mfp_atm = Ptr[float](mfp_atm_p)
    slp_crc = Ptr[float](slp_crc_p)
    radius_moment = Ptr[float](radius_moment_p)

    gamma = (0.56, 0.54, 0.54, 0.56, 0.56, 0.56, 0.50, 0.54, 0.54, 0.54, 0.54)
    alpha = (1.50, 1.20, 1.20, 0.80, 1.00, 0.80, 100.00, 50.00, 2.00, 1.20, 50.00)
    radius_collector = (
        10.00e-03,
        3.50e-03,
        3.50e-03,
        5.10e-03,
        2.00e-03,
        5.00e-03,
        -1.00,
        -1.00,
        10.00e-03,
        3.50e-03,
        -1.00,
    )
    iwet = (-1, -1, -1, -1, -1, -1, 1, -1, 1, -1, -1)

    moment_f = float(moment)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            lnsig = log(sig_part[_idx2(i, k, pcols)])
            radius_moment[_idx2(i, k, pcols)] = min(
                50.0e-6, radius_part[_idx2(i, k, pcols)]
            ) * exp((moment_f - 1.5) * lnsig * lnsig)
            dispersion = exp(2.0 * lnsig * lnsig)

            rho = pmid[_idx2(i, k, pcols)] / rair / t[_idx2(i, k, pcols)]

            vsc_dyn_atm[_idx2(i, k, pcols)] = (
                1.72e-5
                * ((t[_idx2(i, k, pcols)] / 273.0) ** 1.5)
                * 393.0
                / (t[_idx2(i, k, pcols)] + 120.0)
            )
            mfp_atm[_idx2(i, k, pcols)] = (
                2.0
                * vsc_dyn_atm[_idx2(i, k, pcols)]
                / (
                    pmid[_idx2(i, k, pcols)]
                    * sqrt(8.0 / (pi * rair * t[_idx2(i, k, pcols)]))
                )
            )
            vsc_knm_atm[_idx2(i, k, pcols)] = vsc_dyn_atm[_idx2(i, k, pcols)] / rho

            slp_crc[_idx2(i, k, pcols)] = (
                1.0
                + mfp_atm[_idx2(i, k, pcols)]
                * (
                    1.257
                    + 0.4
                    * exp(
                        -1.1
                        * radius_moment[_idx2(i, k, pcols)]
                        / mfp_atm[_idx2(i, k, pcols)]
                    )
                )
                / radius_moment[_idx2(i, k, pcols)]
            )
            vlc_grv[_idx2(i, k, pcols)] = (
                (4.0 / 18.0)
                * radius_moment[_idx2(i, k, pcols)]
                * radius_moment[_idx2(i, k, pcols)]
                * density_part[_idx2(i, k, pcols)]
                * gravit
                * slp_crc[_idx2(i, k, pcols)]
                / vsc_dyn_atm[_idx2(i, k, pcols)]
            )
            vlc_grv[_idx2(i, k, pcols)] = vlc_grv[_idx2(i, k, pcols)] * dispersion

            vlc_dry[_idx2(i, k, pcols)] = vlc_grv[_idx2(i, k, pcols)]

    k = pver
    for i in range(1, ncol + 1):
        dff_aer = (
            boltz
            * t[_idx2(i, k, pcols)]
            * slp_crc[_idx2(i, k, pcols)]
            / (6.0 * pi * vsc_dyn_atm[_idx2(i, k, pcols)] * radius_moment[_idx2(i, k, pcols)])
        )
        shm_nbr = vsc_knm_atm[_idx2(i, k, pcols)] / dff_aer

        wrk2 = 0.0
        wrk3 = 0.0
        for lt in range(1, n_land_type + 1):
            lnd_frc = fraction_landuse[_idx2(i, lt, pcols)]
            if lnd_frc != 0.0:
                brownian = shm_nbr ** (-gamma[lt - 1])
                if radius_collector[lt - 1] > 0.0:
                    stk_nbr = (
                        vlc_grv[_idx2(i, k, pcols)] * fv[i - 1] / (gravit * radius_collector[lt - 1])
                    )
                    interception = 2.0 * (
                        radius_moment[_idx2(i, k, pcols)] / radius_collector[lt - 1]
                    ) ** 2.0
                else:
                    stk_nbr = (
                        vlc_grv[_idx2(i, k, pcols)]
                        * fv[i - 1]
                        * fv[i - 1]
                        / (gravit * vsc_knm_atm[_idx2(i, k, pcols)])
                    )
                    interception = 0.0
                impaction = (stk_nbr / (alpha[lt - 1] + stk_nbr)) ** 2.0

                if iwet[lt - 1] > 0:
                    stickfrac = 1.0
                else:
                    stickfrac = exp(-sqrt(stk_nbr))
                    if stickfrac < 1.0e-10:
                        stickfrac = 1.0e-10
                rss_lmn = 1.0 / (
                    3.0 * fv[i - 1] * stickfrac * (brownian + interception + impaction)
                )
                rss_trb = (
                    ram1[i - 1]
                    + rss_lmn
                    + ram1[i - 1] * rss_lmn * vlc_grv[_idx2(i, k, pcols)]
                )

                wrk1 = 1.0 / rss_trb
                wrk2 = wrk2 + lnd_frc * wrk1
                wrk3 = wrk3 + lnd_frc * (wrk1 + vlc_grv[_idx2(i, k, pcols)])
        vlc_trb[i - 1] = wrk2
        vlc_dry[_idx2(i, k, pcols)] = wrk3


@export
def aero_model_drydep_select_branches_codon(
    apply_srf_drydep: int,
    branch_mask_p: cobj,
):
    branch_mask = Ptr[int](branch_mask_p)

    mask = 0
    if apply_srf_drydep != 0:
        mask |= 1

    branch_mask[0] = mask


def _modal_aero_bcscavcoef_get_core(
    m: int,
    ncol: int,
    pcols: int,
    pver: int,
    ntot_amode: int,
    nimptblgrow_mind: int,
    nimptblgrow_maxd: int,
    dlndg_nimptblgrow: float,
    dgnum_mode: float,
    isprx_mask: Ptr[int],
    dgn_awet: Ptr[float],
    scavimptblnum_mode: Ptr[float],
    scavimptblvol_mode: Ptr[float],
    scavcoefnum: Ptr[float],
    scavcoefvol: Ptr[float],
):
    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            if isprx_mask[_idx2(i, k, pcols)] != 0:
                dumdgratio = dgn_awet[_idx3(i, k, m, pcols, pver)] / dgnum_mode

                if dumdgratio >= 0.99 and dumdgratio <= 1.01:
                    tbl_idx = 0 - nimptblgrow_mind
                    scavimpvol = scavimptblvol_mode[tbl_idx]
                    scavimpnum = scavimptblnum_mode[tbl_idx]
                else:
                    xgrow = log(dumdgratio) / dlndg_nimptblgrow
                    jgrow = int(xgrow)
                    if xgrow < 0.0:
                        jgrow = jgrow - 1
                    if jgrow < nimptblgrow_mind:
                        jgrow = nimptblgrow_mind
                        xgrow = float(jgrow)
                    else:
                        jgrow = min(jgrow, nimptblgrow_maxd - 1)

                    dumfhi = xgrow - jgrow
                    dumflo = 1.0 - dumfhi
                    tbl_idx = jgrow - nimptblgrow_mind

                    scavimpvol = (
                        dumflo * scavimptblvol_mode[tbl_idx]
                        + dumfhi * scavimptblvol_mode[tbl_idx + 1]
                    )
                    scavimpnum = (
                        dumflo * scavimptblnum_mode[tbl_idx]
                        + dumfhi * scavimptblnum_mode[tbl_idx + 1]
                    )

                scavcoefvol[_idx2(i, k, pcols)] = exp(scavimpvol)
                scavcoefnum[_idx2(i, k, pcols)] = exp(scavimpnum)
            else:
                scavcoefvol[_idx2(i, k, pcols)] = 0.0
                scavcoefnum[_idx2(i, k, pcols)] = 0.0


@export
def modal_aero_bcscavcoef_get_codon(
    m: int,
    ncol: int,
    pcols: int,
    pver: int,
    ntot_amode: int,
    nimptblgrow_mind: int,
    nimptblgrow_maxd: int,
    dlndg_nimptblgrow: float,
    dgnum_mode: float,
    isprx_mask_p: cobj,
    dgn_awet_p: cobj,
    scavimptblnum_mode_p: cobj,
    scavimptblvol_mode_p: cobj,
    scavcoefnum_p: cobj,
    scavcoefvol_p: cobj,
):
    _modal_aero_bcscavcoef_get_core(
        m,
        ncol,
        pcols,
        pver,
        ntot_amode,
        nimptblgrow_mind,
        nimptblgrow_maxd,
        dlndg_nimptblgrow,
        dgnum_mode,
        Ptr[int](isprx_mask_p),
        Ptr[float](dgn_awet_p),
        Ptr[float](scavimptblnum_mode_p),
        Ptr[float](scavimptblvol_mode_p),
        Ptr[float](scavcoefnum_p),
        Ptr[float](scavcoefvol_p),
    )


@export
def aero_model_wetdep_inputs_codon(
    cam5_flag: int,
    ncol: int,
    pcols: int,
    pver: int,
    qliq_p: cobj,
    qice_p: cobj,
    icwmrdp_p: cobj,
    icwmrsh_p: cobj,
    rprddp_p: cobj,
    rprdsh_p: cobj,
    sh_frac_p: cobj,
    dp_frac_p: cobj,
    evapcsh_p: cobj,
    evapcdp_p: cobj,
    cldcu_p: cobj,
    evapc_p: cobj,
    cmfdqr_p: cobj,
    conicw_p: cobj,
    totcond_p: cobj,
):
    qliq = Ptr[float](qliq_p)
    qice = Ptr[float](qice_p)
    icwmrdp = Ptr[float](icwmrdp_p)
    icwmrsh = Ptr[float](icwmrsh_p)
    rprddp = Ptr[float](rprddp_p)
    rprdsh = Ptr[float](rprdsh_p)
    sh_frac = Ptr[float](sh_frac_p)
    dp_frac = Ptr[float](dp_frac_p)
    evapcsh = Ptr[float](evapcsh_p)
    evapcdp = Ptr[float](evapcdp_p)
    cldcu = Ptr[float](cldcu_p)
    evapc = Ptr[float](evapc_p)
    cmfdqr = Ptr[float](cmfdqr_p)
    conicw = Ptr[float](conicw_p)
    totcond = Ptr[float](totcond_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            cldcu[idx] = dp_frac[idx] + sh_frac[idx]
            evapc[idx] = evapcsh[idx] + evapcdp[idx]
            cmfdqr[idx] = rprddp[idx] + rprdsh[idx]

            if cam5_flag != 0:
                conicw[idx] = (
                    icwmrdp[idx] * dp_frac[idx] + icwmrsh[idx] * sh_frac[idx]
                ) / max(0.01, sh_frac[idx] + dp_frac[idx])
            else:
                conicw[idx] = icwmrdp[idx] + icwmrsh[idx]

            totcond[idx] = qliq[idx] + qice[idx]


@export
def aero_model_wetdep_precip_mask_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gravit: float,
    pdel_p: cobj,
    prain_p: cobj,
    cmfdqr_p: cobj,
    evapr_p: cobj,
    prec_p: cobj,
    isprx_mask_p: cobj,
):
    pdel = Ptr[float](pdel_p)
    prain = Ptr[float](prain_p)
    cmfdqr = Ptr[float](cmfdqr_p)
    evapr = Ptr[float](evapr_p)
    prec = Ptr[float](prec_p)
    isprx_mask = Ptr[int](isprx_mask_p)

    for i in range(1, ncol + 1):
        prec[i - 1] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            if prec[i - 1] >= 1.0e-7:
                isprx_mask[idx] = 1
            else:
                isprx_mask[idx] = 0

        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            prec[i - 1] = prec[i - 1] + (
                prain[idx] + cmfdqr[idx] - evapr[idx]
            ) * pdel[idx] / gravit


@export
def aero_model_wetdep_f_act_conv_coarse_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    dt: float,
    lcoardust: int,
    lcoarnacl: int,
    state_q_p: cobj,
    ptend_q_p: cobj,
    f_act_conv_coarse_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    ptend_q = Ptr[float](ptend_q_p)
    f_act_conv_coarse = Ptr[float](f_act_conv_coarse_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            f_act_conv_coarse[_idx2(i, k, pcols)] = 0.60

    if lcoardust <= 0 or lcoarnacl <= 0:
        return

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            tmpdust = state_q[_idx3(i, k, lcoardust, pcols, pver)] + dt * ptend_q[
                _idx3(i, k, lcoardust, pcols, pver)
            ]
            if tmpdust < 0.0:
                tmpdust = 0.0
            tmpnacl = state_q[_idx3(i, k, lcoarnacl, pcols, pver)] + dt * ptend_q[
                _idx3(i, k, lcoarnacl, pcols, pver)
            ]
            if tmpnacl < 0.0:
                tmpnacl = 0.0
            if tmpdust + tmpnacl > 1.0e-30:
                f_act_conv_coarse[_idx2(i, k, pcols)] = (
                    0.40 * tmpdust + 0.80 * tmpnacl
                ) / (tmpdust + tmpnacl)


def _aero_model_wetdep_column_flux_core(
    ncol: int,
    pcols: int,
    pver: int,
    gravit: float,
    field: Ptr[float],
    field_offset: int,
    pdel: Ptr[float],
    sflx: Ptr[float],
    sflx_offset: int,
):
    for i in range(1, ncol + 1):
        total = 0.0
        for k in range(1, pver + 1):
            total += field[field_offset + _idx2(i, k, pcols)] * pdel[_idx2(i, k, pcols)] / gravit
        sflx[sflx_offset + i - 1] = total


@export
def aero_model_wetdep_column_flux_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gravit: float,
    field_p: cobj,
    pdel_p: cobj,
    sflx_p: cobj,
):
    _aero_model_wetdep_column_flux_core(
        ncol,
        pcols,
        pver,
        gravit,
        Ptr[float](field_p),
        0,
        Ptr[float](pdel_p),
        Ptr[float](sflx_p),
        0,
    )


@export
def aero_model_wetdep_codon(
    stage: int,
    ncol: int,
    pcols: int,
    pver: int,
    dt: float,
    tmpa: float,
    gravit: float,
    pdel_p: cobj,
    state_tracer_p: cobj,
    ptend_tracer_p: cobj,
    q_tmp_p: cobj,
    dqdt_p: cobj,
    sflx_p: cobj,
    sflx_ics_p: cobj,
    sflx_iss_p: cobj,
    sflx_bcs_p: cobj,
    sflx_bss_p: cobj,
    hygro_sum_old_p: cobj,
    hygro_sum_del_p: cobj,
    qaerwat_p: cobj,
    fldcw_p: cobj,
    icscavt_p: cobj,
    isscavt_p: cobj,
    bcscavt_p: cobj,
    bsscavt_p: cobj,
    aerdep_p: cobj,
):
    pdel = Ptr[float](pdel_p)
    state_tracer = Ptr[float](state_tracer_p)
    ptend_tracer = Ptr[float](ptend_tracer_p)
    q_tmp = Ptr[float](q_tmp_p)
    dqdt = Ptr[float](dqdt_p)
    sflx = Ptr[float](sflx_p)
    sflx_ics = Ptr[float](sflx_ics_p)
    sflx_iss = Ptr[float](sflx_iss_p)
    sflx_bcs = Ptr[float](sflx_bcs_p)
    sflx_bss = Ptr[float](sflx_bss_p)
    hygro_sum_old = Ptr[float](hygro_sum_old_p)
    hygro_sum_del = Ptr[float](hygro_sum_del_p)
    qaerwat = Ptr[float](qaerwat_p)
    fldcw = Ptr[float](fldcw_p)
    icscavt = Ptr[float](icscavt_p)
    isscavt = Ptr[float](isscavt_p)
    bcscavt = Ptr[float](bcscavt_p)
    bsscavt = Ptr[float](bsscavt_p)
    aerdep = Ptr[float](aerdep_p)

    if stage == 1:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, pcols)
                q_tmp[idx] = state_tracer[idx] + ptend_tracer[idx] * dt
        return

    if stage == 2:
        tmpb = tmpa * dt
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, pcols)
                ptend_tracer[idx] = ptend_tracer[idx] + dqdt[idx]
                if tmpa != 0.0:
                    hygro_sum_old[idx] = hygro_sum_old[idx] + tmpa * q_tmp[idx]
                    hygro_sum_del[idx] = hygro_sum_del[idx] + tmpb * dqdt[idx]

        _aero_model_wetdep_column_flux_core(ncol, pcols, pver, gravit, dqdt, 0, pdel, sflx, 0)
        _aero_model_wetdep_column_flux_core(ncol, pcols, pver, gravit, icscavt, 0, pdel, sflx_ics, 0)
        _aero_model_wetdep_column_flux_core(ncol, pcols, pver, gravit, isscavt, 0, pdel, sflx_iss, 0)
        _aero_model_wetdep_column_flux_core(ncol, pcols, pver, gravit, bcscavt, 0, pdel, sflx_bcs, 0)
        _aero_model_wetdep_column_flux_core(ncol, pcols, pver, gravit, bsscavt, 0, pdel, sflx_bss, 0)

        for i in range(1, ncol + 1):
            aerdep[i - 1] = sflx[i - 1]
        return

    if stage == 3:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, pcols)
                water_old = max(0.0, qaerwat[idx])
                hygro_sum_old_ik = max(0.0, hygro_sum_old[idx])
                hygro_sum_new_ik = max(0.0, hygro_sum_old_ik + hygro_sum_del[idx])
                if hygro_sum_new_ik >= 10.0 * hygro_sum_old_ik:
                    water_new = 10.0 * water_old
                else:
                    water_new = water_old * (hygro_sum_new_ik / hygro_sum_old_ik)
                qaerwat[idx] = water_new
        return

    if stage == 4:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, pcols)
                fldcw[idx] = fldcw[idx] + dqdt[idx] * dt

        _aero_model_wetdep_column_flux_core(ncol, pcols, pver, gravit, dqdt, 0, pdel, sflx, 0)
        _aero_model_wetdep_column_flux_core(ncol, pcols, pver, gravit, icscavt, 0, pdel, sflx_ics, 0)
        _aero_model_wetdep_column_flux_core(ncol, pcols, pver, gravit, isscavt, 0, pdel, sflx_iss, 0)
        _aero_model_wetdep_column_flux_core(ncol, pcols, pver, gravit, bcscavt, 0, pdel, sflx_bcs, 0)
        _aero_model_wetdep_column_flux_core(ncol, pcols, pver, gravit, bsscavt, 0, pdel, sflx_bss, 0)

        for i in range(1, ncol + 1):
            aerdep[i - 1] = sflx[i - 1]


@export
def clddiag_codon(
    pcols: int,
    pver: int,
    ncol: int,
    tmelt: float,
    rair: float,
    gravit: float,
    convfw: float,
    t_p: cobj,
    pmid_p: cobj,
    pdel_p: cobj,
    cmfdqr_p: cobj,
    evapc_p: cobj,
    cldt_p: cobj,
    cldcu_p: cobj,
    cldst_p: cobj,
    cme_p: cobj,
    evapr_p: cobj,
    prain_p: cobj,
    cldv_p: cobj,
    cldvcu_p: cobj,
    cldvst_p: cobj,
    rain_p: cobj,
    sumppr_p: cobj,
    sumpppr_p: cobj,
    cldv1_p: cobj,
    sumppr_cu_p: cobj,
    sumpppr_cu_p: cobj,
    cldv1_cu_p: cobj,
    sumppr_st_p: cobj,
    sumpppr_st_p: cobj,
    cldv1_st_p: cobj,
):
    t = Ptr[float](t_p)
    pmid = Ptr[float](pmid_p)
    pdel = Ptr[float](pdel_p)
    cmfdqr = Ptr[float](cmfdqr_p)
    evapc = Ptr[float](evapc_p)
    cldt = Ptr[float](cldt_p)
    cldcu = Ptr[float](cldcu_p)
    cldst = Ptr[float](cldst_p)
    evapr = Ptr[float](evapr_p)
    prain = Ptr[float](prain_p)
    cldv = Ptr[float](cldv_p)
    cldvcu = Ptr[float](cldvcu_p)
    cldvst = Ptr[float](cldvst_p)
    rain = Ptr[float](rain_p)
    sumppr = Ptr[float](sumppr_p)
    sumpppr = Ptr[float](sumpppr_p)
    cldv1 = Ptr[float](cldv1_p)
    sumppr_cu = Ptr[float](sumppr_cu_p)
    sumpppr_cu = Ptr[float](sumpppr_cu_p)
    cldv1_cu = Ptr[float](cldv1_cu_p)
    sumppr_st = Ptr[float](sumppr_st_p)
    sumpppr_st = Ptr[float](sumpppr_st_p)
    cldv1_st = Ptr[float](cldv1_st_p)

    for i in range(1, ncol + 1):
        sumppr[i - 1] = 0.0
        cldv1[i - 1] = 0.0
        sumpppr[i - 1] = 1.0e-36
        sumppr_cu[i - 1] = 0.0
        cldv1_cu[i - 1] = 0.0
        sumpppr_cu[i - 1] = 1.0e-36
        sumppr_st[i - 1] = 0.0
        cldv1_st[i - 1] = 0.0
        sumpppr_st[i - 1] = 1.0e-36

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            cldv[_idx2(i, k, pcols)] = max(
                min(1.0, cldv1[i - 1] / sumpppr[i - 1]) * sumppr[i - 1] / sumpppr[i - 1],
                cldt[_idx2(i, k, pcols)],
            )
            lprec = (
                pdel[_idx2(i, k, pcols)] / gravit
                * (
                    prain[_idx2(i, k, pcols)]
                    + cmfdqr[_idx2(i, k, pcols)]
                    - evapr[_idx2(i, k, pcols)]
                )
            )
            lprecp = max(lprec, 1.0e-30)
            cldv1[i - 1] = cldv1[i - 1] + cldt[_idx2(i, k, pcols)] * lprecp
            sumppr[i - 1] = sumppr[i - 1] + lprec
            sumpppr[i - 1] = sumpppr[i - 1] + lprecp

            cldvcu[_idx2(i, k, pcols)] = max(
                min(1.0, cldv1_cu[i - 1] / sumpppr_cu[i - 1])
                * (sumppr_cu[i - 1] / sumpppr_cu[i - 1]),
                0.0,
            )
            lprec_cu = (
                pdel[_idx2(i, k, pcols)] / gravit
                * (cmfdqr[_idx2(i, k, pcols)] - evapc[_idx2(i, k, pcols)])
            )
            lprecp_cu = max(lprec_cu, 1.0e-30)
            cldv1_cu[i - 1] = cldv1_cu[i - 1] + cldcu[_idx2(i, k, pcols)] * lprecp_cu
            sumppr_cu[i - 1] = sumppr_cu[i - 1] + lprec_cu
            sumpppr_cu[i - 1] = sumpppr_cu[i - 1] + lprecp_cu

            cldvst[_idx2(i, k, pcols)] = max(
                min(1.0, cldv1_st[i - 1] / sumpppr_st[i - 1])
                * (sumppr_st[i - 1] / sumpppr_st[i - 1]),
                0.0,
            )
            lprec_st = (
                pdel[_idx2(i, k, pcols)] / gravit
                * (prain[_idx2(i, k, pcols)] - evapr[_idx2(i, k, pcols)])
            )
            lprecp_st = max(lprec_st, 1.0e-30)
            cldv1_st[i - 1] = cldv1_st[i - 1] + cldst[_idx2(i, k, pcols)] * lprecp_st
            sumppr_st[i - 1] = sumppr_st[i - 1] + lprec_st
            sumpppr_st[i - 1] = sumpppr_st[i - 1] + lprecp_st

            rain[_idx2(i, k, pcols)] = 0.0
            if t[_idx2(i, k, pcols)] > tmelt:
                rho = pmid[_idx2(i, k, pcols)] / (rair * t[_idx2(i, k, pcols)])
                vfall = convfw / sqrt(rho)
                rain[_idx2(i, k, pcols)] = sumppr[i - 1] / (rho * vfall)
                if rain[_idx2(i, k, pcols)] < 1.0e-14:
                    rain[_idx2(i, k, pcols)] = 0.0


def _wetdepa_v2_core(
    pcols: int,
    pver: int,
    ncol: int,
    branch_mode: int,
    gravit: float,
    deltat: float,
    omsm: float,
    sol_facti: float,
    sol_factb: float,
    p: Ptr[float],
    q: Ptr[float],
    pdel: Ptr[float],
    cldt: Ptr[float],
    cldc: Ptr[float],
    cmfdqr: Ptr[float],
    evapc: Ptr[float],
    conicw: Ptr[float],
    precs: Ptr[float],
    conds: Ptr[float],
    evaps: Ptr[float],
    cwat: Ptr[float],
    tracer: Ptr[float],
    tracer_offset: int,
    scavt: Ptr[float],
    scavt_offset: int,
    iscavt: Ptr[float],
    iscavt_offset: int,
    cldvcu: Ptr[float],
    cldvst: Ptr[float],
    dlf: Ptr[float],
    fracis: Ptr[float],
    fracis_offset: int,
    scavcoef: Ptr[float],
    scavcoef_offset: int,
    sol_factic: Ptr[float],
    sol_factic_offset: int,
    qqcw: Ptr[float],
    qqcw_offset: int,
    f_act_conv: Ptr[float],
    f_act_conv_offset: int,
    icscavt: Ptr[float],
    icscavt_offset: int,
    isscavt: Ptr[float],
    isscavt_offset: int,
    bcscavt: Ptr[float],
    bcscavt_offset: int,
    bsscavt: Ptr[float],
    bsscavt_offset: int,
    clds: Ptr[float],
    fracev: Ptr[float],
    fracev_cu: Ptr[float],
    fracp: Ptr[float],
    pdog: Ptr[float],
    rpdog: Ptr[float],
    precabc: Ptr[float],
    precabs: Ptr[float],
    rat: Ptr[float],
    scavab: Ptr[float],
    scavabc: Ptr[float],
    srcc: Ptr[float],
    srcs: Ptr[float],
    srct: Ptr[float],
    fins: Ptr[float],
    finc: Ptr[float],
    conv_scav_ic: Ptr[float],
    conv_scav_bc: Ptr[float],
    st_scav_ic: Ptr[float],
    st_scav_bc: Ptr[float],
    odds: Ptr[float],
    dblchek: Ptr[float],
    trac_qqcw: Ptr[float],
    tracer_incu: Ptr[float],
    tracer_mean: Ptr[float],
    dblchek_hist: Ptr[float],
    dblchek_hist_offset: int,
    srct_hist: Ptr[float],
    srct_hist_offset: int,
    rat_hist: Ptr[float],
    rat_hist_offset: int,
    fracev_hist: Ptr[float],
    fracev_hist_offset: int,
):
    for i in range(1, ncol + 1):
        precabs[i - 1] = 0.0
        precabc[i - 1] = 0.0
        scavab[i - 1] = 0.0
        scavabc[i - 1] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            clds[i - 1] = cldt[_idx2(i, k, pcols)] - cldc[_idx2(i, k, pcols)]
            pdog[i - 1] = pdel[_idx2(i, k, pcols)] / gravit
            rpdog[i - 1] = gravit / pdel[_idx2(i, k, pcols)]
            rdeltat = 1.0 / deltat

            fracev[i - 1] = (
                evaps[_idx2(i, k, pcols)] * pdog[i - 1] / max(1.0e-12, precabs[i - 1])
            )
            fracev[i - 1] = max(0.0, min(1.0, fracev[i - 1]))

            fracev_cu[i - 1] = (
                evapc[_idx2(i, k, pcols)] * pdog[i - 1] / max(1.0e-12, precabc[i - 1])
            )
            fracev_cu[i - 1] = max(0.0, min(1.0, fracev_cu[i - 1]))

            fracp[i - 1] = (
                cmfdqr[_idx2(i, k, pcols)]
                * deltat
                / max(
                    1.0e-12,
                    cldc[_idx2(i, k, pcols)] * conicw[_idx2(i, k, pcols)]
                    + (cmfdqr[_idx2(i, k, pcols)] + dlf[_idx2(i, k, pcols)]) * deltat,
                )
            )
            fracp[i - 1] = max(min(1.0, fracp[i - 1]), 0.0)

            if branch_mode != 0:
                if branch_mode == 1:
                    conv_scav_ic[i - 1] = 0.0
                    conv_scav_bc[i - 1] = 0.0

                    fracp[i - 1] = (
                        precs[_idx2(i, k, pcols)]
                        * deltat
                        / max(
                            1.0e-12,
                            cwat[_idx2(i, k, pcols)] + precs[_idx2(i, k, pcols)] * deltat,
                        )
                    )
                    fracp[i - 1] = max(0.0, min(1.0, fracp[i - 1]))
                    st_scav_ic[i - 1] = (
                        sol_facti
                        * fracp[i - 1]
                        * tracer[tracer_offset + _idx2(i, k, pcols)]
                        * rdeltat
                    )
                    st_scav_bc[i - 1] = 0.0
                else:
                    trac_qqcw[i - 1] = min(
                        qqcw[qqcw_offset + _idx2(i, k, pcols)],
                        tracer[tracer_offset + _idx2(i, k, pcols)]
                        * (clds[i - 1] / max(0.01, 1.0 - clds[i - 1])),
                    )
                    tracer_incu[i - 1] = (
                        f_act_conv[f_act_conv_offset + _idx2(i, k, pcols)]
                        * (tracer[tracer_offset + _idx2(i, k, pcols)] + trac_qqcw[i - 1])
                    )
                    conv_scav_ic[i - 1] = (
                        sol_factic[sol_factic_offset + _idx2(i, k, pcols)]
                        * cldc[_idx2(i, k, pcols)]
                        * fracp[i - 1]
                        * tracer_incu[i - 1]
                        * rdeltat
                    )
                    tracer_mean[i - 1] = (
                        tracer[tracer_offset + _idx2(i, k, pcols)]
                        * (
                            1.0
                            - cldc[_idx2(i, k, pcols)]
                            * f_act_conv[f_act_conv_offset + _idx2(i, k, pcols)]
                        )
                        - cldc[_idx2(i, k, pcols)]
                        * f_act_conv[f_act_conv_offset + _idx2(i, k, pcols)]
                        * trac_qqcw[i - 1]
                    )
                    tracer_mean[i - 1] = max(0.0, tracer_mean[i - 1])

                    odds[i - 1] = (
                        precabc[i - 1]
                        / max(cldvcu[_idx2(i, k, pcols)], 1.0e-5)
                        * scavcoef[scavcoef_offset + _idx2(i, k, pcols)]
                        * deltat
                    )
                    odds[i - 1] = max(min(1.0, odds[i - 1]), 0.0)
                    conv_scav_bc[i - 1] = (
                        sol_factb
                        * cldvcu[_idx2(i, k, pcols)]
                        * odds[i - 1]
                        * tracer_mean[i - 1]
                        * rdeltat
                    )

                    st_scav_ic[i - 1] = 0.0

                    odds[i - 1] = (
                        precabs[i - 1]
                        / max(cldvst[_idx2(i, k, pcols)], 1.0e-5)
                        * scavcoef[scavcoef_offset + _idx2(i, k, pcols)]
                        * deltat
                    )
                    odds[i - 1] = max(min(1.0, odds[i - 1]), 0.0)
                    st_scav_bc[i - 1] = (
                        sol_factb
                        * cldvst[_idx2(i, k, pcols)]
                        * odds[i - 1]
                        * tracer_mean[i - 1]
                        * rdeltat
                    )
            else:
                conv_scav_ic[i - 1] = (
                    sol_factic[sol_factic_offset + _idx2(i, k, pcols)]
                    * cldc[_idx2(i, k, pcols)]
                    * fracp[i - 1]
                    * tracer[tracer_offset + _idx2(i, k, pcols)]
                    * rdeltat
                )

                odds[i - 1] = (
                    precabc[i - 1]
                    / max(cldvcu[_idx2(i, k, pcols)], 1.0e-5)
                    * scavcoef[scavcoef_offset + _idx2(i, k, pcols)]
                    * deltat
                )
                odds[i - 1] = max(min(1.0, odds[i - 1]), 0.0)
                conv_scav_bc[i - 1] = (
                    sol_factb
                    * cldvcu[_idx2(i, k, pcols)]
                    * odds[i - 1]
                    * tracer[tracer_offset + _idx2(i, k, pcols)]
                    * rdeltat
                )

                fracp[i - 1] = (
                    precs[_idx2(i, k, pcols)]
                    * deltat
                    / max(
                        1.0e-12,
                        cwat[_idx2(i, k, pcols)] + precs[_idx2(i, k, pcols)] * deltat,
                    )
                )
                fracp[i - 1] = max(0.0, min(1.0, fracp[i - 1]))
                st_scav_ic[i - 1] = (
                    sol_facti
                    * clds[i - 1]
                    * fracp[i - 1]
                    * tracer[tracer_offset + _idx2(i, k, pcols)]
                    * rdeltat
                )

                odds[i - 1] = (
                    precabs[i - 1]
                    / max(cldvst[_idx2(i, k, pcols)], 1.0e-5)
                    * scavcoef[scavcoef_offset + _idx2(i, k, pcols)]
                    * deltat
                )
                odds[i - 1] = max(min(1.0, odds[i - 1]), 0.0)
                st_scav_bc[i - 1] = (
                    sol_factb
                    * (cldvst[_idx2(i, k, pcols)] * odds[i - 1])
                    * tracer[tracer_offset + _idx2(i, k, pcols)]
                    * rdeltat
                )

            srcc[i - 1] = conv_scav_ic[i - 1] + conv_scav_bc[i - 1]
            finc[i - 1] = conv_scav_ic[i - 1] / (srcc[i - 1] + 1.0e-36)

            srcs[i - 1] = st_scav_ic[i - 1] + st_scav_bc[i - 1]
            fins[i - 1] = st_scav_ic[i - 1] / (srcs[i - 1] + 1.0e-36)

            rat[i - 1] = (
                tracer[tracer_offset + _idx2(i, k, pcols)]
                / max(deltat * (srcc[i - 1] + srcs[i - 1]), 1.0e-36)
            )
            if rat[i - 1] < 1.0:
                srcs[i - 1] = srcs[i - 1] * rat[i - 1]
                srcc[i - 1] = srcc[i - 1] * rat[i - 1]
            srct[i - 1] = (srcc[i - 1] + srcs[i - 1]) * omsm

            fracp[i - 1] = (
                deltat
                * srct[i - 1]
                / max(
                    cldvst[_idx2(i, k, pcols)] * tracer[tracer_offset + _idx2(i, k, pcols)],
                    1.0e-36,
                )
            )
            fracp[i - 1] = max(0.0, min(1.0, fracp[i - 1]))
            fracis[fracis_offset + _idx2(i, k, pcols)] = 1.0 - fracp[i - 1]

            scavt[scavt_offset + _idx2(i, k, pcols)] = -srct[i - 1] + (
                fracev[i - 1] * scavab[i - 1] + fracev_cu[i - 1] * scavabc[i - 1]
            ) * rpdog[i - 1]
            iscavt[iscavt_offset + _idx2(i, k, pcols)] = (
                -(srcc[i - 1] * finc[i - 1] + srcs[i - 1] * fins[i - 1]) * omsm
            )

            icscavt[icscavt_offset + _idx2(i, k, pcols)] = -(srcc[i - 1] * finc[i - 1]) * omsm
            isscavt[isscavt_offset + _idx2(i, k, pcols)] = -(srcs[i - 1] * fins[i - 1]) * omsm
            bcscavt[bcscavt_offset + _idx2(i, k, pcols)] = (
                -(srcc[i - 1] * (1.0 - finc[i - 1])) * omsm
                + fracev_cu[i - 1] * scavabc[i - 1] * rpdog[i - 1]
            )
            bsscavt[bsscavt_offset + _idx2(i, k, pcols)] = (
                -(srcs[i - 1] * (1.0 - fins[i - 1])) * omsm
                + fracev[i - 1] * scavab[i - 1] * rpdog[i - 1]
            )

            dblchek[i - 1] = (
                tracer[tracer_offset + _idx2(i, k, pcols)]
                + deltat * scavt[scavt_offset + _idx2(i, k, pcols)]
            )

            scavab[i - 1] = scavab[i - 1] * (1.0 - fracev[i - 1]) + srcs[i - 1] * pdog[i - 1]
            precabs[i - 1] = precabs[i - 1] + (
                precs[_idx2(i, k, pcols)] - evaps[_idx2(i, k, pcols)]
            ) * pdog[i - 1]
            scavabc[i - 1] = scavabc[i - 1] * (1.0 - fracev_cu[i - 1]) + srcc[i - 1] * pdog[i - 1]
            precabc[i - 1] = precabc[i - 1] + (
                cmfdqr[_idx2(i, k, pcols)] - evapc[_idx2(i, k, pcols)]
            ) * pdog[i - 1]

            dblchek_hist[dblchek_hist_offset + _idx2(i, k, pcols)] = dblchek[i - 1]
            srct_hist[srct_hist_offset + _idx2(i, k, pcols)] = srct[i - 1]
            rat_hist[rat_hist_offset + _idx2(i, k, pcols)] = rat[i - 1]
            fracev_hist[fracev_hist_offset + _idx2(i, k, pcols)] = fracev[i - 1]


@export
def wetdepa_v2_codon(
    pcols: int,
    pver: int,
    ncol: int,
    branch_mode: int,
    gravit: float,
    deltat: float,
    omsm: float,
    sol_facti: float,
    sol_factb: float,
    p_p: cobj,
    q_p: cobj,
    pdel_p: cobj,
    cldt_p: cobj,
    cldc_p: cobj,
    cmfdqr_p: cobj,
    evapc_p: cobj,
    conicw_p: cobj,
    precs_p: cobj,
    conds_p: cobj,
    evaps_p: cobj,
    cwat_p: cobj,
    tracer_p: cobj,
    scavt_p: cobj,
    iscavt_p: cobj,
    cldvcu_p: cobj,
    cldvst_p: cobj,
    dlf_p: cobj,
    fracis_p: cobj,
    scavcoef_p: cobj,
    sol_factic_p: cobj,
    qqcw_p: cobj,
    f_act_conv_p: cobj,
    icscavt_p: cobj,
    isscavt_p: cobj,
    bcscavt_p: cobj,
    bsscavt_p: cobj,
    clds_p: cobj,
    fracev_p: cobj,
    fracev_cu_p: cobj,
    fracp_p: cobj,
    pdog_p: cobj,
    rpdog_p: cobj,
    precabc_p: cobj,
    precabs_p: cobj,
    rat_p: cobj,
    scavab_p: cobj,
    scavabc_p: cobj,
    srcc_p: cobj,
    srcs_p: cobj,
    srct_p: cobj,
    fins_p: cobj,
    finc_p: cobj,
    conv_scav_ic_p: cobj,
    conv_scav_bc_p: cobj,
    st_scav_ic_p: cobj,
    st_scav_bc_p: cobj,
    odds_p: cobj,
    dblchek_p: cobj,
    trac_qqcw_p: cobj,
    tracer_incu_p: cobj,
    tracer_mean_p: cobj,
    dblchek_hist_p: cobj,
    srct_hist_p: cobj,
    rat_hist_p: cobj,
    fracev_hist_p: cobj,
):
    _wetdepa_v2_core(
        pcols,
        pver,
        ncol,
        branch_mode,
        gravit,
        deltat,
        omsm,
        sol_facti,
        sol_factb,
        Ptr[float](p_p),
        Ptr[float](q_p),
        Ptr[float](pdel_p),
        Ptr[float](cldt_p),
        Ptr[float](cldc_p),
        Ptr[float](cmfdqr_p),
        Ptr[float](evapc_p),
        Ptr[float](conicw_p),
        Ptr[float](precs_p),
        Ptr[float](conds_p),
        Ptr[float](evaps_p),
        Ptr[float](cwat_p),
        Ptr[float](tracer_p),
        0,
        Ptr[float](scavt_p),
        0,
        Ptr[float](iscavt_p),
        0,
        Ptr[float](cldvcu_p),
        Ptr[float](cldvst_p),
        Ptr[float](dlf_p),
        Ptr[float](fracis_p),
        0,
        Ptr[float](scavcoef_p),
        0,
        Ptr[float](sol_factic_p),
        0,
        Ptr[float](qqcw_p),
        0,
        Ptr[float](f_act_conv_p),
        0,
        Ptr[float](icscavt_p),
        0,
        Ptr[float](isscavt_p),
        0,
        Ptr[float](bcscavt_p),
        0,
        Ptr[float](bsscavt_p),
        0,
        Ptr[float](clds_p),
        Ptr[float](fracev_p),
        Ptr[float](fracev_cu_p),
        Ptr[float](fracp_p),
        Ptr[float](pdog_p),
        Ptr[float](rpdog_p),
        Ptr[float](precabc_p),
        Ptr[float](precabs_p),
        Ptr[float](rat_p),
        Ptr[float](scavab_p),
        Ptr[float](scavabc_p),
        Ptr[float](srcc_p),
        Ptr[float](srcs_p),
        Ptr[float](srct_p),
        Ptr[float](fins_p),
        Ptr[float](finc_p),
        Ptr[float](conv_scav_ic_p),
        Ptr[float](conv_scav_bc_p),
        Ptr[float](st_scav_ic_p),
        Ptr[float](st_scav_bc_p),
        Ptr[float](odds_p),
        Ptr[float](dblchek_p),
        Ptr[float](trac_qqcw_p),
        Ptr[float](tracer_incu_p),
        Ptr[float](tracer_mean_p),
        Ptr[float](dblchek_hist_p),
        0,
        Ptr[float](srct_hist_p),
        0,
        Ptr[float](rat_hist_p),
        0,
        Ptr[float](fracev_hist_p),
        0,
    )


@export
def aero_model_wetdep_mode_phase_codon(
    m: int,
    lphase: int,
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    ntot_amode: int,
    nslot_max: int,
    nimptblgrow_mind: int,
    nimptblgrow_maxd: int,
    dt: float,
    gravit: float,
    omsm: float,
    dgnum_mode: float,
    dlndg_nimptblgrow: float,
    sol_factb: float,
    sol_facti: float,
    sol_factic_scalar: float,
    base_f_act_scalar: float,
    is_coarse_interstitial: int,
    f_act_conv_coarse_dust: float,
    f_act_conv_coarse_nacl: float,
    pmid_p: cobj,
    q1_p: cobj,
    pdel_p: cobj,
    cldt_p: cobj,
    cldcu_p: cobj,
    cmfdqr_p: cobj,
    evapc_p: cobj,
    conicw_p: cobj,
    prain_p: cobj,
    qme_p: cobj,
    evapr_p: cobj,
    totcond_p: cobj,
    cldvcu_p: cobj,
    cldvst_p: cobj,
    dlf_p: cobj,
    isprx_mask_p: cobj,
    dgnumwet_p: cobj,
    scavimptblnum_mode_p: cobj,
    scavimptblvol_mode_p: cobj,
    state_q_p: cobj,
    ptend_q_p: cobj,
    qaerwat_mode_p: cobj,
    fracis_full_p: cobj,
    f_act_conv_coarse_p: cobj,
    qqcw_mode_phase_p: cobj,
    q_tmp_work_p: cobj,
    hygro_sum_old_p: cobj,
    hygro_sum_del_p: cobj,
    scavcoefnum_p: cobj,
    scavcoefvol_p: cobj,
    scavcoef_work_p: cobj,
    iscavt_work_p: cobj,
    f_act_conv_work_p: cobj,
    sol_factic_work_p: cobj,
    slot_active_p: cobj,
    slot_mm_p: cobj,
    slot_jnv_p: cobj,
    slot_mass_kind_p: cobj,
    slot_hygro_scale_p: cobj,
    diag_dqdt_p: cobj,
    diag_icscavt_p: cobj,
    diag_isscavt_p: cobj,
    diag_bcscavt_p: cobj,
    diag_bsscavt_p: cobj,
    diag_sflx_p: cobj,
    diag_sflx_ics_p: cobj,
    diag_sflx_iss_p: cobj,
    diag_sflx_bcs_p: cobj,
    diag_sflx_bss_p: cobj,
    clds_p: cobj,
    fracev_p: cobj,
    fracev_cu_p: cobj,
    fracp_p: cobj,
    pdog_p: cobj,
    rpdog_p: cobj,
    precabc_p: cobj,
    precabs_p: cobj,
    rat_p: cobj,
    scavab_p: cobj,
    scavabc_p: cobj,
    srcc_p: cobj,
    srcs_p: cobj,
    srct_p: cobj,
    fins_p: cobj,
    finc_p: cobj,
    conv_scav_ic_p: cobj,
    conv_scav_bc_p: cobj,
    st_scav_ic_p: cobj,
    st_scav_bc_p: cobj,
    odds_p: cobj,
    dblchek_p: cobj,
    trac_qqcw_p: cobj,
    tracer_incu_p: cobj,
    tracer_mean_p: cobj,
    fracis_dummy_p: cobj,
    dblchek_hist_p: cobj,
    srct_hist_p: cobj,
    rat_hist_p: cobj,
    fracev_hist_p: cobj,
):
    pmid = Ptr[float](pmid_p)
    q1 = Ptr[float](q1_p)
    pdel = Ptr[float](pdel_p)
    cldt = Ptr[float](cldt_p)
    cldcu = Ptr[float](cldcu_p)
    cmfdqr = Ptr[float](cmfdqr_p)
    evapc = Ptr[float](evapc_p)
    conicw = Ptr[float](conicw_p)
    prain = Ptr[float](prain_p)
    qme = Ptr[float](qme_p)
    evapr = Ptr[float](evapr_p)
    totcond = Ptr[float](totcond_p)
    cldvcu = Ptr[float](cldvcu_p)
    cldvst = Ptr[float](cldvst_p)
    dlf = Ptr[float](dlf_p)
    isprx_mask = Ptr[int](isprx_mask_p)
    dgnumwet = Ptr[float](dgnumwet_p)
    scavimptblnum_mode = Ptr[float](scavimptblnum_mode_p)
    scavimptblvol_mode = Ptr[float](scavimptblvol_mode_p)
    state_q = Ptr[float](state_q_p)
    ptend_q = Ptr[float](ptend_q_p)
    qaerwat_mode = Ptr[float](qaerwat_mode_p)
    fracis_full = Ptr[float](fracis_full_p)
    f_act_conv_coarse = Ptr[float](f_act_conv_coarse_p)
    qqcw_mode_phase = Ptr[float](qqcw_mode_phase_p)
    q_tmp_work = Ptr[float](q_tmp_work_p)
    hygro_sum_old = Ptr[float](hygro_sum_old_p)
    hygro_sum_del = Ptr[float](hygro_sum_del_p)
    scavcoefnum = Ptr[float](scavcoefnum_p)
    scavcoefvol = Ptr[float](scavcoefvol_p)
    scavcoef_work = Ptr[float](scavcoef_work_p)
    iscavt_work = Ptr[float](iscavt_work_p)
    f_act_conv_work = Ptr[float](f_act_conv_work_p)
    sol_factic_work = Ptr[float](sol_factic_work_p)
    slot_active = Ptr[int](slot_active_p)
    slot_mm = Ptr[int](slot_mm_p)
    slot_jnv = Ptr[int](slot_jnv_p)
    slot_mass_kind = Ptr[int](slot_mass_kind_p)
    slot_hygro_scale = Ptr[float](slot_hygro_scale_p)
    diag_dqdt = Ptr[float](diag_dqdt_p)
    diag_icscavt = Ptr[float](diag_icscavt_p)
    diag_isscavt = Ptr[float](diag_isscavt_p)
    diag_bcscavt = Ptr[float](diag_bcscavt_p)
    diag_bsscavt = Ptr[float](diag_bsscavt_p)
    diag_sflx = Ptr[float](diag_sflx_p)
    diag_sflx_ics = Ptr[float](diag_sflx_ics_p)
    diag_sflx_iss = Ptr[float](diag_sflx_iss_p)
    diag_sflx_bcs = Ptr[float](diag_sflx_bcs_p)
    diag_sflx_bss = Ptr[float](diag_sflx_bss_p)
    clds = Ptr[float](clds_p)
    fracev = Ptr[float](fracev_p)
    fracev_cu = Ptr[float](fracev_cu_p)
    fracp = Ptr[float](fracp_p)
    pdog = Ptr[float](pdog_p)
    rpdog = Ptr[float](rpdog_p)
    precabc = Ptr[float](precabc_p)
    precabs = Ptr[float](precabs_p)
    rat = Ptr[float](rat_p)
    scavab = Ptr[float](scavab_p)
    scavabc = Ptr[float](scavabc_p)
    srcc = Ptr[float](srcc_p)
    srcs = Ptr[float](srcs_p)
    srct = Ptr[float](srct_p)
    fins = Ptr[float](fins_p)
    finc = Ptr[float](finc_p)
    conv_scav_ic = Ptr[float](conv_scav_ic_p)
    conv_scav_bc = Ptr[float](conv_scav_bc_p)
    st_scav_ic = Ptr[float](st_scav_ic_p)
    st_scav_bc = Ptr[float](st_scav_bc_p)
    odds = Ptr[float](odds_p)
    dblchek = Ptr[float](dblchek_p)
    trac_qqcw = Ptr[float](trac_qqcw_p)
    tracer_incu = Ptr[float](tracer_incu_p)
    tracer_mean = Ptr[float](tracer_mean_p)
    fracis_dummy = Ptr[float](fracis_dummy_p)
    dblchek_hist = Ptr[float](dblchek_hist_p)
    srct_hist = Ptr[float](srct_hist_p)
    rat_hist = Ptr[float](rat_hist_p)
    fracev_hist = Ptr[float](fracev_hist_p)

    if lphase == 1:
        for k in range(1, pver + 1):
            for i in range(1, pcols + 1):
                idx2 = _idx2(i, k, pcols)
                hygro_sum_old[idx2] = 0.0
                hygro_sum_del[idx2] = 0.0
        _modal_aero_bcscavcoef_get_core(
            m,
            ncol,
            pcols,
            pver,
            ntot_amode,
            nimptblgrow_mind,
            nimptblgrow_maxd,
            dlndg_nimptblgrow,
            dgnum_mode,
            isprx_mask,
            dgnumwet,
            scavimptblnum_mode,
            scavimptblvol_mode,
            scavcoefnum,
            scavcoefvol,
        )

    for slot in range(1, nslot_max + 1):
        slot_offset = (slot - 1) * pcols * pver
        sflx_offset = (slot - 1) * pcols

        for k in range(1, pver + 1):
            for i in range(1, pcols + 1):
                idx2 = _idx2(i, k, pcols)
                diag_dqdt[slot_offset + idx2] = 0.0
                diag_icscavt[slot_offset + idx2] = 0.0
                diag_isscavt[slot_offset + idx2] = 0.0
                diag_bcscavt[slot_offset + idx2] = 0.0
                diag_bsscavt[slot_offset + idx2] = 0.0
        for i in range(1, pcols + 1):
            diag_sflx[sflx_offset + i - 1] = 0.0
            diag_sflx_ics[sflx_offset + i - 1] = 0.0
            diag_sflx_iss[sflx_offset + i - 1] = 0.0
            diag_sflx_bcs[sflx_offset + i - 1] = 0.0
            diag_sflx_bss[sflx_offset + i - 1] = 0.0

        if slot_active[slot - 1] == 0:
            continue

        mm = slot_mm[slot - 1]
        mm_offset = (mm - 1) * pcols * pver
        jnv = slot_jnv[slot - 1]
        hygro_scale = slot_hygro_scale[slot - 1]
        mass_kind = slot_mass_kind[slot - 1]

        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx2 = _idx2(i, k, pcols)
                sol_factic_work[idx2] = sol_factic_scalar
                if jnv == 1:
                    scavcoef_work[idx2] = scavcoefnum[idx2]
                elif jnv == 2:
                    scavcoef_work[idx2] = scavcoefvol[idx2]
                else:
                    scavcoef_work[idx2] = 0.0

                if is_coarse_interstitial != 0:
                    if mass_kind == 1:
                        f_act_conv_work[idx2] = f_act_conv_coarse_dust
                    elif mass_kind == 2:
                        f_act_conv_work[idx2] = f_act_conv_coarse_nacl
                    else:
                        f_act_conv_work[idx2] = f_act_conv_coarse[idx2]
                else:
                    f_act_conv_work[idx2] = base_f_act_scalar

        if lphase == 1:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    idx2 = _idx2(i, k, pcols)
                    q_tmp_work[idx2] = state_q[mm_offset + idx2] + ptend_q[mm_offset + idx2] * dt

            _wetdepa_v2_core(
                pcols,
                pver,
                ncol,
                2,
                gravit,
                dt,
                omsm,
                sol_facti,
                sol_factb,
                pmid,
                q1,
                pdel,
                cldt,
                cldcu,
                cmfdqr,
                evapc,
                conicw,
                prain,
                qme,
                evapr,
                totcond,
                q_tmp_work,
                0,
                diag_dqdt,
                slot_offset,
                iscavt_work,
                0,
                cldvcu,
                cldvst,
                dlf,
                fracis_full,
                mm_offset,
                scavcoef_work,
                0,
                sol_factic_work,
                0,
                qqcw_mode_phase,
                slot_offset,
                f_act_conv_work,
                0,
                diag_icscavt,
                slot_offset,
                diag_isscavt,
                slot_offset,
                diag_bcscavt,
                slot_offset,
                diag_bsscavt,
                slot_offset,
                clds,
                fracev,
                fracev_cu,
                fracp,
                pdog,
                rpdog,
                precabc,
                precabs,
                rat,
                scavab,
                scavabc,
                srcc,
                srcs,
                srct,
                fins,
                finc,
                conv_scav_ic,
                conv_scav_bc,
                st_scav_ic,
                st_scav_bc,
                odds,
                dblchek,
                trac_qqcw,
                tracer_incu,
                tracer_mean,
                dblchek_hist,
                0,
                srct_hist,
                0,
                rat_hist,
                0,
                fracev_hist,
                0,
            )

            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    idx2 = _idx2(i, k, pcols)
                    ptend_q[mm_offset + idx2] = ptend_q[mm_offset + idx2] + diag_dqdt[slot_offset + idx2]
                    if hygro_scale != 0.0:
                        hygro_sum_old[idx2] = hygro_sum_old[idx2] + hygro_scale * q_tmp_work[idx2]
                        hygro_sum_del[idx2] = hygro_sum_del[idx2] + hygro_scale * dt * diag_dqdt[slot_offset + idx2]
        else:
            _wetdepa_v2_core(
                pcols,
                pver,
                ncol,
                1,
                gravit,
                dt,
                omsm,
                sol_facti,
                sol_factb,
                pmid,
                q1,
                pdel,
                cldt,
                cldcu,
                cmfdqr,
                evapc,
                conicw,
                prain,
                qme,
                evapr,
                totcond,
                qqcw_mode_phase,
                slot_offset,
                diag_dqdt,
                slot_offset,
                iscavt_work,
                0,
                cldvcu,
                cldvst,
                dlf,
                fracis_dummy,
                0,
                scavcoef_work,
                0,
                sol_factic_work,
                0,
                qqcw_mode_phase,
                slot_offset,
                f_act_conv_work,
                0,
                diag_icscavt,
                slot_offset,
                diag_isscavt,
                slot_offset,
                diag_bcscavt,
                slot_offset,
                diag_bsscavt,
                slot_offset,
                clds,
                fracev,
                fracev_cu,
                fracp,
                pdog,
                rpdog,
                precabc,
                precabs,
                rat,
                scavab,
                scavabc,
                srcc,
                srcs,
                srct,
                fins,
                finc,
                conv_scav_ic,
                conv_scav_bc,
                st_scav_ic,
                st_scav_bc,
                odds,
                dblchek,
                trac_qqcw,
                tracer_incu,
                tracer_mean,
                dblchek_hist,
                0,
                srct_hist,
                0,
                rat_hist,
                0,
                fracev_hist,
                0,
            )

            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    idx2 = _idx2(i, k, pcols)
                    qqcw_mode_phase[slot_offset + idx2] = (
                        qqcw_mode_phase[slot_offset + idx2] + diag_dqdt[slot_offset + idx2] * dt
                    )

        _aero_model_wetdep_column_flux_core(
            ncol, pcols, pver, gravit, diag_dqdt, slot_offset, pdel, diag_sflx, sflx_offset
        )
        _aero_model_wetdep_column_flux_core(
            ncol, pcols, pver, gravit, diag_icscavt, slot_offset, pdel, diag_sflx_ics, sflx_offset
        )
        _aero_model_wetdep_column_flux_core(
            ncol, pcols, pver, gravit, diag_isscavt, slot_offset, pdel, diag_sflx_iss, sflx_offset
        )
        _aero_model_wetdep_column_flux_core(
            ncol, pcols, pver, gravit, diag_bcscavt, slot_offset, pdel, diag_sflx_bcs, sflx_offset
        )
        _aero_model_wetdep_column_flux_core(
            ncol, pcols, pver, gravit, diag_bsscavt, slot_offset, pdel, diag_sflx_bss, sflx_offset
        )

    if lphase == 1:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx2 = _idx2(i, k, pcols)
                water_old = max(0.0, qaerwat_mode[idx2])
                hygro_sum_old_ik = max(0.0, hygro_sum_old[idx2])
                hygro_sum_new_ik = max(0.0, hygro_sum_old_ik + hygro_sum_del[idx2])
                if hygro_sum_new_ik >= 10.0 * hygro_sum_old_ik:
                    water_new = 10.0 * water_old
                else:
                    water_new = water_old * (hygro_sum_new_ik / hygro_sum_old_ik)
                qaerwat_mode[idx2] = water_new


@export
def calcram_codon(
    ncol: int,
    pcols: int,
    rair: float,
    gravit: float,
    ram1in_p: cobj,
    fvin_p: cobj,
    ram1_p: cobj,
    fv_p: cobj,
    obklen_p: cobj,
    ustar_p: cobj,
    landfrac_p: cobj,
    icefrac_p: cobj,
    ocnfrac_p: cobj,
    t_p: cobj,
    pmid_p: cobj,
    pdel_p: cobj,
):
    ram1in = Ptr[float](ram1in_p)
    fvin = Ptr[float](fvin_p)
    ram1 = Ptr[float](ram1_p)
    fv = Ptr[float](fv_p)
    obklen = Ptr[float](obklen_p)
    ustar = Ptr[float](ustar_p)
    landfrac = Ptr[float](landfrac_p)
    icefrac = Ptr[float](icefrac_p)
    ocnfrac = Ptr[float](ocnfrac_p)
    t = Ptr[float](t_p)
    pmid = Ptr[float](pmid_p)
    pdel = Ptr[float](pdel_p)

    zzocen = 0.0001
    zzsice = 0.0400
    xkar = 0.4

    for i in range(1, ncol + 1):
        z = pdel[i - 1] * rair * t[i - 1] / pmid[i - 1] / gravit / 2.0
        if obklen[i - 1] == 0.0:
            psi = 0.0
            psi0 = 0.0
        else:
            psi = min(max(z / obklen[i - 1], -1.0), 1.0)
            psi0 = min(max(zzocen / obklen[i - 1], -1.0), 1.0)

        temp = z / zzocen
        if icefrac[i - 1] > 0.5:
            if obklen[i - 1] > 0.0:
                psi0 = min(max(zzsice / obklen[i - 1], -1.0), 1.0)
            else:
                psi0 = 0.0
            temp = z / zzsice

        if psi > 0.0:
            ram = 1.0 / xkar / ustar[i - 1] * (log(temp) + 4.7 * (psi - psi0))
        else:
            nu = (1.0 - 15.0 * psi) ** 0.25
            nu0 = (1.0 - 15.0 * psi0) ** 0.25
            if ustar[i - 1] != 0.0:
                ram = 1.0 / xkar / ustar[i - 1] * (
                    log(temp)
                    + log(
                        ((nu0**2.0 + 1.0) * (nu0 + 1.0) ** 2.0)
                        / ((nu**2.0 + 1.0) * (nu + 1.0) ** 2.0)
                    )
                    + 2.0 * (atan(nu) - atan(nu0))
                )
            else:
                ram = 0.0

        if landfrac[i - 1] < 0.000000001:
            fv[i - 1] = ustar[i - 1]
            ram1[i - 1] = ram
        else:
            fv[i - 1] = fvin[i - 1]
            ram1[i - 1] = ram1in[i - 1]

    for i in range(1, ncol + 1):
        if fv[i - 1] == 0.0:
            fv[i - 1] = 1.0e-12


def _dust_cfint2(
    ncol: int,
    pcols: int,
    pverp: int,
    xin_k: int,
    x_p: cobj,
    f_p: cobj,
    fdot_p: cobj,
    xxk_p: cobj,
    fxdot_p: cobj,
    fxdd_p: cobj,
    psistar_p: cobj,
    xins_p: cobj,
    intz_p: cobj,
    status_p: cobj,
    fail_i_p: cobj,
    fail_k_p: cobj,
):
    x = Ptr[float](x_p)
    f = Ptr[float](f_p)
    fdot = Ptr[float](fdot_p)
    xxk = Ptr[float](xxk_p)
    fxdot = Ptr[float](fxdot_p)
    fxdd = Ptr[float](fxdd_p)
    psistar = Ptr[float](psistar_p)
    xins = Ptr[float](xins_p)
    intz = Ptr[int](intz_p)
    status = Ptr[int](status_p)
    fail_i = Ptr[int](fail_i_p)
    fail_k = Ptr[int](fail_k_p)

    for i in range(1, ncol + 1):
        xins[i - 1] = _medan(
            x[_idx2(i, 1, pcols)],
            xxk[_idx2(i, xin_k, pcols)],
            x[_idx2(i, pverp, pcols)],
        )
        intz[i - 1] = 0

    for k in range(1, pverp):
        for i in range(1, ncol + 1):
            if (
                (xins[i - 1] - x[_idx2(i, k, pcols)])
                * (x[_idx2(i, k + 1, pcols)] - xins[i - 1])
            ) >= 0.0:
                intz[i - 1] = k

    for i in range(1, ncol + 1):
        if intz[i - 1] == 0:
            status[0] = 1
            fail_i[0] = i
            fail_k[0] = xin_k
            return

    for i in range(1, ncol + 1):
        k = int(intz[i - 1])
        dx = x[_idx2(i, k + 1, pcols)] - x[_idx2(i, k, pcols)]
        s = (f[_idx2(i, k + 1, pcols)] - f[_idx2(i, k, pcols)]) / dx
        c2 = (3.0 * s - 2.0 * fdot[_idx2(i, k, pcols)] - fdot[_idx2(i, k + 1, pcols)]) / dx
        c3 = (
            fdot[_idx2(i, k, pcols)] + fdot[_idx2(i, k + 1, pcols)] - 2.0 * s
        ) / (dx * dx)
        xx = xins[i - 1] - x[_idx2(i, k, pcols)]
        fxdot[i - 1] = (3.0 * c3 * xx + 2.0 * c2) * xx + fdot[_idx2(i, k, pcols)]
        fxdd[i - 1] = 6.0 * c3 * xx + 2.0 * c2
        cfint = ((c3 * xx + c2) * xx + fdot[_idx2(i, k, pcols)]) * xx + f[_idx2(i, k, pcols)]

        psi1 = f[_idx2(i, k, pcols)] + (
            (f[_idx2(i, k + 1, pcols)] - f[_idx2(i, k, pcols)]) * xx / dx
        )
        if k == 1:
            psi2 = f[_idx2(i, 1, pcols)]
        else:
            psi2 = f[_idx2(i, k, pcols)] + (
                (f[_idx2(i, k, pcols)] - f[_idx2(i, k - 1, pcols)])
                * xx
                / (x[_idx2(i, k, pcols)] - x[_idx2(i, k - 1, pcols)])
            )

        if (k + 1) == pverp:
            psi3 = f[_idx2(i, pverp, pcols)]
        else:
            psi3 = f[_idx2(i, k + 1, pcols)] - (
                (f[_idx2(i, k + 2, pcols)] - f[_idx2(i, k + 1, pcols)])
                * (dx - xx)
                / (x[_idx2(i, k + 2, pcols)] - x[_idx2(i, k + 1, pcols)])
            )

        psim = _medan(psi1, psi2, psi3)
        cfnew = _medan(cfint, psi1, psim)
        psistar[i - 1] = cfnew


def _dust_cfdotmc_pro(
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    x_p: cobj,
    f_p: cobj,
    fdot_p: cobj,
    s_p: cobj,
    sh_p: cobj,
    d_p: cobj,
    dh_p: cobj,
    e_p: cobj,
    eh_p: cobj,
    ppl_p: cobj,
    ppr_p: cobj,
    delxh_p: cobj,
):
    x = Ptr[float](x_p)
    f = Ptr[float](f_p)
    fdot = Ptr[float](fdot_p)
    s = Ptr[float](s_p)
    sh = Ptr[float](sh_p)
    d = Ptr[float](d_p)
    dh = Ptr[float](dh_p)
    e = Ptr[float](e_p)
    eh = Ptr[float](eh_p)
    ppl = Ptr[float](ppl_p)
    ppr = Ptr[float](ppr_p)
    delxh = Ptr[float](delxh_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            delxh[_idx2(i, k, pcols)] = x[_idx2(i, k + 1, pcols)] - x[_idx2(i, k, pcols)]
            sh[_idx2(i, k, pcols)] = (
                f[_idx2(i, k + 1, pcols)] - f[_idx2(i, k, pcols)]
            ) / delxh[_idx2(i, k, pcols)]

        if k >= 2:
            for i in range(1, ncol + 1):
                d[_idx2(i, k, pcols)] = (
                    sh[_idx2(i, k, pcols)] - sh[_idx2(i, k - 1, pcols)]
                ) / (x[_idx2(i, k + 1, pcols)] - x[_idx2(i, k - 1, pcols)])
                s[_idx2(i, k, pcols)] = _minmod(
                    sh[_idx2(i, k, pcols)], sh[_idx2(i, k - 1, pcols)]
                )

    for k in range(2, pver):
        for i in range(1, ncol + 1):
            eh[_idx2(i, k, pcols)] = (
                d[_idx2(i, k + 1, pcols)] - d[_idx2(i, k, pcols)]
            ) / (x[_idx2(i, k + 2, pcols)] - x[_idx2(i, k - 1, pcols)])
            dh[_idx2(i, k, pcols)] = _minmod(
                d[_idx2(i, k, pcols)], d[_idx2(i, k + 1, pcols)]
            )

    for i in range(1, ncol + 1):
        e[_idx2(i, 2, pcols)] = eh[_idx2(i, 2, pcols)]
        e[_idx2(i, pver, pcols)] = eh[_idx2(i, pver - 1, pcols)]

        fdot[_idx2(i, 1, pcols)] = (
            sh[_idx2(i, 1, pcols)]
            - d[_idx2(i, 2, pcols)] * delxh[_idx2(i, 1, pcols)]
            - eh[_idx2(i, 2, pcols)]
            * delxh[_idx2(i, 1, pcols)]
            * (x[_idx2(i, 1, pcols)] - x[_idx2(i, 3, pcols)])
        )
        fdot[_idx2(i, 1, pcols)] = _minmod(
            fdot[_idx2(i, 1, pcols)], 3.0 * sh[_idx2(i, 1, pcols)]
        )

        fdot[_idx2(i, pverp, pcols)] = (
            sh[_idx2(i, pver, pcols)]
            + d[_idx2(i, pver, pcols)] * delxh[_idx2(i, pver, pcols)]
            + eh[_idx2(i, pver - 1, pcols)]
            * delxh[_idx2(i, pver, pcols)]
            * (x[_idx2(i, pverp, pcols)] - x[_idx2(i, pver - 1, pcols)])
        )
        fdot[_idx2(i, pverp, pcols)] = _minmod(
            fdot[_idx2(i, pverp, pcols)], 3.0 * sh[_idx2(i, pver, pcols)]
        )

        fdot[_idx2(i, 2, pcols)] = (
            sh[_idx2(i, 1, pcols)]
            + d[_idx2(i, 2, pcols)] * delxh[_idx2(i, 1, pcols)]
            - eh[_idx2(i, 2, pcols)]
            * delxh[_idx2(i, 1, pcols)]
            * delxh[_idx2(i, 2, pcols)]
        )
        fdot[_idx2(i, 2, pcols)] = _minmod(
            fdot[_idx2(i, 2, pcols)], 3.0 * s[_idx2(i, 2, pcols)]
        )

        fdot[_idx2(i, pver, pcols)] = (
            sh[_idx2(i, pver, pcols)]
            - d[_idx2(i, pver, pcols)] * delxh[_idx2(i, pver, pcols)]
            - eh[_idx2(i, pver - 1, pcols)]
            * delxh[_idx2(i, pver, pcols)]
            * delxh[_idx2(i, pver - 1, pcols)]
        )
        fdot[_idx2(i, pver, pcols)] = _minmod(
            fdot[_idx2(i, pver, pcols)], 3.0 * s[_idx2(i, pver, pcols)]
        )

    for k in range(3, pver):
        for i in range(1, ncol + 1):
            e[_idx2(i, k, pcols)] = _minmod(
                eh[_idx2(i, k, pcols)], eh[_idx2(i, k - 1, pcols)]
            )

    for k in range(3, pver):
        for i in range(1, ncol + 1):
            ppl[_idx2(i, k, pcols)] = (
                sh[_idx2(i, k - 1, pcols)] + dh[_idx2(i, k - 1, pcols)] * delxh[_idx2(i, k - 1, pcols)]
            )
            ppr[_idx2(i, k, pcols)] = (
                sh[_idx2(i, k, pcols)] - dh[_idx2(i, k, pcols)] * delxh[_idx2(i, k, pcols)]
            )

            t = _minmod(ppl[_idx2(i, k, pcols)], ppr[_idx2(i, k, pcols)])

            pp = sh[_idx2(i, k - 1, pcols)] + d[_idx2(i, k, pcols)] * delxh[_idx2(i, k - 1, pcols)]

            fdot[_idx2(i, k, pcols)] = pp - (
                delxh[_idx2(i, k - 1, pcols)]
                * delxh[_idx2(i, k, pcols)]
                * (
                    eh[_idx2(i, k - 1, pcols)] * (x[_idx2(i, k + 2, pcols)] - x[_idx2(i, k, pcols)])
                    + eh[_idx2(i, k, pcols)] * (x[_idx2(i, k, pcols)] - x[_idx2(i, k - 2, pcols)])
                )
                / (x[_idx2(i, k + 2, pcols)] - x[_idx2(i, k - 2, pcols)])
            )

            qpl = sh[_idx2(i, k - 1, pcols)] + delxh[_idx2(i, k - 1, pcols)] * _minmod(
                d[_idx2(i, k - 1, pcols)]
                + e[_idx2(i, k - 1, pcols)] * (x[_idx2(i, k, pcols)] - x[_idx2(i, k - 2, pcols)]),
                d[_idx2(i, k, pcols)] - e[_idx2(i, k, pcols)] * delxh[_idx2(i, k, pcols)],
            )
            qpr = sh[_idx2(i, k, pcols)] + delxh[_idx2(i, k, pcols)] * _minmod(
                d[_idx2(i, k, pcols)] + e[_idx2(i, k, pcols)] * delxh[_idx2(i, k - 1, pcols)],
                d[_idx2(i, k + 1, pcols)]
                + e[_idx2(i, k + 1, pcols)] * (x[_idx2(i, k, pcols)] - x[_idx2(i, k + 2, pcols)]),
            )

            fdot[_idx2(i, k, pcols)] = _medan(fdot[_idx2(i, k, pcols)], qpl, qpr)

            ttt = _minmod(qpl, qpr)
            tmin = _min4(
                0.0,
                3.0 * s[_idx2(i, k, pcols)],
                1.5 * t,
                ttt,
            )
            tmax = _max4(
                0.0,
                3.0 * s[_idx2(i, k, pcols)],
                1.5 * t,
                ttt,
            )

            fdot[_idx2(i, k, pcols)] = fdot[_idx2(i, k, pcols)] + _minmod(
                tmin - fdot[_idx2(i, k, pcols)],
                tmax - fdot[_idx2(i, k, pcols)],
            )


def _dust_getflx(
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    deltat: float,
    xw_p: cobj,
    phi_p: cobj,
    vel_p: cobj,
    flux_p: cobj,
    psi_p: cobj,
    fdot_p: cobj,
    xxk_p: cobj,
    fxdot_p: cobj,
    fxdd_p: cobj,
    psistar_p: cobj,
    xins_p: cobj,
    intz_p: cobj,
    status_p: cobj,
    fail_i_p: cobj,
    fail_k_p: cobj,
    s_p: cobj,
    sh_p: cobj,
    d_p: cobj,
    dh_p: cobj,
    e_p: cobj,
    eh_p: cobj,
    ppl_p: cobj,
    ppr_p: cobj,
    delxh_p: cobj,
):
    xw = Ptr[float](xw_p)
    phi = Ptr[float](phi_p)
    vel = Ptr[float](vel_p)
    flux = Ptr[float](flux_p)
    psi = Ptr[float](psi_p)
    xxk = Ptr[float](xxk_p)
    psistar = Ptr[float](psistar_p)
    status = Ptr[int](status_p)

    for i in range(1, ncol + 1):
        psi[_idx2(i, 1, pcols)] = 0.0
        flux[_idx2(i, 1, pcols)] = 0.0
        flux[_idx2(i, pverp, pcols)] = 0.0

    for k in range(2, pverp + 1):
        for i in range(1, ncol + 1):
            psi[_idx2(i, k, pcols)] = (
                phi[_idx2(i, k - 1, pcols)]
                * (xw[_idx2(i, k, pcols)] - xw[_idx2(i, k - 1, pcols)])
                + psi[_idx2(i, k - 1, pcols)]
            )

    _dust_cfdotmc_pro(
        ncol, pcols, pver, pverp, xw_p, psi_p, fdot_p, s_p, sh_p, d_p, dh_p, e_p, eh_p, ppl_p, ppr_p, delxh_p
    )

    for k in range(2, pver + 1):
        for i in range(1, ncol + 1):
            xxk[_idx2(i, k, pcols)] = xw[_idx2(i, k, pcols)] - vel[_idx2(i, k, pcols)] * deltat

    for k in range(2, pver + 1):
        _dust_cfint2(
            ncol,
            pcols,
            pverp,
            k,
            xw_p,
            psi_p,
            fdot_p,
            xxk_p,
            fxdot_p,
            fxdd_p,
            psistar_p,
            xins_p,
            intz_p,
            status_p,
            fail_i_p,
            fail_k_p,
        )
        if status[0] != 0:
            return
        for i in range(1, ncol + 1):
            flux[_idx2(i, k, pcols)] = psi[_idx2(i, k, pcols)] - psistar[i - 1]


@export
def dust_sediment_tend_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pverp: int,
    dtime: float,
    mxsedfac: float,
    gravit: float,
    pint_p: cobj,
    pdel_p: cobj,
    dustmr_p: cobj,
    pvdust_p: cobj,
    dusttend_p: cobj,
    sfdust_p: cobj,
    fxdust_p: cobj,
    psi_p: cobj,
    fdot_p: cobj,
    xxk_p: cobj,
    fxdot_p: cobj,
    fxdd_p: cobj,
    psistar_p: cobj,
    s_p: cobj,
    sh_p: cobj,
    d_p: cobj,
    dh_p: cobj,
    e_p: cobj,
    eh_p: cobj,
    ppl_p: cobj,
    ppr_p: cobj,
    delxh_p: cobj,
    xins_p: cobj,
    intz_p: cobj,
    status_p: cobj,
    fail_i_p: cobj,
    fail_k_p: cobj,
):
    pdel = Ptr[float](pdel_p)
    dustmr = Ptr[float](dustmr_p)
    pvdust = Ptr[float](pvdust_p)
    fxdust = Ptr[float](fxdust_p)
    dusttend = Ptr[float](dusttend_p)
    sfdust = Ptr[float](sfdust_p)
    status = Ptr[int](status_p)
    fail_i = Ptr[int](fail_i_p)
    fail_k = Ptr[int](fail_k_p)

    status[0] = 0
    fail_i[0] = 0
    fail_k[0] = 0

    for k in range(1, pverp + 1):
        for i in range(1, ncol + 1):
            fxdust[_idx2(i, k, pcols)] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            dusttend[_idx2(i, k, pcols)] = 0.0

    for i in range(1, ncol + 1):
        sfdust[i - 1] = 0.0

    _dust_getflx(
        ncol,
        pcols,
        pver,
        pverp,
        dtime,
        pint_p,
        dustmr_p,
        pvdust_p,
        fxdust_p,
        psi_p,
        fdot_p,
        xxk_p,
        fxdot_p,
        fxdd_p,
        psistar_p,
        xins_p,
        intz_p,
        status_p,
        fail_i_p,
        fail_k_p,
        s_p,
        sh_p,
        d_p,
        dh_p,
        e_p,
        eh_p,
        ppl_p,
        ppr_p,
        delxh_p,
    )
    if status[0] != 0:
        return

    for i in range(1, ncol + 1):
        fxdust[_idx2(i, 1, pcols)] = 0.0
        fxdust[_idx2(i, pverp, pcols)] = (
            dustmr[_idx2(i, pver, pcols)] * pvdust[_idx2(i, pverp, pcols)] * dtime
        )

    for k in range(2, pver + 1):
        for i in range(1, ncol + 1):
            fxdust[_idx2(i, k, pcols)] = max(0.0, fxdust[_idx2(i, k, pcols)])

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            fxdust[_idx2(i, k + 1, pcols)] = min(
                fxdust[_idx2(i, k + 1, pcols)],
                mxsedfac * dustmr[_idx2(i, k, pcols)] * pdel[_idx2(i, k, pcols)],
            )

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            dusttend[_idx2(i, k, pcols)] = (
                fxdust[_idx2(i, k, pcols)] - fxdust[_idx2(i, k + 1, pcols)]
            ) / (dtime * pdel[_idx2(i, k, pcols)])

    for i in range(1, ncol + 1):
        sfdust[i - 1] = fxdust[_idx2(i, pverp, pcols)] / (dtime * gravit)


@export
def qqcw2vmr_codon(
    ncol: int,
    pver: int,
    fldcw_ld1: int,
    mbar_p: cobj,
    fldcw_p: cobj,
    adv_mass: float,
    vmr_p: cobj,
):
    mbar = Ptr[float](mbar_p)
    fldcw = Ptr[float](fldcw_p)
    vmr = Ptr[float](vmr_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            vmr[_idx2(i, k, ncol)] = (
                mbar[_idx2(i, k, ncol)] * fldcw[_idx2(i, k, fldcw_ld1)] / adv_mass
            )


@export
def vmr2qqcw_codon(
    ncol: int,
    pver: int,
    fldcw_ld1: int,
    vmr_p: cobj,
    mbar_p: cobj,
    adv_mass: float,
    fldcw_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    mbar = Ptr[float](mbar_p)
    fldcw = Ptr[float](fldcw_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            fldcw[_idx2(i, k, fldcw_ld1)] = (
                adv_mass * vmr[_idx2(i, k, ncol)] / mbar[_idx2(i, k, ncol)]
            )


@export
def modal_aero_gasaerexch_zero_tendencies_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnstxx: int,
    nsrflx: int,
    dqdt_p: cobj,
    dqqcwdt_p: cobj,
    qsrflx_p: cobj,
    qqcwsrflx_p: cobj,
):
    dqdt = Ptr[float](dqdt_p)
    dqqcwdt = Ptr[float](dqqcwdt_p)
    qsrflx = Ptr[float](qsrflx_p)
    qqcwsrflx = Ptr[float](qqcwsrflx_p)

    for m in range(1, pcnstxx + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                dqdt[_idx3(i, k, m, ncol, pver)] = 0.0
                dqqcwdt[_idx3(i, k, m, ncol, pver)] = 0.0

    for jsrf in range(1, nsrflx + 1):
        for m in range(1, pcnstxx + 1):
            for i in range(1, pcols + 1):
                idx = _idx3(i, m, jsrf, pcols, pcnstxx)
                qsrflx[idx] = 0.0
                qqcwsrflx[idx] = 0.0


@export
def modal_aero_newnuc_zero_tendencies_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnstxx: int,
    pcnst: int,
    nsrflx: int,
    dqdt_p: cobj,
    qsrflx_p: cobj,
):
    dqdt = Ptr[float](dqdt_p)
    qsrflx = Ptr[float](qsrflx_p)

    for m in range(1, pcnstxx + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                dqdt[_idx3(i, k, m, ncol, pver)] = 0.0

    for jsrf in range(1, nsrflx + 1):
        for m in range(1, pcnst + 1):
            for i in range(1, ncol + 1):
                qsrflx[_idx3(i, m, jsrf, pcols, pcnst)] = 0.0


@export
def modal_aero_newnuc_prepare_box_inputs_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    l_h2so4: int,
    l_nh3: int,
    do_nh3: int,
    deltat: float,
    qh2so4_cutoff: float,
    q_p: cobj,
    qv_p: cobj,
    cld_p: cobj,
    qv_sat_p: cobj,
    del_h2so4_gasprod_p: cobj,
    del_h2so4_aeruptk_p: cobj,
    active_mask_p: cobj,
    cldx_p: cobj,
    qh2so4_cur_p: cobj,
    qh2so4_avg_p: cobj,
    qnh3_cur_p: cobj,
    tmp_uptkrate_p: cobj,
    relhumnn_p: cobj,
):
    q = Ptr[float](q_p)
    qv = Ptr[float](qv_p)
    cld = Ptr[float](cld_p)
    qv_sat = Ptr[float](qv_sat_p)
    del_h2so4_gasprod = Ptr[float](del_h2so4_gasprod_p)
    del_h2so4_aeruptk = Ptr[float](del_h2so4_aeruptk_p)
    active_mask = Ptr[int](active_mask_p)
    cldx_out = Ptr[float](cldx_p)
    qh2so4_cur_out = Ptr[float](qh2so4_cur_p)
    qh2so4_avg_out = Ptr[float](qh2so4_avg_p)
    qnh3_cur_out = Ptr[float](qnh3_cur_p)
    tmp_uptkrate_out = Ptr[float](tmp_uptkrate_p)
    relhumnn_out = Ptr[float](relhumnn_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            out_idx = _idx2(i, k, ncol)
            active_mask[out_idx] = 0
            cldx_out[out_idx] = 0.0
            qh2so4_cur_out[out_idx] = 0.0
            qh2so4_avg_out[out_idx] = 0.0
            qnh3_cur_out[out_idx] = 0.0
            tmp_uptkrate_out[out_idx] = 0.0
            relhumnn_out[out_idx] = 0.0

            if k < top_lev:
                continue

            cld_val = cld[_idx2(i, k, ncol)]
            if cld_val >= 0.99:
                continue

            qh2so4_cur = q[_idx3(i, k, l_h2so4, ncol, pver)]
            if qh2so4_cur <= qh2so4_cutoff:
                continue

            tmpa = del_h2so4_gasprod[_idx2(i, k, ncol)]
            if tmpa < 0.0:
                tmpa = 0.0
            tmp_q3 = qh2so4_cur

            tmpb = -del_h2so4_aeruptk[_idx2(i, k, ncol)]
            if tmpb < 0.0:
                tmpb = 0.0
            tmp_q2 = tmp_q3 + tmpb

            if tmp_q2 <= tmp_q3:
                tmpb = 0.0
            else:
                tmpc = tmp_q2 * exp(-20.0)
                if tmp_q3 <= tmpc:
                    tmp_q3 = tmpc
                    tmpb = 20.0
                else:
                    tmpb = log(tmp_q2 / tmp_q3)

            tmp_uptkrate = tmpb / deltat

            if tmpb <= 0.1:
                qh2so4_avg = tmp_q3 * (1.0 + 0.5 * tmpb) - 0.5 * tmpa
            else:
                tmpc = tmpa / tmpb
                qh2so4_avg = (tmp_q3 - tmpc) * ((exp(tmpb) - 1.0) / tmpb) + tmpc
            if qh2so4_avg <= qh2so4_cutoff:
                continue

            if do_nh3 != 0:
                qnh3_cur = q[_idx3(i, k, l_nh3, ncol, pver)]
                if qnh3_cur < 0.0:
                    qnh3_cur = 0.0
            else:
                qnh3_cur = 0.0

            qvswtr = qv_sat[_idx2(i, k, pcols)]
            if qvswtr < 1.0e-20:
                qvswtr = 1.0e-20
            relhumav = qv[_idx2(i, k, pcols)] / qvswtr
            if relhumav < 0.0:
                relhumav = 0.0
            elif relhumav > 1.0:
                relhumav = 1.0

            cldx = cld_val
            if cldx < 0.0:
                cldx = 0.0
            relhum = (relhumav - cldx) / (1.0 - cldx)
            if relhum < 0.0:
                relhum = 0.0
            elif relhum > 1.0:
                relhum = 1.0

            relhumnn = relhum
            if relhumnn < 0.01:
                relhumnn = 0.01
            elif relhumnn > 0.99:
                relhumnn = 0.99

            active_mask[out_idx] = 1
            cldx_out[out_idx] = cldx
            qh2so4_cur_out[out_idx] = qh2so4_cur
            qh2so4_avg_out[out_idx] = qh2so4_avg
            qnh3_cur_out[out_idx] = qnh3_cur
            tmp_uptkrate_out[out_idx] = tmp_uptkrate
            relhumnn_out[out_idx] = relhumnn


@export
def pbl_nuc_wang2008_codon(
    so4vol: float,
    newnuc_method_flagaa: int,
    newnuc_method_flagaa2_p: cobj,
    ratenucl_p: cobj,
    rateloge_p: cobj,
    cnum_tot_p: cobj,
    cnum_h2so4_p: cobj,
    cnum_nh3_p: cobj,
    radius_cluster_p: cobj,
):
    newnuc_method_flagaa2 = Ptr[int](newnuc_method_flagaa2_p)
    ratenucl = Ptr[float](ratenucl_p)
    rateloge = Ptr[float](rateloge_p)
    cnum_tot = Ptr[float](cnum_tot_p)
    cnum_h2so4 = Ptr[float](cnum_h2so4_p)
    cnum_nh3 = Ptr[float](cnum_nh3_p)
    radius_cluster = Ptr[float](radius_cluster_p)

    if newnuc_method_flagaa == 11:
        tmp_ratenucl = 1.0e-6 * so4vol
    elif newnuc_method_flagaa == 12:
        tmp_ratenucl = 1.0e-12 * (so4vol * so4vol)
    else:
        return

    tmp_rateloge = log(tmp_ratenucl)
    if tmp_rateloge <= rateloge[0]:
        return

    rateloge[0] = tmp_rateloge
    ratenucl[0] = tmp_ratenucl
    newnuc_method_flagaa2[0] = newnuc_method_flagaa

    radius_cluster[0] = 0.5
    tmp_diam = radius_cluster[0] * 2.0e-7
    tmp_volu = (tmp_diam * tmp_diam * tmp_diam) * (pi / 6.0)
    tmp_mass = tmp_volu * 1.8
    cnum_h2so4[0] = (tmp_mass / 98.0) * 6.023e23
    cnum_tot[0] = cnum_h2so4[0]
    cnum_nh3[0] = 0.0


@export
def binary_nuc_vehk2002_codon(
    temp: float,
    rh: float,
    so4vol: float,
    ratenucl_p: cobj,
    rateloge_p: cobj,
    cnum_h2so4_p: cobj,
    cnum_tot_p: cobj,
    radius_cluster_p: cobj,
):
    ratenucl = Ptr[float](ratenucl_p)
    rateloge = Ptr[float](rateloge_p)
    cnum_h2so4 = Ptr[float](cnum_h2so4_p)
    cnum_tot = Ptr[float](cnum_tot_p)
    radius_cluster = Ptr[float](radius_cluster_p)

    crit_x = (
        0.740997
        - 0.00266379 * temp
        - 0.00349998 * log(so4vol)
        + 0.0000504022 * temp * log(so4vol)
        + 0.00201048 * log(rh)
        - 0.000183289 * temp * log(rh)
        + 0.00157407 * (log(rh)) ** 2.0
        - 0.0000179059 * temp * (log(rh)) ** 2.0
        + 0.000184403 * (log(rh)) ** 3.0
        - 1.50345e-6 * temp * (log(rh)) ** 3.0
    )

    acoe = 0.14309 + 2.21956 * temp - 0.0273911 * temp**2.0 + 0.0000722811 * temp**3.0 + 5.91822 / crit_x
    bcoe = 0.117489 + 0.462532 * temp - 0.0118059 * temp**2.0 + 0.0000404196 * temp**3.0 + 15.7963 / crit_x
    ccoe = -0.215554 - 0.0810269 * temp + 0.00143581 * temp**2.0 - 4.7758e-6 * temp**3.0 - 2.91297 / crit_x
    dcoe = -3.58856 + 0.049508 * temp - 0.00021382 * temp**2.0 + 3.10801e-7 * temp**3.0 - 0.0293333 / crit_x
    ecoe = 1.14598 - 0.600796 * temp + 0.00864245 * temp**2.0 - 0.0000228947 * temp**3.0 - 8.44985 / crit_x
    fcoe = 2.15855 + 0.0808121 * temp - 0.000407382 * temp**2.0 - 4.01957e-7 * temp**3.0 + 0.721326 / crit_x
    gcoe = 1.6241 - 0.0160106 * temp + 0.0000377124 * temp**2.0 + 3.21794e-8 * temp**3.0 - 0.0113255 / crit_x
    hcoe = 9.71682 - 0.115048 * temp + 0.000157098 * temp**2.0 + 4.00914e-7 * temp**3.0 + 0.71186 / crit_x
    icoe = -1.05611 + 0.00903378 * temp - 0.0000198417 * temp**2.0 + 2.46048e-8 * temp**3.0 - 0.0579087 / crit_x
    jcoe = -0.148712 + 0.00283508 * temp - 9.24619e-6 * temp**2.0 + 5.00427e-9 * temp**3.0 - 0.0127081 / crit_x

    tmpa = (
        acoe
        + bcoe * log(rh)
        + ccoe * (log(rh)) ** 2.0
        + dcoe * (log(rh)) ** 3.0
        + ecoe * log(so4vol)
        + fcoe * (log(rh)) * (log(so4vol))
        + gcoe * ((log(rh)) ** 2.0) * (log(so4vol))
        + hcoe * (log(so4vol)) ** 2.0
        + icoe * log(rh) * ((log(so4vol)) ** 2.0)
        + jcoe * (log(so4vol)) ** 3.0
    )
    rateloge[0] = tmpa
    tmpa = min(tmpa, log(1.0e38))
    ratenucl[0] = exp(tmpa)

    acoe = -0.00295413 - 0.0976834 * temp + 0.00102485 * temp**2.0 - 2.18646e-6 * temp**3.0 - 0.101717 / crit_x
    bcoe = -0.00205064 - 0.00758504 * temp + 0.000192654 * temp**2.0 - 6.7043e-7 * temp**3.0 - 0.255774 / crit_x
    ccoe = 0.00322308 + 0.000852637 * temp - 0.0000154757 * temp**2.0 + 5.66661e-8 * temp**3.0 + 0.0338444 / crit_x
    dcoe = 0.0474323 - 0.000625104 * temp + 2.65066e-6 * temp**2.0 - 3.67471e-9 * temp**3.0 - 0.000267251 / crit_x
    ecoe = -0.0125211 + 0.00580655 * temp - 0.000101674 * temp**2.0 + 2.88195e-7 * temp**3.0 + 0.0942243 / crit_x
    fcoe = -0.038546 - 0.000672316 * temp + 2.60288e-6 * temp**2.0 + 1.19416e-8 * temp**3.0 - 0.00851515 / crit_x
    gcoe = -0.0183749 + 0.000172072 * temp - 3.71766e-7 * temp**2.0 - 5.14875e-10 * temp**3.0 + 0.00026866 / crit_x
    hcoe = -0.0619974 + 0.000906958 * temp - 9.11728e-7 * temp**2.0 - 5.36796e-9 * temp**3.0 - 0.00774234 / crit_x
    icoe = 0.0121827 - 0.00010665 * temp + 2.5346e-7 * temp**2.0 - 3.63519e-10 * temp**3.0 + 0.000610065 / crit_x
    jcoe = 0.000320184 - 0.0000174762 * temp + 6.06504e-8 * temp**2.0 - 1.4177e-11 * temp**3.0 + 0.000135751 / crit_x

    cnum_tot[0] = exp(
        acoe
        + bcoe * log(rh)
        + ccoe * (log(rh)) ** 2.0
        + dcoe * (log(rh)) ** 3.0
        + ecoe * log(so4vol)
        + fcoe * (log(rh)) * (log(so4vol))
        + gcoe * ((log(rh)) ** 2.0) * (log(so4vol))
        + hcoe * (log(so4vol)) ** 2.0
        + icoe * log(rh) * ((log(so4vol)) ** 2.0)
        + jcoe * (log(so4vol)) ** 3.0
    )

    cnum_h2so4[0] = cnum_tot[0] * crit_x
    radius_cluster[0] = exp(-1.6524245 + 0.42316402 * crit_x + 0.3346648 * log(cnum_tot[0]))


def _ternary_nuc_merik2007_core(
    t: float,
    rh: float,
    c2: float,
    c3: float,
):
    j_log = 0.0
    ntot = 0.0
    nacid = 0.0
    namm = 0.0
    r = 0.0
    log_c2 = log(c2)
    log_c3 = log(c3)
    log_rh = log(rh)
    t_sq = t**2
    t_cu = t**3
    c3_cu = c3**3
    log_c2_sq = log_c2**2
    log_c3_sq = log_c3**2
    log_c3_cu = log_c3**3

    t_onset = (
        143.6002929064716
        + 1.0178856665693992 * rh
        + 10.196398812974294 * log_c2
        - 0.1849879416839113 * log_c2_sq
        - 17.161783213150173 * log_c3
        + (109.92469248546053 * log_c3) / log_c2
        + 0.7734119613144357 * log_c2 * log_c3
        - 0.15576469879527022 * log_c3_sq
    )

    if t_onset > t:
        j_log = (
            -12.861848898625231
            + 4.905527742256349 * c3
            - 358.2337705052991 * rh
            - 0.05463019231872484 * c3 * t
            + 4.8630382337426985 * rh * t
            + 0.00020258394697064567 * c3 * t_sq
            - 0.02175548069741675 * rh * t_sq
            - 2.502406532869512e-7 * c3 * t_cu
            + 0.00003212869941055865 * rh * t_cu
            - 4.39129415725234e6 / log_c2_sq
            + (56383.93843154586 * t) / log_c2_sq
            - (239.835990963361 * t_sq) / log_c2_sq
            + (0.33765136625580167 * t_cu) / log_c2_sq
            - (629.7882041830943 * rh) / (c3_cu * log_c2)
            + (7.772806552631709 * rh * t) / (c3_cu * log_c2)
            - (0.031974053936299256 * rh * t_sq) / (c3_cu * log_c2)
            + (0.00004383764128775082 * rh * t_cu) / (c3_cu * log_c2)
            + 1200.472096232311 * log_c2
            - 17.37107890065621 * t * log_c2
            + 0.08170681335921742 * t_sq * log_c2
            - 0.00012534476159729881 * t_cu * log_c2
            - 14.833042158178936 * log_c2_sq
            + 0.2932631303555295 * t * log_c2_sq
            - 0.0016497524241142845 * t_sq * log_c2_sq
            + 2.844074805239367e-6 * t_cu * log_c2_sq
            - 231375.56676032578 * log_c3
            - 100.21645273730675 * rh * log_c3
            + 2919.2852552424706 * t * log_c3
            + 0.977886555834732 * rh * t * log_c3
            - 12.286497122264588 * t_sq * log_c3
            - 0.0030511783284506377 * rh * t_sq * log_c3
            + 0.017249301826661612 * t_cu * log_c3
            + 2.967320346100855e-6 * rh * t_cu * log_c3
            + (2.360931724951942e6 * log_c3) / log_c2
            - (29752.130254319443 * t * log_c3) / log_c2
            + (125.04965118142027 * t_sq * log_c3) / log_c2
            - (0.1752996881934318 * t_cu * log_c3) / log_c2
            + 5599.912337254629 * log_c2 * log_c3
            - 70.70896612937771 * t * log_c2 * log_c3
            + 0.2978801613269466 * t_sq * log_c2 * log_c3
            - 0.00041866525019504 * t_cu * log_c2 * log_c3
            + 75061.15281456841 * log_c3_sq
            - 931.8802278173565 * t * log_c3_sq
            + 3.863266220840964 * t_sq * log_c3_sq
            - 0.005349472062284983 * t_cu * log_c3_sq
            - (732006.8180571689 * log_c3_sq) / log_c2
            + (9100.06398573816 * t * log_c3_sq) / log_c2
            - (37.771091915932004 * t_sq * log_c3_sq) / log_c2
            + (0.05235455395566905 * t_cu * log_c3_sq) / log_c2
            - 1911.0303773001353 * log_c2 * log_c3_sq
            + 23.6903969622286 * t * log_c2 * log_c3_sq
            - 0.09807872005428583 * t_sq * log_c2 * log_c3_sq
            + 0.00013564560238552576 * t_cu * log_c2 * log_c3_sq
            - 3180.5610833308 * log_c3_cu
            + 39.08268568672095 * t * log_c3_cu
            - 0.16048521066690752 * t_sq * log_c3_cu
            + 0.00022031380023793877 * t_cu * log_c3_cu
            + (40751.075322248245 * log_c3_cu) / log_c2
            - (501.66977622013934 * t * log_c3_cu) / log_c2
            + (2.063469732254135 * t_sq * log_c3_cu) / log_c2
            - (0.002836873785758324 * t_cu * log_c3_cu) / log_c2
            + 2.792313345723013 * log_c2_sq * log_c3_cu
            - 0.03422552111802899 * t * log_c2_sq * log_c3_cu
            + 0.00014019195277521142 * t_sq * log_c2_sq * log_c3_cu
            - 1.9201227328396297e-7 * t_cu * log_c2_sq * log_c3_cu
            - 980.923146020468 * log_rh
            + 10.054155220444462 * t * log_rh
            - 0.03306644502023841 * t_sq * log_rh
            + 0.000034274041225891804 * t_cu * log_rh
            + (16597.75554295064 * log_rh) / log_c2
            - (175.2365504237746 * t * log_rh) / log_c2
            + (0.6033215603167458 * t_sq * log_rh) / log_c2
            - (0.0006731787599587544 * t_cu * log_rh) / log_c2
            - 89.38961120336789 * log_c3 * log_rh
            + 1.153344219304926 * t * log_c3 * log_rh
            - 0.004954549700267233 * t_sq * log_c3 * log_rh
            + 7.096309866238719e-6 * t_cu * log_c3 * log_rh
            + 3.1712136610383244 * log_c3_cu * log_rh
            - 0.037822330602328806 * t * log_c3_cu * log_rh
            + 0.0001500555743561457 * t_sq * log_c3_cu * log_rh
            - 1.9828365865570703e-7 * t_cu * log_c3_cu * log_rh
        )

        j = exp(j_log)
        log_j = log(j)
        log_j_sq = log_j**2

        ntot = (
            57.40091052369212
            - 0.2996341884645408 * t
            + 0.0007395477768531926 * t_sq
            - 5.090604835032423 * log_c2
            + 0.011016634044531128 * t * log_c2
            + 0.06750032251225707 * log_c2_sq
            - 0.8102831333223962 * log_c3
            + 0.015905081275952426 * t * log_c3
            - 0.2044174683159531 * log_c2 * log_c3
            + 0.08918159167625832 * log_c3_sq
            - 0.0004969033586666147 * t * log_c3_sq
            + 0.005704394549007816 * log_c3_cu
            + 3.4098703903474368 * log_j
            - 0.014916956508210809 * t * log_j
            + 0.08459090011666293 * log_c3 * log_j
            - 0.00014800625143907616 * t * log_c3 * log_j
            + 0.00503804694656905 * log_j_sq
        )

        r = (
            3.2888553966535506e-10
            - 3.374171768439839e-12 * t
            + 1.8347359507774313e-14 * t_sq
            + 2.5419844298881856e-12 * log_c2
            - 9.498107643050827e-14 * t * log_c2
            + 7.446266520834559e-13 * log_c2_sq
            + 2.4303397746137294e-11 * log_c3
            + 1.589324325956633e-14 * t * log_c3
            - 2.034596219775266e-12 * log_c2 * log_c3
            - 5.59303954457172e-13 * log_c3_sq
            - 4.889507104645867e-16 * t * log_c3_sq
            + 1.3847024107506764e-13 * log_c3_cu
            + 4.141077193427042e-15 * log_j
            - 2.6813110884009767e-14 * t * log_j
            + 1.2879071621313094e-12 * log_c3 * log_j
            - 3.80352446061867e-15 * t * log_c3 * log_j
            - 1.8790172502456827e-14 * log_j_sq
        )

        nacid = (
            -4.7154180661803595
            + 0.13436423483953885 * t
            - 0.00047184686478816176 * t_sq
            - 2.564010713640308 * log_c2
            + 0.011353312899114723 * t * log_c2
            + 0.0010801941974317014 * log_c2_sq
            + 0.5171368624197119 * log_c3
            - 0.0027882479896204665 * t * log_c3
            + 0.8066971907026886 * log_c3_sq
            - 0.0031849094214409335 * t * log_c3_sq
            - 0.09951184152927882 * log_c3_cu
            + 0.00040072788891745513 * t * log_c3_cu
            + 1.3276469271073974 * log_j
            - 0.006167654171986281 * t * log_j
            - 0.11061390967822708 * log_c3 * log_j
            + 0.0004367575329273496 * t * log_c3 * log_j
            + 0.000916366357266258 * log_j_sq
        )

        namm = (
            71.20073903979772
            - 0.8409600103431923 * t
            + 0.0024803006590334922 * t_sq
            + 2.7798606841602607 * log_c2
            - 0.01475023348171676 * t * log_c2
            + 0.012264508212031405 * log_c2_sq
            - 2.009926050440182 * log_c3
            + 0.008689123511431527 * t * log_c3
            - 0.009141180198955415 * log_c2 * log_c3
            + 0.1374122553905617 * log_c3_sq
            - 0.0006253227821679215 * t * log_c3_sq
            + 0.00009377332742098946 * log_c3_cu
            + 0.5202974341687757 * log_j
            - 0.002419872323052805 * t * log_j
            + 0.07916392322884074 * log_c3 * log_j
            - 0.0003021586030317366 * t * log_c3 * log_j
            + 0.0046977006608603395 * log_j_sq
        )
    else:
        j_log = -300.0

    return (j_log, ntot, nacid, namm, r)


@export
def ternary_nuc_merik2007_codon(
    t: float,
    rh: float,
    c2: float,
    c3: float,
    j_log_p: cobj,
    ntot_p: cobj,
    nacid_p: cobj,
    namm_p: cobj,
    r_p: cobj,
):
    j_log = Ptr[float](j_log_p)
    ntot = Ptr[float](ntot_p)
    nacid = Ptr[float](nacid_p)
    namm = Ptr[float](namm_p)
    r = Ptr[float](r_p)

    (
        j_log[0],
        ntot[0],
        nacid[0],
        namm[0],
        r[0],
    ) = _ternary_nuc_merik2007_core(t, rh, c2, c3)


@export
def mer07_veh02_nuc_mosaic_init_state_codon(
    newnuc_method_flagaa: int,
    isize_nuc_p: cobj,
    qnuma_del_p: cobj,
    qso4a_del_p: cobj,
    qnh4a_del_p: cobj,
    qh2so4_del_p: cobj,
    qnh3_del_p: cobj,
    valid_method_p: cobj,
):
    isize_nuc = Ptr[int](isize_nuc_p)
    qnuma_del = Ptr[float](qnuma_del_p)
    qso4a_del = Ptr[float](qso4a_del_p)
    qnh4a_del = Ptr[float](qnh4a_del_p)
    qh2so4_del = Ptr[float](qh2so4_del_p)
    qnh3_del = Ptr[float](qnh3_del_p)
    valid_method = Ptr[int](valid_method_p)

    isize_nuc[0] = 1
    qnuma_del[0] = 0.0
    qso4a_del[0] = 0.0
    qnh4a_del[0] = 0.0
    qh2so4_del[0] = 0.0
    qnh3_del[0] = 0.0
    valid_method[0] = 0

    if (
        newnuc_method_flagaa == 1
        or newnuc_method_flagaa == 2
        or newnuc_method_flagaa == 11
        or newnuc_method_flagaa == 12
    ):
        valid_method[0] = 1


@export
def mer07_veh02_nuc_mosaic_prepare_finalize_inputs_codon(
    rateloge: float,
    ratenuclt_p: cobj,
    ratenuclt_bb_p: cobj,
    continue_flag_p: cobj,
):
    ratenuclt = Ptr[float](ratenuclt_p)
    ratenuclt_bb = Ptr[float](ratenuclt_bb_p)
    continue_flag = Ptr[int](continue_flag_p)

    ratenuclt[0] = 0.0
    ratenuclt_bb[0] = 0.0
    continue_flag[0] = 0

    if rateloge <= -13.82:
        return

    ratenuclt[0] = exp(rateloge)
    ratenuclt_bb[0] = ratenuclt[0] * 1.0e6
    continue_flag[0] = 1


@export
def mer07_veh02_nuc_mosaic_prepare_rates_codon(
    newnuc_method_flagaa: int,
    temp_in: float,
    rh_in: float,
    press_in: float,
    zm_in: float,
    pblh_in: float,
    qh2so4_avg: float,
    qnh3_cur: float,
    rgas: float,
    avogad: float,
    cair_p: cobj,
    so4vol_in_p: cobj,
    nh3ppt_p: cobj,
    ratenuclt_p: cobj,
    rateloge_p: cobj,
    temp_bb_p: cobj,
    rh_bb_p: cobj,
    so4vol_bb_p: cobj,
    nh3ppt_bb_p: cobj,
    newnuc_method_flagaa2_p: cobj,
    use_ternary_rate_p: cobj,
    use_binary_rate_p: cobj,
    do_pbl_rate_p: cobj,
):
    cair = Ptr[float](cair_p)
    so4vol_in = Ptr[float](so4vol_in_p)
    nh3ppt = Ptr[float](nh3ppt_p)
    ratenuclt = Ptr[float](ratenuclt_p)
    rateloge = Ptr[float](rateloge_p)
    temp_bb = Ptr[float](temp_bb_p)
    rh_bb = Ptr[float](rh_bb_p)
    so4vol_bb = Ptr[float](so4vol_bb_p)
    nh3ppt_bb = Ptr[float](nh3ppt_bb_p)
    newnuc_method_flagaa2 = Ptr[int](newnuc_method_flagaa2_p)
    use_ternary_rate = Ptr[int](use_ternary_rate_p)
    use_binary_rate = Ptr[int](use_binary_rate_p)
    do_pbl_rate = Ptr[int](do_pbl_rate_p)

    cair[0] = press_in / (temp_in * rgas)
    so4vol_in[0] = qh2so4_avg * cair[0] * avogad * 1.0e-6
    nh3ppt[0] = qnh3_cur * 1.0e12
    ratenuclt[0] = 1.0e-38
    rateloge[0] = log(ratenuclt[0])
    temp_bb[0] = 0.0
    rh_bb[0] = 0.0
    so4vol_bb[0] = 0.0
    nh3ppt_bb[0] = 0.0
    use_ternary_rate[0] = 0
    use_binary_rate[0] = 0
    do_pbl_rate[0] = 0

    if (newnuc_method_flagaa != 2) and (nh3ppt[0] >= 0.1):
        if so4vol_in[0] >= 5.0e4:
            temp_bb[0] = max(235.0, min(295.0, temp_in))
            rh_bb[0] = max(0.05, min(0.95, rh_in))
            so4vol_bb[0] = max(5.0e4, min(1.0e9, so4vol_in[0]))
            nh3ppt_bb[0] = max(0.1, min(1.0e3, nh3ppt[0]))
            use_ternary_rate[0] = 1
        newnuc_method_flagaa2[0] = 1
    else:
        if so4vol_in[0] >= 1.0e4:
            temp_bb[0] = max(230.15, min(305.15, temp_in))
            rh_bb[0] = max(1.0e-4, min(1.0, rh_in))
            so4vol_bb[0] = max(1.0e4, min(1.0e11, so4vol_in[0]))
            use_binary_rate[0] = 1
        newnuc_method_flagaa2[0] = 2

    if (newnuc_method_flagaa == 11) or (newnuc_method_flagaa == 12):
        if zm_in <= max(pblh_in, 100.0):
            do_pbl_rate[0] = 1


def _mer07_veh02_nuc_mosaic_postprocess_core(
    qnuma_del: float,
    qso4a_del: float,
    qnh4a_del: float,
    deltat: float,
    specmw_so4_amode: float,
    specmw_nh4_amode: float,
    mass1p_aitlo: float,
    mass1p_aithi: float,
):
    qnuma_del = qnuma_del * 1.0e3
    dndt_ait = qnuma_del / deltat
    tmpa = qso4a_del * specmw_so4_amode
    tmpb = tmpa + qnh4a_del * specmw_nh4_amode
    tmp_frso4 = max(tmpa, 1.0e-35) / max(tmpb, 1.0e-35)
    dmdt_ait = max(0.0, tmpb / deltat)

    dndt_aitsv1 = dndt_ait
    dmdt_aitsv1 = dmdt_ait
    dndt_aitsv2 = 0.0
    dmdt_aitsv2 = 0.0
    dndt_aitsv3 = 0.0
    dmdt_aitsv3 = 0.0
    postprocess_code = 0

    if dndt_ait < 1.0e2:
        dndt_ait = 0.0
        dmdt_ait = 0.0
        postprocess_code = 1
    else:
        dndt_aitsv2 = dndt_ait
        dmdt_aitsv2 = dmdt_ait
        postprocess_code = 2
        mass1p = dmdt_ait / dndt_ait
        dndt_aitsv3 = dndt_ait
        dmdt_aitsv3 = dmdt_ait

        if mass1p < mass1p_aitlo:
            dndt_ait = dmdt_ait / mass1p_aitlo
            postprocess_code = 3
        elif mass1p > mass1p_aithi:
            dmdt_ait = dndt_ait * mass1p_aithi
            postprocess_code = 4

    dso4dt_ait = dmdt_ait * tmp_frso4 / specmw_so4_amode
    dnh4dt_ait = dmdt_ait * (1.0 - tmp_frso4) / specmw_nh4_amode
    return (
        qnuma_del,
        dndt_ait,
        dmdt_ait,
        dso4dt_ait,
        dnh4dt_ait,
        dndt_aitsv1,
        dmdt_aitsv1,
        dndt_aitsv2,
        dmdt_aitsv2,
        dndt_aitsv3,
        dmdt_aitsv3,
        postprocess_code,
    )


@export
def mer07_veh02_nuc_mosaic_postprocess_codon(
    qnuma_del_p: cobj,
    qso4a_del: float,
    qnh4a_del: float,
    deltat: float,
    specmw_so4_amode: float,
    specmw_nh4_amode: float,
    mass1p_aitlo: float,
    mass1p_aithi: float,
    dndt_ait_p: cobj,
    dmdt_ait_p: cobj,
    dso4dt_ait_p: cobj,
    dnh4dt_ait_p: cobj,
    dndt_aitsv1_p: cobj,
    dmdt_aitsv1_p: cobj,
    dndt_aitsv2_p: cobj,
    dmdt_aitsv2_p: cobj,
    dndt_aitsv3_p: cobj,
    dmdt_aitsv3_p: cobj,
    postprocess_code_p: cobj,
):
    qnuma_del = Ptr[float](qnuma_del_p)
    dndt_ait = Ptr[float](dndt_ait_p)
    dmdt_ait = Ptr[float](dmdt_ait_p)
    dso4dt_ait = Ptr[float](dso4dt_ait_p)
    dnh4dt_ait = Ptr[float](dnh4dt_ait_p)
    dndt_aitsv1 = Ptr[float](dndt_aitsv1_p)
    dmdt_aitsv1 = Ptr[float](dmdt_aitsv1_p)
    dndt_aitsv2 = Ptr[float](dndt_aitsv2_p)
    dmdt_aitsv2 = Ptr[float](dmdt_aitsv2_p)
    dndt_aitsv3 = Ptr[float](dndt_aitsv3_p)
    dmdt_aitsv3 = Ptr[float](dmdt_aitsv3_p)
    postprocess_code = Ptr[int](postprocess_code_p)

    (
        qnuma_del[0],
        dndt_ait[0],
        dmdt_ait[0],
        dso4dt_ait[0],
        dnh4dt_ait[0],
        dndt_aitsv1[0],
        dmdt_aitsv1[0],
        dndt_aitsv2[0],
        dmdt_aitsv2[0],
        dndt_aitsv3[0],
        dmdt_aitsv3[0],
        postprocess_code[0],
    ) = _mer07_veh02_nuc_mosaic_postprocess_core(
        qnuma_del[0],
        qso4a_del,
        qnh4a_del,
        deltat,
        specmw_so4_amode,
        specmw_nh4_amode,
        mass1p_aitlo,
        mass1p_aithi,
    )


@export
def mer07_veh02_nuc_mosaic_finalize_codon(
    dtnuc: float,
    temp_in: float,
    rh_in: float,
    press_in: float,
    qh2so4_cur: float,
    qnh3_cur: float,
    h2so4_uptkrate: float,
    mw_so4a_host: float,
    nsize: int,
    dplom_sect_p: cobj,
    dphim_sect_p: cobj,
    cnum_h2so4: float,
    cnum_nh3: float,
    radius_cluster: float,
    so4vol_in: float,
    ratenuclt_bb: float,
    pi_c: float,
    rgas: float,
    avogad: float,
    mw_so4a: float,
    mw_nh4a: float,
    isize_nuc_p: cobj,
    qnuma_del_p: cobj,
    qso4a_del_p: cobj,
    qnh4a_del_p: cobj,
    qh2so4_del_p: cobj,
    qnh3_del_p: cobj,
    dens_nh4so4a_p: cobj,
):
    dplom_sect = Ptr[float](dplom_sect_p)
    dphim_sect = Ptr[float](dphim_sect_p)
    isize_nuc = Ptr[int](isize_nuc_p)
    qnuma_del = Ptr[float](qnuma_del_p)
    qso4a_del = Ptr[float](qso4a_del_p)
    qnh4a_del = Ptr[float](qnh4a_del_p)
    qh2so4_del = Ptr[float](qh2so4_del_p)
    qnh3_del = Ptr[float](qnh3_del_p)
    dens_nh4so4a = Ptr[float](dens_nh4so4a_p)

    onethird = 1.0 / 3.0
    accom_coef_h2so4 = 0.65
    dens_ammsulf = 1.770e3
    dens_ammbisulf = 1.770e3
    dens_sulfacid = 1.770e3
    mw_ammsulf = 132.0
    mw_ammbisulf = 114.0
    mw_sulfacid = 96.0

    cair = press_in / (temp_in * rgas)

    tmpa = max(0.10, min(0.95, rh_in))
    wetvol_dryvol = 1.0 - 0.56 / log(tmpa)

    voldry_clus = (max(cnum_h2so4, 1.0) * mw_so4a + cnum_nh3 * mw_nh4a) / (1.0e3 * dens_sulfacid * avogad)
    voldry_clus = voldry_clus * (mw_so4a_host / mw_so4a)
    dpdry_clus = (voldry_clus * 6.0 / pi_c) ** onethird

    dpdry_part = dplom_sect[0]
    if dpdry_clus <= dplom_sect[0]:
        igrow = 1
    elif dpdry_clus >= dphim_sect[nsize - 1]:
        igrow = 0
        isize_nuc[0] = nsize
        dpdry_part = dphim_sect[nsize - 1]
    else:
        igrow = 0
        for i in range(1, nsize + 1):
            if dpdry_clus < dphim_sect[i - 1]:
                isize_nuc[0] = i
                dpdry_part = dpdry_clus
                dpdry_part = min(dpdry_part, dphim_sect[i - 1])
                dpdry_part = max(dpdry_part, dplom_sect[i - 1])
                break

    voldry_part = (pi_c / 6.0) * (dpdry_part**3)

    if igrow <= 0:
        tmp_n1 = 0.0
        tmp_n2 = 0.0
        tmp_n3 = 1.0
    elif qnh3_cur >= qh2so4_cur:
        tmp_n1 = (qnh3_cur / qh2so4_cur) - 1.0
        tmp_n1 = max(0.0, min(1.0, tmp_n1))
        tmp_n2 = 1.0 - tmp_n1
        tmp_n3 = 0.0
    else:
        tmp_n1 = 0.0
        tmp_n2 = qnh3_cur / qh2so4_cur
        tmp_n2 = max(0.0, min(1.0, tmp_n2))
        tmp_n3 = 1.0 - tmp_n2

    tmp_m1 = tmp_n1 * mw_ammsulf
    tmp_m2 = tmp_n2 * mw_ammbisulf
    tmp_m3 = tmp_n3 * mw_sulfacid
    dens_part = (tmp_m1 + tmp_m2 + tmp_m3) / ((tmp_m1 / dens_ammsulf) + (tmp_m2 / dens_ammbisulf) + (tmp_m3 / dens_sulfacid))
    dens_nh4so4a[0] = dens_part
    mass_part = voldry_part * dens_part
    molenh4a_per_moleso4a = 2.0 * tmp_n1 + tmp_n2
    kgaero_per_moleso4a = 1.0e-3 * (tmp_m1 + tmp_m2 + tmp_m3)
    kgaero_per_moleso4a = kgaero_per_moleso4a * (mw_so4a_host / mw_so4a)

    tmpb = 1.0 + molenh4a_per_moleso4a * 17.0 / 98.0
    wet_volfrac_so4a = 1.0 / (wetvol_dryvol * tmpb)

    if igrow <= 0:
        factor_kk = 1.0
    else:
        tmp_spd = 14.7 * sqrt(temp_in)
        gr_kk = 3.0e-9 * tmp_spd * mw_sulfacid * so4vol_in / (dens_part * wet_volfrac_so4a)
        dfin_kk = 1.0e9 * dpdry_part * (wetvol_dryvol**onethird)
        dnuc_kk = 2.0 * radius_cluster
        dnuc_kk = max(dnuc_kk, 1.0)
        gamma_kk = 0.23 * (dnuc_kk**0.2) * (dfin_kk / 3.0) ** 0.075 * (dens_part * 1.0e-3) ** (-0.33) * (temp_in / 293.0) ** (-0.75)
        tmpa = h2so4_uptkrate * 3600.0
        tmpa = max(tmpa, 0.0)
        tmpb = 6.7037e-6 * (temp_in**0.75) / cair
        tmpb = tmpb * 3600.0
        cs_prime_kk = tmpa / (4.0 * pi_c * tmpb * accom_coef_h2so4)
        nu_kk = gamma_kk * cs_prime_kk / gr_kk
        factor_kk = exp((nu_kk / dfin_kk) - (nu_kk / dnuc_kk))

    ratenuclt_kk = ratenuclt_bb * factor_kk

    tmpa = max(0.0, ratenuclt_kk * dtnuc * mass_part)
    tmpe = tmpa / (kgaero_per_moleso4a * cair)
    qmolso4a_del_max = tmpe

    freducea = 1.0
    if qmolso4a_del_max > qh2so4_cur:
        freducea = qh2so4_cur / qmolso4a_del_max

    freduceb = 1.0
    if molenh4a_per_moleso4a >= 1.0e-10:
        qmolnh4a_del_max = qmolso4a_del_max * molenh4a_per_moleso4a
        if qmolnh4a_del_max > qnh3_cur:
            freduceb = qnh3_cur / qmolnh4a_del_max

    freduce = min(freducea, freduceb)
    if freduce * ratenuclt_kk <= 1.0e-12:
        return

    tmpa = 0.9999
    qh2so4_del[0] = min(tmpa * qh2so4_cur, freduce * qmolso4a_del_max)
    qnh3_del[0] = min(tmpa * qnh3_cur, qh2so4_del[0] * molenh4a_per_moleso4a)
    qh2so4_del[0] = -qh2so4_del[0]
    qnh3_del[0] = -qnh3_del[0]
    qso4a_del[0] = -qh2so4_del[0]
    qnh4a_del[0] = -qnh3_del[0]
    qnuma_del[0] = 1.0e-3 * (qso4a_del[0] * mw_so4a + qnh4a_del[0] * mw_nh4a) / mass_part


def _mer07_veh02_nuc_mosaic_1box_core(
    newnuc_method_flagaa: int,
    dtnuc: float,
    temp_in: float,
    rh_in: float,
    press_in: float,
    zm_in: float,
    pblh_in: float,
    qh2so4_cur: float,
    qh2so4_avg: float,
    qnh3_cur: float,
    h2so4_uptkrate: float,
    mw_so4a_host: float,
    nsize: int,
    dplom_sect: Ptr[float],
    dphim_sect: Ptr[float],
    ldiagaa: int,
    pi_c: float,
    rgas: float,
    avogad: float,
    mw_so4a: float,
    mw_nh4a: float,
):
    if ldiagaa > 0:
        return (1, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    isize_nuc_local = 1
    qnuma_del_local = 0.0
    qso4a_del_local = 0.0
    qnh4a_del_local = 0.0
    qh2so4_del_local = 0.0
    qnh3_del_local = 0.0
    dens_nh4so4a_local = 0.0

    valid_method = 0
    if (
        newnuc_method_flagaa == 1
        or newnuc_method_flagaa == 2
        or newnuc_method_flagaa == 11
        or newnuc_method_flagaa == 12
    ):
        valid_method = 1

    if valid_method == 0:
        return (
            0,
            isize_nuc_local,
            qnuma_del_local,
            qso4a_del_local,
            qnh4a_del_local,
            qh2so4_del_local,
            qnh3_del_local,
            dens_nh4so4a_local,
        )

    cair = press_in / (temp_in * rgas)
    so4vol_in = qh2so4_avg * cair * avogad * 1.0e-6
    nh3ppt = qnh3_cur * 1.0e12
    ratenuclt = 1.0e-38
    rateloge = log(ratenuclt)
    temp_bb = 0.0
    rh_bb = 0.0
    so4vol_bb = 0.0
    nh3ppt_bb = 0.0
    crit_x = 0.0
    cnum_tot = 0.0
    cnum_h2so4 = 0.0
    cnum_nh3 = 0.0
    radius_cluster = 0.0
    use_ternary_rate = 0
    use_binary_rate = 0
    do_pbl_rate = 0

    if (newnuc_method_flagaa != 2) and (nh3ppt >= 0.1):
        if so4vol_in >= 5.0e4:
            temp_bb = max(235.0, min(295.0, temp_in))
            rh_bb = max(0.05, min(0.95, rh_in))
            so4vol_bb = max(5.0e4, min(1.0e9, so4vol_in))
            nh3ppt_bb = max(0.1, min(1.0e3, nh3ppt))
            use_ternary_rate = 1
        newnuc_method_flagaa2 = 1
    else:
        if so4vol_in >= 1.0e4:
            temp_bb = max(230.15, min(305.15, temp_in))
            rh_bb = max(1.0e-4, min(1.0, rh_in))
            so4vol_bb = max(1.0e4, min(1.0e11, so4vol_in))
            use_binary_rate = 1
        newnuc_method_flagaa2 = 2

    if (newnuc_method_flagaa == 11) or (newnuc_method_flagaa == 12):
        if zm_in <= max(pblh_in, 100.0):
            do_pbl_rate = 1

    if use_ternary_rate != 0:
        log_c2 = log(so4vol_bb)
        log_c3 = log(nh3ppt_bb)
        log_rh = log(rh_bb)
        temp_sq = temp_bb**2
        temp_cu = temp_bb**3
        nh3ppt_bb_cu = nh3ppt_bb**3
        log_c2_sq = log_c2**2
        log_c3_sq = log_c3**2
        log_c3_cu = log_c3**3

        t_onset = (
            143.6002929064716
            + 1.0178856665693992 * rh_bb
            + 10.196398812974294 * log_c2
            - 0.1849879416839113 * log_c2_sq
            - 17.161783213150173 * log_c3
            + (109.92469248546053 * log_c3) / log_c2
            + 0.7734119613144357 * log_c2 * log_c3
            - 0.15576469879527022 * log_c3_sq
        )

        if t_onset > temp_bb:
            rateloge = (
                -12.861848898625231
                + 4.905527742256349 * nh3ppt_bb
                - 358.2337705052991 * rh_bb
                - 0.05463019231872484 * nh3ppt_bb * temp_bb
                + 4.8630382337426985 * rh_bb * temp_bb
                + 0.00020258394697064567 * nh3ppt_bb * temp_sq
                - 0.02175548069741675 * rh_bb * temp_sq
                - 2.502406532869512e-7 * nh3ppt_bb * temp_cu
                + 0.00003212869941055865 * rh_bb * temp_cu
                - 4.39129415725234e6 / log_c2_sq
                + (56383.93843154586 * temp_bb) / log_c2_sq
                - (239.835990963361 * temp_sq) / log_c2_sq
                + (0.33765136625580167 * temp_cu) / log_c2_sq
                - (629.7882041830943 * rh_bb) / (nh3ppt_bb_cu * log_c2)
                + (7.772806552631709 * rh_bb * temp_bb) / (nh3ppt_bb_cu * log_c2)
                - (0.031974053936299256 * rh_bb * temp_sq) / (nh3ppt_bb_cu * log_c2)
                + (0.00004383764128775082 * rh_bb * temp_cu) / (nh3ppt_bb_cu * log_c2)
                + 1200.472096232311 * log_c2
                - 17.37107890065621 * temp_bb * log_c2
                + 0.08170681335921742 * temp_sq * log_c2
                - 0.00012534476159729881 * temp_cu * log_c2
                - 14.833042158178936 * log_c2_sq
                + 0.2932631303555295 * temp_bb * log_c2_sq
                - 0.0016497524241142845 * temp_sq * log_c2_sq
                + 2.844074805239367e-6 * temp_cu * log_c2_sq
                - 231375.56676032578 * log_c3
                - 100.21645273730675 * rh_bb * log_c3
                + 2919.2852552424706 * temp_bb * log_c3
                + 0.977886555834732 * rh_bb * temp_bb * log_c3
                - 12.286497122264588 * temp_sq * log_c3
                - 0.0030511783284506377 * rh_bb * temp_sq * log_c3
                + 0.017249301826661612 * temp_cu * log_c3
                + 2.967320346100855e-6 * rh_bb * temp_cu * log_c3
                + (2.360931724951942e6 * log_c3) / log_c2
                - (29752.130254319443 * temp_bb * log_c3) / log_c2
                + (125.04965118142027 * temp_sq * log_c3) / log_c2
                - (0.1752996881934318 * temp_cu * log_c3) / log_c2
                + 5599.912337254629 * log_c2 * log_c3
                - 70.70896612937771 * temp_bb * log_c2 * log_c3
                + 0.2978801613269466 * temp_sq * log_c2 * log_c3
                - 0.00041866525019504 * temp_cu * log_c2 * log_c3
                + 75061.15281456841 * log_c3_sq
                - 931.8802278173565 * temp_bb * log_c3_sq
                + 3.863266220840964 * temp_sq * log_c3_sq
                - 0.005349472062284983 * temp_cu * log_c3_sq
                - (732006.8180571689 * log_c3_sq) / log_c2
                + (9100.06398573816 * temp_bb * log_c3_sq) / log_c2
                - (37.771091915932004 * temp_sq * log_c3_sq) / log_c2
                + (0.05235455395566905 * temp_cu * log_c3_sq) / log_c2
                - 1911.0303773001353 * log_c2 * log_c3_sq
                + 23.6903969622286 * temp_bb * log_c2 * log_c3_sq
                - 0.09807872005428583 * temp_sq * log_c2 * log_c3_sq
                + 0.00013564560238552576 * temp_cu * log_c2 * log_c3_sq
                - 3180.5610833308 * log_c3_cu
                + 39.08268568672095 * temp_bb * log_c3_cu
                - 0.16048521066690752 * temp_sq * log_c3_cu
                + 0.00022031380023793877 * temp_cu * log_c3_cu
                + (40751.075322248245 * log_c3_cu) / log_c2
                - (501.66977622013934 * temp_bb * log_c3_cu) / log_c2
                + (2.063469732254135 * temp_sq * log_c3_cu) / log_c2
                - (0.002836873785758324 * temp_cu * log_c3_cu) / log_c2
                + 2.792313345723013 * log_c2_sq * log_c3_cu
                - 0.03422552111802899 * temp_bb * log_c2_sq * log_c3_cu
                + 0.00014019195277521142 * temp_sq * log_c2_sq * log_c3_cu
                - 1.9201227328396297e-7 * temp_cu * log_c2_sq * log_c3_cu
                - 980.923146020468 * log_rh
                + 10.054155220444462 * temp_bb * log_rh
                - 0.03306644502023841 * temp_sq * log_rh
                + 0.000034274041225891804 * temp_cu * log_rh
                + (16597.75554295064 * log_rh) / log_c2
                - (175.2365504237746 * temp_bb * log_rh) / log_c2
                + (0.6033215603167458 * temp_sq * log_rh) / log_c2
                - (0.0006731787599587544 * temp_cu * log_rh) / log_c2
                - 89.38961120336789 * log_c3 * log_rh
                + 1.153344219304926 * temp_bb * log_c3 * log_rh
                - 0.004954549700267233 * temp_sq * log_c3 * log_rh
                + 7.096309866238719e-6 * temp_cu * log_c3 * log_rh
                + 3.1712136610383244 * log_c3_cu * log_rh
                - 0.037822330602328806 * temp_bb * log_c3_cu * log_rh
                + 0.0001500555743561457 * temp_sq * log_c3_cu * log_rh
                - 1.9828365865570703e-7 * temp_cu * log_c3_cu * log_rh
            )

            j = exp(rateloge)
            log_j = log(j)
            log_j_sq = log_j**2

            cnum_tot = (
                57.40091052369212
                - 0.2996341884645408 * temp_bb
                + 0.0007395477768531926 * temp_sq
                - 5.090604835032423 * log_c2
                + 0.011016634044531128 * temp_bb * log_c2
                + 0.06750032251225707 * log_c2_sq
                - 0.8102831333223962 * log_c3
                + 0.015905081275952426 * temp_bb * log_c3
                - 0.2044174683159531 * log_c2 * log_c3
                + 0.08918159167625832 * log_c3_sq
                - 0.0004969033586666147 * temp_bb * log_c3_sq
                + 0.005704394549007816 * log_c3_cu
                + 3.4098703903474368 * log_j
                - 0.014916956508210809 * temp_bb * log_j
                + 0.08459090011666293 * log_c3 * log_j
                - 0.00014800625143907616 * temp_bb * log_c3 * log_j
                + 0.00503804694656905 * log_j_sq
            )

            radius_cluster = (
                3.2888553966535506e-10
                - 3.374171768439839e-12 * temp_bb
                + 1.8347359507774313e-14 * temp_sq
                + 2.5419844298881856e-12 * log_c2
                - 9.498107643050827e-14 * temp_bb * log_c2
                + 7.446266520834559e-13 * log_c2_sq
                + 2.4303397746137294e-11 * log_c3
                + 1.589324325956633e-14 * temp_bb * log_c3
                - 2.034596219775266e-12 * log_c2 * log_c3
                - 5.59303954457172e-13 * log_c3_sq
                - 4.889507104645867e-16 * temp_bb * log_c3_sq
                + 1.3847024107506764e-13 * log_c3_cu
                + 4.141077193427042e-15 * log_j
                - 2.6813110884009767e-14 * temp_bb * log_j
                + 1.2879071621313094e-12 * log_c3 * log_j
                - 3.80352446061867e-15 * temp_bb * log_c3 * log_j
                - 1.8790172502456827e-14 * log_j_sq
            )

            cnum_h2so4 = (
                -4.7154180661803595
                + 0.13436423483953885 * temp_bb
                - 0.00047184686478816176 * temp_sq
                - 2.564010713640308 * log_c2
                + 0.011353312899114723 * temp_bb * log_c2
                + 0.0010801941974317014 * log_c2_sq
                + 0.5171368624197119 * log_c3
                - 0.0027882479896204665 * temp_bb * log_c3
                + 0.8066971907026886 * log_c3_sq
                - 0.0031849094214409335 * temp_bb * log_c3_sq
                - 0.09951184152927882 * log_c3_cu
                + 0.00040072788891745513 * temp_bb * log_c3_cu
                + 1.3276469271073974 * log_j
                - 0.006167654171986281 * temp_bb * log_j
                - 0.11061390967822708 * log_c3 * log_j
                + 0.0004367575329273496 * temp_bb * log_c3 * log_j
                + 0.000916366357266258 * log_j_sq
            )

            cnum_nh3 = (
                71.20073903979772
                - 0.8409600103431923 * temp_bb
                + 0.0024803006590334922 * temp_sq
                + 2.7798606841602607 * log_c2
                - 0.01475023348171676 * temp_bb * log_c2
                + 0.012264508212031405 * log_c2_sq
                - 2.009926050440182 * log_c3
                + 0.008689123511431527 * temp_bb * log_c3
                - 0.009141180198955415 * log_c2 * log_c3
                + 0.1374122553905617 * log_c3_sq
                - 0.0006253227821679215 * temp_bb * log_c3_sq
                + 0.00009377332742098946 * log_c3_cu
                + 0.5202974341687757 * log_j
                - 0.002419872323052805 * temp_bb * log_j
                + 0.07916392322884074 * log_c3 * log_j
                - 0.0003021586030317366 * temp_bb * log_c3 * log_j
                + 0.0046977006608603395 * log_j_sq
            )
        else:
            rateloge = -300.0
    else:
        if use_binary_rate != 0:
            crit_x = (
                0.740997
                - 0.00266379 * temp_bb
                - 0.00349998 * log(so4vol_bb)
                + 0.0000504022 * temp_bb * log(so4vol_bb)
                + 0.00201048 * log(rh_bb)
                - 0.000183289 * temp_bb * log(rh_bb)
                + 0.00157407 * (log(rh_bb)) ** 2.0
                - 0.0000179059 * temp_bb * (log(rh_bb)) ** 2.0
                + 0.000184403 * (log(rh_bb)) ** 3.0
                - 1.50345e-6 * temp_bb * (log(rh_bb)) ** 3.0
            )

            acoe = 0.14309 + 2.21956 * temp_bb - 0.0273911 * temp_bb**2.0 + 0.0000722811 * temp_bb**3.0 + 5.91822 / crit_x
            bcoe = 0.117489 + 0.462532 * temp_bb - 0.0118059 * temp_bb**2.0 + 0.0000404196 * temp_bb**3.0 + 15.7963 / crit_x
            ccoe = -0.215554 - 0.0810269 * temp_bb + 0.00143581 * temp_bb**2.0 - 4.7758e-6 * temp_bb**3.0 - 2.91297 / crit_x
            dcoe = -3.58856 + 0.049508 * temp_bb - 0.00021382 * temp_bb**2.0 + 3.10801e-7 * temp_bb**3.0 - 0.0293333 / crit_x
            ecoe = 1.14598 - 0.600796 * temp_bb + 0.00864245 * temp_bb**2.0 - 0.0000228947 * temp_bb**3.0 - 8.44985 / crit_x
            fcoe = 2.15855 + 0.0808121 * temp_bb - 0.000407382 * temp_bb**2.0 - 4.01957e-7 * temp_bb**3.0 + 0.721326 / crit_x
            gcoe = 1.6241 - 0.0160106 * temp_bb + 0.0000377124 * temp_bb**2.0 + 3.21794e-8 * temp_bb**3.0 - 0.0113255 / crit_x
            hcoe = 9.71682 - 0.115048 * temp_bb + 0.000157098 * temp_bb**2.0 + 4.00914e-7 * temp_bb**3.0 + 0.71186 / crit_x
            icoe = -1.05611 + 0.00903378 * temp_bb - 0.0000198417 * temp_bb**2.0 + 2.46048e-8 * temp_bb**3.0 - 0.0579087 / crit_x
            jcoe = -0.148712 + 0.00283508 * temp_bb - 9.24619e-6 * temp_bb**2.0 + 5.00427e-9 * temp_bb**3.0 - 0.0127081 / crit_x

            tmpa = (
                acoe
                + bcoe * log(rh_bb)
                + ccoe * (log(rh_bb)) ** 2.0
                + dcoe * (log(rh_bb)) ** 3.0
                + ecoe * log(so4vol_bb)
                + fcoe * (log(rh_bb)) * (log(so4vol_bb))
                + gcoe * ((log(rh_bb)) ** 2.0) * (log(so4vol_bb))
                + hcoe * (log(so4vol_bb)) ** 2.0
                + icoe * log(rh_bb) * ((log(so4vol_bb)) ** 2.0)
                + jcoe * (log(so4vol_bb)) ** 3.0
            )
            rateloge = tmpa
            tmpa = min(tmpa, log(1.0e38))
            ratenuclt = exp(tmpa)

            acoe = -0.00295413 - 0.0976834 * temp_bb + 0.00102485 * temp_bb**2.0 - 2.18646e-6 * temp_bb**3.0 - 0.101717 / crit_x
            bcoe = -0.00205064 - 0.00758504 * temp_bb + 0.000192654 * temp_bb**2.0 - 6.7043e-7 * temp_bb**3.0 - 0.255774 / crit_x
            ccoe = 0.00322308 + 0.000852637 * temp_bb - 0.0000154757 * temp_bb**2.0 + 5.66661e-8 * temp_bb**3.0 + 0.0338444 / crit_x
            dcoe = 0.0474323 - 0.000625104 * temp_bb + 2.65066e-6 * temp_bb**2.0 - 3.67471e-9 * temp_bb**3.0 - 0.000267251 / crit_x
            ecoe = -0.0125211 + 0.00580655 * temp_bb - 0.000101674 * temp_bb**2.0 + 2.88195e-7 * temp_bb**3.0 + 0.0942243 / crit_x
            fcoe = -0.038546 - 0.000672316 * temp_bb + 2.60288e-6 * temp_bb**2.0 + 1.19416e-8 * temp_bb**3.0 - 0.00851515 / crit_x
            gcoe = -0.0183749 + 0.000172072 * temp_bb - 3.71766e-7 * temp_bb**2.0 - 5.14875e-10 * temp_bb**3.0 + 0.00026866 / crit_x
            hcoe = -0.0619974 + 0.000906958 * temp_bb - 9.11728e-7 * temp_bb**2.0 - 5.36796e-9 * temp_bb**3.0 - 0.00774234 / crit_x
            icoe = 0.0121827 - 0.00010665 * temp_bb + 2.5346e-7 * temp_bb**2.0 - 3.63519e-10 * temp_bb**3.0 + 0.000610065 / crit_x
            jcoe = 0.000320184 - 0.0000174762 * temp_bb + 6.06504e-8 * temp_bb**2.0 - 1.4177e-11 * temp_bb**3.0 + 0.000135751 / crit_x

            cnum_tot = exp(
                acoe
                + bcoe * log(rh_bb)
                + ccoe * (log(rh_bb)) ** 2.0
                + dcoe * (log(rh_bb)) ** 3.0
                + ecoe * log(so4vol_bb)
                + fcoe * (log(rh_bb)) * (log(so4vol_bb))
                + gcoe * ((log(rh_bb)) ** 2.0) * (log(so4vol_bb))
                + hcoe * (log(so4vol_bb)) ** 2.0
                + icoe * log(rh_bb) * ((log(so4vol_bb)) ** 2.0)
                + jcoe * (log(so4vol_bb)) ** 3.0
            )

            cnum_h2so4 = cnum_tot * crit_x
            radius_cluster = exp(-1.6524245 + 0.42316402 * crit_x + 0.3346648 * log(cnum_tot))

    if do_pbl_rate != 0:
        if newnuc_method_flagaa == 11:
            tmp_ratenucl = 1.0e-6 * so4vol_in
        elif newnuc_method_flagaa == 12:
            tmp_ratenucl = 1.0e-12 * (so4vol_in * so4vol_in)
        else:
            tmp_ratenucl = -1.0

        if tmp_ratenucl > 0.0:
            tmp_rateloge = log(tmp_ratenucl)
            if tmp_rateloge > rateloge:
                rateloge = tmp_rateloge
                ratenuclt = tmp_ratenucl
                newnuc_method_flagaa2 = newnuc_method_flagaa
                radius_cluster = 0.5
                tmp_diam = radius_cluster * 2.0e-7
                tmp_volu = (tmp_diam * tmp_diam * tmp_diam) * (pi / 6.0)
                tmp_mass = tmp_volu * 1.8
                cnum_h2so4 = (tmp_mass / 98.0) * 6.023e23
                cnum_tot = cnum_h2so4
                cnum_nh3 = 0.0

    ratenuclt_bb = 0.0
    continue_flag = 0
    if rateloge > -13.82:
        ratenuclt = exp(rateloge)
        ratenuclt_bb = ratenuclt * 1.0e6
        continue_flag = 1

    if continue_flag == 0:
        return (
            0,
            isize_nuc_local,
            qnuma_del_local,
            qso4a_del_local,
            qnh4a_del_local,
            qh2so4_del_local,
            qnh3_del_local,
            dens_nh4so4a_local,
        )

    onethird = 1.0 / 3.0
    accom_coef_h2so4 = 0.65
    dens_ammsulf = 1.770e3
    dens_ammbisulf = 1.770e3
    dens_sulfacid = 1.770e3
    mw_ammsulf = 132.0
    mw_ammbisulf = 114.0
    mw_sulfacid = 96.0

    tmpa = max(0.10, min(0.95, rh_in))
    wetvol_dryvol = 1.0 - 0.56 / log(tmpa)

    voldry_clus = (max(cnum_h2so4, 1.0) * mw_so4a + cnum_nh3 * mw_nh4a) / (1.0e3 * dens_sulfacid * avogad)
    voldry_clus = voldry_clus * (mw_so4a_host / mw_so4a)
    dpdry_clus = (voldry_clus * 6.0 / pi_c) ** onethird

    dpdry_part = dplom_sect[0]
    if dpdry_clus <= dplom_sect[0]:
        igrow = 1
    elif dpdry_clus >= dphim_sect[nsize - 1]:
        igrow = 0
        isize_nuc_local = nsize
        dpdry_part = dphim_sect[nsize - 1]
    else:
        igrow = 0
        for i in range(1, nsize + 1):
            if dpdry_clus < dphim_sect[i - 1]:
                isize_nuc_local = i
                dpdry_part = dpdry_clus
                dpdry_part = min(dpdry_part, dphim_sect[i - 1])
                dpdry_part = max(dpdry_part, dplom_sect[i - 1])
                break

    voldry_part = (pi_c / 6.0) * (dpdry_part**3)

    if igrow <= 0:
        tmp_n1 = 0.0
        tmp_n2 = 0.0
        tmp_n3 = 1.0
    elif qnh3_cur >= qh2so4_cur:
        tmp_n1 = (qnh3_cur / qh2so4_cur) - 1.0
        tmp_n1 = max(0.0, min(1.0, tmp_n1))
        tmp_n2 = 1.0 - tmp_n1
        tmp_n3 = 0.0
    else:
        tmp_n1 = 0.0
        tmp_n2 = qnh3_cur / qh2so4_cur
        tmp_n2 = max(0.0, min(1.0, tmp_n2))
        tmp_n3 = 1.0 - tmp_n2

    tmp_m1 = tmp_n1 * mw_ammsulf
    tmp_m2 = tmp_n2 * mw_ammbisulf
    tmp_m3 = tmp_n3 * mw_sulfacid
    dens_part = (tmp_m1 + tmp_m2 + tmp_m3) / ((tmp_m1 / dens_ammsulf) + (tmp_m2 / dens_ammbisulf) + (tmp_m3 / dens_sulfacid))
    dens_nh4so4a_local = dens_part
    mass_part = voldry_part * dens_part
    molenh4a_per_moleso4a = 2.0 * tmp_n1 + tmp_n2
    kgaero_per_moleso4a = 1.0e-3 * (tmp_m1 + tmp_m2 + tmp_m3)
    kgaero_per_moleso4a = kgaero_per_moleso4a * (mw_so4a_host / mw_so4a)

    tmpb = 1.0 + molenh4a_per_moleso4a * 17.0 / 98.0
    wet_volfrac_so4a = 1.0 / (wetvol_dryvol * tmpb)

    if igrow <= 0:
        factor_kk = 1.0
    else:
        tmp_spd = 14.7 * sqrt(temp_in)
        gr_kk = 3.0e-9 * tmp_spd * mw_sulfacid * so4vol_in / (dens_part * wet_volfrac_so4a)
        dfin_kk = 1.0e9 * dpdry_part * (wetvol_dryvol**onethird)
        dnuc_kk = 2.0 * radius_cluster
        dnuc_kk = max(dnuc_kk, 1.0)
        gamma_kk = 0.23 * (dnuc_kk**0.2) * (dfin_kk / 3.0) ** 0.075 * (dens_part * 1.0e-3) ** (-0.33) * (temp_in / 293.0) ** (-0.75)
        tmpa = h2so4_uptkrate * 3600.0
        tmpa = max(tmpa, 0.0)
        tmpb = 6.7037e-6 * (temp_in**0.75) / cair
        tmpb1 = tmpb
        tmpb = tmpb * 3600.0
        cs_prime_kk = tmpa / (4.0 * pi_c * tmpb * accom_coef_h2so4)
        cs_kk = cs_prime_kk * 4.0 * pi_c * tmpb1
        nu_kk = gamma_kk * cs_prime_kk / gr_kk
        factor_kk = exp((nu_kk / dfin_kk) - (nu_kk / dnuc_kk))

    ratenuclt_kk = ratenuclt_bb * factor_kk

    tmpa = max(0.0, ratenuclt_kk * dtnuc * mass_part)
    tmpe = tmpa / (kgaero_per_moleso4a * cair)
    qmolso4a_del_max = tmpe

    freducea = 1.0
    if qmolso4a_del_max > qh2so4_cur:
        freducea = qh2so4_cur / qmolso4a_del_max

    freduceb = 1.0
    if molenh4a_per_moleso4a >= 1.0e-10:
        qmolnh4a_del_max = qmolso4a_del_max * molenh4a_per_moleso4a
        if qmolnh4a_del_max > qnh3_cur:
            freduceb = qnh3_cur / qmolnh4a_del_max

    freduce = min(freducea, freduceb)
    if freduce * ratenuclt_kk <= 1.0e-12:
        return (
            0,
            isize_nuc_local,
            qnuma_del_local,
            qso4a_del_local,
            qnh4a_del_local,
            qh2so4_del_local,
            qnh3_del_local,
            dens_nh4so4a_local,
        )

    tmpa = 0.9999
    qh2so4_del_local = min(tmpa * qh2so4_cur, freduce * qmolso4a_del_max)
    qnh3_del_local = min(tmpa * qnh3_cur, qh2so4_del_local * molenh4a_per_moleso4a)
    qh2so4_del_local = -qh2so4_del_local
    qnh3_del_local = -qnh3_del_local
    qso4a_del_local = -qh2so4_del_local
    qnh4a_del_local = -qnh3_del_local
    qnuma_del_local = 1.0e-3 * (qso4a_del_local * mw_so4a + qnh4a_del_local * mw_nh4a) / mass_part

    return (
        0,
        isize_nuc_local,
        qnuma_del_local,
        qso4a_del_local,
        qnh4a_del_local,
        qh2so4_del_local,
        qnh3_del_local,
        dens_nh4so4a_local,
    )


@export
def mer07_veh02_nuc_mosaic_1box_codon(
    newnuc_method_flagaa: int,
    dtnuc: float,
    temp_in: float,
    rh_in: float,
    press_in: float,
    zm_in: float,
    pblh_in: float,
    qh2so4_cur: float,
    qh2so4_avg: float,
    qnh3_cur: float,
    h2so4_uptkrate: float,
    mw_so4a_host: float,
    nsize: int,
    dplom_sect_p: cobj,
    dphim_sect_p: cobj,
    ldiagaa: int,
    pi_c: float,
    rgas: float,
    avogad: float,
    mw_so4a: float,
    mw_nh4a: float,
    isize_nuc_p: cobj,
    qnuma_del_p: cobj,
    qso4a_del_p: cobj,
    qnh4a_del_p: cobj,
    qh2so4_del_p: cobj,
    qnh3_del_p: cobj,
    dens_nh4so4a_p: cobj,
    fallback_required_p: cobj,
):
    dplom_sect = Ptr[float](dplom_sect_p)
    dphim_sect = Ptr[float](dphim_sect_p)
    isize_nuc = Ptr[int](isize_nuc_p)
    qnuma_del = Ptr[float](qnuma_del_p)
    qso4a_del = Ptr[float](qso4a_del_p)
    qnh4a_del = Ptr[float](qnh4a_del_p)
    qh2so4_del = Ptr[float](qh2so4_del_p)
    qnh3_del = Ptr[float](qnh3_del_p)
    dens_nh4so4a = Ptr[float](dens_nh4so4a_p)
    fallback_required = Ptr[int](fallback_required_p)

    (
        fallback_required[0],
        isize_nuc[0],
        qnuma_del[0],
        qso4a_del[0],
        qnh4a_del[0],
        qh2so4_del[0],
        qnh3_del[0],
        dens_nh4so4a[0],
    ) = _mer07_veh02_nuc_mosaic_1box_core(
        newnuc_method_flagaa,
        dtnuc,
        temp_in,
        rh_in,
        press_in,
        zm_in,
        pblh_in,
        qh2so4_cur,
        qh2so4_avg,
        qnh3_cur,
        h2so4_uptkrate,
        mw_so4a_host,
        nsize,
        dplom_sect,
        dphim_sect,
        ldiagaa,
        pi_c,
        rgas,
        avogad,
        mw_so4a,
        mw_nh4a,
    )


@export
def modal_aero_gasaerexch_snapshot_state_codon(
    ncol: int,
    pver: int,
    pcnstxx: int,
    q_p: cobj,
    qqcw_p: cobj,
    dqdt_p: cobj,
    dqqcwdt_p: cobj,
    qold_p: cobj,
    qqcwold_p: cobj,
    dqdtsv1_p: cobj,
    dqqcwdtsv1_p: cobj,
):
    q = Ptr[float](q_p)
    qqcw = Ptr[float](qqcw_p)
    dqdt = Ptr[float](dqdt_p)
    dqqcwdt = Ptr[float](dqqcwdt_p)
    qold = Ptr[float](qold_p)
    qqcwold = Ptr[float](qqcwold_p)
    dqdtsv1 = Ptr[float](dqdtsv1_p)
    dqqcwdtsv1 = Ptr[float](dqqcwdtsv1_p)

    for m in range(1, pcnstxx + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx3(i, k, m, ncol, pver)
                qold[idx] = q[idx]

    for m in range(1, pcnstxx + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx3(i, k, m, ncol, pver)
                qqcwold[idx] = qqcw[idx]

    for m in range(1, pcnstxx + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx3(i, k, m, ncol, pver)
                dqdtsv1[idx] = dqdt[idx]

    for m in range(1, pcnstxx + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx3(i, k, m, ncol, pver)
                dqqcwdtsv1[idx] = dqqcwdt[idx]


@export
def gas_aer_uptkrates_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    ntot_amode: int,
    q_p: cobj,
    t_p: cobj,
    pmid_p: cobj,
    dgncur_awet_p: cobj,
    numptr_p: cobj,
    sigmag_p: cobj,
    mwdry: float,
    rair: float,
    uptkrate_p: cobj,
):
    q = Ptr[float](q_p)
    t = Ptr[float](t_p)
    pmid = Ptr[float](pmid_p)
    dgncur_awet = Ptr[float](dgncur_awet_p)
    numptr = Ptr[int](numptr_p)
    sigmag = Ptr[float](sigmag_p)
    uptkrate = Ptr[float](uptkrate_p)

    tworootpi = 3.5449077
    root2 = 1.4142135
    beta = 2.0
    xghq0 = 0.70710678
    xghq1 = -0.70710678
    wghq0 = 0.88622693
    wghq1 = 0.88622693

    for n in range(1, ntot_amode + 1):
        lnsg = log(sigmag[n - 1])
        beta_lnsg_sq = beta * (lnsg**2.0)
        half_beta_lnsg_sq = 0.5 * ((beta * lnsg) ** 2.0)
        numptr_idx = numptr[n - 1]

        for k in range(top_lev, pver + 1):
            for i in range(1, ncol + 1):
                temp = t[_idx2(i, k, pcols)]
                pmid_ik = pmid[_idx2(i, k, pcols)]
                rhoair = pmid_ik / (rair * temp)
                aircon = rhoair / mwdry
                num_a = q[_idx3(i, k, numptr_idx, ncol, pver)] * aircon

                gasdiffus = 0.557e-4 * (temp**1.75) / pmid_ik
                gasspeed = 1.470e1 * sqrt(temp)
                freepathx2 = 6.0 * gasdiffus / gasspeed

                lndpgn = log(dgncur_awet[_idx3(i, k, n, pcols, pver)])
                const = tworootpi * num_a * exp(beta * lndpgn + half_beta_lnsg_sq)

                lndp = lndpgn + beta_lnsg_sq + root2 * lnsg * xghq0
                dp = exp(lndp)
                knudsen = freepathx2 / dp
                fuchs_sutugin = (0.4875 * (1.0 + knudsen)) / (
                    knudsen * (1.184 + knudsen) + 0.4875
                )
                sumghq = wghq0 * dp * fuchs_sutugin / (dp**beta)

                lndp = lndpgn + beta_lnsg_sq + root2 * lnsg * xghq1
                dp = exp(lndp)
                knudsen = freepathx2 / dp
                fuchs_sutugin = (0.4875 * (1.0 + knudsen)) / (
                    knudsen * (1.184 + knudsen) + 0.4875
                )
                sumghq += wghq1 * dp * fuchs_sutugin / (dp**beta)

                uptkrate[_idx3(n, i, k, ntot_amode, pcols)] = const * gasdiffus * sumghq


@export
def modal_aero_soaexch_codon(
    dtfull: float,
    temp: float,
    pres: float,
    niter_max: int,
    ntot_soamode: int,
    g_soa_in: float,
    a_soa_in_p: cobj,
    a_poa_in_p: cobj,
    xferrate_p: cobj,
    rgas: float,
    a_opoa_p: cobj,
    a_soa_p: cobj,
    beta_p: cobj,
    g_star_p: cobj,
    phi_p: cobj,
    sat_p: cobj,
    niter_p: cobj,
    g_soa_tend_p: cobj,
    a_soa_tend_p: cobj,
):
    a_soa_in = Ptr[float](a_soa_in_p)
    a_poa_in = Ptr[float](a_poa_in_p)
    xferrate = Ptr[float](xferrate_p)
    a_opoa = Ptr[float](a_opoa_p)
    a_soa = Ptr[float](a_soa_p)
    beta = Ptr[float](beta_p)
    g_star = Ptr[float](g_star_p)
    phi = Ptr[float](phi_p)
    sat = Ptr[float](sat_p)
    niter_out = Ptr[int](niter_p)
    g_soa_tend_out = Ptr[float](g_soa_tend_p)
    a_soa_tend = Ptr[float](a_soa_tend_p)

    alpha = 0.05
    g_min1 = 1.0e-20
    opoa_frac = 0.1
    delh_vap_soa = 156.0e3
    p0_soa_298 = 1.0e-10

    g_soa = g_soa_in
    if g_soa < 0.0:
        g_soa = 0.0
    tot_soa = g_soa

    for m in range(1, ntot_soamode + 1):
        a_soa_val = a_soa_in[m - 1]
        if a_soa_val < 0.0:
            a_soa_val = 0.0
        a_soa[m - 1] = a_soa_val
        tot_soa += a_soa_val

        a_opoa_val = opoa_frac * a_poa_in[m - 1]
        if a_opoa_val < 1.0e-20:
            a_opoa_val = 1.0e-20
        a_opoa[m - 1] = a_opoa_val

    p0_soa = p0_soa_298 * exp(
        -(delh_vap_soa / rgas) * ((1.0 / temp) - (1.0 / 298.0))
    )
    g0_soa = 1.01325e5 * p0_soa / pres
    g0_soa = g0_soa * (150.0 / 12.0)

    niter = 0
    tcur = 0.0
    dtcur = 0.0
    for m in range(1, ntot_soamode + 1):
        phi[m - 1] = 0.0
        g_star[m - 1] = 0.0

    while tcur < dtfull - 1.0e-3:
        niter += 1
        if niter > niter_max:
            break

        tmpa = 0.0
        for m in range(1, ntot_soamode + 1):
            sat[m - 1] = g0_soa / (a_soa[m - 1] + a_opoa[m - 1])
            g_star[m - 1] = sat[m - 1] * a_soa[m - 1]
            denom = _max3(g_soa, g_star[m - 1], g_min1)
            phi[m - 1] = (g_soa - g_star[m - 1]) / denom
            tmpa += xferrate[m - 1] * abs(phi[m - 1])

        dtmax = dtfull - tcur
        if dtmax * tmpa <= alpha:
            dtcur = dtmax
            tcur = dtfull
        else:
            dtcur = alpha / tmpa
            tcur += dtcur

        for m in range(1, ntot_soamode + 1):
            beta[m - 1] = dtcur * xferrate[m - 1]
            tmpa = g_soa - g_star[m - 1]
            if tmpa > 0.0:
                a_soa_tmp = a_soa[m - 1] + beta[m - 1] * tmpa
                sat[m - 1] = g0_soa / (a_soa_tmp + a_opoa[m - 1])
                g_star[m - 1] = sat[m - 1] * a_soa_tmp

        tmpa = 0.0
        tmpb = 0.0
        for m in range(1, ntot_soamode + 1):
            denom = 1.0 + beta[m - 1] * sat[m - 1]
            tmpa += a_soa[m - 1] / denom
            tmpb += beta[m - 1] / denom

        g_soa = (tot_soa - tmpa) / (1.0 + tmpb)
        if g_soa < 0.0:
            g_soa = 0.0
        for m in range(1, ntot_soamode + 1):
            a_soa[m - 1] = (a_soa[m - 1] + beta[m - 1] * g_soa) / (
                1.0 + beta[m - 1] * sat[m - 1]
            )

    g_soa_tend_out[0] = (g_soa - g_soa_in) / dtfull
    for m in range(1, ntot_soamode + 1):
        a_soa_tend[m - 1] = (a_soa[m - 1] - a_soa_in[m - 1]) / dtfull
    niter_out[0] = niter


@export
def modal_aero_gasaerexch_sub_codon(
    stage: int,
    ncol: int,
    pcols: int,
    pver: int,
    pcnstxx: int,
    pcnst: int,
    top_lev: int,
    ntot_amode: int,
    maxd_aspectype: int,
    maxspec_pcage: int,
    nsrflx: int,
    jsrflx_gaexch: int,
    jsrflx_rename: int,
    l_so4g: int,
    l_nh4g: int,
    l_msag: int,
    l_soag: int,
    do_nh4g: int,
    do_msag: int,
    do_soag: int,
    method_soa: int,
    ntot_soamode: int,
    modefrm_pcage: int,
    modetoo_pcage: int,
    modeptr_pcarbon: int,
    nspecfrm_pcage: int,
    has_sulfeq: int,
    deltat: float,
    deltatxx: float,
    gravit: float,
    mwdry: float,
    rgas: float,
    fac_m2v_so4: float,
    fac_m2v_nh4: float,
    fac_m2v_soa: float,
    fac_volsfc_pcarbon: float,
    xferfrac_max: float,
    soa_equivso4_factor: float,
    dr_so4_monolayers_pcage: float,
    q_p: cobj,
    qqcw_p: cobj,
    dqdt_p: cobj,
    dqqcwdt_p: cobj,
    qsrflx_p: cobj,
    qqcwsrflx_p: cobj,
    qold_p: cobj,
    qqcwold_p: cobj,
    dqdtsv1_p: cobj,
    dqqcwdtsv1_p: cobj,
    t_p: cobj,
    pmid_p: cobj,
    pdel_p: cobj,
    dgncur_a_p: cobj,
    uptkrate_p: cobj,
    sulfeq_p: cobj,
    troplev_p: cobj,
    ido_so4a_p: cobj,
    ido_nh4a_p: cobj,
    ido_soaa_p: cobj,
    lptr_so4_a_p: cobj,
    lptr_nh4_a_p: cobj,
    lptr_soa_a_p: cobj,
    lptr_pom_a_p: cobj,
    lmassptr_p: cobj,
    nspec_amode_p: cobj,
    lspecfrm_p: cobj,
    lspectoo_p: cobj,
    fac_m2v_pcarbon_p: cobj,
    dqdt_so4_p: cobj,
    dqdt_nh4_p: cobj,
    dqdt_soa_p: cobj,
    fgain_so4_p: cobj,
    fgain_nh4_p: cobj,
    fgain_soa_p: cobj,
    qold_so4_p: cobj,
    qold_nh4_p: cobj,
    qold_soa_p: cobj,
    qold_poa_p: cobj,
    uptkratebb_p: cobj,
    uptkrate_soa_p: cobj,
    a_opoa_soa_p: cobj,
    a_soa_work_p: cobj,
    beta_soa_p: cobj,
    g_star_soa_p: cobj,
    phi_soa_p: cobj,
    sat_soa_p: cobj,
    niter_soa_p: cobj,
    g_soa_tend_p: cobj,
    adv_mass_p: cobj,
    dotend_p: cobj,
    dotendqqcw_p: cobj,
    dotendrn_p: cobj,
    dotendqqcwrn_p: cobj,
):
    q = Ptr[float](q_p)
    qqcw = Ptr[float](qqcw_p)
    dqdt = Ptr[float](dqdt_p)
    dqqcwdt = Ptr[float](dqqcwdt_p)
    qsrflx = Ptr[float](qsrflx_p)
    qqcwsrflx = Ptr[float](qqcwsrflx_p)
    t = Ptr[float](t_p)
    pmid = Ptr[float](pmid_p)
    pdel = Ptr[float](pdel_p)
    dgncur_a = Ptr[float](dgncur_a_p)
    uptkrate = Ptr[float](uptkrate_p)
    sulfeq = Ptr[float](sulfeq_p)
    troplev = Ptr[int](troplev_p)
    ido_so4a = Ptr[int](ido_so4a_p)
    ido_nh4a = Ptr[int](ido_nh4a_p)
    ido_soaa = Ptr[int](ido_soaa_p)
    lptr_so4_a = Ptr[int](lptr_so4_a_p)
    lptr_nh4_a = Ptr[int](lptr_nh4_a_p)
    lptr_soa_a = Ptr[int](lptr_soa_a_p)
    lptr_pom_a = Ptr[int](lptr_pom_a_p)
    lmassptr = Ptr[int](lmassptr_p)
    nspec_amode = Ptr[int](nspec_amode_p)
    lspecfrm = Ptr[int](lspecfrm_p)
    lspectoo = Ptr[int](lspectoo_p)
    fac_m2v_pcarbon = Ptr[float](fac_m2v_pcarbon_p)
    dqdt_so4 = Ptr[float](dqdt_so4_p)
    dqdt_nh4 = Ptr[float](dqdt_nh4_p)
    dqdt_soa = Ptr[float](dqdt_soa_p)
    fgain_so4 = Ptr[float](fgain_so4_p)
    fgain_nh4 = Ptr[float](fgain_nh4_p)
    fgain_soa = Ptr[float](fgain_soa_p)
    qold_so4 = Ptr[float](qold_so4_p)
    qold_nh4 = Ptr[float](qold_nh4_p)
    qold_soa = Ptr[float](qold_soa_p)
    qold_poa = Ptr[float](qold_poa_p)
    uptkratebb = Ptr[float](uptkratebb_p)
    uptkrate_soa = Ptr[float](uptkrate_soa_p)
    g_soa_tend = Ptr[float](g_soa_tend_p)
    adv_mass = Ptr[float](adv_mass_p)
    dotend = Ptr[int](dotend_p)
    dotendqqcw = Ptr[int](dotendqqcw_p)
    dotendrn = Ptr[int](dotendrn_p)
    dotendqqcwrn = Ptr[int](dotendqqcwrn_p)

    if stage == 1:
        modal_aero_gasaerexch_zero_tendencies_codon(
            ncol, pcols, pver, pcnstxx, nsrflx, dqdt_p, dqqcwdt_p, qsrflx_p, qqcwsrflx_p
        )
        return

    if stage == 3:
        modal_aero_gasaerexch_snapshot_state_codon(
            ncol, pver, pcnstxx, q_p, qqcw_p, dqdt_p, dqqcwdt_p, qold_p, qqcwold_p, dqdtsv1_p, dqqcwdtsv1_p
        )
        return

    if stage == 4:
        for l in range(1, pcnstxx + 1):
            if dotend[l - 1] != 0 or dotendrn[l - 1] != 0:
                for k in range(top_lev, pver + 1):
                    for i in range(1, ncol + 1):
                        idx = _idx3(i, k, l, ncol, pver)
                        q[idx] = q[idx] + dqdt[idx] * deltat

            if dotendqqcw[l - 1] != 0 or dotendqqcwrn[l - 1] != 0:
                for k in range(top_lev, pver + 1):
                    for i in range(1, ncol + 1):
                        idx = _idx3(i, k, l, ncol, pver)
                        qqcw[idx] = qqcw[idx] + dqqcwdt[idx] * deltat
        return

    if stage == 5:
        for l in range(1, pcnstxx + 1):
            for jsrf in range(1, nsrflx + 1):
                if jsrf == jsrflx_gaexch:
                    if dotend[l - 1] == 0:
                        continue
                elif jsrf == jsrflx_rename:
                    if dotendrn[l - 1] == 0:
                        continue
                else:
                    continue
                for i in range(1, ncol + 1):
                    idx = _idx3(i, l, jsrf, pcols, pcnstxx)
                    qsrflx[idx] = qsrflx[idx] * (adv_mass[l - 1] / mwdry)

        for l in range(1, pcnstxx + 1):
            if dotendqqcwrn[l - 1] == 0:
                continue
            for i in range(1, ncol + 1):
                idx = _idx3(i, l, jsrflx_rename, pcols, pcnstxx)
                qqcwsrflx[idx] = qqcwsrflx[idx] * (adv_mass[l - 1] / mwdry)
        return

    if stage != 2:
        return

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            sum_uprt_so4 = 0.0
            sum_uprt_nh4 = 0.0
            sum_uprt_soa = 0.0

            for n in range(1, ntot_amode + 1):
                n0 = n - 1
                uptkratebb[n0] = uptkrate[_idx3(n, i, k, ntot_amode, pcols)]

                if ido_so4a[n0] > 0:
                    fgain_so4[n0] = uptkratebb[n0]
                    sum_uprt_so4 = sum_uprt_so4 + fgain_so4[n0]
                    if ido_so4a[n0] == 1:
                        qold_so4[n0] = q[_idx3(i, k, lptr_so4_a[n0], ncol, pver)]
                    else:
                        qold_so4[n0] = 0.0
                else:
                    fgain_so4[n0] = 0.0
                    qold_so4[n0] = 0.0

                if ido_nh4a[n0] > 0:
                    fgain_nh4[n0] = uptkratebb[n0] * 2.08
                    sum_uprt_nh4 = sum_uprt_nh4 + fgain_nh4[n0]
                    if ido_nh4a[n0] == 1:
                        qold_nh4[n0] = q[_idx3(i, k, lptr_nh4_a[n0], ncol, pver)]
                    else:
                        qold_nh4[n0] = 0.0
                else:
                    fgain_nh4[n0] = 0.0
                    qold_nh4[n0] = 0.0

                if ido_soaa[n0] > 0:
                    fgain_soa[n0] = uptkratebb[n0] * 0.81
                    sum_uprt_soa = sum_uprt_soa + fgain_soa[n0]
                    if ido_soaa[n0] == 1:
                        qold_soa[n0] = q[_idx3(i, k, lptr_soa_a[n0], ncol, pver)]
                        l = lptr_pom_a[n0]
                        if l > 0:
                            qold_poa[n0] = q[_idx3(i, k, l, ncol, pver)]
                        else:
                            qold_poa[n0] = 0.0
                    else:
                        qold_soa[n0] = 0.0
                        qold_poa[n0] = 0.0
                else:
                    fgain_soa[n0] = 0.0
                    qold_soa[n0] = 0.0
                    qold_poa[n0] = 0.0
                uptkrate_soa[n0] = fgain_soa[n0]

            if sum_uprt_so4 > 0.0:
                for n in range(1, ntot_amode + 1):
                    fgain_so4[n - 1] = fgain_so4[n - 1] / sum_uprt_so4

            if sum_uprt_nh4 > 0.0:
                for n in range(1, ntot_amode + 1):
                    fgain_nh4[n - 1] = fgain_nh4[n - 1] / sum_uprt_nh4

            if sum_uprt_soa > 0.0:
                for n in range(1, ntot_amode + 1):
                    fgain_soa[n - 1] = fgain_soa[n - 1] / sum_uprt_soa

            avg_uprt_so4 = (1.0 - exp(-deltatxx * sum_uprt_so4)) / deltatxx
            avg_uprt_nh4 = (1.0 - exp(-deltatxx * sum_uprt_nh4)) / deltatxx
            avg_uprt_soa = (1.0 - exp(-deltatxx * sum_uprt_soa)) / deltatxx

            sum_dqdt_so4 = q[_idx3(i, k, l_so4g, ncol, pver)] * avg_uprt_so4
            if do_msag != 0:
                sum_dqdt_msa = q[_idx3(i, k, l_msag, ncol, pver)] * avg_uprt_so4
            else:
                sum_dqdt_msa = 0.0
            if do_nh4g != 0:
                sum_dqdt_nh4 = q[_idx3(i, k, l_nh4g, ncol, pver)] * avg_uprt_nh4
            else:
                sum_dqdt_nh4 = 0.0
            if do_soag != 0:
                sum_dqdt_soa = q[_idx3(i, k, l_soag, ncol, pver)] * avg_uprt_soa
            else:
                sum_dqdt_soa = 0.0

            if has_sulfeq != 0 and k <= troplev[i - 1]:
                tmp_kxt = deltatxx * sum_uprt_so4
                tmp_pxt = 0.0
                for n in range(1, ntot_amode + 1):
                    if ido_so4a[n - 1] <= 0:
                        continue
                    tmp_pxt = tmp_pxt + uptkratebb[n - 1] * sulfeq[_idx3(i, k, n, pcols, pver)]
                tmp_pxt = max(0.0, tmp_pxt * deltatxx)
                tmp_so4g_bgn = q[_idx3(i, k, l_so4g, ncol, pver)]
                if tmp_kxt >= 1.0e-5:
                    tmp_so4g_equ = tmp_pxt / tmp_kxt
                    tmp_so4g_avg = tmp_so4g_equ + (tmp_so4g_bgn - tmp_so4g_equ) * (1.0 - exp(-tmp_kxt)) / tmp_kxt
                else:
                    tmp_so4g_avg = tmp_so4g_bgn * (1.0 - 0.5 * tmp_kxt) + 0.5 * tmp_pxt
                sum_dqdt_so4 = 0.0
                for n in range(1, ntot_amode + 1):
                    n0 = n - 1
                    if ido_so4a[n0] <= 0:
                        continue
                    if ido_so4a[n0] == 1:
                        l = lptr_so4_a[n0]
                        tmp_so4a_bgn = q[_idx3(i, k, l, ncol, pver)]
                    else:
                        tmp_so4a_bgn = 0.0
                    tmp_so4a_end = tmp_so4a_bgn + deltatxx * uptkratebb[n0] * (
                        tmp_so4g_avg - sulfeq[_idx3(i, k, n, pcols, pver)]
                    )
                    tmp_so4a_end = max(0.0, tmp_so4a_end)
                    dqdt_so4[n0] = (tmp_so4a_end - tmp_so4a_bgn) / deltatxx
                    sum_dqdt_so4 = sum_dqdt_so4 + dqdt_so4[n0]
                if do_msag != 0:
                    sum_dqdt_msa = 0.0
            else:
                for n in range(1, ntot_amode + 1):
                    dqdt_so4[n - 1] = fgain_so4[n - 1] * (sum_dqdt_so4 + sum_dqdt_msa)

            sum_dqdt_nh4_b = 0.0
            for n in range(1, ntot_amode + 1):
                dqdt_nh4[n - 1] = 0.0
            if do_nh4g != 0:
                for n in range(1, ntot_amode + 1):
                    n0 = n - 1
                    dqdt_nh4[n0] = fgain_nh4[n0] * sum_dqdt_nh4
                    qnew_nh4 = qold_nh4[n0] + dqdt_nh4[n0] * deltat
                    qnew_so4 = qold_so4[n0] + dqdt_so4[n0] * deltat
                    qmax_nh4 = 2.0 * qnew_so4
                    if qnew_nh4 > qmax_nh4:
                        dqdt_nh4[n0] = (qmax_nh4 - qold_nh4[n0]) / deltatxx
                    sum_dqdt_nh4_b = sum_dqdt_nh4_b + dqdt_nh4[n0]

            if do_soag != 0 and method_soa > 1:
                niter_max = 1000
                for n in range(1, ntot_amode + 1):
                    dqdt_soa[n - 1] = 0.0
                modal_aero_soaexch_codon(
                    deltat,
                    t[_idx2(i, k, pcols)],
                    pmid[_idx2(i, k, pcols)],
                    niter_max,
                    ntot_soamode,
                    q[_idx3(i, k, l_soag, ncol, pver)],
                    qold_soa_p,
                    qold_poa_p,
                    uptkrate_soa_p,
                    rgas,
                    a_opoa_soa_p,
                    a_soa_work_p,
                    beta_soa_p,
                    g_star_soa_p,
                    phi_soa_p,
                    sat_soa_p,
                    niter_soa_p,
                    g_soa_tend_p,
                    dqdt_soa_p,
                )
                sum_dqdt_soa = -g_soa_tend[0]
            elif do_soag != 0:
                for n in range(1, ntot_amode + 1):
                    dqdt_soa[n - 1] = fgain_soa[n - 1] * sum_dqdt_soa
            else:
                for n in range(1, ntot_amode + 1):
                    dqdt_soa[n - 1] = 0.0

            pdel_fac = pdel[_idx2(i, k, pcols)] / gravit
            for n in range(1, ntot_amode + 1):
                n0 = n - 1
                if ido_so4a[n0] == 1:
                    l = lptr_so4_a[n0]
                    dqdt[_idx3(i, k, l, ncol, pver)] = dqdt_so4[n0]
                    qsrflx[_idx3(i, l, jsrflx_gaexch, pcols, pcnstxx)] = (
                        qsrflx[_idx3(i, l, jsrflx_gaexch, pcols, pcnstxx)] + dqdt_so4[n0] * pdel_fac
                    )

                if do_nh4g != 0:
                    if ido_nh4a[n0] == 1:
                        l = lptr_nh4_a[n0]
                        dqdt[_idx3(i, k, l, ncol, pver)] = dqdt_nh4[n0]
                        qsrflx[_idx3(i, l, jsrflx_gaexch, pcols, pcnstxx)] = (
                            qsrflx[_idx3(i, l, jsrflx_gaexch, pcols, pcnstxx)] + dqdt_nh4[n0] * pdel_fac
                        )

                if do_soag != 0:
                    if ido_soaa[n0] == 1:
                        l = lptr_soa_a[n0]
                        dqdt[_idx3(i, k, l, ncol, pver)] = dqdt_soa[n0]
                        qsrflx[_idx3(i, l, jsrflx_gaexch, pcols, pcnstxx)] = (
                            qsrflx[_idx3(i, l, jsrflx_gaexch, pcols, pcnstxx)] + dqdt_soa[n0] * pdel_fac
                        )

            l = l_so4g
            dqdt[_idx3(i, k, l, ncol, pver)] = -sum_dqdt_so4
            qsrflx[_idx3(i, l, jsrflx_gaexch, pcols, pcnstxx)] = (
                qsrflx[_idx3(i, l, jsrflx_gaexch, pcols, pcnstxx)]
                + dqdt[_idx3(i, k, l, ncol, pver)] * pdel_fac
            )

            if do_msag != 0:
                l = l_msag
                dqdt[_idx3(i, k, l, ncol, pver)] = -sum_dqdt_msa
                qsrflx[_idx3(i, l, jsrflx_gaexch, pcols, pcnstxx)] = (
                    qsrflx[_idx3(i, l, jsrflx_gaexch, pcols, pcnstxx)]
                    + dqdt[_idx3(i, k, l, ncol, pver)] * pdel_fac
                )

            if do_nh4g != 0:
                l = l_nh4g
                dqdt[_idx3(i, k, l, ncol, pver)] = -sum_dqdt_nh4_b
                qsrflx[_idx3(i, l, jsrflx_gaexch, pcols, pcnstxx)] = (
                    qsrflx[_idx3(i, l, jsrflx_gaexch, pcols, pcnstxx)]
                    + dqdt[_idx3(i, k, l, ncol, pver)] * pdel_fac
                )

            if do_soag != 0:
                l = l_soag
                dqdt[_idx3(i, k, l, ncol, pver)] = -sum_dqdt_soa
                qsrflx[_idx3(i, l, jsrflx_gaexch, pcols, pcnstxx)] = (
                    qsrflx[_idx3(i, l, jsrflx_gaexch, pcols, pcnstxx)]
                    + dqdt[_idx3(i, k, l, ncol, pver)] * pdel_fac
                )

            if modefrm_pcage > 0:
                n = modeptr_pcarbon
                n0 = n - 1
                vol_shell = deltat * (
                    dqdt_so4[n0] * fac_m2v_so4
                    + dqdt_nh4[n0] * fac_m2v_nh4
                    + dqdt_soa[n0] * fac_m2v_soa * soa_equivso4_factor
                )
                vol_core = 0.0
                for l in range(1, nspec_amode[n0] + 1):
                    vol_core = vol_core + q[_idx3(i, k, lmassptr[_idx2(l, n, maxd_aspectype)], ncol, pver)] * fac_m2v_pcarbon[l - 1]

                tmp1 = vol_shell * dgncur_a[_idx3(i, k, n, pcols, pver)] * fac_volsfc_pcarbon
                tmp2 = max(6.0 * dr_so4_monolayers_pcage * vol_core, 0.0)
                if tmp1 >= tmp2:
                    xferfrac_pcage = xferfrac_max
                else:
                    xferfrac_pcage = min(tmp1 / tmp2, xferfrac_max)

                if xferfrac_pcage > 0.0:
                    for iq in range(1, nspecfrm_pcage + 1):
                        lsfrm = lspecfrm[iq - 1]
                        lstoo = lspectoo[iq - 1]
                        xferrate = (xferfrac_pcage / deltat) * q[_idx3(i, k, lsfrm, ncol, pver)]
                        dqdt[_idx3(i, k, lsfrm, ncol, pver)] = dqdt[_idx3(i, k, lsfrm, ncol, pver)] - xferrate
                        qsrflx[_idx3(i, lsfrm, jsrflx_gaexch, pcols, pcnstxx)] = (
                            qsrflx[_idx3(i, lsfrm, jsrflx_gaexch, pcols, pcnstxx)] - xferrate * pdel_fac
                        )
                        if lstoo > 0 and lstoo <= pcnst:
                            dqdt[_idx3(i, k, lstoo, ncol, pver)] = dqdt[_idx3(i, k, lstoo, ncol, pver)] + xferrate
                            qsrflx[_idx3(i, lstoo, jsrflx_gaexch, pcols, pcnstxx)] = (
                                qsrflx[_idx3(i, lstoo, jsrflx_gaexch, pcols, pcnstxx)] + xferrate * pdel_fac
                            )

                    if ido_so4a[modetoo_pcage - 1] > 0:
                        l = lptr_so4_a[modetoo_pcage - 1]
                        dqdt[_idx3(i, k, l, ncol, pver)] = dqdt[_idx3(i, k, l, ncol, pver)] + dqdt_so4[modefrm_pcage - 1]
                        qsrflx[_idx3(i, l, jsrflx_gaexch, pcols, pcnstxx)] = (
                            qsrflx[_idx3(i, l, jsrflx_gaexch, pcols, pcnstxx)]
                            + dqdt_so4[modefrm_pcage - 1] * pdel_fac
                        )

                    if ido_nh4a[modetoo_pcage - 1] > 0:
                        l = lptr_nh4_a[modetoo_pcage - 1]
                        dqdt[_idx3(i, k, l, ncol, pver)] = dqdt[_idx3(i, k, l, ncol, pver)] + dqdt_nh4[modefrm_pcage - 1]
                        qsrflx[_idx3(i, l, jsrflx_gaexch, pcols, pcnstxx)] = (
                            qsrflx[_idx3(i, l, jsrflx_gaexch, pcols, pcnstxx)]
                            + dqdt_nh4[modefrm_pcage - 1] * pdel_fac
                        )

                    if ido_soaa[modetoo_pcage - 1] > 0:
                        l = lptr_soa_a[modetoo_pcage - 1]
                        dqdt[_idx3(i, k, l, ncol, pver)] = dqdt[_idx3(i, k, l, ncol, pver)] + dqdt_soa[modefrm_pcage - 1]
                        qsrflx[_idx3(i, l, jsrflx_gaexch, pcols, pcnstxx)] = (
                            qsrflx[_idx3(i, l, jsrflx_gaexch, pcols, pcnstxx)]
                            + dqdt_soa[modefrm_pcage - 1] * pdel_fac
                        )


@export
def modal_aero_rename_no_acc_crs_dryvols_codon(
    ncol: int,
    pver: int,
    pcnstxx: int,
    ntot_amode: int,
    maxspec_renamexf: int,
    loffset: int,
    deltat: float,
    idomode_p: cobj,
    q_p: cobj,
    qqcw_p: cobj,
    dqdt_p: cobj,
    dqdt_other_p: cobj,
    dqqcwdt_p: cobj,
    dqqcwdt_other_p: cobj,
    nspec_amode_p: cobj,
    lspectype_amode_p: cobj,
    specmw_amode_p: cobj,
    specdens_amode_p: cobj,
    lmassptr_amode_p: cobj,
    lmassptrcw_amode_p: cobj,
    dryvol_a_p: cobj,
    dryvol_c_p: cobj,
    deldryvol_a_p: cobj,
    deldryvol_c_p: cobj,
):
    idomode = Ptr[int](idomode_p)
    q = Ptr[float](q_p)
    qqcw = Ptr[float](qqcw_p)
    dqdt = Ptr[float](dqdt_p)
    dqdt_other = Ptr[float](dqdt_other_p)
    dqqcwdt = Ptr[float](dqqcwdt_p)
    dqqcwdt_other = Ptr[float](dqqcwdt_other_p)
    nspec_amode = Ptr[int](nspec_amode_p)
    lspectype_amode = Ptr[int](lspectype_amode_p)
    specmw_amode = Ptr[float](specmw_amode_p)
    specdens_amode = Ptr[float](specdens_amode_p)
    lmassptr_amode = Ptr[int](lmassptr_amode_p)
    lmassptrcw_amode = Ptr[int](lmassptrcw_amode_p)
    dryvol_a = Ptr[float](dryvol_a_p)
    dryvol_c = Ptr[float](dryvol_c_p)
    deldryvol_a = Ptr[float](deldryvol_a_p)
    deldryvol_c = Ptr[float](deldryvol_c_p)

    for n in range(1, ntot_amode + 1):
        if idomode[n - 1] > 0:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    dryvol_a[_idx3(i, k, n, ncol, pver)] = 0.0
                    dryvol_c[_idx3(i, k, n, ncol, pver)] = 0.0
                    deldryvol_a[_idx3(i, k, n, ncol, pver)] = 0.0
                    deldryvol_c[_idx3(i, k, n, ncol, pver)] = 0.0

            for l1 in range(1, nspec_amode[n - 1] + 1):
                l2 = lspectype_amode[_idx2(l1, n, maxspec_renamexf)]
                dum_m2v = specmw_amode[l2 - 1] / specdens_amode[l2 - 1]
                dum_m2vdt = dum_m2v * deltat

                la = lmassptr_amode[_idx2(l1, n, maxspec_renamexf)] - loffset
                if la > 0:
                    for k in range(1, pver + 1):
                        for i in range(1, ncol + 1):
                            qold = q[_idx3(i, k, la, ncol, pver)] - deltat * dqdt_other[
                                _idx3(i, k, la, ncol, pver)
                            ]
                            if qold < 0.0:
                                qold = 0.0
                            dryvol_a[_idx3(i, k, n, ncol, pver)] += dum_m2v * qold
                            deldryvol_a[_idx3(i, k, n, ncol, pver)] += (
                                dqdt_other[_idx3(i, k, la, ncol, pver)]
                                + dqdt[_idx3(i, k, la, ncol, pver)]
                            ) * dum_m2vdt

                lc = lmassptrcw_amode[_idx2(l1, n, maxspec_renamexf)] - loffset
                if lc > 0:
                    for k in range(1, pver + 1):
                        for i in range(1, ncol + 1):
                            qqcwold = qqcw[_idx3(i, k, lc, ncol, pver)] - deltat * dqqcwdt_other[
                                _idx3(i, k, lc, ncol, pver)
                            ]
                            if qqcwold < 0.0:
                                qqcwold = 0.0
                            dryvol_c[_idx3(i, k, n, ncol, pver)] += dum_m2v * qqcwold
                            deldryvol_c[_idx3(i, k, n, ncol, pver)] += (
                                dqqcwdt_other[_idx3(i, k, lc, ncol, pver)]
                                + dqqcwdt[_idx3(i, k, lc, ncol, pver)]
                            ) * dum_m2vdt


@export
def modal_aero_rename_no_acc_crs_xferfracs_codon(
    ncol: int,
    pver: int,
    pcnstxx: int,
    ntot_amode: int,
    maxpair_renamexf: int,
    loffset: int,
    npair_renamexf: int,
    q_p: cobj,
    qqcw_p: cobj,
    dryvol_a_p: cobj,
    dryvol_c_p: cobj,
    deldryvol_a_p: cobj,
    deldryvol_c_p: cobj,
    modefrm_renamexf_p: cobj,
    modetoo_renamexf_p: cobj,
    numptr_amode_p: cobj,
    numptrcw_amode_p: cobj,
    dgnum_amode_p: cobj,
    factoraa_p: cobj,
    factoryy_p: cobj,
    dryvol_smallest_p: cobj,
    v2nlorlx_p: cobj,
    v2nhirlx_p: cobj,
    dum3alnsg2_p: cobj,
    dp_cut_p: cobj,
    lndp_cut_p: cobj,
    dp_belowcut_p: cobj,
    onethird: float,
    xferfrac_max: float,
    xferfrac_vol_p: cobj,
    xferfrac_num_p: cobj,
):
    q = Ptr[float](q_p)
    qqcw = Ptr[float](qqcw_p)
    dryvol_a = Ptr[float](dryvol_a_p)
    dryvol_c = Ptr[float](dryvol_c_p)
    deldryvol_a = Ptr[float](deldryvol_a_p)
    deldryvol_c = Ptr[float](deldryvol_c_p)
    modefrm_renamexf = Ptr[int](modefrm_renamexf_p)
    modetoo_renamexf = Ptr[int](modetoo_renamexf_p)
    numptr_amode = Ptr[int](numptr_amode_p)
    numptrcw_amode = Ptr[int](numptrcw_amode_p)
    dgnum_amode = Ptr[float](dgnum_amode_p)
    factoraa = Ptr[float](factoraa_p)
    factoryy = Ptr[float](factoryy_p)
    dryvol_smallest = Ptr[float](dryvol_smallest_p)
    v2nlorlx = Ptr[float](v2nlorlx_p)
    v2nhirlx = Ptr[float](v2nhirlx_p)
    dum3alnsg2 = Ptr[float](dum3alnsg2_p)
    dp_cut = Ptr[float](dp_cut_p)
    lndp_cut = Ptr[float](lndp_cut_p)
    dp_belowcut = Ptr[float](dp_belowcut_p)
    xferfrac_vol = Ptr[float](xferfrac_vol_p)
    xferfrac_num = Ptr[float](xferfrac_num_p)

    for ipair in range(1, maxpair_renamexf + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                xferfrac_vol[_idx3(i, k, ipair, ncol, pver)] = 0.0
                xferfrac_num[_idx3(i, k, ipair, ncol, pver)] = 0.0

    for ipair in range(1, npair_renamexf + 1):
        mfrm = modefrm_renamexf[ipair - 1]
        mtoo = modetoo_renamexf[ipair - 1]

        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                dryvol_t_old = dryvol_a[_idx3(i, k, mfrm, ncol, pver)] + dryvol_c[
                    _idx3(i, k, mfrm, ncol, pver)
                ]
                dryvol_t_del = deldryvol_a[_idx3(i, k, mfrm, ncol, pver)] + deldryvol_c[
                    _idx3(i, k, mfrm, ncol, pver)
                ]
                dryvol_t_new = dryvol_t_old + dryvol_t_del
                dryvol_t_oldbnd = max(dryvol_t_old, dryvol_smallest[mfrm - 1])

                if dryvol_t_new <= dryvol_smallest[mfrm - 1]:
                    continue
                if dryvol_t_del <= 1.0e-6 * dryvol_t_oldbnd:
                    continue

                num_t_old = q[_idx3(i, k, numptr_amode[mfrm - 1] - loffset, ncol, pver)]
                num_t_old += qqcw[
                    _idx3(i, k, numptrcw_amode[mfrm - 1] - loffset, ncol, pver)
                ]
                if num_t_old < 0.0:
                    num_t_old = 0.0

                dryvol_t_oldbnd = max(dryvol_t_old, dryvol_smallest[mfrm - 1])
                num_t_oldbnd = min(dryvol_t_oldbnd * v2nlorlx[mfrm - 1], num_t_old)
                num_t_oldbnd = max(dryvol_t_oldbnd * v2nhirlx[mfrm - 1], num_t_oldbnd)

                dgn_t_new = (dryvol_t_new / (num_t_oldbnd * factoraa[mfrm - 1])) ** onethird
                if dgn_t_new <= dgnum_amode[mfrm - 1]:
                    continue

                lndgn_new = log(dgn_t_new)
                lndgv_new = lndgn_new + dum3alnsg2[ipair - 1]
                yn_tail = (lndp_cut[ipair - 1] - lndgn_new) * factoryy[mfrm - 1]
                yv_tail = (lndp_cut[ipair - 1] - lndgv_new) * factoryy[mfrm - 1]
                tailfr_numnew = 0.5 * erfc(yn_tail)
                tailfr_volnew = 0.5 * erfc(yv_tail)

                dgn_t_old = (dryvol_t_oldbnd / (num_t_oldbnd * factoraa[mfrm - 1])) ** onethird
                if dgn_t_new >= dp_cut[ipair - 1]:
                    dgn_t_old = min(dgn_t_old, dp_belowcut[ipair - 1])

                lndgn_old = log(dgn_t_old)
                lndgv_old = lndgn_old + dum3alnsg2[ipair - 1]
                yn_tail = (lndp_cut[ipair - 1] - lndgn_old) * factoryy[mfrm - 1]
                yv_tail = (lndp_cut[ipair - 1] - lndgv_old) * factoryy[mfrm - 1]
                tailfr_numold = 0.5 * erfc(yn_tail)
                tailfr_volold = 0.5 * erfc(yv_tail)

                dum = tailfr_volnew * dryvol_t_new - tailfr_volold * dryvol_t_old
                if dum <= 0.0:
                    continue

                xferfrac_vol_val = min(dum, dryvol_t_new) / dryvol_t_new
                xferfrac_vol_val = min(xferfrac_vol_val, xferfrac_max)
                xferfrac_num_val = tailfr_numnew - tailfr_numold
                xferfrac_num_val = max(0.0, min(xferfrac_num_val, xferfrac_vol_val))

                xferfrac_vol[_idx3(i, k, ipair, ncol, pver)] = xferfrac_vol_val
                xferfrac_num[_idx3(i, k, ipair, ncol, pver)] = xferfrac_num_val


@export
def modal_aero_rename_no_acc_crs_tendencies_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnstxx: int,
    maxpair_renamexf: int,
    maxspec_renamexf: int,
    loffset: int,
    npair_renamexf: int,
    jsrflx_rename: int,
    nsrflx: int,
    is_dorename_atik: int,
    deltat: float,
    deltatinv: float,
    gravit: float,
    pdel_p: cobj,
    dorename_atik_p: cobj,
    q_p: cobj,
    qqcw_p: cobj,
    dqdt_p: cobj,
    dqqcwdt_p: cobj,
    qsrflx_p: cobj,
    qqcwsrflx_p: cobj,
    xferfrac_vol_p: cobj,
    xferfrac_num_p: cobj,
    nspecfrm_renamexf_p: cobj,
    lspecfrma_renamexf_p: cobj,
    lspecfrmc_renamexf_p: cobj,
    lspectooa_renamexf_p: cobj,
    lspectooc_renamexf_p: cobj,
):
    pdel = Ptr[float](pdel_p)
    dorename_atik = Ptr[int](dorename_atik_p)
    q = Ptr[float](q_p)
    qqcw = Ptr[float](qqcw_p)
    dqdt = Ptr[float](dqdt_p)
    dqqcwdt = Ptr[float](dqqcwdt_p)
    qsrflx = Ptr[float](qsrflx_p)
    qqcwsrflx = Ptr[float](qqcwsrflx_p)
    xferfrac_vol = Ptr[float](xferfrac_vol_p)
    xferfrac_num = Ptr[float](xferfrac_num_p)
    nspecfrm_renamexf = Ptr[int](nspecfrm_renamexf_p)
    lspecfrma_renamexf = Ptr[int](lspecfrma_renamexf_p)
    lspecfrmc_renamexf = Ptr[int](lspecfrmc_renamexf_p)
    lspectooa_renamexf = Ptr[int](lspectooa_renamexf_p)
    lspectooc_renamexf = Ptr[int](lspectooc_renamexf_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            if is_dorename_atik != 0:
                if dorename_atik[_idx2(i, k, ncol)] == 0:
                    continue

            pdel_fac = pdel[_idx2(i, k, pcols)] / gravit

            for ipair in range(1, npair_renamexf + 1):
                xferfrac_vol_local = xferfrac_vol[_idx3(i, k, ipair, ncol, pver)]
                xferfrac_num_local = xferfrac_num[_idx3(i, k, ipair, ncol, pver)]
                if xferfrac_vol_local <= 0.0:
                    continue

                for iq in range(1, nspecfrm_renamexf[ipair - 1] + 1):
                    xfercoef = xferfrac_vol_local * deltatinv
                    if iq == 1:
                        xfercoef = xferfrac_num_local * deltatinv

                    lsfrma = lspecfrma_renamexf[_idx2(iq, ipair, maxspec_renamexf)] - loffset
                    lsfrmc = lspecfrmc_renamexf[_idx2(iq, ipair, maxspec_renamexf)] - loffset
                    lstooa = lspectooa_renamexf[_idx2(iq, ipair, maxspec_renamexf)] - loffset
                    lstooc = lspectooc_renamexf[_idx2(iq, ipair, maxspec_renamexf)] - loffset

                    if lsfrma > 0:
                        xfertend = xfercoef * max(
                            0.0,
                            q[_idx3(i, k, lsfrma, ncol, pver)]
                            + dqdt[_idx3(i, k, lsfrma, ncol, pver)] * deltat,
                        )
                        dqdt[_idx3(i, k, lsfrma, ncol, pver)] -= xfertend
                        qsrflx[_idx3(i, lsfrma, jsrflx_rename, pcols, pcnstxx)] -= (
                            xfertend * pdel_fac
                        )
                        if lstooa > 0:
                            dqdt[_idx3(i, k, lstooa, ncol, pver)] += xfertend
                            qsrflx[_idx3(i, lstooa, jsrflx_rename, pcols, pcnstxx)] += (
                                xfertend * pdel_fac
                            )

                    if lsfrmc > 0:
                        xfertend = xfercoef * max(
                            0.0,
                            qqcw[_idx3(i, k, lsfrmc, ncol, pver)]
                            + dqqcwdt[_idx3(i, k, lsfrmc, ncol, pver)] * deltat,
                        )
                        dqqcwdt[_idx3(i, k, lsfrmc, ncol, pver)] -= xfertend
                        qqcwsrflx[_idx3(i, lsfrmc, jsrflx_rename, pcols, pcnstxx)] -= (
                            xfertend * pdel_fac
                        )
                        if lstooc > 0:
                            dqqcwdt[_idx3(i, k, lstooc, ncol, pver)] += xfertend
                            qqcwsrflx[_idx3(i, lstooc, jsrflx_rename, pcols, pcnstxx)] += (
                                xfertend * pdel_fac
                            )


@export
def modal_aero_rename_no_acc_crs_sub_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnstxx: int,
    ntot_amode: int,
    maxpair_renamexf: int,
    maxspec_renamexf: int,
    loffset: int,
    npair_renamexf: int,
    is_dorename_atik: int,
    jsrflx_rename: int,
    nsrflx: int,
    deltat: float,
    deltatinv: float,
    onethird: float,
    xferfrac_max: float,
    pi: float,
    gravit: float,
    pdel_p: cobj,
    dorename_atik_p: cobj,
    q_p: cobj,
    qqcw_p: cobj,
    dqdt_p: cobj,
    dqdt_other_p: cobj,
    dqqcwdt_p: cobj,
    dqqcwdt_other_p: cobj,
    qsrflx_p: cobj,
    qqcwsrflx_p: cobj,
    modefrm_renamexf_p: cobj,
    modetoo_renamexf_p: cobj,
    nspec_amode_p: cobj,
    lspectype_amode_p: cobj,
    specmw_amode_p: cobj,
    specdens_amode_p: cobj,
    lmassptr_amode_p: cobj,
    lmassptrcw_amode_p: cobj,
    numptr_amode_p: cobj,
    numptrcw_amode_p: cobj,
    dgnum_amode_p: cobj,
    alnsg_amode_p: cobj,
    voltonumblo_amode_p: cobj,
    voltonumbhi_amode_p: cobj,
    nspecfrm_renamexf_p: cobj,
    lspecfrma_renamexf_p: cobj,
    lspecfrmc_renamexf_p: cobj,
    lspectooa_renamexf_p: cobj,
    lspectooc_renamexf_p: cobj,
    idomode_p: cobj,
    dryvol_a_p: cobj,
    dryvol_c_p: cobj,
    deldryvol_a_p: cobj,
    deldryvol_c_p: cobj,
    factoraa_p: cobj,
    factoryy_p: cobj,
    dryvol_smallest_p: cobj,
    v2nlorlx_p: cobj,
    v2nhirlx_p: cobj,
    dum3alnsg2_p: cobj,
    dp_cut_p: cobj,
    lndp_cut_p: cobj,
    dp_belowcut_p: cobj,
    xferfrac_vol_p: cobj,
    xferfrac_num_p: cobj,
    dotendrn_p: cobj,
    dotendqqcwrn_p: cobj,
):
    idomode = Ptr[int](idomode_p)
    modefrm_renamexf = Ptr[int](modefrm_renamexf_p)
    modetoo_renamexf = Ptr[int](modetoo_renamexf_p)
    nspec_amode = Ptr[int](nspec_amode_p)
    numptr_amode = Ptr[int](numptr_amode_p)
    numptrcw_amode = Ptr[int](numptrcw_amode_p)
    dgnum_amode = Ptr[float](dgnum_amode_p)
    alnsg_amode = Ptr[float](alnsg_amode_p)
    voltonumblo_amode = Ptr[float](voltonumblo_amode_p)
    voltonumbhi_amode = Ptr[float](voltonumbhi_amode_p)
    factoraa = Ptr[float](factoraa_p)
    factoryy = Ptr[float](factoryy_p)
    dryvol_smallest = Ptr[float](dryvol_smallest_p)
    v2nlorlx = Ptr[float](v2nlorlx_p)
    v2nhirlx = Ptr[float](v2nhirlx_p)
    dum3alnsg2 = Ptr[float](dum3alnsg2_p)
    dp_cut = Ptr[float](dp_cut_p)
    lndp_cut = Ptr[float](lndp_cut_p)
    dp_belowcut = Ptr[float](dp_belowcut_p)

    frelax = 27.0

    for n in range(1, ntot_amode + 1):
        idomode[n - 1] = 0

    for ipair in range(1, npair_renamexf + 1):
        mfrm = modefrm_renamexf[ipair - 1]
        mtoo = modetoo_renamexf[ipair - 1]
        idomode[mfrm - 1] = 1

        factoraa[mfrm - 1] = (pi / 6.0) * exp(4.5 * (alnsg_amode[mfrm - 1] ** 2))
        factoraa[mtoo - 1] = (pi / 6.0) * exp(4.5 * (alnsg_amode[mtoo - 1] ** 2))
        factoryy[mfrm - 1] = sqrt(0.5) / alnsg_amode[mfrm - 1]
        dryvol_smallest[mfrm - 1] = 1.0e-25
        v2nlorlx[mfrm - 1] = voltonumblo_amode[mfrm - 1] * frelax
        v2nhirlx[mfrm - 1] = voltonumbhi_amode[mfrm - 1] / frelax

        dum3alnsg2[ipair - 1] = 3.0 * (alnsg_amode[mfrm - 1] ** 2)
        dp_cut[ipair - 1] = sqrt(
            dgnum_amode[mfrm - 1] * exp(1.5 * (alnsg_amode[mfrm - 1] ** 2))
            * dgnum_amode[mtoo - 1]
            * exp(1.5 * (alnsg_amode[mtoo - 1] ** 2))
        )
        lndp_cut[ipair - 1] = log(dp_cut[ipair - 1])
        dp_belowcut[ipair - 1] = 0.99 * dp_cut[ipair - 1]

    modal_aero_rename_no_acc_crs_dryvols_codon(
        ncol,
        pver,
        pcnstxx,
        ntot_amode,
        maxspec_renamexf,
        loffset,
        deltat,
        idomode_p,
        q_p,
        qqcw_p,
        dqdt_p,
        dqdt_other_p,
        dqqcwdt_p,
        dqqcwdt_other_p,
        nspec_amode_p,
        lspectype_amode_p,
        specmw_amode_p,
        specdens_amode_p,
        lmassptr_amode_p,
        lmassptrcw_amode_p,
        dryvol_a_p,
        dryvol_c_p,
        deldryvol_a_p,
        deldryvol_c_p,
    )

    modal_aero_rename_no_acc_crs_xferfracs_codon(
        ncol,
        pver,
        pcnstxx,
        ntot_amode,
        maxpair_renamexf,
        loffset,
        npair_renamexf,
        q_p,
        qqcw_p,
        dryvol_a_p,
        dryvol_c_p,
        deldryvol_a_p,
        deldryvol_c_p,
        modefrm_renamexf_p,
        modetoo_renamexf_p,
        numptr_amode_p,
        numptrcw_amode_p,
        dgnum_amode_p,
        factoraa_p,
        factoryy_p,
        dryvol_smallest_p,
        v2nlorlx_p,
        v2nhirlx_p,
        dum3alnsg2_p,
        dp_cut_p,
        lndp_cut_p,
        dp_belowcut_p,
        onethird,
        xferfrac_max,
        xferfrac_vol_p,
        xferfrac_num_p,
    )

    modal_aero_rename_no_acc_crs_tendencies_codon(
        ncol,
        pcols,
        pver,
        pcnstxx,
        maxpair_renamexf,
        maxspec_renamexf,
        loffset,
        npair_renamexf,
        jsrflx_rename,
        nsrflx,
        is_dorename_atik,
        deltat,
        deltatinv,
        gravit,
        pdel_p,
        dorename_atik_p,
        q_p,
        qqcw_p,
        dqdt_p,
        dqqcwdt_p,
        qsrflx_p,
        qqcwsrflx_p,
        xferfrac_vol_p,
        xferfrac_num_p,
        nspecfrm_renamexf_p,
        lspecfrma_renamexf_p,
        lspecfrmc_renamexf_p,
        lspectooa_renamexf_p,
        lspectooc_renamexf_p,
    )

    modal_aero_rename_set_dotend_flags_codon(
        pcnstxx,
        maxpair_renamexf,
        maxspec_renamexf,
        loffset,
        npair_renamexf,
        nspecfrm_renamexf_p,
        lspecfrma_renamexf_p,
        lspecfrmc_renamexf_p,
        lspectooa_renamexf_p,
        lspectooc_renamexf_p,
        dotendrn_p,
        dotendqqcwrn_p,
    )


@export
def modal_aero_rename_acc_crs_dryvols_codon(
    ncol: int,
    pver: int,
    pcnstxx: int,
    maxspec_renamexf: int,
    loffset: int,
    ixferable_all: int,
    nspec_mfrm: int,
    deltat: float,
    q_p: cobj,
    qqcw_p: cobj,
    dqdt_p: cobj,
    dqdt_other_p: cobj,
    dqqcwdt_p: cobj,
    dqqcwdt_other_p: cobj,
    lspectype_mfrm_p: cobj,
    specmw_amode_p: cobj,
    specdens_amode_p: cobj,
    lmassptr_mfrm_p: cobj,
    lmassptrcw_mfrm_p: cobj,
    ixferable_a_p: cobj,
    ixferable_c_p: cobj,
    dryvol_a_p: cobj,
    dryvol_c_p: cobj,
    deldryvol_a_p: cobj,
    deldryvol_c_p: cobj,
    dryvol_a_xfab_p: cobj,
    dryvol_c_xfab_p: cobj,
):
    q = Ptr[float](q_p)
    qqcw = Ptr[float](qqcw_p)
    dqdt = Ptr[float](dqdt_p)
    dqdt_other = Ptr[float](dqdt_other_p)
    dqqcwdt = Ptr[float](dqqcwdt_p)
    dqqcwdt_other = Ptr[float](dqqcwdt_other_p)
    lspectype_mfrm = Ptr[int](lspectype_mfrm_p)
    specmw_amode = Ptr[float](specmw_amode_p)
    specdens_amode = Ptr[float](specdens_amode_p)
    lmassptr_mfrm = Ptr[int](lmassptr_mfrm_p)
    lmassptrcw_mfrm = Ptr[int](lmassptrcw_mfrm_p)
    ixferable_a = Ptr[int](ixferable_a_p)
    ixferable_c = Ptr[int](ixferable_c_p)
    dryvol_a = Ptr[float](dryvol_a_p)
    dryvol_c = Ptr[float](dryvol_c_p)
    deldryvol_a = Ptr[float](deldryvol_a_p)
    deldryvol_c = Ptr[float](deldryvol_c_p)
    dryvol_a_xfab = Ptr[float](dryvol_a_xfab_p)
    dryvol_c_xfab = Ptr[float](dryvol_c_xfab_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            dryvol_a[_idx2(i, k, ncol)] = 0.0
            dryvol_c[_idx2(i, k, ncol)] = 0.0
            deldryvol_a[_idx2(i, k, ncol)] = 0.0
            deldryvol_c[_idx2(i, k, ncol)] = 0.0
            dryvol_a_xfab[_idx2(i, k, ncol)] = 0.0
            dryvol_c_xfab[_idx2(i, k, ncol)] = 0.0

    for l1 in range(1, nspec_mfrm + 1):
        l2 = lspectype_mfrm[l1 - 1]
        tmp_m2v = specmw_amode[l2 - 1] / specdens_amode[l2 - 1]
        tmp_m2vdt = tmp_m2v * deltat

        la = lmassptr_mfrm[l1 - 1] - loffset
        if la > 0:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    dryvol_a[_idx2(i, k, ncol)] += tmp_m2v * max(
                        0.0,
                        q[_idx3(i, k, la, ncol, pver)]
                        - deltat * dqdt_other[_idx3(i, k, la, ncol, pver)],
                    )
                    deldryvol_a[_idx2(i, k, ncol)] += (
                        dqdt_other[_idx3(i, k, la, ncol, pver)]
                        + dqdt[_idx3(i, k, la, ncol, pver)]
                    ) * tmp_m2vdt
                    if ixferable_all <= 0 and ixferable_a[l1 - 1] > 0:
                        dryvol_a_xfab[_idx2(i, k, ncol)] += tmp_m2v * max(
                            0.0,
                            q[_idx3(i, k, la, ncol, pver)]
                            + deltat * dqdt[_idx3(i, k, la, ncol, pver)],
                        )

        lc = lmassptrcw_mfrm[l1 - 1] - loffset
        if lc > 0:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    dryvol_c[_idx2(i, k, ncol)] += tmp_m2v * max(
                        0.0,
                        qqcw[_idx3(i, k, lc, ncol, pver)]
                        - deltat * dqqcwdt_other[_idx3(i, k, lc, ncol, pver)],
                    )
                    deldryvol_c[_idx2(i, k, ncol)] += (
                        dqqcwdt_other[_idx3(i, k, lc, ncol, pver)]
                        + dqqcwdt[_idx3(i, k, lc, ncol, pver)]
                    ) * tmp_m2vdt
                    if ixferable_all <= 0 and ixferable_c[l1 - 1] > 0:
                        dryvol_c_xfab[_idx2(i, k, ncol)] += tmp_m2v * max(
                            0.0,
                            qqcw[_idx3(i, k, lc, ncol, pver)]
                            + deltat * dqqcwdt[_idx3(i, k, lc, ncol, pver)],
                        )


@export
def modal_aero_rename_acc_crs_tendencies_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnstxx: int,
    maxspec_renamexf: int,
    loffset: int,
    nspecfrm_renamexf: int,
    jsrflx_rename: int,
    nsrflx: int,
    is_dorename_atik: int,
    l_dqdt_rnpos: int,
    deltat: float,
    deltatinv: float,
    gravit: float,
    pdel_p: cobj,
    dorename_atik_p: cobj,
    q_p: cobj,
    qqcw_p: cobj,
    dqdt_p: cobj,
    dqqcwdt_p: cobj,
    qsrflx_p: cobj,
    qqcwsrflx_p: cobj,
    xferfrac_vol_p: cobj,
    xferfrac_num_p: cobj,
    lspecfrma_renamexf_p: cobj,
    lspecfrmc_renamexf_p: cobj,
    lspectooa_renamexf_p: cobj,
    lspectooc_renamexf_p: cobj,
    dqdt_rnpos_p: cobj,
):
    pdel = Ptr[float](pdel_p)
    dorename_atik = Ptr[int](dorename_atik_p)
    q = Ptr[float](q_p)
    qqcw = Ptr[float](qqcw_p)
    dqdt = Ptr[float](dqdt_p)
    dqqcwdt = Ptr[float](dqqcwdt_p)
    qsrflx = Ptr[float](qsrflx_p)
    qqcwsrflx = Ptr[float](qqcwsrflx_p)
    xferfrac_vol = Ptr[float](xferfrac_vol_p)
    xferfrac_num = Ptr[float](xferfrac_num_p)
    lspecfrma_renamexf = Ptr[int](lspecfrma_renamexf_p)
    lspecfrmc_renamexf = Ptr[int](lspecfrmc_renamexf_p)
    lspectooa_renamexf = Ptr[int](lspectooa_renamexf_p)
    lspectooc_renamexf = Ptr[int](lspectooc_renamexf_p)
    dqdt_rnpos = Ptr[float](dqdt_rnpos_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            if is_dorename_atik != 0:
                if dorename_atik[_idx2(i, k, ncol)] == 0:
                    continue

            xferfrac_vol_local = xferfrac_vol[_idx2(i, k, ncol)]
            xferfrac_num_local = xferfrac_num[_idx2(i, k, ncol)]
            if xferfrac_vol_local <= 0.0:
                continue

            pdel_fac = pdel[_idx2(i, k, pcols)] / gravit

            for iq in range(1, nspecfrm_renamexf + 1):
                xfercoef = xferfrac_vol_local * deltatinv
                if iq == 1:
                    xfercoef = xferfrac_num_local * deltatinv

                lsfrma = lspecfrma_renamexf[iq - 1] - loffset
                lsfrmc = lspecfrmc_renamexf[iq - 1] - loffset
                lstooa = lspectooa_renamexf[iq - 1] - loffset
                lstooc = lspectooc_renamexf[iq - 1] - loffset

                if lsfrma > 0:
                    xfertend = xfercoef * max(
                        0.0,
                        q[_idx3(i, k, lsfrma, ncol, pver)]
                        + dqdt[_idx3(i, k, lsfrma, ncol, pver)] * deltat,
                    )
                    dqdt[_idx3(i, k, lsfrma, ncol, pver)] -= xfertend
                    qsrflx[_idx3(i, lsfrma, jsrflx_rename, pcols, pcnstxx)] -= (
                        xfertend * pdel_fac
                    )
                    if lstooa > 0:
                        dqdt[_idx3(i, k, lstooa, ncol, pver)] += xfertend
                        qsrflx[_idx3(i, lstooa, jsrflx_rename, pcols, pcnstxx)] += (
                            xfertend * pdel_fac
                        )
                        if l_dqdt_rnpos != 0:
                            dqdt_rnpos[_idx3(i, k, lstooa, ncol, pver)] += xfertend

                if lsfrmc > 0:
                    xfertend = xfercoef * max(
                        0.0,
                        qqcw[_idx3(i, k, lsfrmc, ncol, pver)]
                        + dqqcwdt[_idx3(i, k, lsfrmc, ncol, pver)] * deltat,
                    )
                    dqqcwdt[_idx3(i, k, lsfrmc, ncol, pver)] -= xfertend
                    qqcwsrflx[_idx3(i, lsfrmc, jsrflx_rename, pcols, pcnstxx)] -= (
                        xfertend * pdel_fac
                    )
                    if lstooc > 0:
                        dqqcwdt[_idx3(i, k, lstooc, ncol, pver)] += xfertend
                        qqcwsrflx[_idx3(i, lstooc, jsrflx_rename, pcols, pcnstxx)] += (
                            xfertend * pdel_fac
                        )


@export
def modal_aero_rename_acc_crs_pair_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnstxx: int,
    maxspec_renamexf: int,
    loffset: int,
    is_dorename_atik: int,
    l_dqdt_rnpos: int,
    jsrflx_rename: int,
    nsrflx: int,
    ixferable_all: int,
    nspec_mfrm: int,
    mfrm: int,
    numptr_amode_mfrm: int,
    numptrcw_amode_mfrm: int,
    igrow_shrink: int,
    method_optbb: int,
    flagaa_shrink: int,
    nspecfrm_renamexf: int,
    deltat: float,
    deltatinv: float,
    gravit: float,
    dgnum_amode_mfrm: float,
    factoraa: float,
    factoryy: float,
    dryvol_smallest: float,
    v2nlorlx: float,
    v2nhirlx: float,
    factor_3alnsg2: float,
    dp_cut: float,
    lndp_cut: float,
    dp_belowcut: float,
    dp_xfernone_thresh: float,
    dp_xferall_thresh: float,
    onethird: float,
    xferfrac_max: float,
    troplev_p: cobj,
    pdel_p: cobj,
    dorename_atik_p: cobj,
    q_p: cobj,
    qqcw_p: cobj,
    dqdt_p: cobj,
    dqdt_other_p: cobj,
    dqqcwdt_p: cobj,
    dqqcwdt_other_p: cobj,
    qsrflx_p: cobj,
    qqcwsrflx_p: cobj,
    lspectype_mfrm_p: cobj,
    specmw_amode_p: cobj,
    specdens_amode_p: cobj,
    lmassptr_mfrm_p: cobj,
    lmassptrcw_mfrm_p: cobj,
    ixferable_a_p: cobj,
    ixferable_c_p: cobj,
    lspecfrma_renamexf_p: cobj,
    lspecfrmc_renamexf_p: cobj,
    lspectooa_renamexf_p: cobj,
    lspectooc_renamexf_p: cobj,
    dryvol_a_p: cobj,
    dryvol_c_p: cobj,
    deldryvol_a_p: cobj,
    deldryvol_c_p: cobj,
    dryvol_a_xfab_p: cobj,
    dryvol_c_xfab_p: cobj,
    xferfrac_vol_p: cobj,
    xferfrac_num_p: cobj,
    dqdt_rnpos_p: cobj,
):
    modal_aero_rename_acc_crs_dryvols_codon(
        ncol,
        pver,
        pcnstxx,
        maxspec_renamexf,
        loffset,
        ixferable_all,
        nspec_mfrm,
        deltat,
        q_p,
        qqcw_p,
        dqdt_p,
        dqdt_other_p,
        dqqcwdt_p,
        dqqcwdt_other_p,
        lspectype_mfrm_p,
        specmw_amode_p,
        specdens_amode_p,
        lmassptr_mfrm_p,
        lmassptrcw_mfrm_p,
        ixferable_a_p,
        ixferable_c_p,
        dryvol_a_p,
        dryvol_c_p,
        deldryvol_a_p,
        deldryvol_c_p,
        dryvol_a_xfab_p,
        dryvol_c_xfab_p,
    )

    modal_aero_rename_acc_crs_xferfracs_codon(
        ncol,
        pcols,
        pver,
        pcnstxx,
        loffset,
        mfrm,
        numptr_amode_mfrm,
        numptrcw_amode_mfrm,
        igrow_shrink,
        ixferable_all,
        method_optbb,
        flagaa_shrink,
        dgnum_amode_mfrm,
        factoraa,
        factoryy,
        dryvol_smallest,
        v2nlorlx,
        v2nhirlx,
        factor_3alnsg2,
        dp_cut,
        lndp_cut,
        dp_belowcut,
        dp_xfernone_thresh,
        dp_xferall_thresh,
        onethird,
        xferfrac_max,
        troplev_p,
        q_p,
        qqcw_p,
        dryvol_a_p,
        dryvol_c_p,
        deldryvol_a_p,
        deldryvol_c_p,
        dryvol_a_xfab_p,
        dryvol_c_xfab_p,
        xferfrac_vol_p,
        xferfrac_num_p,
    )

    modal_aero_rename_acc_crs_tendencies_codon(
        ncol,
        pcols,
        pver,
        pcnstxx,
        maxspec_renamexf,
        loffset,
        nspecfrm_renamexf,
        jsrflx_rename,
        nsrflx,
        is_dorename_atik,
        l_dqdt_rnpos,
        deltat,
        deltatinv,
        gravit,
        pdel_p,
        dorename_atik_p,
        q_p,
        qqcw_p,
        dqdt_p,
        dqqcwdt_p,
        qsrflx_p,
        qqcwsrflx_p,
        xferfrac_vol_p,
        xferfrac_num_p,
        lspecfrma_renamexf_p,
        lspecfrmc_renamexf_p,
        lspectooa_renamexf_p,
        lspectooc_renamexf_p,
        dqdt_rnpos_p,
    )


def modal_aero_rename_acc_crs_dryvols_full_codon(
    ncol: int,
    pver: int,
    maxspec_renamexf: int,
    loffset: int,
    ipair: int,
    mfrm: int,
    nspec_mfrm: int,
    ixferable_all: int,
    deltat: float,
    q_p: cobj,
    qqcw_p: cobj,
    dqdt_p: cobj,
    dqdt_other_p: cobj,
    dqqcwdt_p: cobj,
    dqqcwdt_other_p: cobj,
    lspectype_amode_p: cobj,
    specmw_amode_p: cobj,
    specdens_amode_p: cobj,
    lmassptr_amode_p: cobj,
    lmassptrcw_amode_p: cobj,
    ixferable_a_renamexf_p: cobj,
    ixferable_c_renamexf_p: cobj,
    dryvol_a_p: cobj,
    dryvol_c_p: cobj,
    deldryvol_a_p: cobj,
    deldryvol_c_p: cobj,
    dryvol_a_xfab_p: cobj,
    dryvol_c_xfab_p: cobj,
):
    q = Ptr[float](q_p)
    qqcw = Ptr[float](qqcw_p)
    dqdt = Ptr[float](dqdt_p)
    dqdt_other = Ptr[float](dqdt_other_p)
    dqqcwdt = Ptr[float](dqqcwdt_p)
    dqqcwdt_other = Ptr[float](dqqcwdt_other_p)
    lspectype_amode = Ptr[int](lspectype_amode_p)
    specmw_amode = Ptr[float](specmw_amode_p)
    specdens_amode = Ptr[float](specdens_amode_p)
    lmassptr_amode = Ptr[int](lmassptr_amode_p)
    lmassptrcw_amode = Ptr[int](lmassptrcw_amode_p)
    ixferable_a_renamexf = Ptr[int](ixferable_a_renamexf_p)
    ixferable_c_renamexf = Ptr[int](ixferable_c_renamexf_p)
    dryvol_a = Ptr[float](dryvol_a_p)
    dryvol_c = Ptr[float](dryvol_c_p)
    deldryvol_a = Ptr[float](deldryvol_a_p)
    deldryvol_c = Ptr[float](deldryvol_c_p)
    dryvol_a_xfab = Ptr[float](dryvol_a_xfab_p)
    dryvol_c_xfab = Ptr[float](dryvol_c_xfab_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            dryvol_a[_idx2(i, k, ncol)] = 0.0
            dryvol_c[_idx2(i, k, ncol)] = 0.0
            deldryvol_a[_idx2(i, k, ncol)] = 0.0
            deldryvol_c[_idx2(i, k, ncol)] = 0.0
            dryvol_a_xfab[_idx2(i, k, ncol)] = 0.0
            dryvol_c_xfab[_idx2(i, k, ncol)] = 0.0

    for l1 in range(1, nspec_mfrm + 1):
        l2 = lspectype_amode[_idx2(l1, mfrm, maxspec_renamexf)]
        tmp_m2v = specmw_amode[l2 - 1] / specdens_amode[l2 - 1]
        tmp_m2vdt = tmp_m2v * deltat

        la = lmassptr_amode[_idx2(l1, mfrm, maxspec_renamexf)] - loffset
        if la > 0:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    dryvol_a[_idx2(i, k, ncol)] += tmp_m2v * max(
                        0.0,
                        q[_idx3(i, k, la, ncol, pver)]
                        - deltat * dqdt_other[_idx3(i, k, la, ncol, pver)],
                    )
                    deldryvol_a[_idx2(i, k, ncol)] += (
                        dqdt_other[_idx3(i, k, la, ncol, pver)]
                        + dqdt[_idx3(i, k, la, ncol, pver)]
                    ) * tmp_m2vdt
                    if (
                        ixferable_all <= 0
                        and ixferable_a_renamexf[_idx2(l1, ipair, maxspec_renamexf)] > 0
                    ):
                        dryvol_a_xfab[_idx2(i, k, ncol)] += tmp_m2v * max(
                            0.0,
                            q[_idx3(i, k, la, ncol, pver)]
                            + deltat * dqdt[_idx3(i, k, la, ncol, pver)],
                        )

        lc = lmassptrcw_amode[_idx2(l1, mfrm, maxspec_renamexf)] - loffset
        if lc > 0:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    dryvol_c[_idx2(i, k, ncol)] += tmp_m2v * max(
                        0.0,
                        qqcw[_idx3(i, k, lc, ncol, pver)]
                        - deltat * dqqcwdt_other[_idx3(i, k, lc, ncol, pver)],
                    )
                    deldryvol_c[_idx2(i, k, ncol)] += (
                        dqqcwdt_other[_idx3(i, k, lc, ncol, pver)]
                        + dqqcwdt[_idx3(i, k, lc, ncol, pver)]
                    ) * tmp_m2vdt
                    if (
                        ixferable_all <= 0
                        and ixferable_c_renamexf[_idx2(l1, ipair, maxspec_renamexf)] > 0
                    ):
                        dryvol_c_xfab[_idx2(i, k, ncol)] += tmp_m2v * max(
                            0.0,
                            qqcw[_idx3(i, k, lc, ncol, pver)]
                            + deltat * dqqcwdt[_idx3(i, k, lc, ncol, pver)],
                        )


def modal_aero_rename_acc_crs_tendencies_full_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnstxx: int,
    maxspec_renamexf: int,
    loffset: int,
    ipair: int,
    nspecfrm_ipair: int,
    jsrflx_rename: int,
    is_dorename_atik: int,
    l_dqdt_rnpos: int,
    deltat: float,
    deltatinv: float,
    gravit: float,
    pdel_p: cobj,
    dorename_atik_p: cobj,
    q_p: cobj,
    qqcw_p: cobj,
    dqdt_p: cobj,
    dqqcwdt_p: cobj,
    qsrflx_p: cobj,
    qqcwsrflx_p: cobj,
    xferfrac_vol_p: cobj,
    xferfrac_num_p: cobj,
    lspecfrma_renamexf_p: cobj,
    lspecfrmc_renamexf_p: cobj,
    lspectooa_renamexf_p: cobj,
    lspectooc_renamexf_p: cobj,
    dqdt_rnpos_p: cobj,
):
    pdel = Ptr[float](pdel_p)
    dorename_atik = Ptr[int](dorename_atik_p)
    q = Ptr[float](q_p)
    qqcw = Ptr[float](qqcw_p)
    dqdt = Ptr[float](dqdt_p)
    dqqcwdt = Ptr[float](dqqcwdt_p)
    qsrflx = Ptr[float](qsrflx_p)
    qqcwsrflx = Ptr[float](qqcwsrflx_p)
    xferfrac_vol = Ptr[float](xferfrac_vol_p)
    xferfrac_num = Ptr[float](xferfrac_num_p)
    lspecfrma_renamexf = Ptr[int](lspecfrma_renamexf_p)
    lspecfrmc_renamexf = Ptr[int](lspecfrmc_renamexf_p)
    lspectooa_renamexf = Ptr[int](lspectooa_renamexf_p)
    lspectooc_renamexf = Ptr[int](lspectooc_renamexf_p)
    dqdt_rnpos = Ptr[float](dqdt_rnpos_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            if is_dorename_atik != 0:
                if dorename_atik[_idx2(i, k, ncol)] == 0:
                    continue

            xferfrac_vol_local = xferfrac_vol[_idx2(i, k, ncol)]
            xferfrac_num_local = xferfrac_num[_idx2(i, k, ncol)]
            if xferfrac_vol_local <= 0.0:
                continue

            pdel_fac = pdel[_idx2(i, k, pcols)] / gravit

            for iq in range(1, nspecfrm_ipair + 1):
                xfercoef = xferfrac_vol_local * deltatinv
                if iq == 1:
                    xfercoef = xferfrac_num_local * deltatinv

                lsfrma = lspecfrma_renamexf[_idx2(iq, ipair, maxspec_renamexf)] - loffset
                lsfrmc = lspecfrmc_renamexf[_idx2(iq, ipair, maxspec_renamexf)] - loffset
                lstooa = lspectooa_renamexf[_idx2(iq, ipair, maxspec_renamexf)] - loffset
                lstooc = lspectooc_renamexf[_idx2(iq, ipair, maxspec_renamexf)] - loffset

                if lsfrma > 0:
                    xfertend = xfercoef * max(
                        0.0,
                        q[_idx3(i, k, lsfrma, ncol, pver)]
                        + dqdt[_idx3(i, k, lsfrma, ncol, pver)] * deltat,
                    )
                    dqdt[_idx3(i, k, lsfrma, ncol, pver)] -= xfertend
                    qsrflx[_idx3(i, lsfrma, jsrflx_rename, pcols, pcnstxx)] -= (
                        xfertend * pdel_fac
                    )
                    if lstooa > 0:
                        dqdt[_idx3(i, k, lstooa, ncol, pver)] += xfertend
                        qsrflx[_idx3(i, lstooa, jsrflx_rename, pcols, pcnstxx)] += (
                            xfertend * pdel_fac
                        )
                        if l_dqdt_rnpos != 0:
                            dqdt_rnpos[_idx3(i, k, lstooa, ncol, pver)] += xfertend

                if lsfrmc > 0:
                    xfertend = xfercoef * max(
                        0.0,
                        qqcw[_idx3(i, k, lsfrmc, ncol, pver)]
                        + dqqcwdt[_idx3(i, k, lsfrmc, ncol, pver)] * deltat,
                    )
                    dqqcwdt[_idx3(i, k, lsfrmc, ncol, pver)] -= xfertend
                    qqcwsrflx[_idx3(i, lsfrmc, jsrflx_rename, pcols, pcnstxx)] -= (
                        xfertend * pdel_fac
                    )
                    if lstooc > 0:
                        dqqcwdt[_idx3(i, k, lstooc, ncol, pver)] += xfertend
                        qqcwsrflx[_idx3(i, lstooc, jsrflx_rename, pcols, pcnstxx)] += (
                            xfertend * pdel_fac
                        )


@export
def modal_aero_rename_acc_crs_sub_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnstxx: int,
    maxpair_renamexf: int,
    maxspec_renamexf: int,
    loffset: int,
    npair_renamexf: int,
    is_dorename_atik: int,
    l_dqdt_rnpos: int,
    jsrflx_rename: int,
    nsrflx: int,
    modeptr_coarse: int,
    modeptr_accum: int,
    method_optbb: int,
    deltat: float,
    deltatinv: float,
    onethird: float,
    xferfrac_max: float,
    gravit: float,
    troplev_p: cobj,
    pdel_p: cobj,
    dorename_atik_p: cobj,
    q_p: cobj,
    qqcw_p: cobj,
    dqdt_p: cobj,
    dqdt_other_p: cobj,
    dqqcwdt_p: cobj,
    dqqcwdt_other_p: cobj,
    qsrflx_p: cobj,
    qqcwsrflx_p: cobj,
    modefrm_renamexf_p: cobj,
    modetoo_renamexf_p: cobj,
    nspec_amode_p: cobj,
    lspectype_amode_p: cobj,
    specmw_amode_p: cobj,
    specdens_amode_p: cobj,
    lmassptr_amode_p: cobj,
    lmassptrcw_amode_p: cobj,
    numptr_amode_p: cobj,
    numptrcw_amode_p: cobj,
    dgnum_amode_p: cobj,
    factoraa_p: cobj,
    factoryy_p: cobj,
    dryvol_smallest_p: cobj,
    v2nlorlx_p: cobj,
    v2nhirlx_p: cobj,
    factor_3alnsg2_p: cobj,
    dp_cut_p: cobj,
    lndp_cut_p: cobj,
    dp_belowcut_p: cobj,
    dp_xfernone_threshaa_p: cobj,
    dp_xferall_thresh_p: cobj,
    igrow_shrink_renamexf_p: cobj,
    ixferable_all_renamexf_p: cobj,
    ixferable_a_renamexf_p: cobj,
    ixferable_c_renamexf_p: cobj,
    nspecfrm_renamexf_p: cobj,
    lspecfrma_renamexf_p: cobj,
    lspecfrmc_renamexf_p: cobj,
    lspectooa_renamexf_p: cobj,
    lspectooc_renamexf_p: cobj,
    dryvol_a_p: cobj,
    dryvol_c_p: cobj,
    deldryvol_a_p: cobj,
    deldryvol_c_p: cobj,
    dryvol_a_xfab_p: cobj,
    dryvol_c_xfab_p: cobj,
    xferfrac_vol_p: cobj,
    xferfrac_num_p: cobj,
    dotendrn_p: cobj,
    dotendqqcwrn_p: cobj,
    dqdt_rnpos_p: cobj,
):
    modefrm_renamexf = Ptr[int](modefrm_renamexf_p)
    modetoo_renamexf = Ptr[int](modetoo_renamexf_p)
    nspec_amode = Ptr[int](nspec_amode_p)
    numptr_amode = Ptr[int](numptr_amode_p)
    numptrcw_amode = Ptr[int](numptrcw_amode_p)
    dgnum_amode = Ptr[float](dgnum_amode_p)
    factoraa = Ptr[float](factoraa_p)
    factoryy = Ptr[float](factoryy_p)
    dryvol_smallest = Ptr[float](dryvol_smallest_p)
    v2nlorlx = Ptr[float](v2nlorlx_p)
    v2nhirlx = Ptr[float](v2nhirlx_p)
    factor_3alnsg2 = Ptr[float](factor_3alnsg2_p)
    dp_cut = Ptr[float](dp_cut_p)
    lndp_cut = Ptr[float](lndp_cut_p)
    dp_belowcut = Ptr[float](dp_belowcut_p)
    dp_xfernone_threshaa = Ptr[float](dp_xfernone_threshaa_p)
    dp_xferall_thresh = Ptr[float](dp_xferall_thresh_p)
    igrow_shrink_renamexf = Ptr[int](igrow_shrink_renamexf_p)
    ixferable_all_renamexf = Ptr[int](ixferable_all_renamexf_p)
    nspecfrm_renamexf = Ptr[int](nspecfrm_renamexf_p)

    for ipair in range(1, npair_renamexf + 1):
        mfrm = modefrm_renamexf[ipair - 1]
        mtoo = modetoo_renamexf[ipair - 1]
        nspec_mfrm = nspec_amode[mfrm - 1]
        ixferable_all = ixferable_all_renamexf[ipair - 1]
        igrow_shrink = igrow_shrink_renamexf[ipair - 1]

        flagaa_shrink = 0
        if mfrm == modeptr_coarse and mtoo == modeptr_accum:
            flagaa_shrink = 1

        modal_aero_rename_acc_crs_dryvols_full_codon(
            ncol,
            pver,
            maxspec_renamexf,
            loffset,
            ipair,
            mfrm,
            nspec_mfrm,
            ixferable_all,
            deltat,
            q_p,
            qqcw_p,
            dqdt_p,
            dqdt_other_p,
            dqqcwdt_p,
            dqqcwdt_other_p,
            lspectype_amode_p,
            specmw_amode_p,
            specdens_amode_p,
            lmassptr_amode_p,
            lmassptrcw_amode_p,
            ixferable_a_renamexf_p,
            ixferable_c_renamexf_p,
            dryvol_a_p,
            dryvol_c_p,
            deldryvol_a_p,
            deldryvol_c_p,
            dryvol_a_xfab_p,
            dryvol_c_xfab_p,
        )

        modal_aero_rename_acc_crs_xferfracs_codon(
            ncol,
            pcols,
            pver,
            pcnstxx,
            loffset,
            mfrm,
            numptr_amode[mfrm - 1],
            numptrcw_amode[mfrm - 1],
            igrow_shrink,
            ixferable_all,
            method_optbb,
            flagaa_shrink,
            dgnum_amode[mfrm - 1],
            factoraa[mfrm - 1],
            factoryy[mfrm - 1],
            dryvol_smallest[mfrm - 1],
            v2nlorlx[mfrm - 1],
            v2nhirlx[mfrm - 1],
            factor_3alnsg2[ipair - 1],
            dp_cut[ipair - 1],
            lndp_cut[ipair - 1],
            dp_belowcut[ipair - 1],
            dp_xfernone_threshaa[ipair - 1],
            dp_xferall_thresh[ipair - 1],
            onethird,
            xferfrac_max,
            troplev_p,
            q_p,
            qqcw_p,
            dryvol_a_p,
            dryvol_c_p,
            deldryvol_a_p,
            deldryvol_c_p,
            dryvol_a_xfab_p,
            dryvol_c_xfab_p,
            xferfrac_vol_p,
            xferfrac_num_p,
        )

        modal_aero_rename_acc_crs_tendencies_full_codon(
            ncol,
            pcols,
            pver,
            pcnstxx,
            maxspec_renamexf,
            loffset,
            ipair,
            nspecfrm_renamexf[ipair - 1],
            jsrflx_rename,
            is_dorename_atik,
            l_dqdt_rnpos,
            deltat,
            deltatinv,
            gravit,
            pdel_p,
            dorename_atik_p,
            q_p,
            qqcw_p,
            dqdt_p,
            dqqcwdt_p,
            qsrflx_p,
            qqcwsrflx_p,
            xferfrac_vol_p,
            xferfrac_num_p,
            lspecfrma_renamexf_p,
            lspecfrmc_renamexf_p,
            lspectooa_renamexf_p,
            lspectooc_renamexf_p,
            dqdt_rnpos_p,
        )

    modal_aero_rename_set_dotend_flags_codon(
        pcnstxx,
        maxpair_renamexf,
        maxspec_renamexf,
        loffset,
        npair_renamexf,
        nspecfrm_renamexf_p,
        lspecfrma_renamexf_p,
        lspecfrmc_renamexf_p,
        lspectooa_renamexf_p,
        lspectooc_renamexf_p,
        dotendrn_p,
        dotendqqcwrn_p,
    )


@export
def modal_aero_rename_set_dotend_flags_codon(
    pcnstxx: int,
    maxpair_renamexf: int,
    maxspec_renamexf: int,
    loffset: int,
    npair_renamexf: int,
    nspecfrm_renamexf_p: cobj,
    lspecfrma_renamexf_p: cobj,
    lspecfrmc_renamexf_p: cobj,
    lspectooa_renamexf_p: cobj,
    lspectooc_renamexf_p: cobj,
    dotendrn_p: cobj,
    dotendqqcwrn_p: cobj,
):
    nspecfrm_renamexf = Ptr[int](nspecfrm_renamexf_p)
    lspecfrma_renamexf = Ptr[int](lspecfrma_renamexf_p)
    lspecfrmc_renamexf = Ptr[int](lspecfrmc_renamexf_p)
    lspectooa_renamexf = Ptr[int](lspectooa_renamexf_p)
    lspectooc_renamexf = Ptr[int](lspectooc_renamexf_p)
    dotendrn = Ptr[int](dotendrn_p)
    dotendqqcwrn = Ptr[int](dotendqqcwrn_p)

    for l in range(1, pcnstxx + 1):
        dotendrn[l - 1] = 0
        dotendqqcwrn[l - 1] = 0

    for ipair in range(1, npair_renamexf + 1):
        for iq in range(1, nspecfrm_renamexf[ipair - 1] + 1):
            lsfrma = lspecfrma_renamexf[_idx2(iq, ipair, maxspec_renamexf)] - loffset
            lsfrmc = lspecfrmc_renamexf[_idx2(iq, ipair, maxspec_renamexf)] - loffset
            lstooa = lspectooa_renamexf[_idx2(iq, ipair, maxspec_renamexf)] - loffset
            lstooc = lspectooc_renamexf[_idx2(iq, ipair, maxspec_renamexf)] - loffset

            if lsfrma > 0:
                dotendrn[lsfrma - 1] = 1
                if lstooa > 0:
                    dotendrn[lstooa - 1] = 1

            if lsfrmc > 0:
                dotendqqcwrn[lsfrmc - 1] = 1
                if lstooc > 0:
                    dotendqqcwrn[lstooc - 1] = 1


@export
def modal_aero_rename_acc_crs_xferfracs_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnstxx: int,
    loffset: int,
    mfrm: int,
    numptr_amode_mfrm: int,
    numptrcw_amode_mfrm: int,
    igrow_shrink: int,
    ixferable_all: int,
    method_optbb: int,
    flagaa_shrink: int,
    dgnum_amode_mfrm: float,
    factoraa: float,
    factoryy: float,
    dryvol_smallest: float,
    v2nlorlx: float,
    v2nhirlx: float,
    factor_3alnsg2: float,
    dp_cut: float,
    lndp_cut: float,
    dp_belowcut: float,
    dp_xfernone_thresh: float,
    dp_xferall_thresh: float,
    onethird: float,
    xferfrac_max: float,
    troplev_p: cobj,
    q_p: cobj,
    qqcw_p: cobj,
    dryvol_a_p: cobj,
    dryvol_c_p: cobj,
    deldryvol_a_p: cobj,
    deldryvol_c_p: cobj,
    dryvol_a_xfab_p: cobj,
    dryvol_c_xfab_p: cobj,
    xferfrac_vol_p: cobj,
    xferfrac_num_p: cobj,
):
    troplev = Ptr[int](troplev_p)
    q = Ptr[float](q_p)
    qqcw = Ptr[float](qqcw_p)
    dryvol_a = Ptr[float](dryvol_a_p)
    dryvol_c = Ptr[float](dryvol_c_p)
    deldryvol_a = Ptr[float](deldryvol_a_p)
    deldryvol_c = Ptr[float](deldryvol_c_p)
    dryvol_a_xfab = Ptr[float](dryvol_a_xfab_p)
    dryvol_c_xfab = Ptr[float](dryvol_c_xfab_p)
    xferfrac_vol = Ptr[float](xferfrac_vol_p)
    xferfrac_num = Ptr[float](xferfrac_num_p)

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            xferfrac_vol[_idx2(i, k, ncol)] = 0.0
            xferfrac_num[_idx2(i, k, ncol)] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            dryvol_t_old = dryvol_a[_idx2(i, k, ncol)] + dryvol_c[_idx2(i, k, ncol)]
            dryvol_t_del = deldryvol_a[_idx2(i, k, ncol)] + deldryvol_c[_idx2(i, k, ncol)]
            dryvol_t_new = dryvol_t_old + dryvol_t_del
            dryvol_t_oldbnd = max(dryvol_t_old, dryvol_smallest)

            if igrow_shrink > 0:
                if dryvol_t_new <= dryvol_smallest:
                    continue
                if method_optbb != 2:
                    if dryvol_t_del <= 1.0e-6 * dryvol_t_oldbnd:
                        continue

                num_t_old = q[_idx3(i, k, numptr_amode_mfrm - loffset, ncol, pver)]
                num_t_old += qqcw[
                    _idx3(i, k, numptrcw_amode_mfrm - loffset, ncol, pver)
                ]
                num_t_old = max(0.0, num_t_old)
                dryvol_t_oldbnd = max(dryvol_t_old, dryvol_smallest)
                num_t_oldbnd = min(dryvol_t_oldbnd * v2nlorlx, num_t_old)
                num_t_oldbnd = max(dryvol_t_oldbnd * v2nhirlx, num_t_oldbnd)

                dgn_t_new = (dryvol_t_new / (num_t_oldbnd * factoraa)) ** onethird
                if dgn_t_new <= dp_xfernone_thresh:
                    continue

                dgn_t_old = (dryvol_t_oldbnd / (num_t_oldbnd * factoraa)) ** onethird
                dgn_t_oldb = dgn_t_old
                dryvol_t_oldb = dryvol_t_old
                if method_optbb == 2:
                    if dgn_t_old >= dp_cut:
                        dryvol_t_oldb = dryvol_t_old * (dp_belowcut / dgn_t_old) ** 3
                        dgn_t_oldb = dp_belowcut
                    if dgn_t_new < dp_xferall_thresh:
                        if (dryvol_t_new - dryvol_t_oldb) <= 1.0e-6 * dryvol_t_oldbnd:
                            continue
                elif dgn_t_new >= dp_cut:
                    dgn_t_oldb = min(dgn_t_oldb, dp_belowcut)

                lndgn_new = log(dgn_t_new)
                lndgv_new = lndgn_new + factor_3alnsg2
                yn_tail = (lndp_cut - lndgn_new) * factoryy
                yv_tail = (lndp_cut - lndgv_new) * factoryy
                tailfr_numnew = 0.5 * erfc(yn_tail)
                tailfr_volnew = 0.5 * erfc(yv_tail)

                lndgn_old = log(dgn_t_oldb)
                lndgv_old = lndgn_old + factor_3alnsg2
                yn_tail = (lndp_cut - lndgn_old) * factoryy
                yv_tail = (lndp_cut - lndgv_old) * factoryy
                tailfr_numold = 0.5 * erfc(yn_tail)
                tailfr_volold = 0.5 * erfc(yv_tail)

                if method_optbb == 2 and dgn_t_new >= dp_xferall_thresh:
                    dryvol_xferamt = dryvol_t_new
                else:
                    dryvol_xferamt = (
                        tailfr_volnew * dryvol_t_new - tailfr_volold * dryvol_t_oldb
                    )
                if dryvol_xferamt <= 0.0:
                    continue

                xferfrac_vol_val = max(0.0, dryvol_xferamt / dryvol_t_new)
                if method_optbb == 2 and xferfrac_vol_val >= xferfrac_max:
                    xferfrac_vol_val = 1.0
                    xferfrac_num_val = 1.0
                else:
                    xferfrac_vol_val = min(xferfrac_vol_val, xferfrac_max)
                    xferfrac_num_val = tailfr_numnew - tailfr_numold
                    xferfrac_num_val = max(
                        0.0, min(xferfrac_num_val, xferfrac_vol_val)
                    )

                if ixferable_all <= 0:
                    dryvol_t_new_xfab = max(
                        0.0,
                        dryvol_a_xfab[_idx2(i, k, ncol)] + dryvol_c_xfab[_idx2(i, k, ncol)],
                    )
                    dryvol_xferamt = xferfrac_vol_val * dryvol_t_new
                    if dryvol_t_new_xfab >= 0.999999 * dryvol_xferamt:
                        xferfrac_vol_val = min(1.0, dryvol_xferamt / dryvol_t_new_xfab)
                    elif dryvol_t_new_xfab >= 1.0e-7 * dryvol_xferamt:
                        xferfrac_vol_val = 1.0
                        xferfrac_num_val = xferfrac_num_val * (
                            dryvol_t_new_xfab / dryvol_xferamt
                        )
                    else:
                        continue

            else:
                if dryvol_t_old <= dryvol_smallest:
                    continue

                if dryvol_t_del >= -1.0e-6 * dryvol_t_oldbnd:
                    if flagaa_shrink != 0 and k < troplev[i - 1]:
                        flagbb_shrink = 1
                    else:
                        continue
                else:
                    flagbb_shrink = 0

                num_t_old = q[_idx3(i, k, numptr_amode_mfrm - loffset, ncol, pver)]
                num_t_old += qqcw[
                    _idx3(i, k, numptrcw_amode_mfrm - loffset, ncol, pver)
                ]
                num_t_old = max(0.0, num_t_old)
                dryvol_t_oldbnd = max(dryvol_t_old, dryvol_smallest)
                num_t_oldbnd = min(dryvol_t_oldbnd * v2nlorlx, num_t_old)
                num_t_oldbnd = max(dryvol_t_oldbnd * v2nhirlx, num_t_oldbnd)

                dgn_t_new = (dryvol_t_new / (num_t_oldbnd * factoraa)) ** onethird
                if dgn_t_new >= dp_xfernone_thresh:
                    continue
                if flagbb_shrink != 0:
                    if dgn_t_new > dp_cut:
                        continue

                if dgn_t_new <= dp_xferall_thresh:
                    tailfr_numnew = 1.0
                    tailfr_volnew = 1.0
                else:
                    lndgn_new = log(dgn_t_new)
                    lndgv_new = lndgn_new + factor_3alnsg2
                    yn_tail = (lndp_cut - lndgn_new) * factoryy
                    yv_tail = (lndp_cut - lndgv_new) * factoryy
                    tailfr_numnew = 1.0 - 0.5 * erfc(yn_tail)
                    tailfr_volnew = 1.0 - 0.5 * erfc(yv_tail)

                dgn_t_old = (dryvol_t_oldbnd / (num_t_oldbnd * factoraa)) ** onethird
                dgn_t_oldb = dgn_t_old
                dryvol_t_oldb = dryvol_t_old
                tailfr_numold = 0.0
                tailfr_volold = 0.0

                xferfrac_vol_val = tailfr_volnew
                if xferfrac_vol_val <= 0.0:
                    continue
                xferfrac_num_val = tailfr_numnew

                if xferfrac_vol_val >= xferfrac_max:
                    xferfrac_vol_val = 1.0
                    xferfrac_num_val = 1.0
                else:
                    xferfrac_vol_val = min(xferfrac_vol_val, xferfrac_max)
                    xferfrac_num_val = max(xferfrac_num_val, xferfrac_vol_val)
                    xferfrac_num_val = min(xferfrac_max, xferfrac_num_val)

                if ixferable_all <= 0:
                    dryvol_t_new_xfab = max(
                        0.0,
                        dryvol_a_xfab[_idx2(i, k, ncol)] + dryvol_c_xfab[_idx2(i, k, ncol)],
                    )
                    dryvol_xferamt = xferfrac_vol_val * dryvol_t_new
                    if dryvol_t_new_xfab >= 0.999999 * dryvol_xferamt:
                        xferfrac_vol_val = min(1.0, dryvol_xferamt / dryvol_t_new_xfab)
                    elif dryvol_t_new_xfab >= 1.0e-7 * dryvol_xferamt:
                        xferfrac_vol_val = 1.0
                        xferfrac_num_val = xferfrac_num_val * (
                            dryvol_t_new_xfab / dryvol_xferamt
                        )
                    else:
                        continue

            xferfrac_vol[_idx2(i, k, ncol)] = xferfrac_vol_val
            xferfrac_num[_idx2(i, k, ncol)] = xferfrac_num_val


@export
def modal_aero_newnuc_apply_tendencies_codon(
    i_c: int,
    k_c: int,
    ncol_c: int,
    pcols_c: int,
    pver_c: int,
    pcnst_c: int,
    pdel_p: cobj,
    dqdt_p: cobj,
    qsrflx_p: cobj,
    q_p: cobj,
    gravit_c: float,
    cldx_c: float,
    deltat_c: float,
    dso4dt_ait_c: float,
    dndt_ait_c: float,
    dnh4dt_ait_c: float,
    l_h2so4_c: int,
    lso4ait_c: int,
    lnumait_c: int,
    l_nh3_c: int,
    lnh4ait_c: int,
    do_nh3_c: int,
):
    pdel = Ptr[float](pdel_p)
    dqdt = Ptr[float](dqdt_p)
    qsrflx = Ptr[float](qsrflx_p)
    q = Ptr[float](q_p)

    i = i_c
    k = k_c
    ncol = ncol_c
    pcols = pcols_c
    pver = pver_c
    pcnst = pcnst_c

    pdel_fac = pdel[_idx2(i, k, pcols)] / gravit_c

    dqdt[_idx3(i, k, l_h2so4_c, ncol, pver)] = -dso4dt_ait_c * (1.0 - cldx_c)
    qsrflx[_idx3(i, l_h2so4_c, 1, pcols, pcnst)] = (
        qsrflx[_idx3(i, l_h2so4_c, 1, pcols, pcnst)]
        + dqdt[_idx3(i, k, l_h2so4_c, ncol, pver)] * pdel_fac
    )
    q[_idx3(i, k, l_h2so4_c, ncol, pver)] = (
        q[_idx3(i, k, l_h2so4_c, ncol, pver)]
        + dqdt[_idx3(i, k, l_h2so4_c, ncol, pver)] * deltat_c
    )

    dqdt[_idx3(i, k, lso4ait_c, ncol, pver)] = dso4dt_ait_c * (1.0 - cldx_c)
    qsrflx[_idx3(i, lso4ait_c, 1, pcols, pcnst)] = (
        qsrflx[_idx3(i, lso4ait_c, 1, pcols, pcnst)]
        + dqdt[_idx3(i, k, lso4ait_c, ncol, pver)] * pdel_fac
    )
    q[_idx3(i, k, lso4ait_c, ncol, pver)] = (
        q[_idx3(i, k, lso4ait_c, ncol, pver)]
        + dqdt[_idx3(i, k, lso4ait_c, ncol, pver)] * deltat_c
    )

    if lnumait_c > 0:
        dqdt[_idx3(i, k, lnumait_c, ncol, pver)] = dndt_ait_c * (1.0 - cldx_c)
        qsrflx[_idx3(i, lnumait_c, 1, pcols, pcnst)] = (
            qsrflx[_idx3(i, lnumait_c, 1, pcols, pcnst)]
            + dqdt[_idx3(i, k, lnumait_c, ncol, pver)] * pdel_fac
        )
        q[_idx3(i, k, lnumait_c, ncol, pver)] = (
            q[_idx3(i, k, lnumait_c, ncol, pver)]
            + dqdt[_idx3(i, k, lnumait_c, ncol, pver)] * deltat_c
        )

    if do_nh3_c != 0 and dnh4dt_ait_c > 0.0:
        dqdt[_idx3(i, k, l_nh3_c, ncol, pver)] = -dnh4dt_ait_c * (1.0 - cldx_c)
        qsrflx[_idx3(i, l_nh3_c, 1, pcols, pcnst)] = (
            qsrflx[_idx3(i, l_nh3_c, 1, pcols, pcnst)]
            + dqdt[_idx3(i, k, l_nh3_c, ncol, pver)] * pdel_fac
        )
        q[_idx3(i, k, l_nh3_c, ncol, pver)] = (
            q[_idx3(i, k, l_nh3_c, ncol, pver)]
            + dqdt[_idx3(i, k, l_nh3_c, ncol, pver)] * deltat_c
        )

        dqdt[_idx3(i, k, lnh4ait_c, ncol, pver)] = dnh4dt_ait_c * (1.0 - cldx_c)
        qsrflx[_idx3(i, lnh4ait_c, 1, pcols, pcnst)] = (
            qsrflx[_idx3(i, lnh4ait_c, 1, pcols, pcnst)]
            + dqdt[_idx3(i, k, lnh4ait_c, ncol, pver)] * pdel_fac
        )
        q[_idx3(i, k, lnh4ait_c, ncol, pver)] = (
            q[_idx3(i, k, lnh4ait_c, ncol, pver)]
            + dqdt[_idx3(i, k, lnh4ait_c, ncol, pver)] * deltat_c
        )


@export
def modal_aero_newnuc_postprocess_label_codon(
    postprocess_code_c: int,
    tmpch1_code_p: cobj,
):
    tmpch1_code = Ptr[int](tmpch1_code_p)

    tmpch1_code[0] = 32
    if postprocess_code_c == 1:
        tmpch1_code[0] = 65
    elif postprocess_code_c == 2:
        tmpch1_code[0] = 66
    elif postprocess_code_c == 3:
        tmpch1_code[0] = 67
    elif postprocess_code_c == 4:
        tmpch1_code[0] = 69


@export
def modal_aero_newnuc_set_tendency_flags_codon(
    pcnst_c: int,
    lnumait_c: int,
    lso4ait_c: int,
    l_h2so4_c: int,
    lnh4ait_c: int,
    l_nh3_c: int,
    do_nh3_flag_c: int,
    dotend_p: cobj,
    do_nh3_p: cobj,
):
    dotend = Ptr[int](dotend_p)
    do_nh3 = Ptr[int](do_nh3_p)

    for l in range(1, pcnst_c + 1):
        dotend[l - 1] = 0

    dotend[lnumait_c - 1] = 1
    dotend[lso4ait_c - 1] = 1
    dotend[l_h2so4_c - 1] = 1

    if do_nh3_flag_c != 0:
        do_nh3[0] = 1
        dotend[lnh4ait_c - 1] = 1
        dotend[l_nh3_c - 1] = 1
    else:
        do_nh3[0] = 0


def _modal_aero_newnuc_setup_modes_core(
    loffset_c: int,
    pcnst_c: int,
    l_h2so4_sv_c: int,
    l_nh3_sv_c: int,
    lnumait_sv_c: int,
    lso4ait_sv_c: int,
    lptr_nh4_aitken_c: int,
    dgnumlo_aitken_c: float,
    dgnum_aitken_c: float,
    dgnumhi_aitken_c: float,
    specdens_so4_amode_c: float,
    pi_c: float,
):
    l_h2so4 = l_h2so4_sv_c - loffset_c
    l_nh3 = l_nh3_sv_c - loffset_c
    lnumait = lnumait_sv_c - loffset_c
    lso4ait = lso4ait_sv_c - loffset_c
    lnh4ait = lptr_nh4_aitken_c - loffset_c

    dplom_mode_1 = 0.0
    dphim_mode_1 = 0.0
    mass1p_aitlo = 0.0
    mass1p_aithi = 0.0
    do_nh3_flag = 0
    valid_mask = 0

    if l_h2so4 <= 0 or lso4ait <= 0 or lnumait <= 0:
        return (
            l_h2so4,
            l_nh3,
            lnumait,
            lnh4ait,
            lso4ait,
            do_nh3_flag,
            valid_mask,
            dplom_mode_1,
            dphim_mode_1,
            mass1p_aitlo,
            mass1p_aithi,
        )

    valid_mask = 1

    if l_nh3 > 0 and l_nh3 <= pcnst_c and lnh4ait > 0 and lnh4ait <= pcnst_c:
        do_nh3_flag = 1

    dplom_mode_1 = exp(0.67 * log(dgnumlo_aitken_c) + 0.33 * log(dgnum_aitken_c))
    dphim_mode_1 = dgnumhi_aitken_c

    tmpa = specdens_so4_amode_c * pi_c / 6.0
    mass1p_aitlo = tmpa * (dplom_mode_1**3)
    mass1p_aithi = tmpa * (dphim_mode_1**3)
    return (
        l_h2so4,
        l_nh3,
        lnumait,
        lnh4ait,
        lso4ait,
        do_nh3_flag,
        valid_mask,
        dplom_mode_1,
        dphim_mode_1,
        mass1p_aitlo,
        mass1p_aithi,
    )


@export
def modal_aero_newnuc_setup_modes_codon(
    loffset_c: int,
    pcnst_c: int,
    l_h2so4_sv_c: int,
    l_nh3_sv_c: int,
    lnumait_sv_c: int,
    lso4ait_sv_c: int,
    lptr_nh4_aitken_c: int,
    dgnumlo_aitken_c: float,
    dgnum_aitken_c: float,
    dgnumhi_aitken_c: float,
    specdens_so4_amode_c: float,
    pi_c: float,
    l_h2so4_p: cobj,
    l_nh3_p: cobj,
    lnumait_p: cobj,
    lnh4ait_p: cobj,
    lso4ait_p: cobj,
    do_nh3_flag_p: cobj,
    valid_mask_p: cobj,
    dplom_mode_1_p: cobj,
    dphim_mode_1_p: cobj,
    mass1p_aitlo_p: cobj,
    mass1p_aithi_p: cobj,
):
    l_h2so4 = Ptr[int](l_h2so4_p)
    l_nh3 = Ptr[int](l_nh3_p)
    lnumait = Ptr[int](lnumait_p)
    lnh4ait = Ptr[int](lnh4ait_p)
    lso4ait = Ptr[int](lso4ait_p)
    do_nh3_flag = Ptr[int](do_nh3_flag_p)
    valid_mask = Ptr[int](valid_mask_p)
    dplom_mode_1 = Ptr[float](dplom_mode_1_p)
    dphim_mode_1 = Ptr[float](dphim_mode_1_p)
    mass1p_aitlo = Ptr[float](mass1p_aitlo_p)
    mass1p_aithi = Ptr[float](mass1p_aithi_p)

    (
        l_h2so4[0],
        l_nh3[0],
        lnumait[0],
        lnh4ait[0],
        lso4ait[0],
        do_nh3_flag[0],
        valid_mask[0],
        dplom_mode_1[0],
        dphim_mode_1[0],
        mass1p_aitlo[0],
        mass1p_aithi[0],
    ) = _modal_aero_newnuc_setup_modes_core(
        loffset_c,
        pcnst_c,
        l_h2so4_sv_c,
        l_nh3_sv_c,
        lnumait_sv_c,
        lso4ait_sv_c,
        lptr_nh4_aitken_c,
        dgnumlo_aitken_c,
        dgnum_aitken_c,
        dgnumhi_aitken_c,
        specdens_so4_amode_c,
        pi_c,
    )


@export
def modal_aero_newnuc_scale_qsrflx_codon(
    ncol_c: int,
    pcols_c: int,
    pcnst_c: int,
    lmz_c: int,
    adv_mass_c: float,
    mwdry_c: float,
    qsrflx_p: cobj,
):
    qsrflx = Ptr[float](qsrflx_p)

    for i in range(1, ncol_c + 1):
        qsrflx[_idx3(i, lmz_c, 1, pcols_c, pcnst_c)] = (
            qsrflx[_idx3(i, lmz_c, 1, pcols_c, pcnst_c)] * (adv_mass_c / mwdry_c)
        )


@export
def modal_aero_newnuc_sub_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    pcnstxx: int,
    top_lev: int,
    loffset: int,
    deltat: float,
    qh2so4_cutoff: float,
    l_h2so4_sv: int,
    l_nh3_sv: int,
    lnumait_sv: int,
    lso4ait_sv: int,
    lptr_nh4_aitken: int,
    dgnumlo_aitken: float,
    dgnum_aitken: float,
    dgnumhi_aitken: float,
    specdens_so4_amode: float,
    specmw_so4_amode: float,
    specmw_nh4_amode: float,
    pi_c: float,
    gravit_c: float,
    mwdry_c: float,
    adv_mass_p: cobj,
    rgas: float,
    avogad: float,
    mw_so4a: float,
    mw_nh4a: float,
    t_p: cobj,
    pmid_p: cobj,
    pdel_p: cobj,
    zm_p: cobj,
    pblh_p: cobj,
    qv_p: cobj,
    cld_p: cobj,
    q_p: cobj,
    qv_sat_p: cobj,
    del_h2so4_gasprod_p: cobj,
    del_h2so4_aeruptk_p: cobj,
    dqdt_p: cobj,
    qsrflx_p: cobj,
    dplom_mode_p: cobj,
    dphim_mode_p: cobj,
    active_mask_p: cobj,
    cldx_p: cobj,
    qh2so4_cur_p: cobj,
    qh2so4_avg_p: cobj,
    qnh3_cur_p: cobj,
    tmp_uptkrate_p: cobj,
    relhumnn_p: cobj,
    dotend_p: cobj,
    fallback_required_p: cobj,
    ternary_codon_used_p: cobj,
):
    t = Ptr[float](t_p)
    pmid = Ptr[float](pmid_p)
    zm = Ptr[float](zm_p)
    pblh = Ptr[float](pblh_p)
    active_mask = Ptr[int](active_mask_p)
    cldx_work = Ptr[float](cldx_p)
    qh2so4_cur_work = Ptr[float](qh2so4_cur_p)
    qh2so4_avg_work = Ptr[float](qh2so4_avg_p)
    qnh3_cur_work = Ptr[float](qnh3_cur_p)
    tmp_uptkrate_work = Ptr[float](tmp_uptkrate_p)
    relhumnn_work = Ptr[float](relhumnn_p)
    dotend = Ptr[int](dotend_p)
    fallback_required = Ptr[int](fallback_required_p)
    ternary_codon_used = Ptr[int](ternary_codon_used_p)
    dplom_mode = Ptr[float](dplom_mode_p)
    dphim_mode = Ptr[float](dphim_mode_p)
    adv_mass = Ptr[float](adv_mass_p)

    fallback_required[0] = 0
    ternary_codon_used[0] = 0
    for l in range(1, pcnst + 1):
        dotend[l - 1] = 0

    (
        l_h2so4,
        l_nh3,
        lnumait,
        lnh4ait,
        lso4ait,
        do_nh3_flag,
        valid_mask,
        dplom_mode_1,
        dphim_mode_1,
        mass1p_aitlo,
        mass1p_aithi,
    ) = _modal_aero_newnuc_setup_modes_core(
        loffset,
        pcnst,
        l_h2so4_sv,
        l_nh3_sv,
        lnumait_sv,
        lso4ait_sv,
        lptr_nh4_aitken,
        dgnumlo_aitken,
        dgnum_aitken,
        dgnumhi_aitken,
        specdens_so4_amode,
        pi_c,
    )

    dplom_mode[0] = dplom_mode_1
    dphim_mode[0] = dphim_mode_1

    if valid_mask == 0:
        return

    dotend[lnumait - 1] = 1
    dotend[lso4ait - 1] = 1
    dotend[l_h2so4 - 1] = 1

    do_nh3 = 0
    if do_nh3_flag != 0:
        do_nh3 = 1
        dotend[lnh4ait - 1] = 1
        dotend[l_nh3 - 1] = 1

    modal_aero_newnuc_zero_tendencies_codon(
        ncol,
        pcols,
        pver,
        pcnstxx,
        pcnst,
        1,
        dqdt_p,
        qsrflx_p,
    )

    modal_aero_newnuc_prepare_box_inputs_codon(
        ncol,
        pcols,
        pver,
        top_lev,
        l_h2so4,
        l_nh3,
        do_nh3,
        deltat,
        qh2so4_cutoff,
        q_p,
        qv_p,
        cld_p,
        qv_sat_p,
        del_h2so4_gasprod_p,
        del_h2so4_aeruptk_p,
        active_mask_p,
        cldx_p,
        qh2so4_cur_p,
        qh2so4_avg_p,
        qnh3_cur_p,
        tmp_uptkrate_p,
        relhumnn_p,
    )

    newnuc_method_flagaa = 11

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, ncol)
            if active_mask[idx2] == 0:
                continue

            if newnuc_method_flagaa != 2:
                nh3ppt = qnh3_cur_work[idx2] * 1.0e12
                if nh3ppt >= 0.1:
                    cair = pmid[_idx2(i, k, pcols)] / (t[_idx2(i, k, pcols)] * rgas)
                    so4vol_in = qh2so4_avg_work[idx2] * cair * avogad * 1.0e-6
                    if so4vol_in >= 5.0e4:
                        ternary_codon_used[0] = 1

            (
                _fallback_required,
                _isize_nuc,
                qnuma_del,
                qso4a_del,
                qnh4a_del,
                _qh2so4_del,
                _qnh3_del,
                _dens_nh4so4a,
            ) = _mer07_veh02_nuc_mosaic_1box_core(
                newnuc_method_flagaa,
                deltat,
                t[_idx2(i, k, pcols)],
                relhumnn_work[idx2],
                pmid[_idx2(i, k, pcols)],
                zm[_idx2(i, k, pcols)],
                pblh[i - 1],
                qh2so4_cur_work[idx2],
                qh2so4_avg_work[idx2],
                qnh3_cur_work[idx2],
                tmp_uptkrate_work[idx2],
                specmw_so4_amode,
                1,
                dplom_mode,
                dphim_mode,
                -1,
                pi_c,
                rgas,
                avogad,
                mw_so4a,
                mw_nh4a,
            )
            if _fallback_required != 0:
                fallback_required[0] = 1
                return

            (
                _qnuma_del_work,
                dndt_ait,
                _dmdt_ait,
                dso4dt_ait,
                dnh4dt_ait,
                _dndt_aitsv1,
                _dmdt_aitsv1,
                _dndt_aitsv2,
                _dmdt_aitsv2,
                _dndt_aitsv3,
                _dmdt_aitsv3,
                _postprocess_code,
            ) = _mer07_veh02_nuc_mosaic_postprocess_core(
                qnuma_del,
                qso4a_del,
                qnh4a_del,
                deltat,
                specmw_so4_amode,
                specmw_nh4_amode,
                mass1p_aitlo,
                mass1p_aithi,
            )

            modal_aero_newnuc_apply_tendencies_codon(
                i,
                k,
                ncol,
                pcols,
                pver,
                pcnst,
                pdel_p,
                dqdt_p,
                qsrflx_p,
                q_p,
                gravit_c,
                cldx_work[idx2],
                deltat,
                dso4dt_ait,
                dndt_ait,
                dnh4dt_ait,
                l_h2so4,
                lso4ait,
                lnumait,
                l_nh3,
                lnh4ait,
                do_nh3,
            )

    for lmz in range(1, pcnst + 1):
        if dotend[lmz - 1] == 0:
            continue
        modal_aero_newnuc_scale_qsrflx_codon(
            ncol,
            pcols,
            pcnst,
            lmz,
            adv_mass[lmz - 1],
            mwdry_c,
            qsrflx_p,
        )


@export
def modal_aero_coag_sub_codon(
    stage: int,
    ncol: int,
    pcols: int,
    pver: int,
    pcnstxx: int,
    pcnst: int,
    top_lev: int,
    ntot_amode: int,
    maxd_aspectype: int,
    maxpair_acoag: int,
    maxspec_acoag: int,
    pair_option_acoag: int,
    npair_acoag: int,
    ip_aitacc: int,
    ip_pcaacc: int,
    ip_aitpca: int,
    macc: int,
    mait: int,
    mpca: int,
    deltat: float,
    deltatinv_main: float,
    xferfrac_max: float,
    dr_so4_monolayers_pcage: float,
    fac_volsfc_pcarbon: float,
    r_universal: float,
    gravit: float,
    mwdry: float,
    q_p: cobj,
    dqdt_p: cobj,
    qsrflx_p: cobj,
    t_p: cobj,
    pmid_p: cobj,
    pdel_p: cobj,
    dgncur_a_p: cobj,
    ybetaij0_p: cobj,
    ybetaij3_p: cobj,
    ybetaii0_p: cobj,
    ybetajj0_p: cobj,
    xnumbconc_p: cobj,
    xnumbconcavg_p: cobj,
    xnumbconcnew_p: cobj,
    iselfcoagdone_p: cobj,
    modefrm_acoag_p: cobj,
    modetoo_acoag_p: cobj,
    nspecfrm_acoag_p: cobj,
    lspecfrm_acoag_p: cobj,
    lspectoo_acoag_p: cobj,
    mprognum_amode_p: cobj,
    numptr_amode_p: cobj,
    nspec_amode_p: cobj,
    lmassptr_amode_p: cobj,
    lptr_so4_a_amode_p: cobj,
    lptr_nh4_a_amode_p: cobj,
    lptr_soa_a_amode_p: cobj,
    idomode_p: cobj,
    fac_m2v_aitage_p: cobj,
    fac_m2v_pcarbon_p: cobj,
    adv_mass_p: cobj,
    dotend_p: cobj,
):
    q = Ptr[float](q_p)
    dqdt = Ptr[float](dqdt_p)
    qsrflx = Ptr[float](qsrflx_p)
    t = Ptr[float](t_p)
    pmid = Ptr[float](pmid_p)
    pdel = Ptr[float](pdel_p)
    dgncur_a = Ptr[float](dgncur_a_p)
    ybetaij0 = Ptr[float](ybetaij0_p)
    ybetaij3 = Ptr[float](ybetaij3_p)
    ybetaii0 = Ptr[float](ybetaii0_p)
    ybetajj0 = Ptr[float](ybetajj0_p)
    xnumbconc = Ptr[float](xnumbconc_p)
    xnumbconcavg = Ptr[float](xnumbconcavg_p)
    xnumbconcnew = Ptr[float](xnumbconcnew_p)
    iselfcoagdone = Ptr[int](iselfcoagdone_p)
    modefrm_acoag = Ptr[int](modefrm_acoag_p)
    modetoo_acoag = Ptr[int](modetoo_acoag_p)
    nspecfrm_acoag = Ptr[int](nspecfrm_acoag_p)
    lspecfrm_acoag = Ptr[int](lspecfrm_acoag_p)
    lspectoo_acoag = Ptr[int](lspectoo_acoag_p)
    mprognum_amode = Ptr[int](mprognum_amode_p)
    numptr_amode = Ptr[int](numptr_amode_p)
    nspec_amode = Ptr[int](nspec_amode_p)
    lmassptr_amode = Ptr[int](lmassptr_amode_p)
    lptr_so4_a_amode = Ptr[int](lptr_so4_a_amode_p)
    lptr_nh4_a_amode = Ptr[int](lptr_nh4_a_amode_p)
    lptr_soa_a_amode = Ptr[int](lptr_soa_a_amode_p)
    idomode = Ptr[int](idomode_p)
    fac_m2v_aitage = Ptr[float](fac_m2v_aitage_p)
    fac_m2v_pcarbon = Ptr[float](fac_m2v_pcarbon_p)
    adv_mass = Ptr[float](adv_mass_p)
    dotend = Ptr[int](dotend_p)

    if stage == 1:
        lmz = 1
        while lmz <= pcnstxx:
            dotend[lmz - 1] = 0
            i = 1
            while i <= pcols:
                qsrflx[_idx2(i, lmz, pcols)] = 0.0
                i += 1
            k = 1
            while k <= pver:
                i = 1
                while i <= ncol:
                    dqdt[_idx3(i, k, lmz, ncol, pver)] = 0.0
                    i += 1
                k += 1
            lmz += 1
        return

    if stage == 3:
        lmz = 1
        while lmz <= pcnstxx:
            dotend[lmz - 1] = 0
            lmz += 1

        ipair = 1
        while ipair <= npair_acoag:
            modefrm = modefrm_acoag[ipair - 1]
            modetoo = modetoo_acoag[ipair - 1]

            iq = 1
            while iq <= nspecfrm_acoag[ipair - 1]:
                lsfrm = lspecfrm_acoag[_idx2(iq, ipair, maxspec_acoag)]
                lstoo = lspectoo_acoag[_idx2(iq, ipair, maxspec_acoag)]
                if lsfrm > 0 and lsfrm <= pcnstxx:
                    dotend[lsfrm - 1] = 1
                if lstoo > 0 and lstoo <= pcnstxx:
                    dotend[lstoo - 1] = 1
                iq += 1

            if mprognum_amode[modefrm - 1] > 0:
                lsfrm = numptr_amode[modefrm - 1]
                if lsfrm > 0 and lsfrm <= pcnstxx:
                    dotend[lsfrm - 1] = 1
            if mprognum_amode[modetoo - 1] > 0:
                lstoo = numptr_amode[modetoo - 1]
                if lstoo > 0 and lstoo <= pcnstxx:
                    dotend[lstoo - 1] = 1
            ipair += 1

        lmz = 1
        while lmz <= pcnstxx:
            if dotend[lmz - 1] == 0:
                lmz += 1
                continue
            i = 1
            while i <= pcols:
                qsrflx[_idx2(i, lmz, pcols)] = 0.0
                i += 1
            k = top_lev
            while k <= pver:
                i = 1
                while i <= ncol:
                    qsrflx[_idx2(i, lmz, pcols)] = (
                        qsrflx[_idx2(i, lmz, pcols)]
                        + dqdt[_idx3(i, k, lmz, ncol, pver)] * pdel[_idx2(i, k, pcols)]
                    )
                    i += 1
                k += 1
            scale = adv_mass[lmz - 1] / (gravit * mwdry)
            i = 1
            while i <= pcols:
                qsrflx[_idx2(i, lmz, pcols)] = qsrflx[_idx2(i, lmz, pcols)] * scale
                i += 1
            lmz += 1
        return

    if stage != 2:
        return

    k = top_lev
    while k <= pver:
        i = 1
        while i <= ncol:
            aircon = pmid[_idx2(i, k, pcols)] / (r_universal * t[_idx2(i, k, pcols)])

            n = 1
            while n <= ntot_amode:
                idx_mode = _idx3(i, k, n, pcols, pver)
                if idomode[n - 1] > 0:
                    lmz = numptr_amode[n - 1]
                    xnumbconc[idx_mode] = q[_idx3(i, k, lmz, ncol, pver)] * aircon
                    xnumbconc[idx_mode] = max(0.0, xnumbconc[idx_mode])
                else:
                    xnumbconc[idx_mode] = 0.0
                xnumbconcavg[idx_mode] = 0.0
                xnumbconcnew[idx_mode] = 0.0
                iselfcoagdone[idx_mode] = 0
                n += 1

            if pair_option_acoag == 1 or pair_option_acoag == 2:
                ipair = 1
                while ipair <= npair_acoag:
                    modefrm = modefrm_acoag[ipair - 1]
                    modetoo = modetoo_acoag[ipair - 1]
                    idx_pair = _idx3(i, k, ipair, pcols, pver)
                    idx_frm = _idx3(i, k, modefrm, pcols, pver)
                    idx_too = _idx3(i, k, modetoo, pcols, pver)

                    if mprognum_amode[modetoo - 1] > 0 and iselfcoagdone[idx_too] <= 0:
                        iselfcoagdone[idx_too] = 1
                        tmpn = xnumbconc[idx_too]
                        xnumbconcnew[idx_too] = tmpn / (1.0 + deltat * ybetajj0[idx_pair] * tmpn)
                        xnumbconcavg[idx_too] = 0.5 * (xnumbconcnew[idx_too] + tmpn)
                        lstoo = numptr_amode[modetoo - 1]
                        q[_idx3(i, k, lstoo, ncol, pver)] = xnumbconcnew[idx_too] / aircon
                        dqdt[_idx3(i, k, lstoo, ncol, pver)] = (
                            (xnumbconcnew[idx_too] - tmpn) * deltatinv_main / aircon
                        )

                    if mprognum_amode[modefrm - 1] > 0 and iselfcoagdone[idx_frm] <= 0:
                        iselfcoagdone[idx_frm] = 1
                        tmpn = xnumbconc[idx_frm]
                        tmpa = deltat * ybetaij0[idx_pair] * xnumbconcavg[idx_too]
                        tmpb = deltat * ybetaii0[idx_pair]
                        tmpc = tmpa + tmpb * tmpn
                        if abs(tmpc) < 0.01:
                            xnumbconcnew[idx_frm] = tmpn * exp(-tmpc)
                        elif abs(tmpa) < 0.001:
                            xnumbconcnew[idx_frm] = exp(-tmpa) * tmpn / (1.0 + tmpb * tmpn)
                        else:
                            tmpf = tmpb * tmpn / tmpc
                            tmpg = exp(-tmpa)
                            tmph = tmpg * (1.0 - tmpf) / (1.0 - tmpg * tmpf)
                            xnumbconcnew[idx_frm] = tmpn * max(0.0, min(1.0, tmph))
                        xnumbconcavg[idx_frm] = 0.5 * (xnumbconcnew[idx_frm] + tmpn)
                        lsfrm = numptr_amode[modefrm - 1]
                        q[_idx3(i, k, lsfrm, ncol, pver)] = xnumbconcnew[idx_frm] / aircon
                        dqdt[_idx3(i, k, lsfrm, ncol, pver)] = (
                            (xnumbconcnew[idx_frm] - tmpn) * deltatinv_main / aircon
                        )

                    dumloss = ybetaij3[idx_pair] * xnumbconcavg[idx_too]
                    xferfracvol = 1.0 - exp(-dumloss * deltat)
                    xferfracvol = max(0.0, min(xferfrac_max, xferfracvol))

                    iq = 1
                    while iq <= nspecfrm_acoag[ipair - 1]:
                        lsfrm = lspecfrm_acoag[_idx2(iq, ipair, maxspec_acoag)]
                        lstoo = lspectoo_acoag[_idx2(iq, ipair, maxspec_acoag)]
                        if lsfrm > 0:
                            idx_qfrm = _idx3(i, k, lsfrm, ncol, pver)
                            xferamt = q[idx_qfrm] * xferfracvol
                            dqdt[idx_qfrm] = dqdt[idx_qfrm] - xferamt * deltatinv_main
                            q[idx_qfrm] = q[idx_qfrm] - xferamt
                            if lstoo > 0:
                                idx_qtoo = _idx3(i, k, lstoo, ncol, pver)
                                dqdt[idx_qtoo] = dqdt[idx_qtoo] + xferamt * deltatinv_main
                                q[idx_qtoo] = q[idx_qtoo] + xferamt
                        iq += 1
                    ipair += 1

            elif pair_option_acoag == 3:
                idx_macc = _idx3(i, k, macc, pcols, pver)
                idx_mpca = _idx3(i, k, mpca, pcols, pver)
                idx_mait = _idx3(i, k, mait, pcols, pver)
                idx_aitacc = _idx3(i, k, ip_aitacc, pcols, pver)
                idx_pcaacc = _idx3(i, k, ip_pcaacc, pcols, pver)
                idx_aitpca = _idx3(i, k, ip_aitpca, pcols, pver)

                if mprognum_amode[macc - 1] > 0:
                    tmpn = xnumbconc[idx_macc]
                    xnumbconcnew[idx_macc] = tmpn / (1.0 + deltat * ybetajj0[idx_aitacc] * tmpn)
                    xnumbconcavg[idx_macc] = 0.5 * (xnumbconcnew[idx_macc] + tmpn)
                    lstoo = numptr_amode[macc - 1]
                    q[_idx3(i, k, lstoo, ncol, pver)] = xnumbconcnew[idx_macc] / aircon
                    dqdt[_idx3(i, k, lstoo, ncol, pver)] = (
                        (xnumbconcnew[idx_macc] - tmpn) * deltatinv_main / aircon
                    )

                if mprognum_amode[mpca - 1] > 0:
                    tmpn = xnumbconc[idx_mpca]
                    tmpa = deltat * ybetaij0[idx_pcaacc] * xnumbconcavg[idx_macc]
                    tmpb = deltat * ybetaii0[idx_pcaacc]
                    tmpc = tmpa + tmpb * tmpn
                    if abs(tmpc) < 0.01:
                        xnumbconcnew[idx_mpca] = tmpn * exp(-tmpc)
                    elif abs(tmpa) < 0.001:
                        xnumbconcnew[idx_mpca] = exp(-tmpa) * tmpn / (1.0 + tmpb * tmpn)
                    else:
                        tmpf = tmpb * tmpn / tmpc
                        tmpg = exp(-tmpa)
                        tmph = tmpg * (1.0 - tmpf) / (1.0 - tmpg * tmpf)
                        xnumbconcnew[idx_mpca] = tmpn * max(0.0, min(1.0, tmph))
                    xnumbconcavg[idx_mpca] = 0.5 * (xnumbconcnew[idx_mpca] + tmpn)
                    lsfrm = numptr_amode[mpca - 1]
                    q[_idx3(i, k, lsfrm, ncol, pver)] = xnumbconcnew[idx_mpca] / aircon
                    dqdt[_idx3(i, k, lsfrm, ncol, pver)] = (
                        (xnumbconcnew[idx_mpca] - tmpn) * deltatinv_main / aircon
                    )

                if mprognum_amode[mait - 1] > 0:
                    tmpn = xnumbconc[idx_mait]
                    tmpa = deltat * (
                        ybetaij0[idx_aitacc] * xnumbconcavg[idx_macc]
                        + ybetaij0[idx_aitpca] * xnumbconcavg[idx_mpca]
                    )
                    tmpb = deltat * ybetaii0[idx_aitacc]
                    tmpc = tmpa + tmpb * tmpn
                    if abs(tmpc) < 0.01:
                        xnumbconcnew[idx_mait] = tmpn * exp(-tmpc)
                    elif abs(tmpa) < 0.001:
                        xnumbconcnew[idx_mait] = exp(-tmpa) * tmpn / (1.0 + tmpb * tmpn)
                    else:
                        tmpf = tmpb * tmpn / tmpc
                        tmpg = exp(-tmpa)
                        tmph = tmpg * (1.0 - tmpf) / (1.0 - tmpg * tmpf)
                        xnumbconcnew[idx_mait] = tmpn * max(0.0, min(1.0, tmph))
                    xnumbconcavg[idx_mait] = 0.5 * (xnumbconcnew[idx_mait] + tmpn)
                    lsfrm = numptr_amode[mait - 1]
                    q[_idx3(i, k, lsfrm, ncol, pver)] = xnumbconcnew[idx_mait] / aircon
                    dqdt[_idx3(i, k, lsfrm, ncol, pver)] = (
                        (xnumbconcnew[idx_mait] - tmpn) * deltatinv_main / aircon
                    )

                dumloss = (
                    ybetaij3[idx_aitacc] * xnumbconcavg[idx_macc]
                    + ybetaij3[idx_aitpca] * xnumbconcavg[idx_mpca]
                )
                tmpa = ybetaij3[idx_aitpca] * xnumbconcavg[idx_mpca] / max(dumloss, 1.0e-37)
                xferfracvol = 1.0 - exp(-dumloss * deltat)
                xferfracvol = max(0.0, min(xferfrac_max, xferfracvol))
                vol_shell = 0.0

                ipair = ip_aitacc
                iq = 1
                while iq <= nspecfrm_acoag[ipair - 1]:
                    lsfrm = lspecfrm_acoag[_idx2(iq, ipair, maxspec_acoag)]
                    lstoo = lspectoo_acoag[_idx2(iq, ipair, maxspec_acoag)]
                    if lsfrm > 0:
                        idx_qfrm = _idx3(i, k, lsfrm, ncol, pver)
                        xferamt = q[idx_qfrm] * xferfracvol
                        dqdt[idx_qfrm] = dqdt[idx_qfrm] - xferamt * deltatinv_main
                        q[idx_qfrm] = q[idx_qfrm] - xferamt
                        if lstoo > 0:
                            idx_qtoo = _idx3(i, k, lstoo, ncol, pver)
                            dqdt[idx_qtoo] = dqdt[idx_qtoo] + xferamt * deltatinv_main
                            q[idx_qtoo] = q[idx_qtoo] + xferamt
                        vol_shell = vol_shell + xferamt * tmpa * fac_m2v_aitage[iq - 1]
                    iq += 1

                vol_core = 0.0
                l = 1
                while l <= nspec_amode[mpca - 1]:
                    lmz = lmassptr_amode[_idx2(l, mpca, maxd_aspectype)]
                    vol_core = vol_core + q[_idx3(i, k, lmz, ncol, pver)] * fac_m2v_pcarbon[l - 1]
                    l += 1

                tmp1 = vol_shell * dgncur_a[_idx3(i, k, mpca, pcols, pver)] * fac_volsfc_pcarbon
                tmp2 = 6.0 * dr_so4_monolayers_pcage * vol_core
                tmp2 = max(tmp2, 0.0)
                xferfrac_pcage = 0.0
                if tmp1 >= tmp2:
                    xferfrac_pcage = xferfrac_max
                else:
                    xferfrac_pcage = min(tmp1 / tmp2, xferfrac_max)

                dumloss = ybetaij3[idx_pcaacc] * xnumbconcavg[idx_macc]
                xferfracvol = 1.0 - exp(-dumloss * deltat)
                xferfracvol = xferfracvol + xferfrac_pcage
                xferfracvol = max(0.0, min(xferfrac_max, xferfracvol))

                ipair = ip_pcaacc
                iq = 1
                while iq <= nspecfrm_acoag[ipair - 1]:
                    lsfrm = lspecfrm_acoag[_idx2(iq, ipair, maxspec_acoag)]
                    lstoo = lspectoo_acoag[_idx2(iq, ipair, maxspec_acoag)]
                    if lsfrm > 0:
                        idx_qfrm = _idx3(i, k, lsfrm, ncol, pver)
                        xferamt = q[idx_qfrm] * xferfracvol
                        dqdt[idx_qfrm] = dqdt[idx_qfrm] - xferamt * deltatinv_main
                        q[idx_qfrm] = q[idx_qfrm] - xferamt
                        if lstoo > 0:
                            idx_qtoo = _idx3(i, k, lstoo, ncol, pver)
                            dqdt[idx_qtoo] = dqdt[idx_qtoo] + xferamt * deltatinv_main
                            q[idx_qtoo] = q[idx_qtoo] + xferamt
                    iq += 1

                lsfrm = numptr_amode[mpca - 1]
                lstoo = numptr_amode[macc - 1]
                if lsfrm > 0:
                    idx_qfrm = _idx3(i, k, lsfrm, ncol, pver)
                    xferamt = q[idx_qfrm] * xferfrac_pcage
                    dqdt[idx_qfrm] = dqdt[idx_qfrm] - xferamt * deltatinv_main
                    q[idx_qfrm] = q[idx_qfrm] - xferamt
                    if lstoo > 0:
                        idx_qtoo = _idx3(i, k, lstoo, ncol, pver)
                        dqdt[idx_qtoo] = dqdt[idx_qtoo] + xferamt * deltatinv_main
                        q[idx_qtoo] = q[idx_qtoo] + xferamt
            i += 1
        k += 1


def _modal_aero_calcsize_zero_state(
    pcols: int,
    pver: int,
    pcnst: int,
    dqqcwdt: Ptr[float],
    qsrflx: Ptr[float],
    dotend: Ptr[int],
    dotendqqcw: Ptr[int],
):
    l = 1
    while l <= pcnst:
        dotend[l - 1] = 0
        dotendqqcw[l - 1] = 0
        k = 1
        while k <= pver:
            i = 1
            while i <= pcols:
                dqqcwdt[_idx3(i, k, l, pcols, pver)] = 0.0
                i += 1
            k += 1
        jsrflx = 1
        while jsrflx <= 4:
            jac = 1
            while jac <= 2:
                i = 1
                while i <= pcols:
                    qsrflx[_idx4(i, l, jsrflx, jac, pcols, pcnst, 4)] = 0.0
                    i += 1
                jac += 1
            jsrflx += 1
        l += 1


@export
def modal_aero_calcsize_sub_codon(
    ncol: int,
    pcols: int,
    pver: int,
    pcnst: int,
    top_lev: int,
    ntot_amode: int,
    maxd_aspectype: int,
    maxspec_renamexf: int,
    nait: int,
    nacc: int,
    do_adjust: int,
    do_aitacc_transfer: int,
    nspecfrm_pair1: int,
    deltat: float,
    gravit: float,
    q_p: cobj,
    qqcw_p: cobj,
    pdel_p: cobj,
    dqdt_p: cobj,
    dqqcwdt_p: cobj,
    qsrflx_p: cobj,
    dgncur_a_p: cobj,
    dryvol_a_p: cobj,
    dryvol_c_p: cobj,
    drv_a_aitsv_p: cobj,
    num_a_aitsv_p: cobj,
    drv_c_aitsv_p: cobj,
    num_c_aitsv_p: cobj,
    drv_a_accsv_p: cobj,
    num_a_accsv_p: cobj,
    drv_c_accsv_p: cobj,
    num_c_accsv_p: cobj,
    dotend_p: cobj,
    dotendqqcw_p: cobj,
    mprognum_amode_p: cobj,
    numptr_amode_p: cobj,
    numptrcw_amode_p: cobj,
    nspec_amode_p: cobj,
    lspectype_amode_p: cobj,
    lmassptr_amode_p: cobj,
    lmassptrcw_amode_p: cobj,
    dgnum_amode_p: cobj,
    dgnumhi_amode_p: cobj,
    dgnumlo_amode_p: cobj,
    alnsg_amode_p: cobj,
    voltonumb_amode_p: cobj,
    voltonumblo_amode_p: cobj,
    voltonumbhi_amode_p: cobj,
    specdens_amode_p: cobj,
    lspecfrma_pair1_p: cobj,
    lspecfrmc_pair1_p: cobj,
    lspectooa_pair1_p: cobj,
    lspectooc_pair1_p: cobj,
):
    q = Ptr[float](q_p)
    qqcw = Ptr[float](qqcw_p)
    pdel = Ptr[float](pdel_p)
    dqdt = Ptr[float](dqdt_p)
    dqqcwdt = Ptr[float](dqqcwdt_p)
    qsrflx = Ptr[float](qsrflx_p)
    dgncur_a = Ptr[float](dgncur_a_p)
    dryvol_a = Ptr[float](dryvol_a_p)
    dryvol_c = Ptr[float](dryvol_c_p)
    drv_a_aitsv = Ptr[float](drv_a_aitsv_p)
    num_a_aitsv = Ptr[float](num_a_aitsv_p)
    drv_c_aitsv = Ptr[float](drv_c_aitsv_p)
    num_c_aitsv = Ptr[float](num_c_aitsv_p)
    drv_a_accsv = Ptr[float](drv_a_accsv_p)
    num_a_accsv = Ptr[float](num_a_accsv_p)
    drv_c_accsv = Ptr[float](drv_c_accsv_p)
    num_c_accsv = Ptr[float](num_c_accsv_p)
    dotend = Ptr[int](dotend_p)
    dotendqqcw = Ptr[int](dotendqqcw_p)
    mprognum_amode = Ptr[int](mprognum_amode_p)
    numptr_amode = Ptr[int](numptr_amode_p)
    numptrcw_amode = Ptr[int](numptrcw_amode_p)
    nspec_amode = Ptr[int](nspec_amode_p)
    lspectype_amode = Ptr[int](lspectype_amode_p)
    lmassptr_amode = Ptr[int](lmassptr_amode_p)
    lmassptrcw_amode = Ptr[int](lmassptrcw_amode_p)
    dgnum_amode = Ptr[float](dgnum_amode_p)
    dgnumhi_amode = Ptr[float](dgnumhi_amode_p)
    dgnumlo_amode = Ptr[float](dgnumlo_amode_p)
    alnsg_amode = Ptr[float](alnsg_amode_p)
    voltonumb_amode = Ptr[float](voltonumb_amode_p)
    voltonumblo_amode = Ptr[float](voltonumblo_amode_p)
    voltonumbhi_amode = Ptr[float](voltonumbhi_amode_p)
    specdens_amode = Ptr[float](specdens_amode_p)
    lspecfrma_pair1 = Ptr[int](lspecfrma_pair1_p)
    lspecfrmc_pair1 = Ptr[int](lspecfrmc_pair1_p)
    lspectooa_pair1 = Ptr[int](lspectooa_pair1_p)
    lspectooc_pair1 = Ptr[int](lspectooc_pair1_p)

    third = 1.0 / 3.0

    _modal_aero_calcsize_zero_state(pcols, pver, pcnst, dqqcwdt, qsrflx, dotend, dotendqqcw)

    deltatinv = 1.0 / (deltat * (1.0 + 1.0e-15))
    tadj = deltat
    tadj = 86400.0
    tadj = max(tadj, deltat)
    tadjinv = 1.0 / (tadj * (1.0 + 1.0e-15))
    fracadj = deltat * tadjinv
    fracadj = max(0.0, min(1.0, fracadj))
    dumfac = 0.0

    n = 1
    while n <= ntot_amode:
        k = top_lev
        while k <= pver:
            i = 1
            while i <= ncol:
                dgncur_a[_idx3(i, k, n, pcols, pver)] = dgnum_amode[n - 1]
                dryvol_a[_idx2(i, k, pcols)] = 0.0
                dryvol_c[_idx2(i, k, pcols)] = 0.0
                i += 1
            k += 1

        l1 = 1
        while l1 <= nspec_amode[n - 1]:
            lsptype = lspectype_amode[_idx2(l1, n, maxd_aspectype)]
            dummwdens = 1.0 / specdens_amode[lsptype - 1]
            la = lmassptr_amode[_idx2(l1, n, maxd_aspectype)]
            lc = lmassptrcw_amode[_idx2(l1, n, maxd_aspectype)]
            k = top_lev
            while k <= pver:
                i = 1
                while i <= ncol:
                    dryvol_a[_idx2(i, k, pcols)] = (
                        dryvol_a[_idx2(i, k, pcols)]
                        + max(0.0, q[_idx3(i, k, la, pcols, pver)]) * dummwdens
                    )
                    dryvol_c[_idx2(i, k, pcols)] = (
                        dryvol_c[_idx2(i, k, pcols)]
                        + max(0.0, qqcw[_idx3(i, k, lc, pcols, pver)]) * dummwdens
                    )
                    i += 1
                k += 1
            l1 += 1

        lna = numptr_amode[n - 1]
        lnc = numptrcw_amode[n - 1]

        if mprognum_amode[n - 1] <= 0:
            if lna > 0:
                dotend[lna - 1] = 1
                k = top_lev
                while k <= pver:
                    i = 1
                    while i <= ncol:
                        dqdt[_idx3(i, k, lna, pcols, pver)] = (
                            dryvol_a[_idx2(i, k, pcols)] * voltonumb_amode[n - 1]
                            - q[_idx3(i, k, lna, pcols, pver)]
                        ) * deltatinv
                        i += 1
                    k += 1
            if lnc > 0:
                dotendqqcw[lnc - 1] = 1
                k = top_lev
                while k <= pver:
                    i = 1
                    while i <= ncol:
                        dqqcwdt[_idx3(i, k, lnc, pcols, pver)] = (
                            dryvol_c[_idx2(i, k, pcols)] * voltonumb_amode[n - 1]
                            - qqcw[_idx3(i, k, lnc, pcols, pver)]
                        ) * deltatinv
                        i += 1
                    k += 1

        frelaxadj = 27.0
        dumfac = exp(4.5 * (alnsg_amode[n - 1] ** 2.0)) * pi / 6.0
        v2nxx = voltonumbhi_amode[n - 1]
        v2nyy = voltonumblo_amode[n - 1]
        v2nxxrl = v2nxx / frelaxadj
        v2nyyrl = v2nyy * frelaxadj
        dgnxx = dgnumhi_amode[n - 1]
        dgnyy = dgnumlo_amode[n - 1]

        if do_aitacc_transfer != 0:
            if n == nait:
                v2nxx = v2nxx / 1.0e6
            if n == nacc:
                v2nyy = v2nyy * 1.0e6
            v2nxxrl = v2nxx / frelaxadj
            v2nyyrl = v2nyy * frelaxadj

        if do_adjust != 0:
            dotend[lna - 1] = 1
            dotendqqcw[lnc - 1] = 1

        k = top_lev
        while k <= pver:
            i = 1
            while i <= ncol:
                drv_a = dryvol_a[_idx2(i, k, pcols)]
                num_a0 = q[_idx3(i, k, lna, pcols, pver)]
                num_a = max(0.0, num_a0)
                drv_c = dryvol_c[_idx2(i, k, pcols)]
                num_c0 = qqcw[_idx3(i, k, lnc, pcols, pver)]
                num_c = max(0.0, num_c0)

                if do_adjust != 0:
                    if (drv_a <= 0.0) and (drv_c <= 0.0):
                        num_a = 0.0
                        dqdt[_idx3(i, k, lna, pcols, pver)] = -num_a0 * deltatinv
                        num_c = 0.0
                        dqqcwdt[_idx3(i, k, lnc, pcols, pver)] = -num_c0 * deltatinv
                    elif drv_c <= 0.0:
                        num_c = 0.0
                        dqqcwdt[_idx3(i, k, lnc, pcols, pver)] = -num_c0 * deltatinv
                        num_a1 = num_a
                        numbnd = max(drv_a * v2nxx, min(drv_a * v2nyy, num_a1))
                        num_a = num_a1 + (numbnd - num_a1) * fracadj
                        dqdt[_idx3(i, k, lna, pcols, pver)] = (num_a - num_a0) * deltatinv
                    elif drv_a <= 0.0:
                        num_a = 0.0
                        dqdt[_idx3(i, k, lna, pcols, pver)] = -num_a0 * deltatinv
                        num_c1 = num_c
                        numbnd = max(drv_c * v2nxx, min(drv_c * v2nyy, num_c1))
                        num_c = num_c1 + (numbnd - num_c1) * fracadj
                        dqqcwdt[_idx3(i, k, lnc, pcols, pver)] = (num_c - num_c0) * deltatinv
                    else:
                        num_a1 = num_a
                        num_c1 = num_c
                        numbnd = max(drv_a * v2nxxrl, min(drv_a * v2nyyrl, num_a1))
                        delnum_a2 = (numbnd - num_a1) * fracadj
                        num_a2 = num_a1 + delnum_a2
                        numbnd = max(drv_c * v2nxxrl, min(drv_c * v2nyyrl, num_c1))
                        delnum_c2 = (numbnd - num_c1) * fracadj
                        num_c2 = num_c1 + delnum_c2
                        if (delnum_a2 == 0.0) and (delnum_c2 != 0.0):
                            num_a2 = max(
                                drv_a * v2nxxrl,
                                min(drv_a * v2nyyrl, num_a1 - delnum_c2),
                            )
                        elif (delnum_a2 != 0.0) and (delnum_c2 == 0.0):
                            num_c2 = max(
                                drv_c * v2nxxrl,
                                min(drv_c * v2nyyrl, num_c1 - delnum_a2),
                            )
                        drv_t = drv_a + drv_c
                        num_t2 = num_a2 + num_c2
                        delnum_a3 = 0.0
                        delnum_c3 = 0.0
                        if num_t2 < drv_t * v2nxx:
                            delnum_t3 = (drv_t * v2nxx - num_t2) * fracadj
                            if (num_a2 < drv_a * v2nxx) and (num_c2 < drv_c * v2nxx):
                                delnum_a3 = delnum_t3 * (num_a2 / num_t2)
                                delnum_c3 = delnum_t3 * (num_c2 / num_t2)
                            elif num_c2 < drv_c * v2nxx:
                                delnum_c3 = delnum_t3
                            elif num_a2 < drv_a * v2nxx:
                                delnum_a3 = delnum_t3
                        elif num_t2 > drv_t * v2nyy:
                            delnum_t3 = (drv_t * v2nyy - num_t2) * fracadj
                            if (num_a2 > drv_a * v2nyy) and (num_c2 > drv_c * v2nyy):
                                delnum_a3 = delnum_t3 * (num_a2 / num_t2)
                                delnum_c3 = delnum_t3 * (num_c2 / num_t2)
                            elif num_c2 > drv_c * v2nyy:
                                delnum_c3 = delnum_t3
                            elif num_a2 > drv_a * v2nyy:
                                delnum_a3 = delnum_t3
                        num_a = num_a2 + delnum_a3
                        dqdt[_idx3(i, k, lna, pcols, pver)] = (num_a - num_a0) * deltatinv
                        num_c = num_c2 + delnum_c3
                        dqqcwdt[_idx3(i, k, lnc, pcols, pver)] = (num_c - num_c0) * deltatinv

                if drv_a > 0.0:
                    if num_a <= drv_a * v2nxx:
                        dgncur_a[_idx3(i, k, n, pcols, pver)] = dgnxx
                    elif num_a >= drv_a * v2nyy:
                        dgncur_a[_idx3(i, k, n, pcols, pver)] = dgnyy
                    else:
                        dgncur_a[_idx3(i, k, n, pcols, pver)] = (drv_a / (dumfac * num_a)) ** third

                pdel_fac = pdel[_idx2(i, k, pcols)] / gravit
                qsrflx[_idx4(i, lna, 1, 1, pcols, pcnst, 4)] = (
                    qsrflx[_idx4(i, lna, 1, 1, pcols, pcnst, 4)]
                    + max(0.0, dqdt[_idx3(i, k, lna, pcols, pver)]) * pdel_fac
                )
                qsrflx[_idx4(i, lna, 2, 1, pcols, pcnst, 4)] = (
                    qsrflx[_idx4(i, lna, 2, 1, pcols, pcnst, 4)]
                    + min(0.0, dqdt[_idx3(i, k, lna, pcols, pver)]) * pdel_fac
                )
                qsrflx[_idx4(i, lnc, 1, 2, pcols, pcnst, 4)] = (
                    qsrflx[_idx4(i, lnc, 1, 2, pcols, pcnst, 4)]
                    + max(0.0, dqqcwdt[_idx3(i, k, lnc, pcols, pver)]) * pdel_fac
                )
                qsrflx[_idx4(i, lnc, 2, 2, pcols, pcnst, 4)] = (
                    qsrflx[_idx4(i, lnc, 2, 2, pcols, pcnst, 4)]
                    + min(0.0, dqqcwdt[_idx3(i, k, lnc, pcols, pver)]) * pdel_fac
                )

                if do_aitacc_transfer != 0:
                    if n == nait:
                        drv_a_aitsv[_idx2(i, k, pcols)] = drv_a
                        num_a_aitsv[_idx2(i, k, pcols)] = num_a
                        drv_c_aitsv[_idx2(i, k, pcols)] = drv_c
                        num_c_aitsv[_idx2(i, k, pcols)] = num_c
                    elif n == nacc:
                        drv_a_accsv[_idx2(i, k, pcols)] = drv_a
                        num_a_accsv[_idx2(i, k, pcols)] = num_a
                        drv_c_accsv[_idx2(i, k, pcols)] = drv_c
                        num_c_accsv[_idx2(i, k, pcols)] = num_c

                i += 1
            k += 1

        n += 1

    if do_aitacc_transfer != 0:
        iq = 1
        while iq <= nspecfrm_pair1:
            lsfrm = lspecfrma_pair1[iq - 1]
            lstoo = lspectooa_pair1[iq - 1]
            if (lsfrm > 0) and (lstoo > 0):
                dotend[lsfrm - 1] = 1
                dotend[lstoo - 1] = 1
            lsfrm = lspecfrmc_pair1[iq - 1]
            lstoo = lspectooc_pair1[iq - 1]
            if (lsfrm > 0) and (lstoo > 0):
                dotendqqcw[lsfrm - 1] = 1
                dotendqqcw[lstoo - 1] = 1
            iq += 1

        v2nzz = sqrt(voltonumb_amode[nait - 1] * voltonumb_amode[nacc - 1])

        k = top_lev
        while k <= pver:
            i = 1
            while i <= ncol:
                pdel_fac = pdel[_idx2(i, k, pcols)] / gravit
                xfertend_num_11 = 0.0
                xfertend_num_12 = 0.0
                xfertend_num_21 = 0.0
                xfertend_num_22 = 0.0
                xferfrac_num_ait2acc = 0.0
                xferfrac_vol_ait2acc = 0.0
                xferfrac_num_acc2ait = 0.0
                xferfrac_vol_acc2ait = 0.0
                xfertend = 0.0
                xfercoef = 0.0
                ixfer_ait2acc = 0
                xfercoef_num_ait2acc = 0.0
                xfercoef_vol_ait2acc = 0.0

                drv_t = drv_a_aitsv[_idx2(i, k, pcols)] + drv_c_aitsv[_idx2(i, k, pcols)]
                num_t = num_a_aitsv[_idx2(i, k, pcols)] + num_c_aitsv[_idx2(i, k, pcols)]
                if drv_t > 0.0:
                    if num_t < drv_t * v2nzz:
                        ixfer_ait2acc = 1
                        if num_t < drv_t * voltonumb_amode[nacc - 1]:
                            xferfrac_num_ait2acc = 1.0
                            xferfrac_vol_ait2acc = 1.0
                        else:
                            xferfrac_vol_ait2acc = ((num_t / drv_t) - v2nzz) / (
                                voltonumb_amode[nacc - 1] - v2nzz
                            )
                            xferfrac_num_ait2acc = xferfrac_vol_ait2acc * (
                                drv_t * voltonumb_amode[nacc - 1] / num_t
                            )
                            if (xferfrac_num_ait2acc <= 0.0) or (xferfrac_vol_ait2acc <= 0.0):
                                xferfrac_num_ait2acc = 0.0
                                xferfrac_vol_ait2acc = 0.0
                            elif (xferfrac_num_ait2acc >= 1.0) or (xferfrac_vol_ait2acc >= 1.0):
                                xferfrac_num_ait2acc = 1.0
                                xferfrac_vol_ait2acc = 1.0
                        xfercoef_num_ait2acc = xferfrac_num_ait2acc * tadjinv
                        xfercoef_vol_ait2acc = xferfrac_vol_ait2acc * tadjinv
                        xfertend_num_11 = num_a_aitsv[_idx2(i, k, pcols)] * xfercoef_num_ait2acc
                        xfertend_num_12 = num_c_aitsv[_idx2(i, k, pcols)] * xfercoef_num_ait2acc

                ixfer_acc2ait = 0
                xfercoef_num_acc2ait = 0.0
                xfercoef_vol_acc2ait = 0.0
                num_t0 = 0.0
                drv_a_noxf = 0.0
                drv_c_noxf = 0.0

                drv_t = drv_a_accsv[_idx2(i, k, pcols)] + drv_c_accsv[_idx2(i, k, pcols)]
                num_t = num_a_accsv[_idx2(i, k, pcols)] + num_c_accsv[_idx2(i, k, pcols)]
                if drv_t > 0.0:
                    if num_t > drv_t * v2nzz:
                        l1 = 1
                        while l1 <= nspec_amode[nacc - 1]:
                            la = lmassptr_amode[_idx2(l1, nacc, maxd_aspectype)]
                            noxf = 1
                            iq = 1
                            while iq <= nspecfrm_pair1:
                                if lspectooa_pair1[iq - 1] == la:
                                    noxf = 0
                                iq += 1
                            if noxf != 0:
                                lsptype = lspectype_amode[_idx2(l1, nacc, maxd_aspectype)]
                                dummwdens = 1.0 / specdens_amode[lsptype - 1]
                                drv_a_noxf = (
                                    drv_a_noxf
                                    + max(0.0, q[_idx3(i, k, la, pcols, pver)]) * dummwdens
                                )
                                lc = lmassptrcw_amode[_idx2(l1, nacc, maxd_aspectype)]
                                drv_c_noxf = (
                                    drv_c_noxf
                                    + max(0.0, qqcw[_idx3(i, k, lc, pcols, pver)]) * dummwdens
                                )
                            l1 += 1
                        drv_t_noxf = drv_a_noxf + drv_c_noxf
                        num_t_noxf = drv_t_noxf * voltonumblo_amode[nacc - 1]
                        num_t0 = num_t
                        num_t = max(0.0, num_t - num_t_noxf)
                        drv_t = max(0.0, drv_t - drv_t_noxf)

                if drv_t > 0.0:
                    if num_t > drv_t * v2nzz:
                        ixfer_acc2ait = 1
                        if num_t > drv_t * voltonumb_amode[nait - 1]:
                            xferfrac_num_acc2ait = 1.0
                            xferfrac_vol_acc2ait = 1.0
                        else:
                            xferfrac_vol_acc2ait = ((num_t / drv_t) - v2nzz) / (
                                voltonumb_amode[nait - 1] - v2nzz
                            )
                            xferfrac_num_acc2ait = xferfrac_vol_acc2ait * (
                                drv_t * voltonumb_amode[nait - 1] / num_t
                            )
                            if (xferfrac_num_acc2ait <= 0.0) or (xferfrac_vol_acc2ait <= 0.0):
                                xferfrac_num_acc2ait = 0.0
                                xferfrac_vol_acc2ait = 0.0
                            elif (xferfrac_num_acc2ait >= 1.0) or (xferfrac_vol_acc2ait >= 1.0):
                                xferfrac_num_acc2ait = 1.0
                                xferfrac_vol_acc2ait = 1.0
                        duma = 1.0e-37
                        xferfrac_num_acc2ait = xferfrac_num_acc2ait * num_t / max(duma, num_t0)
                        xfercoef_num_acc2ait = xferfrac_num_acc2ait * tadjinv
                        xfercoef_vol_acc2ait = xferfrac_vol_acc2ait * tadjinv
                        xfertend_num_21 = num_a_accsv[_idx2(i, k, pcols)] * xfercoef_num_acc2ait
                        xfertend_num_22 = num_c_accsv[_idx2(i, k, pcols)] * xfercoef_num_acc2ait

                if ixfer_ait2acc + ixfer_acc2ait > 0:
                    duma = (xfertend_num_11 - xfertend_num_21) * deltat
                    num_a = max(0.0, num_a_aitsv[_idx2(i, k, pcols)] - duma)
                    num_a_acc = max(0.0, num_a_accsv[_idx2(i, k, pcols)] + duma)

                    duma = (
                        drv_a_aitsv[_idx2(i, k, pcols)] * xfercoef_vol_ait2acc
                        - (drv_a_accsv[_idx2(i, k, pcols)] - drv_a_noxf) * xfercoef_vol_acc2ait
                    ) * deltat
                    drv_a = max(0.0, drv_a_aitsv[_idx2(i, k, pcols)] - duma)
                    drv_a_acc = max(0.0, drv_a_accsv[_idx2(i, k, pcols)] + duma)

                    duma = (xfertend_num_12 - xfertend_num_22) * deltat
                    num_c = max(0.0, num_c_aitsv[_idx2(i, k, pcols)] - duma)
                    num_c_acc = max(0.0, num_c_accsv[_idx2(i, k, pcols)] + duma)

                    duma = (
                        drv_c_aitsv[_idx2(i, k, pcols)] * xfercoef_vol_ait2acc
                        - (drv_c_accsv[_idx2(i, k, pcols)] - drv_c_noxf) * xfercoef_vol_acc2ait
                    ) * deltat
                    drv_c = max(0.0, drv_c_aitsv[_idx2(i, k, pcols)] - duma)
                    drv_c_acc = max(0.0, drv_c_accsv[_idx2(i, k, pcols)] + duma)

                    if drv_a > 0.0:
                        if num_a <= drv_a * voltonumbhi_amode[nait - 1]:
                            dgncur_a[_idx3(i, k, nait, pcols, pver)] = dgnumhi_amode[nait - 1]
                        elif num_a >= drv_a * voltonumblo_amode[nait - 1]:
                            dgncur_a[_idx3(i, k, nait, pcols, pver)] = dgnumlo_amode[nait - 1]
                        else:
                            dgncur_a[_idx3(i, k, nait, pcols, pver)] = (drv_a / (dumfac * num_a)) ** third
                    else:
                        dgncur_a[_idx3(i, k, nait, pcols, pver)] = dgnum_amode[nait - 1]

                    if drv_a_acc > 0.0:
                        if num_a_acc <= drv_a_acc * voltonumbhi_amode[nacc - 1]:
                            dgncur_a[_idx3(i, k, nacc, pcols, pver)] = dgnumhi_amode[nacc - 1]
                        elif num_a_acc >= drv_a_acc * voltonumblo_amode[nacc - 1]:
                            dgncur_a[_idx3(i, k, nacc, pcols, pver)] = dgnumlo_amode[nacc - 1]
                        else:
                            dgncur_a[_idx3(i, k, nacc, pcols, pver)] = (
                                drv_a_acc / (dumfac * num_a_acc)
                            ) ** third
                    else:
                        dgncur_a[_idx3(i, k, nacc, pcols, pver)] = dgnum_amode[nacc - 1]

                    j = 1
                    while j <= 2:
                        if ((j == 1) and (ixfer_ait2acc > 0)) or ((j == 2) and (ixfer_acc2ait > 0)):
                            jsrflx = j + 2
                            if j == 1:
                                xfercoef = xfercoef_vol_ait2acc
                            else:
                                xfercoef = xfercoef_vol_acc2ait

                            iq = 1
                            while iq <= nspecfrm_pair1:
                                jac = 1
                                while jac <= 2:
                                    lsfrm = 0
                                    lstoo = 0
                                    if j == 1:
                                        if jac == 1:
                                            lsfrm = lspecfrma_pair1[iq - 1]
                                            lstoo = lspectooa_pair1[iq - 1]
                                        else:
                                            lsfrm = lspecfrmc_pair1[iq - 1]
                                            lstoo = lspectooc_pair1[iq - 1]
                                    else:
                                        if jac == 1:
                                            lsfrm = lspectooa_pair1[iq - 1]
                                            lstoo = lspecfrma_pair1[iq - 1]
                                        else:
                                            lsfrm = lspectooc_pair1[iq - 1]
                                            lstoo = lspecfrmc_pair1[iq - 1]

                                    if (lsfrm > 0) and (lstoo > 0):
                                        if jac == 1:
                                            if iq == 1:
                                                if j == 1:
                                                    xfertend = xfertend_num_11
                                                else:
                                                    xfertend = xfertend_num_21
                                            else:
                                                xfertend = max(0.0, q[_idx3(i, k, lsfrm, pcols, pver)]) * xfercoef
                                            dqdt[_idx3(i, k, lsfrm, pcols, pver)] = (
                                                dqdt[_idx3(i, k, lsfrm, pcols, pver)] - xfertend
                                            )
                                            dqdt[_idx3(i, k, lstoo, pcols, pver)] = (
                                                dqdt[_idx3(i, k, lstoo, pcols, pver)] + xfertend
                                            )
                                        else:
                                            if iq == 1:
                                                if j == 1:
                                                    xfertend = xfertend_num_12
                                                else:
                                                    xfertend = xfertend_num_22
                                            else:
                                                xfertend = max(0.0, qqcw[_idx3(i, k, lsfrm, pcols, pver)]) * xfercoef
                                            dqqcwdt[_idx3(i, k, lsfrm, pcols, pver)] = (
                                                dqqcwdt[_idx3(i, k, lsfrm, pcols, pver)] - xfertend
                                            )
                                            dqqcwdt[_idx3(i, k, lstoo, pcols, pver)] = (
                                                dqqcwdt[_idx3(i, k, lstoo, pcols, pver)] + xfertend
                                            )

                                        qsrflx[_idx4(i, lsfrm, jsrflx, jac, pcols, pcnst, 4)] = (
                                            qsrflx[_idx4(i, lsfrm, jsrflx, jac, pcols, pcnst, 4)]
                                            - xfertend * pdel_fac
                                        )
                                        qsrflx[_idx4(i, lstoo, jsrflx, jac, pcols, pcnst, 4)] = (
                                            qsrflx[_idx4(i, lstoo, jsrflx, jac, pcols, pcnst, 4)]
                                            + xfertend * pdel_fac
                                        )

                                    jac += 1
                                iq += 1

                        j += 1

                i += 1
            k += 1

    lc = 1
    while lc <= pcnst:
        if dotendqqcw[lc - 1] != 0:
            k = top_lev
            while k <= pver:
                i = 1
                while i <= ncol:
                    qqcw[_idx3(i, k, lc, pcols, pver)] = max(
                        0.0,
                        qqcw[_idx3(i, k, lc, pcols, pver)]
                        + dqqcwdt[_idx3(i, k, lc, pcols, pver)] * deltat,
                    )
                    i += 1
                k += 1
        lc += 1


@export
def modal_aero_wateruptake_dr_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    nmodes: int,
    maxd_aspectype: int,
    pi_const: float,
    pi43_const: float,
    rhoh2o_const: float,
    rh_p: cobj,
    dgncur_a_p: cobj,
    dgncur_awet_p: cobj,
    qaerwat_p: cobj,
    wetdens_p: cobj,
    nspec_mode_p: cobj,
    sigmag_p: cobj,
    rhcrystal_p: cobj,
    rhdeliques_p: cobj,
    raer_work_p: cobj,
    specdens_work_p: cobj,
    spechygro_work_p: cobj,
    maer_p: cobj,
    hygro_p: cobj,
    naer_p: cobj,
    dryvol_p: cobj,
    drymass_p: cobj,
    dryrad_p: cobj,
    wetrad_p: cobj,
    wetvol_p: cobj,
    wtrvol_p: cobj,
    specdens_1_p: cobj,
    dryvolmr_p: cobj,
):
    rh = Ptr[float](rh_p)
    dgncur_a = Ptr[float](dgncur_a_p)
    dgncur_awet = Ptr[float](dgncur_awet_p)
    qaerwat = Ptr[float](qaerwat_p)
    wetdens = Ptr[float](wetdens_p)
    nspec_mode = Ptr[int](nspec_mode_p)
    sigmag = Ptr[float](sigmag_p)
    rhcrystal = Ptr[float](rhcrystal_p)
    rhdeliques = Ptr[float](rhdeliques_p)
    raer_work = Ptr[float](raer_work_p)
    specdens_work = Ptr[float](specdens_work_p)
    spechygro_work = Ptr[float](spechygro_work_p)
    maer = Ptr[float](maer_p)
    hygro_work = Ptr[float](hygro_p)
    naer = Ptr[float](naer_p)
    dryvol = Ptr[float](dryvol_p)
    drymass = Ptr[float](drymass_p)
    dryrad = Ptr[float](dryrad_p)
    wetrad = Ptr[float](wetrad_p)
    wetvol = Ptr[float](wetvol_p)
    wtrvol = Ptr[float](wtrvol_p)
    specdens_1 = Ptr[float](specdens_1_p)
    dryvolmr = Ptr[float](dryvolmr_p)

    m = 1
    while m <= nmodes:
        specdens_1[m - 1] = 0.0
        spechygro_1 = 0.0

        k = top_lev
        while k <= pver:
            i = 1
            while i <= ncol:
                idx2 = _idx2(i, k, pcols)
                idx3m = _idx3(i, k, m, pcols, pver)
                dryvolmr[idx2] = 0.0
                maer[idx3m] = 0.0
                hygro_work[idx3m] = 0.0
                naer[idx3m] = 0.0
                dryvol[idx3m] = 0.0
                drymass[idx3m] = 0.0
                dryrad[idx3m] = 0.0
                wetrad[idx3m] = 0.0
                wetvol[idx3m] = 0.0
                wtrvol[idx3m] = 0.0
                i += 1
            k += 1

        l = 1
        while l <= nspec_mode[m - 1]:
            specdens = specdens_work[_idx2(l, m, maxd_aspectype)]
            spechygro = spechygro_work[_idx2(l, m, maxd_aspectype)]
            if l == 1:
                specdens_1[m - 1] = specdens
                spechygro_1 = spechygro

            k = top_lev
            while k <= pver:
                i = 1
                while i <= ncol:
                    idx2 = _idx2(i, k, pcols)
                    idx3m = _idx3(i, k, m, pcols, pver)
                    duma = raer_work[_idx4(i, k, l, m, pcols, pver, maxd_aspectype)]
                    maer[idx3m] = maer[idx3m] + duma
                    dumb = duma / specdens
                    dryvolmr[idx2] = dryvolmr[idx2] + dumb
                    hygro_work[idx3m] = hygro_work[idx3m] + dumb * spechygro
                    i += 1
                k += 1
            l += 1

        alnsg = log(sigmag[m - 1])
        k = top_lev
        while k <= pver:
            i = 1
            while i <= ncol:
                idx2 = _idx2(i, k, pcols)
                idx3m = _idx3(i, k, m, pcols, pver)

                if dryvolmr[idx2] > 1.0e-30:
                    hygro_work[idx3m] = hygro_work[idx3m] / dryvolmr[idx2]
                else:
                    hygro_work[idx3m] = spechygro_1

                v2ncur_a = _modal_aero_v2ncur(dgncur_a[idx3m], pi_const, alnsg)
                naer[idx3m] = dryvolmr[idx2] * v2ncur_a

                if maer[idx3m] > 1.0e-31:
                    drydens = maer[idx3m] / dryvolmr[idx2]
                else:
                    drydens = 1.0

                dryvol[idx3m] = 1.0 / v2ncur_a
                drymass[idx3m] = drydens * dryvol[idx3m]
                dryrad[idx3m] = _modal_aero_radius_from_vol(dryvol[idx3m], pi43_const)
                i += 1
            k += 1

        m += 1


@export
def modal_aero_wateruptake_kohler_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    nmodes: int,
    solver_stage: int,
    dryrad_p: cobj,
    hygro_p: cobj,
    rh_p: cobj,
    wetrad_p: cobj,
):
    dryrad = Ptr[float](dryrad_p)
    hygro = Ptr[float](hygro_p)
    rh = Ptr[float](rh_p)
    wetrad = Ptr[float](wetrad_p)

    m = 1
    while m <= nmodes:
        k = top_lev
        while k <= pver:
            i = 1
            while i <= ncol:
                idx2 = _idx2(i, k, pcols)
                idx3m = _idx3(i, k, m, pcols, pver)
                if solver_stage == 1:
                    wetrad[idx3m] = modal_aero_kohler_native_cb(
                        dryrad[idx3m], hygro[idx3m], rh[idx2]
                    )
                elif solver_stage == 3:
                    wetrad[idx3m] = _modal_aero_kohler_scalar_native_roots(
                        dryrad[idx3m], hygro[idx3m], rh[idx2]
                    )
                elif solver_stage == 4:
                    wetrad[idx3m] = _modal_aero_kohler_scalar_sat_native(
                        dryrad[idx3m], hygro[idx3m], rh[idx2]
                    )
                elif solver_stage == 5:
                    wetrad[idx3m] = _modal_aero_kohler_scalar_subsat_native(
                        dryrad[idx3m], hygro[idx3m], rh[idx2]
                    )
                elif solver_stage == 6:
                    wetrad[idx3m] = _modal_aero_kohler_scalar_quartic_native(
                        dryrad[idx3m], hygro[idx3m], rh[idx2]
                    )
                elif solver_stage == 7:
                    wetrad[idx3m] = _modal_aero_kohler_scalar_cubic_native(
                        dryrad[idx3m], hygro[idx3m], rh[idx2]
                    )
                elif solver_stage == 8:
                    wetrad[idx3m] = _modal_aero_kohler_scalar_quartic_sqrt_native(
                        dryrad[idx3m], hygro[idx3m], rh[idx2]
                    )
                elif solver_stage == 9:
                    wetrad[idx3m] = _modal_aero_kohler_scalar_quartic_pow_native(
                        dryrad[idx3m], hygro[idx3m], rh[idx2]
                    )
                elif solver_stage == 10:
                    wetrad[idx3m] = _modal_aero_kohler_scalar_quartic_sqrt_pow_native(
                        dryrad[idx3m], hygro[idx3m], rh[idx2]
                    )
                else:
                    wetrad[idx3m] = _modal_aero_kohler_scalar_all_codon(
                        dryrad[idx3m], hygro[idx3m], rh[idx2]
                    )
                i += 1
            k += 1
        m += 1


@export
def modal_aero_wateruptake_base_guard_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    nmodes: int,
    dryrad_p: cobj,
    wetrad_p: cobj,
):
    dryrad = Ptr[float](dryrad_p)
    wetrad = Ptr[float](wetrad_p)

    m = 1
    while m <= nmodes:
        k = top_lev
        while k <= pver:
            i = 1
            while i <= ncol:
                idx3m = _idx3(i, k, m, pcols, pver)
                if wetrad[idx3m] < dryrad[idx3m]:
                    wetrad[idx3m] = dryrad[idx3m]
                i += 1
            k += 1
        m += 1


@export
def modal_aero_wateruptake_base_wtrvol_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    nmodes: int,
    dryvol_p: cobj,
    wetvol_p: cobj,
    wtrvol_p: cobj,
):
    dryvol = Ptr[float](dryvol_p)
    wetvol = Ptr[float](wetvol_p)
    wtrvol = Ptr[float](wtrvol_p)

    m = 1
    while m <= nmodes:
        k = top_lev
        while k <= pver:
            i = 1
            while i <= ncol:
                idx3m = _idx3(i, k, m, pcols, pver)
                wtrvol[idx3m] = wetvol[idx3m] - dryvol[idx3m]
                if wtrvol[idx3m] < 0.0:
                    wtrvol[idx3m] = 0.0
                i += 1
            k += 1
        m += 1


@export
def modal_aero_wateruptake_base_pow_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    nmodes: int,
    wetrad_p: cobj,
    wetvol_p: cobj,
):
    wetrad = Ptr[float](wetrad_p)
    wetvol = Ptr[float](wetvol_p)

    m = 1
    while m <= nmodes:
        k = top_lev
        while k <= pver:
            i = 1
            while i <= ncol:
                idx3m = _idx3(i, k, m, pcols, pver)
                wetvol[idx3m] = modal_aero_vol_from_radius_native_cb(wetrad[idx3m])
                i += 1
            k += 1
        m += 1


@export
def modal_aero_wateruptake_base_clamp_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    nmodes: int,
    dryvol_p: cobj,
    wetvol_p: cobj,
):
    dryvol = Ptr[float](dryvol_p)
    wetvol = Ptr[float](wetvol_p)

    m = 1
    while m <= nmodes:
        k = top_lev
        while k <= pver:
            i = 1
            while i <= ncol:
                idx3m = _idx3(i, k, m, pcols, pver)
                if wetvol[idx3m] < dryvol[idx3m]:
                    wetvol[idx3m] = dryvol[idx3m]
                i += 1
            k += 1
        m += 1


@export
def modal_aero_wateruptake_base_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    nmodes: int,
    dryrad_p: cobj,
    dryvol_p: cobj,
    wetrad_p: cobj,
    wetvol_p: cobj,
    wtrvol_p: cobj,
):
    dryrad = Ptr[float](dryrad_p)
    dryvol = Ptr[float](dryvol_p)
    wetrad = Ptr[float](wetrad_p)
    wetvol = Ptr[float](wetvol_p)
    wtrvol = Ptr[float](wtrvol_p)

    m = 1
    while m <= nmodes:
        k = top_lev
        while k <= pver:
            i = 1
            while i <= ncol:
                idx3m = _idx3(i, k, m, pcols, pver)

                if wetrad[idx3m] < dryrad[idx3m]:
                    wetrad[idx3m] = dryrad[idx3m]

                wetvol[idx3m] = modal_aero_vol_from_radius_native_cb(wetrad[idx3m])
                if wetvol[idx3m] < dryvol[idx3m]:
                    wetvol[idx3m] = dryvol[idx3m]

                wtrvol[idx3m] = wetvol[idx3m] - dryvol[idx3m]
                if wtrvol[idx3m] < 0.0:
                    wtrvol[idx3m] = 0.0

                i += 1
            k += 1
        m += 1


@export
def modal_aero_wateruptake_hyst_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    nmodes: int,
    pi43_const: float,
    rhcrystal_p: cobj,
    rhdeliques_p: cobj,
    dryrad_p: cobj,
    rh_p: cobj,
    dryvol_p: cobj,
    wetrad_p: cobj,
    wetvol_p: cobj,
    wtrvol_p: cobj,
):
    rhcrystal = Ptr[float](rhcrystal_p)
    rhdeliques = Ptr[float](rhdeliques_p)
    dryrad = Ptr[float](dryrad_p)
    rh = Ptr[float](rh_p)
    dryvol = Ptr[float](dryvol_p)
    wetrad = Ptr[float](wetrad_p)
    wetvol = Ptr[float](wetvol_p)
    wtrvol = Ptr[float](wtrvol_p)

    m = 1
    while m <= nmodes:
        hystfac = 1.0 / max(1.0e-5, rhdeliques[m - 1] - rhcrystal[m - 1])

        k = top_lev
        while k <= pver:
            i = 1
            while i <= ncol:
                idx2 = _idx2(i, k, pcols)
                idx3m = _idx3(i, k, m, pcols, pver)

                if rh[idx2] < rhcrystal[m - 1]:
                    wetrad[idx3m] = dryrad[idx3m]
                    wetvol[idx3m] = dryvol[idx3m]
                    wtrvol[idx3m] = 0.0
                elif rh[idx2] < rhdeliques[m - 1]:
                    wtrvol[idx3m] = wtrvol[idx3m] * hystfac * (rh[idx2] - rhcrystal[m - 1])
                    if wtrvol[idx3m] < 0.0:
                        wtrvol[idx3m] = 0.0
                    wetvol[idx3m] = dryvol[idx3m] + wtrvol[idx3m]
                    wetrad[idx3m] = _modal_aero_radius_from_vol(wetvol[idx3m], pi43_const)

                i += 1
            k += 1
        m += 1


@export
def modal_aero_wateruptake_postpow_wet_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    nmodes: int,
    do_hyst: int,
    pi43_const: float,
    rhcrystal_p: cobj,
    rhdeliques_p: cobj,
    dryrad_p: cobj,
    rh_p: cobj,
    dryvol_p: cobj,
    wetrad_p: cobj,
    wetvol_p: cobj,
    wtrvol_p: cobj,
):
    rhcrystal = Ptr[float](rhcrystal_p)
    rhdeliques = Ptr[float](rhdeliques_p)
    dryrad = Ptr[float](dryrad_p)
    rh = Ptr[float](rh_p)
    dryvol = Ptr[float](dryvol_p)
    wetrad = Ptr[float](wetrad_p)
    wetvol = Ptr[float](wetvol_p)
    wtrvol = Ptr[float](wtrvol_p)

    m = 1
    while m <= nmodes:
        k = top_lev
        while k <= pver:
            i = 1
            while i <= ncol:
                idx3m = _idx3(i, k, m, pcols, pver)

                if wetvol[idx3m] < dryvol[idx3m]:
                    wetvol[idx3m] = dryvol[idx3m]

                wtrvol[idx3m] = wetvol[idx3m] - dryvol[idx3m]
                if wtrvol[idx3m] < 0.0:
                    wtrvol[idx3m] = 0.0

                i += 1
            k += 1
        m += 1

    if do_hyst != 0:
        m = 1
        while m <= nmodes:
            hystfac = 1.0 / max(1.0e-5, rhdeliques[m - 1] - rhcrystal[m - 1])

            k = top_lev
            while k <= pver:
                i = 1
                while i <= ncol:
                    idx2 = _idx2(i, k, pcols)
                    idx3m = _idx3(i, k, m, pcols, pver)

                    if rh[idx2] < rhcrystal[m - 1]:
                        wetrad[idx3m] = dryrad[idx3m]
                        wetvol[idx3m] = dryvol[idx3m]
                        wtrvol[idx3m] = 0.0
                    elif rh[idx2] < rhdeliques[m - 1]:
                        wtrvol[idx3m] = wtrvol[idx3m] * hystfac * (rh[idx2] - rhcrystal[m - 1])
                        if wtrvol[idx3m] < 0.0:
                            wtrvol[idx3m] = 0.0
                        wetvol[idx3m] = dryvol[idx3m] + wtrvol[idx3m]
                        wetrad[idx3m] = _modal_aero_radius_from_vol(wetvol[idx3m], pi43_const)

                    i += 1
                k += 1
            m += 1


@export
def modal_aero_wateruptake_finalize_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    nmodes: int,
    rhoh2o_const: float,
    dgncur_a_p: cobj,
    dgncur_awet_p: cobj,
    qaerwat_p: cobj,
    wetdens_p: cobj,
    naer_p: cobj,
    dryrad_p: cobj,
    drymass_p: cobj,
    wetrad_p: cobj,
    wetvol_p: cobj,
    wtrvol_p: cobj,
    specdens_1_p: cobj,
):
    dgncur_a = Ptr[float](dgncur_a_p)
    dgncur_awet = Ptr[float](dgncur_awet_p)
    qaerwat = Ptr[float](qaerwat_p)
    wetdens = Ptr[float](wetdens_p)
    naer = Ptr[float](naer_p)
    dryrad = Ptr[float](dryrad_p)
    drymass = Ptr[float](drymass_p)
    wetrad = Ptr[float](wetrad_p)
    wetvol = Ptr[float](wetvol_p)
    wtrvol = Ptr[float](wtrvol_p)
    specdens_1 = Ptr[float](specdens_1_p)

    m = 1
    while m <= nmodes:
        k = top_lev
        while k <= pver:
            i = 1
            while i <= ncol:
                idx3m = _idx3(i, k, m, pcols, pver)

                dgncur_awet[idx3m] = dgncur_a[idx3m] * (wetrad[idx3m] / dryrad[idx3m])
                qaerwat[idx3m] = rhoh2o_const * naer[idx3m] * wtrvol[idx3m]

                if wetvol[idx3m] > 1.0e-30:
                    wetdens[idx3m] = (
                        drymass[idx3m] + rhoh2o_const * wtrvol[idx3m]
                    ) / wetvol[idx3m]
                else:
                    wetdens[idx3m] = specdens_1[m - 1]

                i += 1
            k += 1
        m += 1
