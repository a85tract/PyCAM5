@inline
def _nint(x: float) -> int:
    if x >= 0.0:
        return int(x + 0.5)
    else:
        return int(x - 0.5)


@inline
def _idx2(i: int, k: int, pcols: int) -> int:
    return i + k * pcols


@inline
def _idx_vdist(kk: int, itype: int) -> int:
    return (kk - 1) + (itype - 1) * 16


@export
def lightning_no_prod_phase1_codon(
    ncol: int,
    pcols: int,
    pver: int,
    rga: float,
    rearth: float,
    geo_factor: float,
    factor: float,
    phis_p: cobj,
    zm_p: cobj,
    zi_p: cobj,
    t_p: cobj,
    cldtop_p: cobj,
    cldbot_p: cobj,
    landfrac_p: cobj,
    ocnfrac_p: cobj,
    wght_p: cobj,
    flash_freq_p: cobj,
    glob_prod_no_col_p: cobj,
    prod_no_col_p: cobj,
    cldhgt_p: cobj,
    dchgzone_p: cobj,
    cgic_p: cobj,
    flash_energy_p: cobj,
    status_p: cobj,
):
    phis = Ptr[float](phis_p)
    zm = Ptr[float](zm_p)
    zi = Ptr[float](zi_p)
    t = Ptr[float](t_p)
    cldtop = Ptr[float](cldtop_p)
    cldbot = Ptr[float](cldbot_p)
    landfrac = Ptr[float](landfrac_p)
    ocnfrac = Ptr[float](ocnfrac_p)
    wght = Ptr[float](wght_p)

    flash_freq = Ptr[float](flash_freq_p)
    glob_prod_no_col = Ptr[float](glob_prod_no_col_p)
    prod_no_col = Ptr[float](prod_no_col_p)
    cldhgt = Ptr[float](cldhgt_p)
    dchgzone = Ptr[float](dchgzone_p)
    cgic = Ptr[float](cgic_p)
    flash_energy = Ptr[float](flash_energy_p)
    status = Ptr[int](status_p)

    t0 = 273.0
    m2km = 1.0e-3
    secpyr = 365.0 * 8.64e4
    ca = 0.021
    cb = -0.648
    cc = 7.49
    cd = -36.54
    ce = 64.09

    status[0] = 0

    for i0 in range(ncol):
        cldtind = _nint(cldtop[i0])
        cldbind = _nint(cldbot[i0])

        while True:
            if cldbind <= cldtind or t[_idx2(i0, cldbind - 1, pcols)] < t0:
                break
            cldbind -= 1

        if cldtind < pver and cldtind > 0 and cldtind < cldbind:
            zsurf = phis[i0] * rga
            zint_top = zi[_idx2(i0, cldtind - 1, pcols)] + zsurf
            zmid_bot = zm[_idx2(i0, cldbind - 1, pcols)] + zsurf

            cldhgt_i = m2km * max(0.0, zint_top)
            dchgz = cldhgt_i - m2km * zmid_bot

            cldhgt[i0] = cldhgt_i
            dchgzone[i0] = dchgz

            flash_freq_land = 3.44e-5 * (cldhgt_i ** 4.9)
            flash_freq_ocn = 6.40e-4 * (cldhgt_i ** 1.7)
            flash_freq_i = landfrac[i0] * flash_freq_land + ocnfrac[i0] * flash_freq_ocn
            flash_freq[i0] = flash_freq_i

            cgic_i = 1.0 / ((((ca * dchgz + cb) * dchgz + cc) * dchgz + cd) * dchgz + ce)
            if dchgz < 5.5:
                cgic_i = 0.0
            elif dchgz > 14.0:
                cgic_i = 0.02
            cgic[i0] = cgic_i

            flash_energy_i = 6.7e9 * flash_freq_i / 60.0
            flash_energy_i = flash_energy_i * wght[i0] * geo_factor
            flash_energy[i0] = flash_energy_i

            prod_no_col_i = 1.0e17 * flash_energy_i / (1.0e4 * rearth * rearth * wght[i0]) * factor
            prod_no_col[i0] = prod_no_col_i

            glob_prod_no_col[i0] = (
                1.0e17
                * flash_energy_i
                * 14.00674
                * 1.65979e-24
                * 1.0e-12
                * secpyr
                * factor
            )


@export
def lightning_no_prod_phase2_codon(
    ncol: int,
    pcols: int,
    pver: int,
    rga: float,
    lat25: float,
    phis_p: cobj,
    zi_p: cobj,
    cldtop_p: cobj,
    landfrac_p: cobj,
    rlats_p: cobj,
    vdist_p: cobj,
    prod_no_col_p: cobj,
    cldhgt_p: cobj,
    prod_no_p: cobj,
    status_p: cobj,
):
    phis = Ptr[float](phis_p)
    zi = Ptr[float](zi_p)
    cldtop = Ptr[float](cldtop_p)
    landfrac = Ptr[float](landfrac_p)
    rlats = Ptr[float](rlats_p)
    vdist = Ptr[float](vdist_p)
    prod_no_col = Ptr[float](prod_no_col_p)
    cldhgt = Ptr[float](cldhgt_p)
    prod_no = Ptr[float](prod_no_p)
    status = Ptr[int](status_p)

    m2km = 1.0e-3
    km2cm = 1.0e5

    status[0] = 0

    for i0 in range(ncol):
        if prod_no_col[i0] > 0.0 and cldhgt[i0] > 0.0:
            cldtind = _nint(cldtop[i0])
            if cldtind <= 0 or cldtind > pver:
                status[0] = 1
                return

            if abs(rlats[i0]) > lat25:
                itype = 1
            elif _nint(landfrac[i0]) == 1:
                itype = 3
            else:
                itype = 2

            zsurf = phis[i0] * rga
            for k in range(cldtind, pver + 1):
                zlow = (zi[_idx2(i0, k, pcols)] + zsurf) * m2km
                zlow_scal = zlow * 16.0 / cldhgt[i0]
                zlow_ind = max(1, int(zlow_scal) + 1)

                zhigh = (zi[_idx2(i0, k - 1, pcols)] + zsurf) * m2km
                zhigh_scal = zhigh * 16.0 / cldhgt[i0]
                zhigh_ind = max(1, min(16, int(zhigh_scal) + 1))

                accum = 0.0
                for kk in range(zlow_ind, zhigh_ind + 1):
                    wrk = float(kk)
                    wrk1 = float(kk - 1)
                    fraction = min(zhigh_scal, wrk) - max(zlow_scal, wrk1)
                    fraction = max(0.0, min(1.0, fraction))
                    accum += fraction * vdist[_idx_vdist(kk, itype)] * 0.01

                prod_no[_idx2(i0, k - 1, pcols)] = (
                    prod_no_col[i0] * accum / (km2cm * (zhigh - zlow))
                )
