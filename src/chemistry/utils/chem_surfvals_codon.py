@export
def chem_surfvals_readnl_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def chem_surfvals_init_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def prescribed_strataero_adv_codon(active: int) -> int:
    return active


@export
def chem_surfvals_get_codon(
    name_len: int,
    name_ascii_p: cobj,
    mwdry: float,
    mwco2: float,
    co2vmr: float,
    n2ovmr: float,
    ch4vmr: float,
    f11vmr: float,
    f12vmr: float,
    o2mmr: float,
    status_p: cobj,
) -> float:
    name_ascii = Ptr[int](name_ascii_p)
    status = Ptr[int](status_p)

    status[0] = 0
    rmwco2 = mwco2 / mwdry

    if name_len == 6:
        if (
            name_ascii[0] == 67
            and name_ascii[1] == 79
            and name_ascii[2] == 50
            and name_ascii[3] == 86
            and name_ascii[4] == 77
            and name_ascii[5] == 82
        ):
            return co2vmr
        if (
            name_ascii[0] == 67
            and name_ascii[1] == 79
            and name_ascii[2] == 50
            and name_ascii[3] == 77
            and name_ascii[4] == 77
            and name_ascii[5] == 82
        ):
            return rmwco2 * co2vmr
        if (
            name_ascii[0] == 78
            and name_ascii[1] == 50
            and name_ascii[2] == 79
            and name_ascii[3] == 86
            and name_ascii[4] == 77
            and name_ascii[5] == 82
        ):
            return n2ovmr
        if (
            name_ascii[0] == 67
            and name_ascii[1] == 72
            and name_ascii[2] == 52
            and name_ascii[3] == 86
            and name_ascii[4] == 77
            and name_ascii[5] == 82
        ):
            return ch4vmr
        if (
            name_ascii[0] == 70
            and name_ascii[1] == 49
            and name_ascii[2] == 49
            and name_ascii[3] == 86
            and name_ascii[4] == 77
            and name_ascii[5] == 82
        ):
            return f11vmr
        if (
            name_ascii[0] == 70
            and name_ascii[1] == 49
            and name_ascii[2] == 50
            and name_ascii[3] == 86
            and name_ascii[4] == 77
            and name_ascii[5] == 82
        ):
            return f12vmr
        if (
            name_ascii[0] == 79
            and name_ascii[1] == 50
            and name_ascii[2] == 77
            and name_ascii[3] == 77
            and name_ascii[4] == 82
            and name_ascii[5] == 32
        ):
            return o2mmr
    elif name_len == 5:
        if (
            name_ascii[0] == 79
            and name_ascii[1] == 50
            and name_ascii[2] == 77
            and name_ascii[3] == 77
            and name_ascii[4] == 82
        ):
            return o2mmr

    status[0] = 1
    return 0.0


@export
def chem_surfvals_co2_rad_codon(
    vmr_present: int,
    vmr_value: int,
    mwdry: float,
    mwco2: float,
    co2vmr_rad: float,
    co2vmr: float,
) -> float:
    convert_vmr = mwco2 / mwdry
    if vmr_present != 0:
        if vmr_value != 0:
            convert_vmr = 1.0

    if co2vmr_rad > 0.0:
        return convert_vmr * co2vmr_rad

    return convert_vmr * co2vmr


@export
def chem_surfvals_set_all_codon(
    fixYear_ghg: int,
    ghg_yearStart_model: int,
    ghg_yearStart_data: int,
    yr: int,
    calday: float,
    ntim: int,
    yrdata_p: cobj,
    co2_p: cobj,
    ch4_p: cobj,
    n2o_p: cobj,
    f11_p: cobj,
    f12_p: cobj,
    adj_p: cobj,
    co2vmr_p: cobj,
    ch4vmr_p: cobj,
    n2ovmr_p: cobj,
    f11vmr_p: cobj,
    f12vmr_p: cobj,
    status_p: cobj,
):
    yrdata = Ptr[int](yrdata_p)
    co2 = Ptr[float](co2_p)
    ch4 = Ptr[float](ch4_p)
    n2o = Ptr[float](n2o_p)
    f11 = Ptr[float](f11_p)
    f12 = Ptr[float](f12_p)
    adj = Ptr[float](adj_p)

    co2vmr = Ptr[float](co2vmr_p)
    ch4vmr = Ptr[float](ch4vmr_p)
    n2ovmr = Ptr[float](n2ovmr_p)
    f11vmr = Ptr[float](f11vmr_p)
    f12vmr = Ptr[float](f12vmr_p)
    status = Ptr[int](status_p)

    yrmodel = 0
    nyrm = 0

    status[0] = 0

    if ghg_yearStart_model > 0 and ghg_yearStart_data > 0:
        if fixYear_ghg > 0:
            yrmodel = fixYear_ghg
            nyrm = fixYear_ghg - yrdata[0] + 1
        else:
            yearRan_model = yr - ghg_yearStart_model
            if yearRan_model < 0:
                status[0] = 1
                return
            yrmodel = yearRan_model + ghg_yearStart_data
            nyrm = ghg_yearStart_data + yearRan_model - yrdata[0] + 1
    else:
        if fixYear_ghg > 0:
            yrmodel = fixYear_ghg
            nyrm = fixYear_ghg - yrdata[0] + 1
        else:
            yrmodel = yr
            nyrm = yr - yrdata[0] + 1

    nyrp = nyrm + 1

    if nyrm < 1:
        status[0] = 2
        return

    if nyrp > ntim:
        status[0] = 3
        return

    doymodel = float(yrmodel) * 365.0 + calday
    doydatam = float(yrdata[nyrm - 1]) * 365.0 + 1.0
    doydatap = float(yrdata[nyrp - 1]) * 365.0 + 1.0

    if doymodel < 1.0:
        status[0] = 4
        return

    deltat = doydatap - doydatam
    fact1 = (doydatap - doymodel) / deltat
    fact2 = (doymodel - doydatam) / deltat

    if (
        abs(fact1 + fact2 - 1.0) > 1.0e-6
        or fact1 > 1.000001
        or fact1 < -1.0e-6
        or fact2 > 1.000001
        or fact2 < -1.0e-6
        or fact1 != fact1
        or fact2 != fact2
    ):
        status[0] = 5
        return

    co2vmr[0] = (co2[nyrm - 1] * fact1 + co2[nyrp - 1] * fact2) * 1.0e-6
    ch4vmr[0] = (ch4[nyrm - 1] * fact1 + ch4[nyrp - 1] * fact2) * 1.0e-9
    n2ovmr[0] = (n2o[nyrm - 1] * fact1 + n2o[nyrp - 1] * fact2) * 1.0e-9

    cfcscl = adj[nyrm - 1] * fact1 + adj[nyrp - 1] * fact2
    f11vmr[0] = (f11[nyrm - 1] * fact1 + f11[nyrp - 1] * fact2) * 1.0e-12 * (1.0 + cfcscl)
    f12vmr[0] = (f12[nyrm - 1] * fact1 + f12[nyrp - 1] * fact2) * 1.0e-12


@export
def chem_surfvals_set_co2_codon(
    daydiff: float,
    co2_base: float,
    co2_daily_factor: float,
    co2_limit: float,
    co2vmr_p: cobj,
):
    co2vmr = Ptr[float](co2vmr_p)

    if daydiff > 0.0:
        co2vmr[0] = co2_base * (co2_daily_factor ** daydiff)

        if co2_daily_factor < 1.0:
            co2vmr[0] = max(co2vmr[0], co2_limit)
        else:
            co2vmr[0] = min(co2vmr[0], co2_limit)
