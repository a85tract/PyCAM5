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
def chem_surfvals_get_codon(value: float) -> float:
    return value


@export
def chem_surfvals_co2_rad_codon(value: float) -> float:
    return value


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
