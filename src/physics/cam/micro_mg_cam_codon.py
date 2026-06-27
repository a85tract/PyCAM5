from C import wv_sat_svp_water_codon(float, int) -> float
from C import wv_sat_svp_ice_codon(float, int) -> float
from C import wv_sat_svp_to_qsat_codon(float, float, float, float) -> float
from C import pow(float, float) -> float
from math import gamma, exp, log10, pi, sqrt


@inline
def pow_r8(x: float, y: float) -> float:
    return pow(x, y)


@inline
def pow_i(x: float, n: int) -> float:
    if n == 0:
        return 1.0
    if n < 0:
        return 1.0 / pow_i(x, -n)
    out = 1.0
    base = x
    expn = n
    while expn > 0:
        if expn % 2 == 1:
            out *= base
        expn //= 2
        if expn > 0:
            base *= base
    return out


@inline
def pow_i_f90_pracs(x: float, n: int) -> float:
    if n == 2:
        return x * x
    if n == 3:
        x2 = x * x
        return x2 * x
    if n == 4:
        x2 = x * x
        return x2 * x2
    if n == 5:
        x2 = x * x
        x3 = x2 * x
        return x2 * x3
    if n == 6:
        x2 = x * x
        x3 = x2 * x
        return x3 * x3
    return pow_i(x, n)


@inline
def _ascii_ptr_to_str(n: int, ptr_p: cobj) -> str:
    ptr = Ptr[int](ptr_p)
    out = ""
    for i in range(n):
        out += chr(ptr[i])
    return out.strip()


@inline
def _strip_comment(line: str) -> str:
    pos = line.find("!")
    if pos >= 0:
        return line[:pos]
    return line


@inline
def _find_group_end(line: str) -> int:
    quote = ""
    for i in range(len(line)):
        ch = line[i]
        if quote:
            if ch == quote:
                quote = ""
        elif ch == "'" or ch == '"':
            quote = ch
        elif ch == "/":
            return i
    return -1


@inline
def _parse_fortran_float(value: str) -> float:
    return float(value.strip().replace("D", "E").replace("d", "e"))


@inline
def _parse_fortran_int(value: str) -> int:
    return int(value.strip())


@inline
def _parse_fortran_bool(value: str) -> int:
    text = value.strip().lower()
    if text == ".true." or text == "t" or text == "true":
        return 1
    return 0


def _idx2(i: int, k: int, ld1: int):
    return (k - 1) * ld1 + (i - 1)


def _idx3(i: int, k: int, m: int, ld1: int, ld2: int):
    return ((m - 1) * ld2 + (k - 1)) * ld1 + (i - 1)


@export
def micro_mg_cam_readnl_codon(
    path_len: int,
    path_ascii_p: cobj,
    reals_p: cobj,
    ints_p: cobj,
    logicals_p: cobj,
    precip_frac_len: int,
    precip_frac_ascii_p: cobj,
) -> int:
    path = _ascii_ptr_to_str(path_len, path_ascii_p)
    reals = Ptr[float](reals_p)
    ints = Ptr[i32](ints_p)
    logicals = Ptr[i32](logicals_p)
    precip_frac_ascii = Ptr[int](precip_frac_ascii_p)

    f = open(path, "r")
    text = f.read()
    f.close()

    in_group = False
    found_group = False
    assignments = ""

    for raw_line in text.split("\n"):
        line = _strip_comment(raw_line).strip()
        lowered = line.lower()
        if not in_group:
            if lowered.startswith("&micro_mg_nl"):
                in_group = True
                found_group = True
                rest = line[len("&micro_mg_nl") :]
                if rest:
                    assignments += rest + ","
            continue

        slash = _find_group_end(line)
        if slash >= 0:
            assignments += line[:slash] + ","
            break
        assignments += line + ","

    if not found_group:
        return 0

    for item in assignments.split(","):
        if "=" not in item:
            continue
        parts = item.split("=", 1)
        key = parts[0].strip().lower()
        value = parts[1].strip()
        if key == "micro_mg_dcs":
            reals[0] = _parse_fortran_float(value)
        elif key == "micro_mg_berg_eff_factor":
            reals[1] = _parse_fortran_float(value)
        elif key == "micro_mg_version":
            ints[0] = i32(_parse_fortran_int(value))
        elif key == "micro_mg_sub_version":
            ints[1] = i32(_parse_fortran_int(value))
        elif key == "micro_mg_num_steps":
            ints[2] = i32(_parse_fortran_int(value))
        elif key == "micro_mg_do_cldice":
            logicals[0] = i32(_parse_fortran_bool(value))
        elif key == "micro_mg_do_cldliq":
            logicals[1] = i32(_parse_fortran_bool(value))
        elif key == "microp_uniform":
            logicals[2] = i32(_parse_fortran_bool(value))
        elif key == "micro_mg_precip_frac_method":
            text_value = value.strip().strip("'").strip('"')
            for i in range(precip_frac_len):
                if i < len(text_value):
                    precip_frac_ascii[i] = ord(text_value[i])
                else:
                    precip_frac_ascii[i] = 32

    return 0


@export
def micro_mg_cam_implements_cnst_codon(name_len: int, name_ascii_p: cobj) -> int:
    name_ascii = Ptr[int](name_ascii_p)
    if _name_eq8(name_len, name_ascii, 67, 76, 68, 76, 73, 81, 32, 32) != 0:
        return 1
    if _name_eq8(name_len, name_ascii, 67, 76, 68, 73, 67, 69, 32, 32) != 0:
        return 1
    if _name_eq8(name_len, name_ascii, 78, 85, 77, 76, 73, 81, 32, 32) != 0:
        return 1
    if _name_eq8(name_len, name_ascii, 78, 85, 77, 73, 67, 69, 32, 32) != 0:
        return 1
    if _name_eq8(name_len, name_ascii, 82, 65, 73, 78, 81, 77, 32, 32) != 0:
        return 1
    if _name_eq8(name_len, name_ascii, 83, 78, 79, 87, 81, 77, 32, 32) != 0:
        return 1
    if _name_eq8(name_len, name_ascii, 78, 85, 77, 82, 65, 73, 32, 32) != 0:
        return 1
    if _name_eq8(name_len, name_ascii, 78, 85, 77, 83, 78, 79, 32, 32) != 0:
        return 1
    return 0


@export
def micro_mg_cam_register_codon(
    micro_mg_version: int,
    micro_mg_sub_version: int,
    prog_modal_aero: int,
    use_subcol_microp: int,
    do_cldice: int,
    subcol_silhs: int,
    cnst_count_p: cobj,
    cnst_codes_p: cobj,
    pbuf_count_p: cobj,
    pbuf_codes_p: cobj,
    subcol_count_p: cobj,
    subcol_codes_p: cobj,
):
    cnst_count = Ptr[int](cnst_count_p)
    cnst_codes = Ptr[int](cnst_codes_p)
    pbuf_count = Ptr[int](pbuf_count_p)
    pbuf_codes = Ptr[int](pbuf_codes_p)
    subcol_count = Ptr[int](subcol_count_p)
    subcol_codes = Ptr[int](subcol_codes_p)

    n = 0
    cnst_codes[n] = 1
    n += 1
    cnst_codes[n] = 2
    n += 1
    if micro_mg_version == 1 and micro_mg_sub_version == 0:
        cnst_codes[n] = 3
        n += 1
        cnst_codes[n] = 4
        n += 1
    else:
        cnst_codes[n] = 103
        n += 1
        cnst_codes[n] = 104
        n += 1
    if micro_mg_version > 1:
        for code in range(105, 109):
            cnst_codes[n] = code
            n += 1
    cnst_count[0] = n

    n = 0
    for code in range(1, 17):
        pbuf_codes[n] = code
        n += 1
    if prog_modal_aero != 0:
        pbuf_codes[n] = 17
        n += 1
    for code in range(18, 36):
        pbuf_codes[n] = code
        n += 1
    if do_cldice == 0:
        for code in range(36, 39):
            pbuf_codes[n] = code
            n += 1
    for code in range(39, 44):
        pbuf_codes[n] = code
        n += 1
    if subcol_silhs != 0:
        for code in range(44, 48):
            pbuf_codes[n] = code
            n += 1
    pbuf_count[0] = n

    n = 0
    if use_subcol_microp != 0:
        for code in range(1, 24):
            subcol_codes[n] = code
            n += 1
        if prog_modal_aero != 0:
            subcol_codes[n] = 24
            n += 1
        for code in range(25, 33):
            subcol_codes[n] = code
            n += 1
    subcol_count[0] = n


@export
def micro_mg_cam_init_plan_codon(
    micro_mg_version: int,
    micro_mg_sub_version: int,
    init_plan_p: cobj,
):
    init_plan = Ptr[int](init_plan_p)
    if micro_mg_version == 1:
        init_plan[0] = 4
    elif micro_mg_version == 2:
        init_plan[0] = 8
    else:
        init_plan[0] = 4


@export
def micro_mg_cam_init_pbuf_set_plan_codon(
    use_subcol_microp: int,
    qrain_present: int,
    qsnow_present: int,
    nrain_present: int,
    nsnow_present: int,
    count_p: cobj,
    codes_p: cobj,
):
    count = Ptr[int](count_p)
    codes = Ptr[int](codes_p)

    n = 0
    for code in range(1, 18):
        codes[n] = code
        n += 1

    if qrain_present != 0:
        codes[n] = 18
        n += 1
    if qsnow_present != 0:
        codes[n] = 19
        n += 1
    if nrain_present != 0:
        codes[n] = 20
        n += 1
    if nsnow_present != 0:
        codes[n] = 21
        n += 1

    if use_subcol_microp != 0:
        for code in range(101, 109):
            codes[n] = code
            n += 1

    count[0] = n


@export
def micro_mg_cam_init_default_plan_codon(
    micro_mg_version: int,
    micro_mg_sub_version: int,
    history_amwg: int,
    history_budget: int,
    ncnst: int,
    budget_histfile: int,
    count_p: cobj,
    codes_p: cobj,
    files_p: cobj,
):
    count = Ptr[int](count_p)
    codes = Ptr[int](codes_p)
    files = Ptr[int](files_p)

    n = 0
    if history_amwg != 0:
        for code in range(1, 17):
            codes[n] = code
            files[n] = 1
            n += 1
        for m in range(1, ncnst + 1):
            codes[n] = 100 + m
            files[n] = 1
            n += 1

    if history_budget != 0:
        for code in range(201, 210):
            codes[n] = code
            files[n] = budget_histfile
            n += 1
        if micro_mg_version > 1:
            for code in range(210, 212):
                codes[n] = code
                files[n] = budget_histfile
                n += 1
        for code in range(212, 239):
            codes[n] = code
            files[n] = budget_histfile
            n += 1
        for code in range(301, 307):
            codes[n] = code
            files[n] = budget_histfile
            n += 1
        if micro_mg_version > 1:
            for code in range(307, 313):
                codes[n] = code
                files[n] = budget_histfile
                n += 1

    count[0] = n


@export
def p1_codon(n: int) -> int:
    if n >= 0:
        return 1
    return 0


@export
def p2_codon(n1: int, n2: int) -> int:
    if n1 >= 0 and n2 >= 0:
        return 1
    return 0


@inline
def _name_eq8(
    name_len: int,
    name_ascii: Ptr[int],
    c0: int,
    c1: int,
    c2: int,
    c3: int,
    c4: int,
    c5: int,
    c6: int,
    c7: int,
) -> int:
    if name_len > 8:
        i = 8
        while i < name_len:
            if name_ascii[i] != 32:
                return 0
            i += 1
    values = (c0, c1, c2, c3, c4, c5, c6, c7)
    i = 0
    while i < 8:
        left = 32
        if i < name_len:
            left = name_ascii[i]
        if left != values[i]:
            return 0
        i += 1
    return 1


@export
def new_mgpacker_codon(
    pcols: int,
    pver: int,
    mgcols_count: int,
    top_lev: int,
    plan_p: cobj,
):
    plan = Ptr[int](plan_p)
    plan[0] = pcols
    plan[1] = pver
    plan[2] = mgcols_count
    plan[3] = pver - top_lev + 1


@export
def mgpacker_finalize_codon(plan_p: cobj):
    plan = Ptr[int](plan_p)
    plan[0] = 1


@export
def mgfieldpostproc_1d_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def mgfieldpostproc_2d_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def mgfieldpostproc_finalize_codon(plan_p: cobj):
    plan = Ptr[int](plan_p)
    plan[0] = 1


@export
def mgfieldpostproc_process_and_unpack_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def mgfieldpostproc_unpack_only_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def new_mgpostproc_codon(plan_p: cobj):
    plan = Ptr[int](plan_p)
    plan[0] = 1


@export
def mgpostproc_finalize_codon(field_count: int, plan_p: cobj):
    plan = Ptr[int](plan_p)
    plan[0] = 1
    plan[1] = 1
    plan[2] = field_count


@export
def add_field_1d_codon(
    fill_present: int,
    accum_present: int,
    fill_in: float,
    accum_in: int,
    plan_p: cobj,
    fill_out_p: cobj,
):
    plan = Ptr[int](plan_p)
    fill_out = Ptr[float](fill_out_p)
    plan[0] = 1
    if accum_present != 0:
        plan[1] = accum_in
    else:
        plan[1] = 1
    if fill_present != 0:
        fill_out[0] = fill_in
    else:
        fill_out[0] = 0.0


@export
def add_field_2d_codon(
    fill_present: int,
    accum_present: int,
    fill_in: float,
    accum_in: int,
    plan_p: cobj,
    fill_out_p: cobj,
):
    plan = Ptr[int](plan_p)
    fill_out = Ptr[float](fill_out_p)
    plan[0] = 1
    if accum_present != 0:
        plan[1] = accum_in
    else:
        plan[1] = 1
    if fill_present != 0:
        fill_out[0] = fill_in
    else:
        fill_out[0] = 0.0


@export
def mgpostproc_accumulate_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def mgpostproc_process_and_unpack_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def mgpostproc_unpack_only_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def mgpostproc_copy_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@inline
def _micro_mg_data_copy_pack_1d(
    mgncol: int,
    mgcols: Ptr[int],
    src: Ptr[float],
    dst: Ptr[float],
):
    for j in range(1, mgncol + 1):
        dst[j - 1] = src[mgcols[j - 1] - 1]


@inline
def _micro_mg_data_copy_pack_2d(
    pcols: int,
    mgncol: int,
    top_lev: int,
    extent2: int,
    mgcols: Ptr[int],
    src: Ptr[float],
    dst: Ptr[float],
):
    for k in range(1, extent2 + 1):
        src_k = top_lev + k - 1
        for j in range(1, mgncol + 1):
            dst[_idx2(j, k, mgncol)] = src[_idx2(mgcols[j - 1], src_k, pcols)]


@inline
def _micro_mg_data_copy_pack_3d(
    pcols: int,
    pver: int,
    mgncol: int,
    nlev: int,
    top_lev: int,
    extent3: int,
    mgcols: Ptr[int],
    src: Ptr[float],
    dst: Ptr[float],
):
    for m in range(1, extent3 + 1):
        for k in range(1, nlev + 1):
            src_k = top_lev + k - 1
            for j in range(1, mgncol + 1):
                dst[_idx3(j, k, m, mgncol, nlev)] = src[
                    _idx3(mgcols[j - 1], src_k, m, pcols, pver)
                ]


@inline
def _micro_mg_data_fill_1d(pcols: int, fillvalue: float, dst: Ptr[float]):
    for i in range(1, pcols + 1):
        dst[i - 1] = fillvalue


@inline
def _micro_mg_data_fill_2d(
    pcols: int,
    extent2: int,
    fillvalue: float,
    dst: Ptr[float],
):
    for k in range(1, extent2 + 1):
        for i in range(1, pcols + 1):
            dst[_idx2(i, k, pcols)] = fillvalue


@inline
def _micro_mg_data_fill_3d(
    pcols: int,
    pver: int,
    extent3: int,
    fillvalue: float,
    dst: Ptr[float],
):
    for m in range(1, extent3 + 1):
        for k in range(1, pver + 1):
            for i in range(1, pcols + 1):
                dst[_idx3(i, k, m, pcols, pver)] = fillvalue


@inline
def _micro_mg_data_copy_1d(pcols: int, src: Ptr[float], dst: Ptr[float]):
    for i in range(1, pcols + 1):
        dst[i - 1] = src[i - 1]


@inline
def _micro_mg_data_copy_2d(
    pcols: int,
    extent2: int,
    src: Ptr[float],
    dst: Ptr[float],
):
    for k in range(1, extent2 + 1):
        for i in range(1, pcols + 1):
            dst[_idx2(i, k, pcols)] = src[_idx2(i, k, pcols)]


@inline
def _micro_mg_data_copy_3d(
    pcols: int,
    pver: int,
    extent3: int,
    src: Ptr[float],
    dst: Ptr[float],
):
    for m in range(1, extent3 + 1):
        for k in range(1, pver + 1):
            for i in range(1, pcols + 1):
                dst[_idx3(i, k, m, pcols, pver)] = src[
                    _idx3(i, k, m, pcols, pver)
                ]


@inline
def _micro_mg_data_scatter_unpack_1d(
    mgncol: int,
    mgcols: Ptr[int],
    src: Ptr[float],
    dst: Ptr[float],
):
    for j in range(1, mgncol + 1):
        dst[mgcols[j - 1] - 1] = src[j - 1]


@inline
def _micro_mg_data_scatter_unpack_2d(
    pcols: int,
    mgncol: int,
    top_lev: int,
    extent2: int,
    mgcols: Ptr[int],
    src: Ptr[float],
    dst: Ptr[float],
):
    for k in range(1, extent2 + 1):
        dst_k = top_lev + k - 1
        for j in range(1, mgncol + 1):
            dst[_idx2(mgcols[j - 1], dst_k, pcols)] = src[_idx2(j, k, mgncol)]


@inline
def _micro_mg_data_scatter_unpack_3d(
    pcols: int,
    pver: int,
    mgncol: int,
    nlev: int,
    top_lev: int,
    extent3: int,
    mgcols: Ptr[int],
    src: Ptr[float],
    dst: Ptr[float],
):
    for m in range(1, extent3 + 1):
        for k in range(1, nlev + 1):
            dst_k = top_lev + k - 1
            for j in range(1, mgncol + 1):
                dst[_idx3(mgcols[j - 1], dst_k, m, pcols, pver)] = src[
                    _idx3(j, k, m, mgncol, nlev)
                ]


@inline
def _micro_mg_data_accumulate_1d(n: int, src: Ptr[float], dst: Ptr[float]):
    for i in range(1, n + 1):
        dst[i - 1] = dst[i - 1] + src[i - 1]


@export
def micro_mg_init_codon(
    gravit: float,
    rair: float,
    rh2o: float,
    cpair: float,
    rhoh2o: float,
    tmelt_in: float,
    rhmini_in: float,
    micro_mg_berg_eff_factor_in: float,
    latvap: float,
    latice: float,
    micro_mg_dcs: float,
    scalars_p: cobj,
):
    scalars = Ptr[float](scalars_p)

    scalars[0] = gravit
    scalars[1] = rair
    scalars[2] = rh2o
    scalars[3] = cpair
    scalars[4] = rhoh2o
    scalars[5] = tmelt_in
    scalars[6] = rhmini_in
    scalars[7] = micro_mg_berg_eff_factor_in
    scalars[8] = latvap
    scalars[9] = latice
    scalars[10] = scalars[8] + scalars[9]
    scalars[11] = scalars[5]
    scalars[12] = scalars[5] - 5.0
    scalars[13] = 250.0
    scalars[14] = 500.0
    scalars[15] = 1000.0
    scalars[16] = 3.0e7
    scalars[17] = 2.0
    scalars[18] = 11.72
    scalars[19] = 0.41
    scalars[20] = 700.0
    scalars[21] = 1.0
    scalars[22] = 841.99667
    scalars[23] = 0.8
    scalars[24] = 3.1415927
    scalars[25] = scalars[14] * scalars[24] / 6.0
    scalars[26] = 3.0
    scalars[27] = scalars[13] * scalars[24] / 6.0
    scalars[28] = 3.0
    scalars[29] = scalars[15] * scalars[24] / 6.0
    scalars[30] = 3.0
    scalars[31] = 0.86
    scalars[32] = 0.28
    scalars[33] = 0.1
    scalars[34] = 1.0
    scalars[35] = 0.78
    scalars[36] = 0.32
    scalars[37] = micro_mg_dcs
    scalars[38] = 1.0e-18
    scalars[39] = 100.0
    scalars[40] = 0.66
    scalars[41] = 85000.0 / (rair * tmelt_in)
    scalars[42] = 0.1e-6
    scalars[43] = 273.15
    scalars[44] = -30.0
    scalars[45] = 26.0
    scalars[46] = -99.0
    scalars[47] = 1.26e-10
    scalars[48] = 1.0 / 10.0e-6
    scalars[49] = 1.0 / (2.0 * micro_mg_dcs)
    scalars[50] = 1.0 / 20.0e-6
    scalars[51] = 1.0 / 500.0e-6
    scalars[52] = 1.0 / 10.0e-6
    scalars[53] = 1.0 / 2000.0e-6
    scalars[54] = 4.0 / 3.0 * scalars[24] * scalars[14] * (10.0e-6) * (10.0e-6) * (10.0e-6)


@inline
def _micro_mg1_0_column_has_cloud(i: int, ldq: int, nlev: int, top_lev: int,
                                  qsmall: float, qcn: Ptr[float], qin: Ptr[float]) -> bool:
    lev_offset = top_lev - 1
    for k in range(top_lev, nlev + lev_offset + 1):
        idx = _idx2(i, k, ldq)
        if qcn[idx] >= qsmall:
            return True
    for k in range(top_lev, nlev + lev_offset + 1):
        idx = _idx2(i, k, ldq)
        if qin[idx] >= qsmall:
            return True
    return False


@export
def micro_mg1_0_get_cols_count_codon(
    ncol: int,
    ldq: int,
    nlev: int,
    top_lev: int,
    qsmall: float,
    qcn_p: cobj,
    qin_p: cobj,
    mgncol_p: cobj,
):
    qcn = Ptr[float](qcn_p)
    qin = Ptr[float](qin_p)
    mgncol = Ptr[int](mgncol_p)

    count = 0
    for i in range(1, ncol + 1):
        if _micro_mg1_0_column_has_cloud(i, ldq, nlev, top_lev, qsmall, qcn, qin):
            count += 1
    mgncol[0] = count


@export
def micro_mg1_0_get_cols_fill_codon(
    ncol: int,
    ldq: int,
    nlev: int,
    top_lev: int,
    qsmall: float,
    qcn_p: cobj,
    qin_p: cobj,
    mgcols_p: cobj,
):
    qcn = Ptr[float](qcn_p)
    qin = Ptr[float](qin_p)
    mgcols = Ptr[i32](mgcols_p)

    out_i = 0
    for i in range(1, ncol + 1):
        if _micro_mg1_0_column_has_cloud(i, ldq, nlev, top_lev, qsmall, qcn, qin):
            mgcols[out_i] = i32(i)
            out_i += 1


@export
def micro_mg_get_cols_codon(
    stage: int,
    ncol: int,
    ldq: int,
    nlev: int,
    top_lev: int,
    qsmall: float,
    qcn_p: cobj,
    qin_p: cobj,
    mgncol_p: cobj,
    mgcols_p: cobj,
):
    if stage == 0:
        micro_mg1_0_get_cols_count_codon(ncol, ldq, nlev, top_lev, qsmall, qcn_p, qin_p, mgncol_p)
    elif stage == 1:
        micro_mg1_0_get_cols_fill_codon(ncol, ldq, nlev, top_lev, qsmall, qcn_p, qin_p, mgcols_p)


@inline
def _micro_mg_data_accumulate_2d(
    n1: int,
    n2: int,
    src: Ptr[float],
    dst: Ptr[float],
):
    for k in range(1, n2 + 1):
        for i in range(1, n1 + 1):
            idx = _idx2(i, k, n1)
            dst[idx] = dst[idx] + src[idx]


@export
def mgfieldpostproc_accumulate_codon(
    accum_method: int,
    rank: int,
    n1: int,
    n2: int,
    plan_p: cobj,
    packed_p: cobj,
    buffer_p: cobj,
):
    plan = Ptr[int](plan_p)
    plan[1] = 0

    if accum_method == 0:
        return

    if accum_method == 1:
        plan[0] = plan[0] + 1
        if rank == 1:
            _micro_mg_data_accumulate_1d(
                n1, Ptr[float](packed_p), Ptr[float](buffer_p)
            )
        elif rank == 2:
            _micro_mg_data_accumulate_2d(
                n1, n2, Ptr[float](packed_p), Ptr[float](buffer_p)
            )
        else:
            plan[1] = 1
        return

    plan[1] = 2


@inline
def _micro_mg_data_mean_1d(
    n: int,
    num_steps: int,
    src: Ptr[float],
    dst: Ptr[float],
):
    for i in range(1, n + 1):
        dst[i - 1] = src[i - 1] / float(num_steps)


@inline
def _micro_mg_data_mean_2d(
    n1: int,
    n2: int,
    num_steps: int,
    src: Ptr[float],
    dst: Ptr[float],
):
    for k in range(1, n2 + 1):
        for i in range(1, n1 + 1):
            idx = _idx2(i, k, n1)
            dst[idx] = src[idx] / float(num_steps)


@export
def micro_mg_data_pack_unpack_codon(
    mode: int,
    pcols: int,
    pver: int,
    mgncol: int,
    nlev: int,
    top_lev: int,
    extent2: int,
    extent3: int,
    count1: int,
    num_steps: int,
    fillvalue: float,
    mgcols_p: cobj,
    src_p: cobj,
    fill_p: cobj,
    dst_p: cobj,
):
    mgcols = Ptr[int](mgcols_p)
    src = Ptr[float](src_p)
    dst = Ptr[float](dst_p)

    if mode == 1:
        _micro_mg_data_copy_pack_1d(mgncol, mgcols, src, dst)
    elif mode == 2 or mode == 3:
        _micro_mg_data_copy_pack_2d(
            pcols, mgncol, top_lev, extent2, mgcols, src, dst
        )
    elif mode == 4:
        _micro_mg_data_copy_pack_3d(
            pcols, pver, mgncol, nlev, top_lev, extent3, mgcols, src, dst
        )
    elif mode == 5:
        _micro_mg_data_fill_1d(pcols, fillvalue, dst)
        _micro_mg_data_scatter_unpack_1d(mgncol, mgcols, src, dst)
    elif mode == 6:
        _micro_mg_data_copy_1d(pcols, Ptr[float](fill_p), dst)
        _micro_mg_data_scatter_unpack_1d(mgncol, mgcols, src, dst)
    elif mode == 7:
        out_levs = pver + extent2 - nlev
        _micro_mg_data_fill_2d(pcols, out_levs, fillvalue, dst)
        _micro_mg_data_scatter_unpack_2d(
            pcols, mgncol, top_lev, extent2, mgcols, src, dst
        )
    elif mode == 8:
        out_levs = pver + extent2 - nlev
        _micro_mg_data_copy_2d(pcols, out_levs, Ptr[float](fill_p), dst)
        _micro_mg_data_scatter_unpack_2d(
            pcols, mgncol, top_lev, extent2, mgcols, src, dst
        )
    elif mode == 9:
        _micro_mg_data_fill_3d(pcols, pver, extent3, fillvalue, dst)
        _micro_mg_data_scatter_unpack_3d(
            pcols, pver, mgncol, nlev, top_lev, extent3, mgcols, src, dst
        )
    elif mode == 10:
        _micro_mg_data_copy_3d(pcols, pver, extent3, Ptr[float](fill_p), dst)
        _micro_mg_data_scatter_unpack_3d(
            pcols, pver, mgncol, nlev, top_lev, extent3, mgcols, src, dst
        )
    elif mode == 11:
        _micro_mg_data_accumulate_1d(count1, src, dst)
    elif mode == 12:
        _micro_mg_data_accumulate_2d(count1, extent2, src, dst)
    elif mode == 13:
        _micro_mg_data_mean_1d(count1, num_steps, src, dst)
    elif mode == 14:
        _micro_mg_data_mean_2d(count1, extent2, num_steps, src, dst)


@export
def pack_2d_codon(
    mode: int,
    pcols: int,
    pver: int,
    mgncol: int,
    nlev: int,
    top_lev: int,
    extent2: int,
    extent3: int,
    count1: int,
    num_steps: int,
    fillvalue: float,
    mgcols_p: cobj,
    src_p: cobj,
    fill_p: cobj,
    dst_p: cobj,
):
    micro_mg_data_pack_unpack_codon(
        mode, pcols, pver, mgncol, nlev, top_lev, extent2, extent3,
        count1, num_steps, fillvalue, mgcols_p, src_p, fill_p, dst_p
    )


@export
def pack_interface_codon(
    mode: int,
    pcols: int,
    pver: int,
    mgncol: int,
    nlev: int,
    top_lev: int,
    extent2: int,
    extent3: int,
    count1: int,
    num_steps: int,
    fillvalue: float,
    mgcols_p: cobj,
    src_p: cobj,
    fill_p: cobj,
    dst_p: cobj,
):
    micro_mg_data_pack_unpack_codon(
        mode, pcols, pver, mgncol, nlev, top_lev, extent2, extent3,
        count1, num_steps, fillvalue, mgcols_p, src_p, fill_p, dst_p
    )


@export
def pack_3d_codon(
    mode: int,
    pcols: int,
    pver: int,
    mgncol: int,
    nlev: int,
    top_lev: int,
    extent2: int,
    extent3: int,
    count1: int,
    num_steps: int,
    fillvalue: float,
    mgcols_p: cobj,
    src_p: cobj,
    fill_p: cobj,
    dst_p: cobj,
):
    micro_mg_data_pack_unpack_codon(
        mode, pcols, pver, mgncol, nlev, top_lev, extent2, extent3,
        count1, num_steps, fillvalue, mgcols_p, src_p, fill_p, dst_p
    )


@export
def unpack_1d_codon(
    mode: int,
    pcols: int,
    pver: int,
    mgncol: int,
    nlev: int,
    top_lev: int,
    extent2: int,
    extent3: int,
    count1: int,
    num_steps: int,
    fillvalue: float,
    mgcols_p: cobj,
    src_p: cobj,
    fill_p: cobj,
    dst_p: cobj,
):
    micro_mg_data_pack_unpack_codon(
        mode, pcols, pver, mgncol, nlev, top_lev, extent2, extent3,
        count1, num_steps, fillvalue, mgcols_p, src_p, fill_p, dst_p
    )


@export
def unpack_2d_codon(
    mode: int,
    pcols: int,
    pver: int,
    mgncol: int,
    nlev: int,
    top_lev: int,
    extent2: int,
    extent3: int,
    count1: int,
    num_steps: int,
    fillvalue: float,
    mgcols_p: cobj,
    src_p: cobj,
    fill_p: cobj,
    dst_p: cobj,
):
    micro_mg_data_pack_unpack_codon(
        mode, pcols, pver, mgncol, nlev, top_lev, extent2, extent3,
        count1, num_steps, fillvalue, mgcols_p, src_p, fill_p, dst_p
    )


@export
def unpack_2d_array_fill_codon(
    mode: int,
    pcols: int,
    pver: int,
    mgncol: int,
    nlev: int,
    top_lev: int,
    extent2: int,
    extent3: int,
    count1: int,
    num_steps: int,
    fillvalue: float,
    mgcols_p: cobj,
    src_p: cobj,
    fill_p: cobj,
    dst_p: cobj,
):
    micro_mg_data_pack_unpack_codon(
        mode, pcols, pver, mgncol, nlev, top_lev, extent2, extent3,
        count1, num_steps, fillvalue, mgcols_p, src_p, fill_p, dst_p
    )


@inline
def _zero2(ncol: int, pver: int, pcols: int, arr: Ptr[float]):
    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            arr[_idx2(i, k, pcols)] = 0.0


@inline
def _zero2_interface(ncol: int, pverp: int, pcols: int, arr: Ptr[float]):
    for k in range(1, pverp + 1):
        for i in range(1, ncol + 1):
            arr[_idx2(i, k, pcols)] = 0.0


@inline
def _zero2_full_pcols(pver: int, pcols: int, arr: Ptr[float]):
    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            arr[_idx2(i, k, pcols)] = 0.0


@inline
def _zero2_interface_full_pcols(pverp: int, pcols: int, arr: Ptr[float]):
    for k in range(1, pverp + 1):
        for i in range(1, pcols + 1):
            arr[_idx2(i, k, pcols)] = 0.0


@export
def micro_mg1_0_init_fields_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    mincld: float,
    qn_p: cobj,
    tn_p: cobj,
    qc_p: cobj,
    qi_p: cobj,
    nc_p: cobj,
    ni_p: cobj,
    ncai_p: cobj,
    ncal_p: cobj,
    rercld_p: cobj,
    arcld_p: cobj,
    pgamrad_p: cobj,
    lamcrad_p: cobj,
    deffi_p: cobj,
    qcsevap_p: cobj,
    qisevap_p: cobj,
    qvres_p: cobj,
    cmeiout_p: cobj,
    vtrmc_p: cobj,
    vtrmi_p: cobj,
    qcsedten_p: cobj,
    qisedten_p: cobj,
    prao_p: cobj,
    prco_p: cobj,
    mnuccco_p: cobj,
    mnuccto_p: cobj,
    msacwio_p: cobj,
    psacwso_p: cobj,
    bergso_p: cobj,
    bergo_p: cobj,
    melto_p: cobj,
    homoo_p: cobj,
    qcreso_p: cobj,
    prcio_p: cobj,
    praio_p: cobj,
    qireso_p: cobj,
    mnuccro_p: cobj,
    pracso_p: cobj,
    meltsdt_p: cobj,
    frzrdt_p: cobj,
    mnuccdo_p: cobj,
    rflx_p: cobj,
    sflx_p: cobj,
    effc_p: cobj,
    effc_fn_p: cobj,
    effi_p: cobj,
    preo_p: cobj,
    prdso_p: cobj,
    frzro_p: cobj,
    meltso_p: cobj,
    wtfc_p: cobj,
    wtfi_p: cobj,
    wtprelat_p: cobj,
    wtpostlat_p: cobj,
    q_p: cobj,
    t_p: cobj,
    t1_p: cobj,
    q1_p: cobj,
    qc1_p: cobj,
    qi1_p: cobj,
    nc1_p: cobj,
    ni1_p: cobj,
    tlat1_p: cobj,
    qvlat1_p: cobj,
    qctend1_p: cobj,
    qitend1_p: cobj,
    nctend1_p: cobj,
    nitend1_p: cobj,
    qrout_p: cobj,
    qsout_p: cobj,
    nrout_p: cobj,
    nsout_p: cobj,
    dsout_p: cobj,
    drout_p: cobj,
    reff_rain_p: cobj,
    reff_snow_p: cobj,
    nevapr_p: cobj,
    nevapr2_p: cobj,
    evapsnow_p: cobj,
    prain_p: cobj,
    prodsnow_p: cobj,
    cmeout_p: cobj,
    am_evp_st_p: cobj,
    rainrt1_p: cobj,
    cldmax_p: cobj,
    dum2l_p: cobj,
    dum2i_p: cobj,
    prect1_p: cobj,
    preci1_p: cobj,
):
    qn = Ptr[float](qn_p)
    tn = Ptr[float](tn_p)
    qc = Ptr[float](qc_p)
    qi = Ptr[float](qi_p)
    nc = Ptr[float](nc_p)
    ni = Ptr[float](ni_p)
    q = Ptr[float](q_p)
    t = Ptr[float](t_p)
    t1 = Ptr[float](t1_p)
    q1 = Ptr[float](q1_p)
    qc1 = Ptr[float](qc1_p)
    qi1 = Ptr[float](qi1_p)
    nc1 = Ptr[float](nc1_p)
    ni1 = Ptr[float](ni1_p)

    _zero2(ncol, pver, pcols, Ptr[float](ncai_p))
    _zero2(ncol, pver, pcols, Ptr[float](ncal_p))
    _zero2(ncol, pver, pcols, Ptr[float](rercld_p))
    _zero2(ncol, pver, pcols, Ptr[float](arcld_p))
    _zero2(ncol, pver, pcols, Ptr[float](pgamrad_p))
    _zero2(ncol, pver, pcols, Ptr[float](lamcrad_p))
    _zero2(ncol, pver, pcols, Ptr[float](deffi_p))
    _zero2(ncol, pver, pcols, Ptr[float](qcsevap_p))
    _zero2(ncol, pver, pcols, Ptr[float](qisevap_p))
    _zero2(ncol, pver, pcols, Ptr[float](qvres_p))
    _zero2(ncol, pver, pcols, Ptr[float](cmeiout_p))
    _zero2(ncol, pver, pcols, Ptr[float](vtrmc_p))
    _zero2(ncol, pver, pcols, Ptr[float](vtrmi_p))
    _zero2(ncol, pver, pcols, Ptr[float](qcsedten_p))
    _zero2(ncol, pver, pcols, Ptr[float](qisedten_p))
    _zero2(ncol, pver, pcols, Ptr[float](prao_p))
    _zero2(ncol, pver, pcols, Ptr[float](prco_p))
    _zero2(ncol, pver, pcols, Ptr[float](mnuccco_p))
    _zero2(ncol, pver, pcols, Ptr[float](mnuccto_p))
    _zero2(ncol, pver, pcols, Ptr[float](msacwio_p))
    _zero2(ncol, pver, pcols, Ptr[float](psacwso_p))
    _zero2(ncol, pver, pcols, Ptr[float](bergso_p))
    _zero2(ncol, pver, pcols, Ptr[float](bergo_p))
    _zero2(ncol, pver, pcols, Ptr[float](melto_p))
    _zero2(ncol, pver, pcols, Ptr[float](homoo_p))
    _zero2(ncol, pver, pcols, Ptr[float](qcreso_p))
    _zero2(ncol, pver, pcols, Ptr[float](prcio_p))
    _zero2(ncol, pver, pcols, Ptr[float](praio_p))
    _zero2(ncol, pver, pcols, Ptr[float](qireso_p))
    _zero2(ncol, pver, pcols, Ptr[float](mnuccro_p))
    _zero2(ncol, pver, pcols, Ptr[float](pracso_p))
    _zero2(ncol, pver, pcols, Ptr[float](meltsdt_p))
    _zero2(ncol, pver, pcols, Ptr[float](frzrdt_p))
    _zero2(ncol, pver, pcols, Ptr[float](mnuccdo_p))
    _zero2_interface_full_pcols(pver + 1, pcols, Ptr[float](rflx_p))
    _zero2_interface_full_pcols(pver + 1, pcols, Ptr[float](sflx_p))
    _zero2_full_pcols(pver, pcols, Ptr[float](effc_p))
    _zero2_full_pcols(pver, pcols, Ptr[float](effc_fn_p))
    _zero2_full_pcols(pver, pcols, Ptr[float](effi_p))
    _zero2(ncol, pver, pcols, Ptr[float](preo_p))
    _zero2(ncol, pver, pcols, Ptr[float](prdso_p))
    _zero2(ncol, pver, pcols, Ptr[float](frzro_p))
    _zero2(ncol, pver, pcols, Ptr[float](meltso_p))
    _zero2(ncol, pver, pcols, Ptr[float](wtfc_p))
    _zero2(ncol, pver, pcols, Ptr[float](wtfi_p))
    _zero2(ncol, pver, pcols, Ptr[float](wtprelat_p))
    _zero2(ncol, pver, pcols, Ptr[float](wtpostlat_p))

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            q[idx] = qn[idx]
            t[idx] = tn[idx]

    for k in range(1, top_lev):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            qc[idx] = 0.0
            qi[idx] = 0.0
            nc[idx] = 0.0
            ni[idx] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            t1[idx] = t[idx]
            q1[idx] = q[idx]
            qc1[idx] = qc[idx]
            qi1[idx] = qi[idx]
            nc1[idx] = nc[idx]
            ni1[idx] = ni[idx]

    _zero2(ncol, pver, pcols, Ptr[float](tlat1_p))
    _zero2(ncol, pver, pcols, Ptr[float](qvlat1_p))
    _zero2(ncol, pver, pcols, Ptr[float](qctend1_p))
    _zero2(ncol, pver, pcols, Ptr[float](qitend1_p))
    _zero2(ncol, pver, pcols, Ptr[float](nctend1_p))
    _zero2(ncol, pver, pcols, Ptr[float](nitend1_p))
    _zero2(ncol, pver, pcols, Ptr[float](qrout_p))
    _zero2(ncol, pver, pcols, Ptr[float](qsout_p))
    _zero2(ncol, pver, pcols, Ptr[float](nrout_p))
    _zero2(ncol, pver, pcols, Ptr[float](nsout_p))
    _zero2(ncol, pver, pcols, Ptr[float](dsout_p))
    _zero2(ncol, pver, pcols, Ptr[float](drout_p))
    _zero2(ncol, pver, pcols, Ptr[float](reff_rain_p))
    _zero2(ncol, pver, pcols, Ptr[float](reff_snow_p))
    _zero2(ncol, pver, pcols, Ptr[float](nevapr_p))
    _zero2(ncol, pver, pcols, Ptr[float](nevapr2_p))
    _zero2(ncol, pver, pcols, Ptr[float](evapsnow_p))
    _zero2(ncol, pver, pcols, Ptr[float](prain_p))
    _zero2(ncol, pver, pcols, Ptr[float](prodsnow_p))
    _zero2(ncol, pver, pcols, Ptr[float](cmeout_p))
    _zero2(ncol, pver, pcols, Ptr[float](am_evp_st_p))
    _zero2(ncol, pver, pcols, Ptr[float](rainrt1_p))

    cldmax = Ptr[float](cldmax_p)
    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            cldmax[_idx2(i, k, pcols)] = mincld

    _zero2(ncol, pver, pcols, Ptr[float](dum2l_p))
    _zero2(ncol, pver, pcols, Ptr[float](dum2i_p))

    prect1 = Ptr[float](prect1_p)
    preci1 = Ptr[float](preci1_p)
    for i in range(1, ncol + 1):
        prect1[i - 1] = 0.0
        preci1[i - 1] = 0.0


@inline
def _zero2_column(i: int, pver: int, pcols: int, arr: Ptr[float]):
    for k in range(1, pver + 1):
        arr[_idx2(i, k, pcols)] = 0.0


@inline
def _micro_mg1_0_zero_tendency_column(i: int, pver: int, pcols: int,
                                      tlat: Ptr[float], qvlat: Ptr[float],
                                      qctend: Ptr[float], qitend: Ptr[float],
                                      qnitend: Ptr[float], qrtend: Ptr[float],
                                      nctend: Ptr[float], nitend: Ptr[float],
                                      nrtend: Ptr[float], nstend: Ptr[float]):
    _zero2_column(i, pver, pcols, tlat)
    _zero2_column(i, pver, pcols, qvlat)
    _zero2_column(i, pver, pcols, qctend)
    _zero2_column(i, pver, pcols, qitend)
    _zero2_column(i, pver, pcols, qnitend)
    _zero2_column(i, pver, pcols, qrtend)
    _zero2_column(i, pver, pcols, nctend)
    _zero2_column(i, pver, pcols, nitend)
    _zero2_column(i, pver, pcols, nrtend)
    _zero2_column(i, pver, pcols, nstend)


@inline
def _micro_mg1_0_zero_precip_diag_column(i: int, pver: int, pcols: int,
                                         qniic: Ptr[float], qric: Ptr[float],
                                         nsic: Ptr[float], nric: Ptr[float],
                                         rainrt: Ptr[float],
                                         qrtend_copy: Ptr[float],
                                         qnitend_copy: Ptr[float]):
    _zero2_column(i, pver, pcols, qniic)
    _zero2_column(i, pver, pcols, qric)
    _zero2_column(i, pver, pcols, nsic)
    _zero2_column(i, pver, pcols, nric)
    _zero2_column(i, pver, pcols, rainrt)
    _zero2_column(i, pver, pcols, qrtend_copy)
    _zero2_column(i, pver, pcols, qnitend_copy)


@export
def micro_mg1_0_no_cloud_zero_column_codon(
    i: int,
    pcols: int,
    pver: int,
    tlat_p: cobj,
    qvlat_p: cobj,
    qctend_p: cobj,
    qitend_p: cobj,
    qnitend_p: cobj,
    qrtend_p: cobj,
    nctend_p: cobj,
    nitend_p: cobj,
    nrtend_p: cobj,
    nstend_p: cobj,
    prect_p: cobj,
    preci_p: cobj,
    qniic_p: cobj,
    qric_p: cobj,
    nsic_p: cobj,
    nric_p: cobj,
    rainrt_p: cobj,
    qrtend_copy_p: cobj,
    qnitend_copy_p: cobj,
):
    _micro_mg1_0_zero_tendency_column(
        i, pver, pcols,
        Ptr[float](tlat_p), Ptr[float](qvlat_p), Ptr[float](qctend_p),
        Ptr[float](qitend_p), Ptr[float](qnitend_p), Ptr[float](qrtend_p),
        Ptr[float](nctend_p), Ptr[float](nitend_p), Ptr[float](nrtend_p),
        Ptr[float](nstend_p),
    )
    prect = Ptr[float](prect_p)
    preci = Ptr[float](preci_p)
    prect[i - 1] = 0.0
    preci[i - 1] = 0.0
    _micro_mg1_0_zero_precip_diag_column(
        i, pver, pcols,
        Ptr[float](qniic_p), Ptr[float](qric_p), Ptr[float](nsic_p),
        Ptr[float](nric_p), Ptr[float](rainrt_p),
        Ptr[float](qrtend_copy_p), Ptr[float](qnitend_copy_p),
    )


@export
def micro_mg1_0_substep_zero_column_codon(
    i: int,
    pcols: int,
    pver: int,
    tlat_p: cobj,
    qvlat_p: cobj,
    qctend_p: cobj,
    qitend_p: cobj,
    qnitend_p: cobj,
    qrtend_p: cobj,
    nctend_p: cobj,
    nitend_p: cobj,
    nrtend_p: cobj,
    nstend_p: cobj,
    qniic_p: cobj,
    qric_p: cobj,
    nsic_p: cobj,
    nric_p: cobj,
    rainrt_p: cobj,
    qrtend_copy_p: cobj,
    qnitend_copy_p: cobj,
):
    _micro_mg1_0_zero_tendency_column(
        i, pver, pcols,
        Ptr[float](tlat_p), Ptr[float](qvlat_p), Ptr[float](qctend_p),
        Ptr[float](qitend_p), Ptr[float](qnitend_p), Ptr[float](qrtend_p),
        Ptr[float](nctend_p), Ptr[float](nitend_p), Ptr[float](nrtend_p),
        Ptr[float](nstend_p),
    )
    _micro_mg1_0_zero_precip_diag_column(
        i, pver, pcols,
        Ptr[float](qniic_p), Ptr[float](qric_p), Ptr[float](nsic_p),
        Ptr[float](nric_p), Ptr[float](rainrt_p),
        Ptr[float](qrtend_copy_p), Ptr[float](qnitend_copy_p),
    )


@export
def micro_mg1_0_substep_setup_column_codon(
    i: int,
    pcols: int,
    pver: int,
    top_lev: int,
    qsmall: float,
    mincld: float,
    qc_p: cobj,
    qi_p: cobj,
    ni_p: cobj,
    cldm_p: cobj,
    cmei_p: cobj,
    cwml_p: cobj,
    cwmi_p: cobj,
    ums_p: cobj,
    uns_p: cobj,
    umr_p: cobj,
    unr_p: cobj,
    nsubi_p: cobj,
    nsubc_p: cobj,
):
    qc = Ptr[float](qc_p)
    qi = Ptr[float](qi_p)
    ni = Ptr[float](ni_p)
    cldm = Ptr[float](cldm_p)
    cmei = Ptr[float](cmei_p)
    cwml = Ptr[float](cwml_p)
    cwmi = Ptr[float](cwmi_p)
    ums = Ptr[float](ums_p)
    uns = Ptr[float](uns_p)
    umr = Ptr[float](umr_p)
    unr = Ptr[float](unr_p)
    nsubi = Ptr[float](nsubi_p)
    nsubc = Ptr[float](nsubc_p)

    for k in range(top_lev, pver + 1):
        idx = _idx2(i, k, pcols)
        cwml[idx] = qc[idx]
        cwmi[idx] = qi[idx]

        ums[k - 1] = 0.0
        uns[k - 1] = 0.0
        umr[k - 1] = 0.0
        unr[k - 1] = 0.0

        if cmei[idx] < 0.0 and qi[idx] > qsmall and cldm[idx] > mincld:
            nsubi[k - 1] = cmei[idx] / qi[idx] * ni[idx] / cldm[idx]
        else:
            nsubi[k - 1] = 0.0
        nsubc[k - 1] = 0.0


@export
def micro_mg1_0_incloud_activation_prep_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    deltat: float,
    qsmall: float,
    omsm: float,
    cdnl: float,
    do_cldice: int,
    cwml_p: cobj,
    cwmi_p: cobj,
    lcldm_p: cobj,
    icldm_p: cobj,
    nc_p: cobj,
    ni_p: cobj,
    qc_p: cobj,
    qi_p: cobj,
    berg_p: cobj,
    cmei_p: cobj,
    cmeout_p: cobj,
    npccnin_p: cobj,
    rho_p: cobj,
    qcic_p: cobj,
    qiic_p: cobj,
    ncic_p: cobj,
    niic_p: cobj,
    npccn_p: cobj,
    dum2l_p: cobj,
) -> float:
    cwml = Ptr[float](cwml_p)
    cwmi = Ptr[float](cwmi_p)
    lcldm = Ptr[float](lcldm_p)
    icldm = Ptr[float](icldm_p)
    nc = Ptr[float](nc_p)
    ni = Ptr[float](ni_p)
    qc = Ptr[float](qc_p)
    qi = Ptr[float](qi_p)
    berg = Ptr[float](berg_p)
    cmei = Ptr[float](cmei_p)
    cmeout = Ptr[float](cmeout_p)
    npccnin = Ptr[float](npccnin_p)
    rho = Ptr[float](rho_p)
    qcic = Ptr[float](qcic_p)
    qiic = Ptr[float](qiic_p)
    ncic = Ptr[float](ncic_p)
    niic = Ptr[float](niic_p)
    npccn = Ptr[float](npccn_p)
    dum2l = Ptr[float](dum2l_p)

    idx = _idx2(i, k, pcols)
    qcic[idx] = min(cwml[idx] / lcldm[idx], 5.0e-3)
    qiic[idx] = min(cwmi[idx] / icldm[idx], 5.0e-3)
    ncic[idx] = max(nc[idx] / lcldm[idx], 0.0)
    niic[idx] = max(ni[idx] / icldm[idx], 0.0)

    if qc[idx] - berg[idx] * deltat < qsmall:
        qcic[idx] = 0.0
        ncic[idx] = 0.0
        if qc[idx] - berg[idx] * deltat < 0.0:
            berg[idx] = qc[idx] / deltat * omsm

    if do_cldice != 0 and qi[idx] + (cmei[idx] + berg[idx]) * deltat < qsmall:
        qiic[idx] = 0.0
        niic[idx] = 0.0
        if qi[idx] + (cmei[idx] + berg[idx]) * deltat < 0.0:
            cmei[idx] = (-qi[idx] / deltat - berg[idx]) * omsm

    cmeout[idx] = cmeout[idx] + cmei[idx]

    ncmax = 0.0
    if qcic[idx] >= qsmall:
        npccn[k - 1] = max(0.0, npccnin[idx])
        dum2l[idx] = (nc[idx] + npccn[k - 1] * deltat) / lcldm[idx]
        dum2l[idx] = max(dum2l[idx], cdnl / rho[idx])
        ncmax = dum2l[idx] * lcldm[idx]
    else:
        npccn[k - 1] = 0.0
        dum2l[idx] = 0.0
        ncmax = 0.0

    return ncmax


@export
def micro_mg1_0_conservation_limiter_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    deltat: float,
    omsm: float,
    qsmall: float,
    qce: float,
    nce: float,
    qie: float,
    nie: float,
    qrtot: float,
    nrtot: float,
    qstot: float,
    nstot: float,
    do_cldice: int,
    use_hetfrz_classnuc: int,
    lcldm_p: cobj,
    icldm_p: cobj,
    cldmax_p: cobj,
    dz_p: cobj,
    rho_p: cobj,
    prc_p: cobj,
    pra_p: cobj,
    mnuccc_p: cobj,
    mnucct_p: cobj,
    msacwi_p: cobj,
    psacws_p: cobj,
    bergs_p: cobj,
    nprc1_p: cobj,
    npra_p: cobj,
    nnuccc_p: cobj,
    nnucct_p: cobj,
    npsacws_p: cobj,
    nsubc_p: cobj,
    prci_p: cobj,
    prai_p: cobj,
    mnudep_p: cobj,
    nprci_p: cobj,
    nprai_p: cobj,
    nsubi_p: cobj,
    nnudep_p: cobj,
    nsacwi_p: cobj,
    mnuccr_p: cobj,
    pre_p: cobj,
    pracs_p: cobj,
    nsubr_p: cobj,
    npracs_p: cobj,
    nnuccr_p: cobj,
    nragg_p: cobj,
    prds_p: cobj,
    nsubs_p: cobj,
    nsagg_p: cobj,
    nprc_p: cobj,
):
    lcldm = Ptr[float](lcldm_p)
    icldm = Ptr[float](icldm_p)
    cldmax = Ptr[float](cldmax_p)
    dz = Ptr[float](dz_p)
    rho = Ptr[float](rho_p)
    prc = Ptr[float](prc_p)
    pra = Ptr[float](pra_p)
    mnuccc = Ptr[float](mnuccc_p)
    mnucct = Ptr[float](mnucct_p)
    msacwi = Ptr[float](msacwi_p)
    psacws = Ptr[float](psacws_p)
    bergs = Ptr[float](bergs_p)
    nprc1 = Ptr[float](nprc1_p)
    npra = Ptr[float](npra_p)
    nnuccc = Ptr[float](nnuccc_p)
    nnucct = Ptr[float](nnucct_p)
    npsacws = Ptr[float](npsacws_p)
    nsubc = Ptr[float](nsubc_p)
    prci = Ptr[float](prci_p)
    prai = Ptr[float](prai_p)
    mnudep = Ptr[float](mnudep_p)
    nprci = Ptr[float](nprci_p)
    nprai = Ptr[float](nprai_p)
    nsubi = Ptr[float](nsubi_p)
    nnudep = Ptr[float](nnudep_p)
    nsacwi = Ptr[float](nsacwi_p)
    mnuccr = Ptr[float](mnuccr_p)
    pre = Ptr[float](pre_p)
    pracs = Ptr[float](pracs_p)
    nsubr = Ptr[float](nsubr_p)
    npracs = Ptr[float](npracs_p)
    nnuccr = Ptr[float](nnuccr_p)
    nragg = Ptr[float](nragg_p)
    prds = Ptr[float](prds_p)
    nsubs = Ptr[float](nsubs_p)
    nsagg = Ptr[float](nsagg_p)
    nprc = Ptr[float](nprc_p)

    idx = _idx2(i, k, pcols)
    kk = k - 1

    dum = (
        prc[kk]
        + pra[kk]
        + mnuccc[kk]
        + mnucct[kk]
        + msacwi[kk]
        + psacws[kk]
        + bergs[kk]
    ) * lcldm[idx] * deltat

    if dum > qce:
        ratio = (
            qce
            / deltat
            / lcldm[idx]
            / (prc[kk] + pra[kk] + mnuccc[kk] + mnucct[kk] + msacwi[kk] + psacws[kk] + bergs[kk])
            * omsm
        )
        prc[kk] = prc[kk] * ratio
        pra[kk] = pra[kk] * ratio
        mnuccc[kk] = mnuccc[kk] * ratio
        mnucct[kk] = mnucct[kk] * ratio
        msacwi[kk] = msacwi[kk] * ratio
        psacws[kk] = psacws[kk] * ratio
        bergs[kk] = bergs[kk] * ratio

    dum = (nprc1[kk] + npra[kk] + nnuccc[kk] + nnucct[kk] + npsacws[kk] - nsubc[kk]) * lcldm[idx] * deltat

    if dum > nce:
        ratio = nce / deltat / ((nprc1[kk] + npra[kk] + nnuccc[kk] + nnucct[kk] + npsacws[kk] - nsubc[kk]) * lcldm[idx]) * omsm
        nprc1[kk] = nprc1[kk] * ratio
        npra[kk] = npra[kk] * ratio
        nnuccc[kk] = nnuccc[kk] * ratio
        nnucct[kk] = nnucct[kk] * ratio
        npsacws[kk] = npsacws[kk] * ratio
        nsubc[kk] = nsubc[kk] * ratio

    if do_cldice != 0:
        frztmp = -mnuccc[kk] - mnucct[kk] - msacwi[kk]
        if use_hetfrz_classnuc != 0:
            frztmp = -mnuccc[kk] - mnucct[kk] - mnudep[kk] - msacwi[kk]
        dum = (frztmp * lcldm[idx] + (prci[kk] + prai[kk]) * icldm[idx]) * deltat

        if dum > qie:
            frztmp = mnuccc[kk] + mnucct[kk] + msacwi[kk]
            if use_hetfrz_classnuc != 0:
                frztmp = mnuccc[kk] + mnucct[kk] + mnudep[kk] + msacwi[kk]
            ratio = (qie / deltat + frztmp * lcldm[idx]) / ((prci[kk] + prai[kk]) * icldm[idx]) * omsm
            prci[kk] = prci[kk] * ratio
            prai[kk] = prai[kk] * ratio

        frztmp = -nnucct[kk] - nsacwi[kk]
        if use_hetfrz_classnuc != 0:
            frztmp = -nnucct[kk] - nnuccc[kk] - nnudep[kk] - nsacwi[kk]
        dum = (frztmp * lcldm[idx] + (nprci[kk] + nprai[kk] - nsubi[kk]) * icldm[idx]) * deltat

        if dum > nie:
            frztmp = nnucct[kk] + nsacwi[kk]
            if use_hetfrz_classnuc != 0:
                frztmp = nnucct[kk] + nnuccc[kk] + nnudep[kk] + nsacwi[kk]
            ratio = (nie / deltat + frztmp * lcldm[idx]) / ((nprci[kk] + nprai[kk] - nsubi[kk]) * icldm[idx]) * omsm
            nprci[kk] = nprci[kk] * ratio
            nprai[kk] = nprai[kk] * ratio
            nsubi[kk] = nsubi[kk] * ratio

    if ((prc[kk] + pra[kk]) * lcldm[idx] + (-mnuccr[kk] + pre[kk] - pracs[kk]) * cldmax[idx]) * dz[idx] * rho[idx] + qrtot < 0.0:
        if -pre[kk] + pracs[kk] + mnuccr[kk] >= qsmall:
            ratio = (qrtot / (dz[idx] * rho[idx]) + (prc[kk] + pra[kk]) * lcldm[idx]) / ((-pre[kk] + pracs[kk] + mnuccr[kk]) * cldmax[idx]) * omsm
            pre[kk] = pre[kk] * ratio
            pracs[kk] = pracs[kk] * ratio
            mnuccr[kk] = mnuccr[kk] * ratio

    nsubr[kk] = 0.0

    if (nprc[kk] * lcldm[idx] + (-nnuccr[kk] + nsubr[kk] - npracs[kk] + nragg[kk]) * cldmax[idx]) * dz[idx] * rho[idx] + nrtot < 0.0:
        if -nsubr[kk] - nragg[kk] + npracs[kk] + nnuccr[kk] >= qsmall:
            ratio = (nrtot / (dz[idx] * rho[idx]) + nprc[kk] * lcldm[idx]) / ((-nsubr[kk] - nragg[kk] + npracs[kk] + nnuccr[kk]) * cldmax[idx]) * omsm
            nsubr[kk] = nsubr[kk] * ratio
            npracs[kk] = npracs[kk] * ratio
            nnuccr[kk] = nnuccr[kk] * ratio
            nragg[kk] = nragg[kk] * ratio

    if (
        (bergs[kk] + psacws[kk]) * lcldm[idx]
        + (prai[kk] + prci[kk]) * icldm[idx]
        + (pracs[kk] + mnuccr[kk] + prds[kk]) * cldmax[idx]
    ) * dz[idx] * rho[idx] + qstot < 0.0:
        if -prds[kk] >= qsmall:
            ratio = (
                qstot / (dz[idx] * rho[idx])
                + (bergs[kk] + psacws[kk]) * lcldm[idx]
                + (prai[kk] + prci[kk]) * icldm[idx]
                + (pracs[kk] + mnuccr[kk]) * cldmax[idx]
            ) / (-prds[kk] * cldmax[idx]) * omsm
            prds[kk] = prds[kk] * ratio

    nsubs[kk] = 0.0

    if (nprci[kk] * icldm[idx] + (nnuccr[kk] + nsubs[kk] + nsagg[kk]) * cldmax[idx]) * dz[idx] * rho[idx] + nstot < 0.0:
        if -nsubs[kk] - nsagg[kk] >= qsmall:
            ratio = (nstot / (dz[idx] * rho[idx]) + nprci[kk] * icldm[idx] + nnuccr[kk] * cldmax[idx]) / ((-nsubs[kk] - nsagg[kk]) * cldmax[idx]) * omsm
            nsubs[kk] = nsubs[kk] * ratio
            nsagg[kk] = nsagg[kk] * ratio


@export
def micro_mg1_0_process_output_accum_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    qrtend_p: cobj,
    qnitend_p: cobj,
    qrtend_copy_p: cobj,
    qnitend_copy_p: cobj,
    cmei_p: cobj,
    cmeiout_p: cobj,
    prds_p: cobj,
    pre_p: cobj,
    cldmax_p: cobj,
    evapsnow_p: cobj,
    nevapr_p: cobj,
    nevapr2_p: cobj,
    pra_p: cobj,
    prc_p: cobj,
    lcldm_p: cobj,
    pracs_p: cobj,
    mnuccr_p: cobj,
    prain_p: cobj,
    prai_p: cobj,
    prci_p: cobj,
    icldm_p: cobj,
    psacws_p: cobj,
    bergs_p: cobj,
    prodsnow_p: cobj,
    qcsinksum_rate1ord_p: cobj,
    qcsum_rate1ord_p: cobj,
    qc_p: cobj,
    prao_p: cobj,
    prco_p: cobj,
    mnuccc_p: cobj,
    mnucct_p: cobj,
    mnuccd_p: cobj,
    msacwi_p: cobj,
    mnuccco_p: cobj,
    mnuccto_p: cobj,
    mnuccdo_p: cobj,
    msacwio_p: cobj,
    psacwso_p: cobj,
    bergso_p: cobj,
    berg_p: cobj,
    bergo_p: cobj,
    prcio_p: cobj,
    praio_p: cobj,
    mnuccro_p: cobj,
    pracso_p: cobj,
    preo_p: cobj,
    prdso_p: cobj,
):
    qrtend = Ptr[float](qrtend_p)
    qnitend = Ptr[float](qnitend_p)
    qrtend_copy = Ptr[float](qrtend_copy_p)
    qnitend_copy = Ptr[float](qnitend_copy_p)
    cmei = Ptr[float](cmei_p)
    cmeiout = Ptr[float](cmeiout_p)
    prds = Ptr[float](prds_p)
    pre = Ptr[float](pre_p)
    cldmax = Ptr[float](cldmax_p)
    evapsnow = Ptr[float](evapsnow_p)
    nevapr = Ptr[float](nevapr_p)
    nevapr2 = Ptr[float](nevapr2_p)
    pra = Ptr[float](pra_p)
    prc = Ptr[float](prc_p)
    lcldm = Ptr[float](lcldm_p)
    pracs = Ptr[float](pracs_p)
    mnuccr = Ptr[float](mnuccr_p)
    prain = Ptr[float](prain_p)
    prai = Ptr[float](prai_p)
    prci = Ptr[float](prci_p)
    icldm = Ptr[float](icldm_p)
    psacws = Ptr[float](psacws_p)
    bergs = Ptr[float](bergs_p)
    prodsnow = Ptr[float](prodsnow_p)
    qcsinksum_rate1ord = Ptr[float](qcsinksum_rate1ord_p)
    qcsum_rate1ord = Ptr[float](qcsum_rate1ord_p)
    qc = Ptr[float](qc_p)
    prao = Ptr[float](prao_p)
    prco = Ptr[float](prco_p)
    mnuccc = Ptr[float](mnuccc_p)
    mnucct = Ptr[float](mnucct_p)
    mnuccd = Ptr[float](mnuccd_p)
    msacwi = Ptr[float](msacwi_p)
    mnuccco = Ptr[float](mnuccco_p)
    mnuccto = Ptr[float](mnuccto_p)
    mnuccdo = Ptr[float](mnuccdo_p)
    msacwio = Ptr[float](msacwio_p)
    psacwso = Ptr[float](psacwso_p)
    bergso = Ptr[float](bergso_p)
    berg = Ptr[float](berg_p)
    bergo = Ptr[float](bergo_p)
    prcio = Ptr[float](prcio_p)
    praio = Ptr[float](praio_p)
    mnuccro = Ptr[float](mnuccro_p)
    pracso = Ptr[float](pracso_p)
    preo = Ptr[float](preo_p)
    prdso = Ptr[float](prdso_p)

    idx = _idx2(i, k, pcols)
    kk = k - 1

    qrtend_copy[idx] = qrtend_copy[idx] + qrtend[idx]
    qnitend_copy[idx] = qnitend_copy[idx] + qnitend[idx]

    cmeiout[idx] = cmeiout[idx] + cmei[idx]

    evapsnow[idx] = evapsnow[idx] - prds[kk] * cldmax[idx]
    nevapr[idx] = nevapr[idx] - pre[kk] * cldmax[idx]
    nevapr2[idx] = nevapr2[idx] - pre[kk] * cldmax[idx]

    prain[idx] = prain[idx] + (pra[kk] + prc[kk]) * lcldm[idx] + (-pracs[kk] - mnuccr[kk]) * cldmax[idx]
    prodsnow[idx] = prodsnow[idx] + (prai[kk] + prci[kk]) * icldm[idx] + (psacws[kk] + bergs[kk]) * lcldm[idx] + (pracs[kk] + mnuccr[kk]) * cldmax[idx]

    qcsinksum_rate1ord[kk] = qcsinksum_rate1ord[kk] + (pra[kk] + prc[kk] + psacws[kk]) * lcldm[idx]
    qcsum_rate1ord[kk] = qcsum_rate1ord[kk] + qc[idx]

    prao[idx] = prao[idx] + pra[kk] * lcldm[idx]
    prco[idx] = prco[idx] + prc[kk] * lcldm[idx]
    mnuccco[idx] = mnuccco[idx] + mnuccc[kk] * lcldm[idx]
    mnuccto[idx] = mnuccto[idx] + mnucct[kk] * lcldm[idx]
    mnuccdo[idx] = mnuccdo[idx] + mnuccd[kk] * lcldm[idx]
    msacwio[idx] = msacwio[idx] + msacwi[kk] * lcldm[idx]
    psacwso[idx] = psacwso[idx] + psacws[kk] * lcldm[idx]
    bergso[idx] = bergso[idx] + bergs[kk] * lcldm[idx]
    bergo[idx] = bergo[idx] + berg[idx]
    prcio[idx] = prcio[idx] + prci[kk] * icldm[idx]
    praio[idx] = praio[idx] + prai[kk] * icldm[idx]
    mnuccro[idx] = mnuccro[idx] + mnuccr[kk] * cldmax[idx]
    pracso[idx] = pracso[idx] + pracs[kk] * cldmax[idx]
    preo[idx] = preo[idx] + pre[kk] * cldmax[idx]
    prdso[idx] = prdso[idx] + prds[kk] * cldmax[idx]


@export
def micro_mg1_0_post_iter_avg_codon(
    i: int,
    pcols: int,
    pver: int,
    top_lev: int,
    iter_count: int,
    prect1_p: cobj,
    preci1_p: cobj,
    prect_p: cobj,
    preci_p: cobj,
    t1_p: cobj,
    q1_p: cobj,
    qc1_p: cobj,
    qi1_p: cobj,
    nc1_p: cobj,
    ni1_p: cobj,
    t_p: cobj,
    q_p: cobj,
    qc_p: cobj,
    qi_p: cobj,
    nc_p: cobj,
    ni_p: cobj,
    tlat1_p: cobj,
    qvlat1_p: cobj,
    qctend1_p: cobj,
    qitend1_p: cobj,
    nctend1_p: cobj,
    nitend1_p: cobj,
    tlat_p: cobj,
    qvlat_p: cobj,
    qctend_p: cobj,
    qitend_p: cobj,
    nctend_p: cobj,
    nitend_p: cobj,
    rainrt1_p: cobj,
    rainrt_p: cobj,
    rflx1_p: cobj,
    sflx1_p: cobj,
    rflx_p: cobj,
    sflx_p: cobj,
    qrout_p: cobj,
    qsout_p: cobj,
    nrout_p: cobj,
    nsout_p: cobj,
    nevapr_p: cobj,
    nevapr2_p: cobj,
    evapsnow_p: cobj,
    prain_p: cobj,
    prodsnow_p: cobj,
    cmeout_p: cobj,
    cmeiout_p: cobj,
    meltsdt_p: cobj,
    frzrdt_p: cobj,
    prao_p: cobj,
    prco_p: cobj,
    mnuccco_p: cobj,
    mnuccto_p: cobj,
    msacwio_p: cobj,
    psacwso_p: cobj,
    bergso_p: cobj,
    bergo_p: cobj,
    prcio_p: cobj,
    praio_p: cobj,
    mnuccro_p: cobj,
    pracso_p: cobj,
    mnuccdo_p: cobj,
    preo_p: cobj,
    prdso_p: cobj,
    frzro_p: cobj,
    meltso_p: cobj,
    wtprelat_p: cobj,
    prer_evap_p: cobj,
):
    prect1 = Ptr[float](prect1_p)
    preci1 = Ptr[float](preci1_p)
    prect = Ptr[float](prect_p)
    preci = Ptr[float](preci_p)
    t1 = Ptr[float](t1_p)
    q1 = Ptr[float](q1_p)
    qc1 = Ptr[float](qc1_p)
    qi1 = Ptr[float](qi1_p)
    nc1 = Ptr[float](nc1_p)
    ni1 = Ptr[float](ni1_p)
    t = Ptr[float](t_p)
    q = Ptr[float](q_p)
    qc = Ptr[float](qc_p)
    qi = Ptr[float](qi_p)
    nc = Ptr[float](nc_p)
    ni = Ptr[float](ni_p)
    tlat1 = Ptr[float](tlat1_p)
    qvlat1 = Ptr[float](qvlat1_p)
    qctend1 = Ptr[float](qctend1_p)
    qitend1 = Ptr[float](qitend1_p)
    nctend1 = Ptr[float](nctend1_p)
    nitend1 = Ptr[float](nitend1_p)
    tlat = Ptr[float](tlat_p)
    qvlat = Ptr[float](qvlat_p)
    qctend = Ptr[float](qctend_p)
    qitend = Ptr[float](qitend_p)
    nctend = Ptr[float](nctend_p)
    nitend = Ptr[float](nitend_p)
    rainrt1 = Ptr[float](rainrt1_p)
    rainrt = Ptr[float](rainrt_p)
    rflx1 = Ptr[float](rflx1_p)
    sflx1 = Ptr[float](sflx1_p)
    rflx = Ptr[float](rflx_p)
    sflx = Ptr[float](sflx_p)
    qrout = Ptr[float](qrout_p)
    qsout = Ptr[float](qsout_p)
    nrout = Ptr[float](nrout_p)
    nsout = Ptr[float](nsout_p)
    nevapr = Ptr[float](nevapr_p)
    nevapr2 = Ptr[float](nevapr2_p)
    evapsnow = Ptr[float](evapsnow_p)
    prain = Ptr[float](prain_p)
    prodsnow = Ptr[float](prodsnow_p)
    cmeout = Ptr[float](cmeout_p)
    cmeiout = Ptr[float](cmeiout_p)
    meltsdt = Ptr[float](meltsdt_p)
    frzrdt = Ptr[float](frzrdt_p)
    prao = Ptr[float](prao_p)
    prco = Ptr[float](prco_p)
    mnuccco = Ptr[float](mnuccco_p)
    mnuccto = Ptr[float](mnuccto_p)
    msacwio = Ptr[float](msacwio_p)
    psacwso = Ptr[float](psacwso_p)
    bergso = Ptr[float](bergso_p)
    bergo = Ptr[float](bergo_p)
    prcio = Ptr[float](prcio_p)
    praio = Ptr[float](praio_p)
    mnuccro = Ptr[float](mnuccro_p)
    pracso = Ptr[float](pracso_p)
    mnuccdo = Ptr[float](mnuccdo_p)
    preo = Ptr[float](preo_p)
    prdso = Ptr[float](prdso_p)
    frzro = Ptr[float](frzro_p)
    meltso = Ptr[float](meltso_p)
    wtprelat = Ptr[float](wtprelat_p)
    prer_evap = Ptr[float](prer_evap_p)

    iter_real = float(iter_count)

    prect[i - 1] = prect1[i - 1] / iter_real
    preci[i - 1] = preci1[i - 1] / iter_real

    for k in range(top_lev, pver + 1):
        idx = _idx2(i, k, pcols)

        t[idx] = t1[idx]
        q[idx] = q1[idx]
        qc[idx] = qc1[idx]
        qi[idx] = qi1[idx]
        nc[idx] = nc1[idx]
        ni[idx] = ni1[idx]

        tlat[idx] = tlat1[idx] / iter_real
        qvlat[idx] = qvlat1[idx] / iter_real
        qctend[idx] = qctend1[idx] / iter_real
        qitend[idx] = qitend1[idx] / iter_real
        nctend[idx] = nctend1[idx] / iter_real
        nitend[idx] = nitend1[idx] / iter_real

        rainrt[idx] = rainrt1[idx] / iter_real

        idx_flux = _idx2(i, k + 1, pcols)
        rflx[idx_flux] = rflx1[idx_flux] / iter_real
        sflx[idx_flux] = sflx1[idx_flux] / iter_real

        qrout[idx] = qrout[idx] / iter_real
        qsout[idx] = qsout[idx] / iter_real
        nrout[idx] = nrout[idx] / iter_real
        nsout[idx] = nsout[idx] / iter_real

        nevapr[idx] = nevapr[idx] / iter_real
        nevapr2[idx] = nevapr2[idx] / iter_real
        evapsnow[idx] = evapsnow[idx] / iter_real
        prain[idx] = prain[idx] / iter_real
        prodsnow[idx] = prodsnow[idx] / iter_real
        cmeout[idx] = cmeout[idx] / iter_real

        cmeiout[idx] = cmeiout[idx] / iter_real
        meltsdt[idx] = meltsdt[idx] / iter_real
        frzrdt[idx] = frzrdt[idx] / iter_real

        prao[idx] = prao[idx] / iter_real
        prco[idx] = prco[idx] / iter_real
        mnuccco[idx] = mnuccco[idx] / iter_real
        mnuccto[idx] = mnuccto[idx] / iter_real
        msacwio[idx] = msacwio[idx] / iter_real
        psacwso[idx] = psacwso[idx] / iter_real
        bergso[idx] = bergso[idx] / iter_real
        bergo[idx] = bergo[idx] / iter_real
        prcio[idx] = prcio[idx] / iter_real
        praio[idx] = praio[idx] / iter_real

        mnuccro[idx] = mnuccro[idx] / iter_real
        pracso[idx] = pracso[idx] / iter_real

        mnuccdo[idx] = mnuccdo[idx] / iter_real

        preo[idx] = preo[idx] / iter_real
        prdso[idx] = prdso[idx] / iter_real
        frzro[idx] = frzro[idx] / iter_real
        meltso[idx] = meltso[idx] / iter_real
        wtprelat[idx] = tlat[idx]

        nevapr[idx] = nevapr[idx] + evapsnow[idx]
        prer_evap[idx] = nevapr2[idx]
        prain[idx] = prain[idx] + prodsnow[idx]


@export
def micro_mg1_0_phase_change_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    deltat: float,
    cpp: float,
    xlf: float,
    tmelt: float,
    qsmall: float,
    pi_v: float,
    rhow: float,
    do_cldice: int,
    qc_p: cobj,
    qi_p: cobj,
    nc_p: cobj,
    ni_p: cobj,
    t_p: cobj,
    qctend_p: cobj,
    qitend_p: cobj,
    nctend_p: cobj,
    nitend_p: cobj,
    tlat_p: cobj,
    dumc_p: cobj,
    dumi_p: cobj,
    dumnc_p: cobj,
    dumni_p: cobj,
    melto_p: cobj,
    homoo_p: cobj,
    wtpostlat_p: cobj,
):
    qc = Ptr[float](qc_p)
    qi = Ptr[float](qi_p)
    nc = Ptr[float](nc_p)
    ni = Ptr[float](ni_p)
    t = Ptr[float](t_p)
    qctend = Ptr[float](qctend_p)
    qitend = Ptr[float](qitend_p)
    nctend = Ptr[float](nctend_p)
    nitend = Ptr[float](nitend_p)
    tlat = Ptr[float](tlat_p)
    dumc = Ptr[float](dumc_p)
    dumi = Ptr[float](dumi_p)
    dumnc = Ptr[float](dumnc_p)
    dumni = Ptr[float](dumni_p)
    melto = Ptr[float](melto_p)
    homoo = Ptr[float](homoo_p)
    wtpostlat = Ptr[float](wtpostlat_p)

    idx = _idx2(i, k, pcols)

    dumc[idx] = max(qc[idx] + qctend[idx] * deltat, 0.0)
    dumi[idx] = max(qi[idx] + qitend[idx] * deltat, 0.0)
    dumnc[idx] = max(nc[idx] + nctend[idx] * deltat, 0.0)
    dumni[idx] = max(ni[idx] + nitend[idx] * deltat, 0.0)

    if dumc[idx] < qsmall:
        dumnc[idx] = 0.0
    if dumi[idx] < qsmall:
        dumni[idx] = 0.0

    if do_cldice != 0:
        if t[idx] + tlat[idx] / cpp * deltat > tmelt:
            if dumi[idx] > 0.0:
                dum = -dumi[idx] * xlf / cpp
                if t[idx] + tlat[idx] / cpp * deltat + dum < tmelt:
                    dum = (t[idx] + tlat[idx] / cpp * deltat - tmelt) * cpp / xlf
                    dum = dum / dumi[idx] * xlf / cpp
                    dum = max(0.0, dum)
                    dum = min(1.0, dum)
                else:
                    dum = 1.0

                qctend[idx] = qctend[idx] + dum * dumi[idx] / deltat
                melto[idx] = dum * dumi[idx] / deltat
                nctend[idx] = nctend[idx] + 3.0 * dum * dumi[idx] / deltat / (4.0 * pi_v * 5.12e-16 * rhow)
                qitend[idx] = ((1.0 - dum) * dumi[idx] - qi[idx]) / deltat
                nitend[idx] = ((1.0 - dum) * dumni[idx] - ni[idx]) / deltat
                tlat[idx] = tlat[idx] - xlf * dum * dumi[idx] / deltat
                wtpostlat[idx] = wtpostlat[idx] - (xlf * dum * dumi[idx] / deltat)

        if t[idx] + tlat[idx] / cpp * deltat < 233.15:
            if dumc[idx] > 0.0:
                dum = dumc[idx] * xlf / cpp
                if t[idx] + tlat[idx] / cpp * deltat + dum > 233.15:
                    dum = -(t[idx] + tlat[idx] / cpp * deltat - 233.15) * cpp / xlf
                    dum = dum / dumc[idx] * xlf / cpp
                    dum = max(0.0, dum)
                    dum = min(1.0, dum)
                else:
                    dum = 1.0

                qitend[idx] = qitend[idx] + dum * dumc[idx] / deltat
                homoo[idx] = dum * dumc[idx] / deltat
                nitend[idx] = nitend[idx] + dum * 3.0 * dumc[idx] / (4.0 * 3.14 * 1.563e-14 * 500.0) / deltat
                qctend[idx] = ((1.0 - dum) * dumc[idx] - qc[idx]) / deltat
                nctend[idx] = ((1.0 - dum) * dumnc[idx] - nc[idx]) / deltat
                tlat[idx] = tlat[idx] + xlf * dum * dumc[idx] / deltat
                wtpostlat[idx] = wtpostlat[idx] + (xlf * dum * dumc[idx] / deltat)


@export
def micro_mg1_0_number_cleanup_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    deltat: float,
    qsmall: float,
    do_cldice: int,
    qc_p: cobj,
    qi_p: cobj,
    nc_p: cobj,
    ni_p: cobj,
    qctend_p: cobj,
    qitend_p: cobj,
    nctend_p: cobj,
    nitend_p: cobj,
):
    qc = Ptr[float](qc_p)
    qi = Ptr[float](qi_p)
    nc = Ptr[float](nc_p)
    ni = Ptr[float](ni_p)
    qctend = Ptr[float](qctend_p)
    qitend = Ptr[float](qitend_p)
    nctend = Ptr[float](nctend_p)
    nitend = Ptr[float](nitend_p)

    for i in range(1, ncol + 1):
        for k in range(top_lev, pver + 1):
            idx = _idx2(i, k, pcols)
            if qc[idx] + qctend[idx] * deltat < qsmall:
                nctend[idx] = -nc[idx] / deltat
            if do_cldice != 0 and qi[idx] + qitend[idx] * deltat < qsmall:
                nitend[idx] = -ni[idx] / deltat


@export
def micro_mg1_0_reflectivity_flags_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    mindbz: float,
    csmin: float,
    csmax: float,
    refl_p: cobj,
    arefl_p: cobj,
    areflz_p: cobj,
    frefl_p: cobj,
    csrfl_p: cobj,
    acsrfl_p: cobj,
    fcsrfl_p: cobj,
):
    refl = Ptr[float](refl_p)
    arefl = Ptr[float](arefl_p)
    areflz = Ptr[float](areflz_p)
    frefl = Ptr[float](frefl_p)
    csrfl = Ptr[float](csrfl_p)
    acsrfl = Ptr[float](acsrfl_p)
    fcsrfl = Ptr[float](fcsrfl_p)

    for i in range(1, ncol + 1):
        for k in range(top_lev, pver + 1):
            idx = _idx2(i, k, pcols)
            if refl[idx] > mindbz:
                arefl[idx] = refl[idx]
                frefl[idx] = 1.0
            else:
                arefl[idx] = 0.0
                areflz[idx] = 0.0
                frefl[idx] = 0.0

            csrfl[idx] = min(csmax, refl[idx])

            if csrfl[idx] > csmin:
                acsrfl[idx] = refl[idx]
                fcsrfl[idx] = 1.0
            else:
                acsrfl[idx] = 0.0
                fcsrfl[idx] = 0.0


@export
def micro_mg1_0_substep_accum_column_codon(
    i: int,
    pcols: int,
    pver: int,
    top_lev: int,
    deltat: float,
    cpp: float,
    qric_p: cobj,
    qniic_p: cobj,
    nric_p: cobj,
    nsic_p: cobj,
    rho_p: cobj,
    cldmax_p: cobj,
    qrout_p: cobj,
    qsout_p: cobj,
    nrout_p: cobj,
    nsout_p: cobj,
    tlat_p: cobj,
    qvlat_p: cobj,
    qctend_p: cobj,
    qitend_p: cobj,
    nctend_p: cobj,
    nitend_p: cobj,
    tlat1_p: cobj,
    qvlat1_p: cobj,
    qctend1_p: cobj,
    qitend1_p: cobj,
    nctend1_p: cobj,
    nitend1_p: cobj,
    t_p: cobj,
    q_p: cobj,
    qc_p: cobj,
    qi_p: cobj,
    nc_p: cobj,
    ni_p: cobj,
    rainrt_p: cobj,
    rainrt1_p: cobj,
    arcld_p: cobj,
    rercld_p: cobj,
    rflx_p: cobj,
    sflx_p: cobj,
    rflx1_p: cobj,
    sflx1_p: cobj,
    umr_p: cobj,
    ums_p: cobj,
):
    qric = Ptr[float](qric_p)
    qniic = Ptr[float](qniic_p)
    nric = Ptr[float](nric_p)
    nsic = Ptr[float](nsic_p)
    rho = Ptr[float](rho_p)
    cldmax = Ptr[float](cldmax_p)
    qrout = Ptr[float](qrout_p)
    qsout = Ptr[float](qsout_p)
    nrout = Ptr[float](nrout_p)
    nsout = Ptr[float](nsout_p)
    tlat = Ptr[float](tlat_p)
    qvlat = Ptr[float](qvlat_p)
    qctend = Ptr[float](qctend_p)
    qitend = Ptr[float](qitend_p)
    nctend = Ptr[float](nctend_p)
    nitend = Ptr[float](nitend_p)
    tlat1 = Ptr[float](tlat1_p)
    qvlat1 = Ptr[float](qvlat1_p)
    qctend1 = Ptr[float](qctend1_p)
    qitend1 = Ptr[float](qitend1_p)
    nctend1 = Ptr[float](nctend1_p)
    nitend1 = Ptr[float](nitend1_p)
    t = Ptr[float](t_p)
    q = Ptr[float](q_p)
    qc = Ptr[float](qc_p)
    qi = Ptr[float](qi_p)
    nc = Ptr[float](nc_p)
    ni = Ptr[float](ni_p)
    rainrt = Ptr[float](rainrt_p)
    rainrt1 = Ptr[float](rainrt1_p)
    arcld = Ptr[float](arcld_p)
    rercld = Ptr[float](rercld_p)
    rflx = Ptr[float](rflx_p)
    sflx = Ptr[float](sflx_p)
    rflx1 = Ptr[float](rflx1_p)
    sflx1 = Ptr[float](sflx1_p)
    umr = Ptr[float](umr_p)
    ums = Ptr[float](ums_p)

    rflx[_idx2(i, 1, pcols)] = 0.0
    sflx[_idx2(i, 1, pcols)] = 0.0

    for k in range(top_lev, pver + 1):
        idx = _idx2(i, k, pcols)
        qrout[idx] = qrout[idx] + qric[idx] * cldmax[idx]
        qsout[idx] = qsout[idx] + qniic[idx] * cldmax[idx]
        nrout[idx] = nrout[idx] + nric[idx] * rho[idx] * cldmax[idx]
        nsout[idx] = nsout[idx] + nsic[idx] * rho[idx] * cldmax[idx]

        tlat1[idx] = tlat1[idx] + tlat[idx]
        qvlat1[idx] = qvlat1[idx] + qvlat[idx]
        qctend1[idx] = qctend1[idx] + qctend[idx]
        qitend1[idx] = qitend1[idx] + qitend[idx]
        nctend1[idx] = nctend1[idx] + nctend[idx]
        nitend1[idx] = nitend1[idx] + nitend[idx]

        t[idx] = t[idx] + tlat[idx] * deltat / cpp
        q[idx] = q[idx] + qvlat[idx] * deltat
        qc[idx] = qc[idx] + qctend[idx] * deltat
        qi[idx] = qi[idx] + qitend[idx] * deltat
        nc[idx] = nc[idx] + nctend[idx] * deltat
        ni[idx] = ni[idx] + nitend[idx] * deltat

        rainrt1[idx] = rainrt1[idx] + rainrt[idx]

        if arcld[idx] > 0.0:
            rercld[idx] = rercld[idx] / arcld[idx]

        idx_flux = _idx2(i, k + 1, pcols)
        rflx[idx_flux] = qrout[idx] * rho[idx] * umr[k - 1]
        sflx[idx_flux] = qsout[idx] * rho[idx] * ums[k - 1]
        rflx1[idx_flux] = rflx1[idx_flux] + rflx[idx_flux]
        sflx1[idx_flux] = sflx1[idx_flux] + sflx[idx_flux]


@export
def micro_mg1_0_sedimentation_state_codon(
    stage: int,
    i: int,
    k: int,
    pcols: int,
    pver: int,
    deltat: float,
    qsmall: float,
    qc_p: cobj,
    qi_p: cobj,
    nc_p: cobj,
    ni_p: cobj,
    qctend_p: cobj,
    qitend_p: cobj,
    nctend_p: cobj,
    nitend_p: cobj,
    lcldm_p: cobj,
    icldm_p: cobj,
    dumc_p: cobj,
    dumi_p: cobj,
    dumnc_p: cobj,
    dumni_p: cobj,
):
    qc = Ptr[float](qc_p)
    qi = Ptr[float](qi_p)
    nc = Ptr[float](nc_p)
    ni = Ptr[float](ni_p)
    qctend = Ptr[float](qctend_p)
    qitend = Ptr[float](qitend_p)
    nctend = Ptr[float](nctend_p)
    nitend = Ptr[float](nitend_p)
    lcldm = Ptr[float](lcldm_p)
    icldm = Ptr[float](icldm_p)
    dumc = Ptr[float](dumc_p)
    dumi = Ptr[float](dumi_p)
    dumnc = Ptr[float](dumnc_p)
    dumni = Ptr[float](dumni_p)

    idx = _idx2(i, k, pcols)
    if stage == 0:
        dumc[idx] = (qc[idx] + qctend[idx] * deltat) / lcldm[idx]
        dumi[idx] = (qi[idx] + qitend[idx] * deltat) / icldm[idx]
        dumnc[idx] = max((nc[idx] + nctend[idx] * deltat) / lcldm[idx], 0.0)
        dumni[idx] = max((ni[idx] + nitend[idx] * deltat) / icldm[idx], 0.0)
    else:
        dumc[idx] = qc[idx] + qctend[idx] * deltat
        dumi[idx] = qi[idx] + qitend[idx] * deltat
        dumnc[idx] = max(nc[idx] + nctend[idx] * deltat, 0.0)
        dumni[idx] = max(ni[idx] + nitend[idx] * deltat, 0.0)

        if dumc[idx] < qsmall:
            dumnc[idx] = 0.0
        if dumi[idx] < qsmall:
            dumni[idx] = 0.0


@export
def micro_mg1_0_sedimentation_velocity_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    qsmall: float,
    g: float,
    umc: float,
    unc: float,
    umi: float,
    uni: float,
    rho_p: cobj,
    dumc_p: cobj,
    dumi_p: cobj,
    vtrmc_p: cobj,
    vtrmi_p: cobj,
    fi_p: cobj,
    fni_p: cobj,
    fc_p: cobj,
    fnc_p: cobj,
    wtfc_p: cobj,
    wtfi_p: cobj,
):
    rho = Ptr[float](rho_p)
    dumc = Ptr[float](dumc_p)
    dumi = Ptr[float](dumi_p)
    vtrmc = Ptr[float](vtrmc_p)
    vtrmi = Ptr[float](vtrmi_p)
    fi = Ptr[float](fi_p)
    fni = Ptr[float](fni_p)
    fc = Ptr[float](fc_p)
    fnc = Ptr[float](fnc_p)
    wtfc = Ptr[float](wtfc_p)
    wtfi = Ptr[float](wtfi_p)

    idx = _idx2(i, k, pcols)
    kidx = k - 1
    if dumc[idx] >= qsmall:
        vtrmc[idx] = umc
    if dumi[idx] >= qsmall:
        vtrmi[idx] = umi

    fi[kidx] = g * rho[idx] * umi
    fni[kidx] = g * rho[idx] * uni
    fc[kidx] = g * rho[idx] * umc
    fnc[kidx] = g * rho[idx] * unc

    wtfc[idx] = fc[kidx]
    wtfi[idx] = fi[kidx]


@export
def micro_mg1_0_sedimentation_ice_prep_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    dumi_p: cobj,
    dumni_p: cobj,
):
    dumi = Ptr[float](dumi_p)
    dumni = Ptr[float](dumni_p)

    idx = _idx2(i, k, pcols)
    dumni[idx] = min(dumni[idx], dumi[idx] * 1.0e20)


@export
def micro_mg1_0_sedimentation_liq_prep_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    cdnl: float,
    dumc_p: cobj,
    dumnc_p: cobj,
    rho_p: cobj,
):
    dumc = Ptr[float](dumc_p)
    dumnc = Ptr[float](dumnc_p)
    rho = Ptr[float](rho_p)

    idx = _idx2(i, k, pcols)
    dumnc[idx] = min(dumnc[idx], dumc[idx] * 1.0e20)
    dumnc[idx] = max(dumnc[idx], cdnl / rho[idx])


@export
def micro_mg1_0_effrad_state_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    deltat: float,
    max_incloud: float,
    qc_p: cobj,
    qi_p: cobj,
    nc_p: cobj,
    ni_p: cobj,
    qctend_p: cobj,
    qitend_p: cobj,
    nctend_p: cobj,
    nitend_p: cobj,
    lcldm_p: cobj,
    icldm_p: cobj,
    dumc_p: cobj,
    dumi_p: cobj,
    dumnc_p: cobj,
    dumni_p: cobj,
):
    qc = Ptr[float](qc_p)
    qi = Ptr[float](qi_p)
    nc = Ptr[float](nc_p)
    ni = Ptr[float](ni_p)
    qctend = Ptr[float](qctend_p)
    qitend = Ptr[float](qitend_p)
    nctend = Ptr[float](nctend_p)
    nitend = Ptr[float](nitend_p)
    lcldm = Ptr[float](lcldm_p)
    icldm = Ptr[float](icldm_p)
    dumc = Ptr[float](dumc_p)
    dumi = Ptr[float](dumi_p)
    dumnc = Ptr[float](dumnc_p)
    dumni = Ptr[float](dumni_p)

    idx = _idx2(i, k, pcols)
    dumc[idx] = max(qc[idx] + qctend[idx] * deltat, 0.0) / lcldm[idx]
    dumi[idx] = max(qi[idx] + qitend[idx] * deltat, 0.0) / icldm[idx]
    dumnc[idx] = max(nc[idx] + nctend[idx] * deltat, 0.0) / lcldm[idx]
    dumni[idx] = max(ni[idx] + nitend[idx] * deltat, 0.0) / icldm[idx]

    dumc[idx] = min(dumc[idx], max_incloud)
    dumi[idx] = min(dumi[idx], max_incloud)


@export
def micro_mg1_0_effdiam_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    rhoi: float,
    do_cldice: int,
    effi_p: cobj,
    deffi_p: cobj,
):
    effi = Ptr[float](effi_p)
    deffi = Ptr[float](deffi_p)

    idx = _idx2(i, k, pcols)
    if do_cldice != 0:
        deffi[idx] = effi[idx] * rhoi / 917.0 * 2.0
    else:
        deffi[idx] = effi[idx] * 2.0


@export
def micro_mg1_0_effrad_liq_prep_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    deltat: float,
    cdnl: float,
    dumc_p: cobj,
    dumnc_p: cobj,
    rho_p: cobj,
    lcldm_p: cobj,
    nc_p: cobj,
    nctend_p: cobj,
):
    dumc = Ptr[float](dumc_p)
    dumnc = Ptr[float](dumnc_p)
    rho = Ptr[float](rho_p)
    lcldm = Ptr[float](lcldm_p)
    nc = Ptr[float](nc_p)
    nctend = Ptr[float](nctend_p)

    idx = _idx2(i, k, pcols)
    dumnc[idx] = min(dumnc[idx], dumc[idx] * 1.0e20)
    min_droplet = cdnl / rho[idx]
    if dumnc[idx] < min_droplet:
        nctend[idx] = (min_droplet * lcldm[idx] - nc[idx]) / deltat
    dumnc[idx] = max(dumnc[idx], min_droplet)


@export
def micro_mg1_0_effrad_ice_prep_codon(
    i: int,
    k: int,
    pcols: int,
    pver: int,
    dumi_p: cobj,
    dumni_p: cobj,
):
    dumi = Ptr[float](dumi_p)
    dumni = Ptr[float](dumni_p)

    idx = _idx2(i, k, pcols)
    dumni[idx] = min(dumni[idx], dumi[idx] * 1.0e20)


@export
def micro_mg1_0_sedimentation_fallout_codon(
    i: int,
    pcols: int,
    pver: int,
    top_lev: int,
    nstep: int,
    do_cldice: int,
    deltat: float,
    g: float,
    xxlv: float,
    xxls: float,
    pdel_p: cobj,
    lcldm_p: cobj,
    icldm_p: cobj,
    qitend_p: cobj,
    nitend_p: cobj,
    qctend_p: cobj,
    nctend_p: cobj,
    qcsedten_p: cobj,
    qisedten_p: cobj,
    qvlat_p: cobj,
    qisevap_p: cobj,
    qcsevap_p: cobj,
    tlat_p: cobj,
    dumi_p: cobj,
    dumni_p: cobj,
    dumc_p: cobj,
    dumnc_p: cobj,
    fi_p: cobj,
    fni_p: cobj,
    fc_p: cobj,
    fnc_p: cobj,
    falouti_p: cobj,
    faloutni_p: cobj,
    faloutc_p: cobj,
    faloutnc_p: cobj,
    prect_p: cobj,
    preci_p: cobj,
):
    pdel = Ptr[float](pdel_p)
    lcldm = Ptr[float](lcldm_p)
    icldm = Ptr[float](icldm_p)
    qitend = Ptr[float](qitend_p)
    nitend = Ptr[float](nitend_p)
    qctend = Ptr[float](qctend_p)
    nctend = Ptr[float](nctend_p)
    qcsedten = Ptr[float](qcsedten_p)
    qisedten = Ptr[float](qisedten_p)
    qvlat = Ptr[float](qvlat_p)
    qisevap = Ptr[float](qisevap_p)
    qcsevap = Ptr[float](qcsevap_p)
    tlat = Ptr[float](tlat_p)
    dumi = Ptr[float](dumi_p)
    dumni = Ptr[float](dumni_p)
    dumc = Ptr[float](dumc_p)
    dumnc = Ptr[float](dumnc_p)
    fi = Ptr[float](fi_p)
    fni = Ptr[float](fni_p)
    fc = Ptr[float](fc_p)
    fnc = Ptr[float](fnc_p)
    falouti = Ptr[float](falouti_p)
    faloutni = Ptr[float](faloutni_p)
    faloutc = Ptr[float](faloutc_p)
    faloutnc = Ptr[float](faloutnc_p)
    prect = Ptr[float](prect_p)
    preci = Ptr[float](preci_p)
    nstep_f = float(nstep)

    for _n in range(1, nstep + 1):
        for k in range(top_lev, pver + 1):
            idx = _idx2(i, k, pcols)
            kidx = k - 1
            if do_cldice != 0:
                falouti[kidx] = fi[kidx] * dumi[idx]
                faloutni[kidx] = fni[kidx] * dumni[idx]
            else:
                falouti[kidx] = 0.0
                faloutni[kidx] = 0.0

            faloutc[kidx] = fc[kidx] * dumc[idx]
            faloutnc[kidx] = fnc[kidx] * dumnc[idx]

        k = top_lev
        idx_top = _idx2(i, k, pcols)
        kidx_top = k - 1
        faltndi = falouti[kidx_top] / pdel[idx_top]
        faltndni = faloutni[kidx_top] / pdel[idx_top]
        faltndc = faloutc[kidx_top] / pdel[idx_top]
        faltndnc = faloutnc[kidx_top] / pdel[idx_top]

        qitend[idx_top] = qitend[idx_top] - faltndi / nstep_f
        nitend[idx_top] = nitend[idx_top] - faltndni / nstep_f
        qctend[idx_top] = qctend[idx_top] - faltndc / nstep_f
        nctend[idx_top] = nctend[idx_top] - faltndnc / nstep_f

        qcsedten[idx_top] = qcsedten[idx_top] - faltndc / nstep_f
        qisedten[idx_top] = qisedten[idx_top] - faltndi / nstep_f

        dumi[idx_top] = dumi[idx_top] - faltndi * deltat / nstep_f
        dumni[idx_top] = dumni[idx_top] - faltndni * deltat / nstep_f
        dumc[idx_top] = dumc[idx_top] - faltndc * deltat / nstep_f
        dumnc[idx_top] = dumnc[idx_top] - faltndnc * deltat / nstep_f

        for k in range(top_lev + 1, pver + 1):
            idx = _idx2(i, k, pcols)
            idx_prev = _idx2(i, k - 1, pcols)
            kidx = k - 1
            kidx_prev = k - 2

            dum = lcldm[idx] / lcldm[idx_prev]
            dum = min(dum, 1.0)
            dum1 = icldm[idx] / icldm[idx_prev]
            dum1 = min(dum1, 1.0)

            faltndqie = (falouti[kidx] - falouti[kidx_prev]) / pdel[idx]
            faltndi = (falouti[kidx] - dum1 * falouti[kidx_prev]) / pdel[idx]
            faltndni = (faloutni[kidx] - dum1 * faloutni[kidx_prev]) / pdel[idx]
            faltndqce = (faloutc[kidx] - faloutc[kidx_prev]) / pdel[idx]
            faltndc = (faloutc[kidx] - dum * faloutc[kidx_prev]) / pdel[idx]
            faltndnc = (faloutnc[kidx] - dum * faloutnc[kidx_prev]) / pdel[idx]

            qitend[idx] = qitend[idx] - faltndi / nstep_f
            nitend[idx] = nitend[idx] - faltndni / nstep_f
            qctend[idx] = qctend[idx] - faltndc / nstep_f
            nctend[idx] = nctend[idx] - faltndnc / nstep_f

            qcsedten[idx] = qcsedten[idx] - faltndc / nstep_f
            qisedten[idx] = qisedten[idx] - faltndi / nstep_f

            qvlat[idx] = qvlat[idx] - (faltndqie - faltndi) / nstep_f
            qisevap[idx] = qisevap[idx] - (faltndqie - faltndi) / nstep_f
            qvlat[idx] = qvlat[idx] - (faltndqce - faltndc) / nstep_f
            qcsevap[idx] = qcsevap[idx] - (faltndqce - faltndc) / nstep_f

            tlat[idx] = tlat[idx] + (faltndqie - faltndi) * xxls / nstep_f
            tlat[idx] = tlat[idx] + (faltndqce - faltndc) * xxlv / nstep_f

            dumi[idx] = dumi[idx] - faltndi * deltat / nstep_f
            dumni[idx] = dumni[idx] - faltndni * deltat / nstep_f
            dumc[idx] = dumc[idx] - faltndc * deltat / nstep_f
            dumnc[idx] = dumnc[idx] - faltndnc * deltat / nstep_f

            fni[kidx] = max(fni[kidx] / pdel[idx], fni[kidx_prev] / pdel[idx_prev]) * pdel[idx]
            fi[kidx] = max(fi[kidx] / pdel[idx], fi[kidx_prev] / pdel[idx_prev]) * pdel[idx]
            fnc[kidx] = max(fnc[kidx] / pdel[idx], fnc[kidx_prev] / pdel[idx_prev]) * pdel[idx]
            fc[kidx] = max(fc[kidx] / pdel[idx], fc[kidx_prev] / pdel[idx_prev]) * pdel[idx]

        prect[i - 1] = prect[i - 1] + (faloutc[pver - 1] + falouti[pver - 1]) / g / nstep_f / 1000.0
        preci[i - 1] = preci[i - 1] + falouti[pver - 1] / g / nstep_f / 1000.0


@export
def micro_mg1_0_flux_ltrue_init_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    qsmall: float,
    rflx1_p: cobj,
    sflx1_p: cobj,
    rflx_p: cobj,
    sflx_p: cobj,
    qc_p: cobj,
    qi_p: cobj,
    cmei_p: cobj,
    ltrue_p: cobj,
):
    rflx1 = Ptr[float](rflx1_p)
    sflx1 = Ptr[float](sflx1_p)
    rflx = Ptr[float](rflx_p)
    sflx = Ptr[float](sflx_p)
    qc = Ptr[float](qc_p)
    qi = Ptr[float](qi_p)
    cmei = Ptr[float](cmei_p)
    ltrue = Ptr[i32](ltrue_p)

    for i in range(1, ncol + 1):
        rflx1[_idx2(i, 1, pcols)] = 0.0
        sflx1[_idx2(i, 1, pcols)] = 0.0
        for k in range(top_lev, pver + 1):
            rflx1[_idx2(i, k + 1, pcols)] = 0.0
            sflx1[_idx2(i, k + 1, pcols)] = 0.0

    for i in range(1, ncol + 1):
        rflx[_idx2(i, 1, pcols)] = 0.0
        sflx[_idx2(i, 1, pcols)] = 0.0
        for k in range(top_lev, pver + 1):
            rflx[_idx2(i, k + 1, pcols)] = 0.0
            sflx[_idx2(i, k + 1, pcols)] = 0.0

    for i in range(1, ncol + 1):
        ltrue[i - 1] = i32(0)
        for k in range(top_lev, pver + 1):
            idx = _idx2(i, k, pcols)
            if qc[idx] >= qsmall or qi[idx] >= qsmall or cmei[idx] >= qsmall:
                ltrue[i - 1] = i32(1)


@export
def micro_mg1_0_rate1ord_zero_codon(
    ncol: int,
    pcols: int,
    pver: int,
    rate1ord_cw2pr_st_p: cobj,
):
    rate1ord_cw2pr_st = Ptr[float](rate1ord_cw2pr_st_p)

    for i in range(1, ncol + 1):
        for k in range(1, pver + 1):
            rate1ord_cw2pr_st[_idx2(i, k, pcols)] = 0.0


@export
def micro_mg1_0_rate1ord_column_codon(
    i: int,
    pcols: int,
    pver: int,
    top_lev: int,
    qcsinksum_rate1ord_p: cobj,
    qcsum_rate1ord_p: cobj,
    rate1ord_cw2pr_st_p: cobj,
):
    qcsinksum_rate1ord = Ptr[float](qcsinksum_rate1ord_p)
    qcsum_rate1ord = Ptr[float](qcsum_rate1ord_p)
    rate1ord_cw2pr_st = Ptr[float](rate1ord_cw2pr_st_p)

    for k in range(top_lev, pver + 1):
        rate1ord_cw2pr_st[_idx2(i, k, pcols)] = qcsinksum_rate1ord[k - 1] / max(qcsum_rate1ord[k - 1], 1.0e-30)


@export
def micro_mg1_0_tail_activation_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    dum2i_p: cobj,
    dum2l_p: cobj,
    rho_p: cobj,
    ncai_p: cobj,
    ncal_p: cobj,
):
    dum2i = Ptr[float](dum2i_p)
    dum2l = Ptr[float](dum2l_p)
    rho = Ptr[float](rho_p)
    ncai = Ptr[float](ncai_p)
    ncal = Ptr[float](ncal_p)

    for i in range(1, ncol + 1):
        for k in range(top_lev, pver + 1):
            idx = _idx2(i, k, pcols)
            ncai[idx] = dum2i[idx] * rho[idx]
            ncal[idx] = dum2l[idx] * rho[idx]


@export
def micro_mg1_0_tail_avg_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    qrout_p: cobj,
    qsout_p: cobj,
    nrout_p: cobj,
    nsout_p: cobj,
    qrout2_p: cobj,
    qsout2_p: cobj,
    nrout2_p: cobj,
    nsout2_p: cobj,
    drout2_p: cobj,
    dsout2_p: cobj,
    freqs_p: cobj,
    freqr_p: cobj,
):
    qrout = Ptr[float](qrout_p)
    qsout = Ptr[float](qsout_p)
    nrout = Ptr[float](nrout_p)
    nsout = Ptr[float](nsout_p)
    qrout2 = Ptr[float](qrout2_p)
    qsout2 = Ptr[float](qsout2_p)
    nrout2 = Ptr[float](nrout2_p)
    nsout2 = Ptr[float](nsout2_p)
    drout2 = Ptr[float](drout2_p)
    dsout2 = Ptr[float](dsout2_p)
    freqs = Ptr[float](freqs_p)
    freqr = Ptr[float](freqr_p)

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            idx = _idx2(i, k, pcols)
            qrout2[idx] = 0.0
            qsout2[idx] = 0.0
            nrout2[idx] = 0.0
            nsout2[idx] = 0.0
            drout2[idx] = 0.0
            dsout2[idx] = 0.0
            freqs[idx] = 0.0
            freqr[idx] = 0.0

    for i in range(1, ncol + 1):
        for k in range(top_lev, pver + 1):
            idx = _idx2(i, k, pcols)
            if qrout[idx] > 1.0e-7 and nrout[idx] > 0.0:
                qrout2[idx] = qrout[idx]
                nrout2[idx] = nrout[idx]
                freqr[idx] = 1.0
            if qsout[idx] > 1.0e-7 and nsout[idx] > 0.0:
                qsout2[idx] = qsout[idx]
                nsout2[idx] = nsout[idx]
                freqs[idx] = 1.0


@export
def micro_mg1_0_tail_fice_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    deltat: float,
    qsmall: float,
    qc_p: cobj,
    qi_p: cobj,
    qctend_p: cobj,
    qitend_p: cobj,
    qsout_p: cobj,
    qrout_p: cobj,
    dumc_p: cobj,
    dumi_p: cobj,
    nfice_p: cobj,
):
    qc = Ptr[float](qc_p)
    qi = Ptr[float](qi_p)
    qctend = Ptr[float](qctend_p)
    qitend = Ptr[float](qitend_p)
    qsout = Ptr[float](qsout_p)
    qrout = Ptr[float](qrout_p)
    dumc = Ptr[float](dumc_p)
    dumi = Ptr[float](dumi_p)
    nfice = Ptr[float](nfice_p)

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            nfice[_idx2(i, k, pcols)] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx = _idx2(i, k, pcols)
            dumc[idx] = qc[idx] + qctend[idx] * deltat
            dumi[idx] = qi[idx] + qitend[idx] * deltat
            dumfice = qsout[idx] + qrout[idx] + dumc[idx] + dumi[idx]

            if dumfice > qsmall and (qsout[idx] + dumi[idx]) > qsmall:
                nfice[idx] = (qsout[idx] + dumi[idx]) / dumfice

            if nfice[idx] > 1.0:
                nfice[idx] = 1.0


@inline
def _pack2d_mgcols(
    mgncol: int,
    nlev: int,
    psetcols: int,
    top_lev: int,
    mgcols: Ptr[int],
    src: Ptr[float],
    dst: Ptr[float],
):
    for j in range(1, mgncol + 1):
        i = mgcols[j - 1]
        for kk in range(1, nlev + 1):
            k = top_lev + kk - 1
            dst[_idx2(j, kk, mgncol)] = src[_idx2(i, k, psetcols)]


@inline
def _pack_interface_mgcols(
    mgncol: int,
    nlev: int,
    psetcols: int,
    top_lev: int,
    mgcols: Ptr[int],
    src: Ptr[float],
    dst: Ptr[float],
):
    for j in range(1, mgncol + 1):
        i = mgcols[j - 1]
        for kk in range(1, nlev + 2):
            k = top_lev + kk - 1
            dst[_idx2(j, kk, mgncol)] = src[_idx2(i, k, psetcols)]


@inline
def _pack3d_mgcols(
    mgncol: int,
    nlev: int,
    psetcols: int,
    pver: int,
    top_lev: int,
    n3: int,
    mgcols: Ptr[int],
    src: Ptr[float],
    dst: Ptr[float],
):
    for m in range(1, n3 + 1):
        for kk in range(1, nlev + 1):
            k = top_lev + kk - 1
            for j in range(1, mgncol + 1):
                i = mgcols[j - 1]
                dst[_idx3(j, kk, m, mgncol, nlev)] = src[_idx3(i, k, m, psetcols, pver)]


@export
def micro_mg_cam_pack_static_inputs_codon(
    mgncol: int,
    nlev: int,
    psetcols: int,
    pver: int,
    pverp: int,
    top_lev: int,
    rndst_n3: int,
    nacon_n3: int,
    mgcols_p: cobj,
    relvar_p: cobj,
    accre_enhan_p: cobj,
    pmid_p: cobj,
    pdel_p: cobj,
    pint_p: cobj,
    ast_p: cobj,
    alst_mic_p: cobj,
    aist_mic_p: cobj,
    naai_p: cobj,
    npccn_p: cobj,
    rndst_p: cobj,
    nacon_p: cobj,
    packed_relvar_p: cobj,
    packed_accre_enhan_p: cobj,
    packed_p_p: cobj,
    packed_pdel_p: cobj,
    packed_pint_p: cobj,
    packed_cldn_p: cobj,
    packed_liqcldf_p: cobj,
    packed_icecldf_p: cobj,
    packed_naai_p: cobj,
    packed_npccn_p: cobj,
    packed_rndst_p: cobj,
    packed_nacon_p: cobj,
):
    mgcols = Ptr[int](mgcols_p)
    _pack2d_mgcols(mgncol, nlev, psetcols, top_lev, mgcols, Ptr[float](relvar_p), Ptr[float](packed_relvar_p))
    _pack2d_mgcols(
        mgncol, nlev, psetcols, top_lev, mgcols, Ptr[float](accre_enhan_p), Ptr[float](packed_accre_enhan_p)
    )
    _pack2d_mgcols(mgncol, nlev, psetcols, top_lev, mgcols, Ptr[float](pmid_p), Ptr[float](packed_p_p))
    _pack2d_mgcols(mgncol, nlev, psetcols, top_lev, mgcols, Ptr[float](pdel_p), Ptr[float](packed_pdel_p))
    _pack_interface_mgcols(
        mgncol, nlev, psetcols, top_lev, mgcols, Ptr[float](pint_p), Ptr[float](packed_pint_p)
    )
    _pack2d_mgcols(mgncol, nlev, psetcols, top_lev, mgcols, Ptr[float](ast_p), Ptr[float](packed_cldn_p))
    _pack2d_mgcols(
        mgncol, nlev, psetcols, top_lev, mgcols, Ptr[float](alst_mic_p), Ptr[float](packed_liqcldf_p)
    )
    _pack2d_mgcols(
        mgncol, nlev, psetcols, top_lev, mgcols, Ptr[float](aist_mic_p), Ptr[float](packed_icecldf_p)
    )
    _pack2d_mgcols(mgncol, nlev, psetcols, top_lev, mgcols, Ptr[float](naai_p), Ptr[float](packed_naai_p))
    _pack2d_mgcols(mgncol, nlev, psetcols, top_lev, mgcols, Ptr[float](npccn_p), Ptr[float](packed_npccn_p))
    _pack3d_mgcols(
        mgncol, nlev, psetcols, pver, top_lev, rndst_n3, mgcols, Ptr[float](rndst_p), Ptr[float](packed_rndst_p)
    )
    _pack3d_mgcols(
        mgncol, nlev, psetcols, pver, top_lev, nacon_n3, mgcols, Ptr[float](nacon_p), Ptr[float](packed_nacon_p)
    )


@export
def micro_mg_cam_pack_state_inputs_codon(
    mgncol: int,
    nlev: int,
    psetcols: int,
    pver: int,
    pcnst: int,
    top_lev: int,
    ixcldliq: int,
    ixcldice: int,
    ixnumliq: int,
    ixnumice: int,
    mgcols_p: cobj,
    state_t_p: cobj,
    state_q_p: cobj,
    packed_t_p: cobj,
    packed_q_p: cobj,
    packed_qc_p: cobj,
    packed_nc_p: cobj,
    packed_qi_p: cobj,
    packed_ni_p: cobj,
):
    mgcols = Ptr[int](mgcols_p)
    state_t = Ptr[float](state_t_p)
    state_q = Ptr[float](state_q_p)
    packed_t = Ptr[float](packed_t_p)
    packed_q = Ptr[float](packed_q_p)
    packed_qc = Ptr[float](packed_qc_p)
    packed_nc = Ptr[float](packed_nc_p)
    packed_qi = Ptr[float](packed_qi_p)
    packed_ni = Ptr[float](packed_ni_p)

    for j in range(1, mgncol + 1):
        i = mgcols[j - 1]
        for kk in range(1, nlev + 1):
            k = top_lev + kk - 1
            dst_idx = _idx2(j, kk, mgncol)
            packed_t[dst_idx] = state_t[_idx2(i, k, psetcols)]
            packed_q[dst_idx] = state_q[_idx3(i, k, 1, psetcols, pver)]
            packed_qc[dst_idx] = state_q[_idx3(i, k, ixcldliq, psetcols, pver)]
            packed_nc[dst_idx] = state_q[_idx3(i, k, ixnumliq, psetcols, pver)]
            packed_qi[dst_idx] = state_q[_idx3(i, k, ixcldice, psetcols, pver)]
            packed_ni[dst_idx] = state_q[_idx3(i, k, ixnumice, psetcols, pver)]


@export
def micro_mg_cam_pack_three_2d_inputs_codon(
    mgncol: int,
    nlev: int,
    psetcols: int,
    top_lev: int,
    mgcols_p: cobj,
    src1_p: cobj,
    src2_p: cobj,
    src3_p: cobj,
    dst1_p: cobj,
    dst2_p: cobj,
    dst3_p: cobj,
):
    mgcols = Ptr[int](mgcols_p)
    _pack2d_mgcols(mgncol, nlev, psetcols, top_lev, mgcols, Ptr[float](src1_p), Ptr[float](dst1_p))
    _pack2d_mgcols(mgncol, nlev, psetcols, top_lev, mgcols, Ptr[float](src2_p), Ptr[float](dst2_p))
    _pack2d_mgcols(mgncol, nlev, psetcols, top_lev, mgcols, Ptr[float](src3_p), Ptr[float](dst3_p))


@export
def micro_mg_cam_pack_precip_state_inputs_codon(
    mgncol: int,
    nlev: int,
    psetcols: int,
    pver: int,
    top_lev: int,
    ixrain: int,
    ixnumrain: int,
    ixsnow: int,
    ixnumsnow: int,
    mgcols_p: cobj,
    state_q_p: cobj,
    packed_qr_p: cobj,
    packed_nr_p: cobj,
    packed_qs_p: cobj,
    packed_ns_p: cobj,
):
    mgcols = Ptr[int](mgcols_p)
    state_q = Ptr[float](state_q_p)
    packed_qr = Ptr[float](packed_qr_p)
    packed_nr = Ptr[float](packed_nr_p)
    packed_qs = Ptr[float](packed_qs_p)
    packed_ns = Ptr[float](packed_ns_p)

    for j in range(1, mgncol + 1):
        i = mgcols[j - 1]
        for kk in range(1, nlev + 1):
            k = top_lev + kk - 1
            dst_idx = _idx2(j, kk, mgncol)
            packed_qr[dst_idx] = state_q[_idx3(i, k, ixrain, psetcols, pver)]
            packed_nr[dst_idx] = state_q[_idx3(i, k, ixnumrain, psetcols, pver)]
            packed_qs[dst_idx] = state_q[_idx3(i, k, ixsnow, psetcols, pver)]
            packed_ns[dst_idx] = state_q[_idx3(i, k, ixnumsnow, psetcols, pver)]


@inline
def _unpack2d_scalar_to_2d(
    mgncol: int,
    nlev: int,
    psetcols: int,
    pver: int,
    top_lev: int,
    mgcols: Ptr[int],
    packed: Ptr[float],
    dst: Ptr[float],
    fill: float,
):
    for k in range(1, pver + 1):
        for i in range(1, psetcols + 1):
            dst[_idx2(i, k, psetcols)] = fill

    for j in range(1, mgncol + 1):
        i = mgcols[j - 1]
        for kk in range(1, nlev + 1):
            k = top_lev + kk - 1
            dst[_idx2(i, k, psetcols)] = packed[_idx2(j, kk, mgncol)]


@inline
def _unpack2d_scalar_to_qslice(
    mgncol: int,
    nlev: int,
    psetcols: int,
    pver: int,
    top_lev: int,
    qidx: int,
    mgcols: Ptr[int],
    packed: Ptr[float],
    dst_q: Ptr[float],
    fill: float,
):
    for k in range(1, pver + 1):
        for i in range(1, psetcols + 1):
            dst_q[_idx3(i, k, qidx, psetcols, pver)] = fill

    for j in range(1, mgncol + 1):
        i = mgcols[j - 1]
        for kk in range(1, nlev + 1):
            k = top_lev + kk - 1
            dst_q[_idx3(i, k, qidx, psetcols, pver)] = packed[_idx2(j, kk, mgncol)]


@inline
def _unpack2d_statefill_to_qslice(
    mgncol: int,
    nlev: int,
    psetcols: int,
    pver: int,
    top_lev: int,
    qidx: int,
    mgcols: Ptr[int],
    state_q: Ptr[float],
    packed: Ptr[float],
    dst_q: Ptr[float],
    dtstep: float,
):
    for k in range(1, pver + 1):
        for i in range(1, psetcols + 1):
            dst_q[_idx3(i, k, qidx, psetcols, pver)] = -state_q[_idx3(i, k, qidx, psetcols, pver)] / dtstep

    for j in range(1, mgncol + 1):
        i = mgcols[j - 1]
        for kk in range(1, nlev + 1):
            k = top_lev + kk - 1
            dst_q[_idx3(i, k, qidx, psetcols, pver)] = packed[_idx2(j, kk, mgncol)]


@export
def micro_mg_cam_ptend_unpack_codon(
    mgncol: int,
    nlev: int,
    psetcols: int,
    pver: int,
    pcnst: int,
    top_lev: int,
    ixcldliq: int,
    ixcldice: int,
    ixnumliq: int,
    ixnumice: int,
    do_cldice: int,
    dtstep: float,
    mgcols_p: cobj,
    state_q_p: cobj,
    packed_tlat_p: cobj,
    packed_qvlat_p: cobj,
    packed_qctend_p: cobj,
    packed_qitend_p: cobj,
    packed_nctend_p: cobj,
    packed_nitend_p: cobj,
    ptend_s_p: cobj,
    ptend_q_p: cobj,
):
    mgcols = Ptr[int](mgcols_p)
    state_q = Ptr[float](state_q_p)
    ptend_q = Ptr[float](ptend_q_p)

    _unpack2d_scalar_to_2d(
        mgncol, nlev, psetcols, pver, top_lev, mgcols, Ptr[float](packed_tlat_p), Ptr[float](ptend_s_p), 0.0
    )
    _unpack2d_scalar_to_qslice(
        mgncol, nlev, psetcols, pver, top_lev, 1, mgcols, Ptr[float](packed_qvlat_p), ptend_q, 0.0
    )
    _unpack2d_scalar_to_qslice(
        mgncol, nlev, psetcols, pver, top_lev, ixcldliq, mgcols, Ptr[float](packed_qctend_p), ptend_q, 0.0
    )
    _unpack2d_scalar_to_qslice(
        mgncol, nlev, psetcols, pver, top_lev, ixcldice, mgcols, Ptr[float](packed_qitend_p), ptend_q, 0.0
    )
    _unpack2d_statefill_to_qslice(
        mgncol, nlev, psetcols, pver, top_lev, ixnumliq, mgcols, state_q, Ptr[float](packed_nctend_p), ptend_q, dtstep
    )
    if do_cldice != 0:
        _unpack2d_statefill_to_qslice(
            mgncol, nlev, psetcols, pver, top_lev, ixnumice, mgcols, state_q, Ptr[float](packed_nitend_p), ptend_q, dtstep
        )
    else:
        for k in range(1, pver + 1):
            for i in range(1, psetcols + 1):
                ptend_q[_idx3(i, k, ixnumice, psetcols, pver)] = 0.0


@inline
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
def _add_process_rate(
    process_rates,
    i: int,
    k: int,
    isrctype: int,
    idsttype: int,
    rtype: int,
    rate_val: float,
    pcols: int,
    pver: int,
    pwtype: int,
):
    dst_idx = _process_rates_idx(i, k, idsttype, isrctype, rtype, pcols, pver, pwtype)
    process_rates[dst_idx] = process_rates[dst_idx] + rate_val
    if isrctype != idsttype:
        src_idx = _process_rates_idx(i, k, isrctype, idsttype, rtype, pcols, pver, pwtype)
        process_rates[src_idx] = process_rates[src_idx] - rate_val


@export
def micro_mg_cam_premg_diag_codon(
    ncol: int,
    psetcols: int,
    pver: int,
    top_lev: int,
    ixcldliq: int,
    ixcldice: int,
    mincld: float,
    gravit: float,
    state_q_p: cobj,
    state_pdel_p: cobj,
    ast_p: cobj,
    cldo_p: cobj,
    iclwpi_p: cobj,
    iciwpi_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    state_pdel = Ptr[float](state_pdel_p)
    ast = Ptr[float](ast_p)
    cldo = Ptr[float](cldo_p)
    iclwpi = Ptr[float](iclwpi_p)
    iciwpi = Ptr[float](iciwpi_p)

    for i in range(1, psetcols + 1):
        iclwpi[i - 1] = 0.0
        iciwpi[i - 1] = 0.0

    for i in range(1, ncol + 1):
        for k in range(top_lev, pver + 1):
            idx2 = _idx2(i, k, psetcols)
            idx_q_liq = _idx3(i, k, ixcldliq, psetcols, pver)
            idx_q_ice = _idx3(i, k, ixcldice, psetcols, pver)
            iclwpi[i - 1] = iclwpi[i - 1] + min(
                state_q[idx_q_liq] / max(mincld, ast[idx2]), 0.005
            ) * state_pdel[idx2] / gravit
            iciwpi[i - 1] = iciwpi[i - 1] + min(
                state_q[idx_q_ice] / max(mincld, ast[idx2]), 0.005
            ) * state_pdel[idx2] / gravit

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, psetcols)
            cldo[idx2] = ast[idx2]


@export
def micro_mg_cam_wtrc_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    pwtype: int,
    iwtvap: int,
    iwtliq: int,
    iwtice: int,
    iwtstrain: int,
    iwtstsnow: int,
    preo_grid_p: cobj,
    prdso_grid_p: cobj,
    cmeiout_grid_p: cobj,
    meltso_grid_p: cobj,
    qcsedten_grid_p: cobj,
    qisedten_grid_p: cobj,
    mnuccco_grid_p: cobj,
    mnuccto_grid_p: cobj,
    msacwio_grid_p: cobj,
    prao_grid_p: cobj,
    prco_grid_p: cobj,
    psacwso_grid_p: cobj,
    bergo_grid_p: cobj,
    bergso_grid_p: cobj,
    praio_grid_p: cobj,
    prcio_grid_p: cobj,
    pracso_grid_p: cobj,
    mnuccro_grid_p: cobj,
    qcreso_grid_p: cobj,
    qireso_grid_p: cobj,
    homoo_grid_p: cobj,
    melto_grid_p: cobj,
    pre_rates_grid_p: cobj,
    sed_rates_grid_p: cobj,
    post_rates_grid_p: cobj,
    pcmei_grid_p: cobj,
    ncmei_grid_p: cobj,
    pmelts_grid_p: cobj,
    nmelts_grid_p: cobj,
):
    preo_grid = Ptr[float](preo_grid_p)
    prdso_grid = Ptr[float](prdso_grid_p)
    cmeiout_grid = Ptr[float](cmeiout_grid_p)
    meltso_grid = Ptr[float](meltso_grid_p)
    qcsedten_grid = Ptr[float](qcsedten_grid_p)
    qisedten_grid = Ptr[float](qisedten_grid_p)
    mnuccco_grid = Ptr[float](mnuccco_grid_p)
    mnuccto_grid = Ptr[float](mnuccto_grid_p)
    msacwio_grid = Ptr[float](msacwio_grid_p)
    prao_grid = Ptr[float](prao_grid_p)
    prco_grid = Ptr[float](prco_grid_p)
    psacwso_grid = Ptr[float](psacwso_grid_p)
    bergo_grid = Ptr[float](bergo_grid_p)
    bergso_grid = Ptr[float](bergso_grid_p)
    praio_grid = Ptr[float](praio_grid_p)
    prcio_grid = Ptr[float](prcio_grid_p)
    pracso_grid = Ptr[float](pracso_grid_p)
    mnuccro_grid = Ptr[float](mnuccro_grid_p)
    qcreso_grid = Ptr[float](qcreso_grid_p)
    qireso_grid = Ptr[float](qireso_grid_p)
    homoo_grid = Ptr[float](homoo_grid_p)
    melto_grid = Ptr[float](melto_grid_p)
    pre_rates_grid = Ptr[float](pre_rates_grid_p)
    sed_rates_grid = Ptr[float](sed_rates_grid_p)
    post_rates_grid = Ptr[float](post_rates_grid_p)
    pcmei_grid = Ptr[float](pcmei_grid_p)
    ncmei_grid = Ptr[float](ncmei_grid_p)
    pmelts_grid = Ptr[float](pmelts_grid_p)
    nmelts_grid = Ptr[float](nmelts_grid_p)

    for rtype in range(1, pwtype + 1):
        for isrctype in range(1, pwtype + 1):
            for idsttype in range(1, pwtype + 1):
                for k in range(top_lev, pver + 1):
                    for i in range(1, pcols + 1):
                        pre_rates_grid[
                            _process_rates_idx(
                                i, k, idsttype, isrctype, rtype, pcols, pver, pwtype
                            )
                        ] = 0.0
                        post_rates_grid[
                            _process_rates_idx(
                                i, k, idsttype, isrctype, rtype, pcols, pver, pwtype
                            )
                        ] = 0.0

    for m in range(1, pwtype + 1):
        for k in range(top_lev, pver + 1):
            for i in range(1, pcols + 1):
                sed_rates_grid[_idx3(i, k, m, pcols, pver)] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, pcols + 1):
            idx2 = _idx2(i, k, pcols)
            pcmei_grid[idx2] = 0.0
            ncmei_grid[idx2] = 0.0
            pmelts_grid[idx2] = 0.0
            nmelts_grid[idx2] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            if cmeiout_grid[idx2] < 0.0:
                ncmei_grid[idx2] = cmeiout_grid[idx2]
            else:
                pcmei_grid[idx2] = cmeiout_grid[idx2]
            if meltso_grid[idx2] < 0.0:
                nmelts_grid[idx2] = meltso_grid[idx2]
            else:
                pmelts_grid[idx2] = meltso_grid[idx2]
            sed_rates_grid[_idx3(i, k, iwtliq, pcols, pver)] = qcsedten_grid[idx2]
            sed_rates_grid[_idx3(i, k, iwtice, pcols, pver)] = qisedten_grid[idx2]

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            _add_process_rate(
                pre_rates_grid, i, k, iwtvap, iwtice, iwtvap, pcmei_grid[idx2], pcols, pver, pwtype
            )
            _add_process_rate(
                pre_rates_grid, i, k, iwtvap, iwtice, iwtice, ncmei_grid[idx2], pcols, pver, pwtype
            )
            _add_process_rate(
                pre_rates_grid, i, k, iwtvap, iwtstrain, iwtstrain, preo_grid[idx2], pcols, pver, pwtype
            )
            _add_process_rate(
                pre_rates_grid, i, k, iwtvap, iwtstsnow, iwtstsnow, prdso_grid[idx2], pcols, pver, pwtype
            )

            rate_val = mnuccco_grid[idx2] + mnuccto_grid[idx2]
            rate_val = rate_val + msacwio_grid[idx2]
            _add_process_rate(
                pre_rates_grid, i, k, iwtliq, iwtice, iwtliq, rate_val, pcols, pver, pwtype
            )

            rate_val = prao_grid[idx2] + prco_grid[idx2]
            _add_process_rate(
                pre_rates_grid, i, k, iwtliq, iwtstrain, iwtliq, rate_val, pcols, pver, pwtype
            )
            _add_process_rate(
                pre_rates_grid, i, k, iwtliq, iwtstsnow, iwtliq, psacwso_grid[idx2], pcols, pver, pwtype
            )
            _add_process_rate(
                pre_rates_grid, i, k, iwtliq, iwtliq, iwtliq, bergo_grid[idx2], pcols, pver, pwtype
            )
            _add_process_rate(
                pre_rates_grid, i, k, iwtice, iwtice, iwtice, bergso_grid[idx2], pcols, pver, pwtype
            )

            rate_val = praio_grid[idx2] + prcio_grid[idx2]
            _add_process_rate(
                pre_rates_grid, i, k, iwtice, iwtstsnow, iwtice, rate_val, pcols, pver, pwtype
            )

            rate_val = pracso_grid[idx2] + mnuccro_grid[idx2]
            _add_process_rate(
                pre_rates_grid, i, k, iwtstrain, iwtstsnow, iwtstrain, rate_val, pcols, pver, pwtype
            )

            _add_process_rate(
                post_rates_grid, i, k, iwtvap, iwtliq, iwtvap, qcreso_grid[idx2], pcols, pver, pwtype
            )
            _add_process_rate(
                post_rates_grid, i, k, iwtvap, iwtice, iwtvap, qireso_grid[idx2], pcols, pver, pwtype
            )
            _add_process_rate(
                post_rates_grid, i, k, iwtliq, iwtice, iwtliq, homoo_grid[idx2], pcols, pver, pwtype
            )
            _add_process_rate(
                post_rates_grid, i, k, iwtice, iwtliq, iwtice, melto_grid[idx2], pcols, pver, pwtype
            )


def _size_dist_param_basic_codon(
    qsmall: float,
    qic: float,
    nic: float,
    eff_dim: float,
    shape_coef: float,
    lambda_lo: float,
    lambda_hi: float,
    min_mean_mass: float,
):
    lam = 0.0
    nic_out = nic

    if qic > qsmall:
        nic_out = min(nic_out, qic / min_mean_mass)
        lam = (shape_coef * nic_out / qic) ** (1.0 / eff_dim)

        if lam < lambda_lo:
            lam = lambda_lo
            nic_out = lam**eff_dim * qic / shape_coef
        elif lam > lambda_hi:
            lam = lambda_hi
            nic_out = lam**eff_dim * qic / shape_coef

    return lam, nic_out


def _avg_diameter_codon(q: float, n: float, rho_air: float, rho_sub: float):
    return (pi * rho_sub * n / (q * rho_air)) ** (-1.0 / 3.0)


def _size_dist_param_liq_codon(
    qsmall: float,
    qcic: float,
    ncic: float,
    rho_air: float,
    liq_rho: float,
    liq_eff_dim: float,
    liq_min_mean_mass: float,
):
    pgam = -100.0
    lamc = 0.0
    ncic_out = ncic

    if qcic > qsmall:
        pgam = 0.0005714 * (ncic / 1.0e6 * rho_air) + 0.2714
        pgam = 1.0 / (pgam ** 2) - 1.0
        pgam = max(pgam, 2.0)
        pgam = min(pgam, 15.0)

        shape_coef = pi * liq_rho / 6.0 * (gamma(pgam + 1.0 + liq_eff_dim) / gamma(pgam + 1.0))
        lambda_lo = (pgam + 1.0) * 1.0 / 50.0e-6
        lambda_hi = (pgam + 1.0) * 1.0 / 2.0e-6
        lamc, ncic_out = _size_dist_param_basic_codon(
            qsmall,
            qcic,
            ncic_out,
            liq_eff_dim,
            shape_coef,
            lambda_lo,
            lambda_hi,
            liq_min_mean_mass,
        )

    return pgam, lamc, ncic_out


@export
def micro_mg_cam_postmg_diag_codon(
    ncol: int,
    psetcols: int,
    pver: int,
    pverp: int,
    top_lev: int,
    micro_mg_version: int,
    rate1_cw2pr_st_idx: int,
    ixcldliq: int,
    ixcldice: int,
    ixnumliq: int,
    ixnumice: int,
    ixrain: int,
    ixsnow: int,
    cpair: float,
    gravit: float,
    mincld: float,
    qsmall: float,
    state_q_p: cobj,
    state_t_p: cobj,
    state_pmid_p: cobj,
    state_pdel_p: cobj,
    naai_p: cobj,
    naai_hom_p: cobj,
    mnuccdo_p: cobj,
    rflx_p: cobj,
    sflx_p: cobj,
    qrout_p: cobj,
    qsout_p: cobj,
    prect_p: cobj,
    preci_p: cobj,
    rate1cld_p: cobj,
    vtrmc_p: cobj,
    tlat_p: cobj,
    qvlat_p: cobj,
    qcten_p: cobj,
    qiten_p: cobj,
    ncten_p: cobj,
    niten_p: cobj,
    alst_mic_p: cobj,
    cmeliq_p: cobj,
    cmeiout_p: cobj,
    ast_p: cobj,
    cld_p: cobj,
    concld_p: cobj,
    mnuccdohet_p: cobj,
    mgflxprc_p: cobj,
    mgflxsnw_p: cobj,
    mgmrprc_p: cobj,
    mgmrsnw_p: cobj,
    cvreffliq_p: cobj,
    cvreffice_p: cobj,
    rate1ord_cw2pr_st_p: cobj,
    wsedl_p: cobj,
    cc_t_p: cobj,
    cc_qv_p: cobj,
    cc_ql_p: cobj,
    cc_qi_p: cobj,
    cc_nl_p: cobj,
    cc_ni_p: cobj,
    cc_qlst_p: cobj,
    qme_p: cobj,
    prec_pcw_p: cobj,
    snow_pcw_p: cobj,
    prec_sed_p: cobj,
    snow_sed_p: cobj,
    prec_str_p: cobj,
    snow_str_p: cobj,
    icecldf_p: cobj,
    liqcldf_p: cobj,
    icinc_p: cobj,
    icwnc_p: cobj,
    iciwpst_p: cobj,
    iclwpst_p: cobj,
    icswp_p: cobj,
    cldfsnow_p: cobj,
    icimrst_p: cobj,
    icwmrst_p: cobj,
    cldmax_p: cobj,
):
    state_q = Ptr[float](state_q_p)
    state_t = Ptr[float](state_t_p)
    state_pmid = Ptr[float](state_pmid_p)
    state_pdel = Ptr[float](state_pdel_p)
    naai = Ptr[float](naai_p)
    naai_hom = Ptr[float](naai_hom_p)
    mnuccdo = Ptr[float](mnuccdo_p)
    rflx = Ptr[float](rflx_p)
    sflx = Ptr[float](sflx_p)
    qrout = Ptr[float](qrout_p)
    qsout = Ptr[float](qsout_p)
    prect = Ptr[float](prect_p)
    preci = Ptr[float](preci_p)
    rate1cld = Ptr[float](rate1cld_p)
    vtrmc = Ptr[float](vtrmc_p)
    tlat = Ptr[float](tlat_p)
    qvlat = Ptr[float](qvlat_p)
    qcten = Ptr[float](qcten_p)
    qiten = Ptr[float](qiten_p)
    ncten = Ptr[float](ncten_p)
    niten = Ptr[float](niten_p)
    alst_mic = Ptr[float](alst_mic_p)
    cmeliq = Ptr[float](cmeliq_p)
    cmeiout = Ptr[float](cmeiout_p)
    ast = Ptr[float](ast_p)
    cld = Ptr[float](cld_p)
    concld = Ptr[float](concld_p)
    mnuccdohet = Ptr[float](mnuccdohet_p)
    mgflxprc = Ptr[float](mgflxprc_p)
    mgflxsnw = Ptr[float](mgflxsnw_p)
    mgmrprc = Ptr[float](mgmrprc_p)
    mgmrsnw = Ptr[float](mgmrsnw_p)
    cvreffliq = Ptr[float](cvreffliq_p)
    cvreffice = Ptr[float](cvreffice_p)
    rate1ord_cw2pr_st = Ptr[float](rate1ord_cw2pr_st_p)
    wsedl = Ptr[float](wsedl_p)
    cc_t = Ptr[float](cc_t_p)
    cc_qv = Ptr[float](cc_qv_p)
    cc_ql = Ptr[float](cc_ql_p)
    cc_qi = Ptr[float](cc_qi_p)
    cc_nl = Ptr[float](cc_nl_p)
    cc_ni = Ptr[float](cc_ni_p)
    cc_qlst = Ptr[float](cc_qlst_p)
    qme = Ptr[float](qme_p)
    prec_pcw = Ptr[float](prec_pcw_p)
    snow_pcw = Ptr[float](snow_pcw_p)
    prec_sed = Ptr[float](prec_sed_p)
    snow_sed = Ptr[float](snow_sed_p)
    prec_str = Ptr[float](prec_str_p)
    snow_str = Ptr[float](snow_str_p)
    icecldf = Ptr[float](icecldf_p)
    liqcldf = Ptr[float](liqcldf_p)
    icinc = Ptr[float](icinc_p)
    icwnc = Ptr[float](icwnc_p)
    iciwpst = Ptr[float](iciwpst_p)
    iclwpst = Ptr[float](iclwpst_p)
    icswp = Ptr[float](icswp_p)
    cldfsnow = Ptr[float](cldfsnow_p)
    icimrst = Ptr[float](icimrst_p)
    icwmrst = Ptr[float](icwmrst_p)
    cldmax = Ptr[float](cldmax_p)

    for k in range(1, pver + 1):
        for i in range(1, psetcols + 1):
            idx2 = _idx2(i, k, psetcols)
            mnuccdohet[idx2] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, psetcols)
            if naai[idx2] > 0.0:
                mnuccdohet[idx2] = mnuccdo[idx2] - (naai_hom[idx2] / naai[idx2]) * mnuccdo[idx2]

    for k in range(top_lev, pverp + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, psetcols)
            mgflxprc[idx2] = rflx[idx2] + sflx[idx2]
            mgflxsnw[idx2] = sflx[idx2]

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, psetcols)
            mgmrprc[idx2] = qrout[idx2] + qsout[idx2]
            mgmrsnw[idx2] = qsout[idx2]

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, psetcols)
            cvreffliq[idx2] = 9.0
            cvreffice[idx2] = 37.0

    if rate1_cw2pr_st_idx > 0:
        for k in range(top_lev, pver + 1):
            for i in range(1, ncol + 1):
                idx2 = _idx2(i, k, psetcols)
                rate1ord_cw2pr_st[idx2] = rate1cld[idx2]

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, psetcols)
            wsedl[idx2] = vtrmc[idx2]
            cc_t[idx2] = tlat[idx2] / cpair
            cc_qv[idx2] = qvlat[idx2]
            cc_ql[idx2] = qcten[idx2]
            cc_qi[idx2] = qiten[idx2]
            cc_nl[idx2] = ncten[idx2]
            cc_ni[idx2] = niten[idx2]
            cc_qlst[idx2] = qcten[idx2] / max(0.01, alst_mic[idx2])
            qme[idx2] = cmeliq[idx2] + cmeiout[idx2]
            icecldf[idx2] = ast[idx2]
            liqcldf[idx2] = ast[idx2]

    for i in range(1, psetcols + 1):
        idx1 = i - 1
        prec_pcw[idx1] = prect[idx1]
        snow_pcw[idx1] = preci[idx1]
        prec_sed[idx1] = 0.0
        snow_sed[idx1] = 0.0
        prec_str[idx1] = prec_pcw[idx1] + prec_sed[idx1]
        snow_str[idx1] = snow_pcw[idx1] + snow_sed[idx1]

    for k in range(1, pver + 1):
        for i in range(1, psetcols + 1):
            idx2 = _idx2(i, k, psetcols)
            icinc[idx2] = 0.0
            icwnc[idx2] = 0.0
            iciwpst[idx2] = 0.0
            iclwpst[idx2] = 0.0
            icswp[idx2] = 0.0
            cldfsnow[idx2] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, psetcols)
            idx_q_liq = _idx3(i, k, ixcldliq, psetcols, pver)
            idx_q_ice = _idx3(i, k, ixcldice, psetcols, pver)
            idx_n_liq = _idx3(i, k, ixnumliq, psetcols, pver)
            idx_n_ice = _idx3(i, k, ixnumice, psetcols, pver)

            icimrst[idx2] = min(state_q[idx_q_ice] / max(mincld, icecldf[idx2]), 0.005)
            icwmrst[idx2] = min(state_q[idx_q_liq] / max(mincld, liqcldf[idx2]), 0.005)
            icinc[idx2] = state_q[idx_n_ice] / max(mincld, icecldf[idx2]) * state_pmid[idx2] / (287.15 * state_t[idx2])
            icwnc[idx2] = state_q[idx_n_liq] / max(mincld, liqcldf[idx2]) * state_pmid[idx2] / (287.15 * state_t[idx2])
            iciwpst[idx2] = min(state_q[idx_q_ice] / max(mincld, ast[idx2]), 0.005) * state_pdel[idx2] / gravit
            iclwpst[idx2] = min(state_q[idx_q_liq] / max(mincld, ast[idx2]), 0.005) * state_pdel[idx2] / gravit

            cldfsnow[idx2] = cld[idx2]
            if cldfsnow[idx2] > 1.0e-4 and concld[idx2] < 1.0e-4 and state_q[idx_q_liq] < 1.0e-10:
                cldfsnow[idx2] = 0.0
            if cldfsnow[idx2] <= 1.0e-4 and qsout[idx2] > 1.0e-6:
                cldfsnow[idx2] = 0.25

            icswp[idx2] = qsout[idx2] / max(mincld, cldfsnow[idx2]) * state_pdel[idx2] / gravit

    if micro_mg_version > 1:
        for k in range(1, pver + 1):
            for i in range(1, psetcols + 1):
                idx2 = _idx2(i, k, psetcols)
                cldmax[idx2] = max(mincld, ast[idx2])

        for k in range(top_lev + 1, pver + 1):
            for i in range(1, ncol + 1):
                idx2 = _idx2(i, k, psetcols)
                idx2_prev = _idx2(i, k - 1, psetcols)
                idx_rain_prev = _idx3(i, k - 1, ixrain, psetcols, pver)
                idx_snow_prev = _idx3(i, k - 1, ixsnow, psetcols, pver)
                if state_q[idx_rain_prev] >= qsmall or state_q[idx_snow_prev] >= qsmall:
                    cldmax[idx2] = max(cldmax[idx2_prev], cldmax[idx2])


@export
def micro_mg_cam_grid_diag_codon(
    ngrdcol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    minlwp: float,
    gravit: float,
    rhoh2o: float,
    iclwpst_grid_p: cobj,
    cld_grid_p: cobj,
    cmeliq_grid_p: cobj,
    pdel_grid_p: cobj,
    prec_str_grid_p: cobj,
    acgcme_grid_p: cobj,
    acprecl_grid_p: cobj,
    acnum_grid_p: cobj,
    prao_grid_p: cobj,
    prco_grid_p: cobj,
    nc_grid_p: cobj,
    liqcldf_grid_p: cobj,
    icwmrst_grid_p: cobj,
    rel_grid_p: cobj,
    icwnc_grid_p: cobj,
    icecldf_grid_p: cobj,
    icimrst_grid_p: cobj,
    rei_grid_p: cobj,
    icinc_grid_p: cobj,
    nevapr_grid_p: cobj,
    evpsnow_st_grid_p: cobj,
    tgliqwp_grid_p: cobj,
    tgcmeliq_grid_p: cobj,
    pe_grid_p: cobj,
    tpr_grid_p: cobj,
    pefrac_grid_p: cobj,
    vprao_grid_p: cobj,
    vprco_grid_p: cobj,
    racau_grid_p: cobj,
    cnt_grid_p: cobj,
    cdnumc_grid_p: cobj,
    efcout_grid_p: cobj,
    efiout_grid_p: cobj,
    ncout_grid_p: cobj,
    niout_grid_p: cobj,
    freql_grid_p: cobj,
    freqi_grid_p: cobj,
    icwmrst_grid_out_p: cobj,
    icimrst_grid_out_p: cobj,
    fcti_grid_p: cobj,
    fctl_grid_p: cobj,
    ctrel_grid_p: cobj,
    ctrei_grid_p: cobj,
    ctnl_grid_p: cobj,
    ctni_grid_p: cobj,
    evprain_st_grid_p: cobj,
):
    iclwpst_grid = Ptr[float](iclwpst_grid_p)
    cld_grid = Ptr[float](cld_grid_p)
    cmeliq_grid = Ptr[float](cmeliq_grid_p)
    pdel_grid = Ptr[float](pdel_grid_p)
    prec_str_grid = Ptr[float](prec_str_grid_p)
    acgcme_grid = Ptr[float](acgcme_grid_p)
    acprecl_grid = Ptr[float](acprecl_grid_p)
    acnum_grid = Ptr[i32](acnum_grid_p)
    prao_grid = Ptr[float](prao_grid_p)
    prco_grid = Ptr[float](prco_grid_p)
    nc_grid = Ptr[float](nc_grid_p)
    liqcldf_grid = Ptr[float](liqcldf_grid_p)
    icwmrst_grid = Ptr[float](icwmrst_grid_p)
    rel_grid = Ptr[float](rel_grid_p)
    icwnc_grid = Ptr[float](icwnc_grid_p)
    icecldf_grid = Ptr[float](icecldf_grid_p)
    icimrst_grid = Ptr[float](icimrst_grid_p)
    rei_grid = Ptr[float](rei_grid_p)
    icinc_grid = Ptr[float](icinc_grid_p)
    nevapr_grid = Ptr[float](nevapr_grid_p)
    evpsnow_st_grid = Ptr[float](evpsnow_st_grid_p)
    tgliqwp_grid = Ptr[float](tgliqwp_grid_p)
    tgcmeliq_grid = Ptr[float](tgcmeliq_grid_p)
    pe_grid = Ptr[float](pe_grid_p)
    tpr_grid = Ptr[float](tpr_grid_p)
    pefrac_grid = Ptr[float](pefrac_grid_p)
    vprao_grid = Ptr[float](vprao_grid_p)
    vprco_grid = Ptr[float](vprco_grid_p)
    racau_grid = Ptr[float](racau_grid_p)
    cnt_grid = Ptr[i32](cnt_grid_p)
    cdnumc_grid = Ptr[float](cdnumc_grid_p)
    efcout_grid = Ptr[float](efcout_grid_p)
    efiout_grid = Ptr[float](efiout_grid_p)
    ncout_grid = Ptr[float](ncout_grid_p)
    niout_grid = Ptr[float](niout_grid_p)
    freql_grid = Ptr[float](freql_grid_p)
    freqi_grid = Ptr[float](freqi_grid_p)
    icwmrst_grid_out = Ptr[float](icwmrst_grid_out_p)
    icimrst_grid_out = Ptr[float](icimrst_grid_out_p)
    fcti_grid = Ptr[float](fcti_grid_p)
    fctl_grid = Ptr[float](fctl_grid_p)
    ctrel_grid = Ptr[float](ctrel_grid_p)
    ctrei_grid = Ptr[float](ctrei_grid_p)
    ctnl_grid = Ptr[float](ctnl_grid_p)
    ctni_grid = Ptr[float](ctni_grid_p)
    evprain_st_grid = Ptr[float](evprain_st_grid_p)

    for i in range(1, ngrdcol + 1):
        idx1 = i - 1
        tgliqwp_grid[idx1] = 0.0
        tgcmeliq_grid[idx1] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ngrdcol + 1):
            idx1 = i - 1
            idx2 = _idx2(i, k, pcols)
            tgliqwp_grid[idx1] = tgliqwp_grid[idx1] + iclwpst_grid[idx2] * cld_grid[idx2]
            if cmeliq_grid[idx2] > 1.0e-12:
                tgcmeliq_grid[idx1] = tgcmeliq_grid[idx1] + cmeliq_grid[idx2] * (pdel_grid[idx2] / gravit) / rhoh2o

    for i in range(1, ngrdcol + 1):
        idx1 = i - 1
        pe_grid[idx1] = 0.0
        tpr_grid[idx1] = 0.0
        pefrac_grid[idx1] = 0.0

    for i in range(1, ngrdcol + 1):
        idx1 = i - 1
        acgcme_grid[idx1] = acgcme_grid[idx1] + tgcmeliq_grid[idx1]
        acprecl_grid[idx1] = acprecl_grid[idx1] + prec_str_grid[idx1]
        acnum_grid[idx1] = i32(int(acnum_grid[idx1]) + 1)

        if tgliqwp_grid[idx1] < minlwp:
            if acprecl_grid[idx1] > 5.0e-8:
                tpr_grid[idx1] = max(acprecl_grid[idx1] / int(acnum_grid[idx1]), 1.0e-15)
                if acgcme_grid[idx1] > 1.0e-10:
                    pe_grid[idx1] = min(max(acprecl_grid[idx1] / acgcme_grid[idx1], 1.0e-15), 1.0e5)
                    pefrac_grid[idx1] = 1.0

            acprecl_grid[idx1] = 0.0
            acgcme_grid[idx1] = 0.0
            acnum_grid[idx1] = i32(0)

        if int(acnum_grid[idx1]) > 1000:
            acnum_grid[idx1] = i32(0)
            acprecl_grid[idx1] = 0.0
            acgcme_grid[idx1] = 0.0

    for i in range(1, pcols + 1):
        idx1 = i - 1
        vprao_grid[idx1] = 0.0
        cnt_grid[idx1] = i32(0)

    for k in range(top_lev, pver + 1):
        for i in range(1, ngrdcol + 1):
            idx1 = i - 1
            idx2 = _idx2(i, k, pcols)
            vprao_grid[idx1] = vprao_grid[idx1] + prao_grid[idx2]
            if prao_grid[idx2] != 0.0:
                cnt_grid[idx1] = i32(int(cnt_grid[idx1]) + 1)

    for i in range(1, pcols + 1):
        idx1 = i - 1
        if int(cnt_grid[idx1]) > 0:
            vprao_grid[idx1] = vprao_grid[idx1] / int(cnt_grid[idx1])

    for i in range(1, pcols + 1):
        idx1 = i - 1
        vprco_grid[idx1] = 0.0
        cnt_grid[idx1] = i32(0)

    for k in range(top_lev, pver + 1):
        for i in range(1, ngrdcol + 1):
            idx1 = i - 1
            idx2 = _idx2(i, k, pcols)
            vprco_grid[idx1] = vprco_grid[idx1] + prco_grid[idx2]
            if prco_grid[idx2] != 0.0:
                cnt_grid[idx1] = i32(int(cnt_grid[idx1]) + 1)

    for i in range(1, pcols + 1):
        idx1 = i - 1
        if int(cnt_grid[idx1]) > 0:
            vprco_grid[idx1] = vprco_grid[idx1] / int(cnt_grid[idx1])
            racau_grid[idx1] = vprao_grid[idx1] / vprco_grid[idx1]
        else:
            racau_grid[idx1] = 0.0

    for i in range(1, pcols + 1):
        idx1 = i - 1
        racau_grid[idx1] = min(racau_grid[idx1], 1.0e10)

    for i in range(1, ngrdcol + 1):
        idx1 = i - 1
        cdnumc_grid[idx1] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ngrdcol + 1):
            idx1 = i - 1
            idx2 = _idx2(i, k, pcols)
            cdnumc_grid[idx1] = cdnumc_grid[idx1] + nc_grid[idx2] * pdel_grid[idx2] / gravit

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            idx2 = _idx2(i, k, pcols)
            efcout_grid[idx2] = 0.0
            efiout_grid[idx2] = 0.0
            ncout_grid[idx2] = 0.0
            niout_grid[idx2] = 0.0
            freql_grid[idx2] = 0.0
            freqi_grid[idx2] = 0.0
            icwmrst_grid_out[idx2] = 0.0
            icimrst_grid_out[idx2] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ngrdcol + 1):
            idx2 = _idx2(i, k, pcols)
            if liqcldf_grid[idx2] > 0.01 and icwmrst_grid[idx2] > 5.0e-5:
                efcout_grid[idx2] = rel_grid[idx2] * liqcldf_grid[idx2]
                ncout_grid[idx2] = icwnc_grid[idx2] * liqcldf_grid[idx2]
                freql_grid[idx2] = liqcldf_grid[idx2]
                icwmrst_grid_out[idx2] = icwmrst_grid[idx2]
            if icecldf_grid[idx2] > 0.01 and icimrst_grid[idx2] > 1.0e-6:
                efiout_grid[idx2] = rei_grid[idx2] * icecldf_grid[idx2]
                niout_grid[idx2] = icinc_grid[idx2] * icecldf_grid[idx2]
                freqi_grid[idx2] = icecldf_grid[idx2]
                icimrst_grid_out[idx2] = icimrst_grid[idx2]

    for i in range(1, pcols + 1):
        idx1 = i - 1
        fcti_grid[idx1] = 0.0
        fctl_grid[idx1] = 0.0
        ctrel_grid[idx1] = 0.0
        ctrei_grid[idx1] = 0.0
        ctnl_grid[idx1] = 0.0
        ctni_grid[idx1] = 0.0

    for i in range(1, ngrdcol + 1):
        idx1 = i - 1
        for k in range(top_lev, pver + 1):
            idx2 = _idx2(i, k, pcols)
            if liqcldf_grid[idx2] > 0.01 and icwmrst_grid[idx2] > 1.0e-7:
                ctrel_grid[idx1] = rel_grid[idx2] * liqcldf_grid[idx2]
                ctnl_grid[idx1] = icwnc_grid[idx2] * liqcldf_grid[idx2]
                fctl_grid[idx1] = liqcldf_grid[idx2]
                break
            if icecldf_grid[idx2] > 0.01 and icimrst_grid[idx2] > 1.0e-7:
                ctrei_grid[idx1] = rei_grid[idx2] * icecldf_grid[idx2]
                ctni_grid[idx1] = icinc_grid[idx2] * icecldf_grid[idx2]
                fcti_grid[idx1] = icecldf_grid[idx2]
                break

    for k in range(1, pver + 1):
        for i in range(1, ngrdcol + 1):
            idx2 = _idx2(i, k, pcols)
            evprain_st_grid[idx2] = nevapr_grid[idx2] - evpsnow_st_grid[idx2]

    for k in range(top_lev, pver + 1):
        for i in range(1, ngrdcol + 1):
            idx2 = _idx2(i, k, pcols)
            evprain_st_grid[idx2] = max(evprain_st_grid[idx2], 0.0)
            evpsnow_st_grid[idx2] = max(evpsnow_st_grid[idx2], 0.0)


@export
def micro_mg_cam_wtrc_prep_codon(
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    pwtype: int,
    iwtliq: int,
    iwtice: int,
    cmeiout_grid_p: cobj,
    meltso_grid_p: cobj,
    qcsedten_grid_p: cobj,
    qisedten_grid_p: cobj,
    pcmei_grid_p: cobj,
    ncmei_grid_p: cobj,
    pmelts_grid_p: cobj,
    nmelts_grid_p: cobj,
    sed_rates_grid_p: cobj,
):
    cmeiout_grid = Ptr[float](cmeiout_grid_p)
    meltso_grid = Ptr[float](meltso_grid_p)
    qcsedten_grid = Ptr[float](qcsedten_grid_p)
    qisedten_grid = Ptr[float](qisedten_grid_p)
    pcmei_grid = Ptr[float](pcmei_grid_p)
    ncmei_grid = Ptr[float](ncmei_grid_p)
    pmelts_grid = Ptr[float](pmelts_grid_p)
    nmelts_grid = Ptr[float](nmelts_grid_p)
    sed_rates_grid = Ptr[float](sed_rates_grid_p)

    for m in range(1, pwtype + 1):
        for k in range(top_lev, pver + 1):
            for i in range(1, pcols + 1):
                sed_rates_grid[_idx3(i, k, m, pcols, pver)] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, pcols + 1):
            idx2 = _idx2(i, k, pcols)
            pcmei_grid[idx2] = 0.0
            ncmei_grid[idx2] = 0.0
            pmelts_grid[idx2] = 0.0
            nmelts_grid[idx2] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)
            if cmeiout_grid[idx2] < 0.0:
                ncmei_grid[idx2] = cmeiout_grid[idx2]
            else:
                pcmei_grid[idx2] = cmeiout_grid[idx2]

            if meltso_grid[idx2] < 0.0:
                nmelts_grid[idx2] = meltso_grid[idx2]
            else:
                pmelts_grid[idx2] = meltso_grid[idx2]

            sed_rates_grid[_idx3(i, k, iwtliq, pcols, pver)] = qcsedten_grid[idx2]
            sed_rates_grid[_idx3(i, k, iwtice, pcols, pver)] = qisedten_grid[idx2]


@export
def micro_mg_cam_budget_diag_codon(
    mode: int,
    ncol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    qcreso_grid_p: cobj,
    melto_grid_p: cobj,
    mnuccco_grid_p: cobj,
    mnuccto_grid_p: cobj,
    bergo_grid_p: cobj,
    homoo_grid_p: cobj,
    msacwio_grid_p: cobj,
    prao_grid_p: cobj,
    prco_grid_p: cobj,
    psacwso_grid_p: cobj,
    bergso_grid_p: cobj,
    cmeiout_grid_p: cobj,
    qireso_grid_p: cobj,
    prcio_grid_p: cobj,
    praio_grid_p: cobj,
    ftem_grid_p: cobj,
):
    qcreso_grid = Ptr[float](qcreso_grid_p)
    melto_grid = Ptr[float](melto_grid_p)
    mnuccco_grid = Ptr[float](mnuccco_grid_p)
    mnuccto_grid = Ptr[float](mnuccto_grid_p)
    bergo_grid = Ptr[float](bergo_grid_p)
    homoo_grid = Ptr[float](homoo_grid_p)
    msacwio_grid = Ptr[float](msacwio_grid_p)
    prao_grid = Ptr[float](prao_grid_p)
    prco_grid = Ptr[float](prco_grid_p)
    psacwso_grid = Ptr[float](psacwso_grid_p)
    bergso_grid = Ptr[float](bergso_grid_p)
    cmeiout_grid = Ptr[float](cmeiout_grid_p)
    qireso_grid = Ptr[float](qireso_grid_p)
    prcio_grid = Ptr[float](prcio_grid_p)
    praio_grid = Ptr[float](praio_grid_p)
    ftem_grid = Ptr[float](ftem_grid_p)

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            ftem_grid[_idx2(i, k, pcols)] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            idx2 = _idx2(i, k, pcols)

            if mode == 1:
                ftem_grid[idx2] = qcreso_grid[idx2]
            elif mode == 2:
                tmp = melto_grid[idx2] - mnuccco_grid[idx2]
                tmp = tmp - mnuccto_grid[idx2]
                tmp = tmp - bergo_grid[idx2]
                tmp = tmp - homoo_grid[idx2]
                tmp = tmp - msacwio_grid[idx2]
                ftem_grid[idx2] = tmp
            elif mode == 3:
                tmp = -prao_grid[idx2]
                tmp = tmp - prco_grid[idx2]
                tmp = tmp - psacwso_grid[idx2]
                tmp = tmp - bergso_grid[idx2]
                ftem_grid[idx2] = tmp
            elif mode == 4:
                ftem_grid[idx2] = cmeiout_grid[idx2] + qireso_grid[idx2]
            elif mode == 5:
                tmp = -melto_grid[idx2] + mnuccco_grid[idx2]
                tmp = tmp + mnuccto_grid[idx2]
                tmp = tmp + bergo_grid[idx2]
                tmp = tmp + homoo_grid[idx2]
                tmp = tmp + msacwio_grid[idx2]
                ftem_grid[idx2] = tmp
            elif mode == 6:
                ftem_grid[idx2] = -prcio_grid[idx2] - praio_grid[idx2]


def _micro_mg_cam_pbuf_copy_fields(
    ncol: int,
    pcols: int,
    pver: int,
    copy_qrain: int,
    copy_qsnow: int,
    copy_nrain: int,
    copy_nsnow: int,
    qrout_grid_p: cobj,
    qsout_grid_p: cobj,
    nrout_grid_p: cobj,
    nsout_grid_p: cobj,
    qrout_grid_ptr_p: cobj,
    qsout_grid_ptr_p: cobj,
    nrout_grid_ptr_p: cobj,
    nsout_grid_ptr_p: cobj,
):
    qrout_grid = Ptr[float](qrout_grid_p)
    qsout_grid = Ptr[float](qsout_grid_p)
    nrout_grid = Ptr[float](nrout_grid_p)
    nsout_grid = Ptr[float](nsout_grid_p)
    qrout_grid_ptr = Ptr[float](qrout_grid_ptr_p)
    qsout_grid_ptr = Ptr[float](qsout_grid_ptr_p)
    nrout_grid_ptr = Ptr[float](nrout_grid_ptr_p)
    nsout_grid_ptr = Ptr[float](nsout_grid_ptr_p)

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            idx2 = _idx2(i, k, pcols)
            if copy_qrain == 1:
                qrout_grid_ptr[idx2] = qrout_grid[idx2]
            if copy_qsnow == 1:
                qsout_grid_ptr[idx2] = qsout_grid[idx2]
            if copy_nrain == 1:
                nrout_grid_ptr[idx2] = nrout_grid[idx2]
            if copy_nsnow == 1:
                nsout_grid_ptr[idx2] = nsout_grid[idx2]


@export
def micro_mg_cam_pbuf_copy_codon(
    ncol: int,
    pcols: int,
    pver: int,
    copy_qrain: int,
    copy_qsnow: int,
    copy_nrain: int,
    copy_nsnow: int,
    qrout_grid_p: cobj,
    qsout_grid_p: cobj,
    nrout_grid_p: cobj,
    nsout_grid_p: cobj,
    qrout_grid_ptr_p: cobj,
    qsout_grid_ptr_p: cobj,
    nrout_grid_ptr_p: cobj,
    nsout_grid_ptr_p: cobj,
):
    _micro_mg_cam_pbuf_copy_fields(
        ncol,
        pcols,
        pver,
        copy_qrain,
        copy_qsnow,
        copy_nrain,
        copy_nsnow,
        qrout_grid_p,
        qsout_grid_p,
        nrout_grid_p,
        nsout_grid_p,
        qrout_grid_ptr_p,
        qsout_grid_ptr_p,
        nrout_grid_ptr_p,
        nsout_grid_ptr_p,
    )


@export
def micro_mg_cam_tail_pbuf_copy_stage_dispatch_codon(
    ncol: int,
    pcols: int,
    pver: int,
    copy_qrain: int,
    copy_qsnow: int,
    copy_nrain: int,
    copy_nsnow: int,
    qrout_grid_p: cobj,
    qsout_grid_p: cobj,
    nrout_grid_p: cobj,
    nsout_grid_p: cobj,
    qrout_grid_ptr_p: cobj,
    qsout_grid_ptr_p: cobj,
    nrout_grid_ptr_p: cobj,
    nsout_grid_ptr_p: cobj,
):
    _micro_mg_cam_pbuf_copy_fields(
        ncol,
        pcols,
        pver,
        copy_qrain,
        copy_qsnow,
        copy_nrain,
        copy_nsnow,
        qrout_grid_p,
        qsout_grid_p,
        nrout_grid_p,
        nsout_grid_p,
        qrout_grid_ptr_p,
        qsout_grid_ptr_p,
        nrout_grid_ptr_p,
        nsout_grid_ptr_p,
    )

@export
def micro_mg_cam_tail_pbuf_copy_codon(
    ncol: int,
    pcols: int,
    pver: int,
    copy_qrain: int,
    copy_qsnow: int,
    copy_nrain: int,
    copy_nsnow: int,
    qrout_grid_p: cobj,
    qsout_grid_p: cobj,
    nrout_grid_p: cobj,
    nsout_grid_p: cobj,
    qrout_grid_ptr_p: cobj,
    qsout_grid_ptr_p: cobj,
    nrout_grid_ptr_p: cobj,
    nsout_grid_ptr_p: cobj,
):
    micro_mg_cam_tail_pbuf_copy_stage_dispatch_codon(
        ncol,
        pcols,
        pver,
        copy_qrain,
        copy_qsnow,
        copy_nrain,
        copy_nsnow,
        qrout_grid_p,
        qsout_grid_p,
        nrout_grid_p,
        nsout_grid_p,
        qrout_grid_ptr_p,
        qsout_grid_ptr_p,
        nrout_grid_ptr_p,
        nsout_grid_ptr_p,
    )


@export
def micro_mg_cam_reff_calc_codon(
    ngrdcol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    micro_mg_version: int,
    qsmall: float,
    mincld: float,
    liq_rho: float,
    liq_eff_dim: float,
    liq_min_mean_mass: float,
    ice_eff_dim: float,
    ice_shape_coef: float,
    ice_lambda_lo: float,
    ice_lambda_hi: float,
    ice_min_mean_mass: float,
    rhosn: float,
    rhoi: float,
    rhow: float,
    rhows: float,
    mucon: float,
    dcon: float,
    deicon: float,
    rho_grid_p: cobj,
    icwmrst_grid_p: cobj,
    liqcldf_grid_p: cobj,
    nc_grid_p: cobj,
    qr_grid_p: cobj,
    nr_grid_p: cobj,
    qs_grid_p: cobj,
    ns_grid_p: cobj,
    qrout_grid_p: cobj,
    nrout_grid_p: cobj,
    qsout_grid_p: cobj,
    nsout_grid_p: cobj,
    ni_grid_p: cobj,
    icecldf_grid_p: cobj,
    icimrst_grid_p: cobj,
    ast_grid_p: cobj,
    mu_grid_p: cobj,
    lambdac_grid_p: cobj,
    rel_fn_grid_p: cobj,
    ncic_grid_p: cobj,
    rel_grid_p: cobj,
    drout2_grid_p: cobj,
    reff_rain_grid_p: cobj,
    des_grid_p: cobj,
    dsout2_grid_p: cobj,
    reff_snow_grid_p: cobj,
    rei_grid_p: cobj,
    niic_grid_p: cobj,
    dei_grid_p: cobj,
    mgreffrain_grid_p: cobj,
    mgreffsnow_grid_p: cobj,
):
    rho_grid = Ptr[float](rho_grid_p)
    icwmrst_grid = Ptr[float](icwmrst_grid_p)
    liqcldf_grid = Ptr[float](liqcldf_grid_p)
    nc_grid = Ptr[float](nc_grid_p)
    qr_grid = Ptr[float](qr_grid_p)
    nr_grid = Ptr[float](nr_grid_p)
    qs_grid = Ptr[float](qs_grid_p)
    ns_grid = Ptr[float](ns_grid_p)
    qrout_grid = Ptr[float](qrout_grid_p)
    nrout_grid = Ptr[float](nrout_grid_p)
    qsout_grid = Ptr[float](qsout_grid_p)
    nsout_grid = Ptr[float](nsout_grid_p)
    ni_grid = Ptr[float](ni_grid_p)
    icecldf_grid = Ptr[float](icecldf_grid_p)
    icimrst_grid = Ptr[float](icimrst_grid_p)
    ast_grid = Ptr[float](ast_grid_p)
    mu_grid = Ptr[float](mu_grid_p)
    lambdac_grid = Ptr[float](lambdac_grid_p)
    rel_fn_grid = Ptr[float](rel_fn_grid_p)
    ncic_grid = Ptr[float](ncic_grid_p)
    rel_grid = Ptr[float](rel_grid_p)
    drout2_grid = Ptr[float](drout2_grid_p)
    reff_rain_grid = Ptr[float](reff_rain_grid_p)
    des_grid = Ptr[float](des_grid_p)
    dsout2_grid = Ptr[float](dsout2_grid_p)
    reff_snow_grid = Ptr[float](reff_snow_grid_p)
    rei_grid = Ptr[float](rei_grid_p)
    niic_grid = Ptr[float](niic_grid_p)
    dei_grid = Ptr[float](dei_grid_p)
    mgreffrain_grid = Ptr[float](mgreffrain_grid_p)
    mgreffsnow_grid = Ptr[float](mgreffsnow_grid_p)

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            idx2 = _idx2(i, k, pcols)
            mu_grid[idx2] = 0.0
            lambdac_grid[idx2] = 0.0
            rel_fn_grid[idx2] = 10.0
            ncic_grid[idx2] = 1.0e8

    for k in range(top_lev, pver + 1):
        for i in range(1, ngrdcol + 1):
            idx2 = _idx2(i, k, pcols)
            mu_grid[idx2], lambdac_grid[idx2], ncic_grid[idx2] = _size_dist_param_liq_codon(
                qsmall,
                icwmrst_grid[idx2],
                ncic_grid[idx2],
                rho_grid[idx2],
                liq_rho,
                liq_eff_dim,
                liq_min_mean_mass,
            )
            if icwmrst_grid[idx2] > qsmall:
                rel_fn_grid[idx2] = (mu_grid[idx2] + 3.0) / lambdac_grid[idx2] / 2.0 * 1.0e6

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            idx2 = _idx2(i, k, pcols)
            mu_grid[idx2] = 0.0
            lambdac_grid[idx2] = 0.0
            rel_grid[idx2] = 10.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ngrdcol + 1):
            idx2 = _idx2(i, k, pcols)
            ncic_grid[idx2] = nc_grid[idx2] / max(mincld, liqcldf_grid[idx2])
            mu_grid[idx2], lambdac_grid[idx2], ncic_grid[idx2] = _size_dist_param_liq_codon(
                qsmall,
                icwmrst_grid[idx2],
                ncic_grid[idx2],
                rho_grid[idx2],
                liq_rho,
                liq_eff_dim,
                liq_min_mean_mass,
            )
            if icwmrst_grid[idx2] >= qsmall:
                rel_grid[idx2] = (mu_grid[idx2] + 3.0) / lambdac_grid[idx2] / 2.0 * 1.0e6
            else:
                mu_grid[idx2] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            idx2 = _idx2(i, k, pcols)
            drout2_grid[idx2] = 0.0
            reff_rain_grid[idx2] = 0.0
            des_grid[idx2] = 0.0
            dsout2_grid[idx2] = 0.0
            reff_snow_grid[idx2] = 0.0

    if micro_mg_version > 1:
        for k in range(top_lev, pver + 1):
            for i in range(1, ngrdcol + 1):
                idx2 = _idx2(i, k, pcols)
                if qr_grid[idx2] >= 1.0e-7:
                    drout2_grid[idx2] = _avg_diameter_codon(
                        qr_grid[idx2],
                        nr_grid[idx2] * rho_grid[idx2],
                        rho_grid[idx2],
                        rhow,
                    )
                    reff_rain_grid[idx2] = drout2_grid[idx2] * 1.5 * 1.0e6

                if qs_grid[idx2] >= 1.0e-7:
                    dsout2_grid[idx2] = _avg_diameter_codon(
                        qs_grid[idx2],
                        ns_grid[idx2] * rho_grid[idx2],
                        rho_grid[idx2],
                        rhosn,
                    )
                    des_grid[idx2] = dsout2_grid[idx2] * 3.0 * rhosn / rhows
                    reff_snow_grid[idx2] = dsout2_grid[idx2] * 1.5 * 1.0e6
    else:
        for k in range(top_lev, pver + 1):
            for i in range(1, ngrdcol + 1):
                idx2 = _idx2(i, k, pcols)
                if qrout_grid[idx2] >= 1.0e-7:
                    drout2_grid[idx2] = _avg_diameter_codon(
                        qrout_grid[idx2],
                        nrout_grid[idx2] * rho_grid[idx2],
                        rho_grid[idx2],
                        rhow,
                    )
                    reff_rain_grid[idx2] = drout2_grid[idx2] * 1.5 * 1.0e6

                if qsout_grid[idx2] >= 1.0e-7:
                    dsout2_grid[idx2] = _avg_diameter_codon(
                        qsout_grid[idx2],
                        nsout_grid[idx2] * rho_grid[idx2],
                        rho_grid[idx2],
                        rhosn,
                    )
                    des_grid[idx2] = dsout2_grid[idx2] * 3.0 * rhosn / rhows
                    reff_snow_grid[idx2] = dsout2_grid[idx2] * 1.5 * 1.0e6

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            idx2 = _idx2(i, k, pcols)
            rei_grid[idx2] = 25.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ngrdcol + 1):
            idx2 = _idx2(i, k, pcols)
            niic_grid[idx2] = ni_grid[idx2] / max(mincld, icecldf_grid[idx2])
            rei_grid[idx2], niic_grid[idx2] = _size_dist_param_basic_codon(
                qsmall,
                icimrst_grid[idx2],
                niic_grid[idx2],
                ice_eff_dim,
                ice_shape_coef,
                ice_lambda_lo,
                ice_lambda_hi,
                ice_min_mean_mass,
            )
            if icimrst_grid[idx2] >= qsmall:
                rei_grid[idx2] = 1.5 / rei_grid[idx2] * 1.0e6
            else:
                rei_grid[idx2] = 25.0

    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            idx2 = _idx2(i, k, pcols)
            dei_grid[idx2] = rei_grid[idx2] * rhoi / rhows * 2.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ngrdcol + 1):
            idx2 = _idx2(i, k, pcols)
            des_grid[idx2] = des_grid[idx2] * 1.0e6
            if ast_grid[idx2] < 1.0e-4:
                mu_grid[idx2] = mucon
                lambdac_grid[idx2] = (mucon + 1.0) / dcon
                dei_grid[idx2] = deicon

    for k in range(top_lev, pver + 1):
        for i in range(1, ngrdcol + 1):
            idx2 = _idx2(i, k, pcols)
            mgreffrain_grid[idx2] = reff_rain_grid[idx2]
            mgreffsnow_grid[idx2] = reff_snow_grid[idx2]


@export
def micro_mg_cam_rho_grid_codon(
    ncol: int,
    pcols: int,
    pver: int,
    psetcols: int,
    top_lev: int,
    rair: float,
    rho_p: cobj,
    pmid_p: cobj,
    t_p: cobj,
    rho_grid_p: cobj,
):
    rho = Ptr[float](rho_p)
    pmid = Ptr[float](pmid_p)
    t = Ptr[float](t_p)
    rho_grid = Ptr[float](rho_grid_p)

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            rho[_idx2(i, k, psetcols)] = pmid[_idx2(i, k, psetcols)] / (
                rair * t[_idx2(i, k, psetcols)]
            )

    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            rho_grid[_idx2(i, k, pcols)] = rho[_idx2(i, k, psetcols)]


@inline
def _copy2d_full(src: Ptr[float], dst: Ptr[float], src_ld: int, dst_ld: int, pver: int):
    for k in range(1, pver + 1):
        for i in range(1, dst_ld + 1):
            dst[_idx2(i, k, dst_ld)] = src[_idx2(i, k, src_ld)]


def micro_mg_cam_tail_grid_copy_codon(
    psetcols: int,
    pcols: int,
    pver: int,
    copy_mg10: int,
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
    p35: cobj,
    p36: cobj,
    p37: cobj,
    p38: cobj,
    p39: cobj,
    p40: cobj,
    p41: cobj,
    p42: cobj,
    p43: cobj,
    p44: cobj,
    p45: cobj,
    p46: cobj,
    p47: cobj,
    p48: cobj,
    p49: cobj,
    p50: cobj,
    p51: cobj,
    p52: cobj,
    p53: cobj,
    p54: cobj,
    p55: cobj,
    p56: cobj,
):
    if copy_mg10 != 0:
        _copy2d_full(Ptr[float](p1), Ptr[float](p2), psetcols, pcols, pver)

    _copy2d_full(Ptr[float](p3), Ptr[float](p4), psetcols, pcols, pver)
    _copy2d_full(Ptr[float](p5), Ptr[float](p6), psetcols, pcols, pver)
    _copy2d_full(Ptr[float](p7), Ptr[float](p8), psetcols, pcols, pver)
    _copy2d_full(Ptr[float](p9), Ptr[float](p10), psetcols, pcols, pver)
    _copy2d_full(Ptr[float](p11), Ptr[float](p12), psetcols, pcols, pver)
    _copy2d_full(Ptr[float](p13), Ptr[float](p14), psetcols, pcols, pver)
    _copy2d_full(Ptr[float](p15), Ptr[float](p16), psetcols, pcols, pver)
    _copy2d_full(Ptr[float](p17), Ptr[float](p18), psetcols, pcols, pver)
    _copy2d_full(Ptr[float](p19), Ptr[float](p20), psetcols, pcols, pver)
    _copy2d_full(Ptr[float](p21), Ptr[float](p22), psetcols, pcols, pver)
    _copy2d_full(Ptr[float](p23), Ptr[float](p24), psetcols, pcols, pver)
    _copy2d_full(Ptr[float](p25), Ptr[float](p26), psetcols, pcols, pver)
    _copy2d_full(Ptr[float](p27), Ptr[float](p28), psetcols, pcols, pver)
    _copy2d_full(Ptr[float](p29), Ptr[float](p30), psetcols, pcols, pver)
    _copy2d_full(Ptr[float](p31), Ptr[float](p32), psetcols, pcols, pver)
    _copy2d_full(Ptr[float](p33), Ptr[float](p34), psetcols, pcols, pver)
    _copy2d_full(Ptr[float](p35), Ptr[float](p36), psetcols, pcols, pver)
    _copy2d_full(Ptr[float](p37), Ptr[float](p38), psetcols, pcols, pver)
    _copy2d_full(Ptr[float](p39), Ptr[float](p40), psetcols, pcols, pver)
    _copy2d_full(Ptr[float](p41), Ptr[float](p42), psetcols, pcols, pver)
    _copy2d_full(Ptr[float](p43), Ptr[float](p44), psetcols, pcols, pver)
    _copy2d_full(Ptr[float](p45), Ptr[float](p46), psetcols, pcols, pver)
    _copy2d_full(Ptr[float](p47), Ptr[float](p48), psetcols, pcols, pver)
    _copy2d_full(Ptr[float](p49), Ptr[float](p50), psetcols, pcols, pver)
    _copy2d_full(Ptr[float](p51), Ptr[float](p52), psetcols, pcols, pver)
    _copy2d_full(Ptr[float](p53), Ptr[float](p54), psetcols, pcols, pver)
    _copy2d_full(Ptr[float](p55), Ptr[float](p56), psetcols, pcols, pver)


@export
def micro_mg_cam_tail_state_grid_copy_codon(
    psetcols: int,
    pcols: int,
    pver: int,
    pcnst: int,
    ixnumliq: int,
    ixnumice: int,
    state_pdel_p: cobj,
    state_q_p: cobj,
    pdel_grid_p: cobj,
    nc_grid_p: cobj,
    ni_grid_p: cobj,
):
    state_pdel = Ptr[float](state_pdel_p)
    state_q = Ptr[float](state_q_p)
    pdel_grid = Ptr[float](pdel_grid_p)
    nc_grid = Ptr[float](nc_grid_p)
    ni_grid = Ptr[float](ni_grid_p)

    _copy2d_full(state_pdel, pdel_grid, psetcols, pcols, pver)
    for k in range(1, pver + 1):
        for i in range(1, pcols + 1):
            nc_grid[_idx2(i, k, pcols)] = state_q[_idx3(i, k, ixnumliq, psetcols, pver)]
            ni_grid[_idx2(i, k, pcols)] = state_q[_idx3(i, k, ixnumice, psetcols, pver)]


@export
def micro_mg_cam_diag_stage_dispatch_codon(
    ngrdcol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    micro_mg_version: int,
    qsmall: float,
    mincld: float,
    liq_rho: float,
    liq_eff_dim: float,
    liq_min_mean_mass: float,
    ice_eff_dim: float,
    ice_shape_coef: float,
    ice_lambda_lo: float,
    ice_lambda_hi: float,
    ice_min_mean_mass: float,
    rhosn: float,
    rhoi: float,
    rhow: float,
    rhows: float,
    mucon: float,
    dcon: float,
    deicon: float,
    minlwp: float,
    gravit: float,
    rhoh2o: float,
    rho_grid_p: cobj,
    icwmrst_grid_p: cobj,
    liqcldf_grid_p: cobj,
    nc_grid_p: cobj,
    qr_grid_p: cobj,
    nr_grid_p: cobj,
    qs_grid_p: cobj,
    ns_grid_p: cobj,
    qrout_grid_p: cobj,
    nrout_grid_p: cobj,
    qsout_grid_p: cobj,
    nsout_grid_p: cobj,
    ni_grid_p: cobj,
    icecldf_grid_p: cobj,
    icimrst_grid_p: cobj,
    ast_grid_p: cobj,
    mu_grid_p: cobj,
    lambdac_grid_p: cobj,
    rel_fn_grid_p: cobj,
    ncic_grid_p: cobj,
    rel_grid_p: cobj,
    drout2_grid_p: cobj,
    reff_rain_grid_p: cobj,
    des_grid_p: cobj,
    dsout2_grid_p: cobj,
    reff_snow_grid_p: cobj,
    rei_grid_p: cobj,
    niic_grid_p: cobj,
    dei_grid_p: cobj,
    mgreffrain_grid_p: cobj,
    mgreffsnow_grid_p: cobj,
    iclwpst_grid_p: cobj,
    cld_grid_p: cobj,
    cmeliq_grid_p: cobj,
    pdel_grid_p: cobj,
    prec_str_grid_p: cobj,
    acgcme_grid_p: cobj,
    acprecl_grid_p: cobj,
    acnum_grid_p: cobj,
    prao_grid_p: cobj,
    prco_grid_p: cobj,
    icwnc_grid_p: cobj,
    icinc_grid_p: cobj,
    nevapr_grid_p: cobj,
    evpsnow_st_grid_p: cobj,
    tgliqwp_grid_p: cobj,
    tgcmeliq_grid_p: cobj,
    pe_grid_p: cobj,
    tpr_grid_p: cobj,
    pefrac_grid_p: cobj,
    vprao_grid_p: cobj,
    vprco_grid_p: cobj,
    racau_grid_p: cobj,
    cnt_grid_p: cobj,
    cdnumc_grid_p: cobj,
    efcout_grid_p: cobj,
    efiout_grid_p: cobj,
    ncout_grid_p: cobj,
    niout_grid_p: cobj,
    freql_grid_p: cobj,
    freqi_grid_p: cobj,
    icwmrst_grid_out_p: cobj,
    icimrst_grid_out_p: cobj,
    fcti_grid_p: cobj,
    fctl_grid_p: cobj,
    ctrel_grid_p: cobj,
    ctrei_grid_p: cobj,
    ctnl_grid_p: cobj,
    ctni_grid_p: cobj,
    evprain_st_grid_p: cobj,
    qcreso_grid_p: cobj,
    melto_grid_p: cobj,
    mnuccco_grid_p: cobj,
    mnuccto_grid_p: cobj,
    bergo_grid_p: cobj,
    homoo_grid_p: cobj,
    msacwio_grid_p: cobj,
    psacwso_grid_p: cobj,
    bergso_grid_p: cobj,
    cmeiout_grid_p: cobj,
    qireso_grid_p: cobj,
    prcio_grid_p: cobj,
    praio_grid_p: cobj,
    budget_ftem_grid_p: cobj,
    copy_qrain: int,
    copy_qsnow: int,
    copy_nrain: int,
    copy_nsnow: int,
    qrout_grid_ptr_p: cobj,
    qsout_grid_ptr_p: cobj,
    nrout_grid_ptr_p: cobj,
    nsout_grid_ptr_p: cobj,
):
    micro_mg_cam_reff_calc_codon(
        ngrdcol,
        pcols,
        pver,
        top_lev,
        micro_mg_version,
        qsmall,
        mincld,
        liq_rho,
        liq_eff_dim,
        liq_min_mean_mass,
        ice_eff_dim,
        ice_shape_coef,
        ice_lambda_lo,
        ice_lambda_hi,
        ice_min_mean_mass,
        rhosn,
        rhoi,
        rhow,
        rhows,
        mucon,
        dcon,
        deicon,
        rho_grid_p,
        icwmrst_grid_p,
        liqcldf_grid_p,
        nc_grid_p,
        qr_grid_p,
        nr_grid_p,
        qs_grid_p,
        ns_grid_p,
        qrout_grid_p,
        nrout_grid_p,
        qsout_grid_p,
        nsout_grid_p,
        ni_grid_p,
        icecldf_grid_p,
        icimrst_grid_p,
        ast_grid_p,
        mu_grid_p,
        lambdac_grid_p,
        rel_fn_grid_p,
        ncic_grid_p,
        rel_grid_p,
        drout2_grid_p,
        reff_rain_grid_p,
        des_grid_p,
        dsout2_grid_p,
        reff_snow_grid_p,
        rei_grid_p,
        niic_grid_p,
        dei_grid_p,
        mgreffrain_grid_p,
        mgreffsnow_grid_p,
    )

    micro_mg_cam_grid_diag_codon(
        ngrdcol,
        pcols,
        pver,
        top_lev,
        minlwp,
        gravit,
        rhoh2o,
        iclwpst_grid_p,
        cld_grid_p,
        cmeliq_grid_p,
        pdel_grid_p,
        prec_str_grid_p,
        acgcme_grid_p,
        acprecl_grid_p,
        acnum_grid_p,
        prao_grid_p,
        prco_grid_p,
        nc_grid_p,
        liqcldf_grid_p,
        icwmrst_grid_p,
        rel_grid_p,
        icwnc_grid_p,
        icecldf_grid_p,
        icimrst_grid_p,
        rei_grid_p,
        icinc_grid_p,
        nevapr_grid_p,
        evpsnow_st_grid_p,
        tgliqwp_grid_p,
        tgcmeliq_grid_p,
        pe_grid_p,
        tpr_grid_p,
        pefrac_grid_p,
        vprao_grid_p,
        vprco_grid_p,
        racau_grid_p,
        cnt_grid_p,
        cdnumc_grid_p,
        efcout_grid_p,
        efiout_grid_p,
        ncout_grid_p,
        niout_grid_p,
        freql_grid_p,
        freqi_grid_p,
        icwmrst_grid_out_p,
        icimrst_grid_out_p,
        fcti_grid_p,
        fctl_grid_p,
        ctrel_grid_p,
        ctrei_grid_p,
        ctnl_grid_p,
        ctni_grid_p,
        evprain_st_grid_p,
    )

    budget_ftem_grid = Ptr[float](budget_ftem_grid_p)
    qcreso_grid = Ptr[float](qcreso_grid_p)
    melto_grid = Ptr[float](melto_grid_p)
    mnuccco_grid = Ptr[float](mnuccco_grid_p)
    mnuccto_grid = Ptr[float](mnuccto_grid_p)
    bergo_grid = Ptr[float](bergo_grid_p)
    homoo_grid = Ptr[float](homoo_grid_p)
    msacwio_grid = Ptr[float](msacwio_grid_p)
    prao_grid = Ptr[float](prao_grid_p)
    prco_grid = Ptr[float](prco_grid_p)
    psacwso_grid = Ptr[float](psacwso_grid_p)
    bergso_grid = Ptr[float](bergso_grid_p)
    cmeiout_grid = Ptr[float](cmeiout_grid_p)
    qireso_grid = Ptr[float](qireso_grid_p)
    prcio_grid = Ptr[float](prcio_grid_p)
    praio_grid = Ptr[float](praio_grid_p)

    for mode in range(1, 7):
        for k in range(1, pver + 1):
            for i in range(1, pcols + 1):
                budget_ftem_grid[_idx3(i, k, mode, pcols, pver)] = 0.0

    for k in range(top_lev, pver + 1):
        for i in range(1, ngrdcol + 1):
            idx2 = _idx2(i, k, pcols)

            budget_ftem_grid[_idx3(i, k, 1, pcols, pver)] = qcreso_grid[idx2]

            tmp = melto_grid[idx2] - mnuccco_grid[idx2]
            tmp = tmp - mnuccto_grid[idx2]
            tmp = tmp - bergo_grid[idx2]
            tmp = tmp - homoo_grid[idx2]
            tmp = tmp - msacwio_grid[idx2]
            budget_ftem_grid[_idx3(i, k, 2, pcols, pver)] = tmp

            tmp = -prao_grid[idx2]
            tmp = tmp - prco_grid[idx2]
            tmp = tmp - psacwso_grid[idx2]
            tmp = tmp - bergso_grid[idx2]
            budget_ftem_grid[_idx3(i, k, 3, pcols, pver)] = tmp

            budget_ftem_grid[_idx3(i, k, 4, pcols, pver)] = (
                cmeiout_grid[idx2] + qireso_grid[idx2]
            )

            tmp = -melto_grid[idx2] + mnuccco_grid[idx2]
            tmp = tmp + mnuccto_grid[idx2]
            tmp = tmp + bergo_grid[idx2]
            tmp = tmp + homoo_grid[idx2]
            tmp = tmp + msacwio_grid[idx2]
            budget_ftem_grid[_idx3(i, k, 5, pcols, pver)] = tmp

            budget_ftem_grid[_idx3(i, k, 6, pcols, pver)] = (
                -prcio_grid[idx2] - praio_grid[idx2]
            )

    if copy_qrain == 1 or copy_qsnow == 1 or copy_nrain == 1 or copy_nsnow == 1:
        _micro_mg_cam_pbuf_copy_fields(
            ngrdcol,
            pcols,
            pver,
            copy_qrain,
            copy_qsnow,
            copy_nrain,
            copy_nsnow,
            qrout_grid_p,
            qsout_grid_p,
            nrout_grid_p,
            nsout_grid_p,
            qrout_grid_ptr_p,
            qsout_grid_ptr_p,
            nrout_grid_ptr_p,
            nsout_grid_ptr_p,
        )

@export
def micro_mg_cam_diag_shell_codon(
    ngrdcol: int,
    pcols: int,
    pver: int,
    top_lev: int,
    micro_mg_version: int,
    qsmall: float,
    mincld: float,
    liq_rho: float,
    liq_eff_dim: float,
    liq_min_mean_mass: float,
    ice_eff_dim: float,
    ice_shape_coef: float,
    ice_lambda_lo: float,
    ice_lambda_hi: float,
    ice_min_mean_mass: float,
    rhosn: float,
    rhoi: float,
    rhow: float,
    rhows: float,
    mucon: float,
    dcon: float,
    deicon: float,
    minlwp: float,
    gravit: float,
    rhoh2o: float,
    rho_grid_p: cobj,
    icwmrst_grid_p: cobj,
    liqcldf_grid_p: cobj,
    nc_grid_p: cobj,
    qr_grid_p: cobj,
    nr_grid_p: cobj,
    qs_grid_p: cobj,
    ns_grid_p: cobj,
    qrout_grid_p: cobj,
    nrout_grid_p: cobj,
    qsout_grid_p: cobj,
    nsout_grid_p: cobj,
    ni_grid_p: cobj,
    icecldf_grid_p: cobj,
    icimrst_grid_p: cobj,
    ast_grid_p: cobj,
    mu_grid_p: cobj,
    lambdac_grid_p: cobj,
    rel_fn_grid_p: cobj,
    ncic_grid_p: cobj,
    rel_grid_p: cobj,
    drout2_grid_p: cobj,
    reff_rain_grid_p: cobj,
    des_grid_p: cobj,
    dsout2_grid_p: cobj,
    reff_snow_grid_p: cobj,
    rei_grid_p: cobj,
    niic_grid_p: cobj,
    dei_grid_p: cobj,
    mgreffrain_grid_p: cobj,
    mgreffsnow_grid_p: cobj,
    iclwpst_grid_p: cobj,
    cld_grid_p: cobj,
    cmeliq_grid_p: cobj,
    pdel_grid_p: cobj,
    prec_str_grid_p: cobj,
    acgcme_grid_p: cobj,
    acprecl_grid_p: cobj,
    acnum_grid_p: cobj,
    prao_grid_p: cobj,
    prco_grid_p: cobj,
    icwnc_grid_p: cobj,
    icinc_grid_p: cobj,
    nevapr_grid_p: cobj,
    evpsnow_st_grid_p: cobj,
    tgliqwp_grid_p: cobj,
    tgcmeliq_grid_p: cobj,
    pe_grid_p: cobj,
    tpr_grid_p: cobj,
    pefrac_grid_p: cobj,
    vprao_grid_p: cobj,
    vprco_grid_p: cobj,
    racau_grid_p: cobj,
    cnt_grid_p: cobj,
    cdnumc_grid_p: cobj,
    efcout_grid_p: cobj,
    efiout_grid_p: cobj,
    ncout_grid_p: cobj,
    niout_grid_p: cobj,
    freql_grid_p: cobj,
    freqi_grid_p: cobj,
    icwmrst_grid_out_p: cobj,
    icimrst_grid_out_p: cobj,
    fcti_grid_p: cobj,
    fctl_grid_p: cobj,
    ctrel_grid_p: cobj,
    ctrei_grid_p: cobj,
    ctnl_grid_p: cobj,
    ctni_grid_p: cobj,
    evprain_st_grid_p: cobj,
    qcreso_grid_p: cobj,
    melto_grid_p: cobj,
    mnuccco_grid_p: cobj,
    mnuccto_grid_p: cobj,
    bergo_grid_p: cobj,
    homoo_grid_p: cobj,
    msacwio_grid_p: cobj,
    psacwso_grid_p: cobj,
    bergso_grid_p: cobj,
    cmeiout_grid_p: cobj,
    qireso_grid_p: cobj,
    prcio_grid_p: cobj,
    praio_grid_p: cobj,
    budget_ftem_grid_p: cobj,
    copy_qrain: int,
    copy_qsnow: int,
    copy_nrain: int,
    copy_nsnow: int,
    qrout_grid_ptr_p: cobj,
    qsout_grid_ptr_p: cobj,
    nrout_grid_ptr_p: cobj,
    nsout_grid_ptr_p: cobj,
):
    micro_mg_cam_diag_stage_dispatch_codon(
        ngrdcol,
        pcols,
        pver,
        top_lev,
        micro_mg_version,
        qsmall,
        mincld,
        liq_rho,
        liq_eff_dim,
        liq_min_mean_mass,
        ice_eff_dim,
        ice_shape_coef,
        ice_lambda_lo,
        ice_lambda_hi,
        ice_min_mean_mass,
        rhosn,
        rhoi,
        rhow,
        rhows,
        mucon,
        dcon,
        deicon,
        minlwp,
        gravit,
        rhoh2o,
        rho_grid_p,
        icwmrst_grid_p,
        liqcldf_grid_p,
        nc_grid_p,
        qr_grid_p,
        nr_grid_p,
        qs_grid_p,
        ns_grid_p,
        qrout_grid_p,
        nrout_grid_p,
        qsout_grid_p,
        nsout_grid_p,
        ni_grid_p,
        icecldf_grid_p,
        icimrst_grid_p,
        ast_grid_p,
        mu_grid_p,
        lambdac_grid_p,
        rel_fn_grid_p,
        ncic_grid_p,
        rel_grid_p,
        drout2_grid_p,
        reff_rain_grid_p,
        des_grid_p,
        dsout2_grid_p,
        reff_snow_grid_p,
        rei_grid_p,
        niic_grid_p,
        dei_grid_p,
        mgreffrain_grid_p,
        mgreffsnow_grid_p,
        iclwpst_grid_p,
        cld_grid_p,
        cmeliq_grid_p,
        pdel_grid_p,
        prec_str_grid_p,
        acgcme_grid_p,
        acprecl_grid_p,
        acnum_grid_p,
        prao_grid_p,
        prco_grid_p,
        icwnc_grid_p,
        icinc_grid_p,
        nevapr_grid_p,
        evpsnow_st_grid_p,
        tgliqwp_grid_p,
        tgcmeliq_grid_p,
        pe_grid_p,
        tpr_grid_p,
        pefrac_grid_p,
        vprao_grid_p,
        vprco_grid_p,
        racau_grid_p,
        cnt_grid_p,
        cdnumc_grid_p,
        efcout_grid_p,
        efiout_grid_p,
        ncout_grid_p,
        niout_grid_p,
        freql_grid_p,
        freqi_grid_p,
        icwmrst_grid_out_p,
        icimrst_grid_out_p,
        fcti_grid_p,
        fctl_grid_p,
        ctrel_grid_p,
        ctrei_grid_p,
        ctnl_grid_p,
        ctni_grid_p,
        evprain_st_grid_p,
        qcreso_grid_p,
        melto_grid_p,
        mnuccco_grid_p,
        mnuccto_grid_p,
        bergo_grid_p,
        homoo_grid_p,
        msacwio_grid_p,
        psacwso_grid_p,
        bergso_grid_p,
        cmeiout_grid_p,
        qireso_grid_p,
        prcio_grid_p,
        praio_grid_p,
        budget_ftem_grid_p,
        copy_qrain,
        copy_qsnow,
        copy_nrain,
        copy_nsnow,
        qrout_grid_ptr_p,
        qsout_grid_ptr_p,
        nrout_grid_ptr_p,
        nsout_grid_ptr_p,
    )


@export
def micro_mg_cam_stage_dispatch_codon(
    stage: int,
    ncol: int,
    psetcols: int,
    pcols: int,
    pver: int,
    pverp: int,
    top_lev: int,
    micro_mg_version: int,
    rate1_cw2pr_st_idx: int,
    ixcldliq: int,
    ixcldice: int,
    ixnumliq: int,
    ixnumice: int,
    ixrain: int,
    ixsnow: int,
    pwtype: int,
    iwtvap: int,
    iwtliq: int,
    iwtice: int,
    iwtstrain: int,
    iwtstsnow: int,
    mincld: float,
    gravit: float,
    cpair: float,
    qsmall: float,
    rair: float,
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
    p35: cobj,
    p36: cobj,
    p37: cobj,
    p38: cobj,
    p39: cobj,
    p40: cobj,
    p41: cobj,
    p42: cobj,
    p43: cobj,
    p44: cobj,
    p45: cobj,
    p46: cobj,
    p47: cobj,
    p48: cobj,
    p49: cobj,
    p50: cobj,
    p51: cobj,
    p52: cobj,
    p53: cobj,
    p54: cobj,
    p55: cobj,
    p56: cobj,
    p57: cobj,
    p58: cobj,
    p59: cobj,
    p60: cobj,
    p61: cobj,
):
    if stage == 1:
        micro_mg_cam_premg_diag_codon(
            ncol,
            psetcols,
            pver,
            top_lev,
            ixcldliq,
            ixcldice,
            mincld,
            gravit,
            p1,
            p2,
            p3,
            p4,
            p5,
            p6,
        )
    elif stage == 2:
        micro_mg_cam_postmg_diag_codon(
            ncol,
            psetcols,
            pver,
            pverp,
            top_lev,
            micro_mg_version,
            rate1_cw2pr_st_idx,
            ixcldliq,
            ixcldice,
            ixnumliq,
            ixnumice,
            ixrain,
            ixsnow,
            cpair,
            gravit,
            mincld,
            qsmall,
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
            p35,
            p36,
            p37,
            p38,
            p39,
            p40,
            p41,
            p42,
            p43,
            p44,
            p45,
            p46,
            p47,
            p48,
            p49,
            p50,
            p51,
            p52,
            p53,
            p54,
            p55,
            p56,
            p57,
            p58,
            p59,
            p60,
            p61,
        )
    elif stage == 3:
        micro_mg_cam_wtrc_shell_codon(
            ncol,
            pcols,
            pver,
            top_lev,
            pwtype,
            iwtvap,
            iwtliq,
            iwtice,
            iwtstrain,
            iwtstsnow,
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
        )
    elif stage == 4:
        micro_mg_cam_rho_grid_codon(
            ncol,
            pcols,
            pver,
            psetcols,
            top_lev,
            rair,
            p1,
            p2,
            p3,
            p4,
        )
    elif stage == 5:
        micro_mg_cam_tail_grid_copy_codon(
            psetcols,
            pcols,
            pver,
            micro_mg_version,
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
            p35,
            p36,
            p37,
            p38,
            p39,
            p40,
            p41,
            p42,
            p43,
            p44,
            p45,
            p46,
            p47,
            p48,
            p49,
            p50,
            p51,
            p52,
            p53,
            p54,
            p55,
            p56,
        )

@export
def micro_mg_tend_codon(pcols: int, pver: int, ncol: int, top_lev: int, microp_uniform_c: int, do_cldice_c: int, do_clubb_sgs_c: int, use_hetfrz_classnuc_c: int, wv_sat_idx: int, deltatin: float, epsilo: float, omeps: float, ptrs_p: cobj, scalars_p: cobj):
    ptrs = Ptr[cobj](ptrs_p)
    scalars = Ptr[float](scalars_p)
    tn = Ptr[float](ptrs[0])
    qn = Ptr[float](ptrs[1])
    relvar = Ptr[float](ptrs[2])
    accre_enhan = Ptr[float](ptrs[3])
    qc = Ptr[float](ptrs[4])
    qi = Ptr[float](ptrs[5])
    nc = Ptr[float](ptrs[6])
    ni = Ptr[float](ptrs[7])
    p = Ptr[float](ptrs[8])
    pdel = Ptr[float](ptrs[9])
    cldn = Ptr[float](ptrs[10])
    icecldf = Ptr[float](ptrs[11])
    liqcldf = Ptr[float](ptrs[12])
    rate1ord_cw2pr_st = Ptr[float](ptrs[13])
    naai = Ptr[float](ptrs[14])
    npccnin = Ptr[float](ptrs[15])
    rndst = Ptr[float](ptrs[16])
    nacon = Ptr[float](ptrs[17])
    tlat = Ptr[float](ptrs[18])
    qvlat = Ptr[float](ptrs[19])
    qctend = Ptr[float](ptrs[20])
    qitend = Ptr[float](ptrs[21])
    nctend = Ptr[float](ptrs[22])
    nitend = Ptr[float](ptrs[23])
    effc = Ptr[float](ptrs[24])
    effc_fn = Ptr[float](ptrs[25])
    effi = Ptr[float](ptrs[26])
    prect = Ptr[float](ptrs[27])
    preci = Ptr[float](ptrs[28])
    nevapr = Ptr[float](ptrs[29])
    evapsnow = Ptr[float](ptrs[30])
    am_evp_st = Ptr[float](ptrs[31])
    prain = Ptr[float](ptrs[32])
    prodsnow = Ptr[float](ptrs[33])
    cmeout = Ptr[float](ptrs[34])
    deffi = Ptr[float](ptrs[35])
    pgamrad = Ptr[float](ptrs[36])
    lamcrad = Ptr[float](ptrs[37])
    qsout = Ptr[float](ptrs[38])
    dsout = Ptr[float](ptrs[39])
    rflx = Ptr[float](ptrs[40])
    sflx = Ptr[float](ptrs[41])
    qrout = Ptr[float](ptrs[42])
    reff_rain = Ptr[float](ptrs[43])
    reff_snow = Ptr[float](ptrs[44])
    qcsevap = Ptr[float](ptrs[45])
    qisevap = Ptr[float](ptrs[46])
    qvres = Ptr[float](ptrs[47])
    cmeiout = Ptr[float](ptrs[48])
    vtrmc = Ptr[float](ptrs[49])
    vtrmi = Ptr[float](ptrs[50])
    qcsedten = Ptr[float](ptrs[51])
    qisedten = Ptr[float](ptrs[52])
    prao = Ptr[float](ptrs[53])
    prco = Ptr[float](ptrs[54])
    mnuccco = Ptr[float](ptrs[55])
    mnuccto = Ptr[float](ptrs[56])
    msacwio = Ptr[float](ptrs[57])
    psacwso = Ptr[float](ptrs[58])
    bergso = Ptr[float](ptrs[59])
    bergo = Ptr[float](ptrs[60])
    melto = Ptr[float](ptrs[61])
    homoo = Ptr[float](ptrs[62])
    qcreso = Ptr[float](ptrs[63])
    prcio = Ptr[float](ptrs[64])
    praio = Ptr[float](ptrs[65])
    qireso = Ptr[float](ptrs[66])
    mnuccro = Ptr[float](ptrs[67])
    pracso = Ptr[float](ptrs[68])
    meltsdt = Ptr[float](ptrs[69])
    frzrdt = Ptr[float](ptrs[70])
    mnuccdo = Ptr[float](ptrs[71])
    nrout = Ptr[float](ptrs[72])
    nsout = Ptr[float](ptrs[73])
    refl = Ptr[float](ptrs[74])
    arefl = Ptr[float](ptrs[75])
    areflz = Ptr[float](ptrs[76])
    frefl = Ptr[float](ptrs[77])
    csrfl = Ptr[float](ptrs[78])
    acsrfl = Ptr[float](ptrs[79])
    fcsrfl = Ptr[float](ptrs[80])
    rercld = Ptr[float](ptrs[81])
    ncai = Ptr[float](ptrs[82])
    ncal = Ptr[float](ptrs[83])
    qrout2 = Ptr[float](ptrs[84])
    qsout2 = Ptr[float](ptrs[85])
    nrout2 = Ptr[float](ptrs[86])
    nsout2 = Ptr[float](ptrs[87])
    drout2 = Ptr[float](ptrs[88])
    dsout2 = Ptr[float](ptrs[89])
    freqs = Ptr[float](ptrs[90])
    freqr = Ptr[float](ptrs[91])
    nfice = Ptr[float](ptrs[92])
    prer_evap = Ptr[float](ptrs[93])
    nevapr2 = Ptr[float](ptrs[94])
    tnd_qsnow = Ptr[float](ptrs[95])
    tnd_nsnow = Ptr[float](ptrs[96])
    re_ice = Ptr[float](ptrs[97])
    frzimm = Ptr[float](ptrs[98])
    frzcnt = Ptr[float](ptrs[99])
    frzdep = Ptr[float](ptrs[100])
    preo = Ptr[float](ptrs[101])
    prdso = Ptr[float](ptrs[102])
    frzro = Ptr[float](ptrs[103])
    meltso = Ptr[float](ptrs[104])
    wtfc = Ptr[float](ptrs[105])
    wtfi = Ptr[float](ptrs[106])
    wtprelat = Ptr[float](ptrs[107])
    wtpostlat = Ptr[float](ptrs[108])
    t1 = Ptr[float](ptrs[109])
    q1 = Ptr[float](ptrs[110])
    qc1 = Ptr[float](ptrs[111])
    qi1 = Ptr[float](ptrs[112])
    nc1 = Ptr[float](ptrs[113])
    ni1 = Ptr[float](ptrs[114])
    tlat1 = Ptr[float](ptrs[115])
    qvlat1 = Ptr[float](ptrs[116])
    qctend1 = Ptr[float](ptrs[117])
    qitend1 = Ptr[float](ptrs[118])
    nctend1 = Ptr[float](ptrs[119])
    nitend1 = Ptr[float](ptrs[120])
    prect1 = Ptr[float](ptrs[121])
    preci1 = Ptr[float](ptrs[122])
    q = Ptr[float](ptrs[123])
    t = Ptr[float](ptrs[124])
    rho = Ptr[float](ptrs[125])
    dv = Ptr[float](ptrs[126])
    mu = Ptr[float](ptrs[127])
    sc = Ptr[float](ptrs[128])
    kap = Ptr[float](ptrs[129])
    rhof = Ptr[float](ptrs[130])
    cldmax = Ptr[float](ptrs[131])
    cldm = Ptr[float](ptrs[132])
    icldm = Ptr[float](ptrs[133])
    lcldm = Ptr[float](ptrs[134])
    icwc = Ptr[float](ptrs[135])
    calpha = Ptr[float](ptrs[136])
    cbeta = Ptr[float](ptrs[137])
    cbetah = Ptr[float](ptrs[138])
    cgamma = Ptr[float](ptrs[139])
    cgamah = Ptr[float](ptrs[140])
    rcgama = Ptr[float](ptrs[141])
    cmec1 = Ptr[float](ptrs[142])
    cmec2 = Ptr[float](ptrs[143])
    cmec3 = Ptr[float](ptrs[144])
    cmec4 = Ptr[float](ptrs[145])
    cme = Ptr[float](ptrs[146])
    cmei = Ptr[float](ptrs[147])
    cwml = Ptr[float](ptrs[148])
    cwmi = Ptr[float](ptrs[149])
    nnuccd = Ptr[float](ptrs[150])
    mnuccd = Ptr[float](ptrs[151])
    lcldn = Ptr[float](ptrs[152])
    lcldo = Ptr[float](ptrs[153])
    nctend_mixnuc = Ptr[float](ptrs[154])
    qcsinksum_rate1ord = Ptr[float](ptrs[155])
    qcsum_rate1ord = Ptr[float](ptrs[156])
    npccn = Ptr[float](ptrs[157])
    qcic = Ptr[float](ptrs[158])
    qiic = Ptr[float](ptrs[159])
    qniic = Ptr[float](ptrs[160])
    qric = Ptr[float](ptrs[161])
    ncic = Ptr[float](ptrs[162])
    niic = Ptr[float](ptrs[163])
    nsic = Ptr[float](ptrs[164])
    nric = Ptr[float](ptrs[165])
    lami = Ptr[float](ptrs[166])
    n0i = Ptr[float](ptrs[167])
    lamc = Ptr[float](ptrs[168])
    n0c = Ptr[float](ptrs[169])
    lams = Ptr[float](ptrs[170])
    n0s = Ptr[float](ptrs[171])
    lamr = Ptr[float](ptrs[172])
    n0r = Ptr[float](ptrs[173])
    cdist1 = Ptr[float](ptrs[174])
    arcld = Ptr[float](ptrs[175])
    pgam = Ptr[float](ptrs[176])
    mnuccc = Ptr[float](ptrs[177])
    nnuccc = Ptr[float](ptrs[178])
    mnucct = Ptr[float](ptrs[179])
    nnucct = Ptr[float](ptrs[180])
    msacwi = Ptr[float](ptrs[181])
    nsacwi = Ptr[float](ptrs[182])
    prc = Ptr[float](ptrs[183])
    nprc = Ptr[float](ptrs[184])
    nprc1 = Ptr[float](ptrs[185])
    nsagg = Ptr[float](ptrs[186])
    psacws = Ptr[float](ptrs[187])
    npsacws = Ptr[float](ptrs[188])
    uns = Ptr[float](ptrs[189])
    ums = Ptr[float](ptrs[190])
    unr = Ptr[float](ptrs[191])
    umr = Ptr[float](ptrs[192])
    pracs = Ptr[float](ptrs[193])
    npracs = Ptr[float](ptrs[194])
    mnuccr = Ptr[float](ptrs[195])
    nnuccr = Ptr[float](ptrs[196])
    pra = Ptr[float](ptrs[197])
    npra = Ptr[float](ptrs[198])
    nragg = Ptr[float](ptrs[199])
    prci = Ptr[float](ptrs[200])
    nprci = Ptr[float](ptrs[201])
    prai = Ptr[float](ptrs[202])
    nprai = Ptr[float](ptrs[203])
    pre = Ptr[float](ptrs[204])
    prds = Ptr[float](ptrs[205])
    dumc = Ptr[float](ptrs[206])
    dumnc = Ptr[float](ptrs[207])
    dumi = Ptr[float](ptrs[208])
    dumni = Ptr[float](ptrs[209])
    dums = Ptr[float](ptrs[210])
    dumns = Ptr[float](ptrs[211])
    dumr = Ptr[float](ptrs[212])
    dumnr = Ptr[float](ptrs[213])
    fr = Ptr[float](ptrs[214])
    fnr = Ptr[float](ptrs[215])
    fc = Ptr[float](ptrs[216])
    fnc = Ptr[float](ptrs[217])
    fi = Ptr[float](ptrs[218])
    fni = Ptr[float](ptrs[219])
    fs = Ptr[float](ptrs[220])
    fns = Ptr[float](ptrs[221])
    faloutr = Ptr[float](ptrs[222])
    faloutnr = Ptr[float](ptrs[223])
    faloutc = Ptr[float](ptrs[224])
    faloutnc = Ptr[float](ptrs[225])
    falouti = Ptr[float](ptrs[226])
    faloutni = Ptr[float](ptrs[227])
    falouts = Ptr[float](ptrs[228])
    faloutns = Ptr[float](ptrs[229])
    relhum = Ptr[float](ptrs[230])
    csigma = Ptr[float](ptrs[231])
    arn = Ptr[float](ptrs[232])
    asn = Ptr[float](ptrs[233])
    acn = Ptr[float](ptrs[234])
    ain = Ptr[float](ptrs[235])
    nsubi = Ptr[float](ptrs[236])
    nsubc = Ptr[float](ptrs[237])
    nsubs = Ptr[float](ptrs[238])
    nsubr = Ptr[float](ptrs[239])
    dz = Ptr[float](ptrs[240])
    rflx1 = Ptr[float](ptrs[241])
    sflx1 = Ptr[float](ptrs[242])
    tsp = Ptr[float](ptrs[243])
    qsp = Ptr[float](ptrs[244])
    qsphy = Ptr[float](ptrs[245])
    qs = Ptr[float](ptrs[246])
    es = Ptr[float](ptrs[247])
    esl = Ptr[float](ptrs[248])
    esi = Ptr[float](ptrs[249])
    qnitend = Ptr[float](ptrs[250])
    nstend = Ptr[float](ptrs[251])
    qrtend = Ptr[float](ptrs[252])
    nrtend = Ptr[float](ptrs[253])
    berg = Ptr[float](ptrs[254])
    bergs = Ptr[float](ptrs[255])
    drout = Ptr[float](ptrs[256])
    dum2i = Ptr[float](ptrs[257])
    dum2l = Ptr[float](ptrs[258])
    cldmw = Ptr[float](ptrs[259])
    qrtend_copy = Ptr[float](ptrs[260])
    qnitend_copy = Ptr[float](ptrs[261])
    rainrt = Ptr[float](ptrs[262])
    rainrt1 = Ptr[float](ptrs[263])
    mnudep = Ptr[float](ptrs[264])
    nnudep = Ptr[float](ptrs[265])
    ltrue = Ptr[int](ptrs[266])
    microp_uniform = microp_uniform_c != 0
    do_cldice = do_cldice_c != 0
    do_clubb_sgs = do_clubb_sgs_c != 0
    use_hetfrz_classnuc = use_hetfrz_classnuc_c != 0
    g = scalars[0]
    r = scalars[1]
    rv = scalars[2]
    cpp = scalars[3]
    rhow = scalars[4]
    tmelt = scalars[5]
    xxlv = scalars[6]
    xlf = scalars[7]
    xxls = scalars[8]
    rhosn = scalars[9]
    rhoi = scalars[10]
    ac = scalars[11]
    bc = scalars[12]
    as_ = scalars[13]
    bs = scalars[14]
    ai = scalars[15]
    bi = scalars[16]
    ar = scalars[17]
    br = scalars[18]
    ci = scalars[19]
    di = scalars[20]
    cs = scalars[21]
    ds = scalars[22]
    cr = scalars[23]
    dr = scalars[24]
    f1s = scalars[25]
    f2s = scalars[26]
    eii = scalars[27]
    ecr = scalars[28]
    f1r = scalars[29]
    f2r = scalars[30]
    dcs = scalars[31]
    qsmall = scalars[32]
    bimm = scalars[33]
    aimm = scalars[34]
    rhosu = scalars[35]
    mi0 = scalars[36]
    rin = scalars[37]
    pi = scalars[38]
    cons1 = scalars[39]
    cons4 = scalars[40]
    cons5 = scalars[41]
    cons6 = scalars[42]
    cons7 = scalars[43]
    cons8 = scalars[44]
    cons11 = scalars[45]
    cons13 = scalars[46]
    cons14 = scalars[47]
    cons16 = scalars[48]
    cons17 = scalars[49]
    cons22 = scalars[50]
    cons23 = scalars[51]
    cons24 = scalars[52]
    cons25 = scalars[53]
    cons27 = scalars[54]
    cons28 = scalars[55]
    lammini = scalars[56]
    lammaxi = scalars[57]
    lamminr = scalars[58]
    lammaxr = scalars[59]
    lammins = scalars[60]
    lammaxs = scalars[61]
    csmin = scalars[62]
    csmax = scalars[63]
    minrefl = scalars[64]
    mindbz = scalars[65]
    rhmini = scalars[66]
    micro_mg_berg_eff_factor = scalars[67]
    cdnl = 0.0
    cons2 = 0.0
    cons3 = 0.0
    cons9 = 0.0
    cons10 = 0.0
    cons12 = 0.0
    cons15 = 0.0
    cons18 = 0.0
    cons19 = 0.0
    cons20 = 0.0
    deltat = 0.0
    omsm = 0.0
    dto2 = 0.0
    mincld = 0.0
    qtmp = 0.0
    dum = 0.0
    qcld = 0.0
    arg = 0.0
    alpha = 0.0
    dum1 = 0.0
    dum2 = 0.0
    actmp = 0.0
    artmp = 0.0
    lammax = 0.0
    lammin = 0.0
    nacnt = 0.0
    dc0 = 0.0
    ds0 = 0.0
    eci = 0.0
    uni = 0.0
    umi = 0.0
    unc = 0.0
    umc = 0.0
    qvs = 0.0
    qvi = 0.0
    dqsdt = 0.0
    dqsidt = 0.0
    ab = 0.0
    qclr = 0.0
    abi = 0.0
    epss = 0.0
    epsr = 0.0
    qce = 0.0
    qie = 0.0
    nce = 0.0
    nie = 0.0
    ratio = 0.0
    faltndr = 0.0
    faltndnr = 0.0
    faltndc = 0.0
    faltndnc = 0.0
    faltndi = 0.0
    faltndni = 0.0
    faltnds = 0.0
    faltndns = 0.0
    faltndqie = 0.0
    faltndqce = 0.0
    rgvm = 0.0
    mtime = 0.0
    qrtot = 0.0
    nrtot = 0.0
    qstot = 0.0
    nstot = 0.0
    dumnnuc = 0.0
    ninew = 0.0
    qinew = 0.0
    qvl = 0.0
    epsi = 0.0
    prd = 0.0
    bergtsf = 0.0
    rhin = 0.0
    dumfice = 0.0
    ncmax = 0.0
    nimax = 0.0
    qcvar = 0.0
    tcnt = 0.0
    viscosity = 0.0
    mfp = 0.0
    slip1 = 0.0
    slip2 = 0.0
    slip3 = 0.0
    slip4 = 0.0
    ndfaer1 = 0.0
    ndfaer2 = 0.0
    ndfaer3 = 0.0
    ndfaer4 = 0.0
    nslip1 = 0.0
    nslip2 = 0.0
    nslip3 = 0.0
    nslip4 = 0.0
    bbi = 0.0
    cci = 0.0
    ak = 0.0
    iciwc = 0.0
    rvi = 0.0
    tk = 0.0
    deles = 0.0
    aprpr = 0.0
    bprpr = 0.0
    cice = 0.0
    qi0 = 0.0
    crate = 0.0
    qidep = 0.0
    ni_secp = 0.0
    esn = 0.0
    qsn = 0.0
    ttmp = 0.0
    tmp = 0.0
    dmc = 0.0
    ssmc = 0.0
    dstrn = 0.0
    con1 = 0.0
    r3lx = 0.0
    mi0l = 0.0
    frztmp = 0.0
    for _s2 in range(1, pver + 1):
        for _s1 in range(1, ncol + 1):
            ncai[_idx2(_s1, _s2, pcols)] = 0.0
    for _s4 in range(1, pver + 1):
        for _s3 in range(1, ncol + 1):
            ncal[_idx2(_s3, _s4, pcols)] = 0.0
    for _s6 in range(1, pver + 1):
        for _s5 in range(1, ncol + 1):
            rercld[_idx2(_s5, _s6, pcols)] = 0.0
    for _s8 in range(1, pver + 1):
        for _s7 in range(1, ncol + 1):
            arcld[_idx2(_s7, _s8, pcols)] = 0.0
    for _s10 in range(1, pver + 1):
        for _s9 in range(1, ncol + 1):
            pgamrad[_idx2(_s9, _s10, pcols)] = 0.0
    for _s12 in range(1, pver + 1):
        for _s11 in range(1, ncol + 1):
            lamcrad[_idx2(_s11, _s12, pcols)] = 0.0
    for _s14 in range(1, pver + 1):
        for _s13 in range(1, ncol + 1):
            deffi[_idx2(_s13, _s14, pcols)] = 0.0
    for _s16 in range(1, pver + 1):
        for _s15 in range(1, ncol + 1):
            qcsevap[_idx2(_s15, _s16, pcols)] = 0.0
    for _s18 in range(1, pver + 1):
        for _s17 in range(1, ncol + 1):
            qisevap[_idx2(_s17, _s18, pcols)] = 0.0
    for _s20 in range(1, pver + 1):
        for _s19 in range(1, ncol + 1):
            qvres[_idx2(_s19, _s20, pcols)] = 0.0
    for _s22 in range(1, pver + 1):
        for _s21 in range(1, ncol + 1):
            cmeiout[_idx2(_s21, _s22, pcols)] = 0.0
    for _s24 in range(1, pver + 1):
        for _s23 in range(1, ncol + 1):
            vtrmc[_idx2(_s23, _s24, pcols)] = 0.0
    for _s26 in range(1, pver + 1):
        for _s25 in range(1, ncol + 1):
            vtrmi[_idx2(_s25, _s26, pcols)] = 0.0
    for _s28 in range(1, pver + 1):
        for _s27 in range(1, ncol + 1):
            qcsedten[_idx2(_s27, _s28, pcols)] = 0.0
    for _s30 in range(1, pver + 1):
        for _s29 in range(1, ncol + 1):
            qisedten[_idx2(_s29, _s30, pcols)] = 0.0
    for _s32 in range(1, pver + 1):
        for _s31 in range(1, ncol + 1):
            prao[_idx2(_s31, _s32, pcols)] = 0.0
    for _s34 in range(1, pver + 1):
        for _s33 in range(1, ncol + 1):
            prco[_idx2(_s33, _s34, pcols)] = 0.0
    for _s36 in range(1, pver + 1):
        for _s35 in range(1, ncol + 1):
            mnuccco[_idx2(_s35, _s36, pcols)] = 0.0
    for _s38 in range(1, pver + 1):
        for _s37 in range(1, ncol + 1):
            mnuccto[_idx2(_s37, _s38, pcols)] = 0.0
    for _s40 in range(1, pver + 1):
        for _s39 in range(1, ncol + 1):
            msacwio[_idx2(_s39, _s40, pcols)] = 0.0
    for _s42 in range(1, pver + 1):
        for _s41 in range(1, ncol + 1):
            psacwso[_idx2(_s41, _s42, pcols)] = 0.0
    for _s44 in range(1, pver + 1):
        for _s43 in range(1, ncol + 1):
            bergso[_idx2(_s43, _s44, pcols)] = 0.0
    for _s46 in range(1, pver + 1):
        for _s45 in range(1, ncol + 1):
            bergo[_idx2(_s45, _s46, pcols)] = 0.0
    for _s48 in range(1, pver + 1):
        for _s47 in range(1, ncol + 1):
            melto[_idx2(_s47, _s48, pcols)] = 0.0
    for _s50 in range(1, pver + 1):
        for _s49 in range(1, ncol + 1):
            homoo[_idx2(_s49, _s50, pcols)] = 0.0
    for _s52 in range(1, pver + 1):
        for _s51 in range(1, ncol + 1):
            qcreso[_idx2(_s51, _s52, pcols)] = 0.0
    for _s54 in range(1, pver + 1):
        for _s53 in range(1, ncol + 1):
            prcio[_idx2(_s53, _s54, pcols)] = 0.0
    for _s56 in range(1, pver + 1):
        for _s55 in range(1, ncol + 1):
            praio[_idx2(_s55, _s56, pcols)] = 0.0
    for _s58 in range(1, pver + 1):
        for _s57 in range(1, ncol + 1):
            qireso[_idx2(_s57, _s58, pcols)] = 0.0
    for _s60 in range(1, pver + 1):
        for _s59 in range(1, ncol + 1):
            mnuccro[_idx2(_s59, _s60, pcols)] = 0.0
    for _s62 in range(1, pver + 1):
        for _s61 in range(1, ncol + 1):
            pracso[_idx2(_s61, _s62, pcols)] = 0.0
    for _s64 in range(1, pver + 1):
        for _s63 in range(1, ncol + 1):
            meltsdt[_idx2(_s63, _s64, pcols)] = 0.0
    for _s66 in range(1, pver + 1):
        for _s65 in range(1, ncol + 1):
            frzrdt[_idx2(_s65, _s66, pcols)] = 0.0
    for _s68 in range(1, pver + 1):
        for _s67 in range(1, ncol + 1):
            mnuccdo[_idx2(_s67, _s68, pcols)] = 0.0
    for _s70 in range(1, pver + 1 + 1):
        for _s69 in range(1, pcols + 1):
            rflx[_idx2(_s69, _s70, pcols)] = 0.0
    for _s72 in range(1, pver + 1 + 1):
        for _s71 in range(1, pcols + 1):
            sflx[_idx2(_s71, _s72, pcols)] = 0.0
    for _s74 in range(1, pver + 1):
        for _s73 in range(1, pcols + 1):
            effc[_idx2(_s73, _s74, pcols)] = 0.0
    for _s76 in range(1, pver + 1):
        for _s75 in range(1, pcols + 1):
            effc_fn[_idx2(_s75, _s76, pcols)] = 0.0
    for _s78 in range(1, pver + 1):
        for _s77 in range(1, pcols + 1):
            effi[_idx2(_s77, _s78, pcols)] = 0.0
    for _s80 in range(1, pver + 1):
        for _s79 in range(1, ncol + 1):
            preo[_idx2(_s79, _s80, pcols)] = 0.0
    for _s82 in range(1, pver + 1):
        for _s81 in range(1, ncol + 1):
            prdso[_idx2(_s81, _s82, pcols)] = 0.0
    for _s84 in range(1, pver + 1):
        for _s83 in range(1, ncol + 1):
            frzro[_idx2(_s83, _s84, pcols)] = 0.0
    for _s86 in range(1, pver + 1):
        for _s85 in range(1, ncol + 1):
            meltso[_idx2(_s85, _s86, pcols)] = 0.0
    for _s88 in range(1, pver + 1):
        for _s87 in range(1, ncol + 1):
            wtfc[_idx2(_s87, _s88, pcols)] = 0.0
    for _s90 in range(1, pver + 1):
        for _s89 in range(1, ncol + 1):
            wtfi[_idx2(_s89, _s90, pcols)] = 0.0
    for _s92 in range(1, pver + 1):
        for _s91 in range(1, ncol + 1):
            wtprelat[_idx2(_s91, _s92, pcols)] = 0.0
    for _s94 in range(1, pver + 1):
        for _s93 in range(1, ncol + 1):
            wtpostlat[_idx2(_s93, _s94, pcols)] = 0.0
    deltat = deltatin
    omsm = 0.99999
    dto2 = 0.5 * deltat
    mincld = 0.0001
    for _s96 in range(1, pver + 1):
        for _s95 in range(1, ncol + 1):
            q[_idx2(_s95, _s96, pcols)] = qn[_idx2(_s95, _s96, pcols)]
    for _s98 in range(1, pver + 1):
        for _s97 in range(1, ncol + 1):
            t[_idx2(_s97, _s98, pcols)] = tn[_idx2(_s97, _s98, pcols)]
    for _s100 in range(1, top_lev - 1 + 1):
        for _s99 in range(1, ncol + 1):
            qc[_idx2(_s99, _s100, pcols)] = 0.0
    for _s102 in range(1, top_lev - 1 + 1):
        for _s101 in range(1, ncol + 1):
            qi[_idx2(_s101, _s102, pcols)] = 0.0
    for _s104 in range(1, top_lev - 1 + 1):
        for _s103 in range(1, ncol + 1):
            nc[_idx2(_s103, _s104, pcols)] = 0.0
    for _s106 in range(1, top_lev - 1 + 1):
        for _s105 in range(1, ncol + 1):
            ni[_idx2(_s105, _s106, pcols)] = 0.0
    for _s108 in range(1, pver + 1):
        for _s107 in range(1, ncol + 1):
            t1[_idx2(_s107, _s108, pcols)] = t[_idx2(_s107, _s108, pcols)]
    for _s110 in range(1, pver + 1):
        for _s109 in range(1, ncol + 1):
            q1[_idx2(_s109, _s110, pcols)] = q[_idx2(_s109, _s110, pcols)]
    for _s112 in range(1, pver + 1):
        for _s111 in range(1, ncol + 1):
            qc1[_idx2(_s111, _s112, pcols)] = qc[_idx2(_s111, _s112, pcols)]
    for _s114 in range(1, pver + 1):
        for _s113 in range(1, ncol + 1):
            qi1[_idx2(_s113, _s114, pcols)] = qi[_idx2(_s113, _s114, pcols)]
    for _s116 in range(1, pver + 1):
        for _s115 in range(1, ncol + 1):
            nc1[_idx2(_s115, _s116, pcols)] = nc[_idx2(_s115, _s116, pcols)]
    for _s118 in range(1, pver + 1):
        for _s117 in range(1, ncol + 1):
            ni1[_idx2(_s117, _s118, pcols)] = ni[_idx2(_s117, _s118, pcols)]
    for _s120 in range(1, pver + 1):
        for _s119 in range(1, ncol + 1):
            tlat1[_idx2(_s119, _s120, pcols)] = 0.0
    for _s122 in range(1, pver + 1):
        for _s121 in range(1, ncol + 1):
            qvlat1[_idx2(_s121, _s122, pcols)] = 0.0
    for _s124 in range(1, pver + 1):
        for _s123 in range(1, ncol + 1):
            qctend1[_idx2(_s123, _s124, pcols)] = 0.0
    for _s126 in range(1, pver + 1):
        for _s125 in range(1, ncol + 1):
            qitend1[_idx2(_s125, _s126, pcols)] = 0.0
    for _s128 in range(1, pver + 1):
        for _s127 in range(1, ncol + 1):
            nctend1[_idx2(_s127, _s128, pcols)] = 0.0
    for _s130 in range(1, pver + 1):
        for _s129 in range(1, ncol + 1):
            nitend1[_idx2(_s129, _s130, pcols)] = 0.0
    for _s132 in range(1, pver + 1):
        for _s131 in range(1, ncol + 1):
            qrout[_idx2(_s131, _s132, pcols)] = 0.0
    for _s134 in range(1, pver + 1):
        for _s133 in range(1, ncol + 1):
            qsout[_idx2(_s133, _s134, pcols)] = 0.0
    for _s136 in range(1, pver + 1):
        for _s135 in range(1, ncol + 1):
            nrout[_idx2(_s135, _s136, pcols)] = 0.0
    for _s138 in range(1, pver + 1):
        for _s137 in range(1, ncol + 1):
            nsout[_idx2(_s137, _s138, pcols)] = 0.0
    for _s140 in range(1, pver + 1):
        for _s139 in range(1, ncol + 1):
            dsout[_idx2(_s139, _s140, pcols)] = 0.0
    for _s142 in range(1, pver + 1):
        for _s141 in range(1, ncol + 1):
            drout[_idx2(_s141, _s142, pcols)] = 0.0
    for _s144 in range(1, pver + 1):
        for _s143 in range(1, ncol + 1):
            reff_rain[_idx2(_s143, _s144, pcols)] = 0.0
    for _s146 in range(1, pver + 1):
        for _s145 in range(1, ncol + 1):
            reff_snow[_idx2(_s145, _s146, pcols)] = 0.0
    for _s148 in range(1, pver + 1):
        for _s147 in range(1, ncol + 1):
            nevapr[_idx2(_s147, _s148, pcols)] = 0.0
    for _s150 in range(1, pver + 1):
        for _s149 in range(1, ncol + 1):
            nevapr2[_idx2(_s149, _s150, pcols)] = 0.0
    for _s152 in range(1, pver + 1):
        for _s151 in range(1, ncol + 1):
            evapsnow[_idx2(_s151, _s152, pcols)] = 0.0
    for _s154 in range(1, pver + 1):
        for _s153 in range(1, ncol + 1):
            prain[_idx2(_s153, _s154, pcols)] = 0.0
    for _s156 in range(1, pver + 1):
        for _s155 in range(1, ncol + 1):
            prodsnow[_idx2(_s155, _s156, pcols)] = 0.0
    for _s158 in range(1, pver + 1):
        for _s157 in range(1, ncol + 1):
            cmeout[_idx2(_s157, _s158, pcols)] = 0.0
    for _s160 in range(1, pver + 1):
        for _s159 in range(1, ncol + 1):
            am_evp_st[_idx2(_s159, _s160, pcols)] = 0.0
    for _s162 in range(1, pver + 1):
        for _s161 in range(1, ncol + 1):
            rainrt1[_idx2(_s161, _s162, pcols)] = 0.0
    for _s164 in range(1, pver + 1):
        for _s163 in range(1, ncol + 1):
            cldmax[_idx2(_s163, _s164, pcols)] = mincld
    for _s166 in range(1, pver + 1):
        for _s165 in range(1, ncol + 1):
            dum2l[_idx2(_s165, _s166, pcols)] = 0.0
    for _s168 in range(1, pver + 1):
        for _s167 in range(1, ncol + 1):
            dum2i[_idx2(_s167, _s168, pcols)] = 0.0
    for _s169 in range(1, ncol + 1):
        prect1[_s169 - 1] = 0.0
    for _s170 in range(1, ncol + 1):
        preci1[_s170 - 1] = 0.0
    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            _idx = _idx2(i, k, pcols)
            rho[_idx] = p[_idx] / (r * t[_idx])
            dv[_idx] = 8.794E-5 * pow_r8(t[_idx], 1.81) / p[_idx]
            mu[_idx] = 1.496E-6 * pow_r8(t[_idx], 1.5) / (t[_idx] + 120.0)
            sc[_idx] = mu[_idx] / (rho[_idx] * dv[_idx])
            kap[_idx] = 1.414e3 * 1.496e-6 * pow_r8(t[_idx], 1.5) / (t[_idx] + 120.0)
            rhof[_idx] = pow_r8(rhosu / rho[_idx], 0.54)
            arn[_idx] = ar * rhof[_idx]
            asn[_idx] = as_ * rhof[_idx]
            acn[_idx] = ac * rhof[_idx]
            ain[_idx] = ai * rhof[_idx]
            dz[_idx] = pdel[_idx] / (rho[_idx] * g)
    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            es[i - 1] = wv_sat_svp_water_codon(t[_idx2(i, k, pcols)], wv_sat_idx)
            qs[i - 1] = wv_sat_svp_to_qsat_codon(es[i - 1], p[_idx2(i, k, pcols)], epsilo, omeps)
            if qs[i - 1] < 0.0:
                qs[i - 1] = 1.0
                es[i - 1] = p[_idx2(i, k, pcols)]
            esl[_idx2(i, k, pcols)] = wv_sat_svp_water_codon(t[_idx2(i, k, pcols)], wv_sat_idx)
            esi[_idx2(i, k, pcols)] = wv_sat_svp_ice_codon(t[_idx2(i, k, pcols)], wv_sat_idx)
            if t[_idx2(i, k, pcols)] > tmelt:
                esi[_idx2(i, k, pcols)] = esl[_idx2(i, k, pcols)]
            relhum[_idx2(i, k, pcols)] = q[_idx2(i, k, pcols)] / qs[i - 1]
            cldm[_idx2(i, k, pcols)] = max(cldn[_idx2(i, k, pcols)], mincld)
            cldmw[_idx2(i, k, pcols)] = max(cldn[_idx2(i, k, pcols)], mincld)
            icldm[_idx2(i, k, pcols)] = max(icecldf[_idx2(i, k, pcols)], mincld)
            lcldm[_idx2(i, k, pcols)] = max(liqcldf[_idx2(i, k, pcols)], mincld)
            if microp_uniform:
                cldm[_idx2(i, k, pcols)] = mincld
                cldmw[_idx2(i, k, pcols)] = mincld
                icldm[_idx2(i, k, pcols)] = mincld
                lcldm[_idx2(i, k, pcols)] = mincld
                if qc[_idx2(i, k, pcols)] >= qsmall:
                    lcldm[_idx2(i, k, pcols)] = 1.0
                    cldm[_idx2(i, k, pcols)] = 1.0
                    cldmw[_idx2(i, k, pcols)] = 1.0
                if qi[_idx2(i, k, pcols)] >= qsmall:
                    cldm[_idx2(i, k, pcols)] = 1.0
                    icldm[_idx2(i, k, pcols)] = 1.0
            nfice[_idx2(i, k, pcols)] = 0.0
            dumfice = qc[_idx2(i, k, pcols)] + qi[_idx2(i, k, pcols)]
            if dumfice > qsmall and qi[_idx2(i, k, pcols)] > qsmall:
                nfice[_idx2(i, k, pcols)] = qi[_idx2(i, k, pcols)] / dumfice
            if do_cldice and t[_idx2(i, k, pcols)] < tmelt - 5.0:
                dum2 = naai[_idx2(i, k, pcols)]
                dumnnuc = (dum2 - ni[_idx2(i, k, pcols)] / icldm[_idx2(i, k, pcols)]) / deltat * icldm[_idx2(i, k, pcols)]
                dumnnuc = max(dumnnuc, 0.0)
                ninew = ni[_idx2(i, k, pcols)] + dumnnuc * deltat
                qinew = qi[_idx2(i, k, pcols)] + dumnnuc * deltat * mi0
            else:
                ninew = ni[_idx2(i, k, pcols)]
                qinew = qi[_idx2(i, k, pcols)]
            cme[_idx2(i, k, pcols)] = 0.0
            cmei[_idx2(i, k, pcols)] = 0.0
            berg[_idx2(i, k, pcols)] = 0.0
            prd = 0.0
            if icldm[_idx2(i, k, pcols)] > 0.0:
                qiic[_idx2(i, k, pcols)] = qinew / icldm[_idx2(i, k, pcols)]
                niic[_idx2(i, k, pcols)] = ninew / icldm[_idx2(i, k, pcols)]
            else:
                qiic[_idx2(i, k, pcols)] = 0.0
                niic[_idx2(i, k, pcols)] = 0.0
            if do_cldice and t[_idx2(i, k, pcols)] < 273.15:
                if qi[_idx2(i, k, pcols)] > qsmall:
                    bergtsf = 0.0
                    qvi = wv_sat_svp_to_qsat_codon(esi[_idx2(i, k, pcols)], p[_idx2(i, k, pcols)], epsilo, omeps)
                    qvl = wv_sat_svp_to_qsat_codon(esl[_idx2(i, k, pcols)], p[_idx2(i, k, pcols)], epsilo, omeps)
                    dqsidt = xxls * qvi / (rv * pow_i(t[_idx2(i, k, pcols)], 2))
                    abi = 1.0 + dqsidt * xxls / cpp
                    if qiic[_idx2(i, k, pcols)] >= qsmall:
                        lami[k - 1] = pow_r8(cons1 * ci * niic[_idx2(i, k, pcols)] / qiic[_idx2(i, k, pcols)], 1.0 / di)
                        n0i[k - 1] = niic[_idx2(i, k, pcols)] * lami[k - 1]
                        if lami[k - 1] < lammini:
                            lami[k - 1] = lammini
                            n0i[k - 1] = pow_r8(lami[k - 1], di + 1.0) * qiic[_idx2(i, k, pcols)] / (ci * cons1)
                        elif lami[k - 1] > lammaxi:
                            lami[k - 1] = lammaxi
                            n0i[k - 1] = pow_r8(lami[k - 1], di + 1.0) * qiic[_idx2(i, k, pcols)] / (ci * cons1)
                        epsi = 2.0 * pi * n0i[k - 1] * rho[_idx2(i, k, pcols)] * dv[_idx2(i, k, pcols)] / (lami[k - 1] * lami[k - 1])
                        if qc[_idx2(i, k, pcols)] > qsmall:
                            prd = epsi * (qvl - qvi) / abi
                        else:
                            prd = 0.0
                        prd = prd * min(icldm[_idx2(i, k, pcols)], lcldm[_idx2(i, k, pcols)])
                        berg[_idx2(i, k, pcols)] = max(0.0, prd)
                    if berg[_idx2(i, k, pcols)] > 0.0:
                        bergtsf = max(0.0, qc[_idx2(i, k, pcols)] / berg[_idx2(i, k, pcols)] / deltat)
                        if bergtsf < 1.0:
                            berg[_idx2(i, k, pcols)] = max(0.0, qc[_idx2(i, k, pcols)] / deltat)
                    if bergtsf < 1.0 or icldm[_idx2(i, k, pcols)] > lcldm[_idx2(i, k, pcols)]:
                        if qiic[_idx2(i, k, pcols)] >= qsmall:
                            if qc[_idx2(i, k, pcols)] >= qsmall:
                                rhin = (1.0 + relhum[_idx2(i, k, pcols)]) / 2.0
                                if rhin * esl[_idx2(i, k, pcols)] / esi[_idx2(i, k, pcols)] > 1.0:
                                    prd = epsi * (rhin * qvl - qvi) / abi
                                    prd = prd * min(icldm[_idx2(i, k, pcols)], lcldm[_idx2(i, k, pcols)])
                                    cmei[_idx2(i, k, pcols)] = cmei[_idx2(i, k, pcols)] + prd * (1.0 - bergtsf)
                            if qc[_idx2(i, k, pcols)] < qsmall or icldm[_idx2(i, k, pcols)] > lcldm[_idx2(i, k, pcols)]:
                                if qc[_idx2(i, k, pcols)] < qsmall:
                                    dum = 0.0
                                else:
                                    dum = lcldm[_idx2(i, k, pcols)]
                                rhin = relhum[_idx2(i, k, pcols)]
                                if rhin * esl[_idx2(i, k, pcols)] / esi[_idx2(i, k, pcols)] > 1.0:
                                    prd = epsi * (rhin * qvl - qvi) / abi
                                    prd = prd * max(icldm[_idx2(i, k, pcols)] - dum, 0.0)
                                    cmei[_idx2(i, k, pcols)] = cmei[_idx2(i, k, pcols)] + prd
                    if cmei[_idx2(i, k, pcols)] > 0.0 and relhum[_idx2(i, k, pcols)] * esl[_idx2(i, k, pcols)] / esi[_idx2(i, k, pcols)] > 1.0:
                        cmei[_idx2(i, k, pcols)] = min(cmei[_idx2(i, k, pcols)], (q[_idx2(i, k, pcols)] - qs[i - 1] * esi[_idx2(i, k, pcols)] / esl[_idx2(i, k, pcols)]) / abi / deltat)
            if -berg[_idx2(i, k, pcols)] < -qc[_idx2(i, k, pcols)] / deltat:
                berg[_idx2(i, k, pcols)] = max(qc[_idx2(i, k, pcols)] / deltat, 0.0)
            if do_cldice and (relhum[_idx2(i, k, pcols)] * esl[_idx2(i, k, pcols)] / esi[_idx2(i, k, pcols)] < 1.0 and qiic[_idx2(i, k, pcols)] >= qsmall):
                qvi = wv_sat_svp_to_qsat_codon(esi[_idx2(i, k, pcols)], p[_idx2(i, k, pcols)], epsilo, omeps)
                qvl = wv_sat_svp_to_qsat_codon(esl[_idx2(i, k, pcols)], p[_idx2(i, k, pcols)], epsilo, omeps)
                dqsidt = xxls * qvi / (rv * pow_i(t[_idx2(i, k, pcols)], 2))
                abi = 1.0 + dqsidt * xxls / cpp
                lami[k - 1] = pow_r8(cons1 * ci * niic[_idx2(i, k, pcols)] / qiic[_idx2(i, k, pcols)], 1.0 / di)
                n0i[k - 1] = niic[_idx2(i, k, pcols)] * lami[k - 1]
                if lami[k - 1] < lammini:
                    lami[k - 1] = lammini
                    n0i[k - 1] = pow_r8(lami[k - 1], di + 1.0) * qiic[_idx2(i, k, pcols)] / (ci * cons1)
                elif lami[k - 1] > lammaxi:
                    lami[k - 1] = lammaxi
                    n0i[k - 1] = pow_r8(lami[k - 1], di + 1.0) * qiic[_idx2(i, k, pcols)] / (ci * cons1)
                epsi = 2.0 * pi * n0i[k - 1] * rho[_idx2(i, k, pcols)] * dv[_idx2(i, k, pcols)] / (lami[k - 1] * lami[k - 1])
                prd = epsi * (relhum[_idx2(i, k, pcols)] * qvl - qvi) / abi * icldm[_idx2(i, k, pcols)]
                cmei[_idx2(i, k, pcols)] = min(prd, 0.0)
            if cmei[_idx2(i, k, pcols)] < -qi[_idx2(i, k, pcols)] / deltat:
                cmei[_idx2(i, k, pcols)] = -qi[_idx2(i, k, pcols)] / deltat
            if cmei[_idx2(i, k, pcols)] < 0.0 and relhum[_idx2(i, k, pcols)] * esl[_idx2(i, k, pcols)] / esi[_idx2(i, k, pcols)] < 1.0:
                cmei[_idx2(i, k, pcols)] = min(0.0, max(cmei[_idx2(i, k, pcols)], (q[_idx2(i, k, pcols)] - qs[i - 1] * esi[_idx2(i, k, pcols)] / esl[_idx2(i, k, pcols)]) / abi / deltat))
            cmei[_idx2(i, k, pcols)] = cmei[_idx2(i, k, pcols)] * omsm
            if do_cldice and t[_idx2(i, k, pcols)] < tmelt - 5.0:
                dum2i[_idx2(i, k, pcols)] = naai[_idx2(i, k, pcols)]
            else:
                dum2i[_idx2(i, k, pcols)] = 0.0
    for i in range(1, ncol + 1):
        rflx1[_idx2(i, 1, pcols)] = 0.0
        sflx1[_idx2(i, 1, pcols)] = 0.0
        for k in range(top_lev, pver + 1):
            rflx1[_idx2(i, k + 1, pcols)] = 0.0
            sflx1[_idx2(i, k + 1, pcols)] = 0.0
    for i in range(1, ncol + 1):
        rflx[_idx2(i, 1, pcols)] = 0.0
        sflx[_idx2(i, 1, pcols)] = 0.0
        for k in range(top_lev, pver + 1):
            rflx[_idx2(i, k + 1, pcols)] = 0.0
            sflx[_idx2(i, k + 1, pcols)] = 0.0
    for i in range(1, ncol + 1):
        ltrue[i - 1] = 0
        for k in range(top_lev, pver + 1):
            if qc[_idx2(i, k, pcols)] >= qsmall or qi[_idx2(i, k, pcols)] >= qsmall or cmei[_idx2(i, k, pcols)] >= qsmall:
                ltrue[i - 1] = 1
    iter = 2
    deltat = deltat / float(iter)
    mtime = 1.0
    for _s172 in range(1, pver + 1):
        for _s171 in range(1, pcols + 1):
            rate1ord_cw2pr_st[_idx2(_s171, _s172, pcols)] = 0.0
    for i in range(1, ncol + 1):
        if ltrue[i - 1] == 0:
            for _s173 in range(1, pver + 1):
                tlat[_idx2(i, _s173, pcols)] = 0.0
            for _s174 in range(1, pver + 1):
                qvlat[_idx2(i, _s174, pcols)] = 0.0
            for _s175 in range(1, pver + 1):
                qctend[_idx2(i, _s175, pcols)] = 0.0
            for _s176 in range(1, pver + 1):
                qitend[_idx2(i, _s176, pcols)] = 0.0
            for _s177 in range(1, pver + 1):
                qnitend[_idx2(i, _s177, pcols)] = 0.0
            for _s178 in range(1, pver + 1):
                qrtend[_idx2(i, _s178, pcols)] = 0.0
            for _s179 in range(1, pver + 1):
                nctend[_idx2(i, _s179, pcols)] = 0.0
            for _s180 in range(1, pver + 1):
                nitend[_idx2(i, _s180, pcols)] = 0.0
            for _s181 in range(1, pver + 1):
                nrtend[_idx2(i, _s181, pcols)] = 0.0
            for _s182 in range(1, pver + 1):
                nstend[_idx2(i, _s182, pcols)] = 0.0
            prect[i - 1] = 0.0
            preci[i - 1] = 0.0
            for _s183 in range(1, pver + 1):
                qniic[_idx2(i, _s183, pcols)] = 0.0
            for _s184 in range(1, pver + 1):
                qric[_idx2(i, _s184, pcols)] = 0.0
            for _s185 in range(1, pver + 1):
                nsic[_idx2(i, _s185, pcols)] = 0.0
            for _s186 in range(1, pver + 1):
                nric[_idx2(i, _s186, pcols)] = 0.0
            for _s187 in range(1, pver + 1):
                rainrt[_idx2(i, _s187, pcols)] = 0.0
            for _s188 in range(1, pver + 1):
                qrtend_copy[_idx2(i, _s188, pcols)] = 0.0
            for _s189 in range(1, pver + 1):
                qnitend_copy[_idx2(i, _s189, pcols)] = 0.0
            continue
        for _s190 in range(1, pver + 1):
            qcsinksum_rate1ord[_s190 - 1] = 0.0
        for _s191 in range(1, pver + 1):
            qcsum_rate1ord[_s191 - 1] = 0.0
        for it in range(1, iter + 1):
            for _s192 in range(1, pver + 1):
                tlat[_idx2(i, _s192, pcols)] = 0.0
            for _s193 in range(1, pver + 1):
                qvlat[_idx2(i, _s193, pcols)] = 0.0
            for _s194 in range(1, pver + 1):
                qctend[_idx2(i, _s194, pcols)] = 0.0
            for _s195 in range(1, pver + 1):
                qitend[_idx2(i, _s195, pcols)] = 0.0
            for _s196 in range(1, pver + 1):
                qnitend[_idx2(i, _s196, pcols)] = 0.0
            for _s197 in range(1, pver + 1):
                qrtend[_idx2(i, _s197, pcols)] = 0.0
            for _s198 in range(1, pver + 1):
                nctend[_idx2(i, _s198, pcols)] = 0.0
            for _s199 in range(1, pver + 1):
                nitend[_idx2(i, _s199, pcols)] = 0.0
            for _s200 in range(1, pver + 1):
                nrtend[_idx2(i, _s200, pcols)] = 0.0
            for _s201 in range(1, pver + 1):
                nstend[_idx2(i, _s201, pcols)] = 0.0
            for _s202 in range(1, pver + 1):
                qrtend_copy[_idx2(i, _s202, pcols)] = 0.0
            for _s203 in range(1, pver + 1):
                qnitend_copy[_idx2(i, _s203, pcols)] = 0.0
            for _s204 in range(1, pver + 1):
                qniic[_idx2(i, _s204, pcols)] = 0.0
            for _s205 in range(1, pver + 1):
                qric[_idx2(i, _s205, pcols)] = 0.0
            for _s206 in range(1, pver + 1):
                nsic[_idx2(i, _s206, pcols)] = 0.0
            for _s207 in range(1, pver + 1):
                nric[_idx2(i, _s207, pcols)] = 0.0
            for _s208 in range(1, pver + 1):
                rainrt[_idx2(i, _s208, pcols)] = 0.0
            qrtot = 0.0
            nrtot = 0.0
            qstot = 0.0
            nstot = 0.0
            prect[i - 1] = 0.0
            preci[i - 1] = 0.0
            for k in range(top_lev, pver + 1):
                qcvar = relvar[_idx2(i, k, pcols)]
                cons2 = gamma(qcvar + 2.47)
                cons3 = gamma(qcvar)
                cons9 = gamma(qcvar + 2.0)
                cons10 = gamma(qcvar + 1.0)
                cons12 = gamma(qcvar + 1.15)
                cons15 = gamma(qcvar + bc / 3.0)
                cons18 = pow_r8(qcvar, 2.47)
                cons19 = pow_i(qcvar, 2)
                cons20 = pow_r8(qcvar, 1.15)
                cwml[_idx2(i, k, pcols)] = qc[_idx2(i, k, pcols)]
                cwmi[_idx2(i, k, pcols)] = qi[_idx2(i, k, pcols)]
                ums[k - 1] = 0.0
                uns[k - 1] = 0.0
                umr[k - 1] = 0.0
                unr[k - 1] = 0.0
                if k == top_lev:
                    cldmax[_idx2(i, k, pcols)] = cldm[_idx2(i, k, pcols)]
                elif do_clubb_sgs:
                    if qc[_idx2(i, k, pcols)] >= qsmall or qi[_idx2(i, k, pcols)] >= qsmall:
                        cldmax[_idx2(i, k, pcols)] = cldm[_idx2(i, k, pcols)]
                    else:
                        cldmax[_idx2(i, k, pcols)] = cldmax[_idx2(i, k - 1, pcols)]
                elif qric[_idx2(i, k - 1, pcols)] >= qsmall or qniic[_idx2(i, k - 1, pcols)] >= qsmall:
                    cldmax[_idx2(i, k, pcols)] = max(cldmax[_idx2(i, k - 1, pcols)], cldm[_idx2(i, k, pcols)])
                else:
                    cldmax[_idx2(i, k, pcols)] = cldm[_idx2(i, k, pcols)]
                if cmei[_idx2(i, k, pcols)] < 0.0 and qi[_idx2(i, k, pcols)] > qsmall and (cldm[_idx2(i, k, pcols)] > mincld):
                    nsubi[k - 1] = cmei[_idx2(i, k, pcols)] / qi[_idx2(i, k, pcols)] * ni[_idx2(i, k, pcols)] / cldm[_idx2(i, k, pcols)]
                else:
                    nsubi[k - 1] = 0.0
                nsubc[k - 1] = 0.0
                if do_cldice and dum2i[_idx2(i, k, pcols)] > 0.0 and (t[_idx2(i, k, pcols)] < tmelt - 5.0) and (relhum[_idx2(i, k, pcols)] * esl[_idx2(i, k, pcols)] / esi[_idx2(i, k, pcols)] > rhmini + 0.05):
                    nnuccd[k - 1] = (dum2i[_idx2(i, k, pcols)] - ni[_idx2(i, k, pcols)] / icldm[_idx2(i, k, pcols)]) / deltat * icldm[_idx2(i, k, pcols)]
                    nnuccd[k - 1] = max(nnuccd[k - 1], 0.0)
                    nimax = dum2i[_idx2(i, k, pcols)] * icldm[_idx2(i, k, pcols)]
                    mnuccd[k - 1] = nnuccd[k - 1] * mi0
                    cmei[_idx2(i, k, pcols)] = cmei[_idx2(i, k, pcols)] + mnuccd[k - 1] * mtime
                    qvi = wv_sat_svp_to_qsat_codon(esi[_idx2(i, k, pcols)], p[_idx2(i, k, pcols)], epsilo, omeps)
                    dqsidt = xxls * qvi / (rv * pow_i(t[_idx2(i, k, pcols)], 2))
                    abi = 1.0 + dqsidt * xxls / cpp
                    cmei[_idx2(i, k, pcols)] = min(cmei[_idx2(i, k, pcols)], (q[_idx2(i, k, pcols)] - qvi) / abi / deltat)
                    cmei[_idx2(i, k, pcols)] = cmei[_idx2(i, k, pcols)] * omsm
                else:
                    nnuccd[k - 1] = 0.0
                    nimax = 0.0
                    mnuccd[k - 1] = 0.0
                qcic[_idx2(i, k, pcols)] = min(cwml[_idx2(i, k, pcols)] / lcldm[_idx2(i, k, pcols)], 0.005)
                qiic[_idx2(i, k, pcols)] = min(cwmi[_idx2(i, k, pcols)] / icldm[_idx2(i, k, pcols)], 0.005)
                ncic[_idx2(i, k, pcols)] = max(nc[_idx2(i, k, pcols)] / lcldm[_idx2(i, k, pcols)], 0.0)
                niic[_idx2(i, k, pcols)] = max(ni[_idx2(i, k, pcols)] / icldm[_idx2(i, k, pcols)], 0.0)
                if qc[_idx2(i, k, pcols)] - berg[_idx2(i, k, pcols)] * deltat < qsmall:
                    qcic[_idx2(i, k, pcols)] = 0.0
                    ncic[_idx2(i, k, pcols)] = 0.0
                    if qc[_idx2(i, k, pcols)] - berg[_idx2(i, k, pcols)] * deltat < 0.0:
                        berg[_idx2(i, k, pcols)] = qc[_idx2(i, k, pcols)] / deltat * omsm
                if do_cldice and qi[_idx2(i, k, pcols)] + (cmei[_idx2(i, k, pcols)] + berg[_idx2(i, k, pcols)]) * deltat < qsmall:
                    qiic[_idx2(i, k, pcols)] = 0.0
                    niic[_idx2(i, k, pcols)] = 0.0
                    if qi[_idx2(i, k, pcols)] + (cmei[_idx2(i, k, pcols)] + berg[_idx2(i, k, pcols)]) * deltat < 0.0:
                        cmei[_idx2(i, k, pcols)] = (-qi[_idx2(i, k, pcols)] / deltat - berg[_idx2(i, k, pcols)]) * omsm
                cmeout[_idx2(i, k, pcols)] = cmeout[_idx2(i, k, pcols)] + cmei[_idx2(i, k, pcols)]
                if qcic[_idx2(i, k, pcols)] >= qsmall:
                    npccn[k - 1] = max(0.0, npccnin[_idx2(i, k, pcols)])
                    dum2l[_idx2(i, k, pcols)] = (nc[_idx2(i, k, pcols)] + npccn[k - 1] * deltat) / lcldm[_idx2(i, k, pcols)]
                    dum2l[_idx2(i, k, pcols)] = max(dum2l[_idx2(i, k, pcols)], cdnl / rho[_idx2(i, k, pcols)])
                    ncmax = dum2l[_idx2(i, k, pcols)] * lcldm[_idx2(i, k, pcols)]
                else:
                    npccn[k - 1] = 0.0
                    dum2l[_idx2(i, k, pcols)] = 0.0
                    ncmax = 0.0
                if qiic[_idx2(i, k, pcols)] >= qsmall:
                    niic[_idx2(i, k, pcols)] = min(niic[_idx2(i, k, pcols)], qiic[_idx2(i, k, pcols)] * 1e+20)
                    lami[k - 1] = pow_r8(cons1 * ci * niic[_idx2(i, k, pcols)] / qiic[_idx2(i, k, pcols)], 1.0 / di)
                    n0i[k - 1] = niic[_idx2(i, k, pcols)] * lami[k - 1]
                    if lami[k - 1] < lammini:
                        lami[k - 1] = lammini
                        n0i[k - 1] = pow_r8(lami[k - 1], di + 1.0) * qiic[_idx2(i, k, pcols)] / (ci * cons1)
                        niic[_idx2(i, k, pcols)] = n0i[k - 1] / lami[k - 1]
                    elif lami[k - 1] > lammaxi:
                        lami[k - 1] = lammaxi
                        n0i[k - 1] = pow_r8(lami[k - 1], di + 1.0) * qiic[_idx2(i, k, pcols)] / (ci * cons1)
                        niic[_idx2(i, k, pcols)] = n0i[k - 1] / lami[k - 1]
                else:
                    lami[k - 1] = 0.0
                    n0i[k - 1] = 0.0
                if qcic[_idx2(i, k, pcols)] >= qsmall:
                    ncic[_idx2(i, k, pcols)] = min(ncic[_idx2(i, k, pcols)], qcic[_idx2(i, k, pcols)] * 1e+20)
                    ncic[_idx2(i, k, pcols)] = max(ncic[_idx2(i, k, pcols)], cdnl / rho[_idx2(i, k, pcols)])
                    pgam[k - 1] = 0.0005714 * (ncic[_idx2(i, k, pcols)] / 1000000.0 * rho[_idx2(i, k, pcols)]) + 0.2714
                    pgam[k - 1] = 1.0 / pow_i(pgam[k - 1], 2) - 1.0
                    pgam[k - 1] = max(pgam[k - 1], 2.0)
                    pgam[k - 1] = min(pgam[k - 1], 15.0)
                    lamc[k - 1] = pow_r8(pi / 6.0 * rhow * ncic[_idx2(i, k, pcols)] * gamma(pgam[k - 1] + 4.0) / (qcic[_idx2(i, k, pcols)] * gamma(pgam[k - 1] + 1.0)), 1.0 / 3.0)
                    lammin = (pgam[k - 1] + 1.0) / 5e-05
                    lammax = (pgam[k - 1] + 1.0) / 2e-06
                    if lamc[k - 1] < lammin:
                        lamc[k - 1] = lammin
                        ncic[_idx2(i, k, pcols)] = 6.0 * pow_i(lamc[k - 1], 3) * qcic[_idx2(i, k, pcols)] * gamma(pgam[k - 1] + 1.0) / (pi * rhow * gamma(pgam[k - 1] + 4.0))
                    elif lamc[k - 1] > lammax:
                        lamc[k - 1] = lammax
                        ncic[_idx2(i, k, pcols)] = 6.0 * pow_i(lamc[k - 1], 3) * qcic[_idx2(i, k, pcols)] * gamma(pgam[k - 1] + 1.0) / (pi * rhow * gamma(pgam[k - 1] + 4.0))
                    cdist1[k - 1] = ncic[_idx2(i, k, pcols)] / gamma(pgam[k - 1] + 1.0)
                else:
                    lamc[k - 1] = 0.0
                    cdist1[k - 1] = 0.0
                if qcic[_idx2(i, k, pcols)] >= 1e-08:
                    if microp_uniform:
                        prc[k - 1] = 1350.0 * pow_r8(qcic[_idx2(i, k, pcols)], 2.47) * pow_r8(ncic[_idx2(i, k, pcols)] / 1000000.0 * rho[_idx2(i, k, pcols)], -1.79)
                        nprc[k - 1] = prc[k - 1] / (4.0 / 3.0 * pi * rhow * pow_i(2.5e-05, 3))
                        nprc1[k - 1] = prc[k - 1] / (qcic[_idx2(i, k, pcols)] / ncic[_idx2(i, k, pcols)])
                    else:
                        prc[k - 1] = cons2 / (cons3 * cons18) * 1350.0 * pow_r8(qcic[_idx2(i, k, pcols)], 2.47) * pow_r8(ncic[_idx2(i, k, pcols)] / 1000000.0 * rho[_idx2(i, k, pcols)], -1.79)
                        nprc[k - 1] = prc[k - 1] / cons22
                        nprc1[k - 1] = prc[k - 1] / (qcic[_idx2(i, k, pcols)] / ncic[_idx2(i, k, pcols)])
                else:
                    prc[k - 1] = 0.0
                    nprc[k - 1] = 0.0
                    nprc1[k - 1] = 0.0
                dum = 0.45
                dum1 = 0.45
                if k == top_lev:
                    qric[_idx2(i, k, pcols)] = prc[k - 1] * lcldm[_idx2(i, k, pcols)] * dz[_idx2(i, k, pcols)] / cldmax[_idx2(i, k, pcols)] / dum
                    nric[_idx2(i, k, pcols)] = nprc[k - 1] * lcldm[_idx2(i, k, pcols)] * dz[_idx2(i, k, pcols)] / cldmax[_idx2(i, k, pcols)] / dum
                else:
                    if qric[_idx2(i, k - 1, pcols)] >= qsmall:
                        dum = umr[k - 1 - 1]
                        dum1 = unr[k - 1 - 1]
                    if qric[_idx2(i, k - 1, pcols)] >= 1e-09 or qniic[_idx2(i, k - 1, pcols)] >= 1e-09:
                        nprc[k - 1] = 0.0
                    qric[_idx2(i, k, pcols)] = (rho[_idx2(i, k - 1, pcols)] * umr[k - 1 - 1] * qric[_idx2(i, k - 1, pcols)] * cldmax[_idx2(i, k - 1, pcols)] + rho[_idx2(i, k, pcols)] * dz[_idx2(i, k, pcols)] * ((pra[k - 1 - 1] + prc[k - 1]) * lcldm[_idx2(i, k, pcols)] + (pre[k - 1 - 1] - pracs[k - 1 - 1] - mnuccr[k - 1 - 1]) * cldmax[_idx2(i, k, pcols)])) / (dum * rho[_idx2(i, k, pcols)] * cldmax[_idx2(i, k, pcols)])
                    nric[_idx2(i, k, pcols)] = (rho[_idx2(i, k - 1, pcols)] * unr[k - 1 - 1] * nric[_idx2(i, k - 1, pcols)] * cldmax[_idx2(i, k - 1, pcols)] + rho[_idx2(i, k, pcols)] * dz[_idx2(i, k, pcols)] * (nprc[k - 1] * lcldm[_idx2(i, k, pcols)] + (nsubr[k - 1 - 1] - npracs[k - 1 - 1] - nnuccr[k - 1 - 1] + nragg[k - 1 - 1]) * cldmax[_idx2(i, k, pcols)])) / (dum1 * rho[_idx2(i, k, pcols)] * cldmax[_idx2(i, k, pcols)])
                if do_cldice:
                    if t[_idx2(i, k, pcols)] <= 273.15 and qiic[_idx2(i, k, pcols)] >= qsmall:
                        nprci[k - 1] = n0i[k - 1] / (lami[k - 1] * 180.0) * exp(-lami[k - 1] * dcs)
                        prci[k - 1] = pi * rhoi * n0i[k - 1] / (6.0 * 180.0) * (cons23 / lami[k - 1] + 3.0 * cons24 / pow_i(lami[k - 1], 2) + 6.0 * dcs / pow_i(lami[k - 1], 3) + 6.0 / pow_i(lami[k - 1], 4)) * exp(-lami[k - 1] * dcs)
                    else:
                        prci[k - 1] = 0.0
                        nprci[k - 1] = 0.0
                else:
                    prci[k - 1] = tnd_qsnow[_idx2(i, k, pcols)] / cldm[_idx2(i, k, pcols)]
                    nprci[k - 1] = tnd_nsnow[_idx2(i, k, pcols)] / cldm[_idx2(i, k, pcols)]
                dum = asn[_idx2(i, k, pcols)] * cons25
                dum1 = asn[_idx2(i, k, pcols)] * cons25
                if k == top_lev:
                    qniic[_idx2(i, k, pcols)] = prci[k - 1] * icldm[_idx2(i, k, pcols)] * dz[_idx2(i, k, pcols)] / cldmax[_idx2(i, k, pcols)] / dum
                    nsic[_idx2(i, k, pcols)] = nprci[k - 1] * icldm[_idx2(i, k, pcols)] * dz[_idx2(i, k, pcols)] / cldmax[_idx2(i, k, pcols)] / dum
                else:
                    if qniic[_idx2(i, k - 1, pcols)] >= qsmall:
                        dum = ums[k - 1 - 1]
                        dum1 = uns[k - 1 - 1]
                    qniic[_idx2(i, k, pcols)] = (rho[_idx2(i, k - 1, pcols)] * ums[k - 1 - 1] * qniic[_idx2(i, k - 1, pcols)] * cldmax[_idx2(i, k - 1, pcols)] + rho[_idx2(i, k, pcols)] * dz[_idx2(i, k, pcols)] * ((prci[k - 1] + prai[k - 1 - 1] + psacws[k - 1 - 1] + bergs[k - 1 - 1]) * icldm[_idx2(i, k, pcols)] + (prds[k - 1 - 1] + pracs[k - 1 - 1] + mnuccr[k - 1 - 1]) * cldmax[_idx2(i, k, pcols)])) / (dum * rho[_idx2(i, k, pcols)] * cldmax[_idx2(i, k, pcols)])
                    nsic[_idx2(i, k, pcols)] = (rho[_idx2(i, k - 1, pcols)] * uns[k - 1 - 1] * nsic[_idx2(i, k - 1, pcols)] * cldmax[_idx2(i, k - 1, pcols)] + rho[_idx2(i, k, pcols)] * dz[_idx2(i, k, pcols)] * (nprci[k - 1] * icldm[_idx2(i, k, pcols)] + (nsubs[k - 1 - 1] + nsagg[k - 1 - 1] + nnuccr[k - 1 - 1]) * cldmax[_idx2(i, k, pcols)])) / (dum1 * rho[_idx2(i, k, pcols)] * cldmax[_idx2(i, k, pcols)])
                if qniic[_idx2(i, k, pcols)] < qsmall:
                    qniic[_idx2(i, k, pcols)] = 0.0
                    nsic[_idx2(i, k, pcols)] = 0.0
                if qric[_idx2(i, k, pcols)] < qsmall:
                    qric[_idx2(i, k, pcols)] = 0.0
                    nric[_idx2(i, k, pcols)] = 0.0
                nric[_idx2(i, k, pcols)] = max(nric[_idx2(i, k, pcols)], 0.0)
                nsic[_idx2(i, k, pcols)] = max(nsic[_idx2(i, k, pcols)], 0.0)
                if qric[_idx2(i, k, pcols)] >= qsmall:
                    lamr[k - 1] = pow_r8(pi * rhow * nric[_idx2(i, k, pcols)] / qric[_idx2(i, k, pcols)], 1.0 / 3.0)
                    n0r[k - 1] = nric[_idx2(i, k, pcols)] * lamr[k - 1]
                    if lamr[k - 1] < lamminr:
                        lamr[k - 1] = lamminr
                        n0r[k - 1] = pow_i(lamr[k - 1], 4) * qric[_idx2(i, k, pcols)] / (pi * rhow)
                        nric[_idx2(i, k, pcols)] = n0r[k - 1] / lamr[k - 1]
                    elif lamr[k - 1] > lammaxr:
                        lamr[k - 1] = lammaxr
                        n0r[k - 1] = pow_i(lamr[k - 1], 4) * qric[_idx2(i, k, pcols)] / (pi * rhow)
                        nric[_idx2(i, k, pcols)] = n0r[k - 1] / lamr[k - 1]
                    unr[k - 1] = min(arn[_idx2(i, k, pcols)] * cons4 / pow_r8(lamr[k - 1], br), 9.1 * rhof[_idx2(i, k, pcols)])
                    umr[k - 1] = min(arn[_idx2(i, k, pcols)] * cons5 / (6.0 * pow_r8(lamr[k - 1], br)), 9.1 * rhof[_idx2(i, k, pcols)])
                else:
                    lamr[k - 1] = 0.0
                    n0r[k - 1] = 0.0
                    umr[k - 1] = 0.0
                    unr[k - 1] = 0.0
                if qniic[_idx2(i, k, pcols)] >= qsmall:
                    lams[k - 1] = pow_r8(cons6 * cs * nsic[_idx2(i, k, pcols)] / qniic[_idx2(i, k, pcols)], 1.0 / ds)
                    n0s[k - 1] = nsic[_idx2(i, k, pcols)] * lams[k - 1]
                    if lams[k - 1] < lammins:
                        lams[k - 1] = lammins
                        n0s[k - 1] = pow_r8(lams[k - 1], ds + 1.0) * qniic[_idx2(i, k, pcols)] / (cs * cons6)
                        nsic[_idx2(i, k, pcols)] = n0s[k - 1] / lams[k - 1]
                    elif lams[k - 1] > lammaxs:
                        lams[k - 1] = lammaxs
                        n0s[k - 1] = pow_r8(lams[k - 1], ds + 1.0) * qniic[_idx2(i, k, pcols)] / (cs * cons6)
                        nsic[_idx2(i, k, pcols)] = n0s[k - 1] / lams[k - 1]
                    ums[k - 1] = min(asn[_idx2(i, k, pcols)] * cons8 / (6.0 * pow_r8(lams[k - 1], bs)), 1.2 * rhof[_idx2(i, k, pcols)])
                    uns[k - 1] = min(asn[_idx2(i, k, pcols)] * cons7 / pow_r8(lams[k - 1], bs), 1.2 * rhof[_idx2(i, k, pcols)])
                else:
                    lams[k - 1] = 0.0
                    n0s[k - 1] = 0.0
                    ums[k - 1] = 0.0
                    uns[k - 1] = 0.0
                if not use_hetfrz_classnuc:
                    if do_cldice and qcic[_idx2(i, k, pcols)] >= qsmall and (t[_idx2(i, k, pcols)] < 269.15):
                        if microp_uniform:
                            mnuccc[k - 1] = pi * pi / 36.0 * rhow * cdist1[k - 1] * gamma(7.0 + pgam[k - 1]) * bimm * (exp(aimm * (273.15 - t[_idx2(i, k, pcols)])) - 1.0) / pow_i(lamc[k - 1], 3) / pow_i(lamc[k - 1], 3)
                            nnuccc[k - 1] = pi / 6.0 * cdist1[k - 1] * gamma(pgam[k - 1] + 4.0) * bimm * (exp(aimm * (273.15 - t[_idx2(i, k, pcols)])) - 1.0) / pow_i(lamc[k - 1], 3)
                        else:
                            mnuccc[k - 1] = cons9 / (cons3 * cons19) * pi * pi / 36.0 * rhow * cdist1[k - 1] * gamma(7.0 + pgam[k - 1]) * bimm * (exp(aimm * (273.15 - t[_idx2(i, k, pcols)])) - 1.0) / pow_i(lamc[k - 1], 3) / pow_i(lamc[k - 1], 3)
                            nnuccc[k - 1] = cons10 / (cons3 * qcvar) * pi / 6.0 * cdist1[k - 1] * gamma(pgam[k - 1] + 4.0) * bimm * (exp(aimm * (273.15 - t[_idx2(i, k, pcols)])) - 1.0) / pow_i(lamc[k - 1], 3)
                        tcnt = pow_r8(270.16 - t[_idx2(i, k, pcols)], 1.3)
                        viscosity = 1.8e-05 * pow_r8(t[_idx2(i, k, pcols)] / 298.0, 0.85)
                        mfp = 2.0 * viscosity / (p[_idx2(i, k, pcols)] * sqrt(8.0 * 0.02896 / (pi * 8.314409 * t[_idx2(i, k, pcols)])))
                        nslip1 = 1.0 + mfp / rndst[_idx3(i, k, 1, pcols, pver)] * (1.257 + 0.4 * exp(-(1.1 * rndst[_idx3(i, k, 1, pcols, pver)] / mfp)))
                        nslip2 = 1.0 + mfp / rndst[_idx3(i, k, 2, pcols, pver)] * (1.257 + 0.4 * exp(-(1.1 * rndst[_idx3(i, k, 2, pcols, pver)] / mfp)))
                        nslip3 = 1.0 + mfp / rndst[_idx3(i, k, 3, pcols, pver)] * (1.257 + 0.4 * exp(-(1.1 * rndst[_idx3(i, k, 3, pcols, pver)] / mfp)))
                        nslip4 = 1.0 + mfp / rndst[_idx3(i, k, 4, pcols, pver)] * (1.257 + 0.4 * exp(-(1.1 * rndst[_idx3(i, k, 4, pcols, pver)] / mfp)))
                        ndfaer1 = 1.381e-23 * t[_idx2(i, k, pcols)] * nslip1 / (6.0 * pi * viscosity * rndst[_idx3(i, k, 1, pcols, pver)])
                        ndfaer2 = 1.381e-23 * t[_idx2(i, k, pcols)] * nslip2 / (6.0 * pi * viscosity * rndst[_idx3(i, k, 2, pcols, pver)])
                        ndfaer3 = 1.381e-23 * t[_idx2(i, k, pcols)] * nslip3 / (6.0 * pi * viscosity * rndst[_idx3(i, k, 3, pcols, pver)])
                        ndfaer4 = 1.381e-23 * t[_idx2(i, k, pcols)] * nslip4 / (6.0 * pi * viscosity * rndst[_idx3(i, k, 4, pcols, pver)])
                        if microp_uniform:
                            mnucct[k - 1] = (ndfaer1 * (nacon[_idx3(i, k, 1, pcols, pver)] * tcnt) + ndfaer2 * (nacon[_idx3(i, k, 2, pcols, pver)] * tcnt) + ndfaer3 * (nacon[_idx3(i, k, 3, pcols, pver)] * tcnt) + ndfaer4 * (nacon[_idx3(i, k, 4, pcols, pver)] * tcnt)) * pi * pi / 3.0 * rhow * cdist1[k - 1] * gamma(pgam[k - 1] + 5.0) / pow_i(lamc[k - 1], 4)
                            nnucct[k - 1] = (ndfaer1 * (nacon[_idx3(i, k, 1, pcols, pver)] * tcnt) + ndfaer2 * (nacon[_idx3(i, k, 2, pcols, pver)] * tcnt) + ndfaer3 * (nacon[_idx3(i, k, 3, pcols, pver)] * tcnt) + ndfaer4 * (nacon[_idx3(i, k, 4, pcols, pver)] * tcnt)) * 2.0 * pi * cdist1[k - 1] * gamma(pgam[k - 1] + 2.0) / lamc[k - 1]
                        else:
                            mnucct[k - 1] = gamma(qcvar + 4.0 / 3.0) / (cons3 * pow_r8(qcvar, 4.0 / 3.0)) * (ndfaer1 * (nacon[_idx3(i, k, 1, pcols, pver)] * tcnt) + ndfaer2 * (nacon[_idx3(i, k, 2, pcols, pver)] * tcnt) + ndfaer3 * (nacon[_idx3(i, k, 3, pcols, pver)] * tcnt) + ndfaer4 * (nacon[_idx3(i, k, 4, pcols, pver)] * tcnt)) * pi * pi / 3.0 * rhow * cdist1[k - 1] * gamma(pgam[k - 1] + 5.0) / pow_i(lamc[k - 1], 4)
                            nnucct[k - 1] = gamma(qcvar + 1.0 / 3.0) / (cons3 * pow_r8(qcvar, 1.0 / 3.0)) * (ndfaer1 * (nacon[_idx3(i, k, 1, pcols, pver)] * tcnt) + ndfaer2 * (nacon[_idx3(i, k, 2, pcols, pver)] * tcnt) + ndfaer3 * (nacon[_idx3(i, k, 3, pcols, pver)] * tcnt) + ndfaer4 * (nacon[_idx3(i, k, 4, pcols, pver)] * tcnt)) * 2.0 * pi * cdist1[k - 1] * gamma(pgam[k - 1] + 2.0) / lamc[k - 1]
                        if nnuccc[k - 1] * lcldm[_idx2(i, k, pcols)] > nnuccd[k - 1]:
                            dum = nnuccd[k - 1] / (nnuccc[k - 1] * lcldm[_idx2(i, k, pcols)])
                            mnuccc[k - 1] = mnuccc[k - 1] * dum
                            nnuccc[k - 1] = nnuccd[k - 1] / lcldm[_idx2(i, k, pcols)]
                    else:
                        mnuccc[k - 1] = 0.0
                        nnuccc[k - 1] = 0.0
                        mnucct[k - 1] = 0.0
                        nnucct[k - 1] = 0.0
                elif do_cldice and qcic[_idx2(i, k, pcols)] >= qsmall:
                    con1 = 1.0 / pow_r8(1.333 * pi, 0.333)
                    r3lx = con1 * pow_r8(rho[_idx2(i, k, pcols)] * qcic[_idx2(i, k, pcols)] / (rhow * max(ncic[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)], 1000000.0)), 0.333)
                    r3lx = max(4e-06, r3lx)
                    mi0l = 4.0 / 3.0 * pi * rhow * pow_r8(r3lx, 3)
                    nnuccc[k - 1] = frzimm[_idx2(i, k, pcols)] * 1000000.0 / rho[_idx2(i, k, pcols)]
                    mnuccc[k - 1] = nnuccc[k - 1] * mi0l
                    nnucct[k - 1] = frzcnt[_idx2(i, k, pcols)] * 1000000.0 / rho[_idx2(i, k, pcols)]
                    mnucct[k - 1] = nnucct[k - 1] * mi0l
                    nnudep[k - 1] = frzdep[_idx2(i, k, pcols)] * 1000000.0 / rho[_idx2(i, k, pcols)]
                    mnudep[k - 1] = nnudep[k - 1] * mi0
                else:
                    nnuccc[k - 1] = 0.0
                    mnuccc[k - 1] = 0.0
                    nnucct[k - 1] = 0.0
                    mnucct[k - 1] = 0.0
                    nnudep[k - 1] = 0.0
                    mnudep[k - 1] = 0.0
                if qniic[_idx2(i, k, pcols)] >= qsmall and t[_idx2(i, k, pcols)] <= 273.15:
                    nsagg[k - 1] = -1108.0 * asn[_idx2(i, k, pcols)] * eii * pow_r8(pi, (1.0 - bs) / 3.0) * pow_r8(rhosn, (-2.0 - bs) / 3.0) * pow_r8(rho[_idx2(i, k, pcols)], (2.0 + bs) / 3.0) * pow_r8(qniic[_idx2(i, k, pcols)], (2.0 + bs) / 3.0) * pow_r8(nsic[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)], (4.0 - bs) / 3.0) / (4.0 * 720.0 * rho[_idx2(i, k, pcols)])
                else:
                    nsagg[k - 1] = 0.0
                if qniic[_idx2(i, k, pcols)] >= qsmall and t[_idx2(i, k, pcols)] <= tmelt and (qcic[_idx2(i, k, pcols)] >= qsmall):
                    dc0 = (pgam[k - 1] + 1.0) / lamc[k - 1]
                    ds0 = 1.0 / lams[k - 1]
                    dum = dc0 * dc0 * uns[k - 1] * rhow / (9.0 * mu[_idx2(i, k, pcols)] * ds0)
                    eci = dum * dum / ((dum + 0.4) * (dum + 0.4))
                    eci = max(eci, 0.0)
                    eci = min(eci, 1.0)
                    psacws[k - 1] = pi / 4.0 * asn[_idx2(i, k, pcols)] * qcic[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)] * n0s[k - 1] * eci * cons11 / pow_r8(lams[k - 1], bs + 3.0)
                    npsacws[k - 1] = pi / 4.0 * asn[_idx2(i, k, pcols)] * ncic[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)] * n0s[k - 1] * eci * cons11 / pow_r8(lams[k - 1], bs + 3.0)
                else:
                    psacws[k - 1] = 0.0
                    npsacws[k - 1] = 0.0
                if not do_cldice:
                    ni_secp = 0.0
                    nsacwi[k - 1] = 0.0
                    msacwi[k - 1] = 0.0
                elif t[_idx2(i, k, pcols)] < 270.16 and t[_idx2(i, k, pcols)] >= 268.16:
                    ni_secp = 350000000.0 * (270.16 - t[_idx2(i, k, pcols)]) / 2.0 * psacws[k - 1]
                    nsacwi[k - 1] = ni_secp
                    msacwi[k - 1] = min(ni_secp * mi0, psacws[k - 1])
                elif t[_idx2(i, k, pcols)] < 268.16 and t[_idx2(i, k, pcols)] >= 265.16:
                    ni_secp = 350000000.0 * (t[_idx2(i, k, pcols)] - 265.16) / 3.0 * psacws[k - 1]
                    nsacwi[k - 1] = ni_secp
                    msacwi[k - 1] = min(ni_secp * mi0, psacws[k - 1])
                else:
                    ni_secp = 0.0
                    nsacwi[k - 1] = 0.0
                    msacwi[k - 1] = 0.0
                psacws[k - 1] = max(0.0, psacws[k - 1] - ni_secp * mi0)
                if qric[_idx2(i, k, pcols)] >= 1e-08 and qniic[_idx2(i, k, pcols)] >= 1e-08 and (t[_idx2(i, k, pcols)] <= 273.15):
                    pracs[k - 1] = pi * pi * ecr * (pow_r8(pow_i(1.2 * umr[k - 1] - 0.95 * ums[k - 1], 2) + 0.08 * ums[k - 1] * umr[k - 1], 0.5) * rhow * rho[_idx2(i, k, pcols)] * n0r[k - 1] * n0s[k - 1] * (5.0 / (pow_i_f90_pracs(lamr[k - 1], 6) * lams[k - 1]) + 2.0 / (pow_i_f90_pracs(lamr[k - 1], 5) * pow_i_f90_pracs(lams[k - 1], 2)) + 0.5 / (pow_i_f90_pracs(lamr[k - 1], 4) * pow_i_f90_pracs(lams[k - 1], 3))))
                    npracs[k - 1] = pi / 2.0 * rho[_idx2(i, k, pcols)] * ecr * pow_r8(1.7 * pow_i(unr[k - 1] - uns[k - 1], 2) + 0.3 * unr[k - 1] * uns[k - 1], 0.5) * n0r[k - 1] * n0s[k - 1] * (1.0 / (pow_i(lamr[k - 1], 3) * lams[k - 1]) + 1.0 / (pow_i(lamr[k - 1], 2) * pow_i(lams[k - 1], 2)) + 1.0 / (lamr[k - 1] * pow_i(lams[k - 1], 3)))
                else:
                    pracs[k - 1] = 0.0
                    npracs[k - 1] = 0.0
                if t[_idx2(i, k, pcols)] < 269.15 and qric[_idx2(i, k, pcols)] >= qsmall:
                    mnuccr[k - 1] = 20.0 * pi * pi * rhow * nric[_idx2(i, k, pcols)] * bimm * (exp(aimm * (273.15 - t[_idx2(i, k, pcols)])) - 1.0) / pow_i(lamr[k - 1], 3) / pow_i(lamr[k - 1], 3)
                    nnuccr[k - 1] = pi * nric[_idx2(i, k, pcols)] * bimm * (exp(aimm * (273.15 - t[_idx2(i, k, pcols)])) - 1.0) / pow_i(lamr[k - 1], 3)
                else:
                    mnuccr[k - 1] = 0.0
                    nnuccr[k - 1] = 0.0
                if qric[_idx2(i, k, pcols)] >= qsmall and qcic[_idx2(i, k, pcols)] >= qsmall:
                    if microp_uniform:
                        pra[k - 1] = 67.0 * pow_r8(qcic[_idx2(i, k, pcols)] * qric[_idx2(i, k, pcols)], 1.15)
                        npra[k - 1] = pra[k - 1] / (qcic[_idx2(i, k, pcols)] / ncic[_idx2(i, k, pcols)])
                    else:
                        pra[k - 1] = accre_enhan[_idx2(i, k, pcols)] * (cons12 / (cons3 * cons20) * 67.0 * pow_r8(qcic[_idx2(i, k, pcols)] * qric[_idx2(i, k, pcols)], 1.15))
                        npra[k - 1] = pra[k - 1] / (qcic[_idx2(i, k, pcols)] / ncic[_idx2(i, k, pcols)])
                else:
                    pra[k - 1] = 0.0
                    npra[k - 1] = 0.0
                if qric[_idx2(i, k, pcols)] >= qsmall:
                    nragg[k - 1] = -8.0 * nric[_idx2(i, k, pcols)] * qric[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)]
                else:
                    nragg[k - 1] = 0.0
                if do_cldice and qniic[_idx2(i, k, pcols)] >= qsmall and (qiic[_idx2(i, k, pcols)] >= qsmall) and (t[_idx2(i, k, pcols)] <= 273.15):
                    prai[k - 1] = pi / 4.0 * asn[_idx2(i, k, pcols)] * qiic[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)] * n0s[k - 1] * eii * cons11 / pow_r8(lams[k - 1], bs + 3.0)
                    nprai[k - 1] = pi / 4.0 * asn[_idx2(i, k, pcols)] * niic[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)] * n0s[k - 1] * eii * cons11 / pow_r8(lams[k - 1], bs + 3.0)
                else:
                    prai[k - 1] = 0.0
                    nprai[k - 1] = 0.0
                pre[k - 1] = 0.0
                prds[k - 1] = 0.0
                if qcic[_idx2(i, k, pcols)] + qiic[_idx2(i, k, pcols)] < 1e-06 or cldmax[_idx2(i, k, pcols)] > lcldm[_idx2(i, k, pcols)]:
                    if qcic[_idx2(i, k, pcols)] + qiic[_idx2(i, k, pcols)] < 1e-06:
                        dum = 0.0
                    else:
                        dum = lcldm[_idx2(i, k, pcols)]
                    esn = wv_sat_svp_water_codon(t[_idx2(i, k, pcols)], wv_sat_idx)
                    qsn = wv_sat_svp_to_qsat_codon(esn, p[_idx2(i, k, pcols)], epsilo, omeps)
                    esl[_idx2(i, k, pcols)] = esn
                    esi[_idx2(i, k, pcols)] = wv_sat_svp_ice_codon(t[_idx2(i, k, pcols)], wv_sat_idx)
                    if t[_idx2(i, k, pcols)] > tmelt:
                        esi[_idx2(i, k, pcols)] = esl[_idx2(i, k, pcols)]
                    qclr = (q[_idx2(i, k, pcols)] - dum * qsn) / (1.0 - dum)
                    if qric[_idx2(i, k, pcols)] >= qsmall:
                        qvs = wv_sat_svp_to_qsat_codon(esl[_idx2(i, k, pcols)], p[_idx2(i, k, pcols)], epsilo, omeps)
                        dqsdt = xxlv * qvs / (rv * pow_i(t[_idx2(i, k, pcols)], 2))
                        ab = 1.0 + dqsdt * xxlv / cpp
                        epsr = 2.0 * pi * n0r[k - 1] * rho[_idx2(i, k, pcols)] * dv[_idx2(i, k, pcols)] * (f1r / (lamr[k - 1] * lamr[k - 1]) + f2r * pow_r8(arn[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)] / mu[_idx2(i, k, pcols)], 0.5) * pow_r8(sc[_idx2(i, k, pcols)], 1.0 / 3.0) * cons13 / pow_r8(lamr[k - 1], 5.0 / 2.0 + br / 2.0))
                        pre[k - 1] = epsr * (qclr - qvs) / ab
                        pre[k - 1] = min(pre[k - 1] * (cldmax[_idx2(i, k, pcols)] - dum), 0.0)
                        pre[k - 1] = pre[k - 1] / cldmax[_idx2(i, k, pcols)]
                        am_evp_st[_idx2(i, k, pcols)] = max(cldmax[_idx2(i, k, pcols)] - dum, 0.0)
                    if qniic[_idx2(i, k, pcols)] >= qsmall:
                        qvi = wv_sat_svp_to_qsat_codon(esi[_idx2(i, k, pcols)], p[_idx2(i, k, pcols)], epsilo, omeps)
                        dqsidt = xxls * qvi / (rv * pow_i(t[_idx2(i, k, pcols)], 2))
                        abi = 1.0 + dqsidt * xxls / cpp
                        epss = 2.0 * pi * n0s[k - 1] * rho[_idx2(i, k, pcols)] * dv[_idx2(i, k, pcols)] * (f1s / (lams[k - 1] * lams[k - 1]) + f2s * pow_r8(asn[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)] / mu[_idx2(i, k, pcols)], 0.5) * pow_r8(sc[_idx2(i, k, pcols)], 1.0 / 3.0) * cons14 / pow_r8(lams[k - 1], 5.0 / 2.0 + bs / 2.0))
                        prds[k - 1] = epss * (qclr - qvi) / abi
                        prds[k - 1] = min(prds[k - 1] * (cldmax[_idx2(i, k, pcols)] - dum), 0.0)
                        prds[k - 1] = prds[k - 1] / cldmax[_idx2(i, k, pcols)]
                        am_evp_st[_idx2(i, k, pcols)] = max(cldmax[_idx2(i, k, pcols)] - dum, 0.0)
                    qtmp = q[_idx2(i, k, pcols)] - (cmei[_idx2(i, k, pcols)] + (pre[k - 1] + prds[k - 1]) * cldmax[_idx2(i, k, pcols)]) * deltat
                    ttmp = t[_idx2(i, k, pcols)] + (pre[k - 1] * cldmax[_idx2(i, k, pcols)] * xxlv + (cmei[_idx2(i, k, pcols)] + prds[k - 1] * cldmax[_idx2(i, k, pcols)]) * xxls) * deltat / cpp
                    ttmp = max(180.0, min(ttmp, 323.0))
                    esn = wv_sat_svp_water_codon(ttmp, wv_sat_idx)
                    qsn = wv_sat_svp_to_qsat_codon(esn, p[_idx2(i, k, pcols)], epsilo, omeps)
                    if qtmp > qsn:
                        if pre[k - 1] + prds[k - 1] < -1e-20:
                            dum1 = pre[k - 1] / (pre[k - 1] + prds[k - 1])
                            qtmp = q[_idx2(i, k, pcols)] - cmei[_idx2(i, k, pcols)] * deltat
                            ttmp = t[_idx2(i, k, pcols)] + cmei[_idx2(i, k, pcols)] * xxls * deltat / cpp
                            esn = wv_sat_svp_water_codon(ttmp, wv_sat_idx)
                            qsn = wv_sat_svp_to_qsat_codon(esn, p[_idx2(i, k, pcols)], epsilo, omeps)
                            dum = (qtmp - qsn) / (1.0 + cons27 * qsn / (cpp * rv * pow_i(ttmp, 2)))
                            dum = min(dum, 0.0)
                            pre[k - 1] = dum * dum1 / deltat / cldmax[_idx2(i, k, pcols)]
                            esn = wv_sat_svp_ice_codon(ttmp, wv_sat_idx)
                            qsn = wv_sat_svp_to_qsat_codon(esn, p[_idx2(i, k, pcols)], epsilo, omeps)
                            dum = (qtmp - qsn) / (1.0 + cons28 * qsn / (cpp * rv * pow_i(ttmp, 2)))
                            dum = min(dum, 0.0)
                            prds[k - 1] = dum * (1.0 - dum1) / deltat / cldmax[_idx2(i, k, pcols)]
                if qniic[_idx2(i, k, pcols)] >= qsmall and qcic[_idx2(i, k, pcols)] >= qsmall and (t[_idx2(i, k, pcols)] < tmelt):
                    qvi = wv_sat_svp_to_qsat_codon(esi[_idx2(i, k, pcols)], p[_idx2(i, k, pcols)], epsilo, omeps)
                    qvs = wv_sat_svp_to_qsat_codon(esl[_idx2(i, k, pcols)], p[_idx2(i, k, pcols)], epsilo, omeps)
                    dqsidt = xxls * qvi / (rv * pow_i(t[_idx2(i, k, pcols)], 2))
                    abi = 1.0 + dqsidt * xxls / cpp
                    epss = 2.0 * pi * n0s[k - 1] * rho[_idx2(i, k, pcols)] * dv[_idx2(i, k, pcols)] * (f1s / (lams[k - 1] * lams[k - 1]) + f2s * pow_r8(asn[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)] / mu[_idx2(i, k, pcols)], 0.5) * pow_r8(sc[_idx2(i, k, pcols)], 1.0 / 3.0) * cons14 / pow_r8(lams[k - 1], 5.0 / 2.0 + bs / 2.0))
                    bergs[k - 1] = epss * (qvs - qvi) / abi
                else:
                    bergs[k - 1] = 0.0
                qce = qc[_idx2(i, k, pcols)] - berg[_idx2(i, k, pcols)] * deltat
                nce = nc[_idx2(i, k, pcols)] + npccn[k - 1] * deltat * mtime
                qie = qi[_idx2(i, k, pcols)] + (cmei[_idx2(i, k, pcols)] + berg[_idx2(i, k, pcols)]) * deltat
                nie = ni[_idx2(i, k, pcols)] + nnuccd[k - 1] * deltat * mtime
                dum = (prc[k - 1] + pra[k - 1] + mnuccc[k - 1] + mnucct[k - 1] + msacwi[k - 1] + psacws[k - 1] + bergs[k - 1]) * lcldm[_idx2(i, k, pcols)] * deltat
                if dum > qce:
                    ratio = qce / deltat / lcldm[_idx2(i, k, pcols)] / (prc[k - 1] + pra[k - 1] + mnuccc[k - 1] + mnucct[k - 1] + msacwi[k - 1] + psacws[k - 1] + bergs[k - 1]) * omsm
                    prc[k - 1] = prc[k - 1] * ratio
                    pra[k - 1] = pra[k - 1] * ratio
                    mnuccc[k - 1] = mnuccc[k - 1] * ratio
                    mnucct[k - 1] = mnucct[k - 1] * ratio
                    msacwi[k - 1] = msacwi[k - 1] * ratio
                    psacws[k - 1] = psacws[k - 1] * ratio
                    bergs[k - 1] = bergs[k - 1] * ratio
                dum = (nprc1[k - 1] + npra[k - 1] + nnuccc[k - 1] + nnucct[k - 1] + npsacws[k - 1] - nsubc[k - 1]) * lcldm[_idx2(i, k, pcols)] * deltat
                if dum > nce:
                    ratio = nce / deltat / ((nprc1[k - 1] + npra[k - 1] + nnuccc[k - 1] + nnucct[k - 1] + npsacws[k - 1] - nsubc[k - 1]) * lcldm[_idx2(i, k, pcols)]) * omsm
                    nprc1[k - 1] = nprc1[k - 1] * ratio
                    npra[k - 1] = npra[k - 1] * ratio
                    nnuccc[k - 1] = nnuccc[k - 1] * ratio
                    nnucct[k - 1] = nnucct[k - 1] * ratio
                    npsacws[k - 1] = npsacws[k - 1] * ratio
                    nsubc[k - 1] = nsubc[k - 1] * ratio
                if do_cldice:
                    frztmp = -mnuccc[k - 1] - mnucct[k - 1] - msacwi[k - 1]
                    if use_hetfrz_classnuc:
                        frztmp = -mnuccc[k - 1] - mnucct[k - 1] - mnudep[k - 1] - msacwi[k - 1]
                    dum = (frztmp * lcldm[_idx2(i, k, pcols)] + (prci[k - 1] + prai[k - 1]) * icldm[_idx2(i, k, pcols)]) * deltat
                    if dum > qie:
                        frztmp = mnuccc[k - 1] + mnucct[k - 1] + msacwi[k - 1]
                        if use_hetfrz_classnuc:
                            frztmp = mnuccc[k - 1] + mnucct[k - 1] + mnudep[k - 1] + msacwi[k - 1]
                        ratio = (qie / deltat + frztmp * lcldm[_idx2(i, k, pcols)]) / ((prci[k - 1] + prai[k - 1]) * icldm[_idx2(i, k, pcols)]) * omsm
                        prci[k - 1] = prci[k - 1] * ratio
                        prai[k - 1] = prai[k - 1] * ratio
                    frztmp = -nnucct[k - 1] - nsacwi[k - 1]
                    if use_hetfrz_classnuc:
                        frztmp = -nnucct[k - 1] - nnuccc[k - 1] - nnudep[k - 1] - nsacwi[k - 1]
                    dum = (frztmp * lcldm[_idx2(i, k, pcols)] + (nprci[k - 1] + nprai[k - 1] - nsubi[k - 1]) * icldm[_idx2(i, k, pcols)]) * deltat
                    if dum > nie:
                        frztmp = nnucct[k - 1] + nsacwi[k - 1]
                        if use_hetfrz_classnuc:
                            frztmp = nnucct[k - 1] + nnuccc[k - 1] + nnudep[k - 1] + nsacwi[k - 1]
                        ratio = (nie / deltat + frztmp * lcldm[_idx2(i, k, pcols)]) / ((nprci[k - 1] + nprai[k - 1] - nsubi[k - 1]) * icldm[_idx2(i, k, pcols)]) * omsm
                        nprci[k - 1] = nprci[k - 1] * ratio
                        nprai[k - 1] = nprai[k - 1] * ratio
                        nsubi[k - 1] = nsubi[k - 1] * ratio
                if ((prc[k - 1] + pra[k - 1]) * lcldm[_idx2(i, k, pcols)] + (-mnuccr[k - 1] + pre[k - 1] - pracs[k - 1]) * cldmax[_idx2(i, k, pcols)]) * dz[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)] + qrtot < 0.0:
                    if -pre[k - 1] + pracs[k - 1] + mnuccr[k - 1] >= qsmall:
                        ratio = (qrtot / (dz[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)]) + (prc[k - 1] + pra[k - 1]) * lcldm[_idx2(i, k, pcols)]) / ((-pre[k - 1] + pracs[k - 1] + mnuccr[k - 1]) * cldmax[_idx2(i, k, pcols)]) * omsm
                        pre[k - 1] = pre[k - 1] * ratio
                        pracs[k - 1] = pracs[k - 1] * ratio
                        mnuccr[k - 1] = mnuccr[k - 1] * ratio
                nsubr[k - 1] = 0.0
                if (nprc[k - 1] * lcldm[_idx2(i, k, pcols)] + (-nnuccr[k - 1] + nsubr[k - 1] - npracs[k - 1] + nragg[k - 1]) * cldmax[_idx2(i, k, pcols)]) * dz[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)] + nrtot < 0.0:
                    if -nsubr[k - 1] - nragg[k - 1] + npracs[k - 1] + nnuccr[k - 1] >= qsmall:
                        ratio = (nrtot / (dz[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)]) + nprc[k - 1] * lcldm[_idx2(i, k, pcols)]) / ((-nsubr[k - 1] - nragg[k - 1] + npracs[k - 1] + nnuccr[k - 1]) * cldmax[_idx2(i, k, pcols)]) * omsm
                        nsubr[k - 1] = nsubr[k - 1] * ratio
                        npracs[k - 1] = npracs[k - 1] * ratio
                        nnuccr[k - 1] = nnuccr[k - 1] * ratio
                        nragg[k - 1] = nragg[k - 1] * ratio
                if ((bergs[k - 1] + psacws[k - 1]) * lcldm[_idx2(i, k, pcols)] + (prai[k - 1] + prci[k - 1]) * icldm[_idx2(i, k, pcols)] + (pracs[k - 1] + mnuccr[k - 1] + prds[k - 1]) * cldmax[_idx2(i, k, pcols)]) * dz[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)] + qstot < 0.0:
                    if -prds[k - 1] >= qsmall:
                        ratio = (qstot / (dz[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)]) + (bergs[k - 1] + psacws[k - 1]) * lcldm[_idx2(i, k, pcols)] + (prai[k - 1] + prci[k - 1]) * icldm[_idx2(i, k, pcols)] + (pracs[k - 1] + mnuccr[k - 1]) * cldmax[_idx2(i, k, pcols)]) / (-prds[k - 1] * cldmax[_idx2(i, k, pcols)]) * omsm
                        prds[k - 1] = prds[k - 1] * ratio
                nsubs[k - 1] = 0.0
                if (nprci[k - 1] * icldm[_idx2(i, k, pcols)] + (nnuccr[k - 1] + nsubs[k - 1] + nsagg[k - 1]) * cldmax[_idx2(i, k, pcols)]) * dz[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)] + nstot < 0.0:
                    if -nsubs[k - 1] - nsagg[k - 1] >= qsmall:
                        ratio = (nstot / (dz[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)]) + nprci[k - 1] * icldm[_idx2(i, k, pcols)] + nnuccr[k - 1] * cldmax[_idx2(i, k, pcols)]) / ((-nsubs[k - 1] - nsagg[k - 1]) * cldmax[_idx2(i, k, pcols)]) * omsm
                        nsubs[k - 1] = nsubs[k - 1] * ratio
                        nsagg[k - 1] = nsagg[k - 1] * ratio
                qvlat[_idx2(i, k, pcols)] = qvlat[_idx2(i, k, pcols)] - (pre[k - 1] + prds[k - 1]) * cldmax[_idx2(i, k, pcols)] - cmei[_idx2(i, k, pcols)]
                tlat[_idx2(i, k, pcols)] = tlat[_idx2(i, k, pcols)] + (pre[k - 1] * cldmax[_idx2(i, k, pcols)] * xxlv + (prds[k - 1] * cldmax[_idx2(i, k, pcols)] + cmei[_idx2(i, k, pcols)]) * xxls + ((bergs[k - 1] + psacws[k - 1] + mnuccc[k - 1] + mnucct[k - 1] + msacwi[k - 1]) * lcldm[_idx2(i, k, pcols)] + (mnuccr[k - 1] + pracs[k - 1]) * cldmax[_idx2(i, k, pcols)] + berg[_idx2(i, k, pcols)]) * xlf)
                qctend[_idx2(i, k, pcols)] = qctend[_idx2(i, k, pcols)] + (-pra[k - 1] - prc[k - 1] - mnuccc[k - 1] - mnucct[k - 1] - msacwi[k - 1] - psacws[k - 1] - bergs[k - 1]) * lcldm[_idx2(i, k, pcols)] - berg[_idx2(i, k, pcols)]
                if do_cldice:
                    frztmp = mnuccc[k - 1] + mnucct[k - 1] + msacwi[k - 1]
                    if use_hetfrz_classnuc:
                        frztmp = mnuccc[k - 1] + mnucct[k - 1] + mnudep[k - 1] + msacwi[k - 1]
                    qitend[_idx2(i, k, pcols)] = qitend[_idx2(i, k, pcols)] + frztmp * lcldm[_idx2(i, k, pcols)] + (-prci[k - 1] - prai[k - 1]) * icldm[_idx2(i, k, pcols)] + cmei[_idx2(i, k, pcols)] + berg[_idx2(i, k, pcols)]
                qrtend[_idx2(i, k, pcols)] = qrtend[_idx2(i, k, pcols)] + (pra[k - 1] + prc[k - 1]) * lcldm[_idx2(i, k, pcols)] + (pre[k - 1] - pracs[k - 1] - mnuccr[k - 1]) * cldmax[_idx2(i, k, pcols)]
                qnitend[_idx2(i, k, pcols)] = qnitend[_idx2(i, k, pcols)] + (prai[k - 1] + prci[k - 1]) * icldm[_idx2(i, k, pcols)] + (psacws[k - 1] + bergs[k - 1]) * lcldm[_idx2(i, k, pcols)] + (prds[k - 1] + pracs[k - 1] + mnuccr[k - 1]) * cldmax[_idx2(i, k, pcols)]
                qrtend_copy[_idx2(i, k, pcols)] = qrtend_copy[_idx2(i, k, pcols)] + qrtend[_idx2(i, k, pcols)]
                qnitend_copy[_idx2(i, k, pcols)] = qnitend_copy[_idx2(i, k, pcols)] + qnitend[_idx2(i, k, pcols)]
                cmeiout[_idx2(i, k, pcols)] = cmeiout[_idx2(i, k, pcols)] + cmei[_idx2(i, k, pcols)]
                evapsnow[_idx2(i, k, pcols)] = evapsnow[_idx2(i, k, pcols)] - prds[k - 1] * cldmax[_idx2(i, k, pcols)]
                nevapr[_idx2(i, k, pcols)] = nevapr[_idx2(i, k, pcols)] - pre[k - 1] * cldmax[_idx2(i, k, pcols)]
                nevapr2[_idx2(i, k, pcols)] = nevapr2[_idx2(i, k, pcols)] - pre[k - 1] * cldmax[_idx2(i, k, pcols)]
                prain[_idx2(i, k, pcols)] = prain[_idx2(i, k, pcols)] + (pra[k - 1] + prc[k - 1]) * lcldm[_idx2(i, k, pcols)] + (-pracs[k - 1] - mnuccr[k - 1]) * cldmax[_idx2(i, k, pcols)]
                prodsnow[_idx2(i, k, pcols)] = prodsnow[_idx2(i, k, pcols)] + (prai[k - 1] + prci[k - 1]) * icldm[_idx2(i, k, pcols)] + (psacws[k - 1] + bergs[k - 1]) * lcldm[_idx2(i, k, pcols)] + (pracs[k - 1] + mnuccr[k - 1]) * cldmax[_idx2(i, k, pcols)]
                qcsinksum_rate1ord[k - 1] = qcsinksum_rate1ord[k - 1] + (pra[k - 1] + prc[k - 1] + psacws[k - 1]) * lcldm[_idx2(i, k, pcols)]
                qcsum_rate1ord[k - 1] = qcsum_rate1ord[k - 1] + qc[_idx2(i, k, pcols)]
                prao[_idx2(i, k, pcols)] = prao[_idx2(i, k, pcols)] + pra[k - 1] * lcldm[_idx2(i, k, pcols)]
                prco[_idx2(i, k, pcols)] = prco[_idx2(i, k, pcols)] + prc[k - 1] * lcldm[_idx2(i, k, pcols)]
                mnuccco[_idx2(i, k, pcols)] = mnuccco[_idx2(i, k, pcols)] + mnuccc[k - 1] * lcldm[_idx2(i, k, pcols)]
                mnuccto[_idx2(i, k, pcols)] = mnuccto[_idx2(i, k, pcols)] + mnucct[k - 1] * lcldm[_idx2(i, k, pcols)]
                mnuccdo[_idx2(i, k, pcols)] = mnuccdo[_idx2(i, k, pcols)] + mnuccd[k - 1] * lcldm[_idx2(i, k, pcols)]
                msacwio[_idx2(i, k, pcols)] = msacwio[_idx2(i, k, pcols)] + msacwi[k - 1] * lcldm[_idx2(i, k, pcols)]
                psacwso[_idx2(i, k, pcols)] = psacwso[_idx2(i, k, pcols)] + psacws[k - 1] * lcldm[_idx2(i, k, pcols)]
                bergso[_idx2(i, k, pcols)] = bergso[_idx2(i, k, pcols)] + bergs[k - 1] * lcldm[_idx2(i, k, pcols)]
                bergo[_idx2(i, k, pcols)] = bergo[_idx2(i, k, pcols)] + berg[_idx2(i, k, pcols)]
                prcio[_idx2(i, k, pcols)] = prcio[_idx2(i, k, pcols)] + prci[k - 1] * icldm[_idx2(i, k, pcols)]
                praio[_idx2(i, k, pcols)] = praio[_idx2(i, k, pcols)] + prai[k - 1] * icldm[_idx2(i, k, pcols)]
                mnuccro[_idx2(i, k, pcols)] = mnuccro[_idx2(i, k, pcols)] + mnuccr[k - 1] * cldmax[_idx2(i, k, pcols)]
                pracso[_idx2(i, k, pcols)] = pracso[_idx2(i, k, pcols)] + pracs[k - 1] * cldmax[_idx2(i, k, pcols)]
                preo[_idx2(i, k, pcols)] = preo[_idx2(i, k, pcols)] + pre[k - 1] * cldmax[_idx2(i, k, pcols)]
                prdso[_idx2(i, k, pcols)] = prdso[_idx2(i, k, pcols)] + prds[k - 1] * cldmax[_idx2(i, k, pcols)]
                nctend[_idx2(i, k, pcols)] = nctend[_idx2(i, k, pcols)] + npccn[k - 1] * mtime + (-nnuccc[k - 1] - nnucct[k - 1] - npsacws[k - 1] + nsubc[k - 1] - npra[k - 1] - nprc1[k - 1]) * lcldm[_idx2(i, k, pcols)]
                if do_cldice:
                    frztmp = nnucct[k - 1] + nsacwi[k - 1]
                    if use_hetfrz_classnuc:
                        frztmp = nnucct[k - 1] + nnuccc[k - 1] + nnudep[k - 1] + nsacwi[k - 1]
                    nitend[_idx2(i, k, pcols)] = nitend[_idx2(i, k, pcols)] + nnuccd[k - 1] * mtime + frztmp * lcldm[_idx2(i, k, pcols)] + (nsubi[k - 1] - nprci[k - 1] - nprai[k - 1]) * icldm[_idx2(i, k, pcols)]
                nstend[_idx2(i, k, pcols)] = nstend[_idx2(i, k, pcols)] + (nsubs[k - 1] + nsagg[k - 1] + nnuccr[k - 1]) * cldmax[_idx2(i, k, pcols)] + nprci[k - 1] * icldm[_idx2(i, k, pcols)]
                nrtend[_idx2(i, k, pcols)] = nrtend[_idx2(i, k, pcols)] + nprc[k - 1] * lcldm[_idx2(i, k, pcols)] + (nsubr[k - 1] - npracs[k - 1] - nnuccr[k - 1] + nragg[k - 1]) * cldmax[_idx2(i, k, pcols)]
                if nctend[_idx2(i, k, pcols)] > 0.0 and nc[_idx2(i, k, pcols)] + nctend[_idx2(i, k, pcols)] * deltat > ncmax:
                    nctend[_idx2(i, k, pcols)] = max(0.0, (ncmax - nc[_idx2(i, k, pcols)]) / deltat)
                if do_cldice and nitend[_idx2(i, k, pcols)] > 0.0 and (ni[_idx2(i, k, pcols)] + nitend[_idx2(i, k, pcols)] * deltat > nimax):
                    nitend[_idx2(i, k, pcols)] = max(0.0, (nimax - ni[_idx2(i, k, pcols)]) / deltat)
                if qric[_idx2(i, k, pcols)] >= qsmall:
                    if k == top_lev:
                        qric[_idx2(i, k, pcols)] = qrtend[_idx2(i, k, pcols)] * dz[_idx2(i, k, pcols)] / cldmax[_idx2(i, k, pcols)] / umr[k - 1]
                        nric[_idx2(i, k, pcols)] = nrtend[_idx2(i, k, pcols)] * dz[_idx2(i, k, pcols)] / cldmax[_idx2(i, k, pcols)] / unr[k - 1]
                    else:
                        qric[_idx2(i, k, pcols)] = (rho[_idx2(i, k - 1, pcols)] * umr[k - 1 - 1] * qric[_idx2(i, k - 1, pcols)] * cldmax[_idx2(i, k - 1, pcols)] + rho[_idx2(i, k, pcols)] * dz[_idx2(i, k, pcols)] * qrtend[_idx2(i, k, pcols)]) / (umr[k - 1] * rho[_idx2(i, k, pcols)] * cldmax[_idx2(i, k, pcols)])
                        nric[_idx2(i, k, pcols)] = (rho[_idx2(i, k - 1, pcols)] * unr[k - 1 - 1] * nric[_idx2(i, k - 1, pcols)] * cldmax[_idx2(i, k - 1, pcols)] + rho[_idx2(i, k, pcols)] * dz[_idx2(i, k, pcols)] * nrtend[_idx2(i, k, pcols)]) / (unr[k - 1] * rho[_idx2(i, k, pcols)] * cldmax[_idx2(i, k, pcols)])
                else:
                    qric[_idx2(i, k, pcols)] = 0.0
                    nric[_idx2(i, k, pcols)] = 0.0
                if qniic[_idx2(i, k, pcols)] >= qsmall:
                    if k == top_lev:
                        qniic[_idx2(i, k, pcols)] = qnitend[_idx2(i, k, pcols)] * dz[_idx2(i, k, pcols)] / cldmax[_idx2(i, k, pcols)] / ums[k - 1]
                        nsic[_idx2(i, k, pcols)] = nstend[_idx2(i, k, pcols)] * dz[_idx2(i, k, pcols)] / cldmax[_idx2(i, k, pcols)] / uns[k - 1]
                    else:
                        qniic[_idx2(i, k, pcols)] = (rho[_idx2(i, k - 1, pcols)] * ums[k - 1 - 1] * qniic[_idx2(i, k - 1, pcols)] * cldmax[_idx2(i, k - 1, pcols)] + rho[_idx2(i, k, pcols)] * dz[_idx2(i, k, pcols)] * qnitend[_idx2(i, k, pcols)]) / (ums[k - 1] * rho[_idx2(i, k, pcols)] * cldmax[_idx2(i, k, pcols)])
                        nsic[_idx2(i, k, pcols)] = (rho[_idx2(i, k - 1, pcols)] * uns[k - 1 - 1] * nsic[_idx2(i, k - 1, pcols)] * cldmax[_idx2(i, k - 1, pcols)] + rho[_idx2(i, k, pcols)] * dz[_idx2(i, k, pcols)] * nstend[_idx2(i, k, pcols)]) / (uns[k - 1] * rho[_idx2(i, k, pcols)] * cldmax[_idx2(i, k, pcols)])
                else:
                    qniic[_idx2(i, k, pcols)] = 0.0
                    nsic[_idx2(i, k, pcols)] = 0.0
                prect[i - 1] = prect[i - 1] + (qrtend[_idx2(i, k, pcols)] * dz[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)] + qnitend[_idx2(i, k, pcols)] * dz[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)]) / rhow
                preci[i - 1] = preci[i - 1] + qnitend[_idx2(i, k, pcols)] * dz[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)] / rhow
                rainrt[_idx2(i, k, pcols)] = qric[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)] * umr[k - 1] / rhow * 3600.0 * 1000.0
                qrtot = max(qrtot + qrtend[_idx2(i, k, pcols)] * dz[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)], 0.0)
                qstot = max(qstot + qnitend[_idx2(i, k, pcols)] * dz[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)], 0.0)
                nrtot = max(nrtot + nrtend[_idx2(i, k, pcols)] * dz[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)], 0.0)
                nstot = max(nstot + nstend[_idx2(i, k, pcols)] * dz[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)], 0.0)
                if t[_idx2(i, k, pcols)] + tlat[_idx2(i, k, pcols)] / cpp * deltat > 275.15:
                    if qstot > 0.0:
                        dum = -xlf / cpp * qstot / (dz[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)])
                        if t[_idx2(i, k, pcols)] + tlat[_idx2(i, k, pcols)] / cpp * deltat + dum < 275.15:
                            dum = (t[_idx2(i, k, pcols)] + tlat[_idx2(i, k, pcols)] / cpp * deltat - 275.15) * cpp / xlf
                            dum = dum / (xlf / cpp * qstot / (dz[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)]))
                            dum = max(0.0, dum)
                            dum = min(1.0, dum)
                        else:
                            dum = 1.0
                        qric[_idx2(i, k, pcols)] = qric[_idx2(i, k, pcols)] + dum * qniic[_idx2(i, k, pcols)]
                        nric[_idx2(i, k, pcols)] = nric[_idx2(i, k, pcols)] + dum * nsic[_idx2(i, k, pcols)]
                        qniic[_idx2(i, k, pcols)] = (1.0 - dum) * qniic[_idx2(i, k, pcols)]
                        nsic[_idx2(i, k, pcols)] = (1.0 - dum) * nsic[_idx2(i, k, pcols)]
                        meltso[_idx2(i, k, pcols)] = meltso[_idx2(i, k, pcols)] + dum * qstot * g
                        tmp = -xlf * dum * qstot / (dz[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)])
                        meltsdt[_idx2(i, k, pcols)] = meltsdt[_idx2(i, k, pcols)] + tmp
                        tlat[_idx2(i, k, pcols)] = tlat[_idx2(i, k, pcols)] + tmp
                        qrtot = qrtot + dum * qstot
                        nrtot = nrtot + dum * nstot
                        qstot = (1.0 - dum) * qstot
                        nstot = (1.0 - dum) * nstot
                        preci[i - 1] = (1.0 - dum) * preci[i - 1]
                if t[_idx2(i, k, pcols)] + tlat[_idx2(i, k, pcols)] / cpp * deltat < tmelt - 5.0:
                    if qrtot > 0.0:
                        dum = xlf / cpp * qrtot / (dz[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)])
                        if t[_idx2(i, k, pcols)] + tlat[_idx2(i, k, pcols)] / cpp * deltat + dum > tmelt - 5.0:
                            dum = -(t[_idx2(i, k, pcols)] + tlat[_idx2(i, k, pcols)] / cpp * deltat - (tmelt - 5.0)) * cpp / xlf
                            dum = dum / (xlf / cpp * qrtot / (dz[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)]))
                            dum = max(0.0, dum)
                            dum = min(1.0, dum)
                        else:
                            dum = 1.0
                        qniic[_idx2(i, k, pcols)] = qniic[_idx2(i, k, pcols)] + dum * qric[_idx2(i, k, pcols)]
                        nsic[_idx2(i, k, pcols)] = nsic[_idx2(i, k, pcols)] + dum * nric[_idx2(i, k, pcols)]
                        qric[_idx2(i, k, pcols)] = (1.0 - dum) * qric[_idx2(i, k, pcols)]
                        nric[_idx2(i, k, pcols)] = (1.0 - dum) * nric[_idx2(i, k, pcols)]
                        frzro[_idx2(i, k, pcols)] = frzro[_idx2(i, k, pcols)] + dum * qrtot * g
                        tmp = xlf * dum * qrtot / (dz[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)])
                        frzrdt[_idx2(i, k, pcols)] = frzrdt[_idx2(i, k, pcols)] + tmp
                        tlat[_idx2(i, k, pcols)] = tlat[_idx2(i, k, pcols)] + tmp
                        qstot = qstot + dum * qrtot
                        qrtot = (1.0 - dum) * qrtot
                        nstot = nstot + dum * nrtot
                        nrtot = (1.0 - dum) * nrtot
                        preci[i - 1] = preci[i - 1] + dum * (prect[i - 1] - preci[i - 1])
                if qniic[_idx2(i, k, pcols)] < qsmall:
                    qniic[_idx2(i, k, pcols)] = 0.0
                    nsic[_idx2(i, k, pcols)] = 0.0
                if qric[_idx2(i, k, pcols)] < qsmall:
                    qric[_idx2(i, k, pcols)] = 0.0
                    nric[_idx2(i, k, pcols)] = 0.0
                nric[_idx2(i, k, pcols)] = max(nric[_idx2(i, k, pcols)], 0.0)
                nsic[_idx2(i, k, pcols)] = max(nsic[_idx2(i, k, pcols)], 0.0)
                if qric[_idx2(i, k, pcols)] >= qsmall:
                    lamr[k - 1] = pow_r8(pi * rhow * nric[_idx2(i, k, pcols)] / qric[_idx2(i, k, pcols)], 1.0 / 3.0)
                    n0r[k - 1] = nric[_idx2(i, k, pcols)] * lamr[k - 1]
                    if lamr[k - 1] < lamminr:
                        lamr[k - 1] = lamminr
                        n0r[k - 1] = pow_i(lamr[k - 1], 4) * qric[_idx2(i, k, pcols)] / (pi * rhow)
                        nric[_idx2(i, k, pcols)] = n0r[k - 1] / lamr[k - 1]
                    elif lamr[k - 1] > lammaxr:
                        lamr[k - 1] = lammaxr
                        n0r[k - 1] = pow_i(lamr[k - 1], 4) * qric[_idx2(i, k, pcols)] / (pi * rhow)
                        nric[_idx2(i, k, pcols)] = n0r[k - 1] / lamr[k - 1]
                    unr[k - 1] = min(arn[_idx2(i, k, pcols)] * cons4 / pow_r8(lamr[k - 1], br), 9.1 * rhof[_idx2(i, k, pcols)])
                    umr[k - 1] = min(arn[_idx2(i, k, pcols)] * cons5 / (6.0 * pow_r8(lamr[k - 1], br)), 9.1 * rhof[_idx2(i, k, pcols)])
                else:
                    lamr[k - 1] = 0.0
                    n0r[k - 1] = 0.0
                    umr[k - 1] = 0.0
                    unr[k - 1] = 0.0
                if lamr[k - 1] > 0.0:
                    artmp = n0r[k - 1] * pi / (2.0 * pow_r8(lamr[k - 1], 3.0))
                else:
                    artmp = 0.0
                if lamc[k - 1] > 0.0:
                    actmp = cdist1[k - 1] * pi * gamma(pgam[k - 1] + 3.0) / (4.0 * pow_r8(lamc[k - 1], 2.0))
                else:
                    actmp = 0.0
                if actmp > 0 or artmp > 0:
                    rercld[_idx2(i, k, pcols)] = rercld[_idx2(i, k, pcols)] + 3.0 * (qric[_idx2(i, k, pcols)] + qcic[_idx2(i, k, pcols)]) / (4.0 * rhow * (actmp + artmp))
                    arcld[_idx2(i, k, pcols)] = arcld[_idx2(i, k, pcols)] + 1.0
                if qniic[_idx2(i, k, pcols)] >= qsmall:
                    lams[k - 1] = pow_r8(cons6 * cs * nsic[_idx2(i, k, pcols)] / qniic[_idx2(i, k, pcols)], 1.0 / ds)
                    n0s[k - 1] = nsic[_idx2(i, k, pcols)] * lams[k - 1]
                    if lams[k - 1] < lammins:
                        lams[k - 1] = lammins
                        n0s[k - 1] = pow_r8(lams[k - 1], ds + 1.0) * qniic[_idx2(i, k, pcols)] / (cs * cons6)
                        nsic[_idx2(i, k, pcols)] = n0s[k - 1] / lams[k - 1]
                    elif lams[k - 1] > lammaxs:
                        lams[k - 1] = lammaxs
                        n0s[k - 1] = pow_r8(lams[k - 1], ds + 1.0) * qniic[_idx2(i, k, pcols)] / (cs * cons6)
                        nsic[_idx2(i, k, pcols)] = n0s[k - 1] / lams[k - 1]
                    ums[k - 1] = min(asn[_idx2(i, k, pcols)] * cons8 / (6.0 * pow_r8(lams[k - 1], bs)), 1.2 * rhof[_idx2(i, k, pcols)])
                    uns[k - 1] = min(asn[_idx2(i, k, pcols)] * cons7 / pow_r8(lams[k - 1], bs), 1.2 * rhof[_idx2(i, k, pcols)])
                else:
                    lams[k - 1] = 0.0
                    n0s[k - 1] = 0.0
                    ums[k - 1] = 0.0
                    uns[k - 1] = 0.0
                qrout[_idx2(i, k, pcols)] = qrout[_idx2(i, k, pcols)] + qric[_idx2(i, k, pcols)] * cldmax[_idx2(i, k, pcols)]
                qsout[_idx2(i, k, pcols)] = qsout[_idx2(i, k, pcols)] + qniic[_idx2(i, k, pcols)] * cldmax[_idx2(i, k, pcols)]
                nrout[_idx2(i, k, pcols)] = nrout[_idx2(i, k, pcols)] + nric[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)] * cldmax[_idx2(i, k, pcols)]
                nsout[_idx2(i, k, pcols)] = nsout[_idx2(i, k, pcols)] + nsic[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)] * cldmax[_idx2(i, k, pcols)]
                tlat1[_idx2(i, k, pcols)] = tlat1[_idx2(i, k, pcols)] + tlat[_idx2(i, k, pcols)]
                qvlat1[_idx2(i, k, pcols)] = qvlat1[_idx2(i, k, pcols)] + qvlat[_idx2(i, k, pcols)]
                qctend1[_idx2(i, k, pcols)] = qctend1[_idx2(i, k, pcols)] + qctend[_idx2(i, k, pcols)]
                qitend1[_idx2(i, k, pcols)] = qitend1[_idx2(i, k, pcols)] + qitend[_idx2(i, k, pcols)]
                nctend1[_idx2(i, k, pcols)] = nctend1[_idx2(i, k, pcols)] + nctend[_idx2(i, k, pcols)]
                nitend1[_idx2(i, k, pcols)] = nitend1[_idx2(i, k, pcols)] + nitend[_idx2(i, k, pcols)]
                t[_idx2(i, k, pcols)] = t[_idx2(i, k, pcols)] + tlat[_idx2(i, k, pcols)] * deltat / cpp
                q[_idx2(i, k, pcols)] = q[_idx2(i, k, pcols)] + qvlat[_idx2(i, k, pcols)] * deltat
                qc[_idx2(i, k, pcols)] = qc[_idx2(i, k, pcols)] + qctend[_idx2(i, k, pcols)] * deltat
                qi[_idx2(i, k, pcols)] = qi[_idx2(i, k, pcols)] + qitend[_idx2(i, k, pcols)] * deltat
                nc[_idx2(i, k, pcols)] = nc[_idx2(i, k, pcols)] + nctend[_idx2(i, k, pcols)] * deltat
                ni[_idx2(i, k, pcols)] = ni[_idx2(i, k, pcols)] + nitend[_idx2(i, k, pcols)] * deltat
                rainrt1[_idx2(i, k, pcols)] = rainrt1[_idx2(i, k, pcols)] + rainrt[_idx2(i, k, pcols)]
                if arcld[_idx2(i, k, pcols)] > 0.0:
                    rercld[_idx2(i, k, pcols)] = rercld[_idx2(i, k, pcols)] / arcld[_idx2(i, k, pcols)]
                rflx[_idx2(i, 1, pcols)] = 0.0
                sflx[_idx2(i, 1, pcols)] = 0.0
                rflx[_idx2(i, k + 1, pcols)] = qrout[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)] * umr[k - 1]
                sflx[_idx2(i, k + 1, pcols)] = qsout[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)] * ums[k - 1]
                rflx1[_idx2(i, k + 1, pcols)] = rflx1[_idx2(i, k + 1, pcols)] + rflx[_idx2(i, k + 1, pcols)]
                sflx1[_idx2(i, k + 1, pcols)] = sflx1[_idx2(i, k + 1, pcols)] + sflx[_idx2(i, k + 1, pcols)]
            prect1[i - 1] = prect1[i - 1] + prect[i - 1]
            preci1[i - 1] = preci1[i - 1] + preci[i - 1]
        for k in range(top_lev, pver + 1):
            rate1ord_cw2pr_st[_idx2(i, k, pcols)] = qcsinksum_rate1ord[k - 1] / max(qcsum_rate1ord[k - 1], 1e-30)
    deltat = deltat * float(iter)
    for i in range(1, ncol + 1):
        skip_to_500 = False
        if ltrue[i - 1] == 0:
            for k in range(1, top_lev - 1 + 1):
                effc[_idx2(i, k, pcols)] = 0.0
                effi[_idx2(i, k, pcols)] = 0.0
                effc_fn[_idx2(i, k, pcols)] = 0.0
                lamcrad[_idx2(i, k, pcols)] = 0.0
                pgamrad[_idx2(i, k, pcols)] = 0.0
                deffi[_idx2(i, k, pcols)] = 0.0
            for k in range(top_lev, pver + 1):
                effc[_idx2(i, k, pcols)] = 10.0
                effi[_idx2(i, k, pcols)] = 25.0
                effc_fn[_idx2(i, k, pcols)] = 10.0
                lamcrad[_idx2(i, k, pcols)] = 0.0
                pgamrad[_idx2(i, k, pcols)] = 0.0
                deffi[_idx2(i, k, pcols)] = 0.0
            skip_to_500 = True
        if not skip_to_500:
            nstep = 1
            prect[i - 1] = prect1[i - 1] / float(iter)
            preci[i - 1] = preci1[i - 1] / float(iter)
            for k in range(top_lev, pver + 1):
                t[_idx2(i, k, pcols)] = t1[_idx2(i, k, pcols)]
                q[_idx2(i, k, pcols)] = q1[_idx2(i, k, pcols)]
                qc[_idx2(i, k, pcols)] = qc1[_idx2(i, k, pcols)]
                qi[_idx2(i, k, pcols)] = qi1[_idx2(i, k, pcols)]
                nc[_idx2(i, k, pcols)] = nc1[_idx2(i, k, pcols)]
                ni[_idx2(i, k, pcols)] = ni1[_idx2(i, k, pcols)]
                tlat[_idx2(i, k, pcols)] = tlat1[_idx2(i, k, pcols)] / float(iter)
                qvlat[_idx2(i, k, pcols)] = qvlat1[_idx2(i, k, pcols)] / float(iter)
                qctend[_idx2(i, k, pcols)] = qctend1[_idx2(i, k, pcols)] / float(iter)
                qitend[_idx2(i, k, pcols)] = qitend1[_idx2(i, k, pcols)] / float(iter)
                nctend[_idx2(i, k, pcols)] = nctend1[_idx2(i, k, pcols)] / float(iter)
                nitend[_idx2(i, k, pcols)] = nitend1[_idx2(i, k, pcols)] / float(iter)
                rainrt[_idx2(i, k, pcols)] = rainrt1[_idx2(i, k, pcols)] / float(iter)
                rflx[_idx2(i, k + 1, pcols)] = rflx1[_idx2(i, k + 1, pcols)] / float(iter)
                sflx[_idx2(i, k + 1, pcols)] = sflx1[_idx2(i, k + 1, pcols)] / float(iter)
                qrout[_idx2(i, k, pcols)] = qrout[_idx2(i, k, pcols)] / float(iter)
                qsout[_idx2(i, k, pcols)] = qsout[_idx2(i, k, pcols)] / float(iter)
                nrout[_idx2(i, k, pcols)] = nrout[_idx2(i, k, pcols)] / float(iter)
                nsout[_idx2(i, k, pcols)] = nsout[_idx2(i, k, pcols)] / float(iter)
                nevapr[_idx2(i, k, pcols)] = nevapr[_idx2(i, k, pcols)] / float(iter)
                nevapr2[_idx2(i, k, pcols)] = nevapr2[_idx2(i, k, pcols)] / float(iter)
                evapsnow[_idx2(i, k, pcols)] = evapsnow[_idx2(i, k, pcols)] / float(iter)
                prain[_idx2(i, k, pcols)] = prain[_idx2(i, k, pcols)] / float(iter)
                prodsnow[_idx2(i, k, pcols)] = prodsnow[_idx2(i, k, pcols)] / float(iter)
                cmeout[_idx2(i, k, pcols)] = cmeout[_idx2(i, k, pcols)] / float(iter)
                cmeiout[_idx2(i, k, pcols)] = cmeiout[_idx2(i, k, pcols)] / float(iter)
                meltsdt[_idx2(i, k, pcols)] = meltsdt[_idx2(i, k, pcols)] / float(iter)
                frzrdt[_idx2(i, k, pcols)] = frzrdt[_idx2(i, k, pcols)] / float(iter)
                prao[_idx2(i, k, pcols)] = prao[_idx2(i, k, pcols)] / float(iter)
                prco[_idx2(i, k, pcols)] = prco[_idx2(i, k, pcols)] / float(iter)
                mnuccco[_idx2(i, k, pcols)] = mnuccco[_idx2(i, k, pcols)] / float(iter)
                mnuccto[_idx2(i, k, pcols)] = mnuccto[_idx2(i, k, pcols)] / float(iter)
                msacwio[_idx2(i, k, pcols)] = msacwio[_idx2(i, k, pcols)] / float(iter)
                psacwso[_idx2(i, k, pcols)] = psacwso[_idx2(i, k, pcols)] / float(iter)
                bergso[_idx2(i, k, pcols)] = bergso[_idx2(i, k, pcols)] / float(iter)
                bergo[_idx2(i, k, pcols)] = bergo[_idx2(i, k, pcols)] / float(iter)
                prcio[_idx2(i, k, pcols)] = prcio[_idx2(i, k, pcols)] / float(iter)
                praio[_idx2(i, k, pcols)] = praio[_idx2(i, k, pcols)] / float(iter)
                mnuccro[_idx2(i, k, pcols)] = mnuccro[_idx2(i, k, pcols)] / float(iter)
                pracso[_idx2(i, k, pcols)] = pracso[_idx2(i, k, pcols)] / float(iter)
                mnuccdo[_idx2(i, k, pcols)] = mnuccdo[_idx2(i, k, pcols)] / float(iter)
                preo[_idx2(i, k, pcols)] = preo[_idx2(i, k, pcols)] / float(iter)
                prdso[_idx2(i, k, pcols)] = prdso[_idx2(i, k, pcols)] / float(iter)
                frzro[_idx2(i, k, pcols)] = frzro[_idx2(i, k, pcols)] / float(iter)
                meltso[_idx2(i, k, pcols)] = meltso[_idx2(i, k, pcols)] / float(iter)
                wtprelat[_idx2(i, k, pcols)] = tlat[_idx2(i, k, pcols)]
                nevapr[_idx2(i, k, pcols)] = nevapr[_idx2(i, k, pcols)] + evapsnow[_idx2(i, k, pcols)]
                prer_evap[_idx2(i, k, pcols)] = nevapr2[_idx2(i, k, pcols)]
                prain[_idx2(i, k, pcols)] = prain[_idx2(i, k, pcols)] + prodsnow[_idx2(i, k, pcols)]
                dumc[_idx2(i, k, pcols)] = (qc[_idx2(i, k, pcols)] + qctend[_idx2(i, k, pcols)] * deltat) / lcldm[_idx2(i, k, pcols)]
                dumi[_idx2(i, k, pcols)] = (qi[_idx2(i, k, pcols)] + qitend[_idx2(i, k, pcols)] * deltat) / icldm[_idx2(i, k, pcols)]
                dumnc[_idx2(i, k, pcols)] = max((nc[_idx2(i, k, pcols)] + nctend[_idx2(i, k, pcols)] * deltat) / lcldm[_idx2(i, k, pcols)], 0.0)
                dumni[_idx2(i, k, pcols)] = max((ni[_idx2(i, k, pcols)] + nitend[_idx2(i, k, pcols)] * deltat) / icldm[_idx2(i, k, pcols)], 0.0)
                if dumi[_idx2(i, k, pcols)] >= qsmall:
                    dumni[_idx2(i, k, pcols)] = min(dumni[_idx2(i, k, pcols)], dumi[_idx2(i, k, pcols)] * 1e+20)
                    lami[k - 1] = pow_r8(cons1 * ci * dumni[_idx2(i, k, pcols)] / dumi[_idx2(i, k, pcols)], 1.0 / di)
                    lami[k - 1] = max(lami[k - 1], lammini)
                    lami[k - 1] = min(lami[k - 1], lammaxi)
                else:
                    lami[k - 1] = 0.0
                if dumc[_idx2(i, k, pcols)] >= qsmall:
                    dumnc[_idx2(i, k, pcols)] = min(dumnc[_idx2(i, k, pcols)], dumc[_idx2(i, k, pcols)] * 1e+20)
                    dumnc[_idx2(i, k, pcols)] = max(dumnc[_idx2(i, k, pcols)], cdnl / rho[_idx2(i, k, pcols)])
                    pgam[k - 1] = 0.0005714 * (ncic[_idx2(i, k, pcols)] / 1000000.0 * rho[_idx2(i, k, pcols)]) + 0.2714
                    pgam[k - 1] = 1.0 / pow_i(pgam[k - 1], 2) - 1.0
                    pgam[k - 1] = max(pgam[k - 1], 2.0)
                    pgam[k - 1] = min(pgam[k - 1], 15.0)
                    lamc[k - 1] = pow_r8(pi / 6.0 * rhow * dumnc[_idx2(i, k, pcols)] * gamma(pgam[k - 1] + 4.0) / (dumc[_idx2(i, k, pcols)] * gamma(pgam[k - 1] + 1.0)), 1.0 / 3.0)
                    lammin = (pgam[k - 1] + 1.0) / 5e-05
                    lammax = (pgam[k - 1] + 1.0) / 2e-06
                    lamc[k - 1] = max(lamc[k - 1], lammin)
                    lamc[k - 1] = min(lamc[k - 1], lammax)
                else:
                    lamc[k - 1] = 0.0
                if dumc[_idx2(i, k, pcols)] >= qsmall:
                    unc = acn[_idx2(i, k, pcols)] * gamma(1.0 + bc + pgam[k - 1]) / (pow_r8(lamc[k - 1], bc) * gamma(pgam[k - 1] + 1.0))
                    umc = acn[_idx2(i, k, pcols)] * gamma(4.0 + bc + pgam[k - 1]) / (pow_r8(lamc[k - 1], bc) * gamma(pgam[k - 1] + 4.0))
                    vtrmc[_idx2(i, k, pcols)] = umc
                else:
                    umc = 0.0
                    unc = 0.0
                if dumi[_idx2(i, k, pcols)] >= qsmall:
                    uni = ain[_idx2(i, k, pcols)] * cons16 / pow_r8(lami[k - 1], bi)
                    umi = ain[_idx2(i, k, pcols)] * cons17 / (6.0 * pow_r8(lami[k - 1], bi))
                    uni = min(uni, 1.2 * rhof[_idx2(i, k, pcols)])
                    umi = min(umi, 1.2 * rhof[_idx2(i, k, pcols)])
                    vtrmi[_idx2(i, k, pcols)] = umi
                else:
                    umi = 0.0
                    uni = 0.0
                fi[k - 1] = g * rho[_idx2(i, k, pcols)] * umi
                fni[k - 1] = g * rho[_idx2(i, k, pcols)] * uni
                fc[k - 1] = g * rho[_idx2(i, k, pcols)] * umc
                fnc[k - 1] = g * rho[_idx2(i, k, pcols)] * unc
                wtfc[_idx2(i, k, pcols)] = fc[k - 1]
                wtfi[_idx2(i, k, pcols)] = fi[k - 1]
                rgvm = max(fi[k - 1], fc[k - 1], fni[k - 1], fnc[k - 1])
                nstep = max(int(rgvm * deltat / pdel[_idx2(i, k, pcols)] + 1.0), nstep)
                dumc[_idx2(i, k, pcols)] = qc[_idx2(i, k, pcols)] + qctend[_idx2(i, k, pcols)] * deltat
                dumi[_idx2(i, k, pcols)] = qi[_idx2(i, k, pcols)] + qitend[_idx2(i, k, pcols)] * deltat
                dumnc[_idx2(i, k, pcols)] = max(nc[_idx2(i, k, pcols)] + nctend[_idx2(i, k, pcols)] * deltat, 0.0)
                dumni[_idx2(i, k, pcols)] = max(ni[_idx2(i, k, pcols)] + nitend[_idx2(i, k, pcols)] * deltat, 0.0)
                if dumc[_idx2(i, k, pcols)] < qsmall:
                    dumnc[_idx2(i, k, pcols)] = 0.0
                if dumi[_idx2(i, k, pcols)] < qsmall:
                    dumni[_idx2(i, k, pcols)] = 0.0
            micro_mg1_0_sedimentation_fallout_codon(
                i,
                pcols,
                pver,
                top_lev,
                nstep,
                1 if do_cldice else 0,
                deltat,
                g,
                xxlv,
                xxls,
                ptrs[9],
                ptrs[134],
                ptrs[133],
                ptrs[21],
                ptrs[23],
                ptrs[20],
                ptrs[22],
                ptrs[51],
                ptrs[52],
                ptrs[19],
                ptrs[46],
                ptrs[45],
                ptrs[18],
                ptrs[208],
                ptrs[209],
                ptrs[206],
                ptrs[207],
                ptrs[218],
                ptrs[219],
                ptrs[216],
                ptrs[217],
                ptrs[226],
                ptrs[227],
                ptrs[224],
                ptrs[225],
                ptrs[27],
                ptrs[28],
            )
            for k in range(top_lev, pver + 1):
                dumc[_idx2(i, k, pcols)] = max(qc[_idx2(i, k, pcols)] + qctend[_idx2(i, k, pcols)] * deltat, 0.0)
                dumi[_idx2(i, k, pcols)] = max(qi[_idx2(i, k, pcols)] + qitend[_idx2(i, k, pcols)] * deltat, 0.0)
                dumnc[_idx2(i, k, pcols)] = max(nc[_idx2(i, k, pcols)] + nctend[_idx2(i, k, pcols)] * deltat, 0.0)
                dumni[_idx2(i, k, pcols)] = max(ni[_idx2(i, k, pcols)] + nitend[_idx2(i, k, pcols)] * deltat, 0.0)
                if dumc[_idx2(i, k, pcols)] < qsmall:
                    dumnc[_idx2(i, k, pcols)] = 0.0
                if dumi[_idx2(i, k, pcols)] < qsmall:
                    dumni[_idx2(i, k, pcols)] = 0.0
                if do_cldice:
                    if t[_idx2(i, k, pcols)] + tlat[_idx2(i, k, pcols)] / cpp * deltat > tmelt:
                        if dumi[_idx2(i, k, pcols)] > 0.0:
                            dum = -dumi[_idx2(i, k, pcols)] * xlf / cpp
                            if t[_idx2(i, k, pcols)] + tlat[_idx2(i, k, pcols)] / cpp * deltat + dum < tmelt:
                                dum = (t[_idx2(i, k, pcols)] + tlat[_idx2(i, k, pcols)] / cpp * deltat - tmelt) * cpp / xlf
                                dum = dum / dumi[_idx2(i, k, pcols)] * xlf / cpp
                                dum = max(0.0, dum)
                                dum = min(1.0, dum)
                            else:
                                dum = 1.0
                            qctend[_idx2(i, k, pcols)] = qctend[_idx2(i, k, pcols)] + dum * dumi[_idx2(i, k, pcols)] / deltat
                            melto[_idx2(i, k, pcols)] = dum * dumi[_idx2(i, k, pcols)] / deltat
                            nctend[_idx2(i, k, pcols)] = nctend[_idx2(i, k, pcols)] + 3.0 * dum * dumi[_idx2(i, k, pcols)] / deltat / (4.0 * pi * 5.12e-16 * rhow)
                            qitend[_idx2(i, k, pcols)] = ((1.0 - dum) * dumi[_idx2(i, k, pcols)] - qi[_idx2(i, k, pcols)]) / deltat
                            nitend[_idx2(i, k, pcols)] = ((1.0 - dum) * dumni[_idx2(i, k, pcols)] - ni[_idx2(i, k, pcols)]) / deltat
                            tlat[_idx2(i, k, pcols)] = tlat[_idx2(i, k, pcols)] - xlf * dum * dumi[_idx2(i, k, pcols)] / deltat
                            wtpostlat[_idx2(i, k, pcols)] = wtpostlat[_idx2(i, k, pcols)] - xlf * dum * dumi[_idx2(i, k, pcols)] / deltat
                    if t[_idx2(i, k, pcols)] + tlat[_idx2(i, k, pcols)] / cpp * deltat < 233.15:
                        if dumc[_idx2(i, k, pcols)] > 0.0:
                            dum = dumc[_idx2(i, k, pcols)] * xlf / cpp
                            if t[_idx2(i, k, pcols)] + tlat[_idx2(i, k, pcols)] / cpp * deltat + dum > 233.15:
                                dum = -(t[_idx2(i, k, pcols)] + tlat[_idx2(i, k, pcols)] / cpp * deltat - 233.15) * cpp / xlf
                                dum = dum / dumc[_idx2(i, k, pcols)] * xlf / cpp
                                dum = max(0.0, dum)
                                dum = min(1.0, dum)
                            else:
                                dum = 1.0
                            qitend[_idx2(i, k, pcols)] = qitend[_idx2(i, k, pcols)] + dum * dumc[_idx2(i, k, pcols)] / deltat
                            homoo[_idx2(i, k, pcols)] = dum * dumc[_idx2(i, k, pcols)] / deltat
                            nitend[_idx2(i, k, pcols)] = nitend[_idx2(i, k, pcols)] + dum * 3.0 * dumc[_idx2(i, k, pcols)] / (4.0 * 3.14 * 1.563e-14 * 500.0) / deltat
                            qctend[_idx2(i, k, pcols)] = ((1.0 - dum) * dumc[_idx2(i, k, pcols)] - qc[_idx2(i, k, pcols)]) / deltat
                            nctend[_idx2(i, k, pcols)] = ((1.0 - dum) * dumnc[_idx2(i, k, pcols)] - nc[_idx2(i, k, pcols)]) / deltat
                            tlat[_idx2(i, k, pcols)] = tlat[_idx2(i, k, pcols)] + xlf * dum * dumc[_idx2(i, k, pcols)] / deltat
                            wtpostlat[_idx2(i, k, pcols)] = wtpostlat[_idx2(i, k, pcols)] + xlf * dum * dumc[_idx2(i, k, pcols)] / deltat
                    qtmp = q[_idx2(i, k, pcols)] + qvlat[_idx2(i, k, pcols)] * deltat
                    ttmp = t[_idx2(i, k, pcols)] + tlat[_idx2(i, k, pcols)] / cpp * deltat
                    esn = wv_sat_svp_water_codon(ttmp, wv_sat_idx)
                    qsn = wv_sat_svp_to_qsat_codon(esn, p[_idx2(i, k, pcols)], epsilo, omeps)
                    if qtmp > qsn and qsn > 0:
                        dum = (qtmp - qsn) / (1.0 + cons27 * qsn / (cpp * rv * pow_i(ttmp, 2))) / deltat
                        cmeout[_idx2(i, k, pcols)] = cmeout[_idx2(i, k, pcols)] + dum
                        if ttmp > 268.15:
                            dum1 = 0.0
                        elif ttmp < 238.15:
                            dum1 = 1.0
                        else:
                            dum1 = (268.15 - ttmp) / 30.0
                        dum = (qtmp - qsn) / (1.0 + pow_i(xxls * dum1 + xxlv * (1.0 - dum1), 2) * qsn / (cpp * rv * pow_i(ttmp, 2))) / deltat
                        qctend[_idx2(i, k, pcols)] = qctend[_idx2(i, k, pcols)] + dum * (1.0 - dum1)
                        qcreso[_idx2(i, k, pcols)] = dum * (1.0 - dum1)
                        qitend[_idx2(i, k, pcols)] = qitend[_idx2(i, k, pcols)] + dum * dum1
                        qireso[_idx2(i, k, pcols)] = dum * dum1
                        qvlat[_idx2(i, k, pcols)] = qvlat[_idx2(i, k, pcols)] - dum
                        qvres[_idx2(i, k, pcols)] = -dum
                        tlat[_idx2(i, k, pcols)] = tlat[_idx2(i, k, pcols)] + dum * (1.0 - dum1) * xxlv + dum * dum1 * xxls
                        wtpostlat[_idx2(i, k, pcols)] = wtpostlat[_idx2(i, k, pcols)] + (dum * (1.0 - dum1) * xxlv + dum * dum1 * xxls)
                dumc[_idx2(i, k, pcols)] = max(qc[_idx2(i, k, pcols)] + qctend[_idx2(i, k, pcols)] * deltat, 0.0) / lcldm[_idx2(i, k, pcols)]
                dumi[_idx2(i, k, pcols)] = max(qi[_idx2(i, k, pcols)] + qitend[_idx2(i, k, pcols)] * deltat, 0.0) / icldm[_idx2(i, k, pcols)]
                dumnc[_idx2(i, k, pcols)] = max(nc[_idx2(i, k, pcols)] + nctend[_idx2(i, k, pcols)] * deltat, 0.0) / lcldm[_idx2(i, k, pcols)]
                dumni[_idx2(i, k, pcols)] = max(ni[_idx2(i, k, pcols)] + nitend[_idx2(i, k, pcols)] * deltat, 0.0) / icldm[_idx2(i, k, pcols)]
                dumc[_idx2(i, k, pcols)] = min(dumc[_idx2(i, k, pcols)], 0.005)
                dumi[_idx2(i, k, pcols)] = min(dumi[_idx2(i, k, pcols)], 0.005)
                if dumi[_idx2(i, k, pcols)] >= qsmall:
                    dumni[_idx2(i, k, pcols)] = min(dumni[_idx2(i, k, pcols)], dumi[_idx2(i, k, pcols)] * 1e+20)
                    lami[k - 1] = pow_r8(cons1 * ci * dumni[_idx2(i, k, pcols)] / dumi[_idx2(i, k, pcols)], 1.0 / di)
                    if lami[k - 1] < lammini:
                        lami[k - 1] = lammini
                        n0i[k - 1] = pow_r8(lami[k - 1], di + 1.0) * dumi[_idx2(i, k, pcols)] / (ci * cons1)
                        niic[_idx2(i, k, pcols)] = n0i[k - 1] / lami[k - 1]
                        if do_cldice:
                            nitend[_idx2(i, k, pcols)] = (niic[_idx2(i, k, pcols)] * icldm[_idx2(i, k, pcols)] - ni[_idx2(i, k, pcols)]) / deltat
                    elif lami[k - 1] > lammaxi:
                        lami[k - 1] = lammaxi
                        n0i[k - 1] = pow_r8(lami[k - 1], di + 1.0) * dumi[_idx2(i, k, pcols)] / (ci * cons1)
                        niic[_idx2(i, k, pcols)] = n0i[k - 1] / lami[k - 1]
                        if do_cldice:
                            nitend[_idx2(i, k, pcols)] = (niic[_idx2(i, k, pcols)] * icldm[_idx2(i, k, pcols)] - ni[_idx2(i, k, pcols)]) / deltat
                    effi[_idx2(i, k, pcols)] = 1.5 / lami[k - 1] * 1000000.0
                else:
                    effi[_idx2(i, k, pcols)] = 25.0
                if not do_cldice:
                    effi[_idx2(i, k, pcols)] = re_ice[_idx2(i, k, pcols)] * 1000000.0
                if dumc[_idx2(i, k, pcols)] >= qsmall:
                    dumnc[_idx2(i, k, pcols)] = min(dumnc[_idx2(i, k, pcols)], dumc[_idx2(i, k, pcols)] * 1e+20)
                    if dumnc[_idx2(i, k, pcols)] < cdnl / rho[_idx2(i, k, pcols)]:
                        nctend[_idx2(i, k, pcols)] = (cdnl / rho[_idx2(i, k, pcols)] * lcldm[_idx2(i, k, pcols)] - nc[_idx2(i, k, pcols)]) / deltat
                    dumnc[_idx2(i, k, pcols)] = max(dumnc[_idx2(i, k, pcols)], cdnl / rho[_idx2(i, k, pcols)])
                    pgam[k - 1] = 0.0005714 * (ncic[_idx2(i, k, pcols)] / 1000000.0 * rho[_idx2(i, k, pcols)]) + 0.2714
                    pgam[k - 1] = 1.0 / pow_i(pgam[k - 1], 2) - 1.0
                    pgam[k - 1] = max(pgam[k - 1], 2.0)
                    pgam[k - 1] = min(pgam[k - 1], 15.0)
                    lamc[k - 1] = pow_r8(pi / 6.0 * rhow * dumnc[_idx2(i, k, pcols)] * gamma(pgam[k - 1] + 4.0) / (dumc[_idx2(i, k, pcols)] * gamma(pgam[k - 1] + 1.0)), 1.0 / 3.0)
                    lammin = (pgam[k - 1] + 1.0) / 5e-05
                    lammax = (pgam[k - 1] + 1.0) * omsm / 2e-06
                    if lamc[k - 1] < lammin:
                        lamc[k - 1] = lammin
                        ncic[_idx2(i, k, pcols)] = 6.0 * pow_i(lamc[k - 1], 3) * dumc[_idx2(i, k, pcols)] * gamma(pgam[k - 1] + 1.0) / (pi * rhow * gamma(pgam[k - 1] + 4.0))
                        nctend[_idx2(i, k, pcols)] = (ncic[_idx2(i, k, pcols)] * lcldm[_idx2(i, k, pcols)] - nc[_idx2(i, k, pcols)]) / deltat
                    elif lamc[k - 1] > lammax:
                        lamc[k - 1] = lammax
                        ncic[_idx2(i, k, pcols)] = 6.0 * pow_i(lamc[k - 1], 3) * dumc[_idx2(i, k, pcols)] * gamma(pgam[k - 1] + 1.0) / (pi * rhow * gamma(pgam[k - 1] + 4.0))
                        nctend[_idx2(i, k, pcols)] = (ncic[_idx2(i, k, pcols)] * lcldm[_idx2(i, k, pcols)] - nc[_idx2(i, k, pcols)]) / deltat
                    effc[_idx2(i, k, pcols)] = gamma(pgam[k - 1] + 4.0) / gamma(pgam[k - 1] + 3.0) / lamc[k - 1] / 2.0 * 1000000.0
                    lamcrad[_idx2(i, k, pcols)] = lamc[k - 1]
                    pgamrad[_idx2(i, k, pcols)] = pgam[k - 1]
                else:
                    effc[_idx2(i, k, pcols)] = 10.0
                    lamcrad[_idx2(i, k, pcols)] = 0.0
                    pgamrad[_idx2(i, k, pcols)] = 0.0
                if do_cldice:
                    deffi[_idx2(i, k, pcols)] = effi[_idx2(i, k, pcols)] * rhoi / 917.0 * 2.0
                else:
                    deffi[_idx2(i, k, pcols)] = effi[_idx2(i, k, pcols)] * 2.0
                dumnc[_idx2(i, k, pcols)] = 100000000.0
                if dumc[_idx2(i, k, pcols)] >= qsmall:
                    pgam[k - 1] = 0.0005714 * (ncic[_idx2(i, k, pcols)] / 1000000.0 * rho[_idx2(i, k, pcols)]) + 0.2714
                    pgam[k - 1] = 1.0 / pow_i(pgam[k - 1], 2) - 1.0
                    pgam[k - 1] = max(pgam[k - 1], 2.0)
                    pgam[k - 1] = min(pgam[k - 1], 15.0)
                    lamc[k - 1] = pow_r8(pi / 6.0 * rhow * dumnc[_idx2(i, k, pcols)] * gamma(pgam[k - 1] + 4.0) / (dumc[_idx2(i, k, pcols)] * gamma(pgam[k - 1] + 1.0)), 1.0 / 3.0)
                    lammin = (pgam[k - 1] + 1.0) / 5e-05
                    lammax = (pgam[k - 1] + 1.0) / 2e-06
                    if lamc[k - 1] < lammin:
                        lamc[k - 1] = lammin
                    elif lamc[k - 1] > lammax:
                        lamc[k - 1] = lammax
                    effc_fn[_idx2(i, k, pcols)] = gamma(pgam[k - 1] + 4.0) / gamma(pgam[k - 1] + 3.0) / lamc[k - 1] / 2.0 * 1000000.0
                else:
                    effc_fn[_idx2(i, k, pcols)] = 10.0
        for k in range(top_lev, pver + 1):
            if qc[_idx2(i, k, pcols)] + qctend[_idx2(i, k, pcols)] * deltat < qsmall:
                nctend[_idx2(i, k, pcols)] = -nc[_idx2(i, k, pcols)] / deltat
            if do_cldice and qi[_idx2(i, k, pcols)] + qitend[_idx2(i, k, pcols)] * deltat < qsmall:
                nitend[_idx2(i, k, pcols)] = -ni[_idx2(i, k, pcols)] / deltat
    for i in range(1, ncol + 1):
        for k in range(top_lev, pver + 1):
            if qsout[_idx2(i, k, pcols)] > 1e-07 and nsout[_idx2(i, k, pcols)] > 0.0:
                dsout[_idx2(i, k, pcols)] = 3.0 * rhosn / 917.0 * pow_r8(pi * rhosn * nsout[_idx2(i, k, pcols)] / qsout[_idx2(i, k, pcols)], -1.0 / 3.0)
    for i in range(1, ncol + 1):
        for k in range(top_lev, pver + 1):
            if qrout[_idx2(i, k, pcols)] > 1e-07 and nrout[_idx2(i, k, pcols)] > 0.0:
                reff_rain[_idx2(i, k, pcols)] = 1.5 * pow_r8(pi * rhow * nrout[_idx2(i, k, pcols)] / qrout[_idx2(i, k, pcols)], -1.0 / 3.0) * 1000000.0
            if qsout[_idx2(i, k, pcols)] > 1e-07 and nsout[_idx2(i, k, pcols)] > 0.0:
                reff_snow[_idx2(i, k, pcols)] = 1.5 * pow_r8(pi * rhosn * nsout[_idx2(i, k, pcols)] / qsout[_idx2(i, k, pcols)], -1.0 / 3.0) * 1000000.0
    for i in range(1, ncol + 1):
        for k in range(top_lev, pver + 1):
            if qc[_idx2(i, k, pcols)] + qctend[_idx2(i, k, pcols)] * deltat >= qsmall:
                dum = pow_i((qc[_idx2(i, k, pcols)] + qctend[_idx2(i, k, pcols)] * deltat) / lcldm[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)] * 1000.0, 2) / (0.109 * (nc[_idx2(i, k, pcols)] + nctend[_idx2(i, k, pcols)] * deltat) / lcldm[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)] / 1000000.0) * lcldm[_idx2(i, k, pcols)] / cldmax[_idx2(i, k, pcols)]
            else:
                dum = 0.0
            if qi[_idx2(i, k, pcols)] + qitend[_idx2(i, k, pcols)] * deltat >= qsmall:
                dum1 = pow_r8((qi[_idx2(i, k, pcols)] + qitend[_idx2(i, k, pcols)] * deltat) * rho[_idx2(i, k, pcols)] / icldm[_idx2(i, k, pcols)] * 1000.0 / 0.1, 1.0 / 0.63) * icldm[_idx2(i, k, pcols)] / cldmax[_idx2(i, k, pcols)]
            else:
                dum1 = 0.0
            if qsout[_idx2(i, k, pcols)] >= qsmall:
                dum1 = dum1 + pow_r8(qsout[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)] * 1000.0 / 0.1, 1.0 / 0.63)
            refl[_idx2(i, k, pcols)] = dum + dum1
            if rainrt[_idx2(i, k, pcols)] >= 0.001:
                dum = log10(pow_r8(rainrt[_idx2(i, k, pcols)], 6.0)) + 16.0
                dum = pow_r8(10.0, dum / 10.0)
            else:
                dum = 0.0
            refl[_idx2(i, k, pcols)] = refl[_idx2(i, k, pcols)] + dum
            areflz[_idx2(i, k, pcols)] = refl[_idx2(i, k, pcols)]
            if refl[_idx2(i, k, pcols)] > minrefl:
                refl[_idx2(i, k, pcols)] = 10.0 * log10(refl[_idx2(i, k, pcols)])
            else:
                refl[_idx2(i, k, pcols)] = -9999.0
            if refl[_idx2(i, k, pcols)] > mindbz:
                arefl[_idx2(i, k, pcols)] = refl[_idx2(i, k, pcols)]
                frefl[_idx2(i, k, pcols)] = 1.0
            else:
                arefl[_idx2(i, k, pcols)] = 0.0
                areflz[_idx2(i, k, pcols)] = 0.0
                frefl[_idx2(i, k, pcols)] = 0.0
            csrfl[_idx2(i, k, pcols)] = min(csmax, refl[_idx2(i, k, pcols)])
            if csrfl[_idx2(i, k, pcols)] > csmin:
                acsrfl[_idx2(i, k, pcols)] = refl[_idx2(i, k, pcols)]
                fcsrfl[_idx2(i, k, pcols)] = 1.0
            else:
                acsrfl[_idx2(i, k, pcols)] = 0.0
                fcsrfl[_idx2(i, k, pcols)] = 0.0
    for _s210 in range(1, pver + 1):
        for _s209 in range(1, pcols + 1):
            qrout2[_idx2(_s209, _s210, pcols)] = 0.0
    for _s212 in range(1, pver + 1):
        for _s211 in range(1, pcols + 1):
            qsout2[_idx2(_s211, _s212, pcols)] = 0.0
    for _s214 in range(1, pver + 1):
        for _s213 in range(1, pcols + 1):
            nrout2[_idx2(_s213, _s214, pcols)] = 0.0
    for _s216 in range(1, pver + 1):
        for _s215 in range(1, pcols + 1):
            nsout2[_idx2(_s215, _s216, pcols)] = 0.0
    for _s218 in range(1, pver + 1):
        for _s217 in range(1, pcols + 1):
            drout2[_idx2(_s217, _s218, pcols)] = 0.0
    for _s220 in range(1, pver + 1):
        for _s219 in range(1, pcols + 1):
            dsout2[_idx2(_s219, _s220, pcols)] = 0.0
    for _s222 in range(1, pver + 1):
        for _s221 in range(1, pcols + 1):
            freqs[_idx2(_s221, _s222, pcols)] = 0.0
    for _s224 in range(1, pver + 1):
        for _s223 in range(1, pcols + 1):
            freqr[_idx2(_s223, _s224, pcols)] = 0.0
    for i in range(1, ncol + 1):
        for k in range(top_lev, pver + 1):
            if qrout[_idx2(i, k, pcols)] > 1e-07 and nrout[_idx2(i, k, pcols)] > 0.0:
                qrout2[_idx2(i, k, pcols)] = qrout[_idx2(i, k, pcols)]
                nrout2[_idx2(i, k, pcols)] = nrout[_idx2(i, k, pcols)]
                drout2[_idx2(i, k, pcols)] = pow_r8(pi * rhow * nrout[_idx2(i, k, pcols)] / qrout[_idx2(i, k, pcols)], -1.0 / 3.0)
                freqr[_idx2(i, k, pcols)] = 1.0
            if qsout[_idx2(i, k, pcols)] > 1e-07 and nsout[_idx2(i, k, pcols)] > 0.0:
                qsout2[_idx2(i, k, pcols)] = qsout[_idx2(i, k, pcols)]
                nsout2[_idx2(i, k, pcols)] = nsout[_idx2(i, k, pcols)]
                dsout2[_idx2(i, k, pcols)] = pow_r8(pi * rhosn * nsout[_idx2(i, k, pcols)] / qsout[_idx2(i, k, pcols)], -1.0 / 3.0)
                freqs[_idx2(i, k, pcols)] = 1.0
    for i in range(1, ncol + 1):
        for k in range(top_lev, pver + 1):
            ncai[_idx2(i, k, pcols)] = dum2i[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)]
            ncal[_idx2(i, k, pcols)] = dum2l[_idx2(i, k, pcols)] * rho[_idx2(i, k, pcols)]
    for _s226 in range(1, pver + 1):
        for _s225 in range(1, pcols + 1):
            nfice[_idx2(_s225, _s226, pcols)] = 0.0
    for k in range(top_lev, pver + 1):
        for i in range(1, ncol + 1):
            dumc[_idx2(i, k, pcols)] = qc[_idx2(i, k, pcols)] + qctend[_idx2(i, k, pcols)] * deltat
            dumi[_idx2(i, k, pcols)] = qi[_idx2(i, k, pcols)] + qitend[_idx2(i, k, pcols)] * deltat
            dumfice = qsout[_idx2(i, k, pcols)] + qrout[_idx2(i, k, pcols)] + dumc[_idx2(i, k, pcols)] + dumi[_idx2(i, k, pcols)]
            if dumfice > qsmall and qsout[_idx2(i, k, pcols)] + dumi[_idx2(i, k, pcols)] > qsmall:
                nfice[_idx2(i, k, pcols)] = (qsout[_idx2(i, k, pcols)] + dumi[_idx2(i, k, pcols)]) / dumfice
            if nfice[_idx2(i, k, pcols)] > 1.0:
                nfice[_idx2(i, k, pcols)] = 1.0
