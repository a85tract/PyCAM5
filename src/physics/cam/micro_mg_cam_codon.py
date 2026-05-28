from math import gamma, pi


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
def micro_mg_cam_register_plan_codon(
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
def micro_mg_cam_p1_codon(n: int) -> int:
    if n >= 0:
        return 1
    return 0


@export
def micro_mg_cam_p2_codon(n1: int, n2: int) -> int:
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
def new_mgpacker_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def mgpacker_finalize_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


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
def mgfieldpostproc_finalize_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


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
def new_mgpostproc_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def mgpostproc_finalize_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def add_field_1d_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def add_field_2d_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


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
def micro_mg1_0_init_scalars_codon(
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
