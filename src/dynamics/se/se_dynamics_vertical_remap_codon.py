import se_dynamics_common_codon as _common

def vertical_remap_rsplit_prepare_codon(
    np: int,
    nlev: int,
    ps0: float,
    hyai_p: cobj,
    hybi_p: cobj,
    ps_v_p: cobj,
    dp3d_p: cobj,
    dp_p: cobj,
    dp_star_p: cobj,
):
    hyai = Ptr[float](hyai_p)
    hybi = Ptr[float](hybi_p)
    ps_v = Ptr[float](ps_v_p)
    dp3d = Ptr[float](dp3d_p)
    dp = Ptr[float](dp_p)
    dp_star = Ptr[float](dp_star_p)

    for k in range(1, nlev + 1):
        hyai_delta = hyai[_common._hy_idx(k + 1)] - hyai[_common._hy_idx(k)]
        hybi_delta = hybi[_common._hy_idx(k + 1)] - hybi[_common._hy_idx(k)]

        for j in range(1, np + 1):
            for i in range(1, np + 1):
                plane_idx = _common._plane_idx(i, j, np)
                vol_idx = _common._vol_idx(i, j, k, np)
                dp[vol_idx] = hyai_delta * ps0 + hybi_delta * ps_v[plane_idx]
                dp_star[vol_idx] = dp3d[vol_idx]


def vertical_remap_t_scale_codon(
    np: int,
    nlev: int,
    t_p: cobj,
    dp_star_p: cobj,
    ttmp_p: cobj,
):
    t = Ptr[float](t_p)
    dp_star = Ptr[float](dp_star_p)
    ttmp = Ptr[float](ttmp_p)

    for k in range(1, nlev + 1):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                vol_idx = _common._vol_idx(i, j, k, np)
                ttmp[vol_idx] = t[vol_idx] * dp_star[vol_idx]


def vertical_remap_t_unscale_codon(
    np: int,
    nlev: int,
    ttmp_p: cobj,
    dp_p: cobj,
    t_p: cobj,
):
    ttmp = Ptr[float](ttmp_p)
    dp = Ptr[float](dp_p)
    t = Ptr[float](t_p)

    for k in range(1, nlev + 1):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                vol_idx = _common._vol_idx(i, j, k, np)
                t[vol_idx] = ttmp[vol_idx] / dp[vol_idx]


def vertical_remap_v_scale_codon(
    np: int,
    nlev: int,
    v_p: cobj,
    dp_star_p: cobj,
    ttmp_p: cobj,
):
    v = Ptr[float](v_p)
    dp_star = Ptr[float](dp_star_p)
    ttmp = Ptr[float](ttmp_p)

    for k in range(1, nlev + 1):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                vol_idx = _common._vol_idx(i, j, k, np)
                ttmp[_common._field_vol_idx(i, j, k, 1, np, nlev)] = v[_common._v_idx(i, j, 1, k, np)] * dp_star[vol_idx]
                ttmp[_common._field_vol_idx(i, j, k, 2, np, nlev)] = v[_common._v_idx(i, j, 2, k, np)] * dp_star[vol_idx]


def vertical_remap_v_unscale_codon(
    np: int,
    nlev: int,
    ttmp_p: cobj,
    dp_p: cobj,
    v_p: cobj,
):
    ttmp = Ptr[float](ttmp_p)
    dp = Ptr[float](dp_p)
    v = Ptr[float](v_p)

    for k in range(1, nlev + 1):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                vol_idx = _common._vol_idx(i, j, k, np)
                v[_common._v_idx(i, j, 1, k, np)] = ttmp[_common._field_vol_idx(i, j, k, 1, np, nlev)] / dp[vol_idx]
                v[_common._v_idx(i, j, 2, k, np)] = ttmp[_common._field_vol_idx(i, j, k, 2, np, nlev)] / dp[vol_idx]


def vertical_remap_ps_v_update_codon(
    np: int,
    nlev: int,
    hyai1: float,
    ps0: float,
    dp3d_p: cobj,
    ps_v_p: cobj,
):
    dp3d = Ptr[float](dp3d_p)
    ps_v = Ptr[float](ps_v_p)

    for j in range(1, np + 1):
        for i in range(1, np + 1):
            total = 0.0
            for k in range(1, nlev + 1):
                total += dp3d[_common._vol_idx(i, j, k, np)]
            ps_v[_common._plane_idx(i, j, np)] = hyai1 * ps0 + total
