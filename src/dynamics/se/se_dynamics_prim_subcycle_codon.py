import se_dynamics_common_codon as _common

def prim_subcycle_dp3d_init_codon(
    np: int,
    nlev: int,
    ps0: float,
    hyai_p: cobj,
    hybi_p: cobj,
    ps_v_p: cobj,
    dp3d_p: cobj,
):
    hyai = Ptr[float](hyai_p)
    hybi = Ptr[float](hybi_p)
    ps_v = Ptr[float](ps_v_p)
    dp3d = Ptr[float](dp3d_p)

    for k in range(1, nlev + 1):
        hyai_delta = hyai[_common._hy_idx(k + 1)] - hyai[_common._hy_idx(k)]
        hybi_delta = hybi[_common._hy_idx(k + 1)] - hybi[_common._hy_idx(k)]

        for j in range(1, np + 1):
            for i in range(1, np + 1):
                plane_idx = _common._plane_idx(i, j, np)
                vol_idx = _common._vol_idx(i, j, k, np)
                dp3d[vol_idx] = hyai_delta * ps0 + hybi_delta * ps_v[plane_idx]


def prim_subcycle_q_update_codon(
    np: int,
    nlev: int,
    qsize: int,
    ps0: float,
    hyai_p: cobj,
    hybi_p: cobj,
    ps_v_p: cobj,
    dp_np1_p: cobj,
    qdp_p: cobj,
    q_p: cobj,
):
    hyai = Ptr[float](hyai_p)
    hybi = Ptr[float](hybi_p)
    ps_v = Ptr[float](ps_v_p)
    dp_np1 = Ptr[float](dp_np1_p)
    qdp = Ptr[float](qdp_p)
    q_out = Ptr[float](q_p)

    for k in range(1, nlev + 1):
        hyai_delta = hyai[_common._hy_idx(k + 1)] - hyai[_common._hy_idx(k)]
        hybi_delta = hybi[_common._hy_idx(k + 1)] - hybi[_common._hy_idx(k)]

        for j in range(1, np + 1):
            for i in range(1, np + 1):
                plane_idx = _common._plane_idx(i, j, np)
                dp_np1[plane_idx] = hyai_delta * ps0 + hybi_delta * ps_v[plane_idx]

        for qidx in range(1, qsize + 1):
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    plane_idx = _common._plane_idx(i, j, np)
                    q_idx = _common._q_idx(i, j, k, qidx, np, nlev)
                    q_out[q_idx] = qdp[q_idx] / dp_np1[plane_idx]


def qdp_time_avg_codon(
    np: int,
    nlev: int,
    qsize: int,
    rkstage: int,
    n0_qdp: int,
    np1_qdp: int,
    qdp_p: cobj,
):
    qdp = Ptr[float](qdp_p)
    rkstage_minus_one = rkstage - 1

    for qidx in range(1, qsize + 1):
        for k in range(1, nlev + 1):
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    n0_idx = _common._q_tl_idx(i, j, k, qidx, n0_qdp, np, nlev, qsize)
                    np1_idx = _common._q_tl_idx(i, j, k, qidx, np1_qdp, np, nlev, qsize)
                    qdp[np1_idx] = (qdp[n0_idx] + rkstage_minus_one * qdp[np1_idx]) / rkstage
