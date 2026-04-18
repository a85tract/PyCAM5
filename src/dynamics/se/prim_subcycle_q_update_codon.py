@inline
def _hy_idx(klev: int) -> int:
    """hvcoord%hyai / hvcoord%hybi declared as (nlev+1)."""
    return klev - 1


@inline
def _plane_idx(iidx: int, jidx: int, np: int) -> int:
    """ps_v and dp_np1 declared as (np,np)."""
    return (iidx - 1) + (jidx - 1) * np


@inline
def _q_idx(iidx: int, jidx: int, klev: int, qidx: int, np: int, nlev: int) -> int:
    """q and qdp slices declared as (np,np,nlev,qsize)."""
    return (
        (iidx - 1)
        + (jidx - 1) * np
        + (klev - 1) * np * np
        + (qidx - 1) * np * np * nlev
    )


@export
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
        hyai_delta = hyai[_hy_idx(k + 1)] - hyai[_hy_idx(k)]
        hybi_delta = hybi[_hy_idx(k + 1)] - hybi[_hy_idx(k)]

        for j in range(1, np + 1):
            for i in range(1, np + 1):
                plane_idx = _plane_idx(i, j, np)
                dp_np1[plane_idx] = hyai_delta * ps0 + hybi_delta * ps_v[plane_idx]

        for qidx in range(1, qsize + 1):
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    plane_idx = _plane_idx(i, j, np)
                    q_idx = _q_idx(i, j, k, qidx, np, nlev)
                    q_out[q_idx] = qdp[q_idx] / dp_np1[plane_idx]
