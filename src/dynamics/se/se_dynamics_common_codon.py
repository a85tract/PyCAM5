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


@inline
def _q_tl_idx(iidx: int, jidx: int, klev: int, qidx: int, tlidx: int, np: int, nlev: int, qsize: int) -> int:
    """state%Qdp declared as (np,np,nlev,qsize,2)."""
    return _q_idx(iidx, jidx, klev, qidx, np, nlev) + (tlidx - 1) * np * np * nlev * qsize


@inline
def _vol_idx(iidx: int, jidx: int, klev: int, np: int) -> int:
    """dp3d, dp, dp_star declared as (np,np,nlev)."""
    return (iidx - 1) + (jidx - 1) * np + (klev - 1) * np * np


@inline
def _vec2_idx(iidx: int, jidx: int, comp: int, np: int) -> int:
    """v and gv declared as (np,np,2)."""
    return (iidx - 1) + (jidx - 1) * np + (comp - 1) * np * np


@inline
def _vec3_idx(iidx: int, jidx: int, comp: int, np: int) -> int:
    """dum_cart declared as (np,np,3)."""
    return (iidx - 1) + (jidx - 1) * np + (comp - 1) * np * np


@inline
def _mat22_idx(iidx: int, jidx: int, row: int, col: int, np: int) -> int:
    """Dinv and D declared as (np,np,2,2)."""
    return (iidx - 1) + (jidx - 1) * np + (row - 1) * np * np + (col - 1) * np * np * 2


@inline
def _mat32_idx(iidx: int, jidx: int, row: int, col: int, np: int) -> int:
    """vec_sphere2cart declared as (np,np,3,2)."""
    return (iidx - 1) + (jidx - 1) * np + (row - 1) * np * np + (col - 1) * np * np * 3


@inline
def _field_vol_idx(iidx: int, jidx: int, klev: int, fidx: int, np: int, nlev: int) -> int:
    """ttmp declared as (np,np,nlev,2)."""
    return _vol_idx(iidx, jidx, klev, np) + (fidx - 1) * np * np * nlev


@inline
def _v_idx(iidx: int, jidx: int, comp: int, klev: int, np: int) -> int:
    """state%v slice declared as (np,np,2,nlev)."""
    return (iidx - 1) + (jidx - 1) * np + (comp - 1) * np * np + (klev - 1) * np * np * 2


@inline
def _ghost_col_idx(klev: int) -> int:
    """dpo declared as (-1:nlev+2)."""
    return klev + 1


@inline
def _ppm_grid_idx(row: int, jidx: int) -> int:
    """ppmdx declared as (10,0:nlev+1)."""
    return (row - 1) + jidx * 10


@inline
def _ppm_scratch_idx(jidx: int) -> int:
    """ppm_ai declared as (0:nlev), ppm_dma declared as (0:nlev+1)."""
    return jidx


@inline
def _ppm_coef_idx(comp: int, jidx: int) -> int:
    """coefs declared as (0:2,nlev)."""
    return comp + (jidx - 1) * 3


@inline
def _col_idx(klev: int) -> int:
    """pio declared as (nlev+2), pin declared as (nlev+1), z1/z2/kid declared as (nlev)."""
    return klev - 1


@inline
def _lev_q_idx(klev: int, qidx: int, nlev: int) -> int:
    """qmin and qmax declared as (nlev,qsize)."""
    return (klev - 1) + (qidx - 1) * nlev


@inline
def _cell_lev_idx(cell: int, klev: int, np: int) -> int:
    """ptens, dpmass, and workspaces declared as (np*np,nlev)."""
    ncols = np * np
    return (cell - 1) + (klev - 1) * ncols
