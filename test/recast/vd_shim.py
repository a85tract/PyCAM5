"""The Codon vertical_diffusion_ptend_core, as a gateable candidate.

Same adapter role as ``dadadj_shim.py``: loads the library the repository's
own ``codon build`` produces, marshals NumPy arrays as the raw pointers the
exported symbol takes, and carries the signature table the gate generates
inputs from. Nothing numeric happens here.

The routine is pure -- no module state, no calls -- so there is nothing to
initialize and no way for the two sides to disagree about state neither was
given. That is why it was chosen as the first three-way comparison.
"""

from __future__ import annotations

import ctypes
import os
from pathlib import Path

import numpy as np

PCOLS, PVER, PCNST = 8, 30, 4


def _dims(*names: str) -> list[dict[str, str]]:
    return [{"lb": "1", "ub": n} for n in names]


def _arg(name: str, intent: str, *dims: str) -> dict:
    entry = {"name": name, "dtype": "float64", "intent": intent, "optional": False}
    if dims:
        entry["dims"] = _dims(*dims)
    return entry


_SIGNATURES = {
    "vertical_diffusion_ptend_core_native": {
        "kind": "subroutine",
        "args": [
            {"name": "ncol", "dtype": "int32", "intent": "IN", "optional": False},
            {"name": "psetcols_local", "dtype": "int32", "intent": "IN", "optional": False},
            _arg("q_tmp_local", "IN", "pcols", "pver", "pcnst"),
            _arg("s_tmp_local", "IN", "pcols", "pver"),
            _arg("u_tmp_local", "IN", "pcols", "pver"),
            _arg("v_tmp_local", "IN", "pcols", "pver"),
            _arg("state_q_local", "IN", "pcols", "pver", "pcnst"),
            _arg("state_s_local", "IN", "pcols", "pver"),
            _arg("state_u_local", "IN", "pcols", "pver"),
            _arg("state_v_local", "IN", "pcols", "pver"),
            _arg("sl_local", "IN", "pcols", "pver"),
            _arg("qt_local", "IN", "pcols", "pver"),
            _arg("sl_prepbl_local", "IN", "pcols", "pver"),
            _arg("qt_prepbl_local", "IN", "pcols", "pver"),
            _arg("rztodt_local", "IN"),
            _arg("ptend_q_local", "INOUT", "psetcols_local", "pver", "pcnst"),
            _arg("ptend_s_local", "INOUT", "psetcols_local", "pver"),
            _arg("ptend_u_local", "INOUT", "psetcols_local", "pver"),
            _arg("ptend_v_local", "INOUT", "psetcols_local", "pver"),
            _arg("slten_local", "INOUT", "pcols", "pver"),
            _arg("qtten_local", "INOUT", "pcols", "pver"),
        ],
        "result": None,
        "result_dtype": None,
    }
}


def _library() -> ctypes.CDLL:
    path = os.environ.get("RECAST_CODON_LIB") or str(
        Path(__file__).resolve().parent / "build" / "libvertical_diffusion_codon.dylib"
    )
    lib = ctypes.CDLL(path)
    fn = lib.vertical_diffusion_ptend_core_codon
    fn.restype = None
    fn.argtypes = [
        ctypes.c_int64,  # ncol
        ctypes.c_int64,  # pcols
        ctypes.c_int64,  # pver
        ctypes.c_int64,  # pcnst
        ctypes.c_int64,  # psetcols
        ctypes.c_double,  # rztodt
        *([ctypes.c_void_p] * 18),
    ]
    return lib


_LIB = _library()


def _p(a: np.ndarray) -> ctypes.c_void_p:
    return ctypes.c_void_p(a.ctypes.data)


def vertical_diffusion_ptend_core_native(
    ncol,
    psetcols_local,
    q_tmp_local,
    s_tmp_local,
    u_tmp_local,
    v_tmp_local,
    state_q_local,
    state_s_local,
    state_u_local,
    state_v_local,
    sl_local,
    qt_local,
    sl_prepbl_local,
    qt_prepbl_local,
    rztodt_local,
    ptend_q_local,
    ptend_s_local,
    ptend_u_local,
    ptend_v_local,
    slten_local,
    qtten_local,
):
    """Same surface as the Fortran routine; computed by Codon."""
    outs = [
        np.asfortranarray(a, dtype=np.float64)
        for a in (
            ptend_q_local,
            ptend_s_local,
            ptend_u_local,
            ptend_v_local,
            slten_local,
            qtten_local,
        )
    ]
    ins = [
        np.asfortranarray(a, dtype=np.float64)
        for a in (
            q_tmp_local,
            s_tmp_local,
            u_tmp_local,
            v_tmp_local,
            state_q_local,
            state_s_local,
            state_u_local,
            state_v_local,
            sl_local,
            qt_local,
            sl_prepbl_local,
            qt_prepbl_local,
        )
    ]
    _LIB.vertical_diffusion_ptend_core_codon(
        int(ncol),
        PCOLS,
        PVER,
        PCNST,
        int(psetcols_local),
        float(rztodt_local),
        *[_p(a) for a in ins],
        *[_p(a) for a in outs],
    )
    return tuple(outs)
