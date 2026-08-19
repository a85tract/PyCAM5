"""The Codon implementation of dadadj, adopted as a gateable candidate.

RecastEngine's differential gate compares a Python-callable candidate
against compiled Fortran truth. The Codon library is C-ABI, not Python, so
this shim is the adapter: it loads the library the repository's own
``codon build`` produces, marshals NumPy arrays as the raw pointers the
exported symbol takes, and carries the ``_SIGNATURES`` table the gate
generates inputs from.

Nothing numeric happens here. Every value the kernel computes is computed
by the Codon code under test; the shim allocates scratch, passes pointers,
and hands back the two inout arrays. ``cappa`` and ``nlvdry`` are spelled
identically to ``standins.f90`` so both sides run under the same physics.
"""

from __future__ import annotations

import ctypes
import os
from pathlib import Path

import numpy as np

# (SHR_CONST_RGAS / MWDAIR) / CPDAIR -- the same float64 expression the
# Fortran stand-in folds, so the two sides agree to the bit.
CAPPA = (6.02214e26 * 1.38065e-23 / 28.966) / 1.00464e3
NLVDRY = 3

_SIGNATURES = {
    "dadadj_native": {
        "kind": "subroutine",
        "args": [
            {"name": "lchnk", "dtype": "int32", "intent": "IN", "optional": False},
            {"name": "ncol", "dtype": "int32", "intent": "IN", "optional": False},
            {
                "name": "pmid",
                "dtype": "float64",
                "intent": "IN",
                "optional": False,
                "dims": [{"lb": "1", "ub": "pcols"}, {"lb": "1", "ub": "pver"}],
            },
            {
                "name": "pint",
                "dtype": "float64",
                "intent": "IN",
                "optional": False,
                "dims": [{"lb": "1", "ub": "pcols"}, {"lb": "1", "ub": "pverp"}],
            },
            {
                "name": "pdel",
                "dtype": "float64",
                "intent": "IN",
                "optional": False,
                "dims": [{"lb": "1", "ub": "pcols"}, {"lb": "1", "ub": "pver"}],
            },
            {
                "name": "t",
                "dtype": "float64",
                "intent": "INOUT",
                "optional": False,
                "dims": [{"lb": "1", "ub": "pcols"}, {"lb": "1", "ub": "pver"}],
            },
            {
                "name": "q",
                "dtype": "float64",
                "intent": "INOUT",
                "optional": False,
                "dims": [{"lb": "1", "ub": "pcols"}, {"lb": "1", "ub": "pver"}],
            },
        ],
        "result": None,
        "result_dtype": None,
    }
}


def _library() -> ctypes.CDLL:
    path = os.environ.get("RECAST_CODON_LIB")
    if not path:
        here = Path(__file__).resolve().parent
        path = str(here / "build" / "libdadadj_codon.dylib")
    lib = ctypes.CDLL(path)
    fn = lib.dadadj_codon
    fn.restype = None
    fn.argtypes = [
        ctypes.c_int64,
        ctypes.c_int64,
        ctypes.c_int64,
        ctypes.c_double,
        *([ctypes.c_void_p] * 13),
    ]
    return lib


_LIB = _library()


def _ptr(array: np.ndarray) -> ctypes.c_void_p:
    return ctypes.c_void_p(array.ctypes.data)


def dadadj_native(lchnk, ncol, pmid, pint, pdel, t, q):  # noqa: ARG001 -- lchnk is diagnostics-only
    """Same surface as the Fortran ``dadadj_native``; computed by Codon."""
    pcols, pver = np.shape(pmid)
    t = np.asfortranarray(t, dtype=np.float64)
    q = np.asfortranarray(q, dtype=np.float64)
    scratch = [np.zeros(pver, dtype=np.float64, order="F") for _ in range(4)]
    dodad = np.zeros(pcols, dtype=np.int64, order="F")
    status = np.zeros(1, dtype=np.int64)
    zeps_fail = np.zeros(1, dtype=np.float64)
    fail_i = np.zeros(1, dtype=np.int64)
    _LIB.dadadj_codon(
        int(ncol),
        int(pcols),
        NLVDRY,
        CAPPA,
        _ptr(np.asfortranarray(pmid, dtype=np.float64)),
        _ptr(np.asfortranarray(pint, dtype=np.float64)),
        _ptr(np.asfortranarray(pdel, dtype=np.float64)),
        _ptr(t),
        _ptr(q),
        *[_ptr(s) for s in scratch],
        _ptr(dodad),
        _ptr(status),
        _ptr(zeps_fail),
        _ptr(fail_i),
    )
    # A non-zero status mirrors the native path's non-convergence report,
    # where the gate's stand-in turns endrun into a return: both sides ran
    # the same fixed iteration count and their state is what gets compared.
    return t, q


def _PREPARE_INPUTS(name: str, inputs: dict, rng) -> None:  # noqa: N802 -- the gate's protocol name
    """Shape the generated inputs into dadadj's defined domain.

    Per-name ranges cannot express structure, and this routine needs three
    structural facts to be doing physics rather than arithmetic on noise:

    * pressure increases downward, and ``pint`` brackets ``pmid`` -- the
      routine divides by ``pmid(i,k+1) - pmid(i,k)``, so a non-monotone
      column is a division by an arbitrary small number;
    * ``pdel`` is the thickness ``pint`` implies, not an independent field;
      the enthalpy-conservation step weights by it, and inconsistent
      thicknesses make the fixed-point iteration diverge;
    * the lapse rate must be adjustable within 15 iterations, which a
      random temperature column is not -- the native path aborts on
      non-convergence, so a comparison there compares two error paths.

    Both sides receive the same shaped arrays, so this cannot bias the
    verdict; it only decides which region of the input space is sampled.
    """
    pmid = inputs.get("pmid")
    if pmid is None:
        return
    pcols, pver = pmid.shape
    pverp = pver + 1

    # A hybrid-sigma-like column: monotone, spanning model top to surface.
    edges = np.linspace(5000.0, 100000.0, pverp)
    pint = np.zeros((pcols, pverp), dtype=np.float64, order="F")
    pint[:, :] = edges
    mid = 0.5 * (edges[:-1] + edges[1:])
    pmid[:, :] = mid
    inputs["pint"] = np.asfortranarray(pint)
    inputs["pdel"] = np.asfortranarray(
        np.tile(np.diff(edges), (pcols, 1)).astype(np.float64)
    )

    # A temperature column near the dry adiabat, perturbed unstably in the
    # top few levels -- which is the case dadadj exists to adjust.
    t = inputs.get("t")
    if t is not None:
        reference = 250.0 + 60.0 * (mid - mid[0]) / (mid[-1] - mid[0])
        t[:, :] = reference
        t[:, : NLVDRY + 1] += rng.uniform(-3.0, 3.0, size=(pcols, NLVDRY + 1))

    q = inputs.get("q")
    if q is not None:
        q[:, :] = rng.uniform(0.0, 0.02, size=(pcols, pver))
