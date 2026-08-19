#!/usr/bin/env python3
"""Three-way bit-exactness: engine translation, Codon port, native Fortran.

The Codon port's agreement with the native Fortran is already established by
this repository's long-run validation. This asks the other two questions:

    engine translation  vs  native Fortran   -- does the rule-driven
                                                translator reproduce CAM?
    engine translation  vs  Codon port       -- do the two independent
                                                modernizations agree?

All three sides derive from the same routine text.
``vertical_diffusion_ptend_core_native`` is extracted byte-for-byte from
``src/physics/cam/vertical_diffusion.F90`` into ``vd_ptend_core.f90``,
because that file use-imports seventeen framework modules and the routine
needs none of them. This script re-derives the extraction and refuses to run
if the copy has drifted, so the isolation cannot become a rewrite.

The routine is pure -- no module state, no calls -- which is why it was
chosen first: there is nothing to initialize, and no way for the three sides
to disagree about state none of them was given.

Usage:
    clang -c -fPIC -O1 -o test/recast/build/vd_callback_stubs.o \
        test/recast/vd_callback_stubs.c
    codon build --relocation-model=pic -release --lib \
        --linker-flags="$PWD/test/recast/build/vd_callback_stubs.o" \
        -o test/recast/build/libvertical_diffusion_codon.dylib \
        src/physics/cam/vertical_diffusion_codon.py
    python test/recast/run_threeway.py
"""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
CAM = REPO / "src" / "physics" / "cam"
NAME = "vertical_diffusion_ptend_core_native"
DIMS = {"pcols": 8, "pver": 30, "pcnst": 4, "psetcols_local": 8, "ncol": 8}
RANGES = {
    "s_tmp_local": [2.5e5, 3.5e5],
    "state_s_local": [2.5e5, 3.5e5],
    "u_tmp_local": [-60.0, 60.0],
    "state_u_local": [-60.0, 60.0],
    "v_tmp_local": [-60.0, 60.0],
    "state_v_local": [-60.0, 60.0],
    "q_tmp_local": [0.0, 0.03],
    "state_q_local": [0.0, 0.03],
    "sl_local": [2.5e5, 3.5e5],
    "sl_prepbl_local": [2.5e5, 3.5e5],
    "qt_local": [0.0, 0.03],
    "qt_prepbl_local": [0.0, 0.03],
    "rztodt_local": [1.0 / 3600.0, 1.0 / 900.0],
    "ptend_q_local": [0.0, 0.0],
    "ptend_s_local": [0.0, 0.0],
    "ptend_u_local": [0.0, 0.0],
    "ptend_v_local": [0.0, 0.0],
    "slten_local": [0.0, 0.0],
    "qtten_local": [0.0, 0.0],
}


def check_extraction() -> None:
    """The extracted copy must still be the routine it claims to be."""
    original = (CAM / "vertical_diffusion.F90").read_text().splitlines(keepends=True)
    start = next(i for i, l in enumerate(original) if re.match(rf"\s*subroutine\s+{NAME}\s*\(", l))
    end = next(i for i, l in enumerate(original) if re.match(rf"\s*end\s+subroutine\s+{NAME}\b", l))
    body = "".join(original[start : end + 1])
    if body not in (HERE / "vd_ptend_core.f90").read_text():
        raise SystemExit(
            f"vd_ptend_core.f90 no longer contains {NAME} verbatim: re-extract it "
            "rather than comparing against a copy that has drifted"
        )
    print(f"extraction verified: {NAME} byte-identical to vertical_diffusion.F90:{start + 1}")


def main() -> int:
    os.environ["RECAST_CODON_LIB"] = str(HERE / "build" / "libvertical_diffusion_codon.dylib")
    check_extraction()

    from recast.executors.local import LocalExecutor
    from recast.model import Candidate, Confidence, Unit
    from recast.oracle.f2py import F2pyGoldenOracle
    from recast.run import _evidence
    from recast.store.filesystem import FilesystemEvidenceStore
    from recast.verify.bitexact import BitexactVerifier
    from recast_cesm.frontend import factory as cesm_frontend
    from recast_cesm.translate import CamTranslation

    workspace = HERE / "work-threeway"
    workspace.mkdir(exist_ok=True)
    executor = LocalExecutor()

    # --- the engine's own translation, from the extracted routine -----------
    frontend = cesm_frontend()
    unit = next(
        u
        for u in frontend.discover(HERE)
        if u.parent is None and u.sources and str(u.sources[0]) == "vd_ptend_core.f90"
    )
    facts = frontend.analyze(unit, HERE)
    engine = CamTranslation().apply(unit, facts, {"root": HERE})
    print(f"engine translation: {sorted(str(p) for p in engine.files)}, "
          f"{len(engine.deferred)} deferred")
    if engine.deferred:
        for entry in engine.deferred:
            print("  deferred:", entry[:120])

    # --- the Codon port, behind the shim ------------------------------------
    codon = Candidate(
        unit=unit.uid,
        transform="pycam5.adopt.codon",
        files={Path("vd_ptend_core_numpy.py"): (HERE / "vd_shim.py").read_bytes()},
        notes={"adopted": "Codon library from src/physics/cam/vertical_diffusion_codon.py"},
    )

    # --- the native Fortran, compiled -------------------------------------
    config = {
        "root": HERE,
        "subprograms": [NAME],
        # The declared shapes name pcols/pver/pcnst, which the wrapper must
        # be able to fold; psetcols_local is a dummy so it needs no entry.
        "wrapper_parameters": {k: v for k, v in DIMS.items() if k in ("pcols", "pver", "pcnst")},
        "dims": DIMS,
        "ranges": RANGES,
        "trials": 20,
    }
    oracle = F2pyGoldenOracle()
    reference = oracle.materialize(unit, facts, workspace, executor, config)
    print(f"oracle: {reference.key}")

    store = FilesystemEvidenceStore(root=HERE / "evidence")

    class _Recipe:
        name = "pycam5-threeway"

    verdicts = {}
    for label, candidate in (("engine", engine), ("codon", codon)):
        space = workspace / label
        space.mkdir(exist_ok=True)
        verdict = BitexactVerifier().verify(unit, candidate, reference, space, executor, config)
        verdicts[label] = verdict
        metrics = verdict.metrics.get("subprograms", {}).get(NAME, {})
        print(f"\n{label:>7} vs native Fortran: {verdict.confidence.value}")
        print(f"         {verdict.detail}")
        if metrics:
            print(f"         {metrics}")

        class _Run:
            pass

        view = _Run()
        view.candidate = candidate
        view.oracle = reference
        uri = store.put(_evidence(_Recipe(), "local", unit, view, verdict))
        print(f"         evidence: {uri}")

    # --- engine vs codon, directly -----------------------------------------
    print("\n engine vs codon (direct):")
    agree = _compare_candidates(engine, codon, workspace, config)
    print(f"         {agree}")

    both = all(v.confidence is Confidence.BIT_EXACT for v in verdicts.values())
    return 0 if both and agree.startswith("bit-identical") else 1


def _compare_candidates(engine, codon, workspace: Path, config: dict) -> str:
    """Call both modernizations on the same inputs and compare their outputs.

    The two verdicts above already imply this when both are BIT_EXACT, but
    stating it directly is cheap and is the question actually asked: do the
    two independent modernizations of this routine agree with each other?
    """
    import numpy as np

    from recast.verify.bitexact import BitexactVerifier
    from recast.verify.ulp import ulp_audit

    verifier = BitexactVerifier()
    engine_module = verifier._load_candidate(engine, workspace / "engine")
    codon_module = verifier._load_candidate(codon, workspace / "codon")
    signature = engine_module._SIGNATURES[NAME]

    audits = []
    for trial in range(int(config["trials"])):
        rng = np.random.default_rng(
            int.from_bytes(f"{NAME}:{trial}".encode(), "big") % 2**32
        )
        inputs = {}
        for argument in signature["args"]:
            if argument["intent"] == "OUT":
                continue
            inputs[argument["name"]] = verifier._value(
                np, argument, config["dims"], {k.lower(): tuple(v) for k, v in config["ranges"].items()}, rng
            )
        for name in ("ncol", "psetcols_local"):
            inputs[name] = np.int32(config["dims"][name])
        left = engine_module.__dict__[NAME](**{k: _copy(np, v) for k, v in inputs.items()})
        right = codon_module.__dict__[NAME](**{k: _copy(np, v) for k, v in inputs.items()})
        for a, b in zip(left, right, strict=True):
            audits.append(ulp_audit(np.ravel(a).tolist(), np.ravel(b).tolist()))
    points = sum(a["total_points"] for a in audits)
    exact = sum(a["bit_exact"] for a in audits)
    worst = max(a["max_ulp"] for a in audits)
    if exact == points:
        return f"bit-identical on {points} points across {len(audits)} output arrays"
    return f"{points - exact}/{points} points differ, max {worst} ULP"


def _copy(np, value):
    return np.array(value, copy=True, order="F") if hasattr(value, "shape") else value


if __name__ == "__main__":
    sys.exit(main())
