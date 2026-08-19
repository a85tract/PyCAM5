#!/usr/bin/env python3
"""Module-level gate: the Codon dadadj against the native Fortran, to the bit.

PyCAM5's standing validation is whole-model: six-month run pairs comparing
`overall_numeric_equal=True`. This is the complementary check at the other
end of the cost scale -- one routine, seconds, runnable on a laptop, driven
by RecastEngine's verification stages:

    frontend (cesm)         analyze dadadj.F90 as it sits in the tree
    candidate (adopted)     the Codon library behind a ctypes shim
    oracle (f2py-golden)    dadadj_native compiled by gfortran, stand-ins
                            supplying the framework modules
    gate (differential.bitexact)
                            same generated inputs to both, every output
                            compared bit for bit, evidence manifest written

Usage:
    codon build --relocation-model=pic -release --lib \
        -o test/recast/build/libdadadj_codon.dylib src/physics/cam/dadadj_codon.py
    python test/recast/run_gate.py
"""

from __future__ import annotations

import sys
from pathlib import Path

from recast.executors.local import LocalExecutor
from recast.model import Candidate, Unit
from recast.oracle.f2py import F2pyGoldenOracle
from recast.store.filesystem import FilesystemEvidenceStore
from recast.run import _evidence  # the runner's evidence assembly
from recast.verify.bitexact import BitexactVerifier
from recast_cesm.frontend import factory as cesm_frontend

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
CAM = REPO / "src" / "physics" / "cam"

DIMS = {"pcols": 8, "pver": 30, "pverp": 31}
RANGES = {
    "t": [250.0, 320.0],
    "q": [0.0, 0.02],
    "pmid": [5000.0, 100000.0],
    "pint": [5000.0, 100000.0],
    "pdel": [10.0, 5000.0],
}


def main() -> int:
    import os

    os.environ["RECAST_CODON_LIB"] = str(HERE / "build" / "libdadadj_codon.dylib")
    workspace = HERE / "work"
    workspace.mkdir(exist_ok=True)
    executor = LocalExecutor()

    frontend = cesm_frontend()
    unit = Unit(uid="fortran:dadadj", kind="file", sources=(Path("dadadj.F90"),))
    facts = frontend.analyze(unit, CAM)

    candidate = Candidate(
        unit=unit.uid,
        transform="pycam5.adopt.codon",
        files={Path("dadadj_native_numpy.py"): (HERE / "dadadj_shim.py").read_bytes()},
        notes={"adopted": "Codon library built from src/physics/cam/dadadj_codon.py"},
    )

    config = {
        "root": CAM,
        "subprograms": ["dadadj_native"],
        "extra_sources": [str(HERE / "standins.f90")],
        "wrapper_parameters": DIMS,
        "dims": DIMS,
        "ranges": RANGES,
        "trials": 20,
    }
    oracle = F2pyGoldenOracle()
    ref = oracle.materialize(unit, facts, workspace, executor, config)

    verdict = BitexactVerifier().verify(unit, candidate, ref, workspace, executor, config)
    print(f"differential.bitexact: {verdict.confidence.value}")
    print(f"detail: {verdict.detail}")
    for name, metrics in verdict.metrics.get("subprograms", {}).items():
        print(f"  {name}: {metrics}")

    class _Recipe:
        name = "pycam5-module-gate"

    class _Run:
        candidate = None
        oracle = ref

    run_view = _Run()
    run_view.candidate = candidate
    store = FilesystemEvidenceStore(root=HERE / "evidence")
    uri = store.put(_evidence(_Recipe(), "local", unit, run_view, verdict))
    print(f"evidence: {uri}")
    write_summary(
        [
            {
                "routine": "dadadj_native",
                "comparison": "codon vs native Fortran",
                "verifier": verdict.verifier,
                "confidence": verdict.confidence.value,
                "passed": verdict.passed,
                "points": verdict.metrics.get("points"),
                "bit_exact": verdict.metrics.get("bit_exact"),
                "max_ulp": verdict.metrics.get("max_ulp"),
                "oracle": ref.key,
            }
        ]
    )
    return 0 if verdict.passed else 1


def write_summary(entries: list[dict]) -> None:
    """Record what these gates concluded, as a committed file.

    The manifests under ``evidence/`` are the audit trail -- one immutable
    document per verdict per run, and one file per attempt including the
    attempts that failed, which is why they are not committed. This is the
    other record: current state, one entry per routine and verifier,
    regenerated rather than appended, with no wall-clock time or paths in it
    so that two runs over the same revisions produce the same bytes.

    Committing it is what lets a reader -- or the CESM dashboard, which counts
    these -- see what has been verified without owning a Fortran compiler and
    a Codon toolchain.
    """
    import json
    import subprocess

    def revision(path: Path) -> str:
        out = subprocess.run(  # noqa: S603 -- git, in this repository
            ["git", "rev-parse", "HEAD"], cwd=path, capture_output=True, text=True, check=False
        )
        return out.stdout.strip()[:7] if out.returncode == 0 else "unknown"

    target = HERE / "verification.json"
    # Merge rather than overwrite: this file is what the whole harness has
    # concluded, not what the last driver to run concluded. Each (routine,
    # comparison) pair has one entry, and re-running a driver replaces its own
    # entries while leaving the other drivers' alone.
    known: dict[tuple[str, str], dict] = {}
    if target.is_file():
        try:
            for entry in json.loads(target.read_text())["routines"]:
                known[(entry["routine"], entry["comparison"])] = entry
        except (OSError, json.JSONDecodeError, KeyError, TypeError):
            known = {}
    for entry in entries:
        known[(entry["routine"], entry["comparison"])] = entry

    payload = {
        "schema": 1,
        "harness": "test/recast",
        "source_revision": revision(REPO),
        "routines": [known[key] for key in sorted(known)],
    }
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"summary: {target} ({len(payload['routines'])} comparison(s) recorded)")


if __name__ == "__main__":
    sys.exit(main())
