#!/usr/bin/env python3
"""Compare one native/Codon CESM run pair.

The script compares CAM monthly history plus final restart-style outputs:

- cam.h0.*.nc
- cam.r.*.nc
- cam.rh0.*.nc
- cam.rs.*.nc

It compares variable data, not full NetCDF file blobs:

- numeric variables: fixed-format dump -> md5
- character variables: per-variable dump -> md5

It also extracts three timing metrics from CESM GPTL timing output:

- physpkg_st1
- bc_physics
- CPL:RUN_LOOP
"""

from __future__ import annotations

import argparse
import hashlib
import subprocess
import sys
from pathlib import Path

NUMERIC_TYPES = {
    "double",
    "float",
    "int",
    "short",
    "byte",
    "ubyte",
    "ushort",
    "uint",
    "int64",
    "uint64",
}
CHAR_TYPES = {"char", "string"}
INCLUDE_PREFIXES = ("h0.", "r.", "rh0.", "rs.")
DEFAULT_TIMERS = ("physpkg_st1", "bc_physics", "CPL:RUN_LOOP")


def run_text(cmd: list[str]) -> str:
    return subprocess.check_output(cmd, text=True)


def run_bytes(cmd: list[str]) -> bytes:
    return subprocess.check_output(cmd)


def rel_cam_key(path: Path) -> str:
    name = path.name
    marker = ".cam."
    idx = name.find(marker)
    if idx < 0:
        raise ValueError(f"unexpected CAM filename: {name}")
    return name[idx + len(marker) :]


def collect_cam_files(run_dir: Path) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for path in sorted(run_dir.glob("*.cam.*.nc")):
        key = rel_cam_key(path)
        if key.startswith(INCLUDE_PREFIXES):
            out[key] = path
    return out


def list_vars(nc_path: Path) -> tuple[list[str], list[str]]:
    text = run_text(["ncks", "-m", str(nc_path)])
    numeric: list[str] = []
    chars: list[str] = []
    for line in text.splitlines():
        parts = line.strip().split()
        if len(parts) < 2:
            continue
        dtype = parts[0]
        if dtype not in NUMERIC_TYPES and dtype not in CHAR_TYPES:
            continue
        name = parts[1].split("(", 1)[0]
        if dtype in NUMERIC_TYPES:
            numeric.append(name)
        else:
            chars.append(name)
    return numeric, chars


def digest_numeric(nc_path: Path, vars_: list[str]) -> str:
    payload = run_bytes(
        ["ncks", "-C", "-H", "-s", "%+.17g\n", "-v", ",".join(vars_), str(nc_path)]
    )
    return hashlib.md5(payload).hexdigest()


def digest_char_var(nc_path: Path, var: str) -> str:
    text = run_text(["ncdump", "-v", var, str(nc_path)])
    marker = "data:\n"
    idx = text.find(marker)
    if idx < 0:
        raise ValueError(f"failed to locate data section for {var} in {nc_path}")
    payload = text[idx + len(marker) :].encode()
    return hashlib.md5(payload).hexdigest()


def compare_file(native_path: Path, codon_path: Path) -> dict[str, object]:
    num_vars, char_vars = list_vars(native_path)
    native_num = digest_numeric(native_path, num_vars) if num_vars else None
    codon_num = digest_numeric(codon_path, num_vars) if num_vars else None
    char_diffs: list[str] = []
    for var in char_vars:
        if digest_char_var(native_path, var) != digest_char_var(codon_path, var):
            char_diffs.append(var)
    return {
        "numeric_count": len(num_vars),
        "numeric_equal": native_num == codon_num,
        "numeric_md5_native": native_num,
        "numeric_md5_codon": codon_num,
        "char_count": len(char_vars),
        "char_diff_count": len(char_diffs),
        "char_diff_vars": char_diffs,
    }


def case_root_from_run_dir(run_dir: Path) -> Path | None:
    caseroot_file = run_dir / "CASEROOT"
    if not caseroot_file.is_file():
        return None
    text = caseroot_file.read_text().strip()
    if not text:
        return None
    return Path(text)


def timing_dirs_for_run(run_dir: Path) -> list[Path]:
    dirs: list[Path] = []
    run_timing = run_dir / "timing"
    if run_timing.is_dir():
        dirs.append(run_timing)
    case_root = case_root_from_run_dir(run_dir)
    if case_root is not None:
        case_timing = case_root / "timing"
        if case_timing.is_dir():
            dirs.append(case_timing)
    return dirs


def find_timing_stats_file(run_dir: Path) -> Path:
    run_timing = run_dir / "timing" / "cesm_timing_stats"
    if run_timing.is_file():
        return run_timing

    candidates: list[Path] = []
    for timing_dir in timing_dirs_for_run(run_dir):
        candidates.extend(
            p
            for p in timing_dir.iterdir()
            if (
                p.name == "cesm_timing_stats"
                or (
                    p.name.startswith("cesm_timing_stats.")
                    and not p.name.endswith(".gz")
                )
            )
        )
    if not candidates:
        raise FileNotFoundError(
            "no cesm_timing_stats file found for "
            f"run_dir={run_dir}"
        )
    return max(candidates, key=lambda path: path.stat().st_mtime)


def extract_timer(run_dir: Path, timer_name: str) -> float:
    timing_file = find_timing_stats_file(run_dir)
    with timing_file.open() as fh:
        for line in fh:
            if f'"{timer_name}"' in line:
                parts = line.split()
                return float(parts[6])
    raise KeyError(f"{timer_name} not found in {timing_file}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native-run-dir", required=True, type=Path)
    parser.add_argument("--codon-run-dir", required=True, type=Path)
    parser.add_argument(
        "--timer",
        action="append",
        dest="timers",
        help="GPTL timer to compare; may be specified multiple times",
    )
    args = parser.parse_args()
    timers = tuple(args.timers) if args.timers else DEFAULT_TIMERS

    native_files = collect_cam_files(args.native_run_dir)
    codon_files = collect_cam_files(args.codon_run_dir)

    missing_native = sorted(set(codon_files) - set(native_files))
    missing_codon = sorted(set(native_files) - set(codon_files))
    if missing_native or missing_codon:
        print("file set mismatch", file=sys.stderr)
        if missing_native:
            print("missing in native:", missing_native, file=sys.stderr)
        if missing_codon:
            print("missing in codon:", missing_codon, file=sys.stderr)
        return 2

    overall_numeric_ok = True
    overall_char_ok = True
    for key in sorted(native_files):
        result = compare_file(native_files[key], codon_files[key])
        print(f"=== {key} ===")
        print(f"numeric_count={result['numeric_count']}")
        print(f"numeric_equal={result['numeric_equal']}")
        print(f"numeric_md5_native={result['numeric_md5_native']}")
        print(f"numeric_md5_codon={result['numeric_md5_codon']}")
        print(f"char_count={result['char_count']}")
        print(f"char_diff_count={result['char_diff_count']}")
        if result["char_diff_vars"]:
            print("char_diff_vars=" + ",".join(result["char_diff_vars"]))
        if not result["numeric_equal"]:
            overall_numeric_ok = False
        if result["char_diff_count"] != 0:
            overall_char_ok = False
        print()

    print("=== timing ===")
    for timer in timers:
        native_val = extract_timer(args.native_run_dir, timer)
        codon_val = extract_timer(args.codon_run_dir, timer)
        delta_pct = (codon_val / native_val - 1.0) * 100.0
        print(
            f"{timer}: native={native_val:.6f} codon={codon_val:.6f} "
            f"delta_pct={delta_pct:.3f}"
        )

    print()
    print(f"overall_numeric_equal={overall_numeric_ok}")
    print(f"overall_char_equal={overall_char_ok}")
    return 0 if overall_numeric_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
