#!/usr/bin/env python3
"""Track the 36-entry PI-atm standalone-physics campaign.

This phase-2 ledger is intentionally independent from the completed 69-entry
decoupling tracker.  JSON is the only editable source of truth; Markdown is a
generated view.  A process can become BFB only through ``record-run`` with
standalone dependency/build evidence and per-entry execution proof.
"""

from __future__ import annotations

import argparse
import copy
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import tempfile


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_STATUS_FILE = REPO_ROOT / "src/atmos_phys/standalone/standalone_status.json"
DEFAULT_PROGRESS_FILE = REPO_ROOT / "src/atmos_phys/standalone/STANDALONE_PROGRESS.md"

ALLOWED_STATUSES = (
    "planned",
    "in_progress",
    "dependency_pass",
    "build_pass",
    "50step_running",
    "bfb",
    "failed",
)
COMPLETED_STATUSES = ("dependency_pass", "build_pass", "50step_running", "bfb")
EXPECTED_BATCHES = {
    "S01": (
        "compute_vdiff_run",
        "gw_prof_run",
        "gw_oro_src_run",
        "gw_drag_prof_run",
        "energy_change_run",
        "ice_macro_tend_run",
        "rrtmg_sw_run",
        "rrtmg_lw_run",
        "modal_aero_depvel_part_run",
    ),
    "S02": (
        "compute_eddy_diff_run",
        "dadadj_run",
        "zm_convr_run",
        "zm_conv_evap_run",
        "zm_conv_momtran_run",
        "zm_conv_convtran_run",
        "compute_uwshcu_inv_run",
        "compute_uwshcu_run",
        "dust_sediment_tend_run",
    ),
    "S03": (
        "gas_phase_chemdr_run",
        "aero_model_gasaerexch_run",
        "modal_aero_gasaerexch_sub_run",
        "modal_aero_newnuc_sub_run",
        "modal_aero_coag_sub_run",
        "aero_model_drydep_run",
        "modal_aero_calcsize_sub_run",
        "modal_aero_wateruptake_dr_run",
        "wetdepa_v2_run",
    ),
    "S04": (
        "neu_wetdep_tend_run",
        "aero_model_wetdep_run",
        "nucleate_ice_cam_calc_run",
        "dropmixnuc_run",
        "cldfrc_run",
        "mmacro_pcond_run",
        "micro_mg_tend_run",
        "rad_rrtmg_sw_run",
        "rad_rrtmg_lw_run",
    ),
}
EXPECTED_DEPENDENCIES = {
    "S01": (),
    "S02": ("S01",),
    "S03": ("S01",),
    "S04": ("S01", "S03"),
}
REQUIRED_BFB_ARTIFACTS = (
    "dependency_report",
    "standalone_build_report",
    "executable",
    "run_environment",
    "filepath",
    "namelist",
    "compare_report",
)


class StatusError(ValueError):
    """Raised when ledger data or evidence is invalid."""


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z"
    )


def load_status(path: Path) -> dict:
    try:
        with path.open("r", encoding="utf-8") as handle:
            data = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise StatusError(f"cannot read status file {path}: {exc}") from exc
    if not isinstance(data, dict):
        raise StatusError("status root must be a JSON object")
    return data


def write_json_atomic(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(data, indent=2, ensure_ascii=False) + "\n"
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=path.parent, prefix=f".{path.name}.", delete=False
    ) as handle:
        handle.write(payload)
        temporary = handle.name
    os.replace(temporary, path)


def write_text_atomic(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=path.parent, prefix=f".{path.name}.", delete=False
    ) as handle:
        handle.write(value)
        temporary = handle.name
    os.replace(temporary, path)


def process_index(data: dict) -> dict[str, dict]:
    return {item["entry_point"]: item for item in data.get("processes", [])}


def run_index(data: dict) -> dict[str, dict]:
    return {item["run_id"]: item for item in data.get("runs", [])}


def processes_for_batch(data: dict, batch: str) -> list[dict]:
    return [item for item in data.get("processes", []) if item.get("batch") == batch]


def valid_commit(value: object) -> bool:
    return isinstance(value, str) and bool(re.fullmatch(r"[0-9a-f]{40}", value))


def valid_sha256(value: object) -> bool:
    return isinstance(value, str) and bool(re.fullmatch(r"[0-9a-f]{64}", value))


def hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact_from_path(raw_path: str) -> dict:
    path = Path(raw_path).expanduser().resolve()
    if not path.is_file():
        raise StatusError(f"evidence file does not exist: {path}")
    return {"path": str(path), "size": path.stat().st_size, "sha256": hash_file(path)}


def output_from_path(raw_path: str) -> dict:
    path = Path(raw_path).expanduser().resolve()
    if not path.is_file():
        raise StatusError(f"output file does not exist: {path}")
    return {"path": str(path), "size": path.stat().st_size}


def validate_artifact(value: object, label: str, errors: list[str]) -> None:
    if not isinstance(value, dict):
        errors.append(f"{label} must be an artifact object")
        return
    if not isinstance(value.get("path"), str) or not value["path"].strip():
        errors.append(f"{label}.path must be non-empty")
    if not isinstance(value.get("size"), int) or value["size"] < 0:
        errors.append(f"{label}.size must be a non-negative integer")
    if not valid_sha256(value.get("sha256")):
        errors.append(f"{label}.sha256 must be a lowercase SHA-256")


def validate_bfb_run(run: dict, batch_entries: set[str], errors: list[str]) -> None:
    prefix = f"run {run.get('run_id', '<unknown>')}"
    if run.get("steps") != 50:
        errors.append(f"{prefix} must use exactly 50 steps")
    if run.get("numeric_equal") is not True:
        errors.append(f"{prefix} requires numeric_equal=true")
    if run.get("char_equal") is not True:
        errors.append(f"{prefix} requires char_equal=true")
    if not valid_commit(run.get("source_commit")):
        errors.append(f"{prefix}.source_commit must be a full lowercase commit")
    artifacts = run.get("artifacts")
    if not isinstance(artifacts, dict):
        errors.append(f"{prefix}.artifacts must be an object")
    else:
        for name in REQUIRED_BFB_ARTIFACTS:
            validate_artifact(artifacts.get(name), f"{prefix}.artifacts.{name}", errors)
        outputs = artifacts.get("output_files")
        if not isinstance(outputs, list) or not outputs:
            errors.append(f"{prefix} requires at least one output file")
        else:
            for number, output in enumerate(outputs):
                if not isinstance(output, dict):
                    errors.append(f"{prefix}.output_files[{number}] must be an object")
                    continue
                if not isinstance(output.get("path"), str) or not output["path"].strip():
                    errors.append(f"{prefix}.output_files[{number}].path must be non-empty")
                if not isinstance(output.get("size"), int) or output["size"] < 0:
                    errors.append(f"{prefix}.output_files[{number}].size is invalid")
    proof = run.get("execution_proof")
    if not isinstance(proof, dict) or not proof:
        errors.append(f"{prefix} requires execution proof")
        return
    extra = set(proof).difference(batch_entries)
    if extra:
        errors.append(f"{prefix} proves entries outside its batch: {', '.join(sorted(extra))}")
    for entry, evidence in proof.items():
        if not isinstance(evidence, dict):
            errors.append(f"{prefix}.proof[{entry}] must be an object")
            continue
        validate_artifact(evidence.get("source"), f"{prefix}.proof[{entry}].source", errors)
        if not isinstance(evidence.get("line"), str) or not evidence["line"].strip():
            errors.append(f"{prefix}.proof[{entry}].line must be non-empty")


def collect_status_errors(data: dict, repo_root: Path = REPO_ROOT) -> list[str]:
    errors: list[str] = []
    if data.get("schema_version") != 1:
        errors.append("schema_version must be 1")
    if data.get("allowed_statuses") != list(ALLOWED_STATUSES):
        errors.append("allowed_statuses does not match the fixed phase-2 vocabulary")
    validation = data.get("validation", {})
    if validation.get("steps") != 50:
        errors.append("validation.steps must be exactly 50")
    if tuple(validation.get("completed_statuses", ())) != COMPLETED_STATUSES:
        errors.append("validation.completed_statuses is inconsistent with the checker")

    batches = data.get("batches")
    if not isinstance(batches, list):
        errors.append("batches must be a list")
        batches = []
    batch_map = {item.get("id"): item for item in batches if isinstance(item, dict)}
    if set(batch_map) != set(EXPECTED_BATCHES) or len(batch_map) != len(batches):
        errors.append("batches must be exactly S01 through S04 with unique ids")
    for batch, dependencies in EXPECTED_DEPENDENCIES.items():
        if tuple(batch_map.get(batch, {}).get("dependencies", ())) != dependencies:
            errors.append(f"batch {batch} dependencies must be {list(dependencies)}")

    processes = data.get("processes")
    if not isinstance(processes, list):
        errors.append("processes must be a list")
        processes = []
    entries = [item.get("entry_point") for item in processes if isinstance(item, dict)]
    if len(processes) != 36:
        errors.append(f"standalone inventory must contain 36 entries, found {len(processes)}")
    if len(entries) != len(set(entries)):
        errors.append("entry_point values must be unique")
    by_entry = process_index({"processes": processes})
    expected_all = {entry for values in EXPECTED_BATCHES.values() for entry in values}
    if set(by_entry) != expected_all:
        missing = sorted(expected_all.difference(by_entry))
        extra = sorted(set(by_entry).difference(expected_all))
        if missing:
            errors.append(f"missing standalone entries: {', '.join(missing)}")
        if extra:
            errors.append(f"unexpected standalone entries: {', '.join(extra)}")

    for batch, expected_entries in EXPECTED_BATCHES.items():
        actual = sorted(processes_for_batch({"processes": processes}, batch), key=lambda x: x.get("order", 0))
        if tuple(item.get("entry_point") for item in actual) != expected_entries:
            errors.append(f"batch {batch} entries/order do not match the approved inventory")

    runs = data.get("runs")
    if not isinstance(runs, list):
        errors.append("runs must be a list")
        runs = []
    run_ids = [item.get("run_id") for item in runs if isinstance(item, dict)]
    if len(run_ids) != len(set(run_ids)):
        errors.append("run_id values must be unique")
    runs_by_id = run_index({"runs": runs})
    for run in runs:
        if not isinstance(run, dict):
            errors.append("every run must be an object")
            continue
        batch = run.get("batch")
        members = set(EXPECTED_BATCHES.get(batch, ()))
        if not members:
            errors.append(f"run {run.get('run_id')} references an unknown batch")
        if run.get("result") == "bfb":
            validate_bfb_run(run, members, errors)
        elif run.get("result") == "failed":
            if not isinstance(run.get("note"), str) or not run["note"].strip():
                errors.append(f"failed run {run.get('run_id')} requires a note")
        else:
            errors.append(f"run {run.get('run_id')} has invalid result")

    module_re = re.compile(r"(?im)^\s*module\s+(?!procedure\b|subroutine\b|function\b)([a-z_]\w*)")
    for process in processes:
        if not isinstance(process, dict):
            errors.append("every process must be an object")
            continue
        entry = process.get("entry_point", "<unknown>")
        if process.get("status") not in ALLOWED_STATUSES:
            errors.append(f"{entry} has invalid status")
        commit = process.get("commit")
        if commit is not None and not valid_commit(commit):
            errors.append(f"{entry}.commit must be null or a full lowercase commit")
        source_raw = process.get("source")
        if not isinstance(source_raw, str) or not source_raw.startswith("src/atmos_phys/"):
            errors.append(f"{entry}.source must be under src/atmos_phys")
            continue
        source = (repo_root / source_raw).resolve()
        atmos_root = (repo_root / "src/atmos_phys").resolve()
        try:
            source.relative_to(atmos_root)
        except ValueError:
            errors.append(f"{entry}.source escapes src/atmos_phys")
            continue
        if not source.is_file():
            errors.append(f"{entry}.source does not exist: {source_raw}")
            continue
        text = source.read_text(encoding="utf-8", errors="replace")
        declared_modules = {match.lower() for match in module_re.findall(text)}
        if str(process.get("module", "")).lower() not in declared_modules:
            errors.append(f"{entry}.module is not declared by {source_raw}")
        entry_re = re.compile(rf"(?im)^\s*(?:\w+\s+)*(?:subroutine|function)\s+{re.escape(entry)}\b")
        if not entry_re.search(text):
            errors.append(f"{entry} is not declared by {source_raw}")
        if process.get("status") == "bfb":
            run = runs_by_id.get(process.get("bfb_run"))
            if run is None or run.get("result") != "bfb":
                errors.append(f"{entry} is bfb without a valid bfb_run")
            elif entry not in run.get("execution_proof", {}):
                errors.append(f"{entry} has no execution proof in its bfb_run")
    return errors


def assert_valid(data: dict, repo_root: Path = REPO_ROOT) -> None:
    errors = collect_status_errors(data, repo_root)
    if errors:
        raise StatusError("status validation failed:\n- " + "\n- ".join(errors))


def escape(value: object) -> str:
    return str(value if value not in (None, "") else "-").replace("|", "\\|").replace("\n", " ")


def render_markdown(data: dict) -> str:
    processes = data["processes"]
    counts = {status: sum(item["status"] == status for item in processes) for status in ALLOWED_STATUSES}
    lines = [
        "# PI-atm standalone physics phase-2 progress",
        "",
        "> Generated from `standalone_status.json`; do not edit this file manually.",
        "",
        f"- Baseline: `{data['baseline_commit']}`",
        f"- Updated: `{data['updated_at']}`",
        "- Gate: direct and transitive dependency closure, standalone build, then fresh 50-step BFB.",
        "",
        "## Summary",
        "",
        "| Total | BFB | 50step | Build | Dependency | In progress | Failed | Planned |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
        (
            f"| {len(processes)} | {counts['bfb']} | {counts['50step_running']} | "
            f"{counts['build_pass']} | {counts['dependency_pass']} | {counts['in_progress']} | "
            f"{counts['failed']} | {counts['planned']} |"
        ),
        "",
        "## Batches",
        "",
        "| Batch | Title | Dependencies | Complete | Status |",
        "|---|---|---|---:|---|",
    ]
    batch_map = {item["id"]: item for item in data["batches"]}
    for batch in EXPECTED_BATCHES:
        members = processes_for_batch(data, batch)
        bfb = sum(item["status"] == "bfb" for item in members)
        active = sorted({item["status"] for item in members})
        lines.append(
            f"| {batch} | {escape(batch_map[batch]['title'])} | "
            f"{escape(', '.join(batch_map[batch]['dependencies']))} | {bfb}/9 | {escape(', '.join(active))} |"
        )
    lines.extend(
        [
            "",
            "## Processes",
            "",
            "| Batch | # | Family | Entry | Module | Source | Status | Commit | BFB run | Note |",
            "|---|---:|---|---|---|---|---|---|---|---|",
        ]
    )
    for process in sorted(processes, key=lambda item: (item["batch"], item["order"])):
        lines.append(
            f"| {process['batch']} | {process['order']} | {escape(process['family'])} | "
            f"`{process['entry_point']}` | `{process['module']}` | `{process['source']}` | "
            f"`{process['status']}` | {escape(process.get('commit'))} | "
            f"{escape(process.get('bfb_run'))} | {escape(process.get('note'))} |"
        )
    lines.extend(
        [
            "",
            "## 50-step run records",
            "",
            "| Run | Batch | Result | Job | Source commit | Proofs | Numeric | Character | Recorded |",
            "|---|---|---|---|---|---:|---|---|---|",
        ]
    )
    for run in data.get("runs", []):
        lines.append(
            f"| {escape(run.get('run_id'))} | {escape(run.get('batch'))} | {escape(run.get('result'))} | "
            f"{escape(run.get('job_id'))} | {escape(run.get('source_commit'))} | "
            f"{len(run.get('execution_proof', {}))} | {escape(run.get('numeric_equal'))} | "
            f"{escape(run.get('char_equal'))} | {escape(run.get('recorded_at'))} |"
        )
    lines.extend(
        [
            "",
            "## Standard commands",
            "",
            "```bash",
            "python3 tools/tphys_standalone_progress.py check",
            "python3 tools/tphys_standalone_progress.py render --check",
            "python3 tools/check_atmos_phys_standalone.py --scope completed-only --mode both",
            "python3 tools/check_atmos_phys_standalone.py --scope all --mode both",
            "cmake -S src/atmos_phys/standalone -B <fresh-build-dir>",
            "cmake --build <fresh-build-dir> --target check-standalone",
            "```",
            "",
        ]
    )
    return "\n".join(lines)


def update_processes(
    data: dict,
    entries: list[str],
    status: str,
    commit: str | None,
    note: str | None,
) -> dict:
    if status == "bfb":
        raise StatusError("use record-run to mark BFB")
    if status not in ALLOWED_STATUSES:
        raise StatusError(f"invalid status: {status}")
    if commit is not None and not valid_commit(commit):
        raise StatusError("commit must be a full lowercase 40-character SHA")
    result = copy.deepcopy(data)
    by_entry = process_index(result)
    unknown = sorted(set(entries).difference(by_entry))
    if unknown:
        raise StatusError(f"unknown entries: {', '.join(unknown)}")
    for entry in entries:
        process = by_entry[entry]
        process["status"] = status
        process["bfb_run"] = None
        if commit is not None:
            process["commit"] = commit
        if note is not None:
            process["note"] = note
    result["updated_at"] = utc_now()
    assert_valid(result)
    return result


def parse_proof(value: str) -> tuple[str, str, str]:
    if "=" not in value or "::" not in value:
        raise StatusError("proof must use ENTRY=SOURCE_FILE::MATCHED_LINE")
    entry, payload = value.split("=", 1)
    source, line = payload.split("::", 1)
    if not entry.strip() or not source.strip() or not line.strip():
        raise StatusError("proof entry, source and matched line must be non-empty")
    return entry.strip(), source.strip(), line.strip()


def build_run(args: argparse.Namespace) -> dict:
    proof: dict[str, dict] = {}
    for raw in args.proof:
        entry, source, line = parse_proof(raw)
        if entry in proof:
            raise StatusError(f"duplicate proof for {entry}")
        proof[entry] = {"source": artifact_from_path(source), "line": line}
    artifacts: dict[str, object] = {
        "output_files": [output_from_path(path) for path in args.output_file]
    }
    for name in REQUIRED_BFB_ARTIFACTS:
        raw_path = getattr(args, name)
        if raw_path:
            artifacts[name] = artifact_from_path(raw_path)
    return {
        "run_id": args.run_id,
        "batch": args.batch,
        "result": args.result,
        "steps": 50,
        "job_id": args.job_id,
        "source_commit": args.source_commit,
        "case_root": args.case_root,
        "run_dir": args.run_dir,
        "baseline_run_dir": args.baseline_run_dir,
        "numeric_equal": args.numeric_equal,
        "char_equal": args.char_equal,
        "artifacts": artifacts,
        "execution_proof": proof,
        "recorded_at": args.recorded_at or utc_now(),
        "note": args.note or "",
    }


def record_run(data: dict, run: dict) -> dict:
    result = copy.deepcopy(data)
    if run.get("run_id") in run_index(result):
        raise StatusError(f"duplicate run id: {run.get('run_id')}")
    members = processes_for_batch(result, run.get("batch"))
    if not members:
        raise StatusError(f"unknown batch: {run.get('batch')}")
    batch_entries = {item["entry_point"] for item in members}
    if run.get("result") == "bfb":
        errors: list[str] = []
        validate_bfb_run(run, batch_entries, errors)
        if errors:
            raise StatusError("BFB evidence rejected:\n- " + "\n- ".join(errors))
    elif run.get("result") == "failed":
        if not run.get("note"):
            raise StatusError("a failed run requires --note")
    else:
        raise StatusError("run result must be bfb or failed")
    result.setdefault("runs", []).append(run)
    if run["result"] == "bfb":
        proven = set(run["execution_proof"])
        for process in members:
            if process["entry_point"] in proven:
                process["status"] = "bfb"
                process["commit"] = run["source_commit"]
                process["last_run"] = run["run_id"]
                process["bfb_run"] = run["run_id"]
    else:
        for process in members:
            if process["status"] != "bfb":
                process["status"] = "failed"
                process["last_run"] = run["run_id"]
                process["bfb_run"] = None
    result["updated_at"] = utc_now()
    assert_valid(result)
    return result


def paths_from_args(args: argparse.Namespace) -> tuple[Path, Path]:
    return Path(args.status_file).resolve(), Path(args.progress_file).resolve()


def save_status_and_progress(status_path: Path, progress_path: Path, data: dict) -> None:
    write_json_atomic(status_path, data)
    write_text_atomic(progress_path, render_markdown(data))


def command_check(args: argparse.Namespace) -> int:
    status_path, progress_path = paths_from_args(args)
    data = load_status(status_path)
    repo_root = status_path.resolve().parents[3]
    errors = collect_status_errors(data, repo_root)
    if errors:
        raise StatusError("check failed:\n- " + "\n- ".join(errors))
    try:
        current = progress_path.read_text(encoding="utf-8")
    except OSError as exc:
        raise StatusError(f"cannot read generated progress file: {exc}") from exc
    if current != render_markdown(data):
        raise StatusError(f"{progress_path} is stale; run render")
    print("check passed: 36 standalone entries, four batches, generated Markdown current")
    return 0


def command_render(args: argparse.Namespace) -> int:
    status_path, progress_path = paths_from_args(args)
    data = load_status(status_path)
    repo_root = status_path.resolve().parents[3]
    assert_valid(data, repo_root)
    rendered = render_markdown(data)
    if args.check:
        if not progress_path.is_file() or progress_path.read_text(encoding="utf-8") != rendered:
            raise StatusError(f"{progress_path} is stale; run render without --check")
        print(f"render check passed: {progress_path}")
    else:
        write_text_atomic(progress_path, rendered)
        print(f"rendered {progress_path}")
    return 0


def selected_entries(args: argparse.Namespace, data: dict) -> list[str]:
    if args.batch:
        return [item["entry_point"] for item in processes_for_batch(data, args.batch)]
    return list(args.process)


def command_update(args: argparse.Namespace) -> int:
    status_path, progress_path = paths_from_args(args)
    data = load_status(status_path)
    result = update_processes(
        data,
        selected_entries(args, data),
        args.status,
        args.commit,
        args.note,
    )
    save_status_and_progress(status_path, progress_path, result)
    print(f"updated {len(selected_entries(args, data))} entries to {args.status}")
    return 0


def command_record_run(args: argparse.Namespace) -> int:
    status_path, progress_path = paths_from_args(args)
    data = load_status(status_path)
    result = record_run(data, build_run(args))
    save_status_and_progress(status_path, progress_path, result)
    print(f"recorded {args.run_id}: {args.result}")
    return 0


def bool_value(value: str) -> bool:
    return value.lower() == "true"


def add_common_paths(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--status-file", default=str(DEFAULT_STATUS_FILE))
    parser.add_argument("--progress-file", default=str(DEFAULT_PROGRESS_FILE))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    check = subparsers.add_parser("check", help="validate JSON and generated Markdown")
    add_common_paths(check)
    check.set_defaults(handler=command_check)

    render = subparsers.add_parser("render", help="render generated Markdown")
    add_common_paths(render)
    render.add_argument("--check", action="store_true")
    render.set_defaults(handler=command_render)

    update = subparsers.add_parser("update", help="update workflow status")
    add_common_paths(update)
    target = update.add_mutually_exclusive_group(required=True)
    target.add_argument("--process", action="append", default=[])
    target.add_argument("--batch", choices=tuple(EXPECTED_BATCHES))
    update.add_argument("--status", choices=ALLOWED_STATUSES, required=True)
    update.add_argument("--commit")
    update.add_argument("--note")
    update.set_defaults(handler=command_update)

    run = subparsers.add_parser("record-run", help="record a process-scoped 50-step result")
    add_common_paths(run)
    run.add_argument("--run-id", required=True)
    run.add_argument("--batch", choices=tuple(EXPECTED_BATCHES), required=True)
    run.add_argument("--result", choices=("bfb", "failed"), required=True)
    run.add_argument("--job-id", required=True)
    run.add_argument("--source-commit", required=True)
    run.add_argument("--case-root", required=True)
    run.add_argument("--run-dir", required=True)
    run.add_argument("--baseline-run-dir", required=True)
    run.add_argument("--numeric-equal", type=bool_value, choices=(True, False), required=True)
    run.add_argument("--char-equal", type=bool_value, choices=(True, False), required=True)
    run.add_argument("--output-file", action="append", default=[])
    run.add_argument("--proof", action="append", default=[])
    for name in REQUIRED_BFB_ARTIFACTS:
        run.add_argument("--" + name.replace("_", "-"), dest=name)
    run.add_argument("--recorded-at")
    run.add_argument("--note")
    run.set_defaults(handler=command_record_run)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.handler(args)
    except StatusError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
