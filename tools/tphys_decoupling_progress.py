#!/usr/bin/env python3
"""Track and render the tphysac/tphysbc decoupling campaign.

The JSON status file is the only editable source of truth.  The Markdown
dashboard is generated from it.  A process may be marked ``bfb`` only through
an evidence-complete, 50-step batch validation record.
"""

from __future__ import annotations

import argparse
import copy
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import xml.etree.ElementTree as ET


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_STATUS_FILE = REPO_ROOT / "src/atmos_phys/decoupling_status.json"
DEFAULT_PROGRESS_FILE = REPO_ROOT / "src/atmos_phys/DECOUPLING_PROGRESS.md"
DEFAULT_SUITE_DIR = REPO_ROOT / "src/atmos_phys/suites"

ALLOWED_STATUSES = (
    "planned",
    "in_progress",
    "build_pass",
    "50step_running",
    "bfb",
    "failed",
    "deferred_inactive",
)
ACTIVE_STATUSES = ALLOWED_STATUSES[:-1]
EXPECTED_BATCH_SIZES = {
    "B01": 9,
    "B02": 10,
    "B03": 10,
    "B04": 10,
    "B05": 10,
    "B06": 10,
    "B07": 10,
}
EXPECTED_HOSTS = {"tphysac", "tphysbc"}


class StatusError(ValueError):
    """Raised when status data or validation evidence is invalid."""


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
        raise StatusError("status file root must be a JSON object")
    return data


def write_json_atomic(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(data, indent=2, ensure_ascii=False) + "\n"
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=path.parent, prefix=f".{path.name}.", delete=False
    ) as handle:
        handle.write(payload)
        temp_name = handle.name
    os.replace(temp_name, path)


def write_text_atomic(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=path.parent, prefix=f".{path.name}.", delete=False
    ) as handle:
        handle.write(text)
        temp_name = handle.name
    os.replace(temp_name, path)


def processes_for_batch(data: dict, batch_id: str) -> list[dict]:
    return [item for item in data.get("processes", []) if item.get("batch") == batch_id]


def process_index(data: dict) -> dict[str, dict]:
    return {item["entry_point"]: item for item in data.get("processes", [])}


def run_index(data: dict) -> dict[str, dict]:
    return {item["run_id"]: item for item in data.get("runs", [])}


def _valid_sha256(value: object) -> bool:
    if not isinstance(value, str) or len(value) != 64:
        return False
    return all(char in "0123456789abcdef" for char in value)


def _validate_artifact(artifact: object, label: str, errors: list[str]) -> None:
    if not isinstance(artifact, dict):
        errors.append(f"{label} must be an artifact object")
        return
    if not isinstance(artifact.get("path"), str) or not artifact["path"].strip():
        errors.append(f"{label}.path must be non-empty")
    if not _valid_sha256(artifact.get("sha256")):
        errors.append(f"{label}.sha256 must be a lowercase SHA-256 digest")
    if not isinstance(artifact.get("size"), int) or artifact["size"] < 0:
        errors.append(f"{label}.size must be a non-negative integer")


def _validate_bfb_run(run: dict, batch_entries: set[str], errors: list[str]) -> None:
    prefix = f"run {run.get('run_id', '<unknown>')}"
    if run.get("steps") != 50:
        errors.append(f"{prefix} must use exactly 50 steps")
    if run.get("numeric_equal") is not True:
        errors.append(f"{prefix} cannot be BFB without numeric_equal=true")
    if run.get("char_equal") is not True:
        errors.append(f"{prefix} cannot be BFB without char_equal=true")
    if not isinstance(run.get("source_commit"), str) or len(run["source_commit"]) != 40:
        errors.append(f"{prefix}.source_commit must be a full 40-character commit")

    artifacts = run.get("artifacts")
    if not isinstance(artifacts, dict):
        errors.append(f"{prefix}.artifacts must be an object")
    else:
        for name in (
            "executable",
            "run_environment",
            "filepath",
            "namelist",
            "compare_report",
        ):
            _validate_artifact(artifacts.get(name), f"{prefix}.artifacts.{name}", errors)
        outputs = artifacts.get("output_files")
        if not isinstance(outputs, list) or not outputs:
            errors.append(f"{prefix} cannot be BFB with an empty output file set")
        else:
            for number, output in enumerate(outputs):
                if not isinstance(output, dict):
                    errors.append(f"{prefix}.output_files[{number}] must be an object")
                    continue
                if not isinstance(output.get("path"), str) or not output["path"].strip():
                    errors.append(f"{prefix}.output_files[{number}].path must be non-empty")
                if not isinstance(output.get("size"), int) or output["size"] < 0:
                    errors.append(
                        f"{prefix}.output_files[{number}].size must be a non-negative integer"
                    )

    proof = run.get("execution_proof")
    if not isinstance(proof, dict):
        errors.append(f"{prefix}.execution_proof must be an object")
        return
    missing = sorted(batch_entries.difference(proof))
    extra = sorted(set(proof).difference(batch_entries))
    if missing:
        errors.append(f"{prefix} lacks execution proof for: {', '.join(missing)}")
    if extra:
        errors.append(f"{prefix} has proof for entries outside its batch: {', '.join(extra)}")
    for entry_point, evidence in proof.items():
        if not isinstance(evidence, dict):
            errors.append(f"{prefix}.execution_proof[{entry_point}] must be an object")
            continue
        _validate_artifact(evidence.get("source"), f"{prefix}.proof[{entry_point}].source", errors)
        if not isinstance(evidence.get("line"), str) or not evidence["line"].strip():
            errors.append(f"{prefix}.proof[{entry_point}].line must be non-empty")


def collect_status_errors(data: dict) -> list[str]:
    errors: list[str] = []
    if data.get("schema_version") != 1:
        errors.append("schema_version must be 1")
    if data.get("validation", {}).get("steps") != 50:
        errors.append("validation.steps must be exactly 50")
    if data.get("allowed_statuses") != list(ALLOWED_STATUSES):
        errors.append("allowed_statuses does not match the tool's fixed status vocabulary")
    if not isinstance(data.get("updated_at"), str) or not data["updated_at"]:
        errors.append("updated_at must be non-empty")

    batches = data.get("batches")
    if not isinstance(batches, list):
        errors.append("batches must be a list")
        batches = []
    batch_ids = [item.get("id") for item in batches if isinstance(item, dict)]
    if len(batch_ids) != len(set(batch_ids)):
        errors.append("batch ids must be unique")
    if set(batch_ids) != set(EXPECTED_BATCH_SIZES):
        errors.append("batches must be exactly B01 through B07")
    batch_map = {item.get("id"): item for item in batches if isinstance(item, dict)}
    for batch_id, item in batch_map.items():
        dependencies = item.get("dependencies")
        if not isinstance(dependencies, list):
            errors.append(f"batch {batch_id}.dependencies must be a list")
            continue
        unknown = sorted(set(dependencies).difference(batch_map))
        if unknown:
            errors.append(f"batch {batch_id} has unknown dependencies: {', '.join(unknown)}")
        if batch_id in dependencies:
            errors.append(f"batch {batch_id} cannot depend on itself")

    processes = data.get("processes")
    if not isinstance(processes, list):
        errors.append("processes must be a list")
        processes = []
    entries = [item.get("entry_point") for item in processes if isinstance(item, dict)]
    if len(processes) != 69:
        errors.append(f"active process inventory must contain 69 entries, found {len(processes)}")
    if len(entries) != len(set(entries)):
        errors.append("active entry_point values must be unique")

    for batch_id, expected_count in EXPECTED_BATCH_SIZES.items():
        actual_count = sum(1 for item in processes if item.get("batch") == batch_id)
        if actual_count != expected_count:
            errors.append(
                f"batch {batch_id} must contain {expected_count} active entries, found {actual_count}"
            )

    runs = data.get("runs")
    if not isinstance(runs, list):
        errors.append("runs must be a list")
        runs = []
    run_ids = [item.get("run_id") for item in runs if isinstance(item, dict)]
    if len(run_ids) != len(set(run_ids)):
        errors.append("run_id values must be unique")
    runs_by_id = {item.get("run_id"): item for item in runs if isinstance(item, dict)}

    for number, process in enumerate(processes):
        prefix = f"process[{number}]"
        if not isinstance(process, dict):
            errors.append(f"{prefix} must be an object")
            continue
        entry = process.get("entry_point")
        prefix = f"process {entry or number}"
        for field in ("entry_point", "scheme", "lifecycle", "batch", "family"):
            if not isinstance(process.get(field), str) or not process[field].strip():
                errors.append(f"{prefix}.{field} must be non-empty")
        if process.get("lifecycle") not in {"init", "timestep_init", "run", "timestep_final", "final"}:
            errors.append(f"{prefix}.lifecycle is invalid")
        if process.get("batch") not in batch_map:
            errors.append(f"{prefix} references an unknown batch")
        hosts = process.get("hosts")
        if not isinstance(hosts, list) or not hosts:
            errors.append(f"{prefix}.hosts must be a non-empty list")
        elif not set(hosts).issubset(EXPECTED_HOSTS):
            errors.append(f"{prefix}.hosts contains an unknown host")
        status = process.get("status")
        if status not in ACTIVE_STATUSES:
            errors.append(f"{prefix}.status must be an active status")
        for field in ("commit", "last_run", "bfb_run"):
            if process.get(field) is not None and not isinstance(process[field], str):
                errors.append(f"{prefix}.{field} must be a string or null")
        if process.get("commit") is not None and len(process["commit"]) != 40:
            errors.append(f"{prefix}.commit must be a full 40-character commit")
        last_run = process.get("last_run")
        if last_run is not None:
            if last_run not in runs_by_id:
                errors.append(f"{prefix}.last_run references an unknown run")
            elif runs_by_id[last_run].get("batch") != process.get("batch"):
                errors.append(f"{prefix}.last_run belongs to a different batch")
        bfb_run = process.get("bfb_run")
        if status == "bfb":
            if bfb_run not in runs_by_id:
                errors.append(f"{prefix} is bfb without a valid bfb_run")
            else:
                run = runs_by_id[bfb_run]
                if run.get("result") != "bfb" or run.get("batch") != process.get("batch"):
                    errors.append(f"{prefix}.bfb_run is not a BFB run for its batch")
                elif entry not in run.get("execution_proof", {}):
                    errors.append(f"{prefix}.bfb_run lacks execution proof for this entry")
        elif bfb_run is not None:
            errors.append(f"{prefix} is not bfb but still has bfb_run")

    all_entries = set(entries)
    for number, run in enumerate(runs):
        prefix = f"run[{number}]"
        if not isinstance(run, dict):
            errors.append(f"{prefix} must be an object")
            continue
        run_id = run.get("run_id")
        prefix = f"run {run_id or number}"
        batch_id = run.get("batch")
        if batch_id not in batch_map:
            errors.append(f"{prefix} references an unknown batch")
            continue
        if run.get("result") not in {"bfb", "failed"}:
            errors.append(f"{prefix}.result must be bfb or failed")
        proof = run.get("execution_proof", {})
        if isinstance(proof, dict):
            unknown_proof = sorted(set(proof).difference(all_entries))
            if unknown_proof:
                errors.append(f"{prefix} references unknown proof entries: {', '.join(unknown_proof)}")
        if run.get("result") == "bfb":
            batch_entries = {
                item["entry_point"] for item in processes if item.get("batch") == batch_id
            }
            _validate_bfb_run(run, batch_entries, errors)

    deferred = data.get("deferred")
    if not isinstance(deferred, list):
        errors.append("deferred must be a list")
    else:
        deferred_names = [item.get("name") for item in deferred if isinstance(item, dict)]
        if len(deferred_names) != len(set(deferred_names)):
            errors.append("deferred names must be unique")
        for number, item in enumerate(deferred):
            if not isinstance(item, dict):
                errors.append(f"deferred[{number}] must be an object")
                continue
            if item.get("status") != "deferred_inactive":
                errors.append(f"deferred {item.get('name', number)} must use deferred_inactive")
            if not isinstance(item.get("reason"), str) or not item["reason"].strip():
                errors.append(f"deferred {item.get('name', number)} needs a reason")
    return errors


def assert_valid_status(data: dict) -> None:
    errors = collect_status_errors(data)
    if errors:
        raise StatusError("invalid status data:\n- " + "\n- ".join(errors))


def collect_suite_errors(data: dict, suite_dir: Path) -> list[str]:
    errors: list[str] = []
    for host in sorted(EXPECTED_HOSTS):
        path = suite_dir / f"suite_{host}_active.xml"
        try:
            root = ET.parse(path).getroot()
        except (OSError, ET.ParseError) as exc:
            errors.append(f"cannot read suite {path}: {exc}")
            continue
        expected_name = f"{host}_active"
        if root.tag != "suite" or root.get("name") != expected_name:
            errors.append(f"{path.name} must define suite name={expected_name!r}")
        listed = [node.text.strip() for node in root.findall(".//scheme") if node.text and node.text.strip()]
        listed_set = set(listed)
        expected = {
            item["scheme"] for item in data["processes"] if host in item.get("hosts", [])
        }
        missing = sorted(expected.difference(listed_set))
        unknown = sorted(listed_set.difference(expected))
        if missing:
            errors.append(f"{path.name} is missing active schemes: {', '.join(missing)}")
        if unknown:
            errors.append(f"{path.name} contains schemes outside active inventory: {', '.join(unknown)}")
    return errors


def _escape(value: object) -> str:
    if value is None or value == "":
        return "—"
    return str(value).replace("|", "\\|").replace("\n", " ")


def _batch_status(processes: list[dict]) -> str:
    statuses = [item["status"] for item in processes]
    if statuses and all(status == "bfb" for status in statuses):
        return "bfb"
    for status in ("failed", "50step_running", "build_pass", "in_progress"):
        if status in statuses:
            return status
    return "planned"


def render_markdown(data: dict) -> str:
    assert_valid_status(data)
    processes = data["processes"]
    counts = {status: sum(item["status"] == status for item in processes) for status in ACTIVE_STATUSES}
    lines = [
        "<!-- Generated by tools/tphys_decoupling_progress.py; do not edit manually. -->",
        "# `tphysac` / `tphysbc` 解耦进度",
        "",
        f"- JSON 唯一真源：`src/atmos_phys/decoupling_status.json`",
        f"- 基线提交：`{data['baseline_commit']}`",
        f"- 当前 PI 配置：{data['active_case']['description']}",
        f"- 验证门槛：每组 `{data['validation']['steps']}step`，numeric、char 和执行证据必须全部通过",
        f"- 最近更新：`{data['updated_at']}`",
        "",
        "## 总览",
        "",
        "| Active总数 | BFB | 进行中 | Build通过 | 50step运行中 | 失败 | 待开始 |",
        "|---:|---:|---:|---:|---:|---:|---:|",
        (
            f"| {len(processes)} | {counts['bfb']} | {counts['in_progress']} | "
            f"{counts['build_pass']} | {counts['50step_running']} | {counts['failed']} | "
            f"{counts['planned']} |"
        ),
        "",
        "## 批次",
        "",
        "| 批次 | 模块组 | 入口数 | 依赖 | 当前状态 | BFB入口 |",
        "|---|---|---:|---|---|---:|",
    ]
    for batch in data["batches"]:
        members = processes_for_batch(data, batch["id"])
        bfb_count = sum(item["status"] == "bfb" for item in members)
        dependencies = ", ".join(batch["dependencies"]) or "—"
        lines.append(
            f"| {batch['id']} | {_escape(batch['title'])} | {len(members)} | "
            f"{dependencies} | `{_batch_status(members)}` | {bfb_count}/{len(members)} |"
        )

    lines.extend(
        [
            "",
            "## Active入口",
            "",
            "| 批次 | Host | Family | 导出入口 | 状态 | 提交 | 最近运行 | BFB运行 | 备注 |",
            "|---|---|---|---|---|---|---|---|---|",
        ]
    )
    batch_order = {item["id"]: number for number, item in enumerate(data["batches"])}
    for process in sorted(processes, key=lambda item: (batch_order[item["batch"]], item["order"])):
        lines.append(
            "| {batch} | {hosts} | {family} | `{entry}` | `{status}` | {commit} | "
            "{last_run} | {bfb_run} | {note} |".format(
                batch=process["batch"],
                hosts=", ".join(process["hosts"]),
                family=_escape(process["family"]),
                entry=process["entry_point"],
                status=process["status"],
                commit=_escape(process.get("commit")),
                last_run=_escape(process.get("last_run")),
                bfb_run=_escape(process.get("bfb_run")),
                note=_escape(process.get("note")),
            )
        )

    lines.extend(
        [
            "",
            "## 暂缓的非活动分支",
            "",
            "这些分支不计入本阶段69个active入口及BFB完成率；切换到能实际执行该分支的专用case后另组验证。",
            "",
            "| 过程/分支 | Host | 状态 | 原因 |",
            "|---|---|---|---|",
        ]
    )
    for item in data["deferred"]:
        lines.append(
            f"| `{item['name']}` | {_escape(', '.join(item['hosts']))} | "
            f"`{item['status']}` | {_escape(item['reason'])} |"
        )

    lines.extend(
        [
            "",
            "## 50step运行记录",
            "",
            "| Run ID | 批次 | 结果 | Job ID | Commit | Numeric | Char | 执行证据 | 时间 |",
            "|---|---|---|---|---|---|---|---:|---|",
        ]
    )
    if data["runs"]:
        for run in data["runs"]:
            lines.append(
                f"| `{run['run_id']}` | {run['batch']} | `{run['result']}` | "
                f"{_escape(run.get('job_id'))} | {_escape(run.get('source_commit'))} | "
                f"{_escape(run.get('numeric_equal'))} | {_escape(run.get('char_equal'))} | "
                f"{len(run.get('execution_proof', {}))} | {_escape(run.get('recorded_at'))} |"
            )
    else:
        lines.append("| — | — | — | — | — | — | — | 0 | — |")
    lines.append("")
    return "\n".join(lines)


def hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def artifact_from_path(path_text: str) -> dict:
    path = Path(path_text).expanduser().resolve()
    if not path.is_file():
        raise StatusError(f"evidence file does not exist or is not a regular file: {path}")
    return {"path": str(path), "size": path.stat().st_size, "sha256": hash_file(path)}


def output_from_path(path_text: str) -> dict:
    path = Path(path_text).expanduser().resolve()
    if not path.is_file():
        raise StatusError(f"output file does not exist or is not a regular file: {path}")
    return {"path": str(path), "size": path.stat().st_size}


def parse_proof(value: str) -> tuple[str, str, str]:
    if "=" not in value or "::" not in value:
        raise StatusError("proof must use ENTRY_POINT=SOURCE_FILE::MATCHED_LINE")
    entry_point, payload = value.split("=", 1)
    source, line = payload.split("::", 1)
    if not entry_point.strip() or not source.strip() or not line.strip():
        raise StatusError("proof entry point, source file, and matched line must be non-empty")
    return entry_point.strip(), source.strip(), line.strip()


def update_processes(
    data: dict,
    entries: list[str],
    status: str,
    commit: str | None = None,
    note: str | None = None,
) -> dict:
    if status not in ACTIVE_STATUSES:
        raise StatusError(f"invalid active status: {status}")
    result = copy.deepcopy(data)
    by_entry = process_index(result)
    unknown = sorted(set(entries).difference(by_entry))
    if unknown:
        raise StatusError(f"unknown active entries: {', '.join(unknown)}")
    if commit is not None and (len(commit) != 40 or any(c not in "0123456789abcdef" for c in commit)):
        raise StatusError("commit must be a full lowercase 40-character git SHA")
    for entry in entries:
        process = by_entry[entry]
        if status == "bfb":
            run_id = process.get("bfb_run")
            run = run_index(result).get(run_id)
            if run is None or run.get("result") != "bfb" or entry not in run.get("execution_proof", {}):
                raise StatusError(
                    f"cannot mark {entry} bfb without an evidence-complete BFB record; use record-run"
                )
        else:
            process["bfb_run"] = None
        process["status"] = status
        if commit is not None:
            process["commit"] = commit
        if note is not None:
            process["note"] = note
    result["updated_at"] = utc_now()
    assert_valid_status(result)
    return result


def record_run(data: dict, run: dict, verify_paths: bool = True) -> dict:
    result = copy.deepcopy(data)
    if run.get("run_id") in run_index(result):
        raise StatusError(f"run_id already exists: {run.get('run_id')}")
    batch_id = run.get("batch")
    members = processes_for_batch(result, batch_id)
    if not members:
        raise StatusError(f"unknown or empty batch: {batch_id}")
    batch_entries = {item["entry_point"] for item in members}
    if run.get("result") == "bfb":
        evidence_errors: list[str] = []
        _validate_bfb_run(run, batch_entries, evidence_errors)
        if evidence_errors:
            raise StatusError("BFB evidence rejected:\n- " + "\n- ".join(evidence_errors))
        if verify_paths:
            artifacts = run["artifacts"]
            for name in ("executable", "run_environment", "filepath", "namelist", "compare_report"):
                path = Path(artifacts[name]["path"])
                if not path.is_file():
                    raise StatusError(f"recorded {name} evidence no longer exists: {path}")
            for output in artifacts["output_files"]:
                if not Path(output["path"]).is_file():
                    raise StatusError(f"recorded output file no longer exists: {output['path']}")
            for entry, proof in run["execution_proof"].items():
                if not Path(proof["source"]["path"]).is_file():
                    raise StatusError(f"execution proof source for {entry} does not exist")
    elif run.get("result") != "failed":
        raise StatusError("run result must be bfb or failed")
    elif not isinstance(run.get("note"), str) or not run["note"].strip():
        raise StatusError("a failed run requires a non-empty note")

    result.setdefault("runs", []).append(run)
    for process in members:
        process["last_run"] = run["run_id"]
        if run["result"] == "bfb":
            process["status"] = "bfb"
            process["bfb_run"] = run["run_id"]
            process["commit"] = run["source_commit"]
        elif process["status"] != "bfb":
            process["status"] = "failed"
            process["bfb_run"] = None
    result["updated_at"] = utc_now()
    assert_valid_status(result)
    return result


def _bool_value(value: str) -> bool:
    return value.lower() == "true"


def build_run_from_args(args: argparse.Namespace, data: dict) -> dict:
    proof: dict[str, dict] = {}
    for raw in args.proof:
        entry, source, line = parse_proof(raw)
        if entry in proof:
            raise StatusError(f"duplicate execution proof for {entry}")
        proof[entry] = {"source": artifact_from_path(source), "line": line}

    artifacts: dict[str, object] = {"output_files": [output_from_path(path) for path in args.output_file]}
    artifact_args = {
        "executable": args.executable,
        "run_environment": args.run_environment,
        "filepath": args.filepath,
        "namelist": args.namelist,
        "compare_report": args.compare_report,
    }
    for name, path in artifact_args.items():
        if path:
            artifacts[name] = artifact_from_path(path)

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


def _paths_from_args(args: argparse.Namespace) -> tuple[Path, Path]:
    return Path(args.status_file).resolve(), Path(args.progress_file).resolve()


def command_render(args: argparse.Namespace) -> int:
    status_path, progress_path = _paths_from_args(args)
    data = load_status(status_path)
    rendered = render_markdown(data)
    if args.check:
        try:
            current = progress_path.read_text(encoding="utf-8")
        except OSError as exc:
            raise StatusError(f"cannot read generated progress file {progress_path}: {exc}") from exc
        if current != rendered:
            raise StatusError(
                f"{progress_path} is stale; run: {Path(__file__).name} render"
            )
        print(f"render check passed: {progress_path}")
        return 0
    write_text_atomic(progress_path, rendered)
    print(f"rendered {progress_path}")
    return 0


def command_check(args: argparse.Namespace) -> int:
    status_path, progress_path = _paths_from_args(args)
    data = load_status(status_path)
    errors = collect_status_errors(data)
    errors.extend(collect_suite_errors(data, Path(args.suite_dir).resolve()))
    if errors:
        raise StatusError("check failed:\n- " + "\n- ".join(errors))
    expected = render_markdown(data)
    try:
        actual = progress_path.read_text(encoding="utf-8")
    except OSError as exc:
        raise StatusError(f"cannot read generated progress file {progress_path}: {exc}") from exc
    if actual != expected:
        raise StatusError(f"generated progress file is stale: {progress_path}")
    print(
        f"check passed: {len(data['processes'])} active entries, "
        f"{len(data['deferred'])} deferred entries, suites covered"
    )
    return 0


def command_update(args: argparse.Namespace) -> int:
    status_path, progress_path = _paths_from_args(args)
    data = load_status(status_path)
    if args.process:
        entries = args.process
    else:
        entries = [item["entry_point"] for item in processes_for_batch(data, args.batch)]
        if not entries:
            raise StatusError(f"unknown batch: {args.batch}")
    updated = update_processes(data, entries, args.status, args.commit, args.note)
    write_json_atomic(status_path, updated)
    write_text_atomic(progress_path, render_markdown(updated))
    print(f"updated {len(entries)} process(es) to {args.status}")
    return 0


def command_record_run(args: argparse.Namespace) -> int:
    status_path, progress_path = _paths_from_args(args)
    data = load_status(status_path)
    run = build_run_from_args(args, data)
    updated = record_run(data, run)
    write_json_atomic(status_path, updated)
    write_text_atomic(progress_path, render_markdown(updated))
    print(f"recorded {run['run_id']} as {run['result']} for {run['batch']}")
    return 0


def add_common_paths(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--status-file", default=str(DEFAULT_STATUS_FILE))
    parser.add_argument("--progress-file", default=str(DEFAULT_PROGRESS_FILE))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    render_parser = subparsers.add_parser("render", help="render the Markdown dashboard")
    add_common_paths(render_parser)
    render_parser.add_argument("--check", action="store_true", help="fail if Markdown is stale")
    render_parser.set_defaults(func=command_render)

    check_parser = subparsers.add_parser("check", help="validate JSON, suites, and Markdown")
    add_common_paths(check_parser)
    check_parser.add_argument("--suite-dir", default=str(DEFAULT_SUITE_DIR))
    check_parser.set_defaults(func=command_check)

    update_parser = subparsers.add_parser("update", help="update process workflow status")
    add_common_paths(update_parser)
    target = update_parser.add_mutually_exclusive_group(required=True)
    target.add_argument("--process", action="append", help="entry point; may be repeated")
    target.add_argument("--batch", help="update every active entry in a batch")
    update_parser.add_argument("--status", choices=ACTIVE_STATUSES, required=True)
    update_parser.add_argument("--commit")
    update_parser.add_argument("--note")
    update_parser.set_defaults(func=command_update)

    run_parser = subparsers.add_parser("record-run", help="record a batch 50-step result")
    add_common_paths(run_parser)
    run_parser.add_argument("--run-id", required=True)
    run_parser.add_argument("--batch", choices=tuple(EXPECTED_BATCH_SIZES), required=True)
    run_parser.add_argument("--result", choices=("bfb", "failed"), required=True)
    run_parser.add_argument("--job-id", required=True)
    run_parser.add_argument("--source-commit", required=True)
    run_parser.add_argument("--case-root", required=True)
    run_parser.add_argument("--run-dir", required=True)
    run_parser.add_argument("--baseline-run-dir", required=True)
    run_parser.add_argument("--numeric-equal", type=_bool_value, choices=(True, False), required=True)
    run_parser.add_argument("--char-equal", type=_bool_value, choices=(True, False), required=True)
    run_parser.add_argument("--output-file", action="append", default=[])
    run_parser.add_argument("--proof", action="append", default=[], metavar="ENTRY=FILE::LINE")
    run_parser.add_argument("--executable")
    run_parser.add_argument("--run-environment")
    run_parser.add_argument("--filepath")
    run_parser.add_argument("--namelist")
    run_parser.add_argument("--compare-report")
    run_parser.add_argument("--recorded-at")
    run_parser.add_argument("--note")
    run_parser.set_defaults(func=command_record_run)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except StatusError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
