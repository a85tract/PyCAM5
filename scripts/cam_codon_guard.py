#!/usr/bin/env python3
from __future__ import annotations

import argparse
import datetime as _datetime
import fnmatch
import gzip
import hashlib
import json
import os
from pathlib import Path
import py_compile
import re
import subprocess
import sys
import tempfile
from typing import Iterable

try:
    import yaml
except ImportError:  # pragma: no cover - exercised only on systems without PyYAML
    yaml = None


DEFAULT_CONFIG = ".codon_guard.yaml"
DEFAULT_RECEIPT_DIR = ".codon_guard_receipts"
ALLOW_HIGH_RISK_ENV = "CAM_CODON_GUARD_ALLOW_HIGH_RISK"
DEFAULT_NUMERIC_HAZARD_OVERRIDE_ENV = "CAM_CODON_GUARD_ALLOW_NUMERIC_HAZARD"


class GuardError(Exception):
    pass


def run_cmd(
    args: list[str],
    cwd: Path,
    *,
    check: bool = False,
    text: bool = True,
) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(
        args,
        cwd=str(cwd),
        text=text,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if check and result.returncode != 0:
        raise GuardError(
            f"command failed ({result.returncode}): {' '.join(args)}\n"
            f"{result.stdout}{result.stderr}"
        )
    return result


def repo_root() -> Path:
    result = run_cmd(["git", "rev-parse", "--show-toplevel"], Path.cwd())
    if result.returncode != 0:
        raise GuardError("not inside a git repository")
    return Path(result.stdout.strip()).resolve()


def load_config(root: Path, config_arg: str | None) -> dict:
    config_path = Path(config_arg or DEFAULT_CONFIG)
    if not config_path.is_absolute():
        config_path = root / config_path
    if yaml is None:
        raise GuardError("PyYAML is required: module 'yaml' is not available")
    if not config_path.is_file():
        raise GuardError(f"missing guard config: {config_path}")
    with config_path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise GuardError(f"guard config must be a YAML mapping: {config_path}")
    if data.get("version") != 1:
        raise GuardError("guard config version must be 1")
    return data


def config_list(config: dict, key: str) -> list[str]:
    value = config.get(key, [])
    if value is None:
        return []
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise GuardError(f"config key {key!r} must be a list of strings")
    return value


def config_string(config: dict, key: str, default: str) -> str:
    value = config.get(key, default)
    if not isinstance(value, str):
        raise GuardError(f"config key {key!r} must be a string")
    return value


def config_rule_list(config: dict, key: str) -> list[dict]:
    value = config.get(key, [])
    if value is None:
        return []
    if not isinstance(value, list) or not all(isinstance(item, dict) for item in value):
        raise GuardError(f"config key {key!r} must be a list of mappings")
    return value


def config_path(config: dict, key: str, *, required: bool = True) -> Path | None:
    paths = config.get("paths", {})
    if not isinstance(paths, dict):
        raise GuardError("config key 'paths' must be a mapping")
    value = paths.get(key)
    if value is None:
        if required:
            raise GuardError(f"config paths.{key} is required")
        return None
    if not isinstance(value, str):
        raise GuardError(f"config paths.{key} must be a string")
    return Path(value)


def receipt_dir(root: Path, config: dict) -> Path:
    paths = config.get("paths", {})
    value = paths.get("receipt_dir") if isinstance(paths, dict) else None
    path = Path(value) if isinstance(value, str) else Path(DEFAULT_RECEIPT_DIR)
    if not path.is_absolute():
        path = root / path
    return path


def staged_files(root: Path) -> list[str]:
    result = run_cmd(["git", "diff", "--cached", "--name-only", "-z"], root, check=True)
    return [item for item in result.stdout.split("\0") if item]


def current_head(root: Path) -> str:
    result = run_cmd(["git", "rev-parse", "HEAD"], root, check=True)
    return result.stdout.strip()


def staged_diff_text(root: Path) -> str:
    result = run_cmd(
        ["git", "diff", "--cached", "--full-index", "--binary", "--no-ext-diff", "--"],
        root,
        check=True,
    )
    return result.stdout


def staged_diff_hash(root: Path) -> str:
    text = staged_diff_text(root)
    return hashlib.sha256(text.encode("utf-8", errors="surrogateescape")).hexdigest()


def receipt_path(root: Path, config: dict, diff_hash: str) -> Path:
    return receipt_dir(root, config) / f"{diff_hash}.json"


def load_matching_receipt(root: Path, config: dict, diff_hash: str) -> dict | None:
    path = receipt_path(root, config, diff_hash)
    if not path.is_file():
        return None
    try:
        with path.open("r", encoding="utf-8") as handle:
            data = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(data, dict):
        return None
    if data.get("version") != 1:
        return None
    if data.get("staged_diff_hash") != diff_hash:
        return None
    if data.get("git_head") != current_head(root):
        return None
    if data.get("overall_numeric_equal") is not True:
        return None
    return data


def path_matches(path: str, patterns: Iterable[str]) -> bool:
    normalized = path.replace(os.sep, "/")
    for pattern in patterns:
        if fnmatch.fnmatch(normalized, pattern):
            return True
        if Path(normalized).match(pattern):
            return True
    return False


def git_diff_check(root: Path) -> list[str]:
    result = run_cmd(["git", "diff", "--cached", "--check"], root)
    if result.returncode == 0:
        return []
    return [line for line in (result.stdout + result.stderr).splitlines() if line.strip()]


def staged_diff_lines(root: Path) -> list[tuple[str, str, int | None, str]]:
    result = run_cmd(
        ["git", "diff", "--cached", "--unified=0", "--no-ext-diff", "--"],
        root,
        check=True,
    )
    entries: list[tuple[str, str, int | None, str]] = []
    old_path: str | None = None
    new_path: str | None = None
    old_line: int | None = None
    new_line: int | None = None
    hunk_re = re.compile(r"^@@ -(?P<old>\d+)(?:,\d+)? \+(?P<new>\d+)(?:,\d+)? @@")

    for line in result.stdout.splitlines():
        if line.startswith("--- a/"):
            old_path = line[6:]
            continue
        if line.startswith("--- /dev/null"):
            old_path = None
            continue
        if line.startswith("+++ b/"):
            new_path = line[6:]
            continue
        if line.startswith("+++ /dev/null"):
            new_path = None
            continue
        match = hunk_re.match(line)
        if match:
            old_line = int(match.group("old"))
            new_line = int(match.group("new"))
            continue
        if line.startswith("+") and not line.startswith("+++"):
            if new_path:
                entries.append((new_path, "added", new_line, line[1:]))
            if new_line is not None:
                new_line += 1
        elif line.startswith("-") and not line.startswith("---"):
            if old_path:
                entries.append((old_path, "removed", old_line, line[1:]))
            if old_line is not None:
                old_line += 1
        elif not line.startswith("\\"):
            if old_line is not None:
                old_line += 1
            if new_line is not None:
                new_line += 1
    return entries


def numeric_hazard_findings(
    root: Path, config: dict, *, receipt_unlocked: bool = False
) -> tuple[list[str], list[str]]:
    warnings: list[str] = []
    failures: list[str] = []
    override_env = config_string(
        config, "numeric_hazard_override_env", DEFAULT_NUMERIC_HAZARD_OVERRIDE_ENV
    )
    override_enabled = os.environ.get(override_env) == "1"
    allowed_policies = {"warn", "require_override", "hard_fail"}
    allowed_polarities = {"added", "removed"}
    lines = staged_diff_lines(root)

    for rule in config_rule_list(config, "numeric_hazard_rules"):
        name = rule.get("name")
        policy = rule.get("policy")
        polarity = rule.get("polarity")
        file_glob = rule.get("file_glob")
        regex = rule.get("regex")
        if not all(isinstance(item, str) for item in (name, policy, polarity, file_glob, regex)):
            raise GuardError("numeric_hazard_rules entries require string name/policy/polarity/file_glob/regex")
        if policy not in allowed_policies:
            raise GuardError(f"numeric hazard rule {name!r} has invalid policy {policy!r}")
        if polarity not in allowed_polarities:
            raise GuardError(f"numeric hazard rule {name!r} has invalid polarity {polarity!r}")
        try:
            compiled = re.compile(regex, re.IGNORECASE)
        except re.error as exc:
            raise GuardError(f"numeric hazard rule {name!r} has invalid regex: {exc}") from exc

        for path, line_polarity, lineno, text in lines:
            if line_polarity != polarity or not path_matches(path, [file_glob]):
                continue
            if not compiled.search(text):
                continue
            location = f"{path}:{lineno}" if lineno is not None else path
            message = f"{name} ({policy}) at {location}: {text.strip()}"
            if policy == "warn":
                warnings.append(message)
            elif policy == "require_override":
                if receipt_unlocked:
                    warnings.append(message + " [covered by matching BFB receipt]")
                elif override_enabled:
                    warnings.append(message + f" [override {override_env}=1]")
                else:
                    failures.append(message + f" [set {override_env}=1 to allow]")
            else:
                failures.append(message)
    return warnings, failures


def py_compile_staged(root: Path, files: list[str], patterns: list[str]) -> list[str]:
    errors: list[str] = []
    for relpath in files:
        if not path_matches(relpath, patterns):
            continue
        path = root / relpath
        if not path.is_file():
            continue
        try:
            with tempfile.TemporaryDirectory(prefix="cam_codon_guard_pyc_") as tmpdir:
                cfile = Path(tmpdir) / (Path(relpath).name + "c")
                py_compile.compile(str(path), cfile=str(cfile), doraise=True)
        except py_compile.PyCompileError as exc:
            errors.append(f"{relpath}: {exc.msg}")
    return errors


def cache_warnings(root: Path) -> list[str]:
    result = run_cmd(["git", "status", "--porcelain=v1", "-z", "--untracked-files=all"], root)
    warnings: list[str] = []
    for entry in result.stdout.split("\0"):
        if not entry:
            continue
        path = entry[3:] if len(entry) > 3 else entry
        if "__pycache__/" in path or path.endswith(".pyc"):
            warnings.append(path)
    return sorted(set(warnings))


def pre_commit(args: argparse.Namespace) -> int:
    root = repo_root()
    config = load_config(root, args.config)
    files = staged_files(root)
    diff_hash = staged_diff_hash(root)
    receipt = load_matching_receipt(root, config, diff_hash) if files else None
    receipt_unlocked = receipt is not None
    failures: list[str] = []

    forbidden = config_list(config, "forbidden_staged_globs")
    high_risk = config_list(config, "high_risk_staged_globs")
    compile_globs = config_list(config, "python_compile_staged_globs")

    forbidden_hits = [path for path in files if path_matches(path, forbidden)]
    if forbidden_hits:
        failures.append("forbidden staged paths:\n  " + "\n  ".join(forbidden_hits))

    high_risk_hits = [path for path in files if path_matches(path, high_risk)]
    if (
        high_risk_hits
        and not receipt_unlocked
        and os.environ.get(ALLOW_HIGH_RISK_ENV) != "1"
    ):
        failures.append(
            "high-risk staged paths require a matching BFB receipt from "
            "`validate-run`, or explicit override "
            f"{ALLOW_HIGH_RISK_ENV}=1:\n  " + "\n  ".join(high_risk_hits)
        )

    diff_errors = git_diff_check(root)
    if diff_errors:
        failures.append("git diff --cached --check failed:\n  " + "\n  ".join(diff_errors))

    compile_errors = py_compile_staged(root, files, compile_globs)
    if compile_errors:
        failures.append("python syntax checks failed:\n  " + "\n  ".join(compile_errors))

    warnings = cache_warnings(root)
    numeric_warnings, numeric_failures = numeric_hazard_findings(
        root, config, receipt_unlocked=receipt_unlocked
    )
    if numeric_failures:
        failures.append("numeric hazard checks failed:\n  " + "\n  ".join(numeric_failures))
    if warnings:
        print("CAM/Codon guard warning: cache files are present but not staged:", file=sys.stderr)
        for path in warnings[:20]:
            print(f"  {path}", file=sys.stderr)
        if len(warnings) > 20:
            print(f"  ... {len(warnings) - 20} more", file=sys.stderr)
    if numeric_warnings:
        print("CAM/Codon guard warning: numeric hazards found in staged diff:", file=sys.stderr)
        for warning in numeric_warnings[:30]:
            print(f"  {warning}", file=sys.stderr)
        if len(numeric_warnings) > 30:
            print(f"  ... {len(numeric_warnings) - 30} more", file=sys.stderr)
    if receipt_unlocked:
        print(
            "CAM/Codon guard: staged diff has matching BFB receipt "
            f"for job {receipt.get('job', '<unknown>')}.",
            file=sys.stderr,
        )

    if failures:
        print("CAM/Codon guard failed pre-commit checks:", file=sys.stderr)
        for failure in failures:
            print(f"\n{failure}", file=sys.stderr)
        return 1

    print(f"CAM/Codon guard pre-commit passed ({len(files)} staged file(s)).")
    return 0


def find_run_environment(case_root: Path, job: str) -> Path:
    logs_dir = case_root / "logs"
    candidates = sorted(logs_dir.glob(f"run_environment.txt.{job}*"))
    if not candidates and "." in job:
        candidates = sorted(logs_dir.glob(f"run_environment.txt.{job.split('.')[0]}*"))
    if not candidates:
        raise GuardError(f"no run_environment.txt found for job {job} in {logs_dir}")
    return candidates[-1]


def read_text_file(path: Path) -> str:
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        return handle.read()


def read_gzip_file(path: Path) -> str:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        return handle.read()


def proof_candidates(run_dir: Path, proof_files: list[str]) -> list[Path]:
    candidates = sorted(run_dir.glob("atm.log.*.gz"))
    candidates.extend(sorted(run_dir.rglob("*.proof")))
    candidates.extend(sorted(run_dir.rglob("*.proof.*")))
    for item in proof_files:
        candidates.append(Path(item))
    seen: set[Path] = set()
    unique: list[Path] = []
    for path in candidates:
        resolved = path.resolve()
        if resolved in seen or not path.is_file():
            continue
        seen.add(resolved)
        unique.append(path)
    return unique


def contains_needle(path: Path, needle: str) -> bool:
    try:
        text = read_gzip_file(path) if path.suffix == ".gz" else read_text_file(path)
    except OSError:
        return False
    return needle in text


def compare_output_text(root: Path, config: dict, args: argparse.Namespace) -> str:
    if args.compare_output:
        output_path = Path(args.compare_output)
        if not output_path.is_file():
            raise GuardError(f"compare output not found: {output_path}")
        return read_text_file(output_path)

    if args.no_compare:
        return ""

    compare_script = config_path(config, "compare_script")
    native_run_dir = Path(args.native_run_dir) if args.native_run_dir else config_path(config, "pristine_baseline")
    if compare_script is None or native_run_dir is None:
        raise GuardError("compare script and native run dir are required")
    if not compare_script.is_file():
        raise GuardError(f"compare script not found: {compare_script}")
    if not native_run_dir.is_dir():
        raise GuardError(f"native run dir not found: {native_run_dir}")

    result = run_cmd(
        [
            sys.executable,
            str(compare_script),
            "--native-run-dir",
            str(native_run_dir),
            "--codon-run-dir",
            str(args.run_dir),
        ],
        root,
    )
    output = result.stdout + ("\n" + result.stderr if result.stderr else "")
    if result.returncode != 0:
        raise GuardError(f"compare command failed ({result.returncode}):\n{tail(output)}")
    return output


def tail(text: str, lines: int = 80) -> str:
    return "\n".join(text.splitlines()[-lines:])


def parse_char_diff_vars(compare_text: str) -> set[str]:
    vars_seen: set[str] = set()
    for match in re.finditer(r"^char_diff_vars=(.*)$", compare_text, flags=re.MULTILINE):
        value = match.group(1).strip()
        if not value:
            continue
        vars_seen.update(item.strip() for item in value.split(",") if item.strip())
    return vars_seen


def write_receipt(
    root: Path,
    config: dict,
    args: argparse.Namespace,
    *,
    run_env: Path,
    proof_count: int,
    compare_text: str,
) -> tuple[Path | None, str]:
    files = staged_files(root)
    if not files:
        return None, "no staged diff"
    if not args.selector or not args.proof_line:
        return None, "missing selector or proof-line"
    if "overall_numeric_equal=True" not in compare_text:
        return None, "compare did not report overall_numeric_equal=True"

    validation_mtime = run_env.stat().st_mtime
    newer_files = [
        relpath
        for relpath in files
        if (root / relpath).is_file() and (root / relpath).stat().st_mtime > validation_mtime
    ]
    if newer_files:
        return (
            None,
            "staged source file(s) are newer than run_environment for this job: "
            + ", ".join(newer_files[:10]),
        )

    diff_hash = staged_diff_hash(root)
    path = receipt_path(root, config, diff_hash)
    path.parent.mkdir(parents=True, exist_ok=True)
    now = _datetime.datetime.now(_datetime.timezone.utc).isoformat()
    payload = {
        "version": 1,
        "created_at_utc": now,
        "git_head": current_head(root),
        "staged_diff_hash": diff_hash,
        "staged_files": files,
        "job": args.job,
        "run_dir": str(Path(args.run_dir).resolve()),
        "case_root": str((Path(args.case_root) if args.case_root else config_path(config, "validation_case")).resolve()),
        "run_environment": str(run_env.resolve()),
        "selectors": args.selector,
        "proof_lines": args.proof_line,
        "proof_files": args.proof_file,
        "proof_candidates_checked": proof_count,
        "compare_output": str(Path(args.compare_output).resolve()) if args.compare_output else None,
        "native_run_dir": str(Path(args.native_run_dir).resolve()) if args.native_run_dir else str(config_path(config, "pristine_baseline").resolve()),
        "overall_numeric_equal": True,
        "char_diff_vars": sorted(parse_char_diff_vars(compare_text)),
    }
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return path, "written"


def validate_run(args: argparse.Namespace) -> int:
    root = repo_root()
    config = load_config(root, args.config)
    case_root = Path(args.case_root) if args.case_root else config_path(config, "validation_case")
    run_dir = Path(args.run_dir)
    if case_root is None or not case_root.is_dir():
        raise GuardError(f"case root not found: {case_root}")
    if not run_dir.is_dir():
        raise GuardError(f"run dir not found: {run_dir}")

    failures: list[str] = []
    run_env = find_run_environment(case_root, args.job)
    env_lines = set(read_text_file(run_env).splitlines())
    for selector in args.selector:
        if selector not in env_lines:
            failures.append(f"missing selector in {run_env}: {selector}")

    candidates = proof_candidates(run_dir, args.proof_file)
    for needle in args.proof_line:
        if not any(contains_needle(path, needle) for path in candidates):
            failures.append(f"missing proof line in atm/proof files: {needle}")

    compare_text = compare_output_text(root, config, args)
    if compare_text:
        if "overall_numeric_equal=True" not in compare_text:
            failures.append("compare output does not contain overall_numeric_equal=True")
        allowed = set(config_list(config, "allowed_metadata_char_diff_vars"))
        unexpected = parse_char_diff_vars(compare_text) - allowed
        if unexpected:
            failures.append("unexpected char diff vars: " + ",".join(sorted(unexpected)))

    if failures:
        print("CAM/Codon guard failed run validation:", file=sys.stderr)
        for failure in failures:
            print(f"  {failure}", file=sys.stderr)
        return 1

    receipt, receipt_reason = write_receipt(
        root,
        config,
        args,
        run_env=run_env,
        proof_count=len(candidates),
        compare_text=compare_text,
    ) if compare_text else (None, "compare skipped")

    print(f"CAM/Codon guard validate-run passed for job {args.job}.")
    print(f"run_environment={run_env}")
    print(f"proof_candidates={len(candidates)}")
    if compare_text:
        print("overall_numeric_equal=True")
    if receipt:
        print(f"receipt={receipt}")
    elif compare_text:
        print(f"receipt=<not written: {receipt_reason}>")
    return 0


def doctor(args: argparse.Namespace) -> int:
    root = repo_root()
    config_path_arg = Path(args.config or DEFAULT_CONFIG)
    if not config_path_arg.is_absolute():
        config_path_arg = root / config_path_arg
    config = load_config(root, args.config)
    hook_path = root / ".githooks" / "pre-commit"
    hooks = run_cmd(["git", "config", "--get", "core.hooksPath"], root)
    staged = staged_files(root)
    warnings = cache_warnings(root)
    diff_hash = staged_diff_hash(root)
    receipt = load_matching_receipt(root, config, diff_hash) if staged else None

    print(f"repo_root={root}")
    print(f"config={config_path_arg}")
    print(f"config_version={config.get('version')}")
    print(f"receipt_dir={receipt_dir(root, config)}")
    print(f"pre_commit_hook={hook_path}")
    print(f"pre_commit_hook_exists={hook_path.is_file()}")
    print(f"pre_commit_hook_executable={os.access(hook_path, os.X_OK)}")
    print(f"core.hooksPath={hooks.stdout.strip() or '<unset>'}")
    print(f"staged_count={len(staged)}")
    print(f"staged_diff_hash={diff_hash}")
    print(f"matching_receipt={'yes' if receipt else 'no'}")
    print(f"cache_warning_count={len(warnings)}")
    for path in warnings[:20]:
        print(f"cache_warning={path}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Guard checks for CAM/Codon migration work.")
    parser.add_argument("--config", default=DEFAULT_CONFIG, help="guard YAML config path")
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("pre-commit", help="run static staged-file checks")
    subparsers.add_parser("doctor", help="print guard configuration and hook status")

    validate = subparsers.add_parser("validate-run", help="validate a completed CAM/Codon run")
    validate.add_argument("--job", required=True, help="PBS job id prefix, e.g. 6082642")
    validate.add_argument("--run-dir", required=True, help="Codon run directory")
    validate.add_argument("--case-root", help="case root; defaults to config paths.validation_case")
    validate.add_argument("--native-run-dir", help="native baseline; defaults to config paths.pristine_baseline")
    validate.add_argument("--compare-output", help="existing compare output to inspect instead of running compare")
    validate.add_argument("--no-compare", action="store_true", help="skip numeric compare check")
    validate.add_argument("--selector", action="append", default=[], help="required KEY=value line in run_environment")
    validate.add_argument("--proof-line", action="append", default=[], help="required substring in atm.log/proof files")
    validate.add_argument("--proof-file", action="append", default=[], help="additional proof file to search")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "pre-commit":
            return pre_commit(args)
        if args.command == "validate-run":
            return validate_run(args)
        if args.command == "doctor":
            return doctor(args)
        raise GuardError(f"unknown command: {args.command}")
    except GuardError as exc:
        print(f"CAM/Codon guard error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
