#!/usr/bin/env python3
"""Audit direct and transitive Fortran module dependencies for atmos_phys.

The standalone boundary is deliberately strict: every non-intrinsic module
must either be implemented under ``src/atmos_phys`` or belong to the small
external runtime allowlist (ISO/IEEE, ``shr_*``, MPI, OpenMP, NetCDF, PIO).
Modules found elsewhere in CAM are host dependencies and are rejected.

``--scope completed-only`` checks ledger entries at ``dependency_pass`` or a
later state.  ``--scope all`` is the 36-entry closure gate.  Duplicate module
definitions under ``src/atmos_phys`` are always rejected because they make
Fortran module resolution dependent on source ordering.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import re
import sys
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_STATUS_FILE = REPO_ROOT / "src/atmos_phys/standalone/standalone_status.json"
FORTRAN_SUFFIXES = {".f", ".for", ".ftn", ".f90", ".f95", ".f03", ".f08"}
DEFAULT_COMPLETED_STATUSES = {"dependency_pass", "build_pass", "50step_running", "bfb"}

MODULE_DEFINITION_RE = re.compile(
    r"^\s*module\s+(?!procedure\b|subroutine\b|function\b)([a-z_]\w*)",
    re.IGNORECASE,
)
USE_RE = re.compile(
    r"^\s*use(?:\s*,\s*((?:non_)?intrinsic))?\s*(?:::\s*)?([a-z_]\w*)",
    re.IGNORECASE,
)


class DependencyError(ValueError):
    """Raised for malformed inputs rather than dependency violations."""


@dataclass(frozen=True)
class UseDependency:
    module: str
    intrinsic: bool
    line: int


@dataclass(frozen=True)
class SourceInfo:
    path: Path
    modules: frozenset[str]
    uses: tuple[UseDependency, ...]


@dataclass(frozen=True)
class Violation:
    mode: str
    entry_point: str
    source: Path
    module: str
    line: int
    reason: str
    chain: tuple[str, ...]

    def format(self, repo_root: Path) -> str:
        try:
            source = self.source.relative_to(repo_root)
        except ValueError:
            source = self.source
        chain = " -> ".join(self.chain)
        return (
            f"{self.mode}: {self.entry_point}: {source}:{self.line}: use {self.module}: "
            f"{self.reason}; chain={chain}"
        )


@dataclass
class AuditResult:
    mode: str
    selected_entries: tuple[str, ...]
    root_sources: tuple[Path, ...]
    closure_sources: tuple[Path, ...]
    violations: tuple[Violation, ...]
    duplicate_modules: dict[str, tuple[Path, ...]]

    @property
    def passed(self) -> bool:
        return not self.violations and not self.duplicate_modules


def strip_fortran_comment(line: str) -> str:
    """Strip a free-form Fortran comment while respecting quoted strings."""
    quote: str | None = None
    index = 0
    while index < len(line):
        char = line[index]
        if quote is not None:
            if char == quote:
                if index + 1 < len(line) and line[index + 1] == quote:
                    index += 2
                    continue
                quote = None
        elif char in ("'", '"'):
            quote = char
        elif char == "!":
            return line[:index]
        index += 1
    return line


def split_fortran_statements(line: str) -> list[str]:
    """Split semicolon-separated statements while respecting strings."""
    statements: list[str] = []
    start = 0
    quote: str | None = None
    index = 0
    while index < len(line):
        char = line[index]
        if quote is not None:
            if char == quote:
                if index + 1 < len(line) and line[index + 1] == quote:
                    index += 2
                    continue
                quote = None
        elif char in ("'", '"'):
            quote = char
        elif char == ";":
            statements.append(line[start:index])
            start = index + 1
        index += 1
    statements.append(line[start:])
    return statements


def logical_fortran_lines(text: str) -> list[tuple[int, str]]:
    """Return joined free-form logical lines with their first physical line."""
    result: list[tuple[int, str]] = []
    pending = ""
    pending_line = 0
    continuing = False
    for number, raw_line in enumerate(text.splitlines(), start=1):
        clean = strip_fortran_comment(raw_line).rstrip()
        if not clean.strip() or clean.lstrip().startswith("#"):
            continue
        part = clean.lstrip()
        if continuing and part.startswith("&"):
            part = part[1:].lstrip()
        has_continuation = part.endswith("&")
        if has_continuation:
            part = part[:-1].rstrip()
        if not pending:
            pending_line = number
            pending = part
        else:
            pending = f"{pending} {part}"
        continuing = has_continuation
        if not continuing:
            for statement in split_fortran_statements(pending):
                if statement.strip():
                    result.append((pending_line, statement.strip()))
            pending = ""
    if pending:
        for statement in split_fortran_statements(pending):
            if statement.strip():
                result.append((pending_line, statement.strip()))
    return result


def parse_source(path: Path) -> SourceInfo:
    text = path.read_text(encoding="utf-8", errors="replace")
    modules: set[str] = set()
    uses: list[UseDependency] = []
    for line_number, statement in logical_fortran_lines(text):
        module_match = MODULE_DEFINITION_RE.match(statement)
        if module_match:
            modules.add(module_match.group(1).lower())
        use_match = USE_RE.match(statement)
        if use_match:
            nature = (use_match.group(1) or "").lower()
            uses.append(
                UseDependency(
                    module=use_match.group(2).lower(),
                    intrinsic=nature == "intrinsic",
                    line=line_number,
                )
            )
    return SourceInfo(path=path.resolve(), modules=frozenset(modules), uses=tuple(uses))


def fortran_files(root: Path) -> list[Path]:
    if not root.is_dir():
        return []
    return sorted(
        path
        for path in root.rglob("*")
        if path.is_file() and path.suffix.lower() in FORTRAN_SUFFIXES
    )


def build_source_index(root: Path) -> dict[Path, SourceInfo]:
    return {path.resolve(): parse_source(path) for path in fortran_files(root)}


def build_module_index(sources: Iterable[SourceInfo]) -> dict[str, list[Path]]:
    index: dict[str, list[Path]] = {}
    for source in sources:
        for module in source.modules:
            index.setdefault(module, []).append(source.path)
    for providers in index.values():
        providers.sort()
    return index


def is_allowed_external_module(module: str, intrinsic: bool = False) -> bool:
    """Return whether *module* belongs to the explicit standalone allowlist."""
    name = module.lower()
    if name.startswith(("iso_", "ieee_", "shr_", "omp_", "netcdf")):
        return True
    if name == "mpi" or name.startswith("mpi_"):
        return True
    if name == "pio" or name.startswith("pio_"):
        return True
    # Marking an arbitrary module intrinsic must not bypass the allowlist.
    return intrinsic and name in {"iso_fortran_env", "iso_c_binding", "ieee_arithmetic"}


def load_ledger(path: Path) -> dict:
    try:
        with path.open("r", encoding="utf-8") as handle:
            data = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise DependencyError(f"cannot read standalone ledger {path}: {exc}") from exc
    if not isinstance(data, dict) or not isinstance(data.get("processes"), list):
        raise DependencyError("standalone ledger must contain a process list")
    return data


def select_processes(
    data: dict,
    scope: str,
    batches: set[str] | None = None,
    entries: set[str] | None = None,
) -> list[dict]:
    processes = list(data["processes"])
    if scope == "completed-only":
        completed = set(data.get("validation", {}).get("completed_statuses", ()))
        if not completed:
            completed = DEFAULT_COMPLETED_STATUSES
        processes = [item for item in processes if item.get("status") in completed]
    elif scope != "all":
        raise DependencyError(f"unknown scope: {scope}")
    if batches:
        processes = [item for item in processes if item.get("batch") in batches]
    if entries:
        known = {item.get("entry_point") for item in data["processes"]}
        unknown = sorted(entries.difference(known))
        if unknown:
            raise DependencyError(f"unknown entries: {', '.join(unknown)}")
        processes = [item for item in processes if item.get("entry_point") in entries]
    return sorted(processes, key=lambda item: (item.get("batch", ""), item.get("order", 0)))


def outside_module_index(repo_root: Path, atmos_root: Path) -> dict[str, list[Path]]:
    sources: list[SourceInfo] = []
    src_root = repo_root / "src"
    for path in fortran_files(src_root):
        resolved = path.resolve()
        try:
            resolved.relative_to(atmos_root)
            continue
        except ValueError:
            pass
        sources.append(parse_source(resolved))
    return build_module_index(sources)


def audit_dependencies(
    repo_root: Path,
    data: dict,
    processes: list[dict],
    mode: str,
    *,
    internal_sources: dict[Path, SourceInfo] | None = None,
    outside_modules: dict[str, list[Path]] | None = None,
) -> AuditResult:
    if mode not in {"direct", "transitive"}:
        raise DependencyError(f"invalid audit mode: {mode}")
    repo_root = repo_root.resolve()
    atmos_root = (repo_root / "src/atmos_phys").resolve()
    source_index = internal_sources or build_source_index(atmos_root)
    module_index = build_module_index(source_index.values())
    duplicates = {
        module: tuple(paths)
        for module, paths in module_index.items()
        if len(set(paths)) > 1
    }
    outside = outside_modules if outside_modules is not None else outside_module_index(repo_root, atmos_root)

    root_sources: list[Path] = []
    source_entries: dict[Path, list[str]] = {}
    for process in processes:
        source_raw = process.get("source")
        if not isinstance(source_raw, str):
            raise DependencyError(f"{process.get('entry_point')} has no source path")
        source = (repo_root / source_raw).resolve()
        try:
            source.relative_to(atmos_root)
        except ValueError as exc:
            raise DependencyError(f"source escapes src/atmos_phys: {source_raw}") from exc
        if source not in source_index:
            raise DependencyError(f"ledger source is missing or not Fortran: {source_raw}")
        if source not in source_entries:
            root_sources.append(source)
        source_entries.setdefault(source, []).append(process["entry_point"])

    violations: list[Violation] = []
    closure: set[Path] = set(root_sources)
    for root_source in root_sources:
        for entry in source_entries[root_source]:
            queue: list[tuple[Path, tuple[str, ...]]] = [
                (root_source, (entry, root_source.name))
            ]
            visited: set[Path] = set()
            while queue:
                source_path, chain = queue.pop(0)
                if source_path in visited:
                    continue
                visited.add(source_path)
                closure.add(source_path)
                source = source_index[source_path]
                for dependency in source.uses:
                    if is_allowed_external_module(dependency.module, dependency.intrinsic):
                        continue
                    providers = module_index.get(dependency.module, [])
                    if len(providers) == 1:
                        provider = providers[0]
                        closure.add(provider)
                        if mode == "transitive" and provider not in visited:
                            queue.append((provider, chain + (dependency.module, provider.name)))
                        continue
                    if len(providers) > 1:
                        # The global duplicate diagnostic is more useful than one use-site error.
                        continue
                    external_providers = outside.get(dependency.module, [])
                    if external_providers:
                        rendered = ", ".join(
                            str(path.relative_to(repo_root))
                            if path.is_relative_to(repo_root)
                            else str(path)
                            for path in external_providers[:4]
                        )
                        reason = f"CAM module outside src/atmos_phys ({rendered})"
                    else:
                        reason = "unresolved non-allowlisted module"
                    violations.append(
                        Violation(
                            mode=mode,
                            entry_point=entry,
                            source=source_path,
                            module=dependency.module,
                            line=dependency.line,
                            reason=reason,
                            chain=chain + (dependency.module,),
                        )
                    )
                if mode == "direct":
                    break

    unique_violations = {
        (
            item.mode,
            item.entry_point,
            item.source,
            item.module,
            item.line,
            item.reason,
            item.chain,
        ): item
        for item in violations
    }
    return AuditResult(
        mode=mode,
        selected_entries=tuple(item["entry_point"] for item in processes),
        root_sources=tuple(sorted(root_sources)),
        closure_sources=tuple(sorted(closure)),
        violations=tuple(
            sorted(
                unique_violations.values(),
                key=lambda item: (item.entry_point, str(item.source), item.line, item.module),
            )
        ),
        duplicate_modules=duplicates,
    )


def render_result(result: AuditResult, repo_root: Path) -> list[str]:
    status = "PASS" if result.passed else "FAIL"
    lines = [
        (
            f"[{result.mode}] {status}: entries={len(result.selected_entries)} "
            f"roots={len(result.root_sources)} closure_sources={len(result.closure_sources)} "
            f"violations={len(result.violations)} duplicates={len(result.duplicate_modules)}"
        )
    ]
    for module, providers in sorted(result.duplicate_modules.items()):
        rendered = []
        for provider in providers:
            try:
                rendered.append(str(provider.relative_to(repo_root)))
            except ValueError:
                rendered.append(str(provider))
        lines.append(f"duplicate module {module}: {', '.join(rendered)}")
    lines.extend(item.format(repo_root) for item in result.violations)
    return lines


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=str(REPO_ROOT))
    parser.add_argument("--status-file", default=str(DEFAULT_STATUS_FILE))
    parser.add_argument("--scope", choices=("completed-only", "all"), default="completed-only")
    parser.add_argument("--mode", choices=("direct", "transitive", "both"), default="both")
    parser.add_argument("--batch", action="append", choices=("S01", "S02", "S03", "S04"))
    parser.add_argument("--process", action="append")
    parser.add_argument(
        "--list-sources",
        action="store_true",
        help="print the transitive, deduplicated source closure only; fail on violations",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        repo_root = Path(args.repo_root).resolve()
        data = load_ledger(Path(args.status_file).resolve())
        processes = select_processes(
            data,
            args.scope,
            set(args.batch or ()),
            set(args.process or ()),
        )
        atmos_root = repo_root / "src/atmos_phys"
        source_index = build_source_index(atmos_root)
        outside = outside_module_index(repo_root, atmos_root.resolve())
        modes = ("direct", "transitive") if args.mode == "both" else (args.mode,)
        results = [
            audit_dependencies(
                repo_root,
                data,
                processes,
                mode,
                internal_sources=source_index,
                outside_modules=outside,
            )
            for mode in modes
        ]
        passed = all(result.passed for result in results)
        if args.list_sources:
            if args.mode == "direct":
                raise DependencyError("--list-sources requires transitive or both mode")
            if not passed:
                for result in results:
                    print("\n".join(render_result(result, repo_root)), file=sys.stderr)
                return 1
            transitive = next(result for result in results if result.mode == "transitive")
            for source in transitive.closure_sources:
                print(source)
            return 0
        print(
            f"scope={args.scope} selected_entries={len(processes)} "
            f"allowed_external=ISO/IEEE,shr_*,MPI,OpenMP,NetCDF,PIO"
        )
        for result in results:
            print("\n".join(render_result(result, repo_root)))
        return 0 if passed else 1
    except DependencyError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
