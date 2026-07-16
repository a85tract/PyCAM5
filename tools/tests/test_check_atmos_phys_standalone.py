#!/usr/bin/env python3
"""Tests for direct/transitive standalone Fortran dependency auditing."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "check_atmos_phys_standalone.py"
SPEC = importlib.util.spec_from_file_location("check_atmos_phys_standalone", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
checker = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = checker
SPEC.loader.exec_module(checker)


class DependencyCheckerTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        (self.root / "src/atmos_phys/schemes").mkdir(parents=True)
        (self.root / "src/physics/cam").mkdir(parents=True)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def write(self, relative: str, text: str) -> Path:
        path = self.root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
        return path

    def ledger(self, processes: list[dict]) -> dict:
        return {
            "validation": {
                "completed_statuses": [
                    "dependency_pass",
                    "build_pass",
                    "50step_running",
                    "bfb",
                ]
            },
            "processes": processes,
        }

    def process(self, entry: str, source: str, status: str = "planned") -> dict:
        return {
            "entry_point": entry,
            "source": source,
            "batch": "S01",
            "order": 1,
            "status": status,
        }

    def audit(self, data: dict, processes: list[dict], mode: str):
        atmos = self.root / "src/atmos_phys"
        return checker.audit_dependencies(
            self.root,
            data,
            processes,
            mode,
            internal_sources=checker.build_source_index(atmos),
            outside_modules=checker.outside_module_index(self.root, atmos.resolve()),
        )

    def test_parser_handles_continuations_semicolons_and_comments(self) -> None:
        source = self.write(
            "src/atmos_phys/schemes/entry.F90",
            """
module entry_mod
  use, intrinsic :: iso_fortran_env, only: real64
  use shr_kind_mod, only: r8 => shr_kind_r8 ! comment
  use science_mod, &
       only: kernel
contains; subroutine entry_run(); end subroutine entry_run
end module entry_mod
""",
        )
        parsed = checker.parse_source(source)
        self.assertEqual(parsed.modules, {"entry_mod"})
        self.assertEqual(
            {dependency.module for dependency in parsed.uses},
            {"iso_fortran_env", "shr_kind_mod", "science_mod"},
        )

    def test_external_allowlist_is_exactly_supported_families(self) -> None:
        for module in (
            "iso_c_binding",
            "ieee_arithmetic",
            "shr_kind_mod",
            "mpi_f08",
            "omp_lib",
            "netcdf",
            "pio_types",
        ):
            with self.subTest(module=module):
                self.assertTrue(checker.is_allowed_external_module(module))
        self.assertFalse(checker.is_allowed_external_module("physics_types"))
        self.assertFalse(checker.is_allowed_external_module("mpishorthand"))
        self.assertFalse(checker.is_allowed_external_module("random_vendor_module"))

    def test_direct_passes_internal_science_but_transitive_rejects_host(self) -> None:
        root_source = self.write(
            "src/atmos_phys/schemes/entry.F90",
            "module entry_mod\n use science_mod\ncontains\n subroutine entry_run()\n end subroutine\nend module\n",
        )
        science_source = self.write(
            "src/atmos_phys/support/science.F90",
            "module science_mod\n use physics_types\nend module science_mod\n",
        )
        self.write(
            "src/physics/cam/physics_types.F90",
            "module physics_types\nend module physics_types\n",
        )
        process = self.process(
            "entry_run", str(root_source.relative_to(self.root)), "dependency_pass"
        )
        data = self.ledger([process])
        direct = self.audit(data, [process], "direct")
        transitive = self.audit(data, [process], "transitive")
        self.assertTrue(direct.passed)
        self.assertFalse(transitive.passed)
        self.assertIn(science_source.resolve(), transitive.closure_sources)
        self.assertEqual(transitive.violations[0].module, "physics_types")
        self.assertIn("outside src/atmos_phys", transitive.violations[0].reason)

    def test_unresolved_nonallowlisted_module_is_rejected(self) -> None:
        source = self.write(
            "src/atmos_phys/schemes/entry.F90",
            "module entry_mod\n use mystery_mod\ncontains\n subroutine entry_run()\n end subroutine\nend module\n",
        )
        process = self.process("entry_run", str(source.relative_to(self.root)))
        result = self.audit(self.ledger([process]), [process], "direct")
        self.assertFalse(result.passed)
        self.assertEqual(result.violations[0].reason, "unresolved non-allowlisted module")

    def test_duplicate_atmos_phys_module_is_always_rejected(self) -> None:
        source = self.write(
            "src/atmos_phys/schemes/entry.F90",
            "module duplicate_mod\ncontains\n subroutine entry_run()\n end subroutine\nend module\n",
        )
        self.write(
            "src/atmos_phys/support/duplicate.F90",
            "module duplicate_mod\nend module duplicate_mod\n",
        )
        process = self.process("entry_run", str(source.relative_to(self.root)))
        result = self.audit(self.ledger([process]), [process], "direct")
        self.assertFalse(result.passed)
        self.assertIn("duplicate_mod", result.duplicate_modules)

    def test_completed_only_excludes_planned_entries(self) -> None:
        completed = self.process("done_run", "src/atmos_phys/schemes/done.F90", "dependency_pass")
        planned = self.process("todo_run", "src/atmos_phys/schemes/todo.F90", "planned")
        data = self.ledger([completed, planned])
        self.assertEqual(
            [item["entry_point"] for item in checker.select_processes(data, "completed-only")],
            ["done_run"],
        )
        self.assertEqual(len(checker.select_processes(data, "all")), 2)

    def test_transitive_source_closure_is_deduplicated(self) -> None:
        first = self.write(
            "src/atmos_phys/schemes/first.F90",
            "module first_mod\n use shared_mod\ncontains\n subroutine first_run()\n end subroutine\nend module\n",
        )
        second = self.write(
            "src/atmos_phys/schemes/second.F90",
            "module second_mod\n use shared_mod\ncontains\n subroutine second_run()\n end subroutine\nend module\n",
        )
        shared = self.write(
            "src/atmos_phys/support/shared.F90",
            "module shared_mod\n use shr_kind_mod\nend module shared_mod\n",
        )
        processes = [
            self.process("first_run", str(first.relative_to(self.root))),
            self.process("second_run", str(second.relative_to(self.root))),
        ]
        result = self.audit(self.ledger(processes), processes, "transitive")
        self.assertTrue(result.passed)
        self.assertEqual(result.closure_sources.count(shared.resolve()), 1)


if __name__ == "__main__":
    unittest.main()
