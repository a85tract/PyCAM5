#!/usr/bin/env python3
"""Unit tests for the tphysac/tphysbc decoupling progress tracker."""

from __future__ import annotations

import copy
import importlib.util
from pathlib import Path
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "tphys_decoupling_progress.py"
SPEC = importlib.util.spec_from_file_location("tphys_decoupling_progress", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
tracker = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(tracker)


class ProgressTrackerTests(unittest.TestCase):
    def setUp(self) -> None:
        self.data = tracker.load_status(tracker.DEFAULT_STATUS_FILE)
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def _file(self, name: str, content: str = "evidence\n") -> Path:
        path = self.root / name
        path.write_text(content, encoding="utf-8")
        return path

    def _bfb_run(self, batch: str = "B01") -> dict:
        executable = self._file("cesm.exe")
        environment = self._file("run_environment.txt")
        filepath = self._file("Filepath")
        namelist = self._file("atm_in")
        compare = self._file("compare.txt", "overall_numeric_equal=True\noverall_char_equal=True\n")
        output = self._file("history.nc")
        proof_file = self._file("atm.log", "all extracted entries executed\n")
        proof_source = tracker.artifact_from_path(str(proof_file))
        proof = {
            item["entry_point"]: {
                "source": copy.deepcopy(proof_source),
                "line": f"{item['entry_point']} calls = 50",
            }
            for item in tracker.processes_for_batch(self.data, batch)
        }
        return {
            "run_id": f"{batch}-test-run",
            "batch": batch,
            "result": "bfb",
            "steps": 50,
            "job_id": "123.desched1",
            "source_commit": self.data["baseline_commit"],
            "case_root": str(self.root / "case"),
            "run_dir": str(self.root / "run"),
            "baseline_run_dir": str(self.root / "native"),
            "numeric_equal": True,
            "char_equal": True,
            "artifacts": {
                "executable": tracker.artifact_from_path(str(executable)),
                "run_environment": tracker.artifact_from_path(str(environment)),
                "filepath": tracker.artifact_from_path(str(filepath)),
                "namelist": tracker.artifact_from_path(str(namelist)),
                "compare_report": tracker.artifact_from_path(str(compare)),
                "output_files": [tracker.output_from_path(str(output))],
            },
            "execution_proof": proof,
            "recorded_at": "2026-07-15T00:00:00Z",
            "note": "test evidence",
        }

    def test_seed_inventory_and_batch_sizes(self) -> None:
        tracker.assert_valid_status(self.data)
        self.assertEqual(len(self.data["processes"]), 69)
        counts = {
            batch_id: len(tracker.processes_for_batch(self.data, batch_id))
            for batch_id in tracker.EXPECTED_BATCH_SIZES
        }
        self.assertEqual(counts, tracker.EXPECTED_BATCH_SIZES)
        self.assertIn(
            "set_dry_to_wet_run",
            {item["entry_point"] for item in tracker.processes_for_batch(self.data, "B02")},
        )

    def test_suite_files_cover_exact_active_inventory(self) -> None:
        self.assertEqual(
            tracker.collect_suite_errors(self.data, tracker.DEFAULT_SUITE_DIR),
            [],
        )

    def test_render_is_deterministic_and_reports_current_counts(self) -> None:
        first = tracker.render_markdown(self.data)
        second = tracker.render_markdown(self.data)
        self.assertEqual(first, second)
        statuses = [item["status"] for item in self.data["processes"]]
        expected_summary = (
            f"| 69 | {statuses.count('bfb')} | {statuses.count('in_progress')} | "
            f"{statuses.count('build_pass')} | {statuses.count('50step_running')} | "
            f"{statuses.count('failed')} | {statuses.count('planned')} |"
        )
        self.assertIn(expected_summary, first)
        self.assertIn("| B01 | 公共状态更新 | 9 |", first)
        self.assertIn("| B07 | 云、微物理与辐射 | 10 |", first)

    def test_manual_bfb_update_is_rejected_without_evidence(self) -> None:
        with self.assertRaisesRegex(tracker.StatusError, "use record-run"):
            tracker.update_processes(
                self.data,
                ["dadadj_run"],
                "bfb",
            )

    def test_bfb_rejects_empty_output_file_set(self) -> None:
        run = self._bfb_run()
        run["artifacts"]["output_files"] = []
        with self.assertRaisesRegex(tracker.StatusError, "empty output file set"):
            tracker.record_run(self.data, run)

    def test_bfb_rejects_nontrue_numeric_or_char_result(self) -> None:
        for field in ("numeric_equal", "char_equal"):
            with self.subTest(field=field):
                run = self._bfb_run()
                run[field] = False
                with self.assertRaisesRegex(tracker.StatusError, field):
                    tracker.record_run(self.data, run)

    def test_bfb_rejects_missing_execution_proof(self) -> None:
        run = self._bfb_run()
        missing = next(iter(run["execution_proof"]))
        del run["execution_proof"][missing]
        with self.assertRaisesRegex(tracker.StatusError, "lacks execution proof"):
            tracker.record_run(self.data, run)

    def test_complete_bfb_run_marks_entire_batch(self) -> None:
        run = self._bfb_run("B01")
        updated = tracker.record_run(self.data, run)
        members = tracker.processes_for_batch(updated, "B01")
        self.assertTrue(all(item["status"] == "bfb" for item in members))
        self.assertTrue(all(item["bfb_run"] == run["run_id"] for item in members))
        self.assertTrue(all(item["commit"] == run["source_commit"] for item in members))
        tracker.assert_valid_status(updated)

    def test_failed_run_is_visible_but_not_bfb(self) -> None:
        run = {
            "run_id": "B02-failed-test",
            "batch": "B02",
            "result": "failed",
            "steps": 50,
            "job_id": "124.desched1",
            "source_commit": self.data["baseline_commit"],
            "case_root": str(self.root / "case"),
            "run_dir": str(self.root / "run"),
            "baseline_run_dir": str(self.root / "native"),
            "numeric_equal": False,
            "char_equal": True,
            "artifacts": {"output_files": []},
            "execution_proof": {},
            "recorded_at": "2026-07-15T00:00:00Z",
            "note": "numeric comparison failed",
        }
        updated = tracker.record_run(self.data, run)
        members = tracker.processes_for_batch(updated, "B02")
        self.assertTrue(all(item["status"] == "failed" for item in members))
        self.assertTrue(all(item["bfb_run"] is None for item in members))
        tracker.assert_valid_status(updated)


if __name__ == "__main__":
    unittest.main()
