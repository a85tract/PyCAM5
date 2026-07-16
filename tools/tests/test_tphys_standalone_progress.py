#!/usr/bin/env python3
"""Tests for the independent 36-entry standalone progress tracker."""

from __future__ import annotations

import copy
import importlib.util
from pathlib import Path
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "tphys_standalone_progress.py"
SPEC = importlib.util.spec_from_file_location("tphys_standalone_progress", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
tracker = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(tracker)


class StandaloneProgressTests(unittest.TestCase):
    def setUp(self) -> None:
        self.data = tracker.load_status(tracker.DEFAULT_STATUS_FILE)
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def evidence(self, name: str, text: str = "pass\n") -> Path:
        path = self.root / name
        path.write_text(text, encoding="utf-8")
        return path

    def bfb_run(self, batch: str = "S01") -> dict:
        artifacts = {
            name: tracker.artifact_from_path(str(self.evidence(name)))
            for name in tracker.REQUIRED_BFB_ARTIFACTS
        }
        artifacts["output_files"] = [
            tracker.output_from_path(str(self.evidence("history.nc")))
        ]
        proof_file = self.evidence("timing.txt", "timer global_calls=50\n")
        proof_source = tracker.artifact_from_path(str(proof_file))
        proof = {
            item["entry_point"]: {
                "source": copy.deepcopy(proof_source),
                "line": f"ap_{item['entry_point']} global_calls=50",
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
            "artifacts": artifacts,
            "execution_proof": proof,
            "recorded_at": "2026-07-16T00:00:00Z",
            "note": "unit-test evidence",
        }

    def test_inventory_is_exactly_four_groups_of_nine(self) -> None:
        tracker.assert_valid(self.data)
        self.assertEqual(len(self.data["processes"]), 36)
        self.assertEqual(
            {
                batch: tuple(
                    item["entry_point"]
                    for item in sorted(
                        tracker.processes_for_batch(self.data, batch),
                        key=lambda item: item["order"],
                    )
                )
                for batch in tracker.EXPECTED_BATCHES
            },
            tracker.EXPECTED_BATCHES,
        )

    def test_generated_markdown_is_current_and_deterministic(self) -> None:
        rendered = tracker.render_markdown(self.data)
        self.assertEqual(rendered, tracker.render_markdown(self.data))
        self.assertEqual(rendered, tracker.DEFAULT_PROGRESS_FILE.read_text(encoding="utf-8"))
        self.assertIn("| 36 | 0 | 0 | 0 | 0 | 0 | 0 | 36 |", rendered)

    def test_phase_one_tracker_is_a_different_file(self) -> None:
        phase_one = tracker.REPO_ROOT / "tools/tphys_decoupling_progress.py"
        self.assertTrue(phase_one.is_file())
        self.assertNotEqual(phase_one.resolve(), SCRIPT.resolve())

    def test_update_cannot_mark_bfb_without_run_evidence(self) -> None:
        with self.assertRaisesRegex(tracker.StatusError, "record-run"):
            tracker.update_processes(
                self.data,
                ["compute_vdiff_run"],
                "bfb",
                self.data["baseline_commit"],
                None,
            )

    def test_update_dependency_pass_is_visible(self) -> None:
        updated = tracker.update_processes(
            self.data,
            ["compute_vdiff_run"],
            "dependency_pass",
            self.data["baseline_commit"],
            "direct and transitive checks passed",
        )
        process = tracker.process_index(updated)["compute_vdiff_run"]
        self.assertEqual(process["status"], "dependency_pass")
        self.assertIn("direct and transitive", tracker.render_markdown(updated))

    def test_bfb_requires_standalone_reports(self) -> None:
        run = self.bfb_run()
        del run["artifacts"]["standalone_build_report"]
        with self.assertRaisesRegex(tracker.StatusError, "standalone_build_report"):
            tracker.record_run(self.data, run)

    def test_partial_run_marks_only_proven_entry(self) -> None:
        run = self.bfb_run("S01")
        proven = next(iter(run["execution_proof"]))
        run["execution_proof"] = {proven: run["execution_proof"][proven]}
        updated = tracker.record_run(self.data, run)
        by_entry = tracker.process_index(updated)
        self.assertEqual(by_entry[proven]["status"], "bfb")
        self.assertTrue(
            all(
                item["status"] == "planned"
                for item in tracker.processes_for_batch(updated, "S01")
                if item["entry_point"] != proven
            )
        )
        tracker.assert_valid(updated)

    def test_full_bfb_run_marks_nine_entries(self) -> None:
        run = self.bfb_run("S04")
        updated = tracker.record_run(self.data, run)
        members = tracker.processes_for_batch(updated, "S04")
        self.assertTrue(all(item["status"] == "bfb" for item in members))
        self.assertTrue(all(item["bfb_run"] == run["run_id"] for item in members))


if __name__ == "__main__":
    unittest.main()
