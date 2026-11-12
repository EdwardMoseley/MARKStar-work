import csv
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import markstar_scan as scan
import markstar_sequential as sequential
from test_markstar_scan import atom


class SequentialTests(unittest.TestCase):
    def test_epsilon_default(self):
        self.assertEqual(sequential.parse_args([]).epsilon, 0.683)

    def test_select_specific_target_and_histidine_variants(self):
        residues = scan.read_residues(scan.DEFAULT_PDB)
        target, neighbors = scan.select_neighbors(residues, 4, 4, "B2")
        self.assertEqual(target.id, "B2")
        self.assertIn("B1", [r.id for _, r in neighbors])
        self.assertEqual(scan.canonical_identity("HID"), "HIS")
        with self.assertRaisesRegex(ValueError, "not found"):
            scan.select_neighbors(residues, 4, 4, "B999")

    def test_upper_bound_order_and_bad_values(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "scores.csv"
            path.write_text(
                "sequence,identity,lower_bound,upper_bound\n"
                "B1=VAL,VAL,1,1e+90\nB1=ALA,ALA,1,1e+90\n"
                "B1=HID,HID,0,NaN\nB1=GLY,GLY,0,-1\nB1=CYS,CYS,0,0\n"
                "B1=ARG,ARG,0,Infinity\nB1=SER,SER,0,-Infinity\n"
            )
            candidates, warnings = sequential.read_candidates(path, "B1")
            self.assertEqual(sequential.choose_candidate(candidates)["identity"], "ALA")
            self.assertEqual(len(warnings), 5)
            self.assertNotIn("ARG", [row["identity"] for row in candidates])
            self.assertNotIn("SER", [row["identity"] for row in candidates])
            self.assertNotIn("CYS", [row["identity"] for row in candidates])
            exact = [
                {"identity": "ALA", "upper_bound": "1.00000000000000000000000000001"},
                {"identity": "VAL", "upper_bound": "1.00000000000000000000000000002"},
            ]
            self.assertEqual(sequential.choose_candidate(exact)["identity"], "VAL")

    def test_extract_real_first_model(self):
        source = (scan.DEFAULT_PDB.parent / "markstar-19-with-structures" /
                  "structures" / "HID.pdb")
        if not source.exists():
            self.skipTest("Prior real export is unavailable")
        residues = scan.read_residues(scan.DEFAULT_PDB)
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "next.pdb"
            sequential.extract_first_model(source, output, residues, "B1", "HID")
            result = scan.read_residues(output)
            self.assertEqual([r.id for r in result], [r.id for r in residues])
            self.assertEqual(next(r.name for r in result if r.id == "B1"), "HID")
            self.assertNotIn("MODEL ", output.read_text())
            with self.assertRaisesRegex(ValueError, "unexpected identities"):
                sequential.extract_first_model(source, Path(tmp) / "bad.pdb", residues, "B1", "ALA")

    def test_workflow_retries_and_carries_forward(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            input_pdb = base / "input.pdb"
            input_pdb.write_text(
                atom(1, "CA", "ALA", "A", 1, 0, 0, 0)
                + atom(2, "CA", "HIS", "B", 1, 1, 0, 0)
                + atom(3, "CA", "GLU", "B", 2, 2, 0, 0)
                + "END\n"
            )
            calls = []

            def fake_attempt(args, current, target, directory, exclude):
                directory.mkdir()
                calls.append((target, exclude, str(current)))
                log = directory.parent / (directory.name + ".log")
                unsupported = target == "B1" and not exclude
                log.write_text("Mutation from HID to PRO not supported yet\n" if unsupported else "")
                info = {"directory": str(directory), "log": str(log),
                        "exit_code": 1 if unsupported else 0, "exclude_proline": exclude,
                        "proline_unsupported": unsupported,
                        "score_validation_only_failure": target == "B1" and exclude,
                        "manifest_status": "failed" if unsupported else "completed",
                        "score_problems": [], "warnings": []}
                if target == "B1" and exclude:
                    (directory / "sequence_scores.csv").write_text(
                        "sequence,identity,lower_bound,upper_bound\n"
                        "B1=ARG,ARG,0,Infinity\n"
                        "B1=GLY,GLY,0,30\n"
                        "B1=ALA,ALA,1,20\nB1=HID,HID,2,10\n"
                    )
                    structures = directory / "structures"
                    structures.mkdir()
                    (structures / "ALA.pdb").write_text(
                        "MODEL        1\n"
                        + atom(1, "CA", "ALA", "A", 1, 0, 0, 0)
                        + atom(2, "CA", "ALA", "B", 1, 1, 0, 0)
                        + atom(3, "CA", "GLU", "B", 2, 2, 0, 0)
                        + "ENDMDL\nEND\n"
                    )
                return info, []

            output = base / "workflow"
            with patch.object(sequential, "run_attempt", side_effect=fake_attempt):
                result = sequential.main(["--pdb", str(input_pdb), "--output", str(output)])
            self.assertEqual(result, 1)
            self.assertEqual([(t, e) for t, e, _ in calls],
                             [("B1", False), ("B1", True), ("B2", False)])
            with (output / "workflow_manifest.csv").open(newline="") as file:
                rows = list(csv.DictReader(file))
            self.assertEqual([r["status"] for r in rows], ["selected", "skipped"])
            self.assertEqual(rows[0]["selected_identity"], "ALA")
            self.assertEqual(rows[0]["translation_rotation_enabled"], "True")
            self.assertEqual(rows[0]["translation_rotation_max_rotation_degrees"], "10.0")
            self.assertEqual(rows[0]["translation_rotation_max_translation_angstrom"], "2.0")
            self.assertEqual(rows[0]["output_pdb"], rows[1]["input_pdb"])
            self.assertEqual(rows[1]["input_pdb"], rows[1]["output_pdb"])
            self.assertEqual(len(json.loads(rows[0]["attempts"])), 2)
            self.assertTrue(json.loads(rows[0]["attempts"])[1]["score_validation_only_failure"])
            self.assertNotIn("ARG", [c["identity"] for c in json.loads(rows[0]["candidate_scores"])])
            self.assertIn("non-finite upper bound", rows[0]["warnings"])
            self.assertIn("omitted", rows[0]["warnings"])
            self.assertIn("no exported structure", rows[0]["warnings"])
            self.assertEqual(scan.read_residues(Path(rows[1]["input_pdb"]))[1].name, "ALA")
            self.assertEqual(len(json.loads(rows[0]["flexible_residues"])), 3)

    def test_proline_warning_detection(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "log"
            path.write_text("ERROR: Mutation from HID to PRO not supported yet\n"
                            "Bounds are incorrect: 1 > 0\n")
            unsupported, warnings = sequential.log_details(path)
            self.assertTrue(unsupported)
            self.assertEqual(len(warnings), 1)

    def test_only_completed_score_validation_failure_is_reusable(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory / "sequence_scores.csv").write_text("sequence,upper_bound\n")
            structures = directory / "structures"
            structures.mkdir()
            (structures / "manifest.tsv").write_text("sequence\tfile\n")
            manifest = [{"phase": "exporting_structures",
                         "error": "RuntimeError: Incomplete or insufficiently precise results: warning"}]
            self.assertTrue(sequential.is_score_validation_only_failure(1, manifest, directory))
            self.assertFalse(sequential.is_score_validation_only_failure(0, manifest, directory))
            self.assertFalse(sequential.is_score_validation_only_failure(
                1, [{"phase": "scoring", "error": manifest[0]["error"]}], directory))
            (structures / "manifest.tsv").unlink()
            self.assertFalse(sequential.is_score_validation_only_failure(1, manifest, directory))

    def test_resume_retries_failure_and_rejects_partial_retry_results(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            input_pdb = base / "input.pdb"
            input_pdb.write_text(
                atom(1, "CA", "ALA", "A", 1, 0, 0, 0)
                + atom(2, "CA", "HIS", "B", 1, 1, 0, 0)
                + atom(3, "CA", "GLU", "B", 2, 2, 0, 0)
                + "END\n"
            )
            calls = []
            retry_succeeds = False

            def fake_attempt(args, current, target, directory, exclude):
                directory.mkdir()
                calls.append((target, exclude))
                succeeds = target == "B1" or (target == "B2" and exclude and retry_succeeds)
                info = {"directory": str(directory), "log": str(directory) + ".log",
                        "exit_code": 0 if succeeds else 1,
                        "exclude_proline": exclude, "proline_unsupported": False,
                        "manifest_status": "completed", "score_problems": [], "warnings": []}
                if target == "B1":
                    (directory / "sequence_scores.csv").write_text(
                        "sequence,identity,lower_bound,upper_bound\nB1=ALA,ALA,1,20\n")
                    structures = directory / "structures"
                    structures.mkdir()
                    (structures / "ALA.pdb").write_text(
                        "MODEL        1\n"
                        + atom(1, "CA", "ALA", "A", 1, 0, 0, 0)
                        + atom(2, "CA", "ALA", "B", 1, 1, 0, 0)
                        + atom(3, "CA", "GLU", "B", 2, 2, 0, 0)
                        + "ENDMDL\nEND\n")
                else:
                    (directory / "sequence_scores.csv").write_text(
                        "sequence,identity,lower_bound,upper_bound\nB2=ALA,ALA,1,20\n")
                    if succeeds:
                        structures = directory / "structures"
                        structures.mkdir()
                        (structures / "ALA.pdb").write_text(
                            "MODEL        1\n"
                            + atom(1, "CA", "ALA", "A", 1, 0, 0, 0)
                            + atom(2, "CA", "ALA", "B", 1, 1, 0, 0)
                            + atom(3, "CA", "ALA", "B", 2, 2, 0, 0)
                            + "ENDMDL\nEND\n")
                return info, []

            original = base / "original"
            with patch.object(sequential, "run_attempt", side_effect=fake_attempt):
                self.assertEqual(sequential.main(["--pdb", str(input_pdb), "--output", str(original)]), 1)
            original_manifest = original / "workflow_manifest.csv"
            before = original_manifest.read_bytes()
            calls.clear()
            continuation = base / "continuation"
            arguments = ["--pdb", str(input_pdb), "--output", str(continuation),
                         "--resume-manifest", str(original_manifest), "--resume-after", "1"]
            with patch.object(sequential, "run_attempt", side_effect=fake_attempt):
                self.assertEqual(sequential.main(arguments), 1)
            self.assertEqual(calls, [("B2", False), ("B2", True)])
            self.assertEqual(original_manifest.read_bytes(), before)
            with (continuation / "workflow_manifest.csv").open(newline="") as file:
                rows = list(csv.DictReader(file))
            self.assertEqual([row["status"] for row in rows], ["selected", "skipped"])
            self.assertEqual(rows[0]["resumed_from_manifest"], str(original_manifest))
            self.assertEqual(rows[1]["input_pdb"], rows[0]["output_pdb"])
            self.assertEqual(rows[1]["output_pdb"], rows[1]["input_pdb"])
            self.assertEqual(json.loads(rows[1]["candidate_scores"]), [])
            self.assertIn("Both attempts failed", rows[1]["reason"])
            self.assertEqual(len(json.loads(rows[1]["attempts"])), 2)

            retry_succeeds = True
            calls.clear()
            recovered = base / "recovered"
            recovered_args = ["--pdb", str(input_pdb), "--output", str(recovered),
                              "--resume-manifest", str(original_manifest), "--resume-after", "1"]
            with patch.object(sequential, "run_attempt", side_effect=fake_attempt):
                self.assertEqual(sequential.main(recovered_args), 1)
            self.assertEqual(calls, [("B2", False), ("B2", True)])
            with (recovered / "workflow_manifest.csv").open(newline="") as file:
                recovered_rows = list(csv.DictReader(file))
            self.assertEqual(recovered_rows[1]["status"], "selected")
            self.assertEqual(recovered_rows[1]["selected_identity"], "ALA")

            with self.assertRaisesRegex(ValueError, "Resume setting"):
                sequential.main(arguments + ["--radius", "5"])
            with original_manifest.open(newline="") as file:
                original_rows = list(csv.DictReader(file))
            original_rows[0]["translation_rotation_enabled"] = "False"
            with original_manifest.open("w", newline="") as file:
                writer = csv.DictWriter(file, fieldnames=sequential.FIELDS)
                writer.writeheader()
                writer.writerows(original_rows)
            with self.assertRaisesRegex(ValueError, "translation_rotation_enabled"):
                sequential.main(arguments)
            original_manifest.write_bytes(before)
            Path(rows[0]["output_pdb"]).write_text("tampered\n")
            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                sequential.main(arguments)


if __name__ == "__main__":
    unittest.main()
