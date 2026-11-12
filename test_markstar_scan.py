import contextlib
import csv
import io
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import markstar_scan as scan


def atom(serial, name, residue, chain, number, x, y, z, alt=" "):
    return (f"ATOM  {serial:5d} {name:4s}{alt}{residue:3s} {chain}{number:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}           C\n")


class SelectionTests(unittest.TestCase):
    def test_epsilon_default(self):
        self.assertEqual(scan.parse_args([]).epsilon, 0.683)

    def test_real_input_default_and_three_neighbors(self):
        residues = scan.read_residues(scan.DEFAULT_PDB)
        target, selected = scan.select_neighbors(residues, 3, 4)
        self.assertEqual((target.id, target.name), ("B1", "HIS"))
        self.assertEqual([res.id for _, res in selected], ["B2"])
        self.assertAlmostEqual(selected[0][0], 3.799, places=3)
        _, selected = scan.select_neighbors(residues, 3, 6)
        self.assertEqual([res.id for _, res in selected], ["B2", "A206", "A210"])

    def test_inclusive_radius_ties_and_zero(self):
        residues = [
            scan.Residue("B", "1", "HIS", (0, 0, 0)),
            scan.Residue("B", "2", "ALA", (4, 0, 0)),
            scan.Residue("A", "1", "ALA", (0, 4, 0)),
            scan.Residue("A", "2", "ALA", (4.001, 0, 0)),
        ]
        _, selected = scan.select_neighbors(residues, 3, 4)
        self.assertEqual([res.id for _, res in selected], ["A1", "B2"])
        self.assertEqual(scan.select_neighbors(residues, 0, 4)[1], [])

    def test_invalid_inputs(self):
        for options in (["--neighbors", "-1"], ["--radius", "nan"],
                        ["--epsilon", "1"], ["--cpus", "0"], ["--top-k", "0"]):
            with self.subTest(options=options), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit):
                    scan.parse_args(options)

    def test_reject_missing_ca_and_alternates(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.pdb"
            path.write_text(atom(1, "N", "ALA", "A", 1, 0, 0, 0)
                            + atom(2, "CA", "HIS", "B", 1, 0, 0, 0))
            with self.assertRaisesRegex(ValueError, "A1 has no CA"):
                scan.read_residues(path)
            path.write_text(atom(1, "CA", "ALA", "A", 1, 0, 0, 0, alt="A")
                            + atom(2, "CA", "HIS", "B", 1, 0, 0, 0))
            with self.assertRaisesRegex(ValueError, "alternate locations"):
                scan.read_residues(path)

    def test_validate_scores_checks_identity_and_effective_epsilon(self):
        class Sequence:
            def __init__(self, identity):
                self.identity = identity

            def get(self, target):
                self.target = target
                return SimpleNamespace(name=self.identity)

            def __str__(self):
                return "B1=" + self.identity

        def result(identity, complex_delta):
            def pfunc(delta):
                return SimpleNamespace(values=SimpleNamespace(getEffectiveEpsilon=lambda: delta))
            return SimpleNamespace(
                sequence=Sequence(identity),
                score=SimpleNamespace(score=1, protein=pfunc(0), ligand=pfunc(0.5),
                                      complex=pfunc(complex_delta)),
            )

        scores = [result("HID", 0.6), result("ASP", 0.68)]
        self.assertEqual(scan.score_problems(scores, "B1", {"HID", "ASP"}, 0.68), [])
        self.assertIn("Unexpected sequence identities",
                      scan.score_problems(scores, "B1", {"HID", "ALA"}, 0.68)[0])
        scores[1] = result("ASP", 0.9966)
        self.assertIn("B1=ASP complex: effective epsilon",
                      scan.score_problems(scores, "B1", {"HID", "ASP"}, 0.68)[0])

    def test_export_ranks_all_assignments_and_marks_warnings(self):
        energies = {(0, 0): 2.0, (0, 1): 1.0, (1, 0): 1.0, (1, 1): float("inf")}

        class RCs:
            def getNumPos(self):
                return 2

            def get(self, pos):
                return [0, 1]

        class Sequence:
            def makeRCs(self, conf_space):
                return RCs()

            def get(self, target):
                return SimpleNamespace(name="ALA")

            def __str__(self):
                return "B1=ALA"

        class Analyzer:
            def analyze(self, assignments):
                return SimpleNamespace(assignments=assignments,
                                       epmol=SimpleNamespace(energy=energies[tuple(assignments)]))

        class ArrayList(list):
            def add(self, value):
                self.append(value)

        def write_pdb(ensemble, path, comment):
            Path(path).write_text(comment)

        fake_jpype = SimpleNamespace(JInt=int, JArray=lambda _: list,
                                     java=SimpleNamespace(util=SimpleNamespace(ArrayList=ArrayList)))
        fake_osprey = SimpleNamespace(
            ConfAnalyzer=lambda _: Analyzer(),
            c=SimpleNamespace(structure=SimpleNamespace(PDBIO=SimpleNamespace(writeFile=write_pdb))),
        )
        result = SimpleNamespace(sequence=Sequence())
        with tempfile.TemporaryDirectory() as tmp, patch.dict("sys.modules", {"jpype": fake_jpype}):
            scan.export_structures(fake_osprey, [result], "B1",
                                   SimpleNamespace(confSpace=None), Path(tmp), 2, ["bad"])
            pdb = Path(tmp) / "structures" / "ALA.pdb"
            self.assertIn("MARK* SCORE VALIDATION FAILED", pdb.read_text())
            with (pdb.parent / "manifest.tsv").open() as manifest:
                rows = list(csv.DictReader(manifest, delimiter="\t"))
            self.assertEqual([r["rc_assignments"] for r in rows], ["0,1", "1,0"])
            self.assertEqual([r["energy_kcal_per_mol"] for r in rows], ["1.0", "1.0"])
            self.assertEqual([r["finite_conformations"] for r in rows], ["3", "3"])

    def test_run_manifest_records_residues_progress_and_failure(self):
        def fake_run(args, residues, target, neighbors, manifest):
            manifest.update(phase="scoring")
            manifest.exported("HID")
            raise RuntimeError("test failure")

        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "run"
            with patch.object(scan, "run_osprey", side_effect=fake_run):
                with self.assertRaisesRegex(RuntimeError, "test failure"):
                    scan.main(["--output", str(output), "--exclude-proline",
                               "--neighbors", "3", "--radius", "4", "--top-k", "5"])
            with (output / "run_manifest.csv").open(newline="") as file:
                rows = list(csv.DictReader(file))
            self.assertEqual([(row["role"], row["residue_id"]) for row in rows],
                             [("mutation_target", "B1"), ("neighbor", "B2")])
            self.assertEqual(rows[0]["status"], "failed")
            self.assertEqual(rows[0]["phase"], "scoring")
            self.assertEqual(rows[0]["error"], "RuntimeError: test failure")
            self.assertEqual(rows[0]["requested_neighbors"], "3")
            self.assertEqual(rows[0]["selected_neighbors"], "1")
            self.assertEqual(rows[0]["translation_rotation_enabled"], "True")
            self.assertEqual(rows[0]["translation_rotation_max_rotation_degrees"], "10.0")
            self.assertEqual(rows[0]["translation_rotation_max_translation_angstrom"], "2.0")
            self.assertEqual(rows[0]["top_k"], "5")
            self.assertEqual(json.loads(rows[0]["exported_identities"]), ["HID"])
            self.assertNotIn("PRO", json.loads(rows[0]["requested_amino_acids"]))
            self.assertAlmostEqual(float(rows[1]["distance_to_target_ca_angstrom"]), 3.799, places=3)
            self.assertEqual(rows[0]["configuration"], "selected_not_configured")
            self.assertEqual(len(rows[0]["input_sha256"]), 64)

    def test_manifest_configured_and_completed(self):
        def fake_run(args, residues, target, neighbors, manifest):
            class Flexibility:
                def __getitem__(self, key):
                    return SimpleNamespace(wildType={"B1": "hid", "B2": "GLU"}[key])
            manifest.configured({"B": SimpleNamespace(flexibility=Flexibility())})

        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "run"
            with patch.object(scan, "run_osprey", side_effect=fake_run):
                scan.main(["--output", str(output)])
            with (output / "run_manifest.csv").open(newline="") as file:
                rows = list(csv.DictReader(file))
            self.assertEqual([row["osprey_wild_type"] for row in rows], ["hid", "GLU"])
            self.assertTrue(all(row["configuration"] == "configured" for row in rows))
            self.assertTrue(all(row["status"] == "completed" for row in rows))

    def test_manifest_records_keyboard_interrupt(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "run"
            with patch.object(scan, "run_osprey", side_effect=KeyboardInterrupt):
                with self.assertRaises(KeyboardInterrupt):
                    scan.main(["--output", str(output)])
            with (output / "run_manifest.csv").open(newline="") as file:
                row = next(csv.DictReader(file))
            self.assertEqual(row["status"], "interrupted")
            self.assertEqual(row["phase"], "selected_residues")


if __name__ == "__main__":
    unittest.main()
