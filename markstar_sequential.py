#!/usr/bin/env python3
"""Run a sequential ligand scan, carrying one chosen structure to each next step."""

import argparse
import csv
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

import markstar_scan as scan


SCAN_SCRIPT = Path(scan.__file__).resolve()
PROLINE_ERROR = re.compile(r"Mutation from \S+ to PRO not supported yet", re.IGNORECASE)
WARNING_PATTERN = re.compile(
    r"warning|bounds are incorrect|bounds are probably wrong|"
    r"error\. bounds got looser|overcorrected", re.IGNORECASE
)
FIELDS = [
    "index", "target", "status", "started_utc", "updated_utc", "input_pdb",
    "input_sha256", "output_pdb", "output_sha256", "radius_angstrom",
    "max_flex_neighbors", "translation_rotation_enabled",
    "translation_rotation_max_rotation_degrees",
    "translation_rotation_max_translation_angstrom",
    "cpus", "heap_mib", "epsilon", "top_k",
    "initial_pdb", "initial_sha256", "scan_script_sha256", "driver_script_sha256",
    "flexible_residues", "attempts", "candidate_scores", "selected_identity",
    "selected_upper_bound", "selected_pdb", "warnings", "reason", "resumed_from_manifest",
]


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as file:
        for block in iter(lambda: file.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_manifest(path, rows):
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", newline="", encoding="utf-8", dir=path.parent,
            prefix=".workflow.", suffix=".tmp", delete=False,
        ) as file:
            temporary = Path(file.name)
            writer = csv.DictWriter(file, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def log_details(path):
    warnings = []
    proline_unsupported = False
    with path.open(encoding="utf-8", errors="replace") as log:
        for line in log:
            if PROLINE_ERROR.search(line):
                proline_unsupported = True
            if WARNING_PATTERN.search(line):
                text = line.strip()
                if text not in warnings:
                    warnings.append(text)
    return proline_unsupported, warnings


def read_attempt_manifest(path):
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as file:
        return list(csv.DictReader(file))


def is_score_validation_only_failure(exit_code, manifest, directory):
    if exit_code == 0 or not manifest:
        return False
    error = manifest[0].get("error", "")
    return (
        manifest[0].get("phase") == "exporting_structures"
        and error.startswith("RuntimeError: Incomplete or insufficiently precise results:")
        and (directory / "sequence_scores.csv").is_file()
        and (directory / "structures" / "manifest.tsv").is_file()
    )


def read_candidates(path, target):
    if not path.exists():
        return [], []
    candidates = []
    warnings = []
    with path.open(newline="", encoding="utf-8") as file:
        reader = csv.DictReader(file)
        if not {"sequence", "identity", "lower_bound", "upper_bound"}.issubset(reader.fieldnames or []):
            raise ValueError(f"Missing required score columns in {path}")
        seen = set()
        for row in reader:
            sequence = row["sequence"].strip()
            prefix = target + "="
            identity = sequence[len(prefix):].upper() if sequence.startswith(prefix) else ""
            if not identity.isalpha() or identity in seen or identity != row["identity"].upper():
                raise ValueError(f"Invalid or duplicate sequence in {path}: {sequence!r}")
            seen.add(identity)
            text = row["upper_bound"].strip()
            try:
                upper = Decimal(text)
            except InvalidOperation:
                warnings.append(f"{sequence}: invalid upper bound {text!r}")
                continue
            if not upper.is_finite():
                warnings.append(f"{sequence}: non-finite upper bound {text!r}; omitted")
                continue
            if upper <= 0:
                warnings.append(f"{sequence}: non-positive upper bound {text!r}; omitted")
                continue
            candidates.append({"sequence": sequence, "identity": identity,
                               "lower_bound": row["lower_bound"], "upper_bound": text,
                               "score_estimated": row.get("score_estimated", ""),
                               "protein_status": row.get("protein_status", ""),
                               "ligand_status": row.get("ligand_status", ""),
                               "complex_status": row.get("complex_status", ""),
                               "protein_epsilon": row.get("protein_epsilon", ""),
                               "ligand_epsilon": row.get("ligand_epsilon", ""),
                               "complex_epsilon": row.get("complex_epsilon", "")})
    return candidates, warnings


def ordered_candidates(candidates):
    by_identity = sorted(candidates, key=lambda row: row["identity"])
    return sorted(by_identity, key=lambda row: Decimal(row["upper_bound"]), reverse=True)


def choose_candidate(candidates):
    if not candidates:
        return None
    return ordered_candidates(candidates)[0]


def extract_first_model(source, destination, input_residues, target, identity):
    model = []
    in_model = False
    complete = False
    with source.open(encoding="utf-8") as file:
        for line in file:
            if line.startswith("MODEL "):
                if in_model:
                    raise ValueError(f"Nested MODEL records in {source}")
                in_model = True
            elif line.startswith("ENDMDL") and in_model:
                complete = True
                break
            elif in_model and line.startswith(("ATOM  ", "TER   ")):
                model.append(line)
    if not complete or not model:
        raise ValueError(f"No complete first PDB model in {source}")
    temporary = destination.with_suffix(".tmp")
    try:
        with temporary.open("x", encoding="utf-8") as file:
            file.writelines(model)
            file.write("END\n")
        output_residues = scan.read_residues(temporary)
        before = [(r.id, scan.canonical_identity(r.name)) for r in input_residues]
        after = [(r.id, scan.canonical_identity(r.name)) for r in output_residues]
        expected = [(res_id, scan.canonical_identity(identity) if res_id == target else name)
                    for res_id, name in before]
        if after != expected:
            raise ValueError(f"First model of {source} changes residue IDs or unexpected identities.")
        os.replace(temporary, destination)
    finally:
        if temporary.exists():
            temporary.unlink()


def run_attempt(args, input_pdb, target, directory, exclude_proline):
    log_path = directory.parent / (directory.name + ".log")
    command = [
        sys.executable, "-u", str(SCAN_SCRIPT), "--pdb", str(input_pdb),
        "--target", target, "--neighbors", str(args.max_flex), "--radius", str(args.radius),
        "--cpus", str(args.cpus), "--heap-mib", str(args.heap_mib), "--epsilon",
        str(args.epsilon), "--top-k", str(args.top_k), "--output", str(directory),
    ]
    if exclude_proline:
        command.append("--exclude-proline")
    with log_path.open("x", encoding="utf-8") as log:
        completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
    unsupported, warnings = log_details(log_path)
    manifest = read_attempt_manifest(directory / "run_manifest.csv")
    validation_only_failure = is_score_validation_only_failure(
        completed.returncode, manifest, directory
    )
    return {
        "directory": str(directory), "log": str(log_path),
        "exit_code": completed.returncode, "exclude_proline": exclude_proline,
        "proline_unsupported": unsupported,
        "score_validation_only_failure": bool(validation_only_failure),
        "manifest_status": manifest[0]["status"] if manifest else "unavailable",
        "score_problems": json.loads(manifest[0]["score_problems"]) if manifest else [],
        "warnings": warnings,
    }, manifest


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", type=Path, default=scan.DEFAULT_PDB)
    parser.add_argument("--output", type=Path, default=Path("sequential-output"))
    parser.add_argument("--max-flex", type=int, default=4, help="maximum additional flexible residues")
    parser.add_argument("--radius", type=float, default=4.0)
    parser.add_argument("--cpus", type=int, default=4)
    parser.add_argument("--heap-mib", type=int, default=4096)
    parser.add_argument("--epsilon", type=float, default=0.683)
    parser.add_argument("--top-k", type=int, default=10)
    parser.add_argument("--resume-manifest", type=Path,
                        help="prior workflow manifest; write the continuation to a new output directory")
    parser.add_argument("--resume-after", type=int,
                        help="last completed step to reuse from the prior manifest")
    args = parser.parse_args(argv)
    if args.max_flex < 0 or not math.isfinite(args.radius) or args.radius < 0:
        parser.error("--max-flex and --radius must be nonnegative and radius must be finite")
    if args.cpus < 1 or args.heap_mib < 1 or args.top_k < 1:
        parser.error("--cpus, --heap-mib and --top-k must be positive")
    if not math.isfinite(args.epsilon) or not 0 < args.epsilon < 1:
        parser.error("--epsilon must be between 0 and 1 (exclusive)")
    if (args.resume_manifest is None) != (args.resume_after is None):
        parser.error("--resume-manifest and --resume-after must be supplied together")
    if args.resume_after is not None and args.resume_after < 1:
        parser.error("--resume-after must be positive")
    return args


def load_resume(args, targets, initial, initial_hash, scan_hash):
    source = args.resume_manifest.resolve()
    with source.open(newline="", encoding="utf-8") as file:
        prior = list(csv.DictReader(file))
    if len(prior) != len(targets) or args.resume_after > len(targets):
        raise ValueError("Resume manifest does not match the input ligand positions.")
    settings = {
        "radius_angstrom": args.radius, "max_flex_neighbors": args.max_flex,
        "translation_rotation_max_rotation_degrees":
            scan.TRANSLATION_ROTATION_MAX_ROTATION_DEGREES,
        "translation_rotation_max_translation_angstrom":
            scan.TRANSLATION_ROTATION_MAX_TRANSLATION_ANGSTROM,
        "cpus": args.cpus, "heap_mib": args.heap_mib,
        "epsilon": args.epsilon, "top_k": args.top_k,
    }
    previous = initial
    expected_hash = initial_hash
    for index, row in enumerate(prior, 1):
        if row.get("index") != str(index) or row.get("target") != targets[index - 1]:
            raise ValueError(f"Resume manifest target mismatch at step {index}.")
        if index > args.resume_after:
            continue
        if row.get("status") not in ("selected", "skipped"):
            raise ValueError(f"Resume step {index} is not completed.")
        if row.get("initial_sha256") != initial_hash or row.get("scan_script_sha256") != scan_hash:
            raise ValueError(f"Input PDB or scan script differs from resume step {index}.")
        if row.get("translation_rotation_enabled") != "True":
            raise ValueError(f"Resume setting translation_rotation_enabled differs at step {index}.")
        for key, value in settings.items():
            try:
                matches = Decimal(row[key]) == Decimal(str(value))
            except (KeyError, InvalidOperation):
                matches = False
            if not matches:
                raise ValueError(f"Resume setting {key} differs at step {index}.")
        if (row.get("input_pdb") != str(previous)
                or row.get("input_sha256") != expected_hash
                or sha256(previous) != expected_hash):
            raise ValueError(f"Resume input chain or checksum mismatch at step {index}.")
        output = Path(row["output_pdb"])
        expected_hash = row["output_sha256"]
        if not output.is_absolute() or not expected_hash or sha256(output) != expected_hash:
            raise ValueError(f"Resume output checksum mismatch at step {index}.")
        if row["status"] == "skipped" and output != previous:
            raise ValueError(f"Skipped resume step {index} did not carry its input forward.")
        previous = output
    return prior[:args.resume_after], previous, source


def main(argv=None):
    args = parse_args(argv)
    initial = args.pdb.resolve()
    residues = scan.read_residues(initial)
    targets = [residue.id for residue in residues if residue.chain == "B"]
    expected_ids = [r.id for r in residues]
    output = args.output.resolve()
    manifest_path = output / "workflow_manifest.csv"
    initial_hash = sha256(initial)
    scan_hash = sha256(SCAN_SCRIPT)
    driver_hash = sha256(Path(__file__).resolve())
    reused, current, source = [], initial, None
    if args.resume_manifest is not None:
        reused, current, source = load_resume(args, targets, initial, initial_hash, scan_hash)
    output.mkdir(parents=True, exist_ok=False)
    rows = []
    for index, target in enumerate(targets, 1):
        rows.append({
            "index": index, "target": target, "status": "pending", "started_utc": "",
            "updated_utc": timestamp(), "input_pdb": "", "input_sha256": "",
            "output_pdb": "", "output_sha256": "", "radius_angstrom": args.radius,
            "max_flex_neighbors": args.max_flex,
            "translation_rotation_enabled": True,
            "translation_rotation_max_rotation_degrees":
                scan.TRANSLATION_ROTATION_MAX_ROTATION_DEGREES,
            "translation_rotation_max_translation_angstrom":
                scan.TRANSLATION_ROTATION_MAX_TRANSLATION_ANGSTROM,
            "cpus": args.cpus,
            "heap_mib": args.heap_mib, "epsilon": args.epsilon, "top_k": args.top_k,
            "initial_pdb": str(initial), "initial_sha256": initial_hash,
            "scan_script_sha256": scan_hash, "driver_script_sha256": driver_hash,
            "flexible_residues": "[]", "attempts": "[]", "candidate_scores": "[]",
            "selected_identity": "", "selected_upper_bound": "", "selected_pdb": "",
            "warnings": "[]", "reason": "",
            "resumed_from_manifest": "",
        })
    for index, old in enumerate(reused):
        rows[index] = {field: old.get(field, "") for field in FIELDS}
        rows[index]["resumed_from_manifest"] = str(source)
    write_manifest(manifest_path, rows)
    had_warnings = any(row["status"] == "skipped" or json.loads(row["warnings"])
                       for row in reused)
    for row in rows[len(reused):]:
        step = output / f"step_{row['index']:03d}"
        step.mkdir()
        current_residues = scan.read_residues(current)
        if [r.id for r in current_residues] != expected_ids:
            raise ValueError(f"Residue IDs changed in {current}; cannot continue.")
        target, neighbors = scan.select_neighbors(current_residues, args.max_flex,
                                                   args.radius, row["target"])
        row.update(status="running", started_utc=timestamp(), updated_utc=timestamp(),
                   input_pdb=str(current), input_sha256=sha256(current),
                   flexible_residues=json.dumps(
                       [{"id": target.id, "name": target.name, "role": "target", "distance": 0.0}]
                       + [{"id": r.id, "name": r.name, "role": "neighbor", "distance": d}
                          for d, r in neighbors]))
        write_manifest(manifest_path, rows)
        print(f"Step {row['index']}/{len(rows)}: {target.id}; input {current}", flush=True)
        attempts = []
        warnings = []
        try:
            attempts.append({"directory": str(step / "attempt_01"), "status": "running",
                             "exclude_proline": False})
            row.update(attempts=json.dumps(attempts), updated_utc=timestamp())
            write_manifest(manifest_path, rows)
            attempt, attempt_manifest = run_attempt(args, current, target.id, step / "attempt_01", False)
            attempts[-1] = attempt
            row.update(attempts=json.dumps(attempts), updated_utc=timestamp())
            write_manifest(manifest_path, rows)
            if attempt["exit_code"] != 0 or attempt["proline_unsupported"]:
                warnings.append("All-identity attempt failed or reported unsupported PRO; retried without PRO.")
                attempts.append({"directory": str(step / "attempt_02"), "status": "running",
                                 "exclude_proline": True})
                row.update(attempts=json.dumps(attempts), warnings=json.dumps(warnings),
                           updated_utc=timestamp())
                write_manifest(manifest_path, rows)
                attempt, attempt_manifest = run_attempt(args, current, target.id, step / "attempt_02", True)
                attempts[-1] = attempt
                row.update(attempts=json.dumps(attempts), updated_utc=timestamp())
                write_manifest(manifest_path, rows)
            for attempt_info in attempts:
                warnings.extend(attempt_info["warnings"])
                warnings.extend(attempt_info["score_problems"])
                if attempt_info["exit_code"] != 0:
                    if attempt_info.get("score_validation_only_failure"):
                        warnings.append(
                            f"Attempt {attempt_info['directory']} completed scoring and structure "
                            "export but exited after score-validation warnings; finite bounds may be used."
                        )
                    else:
                        warnings.append(
                            f"Attempt {attempt_info['directory']} exited {attempt_info['exit_code']}."
                        )
            score_file = Path(attempt["directory"]) / "sequence_scores.csv"
            usable_attempt = (
                attempt["exit_code"] == 0 or attempt.get("score_validation_only_failure")
            ) and not attempt["proline_unsupported"]
            if not usable_attempt:
                candidates, score_warnings = [], []
                warnings.append("Retry without PRO failed; partial results were not used.")
            else:
                candidates, score_warnings = read_candidates(score_file, target.id)
            warnings.extend(score_warnings)
            row["candidate_scores"] = json.dumps(candidates)
            chosen = None
            for candidate in ordered_candidates(candidates):
                source = Path(attempt["directory"]) / "structures" / f"{candidate['identity']}.pdb"
                if not source.is_file():
                    warnings.append(
                        f"{candidate['sequence']}: finite upper bound but no exported structure; "
                        "trying the next candidate."
                    )
                    continue
                try:
                    destination = step / "next_input.pdb"
                    extract_first_model(source, destination, current_residues,
                                        target.id, candidate["identity"])
                except ValueError as exc:
                    warnings.append(
                        f"{candidate['sequence']}: exported structure is invalid ({exc}); "
                        "trying the next candidate."
                    )
                    continue
                chosen = candidate
                current = destination
                row.update(status="selected", selected_identity=candidate["identity"],
                           selected_upper_bound=candidate["upper_bound"],
                           selected_pdb=str(source))
                break
            if chosen is None:
                if not usable_attempt:
                    reason = "Both attempts failed; input carried forward."
                elif not candidates:
                    reason = "No finite positive K* upper bounds; input carried forward."
                else:
                    reason = "No finite-bound candidate had a valid exported structure; input carried forward."
                row.update(status="skipped", reason=reason)
                print(f"ERROR: {target.id}: {reason}", file=sys.stderr, flush=True)
            if len(neighbors) < args.max_flex:
                warnings.append(f"Selected {len(neighbors)} of {args.max_flex} requested neighbors within radius.")
            row["warnings"] = json.dumps(warnings)
            row["output_pdb"] = str(current)
            row["output_sha256"] = sha256(current)
            if warnings or row["status"] == "skipped":
                had_warnings = True
            row["updated_utc"] = timestamp()
            write_manifest(manifest_path, rows)
        except (Exception, KeyboardInterrupt) as exc:
            row.update(status="interrupted" if isinstance(exc, KeyboardInterrupt) else "failed",
                       reason=f"{type(exc).__name__}: {exc}", attempts=json.dumps(attempts),
                       warnings=json.dumps(warnings), updated_utc=timestamp())
            write_manifest(manifest_path, rows)
            raise
    print(f"Workflow finished; manifest: {manifest_path}", flush=True)
    if had_warnings:
        print("WARNING: one or more steps were skipped or had warnings; inspect the manifest and logs.",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
