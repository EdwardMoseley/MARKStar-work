#!/usr/bin/env python3
"""Score a ligand position with K* using MARK* partition functions."""

import argparse
from collections import OrderedDict
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import heapq
import itertools
import json
import math
import os
from pathlib import Path
import sys
import tempfile


DEFAULT_PDB = Path(__file__).with_name("PKN2_SAV1_forced_pilot_model_0.clean.pdb")
TRANSLATION_ROTATION_MAX_ROTATION_DEGREES = 10.0
TRANSLATION_ROTATION_MAX_TRANSLATION_ANGSTROM = 2.0
AMINO_ACIDS = (
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
)
HISTIDINE_VARIANTS = {"HID", "HIE", "HIP"}


def canonical_identity(name):
    return "HIS" if name.upper() in HISTIDINE_VARIANTS else name.upper()


@dataclass
class Residue:
    chain: str
    number: str
    name: str
    ca: tuple | None = None

    @property
    def id(self):
        return self.chain + self.number


def read_residues(path):
    residues = OrderedDict()
    current_chain = None
    closed_chains = set()
    current_residue = None
    with path.open() as pdb:
        for line_number, line in enumerate(pdb, 1):
            record = line[:6].strip()
            if record in ("MODEL", "ENDMDL"):
                raise ValueError("Multi-model PDBs are not supported; provide one model.")
            if record == "HETATM":
                raise ValueError(f"Line {line_number}: HETATM records are not supported.")
            if record != "ATOM":
                continue
            if len(line) < 54:
                raise ValueError(f"Line {line_number}: incomplete ATOM record.")
            chain = line[21]
            if chain not in ("A", "B"):
                raise ValueError(f"Line {line_number}: expected only chains A and B, got {chain!r}.")
            if chain != current_chain:
                if chain in closed_chains:
                    raise ValueError(f"Chain {chain} is not contiguous in the PDB.")
                if current_chain is not None:
                    closed_chains.add(current_chain)
                current_chain = chain
            if line[16] != " ":
                raise ValueError(f"Line {line_number}: alternate locations are not supported.")
            number = line[22:26].strip() + line[26].strip()
            if not number:
                raise ValueError(f"Line {line_number}: missing residue number.")
            name = line[17:20].strip()
            key = (chain, number)
            if key not in residues:
                residues[key] = Residue(chain, number, name)
            elif key != current_residue:
                raise ValueError(f"Line {line_number}: residue {chain}{number} is not contiguous.")
            current_residue = key
            residue = residues[key]
            if residue.name != name:
                raise ValueError(f"Line {line_number}: conflicting names for {residue.id}.")
            if line[12:16].strip() == "CA":
                if residue.ca is not None:
                    raise ValueError(f"Line {line_number}: duplicate CA for {residue.id}.")
                try:
                    xyz = tuple(float(line[i:i + 8]) for i in (30, 38, 46))
                except ValueError as exc:
                    raise ValueError(f"Line {line_number}: invalid CA coordinates.") from exc
                if not all(math.isfinite(x) for x in xyz):
                    raise ValueError(f"Line {line_number}: non-finite CA coordinates.")
                residue.ca = xyz
    for chain in ("A", "B"):
        if not any(r.chain == chain for r in residues.values()):
            raise ValueError(f"No ATOM residues found in chain {chain}.")
    for residue in residues.values():
        if residue.ca is None:
            raise ValueError(f"Residue {residue.id} has no CA; cannot rank all residues.")
    return list(residues.values())


def select_neighbors(residues, count, radius, target_id=None):
    target = next((r for r in residues if r.chain == "B"
                   and (target_id is None or r.id == target_id)), None)
    if target is None:
        raise ValueError(f"Ligand position {target_id!r} was not found in chain B.")
    candidates = sorted(
        ((math.dist(target.ca, r.ca), r) for r in residues if r is not target),
        key=lambda item: (item[0], item[1].chain, item[1].number),
    )
    return target, [(distance, r) for distance, r in candidates if distance <= radius][:count]


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", type=Path, default=DEFAULT_PDB)
    parser.add_argument("--target", help="chain B PDB residue ID, e.g. B1 or B12A (default: first chain B residue)")
    parser.add_argument("--neighbors", type=int, default=3, help="maximum number of flexible neighbors (default: 3)")
    parser.add_argument("--radius", type=float, default=4.0, help="inclusive CA distance in angstroms (default: 4)")
    parser.add_argument("--cpus", type=int, default=4)
    parser.add_argument("--heap-mib", type=int, default=4096)
    parser.add_argument("--epsilon", type=float, default=0.683)
    parser.add_argument("--top-k", type=int, default=10,
                        help="number of lowest-energy complex conformations to export per sequence (default: 10)")
    parser.add_argument("--output", type=Path, default=Path("markstar-output"))
    parser.add_argument("--exclude-proline", action="store_true",
                        help="omit PRO when the installed OSPREY cannot mutate the N terminus to proline")
    parser.add_argument("--dry-run", action="store_true", help="inspect PDB and select residues without starting OSPREY")
    parser.add_argument("--setup-only", action="store_true", help="build OSPREY conformational spaces but do not calculate scores")
    args = parser.parse_args(argv)
    if args.neighbors < 0:
        parser.error("--neighbors must be nonnegative")
    if not math.isfinite(args.radius) or args.radius < 0:
        parser.error("--radius must be finite and nonnegative")
    if args.cpus < 1 or args.heap_mib < 1:
        parser.error("--cpus and --heap-mib must be positive")
    if args.top_k < 1:
        parser.error("--top-k must be positive")
    if not math.isfinite(args.epsilon) or not 0 < args.epsilon < 1:
        parser.error("--epsilon must be between 0 and 1 (exclusive)")
    if args.dry_run and args.setup_only:
        parser.error("--dry-run and --setup-only cannot be combined")
    return args


class RunManifest:
    """Persist a run's parameters and selected residues in a spreadsheet-friendly CSV."""

    def __init__(self, output, args, target, neighbors):
        digest = hashlib.sha256()
        with args.pdb.open("rb") as pdb:
            for block in iter(lambda: pdb.read(1024 * 1024), b""):
                digest.update(block)
        self.path = output / "run_manifest.csv"
        self.data = {
            "schema_version": 1,
            "started_utc": datetime.now(timezone.utc).isoformat(),
            "updated_utc": "",
            "status": "running",
            "phase": "selected_residues",
            "error": "",
            "input_pdb": str(args.pdb.resolve()),
            "input_sha256": digest.hexdigest(),
            "output_directory": str(output),
            "python_version": sys.version.split()[0],
            "forcefield": "OSPREY default (not yet initialized)",
            "requested_neighbors": args.neighbors,
            "selected_neighbors": len(neighbors),
            "radius_angstrom": args.radius,
            "distance_metric": "CA-to-CA",
            "translation_rotation_enabled": True,
            "translation_rotation_max_rotation_degrees":
                TRANSLATION_ROTATION_MAX_ROTATION_DEGREES,
            "translation_rotation_max_translation_angstrom":
                TRANSLATION_ROTATION_MAX_TRANSLATION_ANGSTROM,
            "cpus": args.cpus,
            "heap_mib": args.heap_mib,
            "epsilon": args.epsilon,
            "top_k": args.top_k,
            "exclude_proline": args.exclude_proline,
            "requested_amino_acids": json.dumps(
                [aa for aa in AMINO_ACIDS if aa != "PRO" or not args.exclude_proline]
            ),
            "stability_filter": "disabled",
            "max_simultaneous_mutations": 1,
            "score_problems": "[]",
            "exported_identities": "[]",
        }
        self.rows = []
        for role, distance, residue in [
            ("mutation_target", 0.0, target),
            *[("neighbor", distance, residue) for distance, residue in neighbors],
        ]:
            self.rows.append({
                "role": role,
                "residue_id": residue.id,
                "chain": residue.chain,
                "residue_number": residue.number,
                "input_residue_name": residue.name,
                "osprey_wild_type": "",
                "ca_x_angstrom": residue.ca[0],
                "ca_y_angstrom": residue.ca[1],
                "ca_z_angstrom": residue.ca[2],
                "distance_to_target_ca_angstrom": distance,
                "flexibility": "continuous",
                "allowed_identities": (
                    "requested_amino_acids" if role == "mutation_target" else "wild_type_only"
                ),
                "configuration": "selected_not_configured",
            })
        self.write()

    def write(self):
        self.data["updated_utc"] = datetime.now(timezone.utc).isoformat()
        temp_path = None
        try:
            with tempfile.NamedTemporaryFile(
                mode="w", newline="", encoding="utf-8", dir=self.path.parent,
                prefix=".run_manifest.", suffix=".tmp", delete=False,
            ) as file:
                temp_path = Path(file.name)
                writer = csv.DictWriter(file, fieldnames=[*self.data, *self.rows[0]])
                writer.writeheader()
                for row in self.rows:
                    writer.writerow({**self.data, **row})
            os.replace(temp_path, self.path)
        finally:
            if temp_path is not None and temp_path.exists():
                temp_path.unlink()

    def update(self, **fields):
        self.data.update(fields)
        self.write()

    def configured(self, strands):
        for row in self.rows:
            flex = strands[row["chain"]].flexibility[row["residue_id"]]
            row["osprey_wild_type"] = str(flex.wildType)
            row["configuration"] = "configured"
        self.update(phase="spaces_built")

    def exported(self, identity):
        identities = json.loads(self.data["exported_identities"])
        identities.append(identity)
        self.update(exported_identities=json.dumps(identities))


def score_problems(scores, target_id, expected_types, epsilon):
    observed_types = [str(result.sequence.get(target_id).name).upper() for result in scores]
    problems = []
    if len(observed_types) != len(expected_types) or set(observed_types) != expected_types:
        problems.append(f"Unexpected sequence identities: expected {sorted(expected_types)}, "
                        f"got {sorted(observed_types)}.")
    for result in scores:
        if result.score.score is None:
            problems.append(f"{result.sequence}: no estimated K* score")
        for state in ("protein", "ligand", "complex"):
            pfunc = getattr(result.score, state)
            delta = pfunc.values.getEffectiveEpsilon()
            if not math.isfinite(delta) or delta > epsilon + 1e-12:
                problems.append(f"{result.sequence} {state}: effective epsilon {delta:.6g} "
                                f"exceeds requested {epsilon:g} or is non-finite")
    return problems


def write_exact_scores(scores, target_id, path):
    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(["sequence", "identity", "lower_bound", "upper_bound",
                         "score_estimated", "protein_status", "ligand_status",
                         "complex_status", "protein_epsilon", "ligand_epsilon",
                         "complex_epsilon"])
        for result in scores:
            states = [getattr(result.score, name) for name in ("protein", "ligand", "complex")]
            identity = str(result.sequence.get(target_id).name).upper()
            writer.writerow([f"{target_id}={identity}", identity, str(result.score.lowerBound),
                             str(result.score.upperBound), result.score.score is not None,
                             *[str(state.status) for state in states],
                             *[state.values.getEffectiveEpsilon() for state in states]])


def export_structures(osprey, scores, target_id, conf_ecalc, output, top_k, warnings,
                      on_export=None):
    import jpype

    structures = output / "structures"
    structures.mkdir()
    (structures / "score_validation.txt").write_text(
        "MARK* score validation failed:\n" + "\n".join(warnings) + "\n"
        if warnings else "MARK* score validation checks passed.\n"
    )
    analyzer = osprey.ConfAnalyzer(conf_ecalc)
    manifest_path = structures / "manifest.tsv"
    with manifest_path.open("w", newline="") as manifest:
        writer = csv.writer(manifest, delimiter="\t")
        writer.writerow(["sequence", "rank", "energy_kcal_per_mol", "rc_assignments",
                         "file", "total_conformations", "finite_conformations"])
        exported = set()
        for result in scores:
            identity = str(result.sequence.get(target_id).name).upper()
            sequence = str(result.sequence)
            if identity in exported or not identity.isalpha():
                raise RuntimeError(f"Cannot safely export duplicate or invalid identity {identity!r}.")
            exported.add(identity)
            rcs = result.sequence.makeRCs(conf_ecalc.confSpace)
            rc_lists = [list(rcs.get(pos)) for pos in range(rcs.getNumPos())]
            total = math.prod(len(values) for values in rc_lists)
            print(f"Enumerating {total} complex conformations for {sequence}...", flush=True)
            best = []
            finite_count = 0
            for index, assignments in enumerate(itertools.product(*rc_lists)):
                if index and index % 100 == 0:
                    print(f"  {sequence}: evaluated {index}/{total} conformations.", flush=True)
                analysis = analyzer.analyze(jpype.JArray(jpype.JInt)(assignments))
                energy = float(analysis.epmol.energy)
                if math.isnan(energy) or energy == -math.inf:
                    raise RuntimeError(f"Invalid energy {energy} for {sequence}, RCs {assignments}.")
                if energy == math.inf:
                    continue
                finite_count += 1
                entry = (-energy, -index, analysis)
                if len(best) < top_k:
                    heapq.heappush(best, entry)
                elif entry[:2] > best[0][:2]:
                    heapq.heapreplace(best, entry)
            if not best:
                raise RuntimeError(f"No finite-energy conformations found for {sequence}.")
            ranked = sorted(best, key=lambda entry: (-entry[0], -entry[1]))
            ensemble = jpype.java.util.ArrayList()
            filename = f"{identity}.pdb"
            rows = []
            for rank, (_, _, analysis) in enumerate(ranked, 1):
                ensemble.add(analysis.epmol)
                rows.append([sequence, rank, repr(float(analysis.epmol.energy)),
                             ",".join(str(rc) for rc in analysis.assignments),
                             filename, total, finite_count])
            comment = (f"{sequence}; exhaustive RC enumeration; {len(ranked)} lowest finite energies"
                       + ("; MARK* SCORE VALIDATION FAILED" if warnings else ""))
            osprey.c.structure.PDBIO.writeFile(ensemble, str(structures / filename), comment)
            writer.writerows(rows)
            manifest.flush()
            if on_export is not None:
                on_export(identity)
            if len(ranked) < top_k:
                print(f"WARNING: {sequence} has only {len(ranked)} finite conformations; "
                      f"requested {top_k}.", file=sys.stderr, flush=True)
            print(f"Wrote {len(ranked)} structures for {sequence} to {filename}.", flush=True)


def run_osprey(args, residues, target, neighbors, manifest=None):
    import osprey

    if manifest:
        manifest.update(phase="starting_osprey")
    osprey.start(heapSizeMiB=args.heap_mib)
    ffparams = osprey.ForcefieldParams()
    if manifest:
        manifest.update(forcefield=str(ffparams.forcefld))
    mol = osprey.readPdb(str(args.pdb.resolve()))
    templates = osprey.TemplateLibrary(ffparams.forcefld)
    strands = {}
    for chain in ("A", "B"):
        chain_residues = [r for r in residues if r.chain == chain]
        strands[chain] = osprey.Strand(
            mol, templateLib=templates,
            residues=[chain_residues[0].id, chain_residues[-1].id],
        )
        if not strands[chain].nonTemplateResNames.isEmpty():
            raise ValueError(f"OSPREY could not match all chain {chain} residues to templates: "
                             f"{strands[chain].nonTemplateResNames}")
    for _, residue in neighbors:
        strands[residue.chain].flexibility[residue.id].setLibraryRotamers(
            osprey.WILD_TYPE
        ).addWildTypeRotamers().setContinuous()
    # Preserve the input's wild-type template, including its protonation.
    target_types = [
        osprey.WILD_TYPE if aa == canonical_identity(target.name) else aa
        for aa in AMINO_ACIDS if aa != "PRO" or not args.exclude_proline
    ]
    target_flex = strands["B"].flexibility[target.id]
    target_flex.setLibraryRotamers(
        *target_types
    ).addWildTypeRotamers().setContinuous()
    expected_types = {
        str(target_flex.wildType).upper() if aa == canonical_identity(target.name) else aa
        for aa in AMINO_ACIDS if aa != "PRO" or not args.exclude_proline
    }

    # Match the flexibility-study setup: bound six-DOF rigid-body motion for all of chain B.
    ligand_translate_rotate = osprey.c.confspace.StrandFlex.TranslateRotate(
        TRANSLATION_ROTATION_MAX_ROTATION_DEGREES,
        TRANSLATION_ROTATION_MAX_TRANSLATION_ANGSTROM,
    )
    ligand_strand_with_pose_dofs = [strands["B"], ligand_translate_rotate]
    protein = osprey.ConfSpace(strands["A"])
    ligand = osprey.ConfSpace([ligand_strand_with_pose_dofs])
    complex_space = osprey.ConfSpace([strands["A"], ligand_strand_with_pose_dofs])
    if manifest:
        manifest.configured(strands)
    if args.setup_only:
        print("OSPREY conformational spaces built successfully.", flush=True)
        return

    output = args.output.resolve()
    original_cwd = Path.cwd()
    try:
        os.chdir(output)
        ecalc = osprey.EnergyCalculator(
            complex_space, ffparams, parallelism=osprey.Parallelism(cpuCores=args.cpus)
        )
    except BaseException:
        os.chdir(original_cwd)
        raise
    try:
        rigid_ecalc = osprey.SharedEnergyCalculator(ecalc, isMinimizing=False)
        kstar = osprey.KStar(
            protein, ligand, complex_space, epsilon=args.epsilon,
            stabilityThreshold=None, maxSimultaneousMutations=1,
            writeSequencesToConsole=True,
            writeSequencesToFile=str(output / "kstar.results.tsv"),
        )
        for info in kstar.confSpaceInfos():
            manifest.update(phase=f"energy_matrices_{info.id}")
            eref = osprey.ReferenceEnergies(info.confSpace, ecalc)
            info.confEcalc = osprey.ConfEnergyCalculator(
                info.confSpace, ecalc, referenceEnergies=eref
            )
            minimized = osprey.EnergyMatrix(info.confEcalc)
            if str(info.id) == "complex":
                complex_conf_ecalc = info.confEcalc
            rigid_conf_ecalc = osprey.ConfEnergyCalculatorCopy(info.confEcalc, rigid_ecalc)
            rigid = osprey.EnergyMatrix(rigid_conf_ecalc)

            def make_pfunc(rcs, conf_ecalc=info.confEcalc, emat=minimized,
                           rigid_conf=rigid_conf_ecalc, rigid_emat=rigid):
                return osprey.MARKStarPfunc(
                    conf_ecalc.confSpace, emat, conf_ecalc, rigid_emat, rigid_conf, rcs
                )

            info.pfuncFactory = osprey.KStar.PfuncFactory(make_pfunc)
        manifest.update(phase="scoring")
        scores = kstar.run(ecalc.tasks)
        print(f"OSPREY returned {len(scores)} scored sequences.", flush=True)
        write_exact_scores(scores, target.id, output / "sequence_scores.csv")
        problems = score_problems(scores, target.id, expected_types, args.epsilon)
        manifest.update(phase="exporting_structures", score_problems=json.dumps(problems))
        if problems:
            print("WARNING: MARK* score validation failed; exported structures are not "
                  "evidence of reliable binding scores:\n  " + "\n  ".join(problems),
                  file=sys.stderr, flush=True)
        export_structures(osprey, scores, target.id, complex_conf_ecalc,
                          output, args.top_k, problems, on_export=manifest.exported)
        if problems:
            raise RuntimeError("Incomplete or insufficiently precise results:\n  " + "\n  ".join(problems))
        print("All requested sequence identities have estimated scores within the requested "
              "partition-function epsilon.", flush=True)
    finally:
        try:
            ecalc.close()
        finally:
            os.chdir(original_cwd)


def main(argv=None):
    args = parse_args(argv)
    residues = read_residues(args.pdb)
    target, neighbors = select_neighbors(residues, args.neighbors, args.radius, args.target)
    if canonical_identity(target.name) not in AMINO_ACIDS:
        raise ValueError(f"Mutation target {target.id} has nonstandard residue name {target.name}; "
                         "expected one of the 20 standard three-letter names.")
    if args.exclude_proline and canonical_identity(target.name) == "PRO":
        raise ValueError("Cannot exclude PRO when the wild-type target is PRO.")
    print(f"Input: {args.pdb.resolve()}")
    print(f"Mutation target: {target.name} {target.id}; CA {target.ca}")
    if args.exclude_proline:
        print("PRO excluded explicitly; this is not a complete 20-identity scan.")
    print(f"Selected {len(neighbors)} of {args.neighbors} requested neighbors within {args.radius:g} Å:")
    for distance, residue in neighbors:
        print(f"  {residue.name} {residue.id}: {distance:.3f} Å")
    if len(neighbors) < args.neighbors:
        print("WARNING: fewer neighbors than requested are within the radius; using all available.",
              file=sys.stderr)
    sys.stdout.flush()
    if not args.dry_run:
        if args.setup_only:
            run_osprey(args, residues, target, neighbors)
            return
        output = args.output.resolve()
        output.mkdir(parents=True, exist_ok=False)
        manifest = RunManifest(output, args, target, neighbors)
        try:
            run_osprey(args, residues, target, neighbors, manifest)
        except (Exception, KeyboardInterrupt) as exc:
            manifest.update(status="interrupted" if isinstance(exc, KeyboardInterrupt) else "failed",
                            error=f"{type(exc).__name__}: {exc}")
            raise
        else:
            manifest.update(status="completed", phase="completed")


if __name__ == "__main__":
    main()
