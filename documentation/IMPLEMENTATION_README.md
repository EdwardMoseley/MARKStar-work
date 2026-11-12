# Implementation README: first-position MARK* scan

## Purpose and important limitation

`markstar_scan.py` sets up K* with MARK* partition-function calculations for
the first residue of peptide/ligand chain B. Chain A is the protein. The target
is flexible and can change identity; up to `n` nearby residues are flexible but
retain their wild-type identities. All other residues retain their input
internal conformations, while the entire chain B ligand also has bounded
rigid-body translation and rotation.

The single-position execution logic, including residue selection, MARK*
configuration, CSV recording, and PDB export, lives in `markstar_scan.py`.
The multi-position driver is `markstar_sequential.py`; its design and usage are
documented in [Sequential implementation](SEQUENTIAL_IMPLEMENTATION.md).
That document also describes checkpoint validation, restart provenance,
automatic retry without PRO, skip behavior, and the November 2026 continuation.
Both depend on the input PDB; the single-position run also depends on the
separately installed OSPREY and JPype packages and Java. Test modules are not
imported by the run scripts.

**The installed OSPREY build cannot complete the requested 20-identity scan for
the supplied input.** It throws `Mutation from HID to PRO not supported yet`
when calculating reference energies: its proline mutation alignment requires a
backbone `H` in the old residue template, whereas the N-terminal histidine
template lacks that atom. This occurs even though conformational-space setup
itself succeeds. To run the other 19 identities, explicitly supply
`--exclude-proline`. This is a partial scan, not a result for proline or a
complete 20-amino-acid comparison. The script does not silently substitute a
different terminus, cap the peptide, or remove proline.

The script without `--exclude-proline` still requests all 20 standard amino
acids, allowing a compatible future OSPREY implementation to attempt the full
scan. With the current environment it is expected to fail during reference
energy calculation. Fixing the N-terminal proline mutation properly requires
changes and validation in OSPREY's mutation implementation, not merely
renaming PDB atoms.

## Input and selection for the supplied structure

The default input is `PKN2_SAV1_forced_pilot_model_0.clean.pdb` next to the
script. By default the script chooses the **first chain B residue in file order**,
not necessarily residue number 1. Use `--target B2` to choose another residue.
For the supplied input the default is HIS B1, with C-alpha
coordinates `(10.926, -0.850, -19.018)` Å.

Selection is based on the input structure before any mutation or minimization:

1. Read C-alpha coordinates for every `ATOM` residue in chains A and B.
2. Measure the Euclidean distance from the target C-alpha to each other
   residue's C-alpha, including residues in either chain.
3. Keep only residues at or within the radius, inclusive.
4. Sort by increasing distance, breaking exact ties by chain and residue-ID
   string, and take at most `--neighbors`.

Defaults are `--neighbors 3 --radius 4.0`. In this PDB, GLU B2 at approximately
3.799 Å is the **only** qualifying neighbor. The next closest residues are GLU
A206 (approximately 5.668 Å) and PRO A210 (approximately 5.850 Å). The script
warns and uses the one available residue; it does not expand the radius. Thus
the default residue-level flexibility is at B1 and B2; chain A has no flexible
residues. Independently, the whole chain B strand has the bounded
translation/rotation described below. A 6 Å radius selects B2, A206, and A210 when
`--neighbors 3`. The reported three-decimal distances are rounded for display;
selection uses unrounded distances.

The PDB parser deliberately rejects multiple models, HETATM records, chains
other than A and B, alternate locations, noncontiguous chain or residue blocks,
duplicate C-alpha atoms, invalid coordinates, and residues without a C-alpha.
It requires both chains and a standard three-letter name at the target. These
restrictions avoid silently discarding data or ranking an incomplete structure.
Prepare a compatible PDB separately if necessary; this script does not prepare,
re-protonate, cap, or modify the input. OSPREY also assigns its own templates
based on the PDB and force field; if it cannot match residues to templates, the
script aborts rather than proceeding with the residues OSPREY would otherwise
automatically delete. A PDB passing the parser is not necessarily a chemically
valid or energetically appropriate input.

## Bounded ligand translation and rotation

Every scan adds
`osprey.c.confspace.StrandFlex.TranslateRotate(10, 2.0)` to the complete
chain B strand when constructing the ligand and complex conformational spaces.
This follows the construction used in `flexibility_study/IAS.py` and
`run_wt_flex_structure_markstar.py`. The protein-only conformational space
does not receive these ligand pose degrees of freedom.

OSPREY represents this as six bounded degrees of freedom: three rotations,
each bounded by ±10 degrees, and three translations, each bounded by ±2.0 Å.
These are whole-strand movements relative to the pose in the input PDB, in
addition to the selected residues' rotamer and continuous local flexibility.
They permit a bounded ligand pose search; they do not enable unrestricted
docking or search arbitrary orientations throughout space.

Translation/rotation changes the ligand and complex partition-function spaces
and therefore changes the MARK* calculation, its cost, and its scores. It is
recorded in both run manifests. Results from a run without these degrees of
freedom are not interchangeable with results from the new setup; restart from
the initial PDB at B1 instead of reusing old no-translation/rotation
checkpoints.

## Environment

The supplied environment has been checked for the `osprey` Python package:

```bash
conda activate osprey-jdk17-cckstar-vtk93
cd /home/users/etm33/src/OSPREY3/MARKStar-work
```

Alternatively invoke its interpreter directly:

```bash
/home/users/etm33/miniconda3/envs/osprey-jdk17-cckstar-vtk93/bin/python markstar_scan.py --dry-run
```

The tested environment reports OSPREY 3.3-dev, Python 3.11 and Java 17.
`--dry-run` uses only Python's standard library and does not import OSPREY or
start the JVM.

## Commands

Inspect the default residue selection:

```bash
python markstar_scan.py --dry-run
```

Validate that OSPREY can read the PDB, assign templates and construct the
protein, ligand, and complex conformational spaces:

```bash
python markstar_scan.py --setup-only --exclude-proline
```

This **does not** calculate reference energies, energy matrices, partition
functions, or K* scores, and therefore does not prove that the calculation
will finish. In particular it does not detect the proline limitation. The
translation/rotation integration was checked with `--setup-only` for B1 using
five requested neighbors, a 4 Å radius and epsilon `0.683`; OSPREY successfully
built the ligand and complex conformational spaces. This confirms that the
installed binding accepts the rigid-body strand configuration, not that the
full MARK* calculation will converge or finish.

Run the explicitly reduced scan:

```bash
python -u markstar_scan.py --exclude-proline --cpus 4 --heap-mib 4096 \
  --output results-19
```

Select up to five neighbors within 8 Å and use a different input and output:

```bash
python -u markstar_scan.py --pdb /path/to/complex.pdb --neighbors 5 \
  --radius 8 --exclude-proline --cpus 4 --heap-mib 8192 \
  --epsilon 0.1 --output results-8A
```

The output directory must **not already exist**; this avoids overwriting
results or accidentally reusing a stale energy matrix. A failed/interrupted
run may leave a partial directory. Choose a new output directory for each
attempt and retain the console output as a record. For example, in bash:

```bash
set -o pipefail
python -u markstar_scan.py --exclude-proline --output results-19 \
  2>&1 | tee results-19.log
```

Do not treat the presence of the output directory or TSV as evidence of a
completed calculation. Check the process exit status, OSPREY's per-sequence
status and bounds, and the final reported count. Some partition functions may
be aborted or otherwise not converge even if a sequence is reported. The
script verifies the identity set, checks that each K* score is estimated, and
checks the effective epsilon of each of its three partition functions; it
prints any problems, exports the available structures, and then exits with an
error if any check fails. If the MARK* calculation itself fails before
returning sequence results, there is nothing to export. The presence of PDB
files is not evidence that the binding score calculation passed.

## Run parameter and flexibility manifest

Each **new full run** automatically creates `run_manifest.csv` in its output
directory immediately after validating the PDB and selecting residues, before
starting OSPREY. Open it in a spreadsheet or read it with Python's `csv`
module. There is one row for the mutable target and one for each selected
flexible neighbor; run-level fields are repeated on each row so each row is
self-contained. For the supplied PDB and default radius, the rows describe
mutable HIS B1 and wild-type-only GLU B2. The count of requested neighbors is
three, but the selected count is one. Nonselected residues receive no local
residue flexibility; the whole chain B still has the bounded pose movement
described above.

The CSV records the input's absolute path and SHA-256, output directory,
UTC start/update times, Python version, initialized force field, requested and
selected neighbor counts, radius and CA-to-CA distance metric, CPU count, heap,
epsilon, top-k, proline-exclusion setting, requested amino-acid identities,
stability-filter and mutation-limit settings. For each flexible position it
includes the role, chain, PDB residue number including insertion code, input
name, CA coordinates and distance from the target CA, flexibility mode and
identity policy. Once OSPREY builds the spaces, `osprey_wild_type` records the
assigned template identity and `configuration` changes from
`selected_not_configured` to `configured`. For B1 this distinguishes input
HIS from OSPREY HID.

`status` is `running`, `completed`, `failed`, or `interrupted`. `phase`
indicates the most recently entered stage, such as energy matrix construction,
scoring, or structure export. `score_problems` and `exported_identities` are
JSON arrays stored in CSV cells; the latter is updated after each identity's
PDB and structure manifest entries are written. After a caught exception or
Ctrl-C, the error and last phase remain in the file and the script still
exits unsuccessfully. Updates replace the CSV atomically, so a reader sees the
previous or the new complete snapshot, not partially written rows. The CSV is
not an event history, a resume file, or proof of numerical validity: a process
killed without Python cleanup, node crash, or power loss can leave `running`
as its last recorded status. Failures before PDB selection or manifest creation
cannot be recorded there. `--dry-run` and `--setup-only` preserve their existing
behavior and do not create an output directory or CSV.

The CSV supplements, but does not capture, OSPREY's console diagnostics;
continue saving the console log as shown above. The previous run predates this
feature, so its CSV is not created retroactively. Its
[`verification.log`](../markstar-19-results/verification.log) does record HIS B1
and GLU B2 near the beginning. `structures/manifest.tsv` is a separate file
containing exported structure ranks, energies and RC assignments.

## Structure export

After MARK* returns, the script writes up to `--top-k` (default 10) complex
structures **per returned amino-acid identity**, not the top 10 identities by
K* score. Each `structures/<identity>.pdb` is a multi-model PDB: model 1 has
the lowest calculated energy for that sequence, and subsequent models are in
ascending energy order. The OSPREY PDB header lists the model energies in
the same order. For this input the wild-type histidine file is `HID.pdb`.
`structures/manifest.tsv` records sequence, rank, energy in kcal/mol, OSPREY
residue-conformation (RC) assignment, filename, total RC combinations, and
number with finite calculated energy. The RC assignment indexes follow the
complex conformational-space position order; they are not PDB residue numbers.
`structures/score_validation.txt` records the script's score-validation
outcome. This report does not capture OSPREY's console warnings, so preserve
and review the console log as well.

To select structures, the script forms the Cartesian product of all allowed
RCs for each returned sequence in the **complex** space, without pruning or
stopping on an energy bound, and asks OSPREY to minimize and evaluate each
combination with the same complex conformation-energy calculator used for
MARK*. It retains the `k` smallest finite energies, breaking equal-energy ties
in enumeration order. Positive-infinite energies are excluded; NaN or negative
infinity aborts the export. If fewer than `k` finite-energy combinations exist,
all of them are exported with a warning; zero is an error. If a sequence's
export fails, previously written files and the manifest may be partial.

This is exhaustive **over the configured discrete RC combinations**, not over
all physically possible structures. Continuous coordinates are minimized by
OSPREY for each combination; numerical minimization and the force field do not
guarantee the true global minimum within a continuous region. Reference energy
terms and the energy-partition conventions of the configured calculator are
included, so compare energies within an identity, not as binding affinities
across different mutations. Different RC assignments may converge to similar
structures. The method intentionally does not rely on the warning-prone MARK*
lower bounds or assume MARK* visited the ten lowest-energy conformations.
Nevertheless the structures are model predictions and do not resolve the
reported K* score or energy-bound warnings.

Exhaustive enumeration scales as the **product** of RC counts at all flexible
positions. Runtime can become very large as `--neighbors` or the amino-acid
rotamer count increases; this additional work is independent of `--top-k`.
The export is performed sequentially and can take substantially longer than
the MARK* score calculation. `--top-k` controls retained structures, not the
number of conformations evaluated. Retaining `k` molecules also increases
memory consumption. `--dry-run` and `--setup-only` do not export structures.

The script runs the OSPREY calculation with the new output directory as its
working directory. OSPREY's `protein.confdb`, `ligand.confdb`, and
`complex.confdb` are therefore isolated there, rather than overwriting
databases in the caller's working directory. They are OSPREY calculation
artifacts, not a supported resume interface for this script.

### Observed result for this PDB

The 19-identity calculation with two CPUs and the default 4 Å radius wrote
[`markstar-19-results/kstar.results.tsv`](../markstar-19-results/kstar.results.tsv).
It contains 19 sequence records, including the wild-type HID and excluding
PRO. **It does not satisfy the requested epsilon for every state**: ASP B1's
complex effective epsilon is approximately 0.99664, above 0.68. Its K* bounds
are correspondingly broad. In this build OSPREY returned an estimated score
despite that discrepancy. The current script explicitly flags this situation
and exits nonzero; the saved TSV remains available for inspection, not as a
fully converged 19-identity scan. This check is based on OSPREY's reported
effective epsilon and does not independently establish the validity of its
energy bounds.

An independent repeat with the current script is recorded in
[`markstar-19-results/verification.log`](../markstar-19-results/verification.log).
In addition to the ASP precision failure, OSPREY printed `Bounds are
incorrect` / `The bounds are probably wrong` while processing ARG B1. This is
an additional reason **not to rely on these scores or their claimed bounds**
without investigating OSPREY's energy-bound calculation. The post-run
identity/epsilon checks do not detect such internal warnings. Review the
complete console log for each new run.

A separate structure-export validation was stopped after completing wild-type
HID and ALA exports because exhaustive enumeration is expensive. See
[`markstar-19-with-structures/INCOMPLETE.txt`](../markstar-19-with-structures/INCOMPLETE.txt);
that directory is **not** a completed scan and its manifest is incomplete.

### Options

| Option | Default | Meaning |
| --- | --- | --- |
| `--pdb` | supplied PKN2/SAV1 PDB | Structure containing protein A and ligand B. |
| `--target` | first chain B residue | Chain B PDB residue ID, including any insertion code. |
| `--neighbors` | 3 | Maximum number of additional flexible residues; zero is allowed. |
| `--radius` | 4.0 | Inclusive C-alpha radius in Å; must be finite and nonnegative. |
| `--exclude-proline` | off | Explicitly omit PRO; cannot be used if the target is itself PRO. |
| `--cpus` | 4 | CPU cores for OSPREY's energy calculator. |
| `--heap-mib` | 4096 | JVM heap limit in MiB; OSPREY reserves part for garbage collection. |
| `--epsilon` | 0.683 | Requested partition-function approximation tolerance, strictly between 0 and 1; smaller is generally more expensive and accurate. |
| `--top-k` | 10 | Maximum number of lowest calculated-energy complex structures exported per returned sequence; must be positive. |
| `--output` | `markstar-output` | New output directory, relative to the current working directory unless absolute. |
| `--dry-run` | off | Parse and print selection only, without starting OSPREY. |
| `--setup-only` | off | Construct OSPREY spaces without calculating energies or scores. |

Use `python markstar_scan.py --help` for the CLI summary. `--dry-run` and
`--setup-only` are mutually exclusive.

## OSPREY calculation details

The implementation follows `examples/python.KStar/markstar.kstar.py` in the
adjacent OSPREY source tree:

* OSPREY's default `ForcefieldParams` and a shared template library are used.
  The input residue range for each chain is inclusive, from its first to last
  residue in file order.
* Each selected neighbor uses OSPREY's `WILD_TYPE` library rotamers, includes
  its input wild-type rotamer, and enables continuous flexibility. The target
  uses library rotamers for all requested standard identities, includes the
  input wild-type rotamer, and enables continuous flexibility. The actual
  wild-type template is used for the target identity, preserving the
  template/protonation OSPREY assigned from the input; for the supplied HIS B1
  OSPREY assigns HID. Other histidine protonation states are not independently
  enumerated. Flexible does not mean unrestricted backbone movement.
* Separate protein, ligand and complex conformational spaces are constructed.
  The ligand and complex spaces include the full chain B strand plus the
  `TranslateRotate(10, 2.0)` rigid-body flexibility; the protein-only space
  does not. The six pose degrees of freedom have independent ±10-degree
  rotation and ±2.0-Å translation bounds.
  A minimizing energy calculator is shared across the spaces; MARK* also gets
  a non-minimizing shared calculator and a rigid energy matrix.
* For each space, the script calculates reference energies, a minimizing
  conformation-energy calculator and matrix, a rigid calculator copy and rigid
  matrix. Each K* partition-function factory constructs
  `osprey.MARKStarPfunc`, matching the OSPREY example. Matrices are not cached;
  reruns recompute them, avoiding stale-cache reuse for a changed PDB or
  configuration.
* K* is explicitly configured with at most one simultaneous mutation. Only B1
  is mutable, so the intended sequences are wild type plus one sequence for
  each other requested identity. The stability filter is disabled
  (`stabilityThreshold=None`) to avoid intentionally filtering identities. It
  cannot guarantee that every partition-function calculation is useful or
  converges. After execution the script requires exactly the expected
  identities, an estimated K* score for each sequence, and finite effective
  epsilon at or below the requested value for each state. These checks do not
  independently validate the model or the bounds.
* The OSPREY output `kstar.results.tsv` contains sequences, score bounds and
  partition-function details in OSPREY's own format. K* is a binding score
  derived from partition functions of complex and unbound states; it is not a
  direct experimental affinity. Compare scores only after checking numerical
  bounds, convergence, modeling assumptions and the protonation/force-field
  choices.
* `sequence_scores.csv` additionally records unrounded K* bounds, each state's
  status and effective epsilon after K* returns. The sequential driver uses
  this file, not the rounded TSV, to rank identities.

This is not a BBK* best-sequence search or an exhaustive scan of every ligand
position. Energy matrices, partition functions, and exhaustive structure
enumeration can be computationally expensive. The added rigid-body degrees of
freedom further increase this cost; no completion-time guarantee is made.

## Tests

From this directory:

```bash
python -m unittest -v test_markstar_scan
```

The standard-library tests exercise the real PDB's 4 Å and 6 Å selections,
inclusive cutoff, tie ordering, zero-neighbor selection, argument validation,
rejection of missing or alternate C-alpha atoms, and rejection of missing
identities or an above-tolerance effective epsilon. They also check manifest
residue records, translation/rotation settings, the epsilon default, updates
and failure/completion states. The separate `--setup-only` smoke test confirms
the installed OSPREY binding constructs the new spaces. These checks do not
replace a MARK* scoring run or validate physical predictions.
