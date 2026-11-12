# Sequential ligand-position MARK* workflow

## What this workflow does

`markstar_sequential.py` runs a series of single-position jobs using
`markstar_scan.py`. It starts with a protein in chain A and a peptide ligand in
chain B, and visits **every chain B residue in the input PDB's file order**,
using the actual PDB residue identifiers, including insertion codes. It does
not assume residue numbering starts at one or is consecutive. For the supplied
input, there are 15 ligand positions. These positions are captured before any
mutations, and each later structure must contain the same residue IDs in the
same order.

The workflow is a **greedy, path-dependent design**, not 15 independent scans:
the selected structure and mutations from one step become the input to the
next. A later flexible neighbor may change its conformation again but retains
its current amino-acid identity. Earlier mutations can therefore accumulate,
even though each individual MARK* run allows at most one mutation relative to
that step's input. Changing ligand order, selection criterion, or starting
structure can change the final result.

## Bounded ligand translation and rotation

Each single-position calculation adds
`osprey.c.confspace.StrandFlex.TranslateRotate(10, 2.0)` to the entire chain B
strand, matching the setup in `flexibility_study/IAS.py`. OSPREY exposes six
rigid-body degrees of freedom: each of the three rotation components is bounded
to ±10 degrees, and each of the three translation components is bounded to
±2.0 Å around the input pose. This is in addition to the target and selected
neighbors' residue-level flexibility. The protein-only conformational space
does not move the protein as a rigid body.

This enables a bounded pose search, not unrestricted docking. The input PDB
still defines the starting relative pose. Because these degrees of freedom
change the ligand and complex partition functions, their inclusion changes
the calculation and may substantially increase its cost. The single-position
and workflow manifests record that translation/rotation is enabled and its
bounds. Checkpoints generated before this setting was added are incompatible;
start a new trajectory from B1 rather than reusing a no-translation/rotation
manifest.

## Running

```bash
conda activate osprey-jdk17-cckstar-vtk93
cd /home/users/etm33/src/OSPREY3/MARKStar-work
python -u markstar_sequential.py --output sequential-results
```

The output directory must not already exist. Defaults:

| Option | Default | Meaning |
| --- | --- | --- |
| `--pdb` | supplied PKN2/SAV1 PDB | Initial complex. |
| `--output` | `sequential-output` | New workflow output directory. |
| `--max-flex` | 4 | Maximum **additional** flexible residues; the target is also flexible. |
| `--radius` | 4.0 | Inclusive CA-to-CA cutoff in Å. |
| `--cpus` | 4 | CPU count per single-position OSPREY calculation. |
| `--heap-mib` | 4096 | JVM heap limit per attempt. |
| `--epsilon` | 0.683 | Requested per-state partition-function approximation tolerance. |
| `--top-k` | 10 | Number of complex models retained per identity, subject to availability. |
| `--resume-manifest` | none | Existing workflow CSV from which to reuse completed steps. |
| `--resume-after` | none | Index of the last completed step to reuse; both resume options are required together. |

For example, `--radius 6 --max-flex 4` uses a 6 Å cutoff with at most five
flexible residues including the target. `--max-flex 0` is allowed. The
per-step job receives the same radius and limit throughout the workflow, but
the actual neighbors may change because selection uses that step's input
coordinates. Neighbors can be in either chain. Ties are resolved by chain and
residue-ID string; if fewer candidates meet the cutoff, all are used and the
shortage is recorded. For the original PDB at B1 only GLU B2 lies within 4 Å.

Each subprocess uses the same Python interpreter as the driver and creates a
new OSPREY JVM. The exhaustive complex-structure calculation scales as the
product of residue-conformation counts and includes the bounded ligand pose
degrees of freedom. Fifteen steps, potentially with two attempts per step and
up to five flexible neighbors, can be **very expensive**. `--max-flex 5`
permits six flexible positions including the target;
`--top-k` reduces saved models, not the number of conformations evaluated.
Do not assume the full workflow will finish quickly. Run long jobs under an
appropriate scheduler and preserve the output directory and logs.

## Restarting from a completed step

Stop the earlier process before resuming so that its manifest and handoff
structures are no longer being modified. Resume into a **new** output
directory, not the existing one. Only resume from a checkpoint produced with
the same input, scan script, and translation/rotation settings. Earlier
checkpoints without rigid-body translation/rotation cannot be resumed by the
current scanner. For example, to retain the first two completed steps of a
compatible run and recalculate B3 onward:

```bash
python -u markstar_sequential.py \
  --pdb PKN2_SAV1_forced_pilot_model_0.clean.pdb \
  --max-flex 5 --radius 4 --cpus 4 --heap-mib 4096 --epsilon 0.683 --top-k 10 \
  --resume-manifest prior-translate-rotate-run/workflow_manifest.csv \
  --resume-after 2 --output sequential-pkn2-flex5-resumed
```

`--resume-after` is a one-based **step index**, not a residue number. The
driver checks the full ordered target list and requires each reused row to be
`selected` or `skipped`. It verifies the original PDB hash, single-position
script hash, numerical settings, input/output SHA-256 hashes, and the handoff
chain for every reused step. A skipped step must carry forward the same input.
Different parameters, translation/rotation settings, or a modified scan
script are rejected rather than silently combining incompatible
calculations. The sequential driver itself may
change: each reused row retains its original driver hash, while new rows record
the current driver hash. These checks do not prove the physical validity of
the earlier results or prevent a later external modification of a file.

The new manifest retains the earlier rows and their original paths, and its
`resumed_from_manifest` field identifies the source for each reused row.
The first new step uses the final reused row's verified output; later source
rows and their artifacts are ignored. New attempts and handoffs are written
only in the new output directory. Original files are not copied, so **retain
the source directory** for provenance and for future resumes. Skipped or
warning-bearing reused rows still count toward the final warning exit status.
Explicitly choosing an earlier checkpoint discards later results from the new
trajectory, without deleting or changing the old ones.

### Checkpoint validation and provenance

The continuation performs these checks before creating its output directory:

1. Read the original PDB specified by `--pdb`, and collect the chain B target
   identifiers in file order. The source manifest must have exactly the same
   number of rows, and every row must have the expected one-based index and
   target, including rows after the chosen checkpoint.
2. For each row through `--resume-after`, require a `selected` or `skipped`
   status. `running`, `pending`, `failed`, and `interrupted` are not completed
   checkpoints. An earlier skipped step is permitted; its unchanged input
   remains part of the reused trajectory. Selecting a later checkpoint after
   such a skip does **not** retry that skipped step.
3. Check that the input PDB and `markstar_scan.py` hashes match those recorded
   in every reused row. Require the same radius, maximum neighbor count, CPU
   count, heap, epsilon, top-k, and ligand translation/rotation settings,
   comparing numeric values rather than their textual formatting.
4. Starting from the original PDB, check each row's input path and hash against
   the preceding output and its recorded hash. Compute the current SHA-256 of
   every referenced input and output file and compare it to its recorded hash.
   The skipped output must be the same path as its input. Missing, changed, or
   inconsistent handoffs cause an error rather than a silent restart.
5. Copy the completed row values into the new manifest and set
   `resumed_from_manifest` on those rows to the immediate source manifest.
   Keep their recorded attempts, warnings, selections, script hashes and
   absolute artifact paths. Create new pending rows for all later targets,
   with the current driver and scan script hashes; these new rows start with
   an empty `resumed_from_manifest`. New execution begins at the next row.

The source directory is not modified or copied, and previously generated
attempts after the checkpoint are not reused. This includes failed attempts:
the first new step starts a new all-identity attempt, even if an old attempt
at that same position failed. For a continuation of a continuation, the
immediate source manifest records the previous link; keep the entire chain
of source manifests and referenced artifacts.

The driver checks the scan script but does **not** hash or verify the OSPREY
JAR, Python interpreter, dependencies, operating system, or Java version.
Use the same environment for scientific comparability. The driver hash can
differ intentionally, so inspect it and the source history when interpreting
a mixed-version trajectory. SHA-256 is an integrity check, not a validation
of molecular geometry or numerical correctness. Avoid modifying either script,
the PDBs, or dependencies during a live run; subprocesses read the scan script
when launched, and the startup hashes do not lock these files.

### Example: November 2026 continuation

For the supplied PKN2/SAV1 input, the original Slurm job `13497969` was
cancelled on November 11, 2026. Its B1 and B2 rows were `selected` (ARG and
HIS respectively). B3 and B4 were `skipped` after null exceptions, and B5
was running when the job was cancelled. A `running` row left by cancellation
is not a completed result. The original directory and all its partial outputs
were preserved.

The first continuation, job `13529861`, reused B1-B2, skipped B3 after the
no-PRO retry failed score validation, selected B4=HIS, B5=ARG, and B6=ARG,
then started exporting B7 structures. B7's score file and partial exports
were written, but the workflow never completed or selected a B7 handoff. A
later Slurm restart attempt failed immediately because the output directory
already existed; the existing continuation directory was preserved.

A second continuation, job `13755368`, was submitted after adding the
finite-positive-bound policy described below. It resumed from the original
completed B2 checkpoint, so it recalculated B3 and subsequent positions; it
did not reuse B3-B7 results or alter the first continuation's files. It used
five **additional** maximum neighbors, a 4 Å cutoff, four CPUs, a 4096 MiB JVM
heap, epsilon 0.68, and top-k 10, but did not include rigid-body
translation/rotation. It was submitted to `compsci` with 4 CPUs, 8 GiB Slurm
memory, and a 30-day limit. The JVM heap is part of, not additional to, the
Slurm memory allocation. The job was canceled on November 12 while its B3
no-PRO retry was in the structure-export phase; B3 had not been selected. Its
partial output was preserved and is not a compatible checkpoint for the
translation/rotation run. Its directory and scheduler log are:

```text
sequential-pkn2-flex5-finite-bounds-20261112/
  workflow_manifest.csv
  step_003/attempt_01.log
  step_003/attempt_02.log          # expected if all-identity attempt fails
  ...
sequential-pkn2-flex5-finite-bounds-20261112.slurm.log
```

The scheduler log and attempt output remain useful for auditing the canceled
run, but are not evidence of a completed B3 or workflow. The replacement run
starts over from the original input at B1.

```bash
sacct -j 13755368 --format=JobID,State,ExitCode,Elapsed
tail -f sequential-pkn2-flex5-finite-bounds-20261112.slurm.log
tail -f sequential-pkn2-flex5-finite-bounds-20261112/step_003/attempt_02.log
```

### Current translation/rotation run

Job `13758376` (`markstar-pkn2-tr`) was submitted on November 12, 2026 after
the OSPREY setup-only check succeeded. It starts at **B1 / position 1** from
the original PDB; it does not resume or reuse any prior position. It enables
the chain B `TranslateRotate(10, 2.0)` flexibility and uses epsilon `0.683`,
five additional flexible neighbors at up to 4 Å, four CPUs, a 4096 MiB JVM
heap, and top-k 10. Slurm allocated 4 CPUs and 8 GiB on `compsci`, with a
30-day limit. Its output and scheduler log are:

```text
sequential-pkn2-flex5-tr10deg-2A-eps0683-20261112/
  workflow_manifest.csv
  step_001/attempt_01/
  step_001/attempt_02/           # PRO-excluded retry if the first attempt fails
sequential-pkn2-flex5-tr10deg-2A-eps0683-20261112.slurm.log
```

At the latest progress check, the job was `RUNNING` on
`compsci-cluster-fitz-37`. Its all-identity B1 attempt failed with the known
`Mutation from HID to PRO not supported yet` error; the automatic PRO-excluded
retry returned 19 sequence scores and was exporting structures. HID and ALA
exports were complete, and ARG export was at 100 of 306 conformations. B1 had
not yet been selected; B2-B15 remained pending. This is a progress snapshot,
not evidence that B1 or the full workflow has completed. Monitor it with:

```bash
squeue -j 13758376
sacct -j 13758376 --format=JobID,State,ExitCode,Elapsed
tail -f sequential-pkn2-flex5-tr10deg-2A-eps0683-20261112.slurm.log
tail -f sequential-pkn2-flex5-tr10deg-2A-eps0683-20261112/step_001/attempt_02.log
```

The scheduler log reports step-level progress and errors; attempt logs contain
OSPREY's detailed output. For the canceled job above, inspect
`step_003/attempt_02.log` to review its retry. For the current job, use the
attempt logs under the current step's directory. Consult that job's workflow
CSV for its actual step, selected residue, flexible residues, parameters, and
outcome; a log file's existence alone is not evidence of success. Job IDs are
scheduler-specific.

## Step-by-step behavior and decisions

1. Validate the original single-model PDB with the single-position parser;
   create `workflow_manifest.csv` with a pending row for each ligand position.
   The initial input is used for step 1.
2. Read the current PDB. Identify the target by its original chain B ID and
   choose the nearest CA neighbors in either chain, excluding the target.
   Configure each selected neighbor for continuous wild-type flexibility and
   configure the target to be flexible and mutable. Add bounded
   `TranslateRotate(10, 2.0)` flexibility to the entire chain B strand in the
   ligand and complex conformational spaces. The single-position manifest
   records OSPREY's assigned templates and rigid-body bounds. Histidine PDB names
   HIS/HID/HIE/HIP are treated as the same amino-acid identity when preserving
   the wild-type template; protonation variants are not separately enumerated.
3. Start an all-20-identity attempt in `attempt_01/`, recording the child
   process output in `attempt_01.log`. If the process exits nonzero, or the
   log reports OSPREY's `Mutation from ... to PRO not supported yet`, preserve
   that attempt and retry once in `attempt_02/` with `--exclude-proline`.
   This policy applies to **any** nonzero exit; it does not imply PRO caused
   the failure. The retry is a reduced scan that does not evaluate PRO at
   that step. If excluding PRO is incompatible with a current wild-type PRO,
   the retry can fail. A failed retry is logged as an error and normally marks the
   step skipped, carrying the prior structure forward. The exception is a
   nonzero exit after completed scoring/export with only the final score-quality
   validation error, as detailed below. Other partial files from failed attempts
   are not used to select a structure. Infrastructure failures may therefore
   cause a step to be skipped too: inspect both logs and warnings.
4. Once K* returns, `markstar_scan.py` writes `sequence_scores.csv` containing
   the full-precision OSPREY K* lower and upper bounds, estimated-score flag,
   state statuses and effective epsilons. This is separate from OSPREY's
   rounded `kstar.results.tsv`. Rank every usable upper bound, **highest
   first**, including finite positive bounds from unvalidated or nonconverged
   scores as requested. Exclude NaN, positive or negative infinity, zero,
   negative, and nonnumeric upper bounds. Excluded scores remain in the raw
   score CSV for audit but are not rankable candidates. Under this workflow's
   interpretation, these rows represent mutations without a usable K* value
   (for example, geometrically or energetically impossible mutations), rather
   than making all the other identities in the scan fail. Zero is also
   excluded because a usable K* ratio is positive.
   Comparisons use decimal arithmetic without converting to floating point;
   exact ties break alphabetically by OSPREY identity name. Finite positive
   upper bounds are eligible
   even when `score_estimated` is false or the requested epsilon was not met;
   those quality issues are retained as warnings. If K* did not return, the
   exact-score CSV is absent and partial TSV rows are not promoted to completed
   scores.
5. From the selected identity's multi-model PDB, choose model 1: the lowest
   calculated complex energy **within that identity**. The structure exporter
   exhaustively minimizes configured RC combinations; this is not a guarantee
   of a continuous global minimum. Do not rank energies across different
   mutations: the calculator includes identity-dependent reference terms.
   Remove `MODEL`/`ENDMDL` records to make a single-model `next_input.pdb`,
   retain atom and TER records, and validate the residue sequence and identity.
   The selected identity is the only allowed sequence change in this handoff.
6. Iterate candidates in descending finite-positive-upper-bound order. If a candidate
   has no exported PDB or its first model fails handoff validation, record a
   warning and try the next-ranked finite-positive candidate. Select the first candidate
   with a valid handoff. If no finite-positive candidate with a valid PDB remains, mark
   the step `skipped` and carry forward its input unchanged.
7. Proceed to the next position using the chosen handoff, or unchanged input
   after a skip. The final row's `output_pdb` is the final structure.

### Retry and failure decision table

| First attempt | Second attempt | Outcome |
| --- | --- | --- |
| Exits zero and no explicit unsupported-PRO message | Not run | Read exact scores and try to select the highest-upper-bound identity and its model. |
| Exits nonzero, including a Java exception or export failure, or reports the explicit unsupported-PRO message | Exits zero, no unsupported-PRO message | Use only the retry's exact scores and models; record both attempts and the retry warning. |
| Exits nonzero or reports unsupported PRO | Exits zero, but has score-quality warnings | Use finite upper bounds from the retry, even if some identities have no score, non-finite bounds, or fail epsilon validation. Preserve the warnings and omit non-finite candidates. |
| Exits nonzero or reports unsupported PRO | Exits nonzero specifically because scoring completed and structure export completed, then the scanner raised its final score-validation `RuntimeError` | Use finite upper bounds and exported structures from that completed attempt. The driver recognizes this narrow case by the attempt manifest phase/error and presence of both score and structure manifests. Preserve the nonzero exit and warnings. |
| Exits nonzero or reports unsupported PRO | Any other nonzero exit, or no usable finite candidate/structure | Log an error, mark the step `skipped`, and carry the step input forward. Do not select from incomplete or failed-at-computation outputs. |

This retry policy does not repair the MARK* implementation and does not
establish why an all-identity attempt failed. It repeats the computation with
PRO excluded even when an unrelated error caused the first failure. If the
current target is already PRO, `markstar_scan.py` rejects `--exclude-proline`
because it cannot remove the target's wild-type identity; if the retry fails,
the step is skipped like any other failed retry.

A finite positive upper bound is a rankable value, not a guarantee of convergence,
accuracy, or physical validity. The workflow retains `score_estimated`, state
statuses, epsilons, and score warnings so users can judge that uncertainty.
Non-finite upper bounds (NaN or either infinity), zero, and negative upper
bounds are not ranked and cannot win merely by sorting above finite values.
Missing structures are also excluded from actual selection: the workflow
continues down the finite-positive-bound ranking until it finds a valid
exported PDB. If none qualifies, it skips the step.

### Non-finite results and finite-bound selection (2026-11-12 update)

The raw `sequence_scores.csv` is preserved as OSPREY returned it. It may
contain `NaN`, `Infinity`, `-Infinity`, or zero upper bounds for individual
identities. The sequential driver parses the bound with `Decimal` and omits
all non-finite, zero, negative, or non-numeric upper bounds from
`candidate_scores` and selection. It records an omission warning; the raw row
is not deleted from the attempt's score file. A non-finite row therefore does
not, by itself, invalidate finite rows for other mutations.

Among remaining candidates, ranking uses the full-precision finite positive
upper bound, descending, with alphabetic identity as the deterministic
tie-break. It does **not** require `score_estimated=True` or an effective
epsilon below the requested threshold. Those are quality/provenance warnings,
not filters on a finite positive bound. Selection still requires an exported,
valid first-model PDB: if a higher-ranked candidate has no PDB or fails
handoff validation, the driver records a warning and tries the next-ranked
candidate. A finite score without a valid model cannot be used as a handoff.
If there is no candidate satisfying both conditions, the step is skipped and
its input is carried forward.

There is a narrow allowance for the single-position scanner's final
score-validation exception. After K* returns, the scanner writes
`sequence_scores.csv`, exports structures, and then raises
`RuntimeError: Incomplete or insufficiently precise results` when any score
fails its quality checks. The subprocess therefore exits nonzero even though
scores and structures may be complete. The sequential driver treats this as
completed output **only** if all of these are true:

1. The subprocess exit code is nonzero.
2. The per-attempt manifest phase is `exporting_structures`.
3. Its error begins with the exact score-validation RuntimeError text.
4. Both `sequence_scores.csv` and `structures/manifest.tsv` exist.

It still records the nonzero exit, the failed child manifest status, and all
score-quality warnings; it does not relabel the scan as validated or
converged. Only finite positive bounds and corresponding valid PDBs can be
selected. Other nonzero exits—including MARK* Java exceptions, energy-matrix
failures, and incomplete structure export—remain failures and their partial
outputs are not used for selection. This is intentionally narrow, so a new or
different failure mode is not mistaken for completed scoring.

In the preserved B3 retry from November 12, 2026,
`step_003/attempt_02/sequence_scores.csv` contained 19 non-PRO results. PHE,
TRP, and TYR had zero upper bounds with `Unstable` complex statuses and no
estimated K* scores; they are excluded by the positive-bound rule. The
highest finite positive upper bound with a present, structurally valid export
was **B3=ILE**, `64973374916819.84491311726991892520585831783635471488907790923085`.
Its `ILE.pdb` first model passed the handoff residue/identity validation in a
read-only check. The original continuation did not record or use this result:
it had the former policy that rejected every nonzero retry exit. Under the
updated policy, a new continuation can resume from the completed B2 checkpoint
and recalculate B3 onward; the existing failed manifest and old artifacts
remain unchanged.

A zero exit is not necessarily a useful result. If no finite positive exact scores are
available, no candidate has a valid exported model, or handoff validation
fails for every candidate, the step is skipped. Unexpected exceptions in the
driver itself, such as a failure to open an attempt log or parse a malformed
score file, are recorded as `failed` when possible and stop the workflow; they
are not caught and converted into a retry or skip. A forced scheduler kill can
prevent even that update.

## Warnings and interpretation

Ranking by the largest upper bound is an explicitly requested **optimistic**
rule. An upper bound alone is not a reliable estimate of binding affinity:
wide finite bounds can win, even if a different sequence has a better
validated score. Non-finite upper bounds are excluded. Continuing after score validation failures or OSPREY
energy-bound warnings is also an explicit choice; it does not make the
selected prediction reliable. Each attempt log, validation problem, exit code,
candidate's state statuses and epsilons, and captured bound warning are recorded
or linked from the manifest. The log-warning match is not exhaustive; inspect
the full logs, including Java exceptions. A nonzero child exit normally
triggers the retry even if it wrote scores or partial structure exports. The
exception is a narrowly detected final score-validation failure after both
K* returned and structure export completed: finite bounds and valid exported
models from that attempt can be used, while all validation problems and its
nonzero exit remain in the manifest. Other failed attempts are never used for
selection. A successful child may still have missing or invalid structures;
the driver then tries the next-ranked finite-bound identity and skips only if
none has a valid handoff.

The supplied environment has already shown an N-terminal HIS-to-PRO exception,
an ASP complex effective epsilon above the requested threshold, and an ARG
energy-bound warning. These are not fixed by the sequential driver. Results
should be considered exploratory until those issues are understood.

The driver exits **1** after completing all steps if any step was skipped or
reported a warning, including a proline retry or shortage of neighbors. It
exits 0 only if all steps finish without recorded warnings. Unexpected errors
or interruptions also exit nonzero. A skipped step is not treated as a
successful structure optimization. Thus a Slurm job can finish its scan and
still report a nonzero exit; inspect the manifest and logs to distinguish
completed-with-warnings from an interrupted or failed workflow.

## Output layout and manifest

```text
sequential-results/
  workflow_manifest.csv
  step_001/
    attempt_01.log
    attempt_01/run_manifest.csv
    attempt_02.log                  # when the first attempt fails or PRO is unsupported
    attempt_02/run_manifest.csv
    attempt_02/sequence_scores.csv
    attempt_02/kstar.results.tsv
    attempt_02/structures/<IDENTITY>.pdb
    attempt_02/structures/manifest.tsv
    next_input.pdb                  # only if a model was selected
  step_002/
  ...
```

The top-level CSV is updated atomically as work proceeds. It has one row per
original ligand position. `status` is `pending`, `running`, `selected`,
`skipped`, `failed`, or `interrupted`. Each row includes the input/output
absolute PDB paths and SHA-256 hashes, target, settings, start/update times,
initial PDB and script hashes, selected identity and unrounded upper bound,
translation/rotation enablement and bounds, selected multi-model source PDB,
and skip/failure reason. Reused rows also identify `resumed_from_manifest`.
`flexible_residues`, `attempts`, `candidate_scores`,
and `warnings` are JSON arrays inside CSV cells. An attempt in progress is
recorded before the subprocess is launched. The per-attempt run manifest has
additional details including exact selected distances and OSPREY templates;
the per-structure manifest maps models to energy and RC assignments.

`attempts` records each attempt's directory, log path, exit code, whether PRO
was excluded, whether the unsupported-PRO text was seen, the child manifest
status and score problems, and captured warning lines. `candidate_scores`
contains only finite, positive upper-bound candidates from the selected
attempt, including finite candidates with unestimated scores or excessive
effective epsilon. It never contains NaN or infinite bounds. The selected
attempt may have a nonzero exit only in the documented completed-output
score-validation case. `reason` explains a skip or driver failure;
`warnings` retains retry, child-exit, score, and selected warning information.
An attempt marked running may have only its directory and exclusion flag,
because its subprocess has not returned yet. `output_pdb` for a skipped step
is the unchanged input path, not a newly created copy.

Atomic replacement keeps each CSV snapshot parseable but is not a journal. A
forced kill or machine failure may leave a step `running` or an attempt
incomplete; `pending` rows are not evidence of success. Neither the CSV nor
the file timestamp guarantees all filesystem data survived a power loss.
Existing output directories are never overwritten. Resumption creates a new
directory and requires an explicit completed checkpoint as described above.

## Verification

```bash
python -m unittest -v test_markstar_scan test_markstar_sequential
```

Tests cover selection at different positions, histidine aliases, unrounded
upper-bound ranking and ties, invalid bounds, retry, skipped-step carry-forward,
manifest updates, checkpoint reuse and hash/parameter rejection, failed retries,
successful non-PRO retries following a generic failure, rejection of partial
scores from a failed retry, and actual PDB first-model extraction when the existing
example export is available. The next-position OSPREY setup was also tested
on the extracted real HID model. A separate real first-step integration run
used `--max-flex 0 --top-k 1`: the all-identity attempt failed on PRO, the
19-identity retry completed, the highest full-precision upper bound selected
ARG, and the resulting PDB built the B2 OSPREY spaces. The test stopped
intentionally before B2 scoring; it was not a full 15-position run or a test of
the default four-neighbor setting. Tests with mocked jobs do not validate the
physical predictions, performance, or correctness of every OSPREY mutation.
The November continuation's real B2 checkpoint, hashes and settings were
validated before submission, and B3 was observed starting; the new retry
policy has not yet been validated through a completed real B3 retry.
