# MARK* source changes and runtime-JAR verification

**Verification date:** 2026-11-11  
**Repository:** `/home/users/etm33/src/OSPREY3`  
**Runtime environment:** `osprey-jdk17-cckstar-vtk93`

This note records which MARK* instrumentation and bound-handling changes are
present in the OSPREY source checkout, whether those changes were present in
the OSPREY runtime used by the MARK* workflow, and the evidence behind that
conclusion. It distinguishes ordinary progress output from new source-level
APIs, which is important when deciding whether a source edit was compiled and
loaded.

## Conclusion

The local source checkout contains uncommitted MARK* edits, but the inspected
runtime JAR does **not** contain those edits. In particular:

- The three new `MARKStarBound` counter/search-space getters are absent from
  the JAR.
- The new `MARKStarNode.setExactEnergyBound()` overloads are absent from the
  JAR.
- The JAR's full-conformation minimization bytecode calls the old
  `setBoundsFromConfLowerAndUpper()` method, not `setExactEnergyBound()`.
- An earlier run log from this environment contains lower-bound inversion
  diagnostics consistent with the old monotonic update behavior.

Therefore, the source-level exact-leaf correction was not active in the
inspected MARK* runtime. The ordinary `MARK* expanded` and `partial
minimizations` progress lines do not change this conclusion; they are existing
engine output, not proof that the new getters were compiled.

This was a verification of the JAR available in the named conda environment.
It is not proof of the identity of every JAR ever used on the cluster. The
workflow scripts launch Python in this environment, and their logs identify
OSPREY 3.3-dev, but future runs should record the actual loaded JAR path and
hash in their manifest to make this provenance explicit.

## Source changes in the checkout

The relevant local source files are:

- [`MARKStarBound.java`](../../src/main/java/edu/duke/cs/osprey/markstar/framework/MARKStarBound.java)
- [`MARKStarNode.java`](../../src/main/java/edu/duke/cs/osprey/markstar/framework/MARKStarNode.java)
- [`TestMARKStarNodeBounds.java`](../../src/test/java/edu/duke/cs/osprey/markstar/TestMARKStarNodeBounds.java)

The MARK* source edits were uncommitted working-tree changes when inspected.
They were not part of the checkout's committed `HEAD` at that time.

### Counter and search-space getters

`MARKStarBound.java` adds:

- `getNumInternalNodesProcessed()`
- `getNumPartialMinimizations()`
- `getNumConfsInSearchSpace()`

These expose existing counters and root-node state. They do not add counting
logic, print output, change search ordering, or change MARK* bounds. They are
API/instrumentation hooks: a caller could query them, but the MARK* console log
does not report them merely because the methods exist.

### Exact-leaf energy-bound replacement

The behavior-changing source edit is in the full-conformation minimization
path. After `ConfAnalyzer` returns a minimized energy, `MARKStarBound` calls
`MARKStarNode.setExactEnergyBound()` instead of the previous monotonic
`setBoundsFromConfLowerAndUpper()` update.

The new method is guarded so it can only be used for a fully assigned leaf
whose conformation count is exactly one and whose assignments contain no
unassigned positions. It replaces both energy endpoints with the accepted
energy and recomputes the single conformation's Boltzmann contribution for
both partition-function endpoints.

The reason for this change was a reported case in a CFKStar flexibility study:

- Previous energy lower bound: `-20.029944216241258` kcal/mol.
- Accepted minimized energy: `-20.04029111337419` kcal/mol.

The minimized energy was lower than the bound that was supposed to be a lower
bound. The pre-existing generic setter only tightened bounds monotonically:
it kept the old energy lower bound but lowered the energy upper bound to the
minimized energy. That produced an inverted energy interval. Since the
Boltzmann transformation decreases with energy, the corresponding pfunc lower
bound could then exceed its upper bound.

The source change permits the fully minimized leaf's energy interval to be
replaced by the accepted minimized value. This restores an ordered interval
for that leaf under the code's convention that this accepted minimization
result is the value assigned to the completed conformation. It does not
explain why the prior A* lower bound was too high, prove every partial-node
bound is sound, or fix other implementations of MARK*. The existing
`checkConfLowerBound()` diagnostic logs the contradiction but does not abort.

The regression test [`TestMARKStarNodeBounds.java`](../../src/test/java/edu/duke/cs/osprey/markstar/TestMARKStarNodeBounds.java)
tests the setter on a constructed complete, single-conformation node. That is
a focused test of the interval-replacement behavior; it is not an end-to-end
MARK* calculation or a proof of the energy model's physical assumptions.

## Runtime JAR inspected

The MARK* Python workflow imports the installed `osprey` package in this
environment. The inspected JAR was:

```text
/home/users/etm33/miniconda3/envs/osprey-jdk17-cckstar-vtk93/lib/python3.11/site-packages/osprey/lib/osprey-3.3.jar
```

At inspection it reported:

| Property | Value |
| --- | --- |
| OSPREY runtime banner | `OSPREY 3.3-dev` |
| Python | `3.11.15` |
| Java | `17.0.18` |
| JAR file timestamp | `2026-03-12 12:51:11 -0400` |
| JAR size | `12,083,594` bytes |
| JAR SHA-256 | `89174c1b2492b18924a725c9b959ed036b3c0cc3f69f001230b64ddc7b54b082` |

The MARK* wrapper in the installed
`osprey/__init__.py` constructs `MARKStarBoundFastQueues`, which extends
`MARKStarBound`. The wrapper is at
`/home/users/etm33/miniconda3/envs/osprey-jdk17-cckstar-vtk93/lib/python3.11/site-packages/osprey/__init__.py`.
The JAR inspection below was performed on the installed classes from that
runtime, not on class files built from the current source checkout.

## How compilation was checked

The class signatures were inspected with `javap` against the installed JAR:

```bash
JAR=/home/users/etm33/miniconda3/envs/osprey-jdk17-cckstar-vtk93/lib/python3.11/site-packages/osprey/lib/osprey-3.3.jar

javap -p -classpath "$JAR" \
  edu.duke.cs.osprey.markstar.framework.MARKStarBound

javap -p -classpath "$JAR" \
  edu.duke.cs.osprey.markstar.framework.MARKStarNode

javap -p -classpath "$JAR" \
  'edu.duke.cs.osprey.markstar.framework.MARKStarNode$Node'
```

The added getters were not listed by `javap`. Neither the outer
`MARKStarNode` class nor its nested `Node` class listed
`setExactEnergyBound()`.

The bytecode for `MARKStarBound`'s full-conformation minimization lambda was
also inspected:

```bash
javap -p -c -classpath "$JAR" \
  edu.duke.cs.osprey.markstar.framework.MARKStarBound
```

That code invokes:

```text
MARKStarNode.setBoundsFromConfLowerAndUpper(double, double)
```

It does not invoke `setExactEnergyBound()`. This is direct evidence that the
inspected JAR retains the old call site.

## What the run logs show—and do not show

The original workflow and its continuation wrote MARK* progress output such
as:

```text
MARK* expanded: ..., partial minimizations: ..., ...
Processing ... internal nodes...
```

These lines report progress already produced by the engine. They are not
output from the newly added getters, and their presence does not indicate
that those getters exist in the loaded bytecode.

The earlier run's B1 reduced-scan log,
`sequential-pkn2-flex5-20261110/step_001/attempt_02.log`, contains diagnostics
of this form:

```text
Bounds are incorrect:-14.75655869454229 > -14.769396140086041
Updating conf lower bound of  -14.75655869454229 with -14.769396140086041, which is lower!?
```

The warning reports a minimized energy below the previous lower bound. The
subsequent “which is lower!?” messages are consistent with the generic
monotonic setter refusing to lower that bound. The behavior agrees with the
old bytecode call found by `javap`, rather than the local exact-leaf setter.

These diagnostics were produced by the earlier attempt in the same workflow
environment. They corroborate the JAR inspection but should not be confused
with a fresh controlled run comparing two built versions.

## What this establishes, and what it does not

Established for the inspected environment:

1. The source checkout had the described local MARK* changes.
2. The installed JAR did not contain the new getters or exact-leaf setter.
3. The installed JAR's full-leaf path called the old setter.
4. Existing run output included the lower-bound diagnostic pattern expected
   from that old behavior.

Not established:

1. That a rebuilt JAR containing the source change would fix all MARK* bound
   problems. The modification addresses an observed full-leaf inversion, not
   the underlying source of an inadmissible earlier lower bound.
2. That the accepted minimized value is always an exact physical minimum.
   The update relies on OSPREY's interpretation of the accepted minimization
   result as the completed leaf's value.
3. That alternate partition-function implementations, such as the separate
   asynchronous implementation, received the same correction. The inspected
   local async path still used the old setter.
4. That all past or future cluster jobs use this exact JAR. A job may use a
   different environment or JAR if its launch configuration changes.
5. That ordinary progress lines alone can establish which source revision
   produced a binary.

The original and resumed MARK* workflows use the named conda environment.
Unless that environment's JAR was rebuilt/replaced, the source correction is
not active in those jobs. A source change only affects a runtime after the
appropriate Java artifacts are rebuilt and the runtime is pointed at the new
artifact.

## Reproducible verification commands

From the repository root, compare the local source diff:

```bash
git diff -- \
  src/main/java/edu/duke/cs/osprey/markstar/framework/MARKStarBound.java \
  src/main/java/edu/duke/cs/osprey/markstar/framework/MARKStarNode.java
```

Check the installed runtime JAR hash and query its methods:

```bash
JAR=/home/users/etm33/miniconda3/envs/osprey-jdk17-cckstar-vtk93/lib/python3.11/site-packages/osprey/lib/osprey-3.3.jar
sha256sum "$JAR"
javap -p -classpath "$JAR" edu.duke.cs.osprey.markstar.framework.MARKStarBound
javap -p -classpath "$JAR" edu.duke.cs.osprey.markstar.framework.MARKStarNode
javap -p -classpath "$JAR" 'edu.duke.cs.osprey.markstar.framework.MARKStarNode$Node'
```

For future provenance, record the resolved JAR path and its SHA-256 in each
workflow's manifest at process startup. A useful stronger check is a small
runtime smoke test that loads the class from the intended environment and
asserts the expected method is present. Neither check replaces scientific
validation of the bound algorithm.
