# PyCAM5 / iCAM5_iHESP CAM Codon Port

## TL;DR

This repository contains an experimental Codon-backed port of selected CAM5
routines inside the isotope-enabled iCESM1.3/iHESP CAM component.

The native Fortran CAM implementation remains the reference. The Codon
implementation is enabled routine by routine through runtime environment
variables such as `*_IMPL=codon`. The main goal is to move computational
kernels into Codon while preserving bit-for-bit (BFB) CAM output against a
pristine native baseline.

As of the 2026-06-16 validation snapshot, the tracked selector set used for
validation contains 675 runtime `*_IMPL` switches, and two 6-month
production-style validations on the Derecho HPC system have achieved
`overall_numeric_equal=True` against matching native baselines.

## Terms Used Below

- CAM5: the Community Atmosphere Model version used by this CESM/iCESM tree.
- Native path: the original CAM Fortran implementation.
- Codon path: the translated implementation compiled from `*_codon.py` into a
  shared library.
- Selector: a runtime environment variable, usually named `*_IMPL`, that chooses
  `native` or `codon` for one entry point.
- BFB: bit-for-bit equality. In this project, BFB means the compare script
  reports `overall_numeric_equal=True`.
- PI and MCO: the internal pre-industrial and Miocene validation cases used for
  long-run testing.

## Project Goals

- Preserve BFB CAM output relative to a pristine native Fortran baseline.
- Move computational CAM kernels into Codon while keeping the CESM/CAM calling
  interface stable.
- Allow routine-by-routine rollout through runtime selectors.
- Track progress with routine-level execution evidence rather than touched-file
  counts.
- Keep the native Fortran path available for comparison, fallback, and
  numerically fragile expression islands.

## Current Status

Snapshot: 2026-06-16 validation runs, after the UWSHCU positive-moisture
expression-order fix.

- Selector coverage: 675 tracked runtime selectors in the validation snapshot.
- Long-run BFB evidence: PI (pre-industrial) and MCO (Miocene) 6-month
  all-Codon runs both compare with `overall_numeric_equal=True` against matching
  pristine native baselines.
- Progress tracking: routine status is maintained in the CAM Codon status
  dashboard, with separate states for complete, partial, in progress, and not
  started routines.
- Reference implementation: native Fortran remains authoritative.

The selector count is a snapshot, not a permanent API. It may change as new CAM
entry points are added or split.

## Architecture: Native CAM + Runtime-Selectable Codon Layer

This port does not replace the CESM case workflow. It adds a runtime-selectable
Codon layer beside the native CAM Fortran implementation.

- Original Fortran routines remain in place.
- Selected routines have Fortran wrappers with `bind(C)` interfaces that call
  Codon shared libraries.
- Codon implementations live in `*_codon.py` files and are compiled into
  `lib*_codon.so`.
- Runtime selectors choose the implementation:

  ```bash
  EXAMPLE_IMPL=native
  EXAMPLE_IMPL=codon
  ```

- CAM logs print execution proof lines such as `implementation = codon` or
  `direct = codon`.
- The BFB rule is strict: a Codon path is accepted only when the compare script
  reports `overall_numeric_equal=True`.

A native island is a minimal Fortran expression or statement block intentionally
kept in the native path because translating that specific operation to Codon was
shown to break proven BFB behavior.

## Repository Layout

Important Codon-port files live beside the original CAM source:

- `src/physics/cam/*_codon.py`: Codon implementations for CAM physics kernels.
- `src/chemistry/*/*_codon.py`: Codon chemistry and aerosol helpers.
- `src/dynamics/se/*_codon.py`: Codon spectral-element dynamics helpers.
- `src/utils/cam_misc_codon.py`: shared CAM utility helpers.
- `cime_config/buildlib`: CIME build hook that compiles Codon sources during
  `case.build`.
- `scripts/cam_codon_guard.py`: validation/pre-commit guard for high-risk
  numerical changes.
- `.codon_guard.yaml`: guard policy, validation metadata, forbidden generated
  artifacts, and metadata-only compare differences.
- `doc/internal_validation.md`: internal dashboard, workspace, and validation
  artifact locations for this project.

## Setup

This CAM tree is intended to sit inside a full iCESM/CESM checkout:

```text
iCESM1.3.1_fzhu/
  cime/
  components/
    cam/        # this repository
    cice/
    clm/
    pop/
```

Requirements:

- A working CESM/CIME environment for the target machine.
- A Fortran compiler and MPI stack supported by the case.
- Codon available as `codon` in `PATH`, or at `~/.codon/bin/codon`.
- Python for helper scripts and compare tooling.
- NetCDF/CIME runtime dependencies required by the parent CESM case.

Build a case normally through CIME:

```bash
cd /path/to/case
./case.setup
./case.build
```

During `./case.build`, `components/cam/cime_config/buildlib` finds the Codon
compiler and builds the CAM Codon shared libraries. If Codon is missing, the
build fails before model execution.

Codon libraries are built with floating-point contraction disabled:

```bash
codon build -release -lib --relocation-model=pic --fp-contract=off --global-ctor=no
```

`--fp-contract=off` is required for BFB because fused multiply-add contraction
changes rounding relative to the native Fortran baseline.

## Quick Start

This assumes you already have a configured CESM case using this CAM source tree.

Build the model:

```bash
cd /path/to/case
./case.setup
./case.build
```

Run one Codon-enabled path:

```bash
env MICROP_DRIVER_IMPL=codon ./case.submit --skip-preview-namelist
```

Check that the Codon path executed:

```bash
zgrep -n 'implementation = codon\|direct = codon' /path/to/run/atm.log.*.gz
```

For BFB validation, compare the run against a pristine native baseline produced
with the same case configuration:

```bash
python /path/to/compare_cesm_runpair.py \
  --native-run-dir /path/to/pristine/native/run \
  --codon-run-dir /path/to/codon/run
```

The expected pass condition is:

```text
overall_numeric_equal=True
```

## Running With Native Vs Codon Implementations

Selectors are ordinary environment variables read by the Fortran wrappers:

```bash
env MICROP_DRIVER_IMPL=codon ./case.submit --skip-preview-namelist
env MICROP_DRIVER_IMPL=native ./case.submit --skip-preview-namelist
```

For all-Codon validation, source or generate a selector file that sets every
known selector to `codon`:

```bash
set -a
source /path/to/env_allcodon_675.sh
set +a
./case.submit --skip-preview-namelist
```

For a mixed run, leave known non-BFB or intentionally native selectors as
`native` and set the remaining selectors to `codon`.

For isolated validation or parallel work, keep each lane separate:

- one CAM source tree or worktree
- one case root
- one `EXEROOT`
- one fresh `RUNDIR`
- lane-local Codon library paths in `Macros.make`

## BFB Validation Workflow

Every validation must compare against a pristine native baseline for the same
case configuration. Case settings that can affect outputs, such as compset,
grid, restart mode, orbit settings, domain files, compiler, and runtime length,
require a matching baseline.

Typical short validation settings:

```bash
./xmlchange \
  STOP_OPTION=nsteps,STOP_N=50,REST_OPTION=nsteps,REST_N=50, \
  BFBFLAG=TRUE,DOUT_S=FALSE,TIMER_DETAIL=6,TIMER_LEVEL=16, \
  CONTINUE_RUN=FALSE
```

After the job completes, prove execution and compare:

```bash
rg -n '<SELECTOR_NAME>' /path/to/case/logs/run_environment.txt.*
zgrep -n 'implementation = codon\|direct = codon' /path/to/run/atm.log.*.gz

python /path/to/compare_cesm_runpair.py \
  --native-run-dir /path/to/pristine/native/run \
  --codon-run-dir /path/to/codon/run
```

Only `overall_numeric_equal=True` is accepted as BFB. Character metadata
differences such as `time_written`, `cpath`, `nfpath`, and `nhfil` are expected
when all numeric variables match.

Do not reuse an old run directory for BFB proof unless old model outputs and
logs have been removed and the cleanup is recorded.

## Progress Tracking

Routine status is tracked outside this repository by the CAM Codon status
dashboard. The internal project deployment is listed in
`doc/internal_validation.md`.

The dashboard tracks which routines are complete, in progress, partial, or not
started, and provides routine pages, formula/equation pages, coverage-case
views, CSV/HTML exports, and REST APIs such as `/api/summary` and
`/api/routines`.

Status meanings:

- `done`: the default active path enters Codon and returns for the same routine.
- `done-native-island`: the active path is Codon except for a minimal native
  expression or statement block retained for proven BFB reasons.
- `processing`: someone is actively editing or validating the routine.
- `partial`: Codon covers helper islands or some branches, but the default
  routine body still has native orchestration.
- `none`: covered by the current case snapshot but no Codon evidence yet.
- `unknown`: parser or evidence conflict that needs manual review.

Use exact `(relpath, routine, kind)` keys when updating status. Example:

```bash
cd /path/to/cam-codon-status
export CAM_STATUS_TOKEN='<token>'

uv run cam-codon-status remote-mark \
  --remote https://<dashboard-host> \
  --relpath src/physics/cam/example.F90 \
  --routine example_subroutine \
  --kind subroutine \
  --status processing \
  --note 'started Codon validation'
```

After validation:

```bash
uv run cam-codon-status remote-mark \
  --remote https://<dashboard-host> \
  --relpath src/physics/cam/example.F90 \
  --routine example_subroutine \
  --kind subroutine \
  --status done \
  --note 'commit <sha>; selector EXAMPLE_IMPL=codon; proof atm.log line; overall_numeric_equal=True'
```

Do not mark a routine `done` just because a helper library exports a related
symbol. The proof must show that the same routine, wrapper, or accepted
same-routine dispatch path actually executed.

## Production-Style Validation Results

The following long validations were run on the Derecho HPC system with GNU
builds after the UWSHCU positive-moisture expression-order fix. Both used all
selectors set to Codon (`675 codon / 0 native`) and compared against matching
pristine native baselines.

| Case | Length | Jobs | Result | Main timing |
| --- | ---: | --- | --- | --- |
| PI pre-industrial case, `ne16_g16`, startup | 6 months | baseline `6467097.desched1`, all-Codon `6467103.desched1` | `overall_numeric_equal=True` | `CPL:RUN_LOOP` 6443.491 -> 7608.488, +18.080% |
| MCO Miocene case, `ne16_g16`, hybrid restart from `MCO/restart/2001-01-01-00000` | 6 months | baseline `6467105.desched1`, all-Codon `6467112.desched1` | `overall_numeric_equal=True` | `CPL:RUN_LOOP` 3826.046 -> 4744.456, +24.004% |

The internal compare outputs and run directories are listed in
`doc/internal_validation.md`.

These runs are validation examples, not a claim that every future compset,
compiler, restart state, or production campaign is automatically BFB.

## Development Rules

- Keep source changes and generated artifacts separate. Do not commit
  `__pycache__`, run directories, build directories, case logs, compare outputs,
  `.so` files, or guard receipts.
- Use fresh run directories for validation proof.
- Record job id, selector counts, run directory, `atm.log`, `run_environment`,
  `END OF MODEL RUN`, compare output, and timing.
- Preserve Fortran floating-point expression order when translating to Codon.
- If a routine is BFB only with a small native expression island, document the
  expression and validation evidence before marking it `done-native-island`.

## Known Limitations

- Native Fortran remains the reference implementation.
- BFB has been proven only for specific case configurations and matching native
  baselines.
- New compsets, grids, compilers, restart states, runtime lengths, and namelist
  changes require separate validation.
- Some numerically fragile expressions intentionally remain in native Fortran.
- The Codon layer is not a general drop-in replacement for all CAM
  configurations.
- Performance is not yet the primary optimization target; the documented
  long-run validations are slower than native Fortran.

## Troubleshooting

- `codon` not found during `case.build`: install Codon or add it to `PATH`.
  The build also checks `~/.codon/bin/codon`.
- Codon selector is set but no proof line appears in `atm.log`: the selected
  path may not execute in this case configuration, or the routine may still
  dispatch through native orchestration.
- Compare reports non-BFB: first verify the baseline uses the same case
  configuration, compiler, restart state, orbit settings, and runtime length.
- Scattered one-ULP differences: check for FMA contraction, changed operation
  order, `pow` lowering, complex `sqrt`, or reduction-order differences.
- A run looks complete but compare is suspicious: confirm output timestamps,
  job ids, `run_environment`, and `END OF MODEL RUN` all match the submitted
  validation run.
- Dashboard counts look inconsistent: use `/api/summary` for totals before
  interpreting paginated `/api/routines` results.
