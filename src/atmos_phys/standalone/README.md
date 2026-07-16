# Standalone physics phase 2

This directory owns infrastructure for making 36 selected PI-atm process
entries independently buildable.  It does not replace or modify the completed
69-entry decoupling campaign.

## Source of truth

- `standalone_status.json` is the only editable progress source.
- `STANDALONE_PROGRESS.md` is generated; never edit it by hand.
- `tools/tphys_standalone_progress.py` validates and updates the ledger.

The four groups contain exactly nine entries each:

| Group | Purpose | Dependencies |
|---|---|---|
| S01 | scientific dependency foundation | none |
| S02 | boundary-layer and convection host extraction | S01 |
| S03 | MAM/chemistry lower-level kernels | S01 |
| S04 | wet removal, cloud, and radiation drivers | S01, S03 |

Statuses advance through:

```text
planned -> in_progress -> dependency_pass -> build_pass -> 50step_running -> bfb
                                  \-> failed
```

`dependency_pass` means both direct and transitive checks pass.  The dependency
checker treats that state and all later states as `completed-only`.

## Dependency boundary

Run the incremental and final gates from the CAM repository root:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 tools/check_atmos_phys_standalone.py \
  --scope completed-only --mode both

PYTHONDONTWRITEBYTECODE=1 python3 tools/check_atmos_phys_standalone.py \
  --scope all --mode both
```

For a worker group or one process, add `--batch S01` or repeat `--process`.
Direct mode inspects only each ledger source.  Transitive mode recursively
follows every module implemented under `src/atmos_phys`.

While a group is still `planned` or `in_progress`, explicitly select it from
the full inventory; otherwise `completed-only` correctly skips it:

```bash
python3 tools/check_atmos_phys_standalone.py \
  --scope all --mode both --batch S01

python3 tools/tphys_standalone_progress.py update \
  --batch S01 --status dependency_pass \
  --commit <full-40-character-source-commit> \
  --note 'direct and transitive standalone dependency reports passed'
```

The only external modules allowed across the standalone boundary are:

- ISO/IEEE intrinsic modules;
- `shr_*` infrastructure;
- MPI (`mpi` or `mpi_*`), OpenMP (`omp_*`), NetCDF, and PIO.

Every other module must have exactly one provider under `src/atmos_phys`.
Modules found elsewhere in CAM, unresolved non-allowlisted modules, and
duplicate module definitions are errors.  In particular `mpishorthand`,
`physics_types`, `ppgrid`, `perf_mod`, history, pbuf, namelist, logging, and
CAM registries are not hidden allowlist exceptions.

The source closure used by CMake can be inspected without compiling:

```bash
python3 tools/check_atmos_phys_standalone.py \
  --scope completed-only --mode transitive --list-sources
```

## Standalone CMake gate

The default configuration needs only CMake and Python.  It checks completed
entries without enabling a Fortran compiler:

```bash
cmake -S src/atmos_phys/standalone -B /path/to/fresh/check-build
cmake --build /path/to/fresh/check-build --target check-standalone
ctest --test-dir /path/to/fresh/check-build --output-on-failure
```

The final 36-entry dependency gate is:

```bash
cmake --build /path/to/fresh/check-build --target check-standalone-all
```

To compile the checked transitive closure in a fresh standalone build tree:

```bash
cmake -S src/atmos_phys/standalone -B /path/to/fresh/fortran-build \
  -DTPHYS_STANDALONE_BUILD_SOURCES=ON \
  -DTPHYS_STANDALONE_SCOPE=completed-only \
  -DTPHYS_STANDALONE_MODULE_INCLUDE_DIRS='/standalone/dependency/modules' \
  -DTPHYS_STANDALONE_EXTERNAL_LIBRARIES='/standalone/dependency/lib/libshr.a'
cmake --build /path/to/fresh/fortran-build --target tphys-standalone-schemes
```

All atmos-physics `.mod` files are emitted into the fresh CMake build's
`modules/` directory.  Include directories containing `bld`, `Buildconf`, or
an ancestor `cesm.exe` are rejected at configure time.  This prevents a
nominally standalone build from silently succeeding by borrowing CAM EXEROOT
modules.  External SHR/MPI/OpenMP/NetCDF/PIO dependencies must come from a
separate installation or standalone dependency build.

## Acceptance sequence

For each group:

1. Move CAM glue to the host adapter and close direct plus transitive module
   dependencies under `src/atmos_phys`.
2. Record `dependency_pass` only after both checker modes pass.
3. Build the closure with this CMake project in a fresh directory and record
   `build_pass` only after its module manifest/build report is saved.
4. Run a fresh matched native/standalone 50-step case and require per-entry
   nonzero execution proof plus numeric and normalized-character BFB.
5. Use `record-run`; direct `update --status bfb` is rejected.  A BFB record
   additionally requires hashed dependency and standalone-build reports.

The evidence-bearing transition has the same process-scoped form as phase 1,
with two additional mandatory artifacts:

```bash
python3 tools/tphys_standalone_progress.py record-run \
  --run-id S01-<date>-job<id> --batch S01 --result bfb \
  --job-id <id>.desched1 --source-commit <full-source-commit> \
  --case-root <isolated-case> --run-dir <fresh-50step-run> \
  --baseline-run-dir <matched-native-run> \
  --numeric-equal true --char-equal true \
  --dependency-report <direct-and-transitive-report> \
  --standalone-build-report <fresh-cmake-build-report> \
  --executable <cesm.exe> --run-environment <run-environment> \
  --filepath <Buildconf/camconf/Filepath> --namelist <atm_in> \
  --compare-report <exact-comparison-report> \
  --output-file <cam-output> \
  --proof '<entry>=<timing-file>::<nonzero timer line>'
```

Standard tracker checks are:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 tools/tphys_standalone_progress.py check
PYTHONDONTWRITEBYTECODE=1 python3 tools/tphys_standalone_progress.py render --check
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest \
  tools.tests.test_tphys_standalone_progress \
  tools.tests.test_check_atmos_phys_standalone
```
