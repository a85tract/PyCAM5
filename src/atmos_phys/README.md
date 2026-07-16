# PI-atm `tphysac` / `tphysbc` modular decoupling

This tree tracks the extraction of the active PI-atm physics path from the two
host drivers in `src/physics/cam/physpkg.F90`.  It follows the CESM atmospheric
physics split without importing the full CCPP runtime:

- a scheme owns array-level scientific calculations and exposes the applicable
  `<scheme>_init`, `<scheme>_timestep_init`, `<scheme>_run`,
  `<scheme>_timestep_final`, and `<scheme>_final` entry points;
- a CAM adapter owns `physics_state`, `physics_ptend`, `pbuf`, namelist, MPI,
  history, and error translation;
- `tphysac` and `tphysbc` retain host orchestration and the original call order;
- suite XML files document the intended active ordering.  Scheme names in XML
  are CCPP-style base names; the lifecycle selects the exact entry point in the
  JSON inventory.

The active campaign has 69 exported entry points in seven dependency-coherent
batches.  B01 contains nine entries and B02--B07 contain ten each.  Branches
that cannot execute in the current PI case are explicitly listed as
`deferred_inactive`; they do not count toward the active BFB denominator.

## Progress source of truth

`decoupling_status.json` is the only editable source of truth.
`DECOUPLING_PROGRESS.md` is generated and must not be edited manually.

The fixed active states are:

```text
planned -> in_progress -> build_pass -> 50step_running -> bfb
                                                  \-> failed
```

Use the standard-library-only tracker from the CAM repository root:

```bash
python3 tools/tphys_decoupling_progress.py check
python3 tools/tphys_decoupling_progress.py render --check

python3 tools/tphys_decoupling_progress.py update \
  --process dadadj_run --status in_progress \
  --note 'worker w06 started the B06 extraction'

python3 tools/tphys_decoupling_progress.py update \
  --batch B06 --status 50step_running \
  --commit <full-40-character-source-commit>
```

`update --status bfb` deliberately fails unless the entries already reference
an evidence-complete BFB run.  The normal BFB transition is `record-run`:

```bash
python3 tools/tphys_decoupling_progress.py record-run \
  --run-id B06-20260715-job123 \
  --batch B06 --result bfb --job-id 123.desched1 \
  --source-commit <full-40-character-source-commit> \
  --case-root <isolated-case-root> \
  --run-dir <fresh-50step-run-dir> \
  --baseline-run-dir <matching-pristine-native-run-dir> \
  --numeric-equal true --char-equal true \
  --executable <cesm.exe> \
  --run-environment <run_environment.txt.job123> \
  --filepath <Buildconf/camconf/Filepath> \
  --namelist <Buildconf/camconf/atm_in> \
  --compare-report <compare-result.txt> \
  --output-file <cam-history-or-restart-file> \
  --proof 'dadadj_run=<atm.log-or-timing-file>::<matched execution line>' \
  --proof 'zm_convr_run=<atm.log-or-timing-file>::<matched execution line>' \
  ...one --proof for every entry in B06...
```

A BFB record is rejected unless all of the following are present:

- an exactly 50-step run;
- at least one output file;
- `numeric_equal=true` and `char_equal=true`;
- execution proof for every entry in the batch;
- executable, run-environment, CAM Filepath, namelist, and comparison report
  files.  Their paths, sizes, and SHA-256 hashes are recorded where applicable.

`char_equal=true` refers to normalized scientific character data.  NetCDF
`date_written` and `time_written` are volatile wall-clock provenance fields and
are excluded; every other character variable must compare exactly.

A failed run can also be recorded with `--result failed`; it requires a clear
`--note`.  This keeps failed 50-step attempts visible without allowing them to
be mistaken for BFB evidence.

## Development checks

```bash
python3 -m unittest discover -s tools/tests -p 'test_tphys_decoupling_progress.py'
python3 tools/tphys_decoupling_progress.py check
python3 tools/tphys_decoupling_progress.py render --check
```

The tracker validates the exact 69-entry inventory, batch sizes, dependency
references, suite coverage, process/run references, and every stored BFB
evidence record.
