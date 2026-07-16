# S03 gas-phase chemistry audit evidence

Commit: `621f922` (`Decouple active gas-phase chemistry orchestration`)

Status: intermediate increment only.  The commit closes the common transforms,
generated rate adjustments, and explicit/implicit solver subgraph, but it does
not yet satisfy the final one-call gas-phase process boundary.  The host still
orchestrates other science providers, so `gas_phase_chemdr_run` must not be
marked completed or `dependency_pass` from this evidence.

Scope: active `src/chemistry/pp_trop_mam3` mechanism only.  This commit does
not move or modify the photolysis support closure assigned to the separate
photolysis lane.

## Host/science boundary

- CAM registration, physics-buffer access, history, timers, logging, and
  fatal handling remain in `src/chemistry/mozart/mo_gas_phase_chemdr.F90`.
- Intrinsic-only common state transforms are in
  `src/atmos_phys/support/chemistry/common/chemistry_common_kernels.F90`.
- Active generated mechanism operations are in
  `src/atmos_phys/support/chemistry/gas_phase/trop_mam3_gas_phase_kernels.F90`.
- Standalone explicit/implicit orchestration is rooted at
  `src/atmos_phys/schemes/mozart_mam/gas_phase_chemdr/gas_phase_chemdr.F90`.

## Active-path order audit

The host keeps the original active PI order:

1. constituent map and `mmr -> vmr` conversion;
2. invariant and stratospheric-state preparation;
3. generated thermal rates (`setrxt` equations);
4. `sulf_interp`, `qsat`, `usrrxt`, and optional GHG rate adjustment;
5. generated reaction-rate adjustment (`adjrxt` equations);
6. column setup and active/table photolysis;
7. `O1D_to_2OH_adj` and the active mechanism `phtadj` no-op;
8. external forcing and heterogeneous/washout rates;
9. explicit solver phase;
10. implicit solver phase;
11. O3S and H2SO4 post-solver updates;
12. modal aerosol gas/aerosol exchange;
13. boundary conditions, `vmr -> mmr`, tendencies, and dry-deposition flux map.

The following equations were compared statement-by-statement with their
active source implementations and preserve operand and assignment order:

- `mo_setrxt:setrxt` -> `trop_mam3_set_thermal_rates`;
- `mo_adjrxt:adjrxt` -> `trop_mam3_adjust_rates`;
- `mo_phtadj:phtadj` -> `trop_mam3_adjust_photolysis` (empty active body);
- `mo_mean_mass:set_mean_mass` -> `chemistry_set_mean_mass`;
- `mo_mass_xforms:mmr2vmr`, `vmr2mmr`, `mmr2vmri`, and `h2o_to_vmr` ->
  corresponding common kernels;
- in-driver constituent map, tendency formation, dry-deposition flux map,
  negative clamp, O3S update, and H2SO4 before/after update -> standalone
  common/gas-phase kernels.

## Dependency closure

Commands:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 tools/check_atmos_phys_standalone.py \
  --scope all --process gas_phase_chemdr_run --mode direct
PYTHONDONTWRITEBYTECODE=1 python3 tools/check_atmos_phys_standalone.py \
  --scope all --process gas_phase_chemdr_run --mode transitive
```

Results:

```text
[direct] PASS: entries=1 roots=1 closure_sources=4 violations=0 duplicates=0
[transitive] PASS: entries=1 roots=1 closure_sources=4 violations=0 duplicates=0
```

## Metadata

Command:

```sh
python3 /glade/u/home/ruitong/code/CAM-SIMA/ccpp_framework/scripts/fortran_tools/offline_check_fortran_vs_metadata.py \
  --directory /glade/derecho/scratch/ruitong/tmp/tphys_standalone_20260716/s03/cam/src/atmos_phys/schemes/mozart_mam/gas_phase_chemdr
```

Result: `All checks passed!`

## Fresh ifx compilation

Compiler:

```text
/glade/u/apps/derecho/25.10/spack/opt/spack/ncarcompilers/1.1.0/intel-oneapi-compilers/2025.2.1/d7a7/bin/ifx
```

Fresh output directory:

```text
/glade/derecho/scratch/ruitong/tmp/tphys_standalone_20260716/s03_check_gas_common_final
```

Standalone compile order, all PASS:

```text
chemistry_common_kernels.F90
trop_mam3_gas_phase_kernels.F90
exp_sol.F90
trop_mam3/imp_sol.F90
gas_phase_chemdr.F90
```

CAM host syntax compile order with the S01 build module tree, all PASS:

```text
src/chemistry/mozart/mo_exp_sol.F90
src/chemistry/pp_trop_mam3/mo_imp_sol.F90
src/chemistry/mozart/mo_gas_phase_chemdr.F90
```

This is compile/static evidence for an intermediate increment only.  It is
not process-completion evidence and is not a 50-step or BFB result; runtime
validation is recorded only after the integrated batch gate runs.
