# Internal Validation Notes

This file records machine-specific validation details for the current project
workspace. It is useful for the internal team, but public forks can replace or
remove it.

## Dashboard

- Status dashboard: `https://cam.idapro.me/`
- Dashboard repository: `/glade/u/home/ruitong/code/cam-codon-status`
- CAM repository: `/glade/u/home/ruitong/code/iCESM1.3.1_fzhu/components/cam`

Use a deployment-specific token when writing status:

```bash
cd /glade/u/home/ruitong/code/cam-codon-status
export CAM_STATUS_TOKEN='<token>'

uv run cam-codon-status remote-mark \
  --remote https://cam.idapro.me \
  --relpath src/physics/cam/example.F90 \
  --routine example_subroutine \
  --kind subroutine \
  --status processing \
  --note 'started Codon validation'
```

## 2026-06-16 PI/MCO 6-Month All-Codon Validation

Validation context:

- Selector snapshot: `675 codon / 0 native`
- Compiler: GNU
- Queue route: Derecho `develop` to `cpudev`
- Source state: after the UWSHCU positive-moisture expression-order fix
- Work root: `/glade/u/home/ruitong/cam_pi_mco_long_bfb_20260615`

Results:

| Case | Length | Jobs | Result | Main timing |
| --- | ---: | --- | --- | --- |
| PI, `ne16_g16`, startup | 6 months | baseline `6467097.desched1`, all-Codon `6467103.desched1` | `overall_numeric_equal=True` | `CPL:RUN_LOOP` 6443.491 -> 7608.488, +18.080% |
| MCO, `ne16_g16`, hybrid restart from `MCO/restart/2001-01-01-00000` | 6 months | baseline `6467105.desched1`, all-Codon `6467112.desched1` | `overall_numeric_equal=True` | `CPL:RUN_LOOP` 3826.046 -> 4744.456, +24.004% |

Evidence:

```text
/glade/u/home/ruitong/cam_pi_mco_long_bfb_20260615/results/20260616_pi_mco_6month_allcodon_after_uwshcu/compare_pi_6month_pristine_vs_allcodon/compare.txt
/glade/u/home/ruitong/cam_pi_mco_long_bfb_20260615/results/20260616_pi_mco_6month_allcodon_after_uwshcu/manual_compare_mco_6month_pristine_vs_allcodon/compare.txt
/glade/u/home/ruitong/cam_pi_mco_long_bfb_20260615/worklog.md
```

Run directories:

```text
/glade/derecho/scratch/ruitong/pi_mco_long_bfb_20260615/runs/pi_6month_pristine_baseline_20260616_pi_mco_6month_allcodon_after_uwshcu/run
/glade/derecho/scratch/ruitong/pi_mco_long_bfb_20260615/runs/pi_6month_allcodon_20260616_pi_mco_6month_allcodon_after_uwshcu/run
/glade/derecho/scratch/ruitong/pi_mco_long_bfb_20260615/runs/mco_6month_pristine_baseline_20260616_pi_mco_6month_allcodon_after_uwshcu/run
/glade/derecho/scratch/ruitong/pi_mco_long_bfb_20260615/runs/mco_6month_allcodon_20260616_pi_mco_6month_allcodon_after_uwshcu/run
```

Notes:

- The automated long-run driver timed out while `mco_6month_allcodon` was still
  running, but the PBS job later completed successfully and was manually
  compared.
- The authoritative MCO result for this validation is the manual compare file
  listed above.
