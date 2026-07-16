# Orographic gravity-wave and surface-resistance schemes

This directory contains the B05 leaf-scheme extraction for five PI-atm
processes:

- `gw_prof_run`
- `gw_oro_src_run`
- `gw_drag_prof_run`
- `energy_change_run`
- `calcram_run`

The original public CAM routines remain as thin host adapters. They retain
their original argument lists, start and stop an independent
`ap_<process>_run` timer, and pass CAM module state or `GWBand` fields to
the corresponding scheme through explicit arguments. This avoids a module
cycle between `gw_common` and the extracted schemes while preserving the
original numerical statement order.

Each scheme has a matching CCPP metadata table. Full B05 acceptance still
requires the integration lane's 50-step BFB comparison.
