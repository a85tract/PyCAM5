# Modal aerosol schemes

This group exposes the B04 modal-aerosol process boundaries while preserving
the existing CAM host APIs:

- `modal_aero_depvel_part_run`
- `modal_aero_calcsize_sub_run`
- `modal_aero_wateruptake_dr_run`
- `aero_model_wetdep_run`
- `wetdepa_v2_run`

The original chemistry modules retain registration, initialization, pbuf index
ownership, history/framework setup, and compatibility entry points. Their run
adapters add `ap_<process>` timers and call these scheme modules. Private
numerical helpers required by a process live with that process.
