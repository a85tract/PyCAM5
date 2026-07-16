# Vertical-diffusion process group

This directory contains the B05 CCPP-style process boundaries extracted from
the PI-atm `tphysac` path.  Each leaf directory owns one `<process>_run`
entry and its metadata:

- `compute_tms`: turbulent mountain-stress numerical kernel.
- `calc_obklen`: surface kinematic fluxes and Obukhov length.
- `compute_eddy_diff`: UW moist-turbulence eddy-diffusivity process.
- `compute_vdiff`: implicit vertical-diffusion solver and selector state.
- `positive_moisture`: post-diffusion moisture positivity correction.

The original CAM modules remain as compatibility adapters.  They preserve the
existing public APIs and own the `ap_<process>` timers; the positivity adapter
also retains CAM logging.  `compute_eddy_diff_run` calls the adapter for
`compute_vdiff`, so both parent and child process timers are exercised on the
iterative UW path.
