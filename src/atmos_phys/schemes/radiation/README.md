# RRTMG radiation schemes

This group owns the four production radiation process boundaries used by the
PI-atm `trop_mam3` configuration:

- `rad_rrtmg_sw_run` and `rrtmg_sw_run` in `rrtmg_shortwave/`;
- `rad_rrtmg_lw_run` and `rrtmg_lw_run` in `rrtmg_longwave/`.

The historical `radsw`, `radlw`, `rrtmg_sw_rad`, and `rrtmg_lw_rad` modules
are compatibility facades.  They re-export the scheme procedures under CAM's
original names, so callers in `radiation.F90` are unchanged.  Moving the full
modules keeps `radsw_init`, `radlw_init`, the saved solar-band irradiance, the
saved longwave top level, and the private RRTMG helpers in exactly one module
instance.

The active low-level implementation is the McICA source formerly located in
`src/physics/rrtmg/ext/rrtmg_mcica`.  The same-basename sources under
`ext/rrtmg_sw` and `ext/rrtmg_lw` are generic alternatives shadowed by the
production Filepath and were deliberately not changed.  Scheme source files
use unique basenames so Filepath lookup cannot silently select an alternative.

Each exported run has one balanced production timer:

- `ap_rad_rrtmg_sw_run`;
- `ap_rad_rrtmg_lw_run`;
- `ap_rrtmg_sw_run`;
- `ap_rrtmg_lw_run`.

The shortwave driver's only early return (`Nday == 0`) is inside its numerical
core.  Its outer run wrapper therefore always reaches the matching stop call.
CCPP metadata is generated from the annotated run interfaces and checked
against the Fortran declarations with the framework's
`offline_check_fortran_vs_metadata.py` utility.  Runtime acceptance still
requires the separate integrated 50-step BFB gate.
