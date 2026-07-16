# Deep-convection schemes

This group contains the numerical implementations extracted from the legacy
CAM physics hosts while preserving their public CAM entry points.

- `dry_adiabatic_adjust/dadadj_scheme.F90`: dry adiabatic adjustment.
- `zhang_mcfarlane/zm_convr.F90`: Zhang-McFarlane deep-convection driver and
  its private numerical helper chain.
- `zhang_mcfarlane/zm_conv_evap.F90`: convective precipitation evaporation.
- `zhang_mcfarlane/zm_conv_convtran.F90`: convective tracer transport.
- `zhang_mcfarlane/zm_conv_momtran.F90`: convective momentum transport.

The compatibility adapters remain in `src/physics/cam/dadadj.F90` and
`src/physics/cam/zm_conv.F90`. Each adapter preserves the original argument
order and wraps its scheme call with an independent `ap_*_run` timer.
