# Atmospheric physics utilities

This leaf directory is reserved for shared, host-independent utilities used by
the decoupled atmospheric-physics schemes.  B01 does not require a shared
Fortran utility yet; its scheme implementations live under `schemes/`.

B01 keeps two host-specific boundaries deliberately narrow:

- `apply_constituent_tendencies_run` receives one active constituent as a
  rank-2 field.  The outer constituent loop and the following CAM minimum
  checks stay in the adapter, preserving their original interleaving.
- `qneg4_run` returns the excess-flux indices and values.  CAM water-tracer
  lookup, `wtrc_ratio`, and isotope correction stay in the adapter because
  those operations depend on CAM's tracer registry; the generic surface-flux
  limiter remains in the scheme.
