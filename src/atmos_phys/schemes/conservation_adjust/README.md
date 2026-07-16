# Conservation-adjustment schemes

This directory contains the process boundaries extracted from the PI-atm
`tphysac`/`tphysbc` path.  CAM-specific orchestration such as physics-buffer
access, global reductions, logging, constituent lookup, and timers remains in
the adapters under `src/physics/cam`.

## PI-atm compatibility entries

`check_energy_scaling_run` and `dycore_energy_consistency_adjust_run` are
explicit lifecycle compatibility entries.  The PI-atm baseline has neither
the newer energy-scaling algorithm nor a separate dycore-energy adjustment at
these call sites.  Therefore scaling returns an identity field and the dycore
entry leaves the physics tendency unchanged while returning zero local
adjustment.  They are exported, called, and timed so the decoupled lifecycle is
complete, but they must not be reported as active PI-atm physical algorithms.

Any change from this compatibility behavior requires its own scientific
review.  Validation status is tracked outside the source tree and may be
marked BFB only after an exact 50-step comparison passes.
