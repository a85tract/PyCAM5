# Chemistry-emission schemes

This group contains the array kernels extracted from the active PI-atm
surface-emission path.  CAM adapters retain derived `physics_state`/`cam_in`
access, tracer lookup, tracer-data pointers, unit strings, grid/time queries,
history output, and performance timers.

The original execution order is unchanged: dust emission is applied first,
the lowest-layer wind is adjusted to the 10 m sea-salt emission factor, sea
salt is accumulated into the shared constituent-flux slots, and prescribed
file emissions are added afterward by the MOZART chemistry adapter.

`sslt_sections_private.F90` is a private numerical dependency of
`seasalt_emis_run`; its 31-section initialization and accumulation order are
copied from the PI-atm `sslt_sections` implementation.  It has no CCPP entry
of its own.
