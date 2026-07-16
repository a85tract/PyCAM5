# UW shallow-convection process group

This directory owns the B06 CCPP-style process boundaries for
`compute_uwshcu_inv_run` and `compute_uwshcu_run`.  The numerical bodies,
private helpers, initialized module state, namelist reader, and initializer
move together so that the original operation order and shared state are
preserved.  `src/physics/cam/uwshcu.F90` remains as a compatibility module
that exports the original CAM names.

The inverse-coordinate process calls the main UW process through its explicit
`_run` boundary, so both `ap_compute_uwshcu_inv_run` and
`ap_compute_uwshcu_run` timers execute on the production path.
