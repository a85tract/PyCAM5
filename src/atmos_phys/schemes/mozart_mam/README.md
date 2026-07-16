# MOZART/MAM schemes

This group contains small portable decisions extracted from MOZART/MAM while
CAM-owned aliases remain in their original modules.  In particular,
`drydep_update_run` returns only whether the land dry-deposition velocity must
be updated.  The pointer association
`lnd(lchnk)%dvel => cam_in%depvel` remains in `mo_drydep.F90`; it is never
replaced by an array copy.

The group also contains the numerical process boundary for the chemistry
solvers used by the PI-atm `pp_trop_mam3` mechanism.

- `exp_sol_run` is the real active explicit-solver entry.  The generated
  mechanism has `clscnt1=0`, so its numerical result is intentionally a no-op;
  the adapter timer and entry counter still prove that the process ran.
- `trop_mam3/imp_sol_run` owns the generated 20-species implicit solve,
  including the original Newton, Jacobian/LU, convergence and adaptive-step
  statement order.  Its generated helper procedures are private to the scheme
  module.

CAM initialization, species/reaction name lookup, history registration and
`outfld`, logging/abort handling, and the existing public chemistry-driver
interfaces remain in `mo_exp_sol` and `mo_imp_sol` adapters.

The `modal_aero` subgroup owns the MAM gas/aerosol and dry-deposition process
boundaries.  Its historical modules under `src/chemistry/modal_aero` remain
thin CAM API adapters:

- `gasaerexch/aero_model_gasaerexch_run` preserves the parent process order
  and calls the separately timed exchange, nucleation, and coagulation scheme
  entries.
- `gasaerexch/modal_aero_gasaerexch_sub_run`,
  `newnuc/modal_aero_newnuc_sub_run`, and
  `coag/modal_aero_coag_sub_run` own their original numerical modules and
  initialization state.
- `drydep/aero_model_drydep_run` owns dry deposition together with its private
  particle deposition-velocity helper.
