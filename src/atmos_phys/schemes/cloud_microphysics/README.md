# Cloud and microphysics process schemes

This group exposes six CAM cloud and microphysics processes through independent
CCPP-style `ap_*_scheme` modules while preserving the legacy public CAM entry
points.

| Legacy entry | Scheme entry | Host facade |
|---|---|---|
| `nucleate_ice_cam_calc` | `nucleate_ice_cam_calc_run` | `nucleate_ice_cam` |
| `dropmixnuc` | `dropmixnuc_run` | `ndrop` |
| `ice_macro_tend` | `ice_macro_tend_run` | `macrop_driver` |
| `cldfrc` | `cldfrc_run` | `cloud_fraction` |
| `mmacro_pcond` | `mmacro_pcond_run` | `cldwat2m_macro` |
| `micro_mg_tend` (MG 1.0) | `micro_mg_tend_run` | `micro_mg1_0` |

The five stateful source modules are relocated as a unit and the original
module name is retained as a thin use-associated facade. This is intentional:
initialization, namelist state, private helper state, and the run routine all
remain in one module instance rather than being split between a host module
and a scheme module.

`ice_macro_tend` remains elemental. Its legacy elemental routine is a scalar
adapter to `ice_macro_tend_run`; the `ap_ice_macro_tend_run` timer therefore
wraps the two array-valued call sites in `macrop_driver` and `clubb_intr`
instead of making the elemental numerical kernel impure or timing every scalar
element.

The active PI configuration dispatches MG version 1.0 and Park macrophysics.
Ice supersaturation adjustment is retained on both the Park and CLUBB paths,
but is inactive when `micro_do_icesupersat` is false. MG 1.5 and MG 2.0 remain
unchanged.

This source checkpoint provides static compile and metadata validation only.
It must be integrated with the rest of its process group and pass the required
50-step BFB run before the group is marked BFB.
