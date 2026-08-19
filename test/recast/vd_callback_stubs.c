/* The three CAM callbacks vertical_diffusion_codon.py declares.
 *
 * They exist so the Codon library links standalone. The routine under test
 * (vertical_diffusion_ptend_core_codon) calls none of them -- they belong to
 * eddy_diff_caleddy and compute_cubic -- so each aborts rather than
 * returning a plausible number: a stub that silently answered would let a
 * future gate on those routines pass on fiction.
 */
#include <stdio.h>
#include <stdlib.h>

static double refuse(const char *name) {
    fprintf(stderr, "recast gate: %s is a CAM callback with no stand-in; "
                    "the gated routine must not call it\n", name);
    abort();
}

double compute_cubic_native_cb(double a, double b, double c) {
    (void)a; (void)b; (void)c;
    return refuse("compute_cubic_native_cb");
}

double eddy_diff_estblf_cb(double t) {
    (void)t;
    return refuse("eddy_diff_estblf_cb");
}

double eddy_diff_svp_to_qsat_cb(double es, double p) {
    (void)es; (void)p;
    return refuse("eddy_diff_svp_to_qsat_cb");
}
