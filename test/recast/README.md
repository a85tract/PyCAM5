# Module-level bit-exactness gates, driven by RecastEngine

This repository's standing validation is whole-model: six-month PI and MCO
run pairs on Derecho, compared with `overall_numeric_equal=True`. That is
the strongest evidence there is, and it is also expensive — it answers "did
anything drift" hours after the fact, for the model as a whole.

This directory adds the complementary check at the other end of the cost
scale: **one routine, seconds, on a laptop**, with the same acceptance bar.
It runs the Codon implementation and the native Fortran side by side on
generated inputs and compares every output bit for bit, using
[RecastEngine](https://github.com/a85tract/RecastEngine)'s verification
stages rather than a script of its own.

The two are not redundant. A long run says the model agrees; a module gate
says *which routine* agrees, catches a regression in the minute it is
introduced, and leaves a machine-readable manifest saying under exactly
which compiler, flags, and inputs the claim holds.

## What runs

    frontend  (recast-cesm `cesm`)   analyzes dadadj.F90 as it sits in the tree
    candidate (adopted)             the Codon library, behind a ctypes shim
    oracle    (f2py-golden)         dadadj_native compiled by gfortran, with
                                    stand-ins supplying the framework modules
    gate      (differential.bitexact)
                                    same inputs to both sides, every output
                                    compared bit for bit, evidence written

Nothing numeric lives here. `dadadj_shim.py` allocates scratch and passes
pointers; every value is computed by the Codon code under test.
`standins.f90` supplies `ppgrid`, `physconst`, and the rest so the routine
compiles standalone — `cappa` is spelled as the same float64 expression on
both sides, because a constant folded two ways is not a bit-exact test of
anything else.

`_PREPARE_INPUTS` in the shim shapes inputs into the routine's defined
domain: pressure monotone with `pint` bracketing `pmid`, `pdel` the
thickness `pint` implies, and a temperature column near the dry adiabat.
Random columns make the fixed-point iteration diverge, and the native path
aborts there — comparing two error paths tests nothing. Both sides receive
the same shaped arrays, so the shaping cannot bias the verdict; it only
chooses which region of the input space is sampled.

## Running it

```bash
# 1. the engine and its CESM extension
pip install -e ../RecastEngine'[fortran,translate,verify]' -e ../recast-cesm

# 2. build the Codon library from the repository's own source
codon build --relocation-model=pic -release --lib \
    -o test/recast/build/libdadadj_codon.dylib src/physics/cam/dadadj_codon.py
cp "$(dirname "$(which codon)")"/../lib/codon/lib{omp,codonrt}.dylib test/recast/build/

# 3. run the gate (needs gfortran)
python test/recast/run_gate.py
```

Expected:

```
differential.bitexact: bit_exact
detail: 9600 points across 1 subprogram(s), all bit-exact
  dadadj_native: {'points': 9600, 'bit_exact': 9600, 'max_ulp': 0, ...}
evidence: file://.../test/recast/evidence/fortran_dadadj/<digest>.json
```

The manifest carries the artifact digest, the oracle's cache key (source
hash + compiler version + flags), the environment, and the full ULP
histogram — enough for someone else to reproduce the same verdict, which is
what separates evidence from a claim.

## Adding a routine

Copy the three files and change what they name: a shim exporting the
routine's surface with a `_SIGNATURES` table, whatever stand-ins the routine
use-imports, and the driver's `subprograms` / `dims` / `ranges`. A routine
whose native path has no side channels needs no `_PREPARE_INPUTS`; one that
iterates to convergence probably does.

Routines whose Codon path is selected at runtime through `*_IMPL` are all
candidates. `dadadj` was chosen first because it is small, iterative, and
numerically fragile — the shape of routine where a port most easily drifts.

## Three-way comparison: engine translation, Codon port, native Fortran

`run_threeway.py` answers a different question from the gate above. This
repository's long runs already establish that the **Codon port** matches the
native Fortran. The two open questions are about the *other* modernization
path — the rule-driven translator in RecastEngine:

    engine translation  vs  native Fortran    does the mechanical translator
                                              reproduce CAM?
    engine translation  vs  Codon port        do the two independent
                                              modernizations agree with each
                                              other?

Result for `vertical_diffusion_ptend_core_native`, 20 trials over 8×30×4
arrays:

```
 engine vs native Fortran: bit_exact   43200/43200 points, 0 ULP
  codon vs native Fortran: bit_exact   43200/43200 points, 0 ULP
 engine vs codon (direct): bit-identical on 43200 points across 120 arrays
```

All three sides derive from the same routine text. The routine is extracted
**byte-for-byte** from `src/physics/cam/vertical_diffusion.F90` into
`vd_ptend_core.f90`, because that file use-imports seventeen framework
modules and this routine needs none of them; the driver re-derives the
extraction on every run and refuses to proceed if the copy has drifted, so
the isolation cannot quietly become a rewrite. Only the environment is
stood in for — `pcols`, `pver`, `pcnst` as parameters.

`vertical_diffusion_ptend_core_native` was chosen as the first three-way
because it is **pure**: no module state read or written, no calls. There is
nothing to initialize, so no way for the three sides to disagree about state
none of them was given. A survey across this repository's 32 Codon-ported
files found 49 such routines, so the same harness extends by changing which
name it points at.

`vd_callback_stubs.c` exists only so the Codon library links standalone:
`vertical_diffusion_codon.py` declares three CAM callbacks
(`compute_cubic_native_cb` and two `eddy_diff_*`) that belong to other
routines. Each stub aborts rather than returning a plausible number — a stub
that answered silently would let a future gate on *those* routines pass on
fiction.

### What this does not claim

The inputs are physically ranged but synthetic, not sampled from a model
run. A bit-exact verdict here says the three implementations agree on the
sampled region of the input space; it does not say the region the model
actually visits was covered. That is what captured dumps are for, and it is
the natural next step for this harness.
