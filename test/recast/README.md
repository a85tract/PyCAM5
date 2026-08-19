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
