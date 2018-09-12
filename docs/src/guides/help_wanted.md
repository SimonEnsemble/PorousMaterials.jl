## Help wanted and needed
* the speed of a GCMC or Henry simulation is determined primarily by how fast `PorousMaterials.jl` can compute the electrostatic and vdw potential energy. Some core functions that can speed up this are:
    * `nearest_image!` in `src/NearestImage.jl`
    * Ewald sums in `src/Electrostatics.jl`. (electrostatics are a huge bottleneck.)
    * `src/VdWEnergetics.jl`
    * The scripts `test/vdw_timing.jl` and `test/ewald_timing.jl` time the functions for benchmarking.

* consolidate `eikar`, `eikbr`, `eikcr` somehow without slowing down the Ewald sum.
* more tests added to `tests/runtests.jl`, `tests/henry_tests.jl`, `tests/gcmc_tests.jl`
* geometric-based pore size calculations (largest free and included spheres), surface area, and porosity calculations that take `Framework`'s as input
* handle .cif's without P1 symmetry. i.e. convert any .cif to P1 symmetry
* extend `gcmc_simulation` to handle mixtures
* better default rules for choosing Ewald sum parameters? alpha, kvectors required...
* Henry coefficient code prints off Ewald sum params 5 times if run with one core...
* set good defaults for `gcmc_simulation` probabilities (as now) but also allow user to change through default arguments to the function
* automatically adjust the translation step `δ` in `gcmc_simulation` during burn cycles to have 50% acceptance of translation moves (online gradient descent?)
* EQEq or other charge equilibration schemes for assinging charges, taking a `Framework` as input.
