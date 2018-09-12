## Input files to describe crystals, molecules, and forcefields

All input files are stored in the path `PorousMaterials.PATH_TO_DATA` (type into Julia). By default, this path is set to be in the present working directory (type `pwd()` into Julia) in a folder `data/`. Go inside `PorousMaterials.jl/test/data` to see example input files for each case below.

There will be example code snippets through the documentation showing how to load in various files. To get a feel for this we have included a Tutorial Mode in `PorousMaterials.jl` that sets the `PorousMaterial.PATH_TO_DATA` to the data folder in our testing directory. To follow along with the examples without downloading your own data simply do the following:

```julia
using PorousMaterials

set_tutorial_mode()
    ┌ Warning: PorousMaterials is now in Tutorial Mode. You have access to the testing data to experiment with PorousMaterials.
    │ To get access to your own data use: reset_path_to_data()
    └ @ PorousMaterials ~/git_files/PorousMaterials.jl/src/PorousMaterials.jl:75
```

#### Atomic masses

Add fancy pseudo-atoms to `data/atomic_masses.csv`.

#### Peng-Robinson gas parameters

Critical temperatures and pressures and acentric factors are stored in `data/PengRobinsonGasProps.csv`.
