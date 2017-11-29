# PorousMaterials.jl

A software package in [Julia](https://julialang.org/) for classical molecular modeling in porous crystals.

For Julia to find the PorousMaterials.jl package, add the directory where the source code 
is found to the `LOAD_PATH` variable in Julia, for example:

```julia
push!(LOAD_PATH, homedir() * "/Dropbox/PorousMaterials.jl/src/")
```

Example use:

```julia
using PorousMaterials

# load abstract data structure containg crystal structure info
framework = read_crystal_structure_file("JAVTAC_clean.cif")

# replicate framework to .xyz for visualization
replicate_to_xyz(framework, "replicated_framework.xyz", repfactors=(2,3,1))

# load abstract data structure Lennard Jones force field (UFF)
forcefield = read_forcefield_file("UFF.csv") 

```

All input files are stored in `PorousMaterials.PATH_TO_DATA`, which by default is 
`pwd() * "/data/"`. The user can change this by editing `PorousMaterials.jl`.

## Crystal

Crystal structure files are stored in `PorousMaterials.PATH_TO_DATA * "crystals/"`. Crystals 
in .cif format must be in P1 symmetry. If your .cif is not in P1 symmetry, our function
`convert_cif_to_P1_symmetry()` calls the Atomic Simulation Environment (ASE) in Python to 
write the .cif in P1 symmetry.

For example:

```
# crystal structure files in `PATH_TO_DATA/crystals`.
convert_cif_to_P1_symmetry("myMOF_notP1.cif", "myMOF_P1.cif")
```

## Forcefield

Forcefield input files are stored in `PorousMaterials.PATH_TO_DATA * "forcefields/"`.

## Molecule/Adsorbate

Molecule input files are stored in `PorousMaterials.PATH_TO_DATA * "molecules/"`.

## Energetics

## TODO
:cactus: document nearest image convention code (why the sign)
:cactus: in `vdw_energy`, ensure that the fractional coordinates are brought into the super cell through `mod()`. [But do not actually modify the molecule object]
:cactus: create `test_structure.cif` and `test_forcefield.csv` and test functions.
:cactus: (for Cory) generate test data from RASPA
:cactus: change `molecule.pos` to `molecule.x` for parallelism with `framework.xf` (`x` is notation for position in math)
:cactus: `translate_to!(molecule::Molecule)`, `perturb_coordinates!(molecule::Molecule)`
