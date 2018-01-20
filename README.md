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

### Lennard-Jones forcefield

Interaction of an adsorbate with the framework is modeled as pair-wise additive and with Lennard-Jones potentials of the form:

`V(r) = 4 * epsilon * [ x ^ 12 - x ^ 6 ], where x = sigma / r`

The Lennard Jones force field input files, e.g. `UFF.csv` contain a list of pure (i.e. X-X, where X is an atom) sigmas and epsilons in units Angstrom and Kelvin, respectively. Note that, in the UFF paper, the Lennard Jones potential is written in a different form and thus parameters need to be converted to correspond to the functional form used in `PorousMaterials.jl`.

## Molecule/Adsorbate

Molecule input files are stored in `PorousMaterials.PATH_TO_DATA * "molecules/"`.

## Energetics

## TODO
* ~~document nearest image convention code (why the sign)~~
* ~~in `vdw_energy`, ensure that the fractional coordinates are brought into the super cell through `mod()`. [But do not actually modify the molecule object]~~
* ~~create `test_structure.cif` and `test_forcefield.csv` and test functions.~~
* ~~(for Cory) generate test data from RASPA~~
* ~~change `molecule.pos` to `molecule.x` for parallelism with `framework.xf` (`x` is notation for position in math)~~
* `translate_to!(molecule::Molecule)`, `perturb_coordinates!(molecule::Molecule)`
* ~~`print(framework::Framework)`, `print(forcefield::LJForcefield)`, `print(molecule::Molecule)`~~
* Read in charges, assign to `Framework` in `read_crystal_structure`
