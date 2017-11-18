# PorousMaterials.jl

A software package in [Julia](https://julialang.org/) for classical molecular modeling in porous crystals.

Example use:

```julia
using PorousMaterials

# load abstract data structure containg crystal structure info
framework = read_crystal_structure_file("cifFiles/JAVTAC_clean.cif")

# replicate framework to .xyz for visualization
replicate_to_xyz(framework, "replicated_framework.xyz", repfactors=(2,3,1))

# load abstract data structure Lennard Jones force field (UFF)
forcefield = read_forcefield_file("UFF.csv") 

```

## Crystal

## Forcefield

## Molecule/Adsorbate

## Energetics
