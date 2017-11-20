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
in .cif format must be in P1 symmetry. If your .cif is not in P1 symmetry, use the Atomic 
Simulation Environment (ASE) in Python to convert to P1 symmetry as follows (Thanks Efrem):

```python
import ase
import ase.io
import ase.build

input_cif_filename='FAU_Si.cif'
output_cif_filename='FAU_Si_unwrapped.cif'

orig_cif = ase.io.read(input_cif_filename, format='cif')
cif_unwrapped = ase.build.make_supercell(orig_cif, [[1,0,0], [0,1,0], [0,0,1]])
ase.io.write(output_cif_filename, cif_unwrapped, format='cif')
```

## Forcefield

Forcefield input files are stored in `PorousMaterials.PATH_TO_DATA * "forcefields/"`.

## Molecule/Adsorbate

Molecule input files are stored in `PorousMaterials.PATH_TO_DATA * "molecules/"`.

## Energetics

## TODO
-- `UFF.csv` epsilon units should be in Kelvin. Also the functional form for LJ potential is different for UFF.
