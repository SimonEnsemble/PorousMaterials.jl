# PorousMaterials.jl

A software package in [Julia](https://julialang.org/) for classical molecular modeling of adsorption in porous crystals.

## Installation

(development version)

Download and install the [Julia programming language](https://julialang.org/).

Loading all functions of PorousMaterials.jl into Julia is as easy as `using PorousMaterials`.

For Julia to find the PorousMaterials.jl package, add the directory where the source code 
is found to the `LOAD_PATH` variable in Julia, for example:

```julia
push!(LOAD_PATH, homedir() * "/Dropbox/PorousMaterials.jl/src/")
```

If you put the above line in your `~/.juliarc` file, the command will run every time you run Julia.

(eventually)

`Pkg.add("PorousMaterials")`

## Quick use cases

If you change your working directory to `PorousMaterials.jl/test`, `PorousMaterials.jl` will find the appropriate crystal structure, molecule, and force field input files to run the following illuminating examples. This should get you started with learning `PorousMaterials.jl`.

#### crystals

```julia
using PorousMaterials

framework = read_crystal_structure_file("SBMOF-1.cif")

replicate_to_xyz(framework, "replicated_framework.xyz", repfactors=(2,3,1))

ρ = crystal_density(framework) # kg/m3

cf = chemical_formula(framework)
```

#### molecules

```julia
co2 = read_molecule_file("CO2")

translate_to!(co2, [1.0, 3.0, 4.0])
translate_by!(co2, [1.0, 1.0, 1.0])
rotate!(co2)
```

#### forcefields

```julia
ljforcefield = read_forcefield_file("UFF.csv")
ljforcefield.ϵ[:Xe][:C] # K
ljforcefield.σ²[:Xe][:C] # Å
```

#### computing potential energies of molecules in crystals

Van der Waals interactions:
```julia
# Unit cell replication factors for application of periodic boundary conditions
repfactors = replication_factors(framework.box, ljforcefield)
energy = vdw_energy(framework, co2, ljforcefield, repfactors) # K
```
Electrostatic potential energy:
```
nu1000 = read_crystal_structure_file("NU-1000_Greg.cif")

k_rep_factors = (11, 11, 9)
α = 0.25
sr_cutoff = 12.0
rep_factors = replication_factors(nu1000, sr_cutoff)
sim_box = replicate_box(nu1000.box, rep_factors)
kvectors = compute_kvectors(sim_box, k_rep_factors, α)

x = [4.2, 23.1, 28.4]
ϕ = electrostatic_potential_energy(nu1000, co2, sim_box, rep_factors, sr_cutoff, kvectors, α) # [K]
```

#### writing energy grids

#### calculating the Henry coefficient of a molecule in a crystal

#### grand-canonical Monte Carlo simulation of a molecule in a crystal
```julia
xe = read_molecule_file("Xe")
sbmof1 = read_crystal_structure_file("SBMOF-1.cif")
ljforcefield = read_forcefield_file("Dreiding.csv")
temperature = 298.0 # K
pressure = 1.0 * 100000 # 1 bar = 100,000 Pa
results = gcmc_simulation(sbmof1, temperature, pressure, xe, ljforcefield,
                          n_burn_cycles=100, n_sample_cycles=100,
                          verbose=true)
```

## From where does PorousMaterials.jl read crystal structure, molecule, and force field files?
All input files are stored in `PorousMaterials.PATH_TO_DATA`, which by default is 
`pwd() * "/data/"`. You can change the variable `PATH_TO_DATA` by editing `PorousMaterials.jl`.

### Crystals

Crystal structure files are stored in `PorousMaterials.PATH_TO_DATA * "crystals/"`. Crystals 
in .cif format must be in P1 symmetry. If your .cif is not in P1 symmetry, our function
`convert_cif_to_P1_symmetry()` calls the Atomic Simulation Environment (ASE) in Python to 
write the .cif in P1 symmetry.

`PorousMaterials.jl` includes a .cif and .cssr reader.

For example:

```
# crystal structure files in `PATH_TO_DATA/crystals`.
convert_cif_to_P1_symmetry("myMOF_notP1.cif", "myMOF_P1.cif")
```

### Forcefields

Forcefield input files are stored in `PorousMaterials.PATH_TO_DATA * "forcefields/"`.

#### Lennard-Jones forcefield

Interaction of an adsorbate with the framework is modeled as pair-wise additive and with Lennard-Jones potentials of the form:

`V(r) = 4 * epsilon * [ x ^ 12 - x ^ 6 ]`, where `x = sigma / r`

The Lennard Jones force field input files, e.g. `UFF.csv` contain a list of pure (i.e. X-X, where X is an atom) sigmas and epsilons in units Angstrom and Kelvin, respectively. Note that, in the UFF paper, the Lennard Jones potential is written in a different form and thus parameters need to be converted to correspond to the functional form used in `PorousMaterials.jl`.

### Molecules/Adsorbates

Molecule input files are stored in `PorousMaterials.PATH_TO_DATA * "molecules/"`. Each molecule possesses its own directory and contains two files: `point_charges.csv` and `lennard_jones_spheres.csv`. Both are .csv files of the point charges and Lennard Jones spheres, respectively, comprising the molecule. Only rigid molecules are currently supported. See `test/data` for examples.

### Energetics

# TODO
* random rotation matrix needs testing and is currently wrong.
* add hook so that tests run before committing
* henry coeff (Melanie)
* fix GCMC to allow complex molecules (Arthur)
* Ewald for molecules in GCMC (Arthur)
* make sure no 100000000000 - 10000000000.01 is going on with electrostatics adding to vdw...
* add guest-guest tests for CO2, other moleucles that are multiple beads
* get test data from Greg on ENTIRE molecule's energies

# Help wanted
* geometric based pore size calculations (largest free and included spheres), surface area, and porosity calculations

# Contribution guidelines

To avoid breaking PorousMaterials.jl upon pushing changes, please copy the `pre-commit` script to the `.git/hooks/` directory inside `PorousMaterials.jl`. Then `chmod +x .git/hooks/pre-commit` to make it executable. Now, when you `git commit`, `tests/runtests.jl` will automatically run to make sure the tests run before you commit a change that could break PorousMaterials.jl.

* keep with spacing and naming conventions used throughout the code. only lower case for variables, upper case for types etc.
* include many comments. include doc strings for your functions
* modularize the code as much as possible by breaking it into small functions
* before you implement a function, check if it already exists; we want to minimize the repeating of code. Less is more!
* no new functions without tests added to `tests/runtests.jl`
