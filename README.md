![PorousMaterials.jl](PMlogo.png)

A pure-[Julia](https://julialang.org/) package for classical molecular modeling of adsorption in porous crystals such as metal-organic frameworks (MOFs).

:hammer: Compute the potential energy of a molecule at particular position and orientation inside of a porous crystal

:hammer: Write a potential energy grid of a molecule inside a porous material to visualize binding sites

:hammer: Compute the Henry coefficient of a gas in a porous crystal

:hammer: Run grand-canonical Monte Carlo simulations of gas adsorption in a porous crystal

Designed for high-throughput computations to minimize input files and querying results from output files. User-friendly. Well-documented (eventually). Easy to install (eventually).

*In development, please contribute, post issues :bug:, and improve!*

## Quick demos

### Henry coefficients

Compute the Henry coefficient of CO2 in CAXVII\_clean at 298 K using the Dreiding force field:

```julia
using PorousMaterials

# read in xtal structure file and populate a Framework data structure
framework = Framework("CAXVII_clean.cif")                                               

# read in Lennard-Jones force field parameters and populate a LJForceField data structure
forcefield = LJForceField("Dreiding.csv", cutoffradius=12.5)                                  

# read in a molecule format file and populate a Molecule data structure
molecule = Molecule("CO2")                                                              

temperature = 298.0 # K

# conduct Widom insertions and compute Henry coefficient, heat of adsorption
result = henry_coefficient(framework, molecule, temperature, forcefield, nb_insertions_per_volume=200)

# ... prints stuff
# results automatically saved to .jld load later in one line of code

# returns dictionary for easy querying
results["Q_st (kJ/mol)"] # -21.0
result["henry coefficient [mol/(kg-Pa)]"] # 2.88e-05
```

### Grand-canonical Monte Carlo simulations

Simulate the adsorption isotherm of CO2 in ZIF-71 at 298 K from 0 to 1 bar using the Dreiding force field:

```julia
using PorousMaterials

# read in xtal structure file and populate a Framework data structure
framework = Framework("ZIF-71.cif")
# remove annoying numbers from atom labels
strip_numbers_from_atom_labels!(framework)

# read in Lennard-Jones force field parameters and populate a LJForceField data structure
forcefield = LJForceField("Dreiding.csv", cutoffradius=12.8)

# read in a molecule format file and populate a Molecule data structure
molecule = Molecule("CO2")

temperature = 298.0 # K
pressures = [0.2, 0.6, 0.8, 1.0] # bar

results = adsorption_isotherm(framework, temperature, pressures, molecule, forcefield,
            n_burn_cycles=5000, n_sample_cycles=5000)
```

## Building blocks

All of the commands below should run if you're in the `PorousMaterials.jl/test` directory. Just fire up Julia and type in:
```julia
using PorousMaterials
```

### Matter

In `PorousMaterials.jl`, crystals and molecules are composed of Lennard-Jones spheres and point charges.

To create a carbon atom at `[0.1, 0.2, 0.5]` fractional coordinates (in the context of some Bravais lattice):
```julia
ljs = LJSphere(:C, [0.1, 0.2, 0.5]) # constructor
ljs.species # :C
ljs.xf # [0.1, 0.2, 0.5]
```

To create a point charge of +1 at `[0.1, 0.2, 0.5]` fractional coordinates (in the context of some Bravais lattice):
```julia
ptc = PtCharge(1.0, [0.1, 0.2, 0.5])
ptc.q # 1.0
ptc.xf # [0.1, 0.2, 0.5]
```

### Bravais lattice

We later apply periodic boundary conditions to mimic a crystal of infinite extent. A `Box` describes a [Bravais lattice](https://en.wikipedia.org/wiki/Bravais_lattice).

To make a 10 by 10 by 10 Å Bravais lattice with right angles:
```julia
box = Box(10.0, 10.0, 10.0, π/2, π/2, π/2)

box.a, box.b, box.c # unit cell dimensions (10.0 Å)
box.α, box.β, box.γ # unit cell angles (1.57... radians)
box.Ω # volume (100 cubic Å)
box.f_to_c # fractional to Cartesian coordinate transformation matrix
box.c_to_f # Cartesian to fractional coordinate transformation matrix
box.reciprocal_lattice # rows are reciprocal lattice vectors
```

You can replicate a box as follows:
```julia
box = replicate(box, (2, 2, 2)) # new box replicated 2 by 2 by 2
box.a # 20 Å
```

### Porous Crystals

```julia
using PorousMaterials

# read in xtal structure file
framework = Framework("SBMOF-1.cif")

# access unit cell, atoms, and charges
framework.box
framework.atoms
framework.charges

# remove annoying numbers on the atom labels
strip_numbers_from_atom_labels!(framework)

# compute crystal density
ρ = crystal_density(framework) # kg/m3

# compute the chemical formula
cf = chemical_formula(framework)
```

### Lennard-Jones forcefields

```julia
# read in Lennard-Jones force field parameters from the Universal Force Field
forcefield = LJForceField("UFF.csv", cutoffradius=14.0, mixing_rules="Lorentz-Berthelot")

# access the Lennard-Jones epsilon & sigma for Xe
forcefield.pure_ϵ[:Xe] # K
forcefield.pure_σ[:Xe] # Å
# access the Lennard-Jones epsilon & sigma for Xe-C interactions
forcefield.ϵ[:Xe][:C] # K                                                                 
forcefield.σ²[:Xe][:C] # Å 
```

### Molecules 

```julia
molecule = Molecule("CO2") # fractional coords in terms of unit cube box

# access Lennard-Jones spheres & point charges that comprise molecule
molecule.atoms
molecule.charges

# translate to [1.0, 2.0, 3.0] fractional coordinates
translate_to!(molecule, [1.0, 2.0, 3.0])

# translate by [0.1, 0.0, 0.0] fractional coordinates
translate_by!(molecule, [0.1, 0.0, 0.0])

# conduct a uniform random rotation
rotate!(molecule, UnitCube()) # b/c now fractional coords defined in context of a unit cube
```

### Potential energies

First, set the fractional coordinates of the molecule in the context of some unit cell box.

```julia
# molecule in a framework
set_fractional_coords!(molecule, framework.box)

# molecule in a 10 by 10 by 10 cube
box = Box(10.0, 10.0, 10.0, π/2, π/2, π/2) # make a box
set_fractional_coords!(molecule, box)
```

#### van der Waals

What is the van der Waals potential energy of a Xe adsorbate inside SBMOF-1 at `[0.0, 1.0, 3.0]` Cartesian coordinates using the UFF as a molecular model?

```julia
using PorousMaterials

framework = Framework("SBMOF-1.cif")

forcefield = LJForceField("UFF.csv")

molecule = Molecule("Xe")
set_fractional_coords!(molecule, framework.box)

translate_to!(molecule, [0.0, 1.0, 0.0], framework.box) # need box b/c we're talking Cartesian

energy = vdw_energy(framework, molecule, forcefield) # K
```

#### Electrostatics

What is the electrostatic potential energy of a CO2 adsorbate inside CAXVII\_clean at `[0.0, 1.0, 0.0]` Cartesian coordinate?

```julia
using PorousMaterials

framework = Framework("CAXVII_clean.cif") # has charges

molecule = Molecule("CO2")
set_fractional_coords!(molecule, framework.box)

translate_to!(molecule, [0.0, 1.0, 0.0], framework.box) # need box b/c we're talking Cartesian

rotate!(molecule, framework.box) # let's give it a random orientation

# this is for speed. pre-compute k-vectors and allocate memory
eparams, kvectors, eikar, eikbr, eikcr = setup_Ewald_sum(12.0, framework.box)

energy = electrostatic_potential_energy(framework, molecule, eparams, kvectors, eikar, eikbr, eikcr)
```

## Input files to describe crystals, molecules, and forcefields

All data is stored in the variable `PorousMaterials.PATH_TO_DATA * "data/"` (type into Julia). By default, the input files are read from your present working directory (type `pwd()` into Julia) in a folder `data/`. Go inside `PorousMaterials.jl/test/data` to see example input files.

#### Crystals

Put your .cif and .cssr crystal structure files in `data/crystals`. `PorousMaterials.jl` only takes crystals in P1 symmetry.

#### Molecules/Adsorbates

Molecule input files are stored in `data/molecules`. Each molecule possesses its own directory and contains two files: `point_charges.csv` and `lennard_jones_spheres.csv`. Both are .csv files of the point charges and Lennard Jones spheres, respectively, comprising the molecule. Only rigid molecules are currently supported. Units of length are in Angstrom here. Charges are in units of electrons.

#### Force field parameters

Lennard-Jones forcefield parameters are stored in .csv format in `data/forcefields/`.

Interaction of an adsorbate with the framework is modeled as pair-wise additive and with Lennard-Jones potentials of the form:

`V(r) = 4 * epsilon * [ x ^ 12 - x ^ 6 ]`, where `x = sigma / r`

The Lennard Jones force field input files, e.g. `UFF.csv` contain a list of pure (i.e. X-X, where X is an atom) sigmas and epsilons in units Angstrom and Kelvin, respectively. Note that, in the UFF paper, the Lennard Jones potential is written in a different form and thus parameters need to be converted to correspond to the functional form used in `PorousMaterials.jl`.

## Installation

1. Download and install the [Julia programming language](https://julialang.org/), v0.6.4.

2. Clone `PorousMaterials.jl` from this repo. Remember the directory where you save it.

3. The package dependencies of `PorousMaterials.jl` are found in the `REQUIRE` file in the `PorousMaterials.jl` directory. Install these package requirements in Julia via e.g. `Pkg.install("CSV")`.

4. For Julia to find the `PorousMaterials.jl` source code, add the directory where you place the cloned Github repo to the `LOAD_PATH` variable in Julia. For example:
```julia
push!(LOAD_PATH, homedir() * "/PorousMaterials.jl/src/")
```
If you put the above line in your `~/.juliarc` file, the command will run silently every time you run Julia!

5. Finally, in Julia, load all functions in `PorousMaterials.jl` into the namespace:

```julia
using PorousMaterials # that's it
```

Note: This package is in development. After stabilized and fully documented, installation will be as easy `Pkg.add("PorousMaterials")`.

## Tests
Please run the test script in `tests/runtests.jl`.

Direct tests for Henry coefficients and grand-canonical Monte Carlo simulations take much longer and are found in `tests/henry_test.jl` and `tests/gcmc_test.jl`.

## Help wanted
* geometric-based pore size calculations (largest free and included spheres), surface area, and porosity calculations that take `Framework`'s as input
* the speed of a GCMC or Henry simulation is determined primarily by how fast `PorousMaterials.jl` can compute the energy. Some core functions that can speed up this are:
    * `nearest_image!`, `nearest_r` in `src/NearestImage.jl`
    * electrostatics take the longest... Ewald sums in `src/Electrostatics.jl`.
    * consolidate `eikar`, `eikbr`, `eikcr` somehow without slowing down the Ewald sum
    The scripts `test/vdw_timing.jl` and `test/ewald_timing.jl` time the functions so you can check if you are speeding it up.
* extend GCMC to mixtures

# Contribution guidelines

To avoid breaking PorousMaterials.jl upon pushing changes, please copy the `pre-commit` script to the `.git/hooks/` directory inside `PorousMaterials.jl`. Then `chmod +x .git/hooks/pre-commit` to make it executable. Now, when you `git commit`, `tests/runtests.jl` will automatically run to make sure the tests run before you commit a change that could break PorousMaterials.jl.

* keep with spacing and naming conventions used throughout the code. only lower case for variables, upper case for types etc.
* include many comments. include doc strings for your functions
* modularize the code as much as possible by breaking it into small functions
* before you implement a function, check if it already exists; we want to minimize the repeating of code. Less is more!
* no new functions without tests added to `tests/runtests.jl`
