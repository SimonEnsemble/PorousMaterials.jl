![PorousMaterials.jl](PMlogo.png)

A pure-[Julia](https://julialang.org/) package for classical molecular modeling of adsorption in porous crystals such as metal-organic frameworks (MOFs).

:hammer: Compute the potential energy of a molecule at particular position and orientation inside of a porous crystal

:hammer: Write a potential energy grid of a molecule inside a porous material to visualize binding sites

:hammer: Compute the Henry coefficient of a gas in a porous crystal

:hammer: Run grand-canonical Monte Carlo simulations of gas adsorption in a porous crystal

Designed for high-throughput computations to minimize input files and querying results from output files. User-friendly. Instructive error messages thrown when they should be. Well-documented (eventually). Easy to install (eventually).

*In development, please contribute, post issues :bug:, and improve!*

## Quick demos

### Henry coefficients

Compute the Henry coefficient of CO<sub>2</sub> in CAXVII\_clean (Fe<sub>2</sub>(dobdc)) at 298 K using the Dreiding force field:

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
results = henry_coefficient(framework, molecule, temperature, forcefield, insertions_per_volume=200)

# ... prints stuff
# results automatically saved to .jld load later in one line of code

# returns dictionary for easy querying
results["Qst (kJ/mol)"] # -21.0
results["henry coefficient [mol/(kg-Pa)]"] # 2.88e-05
```

The simulation is parallelized across a maximum of 5 cores.

### Grand-canonical Monte Carlo simulations

Simulate the adsorption of CO<sub>2</sub> in FIQCEN\_clean\_min\_charges (CuBTC) at 298 K at 1 bar using the Universal Force Field:

```julia
using PorousMaterials

# read in xtal structure file and populate a Framework data structure
framework = Framework("FIQCEN_clean_min_charges.cif")
# remove annoying numbers from atom labels
strip_numbers_from_atom_labels!(framework)

# read in Lennard-Jones force field parameters and populate a LJForceField data structure
forcefield = LJForceField("UFF.csv", cutoffradius=12.8)

# read in a molecule format file and populate a Molecule data structure
molecule = Molecule("CO2")

temperature = 298.0 # K
pressure = 1.0 # bar

# conduct grand-canonical Monte Carlo simulation
results, molecules = gcmc_simulation(framework, molecule, temperature, pressure, forcefield,
            n_burn_cycles=5000, n_sample_cycles=5000)

# ... prints stuff
# results automatically saved to .jld load later in one line of code

# returns dictionary for easy querying
results["⟨N⟩ (molecules/unit cell)"]
results["Q_st (K)"]
```

Or, compute the entire adsorption isotherm at once, parallelized across many cores:
```julia
pressures = [0.2, 0.6, 0.8, 1.0] # bar

# loop over all pressures and compute entire adsorption isotherm in parallel
results = adsorption_isotherm(framework, molecule, temperature, pressures, forcefield,
            n_burn_cycles=5000, n_sample_cycles=5000)
```

Or, compute the adsorption isotherm in a step-wise manner, loading the molecules from the previous simulation to save on burn cycles:
```julia
# loop over all pressures and run GCMC simulations in series. 
# load in the configurations of the molecules from the previous pressure.
results = stepwise_adsorption_isotherm(framework, molecule, temperature, pressures, forcefield,
            n_burn_cycles=1000, n_sample_cycles=5000)
```

### Potential Energy Grid

Superimpose a grid of points about the unit cell of SBMOF-1. Compute the potential energy of xenon at each point and store as a grid.

```julia
using PorousMaterials

framework = Framework("SBMOF-1.cif")
molecule = Molecule("Xe")
forcefield = LJForceField("UFF.csv")

grid = energy_grid(framework, molecule, forcefield, 
    n_pts=(50, 50, 50), units=:kJ_mol) # Grid data structure
```

Write to a .cube volume file to visualize the potential energy contours.
```julia
write_cube(grid, "CH4_in_SBMOF1.cube")
```

## Building blocks

All of the commands below (and above) should run if you're in the `PorousMaterials.jl/test` directory so that `PorousMaterials.jl` can find the right input files. By default, if you `Pkg.clone()`'d `PorousMaterials.jl`, the `test` directory is in `~/.julia/v0.6/PorousMaterials`.
Just fire up Julia and type in:

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
box.Ω # volume (1000.0 Å³)
box.f_to_c # fractional to Cartesian coordinate transformation matrix
box.c_to_f # Cartesian to fractional coordinate transformation matrix
box.reciprocal_lattice # rows are reciprocal lattice vectors
```

Replicate a box as follows:
```julia
box = replicate(box, (2, 2, 2)) # new box replicated 2 by 2 by 2
box.a # 20 Å
```

### Porous Crystals

```julia
using PorousMaterials

# read in xtal structure file
framework = Framework("SBMOF-1.cif")

# access unit cell box
framework.box

# access Lennard-Jones spheres and point charges comprising the crystal
framework.atoms
framework.charges

# remove annoying numbers on the atom labels
strip_numbers_from_atom_labels!(framework)

# compute crystal density
ρ = crystal_density(framework) # kg/m3

# compute the chemical formula
cf = chemical_formula(framework)

# assign charges according to atom type
charges = Dict(:Ca => 3.0, :O => 2.0, :C => -1.0, :S => 7.0, :H => -1.0)
charged_framework = assign_charges(framework, charges)

# replicate & visualize
framework = replicate(framework, (3, 3, 3))
write_to_xyz(framework, "SBMOF-1.xyz")
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
forcefield.σ²[:Xe][:C] # Å (store σ² for faster computation)
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

What is the electrostatic potential energy of a CO<sub>2</sub> adsorbate inside CAXVII\_clean at `[0.0, 1.0, 0.0]` Cartesian coordinate?

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

#### Equations of state

Calculate fugacity, density of methane at 298 K and 65 bar using the Peng-Robinson EOS:
```julia
gas = PengRobinsonGas(:CH4)
props = calculate_properties(gas, 298.0, 65.0) # dictionary of properties
props["fugacity coefficient"] # 0.8729
```

Pass `eos=:PengRobinson` to `gcmc_simulation` to automatically convert pressure to fugacity using the Peng-Robinson equation of state.

## Input files to describe crystals, molecules, and forcefields

All input files are stored in the path `PorousMaterials.PATH_TO_DATA` (type into Julia). By default, this path is set to be in the present working directory (type `pwd()` into Julia) in a folder `data/`. Go inside `PorousMaterials.jl/test/data` to see example input files for each case below.

#### Crystals

Place `.cif` and `.cssr` crystal structure files in `data/crystals`. `PorousMaterials.jl` currently takes crystals in P1 symmetry only.

#### Molecules/Adsorbates

Molecule input files are stored in `data/molecules`. Each molecule possesses its own directory and contains two files: `point_charges.csv` and `lennard_jones_spheres.csv`, comma-separated-value files describing the point charges and Lennard Jones spheres, respectively, comprising the molecule. Only rigid molecules are currently supported. Units of length are in Angstrom; units of charges are electrons.

#### Force field parameters

Lennard-Jones forcefield parameters are stored in comma-separated-value format in `data/forcefields/`.

Interaction of an adsorbate with the framework is modeled as pair-wise additive and with Lennard-Jones potentials of the form:

`V(r) = 4 * epsilon * [ x ^ 12 - x ^ 6 ]`, where `x = sigma / r`

The Lennard Jones force field input files, e.g. `UFF.csv` contain a list of pure (i.e. X-X, where X is an atom) sigmas and epsilons in units Angstrom and Kelvin, respectively. Note that, e.g., in the UFF paper, the Lennard Jones potential is written in a different form and thus parameters need to be converted to correspond to the functional form used in `PorousMaterials.jl`.

#### Atomic masses

Add fancy pseudo-atoms to `data/atomic_masses.csv`.

#### Peng-Robinson gas parameters

Critical temperatures and pressures and acentric factors are stored in `data/PengRobinsonGasProps.csv`.

## Installation

1. Download and install the [Julia programming language](https://julialang.org/), v0.6.4.

2. In Julia, type `Pkg.clone("https://github.com/SimonEnsemble/PorousMaterials.jl.git")` to clone this repository and install Julia package dependencies in `REQUIRE`.

3. In Julia, load all functions in `PorousMaterials.jl` into the namespace:

```julia
using PorousMaterials # that's it
```

Note: This package is in development. After stabilized and fully documented, installation will be as easy as `Pkg.add("PorousMaterials")`.

## Tests
Run the unit-ish tests in the script `tests/runtests.jl` manually or type `Pkg.test("PorousMaterials")` into Julia.

Direct tests for Henry coefficients and grand-canonical Monte Carlo simulations take much longer and are found in `tests/henry_test.jl` and `tests/gcmc_test.jl`.

## FAQ

**How do I type out the math symbols? e.g. `box.α`?**

Julia supports [unicode input](https://docs.julialang.org/en/release-0.4/manual/unicode-input/)! Type `box.\alpha`, then hit tab. Voilà. There is a vim extension for Julia [here](https://github.com/JuliaEditorSupport/julia-vim). 


**How do I run as a script in the command line?**

It is instructive to first run an example in the Julia REPL so you can print out and interact with attributes of your `forcefield`, `framework`, and `molecule` to ensure they are correct. If you want to then run the Julia code in the command line, simply put the commands in a text file with a `.jl` extension and run in terminal as `julia my_script.jl`. For parallelization in `adsorption_isotherm` and `henry_coefficient`, call e.g. 4 cores with `julia -p 4 my_script.jl`.

**Can I use `PorousMaterials.jl` in Jupyter Notebook/ Jupyter Lab?**

Yes! See [here](https://github.com/JuliaLang/IJulia.jl).

**How can I convert my `.cif` into P1 symmetry for `PorousMaterials.jl`?**

We hope someone will contribute this feature to `PorousMaterials.jl` eventually. For now, you can use [OpenBabel](http://openbabel.org/wiki/Main_Page):

```
obabel -icif non-P1.cif -ocif -O P1.cif --fillUC strict
```

## Help wanted and needed
* the speed of a GCMC or Henry simulation is determined primarily by how fast `PorousMaterials.jl` can compute the electrostatic and vdw potential energy. Some core functions that can speed up this are:
    * `nearest_image!`, `nearest_r` in `src/NearestImage.jl`
    * Ewald sums in `src/Electrostatics.jl`. (electrostatics are a huge bottleneck.)
    * `src/VdWEnergetics.jl`
    The scripts `test/vdw_timing.jl` and `test/ewald_timing.jl` time the functions for benchmarking.
* consolidate `eikar`, `eikbr`, `eikcr` somehow without slowing down the Ewald sum.
* more tests added to `tests/runtests.jl`, `tests/henry_tests.jl`, `tests/gcmc_tests.jl`
* code coverage badge
* how to hook up to Travis CI to automatically run tests upon a pull request?
* geometric-based pore size calculations (largest free and included spheres), surface area, and porosity calculations that take `Framework`'s as input
* handle .cif's without P1 symmetry. i.e. convert any .cif to P1 symmetry
* generate a docs website
* extend `gcmc_simulation` to handle mixtures
* better default rules for choosing Ewald sum parameters? alpha, kvectors required...
* Henry coefficient code prints off Ewald sum params 5 times if run with one core...
* set good defaults for `gcmc_simulation` probabilities (as now) but also allow user to change through default arguments to the function
* automatically adjust the translation step `δ` in `gcmc_simulation` during burn cycles to have 50% acceptance of translation moves (online gradient descent?)
* EQEq or other charge equilibration schemes for assinging charges, taking a `Framework` as input.

## Contribution guidelines

Please run `tests/runtests.jl` and assert that the tests run before you submit a pull request.
For substantial changes aside from performance optimizations/bug fixes, please check with us before moving forward. And it's best if you post an issue stating your intentions to implement a feature in case someone else already is.

* keep with spacing and naming conventions used throughout the code. only lower case for variables, upper case for types etc.
* always have type assertions in the function arguments
* include doc strings for your functions that are exposed to the user or comments for internal functions
* modularize the code as much as possible by breaking it into small functions
* before you implement a function, check if it already exists; we want to minimize the repeating of code. Less is more!
* ensure your new function has tests added to `tests/runtests.jl`
