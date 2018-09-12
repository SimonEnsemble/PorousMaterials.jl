## Building blocks

All of the commands below (and above) should run if you're in the `PorousMaterials.jl/test` directory so that `PorousMaterials.jl` can find the right input files. By default, if you `Pkg.add()`'d `PorousMaterials.jl`, the `test` directory is in `~/.julia/packages/UUID/PorousMaterials`.
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

#### Van der Waals

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
