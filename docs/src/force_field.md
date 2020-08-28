# Lennard-Jones Force Fields and Potential Energy
Lennard-Jones force field parameters are stored in comma-separated-value format in `PorousMaterials.PATH_TO_FORCEFIELDS`.

Interaction of an adsorbate with the crystal is modeled as pair-wise additive and with Lennard-Jones potentials of the form:

`V(r) = 4 * ϵ * [ x ^ 12 - x ^ 6 ]`, where `x = σ / r`

The Lennard-Jones force field input files, e.g. `UFF.csv` contain a list of pure (i.e. X-X, where X is an atom) sigmas and epsilons in units Angstrom (Å) and Kelvin (K), respectively. Note that, e.g., in the [UFF paper](https://doi.org/10.1021/ja00051a040), the Lennard-Jones potential is written in a different form; thus, parameters need to be converted to correspond to the functional form used in `PorousMaterials.jl`.

## Building Blocks of PorousMaterials: Lennard-Jones Force Fields

### Loading Force Field Files and Acessing Atributes

Reading in Lennard-Jones force field parameters is made easy with the [`LJForceField`](@ref) function. Let's load in the parameters from the Universal Force Field file (`UFF.csv`):

```julia
using PorousMaterials

# read in Lennard-Jones force field parameters from the Universal Force Field
ljforcefield = LJForceField("UFF", cutoffradius=14.0, mixing_rules="Lorentz-Berthelot")
```

PorousMaterials will then output information about the force field file you just loaded:

```julia
    Force field: UFF
    Number of atoms included: 106
    Cut-off radius (Å) = 14.0
```

This also prints all of the atoms included in the loaded forcefield with their given ϵ and σ. This was excluded because it would use too much space on this page. 

We can access attributes `LJForceField` such as `pure_σ`, `pure_ϵ`, and interaction values:

```julia
# access the Lennard-Jones epsilon & sigma for Xe 
ljforcefield.pure_ϵ[:Xe] # K
ljforcefield.pure_σ[:Xe] # Å

# access the Lennard-Jones epsilon & sigma for Xe-C interactions
ljforcefield.ϵ[:Xe][:C]  # K
ljforcefield.σ²[:Xe][:C] # Å (store σ² for faster computation)
```

### Checking Force Field Coverage

When running simulations, it is necessary to have the force field terms for all of the atoms. This can be checked using [`forcefield_coverage`](@ref):
```julia
using PorousMaterials

ljforcefield = LJForceField("UFF")

xtal = Crystal("SBMOF-1.cif")
strip_numbers_from_atom_labels!(xtal)

# check is the atoms in the crystal are covered
forcefield_coverage(xtal, ljforcefield)

molecule = Molecule("CO2")

# check if the atoms in the molecule are covered
forcefield_coverage(molecule, ljforcefield)
```
### Simulation Box and the Cutoff Radius
Find the replication factors needed to make a supercell big enough to fit a sphere with the specified cutoff radius.In PorousMaterials.jl, rather than replicating the atoms in the home unit cell to build the supercell that serves as a simulation box, we replicate the home unit cell to form the supercell (simulation box) in a for loop.The [`replication_factors`](@ref) function ensures enough replication factors such that the nearest image convention can be applied.

```julia
using PorousMaterials

ljforcefield = LJForceField("UFF")

xtal = Crystal("SBMOF-1.cif")

r_cutoff = 14.0 # Å
repfactors = replication_factors(xtal.box, r_cutoff) 
```

### Potential Energies: Van der Waals

What is the van der Waals potential energy of a Xe adsorbate inside SBMOF-1 at `[0.0, 1.0, 3.0]` Cartesian coordinates using the UFF as a molecular model?

```julia
using PorousMaterials

xtal = Crystal("SBMOF-1.cif")
strip_numbers_from_atom_labels!(xtal)

ljforcefield = LJForceField("UFF")

molecule = Molecule("Xe")

# convert to fractional
molecule = Frac(molecule, xtal.box) 

translate_to!(molecule, Cart([0.0, 1.0, 0.0]), xtal.box) # need box b/c we're talking Cartesian

energy = vdw_energy(xtal, molecule, ljforcefield) # K
```

#### Potential Energies: Electrostatics

What is the electrostatic potential energy of a CO$_2$ adsorbate inside CAXVII\_clean at `[0.0, 1.0, 0.0]` Cartesian coordinate?

```julia
using PorousMaterials

xtal = Crystal("CAXVII_clean.cif") # has charges
strip_numbers_from_atom_labels!(xtal)

# load molecule and convert it to fractional
molecule = Molecule("CO2")
molecule = Frac(molecule, xtal.box)

# let's give it a random orientation
translate_to!(molecule, Cart([0.0, 1.0, 0.0]), xtal.box)
random_rotation!(molecule, xtal.box) 

# this is for speed. pre-compute k-vectors and allocate memory
eparams = setup_Ewald_sum(xtal.box, 12.0)
eikr = Eikr(molecule, eparams)

# now compute the electrostatic potential energy of a molecule inside a crystal
energy = electrostatic_potential_energy(xtal, molecule, eparams, eikr)
```

# detailed docs
## Forcefields
```@docs
    LJForceField
    replication_factors
    forcefield_coverage
```

## Potential Energy
```@docs
    PotentialEnergy
    SystemPotentialEnergy
```

## Nearest Image Conventions
```@docs
    nearest_image!
```

## Electrostatics Energy
```@docs
    Eikr
    total
    electrostatic_potential
    electrostatic_potential_energy
    precompute_kvec_wts
    setup_Ewald_sum
    total_electrostatic_potential_energy
```

## Van der Waals Energy
```@docs
    lennard_jones
    vdw_energy
    vdw_energy_no_PBC
```
