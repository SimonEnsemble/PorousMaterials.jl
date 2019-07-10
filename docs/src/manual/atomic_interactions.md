## Loading in Forcefield Files

Lennard-Jones forcefield parameters are stored in comma-separated-value format in `data/forcefields/`.

Interaction of an adsorbate with the framework is modeled as pair-wise additive and with Lennard-Jones potentials of the form:

`V(r) = 4 * epsilon * [ x ^ 12 - x ^ 6 ]`, where `x = sigma / r`

The Lennard Jones force field input files, e.g. `UFF.csv` contain a list of pure (i.e. X-X, where X is an atom) sigmas and epsilons in units Angstrom and Kelvin, respectively. Note that, e.g., in the UFF paper, the Lennard Jones potential is written in a different form and thus parameters need to be converted to correspond to the functional form used in `PorousMaterials.jl`.

```julia
using PorousMaterials

ljff = LJForceField("UFF.csv")
```

PorousMaterials will the output information about the forcefield file you just loaded:

```julia
    Force field: UFF.csv
    Number of atoms included: 106
    Cut-off radius (Å) = 14.0    
```

This also prints all of the atoms included in the loaded forcefield with their given ϵ and σ. This was excluded because it would use too much space on this page.

## Building Blocks of PorousMaterials: Lennard-Jones forcefields

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

## Building Blocks of PorousMaterials: Potential energies

First, set the fractional coordinates of the molecule in the context of some unit cell box.

```julia
# molecule in a framework
set_fractional_coords!(molecule, framework.box)

# molecule in a 10 by 10 by 10 cube
box = Box(10.0, 10.0, 10.0, π/2, π/2, π/2) # make a box
set_fractional_coords!(molecule, box)
```

#### Potential Energies: Van der Waals

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

#### Potential Energies: Electrostatics

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

#### Potential Energies: Equations of state

Calculate fugacity of methane at 298 K and 65 bar using the Peng-Robinson EOS:
```julia
fluid = PengRobinsonFluid(:CH4)
props = calculate_properties(fluid, 298.0, 65.0) # dictionary of PREOS properties
props["fugacity coefficient"] # 0.8729
```

Calculate compressibility factor of hydrogen at 300 K and 50 bar using van der Waals EOS:
```julia
fluid = VdWFluid(:H2)
props = calculate_properties(fluid, 300.0, 50.0) # dictionary of VdW properties
props["compressibility factor"] # 1.03511
```

Pass `eos=:PengRobinson` to `gcmc_simulation` to automatically convert pressure to fugacity using the Peng-Robinson equation of state.

## Forcefields
```@docs
    LJForceField
    replication_factors
    check_forcefield_coverage
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
