
<a id='Forcefields-1'></a>

## Forcefields

<a id='PorousMaterials.LJForceField' href='#PorousMaterials.LJForceField'>#</a>
**`PorousMaterials.LJForceField`** &mdash; *Type*.



Data structure for a Lennard Jones forcefield.

**Attributes**

  * `name::String`: name of forcefield; correponds to filename
  * `pure_σ::Dict{Symbol, Float64}`: Dictionary that returns Lennard-Jones σ of an X-X interaction, where X is an atom. (units: Angstrom)
  * `pure_ϵ::Dict{Symbol, Float64}`: Dictionary that returns Lennard-Jones ϵ of an X-X interaction, where X is an atom. (units: K)
  * `σ²::Dict{Symbol, Dict{Symbol, Float64}}`: Lennard Jones σ² (units: Angstrom²) for cross-interactions. Example use is `sigmas_squared[:He][:C]`
  * `ϵ::Dict{Symbol, Dict{Symbol, Float64}}`: Lennard Jones ϵ (units: K) for cross-interactions. Example use is `epsilons[:He][:C]`
  * `cutoffradius_squared::Float64`: The square of the cut-off radius beyond which we define the potential energy to be zero (units: Angstrom²). We store σ² to speed up computations, which involve σ², not σ.


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/845505808b4e0fb9638d52f835a7b7cb0cde5b8f/src/Forcefield.jl#L1-L11' class='documenter-source'>source</a><br>

<a id='PorousMaterials.replication_factors' href='#PorousMaterials.replication_factors'>#</a>
**`PorousMaterials.replication_factors`** &mdash; *Function*.



```
repfactors = replication_factors(unitcell, cutoffradius)
```

Find the replication factors needed to make a supercell big enough to fit a sphere with the specified cutoff radius. In PorousMaterials.jl, rather than replicating the atoms in the home unit cell to build the supercell that serves as a simulation box, we replicate the home unit cell to form the supercell (simulation box) in a for loop. This function ensures enough replication factors such that the nearest image convention can be applied.

A non-replicated supercell has 1 as the replication factor in each dimension (`repfactors = (1, 1, 1)`).

**Arguments**

  * `unitcell::Box`: The unit cell of the framework
  * `cutoff_radius::Float64`: Cutoff radius beyond which we define the potential energy to be zero (units: Angstrom)

**Returns**

  * `repfactors::Tuple{Int, Int, Int}`: The replication factors in the a, b, c directions


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/845505808b4e0fb9638d52f835a7b7cb0cde5b8f/src/Forcefield.jl#L100-L116' class='documenter-source'>source</a><br>

<a id='PorousMaterials.check_forcefield_coverage' href='#PorousMaterials.check_forcefield_coverage'>#</a>
**`PorousMaterials.check_forcefield_coverage`** &mdash; *Function*.



```
check_forcefield_coverage(framework, ljforcefield)
check_forcefield_coverage(molecule, ljforcefield)
```

Check that the force field contains parameters for every atom present in a framework or molecule. Will print out which atoms are missing.

**Arguments**

  * `framework::Framework`: The framework containing the crystal structure information
  * `molecule::Molecule`: A molecule object
  * `ljforcefield::LJForceField`: A Lennard Jones forcefield object containing information on atom interactions

**Returns**

  * `all_covered::Bool`: Returns true if all atoms in the `framework` are also included in `ljforcefield`. False otherwise


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/845505808b4e0fb9638d52f835a7b7cb0cde5b8f/src/Forcefield.jl#L161-L175' class='documenter-source'>source</a><br>


<a id='Potential-Energy-1'></a>

## Potential Energy

<a id='PorousMaterials.PotentialEnergy' href='#PorousMaterials.PotentialEnergy'>#</a>
**`PorousMaterials.PotentialEnergy`** &mdash; *Type*.



```
pe = PotentialEnergy()
```

Data structure to store potential energy, partitioned into van der Waals (`energy.vdw`) and electrostatic (`energy.coulomb`) interactions, both `Float64`.

This returns a PotentialEnergy data type where the vdw and coulomb attributes are set to 0.0

**Returns**

  * `pe::PotentialEnergy`: A structure containing van der Waals and electrostatic energies, initialized at 0.0

**Attributes**

  * `vdw::Float64`: The potential energy contributions from Van der Waals interactions
  * `coulomb::Float64`: The potential energy contributions from electrostatic interactions


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/845505808b4e0fb9638d52f835a7b7cb0cde5b8f/src/Energetics_Util.jl#L3-L20' class='documenter-source'>source</a><br>

<a id='PorousMaterials.SystemPotentialEnergy' href='#PorousMaterials.SystemPotentialEnergy'>#</a>
**`PorousMaterials.SystemPotentialEnergy`** &mdash; *Type*.



```
system_potential_energy = SystemPotentialEnergy()
```

Data structure to facilitate storing/partitioning potential energy of a system. It stores the potential energy from guest-host and guest-guest interactions separately.

This initializes guest*host and guest*guest with PotentialEnergy(), so when it is created the total energy recorded is 0.0

**Returns**

  * `system_potential_energy::SystemPotentialEnergy`: A structure containing the potential energy of the system,   broken down into guest-guest and guest-host interactions

**Attributes**

  * `guest_host::PotentialEnergy`: The total potential energy from all guest-host   interactions in the system
  * `guest_guest::PotentialEnergy`: The total potential energy from all guest-guest   interactions in the system


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/845505808b4e0fb9638d52f835a7b7cb0cde5b8f/src/Energetics_Util.jl#L51-L69' class='documenter-source'>source</a><br>


<a id='Nearest-Image-Conventions-1'></a>

## Nearest Image Conventions

<a id='PorousMaterials.nearest_image!' href='#PorousMaterials.nearest_image!'>#</a>
**`PorousMaterials.nearest_image!`** &mdash; *Function*.



```
nearest_image!(dxf)
```

Applies the nearest image convention on a vector `dxf` between two atoms in fractional space; modifies `dxf` for nearest image convention. Fractional coordinates here fall in [0, 1] so that the box is [0, 1]^3 in fractional space.

Warning: this assumes the two molecules are in the box described by fractional coords [0, 1]³.

**Arguments**

  * `dxf::Array{Float64}`: A vector between two atoms in fractional space


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/845505808b4e0fb9638d52f835a7b7cb0cde5b8f/src/NearestImage.jl#L1-L12' class='documenter-source'>source</a><br>


<a id='Electrostatics-Energy-1'></a>

## Electrostatics Energy

<a id='PorousMaterials.Eikr' href='#PorousMaterials.Eikr'>#</a>
**`PorousMaterials.Eikr`** &mdash; *Type*.



```
eikr = Eikr(eikar, eikbr, eikcr)
```

mutable struct for holding the eikr vectors

**Attributes**

  * `eikar::OffsetArray{Complex{Float64}}`: array for storing e^{i * ka ⋅ r}; has indices   0:kreps[1] and corresponds to recip. vectors in a-direction
  * `eikbr::OffsetArray{Complex{Float64}}`: array for storing e^{i * kb ⋅ r}; has indices   -kreps[2]:kreps[2] and corresponds to recip. vectors in b-direction
  * `eikcr::OffsetArray{Complex{Float64}}`: array for storing e^{i * kc ⋅ r}; has indices   -kreps[2]:kreps[1] and corresponds to recip. vectors in c-direction


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/845505808b4e0fb9638d52f835a7b7cb0cde5b8f/src/ElectrostaticEnergetics.jl#L77-L89' class='documenter-source'>source</a><br>

<a id='PorousMaterials.electrostatic_potential_energy' href='#PorousMaterials.electrostatic_potential_energy'>#</a>
**`PorousMaterials.electrostatic_potential_energy`** &mdash; *Function*.



```
ϕ = electrostatic_potential_energy(framework, molecule, eparams, eikr)
```

Compute the electrostatic potential energy of a molecule inside a framework.

The electrostatic potential is created by the point charges assigned to the framework atoms in `framework.charges`. Periodic boundary conditions are applied through the Ewald summation. The spurious self-interaction term is neglected here because we are looking at *differences* in energy in a Monte Carlo simulation.

Warning: it is assumed that the framework is replicated enough such that the nearest image convention can be applied for the short-range cutoff radius supplied in `eparams.sr_cutoff_r`.

**Arguments**

  * `framework::Framework`: Crystal structure (see `framework.charges` for charges)
  * `molecule::Molecule`: The molecule being compared to the atoms in the framework.
  * `eparams::EwaldParams`: data structure containing Ewald summation settings
  * `eikr::Eikr`: Stores the eikar, eikbr, and eikcr OffsetArrays used in this calculation.

**Returns**

  * `pot::EwaldSum`: Electrostatic potential between `framework` and `molecule` (units: K)


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/845505808b4e0fb9638d52f835a7b7cb0cde5b8f/src/ElectrostaticEnergetics.jl#L307-L330' class='documenter-source'>source</a><br>


```
ϕ = electrostatic_potential_energy(molecules, eparams, box, eikr)
```

Compute the electrostatic potential energy of a system comprised of an array of `Molecule`s.

The EWald summation is used here in a double for loop; do not use this function for Monte Carlo simulations because it is computationally expensive.

Returns an `EwaldSum` type containing short-range and long-range contributions to the Ewald sum as well as the spurious self-interaction and intramolecular interactions. Access via (ϕ.sr, ϕ.lr, ϕ.self, ϕ.intra).

Units of energy: Kelvin

**Arguments**

  * `molecules::Array{Molecules, 1}`: array of molecules comprising the system.
  * `eparams::EwaldParams`: data structure containing Ewald summation settings
  * `box::Box`: the box the energy is being computed in
  * `eikr::Eikr`: Stores the eikar, eikbr, and eikcr OffsetArrays used in this calculation.

**Returns**

  * `ϕ::GGEwaldSum`: The total electrostatic potential energy


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/845505808b4e0fb9638d52f835a7b7cb0cde5b8f/src/ElectrostaticEnergetics.jl#L411-L433' class='documenter-source'>source</a><br>

<a id='PorousMaterials.precompute_kvec_wts' href='#PorousMaterials.precompute_kvec_wts'>#</a>
**`PorousMaterials.precompute_kvec_wts`** &mdash; *Function*.



```
kvectors = precompute_kvec_wts(kreps, box, α, max_mag_k_sqrd=Inf)
```

For speed, pre-compute the weights for each reciprocal lattice vector for the Ewald sum in Fourier space. This function takes advantage of the symmetry:     cos(-k⋅(x-xᵢ)) + cos(k⋅(x-xᵢ)) = 2 cos(k⋅(x-xᵢ))

If `max_mag_k_sqrd` is passed, k-vectors with a magnitude greater than `max_mag_k_sqrd` are not included.

**Arguments**

  * `kreps::Tuple{Int, Int, Int}`: number of k-vector replications required in a, b, c
  * `box::Box`: the simulation box containing the reciprocal lattice.
  * `α::Float64`: Ewald sum convergence parameter (units: inverse Å)
  * `max_mag_k_sqrd::Float64`: cutoff for |k|² in Fourier sum; if passed, do not include

k-vectors with magnitude squared greater than this.

**Returns**

  * `kvectors::Array{Kvector, 1}`: array of k-vectors to include in the Fourier sum and their

corresponding weights indicating the contribution to the Fourier sum.


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/845505808b4e0fb9638d52f835a7b7cb0cde5b8f/src/ElectrostaticEnergetics.jl#L177-L197' class='documenter-source'>source</a><br>

<a id='PorousMaterials.setup_Ewald_sum' href='#PorousMaterials.setup_Ewald_sum'>#</a>
**`PorousMaterials.setup_Ewald_sum`** &mdash; *Function*.



```
eparams = setup_Ewald_sum(box, sr_cutoff_r; ϵ=1e-6, verbose=false)
```

Given the short-range cutoff radius and simulation box, automatically compute Ewald convergence parameter and number of k-vector replications in Fourier space required for a given precision. Constructs and returns Ewald parameters data type with this information.

Also, pre-compute weights on k-vector contributions to Ewald sum in Fourier space.

Also, allocate OffsetArrays for storing e^{i * k ⋅ r} where r = x - xⱼ and k is a reciprocal lattice vector.

**Arguments**

  * `box::Box`: the simulation box containing the reciprocal lattice.
  * `sr_cutoff_r::Float64`: cutoff-radius (units: Å) for short-range contributions to Ewald
  * `ϵ::Float64`: desired level of precision. Typical value is 1e-6, but this does not
  * `verbose::Bool`: If `true` will print results

**Returns**

  * `eparams::EwaldParams`: data structure containing Ewald summation settings

corresponding weights indicating the contribution to the Fourier sum.


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/845505808b4e0fb9638d52f835a7b7cb0cde5b8f/src/ElectrostaticEnergetics.jl#L234-L255' class='documenter-source'>source</a><br>

<a id='PorousMaterials.total_electrostatic_potential_energy' href='#PorousMaterials.total_electrostatic_potential_energy'>#</a>
**`PorousMaterials.total_electrostatic_potential_energy`** &mdash; *Function*.



```
total_ϕ = total_electrostatic_potential_energy(molecules, eparams, box, eikr)
```

Calculates the total electrostatic potential energy of an array of `Molecule`s using a Grand Canonical Monte Carlo (GCMC) algorithm. #TODO add to this

**Arguments**

  * `molecules::Array{Molecule, 1}`: The molecules comprising the system.
  * `eparams::EwaldParams`: data structure containing Ewald summation settings
  * `box::Box`: The box the energy is being computed in.
  * `eikr::Eikr`: Stores the eikar, eikbr, and eikcr OffsetArrays used in this calculation.

**Returns**

  * `ϕ::GGEwaldSum`: The total electrostatic potential energy


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/845505808b4e0fb9638d52f835a7b7cb0cde5b8f/src/ElectrostaticEnergetics.jl#L516-L531' class='documenter-source'>source</a><br>


```
total_ϕ = total_electrostatic_potential_energy(framework, molecules, eparams, eikr)
```

Explanation of total*electrostatic*potential_energy that uses framework

**Arguments**

  * `framework::Framework`: Crystal structure (see `framework.charges` for charges)
  * `molecules::Array{Molecule, 1}`: The molecules comprising the system.
  * `eparams::EwaldParams`: data structure containing Ewald summation settings
  * `eikr::Eikr`: Stores the eikar, eikbr, and eikcr OffsetArrays used in this calculation.


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/845505808b4e0fb9638d52f835a7b7cb0cde5b8f/src/ElectrostaticEnergetics.jl#L546-L557' class='documenter-source'>source</a><br>


<a id='Van-der-Waals-Energy-1'></a>

## Van der Waals Energy

<a id='PorousMaterials.lennard_jones' href='#PorousMaterials.lennard_jones'>#</a>
**`PorousMaterials.lennard_jones`** &mdash; *Function*.



```
energy = lennard_jones(r², σ², ϵ)  (units: Kelvin)
```

Calculate the lennard jones potential energy given the square of the radius r between two lennard-jones spheres. σ and ϵ are specific to interaction between two elements. Return the potential energy in units Kelvin (well, whatever the units of ϵ are).

**Arguments**

  * `r²::Float64`: distance between two (pseudo)atoms in question squared (Angstrom²)
  * `σ²::Float64`: sigma parameter in Lennard Jones potential squared (units: Angstrom²)
  * `ϵ::Float64`: epsilon parameter in Lennard Jones potential (units: Kelvin)

**Returns**

  * `energy::Float64`: Lennard Jones potential energy


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/845505808b4e0fb9638d52f835a7b7cb0cde5b8f/src/VdWEnergetics.jl#L4-L18' class='documenter-source'>source</a><br>

<a id='PorousMaterials.vdw_energy' href='#PorousMaterials.vdw_energy'>#</a>
**`PorousMaterials.vdw_energy`** &mdash; *Function*.



```
energy = vdw_energy(framework, molecule, ljforcefield)
```

Calculates the van der Waals interaction energy between a molecule and a framework. Applies the nearest image convention to find the closest replicate of a specific atom.

WARNING: it is assumed that the framework is replicated sufficiently such that the nearest image convention can be applied. See [`replicate`](../boxes_crystals_grids/docs.md#PorousMaterials.replicate).

**Arguments**

  * `framework::Framework`: Crystal structure
  * `molecule::Molecule`: adsorbate (includes position/orientation/atoms)
  * `ljforcefield::LJForceField`: Lennard Jones force field

**Returns**

  * `energy::Float64`: Van der Waals interaction energy


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/845505808b4e0fb9638d52f835a7b7cb0cde5b8f/src/VdWEnergetics.jl#L49-L65' class='documenter-source'>source</a><br>


```
gg_energy = vdw_energy(molecule_id, molecules, ljforcefield, simulation_box)
```

Calculates van der Waals interaction energy of a single adsorbate `molecules[molecule_id]` with all of the other molecules in the system. Periodic boundary conditions are applied, using the nearest image convention.

**Arguments**

  * `molecule_id::Int`: Molecule ID used to determine which molecule in `molecules` we wish to calculate the guest-guest interactions
  * `molecules::Array{Molecule, 1}`: An array of Molecule data structures
  * `ljforcefield::LJForceField`: A Lennard Jones forcefield data structure describing the interactions between different atoms
  * `simulation_box::Box`: The simulation box for the computation.

**Returns**

  * `gg_energy::Float64`: The guest-guest interaction energy of `molecules[molecule_id]` with the other molecules in `molecules`


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/845505808b4e0fb9638d52f835a7b7cb0cde5b8f/src/VdWEnergetics.jl#L70-L85' class='documenter-source'>source</a><br>

<a id='PorousMaterials.vdw_energy_no_PBC' href='#PorousMaterials.vdw_energy_no_PBC'>#</a>
**`PorousMaterials.vdw_energy_no_PBC`** &mdash; *Function*.



Assumes unit cell box is a unit cube and no periodic boundary conditions are applied.


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/845505808b4e0fb9638d52f835a7b7cb0cde5b8f/src/VdWEnergetics.jl#L134-L137' class='documenter-source'>source</a><br>

