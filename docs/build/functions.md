
#Functions This page contains all of the functions exported by PorousMaterials. They are sorted by the .jl files they are found in.


##Box.jl

<a id='PorousMaterials.Box' href='#PorousMaterials.Box'>#</a>
**`PorousMaterials.Box`** &mdash; *Type*.



```
box = Box(a, b, c, α, β, γ, volume, f_to_c, c_to_f, reciprocal_lattice)
box = Box(a, b, c, α, β, γ)
box = Box(f_to_c)
```

Data structure to describe a unit cell box (Bravais lattice) and convert between fractional and Cartesian coordinates.

**Attributes**

  * `a,b,c::Float64`: unit cell dimensions (units: Angstroms)
  * `α,β,γ::Float64`: unit cell angles (units: radians)
  * `Ω::Float64`: volume of the unit cell (units: cubic Angtroms)
  * `f_to_c::Array{Float64,2}`: the 3x3 transformation matrix used to map fractional

coordinates to cartesian coordinates. The columns of this matrix define the unit cell axes. Columns are the vectors defining the unit cell box. units: Angstrom

  * `c_to_f::Array{Float64,2}`: the 3x3 transformation matrix used to map Cartesian

coordinates to fractional coordinates. units: inverse Angstrom

  * `reciprocal_lattice::Array{Float64, 2}`: the *rows* are the reciprocal lattice vectors.

This choice was made (instead of columns) for speed of Ewald Sums.


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/Box.jl#L1-L20' class='documenter-source'>source</a><br>

<a id='PorousMaterials.replicate' href='#PorousMaterials.replicate'>#</a>
**`PorousMaterials.replicate`** &mdash; *Function*.



```
new_box = replicate(original_box, repfactors)
```

Replicates a `Box` in positive directions to construct a new `Box` representing a supercell. The `original_box` is replicated according to the factors in `repfactors`. Note `replicate(original_box, repfactors=(1, 1, 1))` returns same `Box`. The new fractional coordinates as described by `f_to_c` and `c_to_f` still ∈ [0, 1].

**Arguments**

  * `original_box::Box`: The box that you want to replicate
  * `repfactors::Tuple{Int, Int, Int}`: The factor you want to replicate the box by

**Returns**

  * `box::Box`: Fully formed Box object


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/Box.jl#L117-L131' class='documenter-source'>source</a><br>


```
replicated_frame = replicate(framework, repfactors)
```

Replicates the atoms and charges in a `Framework` in positive directions to construct a new `Framework`. Note `replicate(framework, (1, 1, 1))` returns the same `Framework`.

**Arguments**

  * `framework::Framework`: The framework to replicate
  * `repfactors::Tuple{Int, Int, Int}`: The factors by which to replicate the crystal structure in each direction.

**Returns**

  * `replicated_frame::Framework`: Replicated framework


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/Crystal.jl#L221-L233' class='documenter-source'>source</a><br>

<a id='PorousMaterials.UnitCube' href='#PorousMaterials.UnitCube'>#</a>
**`PorousMaterials.UnitCube`** &mdash; *Function*.



```
unit_cube = UnitCube()
```

This function generates a unit cube, each side is 1.0 Angstrom long, and all the corners are right angles.


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/Box.jl#L108-L113' class='documenter-source'>source</a><br>

<a id='PorousMaterials.write_vtk' href='#PorousMaterials.write_vtk'>#</a>
**`PorousMaterials.write_vtk`** &mdash; *Function*.



```
write_vtk(box, "box.vtk"; verbose=true)
write_vtk(framework)
```

Write a `Box` to a .vtk file for visualizing e.g. the unit cell boundary of a crystal. If a `Framework` is passed, the `Box` of that framework is written to a file that is the same as the crystal structure filename but with a .vtk extension.

Appends ".vtk" extension to `filename` automatically if not passed.

**Arguments**

  * `box::Box`: a Bravais lattice
  * `filename::AbstractString`: filename of the .vtk file output (absolute path)


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/Box.jl#L137-L150' class='documenter-source'>source</a><br>


##Crystal.jl

<a id='PorousMaterials.Framework' href='#PorousMaterials.Framework'>#</a>
**`PorousMaterials.Framework`** &mdash; *Type*.



```
framework = Framework(filename, check_charge_neutrality=true,
                      net_charge_tol=0.001, check_atom_and_charge_overlap=true,
                      remove_overlap=false)
framework = Framework(name, box, atoms, charges)
```

Read a crystal structure file (.cif or .cssr) and populate a `Framework` data structure, or construct a `Framework` data structure directly.

**Arguments**

  * `filename::AbstractString`: the name of the crystal structure file (include ".cif" or ".cssr") read from `PorousMaterials.PATH_TO_DATA * "structures/"`.
  * `check_charge_neutrality::Bool`: check for charge neutrality
  * `net_charge_tol::Float64`: when checking for charge neutrality, throw an error if the absolute value of the net charge is larger than this value.
  * `check_atom_and_charge_overlap::Bool`: throw an error if overlapping atoms are detected.
  * `remove_overlap::Bool`: remove identical atoms automatically. Identical atoms are the same element and overlap.

**Returns**

  * `framework::Framework`: A framework containing the crystal structure information

**Attributes**

  * `name::String`: name of crystal structure
  * `box::Box`: unit cell (Bravais Lattice)
  * `atoms::Array{LJSphere, 1}`: list of Lennard-Jones spheres in crystal unit cell
  * `charges::Array{PtCharge, 1}`: list of point charges in crystal unit cell


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/Crystal.jl#L11-L35' class='documenter-source'>source</a><br>

<a id='PorousMaterials.remove_overlapping_atoms_and_charges' href='#PorousMaterials.remove_overlapping_atoms_and_charges'>#</a>
**`PorousMaterials.remove_overlapping_atoms_and_charges`** &mdash; *Function*.



```
new_framework = remove_overlapping_atoms_and_charges(framework, overlap_tol=0.1, verbose=true)
```

Takes in a framework and returns a new framework with where overlapping atoms and overlapping charges were removed. i.e. if there is an overlapping pair, one in the pair is removed. For any atoms or charges to be removed, the species and charge, respectively, must be identical.

**Arguments**

  * `framework::Framework`: The framework containing the crystal structure information
  * `atom_overlap_tol::Float64`: The minimum distance between two atoms that is tolerated
  * `charge_overlap_tol::Float64`: The minimum distance between two charges that is tolerated

**Returns**

  * `new_framework::Framework`: A new framework where identical atoms have been removed.


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/Crystal.jl#L342-L357' class='documenter-source'>source</a><br>

<a id='PorousMaterials.strip_numbers_from_atom_labels!' href='#PorousMaterials.strip_numbers_from_atom_labels!'>#</a>
**`PorousMaterials.strip_numbers_from_atom_labels!`** &mdash; *Function*.



```
strip_numbers_from_atom_labels!(framework)
```

Strip numbers from labels for `framework.atoms`. Precisely, for `atom` in `framework.atoms`, find the first number that appears in `atom`. Remove this number and all following characters from `atom`. e.g. C12 –> C 	 Ba12A_3 –> Ba

**Arguments**

  * `framework::Framework`: The framework containing the crystal structure information


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/Crystal.jl#L445-L456' class='documenter-source'>source</a><br>

<a id='PorousMaterials.chemical_formula' href='#PorousMaterials.chemical_formula'>#</a>
**`PorousMaterials.chemical_formula`** &mdash; *Function*.



```
formula = chemical_formula(framework)
```

Find the irreducible chemical formula of a crystal structure.

**Arguments**

  * `framework::Framework`: The framework containing the crystal structure information

**Returns**

  * `formula::Dict{Symbol, Int}`: A dictionary with the irreducible chemical formula of a crystal structure


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/Crystal.jl#L475-L485' class='documenter-source'>source</a><br>

<a id='PorousMaterials.molecular_weight' href='#PorousMaterials.molecular_weight'>#</a>
**`PorousMaterials.molecular_weight`** &mdash; *Function*.



```
mass_of_framework = molecular_weight(framework)
```

Calculates the molecular weight of a unit cell of the framework in amu using information stored in `data/atomicmasses.csv`.

**Arguments**

  * `framework::Framework`: The framework containing the crystal structure information

**Returns**

  * `mass_of_framework::Float64`: The molecular weight of a unit cell of the framework in amu


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/Crystal.jl#L514-L525' class='documenter-source'>source</a><br>

<a id='PorousMaterials.crystal_density' href='#PorousMaterials.crystal_density'>#</a>
**`PorousMaterials.crystal_density`** &mdash; *Function*.



```
ρ = crystal_density(framework) # kg/m²
```

Compute the crystal density of a framework. Pulls atomic masses from [`read_atomic_masses`](functions.md#PorousMaterials.read_atomic_masses).

**Arguments**

  * `framework::Framework`: The framework containing the crystal structure information

**Returns**

  * `ρ::Float64`: The crystal density of a framework in kg/m³


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/Crystal.jl#L537-L547' class='documenter-source'>source</a><br>

<a id='PorousMaterials.replicate-Tuple{PorousMaterials.Framework,Tuple{Int64,Int64,Int64}}' href='#PorousMaterials.replicate-Tuple{PorousMaterials.Framework,Tuple{Int64,Int64,Int64}}'>#</a>
**`PorousMaterials.replicate`** &mdash; *Method*.



```
replicated_frame = replicate(framework, repfactors)
```

Replicates the atoms and charges in a `Framework` in positive directions to construct a new `Framework`. Note `replicate(framework, (1, 1, 1))` returns the same `Framework`.

**Arguments**

  * `framework::Framework`: The framework to replicate
  * `repfactors::Tuple{Int, Int, Int}`: The factors by which to replicate the crystal structure in each direction.

**Returns**

  * `replicated_frame::Framework`: Replicated framework


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/Crystal.jl#L221-L233' class='documenter-source'>source</a><br>

<a id='PorousMaterials.charged-Tuple{PorousMaterials.Framework}' href='#PorousMaterials.charged-Tuple{PorousMaterials.Framework}'>#</a>
**`PorousMaterials.charged`** &mdash; *Method*.



```
charged_flag = charged(framework, verbose=false) # true or false
```

Determine if a framework has point charges


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/Crystal.jl#L432-L436' class='documenter-source'>source</a><br>

<a id='PorousMaterials.write_cif' href='#PorousMaterials.write_cif'>#</a>
**`PorousMaterials.write_cif`** &mdash; *Function*.



```
write_cif(framework, filename)
```

Write a `framework::Framework` to a .cif file with `filename::String`. If `filename` does not include the .cif extension, it will automatically be added.


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/Crystal.jl#L553-L558' class='documenter-source'>source</a><br>

<a id='PorousMaterials.assign_charges' href='#PorousMaterials.assign_charges'>#</a>
**`PorousMaterials.assign_charges`** &mdash; *Function*.



```
new_framework = assign_charges(framework, charges, net_charge_tol=1e-5)
```

Assign charges to the atoms (represented by `LJSphere`s in `framework.atoms`) present in the framework. Pass a dictionary of charges that place charges according to the species of the atoms or pass an array of charges to assign to each atom, with the order of the array consistent with the order of `framework.atoms`.

If the framework already has charges, the charges are removed and new charges are added accordingly so that `length(framework.atoms) == length(framework.charges)`.

**Examples**

```
charges = Dict(:Ca => 2.0, :C => 1.0, :H => -1.0)
new_framework = assign_charges(framework, charges)
```

```
charges = [4.0, 2.0, -6.0] # framework.atoms is length 3
new_framework = assign_charges(framework, charges)
```

**Arguments**

  * `framework::Framework`: the framework to which we should add charges (not modified in

this function)

  * `charges::Union{Dict{Symbol, Float64}, Array{Float64, 1}}`: a dictionary that returns the

charge assigned to the species of atom or an array of charges to assign, with order consistent with the order in `framework.atoms` (units: electrons).

  * `net_charge_tol::Float64`: the net charge tolerated when asserting charge neutrality of

the resulting framework

**Returns**

  * `new_framework::Framework`: a new framework identical to the one passed except charges

are assigned.


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/Crystal.jl#L607-L641' class='documenter-source'>source</a><br>


##ElectrostaticsEnergetics.jl

<a id='PorousMaterials.electrostatic_potential_energy' href='#PorousMaterials.electrostatic_potential_energy'>#</a>
**`PorousMaterials.electrostatic_potential_energy`** &mdash; *Function*.



```
ϕ = electrostatic_potential_energy(framework, molecule, eparams, kvectors,
                                   eikar, eikbr, eikcr)
```

Compute the electrostatic potential energy of a molecule inside a framework.

The electrostatic potential is created by the point charges assigned to the framework atoms in `framework.charges`. Periodic boundary conditions are applied through the Ewald summation. The spurious self-interaction term is neglected here because we are looking at *differences* in energy in a Monte Carlo simulation.

Warning: it is assumed that the framework is replicated enough such that the nearest image convention can be applied for the short-range cutoff radius supplied in `eparams.sr_cutoff_r`.

**Arguments**

  * `framework::Framework`: Crystal structure (see `framework.charges` for charges)
  * `x::Array{Float64, 1}`: point (Cartesian coordinates) at which we compute the electrostatic   potential
  * `eparams::EwaldParams`: data structure containing Ewald summation settings
  * `kvectors::Array{Kvector, 1}`: array of k-vectors to include in the Fourier sum and their

corresponding weights indicating the contribution to the Fourier sum.

  * `eikra::OffsetArray{Complex{Float64}}`: array for storing e^{i * ka ⋅ r}; has indices   0:kreps[1] and corresponds to recip. vectors in a-direction
  * `eikrb::OffsetArray{Complex{Float64}}`: array for storing e^{i * kb ⋅ r}; has indices   -kreps[2]:kreps[2] and corresponds to recip. vectors in b-direction
  * `eikra::OffsetArray{Complex{Float64}}`: array for storing e^{i * kc ⋅ r}; has indices   -kreps[2]:kreps[1] and corresponds to recip. vectors in c-direction

**Returns**

electrostatic potential inside `framework` at point `x` (units: K)


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/ElectrostaticEnergetics.jl#L261-L292' class='documenter-source'>source</a><br>


```
ϕ = electrostatic_potential_energy(molecules, eparams, kvectors, eikar, eikbr, eikcr)
```

Compute the electrostatic potential energy of a system comprised of an array of `Molecule`s.

The EWald summation is used here in a double for loop; do not use this function for Monte Carlo simulations because it is computationally expensive.

Returns an `EwaldSum` type containing short-range and long-range contributions to the Ewald sum as well as the spurious self-interaction and intramolecular interactions. Access via (ϕ.sr, ϕ.lr, ϕ.self, ϕ.intra).

Units of energy: Kelvin

**Arguments**

  * `molecules::Array{Molecules, 1}`: array of molecules comprising the system.
  * `eparams::EwaldParams`: data structure containing Ewald summation settings
  * `kvectors::Array{Kvector, 1}`: array of k-vectors to include in the Fourier sum and their

corresponding weights indicating the contribution to the Fourier sum.

  * `eikra::OffsetArray{Complex{Float64}}`: array for storing e^{i * ka ⋅ r}; has indices   0:kreps[1] and corresponds to recip. vectors in a-direction
  * `eikrb::OffsetArray{Complex{Float64}}`: array for storing e^{i * kb ⋅ r}; has indices   -kreps[2]:kreps[2] and corresponds to recip. vectors in b-direction
  * `eikra::OffsetArray{Complex{Float64}}`: array for storing e^{i * kc ⋅ r}; has indices   -kreps[2]:kreps[1] and corresponds to recip. vectors in c-direction


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/ElectrostaticEnergetics.jl#L374-L399' class='documenter-source'>source</a><br>

<a id='PorousMaterials.precompute_kvec_wts' href='#PorousMaterials.precompute_kvec_wts'>#</a>
**`PorousMaterials.precompute_kvec_wts`** &mdash; *Function*.



```
kvectors = precompute_kvec_wts(eparams, max_mag_k_sqrd=Inf)
```

For speed, pre-compute the weights for each reciprocal lattice vector for the Ewald sum in Fourier space. This function takes advantage of the symmetry:     cos(-k⋅(x-xᵢ)) + cos(k⋅(x-xᵢ)) = 2 cos(k⋅(x-xᵢ))

If `max_mag_k_sqrd` is passed, k-vectors with a magnitude greater than `max_mag_k_sqrd` are not included.

**Arguments**

  * `eparams::EwaldParams`: data structure containing Ewald summation settings
  * `max_mag_k_sqrd::Float64`: cutoff for |k|² in Fourier sum; if passed, do not include

k-vectors with magnitude squared greater than this.

**Returns**

  * `kvectors::Array{Kvector, 1}`: array of k-vectors to include in the Fourier sum and their

corresponding weights indicating the contribution to the Fourier sum.


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/ElectrostaticEnergetics.jl#L133-L151' class='documenter-source'>source</a><br>

<a id='PorousMaterials.setup_Ewald_sum' href='#PorousMaterials.setup_Ewald_sum'>#</a>
**`PorousMaterials.setup_Ewald_sum`** &mdash; *Function*.



```
eparams, kvecs, eikar, eikbr, eikcr = setup_Ewald_sum(sr_cutoff_r, sim_box, ϵ=1e-6, verbose=false)
```

Given the short-range cutoff radius and simulation box, automatically compute Ewald convergence parameter and number of k-vector replications in Fourier space required for a given precision. Constructs and returns Ewald parameters data type with this information.

Also, pre-compute weights on k-vector contributions to Ewald sum in Fourier space.

Also, allocate OffsetArrays for storing e^{i * k ⋅ r} where r = x - xⱼ and k is a reciprocal lattice vector.

**Arguments**

**Returns**

  * `eparams::EwaldParams`: data structure containing Ewald summation settings
  * `kvectors::Array{Kvector, 1}`: array of k-vectors to include in the Fourier sum and their

corresponding weights indicating the contribution to the Fourier sum.

  * `eikar::OffsetArray{Complex{Float64}}`: array for storing e^{i * ka ⋅ r}; has indices   0:kreps[1] and corresponds to recip. vectors in a-direction
  * `eikbr::OffsetArray{Complex{Float64}}`: array for storing e^{i * kb ⋅ r}; has indices   -kreps[2]:kreps[2] and corresponds to recip. vectors in b-direction
  * `eikcr::OffsetArray{Complex{Float64}}`: array for storing e^{i * kc ⋅ r}; has indices   -kreps[2]:kreps[1] and corresponds to recip. vectors in c-direction


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/ElectrostaticEnergetics.jl#L187-L211' class='documenter-source'>source</a><br>

<a id='PorousMaterials.total_electrostatic_potential_energy' href='#PorousMaterials.total_electrostatic_potential_energy'>#</a>
**`PorousMaterials.total_electrostatic_potential_energy`** &mdash; *Function*.



```
total_ϕ = total_electrostatic_potential_energy(molecules, eparams, kvectors,
                                              eikar, eikbr, eikcr)
```

Explanation of total_electrostatic_potential_energy that doesn't use framework

**Arguments**

  * `molecules::Array{Molecule, 1}`:
  * `eparams::EwaldParams`:
  * `kvectors::Array{Kvector, 1}`:
  * `eikar::OffsetArray{Complex{Float64}}`:
  * `eikbr::OffsetArray{Complex{Float64}}`:
  * `eikcr::OffsetArray{Complex{Float64}}`:


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/ElectrostaticEnergetics.jl#L481-L494' class='documenter-source'>source</a><br>


```
total_ϕ = total_electrostatic_potential_energy(framework, molecules, eparams,
                                              kvectors, eikar, eikbr, eikcr)
```

Explanation of total_electrostatic_potential_energy that uses framework

**Arguments**

  * `framework::Framework`:
  * `molecules::Array{Molecule, 1}`:
  * `eparams::EwaldParams`:
  * `kvectors::Array{Kvector, 1}`:
  * `eikar::OffsetArray{Complex{Float64}}`:
  * `eikbr::OffsetArray{Complex{Float64}}`:
  * `eikcr::OffsetArray{Complex{Float64}}`:


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/ElectrostaticEnergetics.jl#L511-L525' class='documenter-source'>source</a><br>


##Energetics_Util.jl


```
    PotentialEnergy
    SystemPotentialEnergy
```


##EOS.jl


```
    PengRobinsonsGas
    calculate_properties
```


##Forcefield.jl


```
    LJForcefield
    replication_factors
    check_forcefield_coverage
```


##GCMC.jl

<a id='PorousMaterials.gcmc_simulation' href='#PorousMaterials.gcmc_simulation'>#</a>
**`PorousMaterials.gcmc_simulation`** &mdash; *Function*.



```
results, molecules = gcmc_simulation(framework, molecule, temperature, pressure,
                                     ljforcefield; n_sample_cycles=5000,
                                     n_burn_cycles=5000, sample_frequency=5,
                                     verbose=false, molecules=Molecule[],
                                     eos=:ideal)
```

Runs a grand-canonical (μVT) Monte Carlo simulation of the adsorption of a molecule in a framework at a particular temperature and pressure using a Lennard Jones force field.

A cycle is defined as max(20, number of adsorbates currently in the system) Markov chain proposals. Current Markov chain moves implemented are particle insertion/deletion and translation.

**Arguments**

  * `framework::Framework`: the porous crystal in which we seek to simulate adsorption
  * `temperature::Float64`: temperature of bulk gas phase in equilibrium with adsorbed phase   in the porous material. units: Kelvin (K)
  * `pressure::Float64`: pressure of bulk gas phase in equilibrium with adsorbed phase in the   porous material. units: bar
  * `molecule::Molecule`: a template of the adsorbate molecule of which we seek to simulate   the adsorption
  * `ljforcefield::LJForceField`: the molecular model used to describe the   energetics of the adsorbate-adsorbate and adsorbate-host van der Waals interactions.
  * `n_burn_cycles::Int`: number of cycles to allow the system to reach equilibrium before   sampling.
  * `n_sample_cycles::Int`: number of cycles used for sampling
  * `sample_frequency::Int`: during the sampling cycles, sample e.g. the number of adsorbed   gas molecules every this number of Markov proposals.
  * `verbose::Bool`: whether or not to print off information during the simulation.
  * `molecules::Array{Molecule, 1}`: a starting configuration of molecules in the framework.

Note that we assume these coordinates are Cartesian, i.e. corresponding to a unit box.

  * `eos::Symbol`: equation of state to use for calculation of fugacity from pressure. Default

is ideal gas, where fugacity = pressure.


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/GCMC.jl#L198-L233' class='documenter-source'>source</a><br>

<a id='PorousMaterials.adsorption_isotherm' href='#PorousMaterials.adsorption_isotherm'>#</a>
**`PorousMaterials.adsorption_isotherm`** &mdash; *Function*.



```
results = adsorption_isotherm(framework, molecule, temperature, pressures,
                              ljforcefield; n_sample_cycles=100000,
                              n_burn_cycles=10000, sample_frequency=25,
                              verbose=false, molecules=Molecule[],
                              ewald_precision=1e-6, eos=:ideal)
```

Run a set of grand-canonical (μVT) Monte Carlo simulations in parallel. Arguments are the same as [`gcmc_simulation`](functions.md#PorousMaterials.gcmc_simulation), as this is the function run in parallel behind the scenes. The only exception is that we pass an array of pressures. To give Julia access to multiple cores, run your script as `julia -p 4 mysim.jl` to allocate e.g. four cores. See [Parallel Computing](https://docs.julialang.org/en/stable/manual/parallel-computing/#Parallel-Computing-1).


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/GCMC.jl#L159-L171' class='documenter-source'>source</a><br>

<a id='PorousMaterials.stepwise_adsorption_isotherm' href='#PorousMaterials.stepwise_adsorption_isotherm'>#</a>
**`PorousMaterials.stepwise_adsorption_isotherm`** &mdash; *Function*.



```
results = stepwise_adsorption_isotherm(framework, molecule, temperature, pressures,
                              ljforcefield; n_sample_cycles=100000,
                              n_burn_cycles=10000, sample_frequency=10,
                              verbose=true, molecules=Molecule[],
                              ewald_precision=1e-6, eos=:ideal)
```

Run a set of grand-canonical (μVT) Monte Carlo simulations in series. Arguments are the same as [`gcmc_simulation`](functions.md#PorousMaterials.gcmc_simulation), as this is the function run behind the scenes. An exception is that we pass an array of pressures. The adsorption isotherm is computed step- wise, where the ending configuration from the previous simulation (array of molecules) is passed into the next simulation as a starting point. The ordering of `pressures` is honored. By giving each simulation a good starting point, (if the next pressure does not differ significantly from the previous pressure), we can reduce the number of burn cycles required to reach equilibrium in the Monte Carlo simulation. Also see [`adsorption_isotherm`](functions.md#PorousMaterials.adsorption_isotherm) which runs the μVT simulation at each pressure in parallel.


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/GCMC.jl#L121-L137' class='documenter-source'>source</a><br>

<a id='PorousMaterials.gcmc_result_savename' href='#PorousMaterials.gcmc_result_savename'>#</a>
**`PorousMaterials.gcmc_result_savename`** &mdash; *Function*.



```
file_save_name = gcmc_result_savename(framework_name, molecule_species
                                    ljforcefield_name, temperature, pressure,
                                    n_burn_cycles, n_sample_cycles)
```

Determine the name of files saved during the GCMC simulation, be molecule positions or results. It uses many pieces of information from the simulation to ensure the file name accurately describes what it holds.

**Arguments**

  * `framework_name::AbstractString`: The porous crystal being tested
  * `molecule_species::Symbol`: The molecule being tested inside the porous crystal
  * `ljforcefield_name::AbstractString`: The molecular model being used in this   simulation to describe intermolecular Van der Waals interactions
  * `temperature::Float64`: The temperature used in the simulation units: Kelvin (K)
  * `pressure::Float64`: The pressure used in the simulation units: bar
  * `n_burn_cycles::Int`: The number of burn cycles used in this simulation
  * `n_sample_cycles::Int`: The number of sample cycles used in this simulation


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/GCMC.jl#L656-L674' class='documenter-source'>source</a><br>


##Grid.jl


```
    Grid
    apply_periodic_boundary_condition
    write_cube
    read_cube
    energy_grid
```


##Henry.jl


```
    henry_coefficient
    henry_result_savename
```


##Matter.jl

<a id='PorousMaterials.LJSphere' href='#PorousMaterials.LJSphere'>#</a>
**`PorousMaterials.LJSphere`** &mdash; *Type*.



Data structure for a Lennard-Jones sphere, containing its species and position in  fractional coordinates.

**Example use**

```
ljs = LJSphere(:C, [0.0, 0.0, 0.0])
```

**Attributes**

  * `species::Symbol`: atom species name, e.g. `:C`
  * `xf::Array{Float64, 1}`: fractional coordinates, e.g. `[1.0, 0.0, 4.0]`.


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/Matter.jl#L3-L13' class='documenter-source'>source</a><br>

<a id='PorousMaterials.PtCharge' href='#PorousMaterials.PtCharge'>#</a>
**`PorousMaterials.PtCharge`** &mdash; *Type*.



Point charge data structure indicates its charge and position in fractional coordinates.

**Example use**

```
ptc = PtCharge(-0.2, [0.0, 0.0, 0.0])
```

**Attributes**

  * `q::Float64`: signed magnitude of charge (units: electrons), e.g. `1.0`
  * `xf::Array{Float64, 1}`: fractional coordinates, e.g. `[1.0, 0.0, 4.0]`.


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/Matter.jl#L23-L32' class='documenter-source'>source</a><br>


##MChelpers.jl


```
    insert_molecule!
    delete_molecule!
    translate_molecule!
    reinsert_molecule!
    rotatable
```


##Misc.jl

<a id='PorousMaterials.read_xyz' href='#PorousMaterials.read_xyz'>#</a>
**`PorousMaterials.read_xyz`** &mdash; *Function*.



```
atoms, x = read_xyz(filename)
```

Return the list of `atoms` (Array{Symbol, 1}) and their Cartesian coordinates `x::Array{Float64, 2}` as stored in the .xyz file. `x[:, k]` will return Cartesian coords of the kth atom.

**Arguments**

  * `filename::AbstractString`: The filename of the .xyz file

**Returns**

  * `atoms::Array{Symbol, 1}`: An array of atoms stored as symbols e.g. [:H, :H, :O] read

from the .xyz file.

  * `x::Array{Float64, 2}`: The Cartesian coordinates of the atoms. `x[:, k]` will return cartesian coordinates of the k-th atom


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/Misc.jl#L1-L15' class='documenter-source'>source</a><br>

<a id='PorousMaterials.read_cpk_colors' href='#PorousMaterials.read_cpk_colors'>#</a>
**`PorousMaterials.read_cpk_colors`** &mdash; *Function*.



```
atom_colors = read_cpk_colors()
```

Read in CPK color scheme for atoms. Return `atom_colors::Dict{Symbol, Tuple{Int, Int, Int}}` such that `atom_colors[":C"]` gives RGB code for carbon as a tuple, `(144, 144, 144)`. https://en.wikipedia.org/wiki/CPK_coloring

**Returns**

  * `atom_colors::Dict{Symbol, Tuple{Int, Int, Int}}`: A dictionary linking an element symbol to its' corresponding CPK color in RGB


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/Misc.jl#L68-L77' class='documenter-source'>source</a><br>

<a id='PorousMaterials.read_atomic_radii' href='#PorousMaterials.read_atomic_radii'>#</a>
**`PorousMaterials.read_atomic_radii`** &mdash; *Function*.



```
atomic_radii = read_atomic_radii()
```

Return `atomic_radii::Dict{Symbol, Float64}`, where `atom_masses[":C"]` gives the atomic radii of carbon (10.87 Angstrom).

**Returns**

  * `atomic_radii::Dict{Symbol, Float64}`: A dictionary linking an element symbol to its' corresponding atomic radius


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/Misc.jl#L87-L95' class='documenter-source'>source</a><br>

<a id='PorousMaterials.read_atomic_masses' href='#PorousMaterials.read_atomic_masses'>#</a>
**`PorousMaterials.read_atomic_masses`** &mdash; *Function*.



```
atomic_masses = read_atomic_masses()
```

Read the `data/atomicmasses.csv` file to construct a dictionary of atoms and their atomic masses in amu.

**Returns**

  * `atomic_masses::Dict{Symbol, Float64}`: A dictionary containing the atomic masses of each atom stored in `data/atomicmasses.csv`


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/Misc.jl#L107-L115' class='documenter-source'>source</a><br>

<a id='PorousMaterials.write_to_xyz' href='#PorousMaterials.write_to_xyz'>#</a>
**`PorousMaterials.write_to_xyz`** &mdash; *Function*.



```
write_to_xyz(atoms, x, filename; comment="")
write_to_xyz(molecules, box, filename; comment="")
write_to_xyz(framework, filename; comment="")
```

Write a molecule, framework, or array of atoms & positions to an .xyz file.

**Arguments**

  * `atoms::Array{Symbol, 1}`: An array of atoms stored as symbols e.g. [:H, :H, :O]
  * `x::Array{Float64, 2}`: The Cartesian coordinates of the atoms.

`x[:, k]` contains Cartesian coordinates of the k-th atom

  * `molecules::Array{Molecule, 1}`: an array of molecules whose atoms to write to .xyz
  * `framework::Framework`: a crystal structure whose atoms to write to .xyz
  * `filename::AbstractString`: The filename of the .xyz file. (".xyz" appended automatically

if the extension is not provided.) (absolute path)

  * `comment::AbstractString`: comment if you'd like to write to the file.


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/Misc.jl#L35-L51' class='documenter-source'>source</a><br>


##Molecules.jl


```
    translate_to!
    rotate!
    rotation_matrix
    rand_point_on_unit_sphere
    charged(::Molecule, ::Bool)
```


##NearestImage.jl


```
    nearest_image!
    nearest_r²
    nearest_r
```


##VdWEnergetics.jl

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


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/VdWEnergetics.jl#L4-L18' class='documenter-source'>source</a><br>

<a id='PorousMaterials.vdw_energy' href='#PorousMaterials.vdw_energy'>#</a>
**`PorousMaterials.vdw_energy`** &mdash; *Function*.



```
energy = vdw_energy(framework, molecule, ljforcefield)
```

Calculates the van der Waals interaction energy between a molecule and a framework. Applies the nearest image convention to find the closest replicate of a specific atom.

WARNING: it is assumed that the framework is replicated sufficiently such that the nearest image convention can be applied. See [`replicate`](functions.md#PorousMaterials.replicate).

**Arguments**

  * `framework::Framework`: Crystal structure
  * `molecule::Molecule`: adsorbate (includes position/orientation/atoms)
  * `ljforcefield::LJForceField`: Lennard Jones force field

**Returns**

  * `energy::Float64`: Van der Waals interaction energy


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/VdWEnergetics.jl#L36-L52' class='documenter-source'>source</a><br>


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


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/VdWEnergetics.jl#L63-L78' class='documenter-source'>source</a><br>

<a id='PorousMaterials.vdw_energy_no_PBC' href='#PorousMaterials.vdw_energy_no_PBC'>#</a>
**`PorousMaterials.vdw_energy_no_PBC`** &mdash; *Function*.



Assumes unit cell box is a unit cube and no periodic boundary conditions are applied.


<a target='_blank' href='https://github.com/ahyork/PorousMaterials.jl/blob/e979536e56ffe9001b5652a7ade389634560a6c7/src/VdWEnergetics.jl#L135-L138' class='documenter-source'>source</a><br>

