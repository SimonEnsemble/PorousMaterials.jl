
<a id='Boxes-1'></a>

## Boxes

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


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/Box.jl#L1-L20' class='documenter-source'>source</a><br>

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


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/Box.jl#L118-L132' class='documenter-source'>source</a><br>


```
replicated_frame = replicate(framework, repfactors)
```

Replicates the atoms and charges in a `Framework` in positive directions to construct a new `Framework`. Note `replicate(framework, (1, 1, 1))` returns the same `Framework`.

**Arguments**

  * `framework::Framework`: The framework to replicate
  * `repfactors::Tuple{Int, Int, Int}`: The factors by which to replicate the crystal structure in each direction.

**Returns**

  * `replicated_frame::Framework`: Replicated framework


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/Crystal.jl#L222-L234' class='documenter-source'>source</a><br>

<a id='PorousMaterials.UnitCube' href='#PorousMaterials.UnitCube'>#</a>
**`PorousMaterials.UnitCube`** &mdash; *Function*.



```
unit_cube = UnitCube()
```

This function generates a unit cube, each side is 1.0 Angstrom long, and all the corners are right angles.


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/Box.jl#L110-L115' class='documenter-source'>source</a><br>

<a id='PorousMaterials.write_vtk' href='#PorousMaterials.write_vtk'>#</a>
**`PorousMaterials.write_vtk`** &mdash; *Function*.



```
write_vtk(box, filename; verbose=true)
write_vtk(framework)
```

Write a `Box` to a .vtk file for visualizing e.g. the unit cell boundary of a crystal. If a `Framework` is passed, the `Box` of that framework is written to a file that is the same as the crystal structure filename but with a .vtk extension.

Appends ".vtk" extension to `filename` automatically if not passed.

**Arguments**

  * `box::Box`: a Bravais lattice
  * `filename::AbstractString`: filename of the .vtk file output (absolute path)
  * `framework::Framework`: A framework containing the crystal structure information


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/Box.jl#L138-L152' class='documenter-source'>source</a><br>


<a id='Crystals-1'></a>

## Crystals

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

  * `filename::AbstractString`: the name of the crystal structure file (include ".cif" or ".cssr") read from `joinpath(PorousMaterials.PATH_TO_DATA, "structures")`.
  * `check_charge_neutrality::Bool`: check for charge neutrality
  * `net_charge_tol::Float64`: when checking for charge neutrality, throw an error if the absolute value of the net charge is larger than this value.
  * `check_atom_and_charge_overlap::Bool`: throw an error if overlapping atoms are detected.
  * `remove_overlap::Bool`: remove identical atoms automatically. Identical atoms are the same element atoms which overlap.

**Returns**

  * `framework::Framework`: A framework containing the crystal structure information

**Attributes**

  * `name::AbstractString`: name of crystal structure
  * `box::Box`: unit cell (Bravais Lattice)
  * `atoms::Atoms`: list of Atoms in crystal unit cell
  * `charges::Charges`: list of point charges in crystal unit cell


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/Crystal.jl#L11-L35' class='documenter-source'>source</a><br>

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


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/Crystal.jl#L364-L379' class='documenter-source'>source</a><br>

<a id='PorousMaterials.strip_numbers_from_atom_labels!' href='#PorousMaterials.strip_numbers_from_atom_labels!'>#</a>
**`PorousMaterials.strip_numbers_from_atom_labels!`** &mdash; *Function*.



```
strip_numbers_from_atom_labels!(framework)
```

Strip numbers from labels for `framework.atoms`. Precisely, for `atom` in `framework.atoms`, find the first number that appears in `atom`. Remove this number and all following characters from `atom`. e.g. C12 –> C 	 Ba12A_3 –> Ba

**Arguments**

  * `framework::Framework`: The framework containing the crystal structure information


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/Crystal.jl#L477-L488' class='documenter-source'>source</a><br>

<a id='PorousMaterials.chemical_formula' href='#PorousMaterials.chemical_formula'>#</a>
**`PorousMaterials.chemical_formula`** &mdash; *Function*.



```
formula = chemical_formula(framework, verbose=false)
```

Find the irreducible chemical formula of a crystal structure.

**Arguments**

  * `framework::Framework`: The framework containing the crystal structure information
  * `verbose::Bool`: If `true`, will print the chemical formula as well

**Returns**

  * `formula::Dict{Symbol, Int}`: A dictionary with the irreducible chemical formula of a crystal structure


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/Crystal.jl#L505-L516' class='documenter-source'>source</a><br>

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


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/Crystal.jl#L545-L556' class='documenter-source'>source</a><br>

<a id='PorousMaterials.crystal_density' href='#PorousMaterials.crystal_density'>#</a>
**`PorousMaterials.crystal_density`** &mdash; *Function*.



```
ρ = crystal_density(framework) # kg/m²
```

Compute the crystal density of a framework. Pulls atomic masses from [`read_atomic_masses`](../other/docs.md#PorousMaterials.read_atomic_masses).

**Arguments**

  * `framework::Framework`: The framework containing the crystal structure information

**Returns**

  * `ρ::Float64`: The crystal density of a framework in kg/m³


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/Crystal.jl#L568-L578' class='documenter-source'>source</a><br>

<a id='PorousMaterials.replicate-Tuple{Framework,Tuple{Int64,Int64,Int64}}' href='#PorousMaterials.replicate-Tuple{Framework,Tuple{Int64,Int64,Int64}}'>#</a>
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


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/Crystal.jl#L222-L234' class='documenter-source'>source</a><br>

<a id='PorousMaterials.charged-Tuple{Framework}' href='#PorousMaterials.charged-Tuple{Framework}'>#</a>
**`PorousMaterials.charged`** &mdash; *Method*.



```
charged_flag = charged(framework, verbose=false) # true or false
```

Determine if a framework has point charges


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/Crystal.jl#L464-L468' class='documenter-source'>source</a><br>

<a id='PorousMaterials.write_cif' href='#PorousMaterials.write_cif'>#</a>
**`PorousMaterials.write_cif`** &mdash; *Function*.



```
write_cif(framework, filename)
```

Write a `framework::Framework` to a .cif file with `filename::AbstractString`. If `filename` does not include the .cif extension, it will automatically be added.


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/Crystal.jl#L584-L589' class='documenter-source'>source</a><br>

<a id='PorousMaterials.assign_charges' href='#PorousMaterials.assign_charges'>#</a>
**`PorousMaterials.assign_charges`** &mdash; *Function*.



```
new_framework = assign_charges(framework, charges, net_charge_tol=1e-5)
```

Assign charges to the atoms present in the framework. Pass a dictionary of charges that place charges according to the species of the atoms or pass an array of charges to assign to each atom, with the order of the array consistent with the order of `framework.atoms`.

If the framework already has charges, the charges are removed and new charges are added accordingly so that `framework.atoms.n_atoms == framework.charges.n_charges`.

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


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/Crystal.jl#L637-L671' class='documenter-source'>source</a><br>


<a id='Grids-1'></a>

## Grids

<a id='PorousMaterials.Grid' href='#PorousMaterials.Grid'>#</a>
**`PorousMaterials.Grid`** &mdash; *Type*.



Data structure for a regular [equal spacing between points in each coordinate] grid of points superimposed on a unit cell box (`Box`). Each grid point has data, `data`, associated with it, of type `T`, stored in a 3D array.

**Attributes**

  * `box::Box`: describes Bravais lattice over which a grid of points is super-imposed. grid points on all faces are included.
  * `n_pts::Tuple{Int, Int, Int}`: number of grid points in x, y, z directions. 0 and 1 fractional coordinates are included.
  * `data::Array{T, 3}`: three dimensional array conaining data associated with each grid point.
  * `units::Symbol`: the units associated with each data point.
  * `origin::Array{Float64, 1}`: the origin of the grid.


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/Grid.jl#L1-L11' class='documenter-source'>source</a><br>

<a id='PorousMaterials.apply_periodic_boundary_condition!' href='#PorousMaterials.apply_periodic_boundary_condition!'>#</a>
**`PorousMaterials.apply_periodic_boundary_condition!`** &mdash; *Function*.



```
apply_periodic_boundary_condition!(molecule)
```

Check if the `center_of_mass` of a `Molecule` is outside of a `Box`. If so, apply periodic boundary conditions and translate the center of mass of the `Molecule` (and its atoms and point charges) so that it is inside of the `Box`.

**Arguments**

  * `molecule::Molecule`: A molecule we're interested in seeing if its' center of mass falls within `simulation_box`


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/MChelpers.jl#L50-L59' class='documenter-source'>source</a><br>

<a id='PorousMaterials.write_cube' href='#PorousMaterials.write_cube'>#</a>
**`PorousMaterials.write_cube`** &mdash; *Function*.



```
write_cube(grid, filename, verbose=true)
```

Write grid to a .cube file format. This format is described here: http://paulbourke.net/dataformats/cube/ The atoms of the unit cell are not printed in the .cube. Instead, use .xyz files to also visualize atoms.

**Arguments**

  * `grid::Grid`: grid with associated data at each grid point.
  * `filename::AbstractString`: name of .cube file to which we write the grid; this is relative to `PorousMaterials.PATH_TO_DATA`/grids/.
  * `verbose::Bool`: print name of file after writing.


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/Grid.jl#L36-L47' class='documenter-source'>source</a><br>

<a id='PorousMaterials.read_cube' href='#PorousMaterials.read_cube'>#</a>
**`PorousMaterials.read_cube`** &mdash; *Function*.



```
grid = read_cube(filename)
```

Read a .cube file and return a populated `Grid` data structure.

**Arguments**

  * `filename::AbstractString`: name of .cube file to which we write the grid; this is relative to `PorousMaterials.PATH_TO_DATA`grids/.

**Returns**

  * `grid::Grid`: A grid data structure


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/Grid.jl#L89-L99' class='documenter-source'>source</a><br>

<a id='PorousMaterials.energy_grid' href='#PorousMaterials.energy_grid'>#</a>
**`PorousMaterials.energy_grid`** &mdash; *Function*.



```
grid = energy_grid(framework, molecule, ljforcefield; n_pts=(50, 50, 50), temperature=298.0, n_rotations=750)
```

Superimposes a regular grid of points (regularly spaced in fractional coordinates of the `framework.box`) over the unit cell of a crystal, with `n_gridpts` dictating the number of grid points in the a, b, c directions (including 0 and 1 fractional coords). The fractional coordinates 0 and 1 are included in the grid, although they are redundant. Then, at each grid point, calculate the ensemble average potential energy of the molecule when its mass is centered at that point. The average is taken over Boltzmann-weighted rotations.

The ensemble average is a Boltzmann average over rotations:  - R T log ⟨e⁻ᵇᵁ⟩

**Arguments**

  * `framework::Framework`: crystal in which we seek to compute an energy grid for a molecule. `grid.box` will be `framework.box`.
  * `molecule::Molecule`: molecule for which we seek an energy grid
  * `ljforcefield::LJForceField`: molecular model for computing molecule-framework interactions
  * `n_pts::Tuple{Int, Int, Int}=(50,50,50)`: number of grid points in each fractional coordinate dimension, including endpoints (0, 1)
  * `n_rotations::Int`: number of random rotations to conduct in a Monte Carlo simulation for finding the free energy of a molecule centered at a given grid point.

This is only relevant for molecules that are comprised of more than one Lennard Jones sphere.

  * `temperature::Float64`: the temperature at which to compute the free energy for molecules where rotations are required. Lower temperatures overemphasize the minimum potential energy rotational conformation at that point.
  * `units::Symbol`: either `:K` or `:kJ_mol`, the units in which the energy should be stored in the returned `Grid`.
  * `center::Bool`: shift coords of grid so that the origin is the center of the unit cell `framework.box`.
  * `verbose::Bool=true`: print some information.

**Returns**

  * `grid::Grid`: A grid data structure containing the potential energy of the system


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/Grid.jl#L155-L178' class='documenter-source'>source</a><br>

