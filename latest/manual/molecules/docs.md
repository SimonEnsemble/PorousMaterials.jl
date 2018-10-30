
<a id='Molecules-1'></a>

## Molecules

<a id='PorousMaterials.Molecule' href='#PorousMaterials.Molecule'>#</a>
**`PorousMaterials.Molecule`** &mdash; *Type*.



Data structure for a molecule/adsorbate.

**Attributes**

  * `species::Symbol`: Species of molecule, e.g. `:CO2`
  * `atoms::Atoms`: array of Lennard-Jones spheres comprising the molecule
  * `charges::Charges`: array of point charges comprising the molecule
  * `xf_com::Array{Float64, 1}`: center of mass of the molecule in fractional coordinates


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/845505808b4e0fb9638d52f835a7b7cb0cde5b8f/src/Molecules.jl#L1-L9' class='documenter-source'>source</a><br>

<a id='PorousMaterials.translate_to!' href='#PorousMaterials.translate_to!'>#</a>
**`PorousMaterials.translate_to!`** &mdash; *Function*.



```
translate_to!(molecule, xf)
translate_to!(molecule, x, box)
```

Translate a molecule a molecule to point `xf` in fractional coordinate space or to `x` in Cartesian coordinate space. For the latter, a unit cell box is required for context. The molecule is translated such that its center of mass is at `xf`/x`.

**Arguments**

  * `molecule::Molecule`: The molecule which will be translated to `xf`
  * `xf::Array{Float64, 1}`: A vector containing the coordinates of the final destination of the molecule


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/845505808b4e0fb9638d52f835a7b7cb0cde5b8f/src/Molecules.jl#L177-L188' class='documenter-source'>source</a><br>

<a id='PorousMaterials.rotate!' href='#PorousMaterials.rotate!'>#</a>
**`PorousMaterials.rotate!`** &mdash; *Function*.



```
rotate!(molecule, box)
```

Conduct a random rotation of the molecule about its center of mass. The box is needed because the molecule contains only its fractional coordinates.

**Arguments**

  * `molecule::Molecule`: The molecule which will be subject to a random rotation
  * `box::Box`: The molecule only contains fractional coordinates, so the box is needed for a correct rotation


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/845505808b4e0fb9638d52f835a7b7cb0cde5b8f/src/Molecules.jl#L258-L267' class='documenter-source'>source</a><br>

<a id='PorousMaterials.rotation_matrix' href='#PorousMaterials.rotation_matrix'>#</a>
**`PorousMaterials.rotation_matrix`** &mdash; *Function*.



```
r = rotation_matrix()
```

Generate a 3x3 random rotation matrix `r` such that when a point `x` is rotated using this rotation matrix via `r * x`, this point `x` is placed at a uniform random distributed position on the surface of a sphere of radius `norm(x)`. See James Arvo. Fast Random Rotation Matrices.

https://pdfs.semanticscholar.org/04f3/beeee1ce89b9adf17a6fabde1221a328dbad.pdf

**Returns**

  * `r::Array{Float64, 2}`: A 3x3 random rotation matrix


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/845505808b4e0fb9638d52f835a7b7cb0cde5b8f/src/Molecules.jl#L234-L244' class='documenter-source'>source</a><br>

<a id='PorousMaterials.rand_point_on_unit_sphere' href='#PorousMaterials.rand_point_on_unit_sphere'>#</a>
**`PorousMaterials.rand_point_on_unit_sphere`** &mdash; *Function*.



```
u = rand_point_on_unit_sphere()
```

Generate a unit vector with a random orientation.

**Returns**

  * `u::Array{Float64, 1}`: A unit vector with a random orientation


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/845505808b4e0fb9638d52f835a7b7cb0cde5b8f/src/Molecules.jl#L217-L224' class='documenter-source'>source</a><br>

<a id='PorousMaterials.charged-Tuple{Molecule}' href='#PorousMaterials.charged-Tuple{Molecule}'>#</a>
**`PorousMaterials.charged`** &mdash; *Method*.



```
charged_flag = charged(molecule, verbose=false)
```

Determine if a molecule has point charges

**Arguments**

  * `molecule::Molecule`: The molecule which will be checked for charges
  * `verbose::Bool`: Will print result if `true`

**Returns**

  * `charged_flag::Bool`: `true` if molecule is charged, `false` otherwise


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/845505808b4e0fb9638d52f835a7b7cb0cde5b8f/src/Molecules.jl#L340-L351' class='documenter-source'>source</a><br>


<a id='Molecular-Movement-1'></a>

## Molecular Movement

<a id='PorousMaterials.insert_molecule!' href='#PorousMaterials.insert_molecule!'>#</a>
**`PorousMaterials.insert_molecule!`** &mdash; *Function*.



```
insert_molecule!(molecules, box, template)
```

Inserts an additional adsorbate molecule into the simulation box using the template provided. The center of mass of the molecule is chosen at a uniform random position in the simulation box. A uniformly random orientation of the molecule is chosen by rotating about the center of mass.

**Arguments**

  * `molecules::Array{Molecule, 1}`: An array of Molecule objects
  * `box::Box`: The simulation box
  * `template::Molecule`: A template molecule used as reference when inserting molecules


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/845505808b4e0fb9638d52f835a7b7cb0cde5b8f/src/MChelpers.jl#L8-L19' class='documenter-source'>source</a><br>

<a id='PorousMaterials.delete_molecule!' href='#PorousMaterials.delete_molecule!'>#</a>
**`PorousMaterials.delete_molecule!`** &mdash; *Function*.



```
delete_molecule!(molecule_id, molecules)
```

Removes a random molecule from the current molecules in the framework. molecule_id decides which molecule will be deleted, for a simulation, it must be a randomly generated value

**Arguments**

  * `molecule_id::Int`: The molecule ID is used to determine which molecule in `molecules` should be removed
  * `molecules::Array{Molecule, 1}`: An array of Molecule objects


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/845505808b4e0fb9638d52f835a7b7cb0cde5b8f/src/MChelpers.jl#L35-L45' class='documenter-source'>source</a><br>

<a id='PorousMaterials.translate_molecule!' href='#PorousMaterials.translate_molecule!'>#</a>
**`PorousMaterials.translate_molecule!`** &mdash; *Function*.



```
translate_molecule!(molecule, box)
```

Perturbs the Cartesian coordinates of a molecule about its center of mass by a random vector of max length δ. Applies periodic boundary conditions to keep the molecule inside the simulation box. Returns a deep copy of the old molecule in case it needs replaced if the Monte Carlo proposal is rejected.

**Arguments**

  * `molecule::Molecule`: The molecule we want to perturb
  * `box::Box`: The simulation box

**Returns**

  * `old_molecule::Molecule`: The old molecule in case the MC proposal is rejected


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/845505808b4e0fb9638d52f835a7b7cb0cde5b8f/src/MChelpers.jl#L84-L98' class='documenter-source'>source</a><br>

<a id='PorousMaterials.reinsert_molecule!' href='#PorousMaterials.reinsert_molecule!'>#</a>
**`PorousMaterials.reinsert_molecule!`** &mdash; *Function*.



```
reinsert_molecule(molecule, box)
```

Move molecule to a new center of mass randomly distrubted in the unit cell and choose a random orientation for it. Return a deep copy of the starting molecule for possible restoration. This MC move can be viewed as a more aggressive `translate_molecule!`.

**Arguments**

  * `molecule::Molecule`: The molecule we want to perturb
  * `box::Box`: The simulation box


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/845505808b4e0fb9638d52f835a7b7cb0cde5b8f/src/MChelpers.jl#L113-L123' class='documenter-source'>source</a><br>

<a id='PorousMaterials.rotatable' href='#PorousMaterials.rotatable'>#</a>
**`PorousMaterials.rotatable`** &mdash; *Function*.



```
need_to_rotate = rotatable(molecule)
```

Determines whether or not a given molecule needs to be rotated. For example, rotating a single atom isn't necessary.

**Arguments**

  * `molecule::Molecule`: The molecule being tested. This function determines if a   rotation of this molecule will do anything.

**Returns**

  * `is_rotatable::Bool`: A boolean describing whether or not rotating the molecule   will alter its interactions with other molecules


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/845505808b4e0fb9638d52f835a7b7cb0cde5b8f/src/MChelpers.jl#L142-L155' class='documenter-source'>source</a><br>

