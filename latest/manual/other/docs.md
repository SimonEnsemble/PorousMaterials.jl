
<a id='PATH_TO_DATA-Control-1'></a>

## PATH_TO_DATA Control

<a id='PorousMaterials.set_path_to_data' href='#PorousMaterials.set_path_to_data'>#</a>
**`PorousMaterials.set_path_to_data`** &mdash; *Function*.



```
set_path_to_data("user/path/to/data")
set_path_to_data()
```

Sets PorousMaterials `PATH_TO_DATA` variable which dictates where crystal, forcefield, and molecule files are loaded from. This function allows the user to set `PATH_TO_DATA` manually to any directory or to a "/data/" folder within their current directory. This function WILL change the `PATH_TO_DATA` regardless of whether or not the path exists, but will give a warning alerting the user that PorousMaterials cannot load files from the chosen path.

**Arguments**

  * `new_path_to_data::String`: The desired `PATH_TO_DATA` in string form.


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/PorousMaterials.jl#L30-L43' class='documenter-source'>source</a><br>

<a id='PorousMaterials.set_tutorial_mode' href='#PorousMaterials.set_tutorial_mode'>#</a>
**`PorousMaterials.set_tutorial_mode`** &mdash; *Function*.



```
set_tutorial_mode()
```

Places PorousMaterials in "Tutorial Mode". It changes the `PATH_TO_DATA` variable to the directory where the PorousMaterials test data is stored. It can be used to follow examples shown in the README. It displays a warning so that the user knows They are no longer using their own data.


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/PorousMaterials.jl#L60-L67' class='documenter-source'>source</a><br>


<a id='Reading-in-Atomic-Values-1'></a>

## Reading in Atomic Values

<a id='PorousMaterials.read_atomic_radii' href='#PorousMaterials.read_atomic_radii'>#</a>
**`PorousMaterials.read_atomic_radii`** &mdash; *Function*.



```
atomic_radii = read_atomic_radii()
```

Return `atomic_radii::Dict{Symbol, Float64}`, where `atom_masses[":C"]` gives the atomic radii of carbon (10.87 Angstrom).

**Returns**

  * `atomic_radii::Dict{Symbol, Float64}`: A dictionary linking an element symbol to its' corresponding atomic radius


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/Misc.jl#L88-L96' class='documenter-source'>source</a><br>

<a id='PorousMaterials.read_atomic_masses' href='#PorousMaterials.read_atomic_masses'>#</a>
**`PorousMaterials.read_atomic_masses`** &mdash; *Function*.



```
atomic_masses = read_atomic_masses()
```

Read the `data/atomicmasses.csv` file to construct a dictionary of atoms and their atomic masses in amu.

**Returns**

  * `atomic_masses::Dict{Symbol, Float64}`: A dictionary containing the atomic masses of each atom stored in `data/atomicmasses.csv`


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/Misc.jl#L108-L116' class='documenter-source'>source</a><br>

<a id='PorousMaterials.read_cpk_colors' href='#PorousMaterials.read_cpk_colors'>#</a>
**`PorousMaterials.read_cpk_colors`** &mdash; *Function*.



```
atom_colors = read_cpk_colors()
```

Read in CPK color scheme for atoms. Return `atom_colors::Dict{Symbol, Tuple{Int, Int, Int}}` such that `atom_colors[":C"]` gives RGB code for carbon as a tuple, `(144, 144, 144)`. https://en.wikipedia.org/wiki/CPK_coloring

**Returns**

  * `atom_colors::Dict{Symbol, Tuple{Int, Int, Int}}`: A dictionary linking an element symbol to its' corresponding CPK color in RGB


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/Misc.jl#L69-L78' class='documenter-source'>source</a><br>


<a id='Using-.xyz-files-1'></a>

## Using .xyz files

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


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/Misc.jl#L1-L15' class='documenter-source'>source</a><br>

<a id='PorousMaterials.write_xyz' href='#PorousMaterials.write_xyz'>#</a>
**`PorousMaterials.write_xyz`** &mdash; *Function*.



```
write_xyz(atoms, x, filename; comment="")
write_xyz(molecules, box, filename; comment="")
write_xyz(framework, filename; comment="", center=false)
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
  * `center::Bool`: shift atoms so that origin is the center of the `framework.box`


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/Misc.jl#L35-L52' class='documenter-source'>source</a><br>

