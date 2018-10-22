
<a id='Matter-1'></a>

## Matter

<a id='PorousMaterials.Atoms' href='#PorousMaterials.Atoms'>#</a>
**`PorousMaterials.Atoms`** &mdash; *Type*.



Data structure holds a set of atom species and their positions in fractional coordinates.

Fractional coords of atom `i` is `charges.xf[:, i]`.

**Example use**

```
atoms = Atoms(2, [:C, :F], [0.0 1.0; 2.0 3.0; 4.0 5.0])
```

**Attributes**

  * `n_atoms::Int`: number of atoms
  * `species::Array{Symbol, 1}`: atom species
  * `xf::Array{Float64, 2}`: fractional coordinates in the columns


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/acc15df28b174dd79344a2142c46f7a0e16e948b/src/Matter.jl#L1-L13' class='documenter-source'>source</a><br>

<a id='PorousMaterials.Charges' href='#PorousMaterials.Charges'>#</a>
**`PorousMaterials.Charges`** &mdash; *Type*.



Data structure holds a set of point charges and their positions in fractional coordinates.

Fractional coords of charge `i` is `charges.xf[:, i]`.

**Example use**

```
charges = Charges(2, [-1.0, 1.0], [0.0 1.0; 2.0 3.0; 4.0 5.0])
```

**Attributes**

  * `n_charges::Int`: number of charges
  * `q::Array{Float64, 1}`: signed magnitude of charges (units: electrons)
  * `xf::Array{Float64, 2}`: fractional coordinates in the columns


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/acc15df28b174dd79344a2142c46f7a0e16e948b/src/Matter.jl#L25-L37' class='documenter-source'>source</a><br>

