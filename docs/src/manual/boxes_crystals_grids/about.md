# Boxes, Crystals, and Grids

## Loading in Crystal Structure Files

Place `.cif` and `.cssr` crystal structure files in `data/crystals`. `PorousMaterials.jl` currently takes crystals in P1 symmetry only. From here you can start julia and do the following to load a framework and start working with it.

```
julia> using PorousMaterials

julia> f = Framework("SBMOF-1.cif")
Name: SBMOF-1.cif
Bravais unit cell of a crystal.
        Unit cell angles α = 90.000000 deg. β = 100.897000 deg. γ = 90.000000 deg.
        Unit cell dimensions a = 11.619300 Å. b = 5.566700 Å, c = 22.931200 Å
        Volume of unit cell: 1456.472102 Å³

Number of atoms = 120
Number of charges = 0
Chemical formula: Dict(:H=>8,:S=>1,:Ca=>1,:O=>6,:C=>14)
```
