```@meta
DocTestSetup = quote
  using PorousMaterials
end
```

# Matter and Coordinates

`Atoms` and `Charges` are the building blocks of `Crystal`s and `Molecule`s in `PorousMaterials.jl`. Each have coordinates in both `Cart`esian and `Frac`tional space (associated with unit cell information, i.e., a `Box`).

## Coordinates

we store coordinates as an abstract `Coords` type that has two subtypes: `Cart` and `Frac` for Cartesian and Fractional, respectively.
see the [Wikipedia](https://en.wikipedia.org/wiki/Fractional_coordinates) page on fractional coordinates, which are defined in the context of a periodic system, e.g. within a crystal.

construct coordinates of `n` particles by passing a `n` by `3` array 

```jldoctest; output=false
coord = Cart([1.0, 2.0, 5.0])  # construct cartesian coordinate of a particle
coord.x                        # 3 x 1 array, [1, 2, 3]

coord = Frac([0.1, 0.2, 0.5])  # construct fractional coordinate of a particle
coord.xf                       # 3 x 1 array, [0.1, 0.2, 0.3]
# output
3×1 Matrix{Float64}:
 0.1
 0.2
 0.5
```

the coordinates of multiple particles are stored column-wise:

```jldoctest matter; output=false
coords = Cart([
  0.651027   0.274176   0.386178  0.371651  0.619131;
  0.0681196  0.0267313  0.836004  0.819681  0.585807;
  0.667704   0.825569   0.780142  0.606194  0.572355
])      # five particles at uniform random coordinates
# output
Cart([0.651027 0.274176 … 0.371651 0.619131; 0.0681196 0.0267313 … 0.819681 0.585807; 0.667704 0.825569 … 0.606194 0.572355])
```

many `Array` operations work on `Coords`, such as:

```jldoctest matter; output=false
coords[2]                      # coordinate of 2nd particle
coords[2:3]                    # (slicing by index) coords of particles 2 and 3
coords[[1, 2, 5]]              # (slicing by index) coords of particles 1, 2, and 5
coords[rand(Bool, 5)]          # (boolean slicing) coords, selected at random
length(coords)                 # number of particles, (5)
# output
5
```

### manipulating coordinates

`Coords` are immutable:

```jldoctest matter
coords.x = rand(3, 5)
# output
ERROR: setfield!: immutable struct of type Cart cannot be changed
```
but we can manipulate the values of `Array{Float64, 2}` where coordinates (through `coords.x` or `coords.xf`) are stored:

```jldoctest matter; output=false
coords.x[2, 3] = 100.0
coords.x[:] = [ # need the [:] to say "overwrite all of the elements"
  0.496997  0.560528   0.496615  0.213062  0.21751;
  0.772372  0.697443   0.133055  0.211073  0.02676;
  0.230555  0.0988727  0.592699  0.193649  0.16536
]
# output
3×5 Matrix{Float64}:
 0.496997  0.560528   0.496615  0.213062  0.21751
 0.772372  0.697443   0.133055  0.211073  0.02676
 0.230555  0.0988727  0.592699  0.193649  0.16536
```

fractional coordinates can be wrapped to be inside the unit cell box:

```jldoctest
coords = Frac([1.2, -0.3, 0.9])
wrap!(coords)
# output
3×1 Matrix{Float64}:
 0.19999999999999996
 0.7
 0.9
```

we can translate coordinates by a vector `dx`:

```jldoctest
dx = Cart([1.0, 2.0, 3.0])
coords = Cart([1.0, 0.0, 0.0])  
translate_by!(coords, dx)
# output
3×1 Matrix{Float64}:
 2.0
 2.0
 3.0
```

if `dx::Frac` and `coords::Cart`, `translate_by!` requires a `Box` to convert between fractional and cartesian, as the last argument:

```jldoctest
dx = Frac([0.1, 0.2, 0.3])
box = unit_cube()
coords = Cart([1.0, 0.0, 0.0])
translate_by!(coords, dx, box)
# output
3×1 Matrix{Float64}:
 1.1
 0.20000000000000004
 0.3
```

## Atoms

an atom is specified by its coordinates and atomic species. we can construct a set of atoms (perhaps, comprising a molecule or crystal) as follows.

```jldoctest matter
species = [:O, :H, :H]            # atomic species are represnted with Symbols
coords = Cart([0.0 0.757 -0.757;  # coordinates of each
               0.0 0.586  0.586; 
               0.0 0.0    0.0   ]
             )
atoms = Atoms(species, coords)    # 3 atoms comprising water
atoms.n                           # number of atoms, 3
atoms.coords                      # coordinates; atoms.coords.x gives the array of coords
atoms.species                     # array of species
atoms::Atoms{Cart}                # successful type assertion, as opposed to atoms::Atoms{Frac}
# output
Atoms{Cart}(3, [:O, :H, :H], Cart([0.0 0.757 -0.757; 0.0 0.586 0.586; 0.0 0.0 0.0]))
```

the last line illustrates the two subtypes of `Atoms`, depending on whether the `Coords` are stored as `Frac`tional or `Cart`esian.

we can slice atoms, such as:

```jldoctest matter
atoms[2:3]
# output
Atoms{Cart}(2, [:H, :H], Cart([0.757 -0.757; 0.586 0.586; 0.0 0.0]))
```

and combine them:

```jldoctest matter
atoms_combined = atoms[1] + atoms[2:3]
isapprox(atoms, atoms_combined)
# output
true
```

## Charges

`Charges`, well, point charges, work analogously to atoms, except instead of `species`, the values of the point charges are stored in an array, `q`.

```jldoctest matter
q = [-1.0, 0.5, 0.5]              # values of point charges, units: electrons
coords = Cart([0.0 0.757 -0.757;  # coordinates of the point charges
               0.0 0.586  0.586; 
               0.0 0.0    0.0   ]
             )
charges = Charges(q, coords)      # 3 point charges
charges.n                         # number of charges, 3
charges.coords                    # retreive coords
charges.q                         # retreive q
charges::Charges{Cart}            # successful type assertion, as opposed to charges::Charges{Frac}
# output
Charges{Cart}(3, [-1.0, 0.5, 0.5], Cart([0.0 0.757 -0.757; 0.0 0.586 0.586; 0.0 0.0 0.0]))
```

we can determine if the set of point charges comprise a charge-neutral system by:
```jldoctest matter; output=false
net_charge(charges)                 # 0.0
neutral(charges)                    # true
# output
true
```

# detailed docs

```@docs
    Coords
    Frac
    Cart
    Atoms
    Charges
    net_charge
    neutral
    translate_by!
```
