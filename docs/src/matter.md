# Matter and Coordinates

`Atoms` and `Charges` are the building blocks of `Crystal`s and `Molecule`s in `PorousMaterials.jl`. Each have coordinates in both `Cart`esian and `Frac`tional space (associated with unit cell information, i.e., a `Box`).

## Coordinates

we store coordinates as an abstract `Coords` type that has two subtypes: `Cart` and `Frac` for Cartesian and Fractional, respectively.
see the [Wikipedia](https://en.wikipedia.org/wiki/Fractional_coordinates) page on fractional coordinates, which are defined in the context of a periodic system, e.g. within a crystal.

construct coordinates of `n` particles by passing a `n` by `3` array 
```julia
coord = Cart([1.0, 2.0, 5.0])  # construct cartesian coordinate of a particle
coord.x                        # 3 x 1 array, [1, 2, 3]

coord = Frac([0.1, 0.2, 0.5])  # construct fractional coordinate of a particle
coord.xf                       # 3 x 1 array, [0.1, 0.2, 0.3]
```

the coordinates of multiple particles are stored column-wise:
```julia
coords = Cart(rand(3, 5))      # five particles at uniform random coordinates
```

many `Array` operations work on `Coords`, such as:
```julia
coords[2]                      # coordinate of 2nd particle
coords[2:3]                    # (slicing by index) coords of particles 2 and 3
coords[[1, 2, 5]]              # (slicing by index) coords of particles 1, 2, and 5
coords[rand(Bool, 5)]          # (boolean slicing) coords, selected at random
length(coords)                 # number of particles, (5)
```

### manipulating coordinates

`Coords` are immutable:
```julia
coords.x = rand(3, 3)          # fails! coordinates are immutable
```
but we can manipulate the values of `Array{Float64, 2}` where coordinates (through `coords.x` or `coords.xf`) are stored:

```julia
coords.x[2, 3] = 100.0         # successful!
coords.x[:] = rand(3, 3)       # successful! (achieves the above, but need the [:] to say "overwrite all of the elements"
```

fractional coordinates can be wrapped to be inside the unit cell box:
```julia
coords = Frac([1.2, -0.3, 0.9])
wrap!(coords)
coords.xf                      # [0.2, 0.7, 0.9]
```

we can translate coordinates by a vector `dx`:
```julia
dx = Cart([1.0, 2.0, 3.0])
coords = Cart([1.0, 0.0, 0.0])  
translate_by!(coords, dx)
coords.x                        # [2.0, 2.0, 3.0]
```

if `dx::Frac` and `coords::Cart`, `translate_by!` requires a `Box` to convert between fractional and cartesian, as the last argument:
```julia
dx = Frac([0.1, 0.2, 0.3])
box = unit_cube()
coords = Cart([1.0, 0.0, 0.0])
translate_by!(coords, dx)       # fails! need to know Box...
translate_by!(coords, dx, box)
coords.x                        # [1.1, 0.2, 0.3]
```

## Atoms

an atom is specified by its coordinates and atomic species. we can construct a set of atoms (perhaps, comprising a molecule or crystal) as follows.

```julia
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
```

the last line illustrates the two subtypes of `Atoms`, depending on whether the `Coords` are stored as `Frac`tional or `Cart`esian.

we can slice atoms, such as:
```julia
atoms[1]                         # 1st atom
atoms[2:3]                       # 2nd and 3rd atom
```

and combine them:
```julia
atoms_combined = atoms[1] + atoms[2:3]   # combine atoms 1, 2, and 3
isapprox(atoms, atoms_combined)          # true
```

## Charges

`Charges`, well, point charges, work analogously to atoms, except instead of `species`, the values of the point charges are stored in an array, `q`.

```
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
```

we can determine if the set of point charges comprise a charge-neutral system by:
```julia
net_charge(charges)                 # 0.0
neutral(charges)                    # true
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
