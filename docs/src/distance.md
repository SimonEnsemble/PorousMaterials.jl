```@meta
DocTestSetup = quote
  using PorousMaterials
end
```

# Distances

The distance between two `Atoms` in a `Crystal` is central to many operations
within `PorousMaterials.jl`.  The `distance` function calculates the `Cart`esian
displacement between the `Coords` (`Cart` or `Frac`) of two points, `i` and `j`,
within a given `box`.

```jldoctest distance
xtal = Crystal("SBMOF-1.cif")
distance(xtal.atoms.coords, xtal.box, 1, 10, false) # Cartesian distance within the unit cell
# output
4.962373067546231
```

The `apply_pbc` argument allows for calculation of distances
across the periodic boundaries of the `box`.

```jldoctest distance
distance(xtal.atoms.coords, xtal.box, 1, 10, true) # Cartesian distance accounting for periodic boundary
# output
4.143597209982431
```

`distance` also works on `Atoms` and `Charges`.

```jldoctest distance
distance(xtal.atoms, xtal.box, 3, 5, true)
# output
10.244292605252747
```

# docs

```@docs
    distance
```
