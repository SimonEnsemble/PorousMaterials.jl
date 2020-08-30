# Distances

The distance between two `Atoms` in a `Crystal` is central to many operations
within `PorousMaterials.jl`.  The `distance` function calculates the `Cart`esian
displacement between the `Coords` (`Cart` or `Frac`) of two points, `i` and `j`,
within a given `box`.

```julia
xtal = Crystal("Co-MOF-74.cif")
distance(xtal.atoms.coords, xtal.box, 1, 2, false) # 23.2 Å
```

The `apply_pbc` argument allows for calculation of distances
across the periodic boundaries of the `box`.

```julia
distance(xtal.atoms.coords, xtal.box, 1, 2, true) # 3.34 Å
```

`distance` also works on `Atoms` and `Charges`.

```julia
distance(xtal.atoms, xtal.box, 3, 5, true)
```

# docs

```@docs
    distance
```
