# Crystals

`PorousMaterials.jl` maintains a data structure `Crystal` that stores information about a crystal structure file.

## reading in a crystal structure file

read in a crystal as:
```julia
xtal = Crystal("IRMOF-1.cif")
```

## retreiving attributes of a `Crystal`

### unit cell box

`xtal.box`

### atoms

`xtal.atoms`

### charges

`xtal.charges`

## crystal density

calculate the crystal density with the function [`crystal_density`](@ref).

```julia
ρ = crystal_density(xtal) # kg/m³
```

# details

```@docs
    Crystal
    crystal_density
```
