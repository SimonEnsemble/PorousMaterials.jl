Matter is the foundation for every data type defined in PorousMaterials. Two structs are used to define all atomic structures: `Atoms` and `Charges`. Every atom, molecule, or crystal structure can be simulated by understanding its atoms and its point charges.

We chose to store each collection as a single object (Atoms) rather than an array of objects (Array{Atom}) because it allows us to store the locations in contiguous memory. We found significant speed increases when storing each collection as a single object.

One array of locations also allowed us to take advantage of broadcasting. For example, it is useful when calculating the energy inside a framework. Instead of looping through every atom, we can instead run this line:

```
dxf = broadcast(-, framework.atoms.xf, molecule.atoms.xf[i])
```

This calculates the distance between one atom in a molecule and every atom in the framework.

## Building Blocks of PorousMaterials: Matter

In `PorousMaterials.jl`, crystals and molecules are composed of `Atoms` and `Charges`

To create a carbon atom at `[0.1, 0.2, 0.5]` fractional coordinates (in the context of some Bravais lattice):
```julia
xf = Array{Float64, 2}(undef, 3, 0)
xf = [xf [0.1, 0.2, 0.5]]
atoms = Atoms([:C], xf) # constructor
atoms.species[1] # :C
atoms.xf[:, 1] # [0.1, 0.2, 0.5]
```

To create a point charge of +1 at `[0.1, 0.2, 0.5]` fractional coordinates (in the context of some Bravais lattice):
```julia
xf = Array{Float64, 2}(undef, 3, 0)
xf = [xf [0.1, 0.2, 0.5]]
charges = Charges([1.0], xf)
charges.q[1] # 1.0
charges.xf[:, 1] # [0.1, 0.2, 0.5]
```

## Matter
```@docs
    Atoms
    Charges
```
