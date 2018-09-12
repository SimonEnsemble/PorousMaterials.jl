# Matter

Matter is the foundation for every data type defined in PorousMaterials. Two structs are used to define all atomic structures: `Atoms` and `Charges`. Every atom, molecule, or crystal structure can be simulated by understanding its atoms and its point charges.

We chose to store each collection as a single object (Atoms) rather than an array of objects (Array{Atom}) because it allows us to store the locations in contiguous memory. We found significant speed increases when storing each collection as a single object.

One array of locations also allowed us to take advantage of broadcasting. For example, it is useful when calculating the energy inside a framework. Instead of looping through every atom, we can instead run this line:

```
dxf = broadcast(-, framework.atoms.xf, molecule.atoms.xf[i])
```

This calculates the distance between one atom in a molecule and every atom in the framework.

## Atoms

The `Atoms` struct has three attributes: `species`, `xf`, and `n_atoms`. The `species` field is a 1D array of `Symbol`s that stores the species of the different atoms in this collection. The `xf` field is a 2D array of `Float64`s that stores the location of atoms in 3D space. Each column represents the x, y, and z coordinates of an atom in this collection. The `n_atoms` field is an `Int` that is used for looping through the atoms deifned by this object. It does not need to be passed in when creating and `Atoms` object, instead it is calculated.

The `xf` and `species` arrays line up, so each species corresponds to a column in `xf` and thus has a location in 3D space.

## Charges

The `Charges` struct is similar to the `Atoms` struct because it also has three attributes: `q`, `xf`, and `n_charges`. `q` is a 1D array of `Float64`s that stores the charge values for different point charges. The `xf` field is the same as the one in `Atoms`, it is a 2D array that stores the location of the charges in 3D space. `n_charges` is an `Int` that is used for looping through the point charges defined by this object. As with `Atoms` the `n_charges` value doesn't need to be passed in, it will be automatically calculated by a constructor.

The `xf` and `q` arrays line up, so each charge corresponds to a column in `xf`, and thus has a location in 3D space.
