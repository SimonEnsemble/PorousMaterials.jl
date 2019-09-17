## Loading in Molecule Files

Molecule input files are stored in `PorousMaterials.PATH_TO_MOLECULES`. Each molecule possesses its own directory and contains two files: `point_charges.csv` and `lennard_jones_spheres.csv`, comma-separated-value files describing the point charges and Lennard Jones spheres, respectively, comprising the molecule. Only rigid molecules are currently supported. Units of length are in Angstrom; units of charges are electrons.

```julia
using PorousMaterials

m = Molecule("CO2")
```

PorousMaterials will then output information about the molecule you just loaded:

```julia
Molecule species: CO2
Center of mass (fractional coords): [0.0, 0.0, 0.0]
Atoms:

        atom = C_CO2, xf = [0.000, 0.000, 0.000]
        atom = O_CO2, xf = [-1.160, 0.000, 0.000]
        atom = O_CO2, xf = [1.160, 0.000, 0.000]
Point charges:
        charge = 0.700000, xf = [0.000, 0.000, 0.000]
        charge = -0.350000, xf = [-1.160, 0.000, 0.000]
        charge = -0.350000, xf = [1.160, 0.000, 0.000]
```

## Building Blocks of PorousMaterials: Molecules

```julia
molecule = Molecule("CO2") # fractional coords in terms of unit cube box

# access Lennard-Jones spheres & point charges that comprise molecule
molecule.atoms
molecule.charges

# translate to [1.0, 2.0, 3.0] fractional coordinates
translate_to!(molecule, [1.0, 2.0, 3.0])

# translate by [0.1, 0.0, 0.0] fractional coordinates
translate_by!(molecule, [0.1, 0.0, 0.0])

# conduct a uniform random rotation
rotate!(molecule, UnitCube()) # b/c now fractional coords defined in context of a unit cube
```

## Molecules
```@docs
    Molecule
    n_atoms
    translate_to!
    rotate!
    rotation_matrix()
    rotation_matrix(::Float64, ::Array{Float64, 1}; ::Bool)
    rand_point_on_unit_sphere
    charged(::Molecule; ::Bool)
```

## Molecular Movement
```@docs
    insert_molecule!
    delete_molecule!
    translate_molecule!
    reinsert_molecule!
    rotatable
```
