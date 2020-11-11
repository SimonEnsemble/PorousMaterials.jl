# Molecules

## Loading Molecule Files

Molecule input files are stored in `PorousMaterials.PATH_TO_MOLECULES`. Each molecule possesses its own directory containing two files: `charges.csv` and `atoms.csv`, comma-separated-value files, which describe the point charges and Lennard Jones spheres, respectively, that compose the molecule. Only rigid molecules are currently supported. Units of length are in Angstroms ($\AA$); units of charges are electrons.

```julia
using PorousMaterials

molecule = Molecule("CO2")
```

PorousMaterials will then output information about the molecule you just loaded:

```julia
Molecule species: CO2
Center of mass (fractional coords): Cart([0.0; 0.0; 0.0])
Atoms:

	atom = C_CO2, x = [0.000, 0.000, 0.000]
	atom = O_CO2, x = [-1.160, 0.000, 0.000]
	atom = O_CO2, x = [1.160, 0.000, 0.000]
Point charges: 
	charge = 0.700000, x = [0.000, 0.000, 0.000]
	charge = -0.350000, x = [-1.160, 0.000, 0.000]
	charge = -0.350000, x = [1.160, 0.000, 0.000]
```

## Building Blocks of PorousMaterials: Molecules

```julia
molecule = Molecule("CO2") # fractional coords in terms of unit cube box

# access the attributes that comprise the molecule object
molecule.species   # molecule species
molecule.com       # center-of-mass
molecule.atoms     # Lennard-Jones spheres
molecule.charges   # point charges 
```

To see specific information about the atoms and charges attributes of the molecule see [`Atoms`](@ref) and [`Charges`](@ref).

## Moving Molecules
We can translate and roatate a molecule:

```julia
using PorousMaterials

# load a molecule and convert it to Molecule{Frac}
molecule = Molecule("CO2")
molecule = Frac(molecule, unit_cube())

# translate center-of-mass to [1.0, 2.0, 3.0] in fractional coordinates
x = Cart([1.0, 2.0, 3.0])
translate_to!(molecule, x, unit_cube())

# translate by [0.1, 0.0, 0.0] in fractional coordinates
dx = Cart([0.1, 0.0, 0.0])
translate_by!(molecule, dx, unit_cube())

# conduct a uniform random rotation about the center-of-mass
random_rotation!(molecule, unit_cube()) 
```

# detailed docs

## Molecules
```@docs
    Molecule
    translate_to!
    random_rotation!
    random_rotation_matrix()
    ion
    distortion
```

## Molecular Movement
```@docs
    apply_periodic_boundary_condition!
    random_insertion!
    remove_molecule!
    random_translation!
    random_reinsertion!
    needs_rotations
```
