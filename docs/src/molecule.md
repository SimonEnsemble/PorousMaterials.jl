```@meta
DocTestSetup = quote
  using PorousMaterials
end
```

# Molecules

## Loading Molecule Files

Molecule input files are stored in the path indicated by the `Xtals`/`PorousMaterials` global dictionary `rc`. 

```julia
rc[:paths][:molecules] # absolute or relative path
```

Each molecule possesses its own directory containing two files: `charges.csv` and `atoms.csv`, comma-separated-value files, which describe the point charges and Lennard Jones spheres, respectively, that compose the molecule. Only rigid molecules are currently supported. Units of length are in Angstroms ($\AA$); units of charges are electrons.

```jldoctest molecule
molecule = Molecule("CO2")
# output
Molecule species: CO2
Center of mass (fractional coords): Cart([0.0; 0.0; 0.0;;])
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

```jldoctest molecule; output=false
# access the attributes that comprise the molecule object
molecule.species   # molecule species
molecule.com       # center-of-mass
molecule.atoms     # Lennard-Jones spheres
molecule.charges   # point charges 
# output
Charges{Cart}(3, [0.7, -0.35, -0.35], Cart([0.0 -1.16 1.16; 0.0 0.0 0.0; 0.0 0.0 0.0]))
```

To see specific information about the atoms and charges attributes of the molecule see [`Atoms`](@ref) and [`Charges`](@ref).

## Pseudo-Atoms

Sometimes, [e.g. modeling quadrupolar molecules](https://github.com/SimonEnsemble/PorousMaterials.jl/issues/236), it is desirable to add a point charge in a location outside of any atomic nucleus.
In order to do this, a pseudo-atom label must be chosen (consistent with input data) and added to the global list of atomic masses.

Take the example of using [TraPPE](https://doi.org/10.1002/aic.690470719)-derived charges for the dinitrogen molecule.
Partial negative charges are assigned to the nitrogen atoms and a corresponding positive charge is placed at a pseudo-atom "dummy site" at the center of mass.
The atoms and charges input data for this "N2_TraPPE" molecule are like so:

```jldoctest molecule
dir = joinpath(rc[:paths][:molecules], "N2_TraPPE")
for file in readdir(dir)
  @info file data=String(read(joinpath(dir, file)))
end
# output
┌ Info: atoms.csv
└   data = "atom,x,y,z\r\nN_in_N2,0,0,0.55\r\nPSEUDOATOM_LABEL,0,0,0\r\nN_in_N2,0,0,-0.55"
┌ Info: charges.csv
└   data = "q,x,y,z\r\n-0.482,0,0,0.55\r\n0.964,0,0,0\r\n-0.482,0,0,-0.55"
```

**atoms.csv**

```
atom,x,y,z
N_in_N2,0,0,0.55
PSEUDOATOM_LABEL,0,0,0
N_in_N2,0,0,-0.55
```

**charges.csv**

```
q,x,y,z
-0.482,0,0,0.55
0.964,0,0,0
-0.482,0,0,-0.55
```

In this example, `PSEUDOATOM_LABEL` is the label for the pseudo-atom; there are several common choices, so make sure you know which one is used in your particular data!

Before loading these data by calling `Molecule`, we need to add our (massless) pseudo-atom to the mass dictionary:

```juliadoctest
rc[:atomic_masses][:PSEUDOATOM_LABEL] = 0.
Molecule("N2_TraPPE")
# ouput
asdf
```

## Moving Molecules
We can translate and roatate a molecule:

```julia
# convert to Molecule{Frac}
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
