## Loading in Crystal Structure Files

Place `.cif` and `.cssr` crystal structure files in `PorousMaterials.PATH_To_CRYSTALS`. `PorousMaterials.jl` can load in `.cif` files of any symmetry as long as the symmetry operations are included. From here you can start julia and do the following to load a framework and start working with it.

```julia
using PorousMaterials

f = Framework("SBMOF-1.cif")
```

PorousMaterials will then output information about the framework you just loaded:

```julia
Name: SBMOF-1.cif
Bravais unit cell of a crystal.
        Unit cell angles α = 90.000000 deg. β = 100.897000 deg. γ = 90.000000 deg.
        Unit cell dimensions a = 11.619300 Å. b = 5.566700 Å, c = 22.931200 Å
        Volume of unit cell: 1456.472102 Å³

Number of atoms = 120
Number of charges = 0
Chemical formula: Dict(:H=>8,:S=>1,:Ca=>1,:O=>6,:C=>14)
```

If the file is not in P1 symmetry, it will be converted within the framework reader and this message will be displayed.

```julia
┌ Warning: Name_of_file.cif is not in P1 symmetry. It is being converted to P1 for use in PorousMaterials.jl.
└ @ PorousMaterials ~/osu_undergrad/simon_ensemble/PorousMaterials.jl/src/Crystal.jl:284
```

PorousMaterials also gives the option to read in structures in the lower level
symmetries and convert them to P1 before simulation.

```julia
framework = Framework("ORIVOC_clean.cif"; remove_overlap=true, convert_to_p1=false)
# the remove_overlap argument is specific to this structure, not all frameworks need it

###
    Perform any operation on the structure while it is not in P1
###

framework = apply_symmetry_rules(framework)
# the framework is now in P1 and it can be used in simulations
```

## Generating Bond Information for Frameworks

The bonds are stored in a `SimpleGraph` from the `LightGraphs.jl` package, and
can be accessed through the `bonds` attribute.  

### Reading from a file

`PorousMaterials` can read in bonds from `.cif` files if they have the tags
`_geom_bond_atom_site_label_1` and `_geom_bond_atom_site_label_2`. To choose to
read bonds from a file, pass `read_bonds_from_file=true` to the `Framework`
constructor.

```julia
using PorousMaterials

f = Framework("KAXQIL_clean.cif"; read_bonds_from_file=true, convert_to_p1=false)

f.bonds
```

This example uses a structure that is not in P1 symmetry. `PorousMaterials`
cannot replicate a structure or apply symmetry rules if it currently has bonds.
However, this structure can be converted to P1 without bonds, and then bonds can
be inferred for the full P1 structure.

### Inferring bonds using `BondingRule`s

`PorousMaterials` can infer bonds for a structure and populate the bond graph by
using `BondingRule`s. Each `BondingRule` has two species of atoms that it works
for. It also has a minimum and maximum distance that a bond can be defined for
the two atoms.

```julia
using PorousMaterials

f = Framework("SBMOF-1.cif")

# define an array of BondingRule's that will be used to define bonds in the
#   framework. These need to be in the order that they are applied
bonding_rules = [BondingRule(:H, :*, 0.4, 1.2),
                 BondingRule(:*, :*, 0.4, 1.9)]

# infer the bonds for the framework f
infer_bonds!(f, bonding_rules)

# redefine bonding_rules to account for edge cases between Ca and O atoms
bonding_rules = [BondingRule(:H, :*, 0.4, 1.2),
                 BondingRule(:Ca, :O, 0.4, 2.5),
                 BondingRule(:*, :*, 0.4, 1.9)]

# remove old bonds from framework before inferring bonds with new rules
remove_bonds!(f)

# re-infer bonds
infer_bonds!(f, bonding_rules)

# output the bond information to visualize it and double check
write_bond_information(f, "SBMOF-1_bonds.vtk")
```


## Building Blocks of PorousMaterials: Bravais lattice

We later apply periodic boundary conditions to mimic a crystal of infinite extent. A `Box` describes a [Bravais lattice](https://en.wikipedia.org/wiki/Bravais_lattice).

To make a 10 by 10 by 10 Å Bravais lattice with right angles:
```julia
box = Box(10.0, 10.0, 10.0, π/2, π/2, π/2)

box.a, box.b, box.c # unit cell dimensions (10.0 Å)
box.α, box.β, box.γ # unit cell angles (1.57... radians)
box.Ω # volume (1000.0 Å³)
box.f_to_c # fractional to Cartesian coordinate transformation matrix
box.c_to_f # Cartesian to fractional coordinate transformation matrix
box.reciprocal_lattice # rows are reciprocal lattice vectors
```

Replicate a box as follows:
```julia
box = replicate(box, (2, 2, 2)) # new box replicated 2 by 2 by 2
box.a # 20 Å
```

## Building Blocks of PorousMaterials: Porous Crystals

```julia
using PorousMaterials

# read in xtal structure file
framework = Framework("SBMOF-1.cif")

# access unit cell box
framework.box

# access Lennard-Jones spheres and point charges comprising the crystal
framework.atoms
framework.charges

# remove annoying numbers on the atom labels
strip_numbers_from_atom_labels!(framework)

# compute crystal density
ρ = crystal_density(framework) # kg/m3

# compute the chemical formula
cf = chemical_formula(framework)

# assign charges according to atom type
charges = Dict(:Ca => 3.0, :O => 2.0, :C => -1.0, :S => 7.0, :H => -1.0)
charged_framework = assign_charges(framework, charges)

# replicate & visualize
framework = replicate(framework, (3, 3, 3))
write_to_xyz(framework, "SBMOF-1.xyz")
```

## Demo of Potential Energy Grid

Superimpose a grid of points about the unit cell of SBMOF-1. Compute the potential energy of xenon at each point and store as a grid.

```julia
using PorousMaterials

framework = Framework("SBMOF-1.cif")
molecule = Molecule("Xe")
forcefield = LJForceField("UFF.csv")

grid = energy_grid(framework, molecule, forcefield,
    n_pts=(50, 50, 50), units=:kJ_mol) # Grid data structure
```

Write to a .cube volume file to visualize the potential energy contours.
```julia
write_cube(grid, "CH4_in_SBMOF1.cube")
```

## Boxes
```@docs
    Box
    replicate
    UnitCube
    write_vtk
    inside
```

## Crystals
```@docs
    Framework
    remove_overlapping_atoms_and_charges
    strip_numbers_from_atom_labels!
    wrap_atoms_to_unit_cell
    chemical_formula
    molecular_weight
    crystal_density
    replicate(::Framework, ::Tuple{Int, Int, Int})
    charged(::Framework; ::Bool)
    assign_charges
    apply_symmetry_rules
    is_symmetry_equal
    write_cif
    BondingRule
    write_bond_information
    infer_bonds!
    remove_bonds!
    compare_bonds_in_framework
```

## Grids
```@docs
    Grid
    xf_to_id
    update_density!
    apply_periodic_boundary_condition!
    write_cube
    read_cube
    energy_grid
```
