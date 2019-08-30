## Loading in Crystal Structure Files

Place `.cif` and `.cssr` crystal structure files in `PorousMaterials.PATH_TO_CRYSTALS`. `PorousMaterials.jl` currently takes crystals in P1 symmetry only. From here you can start julia and do the following to load a framework and start working with it.

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
    chemical_formula
    molecular_weight
    crystal_density
    replicate(::Framework, ::Tuple{Int, Int, Int})
    charged(::Framework; ::Bool)
    write_cif
    assign_charges
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
