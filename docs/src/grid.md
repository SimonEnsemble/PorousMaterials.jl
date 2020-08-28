# Potential Energy Grid

`PorousMaterials.jl` allows us to calculate and store the ensemble average potential energy of a molecule inside a crystal. This is done by using the molecule as a probe to measure the potential energy on a grid of points superimposed on the unit cell of the crystal.

### Calculating Potential Energy Grids

Superimpose a grid of points about the unit cell of SBMOF-1, compute the potential energy of xenon at each point, and store the data in a [`Grid`](@ref) object using [`energy_grid`](@ref).

```julia
using PorousMaterials

xtal = Crystal("SBMOF-1.cif")
strip_numbers_from_atom_labels!(xtal)
molecule = Molecule("Xe")
ljforcefield = LJForceField("UFF")

# Grid data structure
grid = energy_grid(xtal, molecule, ljforcefield, 
                   n_pts=(50, 50, 50), units=:kJ_mol) 
```

The [`Grid`](@ref) object has the following attributes:

```julia
grid.box      # Bravais lattice over which a grid of points is superimposed
grid.data     # 3 dim array containing data for each point
grid.n_pts    # number of grid points in x, y, z
grid.origin   # the origin of the grid
grid.units    # units associated with each data point
``` 

### Saving and Retrieving Grids

Write to a `.cube` volume file to visualize the potential energy contours. The output file location is determined by `PorousMaterials.PATH_TO_GRIDS`.

```julia
write_cube(grid, "CH4_in_SBMOF1.cube")
```

Likewise, we can read a `.cube` file to populate a `Grid` object:

```julia
filename = joinpath(PorousMaterials.PATH_TO_GRIDS, "CH4_in_SBMOF1.cube")
grid = read_cube(filename)
```

# detailed docs
## Grids
```@docs
    Grid
    energy_grid
    write_cube
    read_cube
    required_n_pts
    xf_to_id
    id_to_xf
    update_density!
    compute_accessibility_grid
    accessible
```
