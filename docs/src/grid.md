```@meta
DocTestSetup = quote
  using PorousMaterials
end
```

# Potential Energy Grid

`PorousMaterials.jl` allows us to calculate and store the ensemble average potential energy of a molecule inside a crystal. This is done by using the molecule as a probe to measure the potential energy on a grid of points superimposed on the unit cell of the crystal.

### Calculating Potential Energy Grids

Superimpose a grid of points about the unit cell of SBMOF-1, compute the potential energy of xenon at each point, and store the data in a [`Grid`](@ref) object using [`energy_grid`](@ref).

```jldoctest grid
xtal = Crystal("SBMOF-1.cif")
strip_numbers_from_atom_labels!(xtal)
molecule = Molecule("Xe")
ljforcefield = LJForceField("UFF")
grid = energy_grid(xtal, molecule, ljforcefield; resolution=0.5, units=:kJ_mol)

# output

Computing energy grid of Xe in SBMOF-1.cif
	Regular grid (in fractional space) of 25 by 13 by 47 points superimposed over the unit cell.
Regular grid of 25 by 13 by 47 points superimposed over a unit cell and associated data.
	units of data attribute: kJ_mol
	origin: [0.000000, 0.000000, 0.000000]
```

The [`Grid`](@ref) object has the following attributes:

```jldoctest grid; output=false
grid.box      # Bravais lattice over which a grid of points is superimposed
grid.data     # 3 dim array containing data for each point
grid.n_pts    # number of grid points in x, y, z
grid.origin   # the origin of the grid
grid.units    # units associated with each data point

# output

:kJ_mol
```

### Saving and Retrieving Grids

Write to a `.cube` volume file to visualize the potential energy contours. The output file location is determined by `rc[:paths][:grids]`.

```julia
write_cube(grid, "CH4_in_SBMOF1.cube")
```

Likewise, we can read a `.cube` file to populate a `Grid` object:

```julia
filename = joinpath(rc[:paths][:grids], "CH4_in_SBMOF1.cube")
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
