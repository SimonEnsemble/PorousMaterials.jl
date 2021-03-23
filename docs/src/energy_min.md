# Find Potential Energy Minimum

`PorousMaterials.jl` allows us to find the potential energy minimum of a probe molecule inside a crystal structure for a given force field. 

### Example
```julia
# Perform an energy grid calculation on a course grid to get initial estimate.
grid = energy_grid(xtal, molecule, ljff; n_pts=(50,50,50))

# Store the minimum value of the grid data.
E_min_estimate = findmin(grid.data)

# Get the Cartesian index of the voxel with the minimum energy estimate.
vox_id = Tuple([E_min_estimate[2][i] for i in 1:3])

# Get the fractional coordinates of the voxel with the minimum energy estimate.
#     These are in fractional coordinates scaled to the size of xtal.box, but still need to be cast as Frac.
xf_minE = id_to_xf(vox_id, grid.n_pts) 

##
# Apply grid search optimization to find the minimum energy starting from input location.
#     If xf₀ is not specified, the previous steps are done automatically.
##
res = find_energy_minimum(xtal, Frac(molecule, xtal.box), ljff; xf₀=Frac(xf_minE))

# The minimum value found by the opimizer
#     units: kJ/mol
energy_minimum = res.minimum

# The location of the minimum value.
#     These are scaled to the size required by the force field cutoff radius
xf_minimum = res.minimizer

# Rescale to input xtal.box (if desired)
xf_minimum_rescale = xf_minimum .* replication_factors(xtal, sqrt(ffield.r²_cutoff))
```

# detailed docs
```@docs
    find_energy_minimum
    find_energy_minimum_gridsearch
```
