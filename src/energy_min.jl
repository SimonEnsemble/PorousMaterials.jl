"""
    res = find_energy_minimum(xtal::Crystal, molecule::Molecule, ljff::LJForceField; xf₀::Union{Frac, Nothing}=nothing, temperature::Float64=NaN)
    res = find_energy_minimum(xtal::Crystal, molecule::Molecule, ljff::LJForceField, x₀::Cart; temperature::Float64=NaN)

Find the minimum energy of a given `molecule` inside the `xtal` using a grid search optimization routine. 
The optimizer needs an initial estimate xf₀. 
One way to obtain this estimate is by running an [`energy_grid`](@ref) calculation (see example). 
If xf₀==nothing (default), the `energy_grid` calculation is automatically performed with the default function values.  

# Arguments
- `xtal::Crystal`: The crystal being investigated
- `molecule::Molecule`: The Molecule used to probe energy surface
- `ljff::LJForceField`: The force field used to calculate interaction energies
- `xf₀::Union{Frac, Nothing}=nothing`: Initial estimate used by the optimizer. 
                                       If xf₀==nothing, an estimate will be automatically generated.
- `temperature::Float64`: The temperature at which to compute the free energy for molecules where rotations are required. 
                          Only used if xf₀==nothing.


# Returns
- `res::MultivariateOptimizationResults`: contains results of optimization

To view the candidate solution, look at the arrtibutes:
- `res.minimum`: The estimated minimum energy
- `res.minimizer`: The location of `res.minimum` in fractional coordinates

# Example
```
# Perform an energy grid calculation on a course grid to get initial estimate.
grid = energy_grid(xtal, molecule, ljff; n_pts=(50, 50, 50), temperature=298.0, n_rotations=750)

# Store the minimum value of the grid data to be unpacked.
E_min_estimate = findmin(grid.data)

# Get the Cartesian index of the voxel with the minimum energy estimate.
vox_id = Tuple([E_min_estimate[2][i] for i in 1:3])

# Get the fractional coordinates of the voxel with the minimum energy estimate.
#     These are in fractional coordinates scaled to the size of xtal.box, but still need to be cast as Frac.
xf_minE = id_to_xf(vox_id, grid.n_pts) # fractional coords of min

# Apply grid search optimization to find the minimum energy starting from input location.
res = find_energy_minimum(xtal, Frac(molecule, xtal.box), ffield, Frac(xf_minE))

# The minimum value found by the opimizer
#     units: kJ/mol
energy_minimum = res.minimum

# The location of the minimum value.
#     These are scaled to the size required by the forcefield cutoff radius
xf_minimum = res.minimizer

# Rescale to input xtal.box (if desired)
xf_minimum_rescale = xf_minimum .* replication_factors(xtal, sqrt(ffield.r²_cutoff))
```
"""
function find_energy_minimum(xtal::Crystal, molecule::Molecule, ljff::LJForceField; 
                             xf₀::Union{Frac, Nothing}=nothing, n_pts::Tuple{Int, Int, Int}=(50,50,50))
    @assert ! needs_rotations(molecule) "can't handle molecules that need rotations"

    if isnothing(xf₀)
        xf₀ = find_energy_minimum_gridsearch(xtal, Cart(molecule, xtal.box), ljff, n_pts=n_pts)
    end

    # get replication factors
    rep_factors = replication_factors(xtal, sqrt(ljff.r²_cutoff))
    # replicate xtal to meet requirements for VDW interactions
    xtal = replicate(xtal, rep_factors)

    # scale fractional coords to new box size
    xf₀.xf .= xf₀.xf ./ rep_factors

    # define input function for optimizer
    function energy(xf)
        # make sure the coords are fractional
        xf = mod.(xf, 1.0)
        # move probe molecule
        translate_to!(molecule, Frac(xf))
        # calculate the guest-host VDW interaction and return
        return vdw_energy(xtal, molecule, ljff) * 8.314 / 1000 # units: kJ/mol
    end

    opt_res = optimize(energy, xf₀.xf)
    if any(opt_res.minimizer .> 1.0) || any(opt_res.minimizer .< 0.0)
        @warn "minimizer is outside fractional box. Still a valid result, but wrap it to the box."
    end

    return opt_res
end

find_energy_minimum(xtal::Crystal, molecule::Molecule, ljff::LJForceField, x₀::Cart; n_pts::Tuple{Int, Int, Int}=(50,50,50)) = find_energy_minimum(xtal, molecule, ljff; xf₀=Frac(x₀, xtal.box); n_pts=n_pts)

# automate find_energy_minimum  calculation if xf₀==nothing
function find_energy_minimum_gridsearch(xtal::Crystal, molecule::Molecule, ljff::LJForceField; n_pts::Tuple{Int, Int, Int}=(50,50,50))
    # Perform an energy grid calculation on a course grid to get initial estimate.
    grid = energy_grid(xtal, molecule, ljff; n_pts=n_pts)

    # Store the minimum value of the grid data to be unpacked.
    E_min_estimate = findmin(grid.data)

    # Get the Cartesian index of the voxel with the minimum energy estimate.
    vox_id = Tuple([E_min_estimate[2][i] for i in 1:3])

    # Get the fractional coordinates of the voxel with the minimum energy estimate.
    #     These are in fractional coordinates scaled to the size of xtal.box, but still need to be cast as Frac.
    xf_minE = id_to_xf(vox_id, grid.n_pts) # fractional coords of min
    return Frac(xf_minE)
end

