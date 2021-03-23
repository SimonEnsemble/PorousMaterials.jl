"""
    res = find_energy_minimum(xtal::Crystal, molecule::Molecule{Frac}, ljff::LJForceField; xf₀::Union{Frac, Nothing}=nothing, n_pts::Tuple{Int, Int, Int}=(50,50,50))
    res = find_energy_minimum(xtal::Crystal, molecule::Molecule{Frac}, ljff::LJForceField, x₀::Cart; n_pts::Tuple{Int, Int, Int}=(50,50,50))

Find the minimum energy of a given `molecule` inside the `xtal`. 
The optimizer needs an initial estimate xf₀. 
If xf₀==nothing (default), automatically generate initial estimate via [`find_energy_minimum_gridsearch`](@ref).
Otherwise, the initial estimate will have to be provided.

# Arguments
- `xtal::Crystal`: The crystal being investigated
- `molecule::Molecule{Frac}`: The molecule used to probe energy surface
- `ljff::LJForceField`: The force field used to calculate interaction energies
- `xf₀::Union{Frac, Nothing}=nothing`: Initial estimate used by the optimizer. If xf₀==nothing, an estimate will be automatically generated.
- `n_pts::Tuple{Int, Int, Int}=(50,50,50)`: Number of grid points in each fractional coordinate dimension, including endpoints (0, 1)

# Returns
- `res::MultivariateOptimizationResults`: contains results of optimization

To view the candidate solution, look at the arrtibutes:
- `res.minimum`: The estimated minimum energy, units: kJ/mol
- `res.minimizer`: The location of `res.minimum` in fractional coordinates
"""
function find_energy_minimum(xtal::Crystal,
                             molecule::Molecule{Frac},
                             ljff::LJForceField;
                             xf₀::Union{Frac, Nothing}=nothing,
                             n_pts::Tuple{Int, Int, Int}=(50,50,50)
                            )
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

    # run the optimizer
    res = optimize(energy, xf₀.xf)

    # check that the solution falls within the xtal.box
    if any(res.minimizer .> 1.0) || any(res.minimizer .< 0.0)
        @warn "minimizer is outside fractional box. Still a valid result, but wrap it to the box."
    end

    return res
end

find_energy_minimum(xtal::Crystal, molecule::Molecule{Frac}, ljff::LJForceField, x₀::Cart; n_pts::Tuple{Int, Int, Int}=(50,50,50)) = find_energy_minimum(xtal, molecule, ljff; xf₀=Frac(x₀, xtal.box), n_pts=n_pts)


"""
    xf₀ = find_energy_minimum_gridsearch(xtal::Crystal,
                                         molecule::Molecule{Cart},
                                         ljff::LJForceField;
                                         n_pts::Tuple{Int, Int, Int}=(50,50,50)
                                        )

Perform an [`energy_grid`](@ref) calculation and return the location of the minimum energy in fractional coordinates.

# Arguments
- `xtal::Crystal`: The crystal being investigated
- `molecule::Molecule{Cart}`: The molecule used to probe energy surface
- `ljff::LJForceField`: The force field used to calculate interaction energies
- `n_pts::Tuple{Int, Int, Int}=(50,50,50)`: Number of grid points in each fractional coordinate dimension, including endpoints (0, 1)

# Returns
- `xf₀::Frac`: location of minimum energy
"""
function find_energy_minimum_gridsearch(xtal::Crystal,
                                        molecule::Molecule{Cart},
                                        ljff::LJForceField;
                                        n_pts::Tuple{Int, Int, Int}=(50,50,50)
                                        )
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
