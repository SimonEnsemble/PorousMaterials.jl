"""
    minimized_molecule, min_energy  = find_energy_minimum(xtal, molecule, ljff) # molecule set at initial guess

find the minimum energy position, and associated minimum energy, of a molecule in a crystal.
n.b. currently works only for molecules with one atom.
the optimizer needs an initial estimate of the minimum energy position. 
pass molecule with good initial position.
if you don't have a good initial position, use [`find_energy_minimum_gridsearch`](@ref).

# Arguments
- `xtal::Crystal`: the crystal
- `molecule::Molecule`: the molecule, whose position we seek to tune until we reach a local minimum. must start at a good initial position close to the minimum.
- `ljff::LJForceField`: the force field used to calculate crystal-molecule interaction energies

# Returns
- `minimized_molecule::Molecule{Frac}`: the molecule at its minimum energy position
- `min_energy::Float64`: the associated minimum molecule-crystal interaciton energy (kJ/mol)
"""
function find_energy_minimum(xtal::Crystal,
                             molecule::Molecule,
                             ljff::LJForceField,
                             )
    if needs_rotations(molecule)
        @warn "needs rotations. does not optimize over configurations, only over center of mass"
    end
    
    # make sure replication factors sufficient
    rep_factors = replication_factors(xtal, sqrt(ljff.r²_cutoff))
    xtal = replicate(xtal, rep_factors)
    
    # make sure molecule is in fractional coords
    if isa(molecule, Molecule{Cart})
        molecule = Frac(molecule, xtal.box)
    else
        translate_to!(molecule, Frac(molecule.com.xf ./ rep_factors))
    end
    
    # xf::Array{Float64, 1}
    function energy(xf)
        # make sure the coords are fractional
        xf = mod.(xf, 1.0)
        # move probe molecule
        translate_to!(molecule, Frac(xf))
        # calculate the guest-host VDW interaction and return
        return vdw_energy(xtal, molecule, ljff) * 8.314 / 1000 # units: kJ/mol
    end

    # find minimum energy position
    res = optimize(energy, molecule.com.xf)
    xf_min = mod.(mod.(res.minimizer, 1.0) .* rep_factors, 1.0)
    # translate molecule to min energy position to return it.
    translate_to!(molecule, Frac(xf_min))
    
    return molecule, res.minimum
end

"""
    xf₀ = find_energy_minimum_gridsearch(xtal, molecule, ljff; n_pts=(50, 50, 50))

perform an [`energy_grid`](@ref) calculation and, via a grid search, find the minimum energy position of a molecule.

# Arguments
- `xtal::Crystal`: The crystal being investigated
- `molecule::Molecule{Cart}`: The molecule used to probe energy surface
- `ljff::LJForceField`: The force field used to calculate interaction energies
- `n_pts::Tuple{Int, Int, Int}=(50,50,50)`: Number of grid points in each fractional coordinate dimension, including endpoints (0, 1)

# Returns
- `minimized_molecule::Molecule{Frac}`: the molecule at its minimum energy position
- `min_energy::Float64`: the associated minimum molecule-crystal interaciton energy (kJ/mol)
"""
function find_energy_minimum_gridsearch(xtal::Crystal,
                                        molecule::Molecule{Cart},
                                        ljff::LJForceField;
                                        n_pts::Tuple{Int, Int, Int}=(50, 50, 50)
                                        )
    # Perform an energy grid calculation on a course grid to get initial estimate.
    grid = energy_grid(xtal, molecule, ljff; n_pts=n_pts)
    E_min, idx_min = findmin(grid.data)

    xf_minE = id_to_xf(Tuple(idx_min), grid.n_pts) # fractional coords of min
    molecule = Frac(molecule, xtal.box)
    translate_to!(molecule, Frac(xf_minE))
    return molecule, E_min
end
