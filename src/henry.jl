"""
```
result = henry_coefficient(crystal, molecule, temperature, ljforcefield,
                            nb_insertions=1e6, verbose=true, ewald_precision=1e-6,
                            autosave=true)
```

Conduct particle insertions to compute the Henry coefficient Kₕ of a molecule in a crystal.
Also, for free, the heat of adsorption and ensemble average energy of adsorption is computed.
The Henry coefficient is a model for adsorption at infinite dilution (low coverage):
⟨N⟩ = Kₕ P, where P is pressure and Kₕ is the Henry coefficient.

Kₕ = β ⟨e^{-β U}⟩, where the average is over positions and orientations of the molecule
in the crystal.

# Arguments
- `crystal::Crystal`: the porous crystal in which we seek to simulate adsorption
- `molecule::Molecule`: the adsorbate molecule
- `temperature::Float64`: temperature of bulk gas phase in equilibrium with adsorbed phase
    in the porous material. units: Kelvin (K)
- `ljforcefield::LJForceField`: the molecular model used to describe the
    energetics of the adsorbate-adsorbate and adsorbate-host van der Waals interactions.
- `insertions_per_volume::Int`: number of Widom insertions to perform for computing the
average, per unit cell volume (Å³)
- `verbose::Bool`: whether or not to print off information during the simulation.
- `ewald_precision::Float64`: desired precision for Ewald summations; used to determine
the replication factors in reciprocal space.
- `autosave::Bool`: save results file as a .jld in rc[:paths][:data] * `sims`
- `filename_comment::AbstractString`: An optional comment that will be appended to the name of the saved file.

# Returns
- `result::Dict{String, Float64}`: A dictionary containing all the results from the Henry coefficient simulation
"""
function henry_coefficient(crystal::Crystal, molecule::Molecule, temperature::Float64,
                           ljforcefield::LJForceField; insertions_per_volume::Union{Int, Float64}=200,
                           verbose::Bool=true, ewald_precision::Float64=1e-6,
                           autosave::Bool=true, filename_comment::AbstractString="",
                           accessibility_grid::Union{Nothing, Grid{Bool}}=nothing)

    # simulation only works if crystal is in P1
    assert_P1_symmetry(crystal)

    time_start = time()
    if verbose
        print("Simulating Henry coefficient of ")
        printstyled(molecule.species; color=:green)
        print(" in ")
        printstyled(crystal.name; color=:green)
        print(" at ")
        printstyled(temperature; color=:green)
        print(" K, using ")
        printstyled(ljforcefield.name; color=:green)
        print(" force field with ")
        printstyled(insertions_per_volume; color=:green)
        println(" insertions per Å³.")

        if !isnothing(accessibility_grid)
            println("Using provided accessibility grid to block inaccessible pockets")
            if ! isapprox(accessibility_grid.box, crystal.box)
                error(@sprintf("accessibility grid box does not match box of %s.\n",
                    crystal.name))
            end
        end
    end

    # determine the number of insertions based on the unit cell volume of the crystal (BEFORE replication)
    nb_insertions = max(N_BLOCKS, ceil(Int, insertions_per_volume * crystal.box.Ω))
    if verbose
        @printf("\t%d total # particle insertions\n", nb_insertions)
    end

    # partition total insertions among blocks.
    if nprocs() > N_BLOCKS
        error("Use $N_BLOCKS cores or less for Henry coefficient calculations to match the number of blocks")
    end

    nb_insertions_per_block = ceil(Int, nb_insertions / N_BLOCKS)

    # replication factors for applying nearest image convention for short-range interactions
    repfactors = replication_factors(crystal.box, ljforcefield)
    if verbose
        @printf("\tReplicating crystal %d by %d by %d for short-range cutoff %.2f\n",
                repfactors[1], repfactors[2], repfactors[3], sqrt(ljforcefield.r²_cutoff))
    end
    # replicate the crystal atoms so fractional coords are in [0, 1] spanning the simulation box
    crystal = replicate(crystal, repfactors)

    # adjust molecule's fractional coordinates according to the replicated crystal box.
    molecule = Frac(molecule, crystal.box)

    # Bool's of whether to compute guest-host and/or guest-guest electrostatic energies
    #   there is no point in going through the computations if all charges are zero!
    charged_system = has_charges(crystal) & has_charges(molecule)

    # get xtal density for conversion to per mass units (up here in case fails due to missing atoms in atomicmasses.csv)
    ρ = crystal_density(crystal) # kg/m³

    # conduct Monte Carlo insertions for less than 5 cores using Julia pmap function
    # set up function to take a tuple of arguments, the number of insertions to
    # perform and the molecule to move around/rotate. each core needs a different
    # molecule because it will change its attributes in the simulation
    # x = (nb_insertions, molecule, block_no) for that core
    henry_loop(x::Tuple{Int, Molecule, Int}) = _conduct_Widom_insertions(crystal, x[2],
                                            temperature, ljforcefield, x[1],
                                            charged_system, x[3], ewald_precision, verbose,
                                            accessibility_grid, repfactors)

    # parallelize insertions across the cores; keep nb_insertions_per_block same
    res = pmap(henry_loop, [(nb_insertions_per_block, deepcopy(molecule), b) for b = 1:N_BLOCKS])

    # unpack the boltzmann factor sum and weighted energy sum from each block
    boltzmann_factor_sums = [res[b][1] for b = 1:N_BLOCKS] # Σᵢ e^(-βEᵢ) for that core
    wtd_energy_sums = [res[b][2] for b = 1:N_BLOCKS] # Σᵢ Eᵢe^(-βEᵢ) for that core

    # compute block ⟨U⟩, Kₕ
    #   ⟨U⟩ = Σ Uᵢ e ^(βUᵢ) / [ ∑ e^(βUᵢ) ]
    #   Kₕ = β Σ e ^(βUᵢ) / nb_insertions_per_block
    #    (these N_BLOCKS-long arrays)
    average_energies = wtd_energy_sums ./ boltzmann_factor_sums # K
    henry_coefficients = boltzmann_factor_sums / (UNIV_GAS_CONST * temperature * nb_insertions_per_block) # mol/(m³-bar)

    if verbose
        for b = 1:N_BLOCKS
            printstyled(@sprintf("\tBlock  %d/%d statistics:\n", b, N_BLOCKS); color=:yellow)
            println("\tHenry coeff. [mmol/(g-bar)]: ", henry_coefficients[b] / ρ)
            println("\t⟨U, vdw⟩ (K): ", average_energies[b].vdw)
            println("\t⟨U, es⟩ (K): ", average_energies[b].es)
        end
    end

    time_stop = time()
    elapsed_time = (time_stop - time_start) # seconds

    # compute error estimates
    err_kh = 2.0 * std(henry_coefficients) / sqrt(N_BLOCKS)
    err_energy = PotentialEnergy()
    err_energy.vdw = 2.0 * std([average_energies[b].vdw for b = 1:N_BLOCKS]) / sqrt(N_BLOCKS)
    err_energy.es = 2.0 * std([average_energies[b].es for b = 1:N_BLOCKS]) / sqrt(N_BLOCKS)

    results = Dict{String, Any}()
    results["henry coefficient [mol/(m³-bar)]"] = mean(henry_coefficients)
    results["xtal"] = crystal.name
    results["henry coefficient [mmol/(g-bar)]"] = results["henry coefficient [mol/(m³-bar)]"] / ρ
    results["err henry coefficient [mmol/(g-bar)]"] = err_kh / ρ
    results["henry coefficient [mol/(kg-Pa)]"] = results["henry coefficient [mmol/(g-bar)]"] / 100000.0
    # note assumes same # insertions per core.
    results["⟨U, vdw⟩ (K)"] = mean([average_energies[b].vdw for b = 1:N_BLOCKS])
    results["⟨U, es⟩ (K)"] = mean([average_energies[b].es for b = 1:N_BLOCKS])
    results["⟨U⟩ (K)"] = results["⟨U, vdw⟩ (K)"] + results["⟨U, es⟩ (K)"]

    results["⟨U⟩ (kJ/mol)"] = results["⟨U⟩ (K)"] * K_TO_KJ_PER_MOL
    results["⟨U, vdw⟩ (kJ/mol)"] = results["⟨U, vdw⟩ (K)"] * K_TO_KJ_PER_MOL
    results["err ⟨U, vdw⟩ (kJ/mol)"] = err_energy.vdw * K_TO_KJ_PER_MOL
    results["⟨U, es⟩ (kJ/mol)"] = results["⟨U, es⟩ (K)"] * K_TO_KJ_PER_MOL
    results["err ⟨U, es⟩ (kJ/mol)"] = err_energy.es * K_TO_KJ_PER_MOL
    results["Qst (kJ/mol)"] = -results["⟨U⟩ (kJ/mol)"] + temperature * K_TO_KJ_PER_MOL
    results["err Qst (kJ/mol)"] = sum(err_energy) * K_TO_KJ_PER_MOL

    results["elapsed time (min)"] = elapsed_time / 60

    if autosave
        if ! isdir(PATH_TO_SIMS)
            mkdir(PATH_TO_SIMS)
        end
        savename = joinpath(PATH_TO_SIMS, henry_result_savename(crystal, molecule, temperature,
                               ljforcefield, insertions_per_volume, comment=filename_comment))
        @save savename results
        if verbose
            println("\tResults saved in: ", savename)
        end
    end

    if verbose
        println("\tElapsed time (min): ", results["elapsed time (min)"])
        printstyled("\t----- final results ----\n"; color=:green)
        for key in ["henry coefficient [mmol/(g-bar)]", "⟨U, vdw⟩ (kJ/mol)", "⟨U, es⟩ (kJ/mol)", "Qst (kJ/mol)"]
            @printf("\t%s = %f +/- %f\n", key, results[key], results["err " * key])
        end
    end
    return results
end

# assumed crystal is already replicated sufficiently for short-range interactions
# to facilitate parallelization
function _conduct_Widom_insertions(crystal::Crystal, molecule::Molecule,
                                   temperature::Float64, ljforcefield::LJForceField,
                                   nb_insertions::Int, charged_system::Bool, which_block::Int,
                                   ewald_precision::Float64, verbose::Bool,
                                   accessibility_grid::Union{Nothing, Grid{Bool}},
                                   repfactors::Tuple{Int, Int, Int})
                                   
    # copy the molecule in case we need to reset it when bond lengths drift
    bond_length_drift_check_frequency = 5000 # every how many insertions check for drift
    ref_molecule = deepcopy(molecule)

    # define Ewald summation params; this must be done on each core since they cannot share eikar for race condition possibility
    # pre-compute weights on k-vector contributions to long-rage interactions in
    #   Ewald summation for electrostatics
    #   allocate memory for exp^{i * n * k ⋅ r}
    eparams = setup_Ewald_sum(crystal.box, sqrt(ljforcefield.r²_cutoff),
                              verbose=(verbose & (myid() == 1)  & charged_system),
                              ϵ=ewald_precision)
    eikr = Eikr(crystal, eparams)

    # to be Σᵢ Eᵢe^(-βEᵢ)
    wtd_energy_sum = PotentialEnergy(0.0, 0.0)
    # to be Σᵢ e^(-βEᵢ)
    boltzmann_factor_sum = 0.0
    
    insertion_start = 1 # which insertion number to start on
    
    for i = insertion_start:nb_insertions
        # determine uniform random center of mass
        xf = rand(3)
        # translate molecule to the new center of mass
        translate_to!(molecule, Frac(xf))
        # rotate randomly
        if needs_rotations(molecule)
            random_rotation!(molecule, crystal.box)
        end

        # calculate potential energy of molecule at that position and orientation
        energy = PotentialEnergy(0.0, 0.0)
        # actually compute if accessibilty grid not available; or if is avail and accessible
        if isnothing(accessibility_grid) || (accessible(accessibility_grid, xf, repfactors))
            energy.vdw = vdw_energy(crystal, molecule, ljforcefield)
            if charged_system
                energy.es = total(electrostatic_potential_energy(crystal, molecule, eparams, eikr))
            end
        else
            energy.vdw = Inf
            energy.es = Inf
        end

        # calculate Boltzmann factor e^(-βE)
        boltzmann_factor = exp(-sum(energy) / temperature)

        boltzmann_factor_sum += boltzmann_factor

        # to avoid NaN; contribution to wtd energy sum when energy = Inf is zero.
        if (! isinf(energy.vdw)) && (! isinf(energy.es))
            wtd_energy_sum += boltzmann_factor * energy
        end
        # else add 0.0 b/c lim E --> ∞ of E exp(-b * E) is zero.

        # check for drift in bond lengths
        if i % bond_length_drift_check_frequency == 0
            # assert bond length drift not so severe as to not trust results
            if distortion(molecule, ref_molecule, crystal.box,
                    atol=1e-5, throw_warning=true)
                error("significant bond length drift observed after $i insertions.
                Change `bond_length_drift_check_frequency` to restart the molecule configuration
                more frequently.")
            end
            # reset molecule if bond length drift observed but not too severe
            if distortion(molecule, ref_molecule, crystal.box,
                    atol=1e-12, throw_warning=true)
                @warn "... resetting molecule coordinates with fresh copy after $i insertions"
                molecule = deepcopy(ref_molecule)
            end
        end
    end #insertions

    return boltzmann_factor_sum, wtd_energy_sum
end

"""
    save_name = henry_result_savename(crystal, molecule, temperature,
                                   ljforcefield, insertions_per_volume;
                                   comment="")

Determine the name of files saved while calculating the henry coefficient. It uses
many pieces of information from the simulation to ensure the file name accurately
describes what it holds.

# Arguments
- `crystal::Crystal`: The porous crystal being tested
- `molecule::Molecule`: The molecule being tested inside the crystal
- `temperature::Float64`: The temperature used in the simulation units: Kelvin (K)
- `ljforcefield::LJForceField`: The molecular model being used in the simulation
    to describe the intermolecular Van der Waals forces
- `insertions_per_volume::Union{Int, Float64}`: The number of widom insertions per unit volume.
    Will be scaled according to the crystal we're working with
- `comment::AbstractString`: An optional comment that will be appended to the filename
"""
function henry_result_savename(crystal::Crystal, molecule::Molecule, temperature::Float64,
                               ljforcefield::LJForceField, insertions_per_volume::Union{Int, Float64}; comment::AbstractString="")
    if comment != "" && comment[1] != '_'
        comment = "_" * comment
    end
    return @sprintf("henry_sim_%s_in_%s_%fK_%s_ff_%d_insertions_per_volume%s.jld2",
                    molecule.species, crystal.name, temperature, ljforcefield.name,
                    insertions_per_volume, comment)
end
