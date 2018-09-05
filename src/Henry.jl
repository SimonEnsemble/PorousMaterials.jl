const universal_gas_constant = 8.3144598e-5 # m³-bar/(K-mol)
const K_to_kJ_mol = 8.3144598 / 1000.0 # kJ/(mol-K)
using Distributed

"""
   result = henry_coefficient(framework, molecule, temperature, ljforcefield,
                             nb_insertions=1e6, verbose=true, ewald_precision=1e-6,
                             autosave=true)

Conduct particle insertions to compute the Henry coefficient Kₕ of a molecule in a framework.
Also, for free, the heat of adsorption and ensemble average energy of adsorption is computed.
The Henry coefficient is a model for adsorption at infinite dilution (low coverage):
⟨N⟩ = Kₕ P, where P is pressure and Kₕ is the Henry coefficient.

Kₕ = β ⟨e^{-β U}⟩, where the average is over positions and orientations of the molecule
in the framework.

Returns a dictionary of results.

# Arguments
- `framework::Framework`: the porous crystal in which we seek to simulate adsorption
- `molecule::molecule`: the adsorbate molecule
- `temperature::Float64`: temperature of bulk gas phase in equilibrium with adsorbed phase
    in the porous material. units: Kelvin (K)
- `ljforcefield::LJForceField`: the molecular model used to describe the
    energetics of the adsorbate-adsorbate and adsorbate-host van der Waals interactions.
- `insertions_per_volume::Int`: number of Widom insertions to perform for computing the
average, per unit cell volume (Å³)
- `verbose::Bool`: whether or not to print off information during the simulation.
- `ewald_precision::Float64`: desired precision for Ewald summations; used to determine
the replication factors in reciprocal space.
- `autosave::Bool`: save results file as a .jld in PATH_TO_DATA * `sims`
- `filename_comment::AbstractString`: An optional comment that will be appended to the name of the saved file.
"""
function henry_coefficient(framework::Framework, molecule_::Molecule, temperature::Float64,
                           ljforcefield::LJForceField; insertions_per_volume::Int=200,
                           verbose::Bool=true, ewald_precision::Float64=1e-6,
                           autosave::Bool=true, filename_comment::AbstractString="")
    time_start = time()
    if verbose
        print("Simulating Henry coefficient of ")
        printstyled(molecule_.species; color=:green)
        print(" in ")
        printstyled(framework.name; color=:green)
        print(" at ")
        printstyled(temperature; color=:green)
        print(" K, using ")
        printstyled(ljforcefield.name; color=:green)
        print(" force field with ")
        printstyled(insertions_per_volume; color=:green)
        println(" insertions per Å³.")
    end

    # determine the number of insertions based on the unit cell volume of the crystal (BEFORE replication)
    nb_insertions = ceil(Int, insertions_per_volume * framework.box.Ω)
    if verbose
        @printf("\t%d total # particle insertions\n", nb_insertions)
    end

    # partition total insertions among blocks.
    if nprocs() > N_BLOCKS
        error("Use $N_BLOCKS cores or less for Henry coefficient calculations to match the number of blocks")
    end
    nb_insertions_per_block = ceil(Int, nb_insertions / N_BLOCKS)

    # replication factors for applying nearest image convention for short-range interactions
    repfactors = replication_factors(framework.box, ljforcefield)
    if verbose
        @printf("\tReplicating framework %d by %d by %d for short-range cutoff %.2f\n",
                repfactors..., sqrt(ljforcefield.cutoffradius_squared))
    end
    # replicate the framework atoms so fractional coords are in [0, 1] spanning the simulation box
    framework = replicate(framework, repfactors)

    # make a deep copy of the molecule passed so we don't change its coords
    molecule = deepcopy(molecule_)

    # adjust molecule's fractional coordinates according to the replicated framework box.
    set_fractional_coords!(molecule, framework.box)

    # Bool's of whether to compute guest-host and/or guest-guest electrostatic energies
    #   there is no point in going through the computations if all charges are zero!
    charged_system = charged(framework, verbose=verbose) & charged(molecule, verbose=verbose)

    # get xtal density for conversion to per mass units (up here in case fails due to missing atoms in atomicmasses.csv)
    ρ = crystal_density(framework) # kg/m³

    # conduct Monte Carlo insertions for less than 5 cores using Julia pmap function
    # set up function to take a tuple of arguments, the number of insertions to
    # perform and the molecule to move around/rotate. each core needs a different
    # molecule because it will change its attributes in the simulation
    # x = (nb_insertions, molecule) for that core
    henry_loop(x::Tuple{Int, Molecule}) = _conduct_Widom_insertions(framework, x[2],
                                            temperature, ljforcefield, x[1],
                                            charged_system, ewald_precision, verbose)

    # parallelize insertions across the cores; keep nb_insertions_per_block same
    res = pmap(henry_loop, [(nb_insertions_per_block, deepcopy(molecule)) for b = 1:N_BLOCKS])

    # unpack the boltzmann factor sum and weighted energy sum from each block
    boltzmann_factor_sums = [res[b][1] for b = 1:N_BLOCKS] # Σᵢ e^(-βEᵢ) for that core
    wtd_energy_sums = [res[b][2] for b = 1:N_BLOCKS] # Σᵢ Eᵢe^(-βEᵢ) for that core

    # compute block ⟨U⟩, Kₕ
    #   ⟨U⟩ = Σ Uᵢ e ^(βUᵢ) / [ ∑ e^(βUᵢ) ]
    #   Kₕ = β Σ e ^(βUᵢ) / nb_insertions_per_block
    #    (these N_BLOCKS-long arrays)
    average_energies = wtd_energy_sums ./ boltzmann_factor_sums # K
    henry_coefficients = boltzmann_factor_sums / (universal_gas_constant * temperature * nb_insertions_per_block) # mol/(m³-bar)

    if verbose
        for b = 1:N_BLOCKS
            printstyled(@sprintf("\tBlock  %d/%d statistics:\n", b, N_BLOCKS); color=:yellow)
            println("\tHenry coeff. [mmol/(g-bar)]: ", henry_coefficients[b] / ρ)
            println("\t⟨U, vdw⟩ (K): ", average_energies[b].vdw)
            println("\t⟨U, Coulomb⟩ (K): ", average_energies[b].coulomb)
        end
    end

    time_stop = time()
    elapsed_time = (time_stop - time_start) # seconds

    # compute error estimates
    err_kh = 2.0 * std(henry_coefficients) / sqrt(N_BLOCKS)
    err_energy = PotentialEnergy()
    err_energy.vdw = 2.0 * std([average_energies[b].vdw for b = 1:N_BLOCKS]) / sqrt(N_BLOCKS)
    err_energy.coulomb = 2.0 * std([average_energies[b].coulomb for b = 1:N_BLOCKS]) / sqrt(N_BLOCKS)

    result = Dict{String, Float64}()
    result["henry coefficient [mol/(m³-bar)]"] = mean(henry_coefficients)
    result["henry coefficient [mmol/(g-bar)]"] = result["henry coefficient [mol/(m³-bar)]"] / ρ
    result["err henry coefficient [mmol/(g-bar)]"] = err_kh / ρ
    result["henry coefficient [mol/(kg-Pa)]"] = result["henry coefficient [mmol/(g-bar)]"] / 100000.0
    # note assumes same # insertions per core.
    result["⟨U, vdw⟩ (K)"] = mean([average_energies[b].vdw for b = 1:N_BLOCKS])
    result["⟨U, Coulomb⟩ (K)"] = mean([average_energies[b].coulomb for b = 1:N_BLOCKS])
    result["⟨U⟩ (K)"] = result["⟨U, vdw⟩ (K)"] + result["⟨U, Coulomb⟩ (K)"]

    result["⟨U⟩ (kJ/mol)"] = result["⟨U⟩ (K)"] * K_to_kJ_mol
    result["⟨U, vdw⟩ (kJ/mol)"] = result["⟨U, vdw⟩ (K)"] * K_to_kJ_mol
    result["err ⟨U, vdw⟩ (kJ/mol)"] = err_energy.vdw * K_to_kJ_mol
    result["⟨U, Coulomb⟩ (kJ/mol)"] = result["⟨U, Coulomb⟩ (K)"] * K_to_kJ_mol
    result["err ⟨U, Coulomb⟩ (kJ/mol)"] = err_energy.coulomb * K_to_kJ_mol
    result["Qst (kJ/mol)"] = -result["⟨U⟩ (kJ/mol)"] + temperature * K_to_kJ_mol
    result["err Qst (kJ/mol)"] = sum(err_energy) * K_to_kJ_mol

    result["elapsed time (min)"] = elapsed_time / 60

    if autosave
        if ! isdir(PATH_TO_DATA * "henry_sims")
            mkdir(PATH_TO_DATA * "henry_sims")
        end
        savename = PATH_TO_DATA * "henry_sims/" * henry_result_savename(framework, molecule, temperature,
                               ljforcefield, insertions_per_volume, comment=filename_comment)
        @save savename result
        if verbose
            println("\tResults saved in: ", savename)
        end
    end

    if verbose
        println("\tElapsed time (min): ", result["elapsed time (min)"])
        printstyled("\t----- final results ----\n"; color=:green)
        for key in ["henry coefficient [mmol/(g-bar)]", "⟨U, vdw⟩ (kJ/mol)", "⟨U, Coulomb⟩ (kJ/mol)", "Qst (kJ/mol)"]
            @printf("\t%s = %f +/- %f\n", key, result[key], result["err " * key])
        end
    end
    return result
end

# assumed framework is already replicated sufficiently for short-range interactions
# to facilitate parallelization
function _conduct_Widom_insertions(framework::Framework, molecule::Molecule,
                                   temperature::Float64, ljforcefield::LJForceField,
                                   nb_insertions::Int, charged_system::Bool,
                                   ewald_precision::Float64, verbose::Bool)
    # compute pairwise atom & charge distances to later ensure these do not drift in the
    #   course of the simulation (rotations and translations)
    atom_distances = pairwise_atom_distances(molecule, framework.box)
    charge_distances = pairwise_charge_distances(molecule, framework.box)

    # define Ewald summation params; this must be done on each core since they cannot share eikar for race condition possibility
    # pre-compute weights on k-vector contributions to long-rage interactions in
    #   Ewald summation for electrostatics
    #   allocate memory for exp^{i * n * k ⋅ r}
    eparams = setup_Ewald_sum(framework.box, sqrt(ljforcefield.cutoffradius_squared),
                              verbose=(verbose & (myid() == 1)  & charged_system),
                              ϵ=ewald_precision)
    eikr = Eikr(framework, eparams)

    # to be Σᵢ Eᵢe^(-βEᵢ)
    wtd_energy_sum = PotentialEnergy(0.0, 0.0)
    # to be Σᵢ e^(-βEᵢ)
    boltzmann_factor_sum = 0.0

    for i = 1:nb_insertions
        # determine uniform random center of mass
        xf = rand(3)
        # translate molecule to the new center of mass
        translate_to!(molecule, xf)
        # rotate randomly
        if rotatable(molecule)
            rotate!(molecule, framework.box)
        end

        # calculate potential energy of molecule at that position and orientation
        energy = PotentialEnergy(0.0, 0.0)
        energy.vdw = vdw_energy(framework, molecule, ljforcefield)
        if charged_system
            energy.coulomb = total(electrostatic_potential_energy(framework, molecule, eparams, eikr))
        end

        # calculate Boltzmann factor e^(-βE)
        boltzmann_factor = exp(-sum(energy) / temperature)

        boltzmann_factor_sum += boltzmann_factor

        # to avoid NaN; contribution to wtd energy sum when energy = Inf is zero.
        if (! isinf(energy.vdw)) && (! isinf(energy.coulomb))
            wtd_energy_sum += boltzmann_factor * energy
        end
        # else add 0.0 b/c lim E --> ∞ of E exp(-b * E) is zero.
    end

    # assert no significant drift in atom distances (e.g. bond lengths) upon rotations
    if ! isapprox(atom_distances, pairwise_atom_distances(molecule, framework.box), atol=1e-10) ||
       ! isapprox(charge_distances, pairwise_charge_distances(molecule, framework.box), atol=1e-10)
        error("significant drift in bond lenghts during rotations/translations")
    end

    return boltzmann_factor_sum, wtd_energy_sum
end

function henry_result_savename(framework::Framework, molecule::Molecule, temperature::Float64,
                               ljforcefield::LJForceField, insertions_per_volume::Union{Int, Float64}; comment::AbstractString="")
    if comment != "" && comment[1] != '_'
        comment = "_" * comment
    end
    return @sprintf("henry_sim_%s_in_%s_%fK_%s_ff_%d_insertions_per_volume%s.jld2",
                    molecule.species, framework.name, temperature, ljforcefield.name,
                    insertions_per_volume, comment)
end
