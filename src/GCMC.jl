const KB = 1.38064852e7 # Boltmann constant (Pa-m3/K --> Pa-A3/K)

# define Markov chain proposals here.
const PROPOSAL_ENCODINGS = Dict(1 => "insertion", 2 => "deletion", 3 => "translation")
const N_PROPOSAL_TYPES = length(keys(PROPOSAL_ENCODINGS))
const INSERTION   = Dict([v => k for (k, v) in PROPOSAL_ENCODINGS])["insertion"]
const DELETION    = Dict([v => k for (k, v) in PROPOSAL_ENCODINGS])["deletion"]
const TRANSLATION = Dict([v => k for (k, v) in PROPOSAL_ENCODINGS])["translation"]

"""
Data structure to keep track of statistics collected during a grand-canonical Monte Carlo
simulation.

* `n` is the number of molecules in the simulation box.
* `V` is the potential energy.
* `g` refers to guest (the adsorbate molecule).
* `h` refers to host (the crystalline framework).
"""
type GCMCstats
    n_samples::Int

    n::Int
    n²::Int

    U_gh::Float64
    U_gh²::Float64

    U_gg::Float64
    U_gg²::Float64

    U_ggU_gh::Float64
    Un::Float64 # ⟨U n⟩
end

"""
Keep track of Markov chain transitions (proposals and acceptances) during a grand-canonical
Monte Carlo simulation. Entry `i` of these arrays corresponds to PROPOSAL_ENCODINGS[i].
"""
type MarkovCounts
    n_proposed::Array{Int, 1}
    n_accepted::Array{Int, 1}
end

# TODO move this to MC helpers? but not sure if it will inline. so wait after test with @time
@inline function potential_energy(molecule_id::Int, molecules::Array{Molecule, 1}, framework::Framework,
                                  ljforcefield::LennardJonesForceField, simulation_box::Box,
                                  repfactors::Tuple{Int, Int, Int}, sr_cutoff_radius::Float64,
                                  kvectors::Array{Kvector, 1}, α::Float64, charged_molecules::Bool, charged_framework::Bool)
    energy = PotentialEnergy(0.0, 0.0, 0.0, 0.0)
    energy.vdw_gg = vdw_energy(molecule_id, molecules, ljforcefield, simulation_box)
    energy.vdw_gh = vdw_energy(framework, molecules[molecule_id], ljforcefield, repfactors)
    if charged_molecules
        energy.electro_gg = electrostatic_potential_energy(molecules, molecule_id, simulation_box, sr_cutoff_radius, kvectors, α)
        if charged_framework
            energy.electro_gh = electrostatic_potential_energy(framework, molecules[molecule_id], simulation_box, repfactors, sr_cutoff_radius, kvectors, α)
        end
    end
    return energy
end

function adsorption_isotherm(framework::Framework, temperature::Float64,
                        fugacities::Array{Float64}, molecule::Molecule,
                        ljforcefield::LennardJonesForceField;
                        n_burn_cycles::Int=10000, n_sample_cycles::Int=100000,
                        sample_frequency::Int=25, verbose::Bool=false)

    function run_fugacity(fugacity::Float64)
        return gcmc_simulation(framework, temperature, fugacity, molecule,
                        ljforcefield, n_burn_cycles=n_burn_cycles, n_sample_cycles=n_sample_cycles,
                        sample_frequency=sample_frequency, verbose=verbose)
    end

    # for load balancing, longer computation time goes first
    ids = sortperm(fugacities, rev=true)

    # gcmc simulations in parallel
    results = pmap(run_fugacity, fugacities[ids])

    # put results in same order as original fugacity
    return results[[find(x -> x==i, ids)[1] for i = 1:length(ids)]]
end


"""
    results = gcmc_simulation(framework, temperature, fugacity, molecule, ljforcefield;
                              n_sample_cycles=100000, n_burn_cycles=10000,
                              sample_frequency=25, verbose=false)

Runs a grand-canonical (μVT) Monte Carlo simulation of the adsorption of a molecule in a
framework at a particular temperature and fugacity (= pressure for an ideal gas) using a
Lennard Jones force field.

A cycle is defined as max(20, number of adsorbates currently in the system) Markov chain
proposals. Current Markov chain moves implemented are particle insertion/deletion and
translation.

# Arguments
- `framework::Framework`: the porous crystal in which we seek to simulate adsorption
- `temperature::Float64`: temperature of bulk gas phase in equilibrium with adsorbed phase
    in the porous material. units: Kelvin (K)
- `fugacity::Float64`: fugacity of bulk gas phase in equilibrium with adsorbed phase in the
    porous material. Equal to pressure for an ideal gas. units: Pascal (Pa)
- `molecule::Molecule`: a template of the adsorbate molecule of which we seek to simulate
    the adsorption
- `ljforcefield::LennardJonesForceField`: the molecular model used to describe the
    energetics of the adsorbate-adsorbate and adsorbate-host van der Waals interactions.
- `n_burn_cycles::Int`: number of cycles to allow the system to reach equilibrium before
    sampling.
- `n_sample_cycles::Int`: number of cycles used for sampling
- `sample_frequency::Int`: during the sampling cycles, sample e.g. the number of adsorbed
    gas molecules every this number of Markov proposals.
- `verbose::Bool`: whether or not to print off information during the simulation.
"""
function gcmc_simulation(framework::Framework, temperature::Float64, fugacity::Float64,
                         molecule::Molecule, ljforcefield::LennardJonesForceField;
                         n_burn_cycles::Int=10000, n_sample_cycles::Int=100000,
                         sample_frequency::Int=25, verbose::Bool=false)
    if verbose
        pretty_print(molecule.species, framework.name, temperature, fugacity)
    end

    const repfactors = replication_factors(framework.box, ljforcefield)
    const simulation_box = replicate_box(framework.box, repfactors)
    # TODO: assert center of mass is origin and make rotate! take optional argument to assume com is at origin?
    const molecule_template = deepcopy(molecule)

    if ! (check_forcefield_coverage(framework, ljforcefield) & check_forcefield_coverage(molecule, ljforcefield))
        error("Missing atoms from forcefield")
    end

    # Bool's of whether to compute guest-host and/or guest-guest electrostatic energies
    #   there is no point in going through the computations if all charges are zero!
    const charged_framework = charged(framework)
    const charged_molecules = charged(molecule)

    # define Ewald summation params
    #  TODO do this automatically with rules!
    const k_rep_factors = (11, 11, 9)
    const α = 0.265058
    const sr_cutoff_radius = sqrt(ljforcefield.cutoffradius_squared)

    # pre-compute k-vectors for Ewald summation for electrostatics
    const kvectors = compute_kvectors(simulation_box, k_rep_factors, α)

    # TODO in adsorption isotherm dump coords from previous pressure and load them in here.
    # only true if starting with 0 molecules
    system_energy = PotentialEnergy(0.0, 0.0, 0.0, 0.0)
    gcmc_stats = GCMCstats(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    molecules = Molecule[]

    markov_counts = MarkovCounts(zeros(Int, length(PROPOSAL_ENCODINGS)),
                                 zeros(Int, length(PROPOSAL_ENCODINGS)))

    # (n_burn_cycles + n_sample_cycles) is number of outer cycles.
    #   for each outer cycle, peform max(20, # molecules in the system) MC proposals.
    markov_chain_time = 0
    for outer_cycle = 1:(n_burn_cycles + n_sample_cycles), inner_cycle = 1:max(20, length(molecules))
        markov_chain_time += 1

        # choose proposed move randomly; keep track of proposals
        which_move = rand(1:N_PROPOSAL_TYPES)
        markov_counts.n_proposed[which_move] += 1

        if which_move == INSERTION
            insert_molecule!(molecules, simulation_box, molecule_template)

            # compute the potential energy of the inserted molecule
            energy = potential_energy(length(molecules), molecules, framework, ljforcefield,
                                      simulation_box, repfactors, sr_cutoff_radius,
                                      kvectors, α, charged_molecules, charged_framework)

            # Metropolis Hastings Acceptance for Insertion
            if rand() < fugacity * simulation_box.Ω / (length(molecules) * KB *
                    temperature) * exp(-sum(energy) / temperature)
                # accept the move, adjust current_energy
                markov_counts.n_accepted[which_move] += 1

                system_energy += energy
            else
                # reject the move, remove the inserted molecule
                pop!(molecules)
            end
        elseif (which_move == DELETION) && (length(molecules) != 0)
            # propose which molecule to delete
            molecule_id = rand(1:length(molecules))

            # compute the potential energy of the molecule we propose to delete
            energy = potential_energy(molecule_id, molecules, framework, ljforcefield,
                                      simulation_box, repfactors, sr_cutoff_radius,
                                      kvectors, α, charged_molecules, charged_framework)

            # Metropolis Hastings Acceptance for Deletion
            if rand() < length(molecules) * KB * temperature / (fugacity *
                    simulation_box.Ω) * exp(sum(energy) / temperature)
                # accept the deletion, delete molecule, adjust current_energy
                markov_counts.n_accepted[which_move] += 1

                delete_molecule!(molecule_id, molecules)

                system_energy -= energy
            end
        elseif (which_move == TRANSLATION) && (length(molecules) != 0)
            # propose which molecule whose coordinates we should perturb
            molecule_id = rand(1:length(molecules))

            # energy of the molecule before it was translated
            energy_old = potential_energy(molecule_id, molecules, framework, ljforcefield,
                                          simulation_box, repfactors, sr_cutoff_radius,
                                          kvectors, α, charged_molecules, charged_framework)

            old_molecule = translate_molecule!(molecules[molecule_id], simulation_box)

            # energy of the molecule after it is translated
            energy_new = potential_energy(molecule_id, molecules, framework, ljforcefield,
                                          simulation_box, repfactors, sr_cutoff_radius,
                                          kvectors, α, charged_molecules, charged_framework)

            # Metropolis Hastings Acceptance for translation
            if rand() < exp(-(sum(energy_new) - sum(energy_old)) / temperature)
                # accept the move, adjust current energy
                markov_counts.n_accepted[which_move] += 1

                system_energy += energy_new - energy_old
            else
                # reject the move, reset the molecule at molecule_id
                molecules[molecule_id] = deepcopy(old_molecule)
            end
        end # which move the code executes

        # TODO remove after testing.
        for molecule in molecules
            @assert(! outside_box(molecule, simulation_box), "molecule outside box!")
        end

        # sample the current configuration
        if (outer_cycle > n_burn_cycles) && (markov_chain_time % sample_frequency == 0)
            current_energy_gg = system_energy.vdw_gg + system_energy.electro_gg
            current_energy_gh = system_energy.vdw_gh + system_energy.electro_gh

            gcmc_stats.n_samples += 1

            gcmc_stats.n += length(molecules)
            gcmc_stats.n² += length(molecules) ^ 2

            gcmc_stats.U_gh += current_energy_gh
            gcmc_stats.U_gh² += current_energy_gh ^ 2

            gcmc_stats.U_gg += current_energy_gg
            gcmc_stats.U_gg² += current_energy_gg ^ 2

            gcmc_stats.U_ggU_gh += current_energy_gg * current_energy_gh

            gcmc_stats.Un += (current_energy_gg + current_energy_gh) * length(molecules)
        end
    end # finished markov chain proposal moves

    # compute total energy, compare to `current_energy*` variables where were incremented
    # TODO need function for total electrostatic energy...
    total_U_gh = total_guest_host_vdw_energy(framework, molecules, ljforcefield, repfactors)
    total_U_gg = total_guest_guest_vdw_energy(molecules, ljforcefield, simulation_box)
    incremented_U_gg = system_energy.vdw_gg + system_energy.electro_gg
    incremented_U_gh = system_energy.vdw_gh + system_energy.electro_gh

    if ! isapprox(total_U_gh, incremented_U_gh, atol=0.01)
        println("U_gh, incremented = ", incremented_U_gh)
        println("U_gh, computed at end of simulation =", total_U_gh)
        error("guest-host energy incremented improperly")
    end
    if ! isapprox(total_U_gg, incremented_U_gg, atol=0.01)
        println("U_gg, incremented = ", incremented_U_gg)
        println("U_gg, computed at end of simulation =", total_U_gg)
        error("guest-guest energy incremented improperly")
    end

    @assert(markov_chain_time == sum(markov_counts.n_proposed))

    results = Dict{String, Any}()
    results["crystal"] = framework.name
    results["adsorbate"] = molecule.species
    results["forcefield"] = ljforcefield.name
    results["fugacity (Pa)"] = fugacity
    results["temperature (K)"] = temperature
    results["repfactors"] = repfactors

    results["# sample cycles"] = n_sample_cycles
    results["# burn cycles"] = n_burn_cycles

    results["# samples"] = gcmc_stats.n_samples

    results["# samples"] = gcmc_stats.n_samples
    results["⟨N⟩ (molecules)"] = gcmc_stats.n / gcmc_stats.n_samples
    results["⟨N⟩ (molecules/unit cell)"] = results["⟨N⟩ (molecules)"] /
        (repfactors[1] * repfactors[2] * repfactors[3])
    # (molecules/unit cell) * (mol/6.02 * 10^23 molecules) * (1000 mmol/mol) *
    #    (unit cell/framework amu) * (amu/ 1.66054 * 10^-24)
    results["⟨N⟩ (mmol/g)"] = results["⟨N⟩ (molecules/unit cell)"] * 1000 /
        (6.022140857e23 * molecular_weight(framework) * 1.66054e-24)
    results["⟨U_gg⟩ (K)"] = gcmc_stats.U_gg / gcmc_stats.n_samples
    results["⟨U_gh⟩ (K)"] = gcmc_stats.U_gh / gcmc_stats.n_samples
    results["⟨Energy⟩ (K)"] = (gcmc_stats.U_gg + gcmc_stats.U_gh) /
        gcmc_stats.n_samples
    #variances
    results["var(N)"] = (gcmc_stats.n² / gcmc_stats.n_samples) -
        (results["⟨N⟩ (molecules)"] ^ 2)
    results["Q_st (K)"] = temperature - (gcmc_stats.Un / gcmc_stats.n_samples - results["⟨Energy⟩ (K)"] * results["⟨N⟩ (molecules)"]) / results["var(N)"]
    results["var(U_gg)"] = (gcmc_stats.U_gg² / gcmc_stats.n_samples) -
        (results["⟨U_gg⟩ (K)"] ^ 2)
    results["var⟨U_gh⟩"] = (gcmc_stats.U_gh² / gcmc_stats.n_samples) -
        (results["⟨U_gh⟩ (K)"] ^ 2)
    results["var(Energy)"] = ((gcmc_stats.U_gg² + gcmc_stats.U_gh² + 2 *
        gcmc_stats.U_ggU_gh) / gcmc_stats.n_samples) -
        (results["⟨Energy⟩ (K)"] ^ 2)
    # Markov stats
    for (proposal_id, proposal_description) in PROPOSAL_ENCODINGS
        results[@sprintf("Total # %s proposals", proposal_description)] = markov_counts.n_proposed[proposal_id]
        results[@sprintf("Fraction of %s proposals accepted", proposal_description)] = markov_counts.n_accepted[proposal_id] / markov_counts.n_proposed[proposal_id]
    end

    if verbose
        print_results(results)
    end

    return results
end # gcmc_simulation

function print_results(results::Dict)
    @printf("GCMC simulation of %s in %s at %f K and %f Pa = %f bar fugacity.\n\n",
            results["adsorbate"], results["crystal"], results["temperature (K)"],
            results["fugacity (Pa)"], results["fugacity (Pa)"] / 100000.0)

    @printf("Unit cell replication factors: %d %d %d\n\n", results["repfactors"][1],
                                                         results["repfactors"][2],
                                                         results["repfactors"][3])
    # Markov stats
    for (proposal_id, proposal_description) in PROPOSAL_ENCODINGS
        for key in [@sprintf("Total # %s proposals", proposal_description),
                    @sprintf("Fraction of %s proposals accepted", proposal_description)]
            println(key * ": ", results[key])
        end
    end

    println("")
    for key in ["# sample cycles", "# burn cycles", "# samples"]
        println(key * ": ", results[key])
    end


    println("")
    for key in ["⟨N⟩ (molecules)", "⟨N⟩ (molecules/unit cell)",
                "⟨N⟩ (mmol/g)", "⟨U_gg⟩ (K)", "⟨U_gh⟩ (K)", "⟨Energy⟩ (K)",
                "var(N)", "var(U_gg)", "var⟨U_gh⟩", "var(Energy)"]
        println(key * ": ", results[key])
    end

    @printf("Q_st (K) = %f = %f kJ/mol\n", results["Q_st (K)"], results["Q_st (K)"] * 8.314 / 1000.0)
    return
end

function pretty_print(adsorbate::Symbol, frameworkname::String, temperature::Float64, fugacity::Float64)
    print("Simulating adsorption of ")
    print_with_color(:green, adsorbate)
    print(" in ")
    print_with_color(:green, frameworkname)
    print(" at ")
    print_with_color(:green, @sprintf("%f K", temperature))
    print(" and ")
    print_with_color(:green, @sprintf("%f Pa", fugacity))
    println(" (fugacity).")
end
