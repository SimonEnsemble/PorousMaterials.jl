import Base: +, /

const KB = 1.38064852e7 # Boltmann constant (Pa-m3/K --> Pa-A3/K)

# define Markov chain proposals here.
const PROPOSAL_ENCODINGS = Dict(1 => "insertion", 2 => "deletion", 3 => "translation", 4 => "rotation", 5 => "reinsertion")
const N_PROPOSAL_TYPES = length(keys(PROPOSAL_ENCODINGS))
const INSERTION   = Dict([v => k for (k, v) in PROPOSAL_ENCODINGS])["insertion"]
const DELETION    = Dict([v => k for (k, v) in PROPOSAL_ENCODINGS])["deletion"]
const TRANSLATION = Dict([v => k for (k, v) in PROPOSAL_ENCODINGS])["translation"]
const ROTATION = Dict([v => k for (k, v) in PROPOSAL_ENCODINGS])["rotation"]
const REINSERTION = Dict([v => k for (k, v) in PROPOSAL_ENCODINGS])["reinsertion"]

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

    U::SystemPotentialEnergy
    U²::SystemPotentialEnergy

    Un::Float64 # ⟨U n⟩
end
GCMCstats() = GCMCstats(0, 0, 0, SystemPotentialEnergy(), SystemPotentialEnergy(), 0.0)
+(s1::GCMCstats, s2::GCMCstats) = GCMCstats(s1.n_samples + s2.n_samples,
                                            s1.n         + s2.n,
                                            s1.n²        + s2.n²,
                                            s1.U         + s2.U,
                                            s1.U²        + s2.U², 
                                            s1.Un        + s2.Un)
function Base.sum(gcmc_stats::Array{GCMCstats, 1})
    sum_stats = GCMCstats()
    for gs in gcmc_stats
        sum_stats += gs
    end
    return sum_stats
end

"""
    avg_n, err_n, avg_U, err_U = mean_stderr_n_U(gcmc_stats)

Compute average and standard error of the number of molecules and potential energy from an
array of `GCMCstats`, each corresponding to statitics from an independent block/bin during
the simulation. The average from each bin is treated as an independent sample and used to
estimate the error in the estimate as a confidence interval.
"""
function mean_stderr_n_U(gcmc_stats::Array{GCMCstats, 1})
    # ⟨N⟩, ⟨U⟩
    avg_n = sum(gcmc_stats).n / sum(gcmc_stats).n_samples
    avg_U = sum(gcmc_stats).U / (1.0 * sum(gcmc_stats).n_samples)

    avg_n_blocks = [gs.n / gs.n_samples for gs in gcmc_stats]
    err_n = 2.0 * std(avg_n_blocks) / sqrt(length(gcmc_stats))

    err_U = SystemPotentialEnergy()
    for gs in gcmc_stats
        avg_U_this_block = gs.U / (1.0 * gs.n_samples)
        err_U += square(avg_U_this_block - avg_U)
    end
    err_U = sqrt(err_U) / sqrt(length(gcmc_stats) - 1) # std(U) at this pt.
    err_U = err_U * 2.0 / sqrt(length(gcmc_stats))
    return avg_n, err_n, avg_U, err_U
end

function Base.print(gcmc_stats::GCMCstats)
    println("\t# samples: ", gcmc_stats.n_samples)
    println("\t⟨N⟩ (molecules) = ", gcmc_stats.n / gcmc_stats.n_samples)

    println("\t⟨U_gh, vdw⟩ (K) = ",     gcmc_stats.U.guest_host.vdw / gcmc_stats.n_samples)
    println("\t⟨U_gh, Coulomb⟩ (K) = ", gcmc_stats.U.guest_host.coulomb / gcmc_stats.n_samples)
    println("\t⟨U_gg, vdw⟩ (K) = ",     gcmc_stats.U.guest_guest.vdw / gcmc_stats.n_samples)
    println("\t⟨U_gg, Coulomb⟩ (K) = ", gcmc_stats.U.guest_guest.coulomb / gcmc_stats.n_samples)

    println("\t⟨U⟩ (K) = ", sum(gcmc_stats.U) / gcmc_stats.n_samples)
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
# potential energy change after inserting/deleting/perturbing coordinates of molecules[molecule_id]
@inline function potential_energy(molecule_id::Int,
                                  molecules::Array{Molecule, 1},
                                  framework::Framework,
                                  ljforcefield::LJForceField,
                                  eparams::EwaldParams,
                                  kvectors::Array{Kvector, 1},
                                  eikar::OffsetArray{Complex{Float64}},
                                  eikbr::OffsetArray{Complex{Float64}},
                                  eikcr::OffsetArray{Complex{Float64}},
                                  charged_molecules::Bool,
                                  charged_framework::Bool)
    energy = SystemPotentialEnergy()
    energy.guest_guest.vdw = vdw_energy(molecule_id, molecules, ljforcefield, framework.box)
    energy.guest_host.vdw = vdw_energy(framework, molecules[molecule_id], ljforcefield)
    if charged_molecules
        energy.guest_guest.coulomb = total(electrostatic_potential_energy(molecules, molecule_id, eparams, kvectors, eikar, eikbr, eikcr))
        if charged_framework
            energy.guest_host.coulomb = electrostatic_potential_energy(framework, molecules[molecule_id], eparams, kvectors, eikar, eikbr, eikcr)
        end
    end
    return energy
end

"""
    results = stepwise_adsorption_isotherm(framework, molecule, temperature, pressures,
                                  ljforcefield; n_sample_cycles=100000,
                                  n_burn_cycles=10000, sample_frequency=10,
                                  verbose=true, molecules=Molecule[],
                                  ewald_precision=1e-6, eos=:ideal, checkpoint=false)

Run a set of grand-canonical (μVT) Monte Carlo simulations in series. Arguments are the
same as [`gcmc_simulation`](@ref), as this is the function run behind the scenes. An
exception is that we pass an array of pressures. The adsorption isotherm is computed step-
wise, where the ending configuration from the previous simulation (array of molecules) is
passed into the next simulation as a starting point. The ordering of `pressures` is
honored. By giving each simulation a good starting point, (if the next pressure does not
differ significantly from the previous pressure), we can reduce the number of burn cycles
required to reach equilibrium in the Monte Carlo simulation. Also see
[`adsorption_isotherm`](@ref) which runs the μVT simulation at each pressure in parallel.
"""
function stepwise_adsorption_isotherm(framework::Framework, molecule::Molecule, 
    temperature::Float64, pressures::Array{Float64, 1}, ljforcefield::LJForceField;
    n_burn_cycles::Int=1000, n_sample_cycles::Int=5000, sample_frequency::Int=5, 
    verbose::Bool=true, ewald_precision::Float64=1e-6, eos::Symbol=:ideal,
    checkpoint::Bool=false, checkpoint_interval::Int=50)

    results = Dict{String, Any}[] # push results to this array
    molecules = Molecule[] # initiate with empty framework
    for (i, pressure) in enumerate(pressures)
        result, molecules = gcmc_simulation(framework, molecule, temperature, pressure, 
                                            ljforcefield, 
                                            n_burn_cycles=n_burn_cycles,
                                            n_sample_cycles=n_sample_cycles,
                                            sample_frequency=sample_frequency,
                                            verbose=verbose, molecules=molecules,
                                            ewald_precision=ewald_precision, eos=eos,
                                            checkpoint=checkpoint)
        push!(results, result)
    end

    return results
end

"""
    results = adsorption_isotherm(framework, molecule, temperature, pressures,
                                  ljforcefield; n_sample_cycles=100000,
                                  n_burn_cycles=10000, sample_frequency=25,
                                  verbose=false, molecules=Molecule[],
                                  ewald_precision=1e-6, eos=:ideal, checkpoint=false)

Run a set of grand-canonical (μVT) Monte Carlo simulations in parallel. Arguments are the
same as [`gcmc_simulation`](@ref), as this is the function run in parallel behind the scenes.
The only exception is that we pass an array of pressures. To give Julia access to multiple
cores, run your script as `julia -p 4 mysim.jl` to allocate e.g. four cores. See
[Parallel Computing](https://docs.julialang.org/en/stable/manual/parallel-computing/#Parallel-Computing-1).
"""
function adsorption_isotherm(framework::Framework, molecule::Molecule, temperature::Float64,
    pressures::Array{Float64, 1}, ljforcefield::LJForceField; n_burn_cycles::Int=5000, 
    n_sample_cycles::Int=5000, sample_frequency::Int=5, verbose::Bool=true,
    ewald_precision::Float64=1e-6, eos::Symbol=:ideal,
    checkpoint::Bool=false, checkpoint_interval::Int=50)
    # make a function of pressure only to facilitate uses of `pmap`
    run_pressure(pressure::Float64) = gcmc_simulation(framework, molecule, temperature, 
                                                      pressure, ljforcefield,
                                                      n_burn_cycles=n_burn_cycles,
                                                      n_sample_cycles=n_sample_cycles,
                                                      sample_frequency=sample_frequency,
                                                      verbose=verbose,
                                                      ewald_precision=ewald_precision,
                                                      eos=eos, checkpoint = checkpoint)[1] # only return results

    # for load balancing, larger pressures with longer computation time goes first
    ids = sortperm(pressures, rev=true)

    # run gcmc simulations in parallel using Julia's pmap parallel computing function
    results = pmap(run_pressure, pressures[ids])

    # return results in same order as original pressure even though we permuted them for
    #  better load balancing.
    return results[[find(x -> x==i, ids)[1] for i = 1:length(ids)]]
end


"""
    results, molecules = gcmc_simulation(framework, molecule, temperature, pressure,
                                         ljforcefield; n_sample_cycles=5000,
                                         n_burn_cycles=5000, sample_frequency=5,
                                         verbose=false, molecules=Molecule[],
                                         eos=:ideal, checkpoint=checkpoint)

Runs a grand-canonical (μVT) Monte Carlo simulation of the adsorption of a molecule in a
framework at a particular temperature and pressure using a
Lennard Jones force field.

A cycle is defined as max(20, number of adsorbates currently in the system) Markov chain
proposals. Current Markov chain moves implemented are particle insertion/deletion and
translation.

# Arguments
- `framework::Framework`: the porous crystal in which we seek to simulate adsorption
- `temperature::Float64`: temperature of bulk gas phase in equilibrium with adsorbed phase
    in the porous material. units: Kelvin (K)
- `pressure::Float64`: pressure of bulk gas phase in equilibrium with adsorbed phase in the
    porous material. units: bar
- `molecule::Molecule`: a template of the adsorbate molecule of which we seek to simulate
    the adsorption
- `ljforcefield::LJForceField`: the molecular model used to describe the
    energetics of the adsorbate-adsorbate and adsorbate-host van der Waals interactions.
- `n_burn_cycles::Int`: number of cycles to allow the system to reach equilibrium before
    sampling.
- `n_sample_cycles::Int`: number of cycles used for sampling
- `sample_frequency::Int`: during the sampling cycles, sample e.g. the number of adsorbed
    gas molecules every this number of Markov proposals.
- `verbose::Bool`: whether or not to print off information during the simulation.
- `molecules::Array{Molecule, 1}`: a starting configuration of molecules in the framework.
Note that we assume these coordinates are Cartesian, i.e. corresponding to a unit box.
- `eos::Symbol`: equation of state to use for calculation of fugacity from pressure. Default
is ideal gas, where fugacity = pressure.
"""
function gcmc_simulation(framework::Framework, molecule_::Molecule, temperature::Float64, 
    pressure::Float64, ljforcefield::LJForceField; n_burn_cycles::Int=25000, 
    n_sample_cycles::Int=25000, sample_frequency::Int=5, verbose::Bool=true,
    molecules::Array{Molecule, 1}=Molecule[], ewald_precision::Float64=1e-6, 
    eos::Symbol=:ideal, autosave::Bool=true, progressbar::Bool=false,
    load_checkpoint::Union{Bool, AbstractString}=false, write_checkpoint::Bool=false, checkpoint_interval::Int=50)
  
    # Initialize relevant checkpoint data if we either want to load or write cheackpoints
    restarted = false
    checkpoint_path = ""
    if load_checkpoint != false

        if typeof(load_checkpoint) == Bool
            checkpoint_path = PATH_TO_DATA * "gcmc_checkpoints/" * gcmc_result_savename(framework.name, 
                    molecule_.species, ljforcefield.name, temperature, pressure, n_burn_cycles, n_sample_cycles, "_checkpoint")
        elseif typeof(load_checkpoint) == AbstractString || typeof(load_checkpoint) == String
            if contains(load_checkpoint, "gcmc_checkpoints/")
                checkpoint_path = load_checkpoint
            else
                checkpoint_path = "data/gcmc_checkpoints/" * load_checkpoint
            end
        end

        if isfile(checkpoint_path)
            @printf("Restarting simulation from a previous job\n")
            mccheckpoint = load(checkpoint_path, "checkpoint")
            molecules = mccheckpoint.molecules
            restarted = true
        else
            @printf("No checkpoint file found. Starting a new job\n")
        end

    end

    if write_checkpoint && !restarted
        checkpoint_path = PATH_TO_DATA * "gcmc_checkpoints/" * gcmc_result_savename(framework.name, 
                        molecule_.species, ljforcefield.name, temperature, pressure, n_burn_cycles, n_sample_cycles, "_checkpoint")
        mccheckpoint = MCCheckpoint()
    end

    start_time = time()
    # to avoid changing the outside object `molecule_` inside this function, we make
    #  a deep copy of it here. this serves as a template to copy when we insert a new molecule.
    molecule = deepcopy(molecule_)

    if verbose
        pretty_print(molecule.species, framework.name, temperature, pressure, ljforcefield)
    end

    # convert pressure to fugacity using an equation of state
    fugacity = NaN
    if eos == :ideal
       fugacity = pressure * 100000.0 # bar --> Pa
    elseif eos == :PengRobinson
        prgas = PengRobinsonGas(molecule.species)
        gas_props = calculate_properties(prgas, temperature, pressure, verbose=false)
        fugacity = gas_props["fugacity (bar)"] * 100000.0 # bar --> Pa
    else
        error("eos=:ideal and eos=:PengRobinson are only valid options for equation of state.")
    end
    if verbose
        @printf("\t%s EOS fugacity = %f bar\n", eos, fugacity / 100000.0)
    end

    # replication factors for applying nearest image convention for short-range interactions
    repfactors = replication_factors(framework.box, ljforcefield)
    # replicate the framework atoms so fractional coords are in [0, 1] spanning the simulation box
    framework = replicate(framework, repfactors)
    # adjust fractional coords of molecule according to *replicated* framework
    set_fractional_coords!(molecule, framework.box)
    if verbose
        @printf("\tFramework replicated (%d,%d,%d) for short-range cutoff of %f Å\.\n",
            repfactors..., sqrt(ljforcefield.cutoffradius_squared))
        println("\tFramework crystal density: ", crystal_density(framework))
        println("\tFramework chemical formula: ", chemical_formula(framework))
        println("\tTotal number of LJ Spheres: ", length(framework.atoms))
        println("\tTotal number of point charges: ", length(framework.charges))
    end

    # TODO: assert center of mass is origin and make rotate! take optional argument to assume com is at origin?
    const molecule_template = deepcopy(molecule)
    if ! (total_charge(molecule_template) ≈ 0.0)
        error(@sprintf("Molecule %s is not charge neutral!\n", molecule_template.species))
    end

    if ! (check_forcefield_coverage(framework, ljforcefield) & check_forcefield_coverage(molecule_template, ljforcefield))
        error("Missing atoms from forcefield.")
    end

    # Bool's of whether to compute guest-host and/or guest-guest electrostatic energies
    #   there is no point in going through the computations if all charges are zero!
    const charged_framework = charged(framework, verbose=verbose)
    const charged_molecules = charged(molecule_template, verbose=verbose)

    # define Ewald summation params
    # pre-compute weights on k-vector contributions to long-rage interactions in
    #   Ewald summation for electrostatics
    #   allocate memory for exp^{i * n * k ⋅ r}
    eparams, kvectors, eikar, eikbr, eikcr = setup_Ewald_sum(sqrt(ljforcefield.cutoffradius_squared), framework.box,
                        verbose=verbose & (charged_framework || charged_molecules), 
                        ϵ=ewald_precision)

    # initiate system energy to which we increment when MC moves are accepted
    system_energy = load_checkpoint != false && restarted ? mccheckpoint.system_energy : SystemPotentialEnergy()
    # if we don't start with an emtpy framework, compute energy of starting configuration
    #  (n=0 corresponds to zero energy)
    if length(molecules) != 0
        # ensure molecule template matches species of starting molecules.
        assert(all([m.species == molecule_template.species for m in molecules]))
        
        # set fractional coords of these molecules 
        set_fractional_coords!.(molecules, framework.box)

        # assert that the bond lengths are equal between the template and array to make
        # sure the right fractional coords were used
        if length(molecule.atoms) > 1
            template_bond_length = norm(framework.box.f_to_c * (molecule_template.atoms[1].xf - molecule_template.atoms[2].xf))
            for m in molecules
                bond_length = norm(framework.box.f_to_c * (m.atoms[1].xf - m.atoms[2].xf))
                if ! isapprox(bond_length, template_bond_length, atol=1e-6)
                    error("A bond length between atoms in a molecule in `molecules` passed 
                    in as an initial configuration is not equal to the molecule template
                    passed.")
                end
            end
a       end

        # assert that the molecules are inside the simulation box
        for m in molecules
            if outside_box(m)
                error("A molecule in `molecules` passed to `gcmc_simulation` as a starting
                configuation of molecules is outside of the framework box!")
            end
        end

        if !(load_checkpoint != false && restarted)
            system_energy.guest_host.vdw = total_vdw_energy(framework, molecules, ljforcefield)
            system_energy.guest_guest.vdw = total_vdw_energy(molecules, ljforcefield, framework.box)
            system_energy.guest_host.coulomb = total_electrostatic_potential_energy(framework, molecules,
                                                        eparams, kvectors, eikar, eikbr, eikcr)
            system_energy.guest_guest.coulomb = total(electrostatic_potential_energy(molecules,
                                                eparams, kvectors, eikar, eikbr, eikcr))
        end
    end

    if progressbar
        progress_bar = Progress(n_burn_cycles + n_sample_cycles, 1)
    end

    # define probabilty of proposing each type of MC move here.
    mc_proposal_probabilities = [0.0 for p = 1:N_PROPOSAL_TYPES]
    # set defaults
    mc_proposal_probabilities[INSERTION] = 0.35
    mc_proposal_probabilities[DELETION] = mc_proposal_probabilities[INSERTION] # must be equal
    mc_proposal_probabilities[REINSERTION] = 0.05
    if rotatable(molecule)
        mc_proposal_probabilities[TRANSLATION] = 0.125
        mc_proposal_probabilities[ROTATION] = 0.125
    else
        mc_proposal_probabilities[TRANSLATION] = 0.25
        mc_proposal_probabilities[ROTATION] = 0.0
    end
    mc_proposal_probabilities /= sum(mc_proposal_probabilities) # normalize
    # StatsBase.jl functionality for sampling
    mc_proposal_probabilities = ProbabilityWeights(mc_proposal_probabilities)
    if verbose
        for p = 1:N_PROPOSAL_TYPES
            @printf("\tProbability of %s: %f\n", PROPOSAL_ENCODINGS[p], mc_proposal_probabilities[p])
        end
    end

    # initiate GCMC statistics for each block # break simulation into `N_BLOCKS` blocks to gauge convergence
    gcmc_stats = load_checkpoint != false && restarted ? mccheckpoint.blocks : [GCMCstats() for block_no = 1:N_BLOCKS]
    current_block = load_checkpoint != false && restarted ? mccheckpoint.current_block : 1
    # make sure the number of sample cycles is at least equal to N_BLOCKS
    if n_sample_cycles < N_BLOCKS
        n_sample_cycles = N_BLOCKS
        warn(@sprintf("# sample cycles set to minimum %d, which is number of blocks.", N_BLOCKS))
    end
    const N_CYCLES_PER_BLOCK = floor(Int, n_sample_cycles / N_BLOCKS)

    markov_counts = load_checkpoint != false && restarted ? mccheckpoint.markov_counts : MarkovCounts(zeros(Int, length(PROPOSAL_ENCODINGS)), zeros(Int, length(PROPOSAL_ENCODINGS)))
    #println(rand())

    # (n_burn_cycles + n_sample_cycles) is number of outer cycles.
    #   for each outer cycle, peform max(20, # molecules in the system) MC proposals.
    markov_chain_time = load_checkpoint != false && restarted ? mccheckpoint.markov_chain_time : 0
    outer_start = load_checkpoint != false && restarted ? mccheckpoint.n_cycle : 1
    for outer_cycle = outer_start:(n_burn_cycles + n_sample_cycles)
        if progressbar
            next!(progress_bar; showvalues=[(:cycle, outer_cycle), (:number_of_molecules, length(molecules))])
        end
        end_inner_cycle = max(20, length(molecules))
        for inner_cycle = 1:end_inner_cycle
            markov_chain_time += 1

            # choose proposed move randomly; keep track of proposals
            which_move = sample(1:N_PROPOSAL_TYPES, mc_proposal_probabilities) # StatsBase.jl
            markov_counts.n_proposed[which_move] += 1

            if which_move == INSERTION
                insert_molecule!(molecules, framework.box, molecule_template)

                energy = potential_energy(length(molecules), molecules, framework, 
                                                ljforcefield, eparams, kvectors, eikar, eikbr, 
                                                eikcr, charged_molecules, charged_framework)

                # Metropolis Hastings Acceptance for Insertion
                if rand() < fugacity * framework.box.Ω / (length(molecules) * KB *
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
                                          eparams, kvectors, eikar, eikbr, eikcr, 
                                          charged_molecules, charged_framework)

                # Metropolis Hastings Acceptance for Deletion
                if rand() < length(molecules) * KB * temperature / (fugacity *
                        framework.box.Ω) * exp(sum(energy) / temperature)
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
                                          eparams,
                                          kvectors, eikar, eikbr, eikcr, charged_molecules, charged_framework)

                old_molecule = translate_molecule!(molecules[molecule_id], framework.box)

                # energy of the molecule after it is translated
                energy_new = potential_energy(molecule_id, molecules, framework, ljforcefield,
                                              eparams, kvectors, eikar, eikbr, eikcr, 
                                              charged_molecules, charged_framework)

                # Metropolis Hastings Acceptance for translation
                if rand() < exp(-(sum(energy_new) - sum(energy_old)) / temperature)
                    # accept the move, adjust current energy
                    markov_counts.n_accepted[which_move] += 1

                    system_energy += energy_new - energy_old
                else
                    # reject the move, put back the old molecule
                    molecules[molecule_id] = deepcopy(old_molecule)
                end
            elseif (which_move == ROTATION) && (length(molecules) != 0)
                # propose which molecule to rotate
                molecule_id = rand(1:length(molecules))
                
                # energy of the molecule before we rotate it
                energy_old = potential_energy(molecule_id, molecules, framework, ljforcefield,
                                              eparams, kvectors, eikar, eikbr, eikcr, 
                                              charged_molecules, charged_framework)
                
                # store old molecule to restore old position in case move is rejected
                old_molecule = deepcopy(molecules[molecule_id])

                # conduct a random rotation
                rotate!(molecules[molecule_id], framework.box)

                # energy of the molecule after it is translated
                energy_new = potential_energy(molecule_id, molecules, framework, ljforcefield,
                                              eparams, kvectors, eikar, eikbr, eikcr, 
                                              charged_molecules, charged_framework)

                # Metropolis Hastings Acceptance for rotation 
                if rand() < exp(-(sum(energy_new) - sum(energy_old)) / temperature)
                    # accept the move, adjust current energy
                    markov_counts.n_accepted[which_move] += 1

                    system_energy += energy_new - energy_old
                else
                    # reject the move, put back the old molecule
                    molecules[molecule_id] = deepcopy(old_molecule)
                end
            elseif (which_move == REINSERTION) && (length(molecules) != 0)
                # propose which molecule to re-insert
                molecule_id = rand(1:length(molecules))

                # compute the potential energy of the molecule we propose to re-insert
                energy_old = potential_energy(molecule_id, molecules, framework, ljforcefield,
                                             eparams,
                                             kvectors, eikar, eikbr, eikcr, charged_molecules,
                                             charged_framework)

                # reinsert molecule; store old configuration of the molecule in case proposal is rejected
                old_molecule = reinsert_molecule!(molecules[molecule_id], framework.box)

                # compute the potential energy of the molecule in its new configuraiton
                energy_new = potential_energy(molecule_id, molecules, framework, ljforcefield,
                                              eparams,
                                             kvectors, eikar, eikbr, eikcr, charged_molecules,
                                             charged_framework)

                # Metropolis Hastings Acceptance for reinsertion
                if rand() < exp(-(sum(energy_new) - sum(energy_old)) / temperature)
                    # accept the move, adjust current energy
                    markov_counts.n_accepted[which_move] += 1

                    system_energy += energy_new - energy_old
                else
                    # reject the move, put back old molecule
                    molecules[molecule_id] = deepcopy(old_molecule)
                end
            end # which move the code executes

            # TODO remove after testing.
            for m in molecules
                @assert(! outside_box(m), "molecule outside box!")
            end

            # if we're in the production MC cycles. i.e. we've done all burn cycles...
            if outer_cycle > n_burn_cycles
                # take a sample.
                if markov_chain_time % sample_frequency == 0
                    gcmc_stats[current_block].n_samples += 1

                    gcmc_stats[current_block].n += length(molecules)
                    gcmc_stats[current_block].n² += length(molecules) ^ 2

                    gcmc_stats[current_block].U += system_energy
                    gcmc_stats[current_block].U² += square(system_energy)

                    gcmc_stats[current_block].Un += sum(system_energy) * length(molecules)
                end

                # print block statistics if last inner cycle of the right outer cycle.
                if (inner_cycle == end_inner_cycle) && ((outer_cycle - n_burn_cycles) % N_CYCLES_PER_BLOCK == 0) &&
                    (current_block != N_BLOCKS || (outer_cycle == (n_burn_cycles + n_sample_cycles)))
                    # print statistics for this block
                    if verbose 
                        print_with_color(:yellow, @sprintf("\tBlock  %d/%d statistics:\n", current_block, N_BLOCKS))
                        print(gcmc_stats[current_block])
                    end
                    # move onto new block unless current_block is N_BLOCKS;
                    # then just keep adding stats to the last block.
                    # this only occurs if sample_cycles not divisible by N_BLOCKS
                    if current_block != N_BLOCKS
                        current_block += 1
                    end
                end
            end # end sampling code
        end # inner
    

    if write_checkpoint && outer_cycle % checkpoint_interval == 0
        mccheckpoint.n_cycle = outer_cycle + 1
        mccheckpoint.current_block = current_block
        mccheckpoint.blocks = gcmc_stats
        mccheckpoint.molecules = get_fractional_coords_to_unit_cube.(molecules, framework.box)
        mccheckpoint.system_energy = system_energy
        mccheckpoint.markov_counts = markov_counts
        mccheckpoint.markov_chain_time = markov_chain_time
        mccheckpoint.time = time() - start_time
        JLD.save(checkpoint_path, "checkpoint", mccheckpoint)
    end

    end # outer cycles
    # finished MC moves at this point.

    # compute total energy, compare to `current_energy*` variables where were incremented
    system_energy_end = SystemPotentialEnergy()
    system_energy_end.guest_host.vdw = total_vdw_energy(framework, molecules, ljforcefield)
    system_energy_end.guest_guest.vdw = total_vdw_energy(molecules, ljforcefield, framework.box)
    system_energy_end.guest_host.coulomb = total_electrostatic_potential_energy(framework, molecules,
                                       eparams, kvectors, eikar, eikbr, eikcr)
    system_energy_end.guest_guest.coulomb = total(total_electrostatic_potential_energy(molecules,
                                        eparams, kvectors, eikar, eikbr, eikcr))

    # see Energetics_Util.jl for this function, overloaded isapprox to print mismatch
    if ! isapprox(system_energy, system_energy_end, verbose=true, atol=0.01)
        error("energy incremented improperly during simulation...")
    end

    @assert(markov_chain_time == sum(markov_counts.n_proposed))
    @printf("Estimated elapsed time: %d seconds\n", load_checkpoint != false && restarted ? time() - start_time + mccheckpoint.time : time() - start_time)

    # build dictionary containing summary of simulation results for easy querying
    results = Dict{String, Any}()
    results["crystal"] = framework.name
    results["adsorbate"] = molecule.species
    results["forcefield"] = ljforcefield.name
    results["pressure (bar)"] = pressure
    results["fugacity (bar)"] = fugacity / 100000.0
    results["temperature (K)"] = temperature
    results["repfactors"] = repfactors

    results["# sample cycles"] = n_sample_cycles
    results["# burn cycles"] = n_burn_cycles
    results["# samples"] = sum(gcmc_stats).n_samples

    # statistics from samples during simulation
    # see here: https://cs.nyu.edu/courses/fall06/G22.2112-001/MonteCarlo.pdf for how
    # error bars are computed; simulation broken into N_BLOCKS and each average from the
    # block is treated as an independent sample.
    avg_n, err_n, avg_U, err_U = mean_stderr_n_U(gcmc_stats)

    # averages
    results["⟨N⟩ (molecules)"]     = avg_n
    results["⟨U_gh, vdw⟩ (K)"]     = avg_U.guest_host.vdw
    results["⟨U_gh, electro⟩ (K)"] = avg_U.guest_host.coulomb
    results["⟨U_gg, vdw⟩ (K)"]     = avg_U.guest_guest.vdw
    results["⟨U_gg, electro⟩ (K)"] = avg_U.guest_guest.coulomb
    results["⟨U⟩ (K)"] = sum(avg_U)

    # variances
    results["var(N)"] = (sum(gcmc_stats).n² / sum(gcmc_stats).n_samples) -
        (results["⟨N⟩ (molecules)"] ^ 2)
    # isosteric heat of adsorption TODO stdev of this too.
    results["Q_st (K)"] = temperature - (sum(gcmc_stats).Un / sum(gcmc_stats).n_samples - results["⟨U⟩ (K)"] * results["⟨N⟩ (molecules)"]) / results["var(N)"]

    # error bars (confidence intervals)
    results["err ⟨N⟩ (molecules)"]     = err_n
    results["err ⟨U_gh, vdw⟩ (K)"]     = err_U.guest_host.vdw
    results["err ⟨U_gh, electro⟩ (K)"] = err_U.guest_host.coulomb
    results["err ⟨U_gg, vdw⟩ (K)"]     = err_U.guest_guest.vdw
    results["err ⟨U_gg, electro⟩ (K)"] = err_U.guest_guest.coulomb
    results["err ⟨U⟩ (K)"] = sum(err_U)


    # average N in more common units
    results["⟨N⟩ (molecules/unit cell)"] = avg_n / (repfactors[1] * repfactors[2] * repfactors[3])
    results["err ⟨N⟩ (molecules/unit cell)"] = err_n / (repfactors[1] * repfactors[2] * repfactors[3])
    # (molecules/unit cell) * (mol/6.02 * 10^23 molecules) * (1000 mmol/mol) *
    #    (unit cell/framework amu) * (amu/ 1.66054 * 10^-24)
    results["⟨N⟩ (mmol/g)"] = results["⟨N⟩ (molecules/unit cell)"] * 1000 /
        (6.022140857e23 * molecular_weight(framework) * 1.66054e-24) * (repfactors[1] * repfactors[2] * repfactors[3])
    results["err ⟨N⟩ (mmol/g)"] = results["err ⟨N⟩ (molecules/unit cell)"] * 1000 /
        (6.022140857e23 * molecular_weight(framework) * 1.66054e-24) * (repfactors[1] * repfactors[2] * repfactors[3])

    # Markov stats
    for (proposal_id, proposal_description) in PROPOSAL_ENCODINGS
        results[@sprintf("Total # %s proposals", proposal_description)] = markov_counts.n_proposed[proposal_id]
        results[@sprintf("Fraction of %s proposals accepted", proposal_description)] = markov_counts.n_accepted[proposal_id] / markov_counts.n_proposed[proposal_id]
    end

    if verbose
        print_results(results, print_title=false)
    end

    # before returning molecules, convert coords back to Cartesian.
    set_fractional_coords_to_unit_cube!.(molecules, framework.box)

    if autosave
        if ! isdir(PATH_TO_DATA * "gcmc_sims")
            mkdir(PATH_TO_DATA * "gcmc_sims")
        end
    
        save_results_filename = PATH_TO_DATA * "gcmc_sims/" * gcmc_result_savename(framework.name, 
            molecule.species, ljforcefield.name, temperature, pressure, n_burn_cycles, n_sample_cycles)

        JLD.save(save_results_filename, "results", results)
        if verbose
            println("\tResults dictionary saved in ", save_results_filename)
        end
    end

    return results, molecules # summary of statistics and ending configuration of molecules
end # gcmc_simulation

"""
Determine the name of files saved during the GCMC simulation, be molecule positions or results.
"""
function gcmc_result_savename(framework_name::AbstractString,
                            molecule_species::Symbol,
                            ljforcefield_name::AbstractString,
                            temperature::Float64,
                            pressure::Float64,
                            n_burn_cycles::Int,
                            n_sample_cycles::Int,
                            comment::AbstractString="")
        framework_name = split(framework_name, ".")[1] # remove file extension
        ljforcefield_name = split(ljforcefield_name, ".")[1] # remove file extension
        if comment != "" && comment[1] != '_'
            comment = "_" * comment
        end
        return @sprintf("gcmc_%s_%s_T%f_P%f_%s_%dburn_%dsample%s.jld", framework_name,
                    molecule_species, temperature, pressure, ljforcefield_name,
                    n_burn_cycles, n_sample_cycles, comment)
end

function print_results(results::Dict; print_title::Bool=true)
    if print_title
        # already print in GCMC tests...
        @printf("GCMC simulation of %s in %s at %f K and %f bar pressure, %f bar fugacity using %s forcefield.\n\n",
                results["adsorbate"], results["crystal"], results["temperature (K)"],
                results["pressure (bar)"], results["fugacity (bar)"] / 100000.0, results["forcefield"])
    end

    @printf("\nUnit cell replication factors: %d %d %d\n\n", results["repfactors"][1],
                                                             results["repfactors"][2],
                                                             results["repfactors"][3])
    # Markov stats
    println("")
    for key in ["# sample cycles", "# burn cycles", "# samples"]
        println(key * ": ", results[key])
    end

    for (proposal_id, proposal_description) in PROPOSAL_ENCODINGS
        total_proposals = results[@sprintf("Total # %s proposals", proposal_description)]
        fraction_accepted = results[@sprintf("Fraction of %s proposals accepted", proposal_description)]
        if total_proposals > 0
            print_with_color(:yellow, proposal_description)
            @printf("\t%d total proposals.\n", total_proposals)
            @printf("\t%f %% proposals accepted.\n", 100.0 * fraction_accepted)
        end
    end

    println("")
    for key in ["⟨N⟩ (molecules)", "⟨N⟩ (molecules/unit cell)", "⟨N⟩ (mmol/g)",
                "⟨U_gg, vdw⟩ (K)", "⟨U_gh, vdw⟩ (K)", "⟨U_gg, electro⟩ (K)",
                "⟨U_gh, electro⟩ (K)", "⟨U⟩ (K)"]
        @printf("%s: %f +/- %f\n", key, results[key], results["err " * key])
        if key == "⟨N⟩ (mmol/g)"
            println("")
        end
    end

    @printf("\nQ_st (K) = %f = %f kJ/mol\n\n", results["Q_st (K)"], results["Q_st (K)"] * 8.314 / 1000.0)
    return
end

function pretty_print(adsorbate::Symbol, frameworkname::String, temperature::Float64, 
                      pressure::Float64, ljff::LJForceField)
    print("Simulating ")
    print_with_color(:yellow, "(μVT)")
    print(" adsorption of ")
    print_with_color(:green, adsorbate)
    print(" in ")
    print_with_color(:green, frameworkname)
    print(" at ")
    print_with_color(:green, @sprintf("%f K", temperature))
    print(" and ")
    print_with_color(:green, @sprintf("%f bar", pressure))
    print(" (bar) with ")
    print_with_color(:green, split(ljff.name, ".")[1])
    println(" force field.")
end

mutable struct MCCheckpoint
    n_cycle::Int
    current_block::Int
    blocks::Array{GCMCstats, 1}
    molecules::Array{Molecule, 1}
    system_energy::SystemPotentialEnergy
    markov_counts::MarkovCounts
    markov_chain_time::Int
    time::Float64
end

MCCheckpoint() = MCCheckpoint(0, 0, Array{GCMCstats, 1}(), Array{Molecule, 1}(), SystemPotentialEnergy(), MarkovCounts(zeros(Int, 1), zeros(Int, 1)), 0, 0.0) 

function get_fractional_coords_to_unit_cube(molecule::Molecule, box::Box)
    new_mol = deepcopy(molecule)
    for ljsphere in new_mol.atoms
        ljsphere.xf[:] = box.f_to_c * ljsphere.xf
    end
    for charge in new_mol.charges
        charge.xf[:] = box.f_to_c * charge.xf
    end
    new_mol.xf_com[:] = box.f_to_c * new_mol.xf_com
    return new_mol
end
