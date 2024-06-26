import Base: +, /

const KB = 1.38064852e7 # Boltmann constant (Pa-m3/K --> Pa-A3/K)

###
#   Markov chain proposals
###
const PROPOSAL_ENCODINGS = Dict(1 => "translation",
                                2 => "rotation",
                                3 => "reinsertion"
                                ) # helps with printing later
const N_PROPOSAL_TYPES = length(keys(PROPOSAL_ENCODINGS))
# each proposal type gets an Int for clearer code
const TRANSLATION = Dict([v => k for (k, v) in PROPOSAL_ENCODINGS])["translation"]
const ROTATION    = Dict([v => k for (k, v) in PROPOSAL_ENCODINGS])["rotation"]
const REINSERTION = Dict([v => k for (k, v) in PROPOSAL_ENCODINGS])["reinsertion"]

 # count proposed/accepted for each subtype
mutable struct MarkovCounts
    n_proposed::Array{Int, 1}
    n_accepted::Array{Int, 1}
end

# TODO move this to MC helpers? but not sure if it will inline. so wait after test with @time
# potential energy change after inserting/deleting/perturbing coordinates of molecules[molecule_id]
@inline function potential_energy(molecule_id::Int,
                                  molecules::Array{Molecule{Frac}, 1},
                                  xtal::Crystal,
                                  ljff::LJForceField,
                                  eparams::EwaldParams,
                                  eikr_gh::Eikr,
                                  eikr_gg::Eikr)
    energy = SystemPotentialEnergy()
    # van der Waals interactions
    energy.gg.vdw = vdw_energy(molecule_id, molecules, ljff, xtal.box) # guest-guest
    energy.gh.vdw = vdw_energy(xtal, molecules[molecule_id], ljff)     # guest-host
    # electrostatic interactions
    if has_charges(molecules[molecule_id])
        # guest-guest
        energy.gg.es = total(electrostatic_potential_energy(molecules, molecule_id, eparams, xtal.box, eikr_gg))
        # guest-host
        if has_charges(xtal)
            energy.gh.es = total(electrostatic_potential_energy(xtal, molecules[molecule_id], eparams, eikr_gh))
        end
    end
    return energy
end

function NVT_sim(xtal::Crystal,
                 molecule_template::Molecule{Cart},
                 molecules::Array{Molecule{Cart}, 1},
                 temperature::Float64,
                 ljff::LJForceField;
                 n_burn_cycles::Int=5000,
                 n_sample_cycles::Int=5000,
                 verbose=true,
                 ewald_precision::Float64=1e-6
                )
    assert_P1_symmetry(xtal)

    start_time = time()
    # # to avoid changing the outside object `molecule_` inside this function, we make
    # #  a deep copy of it here. this serves as a template to copy when we insert a new molecule.
    # molecule = deepcopy(molecule_)

    ###
    #  NON-RIGOROUS PLACEHOLDER -- DELETE
    ###
    molecule = deepcopy(molecule_template)
    
    ###
    #   replicate xtal so that nearest image convention can be applied for short-range interactions
    ###
    repfactors = replication_factors(xtal.box, ljff)
    xtal = replicate(xtal, repfactors) # frac coords still in [0, 1]
    
    # convert molecules array to fractional using this box.
    molecules = Frac.(molecules, xtal.box)

    if ! neutral(molecule.charges)
        error(@sprintf("Molecule %s is not charge neutral!\n", molecule.species))
    end

    if ! (forcefield_coverage(xtal.atoms, ljff) & forcefield_coverage(molecule.atoms, ljff))
        error("Missing atoms from forcefield.")
    end

    # define Ewald summation params
    # pre-compute weights on k-vector contributions to long-rage interactions in
    #   Ewald summation for electrostatics
    #   allocate memory for exp^{i * n * k ⋅ r}
    eparams = setup_Ewald_sum(xtal.box, sqrt(ljff.r²_cutoff),
                        verbose=verbose & (has_charges(molecule) || has_charges(xtal)),
                        ϵ=ewald_precision)
    eikr_gh = Eikr(xtal, eparams)
    eikr_gg = Eikr(molecule, eparams)

    # initiate system energy to which we increment when MC moves are accepted
    system_energy = SystemPotentialEnergy()
    # if we don't start with an emtpy xtal, compute energy of starting configuration
    #  (n=0 corresponds to zero energy)
    if length(molecules) != 0
        # some checks
        for m in molecules
            # ensure molecule template matches species of starting molecules.
            @assert m.species == molecule.species "initializing with wrong molecule species"
            # assert that the molecules are inside the simulation box
            @assert inside(m) "initializing with molecules outside simulation box!"
            # ensure pair-wise bond distances match template
            @assert ! distortion(m, Frac(molecule, xtal.box), xtal.box)
        end

        system_energy.gh.vdw = total_vdw_energy(xtal, molecules, ljff)
        system_energy.gg.vdw = total_vdw_energy(molecules, ljff, xtal.box)
        system_energy.gh.es  = total(total_electrostatic_potential_energy(xtal, molecules, eparams, eikr_gh))
        system_energy.gg.es  = total(electrostatic_potential_energy(molecules, eparams, xtal.box, eikr_gg))
    end

    ####
    #   proposal probabilities
    ###
    mc_proposal_probabilities = [0.0 for p = 1:N_PROPOSAL_TYPES]
    # set defaults
    mc_proposal_probabilities[REINSERTION] = 0.05
    if needs_rotations(molecule)
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
        println("\tMarkov chain proposals:")
        for p = 1:N_PROPOSAL_TYPES
            @printf("\t\tprobability of %s: %f\n", PROPOSAL_ENCODINGS[p], mc_proposal_probabilities[p])
        end
    end

    markov_counts = MarkovCounts(zeros(Int, length(PROPOSAL_ENCODINGS)), zeros(Int, length(PROPOSAL_ENCODINGS)))

    # (n_burn_cycles + n_sample_cycles) is number of outer cycles.
    #   for each outer cycle, peform max(20, # molecules in the system) MC proposals.
    markov_chain_time = 0
    outer_cycle_start = 1
    for outer_cycle = outer_cycle_start:(n_burn_cycles + n_sample_cycles)
        for inner_cycle = 1:max(20, length(molecules))
            markov_chain_time += 1

            # choose proposed move randomly; keep track of proposals
            which_move = sample(1:N_PROPOSAL_TYPES, mc_proposal_probabilities) # StatsBase.jl
            markov_counts.n_proposed[which_move] += 1

            if which_move == TRANSLATION
                # propose which molecule whose coordinates we should perturb
                molecule_id = rand(1:length(molecules))

                # energy of the molecule before it was translated
                energy_old = potential_energy(molecule_id, molecules, xtal, ljff,
                                              eparams, eikr_gh, eikr_gg)

                old_molecule = random_translation!(molecules[molecule_id], xtal.box)

                # energy of the molecule after it is translated
                energy_new = potential_energy(molecule_id, molecules, xtal, ljff,
                                              eparams, eikr_gh, eikr_gg)

                # Metropolis Hastings Acceptance for translation
                if rand() < exp(-(sum(energy_new) - sum(energy_old)) / temperature)
                    # accept the move, adjust current energy
                    markov_counts.n_accepted[which_move] += 1

                    system_energy += energy_new - energy_old
                else
                    # reject the move, put back the old molecule
                    molecules[molecule_id] = deepcopy(old_molecule)
                end
            elseif which_move == ROTATION
                # propose which molecule to rotate
                molecule_id = rand(1:length(molecules))

                # energy of the molecule before we rotate it
                energy_old = potential_energy(molecule_id, molecules, xtal, ljff,
                                              eparams, eikr_gh, eikr_gg)

                # store old molecule to restore old position in case move is rejected
                old_molecule = deepcopy(molecules[molecule_id])

                # conduct a random rotation
                random_rotation!(molecules[molecule_id], xtal.box)

                # energy of the molecule after it is translated
                energy_new = potential_energy(molecule_id, molecules, xtal, ljff,
                                              eparams, eikr_gh, eikr_gg)

                # Metropolis Hastings Acceptance for rotation
                if rand() < exp(-(sum(energy_new) - sum(energy_old)) / temperature)
                    # accept the move, adjust current energy
                    markov_counts.n_accepted[which_move] += 1

                    system_energy += energy_new - energy_old
                else
                    # reject the move, put back the old molecule
                    molecules[molecule_id] = deepcopy(old_molecule)
                end
            elseif which_move == REINSERTION
                # propose which molecule to re-insert
                molecule_id = rand(1:length(molecules))

                # compute the potential energy of the molecule we propose to re-insert
                energy_old = potential_energy(molecule_id, molecules, xtal, ljff,
                                             eparams, eikr_gh, eikr_gg)

                # reinsert molecule; store old configuration of the molecule in case proposal is rejected
                old_molecule = random_reinsertion!(molecules[molecule_id], xtal.box)

                # compute the potential energy of the molecule in its new configuraiton
                energy_new = potential_energy(molecule_id, molecules, xtal, ljff,
                                              eparams, eikr_gh, eikr_gg)

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
        end # inner cycles
    end # outer cycles
    # finished MC moves at this point.

    for m in molecules
        @assert inside(m) "molecule outside box!"
        @assert ! distortion(m, Frac(molecule, xtal.box), xtal.box, atol=1e-10) "molecules distorted"
    end

    # compute total energy, compare to `current_energy*` variables where were incremented
    system_energy_end = SystemPotentialEnergy()
    system_energy_end.gh.vdw = total_vdw_energy(xtal, molecules, ljff)
    system_energy_end.gg.vdw = total_vdw_energy(molecules, ljff, xtal.box)
    system_energy_end.gh.es = total(total_electrostatic_potential_energy(xtal, molecules,
                                                 eparams, eikr_gh))
    system_energy_end.gg.es = total(total_electrostatic_potential_energy(molecules,
                                        eparams, xtal.box, eikr_gg))

    # see Energetics_Util.jl for this function, overloaded isapprox to print mismatch
    if ! isapprox(system_energy, system_energy_end, verbose=true, atol=0.01)
        error("energy incremented improperly during simulation...")
    end

    @assert (markov_chain_time == sum(markov_counts.n_proposed))

    # return molecules in Cartesian format
    molecules = Cart.(molecules, xtal.box)

    return molecules # summary of statistics and ending configuration of molecules
end # μVT_sim
