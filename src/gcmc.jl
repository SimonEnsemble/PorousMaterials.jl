import Base: +, /

const KB = 1.38064852e7 # Boltmann constant (Pa-m3/K --> Pa-A3/K)

###
#   Markov chain proposals
###
const PROPOSAL_ENCODINGS = Dict(1 => "insertion", 
                                2 => "deletion",
                                3 => "translation",
                                4 => "rotation",
                                5 => "reinsertion"
                                ) # helps with printing later
const N_PROPOSAL_TYPES = length(keys(PROPOSAL_ENCODINGS))
# each proposal type gets an Int for clearer code
const INSERTION   = Dict([v => k for (k, v) in PROPOSAL_ENCODINGS])["insertion"]
const DELETION    = Dict([v => k for (k, v) in PROPOSAL_ENCODINGS])["deletion"]
const TRANSLATION = Dict([v => k for (k, v) in PROPOSAL_ENCODINGS])["translation"]
const ROTATION    = Dict([v => k for (k, v) in PROPOSAL_ENCODINGS])["rotation"]
const REINSERTION = Dict([v => k for (k, v) in PROPOSAL_ENCODINGS])["reinsertion"]

 # count proposed/accepted for each subtype
mutable struct MarkovCounts
    n_proposed::Array{Int, 1}
    n_accepted::Array{Int, 1}
end

###
#   collecting statistics
###
mutable struct GCMCstats
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

function Base.print(gcmc_stats::GCMCstats)
    println("\t# samples: ", gcmc_stats.n_samples)
    println("\t⟨N⟩ (molecules) = ", gcmc_stats.n / gcmc_stats.n_samples)

    println("\t⟨U_gh, vdw⟩ (K) = ",     gcmc_stats.U.gh.vdw / gcmc_stats.n_samples)
    println("\t⟨U_gh, Coulomb⟩ (K) = ", gcmc_stats.U.gh.es / gcmc_stats.n_samples)
    println("\t⟨U_gg, vdw⟩ (K) = ",     gcmc_stats.U.gg.vdw / gcmc_stats.n_samples)
    println("\t⟨U_gg, Coulomb⟩ (K) = ", gcmc_stats.U.gg.es / gcmc_stats.n_samples)

    println("\t⟨U⟩ (K) = ", sum(gcmc_stats.U) / gcmc_stats.n_samples)
end

# Compute average and standard error of the number of molecules and potential
# energy from an array of `GCMCstats`, each corresponding to statitics from an
# independent block/bin during the simulation. The average from each bin is
# treated as an independent sample and used to estimate the error in the
# estimate as a confidence interval.
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

# TODO move this to MC helpers? but not sure if it will inline. so wait after test with @time
# potential energy change after inserting/deleting/perturbing coordinates of molecules[which_species][molecule_id]
@inline function potential_energy(which_species::Int,
								  molecule_id::Int,
                                  molecules::Array{Array{Molecule{Frac}, 1}, 1},
                                  xtal::Crystal,
                                  ljff::LJForceField)
    energy = SystemPotentialEnergy()
    # van der Waals interactions
    energy.gg.vdw = vdw_energy(which_species, molecule_id, molecules, ljff, xtal.box) # guest-guest
    energy.gh.vdw = vdw_energy(xtal, molecules[which_species][molecule_id], ljff)     # guest-host
    # TODO electrostatic interactions
    return energy
end

function n_this_species(which_species::Int, molecules::Array{Array{Molecule, 1}, 1}, molecular_species::Symbol)
	if length(molecules) == 0
		return 0
	end

	n_i = length(molecules[which_species])
end

"""
    results, molecules = μVT_sim(xtal, molecule_templates, temperature, pressure,
                                 ljff; molecules=Array{Molecule, 1}[], settings=settings)

Runs a grand-canonical (μVT) Monte Carlo simulation of the adsorption of a molecule in a
xtal at a particular temperature and pressure using a
Lennard Jones force field.

A cycle is defined as max(20, number of adsorbates currently in the system) Markov chain
proposals. Current Markov chain moves implemented are particle insertion/deletion and
translation.

# Arguments
- `xtal::Crystal`: the porous xtal in which we seek to simulate adsorption
- `molecule_templates::Array{Molecule, 1}`: an array of the templates of unique adsorbate molecules of which we seek to simulate
- `temperature::Float64`: temperature of bulk gas phase in equilibrium with adsorbed phase
    in the porous material. units: Kelvin (K)
- `pressures::Array{Float64, 1}`: pressure of bulk gas phase in equilibrium with adsorbed phase in the
    porous material for each adsorbate. units: bar
    the adsorption
- `ljff::LJForceField`: the molecular model used to describe the
- `molecules::Array{Array{Molecule{Cart}, 1}, 1}`: a starting configuration of molecules in the xtal with an array per species.
- `n_burn_cycles::Int`: number of cycles to allow the system to reach
    equilibrium before sampling.
- `n_sample_cycles::Int`: number of cycles used for sampling
- `sample_frequency::Int`: during the sampling cycles, sample e.g. the number of
 adsorbed gas molecules every this number of Markov proposals
- `verbose::Bool`: whether or not to print off information during the simulation
- `ewald_precision::Float64`: desired precision for the long range Ewald summation
- `eos::Symbol`: equation of state to use for calculation of fugacity from pressure
- `write_adsorbate_snapshots::Bool`: whether the simulation will create and save a snapshot file
- `snapshot_frequency::Int`: the number of cycles taken between each snapshot (after burn cycle completion)
- `calculate_density_grid::Bool`: whether the simulation will keep track of a density grid for adsorbates
- `density_grid_dx::Float64`: The (approximate) space between voxels (in Angstroms) in the density grid. The number of voxels in the simulation box is computed automatically by [`required_n_pts`](@ref).
- `density_grid_molecular_species::Symbol`: the adsorbate for which we will make a density grid of its position (center).
- `density_grid_sim_box::Bool`: `true` if we wish for the density grid to be over the 
entire simulation box as opposed to the box of the crystal passed in. `false` if we wish the
density grid to be over the original `xtal.box`, before replication, passed in.
- `results_filename_comment::AbstractString`: An optional comment that will be appended to the name of the saved file (if autosaved)
"""
function μVT_sim(xtal::Crystal, 
                 molecule_templates::Array{Molecule, 1},
                 temperature::Float64,
                 pressures::Array{Float64, 1}, 
                 ljff::LJForceField; 
                 molecules::Array{Array{Molecule{Cart}, 1}, 1}=Array{Molecule{Cart}, 1}[], 
                 n_burn_cycles::Int=5000,
                 n_sample_cycles::Int=5000,
                 sample_frequency::Int=1,
                 verbose::Bool=true,
                 ewald_precision::Float64=1e-6,
                 eos::Symbol=:ideal,
                 autosave::Bool=true,
                 show_progress_bar::Bool=false,
                 write_adsorbate_snapshots::Bool=false,
                 snapshot_frequency::Int=1,
                 calculate_density_grid::Bool=false,
                 density_grid_dx::Float64=1.0,
                 density_grid_molecular_species::Union{Nothing, Symbol}=nothing,
                 density_grid_sim_box::Bool=true,
                 results_filename_comment::String=""
                )
    assert_P1_symmetry(xtal)

    start_time = time()
    # # to avoid changing the outside object `molecule_` inside this function, we make
    # #  a deep copy of it here. this serves as a template to copy when we insert a new molecule.
    # molecule = deepcopy(molecule_)

    nb_species = length(molecule_templates)
    molecular_species = [mt.species for mt in molecule_templates]

    if nb_species != length(pressures)
        error("# molecules in simulation: $nb_species
        # partial pressures: $(length(pressures))
        these should be equal!")
    end

    if verbose
        pretty_print(xtal, molecule_templates, temperature, pressures, ljff)
        println("\t# burn cycles: ", n_burn_cycles)
        println("\t# sample cycles: ", n_sample_cycles)
    end
    
    ###
    #  xyz file for storing snapshots of adsorbate positions
    ###
    num_snapshots = 0
    xyz_snapshots_filename = μVT_output_filename(xtal, molecule_templates, temperature, pressures, ljff, 
                                                 n_burn_cycles, n_sample_cycles, extension=".xyz")
    xyz_snapshot_file = IOStream(xyz_snapshots_filename) # declare a variable outside of scope so we only open a file if we want to snapshot
    if write_adsorbate_snapshots
        xyz_snapshot_file = open(xyz_snapshot_filename, "w")
    end

    ###
    #  Convert pressure to fugacity (units: Pascal) using an equation of state
    ###
    # TODO make PengRobinson work for multiple species (ignore for now)
    fugacities = [NaN for s = 1:nb_species] # Pa
    if eos == :ideal
       fugacities = pressures * 100000.0 # bar --> Pa
    elseif eos == :PengRobinson
        if nb_species != 1
            error("Peng Robinson EOS not implemented for >1 species")
        end
        prfluid = PengRobinsonFluid(molecule_templates[1].species)
        gas_props = calculate_properties(prfluid, temperature, pressures[1], verbose=false)
        fugacities[1] = gas_props["fugacity (bar)"] * 100000.0 # bar --> Pa
    else
        error("eos=:ideal and eos=:PengRobinson are the only valid options for an equation of state.")
    end
    if verbose
        for s = 1:nb_species
            @printf("\t%s equation of state, %s partial fugacity = %f bar\n", eos, 
                    molecule_templates[s], fugacities[s] / 100000.0)
        end
    end
    
    ###
    #   replicate xtal so that nearest image convention can be applied for short-range interactions
    ###
    repfactors = replication_factors(xtal.box, ljff)
    original_xtal_box = deepcopy(xtal.box)
    xtal = replicate(xtal, repfactors) # frac coords still in [0, 1]
    
    # convert molecules array to fractional using this box.
    molecules = [Frac.(mols, xtal.box) for mols in molecules]

    ###
    #   Density grid for adsorbate
    #   (if more than one adsorbate, user must specify which adsorbate species to make
    #     the grid for)
    ###
    if calculate_density_grid && isnothing(density_grid_molecular_species)
        if nb_species == 1
            # obviously we are keeping track of the only atom in the adsorbate.
            density_grid_molecular_species = molecule_templates[1].species
        else
            # cannot proceed if we do not know which molecule to keep track of!
            error(@printf("Passed `calculate_density_grid=true` but there is %d
                 different adsorbates in the system. Pass `density_grid_species` to specify
                 which adsorbate to make a grid for.", nb_species)
                 )
        end
    end
    # Initialize a density grid based on the *simulation box* (not xtal box passed in) and the passed in density_grid_dx
    # Calculate `n_pts`, number of voxels in grid, based on the sim box and specified voxel spacing
    n_pts = (0, 0, 0) # don't store a huge grid if we aren't tracking a density grid
    if calculate_density_grid
        if density_grid_sim_box
            n_pts = required_n_pts(xtal.box, density_grid_dx)
        else
            n_pts = required_n_pts(original_xtal_box, density_grid_dx)
        end
    end
    density_grid = Grid(density_grid_sim_box ? xtal.box : original_xtal_box, n_pts, 
        zeros(n_pts...), :inverse_A3, [0.0, 0.0, 0.0])

    if verbose
        println("\tthe crystal:")
        @printf("\t\treplicated (%d,%d,%d) for short-range cutoff of %f Å\n",
                repfactors[1], repfactors[2], repfactors[3],
                sqrt(ljff.r²_cutoff))
        println("\t\tdensity [kg/m³]: ", crystal_density(xtal))
        println("\t\tchemical formula: ", chemical_formula(xtal))
        println("\t\t# atoms: ", xtal.atoms.n)
        println("\t\t# point charges: ", xtal.charges.n)
        println("\tthe molecules:")
        for s = 1:nb_species
            println("\t\tunique species: ", unique(molecule_templates[s].atoms.species))
            println("\t\t# atoms: ", molecule_templates[s].atoms.n)
            println("\t\t# point charges: ", molecule_templates[s].charges.n)
        end
        if write_adsorbate_snapshots
            @printf("\tWriting snapshots of adsorption positions every %d cycles (after burn cycles)\n", snapshot_frequency)
            @printf("\t\tWriting to file: %s\n", xyz_filename)
        end
        if calculate_density_grid
            @printf("\tTracking adsorbate spatial probability density grid of adsorbate %s, updated every %d cycles (after burn cycles)\n", density_grid_molecular_species, snapshot_frequency)
            @printf("\t\tdensity grid voxel spacing specified as %.3f Å => %d by %d by %d voxels\n", density_grid_dx, n_pts...)
            if density_grid_sim_box
                @printf("\t\tdensity grid is over simulation box\n")
            else
                @printf("\t\tdensity grid is over original crystal box\n")
            end
        end
    end
    if ! forcefield_coverage(xtal.atoms, ljff)
        error("crystal $(xtal.name) not covered by the force field $(ljff.name)")
    end
    
    for s = 1:nb_species
        if ! neutral(molecule_templates[s].charges)
            error(@sprintf("Molecule %s is not charge neutral!\n", molecule_templates[s].species))
        end
        if ! forcefield_coverage(molecule_templates[s].atoms, ljff)
            error("molecule $(molecule_templates[s].species) not covered by the force field $(ljff.name)")
        end
    end

    # TODO electrostatics 
    electrostatics_flag = has_charges(xtal) && any(has_charges.(molecule_templates))
    if electrostatics_flag
        error("sorry, electrostatics not supported")
    end

    # initiate system energy to which we increment when MC moves are accepted
    system_energy = SystemPotentialEnergy()
    # if we don't start with an emtpy xtal, compute energy of starting configuration
    #  (n=0 corresponds to zero energy)
    if length(molecules) != 0
        # some checks
#        for m in molecules 
#            # ensure molecule template matches species of starting molecules.
#            @assert m.species in [mt.species for mt in molecule_templates] "initializing with wrong molecule species"
#            # assert that the molecules are inside the simulation box
#            @assert inside(m) "initializing with molecules outside simulation box!"
#            # ensure pair-wise bond distances match template
#            id_mt = findfirst([mt.species for mt in molecule_templates] .== m.species)
#            @assert ! distortion(m, Frac(molecule_template[id_mt], xtal.box), xtal.box)
#        end

		for m in molecules # these are arrays of molecule structs
			# ensuer molecule tempale matches species of starting molecules.
			@assert all(i -> (i.species in [mt.species for mt in molecule_templates]), m) "initializing with wrong molecule species"
			# assert that the molecules are inside the simulation box
			@assert all(i -> inside(i), m) "initializing with molecules outside simulation box!"
			# ensure pair-wise bond distance match template
			ids_mt = [findfirst([mt.species for mt in molecule_templates] .== i.species) for i in m] # collection of ids
			@assert all(id -> [! distortion(i, Frac(molecule_template[id], xtal.box)) for i in m], ids_mt) "initializing with distorted molecules"

		end

        system_energy.gh.vdw = total_vdw_energy(xtal, molecules, ljff)
        system_energy.gg.vdw = total_vdw_energy(molecules, ljff, xtal.box)
        system_energy.gh.es = 0.0 #total(total_electrostatic_potential_energy(xtal, molecules, eparams, eikr_gh))
        system_energy.gg.es = 0.0 #total(electrostatic_potential_energy(molecules, eparams, xtal.box, eikr_gg))
    end

    if show_progress_bar
        progress_bar = Progress(n_burn_cycles + n_sample_cycles, 1)
    end

    ####
    #   proposal probabilities
    ###
    mc_proposal_probabilities = [[0.0 for p = 1:N_PROPOSAL_TYPES] for sp = 1:nb_species]
    # set defaults
    for s in 1:nb_species
		mc_proposal_probabilities[s][INSERTION] = 0.35
    	mc_proposal_probabilities[s][DELETION] = mc_proposal_probabilities[s][INSERTION] # must be equal
    	mc_proposal_probabilities[s][REINSERTION] = 0.05
    	if needs_rotations(molecule_templates[s])
        	mc_proposal_probabilities[s][TRANSLATION] = 0.125
        	mc_proposal_probabilities[s][ROTATION] = 0.125
    	else
       		mc_proposal_probabilities[s][TRANSLATION] = 0.25
        	mc_proposal_probabilities[s][ROTATION] = 0.0
    	end
		# Are these normalized across all species?
    	mc_proposal_probabilities[s] /= sum(mc_proposal_probabilities[s]) # normalize
    	# StatsBase.jl functionality for sampling
    	mc_proposal_probabilities[s] = ProbabilityWeights(mc_proposal_probabilities[s])
    	if verbose
        	println("\tMarkov chain proposals:")
        	for p = 1:N_PROPOSAL_TYPES
            	@printf("For species %s\n\t\tprobability of %s: %f\n", molecule_templates[s].species, 
						PROPOSAL_ENCODINGS[p], mc_proposal_probabilities[p])
        	end
    	end
	end

    # initiate GCMC statistics for each block # break simulation into `N_BLOCKS` blocks to gauge convergence
    gcmc_stats = [GCMCstats() for block_no = 1:N_BLOCKS]
    current_block = 1
    # make sure the number of sample cycles is at least equal to N_BLOCKS
    if n_sample_cycles < N_BLOCKS
        n_sample_cycles = N_BLOCKS
        @warn @sprintf("# sample cycles set to minimum %d, which is number of blocks.", N_BLOCKS)
    end
    N_CYCLES_PER_BLOCK = floor(Int, n_sample_cycles / N_BLOCKS)

    markov_counts = MarkovCounts(zeros(Int, length(PROPOSAL_ENCODINGS)), zeros(Int, length(PROPOSAL_ENCODINGS)))

    # (n_burn_cycles + n_sample_cycles) is number of outer cycles.
    #   for each outer cycle, peform max(20, # molecules in the system) MC proposals.
    markov_chain_time = 0
    outer_cycle_start = 1
    for outer_cycle = outer_cycle_start:(n_burn_cycles + n_sample_cycles)
        if show_progress_bar
            next!(progress_bar; showvalues=[(:cycle, outer_cycle), (:number_of_molecules, length(molecules))])
        end
        for inner_cycle = 1:max(20, length(molecules))
            markov_chain_time += 1
            
            # choose a species
            which_species = rand(1:nb_species)
			# number of that species in molecules array
			# TODO: write this function
			n_i = n_this_species(molecules, molecular_species[which_species])

			# choose proposed move randomly; keep track of proposals
            which_move = sample(1:N_PROPOSAL_TYPES, mc_proposal_probabilities[which_species]) # StatsBase.jl
            markov_counts.n_proposed[which_move] += 1

			# TODO: What to do about length(molecules) != 0 statements if the sub-array for which_species hasn't been created?
			# 			i.e. if a sub-array for one
            if which_move == INSERTION
				# TODO: overload for new molecules array
                random_insertion!(molecules[which_species], xtal.box, molecule_templates[which_species])

				# inserted molecule pushed to end of molecules[which_species] array
                ΔE = potential_energy(which_species, length(molecules[which_species]), molecules, xtal, ljff)

                # Metropolis Hastings Acceptance for Insertion
                if rand() < fugacities[which_species] * xtal.box.Ω / (n_i * KB *
                        temperature) * exp(-sum(ΔE) / temperature)
                    # accept the move, adjust current_energy
                    markov_counts.n_accepted[which_move] += 1

                    system_energy += ΔE
                else
                    # reject the move, remove the inserted molecule
                    pop!(molecules)
                end
            elseif (which_move == DELETION) && (length(molecules) != 0)
                # get ids that belong to this species in the molecules array

				#	ids_this_species = [i_m if (molecules[which_species][i_m].species == molecular_species[which_species]) for i_m = 1:length(molecules[which_species])]

				ids_this_species = [i_m for i_m in 1:length(molecules[which_species])]
				@assert all(i_m -> molecules[which_species][i_m].species == moleculear_species[which_species], ids_this_species) "ids do not match chosen species in molecules"

				# propose which molecule to delete
                # molecule_id = rand(1:length(molecules))
				molecule_id = rand(ids_this_species)

                # compute the potential energy of the molecule we propose to delete
                ΔE = potential_energy(which_species, molecule_id, molecules, xtal, ljff)

                # Metropolis Hastings Acceptance for Deletion
                if rand() < length(molecules) * KB * temperature / (fugacity *
                        xtal.box.Ω) * exp(sum(ΔE) / temperature)
                    # accept the deletion, delete molecule, adjust current_energy
                    markov_counts.n_accepted[which_move] += 1

                    remove_molecule!(molecule_id, molecules) # overload for new molecules array

                    system_energy -= ΔE
                end
            elseif (which_move == TRANSLATION) && (length(molecules) != 0)
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
            elseif (which_move == ROTATION) && (length(molecules) != 0)
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
            elseif (which_move == REINSERTION) && (length(molecules) != 0)
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

            # if we've done all burn cycles, take samples for statistics
            if outer_cycle > n_burn_cycles
                if markov_chain_time % sample_frequency == 0
                    gcmc_stats[current_block].n_samples += 1

                    gcmc_stats[current_block].n += length(molecules)
                    gcmc_stats[current_block].n² += length(molecules) ^ 2

                    gcmc_stats[current_block].U += system_energy
                    gcmc_stats[current_block].U² += square(system_energy)

                    gcmc_stats[current_block].Un += sum(system_energy) * length(molecules)
                end
            end # sampling
        end # inner cycles

        # print block statistics / increment block
        if (outer_cycle > n_burn_cycles) && (current_block != N_BLOCKS) && (
            (outer_cycle - n_burn_cycles) % N_CYCLES_PER_BLOCK == 0)
            # move onto new block unless current_block is N_BLOCKS;
            # then just keep adding stats to the last block.
            # this only occurs if sample_cycles not divisible by N_BLOCKS
            # print GCMC stats later and do not increment block if we are in last block.
            # print statistics for this block
            if verbose
                printstyled(@sprintf("\tBlock  %d/%d statistics:\n", current_block, N_BLOCKS); color=:yellow)
                print(gcmc_stats[current_block])
            end
            current_block += 1
        end
        # print the last cycle in the last block
        if outer_cycle == (n_sample_cycles + n_burn_cycles)
            if verbose
                printstyled(@sprintf("\tBlock  %d/%d statistics:\n", current_block, N_BLOCKS); color=:yellow)
                print(gcmc_stats[current_block])
            end
        end

        # snapshot cycle
        if (outer_cycle > n_burn_cycles) && (outer_cycle % snapshot_frequency == 0)
            if write_adsorbate_snapshots
                # have a '\n' for every new set of atoms, leaves no '\n' at EOF
                if num_snapshots > 0
                    @printf(xyz_snapshot_file, "\n")
                end
                write_xyz(xtal.box, molecules, xyz_snapshot_file)
            end
            if calculate_density_grid
                if density_grid_sim_box
                    update_density!(density_grid, molecules, density_grid_species)
                else
                    update_density!(density_grid, Cart.(molecules, xtal.box), density_grid_species)
                end
            end
            num_snapshots += 1
        end
    end # outer cycles
    # finished MC moves at this point.

    # close snapshot xyz file
    close(xyz_snapshot_file)

    if calculate_density_grid
        # divide number of molecules in a given voxel by total snapshots
        density_grid.data ./= num_snapshots
    end

    # checks
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
    elapsed_time = time() - start_time
    if verbose
        @printf("\tEstimated elapsed time: %d seconds\n", elapsed_time)
        println("\tTotal # MC steps: ", markov_chain_time)
    end

    # build dictionary containing summary of simulation results for easy querying
    results = Dict{String, Any}()
    results["xtal"] = xtal.name
    results["adsorbate"] = molecule.species
    results["forcefield"] = ljff.name
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
    results["⟨U_gh, vdw⟩ (K)"]     = avg_U.gh.vdw
    results["⟨U_gh, electro⟩ (K)"] = avg_U.gh.es
    results["⟨U_gg, vdw⟩ (K)"]     = avg_U.gg.vdw
    results["⟨U_gg, electro⟩ (K)"] = avg_U.gg.es
    results["⟨U⟩ (K)"] = sum(avg_U)

    # variances
    results["var(N)"] = (sum(gcmc_stats).n² / sum(gcmc_stats).n_samples) -
        (results["⟨N⟩ (molecules)"] ^ 2)
    # isosteric heat of adsorption TODO stdev of this too.
    results["Q_st (K)"] = temperature - (sum(gcmc_stats).Un / sum(gcmc_stats).n_samples - results["⟨U⟩ (K)"] * results["⟨N⟩ (molecules)"]) / results["var(N)"]

    # error bars (confidence intervals)
    results["err ⟨N⟩ (molecules)"]     = err_n
    results["err ⟨U_gh, vdw⟩ (K)"]     = err_U.gh.vdw
    results["err ⟨U_gh, electro⟩ (K)"] = err_U.gh.es
    results["err ⟨U_gg, vdw⟩ (K)"]     = err_U.gg.vdw
    results["err ⟨U_gg, electro⟩ (K)"] = err_U.gg.es
    results["err ⟨U⟩ (K)"] = sum(err_U)


    # average N in more common units
    results["⟨N⟩ (molecules/unit cell)"] = avg_n / (repfactors[1] * repfactors[2] * repfactors[3])
    results["err ⟨N⟩ (molecules/unit cell)"] = err_n / (repfactors[1] * repfactors[2] * repfactors[3])
    # (molecules/unit cell) * (mol/6.02 * 10^23 molecules) * (1000 mmol/mol) *
    #    (unit cell/xtal amu) * (amu/ 1.66054 * 10^-24)
    results["⟨N⟩ (mmol/g)"] = results["⟨N⟩ (molecules/unit cell)"] * 1000 /
        (6.022140857e23 * molecular_weight(xtal) * 1.66054e-24) * (repfactors[1] * repfactors[2] * repfactors[3])
    results["err ⟨N⟩ (mmol/g)"] = results["err ⟨N⟩ (molecules/unit cell)"] * 1000 /
        (6.022140857e23 * molecular_weight(xtal) * 1.66054e-24) * (repfactors[1] * repfactors[2] * repfactors[3])

    # Markov stats
    for (proposal_id, proposal_description) in PROPOSAL_ENCODINGS
        results[@sprintf("Total # %s proposals", proposal_description)] = markov_counts.n_proposed[proposal_id]
        results[@sprintf("Fraction of %s proposals accepted", proposal_description)] = markov_counts.n_accepted[proposal_id] / markov_counts.n_proposed[proposal_id]
    end

    # Snapshot information
    results["density grid"] = deepcopy(density_grid)
    results["num snapshots"] = num_snapshots

    if verbose
        print_results(results, print_title=false)
    end

    # return molecules in Cartesian format
    molecules = Cart.(molecules, xtal.box)

    if autosave
        if ! isdir(PATH_TO_SIMS)
            mkdir(PATH_TO_SIMS)
        end

        save_results_filename = joinpath(PATH_TO_SIMS,
                                         μVT_output_filename(xtal, molecule_templates, temperature, pressures, 
                                                             ljff, n_burn_cycles, n_sample_cycles, 
                                                             comment=results_filename_comment
                                                             )
                                        )

        @save save_results_filename results
        if verbose
            println("\tresults dictionary saved in ", save_results_filename)
        end
    end

    return results, molecules # summary of statistics and ending configuration of molecules
end # μVT_sim

"""
    filename = μVT_output_filename(xtal, molecule_templates, temperature, 
                                   pressures, ljff, n_burn_cycles, 
                                   n_sample_cycles; comment="", extension=".jld2")

This is the function that establishes the file naming convention used by [μVT_sim](@ref).

# Arguments
- `xtal::Crystal`: porous xtal used in adsorption simulation
- `molecule_templates::Array{Molecule, 1}`: template of the adsorbate molecules used in adsorption simulation
- `temperature::Float64`:temperature of bulk gas phase in equilibrium with adsorbed phase in
        the porous material. units: Kelvin (K) 
- `pressures::Array{Float64, 1}`: partial pressures of bulk gas phase in equilibrium with adsorbed phase in the
        porous material. units: bar
- `ljff::LJForceField`: the molecular model used in adsorption simulation
- `n_burn_cycles::Int`: number of cycles to allow the system to reach equilibrium before sampling.    
- `n_sample_cycles::Int`: number of cycles used for sampling
- `comment::String=""`: remarks to be included in the filename
- `extension::String=".jld2"`: the file extension

# Returns
- `filename::String`: the name of the specific `.jld2` simulation file
"""
function μVT_output_filename(xtal::Crystal, molecule_templates::Array{Molecule, 1}, temperature::Float64,
                             pressures::Array{Float64, 1}, ljff::LJForceField, n_burn_cycles::Int, 
                             n_sample_cycles::Int; comment::String="", extension::String=".jld2")
    filename = @sprintf("muVT_xtal_%s_T_%.3fK", xtal.name, temperature)
    for m = 1:length(molecule_templates)
        filename *= @sprintf("_%s_P_%.6fbar", molecule_templates[m].species, pressures[m])
    end
    filename *= @sprintf("_%s_%dburn_%dsample", ljff.name, n_burn_cycles, n_sample_cycles)
    return filename * comment * extension
end

function print_results(results::Dict; print_title::Bool=true)
    if print_title
        # already print in GCMC tests...
        @printf("GCMC simulation of %s in %s at %f K and %f bar pressure, %f bar fugacity using %s forcefield.\n\n",
                results["adsorbate"], results["xtal"], results["temperature (K)"],
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
            printstyled(proposal_description; color=:yellow)
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

function pretty_print(xtal::Crystal, 
                      molecule_templates::Array{Molecule, 1}, 
                      temperature::Float64, 
                      pressures::Array{Float64, 1},
                      ljff::LJForceField)
    printstyled("μVT simulation\n"; color=:yellow)

    print("crystal: ")
    printstyled(xtal.name; color=:green)

    printstyled(@sprintf("temperature = %f K", temperature); color=:green)

    println("adsorbates:")
    for m = 1:length(molecule_templates)
        printstyled(molecule_templates[m].species; color=:green)
        printstyled(@sprintf("\tpartial pressure = %f bar", pressures[m]); color=:green)
    end

    printstyled(split(ljff.name, ".")[1]; color=:green)
    println(" force field.")
end

"""
    results = stepwise_adsorption_isotherm(xtal, molecule_templates, temperature, pressures,
                                  ljff; n_sample_cycles=5000,
                                  n_burn_cycles=5000, sample_frequency=1,
                                  verbose=true, ewald_precision=1e-6, eos=:ideal,
                                  write_adsorbate_snapshots=false,
                                  snapshot_frequency=1, calculate_density_grid=false,
                                  density_grid_dx=1.0, density_grid_species=nothing,
                                  density_grid_sim_box::Bool=true,
                                  results_filename_comment="", show_progress_bar=false)

Run a set of grand-canonical (μVT) Monte Carlo simulations in series. Arguments are the
same as [`μVT_sim`](@ref), as this is the function run behind the scenes. An
exception is that we pass an array of pressures. The adsorption isotherm is computed step-
wise, where the ending configuration from the previous simulation (array of molecules) is
passed into the next simulation as a starting point. The ordering of `pressures` is
honored. By giving each simulation a good starting point, (if the next pressure does not
differ significantly from the previous pressure), we can reduce the number of burn cycles
required to reach equilibrium in the Monte Carlo simulation. Also see
[`adsorption_isotherm`](@ref) which runs the μVT simulation at each pressure in parallel.
"""
function stepwise_adsorption_isotherm(xtal::Crystal,
                                      molecule_templates::Array{Molecule, 1},
                                      temperature::Float64,
                                      pressures::Array{Float64, 1},
                                      ljff::LJForceField;
                                      n_burn_cycles::Int=5000, n_sample_cycles::Int=5000,
                                      sample_frequency::Int=1, verbose::Bool=true,
                                      ewald_precision::Float64=1e-6, eos::Symbol=:ideal,
                                      show_progress_bar::Bool=false,
                                      write_adsorbate_snapshots::Bool=false,
                                      snapshot_frequency::Int=1, 
                                      calculate_density_grid::Bool=false,
                                      density_grid_dx::Float64=1.0,
                                      density_grid_species::Union{Nothing, Symbol}=nothing,
                                      density_grid_sim_box::Bool=true,
                                      results_filename_comment::AbstractString="")

    # simulation only works if xtal is in P1
    assert_P1_symmetry(xtal)

    results = Dict{String, Any}[] # push results to this array
    molecules = Array{Molecule{Cart}, 1}[] # initiate with empty xtal
    for (i, pressure) in enumerate(pressures)
        result, molecules = μVT_sim(xtal, molecule_templates, temperature, pressure,
                                            ljff,
                                            n_burn_cycles=n_burn_cycles,
                                            n_sample_cycles=n_sample_cycles,
                                            sample_frequency=sample_frequency,
                                            verbose=verbose, molecules=molecules, # essential step here
                                            ewald_precision=ewald_precision, eos=eos,
                                            write_adsorbate_snapshots=write_adsorbate_snapshots,
                                            snapshot_frequency=snapshot_frequency,
                                            calculate_density_grid=calculate_density_grid,
                                            density_grid_dx=density_grid_dx,
                                            density_grid_species=density_grid_species,
                                            density_grid_sim_box=density_grid_sim_box,
                                            results_filename_comment=results_filename_comment)
        push!(results, result)
    end

    return results
end


"""
    results = adsorption_isotherm(xtal, molecule_templates, temperature, 
                                  partial_pressures, sim_pressures,
                                  ljff; n_sample_cycles=5000,
                                  n_burn_cycles=5000, sample_frequency=1,
                                  verbose=true, ewald_precision=1e-6, eos=:ideal,
                                  write_adsorbate_snapshots=false,
                                  snapshot_frequency=1, calculate_density_grid=false,
                                  density_grid_dx=1.0, density_grid_species=nothing,
                                  density_grid_sim_box=true,
                                  results_filename_comment="", show_progress_bar=false)

Run a set of grand-canonical (μVT) Monte Carlo simulations in parallel. Arguments are the
same as [`μVT_sim`](@ref), as this is the function run in parallel behind the scenes.
The only exception is that we pass an array of pressures, and we only consider a single species.
To give Julia access to multiple cores, run your script as `julia -p 4 mysim.jl` to allocate e.g. four cores. See
[Parallel Computing](https://docs.julialang.org/en/stable/manual/parallel-computing/#Parallel-Computing-1).
"""
function adsorption_isotherm(xtal::Crystal,
                             molecule_templates::Array{Molecule, 1},
                             temperature::Float64,
                             pressures::Array{Float64, 1},
                             ljff::LJForceField;
                             n_burn_cycles::Int=5000, n_sample_cycles::Int=5000,
                             sample_frequency::Int=1, verbose::Bool=true,
                             ewald_precision::Float64=1e-6, eos::Symbol=:ideal,
                             show_progress_bar::Bool=false,
                             write_adsorbate_snapshots::Bool=false,
                             snapshot_frequency::Int=1, calculate_density_grid::Bool=false,
                             density_grid_dx::Float64=1.0,
                             density_grid_species::Union{Nothing, Symbol}=nothing,
                             density_grid_sim_box::Bool=true,
                             results_filename_comment::AbstractString="")

    # simulation only works if xtal is in P1
    assert_P1_symmetry(xtal)

    # we only consider a pure gas isotherm
    if length(molecules) != 1
        error("We only calculate pure gas isotherms")
    end

    # make a function of pressure only to facilitate uses of `pmap`
    run_pressure(partial_pressures::Array{Float64, 1}) = μVT_sim(xtal, molecule_templates, temperature,
                                                                 partial_pressures, ljff,
                                                                 n_burn_cycles=n_burn_cycles,
                                                                 n_sample_cycles=n_sample_cycles,
                                                                 sample_frequency=sample_frequency,
                                                                 verbose=verbose,
                                                                 ewald_precision=ewald_precision,
                                                                 eos=eos, 
                                                                 show_progress_bar=show_progress_bar,
                                                                 write_adsorbate_snapshots=write_adsorbate_snapshots,
                                                                 snapshot_frequency=snapshot_frequency,
                                                                 calculate_density_grid=calculate_density_grid,
                                                                 density_grid_dx=density_grid_dx,
                                                                 density_grid_species=density_grid_species,
                                                                 density_grid_sim_box=density_grid_sim_box,
                                                                 results_filename_comment=results_filename_comment)[1] # only return results

    # for load balancing, larger pressures with longer computation time goes first
    ids = sortperm(sim_pressures, rev=true)

    # define pressure = Array{Array{Float64, 1}, 1}
    #

    # run gcmc simulations in parallel using Julia's pmap parallel computing function
    results = pmap(run_pressure, pressures[ids])

    # return results in same order as original pressure even though we permuted them for
    #  better load balancing.
    return results[[findall(x -> x==i, ids)[1] for i = 1:length(ids)]]
end

"""
    dataframe = isotherm_sim_results_to_dataframe(desired_props, xtal,
                                                  molecule, temperature,
                                                  pressures, ljff,
                                                  n_burn_cycles, n_sample_cycles;
                                                  where_are_jld_files=nothing)`

convert the `.jld2` results output files auto-saved from [`μVT_sim`](@ref) into a `DataFrame`.
each row of the `DataFrame` corresponds to a pressure in the adsorption isotherm.
`desired_props` is an array of desired properties from the results. 
to locate the requested output files, this function calls [`μVT_output_filename`](@ref) to retrieve the file names.

# Arguments
- `desired_props::Array{String, 1}`: An array containing names of properties to be retrieved from the `.jld2` file
- `xtal::Crystal`: The porous crystal
- `molecule::Molecule{Cart}`: The adsorbate molecule
- `temperature::Float64`: The temperature in the simulation, units: Kelvin (K)
- `pressures::Array{Float64}`: The pressures in the adsorption isotherm simulation(s), units: bar
- `ljff::LJForceField`: The molecular model being used in the simulation to describe intermolecular Van
        der Waals interactions
- `n_burn_cycles::Int64`: The number of burn cycles used in this simulation
- `n_sample_cycles::Int64`: The number of sample cycles used in the simulations
- `where_are_jld_files::Union{String, Nothing}=nothing`: The location to the simulation files. defaults to
    `PorousMaterials.PATH_TO_SIMS`.
- `comment::String=""`: comment appended to outputfilename

# Returns
- `dataframe::DataFrame`: A `DataFrame` containing the simulated data along the adsorption isotherm, whose
    columns are for the specified properties

# Note
A range of pressures can be used to select a batch of simulation files to be included in the `DataFrame`.

# Example

```julia
xtal = Crystal("SBMOF-1.cif")
molecule = Molecule("Xe")
ljff = LJForceField("UFF", mixing_rules="Lorentz-Berthelot")
temperature = 298.0 # K
pressures = 10 .^ range(-2, stop=2, length=15)

dataframe = isotherm_sim_results_to_dataframe(["pressure (bar)", "⟨N⟩ (mmol/g)"], 
                                              xtal, molecule, temperature,
                                              pressures, ljff, 10000, 10000)
dataframe[Symbol("pressure (bar)")] # pressures
dataframe[Symbol("⟨N⟩ (mmol/g)")] # adsorption at corresponding pressures
```
"""
function isotherm_sim_results_to_dataframe(desired_props::Array{String, 1},
                                           xtal::Crystal, 
                                           molecule::Molecule,
                                           temperature::Float64,
                                           pressures::Array{Float64, 1},
                                           ljff::LJForceField,
                                           n_burn_cycles::Int64, 
                                           n_sample_cycles::Int64;
                                           comment::String="",
                                           where_are_jld_files::Union{String, Nothing}=nothing)
    # determine the location of the data files
    if isnothing(where_are_jld_files)
        where_are_jld_files = PATH_TO_SIMS
    end
    # prepare dataframe to populate
    df = DataFrame()
    # loop over pressures and populate dataframe
    for (i, pressure) in enumerate(pressures)
        jld2_filename = μVT_output_filename(xtal, molecule_templates, temperature, 
                                            pressures, ljff, n_burn_cycles, 
                                            n_sample_cycles, comment=comment)
        # load in the results as a dictionary
        @load joinpath(where_are_jld_files, jld2_filename) results
        
        # population column names, taking into account types
        if i == 1
            for col in desired_props
                insertcols!(df, length(names(df)) + 1, Symbol(col) => typeof(results[col])[])
            end
        end
        push!(df, [results[prop] for prop in desired_props])
    end
    return df
end
