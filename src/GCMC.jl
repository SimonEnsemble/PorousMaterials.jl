#Arbitrary value for delta for now, can change or make more complex later
const δ = 0.1


type GCMCstats
    n::Int
    n²::Int
    U_gh::Float64
    U_gh²::Float64
    U_gg::Float64
    U_gg²::Float64
    U_ggU_gh::Float64
    num_samples::Int
end

type MarkovCounts
    num_moves::Int

    insert_proposal::Int
    insert_acceptance::Int

    translate_proposal::Int
    translate_acceptance::Int

    delete_proposal::Int
    delete_acceptance::Int
end

"""
    insert_molecule!(molecules, simulation_box, adsorbate)

Inserts an additional molecule into the system and then checks the metropolis
hastings acceptance rules to find whether the proposed insertion is accepted or
rejected. Function then returns the new list of molecules
"""
function insert_molecule!(molecules::Array{Molecule}, simulation_box::Box,
        adsorbate::String)
    #TODO create template to handle more complex molecules
    x_new = simulation_box.f_to_c * [rand(), rand(), rand()]
#    new_molecule = Molecule(1, [adsorbate], x_new, [0.0])
    new_molecule = Molecule(1, [adsorbate], x_new[:, :], [0.0])
    push!(molecules, new_molecule)
end

"""
    delete_molecule!(molecule_id, molecules)

Removes a random molecule from the current molecules in the framework.
molecule_id decides which molecule will be deleted, for a simulation, it must
    be a randomly generated value
"""
function delete_molecule!(molecule_id::Int, molecules::Array{Molecule})
    deleteat!(molecules, molecule_id)
end

"""
    translate_molecule(molecule_id, molecules, simulation_box)

Translates a random molecule a random amount in the given array and simulation
box
reflects if entire molecule goes outside in any cartesian direction
"""
function translate_molecule!(molecule_id::Int, molecules::Array{Molecule},
        simulation_box::Box)

    old_coords = molecules[molecule_id].x

    dx = [δ * (rand() - 0.5) for i = 1:3]
    molecules[molecule_id].x .+= dx
    xf_molecule = simulation_box.c_to_f * molecules[molecule_id].x
    for coords = 1:3
        if sum(xf_molecule[coords, :] .< 1.0) == 0 #if all atoms of the molecule
                    #are past the upper end of the simulation box
            xf_molecule[coords, :] -= 1.0
        elseif sum(xf_molecules[coords, :] .> 0.0) == 0 #if all atoms of the
                    #molecule are past the lower end of the simulation box
            xf_molecule[coords, :] += 1.0
        end #if statement that checks for reflection
    end #for loop to go over x, y, and, z coordinate systems

    return old_coords
end

"""
    proposal_energy = guest_guest_vdw_energy(molecule_id, molecules,
        ljforcefield, simulation_box)

Calculates energy of a single adsorbate in the system. This can be used to find
the change in energy after an accepted proposal

Code copied from Arni's Energetics.jl with minor adjustments to calculate
    interactions between the adsorbates as well as the framework
"""
function guest_guest_vdw_energy(molecule_id::Int, molecules::Array{Molecule},
                        ljforcefield::LennardJonesForceField,
                        simulation_box::Box)
    #start energy at 0 and add to it as comparisons are made
    energy = 0.0
    #Loop over all atoms in the given molecule
    for atom_id = 1:molecules[molecule_id].n_atoms
        xf_molecule_atom = mod.(simulation_box.c_to_f *
            molecules[molecule_id].x[:, atom_id], 1.0)

        #loop over other adsorbate atoms leading up to central atom in array
        for other_molecule_id = 1:length(molecules)
            #molecule cannot interact with itself
            if other_molecule_id == molecule_id
                continue
            end #if statement prevents molecule from interacting with itself
            #loop over every atom in second molecules
            for other_atom_id = 1:molecules[other_molecule_id].n_atoms
                xf_other_molecule_atom = mod.(simulation_box.c_to_f *
                    molecules[other_molecule_id].x[:, other_atom_id], 1.0)
                dxf = xf_molecule_atom - xf_other_molecule_atom

                #runs nearest image convention and adjusts accordingly
                nearest_image!(dxf, 1.0)
                #converts fractional distance to cartesian distance
                dx = simulation_box.f_to_c * dxf

                r_squared = dot(dx, dx)
                if r_squared < R_OVERLAP_squared
                    #the atoms are overlapping, return inf energy
                    return Inf
                elseif r_squared < ljforcefield.cutoffradius_squared
                    #atoms can interact and won't produce Inf energy
                    #TODO test whether it is more efficient to store this as a variable up top
                    energy += lennard_jones(r_squared,
                        ljforcefield.sigmas_squared[molecules[other_molecule_id].atoms[other_atom_id]][molecules[molecule_id].atoms[atom_id]],
                        ljforcefield.epsilons[molecules[other_molecule_id].atoms[other_atom_id]][molecules[molecule_id].atoms[atom_id]])
                end #if statement for whether energy will be calculated
            end #for loop to check all atoms in other molecule
        end #for loop to go over all other molecules
    end #for loop for every atom in test molecule
    return energy #units are the same as in ϵ for forcefield (Kelvin)
end

"""
    gcmc_sim(framework, temperature, pressure, adsorbate, ljforcefield)

runs a monte carlo simulation using the given framework and adsorbate.
runs at the given temperature and pressure
"""
#will pass in molecules::Array{Molecule} later
function gcmc_sim(framework::Framework, temperature::Float64, pressure::Float64,
                        adsorbate::String, ljforcefield::LennardJonesForceField,
                        fugacity::Float64; num_mc_trials::Int=100000,
                        num_burn_cycles::Int=10000, sample_frequency::Int=25)

    # Boltmann constant (Pa-m3/K --> Pa-A3/K)
     # const KB = 1.38064852e-23 * 1e30
    const KB = 1.38064852e7
    const repfactors = replication_factors(framework.box, ljforcefield)
    const simulation_box = replicate_box(framework.box, repfactors)

    current_energy_gg = 0.0
    current_energy_gh = 0.0
    gcmc_stats = GCMCstats(0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0)
    markov_counts = MarkovCounts([0 for i = 0:7]...)

    molecules = Molecule[]
"""
#will need this later for when starting with a list of molecules
    for molecule_id = 1:length(molecules)
        current_energy += guest_guest_vdw_energy(molecule_id, molecules,
            ljforcefield, framework.box)/2
        current_energy += vdw_energy(framework, molecules[molecule_id],
            ljforcefield, repfactors)
    end
"""

    for t = 1:num_mc_trials
        markov_counts.num_moves += 1
        proposal_id = rand(1:3)
        if proposal_id == 1 #insert
            markov_counts.insert_proposal += 1
            insert_molecule!(molecules, framework.box)
            U_gg = guest_guest_vdw_energy(length(molecules), molecules,
                ljforcefield, simulation_box)
            U_gh = vdw_energy(framework, molecules[end],
                ljforcefield, repfactors)
            #Metropolis Hastings Acceptance for Insertion
            if rand() < fugacity * sim_box.Ω / (length(molecules) * KB *
                    temperature) * exp(-(U_gh + U_gg) / temperature)
                #accept the move, adjust current_energy
                markov_counts.insert_acceptance += 1
                current_energy_gg += U_gg
                current_energy_gh += U_gh
            else
                #reject the move, remove the added molecule
                pop!(molecules)
            end
        elseif proposal_id == 2 #delete
            markov_counts.delete_proposal += 1
            molecule_id = rand(1:length(molecules))
            U_gg = guest_guest_vdw_energy(molecule_id, molecules, ljforcefield,
                simulation_box)
            U_gh = vdw_energy(framework, molecules[molecule_id], ljforcefield,
                repfactors)
            #Metropolis Hastings Acceptance for Deletion
            if rand() < fugacity * sim_box.Ω / (length(molecules) * KB *
                    temperature) * exp((U_gh + U_gg) / temperature)
                #accept the deletion, delete molecule, adjust current_energy
                markov_counts.delete_acceptance += 1
                delete_molecule(molecule_id, molecules)
                current_energy_gg += U_gg
                current_energy_gh += U_gh
            end

        else #translate
            markov_counts.translate_proposal += 1
            molecule_id = rand(1:length(molecules))
            #energy of the molecule before it was translated
            U_gg_old = guest_guest_vdw_energy(molecule_id, molecules,
                ljforcefield, simulation_box)
            U_gh_old = vdw_energy(framework, molecules[molecule_id],
                ljforcefield, repfactors)

            old_molecule = translate_molecule!(molecule_id, molecules,
                simulation_box)

            #energy of the molecule after it is translated
            U_gg_new = guest_guest_vdw_energy(molecule_id, molecules,
                ljforcefield, simulation_box)
            U_gh_new = vdw_energy(framework, molecules[molecule_id],
                ljforcefield, rep_factors)

            #Metropolis Hastings Acceptance for translation
            if rand() < exp(-((U_gg_new + U_gh_new) - (U_gg_old + U_gh_old))
                / temperature)
                #accept the move, adjust current energy
                markov_counts.translate_acceptance += 1
                current_energy_gg += U_gg_new - U_gg_old
                current_energy_gh += U_gh_new - U_gh_old
            else
                #reject the move, reset the molecule at molecule_id
                molecules[molecule_id] = old_molecule
            end

        end #which move the code executes
        #sample the current state space
        if(t > num_burn_cycles && t % sample_frequency == 0)
            gcmc_stats.n += length(molecules)
            gcmc_stats.n² += length(molecules) ^ 2
            gcmc_stats.U_gh += current_energy_gh
            gcmc_stats.U_gh² += current_energy_gh ^ 2
            gcmc_stats.U_gg += current_energy_gg
            gcmc_stats.U_gg² += current_energy_gg ^ 2
            gcmc_stats.U_ggU_gh += current_energy_gg * current_energy_gh
            gcmc_stats.num_samples += 1
        end

    end #finished markov chain proposal moves

    results = Dict{String, Union{Int, Float64}}()
    results["⟨N⟩ (molecules)"] = gcmc_stats.n / gcmc_stats.num_samples
    results["⟨N⟩ (molecules/unit cell)"] = results["⟨N⟩ (molecules)"] /
        (repfactors[1] * repfactors[2] * repfactors[3])
    #(molecules/unit cell) * (mol/6.02 * 10^23 molecules) * (1000 mmol/mol) *
    #    (unit cell/framework amu) * (amu/ 1.66054 * 10^-24)
    results["⟨N⟩ (mmol/g)"] = results["⟨N⟩ (molecules/unit cell)"] * 1000 /
        (6.022140857e23 * molecular_weight(framework) * 1.66054e-24)
    results["⟨U_gg⟩ (K)"] = gcmc_stats.U_gg / gcmc_stats.num_samples
    results["⟨U_gh⟩ (K)"] = gcmc_stats.U_gh / gcmc_stats.num_samples
    results["⟨Energy⟩ (K)"] = (gcmc_stats.U_gg + gcmc_stats.U_gh) /
        gcmc_stats.num_samples
    #variances
    results["var(N)"] = (gcmc_stats.n² / gcmc_stats.num_samples) -
        (results["⟨N⟩ (molecules)"] ^ 2)
    results["var(U_gg)"] = (gcmc_stats.U_gg² / gcmc_stats.num_samples) -
        (results["⟨U_gg⟩ (K)"] ^ 2)
    results["var⟨U_gh⟩"] = (gcmc_stats.U_gh² / gcmc_stats.num_samples) -
        (results["⟨U_gh⟩ (K)"] ^ 2)
    results["var(Energy)"] = ((gcmc_stats.U_gg² + gcmc_stats.U_gh² + 2 *
        gcmc_stats.U_ggU_gh) / gcmc_stats.num_samples) -
        (results["⟨Energy⟩ (K)"] ^ 2)
    #insertions stats
    results["Number of attempted Insertions"] = markov_counts.insert_proposal
    results["Percentage of accepted Insertions"] =
        markov_counts.insert_acceptance / markov_counts.insert_proposal
    #deletion stats
    results["Number of attempted Deletions"] = markov_counts.delete_proposal
    results["Percentage of accepted Deletions"] =
        markov_counts.delete_acceptance / markov_counts.delete_proposal
    #translation stats
    results["Number of attempted Translations"] = markov_counts.translate_proposal
    results["Percentage of accepted Translations"] =
        markov_counts.translate_acceptance / markov_counts.translate_proposal

    return results

end #gcmc_sim
