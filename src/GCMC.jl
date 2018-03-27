# δ (units: Angstrom) is the maximal distance a particle is perturbed in a given coordinate
#  during particle translations
const δ = 0.35
const KB = 1.38064852e7 # Boltmann constant (Pa-m3/K --> Pa-A3/K)

# Markov chain proposals
const PROPOSAL_ENCODINGS = Dict(1 => "insertion", 2 => "deletion", 3 => "translation")
const N_PROPOSAL_TYPES = length(keys(PROPOSAL_ENCODINGS))
const INSERTION = Dict([value => key for (key, value) in PROPOSAL_ENCODINGS])["insertion"]
const DELETION = Dict([value => key for (key, value) in PROPOSAL_ENCODINGS])["deletion"]
const TRANSLATION = Dict([value => key for (key, value) in PROPOSAL_ENCODINGS])["translation"]

"""
Keep track of statistics during a grand-canonical Monte Carlo simultion

* `n` is the number of molecules
* `U` is the potential energy
* `g` refers to guest (the adsorbate molecule)
* `h` refers to host (the crystalline framework)
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
Monte Carlo simulation.
"""
type MarkovCounts
    n_proposed::Array{Int, 1}
    n_accepted::Array{Int, 1}
end

"""
    insert_molecule!(molecules, simulation_box, adsorbate)

Inserts an additional `adsorbate` molecule into the system at random coordinates inside 
the `simulation_box`.
"""
function insert_molecule!(molecules::Array{Molecule, 1}, box::Box, adsorbate::Symbol)
    #TODO create template to handle more complex molecules
    x = box.f_to_c * rand(3, 1)
    molecule = Molecule(1, [adsorbate], x[:, :], [0.0])
    push!(molecules, molecule)
end

"""
    delete_molecule!(molecule_id, molecules)

Removes a random molecule from the current molecules in the framework.
molecule_id decides which molecule will be deleted, for a simulation, it must
    be a randomly generated value
"""
function delete_molecule!(molecule_id::Int, molecules::Array{Molecule, 1})
    splice!(molecules, molecule_id)
end

"""
    bring_molecule_inside_box!(molecule, simulation_box)

Apply periodic boundary conditions to bring a molecule inside the a box if it is
completely outside of the box.
"""
function bring_molecule_inside_box!(molecule::Molecule, box::Box)
    outside_box = false # do nothing if not outside the box

    # compute its fractional coordinates
    xf = box.c_to_f * molecule.x

    # apply periodic boundary conditions if any of its x, y, z fractional coords are
    #  greater than 1.0
    for xyz = 1:3 # loop over x, y, z coordinates
        # if all atoms of the molecule have x, y, or z > 1.0, shift down
        if sum(xf[xyz, :] .<= 1.0) == 0
            outside_box = true
            xf[xyz, :] -= 1.0
        # if all atoms of the molecule have x, y, or z < 0.0, shift up
        elseif sum(xf[xyz, :] .> 0.0) == 0
            outside_box = true
            xf[xyz, :] += 1.0
        end
    end

    # update its Cartesian coordinates if it was outside of the box.
    if outside_box
        molecule.x = box.f_to_c * xf
    end
end

"""
    translate_molecule!(molecule, simulation_box)

Perturbs the Cartesian coordinates of a molecule by a random vector of max length δ.
Applies periodic boundary conditions to keep the molecule inside the simulation box.
"""
function translate_molecule!(molecule::Molecule, simulation_box::Box)
    # store old coordinates and return at the end for possible restoration of old coords
    x_old = deepcopy(molecule.x)
    # peturb in Cartesian coords in a random cube centered at current coords.
    dx = δ * (rand(3, 1) - 0.5) # move every atom of the molecule by the same vector.
    # change coordinates of molecule
    molecule.x .+= dx
    # done, unless the molecule has moved completely outside of the box...
    bring_molecule_inside_box!(molecule, simulation_box)

    return x_old # in case we need to restore
end

"""
    proposal_energy = guest_guest_vdw_energy(molecule_id, molecules,
        ljforcefield, simulation_box)

Calculates energy of a single adsorbate in the system. This can be used to find
the change in energy after an accepted proposal

Code copied from Arni's Energetics.jl with minor adjustments to calculate
    interactions between the adsorbates as well as the framework
"""
function guest_guest_vdw_energy(molecule_id::Int, molecules::Array{Molecule, 1},
                                ljforcefield::LennardJonesForceField, simulation_box::Box)
    energy = 0.0 # energy is pair-wise additive
    # fractional coordinates of molecule
    xf_molecule = simulation_box.c_to_f * molecules[molecule_id].x
    # Look at interaction with all other molecules in the system
    for other_molecule_id = 1:length(molecules)
        # molecule cannot interact with itself
        if other_molecule_id == molecule_id
            continue
        end
        # fractional coodinates of other molecule
        xf_other_molecule = simulation_box.c_to_f * molecules[other_molecule_id].x
        # Loop over all atoms in the given molecule
        for atom_id = 1:molecules[molecule_id].n_atoms
            # loop over every atom in the other molecule
            for other_atom_id = 1:molecules[other_molecule_id].n_atoms
                # compute vector between molecules in fractional coordinates
                dxf = xf_molecule[:, atom_id] - xf_other_molecule[:, other_atom_id]

                nearest_image!(dxf, (1, 1, 1))
                # converts fractional distance to cartesian distance
                dx = simulation_box.f_to_c * dxf

                r² = dot(dx, dx)

                if r² < R_OVERLAP_squared
                    return Inf
                elseif r² < ljforcefield.cutoffradius_squared
                    # TODO test whether it is more efficient to store this as a variable up top
                    energy += lennard_jones(r²,
                        ljforcefield.σ²[molecules[other_molecule_id].atoms[other_atom_id]][molecules[molecule_id].atoms[atom_id]],
                        ljforcefield.ϵ[molecules[other_molecule_id].atoms[other_atom_id]][molecules[molecule_id].atoms[atom_id]])
                end
            end # loop over all atoms in other molecule
        end # loop over all atoms of molecule_id of interest
    end # loop over all other molecules
    return energy # units are the same as in ϵ for forcefield (Kelvin)
end

"""
Compute total guest-host interaction energy (sum over all adsorbates).
"""
function total_guest_host_vdw_energy(framework::Framework, molecules::Array{Molecule, 1}, ljforcefield::LennardJonesForceField, repfactors::Tuple{Int, Int, Int})
    total_energy = 0.0
    for molecule_id = 1:length(molecules)
        total_energy += vdw_energy(framework, molecules[molecule_id], ljforcefield, repfactors)
    end
    return total_energy
end

"""
Compute sum of all guest-guest interaction energy from vdW interactions.
"""
function total_guest_guest_vdw_energy(molecules::Array{Molecule, 1}, ljforcefield::LennardJonesForceField, simulation_box::Box)
    total_energy = 0.0
    for molecule_id = 1:length(molecules)
        total_energy += guest_guest_vdw_energy(molecule_id, molecules, ljforcefield, simulation_box)
    end
    return total_energy / 2.0 # avoid double-counting pairs
end

"""
    results = gcmc_simulation(framework, temperature, pressure, adsorbate, ljforcefield)

runs a monte carlo simulation using the given framework and adsorbate.
runs at the given temperature and pressure
"""
#will pass in molecules::Array{Molecule} later
function gcmc_simulation(framework::Framework, temperature::Float64, fugacity::Float64,
                         adsorbate::Symbol, ljforcefield::LennardJonesForceField; n_sample_cycles::Int=100000,
                         n_burn_cycles::Int=10000, sample_frequency::Int=25, verbose::Bool=false)
    if verbose
        pretty_print(adsorbate, framework.name, temperature, fugacity)
    end

    const repfactors = replication_factors(framework.box, ljforcefield)
    const simulation_box = replicate_box(framework.box, repfactors)

    current_energy_gg = 0.0 # only true if starting with 0 molecules
    current_energy_gh = 0.0
    gcmc_stats = GCMCstats(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    molecules = Molecule[]
    
    markov_counts = MarkovCounts([0 for i = 1:length(PROPOSAL_ENCODINGS)], [0 for i = 1:length(PROPOSAL_ENCODINGS)])
    
    # (n_burn_cycles + n_sample_cycles) is # of outer cycles; for each outer cycle, peform max(20, # molecules in the system) 
    #  Markov chain proposals.
    markov_chain_time = 0
    for outer_cycle = 1:(n_burn_cycles + n_sample_cycles), inner_cycle = 1:max(20, length(molecules))
        markov_chain_time += 1
        
        # choose move randomly; keep track of proposals
        which_move = rand(1:N_PROPOSAL_TYPES)
        markov_counts.n_proposed[which_move] += 1
        
        if which_move == INSERTION
            insert_molecule!(molecules, simulation_box, adsorbate)

            U_gg = guest_guest_vdw_energy(length(molecules), molecules, ljforcefield, simulation_box)
            U_gh = vdw_energy(framework, molecules[end], ljforcefield, repfactors)

            # Metropolis Hastings Acceptance for Insertion
            if rand() < fugacity * simulation_box.Ω / (length(molecules) * KB *
                    temperature) * exp(-(U_gh + U_gg) / temperature)
                # accept the move, adjust current_energy
                markov_counts.n_accepted[which_move] += 1

                current_energy_gg += U_gg
                current_energy_gh += U_gh
            else
                # reject the move, remove the inserted molecule
                pop!(molecules)
            end
        elseif (which_move == DELETION) && (length(molecules) != 0)
            # propose which molecule to delete
            molecule_id = rand(1:length(molecules))

            U_gg = guest_guest_vdw_energy(molecule_id, molecules, ljforcefield,
                simulation_box)
            U_gh = vdw_energy(framework, molecules[molecule_id], ljforcefield,
                repfactors)

            # Metropolis Hastings Acceptance for Deletion
            if rand() < length(molecules) * KB * temperature / (fugacity * 
                    simulation_box.Ω) * exp((U_gh + U_gg) / temperature)
                # accept the deletion, delete molecule, adjust current_energy
                markov_counts.n_accepted[which_move] += 1

                delete_molecule!(molecule_id, molecules)

                current_energy_gg -= U_gg
                current_energy_gh -= U_gh
            end
        elseif (which_move == TRANSLATION) && (length(molecules) != 0)
            # propose which molecule whose coordinates we should perturb
            molecule_id = rand(1:length(molecules))

            # energy of the molecule before it was translated
            U_gg_old = guest_guest_vdw_energy(molecule_id, molecules,
                ljforcefield, simulation_box)
            U_gh_old = vdw_energy(framework, molecules[molecule_id],
                ljforcefield, repfactors)

            x_old = translate_molecule!(molecules[molecule_id], simulation_box)

            # energy of the molecule after it is translated
            U_gg_new = guest_guest_vdw_energy(molecule_id, molecules,
                ljforcefield, simulation_box)
            U_gh_new = vdw_energy(framework, molecules[molecule_id],
                ljforcefield, repfactors)

            # Metropolis Hastings Acceptance for translation
            if rand() < exp(-((U_gg_new + U_gh_new) - (U_gg_old + U_gh_old))
                / temperature)
                # accept the move, adjust current energy
                markov_counts.n_accepted[which_move] += 1

                current_energy_gg += U_gg_new - U_gg_old
                current_energy_gh += U_gh_new - U_gh_old
            else
                # reject the move, reset the molecule at molecule_id
                molecules[molecule_id].x = deepcopy(x_old)
            end
        end # which move the code executes

        # TODO remove after testing.
        for m = 1:length(molecules)
            @assert(! completely_outside_box(molecules[m], simulation_box), "molecule outside box!")
        end

        # sample the current configuration
        if (outer_cycle > n_burn_cycles) && (markov_chain_time % sample_frequency == 0)
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
    end #finished markov chain proposal moves

    # compute total energy, compare to `current_energy*` variables where were incremented
    total_U_gh = total_guest_host_vdw_energy(framework, molecules, ljforcefield, repfactors)
    total_U_gg = total_guest_guest_vdw_energy(molecules, ljforcefield, simulation_box)
    if ! isapprox(total_U_gh, current_energy_gh, atol=0.01)
        println("U_gh, incremented = ", current_energy_gh)
        println("U_gh, computed at end of simulation =", total_U_gh)
        error("guest-host energy incremented improperly")
    end
    if ! isapprox(total_U_gg, current_energy_gg, atol=0.01)
        println("U_gg, incremented = ", current_energy_gg)
        println("U_gg, computed at end of simulation =", total_U_gg)
        error("guest-guest energy incremented improperly")
    end

    @assert(markov_chain_time == sum(markov_counts.n_proposed))

    results = Dict{String, Any}()
    results["crystal"] = framework.name
    results["adsorbate"] = adsorbate
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
