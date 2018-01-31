module GCMC

#TODO how to include other porous materials files
include ("Crystal.jl")
using Energetics_Util.jl
export gcmc_sim

#Arbitrary value for delta for now, can change or make more complex later
const δ = 0.1

"""
    insert_molecule!(molecules, simulation_box, adsorbate)

Inserts an additional molecule into the system and then checks the metropolis
hastings acceptance rules to find whether the proposed insertion is accepted or
rejected. Function then returns the new list of molecules
"""
function insert_molecule!(molecules::Array{Molecule}, simulation_box::Box,
        adsorbate::String)
    #TODO how do I randomly create a Molecule? It should be a given type, but
    #     then how do I generate locations for the atoms in the molecule
    #TODO create template to handle more complex molecules
    x_new = simulation_box.f_to_c * [rand(), rand(), rand()]
    new_molecule = Molecule(1, [adsorbate], x_new, [0.0])
    push!(molecules, new_molecule)
end

"""
    delete_molecule!(molecule_id, molecules)

Removes a random molecule from the current molecules in the framework.
molecule_id decides which molecule will be deleted, for a simulation, it must
    be a randomly generated value
"""
function delete_molecule!(molecule_id::Int, molecules::Array{Molecule})
    # could also generate a value here, would it work to return two values?
    #molecule_id = rand(0:length(molecules))
    deleted_molecule = molecules[molecule_id]
    deleteat!(molecules, molecule_id)
    return deleted_molecule
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

Code copied from Arni's Energetics.jl with minor adjustments to calculate interactions
between the adsorbates as well as the framework
"""
function guest_guest_vdw_energy(molecule_id::Int, molecules::Array{Molecule},
                        ljforcefield::LennardJonesForceField, simulation_box::Box)
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
    gcmc_sim(framework, molecules, temperature, pressure)
"""
#will pass in molecules::Array{Molecule} later
function gcmc_sim(framework::Framework,
                        temperature::Float64, pressure::Float64,
                        adsorbate::String, ljforcefield::LennardJonesForceField;
                        number_mc_trials::Int64=100000)
    const repfactors = replication_factors(framework.box, ljforcefield)
    const simulation_box = replicate_box(framework.box, repfactors)

    #used to calculate average energy of the system
    total_energy = 0.0
    #used to add accurate values to total_energy
    current_energy = 0.0
    #used to calculate average number of molecules
    total_n = 0.0
    #start with 0 molecules
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

    for t = 1:number_mc_trials
        proposal_id = rand(1:3)
        if proposal_id == 1
            insert_molecule!(molecules, framework.box)
            U_gg = guest_guest_vdw_energy(length(molecules), molecules,
                ljforcefield, simulation_box)
            U_gh = vdw_energy(framework, molecules[end],
                ljforcefield, repfactors)
            if rand() < exp(-(U_gg + U_gh) / temperature)
                #accept the move, add the energy of the molecule to the total energy of the system
                current_energy += U_gg + U_gh
            else
                #reject the move, remove the added molecule
                pop!(molecules)
            end
        elseif proposal_id == 2
            molecule_id = rand(1:length(molecules))
            U_gg = guest_guest_vdw_energy(molecule_id, molecules, ljforcefield,
                simulation_box)
            U_gh = vdw_energy(framework, molecules[molecule_id], ljforcefield,
                repfactors)

        else
            molecule_id = rand(1:length(molecules))
            old_molecule = translate_molecule(molecule_id, molecules,
                simulation_box)
            
        end
        total_energy += current_energy
    end

end #gcmc_sim

end
