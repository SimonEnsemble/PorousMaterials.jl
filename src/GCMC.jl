module GCMC
using PorousMaterials
using Energetics_Util
export gcmc_sim

#Arbitrary value for delta for now, can change or make more complex later
const δ = 0.1

"""
    insert_random_molecule!(molecules, simulation_box)

Inserts an additional molecule into the system and then checks the metropolis
hastings acceptance rules to find whether the proposed insertion is accepted or
rejected. Function then returns the new list of molecules
"""
function insert_random_molecule!(molecules::Array{Molecule}, simulation_box::Box)
    #TODO how do I randomly create a Molecule? It should be a given type, but
    #     then how do I generate locations for the atoms in the molecule
    new_molecule = Molecule(n_atoms, atoms, new_molecule_position, charges)
    push!(molecules, new_molecule)
end #insert_random_molecule!

"""
    delete_random_molecule!(molecules, simulation_box)

Removes a random molecule from the current molecules in the framework.
"""
function delete_random_molecule!(molecules::Array{Molecule}, simulation_box::Box)
    molecule_id = rand(1:length(molecules))
    deleted_molecule = molecules[molecule_id]
    deleteat!(molecules, molecule_id)
    return deleted_molecule
end #delete_random_molecule!

"""
    translate_random_molecule(molecules)

Translates a random molecule a random amount in the given array and simulation
box
reflects if entire molecule goes outside in any cartesian direction
"""
function translate_random_molecule!(molecules::Array{Molecule}, simulation_box::Box)
    molecule_id = rand(1:length(molecules))
    dx = [δ * (rand() - 0.5) for i = 1:3]
    molecules[molecule_id].x .+= dx
    xf_molecule = simulation_box.c_to_f * molecules[molecule_id]
    for coords = 1:3
        if sum(xf_molecule[coords, :] .< 1.0) == 0
            xf_molecule[coords, :] -= 1.0
        elseif sum(xf_molecules[coords, :] .>) == 0
            xf_molecule[coords, :] += 1.0
        end #if statement that checks for reflection
    end #for loop to go over x, y, and, z coordinate systems
end #translate_random_molecule!

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
end #function

"""
    gcmc_sim(framework, molecules, temperature, pressure)
"""
function gcmc_sim(framework::Framework, molecules::Array{Molecule},
                        temperature::Float64, pressure::Float64)
    const NUMBER_SIMULATIONS = 1000000 #one million

end #gcmc_sim

end
