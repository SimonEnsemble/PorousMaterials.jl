module GCMC
using PorousMaterials
using Energetics_Util
export gcmc_sim

#Arbitrary value for delta for now, can change or make more complex later
const δ = 0.1

"""
    insert(molecules::Array, framework::framework, repfactor::Tuple)

Inserts an additional molecule into the system and then checks the metropolis
hastings acceptance rules to find whether the proposed insertion is accepted or
rejected. Function then returns the new list of molecules
"""
function insert(molecules::Array{Molecule}, framework::framework, sim_box::Box)
    U_old = gg_energy(molecules, simulation_box)
    new_mol = Molecule()
    prop_energy = get_energy(proposal, framework)
    #Metropolis Hastings Acceptance rules for inserting a random particle
    if(rand() < metropolis_hastings_acceptance)
        return proposal
    else
        return molecules
    end
end

"""
    delete(molecules,framework,repfactor)

Removes a random molecule from the current molecules in the framework. Checks
metropolis hastings to see whether this is approved or rejected, and returns the
new list of molecules
"""
#TODO add guest-guest interactions, using this as a baseline to make sure it works
function delete(molecules::Array{Molecule}, framework::Framework)
    U_old = gg_energy(molecules, framework.box) + gh_energy(molecules, framework)

    #Metropolis Hastings Acceptance rules for deleting a random particle
    if(rand() > metropolis_hastings_acceptance)

    end
end

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
            xf_molecules[coords, :] += 1.0
        end
    end
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
            end
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
                end
            end
        end
    end #TODO label what these ends end
    return energy #units are the same as in ϵ for forcefield (Kelvin)
end #function



"""

Need sim_box in order to do nearest image convention
"""
function gg_energy(molecules::Array{Molecule},sim_box::Box)
#    for every atom in the molecule
#        for every molecule in the sim
#            for every atom in the pairing molecule
#                calculate lennard jones energy, add to sum
#    return sum of all energy
end

"""

Copied Arni's vdw-energy function, and added a line to calculate the energy of
multiple molecules in the MOF interacting with the framework
"""
function gh_energy(framework::Framework, molecules::Array{Molecule},
                    ljforcefield::LennardJonesForceField, repfactors::Tuple{Int64, Int64, Int64})
    energy = 0.0
    #loop over replications of the home unit cell to build the supercell
    for nA = 0:repfactors[1]-1, nB = 0:repfactors[2]-1, nC = 0:repfactors[3]-1
        #loop over all molecules in the framework
        for l = 1:length(molecules)
            #loop over atoms in the current molecule/adosrbate
            for i = 1:molecules[l].n_atoms
                xf_molecule_atom = mod.(simulation_box.c_to_f * molecule.x[:,i],repfactors)
                #loop over framework atoms in the current cell

end

function gcmc_sim(framework::Framework, molecule::Molecule, T, P)

end

end
