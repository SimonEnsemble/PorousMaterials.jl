module GCMC
using PorousMaterials

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
    U_old = gg_energy(molecules, framework.box)
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
    move(molecules, framework, repfactor)

Translates a random molecule to a nearby coordinate. Then uses metropolis hastings
to approve or reject the move. This then returns the new list of molecules
"""
function move(molecules::Array{Molecule}, framework::framework)
    mol_k = rand(1:length(molecules))
    U_old = one_mol_energy(k, molecules, framework)
    old_mol = molecules[k]
    new_mol = old_mol
    new_mol.x[1,:] += δ * (rand() - 0.5)
    new_mol.x[2,:] += δ * (rand() - 0.5)
    new_mol.x[3,:] += δ * (rand() - 0.5)
    molecules[k] = new_mol
    U_new = one_mol_energy(k, molecules, framework)
    #Metropolis Hastings Acceptance rules for moving a random particle
    if(rand() > metropolis_hastings_acceptance)
        molecules[k] = old_mol
    end
end

"""
    num_molecules = get_N(molecules)

returns the number of molecules currently in the framework/simulation
"""
function get_N(molecules::Array{Molecule})
    return length(molecules)
end

"""
    current_energy = one_mol_energy(k, molecules, framework)

returns the current energy in the system. Calculates energy based on each molecule
interacting with the framework and every other molecule.

Code copied from Arni's Energetics.jl with minor adjustments to calculate interactions
between the adsorbates as well as the framework
"""
function one_mol_energy(m::Int, molecules::Array{Molecule}, framework::framework
                        ljforcefield::LennardJonesForceField, repfactors::Tuple{Int64,Int64,Int64})
    #TODO Use vdw_energy to find energy initially
    #TODO Develop guest-guest interactions equation
    #TODO find a way to get ljforcefield quickly, get initial input from higher function?
    energy = 0.0
    # loop over replications of the home unit cell to build the supercell (simulation box)
    for nA = 0:repfactors[1]-1, nB = 0:repfactors[2]-1, nC = 0:repfactors[3]-1
        #Loop over all atoms in the given molecule
        for i = 1:molecules[m].n_atoms
            xf_molecule_atom = mod.(framework.box.c_to_f * molecules[m].x[:,i],repfactors)
            # loop over framework atoms in the home unit cell
			for k = 1:framework.n_atoms
                dxf = xf_molecule_atom - (framework.xf[:,k] + [nA, nB, nC])

                #runs nearest image convention and adjusts accordingly
                nearest_image(dxf,repfactors)
                #converts fractional distance to cartesian distance
                dx = framework.box.f_to_c * dxf

                r_squared = sum(dx .* dx)
                if r_squared < R_OVERLAP_squared
                    #the atoms are overlapping, return inf energy
                    return Inf
                elseif r_squared < ljforcefield.cutoffradius_squared
                    #atoms can interact and won't produce Inf energy
                    #I removed the sigma_squared and ϵ variables because they are only used once and would
                    #have to be initialized for every energy pair
                    energy += lennard_jones(r_squared,
                                            ljforcefield.sigmas_squared[framework.atoms[k]][molecule.atoms[i]],
                                            ϵ = ljforcefield.epsilons[framework.atoms[k]][molecule.atoms[i]])
                end
            end

            #loop over other adsorbate atoms leading up to central atom in array
            for l = 1:m-1
                #loop over every atom in second molecules
                for k = 1:molecules[l].n_atoms
                    dxf = xf_molecule_atom - molecules[l].x[:,k]

                    #runs nearest image convention and adjusts accordingly
                    nearest_image(dxf,repfactors)
                    #converts fractional distance to cartesian distance
                    dx = framework.box.f_to_c * dxf

                    r_squared = sum(dx .* dx)
                    if r_squared < R_OVERLAP_squared
                        #the atoms are overlapping, return inf energy
                        return Inf
                    elseif r_squared < ljforcefield.cutoffradius_squared
                        #atoms can interact and won't produce Inf energy
                        #I removed the sigma_squared and ϵ variables because they are only used once and would
                        #have to be initialized for every energy pair
                        energy += lennard_jones(r_squared,
                                                ljforcefield.sigmas_squared[framework.atoms[k]][molecule.atoms[i]],
                                                ϵ = ljforcefield.epsilons[framework.atoms[k]][molecule.atoms[i]])
                    end
                end
            end

            #Loop over other adsorbate atoms following central atom in array
            for l = (m+1):length(molecules)
                #loop over every atom in the second molecule
                for k = 1:molecules[l].n_atoms
                    dxf = xf_molecule_atom - molecule[l].x[:,k]

                    #runs nearest image convention and adjusts accordingly
                    nearest_image(dxf,repfactors)
                    #converts fractional distance to cartesian distance
                    dx = framework.box.f_to_c * dxf

                    r_squared = sum(dx .* dx)
                    if r_squared < R_OVERLAP_squared
                        #the atoms are overlapping, return inf energy
                        return Inf
                    elseif r_squared < ljforcefield.cutoffradius_squared
                        #atoms can interact and won't produce Inf energy
                        #I removed the sigma_squared and ϵ variables because they are only used once and would
                        #have to be initialized for every energy pair
                        energy += lennard_jones(r_squared,
                                                ljforcefield.sigmas_squared[framework.atoms[k]][molecule.atoms[i]],
                                                ϵ = ljforcefield.epsilons[framework.atoms[k]][molecule.atoms[i]])
                    end
                end
            end
        end
    end
end

"""
    nearest_image(distance,repfactors)
"""
function nearest_image(distance,repfactors)
    for j = 1:3
        if abs(distance[j]) > repfactors[j] / 2
            distance[j] -= sign(distance[j]) * repfactors[j]
        end
    end
end

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
function gh_energy(framework::Framework, )

end

function gcmc_sim(framework::Framework, molecule::Molecule, T, P)

end

end
