module GCMC
using PorousMaterials

export gcmc_sim, sim_for_N_equilibrium, sim_for_energy_equilibrium

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
function delete(molecules::Array{Molecule}, framework::framework, sim_box::Box)
    proposal = [molecules[i] for i = 1:length(molecules)]
    deleteat!(proposal,rand(1:length(proposal)))
    base_energy = get_energy(molecules, framework)
    prop_energy = get_energy(proposal, framework)
    #Metropolis Hastings Acceptance rules for deleting a random particle
    if(rand() < metropolis_hastings_acceptance)
        return proposal
    else
        return molecules
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
"""
function one_mol_energy(k::Int, molecules::Array{Molecule}, framework::framework, sim_box::Box)
    #TODO Use vdw_energy to find energy initially
    #TODO Develop guest-guest interactions equation
    #TODO find a way to get ljforcefield quickly, get initial input from higher function?
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

"""
"
    N_equilibrium = sim_for_N_equilibrium()

runs a MC simulation to calculate the average value of N at the given:
Pressure
Temperature
The molecule being inserted
The framework being tested
"
function sim_for_N_equilibrium()

end

"
    energy_equilibrium = sim_for_energy_equilibrium()

runs a MC simulation to calculate the average value of energy at the given:
Pressure
Temperature
Based on:
The molecule being inserted
The framework being tested
"
function sim_for_energy_equilibrium()

end
"""
