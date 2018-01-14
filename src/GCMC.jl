module GCMC
using PorousMaterials

export gcmc_sim

"""
    insert(molecule_pos::Array, framework::framework, repfactor::Tuple)

Inserts an additional molecule into the system and then checks the metropolis
hastings acceptance rules to find whether the proposed insertion is accepted or
rejected. Function then returns the new list of molecules
"""
function insert(molecule_pos::Array{Tuple{Float64,Float64,Float64}}, framework::framework, repfactor::Tuple{Int64,Int64,Int64})
    proposal = [molecule_pos[i] for i = 1:length(molecule_pos)]
    push!(proposal, (rand()*repfactor[1],rand()*repfactor[2],rand()*repfactor[3]))
    base_energy = get_energy(molecule_pos, framework)
    prop_energy = get_energy(proposal, framework)
    #Metropolis Hastings Acceptance rules for inserting a random particle
    if(rand() < metropolis_hastings_acceptance)
        return proposal
    else
        return molecule_pos
    end
end

"""
    delete(molecule_pos,framework,repfactor)

Removes a random molecule from the current molecules in the framework. Checks
metropolis hastings to see whether this is approved or rejected, and returns the
new list of molecules
"""
#TODO add guest-guest interactions, using this as a baseline to make sure it works
function delete(molecule_pos::Array{Tuple{Float64,Float64,Float64}}, framework::framework, repfactor::Tuple{Int64,Int64,Int64})
    proposal = [molecule_pos[i] for i = 1:length(molecule_pos)]
    deleteat!(proposal,rand(1:length(proposal)))
    base_energy = get_energy(molecule_pos, framework)
    prop_energy = get_energy(proposal, framework)
    #Metropolis Hastings Acceptance rules for deleting a random particle
    if(rand() < metropolis_hastings_acceptance)
        return proposal
    else
        return molecule_pos
    end
end

"""
    move(molecule_pos, framework, repfactor)

Translates a random molecule to a nearby coordinate. Then uses metropolis hastings
to approve or reject the move. This then returns the new list of molecules
"""
function move(molecule_pos::Array{Tuple{Float64,Float64,Float64}}, framework::framework, repfactor::Tuple{Int64,Int64,Int64})
    proposal = [molecule_pos[i] for i = 1:length(molecule_pos)]
    #Arbitray value for delta for now, can change or make more complex later
    δ = 0.1
    move_spot = rand(1:length(proposal))
    proposal[move_spot] = (δ*(proposal[move_spot])[1], δ*(proposal[move_spot])[2], δ*(proposal[move_spot])[3])
    base_energy = get_energy(molecule_pos, framework)
    prop_energy = get_energy(proposal, framework)
    #Metropolis Hastings Acceptance rules for moving a random particle
    if(rand() < metropolis_hastings_acceptance)
        return proposal
    else
        return molecule_pos
    end
end

"""
    num_molecules = get_N(molecule_pos)

returns the number of molecules currently in the framework/simulation
"""
function get_N(molecule_pos)
    return length(molecule_pos)
end

"""
    current_energy = get_energy(molecule_pos, framework)

returns the current energy in the system. Calculates energy based on each molecule
interacting with the framework and every other molecule.
"""
function get_energy(molecule_pos::Array(Tuple), framework::framework)
    #TODO Use vdw_energy to find energy initially
    #TODO Develop guest-guest interactions equation
    #TODO find a way to get ljforcefield quickly, get initial input from higher function?
end

"""
    N_equilibrium = sim_for_N_equilibrium()

runs a MC simulation to calculate the average value of N at the given:
Pressure
Temperature
The molecule being inserted
The framework being tested
"""
function sim_for_N_equilibrium()

end

"""
    energy_equilibrium = sim_for_energy_equilibrium()

runs a MC simulation to calculate the average value of energy at the given:
Pressure
Temperature
Based on:
The molecule being inserted
The framework being tested
"""
function sim_for_energy_equilibrium()

end

function gcmc_sim(framework::framework, molecule::molecule, T, P)

end

end
