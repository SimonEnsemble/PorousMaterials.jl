module GCMC
using PorousMaterials

export gcmc_sim

"""
    insert(molecule_pos::Array, framework::framework, repfactor::Tuple)

Inserts an additional molecule into the system and then checks the metropolis
hastings acceptance rules to find whether the proposed insertion is accepted or rejected.
"""
function insert(molecule_pos::Array{Tuple{Float64,Float64,Float64}}, framework::framework, repfactor::Tuple{Int64,Int64,Int64})
    proposal = [molecule_pos[i] for i = 1:length(molecule_pos)]
    push!(proposal, (rand()*repfactor[1],rand()*repfactor[2],rand()*repfactor[3]))
    base_energy = get_energy(molecule_pos, framework)
    prop_energy = get_energy(proposal, framework)
    #Metropolis Hastings Acceptance rules for inserting a random particle
end

#TODO add guest-guest interactions, using this as a baseline to make sure it works
function delete(molecule_pos::Array{Tuple{Float64,Float64,Float64}}, framework::framework, repfactor::Tuple{Int64,Int64,Int64})
    proposal = [molecule_pos[i] for i = 1:length(molecule_pos)]
    deleteat!(proposal,rand(1:length(proposal)))
    base_energy = get_energy(molecule_pos, framework)
    prop_energy = get_energy(proposal, framework)
    #Metropolis Hastings Acceptance rules for deleting a random particle
end

function move(molecule_pos::Array{Tuple{Float64,Float64,Float64}}, framework::framework, repfactor::Tuple{Int64,Int64,Int64})
    proposal = [molecule_pos[i] for i = 1:length(molecule_pos)]
    #Arbitray value for delta for now, can change or make more complex later
    δ = 0.1
    move_spot = rand(1:length(proposal))
    proposal[move_spot] = (δ*(proposal[move_spot])[1], δ*(proposal[move_spot])[2], δ*(proposal[move_spot])[3])
    base_energy = get_energy(molecule_pos, framework)
    prop_energy = get_energy(proposal, framework)
    #Metropolis Hastings Acceptance rules for moving a random particle
end

function get_N(molecule_pos)
    return length(molecule_pos)
end

function get_energy(molecule_pos::Array(Tuple), framework::framework)
    #TODO Use vdw_energy to find energy initially
    #TODO Develop guest-guest interactions equation
    #TODO find a way to get ljforcefield quickly, get initial input from higher function?
end

function sim_for_N_equilibrium()

end

function sim_for_energy_equilibrium()

end

function gcmc_sim(framework::framework, molecule::molecule, T, P)

end

end
