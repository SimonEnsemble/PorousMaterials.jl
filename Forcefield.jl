"""All things energy related"""
module Forcefield

export LennardJonesForceField, ljffConstruct, lennard_jones

const NA = 6.022e23
const kcal_to_kJ = 4.184

"""
	ljforcefield = LennardJonesForceField(cutoffradius, epsilon_dict, sigma_dict, atom_to_id, epsilons, sigmas)

Data structure for a Lennard Jones forcefield, read from a file containing UFF parameters.

# Arguments
- `cutoffradius::Float64`: cut-off radius to define where the potential energy goes to zero, rather than taking the limit as r → ∞ (units: Angstrom)
- `epsilon_dict::Dict{String, Float64}`: Dictionary that connects element acronyms to an ϵ, which is the depth of a Lennard Jones potential well
- `sigma_dict::Dict{String, Float64}`: Dictionary that connects element acronyms to a σ, which is the finite distance where the potential between atoms goes to zero
- `atom_to_id::Dict{String, Int64}`: Dictionary that connects element acronyms to a unique integer value.
- `epsilons::Array{Float64,2}`: Two dimensional matrix that contains the ϵ interacting values between two elements. Row/Column number correspond to id's from atom_to_id (units: kcal/mol)
- `sigmas::Array{Float64,2}`: Two dimensional matrix that contains the σ interacting values between two elements. Row/Column number correspond to id's from atom_to_id (units: Angstrom)
"""
struct LennardJonesForceField
	cutoffradius::Float64
	epsilon_dict::Dict{String, Float64}
	sigma_dict::Dict{String, Float64}
	atom_to_id::Dict{String, Int64}
	epsilons::Array{Float64, 2}
	sigmas::Array{Float64,2}
end

"""
	ljforcefield = ljffConstruct("filename.csv")

Read a .csv file containing UFF parameters (with the following column order: [Element, σ, ϵ]) and constructs a LennardJonesForceField object
"""
function ljffConstruct(filename::String)
	f = open(filename,"r")
	lines = readlines(f)

	n_ele = length(lines)-1
	
	cutoffradius = 14.0
	epsilon_dict = Dict{String, Float64}()
	sigma_dict = Dict{String, Float64}()
	atom_to_id = Dict{String, Int64}()
	elements = Array{String}(n_ele)
	epsilons = Array{Float64,2}(n_ele,n_ele)
	sigmas = similar(epsilons)

	for (i,line) in enumerate(lines)
		if (i > 1)
			str = split(line,",")
			σ,ϵ = map(x->parse(Float64, x), str[2:3])
			epsilon_dict[str[1]] = ϵ
			sigma_dict[str[1]] = σ	
			atom_to_id[str[1]] = i-1
			elements[i-1] = str[1]
		end
	end

	for (i,ele1) in enumerate(elements)
		for (k,ele2) in enumerate(elements[i:end])
			epsilons[i,k+i-1] = sqrt(epsilon_dict[ele1]*epsilon_dict[ele2])
			epsilons[k+i-1,i] = epsilons[i,k+i-1]
			sigmas[i,k+i-1] = (sigma_dict[ele1]+sigma_dict[ele2])/2
			sigmas[k+i-1,i] = sigmas[i,k+i-1]
		end
	end

#	@printf("%s\t%s\t%s\t%s\t%s\t%s",typeof(cutoffradius),typeof(epsilon_dict),typeof(sigma_dict),typeof(atom_to_id),typeof(epsilons),typeof(sigmas))
	return LennardJonesForceField(cutoffradius, epsilon_dict, sigma_dict, atom_to_id, epsilons, sigmas)

end # function end

"""
	V = lennard_jones_potential_energy(r, σ, ϵ)

Calculate the lennard jones potential energy given a radius r between two molecules. σ and ϵ are specific to interaction between two elements
"""
function lennard_jones(r::Float64, σ::Float64, ϵ::Float64 )
	ϵ = ϵ*kcal_to_kJ*NA
	ratio = (σ/r)^2
	return 4*ϵ*(ratio^2 - ratio)
end # function end

function vdw_energy(ljforcefield::LennardJonesForceField, molecule::Molecule, pos)
	@pass
end # function end


end # end module
