""" """
module Forcefield

export LennardJonesForceField, ljffConstruct, lennard_jones_potential_energy

struct LennardJonesForceField
	cutoffradius::Array{Float64}
	epsilon_dict::Dict{String, Float64}
	sigma_dict::Dict{String, Float64}
	atom_to_id::Dict{String, Int64}
	epsilons::Array{Float64, 2}
	sigmas::Array{Float64,2}
end

function ljffConstruct(filename::String)
	f = open(filename,"r")
	lines = readlines(f)

	n_ele = length(lines)-1
	
	cutoffradius = Array{Float64}(n_ele)
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
			cutoffradius[i-1] = 2.5*σ
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

	@printf("%s\t%s\t%s\t%s\t%s\t%s",typeof(cutoffradius),typeof(epsilon_dict),typeof(sigma_dict),typeof(atom_to_id),typeof(epsilons),typeof(sigmas))
	return LennardJonesForceField(cutoffradius, epsilon_dict, sigma_dict, atom_to_id, epsilons, sigmas)

end # function end

function lennard_jones_potential_energy(r::Float64, ljforcefield::LennardJonesForceField, ele1::String, ele2::String)
	σ = ljforcefield.sigmas[ljforcefield.atom_to_id[ele1], ljforcefield.atom_to_id[ele2]]
	if (r > 2.5*σ)
		NA = 6.022e23
		kcal_to_kJ = 4.184
		ratio = (σ/r)^2
		epsilon = ljforcefield.epsilons[ljforcefield.atom_to_id[ele1], ljforcefield.atom_to_id[ele2]]*kcal_to_kJ*NA
		return 4*epsilon*(ratio^2 - ratio)
	else
		return 0
	end
end # function end



end # end module
