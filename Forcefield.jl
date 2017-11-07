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
	close(f)
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
	ratio = (σ/r)^6
	return 4*ϵ*(ratio^2 - ratio)
end # function end

"""
	repfactors = rep_factors(frame,cutoff)

Find the repetition factors needed to make a supercell big enough to fit a sphere with the specified cutoff radius.
Rather than adding all the new atoms to the coordinate matrix, we only keep check on how many times they're repeated.
"""

function rep_factors(frame::Framework, cutoff::Float64)
	# Unit vectors used to transform from fractional coordinates to cartesian coordinates. We'll be 
	a = frame.f_to_C[:,1]
	b = frame.f_to_C[:,2]
	c = frame.f_to_C[:,3]
	
	n_ab = cross(a,b)
	n_ac = cross(a,c)
	n_bc = cross(b,c)
	
	# c0 defines a center in the unit cell
	c0 = [a b c] * [.5, .5, .5]

	rep = [1, 1, 1]

	# Repeat for `a`
	while abs(dot(n_bc, c0)) / vecnorm(n_bc) < cutoff
		rep[1] += 1
		a += frame.f_to_C[:,1]
		c0 = [a b c] * [.5, .5, .5]
	end
	
	# Repeat for `b`
	while abs(dot(n_ac, c0)) / vecnorm(n_ac) < cutoff
		rep[2] += 1
		b += frame.f_to_C[:,2]
		c0 = [a b c] * [.5, .5, .5]
	end
	
	# Repeat for `c`
	while abs(dot(n_ab, c0)) / vecnorm(n_ab) < cutoff
		rep[3] += 1
		c += frame.f_to_C[:,3]
		c0 = [a b c] * [.5, .5, .5]
	end

	return rep
end


function vdw_energy(frame::Framework, molecule::Molecule, ljforcefield::LennardJonesForceField, pos)
	r = 0.0
	σ = 0.0
	ϵ = 0.0
	sum = 0
	for i=1:molecule.n_atoms
		for k=1:frame.n_atoms
			r = norm(molecule.x[i,:]-(frame.f_coords*frame.f_to_C)[k,:])
			σ = ljforcefield.sigmas[ljforcefield.atom_to_id[frame.atoms[k]],ljforcefield.atom_to_id[molecule.atoms[i]]]
			ϵ = ljforcefield.epsilons[ljforcefield.atom_to_id[frame.atoms[k]], ljforcefield.atom_to_id[molecule.atoms[i]]]
			if (r < ljforcefield.cutoffradius && r > 0.1)
				sum += lennard_jones(r,σ,ϵ)	
			end
		end
	end	
	return sum
end # function end

end # end module
