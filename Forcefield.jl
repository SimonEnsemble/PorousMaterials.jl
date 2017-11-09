"""All things energy related"""
module Forcefield

using Crystal
using Mols

export LennardJonesForceField, ljffConstruct, lennard_jones, readElementProps, rep_factors, centerOfMass, vdw_energy

const NA = 6.022e23 # 1/mol
const R = 1.9872036e-3 # kcal/(mol K)

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
			epsilon_dict[str[1]] = ϵ/(NA*R)
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
	ratio = (σ/r)^6
	return 4*ϵ*(ratio^2 - ratio)
end # function end

"""
	repfactors = rep_factors(frame,cutoff)

Find the replication factors needed to make a supercell big enough to fit a sphere with the specified cutoff radius.
Rather than adding all the new atoms to the coordinate matrix, we only keep check on how many times they're replicated.
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
	# |n_bc ⋅ c0|/|n_bc| defines the distance from the end of the supercell and the center. As long as that distance is less than the cutoff radius, we need to increase it
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

"""
	readElementProps("filepath")

Reads a .csv file with the following header: Element name, Molecular mass, Atomic radius, Ionic radius and returns a dictionary with an element name pointing at an array of values,
EleProps[Element] => [Molecular mass(amu), Atomic radius(Angstrom), Ionic radius(Angstrom)]
"""
function readElementProps(filename::String)
	EleProps = Dict{String,Array{Float64}}()
	f = open(filename,"r")
	lines = readlines(f)

	for (i,line) in enumerate(lines)
		if (i>1)
			str = split(line,",")
			temp = zeros(3)
			for (k,val) in enumerate(str[2:4])
				if val!="NA"
					temp[k] = parse(Float64,val)
				else
					temp[k] = 0.0
				end
			end
			EleProps[str[1]] = temp
		end
	end
	close(f)
	return EleProps
end


"""
	centerOfMass(frame,"~/example/properties.csv")

Uses `readElementProps` to get a dictionary of element properties and uses that to calculate the center of mass of a supercell made from Framework and replication factors. (See readElementProps for more info on that function)
Calculates the center of mass according to r_cm = ∑r_i*m_i/m_tot in fractional coordinates but returns in cartesian coordinates
"""
function centerOfMass(frame::Framework, filename::String, rep_factors::Array{Int64})
	repA = rep_factors[1]
	repB = rep_factors[2]
	repC = rep_factors[3]

	EleProps = readElementProps(filename)
	rvec = zeros(3)
	mtot = 0.0
	for nA=1:repA, nB=1:repB, nC=1:repC
		for i=1:frame.n_atoms
			rvec += (frame.f_coords[i,:]+[nA-1, nB-1, nC-1])*EleProps[frame.atoms[i]][1]
			mtot += EleProps[frame.atoms[i]][1]
		end
	end
	return frame.f_to_C*(rvec/mtot)
end

"""
	vdw_energy(frame, molecule, ljforcefield, pos, repfactors)

Calculates the van der Waals energy for a molecule locates at a specific position in a MOF framework.
"""
function vdw_energy(frame::Framework, molecule::Molecule, ljforcefield::LennardJonesForceField, pos::Array{Float64}, repfactors::Array{Int64})
	repA = repfactors[1]
	repB = repfactors[2]
	repC = repfactors[3]
	repvec = Array{Int64}(3)
	σ = 0.0
	r = 0.0
	ϵ = 0.0
	potsum = 0
	for nA=0:repA-1, nB=0:repB-1, nC=0:repC-1
		for i=1:molecule.n_atoms
			for k=1:frame.n_atoms
				# Nearest image convention. If the interaction between the probe molecule and atom k is being looked at, we'll only look at the interaction between the probe molecule and the closest replication of atom k. This is done with fractional coordinates for simplication and transformation to cartesian is done later.
				repvec = [nA, nB, nC]
				if frame.f_coords[k,1]+nA > repA/2
					repvec -= [1,0,0]
				end
				if frame.f_coords[k,2]+nB > repB/2
					repvec -= [0,1,0]
				end
				if frame.f_coords[k,3]+nC > repC/2
					repvec -= [0,0,1]
				end

				r = norm((pos+molecule.x[i,:])-frame.f_to_C*(frame.f_coords[k,:]+repvec))
				σ = ljforcefield.sigmas[ljforcefield.atom_to_id[frame.atoms[k]],ljforcefield.atom_to_id[molecule.atoms[i]]]
				ϵ = ljforcefield.epsilons[ljforcefield.atom_to_id[frame.atoms[k]], ljforcefield.atom_to_id[molecule.atoms[i]]]
				if (r < ljforcefield.cutoffradius)
					potsum += lennard_jones(r,σ,ϵ)
				end
			end
		end
	end
	return potsum
end # function end

end # end module
