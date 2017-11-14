"""All things energy related"""
module Forcefield

using Crystal
using Mols

export LennardJonesForceField, constructforcefield, lennard_jones, readproperties, rep_factors, centerofmass, vdw_energy, rotate, exploreframe

const NA = 6.022e23 # 1/mol
const R = 1.9872036e-3 # kcal/(mol K)

"""
	ljforcefield = LennardJonesForceField(cutoffradius, epsilon_dict, sigma_dict, atom_to_id, epsilons, sigmas)

Data structure for a Lennard Jones forcefield, read from a file containing UFF parameters.

# Arguments
- `cutoffradius::Float64`: cut-off radius to define where the potential energy goes to zero, rather than taking the limit as r → ∞ (units: Angstrom)
- `epsilon_dict::Dict{AbstractString, Float64}`: Dictionary that connects element acronyms to an ϵ, which is the depth of a Lennard Jones potential well
- `sigma_dict::Dict{AbstractString, Float64}`: Dictionary that connects element acronyms to a σ, which is the finite distance where the potential between atoms goes to zero
- `atom_to_id::Dict{AbstractString, Int64}`: Dictionary that connects element acronyms to a unique integer value.
- `epsilons::Array{Float64,2}`: Two dimensional matrix that contains the ϵ interacting values between two elements. Row/Column number correspond to id's from atom_to_id (units: kcal/mol)
- `sigmas::Array{Float64,2}`: Two dimensional matrix that contains the σ interacting values between two elements. Row/Column number correspond to id's from atom_to_id (units: Angstrom)
"""
struct LennardJonesForceField
	cutoffradius::Float64
	epsilon_dict::Dict{AbstractString, Float64}
	sigma_dict::Dict{AbstractString, Float64}
	atom_to_id::Dict{AbstractString, Int64}
	epsilons::Array{Float64, 2}
	sigmas::Array{Float64,2}
end

"""
	ljforcefield = constructforcefield("filename.csv")

Read a .csv file containing UFF parameters (with the following column order: [Element, σ, ϵ]) and constructs a LennardJonesForceField object
"""
function constructforcefield(filename::AbstractString; cutoffradius::Float64=14.0)
	f = open(filename,"r")
	lines = readlines(f)

	n_ele = length(lines)-1

	epsilon_dict = Dict{AbstractString, Float64}()
	sigma_dict = Dict{AbstractString, Float64}()
	atom_to_id = Dict{AbstractString, Int64}()
	elements = Array{AbstractString}(n_ele)
	epsilons = Array{Float64,2}(n_ele,n_ele)
	sigmas = similar(epsilons)

	for (i,line) in enumerate(lines)
		if i > 1
			str = split(line,",")
			σ,ϵ = map(x->parse(Float64, x), str[2:3])
			epsilon_dict[str[1]] = ϵ/R
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

	return LennardJonesForceField(cutoffradius, epsilon_dict, sigma_dict, atom_to_id, epsilons, sigmas)

end # constructforcefield end

"""
	V = lennard_jones_potential_energy(r::Float64, σ::Float64, ϵ::Float64)

Calculate the lennard jones potential energy given a radius r between two molecules. σ and ϵ are specific to interaction between two elements
"""
function lennard_jones(r::Float64, σ::Float64, ϵ::Float64 )
	ratio = (σ/r)^6
#	if σ > r
#		@printf("σ = %f, r = %f\n",σ,r)
#	end
	return 4*ϵ*(ratio^2 - ratio)
end # lennard_jones end

"""
	repfactors = rep_factors(frame::Framework,cutoff::Float64)

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
end # end rep_factors

"""
	EleProps = readproperties("filepath")

Reads a .csv file with the following header: Element name, Molecular mass, Atomic radius, Ionic radius and returns a dictionary with an element name pointing at an array of values,
EleProps[Element] => [Molecular mass(amu), Atomic radius(Angstrom), Ionic radius(Angstrom)]
"""
function readproperties(filename::AbstractString)
	EleProps = Dict{AbstractString,Array{Float64}}()
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
end # end readproperties


"""
	c0 = centerofmass(frame,"~/example/properties.csv")

Uses `readproperties` to get a dictionary of element properties and uses that to calculate the center of mass of a supercell made from Framework and replication factors. (See readElementProps for more info on that function)
Calculates the center of mass according to r_cm = ∑r_i*m_i/m_tot in fractional coordinates but returns in cartesian coordinates
"""
function centerofmass(frame::Framework, filename::AbstractString, rep_factors::Array{Int64})
	repA = rep_factors[1]
	repB = rep_factors[2]
	repC = rep_factors[3]

	EleProps = readproperties(filename)
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
	V = vdw_energy(frame, molecule, ljforcefield, pos, repfactors)

Calculates the van der Waals energy for a molecule locates at a specific position in a MOF supercell. Uses the nearest image convention to find the closest replicate of a specific atom
"""
function vdw_energy(frame::Framework, molecule::Molecule, ljforcefield::LennardJonesForceField, pos::Array{Float64}, repfactors::Array{Int64})
	repA = repfactors[1]
	repB = repfactors[2]
	repC = repfactors[3]
	repvec = Array{Int64}(3)
	σ = 0.0
	r = 0.0
	ϵ = 0.0
	fpos = frame.C_to_f*pos
	potsum = 0
	for nA=0:repA-1, nB=0:repB-1, nC=0:repC-1
		for i=1:molecule.n_atoms
			for k=1:frame.n_atoms
				# Nearest image convention. If the interaction between the probe molecule and atom k is being looked at, we'll only look at the interaction between the probe molecule and the closest replication of atom k. This is done with fractional coordinates for simplication and transformation to cartesian is done later.
				repvec = [nA, nB, nC]
				if abs(fpos[1]+molecule.x[i,1]-(frame.f_coords[k,1]+nA)) > repA/2
					repvec -= [repA,0,0]
				end
				if abs(fpos[2]+molecule.x[i,2]-(frame.f_coords[k,2]+nB)) > repB/2
					repvec -= [0,repB,0]
				end
				if abs(fpos[3]+molecule.x[i,3]-(frame.f_coords[k,3]+nC)) > repC/2
					repvec -= [0,0,repC]
				end
#				println(repvec)
#				println("==========================\n")

				temp = frame.f_to_C*(frame.f_coords[k,:]+repvec)
				r = vecnorm( (pos+molecule.x[i,:]) - frame.f_to_C*(frame.f_coords[k,:]+repvec) )
				σ = ljforcefield.sigmas[ ljforcefield.atom_to_id[ frame.atoms[k] ] , ljforcefield.atom_to_id[ molecule.atoms[i] ] ]
				ϵ = ljforcefield.epsilons[ ljforcefield.atom_to_id[ frame.atoms[k] ] , ljforcefield.atom_to_id[ molecule.atoms[i] ] ]
				if (r < ljforcefield.cutoffradius)
#					@printf("Calling lennard_jones(%f,%f,%f)\n",r,σ,ϵ)
					if lennard_jones(r,σ,ϵ) > 0
#						@printf("%s-%d (nA = %d, nB = %d, nC = %d) and %s-%d -> r = %f | pos = [%f,%f,%f]\n",frame.atoms[k],k,nA,nB,nC,molecule.atoms[i],i,r,pos[1],pos[2],pos[3])
					@printf("pos = [%f,%f,%f], frame.f_to_C*(frame.f_coords[k,:]+repvec) = [%f,%f,%f]\n",pos[1],pos[2],pos[3],temp[1],temp[2],temp[3])
					end
					potsum += lennard_jones(r,σ,ϵ)
				end
			end
		end
	end
	return potsum
end # vdw_energy end

"""
	rotate!(rotMatrix::Array{Float64,2},θ::Float64,ϕ::Float64,ψ::Float64)

Changes the rotation matrix `rotMatrix` to a rotation defined by the angles θ,ϕ,ψ (see http://mathworld.wolfram.com/EulerAngles.html)
"""
function rotate(rotMatrix::Array{Float64,2},θ::Float64,ϕ::Float64,ψ::Float64)
	B = [cos(ψ) sin(ψ) 0; -sin(ψ) cos(ψ) 0; 0 0 1]
	C = [1 0 0; 0 cos(θ) sin(θ); 0 -sin(θ) cos(θ)]
	D = [cos(ϕ) sin(ϕ) 0; -sin(ϕ) cos(ϕ) 0; 0 0 1]
	rotMatrix = B*C*D
end

function exploreframe(frame::Framework, molecule::Molecule, ljforcefield::LennardJonesForceField, repfactors::Array{Int64}; rotation=false, mesh=101)
	if rotation
		rotMatrix = zeros(3,3)
		EnergyMatrix = zeros(mesh,mesh,mesh,mesh)
		# How should I determine what rotations to use?
		rotcnt = 10
	else
		EnergyMatrix = zeros(mesh,mesh,mesh)
		rotcnt = 1
	end
	coordmatrix = Dict{AbstractString,Array{Float64,1}}()
	frac_range_x = linspace(0,repfactors[1],mesh)	
	frac_range_y = linspace(0,repfactors[2],mesh)	
	frac_range_z = linspace(0,repfactors[3],mesh)	

	for (i,xf) in enumerate(frac_range_x), (j,yf) in enumerate(frac_range_y), (k,zf) in enumerate(frac_range_z)
		pos = frame.f_to_C*[xf, yf, zf]
#		@printf("pos = [%f, %f, %f]\n",pos[1],pos[2],pos[3])
#		@printf("i = %d, j = %d, k = %d\n",i,j,k)
		EnergyMatrix[i,j,k] = vdw_energy(frame, molecule, ljforcefield, pos, repfactors)	
		tempstr = @sprintf("%d-%d-%d",i,j,k)
		coordmatrix[tempstr]=pos
	end
	println(typeof(coordmatrix))
	return (EnergyMatrix,coordmatrix)
end
end # end module
