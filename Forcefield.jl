"""All things energy related"""
module Forcefield

using DataFrames
using Crystal
using Mols

export LennardJonesForceField, read_forcefield_file, lennard_jones, readproperties, rep_factors, centerofmass, vdw_energy, rotate, exploreframe

"""
	ljforcefield = LennardJonesForceField(cutoffradius, epsilon_dict, sigma_dict, atom_to_id, epsilons, sigmas)

Data structure for a Lennard Jones forcefield, read from a file containing UFF parameters.

# Arguments
- `pure_sigmas::Dict{AbstractString, Float64}`: Dictionary that connects element acronyms to a σ, which is the finite distance where the potential between atoms goes to zero
- `pure_epsilons::Dict{AbstractString, Float64}`: Dictionary that connects element acronyms to an ϵ, which is the depth of a Lennard Jones potential well
- `epsilons::Dict{AbstractString, Dict{AbstractString, Float64}}`: Lennard Jones ϵ (units: K) for cross-interactions. Example use is `epsilons["He"]["C"]`
- `sigmas::Dict{AbstractString, Dict{AbstractString, Float64}}`: Lennard Jones σ (units: A) for cross-interactions. Example use is `sigmas["He"]["C"]`
- `cutoffradius::Float64`: cut-off radius beyond which we define the potential energy to be zero (units: Angstrom)
"""
struct LennardJonesForceField
	pure_sigmas::Dict{AbstractString, Float64}
	pure_epsilons::Dict{AbstractString, Float64}

	sigmas::Dict{AbstractString, Dict{AbstractString, Float64}}
	epsilons::Dict{AbstractString, Dict{AbstractString, Float64}}

	cutoffradius::Float64
end

"""
	ljforcefield = read_forcefield_file("filename.csv")

Read a .csv file containing Lennard Jones parameters (with the following columns: `atom,sigma,epsilon` and constructs a LennardJonesForceField object.
"""
function read_forcefield_file(filename::AbstractString; cutoffradius::Float64=14.0, mixing_rules::AbstractString="Lorentz-Berthelot")
    if ! (mixing_rules in ["Lorentz-Berthelot"])
        error(@sprintf("%s mixing rules not implemented...\n", mixing_rules))
    end

    df = readtable(filename, allowcomments=true)

    pure_sigmas = Dict{AbstractString, Float64}()
    pure_epsilons = Dict{AbstractString, Float64}()
    for row in eachrow(df)
        pure_sigmas[row[:atom]] = row[:sigma]
        pure_epsilons[row[:atom]] = row[:epsilon]
    end
    
    # cross interactions
    epsilons = Dict{AbstractString, Dict{AbstractString, Float64}}()
    sigmas = Dict{AbstractString, Dict{AbstractString, Float64}}()
	for atom in keys(pure_sigmas)
        epsilons[atom] = Dict{AbstractString, Float64}()
        sigmas[atom] = Dict{AbstractString, Float64}()
        for other_atom in keys(pure_sigmas)
			epsilons[atom][other_atom] = sqrt(pure_epsilons[atom] * pure_epsilons[other_atom])
			sigmas[atom][other_atom] = (pure_sigmas[atom] + pure_sigmas[other_atom]) / 2.0
		end
	end

	return LennardJonesForceField(pure_sigmas, pure_epsilons, sigmas, epsilons, cutoffradius)
end # constructforcefield end

"""
	V = lennard_jones_potential_energy(r::Float64, σ::Float64, ϵ::Float64) # returns units Kelvin

Calculate the lennard jones potential energy given a radius r between two molecules. σ and ϵ are specific to interaction between two elements.
returns potential energy in units Kelvin.
# Arguments
- `r::Float64`: distance between two (pseudo)atoms in question (Angstrom)
- `σ::Float64`: sigma parameter in Lennard Jones potential (units: Angstrom)
- `ϵ::Float64`: epsilon parameter in Lennard Jones potential (units: Kelvin)
"""
function lennard_jones(r::Float64, σ::Float64, ϵ::Float64)
	ratio = (σ / r) ^ 6
	return 4 * ϵ * (ratio ^ 2 - ratio)
end

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
	eleprops = Dict{AbstractString,Array{Float64}}()
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
			eleprops[str[1]] = temp
		end
	end
	close(f)
	return eleprops
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

	eleprops = readproperties(filename)
	rvec = zeros(3)
	mtot = 0.0
	for nA=1:repA, nB=1:repB, nC=1:repC
		for i=1:frame.n_atoms
			rvec += (frame.f_coords[:,i]+[nA-1, nB-1, nC-1])*eleprops[frame.atoms[i]][1]
			mtot += eleprops[frame.atoms[i]][1]
		end
	end
	return frame.f_to_C*(rvec/mtot)
end


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
	coordmatrix = Dict{AbstractString,Array{Float64,2}}()
	frac_range_x = linspace(0,repfactors[1],mesh)	
	frac_range_y = linspace(0,repfactors[2],mesh)	
	frac_range_z = linspace(0,repfactors[3],mesh)	

	for (i,xf) in enumerate(frac_range_x), (j,yf) in enumerate(frac_range_y), (k,zf) in enumerate(frac_range_z)
		molecule.pos = (frame.f_to_C*[xf,yf,zf])[:,:]
#		@printf("pos = [%f, %f, %f]\n",pos[1],pos[2],pos[3])
#		@printf("i = %d, j = %d, k = %d\n",i,j,k)
		EnergyMatrix[i,j,k] = vdw_energy(frame, molecule, ljforcefield, repfactors)	
		tempstr = @sprintf("%d-%d-%d",i,j,k)
		coordmatrix[tempstr] = molecule.pos
	end
	println(typeof(coordmatrix))
	return (EnergyMatrix,coordmatrix)
end
end # end module
