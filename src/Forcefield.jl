using DataFrames

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

    df = readtable(PATH_TO_DATA * "forcefields/" * filename, allowcomments=true) # from DataFrames
    
    # pure X-X interactions (X = (pseudo)atom)
    pure_sigmas = Dict{AbstractString, Float64}()
    pure_epsilons = Dict{AbstractString, Float64}()
    for row in eachrow(df)
        pure_sigmas[row[:atom]] = row[:sigma]
        pure_epsilons[row[:atom]] = row[:epsilon]
    end
    
    # cross X-Y interactions (X, Y = generally different (pseduo)atoms)
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
end # read_forcefield_file end

"""
	repfactors = replication_factors(framework::Framework, cutoff::Float64)

Find the replication factors needed to make a supercell big enough to fit a sphere with the specified cutoff radius.
In PorousMaterials.jl, rather than replicating the atoms in the home unit cell to build the supercell that
serves as a simulation box, we replicate the home unit cell to form the supercell (simulation box) in a for loop.
This function ensures enough replication factors such that the nearest image convention can be applied.

Returns tuple of replication factors in the a, b, c directions.

# TODO comment on whether it starts at 0 or 1.. like, repfactors = [0, 0, 0] is that possible?
"""
function replication_factors(framework::Framework, cutoff::Float64)
	# Unit vectors used to transform from fractional coordinates to cartesian coordinates. We'll be
	a = framework.f_to_C[:, 1]
	b = framework.f_to_C[:, 2]
	c = framework.f_to_C[:, 3]

	n_ab = cross(a, b)
	n_ac = cross(a, c)
	n_bc = cross(b, c)

	# c0 defines a center in the unit cell
	c0 = [a b c] * [.5, .5, .5]

	rep = [1, 1, 1]

	# Repeat for `a`
	# |n_bc ⋅ c0|/|n_bc| defines the distance from the end of the supercell and the center. As long as that distance is less than the cutoff radius, we need to increase it
	while abs(dot(n_bc, c0)) / vecnorm(n_bc) < cutoff
		rep[1] += 1
		a += framework.f_to_C[:,1]
		c0 = [a b c] * [.5, .5, .5]
	end

	# Repeat for `b`
	while abs(dot(n_ac, c0)) / vecnorm(n_ac) < cutoff
		rep[2] += 1
		b += framework.f_to_C[:,2]
		c0 = [a b c] * [.5, .5, .5]
	end

	# Repeat for `c`
	while abs(dot(n_ab, c0)) / vecnorm(n_ab) < cutoff
		rep[3] += 1
		c += framework.f_to_C[:,3]
		c0 = [a b c] * [.5, .5, .5]
	end

	return (rep[1], rep[2], rep[3])::Tuple{Int, Int, Int}
end # end rep_factors
