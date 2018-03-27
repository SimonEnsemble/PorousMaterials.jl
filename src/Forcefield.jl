using DataFrames
using CSV

"""
Data structure for a Lennard Jones forcefield.

# Attributes
- `name::String`: name of forcefield; correponds to filename
- `pure_σ::Dict{Symbol, Float64}`: Dictionary that returns Lennard-Jones σ of an X-X interaction, where X is an atom. (units: Angstrom)
- `pure_ϵ::Dict{Symbol, Float64}`: Dictionary that returns Lennard-Jones ϵ of an X-X interaction, where X is an atom. (units: K)
- `ϵ::Dict{Symbol, Dict{Symbol, Float64}}`: Lennard Jones ϵ (units: K) for cross-interactions. Example use is `epsilons[:He][:C]`
- `σ²::Dict{Symbol, Dict{Symbol, Float64}}`: Lennard Jones σ² (units: Angstrom²) for cross-interactions. Example use is `sigmas_squared[:He][:C]`
- `cutoffradius_squared::Float64`: The square of the cut-off radius beyond which we define the potential energy to be zero (units: Angstrom²). We store σ² to speed up computations, which involve σ², not σ.
"""
struct LennardJonesForceField
    name::String
    
	pure_σ::Dict{Symbol, Float64}
	pure_ϵ::Dict{Symbol, Float64}

	σ²::Dict{Symbol, Dict{Symbol, Float64}}
	ϵ::Dict{Symbol, Dict{Symbol, Float64}}

	cutoffradius_squared::Float64
end


"""
	ljforcefield = read_forcefield_file("filename.csv", cutoffradius::Float64=14.0, mixing_rules::AbstractString="Lorentz-Berthelot")

Read a .csv file containing Lennard Jones parameters (with the following columns: `atom,sigma,epsilon` and constructs a LennardJonesForceField object.
"""
function read_forcefield_file(filename::AbstractString; cutoffradius::Float64=14.0, mixing_rules::AbstractString="Lorentz-Berthelot")
    if ! (mixing_rules in ["Lorentz-Berthelot"])
        error(@sprintf("%s mixing rules not implemented...\n", mixing_rules))
    end

    df = CSV.read(PATH_TO_DATA * "forcefields/" * filename) # from DataFrames
    # assert that all atoms in the force field are unique (i.e. no duplicates)
    @assert(length(unique(df[:atom])) == size(df, 1), 
        @sprintf("Duplicate atoms found in force field file %s\n", filename))
    
    ljff = LennardJonesForceField(filename, Dict(), Dict(), Dict(), Dict(), cutoffradius ^ 2)
    
    # pure X-X interactions (X = (pseudo)atom)
    for row in eachrow(df)
		ljff.pure_σ[Symbol(row[:atom])] = row[Symbol("sigma(A)")]
		ljff.pure_ϵ[Symbol(row[:atom])] = row[Symbol("epsilon(K)")]
    end
    
    # cross X-Y interactions (X, Y = generally different (pseduo)atoms)
	for atom in [Symbol(atom) for atom in keys(ljff.pure_σ)]
        ljff.ϵ[atom] = Dict{Symbol, Float64}()
        ljff.σ²[atom] = Dict{Symbol, Float64}()
		for other_atom in [Symbol(other_atom) for other_atom in keys(ljff.pure_σ)]
			ljff.ϵ[atom][other_atom] = sqrt(ljff.pure_ϵ[atom] * ljff.pure_ϵ[other_atom])
			ljff.σ²[atom][other_atom] = ((ljff.pure_σ[atom] + ljff.pure_σ[other_atom]) / 2.0) ^ 2
		end
	end

	return ljff
end

"""
	repfactors = replication_factors(unitcell::Box, cutoffradius::Float64)

Find the replication factors needed to make a supercell big enough to fit a sphere with the specified cutoff radius.
In PorousMaterials.jl, rather than replicating the atoms in the home unit cell to build the supercell that
serves as a simulation box, we replicate the home unit cell to form the supercell (simulation box) in a for loop.
This function ensures enough replication factors such that the nearest image convention can be applied.

Returns tuple of replication factors in the a, b, c directions.

A non-replicated supercell has 1 as the replication factor in each dimension (`repfactors = [1, 1, 1]`).
#TODO comment on whether it starts at 0 or 1.. like, repfactors = [0, 0, 0] is that possible?
"""
function replication_factors(unitcell::Box, ljforcefield::LennardJonesForceField)
	# Unit vectors used to transform from fractional coordinates to cartesian coordinates. We'll be
	a = unitcell.f_to_c[:, 1]
	b = unitcell.f_to_c[:, 2]
	c = unitcell.f_to_c[:, 3]

	n_ab = cross(a, b)
	n_ac = cross(a, c)
	n_bc = cross(b, c)

	# c0 defines a center in the unit cell
	c0 = [a b c] * [.5, .5, .5]

	rep = [1, 1, 1]
	cutoff = sqrt(ljforcefield.cutoffradius_squared)

	# Repeat for `a`
	# |n_bc ⋅ c0|/|n_bc| defines the distance from the end of the supercell and the center. As long as that distance is less than the cutoff radius, we need to increase it
	while abs(dot(n_bc, c0)) / vecnorm(n_bc) < cutoff
		rep[1] += 1
		a += unitcell.f_to_c[:,1]
		c0 = [a b c] * [.5, .5, .5]
	end

	# Repeat for `b`
	while abs(dot(n_ac, c0)) / vecnorm(n_ac) < cutoff
		rep[2] += 1
		b += unitcell.f_to_c[:,2]
		c0 = [a b c] * [.5, .5, .5]
	end

	# Repeat for `c`
	while abs(dot(n_ab, c0)) / vecnorm(n_ab) < cutoff
		rep[3] += 1
		c += unitcell.f_to_c[:,3]
		c0 = [a b c] * [.5, .5, .5]
	end

	return (rep[1], rep[2], rep[3])::Tuple{Int64, Int64, Int64}
end # end rep_factors

"""
    check_forcefield_coverage(framework::Framework, ljforcefield::LennardJonesForceField; verbose::Bool=true)

Check that the force field contains parameters for every atom present in the framework.
returns true or false; prints which atoms are missing by default if `verbose=true`.
"""
function check_forcefield_coverage(framework::Framework, ljforcefield::LennardJonesForceField; verbose::Bool=true)
    framework_atoms = unique(framework.atoms)
    forcefield_atoms = keys(ljforcefield.pure_ϵ)

    full_coverage = true

    for atom in framework_atoms
        if !(atom in forcefield_atoms)
            if verbose
                @printf("%s framework atom \"%s\" is not covered by the forcefield %s.\n", framework.name, atom, ljforcefield.name)
            end
            full_coverage = false
        end
    end
    return full_coverage
end

import Base.print
function print(io::IO, ljforcefield::LennardJonesForceField)
	println(io, "Amount of atoms included: ",length(ljforcefield.pure_σ))
	print(io, "Cut-off radius = ",sqrt(ljforcefield.cutoffradius_squared))
end

import Base.show
function show(io::IO, ljforcefield::LennardJonesForceField) 
	print(io, ljforcefield)
end
