"""
Data structure for a Lennard Jones forcefield.

# Attributes
- `name::String`: name of forcefield; correponds to filename
- `pure_σ::Dict{Symbol, Float64}`: Dictionary that returns Lennard-Jones σ of an X-X interaction, where X is an atom. (units: Angstrom)
- `pure_ϵ::Dict{Symbol, Float64}`: Dictionary that returns Lennard-Jones ϵ of an X-X interaction, where X is an atom. (units: K)
- `σ²::Dict{Symbol, Dict{Symbol, Float64}}`: Lennard Jones σ² (units: Angstrom²) for cross-interactions. Example use is `sigmas_squared[:He][:C]`
- `ϵ::Dict{Symbol, Dict{Symbol, Float64}}`: Lennard Jones ϵ (units: K) for cross-interactions. Example use is `epsilons[:He][:C]`
- `cutoffradius_squared::Float64`: The square of the cut-off radius beyond which we define the potential energy to be zero (units: Angstrom²). We store σ² to speed up computations, which involve σ², not σ.
"""
struct LJForceField
    name::String

	pure_σ::Dict{Symbol, Float64}
	pure_ϵ::Dict{Symbol, Float64}

	σ²::Dict{Symbol, Dict{Symbol, Float64}}
	ϵ::Dict{Symbol, Dict{Symbol, Float64}}

	cutoffradius_squared::Float64
end

Base.broadcastable(ljff::LJForceField) = Ref(ljff)

"""
	ljforcefield = ForceField(forcefieldfile; cutoffradius=14.0, mixing_rules="Lorentz-Berthelot")

Read a .csv file containing Lennard Jones parameters (with the following columns: `atom,sigma,epsilon` and constructs a LJForceField object.

The following mixing rules are implemented:
* Kong mixing rules: DOI 10.1063/1.1680358
* Lorentz-Berthelot: https://en.wikipedia.org/wiki/Combining_rules#Lorentz-Berthelot_rules
* Geometric

# Arguments
- `forcefieldfile::AbstractString`: Name of the forcefield csv file
- `cutoffradius::Float64`: Cutoff radius beyond which we define the potential energy to be zero (units: Angstrom)
- `mixing_rules::AbstractString`: The mixing rules used to compute the cross-interaction terms of the forcefield

# Returns
- `ljforcefield::LJForceField`: The data structure containing the forcefield parameters (pure σ, ϵ and cross interaction terms as well)
"""
function LJForceField(forcefieldfile::AbstractString; cutoffradius::Float64=14.0,
                      mixing_rules::AbstractString="Lorentz-Berthelot")
    if ! (lowercase(mixing_rules) in ["lorentz-berthelot", "kong", "geometric"])
        error(@sprintf("%s mixing rules not implemented...\n", mixing_rules))
    end

    forcefield_file_path = joinpath(PATH_TO_DATA, "forcefields", forcefieldfile)

    df = CSV.read(forcefield_file_path) # from DataFrames

    ljff = LJForceField(forcefieldfile, Dict(), Dict(), Dict(), Dict(), cutoffradius ^ 2)

    # pure X-X interactions (X = (pseudo)atom)
    for row in eachrow(df)
        atom_species = Symbol(row[:atom])
        # if atom already recorded, we have a duplicate. this is dangerous to overwrite.
        if atom_species in keys(ljff.pure_σ)
            error(@sprintf("Atom %s listed at least twice in %s.\n", atom_species,
                forcefield_file_path))
        end
		ljff.pure_σ[atom_species] = row[Symbol("sigma(A)")]
		ljff.pure_ϵ[atom_species] = row[Symbol("epsilon(K)")]
    end

    # cross X-Y interactions (X, Y = generally different (pseduo)atoms)
	for atom in [Symbol(atom) for atom in keys(ljff.pure_σ)]
        ljff.ϵ[atom] = Dict{Symbol, Float64}()
        ljff.σ²[atom] = Dict{Symbol, Float64}()
		for other_atom in [Symbol(other_atom) for other_atom in keys(ljff.pure_σ)]
            if lowercase(mixing_rules) == "lorentz-berthelot"
                ϵ_ij = sqrt(ljff.pure_ϵ[atom] * ljff.pure_ϵ[other_atom])
                σ_ij² = ((ljff.pure_σ[atom] + ljff.pure_σ[other_atom]) / 2.0) ^ 2
            elseif lowercase(mixing_rules) == "kong"
                ϵ_iiσ_ii⁶ = ljff.pure_ϵ[atom]       * ljff.pure_σ[atom] ^ 6
                ϵ_jjσ_jj⁶ = ljff.pure_ϵ[other_atom] * ljff.pure_σ[other_atom] ^ 6

                ϵ_iiσ_ii¹² = ϵ_iiσ_ii⁶ * ljff.pure_σ[atom] ^ 6
                ϵ_jjσ_jj¹² = ϵ_jjσ_jj⁶ * ljff.pure_σ[other_atom] ^ 6

                ϵ_ijσ_ij⁶ = sqrt(ϵ_iiσ_ii⁶ * ϵ_jjσ_jj⁶)
                ϵ_ijσ_ij¹² = ((ϵ_iiσ_ii¹² ^ (1/13) + ϵ_jjσ_jj¹² ^ (1/13)) / 2) ^ 13

                ϵ_ij = ϵ_ijσ_ij⁶ ^ 2 / ϵ_ijσ_ij¹²
                σ_ij² = (ϵ_ijσ_ij¹² / ϵ_ijσ_ij⁶) ^ (1/3)
            elseif lowercase(mixing_rules) == "geometric"
                ϵ_ij = sqrt(ljff.pure_ϵ[atom] * ljff.pure_ϵ[other_atom])
                σ_ij² = ljff.pure_σ[atom] * ljff.pure_σ[other_atom] # √(σ_i σ_j)²
            end
            ljff.ϵ[atom][other_atom] = ϵ_ij
            ljff.σ²[atom][other_atom] = σ_ij²
		end
	end

	return ljff
end

"""
	repfactors = replication_factors(unitcell, cutoffradius)

Find the replication factors needed to make a supercell big enough to fit a sphere with the specified cutoff radius.
In PorousMaterials.jl, rather than replicating the atoms in the home unit cell to build the supercell that
serves as a simulation box, we replicate the home unit cell to form the supercell (simulation box) in a for loop.
This function ensures enough replication factors such that the nearest image convention can be applied.

A non-replicated supercell has 1 as the replication factor in each dimension (`repfactors = (1, 1, 1)`).

# Arguments
- `unitcell::Box`: The unit cell of the framework
- `cutoff_radius::Float64`: Cutoff radius beyond which we define the potential energy to be zero (units: Angstrom)

# Returns
- `repfactors::Tuple{Int, Int, Int}`: The replication factors in the a, b, c directions
"""
function replication_factors(unitcell::Box, cutoff_radius::Float64)
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

	# Repeat for `a`
	# |n_bc ⋅ c0|/|n_bc| defines the distance from the end of the supercell and the center. As long as that distance is less than the cutoff radius, we need to increase it
	while abs(dot(n_bc, c0)) / norm(n_bc) < cutoff_radius
		rep[1] += 1
		a += unitcell.f_to_c[:,1]
		c0 = [a b c] * [.5, .5, .5]
	end

	# Repeat for `b`
	while abs(dot(n_ac, c0)) / norm(n_ac) < cutoff_radius
		rep[2] += 1
		b += unitcell.f_to_c[:,2]
		c0 = [a b c] * [.5, .5, .5]
	end

	# Repeat for `c`
	while abs(dot(n_ab, c0)) / norm(n_ab) < cutoff_radius
		rep[3] += 1
		c += unitcell.f_to_c[:,3]
		c0 = [a b c] * [.5, .5, .5]
	end

	return (rep[1], rep[2], rep[3])::Tuple{Int, Int, Int}
end

replication_factors(unitcell::Box, ljforcefield::LJForceField) = replication_factors(unitcell, sqrt(ljforcefield.cutoffradius_squared))
replication_factors(framework::Framework, cutoff_radius::Float64) = replication_factors(framework.box, cutoff_radius)
replication_factors(framework::Framework, ljforcefield::LJForceField) = replication_factors(framework.box, sqrt(ljforcefield.cutoffradius_squared))

"""
    check_forcefield_coverage(framework, ljforcefield)
    check_forcefield_coverage(molecule, ljforcefield)

Check that the force field contains parameters for every atom present in a framework or molecule.
Will print out which atoms are missing.

# Arguments
- `framework::Framework`: The framework containing the crystal structure information
- `molecule::Molecule`: A molecule object
- `ljforcefield::LJForceField`: A Lennard Jones forcefield object containing information on atom interactions

# Returns
- `all_covered::Bool`: Returns true if all atoms in the `framework` are also included in `ljforcefield`. False otherwise
"""
function check_forcefield_coverage(f_or_m::Union{Framework, Molecule}, ljff::LJForceField)
    unique_species = unique(f_or_m.atoms.species)
    all_covered = true
    for species in unique_species
        if !(species in keys(ljff.pure_ϵ))
            @warn @sprintf("\t%s in %s missing from %s force field.", species,
                isa(f_or_m, Framework) ? f_or_m.name : f_or_m.species, ljff.name)
            all_covered = false
        end
    end
    return all_covered
end

function Base.show(io::IO, ff::LJForceField)
    println(io, "Force field: ", ff.name)
	println(io, "Number of atoms included: ", length(ff.pure_σ))
	println(io, "Cut-off radius (Å) = ", sqrt(ff.cutoffradius_squared))
    for atom in keys(ff.pure_σ)
        @printf(io, "%5s-%5s ϵ = %10.5f K, σ = %10.5f Å\n", atom, atom, ff.pure_ϵ[atom], ff.pure_σ[atom])
    end
end
