# Matter is composed of Lennard-Jones spheres and point charges.

"""
Data structure for a Lennard-Jones sphere, containing its species and position in
fractional coordinates.

# Example use
    ljs = LJSphere(:C, [0.0, 0.0, 0.0])

# Attributes
- `species::Symbol`: atom species name, e.g. `:C`
- `xf::Array{Float64, 1}`: fractional coordinates, e.g. `[1.0, 0.0, 4.0]`.
"""
struct Atoms
    n_atoms::Int
    species::Array{Symbol, 1}
    xf::Array{Float64, 2}
end

Atoms(species::Array{Symbol, 1}, xf::Array{Float64, 2}) = Atoms(size(xf, 2), species, xf)

function Base.isapprox(atoms1::Atoms, atoms2::Atoms)
    return all((atoms1.species[i] == atoms2.species[i]) for i = 1:length(atoms1.species)) &&
            isapprox(atoms1.xf, atoms2.xf)
end

"""
Point charge data structure indicates its charge and position in fractional coordinates.

# Example use
    ptc = PtCharge(-0.2, [0.0, 0.0, 0.0])

# Attributes
- `q::Float64`: signed magnitude of charge (units: electrons), e.g. `1.0`
- `xf::Array{Float64, 1}`: fractional coordinates, e.g. `[1.0, 0.0, 4.0]`.
"""
struct Charges
    n_charges::Int
    q::Array{Float64, 1}
    xf::Array{Float64, 2}
end

Charges(q::Array{Float64, 1}, xf::Array{Float64, 2}) = Charges(size(xf, 2), q, xf)

function Base.isapprox(c1::Charges, c2::Charges)
    return (isapprox(c1.q, c2.q) && isapprox(c1.xf, c2.xf))
end
