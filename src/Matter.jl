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
mutable struct LJSpheres
    species::Array{Symbol, 1}
    xf::Array{Float64, 2}
end

function LJSpheres(species::Array{AbstractString, 1}, xf::Array{Float64, 2})
    return LJSpheres([Symbol(s) for s in species], xf)
end

function push!(atoms::LJSpheres, new_col::Array{Float64, 1}, new_species::Symbol)
    Base.push!(atoms.species, new_species)
    atoms.xf = [atoms.xf new_col]
end

function Base.isapprox(ljs1::LJSpheres, ljs2::LJSpheres)
    return all((ljs1.species[i] == ljs2.species[i]) for i = 1:length(ljs1.species)) &&
            isapprox(ljs1.xf, ljs2.xf)
end

#function Base.isapprox(ljs1::LJSphere, ljs2::LJSphere)
#    return ((ljs1.species == ljs2.species) && isapprox(ljs1.xf, ljs2.xf))
#end

"""
Point charge data structure indicates its charge and position in fractional coordinates.

# Example use
    ptc = PtCharge(-0.2, [0.0, 0.0, 0.0])

# Attributes
- `q::Float64`: signed magnitude of charge (units: electrons), e.g. `1.0`
- `xf::Array{Float64, 1}`: fractional coordinates, e.g. `[1.0, 0.0, 4.0]`.
"""
mutable struct Charges
    q::Array{Float64, 1}
    xf::Array{Float64, 2}
end
Base

function push!(charges::Charges, new_col::Array{Float64, 1}, new_q::Float64)
    Base.push!(charges.q, new_q)
    charges.xf = [charges.xf new_col]
end

#function Base.isapprox(c1::PtCharge, c2::PtCharge)
#    return (isapprox(c1.q, c2.q) && isapprox(c1.xf, c2.xf))
#end
