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
struct LJSphere
    species::Symbol
    xf::Array{Float64, 1}
end

function Base.isapprox(ljs1::LJSphere, ljs2::LJSphere)
    return ((ljs1.species == ljs2.species) && isapprox(ljs1.xf, ljs2.xf))
end

"""
Point charge data structure indicates its charge and position in fractional coordinates.

# Example use
    ptc = PtCharge(-0.2, [0.0, 0.0, 0.0])

# Attributes
- `q::Float64`: signed magnitude of charge (units: electrons), e.g. `1.0`
- `xf::Array{Float64, 1}`: fractional coordinates, e.g. `[1.0, 0.0, 4.0]`.
"""
struct PtCharge
    q::Float64
    xf::Array{Float64, 1}
end

function Base.isapprox(c1::PtCharge, c2::PtCharge)
    return (isapprox(c1.q, c2.q) && isapprox(c1.xf, c2.xf))
end
