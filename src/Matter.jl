# Matter is composed of Lennard-Jones spheres and point charges.

"""
Lennard-Jones sphere data structure indicates its species and its Cartesian coordinates.

# Attributes
- `atom::Symbol`: atom species name, e.g. `:C`
- `x::Array{Float64, 1}`: Cartesian coordinates (units: Angstrom), e.g. `[1.0, 0.0, 4.0]`.
"""
struct LJSphere
    atom::Symbol
    x::Array{Float64, 1}
end

function Base.isapprox(ljs1::LJSphere, ljs2::LJSphere)
    return ((ljs1.atom == ljs2.atom) && isapprox(ljs1.x, ljs2.x))
end

"""
Point charge data structure indicates its charge and its position in Cartesian coordinates.

# Attributes
- `q::Float64`: signed magnitude of charge (units: electrons)
- `x::Array{Float64, 1}`: Cartesian coordinates (units: Angstrom), e.g. `[1.0, 0.0, 4.0]`.
"""
struct PointCharge
    q::Float64
    x::Array{Float64, 1}
end

function Base.isapprox(c1::PointCharge, c2::PointCharge)
    return (isapprox(c1.q, c2.q) && isapprox(c1.x, c2.x))
end
