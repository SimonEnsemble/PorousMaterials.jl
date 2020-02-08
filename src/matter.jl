###
#   Coordinates: Fractional and Cartesian
###
abstract type Coords end

struct Frac<:Coords
    xf::Array{Float64, 2}
end
Base.isapprox(c1::Frac, c2::Frac) = isapprox(c1.xf, c2.xf)
Base.hcat(c1::Frac, c2::Frac) = Frac(hcat(c1.xf, c2.xf))

struct Cart<:Coords
    x::Array{Float64, 2}
end
Base.isapprox(c1::Cart, c2::Cart) = isapprox(c1.x, c2.x)
Base.hcat(c1::Cart, c2::Cart) = Cart(hcat(c1.x, c2.x))


###
#   Atoms
###
# c is either Cart or Frac
struct Atoms{T}
    n::Int # how many atoms?
    species::Array{Symbol, 1} # list of species
    coords::T # coordinates
end
Atoms(species::Array{Symbol, 1}, coords::Coords) = Atoms(length(species), species, coords)

Base.isapprox(a1::Atoms, a2::Atoms) = (a1.species == a2.species) && isapprox(a1.coords, a2.coords)
Base.:+(a1::Atoms, a2::Atoms) = Atoms(a1.n + a2.n, [a1.species; a2.species], hcat(a1.coords, a2.coords))


###
#   Charges
###
struct Charges{T}
    n::Int
    q::Array{Float64, 1}
    coords::T
end
Charges(q::Array{Float64, 1}, coords::Coords) = Charges(length(q), q, coords)

Base.isapprox(c1::Charges, c2::Charges) = isapprox(c1.q, c2.q) && isapprox(c1.coords, c2.coords)
Base.:+(c1::Charges, c2::Charges) = Charges(c1.n + c2.n, [c1.q; c2.q], hcat(c1.coords, c2.coords))
