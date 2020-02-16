###
#   Coordinates: Fractional and Cartesian
###
abstract type Coords end

struct Frac<:Coords
    xf::Array{Float64, 2}
end
Base.isapprox(c1::Frac, c2::Frac; atol::Real=0) = isapprox(c1.xf, c2.xf, atol=atol)
Base.hcat(c1::Frac, c2::Frac) = Frac(hcat(c1.xf, c2.xf))

struct Cart<:Coords
    x::Array{Float64, 2}
end
Base.isapprox(c1::Cart, c2::Cart; atol::Real=0) = isapprox(c1.x, c2.x, atol=atol)
Base.hcat(c1::Cart, c2::Cart) = Cart(hcat(c1.x, c2.x))

# helpers for single vectors
Frac(xf::Array{Float64, 1}) = Frac(reshape(xf, (3, 1)))
Cart(x::Array{Float64, 1}) = Cart(reshape(x, (3, 1)))

Base.getindex(coords::Frac, ids) = Frac(coords.xf[:, ids])
Base.getindex(coords::Cart, ids) = Cart(coords.x[:, ids])

Base.length(coords::Cart) = size(coords.x, 2)
Base.length(coords::Frac) = size(coords.xf, 2)

"""
    wrap!(f::Frac)
    wrap!(crystal::Crystal)

wrap fractional coordinates to [0, 1] via `mod(⋅, 1.0)`.
e.g. -0.1 --> 0.9 and 1.1 -> 0.1
"""
function wrap!(f::Frac)
    f.xf .= mod.(f.xf, 1.0)
end

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
Atoms{Frac}(n::Int) = Atoms([:_ for a = 1:n], Frac([NaN for i = 1:3, a = 1:n])) # safe pre-allocation
Atoms{Cart}(n::Int) = Atoms([:_ for a = 1:n], Cart([NaN for i = 1:3, a = 1:n])) # safe pre-allocation

Base.isapprox(a1::Atoms, a2::Atoms; atol::Real=0) = (a1.species == a2.species) && isapprox(a1.coords, a2.coords, atol=atol)
Base.:+(a1::Atoms, a2::Atoms) = Atoms(a1.n + a2.n, [a1.species; a2.species], hcat(a1.coords, a2.coords))

Base.getindex(atoms::Atoms, ids) = Atoms(atoms.species[ids], atoms.coords[ids])

###
#   Charges
###
struct Charges{T}
    n::Int
    q::Array{Float64, 1}
    coords::T
end
Charges(q::Array{Float64, 1}, coords::Coords) = Charges(length(q), q, coords)
Charges{Frac}(n::Int) = Charges([NaN for c = 1:n], Frac([NaN for i = 1:3, c = 1:n])) # safe pre-allocation
Charges{Cart}(n::Int) = Charges([NaN for c = 1:n], Cart([NaN for i = 1:3, c = 1:n])) # safe pre-allocation

Base.isapprox(c1::Charges, c2::Charges; atol::Real=0) = isapprox(c1.q, c2.q) && isapprox(c1.coords, c2.coords, atol=atol)
Base.:+(c1::Charges, c2::Charges) = Charges(c1.n + c2.n, [c1.q; c2.q], hcat(c1.coords, c2.coords))

Base.getindex(charges::Charges, ids) = Charges(charges.q[ids], charges.coords[ids])

"""
    nc = net_charge(charges)
    nc = net_charge(crystal)

find the sum of charges in `charges::Charges` or `crystal.charges` where `crystal::Crystal`.
"""
function net_charge(charges::Charges)
    if charges.n == 0
        return 0.0
    else
        return sum(charges.q)
    end
end

"""
    neutral(charges, tol) # true or false. default tol = 1e-5
    neutral(crystal, tol) # true or false. default tol = 1e-5

determine if a set of `charges::Charges` (`charges.q`) sum to an absolute value less than `tol::Float64`. if `crystal::Crystal` is passed, the function looks at the `crystal.charges`. i.e. determine the absolute value of the net charge is less than `tol`.
"""
neutral(charges::Charges, tol::Float64=1e-5) = abs(net_charge(charges)) < tol

"""
    multiequal_sets(atoms_i, atoms_j, digits=5) # true or false
    multiequal_sets(charges_i, charges_j, digits=5) # true or false

Are two multisets of atoms/charges equal? Think of this as a permutation-invariant `isapprox`.
`digits` is up to how many digits we consider when considering the coordinates/charges approximately equal.
"""
function equal_sets(a1::Atoms{Frac}, a2::Atoms{Frac}; digits::Int=5)
    if a1.n != a2.n
        return false
    end
    a1_s = Set()
    a2_s = Set()
    for a = 1:a1.n
        push!(a1_s, (a1.species[a], round.(a1.coords.xf[:, a], digits=digits)))
        push!(a2_s, (a2.species[a], round.(a2.coords.xf[:, a], digits=digits)))
    end
    return a1_s == a2_s
end

function equal_sets(c1::Charges{Frac}, c2::Charges{Frac}; digits::Int=5)
    if c1.n != c2.n
        return false
    end
    c1_s = Set()
    c2_s = Set()
    for c = 1:c1.n
        push!(c1_s, (round(c1.q[c], digits=digits), round.(c1.coords.xf[:, c], digits=digits)))
        push!(c2_s, (round(c2.q[c], digits=digits), round.(c2.coords.xf[:, c], digits=digits)))
    end
    return c1_s == c2_s
end
