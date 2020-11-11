###
#   coordinates: Fractional and Cartesian
###
"""
abstract type for coordinates.
"""
abstract type Coords end

#Base.IndexStyle(::Type{<:Coords}) = IndexLinear()

"""
fractional coordinates, a subtype of `Coords`.

construct by passing an `Array{Float64, 2}` whose columns are the coordinates.

generally, fractional coordinates should be in [0, 1] and are implicitly
associated with a `Box` to represent a periodic coordinate system.

e.g.
```julia
f_coords = Frac(rand(3, 2))  # 2 particles
f_coords.xf                  # retreive fractional coords
```
"""
struct Frac<:Coords
    xf::Array{Float64, 2}
end
Base.isapprox(c1::Frac, c2::Frac; atol::Real=0) = isapprox(c1.xf, c2.xf, atol=atol)
Base.hcat(c1::Frac, c2::Frac) = Frac(hcat(c1.xf, c2.xf))

"""
cartesian coordinates, a subtype of `Coords`.

construct by passing an `Array{Float64, 2}` whose columns are the coordinates.

e.g.
```julia
c_coords = Cart(rand(3, 2))  # 2 particles
c_coords.x                   # retreive cartesian coords
```
"""
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
Base.lastindex(coords::Coords) = length(coords)

Base.setindex!(coords::Frac, val::Array{Float64, 1}, ids) = (coords.xf[:, ids] = val)
Base.setindex!(coords::Cart, val::Array{Float64, 1}, ids) = (coords.x[:, ids] = val)

Base.size(coords::Cart) = size(coords.x)
Base.size(coords::Frac) = size(coords.xf)

origin(T::DataType) = T([0.0, 0.0, 0.0])

"""
    translate_by!(coords, dx)
    translate_by!(coords, dx, box)
    translate_by!(molecule, dx)
    translate_by!(molecule, dx, box)

translate `coords` by the vector `dx`. that is, add the vector `dx`.

this works for any combination of `Frac` and `Cart` coords.

modifies coordinates in place.

`box` is needed when mixing `Frac` and `Cart` coords.

note that periodic boundary conditions are *not* subsequently applied here.

if applied to a `molecule::Molecule`, the coords of atoms, charges, and center of mass
are all translated.
"""
function translate_by!(coords::Cart, dx::Cart)
    coords.x .= broadcast(+, coords.x, dx.x)
end

function translate_by!(coords::Frac, dxf::Frac)
    coords.xf .= broadcast(+, coords.xf, dxf.xf)
end

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
# T is either Cart or Frac
"""
used to represent a set of atoms in space (their atomic species and coordinates).

```julia
struct Atoms{T<:Coords} # enforce that the type specified is `Coords`
    n::Int # how many atoms?
    species::Array{Symbol, 1} # list of species
    coords::T # coordinates
end
```

here, `T` is `Frac` or `Cart`.

helper constructor (infers `n`):
```julia
species = [:H, :H]
coords = Cart(rand(3, 2))
atoms = Atoms(species, coords)
```
"""
struct Atoms{T<:Coords} # enforce that the type specified is `Coords`
    n::Int # how many atoms?
    species::Array{Symbol, 1} # list of species
    coords::T # coordinates
end

function Atoms(species::Array{Symbol, 1}, coords::Coords)
    @assert length(species) == size(coords)[2]
    @assert size(coords)[1] == 3
    return Atoms(length(species), species, coords)
end

Atoms(species::Symbol, coords::Coords) = Atoms([species], coords)
Atoms{Frac}(n::Int) = Atoms([:_ for a = 1:n], Frac([NaN for i = 1:3, a = 1:n])) # safe pre-allocation
Atoms{Cart}(n::Int) = Atoms([:_ for a = 1:n], Cart([NaN for i = 1:3, a = 1:n])) # safe pre-allocation

Base.isapprox(a1::Atoms, a2::Atoms; atol::Real=0) = (a1.species == a2.species) && isapprox(a1.coords, a2.coords, atol=atol)
Base.:+(a1::Atoms, a2::Atoms) = Atoms(a1.n + a2.n, [a1.species; a2.species], hcat(a1.coords, a2.coords))

Base.getindex(atoms::Atoms, ids) = Atoms(atoms.species[ids], atoms.coords[ids])
Base.lastindex(atoms::Atoms) = atoms.n

###
#   point charges
###
"""
used to represent a set of partial point charges in space (their charges and coordinates).

```julia
struct Charges{T<:Coords} # enforce that the type specified is `Coords`
    n::Int
    q::Array{Float64, 1}
    coords::T
end
```

here, `T` is `Frac` or `Cart`.

helper constructor (infers `n`):
```julia
q = [0.1, -0.1]
coords = Cart(rand(3, 2))
charges = Charges(q, coords)
```
"""
struct Charges{T<:Coords} # enforce that the type specified is `Coords`
    n::Int
    q::Array{Float64, 1}
    coords::T
end
function Charges(q::Array{Float64, 1}, coords::Coords)
    @assert length(q) == size(coords)[2]
    @assert size(coords)[1] == 3
    return Charges(length(q), q, coords)
end

Charges(q::Float64, coords::Coords) = Charges([q], coords) # for one charge
Charges{Frac}(n::Int) = Charges([NaN for c = 1:n], Frac([NaN for i = 1:3, c = 1:n])) # safe pre-allocation
Charges{Cart}(n::Int) = Charges([NaN for c = 1:n], Cart([NaN for i = 1:3, c = 1:n])) # safe pre-allocation

Base.isapprox(c1::Charges, c2::Charges; atol::Real=0) = isapprox(c1.q, c2.q) && isapprox(c1.coords, c2.coords, atol=atol)
Base.:+(c1::Charges, c2::Charges) = Charges(c1.n + c2.n, [c1.q; c2.q], hcat(c1.coords, c2.coords))

Base.getindex(charges::Charges, ids) = Charges(charges.q[ids], charges.coords[ids])
Base.lastindex(charges::Charges) = charges.n

"""
    nc = net_charge(charges)
    nc = net_charge(crystal)
    nc = net_charge(molecule)

find the sum of charges in `charges::Charges` or charges in `crystal::Crystal` or `molecule::Molecule`.
(if there are no charges, the net charge is zero.)
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
