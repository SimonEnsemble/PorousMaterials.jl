"""
Data structure holds a set of atom species and their positions in fractional coordinates.

Fractional coords of atom `i` is `charges.xf[:, i]`.

# Example use
    atoms = Atoms(2, [:C, :F], [0.0 1.0; 2.0 3.0; 4.0 5.0])

# Attributes
- `n_atoms::Int`: number of atoms
- `species::Array{Symbol, 1}`: atom species
- `xf::Array{Float64, 2}`: fractional coordinates in the columns
"""
struct Atoms
    n_atoms::Int
    species::Array{Symbol, 1}
    xf::Array{Float64, 2}
end

# compute n_species automatically from array sizes
Atoms(species::Array{Symbol, 1}, xf::Array{Float64, 2}) = Atoms(size(xf, 2), species, xf)

Base.isapprox(a1::Atoms, a2::Atoms) = issetequal(Set(a1.species), Set(a2.species)) &&
    issetequal(Set([round.(a1.xf[:, i], digits=5) for i in 1:a1.n_atoms]), Set([round.(a2.xf[:, i], digits=5) for i in 1:a2.n_atoms])) 

function Base.:+(a1::Atoms, a2::Atoms)
    new_species = [a1.species; a2.species]
    new_xf = [a1.xf a2.xf]

    return Atoms(a1.n_atoms + a2.n_atoms, new_species, new_xf)
end

"""
Data structure holds a set of point charges and their positions in fractional coordinates.

Fractional coords of charge `i` is `charges.xf[:, i]`.

# Example use
    charges = Charges(2, [-1.0, 1.0], [0.0 1.0; 2.0 3.0; 4.0 5.0])

# Attributes
- `n_charges::Int`: number of charges
- `q::Array{Float64, 1}`: signed magnitude of charges (units: electrons)
- `xf::Array{Float64, 2}`: fractional coordinates in the columns
"""
struct Charges
    n_charges::Int
    q::Array{Float64, 1}
    xf::Array{Float64, 2}
end

# compute n_charges automatically from array sizes
Charges(q::Array{Float64, 1}, xf::Array{Float64, 2}) = Charges(size(xf, 2), q, xf)

Base.isapprox(c1::Charges, c2::Charges) = issetequal(Set(c1.q), Set(c2.q)) &&
    issetequal(Set([round.(c1.xf[:, i], digits=5) for i in 1:c1.n_charges]), Set([round.(c2.xf[:, i], digits=5) for i in 1:c2.n_charges]))

function Base.:+(c1::Charges, c2::Charges)
    new_q = [c1.q; c2.q]
    new_xf = [c1.xf c2.xf]
    
    return Charges(c1.n_charges + c2.n_charges, new_q, new_xf)
end
