"""
	molecule = Molecule(n_atoms::Int64, atoms::Array{String}, x::Array{Float64,2}, charges::Array{Float64})

Data structure for a molecule/adsorbate.

# Arguments
- `n_atoms::Int64`: Number of atoms in the molecule
- `atoms::Array{AbstractString}`: List of (pseudo)atoms
- `x::Array{Float64,2}`: A matrix of cartesian coordinates characterizing the position of the atoms in the molecule, stored column-wise so that x[:, i] is the coordinate for atom `atoms[i]`.
- `charges::Array{Float64}`: An array of charges for each specific atom.
"""
mutable struct Molecule
	n_atoms::Int64
	atoms::Array{AbstractString}
	x::Array{Float64, 2}
	charges::Array{Float64}
end

import Base.print
function print(io::IO, molecule::Molecule)
	println("Number of atoms in molecule: ",molecule.n_atoms)
	print("Position of atoms: ")
	for i=1:molecule.n_atoms
		@printf("\n%s:\t[%.3f, %.3f, %.3f]", molecule.atoms[i], molecule.x[1,i], molecule.x[2,i], molecule.x[3,i])
	end
end

import Base.show
function show(io::IO, molecule::Molecule)
	print(io, molecule)
end
