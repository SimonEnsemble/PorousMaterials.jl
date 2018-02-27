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

"""
	molecule = read_molecule_file("~/example/filename.mol")

Reads a .mol file and gathers the relevant information to construct a Molecule.
# TODO describe the format of the input file
"""
function read_molecule_file(mol_filename::AbstractString)
	extension = split(mol_filename, ".")[end]
	if ! (extension in ["mol"])
		error("PorousMaterials.jl can only read .mol molecule files.")
	end

	f = open(mol_filename, "r")
	lines = readlines(f)
	close(f)

	nAtoms = 0
    # TODO n_atoms, not nAtoms, no capital letters in var names = convention
	pos = Array{Float64,2}(0,0)
	Atoms = Array{AbstractString}(0)
	Charge = Array{Float64}(0)

	for (i, line) in enumerate(lines)
		if (i == 4)
			vals = split(line)
			nAtoms = parse(Int64, vals[1])
			pos = Array{Float64,2}(3, nAtoms)
			Atoms = Array{AbstractString}(nAtoms)
			Charge = Array{Float64}(nAtoms)
		elseif (i > 4 && i < 5 + nAtoms)
			vals = split(line)
			pos[:, i-4] = map(x -> parse(Float64, x), vals[1:3])
			Atoms[i-4] = vals[4]
            # TODO Atoms --> atoms, Charge --> charge. Capital letters reserved for structs, types!
			Charge[i-4] = readcharge(vals[6])
		end
	end

	return Molecule(nAtoms, Atoms, pos, Charge)
end

"""
	charge::Int64 = readcharge(val::String)

Converts the .mol charge convention into real charges. Grabs a string and returns an integer value
"""
function readcharge(val::AbstractString)
	num = parse(Int64,val)
	if (num == 7)
		return -3
	elseif (num == 6)
		return -2
	elseif (num == 5)
		return -1
	elseif (num == 0)
		return 0
	elseif (num == 3)
		return 1
	elseif (num == 2)
		return 2
	elseif (num == 1)
		return 3
	else
		error("Charge not valid")
	end # end if/if-else/else
end # end readcharge
# TODO remove this? what is this doing? charges can be floats... so why not just read in 
#  the floats that are the charges assigned to that bead of the molecule?
