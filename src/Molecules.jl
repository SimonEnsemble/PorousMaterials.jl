"""
	molecule = Molecule(n_atoms::Int64, atoms::Array{String}, x::Array{Float64,2}, charges::Array{Float64})

Data structure for a molecule/adsorbate.

# Arguments
- `n_atoms::Int64`: Number of atoms in the molecule
- `atoms::Array{String}`: List of (pseudo)atoms
- `pos::Array{Float64,2}`: A matrix of cartesian coordinates characterizing the position of the atoms in the molecule, stored column-wise so that x[:, i] is the coordinate for atom `atoms[i]`.
- `charges::Array{Float64}`: An array of charges for each specific atom.
"""
type Molecule
	n_atoms::Int64
	atoms::Array{String}
	pos::Array{Float64, 2}
	charges::Array{Float64}
end

"""
	molecule = read_molecule_file("~/example/filename.mol")

Reads a .mol file and gathers the relevant information to construct a Molecule.
# TODO describe the format of the input file
"""
function read_molecule_file(mol_filename::String)
	f = open(mol_filename, "r")
	lines = readlines(f)
	close(f)

	nAtoms = 0
	pos = Array{Float64,2}(0,0)
	Atoms = Array{String}(0)
	Charge = Array{Float64}(0)

	for (i, line) in enumerate(lines)
		if (i == 4)
			vals = split(line)
			nAtoms = parse(Int64, vals[1])
			pos = Array{Float64,2}(3, nAtoms)
			Atoms = Array{String}(nAtoms)
			Charge = Array{Float64}(nAtoms)
		elseif (i > 4 && i < 5 + nAtoms)
			vals = split(line)
			pos[:, i-4] = map(x -> parse(Float64, x), vals[1:3])
			Atoms[i-4] = vals[4]
			Charge[i-4] = readcharge(vals[6])
		end
	end

	return Molecule(nAtoms, Atoms, pos, Charge)
end

"""
	charge::Int64 = readcharge(val::String)

Converts the .mol charge convention into real charges. Grabs a string and returns an integer value
"""
function readcharge(val::String)
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
