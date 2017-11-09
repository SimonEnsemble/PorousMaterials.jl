module Mols

export Molecule, constructmolecule, readcharge

"""
	molecule = Molecule(n_atoms::Int64, atoms::Array{String}, x::Array{Float64,2}, charges::Array{Float64})

Construct a molecule used to probe a MOF framework

# Arguments
- `n_atoms::Int64`: Number of atoms in the molecule
- `atoms::Array{String}`: List of element abbreviations
- `x::Array{Float64,2}`: A matrix of cartesian coordinates for atoms in the molecule. x[1,:] is the 3D coordinate for atom 1. The order is the same as in `atoms`
- `charges::Array{Float64}`: An array of charges for each specific atom.
"""
type Molecule
	n_atoms::Int64
	atoms::Array{String}
	x::Array{Float64,2}
	charges::Array{Float64}
end # end Molecule

"""
	Mol = constructmolecule("~/example/filename.mol")

Reads a .mol file and gathers the relevant information to construct a Molecule.
"""
function constructmolecule(mol_filename::String)
	f = open(mol_filename,"r")
	lines = readlines(f)

	nAtoms = 0
	pos = Array{Float64,2}(0,0)
	Atoms = Array{String}(0)
	Charge = Array{Float64}(0)

	for (i,line) in enumerate(lines)
		if (i==4)
			vals = split(line)
			nAtoms = parse(Int64,vals[1])
			pos = Array{Float64,2}(nAtoms,3)
			Atoms = Array{String}(nAtoms)
			Charge = Array{Float64}(nAtoms)
		elseif (i>4 && i<5+nAtoms)
			vals = split(line)
			pos[i-4,:] = map(x->parse(Float64,x), vals[1:3])
			Atoms[i-4] = vals[4]
			Charge[i-4] = readcharge(vals[6])
		end
	end
	close(f)

	return Molecule(nAtoms, Atoms, pos, Charge)
end # end constructmolecule

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

end # End module
