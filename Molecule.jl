module Mols

export Molecule, ConstructMolecule

type Molecule
	n_atoms::Int64
	atoms::Array{String}
	x::Array{Float64,2}
	charges::Array{Float64}
end

function ConstructMolecule(filename::String)
	f = open(filename,"r")
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

	return Molecule(nAtoms, Atoms, pos)
end # End function

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
	end # End if/if-else/else
end # End function

end # End module
