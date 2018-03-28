cd("../PorousMaterials.jl")
using PorousMaterials

@printf("Starting script:\n")
files = readdir("data/crystals")

ljforcefield = read_forcefield_file("UFF.csv", cutoffradius = 12.5)

@printf("Iterating through %d crystal file structures\n",length(files))

for filename in files[66:67]
	try
		frame = read_crystal_structure_file(filename)
		strip_numbers_from_atom_labels!(frame)
		mol = Molecule(1, ["He"], zeros(3,1), [0.0])

		reps = replication_factors(frame.box, ljforcefield)

		snapshot = Snapshot(frame, mol, ljforcefield, [10., 10., 10.], [0., 0., 0.], reps, 1)
		occ = snap(snapshot, "occupancy")
		write_to_npy(occ, split(filename, ".")[1] * "drasl_r02.npy")
	catch e
		@printf("%s caused the following error:\n%s\n",filename,e)
	end
end
