using PorousMaterials

@printf("Starting script:\n")
files = readdir("data/crystals")

ljforcefield = read_forcefield_file("UFF.csv", cutoffradius = 12.5)

@printf("Iterating through %d crystal file structures\n",length(files))

for filename in files[1:1]
    frame = read_crystal_structure_file(filename)
    strip_numbers_from_atom_labels!(frame)
    mol = Molecule(1, ["He"], zeros(3,1), [0.0])

    reps = replication_factors(frame.box, ljforcefield)

    snapshot = Snapshot(frame, mol, ljforcefield, [15., 15., 15.], [0., 0., 0.], reps, 0.4)
    occ = snap(snapshot, "occupancy")
#    write_to_npy(occ)
#	write_snapshot_to_vtk(snapshot, "snapp.vtk")
	print(occ)
end
