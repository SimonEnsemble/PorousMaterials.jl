using PorousMaterials

@printf("Starting script:\n")
filename = "SBMOF-1_cory.cif"

ljforcefield = read_forcefield_file("UFF.csv", cutoffradius = 12.5)


frame = read_crystal_structure_file(filename)
strip_numbers_from_atom_labels!(frame)
mol = Molecule(1, ["He"], zeros(3,1), [0.0])

reps = replication_factors(frame.box, ljforcefield)

snapshot = Snapshot(frame, mol, ljforcefield, [10., 10., 10.], [0., 0., 0.], reps, 1)
occ = snap(snapshot, "occupancy")
write_to_npy(occ, split(filename, ".")[1] * "imgposter.npy")
