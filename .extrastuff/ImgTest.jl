using PorousMaterials

filename = "SBMOF-1_cory.cif"

ljforcefield = read_forcefield_file("UFF.csv", cutoffradius = 12.5)

frame = read_crystal_structure_file(filename)
strip_numbers_from_atom_labels!(frame)

mol = Molecule(1, ["He"], zeros(3,1), [0.0])

reps = replication_factors(frame.box, ljforcefield)

snapshot = Snapshot(frame, mol, ljforcefield, [30., 30., 30.], [0., 0., 0.], reps, 0.5)
#energy = snap(snapshot, "energy")
energy = fractional_energy(snapshot)

replicate_to_xyz(frame,"imgposter.xyz", repfactors = reps)

grid = Grid(frame.box, (75, 75, 75), energy, "Float64", snapshot.repfactors)

write_to_cube(grid, "imgposter2.cube")
