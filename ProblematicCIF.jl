using PorousMaterials

# Convert to P1

#badfile = "mof3.cif"
#convert_cif_to_P1_symmetry(badfile, "mof3_2.cif")

# Read files + initialize core components

file = "SBMOF-1_cory.cif"
frame = read_crystal_structure_file(file)
strip_numbers_from_atom_labels!(frame)
mol = Molecule(1, ["Xe"], zeros(3,1), [0.0])
ljforcefield = read_forcefield_file("UFF.csv", cutoffradius = 12.5)
reps = replication_factors(frame.box, ljforcefield)

# Define a Snapshot. I define a simulation box where I sample the energy and return it

#snapshot = Snapshot(frame::Framework, molecule::Molecule,
#			ljforcefield::LennardJonesForceField, dimensions::Array{Float64},
#			startpoint::Array{Float64}, repfactors::Tuple{Int64, Int64, Int64},
#			roughness::Float64)
#
#			Forms a Snapshot object to take a snapshot for with dimensions dimensions,
#			starting at startpoint. The mesh is calculated using roughness to measure
#			 distance between points in all dimensions.
									

#snapshot = Snapshot(frame, mol, ljforcefield, [30.,30.,30.], [0.,0.,0.], reps, 0.5)
#energy = snap(snapshot, "energy")
energy = unitcell_energy(frame, mol, ljforcefield, reps, mesh=[100,100,100])

grid = Grid(frame.box, size(energy), energy.*8.3144598, "Float64")

#Write files to visualize with VisIt

replicate_to_xyz(frame, "molecules_sbmof1.xyz", repfactors=reps)
write_unitcell_boundary_vtk(frame, "frame_sbmof1.vtk")
#write_snapshot_to_vtk(snapshot, "snapshot_frame.vtk")
write_to_cube(grid, "energy_sbmof1.cube")
