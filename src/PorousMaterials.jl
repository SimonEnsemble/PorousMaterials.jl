module PorousMaterials

# this is the directory where crystal structures, forcefields, and molecules data is stored
global PATH_TO_DATA = pwd() * "/data/"

include("Crystal.jl")
include("Forcefield.jl")
include("Molecules.jl")
include("Energetics.jl")
include("Misc.jl")
include("Grid.jl")
#TODO Remove Snapshot.jl from PorousMaterials.jl (This is from Arni's project)
include("Snapshot.jl")
include("GCMC.jl")

export Framework, read_crystal_structure_file, replicate_to_xyz,
       strip_numbers_from_atom_labels!, write_unitcell_boundary_vtk, chemical_formula,
       convert_cif_to_P1_symmetry, # Crystal.jl
       LennardJonesForceField, read_forcefield_file, replication_factors, check_forcefield_coverage, # Forcefield.jl
       Molecule, read_molecule_file, readcharge, # Molecules.jl
       lennard_jones, vdw_energy, # Energetics.jl
       read_xyz, read_cpk_colors, read_atomic_masses, read_atomic_radii, # Misc.jl
	   Grid, write_to_cube, # Grid.jl
	   Snapshot, snap, write_snapshot_to_vtk, write_to_npy # Snapshot.jl #TODO Remove this line later
end
