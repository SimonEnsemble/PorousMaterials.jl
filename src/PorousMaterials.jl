module PorousMaterials

# this is the directory where crystal structures, forcefields, and molecules data is stored
global PATH_TO_DATA = pwd() * "/data/"

include("Crystal.jl")
include("Forcefield.jl")
include("Molecules.jl")
include("Energetics_Util.jl")
include("Energetics.jl")
include("Misc.jl")
include("Grid.jl")
include("GCMC.jl")

export Box, Framework, read_crystal_structure_file, replicate_to_xyz,
       strip_numbers_from_atom_labels!, write_unitcell_boundary_vtk, chemical_formula, molecular_weight, crystal_density,
       convert_cif_to_P1_symmetry, construct_box, replicate_box, # Crystal.jl
       LennardJonesForceField, read_forcefield_file, replication_factors, check_forcefield_coverage, # Forcefield.jl
       Molecule, read_molecule_file, readcharge, # Molecules.jl
       lennard_jones, vdw_energy, # Energetics.jl
       completely_outside_box,
       read_xyz, read_cpk_colors, read_atomic_masses, read_atomic_radii, # Misc.jl
	Snapshot, snap, write_snapshot_to_vtk, write_to_npy, write_snapshot_to_cube, # Snapshot.jl #TODO Remove this line later
	GCMCstats, MarkovCounts, insert_molecule!, delete_molecule!, translate_molecule!,
	guest_guest_vdw_energy, gcmc_simulation, bring_molecule_inside_box! # GCMC.jl
end
