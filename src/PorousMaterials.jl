module PorousMaterials

# this is the directory where crystal structures, forcefields, and molecules data is stored
global PATH_TO_DATA = pwd() * "/data/"

include("Crystal.jl")
include("Molecules.jl")
include("Forcefield.jl")
include("Energetics_Util.jl")
include("VdWEnergetics.jl")
include("ElectrostaticEnergetics.jl")
include("Misc.jl")
include("Grid.jl")
include("MChelpers.jl")
include("GCMC.jl")

export Box, Framework, read_crystal_structure_file, replicate_to_xyz, remove_overlapping_atoms,
       strip_numbers_from_atom_labels!, write_unitcell_boundary_vtk, chemical_formula, molecular_weight, crystal_density,
       convert_cif_to_P1_symmetry, construct_box, replicate_box, read_atomic_masses, charged, # Crystal.jl
       LennardJonesForceField, read_forcefield_file, replication_factors, check_forcefield_coverage, # Forcefield.jl
       Molecule, read_molecule_file, translate_by!, translate_to!, rotate!, rotation_matrix, rand_point_on_unit_sphere, charged, # Molecules.jl
       outside_box, write_to_xyz,
       nearest_image!, V_vdw, V_electro, PotentialEnergy, guest_guest, guest_host, # Energetics_Util.jl
       lennard_jones, vdw_energy, # VdWEnergetics.jl
       read_xyz, read_cpk_colors, read_atomic_radii, # Misc.jl
       Grid, write_cube, read_cube, energy_grid, # Grid.jl
       insert_molecule!, delete_molecule!, translate_molecule!, # MChelpers.jl
       apply_periodic_boundary_condition!,
       guest_guest_vdw_energy, gcmc_simulation,  # GCMC.jl
       electrostatic_potential, electrostatic_potential_energy, precompute_kvec_wts, allocate_eikr, setup_Ewald_sum # ElectrostaticEnergetics.jl
end
