module PorousMaterials

using CSV
using DataFrames

# this is the directory where crystal structures, forcefields, and molecules data is stored
global PATH_TO_DATA = pwd() * "/data/"

include("Box.jl")
include("Matter.jl")
include("NearestImage.jl")
include("Misc.jl")
include("Crystal.jl")
include("Molecules.jl")
include("Forcefield.jl")
 # include("Energetics_Util.jl")
 # include("VdWEnergetics.jl")
 # include("ElectrostaticEnergetics.jl")
 # include("Misc.jl")
 # include("Grid.jl")
 # include("MChelpers.jl")
 # include("EOS.jl")
 # include("Henry.jl")
 # include("GCMC.jl")
 # 
export Box, replicate, # Box.jl
       LJSphere, PointCharge, # Matter.jl 
        Framework, read_crystal_structure_file, replicate_to_xyz, remove_overlapping_atoms,
       strip_numbers_from_atom_labels!, write_unitcell_boundary_vtk, chemical_formula, molecular_weight, crystal_density,
       convert_cif_to_P1_symmetry, construct_box, replicate, read_atomic_masses, charged, write_cif, # Crystal.jl
       Molecule, PointCharge, LJSphere, translate_by!,
       translate_to!, rotate!, rotation_matrix, rand_point_on_unit_sphere, charged, # Molecules.jl
       LennardJonesForceField, ForceField, replication_factors, check_forcefield_coverage # Forcefield.jl
 #        outside_box, write_to_xyz,
 #        nearest_image!, PotentialEnergy, SystemPotentialEnergy, # Energetics_Util.jl
 #        lennard_jones, vdw_energy, vdw_energy_no_PBC, # VdWEnergetics.jl
 #        read_xyz, read_cpk_colors, read_atomic_radii, # Misc.jl
 #        Grid, write_cube, read_cube, energy_grid, # Grid.jl
 #        insert_molecule!, delete_molecule!, translate_molecule!, reinsert_molecule!, rotatable, # MChelpers.jl
 #        apply_periodic_boundary_condition!,
 #        gcmc_simulation, adsorption_isotherm, stepwise_adsorption_isotherm,  # GCMC.jl
 #        electrostatic_potential, electrostatic_potential_energy, precompute_kvec_wts, setup_Ewald_sum, total, # ElectrostaticEnergetics.jl
 #        calculate_properties, PengRobinsonGas, # EOS.jl
 #        henry_coefficient # Henry.jl
end
