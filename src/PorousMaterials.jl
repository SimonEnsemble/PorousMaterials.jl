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
include("Energetics_Util.jl")
include("VdWEnergetics.jl")
include("ElectrostaticEnergetics.jl")
include("MChelpers.jl")
include("Grid.jl")
include("EOS.jl")
include("Henry.jl")
include("GCMC.jl")

export Box, replicate, UnitCube, # Box.jl
       LJSphere, PtCharge, # Matter.jl 
       nearest_image!, nearest_r², nearest_r, # NearestImage.jl
       read_xyz, read_cpk_colors, read_atomic_radii, write_to_xyz, # Misc.jl
       Framework, read_crystal_structure_file, replicate_to_xyz, remove_overlapping_atoms,
       strip_numbers_from_atom_labels!, write_unitcell_boundary_vtk, chemical_formula, molecular_weight, crystal_density,
       construct_box, replicate, read_atomic_masses, charged, write_cif, # Crystal.jl
       Molecule, set_fractional_coords!, translate_by!, outside_box, set_fractional_coords_to_unit_cube!,
       translate_to!, rotate!, rotation_matrix, rand_point_on_unit_sphere, charged, # Molecules.jl
       LJForceField, replication_factors, check_forcefield_coverage, # Forcefield.jl
       PotentialEnergy, SystemPotentialEnergy, # Energetics_Util.jl
       lennard_jones, vdw_energy, vdw_energy_no_PBC, # VdWEnergetics.jl
       electrostatic_potential, electrostatic_potential_energy, precompute_kvec_wts, setup_Ewald_sum, total, # ElectrostaticEnergetics.jl
       insert_molecule!, delete_molecule!, translate_molecule!, reinsert_molecule!, rotatable, # MChelpers.jl
       apply_periodic_boundary_condition!,
       Grid, write_cube, read_cube, energy_grid, # Grid.jl
       calculate_properties, PengRobinsonGas, # EOS.jl
       gcmc_simulation, adsorption_isotherm, stepwise_adsorption_isotherm,  # GCMC.jl
       henry_coefficient # Henry.jl
end
