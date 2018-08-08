__precompile__()
module PorousMaterials

using CSV
using DataFrames
using Roots # for fzero
using OffsetArrays # used for Ewald sum
using SpecialFunctions # for erfc
using StatsBase
using ProgressMeter
using JLD
using Polynomials

# this runs everytime PorousMaterials is loaded, so if the user changes directory
#   then the PATH_TO_DATA will change as well
function __init__()
    # this is the directory where crystal structures, forcefields, and molecules data is stored
    global PATH_TO_DATA = pwd() * "/data/"
end

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

export Box, replicate, UnitCube, write_vtk, # Box.jl
       LJSphere, PtCharge, # Matter.jl
       nearest_image!, nearest_r², nearest_r, # NearestImage.jl
       read_xyz, read_cpk_colors, read_atomic_radii, write_to_xyz, # Misc.jl
       Framework, read_crystal_structure_file, remove_overlapping_atoms,
<<<<<<< HEAD
       strip_numbers_from_atom_labels!, write_unitcell_boundary_vtk, chemical_formula, molecular_weight, crystal_density,
       convert_cif_to_P1_symmetry, construct_box, replicate, read_atomic_masses, charged, write_cif, assign_charges_to_zeolites!, # Crystal.jl
       LennardJonesForceField, read_forcefield_file, replication_factors, check_forcefield_coverage, # Forcefield.jl
       Molecule, PointCharge, LennardJonesSphere, read_molecule_file, translate_by!,
=======
       strip_numbers_from_atom_labels!, chemical_formula, molecular_weight, crystal_density,
       construct_box, replicate, read_atomic_masses, charged, write_cif, assign_charges, # Crystal.jl
       Molecule, set_fractional_coords!, translate_by!, outside_box, set_fractional_coords_to_unit_cube!,
>>>>>>> bde520c00e87fdfdf7777166b08c7aeb50d8ad96
       translate_to!, rotate!, rotation_matrix, rand_point_on_unit_sphere, charged, # Molecules.jl
       LJForceField, replication_factors, check_forcefield_coverage, # Forcefield.jl
       PotentialEnergy, SystemPotentialEnergy, # Energetics_Util.jl
       lennard_jones, vdw_energy, vdw_energy_no_PBC, # VdWEnergetics.jl
       electrostatic_potential, electrostatic_potential_energy, precompute_kvec_wts, setup_Ewald_sum, total, # ElectrostaticEnergetics.jl
       insert_molecule!, delete_molecule!, translate_molecule!, reinsert_molecule!, rotatable, # MChelpers.jl
       apply_periodic_boundary_condition!,
       gcmc_simulation, adsorption_isotherm, stepwise_adsorption_isotherm, root_save_filename,  # GCMC.jl
       electrostatic_potential, electrostatic_potential_energy, precompute_kvec_wts, setup_Ewald_sum, total, # ElectrostaticEnergetics.jl
       calculate_properties, PengRobinsonGas, # EOS.jl
       gcmc_simulation, adsorption_isotherm, stepwise_adsorption_isotherm,
       gcmc_result_savename, # GCMC.jl
       henry_coefficient, henry_result_savename # Henry.jl
end
