module PorousMaterials

using CSV
using DataFrames
 # using Roots # for fzero
 # using OffsetArrays # used for Ewald sum
 # using SpecialFunctions # for erfc
 # using StatsBase
 # using ProgressMeter
 # using Polynomials
 # using JLD2
 # using Statistics
using Printf
using LinearAlgebra
using LightGraphs
 # using Distributed
using Optim
 # import Base.push!
 # 

# atoms are considered to overlap if this close.
const R²_OVERLAP = 0.1 # Units: Angstrom²

"""
    print_file_paths()

print off paths where PorousMaterials.jl looks for input files and writes output files.
"""
function print_file_paths()
    println("general data folder: ", PATH_TO_DATA)
    println("\tcrystal structures (.cif, .cssr): ", PATH_TO_CRYSTALS)
    println("\tforce field files (.csv): ", PATH_TO_FORCEFIELDS)
    println("\tmolecule input files: ", PATH_TO_MOLECULES)
    println("\tgrids (.cube): ", PATH_TO_GRIDS)
end

"""
    set_default_file_paths(print_paths=true)

sets the default paths for where input files and some output files are stored.
to see current set up, call [`print_file_paths`](@ref)
"""
function set_default_file_paths(;print_paths::Bool=true)
    # this is the main directory where crystal structures, forcefields, and molecules data is stored
    global PATH_TO_DATA = joinpath(pwd(), "data")

    global PATH_TO_CRYSTALS = joinpath(PATH_TO_DATA, "crystals")
    global PATH_TO_FORCEFIELDS = joinpath(PATH_TO_DATA, "forcefields")
    global PATH_TO_MOLECULES = joinpath(PATH_TO_DATA, "molecules")
    global PATH_TO_GRIDS = joinpath(PATH_TO_DATA, "grids")
    
    if print_paths
        print_file_paths()
    end
end

# this runs everytime porousmaterials is loaded, so if the user changes directory
#   then the path_to_data will change as well
function __init__()
    set_default_file_paths(print_paths=false)
end

 # """
 #     set_tutorial_mode()
 # 
 # places porousmaterials in "tutorial mode". it changes the `path_to_data` variable to
 # the directory where the porousmaterials test data is stored. it can be used to
 # follow examples shown in the readme. it displays a warning so that the user knows
 # they are no longer using their own data.
 # """
 # function set_tutorial_mode()
 #     new_path = joinpath(dirname(pathof(porousmaterials)), "..", "test", "data")
 #     if ! isdir(new_path)
 #         @error @sprintf("directory for testing data %s does not exist.\nnot entering tutorial mode.\n", new_path)
 #     else
 #         global path_to_data = new_path
 #         global path_to_crystals = joinpath(path_to_data, "crystals")
 #         global path_to_forcefields = joinpath(path_to_data, "forcefields")
 #         global path_to_molecules = joinpath(path_to_data, "molecules")
 #         global path_to_grids = joinpath(path_to_data, "grids")
 #         @warn "porousmaterials is now in tutorial mode. you have access to the testing data to experiment with porousmaterials.\nto reset to default file paths, call `set_default_file_paths()`\n"
 #     end
 # end
 # 
include("matter.jl")
include("box.jl")
include("distance.jl")
include("misc.jl")
include("isotherm_fitting.jl")
include("crystal.jl")
include("bonds.jl")
include("forcefield.jl")
include("molecule.jl")
 # include("energetics_util.jl")
include("vdwenergetics.jl")
 # include("electrostaticenergetics.jl")
 # include("mchelpers.jl")
 # include("grid.jl")
 # include("eos.jl")
 # include("henry.jl")
 # include("gcmc.jl")
 # include("generic_rotations.jl")

export
 #     # porousmaterials.jl
 #     set_default_file_paths, print_file_paths, set_tutorial_mode,
 # 
    # matter.jl
    Coords, Frac, Cart, Atoms, Charges, wrap!, neutral, net_charge, translate_by!,
    
    # box.jl
    Box, replicate, unit_cube, write_vtk, inside, fractional_coords, cartesian_coords,

    # distance.jl
    nearest_image!, distance, overlap, remove_duplicates,

    # misc.jl
    read_xyz, read_cpk_colors, write_xyz, read_atomic_masses,
    
    # isotherm_fitting.jl
    fit_adsorption_isotherm,

    # crystal.jl
    Crystal, strip_numbers_from_atom_labels!, assign_charges,
    chemical_formula, molecular_weight, crystal_density, write_cif, has_charges,
    apply_symmetry_operations, 

    # bonds.jl
    infer_bonds!, write_bond_information, BondingRule, bond_sanity_check, remove_bonds!,

 #     construct_box,
 #     replicate, read_atomic_masses, charged, write_cif, assign_charges,
 #     is_symmetry_equal, apply_symmetry_rules, assert_p1_symmetry, infer_bonds!,
 #     remove_bonds!, compare_bonds_in_framework, wrap_atoms_to_unit_cell!,
 #     write_bond_information, is_bonded, default_bondingrules, has_same_sets_of_atoms_and_charges, 
 #     distance, bond_sanity_check,
 #
 #     # FrameworkOperations.jl
 #     partition_framework, subtract_atoms,
 # 
    # molecule.jl
    Molecule, translate_by!, translate_to!, rotate!, rotation_matrix,
    rand_point_on_unit_sphere, ion,
  
    # Forcefield.jl
    LJForceField, replication_factors, check_forcefield_coverage,
 # 
 #     # Energetics_Util.jl
 #     PotentialEnergy, SystemPotentialEnergy,
 # 
      # VdWEnergetics.jl
      lennard_jones, vdw_energy#, vdw_energy_no_PBC
 # 
 #     # ElectrostaticEnergetics.jl
 #     electrostatic_potential, electrostatic_potential_energy, precompute_kvec_wts,
 #     setup_Ewald_sum, total, Eikr, total_electrostatic_potential_energy,
 # 
 #     # MChelpers.jl
 #     insert_molecule!, delete_molecule!, translate_molecule!, reinsert_molecule!, rotatable,
 # 
 #     # Grid.jl
 #     apply_periodic_boundary_condition!,
 #     Grid, write_cube, read_cube, energy_grid, compute_accessibility_grid, accessible,
 #     required_n_pts, xf_to_id, id_to_xf, update_density!,
 # 
 #     # EOS.jl
 #     calculate_properties, PengRobinsonFluid, VdWFluid,
 # 
 #     # GCMC.jl
 #     gcmc_simulation, adsorption_isotherm, stepwise_adsorption_isotherm,
 #     gcmc_result_savename, GCMCstats, MarkovCounts,
 # 
 #     # Henry.jl
 #     henry_coefficient, henry_result_savename,
 # 
 #     # generic_rotations.jl
 #     rotation_matrix
end
