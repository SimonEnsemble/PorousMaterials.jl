module PorousMaterials

using CSV
using DataFrames
using Roots # for fzero
using OffsetArrays # used for Ewald sum
using SpecialFunctions # for erfc
using StatsBase
using ProgressMeter
using Polynomials
using JLD2
using Statistics
using Printf
using LinearAlgebra
using LightGraphs
using Distributed
using Optim
import Base.push!


"""
    print_file_paths()

Print off paths where PorousMaterials.jl looks for input files and writes output files.
"""
function print_file_paths()
    println("General data folder: ", PATH_TO_DATA)
    println("\tcrystal structures (.cif, .cssr): ", PATH_TO_CRYSTALS)
    println("\tforce field files (.csv): ", PATH_TO_FORCEFIELDS)
    println("\tmolecule input files: ", PATH_TO_MOLECULES)
    println("\tgrids (.cube): ", PATH_TO_GRIDS)
end

"""
    set_default_file_paths(print_paths=true)

Sets the default paths for where input files and some output files are stored.
To see current set up, call [`print_file_paths`](@ref)
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

# this runs everytime PorousMaterials is loaded, so if the user changes directory
#   then the PATH_TO_DATA will change as well
function __init__()
    set_default_file_paths(print_paths=false)
end

"""
    set_tutorial_mode()

Places PorousMaterials in "Tutorial Mode". It changes the `PATH_TO_DATA` variable to
the directory where the PorousMaterials test data is stored. It can be used to
follow examples shown in the README. It displays a warning so that the user knows
They are no longer using their own data.
"""
function set_tutorial_mode()
    new_path = joinpath(dirname(pathof(PorousMaterials)), "..", "test", "data")
    if ! isdir(new_path)
        @error @sprintf("Directory for testing data %s does not exist.\nNot entering Tutorial Mode.\n", new_path)
    else
        global PATH_TO_DATA = new_path
        global PATH_TO_CRYSTALS = joinpath(PATH_TO_DATA, "crystals")
        global PATH_TO_FORCEFIELDS = joinpath(PATH_TO_DATA, "forcefields")
        global PATH_TO_MOLECULES = joinpath(PATH_TO_DATA, "molecules")
        global PATH_TO_GRIDS = joinpath(PATH_TO_DATA, "grids")
        @warn "PorousMaterials is now in Tutorial Mode. You have access to the testing data to experiment with PorousMaterials.\nTo reset to default file paths, call `set_default_file_paths()`\n"
    end
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
include("generic_rotations.jl")

export
    # PorousMaterials.jl
    set_default_file_paths, print_file_paths, set_tutorial_mode,

    # Box.jl
    Box, replicate, UnitCube, write_vtk, inside,

    # Matter.jl
    Atoms, Charges,

    # NearestImage.jl
    nearest_image!, nearest_r², nearest_r,

    # Misc.jl
    read_xyz, read_cpk_colors, read_atomic_radii, write_xyz, fit_adsorption_isotherm,

    # Crystal.jl
    Framework, read_crystal_structure_file, remove_overlapping_atoms_and_charges,
    strip_numbers_from_atom_labels!, chemical_formula, molecular_weight, crystal_density,
    construct_box, replicate, read_atomic_masses, charged, write_cif, assign_charges,

    # Molecules.jl
    Molecule, n_atoms, set_fractional_coords!, translate_by!, outside_box,
    set_fractional_coords_to_unit_cube!, translate_to!, rotate!, rotation_matrix,
    rand_point_on_unit_sphere, charged, pairwise_atom_distances,
    pairwise_charge_distances, Ion, bond_length_drift,

    # Forcefield.jl
    LJForceField, replication_factors, check_forcefield_coverage,

    # Energetics_Util.jl
    PotentialEnergy, SystemPotentialEnergy,

    # VdWEnergetics.jl
    lennard_jones, vdw_energy, vdw_energy_no_PBC,

    # ElectrostaticEnergetics.jl
    electrostatic_potential, electrostatic_potential_energy, precompute_kvec_wts,
    setup_Ewald_sum, total, Eikr, total_electrostatic_potential_energy,

    # MChelpers.jl
    insert_molecule!, delete_molecule!, translate_molecule!, reinsert_molecule!, rotatable,

    # Grid.jl
    apply_periodic_boundary_condition!,
    Grid, write_cube, read_cube, energy_grid, compute_accessibility_grid, accessible,
    required_n_pts, xf_to_id, id_to_xf, update_density!,

    # EOS.jl
    calculate_properties, PengRobinsonFluid, VdWFluid,

    # GCMC.jl
    gcmc_simulation, adsorption_isotherm, stepwise_adsorption_isotherm,
    gcmc_result_savename, GCMCstats, MarkovCounts,

    # Henry.jl
    henry_coefficient, henry_result_savename,

    # generic_rotations.jl
    rotation_matrix
end
