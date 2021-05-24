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
using PyCall
using Reexport
@reexport using Xtals


import Xtals.print_file_paths
"""
    print_file_paths()

print off paths where PorousMaterials.jl looks for input files and writes output files.
"""
function print_file_paths()
    println("general data folder: ", PATH_TO_DATA)
    println("\tcrystal structures (.cif, .cssr): ", Xtals.PATH_TO_CRYSTALS)
    println("\tforce field files (.csv): ", PATH_TO_FORCEFIELDS)
    println("\tmolecule input files: ", PATH_TO_MOLECULES)
    println("\tsimulation output files: ", PATH_TO_SIMS)
    println("\tgrids (.cube): ", PATH_TO_GRIDS)
end


import Xtals.set_path_to_data
"""
    set_path_to_data(path; print_paths=true)

Set the path variables: `PATH_TO_DATA`, `PATH_TO_CRYSTALS`, `PATH_TO_FORCEFIELDS`, `PATH_TO_MOLECULES`,
`PATH_TO_GRIDS`, and `PATH_TO_SIMS`.  The latter five paths are set relative to the root data path.

# Arguments
- `path::String`: the absolute path to the root of the data directory.
- `print_paths::Bool`: set false to suppress printing of path values.
"""
function set_path_to_data(path::String; print_paths::Bool=true)
    global PATH_TO_DATA = path
    set_path_to_crystals(joinpath(PATH_TO_DATA, "crystals"))
    global PATH_TO_FORCEFIELDS = joinpath(PATH_TO_DATA, "forcefields")
    global PATH_TO_MOLECULES = joinpath(PATH_TO_DATA, "molecules")
    global PATH_TO_GRIDS = joinpath(PATH_TO_DATA, "grids")
    global PATH_TO_SIMS = joinpath(PATH_TO_DATA, "simulations")

    if print_paths
        print_file_paths()
    end
end


# this runs everytime porousmaterials is loaded, so if the user changes directory
#   then the path_to_data will change as well
function __init__()
    global PATH_TO_DATA = joinpath(pwd(), "data")
    set_path_to_crystals(joinpath(PATH_TO_DATA, "crystals"))
    global PATH_TO_FORCEFIELDS = joinpath(PATH_TO_DATA, "forcefields")
    global PATH_TO_MOLECULES = joinpath(PATH_TO_DATA, "molecules")
    global PATH_TO_GRIDS = joinpath(PATH_TO_DATA, "grids")
    global PATH_TO_SIMS = joinpath(PATH_TO_DATA, "simulations")
end


include("misc.jl")
include("isotherm_fitting.jl")
include("forcefield.jl")
include("molecule.jl")
include("energy_utilities.jl")
include("vdw_energetics.jl")
include("electrostatics.jl")
include("mc_helpers.jl")
include("grid.jl")
include("eos.jl")
include("henry.jl")
include("gcmc.jl")
include("energy_min.jl")


export
    # porousmaterials.jl
    print_file_paths, set_path_to_data,

    # isotherm_fitting.jl
    fit_adsorption_isotherm,

    # molecule.jl
    Molecule, translate_by!, translate_to!, random_rotation!, random_rotation_matrix, ion, distortion,

    # forcefield.jl
    LJForceField, replication_factors, forcefield_coverage,

    # energy_utilities.jl
    PotentialEnergy, SystemPotentialEnergy,

    # vdw_energetics.jl
    lennard_jones, vdw_energy, vdw_energy_no_PBC,

    # electrostatics.jl
    electrostatic_potential, electrostatic_potential_energy, precompute_kvec_wts,
    setup_Ewald_sum, total, Eikr, total_electrostatic_potential_energy,

    # mc_helpers.jl
    random_insertion!, remove_molecule!, random_translation!, random_reinsertion!, needs_rotations,

    # Grid.jl
    apply_periodic_boundary_condition!,
    Grid, write_cube, read_cube, energy_grid, compute_accessibility_grid, accessible,
    required_n_pts, xf_to_id, id_to_xf, update_density!, find_energy_minimum,
 
    # EOS.jl
    calculate_properties, PengRobinsonFluid, VdWFluid,

    # gcmc.jl
    μVT_sim, adsorption_isotherm, stepwise_adsorption_isotherm,
    μVT_output_filename, GCMCstats, MarkovCounts, isotherm_sim_results_to_dataframe,
 
    # henry.jl
    henry_coefficient, henry_result_savename,

    # energy_min.jl
    find_energy_minimum, find_energy_minimum_gridsearch
end
