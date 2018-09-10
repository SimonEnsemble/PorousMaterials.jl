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
import Base.push!

global const PATH_TO_PACKAGE_DATA = dirname(pathof(PorousMaterials)) * "/../data/"

# this runs everytime PorousMaterials is loaded, so if the user changes directory
#   then the PATH_TO_DATA will change as well
function __init__()
    # every time PorousMaterials is started, the PATH_TO_DATA defaults to using the
    #   packages data, but can be changed with the set_path_to_data and reset_path_to_data function
    global PATH_TO_DATA = PATH_TO_PACKAGE_DATA
    global USING_PACKAGE_DATA = true
end

function reset_path_to_data()
    global PATH_TO_DATA
    global USING_PACKAGE_DATA
    USING_PACKAGE_DATA = true
    PATH_TO_DATA = PATH_TO_PACKAGE_DATA
    @printf("Current PATH_TO_DATA: %s\n", PATH_TO_DATA)
end

function set_path_to_data()
    global PATH_TO_DATA
    global USING_PACKAGE_DATA
    USING_PACKAGE_DATA = false
    PATH_TO_DATA = pwd() * "/"
    @printf("Current PATH_TO_DATA: %s\n", PATH_TO_DATA)
end

function set_path_to_data(new_path::AbstractString)
    global PATH_TO_DATA
    global USING_PACKAGE_DATA
    if isdir(new_path)
        if new_path[end] != "/"
            new_path = new_path * "/"
        end
        USING_PACKAGE_DATA = false
        PATH_TO_DATA = new_path
    else
        error(@sprintf("The path %s does not exist", new_path))
    end
    @printf("Current PATH_TO_DATA: %s\n", PATH_TO_DATA)
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

export
    # PorousMaterials.jl
    reset_path_to_data, set_path_to_data,

    # Box.jl
    Box, replicate, UnitCube, write_vtk,

    # Matter.jl
    Atoms, Charges,

    # NearestImage.jl
    nearest_image!, nearest_r², nearest_r,

    # Misc.jl
    read_xyz, read_cpk_colors, read_atomic_radii, write_xyz,

    # Crystal.jl
    Framework, read_crystal_structure_file, remove_overlapping_atoms,
    strip_numbers_from_atom_labels!, chemical_formula, molecular_weight, crystal_density,
    construct_box, replicate, read_atomic_masses, charged, write_cif, assign_charges,

    # Molecules.jl
    Molecule, set_fractional_coords!, translate_by!, outside_box, set_fractional_coords_to_unit_cube!,
    translate_to!, rotate!, rotation_matrix, rand_point_on_unit_sphere, charged,
    pairwise_atom_distances, pairwise_charge_distances, Ion,

    # Forcefield.jl
    LJForceField, replication_factors, check_forcefield_coverage,

    # Energetics_Util.jl
    PotentialEnergy, SystemPotentialEnergy,

    # VdWEnergetics.jl
    lennard_jones, vdw_energy, vdw_energy_no_PBC,

    # ElectrostaticEnergetics.jl
    electrostatic_potential, electrostatic_potential_energy, precompute_kvec_wts, setup_Ewald_sum, total, Eikr,

    # MChelpers.jl
    insert_molecule!, delete_molecule!, translate_molecule!, reinsert_molecule!, rotatable,

    # Grid.jl
    apply_periodic_boundary_condition!,
    Grid, write_cube, read_cube, energy_grid,

    # EOS.jl
    calculate_properties, PengRobinsonGas,

    # GCMC.jl
    gcmc_simulation, adsorption_isotherm, stepwise_adsorption_isotherm,
    gcmc_result_savename, GCMCstats, MarkovCounts,

    # Henry.jl
    henry_coefficient, henry_result_savename
end
