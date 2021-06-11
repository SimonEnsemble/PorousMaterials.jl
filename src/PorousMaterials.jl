module PorousMaterials

using Roots, OffsetArrays, SpecialFunctions, StatsBase, ProgressMeter, Polynomials,
JLD2, Statistics, Distributed, Optim, Printf, DataFrames, LightGraphs, CSV, LinearAlgebra

# extend Xtals
using Reexport
@reexport using Xtals
import Xtals.Cart, Xtals.Frac, Xtals.write_xyz

# physical constants
const univ_gas_const = 8.3144598e-5 # m³-bar/(K-mol)
const K_to_kJ_per_mol = 8.3144598e-3 # kJ/(mol-K)
const boltzmann = 1.38064852e7 # Boltmann constant (Pa-m3/K --> Pa-A3/K)

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
include("atomic_masses.jl")

function __init__()
    rc[:paths][:forcefields] = joinpath(rc[:paths][:data], "forcefields")
    rc[:paths][:molecules] = joinpath(rc[:paths][:data], "molecules")
    rc[:paths][:grids] = joinpath(rc[:paths][:data], "grids")
    rc[:paths][:sims] = joinpath(rc[:paths][:data], "simulations")
    append_atomic_masses()
    @debug "Environment variables" rc
end


export
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
    required_n_pts, xf_to_id, id_to_xf, update_density!, find_energy_minimum, origin,
 
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
