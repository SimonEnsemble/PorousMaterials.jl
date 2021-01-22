using Distributed
using PorousMaterials
using Test
using CSV
using DataFrames
using JLD2
using Printf
using Statistics

# put quick tests here that can run on Travis in a reasonable amount of time.

@testset "GCMC (quick) tests" begin
    ###
    #   isotherm_sim_results_to_dataframe test
    ###
    desired_props = ["xtal", "adsorbate", "pressure (bar)", 
                     "repfactors", "⟨N⟩ (mmol/g)", "density grid"]

    xtal = Crystal("SBMOF-1.cif")
    molecule = Molecule("Xe") # adsorbate
    ljff = LJForceField("UFF", mixing_rules="Lorentz-Berthelot")
    temperature= 298.0 # temperature: K
    n_sample_cycles = 5000
    n_burn_cycles = 5000

    pressures = 10 .^ range(-2.0, stop=0.0, length=3)

    ### manually checked test data ###
    test_mmol_g = [1.088001225436118, 1.3978990363878692, 1.4456497503632355]

    # read output files into dataframe
    df_data =  isotherm_sim_results_to_dataframe(desired_props, xtal, 
                                                 molecule, temperature, pressures, 
                                                 ljff, n_burn_cycles, n_sample_cycles)

    # check the dataframe entries against manualy checked output
    @test isapprox(df_data[:, Symbol("pressure (bar)")],  pressures)
    @test isapprox(df_data[:, Symbol("⟨N⟩ (mmol/g)")], test_mmol_g)
    @test all(df_data[:, :xtal] .== xtal.name)
    @test all(df_data[:, :adsorbate] .== molecule.species)

    # check the data type of specific columns
    @test all(typeof.(df_data[:, :repfactors]) .== Tuple{Int64,Int64,Int64})
    @test all(typeof.(df_data[:, Symbol("density grid")]) .== Grid{Float64})
end

@warn "GCMC simulations not included in runtests.jl; run gcmc_long.jl to test GCMC simulations (will take hours)"
