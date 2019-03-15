module GCMC_Checkpoints_Test

using PorousMaterials
using OffsetArrays
using LinearAlgebra
using Test
using JLD2
using Statistics
using Random

@testset "GCMC Checkpoint Tests" begin
    # idea here: run a simulation with 20 burn, 5 sample shld hv same result as 3 simulations:
    #  (1)                    5 sample, 10 burn, dump checkpoint
    #  (2)   load checkpoint, 5 sample, 15 burn, dump checkpoint
    #  (3)   load checkpoint, 5 sample, 20 burn
    # thus we set random number seed before each simulation to ensure exaclty same moves are
    # conducted.
    framework = Framework("SBMOF-1.cif")
    co2 = Molecule("CO2")
    atom_distances = pairwise_atom_distances(co2, UnitCube())
    charge_distances = pairwise_charge_distances(co2, UnitCube())
    ljff = LJForceField("UFF.csv")
    temp = 298.0
    pressure = 0.5
    molecules = Array{Molecule, 1}

    results = Dict()
    Random.seed!(1234)
    for i = 1:3
        if i == 1
            checkpoint = Dict()
        else
            Core.eval(Main, :(import PorousMaterials)) # used to be able to load in objects from Porous Materials
            # This was found here: https://github.com/JuliaIO/JLD.jl/issues/216
            @load (joinpath(PorousMaterials.PATH_TO_DATA, "gcmc_checkpoints", gcmc_result_savename(framework.name, co2.species, ljff.name, temp, pressure, 5, 5 * i, comment = "_checkpoint", extension=".jld2"))) checkpoint
        end
        results, molecules = gcmc_simulation(framework, deepcopy(co2), temp, pressure, ljff,
                                             n_burn_cycles=5, n_sample_cycles=5 * (i + 1),
                                             verbose=false, sample_frequency=1, eos=:PengRobinson,
                                             autosave=false, write_checkpoints=true,
                                             checkpoint=checkpoint, checkpoint_frequency=1)
        @test all([isapprox(atom_distances, pairwise_atom_distances(molec, UnitCube())) for molec in molecules])
        @test all([isapprox(charge_distances, pairwise_charge_distances(molec, UnitCube())) for molec in molecules])
    end

    Random.seed!(1234)
    results2, molecules2 = gcmc_simulation(framework, deepcopy(co2), temp, pressure, ljff,
                                         n_burn_cycles=5, n_sample_cycles=20,
                                         verbose=false, sample_frequency=1, eos=:PengRobinson,
                                         autosave=false, write_checkpoints=false,
                                         load_checkpoint_file=false, checkpoint_frequency=1)
    @test all([isapprox(atom_distances, pairwise_atom_distances(molec, UnitCube())) for molec in molecules])
    @test all([isapprox(charge_distances, pairwise_charge_distances(molec, UnitCube())) for molec in molecules])
    @test all([isapprox(molecules[i], molecules2[i]) for i = 1:length(molecules)])
    @test length(molecules) == length(molecules2)
    @test isapprox(results["Q_st (K)"], results2["Q_st (K)"])
    @test isapprox(results["⟨N⟩ (mmol/g)"], results2["⟨N⟩ (mmol/g)"])
end
end
