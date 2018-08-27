module GCMC_Checkpoints_Test

using PorousMaterials
using OffsetArrays
using LinearAlgebra
using Test
using JLD2
using Statistics
using Random

@testset "GCMC Checkpoint Tests" begin
    using PorousMaterials
    # idea here: run a simulation with 20 burn, 5 sample shld hv same result as 3 simulations:
    #  (1)                    5 sample, 10 burn, dump checkpoint
    #  (2)   load checkpoint, 5 sample, 15 burn, dump checkpoint
    #  (3)   load checkpoint, 5 sample, 20 burn
    # thus we set random number seed before each simulation to ensure exaclty same moves are
    # conducted.
    framework = Framework("SBMOF-1.cif")
    co2 = Molecule("CO2")
    co_bond_length = norm(co2.atoms[1].xf - co2.atoms[2].xf)
    ljff = LJForceField("UFF.csv")
    temp = 298.0
    pressure = 0.5
    molecules = Array{Molecule, 1}

    results = Dict()
    Random.seed!(1234)
    for i = 1:3
        using PorousMaterials
        if i == 1
            checkpoint = Dict()
        else
            using PorousMaterials
            println("LOADING...")
            @load (PorousMaterials.PATH_TO_DATA * "/gcmc_checkpoints/" * gcmc_result_savename(framework.name, co2.species, ljff.name, temp, pressure, 5, 5 * i, comment = "_checkpoint")) checkpoint
        end
        results, molecules = gcmc_simulation(framework, deepcopy(co2), temp, pressure, ljff,
                                             n_burn_cycles=5, n_sample_cycles=5 * (i + 1),
                                             verbose=true, sample_frequency=1, eos=:PengRobinson,
                                             autosave=false, write_checkpoints=true,
                                             checkpoint=checkpoint, checkpoint_frequency=1)
        @test isapprox(norm(molecules[1].atoms[1].xf - molecules[1].atoms[2].xf), co_bond_length)
    end

    Random.seed!(1234)
    results2, molecules2 = gcmc_simulation(framework, deepcopy(co2), temp, pressure, ljff,
                                         n_burn_cycles=5, n_sample_cycles=20,
                                         verbose=true, sample_frequency=1, eos=:PengRobinson,
                                         autosave=false, write_checkpoints=false,
                                         load_checkpoint_file=false, checkpoint_frequency=1)
    @test isapprox(norm(molecules2[1].atoms[1].xf - molecules2[1].atoms[2].xf), co_bond_length)
    @test length(molecules) == length(molecules2)
    @test all(isapprox.(molecules, molecules2))
    @test isapprox(results["Q_st (K)"], results2["Q_st (K)"])
    @test isapprox(results["⟨N⟩ (mmol/g)"], results2["⟨N⟩ (mmol/g)"])
end
end
