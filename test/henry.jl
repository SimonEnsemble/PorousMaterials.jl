using PorousMaterials
using Test
using Distributed

@warn "This will take a while..."

insertions_per_volume = 500

@testset "Henry coefficient tests" begin
    ###
    #  Henry test 1: Xe in SBMOF-1
    ###
    crystal = Crystal("SBMOF-1.cif")
    ljff = LJForceField("Dreiding", r_cutoff=12.5)
    molecule = Molecule("Xe")
    temperature = 298.0

    result = henry_coefficient(crystal, molecule, temperature, ljff,
                               insertions_per_volume=insertions_per_volume, verbose=true)
    @test isapprox(result["henry coefficient [mol/(kg-Pa)]"], 0.00985348, rtol=0.025)
    @test isapprox(result["⟨U⟩ (kJ/mol)"], -39.6811, rtol=0.025)

    ###
    #  Henry test 2: CO2 in CAXVII_clean.cif
    ###
    crystal = Crystal("CAXVII_clean.cif")
    ljff = LJForceField("Dreiding", r_cutoff=12.5)
    molecule = Molecule("CO2")
    temperature = 298.0

    result = henry_coefficient(crystal, molecule, temperature, ljff,
                               insertions_per_volume=insertions_per_volume, verbose=true)
    @test isapprox(result["henry coefficient [mol/(kg-Pa)]"], 2.88317e-05, rtol=0.025)
    @test isapprox(result["⟨U⟩ (kJ/mol)"], -18.69582223, rtol=0.025)
    # should not change molecule passed...
    @test isapprox(molecule, Molecule("CO2"))

    ###
    #   Blocked accessible pockets vs. not.
    ###
    crystal = Crystal("LTA.cif")
    molecule = Molecule("CH4")
    ljff = LJForceField("UFF")
    temperature = 298.0
    result = henry_coefficient(crystal, molecule, temperature, ljff,
                               insertions_per_volume=insertions_per_volume, verbose=true)

    accessibility_grid, nb_blocked_pockets = compute_accessibility_grid(crystal, 
        molecule, ljff, n_pts=(50, 50, 50), energy_tol=3.0 * temperature, verbose=true, 
        write_b4_after_grids=false)
    @test nb_blocked_pockets > 0
    
    result = henry_coefficient(crystal, molecule, temperature, ljff,
                               insertions_per_volume=insertions_per_volume, verbose=true,
                               accessibility_grid=accessibility_grid)
    
    crystal = Crystal("LTA_manually_blocked.cif")
    result_manual_block = henry_coefficient(crystal, molecule, temperature, ljff,
                               insertions_per_volume=insertions_per_volume, verbose=true,
                               accessibility_grid=nothing)
    @test isapprox(result["⟨U⟩ (kJ/mol)"], result_manual_block["⟨U⟩ (kJ/mol)"], rtol=0.01)
    @test isapprox(result["henry coefficient [mmol/(g-bar)]"], result_manual_block["henry coefficient [mmol/(g-bar)]"], rtol=0.01)
end
