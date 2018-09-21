using PorousMaterials
using Test
using Distributed

@warn "This will take a while..."

insertions_per_volume = 1000

@testset "Henry coefficient tests" begin
    ###
    #  Henry test 1: Xe in SBMOF-1
    ###
    framework = Framework("SBMOF-1.cif")
    ljff = LJForceField("Dreiding.csv", cutoffradius=12.5)
    molecule = Molecule("Xe")
    temperature = 298.0

    result = henry_coefficient(framework, molecule, temperature, ljff,
                               insertions_per_volume=insertions_per_volume, verbose=true)
    @test isapprox(result["henry coefficient [mol/(kg-Pa)]"], 0.00985348, atol=0.0002)
    @test isapprox(result["⟨U⟩ (kJ/mol)"], -39.6811, atol=0.1)

    ###
    #  Henry test 2: CO2 in CAXVII_clean.cif
    ###
    framework = Framework("CAXVII_clean.cif")
    ljff = LJForceField("Dreiding.csv", cutoffradius=12.5)
    molecule = Molecule("CO2")
    temperature = 298.0

    result = henry_coefficient(framework, molecule, temperature, ljff,
                               insertions_per_volume=insertions_per_volume, verbose=true)
    @test isapprox(result["henry coefficient [mol/(kg-Pa)]"], 2.88317e-05, atol=1.5e-7)
    @test isapprox(result["⟨U⟩ (kJ/mol)"], -18.69582223, atol=0.1)
    # should not change molecule passed...
    @test isapprox(molecule, Molecule("CO2"))

    ###
    #   Blocked accessible pockets vs. not.
    ###
    framework = Framework("LTA.cif")
    molecule = Molecule("CH4")
    ljff = LJForceField("UFF.csv")
    temperature = 298.0
    result = henry_coefficient(framework, molecule, temperature, ljff,
                               insertions_per_volume=insertions_per_volume, verbose=true)

    accessibility_grid, some_pockets_were_blocked = compute_accessibility_grid(framework, 
        molecule, ljff, n_pts=(50, 50, 50), energy_tol=3.0 * temperature, verbose=true, 
        write_b4_after_grids=false)
    @test some_pockets_were_blocked
    
    result = henry_coefficient(framework, molecule, temperature, ljff,
                               insertions_per_volume=insertions_per_volume, verbose=true,
                               accessibility_grid=accessibility_grid)
    
    framework = Framework("LTA_manually_blocked.cif")
    result_manual_block = henry_coefficient(framework, molecule, temperature, ljff,
                               insertions_per_volume=insertions_per_volume, verbose=true,
                               accessibility_grid=nothing)
    @test isapprox(result["⟨U⟩ (kJ/mol)"], result_manual_block["⟨U⟩ (kJ/mol)"], rtol=0.01)
    @test isapprox(result["henry coefficient [mmol/(g-bar)]"], result_manual_block["henry coefficient [mmol/(g-bar)]"], rtol=0.01)
end
