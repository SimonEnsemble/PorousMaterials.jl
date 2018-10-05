using PorousMaterials
using JLD2
using FileIO
using Test

@testset "Henry Checkpoint Tests" begin
    # Write checkpoint files during a henry simulation. Then load it in and compare the results
    framework = Framework("FIQCEN_clean_min_charges.cif")
    strip_numbers_from_atom_labels!(framework)
    co2 = Molecule("CO2")
    ljff = LJForceField("UFF.csv")
    temp = 298.0
    insertions_per_volume = 750.0 / framework.box.Ω
    results = Dict()

    results = henry_coefficient(framework, co2, temp, ljff, 
                                insertions_per_volume=insertions_per_volume, 
                                write_checkpoint=true, checkpoint_frequency=25,
                                load_checkpoint=false)

    results2 = henry_coefficient(framework, co2, temp, ljff,
                                 insertions_per_volume=insertions_per_volume,
                                 load_checkpoint=true)

    @test(isapprox(results["⟨U⟩ (K)"], results2["⟨U⟩ (K)"]))
    @test(isapprox(results["⟨U, vdw⟩ (K)"], results2["⟨U, vdw⟩ (K)"]))
    @test(isapprox(results["⟨U, Coulomb⟩ (K)"], results2["⟨U, Coulomb⟩ (K)"]))
    @test(isapprox(results["Qst (kJ/mol)"], results2["Qst (kJ/mol)"]))
    @test(isapprox(results["henry coefficient [mol/(m³-bar)]"], results2["henry coefficient [mol/(m³-bar)]"]))
end
