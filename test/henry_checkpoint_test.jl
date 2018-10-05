using PorousMaterials
using JLD2
using FileIO
using Test

@testset "Henry Checkpoint Tests" begin
    # Write checkpoint files during a henry simulation. Then load it in and compare the results
    framework = Framework("SBMOF-1.cif")
    co2 = Molecule("CO2")
    ljff = LJForceField("UFF.csv")
    temp = 298.0
    insertions_per_volume = 750.0/framework.box.Ω
    ins_per_vol = [150.0, 300.0, 450.0, 600.0, 750.0] ./ framework.box.Ω
    results = Dict()

    for i = 1:5
        #global results
        results = henry_coefficient(framework, co2, temp, ljff, insertions_per_volume = ins_per_vol[i], write_checkpoint=true, checkpoint_frequency=25, load_checkpoint=false)
    end

    results2 = henry_coefficient(framework, co2, temp, ljff, insertions_per_volume = insertions_per_volume, load_checkpoint=true)

    @test(isapprox(results["⟨U⟩ (K)"], results2["⟨U⟩ (K)"]))
    @test(isapprox(results["Qst (kJ/mol)"], results2["Qst (kJ/mol)"]))
    @test(isapprox(results["henry coefficient [mol/(m³-bar)]"], results2["henry coefficient [mol/(m³-bar)]"]))
end
