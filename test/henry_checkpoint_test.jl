using PorousMaterials
using JLD2
using FileIO
using Test

@testset "Henry Checkpoint Tests" begin
    framework = Framework("FIQCEN_clean_min_charges.cif")
    strip_numbers_from_atom_labels!(framework)
    co2 = Molecule("CO2")
    ljff = LJForceField("UFF.csv")
    temp = 298.0

    # trick to ensure there are 750 total insertions
    #  and thus 150 insertions per block, divisible by 25 which is the
    #  check point writing frequency below.
    insertions_per_volume = 750.0 / framework.box.Ω
    
    # the first simulation writes checkpoints but does not load them
    # we allow it to completely finish. then the checkpoints should be saved.
    results = henry_coefficient(framework, co2, temp, ljff, 
                                insertions_per_volume=insertions_per_volume, 
                                write_checkpoint=true, checkpoint_frequency=25,
                                load_checkpoint=false)
    
    # the second simulation is the same as the first except we load in the checkpoint.
    #   from the first simulation.
    #   since the first simulation was allowed to finish, this just returns the loaded results.
    #   so we're making sure this reads and writes the checkpoint correctly.
    #   TODO maybe a more rigorous test where we actually carry off from an interrupted sim
    #    but this is difficult b/c the random numbers need to correspond and the history
    #    does implicitly matter b/c we rotate about a starting config.
    results2 = henry_coefficient(framework, co2, temp, ljff,
                                 insertions_per_volume=insertions_per_volume,
                                 load_checkpoint=true)

    # we expect the last simulation to have loaded in the checkpoint from the first and thus
    #  be equivalent. (this would not be case if insertions per block were not
    #  exactly divisible by checkpoint_frequency!
    @test(isapprox(results["⟨U⟩ (K)"], results2["⟨U⟩ (K)"]))
    @test(isapprox(results["⟨U, vdw⟩ (K)"], results2["⟨U, vdw⟩ (K)"]))
    @test(isapprox(results["⟨U, Coulomb⟩ (K)"], results2["⟨U, Coulomb⟩ (K)"]))
    @test(isapprox(results["Qst (kJ/mol)"], results2["Qst (kJ/mol)"]))
    @test(isapprox(results["henry coefficient [mol/(m³-bar)]"], results2["henry coefficient [mol/(m³-bar)]"]))
end
