using PorousMaterials
using Random
using JLD2
using FileIO
using Test

framework = Framework("SBMOF-1.cif")
co2 = Molecule("CO2")
ljff = LJForceField("UFF.csv")
temp = 298.0
insertions_per_volume = 1500.0/framework.box.Ω
results = Dict()

Random.seed!(1234)
for i = 1:5
    if i == 1
        checkpoint = Dict()
    else
        checkpoint = load(joinpath(PorousMaterials.PATH_TO_DATA, "henry_checkpoints", henry_result_savename(framework, co2, temp, ljff, insertions_per_volume * (i - 1), comment="checkpoint")))
    end

    results = henry_coefficient(framework, co2, temp, ljff, insertions_per_volume = insertions_per_volume * i, write_checkpoint=true, checkpoint_frequency=10, checkpoint=checkpoint)
    print(results)
end

Random.seed!(1234)
results2 = henry_coefficient(framework, co2, temp, ljff, insertions_per_volume = 5 * insertions_per_volume)
