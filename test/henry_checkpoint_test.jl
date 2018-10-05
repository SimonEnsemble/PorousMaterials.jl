using PorousMaterials
using Random
using JLD2
using FileIO
using Test

framework = Framework("SBMOF-1.cif")
co2 = Molecule("CO2")
ljff = LJForceField("UFF.csv")
temp = 298.0
insertions_per_volume = 750.0/framework.box.Ω
ins_per_vol = [150.0, 300.0, 450.0, 600.0, 750.0] ./ framework.box.Ω
results = Dict()

Random.seed!(1234)
for i = 1:5
    global results
    results = henry_coefficient(framework, co2, temp, ljff, insertions_per_volume = ins_per_vol[i], write_checkpoint=true, checkpoint_frequency=10, load_checkpoint=true)
    print(results)
end

Random.seed!(1234)
results2 = henry_coefficient(framework, co2, temp, ljff, insertions_per_volume = insertions_per_volume)
