using PorousMaterials
using Base.Test
using Base.Random
using CSV
using DataFrames
using JLD


frame = Framework("FIQCEN_clean_min_charges.cif")
strip_numbers_from_atom_labels!(frame)
co2 = Molecule("CO2")
co_bond_length = norm(co2.atoms[1].xf - co2.atoms[2].xf)
ljff = LJForceField("UFF.csv")
temp = 298.0
pressure = 0.5
molecules = Molecule[]

results = Dict()
srand(1234)
println(" a random no: ", rand())

for i = 1:3
    @printf("i = %d\n",i)
    if i == 1
        checkpoint = Dict()
    else
        checkpoint = JLD.load(PorousMaterials.PATH_TO_DATA * "/gcmc_checkpoints/" * 
            gcmc_result_savename(frame.name, co2.species, ljff.name, temp, pressure, 5, 5 * i, "_checkpoint"), "checkpoint")
    end
    results, molecules = gcmc_simulation(frame, co2, temp, pressure, ljff,
                                         n_burn_cycles=5, n_sample_cycles=5 * (i + 1),
                                         verbose=true, sample_frequency=1, eos=:PengRobinson,
                                         autosave=false, write_checkpoints=true, 
                                         checkpoint=checkpoint, checkpoint_frequency=1)
    @assert isapprox(norm(molecules[1].atoms[1].xf - molecules[1].atoms[2].xf), co_bond_length)
end

println("---------------------------------")
println(" a random no at end: ", rand())

srand(1234)
println(" a random no: ", rand())

results2, molecules2 = gcmc_simulation(frame, co2, temp, pressure, ljff,
                          n_burn_cycles=5, n_sample_cycles=20,
                          verbose=true, sample_frequency=1, eos=:PengRobinson,
                          autosave=false, write_checkpoints=false, load_checkpoint=false, checkpoint_frequency=1, progressbar=true)
@assert isapprox(norm(molecules2[1].atoms[1].xf - molecules2[1].atoms[2].xf), co_bond_length)
 # @test all(isapprox.(molecules, molecules2))
println(" a random no at end: ", rand())
