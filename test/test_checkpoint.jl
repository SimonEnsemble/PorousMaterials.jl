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

results = Dict()
srand(1234)

for i = 1:3
    @printf("i = %d\n",i)
    checkpoint_file = i == 1 ? false : gcmc_result_savename(frame.name, co2.species, ljff.name, temp, pressure, 5, 5*(i-1), "_checkpoint")
    results, molecules = gcmc_simulation(frame, co2, temp, pressure, ljff,
                                         n_burn_cycles = 5, n_sample_cycles=5*i,
                                         verbose=true, sample_frequency=1, eos=:PengRobinson,
                                         autosave=false, write_checkpoints=true, load_checkpoint=checkpoint_file, checkpoint_frequency=1)
    @test isapprox(norm(molecules[1].atoms[1].xf - molecules[1].atoms[2].xf), co_bond_length)
end

println("---------------------------------")

srand(1234)

results2, molecules2 = gcmc_simulation(frame, co2, temp, pressure, ljff,
                          n_burn_cycles=5, n_sample_cycles = 15,
                          verbose=true, sample_frequency=1, eos=:PengRobinson,
                          autosave=false, write_checkpoints=false, load_checkpoint=false, checkpoint_frequency=1, progressbar=true)
 # @test all(isapprox.(molecules, molecules2))
