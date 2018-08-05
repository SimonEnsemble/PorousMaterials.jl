using PorousMaterials
using Base.Test
using Base.Random
using CSV
using DataFrames
using JLD


frame = Framework("FIQCEN_clean.cif")
strip_numbers_from_atom_labels!(frame)
co2 = Molecule("CO2")
ljff = LJForceField("UFF.csv")
n_burn_cycles = 51
n_sample_cycles = 51
temp = 298.0
pressure = 0.5

srand(1234)

try
    results, molecules = gcmc_simulation(frame, co2, temp, pressure, ljff,
                              n_burn_cycles=51, n_sample_cycles = 51,
                              verbose=true, sample_frequency=1, eos=:PengRobinson,
                              autosave=false, write_checkpoint=true, read_checkpoint=false, checkpoint_interval=1, troubleshoot=true)
end
results, molecules = gcmc_simulation(frame, co2, temp, pressure, ljff,
                          n_burn_cycles=51, n_sample_cycles = 51,
                          verbose=true, sample_frequency=1, eos=:PengRobinson,
                          autosave=false, write_checkpoint=true, read_checkpoint=true, checkpoint_interval=1, troubleshoot=false)
println("---------------------------------")

srand(1234)

results2, molecules2 = gcmc_simulation(frame, co2, temp, pressure, ljff,
                          n_burn_cycles=n_burn_cycles, n_sample_cycles = n_sample_cycles,
                          verbose=true, sample_frequency=1, eos=:PengRobinson,
                          autosave=false, write_checkpoint=false, read_checkpoint=false, checkpoint_interval=1)
