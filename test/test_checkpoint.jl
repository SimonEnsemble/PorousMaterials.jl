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
n_burn_cycles = 100
n_sample_cycles = 100
temp = 298.0
pressure = 0.5

srand(1234)
try
    results = gcmc_simulation(frame, co2, temp, pressure, ljff,
                              n_burn_cycles=n_burn_cycles, n_sample_cycles = n_sample_cycles,
                              verbose=true, sample_frequency=1, eos=:PengRobinson,
                              autosave=false, checkpoint=true, checkpoint_interval=1,
                              troubleshoot_checkpoint=true)
end

    
results, molecules = gcmc_simulation(frame, co2, temp, pressure, ljff,
                          n_burn_cycles=n_burn_cycles, n_sample_cycles = n_sample_cycles,
                          verbose=true, sample_frequency=1, eos=:PengRobinson,
                          autosave=false, checkpoint=true, checkpoint_interval=1,
                          troubleshoot_checkpoint=false)

srand(1234)

results2, molecules2 = gcmc_simulation(frame, co2, temp, pressure, ljff,
                          n_burn_cycles=n_burn_cycles, n_sample_cycles = n_sample_cycles,
                          verbose=true, sample_frequency=1, eos=:PengRobinson,
                          autosave=false, checkpoint=false, checkpoint_interval=1,
                          troubleshoot_checkpoint=false)
