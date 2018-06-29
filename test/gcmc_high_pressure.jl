@everywhere using PorousMaterials
@everywhere using Base.Test
using CSV
using PyPlot
using DataFrames
using JLD
    
    
zif71 = read_crystal_structure_file("zif71_bogus_charges.cif")
strip_numbers_from_atom_labels!(zif71)
ff = read_forcefield_file("Greg_bogus_ZIF71.csv", cutoffradius=12.8)
co2 = read_molecule_file("CO2EPM2")

# simulate with PorousMaterials.jl in parallel
results, molecules = gcmc_simulation(zif71, 298.0, 1786780.57, co2, ff, n_burn_cycles=15000, n_sample_cycles=15000, verbose=true, sample_frequency=1, ewald_precision=1e-6)
JLD.save("molecules.jld", "molecules", molecules)
