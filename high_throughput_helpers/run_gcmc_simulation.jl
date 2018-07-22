# example use, e.g. to run μVT sim of CH4 in CAXVII_clean.cif at 1 bar, 298 K:
#   julia run_gcmc_simulation.jl CAXVII_clean CH4 298.0 1.0

using PorousMaterials
using JLD

# read in command line arguments
if length(ARGS) != 4
    error("Run as: julia run_gcmc_simulation.jl structure gas temperature(K) pressure(bar)")
end

structurename = ARGS[1]
gasname = ARGS[2]
temperature = parse(Float64, ARGS[3])
pressure = parse(Float64, ARGS[4])

# simulation set-up
ff = LJForceField("Greg_bogus_ZIF71.csv", cutoffradius=12.8)
n_burn_cycles = 5000
n_sample_cycles = 5000

xtal = Framework(structurename * ".cif")
strip_numbers_from_atom_labels!(xtal)
adsorbate = Molecule(gasname, xtal.box)

# run the simulation
results = gcmc_simulation(xtal, temperature, pressure, adsorbate, ljff, 
                          n_burn_cycles=n_burn_cycles, n_sample_cycles=n_sample_cycles, 
                          verbose=true, sample_frequency=1, eos=:PengRobinson, autosave=true)

# results dictionary autosaved in data/gcmc_sims
