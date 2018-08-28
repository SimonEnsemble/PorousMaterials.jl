# example use, e.g. to run μVT sim of CH4 in CAXVII_clean.cif at 1 bar, 298 K:
#   julia run_gcmc_simulation.jl CAXVII_clean CH4 298.0 1.0

# read in command line arguments
if length(ARGS) != 4
    error("Run as: julia run_gcmc_simulation.jl structure gas temperature(K) pressure(bar)")
end

structurename = ARGS[1]
gasname = ARGS[2]
temperature = parse(Float64, ARGS[3])
pressure = parse(Float64, ARGS[4])

include("jobs_to_run.jl") # contains ljff and # samples

# simulation set-up
ljff = jobs_to_run[structurename]["forcefield"]
n_burn_cycles = jobs_to_run[structurename]["n_burn_cycles"]
n_sample_cycles = jobs_to_run[structurename]["n_sample_cycles"]

xtal = Framework(structurename * ".cif")
strip_numbers_from_atom_labels!(xtal)
adsorbate = Molecule(gasname)

# run the simulation
results = gcmc_simulation(xtal, temperature, pressure, adsorbate, ljff, 
                          n_burn_cycles=n_burn_cycles, n_sample_cycles=n_sample_cycles, 
                          verbose=true, sample_frequency=1, eos=:PengRobinson, autosave=true)

# results dictionary autosaved in data/gcmc_sims
