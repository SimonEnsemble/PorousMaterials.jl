## Grand-canonical Monte Carlo Simulations

Simulate the adsorption of CO$_2$ in FIQCEN\_clean\_min\_charges (CuBTC) at 298 K at 1 bar using the Universal Force Field:

```julia
using PorousMaterials

# read in xtal structure file and populate a Framework data structure
xtal = Crystal("FIQCEN_clean_min_charges.cif")
# remove numbers from atom labels
strip_numbers_from_atom_labels!(xtal)

# read in Lennard-Jones force field parameters and populate a LJForceField data structure
ljforcefield = LJForceField("UFF", r_cutoff=12.8)

# read in a molecule format file and populate a Molecule data structure
molecule = Molecule("CO2")

temperature = 298.0 # K
pressure = 1.0 # bar

# conduct Grand-Canonical Monte Carlo simulation
results, molecules = μVT_sim(xtal, molecule, temperature, pressure, forcefield,
            n_burn_cycles=5000, n_sample_cycles=5000)

# ... prints stuff
# results automatically saved to .jld2 load later in one line of code

# returns dictionary for easy querying
results["⟨N⟩ (molecules/unit cell)"]   
results["err ⟨N⟩ (molecules)"]         
results["Q_st (K)"]                     
```

Or, compute the entire adsorption isotherm at once, parallelized across many cores (this works by cleverly queuing a [`μVT_sim`](@ref) for each pressue across the specified number of cores for optimal efficiency):
```julia
pressures = [0.2, 0.6, 0.8, 1.0] # bar

# loop over all pressures and compute entire adsorption isotherm in parallel
results = adsorption_isotherm(xtal, molecule, temperature, pressures, forcefield,
            n_burn_cycles=5000, n_sample_cycles=5000)
```
 
Or, compute the adsorption isotherm in a step-wise manner, loading the molecules from the previous simulation to save on burn cycles:
```julia
# loop over all pressures and run GCMC simulations in series.
# load in the configurations of the molecules from the previous pressure.
results = stepwise_adsorption_isotherm(xtal, molecule, temperature, pressures, forcefield,
               n_burn_cycles=1000, n_sample_cycles=5000)
```
# detailed docs

## Grand-Canonical Monte Carlo Simulations
```@docs
    μVT_sim
    adsorption_isotherm
    stepwise_adsorption_isotherm
    μVT_output_filename
    isotherm_sim_results_to_dataframe
```

