## Grand-canonical Monte Carlo Simulations

Simulate the adsorption of CO$_2$ in FIQCEN\_clean\_min\_charges (CuBTC) at 298 K at 1 bar using the Universal Force Field:

```julia
xtal = Crystal("FIQCEN_clean.cif") # load a crystal structure
strip_numbers_from_atom_labels!(xtal) # clean up the atom labels
ljforcefield = LJForceField("UFF", r_cutoff=12.8) # load the UFF forcefield
molecule = Molecule("CO2") # load the CO2 molecule
temperature = 298.0 # K
pressure = 1.0 # bar
# conduct Grand-Canonical Monte Carlo simulation (VERY short, should use thousands of cycles!)
results, molecules = μVT_sim(xtal, molecule, temperature, pressure, ljforcefield,
            n_burn_cycles=50, n_sample_cycles=50)
```

Or, compute the entire adsorption isotherm at once, parallelized across many cores (this works by cleverly queuing a [`μVT_sim`](@ref) for each pressue across the specified number of cores for optimal efficiency):

```julia
pressures = [0.2, 0.6, 0.8, 1.0] # bar
# loop over all pressures and compute entire adsorption isotherm in parallel
results = adsorption_isotherm(xtal, molecule, temperature, pressures, ljforcefield,
            n_burn_cycles=50, n_sample_cycles=50)
```
 
Or, compute the adsorption isotherm in a step-wise manner, loading the molecules from the previous simulation to save on burn cycles:

```julia
results = stepwise_adsorption_isotherm(xtal, molecule, temperature, pressures, forcefield,
               n_burn_cycles=50, n_sample_cycles=50)
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

