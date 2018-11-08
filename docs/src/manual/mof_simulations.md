## Demo of Henry coefficients Calculation

Compute the Henry coefficient of CO<sub>2</sub> in CAXVII\_clean (Fe<sub>2</sub>(dobdc)) at 298 K using the Dreiding force field:

```julia
using PorousMaterials

# read in xtal structure file and populate a Framework data structure
framework = Framework("CAXVII_clean.cif")                                               

# read in Lennard-Jones force field parameters and populate a LJForceField data structure
forcefield = LJForceField("Dreiding.csv", cutoffradius=12.5)                                  

# read in a molecule format file and populate a Molecule data structure
molecule = Molecule("CO2")                                                              

temperature = 298.0 # K

# conduct Widom insertions and compute Henry coefficient, heat of adsorption
results = henry_coefficient(framework, molecule, temperature, forcefield, insertions_per_volume=200)

# ... prints stuff
# results automatically saved to .jld load later in one line of code

# returns dictionary for easy querying
results["Qst (kJ/mol)"] # -21.0
results["henry coefficient [mol/(kg-Pa)]"] # 2.88e-05
```

The simulation is parallelized across a maximum of 5 cores.

## Demo of Grand-canonical Monte Carlo Simulations

Simulate the adsorption of CO<sub>2</sub> in FIQCEN\_clean\_min\_charges (CuBTC) at 298 K at 1 bar using the Universal Force Field:

```julia
using PorousMaterials

# read in xtal structure file and populate a Framework data structure
framework = Framework("FIQCEN_clean_min_charges.cif")
# remove annoying numbers from atom labels
strip_numbers_from_atom_labels!(framework)

# read in Lennard-Jones force field parameters and populate a LJForceField data structure
forcefield = LJForceField("UFF.csv", cutoffradius=12.8)

# read in a molecule format file and populate a Molecule data structure
molecule = Molecule("CO2")

temperature = 298.0 # K
pressure = 1.0 # bar

# conduct grand-canonical Monte Carlo simulation
results, molecules = gcmc_simulation(framework, molecule, temperature, pressure, forcefield,
            n_burn_cycles=5000, n_sample_cycles=5000)

# ... prints stuff
# results automatically saved to .jld load later in one line of code

# returns dictionary for easy querying
results["⟨N⟩ (molecules/unit cell)"]
results["Q_st (K)"]
```

Or, compute the entire adsorption isotherm at once, parallelized across many cores:
```julia
pressures = [0.2, 0.6, 0.8, 1.0] # bar

# loop over all pressures and compute entire adsorption isotherm in parallel
results = adsorption_isotherm(framework, molecule, temperature, pressures, forcefield,
            n_burn_cycles=5000, n_sample_cycles=5000)
```

Or, compute the adsorption isotherm in a step-wise manner, loading the molecules from the previous simulation to save on burn cycles:
```julia
# loop over all pressures and run GCMC simulations in series.
# load in the configurations of the molecules from the previous pressure.
results = stepwise_adsorption_isotherm(framework, molecule, temperature, pressures, forcefield,
            n_burn_cycles=1000, n_sample_cycles=5000)
```

## Henry Coefficient Calculations
```@docs
    henry_coefficient
    henry_result_savename
```

## Grand-Canonical Monte Carlo Simulations
```@docs
    gcmc_simulation
    adsorption_isotherm
    stepwise_adsorption_isotherm
    gcmc_result_savename
```
