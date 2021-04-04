# Find Potential Energy Minimum

Here we show how to find the minimum energy (acc. to a force field) position of a molecule in a crystal.

### Example
for example, we wish to find the minimum energy (acc. to the UFF) position of a xenon adsorbate in SBMOF-1. 

```julia
xtal = Crystal("SBMOF-1.cif")
molecule  = Molecule("Xe")
ljff = LJForceField("UFF")

# grid search to find min energy position.
#  gives good starting guess for optimization algorithm to fine tune.
n_pts = (15, 15, 15) # resolution of grid points
minimized_molecule, min_E = find_energy_minimum_gridsearch(xtal, molecule, ljff, n_pts=n_pts)
# minimized_molecule: xenon at its min energy position
# min_E: associated minimum energy of xenon (kJ/mol)

# fine tune the minimum energy position according to the grid search.
minimized_molecule, min_E = find_energy_minimum(xtal, minimized_molecule, ljff)
```

# detailed docs
```@docs
    find_energy_minimum
    find_energy_minimum_gridsearch
```
