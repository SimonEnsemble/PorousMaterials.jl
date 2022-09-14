```@meta
DocTestSetup = quote
  using PorousMaterials
end
```

# Find Potential Energy Minimum

Here we show how to find the minimum energy (acc. to a force field) position of a molecule in a crystal.

### Example

For example, we wish to find the minimum energy (acc. to the UFF) position of a xenon adsorbate in SBMOF-1. 

```jldoctest
xtal = Crystal("SBMOF-1.cif")
molecule  = Molecule("Xe")
ljff = LJForceField("UFF")

# grid search to find min energy position.
#  gives good starting guess for optimization algorithm to fine tune.
resolution = 1.0 # resolution of grid points in Å
minimized_molecule, min_E = find_energy_minimum_gridsearch(xtal, molecule, ljff, resolution=resolution)
# minimized_molecule: xenon at its min energy position
# min_E: associated minimum energy of xenon (kJ/mol)

# fine tune the minimum energy position according to the grid search.
minimized_molecule, min_E = find_energy_minimum(xtal, minimized_molecule, ljff)
# output
Computing energy grid of Xe in SBMOF-1.cif
	Regular grid (in fractional space) of 13 by 7 by 24 points superimposed over the unit cell.
(Molecule species: Xe
Center of mass (fractional coords): Frac([0.01749943805846959; 0.9372916114895011; 0.011192272400742498;;])
Atoms:

	atom = Xe, xf = [0.017, 0.937, 0.011], -37.69376112588296)
```

# detailed docs
```@docs
    find_energy_minimum
    find_energy_minimum_gridsearch
```
