using PorousMaterials
using Test

# set up sim
xtal            = Crystal("SBMOF-1.cif")
molecules       = [Molecule("Kr"), Molecule("Xe")]
ljff            = LJForceField("UFF", mixing_rules="Lorentz-Berthelot")
temperature     = 298.0
n_sample_cycles = 5000
n_burn_cycles   = 5000
pressures       = [0.2, 0.8]

results, molecules = μVT_sim(xtal, molecules, temperature, pressures, ljff, n_burn_cycles=n_burn_cycles, n_sample_cycles=n_sample_cycles)
