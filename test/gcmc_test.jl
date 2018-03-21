using PorousMaterials
f = read_crystal_structure_file("SBMOF-1.cif")
ff = read_forcefield_file("test_forcefield.csv")
results = gcmc_simulation(f, 298.0, 10000.0, "Xe", ff, n_burn_cycles=20, n_mc_proposals=200, verbose=true)
