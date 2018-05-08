using PorousMaterials

co2 = read_molecule_file("CO2")
xtal = read_crystal_structure_file("ZnCo-ZIF-71_atom_relax_RESP.cif")
strip_numbers_from_atom_labels!(xtal)
ff = read_forcefield_file("Greg_CO2_GCMCtest_ff.csv", cutoffradius=12.8)
results = gcmc_simulation(xtal, 313.0, 20000.0, co2, ff, n_burn_cycles=1000, n_sample_cycles=10000, verbose=true)
