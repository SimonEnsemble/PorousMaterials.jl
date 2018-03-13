using PorousMaterials
f = read_crystal_structure_file("SBMOF-1.cif")
ff = read_forcefield_file("test_forcefield.csv")
gcmc_sim(f, 298.0, 100.0, "Xe", ff)
