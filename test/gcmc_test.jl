using PorousMaterials
sbmof1 = read_crystal_structure_file("SBMOF-1.cif")
dreiding_forcefield = read_forcefield_file("test_forcefield.csv")
 
# very quick test
results = gcmc_simulation(sbmof1, 298.0, 2300.0, "Xe", dreiding_forcefield, n_burn_cycles=10, n_sample_cycles=10, verbose=true)

test_fugacities = [20.0, 200.0, 2000.0]
test_mmol_g = [0.1904, 1.007, 1.4007]
test_molec_unit_cell = [0.262, 1.388, 1.929]

for (i, fugacity) in enumerate(test_fugacities)
    results = gcmc_simulation(sbmof1, 298.0, fugacity, "Xe", dreiding_forcefield, n_burn_cycles=5000, n_sample_cycles=5000, verbose=true)
    isapprox(results["⟨N⟩ (molecules/unit cell)"], test_molec_unit_cell[i], rtol=0.005)
    isapprox(results["⟨N⟩ (mmol/g)"], test_mmol_g[i], rtol=0.005)
end
