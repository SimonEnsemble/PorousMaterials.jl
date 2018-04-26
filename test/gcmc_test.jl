using PorousMaterials
using Base.Test

run_ig_tests = true

#
# Ideal gas tests.
#  run GCMC in empty box; should get ideal gas law.
#  "ig" in test_forcefield.csv has sigma tiny and epsilon 0.0 to match ideal gas.
#  basically, this tests the acceptance rules when energy is always zero.
#
if run_ig_tests
    empty_space = read_crystal_structure_file("empty_box.cssr") # zero atoms!
    ideal_gas = read_molecule_file("IG")
    @assert(empty_space.n_atoms == 0)
    forcefield = read_forcefield_file("Dreiding.csv")
    temperature = 298.0
    fugacity = 10.0 .^ [0.1, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    # according to ideal gas law, number of molecules in box should be:
    n_ig = fugacity * empty_space.box.Ω / (PorousMaterials.KB * temperature)
    n_sim = similar(n_ig)
    for i = 1:length(fugacity)
        results = gcmc_simulation(empty_space, temperature, fugacity[i], ideal_gas, forcefield,
                    n_burn_cycles=100000, n_sample_cycles=100000)
        n_sim[i] = results["⟨N⟩ (molecules/unit cell)"]
        @printf("fugacity = %f Pa; n_ig = %e; n_sim = %e\n", fugacity[i], n_ig[i], n_sim[i])
    end
end

## SBMOF-1 tests
sbmof1 = read_crystal_structure_file("SBMOF-1.cif")
dreiding_forcefield = read_forcefield_file("Dreiding.csv", cutoffradius=12.5)

##Getting template for use with the sim
molecule = read_molecule_file("Xe")

# very quick test
results = gcmc_simulation(sbmof1, 298.0, 2300.0, molecule, dreiding_forcefield, n_burn_cycles=10, n_sample_cycles=10, verbose=true)


test_fugacities = [20.0, 200.0, 2000.0]
test_mmol_g = [0.1931, 1.007, 1.4007]
test_molec_unit_cell = [0.266, 1.388, 1.929]

for (i, fugacity) in enumerate(test_fugacities)
    @time results = gcmc_simulation(sbmof1, 298.0, fugacity, molecule, dreiding_forcefield, n_burn_cycles=10000, n_sample_cycles=10000, verbose=true)
    println("Should be (molec/unit cell: ", test_molec_unit_cell[i])
    @test isapprox(results["⟨N⟩ (molecules/unit cell)"], test_molec_unit_cell[i], rtol=0.01)
    @test isapprox(results["⟨N⟩ (mmol/g)"], test_mmol_g[i], rtol=0.01)
end
