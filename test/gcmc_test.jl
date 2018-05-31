@everywhere using PorousMaterials
@everywhere using Base.Test
using CSV
using DataFrames

ig_tests = false
xe_in_sbmof1_tests = false
co2_tests = false

#
# Ideal gas tests.
#  run GCMC in empty box; should get ideal gas law.
#  "ig" in test_forcefield.csv has sigma tiny and epsilon 0.0 to match ideal gas.
#  basically, this tests the acceptance rules when energy is always zero.
#
if ig_tests
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

## Xe adsorption in SBMOF-1 tests
if xe_in_sbmof1_tests
    sbmof1 = read_crystal_structure_file("SBMOF-1.cif")
    dreiding_forcefield = read_forcefield_file("Dreiding.csv", cutoffradius=12.5)
    molecule = read_molecule_file("Xe")

    test_fugacities = [20.0, 200.0, 2000.0]
    test_mmol_g = [0.1931, 1.007, 1.4007]
    test_molec_unit_cell = [0.266, 1.388, 1.929]

    results = adsorption_isotherm(sbmof1, 298.0, test_fugacities, molecule, dreiding_forcefield, n_burn_cycles=25000, n_sample_cycles=25000, verbose=true)

    for i = 1:length(test_fugacities)
        @test isapprox(results[i]["⟨N⟩ (molecules/unit cell)"], test_molec_unit_cell[i], rtol=0.025)
        @test isapprox(results[i]["⟨N⟩ (mmol/g)"], test_mmol_g[i], rtol=0.025)
    end
end

if co2_tests
    # load test data, adsorption isotherm at 313 K simulated by Greg
    df = CSV.read("greg_chung/ZnCo-ZIF-71_atom_relax_RESP_CO2_adsorption_isotherm313K_test.csv")

    co2 = read_molecule_file("CO2")
    f = read_crystal_structure_file("ZnCo-ZIF-71_atom_relax_RESP.cif")
    strip_numbers_from_atom_labels!(f)
    ff = read_forcefield_file("Greg_CO2_GCMCtest_ff.csv", cutoffradius=12.5)
    
    # simulate with PorousMaterials.jl in parallel
    results = adsorption_isotherm(f, 313.0, convert(Array{Float64, 1}, df[:P_Pa]), co2, ff, n_burn_cycles=5000, n_sample_cycles=5000, verbose=true)
    n_sim = [result["⟨N⟩ (molecules/unit cell)"] for result in results]
    
    # plot comparison
    using PyPlot
    figure()
    xlabel("Pressure (bar)")
    ylabel("Molecules/unit cell")
    scatter(df[:P_Pa], df[:molecules_per_unit_cell], label="Greg")
    scatter(df[:P_Pa], n_sim, label="PorousMaterials.jl")
    legend()
    savefig("ZnCo-ZIF-71_atom_relax_RESP_CO2_adsorption_isotherm313K_test.png", format="png", dpi=300)
end
        
co2 = read_molecule_file("CO2")
f = read_crystal_structure_file("ZnCo-ZIF-71_atom_relax_RESP.cif")
strip_numbers_from_atom_labels!(f)
ff = read_forcefield_file("Greg_CO2_GCMCtest_ff.csv", cutoffradius=12.5)

results = gcmc_simulation(f, 313.0, 20.0*100000, co2, ff,
            n_burn_cycles=1, n_sample_cycles=1, verbose=true)
@time results = gcmc_simulation(f, 313.0, 1816566.334, co2, ff,
            n_burn_cycles=10000, n_sample_cycles=10000, verbose=true)
