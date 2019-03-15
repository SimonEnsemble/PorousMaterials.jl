using Distributed
@everywhere using PorousMaterials
@everywhere using Test
using CSV
using DataFrames
using JLD2
using Printf
using Statistics

ig_tests = false
xe_in_sbmof1_tests = true
co2_tests = true
density_grid_test = false

#
# Ideal gas tests.
#  run GCMC in empty box; should get ideal gas law.
#  "ig" in test_forcefield.csv has sigma tiny and epsilon 0.0 to match ideal gas.
#  basically, this tests the acceptance rules when energy is always zero.
#
if ig_tests
    empty_space = Framework("empty_box.cssr") # zero atoms!
    ideal_gas = Molecule("IG")
    @assert empty_space.atoms.n_atoms == 0
    forcefield = LJForceField("Dreiding.csv")
    temperature = 298.0
    fugacity = 10.0 .^ [0.1, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0] / 100000.0 # bar
    # according to ideal gas law, number of molecules in box should be:
    n_ig = fugacity * empty_space.box.Ω / (PorousMaterials.KB * temperature)
    n_sim = similar(n_ig)
    for i = 1:length(fugacity)
        results, molecules = gcmc_simulation(empty_space, ideal_gas, temperature, fugacity[i], forcefield,
                    n_burn_cycles=25000, n_sample_cycles=25000)
#                    n_burn_cycles=100000, n_sample_cycles=100000)
        n_sim[i] = results["⟨N⟩ (molecules/unit cell)"]
        @printf("fugacity = %f Pa; n_ig = %e; n_sim = %e\n", fugacity[i], n_ig[i], n_sim[i])
    end
end

if density_grid_test
    # test denstiy grid w./ ideal gas
    n_pts = (4, 4, 4) # testing a grid with 4x4x4 voxels
    empty_space = Framework("empty_box.cssr") # zero atoms!
    ideal_gas = Molecule("IG")
    forcefield = LJForceField("Dreiding.csv", cutoffradius=4.0)
    temperature = 25.0
    fugacity = 500.0 # mk sure lots of molecules
    results, molecules = gcmc_simulation(empty_space, ideal_gas, temperature, fugacity, forcefield,
                n_burn_cycles=1, n_sample_cycles=1000, calculate_density_grid=true, density_grid_dx=10.1,
                show_progress_bar=true)
    mean_density = mean(results["density grid"].data)
    println(results["density grid"].data)
    @test all(isapprox.(results["density grid"].data, mean_density, rtol=0.05))
end

## Xe adsorption in SBMOF-1 tests
if xe_in_sbmof1_tests
    sbmof1 = Framework("SBMOF-1.cif")
    dreiding_forcefield = LJForceField("Dreiding.csv", cutoffradius=12.5)
    molecule = Molecule("Xe")

    test_fugacities = [20.0, 200.0, 2000.0] / 100000.0 # bar
    test_mmol_g = [0.18650, 1.00235, 1.39812]
    test_molec_unit_cell = [0.2568, 1.3806, 1.9257]

 #     results = adsorption_isotherm(sbmof1, 298.0, test_fugacities, molecule, dreiding_forcefield, n_burn_cycles=20, n_sample_cycles=20, verbose=true)
    results = stepwise_adsorption_isotherm(sbmof1, molecule, 298.0, test_fugacities, dreiding_forcefield,
                        n_burn_cycles=25000, n_sample_cycles=25000, verbose=true, sample_frequency=1, show_progress_bar=true)

    for i = 1:length(test_fugacities)
        @test isapprox(results[i]["⟨N⟩ (molecules/unit cell)"], test_molec_unit_cell[i], rtol=0.025)
        @test isapprox(results[i]["⟨N⟩ (mmol/g)"], test_mmol_g[i], rtol=0.025)
    end
end

if co2_tests
    plot_results = false
    run_sims = true
    ###
    #  Test isotherm 1: by greg chung. co2 at 313 k
    ###
 #     f = Framework("ZnCo-ZIF-71_atom_relax_RESP.cif")
 #     co2 = Molecule("CO2")

 #     strip_numbers_from_atom_labels!(f)
 #     ff = LJForceField("Greg_CO2_GCMCtest_ff.csv", cutoffradius=12.5)
 #
 #     # load in test data
 #     df = CSV.read("greg_chung/ZnCo-ZIF-71_atom_relax_RESP_CO2_adsorption_isotherm313K_test.csv")
 #
 #     # simulate with PorousMaterials.jl in parallel
 #     if run_sims
 #         results = adsorption_isotherm(f, co2, 313.0, convert(Array{Float64, 1}, df[:fugacity_Pa] / 100000.0), ff,
 #             n_burn_cycles=10000, n_sample_cycles=10000, verbose=true, sample_frequency=5)
 #         JLD.save("ZnCo-ZIF-71_atom_relax_RESP_co2_simulated_isotherm.jld", "results", results)
 #     else
 #         results = JLD.load("ZnCo-ZIF-71_atom_relax_RESP_co2_simulated_isotherm.jld")["results"]
 #     end
 #     n_sim = [result["⟨N⟩ (molecules/unit cell)"] for result in results]
 #
 #     # plot comparison
 #     if plot_results
 #         figure()
 #         xlabel("Fugacity (bar)")
 #         ylabel("Molecules/unit cell")
 #         scatter(df[:fugacity_Pa] / 100000.0, df[:molecules_per_unit_cell], label="Greg")
 #         scatter(df[:fugacity_Pa] / 100000.0, n_sim, label="PorousMaterials.jl")
 #         legend()
 #         savefig("ZnCo-ZIF-71_atom_relax_RESP_CO2_adsorption_isotherm313K_test.png", format="png", dpi=300)
 #     end

    ###
    #  Test isotherm 2: by greg chung. co2 at 298 K
    ###
    zif71 = Framework("zif71_bogus_charges.cif")
    strip_numbers_from_atom_labels!(zif71)
    ff = LJForceField("Greg_bogus_ZIF71.csv", cutoffradius=12.8)
    co2 = Molecule("CO2EPM2")

    # make sure bond lenghts are preserved
    bls = pairwise_atom_distances(co2, UnitCube())
    results, molecules = gcmc_simulation(zif71, co2, 298.0, 1.0, ff,
                        n_burn_cycles=25, n_sample_cycles=25, verbose=false)
    @printf("Testing that bond lenghts are preserved for %d molecules.\n", length(molecules))
    for m in molecules
        # bond lengths preserved?
        @assert isapprox(bls, pairwise_atom_distances(m, UnitCube()))
        # charges hv same coords as atoms?
        for i = 1:3
            @assert isapprox(m.atoms.xf[:, i], m.charges.xf[:, i])
        end
        @assert isapprox(m.xf_com, m.atoms.xf[:, 1]) # C atom is center
    end


    # load in test data
    df = CSV.read("greg_chung/zif_71_co2_isotherm_w_preos_fugacity.csv")

    # simulate with PorousMaterials.jl in parallel
    if run_sims
        results = stepwise_adsorption_isotherm(zif71, co2, 298.0, convert(Array{Float64, 1}, df[:fugacity_Pa] / 100000.0), ff,
            n_burn_cycles=2000, n_sample_cycles=5000, verbose=true, sample_frequency=1, ewald_precision=1e-6)
        @save "ZIF71_bogus_charges_co2_simulated_isotherm.jld" results
    else
        @load "ZIF71_bogus_charges_co2_simulated_isotherm.jld" results
    end
    n_sim = [result["⟨N⟩ (mmol/g)"] for result in results]

    # plot comparison
    if plot_results
        using PyPlot
        figure()
        xlabel("Fugacity (bar)")
        ylabel("Molecules/unit cell")
        scatter(df[:fugacity_Pa] / 100000.0, df[:L_mmol_g], label="Greg")
        scatter(df[:fugacity_Pa] / 100000.0, n_sim, label="PorousMaterials.jl")
        legend()
        savefig("Greg_bogus_ZIF71_298K_co2_isotherm_test.png", format="png", dpi=300)
    end
end

 # co2 = Molecule("CO2")
 # f = Framework("ZnCo-ZIF-71_atom_relax_RESP.cif")
 # strip_numbers_from_atom_labels!(f)
 # ff = LJForceField("Greg_CO2_GCMCtest_ff.csv", cutoffradius=12.5)
 #
 # results = gcmc_simulation(f, 313.0, 20.0*100000, co2, ff,
 #             n_burn_cycles=1, n_sample_cycles=1, verbose=true)
 # @time results = gcmc_simulation(f, 313.0, 1816566.334, co2, ff,
 #             n_burn_cycles=10000, n_sample_cycles=10000, verbose=true)
