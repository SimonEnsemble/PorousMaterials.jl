using PorousMaterials
using Test
using Printf
using CSV
using DataFrames
using JLD2

tests_to_run = Dict("Kr/Xe in SBMOF-1"  => true,
                    "run mixture sims"  => false,
                    "run pure gas sims" => false,
                    "ideal gas mixture" => false
                    )

###
#  Simulate Kr/Xe mixture in SBMOF-1
###
if tests_to_run["Kr/Xe in SBMOF-1"]
    @info "running Kr/Xe in SBMOF-1 test"
    # set up sim
    xtal            = Crystal("SBMOF-1.cif")
    adsorbates      = ["Kr", "Xe"]
    mol_templates   = Molecule.(adsorbates)
    ljff            = LJForceField("UFF", mixing_rules="Lorentz-Berthelot")
    temperature     = 298.0
    # n_sample_cycles and n_burn_cycles are on the lower end of minimum requirment to get agreement with RASPA
    n_sample_cycles = 250000
    n_burn_cycles   = 250000
    mol_fraction    = [0.9, 0.1]
    # load RASPA (benchmark) data
    raspa_data = Dict(:Kr  => CSV.read(joinpath(pwd(), "data", "raspa_data/Kr.csv"), DataFrame),
                      :Xe  => CSV.read(joinpath(pwd(), "data", "raspa_data/Xe.csv"), DataFrame),
                      :mix => CSV.read(joinpath(pwd(), "data", "raspa_data/Kr_Xe.csv"), DataFrame))

    # run sim over range of total pressures
    for (i, p) in enumerate([0.01, 0.1, 0.5])
        if tests_to_run["run mixture sims"]
            @warn "Mixture simulation set will take about a day to run on a serial machine."
            results, molecules = μVT_sim(xtal,
                                         mol_templates,
                                         temperature,
                                         p *  mol_fraction,
                                         ljff,
                                         n_burn_cycles=n_burn_cycles,
                                         n_sample_cycles=n_sample_cycles)
        else
            # use the simulation files that are provided
            filename = μVT_output_filename(xtal, mol_templates, temperature, p * mol_fraction, ljff,
                                           n_burn_cycles, n_sample_cycles)
            @load joinpath(rc[:paths][:simulations], filename) results
        end

        # evaluate mixture results
        for (j, sp) in enumerate(adsorbates)
            # make sure the pressures are correct
            @test results["pressure (bar)"][j] ≈ raspa_data[:mix][i, Symbol(sp * " pressure (bar)")]
            # evaluate adsorption uptake
            std_err = raspa_data[:mix][i, Symbol("err $sp loading (mmol/g)")] + results["err ⟨N⟩ (mmol/g)"][j]
            @test isapprox(raspa_data[:mix][i, Symbol(sp * " loading (mmol/g)")], results["⟨N⟩ (mmol/g)"][j], atol=std_err)
        end

        # evaluate pure gas results, keep separate sinse we re-assign results dictionary
        for mol in mol_templates
            if  tests_to_run["run pure gas sims"]
                @warn "Running pure component simulation set will take about a day to run on a serial machine."
                results, molecules = μVT_sim(xtal,
                                             [mol],
                                             temperature,
                                             [p],
                                             ljff,
                                             n_burn_cycles=50000,
                                             n_sample_cycles=50000)
            else
                # load pure gas dictionary
                pure_gas_filename = μVT_output_filename(xtal, [mol], temperature, [p], ljff, 50000, 50000)
                @load joinpath(rc[:paths][:simulations],  pure_gas_filename) results
            end
            # make sure the pressures are correct
            @test results["pressure (bar)"][1] ≈ raspa_data[mol.species][i, Symbol(String(mol.species) * " pressure (bar)")]
            # evaluate adsorption uptake
            std_err = raspa_data[mol.species][i, Symbol("err $(String(mol.species)) loading (mmol/g)")] + results["err ⟨N⟩ (mmol/g)"][1]
            @test isapprox(raspa_data[mol.species][i, Symbol(String(mol.species) * " loading (mmol/g)")], results["⟨N⟩ (mmol/g)"][1], atol=std_err)
        end
    end
end

###
#  ideal gas mixture tests
###
#  run GCMC in empty box; should get ideal gas law.
#  "ig" in test_forcefield.csv has sigma tiny and epsilon 0.0 to match ideal gas.
#  basically, this tests the acceptance rules when energy is always zero.
if tests_to_run["ideal gas mixture"]
    @info "running ideal gas mixture test"
    empty_space  = replicate(Crystal("empty_box.cssr"), (3,3,3)) # zero atoms!
    ideal_gas1   = Molecule("ig")
    ideal_gas2   = Molecule("ig")
    @assert empty_space.atoms.n == 0
    forcefield   = LJForceField("Dreiding")
    temperature  = 298.0
    mol_fraction = [0.70, 0.30] 
    fugacities   = 10.0 .^ [2.0, 4.0, 5.0, 6.0, 7.0] / 100000.0 # bar
    partial_fugacities = [mol_fraction * f for f in fugacities]

    # according to ideal gas law (P_i*V=n_i*k_b*T), number of molecules per species in box should be
    n_ig  = partial_fugacities * empty_space.box.Ω / (PorousMaterials.BOLTZMANN * temperature) * 100000.0
    n_sim = similar(n_ig)
    for i in 1:length(partial_fugacities)
        @info partial_fugacities[i]
        results, molecules = μVT_sim(empty_space, 
                                     [ideal_gas1, ideal_gas2], 
                                     temperature, 
                                     partial_fugacities[i], 
                                     forcefield,
                                     n_burn_cycles=1000000, 
                                     n_sample_cycles=1000000, 
                                     verbose=false)
        # print fugacities 
        for s in 1:2
            @printf("\nfugacities %f, n_ig = %f, n_sim = %f\n", 
                    partial_fugacities[i][s], n_ig[i][s], results["⟨N⟩ (molecules/unit cell)"][s])
        end
        # test results
        @test isapprox(sum(results["⟨N⟩ (molecules/unit cell)"]), sum(n_ig[i]), rtol=0.05)
        @test isapprox(results["⟨N⟩ (molecules/unit cell)"][1], n_ig[i][1], rtol=0.05)
        @test isapprox(results["⟨N⟩ (molecules/unit cell)"][2], n_ig[i][2], rtol=0.05)
    end
end
