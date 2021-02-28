using PorousMaterials
using Test
using Printf

tests_to_run = Dict("Kr/Xe in SBMOF-1" => false,
                    "ideal gas mixture" => true
                    )
###
#  Simulate Kr/Xe mixture in SBMOF-1
###
if tests_to_run["Kr/Xe in SBMOF-1"]
    @info "running Kr/Xe in SBMOF-1 test"
    # set up sim
    xtal            = Crystal("SBMOF-1.cif")
    molecules       = [Molecule("Kr"), Molecule("Xe")]
    ljff            = LJForceField("UFF", mixing_rules="Lorentz-Berthelot")
    temperature     = 298.0
    n_sample_cycles = 50000
    n_burn_cycles   = 50000
    pressures       = [0.8, 0.2]
    # run sim
    results, molecules = μVT_sim(xtal, molecules, temperature, pressures, ljff, n_burn_cycles=n_burn_cycles, n_sample_cycles=n_sample_cycles)
    # TODO: evaluate results.
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
    ideal_gas1   = Molecule("IG")
    ideal_gas2   = Molecule("IG")
    @assert empty_space.atoms.n == 0
    forcefield   = LJForceField("Dreiding")
    temperature  = 298.0
    mol_fraction = [0.70, 0.30] 
    fugacities    = 10.0 .^ [2.0, 4.0, 5.0, 6.0, 7.0] / 100000.0 # bar
    partial_fugacities = [mol_fraction * f for f in fugacities]

    # according to ideal gas law (P_i*V=n_i*k_b*T), number of molecules per species in box should be
    n_ig  = partial_fugacities * empty_space.box.Ω / (PorousMaterials.KB * temperature) * 100000.0
    n_sim = similar(n_ig)
    for i in 1:length(partial_fugacities)
        results, molecules = μVT_sim(empty_space, [ideal_gas1, ideal_gas2], temperature, partial_fugacities[i], forcefield,
                                     n_burn_cycles=500000, n_sample_cycles=500000, verbose=true)

        for s in 1:2
            @printf("\nfugacities %f, n_ig = %f, n_sim = %f\n", 
                    partial_fugacities[i][s], n_ig[i][s], results["⟨N⟩ (molecules/unit cell)"][s])
        end
        @test isapprox(sum(results["⟨N⟩ (molecules/unit cell)"]), sum(n_ig[i]), rtol=0.05)
        @test isapprox(results["⟨N⟩ (molecules/unit cell)"][1], n_ig[i][1], rtol=0.05)
        @test isapprox(results["⟨N⟩ (molecules/unit cell)"][2], n_ig[i][2], rtol=0.05)
    end
end

