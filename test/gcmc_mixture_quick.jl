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

# TODO
###
#  ideal gas mixture tests
###

# To demonstrate how to set up an ideal gas. basically it's a LJ particle with epsilon zero so that it doesn't interact (ideal gas particles). 
# Instead of a MOF, it inserts/deletes/moves particles in an empty box to simulate the idea gas.
# For the mixture, you'd do the same thing but x pv=n_xRT. where x is the mole fraction of the component and n_x is the amount "adsorbed" of that component in the empty box.


#  run GCMC in empty box; should get ideal gas law.
#  "ig" in test_forcefield.csv has sigma tiny and epsilon 0.0 to match ideal gas.
#  basically, this tests the acceptance rules when energy is always zero.

if tests_to_run["ideal_gas"]
    @info "running ideal gas test"
    empty_space = replicate(Crystal("empty_box.cssr"), (3,3,3)) # zero atoms!
    ideal_gas = Molecule("IG")
    @assert empty_space.atoms.n == 0
    forcefield = LJForceField("Dreiding")
    temperature = 298.0
    fugacity = 10.0 .^ [2.0, 4.0, 5.0, 6.0, 7.0, 8.0] / 100000.0 # bar
    # according to ideal gas law, number of molecules in box should be:
    n_ig = fugacity * empty_space.box.Ω / (PorousMaterials.KB * temperature) * 100000.0
    n_sim = similar(n_ig)
    for i = 1:length(fugacity)
        results, molecules = μVT_sim(empty_space, ideal_gas, temperature, fugacity[i], forcefield,
                    n_burn_cycles=100000, n_sample_cycles=100000, verbose=false)
#                    n_burn_cycles=100000, n_sample_cycles=100000)
        @printf("fugacity %f, n_ig = %f, n_sim = %f\n", fugacity[i], n_ig[i], results["⟨N⟩ (molecules/unit cell)"])
        @test isapprox(results["⟨N⟩ (molecules/unit cell)"], n_ig[i], rtol=0.05)
    end
end 

