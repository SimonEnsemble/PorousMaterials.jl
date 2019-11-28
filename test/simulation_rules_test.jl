module Simulation_Rules_Test

using PorousMaterials
using Test

# Test set for making sure simulations meet certain rules
#   i.e. frameworks musts be in P1 symmetry to be used in GCMC or Henry
#   coefficient test
@testset "Simulation Rules Tests" begin
    non_P1_framework = Framework("symmetry_test_structure.cif", convert_to_p1=false)
    molecule = Molecule("CO2")
    ljff = LJForceField("UFF.csv")
    temp = 298.0
    pressure = 5.0

    # make sure that the framework is not in P1 before running tests

    # assertion errors should be thrown when trying to run gcmc simulations
    #   with non-P1 frameworks
    @test_throws AssertionError gcmc_simulation(non_P1_framework, molecule,
                                               temp, pressure, ljff)
    pressures = [1.0, 5.0, 10.0, 15.0, 20.0]
    @test_throws AssertionError adsorption_isotherm(non_P1_framework, molecule,
                                                   temp, pressures, ljff)
    @test_throws AssertionError stepwise_adsorption_isotherm(non_P1_framework,
                                            molecule, temp, pressures, ljff)

    # assertion error should be thrown when attempting to run a henry
    #   coefficient with a non-P1 framework
    @test_throws AssertionError henry_coefficient(non_P1_framework, molecule,
                                                 temp, ljff)

    # Test that an assertion is thrown when a non-P1 structure is passed into
    #   energy_grid(), and compute_accessibility_grid()
    @test_throws AssertionError energy_grid(non_P1_framework, molecule, ljff)
    @test_throws AssertionError compute_accessibility_grid(non_P1_framework, molecule, ljff)

    # Test that an assertion is thrown when trying to replicate a non-P1
    #   structure
    @test_throws AssertionError replicate(non_P1_framework, (2, 2, 2))
    @test_throws AssertionError replicate(non_P1_framework, (1, 1, 1))

end
end
