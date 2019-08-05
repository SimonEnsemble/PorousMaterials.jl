module Simulation_Rules_Test

using PorousMaterials
using Test

# Test set for making sure simulations meet certain rules
#   i.e. frameworks musts be in P1 symmetry to be used in GCMC or Henry
#   coefficient test
@testset "Simulation Rules Tests" begin
    non_P1_framework = Framework("ORIVOC_clean.cif", convert_to_p1=false)
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
end
end
