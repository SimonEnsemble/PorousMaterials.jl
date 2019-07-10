module EOS_Test

using PorousMaterials
using OffsetArrays
using LinearAlgebra
using Test
using JLD2
using Statistics
using Random

@testset "EOS Tests" begin
    # Peng-Robinsion EOS test for methane.
    fluid = PengRobinsonFluid(:CH4)
    props = calculate_properties(fluid, 298.0, 65.0, verbose=false)
    @test isapprox(props["compressibility factor"], 0.874496226625811, atol=1e-4)
    @test isapprox(props["fugacity coefficient"], 0.8729028157628362, atol=1e-4)
    @test isapprox(props["fugacity (bar)"], 65.0 * 0.8729028157628362, atol=1e-4)
    @test isapprox(props["density (mol/m³)"], 3000.054418, atol=0.2)
    @test isapprox(props["molar volume (L/mol)"], 0.333327, atol=1e-4)

    # Van der Waals EOS test for CO2 compressibility, molar volume, and density.
    # https://www.webqc.org/van_der_waals_gas_law.html
    fluid = VdWFluid(:CO2_test)
    props = calculate_properties(fluid, 150., 50., verbose=false)
    @test isapprox(props["compressibility factor"], 0.20609717091, atol=1e-4)
    @test isapprox(props["density (mol/m³)"], 19452.3715854, atol=0.2)
    @test isapprox(props["molar volume (L/mol)"], 0.051407613493789, atol=1e-4)

    # Van der Waals EOS test for nitrogen fugacity and fugacity coefficient.
    # https://mabdelsalam.kau.edu.sa/Files/0053615/files/16538_Lecture%207%20phase%20equilibrium.pdf
    fluid = VdWFluid(:N2_test)
    props = calculate_properties(fluid, 298., 50., verbose=false)
    @test isapprox(props["fugacity coefficient"], 0.964, atol = 1e-2)
    @test isapprox(props["fugacity (bar)"], 48.2, atol = 0.1)
end
end
