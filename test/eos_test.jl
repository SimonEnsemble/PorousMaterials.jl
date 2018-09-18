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
    gas = PengRobinsonFluid(:CH4)
    props = calculate_properties(gas, 298.0, 65.0, verbose=false)
    @test isapprox(props["compressibility factor"], 0.874496226625811, atol=1e-4)
    @test isapprox(props["fugacity coefficient"], 0.8729028157628362, atol=1e-4)
    @test isapprox(props["fugacity (bar)"], 65.0 * 0.8729028157628362, atol=1e-4)
    @test isapprox(props["density (mol/m³)"], 3000.054418, atol=0.2)
    @test isapprox(props["molar volume (L/mol)"], 0.333327, atol=1e-4)

    #Van der Waals EOS test for Hydrogen
    gas = VDWFluid(:H2)
    props = calculate_properties(gas, 300.0, 65.0)
    @test isapprox(props["Compressibility Factor"], 1.0462045 , atol=0.001)
    @test isapprox(props["Fugacity Coefficient"], 67.98203133 / 65.0 , atol=1e-4)
    @test isapprox(props["Fugacity (bar)"], 67.98203133 , atol=1e-4)
    @test isapprox(props["Density (mol/m³)"], 2490.815 , atol=0.2)
    @test isapprox(props["Molar Volume (L/mol)"], 0.401475 , atol=1e-4)
end
end
