module Misc_Test

using PorousMaterials
using Test
using DataFrames
using CSV
using Optim

@testset "Misc Tests" begin
    ###
    # adsorption isotherm tests
    ###
    # Henry
    df = DataFrame(P = [0.0, 0.56, 1.333], N = [0.0, 0.534, 1.295])
    henry = fit_adsorption_isotherm(df, :P, :N, :henry)["H"]
    @test isapprox(henry, 0.9688044280548711)
    
    # Langmuir
    P = range(0, stop=1, length=100)
    M = 23.0
    K = 11.9
    N = (M * K .* P) ./ (1 .+ K.*P)
    ids_fit = [1, 6, 11, 21, 31, 36, 41, 51, 61, 71, 76, 81, 91]
    df = DataFrame(P = P[ids_fit], N = N[ids_fit])
    x = fit_adsorption_isotherm(df, :P, :N, :langmuir)
    M2, K2 = x["M"], x["K"]
    @test isapprox(M, M2, rtol=1e-6)
    @test isapprox(K, K2, rtol=1e-6)
    df = CSV.read("Ni-MOF-74_isotherm_test.csv", copycols=true)
    x = fit_adsorption_isotherm(df, Symbol("P(bar)"), Symbol("mol_m3"), :langmuir)
    M3, K3 = x["M"], x["K"]
    @test isapprox(M3, 8546.37534030619)
    @test isapprox(K3, 1.51701035589)
end


end
