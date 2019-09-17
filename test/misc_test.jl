module Misc_Test

using PorousMaterials
using Test
using DataFrames
using CSV
using Optim
using Random

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
    shuffle!(ids_fit)
    df = DataFrame(P = P[ids_fit], N = N[ids_fit])
    x = fit_adsorption_isotherm(df, :P, :N, :langmuir)
    M_opt, K_opt = x["M"], x["K"]
    @test isapprox(M, M_opt, rtol=1e-6)
    @test isapprox(K, K_opt, rtol=1e-6)

    df = CSV.read("Ni-MOF-74_isotherm_test.csv", copycols=true)
    x = fit_adsorption_isotherm(df, Symbol("P(bar)"), Symbol("mol_m3"), :langmuir)
    M_opt, K_opt = x["M"], x["K"]
    @test isapprox(M_opt, 8546.37534030619, atol=1e-5)
    @test isapprox(K_opt, 1.51701035589, atol=1e-5)
end


end
