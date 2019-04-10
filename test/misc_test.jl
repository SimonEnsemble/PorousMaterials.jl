module Misc_Test

using PorousMaterials
using Test
using DataFrames
using CSV
using Optim

@testset "Misc Tests" begin
    df = CSV.read("isotherm.csv")
    p_col_name, l_col_name = names(df)
    henry = fit_isotherm(df[1:3,:], p_col_name, l_col_name, :henry)["K"]
    @test isapprox(henry, 0.9688044280548711)

    P = range(0, stop=1, length=100)
    M = 23.0
    K = 11.9
    N = (M*K.*P)./(1 .+ K.*P)
    N = [N[1], N[11], N[21], N[31], N[41], N[51], N[61], N[71], N[81]]
    P = [P[1], P[11], P[21], P[31], P[41], P[51], P[61], P[71], P[81]]
    new_df = DataFrame(Pressure = P, Loading = N)
    x = fit_isotherm(new_df, Symbol("Pressure"), Symbol("Loading"), :langmuir)
    M2, K2 = x["M"], x["K"]
    @test isapprox(M, M2)
    @test isapprox(K, K2)
end


end
