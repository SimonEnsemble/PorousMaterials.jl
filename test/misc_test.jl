module Misc_Test

using PorousMaterials
using Test
using DataFrames
using CSV
using MultivariateStats

@testset "Misc Tests" begin
    df = CSV.read("isotherm.csv")
    p_col_name, l_col_name = names(df)
    henry = extract_henry_coefficient(df, p_col_name, l_col_name, 4)
    @test isapprox(henry, 0.9404942451976929)
end

end
