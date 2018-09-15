module Path_Test

using PorousMaterials
using OffsetArrays
using LinearAlgebra
using Test
using JLD2
using Statistics
using Random

@testset "Path Tests" begin
    set_tutorial_mode()
    @test PorousMaterials.PATH_TO_DATA == joinpath(dirname(pathof(PorousMaterials)),"..","test","data")
    set_path_to_data(pwd())
    @test PorousMaterials.PATH_TO_DATA == pwd()
    set_path_to_data()
    @test PorousMaterials.PATH_TO_DATA == joinpath(pwd(),"data")
end
end
