module Nearest_Image_Test

using PorousMaterials
using OffsetArrays
using LinearAlgebra
using Test
using JLD2
using Statistics
using Random

@testset "NearestImage Tests" begin
    dxf = [-0.8, -0.4, 0.7]
    nearest_image!(dxf)
    @test isapprox(dxf, [0.2, -0.4, -0.3])

    dxf = [-0.3, -0.1, -0.9]
    nearest_image!(dxf)
    @test isapprox(dxf, [-0.3, -0.1, 0.1])

    box = UnitCube()
    x = [0.9, 0.1, 0.1]
    y = [0.0, 0.0, 0.0]
    @test nearest_r(x, y, box) ≈ norm([-0.1, 0.1, 0.1])
    x = [0.5, 0.5, 0.5]
    y = [0.4, 0.4, 0.35]
    @test nearest_r(x, y, box) ≈ norm(x-y)
end
end
