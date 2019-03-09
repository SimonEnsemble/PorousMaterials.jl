module Generic_Rotation_Test

using PorousMaterials
using Test
using LinearAlgebra

@testset "Generic Rotation Tests" begin
    # direction and angle about which to rotate
    u = randn(3)
    u = u / norm(u)
    θ = π / 7

    # test unit vector normalization
    @test isapprox(rotation_matrix(θ, u), rotation_matrix(θ, u * 10))
    R = rotation_matrix(θ, u)
    @test isapprox(transpose(R) * R, diagm(0 => ones(3)))
    # rotate θ, then -θ, will get back...
    @test isapprox(rotation_matrix(-θ, u) * rotation_matrix(θ, u * 10), diagm(0 => ones(3)))
    
    # test rotating about x, y, z axes
    @test isapprox(rotation_matrix(θ, [1.0, 0.0, 0.0]), rotation_matrix(θ, 1))
    @test isapprox(rotation_matrix(θ, [0.0, 1.0, 0.0]), rotation_matrix(θ, 2))
    @test isapprox(rotation_matrix(θ, [0.0, 0.0, 1.0]), rotation_matrix(θ, 3))
    
    # determinant should be 1.0
    @test isapprox(det(rotation_matrix(rand() * 2 * π, randn(3))), 1.0)
end
end
