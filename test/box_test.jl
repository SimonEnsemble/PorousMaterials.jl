module Box_Test

using PorousMaterials
using OffsetArrays
using LinearAlgebra
using Test
using JLD2
using Statistics
using Random

@testset "Box Tests" begin
    framework = Framework("SBMOF-1_cory.cif")
    @test isapprox(framework.box, Box(framework.box.f_to_c))
    @test framework.box.f_to_c * framework.box.c_to_f ≈ Matrix{Float64}(I, 3, 3)
    @test isapprox(framework.box.reciprocal_lattice, 2 * π * inv(framework.box.f_to_c))
    @test isapprox(framework.box, Box(framework.box.a, framework.box.b, framework.box.c,
                                      framework.box.α, framework.box.β, framework.box.γ))
    @test isapprox(replicate(framework.box, (1, 1, 1)), framework.box)
    box = UnitCube()
    @test box.Ω ≈ 1.0
    @test isapprox(replicate(box, (3, 5, 4)), Box(3.0, 5.0, 4.0, π/2, π/2, π/2))
    @test framework.box.Ω ≈ det(framework.box.f_to_c)
    # test alternative constructor using f_to_c matrix
    @test isapprox(framework.box, Box(framework.box.f_to_c))
    
    # inside box function
    box = Box(1.0, 20.0, 30.0, π/2, π/2, π/2)
    @test inside([0.001, 0.0001*20, 0.0001*30], box)
    @test inside([0.99, 0.99*20, 0.99*30], box)
    @test ! inside([-0.01, 0.99*20, 0.99*30], box)
    @test ! inside([0.01, 1.01*20, 0.99*30], box)
    @test inside(framework.box.f_to_c * [0.1, 0.9, 0.7], framework.box)
    @test ! inside(framework.box.f_to_c * [0.1, -0.9, 0.7], framework.box)
    @test ! inside(framework.box.f_to_c * [0.1, 1.2, 0.7], framework.box)
end
end
