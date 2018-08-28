module Energetics_Util_Test

using PorousMaterials
using OffsetArrays
using LinearAlgebra
using Test
using JLD2
using Statistics
using Random

@testset "Energetics_Util Tests" begin
    # data types for potential energies
    u = PotentialEnergy(10.0, 30.0)
    v = PotentialEnergy(3.0, 4.0)
    @test ! isapprox(v, PotentialEnergy(3.0, 1.2), verbose=false) # isapprox
    @test isapprox(sum(v), 7.0) # sum
    @test isapprox(u + v, PotentialEnergy(13.0, 34.0)) # +
    @test isapprox(u - v, PotentialEnergy(7.0, 26.0)) # -
    @test isapprox(2.0 * v, PotentialEnergy(6.0, 8.0)) # *
    @test isapprox(v * 2.0, PotentialEnergy(6.0, 8.0)) # *
    @test isapprox(u / 2.0, PotentialEnergy(5.0, 15.0)) # /
    @test isapprox(sqrt(PotentialEnergy(4.0, 16.0)), PotentialEnergy(2.0, 4.0)) # sqrt
    @test isapprox(PorousMaterials.square(PotentialEnergy(2.0, 4.0)), PotentialEnergy(4.0, 16.0)) # square

    t = PotentialEnergy(1.0, 2.0)
    s = PotentialEnergy(300.0, 100.0)
    us = SystemPotentialEnergy(u, v)
    vs = SystemPotentialEnergy(s, t)
    @test isapprox(sum(vs), 403.0) # sum
    @test isapprox(us - vs, SystemPotentialEnergy(u - s, v - t)) # -
    @test isapprox(us + vs, SystemPotentialEnergy(u + s, v + t)) # -
    @test isapprox(2.0 * us, SystemPotentialEnergy(2.0 * u, 2.0 * v)) # *
    @test isapprox(2.0 * us, SystemPotentialEnergy(2.0 * u, 2.0 * v)) # *
    @test isapprox(us * 2.0, SystemPotentialEnergy(2.0 * u, 2.0 * v)) # *
    @test isapprox(us / 2.0, SystemPotentialEnergy(u / 2.0, v / 2.0)) # /
    @test isapprox(sqrt(us), SystemPotentialEnergy(sqrt(u), sqrt(v))) # sqrt
    @test isapprox(PorousMaterials.square(us), SystemPotentialEnergy(PorousMaterials.square(u), PorousMaterials.square(v))) # square
end
end
