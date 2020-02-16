module Distance_Test

using PorousMaterials
using Test
using LinearAlgebra

@testset "Distance Tests" begin
    dxf = [-0.8, -0.4, 0.7]
    nearest_image!(dxf)
    @test isapprox(dxf, [0.2, -0.4, -0.3])

    dxf = [-0.3, -0.1, -0.9]
    nearest_image!(dxf)
    @test isapprox(dxf, [-0.3, -0.1, 0.1])
    
    # distance (fractional)
    f = Frac([0.2 0.4;
              0.1 0.8;
              0.8 0.6]
              )
    box = unit_cube()
    @test isapprox(distance(f, box, 1, 2, false),  norm(f.xf[:, 1] - f.xf[:, 2]))
    @test isapprox(distance(f, box, 2, 2, false), 0.0)
    @test isapprox(distance(f, box, 2, 1, false), distance(f, box, 1, 2, false))
    @test isapprox(distance(f, box, 1, 2, true),  norm(f.xf[:, 1] - [0.4, -0.2, 0.6]))
    box = Box(1.0, 10.0, 100.0)
    @test isapprox(distance(f, box, 1, 2, false),  norm([0.2, 1.0, 80.0] - [0.4, 8.0, 60.0]))
    @test isapprox(distance(f, box, 2, 1, true),  norm([0.2, 1.0, 80.0] - [0.4, -2.0, 60.0]))
    
    # distance (Cartesian)
    c = Cart([0.2 0.4;
              0.1 0.8;
              0.8 0.6]
              )
    box = unit_cube()
    @test isapprox(distance(c, box, 1, 2, false),  norm(c.x[:, 1] - c.x[:, 2]))
    @test isapprox(distance(c, box, 2, 2, false), 0.0)
    @test isapprox(distance(c, box, 2, 1, false), distance(c, box, 1, 2, false))
    @test isapprox(distance(c, box, 1, 2, true),  norm(c.x[:, 1] - [0.4, -0.2, 0.6]))
    box = Box(10.0, 1.0, 100.0)
    @test isapprox(distance(c, box, 1, 2, false), norm(c.x[:, 1] - c.x[:, 2]))
    @test isapprox(distance(c, box, 2, 1, true),  norm([0.2, 0.1, 0.8] - [0.4, -0.2, 0.6]))

    # distance (Tests from Avogadro measurement tool)
    crystal = Crystal("distance_test.cif")
    @test isapprox(distance(crystal.atoms, crystal.box, 1, 3, false), 12.334, atol=0.001) # C- Ca
    @test isapprox(distance(crystal.atoms, crystal.box, 1, 2, false), 14.355, atol=0.001) # C-S
    @test isapprox(distance(crystal.atoms, crystal.box, 1, 2, true), 8.841, atol=0.001) # C-S
    @test isapprox(distance(crystal.atoms, crystal.box, 2, 3, false), 6.292, atol=0.001) # S-Ca

    # overlap
    box = unit_cube()
    f = Frac([0.2 0.4 0.6 0.401;
              0.1 0.8 0.7 0.799;
              0.8 0.6 0.5 0.602]
              )
    o_flag, o_ids = overlap(f, box, true, tol=0.01)
    @test o_flag
    @test o_ids == [(2, 4)]
    
    f = Frac([0.2 0.4 0.6 0.2;
              0.1 0.8 0.7 0.1;
              0.99 0.6 0.5 0.01]
              )
    o_flag, o_ids = overlap(f, box, true, tol=0.03)
    @test o_flag
    @test o_ids == [(1, 4)]
    o_flag, o_ids = overlap(f, box, false, tol=0.03)
    @test !o_flag
    @test o_ids == []
    
    # test distance function (via Avogadro)
    crystal = Crystal("simple_test.cif")
    @test distance(crystal.atoms, crystal.box, 1, 1, true) == 0.0
    @test isapprox(distance(crystal.atoms, crystal.box, 2, 5, true), 4.059, atol=0.001)
    @test isapprox(distance(crystal.atoms, crystal.box, 2, 5, false), 4.059, atol=0.001)
    @test isapprox(distance(crystal.atoms, crystal.box, 1, 5, false), 17.279, atol=0.001)
    @test isapprox(distance(crystal.atoms, crystal.box, 1, 5, true), 1.531, atol=0.001)

end
end
