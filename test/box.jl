module Box_Test

using PorousMaterials
using LinearAlgebra
using Test

@testset "Box Tests" begin
    box = Box(11.6, 5.5, 22.9, 90.0 * π / 180, 100.8 * π / 180.0, 90.0 * π / 180.0)
    @test isapprox(box, Box(box.f_to_c)) # alt. constructor
    @test isapprox(box, Box(box.a, box.b, box.c,
                            box.α, box.β, box.γ)) # alt. constructor
    @test box.f_to_c * box.c_to_f ≈ Matrix{Float64}(I, 3, 3)
    @test isapprox(box.reciprocal_lattice, 2 * π * inv(box.f_to_c)) #TODO just use this method to construct recip. lattice.
    @test isapprox(replicate(box, (1, 1, 1)), box)
    @test box.Ω ≈ det(box.f_to_c)
    
    # converters
    box = Box(11.6, 5.5, 22.9, 90.0 * π / 180, 100.8 * π / 180.0, 92.0 * π / 180.0)
    f = Frac([0.05, 0.9, 0.99])
    c = Cart(f, box)
    f_conv = Frac(c, box)
    @test isapprox(c.x, box.f_to_c * f.xf)
    @test isapprox(f, f_conv)

    box = unit_cube()
    @test isapprox(Box(1.0, 1.0, 1.0), unit_cube()) # alt. constructor, right angles assumed
    @test box.Ω ≈ 1.0
    @test isapprox(replicate(box, (3, 5, 4)), Box(3.0, 5.0, 4.0))
    
    box = Box(10.0, 10.0, 10.0)
    f = Frac([0.1 0.1;
              0.3 0.2;
              0.5 0.6]
             )
    atoms = Atoms([:C, :O], f)
    atoms_c = Cart(atoms, box)
    @test atoms_c.species == atoms.species
    @test isapprox(atoms_c.coords.x, box.f_to_c * atoms.coords.xf)
    @test isapprox(atoms, Frac(atoms_c, box))
    charges = Charges([1.0, -1.0], f)
    charges_c = Cart(charges, box)
    @test isapprox(charges_c.q, charges.q)
    @test isapprox(charges_c.coords.x, box.f_to_c * charges.coords.xf)
    @test isapprox(charges, Frac(charges_c, box))
    
    # inside box function
    box = Box(1.0, 20.0, 30.0)
    f = Frac([0.1 0.6;
               0.9 0.8;
               0.8 0.9]
              )
    @test inside(f)
    f.xf[2, 2] = -0.1
    @test ! inside(f)
    f.xf[2, 2] = 1.2
    @test ! inside(f)

    c = Cart([0.1 0.4;
              4.5 13.0;
              22. 10.1]
             )
    @test inside(c, box)
    c.x[2, 2] = -0.1
    @test !inside(c, box)
    c.x[2, 2] = 20.1
    @test !inside(c, box)

    box = Box(11.6, 5.5, 22.9, 90.0 * π / 180, 100.8 * π / 180.0, 90.0 * π / 180.0)
    f = Frac([0.1, 0.8, 0.9])
    @test inside(f)
    @test inside(Cart(f, box), box)

    # translate_by!
    box = Box(1.0, 10.0, 40.0)
    f = Frac([0.1 0.2;
              0.2 0.5;
              0.3 0.8]
            )
    dx = Cart([0.2, 1.0, 4.0])
    translate_by!(f, dx, box)
    @test isapprox(f, Frac([0.3 0.4;
                             0.3 0.6;
                             0.4 0.9]
                            )
                  )

    c = Cart([1.0 0.0;
              10.0 0.0;
              40.0 0.0]) #corner of box
    dxf = Frac([-1.0, -1.0, -1.0])
    translate_by!(c, dxf, box)
    @test isapprox(c, Cart([0.0 -1.0;
                            0.0 -10.0;
                            0.0 -40.0])
                  )
end
end
