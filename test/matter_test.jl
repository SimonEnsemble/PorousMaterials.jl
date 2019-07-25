module Matter_Test

using PorousMaterials
using Random
using Test

@testset "Matter Tests" begin
    a1 = Atoms([:a, :b], [1.0 4.0;
                          2.0 5.0;
                          3.0 6.0])
    a2 = Atoms([:c, :d], [7.0 10.0;
                          8.0 11.0;
                          9.0 12.0])
    a3 = a1 + a2
    @test a3.n_atoms == a1.n_atoms + a2.n_atoms
    @test issubset(Set(a1.xf), Set(a3.xf))
    @test issubset(Set(a2.xf), Set(a3.xf))
    @test issubset(Set(a1.species), Set(a3.species))
    @test issubset(Set(a2.species), Set(a3.species))

    c1 = Charges([0.1, 0.2], [1.0 4.0;
                              2.0 5.0;
                              3.0 6.0])
    c2 = Charges([0.3, 0.4], [7.0 10.0;
                              8.0 11.0;
                              9.0 12.0])
    c3 = c1 + c2
    @test c3.n_charges == c1.n_charges + c2.n_charges
    @test issubset(Set(c1.xf), Set(c3.xf))
    @test issubset(Set(c2.xf), Set(c3.xf))
    @test issubset(Set(c1.q), Set(c3.q))
    @test issubset(Set(c2.q), Set(c3.q))

end
end
