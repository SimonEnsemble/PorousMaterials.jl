module Matter_Test

using PorousMaterials
using Random
using Test

@testset "Matter Tests" begin
    # testing addition for atoms type
    a1 = Atoms([:a, :b], [1.0 4.0;
                          2.0 5.0;
                          3.0 6.0])
    a2 = Atoms([:c, :d], [7.0 10.0;
                          8.0 11.0;
                          9.0 12.0])
    a3 = a1 + a2
    # total atoms in addition is number of atoms combined
    @test a3.n_atoms == a1.n_atoms + a2.n_atoms
    # atoms from a1 & a2 should be within the atoms for a3
    @test issubset(Set(a1.xf), Set(a3.xf))
    @test issubset(Set(a2.xf), Set(a3.xf))
    @test issubset(Set(a1.species), Set(a3.species))
    @test issubset(Set(a2.species), Set(a3.species))
    # no atoms that weren't in a1 or a2 are in a3
    @test issetequal(setdiff(Set(a3.xf), union(Set(a1.xf), Set(a2.xf))), Set())
    @test issetequal(setdiff(Set(a3.species), union(Set(a1.species), Set(a2.species))), Set())
    a1_mismatch = Atoms([:b, :a], [1.0 4.0;
                                   2.0 5.0;
                                   3.0 6.0])
    @test ! has_same_set_of_atoms(a1, a1_mismatch)

    # testing addition for charges type
    c1 = Charges([0.1, 0.2], [1.0 4.0;
                              2.0 5.0;
                              3.0 6.0])
    c2 = Charges([0.3, 0.4], [7.0 10.0;
                              8.0 11.0;
                              9.0 12.0])
    c3 = c1 + c2
    # total charges in c3 must be the number of charges in c1 and c2
    @test c3.n_charges == c1.n_charges + c2.n_charges
    # charges from c1 & c2 should be within c3
    @test issubset(Set(c1.xf), Set(c3.xf))
    @test issubset(Set(c2.xf), Set(c3.xf))
    @test issubset(Set(c1.q), Set(c3.q))
    @test issubset(Set(c2.q), Set(c3.q))
    # no charges are added to c3 that weren't within c1 or c2
    @test issetequal(setdiff(Set(c3.xf), union(Set(c1.xf), Set(c2.xf))), Set())
    @test issetequal(setdiff(Set(c3.q), union(Set(c1.q), Set(c2.q))), Set())
    c1_mismatch = Charges([0.2, 0.1], [1.0 4.0;
                                       2.0 5.0;
                                       3.0 6.0])
    @test ! has_same_set_of_charges(c1, c1_mismatch)

end
end
