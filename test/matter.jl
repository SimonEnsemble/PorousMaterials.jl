module Matter_Test

using PorousMaterials
using Test

@testset "Matter Tests" begin
    f1 = Frac([1.0 4.0;
               2.0 5.0;
               3.0 6.0]
              )
    f2 = Frac([7.0 10.0;
               8.0 11.0;
               9.0 12.0]
              )

    c1 = Cart([1.0 4.0;
               2.0 5.0;
               3.0 6.0]
              )
    c2 = Cart([7.0 10.0;
               8.0 11.0;
               9.0 12.0]
              )

    s1 = [:a, :b]
    s2 = [:c, :d]

    ## Atoms
    a1_f = Atoms(s1, f1)
    a2_f = Atoms(s2, f2)

    a1_c = Atoms(s1, c1)
    a2_c = Atoms(s2, c2)

    @test ! isapprox(a1_c, a2_c)
    @test ! isapprox(a1_f, a2_f)
    @test ! isapprox(Atoms(s1, f1), Atoms(s1, f2))
    @test ! isapprox(Atoms(s1, f1), Atoms(s2, f1))
    @test ! isapprox(Atoms(s1, c1), Atoms(s1, c2))
    @test ! isapprox(Atoms(s1, c1), Atoms(s2, c1))

    a3_c = a1_c + a2_c
    @test a3_c.n == a1_c.n + a2_c.n
    @test a3_c.species[1:2] == s1
    @test a3_c.species[3:4] == s2
    @test a3_c.coords.x[:, 1:2] == c1.x
    @test a3_c.coords.x[:, 3:4] == c2.x
    
    a3_f = a1_f + a2_f
    @test a3_f.n == a1_f.n + a2_f.n
    @test a3_f.species[1:2] == s1
    @test a3_f.species[3:4] == s2
    @test a3_f.coords.xf[:, 1:2] == f1.xf
    @test a3_f.coords.xf[:, 3:4] == f2.xf

    ## Charges
    q_vals_1 = [1.0, 4.0]
    q_vals_2 = [0.4, 0.3]

    q1_c = Charges(q_vals_1, c1)
    q2_c = Charges(q_vals_2, c2)
    @test ! neutral(q1_c)

    q1_f = Charges(q_vals_1, f1)
    q2_f = Charges(q_vals_2, f2)
    
    @test ! isapprox(q1_c, q2_c)
    @test ! isapprox(q1_f, q2_f)
    @test ! isapprox(Charges(q_vals_1, f1), Charges(q_vals_1, f2))
    @test ! isapprox(Charges(q_vals_1, f1), Charges(q_vals_2, f1))
    @test ! isapprox(Charges(q_vals_1, c1), Charges(q_vals_1, c2))
    @test ! isapprox(Charges(q_vals_1, c1), Charges(q_vals_2, c1))

    q3_c = q1_c + q2_c
    @test q3_c.n == q1_c.n + q2_c.n
    @test q3_c.q[1:2] == q_vals_1
    @test q3_c.q[3:4] == q_vals_2
    @test q3_c.coords.x[:, 1:2] == c1.x
    @test q3_c.coords.x[:, 3:4] == c2.x
    
    q3_f = q1_f + q2_f
    @test q3_f.n == q1_f.n + q2_f.n
    @test q3_f.q[1:2] == q_vals_1
    @test q3_f.q[3:4] == q_vals_2
    @test q3_f.coords.xf[:, 1:2] == f1.xf
    @test q3_f.coords.xf[:, 3:4] == f2.xf

    # wrap
    f = Frac([0.1 -0.1;
              2.2 0.56;
             -0.4 1.2]
             )
    fw = Frac([0.1 0.9;
               0.2 0.56;
               0.6 0.2]
             )
    wrap!(f)
    @test isapprox(f, fw)
    
    q = Charges([-0.1, 0.2, -0.1], Frac(zeros(3, 3)))
    @test neutral(q)
    @test isapprox(net_charge(q), 0.0)
    q.q[2] = 0.0
    @test ! neutral(q)
    @test isapprox(net_charge(q), -0.2)
    q = Charges{Frac}(0)
    @test net_charge(q) == 0.0
    q = Charges{Cart}(0)
    @test net_charge(q) == 0.0
    
    # safe constructors
    atoms = Atoms{Frac}(10)
    @test all(isnan.(atoms.coords.xf))
    @test all(atoms.species .== :_)
    atoms = Atoms{Cart}(10)
    @test all(isnan.(atoms.coords.x))
    @test all(atoms.species .== :_)

    charges = Charges{Frac}(10)
    @test all(isnan.(charges.coords.xf))
    @test all(isnan.(charges.q))
    charges = Charges{Cart}(10)
    @test all(isnan.(charges.coords.x))
    @test all(isnan.(charges.q))

    # set equal
    f = rand(3, 6)
    s = [:a, :b, :c, :d, :e, :f]
    q = rand(6)
    ids = [3, 6, 5, 4, 2, 1]

    a = Atoms(s, Frac(f))
    as = Atoms(s[ids], Frac(f[:, ids]))
    @test equal_sets(a, as)
    @test ! isapprox(a, as)
    a.coords.xf[2, 3] += 0.005
    @test ! equal_sets(a, as)
    @test equal_sets(a, as, digits=1)

    c = Charges(q, Frac(f))
    cs = Charges(q[ids], Frac(f[:, ids]))
    @test equal_sets(c, cs)
    @test ! isapprox(c, cs)
    c.coords.xf[2, 3] += 0.005
    @test ! equal_sets(c, cs)
    @test equal_sets(c, cs, digits=1)
    
    # getindex (array indexing!)
    atoms = Crystal("SBMOF-1.cif").atoms
    @test isapprox(atoms[3:4].coords, atoms.coords[3:4])
    @test isapprox(atoms[3:4].coords.xf, atoms.coords.xf[:, 3:4])
    @test atoms[3:4].species == atoms.species[3:4]
    @test length(atoms[6:10].species) == 5 == atoms[6:10].n

    charges = Crystal("ATIBOU01_clean.cif").charges
    @test isapprox(charges[3:4].coords, charges.coords[3:4])
    @test isapprox(charges[3:4].coords.xf, charges.coords.xf[:, 3:4])
    @test isapprox(charges[3:4].q, charges.q[3:4])
    @test length(charges[6:10].q) == 5 == charges[6:10].n
end
end
