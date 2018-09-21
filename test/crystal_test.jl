module Crystal_Test

using PorousMaterials
using OffsetArrays
using LinearAlgebra
using Test
using JLD2
using Statistics
using Random

@testset "Crystal Tests" begin
    framework = Framework("test_structure2.cif")
    strip_numbers_from_atom_labels!(framework)
    @test framework.name == "test_structure2.cif"
    @test isapprox(framework.box, Box(10.0, 20.0, 30.0, 90*π/180, 45*π/180, 120*π/180))
    @test framework.atoms.n_atoms == 2
    @test isapprox(framework.atoms, Atoms([:Ca, :O], [0.2 0.6; 0.5 0.3; 0.7 0.1]))
    @test isapprox(framework.charges, Charges([1.0, -1.0], [0.2 0.6; 0.5 0.3; 0.7 0.1]))
    new_frame = assign_charges(framework, Dict(:Ca => -2.0, :O => 2.0))
    @test isapprox(new_frame.charges, Charges([-2.0, 2.0], [0.2 0.6; 0.5 0.3; 0.7 0.1]))
    new_frame = assign_charges(framework, [4.0, -4.0])
    @test isapprox(new_frame.charges, Charges([4.0, -4.0], [0.2 0.6; 0.5 0.3; 0.7 0.1]))
    @test charged(framework)
    @test chemical_formula(framework) == Dict(:Ca => 1, :O => 1)
    @test molecular_weight(framework) ≈ 15.9994 + 40.078
    # same as test_structure.cif but with overlapping atoms.
    framework2 = Framework("test_structure2B.cif", remove_overlap=true, check_charge_neutrality=false)
    strip_numbers_from_atom_labels!(framework2)
    @test isapprox(framework.atoms, framework2.atoms) && isapprox(framework.charges, framework2.charges)

    # test .cif writer; write, read in, assert equal
    write_cif(framework, joinpath("data", "crystals", "rewritten_test_structure2.cif"))
    framework_rewritten = Framework("rewritten_test_structure2.cif")
    @test isapprox(framework, framework_rewritten)

    # test .cssr reader too; test_structure2.{cif,cssr} designed to be the same.
    framework_from_cssr = Framework("test_structure2.cssr")
    strip_numbers_from_atom_labels!(framework_from_cssr)
    @test isapprox(framework_from_cssr, framework, checknames=false)

    # test replicate framework
    sbmof = Framework("SBMOF-1.cif")
    replicated_sbmof = replicate(sbmof, (1, 1, 1))
    @test isapprox(sbmof, replicated_sbmof)

    repfactors = replication_factors(sbmof.box, 14.0)
    replicated_sbmof = replicate(sbmof, repfactors)
    @test replication_factors(replicated_sbmof.box, 14.0) == (1, 1, 1)
    @test isapprox(sbmof.atoms.xf[:, 1] ./ repfactors, replicated_sbmof.atoms.xf[:, 1])
    @test isapprox(replicated_sbmof.box.reciprocal_lattice, 2 * π * inv(replicated_sbmof.box.f_to_c))
    @test chemical_formula(sbmof) == chemical_formula(replicated_sbmof)
    @test isapprox(crystal_density(sbmof), crystal_density(replicated_sbmof), atol=1e-7)

    # more xtal tests
    sbmof1 = Framework("SBMOF-1.cif")
    @test !charged(sbmof1)
    @test isapprox(sbmof1.box.reciprocal_lattice, 2 * π * inv(sbmof1.box.f_to_c))
    @test sbmof1.box.Ω ≈ det(sbmof1.box.f_to_c) # sneak in crystal test
    @test isapprox(crystal_density(sbmof1), 1570.4, atol=0.5) # kg/m3

    # replicating the unit cell to construct simulation box
    sbmof1 = Framework("SBMOF-1.cif")
    rbox = replicate(sbmof1.box, (2, 3, 4))
    @test rbox.Ω ≈ sbmof1.box.Ω * 2 * 3 * 4
    @test all(rbox.c_to_f * sbmof1.box.f_to_c * [1.0, 1.0, 1.0] .≈ [1/2, 1/3, 1/4])

#=    ## .cssr reader test
    # test replicate framework
    sbmof = Framework("SBMOF-1.cif")
    replicated_sbmof = replicate(sbmof, (1, 1, 1))
    @test isapprox(sbmof, replicated_sbmof)

    repfactors = replication_factors(sbmof.box, 14.0)
    replicated_sbmof = replicate(sbmof, repfactors)
    @test replication_factors(replicated_sbmof.box, 14.0) == (1, 1, 1)
    @test isapprox(sbmof.atoms.xf[:, 1] ./ repfactors, replicated_sbmof.atoms.xf[:, 1])
    @test isapprox(replicated_sbmof.box.reciprocal_lattice, 2 * π * inv(replicated_sbmof.box.f_to_c))
    @test chemical_formula(sbmof) == chemical_formula(replicated_sbmof)
    @test isapprox(crystal_density(sbmof), crystal_density(replicated_sbmof), atol=1e-7)

    # more xtal tests
    sbmof1 = Framework("SBMOF-1.cif")
    @test !charged(sbmof1)
    @test isapprox(sbmof1.box.reciprocal_lattice, 2 * π * inv(sbmof1.box.f_to_c))
    @test sbmof1.box.Ω ≈ det(sbmof1.box.f_to_c) # sneak in crystal test
    @test isapprox(crystal_density(sbmof1), 1570.4, atol=0.5) # kg/m3

    # replicating the unit cell to construct simulation box
    sbmof1 = Framework("SBMOF-1.cif")
    rbox = replicate(sbmof1.box, (2, 3, 4))
    @test rbox.Ω ≈ sbmof1.box.Ω * 2 * 3 * 4
    @test all(rbox.c_to_f * sbmof1.box.f_to_c * [1.0, 1.0, 1.0] .≈ [1/2, 1/3, 1/4])=#
end
end
