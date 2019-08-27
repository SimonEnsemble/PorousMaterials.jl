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
    #   more write_cif tests below with symmetry tests
    write_cif(framework, joinpath("data", "crystals", "rewritten_test_structure2.cif"))
    framework_rewritten = Framework("rewritten_test_structure2.cif")
    write_cif(framework, joinpath("data", "crystals", "rewritten_test_structure2_cartn.cif"); fractional=false)
    framework_rewritten_cartn = Framework("rewritten_test_structure2_cartn.cif")
    @test isapprox(framework, framework_rewritten)
    @test isapprox(framework, framework_rewritten_cartn)

    # test .cif reader for non-P1 symmetry
    #   no atoms should overlap
    #   should place atoms in the same positions as the P1 conversion using
    #       openBabel
    non_P1_framework = Framework("ORIVOC_clean_fract.cif", remove_overlap=true)
    non_P1_cartesian = Framework("ORIVOC_clean.cif", remove_overlap=true)
    P1_framework = Framework("ORIVOC_clean_P1.cif", remove_overlap=true)

    # wrap all atoms and charges to be within the unit cell
    non_P1_framework.atoms.xf .= mod.(non_P1_framework.atoms.xf, 1.0)
    non_P1_framework.charges.xf .= mod.(non_P1_framework.charges.xf, 1.0)

    non_P1_cartesian.atoms.xf .= mod.(non_P1_cartesian.atoms.xf, 1.0)
    non_P1_cartesian.charges.xf .= mod.(non_P1_cartesian.charges.xf, 1.0)

    P1_framework.atoms.xf .= mod.(P1_framework.atoms.xf, 1.0)
    P1_framework.charges.xf .= mod.(P1_framework.charges.xf, 1.0)

    @test isapprox(non_P1_framework, P1_framework; atol=1e-2) 
    # test that fractional and cartesian produce same results
    @test isapprox(non_P1_framework, non_P1_cartesian; atol=1e-2)
    # test that cartesian and P1 produce same results
    @test isapprox(non_P1_cartesian, P1_framework; atol=1e-2)
    # test that incorrect file formats throw proper errors
    @test_throws ErrorException Framework("non_P1_no_symmetry.cif")
    # test that a file with no atoms throws error
    @test_throws ErrorException Framework("no_atoms.cif")
    # test that these carry the is_p1 flag
    @test non_P1_framework.is_p1
    @test non_P1_cartesian.is_p1
    @test P1_framework.is_p1

    # test reading in non-P1 then applying symmetry later
    # read in the same files as above, then convert to P1, then compare
    non_P1_framework_symmetry = Framework("ORIVOC_clean_fract.cif", convert_to_p1=false)
    non_P1_cartesian_symmetry = Framework("ORIVOC_clean.cif", convert_to_p1=false)

    # make sure these frameworks are not in P1 symmetry when convert_to_p1 is
    #   set to false
    @test ! non_P1_framework_symmetry.is_p1
    @test ! non_P1_cartesian_symmetry.is_p1

    # test write_cif in non_p1 symmetry
    write_cif(non_P1_framework_symmetry, joinpath("data", "crystals", "rewritten_ORIVOC_clean_fract.cif"))
    # keep this in cartesian to test both
    write_cif(non_P1_cartesian_symmetry, joinpath("data", "crystals", "rewritten_ORIVOC_clean.cif"), fractional=false)
    rewritten_non_p1_fractional = Framework("rewritten_ORIVOC_clean_fract.cif"; convert_to_p1=false)
    rewritten_non_p1_cartesian = Framework("rewritten_ORIVOC_clean.cif"; convert_to_p1=false)

    @test isapprox(rewritten_non_p1_fractional, non_P1_framework_symmetry)
    @test isapprox(rewritten_non_p1_cartesian, non_P1_cartesian_symmetry)

    non_P1_framework_symmetry = apply_symmetry_rules(non_P1_framework_symmetry, remove_overlap=true)
    non_P1_cartesian_symmetry = apply_symmetry_rules(non_P1_cartesian_symmetry, remove_overlap=true)

    # wrap all atoms and charges to be within the unit cell
    non_P1_framework_symmetry.atoms.xf .= mod.(non_P1_framework_symmetry.atoms.xf, 1.0)
    non_P1_framework_symmetry.charges.xf .= mod.(non_P1_framework_symmetry.charges.xf, 1.0)

    non_P1_cartesian_symmetry.atoms.xf .= mod.(non_P1_cartesian_symmetry.atoms.xf, 1.0)
    non_P1_cartesian_symmetry.charges.xf .= mod.(non_P1_cartesian_symmetry.charges.xf, 1.0)

    # make sure frameworks are now recorded as being in P1 with their is_p1 flag
    @test non_P1_framework_symmetry.is_p1
    @test non_P1_cartesian_symmetry.is_p1

    # test that same structure is created when reading and converting to P1 and
    #   when reading then converting to P1
    @test isapprox(non_P1_framework, non_P1_framework_symmetry)
    @test isapprox(non_P1_cartesian, non_P1_cartesian_symmetry)

    # test .cssr reader too; test_structure2.{cif,cssr} designed to be the same.
    framework_from_cssr = Framework("test_structure2.cssr")
    strip_numbers_from_atom_labels!(framework_from_cssr)
    @test isapprox(framework_from_cssr, framework, checknames=false)

    # test replicate framework
    sbmof = Framework("SBMOF-1.cif")
    replicated_sbmof = replicate(sbmof, (1, 1, 1))
    @test isapprox(sbmof, replicated_sbmof)
    # test replication no bonds assertion
    sbmof_bonds = Framework("SBMOF-1.cif")
    bonding_rules = [BondingRule(:H, :*, 0.4, 1.2),
                     BondingRule(:Ca, :*, 0.4, 2.3),
                     BondingRule(:*, :*, 0.4, 1.9)]
    infer_bonds!(sbmof_bonds, bonding_rules)
    @test_throws AssertionError 

    repfactors = replication_factors(sbmof.box, 14.0)
    replicated_sbmof = replicate(sbmof, repfactors)
    @test replication_factors(replicated_sbmof.box, 14.0) == (1, 1, 1)
    @test isapprox(sbmof.atoms.xf[:, 1] ./ repfactors, replicated_sbmof.atoms.xf[:, 1])
    @test isapprox(replicated_sbmof.box.reciprocal_lattice, 2 * π * inv(replicated_sbmof.box.f_to_c))
    @test chemical_formula(sbmof) == chemical_formula(replicated_sbmof)
    @test isapprox(crystal_density(sbmof), crystal_density(replicated_sbmof), atol=1e-7)

    # test symmetry_rules
    # define default symmetry_rules
    symmetry_rules = [Array{AbstractString, 2}(undef, 3, 0) ["x", "y", "z"]]
    other_symmetry_rules = [Array{AbstractString, 2}(undef, 3, 0) ["y + z", "x + z", "x + y"]]
    symmetry_rules_two = [Array{AbstractString, 2}(undef, 3, 0) ["x" "y + z";
                                                                 "y" "x + z";
                                                                 "z" "x + y"]]
    symmetry_rules_two_cpy = deepcopy(symmetry_rules_two)
    @test ! is_symmetry_equal(symmetry_rules, symmetry_rules_two)
    @test ! is_symmetry_equal(symmetry_rules, other_symmetry_rules)
    @test is_symmetry_equal(symmetry_rules, symmetry_rules)
    @test is_symmetry_equal(symmetry_rules_two, symmetry_rules_two_cpy)

    # test framework addition
    f1 = Framework("framework 1", UnitCube(), Atoms([:a, :b],
                                                    [1.0 4.0;
                                                     2.0 5.0;
                                                     3.0 6.0]),
                                              Charges([0.1, 0.2],
                                                      [1.0 4.0;
                                                       2.0 5.0;
                                                       3.0 6.0]))
    f2 = Framework("framework 2", UnitCube(), Atoms([:c, :d],
                                                    [7.0 10.0;
                                                     8.0 11.0;
                                                     9.0 12.0]),
                                              Charges([0.3, 0.4],
                                                      [7.0 10.0;
                                                       8.0 11.0;
                                                       9.0 12.0]))
    f3 = f1 + f2
    @test_throws AssertionError f1 + sbmof # only allow frameworks with same box
    @test isapprox(f1.box, f3.box)
    @test isapprox(f2.box, f3.box)
    @test isapprox(f1.atoms + f2.atoms, f3.atoms)
    @test isapprox(f1.charges + f2.charges, f3.charges)

    # test infinite framework_addition
    f4 = +(f1, f2, f3; check_overlap=false)
    @test isapprox(f3.box, f4.box)
    @test isapprox(f4.atoms, f1.atoms + f2.atoms + f3.atoms)
    @test isapprox(f4.charges, f1.charges + f2.charges + f3.charges)

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
