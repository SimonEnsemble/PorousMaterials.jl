module Crystal_Test

using PorousMaterials
using Test

@testset "Crystal Tests" begin
    # cif reader
    xtal = Crystal("test_structure2.cif")
    strip_numbers_from_atom_labels!(xtal)
    @test xtal.name == "test_structure2.cif"
    @test isapprox(xtal.box, Box(10.0, 20.0, 30.0, 90*π/180, 45*π/180, 120*π/180))
    @test xtal.atoms.n == 2
    @test isapprox(xtal.atoms, Atoms([:Ca, :O], Frac([0.2 0.6; 0.5 0.3; 0.7 0.1])))
    @test isapprox(xtal.charges, Charges([1.0, -1.0], Frac([0.2 0.6; 0.5 0.3; 0.7 0.1])))
    @test xtal.symmetry.is_p1
    
    # assign charges function
    xtal = Crystal("test_structure3.cif")
    @test xtal.charges.n == 0
    strip_numbers_from_atom_labels!(xtal)
    xtal2 = assign_charges(xtal, Dict(:Ca => -2.0, :O => 2.0))
    @test isapprox(xtal2.charges, Charges([-2.0, 2.0], Frac([0.2 0.6; 0.5 0.3; 0.7 0.1])))
    @test ! neutral(assign_charges(xtal, Dict(:Ca => 2.0, :O => 2.0), 100.0)) # not charge neutral
    @test_throws ErrorException assign_charges(xtal, Dict(:Ca => 2.0, :O => 2.0)) # not charge neutral
    @test chemical_formula(xtal) == Dict(:Ca => 1, :O => 1)
    @test molecular_weight(xtal) ≈ 15.9994 + 40.078

    # overlap checker, charge neutrality checker, wrap coords checker
    @test_throws ErrorException Crystal("test_structure2B.cif", check_overlap=true, check_neutrality=false)
    @test_throws ErrorException Crystal("test_structure2B.cif", check_overlap=false, check_neutrality=true)
    xtal = Crystal("test_structure2B.cif", check_overlap=false, check_neutrality=false, wrap_coords=true)
    f = Frac([0.2 0.5 0.7;
              0.2 0.5 0.7;
              0.6 0.3 0.1;
              0.9 0.24 0.15]')

    @test isapprox(xtal.atoms, Atoms([:Ca, :Ca, :O, :C], f))
    @test isapprox(net_charge(xtal), 2.0)
    @test chemical_formula(xtal) == Dict(:Ca => 2, :O => 1, :C => 1)

    # xyz writer
    xtal = Crystal("SBMOF-1.cif")
    xtal_xyz_filename = replace(replace(xtal.name, ".cif" => ""), ".cssr" => "") * ".xyz"
    write_xyz(xtal)
    atoms_read = read_xyz(xtal_xyz_filename) # Atoms{Cart}
    atoms_read_f = Frac(atoms_read, xtal.box) # Atoms{Frac}
    @test isapprox(atoms_read_f, xtal.atoms, atol=0.001)
    write_xyz(xtal, center=true, comment="blah") #center coords
    atoms_read = read_xyz(xtal_xyz_filename) # Atoms{Cart}
    atoms_read_f = Frac(atoms_read, xtal.box) # Atoms{Frac}
    @test isapprox(atoms_read_f.coords.xf, xtal.atoms.coords.xf .- [0.5, 0.5, 0.5], atol=0.001)
    rm(xtal_xyz_filename) # clean up

    # test .cif writer; write, read in, assert equal
    crystal = Crystal("SBMOF-1.cif")
    rewrite_filename = "rewritten.cif"
    write_cif(crystal, joinpath(PorousMaterials.PATH_TO_CRYSTALS, rewrite_filename))
    crystal_reloaded = Crystal(rewrite_filename)
    strip_numbers_from_atom_labels!(crystal_reloaded) # TODO Arthur remove this necessity from write_cif
    @test isapprox(crystal, crystal_reloaded, atol=0.0001)
    write_cif(crystal, joinpath(PorousMaterials.PATH_TO_CRYSTALS, rewrite_filename), fractional_coords=false) # cartesian
    crystal_reloaded = Crystal(rewrite_filename)
    strip_numbers_from_atom_labels!(crystal_reloaded) # TODO Arthur remove this necessity from write_cif
    @test isapprox(crystal, crystal_reloaded, atol=0.0001)
 #     rm(rewrite_filename) # clean up.

    # test .cif reader for non-P1 symmetry
    #   no atoms should overlap
    #   should place atoms in the same positions as the P1 conversion using
    #       openBabel
    non_P1_crystal = Crystal("symmetry_test_structure.cif")
    strip_numbers_from_atom_labels!(non_P1_crystal)
    non_P1_cartesian = Crystal("symmetry_test_structure_cartn.cif")
    strip_numbers_from_atom_labels!(non_P1_cartesian)
    P1_crystal = Crystal("symmetry_test_structure_P1.cif")
    strip_numbers_from_atom_labels!(P1_crystal)

    @test isapprox(non_P1_crystal, P1_crystal, atol=0.0001)
    # test that fractional and cartesian produce same results
    @test isapprox(non_P1_crystal, non_P1_cartesian, atol=0.0001)
    # test that cartesian and P1 produce same results
    @test isapprox(non_P1_cartesian, P1_crystal, atol=0.0001)
    # test that incorrect file formats throw proper errors
    @test_throws ErrorException Crystal("non_P1_no_symmetry.cif")
    # test that a file with no atoms throws error
    @test_throws ErrorException Crystal("no_atoms.cif")
    # should all be in P1 now.
    @test non_P1_crystal.symmetry.is_p1
    @test non_P1_cartesian.symmetry.is_p1
    @test P1_crystal.symmetry.is_p1

    # test reading in non-P1 then applying symmetry later
    # read in the same files as above, then convert to P1, then compare
    non_P1_crystal_symmetry = Crystal("symmetry_test_structure.cif", convert_to_p1=false)
    strip_numbers_from_atom_labels!(non_P1_crystal_symmetry)
    non_P1_cartesian_symmetry = Crystal("symmetry_test_structure_cartn.cif", convert_to_p1=false)
    strip_numbers_from_atom_labels!(non_P1_cartesian_symmetry)

    # make sure these crystals are not in P1 symmetry when convert_to_p1 is
    #   set to false
    @test ! non_P1_crystal_symmetry.symmetry.is_p1
    @test ! non_P1_cartesian_symmetry.symmetry.is_p1

    # test write_cif in non_p1 symmetry
    write_cif(non_P1_crystal_symmetry, joinpath("data", "crystals", "rewritten_symmetry_test_structure.cif"))
    # keep this in cartesian to test both
    write_cif(non_P1_cartesian_symmetry, joinpath("data", "crystals", "rewritten_symmetry_test_structure_cartn.cif"), fractional_coords=false)
    rewritten_non_p1_fractional = Crystal("rewritten_symmetry_test_structure.cif"; convert_to_p1=false)
    strip_numbers_from_atom_labels!(rewritten_non_p1_fractional) # TODO remove
    rewritten_non_p1_cartesian = Crystal("rewritten_symmetry_test_structure_cartn.cif"; convert_to_p1=false)
    strip_numbers_from_atom_labels!(rewritten_non_p1_cartesian)

    @test isapprox(rewritten_non_p1_fractional, non_P1_crystal_symmetry, atol=0.0001)
    @test isapprox(rewritten_non_p1_cartesian, non_P1_cartesian_symmetry, atol=0.0001)

    non_P1_crystal_symmetry = apply_symmetry_operations(non_P1_crystal_symmetry)
    non_P1_cartesian_symmetry = apply_symmetry_operations(non_P1_cartesian_symmetry)
    @test non_P1_crystal_symmetry.symmetry.is_p1
    @test non_P1_cartesian_symmetry.symmetry.is_p1

    # test that same structure is created when reading and converting to P1 and
    #   when reading then converting to P1
    @test isapprox(non_P1_crystal, non_P1_crystal_symmetry)
    @test isapprox(non_P1_cartesian, non_P1_cartesian_symmetry)
 # 
 #     # test .cssr reader too; test_structure2.{cif,cssr} designed to be the same.
 #     crystal_from_cssr = Crystal("test_structure2.cssr")
 #     strip_numbers_from_atom_labels!(crystal_from_cssr)
 #     @test has_same_sets_of_atoms_and_charges(crystal_from_cssr, crystal, checknames=false)
 # 
 #     # test replicate crystal
 #     sbmof = Crystal("SBMOF-1.cif")
 #     strip_numbers_from_atom_labels!(sbmof)
 #     replicated_sbmof = replicate(sbmof, (1, 1, 1))
 #     @test isapprox(sbmof, replicated_sbmof)
 #     # test replication no bonds assertion
 #     sbmof_bonds = Crystal("SBMOF-1.cif")
 #     strip_numbers_from_atom_labels!(sbmof_bonds)
 #     bonding_rules = [BondingRule(:H, :*, 0.4, 1.2),
 #                      BondingRule(:Ca, :O, 0.4, 2.5),
 #                      BondingRule(:*, :*, 0.4, 1.9)]
 #     infer_bonds!(sbmof_bonds, true, bonding_rules)
 #     @test_throws AssertionError replicate(sbmof_bonds, (2, 2, 2))
 #     # write out and compare to inferred bonds
 #     write_cif(sbmof_bonds, joinpath(pwd(), "data", "crystals", "SBMOF-1_inferred_bonds.cif"))
 #     sbmof_inferred_bonds = Crystal("SBMOF-1_inferred_bonds.cif"; read_bonds_from_file=true)
 #     strip_numbers_from_atom_labels!(sbmof_inferred_bonds)
 #     @test compare_bonds_in_crystal(sbmof_bonds, sbmof_inferred_bonds)
 #     # other bond info tests
 #     # TODO find more robust test/confirm these are the correct numbers
 #     # replacing this test with the one below comparing pdb bond info to inferred
 #     #   bond info
 #     sbmof_bonds_copy = Crystal("SBMOF-1.cif")
 #     strip_numbers_from_atom_labels!(sbmof_bonds_copy)
 #     # reverse the order of the atoms and bond info should still be the same
 #     sbmof_bonds_copy.atoms.xf .= reverse(sbmof_bonds_copy.atoms.xf; dims=2)
 #     sbmof_bonds_copy.atoms.species .= reverse(sbmof_bonds_copy.atoms.species)
 #     infer_bonds!(sbmof_bonds_copy, true, bonding_rules)
 #     @test compare_bonds_in_crystal(sbmof_bonds, sbmof_bonds_copy)
 #     remove_bonds!(sbmof_bonds)
 #     @test ne(sbmof_bonds.bonds) == 0
 #     @test !compare_bonds_in_crystal(sbmof_bonds, sbmof_bonds_copy)
 # 
 #     # test reading in bonds as part of `Crystal()`
 #     sbmof_read_bonds = Crystal("test_bond_viz.cif"; read_bonds_from_file=true, check_atom_and_charge_overlap=false)
 #     strip_numbers_from_atom_labels!(sbmof_read_bonds)
 #     @test ne(sbmof_read_bonds.bonds) == 5
 #     write_cif(sbmof_read_bonds, joinpath(pwd(), "data", "crystals", "rewritten_sbmof_read_bonds.cif"))
 #     reloaded_sbmof_read_bonds = Crystal("rewritten_sbmof_read_bonds.cif"; read_bonds_from_file=true, check_atom_and_charge_overlap=false)
 #     strip_numbers_from_atom_labels!(reloaded_sbmof_read_bonds)
 #     @test compare_bonds_in_crystal(sbmof_read_bonds, reloaded_sbmof_read_bonds)
 # 
 #     # Test that reading in bonding information is the same as inferring the
 #     #   bonding info
 #     # Bonding information is from a pdb file saved by avogadro
 #     # using non-p1 because it meant copying over fewer bonds
 #     read_bonds = Crystal("KAXQIL_clean_cartn.cif"; convert_to_p1=false, read_bonds_from_file=true)
 #     strip_numbers_from_atom_labels!(read_bonds)
 #     inferred_bonds = Crystal("KAXQIL_clean.cif"; convert_to_p1=false)
 #     strip_numbers_from_atom_labels!(inferred_bonds)
 #     # Using same bonding rules as above
 #     infer_bonds!(inferred_bonds, true, bonding_rules)
 #     @test compare_bonds_in_crystal(read_bonds, inferred_bonds; atol=1e-6)
 # 
 #     repfactors = replication_factors(sbmof.box, 14.0)
 #     replicated_sbmof = replicate(sbmof, repfactors)
 #     @test replication_factors(replicated_sbmof.box, 14.0) == (1, 1, 1)
 #     @test isapprox(sbmof.atoms.xf[:, 1] ./ repfactors, replicated_sbmof.atoms.xf[:, 1])
 #     @test isapprox(replicated_sbmof.box.reciprocal_lattice, 2 * π * inv(replicated_sbmof.box.f_to_c))
 #     @test chemical_formula(sbmof) == chemical_formula(replicated_sbmof)
 #     @test isapprox(crystal_density(sbmof), crystal_density(replicated_sbmof), atol=1e-7)
 # 
 #     # test symmetry_rules
 #     # define default symmetry_rules
 #     symmetry_rules = [Array{AbstractString, 2}(undef, 3, 0) ["x", "y", "z"]]
 #     other_symmetry_rules = [Array{AbstractString, 2}(undef, 3, 0) ["y + z", "x + z", "x + y"]]
 #     symmetry_rules_two = [Array{AbstractString, 2}(undef, 3, 0) ["x" "y + z";
 #                                                                  "y" "x + z";
 #                                                                  "z" "x + y"]]
 #     symmetry_rules_two_cpy = deepcopy(symmetry_rules_two)
 #     @test ! is_symmetry_equal(symmetry_rules, symmetry_rules_two)
 #     @test ! is_symmetry_equal(symmetry_rules, other_symmetry_rules)
 #     @test is_symmetry_equal(symmetry_rules, symmetry_rules)
 #     @test is_symmetry_equal(symmetry_rules_two, symmetry_rules_two_cpy)
 # 
 #     # test crystal addition
 #     f1 = Crystal("crystal 1", UnitCube(), Atoms([:a, :b],
 #                                                     [1.0 4.0;
 #                                                      2.0 5.0;
 #                                                      3.0 6.0]),
 #                                               Charges([0.1, 0.2],
 #                                                       [1.0 4.0;
 #                                                        2.0 5.0;
 #                                                        3.0 6.0]))
 #     f2 = Crystal("crystal 2", UnitCube(), Atoms([:c, :d],
 #                                                     [7.0 10.0;
 #                                                      8.0 11.0;
 #                                                      9.0 12.0]),
 #                                               Charges([0.3, 0.4],
 #                                                       [7.0 10.0;
 #                                                        8.0 11.0;
 #                                                        9.0 12.0]))
 #     f3 = f1 + f2
 #     addition_bonding_rules = [BondingRule(:a, :b, 4.5, 5.3),
 #                               BondingRule(:c, :d, 4.5, 5.3)]
 #     @test is_bonded(f1, 1, 2, [BondingRule(:a, :b, 1.0, 5.5)]; include_bonds_across_periodic_boundaries=false)
 #     @test ! is_bonded(f2, 1, 2, [BondingRule(:c, :d, 1.0, 4.5)]; include_bonds_across_periodic_boundaries=false)
 #     infer_bonds!(f1, false, addition_bonding_rules)
 #     infer_bonds!(f2, false, addition_bonding_rules)
 #     @test ! compare_bonds_in_crystal(f1 + f2, f3)
 #     infer_bonds!(f3, false, addition_bonding_rules)
 #     @test compare_bonds_in_crystal(f1 + f2, f3)
 #     @test_throws AssertionError f1 + sbmof # only allow crystals with same box
 #     @test isapprox(f1.box, f3.box)
 #     @test isapprox(f2.box, f3.box)
 #     @test isapprox(f1.atoms + f2.atoms, f3.atoms)
 #     @test isapprox(f1.charges + f2.charges, f3.charges)
 # 
 #     # test infinite crystal_addition
 #     f4 = +(f1, f2, f3; check_overlap=false)
 #     @test isapprox(f3.box, f4.box)
 #     @test isapprox(f4.atoms, f1.atoms + f2.atoms + f3.atoms)
 #     @test isapprox(f4.charges, f1.charges + f2.charges + f3.charges)
 # 
 #     # more xtal tests
 #     sbmof1 = Crystal("SBMOF-1.cif")
 #     @test !charged(sbmof1)
 #     @test isapprox(sbmof1.box.reciprocal_lattice, 2 * π * inv(sbmof1.box.f_to_c))
 #     @test sbmof1.box.Ω ≈ det(sbmof1.box.f_to_c) # sneak in crystal test
 #     @test isapprox(crystal_density(sbmof1), 1570.4, atol=0.5) # kg/m3
 # 
 #     # replicating the unit cell to construct simulation box
 #     sbmof1 = Crystal("SBMOF-1.cif")
 #     rbox = replicate(sbmof1.box, (2, 3, 4))
 #     @test rbox.Ω ≈ sbmof1.box.Ω * 2 * 3 * 4
 #     @test all(rbox.c_to_f * sbmof1.box.f_to_c * [1.0, 1.0, 1.0] .≈ [1/2, 1/3, 1/4])
 # 
 #     # write bond information to manually inspect if bonds are in order
 #     hkust1 = Crystal("FIQCEN_clean.cif")
 #     strip_numbers_from_atom_labels!(hkust1)
 #     br = default_bondingrules()
 #     pushfirst!(br, BondingRule(:Cu, :O, 0.3, 2.4))
 #     infer_bonds!(hkust1, true, br)
 #     reference_bonds = loadgraph("hkust1_reference_bonds.lgz")
 #     @test reference_bonds == hkust1.bonds
 # 
 #     # test distance function (via Avogadro)
 #     frame = Crystal("simple_test.cif")
 #     @test distance(frame, 1, 1, true) == 0.0
 #     @test isapprox(distance(frame, 2, 5, true), 4.059, atol=0.001)
 #     @test isapprox(distance(frame, 2, 5, false), 4.059, atol=0.001)
 #     @test isapprox(distance(frame, 1, 5, false), 17.279, atol=0.001)
 #     @test isapprox(distance(frame, 1, 5, true), 1.531, atol=0.001)
 # 
    # test reading crystals and include 0.0 charges
    frame1 = Crystal("ATIBOU01_clean.cif")
    @test ! any(frame1.charges.q .== 0.0)
    @test frame1.charges.n == frame1.atoms.n - 4 # 4 charges are zero
    frame2 = Crystal("ATIBOU01_clean.cif"; include_zero_charges=true)
    @test frame2.charges.n == frame2.atoms.n
    @test isapprox(frame2.charges.coords, frame2.atoms.coords)
 # 
 # #=    ## .cssr reader test
 #     # test replicate crystal
 #     sbmof = Crystal("SBMOF-1.cif")
 #     replicated_sbmof = replicate(sbmof, (1, 1, 1))
 #     @test isapprox(sbmof, replicated_sbmof)
 # 
 #     repfactors = replication_factors(sbmof.box, 14.0)
 #     replicated_sbmof = replicate(sbmof, repfactors)
 #     @test replication_factors(replicated_sbmof.box, 14.0) == (1, 1, 1)
 #     @test isapprox(sbmof.atoms.xf[:, 1] ./ repfactors, replicated_sbmof.atoms.xf[:, 1])
 #     @test isapprox(replicated_sbmof.box.reciprocal_lattice, 2 * π * inv(replicated_sbmof.box.f_to_c))
 #     @test chemical_formula(sbmof) == chemical_formula(replicated_sbmof)
 #     @test isapprox(crystal_density(sbmof), crystal_density(replicated_sbmof), atol=1e-7)
 # 
 #     # more xtal tests
 #     sbmof1 = Crystal("SBMOF-1.cif")
 #     @test !charged(sbmof1)
 #     @test isapprox(sbmof1.box.reciprocal_lattice, 2 * π * inv(sbmof1.box.f_to_c))
 #     @test sbmof1.box.Ω ≈ det(sbmof1.box.f_to_c) # sneak in crystal test
 #     @test isapprox(crystal_density(sbmof1), 1570.4, atol=0.5) # kg/m3
 # 
 #     # replicating the unit cell to construct simulation box
 #     sbmof1 = Crystal("SBMOF-1.cif")
 #     rbox = replicate(sbmof1.box, (2, 3, 4))
 #     @test rbox.Ω ≈ sbmof1.box.Ω * 2 * 3 * 4
 #     @test all(rbox.c_to_f * sbmof1.box.f_to_c * [1.0, 1.0, 1.0] .≈ [1/2, 1/3, 1/4])=#
end
end
