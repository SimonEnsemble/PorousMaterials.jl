#!/usr/bin/env julia

# Details from http://www.stochasticlifestyle.com/finalizing-julia-package-documentation-testing-coverage-publishing/
# Start Test Script
using PorousMaterials
using Base.Test

# Run Tests

@printf("------------------------------\nTesting Crystal.jl\n\n")
framework = read_crystal_structure_file("test_structure2.cif")
strip_numbers_from_atom_labels!(framework)
@testset "Crystal Tests" begin
    @test framework.name == "test_structure2.cif"
	@test framework.box.f_to_c * framework.box.c_to_f ≈ eye(3)
	@test framework.box.a ≈ 10.0
	@test framework.box.b ≈ 20.0
	@test framework.box.c ≈ 30.0
	@test framework.box.α ≈ 90.0 * (pi / 180)
	@test framework.box.β ≈ 45.0 * (pi / 180)
	@test framework.box.γ ≈ 120.0 * (pi / 180)
	@test framework.box.Ω ≈ det(framework.box.f_to_c)
    @test framework.n_atoms == 2
	@test framework.atoms == ["Ca", "O"]
	@test framework.xf == [0.2 0.6; 0.5 0.3; 0.7 0.1]
    @test framework.charges == [1.0, -1.0]
    @test chemical_formula(framework) == Dict("Ca" => 1, "O" => 1)
    @test molecular_weight(framework) ≈ 15.9994 + 40.078

    # test .cssr reader too; test_structure2.{cif,cssr} designed to be the same.
    framework_from_cssr = read_crystal_structure_file("test_structure2.cif")
    strip_numbers_from_atom_labels!(framework_from_cssr)
    @test all(framework_from_cssr.xf .== framework.xf)
    @test all(framework_from_cssr.charges .== framework.charges)
    @test all(framework_from_cssr.atoms == framework.atoms)
    @test framework_from_cssr.n_atoms == framework.n_atoms
    @test all(framework_from_cssr.box.f_to_c .== framework.box.f_to_c)
    @test all(framework_from_cssr.box.c_to_f .== framework.box.c_to_f)
end;

@printf("------------------------------\nTesting Forcefield.jl\n\n")
const ljforcefield = read_forcefield_file("test_forcefield.csv", cutoffradius=12.5, mixing_rules="Lorentz-Berthelot") # Dreiding
frame = read_crystal_structure_file("test_structure.cif") # .cif
strip_numbers_from_atom_labels!(frame)
rep_factors = replication_factors(frame.box, ljforcefield)
@testset "Forcefield Tests" begin
	@test ljforcefield.pure_σ["He"] == 1.0
	@test ljforcefield.pure_ϵ["Zn"] == 12.0
	@test ljforcefield.σ²["Zn"]["He"] == ((1.0 + 3.0) / 2) ^ 2
	@test ljforcefield.ϵ["He"]["Zn"] == sqrt(12.0 * 3.0)
	@test ljforcefield.ϵ["He"]["Zn"] == ljforcefield.ϵ["Zn"]["He"] # symmetry
	@test ljforcefield.σ²["He"]["Zn"] == ljforcefield.σ²["Zn"]["He"] # symmetry
	@test ljforcefield.cutoffradius_squared == 12.5 ^ 2
	@test rep_factors == (25, 25, 25)
end;

@printf("------------------------------\nTesting Energetics.jl\n\n")
@testset "Energetics Tests" begin
    # test Periodic boundary conditions
    molecule1 = Molecule(1, ["He"], [0.5, 0.5, 0.5][:,:], [0.0])
    molecule2 = Molecule(1, ["He"], [0.5 + rep_factors[1], 0.5 + rep_factors[2], 0.5 + rep_factors[3]][:,:], [0.0])
	@test vdw_energy(frame, molecule1, ljforcefield, rep_factors) ≈ vdw_energy(frame, molecule2, ljforcefield, rep_factors)
	@test vdw_energy(frame, molecule1, ljforcefield, (1,1,1)) ≈ 4 * ljforcefield.ϵ["He"]["Zn"] * ((ljforcefield.σ²["Zn"]["He"] / 0.75) ^ 6 - (ljforcefield.σ²["Zn"]["He"] / 0.75) ^ 3)
    # Xe in SBMOF-1 tests, comparing to RASPA
    sbmof1 = read_crystal_structure_file("SBMOF-1.cif")
    rep_factors_sbmof1 = replication_factors(sbmof1.box, ljforcefield)
    xenon = Molecule(1, ["Xe"], zeros(3, 1), [0.0])
    energy = vdw_energy(sbmof1, xenon, ljforcefield, rep_factors_sbmof1)
	@test isapprox(energy, -5041.58, atol = 0.005)
    xenon.x[1] = 0.494265; xenon.x[2] = 2.22668; xenon.x[3] = 0.450354;
    energy = vdw_energy(sbmof1, xenon, ljforcefield, rep_factors_sbmof1)
	@test isapprox(energy, 12945.838, atol = 0.005)

    # test PBC, rep factors are (3, 5, 2)
    xf = [0.05, 0.4, 0.02][:, :]
    xenon.x = sbmof1.box.f_to_c * xf
    energy1 = vdw_energy(sbmof1, xenon, ljforcefield, rep_factors_sbmof1)
    xf = [1.05, 0.4, 0.02][:, :]
    xenon.x = sbmof1.box.f_to_c * xf
    energy2 = vdw_energy(sbmof1, xenon, ljforcefield, rep_factors_sbmof1)
	@test isapprox(energy1, energy2, atol = 0.00001)
    xf = [1.05, 4.4, 1.02][:, :]
    xenon.x = sbmof1.box.f_to_c * xf
    energy3 = vdw_energy(sbmof1, xenon, ljforcefield, rep_factors_sbmof1)
	@test isapprox(energy1, energy3, atol = 0.00001)
    # outside box
    xf = [4.05, 5.4, 2.02][:, :]
    xenon.x = sbmof1.box.f_to_c * xf
    energy4 = vdw_energy(sbmof1, xenon, ljforcefield, rep_factors_sbmof1)
	@test isapprox(energy1, energy4, atol = 0.00001)
    xf = [-0.95, 5.4, 2.02][:, :]
    xenon.x = sbmof1.box.f_to_c * xf
    energy5 = vdw_energy(sbmof1, xenon, ljforcefield, rep_factors_sbmof1)
	@test isapprox(energy1, energy5, atol = 0.00001)
    xf = [-0.95, -0.6, -0.98][:, :]
    xenon.x = sbmof1.box.f_to_c * xf
    energy6 = vdw_energy(sbmof1, xenon, ljforcefield, rep_factors_sbmof1)
	@test isapprox(energy1, energy6, atol = 0.00001)
end;

@printf("------------------------------\n")

@printf("------------------------------\nTesting GCMC.jl\n\n")
@testset "Monte Carlo Functions Tests" begin
    #
    #INSERTION TESTS
    #
    insertion_inside_box = true
    insertion_at_random_coords = true
    insertion_adds_molecule = true

    molecules = Array{Molecule}(0)
    repfactors = replication_factors(frame.box, ljforcefield)
    sim_box = replicate_box(frame.box, repfactors)

    for i = 1:100
        insert_molecule!(molecules, sim_box, "C")
        if completely_outside_box(molecules[i], sim_box)
            insertion_inside_box = false
        end
        if ! (length(molecules) == i)
            insertion_adds_molecule = false
        end
        if i > 1
            if sum(molecules[i - 1].x .≈ molecules[i].x) != 0
                # by chance this could fail but highly unlikely!
                insertion_at_random_coords = false
            end
        end
    end
    @test insertion_inside_box
    @test insertion_at_random_coords
    @test insertion_adds_molecule

    #
    #DELETION TESTS
    #
    deletion_removes_a_molecule = true
    for i = 1:100
        delete_molecule!(rand(1:length(molecules)), molecules)
        if length(molecules) != 100 - i
            deletion_removes_a_molecule = false
        end
    end
    @test deletion_removes_a_molecule

    #
    #TRANSLATION TESTS
    #
    translation_old_coords_stored_properly = true
    translation_coords_changed = true
    translation_inside_box = true
    molecules = [Molecule(1, ["C"], sim_box.f_to_c * [0.99, 0.99, 0.01][:, :], [0.0]),
                 Molecule(1, ["F"], sim_box.f_to_c * [0.01, 0.01, 0.99][:, :], [0.0])]
    x_old = translate_molecule!(1, molecules, sim_box)
    if ! all(x_old .== sim_box.f_to_c * [0.99, 0.99, 0.01][:, :])
        translation_old_coords_stored_properly = false
    end
    if ! all(molecules[1].x .!= x_old)
        translation_coords_changed = false
    end

    for i = 1:100000
        which_molecule = rand(1:2) # choose molecule to move
        xf_old_should_be = deepcopy(molecules[which_molecule].x)
        xf_old = translate_molecule!(which_molecule, molecules, sim_box)
        if ! all(molecules[which_molecule].x .!= x_old)
            translation_coords_changed = false
        end
        if completely_outside_box(molecules[which_molecule], sim_box)
            translation_inside_box = false
        end
    end
    @test translation_old_coords_stored_properly
    @test translation_coords_changed
    @test translation_inside_box
end
@printf("------------------------------\n")
@testset "Guest-guest Energetics Tests" begin
    sim_box = construct_box(25.0, 25.0, 25.0, π/2, π/2, π/2)
    # distance of 6.0 away
    molecules = [Molecule(1, ["C"], [5.0, 12.0, 12.0][:, :], [0.0]),
                 Molecule(1, ["O"], [11.0, 12.0, 12.0][:, :], [0.0])]
    r² = (11.0 - 5.0) ^ 2 # duh
    energy = lennard_jones(r², ljforcefield.σ²["C"]["O"], ljforcefield.ϵ["O"]["C"])
    @test energy ≈ guest_guest_vdw_energy(1, molecules, ljforcefield, sim_box)
    # symmetry
    @test energy ≈ guest_guest_vdw_energy(2, molecules, ljforcefield, sim_box)
    
    # via PBC, a distance (24.0 - 5.0) > (1+5)
    molecules[2] = Molecule(1, ["O"], [24.0, 12.0, 12.0][:, :], [0.0])
    r² = (1.0 + 5.0) ^ 2 # PBC
    energy = lennard_jones(r², ljforcefield.σ²["C"]["O"], ljforcefield.ϵ["C"]["O"])
    @test energy ≈ guest_guest_vdw_energy(2, molecules, ljforcefield, sim_box)
    @test energy ≈ guest_guest_vdw_energy(1, molecules, ljforcefield, sim_box) # symmetry again.
    
    # put a molecule on top of first one.
    push!(molecules, Molecule(1, ["C"], [5.0, 12.0, 12.0][:, :], [0.0]))
    @test guest_guest_vdw_energy(2, molecules, ljforcefield, sim_box) ≈ 2 * energy
    
    @test guest_guest_vdw_energy(1, molecules, ljforcefield, sim_box) == Inf
    @test guest_guest_vdw_energy(3, molecules, ljforcefield, sim_box) == Inf
    
    molecules_a = [Molecule(1, ["C"], [11.0, 1.0, 12.0][:, :], [0.0]),
                   Molecule(1, ["O"], [11.0, 4.0, 12.0][:, :], [0.0])]
    molecules_b = [Molecule(1, ["C"], [11.0, 1.0, 12.0][:, :], [0.0]),
                   Molecule(1, ["O"], [11.0, 23.0, 12.0][:, :], [0.0])]
    @test guest_guest_vdw_energy(1, molecules_a, ljforcefield, sim_box) ≈ guest_guest_vdw_energy(1, molecules_b, ljforcefield, sim_box)

end
