#!/usr/bin/env julia

#TODO split these into different files for different tests so e.g. framework in one suite is not used in another suite so it's clear what's goign on.

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
    @test isapprox(framework.box, construct_box(framework.box.f_to_c))
	@test framework.box.a ≈ 10.0
	@test framework.box.b ≈ 20.0
	@test framework.box.c ≈ 30.0
	@test framework.box.α ≈ 90.0 * (pi / 180)
	@test framework.box.β ≈ 45.0 * (pi / 180)
	@test framework.box.γ ≈ 120.0 * (pi / 180)
	@test framework.box.Ω ≈ det(framework.box.f_to_c)
    @test framework.n_atoms == 2
	@test framework.atoms == [:Ca, :O]
	@test framework.xf == [0.2 0.6; 0.5 0.3; 0.7 0.1]
    @test framework.charges == [1.0, -1.0]
    @test charged(framework)
    @test chemical_formula(framework) == Dict(:Ca => 1, :O => 1)
    @test molecular_weight(framework) ≈ 15.9994 + 40.078
    @test isapprox(transpose(framework.box.reciprocal_lattice), 2 * π * inv(framework.box.f_to_c))

    # test .cssr reader too; test_structure2.{cif,cssr} designed to be the same.
    framework_from_cssr = read_crystal_structure_file("test_structure2.cif")
    strip_numbers_from_atom_labels!(framework_from_cssr)
    @test isapprox(framework_from_cssr, framework, checknames=false)
end;

@printf("------------------------------\nTesting Forcefield.jl\n\n")
const ljforcefield = read_forcefield_file("Dreiding.csv", cutoffradius=12.5, mixing_rules="Lorentz-Berthelot") # Dreiding
frame = read_crystal_structure_file("test_structure.cif") # .cif
strip_numbers_from_atom_labels!(frame)
rep_factors = replication_factors(frame.box, ljforcefield)
@testset "Forcefield Tests" begin
	@test ljforcefield.pure_σ[:He] == 1.0
	@test ljforcefield.pure_ϵ[:Zn] == 12.0
	@test ljforcefield.σ²[:Zn][:He] == ((1.0 + 3.0) / 2) ^ 2
	@test ljforcefield.ϵ[:He][:Zn] == sqrt(12.0 * 3.0)
	@test ljforcefield.ϵ[:He][:Zn] == ljforcefield.ϵ[:Zn][:He] # symmetry
	@test ljforcefield.σ²[:He][:Zn] == ljforcefield.σ²[:Zn][:He] # symmetry
	@test ljforcefield.cutoffradius_squared == 12.5 ^ 2
	@test rep_factors == (25, 25, 25)
    # force field coverage function
    framework10 = read_crystal_structure_file("SBMOF-1.cif")
    @test isapprox(framework10.box, construct_box(framework10.box.f_to_c))
    @test check_forcefield_coverage(framework10, ljforcefield)
    push!(framework10.atoms, :bogus_atom)
    @test !check_forcefield_coverage(framework10, ljforcefield)
end;

@printf("------------------------------\nTesting Molecules.jl\n\n")
@testset "Molecules Tests" begin
    # test reader
    molecule = read_molecule_file("CO2")
    @test check_forcefield_coverage(molecule, ljforcefield)
    @test charged(molecule)
    atomic_masses = read_atomic_masses()
    @test molecule.species == :CO2
    @test length(molecule.ljspheres) == 3
    @test molecule.ljspheres[1].atom == :C_CO2
    @test molecule.ljspheres[2].atom == :O_CO2
    @test molecule.ljspheres[3].atom == :O_CO2
    @test all(molecule.ljspheres[1].x .≈ [0.0, 0.0, 0.0])
    @test all(molecule.ljspheres[2].x .≈ [-1.16, 0.0, 0.0])
    @test all(molecule.ljspheres[3].x .≈ [1.16, 0.0, 0.0])
    @test all(molecule.center_of_mass .≈ [0.0, 0.0, 0.0])
    @test length(molecule.charges) == 3
    @test molecule.charges[1].q ≈ 0.7
    @test molecule.charges[2].q ≈ -0.35
    @test molecule.charges[3].q ≈ -0.35
    for i = 1:3
        @test all(molecule.charges[i].x ≈ molecule.ljspheres[i].x)
    end

    # test translate
    m1 = read_molecule_file("CO2")
    m2 = read_molecule_file("CO2")
    @test isapprox(m1, m2) # overloaded this function for molecules
    translate_by!(m2, [0.0, 0.0, 0.0])
    @test isapprox(m1, m2)
    translate_by!(m2, [0.0, 10.0, 0.0])
    @test ! isapprox(m1, m2)
    translate_to!(m2, m1.center_of_mass)
    @test isapprox(m1, m2)
    translate_to!(m2, [50.0, 100.0, 150.0])
    @test isapprox(m2.center_of_mass, [50.0, 100.0, 150.0])
    for i = 1:200
        translate_by!(m2, [randn(), randn(), randn()])
    end
    @test norm(m2.charges[1].x - m2.charges[2].x) ≈ norm(m1.charges[1].x - m1.charges[2].x)
    @test norm(m2.ljspheres[1].x - m2.ljspheres[2].x) ≈ norm(m1.ljspheres[1].x - m1.ljspheres[2].x)

    # test unit vector on sphere generator
    ms = [read_molecule_file("He") for i = 1:10000]
    for m in ms
        translate_to!(m, rand_point_on_unit_sphere())
    end
    @test all(isapprox.([norm(m.ljspheres[1].x) for m in ms], 1.0))
    write_to_xyz(ms, "random_vectors_on_sphere")
    println("See random_vectors_on_sphere")
    
    # rotation matrix shld be orthogonal
    r_orthogonal = true
    r_det_1 = true
    for i = 1:300
        r = rotation_matrix()
        if ! isapprox(r * transpose(r), eye(3))
            r_orthogonal = false
        end
        if ! isapprox(det(r), 1.0)
            r_det_1 = false
        end
    end
    @test r_orthogonal
    @test r_det_1
    
    # test rotate function
    translate_to!(m2, [50.0, 100.0, 150.0])
    for i = 1:2000
        rotate!(m2)
    end
    @test isapprox(m2.center_of_mass, [50.0, 100.0, 150.0])
    @test norm(m2.charges[1].x - m2.charges[2].x) ≈ norm(m1.charges[1].x - m1.charges[2].x)
    @test norm(m2.ljspheres[1].x - m2.ljspheres[2].x) ≈ norm(m1.ljspheres[1].x - m1.ljspheres[2].x)
    m2_old = deepcopy(m2)
    rotate!(m2)
    @test ! isapprox(m2_old, m2)
    # visually inspect
    ms = [read_molecule_file("CO2") for i = 1:1000]
    for m in ms
       rotate!(m)
    end
    write_to_xyz(ms, "co2s")
    println("see co2s.xyz for dist'n of rotations")


end

@printf("------------------------------\nTesting Energetics.jl\n\n")
@testset "Energetics Tests" begin
    # test Periodic boundary conditions
    molecule1 = read_molecule_file("He")
    molecule1.ljspheres[1].x[:] = [0.5, 0.5, 0.5]
    molecule2 = read_molecule_file("He")
    molecule2.ljspheres[1].x[:] = [0.5 + rep_factors[1], 0.5 + rep_factors[2], 0.5 + rep_factors[3]]

	@test vdw_energy(frame, molecule1, ljforcefield, rep_factors) ≈ vdw_energy(frame, molecule2, ljforcefield, rep_factors)
	@test vdw_energy(frame, molecule1, ljforcefield, (1,1,1)) ≈ 4 * ljforcefield.ϵ[:He][:Zn] * ((ljforcefield.σ²[:Zn][:He] / 0.75) ^ 6 - (ljforcefield.σ²[:Zn][:He] / 0.75) ^ 3)
    # the position of a molecule should not change inside vdw_energy.
    @assert(all(molecule1.ljspheres[1].x .≈ [0.5, 0.5, 0.5]))
    @assert(all(molecule2.ljspheres[1].x .≈ [0.5 + rep_factors[1], 0.5 + rep_factors[2], 0.5 + rep_factors[3]]))

    # Xe in SBMOF-1 tests, comparing to RASPA
    sbmof1 = read_crystal_structure_file("SBMOF-1.cif")
    @test !charged(sbmof1)
    @test isapprox(transpose(sbmof1.box.reciprocal_lattice), 2 * π * inv(sbmof1.box.f_to_c))
    @test sbmof1.box.Ω ≈ det(sbmof1.box.f_to_c) # sneak in crystal test
    @test isapprox(crystal_density(sbmof1), 1570.4, atol=0.5) # kg/m3
    rep_factors_sbmof1 = replication_factors(sbmof1.box, ljforcefield)
    xenon = read_molecule_file("Xe")
    @test ! charged(xenon)
    xenon.ljspheres[1].x[:] = zeros(3)
    energy = vdw_energy(sbmof1, xenon, ljforcefield, rep_factors_sbmof1)
	@test isapprox(energy, -5041.58, atol = 0.005)
    xenon.ljspheres[1].x[1] = 0.494265; xenon.ljspheres[1].x[2] = 2.22668; xenon.ljspheres[1].x[3] = 0.450354;
    energy = vdw_energy(sbmof1, xenon, ljforcefield, rep_factors_sbmof1)
	@test isapprox(energy, 12945.838, atol = 0.005)

    # test PBC, rep factors are (3, 5, 2)
    xf = [0.05, 0.4, 0.02][:, :]
    xenon.ljspheres[1].x[:] = sbmof1.box.f_to_c * xf
    energy1 = vdw_energy(sbmof1, xenon, ljforcefield, rep_factors_sbmof1)
    xf = [1.05, 0.4, 0.02][:, :]
    xenon.ljspheres[1].x[:] = sbmof1.box.f_to_c * xf
    energy2 = vdw_energy(sbmof1, xenon, ljforcefield, rep_factors_sbmof1)
	@test isapprox(energy1, energy2, atol = 0.00001)
    xf = [2.05, 4.4, 1.02][:, :]
    xenon.ljspheres[1].x[:] = sbmof1.box.f_to_c * xf
    energy3 = vdw_energy(sbmof1, xenon, ljforcefield, rep_factors_sbmof1)
	@test isapprox(energy1, energy3, atol = 0.00001)
    # outside box
    xf = [4.05, 5.4, 2.02][:, :]
    xenon.ljspheres[1].x[:] = sbmof1.box.f_to_c * xf
    energy4 = vdw_energy(sbmof1, xenon, ljforcefield, rep_factors_sbmof1)
	@test isapprox(energy1, energy4, atol = 0.00001)
    xf = [-0.95, 5.4, 2.02][:, :]
    xenon.ljspheres[1].x[:] = sbmof1.box.f_to_c * xf
    energy5 = vdw_energy(sbmof1, xenon, ljforcefield, rep_factors_sbmof1)
	@test isapprox(energy1, energy5, atol = 0.00001)
    xf = [-0.95, -0.6, -0.98][:, :]
    xenon.ljspheres[1].x[:] = sbmof1.box.f_to_c * xf
    energy6 = vdw_energy(sbmof1, xenon, ljforcefield, rep_factors_sbmof1)
	@test isapprox(energy1, energy6, atol = 0.00001)

    # data types for potential energies
    v = PotentialEnergy(3.0, 4.0, 5.0, 6.0)
    @test isapprox(sum(v), 18.0)
    @test isapprox(v, PotentialEnergy(3.0, 4.0, 5.0, 6.0))
    @test ! isapprox(v, PotentialEnergy(3.1, 4.0, 5.0, 6.0))
    @test isapprox(v + PotentialEnergy(2.0, 3.0, 4.0, 5.0), PotentialEnergy(5.0, 7.0, 9.0, 11.0))
    @test isapprox(v - PotentialEnergy(2.0, 3.0, 4.0, 5.0), PotentialEnergy(1.0, 1.0, 1.0, 1.0))
end;
@printf("------------------------------\n")


@printf("------------------------------\nTesting Electrostatics\n\n")
framework = read_crystal_structure_file("NU-1000_Greg.cif")

kreps = (11, 11, 9)
α = 0.265058
sr_cutoff_r = 12.5
rep_factors = replication_factors(framework, sr_cutoff_r)
sim_box = replicate_box(framework.box, rep_factors)
eparams, kvec_wts, eikar, eikbr, eikcr = setup_Ewald_sum(kreps, α, sr_cutoff_r, sim_box)

q_test = 0.8096

@testset "Ewald summation Tests" begin
    x = [9.535619863743, 20.685576379935, 0.127344239990]
    ϕ = electrostatic_potential(framework, x, rep_factors, eparams, kvec_wts, eikar, eikbr, eikcr)
    @test isapprox(ϕ * q_test, 111373.38, atol=2.5)

    x = [4.269654927228, 23.137319129548, 28.352847101096]
    ϕ = electrostatic_potential(framework, x, rep_factors, eparams, kvec_wts, eikar, eikbr, eikcr)
    @test isapprox(ϕ * q_test, -531.0, atol=0.5)

    x = [-0.047382031804, 7.209555961450, 5.158180463556]
    ϕ = electrostatic_potential(framework, x, rep_factors, eparams, kvec_wts, eikar, eikbr, eikcr)
    @test isapprox(ϕ * q_test, -2676.8230141, atol=0.5)
end

@printf("------------------------------\n")

@printf("------------------------------\nTesting GCMC.jl\n\n")
@testset "Monte Carlo Functions Tests" begin
    # replicating the unit cell to construct simulation box
    sbmof1 = read_crystal_structure_file("SBMOF-1.cif")
    sim_box = replicate_box(sbmof1.box, (1, 1, 1))
    @test isapprox(sim_box, sbmof1.box)
    sim_box = replicate_box(sbmof1.box, (2, 3, 4))
    @test sim_box.Ω ≈ sbmof1.box.Ω * 2 * 3 * 4
    @test all(sim_box.c_to_f * sbmof1.box.f_to_c * [1.0, 1.0, 1.0] .≈ [1/2, 1/3, 1/4])

    #
    #INSERTION TESTS
    #
    insertion_inside_box = true
    insertion_at_random_coords = true
    insertion_adds_molecule = true

    molecules = Array{Molecule}(0)
    repfactors = replication_factors(frame.box, ljforcefield)
    sim_box = replicate_box(frame.box, repfactors)
    
    m = read_molecule_file("He")
    for i = 1:100
        insert_molecule!(molecules, sim_box, m)
        if outside_box(molecules[i], sim_box)
            insertion_inside_box = false
        end
        if ! (length(molecules) == i)
            insertion_adds_molecule = false
        end
        if i > 1
            if isapprox(molecules[i - 1], molecules[i])
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
    # first, test function to bring molecule inside a box.
    box = construct_box(25.0, 25.0, 25.0, π/2, π/2, π/2)
    molecule = read_molecule_file("He")
    translate_to!(molecule, [26.0, -0.2, 12.])
    apply_periodic_boundary_condition!(molecule, box)
    @test isapprox(molecule.center_of_mass, [1.0, 24.8, 12.0])
    @test isapprox(molecule.ljspheres[1].x, [1.0, 24.8, 12.0])

    translation_old_molecule_stored_properly = true
    translation_coords_changed = true
    translation_inside_box = true
    molecules = [read_molecule_file("He"), read_molecule_file("He")]
    translate_to!(molecules[1], sim_box.f_to_c * [0.99, 0.99, 0.01])
    translate_to!(molecules[2], sim_box.f_to_c * [0.99, 0.99, 0.01])
    old_molecule = translate_molecule!(molecules[1], sim_box)
    if ! isapprox(old_molecule, molecules[2]) # constructed to be identitical!
        translation_old_molecule_stored_properly = false
    end
    if isapprox(molecules[1], molecules[2])
        translation_coords_changed = false
    end

    for i = 1:100000
        which_molecule = rand(1:2) # choose molecule to move
        old_molecule_should_be = deepcopy(molecules[which_molecule])
        old_molecule = translate_molecule!(molecules[which_molecule], sim_box)
        if ! isapprox(old_molecule, old_molecule_should_be)
            translation_coords_changed = false
        end
        if outside_box(molecules[which_molecule], sim_box)
            translation_inside_box = false
        end
    end
    @test translation_old_molecule_stored_properly
    @test translation_coords_changed
    @test translation_inside_box
end
@printf("------------------------------\n")
@testset "Guest-guest Energetics Tests" begin
    repfactors = (1, 2, 3)
    dxf = [1.2, -0.2, 2.4]
    nearest_image!(dxf, repfactors)
    @test isapprox(dxf, [0.2, -0.2, -0.6])
    
    dxf = [1.2 -0.3; -0.2 -0.1; 2.4 3.4]
    nearest_image!(dxf, repfactors)
    @test isapprox(dxf, [0.2 -0.3; -0.2 -0.1; -0.6 0.4])

    sim_box = construct_box(25.0, 25.0, 25.0, π/2, π/2, π/2)
    # a He and Xe a distance of 6.0 away
    xe = read_molecule_file("Xe")
    he = read_molecule_file("He")
    translate_to!(xe, [5.0, 12.0, 12.0])
    translate_to!(he, [11.0, 12.0, 12.0])
    molecules = [xe, he]
    r² = (11.0 - 5.0) ^ 2 # duh
    energy = lennard_jones(r², ljforcefield.σ²[:Xe][:He], ljforcefield.ϵ[:Xe][:He])
    @test energy ≈ vdw_energy(1, molecules, ljforcefield, sim_box)
    @test energy ≈ vdw_energy(2, molecules, ljforcefield, sim_box) # symmetry
    
    # via PBC, a distance (24.0 - 5.0) > (1+5)
    translate_to!(molecules[2], [24.0, 12.0, 12.0])
    r² = (1.0 + 5.0) ^ 2 # PBC
    energy = lennard_jones(r², ljforcefield.σ²[:Xe][:He], ljforcefield.ϵ[:He][:Xe])
    @test energy ≈ vdw_energy(2, molecules, ljforcefield, sim_box)
    @test energy ≈ vdw_energy(1, molecules, ljforcefield, sim_box) # symmetry again.
    
    # put a molecule on top of first one.
    push!(molecules, deepcopy(molecules[1]))
    @test vdw_energy(2, molecules, ljforcefield, sim_box) ≈ 2 * energy
    
    @test vdw_energy(1, molecules, ljforcefield, sim_box) == Inf
    @test vdw_energy(3, molecules, ljforcefield, sim_box) == Inf
    
    # interaction energy between first and second should be same via PBC
    molecules_a = [read_molecule_file("Xe"), read_molecule_file("He")]
    translate_to!(molecules_a[1], [11.0, 1.0, 12.0])
    translate_to!(molecules_a[2], [11.0, 4.0, 12.0])
    molecules_b = [read_molecule_file("Xe"), read_molecule_file("He")]
    translate_to!(molecules_b[1], [11.0, 1.0, 12.0])
    translate_to!(molecules_b[2], [11.0, 23.0, 12.0])
    @test vdw_energy(1, molecules_a, ljforcefield, sim_box) ≈ vdw_energy(1, molecules_b, ljforcefield, sim_box)
    
    # another PBC one where three coords are different.
    molecules = [read_molecule_file("Xe"), read_molecule_file("He")]
    translate_to!(molecules[1], [24.0, 23.0, 11.0])
    translate_to!(molecules[2], [22.0, 2.0, 12.0])
    r² = 4.0^2 + 2.0^2 + 1.0^2
    energy = lennard_jones(r², ljforcefield.σ²[:He][:Xe], ljforcefield.ϵ[:He][:Xe])
    @test vdw_energy(1, molecules, ljforcefield, sim_box) ≈ energy
    @test vdw_energy(2, molecules, ljforcefield, sim_box) ≈ energy

    # test cutoff radius. molecules here are too far to interact
    translate_to!(molecules[1], [0.0, 0.0, 0.0])
    translate_to!(molecules[2], [12.0, 12.0, 12.0])
    @test vdw_energy(1, molecules, ljforcefield, sim_box) ≈ 0.0
    @test vdw_energy(2, molecules, ljforcefield, sim_box) ≈ 0.0
    # the position of a molecule should not change inside vdw_energy.
    @test all(molecules[1].ljspheres[1].x .== [0.0, 0.0, 0.0])
    @test all(molecules[2].ljspheres[1].x .== [12.0, 12.0, 12.0])
    # TODO write tests for CO2 where there are more than one beads
end
