#!/usr/bin/env julia

#TODO split these into different files for different tests so e.g. framework in one suite is not used in another suite so it's clear what's goign on.

# Details from http://www.stochasticlifestyle.com/finalizing-julia-package-documentation-testing-coverage-publishing/
# Start Test Script
using PorousMaterials
using Base.Test
using OffsetArrays

# Run Tests

@printf("\n------------------------------\nTesting Crystal.jl\n\n")
framework = read_crystal_structure_file("test_structure2.cif")
framework2 = read_crystal_structure_file("test_structure2B.cif", remove_overlap = true)
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
    @test framework.atoms == framework2.atoms && framework.xf == framework2.xf && framework.charges == framework2.charges

    # test .cif writer; write, read in, assert equal
    write_cif(framework, "data/crystals/rewritten_test_structure2.cif")
    framework_rewritten = read_crystal_structure_file("rewritten_test_structure2.cif")
    @test isapprox(framework, framework_rewritten)

    # test .cssr reader too; test_structure2.{cif,cssr} designed to be the same.
    framework_from_cssr = read_crystal_structure_file("test_structure2.cif")
    strip_numbers_from_atom_labels!(framework_from_cssr)
    @test isapprox(framework_from_cssr, framework, checknames=false)

    # test replicate framework
    sbmof = read_crystal_structure_file("SBMOF-1.cif")
    replicated_sbmof = replicate(sbmof, (1, 1, 1))
    @test isapprox(sbmof, replicated_sbmof)
    repfactors = replication_factors(sbmof.box, 14.0)
    replicated_sbmof = replicate(sbmof, repfactors)
    @test replication_factors(replicated_sbmof.box, 14.0) == (1, 1, 1)
    @test isapprox(sbmof.box.f_to_c * sbmof.xf[:, 1], replicated_sbmof.box.f_to_c * replicated_sbmof.xf[:, 1])
    @test isapprox(transpose(replicated_sbmof.box.reciprocal_lattice), 2 * π * inv(replicated_sbmof.box.f_to_c))
    @test chemical_formula(sbmof) == chemical_formula(replicated_sbmof)
    @test isapprox(crystal_density(sbmof), crystal_density(replicated_sbmof), atol=1e-7)
end;

@printf("\n------------------------------\nTesting Forcefield.jl\n\n")
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

@printf("\n------------------------------\nTesting Molecules.jl\n\n")
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

    # Test to see if rotation_matrix() is random and uniform on sphere surface
    N = 1000000
    points = Array{Float64, 2}(3,N)
    for i = 1:N
        points[:,i] = rotation_matrix() * [0., 0., 1.]
    end

    for i = 1:3
        r = rand()
        count = zeros(10)
        for j = 1:10
            for k = 1:N
                if points[1,k] > 0 && points[2,k]^2 + points[3,k]^2 <= r^2
                    count[j] += 1
                end
            end
            points = rotation_matrix() * points
        end
        @test (maximum(count) - minimum(count)) / N < 0.01
    end

    # rotation matrix should be orthogonal
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

    # make sure rotation, translate does not chage bond lengths or mess up center of mass
    co2 = read_molecule_file("CO2")
    bond_length = norm(co2.charges[1].x - co2.charges[2].x)
    for i = 1:100000
        translate_to!(co2, 10.0 * [rand(), rand(), rand()])
        translate_by!(co2, 5.0 * [randn(), randn(), randn()])
        rotate!(co2)
    end
    @test isapprox(norm(co2.charges[1].x - co2.charges[2].x), bond_length, atol=1e-12)
    @test isapprox(norm(co2.ljspheres[1].x - co2.ljspheres[2].x), bond_length, atol=1e-12)
    @test isapprox(co2.center_of_mass, co2.ljspheres[1].x, atol=1e-12) # should be on carbon
    # ljspheres and charges shld have same coords still
    @test all([isapprox(co2.ljspheres[k].x, co2.charges[k].x, atol=1e-12) for k = 1:3])
    # bond angles preserved.
    co_vector1 = co2.ljspheres[2].x - co2.ljspheres[1].x
    co_vector2 = co2.ljspheres[3].x - co2.ljspheres[1].x
    @test isapprox(dot(co_vector1, co_vector2), -bond_length^2, atol=1e-12)
end

@printf("\n------------------------------\nTesting Energetics.jl\n\n")
@testset "Energetics Tests" begin
    # test Periodic boundary conditions
    molecule1 = read_molecule_file("He")
    molecule1.ljspheres[1].x[:] = [0.5, 0.5, 0.5]
    molecule2 = read_molecule_file("He")
    molecule2.ljspheres[1].x[:] = [0.5 + rep_factors[1], 0.5 + rep_factors[2], 0.5 + rep_factors[3]]

	@test vdw_energy(frame, molecule1, ljforcefield) ≈ vdw_energy(frame, molecule2, ljforcefield)
	@test vdw_energy(frame, molecule1, ljforcefield) ≈ 4 * ljforcefield.ϵ[:He][:Zn] * ((ljforcefield.σ²[:Zn][:He] / 0.75) ^ 6 - (ljforcefield.σ²[:Zn][:He] / 0.75) ^ 3)
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
    sbmof1 = replicate(sbmof1, rep_factors_sbmof1) # replicate so nearest image convention can be applied

    xenon = read_molecule_file("Xe")
    @test ! charged(xenon)
    xenon.ljspheres[1].x[:] = zeros(3)
    energy = vdw_energy(sbmof1, xenon, ljforcefield)
	@test isapprox(energy, -5041.58, atol = 0.005)
    xenon.ljspheres[1].x[1] = 0.494265; xenon.ljspheres[1].x[2] = 2.22668; xenon.ljspheres[1].x[3] = 0.450354;
    energy = vdw_energy(sbmof1, xenon, ljforcefield)
	@test isapprox(energy, 12945.838, atol = 0.005)

    # test PBC, rep factors are (3, 5, 2)
    xf = [0.05, 0.4, 0.02][:, :]
    xenon.ljspheres[1].x[:] = sbmof1.box.f_to_c * xf
    energy1 = vdw_energy(sbmof1, xenon, ljforcefield)
    xf = [1.05, 0.4, 0.02][:, :]
    xenon.ljspheres[1].x[:] = sbmof1.box.f_to_c * xf
    energy2 = vdw_energy(sbmof1, xenon, ljforcefield)
	@test isapprox(energy1, energy2, atol = 0.00001)
    xf = [2.05, 4.4, 1.02][:, :]
    xenon.ljspheres[1].x[:] = sbmof1.box.f_to_c * xf
    energy3 = vdw_energy(sbmof1, xenon, ljforcefield)
	@test isapprox(energy1, energy3, atol = 0.00001)
    # outside box
    xf = [4.05, 5.4, 2.02][:, :]
    xenon.ljspheres[1].x[:] = sbmof1.box.f_to_c * xf
    energy4 = vdw_energy(sbmof1, xenon, ljforcefield)
	@test isapprox(energy1, energy4, atol = 0.00001)
    xf = [-0.95, 5.4, 2.02][:, :]
    xenon.ljspheres[1].x[:] = sbmof1.box.f_to_c * xf
    energy5 = vdw_energy(sbmof1, xenon, ljforcefield)
	@test isapprox(energy1, energy5, atol = 0.00001)
    xf = [-0.95, -0.6, -0.98][:, :]
    xenon.ljspheres[1].x[:] = sbmof1.box.f_to_c * xf
    energy6 = vdw_energy(sbmof1, xenon, ljforcefield)
	@test isapprox(energy1, energy6, atol = 0.00001)

    # data types for potential energies
    u = PotentialEnergy(10.0, 30.0)
    v = PotentialEnergy(3.0, 4.0)
    @test ! isapprox(v, PotentialEnergy(3.0, 1.2), verbose=false) # isapprox
    @test isapprox(sum(v), 7.0) # sum
    @test isapprox(u + v, PotentialEnergy(13.0, 34.0)) # +
    @test isapprox(u - v, PotentialEnergy(7.0, 26.0)) # -
    @test isapprox(2.0 * v, PotentialEnergy(6.0, 8.0)) # *
    @test isapprox(v * 2.0, PotentialEnergy(6.0, 8.0)) # *
    @test isapprox(u / 2.0, PotentialEnergy(5.0, 15.0)) # / 
    @test isapprox(sqrt(PotentialEnergy(4.0, 16.0)), PotentialEnergy(2.0, 4.0)) # sqrt
    @test isapprox(PorousMaterials.square(PotentialEnergy(2.0, 4.0)), PotentialEnergy(4.0, 16.0)) # square
    t = PotentialEnergy(1.0, 2.0)
    s = PotentialEnergy(300.0, 100.0)
    us = SystemPotentialEnergy(u, v)
    vs = SystemPotentialEnergy(s, t)
    @test isapprox(sum(vs), 403.0) # sum
    @test isapprox(us - vs, SystemPotentialEnergy(u - s, v - t)) # - 
    @test isapprox(us + vs, SystemPotentialEnergy(u + s, v + t)) # - 
    @test isapprox(2.0 * us, SystemPotentialEnergy(2.0 * u, 2.0 * v)) # *
    @test isapprox(2.0 * us, SystemPotentialEnergy(2.0 * u, 2.0 * v)) # *
    @test isapprox(us * 2.0, SystemPotentialEnergy(2.0 * u, 2.0 * v)) # *
    @test isapprox(us / 2.0, SystemPotentialEnergy(u / 2.0, v / 2.0)) # /
    @test isapprox(sqrt(us), SystemPotentialEnergy(sqrt(u), sqrt(v))) # sqrt
    @test isapprox(PorousMaterials.square(us), SystemPotentialEnergy(PorousMaterials.square(u), PorousMaterials.square(v))) # square

    # NIST data to test LJ potentials
    # data from here: https://www.nist.gov/mml/csd/chemical-informatics-research-group/lennard-jones-fluid-reference-calculations
    # created bogus atom X for this purpose.
    ljff = read_forcefield_file("NIST.csv", cutoffradius=3.0)
    energies_should_be = [-4.3515E+03, -6.9000E+02, -1.1467E+03, -1.6790E+01]
    for c = 1:4 # four configurations
        # read in positions of atoms provided by NIST ("X" atoms)
        posfile = open("nist/lennardjones/lj_sample_config_periodic$c.txt")
        lines = readlines(posfile)
        # first line is dims of unit cell box
        dims = parse.(Float64, split(lines[1]))
        box = construct_box(dims..., π/2, π/2, π/2)
        # second line is # of molecules
        n = parse(Int, lines[2])

        # read in molecule positions, construct them
        ms = Molecule[]
        for i = 1:n
            xyz = split(lines[2+i])[2:end]
            x = parse.(Float64, xyz)
            m = read_molecule_file("X")
            translate_to!(m, x)
            push!(ms, m)
        end
        close(posfile)

        # compute energy of the configuration
        energy = PorousMaterials.total_vdw_energy(ms, ljff, box)
        @test isapprox(energy, energies_should_be[c], atol=0.1)
    end

    ###
    #  Greg Chung's ZIF-71 w/ bogus charges tests
    ###
    zif71 = read_crystal_structure_file("zif71_bogus_charges.cif")
    strip_numbers_from_atom_labels!(zif71)
    ff = read_forcefield_file("Greg_bogus_ZIF71.csv", cutoffradius=12.8)
    co2 = read_molecule_file("CO2EPM2")

    # test 1: guest-host
    @assert(co2.ljspheres[1].atom == :C_CO2) # assumed C is first in input file...
    @assert(isapprox(co2.charges[1].q, 0.7)) # assumed C is first in input file...
    # load in coordinates of CO2 at test location
    co2.ljspheres[1].x[:] =  zif71.box.f_to_c * [0.50543, 0.57349, 0.50788] # C
    co2.ljspheres[2].x[:] =  zif71.box.f_to_c * [0.46884, 0.57393, 0.52461] # O
    co2.ljspheres[3].x[:] =  zif71.box.f_to_c * [0.54203, 0.57305, 0.49116] # O
    co2.charges[1].x[:] =  zif71.box.f_to_c * [0.50543, 0.57349, 0.50788] # C
    co2.charges[2].x[:] =  zif71.box.f_to_c * [0.46884, 0.57393, 0.52461] # O
    co2.charges[3].x[:] =  zif71.box.f_to_c * [0.54203, 0.57305, 0.49116] # O
    # test vdW energy
	@test isapprox(vdw_energy(zif71, co2, ff), -132.56, atol=0.01)
    # test electrostatics
    eparams, kvecs, eikar, eikbr, eikcr = setup_Ewald_sum(12.0, zif71.box, verbose=false, ϵ=1e-6)
    ϕ = electrostatic_potential_energy(zif71, co2, eparams, kvecs, eikar, eikbr, eikcr)
    @test isapprox(ϕ, -9.37846564, atol=0.1)

    # test 2: guest-guest
    co2.ljspheres[1].x[:] =  zif71.box.f_to_c * [0.50543, 0.57349, 0.50788]
    co2.ljspheres[2].x[:] =  zif71.box.f_to_c * [0.54203, 0.57305, 0.49116]
    co2.ljspheres[3].x[:] =  zif71.box.f_to_c * [0.46884, 0.57393, 0.52461]
    co2.charges[1].x[:] =  zif71.box.f_to_c * [0.50543, 0.57349, 0.50788]
    co2.charges[2].x[:] =  zif71.box.f_to_c * [0.54203, 0.57305, 0.49116]
    co2.charges[3].x[:] =  zif71.box.f_to_c * [0.46884, 0.57393, 0.52461]

    co2_2 = read_molecule_file("CO2EPM2")
    co2_2.ljspheres[1].x[:] =  zif71.box.f_to_c * [0.50680, 0.38496, 0.50788]
    co2_2.ljspheres[2].x[:] =  zif71.box.f_to_c * [0.54340, 0.38451, 0.49116]
    co2_2.ljspheres[3].x[:] =  zif71.box.f_to_c * [0.47020, 0.38540, 0.52461]
    co2_2.charges[1].x[:] =  zif71.box.f_to_c * [0.50680, 0.38496, 0.50788]
    co2_2.charges[2].x[:] =  zif71.box.f_to_c * [0.54340, 0.38451, 0.49116]
    co2_2.charges[3].x[:] =  zif71.box.f_to_c * [0.47020, 0.38540, 0.52461]
    @test isapprox(PorousMaterials.total_vdw_energy(zif71, [co2, co2_2], ff), -311.10392551, atol=0.1)
    @test isapprox(PorousMaterials.total_vdw_energy([co2, co2_2], ff, zif71.box), -50.975, atol=0.1)
    @test isapprox(PorousMaterials.total_electrostatic_potential_energy(zif71, [co2, co2_2], eparams, kvecs, eikar, eikbr, eikcr), -36.00, atol=0.3)
    ϕ = electrostatic_potential_energy([co2, co2_2], eparams, kvecs, eikar, eikbr, eikcr)
    @test isapprox(total(ϕ), 59.3973, atol=0.05)

    # MC function for electrostatic potential should be equal to difference in two systems
    #  one with both co2's and one with just one of the co2's
    ϕ_both = electrostatic_potential_energy([co2, co2_2], eparams, kvecs, eikar, eikbr, eikcr)
    ϕ_one = electrostatic_potential_energy([co2], eparams, kvecs, eikar, eikbr, eikcr)
    ϕ_for_MC = electrostatic_potential_energy([co2, co2_2], 2, eparams, kvecs, eikar, eikbr, eikcr)
    @test isapprox(total(ϕ_both) - total(ϕ_one), total(ϕ_for_MC))
    # difference in potential energy when adding 1 co2 should be that of the system with one co2
    ϕ_for_MC = electrostatic_potential_energy([co2], 1, eparams, kvecs, eikar, eikbr, eikcr)
    @test isapprox(total(ϕ_for_MC), total(ϕ_one))

    # assert total_electrostatic_potential function incrementer works
    @test isapprox(total(PorousMaterials.total_electrostatic_potential_energy([co2, co2_2], eparams, kvecs, eikar, eikbr, eikcr)),
                   total(electrostatic_potential_energy([co2, co2_2], eparams, kvecs, eikar, eikbr, eikcr)))

    # test vdw_energy_no_PBC, which is the vdw_energy function when no PBCs are applied.
    #  The following "framework" is a cage floating in space so no atoms are near the boundary
    #   of the unit cell box. So with cutoff should get same with or without PBCs.
    co2 = read_molecule_file("CO2")
    translate_to!(co2, [50.0, 50.0, 50.0])
    atoms, x = read_xyz("data/crystals/CB5.xyz") # raw .xyz of cage
    f = read_crystal_structure_file("cage_in_space.cif") # same cage, but shifted to [50, 50, 50] in unit cell box 100 by 100 by 100.
    ljff = read_forcefield_file("UFF.csv")
    @test isapprox(vdw_energy(f, co2, ljff), vdw_energy_no_PBC(co2, atoms, x .+ [50.0, 50.0, 50.0], ljff))

end;
#@printf("------------------------------\n")

# TODO energetics and electrostatics can be grouped together...
@printf("------------------------------\nTesting Electrostatics\n\n")
framework = read_crystal_structure_file("NU-1000_Greg.cif")

 # kreps = (11, 11, 9)
 # α = 0.265058
sr_cutoff_r = 12.5
rep_factors = replication_factors(framework, sr_cutoff_r)
sim_box = replicate(framework.box, rep_factors)
framework = replicate(framework, rep_factors)
eparams, kvecs, eikar, eikbr, eikcr = setup_Ewald_sum(sr_cutoff_r, sim_box, verbose=false, ϵ=1e-6)

q_test = 0.8096

@testset "Ewald summation Tests" begin
    # ensure getting right Ewald settings
    #  note there are differnet method to choose
    #  these params for a givne precision so if you changed
    #  `determine_ewald_params` that may be ok if you still get the
    #  right electrostatic potential...
    @test eparams.kreps == (9, 9, 9)
    @test isapprox(eparams.α, 0.2471, atol=0.05)
    # construct box so recip. lattice is dimension (2, 10, 5)
    box = construct_box(0.5*2*π, 0.1*2*π, 0.2*2*π, π/2, π/2, π/2)
    @test PorousMaterials.required_kreps(box, 2.1^2) == (1, 0, 0)
    @test PorousMaterials.required_kreps(box, 5.1^2) == (2, 0, 1)
    @test PorousMaterials.required_kreps(box, 10.1^2) == (5, 1, 2)

    x = [9.535619863743, 20.685576379935, 0.127344239990]
    pc = PointCharge(q_test, x)
    m = Molecule(:pt_charge, LennardJonesSphere[], [pc], x)
    ϕ = electrostatic_potential_energy(framework, m, eparams, kvecs, eikar, eikbr, eikcr)
    @test isapprox(ϕ, 111373.38, atol=2.5)

    x = [4.269654927228, 23.137319129548, 28.352847101096]
    pc = PointCharge(q_test, x)
    m = Molecule(:pt_charge, LennardJonesSphere[], [pc], x)
    ϕ = electrostatic_potential_energy(framework, m, eparams, kvecs, eikar, eikbr, eikcr)
    @test isapprox(ϕ, -531.0, atol=0.5)

    x = [-0.047382031804, 7.209555961450, 5.158180463556]
    pc = PointCharge(q_test, x)
    m = Molecule(:pt_charge, LennardJonesSphere[], [pc], x)
    ϕ = electrostatic_potential_energy(framework, m, eparams, kvecs, eikar, eikbr, eikcr)
    @test isapprox(ϕ, -2676.8230141, atol=0.5)

    # NIST data to test Ewald sums
    # data from here:  https://www.nist.gov/mml/csd/chemical-informatics-research-group/spce-water-reference-calculations-10%C3%A5-cutoff
    # what the energies should be for all four configurations provided by NIST
    energies_should_be = [Dict(["real"=> -5.58889e05, "fourier"=> 6.27009e03, "self"=> -2.84469e06, "intra" => 2.80999e06]),
                          Dict(["real"=> -1.19295e06, "fourier"=> 6.03495e03, "self"=> -5.68938e06, "intra" => 5.61998e06]),
                          Dict(["real"=> -1.96297e06, "fourier"=> 5.24461e03, "self"=> -8.53407e06, "intra" => 8.42998e06]),
                          Dict(["real"=> -3.57226e06, "fourier"=> 7.58785e03, "self"=> -1.42235e07, "intra" => 1.41483e07])]
    # loop over all four configurations provided by NIST
    for c = 1:length(energies_should_be)
        # read in positions of atoms provided by NIST ("X" atoms)
        posfile = open("nist/electrostatics/spce_sample_config_periodic$c.txt")
        lines = readlines(posfile)
        # first line is dims of unit cell box
        dims = parse.(Float64, split(lines[1]))
        box = construct_box(dims..., π/2, π/2, π/2)
        # second line is # of molecules
        n = parse(Int, lines[2]) * 3 # 2H, 1O per n

        # read in molecule positions, construct them
        ms = Molecule[]
        q_H = 0.42380 # on H, -2q on O
        qs = PorousMaterials.PointCharge[]
        for i = 1:n
            if i % 3 == 1 # new water molecule
                qs = PorousMaterials.PointCharge[]
            end
            # get x position
            xyz = split(lines[2+i])[2:4]
            x = parse.(Float64, xyz)
            # get species
            O_or_H = split(lines[2+i])[end]
            q = O_or_H == "O" ? -2 * q_H : q_H
            # add to charges
            push!(qs, PorousMaterials.PointCharge(q, x))
            # construct molecule
            if i % 3 == 0
                com = [0.0, 0.0, 0.0]
                for q in qs
                   com += q.x
                end
                com /= 3
                m = Molecule(:H2O, [], qs, com)
                @assert(PorousMaterials.total_charge(m) == 0.0)
                push!(ms, m)
                @assert(isapprox(PorousMaterials.total_charge(m), 0.0, rtol=0.001))
            end
        end
        @assert(length(ms) == n/3)
        close(posfile)

        # compute energy of the configuration
        sr_cutoff_r = 10.0
        # use PorousMaterials.jl settings
     #     eparams, kvecs, eikar, eikbr, eikcr = setup_Ewald_sum(sr_cutoff_r, box, verbose=true, ϵ=1e-6)
        # use NIST reported settings
        kreps = (5, 5, 5)
        eparams = PorousMaterials.EwaldParams(kreps, 5.6/box.a, sr_cutoff_r, box)
        kvecs = PorousMaterials.precompute_kvec_wts(eparams)
        # only include kvecs with <27 acc to NIST website
        kvec_keep = [kvec.ka ^ 2 + kvec.kb ^2 + kvec.kc ^2 < 27 for kvec in kvecs]
        kvecs = kvecs[kvec_keep]
        eikar = OffsetArray(Complex{Float64}, 0:kreps[1])
        eikbr = OffsetArray(Complex{Float64}, -kreps[2]:kreps[2])
        eikcr = OffsetArray(Complex{Float64}, -kreps[3]:kreps[3])
        ϕ = electrostatic_potential_energy(ms, eparams, kvecs, eikar, eikbr, eikcr)
        @test isapprox(ϕ.self, energies_should_be[c]["self"], rtol=0.00001)
        @test isapprox(ϕ.sr, energies_should_be[c]["real"], rtol=0.00001)
        @test isapprox(ϕ.lr + ϕ.lr_own_images, energies_should_be[c]["fourier"],  rtol=0.00001)
        @test isapprox(ϕ.intra, energies_should_be[c]["intra"],  rtol=0.0001)
        # test incremented potential energy
        @test isapprox(total(PorousMaterials.total_electrostatic_potential_energy(ms, eparams, kvecs, eikar, eikbr, eikcr)), total(ϕ), atol=0.01)
        # test potential energy function for MC sims
        ϕ_minus_molecule = electrostatic_potential_energy(ms[2:end], eparams, kvecs, eikar, eikbr, eikcr)
        ϕ_MC =  electrostatic_potential_energy(ms, 1, eparams, kvecs, eikar, eikbr, eikcr)
        @test(isapprox(total(ϕ) - total(ϕ_minus_molecule), total(ϕ_MC), atol=0.01))
    end
end
@printf("------------------------------\n")

@printf("\n------------------------------\nTesting GCMC.jl\n\n")
@testset "Monte Carlo Functions Tests" begin
    # replicating the unit cell to construct simulation box
    sbmof1 = read_crystal_structure_file("SBMOF-1.cif")
    sim_box = replicate(sbmof1.box, (1, 1, 1))
    @test isapprox(sim_box, sbmof1.box)
    sim_box = replicate(sbmof1.box, (2, 3, 4))
    @test sim_box.Ω ≈ sbmof1.box.Ω * 2 * 3 * 4
    @test all(sim_box.c_to_f * sbmof1.box.f_to_c * [1.0, 1.0, 1.0] .≈ [1/2, 1/3, 1/4])

    #
    #INSERTION TESTS
    #
    insertion_inside_box = true
    insertion_at_random_coords = true
    insertion_adds_molecule = true

    molecules = Array{Molecule}(0)
    # TODO rm these?
    repfactors = replication_factors(frame.box, ljforcefield)
    sim_box = replicate(frame.box, repfactors)

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

    #
    #REINSERTION TESTS
    #
    box = construct_box(25.0, 25.0, 25.0, π/2, π/2, π/2)
    molecules = [read_molecule_file("He"), read_molecule_file("CO2")]
    old_he = reinsert_molecule!(molecules[1], box)
    old_co2 = reinsert_molecule!(molecules[2], box)
    @test isapprox(old_he, read_molecule_file("He"))
    @test isapprox(old_co2, read_molecule_file("CO2"))
    @test ! isapprox(molecules[1].center_of_mass, read_molecule_file("He").center_of_mass)
    @test ! isapprox(molecules[2].center_of_mass, read_molecule_file("CO2").center_of_mass)
end
@printf("------------------------------\n")
@testset "Guest-guest Energetics Tests" begin
    dxf = [-0.8, -0.4, 0.7]
    nearest_image!(dxf)
    @test isapprox(dxf, [0.2, -0.4, -0.3])

    dxf = [0.8 -0.3; -0.2 -0.1; 0.9 -0.9]
    nearest_image!(dxf)
    @test isapprox(dxf, [-0.2 -0.3; -0.2 -0.1; -0.1 0.1])

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

    # Molecules with more than one ljsphere

    # two CO2 molecules 6.0 units apart
    molecules_co2 = [read_molecule_file("CO2"), read_molecule_file("CO2")]
    translate_to!(molecules_co2[1], [12.0, 9.0, 12.0])
    translate_to!(molecules_co2[2], [12.0, 15.0, 12.0])
    # because the molecules have not been rotated, all corresponding beads are same
    #   distance apart when they are separated along the y-axis
    r²_com = (15.0 - 9.0)^2
    # distance between teh central carbon and an oxygen in one molecule this
    #   takes advantage of the fact that the carbon is the central atom, and that
    #   all three atoms are in a line
    r²_co = 1.16^2
    # distance between the two oxygens in one molecule
    r²_oo = (2.0 * 1.16)^2
    energy = (2.0 * lennard_jones(r²_com, ljforcefield.σ²[:O_CO2][:O_CO2], ljforcefield.ϵ[:O_CO2][:O_CO2])
        + 4.0 * lennard_jones(r²_com + r²_co, ljforcefield.σ²[:O_CO2][:C_CO2], ljforcefield.ϵ[:O_CO2][:C_CO2])
        + 2.0 * lennard_jones(r²_com + r²_oo, ljforcefield.σ²[:O_CO2][:O_CO2], ljforcefield.ϵ[:O_CO2][:O_CO2])
        + lennard_jones(r²_com, ljforcefield.σ²[:C_CO2][:C_CO2], ljforcefield.ϵ[:C_CO2][:C_CO2]))
    @test vdw_energy(1, molecules_co2, ljforcefield, sim_box) ≈ energy
    @test vdw_energy(2, molecules_co2, ljforcefield, sim_box) ≈ energy
    # TODO this didn't work when the symbols were regular (:O and :C) but it did work when they were :C_CO2 and :O_CO2

    # PBC placing one at 2.0 and the other at 21.0
    translate_to!(molecules_co2[1], [12.0, 2.0, 12.0])
    translate_to!(molecules_co2[2], [12.0, 21.0, 12.0])
    r²_com = (4.0 + 2.0)^2
    energy = (2.0 * lennard_jones(r²_com, ljforcefield.σ²[:O_CO2][:O_CO2], ljforcefield.ϵ[:O_CO2][:O_CO2])
        + 4.0 * lennard_jones(r²_com + r²_co, ljforcefield.σ²[:O_CO2][:C_CO2], ljforcefield.ϵ[:O_CO2][:C_CO2])
        + 2.0 * lennard_jones(r²_com + r²_oo, ljforcefield.σ²[:O_CO2][:O_CO2], ljforcefield.ϵ[:O_CO2][:O_CO2])
        + lennard_jones(r²_com, ljforcefield.σ²[:C_CO2][:C_CO2], ljforcefield.ϵ[:C_CO2][:C_CO2]))
    @test vdw_energy(1, molecules_co2, ljforcefield, sim_box) ≈ energy
    @test vdw_energy(2, molecules_co2, ljforcefield, sim_box) ≈ energy

    # testing cutoff radius, so only one oxygen from each will be able to interact
    # making a larger sim_box so that only a few ljspheres from each CO2 will be able to interact
    sim_box_large = construct_box(50.0, 50.0, 50.0, π/2, π/2, π/2)
    # placed 12.6 units apart so the C atoms will be outside the cutoff radius,
    #   but one O atom from each will be inside, so these will interact
    translate_to!(molecules_co2[1], [0.0, 0.0, 0.0])
    translate_to!(molecules_co2[2], [13.0, 0.0, 0.0])
    r²_com = (13.0)^2
    r²_o = (13.0 - (2.0 * 1.16))^2
    r²_co = (13.0 - 1.16)^2
    energy = (lennard_jones(r²_o, ljforcefield.σ²[:O_CO2][:O_CO2], ljforcefield.ϵ[:O_CO2][:O_CO2])
        + 2 * lennard_jones(r²_co, ljforcefield.σ²[:O_CO2][:C_CO2], ljforcefield.ϵ[:O_CO2][:C_CO2]))
    @test vdw_energy(1, molecules_co2, ljforcefield, sim_box_large) ≈ energy
    @test vdw_energy(2, molecules_co2, ljforcefield, sim_box_large) ≈ energy
end

@printf("------------------------------\n")
@testset "Equation of State (EOS.jl) Tests" begin
    # Peng-Robinsion EOS test for methane.
    gas = PengRobinsonGas(:CH4)
    props = calculate_properties(gas, 298.0, 65.0)
    @test isapprox(props["compressibility factor"], 0.874910311935475, atol=1e-4)
    @test isapprox(props["fugacity coefficient"], 0.8732765620181274, atol=1e-4)
    @test isapprox(props["fugacity (bar)"], 65.0 * 0.8732765620181274, atol=1e-4)
    @test isapprox(props["density (mol/m³)"], 2998.634526, atol=0.2)
    @test isapprox(props["molar volume (L/mol)"], 0.333485, atol=1e-4)

    #Van der Waals EOS test for Hydrogen
    gas = VDWGasGas(:Hydrogen)
    props = calculate_properties(gas, 300.0, 65.0)
    @test isapprox(props["compressibility factor"], 1.0462045 , atol=0.001)
    @test isapprox(props["fugacity coefficient"], 67.98166172 / 65.0 , atol=1e-4)
    @test isapprox(props["fugacity (bar)"], 67.98166172 , atol=1e-4)
    @test isapprox(props["density (mol/m³)"], 2490.815 , atol=0.2)
    @test isapprox(props["molar volume (L/mol)"], 0.401475 , atol=1e-4)

end
