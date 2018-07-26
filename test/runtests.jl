#TODO split these into different files for different tests so e.g. framework in one suite is not used in another suite so it's clear what's goign on.

# Details from http://www.stochasticlifestyle.com/finalizing-julia-package-documentation-testing-coverage-publishing/
# Start Test Script
using PorousMaterials
using Base.Test
using OffsetArrays

@testset "Box Tests" begin
    framework = Framework("SBMOF-1.cif")
    @test isapprox(framework.box, Box(framework.box.f_to_c))
    @test framework.box.f_to_c * framework.box.c_to_f ≈ eye(3)
    @test isapprox(framework.box.reciprocal_lattice, 2 * π * inv(framework.box.f_to_c))
    @test isapprox(framework.box, Box(framework.box.a, framework.box.b, framework.box.c,
                                      framework.box.α, framework.box.β, framework.box.γ))
    @test isapprox(replicate(framework.box, (1, 1, 1)), framework.box)
    box = UnitCube()
    @test box.Ω ≈ 1.0
    @test isapprox(replicate(box, (3, 5, 4)), Box(3.0, 5.0, 4.0, π/2, π/2, π/2))
    @test framework.box.Ω ≈ det(framework.box.f_to_c)
    # test alternative constructor using f_to_c matrix
    @test isapprox(framework.box, Box(framework.box.f_to_c))
end

@testset "Crystal Tests" begin
    framework = Framework("test_structure2.cif")
    strip_numbers_from_atom_labels!(framework)
    @test framework.name == "test_structure2.cif"
    @test isapprox(framework.box, Box(10.0, 20.0, 30.0, 90*π/180, 45*π/180, 120*π/180))
    @test length(framework.atoms) == 2
    @test isapprox(framework.atoms[1], LJSphere(:Ca, [0.2, 0.5, 0.7]))
    @test isapprox(framework.atoms[2], LJSphere(:O, [0.6, 0.3, 0.1]))
    @test isapprox(framework.charges[1], PtCharge(1.0, [0.2, 0.5, 0.7]))
    @test isapprox(framework.charges[2], PtCharge(-1.0, [0.6, 0.3, 0.1]))
    new_frame = assign_charges(framework, Dict(:Ca => -2.0, :O => 2.0))
    @test isapprox(new_frame.charges[1], PtCharge(-2.0, [0.2, 0.5, 0.7]))
    @test isapprox(new_frame.charges[2], PtCharge(2.0, [0.6, 0.3, 0.1]))
    new_frame = assign_charges(framework, [4.0, -4.0])
    @test isapprox(new_frame.charges[1], PtCharge(4.0, [0.2, 0.5, 0.7]))
    @test isapprox(new_frame.charges[2], PtCharge(-4.0, [0.6, 0.3, 0.1]))
    @test charged(framework)
    @test chemical_formula(framework) == Dict(:Ca => 1, :O => 1)
    @test molecular_weight(framework) ≈ 15.9994 + 40.078
    # same as test_structure.cif but with overlapping atoms.
    framework2 = Framework("test_structure2B.cif", remove_overlap=true, check_charge_neutrality=false)
    strip_numbers_from_atom_labels!(framework2)
    @test all(isapprox.(framework.atoms, framework2.atoms)) && all(isapprox.(framework.charges, framework2.charges))

    # test .cif writer; write, read in, assert equal
    write_cif(framework, "data/crystals/rewritten_test_structure2.cif")
    framework_rewritten = Framework("rewritten_test_structure2.cif")
    @test isapprox(framework, framework_rewritten)

    # test .cssr reader too; test_structure2.{cif,cssr} designed to be the same.
    framework_from_cssr = Framework("test_structure2.cif")
    strip_numbers_from_atom_labels!(framework_from_cssr)
    @test isapprox(framework_from_cssr, framework, checknames=false)

    # test replicate framework
    sbmof = Framework("SBMOF-1.cif")
    replicated_sbmof = replicate(sbmof, (1, 1, 1))
    @test isapprox(sbmof, replicated_sbmof)

    repfactors = replication_factors(sbmof.box, 14.0)
    replicated_sbmof = replicate(sbmof, repfactors)
    @test replication_factors(replicated_sbmof.box, 14.0) == (1, 1, 1)
    @test isapprox(sbmof.atoms[1].xf ./ repfactors, replicated_sbmof.atoms[1].xf)
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
end;

@testset "Molecules Tests" begin
    # test reader
    molecule = Molecule("CO2")
    @test charged(molecule)
    atomic_masses = read_atomic_masses()
    @test molecule.species == :CO2
    @test length(molecule.atoms) == 3
    @test molecule.atoms[1].species == :C_CO2
    @test molecule.atoms[2].species == :O_CO2
    @test molecule.atoms[3].species == :O_CO2
    @test all(molecule.atoms[1].xf .≈ [0.0, 0.0, 0.0])
    @test all(molecule.atoms[2].xf .≈ [-1.16, 0.0, 0.0])
    @test all(molecule.atoms[3].xf .≈ [1.16, 0.0, 0.0])
    @test all(molecule.xf_com .≈ [0.0, 0.0, 0.0])
    @test length(molecule.charges) == 3
    @test molecule.charges[1].q ≈ 0.7
    @test molecule.charges[2].q ≈ -0.35
    @test molecule.charges[3].q ≈ -0.35
    for i = 1:3
        @test all(molecule.charges[i].xf ≈ molecule.atoms[i].xf)
    end
    
    m = Molecule("CO2")
    box = Framework("SBMOF-1.cif").box
    set_fractional_coords!.(m, box)
    set_fractional_coords_to_unit_cube!.(m, box)
    @test isapprox(m, Molecule("CO2")) # should restore.
    set_fractional_coords!(m, box)
    for i = 1:200
        translate_by!(m, [randn(), randn(), randn()])
        translate_by!(m, [randn(), randn(), randn()], box)
        translate_to!(m, [randn(), randn(), randn()])
        translate_to!(m, [randn(), randn(), randn()], box)
    end
    set_fractional_coords_to_unit_cube!.(m, box)
    fresh_m = Molecule("CO2")
    translate_to!(fresh_m, m.xf_com)
    @test isapprox(m, fresh_m) # should restore.

    box = UnitCube()
    m = Molecule("CO2")
    set_fractional_coords!(m, box)
    @test isapprox(m, Molecule("CO2"))

    # test translate_to, translate_by
    box = Box(1.23, 0.4, 6.0, π/2, π/2, π/2)
    ms = [Molecule("CO2") for i = 1:2]
    set_fractional_coords!.(ms, box)
    @test isapprox(ms[1], ms[2])
    translate_by!(ms[2], [0.0, 0.0, 0.0])
    @test isapprox(ms[1], ms[2])
    translate_by!(ms[2], [0.0, 1.2, 0.0])
    @test ! isapprox(ms[1], ms[2])
    translate_to!(ms[2], ms[1].xf_com)
    @test isapprox(ms[1], ms[2])
    translate_to!(ms[2], [50.0, 100.0, 150.0], box)
    @test isapprox(box.f_to_c * ms[2].xf_com, [50.0, 100.0, 150.0])
    # make sure bond lenghts are not perturbed by translate
    for i = 1:200
        translate_by!(ms[2], [randn(), randn(), randn()])
        translate_by!(ms[2], [randn(), randn(), randn()], box)
        translate_to!(ms[2], [randn(), randn(), randn()])
        translate_to!(ms[2], [randn(), randn(), randn()], box)
    end
    @test isapprox(norm(box.f_to_c * (ms[2].atoms[2].xf - ms[2].atoms[1].xf)), 
                   norm(box.f_to_c * (ms[1].atoms[2].xf - ms[1].atoms[1].xf))) # shldn't change bond lengths
    @test isapprox(norm(box.f_to_c * (ms[2].charges[2].xf - ms[2].charges[1].xf)), 
                   norm(box.f_to_c * (ms[1].charges[2].xf - ms[1].charges[1].xf))) # shldn't change bond lengths
    translate_to!(ms[1], [0.1, 0.2, 1.4])
    translate_to!(ms[2], box.f_to_c * [0.1, 0.2, 1.4], box)
    @test isapprox(ms[1], ms[2])
    @test outside_box(ms[1])
    translate_by!(ms[1], [-0.1, -0.2, -1.1])
    translate_by!(ms[2], box.f_to_c * [-0.1, -0.2, -1.1], box)
    @test isapprox(ms[1], ms[2])

    # test unit vector on sphere generator
    ms = [Molecule("He") for i = 1:10000]
    for m in ms
        translate_to!(m, rand_point_on_unit_sphere())
    end
    @test all(isapprox.([norm(m.atoms[1].xf) for m in ms], 1.0))
    write_to_xyz(ms, box, "random_vectors_on_sphere")
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
    
    # test translate_by for fractional and cartesian
    box = Framework("SBMOF-1.cif").box
    ms = [Molecule("CO2") for i = 1:2]
    set_fractional_coords!.(ms, box)
    for i = 1:200
        translate_by!(ms[2], [randn(), randn(), randn()])
        translate_by!(ms[2], [randn(), randn(), randn()], box)
        translate_to!(ms[2], [randn(), randn(), randn()])
        translate_to!(ms[2], [randn(), randn(), randn()], box)
    end
    @test isapprox(norm(box.f_to_c * (ms[2].atoms[2].xf - ms[2].atoms[1].xf)), 
                   norm(box.f_to_c * (ms[1].atoms[2].xf - ms[1].atoms[1].xf))) # shldn't change bond lengths
    @test isapprox(norm(box.f_to_c * (ms[2].charges[2].xf - ms[2].charges[1].xf)), 
                   norm(box.f_to_c * (ms[1].charges[2].xf - ms[1].charges[1].xf))) # shldn't change bond lengths

    # test fractional, cartesian translates
    ms = [Molecule("CO2") for i = 1:2]
    translate_by!(ms[2], [0.1, 0.2, 0.3]) # fractional
    translate_by!(ms[1], box.f_to_c * [0.1, 0.2, 0.3], box)
    @test isapprox(ms[1], ms[2])
    translate_to!(ms[1], [0.5, 0.6, 0.4])
    translate_to!(ms[2], box.f_to_c * [0.5, 0.6, 0.4], box)
    @test isapprox(ms[1], ms[2])

    # test translate to
    translate_to!(ms[1], [rand(), rand(), rand()])
    translate_to!(ms[1], [0.2, 0.4, 0.6])
    @test isapprox(ms[1].xf_com, [0.2, 0.4, 0.6])
    translate_to!(ms[2], box.f_to_c * [0.2, 0.4, 0.6], box)
    @test isapprox(ms[1], ms[2])
    rotate!(ms[1], box)
    rotate!(ms[2], box)
    @test ! isapprox(ms[1], ms[2])

    # test rotate function; bond lengths must preserve, center of mass must preserve.
    box = Framework("SBMOF-1.cif").box
    m1 = Molecule("CO2")
    m2 = Molecule("CO2")
    set_fractional_coords!(m1, box)
    set_fractional_coords!(m2, box)
    # test fractional, cartesian translates
    translate_to!(m2, [5.0, 10.0, 15.0], box)
    for i = 1:2000
        rotate!(m2, box)
    end
    @test isapprox(m2.xf_com, box.c_to_f * [5.0, 10.0, 15.0])
    @test isapprox(norm(box.f_to_c * (m2.charges[1].xf - m2.charges[2].xf)),
                   norm(box.f_to_c * (m1.charges[1].xf - m1.charges[2].xf)))
    @test isapprox(norm(box.f_to_c * (m2.atoms[1].xf - m2.atoms[2].xf)),
                   norm(box.f_to_c * (m1.atoms[1].xf - m1.atoms[2].xf)))
    m2_old = deepcopy(m2)
    rotate!(m2, box)
    @test ! isapprox(m2_old, m2)
    # visually inspect
    ms = [Molecule("CO2") for i = 1:1000]
    set_fractional_coords!.(ms, box)
    for m in ms
       rotate!(m, box)
    end
    write_to_xyz(ms, box, "co2s")
    println("see co2s.xyz for dist'n of rotations")

    # make sure rotation, translate does not chage bond lengths or mess up center of mass
    co2 = Molecule("CO2")
    bond_length = norm(co2.charges[1].xf - co2.charges[2].xf)
    set_fractional_coords!(co2, box)
    @test isapprox(bond_length, norm(box.f_to_c * (co2.charges[1].xf - co2.charges[2].xf)))
    for i = 1:100000
        translate_to!(co2, [rand(), rand(), rand()])
        translate_by!(co2, [randn(), randn(), randn()])
        translate_to!(co2, 4.0 * [rand(), rand(), rand()], box)
        translate_by!(co2, 4.0 * [rand(), rand(), rand()], box)
        rotate!(co2, box)
    end
    @test isapprox(norm(box.f_to_c * (co2.charges[1].xf - co2.charges[2].xf)), 
                bond_length, atol=1e-12)
    @test isapprox(norm(box.f_to_c * (co2.atoms[1].xf - co2.atoms[2].xf)), 
                bond_length, atol=1e-12)
    @test isapprox(co2.xf_com, co2.atoms[1].xf, atol=1e-12) # should be on carbon
    #.atoms and charges shld have same coords still
    @test all([isapprox(co2.atoms[k].xf, co2.charges[k].xf, atol=1e-12) for k = 1:3]) 
    # bond angles preserved.
    co_vector1 = box.f_to_c * (co2.atoms[2].xf - co2.atoms[1].xf)
    co_vector2 = box.f_to_c * (co2.atoms[3].xf - co2.atoms[1].xf)
    @test isapprox(dot(co_vector1, co_vector2), -bond_length^2, atol=1e-12)
end

@testset "NearestImage Tests" begin
    dxf = [-0.8, -0.4, 0.7]
    nearest_image!(dxf)
    @test isapprox(dxf, [0.2, -0.4, -0.3])

    dxf = [-0.3, -0.1, -0.9]
    nearest_image!(dxf)
    @test isapprox(dxf, [-0.3, -0.1, 0.1])

    box = UnitCube()
    x = [0.9, 0.1, 0.1]
    y = [0.0, 0.0, 0.0]
    @test nearest_r(x, y, box) ≈ norm([-0.1, 0.1, 0.1])
    x = [0.5, 0.5, 0.5]
    y = [0.4, 0.4, 0.35]
    @test nearest_r(x, y, box) ≈ norm(x-y)
end

@testset "Forcefield Tests" begin
    ljforcefield = LJForceField("Dreiding.csv", cutoffradius=12.5,
        mixing_rules="Lorentz-Berthelot") # Dreiding
    # test reading of force field
    @test ljforcefield.pure_σ[:He] == 1.0
    @test ljforcefield.pure_ϵ[:Zn] == 12.0
    @test ljforcefield.σ²[:Zn][:He] == ((1.0 + 3.0) / 2) ^ 2
    @test ljforcefield.ϵ[:He][:Zn] == sqrt(12.0 * 3.0)
    @test ljforcefield.ϵ[:He][:Zn] == ljforcefield.ϵ[:Zn][:He] # symmetry
    @test ljforcefield.σ²[:He][:Zn] == ljforcefield.σ²[:Zn][:He] # symmetry
    @test ljforcefield.cutoffradius_squared == 12.5 ^ 2
    
    # test calculation of replication factors required
    frame = Framework("test_structure.cif") # .cif
    strip_numbers_from_atom_labels!(frame)
    rep_factors = replication_factors(frame.box, ljforcefield)
    @test rep_factors == (25, 25, 25)

    # test check for force field coverage
    @test check_forcefield_coverage(Molecule("CO2"), ljforcefield)
    framework = Framework("SBMOF-1.cif")
    @test check_forcefield_coverage(framework, ljforcefield)
    push!(framework.atoms, LJSphere(:bogus_atom, [rand(), rand(), rand()]))
    @test !check_forcefield_coverage(framework, ljforcefield)
end;

@testset "VdwEnergetics Tests" begin
    # Xe in SBMOF-1 tests, comparing to RASPA
    ljforcefield = LJForceField("Dreiding.csv", cutoffradius=12.5, mixing_rules="Lorentz-Berthelot") # Dreiding
    sbmof1 = Framework("SBMOF-1.cif")
    rep_factors_sbmof1 = replication_factors(sbmof1.box, ljforcefield)
    sbmof1 = replicate(sbmof1, rep_factors_sbmof1)
    write_to_xyz(sbmof1, "replicated_sbmof1")
    xenon = Molecule("Xe")
    set_fractional_coords!(xenon, sbmof1.box)
    @test ! charged(xenon)
    xenon.atoms[1].xf[:] = sbmof1.box.c_to_f * zeros(3)
    energy = vdw_energy(sbmof1, xenon, ljforcefield)
    @test isapprox(energy, -5041.58, atol = 0.005)
    xenon.atoms[1].xf[:] = sbmof1.box.c_to_f * [0.494265, 2.22668, 0.450354]
    energy = vdw_energy(sbmof1, xenon, ljforcefield)
    @test isapprox(energy, 12945.838, atol = 0.005)

    # NIST data to test LJ potentials
    # data from here: https://www.nist.gov/mml/csd/chemical-informatics-research-group/lennard-jones-fluid-reference-calculations
    # created bogus atom X for this purpose.
    ljff = LJForceField("NIST.csv", cutoffradius=3.0)
    energies_should_be = [-4.3515E+03, -6.9000E+02, -1.1467E+03, -1.6790E+01]
    for c = 1:4 # four configurations
        # read in positions of atoms provided by NIST ("X" atoms)
        posfile = open("nist/lennardjones/lj_sample_config_periodic$c.txt")
        lines = readlines(posfile)
        # first line is dims of unit cell box
        dims = parse.(Float64, split(lines[1]))
        box = Box(dims..., π/2, π/2, π/2)
        # second line is # of molecules
        n = parse(Int, lines[2])

        # read in molecule positions, construct them
        ms = Molecule[]
        for i = 1:n
            xyz = split(lines[2+i])[2:end]
            x = parse.(Float64, xyz)
            m = Molecule("X")
            set_fractional_coords!(m, box)
            translate_to!(m, box.c_to_f * x)
            push!(ms, m)
        end
        close(posfile)

        # compute energy of the configuration
        energy = PorousMaterials.total_vdw_energy(ms, ljff, box)
        @test isapprox(energy, energies_should_be[c], atol=0.1)
    end
    # test vdw_energy_no_PBC, which is the vdw_energy function when no PBCs are applied.
    #  The following "framework" is a cage floating in space so no atoms are near the boundary
    #   of the unit cell box. So with cutoff should get same with or without PBCs.
    box = Box(100.0, 100.0, 100.0, π/2, π/2, π/2)
    co2 = Molecule("CO2")
    set_fractional_coords!(co2, box)
    translate_to!(co2, [50.0, 50.0, 50.0], box)
    f = Framework("cage_in_space.cif") # same cage, but shifted to [50, 50, 50] in unit cell box 100 by 100 by 100.
    ljff = LJForceField("UFF.csv")
    # energy with PBC but padded so effetive periodic interactions are zero, bc beyond cutoff
    energy = vdw_energy(f, co2, ljff)
    atoms, x = read_xyz("data/crystals/CB5.xyz") # raw .xyz of cage
    ljspheres = [LJSphere(atoms[i], x[:, i] + [50.0, 50.0, 50.0]) for i = 1:length(atoms)] # cage as LJSphere array
    co2 = Molecule("CO2")
    translate_to!(co2, [50.0, 50.0, 50.0])
    @test isapprox(energy, vdw_energy_no_PBC(co2, ljspheres, ljff))
end

@testset "Energetics_Until Tests" begin
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
end

@testset "ElectrostaticEnergetics Tests" begin
    framework = Framework("NU-1000_Greg.cif")
    sr_cutoff_r = 12.5
    rep_factors = replication_factors(framework, sr_cutoff_r)
    sim_box = replicate(framework.box, rep_factors)
    framework = replicate(framework, rep_factors)
    eparams, kvecs, eikar, eikbr, eikcr = setup_Ewald_sum(sr_cutoff_r, sim_box, verbose=false, ϵ=1e-6)
    q_test = 0.8096
    # ensure getting right Ewald settings
    #  note there are differnet method to choose
    #  these params for a givne precision so if you changed
    #  `determine_ewald_params` that may be ok if you still get the
    #  right electrostatic potential...
    @test eparams.kreps == (9, 9, 9)
    @test isapprox(eparams.α, 0.2471, atol=0.05)
    # construct box so recip. lattice is dimension (2, 10, 5)
    box = Box(0.5*2*π, 0.1*2*π, 0.2*2*π, π/2, π/2, π/2)
    @test PorousMaterials.required_kreps(box, 2.1^2) == (1, 0, 0)
    @test PorousMaterials.required_kreps(box, 5.1^2) == (2, 0, 1)
    @test PorousMaterials.required_kreps(box, 10.1^2) == (5, 1, 2)

    xf = framework.box.c_to_f * [9.535619863743, 20.685576379935, 0.127344239990]
    pc = PtCharge(q_test, xf)
    m = Molecule(:pt_charge, LJSphere[], [pc], xf)
    ϕ = electrostatic_potential_energy(framework, m, eparams, kvecs, eikar, eikbr, eikcr)
    @test isapprox(ϕ, 111373.38, atol=2.5)

    xf = framework.box.c_to_f * [4.269654927228, 23.137319129548, 28.352847101096]
    pc = PtCharge(q_test, xf)
    m = Molecule(:pt_charge, LJSphere[], [pc], xf)
    ϕ = electrostatic_potential_energy(framework, m, eparams, kvecs, eikar, eikbr, eikcr)
    @test isapprox(ϕ, -531.0, atol=0.5)

    xf = framework.box.c_to_f * [-0.047382031804, 7.209555961450, 5.158180463556]
    pc = PtCharge(q_test, xf)
    m = Molecule(:pt_charge, LJSphere[], [pc], xf)
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
        box = Box(dims..., π/2, π/2, π/2)
        # second line is # of molecules
        n = parse(Int, lines[2]) * 3 # 2H, 1O per n

        # read in molecule positions, construct them
        ms = Molecule[]
        q_H = 0.42380 # on H, -2q on O
        qs = PorousMaterials.PtCharge[]
        for i = 1:n
            if i % 3 == 1 # new water molecule
                qs = PorousMaterials.PtCharge[]
            end
            # get x position
            xyz = split(lines[2+i])[2:4]
            x = parse.(Float64, xyz)
            # get species
            O_or_H = split(lines[2+i])[end]
            q = O_or_H == "O" ? -2 * q_H : q_H
            # add to charges
            push!(qs, PorousMaterials.PtCharge(q, box.c_to_f * x))
            # construct molecule
            if i % 3 == 0
                com = [0.0, 0.0, 0.0]
                for q in qs
                   com += box.f_to_c * q.xf
                end
                com /= 3
                m = Molecule(:H2O, LJSphere[], qs, box.c_to_f * com)
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

    ###
    #  Greg Chung's ZIF-71 w/ bogus charges tests
    ###
    zif71 = Framework("zif71_bogus_charges.cif")
    strip_numbers_from_atom_labels!(zif71)
    ff = LJForceField("Greg_bogus_ZIF71.csv", cutoffradius=12.8)
    co2 = Molecule("CO2EPM2")
    set_fractional_coords!(co2, zif71.box)

    # test 1: guest-host
    @assert(co2.atoms[1].species == :C_CO2) # assumed C is first in input file...
    @assert(isapprox(co2.charges[1].q, 0.7)) # assumed C is first in input file...
    # load in coordinates of CO2 at test location
    co2.atoms[1].xf[:] = [0.50543, 0.57349, 0.50788] # C
    co2.atoms[2].xf[:] = [0.46884, 0.57393, 0.52461] # O
    co2.atoms[3].xf[:] = [0.54203, 0.57305, 0.49116] # O
    co2.charges[1].xf[:] = [0.50543, 0.57349, 0.50788] # C
    co2.charges[2].xf[:] = [0.46884, 0.57393, 0.52461] # O
    co2.charges[3].xf[:] = [0.54203, 0.57305, 0.49116] # O
    # test vdW energy
    @test isapprox(vdw_energy(zif71, co2, ff), -132.56, atol=0.01)
    # test electrostatics
    eparams, kvecs, eikar, eikbr, eikcr = setup_Ewald_sum(12.0, zif71.box, verbose=false, ϵ=1e-6)
    ϕ = electrostatic_potential_energy(zif71, co2, eparams, kvecs, eikar, eikbr, eikcr)
    @test isapprox(ϕ, -9.37846564, atol=0.1)

    # test 2: guest-guest
    co2.atoms[1].xf[:] = [0.50543, 0.57349, 0.50788]
    co2.atoms[2].xf[:] = [0.54203, 0.57305, 0.49116]
    co2.atoms[3].xf[:] = [0.46884, 0.57393, 0.52461]
    co2.charges[1].xf[:] = [0.50543, 0.57349, 0.50788]
    co2.charges[2].xf[:] = [0.54203, 0.57305, 0.49116]
    co2.charges[3].xf[:] = [0.46884, 0.57393, 0.52461]

    co2_2 = Molecule("CO2EPM2")
    set_fractional_coords!(co2_2, zif71.box)
    co2_2.atoms[1].xf[:] = [0.50680, 0.38496, 0.50788]
    co2_2.atoms[2].xf[:] = [0.54340, 0.38451, 0.49116]
    co2_2.atoms[3].xf[:] = [0.47020, 0.38540, 0.52461]
    co2_2.charges[1].xf[:] = [0.50680, 0.38496, 0.50788]
    co2_2.charges[2].xf[:] = [0.54340, 0.38451, 0.49116]
    co2_2.charges[3].xf[:] = [0.47020, 0.38540, 0.52461]
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

end

@testset "MCHelpers Tests" begin
    sim_box = Framework("SBMOF-1.cif").box

    #
    #INSERTION TESTS
    #
    insertion_inside_box = true
    insertion_at_random_coords = true
    insertion_adds_molecule = true

    molecules = Array{Molecule}(0)

    m = Molecule("He")
    set_fractional_coords!(m, sim_box)
    for i = 1:100
        insert_molecule!(molecules, sim_box, m)
        if outside_box(molecules[i])
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
    box = Box(25.0, 25.0, 25.0, π/2, π/2, π/2)
    molecule = Molecule("He")
    set_fractional_coords!(molecule, box)
    translate_to!(molecule, [26.0, -0.2, 12.], box)
    apply_periodic_boundary_condition!(molecule)
    @test isapprox(box.f_to_c * molecule.xf_com, [1.0, 24.8, 12.0])
    @test isapprox(box.f_to_c * molecule.atoms[1].xf, [1.0, 24.8, 12.0])

    translation_old_molecule_stored_properly = true
    translation_coords_changed = true
    translation_inside_box = true
    molecules = [Molecule("He"), Molecule("He")]
    set_fractional_coords!.(molecules, box)
    translate_to!(molecules[1], [0.99, 0.99, 0.01])
    translate_to!(molecules[2], box.f_to_c * [0.99, 0.99, 0.01], box)
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
        if outside_box(molecules[which_molecule])
            translation_inside_box = false
        end
    end
    @test translation_old_molecule_stored_properly
    @test translation_coords_changed
    @test translation_inside_box

    #
    #REINSERTION TESTS
    #
    box = Box(25.0, 25.0, 25.0, π/2, π/2, π/2)
    molecules = [Molecule("He"), Molecule("CO2"), Molecule("He"), Molecule("CO2")]
    set_fractional_coords!.(molecules, box)
    old_he = reinsert_molecule!(molecules[1], box)
    old_co2 = reinsert_molecule!(molecules[2], box)
    @test ! (outside_box(molecules[1]) | outside_box(molecules[2]))
    @test isapprox(old_he, molecules[3])
    @test isapprox(old_co2, molecules[4])
    @test ! isapprox(molecules[1].xf_com, molecules[3].xf_com)
    @test ! isapprox(molecules[2].xf_com, molecules[4].xf_com)
end

@testset "Guest-guest Energetics Tests" begin
    # TODO put these with VdWEnergetics
    ljforcefield = LJForceField("Dreiding.csv", cutoffradius=12.5)
    sim_box = Box(25.0, 25.0, 25.0, π/2, π/2, π/2)
    # a He and Xe a distance of 6.0 away
    xe = Molecule("Xe")
    he = Molecule("He")
    set_fractional_coords!(xe, sim_box)
    set_fractional_coords!(he, sim_box)
    translate_to!(xe, [5.0, 12.0, 12.0], sim_box)
    translate_to!(he, [11.0, 12.0, 12.0], sim_box)
    molecules = [xe, he]
    r² = (11.0 - 5.0) ^ 2 # duh
    energy = lennard_jones(r², ljforcefield.σ²[:Xe][:He], ljforcefield.ϵ[:Xe][:He])
    @test energy ≈ vdw_energy(1, molecules, ljforcefield, sim_box)
    @test energy ≈ vdw_energy(2, molecules, ljforcefield, sim_box) # symmetry

    # via PBC, a distance (24.0 - 5.0) > (1+5)
    translate_to!(molecules[2], [24.0, 12.0, 12.0], sim_box)
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
    molecules_a = [Molecule("Xe"), Molecule("He")]
    set_fractional_coords!.(molecules_a, sim_box)
    translate_to!(molecules_a[1], [11.0, 1.0, 12.0], sim_box)
    translate_to!(molecules_a[2], [11.0, 4.0, 12.0], sim_box)
    molecules_b = [Molecule("Xe"), Molecule("He")]
    set_fractional_coords!.(molecules_b, sim_box)
    translate_to!(molecules_b[1], [11.0, 1.0, 12.0], sim_box)
    translate_to!(molecules_b[2], [11.0, 23.0, 12.0], sim_box)
    @test vdw_energy(1, molecules_a, ljforcefield, sim_box) ≈ vdw_energy(1, molecules_b, ljforcefield, sim_box)

    # another PBC one where three coords are different.
    molecules = [Molecule("Xe"), Molecule("He")]
    set_fractional_coords!.(molecules, sim_box)
    translate_to!(molecules[1], [24.0, 23.0, 11.0], sim_box)
    translate_to!(molecules[2], [22.0, 2.0, 12.0], sim_box)
    r² = 4.0^2 + 2.0^2 + 1.0^2
    energy = lennard_jones(r², ljforcefield.σ²[:He][:Xe], ljforcefield.ϵ[:He][:Xe])
    @test vdw_energy(1, molecules, ljforcefield, sim_box) ≈ energy
    @test vdw_energy(2, molecules, ljforcefield, sim_box) ≈ energy

    # test cutoff radius. molecules here are too far to interact
    translate_to!(molecules[1], [0.0, 0.0, 0.0], sim_box)
    translate_to!(molecules[2], [12.0, 12.0, 12.0], sim_box)
    @test vdw_energy(1, molecules, ljforcefield, sim_box) ≈ 0.0
    @test vdw_energy(2, molecules, ljforcefield, sim_box) ≈ 0.0
    # the position of a molecule should not change inside vdw_energy.
    @test all(molecules[1].atoms[1].xf .== sim_box.c_to_f * [0.0, 0.0, 0.0])
    @test all(molecules[2].atoms[1].xf .== sim_box.c_to_f * [12.0, 12.0, 12.0])
    # TODO write tests for CO2 where there are more than one beads

    # Molecules with more than one ljsphere

    # two CO2 molecules 6.0 units apart
    molecules_co2 = [Molecule("CO2"), Molecule("CO2")]
    set_fractional_coords!.(molecules_co2, sim_box)
    translate_to!(molecules_co2[1], [12.0, 9.0, 12.0], sim_box)
    translate_to!(molecules_co2[2], [12.0, 15.0, 12.0], sim_box)
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

    # PBC placing one at 2.0 and the other at 21.0
    translate_to!(molecules_co2[1], [12.0, 2.0, 12.0], sim_box)
    translate_to!(molecules_co2[2], [12.0, 21.0, 12.0], sim_box)
    r²_com = (4.0 + 2.0)^2
    energy = (2.0 * lennard_jones(r²_com, ljforcefield.σ²[:O_CO2][:O_CO2], ljforcefield.ϵ[:O_CO2][:O_CO2])
        + 4.0 * lennard_jones(r²_com + r²_co, ljforcefield.σ²[:O_CO2][:C_CO2], ljforcefield.ϵ[:O_CO2][:C_CO2])
        + 2.0 * lennard_jones(r²_com + r²_oo, ljforcefield.σ²[:O_CO2][:O_CO2], ljforcefield.ϵ[:O_CO2][:O_CO2])
        + lennard_jones(r²_com, ljforcefield.σ²[:C_CO2][:C_CO2], ljforcefield.ϵ[:C_CO2][:C_CO2]))
    @test vdw_energy(1, molecules_co2, ljforcefield, sim_box) ≈ energy
    @test vdw_energy(2, molecules_co2, ljforcefield, sim_box) ≈ energy

    # testing cutoff radius, so only one oxygen from each will be able to interact
    # making a larger sim_box so that only a few.atoms from each CO2 will be able to interact
    sim_box_large = Box(50.0, 50.0, 50.0, π/2, π/2, π/2)
    molecules_co2 = [Molecule("CO2"), Molecule("CO2")]
    set_fractional_coords!.(molecules_co2, sim_box_large)
    # placed 12.6 units apart so the C atoms will be outside the cutoff radius,
    #   but one O atom from each will be inside, so these will interact
    translate_to!(molecules_co2[1], [0.0, 0.0, 0.0], sim_box_large)
    translate_to!(molecules_co2[2], [13.0, 0.0, 0.0], sim_box_large)
    r²_com = (13.0)^2
    r²_o = (13.0 - (2.0 * 1.16))^2
    r²_co = (13.0 - 1.16)^2
    energy = (lennard_jones(r²_o, ljforcefield.σ²[:O_CO2][:O_CO2], ljforcefield.ϵ[:O_CO2][:O_CO2])
        + 2 * lennard_jones(r²_co, ljforcefield.σ²[:O_CO2][:C_CO2], ljforcefield.ϵ[:O_CO2][:C_CO2]))
    @test vdw_energy(1, molecules_co2, ljforcefield, sim_box_large) ≈ energy
    @test vdw_energy(2, molecules_co2, ljforcefield, sim_box_large) ≈ energy
end

@testset "EOS Tests" begin
    # Peng-Robinsion EOS test for methane.
    gas = PengRobinsonGas(:CH4)
    props = calculate_properties(gas, 298.0, 65.0, verbose=false)
    @test isapprox(props["compressibility factor"], 0.874496226625811, atol=1e-4)
    @test isapprox(props["fugacity coefficient"], 0.8729028157628362, atol=1e-4)
    @test isapprox(props["fugacity (bar)"], 65.0 * 0.8729028157628362, atol=1e-4)
    @test isapprox(props["density (mol/m³)"], 3000.054418, atol=0.2)
    @test isapprox(props["molar volume (L/mol)"], 0.333327, atol=1e-4)
end
