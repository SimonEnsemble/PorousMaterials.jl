module Molecule_Test

using PorousMaterials
using OffsetArrays
using LinearAlgebra
using Test
using JLD2
using Statistics
using Random

function pairwise_atom_distances(molecule::Molecule, box::Box; verbose::Bool=false)
    bond_lengths = zeros(molecule.atoms.n, molecule.atoms.n)
    for i = 1:molecule.atoms.n
        for j = (i+1):molecule.atoms.n
            dx = 0.0
            if typeof(molecule.atoms.coords) == Frac
                dx = box.f_to_c * (molecule.atoms.coords.xf[:, i] - molecule.atoms.coords.xf[:, j])
                if verbose
                    println(dx)
                end
            else
                dx = molecule.atoms.coords.x[:, i] - molecule.atoms.coords.x[:, j]
            end
            bond_lengths[i, j] = norm(dx)
            bond_lengths[j, i] = bond_lengths[i, j]
        end
    end
    return bond_lengths
end

function pairwise_charge_distances(molecule::Molecule, box::Box)
    bond_lengths = zeros(molecule.charges.n, molecule.charges.n)
    for i = 1:molecule.charges.n
        for j = (i+1):molecule.charges.n
            dx = 0.0
            if typeof(molecule.charges.coords) == Frac
                dx = box.f_to_c * (molecule.charges.coords.xf[:, i] - molecule.charges.coords.xf[:, j])
            else
                dx = molecule.charges.coords.x[:, i] - molecule.charges.coords.x[:, j]
            end
            bond_lengths[i, j] = norm(dx)
            bond_lengths[j, i] = bond_lengths[i, j]
        end
    end
    return bond_lengths
end

@testset "Molecules Tests" begin
    molecule = Molecule("CO2")
    rotate!(molecule)

    # test reader
    molecule = Molecule("CO2")
    @test has_charges(molecule)
    atomic_masses = read_atomic_masses()
    @test molecule.species == :CO2
    @test molecule.atoms.n == 3
    @test molecule.atoms.species[1] == :C_CO2
    @test molecule.atoms.species[2] == :O_CO2
    @test molecule.atoms.species[3] == :O_CO2
    @test all(molecule.atoms.coords.x[:,1] .≈ [0.0, 0.0, 0.0])
    @test all(molecule.atoms.coords.x[:,2] .≈ [-1.16, 0.0, 0.0])
    @test all(molecule.atoms.coords.x[:,3] .≈ [1.16, 0.0, 0.0])
    @test all(molecule.com.x .≈ [0.0, 0.0, 0.0])
    @test molecule.charges.n == 3
    @test molecule.charges.q[1] ≈ 0.7
    @test molecule.charges.q[2] ≈ -0.35
    @test molecule.charges.q[3] ≈ -0.35
    for i = 1:3
        @test all(molecule.charges.coords.x[i] ≈ molecule.atoms.coords.x[i])
    end

    m = Molecule("CO2")
    box = Crystal("SBMOF-1.cif").box
    m = Frac(m, box)
    m = Cart(m, box)
    @test isapprox(m, Molecule("CO2")) # should restore.
    m = Frac(m, box)
    for i = 1:200
        translate_by!(m, Frac([randn(), randn(), randn()]))
        translate_by!(m, Cart([randn(), randn(), randn()]), box)
        translate_to!(m, Frac([randn(), randn(), randn()]))
        translate_to!(m, Cart([randn(), randn(), randn()]), box)
    end
    m = Cart(m, box)
    fresh_m = Molecule("CO2")
    translate_to!(fresh_m, m.com)
    @test isapprox(m, fresh_m) # should restore.

    box = unit_cube()
    m = Molecule("CO2")
    m = Frac(m, box)
    @test isapprox(Cart(m, box), Molecule("CO2"))

    # test translate_to, translate_by
    box = Crystal("SBMOF-1.cif").box
    m1 = Molecule("H2S")
    m2 = Molecule("H2S")
    m1 = Frac(m1, box)
    m2 = Frac(m2, box)
    @test isapprox(m1, m2)
    translate_by!(m2, Frac([0.0, 0.0, 0.0]))
    @test isapprox(m1, m2)
    translate_by!(m2, Frac([0.0, 1.2, 0.0]))
    @test ! isapprox(m1, m2)
    translate_to!(m2, m1.com)
    @test isapprox(m1, m2)
    translate_to!(m2, Cart([50.0, 100.0, 150.0]), box)
    @test isapprox(box.f_to_c * m2.com.xf, [50.0, 100.0, 150.0])

    translate_to!(m1, Frac([0.1, 0.2, 1.4]))
    translate_to!(m2, Cart(box.f_to_c * [0.1, 0.2, 1.4]), box)
    @test isapprox(m1, m2)
    
    translate_by!(m1, Frac([-0.1, -0.2, -1.1]))
    translate_by!(m2, Cart(box.f_to_c * [-0.1, -0.2, -1.1]), box)
    @test isapprox(m1, m2)
    rotate!(m2, box)
    rotate!(m1, box)

    m1 = Molecule("H2S")
    m2 = Molecule("H2S")
    for i = 1:1000
        if rand() > 0.5
            rotate!(m1)
        else
            rotate!(m2)
        end
    end
    pairwise_dist1 = pairwise_atom_distances(m1, box)
    pairwise_dist2 = pairwise_atom_distances(m2, box)
    ref_molecule = Molecule("H2S")
    ref_pairwise_dist = pairwise_atom_distances(ref_molecule, box)
    @test isapprox(pairwise_dist1, ref_pairwise_dist)
    @test isapprox(pairwise_dist2, ref_pairwise_dist)
    @test isapprox(m1.com.x, m2.com.x, atol=1e-12)
    
    m1 = Molecule("H2S")
    m2 = Molecule("H2S")
    box = Crystal("SBMOF-1.cif").box
    translate_to!(m1, Cart([1.0, 2.0, 3.0]))
    translate_to!(m2, Cart([1.0, 2.0, 3.0]))
    m1 = Frac(m1, box)
    m2 = Frac(m2, box)
    translate_to!(m1, Frac([3.0, 4.0, 5.0]))
    translate_to!(m2, Frac([3.0, 4.0, 5.0]))
    for i = 1:1000
        if rand() > 0.5
            rotate!(m1, box)
        else
            rotate!(m2, box)
        end
    end
    ref_molecule = Molecule("H2S")
    ref_molecule = Frac(ref_molecule, box)
    pairwise_dist1 = pairwise_atom_distances(m1, box)
    pairwise_dist2 = pairwise_atom_distances(m2, box)
    ref_pairwise_dist = pairwise_atom_distances(ref_molecule, box)
    @test isapprox(pairwise_dist1, ref_pairwise_dist, atol=1e-12)
    @test isapprox(pairwise_dist2, ref_pairwise_dist, atol=1e-12)
    @test isapprox(m1.com.xf, m2.com.xf, atol=1e-12)
    
    # bond drift
    m1 = Molecule("CO2")
    m2 = Molecule("CO2")
    bond_lengths1 = [norm(m1.atoms.coords.x[:,i]) for i = 1:m1.atoms.n]
    bond_lengths2 = [norm(m2.atoms.coords.x[:,i]) for i = 1:m2.atoms.n]
    @test isapprox(bond_lengths1, bond_lengths2, atol=1e-15)
    rotate!(m1)
    rotate!(m2)
    bond_lengths1 = [norm(m1.atoms.coords.x[:,i]) for i = 1:m1.atoms.n]
    bond_lengths2 = [norm(m2.atoms.coords.x[:,i]) for i = 1:m2.atoms.n]
    @test isapprox(bond_lengths1, bond_lengths2, atol=1e-14)
    m1.atoms.coords.x[:, 1] += [1e-10, 1e-9, 1e-10]
    bond_lengths1 = [norm(m1.atoms.coords.x[:,i]) for i = 1:m1.atoms.n]
    bond_lengths2 = [norm(m2.atoms.coords.x[:,i]) for i = 1:m2.atoms.n]
    @test ! isapprox(bond_lengths1, bond_lengths2, atol=1e-14)
    m1.atoms.coords.x[:, 1] -= [1e-10, 1e-9, 1e-10]
    bond_lengths1 = [norm(m1.atoms.coords.x[:,i]) for i = 1:m1.atoms.n]
    bond_lengths2 = [norm(m2.atoms.coords.x[:,i]) for i = 1:m2.atoms.n]
    @test isapprox(bond_lengths1, bond_lengths2, atol=1e-14)

    m1 = Molecule("Xe")
    m2 = Molecule("Xe")
    translate_to!(m1, Cart(rand(3)))
    bond_lengths1 = [norm(m1.atoms.coords.x[:,i] .- m1.com.x) for i = 1:m1.atoms.n]
    bond_lengths2 = [norm(m2.atoms.coords.x[:,i] .- m2.com.x) for i = 1:m2.atoms.n]
    @test isapprox(bond_lengths1, bond_lengths2, atol=1e-14)

    # test unit vector on sphere generator
    ms = [Molecule("He") for i = 1:10000]
    for m in ms
        translate_to!(m, Cart(rand_point_on_unit_sphere()))
    end
    @test all(isapprox.([norm(m.atoms.coords.x[:, 1]) for m in ms], 1.0))
    write_xyz(ms, unit_cube(), "random_vectors_on_sphere")
    println("See random_vectors_on_sphere")

    # Test to see if rotation_matrix() is random and uniform on sphere surface
    N = 1000000
    points = Array{Float64, 2}(undef, 3, N)
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
        if ! isapprox(r * transpose(r), Matrix{Float64}(I, 3, 3))
            r_orthogonal = false
        end
        if ! isapprox(det(r), 1.0)
            r_det_1 = false
        end
    end
    @test r_orthogonal
    @test r_det_1

    # test translate_by for fractional and cartesian
    box = Crystal("SBMOF-1.cif").box
    m1 = Molecule("CO2")
    m2 = Molecule("CO2")
    m1 = Frac(m1, box)
    m2 = Frac(m2, box)
    for i = 1:200
        translate_by!(m2, Frac([randn(), randn(), randn()]))
        translate_by!(m2, Cart([randn(), randn(), randn()]), box)
        translate_to!(m2, Frac([randn(), randn(), randn()]))
        translate_to!(m2, Cart([randn(), randn(), randn()]), box)
        rotate!(m2, box)
    end
    ref_molecule = Molecule("CO2")
    ref_pairwise_atoms_dist = pairwise_atom_distances(ref_molecule, box)
    ref_pairwise_charges_dist = pairwise_charge_distances(ref_molecule, box)
    pairwise_atoms_dist1 = pairwise_atom_distances(m1, box)
    pairwise_atoms_dist2 = pairwise_atom_distances(m2, box)
    pairwise_charges_dist1 = pairwise_charge_distances(m1, box)
    pairwise_charges_dist2 = pairwise_charge_distances(m2, box)
    @test isapprox(pairwise_atoms_dist1, ref_pairwise_atoms_dist)
    @test isapprox(pairwise_atoms_dist2, ref_pairwise_atoms_dist)
    @test isapprox(pairwise_charges_dist1, ref_pairwise_charges_dist)
    @test isapprox(pairwise_charges_dist2, ref_pairwise_charges_dist)

    # test rotate function; bond lengths must preserve, center of mass must preserve.
    box = Crystal("SBMOF-1.cif").box
    m1 = Molecule("CO2")
    m2 = Molecule("CO2")
    m1 = Frac(m1, box)
    m2 = Frac(m2, box)
    # test fractional, cartesian translates
    translate_to!(m2, Cart([5.0, 10.0, 15.0]), box)
    for i = 1:2000
        rotate!(m2, box)
    end
    @test isapprox(m2.com.xf, box.c_to_f * [5.0, 10.0, 15.0])
    @test isapprox(norm(box.f_to_c * (m2.charges.coords.xf[:, 1] - m2.charges.coords.xf[:, 2])),
                   norm(box.f_to_c * (m1.charges.coords.xf[:, 1] - m1.charges.coords.xf[:, 2])))
    @test isapprox(norm(box.f_to_c * (m2.atoms.coords.xf[:, 1] - m2.atoms.coords.xf[:, 2])),
                   norm(box.f_to_c * (m1.atoms.coords.xf[:, 1] - m1.atoms.coords.xf[:, 2])))
    m2_old = deepcopy(m2)
    rotate!(m2, box)
    @test ! isapprox(m2_old, m2)
    @test isapprox(m2_old.com, m2.com) # center of mass shld not change
    # visually inspect
    ms = [Molecule("CO2") for i = 1:1000]
    for m in ms
        rotate!(m)
    end
    write_xyz(ms, box, "co2s")
    println("see co2s.xyz for dist'n of rotations")

    # make sure rotation, translate does not chage bond lengths or mess up center of mass
    co2 = Molecule("CO2")
    atom_distances = pairwise_atom_distances(co2, unit_cube())
    charge_distances = pairwise_charge_distances(co2, unit_cube())
    co2 = Frac(co2, box)
    @test isapprox(atom_distances, pairwise_atom_distances(co2, box; verbose=true))
    @test isapprox(charge_distances, pairwise_charge_distances(co2, box))
    for i = 1:100000
        translate_to!(co2, Frac([rand(), rand(), rand()]))
        translate_by!(co2, Frac([randn(), randn(), randn()]))
        translate_to!(co2, Frac(4.0 * [rand(), rand(), rand()]))
        translate_by!(co2, Frac(4.0 * [rand(), rand(), rand()]))
        rotate!(co2, box)
    end
    println("atom dist ", pairwise_atom_distances(co2, box; verbose=true))
    println("charge dist ", pairwise_charge_distances(co2, box))
    println("com ", co2.com.xf)
    println("atom coords ", co2.atoms.coords.xf)
    println("charge coords ", co2.charges.coords.xf)
    @test isapprox(atom_distances, pairwise_atom_distances(co2, box), atol=1e-10)
    @test isapprox(charge_distances, pairwise_charge_distances(co2, box), atol=1e-10)
    @test isapprox(co2.com.xf, co2.atoms.coords.xf[:, 1], atol=1e-12) # should be on carbon
    #.atoms and charges shld have same coords still (this is just for CO2...
    @test all([isapprox(co2.atoms.coords.xf[:, k], co2.charges.coords.xf[:, k], atol=1e-12) for k = 1:3])
    # shld still be linear...
    co_vector1 = box.f_to_c * (co2.atoms.coords.xf[:, 2] - co2.atoms.coords.xf[:, 1])
    co_vector2 = box.f_to_c * (co2.atoms.coords.xf[:, 3] - co2.atoms.coords.xf[:, 1])
    @test isapprox(dot(co_vector1, co_vector2), -norm(co_vector1)^2, atol=1e-10)
end
end
