module Molecule_Test

using PorousMaterials
using OffsetArrays
using LinearAlgebra
using Test
using JLD2
using Statistics
using Random

@testset "Molecules Tests" begin
    molecule = Molecule("CO2")
    rotate!(molecule, UnitCube())
    @test isapprox(pairwise_atom_distances(molecule, UnitCube()),
        [0 1.16 1.16; 1.16  0.0   1.16*2; 1.16 1.16*2 0])
    @test isapprox(pairwise_charge_distances(molecule, UnitCube()),
        [0 1.16 1.16; 1.16  0.0   1.16*2; 1.16 1.16*2 0])

    # test reader
    molecule = Molecule("CO2")
    @test charged(molecule)
    atomic_masses = read_atomic_masses()
    @test molecule.species == :CO2
    @test molecule.atoms.n_atoms == 3
    @test molecule.atoms.species[1] == :C_CO2
    @test molecule.atoms.species[2] == :O_CO2
    @test molecule.atoms.species[3] == :O_CO2
    @test all(molecule.atoms.xf[:, 1] .≈ [0.0, 0.0, 0.0])
    @test all(molecule.atoms.xf[:, 2] .≈ [-1.16, 0.0, 0.0])
    @test all(molecule.atoms.xf[:, 3] .≈ [1.16, 0.0, 0.0])
    @test all(molecule.xf_com .≈ [0.0, 0.0, 0.0])
    @test molecule.charges.n_charges == 3
    @test molecule.charges.q[1] ≈ 0.7
    @test molecule.charges.q[2] ≈ -0.35
    @test molecule.charges.q[3] ≈ -0.35
    for i = 1:3
        @test all(molecule.charges.xf[:, i] ≈ molecule.atoms.xf[:, i])
    end

    m = Molecule("CO2")
    box = Framework("SBMOF-1.cif").box
    set_fractional_coords!(m, box)
    set_fractional_coords_to_unit_cube!(m, box)
    @test isapprox(m, Molecule("CO2")) # should restore.
    set_fractional_coords!(m, box)
    for i = 1:200
        translate_by!(m, [randn(), randn(), randn()])
        translate_by!(m, [randn(), randn(), randn()], box)
        translate_to!(m, [randn(), randn(), randn()])
        translate_to!(m, [randn(), randn(), randn()], box)
    end
    set_fractional_coords_to_unit_cube!(m, box)
    fresh_m = Molecule("CO2")
    translate_to!(fresh_m, m.xf_com)
    @test isapprox(m, fresh_m) # should restore.

    box = UnitCube()
    m = Molecule("CO2")
    set_fractional_coords!(m, box)
    @test isapprox(m, Molecule("CO2"))

    # test translate_to, translate_by
    box = Framework("SBMOF-1.cif").box
    ms = [Molecule("H2S") for i = 1:2]
    for m in ms
        set_fractional_coords!(m, box)
    end
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
    # bond lengths preserved?
    for m in ms
        @test isapprox(pairwise_atom_distances(m, box),
            pairwise_atom_distances(Molecule("H2S"), UnitCube()))
        @test isapprox(pairwise_charge_distances(m, box),
            pairwise_charge_distances(Molecule("H2S"), UnitCube()))
    end
    translate_to!(ms[1], [0.1, 0.2, 1.4])
    translate_to!(ms[2], box.f_to_c * [0.1, 0.2, 1.4], box)
    @test isapprox(ms[1], ms[2])
    @test outside_box(ms[1])
    translate_by!(ms[1], [-0.1, -0.2, -1.1])
    translate_by!(ms[2], box.f_to_c * [-0.1, -0.2, -1.1], box)
    @test isapprox(ms[1], ms[2])
    rotate!(ms[2], box)
    rotate!(ms[1], box)
    @test isapprox(pairwise_atom_distances(ms[2], box),
                   pairwise_atom_distances(Molecule("H2S"), UnitCube()))
    @test isapprox(pairwise_charge_distances(ms[2], box),
                   pairwise_charge_distances(Molecule("H2S"), UnitCube()))

    ms = [Molecule("H2S") for i = 1:2]
    for i = 1:1000
        rotate!(ms[rand(1:2)], UnitCube())
    end
    for m in ms
        @test ! bond_length_drift(m, Molecule("H2S"), UnitCube(), atol=1e-12)
    end
    @test isapprox(ms[1].xf_com, ms[2].xf_com, atol=1e-12)
    
    ms = [Molecule("H2S") for i = 1:2]
    box = Framework("SBMOF-1.cif").box
    for m in ms
        translate_to!(m, [1.0, 2.0, 3.0])
        set_fractional_coords!(m, box)
        translate_to!(m, [3.0, 4.0, 5.0])
    end
    for i = 1:1000
        rotate!(ms[rand(1:2)], box)
    end
    ref_molecule = Molecule("H2S")
    set_fractional_coords!(ref_molecule, box)
    for m in ms
        @test ! bond_length_drift(m, ref_molecule, box, atol=1e-12)
    end
    @test isapprox(ms[1].xf_com, ms[2].xf_com, atol=1e-12)
    
    # bond drift
    m1 = Molecule("CO2")
    m2 = Molecule("CO2")
    @test ! bond_length_drift(m1, m2, UnitCube(), atol=1e-15)
    rotate!(m1, UnitCube())
    rotate!(m2, UnitCube())
    @test ! bond_length_drift(m1, m2, UnitCube(), atol=1e-14)
    m1.atoms.xf[:, 1] += [1e-10, 1e-9, 1e-10]
    @test bond_length_drift(m1, m2, UnitCube(), atol=1e-14, throw_warning=false)
    m1.atoms.xf[:, 1] -= [1e-10, 1e-9, 1e-10]
    @test ! bond_length_drift(m1, m2, UnitCube(), atol=1e-14)
    m2.charges.xf[:, 1] += [1e-10, 1e-9, 1e-10]
    @test bond_length_drift(m1, m2, UnitCube(), atol=1e-14, throw_warning=false)
    m1 = Molecule("Xe")
    m2 = Molecule("Xe")
    translate_to!(m1, rand(3))
    @test ! bond_length_drift(m1, m2, UnitCube(), atol=1e-14)

    # test unit vector on sphere generator
    ms = [Molecule("He") for i = 1:10000]
    for m in ms
        translate_to!(m, rand_point_on_unit_sphere(), UnitCube())
    end
    @test all(isapprox.([norm(m.atoms.xf[:, 1]) for m in ms], 1.0))
    write_xyz(ms, UnitCube(), "random_vectors_on_sphere")
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
    box = Framework("SBMOF-1.cif").box
    ms = [Molecule("CO2") for i = 1:2]
    for m in ms
        set_fractional_coords!(m, box)
    end
    for i = 1:200
        translate_by!(ms[2], [randn(), randn(), randn()])
        translate_by!(ms[2], [randn(), randn(), randn()], box)
        translate_to!(ms[2], [randn(), randn(), randn()])
        translate_to!(ms[2], [randn(), randn(), randn()], box)
        rotate!(ms[2], box)
    end
    for m in ms
        @test isapprox(pairwise_atom_distances(m,               box),
                       pairwise_atom_distances(Molecule("CO2"), UnitCube()))
        @test isapprox(pairwise_charge_distances(m,               box),
                       pairwise_charge_distances(Molecule("CO2"), UnitCube()))
    end
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
    @test isapprox(norm(box.f_to_c * (m2.charges.xf[:, 1] - m2.charges.xf[:, 2])),
                   norm(box.f_to_c * (m1.charges.xf[:, 1] - m1.charges.xf[:, 2])))
    @test isapprox(norm(box.f_to_c * (m2.atoms.xf[:, 1] - m2.atoms.xf[:, 2])),
                   norm(box.f_to_c * (m1.atoms.xf[:, 1] - m1.atoms.xf[:, 2])))
    m2_old = deepcopy(m2)
    rotate!(m2, box)
    @test ! isapprox(m2_old, m2)
    @test isapprox(m2_old.xf_com, m2.xf_com) # center of mass shld not change
    # visually inspect
    ms = [Molecule("CO2") for i = 1:1000]
    for m in ms
        set_fractional_coords!(m, box)
    end
    for m in ms
       rotate!(m, box)
    end
    write_xyz(ms, box, "co2s")
    println("see co2s.xyz for dist'n of rotations")

    # make sure rotation, translate does not chage bond lengths or mess up center of mass
    co2 = Molecule("CO2")
    atom_distances = pairwise_atom_distances(co2, UnitCube())
    charge_distances = pairwise_charge_distances(co2, UnitCube())
    set_fractional_coords!(co2, box)
    @test isapprox(atom_distances, pairwise_atom_distances(co2, box))
    @test isapprox(charge_distances, pairwise_charge_distances(co2, box))
    for i = 1:100000
        translate_to!(co2, [rand(), rand(), rand()])
        translate_by!(co2, [randn(), randn(), randn()])
        translate_to!(co2, 4.0 * [rand(), rand(), rand()], box)
        translate_by!(co2, 4.0 * [rand(), rand(), rand()], box)
        rotate!(co2, box)
    end
    println("atom dist ", pairwise_atom_distances(co2, box))
    println("charge dist ", pairwise_charge_distances(co2, box))
    @test isapprox(atom_distances, pairwise_atom_distances(co2, box), atol=1e-10)
    @test isapprox(charge_distances, pairwise_charge_distances(co2, box), atol=1e-10)
    @test isapprox(co2.xf_com, co2.atoms.xf[:, 1], atol=1e-12) # should be on carbon
    #.atoms and charges shld have same coords still (this is just for CO2...
    @test all([isapprox(co2.atoms.xf[:, k], co2.charges.xf[:, k], atol=1e-12) for k = 1:3])
    # shld still be linear...
    co_vector1 = box.f_to_c * (co2.atoms.xf[:, 2] - co2.atoms.xf[:, 1])
    co_vector2 = box.f_to_c * (co2.atoms.xf[:, 3] - co2.atoms.xf[:, 1])
    @test isapprox(dot(co_vector1, co_vector2), -norm(co_vector1)^2, atol=1e-10)
end
end
