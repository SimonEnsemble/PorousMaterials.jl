module Molecule_Test

using PorousMaterials
using OffsetArrays
using LinearAlgebra
using Test
using JLD2
using Statistics
using Random

function rand_point_on_unit_sphere()
    u = randn(3)
    u_norm = norm(u)
    if u_norm < 1e-6 # avoid numerical error in division
        return rand_point_on_unit_sphere()
    end
    return u / u_norm
end

function pairwise_distances(coords::Frac, box::Box)
    n = size(coords)[2]
    pad = zeros(n, n)
    for i = 1:n
        for j = 1:n
            pad[i, j] = distance(coords, box, i, j, false)
        end
    end
    return pad
end

function pairwise_distances(coords::Cart)
    n = size(coords)[2]
    pad = zeros(n, n)
    for i = 1:n
        for j = 1:n
            pad[i, j] = distance(coords, unit_cube(), i, j, false)
        end
    end
    return pad
end

pairwise_distances(m::Molecule{Cart}) = [pairwise_distances(m.atoms.coords), pairwise_distances(m.charges.coords)]
pairwise_distances(m::Molecule{Frac}, box::Box) = [pairwise_distances(m.atoms.coords, box), pairwise_distances(m.charges.coords, box)]

@testset "Molecules Tests" begin
    ###
    #   molecule file reader
    ###
    molecule = Molecule("CO2")
    @test needs_rotations(molecule)
    @test has_charges(molecule)
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
    
    ###
    #  extremely basic: make sure these functions change the molecule and don't screw up bond distances
    ###
    m = Molecule("CO2")
    pad = pairwise_distances(m)
    dx = 4 * randn(3)
    translate_by!(m, Cart(dx))
    @test ! isapprox(m.atoms, Molecule("CO2").atoms)
    @test ! isapprox(m.charges, Molecule("CO2").charges)
    @test ! isapprox(m.com, Molecule("CO2").com)
    @test isapprox(pad, pairwise_distances(m))
    translate_by!(m, Cart(-1 * dx))
    @test isapprox(m, Molecule("CO2"))
    
    m = Molecule("CO2")
    x_new_com = 4 * randn(3)
    new_com = Cart(x_new_com)
    translate_to!(m, new_com)
    @test ! isapprox(m.atoms, Molecule("CO2").atoms)
    @test ! isapprox(m.charges, Molecule("CO2").charges)
    @test ! isapprox(m.com, Molecule("CO2").com)
    @test isapprox(m.com, new_com)
    @test isapprox(pad, pairwise_distances(m))
    
    m = Molecule("CO2")
    random_rotation!(m)
    @test ! isapprox(m.atoms, Molecule("CO2").atoms)
    @test ! isapprox(m.charges, Molecule("CO2").charges)
    @test isapprox(m.com, Molecule("CO2").com)
    @test isapprox(pad, pairwise_distances(m))
        
    # ... in frac coords
    m = Molecule("H2S")
    pad = pairwise_distances(m)
    box = Crystal("SBMOF-1.cif").box
    m = Frac(m, box)
    for i = 1:200
        translate_by!(m, Frac([randn(), randn(), randn()]))
        translate_by!(m, Cart([randn(), randn(), randn()]), box)
        translate_to!(m, Frac([randn(), randn(), randn()]))
        translate_to!(m, Cart([randn(), randn(), randn()]), box)
        random_rotation!(m, box)
    end
    pad_after = pairwise_distances(m, box)
    @test isapprox(pad, pad_after)
    
    ###
    #   Frac to Cart converter
    ###
    m = Molecule("CO2")
    box = Crystal("SBMOF-1.cif").box
    m = Frac(m, box)
    m = Cart(m, box)
    @test isapprox(m, Molecule("CO2")) # should restore.
    
    box = unit_cube()
    m = Molecule("CO2")
    m = Frac(m, box)
    @test isapprox(Cart(m, box), Molecule("CO2"))
    
    ###
    # make sure translate by/to preserve orientation 
    ###
    # ... in frac coords
    m = Molecule("CO2")
    m = Frac(m, box)
    for i = 1:10
        translate_by!(m, Frac([randn(), randn(), randn()]))
        translate_by!(m, Cart([randn(), randn(), randn()]), box)
        translate_to!(m, Frac([randn(), randn(), randn()]))
        translate_to!(m, Cart([randn(), randn(), randn()]), box)
    end
    fresh_m = Molecule("CO2")
    m = Cart(m, box)
    translate_to!(fresh_m, m.com)
    @test isapprox(m, fresh_m) # should restore and preserve orientation
    
    # ... in cart coords
    m = Molecule("CO2")
    for i = 1:10
        translate_by!(m, Cart([randn(), randn(), randn()]))
        translate_by!(m, Frac([randn(), randn(), randn()]), box)
        translate_to!(m, Cart([randn(), randn(), randn()]))
        translate_to!(m, Frac([randn(), randn(), randn()]), box)
    end
    fresh_m = Molecule("CO2")
    translate_to!(fresh_m, m.com)
    @test isapprox(m, fresh_m) # should restore and preserve orientation

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
    random_rotation!(m2, box)
    random_rotation!(m1, box)
    
    ###
    # make sure random rotations preserve bond distances.
    ###
    #  ... rotations in cartesian space
    m = Molecule("H2S")
    com = Cart(10.0 * randn(3))
    translate_to!(m, com)
    pad = pairwise_distances(m)
    for i = 1:1000
        random_rotation!(m)
    end
    @test isapprox(pad, pairwise_distances(m))
    @test isapprox(m.com, com, atol=1e-12) # com should not change
    
    #  ... rotations in frac space
    m = Molecule("H2S")
    pad = pairwise_distances(m)
    box = Crystal("SBMOF-1.cif").box
    m = Frac(m, box)
    com = Frac(10.0 * randn(3))
    translate_to!(m, com)
    for i = 1:1000
        random_rotation!(m, box)
    end
    @test isapprox(pad, pairwise_distances(m, box))
    @test isapprox(m.com, com, atol=1e-12)

    # test unit vector on sphere generator
    ms = [Molecule("He") for i = 1:10000]
    for m in ms
        translate_to!(m, Cart(rand_point_on_unit_sphere()))
    end
    @test all(isapprox.([norm(m.atoms.coords.x[:, 1]) for m in ms], 1.0))
    write_xyz(ms, "random_vectors_on_sphere")
    println("See random_vectors_on_sphere")

    # Test to see if random_rotation_matrix() is random and uniform on sphere surface
    N = 1000000
    points = Array{Float64, 2}(undef, 3, N)
    for i = 1:N
        points[:, i] = random_rotation_matrix() * [0., 0., 1.]
    end

    for i = 1:3
        r = rand()
        count = zeros(10)
        for j = 1:10
            for k = 1:N
                if points[1, k] > 0 && points[2, k] ^ 2 + points[3, k] ^ 2 <= r ^ 2
                    count[j] += 1
                end
            end
            points = random_rotation_matrix() * points
        end
        @test (maximum(count) - minimum(count)) / N < 0.01
    end

    ###
    #   random rotation should place H in each half-quadrant equal amount of time...
    ###
    m = Molecule("H2S")
    translate_to!(m, origin(Cart))
    @assert m.atoms.species[2] == :H_H2S
    N = 10000000
    up = 0 # H in upper half-sphere
    left = 0 # H in left half-sphere
    back = 0 # H in back half-sphere
    for i = 1:N
        random_rotation!(m)
        if m.atoms.coords.x[3, 2] > 0.0
            up += 1
        end
        if m.atoms.coords.x[2, 2] < 0.0
            back += 1
        end
        if m.atoms.coords.x[1, 2] < 0.0
            left += 1
        end
    end
    @test isapprox(up / N, 0.5, atol=0.01)
    @test isapprox(left / N, 0.5, atol=0.01)
    @test isapprox(back / N, 0.5, atol=0.01)
    
    ###
    #   rotation matrix should be orthogonal
    ###
    r_orthogonal = true
    r_det_1 = true
    for i = 1:300
        r = random_rotation_matrix()
        if ! isapprox(r * transpose(r), Matrix{Float64}(I, 3, 3))
            r_orthogonal = false
        end
        if ! isapprox(det(r), 1.0)
            r_det_1 = false
        end
    end
    @test r_orthogonal
    @test r_det_1

    # visually inspection that rotations are random
    ms = [Molecule("CO2") for i = 1:1000]
    for m in ms
        random_rotation!(m)
    end
    write_xyz(ms, "co2s")
    @info "see co2s.xyz for dist'n of rotations"

    @test ! needs_rotations(Molecule("Xe"))
    
    # frac and cart
    m = Molecule("H2S")
    box = Box(2.0, 4.0, 8.0)
    m_f = Frac(m, box)
    @test isapprox(m_f.atoms.coords.xf[1, :], m.atoms.coords.x[1, :] / 2.0)
    @test isapprox(m_f.atoms.coords.xf[2, :], m.atoms.coords.x[2, :] / 4.0)
    @test isapprox(m_f.charges.coords.xf[1, :], m.charges.coords.x[1, :] / 2.0)
    @test isapprox(m_f.charges.coords.xf[2, :], m.charges.coords.x[2, :] / 4.0)
    m_f_c = Cart(m_f, box)
    @test isapprox(m_f_c.charges, m.charges)
    @test isapprox(m_f_c.atoms, m.atoms)

    # distorted
    adaptive_δ = AdaptiveTranslationStepSize(2.0) # default is 2 Å
    box = Crystal("SBMOF-1.cif").box
    m = Frac(Molecule("H2S"), box)
    m_ref = Frac(Molecule("H2S"), box)
    random_translation!(m, box, adaptive_δ)
    random_translation!(m_ref, box, adaptive_δ)
    random_rotation!(m, box)
    random_rotation!(m_ref, box)
    @test ! PorousMaterials.distortion(m, m_ref, box)
    m.atoms.coords.xf[:, 1] += randn(3)
    @test PorousMaterials.distortion(m, m_ref, box)

    # "center of mass" for massless molecules
    rc[:atomic_masses][:null] = 0.0
    molecule = Molecule("com_test")
    @test isapprox(molecule.com, Cart([0.;0.;0.]))
end
end
