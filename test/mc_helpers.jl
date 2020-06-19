module MCHelpers_Test

using PorousMaterials
using OffsetArrays
using LinearAlgebra
using Test
using JLD2
using Statistics
using Random

@testset "MCHelpers Tests" begin
    sim_box = Crystal("SBMOF-1.cif").box

    #
    #INSERTION TESTS
    #
    insertion_inside_box = true
    insertion_at_random_coords = true
    insertion_adds_molecule = true

    molecules = Array{Molecule{Frac}}(undef, 0)

    m = Molecule("He")
    for i = 1:100
        random_insertion!(molecules, sim_box, m)
        if ! inside(molecules[i])
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
        remove_molecule!(rand(1:length(molecules)), molecules)
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
    molecule = Frac(molecule, box)
    translate_to!(molecule, Cart([26.0, -0.2, 12.]), box)
    @test ! inside(molecule)
    apply_periodic_boundary_condition!(molecule)
    @test isapprox(box.f_to_c * molecule.com.xf, [1.0, 24.8, 12.0])
    @test isapprox(box.f_to_c * molecule.atoms.coords.xf[:, 1], [1.0, 24.8, 12.0])
    @test inside(molecule)

    translation_old_molecule_stored_properly = true
    translation_coords_changed = true
    translation_inside_box = true
    molecules = Frac.([Molecule("He"), Molecule("He")], box)
    translate_to!(molecules[1], Frac([0.99, 0.99, 0.01]))
    translate_to!(molecules[2], Cart(box.f_to_c * [0.99, 0.99, 0.01]), box)
    old_molecule = random_translation!(molecules[1], sim_box)
    if ! isapprox(old_molecule, molecules[2]) # constructed to be identitical!
        translation_old_molecule_stored_properly = false
    end
    if isapprox(molecules[1], molecules[2])
        translation_coords_changed = false
    end

    for i = 1:100000
        which_molecule = rand(1:2) # choose molecule to move
        old_molecule_should_be = deepcopy(molecules[which_molecule])
        old_molecule = random_translation!(molecules[which_molecule], sim_box)
        if ! isapprox(old_molecule, old_molecule_should_be)
            translation_coords_changed = false
        end
        if ! inside(molecules[which_molecule])
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
    molecules = Frac.([Molecule("He"), Molecule("CO2"), Molecule("He"), Molecule("CO2")], box)
    old_he = random_reinsertion!(molecules[1], box)
    old_co2 = random_reinsertion!(molecules[2], box)
    @test (inside(molecules[1]) && inside(molecules[2]))
    @test isapprox(old_he, molecules[3])
    @test isapprox(old_co2, molecules[4])
    @test ! isapprox(molecules[1].com, molecules[3].com)
    @test ! isapprox(molecules[2].com, molecules[4].com)
end
end
