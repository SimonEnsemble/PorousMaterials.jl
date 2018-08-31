module MCHelpers_Test

using PorousMaterials
using OffsetArrays
using LinearAlgebra
using Test
using JLD2
using Statistics
using Random

@testset "MCHelpers Tests" begin
    sim_box = Framework("SBMOF-1.cif").box

    #
    #INSERTION TESTS
    #
    insertion_inside_box = true
    insertion_at_random_coords = true
    insertion_adds_molecule = true

    molecules = Array{Molecule}(undef, 0)

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
    @test isapprox(box.f_to_c * molecule.atoms.xf[:, 1], [1.0, 24.8, 12.0])

    translation_old_molecule_stored_properly = true
    translation_coords_changed = true
    translation_inside_box = true
    molecules = [Molecule("He"), Molecule("He")]
    for molecule in molecules
        set_fractional_coords!(molecule, box)
    end
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
    for molecule in molecules
        set_fractional_coords!(molecule, box)
    end
    old_he = reinsert_molecule!(molecules[1], box)
    old_co2 = reinsert_molecule!(molecules[2], box)
    @test ! (outside_box(molecules[1]) | outside_box(molecules[2]))
    @test isapprox(old_he, molecules[3])
    @test isapprox(old_co2, molecules[4])
    @test ! isapprox(molecules[1].xf_com, molecules[3].xf_com)
    @test ! isapprox(molecules[2].xf_com, molecules[4].xf_com)
end
end
