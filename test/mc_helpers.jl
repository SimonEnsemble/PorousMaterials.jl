module MCHelpers_Test

using PorousMaterials
using OffsetArrays
using LinearAlgebra
using Test
using JLD2
using Statistics
using Random

SBMOF1 = Crystal("SBMOF-1.cif")
He = Molecule("He")
CO2 = Molecule("CO2")

@testset "MCHelpers Tests" begin
    box = SBMOF1.box

    ###
    #   random insertions
    ###
    insertion_inside_box = true # make sure insertions inside box
    insertion_at_random_coords = true # make sure coords are different for each molecule
    insertion_adds_molecule = true # make sure molecule vector increases in length

    molecules = Array{Molecule{Frac}}(undef, 0)

    m = deepcopy(He)
    for i = 1:100
        random_insertion!(molecules, box, m)
        if ! inside(molecules[i])
            insertion_inside_box = false
        end
        if ! (length(molecules) == i)
            insertion_adds_molecule = false
        end
        if i > 1
            for j = 1:(i-1)
                if isapprox(molecules[j], molecules[i])
                    # by chance this could fail but highly unlikely!
                    insertion_at_random_coords = false
                end
            end
        end
    end
    @test insertion_inside_box
    @test insertion_at_random_coords
    @test insertion_adds_molecule

    ###
    #   particle deletions
    ###
    molecule_24 = deepcopy(molecules[24])
    molecule_26 = deepcopy(molecules[26])
    remove_molecule!(25, molecules)
    @test length(molecules) == 99
    @test isapprox(molecules[24], molecule_24)
    @test isapprox(molecules[25], molecule_26)
    remove_molecule!(24, molecules)
    @test length(molecules) == 98
    @test isapprox(molecules[24], molecule_26)

    ###
    #  random translations
    ###
    adaptive_δ = AdaptiveTranslationStepSize(2.0) # default is 2 Å

    # first, test function to bring molecule inside a box.
    box = Box(25.0, 25.0, 25.0, π/2, π/2, π/2)
    molecule = deepcopy(He)
    molecule = Frac(molecule, box)
    translate_to!(molecule, Cart([26.0, -0.2, 12.]), box)
    @test ! inside(molecule)
    apply_periodic_boundary_condition!(molecule)
    @test isapprox(box.f_to_c * molecule.com.xf, [1.0, 24.8, 12.0])
    @test isapprox(box.f_to_c * molecule.atoms.coords.xf[:, 1], [1.0, 24.8, 12.0])
    @test inside(molecule)

    molecule = deepcopy(CO2)
    translate_to!(molecule, Cart([0.0, 0.0, 0.0]))
    @test isapprox(molecule.atoms.coords, molecule.charges.coords)
    box = Box(25.0, 5.0, 10.0, π/2, π/2, π/2)
    molecule = Frac(molecule, box)
    translate_to!(molecule, Cart([0.1, 7.0, -1.0]), box)
    apply_periodic_boundary_condition!(molecule)
    @test isapprox(molecule.com, Frac(Cart([0.1, 2.0, 9.0]), box))

    translation_old_molecule_stored_properly = true
    translation_coords_changed = true
    translation_inside_box = true
    molecules = Frac.([deepcopy(He), deepcopy(He)], box)
    translate_to!(molecules[1], Frac([0.99, 0.99, 0.01]))
    translate_to!(molecules[2], Cart(box.f_to_c * [0.99, 0.99, 0.01]), box)
    old_molecule = random_translation!(molecules[1], box, adaptive_δ)
    if ! isapprox(old_molecule, molecules[2]) # constructed to be identitical!
        translation_old_molecule_stored_properly = false
    end
    if isapprox(molecules[1], molecules[2])
        translation_coords_changed = false
    end

    for i = 1:100
        which_molecule = rand(1:2) # choose molecule to move
        old_molecule_should_be = deepcopy(molecules[which_molecule])
        old_molecule = random_translation!(molecules[which_molecule], box, adaptive_δ)
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
    
    # now try for a molecule with both charges and atoms.
    molecules = Molecule{Frac}[]
    box = SBMOF1.box
    m = Molecule("H2S")
    @test needs_rotations(m)
    for i = 1:10
        random_insertion!(molecules, box, m)
    end
    @test length(molecules) == 10
    other_molecules_untouched = true
    translation_old_molecule_stored_properly = true
    position_changed = true
    translation_inside_box = true
    charges_atoms_both_move_same = true
    # many random translations...
    for i = 1:5000
        which_molecule = rand(1:length(molecules))
        old_molecules = deepcopy(molecules)
        old_molecule = random_translation!(molecules[which_molecule], box, adaptive_δ)
        # do atoms and charges move by same vector?
        dx_a = molecules[which_molecule].atoms.coords.xf - old_molecules[which_molecule].atoms.coords.xf
        dx_c = molecules[which_molecule].charges.coords.xf - old_molecules[which_molecule].charges.coords.xf
        dx_c = dx_c .- dx_a[:, 1]
        dx_a = dx_a .- dx_a[:, 1]
        if ! (isapprox(dx_a, zeros(3, 3), atol=1e-8) && isapprox(dx_c, zeros(3, 3), atol=1e-8))
            charges_atoms_both_move_same = false
        end
        # is old molecule config stored properly?
        if ! isapprox(old_molecule, old_molecules[which_molecule])
            translation_old_molecule_stored_properly = false
        end
        # have the coords changed?
        if isapprox(old_molecule, molecules[which_molecule]) || isapprox(
                old_molecule.com, molecules[which_molecule].com) || isapprox(
                old_molecule.charges, molecules[which_molecule].charges) || isapprox(
                old_molecule.atoms, molecules[which_molecule].atoms)
            position_changed = false
        end
        # is it inside the box? have all other molecule been untouched?
        apply_periodic_boundary_condition!(molecules[which_molecule])
        for j = 1:length(molecules)
            if j == which_molecule
                continue
            end
            if ! isapprox(molecules[j], old_molecules[j])
                other_molecules_untouched = false
            end
        end
        for m in molecules
            if ! inside(m)
                translation_inside_box = false
            end
        end
    end
    @test other_molecules_untouched
    @test translation_old_molecule_stored_properly
    @test position_changed
    @test translation_inside_box
    @test charges_atoms_both_move_same

    ###
    #  random reinsertions
    ###
    box = Box(25.0, 25.0, 25.0, π/2, π/2, π/2)
    molecules = Frac.([deepcopy(He), deepcopy(CO2), deepcopy(He), deepcopy(CO2)], box)
    old_he = random_reinsertion!(molecules[1], box)
    old_co2 = random_reinsertion!(molecules[2], box)
    @test (inside(molecules[1]) && inside(molecules[2]))
    @test isapprox(old_he, molecules[3])
    @test isapprox(old_co2, molecules[4])
    @test ! isapprox(molecules[1].com, molecules[3].com)
    @test ! isapprox(molecules[2].com, molecules[4].com)
    @test ! isapprox(molecules[2].charges.coords, molecules[4].charges.coords)
    @test ! isapprox(molecules[2].atoms.coords, molecules[4].atoms.coords)
end
end
