module Guest_Guest_Energetics_Test

using PorousMaterials
using OffsetArrays
using LinearAlgebra
using Test
using JLD2
using Statistics
using Random

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
    for m in molecules_a
        set_fractional_coords!(m, sim_box)
    end
    translate_to!(molecules_a[1], [11.0, 1.0, 12.0], sim_box)
    translate_to!(molecules_a[2], [11.0, 4.0, 12.0], sim_box)
    molecules_b = [Molecule("Xe"), Molecule("He")]
    for m in molecules_b
        set_fractional_coords!(m, sim_box)
    end
    translate_to!(molecules_b[1], [11.0, 1.0, 12.0], sim_box)
    translate_to!(molecules_b[2], [11.0, 23.0, 12.0], sim_box)
    @test vdw_energy(1, molecules_a, ljforcefield, sim_box) ≈ vdw_energy(1, molecules_b, ljforcefield, sim_box)

    # another PBC one where three coords are different.
    molecules = [Molecule("Xe"), Molecule("He")]
    for m in molecules
        set_fractional_coords!(m, sim_box)
    end
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
    @test all(molecules[1].atoms.xf[:, 1] .== sim_box.c_to_f * [0.0, 0.0, 0.0])
    @test all(molecules[2].atoms.xf[:, 1] .== sim_box.c_to_f * [12.0, 12.0, 12.0])
    # TODO write tests for CO2 where there are more than one beads

    # Molecules with more than one ljsphere

    # two CO2 molecules 6.0 units apart
    molecules_co2 = [Molecule("CO2"), Molecule("CO2")]
    for m in molecules_co2
        set_fractional_coords!(m, sim_box)
    end
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
    for m in molecules_co2
        set_fractional_coords!(m, sim_box_large)
    end
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
end
