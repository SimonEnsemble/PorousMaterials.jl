module Vdw_Energetics_Test

using PorousMaterials
using OffsetArrays
using LinearAlgebra
using Test
using JLD2
using Statistics
using Random

@testset "VdwEnergetics Tests" begin
    # lennard jones function
    σ² = 1.0
    ϵ = 3.0
    @test isapprox(PorousMaterials.lennard_jones(σ², σ², ϵ), 0.0)
    @test isapprox(PorousMaterials.lennard_jones(2 ^ (2/6) * σ², σ², ϵ), -1.0 * ϵ)

    ###
    #  Xe in SBMOF-1 tests, comparing to RASPA
    ###
    ljff = LJForceField("Dreiding", r_cutoff=12.5, mixing_rules="Lorentz-Berthelot")
    xtal = Crystal("SBMOF-1.cif")
    rep_factors = replication_factors(xtal.box, ljff)
    xtal = replicate(xtal, rep_factors)
    mol = Molecule("Xe")
    mol = Frac(mol, xtal.box)
    @test ! has_charges(mol)
    # point #1
    translate_to!(mol, origin(Frac))
    energy = vdw_energy(xtal, mol, ljff)
    @test isapprox(energy, -5041.58, atol = 0.005)
    # point #2
    translate_to!(mol, Frac(Cart([0.494265, 2.22668, 0.450354]), xtal.box))
 #     xenon.atoms.xf[:, 1] = sbmof1.box.c_to_f * [0.494265, 2.22668, 0.450354]
    energy = vdw_energy(xtal, mol, ljff)
    @test isapprox(energy, 12945.838, atol = 0.005)
    
    random_rotation!(mol, xtal.box)
    energy = vdw_energy(xtal, mol, ljff)
    mol = Cart(mol, xtal.box)
    @test isapprox(energy, vdw_energy(xtal, mol, ljff))
    
    ###
    # NIST data to test LJ potentials
    # data from here: https://www.nist.gov/mml/csd/chemical-informatics-research-group/lennard-jones-fluid-reference-calculations
    # created bogus atom X for this purpose.
    ###
    ljff = LJForceField("NIST", r_cutoff=3.0)
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
            m = Frac(m, box)
            translate_to!(m, Frac(box.c_to_f * x))
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
    translate_to!(co2, Cart([50.0, 50.0, 50.0]))
    f = Crystal("cage_in_space.cif") # same cage, but shifted to [50, 50, 50] in unit cell box 100 by 100 by 100.
    ljff = LJForceField("UFF")
    # energy with PBC but padded so effetive periodic interactions are zero, bc beyond cutoff
    energy = vdw_energy(f, co2, ljff)
    atoms = read_xyz("data/crystals/CB5.xyz") # raw .xyz of cage
    co2 = Molecule("CO2")
    @test isapprox(energy, vdw_energy_no_PBC(atoms, co2.atoms, ljff))

    ###
    #   guest-guest energetics
    ###
    ljff = LJForceField("Dreiding", r_cutoff=12.5)
    box = Box(25.0, 25.0, 25.0, π/2, π/2, π/2)
    # a He and Xe a distance of 6.0 away
    xe = Frac(Molecule("Xe"), box)
    he = Frac(Molecule("He"), box)
    translate_to!(xe, Cart([5.0, 12.0, 12.0]), box)
    translate_to!(he, Cart([11.0, 12.0, 12.0]), box)
    molecules = [xe, he]
    r² = (11.0 - 5.0) ^ 2 # duh
    energy = lennard_jones(r², ljff.σ²[:Xe][:He], ljff.ϵ[:Xe][:He])
    @test energy ≈ vdw_energy(1, molecules, ljff, box)
    @test energy ≈ vdw_energy(2, molecules, ljff, box) # symmetry

    # via PBC, a distance (24.0 - 5.0) > (1+5)
    translate_to!(molecules[2], Cart([24.0, 12.0, 12.0]), box)
    r² = (1.0 + 5.0) ^ 2 # PBC
    energy = lennard_jones(r², ljff.σ²[:Xe][:He], ljff.ϵ[:He][:Xe])
    @test energy ≈ vdw_energy(2, molecules, ljff, box)
    @test energy ≈ vdw_energy(1, molecules, ljff, box) # symmetry again.

    # put a molecule on top of first one.
    push!(molecules, deepcopy(molecules[1]))
    @test vdw_energy(2, molecules, ljff, box) ≈ 2 * energy

    @test vdw_energy(1, molecules, ljff, box) == Inf
    @test vdw_energy(3, molecules, ljff, box) == Inf

    # interaction energy between first and second should be same via PBC
    molecules_a = Frac.([Molecule("Xe"), Molecule("He")], box)
    translate_to!(molecules_a[1], Cart([11.0, 1.0, 12.0]), box)
    translate_to!(molecules_a[2], Cart([11.0, 4.0, 12.0]), box)
    molecules_b = Frac.([Molecule("Xe"), Molecule("He")], box)
    translate_to!(molecules_b[1], Cart([11.0, 1.0, 12.0]), box)
    translate_to!(molecules_b[2], Cart([11.0, 23.0, 12.0]), box)
    @test vdw_energy(1, molecules_a, ljff, box) ≈ vdw_energy(1, molecules_b, ljff, box)

    # another PBC one where three coords are different.
    molecules = Frac.([Molecule("Xe"), Molecule("He")], box)
    translate_to!(molecules[1], Cart([24.0, 23.0, 11.0]), box)
    translate_to!(molecules[2], Cart([22.0, 2.0, 12.0]), box)
    r² = 4.0^2 + 2.0^2 + 1.0^2
    energy = lennard_jones(r², ljff.σ²[:He][:Xe], ljff.ϵ[:He][:Xe])
    @test vdw_energy(1, molecules, ljff, box) ≈ energy
    @test vdw_energy(2, molecules, ljff, box) ≈ energy

    # test cutoff radius. molecules here are too far to interact
    translate_to!(molecules[1], Cart([0.0, 0.0, 0.0]), box)
    translate_to!(molecules[2], Cart([12.0, 12.0, 12.0]), box)
    @test vdw_energy(1, molecules, ljff, box) ≈ 0.0
    @test vdw_energy(2, molecules, ljff, box) ≈ 0.0
    # the position of a molecule should not change inside vdw_energy.
    @test all(molecules[1].atoms.coords.xf[:, 1] .== box.c_to_f * [0.0, 0.0, 0.0])
    @test all(molecules[2].atoms.coords.xf[:, 1] .== box.c_to_f * [12.0, 12.0, 12.0])
    # TODO write tests for CO2 where there are more than one beads

    # Molecules with more than one ljsphere

    # two CO2 molecules 6.0 units apart
    molecules_co2 = Frac.([Molecule("CO2"), Molecule("CO2")], box)
    translate_to!(molecules_co2[1], Cart([12.0, 9.0, 12.0]), box)
    translate_to!(molecules_co2[2], Cart([12.0, 15.0, 12.0]), box)
    # because the molecules have not been rotated, all corresponding beads are same
    #   distance apart when they are separated along the y-axis
    r²_com = (15.0 - 9.0)^2
    # distance between teh central carbon and an oxygen in one molecule this
    #   takes advantage of the fact that the carbon is the central atom, and that
    #   all three atoms are in a line
    r²_co = 1.16^2
    # distance between the two oxygens in one molecule
    r²_oo = (2.0 * 1.16)^2
    energy = (2.0 * lennard_jones(r²_com, ljff.σ²[:O_CO2][:O_CO2], ljff.ϵ[:O_CO2][:O_CO2])
        + 4.0 * lennard_jones(r²_com + r²_co, ljff.σ²[:O_CO2][:C_CO2], ljff.ϵ[:O_CO2][:C_CO2])
        + 2.0 * lennard_jones(r²_com + r²_oo, ljff.σ²[:O_CO2][:O_CO2], ljff.ϵ[:O_CO2][:O_CO2])
        + lennard_jones(r²_com, ljff.σ²[:C_CO2][:C_CO2], ljff.ϵ[:C_CO2][:C_CO2]))
    @test vdw_energy(1, molecules_co2, ljff, box) ≈ energy
    @test vdw_energy(2, molecules_co2, ljff, box) ≈ energy

    # PBC placing one at 2.0 and the other at 21.0
    translate_to!(molecules_co2[1], Cart([12.0, 2.0, 12.0]), box)
    translate_to!(molecules_co2[2], Cart([12.0, 21.0, 12.0]), box)
    r²_com = (4.0 + 2.0)^2
    energy = (2.0 * lennard_jones(r²_com, ljff.σ²[:O_CO2][:O_CO2], ljff.ϵ[:O_CO2][:O_CO2])
        + 4.0 * lennard_jones(r²_com + r²_co, ljff.σ²[:O_CO2][:C_CO2], ljff.ϵ[:O_CO2][:C_CO2])
        + 2.0 * lennard_jones(r²_com + r²_oo, ljff.σ²[:O_CO2][:O_CO2], ljff.ϵ[:O_CO2][:O_CO2])
        + lennard_jones(r²_com, ljff.σ²[:C_CO2][:C_CO2], ljff.ϵ[:C_CO2][:C_CO2]))
    @test vdw_energy(1, molecules_co2, ljff, box) ≈ energy
    @test vdw_energy(2, molecules_co2, ljff, box) ≈ energy

    # testing cutoff radius, so only one oxygen from each will be able to interact
    # making a larger box so that only a few.atoms from each CO2 will be able to interact
    box_large = Box(50.0, 50.0, 50.0, π/2, π/2, π/2)
    molecules_co2 = Frac.([Molecule("CO2"), Molecule("CO2")], box_large)
    # placed 12.6 units apart so the C atoms will be outside the cutoff radius,
    #   but one O atom from each will be inside, so these will interact
    translate_to!(molecules_co2[1], Cart([0.0, 0.0, 0.0]), box_large)
    translate_to!(molecules_co2[2], Cart([13.0, 0.0, 0.0]), box_large)
    r²_com = (13.0)^2
    r²_o = (13.0 - (2.0 * 1.16))^2
    r²_co = (13.0 - 1.16)^2
    energy = (lennard_jones(r²_o, ljff.σ²[:O_CO2][:O_CO2], ljff.ϵ[:O_CO2][:O_CO2])
        + 2 * lennard_jones(r²_co, ljff.σ²[:O_CO2][:C_CO2], ljff.ϵ[:O_CO2][:C_CO2]))
    @test vdw_energy(1, molecules_co2, ljff, box_large) ≈ energy
    @test vdw_energy(2, molecules_co2, ljff, box_large) ≈ energy


    ###
    #  Mixture Tests
    ###
    # test the vdW guest-guest interation between two species
    # set up the system
    molecules = [[Molecule("Xe"), Molecule("Xe")], [Molecule("Kr")]]
    ljff      = LJForceField("UFF", r_cutoff=5.1)
    box       = Box(10.0, 10.0, 10.0)
    # convert molecules array to fractional using this box.
    molecules = [Frac.(mols, box) for mols in molecules]
    # position the molecules
    translate_to!( molecules[1][1], Frac([1.0, 1.0, 0.1]) )
    translate_to!( molecules[1][2], Frac([1.0, 5.0, 0.1]) )
    translate_to!( molecules[2][1], Frac([5.0, 1.0, 0.1]) )
    # calculate vdW_energy interaction
    r12 = sum( (molecules[1][2].atoms.coords.xf - molecules[1][1].atoms.coords.xf) .^ 2) # 16.0
    r13 = sum( (molecules[1][2].atoms.coords.xf - molecules[1][1].atoms.coords.xf) .^ 2) # 16.0
    energy_12 = 4.0 * ljff.ϵ[:Xe][:Xe] * ((ljff.σ²[:Xe][:Xe] / (r12)) ^ 6 - (ljff.σ²[:Xe][:Xe] / (r12)) ^ 3)
    energy_13 = 4.0 * ljff.ϵ[:Xe][:Kr] * ((ljff.σ²[:Xe][:Kr] / (r13)) ^ 6 - (ljff.σ²[:Xe][:Kr] / (r13)) ^ 3)
    @test (energy_12 + energy_13) ≈ vdw_energy(1, 1, molecules, ljff, box)
end
end
