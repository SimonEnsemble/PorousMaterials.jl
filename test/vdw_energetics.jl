module Vdw_Energetics_Test

using PorousMaterials
using OffsetArrays
using LinearAlgebra
using Test
using JLD2
using Statistics
using Random

@testset "VdwEnergetics Tests" begin
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
end
end
