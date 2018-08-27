module Electrostatics_Energetics_Test

using PorousMaterials
using OffsetArrays
using LinearAlgebra
using Test
using JLD2
using Statistics
using Random

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
                @assert (PorousMaterials.total_charge(m) == 0.0)
                push!(ms, m)
                @assert (isapprox(PorousMaterials.total_charge(m), 0.0, rtol=0.001))
            end
        end
        @assert (length(ms) == n/3)
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
        eikar = OffsetArray{Complex{Float64}}(undef, 0:kreps[1])
        eikbr = OffsetArray{Complex{Float64}}(undef, -kreps[2]:kreps[2])
        eikcr = OffsetArray{Complex{Float64}}(undef, -kreps[3]:kreps[3])
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
    @assert (co2.atoms[1].species == :C_CO2) # assumed C is first in input file...
    @assert (isapprox(co2.charges[1].q, 0.7)) # assumed C is first in input file...
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
end
