# NIST data to test LJ potentials
# data from here: https://www.nist.gov/mml/csd/chemical-informatics-research-group/lennard-jones-fluid-reference-calculations
# created bogus molecule X for this purpose.
using PorousMaterials
using OffsetArrays
ljff = read_forcefield_file("NIST.csv", cutoffradius=3.0)

energies_should_be = [-4.3515E+03, -6.9000E+02, -1.1467E+03, -1.6790E+01]

for c = 1:4
    # read in positions of atoms provided by NIST ("X" atoms)
    posfile = open("nist/lennardjones/lj_sample_config_periodic$c.txt")
    lines = readlines(posfile)
    # first line is dims of unit cell box
    dims = parse.(Float64, split(lines[1]))
    box = construct_box(dims..., π/2, π/2, π/2)
    # second line is # of molecules
    n = parse(Int, lines[2])

    # read in molecule positions, construct them
    ms = Molecule[]
    for i = 1:n
        xyz = split(lines[2+i])[2:end]
        x = parse.(Float64, xyz)
        m = read_molecule_file("X")
        translate_to!(m, x)
        push!(ms, m)
    end
    close(posfile)
    
    # compute energy of the configuration
    repfactors = (1, 1, 1)
    energy = 0.0
    for i = 1:length(ms)
        energy += vdw_energy(i, ms, ljff, box)
    end
    energy /= 2
    @assert(isapprox(energy, energies_should_be[c], atol=1.0))
end
println("LJ tests pass")

# NIST data to test Ewald sums
# data from here:  https://www.nist.gov/mml/csd/chemical-informatics-research-group/spce-water-reference-calculations-10%C3%A5-cutoff
energies_should_be = [-5.58889E+05+6.27009E+03+-2.84469E+06+2.80999E+06,
                     -1.19295E+06+6.03495E+03+-5.68938E+06+5.61998E+06,
                     -1.96297E+06+5.24461E+03+-8.53407E+06+8.42998E+06,
                     -3.57226E+06+7.58785E+03+-1.42235E+07+1.41483E+07]

for c = 1:4
    # read in positions of atoms provided by NIST ("X" atoms)
    posfile = open("nist/electrostatics/spce_sample_config_periodic$c.txt")
    lines = readlines(posfile)
    # first line is dims of unit cell box
    dims = parse.(Float64, split(lines[1]))
    box = construct_box(dims..., π/2, π/2, π/2)
    # second line is # of molecules
    n = parse(Int, lines[2]) * 3 # 2H, 1O per n

    # read in molecule positions, construct them
    ms = Molecule[]
    q_H = 0.42380 # on H, -2q on O
    qs = PorousMaterials.PointCharge[]
    for i = 1:n
        if i % 3 == 1 # new water molecule
            qs = PorousMaterials.PointCharge[]
        end
        # get x position
        xyz = split(lines[2+i])[2:4]
        x = parse.(Float64, xyz)
        # get species
        O_or_H = split(lines[2+i])[end]
        q = O_or_H == "O" ? -2 * q_H : q_H
        # add to charges
        push!(qs, PorousMaterials.PointCharge(q, x))
        # construct molecule
        if i % 3 == 0
            com = [0.0, 0.0, 0.0]
            for q in qs
               com += q.x
            end
            com /= 3
            m = Molecule(:H2O, [], qs, com)
            push!(ms, m)
            @assert(isapprox(PorousMaterials.total_charge(m), 0.0, rtol=0.001))
        end
    end
    @assert(length(ms) == n/3)
    close(posfile)
    
    # compute energy of the configuration
    sr_cutoff_r = 10.0
    # use PorousMaterials.jl settings
    eparams, kvecs, eikar, eikbr, eikcr = setup_Ewald_sum(sr_cutoff_r, box, verbose=true, ϵ=1e-6)
    # use NIST reported settings
 #     kreps = (5, 5, 5)
 #     eparams = PorousMaterials.EwaldParams((5,5,5), 5.6/box.a, sr_cutoff_r, box)
 #     kvecs = PorousMaterials.precompute_kvec_wts(eparams, 27.0)
 #     eikar = OffsetArray(Complex{Float64}, 0:kreps[1])
 #     eikbr = OffsetArray(Complex{Float64}, -kreps[2]:kreps[2])
 #     eikcr = OffsetArray(Complex{Float64}, -kreps[3]:kreps[3])
    energy = PorousMaterials.total_electrostatic_potential_energy(ms, eparams, kvecs, eikar, eikbr, eikcr)
    @assert(isapprox(energy, energies_should_be[c], rtol=0.01))
    println(energy)
    println(energies_should_be[c])
end
