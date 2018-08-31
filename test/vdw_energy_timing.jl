
using BenchmarkTools, Compat
using PorousMaterials
using Test
using Random

ljforcefield = LJForceField("Dreiding.csv", cutoffradius=12.5, mixing_rules="Lorentz-Berthelot") # Dreiding

# Xe in SBMOF-1 tests, comparing to RASPA
sbmof1 = Framework("SBMOF-1.cif")
rep_factors_sbmof1 = replication_factors(sbmof1.box, ljforcefield)
sbmof1 = replicate(sbmof1, rep_factors_sbmof1)
xenon = Molecule("Xe")
set_fractional_coords!(xenon, sbmof1.box)

@test ! charged(xenon)
xenon.atoms.xf[:, 1] = sbmof1.box.c_to_f * zeros(3)
energy = vdw_energy(sbmof1, xenon, ljforcefield)
@test isapprox(energy, -5041.58, atol = 0.005)
xenon.atoms.xf[1, 1] = 0.494265; xenon.atoms.xf[2, 1] = 2.22668; xenon.atoms.xf[3, 1] = 0.450354;
xenon.atoms.xf[:, 1] = sbmof1.box.c_to_f * xenon.atoms.xf[:, 1]
energy = vdw_energy(sbmof1, xenon, ljforcefield)
@test isapprox(energy, 12945.838, atol = 0.005)

# guest-host vdw_energy timing
vdw_energy(sbmof1, xenon, ljforcefield)
@btime vdw_energy(sbmof1, xenon, ljforcefield)

# guest-guest vdw_energy timing
Random.seed!(1234) # so that every trial will be the same
box = UnitCube()
co2 = Molecule("CO2")
ms = Array{Molecule, 1}()
for i = 1:100 # give a decent number of molecules to test
    insert_molecule!(ms, box, co2)
end
vdw_energy(1, ms, ljforcefield, box)
@btime vdw_energy(1, ms, ljforcefield, box)
@btime begin
    for i = 1:length(ms)
        vdw_energy(i, ms, ljforcefield, box)
    end
end
