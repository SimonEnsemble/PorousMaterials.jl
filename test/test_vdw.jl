
using BenchmarkTools, Compat
using PorousMaterials
using Base.Test
    
ljforcefield = LJForceField("Dreiding.csv", cutoffradius=12.5, mixing_rules="Lorentz-Berthelot") # Dreiding

# Xe in SBMOF-1 tests, comparing to RASPA
sbmof1 = Framework("SBMOF-1.cif")
rep_factors_sbmof1 = replication_factors(sbmof1.box, ljforcefield)
sbmof1 = replicate(sbmof1, rep_factors_sbmof1)
xenon = Molecule("Xe")
set_fractional_coords!(xenon, sbmof1.box)

@test ! charged(xenon)
xenon.atoms[1].xf[:] = sbmof1.box.c_to_f * zeros(3)
energy = vdw_energy(sbmof1, xenon, ljforcefield)
@test isapprox(energy, -5041.58, atol = 0.005)
xenon.atoms[1].xf[1] = 0.494265; xenon.atoms[1].xf[2] = 2.22668; xenon.atoms[1].xf[3] = 0.450354;
xenon.atoms[1].xf[:] = sbmof1.box.c_to_f * xenon.atoms[1].xf[:]
energy = vdw_energy(sbmof1, xenon, ljforcefield)
@test isapprox(energy, 12945.838, atol = 0.005)

vdw_energy(sbmof1, xenon, ljforcefield)
@btime vdw_energy(sbmof1, xenon, ljforcefield)
