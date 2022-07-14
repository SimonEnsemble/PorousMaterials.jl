using BenchmarkTools
using PorousMaterials
using Test
using Random
using Printf

ljforcefield = LJForceField("Dreiding", r_cutoff=12.5, mixing_rules="Lorentz-Berthelot") # Dreiding

# Xe in SBMOF-1 tests, comparing to RASPA
sbmof1 = Crystal("SBMOF-1.cif")
rep_factors_sbmof1 = replication_factors(sbmof1.box, ljforcefield)
sbmof1 = replicate(sbmof1, rep_factors_sbmof1)
xenon = Molecule("Xe")
xenon = Frac(xenon, sbmof1.box)

@test ! has_charges(xenon)
xenon.atoms.coords.xf[:, 1] = sbmof1.box.c_to_f * zeros(3)
energy = vdw_energy(sbmof1, xenon, ljforcefield)
@test isapprox(energy, -5041.58, atol = 0.005)
xenon.atoms.coords.xf[1, 1] = 0.494265; xenon.atoms.coords.xf[2, 1] = 2.22668; xenon.atoms.coords.xf[3, 1] = 0.450354;
xenon.atoms.coords.xf[:, 1] = sbmof1.box.c_to_f * xenon.atoms.coords.xf[:, 1]
energy = vdw_energy(sbmof1, xenon, ljforcefield)
@test isapprox(energy, 12945.838, atol = 0.005)

# guest-host vdw_energy timing
vdw_energy(sbmof1, xenon, ljforcefield)
println("Guest-Host Van der Waals energy computation:")
@btime vdw_energy(sbmof1, xenon, ljforcefield)

# guest-guest vdw_energy timing
Random.seed!(1234) # so that every trial will be the same
box = Box(100.0, 100.0, 100.0, π/2, π/2, π/2)
co2 = Molecule("CO2")
ms = Array{Molecule{Frac}, 1}()
num_molecules = 100
while length(ms) < num_molecules # give a decent number of molecules to test
    random_insertion!(ms, box, co2)
    # Ensure that no molecules have infinite lj energy, so that it calculates
    #   the interaction between all molecules
    if vdw_energy(length(ms), ms, ljforcefield, box) >= Inf
        pop!(ms)
    end
end
vdw_energy(1, ms, ljforcefield, box)
println("Guest-Guest Van der Waals energy computation for a single molecule:")
println("\t# of molecules: ", length(ms))
@btime vdw_energy(1, ms, ljforcefield, box)
@printf("Guest-Guest Van der Waals energy computation for all %i molecules:\n", num_molecules)
@btime begin
    for i = eachindex(ms)
        vdw_energy(i, ms, ljforcefield, box)
    end
end
