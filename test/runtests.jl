#!/usr/bin/env julia

# Details from http://www.stochasticlifestyle.com/finalizing-julia-package-documentation-testing-coverage-publishing/
# Start Test Script
using PorousMaterials
using Base.Test

# Run Tests

@printf("------------------------------\nTesting Crystal.jl\n\n")
frame = read_crystal_structure_file("test_structure.cif")
@testset "Crystal Tests" begin
	@test frame.box.f_to_c * frame.box.c_to_f ≈ eye(3)
	@test frame.box.a == 1
	@test frame.box.β == 90 * (pi / 180)
	@test frame.box.Ω ≈ 1
	@test frame.xf == [0.0, 0.0, 0.0][:,:]
	@test frame.atoms == ["Zn"]
end;

@printf("------------------------------\nTesting Forcefield.jl\n\n") 
ljforcefield = read_forcefield_file("test_forcefield.csv", cutoffradius=12.5, mixing_rules="Lorentz-Berthelot")
rep_factors = replication_factors(frame.box, ljforcefield)
@testset "Forcefield Tests" begin
	@test ljforcefield.pure_sigmas["He"] == 1.0
	@test ljforcefield.pure_epsilons["Zn"] == 12.0
	@test ljforcefield.sigmas_squared["Zn"]["He"] == ( (1.0 + 3.0) / 2 ) ^ 2
	@test ljforcefield.epsilons["He"]["Zn"] == sqrt(12.0 * 3.0)
	@test ljforcefield.cutoffradius_squared == 12.5 ^ 2
	@test rep_factors == (25, 25, 25)
end;

@printf("------------------------------\nTesting Energetics.jl\n\n")
molecule1 = Molecule(1,["He"], [0.5, 0.5, 0.5][:,:], [0.0])
molecule2 = Molecule(1,["He"], [0.5 + rep_factors[1], 0.5 + rep_factors[2], 0.5 + rep_factors[3]][:,:], [0.0])
frame2 = read_crystal_structure_file("SBMOF-1.cif")
rep_factors2 = replication_factors(frame2.box, ljforcefield)
molecule3 = Molecule(1,["Xe"], zeros(3,1), [0.0])
energy1 = vdw_energy(frame2, molecule3, ljforcefield, rep_factors2)
molecule3.x[1] = 0.494265; molecule3.x[2] = 2.22668; molecule3.x[3] = 0.450354;
energy2 = vdw_energy(frame2, molecule3, ljforcefield, rep_factors2)
@testset "Energetics Tests" begin
	@test vdw_energy(frame, molecule1, ljforcefield, rep_factors) ≈ vdw_energy(frame, molecule2, ljforcefield, rep_factors)
	@test vdw_energy(frame, molecule1, ljforcefield, (1,1,1)) ≈ 4 * ljforcefield.epsilons["He"]["Zn"] * ( (ljforcefield.sigmas_squared["Zn"]["He"] / 0.75) ^ 6 - (ljforcefield.sigmas_squared["Zn"]["He"] / 0.75) ^ 3 )
	@test isapprox(energy1, -5041.58, atol = 0.005)
	@test isapprox(energy2, 12945.838, atol = 0.005)
end;

@printf("------------------------------\n")

