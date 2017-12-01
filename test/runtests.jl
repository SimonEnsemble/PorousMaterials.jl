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
ljforcefield = read_forcefield_file("test_forcefield.csv", cutoffradius=5.0, mixing_rules="Lorentz-Berthelot")
rep_factors = replication_factors(frame.box, 0.25)
@testset "Forcefield Tests" begin
	@test ljforcefield.pure_sigmas["He"] == 1.0
	@test ljforcefield.pure_epsilons["Zn"] == 12.0
	@test ljforcefield.sigmas_squared["Zn"]["He"] == ( (1.0 + 3.0) / 2 ) ^ 2
	@test ljforcefield.epsilons["He"]["Zn"] == sqrt(12.0 * 3.0)
	@test ljforcefield.cutoffradius_squared == 5.0 ^ 2
	@test rep_factors == (1, 1, 1)
end;

@printf("------------------------------\nTesting Energetics.jl\n\n")
molecule1 = Molecule(1,["He"], [0.5, 0.5, 0.5][:,:], [0.0])
molecule2 = Molecule(1,["He"], [0.5 + rep_factors[1], 0.5 + rep_factors[2], 0.5 + rep_factors[3]][:,:], [0.0])
@testset "Energetics Tests" begin
	@test vdw_energy(frame, molecule1, ljforcefield, rep_factors) ≈ vdw_energy(frame, molecule2, ljforcefield, rep_factors)
	@test vdw_energy(frame, molecule1, ljforcefield, rep_factors) ≈ 4 * ljforcefield.epsilons["He"]["Zn"] * ( (ljforcefield.sigmas_squared["Zn"]["He"] / 0.75) ^ 6 - (ljforcefield.sigmas_squared["Zn"]["He"] / 0.75) ^ 3 )
end;
@printf("------------------------------\n")

