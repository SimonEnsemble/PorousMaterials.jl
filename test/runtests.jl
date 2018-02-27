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
const ljforcefield = read_forcefield_file("test_forcefield.csv", cutoffradius=12.5, mixing_rules="Lorentz-Berthelot")
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

#
#GCMC Tests
#

"""

    use_nearest_image = completely_outside_box(molecule, sim_box)

returns true if the given molecule is completely outside of a given box
"""
function completely_outside_box(molecule::Molecule, box::Box)
    xf = box.c_to_f * molecule.x
    for coords = 1:3
        if sum(xf[coords, :] .<= 1.0) == 0
            return false
        else if sum(xf[coords, :] .>= 0.0) == 0
            return false
        end
    end
    return true
end

#
#INSERTION TESTS
#
molecules = Array{Molecule}(0)
repfactors = replication_factors(frame.box, ljforcefield)
sim_box = replicate_box(frame.box, repfactors)

for i = 1:100
    insert_molecule!(molecules, sim_box, "CH4")
    @assert(completely_outside_box(molecules[i], sim_box), "Molecule outside of simulation box")
    @assert(length(molecules) == i, "Molecules array not modified")
    if i > 1
        @assert(sum(molecules[i - 1].x .≈ molecules[i].x) == 0, "Molecules not being inserted at random coordinates")
    end
end

#
#DELETION TESTS
#
for i = 1:100
    delete_molecule!(molecules, rand(1:length(molecules)))
    @assert(length(molecules) == 100 - i, "Molecule not deleted")
end

#
#TRANSLATION TESTS
#
molecules = [Molecule(1, "CH4", sim_box.c_to_f * [0.99, 0.99, 0.01], [0.0])]
for i = 1:10000
    old_coords = molecules[1].x
    translate_molecule!(1, molecules, sim_box)
    @assert(sum(molecules[1].x .≈ old_coords) == 0, "Molecule not moved")
    @assert(completely_outside_box(molecules[1], sim_box),
        "Molecule completely outside of simulation box")
end

#
#VAN DER WAALS GUEST-GUEST ENERGY TESTS
#
#TODO make sim box of 25 25 25

sim_box = construct_box(25, 25, 25, π/2, π/2, π/2)
#TODO add CH4 to UFF file to have sigmas and epsilons ready
molecules = [Molecule(1, "C", [5.0, 12.0, 12.0], [0.0],
    Molecule(1, "C", [11.0, 12.0, 12.0], [0.0])]
@assert(guest_guest_vdw_energy(1, molecules, ljforcefield, sim_box) == -6.11807,
    "Did not calculate energy correctly between two molecules")
delete_molecule!(2, molecules)
push!(molecules, Molecule(1, "C", [24.0, 12.0, 12.0], [0.0]))
@assert(guest_guest_vdw_energy(1, molecules, ljforcefield, sim_box) == -6.11087,
    "Did not use nearest image convention")
push!(molecules, Molecule(1, "C", [11.0, 12.0, 12.0], [0.0]))
@assert(guest_guest_vdw_energy(1, molecules, ljforcefield, sim_box) == -6.11087 * 2,
    "Did not calculate guest-guest energy correctly for more than two molecules")
