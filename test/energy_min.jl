using PorousMaterials
using OffsetArrays
using LinearAlgebra
using Test
using JLD2
using Statistics
using Random

FIQCEN_clean = Crystal("FIQCEN_clean.cif")
Xe = Molecule("Xe")

@testset "Energy Minimizaiton Tests" begin
    xtal = deepcopy(FIQCEN_clean)
    strip_numbers_from_atom_labels!(xtal)
    molecule = deepcopy(Xe)
    ljff = LJForceField("UFF")
    n_pts = (15, 15, 15)

    # grid search (gs)
    min_mol_gs, min_E_gs = find_energy_minimum_gridsearch(xtal, molecule, ljff; n_pts=n_pts)
    @test isapprox(min_E_gs, -25.1, atol=1.0)
    @test isapprox(min_mol_gs.com.xf, [0.785, 0.785, 0.785], atol=0.05)

    # fine tune it.
    min_mol, min_E = find_energy_minimum(xtal, min_mol_gs, ljff)
    @test isapprox(min_E, -27.0269, atol=0.01)
    @test min_E < min_E_gs
    @test isapprox(min_mol.com.xf, [0.75, 0.75, 0.75], atol=0.001)
    # test replications handled correctly
    min_mol = Cart(min_mol, xtal.box)
    rep_factors = replication_factors(xtal, sqrt(ljff.r²_cutoff))
    xtal = replicate(xtal, rep_factors)
    @test isapprox(vdw_energy(xtal, min_mol, ljff) * 8.314 / 1000, min_E, atol=0.1)

    # make sure it works with cart coords too.
    xtal = deepcopy(FIQCEN_clean)
    strip_numbers_from_atom_labels!(xtal)
    molecule = deepcopy(Xe)
    translate_to!(molecule, Cart(Frac([0.785, 0.785, 0.785]), xtal.box))
    min_mol, min_E = find_energy_minimum(xtal, molecule, ljff)
    @test isapprox(min_E, -27.0269, atol=0.01)
    @test isapprox(min_mol.com.xf, [0.75, 0.75, 0.75], atol=0.001)

    # one more, know the min is in the mid b/c this is a centered cage.
    xtal = Crystal("cage_in_space.cif")
    min_mol, min_E = find_energy_minimum_gridsearch(xtal, molecule, ljff; n_pts=n_pts)
    @test isapprox(min_mol.com.xf, [0.5, 0.5, 0.5], atol=0.1)
    min_mol, min_E = find_energy_minimum(xtal, min_mol, ljff)
    @test isapprox(min_E, -38.55, atol=0.1)
    @test isapprox(min_mol.com.xf, [0.5, 0.5, 0.5], atol=0.1)
    @test isapprox(vdw_energy(xtal, min_mol, ljff) * 8.314 / 1000, min_E, atol=0.1)
end
