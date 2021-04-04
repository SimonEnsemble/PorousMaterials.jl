using PorousMaterials
using OffsetArrays
using LinearAlgebra
using Test
using JLD2
using Statistics
using Random

@testset "Energy Minimizaiton Tests" begin
    xtal = Crystal("SBMOF-1.cif")
    molecule  = Molecule("Xe")
    ljff = LJForceField("UFF")
    n_pts = (15, 15, 15)
    
    # grid search (gs)
    minimized_molecule, min_E = find_energy_minimum_gridsearch(xtal, molecule, ljff, n_pts=n_pts)
    @test isapprox(min_E, -37.0, atol=1.0)

    # fine tune it.
    minimized_molecule, min_E = find_energy_minimum(xtal, minimized_molecule, ljff)
    @test isapprox(min_E, -37.6938, atol=0.001)
    @test isapprox(minimized_molecule.com.xf, [0.00583315, 0.156215, 0.505596], atol=0.001)
    
    xtal = Crystal("cage_in_space.cif")
    min_mol, E_min = find_energy_minimum_gridsearch(xtal, molecule, ljff, n_pts=n_pts)
    @test isapprox(min_mol.com.xf, [0.5, 0.5, 0.5], atol=0.1)
    min_mol, E_min = find_energy_minimum(xtal, min_mol, ljff)
    @test isapprox(E_min, -38.55, atol=0.1)
    @test isapprox(min_mol.com.xf, [0.5, 0.5, 0.5], atol=0.1)
end
