module Grid_Test

using PorousMaterials
using OffsetArrays
using LinearAlgebra
using Test
using JLD2
using Statistics
using Random

@testset "Grid Tests" begin
    # test read and write
    grid = Grid(Box(0.7, 0.8, 0.9, 1.5, 1.6, 1.7), (3, 3, 3), rand(Float64, (3, 3, 3)),
        :kJ_mol, [1., 2., 3.])
    write_cube(grid, "test_grid.cube")
    grid2 = read_cube("test_grid.cube")
    @test isapprox(grid, grid2)
    
    framework = Framework("LTA.cif")
    molecule = Molecule("Xe")
    forcefield = LJForceField("UFF.csv")
    grid = energy_grid(framework, molecule, forcefield, n_pts=(10, 10, 10))
end
end
