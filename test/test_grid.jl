using PorousMaterials
using Test

@testset "Grid Tests" begin
    framework = Framework("SBMOF-1.cssr")
    forcefield = LJForceField("UFF.csv")
    molecule = Molecule("He")
    grid = energy_grid(framework, molecule, forcefield, n_pts=(13,14,15), verbose=false)
    write_cube(grid, "test")
    grid_reconstructed = read_cube("test.cube")
    @test isapprox(grid, grid_reconstructed)
end
