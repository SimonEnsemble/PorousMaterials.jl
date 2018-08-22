using PorousMaterials
framework = Framework("SBMOF-1.cssr")
forcefield = LJForceField("UFF.csv")
molecule = Molecule("He")
grid = energy_grid(framework, molecule, forcefield, n_pts=(13,14,15))
write_cube(grid, "test")
grid_reconstructed = read_cube("test.cube")
@assert(isapprox(grid, grid_reconstructed))
