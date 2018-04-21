using PorousMaterials
framework = read_crystal_structure_file("SBMOF-1.cssr")
forcefield = read_forcefield_file("test_forcefield.csv")
molecule = read_molecule_file("He")
grid = energy_grid(framework, molecule, forcefield, n_pts=(13,14,15))
write_cube(grid, "test")
grid_reconstructed = read_cube("test.cube")
@assert(isapprox(grid, grid_reconstructed))
