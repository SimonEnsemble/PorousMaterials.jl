using PorousMaterials

framework = Framework("LTA.cif")
molecule = Molecule("CH4")
forcefield = LJForceField("UFF.csv")
grid = energy_grid(framework, molecule, forcefield, n_pts=(20, 20, 20))

segmented_grid = PorousMaterials._segment_grid(grid, energy_tol=0.0, verbose=true)
write_cube(segmented_grid, "segmented_grid_LTA_b4.cube")
graph = PorousMaterials._build_connectivity_graph(segmented_grid)
segment_classifications = PorousMaterials._classify_segments(segmented_grid, graph)
PorousMaterials._assign_inaccessible_pockets_minus_one!(segmented_grid, segment_classifications)
write_cube(segmented_grid, "segmented_grid_LTA_after.cube")

accessibility_grid, some_pockets_were_blocked = compute_accessibility_grid(framework, molecule, forcefield, n_pts=(20, 20, 20), energy_tol=0.0, verbose=true, write_b4_after_grids=true)
