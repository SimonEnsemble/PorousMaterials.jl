using PorousMaterials

framework = Framework("LTA.cif")
write_xyz(framework)
molecule = Molecule("CH4")
forcefield = LJForceField("UFF.csv")
grid = energy_grid(framework, molecule, forcefield, n_pts=(20, 20, 20))

segmented_grid = PorousMaterials._segment_grid(grid, 298.0, true)
verbose = true
write_cube(segmented_grid, "segmented_grid_LTA_b4.cube")
connections = PorousMaterials._build_list_of_connections(segmented_grid)
graph, vertex_to_direction = PorousMaterials._translate_into_graph(segmented_grid, connections)
segment_classifications = PorousMaterials._classify_segments(segmented_grid, graph, vertex_to_direction, verbose)
PorousMaterials._assign_inaccessible_pockets_minus_one!(segmented_grid, segment_classifications)
write_cube(segmented_grid, "segmented_grid_LTA_after.cube")

accessibility_grid, some_pockets_were_blocked = compute_accessibility_grid(framework, molecule, forcefield, n_pts=(20, 20, 20), energy_tol=0.0, verbose=true, write_b4_after_grids=true)
