using PorousMaterials

framework = Framework("LTA.cif")
molecule = Molecule("Xe")
forcefield = LJForceField("UFF.csv")
grid = energy_grid(framework, molecule, forcefield, n_pts=(10, 10, 10))

segmented_grid = PorousMaterials._segment_grid(grid, energy_tol=0.0, verbose=true)
PorousMaterials._merge_segments_connected_across_periodic_boundary!(segmented_grid)
