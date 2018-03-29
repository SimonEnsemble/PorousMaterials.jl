using PorousMaterials

@printf("=================================\n")

filename = "SBMOF-1_cory.cif"
# Make a framework for the MOF

frame = read_crystal_structure_file(filename);
#frame = constructframework("small.cif")
@printf("\n Framework successfully made! \n=================================\n")



# Hard coding a He atom in
temp = zeros(3,1);
mol = Molecule(1,["Xe"], temp, [0.0]);
@printf("\n Molecule successfully made! \n=================================\n")

# Finding the center of mass in the supercell
# Read parameters for universal force fields, used for the lennard jones calcs
ljforcefield = read_forcefield_file("test_forcefield.csv", cutoffradius=12.5);


@printf("\n LJForceField successfully made! \n=================================\n")

# Finding replication factors for a supercell that can contain a sphere with a 14 Angstrom radius
reps = replication_factors(frame.box, ljforcefield);
@printf("%d,%d,%d\n",reps[1],reps[2],reps[3])


# Try to calculate energy at one position
#energy = vdw_energy(frame, mol, ljforcefield, reps);
#@printf("\n Energy calculated! Energy = %f\n=================================\n",energy)


(EnergyM, CoordM, cnt) = exploreframe(frame, mol, ljforcefield, (1, 1, 1))
# Try to explore frame
#occupancy = takesnapshot(frame, mol, ljforcefield, [10.,10.,10.] , reps, mesh=[21,21,21], startpoint=[15.,15.,15]);
#@printf("Occupancy Calculated! See variable `occupancy`\n=================================\n")

replicate_to_xyz(frame, "imgposter.xyz", repfactors=(1,1,1), negative_replications = true)

write_unitcell_boundary_vtk(frame, "imgposter.vtk")

grid = Grid(frame.box, cnt, EnergyM, "Float64", (1, 1, 1))

write_to_cube(grid, "imgposter.cube")
