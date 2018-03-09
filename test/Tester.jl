using PorousMaterials

@printf("=================================\n")

filename = "SBMOF-1_cory.cif"
# Make a framework for the MOF

frame = read_crystal_structure_file(filename);
#frame = constructframework("small.cif")
@printf("\n Framework successfully made! \n=================================\n")


# Finding replication factors for a supercell that can contain a sphere with a 14 Angstrom radius
reps = replication_factors(frame.box, 12.5);
@printf("%d,%d,%d\n",reps[1],reps[2],reps[3])

# Hard coding a He atom in
temp = zeros(3,1);
mol = Molecule(1,["Xe"], temp, [0.0]);
@printf("\n Molecule successfully made! \n=================================\n")

# Finding the center of mass in the supercell
# Read parameters for universal force fields, used for the lennard jones calcs
ljforcefield = read_forcefield_file("test_forcefield.csv", cutoffradius=12.5);


@printf("\n LJForceField successfully made! \n=================================\n")


energy1 = vdw_energy(frame, mol, ljforcefield, reps)
@printf("Xe at [%f, %f, %f]:\t Energy = %f\n",mol.x[1],mol.x[2],mol.x[3],energy1)

mol.x[1] = 0.494265; mol.x[2] = 2.22668; mol.x[3] = 0.450354;

energy2 = vdw_energy(frame, mol, ljforcefield, reps)
@printf("Xe at [%f, %f, %f]:\t Energy = %f\n",mol.x[1],mol.x[2],mol.x[3],energy2)

# Try to calculate energy at one position
#energy = vdw_energy(frame, mol, ljforcefield, reps);
#@printf("\n Energy calculated! Energy = %f\n=================================\n",energy)


# Try to explore frame
#occupancy = takesnapshot(frame, mol, ljforcefield, [10.,10.,10.] , reps, mesh=[21,21,21], startpoint=[15.,15.,15]);
#@printf("Occupancy Calculated! See variable `occupancy`\n=================================\n")

#replicate_to_xyz(frame,repfactors=reps)

#grid = Grid(frame.box, (30, 30, 30), (occupancy+0), "Int64")

#write_to_cube(grid, "ZASJAG_test1.cube")
