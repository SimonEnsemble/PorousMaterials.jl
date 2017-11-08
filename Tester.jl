using Crystal
using Forcefield
using Mols

@printf("=================================\n")
# Path to cssr files
strucPath = homedir()*"/Dropbox/Code/PorousMaterials.jl/cssrFiles/"


# Make a framework for the MOF
frame = readcssr(strucPath*"ABUWOJ_clean_min_charges.cssr")
@printf("\n Framework successfully made! \n=================================\n")


# Finding replication factors for a supercell that can contain a sphere with a 14 Angstrom radius
reps = rep_factors(frame, 14.0)

# Hard coding a He atom in
mol = Molecule(1,["He"], [0 0 0], [0.0])
@printf("\n Molecule successfully made! \n=================================\n")

# Finding the center of mass in the supercell
c0 = centerOfMass(frame, homedir()*"/Dropbox/Code/PorousMaterials.jl/MolProps.csv", reps)
# Read parameters for universal force fields, used for the lennard jones calcs
ljforcefield = ljffConstruct(homedir()*"/Dropbox/Code/PorousMaterials.jl/UFF.csv")
@printf("\n LJForceField successfully made! \n=================================\n")

# Try to calculate energy at one position
energy = vdw_energy(frame, mol, ljforcefield, c0, reps)
@printf("\n Energy calculated! Energy = %f\n=================================\n",energy/6.022e23)
