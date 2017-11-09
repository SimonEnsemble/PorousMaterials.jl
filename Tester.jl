using Crystal
using Forcefield
using Mols

@printf("=================================\n")
# Path to cssr files
strucPath = pwd()*"/../"


# Make a framework for the MOF
frame = readcssr(strucPath*"SBMOF-1.cssr")
@printf("\n Framework successfully made! \n=================================\n")


# Finding replication factors for a supercell that can contain a sphere with a 14 Angstrom radius
reps = rep_factors(frame, 12.5)

# Hard coding a He atom in
mol = Molecule(1,["Xe"], [0 0 0], [0.0])
@printf("\n Molecule successfully made! \n=================================\n")

# Finding the center of mass in the supercell
c0 = centerOfMass(frame, homedir()*"/Dropbox/Code/PorousMaterials.jl/MolProps.csv", reps)
# Read parameters for universal force fields, used for the lennard jones calcs
ljforcefield = ljffConstruct(homedir()*"/Dropbox/Code/PorousMaterials.jl/UFF.csv")
@printf("\n LJForceField successfully made! \n=================================\n")

# Try to calculate energy at one position
energy = vdw_energy(frame, mol, ljforcefield, [0,0,0], reps)
@printf("\n Energy calculated! Energy = %f\n=================================\n",energy)
