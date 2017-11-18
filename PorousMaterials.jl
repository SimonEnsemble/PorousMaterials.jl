module PorousMaterials

include("Crystal.jl")
include("Forcefield.jl")
include("Molecules.jl")
include("Energetics.jl")

global PATH_TO_STRUCTURE_FILES = homedir() * "/Dropbox/Code/PorousMaterials.jl/cssrFiles"
# TODO somehow let user adjust this variable and look for files outside of the directory? although this is fine for now, it could also be confusing to try to find out where all the files should be placed.
 
export Framework, read_crystal_structure_file, replicate_to_xyz, # Crystal.jl
       LennardJonesForceField, read_forcefield_file, replication_factors, # Forcefield.jl
       Molecule, constructmolecule, readcharge, # Molecules.jl
       lennard_jones, vdw_energy # Energetics.jl

end
