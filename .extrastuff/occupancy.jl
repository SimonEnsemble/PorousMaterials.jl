using PorousMaterials
#using JLD
#using FileIO
#using Pandas
using NPZ

@printf("=================================\n")

files = readdir("data/crystals")

ljforcefield = read_forcefield_file("UFF.csv", cutoffradius = 12.5)

for filename in files[3:5]
    frame = read_crystal_structure_file(filename)
    strip_numbers_from_atom_labels!(frame)
    mol = Molecule(1, ["He"], zeros(3,1), [0.0])

    reps = replication_factors(frame.box, ljforcefield)

    occupancy = takesnapshot(frame, mol, ljforcefield, [10., 10., 10.], reps, startpoint = [0., 0., 0.])
    str = split(filename, ".")[1]
    npzwrite("data/db/" * str * ".npy", convert.(UInt8, occupancy))
#    save(File(format"JLD","data/db/" * str * ".jld"),"occ",occupancy)
#    df = DataFrame((occupancy+0)[:])
#    to_json(df, "data/db/" * str * ".json")
    @printf("See data/db/%s.npy\n", str)
end
