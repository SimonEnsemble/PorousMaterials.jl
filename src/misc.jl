function add_extension(filename::String, extension::String)
    if ! occursin(extension, filename)
        filename *= extension
    end
    return filename
end

"""
    atoms = read_xyz(filename)

read a list of atomic species and their corresponding coordinates from an .xyz file.

# Arguments
- `filename::AbstractString`: The filename of the .xyz file

# Returns
- `atoms::Atoms{Cart}`: the set of atoms read from the .xyz file.
"""
function read_xyz(filename::AbstractString)
    f = open(filename)
    lines = readlines(f)
    if length(lines) == 0
        return Symbol[], Float64[]
    end
    n = parse(Int, lines[1]) # get number of atoms
    species = Symbol[]
    x = zeros(Float64, 3, n)
    for i = 1:n
        push!(species, Symbol(split(lines[i + 2])[1]))
        for j = 1:3
            x[j, i] = parse(Float64, split(lines[i + 2])[1 + j])
        end
    end
    close(f)
    return Atoms(species, Cart(x))
end

"""
    write_xyz(atoms, filename; comment="")
    write_xyz(crystal; comment="", center_at_origin=false)
    write_xyz(molecules, box, filename; comment="") # fractional
    write_xyz(molecules, box, filename; comment="") # Cartesian

write atoms to an .xyz file.

# Arguments
- `atoms::Atoms`: the set of atoms.
- `filename::AbstractString`: the filename (absolute path) of the .xyz file. (".xyz" appended automatically
if the extension is not provided.)
- `comment::AbstractString`: comment if you'd like to write to the file.
- `center_at_origin::Bool`: (for crystal only) if `true`, translate all coords such that the origin is the center of the unit cell.
"""
function write_xyz(atoms::Atoms{Cart}, filename::AbstractString; comment::AbstractString="")
    filename = add_extension(filename, ".xyz")

    xyzfile = open(filename, "w")
    @printf(xyzfile, "%d\n%s\n", atoms.n, comment)
    for i = 1:atoms.n
		@printf(xyzfile, "%s\t%.4f\t%.4f\t%.4f\n", atoms.species[i], 
            atoms.coords.x[1, i], atoms.coords.x[2, i], atoms.coords.x[3, i])
    end
    close(xyzfile)
    return nothing
end

"""
    atom_colors = read_cpk_colors()

Read in CPK color scheme for atoms. Return `atom_colors::Dict{Symbol, Tuple{Int, Int, Int}}` such that
`atom_colors[":C"]` gives RGB code for carbon as a tuple, `(144, 144, 144)`.
https://en.wikipedia.org/wiki/CPK_coloring

# Returns
- `atom_colors::Dict{Symbol, Tuple{Int, Int, Int}}`: A dictionary linking an element symbol to its' corresponding CPK color in RGB
"""
function read_cpk_colors()
    atom_colors = Dict{Symbol, Tuple{Int, Int, Int}}()
    df_colors = CSV.read(joinpath(PATH_TO_DATA, "cpk_atom_colors.csv"))
    for row in eachrow(df_colors)
        atom_colors[Symbol(row[:atom])] = (row[:R], row[:G], row[:B])
    end
    return atom_colors
end

"""
    atomic_masses = read_atomic_masses()

Read the `data/atomicmasses.csv` file to construct a dictionary of atoms and their atomic
masses in amu.

# Returns
- `atomic_masses::Dict{Symbol, Float64}`: A dictionary containing the atomic masses of each atom stored in `data/atomicmasses.csv`
"""
function read_atomic_masses()
    if ! isfile(joinpath(PATH_TO_DATA, "atomicmasses.csv"))
        error("Cannot find atomicmasses.csv file in your data folder\n")
    end

    df_am = CSV.read(joinpath(PATH_TO_DATA, "atomicmasses.csv"))

    atomic_masses = Dict{Symbol, Float64}()

    for row in eachrow(df_am)
		atomic_masses[Symbol(row[:atom])] = row[:mass]
    end

    return atomic_masses
end
