"""
    atoms, x = read_xyz(filename)

Return the list of `atoms` (Array{Symbol, 1}) and their Cartesian coordinates
`x::Array{Float64, 2}` as stored in the .xyz file. `x[:, k]` will return Cartesian
coords of the kth atom.

# Arguments
- `filename::AbstractString`: The filename of the .xyz file

# Returns
- `atoms::Array{Symbol, 1}`: An array of atoms stored as symbols e.g. [:H, :H, :O] read
from the .xyz file.
- `x::Array{Float64, 2}`: The Cartesian coordinates of the atoms. `x[:, k]` will return cartesian coordinates of the k-th atom
"""
function read_xyz(filename::AbstractString)
    f = open(filename)
    lines = readlines(f)
    if length(lines) == 0
        return Symbol[], Float64[]
    end
    n = parse(Int, lines[1]) # get number of atoms
    atoms = Symbol[]
    x = zeros(Float64, 3, n)
    for i = 1:n
        push!(atoms, Symbol(split(lines[i + 2])[1]))
        for j = 1:3
            x[j, i] = parse(Float64, split(lines[i + 2])[1 + j])
        end
    end
    close(f)
    return atoms, x
end

"""
    write_to_xyz(atoms, x, filename; comment="")
    write_to_xyz(molecules, box, filename; comment="")
    write_to_xyz(framework, filename; comment="")

Write a molecule, framework, or array of atoms & positions to an .xyz file.

# Arguments
- `atoms::Array{Symbol, 1}`: An array of atoms stored as symbols e.g. [:H, :H, :O]
- `x::Array{Float64, 2}`: The Cartesian coordinates of the atoms.
 `x[:, k]` contains Cartesian coordinates of the k-th atom
- `molecules::Array{Molecule, 1}`: an array of molecules whose atoms to write to .xyz
- `framework::Framework`: a crystal structure whose atoms to write to .xyz
- `filename::AbstractString`: The filename of the .xyz file. (".xyz" appended automatically
if the extension is not provided.) (absolute path)
- `comment::AbstractString`: comment if you'd like to write to the file.
"""
function write_to_xyz(atoms::Array{Symbol, 1}, x::Array{Float64, 2}, 
    filename::AbstractString; comment::AbstractString="")

    if ! contains(filename, ".xyz")
        filename *= ".xyz"
    end

    xyzfile = open(filename, "w")
    @printf(xyzfile, "%d\n%s\n", length(atoms), comment)
    for i = 1:length(atoms)
		@printf(xyzfile, "%s\t%.4f\t%.4f\t%.4f\n", atoms[i], x[1, i], x[2, i], x[3, i])
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
    df_colors = CSV.read(PATH_TO_DATA * "cpk_atom_colors.csv")
    for row in eachrow(df_colors)
        atom_colors[Symbol(row[:atom])] = (row[:R], row[:G], row[:B])
    end
    return atom_colors
end

"""
    atomic_radii = read_atomic_radii()

Return `atomic_radii::Dict{Symbol, Float64}`, where `atom_masses[":C"]` gives
the atomic radii of carbon (10.87 Angstrom).

# Returns
- `atomic_radii::Dict{Symbol, Float64}`: A dictionary linking an element symbol to its' corresponding atomic radius
"""
function read_atomic_radii()
    atomic_radii = Dict{Symbol, Float64}()
    df_props = CSV.read(PATH_TO_DATA * "atom_properties.csv")
    for row in eachrow(df_props)
        if ! ismissing(row[Symbol("atomicradius[Angstrom]")])
            atomic_radii[Symbol(row[:atom])] = row[Symbol("atomicradius[Angstrom]")]
        end
    end
    return atomic_radii
end

"""
    atomic_masses = read_atomic_masses()

Read the `data/atomicmasses.csv` file to construct a dictionary of atoms and their atomic
masses in amu.

# Returns
- `atomic_masses::Dict{Symbol, Float64}`: A dictionary containing the atomic masses of each atom stored in `data/atomicmasses.csv`
"""
function read_atomic_masses()
    if ! isfile(PATH_TO_DATA * "atomicmasses.csv")
        error("Cannot find atomicmasses.csv file in your data folder\n")
    end

    df_am = CSV.read(PATH_TO_DATA * "atomicmasses.csv")

    atomic_masses = Dict{Symbol, Float64}()

    for row in eachrow(df_am)
		atomic_masses[Symbol(row[:atom])] = row[:mass]
    end

    return atomic_masses
end
