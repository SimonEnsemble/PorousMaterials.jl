"""
    atoms, x = read_xyz(filename)

Return the list of `atoms` (Array{Symbol, 1}) and their Cartesian coordinates 
`x::Array{Float64, 2}` as stored in the .xyz file. `x[:, k]` will return Cartesian
coords of the kth atom.

# Arguments
- `filename::AbstractString`: The filename of the .xyz file

# Returns
- `atoms::Array{Symbol, 1}`: A symbolic array of the atoms in the .xyz file
- `x::Array{Float64, 2}`: The cartesian coordinates of the atoms. `x[:, k]` will return cartesian coordinates of the k-th atom
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
    df_props = CSV.read(PATH_TO_DATA * "atom_properties.csv", nullable=true)
    for row in eachrow(df_props)
        if ! ismissing(row[Symbol("atomicradius[Angstrom]")])
            atomic_radii[Symbol(row[:atom])] = row[Symbol("atomicradius[Angstrom]")]
        end
    end
    return atomic_radii
end
