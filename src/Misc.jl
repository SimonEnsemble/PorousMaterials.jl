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
    write_xyz(atoms, x, filename; comment="")
    write_xyz(molecules, box, filename; comment="")
    write_xyz(framework, filename; comment="", center=false)

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
- `center::Bool`: shift atoms so that origin is the center of the `framework.box`
"""
function write_xyz(atoms::Array{Symbol, 1}, x::Array{Float64, 2},
    filename::AbstractString; comment::AbstractString="")

    if ! occursin(".xyz", filename)
        filename *= ".xyz"
    end

    if length(atoms) != size(x)[2]
        error("Number of atoms does not match number of coordinates provided!\n")
    end
    if size(x)[1] != 3
        error("x should be 3 by number of atoms.\n")
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
    df_colors = CSV.read(joinpath(PATH_TO_DATA, "cpk_atom_colors.csv"))
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
    df_props = CSV.read(joinpath(PATH_TO_DATA, "atom_properties.csv"))
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

"""
    henry, rmse = extract_henry_coefficient(df, pressure_col_name, loading_col_name, nb_pts_to_include)

Extract the Henry coefficient from an isotherm contained in `df`. Will use `nb_pts_to_include` points to find the slope (which is equivalent
to the Henry coefficient). Does not convert any units.

# Arguments
- `df::DataFrame`: A DataFrame object that contains the isotherm we wish to extract the Henry coefficient from
- `pressure_col_name::Symbol`: The name of the pressure column. Takes in symbols due to DataFrame structure
- `loading_col_name::Symbol`: The name of the loading/adsorption column. Takes in symbols due to DataFrame structure
- `nb_pts_to_include::Int`: The number of points used to fit the best Henry coefficient

# Returns
- `henry::Float64`: The henry coefficient from the isotherm in `df`
- `rmse::Float64`: The Root Mean Squared Error of the fit
"""
function extract_henry_coefficient(df::DataFrame, pressure_col_name::Symbol, loading_col_name::Symbol, nb_pts_to_include::Int)
    sort!(df, [pressure_col_name])
    P = zeros(nb_pts_to_include, 1)
    n = zeros(nb_pts_to_include)
    for (i, row) in enumerate(eachrow(df))
        if i > nb_pts_to_include
            break
        end
        P[i, 1] = row[pressure_col_name]
        n[i] = row[loading_col_name]
    end
    henry = llsq(P, n; bias=false)[1]
    prediction = henry * P
    rmse = sqrt(sum(abs2.(prediction - n))/length(prediction))
    return henry, rmse
end

"""
    M, K = fit_langmuir(df, pressure_col_name, loading_col_name, M0=1.0, K0=1.0)

Takes in a DataFrame `df` containing an isotherm. Will try to fit a Langmuir model to the data.
The Langmuir model is in the following form:
N = (MKP)/(1+KP)
where N is the total adsorption, M is the maximum monolayer coverage, K is the Langmuir constant and P is the partial pressure of the gas.
A minimization problem is solved using `M0` and `K0` as initial guesses for the maximum monolayer coverage and Langmuir constant, respectively.

# Arguments
- `df::DataFrame`: The DataFrame containing the pressure and adsorption data for the isotherm
- `pressure_col_name::Symbol`: The header of the pressure column. Can be found with `names(df)`
- `loading_col_name::Symbol`: The header of the loading/adsorption column. Can be found with `names(df)`
- `M0::Float64`: An initial guess for the maximum monolayer coverage
- `K0::Float64`: An initial guess for the Langmuir constant

# Returns
- `M::Float64`: The maximum monolayer coverage according to Langmuir's theory, N = (MKP)/(1+KP)
- `K::Float64`: The Langmuir constant
"""
function fit_langmuir(df::DataFrame, pressure_col_name::Symbol, loading_col_name::Symbol, method::Symbol=:nonlinear)
    sort!(df, [pressure_col_name])
    if method == :nonlinear
        n = df[loading_col_name]
        p = df[pressure_col_name]
        #M0 = n[end] * 1.1
        #K0 = n[2] / (p[2] * (M0-n[2]))
        M0, K0 = fit_langmuir(df, pressure_col_name, loading_col_name, :linear)
        function_to_minimize(θ) = return sum([(n[i] - θ[1] * θ[2] * p[i]/(1 + θ[2] * p[i]))^2 for i = 1:length(n)])
        θ0 = [M0, K0]
        res = optimize(function_to_minimize, θ0, LBFGS())
        M, K = res.minimizer
    elseif method == :linear
        inv_N = (df[loading_col_name]).^(-1)
        inv_P = (df[pressure_col_name]).^(-1)
        if inv_P[1] == Inf || inv_N[1] == Inf
            slope, intercept = llsq(inv_P[2:end,:], inv_N[2:end])
        else
            slope, intercept = llsq(inv_P, inv_N)
        end
        M = 1/intercept
        K = 1/(M * slope)
    else
        error("Method not available. Try `:linear` or `:nonlinear`")
    end
    return M, K
end
