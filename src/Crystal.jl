using Base.Test
using PyCall # see pyimport, requires ASE

"""
Data structure to describe a unit cell box (Bravais lattice) and convert between
fractional and Cartesian coordinates.

# Attributes
- `a,b,c::Float64`: unit cell dimensions (units: Angstroms)
- `α,β,γ::Float64`: unit cell angles (units: radians)
- `Ω::Float64`: volume of the unit cell (units: cubic Angtroms)
- `f_to_c::Array{Float64,2}`: the 3x3 transformation matrix used to map fractional
coordinates to cartesian coordinates. The columns of this matrix define the unit cell
axes.
- `c_to_f::Array{Float64,2}`: the 3x3 transformation matrix used to map Cartesian
coordinates to fractional coordinates
- `reciprocal_lattice::Array{Float64, 2}`: the columns are the reciprocal lattice vectors
"""
struct Box
    a::Float64
    b::Float64
    c::Float64

    α::Float64
    β::Float64
    γ::Float64

    Ω::Float64

    f_to_c::Array{Float64, 2}
    c_to_f::Array{Float64, 2}

    reciprocal_lattice::Array{Float64, 2}
end

"""
    r = reciprocal_lattice(a₁, a₂, a₃)

Given the unit cell vectors defining the Bravais lattice, a₁, a₂, a₃, compute the reciprocal lattice vectors.

# Arguments
- `a₁::Array{Float64, 1}`: The first unit cell vector of a Bravais lattice
- `a₂::Array{Float64, 1}`: The second unit cell vector of a Bravais lattice
- `a₃::Array{Float64, 1}`: The third unit cell vector of a Bravais lattice

# Returns
- `r::Array{Float64, 2}`: Reciprocal lattice vectors in a matrix format, where the columns are the reciprocal lattice vectors.
"""
function reciprocal_lattice(a₁::Array{Float64, 1}, a₂::Array{Float64, 1}, a₃::Array{Float64, 1})
    r = zeros(Float64, 3, 3)
    r[:, 1] = 2 * π * cross(a₂, a₃) / dot(a₁, cross(a₂, a₃))
    r[:, 2] = 2 * π * cross(a₃, a₁) / dot(a₂, cross(a₃, a₁))
    r[:, 3] = 2 * π * cross(a₁, a₂) / dot(a₃, cross(a₁, a₂))
    return r
end

"""
    box = construct_box(a, b, c, α, β, γ)

Constructs a `Box` with unit cell dimensions a, b, and c and angles α, β, and γ.
Automatically calculates Ω, f_to_c, and c_to_f for `Box` data structure and returns a `Box`.

# Arguments
- `a,b,c::Float64`: unit cell dimensions (units: Angstroms)
- `α,β,γ::Float64`: unit cell angles (units: radians)

# Returns
- `box::Box`: Fully formed Box object
"""
function construct_box(a::Float64, b::Float64, c::Float64,
                       α::Float64, β::Float64, γ::Float64)
    # unit cell volume (A³)
    Ω = a * b * c * sqrt(1 - cos(α) ^ 2 - cos(β) ^ 2 - cos(γ) ^ 2 + 2 * cos(α) * cos(β) * cos(γ))
    # matrices to map fractional coords <--> Cartesian coords
    f_to_c = [[a, 0, 0] [b * cos(γ), b * sin(γ), 0] [c * cos(β), c * (cos(α) - cos(β) * cos(γ)) / sin(γ), Ω / (a * b * sin(γ))]]
    c_to_f = [[1/a, 0, 0] [-cos(γ) / (a * sin(γ)), 1 / (b * sin(γ)), 0] [b * c * (cos(α) * cos(γ) - cos(β)) / (Ω * sin(γ)), a * c * (cos(β) * cos(γ) - cos(α)) / (Ω * sin(γ)), a * b * sin(γ) / Ω]]
    # the columns of f_to_c are the unit cell axes
    r = reciprocal_lattice(f_to_c[:, 1], f_to_c[:, 2], f_to_c[:, 3])

    @test f_to_c * c_to_f ≈ eye(3)
    @test isapprox(transpose(r), 2.0 * π * inv(f_to_c))

    return Box(a, b, c, α, β, γ, Ω, f_to_c, c_to_f, r)
end

"""
    box = construct_box(f_to_c)

Reconstructs a `Box` from the fractional to Cartesian coordinate transformation matrix.
The columns of this matrix are the unit cell axes. Units must be in Å.

# Arguments
- `f_to_c::Array{Float64, 2}`: fractional to Cartesian coordinate transformation matrix.

# Returns
- `box::Box`: Fully formed Box object
"""
function construct_box(f_to_c::Array{Float64, 2})
    # unit cell volume (A³) is the determinant of the matrix
    Ω = det(f_to_c)

    # unit cell dimensions are lengths of the unit cell axes (Å)
    a = norm(f_to_c[:, 1])
    b = norm(f_to_c[:, 2])
    c = norm(f_to_c[:, 3])

    # c_to_f is the inverse
    c_to_f = inv(f_to_c)

    # angles (radians)
    α = acos(dot(f_to_c[:, 2], f_to_c[:, 3]) / (b * c))
    β = acos(dot(f_to_c[:, 1], f_to_c[:, 3]) / (a * c))
    γ = acos(dot(f_to_c[:, 1], f_to_c[:, 2]) / (a * b))

    # the columns of f_to_c are the unit cell axes
    r = reciprocal_lattice(f_to_c[:, 1], f_to_c[:, 2], f_to_c[:, 3])

    return Box(a, b, c, α, β, γ, Ω, f_to_c, c_to_f, r)
end

"""
    new_box = replicate_box(original_box, repfactors)

Replicates a `Box` in positive directions to construct a new `Box` representing a supercell.
The `original_box` is replicated according to the factors in `repfactors`.
Note `replicate_box(original_box, repfactors=(1, 1, 1))` returns same `Box`.
The new fractional coordinates as described by `f_to_c` and `c_to_f` still ∈ [0, 1].

# Arguments
- `original_box::Box`: The box that you want to replicate
- `repfactors::Tuple{Int, Int, Int}`: The factor you want to replicate the box by

# Returns
- `box::Box`: Fully formed Box object
"""
function replicate_box(box::Box, repfactors::Tuple{Int, Int, Int})
    #because this uses construct_box, its fractional coords still go 0 - 1
    return construct_box(box.a * repfactors[1], box.b * repfactors[2], box.c * repfactors[3],
                         box.α, box.β, box.γ)
end

"""
Data structure for a 3D crystal structure.

# Attributes
- `name::String`: corresponds to crystal structure filename from which it was read.
- `box::Box`: description of unit cell (Bravais Lattice); see `Box` struct.
- `n_atoms::Int64`: number of atoms in the unit cell
- `atoms::Array{Symbol,1}`: list of (pseudo)atoms (e.g. elements) composing crystal unit cell, in strict order.
- `xf::Array{Float64,2}`: fractional coordinates of the atoms, in strict order corresponding to `atoms`, stored column-wise so that `xf[:, i]` possesses the fractional coordinates of atom `i`.
- `charges::Array{Float64,1}`: the point charges of the atoms in corresponding order as `atoms`.
"""
struct Framework
    name::String

    box::Box

    n_atoms::Int64
    atoms::Array{Symbol, 1}
    xf::Array{Float64, 2}
    charges::Array{Float64, 1}
end

"""
    framework = read_crystal_structure_file("filename.cssr"; run_checks=true, net_charge_tol=0.001, remove_overlap = false)

Read a crystal structure file (.cif or .cssr) and construct a Framework object.
If `run_checks=True`, checks for atom overlap and charge neutrality. `net_charge_tol` is the tolerance for net charge.

# Arguments
- `filename::AbstractString`: the name of the crystal structure file
- `run_checks::Bool`: Will run overlap check and charge neutrality check if `true`
- `net_charge_tol::Float64`: Charge tolerance for charge neutrality check
- `remove_overlap::Bool`: Will remove "identical" atoms if true. Identical atoms are of the same element, occupying the same space.

# Returns
- `framework::Framework`: A framework containing the crystal structure information
"""
function read_crystal_structure_file(filename::AbstractString; run_checks::Bool=true, net_charge_tol::Float64=0.001, remove_overlap = false)
    # Read file extension. Ensure we can read the file type
    extension = split(filename, ".")[end]
    if ! (extension in ["cif", "cssr"])
        error("PorousMaterials.jl can only read .cif or .cssr crystal structure files.")
    end

    path_to_file = PATH_TO_DATA * "crystals/" * filename
    # Read in crystal structure file
    if ! isfile(path_to_file)
        error(@sprintf("Could not open crystal structure file %s\n",
                       path_to_file))
    end

    f = open(path_to_file, "r")
    lines = readlines(f)
    close(f)

    # Initialize arrays. We'll populate them when reading through the crystal structure file.
    charges = Array{Float64, 1}()
    xf = Array{Float64, 1}()
    yf = Array{Float64, 1}()
    zf = Array{Float64, 1}()
    atoms = Array{Symbol, 1}()


    # Start of .cif reader
    if extension == "cif"
        data = Dict{AbstractString, Float64}()
        loop_starts = -1
        for (i, line) in enumerate(lines)
            line = split(line)
            # Skip empty lines
            if length(line) == 0
                continue
            end

            # Make sure the space group is P1
            if line[1] == "_symmetry_space_group_name_H-M"
                if length(line) == 3
                    @assert(contains(line[2] * line[3], "P1") ||
                            contains(line[2] * line[3], "P 1") ||
                            contains(line[2] * line[3], "P-1") ||
                            contains(line[2] * line[3], "-P1"),
                            "cif must have P1 symmetry.\n")
                elseif length(line) == 2
                    @assert(contains(line[2], "P1") ||
                            contains(line[2], "P 1") ||
                            contains(line[2], "P -1") ||
                            contains(line[2], "-P1"),
                            "cif must have P1 symmetry.\n")
                else
                    println(line)
                    error("Does this .cif have P1 symmetry? Use `convert_cif_to_P1_symmetry` to convert to P1 symmetry")
                end
            end

            for axis in ["a", "b", "c"]
                if line[1] == @sprintf("_cell_length_%s", axis)
                    data[axis] = parse(Float64, split(line[2],'(')[1])
                end
            end

            for angle in ["alpha", "beta", "gamma"]
                if line[1] == @sprintf("_cell_angle_%s", angle)
                    data[angle] = parse(Float64, split(line[2],'(')[1]) * pi / 180.0
                end
            end

            # As soon as we reach the coordinate loop, we break this for-loop
            # and replace it with a while-loop further down
            if line[1] == "loop_"
                next_line = split(lines[i+1])
                if contains(next_line[1], "_atom_site")
                    loop_starts = i + 1
                    break
                end
            end

        end # End loop over lines

        if loop_starts == -1
            error("Could not find _atom_site* after loop_ in .cif file\n")
        end


        atom_column_name = ""
        # name_to_column is a dictionary that e.g. returns which column contains x fractional coord
        #   use example: name_to_column["_atom_site_fract_x"] gives 3
        name_to_column = Dict{AbstractString, Int}()

        i = loop_starts
        while length(split(lines[i])) == 1
            if i == loop_starts
                atom_column_name = split(lines[i])[1]
            end
            name_to_column[split(lines[i])[1]] = i + 1 - loop_starts
            i += 1
        end


        for i = loop_starts+length(name_to_column):length(lines)
            # Skip identical atoms if they're present

            line = split(lines[i])
            if length(line) != length(name_to_column)
                break
            end

            push!(atoms, line[name_to_column[atom_column_name]])
            push!(xf, mod(parse(Float64, line[name_to_column["_atom_site_fract_x"]]), 1.0))
            push!(yf, mod(parse(Float64, line[name_to_column["_atom_site_fract_y"]]), 1.0))
            push!(zf, mod(parse(Float64, line[name_to_column["_atom_site_fract_z"]]), 1.0))
            # If charges present, import them
            if haskey(name_to_column, "_atom_site_charge")
                push!(charges, parse(Float64, line[name_to_column["_atom_site_charge"]]))
            else
                push!(charges, 0.0)
            end
        end

        n_atoms = length(xf)
        a = data["a"]
        b = data["b"]
        c = data["c"]
        α = data["alpha"]
        β = data["beta"]
        γ = data["gamma"]

    # Start of cssr reader #TODO make sure this works for different .cssr files!
    elseif extension == "cssr"
        # First line contains unit cell lenghts
        line = split(lines[1])
        a = parse(Float64, line[1])
        b = parse(Float64, line[2])
        c = parse(Float64, line[3])

        # Second line contains unit cell angles
        line = split(lines[2])
        α = parse(Float64, line[1]) * pi / 180.0
        β = parse(Float64, line[2]) * pi / 180.0
        γ = parse(Float64, line[3]) * pi / 180.0

        n_atoms = parse(Int, split(lines[3])[1])

        # Read in atoms and fractional coordinates
        for i = 1:n_atoms
            line = split(lines[4 + i])
            push!(atoms, line[2])

            push!(xf, mod(parse(Float64, line[3]), 1.0)) # Wrap to [0,1]
            push!(yf, mod(parse(Float64, line[4]), 1.0)) # Wrap to [0,1]
            push!(zf, mod(parse(Float64, line[5]), 1.0)) # Wrap to [0,1]

            push!(charges, parse(Float64, line[14]))
        end
    end

    # Construct the unit cell box
    box = construct_box(a, b, c, α, β, γ)

    fractional_coords = Array{Float64, 2}(3, n_atoms)
    fractional_coords[1, :] = xf[:]; fractional_coords[2,:] = yf[:]; fractional_coords[3,:] = zf[:]

    # And finally construct the framework
    framework = Framework(filename, box, n_atoms, atoms, fractional_coords, charges)

    if run_checks && !remove_overlap
        if atom_overlap(framework)
            error(@sprintf("At least one pair of atoms overlap in %s\nConsider calling this function again with `remove_overlap = true`\n", framework.name))
        end
        if ! charge_neutral(framework, net_charge_tol)
            error(@sprintf("Framework %s is not charge neutral; sum of charges is : %f e\n", framework.name, sum(framework.charges)))
        end
    end

    if remove_overlap
        framework = remove_overlapping_atoms(framework, run_checks = run_checks)
    end

    return framework
end


"""
    replicate_to_xyz(framework, xyzfilename=nothing; comment="", repfactors=(0, 0, 0),
                     negative_replications=false)

Write a .xyz file representation of a crystal structure from a `Framework` type.
Write an optional `comment` to the .xyz file if desired.
Replicate the structure in the x-,y- and/or z-direction by changing the tuple `repfactors`.
A value of 1 replicates the structure once in the desired direction, so
`repfactors=(0, 0, 0)` includes only the "home" unit cell. Pass `negative_replications=true`
if home unit cell should be replicated in the negative directions too.

# Arguments
- `framework::Framework`: The framework containing the crystal structure information
- `xyzfilename::Union{AbstractString, Void}`: Name of the output file. If left blank, it will be named using the framework's name
- `comment::AbstractString`: An optional comment for the xyz file
- `repfactors::Tuple{Int, Int, Int}`: The replication factors used for the xyz file
- `negative_replications::Bool`: If true, the function will replicate the framework in both directions
"""
function replicate_to_xyz(framework::Framework, xyzfilename::Union{AbstractString, Void}=nothing;
                          comment::AbstractString="", repfactors::Tuple{Int, Int, Int}=(1, 1, 1),
                          negative_replications::Bool=false)
    # pre-calculate # of total atoms in .xyz
    if negative_replications
        n_atoms = 2 ^ 3 * framework.n_atoms * repfactors[1] * repfactors[2] * repfactors[3]
        neg_repfactors = (-repfactors[1] + 1, -repfactors[2] + 1, -repfactors[3] + 1)
    else
        n_atoms = framework.n_atoms * repfactors[1] * repfactors[2] * repfactors[3]
        neg_repfactors = (1, 1, 1)
    end

    # if no filename given, use framework's name
    if xyzfilename == nothing
        xyzfilename = split(framework.name, ".")[1] * ".xyz"
    end

    if ! contains(xyzfilename, ".xyz")
        xyzfilename *= ".xyz"
    end

    xyzfile = open(xyzfilename, "w")
    @printf(xyzfile, "%d\n%s\n", n_atoms, comment)

    for i = neg_repfactors[1]:repfactors[1], j = neg_repfactors[2]:repfactors[2], k = neg_repfactors[3]:repfactors[3]
        for a = 1:framework.n_atoms
            xf = framework.xf[:, a] + [i - 1.0, j - 1.0, k - 1.0]
            x = framework.box.f_to_c * xf
			@printf(xyzfile, "%s\t%.4f\t%.4f\t%.4f\n", string(framework.atoms[a]), x[1], x[2], x[3])
        end
    end
    close(xyzfile)

    println("See ", xyzfilename)
    return
end # replicate_to_xyz end


"""
    is_overlap = atom_overlap(framework; hard_diameter=0.1, verbose=false)

Return true iff any two atoms in the crystal overlap by calculating the distance
between every pair of atoms and ensuring distance is greater than
`hard_diameter`. If verbose, print the pair of atoms which are culprits.
# TODO include different radii for different atoms?

# Arguments
- `framework::Framework`: The framework containing the crystal structure information
- `hard_diameter::Float64`: The minimum distance between two atoms without them overlapping
- `verbose:Bool`: If true, will print out extra information as it's running

# Returns
- `is_overlap::Bool`: A Boolean telling us if any two atoms in the framework are overlapping
"""
function atom_overlap(framework::Framework; hard_diameter::Float64=0.1, verbose::Bool=false)
    overlap = false
    # loop over pairs of atoms
    for i = 1:framework.n_atoms
        for j = (i + 1):framework.n_atoms
            dxf = framework.xf[:, i] - framework.xf[:, j]
            # vector pointing from atom j to atom i in cartesian coords
            dx = framework.box.f_to_c * dxf
            if norm(dx) < hard_diameter
                overlap = true
                if verbose
                    @sprintf("atoms %d and %d are too close, distance %f Å < %f Å threshold\n", i, j, norm(dx), hard_diameter)
                end
            end
        end
    end
    return overlap
end


"""
    new_framework = remove_overlapping_atoms(framework; hard_diameter=0.1, verbose=false, run_checks = true, net_charge_tol=0.001)

Takes in a framework and checks to see if there are any "identical" atoms and promptly removes them. Identical atoms are two atoms of the same element, occupying the same space.

# Arguments
- `framework::Framework`: The framework containing the crystal structure information
- `hard_diameter::Float64`: The minimum distance between two atoms without them overlapping
- `run_checks::Bool`: Will run overlap check and charge neutrality check if `true` on new framework
- `net_charge_tol::Float64`: Charge tolerance for charge neutrality check

# Returns
- `new_framework::Framework`: A new framework where identical atoms have been removed.
"""
function remove_overlapping_atoms(framework::Framework; hard_diameter::Float64=0.1, run_checks = true, net_charge_tol::Float64=0.001)
    if !atom_overlap(framework)
        @printf("No atoms were overlapping.\nReturning original framework...\n")
        return framework
    end
    atoms_to_keep = trues(framework.n_atoms)

    for i = 1:framework.n_atoms
        for j = (i + 1):framework.n_atoms
            dxf = framework.xf[:, i] - framework.xf[:, j]
            dx = framework.box.f_to_c * dxf

            if norm(dx) < hard_diameter && framework.atoms[i] == framework.atoms[j]
                atoms_to_keep[j] = false
            end
        end
    end
    n_atoms = sum(atoms_to_keep)

    new_framework = Framework(framework.name, framework.box,
                          n_atoms, framework.atoms[atoms_to_keep],
                          framework.xf[:,atoms_to_keep], framework.charges[atoms_to_keep])

    if run_checks
        if atom_overlap(new_framework)
            error(@sprintf("At least one pair of atoms still overlaps in %s after we removed identical atoms. Consider checking if different elements are overlapping.\n", new_framework.name))
        end
        if ! charge_neutral(new_framework, net_charge_tol)
            error(@sprintf("Framework %s is not charge neutral; sum of charges is : %f e\n For comparison, the sum of charges for the old framework was : %f e\n", new_framework.name, sum(new_framework.charges), sum(framework.charges)))
        end
    end
    return new_framework
end


"""
    cn = charge_neutral(framework, net_charge_tol)

# Arguments
- `framework::Framework`: The framework containing the crystal structure information
- `net_charge_tol::Float64`: Charge tolerance for charge neutrality check

# Returns
- `cn::Bool`: True iff the framework is charge neutral, i.e. if the absolute value of the sum of the point charges assigned to its atoms is less than `net_charge_tol`. False otherwise.
"""
function charge_neutral(framework::Framework, net_charge_tol::Float64)
    sum_of_charges = sum(framework.charges)
    if abs(sum_of_charges) > net_charge_tol
        return false
    else
        return true
    end
end

"""
    charged_flag = charged(framework)

# Arguments
- `framework::Framework`: The framework containing the crystal structure information

# Returns
- `charged_flag::Bool`: True iff any of the atoms of the framework are assigned nonzero point charge. False otherwise.
"""
function charged(framework::Framework)
    zero_charge_flags = isapprox.(abs.(framework.charges), 0.0)
    return ! all(zero_charge_flags)
end

"""
    strip_numbers_from_atom_labels!(framework)

Strip numbers from labels for `framework.atoms`.
Precisely, for `atom` in `framework.atoms`, find the first number that appears in `atom`. Remove this number and all following characters from `atom`.
e.g. C12 --> C
	 Ba12A_3 --> Ba

# Arguments
- `framework::Framework`: The framework containing the crystal structure information
"""
function strip_numbers_from_atom_labels!(framework::Framework)
	for i = 1:framework.n_atoms
		atom_string = string(framework.atoms[i])
		for j = 1:length(atom_string)
			if !isalpha(atom_string[j])
				atom_string = atom_string[1:j-1]
				break
			end
		end
		framework.atoms[i] = Symbol(atom_string)
	end
	return
end

"""
    write_unitcell_boundary_vtk(box, filename)

Write unit cell boundary as a .vtk file for visualizing the unit cell boundary.
#TODO Is this function working properly?

# Arguments
- `box::Box`: The data structure describing the unit cell box (Bravais lattice)
- `filename::AbstractString`: The filename of the .vtk file
"""
function write_unitcell_boundary_vtk(box::Box, filename::AbstractString)
    # if no filename given, use framework's name
    if filename == nothing
        filename = split(framework.name, ".")[1] * ".vtk"
    end

    vtk_file = open(filename, "w")

    @printf(vtk_file, "# vtk DataFile Version 2.0\nunit cell boundary\n
                       ASCII\nDATASET POLYDATA\nPOINTS 8 double\n")

    # write points on boundary of unit cell
    for i = 0:1
        for j = 0:1
            for k = 0:1
                xf = [i, j, k] # fractional coordinates of corner
                cornerpoint = box.f_to_c * xf
                @printf(vtk_file, "%.3f %.3f %.3f\n",
                        cornerpoint[1], cornerpoint[2], cornerpoint[3])
            end
        end
    end

    # define connections
    @printf(vtk_file, "LINES 12 36\n2 0 1\n2 0 2\n2 1 3\n2 2 3\n2 4 5\n2 4 6\n2 5 7\n2 6 7\n2 0 4\n2 1 5\n2 2 6\n2 3 7\n")
    close(vtk_file)
    println("See ", filename)
    return
end

write_unitcell_boundary_vtk(framework::Framework) = write_unitcell_boundary_vtk(framework.box, split(framework.name, ".")[1] * ".vtk")

"""
    formula = chemical_formula(framework)

Find the irreducible chemical formula of a crystal structure.

# Arguments
- `framework::Framework`: The framework containing the crystal structure information

# Returns
- `formula::Dict{Symbol, Int}`: A dictionary with the irreducible chemical formula of a crystal structure
"""
function chemical_formula(framework::Framework; verbose::Bool=false)
    # use dictionary to count atom types
    atom_counts = Dict(zip(unique(framework.atoms), zeros(Int, length(unique(framework.atoms)))))
    for i = 1:framework.n_atoms
        atom_counts[framework.atoms[i]] += 1
    end

    # get greatest common divisor
    gcd_ = gcd([k for k in values(atom_counts)]...)

    # turn into irreducible chemical formula
    for atom in keys(atom_counts)
        atom_counts[atom] = atom_counts[atom] / gcd_
    end

    # print result
    if verbose
        @printf("Chemical formula of %s:\n\t", framework.name)
        for (atom, formula_unit) in atom_counts
			@printf("%s_%d", string(atom), formula_unit)
        end
        @printf("\n")
    end

    return atom_counts
end

"""
    convert_cif_to_P1_symmetry(filename, outputfilename; verbose=true)

Use Atomic Simulation Environment (ASE) Python package to convert .cif file in non-P1
symmetry to P1 symmetry. Writes .cif with P1 symmetry to `outputfilename1`.
Filenames correspond to files in `PATH_TO_DATA/crystals`.
# TODO add to README that ase is required. Later look at Crystals.jl

# Arguments
- `filename::AbstractString`: The name of the input file
- `outputfilename::AbstractString`: The name of the output file
- `verbose::Bool`: Will print updates as the function's running if true.
"""
function convert_cif_to_P1_symmetry(filename::AbstractString, outputfilename::AbstractString; verbose::Bool=true)
    # Import Atomic Simulation Environment Python package
    @pyimport ase
    @pyimport ase.io as aseio
    @pyimport ase.build as asebuild

    non_p1_cif_location = PATH_TO_DATA * "crystals/" * filename
    println(non_p1_cif_location)
    non_p1_cif = aseio.read(non_p1_cif_location, format="cif")

    p1_cif = asebuild.make_supercell(non_p1_cif, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    p1_cif_location = PATH_TO_DATA * "crystals/" * outputfilename
    aseio.write(p1_cif_location, p1_cif, format="cif")

    if verbose
        @printf("Converting to P1 symmetry using ASE.\n\t%s\n\t\t--->\n\t%s\n\n", non_p1_cif_location, p1_cif_location)
    end

    return
end


"""
    atomic_masses = read_atomic_masses()

Read the `data/atomicmasses.csv` file to construct a dictionary of atoms and their atomic
masses in amu.

# Returns
- `atomic_masses::Dict{Symbol, Float64}`: A dictionary containing the atomic masses of each atom stored in `data/atomicmasses.csv`
"""
function read_atomic_masses()
    if ! isfile("data/atomicmasses.csv")
        error("Cannot find atomicmasses.csv file in your data folder\n")
    end

    df_am = CSV.read("data/atomicmasses.csv")

    atomic_masses = Dict{Symbol, Float64}()

    for row in eachrow(df_am)
		atomic_masses[Symbol(row[:atom])] = row[:mass]
    end

    return atomic_masses
end

"""

    mass_of_framework = molecular_weight(framework)

Calculates the molecular weight of a unit cell of the framework in amu using information stored in `data/atomicmasses.csv`.

# Arguments
- `framework::Framework`: The framework containing the crystal structure information

# Returns
- `mass_of_framework::Float64`: The molecular weight of a unit cell of the framework in amu
"""
function molecular_weight(framework::Framework)
    mass = 0.0
    atomic_masses = read_atomic_masses()

	for atom in framework.atoms
        mass += atomic_masses[atom]
    end

    return mass #amu
end

"""
    ρ = crystal_density(framework) # kg / m3

Compute the crystal density of a framework in units kg/m3.

# Arguments
- `framework::Framework`: The framework containing the crystal structure information

# Returns
- `ρ::Float64`: The crystal density of a framework in kg/m³
"""
function crystal_density(framework::Framework)
    mw = molecular_weight(framework)
    return mw / framework.box.Ω * 1660.53892  # --> kg/m3
end

"""
    com = center_of_mass(framework::Framework; fractional_coordinates=false)
    com = center_of_mass(molecule::Molecule)

Computes the center of mass of a crystal or molecule. Atomic masses are retreived from
[`read_atomic_masses`](@ref).

For the `Framework`, periodic boundaries are taken into consideration. The unit cell
fractional coordinates are mapped to a circle, then the periodic nature of the circle is
used to calculate the center of mass.

# Arguments
- `framework::Framework`: The framework containing the crystal structure information
- `molecule::Molecule`: A molecule object
- `fractional_coordinates::Bool`: optional argument to return center of mass in fractional
coordinates instead of Cartesian. (Not an option for `Molecule`s.)

# Returns
- `com::Array{Float64, 1}`: center of mass coordinates
"""
function center_of_mass(frame::Framework; fractional_coordinates::Bool=false)
    atomic_masses = read_atomic_masses()
    xf_com = Array{Float64, 1}(3) # center of mass in fractional coordinates
    for axis in 1:3
        ξ_bar = 0.0
        ζ_bar = 0.0
        total_mass = 0.0
        for i in 1:frame.n_atoms
            θ = frame.xf[axis, i] * 2 * pi
            ξ = cos(θ) * atomic_masses[frame.atoms[i]] / (2 * pi)
            ζ = sin(θ) * atomic_masses[frame.atoms[i]] / (2 * pi)
            ξ_bar += ξ
            ζ_bar += ζ
            total_mass += atomic_masses[frame.atoms[i]]
        end
        ξ_bar /= total_mass
        ζ_bar /= total_mass
        θ_bar = atan2(-ζ_bar, -ξ_bar) + pi
        xf_com[axis] = θ_bar / (2 * pi)
    end
    if fractional_coordinates
        return xf_com
    end
    return frame.box.f_to_c * xf_com
end

function Base.show(io::IO, framework::Framework)
	@printf(io, "a = %.3f Angstrom\n", framework.box.a)
	@printf(io, "b = %.3f Angstrom\n", framework.box.b)
	@printf(io, "c = %.3f Angstrom\n", framework.box.c)

	@printf(io, "α = %.3f radians\n", framework.box.α)
	@printf(io, "β = %.3f radians\n", framework.box.β)
	@printf(io, "γ = %.3f radians\n", framework.box.γ)

	@printf(io, "Ω = %.3f Angstrom³\n", framework.box.Ω)

	@printf(io, "Number of atoms = %d", framework.n_atoms)
end


function Base.show(io::IO, box::Box)
    println(io, "Bravais unit cell of a crystal.")
    @printf(io, "\tUnit cell angles α = %f deg. β = %f deg. γ = %f deg.\n",
        box.α * 180.0 / π, box.β * 180.0 / π, box.γ * 180.0 / π)
    @printf(io, "\tUnit cell dimensions a = %f Å. b = %f Å, c = %f Å\n",
        box.a, box.b, box.c)
    @printf(io, "\tVolume of unit cell: %f Å³\n", box.Ω)
end

function Base.isapprox(box1::Box, box2::Box; rtol::Real=sqrt(eps()))
    return (isapprox(box1.a, box2.a, rtol=rtol) &&
            isapprox(box1.b, box2.b, rtol=rtol) &&
            isapprox(box1.c, box2.c, rtol=rtol) &&
            isapprox(box1.α, box2.α, rtol=rtol) &&
            isapprox(box1.β, box2.β, rtol=rtol) &&
            isapprox(box1.γ, box2.γ, rtol=rtol) &&
            isapprox(box1.Ω, box2.Ω, rtol=rtol) &&
            isapprox(box1.f_to_c, box2.f_to_c, rtol=rtol) &&
            isapprox(box1.c_to_f, box2.c_to_f, rtol=rtol) &&
            isapprox(box1.reciprocal_lattice, box2.reciprocal_lattice, rtol=rtol))
end

function Base.isapprox(f1::Framework, f2::Framework; checknames::Bool=false)
    names = f1.name == f2.name
    box = isapprox(f1.box, f2.box)
    n_atoms = f1.n_atoms == f2.n_atoms
    if checknames && n_atoms
        return names && box && n_atoms && (f1.atoms == f2.atoms) && isapprox(f1.xf, f2.xf) && isapprox(f1.charges, f2.charges)
    else
        return box && n_atoms && (f1.atoms == f2.atoms) && isapprox(f1.xf, f2.xf) && isapprox(f1.charges, f2.charges)
    end
end
