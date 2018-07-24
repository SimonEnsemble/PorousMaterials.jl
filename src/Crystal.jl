using Base.Test

"""
Data structure for a 3D crystal structure.

# Attributes
- `name::String`: name of crystal structure
- `box::Box`: unit cell (Bravais Lattice)
- `atoms::Array{LJSphere, 1}`: list of Lennard-Jones spheres in crystal unit cell
- `charges::Array{PtCharge, 1}`: list of point charges in crystal unit cell
"""
struct Framework
    name::String
    box::Box
    atoms::Array{LJSphere, 1}
    charges::Array{PtCharge, 1}
end

"""
    framework = Framework("filename.cssr"; run_checks=true, 
                          net_charge_tol=0.001, remove_overlap = false)

Read a crystal structure file (.cif or .cssr) and construct a Framework object.
If `run_checks=True`, checks for atom overlap and charge neutrality. `net_charge_tol` is the tolerance for net charge.

# Arguments
- `filename::AbstractString`: the name of the crystal structure file
- `check_charge_neutrality::Bool`: check for charge neutrality
- `net_charge_tol::Float64`: Charge tolerance for charge neutrality check
- `remove_overlap::Bool`: Will remove "identical" atoms if true. Identical atoms are of the same element, occupying the same space.

# Returns
- `framework::Framework`: A framework containing the crystal structure information
"""
function Framework(filename::AbstractString; check_charge_neutrality::Bool=true,
                   net_charge_tol::Float64=0.001, check_atom_overlap::Bool=true,
                   remove_overlap::Bool=false)
    # Read file extension. Ensure we can read the file type
    extension = split(filename, ".")[end]
    if ! (extension in ["cif", "cssr"])
        error("PorousMaterials.jl can only read .cif or .cssr crystal structure files.")
    end

    # read file
    f = open(PATH_TO_DATA * "crystals/" * filename, "r")
    lines = readlines(f)
    close(f)

    # Initialize arrays. We'll populate them when reading through the crystal structure file.
    charge_values = Array{Float64, 1}()
    xf = Array{Float64, 1}()
    yf = Array{Float64, 1}()
    zf = Array{Float64, 1}()
    species = Array{Symbol, 1}()

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
            
            # pick up unit cell lengths
            for axis in ["a", "b", "c"]
                if line[1] == @sprintf("_cell_length_%s", axis)
                    data[axis] = parse(Float64, split(line[2],'(')[1])
                end
            end
            
            # pick up unit cell angles
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
            line = split(lines[i])
            if length(line) != length(name_to_column)
                break
            end

            push!(species, line[name_to_column[atom_column_name]])
            push!(xf, mod(parse(Float64, line[name_to_column["_atom_site_fract_x"]]), 1.0))
            push!(yf, mod(parse(Float64, line[name_to_column["_atom_site_fract_y"]]), 1.0))
            push!(zf, mod(parse(Float64, line[name_to_column["_atom_site_fract_z"]]), 1.0))
            # If charges present, import them
            if haskey(name_to_column, "_atom_site_charge")
                push!(charge_values, parse(Float64, line[name_to_column["_atom_site_charge"]]))
            else
                push!(charge_values, 0.0)
            end
        end

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
            push!(species, line[2])

            push!(xf, mod(parse(Float64, line[3]), 1.0)) # Wrap to [0,1]
            push!(yf, mod(parse(Float64, line[4]), 1.0)) # Wrap to [0,1]
            push!(zf, mod(parse(Float64, line[5]), 1.0)) # Wrap to [0,1]

            push!(charge_values, parse(Float64, line[14]))
        end
    end

    # Construct the unit cell box
    box = Box(a, b, c, α, β, γ)

    atoms = LJSphere[]
    charges = PtCharge[]
    for a = 1:length(species)
        frac_coord = [xf[a], yf[a], zf[a]]
        push!(atoms, LJSphere(species[a], frac_coord))
        if abs(charge_values[a]) > 0.0
            push!(charges, PtCharge(charge_values[a], frac_coord))
        end
    end

    framework = Framework(filename, box, atoms, charges)

    if check_charge_neutrality
        if ! charge_neutral(framework, net_charge_tol)
            error(@sprintf("Framework %s is not charge neutral; net charge is %f e. Ignore 
            this error message by passing check_charge_neutrality=false or increasing the
            net charge tolerance `net_charge_tol`\n", 
                            framework.name, total_charge(framework)))
        end
    end

    if remove_overlap
        return remove_overlapping_atoms(framework)
    end

    if check_atom_overlap
        if atom_overlap(framework)
            error(@sprintf("At least one pair of atoms overlap in %s\nConsider passing `remove_overlap=true`\n", framework.name))
        end
    end

    return framework
end

"""
    replicated_frame = replicate(framework, repfactors)

Replicates the atoms and charges in a `Framework` in positive directions to 
construct a new `Framework`. Note `replicate(framework, (1, 1, 1))` returns the same `Framework`.

# Arguments
- `framework::Framework`: The framework to replicate
- `repfactors::Tuple{Int, Int, Int}`: The factors by which to replicate the crystal structure in each direction.

# Returns
- `replicated_frame::Framework`: Replicated framework
"""
function replicate(framework::Framework, repfactors::Tuple{Int, Int, Int})
    # determine number of atoms in replicated framework
    n_atoms = length(framework.atoms) * repfactors[1] * repfactors[2] * repfactors[3]

    # replicate box
    new_box = replicate(framework.box, repfactors)

    # replicate atoms and charges
    new_charges = PtCharge[]
    new_atoms = LJSphere[]
    for ra = 0:(repfactors[1] - 1), rb = 0:(repfactors[2] - 1), rc = 0:(repfactors[3] - 1)
        for atom in framework.atoms
            xf = atom.xf + 1.0 * [ra, rb, rc]
            # scale fractional coords
            xf = xf ./ repfactors
            push!(new_atoms, LJSphere(atom.species, xf))
        end
        for charge in framework.charges
            xf = charge.xf + 1.0 * [ra, rb, rc]
            # scale fractional coords
            xf = xf ./ repfactors
            push!(new_charges, PtCharge(charge.q, xf))
        end
    end
    @assert(length(new_charges) == length(framework.charges) * prod(repfactors))
    @assert(length(new_atoms) == length(framework.atoms) * prod(repfactors))
    return Framework(framework.name, new_box, new_atoms, new_charges)
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
- `repfactors::Tuple{Int, Int, Int}`: How many times to replicate the framework in each direction.
- `negative_replications::Bool`: If true, the function will replicate the framework in both directions
"""
function replicate_to_xyz(framework::Framework, xyzfilename::Union{AbstractString, Void}=nothing;
                          comment::AbstractString="", repfactors::Tuple{Int, Int, Int}=(1, 1, 1),
                          negative_replications::Bool=false)
    # pre-calculate # of total atoms in .xyz
    if negative_replications
        n_atoms = 2 ^ 3 * length(framework.atoms) * repfactors[1] * repfactors[2] * repfactors[3]
        neg_repfactors = (-repfactors[1] + 1, -repfactors[2] + 1, -repfactors[3] + 1)
    else
        n_atoms = length(framework.atoms) * repfactors[1] * repfactors[2] * repfactors[3]
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
        for atom in framework.atoms
            xf = atom.xf + [i - 1.0, j - 1.0, k - 1.0]
            x = framework.box.f_to_c * xf
			@printf(xyzfile, "%s\t%.4f\t%.4f\t%.4f\n", string(atom.species), x[1], x[2], x[3])
        end
    end
    close(xyzfile)

    println("See ", xyzfilename)
    return
end # replicate_to_xyz end


"""
    is_overlap = atom_overlap(framework; hard_diameter=0.1, verbose=false)

Return true iff any two `LJSphere`'s in the crystal overlap by calculating the distance
between every pair of atoms and ensuring distance is greater than
`hard_diameter`. If verbose, print the pair of atoms which are culprits.

# Arguments
- `framework::Framework`: The framework containing the crystal structure information
- `hard_diameter::Float64`: The minimum distance between two atoms without them overlapping
- `verbose:Bool`: If true, will print out extra information as it's running

# Returns
- `overlap::Bool`: A Boolean telling us if any two atoms in the framework are overlapping
"""
function atom_overlap(framework::Framework; hard_diameter::Float64=0.1, verbose::Bool=false)
    overlap = false
    for (i, atom_i) in enumerate(framework.atoms)
        for (j, atom_j) in enumerate(framework.atoms)
            if j >= i
                continue
            end
            dxf = atom_i.xf - atom_j.xf
            nearest_image!(dxf)
            r = norm(framework.box.f_to_c * dxf)
            if r < hard_diameter
                overlap = true
                if verbose
                    @sprintf("The distance between atom %d and atom %d in %s are too close.
                              Distance %f Å < %f Å threshold\n", i, j, framework.name, r, 
                              hard_diameter)
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
function remove_overlapping_atoms(framework::Framework; hard_diameter::Float64=0.1)
    if !atom_overlap(framework)
        @printf("No atoms were overlapping.\nReturning original framework...\n")
        return framework
    end
    atoms_to_keep = trues(length(framework.atoms))

    if charged(framework) && (length(framework.atoms) != length(framework.charges))
        # TODO hv option to read in zero charges?
        error("length of charges not equal to length of atoms; must remove overlapping
        atoms manually. (some charges are zero...)")
    end
    
    for (i, atom_i) in enumerate(framework.atoms)
        for (j, atom_j) in enumerate(framework.atoms)
            if j >= i
                continue
            end

            dxf = atom_i.xf - atom_j.xf
            nearest_image!(dxf)
            r = norm(framework.box.f_to_c * dxf)

            if r < hard_diameter
                atoms_to_keep[i] = false
            end
        end
    end

    new_framework = Framework(framework.name, framework.box,
                              framework.atoms[atoms_to_keep],
                              charged(framework) ? framework.charges[atoms_to_keep] : PtCharge[])
    
    @assert (! atom_overlap(new_framework))
    return new_framework
end
 #

function total_charge(framework::Framework)
    q = 0.0
    for charge in framework.charges
        q += charge.q
    end
    return q
end

"""
    charge_neutral_flag = charge_neutral(framework, net_charge_tol) # true or false

Determine if the absolute value of the net charge in `framework` is less than `net_charge_tol`.
"""
function charge_neutral(framework::Framework, net_charge_tol::Float64)
    q = total_charge(framework)
    return abs(q) < net_charge_tol
end

"""
    charged_flag = charged(framework, verbose=false) # true or false

Determine if a framework has point charges
"""
function charged(framework::Framework; verbose::Bool=false)
    charged_flag = length(framework.charges) > 0
    if verbose
        @printf("\tFramework atoms of %s have charges? %s\n", framework.name, charged_flag)
    end
    return charged_flag
end

"""
    strip_numbers_from_atom_labels!(framework)

Strip numbers from labels for `framework.atoms`.
Precisely, for `atom` in `framework.atoms`, find the first number that appears in `atom`. 
Remove this number and all following characters from `atom`.
e.g. C12 --> C
	 Ba12A_3 --> Ba

# Arguments
- `framework::Framework`: The framework containing the crystal structure information
"""
function strip_numbers_from_atom_labels!(framework::Framework)
    for (a, atom) in enumerate(framework.atoms)
        # atom species in string format
		species = string(atom.species)
		for j = 1:length(species)
			if ! isalpha(species[j])
                framework.atoms[a] = LJSphere(species[1:j-1], atom.xf)
				break
			end
		end
	end
    return
end

"""
    write_unitcell_boundary_vtk(box, filename, verbose=true)

Write unit cell boundary as a .vtk file for visualizing the unit cell boundary.

Appends ".vtk" extension to `filename` automatically if not passed.

# Arguments
- `box::Box`: the unit cell box (Bravais lattice)
- `filename::AbstractString`: filename of the .vtk file output (absolute path)
"""
function write_unitcell_boundary_vtk(box::Box, filename::AbstractString; verbose::Bool=true)
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
    if verbose
        println("See ", filename)
    end
    return
end

write_unitcell_boundary_vtk(framework::Framework) = write_unitcell_boundary_vtk(framework.box, split(framework.name, ".")[1])

"""
    formula = chemical_formula(framework)

Find the irreducible chemical formula of a crystal structure.

# Arguments
- `framework::Framework`: The framework containing the crystal structure information

# Returns
- `formula::Dict{Symbol, Int}`: A dictionary with the irreducible chemical formula of a crystal structure
"""
function chemical_formula(framework::Framework; verbose::Bool=false)
    unique_atoms = unique([atom.species for atom in framework.atoms])
    # use dictionary to count atom types
    atom_counts = Dict{Symbol, Int}([a => 0 for a in unique_atoms])
    for atom in framework.atoms
        atom_counts[atom.species] += 1
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

    mass_of_framework = molecular_weight(framework)

Calculates the molecular weight of a unit cell of the framework in amu using information stored in `data/atomicmasses.csv`.

# Arguments
- `framework::Framework`: The framework containing the crystal structure information

# Returns
- `mass_of_framework::Float64`: The molecular weight of a unit cell of the framework in amu
"""
function molecular_weight(framework::Framework)
    atomic_masses = read_atomic_masses()

    mass = 0.0
	for atom in framework.atoms
        mass += atomic_masses[atom.species]
    end

    return mass # amu
end

"""
    ρ = crystal_density(framework) # kg/m²

Compute the crystal density of a framework. Pulls atomic masses from [`read_atomic_masses`](@ref).

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
    write_cif(framework, filename)

Write a `framework::Framework` to a .cif file with `filename::String`. If `filename` does
not include the .cif extension, it will automatically be added.
"""
function write_cif(framework::Framework, filename::String)
    if charged(framework) && (length(framework.atoms) != length(framework.charges))
        error("write_cif assumes equal numbers of Charges and LJSpheres (or zero charges)")
    end
    # append ".cif" to filename if it doesn't already have the extension
    if ! contains(filename, ".cif")
        filename *= ".cif"
    end
    cif_file = open(filename, "w")
    @printf(cif_file, "_symmetry_space_group_name_H-M   'P 1'\n")

    @printf(cif_file, "_cell_length_a %f\n", framework.box.a)
    @printf(cif_file, "_cell_length_b %f\n", framework.box.b)
    @printf(cif_file, "_cell_length_c %f\n", framework.box.c)

    @printf(cif_file, "_cell_angle_alpha %f\n", framework.box.α * 180.0 / pi)
    @printf(cif_file, "_cell_angle_beta %f\n", framework.box.β * 180.0 / pi)
    @printf(cif_file, "_cell_angle_gamma %f\n", framework.box.γ * 180.0 / pi)

    @printf(cif_file, "_symmetry_Int_Tables_number 1\n\n")
    @printf(cif_file, "loop_\n_symmetry_equiv_pos_as_xyz\n 'x, y, z'\n\n")

    @printf(cif_file, "loop_\n_atom_site_label\n")
    @printf(cif_file, "_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n")
    @printf(cif_file, "_atom_site_charge\n")

    for (a, atom) in enumerate(framework.atoms)
        q = 0.0
        if charged(framework)
            charge = framework.charges[a]
            q = charge.q
            if ! isapprox(charge.xf, atom.xf)
                error("write_cif assumes charges correspond to LJspheres")
            end
        end
        @printf(cif_file, "%s %f %f %f %f\n", atom.species, atom.xf..., q)
     end
     close(cif_file)
end

function Base.show(io::IO, framework::Framework)
    println(io, "Name: ", framework.name)
    println(io, framework.box)
	@printf(io, "Number of atoms = %d\n", length(framework.atoms))
	@printf(io, "Number of charges = %d\n", length(framework.charges))
    println(io, "Chemical formula: ", chemical_formula(framework))
end

function Base.isapprox(f1::Framework, f2::Framework; checknames::Bool=false)
    names_flag = f1.name == f2.name
    if checknames && (! names_flag)
        return false
    end
    box_flag = isapprox(f1.box, f2.box)
    if length(f1.charges) != length(f2.charges)
        return false
    end
    if length(f1.atoms) != length(f2.atoms)
        return false
    end
    charges_flag = all(isapprox.(f1.charges, f2.charges))
    atoms_flag = all(isapprox.(f1.atoms, f2.atoms))
    return box_flag && charges_flag && atoms_flag
end
