using Base.Test

# Data structure for a framework; user-friendly constructor below
struct Framework
    name::String
    box::Box
    atoms::Array{LJSphere, 1}
    charges::Array{PtCharge, 1}
end

"""
    framework = Framework(filename, check_charge_neutrality=true,
                          net_charge_tol=0.001, check_atom_and_charge_overlap=true,
                          remove_overlap=false)
    framework = Framework(name, box, atoms, charges)

Read a crystal structure file (.cif or .cssr) and populate a `Framework` data structure,
or construct a `Framework` data structure directly.

# Arguments
- `filename::AbstractString`: the name of the crystal structure file (include ".cif" or ".cssr") read from `PorousMaterials.PATH_TO_DATA * "structures/"`.
- `check_charge_neutrality::Bool`: check for charge neutrality
- `net_charge_tol::Float64`: when checking for charge neutrality, throw an error if the absolute value of the net charge is larger than this value.
- `check_atom_and_charge_overlap::Bool`: throw an error if overlapping atoms are detected.
- `remove_overlap::Bool`: remove identical atoms automatically. Identical atoms are the same element and overlap.

# Returns
- `framework::Framework`: A framework containing the crystal structure information

# Attributes
- `name::String`: name of crystal structure
- `box::Box`: unit cell (Bravais Lattice)
- `atoms::Array{LJSphere, 1}`: list of Lennard-Jones spheres in crystal unit cell
- `charges::Array{PtCharge, 1}`: list of point charges in crystal unit cell
"""
function Framework(filename::AbstractString; check_charge_neutrality::Bool=true,
                   net_charge_tol::Float64=0.001, check_atom_and_charge_overlap::Bool=true,
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
                            contains(line[2] * line[3], "-P1"),
                            "cif must have P1 symmetry.\n")
                elseif length(line) == 2
                    @assert(contains(line[2], "P1") ||
                            contains(line[2], "P 1") ||
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
        return remove_overlapping_atoms_and_charges(framework)
    end

    if check_atom_and_charge_overlap
        if atom_overlap(framework) | charge_overlap(framework)
            error(@sprintf("At least one pair of atoms/charges overlap in %s.
            Consider passing `remove_overlap=true`\n", framework.name))
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

# doc string in Misc.jl
function write_to_xyz(framework::Framework, filename::AbstractString;
                      comment::AbstractString="")
    atoms = [atom.species for atom in framework.atoms]
    x = zeros(Float64, 3, length(framework.atoms))
    for (a, atom) in enumerate(framework.atoms)
        x[:, a] = framework.box.f_to_c * atom.xf
    end

    write_to_xyz(atoms, x, filename, comment=comment)
end 


"""
    is_overlap = atom_overlap(framework; overlap_tol=0.1, verbose=false)

Return true iff any two `LJSphere`'s in the crystal overlap by calculating the distance
between every pair of atoms and ensuring distance is greater than
`overlap_tol`. If verbose, print the pair of atoms which are culprits.

# Arguments
- `framework::Framework`: The framework containing the crystal structure information
- `overlap_tol::Float64`: The minimum distance between two atoms without them overlapping
- `verbose:Bool`: If true, will print out extra information as it's running

# Returns
- `overlap::Bool`: A Boolean telling us if any two atoms in the framework are overlapping
"""
function atom_overlap(framework::Framework; overlap_tol::Float64=0.1, verbose::Bool=true)
    overlap = false
    for (i, atom_i) in enumerate(framework.atoms)
        for (j, atom_j) in enumerate(framework.atoms)
            if j >= i
                continue
            end
            if _overlap(atom_i, atom_j, framework.box, overlap_tol)
                overlap = true
                if verbose
                    warn(@sprintf("Atoms %d and %d in %s are less than %d Å apart.", i, j,
                        framework.name, overlap_tol))
                end
            end
        end
    end
    return overlap
end

function charge_overlap(framework::Framework; overlap_tol::Float64=0.1, verbose::Bool=true)
    overlap = false
    for (i, charge_i) in enumerate(framework.charges)
        for (j, charge_j) in enumerate(framework.charges)
            if j >= i
                continue
            end
            if _overlap(charge_i, charge_j, framework.box, overlap_tol)
                overlap = true
                if verbose
                    warn(@sprintf("Charges %d and %d in %s are less than %d Å apart.", i, j,
                        framework.name, overlap_tol))
                end
            end
        end
    end
    return overlap
end

# determine if two lennard-jones spheres overlap
function _overlap(a1::Union{PtCharge, LJSphere}, a2::Union{PtCharge, LJSphere}, box::Box, 
                  overlap_tol::Float64)
    dxf = a1.xf - a2.xf
    nearest_image!(dxf)
    if norm(box.f_to_c * dxf) < overlap_tol
        return true
    else
        return false
    end
end

#TODO write tests for this! one with diff elements
"""
    new_framework = remove_overlapping_atoms_and_charges(framework, overlap_tol=0.1, verbose=true)

Takes in a framework and returns a new framework with where overlapping atoms and overlapping
charges were removed. i.e. if there is an overlapping pair, one in the pair is removed. 
For any atoms or charges to be removed, the species and charge, respectively,
must be identical.

# Arguments
- `framework::Framework`: The framework containing the crystal structure information
- `atom_overlap_tol::Float64`: The minimum distance between two atoms that is tolerated
- `charge_overlap_tol::Float64`: The minimum distance between two charges that is tolerated

# Returns
- `new_framework::Framework`: A new framework where identical atoms have been removed.
"""
function remove_overlapping_atoms_and_charges(framework::Framework; 
    atom_overlap_tol::Float64=0.1, charge_overlap_tol::Float64=0.1, verbose::Bool=true)

    atoms_to_keep = trues(length(framework.atoms))
    charges_to_keep = trues(length(framework.atoms))

    for (i, atom_i) in enumerate(framework.atoms)
        for (j, atom_j) in enumerate(framework.atoms)
            if j >= i
                continue
            end
            if _overlap(atom_i, atom_j, framework.box, atom_overlap_tol)
                if atom_i.species != atom_j.species
                    error(@sprintf("Atom %d, %s and atom %d, %s overlap but are not the 
                    same element so we will not automatically remove one in the pair.\n", 
                    i, atom_i.species, j, atom_j.species))
                else
                    atoms_to_keep[i] = false
                end
            end
        end
    end
    if verbose
        println("# atoms removed: ", sum(.! atoms_to_keep))
    end

    for (i, charge_i) in enumerate(framework.charges)
        for (j, charge_j) in enumerate(framework.charges)
            if j >= i
                continue
            end
            if _overlap(charge_i, charge_j, framework.box, charge_overlap_tol)
                if ! isapprox(charge_i.q, charge_j.q)
                    error(@sprintf("charge %d of %f and charge %d of %f overlap but are 
                    not the same charge so we will not automatically remove one in the pair.\n", 
                    i, charge_i.q, j, charge_j.q))
                else
                    charges_to_keep[i] = false
                end
            end
        end
    end
    if verbose
        println("# charges removed: ", sum(.! charges_to_keep))
    end

    new_framework = Framework(framework.name, framework.box, 
        framework.atoms[atoms_to_keep],
        charged(framework) ? framework.charges[charges_to_keep] : PtCharge[])

    @assert (! atom_overlap(new_framework, overlap_tol=atom_overlap_tol))
    @assert (! charge_overlap(new_framework, overlap_tol=charge_overlap_tol))
    
    return new_framework
end

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

write_vtk(framework::Framework) = write_vtk(framework.box, split(framework.name, ".")[1])

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
    # first line should be data_xtalname_PM
    if framework.name == ""
        @printf(cif_file, "data_PM\n")
    else
        # don't include file extension!
        @printf(cif_file, "data_%s_PM\n", split(framework.name, ".")[1])
    end

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

"""
    new_framework = assign_charges(framework, charges, net_charge_tol=1e-5)

Assign charges to the atoms (represented by `LJSphere`s in `framework.atoms`) present in 
the framework. Pass a dictionary of charges that place charges according to the species
of the atoms or pass an array of charges to assign to each atom, with the order of the 
array consistent with the order of `framework.atoms`.

If the framework already has charges, the charges are removed and new charges are added
accordingly so that `length(framework.atoms) == length(framework.charges)`.

# Examples
```
charges = Dict(:Ca => 2.0, :C => 1.0, :H => -1.0)
new_framework = assign_charges(framework, charges)
```

```
charges = [4.0, 2.0, -6.0] # framework.atoms is length 3
new_framework = assign_charges(framework, charges)
```

# Arguments
- `framework::Framework`: the framework to which we should add charges (not modified in 
this function)
- `charges::Union{Dict{Symbol, Float64}, Array{Float64, 1}}`: a dictionary that returns the
charge assigned to the species of atom or an array of charges to assign, with order
consistent with the order in `framework.atoms` (units: electrons).
- `net_charge_tol::Float64`: the net charge tolerated when asserting charge neutrality of
the resulting framework

# Returns
- `new_framework::Framework`: a new framework identical to the one passed except charges 
are assigned.
"""
function assign_charges(framework::Framework, charges::Union{Dict{Symbol, Float64}, Array{Float64, 1}},
    net_charge_tol::Float64=1e-5)
    # if charges are already present, may make little sense to assign charges to atoms again
    if length(framework.charges) != 0
        warn(@sprintf("Charges are already present in %s. Removing the current charges on the
        framework and adding new ones...\n", framework.name))
    end
    
    # build the array of point charges according to atom species
    pt_charges = PtCharge[]
    for (a, atom) in enumerate(framework.atoms)
        if isa(charges, Dict{Symbol, Float64})
            if ! (atom.species in keys(charges))
                error(@sprintf("Atom %s is not present in the charge dictionary passed to 
                `assign_charges` for %s\n", atom.species, framework.name))
            end
            push!(pt_charges, PtCharge(charges[atom.species], atom.xf))
        else
            if length(charges) != length(framework.atoms)
                error(@sprintf("Length of `charges` array passed to `assign_charges` is not
                equal to the number of atoms in %s = %d\n", framework.name, 
                length(framework.atoms)))
            end
            push!(pt_charges, PtCharge(charges[a], atom.xf))
        end
    end

    # construct new framework
    new_framework = Framework(framework.name, framework.box, framework.atoms, pt_charges)

    # check for charge neutrality
    if abs(total_charge(new_framework)) > net_charge_tol
        error(@sprintf("Net charge of framework %s = %f > net charge tolerance %f. If
        charge neutrality is not a problem, pass `net_charge_tol=Inf`\n", framework.name,
        total_charge(new_framework), net_charge_tol))
    end

    return new_framework
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
