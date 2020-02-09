"""
    SymmetryInfo(symmetry, space_group, is_p1)

# Attributes
- `symmetry_ops::Array{Function, 2}`: 2D array of anonymous functions that represent
    the symmetry operations. If the structure is in P1 there will be one
    symmetry operation.
- `space_group::AbstractString`: The name of the space group. This is stored
    so that it can be written out again in the write_cif function. The space
    group is not used to verify the symmetry rules.
- `is_p1::Bool`: Stores whether the crystal is currently in P1 symmetry. 
"""
struct SymmetryInfo
    symmetry_ops::Array{String, 2}
    space_group::String
    is_p1::Bool
end
SymmetryInfo() = SymmetryInfo(["x", "y", "z"], "P1", true) # default

struct Crystal
    name::String
    box::Box
    atoms::Atoms{Frac}
    charges::Charges{Frac}
    bonds::SimpleGraph
    symmetry::SymmetryInfo
end

# default constructor without bond info or symmetry info
Crystal(name::String, box::Box, atoms::Atoms{Frac}, charges::Charges{Frac}) = Crystal(
    name, box, atoms, charges, SimpleGraph(atoms.n), SymmetryInfo())

"""
    crystal = Crystal(filename;
        check_neutrality=true, net_charge_tol=0.001, 
        check_overlap=true, overlap_tol=0.1,
        convert_to_p1=true, read_bonds_from_file=false, wrap_coords=true,
        include_zero_charges=false) # read from file

    crystal = Crystal(name, box, atoms, charges) # construct from matter, no bonds, P1-symmetry assumed

Read a crystal structure file (.cif or .cssr) and populate a `Crystal` data structure,
or construct a `Crystal` data structure directly.

# Arguments
- `filename::String`: the name of the crystal structure file (include ".cif" or ".cssr") read from `PATH_TO_CRYSTALS`.
- `check_neutrality::Bool`: check for charge neutrality
- `net_charge_tol::Float64`: when checking for charge neutrality, throw an error if the absolute value of the net charge is larger than this value.
- `check_overlap::Bool`: throw an error if overlapping atoms are detected.
- `convert_to_p1::Bool`: If the structure is not in P1 it will be converted to
    P1 symmetry using the symmetry rules from the `_symmetry_equiv_pos_as_xyz` list in the .cif file. 
    (We do not use the space groups name to look up symmetry rules).
- `read_bonds_from_file::Bool`: Whether or not to read bonding information from
    cif file. If false, the bonds can be inferred later. note that, if the crystal is not in P1 symmetry, we cannot *both* read bonds and convert to P1 symmetry.
- `wrap_coords::Bool`: if `true`, enforce that fractional coords of atoms and charges are in [0,1]³ by mod(x, 1)
- `include_zero_charges::Bool`: if `false`, do not include in `crystal.charges` atoms which have zero charges, in order to speed up the electrostatic calculations.
    If `true,` include the atoms in `crystal.charges` that have zero charge, ensuring that the number of atoms is equal to the number of charges and that `crystal.charges.xf` and crystal.atoms.xf` are the same.

# Returns
- `crystal::Crystal`: A crystal containing the crystal structure information

# Attributes
- `name::AbstractString`: name of crystal structure
- `box::Box`: unit cell (Bravais Lattice)
- `atoms::Atoms`: list of Atoms in crystal unit cell
- `charges::Charges`: list of point charges in crystal unit cell
- `bonds::SimpleGraph`: Unweighted, undirected graph showing all of the atoms
    that are bonded within the crystal
- `symmetry::SymmetryInfo`: symmetry inforomation
"""
function Crystal(filename::String; 
                 check_neutrality::Bool=true, net_charge_tol::Float64=0.001, 
                 check_overlap::Bool=true, overlap_tol::Float64=0.1,
                 convert_to_p1::Bool=true,
                 read_bonds_from_file::Bool=false, wrap_coords::Bool=true,
                 include_zero_charges::Bool=false)
    # Read file extension. Ensure we can read the file type
    extension = split(filename, ".")[end]
    if ! (extension in ["cif", "cssr"])
        error("I can only read .cif or .cssr crystal structure files.")
    end

    # read all lines of crystal structure file
    _f = open(joinpath(PATH_TO_CRYSTALS, filename), "r")
    lines = readlines(_f)
    close(_f)

    # Initialize arrays. We'll populate them when reading through the crystal structure file.
    charge_values = Array{Float64, 1}()
    species = Array{Symbol, 1}()
    xf = Array{Float64, 1}()
    yf = Array{Float64, 1}()
    zf = Array{Float64, 1}()
    coords = Array{Float64, 2}(undef, 3, 0)
    # default for symmetry rules is P1.
    # These will be overwritten if the user chooses to read in non-P1
    symmetry_ops = Array{String, 2}(undef, 3, 0)
    # creating empty SimpleGraph, might not have any information read in
    bonds = SimpleGraph()
    # used for remembering whether fractional/cartesian coordinates are read in
    # placed here so it will be defined for the if-stmt after the box is defined
    fractional = false
    cartesian = false
    # used for determining if the crystal is in P1 symmetry for simulations
    is_p1 = false
    space_group = ""


    # Start of .cif reader
    ###################################
    # CIF READER
    ###################################
    if extension == "cif"
        data = Dict{AbstractString, Float64}()
        loop_starts = -1
        i = 1
        # used for reading in symmetry options and replications
        symmetry_info = false
        atom_info = false
        label_num_to_idx = Dict{AbstractString, Int}()
        fractional = false
        cartesian = false
        while i <= length(lines)
            line = split(lines[i])
            # Skip empty lines
            if length(line) == 0
                i += 1
                continue
            end

            # Make sure the space group is P1
            if line[1] == "_symmetry_space_group_name_H-M"
                # use anonymous function to combine all terms past the first
                #   to extract space group name
                space_group = reduce((x, y) -> x * " " * y, line[2:end])
                space_group = split(space_group, [''', '"'], keepempty=false)[1]
                if space_group == "P1" || space_group == "P 1" ||
                        space_group == "-P1"
                    # simplify by only having one P1 space_group name
                    space_group = "P1"
                    is_p1 = true
                end
            end

            # checking for information about atom sites and symmetry
            if line[1] == "loop_"
                # creating dictionary of column names to determine what should be done
                atom_column_name = ""
                # name_to_column is a dictionary that e.g. returns which column contains x fractional coord
                #   use example: name_to_column["_atom_site_fract_x"] gives 3
                name_to_column = Dict{AbstractString, Int}()

                i += 1
                loop_starts = i
                while length(split(lines[i])) == 1 && split(lines[i])[1][1] == '_'
                    if i == loop_starts
                        atom_column_name = split(lines[i])[1]
                    end
                    name_to_column[split(lines[i])[1]] = i + 1 - loop_starts
                    # iterate to next line in file
                    i += 1
                end

                fractional = fractional || haskey(name_to_column, "_atom_site_fract_x") &&
                                haskey(name_to_column, "_atom_site_fract_y") &&
                                haskey(name_to_column, "_atom_site_fract_z")
                # if the file provides cartesian coordinates
                cartesian = cartesian || ! fractional && haskey(name_to_column, "_atom_site_Cartn_x") &&
                                haskey(name_to_column, "_atom_site_Cartn_y") &&
                                haskey(name_to_column, "_atom_site_Cartn_z")
                                             # if both are provided, will default
                                             #  to using fractional, so keep cartesian
                                             #  false

                # =====================
                # SYMMETRY READER
                # =====================
                if haskey(name_to_column, "_symmetry_equiv_pos_as_xyz")
                    symmetry_info = true

                    symmetry_count = 0
                    # CSD stores symmetry as one column in a string that ends
                    #   up getting split on the spaces between commas (i.e. its
                    #   not really one column) the length(name_to_column) + 2
                    #   should catch this hopefully there aren't other weird
                    #   ways of writing cifs...
                    while i <= length(lines) && length(lines[i]) > 0 && lines[i][1] != '_' && !occursin("loop_", lines[i])
                        symmetry_count += 1
                        line = lines[i]
                        sym_funcs = split(line, [' ', ',', ''', '"'], keepempty=false)

                        # store as strings so it can be written out later
                        new_sym_rule = Array{AbstractString, 1}(undef, 3)

                        sym_start = name_to_column["_symmetry_equiv_pos_as_xyz"] - 1
                        for j = 1:3
                            new_sym_rule[j] = sym_funcs[j + sym_start]
                        end

                        symmetry_ops = [symmetry_ops new_sym_rule]

                        i += 1
                    end

                    @assert symmetry_count == size(symmetry_ops, 2) "number of symmetry rules must match the count"

                    # finish reading in symmetry information, skip to next
                    #   iteration of outer while-loop
                    continue
                # =====================
                # FRACTIONAL READER
                # =====================
                elseif fractional && ! atom_info
                    atom_info = true
                    atom_column_name = [name for (name, column) in name_to_column if column == 1][end]
                    
                    while i <= length(lines) && length(split(lines[i])) == length(name_to_column)
                        line = split(lines[i])

                        push!(species, Symbol(line[name_to_column[atom_column_name]]))
                        coords = [coords [mod(parse(Float64, split(line[name_to_column["_atom_site_fract_x"]], '(')[1]), 1.0),
                                mod(parse(Float64, split(line[name_to_column["_atom_site_fract_y"]], '(')[1]), 1.0),
                                mod(parse(Float64, split(line[name_to_column["_atom_site_fract_z"]], '(')[1]), 1.0)]]
                        # If charges present, import them
                        if haskey(name_to_column, "_atom_site_charge")
                            push!(charge_values, parse(Float64, line[name_to_column["_atom_site_charge"]]))
                        else
                            push!(charge_values, 0.0)
                        end
                        # add to label_num_to_idx so that bonds can be converted later
                        if read_bonds_from_file
                            label_num_to_idx[line[name_to_column["_atom_site_label"]]] = length(species)
                        end
                        # iterate to next line in file
                        i += 1
                    end
                    # set up graph of correct size
                    bonds = SimpleGraph(length(species))
                    # finish reading in atom_site information, skip to next
                    #   iteration of outer while-loop
                    # prevents skipping a line after finishing reading atoms
                    continue
                # =====================
                # CARTESIAN READER
                # =====================
                elseif cartesian && ! atom_info
                    atom_info = true
                    atom_column_name = ""
                    for (name, column) in name_to_column
                        if column == 1
                            atom_column_name = name
                        end
                    end

                    while i <= length(lines) && length(split(lines[i])) == length(name_to_column)
                        line = split(lines[i])

                        push!(species, Symbol(line[name_to_column[atom_column_name]]))
                        coords = [coords [parse(Float64, split(line[name_to_column["_atom_site_Cartn_x"]], '(')[1]),
                                parse(Float64, split(line[name_to_column["_atom_site_Cartn_y"]], '(')[1]),
                                parse(Float64, split(line[name_to_column["_atom_site_Cartn_z"]], '(')[1])]]
                        # If charges present, import them
                        if haskey(name_to_column, "_atom_site_charge")
                            push!(charge_values, parse(Float64, line[name_to_column["_atom_site_charge"]]))
                        else
                            push!(charge_values, 0.0)
                        end
                        # add to label_num_to_idx so that bonds can be converted later
                        if read_bonds_from_file
                            label_num_to_idx[line[name_to_column["_atom_site_label"]]] = length(species)
                        end
                        # iterate to next line in file
                        i += 1
                    end
                    # set up graph of correct size
                    bonds = SimpleGraph(length(species))
                    # finish reading in atom_site information, skip to next
                    #   iteration of outer while-loop
                    # prevents skipping a line after finishing reading atoms
                    continue
                # =====================
                # BOND READER
                # =====================
                elseif read_bonds_from_file &&
                       haskey(name_to_column, "_geom_bond_atom_site_label_1") &&
                       haskey(name_to_column, "_geom_bond_atom_site_label_2")
                    while i <= length(lines) && length(split(lines[i])) == length(name_to_column)
                        line = split(lines[i])

                        atom_one_idx = label_num_to_idx[line[name_to_column["_geom_bond_atom_site_label_1"]]]
                        atom_two_idx = label_num_to_idx[line[name_to_column["_geom_bond_atom_site_label_2"]]]
                        add_edge!(bonds, atom_one_idx, atom_two_idx) 

                        # iterate to next line in file
                        i += 1
                    end

                    # skip to next iteration in outer while loop
                    continue
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

            i += 1
        end # End loop over lines

        if !atom_info
            error("Could not find _atom_site* after loop_ in .cif file\n")
        end

        # Structure must either be in P1 symmetry or have replication information
        if ! is_p1 && !symmetry_info
            error(@sprintf("%s is not in P1 symmetry and the .cif does not have a _symmetry_equiv_pos_as_xyz column 
            for us to apply symmetry operations to convert into P1 symmetry.", filename))
        end

        a = data["a"]
        b = data["b"]
        c = data["c"]
        α = data["alpha"]
        β = data["beta"]
        γ = data["gamma"]

        # redo coordinates if they were read in cartesian
        if cartesian && ! fractional
            coords = Box(a, b, c, α, β, γ).c_to_f * coords
        end

    # Start of cssr reader #TODO make sure this works for different .cssr files!
    ###################################
    # CSSR READER
    ###################################
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
            push!(species, Symbol(line[2]))

            push!(xf, parse(Float64, line[3]))
            push!(yf, parse(Float64, line[4]))
            push!(zf, parse(Float64, line[5]))

            push!(charge_values, parse(Float64, line[14]))
        end

        for i = 1:n_atoms
            coords = [ coords [xf[i], yf[i], zf[i]] ]
        end

        # add P1 symmetry rules for consistency
        symmetry_ops = [symmetry_ops ["x", "y", "z"]]
        is_p1 = true
        space_group = "P1"
    end

    # Construct the unit cell box
    box = Box(a, b, c, α, β, γ)
    # construct atoms attribute of crystal
    atoms = Atoms(species, Frac(coords))
    # construct charges attribute of crystal
    if ! include_zero_charges
        # include only nonzero charges
        idx_nz = charge_values .!= 0.0
        charges = Charges(charge_values[idx_nz], Frac(coords[:, idx_nz]))
    else
        # include all charges, even if some are zero.
        charges = Charges(charge_values, Frac(coords))
    end
    
    symmetry = SymmetryInfo(symmetry_ops, space_group, is_p1)
    crystal = Crystal(filename, box, atoms, charges, bonds, symmetry)

    if check_neutrality
        if ! neutral(crystal, net_charge_tol)
            error(@sprintf("Crystal %s is not charge neutral; net charge is %f e. Ignore
            this error message by passing check_charge_neutrality=false or increasing the
            net charge tolerance `net_charge_tol`\n",
                            crystal.name, total_charge(crystal)))
        end
    end

    if convert_to_p1 && ! is_p1 && ! read_bonds_from_file
        @warn @sprintf("Crystal %s has %s space group. We are converting it to P1 symmetry for use in molecular simulations.
        To afrain from this, pass `convert_to_p1=false` to the `Crystal` constructor.\n",
            crystal.name, crystal.space_group)
        return apply_symmetry_ops(crystal;
                        check_charge_neutrality=check_charge_neutrality,
                        check_atom_and_charge_overlap=check_atom_and_charge_overlap)
    end
    
    if wrap_coords
        wrap!(crystal)
    end

    if check_overlap
        overlap_flag, overlapping_pairs = overlap(crystal.atoms.coords, crystal.box, true)
        if overlap_flag
            for p in overlapping_pairs
                i, j = p
                @warn @sprintf("atom %d (%s) and %d (%s) are overlapping\n", i, crystal.atoms.species[i],
                    j, crystal.atoms.species[j])
            end
            error(crystal.name * " has overlapping pairs of atoms.
            pass `check_overlap=false` then run `overlap(crystal.atoms.coords, crystal.box, true)` 
            to obtain pairs of overlapping atoms.`\n")
        end
    end

    return crystal
end

# documented in matter.jl
function wrap!(crystal::Crystal)
    wrap!(crystal.atoms.coords)
    wrap!(crystal.charges.coords)
end

# documented in matter.jl
net_charge(crystal::Crystal) = net_charge(crystal.charges)

 # """
 #     replicated_frame = replicate(crystal, repfactors)
 # 
 # Replicates the atoms and charges in a `Crystal` in positive directions to
 # construct a new `Crystal`. Note `replicate(crystal, (1, 1, 1))` returns the same `Crystal`.
 # 
 # # Arguments
 # - `crystal::Crystal`: The crystal to replicate
 # - `repfactors::Tuple{Int, Int, Int}`: The factors by which to replicate the crystal structure in each direction.
 # 
 # # Returns
 # - `replicated_frame::Crystal`: Replicated crystal
 # """
 # function replicate(crystal::Crystal, repfactors::Tuple{Int, Int, Int})
 #     @assert ne(crystal.bonds) == 0 @sprintf("The crystal %s has bonds within it. Remove the bonds with `remove_bonds!` to replicate, and then use `infer_bonds(crystal)` to recalculate bond information", crystal.name)
 #     assert_P1_symmetry(crystal)
 #     # determine number of atoms in replicated crystal
 #     n_atoms = size(crystal.atoms.xf, 2) * repfactors[1] * repfactors[2] * repfactors[3]
 # 
 #     # replicate box
 #     new_box = replicate(crystal.box, repfactors)
 # 
 #     # replicate atoms and charges
 #     charge_coords = Array{Float64, 2}(undef, 3, 0)
 #     charge_vals = Array{Float64, 1}()
 #     atom_coords = Array{Float64, 2}(undef, 3, 0)
 #     species = Array{Symbol, 1}()
 #     for ra = 0:(repfactors[1] - 1), rb = 0:(repfactors[2] - 1), rc = 0:(repfactors[3] - 1)
 #         for i = 1:crystal.atoms.n
 #             xf = crystal.atoms.xf[:, i] + 1.0 * [ra, rb, rc]
 #             # scale fractional coords
 #             xf = xf ./ repfactors
 #             atom_coords = [atom_coords xf]
 #             push!(species, Symbol(crystal.atoms.species[i]))
 #         end
 #         for j = 1:crystal.charges.n_charges
 #             xf = crystal.charges.xf[:, j] + 1.0 * [ra, rb, rc]
 #             # scale fractional coords
 #             xf = xf ./ repfactors
 #             charge_coords = [charge_coords xf]
 #             push!(charge_vals, crystal.charges.q[j])
 #         end
 #     end
 # 
 #     new_atoms = Atoms(species, atom_coords)
 #     new_charges = Charges(charge_vals, charge_coords)
 # 
 #     @assert (new_charges.n_charges == crystal.charges.n_charges * prod(repfactors))
 #     @assert (new_atoms.n == crystal.atoms.n * prod(repfactors))
 #     return Crystal(crystal.name, new_box, new_atoms, new_charges,
 #                      symmetry=deepcopy(crystal.symmetry),
 #                      space_group=crystal.space_group, is_p1=crystal.is_p1)
 # end
 # 
 # # doc string in Misc.jl
 # function write_xyz(crystal::Crystal, filename::AbstractString;
 #                       comment::AbstractString="", center::Bool=false)
 #     atoms = crystal.atoms.species
 #     x = zeros(Float64, 3, crystal.atoms.n)
 #     for i = 1:crystal.atoms.n
 #         x[:, i] = crystal.box.f_to_c * crystal.atoms.xf[:, i]
 #     end
 #     if center
 #         center_of_box = crystal.box.f_to_c * [0.5, 0.5, 0.5]
 #         for i = 1:crystal.atoms.n
 #             x[:, i] -= center_of_box
 #         end
 #     end
 # 
 #     write_xyz(atoms, x, filename, comment=comment)
 # end
 # write_xyz(crystal::Crystal; comment::AbstractString="", center::Bool=false) = write_xyz(
 #     crystal,
 #     replace(replace(crystal.name, ".cif" => ""), ".cssr" => "") * ".xyz",
 #     comment=comment, center=center)
 # 
"""
    overlap = overlap(crystal; tol=0.1, verbose=false)

To detect duplicate atoms and duplicate charges, go through the atoms and charges and make sure no atoms
overlap with each other and no charges overlap with each other. here, overlap is defined when the distance
is less than `tol`.

# Arguments
- `crystal::Crystal`: the crystal
- `tol::Float64`: if distance of particles is less than `tol`, they are classified as overlapping.

# Returns
- `overlap::Bool`: `true` if *any* particles in the crystal overlap.
"""
function overlap(crystal::Crystal; tol::Float64=0.1, verbose::Bool=true)
    overlap = false
    for i = 1:crystal.atoms.n
        for j = 1:crystal.atoms.n
            if j >= i
                continue
            end
            if _overlap(crystal.atoms.xf[:, i], crystal.atoms.xf[:, j],
                        crystal.box, overlap_tol)
                overlap = true
                if verbose
                    @warn @sprintf("Atoms %d and %d in %s are less than %f Å apart.", i, j,
                        crystal.name, overlap_tol)
                end
            end
        end
    end
    return overlap
end

 # function charge_overlap(crystal::Crystal; overlap_tol::Float64=0.1, verbose::Bool=true)
 #     overlap = false
 #     for i = 1:crystal.charges.n_charges
 #         for j = 1:crystal.charges.n_charges
 #             if j >= i
 #                 continue
 #             end
 #             if _overlap(crystal.charges.xf[:, i], crystal.charges.xf[:, j],
 #                         crystal.box, overlap_tol)
 #                 overlap = true
 #                 if verbose
 #                     @warn @sprintf("Charges %d and %d in %s are less than %f Å apart.", i, j,
 #                         crystal.name, overlap_tol)
 #                 end
 #             end
 #         end
 #     end
 #     return overlap
 # end
 # 
 # # determine if two atoms overlap, returns the number of Atoms that
 # #   do overlap, and can then use that number to determine if they overlap or are repeats
 # function _overlap(xf_1::Array{Float64, 1}, xf_2::Array{Float64, 1},
 #                   box::Box, overlap_tol::Float64)
 #     dxf = mod.(xf_1, 1.0) .- mod.(xf_2, 1.0)
 #     nearest_image!(dxf)
 #     dxc = box.f_to_c * dxf
 #     return norm(dxc) < overlap_tol
 # end
 # 
 # function _overlap(xf::Union{Charges, Atoms}, box::Box, overlap_tol::Float64)
 #     return _overlap(xf, xf, box, overlap_tol)
 # end

# docstring in matter.jl
neutral(crystal::Crystal, tol::Float64=1e-5) = neutral(crystal.charges, tol)

 # 
 # """
 #     charged_flag = charged(crystal, verbose=false) # true or false
 # 
 # Determine if a crystal has point charges
 # """
 # function charged(crystal::Crystal; verbose::Bool=false)
 #     charged_flag = crystal.charges.n_charges > 0
 #     if verbose
 #         @printf("\tCrystal atoms of %s have charges? %s\n", crystal.name, charged_flag)
 #     end
 #     return charged_flag
 # end
 # 
"""
    strip_numbers_from_atom_labels!(crystal)

Strip numbers from labels for `crystal.atoms`.
Precisely, for `atom` in `crystal.atoms`, find the first number that appears in `atom`.
Remove this number and all following characters from `atom`.
e.g. C12 --> C
	 Ba12A_3 --> Ba

# Arguments
- `crystal::Crystal`: The crystal containing the crystal structure information
"""
function strip_numbers_from_atom_labels!(crystal::Crystal)
    for i = 1:crystal.atoms.n
        # atom species in string format
		species = string(crystal.atoms.species[i])
		for j = 1:length(species)
			if ! isletter(species[j])
                crystal.atoms.species[i] = Symbol(species[1:j-1])
				break
			end
		end
	end
end

 # write_vtk(crystal::Crystal) = write_vtk(crystal.box, split(crystal.name, ".")[1])
 # 
 # """
 #     formula = chemical_formula(crystal, verbose=false)
 # 
 # Find the irreducible chemical formula of a crystal structure.
 # 
 # # Arguments
 # - `crystal::Crystal`: The crystal containing the crystal structure information
 # - `verbose::Bool`: If `true`, will print the chemical formula as well
 # 
 # # Returns
 # - `formula::Dict{Symbol, Int}`: A dictionary with the irreducible chemical formula of a crystal structure
 # """
 # function chemical_formula(crystal::Crystal; verbose::Bool=false)
 #     unique_atoms = unique(crystal.atoms.species)
 #     # use dictionary to count atom types
 #     atom_counts = Dict{Symbol, Int}([a => 0 for a in unique_atoms])
 #     for i = 1:crystal.atoms.n
 #         atom_counts[crystal.atoms.species[i]] += 1
 #     end
 # 
 #     # get greatest common divisor
 #     gcd_ = gcd([k for k in values(atom_counts)])
 # 
 #     # turn into irreducible chemical formula
 #     for atom in keys(atom_counts)
 #         atom_counts[atom] = atom_counts[atom] / gcd_
 #     end
 # 
 #     # print result
 #     if verbose
 #         @printf("Chemical formula of %s:\n\t", crystal.name)
 #         for (atom, formula_unit) in atom_counts
 # 			@printf("%s_%d", string(atom), formula_unit)
 #         end
 #         @printf("\n")
 #     end
 # 
 #     return atom_counts
 # end
 # 
 # """
 # 
 #     mass_of_crystal = molecular_weight(crystal)
 # 
 # Calculates the molecular weight of a unit cell of the crystal in amu using information stored in `data/atomicmasses.csv`.
 # 
 # # Arguments
 # - `crystal::Crystal`: The crystal containing the crystal structure information
 # 
 # # Returns
 # - `mass_of_crystal::Float64`: The molecular weight of a unit cell of the crystal in amu
 # """
 # function molecular_weight(crystal::Crystal)
 #     atomic_masses = read_atomic_masses()
 # 
 #     mass = 0.0
 # 	for i = 1:crystal.atoms.n
 #         mass += atomic_masses[crystal.atoms.species[i]]
 #     end
 # 
 #     return mass # amu
 # end
 # 
 # """
 #     ρ = crystal_density(crystal) # kg/m²
 # 
 # Compute the crystal density of a crystal. Pulls atomic masses from [`read_atomic_masses`](@ref).
 # 
 # # Arguments
 # - `crystal::Crystal`: The crystal containing the crystal structure information
 # 
 # # Returns
 # - `ρ::Float64`: The crystal density of a crystal in kg/m³
 # """
 # function crystal_density(crystal::Crystal)
 #     mw = molecular_weight(crystal)
 #     return mw / crystal.box.Ω * 1660.53892  # --> kg/m3
 # end
 # 
 # """
 #     simulation_ready_crystal = apply_symmetry_ops(non_p1_crystal;
 #                                                 check_charge_neutrality=true,
 #                                                 net_charge_tol=0.001,
 #                                                 check_atom_and_charge_overlap=true,
 #                                                 remove_overlap=false,
 #                                                 wrap_to_unit_cell=true)
 # 
 # Convert a crystal to P1 symmetry based on internal symmetry rules. This will
 # return the new crystal.
 # 
 # # Arguments
 # - `f::Crystal`: The crystal to be converted to P1 symmetry
 # - `check_charge_neutrality::Bool`: check for charge neutrality
 # - `net_charge_tol::Float64`: when checking for charge neutrality, throw an error if the absolute value of the net charge is larger than this value.
 # - `check_atom_and_charge_overlap::Bool`: throw an error if overlapping atoms are detected.
 # - `remove_overlap::Bool`: remove identical atoms automatically. Identical atoms are the same element atoms which overlap.
 # - `wrap_to_unit_cell::Bool`: if true, enforce that fractional coords of atoms/charges are in [0,1]³ by mod(x, 1)
 # 
 # # Returns
 # - `P1_crystal::Crystal`: The crystal after it has been converted to P1
 #     symmetry. The new symmetry rules will be the P1 symmetry rules
 # """
 # function apply_symmetry_ops(crystal::Crystal; check_charge_neutrality::Bool=true,
 #                               net_charge_tol::Float64=0.001, check_atom_and_charge_overlap::Bool=true,
 #                               remove_overlap::Bool=false, wrap_to_unit_cell::Bool=true)
 #     if crystal.is_p1
 #         return crystal
 #     end
 #     new_atom_xfs = Array{Float64, 2}(undef, 3, crystal.atoms.n * size(crystal.symmetry, 2))
 #     new_charge_xfs = Array{Float64, 2}(undef, 3, crystal.charges.n_charges * size(crystal.symmetry, 2))
 #     new_atom_species = Array{Symbol, 1}(undef, 0)
 #     new_charge_qs = Array{Float64, 1}(undef, 0)
 # 
 #     # for each symmetry rule
 #     for i in 1:size(crystal.symmetry, 2)
 #         # loop over all atoms in lower level symmetry
 #         sym_rule = eval.(Meta.parse.("(x, y, z) -> " .* crystal.symmetry[:, i]))
 #         for j in 1:size(crystal.atoms.xf, 2)
 #             # apply current symmetry rule to current atom for x, y, and z coordinates
 #             current_atom_idx = (i - 1) * crystal.atoms.n + j
 #             new_atom_xfs[:, current_atom_idx] .= [Base.invokelatest.(
 #                         sym_rule[k], crystal.atoms.xf[:, j]...) for k in 1:3]
 #         end
 #         # loop over all charges in lower level symmetry
 #         for j in 1:size(crystal.charges.xf, 2)
 #             # apply current symmetry rule to current atom for x, y, and z coordinates
 #             current_charge_idx = (i - 1) * crystal.charges.n_charges + j
 #             new_charge_xfs[:, current_charge_idx] .= [Base.invokelatest.(
 #                         sym_rule[k], crystal.charges.xf[:, j]...) for k in 1:3]
 #         end
 #         # repeat charge_qs and atom_species for every symmetry applied
 #         new_atom_species = [new_atom_species; crystal.atoms.species]
 #         new_charge_qs = [new_charge_qs; crystal.charges.q]
 #     end
 # 
 #     new_crystal = Crystal(crystal.name, crystal.box,
 #         Atoms(new_atom_species, new_atom_xfs),
 #         Charges(new_charge_qs, new_charge_xfs))
 # 
 #     if check_charge_neutrality
 #         if ! charge_neutral(new_crystal, net_charge_tol)
 #             error(@sprintf("Crystal %s is not charge neutral; net charge is %f e. Ignore
 #             this error message by passing check_charge_neutrality=false or increasing the
 #             net charge tolerance `net_charge_tol`\n",
 #                             new_crystal.name, total_charge(new_crystal)))
 #         end
 #     end
 #         
 #     if wrap_to_unit_cell
 #         wrap_atoms_to_unit_cell!(crystal)
 #     end
 # 
 #     if remove_overlap
 #         return remove_overlapping_atoms_and_charges(new_crystal)
 #     end
 # 
 #     if check_atom_and_charge_overlap
 #         if atom_overlap(new_crystal) | charge_overlap(new_crystal)
 #             error(@sprintf("At least one pair of atoms/charges overlap in %s.
 #             Consider passing `remove_overlap=true`\n", new_crystal.name))
 #         end
 #     end
 # 
 #     return new_crystal
 # end
 # 
 # """
 #     symmetry_equal = is_symmetry_equal(crystal1.symmetry, crystal2.symmetry)
 # 
 # Returns true if both symmetry rules can create the same set from the same set
 # of coordinates. Returns false if they don't contain the same number of rules or
 # if they create different sets of points.
 # 
 # # Arguments
 # - `sym1::Array{AbstractString, 2}`: Array of strings that represent
 #     symmetry operations
 # - `sym2::Array{AbstractString, 2}`: Array of strings that represent
 #     symmetry operations
 # 
 # # Returns
 # - `is_equal::Bool`: True if they are the same set of symmetry rules
 #     False if they are different
 # """
 # function is_symmetry_equal(sym1::Array{AbstractString, 2}, sym2::Array{AbstractString, 2})
 #     # need same number of symmetry operations
 #     if size(sym1, 2) != size(sym2, 2)
 #         return false
 #     end
 #     # define a test array that operations will be performed on
 #     test_array = [0.0 0.25 0.0  0.0  0.0  0.25 0.25 0.25;
 #                   0.0 0.0  0.25 0.0  0.25 0.0  0.25 0.25;
 #                   0.0 0.0  0.0  0.25 0.25 0.25 0.25 0.25]
 #     # set up both arrays for storing replicated coords
 #     sym1_applied_to_test = Array{Float64, 2}(undef, 3, 0)
 #     sym2_applied_to_test = Array{Float64, 2}(undef, 3, 0)
 # 
 #     # loop over all positions in the test_array
 #     for i in 1:size(test_array, 2)
 #         # loop over f1 symmetry rules
 #         for j in 1:size(sym1, 2)
 #             sym1_applied_to_test = [sym1_applied_to_test [Base.invokelatest.(
 #                 eval(Meta.parse("(x, y, z) -> " * sym1[k, j])), test_array[:, i]...) for k in 1:3]]
 #         end
 #         # loop over f2 symmetry rules
 #         for j in 1:size(sym2, 2)
 #             sym2_applied_to_test = [sym2_applied_to_test [Base.invokelatest.(
 #                 eval(Meta.parse("(x, y, z) -> " * sym2[k, j])), test_array[:, i]...) for k in 1:3]]
 #         end
 #     end
 # 
 #     # convert to sets for using issetequal, symmetry rules might be in a a different order
 #     sym1_set = Set([sym1_applied_to_test[:, i] for i in 1:size(sym1_applied_to_test, 2)])
 #     sym2_set = Set([sym2_applied_to_test[:, i] for i in 1:size(sym2_applied_to_test, 2)])
 # 
 #     # return if the sets of coords are equal
 #     return issetequal(sym1_set, sym2_set)
 # end
 # 
 # """
 #     assert_P1_symmetry(crystal::Crystal)
 # 
 # Throw an error if and only if the crystal is not in P1 symmetry.
 # """
 # function assert_P1_symmetry(crystal::Crystal)
 #     if ! crystal.symmetry.is_p1 
 #         error("The crystal %s is not in P1 symmetry.\n
 #                To convert to P1 symmetry, try:\n
 #                \tcrystal_p1 = apply_symmetry_ops(crystal)", crystal.name)
 # end
 # 
 # """
 #     write_cif(crystal, filename; fractional=true)
 # 
 # Write a `crystal::Crystal` to a .cif file with `filename::AbstractString`. If `filename` does
 # not include the .cif extension, it will automatically be added.
 # """
 # function write_cif(crystal::Crystal, filename::AbstractString; fractional::Bool=true)
 #     if charged(crystal) && (crystal.atoms.n != crystal.charges.n_charges)
 #         error("write_cif assumes equal numbers of Charges and Atoms (or zero charges)")
 #     end
 # 
 #     # create dictionary for tracking label numbers
 #     label_numbers = Dict{Symbol, Int}()
 #     for atom in crystal.atoms.species
 #         if !haskey(label_numbers, atom)
 #             label_numbers[atom] = 1
 #         end
 #     end
 # 
 #     # append ".cif" to filename if it doesn't already have the extension
 #     if ! occursin(".cif", filename)
 #         filename *= ".cif"
 #     end
 #     cif_file = open(filename, "w")
 #     # first line should be data_xtalname_PM
 #     if crystal.name == ""
 #         @printf(cif_file, "data_PM\n")
 #     else
 #         # don't include file extension!
 #         @printf(cif_file, "data_%s_PM\n", split(crystal.name, ".")[1])
 #     end
 # 
 #     @printf(cif_file, "_symmetry_space_group_name_H-M\t'%s'\n", crystal.space_group)
 # 
 #     @printf(cif_file, "_cell_length_a\t%f\n", crystal.box.a)
 #     @printf(cif_file, "_cell_length_b\t%f\n", crystal.box.b)
 #     @printf(cif_file, "_cell_length_c\t%f\n", crystal.box.c)
 # 
 #     @printf(cif_file, "_cell_angle_alpha\t%f\n", crystal.box.α * 180.0 / pi)
 #     @printf(cif_file, "_cell_angle_beta\t%f\n", crystal.box.β * 180.0 / pi)
 #     @printf(cif_file, "_cell_angle_gamma\t%f\n", crystal.box.γ * 180.0 / pi)
 # 
 #     @printf(cif_file, "_symmetry_Int_Tables_number 1\n\n")
 #     @printf(cif_file, "loop_\n_symmetry_equiv_pos_as_xyz\n")
 #     for i in 1:size(crystal.symmetry, 2)
 #         @printf(cif_file, "'%s,%s,%s'\n", crystal.symmetry[:, i]...)
 #     end
 #     @printf(cif_file, "\n")
 # 
 #     @printf(cif_file, "loop_\n_atom_site_label\n_atom_site_type_symbol\n")
 #     if fractional
 #         @printf(cif_file, "_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n")
 #     else
 #         @printf(cif_file, "_atom_site_Cartn_x\n_atom_site_Cartn_y\n_atom_site_Cartn_z\n")
 #     end
 #     @printf(cif_file, "_atom_site_charge\n")
 # 
 #     idx_to_label = Array{AbstractString, 1}(undef, crystal.atoms.n)
 #     for i = 1:crystal.atoms.n
 #         q = 0.0
 #         if charged(crystal)
 #             q = crystal.charges.q[i]
 #             if ! isapprox(crystal.charges.xf[:, i], crystal.atoms.xf[:, i])
 #                 error("write_cif assumes charges correspond to LJspheres")
 #             end
 #         end
 #         # print label and type symbol
 #         @printf(cif_file, "%s\t%s\t", string(crystal.atoms.species[i]) *
 #                 string(label_numbers[crystal.atoms.species[i]]),
 #                 crystal.atoms.species[i])
 #         # store label for this atom idx
 #         idx_to_label[i] = string(crystal.atoms.species[i]) *
 #                     string(label_numbers[crystal.atoms.species[i]])
 #         # increment label
 #         label_numbers[crystal.atoms.species[i]] += 1
 #         if fractional
 #             @printf(cif_file, "%f\t%f\t%f\t%f\n", crystal.atoms.xf[:, i]..., q)
 #         else
 #             
 #             @printf(cif_file, "%f\t%f\t%f\t%f\n", (crystal.box.f_to_c * crystal.atoms.xf[:, i])..., q)
 #         end
 #     end
 # 
 #     # only print bond information if it is in the crystal
 #     if ne(crystal.bonds) > 0
 #         # print column names for bond information
 #         @printf(cif_file, "\nloop_\n_geom_bond_atom_site_label_1\n_geom_bond_atom_site_label_2\n_geom_bond_distance\n")
 # 
 #         for edge in collect(edges(crystal.bonds))
 #             dxf = crystal.atoms.xf[:, edge.src] - crystal.atoms.xf[:, edge.dst]
 #             nearest_image!(dxf)
 #             @printf(cif_file, "%s\t%s\t%0.5f\n", idx_to_label[edge.src], idx_to_label[edge.dst],
 #                     norm(dxf))
 #         end
 #     end
 #     close(cif_file)
 # end
 # 
 # """
 #     write_bond_information(crystal, filename)
 #     write_bond_information(crystal)
 # 
 # Writes the bond information from a crystal to the selected filename.
 # 
 # # Arguments
 # -`crystal::Crystal`: The crystal to have its bonds written to a vtk file
 # -`filename::AbstractString`: The filename the bond information will be saved to. If left out, will default to crystal name.
 # """
 # function write_bond_information(crystal::Crystal, filename::AbstractString)
 #     if ne(crystal.bonds) == 0
 #         @warn("Crystal %s has no bonds present. To get bonding information for this crystal run `infer_bonds!` with an array of bonding rules\n", crystal.name)
 #     end
 #     if ! occursin(".vtk", filename)
 #         filename *= ".vtk"
 #     end
 # 
 #     vtk_file = open(filename, "w")
 # 
 #     @printf(vtk_file, "# vtk DataFile Version 2.0\n%s bond information\nASCII\nDATASET POLYDATA\nPOINTS %d double\n", crystal.name, nv(crystal.bonds))
 # 
 #     for i = 1:crystal.atoms.n
 #         @printf(vtk_file, "%0.5f\t%0.5f\t%0.5f\n", (crystal.box.f_to_c * crystal.atoms.xf[:, i])...)
 #     end
 #     @printf(vtk_file, "\nLINES %d %d\n", ne(crystal.bonds), 3 * ne(crystal.bonds))
 #     for edge in collect(edges(crystal.bonds))
 #         @printf(vtk_file, "2\t%d\t%d\n", edge.src - 1, edge.dst - 1)
 #     end
 #     close(vtk_file)
 #     @printf("Saving bond information for crystal %s to %s.\n", crystal.name, joinpath(pwd(), filename))
 # end
 # 
 # write_bond_information(crystal::Crystal) = write_bond_information(crystal, crystal.name * "_bonds.vtk")
 # 
"""
    crystal_with_charges = assign_charges(crystal, species_to_charge, net_charge_tol=1e-5)

assign charges to the atoms present in the crystal based on atom type.
pass a dictionary `species_to_charge` that maps atomic species to a charge.

if the crystal already has charges, the charges are removed and new charges are added. a warning is thrown if this is the case.

checks for charge neutrality in the end.

returns a new crystal.

# Examples

```
species_to_charge = Dict(:Ca => 2.0, :C => 1.0, :H => -1.0)
crystal_with_charges = assign_charges(crystal, species_to_charge, 1e-7)
crystal_with_charges = assign_charges(crystal, species_to_charge) # tol 1e-5 default
```

# Arguments
- `crystal::Crystal`: the crystal 
- `species_to_charge::Dict{Symbol, Float64}`: a dictionary that maps atomic species to charge
- `net_charge_tol::Float64`: the net charge tolerated when asserting charge neutrality of
the resulting crystal
"""
function assign_charges(crystal::Crystal, species_to_charge::Dict{Symbol, Float64}, net_charge_tol::Float64=1e-5)
    if crystal.charges.n != 0
        @warn @sprintf("Charges are already present in %s. Removing the current charges on the
        crystal and adding new ones...\n", crystal.name)
    end

    q = [species_to_charge[s] for s in crystal.atoms.species]
    charges = Charges(q, crystal.atoms.coords)
    
    new_crystal = Crystal(crystal.name, crystal.box, crystal.atoms, charges, crystal.bonds, crystal.symmetry)

    # check for charge neutrality
    if ! neutral(new_crystal, net_charge_tol)
        error(@sprintf("Net charge of crystal %s = %f > net charge tolerance %f. If
        charge neutrality is not a problem, pass `net_charge_tol=Inf`\n", crystal.name,
        net_charge(new_crystal), net_charge_tol))
    end

    return new_crystal
end

 # """
 #     r = distance(crystal, i, j, apply_pbc)
 # 
 # Calculate the (Cartesian) distance between atoms `i` and `j` in a crystal.
 # `if apply_pbc`, use the nearest image convention to apply periodic boundary conditions.
 # `if ! apply_pbc`, do not apply periodic boundary conditions
 # 
 # # Arguments
 # -`crystal::Crystal`: the crystal structure
 # -`i::Int`: index of the first atom
 # -`j::Int`: Index of the second atom
 # - `apply_pbc::Bool`: `true` if we wish to apply periodic boundary conditions, `false` otherwise
 # """
 # function distance(crystal::Crystal, i::Int, j::Int, apply_pbc::Bool)
 #     dxf = crystal.atoms.xf[:, i] - crystal.atoms.xf[:, j]
 #     if apply_pbc
 #         nearest_image!(dxf)
 #     end
 #     return norm(crystal.box.f_to_c * dxf)
 # end
 # 
 # function Base.show(io::IO, crystal::Crystal)
 #     println(io, "Name: ", crystal.name)
 #     println(io, crystal.box)
 # 	@printf(io, "Number of atoms = %d\n", crystal.atoms.n)
 # 	@printf(io, "Number of charges = %d\n", crystal.charges.n)
 #     println(io, "Chemical formula: ", chemical_formula(crystal))
 #     @printf(io, "Space Group: %s\n", crystal.space_group)
 #     @printf(io, "Symmetry Operations:\n")
 #     for i in 1:size(crystal.symmetry, 2)
 #         @printf(io, "\t'%s, %s, %s'\n", crystal.symmetry[:, i]...)
 #     end
 # end
 # 
 # function has_same_sets_of_atoms_and_charges(f1::Crystal, f2::Crystal; atol::Float64=1e-6, checknames::Bool=false)
 #     names_flag = f1.name == f2.name
 #     if checknames && (! names_flag)
 #         return false
 #     end
 #     box_flag = isapprox(f1.box, f2.box)
 #     if f1.charges.n_charges != f2.charges.n_charges
 #         return false
 #     end
 #     if f1.atoms.n != f2.atoms.n
 #         return false
 #     end
 #     charges_flag = has_same_set_of_charges(f1.charges, f2.charges; atol=atol)
 #     atoms_flag = has_same_set_of_atoms(f1.atoms, f2.atoms; atol=atol)
 #     symmetry_flag = is_symmetry_equal(f1.symmetry, f2.symmetry)
 #     return box_flag && charges_flag && atoms_flag && symmetry_flag
 # end
 # 
 # 
 # function Base.isapprox(f1::Crystal, f2::Crystal)
 #     box_flag = isapprox(f1.box, f2.box)
 #     if f1.charges.n_charges != f2.charges.n_charges
 #         return false
 #     end
 #     if f1.atoms.n != f2.atoms.n
 #         return false
 #     end
 #     charges_flag = isapprox(f1.charges, f2.charges)
 #     atoms_flag = isapprox(f1.atoms, f2.atoms)
 #     symmetry_flag = is_symmetry_equal(f1.symmetry, f2.symmetry)
 #     return box_flag && charges_flag && atoms_flag && symmetry_flag
 # end
 # 
 # function Base.:+(crystals::Crystal...; check_overlap::Bool=true)
 #     new_crystal = deepcopy(crystals[1])
 #     for (i, f) in enumerate(crystals)
 #         if i == 1
 #             continue
 #         end
 #         @assert isapprox(new_crystal.box, f.box) @sprintf("Crystal %s has a different box\n", f.name)
 #         @assert is_symmetry_equal(new_crystal.symmetry, f.symmetry) @sprintf("Crystal %s has different symmetry rules\n", f.name)
 #         @assert new_crystal.space_group == f.space_group
 # 
 #         new_atoms = new_crystal.atoms + f.atoms
 #         new_charges = new_crystal.charges + f.charges
 # 
 #         nf_n_atoms = new_crystal.atoms.n
 #         add_vertices!(new_crystal.bonds, nf_n_atoms)
 #         for edge in collect(edges(f.bonds))
 #             add_edge!(new_crystal.bonds, nf_n_atoms + edge.src, nf_n_atoms + edge.dst)
 #         end
 # 
 #         new_crystal = Crystal(split(new_crystal.name, ".")[1] * "_" * split(f.name, ".")[1],
 #                                  new_crystal.box, new_atoms, new_charges,
 #                                  symmetry=new_crystal.symmetry,space_group=new_crystal.space_group,
 #                                  is_p1=new_crystal.is_p1, bonds=new_crystal.bonds)
 #     end
 #     if check_overlap
 #         if atom_overlap(new_crystal)
 #             @warn "This new crystal has overlapping atoms, use:\n`remove_overlapping_atoms_and_charges(crystal)`\nto remove them"
 #         end
 # 
 #         if charge_overlap(new_crystal)
 #             @warn "This new crystal has overlapping charges, use:\n`remove_overlapping_atoms_and_charges(crystal)`\nto remove them"
 #         end
 #     end
 # 
 #     return new_crystal
 # end
