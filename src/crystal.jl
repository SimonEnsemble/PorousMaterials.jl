"""
    SymmetryInfo(symmetry, space_group, is_p1)

# Attributes
- `operations::Array{Function, 2}`: 2D array of anonymous functions that represent
    the symmetry operations. If the structure is in P1 there will be one
    symmetry operation.
- `space_group::AbstractString`: The name of the space group. This is stored
    so that it can be written out again in the write_cif function. The space
    group is not used to verify the symmetry rules.
- `is_p1::Bool`: Stores whether the crystal is currently in P1 symmetry.
"""
struct SymmetryInfo
    operations::Array{String, 2}
    space_group::String
    is_p1::Bool
end
SymmetryInfo() = SymmetryInfo([Array{String, 2}(undef, 3, 0) ["x", "y", "z"]], "P1", true) # default

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

const NET_CHARGE_TOL = 1e-4 # net charge tolerance

"""
    crystal = Crystal(filename;
        check_neutrality=true, net_charge_tol=1e-4,
        check_overlap=true, overlap_tol=0.1,
        convert_to_p1=true, read_bonds_from_file=false, wrap_coords=true,
        include_zero_charges=false,
        remove_duplicates=false,
        species_col=["_atom_site_label", "_atom_site_type_symbol"]
        ) # read from file

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
    If `true,` include the atoms in `crystal.charges` that have zero charge, ensuring that the number of atoms is equal to the number of charges and that `crystal.charges.coords.xf` and `crystal.atoms.coords.xf` are the same.
- `remove_duplicates::Bool`: remove duplicate atoms and charges. an atom is duplicate only if it is the same species.

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
- `species_col::Array{String}`: which column to use for species identification for `crystal.atoms.species`. we use a priority list:
    we check for the first entry of `species_col` in the .cif file; if not present, we then use the second entry, and so on.
"""
function Crystal(filename::String;
                 check_neutrality::Bool=true, 
                 net_charge_tol::Float64=NET_CHARGE_TOL,
                 check_overlap::Bool=true,
                 overlap_tol::Float64=0.1,
                 convert_to_p1::Bool=true,
                 read_bonds_from_file::Bool=false, 
                 wrap_coords::Bool=true,
                 include_zero_charges::Bool=false,
                 remove_duplicates::Bool=false,
                 species_col::Array{String, 1}=["_atom_site_label", "_atom_site_type_symbol"])
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
    operations = Array{String, 2}(undef, 3, 0)
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
        # find later.
        species_column = "null"

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
                #atom_column_name = ""
                # name_to_column is a dictionary that e.g. returns which column contains x fractional coord
                #   use example: name_to_column["_atom_site_fract_x"] gives 3
                name_to_column = Dict{AbstractString, Int}()

                i += 1
                loop_starts = i
                while length(split(lines[i])) == 1 && split(lines[i])[1][1] == '_'
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

                # Assign species_column by matching to priority list
                if haskey(name_to_column, "_atom_site_Cartn_x") || haskey(name_to_column, "_atom_site_fract_x") # to have entered _atom_site loop
                    found_species_col = false
                    for col in species_col
                        if col ∈ keys(name_to_column)
                            species_column = col
                            found_species_col = true
                            break
                        end
                    end
                    if ! found_species_col
                        @error "could not find species_col=$(species_col) columns in the .cif file for atomic species labels for $(filename)."
                    end
                end

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
                        line = lines[i]
                        sym_funcs = split(line, [' ', ',', ''', '"'], keepempty=false)
                        if length(sym_funcs) == 0
                            break
                        end
                        symmetry_count += 1

                        # store as strings so it can be written out later
                        new_sym_rule = Array{AbstractString, 1}(undef, 3)

                        sym_start = name_to_column["_symmetry_equiv_pos_as_xyz"] - 1
                        for j = 1:3
                            new_sym_rule[j] = sym_funcs[j + sym_start]
                        end

                        operations = [operations new_sym_rule]

                        i += 1
                    end

                    @assert symmetry_count == size(operations, 2) "number of symmetry rules must match the count"

                    # finish reading in symmetry information, skip to next
                    #   iteration of outer while-loop
                    continue
                # =====================
                # FRACTIONAL READER
                # =====================
                elseif fractional && ! atom_info
                    atom_info = true

                    while i <= length(lines) && length(split(lines[i])) == length(name_to_column)
                        line = split(lines[i])

                        push!(species, Symbol(line[name_to_column[species_column]]))
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
                            label_num_to_idx[line[name_to_column[species_column]]] = length(species)
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

                    while i <= length(lines) && length(split(lines[i])) == length(name_to_column)
                        line = split(lines[i])

                        push!(species, Symbol(line[name_to_column[species_column]]))
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
                            label_num_to_idx[line[name_to_column[species_column]]] = length(species)
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
        bonds = SimpleGraph(n_atoms)

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
        operations = [operations ["x", "y", "z"]]
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

    symmetry = SymmetryInfo(operations, space_group, is_p1)
    crystal = Crystal(filename, box, atoms, charges, bonds, symmetry)

    if convert_to_p1 && ! is_p1 && ! read_bonds_from_file
        @info @sprintf("Crystal %s has %s space group. I am converting it to P1 symmetry.
        To afrain from this, pass `convert_to_p1=false` to the `Crystal` constructor.\n",
            crystal.name, crystal.symmetry.space_group)
        crystal = apply_symmetry_operations(crystal)
    end
    
    if wrap_coords
        wrap!(crystal) # do before checking overlap!
    end
    
    if remove_duplicates
        crystal = remove_duplicate_atoms_and_charges(crystal)
    end

    if check_neutrality
        if ! neutral(crystal, net_charge_tol)
            error(@sprintf("Crystal %s is not charge neutral; net charge is %f e. Ignore
            this error message by passing check_charge_neutrality=false or increasing the
            net charge tolerance `net_charge_tol`\n", crystal.name, net_charge(crystal)))
        end
    end

    if check_overlap
        _check_overlap(crystal)
    end

    return crystal
end

overlap(crystal::Crystal) = overlap(crystal.atoms.coords, crystal.box, true)

function _check_overlap(crystal::Crystal)
    overlap_flag, overlapping_pairs = overlap(crystal)
    if overlap_flag
        for (i, j) in overlapping_pairs
            @warn @sprintf("atom %d (%s) and %d (%s) are overlapping\n", i, crystal.atoms.species[i],
                j, crystal.atoms.species[j])
        end
        error(crystal.name * " has overlapping pairs of atoms! 
            (this occurs often when applying symmetry rules, if your structure was not in P1 symmetry to begin with.)
            pass `check_overlap=false` then run `overlap(crystal)` to find pairs of overlapping atoms for inspection. 
            or pass `remove_duplicates=true` to the `Crystal` constructor to automatically remove the duplicate atoms and charges.
            you should then visualize your .cif to make sure it is proper.`\n")
    end
end

# documented in matter.jl
function wrap!(crystal::Crystal)
    wrap!(crystal.atoms.coords)
    wrap!(crystal.charges.coords)
end

# documented in matter.jl
net_charge(crystal::Crystal) = net_charge(crystal.charges)

"""
    replicated_crystal = replicate(crystal, repfactors)

replicate the atoms and charges in a `Crystal` in positive directions to construct a new
`Crystal`. Note `replicate(crystal, (1, 1, 1))` returns the same `Crystal`. the fractional
coordinates will be rescaled to be in [0, 1].

# arguments
- `crystal::Crystal`: The crystal to replicate
- `repfactors::Tuple{Int, Int, Int}`: The factors by which to replicate the crystal structure in each crystallographic direction (a, b, c).

# returns
- `replicated_frame::Crystal`: replicated crystal
"""
function replicate(crystal::Crystal, repfactors::Tuple{Int, Int, Int})
    if ne(crystal.bonds) != 0
        error("the crystal " * crystal.name * " has assigned bonds. to replicate, remove
        its bonds with `remove_bonds!(crystal)`. then use `infer_bonds(crystal)` to
        reassign the bonds")
    end

    assert_P1_symmetry(crystal)

    n_atoms = crystal.atoms.n * prod(repfactors)
    n_charges = crystal.charges.n * prod(repfactors)

    box = replicate(crystal.box, repfactors)
    atoms = Atoms{Frac}(n_atoms)
    charges = Charges{Frac}(n_charges)

    atom_counter = 0
    charge_counter = 0
    for ra = 0:(repfactors[1] - 1), rb = 0:(repfactors[2] - 1), rc = 0:(repfactors[3] - 1)
        xf_shift = 1.0 * [ra, rb, rc]

        # replicate atoms
        for i = 1:crystal.atoms.n
            atom_counter += 1

            atoms.species[atom_counter] = crystal.atoms.species[i]

            xf = crystal.atoms.coords.xf[:, i] + xf_shift
            atoms.coords.xf[:, atom_counter] = xf ./ repfactors
        end

        # replicate charges
        for i = 1:crystal.charges.n
            charge_counter += 1

            charges.q[charge_counter] = crystal.charges.q[i]

            xf = crystal.charges.coords.xf[:, i] + xf_shift
            charges.coords.xf[:, charge_counter] = xf ./ repfactors
        end
    end

    return Crystal(crystal.name, box, atoms, charges, SimpleGraph(n_atoms), crystal.symmetry)
end

# doc string in Misc.jl
xyz_filename(crystal::Crystal) = replace(replace(crystal.name, ".cif" => ""), ".cssr" => "") * ".xyz"
function write_xyz(crystal::Crystal; comment::AbstractString="", center_at_origin::Bool=false)
    filename = xyz_filename(crystal)
    atoms = Atoms(crystal.atoms.species,
                  Cart(crystal.atoms.coords, crystal.box)
                  ) # put in Cartesian
    if center_at_origin
        x_c = crystal.box.f_to_c * [0.5, 0.5, 0.5]
        atoms.coords.x .-= x_c
        write_xyz(atoms, filename, comment=comment)
    else
        write_xyz(atoms, filename, comment=comment)
    end
end

# docstring in matter.jl
neutral(crystal::Crystal, tol::Float64=1e-5) = neutral(crystal.charges, tol)


"""
    charged = has_charges(crystal) # true or false
    charged = has_charges(molecule) # true or false

`true` if any only if `crystal::Crystal`/`molecule::Molecule` has point charges.
"""
has_charges(crystal::Crystal) = crystal.charges.n > 0

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

vtk_filename(crystal::Crystal) = replace(replace(crystal.name, ".cif" => ""), ".cssr" => "") * ".vtk"
write_vtk(crystal::Crystal; center_at_origin::Bool=false) = write_vtk(crystal.box, vtk_filename(crystal), center_at_origin=center_at_origin)

"""
    formula = chemical_formula(crystal, verbose=false)

Find the irreducible chemical formula of a crystal structure.

# Arguments
- `crystal::Crystal`: The crystal containing the crystal structure information
- `verbose::Bool`: If `true`, will print the chemical formula as well

# Returns
- `formula::Dict{Symbol, Int}`: A dictionary with the irreducible chemical formula of a crystal structure
"""
function chemical_formula(crystal::Crystal; verbose::Bool=false)
    unique_atoms = unique(crystal.atoms.species)
    # use dictionary to count atom types
    atom_counts = Dict{Symbol, Int}([a => 0 for a in unique_atoms])
    for i = 1:crystal.atoms.n
        atom_counts[crystal.atoms.species[i]] += 1
    end

    # get greatest common divisor
    gcd_ = gcd([k for k in values(atom_counts)])

    # turn into irreducible chemical formula
    for atom in keys(atom_counts)
        atom_counts[atom] = atom_counts[atom] / gcd_
    end

    # print result
    if verbose
        @printf("Chemical formula of %s:\n\t", crystal.name)
        for (atom, formula_unit) in atom_counts
			@printf("%s_%d", string(atom), formula_unit)
        end
        @printf("\n")
    end

    return atom_counts
end

"""

    mass_of_crystal = molecular_weight(crystal)

Calculates the molecular weight of a unit cell of the crystal in amu using information stored in `data/atomicmasses.csv`.

# Arguments
- `crystal::Crystal`: The crystal containing the crystal structure information

# Returns
- `mass_of_crystal::Float64`: The molecular weight of a unit cell of the crystal in amu
"""
function molecular_weight(crystal::Crystal)
    atomic_masses = read_atomic_masses()

    mass = 0.0
	for i = 1:crystal.atoms.n
        mass += atomic_masses[crystal.atoms.species[i]]
    end

    return mass # amu
end

"""
    ρ = crystal_density(crystal) # kg/m²

Compute the crystal density of a crystal. Pulls atomic masses from [`read_atomic_masses`](@ref).

# Arguments
- `crystal::Crystal`: The crystal containing the crystal structure information

# Returns
- `ρ::Float64`: The crystal density of a crystal in kg/m³
"""
crystal_density(crystal::Crystal) = molecular_weight(crystal) / crystal.box.Ω * 1660.53892  # --> kg/m3

"""
    simulation_ready_crystal = apply_symmetry_operations(non_p1_crystal)

Convert a crystal to P1 symmetry based on internal symmetry rules. This will return a new crystal.

# Arguments
- `non_p1_crystal::Crystal`: The crystal to be converted to P1 symmetry

# Returns
- `P1_crystal::Crystal`: The crystal after it has been converted to P1
    symmetry. The new symmetry rules will be the P1 symmetry rules
"""
function apply_symmetry_operations(crystal::Crystal)
    if crystal.symmetry.is_p1
        return crystal
    end

    nb_symmetry_ops = size(crystal.symmetry.operations, 2)

    n_atoms = crystal.atoms.n * nb_symmetry_ops
    atoms = Atoms{Frac}(n_atoms)

    n_charges = crystal.charges.n * nb_symmetry_ops
    charges = Charges{Frac}(n_charges)

    # for each symmetry rule
    for sr in 1:nb_symmetry_ops
        # loop over all atoms in lower level symmetry
        sym_rule = eval.(Meta.parse.("(x, y, z) -> " .* crystal.symmetry.operations[:, sr]))
        for a in 1:crystal.atoms.n
            # apply current symmetry rule to current atom for x, y, and z coordinates
            atom_id = (sr - 1) * crystal.atoms.n + a
            atoms.species[atom_id] = crystal.atoms.species[a]
            atoms.coords.xf[:, atom_id] .= [Base.invokelatest.(
                        sym_rule[k], crystal.atoms.coords.xf[:, a]...) for k in 1:3]
        end

        # loop over all charges in lower level symmetry
        for c in 1:crystal.charges.n
            # apply current symmetry rule to current atom for x, y, and z coordinates
            charge_id = (sr - 1) * crystal.charges.n + c
            charges.q[c] = crystal.charges.q[c]
            charges.coords.xf[:, charge_id] .= [Base.invokelatest.(
                        sym_rule[k], crystal.charges.coords.xf[:, c]...) for k in 1:3]
        end
    end

    return Crystal(crystal.name, crystal.box, atoms, charges)
end

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
"""
    assert_P1_symmetry(crystal::Crystal)

Throw an error if and only if the crystal is not in P1 symmetry.
"""
function assert_P1_symmetry(crystal::Crystal)
    if ! crystal.symmetry.is_p1
        error("the crystal " * crystal.name * " is not in P1 symmetry.\n
               To convert to P1 symmetry, try:\n
               \tcrystal_p1 = apply_symmetry_operations(crystal)")
    end
end

"""
    write_cif(crystal, filename; fractional_coords=true, number_atoms=true)

Write a `crystal::Crystal` to a .cif file with `filename::AbstractString`. If `filename` does
not include the .cif extension, it will automatically be added. the `fractional_coords` flag
allows us to write either fractional or Cartesian coordinates.
"""
function write_cif(crystal::Crystal, filename::AbstractString; fractional_coords::Bool=true,
		   number_atoms::Bool=true)
    if has_charges(crystal)
        if crystal.atoms.n != crystal.charges.n
            error("write_cif assumes equal numbers of Charges and Atoms (or zero charges)")
        end
        if ! isapprox(crystal.charges.coords, crystal.atoms.coords)
            error("write_cif needs coords of atoms and charges to correspond.")
        end
    end

    # TODO is this labeling necessary for the bonds, arthur?
    # create dictionary for tracking label numbers
    label_numbers = Dict{Symbol, Int}()
    for atom in crystal.atoms.species
        if !haskey(label_numbers, atom)
            label_numbers[atom] = 1
        end
    end

    # append ".cif" to filename if it doesn't already have the extension
    if ! occursin(".cif", filename)
        filename *= ".cif"
    end
    cif_file = open(filename, "w")
    # first line should be data_xtalname_PM
    if crystal.name == ""
        @printf(cif_file, "data_PM\n")
    else
        # don't include file extension!
        @printf(cif_file, "data_%s_PM\n", split(crystal.name, ".")[1])
    end

    @printf(cif_file, "_symmetry_space_group_name_H-M\t'%s'\n", crystal.symmetry.space_group)

    @printf(cif_file, "_cell_length_a\t%f\n", crystal.box.a)
    @printf(cif_file, "_cell_length_b\t%f\n", crystal.box.b)
    @printf(cif_file, "_cell_length_c\t%f\n", crystal.box.c)

    @printf(cif_file, "_cell_angle_alpha\t%f\n", crystal.box.α * 180.0 / pi)
    @printf(cif_file, "_cell_angle_beta\t%f\n", crystal.box.β * 180.0 / pi)
    @printf(cif_file, "_cell_angle_gamma\t%f\n", crystal.box.γ * 180.0 / pi)

    @printf(cif_file, "_symmetry_Int_Tables_number 1\n\n")
    @printf(cif_file, "loop_\n_symmetry_equiv_pos_as_xyz\n")
    for i in 1:size(crystal.symmetry.operations, 2)
        @printf(cif_file, "'%s,%s,%s'\n", crystal.symmetry.operations[:, i]...)
    end
    @printf(cif_file, "\n")

    @printf(cif_file, "loop_\n_atom_site_label\n_atom_site_type_symbol\n")
    if fractional_coords
        @printf(cif_file, "_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n")
    else
        @printf(cif_file, "_atom_site_Cartn_x\n_atom_site_Cartn_y\n_atom_site_Cartn_z\n")
    end
    high_precision_charges = false # if, for neutrality, need high-precision charges
    if has_charges(crystal)
        @printf(cif_file, "_atom_site_charge\n")
        # if crystal will not be charge neutral to a 1e-5 tolerance when loading it
        #    into PorousMaterials.jl, then write higher-precision charges
        if abs(sum(round.(crystal.charges.q, digits=6))) > NET_CHARGE_TOL
            @info "writing high-precision charges for " * filename * ".\n"
            high_precision_charges = true
        end
    end

    idx_to_label = Array{AbstractString, 1}(undef, crystal.atoms.n)
    for i = 1:crystal.atoms.n
        # print label and type symbol
        @printf(cif_file, "%s\t%s\t", string(crystal.atoms.species[i]) *
                (number_atoms ? string(label_numbers[crystal.atoms.species[i]]) : ""),
                crystal.atoms.species[i])
        # store label for this atom idx
        idx_to_label[i] = string(crystal.atoms.species[i]) *
                    string(label_numbers[crystal.atoms.species[i]])
        # increment label
        label_numbers[crystal.atoms.species[i]] += 1
        if fractional_coords
            @printf(cif_file, "%f\t%f\t%f", crystal.atoms.coords.xf[:, i]...)
        else
            @printf(cif_file, "%f\t%f\t%f", (crystal.box.f_to_c * crystal.atoms.coords.xf[:, i])...)
        end
        if has_charges(crystal)
            if high_precision_charges
                @printf(cif_file, "\t%.10f\n", crystal.charges.q[i])
            else
                @printf(cif_file, "\t%f\n", crystal.charges.q[i])
            end
        else
            @printf(cif_file, "\n")
        end
    end

    # only print bond information if it is in the crystal
    if ne(crystal.bonds) > 0
        if ! number_atoms
             error("must label atoms with numbers to write bond information.\n")
        end
        # print column names for bond information
        @printf(cif_file, "\nloop_\n_geom_bond_atom_site_label_1\n_geom_bond_atom_site_label_2\n_geom_bond_distance\n")

        for edge in collect(edges(crystal.bonds))
            dxf = crystal.atoms.coords.xf[:, edge.src] - crystal.atoms.coords.xf[:, edge.dst]
            nearest_image!(dxf)
            @printf(cif_file, "%s\t%s\t%0.5f\n", idx_to_label[edge.src], idx_to_label[edge.dst],
                    norm(dxf))
        end
    end
    close(cif_file)
end

write_cif(crystal::Crystal) = write_cif(crystal, crystal.name)

"""
    crystal = remove_duplicate_atoms_and_charges(crystal, r_tol=0.1, q_tol=0.0001, verbose=false)

Remove duplicate atoms and charges from a crystal. See [`remove_duplicates`](@ref).

# arguments
- `crystal::Crystal`: a crystal
- `r_tol::Float64`: atoms/charges are overlapping if within `r_tol` distance (PBC applied)
- `q_tol::Float64`: charges have the same charge value if their charges are within `q_tol` of each other
- `verbose::Bool`: print off which dupicates are removed

# returns
- `crystal::Crystal`: crystal with duplicate atoms and charges removed.
"""
function remove_duplicate_atoms_and_charges(crystal::Crystal, r_tol::Float64=0.1, q_tol::Float64=0.0001)
    if ne(crystal.bonds) != 0
        error("cannot remove duplicates with bonds assigned")
    end
    atoms = remove_duplicates(crystal.atoms, crystal.box, true, r_tol=r_tol, q_tol=q_tol)
    charges = remove_duplicates(crystal.charges, crystal.box, true, r_tol=r_tol, q_tol=q_tol)
    return Crystal(crystal.name, crystal.box, atoms, charges, SimpleGraph(crystal.atoms.n), crystal.symmetry)
end

inside(crystal::Crystal) = inside(crystal.atoms.coords) && inside(crystal.charges.coords)

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

function Base.show(io::IO, crystal::Crystal)
    println(io, "Name: ", crystal.name)
    println(io, crystal.box)
	@printf(io, "\t# atoms = %d\n", crystal.atoms.n)
	@printf(io, "\t# charges = %d\n", crystal.charges.n)
    println(io, "\tchemical formula: ", chemical_formula(crystal))
    @printf(io, "\tspace Group: %s\n", crystal.symmetry.space_group)
    @printf(io, "\tsymmetry Operations:\n")
    for i in 1:size(crystal.symmetry.operations, 2)
        @printf(io, "\t\t'%s, %s, %s'\n", crystal.symmetry.operations[:, i]...)
    end
end

function Base.isapprox(c1::Crystal, c2::Crystal; atol::Real=0.0)
    # name not included
    box_flag = isapprox(c1.box, c2.box)
    if c1.charges.n != c2.charges.n
        return false
    end
    if c1.atoms.n != c2.atoms.n
        return false
    end
    charges_flag = isapprox(c1.charges, c2.charges, atol=atol)
    atoms_flag = isapprox(c1.atoms, c2.atoms, atol=atol)
    sym1 = c1.symmetry
    sym2 = c2.symmetry
    symmetry_flag = (sym1.is_p1 == sym2.is_p1) && (sym1.space_group == sym1.space_group) && (sym1.operations == sym2.operations)
    return box_flag && charges_flag && atoms_flag && symmetry_flag
end

function Base.getindex(crystal::Crystal, ids::Union{BitArray{1}, Array{Int, 1},
                                                    UnitRange{Int}})
#    if ne(crystal.bonds) != 0
#        error("indexing a crystal with bonds is not supported.")
#    end

    old_to_new = Dict(ids[i] => i for i in 1:length(ids))
    bonds = SimpleGraph(length(ids))
    for edge in collect(edges(crystal.bonds))
        if edge.src in ids && edge.dst in ids
            add_edge!(bonds, old_to_new[edge.src], old_to_new[edge.dst])
        end
    end

    if crystal.charges.n == 0
        return Crystal(crystal.name, crystal.box, crystal.atoms[ids],
            crystal.charges, bonds, crystal.symmetry)
    elseif (crystal.charges.n == crystal.atoms.n) && isapprox(crystal.charges.coords,
                                                              crystal.atoms.coords)
        return Crystal(crystal.name, crystal.box, crystal.atoms[ids],
            crystal.charges[ids], bonds, crystal.symmetry)
    else
        error("for getindex(crystal), crystal must have 0 charges or an equal number of charges and atoms
        that share coordinates")
    end

    return crystal
end

function Base.lastindex(crystal::Crystal)
    if (crystal.atoms.n == crystal.charges.n) || (crystal.charges.n == 0)
        return crystal.atoms.n
    else
        error("to index the crystal, it must have 0 charges or an equal number of charges and atoms")
    end
end

function Base.:+(crystals::Crystal...; check_overlap::Bool=true)
    crystal = deepcopy(crystals[1])
    for (i, f) in enumerate(crystals)
        if i == 1
            continue
        end
        @assert isapprox(crystal.box, f.box) @sprintf("Crystal %s has a different box\n", f.name)
 #         @assert is_symmetry_equal(crystal.symmetry, f.symmetry) @sprintf("Crystal %s has different symmetry rules\n", f.name)
        @assert crystal.symmetry != f.symmetry @sprintf("Crystal %s has different symmetry rules\n", f.name)
        @assert crystal.symmetry.space_group == f.symmetry.space_group

        atoms = crystal.atoms + f.atoms
        charges = crystal.charges + f.charges

        nf_n_atoms = crystal.atoms.n
        add_vertices!(crystal.bonds, nf_n_atoms)
        for edge in collect(edges(f.bonds))
            add_edge!(crystal.bonds, nf_n_atoms + edge.src, nf_n_atoms + edge.dst)
        end

        crystal = Crystal(split(crystal.name, ".")[1] * "_" * split(f.name, ".")[1],
                                 crystal.box, atoms, charges, crystal.bonds,
                                 crystal.symmetry)
    end

    if check_overlap
        _check_overlap(crystal)
    end

    return crystal
end
