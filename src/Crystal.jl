using Base.Test

"""
    framework = Framework(name, a, b, c, α, β, γ, Ω, n_atoms, atoms, xf, f_to_c, c_to_f)

Data structure for a 3D crystal structure.

# Arguments
- `name::String`: corresponds to crystal structure filename from which it was read.
- `a,b,c::Float64`: unit cell dimensions (units: Angstroms)
- `α,β,γ::Float64`: unit cell angles (units: radians)
- `n_atoms::Int64`: number of atoms in the unit cell
- `Ω::Float64`: volume of the unit cell (units: cubic Angtroms)
- `atoms::Array{String,1}`: list of (pseudo)atoms (e.g. elements) composing crystal unit 
cell, in strict order
- `xf::Array{Float64,2}`: a 2D array of fractional coordinates of the atoms, in strict 
order corresponding to `atoms`, stored column-wise so that xf[:, 1] is first atom's 
fractional coordinates.
- `f_to_c::Array{Float64,2}`: a 3x3 matrix used to convert fractional coordinates to 
cartesian coordinates
- `c_to_f::Array{Float64,2}`: a 3x3 matrix used to convert cartesian coordinates to 
fractional coordinates
"""
struct Framework
    name::String

    a::Float64
    b::Float64
    c::Float64

    α::Float64
    β::Float64
    γ::Float64

    Ω::Float64

    n_atoms::Int64
    atoms::Array{String, 1}
    xf::Array{Float64, 2}

    f_to_C::Array{Float64, 2}
    C_to_f::Array{Float64, 2}
end

"""
    framework = read_crystal_structure_file("filename.cssr"; run_checks=true)

Read a crystal structure file (.cif or .cssr) and construct a Framework object.
If `run_checks=True`, ensures no atom overlap or atoms outside of the unit cell box.
"""
function read_crystal_structure_file(filename::String; run_checks::Bool=true)
    # TODO add charges
    # read file extension, ensure reader implemented.
    extension = split(filename, ".")[end]
    if ! (extension in ["cif", "cssr"])
        error("PorousMaterials.jl can only read .cif or .cssr crystal structure files.")
    end
    
    # read in crystal structure file
    if ! isfile(PATH_TO_DATA * "crystals/" * filename)
        error(@sprintf("Could not open crystal structure file %s\n", 
                       PATH_TO_DATA * "crystals/" * filename))
    end
    f = open(PATH_TO_DATA * "crystals/" * filename, "r")
    lines = readlines(f)
    close(f)
    # TODO This reader is a mishmash of a lot of things. This needs to be cleaned up!

    if extension == "cssr"
        # Initialize variables
        charges = Float64[]
        n_atoms = length(lines) - 5
        a, b, c, α, β, γ = Array{Float64}(6)
        x = Array{Float64}(n_atoms) # fractional
        y = similar(x)
        z = similar(x)
        atoms = Array{String}(n_atoms)

        # Make a boolean variable to fix discrepency with fractional and cartesian 
        #  coordinates in cssr files
        corr = false
        for line in lines[6:6+16]
            str = split(line)
            for val in str[3:5]
                # Check if values are fractional or not. I only check 10 lines, 
                #  which is probably enough.
                if parse(Float64,val)>1
                    corr = true
                    break
                end
            end
        end

        # Iterate through the lines of the cssr file.
        for (i,line) in enumerate(lines)
            str = split(line)
            # Unit cell dimension line
            if i == 1
                a, b, c = map(x->parse(Float64, x), str[end-2:end])
            # Unit cell angle line
            elseif i == 2
                temp = zeros(3)
                cnt = 1
                for val in str
                    try
                        temp[cnt] = parse(Float64,val)*π/180
                        cnt += 1
                    end
                end
                α, β, γ = temp[1:3]
            # Atom lines
            elseif i > 5
                try # Fix faulty cssr files where columns one and two merge for some reason
                    parse(Float64,str[1])
                catch
                    str = correct_cssr_line(str)
                end

                # Fix element column where some cssr files have a number concated to the 
                #  end of the element tag (e.g. H3)
                tempch = ""
                for ch in str[2]
                    if !isdigit(ch)
                        tempch = string(tempch,ch)
                    end
                end
                atoms[i - 5] = tempch
                # Use the correction boolean to fix fractional/cartesian discrepency
                if corr
                    x[i - 5], y[i - 5], z[i - 5] = map(x->parse(Float64, x), str[3:5])./[a, b, c]
                    if i == 6
                        @printf("x = %f, y = %f, z = %f\n",x[i-5],y[i-5],z[i-5])
                    end
                else
                    x[i - 5], y[i - 5], z[i - 5] = map(x->parse(Float64, x), str[3:5])
                end
            end
        end
    end

    # Cif reader from Cory Simon
    if extension == "cif"
        data = Dict()
        charges = Float64[]
        x = Float64[]
        y = Float64[]
        z = Float64[]
        atoms = String[]

        loop_starts = -1
        for i = 1:length(lines)
            line = split(lines[i])
            if length(line) == 0
                continue
            end
                        
            if line[1] == "_symmetry_space_group_name_H-M"
                # TODO long-term write our own replicator?
                if length(line) == 3
                    @assert(contains(line[2] * line[3], "P1") || contains(line[2] * line[3], "P 1"), ".cif must have P1 symmetry.\n")
                elseif length(line) == 2
                    @assert(contains(line[2], "P1") || contains(line[2], "P 1"), ".cif must have P1 symmetry\n")
                else
                    println(line)
                    error("Does this .cif have P1 symmetry? Use `convert_cif_to_P1_symmetry` to convert to P1 symmetry")
                end
            end

            for axis in ["a", "b", "c"]
                if line[1] == @sprintf("_cell_length_%s", axis)
                    data[axis] = parse(Float64, line[2])
                end
            end

            for angle in ["alpha", "beta", "gamma"]
                if line[1] == @sprintf("_cell_angle_%s", angle)
                    data[angle] = parse(Float64, line[2]) * pi / 180.0
                end
            end

            if (line[1] == "loop_")
                next_line = split(lines[i+1])
                if contains(next_line[1], "_atom_site")
                    loop_starts = i + 1
                    break
                end
            end
        end  # end loop over lines

        if loop_starts == -1
            error("Could not find _atom_site* after loop_ if .cif file\n")
        end

        # broke the loop. so loop_starts is line where "_loop" first starts
        # name_to_column is a dictionary that e.g. returns which column contains x fractional coord
        #   use example: name_to_column["_atom_site_fract_x"] gives 3
        atom_column_name = "_atom_site_type_symbol"  # this can be different for different .cifs...
        name_to_column = Dict{AbstractString, Int}()
        i = loop_starts
        while length(split(lines[i])) == 1
            if i == loop_starts
                atom_column_name = split(lines[i])[1]
            end
            name_to_column[split(lines[i])[1]] = i + 1 - loop_starts
            i += 1
        end

        # now extract fractional coords of atoms and their charges
        for i = loop_starts+length(name_to_column):length(lines)
            line = split(lines[i])
            if length(line) != length(name_to_column)
                break
            end
            push!(atoms, line[name_to_column[atom_column_name]])
            push!(x, mod(parse(Float64, line[name_to_column["_atom_site_fract_x"]]), 1.0))
            push!(y, mod(parse(Float64, line[name_to_column["_atom_site_fract_y"]]), 1.0))
            push!(z, mod(parse(Float64, line[name_to_column["_atom_site_fract_z"]]), 1.0))
            # if charges present, import them
            if haskey(name_to_column, "_atom_site_charge")
                push!(charges, parse(Float64, line[name_to_column["_atom_site_charge"]]))
            else
                push!(charges, 0.0)
            end
        end
        n_atoms = length(x)
        a = data["a"]
        b = data["b"]
        c = data["c"]
        α = data["alpha"]
        β = data["beta"]
        γ = data["gamma"]
    end

    Ω = a * b * c * sqrt(1 - cos(α) ^ 2 - cos(β) ^ 2 - cos(γ) ^ 2 + 2 * cos(α) * cos(β) * cos(γ)) # unit cell volume
    f_to_C = [[a, 0, 0] [b * cos(γ), b * sin(γ), 0] [c * cos(β), c * (cos(α) - cos(β) * cos(γ)) / sin(γ), Ω / (a * b * sin(γ))]]
    C_to_f = [[1/a, 0, 0] [-cos(γ) / (a * sin(γ)), 1 / (b * sin(γ)), 0] [b * c * (cos(α) * cos(γ) - cos(β)) / (Ω * sin(γ)), a * c * (cos(β) * cos(γ) - cos(α)) / (Ω * sin(γ)), a * b * sin(γ) / Ω]]
    
    @test f_to_C * C_to_f ≈ eye(3)

    fractional_coords = Array{Float64, 2}(3, length(x))
    fractional_coords[1, :] = x[:]; fractional_coords[2, :] = y[:]; fractional_coords[3, :] = z[:]
    
    # finally construct the framework
    framework = Framework(filename, a, b, c, α, β, γ, Ω, n_atoms, atoms, fractional_coords, f_to_C, C_to_f)
    
    if run_checks
        check_for_atom_overlap(framework)
        check_for_atoms_out_of_unit_cell_box(framework)
    end

    return framework
end # constructframework end

"""
    replicate_to_xyz(framework, xyzfilename; comment="", repfactors=(0, 0, 0),
                     negative_replications=false)

Write a .xyz file representation of a crystal structure from a Framework type.
Write an optional `comment` to the .xyz file if desired.
Extend the structure in the x-,y- or z-direction by changing the tuple `repfactors`.
A value of 1 replicates the structure once in the desired direction, so 
`repfactors=(0, 0, 0)` includes only the "home" unit cell. Ensure `negative_replications`
is true if home unit cell should be replicated in the negative directions too.
"""
function replicate_to_xyz(framework::Framework, xyzfilename::String; 
                          comment::String="", repfactors::Tuple{Int, Int, Int}=(0, 0, 0),
                          negative_replications::Bool=false)
    # pre-calculate # of total atoms in .xyz
    if negative_replications
        n_atoms = framework.n_atoms * (2 * repfactors[1] + 1) * (2 * repfactors[2] + 1) * (2 * repfactors[3] + 1)
        neg_repfactors = (-repfactors[1], -repfactors[2], -repfactors[3])
    else
        n_atoms = framework.n_atoms * (repfactors[1] + 1) * (repfactors[2] + 1) * (repfactors[3] + 1)
        neg_repfactors = (0, 0, 0)
    end

    f = open(xyzfilename, "w")
    @printf(f, "%d\n%s\n", n_atoms, comment)

    for i = neg_repfactors[1]:repfactors[1], j = neg_repfactors[2]:repfactors[2], k = neg_repfactors[3]:repfactors[3]
        xf = framework.xf .+ [i, j, k]
        c_coords = framework.f_to_C * xf
        for ii = 1:size(c_coords, 2)
            @printf(f, "%s\t%.4f\t%.4f\t%.4f\n", framework.atoms[ii], c_coords[1, ii], c_coords[2, ii], c_coords[3, ii])
        end
    end
    close(f)

    println("See ", xyzfilename)
    return
end # replicate_to_xyz end

# TODO: remove this once we figure out how to read .cssr another way.
"""
    corrected_line = correct_cssr_line(str)

Take an array of string values and correct an error from the openbabel python module.
The error merges two values together if the element abbreviation contains two letters,
such as Zn, Cl and so on.
"""
function correct_cssr_line(str)
    tempbool = Array{Bool}(length(str[1]),1)
    for (k,ch) in enumerate(str[1])
        tempbool[k] = isdigit(ch)
    end

    ind = findfirst(tempbool,false);
    unshift!(str,str[1][1:ind-1])
    str[2] = str[2][ind:end]

    for (k,ch) in enumerate(str[2])
        if isdigit(ch)
            str[2] = str[2][1:k-1]
            break
        end
    end
    return str
end #correct_cssr_line end

"""
    check_for_atom_overlap(framework; threshold_distance_for_overlap=0.1, verbose=false)

Check if any two atoms are lying on top of each other by calculating the 2-norm distance
between every pair of atoms and ensuring distance is greater than a threshold.
Throw error if atoms overlap and tell which atoms are culprit.
"""
function check_for_atom_overlap(framework::Framework; 
                                threshold_distance_for_overlap::Float64=0.1,
                                verbose::Bool=false)
    for i = 1:framework.n_atoms
        xf_i = framework.xf[:, i]
        for j = i+1:framework.n_atoms
            xf_j = framework.xf[:, j]
            # vector pointing from atom j to atom i in carteisan coords
            dx = framework.f_to_C * (xf_i - xf_j)
            if (norm(dx) < threshold_distance_for_overlap)
                error(@sprintf("Atoms %d and %d are too close, distance %f Å) < %f Å threshold\n", i, j, norm(dx), threshold_distance_for_overlap))
            end
        end
    end
    if verbose
        @printf("No atoms are on top of each other!")
    end
    return
end

"""
    check_for_atoms_out_of_unit_cell_box(framework)

Check if any atoms lie outside of the unit cell box and print a warning if so.
"""
function check_for_atoms_out_of_unit_cell_box(framework::Framework)
    for i = 1:framework.n_atoms
        if sum(framework.xf[:, i] .< 0.0) != 0
            @printf("WARNING: Atom %d has fractional coordinate less than 0.0\n", i)
        end
        if sum(framework.xf[:, i] .> 1.0) != 0
            @printf("WARNING: Atom %d has fractional coordinate greater than 1.0\n", i)
        end
    end
    return true
end

"""
    strip_numbers_from_atom_labels(framework::Framework)

Strip numbers from labels for `framework.atoms`.
e.g. C12 --> C
"""
function strip_numbers_from_atom_labels!(framework::Framework)
    for i = 1:framework.n_atoms
        while !isalpha(framework.atoms[i][end])
            framework.atoms[i] = chop(framework.atoms[i])
        end
    end
    return 
end

"""       
	write_unitcell_boundary_vtk(framework::Framework, filename::String)

Write unit cell boundary as a .vtk file for visualizing the unit cell boundary.    
"""                                                                             
function write_unitcell_boundary_vtk(framework::Framework, filename::String)
	vtk_file = open(filename, "w")                          
																					
	@printf(vtk_file, "# vtk DataFile Version 2.0\nunit cell boundary\n             
					   ASCII\nDATASET POLYDATA\nPOINTS 8 double\n")                 
																					
	# write points on boundary of unit cell                                         
	for i = 0:1                                                                     
		for j = 0:1                                                                 
			for k = 0:1                                                             
				xf = [i, j, k] # fractional coordinate                               
				cellpoint = framework.f_to_C * xf                    
				@printf(vtk_file, "%.3f %.3f %.3f\n", 
						cellpoint[1], cellpoint[2], cellpoint[3])
			end                                                                     
		end                                                                         
	end                                                                             
																					
	# define connections                                                            
	@printf(vtk_file, "LINES 12 36\n2 0 1\n2 0 2\n2 1 3\n2 2 3\n2 4 5\n             
					   2 4 6\n2 5 7\n2 6 7\n2 0 4\n2 1 5\n2 2 6\n2 3 7\n")          
	close(vtk_file)                                                                 
	println("See ", filename)
	return
end              

    
function chemical_formula(framework::Framework)
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
    @printf("Chemical formula of %s:\n\t", framework.name)
    for (atom, formula_unit) in atom_counts
        @printf("%s_%d", atom, formula_unit)
    end
    @printf("\n")

    return atom_counts
end

"""
    convert_cif_to_P1_symmetry(filename::String, outputfilename::String)

Use Atomic Simulation Environment (ASE) Python package to convert .cif file in non-P1 
symmetry to P1 symmetry. Writes .cif with P1 symmetry to `outputfilename1`.
"""
function convert_cif_to_P1_symmetry(filename::String, outputfilename::String)               
    # Import Atomic Simulation Environment Python package                                   
    @pyimport ase
    @pyimport ase.io as aseio
    @pyimport ase.build as asebuild
                                                                                            
    non_p1_cif = aseio.read(filename, format="cif")
    p1_cif = asebuild.make_supercell(non_p1_cif, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    aseio.write(outputfilename, p1_cif, format="cif")

    return
end 
