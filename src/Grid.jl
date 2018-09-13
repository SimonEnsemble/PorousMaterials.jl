"""
Data structure for a regular [equal spacing between points in each coordinate] grid of points superimposed on a unit cell box (`Box`).
Each grid point has data, `data`, associated with it, of type `T`, stored in a 3D array.

# Attributes
- `box::Box`: describes Bravais lattice over which a grid of points is super-imposed. grid points on all faces are included.
- `n_pts::Tuple{Int, Int, Int}`: number of grid points in x, y, z directions. 0 and 1 fractional coordinates are included.
- `data::Array{T, 3}`: three dimensional array conaining data associated with each grid point.
- `units::Symbol`: the units associated with each data point.
- `origin::Array{Float64, 1}`: the origin of the grid.
"""
struct Grid{T}
    box::Box
    n_pts::Tuple{Int, Int, Int}
    data::Array{T, 3}
    units::Symbol
    origin::Array{Float64, 1}
end

"""
    write_cube(grid, filename, verbose=true)

Write grid to a .cube file format. This format is described here:
http://paulbourke.net/dataformats/cube/
The atoms of the unit cell are not printed in the .cube. Instead, use .xyz files to also visualize atoms.

# Arguments
- `grid::Grid`: grid with associated data at each grid point.
- `filename::AbstractString`: name of .cube file to which we write the grid; this is relative to `PorousMaterials.PATH_TO_DATA`/grids/.
- `verbose::Bool`: print name of file after writing.
"""
function write_cube(grid::Grid, filename::AbstractString; verbose::Bool=true)
    if ! isdir(PATH_TO_DATA * "grids/")
        mkdir(PATH_TO_DATA * "grids/")
    end

    if ! occursin(".cube", filename)
        filename = filename * ".cube"
    end
    cubefile = open(PATH_TO_DATA * "grids/" * filename, "w")

    @printf(cubefile, "Units of data: %s\nLoop order: x, y, z\n", grid.units)

    # the integer refers to 0 atoms (just use .xyz to visualize atoms)
    # the next three floats correspond to the origin
    @printf(cubefile, "%d %f %f %f\n" , 0, grid.origin...)
    for k = 1:3
        # these are the vectors that form the parallelogram comprising the voxels
        # 0 and 1 fractional coords were included. so voxel vector is unit cell axis divided by # grid pts - 1
        voxel_vector = grid.box.f_to_c[:, k] / (grid.n_pts[k] - 1)
        @printf(cubefile, "%d %f %f %f\n" , grid.n_pts[k],
            voxel_vector[1], voxel_vector[2], voxel_vector[3])
    end

    for i = 1:grid.n_pts[1]
        for j = 1:grid.n_pts[2]
            for k = 1:grid.n_pts[3]
                @printf(cubefile, "%e ", grid.data[i, j, k])
                if (k % 6) == 0
                    @printf(cubefile, "\n")
                end
            end # loop over z points
            @printf(cubefile, "\n")
        end # loop over y points
    end # loop over x points
    close(cubefile)
    if verbose
        println("See ", PATH_TO_DATA * "grids/" * filename)
    end
    return
end

"""
    grid = read_cube(filename)

Read a .cube file and return a populated `Grid` data structure.

# Arguments
- `filename::AbstractString`: name of .cube file to which we write the grid; this is relative to `PorousMaterials.PATH_TO_DATA`grids/.

# Returns
- `grid::Grid`: A grid data structure
"""
function read_cube(filename::AbstractString)
    if ! occursin(".cube", filename)
        filename *= ".cube"
    end

    cubefile = open(PATH_TO_DATA * "grids/" * filename)

    # waste two lines
    line = readline(cubefile)
    units = Symbol(split(line)[4])
    readline(cubefile)

    # read origin
    line = readline(cubefile)
    origin = [parse(Float64, split(line)[1 + i]) for i = 1:3]

    # read box info
    box_lines = [readline(cubefile) for i = 1:3]

    # number of grid pts
    n_pts = Tuple([parse(Int, split(box_lines[i])[1]) for i = 1:3])

    # f_to_c matrix (given by voxel vectors)
    f_to_c = zeros(Float64, 3, 3)
    for i = 1:3, j=1:3
        f_to_c[j, i] = parse(Float64, split(box_lines[i])[j + 1])
    end
    for k = 1:3
        f_to_c[:, k] = f_to_c[:, k] * (n_pts[k] - 1.0)
    end

    # reconstruct box from f_to_c matrix
    box = Box(f_to_c)

    # read in data
    data = zeros(Float64, n_pts...)
    line = readline(cubefile)
    for i = 1:n_pts[1]
        for j = 1:n_pts[2]
            read_count = 0
            for k = 1:n_pts[3]
                data[i, j, k] = parse(Float64, split(line)[read_count+1])
                read_count += 1
                if (k % 6 == 0)
                    line = readline(cubefile)
                    read_count = 0
                end
            end # loop over z points
            line = readline(cubefile)
        end # loop over y points
    end # loop over x points

    return Grid(box, n_pts, data, units, origin)
end

"""
	grid = energy_grid(framework, molecule, ljforcefield; n_pts=(50, 50, 50), temperature=298.0, n_rotations=750)

Superimposes a regular grid of points (regularly spaced in fractional coordinates of the `framework.box`) over the unit cell of a crystal, with `n_gridpts` dictating the number of grid points in the a, b, c directions (including 0 and 1 fractional coords).
The fractional coordinates 0 and 1 are included in the grid, although they are redundant.
Then, at each grid point, calculate the ensemble average potential energy of the molecule when its mass is centered at that point. The average is taken over Boltzmann-weighted rotations.

The ensemble average is a Boltzmann average over rotations:  - R T log ⟨e⁻ᵇᵁ⟩

# Arguments
- `framework::Framework`: crystal in which we seek to compute an energy grid for a molecule. `grid.box` will be `framework.box`.
- `molecule::Molecule`: molecule for which we seek an energy grid
- `ljforcefield::LJForceField`: molecular model for computing molecule-framework interactions
- `n_pts::Tuple{Int, Int, Int}=(50,50,50)`: number of grid points in each fractional coordinate dimension, including endpoints (0, 1)
- `n_rotations::Int`: number of random rotations to conduct in a Monte Carlo simulation for finding the free energy of a molecule centered at a given grid point.
This is only relevant for molecules that are comprised of more than one Lennard Jones sphere.
- `temperature::Float64`: the temperature at which to compute the free energy for molecules where rotations are required. Lower temperatures overemphasize the minimum potential energy rotational conformation at that point.
- `units::Symbol`: either `:K` or `:kJ_mol`, the units in which the energy should be stored in the returned `Grid`.
- `center::Bool`: shift coords of grid so that the origin is the center of the unit cell `framework.box`.
- `verbose::Bool=true`: print some information.

# Returns
- `grid::Grid`: A grid data structure containing the potential energy of the system
"""
function energy_grid(framework::Framework, molecule::Molecule, ljforcefield::LJForceField;
                     n_pts::Tuple{Int, Int, Int}=(50,50,50), n_rotations::Int=1000, 
                     temperature::Float64=NaN, units::Symbol=:kJ_mol, center::Bool=false,
                     verbose::Bool=true)
    if ! (units in [:kJ_mol, :K])
        error("Pass :kJ_mol or :K for units of kJ/mol or K, respectively.")
    end

    rotations_required = rotatable(molecule)
    charged_system = charged(framework) && charged(molecule)
    if rotations_required & isnan(temperature)
        error("Must pass temperature (K) for Boltzmann weighted rotations.\n")
    end

    eparams = setup_Ewald_sum(framework.box, sqrt(ljforcefield.cutoffradius_squared),
                              verbose=verbose & charged_system)
    eikr = Eikr(framework, eparams)

    # grid of voxel centers (each axis at least).
    grid_pts = [collect(range(0.0; stop=1.0, length=n_pts[i])) for i = 1:3]

    grid = Grid(framework.box, n_pts, zeros(Float64, n_pts...), units, 
        center ? framework.box.f_to_c * [-0.5, -0.5, -0.5] : [0.0, 0.0, 0.0])

    if verbose
        @printf("Computing energy grid of %s in %s\n", molecule.species, framework.name)
        @printf("\tRegular grid (in fractional space) of %d by %d by %d points superimposed over the unit cell.\n", n_pts...)
        if rotations_required
            @printf("\t%d molecule rotations per grid point with temperature %f K.\n", n_rotations, temperature)
        end
    end
    
    # compute replication factors required for short-range interactions & cutoff
    repfactors = replication_factors(framework.box, ljforcefield)
    framework = replicate(framework, repfactors)

	for (i, xf) in enumerate(grid_pts[1]), (j, yf) in enumerate(grid_pts[2]), (k, zf) in enumerate(grid_pts[3])
        # must account for fact that framework is now replicated; use coords in home box
        translate_to!(molecule, [xf, yf, zf] ./ repfactors)
        if ! rotations_required
            ensemble_average_energy = vdw_energy(framework, molecule, ljforcefield)
        else
            boltzmann_factor_sum = 0.0
            for r = 1:n_rotations
                rotate!(molecule, framework.box)

                energy = PotentialEnergy(0.0, 0.0)
                energy.vdw = vdw_energy(framework, molecule, ljforcefield)
                if charged_system
                    energy.coulomb = total(electrostatic_potential_energy(framework, molecule, 
                                                                         eparams, eikr))
                end
                boltzmann_factor_sum += exp(-sum(energy) / temperature)
            end
            ensemble_average_energy = -temperature * log(boltzmann_factor_sum / n_rotations)
        end
		grid.data[i, j, k] = ensemble_average_energy # K
	end

    if units == :kJ_mol # K - kJ/mol
        grid.data[:] = grid.data[:] * 8.314 / 1000.0
    end

	return grid
end

function Base.show(io::IO, grid::Grid)
    @printf(io, "Regular grid of %d by %d by %d points superimposed over a unit cell and associated data.\n", grid.n_pts...)
    @printf(io, "\tunits of data attribute: %s\n", grid.units)
    @printf(io, "\torigin: [%f, %f, %f]\n", grid.origin...)
end

# comparing very large numbers in grid.data, so increase rtol to account
#  for loss of precision when writing grid.data to a cube file.
function Base.isapprox(g1::Grid, g2::Grid; rtol::Float64=0.000001)
    return (isapprox(g1.box, g2.box, rtol=rtol) &&
            (g1.n_pts == g2.n_pts) &&
            isapprox(g1.data, g2.data, rtol=rtol) &&
            (g1.units == g2.units) &&
            isapprox(g1.origin, g2.origin))
end

function _flood_fill!(grid::Grid, segmented_grid::Grid, i::Int, j::Int, k::Int, energy_tol::Float64)
    # look left, right, up, down "_s" for shift
    for i_s = -1:1, j_s = -1:1, k_s = -1:1
        ind = (i + i_s, j + j_s, k + k_s)
        # well already know segmented_grid[i, j, k]
        if ind == (0, 0, 0)
            continue
        end
        # avoid out of bounds
        if any(ind .<= (0, 0, 0)) || any(ind .> grid.n_pts)
            continue
        end
        # if already assigned, continue
        if segmented_grid.data[ind...] != 0
            continue
        end
        # if accessible, assign same segment ID and keep looking to left and right.
        if grid.data[ind...] < energy_tol
            segmented_grid.data[ind...] = segmented_grid.data[i, j, k]
            # recursively look to left, right, up, down
            _flood_fill!(grid, segmented_grid, ind..., energy_tol)
        end
    end
end

function _segment_grid(grid::Grid; energy_tol::Float64=0.0, verbose::Bool=true)
    # grid of Int's corresponding to each original grid point.
    # let "0" be "unsegmented"
    # let "-1" be "not accessible"
    segmented_grid = Grid(grid.box, grid.n_pts, zeros(Int, grid.n_pts...), :Segment_No, grid.origin)
    segment_no = 0
    for i = 1:grid.n_pts[1], j = 1:grid.n_pts[2], k = 1:grid.n_pts[3]
        if segmented_grid.data[i, j, k] == 0 # not yet assigned segment 
            if grid.data[i, j, k] > energy_tol # not accessible
                segmented_grid.data[i, j, k] = -1
            else # accessible
                # start a new segment!
                segment_no += 1
                segmented_grid.data[i, j, k] = segment_no
                _flood_fill!(grid, segmented_grid, i, j, k, energy_tol)
            end
        else # already segmented
            nothing
        end
    end
    @assert sum(segmented_grid.data == 0) == 0 "some points were unsegmented!"
    verbose ? @printf("Found %d segments\n", segment_no) : nothing
    return segmented_grid
end

function _check_and_merge_if_applicable!(segmented_grid::Grid, left_segment_no::Int,
    right_segment_no::Int)
    # if unaccessible, do no merging.
    if (left_segment_no == -1) || (right_segment_no == -1)
        return nothing
    end
    # at this point, we're only looking at accessible points.

    #  if not the same segment... and accessible... these two segments are connected
    #     and should be merged into one segment with the same segment ID
    if left_segment_no != right_segment_no
        # doesn't matter which we overwrite.
        @printf("Merging segments %d and %d (%d --> %d)\n", left_segment_no, right_segment_no,
            left_segment_no, right_segment_no)
        segmented_grid.data[segmented_grid.data .== left_segment_no] .= right_segment_no
    end
end

function _rename_accessible_segments!(segmented_grid::Grid)
    unique_segments = unique(segmented_grid.data)
    nb_unique_segments = length(unique_segments)
    sort!(unique_segments)
    @printf("Now %d unique segments.\n", nb_unique_segments)
    println("\tUnique segments: ", unique_segments)

    nb_accessible_segments = length(unique_segments)
    # exclude -1
    if -1 in unique_segments
        nb_accessible_segments -= 1
    end

    # rename so that accessible segments are 1, 2, ..., n_unique_segments
    for s = 1:nb_accessible_segments
        if s in unique_segments
            # segment already there, no need to rename
            continue
        end
        # pop largest segment from stack
        replace_segment = pop!(unique_segments)
        # rename this segment 
        @printf("... renaming segment %d to %d\n", replace_segment, s)
        segmented_grid.data[segmented_grid.data .== replace_segment] .= s
    end
    
    new_unique_segments = unique(segmented_grid.data)
    println("new unique segments: ", new_unique_segments)
    @assert(maximum(new_unique_segments) == nb_accessible_segments)
    @assert(length(new_unique_segments) == nb_unique_segments)
end

# find all segment pairs connected across the unit cell boundary
function _merge_segments_connected_across_periodic_boundary!(segmented_grid::Grid)
    # loop over faces of unit cell.
    for i = 1:segmented_grid.n_pts[1], j = 1:segmented_grid.n_pts[2]
        _check_and_merge_if_applicable!(segmented_grid, 
            segmented_grid.data[i, j, 1], segmented_grid.data[i, j, end])
    end
    for j = 1:segmented_grid.n_pts[2], k = 1:segmented_grid.n_pts[3]
        _check_and_merge_if_applicable!(segmented_grid, 
            segmented_grid.data[1, j, k], segmented_grid.data[end, j, k])
    end
    for i = 1:segmented_grid.n_pts[1], k = 1:segmented_grid.n_pts[3]
        _check_and_merge_if_applicable!(segmented_grid, 
            segmented_grid.data[i, 1, k], segmented_grid.data[i, end, k])
    end
    
    # so 1 2 3 4 not e.g. 6 19
    _rename_accessible_segments!(segmented_grid)

    return segmented_grid
end
