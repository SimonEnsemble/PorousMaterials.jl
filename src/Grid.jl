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
    voxel_id = xf_to_id(n_pts, xf)

Returns the indices of the voxel in which it falls when a unit cube is
partitioned into a regular grid of `n_pts[1]` by `n_pts[2]` by `n_pts[3]` voxels.
Periodic boundary conditions are applied.

# Arguments
 - `n_pts::Tuple{Int, Int, Int}`: The number of points for each axis in the `Grid`
 - `xf::Array{Float64, 1}`: The fractional coordinates to be converted to an id
# Returns
 - `id::Array{Int, 1}`: The array indices for storing this point in space
"""
function xf_to_id(n_pts::Tuple{Int, Int, Int}, xf::Array{Float64, 1})
    voxel_id = floor.(Int, xf .* n_pts) .+ 1
    for k = 1:3
        if voxel_id[k] <= 0
            voxel_id[k] += n_pts[k]
        elseif voxel_id[k] > n_pts[k]
            voxel_id[k] -= n_pts[k]
        end
    end
    return voxel_id
end

"""
    xf = id_to_xf(voxel_id, n_pts)

Given a `voxel_id` in a `Grid`, return the fractional coordinates to which this voxel
corresponds.

# Arguments
 - `n_pts::Tuple{Int, Int, Int}`: The number of voxels along each axis in the `Grid`
 - `voxel_id::Array{Int, 1}`: the voxel coordinates in `grid.data`
# Returns
 - `xf::Array{Float64, 1}`: The fractional coordinates corresponding to the grid voxel
"""
function id_to_xf(id::Tuple{Int, Int, Int}, n_pts::Tuple{Int, Int, Int})
    return [(id[k] - 1.0) / (n_pts[k] - 1.0) for k = 1:3]
end 

"""
    update_density!(grid, molecule, species)

updates the density grid based on an array of molecules. If a molecule doesn't
match the specified species it won't be added to the density grid. This function
doesn't calculate the actual densities, it will need a `./ = num_snapshots`
at the end of the GCMC simulation.

# Arguments
 - `grid::Grid`: the grid to be updated
 - `molecules::Array{Molecule, 1}`: An array of molecules whose positions will
    be added to the grid
 - `species::Symbol`: The species of atom that can be added to this density grid
"""
function update_density!(grid::Grid, molecules::Array{Molecule, 1}, species::Symbol)
    for molecule in molecules
        for i = 1:molecule.atoms.n_atoms
            if molecule.atoms.species[i] == species
                voxel_id = xf_to_id(grid.n_pts, molecule.atoms.xf[:, i])
                grid.data[voxel_id...] += 1
            end
        end
    end
end

"""
    n_pts = required_n_pts(box, dx)

Calculate the required number of grid pts in a, b, c unit cell directions required to keep
distances between grid points less than `dx` apart, where `dx` is in units of Angstrom.
"""
function required_n_pts(box::Box, dx::Float64)
    # columns of f_to_c are the unit cell lattice vectors.
    cell_vector_norms = [norm(box.f_to_c[:, i]) for i = 1:3]
    n_pts = zeros(Int, 3)
    for i = 1:3
        n_pts[i] = ceil(Int, cell_vector_norms[i] / dx) + 1
    end
    return Tuple(n_pts)
end

"""
    write_cube(grid, filename, verbose=true)

Write grid to a .cube file format. This format is described here:
http://paulbourke.net/dataformats/cube/
The atoms of the unit cell are not printed in the .cube. Instead, use .xyz files to also visualize atoms.

# Arguments
- `grid::Grid`: grid with associated data at each grid point.
- `filename::AbstractString`: name of .cube file to which we write the grid; this is relative to `PATH_TO_GRIDS`.
- `verbose::Bool`: print name of file after writing.
"""
function write_cube(grid::Grid, filename::AbstractString; verbose::Bool=true)
    if ! isdir(PATH_TO_GRIDS)
        mkdir(PATH_TO_GRIDS)
    end

    if ! occursin(".cube", filename)
        filename = filename * ".cube"
    end
    cubefile = open(joinpath(PATH_TO_GRIDS, filename), "w")

    @printf(cubefile, "Units of data: %s\nLoop order: x, y, z\n", grid.units)

    # the integer refers to 0 atoms (just use .xyz to visualize atoms)
    # the next three floats correspond to the origin
    @printf(cubefile, "%d %f %f %f\n" , 0, grid.origin[1], grid.origin[2], grid.origin[3])
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
        println("\tSee ", joinpath(PATH_TO_GRIDS, filename))
    end
    return
end

"""
    grid = read_cube(filename)

Read a .cube file and return a populated `Grid` data structure.

# Arguments
- `filename::AbstractString`: name of .cube file to which we write the grid; this is relative to `PATH_TO_GRIDS`

# Returns
- `grid::Grid`: A grid data structure
"""
function read_cube(filename::AbstractString)
    if ! occursin(".cube", filename)
        filename *= ".cube"
    end

    cubefile = open(joinpath(PATH_TO_GRIDS, filename))

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
        @printf("\tRegular grid (in fractional space) of %d by %d by %d points superimposed over the unit cell.\n", n_pts[1], n_pts[2], n_pts[3])
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
    @printf(io, "Regular grid of %d by %d by %d points superimposed over a unit cell and associated data.\n", grid.n_pts[1], grid.n_pts[2], grid.n_pts[3])
    @printf(io, "\tunits of data attribute: %s\n", grid.units)
    @printf(io, "\torigin: [%f, %f, %f]\n", grid.origin[1], grid.origin[2], grid.origin[3])
end

# comparing very large numbers in grid.data, so increase rtol to account
#  for loss of precision when writing grid.data to a cube file.
function Base.isapprox(g1::Grid, g2::Grid; atol::Real=0.0, rtol::Real=atol > 0.0 ? 0.0 : sqrt(eps()))
    return (isapprox(g1.box, g2.box, atol=atol, rtol=rtol) &&
            (g1.n_pts == g2.n_pts) &&
            isapprox(g1.data, g2.data, atol=atol, rtol=rtol) &&
            (g1.units == g2.units) &&
            isapprox(g1.origin, g2.origin, atol=atol, rtol=rtol))
end

function _flood_fill!(grid::Grid, segmented_grid::Grid, 
            queue_of_grid_pts::Array{Tuple{Int, Int, Int}, 1},
            i::Int, j::Int, k::Int, energy_tol::Float64)
    # look left, right, up, down "_s" for shift
    for i_s = -1:1, j_s = -1:1, k_s = -1:1
        # well already know segmented_grid[i, j, k]
        if (i_s, j_s, k_s) == (0, 0, 0)
            continue
        end
        
        # index of a neighboring point
        ind = (i + i_s, j + j_s, k + k_s)
        
        # avoid out of bounds
        if any(ind .<= (0, 0, 0)) || any(ind .> grid.n_pts)
            continue
        end
        # if neighbor already assigned or un-occupiable, continue
        if segmented_grid.data[ind...] != 0
            continue
        end

        # if accessible, assign same segment ID and keep looking at surrounding points
        if grid.data[ind...] < energy_tol
            segmented_grid.data[ind...] = segmented_grid.data[i, j, k]
            # add grid pt. to queue to look at ITS surrounding neighbors
            # originally had recursive algo but get stack overflow
            push!(queue_of_grid_pts, ind)
        else
            # if not accessible, assign -1 and don't push any points to the queue
            #  (this is how the algo eventually dies out)
            segmented_grid.data[ind...] = -1
        end
    end
    return nothing
end

function _segment_grid(grid::Grid, energy_tol::Float64, verbose::Bool)
    # grid of Int's corresponding to each original grid point.
    # let "0" be "unsegmented"
    # let "-1" be "not occupiable" and, eventually, "not accessible"
    segmented_grid = Grid(grid.box, grid.n_pts, zeros(Int, grid.n_pts...),
        :Segment_No, grid.origin)
    segment_no = 0
    for i = 1:grid.n_pts[1], j = 1:grid.n_pts[2], k = 1:grid.n_pts[3]
        if segmented_grid.data[i, j, k] == 0 # not yet assigned segment 
            if grid.data[i, j, k] > energy_tol # not accessible
                segmented_grid.data[i, j, k] = -1
            else # accessible
                # start a new segment!
                segment_no += 1
                # initiate a queue of grid points
                queue_of_grid_pts = [(i, j, k)]
                while length(queue_of_grid_pts) != 0
                    # take first index in line
                    id = queue_of_grid_pts[1]
                    # assign segment number
                    segmented_grid.data[id...] = segment_no
                    # look at surroudning points, add to queue if they are also accessible
                    _flood_fill!(grid, segmented_grid, queue_of_grid_pts,
                                 id..., energy_tol)
                    # handled first one in queue, remove from queue
                    deleteat!(queue_of_grid_pts, 1)
                end
            end
        else # already segmented
            nothing
        end
    end
    @assert sum(segmented_grid.data == 0) == 0 "some points were unsegmented!"
    verbose ? @printf("Found %d segments\n", segment_no) : nothing
    return segmented_grid
end

# this returns number of segments in a segmented grid.
#  it excludes the inaccessible portions -1.
function _count_segments(segmented_grid::Grid)
    unique_segments = unique(segmented_grid.data)
    nb_segments = length(unique_segments)
    # not quite; if -1 present (unoccupiable), this shouldn't count as a segment
    if -1 in unique_segments
        nb_segments -= 1
    end
    if nb_segments != 0
        @assert(nb_segments == maximum(segmented_grid.data))
    end
    return nb_segments
end

# struct describing directed edge describing connection between two segments of a flood-
#  filled grid. the direction (which unit cell face is traversed in the connection) is
#  also included.
struct SegmentConnection
    src::Int # source segment
    dst::Int # destination segment
    direction::Tuple{Int, Int, Int} # unit cell face traversed
end

function _note_connection!(segment_1::Int, segment_2::Int, connections::Array{SegmentConnection, 1},
                           direction::Tuple{Int, Int, Int}, verbose::Bool=true)
    # ignore connections between un-occupiable regions
    if (segment_1 == -1) || (segment_2 == -1)
        return nothing
    end
    
    # construct edge; if not in list of edges, push it and note connection
    segment_connection = SegmentConnection(segment_1, segment_2, direction)
    if ! (segment_connection in connections)
        push!(connections, segment_connection)
        if verbose
            @printf("Noted seg. %d --> %d connection in (%d, %d, %d) direction.\n", 
                    segment_1, segment_2, direction[1], direction[2], direction[3])
        end
        # also add opposite direction to take into account symmetry
        push!(connections, SegmentConnection(segment_2, segment_1, segment_connection.direction .* -1))
    end

    return nothing
end

# obtain set of Segment Connections
function _build_list_of_connections(segmented_grid::Grid; verbose::Bool=true)
    # initialize empty list of edges
    connections = SegmentConnection[]

    # loop over faces of unit cell
    for i = 1:segmented_grid.n_pts[1], j = 1:segmented_grid.n_pts[2]
        _note_connection!(segmented_grid.data[i, j, end], segmented_grid.data[i, j, 1], 
            connections, (0, 0, 1), verbose)
    end

    for j = 1:segmented_grid.n_pts[2], k = 1:segmented_grid.n_pts[3]
        _note_connection!(segmented_grid.data[end, j, k], segmented_grid.data[1, j, k],
            connections, (1, 0, 0), verbose)
    end
    
    for i = 1:segmented_grid.n_pts[1], k = 1:segmented_grid.n_pts[3]
        _note_connection!(segmented_grid.data[i, end, k], segmented_grid.data[i, 1, k],
            connections, (0, 1, 0), verbose)
    end
    
    return connections
end

function _translate_into_graph(segmented_grid::Grid, connections::Array{SegmentConnection, 1})
    nb_segments = _count_segments(segmented_grid)
        
    # directed graph
    graph = DiGraph(nb_segments)

    # so each edge REALLY has some data associated with it: which unit cell face is traversed
    #  when traveling from the source segment to the destination segment.
    #  goal: include this info about an edge indirectly by having an artificial vertex 
    #  indicating the direction. so all vertices have a direction. if the vertex is actually
    #  for representing a segment, it doens't have a direction, so we assign it a 
    #  direction `nothing`.
    vertex_to_direction = Dict{Int, Union{Nothing, Tuple{Int, Int, Int}}}()
    for v in vertices(graph)
        vertex_to_direction[v] = nothing
    end
    
    # loop over edges, add edges to graph but with intermediate dummy vertex corresponding
    #  to direction
    for segment_connection in connections
        # add vertex for the direction
        add_vertex!(graph)
        # id of vertex indicating connection
        connection_vertex_id = length(vertices(graph))
        # record direction of unit cell traversal in this edge
        vertex_to_direction[connection_vertex_id] = segment_connection.direction
        # add edge of segment connectivity, first going through direction segment
        #   src segment --> connection_vertex_id --> dst segment
        add_edge!(graph, segment_connection.src, connection_vertex_id)
        add_edge!(graph, connection_vertex_id, segment_connection.dst)
    end

    return graph, vertex_to_direction
end

function _classify_segments(segmented_grid::Grid, graph::SimpleDiGraph{Int64}, 
                            vertex_to_direction::Dict{Int, Union{Nothing, Tuple{Int, Int, Int}}},
                            verbose::Bool=true)
    nb_segments = _count_segments(segmented_grid)

    # 0 = inaccessible, 1 = accessible. start out with assuming inaccessible.
    segment_classifiction = zeros(Int, nb_segments)
    
    # look for simple cycles in the graph
    all_cycles = simplecycles(graph)
    if verbose
        @printf("\tFound %d simple cycles in segment connectivity graph.\n", length(all_cycles))
    end
    
    # all edges involved in a cycle are connected together and are accessible channels
    for cyc in all_cycles
        @assert vertex_to_direction[cyc[1]] == nothing "shld start with segment."
        @assert vertex_to_direction[cyc[end]] != nothing "shld end with direction vertex."
        # travel through cycle, keep track of which unit cell we're in
        #  start in home unit cell
        unit_cell = (0, 0, 0)
        # loop over segments in the cycle (ordered)
        for s in cyc
            # if an intermediate unit cell boundary vertex, then keep track of unit cell
            #  we're in by looking at boundary we traversed.
            if vertex_to_direction[s] != nothing
                unit_cell = unit_cell .+ vertex_to_direction[s]
            end
        end

        # if we traversed the cycle and ended up in a different unit cell, these all
        #   must be accessible pockets!
        if unit_cell != (0, 0, 0)
            if verbose
                println("\t...found a cycle of accessible segments!")
            end
            for s in cyc
                if vertex_to_direction[s] == nothing # i.e. if this is not a direciton vertex
                    segment_classifiction[s] = 1
                end
            end
        end
    end

    # now could have something like 1 --> 2 ---> 4 ---> 2. then 1 is also accessible
    # loop over all segments
    for s = 1:nb_segments
        # if this segment was involved in a cycle that resulted in ending up in a different
        #  unit cell
        if segment_classifiction[s] == 1
            # look for edges connected to this segment; these must be accessible
            for t = 1:nb_segments
                if Edge(t, s) in edges(graph)
                    segment_classifiction[t] = 1
                end
            end
        end
    end

    return segment_classifiction
end

function _assign_inaccessible_pockets_minus_one!(segmented_grid::Grid, 
            segment_classifiction::Array{Int, 1}; verbose::Bool=true)
    for s = 1:length(segment_classifiction)
        if segment_classifiction[s] == 1
            if verbose
                @printf("Segment %s classified as accessible channel.\n", s)
            end
        else
            if verbose
                @printf("Segment %s classified as inaccessible pocket. Overwriting.\n", s)
            end
            segmented_grid.data[segmented_grid.data .== s] .= -1
        end
    end
end

"""
    accessibility_grid, nb_segments_blocked, porosity = compute_accessibility_grid(framework, 
    probe_molecule, ljforcefield; n_pts=(20, 20, 20), energy_tol=10.0, energy_unit=:kJ_mol,
    verbose=true, write_b4_after_grids=false, block_inaccessible_pockets=true)

Overlay a grid of points about the unit cell. Compute the potential energy of a probe
molecule at each point. If the potential energy is less than `energy_tol`, the grid point
is declared as accessible to an adsorbate; otherwise inaccessible.

If `block_pockets` is true:
Then perform a flood fill algorithm to label disparate (unconnected) segments in the grid.

Then build a graph whose vertices are the unconnected segments in the flood-filled grid and
whose edges are the connections between the segments across the periodic boundary.

Then find any simple cycles in the grid. Any vertex that is involved in a simple cycle is
considered accessible since a molecule can travel from that segment in the home unit cell
to the same segment but in a different unit cell. If any vertex is not involved in a cycle,
the segment is declared as inaccessible and all grid points in this segment are re-labeled
as inaccessible.

Returns `accessibility_grid::Grid{Bool}` and `nb_segments_blocked`, the latter the number
of segments that were blocked because they were determined to be inaccessible.

# Arguments
* `framework::Framework`: the crystal for which we seek to compute an accessibility grid.
* `probe_molecule::Molecule` a molecule serving as a probe to determine whether a given 
point can be occupied and accessed.
* `LJForceField::LJForceField`: the force field used to compute the potential energy of 
the probe molecule
* `n_pts::Tuple{Int, Int, Int}`: number of grid points in a, b, c directions
* `energy_tol::Float64`: if the computed potential energy is less than this, we declare the
grid point to be occupiable. Also this is the energy barrier beyond which we assume the
probe adsorbate cannot pass. Units given by `energy_units` argument
* `energy_units::Symbol`: units of energy (`:kJ_mol` or `:K`) to be used in determining
threshold for occupiability and whether molecule can percolate over barrier in channel.
(see `energy_tol`)
* `write_b4_after_grids::Bool`: write a .cube file of occupiability for visualization both
before and after flood fill/blocking inaccessible pockets
"""
function compute_accessibility_grid(framework::Framework, probe::Molecule, forcefield::LJForceField;
                                    n_pts::Tuple{Int, Int, Int}=(20, 20, 20), 
                                    energy_tol::Float64=10.0,  energy_units::Symbol=:kJ_mol,
                                    verbose::Bool=true, write_b4_after_grids::Bool=true, 
                                    block_inaccessible_pockets::Bool=true)
    if verbose
        printstyled(@sprintf("Computing accessibility grid of %s using %f %s potential energy tol and %s probe...\n",
            framework.name, energy_tol, energy_units, probe.species), color=:green)
    end
    
    # write potential energy grid
    grid = energy_grid(framework, probe, forcefield, n_pts=n_pts, verbose=verbose, 
                       units=energy_units)

    if ! block_inaccessible_pockets
        accessibility_grid = Grid{Bool}(grid.box, grid.n_pts, grid.data .< energy_tol, 
                                        :accessibility, grid.origin)

        porosity = sum(accessibility_grid.data) / length(accessibility_grid.data)

        return accessibility_grid, 0, porosity
    end
    
    # flood fill and label segments
    segmented_grid = _segment_grid(grid, energy_tol, verbose)

    if write_b4_after_grids
        _segmented_grid = deepcopy(segmented_grid)
        _segmented_grid.data[segmented_grid.data .!= -1] .= 1
        gridfilename = @sprintf("%s_in_%s_%s_b4_pocket_blocking.cube", probe.species,
            replace(replace(framework.name, ".cif" => ""), ".cssr" => ""),
            replace(forcefield.name, ".csv" => ""))
        write_cube(_segmented_grid, gridfilename)
    end
    
    # get list of edges describing connectivity of segments across unit cell boundaries
    connections = _build_list_of_connections(segmented_grid, verbose=verbose)

    # translate this into a directed LightGraph. Note that these include directions as an
    #  artificial vertex to keep track of which unit cell boundary is traversed.
    graph, vertex_to_direction = _translate_into_graph(segmented_grid, connections)

    # get classifications of the segments
    segment_classifications = _classify_segments(segmented_grid, graph, vertex_to_direction, verbose)

    # assign inaccessible pockets minus one if cycle not found in graph
    _assign_inaccessible_pockets_minus_one!(segmented_grid, segment_classifications, verbose=verbose)
    
    # -1 for not accessible, 1 for accessible
    segmented_grid.data[segmented_grid.data .!= -1] .= 1
    
    if write_b4_after_grids
        gridfilename = @sprintf("%s_in_%s_%s_after_pocket_blocking.cube", probe.species,
            replace(replace(framework.name, ".cif" => ""), ".cssr" => ""),
            replace(forcefield.name, ".csv" => ""))
        write_cube(segmented_grid, gridfilename)
    end
    
    # print warning if there is no use in running a simulation since all pockets are inaccessible
    if all(segmented_grid.data .== -1)
        @warn @sprintf("%s cannot enter the pores of %s with %f K energy tolerance.",
            probe.species, framework.name, energy_tol)
    end
    
    nb_segments_blocked = sum(segment_classifications .== 0)

    # instead of returning grid of int's return a boolean grid for storage efficiency
    accessibility_grid = Grid{Bool}(segmented_grid.box, segmented_grid.n_pts,
        segmented_grid.data .== 1, :accessibility, segmented_grid.origin)

    # compute porosity before and after blocking
    porosity = Dict(:b4_blocking    => sum(grid.data .< energy_tol) / length(grid.data),
                    :after_blocking => sum(accessibility_grid.data) / length(accessibility_grid.data))

    if nb_segments_blocked != 0
        printstyled(@sprintf("\t%d pockets in %s were found to be inaccessible to %s and blocked.\n",
            sum(segment_classifications .== 0), framework.name, probe.species), color=:yellow)
        @printf("\tPorosity of %s b4 pocket blocking is %f\n", framework.name, porosity[:b4_blocking])
        @printf("\tPorosity of %s after pocket blocking is %f\n", framework.name, porosity[:after_blocking])
        @printf("\tpotential energy barrier used to determine accessibility: %f %s\n",
            energy_tol, energy_units)
    end

    return accessibility_grid, nb_segments_blocked, porosity
end

# return ID that is the nearest neighbor.
function _arg_nearest_neighbor(n_pts::Tuple{Int, Int, Int}, xf::Array{Float64, 1})
    return 1 .+ round.(Int, xf .* (n_pts .- 1))
end

function _apply_pbc_to_index!(id::Array{Int, 1}, n_pts::Tuple{Int, Int, Int})
    for k = 1:3
        if id[k] == 0
            id[k] = n_pts[k]
        end
        if id[k] == n_pts[k] + 1
            id[k] = 1
        end
    end
    return nothing
end

"""
    accessible(accessibility_grid, xf)
    accessible(accessibility_grid, xf, repfactors)

Using `accessibility_grid`, determine if fractional coordinate `xf` (relative to 
`accessibility_grid.box` is accessible or not. Here, we search for the nearest grid point.
We then look at the accessibility of this nearest grid point and all surroudning 9 other 
grid points. The point `xf` is declared inaccessible if and only if all 10 of these grid
points are inaccessible. We take this approach because, if the grid is coarse, we can allow
energy computations to automatically determine accessibility at the boundary of 
accessibility e.g. during a molecular simulation where inaccessible pockets are blocked.

If a tuple of replication factors are also passed, it is assumed that the passed `xf` is 
relative to a replicated `accessibility_grid.box` so that `xf` is scaled by these rep. 
factors. So `xf = [0.5, 0.25, 0.1]` with `repfactors=(1, 2, 4)` actually is, relative to
`accessibility_grid.box`, fractional coordinate `[0.5, 0.5, 0.4]`.
"""
function accessible(accessibility_grid::Grid{Bool}, xf::Array{Float64, 1})
    # find ID of nearest neighbor
    id_nearest_neighbor = _arg_nearest_neighbor(accessibility_grid.n_pts, xf)
    # the point is inaccessible if and only if ALL surrounding points are inaccessible;
    #   if at the boundary, let the energy be computed and automatically address
    #   accessibility during a simulation
    for i = -1:1, j = -1:1, k = -1:1
        id_neighbor = id_nearest_neighbor .+ [i, j, k]
        _apply_pbc_to_index!(id_neighbor, accessibility_grid.n_pts)
        # again, if ANY surrounding point is accessible, let the energy computation go on
        if accessibility_grid.data[id_neighbor...]
            return true
        end
    end
    # if we made it to here, there is no surrounding point that is accessible.
    return false
end

# when rep factors are used in the simulation so fractional coord system in grid does not match
#   that in the simulation. so `xf` is repect to simulation box which is replicated box in
#   the accessibility grid.
function accessible(accessibility_grid::Grid{Bool}, xf::Array{Float64, 1}, repfactors::Tuple{Int, Int, Int})
    return accessible(accessibility_grid, mod.(xf .* repfactors, 1.0))
end
