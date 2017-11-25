"""
    Grid(box, nb_grid_pts, data, end_included)

Data structure for a regular [equal spacing between points] grid and associated data
with each grid point/voxel.
# TODO: define if data is value at center of voxel or at corner...

# Arguments
- `box::Box`: describes Bravais lattice onto which the grid is superimposed.
- `nb_grid_pts::Tuple{Int, Int, Int}`: number of grid points in x, y, z directions.
- `data::Array{T, 3}`: three dimensional array conaining data associated with each grid point.
- `units::AbstractString`: the units associated with each data point.
"""
struct Grid{T}
    box::Box
    nb_grid_pts::Tuple{Int, Int, Int}
    data::Array{T, 3}
    units::AbstractString
end

"""
    write_to_cube(grid::Grid, filename::AbstractString)

Write volume data to a .cube file format. This format is described here:
http://paulbourke.net/dataformats/cube/
The origin is assumed to be (0, 0, 0).
The atoms of the unit cell are not printed in the .cube. Instead, use .xyz files.

# Arguments
- `grid::Grid`: grid with associated volume data (see Grid struct)
- `filename::AbstractString`: name of .cube file to which we write the grid; this is relative to `PorousMaterials.PATH_TO_DATA`/grids/.
"""
function write_to_cube(grid::Grid, filename::AbstractString)
    cubefile = open(PATH_TO_DATA * "grids/" * filename, "w")
    
    @printf(cubefile, "Units of data: %s\nLoop order: x, y, z\n", grid.data)
    # the integer refers to 0 atoms (just use .xyz to visualize atoms)
    # the next three floats correspond to the origin, assumed to be (0,0,0)
    @printf(cubefile, "%d %f %f %f\n" , 0, 0.0, 0.0, 0.0)

    # these are the vectors that form the voxels of each grid parallelogram.
    for k = 1:3
        # TODO re-evaluate this. depending on how you define grid, may be
        #   / (grid.nb_grid_pts[k] - 1) e.g. if you 
        voxel_vector = box.f_to_c[:, k] / grid.nb_grid_pts[k]
        @printf(cubefile, "%d %f %f %f\n" , grid.nb_grid_pts[k],
            voxel_vector[1], voxel_vector[2], voxel_vector[3])
    end

    for i = 1:grid.nx
        for j = 1:grid.ny
            for k = 1:grid.nz
                @printf(cubefile, "%e", grid.data[i, j, k])
                if (k % 6) == 0
                    @printf(cubefile, "\n")
                end
            end # loop over z points
            @printf(cubefile, "\n")
        end # loop over y points
    end # loop over x points
    close(cubefile)
    println("See ", PATH_TO_DATA * "grids/" * filename)
    return 
end

# TODO: function to read .cube and return grid.
 # function read_cube(filename::AbstractString)
 #     cubefile = open(filename)
 #     
 #     grid = Grid(box, nb_grid_pts, data_matrix, units)
 #     return grid
 # end
