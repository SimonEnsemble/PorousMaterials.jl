"""
Data structure for a regular [equal spacing between points in each coordinate] grid of voxels and associated data with each voxel.
A unit cell (`Box`) is partitioned into voxels.

# Arguments
- `box::Box`: describes Bravais lattice which is partitioned into voxels in fractional coordinates
- `nb_voxels::Tuple{Int, Int, Int}`: number of grid points in x, y, z directions.
- `data::Array{T, 3}`: three dimensional array conaining data associated with each grid point.
- `units::Symbol`: the units associated with each data point.
"""
struct Grid{T}
    box::Box
    nb_voxels::Tuple{Int64, Int64, Int64}
    data::Array{T, 3}
    units::Symbol
end

function Base.show(io::IO, grid::Grid)
    @printf(io, "Grid of %d by %d by %d voxels comprising a unit cell and associated data.\n", grid.nb_voxels...)
    @printf(io, "\tunits of data attribute: %s\n", grid.units)
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
    
    @printf(cubefile, "Units of data: %s\nLoop order: x, y, z\n", grid.units)
    # the integer refers to 0 atoms (just use .xyz to visualize atoms)
    # the next three floats correspond to the origin, assumed to be (0,0,0)
    @printf(cubefile, "%d %f %f %f\n" , 0, 0.0, 0.0, 0.0)

    # these are the vectors that form the voxels of each grid parallelogram.
    for k = 1:3
        # TODO re-evaluate this. depending on how you define grid, may be
        #   / (grid.nb_grid_pts[k] - 1) e.g. if you 
        voxel_vector = grid.box.f_to_c[:, k] / grid.nb_grid_pts[k]
        @printf(cubefile, "%d %f %f %f\n" , grid.nb_grid_pts[k],
            voxel_vector[1], voxel_vector[2], voxel_vector[3])
    end

    for i = 1:grid.nb_grid_pts[1]
        for j = 1:grid.nb_grid_pts[2]
            for k = 1:grid.nb_grid_pts[3]
                @printf(cubefile, "%e ", grid.data[i, j, k])
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

"""
	grid = energy_grid(framework::Framework, moleclule::Molecule, ljforcefield::LennardJonesForcefield; n_gridpts::Tuple{Int, Int, Int}=(50,50,50))

Partitions the unit cell of a crystal into a grid of voxels (in fractional coordinates), with `n_gridpts` dictating the number of voxels in the a, b, c directions.
Then, for each grid point, calculate the ensemble average potential energy of the molecule when its mass is centered at that point. The average is taken over Boltzmann-weighted rotations.

The ensemble average is a Boltzmann average over rotations:  - R T log ⟨e⁻ᵇᵁ⟩
"""
function energy_grid(framework::Framework, molecule::Molecule, ljforcefield::LennardJonesForceField;
                     n_voxels::Tuple{Int, Int, Int}=(50,50,50), n_rotations::Int=1000, temperature::Float64=NaN, units::Symbol=:kJ_mol, verbose::Bool=true)
    if ! (units in [:kJ_mol, :K])
        error("Pass :kJ_mol or :K for units of kJ/mol or K, respectively.")
    end

    # TODO electrostatics
    if length(molecule.charges) != 0
        error("Electrostatics not implemented yet.")
    end

    const rotations_required = ((length(molecule.ljspheres) > 1) | (length(molecule.charges) > 1))
    if rotations_required & isnan(temperature)
        error("Must pass temperature (K) for Boltzmann weighted rotations.\n")
    end

    const repfactors = replication_factors(framework.box, ljforcefield)
    
    # grid of voxel centers (each axis at least).
    voxel_centers = [collect(linspace(0.0, 1.0, n_voxels[i] + 1)) for i = 1:3]
    for i = 1:3
        dxf = voxel_centers[i][2] - voxel_centers[i][1]
        pop!(voxel_centers[i])
        voxel_centers[i] += dxf / 2.0
    end

    grid = Grid(framework.box, n_voxels, zeros(Float64, n_voxels...), units)

    if verbose
        @printf("Computing energy grid of %s in %s\n", molecule.species, framework.name)
        @printf("\tUnit cell partitioned into %d by %d by %d voxels.\n", n_voxels...)
        if rotations_required
            @printf("\t%d rotations per voxel.\n", n_rotations)
        end
    end

	for (i, xf) in enumerate(voxel_centers[1]), (j, yf) in enumerate(voxel_centers[2]), (k, zf) in enumerate(voxel_centers[3])
        translate_to!(molecule, framework.box.f_to_c * [xf, yf, zf])
        if ! rotations_required
            ensemble_average_energy = vdw_energy(framework, molecule, ljforcefield, repfactors)
        else
            boltzmann_factor_sum = 0.0
            for r = 1:n_rotations
                rotate!(molecule)
                energy = vdw_energy(framework, molecule, ljforcefield, repfactors)
                boltzmann_factor_sum += exp(-energy / temperature)
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
