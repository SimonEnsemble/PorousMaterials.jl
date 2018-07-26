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

    if ! contains(filename, ".cube")
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
    if ! contains(filename, ".cube")
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
    box = construct_box(f_to_c)

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
- `verbose::Bool=true`: print some information.

# Returns
- `grid::Grid`: A grid data structure containing the potential energy of the system
"""
function energy_grid(framework::Framework, molecule::Molecule, ljforcefield::LJForceField;
                     n_pts::Tuple{Int, Int, Int}=(50,50,50), n_rotations::Int=1000, temperature::Float64=NaN, units::Symbol=:kJ_mol, verbose::Bool=true)
    if ! (units in [:kJ_mol, :K])
        error("Pass :kJ_mol or :K for units of kJ/mol or K, respectively.")
    end

    const rotations_required = rotatable(molecule)
    const charged_system = (length(framework.charges) > 1) && (length(molecule.charges) > 1)
    if rotations_required & isnan(temperature)
        error("Must pass temperature (K) for Boltzmann weighted rotations.\n")
    end
    
    eparams, kvecs, eikar, eikbr, eikcr = setup_Ewald_sum(sqrt(ljforcefield.cutoffradius_squared),
                                                          framework.box,
                                                          verbose=verbose & charged_system)


    repfactors = replication_factors(framework.box, ljforcefield)
    framework = replicate(framework, repfactors)
    
    # grid of voxel centers (each axis at least).
    grid_pts = [collect(linspace(0.0, 1.0, n_pts[i])) for i = 1:3]

    grid = Grid(framework.box, n_pts, zeros(Float64, n_pts...), units, [0.0, 0.0, 0.0])

    if verbose
        @printf("Computing energy grid of %s in %s\n", molecule.species, framework.name)
        @printf("\tRegular grid (in fractional space) of %d by %d by %d points superimposed over the unit cell.\n", n_pts...)
        if rotations_required
            @printf("\t%d molecule rotations per grid point with temperature %f K.\n", n_rotations, temperature)
        end
    end

	for (i, xf) in enumerate(grid_pts[1]), (j, yf) in enumerate(grid_pts[2]), (k, zf) in enumerate(grid_pts[3])
        translate_to!(molecule, [xf, yf, zf])
        if ! rotations_required
            ensemble_average_energy = vdw_energy(framework, molecule, ljforcefield)
        else
            boltzmann_factor_sum = 0.0
            for r = 1:n_rotations
                rotate!(molecule, framework.box)

                energy = PotentialEnergy(0.0, 0.0)
                energy.vdw = vdw_energy(framework, molecule, ljforcefield)
                if charged_system
                    energy.coulomb = electrostatic_potential_energy(framework, molecule, eparams,
                                                                    kvecs, eikar, eikbr, eikcr)
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
