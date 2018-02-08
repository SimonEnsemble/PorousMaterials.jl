using NPZ

"""
snapshot = Snapshot(frame::Framework, molecule::Molecule, ljforcefield::LennardJonesForceField, dimensions::Array{Float64}, startpoint::Array{Float64}, repfactors::Tuple{Int64, Int64, Int64}, roughness::Float64)

Forms a Snapshot object to take a snapshot for with dimensions `dimensions`, starting at `startpoint`. 
The mesh is calculated using `roughness` to measure distance between points in all dimensions.

# Arguments
- `frame::Framework`: Crystal structure
- `molecule::Molecule`: adsorbate (includes position/orintation)
- `ljforcefield::LennardJonesForceField`: Lennard Jones force field
- `dimensions::Array{Float64}`: dimensions of the 3D snapshot being taken (Units: Angstrom)
- `startpoint::Array{Float64}`: Startpoint for the snapshot
- `repfactors::Tuple{Int64, Int64, Int64}`: Replication factors used to simulate the simulation box
- `roughness::Float64`: Distance between points in all dimensions. Used to calculate the meshgrid
"""
type Snapshot
	frame::Framework
	molecule::Molecule
	ljforcefield::LennardJonesForceField
	dimensions::Array{Float64}
	startpoint::Array{Float64}
	repfactors::Tuple{Int64, Int64, Int64}
	roughness::Float64
end

"""
	{occupancy, energy} = snap(snapshot::Snapshot, output_type::AbstractString="occupancy")

Returns a occupancy matrix or an energy matrix depending on which output type you choose. Supported output types are:
	"`occupancy`": Returns a bitArray corresponding to the van Der Waals energy in the Framework. The bitArray will be a three dimensional matrix with dimensions described in `Snapshot.dimensions` in cartesian coordinates. The `Snapshot.startpoint` will determine where the iteration starts at. If the van Der Waals interaction between the adsorbate `molecule` and the Framework `frame` greater than zero, that corresponds to a true (occupied) value in the bitArray. If the interaction is less than zero (not occupied) it will correspond to a false value in the bitArray.
	"`energy`": Returns a three dimensional Float64 Matrix with the van der Waals energy at every point in the meshgrid (the same as occupancy except there's no check to see if the energy is greater than or less than zero.
	"`both`": Returns both matrices (ex. (occupancy, energy) = snap(snapshot, "both")
"""
function snap(snapshot::Snapshot, output_type::AbstractString="occupancy")
	if output_type=="occupancy"
		return occupancy_snap(snapshot)
	elseif output_type=="energy"
		return energy_snap(snapshot)
	elseif output_type=="both"
		return both_snap(snapshot)
	else
		error("Output_type not supported. Only supported output_types are 'energy' and 'occupancy'")
	end
end

function energy_snap(snapshot::Snapshot)
	x_cnt, y_cnt, z_cnt = floor.(Int64, snapshot.dimensions ./ snapshot.roughness)
	energy = zeros(x_cnt, y_cnt, z_cnt)

	# Cartesian range in x,y and z-dimensions
	cart_range_x = snapshot.startpoint[1] + linspace(0,snapshot.dimensions[1],x_cnt)
	cart_range_y = snapshot.startpoint[2] + linspace(0,snapshot.dimensions[2],y_cnt)
	cart_range_z = snapshot.startpoint[3] + linspace(0,snapshot.dimensions[3],z_cnt)

	for (i,x) in enumerate(cart_range_x), (j,y) in enumerate(cart_range_y), (k,z) in enumerate(cart_range_z)
		snapshot.molecule.x = ([x,y,z])[:,:]
		energy[i,j,k] = vdw_energy(snapshot.frame, snapshot.molecule, snapshot.ljforcefield, snapshot.repfactors)
	end

	return energy
end

function occupancy_snap(snapshot::Snapshot)
	x_cnt, y_cnt, z_cnt = floor.(Int64, snapshot.dimensions ./ snapshot.roughness)
	occupancy = trues(x_cnt, y_cnt, z_cnt)

	# Cartesian range in x,y and z-dimensions
	cart_range_x = snapshot.startpoint[1] + linspace(0,snapshot.dimensions[1],x_cnt)
	cart_range_y = snapshot.startpoint[2] + linspace(0,snapshot.dimensions[2],y_cnt)
	cart_range_z = snapshot.startpoint[3] + linspace(0,snapshot.dimensions[3],z_cnt)

	for (i,x) in enumerate(cart_range_x), (j,y) in enumerate(cart_range_y), (k,z) in enumerate(cart_range_z)
		snapshot.molecule.x = ([x,y,z])[:,:]
		if vdw_energy(snapshot.frame, snapshot.molecule, snapshot.ljforcefield, snapshot.repfactors) < 0
			occupancy[i,j,k] = false
		end
	end

	return occupancy
end

function both_snap(snapshot::Snapshot)
	x_cnt, y_cnt, z_cnt = floor.(Int64, snapshot.dimensions ./ snapshot.roughness)
	occupancy = trues(x_cnt, y_cnt, z_cnt)
	energy = zeros(x_cnt, y_cnt, z_cnt)

	# Cartesian range in x,y and z-dimensions
	cart_range_x = snapshot.startpoint[1] + linspace(0,snapshot.dimensions[1],x_cnt)
	cart_range_y = snapshot.startpoint[2] + linspace(0,snapshot.dimensions[2],y_cnt)
	cart_range_z = snapshot.startpoint[3] + linspace(0,snapshot.dimensions[3],z_cnt)

	for (i,x) in enumerate(cart_range_x), (j,y) in enumerate(cart_range_y), (k,z) in enumerate(cart_range_z)
		snapshot.molecule.x = ([x,y,z])[:,:]
		energy[i,j,k] = vdw_energy(snapshot.frame, snapshot.molecule, snapshot.ljforcefield, snapshot.repfactors) 
		if energy[i,j,k] < 0
			occupancy[i,j,k] = false
		end
	end

	return (occupancy, energy)
end

"""
	write_snapshot_to_vtk(snapshot::Snapshot, filename::Union{Void, AbstractString}=nothing)

Writes the dimensions of the snapshot to a vtk for visualization.
If no filename is specified this defaults to a `frame.name + .vtk` filename
"""
function write_snapshot_to_vtk(snapshot::Snapshot, filename::Union{Void, AbstractString}=nothing)
	if filename == nothing
		filename = split(snapshot.frame.name, ".")[1] * ".vtk"
	end

	endx = [snapshot.startpoint[1], snapshot.startpoint[1] + snapshot.dimensions[1]]
	endy = [snapshot.startpoint[2], snapshot.startpoint[2] + snapshot.dimensions[2]]
	endz = [snapshot.startpoint[3], snapshot.startpoint[3] + snapshot.dimensions[3]]

	vtk_file = open(filename, "w")

	@printf(vtk_file, "# vtk DataFile Version 2.0\nunit cell boundary\n
					   ASCII\nDATASET POLYDATA\nPOINTS 8 double\n")


    # write points on boundary of unit cell
    for i = 1:2
        for j = 1:2
            for k = 1:2
				cornerpoint = [endx[i], endy[j], endz[k]] # fractional coordinates of corner
                @printf(vtk_file, "%.3f %.3f %.3f\n",
                        cornerpoint[1], cornerpoint[2], cornerpoint[3])
            end
        end
    end

    # define connections
    @printf(vtk_file, "LINES 12 36\n2 0 1\n2 0 2\n2 1 3\n2 2 3\n2 4 5\n2 4 6\n2 5 7\n2 6 7\n2 0 4\n2 1 5\n2 2 6\n2 3 7\n")
    close(vtk_file)
    println("See ", filename)
end


"""
	write_to_npy(matrix, "my_file.npy")

writes the information in `matrix` to a .npy file. The .npy file is able to be read from Python 3 using the numpy package
"""
function write_to_npy(matrix::Union{Array{Float64,3}, BitArray{3}}, filename::Union{Void, AbstractString}=nothing)
	if filename == nothing
		filename = "numpy_data.npy"
	end

	if typeof(matrix) == BitArray{3}
		npzwrite("data/db/" * filename, convert.(UInt8, matrix))
		@printf("BitArray converted to UInt8 and saved in %s\n","data/db/" * filename)
	else
		npzwrite("data/db/" * filename, matrix)
		@printf("Array{Float64,3} saved in %s\n","data/db/" * filename)
	end
end

import Base.print
function print(io::IO, snapshot::Snapshot)
	 @printf(io, "Using frame = %s\n",snapshot.frame.name)
	 @printf(io, "Corners of snapshot:\n")
	 @printf(io, "\t%.3f\t%.3f\t%.3f\n", snapshot.startpoint[1], snapshot.startpoint[2], snapshot.startpoint[3])
	 @printf(io, "\t%.3f\t%.3f\t%.3f\n", snapshot.dimensions[1] + snapshot.startpoint[1], snapshot.dimensions[2] + snapshot.startpoint[2], snapshot.dimensions[3] + snapshot.startpoint[3])
	 @printf(io, "Roughness = %.3f", snapshot.roughness)
 end

 import Base.show
 function show(io::IO, snapshot::Snapshot)
	 print(io, snapshot)
 end
