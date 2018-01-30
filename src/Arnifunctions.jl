"""
	EleProps = readproperties("filepath")

Reads a .csv file with the following header: Element name, Molecular mass, Atomic radius, Ionic radius and returns a dictionary with an element name pointing at an array of values,
EleProps[Element] => [Molecular mass(amu), Atomic radius(Angstrom), Ionic radius(Angstrom)]
"""
function readproperties(filename::AbstractString)
	eleprops = Dict{AbstractString,Array{Float64}}()
	f = open(filename,"r")
	lines = readlines(f)

	for (i,line) in enumerate(lines)
		if (i>1)
			str = split(line,",")
			temp = zeros(3)
			for (k,val) in enumerate(str[2:4])
				if val!="NA"
					temp[k] = parse(Float64,val)
				else
					temp[k] = 0.0
				end
			end
			eleprops[str[1]] = temp
		end
	end
	close(f)
	return eleprops
end # end readproperties


"""
	rotate(rotMatrix::Array{Float64,2},θ::Float64,ϕ::Float64,ψ::Float64)

Changes the rotation matrix `rotMatrix` to a rotation defined by the angles θ,ϕ,ψ (see http://mathworld.wolfram.com/EulerAngles.html)
"""
#TODO See if rotate! changes the inputting argument rather than outputting a new one
function rotate(rotMatrix::Array{Float64,2},θ::Float64,ϕ::Float64,ψ::Float64)
	B = [cos(ψ) sin(ψ) 0; -sin(ψ) cos(ψ) 0; 0 0 1]
	C = [1 0 0; 0 cos(θ) sin(θ); 0 -sin(θ) cos(θ)]
	D = [cos(ϕ) sin(ϕ) 0; -sin(ϕ) cos(ϕ) 0; 0 0 1]
	rotMatrix = B*C*D
end


function exploreframe(frame::Framework, molecule::Molecule, ljforcefield::LennardJonesForceField, repfactors::Tuple{Int64, Int64, Int64}; roughness = 0.01)
	coordmatrix = Dict{AbstractString,Array{Float64,2}}()
	a_cnt, b_cnt, c_cnt = [repfactors[1], repfactors[2], repfactors[3]] / roughness
	a_cnt, b_cnt, c_cnt = [convert(Int64, floor(a_cnt)), convert(Int64, floor(b_cnt)), convert(Int64, floor(c_cnt))]
	cnt = (a_cnt, b_cnt, c_cnt)
	@printf("%f, %f, %f\n",a_cnt, b_cnt, c_cnt)
	EnergyMatrix = zeros(a_cnt,b_cnt,c_cnt)
	frac_range_x = range(0.,roughness,a_cnt)
	frac_range_y = range(0.,roughness,b_cnt)
	frac_range_z = range(0.,roughness,c_cnt)
	@printf("%f, %f, %f\n",length(frac_range_x), length(frac_range_y), length(frac_range_z))

	for (i,xf) in enumerate(frac_range_x), (j,yf) in enumerate(frac_range_y), (k,zf) in enumerate(frac_range_z)
		molecule.x = (frame.box.f_to_c*[xf,yf,zf])[:,:]
		EnergyMatrix[i,j,k] = vdw_energy(frame, molecule, ljforcefield, repfactors)	
		tempstr = @sprintf("%d-%d-%d",i,j,k)
		coordmatrix[tempstr] = molecule.x
	end
	return (EnergyMatrix,coordmatrix, cnt)
end
"""
occupancy = takesnapshot(frame::Framework, molecule::Molecule, ljforcefield::LennardJonesForceField, dimensions::Array{Float64}, repfactors::Tuple{Int64, Int64, Int64}; startpoint::Array{Float64}=[0.0,0.0,0.0], roughness = 0.1) 

Forms a bitArray corresponding to the van Der Waals energy in the Framework. bitArray will be a three dimensional matrix with dimensions described in `dimensions` in cartesian coordinates. The startpoint will determine where the iteration starts at. If the van Der Waals interaction between the adsorbate `molecule` and the Framework `frame` is positive, that corresponds to a true (occupied) value in the bitArray. If the interaction is negative (not occupied) it will correspond to a false value in the bitArray.

`mesh` is an array describing how fine we want our grid (how big of a matrix). The three numbers in the array correspond to the x-, y- and z-direction respectively. 
"""
function takesnapshot(frame::Framework, molecule::Molecule, ljforcefield::LennardJonesForceField, dimensions::Array{Float64}, repfactors::Tuple{Int64, Int64, Int64}; startpoint::Array{Float64}=[0.0,0.0,0.0], roughness = 0.1) 
	# If vdw_energy > 0 at a certain point, it will correspond to a 'true' in the bitArray. If vdw_energy < 0 (void space) it will correspond to a 'false' in the bitArray.	
	x_cnt, y_cnt, z_cnt = floor.(Int64, dimensions ./ roughness)
	@printf("%f, %f, %f\n",x_cnt,y_cnt,z_cnt)
	@printf("%s, %s, %s\n",typeof(x_cnt),typeof(y_cnt),typeof(z_cnt))
	occupancy = trues(x_cnt, y_cnt, z_cnt)

	# Cartesian range in x,y and z-dimensions
	cart_range_x = startpoint[1] + linspace(0,dimensions[1],x_cnt)
	cart_range_y = startpoint[2] + linspace(0,dimensions[2],y_cnt)
	cart_range_z = startpoint[3] + linspace(0,dimensions[3],z_cnt)
	endpoints = [cart_range_x[1] cart_range_x[end]; cart_range_y[1] cart_range_y[end]; cart_range_z[1] cart_range_z[end]]

	for (i,x) in enumerate(cart_range_x), (j,y) in enumerate(cart_range_y), (k,z) in enumerate(cart_range_z)
		molecule.x = ([x,y,z])[:,:]
		if vdw_energy(frame, molecule, ljforcefield, repfactors) < 0
			occupancy[i,j,k] = false
		end
	end
	return occupancy
end
