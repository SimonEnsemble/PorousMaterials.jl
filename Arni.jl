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

# TODO: Note sure if you will actually use this function, but consider that computing
# TODO   center of mass for a periodic image is tricky. https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions
# TODO   is this function invariant to shifting the home unit cell?
"""
	c0 = centerofmass(frame,"~/example/properties.csv")

Uses `readproperties` to get a dictionary of element properties and uses that to calculate the center of mass of a supercell made from Framework and replication factors. (See readElementProps for more info on that function)
Calculates the center of mass according to r_cm = ∑r_i*m_i/m_tot in fractional coordinates but returns in cartesian coordinates
"""
function centerofmass(frame::Framework, filename::AbstractString, rep_factors::Array{Int64})
	repA = rep_factors[1]
	repB = rep_factors[2]
	repC = rep_factors[3]

	eleprops = readproperties(filename)
	rvec = zeros(3)
	mtot = 0.0
	for nA=1:repA, nB=1:repB, nC=1:repC
		for i=1:frame.n_atoms
			rvec += (frame.f_coords[:,i]+[nA-1, nB-1, nC-1])*eleprops[frame.atoms[i]][1]
			mtot += eleprops[frame.atoms[i]][1]
		end
	end
	return frame.f_to_C*(rvec/mtot)
end


"""
	rotate!(rotMatrix::Array{Float64,2},θ::Float64,ϕ::Float64,ψ::Float64)

Changes the rotation matrix `rotMatrix` to a rotation defined by the angles θ,ϕ,ψ (see http://mathworld.wolfram.com/EulerAngles.html)
"""
function rotate(rotMatrix::Array{Float64,2},θ::Float64,ϕ::Float64,ψ::Float64)
	B = [cos(ψ) sin(ψ) 0; -sin(ψ) cos(ψ) 0; 0 0 1]
	C = [1 0 0; 0 cos(θ) sin(θ); 0 -sin(θ) cos(θ)]
	D = [cos(ϕ) sin(ϕ) 0; -sin(ϕ) cos(ϕ) 0; 0 0 1]
	rotMatrix = B*C*D
end

function exploreframe(frame::Framework, molecule::Molecule, ljforcefield::LennardJonesForceField, repfactors::Array{Int64}; rotation=false, mesh=101)
	if rotation
		rotMatrix = zeros(3,3)
		EnergyMatrix = zeros(mesh,mesh,mesh,mesh)
		# How should I determine what rotations to use?
		rotcnt = 10
	else
		EnergyMatrix = zeros(mesh,mesh,mesh)
		rotcnt = 1
	end
	coordmatrix = Dict{AbstractString,Array{Float64,2}}()
	frac_range_x = linspace(0,repfactors[1],mesh)	
	frac_range_y = linspace(0,repfactors[2],mesh)	
	frac_range_z = linspace(0,repfactors[3],mesh)	

	for (i,xf) in enumerate(frac_range_x), (j,yf) in enumerate(frac_range_y), (k,zf) in enumerate(frac_range_z)
		molecule.pos = (frame.f_to_C*[xf,yf,zf])[:,:]
#		@printf("pos = [%f, %f, %f]\n",pos[1],pos[2],pos[3])
#		@printf("i = %d, j = %d, k = %d\n",i,j,k)
		EnergyMatrix[i,j,k] = vdw_energy(frame, molecule, ljforcefield, repfactors)	
		tempstr = @sprintf("%d-%d-%d",i,j,k)
		coordmatrix[tempstr] = molecule.pos
	end
	println(typeof(coordmatrix))
	return (EnergyMatrix,coordmatrix)
end
