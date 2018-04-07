# Molecules are composed of Lennard-Jones spheres and point charges.
struct LennardJonesSphere
    "atom name corresponding to force field"
    atom::Symbol
    "Cartesian coordinates (units: A)"
    x::Array{Float64, 1}
end

function Base.isapprox(ljs1::LennardJonesSphere, ljs2::LennardJonesSphere)
    return ((ljs1.atom == ljs2.atom) & isapprox(ljs1.x, ljs2.x))
end

"A point charge"
struct PointCharge
    "charge magnitude (electrons)"
    q::Float64
    "Cartesian coordinates (units: A)"
    x::Array{Float64, 1}
end

function Base.isapprox(c1::PointCharge, c2::PointCharge)
    return ((c1.q == c2.q) & isapprox(c1.x, c2.x))
end

"""
Data structure for a molecule/adsorbate.

# Attributes
- `species::Symbol`: Species of molecule, e.g. `CO2` or `ethane`
- `ljspheres::Array{LennardJonesSphere, 1}`: array of Lennard-Jones spheres comprising the molecule
- `charges::Array{PointCharges, 1}`: array of point charges comprising the molecule
- `center_of_mass::Array{Float64, 1}`: center of mass of the molecule
"""
struct Molecule
    species::Symbol
    ljspheres::Array{LennardJonesSphere, 1}
    charges::Array{PointCharge, 1}
    center_of_mass::Array{Float64, 1}
end

function Base.isapprox(m1::Molecule, m2::Molecule)
    species = m1.species == m2.species
    centers_of_mass = isapprox(m1.center_of_mass, m2.center_of_mass)
    ljs_counts = length(m1.ljspheres) == length(m2.ljspheres)
    pc_counts = length(m1.charges) == length(m2.charges)
    if (species & centers_of_mass & ljs_counts & pc_counts)
        ljspheres = all([isapprox(m1.ljspheres[i], m2.ljspheres[i]) for i = 1:length(m1.ljspheres)])
        charges = all([isapprox(m1.charges[i], m2.charges[i]) for i = 1:length(m1.charges)])
    else
        return false
    end
end

"""
    molecule = read_molecule_file("Xe")

Reads molecule files in the directory `PorousMaterials.PATH_TO_DATA * "/molecule/" * species * "/"`.
Center of mass assigned using atomic masses from `read_atomic_masses()`
"""
function read_molecule_file(species::String)
    if ! isdir(PATH_TO_DATA * "molecules/" * species)
        error(@sprintf("No directory created for %s in %s\n", species,
                       PATH_TO_DATA * "molecules/"))
    end

    # Read in atomic masses
    if ! isfile(PATH_TO_DATA * "atomicmasses.csv")
        error(@sprintf("No file 'atomicmasses.csv' exists. This file is needed to calculate the center of mass"))
    end
    atomic_masses = read_atomic_masses()
    
    # Read in Lennard Jones spheres
    ljspheresfilename = PATH_TO_DATA * "molecules/" * species * "/lennard_jones_spheres.csv"
    if ! isfile(ljspheresfilename)
        error(@sprintf("No file %s exists. Even if there are no Lennard Jones spheres in your molecule, include a .csv file with the proper headers but no rows.",
                       ljspheresfilename))
    end
    df_lj = CSV.read(ljspheresfilename)
    
    COM = [0.0, 0.0, 0.0]
    total_mass = 0.0
    ljspheres = LennardJonesSphere[]
    for row in eachrow(df_lj)
        x = [row[:x], row[:y], row[:z]]
        atom = Symbol(row[:atom])
        push!(ljspheres, LennardJonesSphere(atom, x))
        if ! (atom in keys(atomic_masses))
            error(@sprintf("Atomic mass of %s not found. See `read_atomic_masses()`\n", atom))
        end
        total_mass += atomic_masses[atom]
        COM += atomic_masses[atom] .* x
    end
    COM /= total_mass

    # Read in point charges
    chargesfilename = PATH_TO_DATA * "molecules/" * species * "/point_charges.csv"
    if ! isfile(chargesfilename)
        error(@sprintf("No file %s exists. Even if there are no point charges in your molecule, include a .csv file with the proper headers but no rows.",
                       chargesfilename))
    end
    df_c = CSV.read(chargesfilename)

    charges = PointCharge[]
    for row in eachrow(df_c)
        x = [row[:x], row[:y], row[:z]]
        push!(charges, PointCharge(row[:q], x))
    end

    # check for charge neutrality
    if length(charges) > 0
        if ! (sum([charges[i].q for i = 1:length(charges)]) ≈ 0.0)
            error(@sprintf("Molecule %s is not charge neutral!", species))
        end
    end
    
    return Molecule(Symbol(species), ljspheres, charges, COM)
end

"""
    translate!(molecule::Molecule, dx::Array{Float64, 1})

Translate a molecule by Cartesian vector `dx`. Adds `dx` to the coordinates of each `LennardJonesSphere` and `Charge` in the `Molecule` and shifts the center of mass as well.
This function does *not* account for periodic boundary conditions applied in the context of a `Box`.

# Arguments
- `molecule::Molecule`: the molecule to translate
- `dx::Array{Float64, 1}`: A 1-D array with 3 elements, corresponding to the Cartesian x-, y- and z- directions, that dictate by how much the molecule is perturbed in that coordinate.
"""
function translate_by!(molecule::Molecule, dx::Array{Float64, 1})
    for ljsphere in molecule.ljspheres
        ljsphere.x[:] += dx
    end

    for charge in molecule.charges
        charge.x[:] += dx
    end
    molecule.center_of_mass[:] += dx
end

"""
    translate_to!(molecule::Molecule, x::Array{Float64, 1})

Translate a molecule so that its center of mass is at Cartesian vector `x`. Charges the coordinates of each `LennardJonesSphere` and `Charge` in the `Molecule` and its center of mass as well.
This function does *not* account for periodic boundary conditions applied in the context of a `Box`.

# Arguments
- `molecule::Molecule`: the molecule to translate
- `x::Array{Float64, 1}`: A 1-D array with 3 elements, corresponding to the Cartesian x-, y- and z- directions, that dictates where its new center of mass will be.
"""
function translate_to!(molecule::Molecule, x::Array{Float64, 1})
    translate_by!(molecule, x - molecule.center_of_mass)
end

function Base.show(io::IO, molecule::Molecule)
    println(io, "Molecule species: ", molecule.species)
    println(io, "Center of mass: ", molecule.center_of_mass)
    if length(molecule.ljspheres) > 0
        print(io, "Lennard-Jones spheres: ")
        for ljsphere in molecule.ljspheres
            @printf(io, "\n\tatom = %s, x = [%.3f, %.3f, %.3f]", ljsphere.atom,
                    ljsphere.x[1], ljsphere.x[2], ljsphere.x[3])
        end
    end
    if length(molecule.charges) > 0
        print(io, "\nPoint charges: ")
        for charge in molecule.charges
            @printf(io, "\n\tq = %f, x = [%.3f, %.3f, %.3f]", charge.q, 
                    charge.x[1], charge.x[2], charge.x[3])
        end
    end
end

"""
    u = rand_point_on_unit_sphere()
    
Generate a unit vector with a random orientation.
"""
function rand_point_on_unit_sphere()
    u = [randn(), randn(), randn()]
    u_norm = norm(u)
    if u_norm < 1e-6 # avoid numerical error in division
        return rand_point_on_unit_sphere()
    end
    return u / u_norm
end

"""
    r = rotation_matrix(u::Array{Float64, 1}, θ::Float64)

Construct a matrix to perform a rotation of `θ` about the unit vector `u`.
"""
function rotation_matrix(u::Array{Float64, 1}, θ::Float64)
    @assert(norm(u) ≈ 1.0, "pass a unit vector!")
    c = cos(θ)
    s = sin(θ)
    r = [c + (1.0 - c) * u[1] ^ 2            (1.0 - c) * u[2] * u[1] - s * u[3]   (1.0 - c) * u[3] * u[1] + s * u[2];
         (1.0 - c) * u[1] * u[2] + s * u[3]  c + (1.0 - c) * u[2] ^ 2             (1.0 - c) * u[3] * u[2] - s * u[1];
         (1.0 - c) * u[1] * u[3] - s * u[2]  (1.0 - c) * u[2] * u[3] + s * u[1]   c + (1.0 - c) * u[3] ^ 2           ]
   return r 
end

"""
    rotate!(molecule::Molecule; assume_center_of_mass_at_origin::Bool=false) # molecule has rotated about its center of mass

Conduct a random rotation of a molecule about its center of mass.
"""
function rotate!(molecule::Molecule; assume_center_of_mass_at_origin::Bool=false)
    # generate a randomly oriented vector. this is the axis about which we rotate
    u = rand_point_on_unit_sphere()
    # generate a rotation matrix for rotating about u for a random θ
    r = rotation_matrix(u, rand() * 2 * π)
    # store the center of mass
    center_of_mass = deepcopy(molecule.center_of_mass)
    # move the molecule to the origin
    translate_by!(molecule, -center_of_mass)
    # conduct the rotation
    for ljsphere in molecule.ljspheres
        ljsphere.x[:] = r * ljsphere.x
    end

    for charge in molecule.charges
        charge.x[:] = r * charge.x
    end
    # no need to rotate center of mass since it is now the origin.
    # translate back to the center of mass.
    translate_by!(molecule, center_of_mass)
end

 # function move_arni!(molecule::Molecule, translation::Array{Float64, 1})
 #     T = construct_translational_matrix(translation[1], translation[2], translation[3])
 #     for ljsphere in molecule.ljspheres
 #         ljsphere.x[:] = matmul(ljsphere.x, T)        
 #     end
 # 
 #     for charge in molecule.charges
 #         charge.x[:] = matmul(charge.x, T)
 #     end
 #     molecule.center_of_mass[:] = matmul(molecule.center_of_mass, T)
 # end
 # 
 # 
 # """
 # Used to adjust to 4x4 translational/rotational matrices
 # """
 # function matmul(vector::Array{Float64,1}, matrix::Array{Float64,2})
 #     vec = deepcopy(vector)
 #     push!(vec, 1)
 #     vec = matrix * vec
 #     return vec[1:3]
 # end
 # 
 # 
 # """
 #     rotate!(molecule::Molecule, θ::Real, ϕ::Real, γ::Real, order_of_rotation::Array{Int64, 1} = [1, 2, 3])
 # 
 # Used to rotate a molecule by the center of mass. Because this is a three dimensional space, the order of rotation matters.
 # 
 # # Arguments
 # - `molecule::Molecule`: The molecule being rotated
 # - `θ::Real`: The rotation about the x-axis in radians
 # - `ϕ::Real`: The rotation about the y-axis in radians
 # - `γ::Real`: The rotation about the z-axis in radians
 # - `order_of_rotation::Array{Int64, 1}`: The order in which the molecule is rotated. Will default to [1,2,3]. That is, it will rotate about the x-axis first, then the y-axis and finally the z-axis.
 # """
 # function rotate!(molecule::Molecule, θ::Real, ϕ::Real, γ::Real, order_of_rotation::Array{Int64, 1} = [1,2,3])
 #     θ = Float64(θ); ϕ = Float64(ϕ); γ = Float64(γ);
 #     #TODO This function is almost identical to the move! function. We could merge them maybe?
 #     R = construct_rot_matrix(molecule.center_of_mass, order_of_rotation, θ, ϕ, γ)
 #     for ljsphere in molecule.ljspheres
 #         ljsphere.x[:] = matmul(ljsphere.x, R)
 #     end
 # 
 #     for charge in molecule.charges
 #         charge.x[:] = matmul(charge.x, R)
 #     end
 # end
 # 
 # 
 # """
 #     rotate_x!(molecule::Molecule, θ::Real)
 # 
 # Used to rotate molecule about the x-axis (see rotate!)
 # """
 # rotate_x!(molecule::Molecule, θ::Real) = rotate!(molecule, θ, 0, 0)
 # 
 # 
 # """
 #     rotate_y!(molecule::Molecule, θ::Real)
 # 
 # Used to rotate molecule about the y-axis (see rotate!)
 # """
 # rotate_y!(molecule::Molecule, θ::Real) = rotate!(molecule, 0, θ, 0)
 # 
 # 
 # """
 #     rotate_z!(molecule::Molecule, θ::Real)
 # 
 # Used to rotate molecule about the z-axis (see rotate!)
 # """
 # rotate_z!(molecule::Molecule, θ::Real) = rotate!(molecule, 0, 0, θ)
 # 
 # """
 # Constructs a 4x4 translational matrix. We use a 4x4 matrix to be able to multiply it with the rotational matrices and combine all translational and rotational operations into a single matrix `R`
 # """
 # function construct_translational_matrix(x::Real, y::Real, z::Real)
 #     return [1 0 0 x;
 #             0 1 0 y;
 #             0 0 1 z;
 #             0 0 0 1]
 # end
 # 
 # """
 # Constructs a 4x4 rotational matrix. Because we're in a three dimensional space the order of rotations matters. We use a 4x4 matrix to be able to multiply it with the rotational matrices and combine all translational and rotational operations into a single matrix `R`
 # """
 # function construct_rot_matrix(center::Array{Float64, 1}, order_of_rotation::Array{Int64, 1}, θ::Real, ϕ::Real, γ::Real)
 #     #TODO optimize
 #     # Translate to origin
 #     R = construct_translational_matrix(-center[1], -center[2], -center[3])
 #     RR = Array{Array{Float64, 2}, 1}(3)
 #     RR[1] = [1 0 0 0; 0 cos(θ) -sin(θ) 0; 0 sin(θ) cos(θ) 0; 0 0 0 1]
 #     RR[2] = [cos(ϕ) 0 sin(ϕ) 0; 0 1 0 0; -sin(ϕ) 0 cos(ϕ) 0; 0 0 0 1]
 #     RR[3] = [cos(γ) -sin(γ) 0 0; sin(γ) cos(γ) 0 0; 0 0 1 0; 0 0 0 1]
 #     for i in order_of_rotation
 #         R = RR[i] * R
 #     end
 #     return construct_translational_matrix(center[1], center[2], center[3]) * R
 # end
