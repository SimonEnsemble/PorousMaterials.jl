function _define_atomic_mass()
    global ATOMIC_MASS = read_atomic_masses() # for center-of-mass calcs
end

"""
Data structure for a molecule/adsorbate.

# Attributes
- `species::Symbol`: Species of molecule, e.g. `:CO2`
- `atoms::Atoms`: array of Lennard-Jones spheres comprising the molecule
- `charges::Charges`: array of point charges comprising the molecule
- `com::Coords`: center of mass
"""
struct Molecule{T} # T = Frac or Cart
    species::Symbol
    atoms::Atoms{T}
    charges::Charges{T}
    com::T # center of mass
end

function Base.isapprox(m1::Molecule, m2::Molecule)
    return (m1.species == m2.species) && isapprox(m1.com, m2.com) &&
         isapprox(m1.atoms, m2.atoms) && isapprox(m1.charges, m2.charges)
end

function center_of_mass(molecule::Molecule{Cart})
    if ! @isdefined ATOMIC_MASS
        _define_atomic_mass()
    end
    total_mass = 0.0 # total mass
    x_com = [0.0, 0.0, 0.0] # center of mass
    for a = 1:molecule.atoms.n
        m = ATOMIC_MASS[molecule.atoms.species[a]]
        total_mass += m 
        x_com += m * molecule.atoms.coords.x[:, a]
    end
    return Cart(x_com / total_mass)
end
    
"""
    molecule = Molecule(species, assert_charge_neutrality=true)

construt a `Molecule` from input files read from:
`joinpath(PATH_TO_MOLECULES, species)`.

center of mass assigned using atomic masses from `read_atomic_masses()`. 

# Arguments
- `species::String`: name of the molecule
- `check_neutrality::Bool`: assert the molecule is charge neutral for safety.

# Returns
- `molecule::Molecule{Cart}`: a molecule in Cartesian coordinates
"""
function Molecule(species::String; check_neutrality::Bool=true)
    ###
    #  Read in Lennard Jones spheres
    ###
    df_lj = CSV.read(joinpath(PATH_TO_MOLECULES, species, "atoms.csv"))
    atoms = Atoms{Cart}(nrow(df_lj)) # pre-allocate atoms

    for (a, row) in enumerate(eachrow(df_lj))
        atoms.species[a] = Symbol(row[:atom])
        atoms.coords[a] = [row[:x], row[:y], row[:z]]
    end

    ###
    #  Read in point charges
    ###
    df_c = CSV.read(joinpath(PATH_TO_MOLECULES, species, "charges.csv"))
    charges = Charges{Cart}(nrow(df_c)) # pre-allocate charges

    for (c, row) in enumerate(eachrow(df_c))
        charges.q[c] = row[:q]
        charges.coords.x[:, c] = [row[:x], row[:y], row[:z]]
    end
    
    molecule = Molecule(Symbol(species), atoms, charges, Cart([NaN, NaN, NaN]))

    # compute center of mass
    molecule.com.x .= center_of_mass(molecule).x

    # check for charge neutrality
    if (! neutral(molecule.charges)) && check_neutrality
        error(@sprintf("Molecule %s is not charge neutral! Pass
        `check_neutrality=false` to ignore this error message.", species))
    end

    return molecule
end

# documented in matter.jl
function wrap!(molecule::Molecule{Frac})
    wrap!(molecule.atoms.coords)
    wrap!(molecule.charges.coords)
    wrap!(molecule.com)
end

# documented in matter.jl
net_charge(molecule::Molecule) = net_charge(molecule.charges)

# convert a molecule to fractional coordinates
function Frac(molecule::Molecule{Cart}, box::Box)
    return Molecule(molecule.species,
        Frac(molecule.atoms, box),
        Frac(molecule.charges, box),
        Frac(molecule.com, box)
        )
end

# convert a molecule to cartesian coordinates
function Cart(molecule::Molecule{Frac}, box::Box)
    return Molecule(molecule.species,
        Cart(molecule.atoms, box),
        Cart(molecule.charges, box),
        Cart(molecule.com, box)
       )
end

 # """
 #     write_xyz(molecules, xyz_file)
 # 
 # Writes the coordinates of all atoms in molecules to the given xyz_file file object
 # passing a file object around is faster for simulation because it can be opened
 # once at the beginning of the simulation and closed at the end.
 # 
 # This writes the coordinates of the molecules in cartesian coordinates, so the
 # box is needed for the conversion.
 # 
 # # Arguments
 #  - `box::Box`: The box the molecules are in, to convert molecule positions
 #         to cartesian coordinates
 #  - `molecules::Array{Molecule, 1}`: The array of molecules to be written to the file
 #  - `xyz_file::IOStream`: The open 'write' file stream the data will be saved to
 # """
 # function write_xyz(box::Box, molecules::Array{Molecule, 1}, xyz_file::IOStream)
 #     num_atoms = n_atoms(molecules)
 #     @printf(xyz_file, "%s\n", num_atoms)
 #     for molecule in molecules
 #         for i = 1:molecule.atoms.n_atoms
 #             cartesian_coords = box.f_to_c * molecule.atoms.xf[:, i]
 #             @printf(xyz_file, "\n%s %f %f %f", molecule.atoms.species[i],
 #                     cartesian_coords...)
 #         end
 #     end
 # end

# documented in matter.jl
function translate_by!(molecule::Molecule{Cart}, dx::Cart)
    translate_by!(molecule.atoms.coords,   dx)
    translate_by!(molecule.charges.coords, dx)
    translate_by!(molecule.com,            dx)
end

function translate_by!(molecule::Molecule{Frac}, dxf::Frac)
    translate_by!(molecule.atoms.coords,   dxf)
    translate_by!(molecule.charges.coords, dxf)
    translate_by!(molecule.com,            dxf)
end

function translate_by!(molecule::Molecule{Cart}, dxf::Frac, box::Box)
    translate_by!(molecule, Cart(dxf, box))
end

function translate_by!(molecule::Molecule{Frac}, dx::Cart, box::Box)
    translate_by!(molecule, Frac(dx, box))
end

"""
    translate_to!(molecule, xf)
    translate_to!(molecule, x)
    translate_to!(molecule, xf, box)
    translate_to!(molecule, x, box)

Translate a molecule a molecule to point `xf` in fractional coordinate space or to `x` in
Cartesian coordinate space. For the latter, a unit cell box is required for context. The
molecule is translated such that its center of mass is at `xf`/x`.

# Arguments
- `molecule::Molecule`: The molecule which will be translated to `xf`
- `xf::Array{Float64, 1}`: A vector containing the coordinates of the final destination of the molecule
"""
function translate_to!(molecule::Molecule{Cart}, x::Cart)
    dx = Cart(broadcast(-, x.x, molecule.com.x))
    translate_by!(molecule, dx)
end

function translate_to!(molecule::Molecule{Frac}, xf::Frac)
    dxf = Frac(broadcast(-, xf.xf, molecule.com.xf))
    translate_by!(molecule, dxf)
end

translate_to!(molecule::Molecule{Frac}, x::Cart, box::Box) = translate_to!(molecule, Frac(box.c_to_f * x.x))

translate_to!(molecule::Molecule{Cart}, xf::Frac, box::Box) = translate_to!(molecule, Cart(box.f_to_c * xf.xf))

function Base.show(io::IO, molecule::Molecule)
    println(io, "Molecule species: ", molecule.species)
    println(io, "Center of mass (fractional coords): ", molecule.com)
    if molecule.atoms.n > 0
        print(io, "Atoms:\n")
        if typeof(molecule.atoms.coords) == Frac
            for i = 1:molecule.atoms.n
                @printf(io, "\n\tatom = %s, xf = [%.3f, %.3f, %.3f]", molecule.atoms.species[i],
                        molecule.atoms.coords[i].xf...)
            end
        elseif typeof(molecule.atoms.coords) == Cart
            for i = 1:molecule.atoms.n
                @printf(io, "\n\tatom = %s, x = [%.3f, %.3f, %.3f]", molecule.atoms.species[i],
                        molecule.atoms.coords[i].x...)
            end
        end
    end
    if molecule.charges.n > 0
        print(io, "\nPoint charges: ")
        if typeof(molecule.charges.coords) == Frac
            for i = 1:molecule.charges.n
                @printf(io, "\n\tcharge = %f, xf = [%.3f, %.3f, %.3f]", molecule.charges.q[i],
                        molecule.charges.coords[i].xf...)
            end
        elseif typeof(molecule.charges.coords) == Cart
            for i = 1:molecule.charges.n
                @printf(io, "\n\tcharge = %f, x = [%.3f, %.3f, %.3f]", molecule.charges.q[i],
                        molecule.charges.coords[i].x...)
            end
        end
    end
end

"""
    u = rand_point_on_unit_sphere()

Generate a unit vector with a random orientation.

# Returns
- `u::Array{Float64, 1}`: A unit vector with a random orientation
"""
function rand_point_on_unit_sphere()
    u = randn(3)
    u_norm = norm(u)
    if u_norm < 1e-6 # avoid numerical error in division
        return rand_point_on_unit_sphere()
    end
    return u / u_norm
end

"""
    r = rotation_matrix()

Generate a 3x3 random rotation matrix `r` such that when a point `x` is rotated using this rotation matrix via `r * x`, this point `x` is placed at a uniform random distributed position on the surface of a sphere of radius `norm(x)`.
See James Arvo. Fast Random Rotation Matrices.

https://pdfs.semanticscholar.org/04f3/beeee1ce89b9adf17a6fabde1221a328dbad.pdf

# Returns
- `r::Array{Float64, 2}`: A 3x3 random rotation matrix
"""
function rotation_matrix()
    # random rotation about the z-axis
    u₁ = rand() * 2.0 * π
    r = [cos(u₁) sin(u₁) 0.0; -sin(u₁) cos(u₁) 0.0; 0.0 0.0 1.0]

    # househoulder matrix
    u₂ = 2.0 * π * rand()
    u₃ = rand()
    v = [cos(u₂) * sqrt(u₃), sin(u₂) * sqrt(u₃), sqrt(1.0 - u₃)]
    h = Matrix{Float64}(I, 3, 3) - 2 * v * transpose(v)
    return - h * r
end

"""
    rotate!(molecule, box)
    rotate!(molecule)

Conduct a random rotation of the molecule about its center of mass.
The box is needed when the molecule contains fractional coordinates,
but if cartesian coordinates are used, only the molecule argument suffices.

# Arguments
- `molecule::Molecule`: The molecule which will be subject to a random rotation
- `box::Box`: The molecule only contains fractional coordinates, so the box is needed for a correct rotation
"""
function rotate!(molecule::Molecule{Frac}, box::Box)
    #TODO make this work with Cart and Frac
    # generate a random rotation matrix
    #    but use c_to_f, f_to_c for fractional
    r = rotation_matrix()
    r = box.c_to_f * r * box.f_to_c
    # conduct the rotation
    # shift to origin
    molecule.atoms.coords.xf[:,:] = broadcast(-, molecule.atoms.coords.xf, molecule.com.xf)
    molecule.charges.coords.xf[:,:] = broadcast(-, molecule.charges.coords.xf, molecule.com.xf)
    # conduct the rotation
    molecule.atoms.coords.xf[:,:] = r * molecule.atoms.coords.xf
    molecule.charges.coords.xf[:,:] = r * molecule.charges.coords.xf
    # shift back to center of mass
    molecule.atoms.coords.xf[:,:] = broadcast(+, molecule.atoms.coords.xf, molecule.com.xf)
    molecule.charges.coords.xf[:,:] = broadcast(+, molecule.charges.coords.xf, molecule.com.xf)
    return nothing
end

function rotate!(molecule::Molecule{Cart})
    r = rotation_matrix()
    # conduct the rotation
    # shift to origin
    molecule.atoms.coords.x[:,:] = broadcast(-, molecule.atoms.coords.x, molecule.com.x)
    molecule.charges.coords.x[:,:] = broadcast(-, molecule.charges.coords.x, molecule.com.x)
    # conduct the rotation
    molecule.atoms.coords.x[:,:] = r * molecule.atoms.coords.x
    molecule.charges.coords.x[:,:] = r * molecule.charges.coords.x
    # shift back to center of mass
    molecule.atoms.coords.x[:,:] = broadcast(+, molecule.atoms.coords.x, molecule.com.x)
    molecule.charges.coords.x[:,:] = broadcast(+, molecule.charges.coords.x, molecule.com.x)
    return nothing
end



inside(molecule::Molecule{Cart}, box::Box) = inside(molecule.atoms.coords, box) && inside(molecule.charges.coords, box)
inside(molecule::Molecule{Frac}) = inside(molecule.atoms.coords) && inside(molecule.charges.coords)

# docstring in Misc.jl
function write_xyz(molecules::Array{<:Molecule, 1}, box::Box, filename::AbstractString;
    comment::AbstractString="")
    
    # append all atoms of the molecule together
    atoms = sum([molecule.atoms for molecule in molecules])
    
    if isa(atoms.coords, Frac)
        # convert to Cartesian
        atoms = Cart(atoms, box)
    end
    
    # send to write_xyz for writing atoms in Cartesian coords.
    write_xyz(atoms, filename, comment=comment) # Misc.jl
end

# documented in crystal.jl
has_charges(molecule::Molecule) = molecule.charges.n > 0

# facilitate constructing a point charge
function Ion(q::Float64, coords::Frac)
    @assert size(coords.xf, 2) == 1
    return Molecule(:ion, Charges([q], coords), Atoms{Frac}(0), coords)
end
