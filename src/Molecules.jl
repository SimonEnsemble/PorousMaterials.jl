"""
Lennard-Jones spheres containing an atom and its' cartesian coordinates

# Attributes
- `atom::Symbol`: Atom name corresponding to force field
- `x::Array{Float64, 1}`: Cartesian coordinates in Angstrom
"""
struct LennardJonesSphere
    atom::Symbol
    x::Array{Float64, 1}
end

function Base.isapprox(ljs1::LennardJonesSphere, ljs2::LennardJonesSphere)
    return ((ljs1.atom == ljs2.atom) & isapprox(ljs1.x, ljs2.x))
end

"""
Point charge data structure containing a charge and the charges' cartesian coordinates

# Attributes
- `q::Float64`: Charge magnitude (electrons)
- `x::Array{Float64, 1}`: Cartesian coordinates in Angstrom
"""
struct PointCharge
    q::Float64
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
    molecule = read_molecule_file("Xe", assert_charge_neutrality=true)

Reads molecule files in the directory `PorousMaterials.PATH_TO_DATA * "/molecule/" * species * "/"`.
Center of mass assigned using atomic masses from `read_atomic_masses()`

# Arguments
- `species::AbstractString`: Name of the molecule
- `assert_charge_neutrality::Bool`: assert the molecule is charge neutral for safety.

# Returns
- `molecule::Molecule`: A fully constructed molecule data structure
"""
function read_molecule_file(species::AbstractString; assert_charge_neutrality::Bool=true)
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

    molecule = Molecule(Symbol(species), ljspheres, charges, COM)

    # check for charge neutrality
    if length(charges) > 0
        if ! (total_charge(molecule) ≈ 0.0)
            if assert_charge_neutrality
                error(@sprintf("Molecule %s is not charge neutral! Pass `assert_charge_neutrality=false` to ignore this error message.", species))
            end
        end
    end

    return molecule
end

"""
    translate!(molecule, dx)

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
    translate_to!(molecule, x)

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
    h = eye(3) - 2 * v * transpose(v)
    return - h * r
end

"""
    rotate!(molecule)

Conduct a uniform random rotation of a molecule about its center of mass.

# Arguments
- `molecule::Molecule`: The molecule which will be rotated
"""
function rotate!(molecule::Molecule)
    # generate a random rotation matrix
    r = rotation_matrix()
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

"""
    outside_box = completely_outside_box(molecule, box)

Checks if a Molecule object is within the boundaries of a Box unitcell.

# Arguments
- `molecule::Molecule`: The molecule object
- `box::Box`: The unit cell object

# Returns
- `outside_box::Bool`: True if the center of mass of `molecule` is outisde of `box`. False otherwise
"""
function outside_box(molecule::Molecule, box::Box)
    xf = box.c_to_f * molecule.center_of_mass
    for k = 1:3
        if (xf[k] > 1.0) | (xf[k] < 0.0)
            return true
        end
    end
    return false
end

"""
    write_to_xyz(molecules, filename; comment="")

Write an array of molecules to an .xyz file. Write only the Lennard-Jones spheres to file (not charges).

# Arguments
- `molecules::Array{Molecule, 1}`: An array of molecules
- `filename::AbstractString`: Name of the output file
- `comment::AbstractString`: A comment that will be printed in the xyz file
"""
function write_to_xyz(molecules::Array{Molecule, 1}, filename::AbstractString; comment::AbstractString="")
    if ! contains(filename, ".xyz")
        filename *= ".xyz"
    end

    n_atoms = sum([length(molecules[i].ljspheres) for i = 1:length(molecules)])

    xyzfile = open(filename, "w")
    @printf(xyzfile, "%d\n%s\n", n_atoms, comment)
    for molecule in molecules
        for ljsphere in molecule.ljspheres
			@printf(xyzfile, "%s\t%.4f\t%.4f\t%.4f\n", string(ljsphere.atom), ljsphere.x[1], ljsphere.x[2], ljsphere.x[3])
        end
    end
    close(xyzfile)
end

"""
    total_charge = total_charge(molecule)

Sum up point charges on a molecule.

# Arguments
- `molecule::Molecule`: the molecule we wish to calculate the total charge of

# Returns
- `total_charge::Float64`: The sum of the point charges of `molecule`
"""
function total_charge(molecule::Molecule)
    total_charge = 0.0
    for charge in molecule.charges
        total_charge += charge.q
    end
    return total_charge
end

function charged(molecule::Molecule; verbose::Bool=false)
    charged_flag = true
    if length(molecule.charges) == 0
        charged_flag = false
    end
    if verbose
        @printf("\tMolecule %s has point charges? %s\n", molecule.species, charged_flag)
    end
    return charged_flag
end
