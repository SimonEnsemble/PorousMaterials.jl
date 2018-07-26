"""
Data structure for a molecule/adsorbate.

# Attributes
- `species::Symbol`: Species of molecule, e.g. `:CO2`
- `atoms::Array{LJSphere, 1}`: array of Lennard-Jones spheres comprising the molecule
- `charges::Array{PtCharge, 1}`: array of point charges comprising the molecule
- `xf_com::Array{Float64, 1}`: center of mass of the molecule in fractional coordinates
"""
struct Molecule
    species::Symbol
    atoms::Array{LJSphere, 1}
    charges::Array{PtCharge, 1}
    xf_com::Array{Float64, 1}
end

function Base.isapprox(m1::Molecule, m2::Molecule)
    return (m1.species == m2.species) && isapprox(m1.xf_com, m2.xf_com) &&
        (length(m1.atoms) == length(m2.atoms)) && 
        (length(m1.charges) == length(m2.charges)) && 
        all([isapprox(m1.atoms[i], m2.atoms[i]) for i = 1:length(m1.atoms)]) &&
        all([isapprox(m1.charges[i], m2.charges[i]) for i = 1:length(m1.charges)])
end

"""
    molecule = Molecule(species, assert_charge_neutrality=true)

Reads molecule files in the directory `PorousMaterials.PATH_TO_DATA * "/molecule/" * species * "/"`.
Center of mass assigned using atomic masses from `read_atomic_masses()`. The fractional
coordinates are determined assuming a unit cube box. These must be adjusted later for
simulations using `set_fractional_coords!(molecule, box)`.

# Arguments
- `species::AbstractString`: Name of the molecule
- `assert_charge_neutrality::Bool`: assert the molecule is charge neutral for safety.

# Returns
- `molecule::Molecule`: A fully constructed molecule data structure
"""
function Molecule(species::AbstractString; assert_charge_neutrality::Bool=true)
    if ! isdir(PATH_TO_DATA * "molecules/" * species)
        error(@sprintf("No directory created for %s in %s\n", species,
                       PATH_TO_DATA * "molecules/"))
    end
    
    # Read in Lennard Jones spheres
    atomsfilename = PATH_TO_DATA * "molecules/" * species * "/lennard_jones_spheres.csv"
    if ! isfile(atomsfilename)
        error(@sprintf("No file %s exists. Even if there are no Lennard Jones spheres in 
        %s, include a .csv file with the proper headers but no rows.", species, atomsfilename))
    end
    df_lj = CSV.read(atomsfilename)
    
    atomic_masses = read_atomic_masses() # for center-of-mass calcs

    x_com = [0.0, 0.0, 0.0]
    total_mass = 0.0

    atoms = LJSphere[]
    for row in eachrow(df_lj)
        x = [row[:x], row[:y], row[:z]]
        atom = Symbol(row[:atom])
        # assume a unit cube as a box for now.
        push!(atoms, LJSphere(atom, x))
        if ! (atom in keys(atomic_masses))
            error(@sprintf("Atomic mass of %s not found. See `read_atomic_masses()`\n", atom))
        end
        total_mass += atomic_masses[atom]
        x_com += atomic_masses[atom] .* x
    end
    x_com /= total_mass

    # Read in point charges
    chargesfilename = PATH_TO_DATA * "molecules/" * species * "/point_charges.csv"
    if ! isfile(chargesfilename)
        error(@sprintf("No file %s exists. Even if there are no point charges in %s, 
        include a .csv file with the proper headers but no rows.", species,
                       chargesfilename))
    end
    df_c = CSV.read(chargesfilename)

    charges = PtCharge[]
    for row in eachrow(df_c)
        # assume unit cube as box for now.
        x = [row[:x], row[:y], row[:z]]
        push!(charges, PtCharge(row[:q], x))
    end

    molecule = Molecule(Symbol(species), atoms, charges, x_com)

    # check for charge neutrality
    if (length(charges) > 0) && (! (total_charge(molecule) ≈ 0.0))
        if assert_charge_neutrality
            error(@sprintf("Molecule %s is not charge neutral! Pass 
            `assert_charge_neutrality=false` to ignore this error message.", species))
        end
    end

    return molecule
end

"""
    set_fractional_coords!(molecule, box)

After a molecule is freshly constructed, its fractional coords are assumed to correspond
to a unit cell box that is a unit cube. This function adjusts the fractional coordinates
of the molecule to be consistent with a different box.
"""
function set_fractional_coords!(molecule::Molecule, box::Box)
    for ljsphere in molecule.atoms
        ljsphere.xf[:] = box.c_to_f * ljsphere.xf
    end
    for charge in molecule.charges
        charge.xf[:] = box.c_to_f * charge.xf
    end
    molecule.xf_com[:] = box.c_to_f * molecule.xf_com
    return nothing
end

"""
    set_fractional_coords_to_unit_cube!(molecule, box)

Change fractional coordinates of a molecule in the context of a given box to Cartesian, 
i.e. to correspond to fractional coords in a unit cube box.
"""
function set_fractional_coords_to_unit_cube!(molecule::Molecule, box::Box)
    for ljsphere in molecule.atoms
        ljsphere.xf[:] = box.f_to_c * ljsphere.xf
    end
    for charge in molecule.charges
        charge.xf[:] = box.f_to_c * charge.xf
    end
    molecule.xf_com[:] = box.f_to_c * molecule.xf_com
    return nothing
end

"""
    translate_by!(molecule, dxf)
    translate_by!(molecule, dx, box)

Translate a molecule by vector `dxf` in fractional coordinate space or by vector `dx` in
Cartesian coordinate space. For the latter, a unit cell box is required for context.
"""
function translate_by!(molecule::Molecule, dxf::Array{Float64, 1})
    # move LJSphere's
    for ljsphere in molecule.atoms
        ljsphere.xf[:] += dxf
    end
    # move PtCharge's
    for charge in molecule.charges
        charge.xf[:] += dxf
    end
    # adjust center of mass
    molecule.xf_com[:] += dxf
end

function translate_by!(molecule::Molecule, dx::Array{Float64, 1}, box::Box)
    # determine shift in fractional coordinate space
    dxf = box.c_to_f * dx
    translate_by!(molecule, dxf)
end

"""
    translate_to!(molecule, xf)
    translate_to!(molecule, x, box)

Translate a molecule a molecule to point `xf` in fractional coordinate space or to `x` in
Cartesian coordinate space. For the latter, a unit cell box is required for context. The
molecule is translated such that its center of mass is at `xf`/x`.
"""
function translate_to!(molecule::Molecule, xf::Array{Float64, 1})
    dxf = xf - molecule.xf_com
    translate_by!(molecule, dxf)
end

function translate_to!(molecule::Molecule, x::Array{Float64, 1}, box::Box)
    translate_to!(molecule, box.c_to_f * x)
end

function Base.show(io::IO, molecule::Molecule)
    println(io, "Molecule species: ", molecule.species)
    println(io, "Center of mass (fractional coords): ", molecule.xf_com)
    if length(molecule.atoms) > 0
        print(io, "Lennard-Jones spheres: ")
        for ljsphere in molecule.atoms
            @printf(io, "\n\tatom = %s, xf = [%.3f, %.3f, %.3f]", ljsphere.species,
                    ljsphere.xf[1], ljsphere.xf[2], ljsphere.xf[3])
        end
    end
    if length(molecule.charges) > 0
        print(io, "\nPoint charges: ")
        for charge in molecule.charges
            @printf(io, "\n\tq = %f, xf = [%.3f, %.3f, %.3f]", charge.q, 
                    charge.xf[1], charge.xf[2], charge.xf[3])
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
    rotate!(molecule, box)

Conduct a random rotation of the molecule about its center of mass. 
The box is needed because the molecule contains only its fractional coordinates.
"""
function rotate!(molecule::Molecule, box::Box)
    # generate a random rotation matrix
    #    but use c_to_f, f_to_c for fractional
    r = rotation_matrix()
    r = box.c_to_f * r * box.f_to_c
    # conduct the rotation
    for ljsphere in molecule.atoms
        ljsphere.xf[:] = molecule.xf_com + r * (ljsphere.xf - molecule.xf_com)
    end
    for charge in molecule.charges
        charge.xf[:] = molecule.xf_com + r * (charge.xf - molecule.xf_com)
    end
end

"""
    outside_box = completely_outside_box(molecule)

Checks if a Molecule object is within the boundaries of a Box unitcell.

# Arguments
- `molecule::Molecule`: The molecule object
- `box::Box`: The unit cell object

# Returns
- `outside_box::Bool`: True if the center of mass of `molecule` is outisde of `box`. False otherwise
"""
function outside_box(molecule::Molecule)
    for k = 1:3
        if (molecule.xf_com[k] > 1.0) || (molecule.xf_com[k] < 0.0)
            return true
        end
    end
    return false
end

# docstring in Misc.jl
function write_to_xyz(molecules::Array{Molecule, 1}, box::Box, filename::AbstractString; 
    comment::AbstractString="")

    n_atoms = sum([length(molecule.atoms) for molecule in molecules])
    atoms = Array{Symbol}(n_atoms)
    x = zeros(Float64, 3, n_atoms)

    atom_counter = 0
    for molecule in molecules
        for atom in molecule.atoms
            atom_counter += 1
            x[:, atom_counter] = box.f_to_c * atom.xf
            atoms[atom_counter] = atom.species
        end
    end
    @assert(atom_counter == n_atoms)
    
    write_to_xyz(atoms, x, filename, comment=comment) # Misc.jl
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
    charged_flag = length(molecule.charges) > 0
    if verbose
        @printf("\tMolecule %s has point charges? %s\n", molecule.species, charged_flag)
    end
    return charged_flag
end

function bond_lengths(molecule::Molecule, box::Box)
    bond_lengths = Float64[]
    for i = 1:length(molecule.atoms)
        for j = (i+1):length(molecule.atoms)
            dx = box.f_to_c * (molecule.atoms[i].xf - molecule.atoms[j].xf)
            push!(bond_lengths, norm(dx))
        end
    end
    return bond_lengths
end
