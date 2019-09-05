"""
Data structure for a molecule/adsorbate.

# Attributes
- `species::Symbol`: Species of molecule, e.g. `:CO2`
- `atoms::Atoms`: array of Lennard-Jones spheres comprising the molecule
- `charges::Charges`: array of point charges comprising the molecule
- `xf_com::Array{Float64, 1}`: center of mass of the molecule in fractional coordinates
"""
struct Molecule
    species::Symbol
    atoms::Atoms
    charges::Charges
    xf_com::Array{Float64, 1}
end

function Base.isapprox(m1::Molecule, m2::Molecule)
    return (m1.species == m2.species) && isapprox(m1.xf_com, m2.xf_com) &&
        (isapprox(m1.atoms, m2.atoms)) &&
        (isapprox(m1.charges, m2.charges))
end

"""
    molecule = Molecule(species, assert_charge_neutrality=true)

Reads molecule files in the directory `joinpath(PATH_TO_MOLECULES, species)`.
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
    if ! isdir(joinpath(PATH_TO_MOLECULES, species))
        error(@sprintf("No directory created for %s in %s\n", species,
                       joinpath(PATH_TO_MOLECULES)))
    end

    ###
    #  Read in Lennard Jones spheres
    ###
    atomsfilename = joinpath(PATH_TO_MOLECULES, species, "lennard_jones_spheres.csv")
    if ! isfile(atomsfilename)
        error(@sprintf("No file %s exists. Even if there are no Lennard Jones spheres in
        %s, include a .csv file with the proper headers but no rows.", species, atomsfilename))
    end
    df_lj = CSV.read(atomsfilename)

    atomic_masses = read_atomic_masses() # for center-of-mass calcs

    x_com = [0.0, 0.0, 0.0]
    total_mass = 0.0

    atom_species = Array{Symbol, 1}()
    atom_coords = Array{Float64, 2}(undef, 3, 0)
    for row in eachrow(df_lj)
        x = [row[:x], row[:y], row[:z]]
        atom = Symbol(row[:atom])
        # assume a unit cube as a box for now.
        if ! (atom in keys(atomic_masses))
            error(@sprintf("Atomic mass of %s not found. See `read_atomic_masses()`\n", atom))
        end
        total_mass += atomic_masses[atom]
        x_com += atomic_masses[atom] .* x
        atom_coords = [atom_coords x]
        push!(atom_species, Symbol(atom))
    end
    x_com /= total_mass
    # construct atoms attribute of molecule
    atoms = Atoms(atom_species, atom_coords)

    ###
    #  Read in point charges
    ###
    chargesfilename = joinpath(PATH_TO_MOLECULES, species, "point_charges.csv")
    if ! isfile(chargesfilename)
        error(@sprintf("No file %s exists. Even if there are no point charges in %s,
        include a .csv file with the proper headers but no rows.", species,
                       chargesfilename))
    end
    df_c = CSV.read(chargesfilename)

    charge_vals = Array{Float64, 1}()
    charge_coords = Array{Float64, 2}(undef, 3, 0)
    for row in eachrow(df_c)
        # assume unit cube as box for now.
        charge_coords = [charge_coords [row[:x], row[:y], row[:z]]]
        push!(charge_vals, row[:q])
    end
    # construct charges attribute of molecule
    charges = Charges(charge_vals, charge_coords)

    # construct molecule
    molecule = Molecule(Symbol(species), atoms, charges, x_com)

    # check for charge neutrality
    if (length(charges.q) > 0) && (! (total_charge(molecule) ≈ 0.0))
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

# Arguments
- `molecule::Molecule`: The molecule which will have its' coordinates set to fractional coordinates
- `box::Box`: The box which will serve as a reference for the fractional coordinates
"""
function set_fractional_coords!(molecule::Molecule, box::Box)
    molecule.atoms.xf[:] = box.c_to_f * molecule.atoms.xf
    molecule.charges.xf[:] = box.c_to_f * molecule.charges.xf
    molecule.xf_com[:] = box.c_to_f * molecule.xf_com
    return nothing
end

"""
    set_fractional_coords_to_unit_cube!(molecule, box)

Change fractional coordinates of a molecule in the context of a given box to Cartesian,
i.e. to correspond to fractional coords in a unit cube box.

# Arguments
- `molecule::Molecule`: The molecule which will have its' coordinates set to fractional coordinates of a unit cube
- `box::Box`: The box which will serve as a reference for the fractional coordinates
"""
function set_fractional_coords_to_unit_cube!(molecule::Molecule, box::Box)
    for i = 1:molecule.atoms.n_atoms
        molecule.atoms.xf[:, i] = box.f_to_c * molecule.atoms.xf[:, i]
    end
    for i = 1:molecule.charges.n_charges
        molecule.charges.xf[:, i] = box.f_to_c * molecule.charges.xf[:, i]
    end
    molecule.xf_com[:] = box.f_to_c * molecule.xf_com
    return nothing
end

"""
    num_atoms = n_atoms(molecules)

calculates the total number of atoms in an array of molecules

# Arguments
 - `molecule::Array{Molecule, 1}`: The molecules to count the number of atoms in
# Returns
 - The total number of atoms in the molecules passed in
"""
function n_atoms(molecules::Array{Molecule, 1})
    num_atoms = 0
    for molecule in molecules
        num_atoms += molecule.atoms.n_atoms
    end
    return num_atoms
end

"""
    write_xyz(box, molecules, xyz_file)

Writes the coordinates of all atoms in molecules to the given xyz_file file object
passing a file object around is faster for simulation because it can be opened
once at the beginning of the simulation and closed at the end.

This writes the coordinates of the molecules in cartesian coordinates, so the
box is needed for the conversion.

# Arguments
 - `box::Box`: The box the molecules are in, to convert molecule positions
        to cartesian coordinates
 - `molecules::Array{Molecule, 1}`: The array of molecules to be written to the file
 - `xyz_file::IOStream`: The open 'write' file stream the data will be saved to
"""
function write_xyz(box::Box, molecules::Array{Molecule, 1}, xyz_file::IOStream)
    num_atoms = n_atoms(molecules)
    @printf(xyz_file, "%s\n", num_atoms)
    for molecule in molecules
        for i = 1:molecule.atoms.n_atoms
            cartesian_coords = box.f_to_c * molecule.atoms.xf[:, i]
            @printf(xyz_file, "\n%s %f %f %f", molecule.atoms.species[i],
                    cartesian_coords...)
        end
    end
end

"""
    translate_by!(molecule, dxf)
    translate_by!(molecule, dx, box)

Translate a molecule by vector `dxf` in fractional coordinate space or by vector `dx` in
Cartesian coordinate space. For the latter, a unit cell box is required for context.

# Arguments
- `molecule::Molecule`: The molecule which will be translated by `dxf`
- `dxf::Array{Float64, 1}`: A vector containing the translation in each dimension
"""
function translate_by!(molecule::Molecule, dxf::Array{Float64, 1})
    # translate both atoms and charges
    molecule.atoms.xf[:] = broadcast(+, molecule.atoms.xf, dxf)
    molecule.charges.xf[:] = broadcast(+, molecule.charges.xf, dxf)
    # adjust center of mass
    molecule.xf_com[:] += dxf
    return nothing
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

# Arguments
- `molecule::Molecule`: The molecule which will be translated to `xf`
- `xf::Array{Float64, 1}`: A vector containing the coordinates of the final destination of the molecule
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
    if molecule.atoms.n_atoms > 0
        print(io, "Atoms:\n")
        for i = 1:molecule.atoms.n_atoms
            @printf(io, "\n\tatom = %s, xf = [%.3f, %.3f, %.3f]", molecule.atoms.species[i],
                    molecule.atoms.xf[1, i], molecule.atoms.xf[2, i], molecule.atoms.xf[3, i])
        end
    end
    if molecule.charges.n_charges > 0
        print(io, "\nPoint charges: ")
        for i = 1:molecule.charges.n_charges
            @printf(io, "\n\tcharge = %f, xf = [%.3f, %.3f, %.3f]", molecule.charges.q[i],
                    molecule.charges.xf[1, i], molecule.charges.xf[2, i], molecule.charges.xf[3, i])
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

Conduct a random rotation of the molecule about its center of mass.
The box is needed because the molecule contains only its fractional coordinates.

# Arguments
- `molecule::Molecule`: The molecule which will be subject to a random rotation
- `box::Box`: The molecule only contains fractional coordinates, so the box is needed for a correct rotation
"""
function rotate!(molecule::Molecule, box::Box)
    # generate a random rotation matrix
    #    but use c_to_f, f_to_c for fractional
    r = rotation_matrix()
    r = box.c_to_f * r * box.f_to_c
    # conduct the rotation
    # shift to origin
    molecule.atoms.xf[:] = broadcast(-, molecule.atoms.xf, molecule.xf_com)
    molecule.charges.xf[:] = broadcast(-, molecule.charges.xf, molecule.xf_com)
    # conduct the rotation
    molecule.atoms.xf[:] = r * molecule.atoms.xf
    molecule.charges.xf[:] = r * molecule.charges.xf
    # shift back to center of mass
    molecule.atoms.xf[:] = broadcast(+, molecule.atoms.xf, molecule.xf_com)
    molecule.charges.xf[:] = broadcast(+, molecule.charges.xf, molecule.xf_com)
    return nothing
end

"""
    outside_box = outside_box(molecule)

Checks if a Molecule object is within the boundaries of a Box unitcell.

# Arguments
- `molecule::Molecule`: The molecule object

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
function write_xyz(molecules::Array{Molecule, 1}, box::Box, filename::AbstractString;
    comment::AbstractString="")

    n_atoms = sum([molecule.atoms.n_atoms for molecule in molecules])
    atoms = Array{Symbol}(undef, n_atoms)
    x = zeros(Float64, 3, n_atoms)

    atom_counter = 0
    for molecule in molecules
        for i = 1:molecule.atoms.n_atoms
            atom_counter += 1
            x[:, atom_counter] = box.f_to_c * molecule.atoms.xf[:, i]
            atoms[atom_counter] = molecule.atoms.species[i]
        end
    end
    @assert(atom_counter == n_atoms)

    write_xyz(atoms, x, filename, comment=comment) # Misc.jl
end

"""
    total_charge = total_charge(molecule)

Sum up point charges on a molecule.

# Arguments
- `molecule::Molecule`: the molecule we wish to calculate the total charge of

# Returns
- `total_charge::Float64`: The sum of the point charges of `molecule`
"""
total_charge(molecule::Molecule) = (molecule.charges.n_charges == 0) ? 0.0 : sum(molecule.charges.q)

"""
    charged_flag = charged(molecule, verbose=false)

Determine if a molecule has point charges

# Arguments
- `molecule::Molecule`: The molecule which will be checked for charges
- `verbose::Bool`: Will print result if `true`

# Returns
- `charged_flag::Bool`: `true` if molecule is charged, `false` otherwise
"""
function charged(molecule::Molecule; verbose::Bool=false)
    charged_flag = molecule.charges.n_charges > 0
    if verbose
        @printf("\tMolecule %s has point charges? %s\n", molecule.species, charged_flag)
    end
    return charged_flag
end

"""
    bond_lenghts = pairwise_atom_distances(molecule, box) 

Loop over all pairs of `Atoms`'s in `molecule.atoms`. Return a matrix whose `(i, j)`
element is the distance between atom `i` and atom `j` in the molecule.

# Arguments
- `molecule::Molecule`: The molecules containing the atoms we wish to calculate the distances between
- `box::Box`: The box we use to get the cartesian coordinates rather than the fractional coordinates stored in the Molecule object

# Returns
- `bond_lenghts::Array{Float64, 2}`: A `n_atoms` by `n_atoms` symmetric matrix containing the bond lengths between atoms in `molecule`
"""
function pairwise_atom_distances(molecule::Molecule, box::Box)
    nb_atoms = molecule.atoms.n_atoms
    bond_lengths = zeros(nb_atoms, nb_atoms)
    for i = 1:nb_atoms
        for j = (i+1):nb_atoms
            dx = box.f_to_c * (molecule.atoms.xf[:, i] - molecule.atoms.xf[:, j])
            bond_lengths[i, j] = norm(dx)
            bond_lengths[j, i] = bond_lengths[i, j]
        end
    end
    return bond_lengths
end

"""
    charge_bond_lengths = pairwise_charge_distances(molecule, box)

Loop over all pairs of `Charges`'s in `molecule.charges`. Return a matrix whose `(i, j)`
element is the distance between point charge `i` and point charge `j` in the molecule.

# Arguments
- `molecule::Molecule`: The molecules containing the atoms we wish to calculate the distances between
- `box::Box`: The box we use to get the cartesian coordinates rather than the fractional coordinates stored in the Molecule object

# Returns
- `charge_bond_lenghts::Array{Float64, 2}`: A `n_charges` by `n_charges` symmetric matrix containing the point charge bond lengths between atoms in `molecule`
"""
function pairwise_charge_distances(molecule::Molecule, box::Box)
    nb_charges = molecule.charges.n_charges
    bond_lengths = zeros(nb_charges, nb_charges)
    for i = 1:nb_charges
        for j = (i+1):nb_charges
            dx = box.f_to_c * (molecule.charges.xf[:, i] - molecule.charges.xf[:, j])
            bond_lengths[i, j] = norm(dx)
            bond_lengths[j, i] = bond_lengths[i, j]
        end
    end
    return bond_lengths
end

"""
    drift = bond_length_drift(molecule, reference_molecule, box, atol=1e-14, throw_warning=true)

Compute pairwise atom & charge distances of `molecule` and compare to those in a reference 
molecule to determine if the pairwise atom & charge distances differ within a tolerance 
`atol`. This is useful for checking for drift in the course of the simulation, during 
which rotations and translations are performed.

# Arguments
- `molecule::Molecule`: The molecule which is subject to bond length drift check
- `reference_molecule::Molecule`: The molecule we will compare `molecule` with
- `box::Box`: The box structure we use to transform to cartesian coordinates
- `atol::Float64`: The tolerance for the bond length drift check
- `throw_warning::Bool`: Will throw warning if `true`

# Returns
- `drift::Bool`: If `true` we have observed a bond length drift. No bond length drift otherwise
"""
function bond_length_drift(molecule::Molecule, ref_molecule::Molecule, box::Box;
    atol::Float64=1e-14, throw_warning::Bool=true)
    @assert molecule.species == ref_molecule.species

    we_have_drift = false

    if ! isapprox(pairwise_atom_distances(ref_molecule, box),
                  pairwise_atom_distances(molecule, box), atol=atol)
        we_have_drift = true
    end

    if ! isapprox(pairwise_charge_distances(ref_molecule, box),
                  pairwise_charge_distances(molecule, box), atol=atol)
        we_have_drift = true
    end

    if throw_warning && we_have_drift
        @warn @sprintf("pairwise charge/atom distances in %s drifted beyond %e tolerance.\n",
            molecule.species, atol)
    end
    return we_have_drift # if made it this far, no bond length drift
end

# facilitate constructing a point charge
Ion(q::Float64, xf::Array{Float64, 1}, species::Symbol=:ion) = Molecule(
    species,
    Atoms(0, Symbol[], zeros(0, 0)),
    Charges(1, [q], reshape(xf, (3, 1))),
    xf)
