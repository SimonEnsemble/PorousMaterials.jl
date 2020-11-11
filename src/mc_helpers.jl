# δ is half the maximal distance a particle is perturbed in a given coordinate
#  during particle translations
const δ = 2.0 # Å

"""
    needs_rotations(molecule) # true or false

Determine whether a molecule needs to undergo rotation.

    `true` if molecule.atoms.n + molecule.charges.n > 1
"""
needs_rotations(molecule::Molecule) = molecule.atoms.n + molecule.charges.n > 1

# break collection of statistics into blocks to gauge convergence and compute standard err
const N_BLOCKS = 5

"""
    random_insertion!(molecules::Array{Molecule{Frac}, 1}, box::Box, template::Molecule{Cart})

Insert a molecule into the simulation box and perform a random rotation if needed.

# Arguments
- `molecules::Array{Molecule{Frac}, 1}`: array containing the molecules in the simulation
- `box::Box`: the box  used for fractional coordinats
- `template::Molecule{Cart}`: reference molecule of the type inserted
"""
function random_insertion!(molecules::Array{Molecule{Frac}, 1}, box::Box, template::Molecule{Cart})
    # copy template
    molecule = deepcopy(template)
    # rotate
    if needs_rotations(molecule)
        random_rotation!(molecule)
    end
    # convert to fractional
    molecule = Frac(molecule, box)
    # translate to uniform random fractional coords in the box
    com = Frac(rand(3))
    translate_to!(molecule, com)
    # add the molecule to the array of molecules
    push!(molecules, molecule)
end

"""
    remove_molecule!(molecule_id, molecules)

Remove a molecule from the array of molecules.

# Arguments
- `molecule_id::Int`: the ID of the molecule to be removed
- `molecules::Array{<:Molecule, 1}`: array of molecules to be modified
"""
function remove_molecule!(molecule_id::Int, molecules::Array{<:Molecule, 1})
    splice!(molecules, molecule_id)
end

"""
    apply_periodic_boundary_condition!(molecule::Molecule{Frac})

Check if each of a molecule's center-of-mass coordinates is within the bounds (0.0, 1.0) in fractional coordinates, and translate if needed.  

# Arguments
- `molecule::Molecule{Frac}`: the molecule to be checked

# Returns
- `nothing`, if the molecule is within the boundary; ohterwise, the coordinates of the input molecule will be modified
"""
function apply_periodic_boundary_condition!(molecule::Molecule{Frac})
    # based on center of mass
    if inside(molecule.com)
        return nothing
    end

    # current center of mass in fractional coordinates; adjust inside loop
    new_com = deepcopy(molecule.com)

    # apply periodic boundary conditions
    for k = 1:3 # loop over xf, yf, zf components
        # if > 1.0, shift down
        if new_com.xf[k] >= 1.0
            new_com.xf[k] -= 1.0
        elseif new_com.xf[k] < 0.0
            new_com.xf[k] += 1.0
        end
    end

    # translate molecule to new center of mass if it was found to be outside of the box
    translate_to!(molecule, new_com)
end

"""
    random_reinsertion!(molecule, box)

Perform a translational perturbation in Cartesian coordinates on a molecule, apply the periodic boundry conditions, and keep a copy of the original in case it needs to be restored.

# Arguements
- `molecule::Molecule{Frac}`: molecule to be translated
- `box::Box`: the box used in rotation

# Returns
- `old_molecule::Molecule{Frac}`: a copy of the original molecule
"""
function random_translation!(molecule::Molecule{Frac}, box::Box)
    # store old molecule and return at the end for possible restoration
    old_molecule = deepcopy(molecule)

    # peturb in Cartesian coords in a random cube centered at current coords.
    dx = Cart(δ * (rand(3) .- 0.5)) # move every atom of the molecule by the same vector.
    translate_by!(molecule, dx, box)
    
    # done, unless the molecule has moved outside of the box, then apply PBC
    apply_periodic_boundary_condition!(molecule)

    return old_molecule # in case we need to restore coords
end

"""
    random_reinsertion!(molecule, box)

Perform a translation and rotated (if needed) on a molecule, and keep a copy of the original in case it needs to be restored.

# Arguements
- `molecule::Molecule{Frac}`: molecule to be translated and rotated (if needed)
- `box::Box`: the box used in rotation

# Returns
- `old_molecule::Molecule{Frac}`: a copy of the original molecule
"""
function random_reinsertion!(molecule::Molecule{Frac}, box::Box)
    # store old molecule and return at the end for possible restoration
    old_molecule = deepcopy(molecule)

    # translate molecule to a new center of mass
    translate_to!(molecule, Frac(rand(3)))

    # conduct a rotation
    if needs_rotations(molecule)
        random_rotation!(molecule, box)
    end

    # no need to apply BCs b/c by construction we inserted it in the sim box.

    return old_molecule # in case we need to restore
end
