# δ is half the maximal distance a particle is perturbed in a given coordinate
#  during particle translations
const δ = 2.0 # Å

needs_rotations(molecule::Molecule) = molecule.atoms.n + molecule.charges.n > 1

# break collection of statistics into blocks to gauge convergence and compute standard err
const N_BLOCKS = 5

function random_insertion!(molecules::Array{Molecule{Frac}, 1}, box::Box, template::Molecule{Cart})
    molecule = deepcopy(template)

    if needs_rotations(molecule)
        random_rotation!(molecule)
    end
    molecule = Frac(molecule, box)
    com = Frac(rand(3))
    translate_to!(molecule, com)
    push!(molecules, molecule)
end

function remove_molecule!(molecule_id::Int, molecules::Array{<:Molecule, 1})
    splice!(molecules, molecule_id)
end

# based on center of mass
function apply_periodic_boundary_condition!(molecule::Molecule{Frac})
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

    @assert inside(new_com)

    # translate molecule to new center of mass if it was found to be outside of the box
    translate_to!(molecule, new_com)
end

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
