# δ is half the maximal distance a particle is perturbed in a given coordinate
#  during particle translations
const δ = 2.0 # Å

"""
    insert_molecule!(molecules, simulation_box, template)

Inserts an additional adsorbate molecule into the simulation box using the template provided.
The center of mass of the molecule is chosen at a uniform random position in the simulation box.
A uniformly random orientation of the molecule is chosen by rotating about the center of mass.

# Arguments
- `molecules::Array{Molecule, 1}`: An array of Molecule objects
- `simulation_box::Box`: The simulation box
- `template::Molecule`: A template molecule used as reference when inserting molecules
"""
function insert_molecule!(molecules::Array{Molecule, 1}, simulation_box::Box, template::Molecule)
    # choose center of mass
    x = simulation_box.f_to_c * rand(3)
    # copy the template
    molecule = deepcopy(template)
    # conduct a rotation
    if (length(molecule.ljspheres) + length(molecule.charges) > 1)
        rotate!(molecule)
    end
    # translate molecule to its new center of mass
    translate_to!(molecule, x)
    # push molecule to array.
    push!(molecules, molecule)
end

"""
    delete_molecule!(molecule_id, molecules)

Removes a random molecule from the current molecules in the framework.
molecule_id decides which molecule will be deleted, for a simulation, it must
be a randomly generated value

# Arguments
- `molecule_id::Int`: The molecule ID is used to determine which molecule in `molecules` should be removed
- `molecules::Array{Molecule, 1}`: An array of Molecule objects
"""
function delete_molecule!(molecule_id::Int, molecules::Array{Molecule, 1})
    splice!(molecules, molecule_id)
end

"""
    apply_periodic_boundary_condition!(molecule, simulation_box)

Check if the `center_of_mass` of a `Molecule` is outside of a `Box`. If so, apply periodic
boundary conditions and translate the center of mass of the `Molecule` (and its atoms
and point charges) so that it is inside of the `Box`.

# Arguments
- `molecule::Molecule`: A molecule we're interested in seeing if its' center of mass falls within `simulation_box`
- `simulation_box::Box`: The simulation box
"""
function apply_periodic_boundary_condition!(molecule::Molecule, box::Box)
    outside_box = false # do nothing if not outside the box

    # compute its center of mass in fractional coordinates
    xf = box.c_to_f * molecule.center_of_mass

    # apply periodic boundary conditions
    for k = 1:3 # loop over xf, yf, zf components
        # if > 1.0, shift down
        if xf[k] >= 1.0
            outside_box = true
            xf[k] -= 1.0
        elseif xf[k] < 0.0
            outside_box = true
            xf[k] += 1.0
        end
    end

    # translate molecule to new center of mass if it was found to be outside of the box
    if outside_box
        new_center_of_mass = box.f_to_c * xf
        translate_to!(molecule, new_center_of_mass)
    end
end

"""
    translate_molecule!(molecule, simulation_box)

Perturbs the Cartesian coordinates of a molecule about its center of mass by a random
vector of max length δ. Applies periodic boundary conditions to keep the molecule inside
the simulation box. Returns a deep copy of the old molecule in case it needs replaced
if the Monte Carlo proposal is rejected.

# Arguments
- `molecule::Molecule`: The molecule we want to perturb
- `simulation_box::Box`: The simulation box

# Returns
- `old_molecule::Molecule`: The old molecule in case the MC proposal is rejected
"""
function translate_molecule!(molecule::Molecule, simulation_box::Box)
    # store old molecule and return at the end for possible restoration
    old_molecule = deepcopy(molecule)
    # peturb in Cartesian coords in a random cube centered at current coords.
    dx = δ * (rand(3) - 0.5) # move every atom of the molecule by the same vector.
    translate_by!(molecule, dx)
    # done, unless the molecule has moved outside of the box, then apply PBC
    apply_periodic_boundary_condition!(molecule, simulation_box)

    return old_molecule # in case we need to restore
end

"""
    reinsert_molecule(molecule, simulation_box)

Move molecule to a new center of mass randomly distrubted in the unit cell and choose
a random orientation for it. Return a deep copy of the starting molecule for possible
restoration. This MC move can be viewed as a more aggressive `translate_molecule!`.
"""
function reinsert_molecule!(molecule::Molecule, simulation_box::Box)
    # store old molecule and return at the end for possible restoration
    old_molecule = deepcopy(molecule)

    # choose new center of mass
    x = simulation_box.f_to_c * rand(3)

    # translate molecule to its new center of mass
    translate_to!(molecule, x)

    # conduct a rotation
    if (length(molecule.ljspheres) + length(molecule.charges) > 1)
        rotate!(molecule)
    end

    # no need to apply BCs b/c by construction we inserted it in the sim box.

    return old_molecule # in case we need to restore
end

# do we need to conduct a rotation or not?
need_rotation(molecule::Molecule) = (length(molecule.ljspheres) + length(molecule.charges) > 1)
