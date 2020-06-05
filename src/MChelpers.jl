# δ is half the maximum distance a particle can be displaced per dimension during a MC translation
# move. δ is adaptively adjusted, but the range that can be adjusted is bounded by maxδ and minδ.
const maxδ = 2.0 # Å
const minδ = 0.01 # Å

# 2π*φ is the  maximum a particle can be rotated during a MC rotation move. If φ = 1, then the
# maximum rotation is 2π hence any rotation is valid.  For φ < 1.0, the rotations are limited to
# a range near the starting point, e.g. if φ = 0.5, the possible rotations would represent one
# half the rotation sphere. φ is adaptively adjusted, but the range that can be adjusted is bounded
# by maxφ and minφ.
const maxφ = 1.0 # all rotations on sphere valid; equivalent to 2π radians
const minφ = 0.0001 # move only 0.01% of full spherical rotation; equivalent to 2π/1000 radians

# δ and φ will be adaptively adjusted until the actual acceptance rate equals the desired acceptance
# rate, or the boundaries to δ and φ defined above are reached.
const desired_acceptance_rate = 0.50

# break collection of statistics into blocks to gauge convergence and compute standard err
const N_BLOCKS = 5

"""
    insert_molecule!(molecules, box, template)

Inserts an additional adsorbate molecule into the simulation box using the template provided.
The center of mass of the molecule is chosen at a uniform random position in the simulation box.
A uniformly random orientation of the molecule is chosen by rotating about the center of mass.

# Arguments
- `molecules::Array{Molecule, 1}`: An array of Molecule objects
- `box::Box`: The simulation box
- `template::Molecule`: A template molecule used as reference when inserting molecules
"""
function insert_molecule!(molecules::Array{Molecule, 1}, box::Box, template::Molecule)
    # choose center of mass
    xf = rand(3)
    # copy the template
    molecule = deepcopy(template)
    # conduct a rotation if needed
    if rotatable(molecule)
        rotate!(molecule, box)
    end
    # translate molecule to its new center of mass
    translate_to!(molecule, xf)
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
    apply_periodic_boundary_condition!(molecule)

Check if the `center_of_mass` of a `Molecule` is outside of a `Box`. If so, apply periodic
boundary conditions and translate the center of mass of the `Molecule` (and its atoms
and point charges) so that it is inside of the `Box`.

# Arguments
- `molecule::Molecule`: A molecule we're interested in seeing if its' center of mass falls within `simulation_box`
"""
function apply_periodic_boundary_condition!(molecule::Molecule)
    outside_box = false # do nothing if not outside the box

    # current center of mass in fractional coordinates; adjust inside loop
    xf = deepcopy(molecule.xf_com)

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
        translate_to!(molecule, xf)
    end
end


"""
    new_δ = adjust_δ(δ, accepted, proposed)

Compares the current acceptance rate (accepted / proposed) to the desired acceptance rate and then
adjusts δ up or down to get closer to the desired acceptance rate. This new δ is bounded by the
maxδ and minδ consts before being returned.

# Arguments
- `δ::Float64`: current value of the translation scale parameter δ.
- `accepted::Int`: number of translation moves accepted.
- `proposed::Int`: number of translation moves proposed.
"""
function adjust_δ(δ::Float64, accepted::Int, proposed::Int)
    recent_acceptance_rate = accepted / proposed
    new_δ = max(min(δ * recent_acceptance_rate / desired_acceptance_rate, maxδ), minδ)
    @printf("adjusting: δ = %f => %f [acceptance rate: %f]\n", δ, new_δ, recent_acceptance_rate)
    return new_δ
end

"""
    new_φ = adjust_φ(φ, accepted, proposed)

Compares the current acceptance rate (accepted / proposed) to the desired acceptance rate and then
adjusts φ up or down to get closer to the desired acceptance rate. This new φ is bounded by the
maxφ and minφ consts before being returned.

# Arguments
- `φ::Float64`: current value of the rotation scale parameter φ.
- `accepted::Int`: number of rotation moves accepted.
- `proposed::Int`: number of rotation moves proposed.
"""
function adjust_φ(φ::Float64, accepted, proposed)
    recent_acceptance_rate = accepted / proposed
    new_φ = max(min(φ * recent_acceptance_rate / desired_acceptance_rate, maxφ), minφ)
    @printf("adjusting: φ = %f => %f [acceptance rate: %f]\n", φ, new_φ, recent_acceptance_rate)
    return new_φ
end

"""
    translate_molecule!(molecule, box)

Perturbs the Cartesian coordinates of a molecule about its center of mass by a random
vector of max length δ. Applies periodic boundary conditions to keep the molecule inside
the simulation box. Returns a deep copy of the old molecule in case it needs replaced
if the Monte Carlo proposal is rejected.

# Arguments
- `molecule::Molecule`: The molecule we want to perturb
- `box::Box`: The simulation box
- `scale::Float64`: the translation scale factor δ.

# Returns
- `old_molecule::Molecule`: The old molecule in case the MC proposal is rejected
"""
function translate_molecule!(molecule::Molecule, box::Box; scale::Float64)
    # store old molecule and return at the end for possible restoration
    old_molecule = deepcopy(molecule)

    # perturb in Cartesian coords in a random cube centered at current coords.
    dx = scale * (rand(3) .- 0.5) # move every atom of the molecule by the same vector.
    translate_by!(molecule, dx, box)

    # done, unless the molecule has moved outside of the box, then apply PBC
    apply_periodic_boundary_condition!(molecule)

    return old_molecule # in case we need to restore coords
end

"""
    reinsert_molecule(molecule, box)

Move molecule to a new center of mass randomly distrubted in the unit cell and choose
a random orientation for it. Return a deep copy of the starting molecule for possible
restoration. This MC move can be viewed as a more aggressive `translate_molecule!`.

# Arguments
- `molecule::Molecule`: The molecule we want to perturb
- `box::Box`: The simulation box
"""
function reinsert_molecule!(molecule::Molecule, box::Box)
    # store old molecule and return at the end for possible restoration
    old_molecule = deepcopy(molecule)

    # translate molecule to a new center of mass
    translate_to!(molecule, rand(3))

    # conduct a rotation
    if rotatable(molecule)
        rotate!(molecule, box)
    end

    # no need to apply BCs b/c by construction we inserted it in the sim box.

    return old_molecule # in case we need to restore
end

# do we need to conduct a rotation or not? # TODO what if it is an ion? No need to rotate...
"""
    need_to_rotate = rotatable(molecule)

Determines whether or not a given molecule needs to be rotated. For example,
rotating a single atom isn't necessary.

# Arguments
- `molecule::Molecule`: The molecule being tested. This function determines if a
    rotation of this molecule will do anything.

# Returns
- `is_rotatable::Bool`: A boolean describing whether or not rotating the molecule
    will alter its interactions with other molecules 
"""
rotatable(molecule::Molecule) = (molecule.atoms.n_atoms + molecule.charges.n_charges > 1)::Bool
