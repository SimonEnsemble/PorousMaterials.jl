const R_OVERLAP_squared = 0.1 # Units: Angstrom²
#TODO Keep consistant with `check_for_atom_overlap` in src/Crystal.jl? `check_for_atom_overlap` uses the threshold 0.1 Angstrom (0.01 Angstrom²). Which one to use?

"""
	energy = lennard_jones(r², σ², ϵ)  (units: Kelvin)

Calculate the lennard jones potential energy given the square of the radius r between two
lennard-jones spheres. σ and ϵ are specific to interaction between two elements. Return
the potential energy in units Kelvin (well, whatever the units of ϵ are).

# Arguments
- `r²::Float64`: distance between two (pseudo)atoms in question squared (Angstrom²)
- `σ²::Float64`: sigma parameter in Lennard Jones potential squared (units: Angstrom²)
- `ϵ::Float64`: epsilon parameter in Lennard Jones potential (units: Kelvin)

# Returns
- `energy::Float64`: Lennard Jones potential energy
"""
function lennard_jones(r²::Float64, σ²::Float64, ϵ::Float64)
	ratio = (σ² / r²) ^ 3
	return 4.0 * ϵ * (ratio ^ 2 - ratio)
end


"""
    energy = vdw_energy(framework, molecule, ljforcefield)

Calculates the van der Waals interaction energy between a molecule and a framework.
Applies the nearest image convention to find the closest replicate of a specific atom.

WARNING: it is assumed that the framework is replicated sufficiently such that the nearest
image convention can be applied. See [`replicate`](@ref).

# Arguments
- `framework::Framework`: Crystal structure
- `molecule::Molecule`: adsorbate (includes position/orientation/atoms)
- `ljforcefield::LennardJonesForceField`: Lennard Jones force field

# Returns
- `energy::Float64`: Van der Waals interaction energy
"""
function vdw_energy(framework::Framework, molecule::Molecule,
                    ljforcefield::LennardJonesForceField)
	energy = 0.0
    for ljsphere in molecule.ljspheres
        energy += vdw_energy(framework, ljsphere, ljforcefield)
    end
	return energy
end

function vdw_energy(framework::Framework, ljsphere::LennardJonesSphere,
                    ljforcefield::LennardJonesForceField)
	energy = 0.0
    # compute fractional coordinate of the Lennard-Jones sphere and apply PBC to bring it into the
    #  home unit cell of the crystal
    xf_molecule = mod.(framework.box.c_to_f * ljsphere.x, 1.0)

    # loop over replications of the home unit cell to build the supercell
    # distance in fractional coordinate space.
    dxf = broadcast(-, xf_molecule, framework.xf)

    nearest_image!(dxf)

    # Distance in cartesian coordinate space
    dx = framework.box.f_to_c * dxf

    # loop over atoms of the framework and compute its contribution to the vdW energy.
    for i = 1:framework.n_atoms
        @inbounds r² = dx[1, i] * dx[1, i] + dx[2, i] * dx[2, i] + dx[3, i] * dx[3, i]

        if r² < R_OVERLAP_squared
            # if adsorbate atom overlaps with an atom, return Inf (R_OVERLAP is defined as 0.01 Angstrom, or `R_OVERLAP_squared = 0.0001 Angstrom²)
            return Inf
        elseif r² < ljforcefield.cutoffradius_squared
            # add pairwise contribution to potential energy
            @inbounds energy += lennard_jones(r²,
                ljforcefield.σ²[framework.atoms[i]][ljsphere.atom],
                ljforcefield.ϵ[framework.atoms[i]][ljsphere.atom])
        end # if-elseif-end
    end # framework atoms end
	return energy
end

"""
    gg_energy = vdw_energy(molecule_id, molecules, ljforcefield, simulation_box)

Calculates van der Waals interaction energy of a single adsorbate `molecules[molecule_id]`
with all of the other molecules in the system. Periodic boundary conditions are applied,
using the nearest image convention.

# Arguments
- `molecule_id::Int`: Molecule ID used to determine which molecule in `molecules` we wish to calculate the guest-guest interactions
- `molecules::Array{Molecule, 1}`: An array of Molecule data structures
- `ljforcefield::LennardJonesForceField`: A Lennard Jones forcefield data structure describing the interactions between different atoms
- `simulation_box::Box`: The simulation box for the computation.

# Returns
- `gg_energy::Float64`: The guest-guest interaction energy of `molecules[molecule_id]` with the other molecules in `molecules`
"""
function vdw_energy(molecule_id::Int, molecules::Array{Molecule, 1}, ljforcefield::LennardJonesForceField, simulation_box::Box)
    energy = 0.0 # energy is pair-wise additive
    # Look at interaction with all other molecules in the system
    for this_ljsphere in molecules[molecule_id].ljspheres
        # Loop over all atoms in the given molecule
        for other_molecule_id = 1:length(molecules)
            # molecule cannot interact with itself
            if other_molecule_id == molecule_id
                continue
            end
            # loop over every ljsphere (atom) in the other molecule
            for other_ljsphere in molecules[other_molecule_id].ljspheres
                # compute vector between molecules in fractional coordinates
                dxf = simulation_box.c_to_f * (this_ljsphere.x - other_ljsphere.x)

                # simulation box has fractional coords [0, 1] by construction
                nearest_image!(dxf)

                # converts fractional distance to cartesian distance
                dx = simulation_box.f_to_c * dxf

                r² = dx[1] * dx[1] + dx[2] * dx[2] + dx[3] * dx[3]

                if r² < R_OVERLAP_squared
                    return Inf
                elseif r² < ljforcefield.cutoffradius_squared
                    energy += lennard_jones(r²,
                        ljforcefield.σ²[this_ljsphere.atom][other_ljsphere.atom],
                        ljforcefield.ϵ[this_ljsphere.atom][other_ljsphere.atom])
                end
            end # loop over all ljspheres in other molecule
        end # loop over all other molecules
    end # loop over all ljspheres in this molecule
    return energy # units are the same as in ϵ for forcefield (Kelvin)
end

"""
    total_gh_energy = total_vdw_energy(framework, molecules, ljforcefield) # guest-host
    total_gg_energy = total_vdw_energy(molecules, ljforcefield, simulation_box) # guest-guest

Compute total guest-host (gh) or guest-guest (gg) interaction energy, i.e. the contribution
from all adsorbates in `molecules`.

WARNING: it is assumed that the framework is replicated sufficiently such that the nearest
image convention can be applied. See [`replicate`](@ref).

# Arguments
- `framework::Framework`: The framework containing the crystal structure information
- `molecules::Array{Molecule, 1}`: An array of Molecule data structures
- `ljforcefield::LennardJonesForceField`: A Lennard Jones forcefield data structure describing the interactions between different atoms
- `simulation_box::Box`: The simulation box for application of PBCs.

# Returns
- `total_energy::Float64`: The total guest-host or guest-guest van der Waals energy
"""
function total_vdw_energy(framework::Framework,
                          molecules::Array{Molecule, 1},
                          ljforcefield::LennardJonesForceField)
    total_energy = 0.0
    for molecule in molecules
        total_energy += vdw_energy(framework, molecule, ljforcefield)
    end
    return total_energy
end

function total_vdw_energy(molecules::Array{Molecule, 1},
                         ljforcefield::LennardJonesForceField,
                         simulation_box::Box)
    total_energy = 0.0
    for molecule_id = 1:length(molecules)
        total_energy += vdw_energy(molecule_id, molecules, ljforcefield, simulation_box)
    end
    return total_energy / 2.0 # avoid double-counting pairs
end

"""
    energy = vdw_energy_no_PBC(molecule, atoms, x, ljff)

Calculates the van der Waals interaction energy between a molecule and a list of `atoms`
at Cartesian positions `x` using a lennard jones force field `ljff`
No periodic boundary conditions are applied.
"""
function vdw_energy_no_PBC(molecule::Molecule, atoms::Array{Symbol, 1}, x::Array{Float64, 2},
                           ljforcefield::LennardJonesForceField)
	energy = 0.0
    # loop over lennard-jones spheres in the molecule
    for ljs in molecule.ljspheres
        # loop over atoms
        for i = 1:length(atoms)
            dx = x[:, i] - ljs.x
            r² = dx[1] ^ 2 + dx[2] ^ 2 + dx[3] ^ 2

            if r² < R_OVERLAP_squared
                return Inf
            elseif r² < ljforcefield.cutoffradius_squared
                # add pairwise contribution to potential energy
                energy += lennard_jones(r²,
                    ljforcefield.σ²[atoms[i]][ljs.atom],
                    ljforcefield.ϵ[atoms[i]][ljs.atom])
            end
        end # loop over atoms
    end # loop over ljspheres
	return energy
end
