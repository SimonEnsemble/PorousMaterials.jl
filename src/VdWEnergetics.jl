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
@inline function lennard_jones(r²::Float64, σ²::Float64, ϵ::Float64)
	ratio = (σ² / r²) ^ 3
	return 4.0 * ϵ * ratio * (ratio - 1.0)
end

# generic atoms-atoms van der waals energy function; pass bigger list of atoms second for speed
@inline function vdw_energy(atoms_i::Atoms, atoms_j::Atoms, box::Box, ljff::LJForceField)
	energy = 0.0
    for i = 1:atoms_i.n_atoms
        # vectors from atom i to atoms_j in fractional space
        @inbounds dxf = broadcast(-, atoms_j.xf, atoms_i.xf[:, i])
        nearest_image!(dxf)
        # convert to Cartesian coords
        @inbounds dxf .= box.f_to_c * dxf
        for j = 1:atoms_j.n_atoms
            @inbounds r² = dxf[1, j] ^ 2 + dxf[2, j] ^ 2 + dxf[3, j] ^ 2
            if r² < ljff.cutoffradius_squared # within cutoff radius
                if r² < R_OVERLAP_squared # overlapping atoms
                    return Inf # no sense in continuing to next atom and adding Inf
                else # not overlapping
                    @inbounds energy += lennard_jones(r²,
                            ljff.σ²[atoms_i.species[i]][atoms_j.species[j]],
                             ljff.ϵ[atoms_i.species[i]][atoms_j.species[j]])  
                end
            end
        end
    end
	return energy
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
- `ljforcefield::LJForceField`: Lennard Jones force field

# Returns
- `energy::Float64`: Van der Waals interaction energy
"""
@inline function vdw_energy(framework::Framework, molecule::Molecule, ljff::LJForceField)
    return vdw_energy(molecule.atoms, framework.atoms, framework.box, ljff)
end

"""
    gg_energy = vdw_energy(molecule_id, molecules, ljforcefield, simulation_box)

Calculates van der Waals interaction energy of a single adsorbate `molecules[molecule_id]`
with all of the other molecules in the system. Periodic boundary conditions are applied,
using the nearest image convention.

# Arguments
- `molecule_id::Int`: Molecule ID used to determine which molecule in `molecules` we wish to calculate the guest-guest interactions
- `molecules::Array{Molecule, 1}`: An array of Molecule data structures
- `ljforcefield::LJForceField`: A Lennard Jones forcefield data structure describing the interactions between different atoms
- `simulation_box::Box`: The simulation box for the computation.

# Returns
- `gg_energy::Float64`: The guest-guest interaction energy of `molecules[molecule_id]` with the other molecules in `molecules`
"""
function vdw_energy(molecule_id::Int, molecules::Array{Molecule, 1}, ljff::LJForceField, box::Box)
    energy = 0.0
    # loop over all other molecules
    for other_molecule_id = 1:length(molecules)
        # molecule cannot interact with itself
        if other_molecule_id == molecule_id
            continue
        end
        @inbounds energy += vdw_energy(molecules[molecule_id].atoms, molecules[other_molecule_id].atoms, box, ljff)
    end
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
- `ljforcefield::LJForceField`: A Lennard Jones forcefield data structure describing the interactions between different atoms
- `simulation_box::Box`: The simulation box for application of PBCs.

# Returns
- `total_energy::Float64`: The total guest-host or guest-guest van der Waals energy
"""
function total_vdw_energy(framework::Framework, molecules::Array{Molecule, 1}, ljff::LJForceField)
    total_energy = 0.0
    for molecule in molecules
        total_energy += vdw_energy(framework, molecule, ljff)
    end
    return total_energy
end

function total_vdw_energy(molecules::Array{Molecule, 1}, ljff::LJForceField, box::Box)
    total_energy = 0.0
    for i = 1:length(molecules)
        total_energy += vdw_energy(i, molecules, ljff, box)
    end
    return total_energy / 2.0 # avoid double-counting pairs
end

"""
    pot_energy = vdw_energy_no_PBC(atoms, molecule, ljff)

Assumes unit cell box is a unit cube and no periodic boundary conditions
are applied.
"""
function vdw_energy_no_PBC(atoms::Atoms, molecule::Molecule, ljff::LJForceField)
    energy = 0.0
    for i = 1:molecule.atoms.n_atoms # loop over all atoms in molecule
        for j = 1:atoms.n_atoms
            dx = molecule.atoms.xf[:, i] - atoms.xf[:, j]
            r² = dx[1] * dx[1] + dx[2] * dx[2] + dx[3] * dx[3]
            energy += lennard_jones(r², 
                ljff.σ²[molecule.atoms.species[i]][atoms.species[j]],
                 ljff.ϵ[molecule.atoms.species[i]][atoms.species[j]])
        end
    end
	return energy
end
