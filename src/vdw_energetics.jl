# atoms are considered to overlap and thus have Inf energy if distance squared is less than this.
const R²_OVERLAP = 0.1 # Units: Angstrom²

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
    return 4.0 * ϵ * ratio * (ratio - 1.0)
end

# generic atoms-atoms van der waals energy function; pass bigger list of atoms second for speed
function vdw_energy(atoms_i::Atoms{Frac}, atoms_j::Atoms{Frac}, box::Box, ljff::LJForceField)
    energy = 0.0
    for i = 1:atoms_i.n
        # vectors from atom i to atoms_j in fractional space
        @inbounds dx = broadcast(-, atoms_j.coords.xf, atoms_i.coords.xf[:, i])
        nearest_image!(dx)
        # convert to Cartesian coords
        @inbounds dx .= box.f_to_c * dx
        for j = 1:atoms_j.n
            @inbounds r² = dx[1, j] ^ 2 + dx[2, j] ^ 2 + dx[3, j] ^ 2
            if r² < ljff.r²_cutoff # within cutoff radius
                if r² < R²_OVERLAP # overlapping atoms
                    return Inf # no sense in continuing to next atom and adding Inf
                else # not overlapping
                    @inbounds energy += lennard_jones(r²,
                             ljff.σ²[atoms_i.species[i]][atoms_j.species[j]],
                             ljff.ϵ[ atoms_i.species[i]][atoms_j.species[j]])
                end
            end
        end
    end
    return energy
end

vdw_energy(atoms_i::Atoms{Cart}, atoms_j::Atoms{Cart}, box::Box, ljff::LJForceField) = vdw_energy(Frac(atoms_i, box), Frac(atoms_j, box), box, ljff)
vdw_energy(atoms_i::Atoms{Frac}, atoms_j::Atoms{Cart}, box::Box, ljff::LJForceField) = vdw_energy(atoms_i, Frac(atoms_j, box), box, ljff)
vdw_energy(atoms_i::Atoms{Cart}, atoms_j::Atoms{Frac}, box::Box, ljff::LJForceField) = vdw_energy(Frac(atoms_i, box), atoms_j, box, ljff)

"""
   energy = vdw_energy(crystal, molecule, ljforcefield)

Calculates the van der Waals interaction energy between a molecule and a crystal.
Applies the nearest image convention to find the closest replicate of a specific atom.

WARNING: it is assumed that the framework is replicated sufficiently such that the nearest
image convention can be applied. See [`replicate`](@ref) and [`replication_factors`](@ref).
"""
function vdw_energy(crystal::Crystal, molecule::Molecule, ljff::LJForceField)
   return vdw_energy(molecule.atoms, crystal.atoms, crystal.box, ljff)
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
function vdw_energy(molecule_id::Int, molecules::Array{<:Molecule, 1}, ljff::LJForceField, box::Box)
   energy = 0.0
   # loop over all other molecule
   for other_molecule_id = 1:length(molecules)
       # molecule cannot interact with itself
       if other_molecule_id == molecule_id
           continue
       end
       energy += vdw_energy(molecules[molecule_id].atoms, molecules[other_molecule_id].atoms, box, ljff)
   end
   return energy # units are the same as in ϵ for forcefield (Kelvin)
end

###
#  for mixture simulations.
#  compute potential energy of molecules[which_spcies][molecule_id] with all other molecules.
###
function vdw_energy(species_id::Int, molecule_id::Int, molecules::Array{Array{Molecule{Frac}, 1}, 1}, ljff::LJForceField, box::Box)
    energy = 0.0
    # loop over species
    for s in 1:length(molecules)
        # loop over molecules of this species, in molecules[s]
        for m in 1:length(molecules[s])
            # molecule cannot interact with itself (only relevant for molecules of the same species)
            if (s == species_id) && (m == molecule_id)
                continue
            end

            energy += vdw_energy(molecules[species_id][molecule_id].atoms, molecules[s][m].atoms, box, ljff)
        end
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
function total_vdw_energy(crystal::Crystal, molecules::Array{<:Molecule, 1}, ljff::LJForceField)
   total_energy = 0.0
   for molecule in molecules
       total_energy += vdw_energy(crystal, molecule, ljff)
   end
   return total_energy
end

function total_vdw_energy(molecules::Array{<:Molecule, 1}, ljff::LJForceField, box::Box)
   total_energy = 0.0
   for i = 1:length(molecules)
       total_energy += vdw_energy(i, molecules, ljff, box)
   end
   return total_energy / 2.0 # avoid double-counting pairs
end

"""
   pot_energy = vdw_energy_no_PBC(atoms_i, atoms_j , ljff)

compute vdw potential energy without periodic boundary conditions
"""
function vdw_energy_no_PBC(atoms_i::Atoms{Cart}, atoms_j::Atoms{Cart}, ljff::LJForceField)
    energy = 0.0
    for i = 1:atoms_i.n
        # vectors from atom i to atoms_j in fractional space
        @inbounds dx = broadcast(-, atoms_j.coords.x, atoms_i.coords.x[:, i])
        for j = 1:atoms_j.n
            @inbounds r² = dx[1, j] ^ 2 + dx[2, j] ^ 2 + dx[3, j] ^ 2
            if r² < ljff.r²_cutoff # within cutoff radius
                if r² < R²_OVERLAP # overlapping atoms
                    return Inf # no sense in continuing to next atom and adding Inf
                else # not overlapping
                    @inbounds energy += lennard_jones(r²,
                             ljff.σ²[atoms_i.species[i]][atoms_j.species[j]],
                             ljff.ϵ[ atoms_i.species[i]][atoms_j.species[j]])
                end
            end
        end
    end
    return energy
end

function total_vdw_energy(crystal::Crystal, molecules::Array{Array{Molecule{Frac}, 1}, 1}, ljff::LJForceField)
    total_energy = 0.0
    for molecules_species in molecules
        total_energy += total_vdw_energy(crystal, molecules_species, ljff)
    end
    return total_energy
end

function total_vdw_energy(molecules::Array{Array{Molecule{Frac}, 1}, 1}, ljff::LJForceField, box::Box)
    total_energy = 0.0 
    for (species_id, molecules_species) in enumerate(molecules)
        for molecule_id in 1:length(molecules_species)
            total_energy += vdw_energy(species_id, molecule_id, molecules, ljff, box)
        end
    end
    return total_energy / 2.0 # avoid double-counting pairs 
end
