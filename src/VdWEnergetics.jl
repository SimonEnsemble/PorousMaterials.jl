const R_OVERLAP_squared = 0.0001 # Units: Angstrom²
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
    energy = vdw_energy(framework, molecule, ljforcefield, repfactors)

Calculates the van der Waals interaction energy between a molecule and a framework.
Applies the nearest image convention to find the closest replicate of a specific atom.

# Arguments
- `framework::Framework`: Crystal structure
- `molecule::Molecule`: adsorbate (includes position/orientation/atoms)
- `ljforcefield::LennardJonesForceField`: Lennard Jones force field
- `repfactors::Tuple{Int, Int, Int}`: replication factors of the home unit cell to build
the supercell such that the nearest image convention can be applied in this function.

# Returns
- `energy::Float64`: Van der Waals interaction energy
"""
function vdw_energy(framework::Framework, molecule::Molecule,
                    ljforcefield::LennardJonesForceField, repfactors::Tuple{Int, Int, Int})
	energy = 0.0
    for ljsphere in molecule.ljspheres
        energy += vdw_energy(framework, ljsphere, ljforcefield, repfactors)
    end
	return energy
end

function vdw_energy(framework::Framework, ljsphere::LennardJonesSphere,
                    ljforcefield::LennardJonesForceField, repfactors::Tuple{Int, Int, Int})
	energy = 0.0
    # compute fractional coordinate of the Lennard-Jones sphere and apply PBC to bring it into the
    #  home unit cell of the crystal
    xf_molecule = mod.(framework.box.c_to_f * ljsphere.x, 1.0)

    # loop over replications of the home unit cell to build the supercell
	for ra = 0:1.0:(repfactors[1] - 1), rb = 0:1.0:(repfactors[2] - 1), rc = 0:1.0:(repfactors[3] - 1)
        # distance in fractional coordinate space. same as xf_molecule - (framework.xf + [ra rb rc])
        dxf = broadcast(-, xf_molecule - [ra, rb, rc], framework.xf)

        nearest_image!(dxf, repfactors)
            
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
	end # repfactor end
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
                nearest_image!(dxf, (1, 1, 1))

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
    total_energy = total_guest_host_vdw_energy(framework, molecules, ljforcefield, repfactors)

Compute total guest-host interaction energy (sum over all adsorbates).

# Arguments
- `framework::Framework`: The framework containing the crystal structure information
- `molecules::Array{Molecule, 1}`: An array of Molecule data structures
- `ljforcefield::LennardJonesForceField`: A Lennard Jones forcefield data structure describing the interactions between different atoms
- `repfactors::Tuple{Int, Int, Int}`: The replication factors we use to replicate our unit cell

# Returns
- `total_energy::Float64`: The total guest-host van der Waals energy
"""
function total_guest_host_vdw_energy(framework::Framework,
                                     molecules::Array{Molecule, 1},
                                     ljforcefield::LennardJonesForceField,
                                     repfactors::Tuple{Int, Int, Int})
    total_energy = 0.0
    for molecule in molecules
        total_energy += vdw_energy(framework, molecule, ljforcefield, repfactors)
    end
    return total_energy
end

"""
    total_energy = total_guest_guest_vdw_energy(molecules, ljforcefield, simulation_box)

Compute sum of all guest-guest interaction energy from vdW interactions.

# Arguments
- `molecules::Array{Molecule, 1}`: An array of Molecule data structures
- `ljforcefield::LennardJonesForceField`: A Lennard Jones forcefield data structure describing the interactions between different atoms
- `simulation_box::Box`: The simulation box for the computation.

# Returns
- `total_energy::Float64`: The total guest-guest van der Waals energy
"""
function total_guest_guest_vdw_energy(molecules::Array{Molecule, 1},
                                      ljforcefield::LennardJonesForceField,
                                      simulation_box::Box)
    total_energy = 0.0
    for molecule_id = 1:length(molecules)
        total_energy += vdw_energy(molecule_id, molecules, ljforcefield, simulation_box)
    end
    return total_energy / 2.0 # avoid double-counting pairs
end
