const R_OVERLAP_squared = 0.0001 # Units: Angstrom²
#TODO Keep consistant with `check_for_atom_overlap` in src/Crystal.jl? `check_for_atom_overlap` uses the threshold 0.1 Angstrom (0.01 Angstrom²). Which one to use?

"""
	V = lennard_jones(r_squared::Float64, σ_squared::Float64, ϵ::Float64)  (units: Kelvin)

Calculate the lennard jones potential energy given a radius r between two molecules.
σ and ϵ are specific to interaction between two elements.
returns potential energy in units Kelvin.

# Arguments
- `r²::Float64`: distance between two (pseudo)atoms in question squared (Angstrom²)
- `σ²::Float64`: sigma parameter in Lennard Jones potential squared (units: Angstrom²)
- `ϵ::Float64`: epsilon parameter in Lennard Jones potential (units: Kelvin)
"""
function lennard_jones(r²::Float64, σ²::Float64, ϵ::Float64)
	ratio = (σ² / r²) ^ 3
	return 4.0 * ϵ * (ratio ^ 2 - ratio)
end


"""
    V = vdw_energy(framework::Framework, molecule::Molecule, ljforcefield::LennardJonesForceField, repfactors::Tuple{Int64, Int64, Int64})

Calculates the van der Waals interaction energy between a molecule and a framework.
Applies the nearest image convention to find the closest replicate of a specific atom.

# Arguments
- `framework::Framework`: Crystal structure
- `molecule::Molecule`: adsorbate (includes position/orientation/atoms)
- `ljforcefield::LennardJonesForceField`: Lennard Jones force field
- `repfactors::Tuple{Int64, Int64, Int64}`: replication factors of the home unit cell to build
the supercell such that the nearest image convention can be applied in this function.
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
