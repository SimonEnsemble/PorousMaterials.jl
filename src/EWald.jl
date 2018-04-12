# Vacuum permittivity eps0 = 8.854187817e-12 # C^2/(J-m)
# 1 m = 1e10 A, 1 e = 1.602176e-19 C, kb = 1.3806488e-23 J/K
# 8.854187817e-12 C^2/(J-m) [1 m / 1e10 A] [1 e / 1.602176e-19 C]^2 [kb = 1.3806488e-23 J/K]
const ϵ₀ = 4.7622424954949676e-7  # \epsilon_0 vacuum permittivity units: electron charge^2 /(A - K)
using SpecialFunctions # for erfc

# TODO use NIST test data https://www.nist.gov/mml/csd/chemical-informatics-research-group/spce-water-reference-calculations-10%C3%A5-cutoff
#   to untar: tar zxvf spce_sample_configurations-tar.gz

struct EWaldParams
    α::Float64
    k_repfactors::Tuple{Int, Int, Int}
    sr_cutoff_radius::Float64
    sr_repfactors::Tuple{Int, Int, Int}
end

"""
# Arguments
- `framework::Framework`: Crystal structure (includes charge)
- `molecule::Molecule`: adsorbate (includes posotion/orientation and charge)
- `cutoffradius_squared::Float64`: cutoff radius squared
- `repfactors::Tuple{Int64, Int64, Int64}`: The replication factors used to form the supercell
- `α::Float64`: Ewald sum convergence parameter
- `cutoff_k::Int64`: The cutoff value for each k-vector (ku, kv, kw)
"""
function electrostatic_potential(framework::Framework, x::Array{Float64, 1}, ewald_params::EWaldParams)
    # fractional coordinate (wrap to [0, 1.0])
    xf = mod.(framework.box.c_to_f * x, 1.0)
	
    #
    # Short-range contribution
    #
	sr_potential = 0.0
    for ra = 0:(ewald_params.sr_repfactors[1] - 1), rb = 0:(ewald_params.sr_repfactors[2] - 1), rc = 0:(ewald_params.sr_repfactors[3] - 1)
        for a = 1:framework.n_atoms
            # distance in fractional coords
            dx = xf - (framework.xf[:, a] + [ra, rb, rc])
            # apply nearest image convention for periodic BCs
            nearest_image!(dx, ewald_params.sr_repfactors)
            # convert distance vector to Cartesian coordinates
            dx = framework.box.f_to_c * dx
            r = norm(dx)
            if r < ewald_params.sr_cutoff_radius
                sr_potential += framework.charges[a] / r * erfc(r * ewald_params.α)
            end
        end
    end
    sr_potential /= 4.0 * π * ϵ₀

    #
    # Long-range contribution
    #
	lr_potential = 0.0
    for ka = 0:ewald_params.k_repfactors[1], kb = 0:ewald_params.k_repfactors[2], kc=0:ewald_params.k_repfactors[3]
		if (ka == 0) && (kb == 0) && (kc == 0)
			continue
		end
        # reciprocal vector, k
        k = framework.reciprocal_lattice * [ka, kb, kc]
        # |k|²
        mag_k² = dot(k, k)

        lr_potential_this_k = 0.0
        for a = 1:framework.n_atoms
            dx = framework.f_to_c * (xf - framework.xf[:, a]) # TODO can precompute this...
            # k ⋅ (x - x_a)
            k_dot_dx = dot(k, dx)
            lr_potential_this_k += framework.charges[a] * cos(k_dot_dx)
        end
        lr_potential += 2.0 * lr_potential_this_k / mag_k² * exp(- mag_k² / (4.0 * ewald_params.α ^ 2)) # TODO precompute 4.0 α^2
    end
    # TODO this Ω must be consistent with whatehver reciprocal vectors we use.
    lr_potential /= ϵ₀ * framework.Ω
    return lr_potential + sr_potential
end		
