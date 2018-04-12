# Vacuum permittivity eps0 = 8.854187817e-12 # C^2/(J-m)
# 1 m = 1e10 A, 1 e = 1.602176e-19 C, kb = 1.3806488e-23 J/K
# 8.854187817e-12 C^2/(J-m) [1 m / 1e10 A] [1 e / 1.602176e-19 C]^2 [kb = 1.3806488e-23 J/K]
const ϵ₀ = 4.7622424954949676e-7  # \epsilon_0 vacuum permittivity units: electron charge^2 /(A - K)
using SpecialFunctions # for erfc

# TODO use NIST test data https://www.nist.gov/mml/csd/chemical-informatics-research-group/spce-water-reference-calculations-10%C3%A5-cutoff
#   to untar: tar zxvf spce_sample_configurations-tar.gz

 # struct EWaldParams
 #     α::Float64
 #     k_repfactors::Tuple{Int, Int, Int}
 #     sr_cutoff_radius::Float64
 # end

"""
# Arguments
- `framework::Framework`: Crystal structure (includes charge)
- `molecule::Molecule`: adsorbate (includes posotion/orientation and charge)
- `cutoffradius_squared::Float64`: cutoff radius squared
- `repfactors::Tuple{Int64, Int64, Int64}`: The replication factors used to form the supercell
- `α::Float64`: Ewald sum convergence parameter
- `cutoff_k::Int64`: The cutoff value for each k-vector (ku, kv, kw)
"""
function electrostatic_potential(framework::Framework, x::Array{Float64, 1}, sim_box::Box, rep_factors::Tuple{Int, Int, Int}, k_rep_factors::Tuple{Int, Int, Int}, α::Float64)
    # fractional coordinate of point at which we're computing the electrostatic potential (wrap to [0, 1.0])
    xf = mod.(framework.box.c_to_f * x, 1.0)
	
    #
    # Short-range contribution
    #
	sr_potential = 0.0 # short-range contribution
	lr_potential = 0.0 # long-range contribution
    # build the super cell with loop over (ra, rb, rc) and atoms in the primitive unit cell
    for ra = 0:(rep_factors[1] - 1), rb = 0:(rep_factors[2] - 1), rc = 0:(rep_factors[3] - 1)
        for i = 1:framework.n_atoms
            # x - x_i vector in fractional coords
            dx = xf - (framework.xf[:, i] + [ra, rb, rc])

            ###
            #  Long-range contribution
            ###
            for ka = 0:k_rep_factors[1], kb = 0:k_rep_factors[2], kc=0:k_rep_factors[3]
                if (ka == 0) && (kb == 0) && (kc == 0)
                    continue
                end

                # reciprocal vector, k
                k = sim_box.reciprocal_lattice * [ka, kb, kc]
                # |k|²
                mag_k² = dot(k, k)
                # k ⋅ (x - x_a)
                k_dot_dx = dot(k, dx)
                lr_potential += 2.0 * framework.charges[i] * cos(k_dot_dx) / mag_k² * exp(- mag_k² / (4.0 * α ^ 2)) # TODO precompute some of this.
            end

            ###
            #  Short range contribution
            ###
            # apply nearest image convention for periodic BCs
            nearest_image!(dx, rep_factors)

            # convert distance vector to Cartesian coordinates
            dx = framework.box.f_to_c * dx

            r = norm(dx)

            if r < sr_cutoff_radius
                sr_potential += framework.charges[i] / r * erfc(r * α)
            end
        end
    end
    sr_potential /= 4.0 * π * ϵ₀
    lr_potential /= ϵ₀ * sim_box.Ω

    return lr_potential + sr_potential
end		
