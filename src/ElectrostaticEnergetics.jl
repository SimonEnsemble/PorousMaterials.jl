using SpecialFunctions
using OffsetArrays
# Vacuum permittivity eps0 = 8.854187817e-12 # C^2/(J-m)
# 1 m = 1e10 A, 1 e = 1.602176e-19 C, kb = 1.3806488e-23 J/K
# 8.854187817e-12 C^2/(J-m) [1 m / 1e10 A] [1 e / 1.602176e-19 C]^2 [kb = 1.3806488e-23 J/K]
const ϵ₀ = 4.7622424954949676e-7  # \epsilon_0 vacuum permittivity units: electron charge^2 /(A - K)

# TODO use NIST test data https://www.nist.gov/mml/csd/chemical-informatics-research-group/spce-water-reference-calculations-10%C3%A5-cutoff
#   to untar: tar zxvf spce_sample_configurations-tar.gz
struct Kvector
    k::Array{Float64, 1}
    wt::Float64
end

"""
    kvectors = compute_kvectors(sim_box, k_repfactors, α)

For speed, pre-compute the reciprocal lattice k-vectors for the long-range Ewald summation 
in Fourier space. This function takes advantage of the symmetry:
    cos(-k⋅(x-xᵢ)) + cos(k⋅(x-xᵢ)) = 2 cos(k⋅(x-xᵢ))
to reduce the number of k-vectors in the long-range computation. This reduction is 
accounted for in the `weight` attribute of the `Kvector`. Returns an array of `Kvector`'s.

# Arguments
- `sim_box::Box`: the simulation box that consits of replicated primitive unit cells of the framework
- `k_repfactors::Tuple{Int, Int, Int}`: The number of reciprocal lattice space replications to include in the long-range Ewald sum.
- `α::Float64`: Ewald sum convergence parameter. Units: inverse Å.
"""
function precompute_kvec_wts(sim_box::Box, kreps::Tuple{Int, Int, Int}, α::Float64)
    kvec_wts = OffsetArray(zeros(Float64, kreps[1] + 1, 2 * kreps[2] + 1, 2 * kreps[3] + 1),
                           0:kreps[1], -kreps[2]:kreps[2], -kreps[3]:kreps[3])
    # take advantage of symmetry. cos(k ⋅ dx) = cos( (-k) ⋅ dx)
    #   don't include both [ka kb kc] [-ka -kb -kc] for all kb, kc
    #   hence ka goes from 0:k_repfactors[3]
    for ka = 0:kreps[1], kb = -kreps[2]:kreps[2], kc=-kreps[3]:kreps[3]
        # don't include home unit cell
        if (ka == 0) && (kb == 0) && (kc == 0)
            kvec_wts[ka, kb, kc] = NaN
            continue
        end
        # if ka == 0, don't include both [0 1 x] and [0 -1 -x]
        #  but need [0 1 x] [0 1 -x]
        if (ka == 0) && (kb < 0)
            kvec_wts[ka, kb, kc] = NaN
            continue
        end
        # don't include both [0 0 1] and [0 0 -1]
        if (ka == 0) && (kb == 0) && (kc < 0)
            kvec_wts[ka, kb, kc] = NaN
            continue
        end

        # reciprocal vector, k
        k = sim_box.reciprocal_lattice * [ka, kb, kc]
        # |k|²
        norm_squared = dot(k, k)
        
        # factor of 2 from cos(-k⋅(x-xᵢ)) + cos(k⋅(x-xᵢ)) = 2 cos(k⋅(x-xᵢ))
        #  and how we include ka>=0 only and the two if statements above
        kvec_wts[ka, kb, kc] = 2 * exp(- norm_squared / (4.0 * α ^ 2)) / norm_squared / ϵ₀
    end
    return kvec_wts
end
        
function allocate_eikx(kreps::Tuple{Int, Int, Int})
    eika = OffsetArray(Complex{Float64}, 0:kreps[1])
    eikb = OffsetArray(Complex{Float64}, -kreps[2]:kreps[2])
    eikc = OffsetArray(Complex{Float64}, -kreps[3]:kreps[3])
    return eika, eikb, eikc
end

@unsafe function fill_eikx!(eikx::OffsetArray{Complex{Float64}},
                    kx_dot_dx::Float64,
                    krep::Int, include_neg_reps::Bool)
    # explicitly compute for k = 0 and k = 1
    eikx[0] = exp(0.0 * im)
    @fastmath eikx[1] = exp(im * kx_dot_dx)

    # recursion relation for higher frequencies
    for k = 2:krep
        eikx[k] = eikx[k - 1] * eikx[1]
    end

    # negative kreps are complex conjugate 
    if include_neg_reps
        for k = -krep:-1
            eikx[k] = conj(eikx[-k])
        end
    end
end

"""
    ϕ = electrostatic_potential(framework, x, simulation_box, repfactors, sr_cutoff_radius,
                                kvectors, α)

Compute the electrostatic potential at a point x inside of a framework. The electrostatic 
potential is created by the point charges assigned to the framework atoms. Periodic boundary 
conditions are applied through the Ewald summation.

# Arguments
- `framework::Framework`: Crystal structure (see `framework.charges` for charges)
- `x::Array{Float64, 1}`: point (Cartesian coordinates) at which we compute the electrostatic
    potential
- `sim_box::Box`: the simulation box for the Ewald summation. This must be consistent with 
    the repfactors applied to `framework.box`. We could compute it in here but we pass it 
    to the function instead for speed.
- `repfactors::Tuple{Int64, Int64, Int64}`: The replication factors used to form the super
    cell defined by `sim_box` from the home unit cell in `framework`.
- `sr_cutoff_radius::Float64`: the short-range cutoff radius for the real space sum in the
    Ewald summation.
- `kvectors::Array{Kvector, 1}`: pre-computed (for speed) reciprocal lattice vectors and
    weights associated with them for the Fourier space sum in the Ewald summation.
- `α::Float64`: Ewald sum convergence parameter. Units: inverse Å.
"""
function electrostatic_potential(framework::Framework, x::Array{Float64, 1}, 
                                 sim_box::Box, repfactors::Tuple{Int, Int, Int}, sr_cutoff_radius::Float64,
                                 eika::OffsetArray{Complex{Float64}},
                                 eikb::OffsetArray{Complex{Float64}},
                                 eikc::OffsetArray{Complex{Float64}},
                                 kreps::Tuple{Int, Int, Int}, 
                                 kvec_wts::OffsetArray{Float64,3,Array{Float64,3}}, α::Float64)
    # fractional coordinate of point at which we're computing the electrostatic potential (wrap to [0, 1.0])
    xf = mod.(framework.box.c_to_f * x, 1.0)

	sr_potential = 0.0 # short-range contribution
	lr_potential::Float64 = 0.0 # long-range contribution
    # build the super cell with loop over (ra, rb, rc) and atoms in the primitive unit cell
    for ra = 0:(repfactors[1] - 1), rb = 0:(repfactors[2] - 1), rc = 0:(repfactors[3] - 1)
        # x - x_i vector in fractional coords; i corresponds to atom of framework.
        #  same as: dxf = xf - (framework.xf[:, i] + [ra, rb, rc]) = (xf - [ra, rb, rc]) - framework.xf[:, i]
        @inbounds dxf = broadcast(-, xf - [ra, rb, rc], framework.xf)
        # convert distance vector to Cartesian coordinates
        @inbounds dx = framework.box.f_to_c * dxf
        # compute dot product with each of the three reciprocal lattice vectors
        #   k_dot_dx[i, j] = kᵢ ⋅ (x - xⱼ)
        @inbounds k_dot_dx = transpose(sim_box.reciprocal_lattice) * dx
        
        ###
        #  Long-range contribution
        ###
        for i = 1:framework.n_atoms
            fill_eikx!(eika, k_dot_dx[1, i], kreps[1], false)
            fill_eikx!(eikb, k_dot_dx[2, i], kreps[2], true)
            fill_eikx!(eikc, k_dot_dx[3, i], kreps[3], true)

            for kc=-kreps[3]:kreps[3], kb = -kreps[2]:kreps[2], ka = 0:kreps[1]
                if ka == 0
                    if kb < 0
                        continue
                    end
                    if kb == 0 && kc <= 0
                        continue
                    end
                end
                @inbounds lr_potential += framework.charges[i] * kvec_wts[ka, kb, kc] * real(eika[ka] * eikb[kb] * eikc[kc])
            end
        end

        ###
        #  Short range contribution
        ###
        # apply nearest image convention for periodic BCs
        nearest_image!(dxf, repfactors)

        # convert distance vector to Cartesian coordinates
        @inbounds dx = framework.box.f_to_c * dxf

        for i = 1:framework.n_atoms
            @inbounds @fastmath r = sqrt(dx[1, i] * dx[1, i] + dx[2, i] * dx[2, i] + dx[3, i] * dx[3, i])
            if r < sr_cutoff_radius
                @inbounds @fastmath sr_potential += framework.charges[i] / r* erfc(r * α)
            end
        end
    end
    sr_potential /= 4.0 * π * ϵ₀
    lr_potential /= sim_box.Ω

    return (lr_potential + sr_potential)::Float64
end		

function electrostatic_potential_energy(framework::Framework, molecule::Molecule, 
                                        sim_box::Box, repfactors::Tuple{Int, Int, Int}, sr_cutoff_radius::Float64,
                                        kvectors::Array{Kvector, 1}, α::Float64)
    ϕ = 0.0
    for charge in molecule.charges
        ϕ += charge.q * electrostatic_potential(framework, charge.x, sim_box, repfactors, sr_cutoff_radius, kvectors, α)
    end
    return ϕ
end

"""
    ϕ = electrostatic_potential(molecules, exclude_molecule_id, x, sim_box,
                                sr_cutoff_radius, kvectors, α)

Compute the electrostatic potential generated by a set of molecules (in `molecules`) in a 
simulation box, excluding the contribution of `molecules[exclude_molecule_id]`.
"""
 # #TODO TEST THIS!
 # function electrostatic_potential(molecules::Array{Molecule, 1}, exclude_molecule_id::Int,
 #                                  x::Array{Float64, 1}, sim_box::Box,
 #                                  sr_cutoff_radius::Float64,
 #                                  kvectors::Array{Kvector, 1}, α::Float64)
 # 	sr_potential = 0.0 # short-range contribution
 # 	lr_potential = 0.0 # long-range contribution
 #     for (i, molecule) in enumerate(molecules)
 #         if i == exclude_molecule_id
 #             continue
 #         end
 #         
 #         for charge in molecule.charges
 #             # vector from pt charge to pt of interest x in Cartesian coordinates
 #             dx = x - charge.x
 #             
 #             ###
 #             #  Long-range contribution
 #             ###
 #             @simd for kvector in kvectors
 #                 # k ⋅ (x - x_a)
 #                 @inbounds k_dot_dx = kvector.k[1] * dx[1] + kvector.k[2] * dx[2] + kvector.k[3] * dx[3]
 #                 @inbounds @fastmath lr_potential += charge.q * cos(k_dot_dx) * kvector.weight
 #             end
 # 
 #             ###
 #             #  Short range contribution
 #             ###
 #             dxf = sim_box.c_to_f * dx
 #             # apply nearest image convention for periodic BCs
 #             nearest_image!(dxf, (1, 1, 1))
 # 
 #             # convert distance vector to Cartesian coordinates
 #             dx = sim_box.f_to_c * dxf
 # 
 #             @inbounds @fastmath r = sqrt(dx[1] * dx[1] + dx[2] * dx[2] + dx[3] * dx[3])
 #             if r < sr_cutoff_radius
 #                 @inbounds @fastmath sr_potential += charge.q / r* erfc(r * α)
 #             end
 #         end
 #     end
 #     sr_potential /= 4.0 * π * ϵ₀
 #     lr_potential /= sim_box.Ω
 # 
 #     return (lr_potential + sr_potential)::Float64
 # end
 # 
 # function electrostatic_potential_energy(molecules::Array{Molecule, 1}, molecule_id::Int, 
 #                                         sim_box::Box, sr_cutoff_radius::Float64,
 #                                         kvectors::Array{Kvector, 1}, α::Float64)
 #     ϕ = 0.0
 #     for charge in molecules[molecule_id].charges
 #         ϕ += charge.q * electrostatic_potential(molecules, molecule_id, charge.x, sim_box,
 #                                                 sr_cutoff_radius, kvectors, α)
 #     end
 #     return ϕ
 # end
