using SpecialFunctions
using OffsetArrays
# Vacuum permittivity eps0 = 8.854187817e-12 # C^2/(J-m)
# 1 m = 1e10 A, 1 e = 1.602176e-19 C, kb = 1.3806488e-23 J/K
# 8.854187817e-12 C^2/(J-m) [1 m / 1e10 A] [1 e / 1.602176e-19 C]^2 [kb = 1.3806488e-23 J/K]
const ϵ₀ = 4.7622424954949676e-7  # \epsilon_0 vacuum permittivity units: electron charge^2 /(A - K)

# TODO use NIST test data https://www.nist.gov/mml/csd/chemical-informatics-research-group/spce-water-reference-calculations-10%C3%A5-cutoff
#   to untar: tar zxvf spce_sample_configurations-tar.gz

struct EwaldParams
    kreps::Tuple{Int, Int, Int}
    α::Float64
    sr_cutoff_r::Float64
    box::Box
end

"""
    kvectors = precompute_kvec_wts(sim_box, kreps, α)

For speed, pre-compute the weights for each reciprocal lattice vector for the Ewald sum in
Fourier space. This function takes advantage of the symmetry:
    cos(-k⋅(x-xᵢ)) + cos(k⋅(x-xᵢ)) = 2 cos(k⋅(x-xᵢ))
Returns an OffsetArray `kvec_wts` so e.g. `kvec_wts[0, -5, 4]` corresponds to the weight 
with ka=0, kb=-5, kc=4.

# Arguments
- `sim_box::Box`: the simulation box that consits of replicated primitive unit cells of the framework
- `kreps::Tuple{Int, Int, Int}`: The number of reciprocal lattice space replications to include in the long-range Ewald sum.
- `α::Float64`: Ewald sum convergence parameter. Units: inverse Å.
"""
function precompute_kvec_wts(eparams::EwaldParams)
    kvec_wts = OffsetArray(Float64, 0:eparams.kreps[1], -eparams.kreps[2]:eparams.kreps[2], -eparams.kreps[3]:eparams.kreps[3])
    # take advantage of symmetry. cos(k ⋅ dx) = cos( (-k) ⋅ dx)
    #   don't include both [ka kb kc] [-ka -kb -kc] for all kb, kc
    #   hence ka goes from 0:k_repfactors[3]
    for ka = 0:eparams.kreps[1], kb = -eparams.kreps[2]:eparams.kreps[2], kc=-eparams.kreps[3]:eparams.kreps[3]
        # don't include the home unit cell
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
        k = eparams.box.reciprocal_lattice * [ka, kb, kc]
        # |k|²
        norm_squared = dot(k, k)
        
        # factor of 2 from cos(-k⋅(x-xᵢ)) + cos(k⋅(x-xᵢ)) = 2 cos(k⋅(x-xᵢ))
        #  and how we include ka>=0 only and the two if statements above
        kvec_wts[ka, kb, kc] = 2 * exp(- norm_squared / (4.0 * eparams.α ^ 2)) / norm_squared / ϵ₀
    end
    return kvec_wts
end

"""
    eikr = allocate_eikr(kreps)

Allocate OffsetArrays for storing e^{i * k ⋅ r} where r = x - xⱼ and k is a reciprocal
lattice vector.
eikra has indices 0:kreps[1] and corresponds to recip. vector in a-direction
eikrb has indices -kreps[2]:kreps[2] and corresponds to recip. vector in b-direction
eikrc has indices -kreps[3]:kreps[3] and corresponds to recip. vector in c-direction
"""
function setup_Ewald_sum(kreps::Tuple{Int, Int, Int}, α::Float64, sr_cutoff_r::Float64, sim_box::Box)
    eparams = EwaldParams(kreps, α, sr_cutoff_r, sim_box)
    kvec_wts = precompute_kvec_wts(eparams)
    eikar = OffsetArray(Complex{Float64}, 0:kreps[1])
    eikbr = OffsetArray(Complex{Float64}, -kreps[2]:kreps[2])
    eikcr = OffsetArray(Complex{Float64}, -kreps[3]:kreps[3])
    return eparams, kvec_wts, eikar, eikbr, eikcr
end

"""
    fill_eikr!(eikr, k_dot_dx, krep, include_neg_reps)

Given k ⋅ r, where r = x - xⱼ, compute e^{i n k ⋅ r} for n = 0:krep to fill the OffsetArray
`eikr`. If `include_neg_reps` is true, also compute n=-krep:-1.
"""
@unsafe function fill_eikr!(eikr::OffsetArray{Complex{Float64}}, k_dot_r::Float64, 
                            krep::Int, include_neg_reps::Bool)
    # explicitly compute for k = 0 and k = 1
    eikr[0] = exp(0.0 * im)
    @fastmath eikr[1] = exp(im * k_dot_r)

    # recursion relation for higher frequencies to avoid expensive computing of cosine.
    #  e^{3 * i * k_dot_r} = e^{2 * i * k_dot_r} * e^{ i * k_dot_r}
    for k = 2:krep
        eikr[k] = eikr[k - 1] * eikr[1]
    end

    # negative kreps are complex conjugate of positive ones.
    #  e^{2 * i * k_dot_r} = conj(e^{-2 * i * k_dot_dr})
    if include_neg_reps
        for k = -krep:-1
            eikr[k] = conj(eikr[-k])
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
function electrostatic_potential(framework::Framework,
                                 x::Array{Float64, 1},
                                 repfactors::Tuple{Int, Int, Int},
                                 eparams::EwaldParams,
                                 kvec_wts::OffsetArray{Float64,3,Array{Float64,3}},
                                 eikar::OffsetArray{Complex{Float64}},
                                 eikbr::OffsetArray{Complex{Float64}},
                                 eikcr::OffsetArray{Complex{Float64}})
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
        # TODO make recip lattice transposed...
        @inbounds k_dot_dx = transpose(eparams.box.reciprocal_lattice) * dx
        
        ###
        #  Long-range contribution
        ###
        for i = 1:framework.n_atoms
            # compute e^{ i * n * k * (x - xᵢ)} for n = -krep:krep
            fill_eikr!(eikar, k_dot_dx[1, i], eparams.kreps[1], false) # via symmetry only need +ve
            fill_eikr!(eikbr, k_dot_dx[2, i], eparams.kreps[2], true)
            fill_eikr!(eikcr, k_dot_dx[3, i], eparams.kreps[3], true)
            
            # loop over kevectors
            for kc=-eparams.kreps[3]:eparams.kreps[3], kb = -eparams.kreps[2]:eparams.kreps[2], ka = 0:eparams.kreps[1]
                # same logic as in `precompute_kvec_wts` but faster.
                if ka == 0
                    if kb < 0
                        continue
                    end
                    if kb == 0 && kc <= 0
                        continue
                    end
                end
                # cos( i * this_k * r) = real part of:
                #     e^{i ka r * vec(ka)} * 
                #     e^{i kb r * vec(kb)} * 
                #     e^{i kb r * vec(kc)} where r = x - xᵢ
                #   and eikar[ka], eikbr[kb], eikcr[kc] contain exactly the above components.
                @unsafe @inbounds lr_potential += framework.charges[i] * kvec_wts[ka, kb, kc] * real(eikar[ka] * eikbr[kb] * eikcr[kc])
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
            if r < eparams.sr_cutoff_r
                @inbounds @fastmath sr_potential += framework.charges[i] / r * erfc(r * eparams.α)
            end
        end
    end
    sr_potential /= 4.0 * π * ϵ₀
    lr_potential /= eparams.box.Ω

    return (lr_potential + sr_potential)::Float64
end		

function electrostatic_potential_energy(framework::Framework,
                                     molecule::Molecule,
                                     repfactors::Tuple{Int, Int, Int},
                                     eparams::EwaldParams,
                                     kvec_wts::OffsetArray{Float64,3,Array{Float64,3}},
                                     eikar::OffsetArray{Complex{Float64}},
                                     eikbr::OffsetArray{Complex{Float64}},
                                     eikcr::OffsetArray{Complex{Float64}})
    ϕ = 0.0
    for charge in molecule.charges
        ϕ += charge.q * electrostatic_potential(framework, charge.x, repfactors, eparams, kvec_wts, eikar, eikbr, eikcr)
    end
    return ϕ
end

"""
    ϕ = electrostatic_potential(molecules, exclude_molecule_id, x, sim_box,
                                sr_cutoff_radius, kvectors, α)

Compute the electrostatic potential generated by a set of molecules (in `molecules`) in a 
simulation box, excluding the contribution of `molecules[exclude_molecule_id]`.
"""
#TODO TEST THIS!
function electrostatic_potential(molecules::Array{Molecule, 1},
                                 exclude_molecule_id::Int,
                                 x::Array{Float64, 1},
                                 eparams::EwaldParams,
                                 kvec_wts::OffsetArray{Float64,3,Array{Float64,3}},
                                 eikar::OffsetArray{Complex{Float64}},
                                 eikbr::OffsetArray{Complex{Float64}},
                                 eikcr::OffsetArray{Complex{Float64}})
	sr_potential = 0.0 # short-range contribution
	lr_potential::Float64 = 0.0 # long-range contribution
    for (i, molecule) in enumerate(molecules)
        if i == exclude_molecule_id
            continue
        end
        
        for charge in molecule.charges
            # vector from pt charge to pt of interest x in Cartesian coordinates
            dx = x - charge.x
            # reciprocal lattice vectors ⋅ dx
            k_dot_dx = transpose(eparams.box.reciprocal_lattice) * dx
        
            ###
            #  Long-range contribution
            ###
            # compute e^{ i * n * k * (x - charge.x)} for n = -krep:krep
            fill_eikr!(eikar, k_dot_dx[1], eparams.kreps[1], false) # via symmetry only need +ve
            fill_eikr!(eikbr, k_dot_dx[2], eparams.kreps[2], true)
            fill_eikr!(eikcr, k_dot_dx[3], eparams.kreps[3], true)
            
            # loop over kevectors
            for kc=-eparams.kreps[3]:eparams.kreps[3], kb = -eparams.kreps[2]:eparams.kreps[2], ka = 0:eparams.kreps[1]
                # same logic as in `precompute_kvec_wts` but faster.
                if ka == 0
                    if kb < 0
                        continue
                    end
                    if kb == 0 && kc <= 0
                        continue
                    end
                end
                # cos( i * this_k * r) = real part of:
                #     e^{i ka r * vec(ka)} * 
                #     e^{i kb r * vec(kb)} * 
                #     e^{i kb r * vec(kc)} where r = x - xᵢ
                #   and eikar[ka], eikbr[kb], eikcr[kc] contain exactly the above components.
                @unsafe @inbounds lr_potential += charge.q * kvec_wts[ka, kb, kc] * real(eikar[ka] * eikbr[kb] * eikcr[kc])
            end
            
            ###
            #  Short range contribution
            ###
            dxf = eparams.box.c_to_f * dx
            # apply nearest image convention for periodic BCs
            nearest_image!(dxf, (1, 1, 1))

            # convert distance vector to Cartesian coordinates
            dx = eparams.box.f_to_c * dxf

            @inbounds @fastmath r = sqrt(dx[1] * dx[1] + dx[2] * dx[2] + dx[3] * dx[3])
            if r < eparams.sr_cutoff_r
                @inbounds @fastmath sr_potential += charge.q / r * erfc(r * eparams.α)
            end
        end
    end
    sr_potential /= 4.0 * π * ϵ₀
    lr_potential /= eparams.box.Ω

    return (lr_potential + sr_potential)::Float64
end

function electrostatic_potential_energy(molecules::Array{Molecule, 1},
                                        molecule_id::Int,
                                     eparams::EwaldParams,
                                     kvec_wts::OffsetArray{Float64,3,Array{Float64,3}},
                                     eikar::OffsetArray{Complex{Float64}},
                                     eikbr::OffsetArray{Complex{Float64}},
                                     eikcr::OffsetArray{Complex{Float64}})
    ϕ = 0.0
    for charge in molecules[molecule_id].charges
        ϕ += charge.q * electrostatic_potential(molecules, molecule_id, charge.x, eparams,
                                                kvec_wts, eikar, eikbr, eikcr)
    end
    return ϕ
end
