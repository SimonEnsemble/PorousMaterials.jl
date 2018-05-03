using SpecialFunctions # for erfc
using OffsetArrays
using Roots # for fzero
# Vacuum permittivity eps0 = 8.854187817e-12 # C^2/(J-m)
# 1 m = 1e10 A, 1 e = 1.602176e-19 C, kb = 1.3806488e-23 J/K
# 8.854187817e-12 C^2/(J-m) [1 m / 1e10 A] [1 e / 1.602176e-19 C]^2 [kb = 1.3806488e-23 J/K]
const ϵ₀ = 4.7622424954949676e-7  # \epsilon_0 vacuum permittivity units: electron charge^2 /(A - K)

# TODO use NIST test data https://www.nist.gov/mml/csd/chemical-informatics-research-group/spce-water-reference-calculations-10%C3%A5-cutoff
#   to untar: tar zxvf spce_sample_configurations-tar.gz

"Data structure for storing Ewald summation settings"
struct EwaldParams
    "number of replications in k-space in a, b, c directions"
    kreps::Tuple{Int, Int, Int}
    "Ewald summation convergence parameter (units: inverse Å)"
    α::Float64
    "short-range cutoff radius (units: Å)"
    sr_cutoff_r::Float64
    "simulation box in which Ewald sum is conducted (contains reciprocal lattice)"
    box::Box
end

function Base.print(io::IO, eparams::EwaldParams)
    @printf(io, "Ewald summation parameters:\n")
    @printf(io, "\tk-replication factors: %d %d %d\n", eparams.kreps[1], eparams.kreps[2], eparams.kreps[3])
    @printf(io, "\tEwald convergence param. α = %f\n", eparams.α)
    @printf(io, "\tshort-range cutoff radius (Å): %f\n", eparams.sr_cutoff_r)
end

"""
Data structure for a k-vector replication, although the reciprocal lattice is needed to
specify the exact k-vector.
Stores the weight determining the contribution of kvector (ka, kb, kc) to the Fourier
component of the Ewald sum. Ewald sums work with an array of `Kvec`'s.
"""
struct Kvector
    ka::Int
    kb::Int
    kc::Int

    wt::Float64
end

"""
    α, max_mag_k_sqrd = determine_ewald_params(sr_cutoff_r, ϵ)

For a given precision `ϵ` and short-range cutoff radius, compute the Ewald summation 
convergence parameter α and the cutoff value for |k|² in the reciprocal space.
"""
function determine_ewald_params(sr_cutoff_r::Float64, ϵ::Float64)
    # α is Ewald sum convergence parameter. It determines how fast the long- and short-
    #  range interactions decay with r in the Fourier and real space sum, respectively.
    #    erfc(α * r) / r dependence for short-range
    #    exp(-mag_k_squared / (4α²)) / mag_k_squared dependence for long-range
    #   as α increases, short-range contributions decays with r more quickly;
    #   thus need smaller sr_cutoff. but long-range decays more slowly with r with larger
    #   α. thus would need more kevectors.
    #      this approach is similar to DLPoly according to here:
    #       http://www.cse.scitech.ac.uk/ccg/software/DL_POLY_CLASSIC/FAQ/FAQ2.shtml
    # solve for α to allow short-range interactions to decay beyond the cutoff radius. 
    sr_error(α) = erfc(sr_cutoff_r * α) / sr_cutoff_r - ϵ
    α = fzero(sr_error, 0.0, 25.0) # 25.0 should be safe bracket.
    # solve for the max squre norm of the k-vectors required to allow long-range 
    #   interactions in Fourier space to decay beyond for k-vectors longer than this. 
    lr_error(mag_k_sqrd) = 2 * exp(- mag_k_sqrd / (4.0 * α ^ 2)) / mag_k_sqrd - ϵ
    max_mag_k_sqrd = fzero(lr_error, 0.0, 100.0)
    return α, max_mag_k_sqrd
end

"""
Determine replications of k-vectors in Fourier sum, `kreps`, required to assert all 
k-vectors with square magnitude less than |k|² (`max_mag_k_sqrd`) are included.
"""
function required_kreps(box::Box, max_mag_k_sqrd::Float64)
    # fill in later
    kreps = [0, 0, 0]
    for abc = 1:3
        n = [0, 0, 0]
        while true
            n[abc] += 1
            # compute reciprocal vector, k
            k = box.reciprocal_lattice * n
            if dot(k, k) > max_mag_k_sqrd
                kreps[abc] = n[abc] - 1
                break
            end
        end
    end
    return (kreps[1], kreps[2], kreps[3])
end

"""
    kvectors = precompute_kvec_wts(eparams, max_mag_k_sqrd=Inf)

For speed, pre-compute the weights for each reciprocal lattice vector for the Ewald sum in
Fourier space. This function takes advantage of the symmetry:
    cos(-k⋅(x-xᵢ)) + cos(k⋅(x-xᵢ)) = 2 cos(k⋅(x-xᵢ))
Returns an array of `Kvector`'s.

If `max_mag_k_sqrd` is passed, k-vectors with a magnitude greater than `max_mag_k_sqrd` are
not included.
"""
function precompute_kvec_wts(eparams::EwaldParams, max_mag_k_sqrd::Float64=Inf)
    kvecs = Kvector[]
    # take advantage of symmetry. cos(k ⋅ dx) = cos( (-k) ⋅ dx)
    #   don't include both [ka kb kc] [-ka -kb -kc] for all kb, kc
    #   hence ka goes from 0:k_repfactors[1]
    for ka = 0:eparams.kreps[1], kb = -eparams.kreps[2]:eparams.kreps[2], kc=-eparams.kreps[3]:eparams.kreps[3]
        # don't include the home unit cell
        if (ka == 0) && (kb == 0) && (kc == 0)
            continue
        end
        # if ka == 0, don't include both [0 1 x] and [0 -1 -x]
        #  but need [0 1 x] [0 1 -x]
        if (ka == 0) && (kb < 0)
            continue
        end
        # don't include both [0 0 1] and [0 0 -1]
        if (ka == 0) && (kb == 0) && (kc < 0)
            continue
        end

        # compute reciprocal vector, k
        k = eparams.box.reciprocal_lattice * [ka, kb, kc]
        # |k|²
        mag_k_sqrd = dot(k, k)
        
        # factor of 2 from cos(-k⋅(x-xᵢ)) + cos(k⋅(x-xᵢ)) = 2 cos(k⋅(x-xᵢ))
        #  and how we include ka>=0 only and the two if statements above
        if mag_k_sqrd < max_mag_k_sqrd
            wt = 2 * exp(- mag_k_sqrd / (4.0 * eparams.α ^ 2)) / mag_k_sqrd / ϵ₀
            push!(kvecs, Kvector(ka, kb, kc, wt))
        end
    end
    return kvecs
end

"""
    eparams, kvecs, eikar, eikbr, eikcr = setup_Ewald_sum(sr_cutoff_r, sim_box, ϵ=1e-6, verbose=false)

Given the short-range cutoff radius `sr_cutoff_r`, automatically compute Ewald convergence 
parameter `α` and number of k-vector replications in Fourier space, `kreps`.
Construct Ewald parameters data type.
Pre-compute weights on k-vector contributions to Ewald sum in Fourier space, return as 
OffsetArray.
Allocate OffsetArrays for storing e^{i * k ⋅ r} where r = x - xⱼ and k is a reciprocal
lattice vector.
eikra has indices 0:kreps[1] and corresponds to recip. vector in a-direction
eikrb has indices -kreps[2]:kreps[2] and corresponds to recip. vector in b-direction
eikrc has indices -kreps[3]:kreps[3] and corresponds to recip. vector in c-direction
"""
function setup_Ewald_sum(sr_cutoff_r::Float64, sim_box::Box; ϵ::Float64=1e-6, verbose::Bool=false)
    # determine Ewald convergence parameter and cutoff for |k|² for sum in Fourier space
    α, max_mag_k_sqrd = determine_ewald_params(sr_cutoff_r, ϵ)
    # determine kreps required to assert all k-vectors with magnitude less than |k|² are included.
    kreps = required_kreps(sim_box, max_mag_k_sqrd)
    # construct Ewald parameters
    eparams = EwaldParams(kreps, α, sr_cutoff_r, sim_box)
    if verbose
        print(eparams)
    end
    # pre-compute weights associated with each k-vector in sum in Fourier space.
    #  discards k-vector replications whose |k|² > max_mag_k_sqrd for efficiency.
    kvectors = precompute_kvec_wts(eparams, max_mag_k_sqrd)
    # pre-allocate memory for e^{i k vec(k) ⋅ r}. Don't put these in EwaldParams for speed,
    #  so they are passed as reference. These are OffsetArrays, which changes indexing scheme.
    eikar = OffsetArray(Complex{Float64}, 0:kreps[1])
    eikbr = OffsetArray(Complex{Float64}, -kreps[2]:kreps[2])
    eikcr = OffsetArray(Complex{Float64}, -kreps[3]:kreps[3])
    return eparams, kvectors, eikar, eikbr, eikcr
end

"""
    fill_eikr!(eikr, k_dot_dx, krep, include_neg_reps)
Given k ⋅ r:w
, where r = x - xⱼ, compute e^{i n k ⋅ r} for n = 0:krep to fill the OffsetArray
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
                                 kvectors::Array{Kvector, 1},
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
        # TODO should this be over primitive unit cell or simulation box?
        for i = 1:framework.n_atoms
            # compute e^{ i * n * k * (x - xᵢ)} for n = -krep:krep
            fill_eikr!(eikar, k_dot_dx[1, i], eparams.kreps[1], false) # via symmetry only need +ve
            fill_eikr!(eikbr, k_dot_dx[2, i], eparams.kreps[2], true)
            fill_eikr!(eikcr, k_dot_dx[3, i], eparams.kreps[3], true)
            
            # loop over kevectors
            for kvector in kvectors
                # cos( i * this_k * r) = real part of:
                #     e^{i ka r * vec(ka)} * 
                #     e^{i kb r * vec(kb)} * 
                #     e^{i kb r * vec(kc)} where r = x - xᵢ
                #   and eikar[ka], eikbr[kb], eikcr[kc] contain exactly the above components.
                @unsafe @inbounds lr_potential += framework.charges[i] * kvector.wt * real(
                        eikar[kvector.ka] * eikbr[kvector.kb] * eikcr[kvector.kc])
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
                                 kvectors::Array{Kvector, 1},
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
            for kvector in kvectors
                # cos( i * this_k * r) = real part of:
                #     e^{i ka r * vec(ka)} * 
                #     e^{i kb r * vec(kb)} * 
                #     e^{i kb r * vec(kc)} where r = x - xᵢ
                #   and eikar[ka], eikbr[kb], eikcr[kc] contain exactly the above components.
                @unsafe @inbounds lr_potential += charge.q * kvector.wt * real(
                        eikar[kvector.ka] * eikbr[kvector.kb] * eikcr[kvector.kc])
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
