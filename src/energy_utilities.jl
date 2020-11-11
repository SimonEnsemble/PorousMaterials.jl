import Base: +, -, /, *

mutable struct PotentialEnergy
    vdw::Float64     # contribution from van der Waals interactions
    es::Float64      # contribution from electrostatic interactions
end

PotentialEnergy() = PotentialEnergy(0.0, 0.0) # constructor

Base.sum(energy::PotentialEnergy) = energy.vdw + energy.es
+(u::PotentialEnergy, v::PotentialEnergy) = PotentialEnergy(u.vdw + v.vdw, u.es + v.es)
-(u::PotentialEnergy, v::PotentialEnergy) = PotentialEnergy(u.vdw - v.vdw, u.es - v.es)
*(a::Float64, energy::PotentialEnergy) = PotentialEnergy(a * energy.vdw, a * energy.es)
*(energy::PotentialEnergy, a::Float64) = *(a, energy)
/(energy::PotentialEnergy, a::Float64) = PotentialEnergy(energy.vdw / a, energy.es / a)
square(u::PotentialEnergy) = PotentialEnergy(u.vdw ^ 2, u.es ^ 2)
Base.sqrt(u::PotentialEnergy) = PotentialEnergy(sqrt(u.vdw), sqrt(u.es))

function Base.isapprox(u::PotentialEnergy, v::PotentialEnergy; verbose::Bool=true, atol::Float64=1e-6)
    if ! isapprox(u.vdw, v.vdw, atol=atol)
        if verbose
            @warn "vdw energy mismatch"
        end
        return false
    end
    if ! isapprox(u.es, v.es, atol=atol)
        if verbose
            @warn "es energy mismatch"
        end
        return false
    end
    return true
end

mutable struct SystemPotentialEnergy
    gh::PotentialEnergy
    gg::PotentialEnergy
end
SystemPotentialEnergy() = SystemPotentialEnergy(PotentialEnergy(), PotentialEnergy()) # constructor
Base.sum(v::SystemPotentialEnergy) = v.gg.vdw + v.gg.es +
                                     v.gh.vdw  + v.gh.es
+(u::SystemPotentialEnergy, v::SystemPotentialEnergy) = SystemPotentialEnergy(u.gh  + v.gh,
                                                                              u.gg + v.gg)
-(u::SystemPotentialEnergy, v::SystemPotentialEnergy) = SystemPotentialEnergy(u.gh  - v.gh,
                                                                              u.gg - v.gg)
*(u::SystemPotentialEnergy, a::Float64) = SystemPotentialEnergy(a * u.gh, a * u.gg)
*(a::Float64, u::SystemPotentialEnergy) = *(u::SystemPotentialEnergy, a::Float64)
/(u::SystemPotentialEnergy, a::Float64) = SystemPotentialEnergy(u.gh / a, u.gg / a)
Base.sqrt(u::SystemPotentialEnergy) = SystemPotentialEnergy(sqrt(u.gh), sqrt(u.gg))
square(u::SystemPotentialEnergy) = SystemPotentialEnergy(square(u.gh), square(u.gg))

function Base.isapprox(u::SystemPotentialEnergy, v::SystemPotentialEnergy;
                       verbose::Bool=true, atol::Float64=1e-6)
    if ! isapprox(u.gh, v.gh, verbose=verbose, atol=atol)
        @warn "(guest-host mismatch)"
        return false
    end
    if ! isapprox(u.gg, v.gg, verbose=verbose, atol=atol)
        @warn "(guest-guest mismatch)"
        return false
    end
    return true
end
