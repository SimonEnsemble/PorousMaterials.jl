import Base: +, -, /, *

"""
Data structure to store potential energy, partitioned into van der Waals (`energy.vdw`) 
and electrostatic (`energy.coulomb`) interactions, both `Float64`.
"""
mutable struct PotentialEnergy
    vdw::Float64 # contribution from van der Waals interactions
    coulomb::Float64 # contribution from electrostatic interactions
end
PotentialEnergy() = PotentialEnergy(0.0, 0.0) # constructor
Base.sum(energy::PotentialEnergy) = energy.vdw + energy.coulomb
+(u::PotentialEnergy, v::PotentialEnergy) = PotentialEnergy(u.vdw + v.vdw, u.coulomb + v.coulomb)
-(u::PotentialEnergy, v::PotentialEnergy) = PotentialEnergy(u.vdw - v.vdw, u.coulomb - v.coulomb)
*(a::Float64, energy::PotentialEnergy) = PotentialEnergy(a * energy.vdw, a * energy.coulomb)
*(energy::PotentialEnergy, a::Float64) = *(a, energy)
/(energy::PotentialEnergy, a::Float64) = PotentialEnergy(energy.vdw / a, energy.coulomb / a)
square(u::PotentialEnergy) = PotentialEnergy(u.vdw ^ 2, u.coulomb ^ 2)
Base.sqrt(u::PotentialEnergy) = PotentialEnergy(sqrt(u.vdw), sqrt(u.coulomb))

function Base.isapprox(u::PotentialEnergy, v::PotentialEnergy; verbose::Bool=true, atol::Float64=1e-6)
    if ! isapprox(u.vdw, v.vdw, atol=atol)
        if verbose
            warn("vdw energy mismatch")
        end
        return false
    end
    if ! isapprox(u.coulomb, v.coulomb, atol=atol)
        if verbose
            warn("coulomb energy mismatch")
        end
        return false
    end
    return true
end

# data structures to facilitate storing/partitioning potential energy of a system
mutable struct SystemPotentialEnergy
    guest_host::PotentialEnergy
    guest_guest::PotentialEnergy
end
SystemPotentialEnergy() = SystemPotentialEnergy(PotentialEnergy(), PotentialEnergy()) # constructor
Base.sum(v::SystemPotentialEnergy) = v.guest_guest.vdw + v.guest_guest.coulomb + 
                                     v.guest_host.vdw  + v.guest_host.coulomb
+(u::SystemPotentialEnergy, v::SystemPotentialEnergy) = SystemPotentialEnergy(u.guest_host  + v.guest_host, 
                                                                              u.guest_guest + v.guest_guest)
-(u::SystemPotentialEnergy, v::SystemPotentialEnergy) = SystemPotentialEnergy(u.guest_host  - v.guest_host, 
                                                                              u.guest_guest - v.guest_guest)
*(u::SystemPotentialEnergy, a::Float64) = SystemPotentialEnergy(a * u.guest_host, a * u.guest_guest)
*(a::Float64, u::SystemPotentialEnergy) = *(u::SystemPotentialEnergy, a::Float64)
/(u::SystemPotentialEnergy, a::Float64) = SystemPotentialEnergy(u.guest_host / a, u.guest_guest / a)
Base.sqrt(u::SystemPotentialEnergy) = SystemPotentialEnergy(sqrt(u.guest_host), sqrt(u.guest_guest))
square(u::SystemPotentialEnergy) = SystemPotentialEnergy(square(u.guest_host), square(u.guest_guest))

function Base.isapprox(u::SystemPotentialEnergy, v::SystemPotentialEnergy; 
                       verbose::Bool=true, atol::Float64=1e-6)
    if ! isapprox(u.guest_host, v.guest_host, verbose=verbose, atol=atol)
        warn("(guest-host mismatch)")
        return false
    end
    if ! isapprox(u.guest_guest, v.guest_guest, verbose=verbose, atol=atol)
        warn("(guest-guest mismatch)")
        return false
    end
    return true
end
