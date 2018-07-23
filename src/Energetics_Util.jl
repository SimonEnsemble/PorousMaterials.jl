import Base: +, -, /, *

"""
    nearest_image!(dxf)

Applies the nearest image convention on a vector `dxf` between two atoms in fractional 
space; modifies `dxf` for nearest image convention. Fractional coordinates here fall in
[0, 1] so that the box is [0, 1]^3 in fractional space.

See comments in vdw_energy for more description.

# Arguments
- `dxf::Array{Float64, 1}`: A vector between two atoms in fractional coordinates
"""
function nearest_image!(dxf::Array{Float64, 1})
    for k = 1:3 # loop over components
        @inbounds if abs(dxf[k]) > 0.5
            @inbounds dxf[k] -= sign(dxf[k])
        end
    end
end

function nearest_image!(dxf::Array{Float64, 2})
    for a = 1:size(dxf)[2] # loop over atoms
        for k = 1:3 # loop over components
            @inbounds if abs(dxf[k, a]) > 0.5
                @inbounds dxf[k, a] -= sign(dxf[k, a])
            end
        end
    end
end

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

# Arni's notes on Nearest image convention.
#  If the interaction between the adsorbate molecule and atom k is being looked
#  at, we'll only look at the interaction between the adsorbate molecule and
#  the closest replication of atom k. This is done with fractional
#  coordinates for simplication and transformation to cartesian is done
#  later.


# NIC condensed into a for-loop
#
# If the absolute value of the distance between the adsorbate atom and the
# framework atom is greater than half the replication factor, we know that
# there is a closer replication of the framework atom.
#
# {Replicat.} ||{Supercell}||{Replicat.}
# |-----|----o||--x--|----o||-----|----o|
#				  |--dxf--|
#
# x = adsorbate atom, o = framework atom
#
# dxf is `x_adsorbate - x_framework` so when the adsorbate atom is to the left of
# the framework atom, dxf is negative.
# When correcting for the position of the framework atom with the Nearest Image Convention
# we use `sign(dxf[j]) * repfactors[j]` to change the distance dxf so it gives the distance
# between the adsorbate atom and the closest replication of the framework atom.
#
# In the above example, the framework atom is to the right of the adsorbate atom, and dxf < 0
# We see that the left replication of `o` is closer to `x`, and we should be calculating the
# distance between that atom and the adsorbate. So by subtracting `sign(dxf[j]) * repfactors[j]`
# (remember that `sign(dxf[j]) = -1` in this case) we're adding `repfactors[j]` to the dxf[j] value
# and dxf[j] becomes positive (which makes sense because we're calculating `x_ads - x_framework`)
# When checking the other case (where the adsorbate atom is to the right of the framework atom),
# we see that the same equation holds (because now `sign(dxf[j]) = 1`)
