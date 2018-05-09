import Base: +, -

"""
    nearest_image!(dxf, repfactors)

Applies the nearest image convention on a vector `dxf` between two atoms
in fractional space; modifies `dxf` for nearest image convention.

See comments in vdw_energy for more description.

# Arguments
- `dxf::Array{Float64, 1}`: A vector between two atoms in fractional coordinates
- `repfactors::Tuple{Int, Int, Int}`: Replication factors used to determine how many times the unit cell should be replicated
"""
function nearest_image!(dxf::Array{Float64, 1}, repfactors::Tuple{Int, Int, Int})
    for k = 1:3 # loop over components
        if abs(dxf[k]) > repfactors[k] / 2.0
            @inbounds dxf[k] -= sign(dxf[k]) * repfactors[k]
        end
    end
end

function nearest_image!(dxf::Array{Float64, 2}, repfactors::Tuple{Int, Int, Int})
    for a = 1:size(dxf)[2] # loop over atoms
        for k = 1:3 # loop over components
            if abs(dxf[k, a]) > repfactors[k] / 2.0
                @inbounds dxf[k, a] -= sign(dxf[k, a]) * repfactors[k]
            end
        end
    end
end

"""
Data structure containing the guest-host and guest-guest potential and electrostatic energy

# Attributes
- `vdw_gh::Float64`: Guest-host van der Waals energy
- `vdw_gg::Float64`: Guest-guest van der Waals energy
- `electro_gh::Float64`: Guest-host electrostatic energy
- `electro_gg::Float64`: Guest-guest electrostatic energy
"""
type PotentialEnergy
    vdw_gh::Float64
    vdw_gg::Float64
    electro_gh::Float64
    electro_gg::Float64
end

Base.sum(v::PotentialEnergy) = v.vdw_gh + v.vdw_gg + v.electro_gh + v.electro_gg
+(u::PotentialEnergy, v::PotentialEnergy) = PotentialEnergy(u.vdw_gh     + v.vdw_gh,
                                                            u.vdw_gg     + v.vdw_gg,
                                                            u.electro_gh + v.electro_gh,
                                                            u.electro_gg + v.electro_gg)
-(u::PotentialEnergy, v::PotentialEnergy) = PotentialEnergy(u.vdw_gh     - v.vdw_gh,
                                                            u.vdw_gg     - v.vdw_gg,
                                                            u.electro_gh - v.electro_gh,
                                                            u.electro_gg - v.electro_gg)
Base.isapprox(u::PotentialEnergy, v::PotentialEnergy) = (isapprox(u.vdw_gh    , v.vdw_gh) &&
                                                         isapprox(u.vdw_gg    , v.vdw_gg) &&
                                                         isapprox(u.electro_gh, v.electro_gh) &&
                                                         isapprox(u.electro_gg, v.electro_gg))

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
