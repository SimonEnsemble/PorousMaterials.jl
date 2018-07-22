"""
    nearest_image!(dxf)

Applies the nearest image convention on a vector `dxf` between two atoms in fractional 
space; modifies `dxf` for nearest image convention. Fractional coordinates here fall in
[0, 1] so that the box is [0, 1]^3 in fractional space.

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

@inline function nearest_r²(xf::Array{Float64, 1}, yf::Array{Float64, 1}, box::Box)
    # vector from y to x in fractional coordinate space
    @inbounds dxf = xf - yf
    # apply nearest image convention for periodic boundary conditions
    nearest_image!(dxf)
    # convert to cartesian
    @inbounds dxf .= box.f_to_c * dxf
    # return r²
    @inbounds return dxf[1] * dxf[1] + dxf[2] * dxf[2] + dxf[3] * dxf[3]
end

@inline function nearest_r(xf::Array{Float64, 1}, yf::Array{Float64, 1}, box::Box)
    return sqrt(nearest_r²(xf, yf, box))
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
#        |--dxf--|
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
