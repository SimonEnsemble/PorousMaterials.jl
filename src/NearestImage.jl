"""
    nearest_image!(dxf)

Applies the nearest image convention on a vector `dxf` between two atoms in fractional
space; modifies `dxf` for nearest image convention. Fractional coordinates here fall in
[0, 1] so that the box is [0, 1]^3 in fractional space.

Warning: this assumes the two molecules are in the box described by fractional coords [0, 1]³.

# Arguments
- `dxf::Array{Float64, 1}`: A vector between two atoms in fractional space
"""
@inline function nearest_image!(dxf::Array{Float64, 1})
    for k = 1:3 # loop over components
        @inbounds if abs(dxf[k]) > 0.5
            @inbounds dxf[k] -= sign(dxf[k])
        end
    end
    return nothing
end

"""
    shortest_r² = nearest_r²(xf, yf, box)

Given two points in fractional space, this determines the square of the closest
they can be if the nearest image convention is being used. If the nearest image
is closer than the original it will return the distance between xf and image of yf.

# Arguments
- `xf::Array{Float64, 1}`: The fractional coordinates of point x in 3-space
- `yf::Array{Float64, 1}`: The fractional coordinates of point y in 3-space
- `box::Box`: The box these points are being compared in
"""
@inline function nearest_r²(xf::Array{Float64, 1}, yf::Array{Float64, 1}, box::Box)
    # vector from y to x in fractional coordinate space
    @inbounds dxf = xf - yf
    # apply nearest image convention for periodic boundary conditions
    nearest_image!(dxf)
    # convert to cartesian
    @inbounds dx = box.f_to_c * dxf
    # return r²
    @inbounds return dx[1] * dx[1] + dx[2] * dx[2] + dx[3] * dx[3]
end

"""
    shortest_r = nearest_r(xf, yf, box)

Given two points in fractional space, this determines the closest they can be if
the nearest image convention is being used. If the nearest image is closer than
the original it will return the distance between xf and image of yf.

# Arguments
- `xf::Array{Float64, 1}`: The fractional coordinates of point x in 3-space
- `yf::Array{Float64, 1}`: The fractional coordinates of point y in 3-space
- `box::Box`: The box these points are being compared in
"""
@inline function nearest_r(xf::Array{Float64, 1}, yf::Array{Float64, 1}, box::Box)
    return sqrt(nearest_r²(xf, yf, box))
end
