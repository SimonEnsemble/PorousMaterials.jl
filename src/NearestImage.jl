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

@inline function nearest_r(xf::Array{Float64, 1}, yf::Array{Float64, 1}, box::Box)
    return sqrt(nearest_r²(xf, yf, box))
end
