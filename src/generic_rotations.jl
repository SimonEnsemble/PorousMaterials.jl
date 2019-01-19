"""
    R = rotation_matrix(θ, u, assume_unit_vector=false) # 3 by 3 rotation matrix, angle θ about vector u
    R = rotation_matrix(θ, dim) # 3 by 3 rotation matrix, angle θ about axis `dim`

Determine the 3D rotation matrix to rotate an angle θ (radians) about axis `u`.

See [Wikipedia](https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle).

# Arguments
- `θ::Float64`: angle to rotate about an axis, in radians
- `u::Array{Float64, 1}`: axis about which to rotate
- `dim::Int`: 1, 2, 3 for rotation about x-, y-, or z-axis, respectively.
- `assume_unit_vector::Bool`: assume `u` is a unit vector; otherwise, `u` will be normalized
internal to this function.

# Returns
- `R::Array{Float64, 2}`: 3D rotation matrix. so `R * x` will rotate vector `x` as desired.
"""
function rotation_matrix(θ::Float64, u::Array{Float64, 1}; assume_unit_vector::Bool=false)
    if ! assume_unit_vector
        u = u / norm(u)
    end

    c = cos(θ) # for speed pre-compute these
    s = sin(θ)

    R = [c + u[1] ^ 2 * (1.0 - c)            u[1] * u[2] * (1.0 - c) - u[3] * s   u[1] * u[3] * (1.0 - c) + u[2] * s;
         u[2] * u[1] * (1.0 - c) + u[3] * s  c + u[2] ^ 2 * (1.0 - c)             u[2] * u[3] * (1.0 - c) - u[1] * s;
         u[3] * u[1] * (1.0 - c) - u[2] * s  u[3] * u[2] * (1.0 - c) + u[1] * s   c + u[3] ^ 2 * (1.0 - c)]

    return R
end

function rotation_matrix(θ::Float64, dim::Int)
    c = cos(θ) # for speed pre-compute these
    s = sin(θ)
    
    # see https://en.wikipedia.org/wiki/Rotation_matrix#Basic_rotations
    if dim == 1
        return [1.0 0.0 0.0;
                0.0  c  -s;
                0.0  s   c]
    elseif dim == 2
        return [ c  0.0  s;
                0.0 1.0 0.0;
                -s  0.0   c]
    elseif dim == 3
        return [ c  -s  0.0;
                 s   c  0.0;
                0.0 0.0 1.0]
    else
        error("dim must be 1, 2, or 3\n")
    end
end
