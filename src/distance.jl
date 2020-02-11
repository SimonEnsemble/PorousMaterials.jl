"""
    nearest_image!(dxf)

Applies the nearest image convention on a vector `dxf` between two atoms in fractional
space; modifies `dxf` for nearest image convention. Fractional coordinates here fall in
[0, 1] so that the box is [0, 1]^3 in fractional space.

Warning: this assumes the two molecules are in the box described by fractional coords [0, 1]³.

# Arguments
- `dxf::Array{Float64}`: A vector between two atoms in fractional space
"""
@inline function nearest_image!(dxf::Array{Float64})
    for i in eachindex(dxf)
        @inbounds if abs(dxf[i]) > 0.5
            @inbounds dxf[i] -= sign(dxf[i])
        end
    end
end

"""
    r = distance(coords, box, i, j, apply_pbc)
    r = distance(atoms, box, i, j, apply_pbc) # atoms i and j
    r = distance(charges, box, i, j, apply_pbc) # atoms i and j

calculate the (Cartesian) distance between particles `i` and `j`.

apply periodic boundary conditions if and only if `apply_pbc` is `true`.

# arguments
- `coords::Coords`: the coordinates (`Frac>:Coords` or `Cart>:Coords`)
- `atoms::Atoms`: atoms
- `charges::charges`: atoms
- `box::Box`: unit cell information
- `i::Int`: index of the first particle
- `j::Int`: Index of the second particle
- `apply_pbc::Bool`: `true` if we wish to apply periodic boundary conditions, `false` otherwise
"""
function distance(coords::Frac, box::Box, i::Int, j::Int, apply_pbc::Bool)
    dxf = coords.xf[:, i] - coords.xf[:, j]
    if apply_pbc
        nearest_image!(dxf)
    end
    return norm(box.f_to_c * dxf)
end

function distance(coords::Cart, box::Box, i::Int, j::Int, apply_pbc::Bool)
    dx = coords.x[:, i] - coords.x[:, j]
    if apply_pbc
        dxf = box.c_to_f * dx
        nearest_image!(dxf)
        return norm(box.f_to_c * dxf)
    else
        return norm(dx)
    end
end

distance(atoms::Atoms, box::Box, i::Int, j::Int, apply_pbc::Bool) = distance(atoms.coords, box, i, j, apply_pbc)
distance(charges::Charges, box::Box, i::Int, j::Int, apply_pbc::Bool) = distance(charges.coords, box, i, j, apply_pbc)

"""
    overlap_flag, overlap_pairs = overlap(frac_coords, box, apply_pbc; tol=0.1)
    overlap_flag, overlap_pairs = overlap(crystal)

determine if any coordinates overlap. here, two coordinates are defined to overlap if their (Cartesian) distance
is less than `tol`.

# Arguments
- `coords::Frac`: the fractional coordinates (`Frac>:Coords`)
- `box::Box`: unit cell information
- `apply_pbc::Bool`: `true` if we wish to apply periodic boundary conditions, `false` otherwise
- `tol::Float64`: tolerance for overlap; if distance between particles less than this, overlap occurs

# Returns
* `overlap_flag::Bool`: `true` if overlap, `false` otherwise
* `overlap_ids::Array{Tuple{Int, Int}, 1}`: ids of coordinate pairs that are overlapping e.g. `[(4, 5), (7, 8)]`
"""
function overlap(coords::Frac, box::Box, apply_pbc::Bool; tol::Float64=0.1)
    overlap_flag = false
    overlap_ids = Array{Tuple{Int, Int}, 1}()

    n = size(coords.xf)[2] # number of coords
    for i = 1:n
        for j = i+1:n
            r = distance(coords, box, i, j, apply_pbc)
            if r < tol
                push!(overlap_ids, (i, j))
                overlap_flag = true
            end
        end
    end
    return overlap_flag, overlap_ids
end
