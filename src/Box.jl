"""
    box = Box(a, b, c, α, β, γ, volume, f_to_c, c_to_f, reciprocal_lattice)
    box = Box(a, b, c, α, β, γ)
    box = Box(f_to_c)

Data structure to describe a unit cell box (Bravais lattice) and convert between
fractional and Cartesian coordinates.

# Attributes
- `a,b,c::Float64`: unit cell dimensions (units: Angstroms)
- `α,β,γ::Float64`: unit cell angles (units: radians)
- `Ω::Float64`: volume of the unit cell (units: cubic Angtroms)
- `f_to_c::Array{Float64,2}`: the 3x3 transformation matrix used to map fractional
coordinates to cartesian coordinates. The columns of this matrix define the unit cell
axes. Columns are the vectors defining the unit cell box. units: Angstrom
- `c_to_f::Array{Float64,2}`: the 3x3 transformation matrix used to map Cartesian
coordinates to fractional coordinates. units: inverse Angstrom
- `reciprocal_lattice::Array{Float64, 2}`: the *rows* are the reciprocal lattice vectors.
This choice was made (instead of columns) for speed of Ewald Sums.
"""
struct Box
    a::Float64
    b::Float64
    c::Float64

    α::Float64
    β::Float64
    γ::Float64

    Ω::Float64

    f_to_c::Array{Float64, 2}
    c_to_f::Array{Float64, 2}

    reciprocal_lattice::Array{Float64, 2}
end

Base.broadcastable(b::Box) = Ref(b)

"""
    r = reciprocal_lattice(f_to_c)

Given the fractional to Cartesian transformation matrix, compute the reciprocal lattice
vectors, which are the rows of the returned matrix `r`.

# Arguments
- `f_to_c::Array{Float64,2}`: the 3x3 transformation matrix used to map fractional
coordinates to cartesian coordinates. The columns of this matrix define the unit cell
axes. units: Angstrom

# Returns
- `r::Array{Float64, 2}`: Reciprocal lattice vectors in a matrix format, where the rows
are the reciprocal lattice vectors.
"""
function reciprocal_lattice(f_to_c::Array{Float64, 2})
    # the unit cell vectors are the columns of f_to_c
    a₁ = f_to_c[:, 1]
    a₂ = f_to_c[:, 2]
    a₃ = f_to_c[:, 3]

    r = zeros(Float64, 3, 3)
    r[1, :] = 2 * π * cross(a₂, a₃) / dot(a₁, cross(a₂, a₃))
    r[2, :] = 2 * π * cross(a₃, a₁) / dot(a₂, cross(a₃, a₁))
    r[3, :] = 2 * π * cross(a₁, a₂) / dot(a₃, cross(a₁, a₂))
    return r
end

# Automatically calculates Ω, f_to_c, and c_to_f for `Box` data structure based on axes lenghts and angles.
function Box(a::Float64, b::Float64, c::Float64,
             α::Float64, β::Float64, γ::Float64)
    # unit cell volume (A³)
    Ω = a * b * c * sqrt(1 - cos(α) ^ 2 - cos(β) ^ 2 - cos(γ) ^ 2 + 2 * cos(α) * cos(β) * cos(γ))
    # matrices to map fractional coords <--> Cartesian coords
    f_to_c = [[a, 0, 0] [b * cos(γ), b * sin(γ), 0] [c * cos(β), c * (cos(α) - cos(β) * cos(γ)) / sin(γ), Ω / (a * b * sin(γ))]]
    c_to_f = [[1/a, 0, 0] [-cos(γ) / (a * sin(γ)), 1 / (b * sin(γ)), 0] [b * c * (cos(α) * cos(γ) - cos(β)) / (Ω * sin(γ)), a * c * (cos(β) * cos(γ) - cos(α)) / (Ω * sin(γ)), a * b * sin(γ) / Ω]]
    # the columns of f_to_c are the unit cell axes
    r = reciprocal_lattice(f_to_c)

    @assert f_to_c * c_to_f ≈ Matrix{Float64}(I, 3, 3)
    @assert isapprox(r, 2.0 * π * inv(f_to_c))

    return Box(a, b, c, α, β, γ, Ω, f_to_c, c_to_f, r)
end

# Constructs a `Box` from the fractional to Cartesian coordinate transformation matrix.
# The columns of this matrix are the unit cell axes. Units must be in Å.
function Box(f_to_c::Array{Float64, 2})
    # unit cell volume (A³) is the determinant of the matrix
    Ω = det(f_to_c)

    # unit cell dimensions are lengths of the unit cell axes (Å)
    a = norm(f_to_c[:, 1])
    b = norm(f_to_c[:, 2])
    c = norm(f_to_c[:, 3])

    # c_to_f is the inverse
    c_to_f = inv(f_to_c)

    # angles (radians)
    α = acos(dot(f_to_c[:, 2], f_to_c[:, 3]) / (b * c))
    β = acos(dot(f_to_c[:, 1], f_to_c[:, 3]) / (a * c))
    γ = acos(dot(f_to_c[:, 1], f_to_c[:, 2]) / (a * b))

    # the columns of f_to_c are the unit cell axes
    r = reciprocal_lattice(f_to_c)

    return Box(a, b, c, α, β, γ, Ω, f_to_c, c_to_f, r)
end

"""
    unit_cube = UnitCube()

This function generates a unit cube, each side is 1.0 Angstrom long, and all the
corners are right angles.
"""
UnitCube() = Box(1.0, 1.0, 1.0, π/2, π/2, π/2)

"""
    new_box = replicate(original_box, repfactors)

Replicates a `Box` in positive directions to construct a new `Box` representing a supercell.
The `original_box` is replicated according to the factors in `repfactors`.
Note `replicate(original_box, repfactors=(1, 1, 1))` returns same `Box`.
The new fractional coordinates as described by `f_to_c` and `c_to_f` still ∈ [0, 1].

# Arguments
- `original_box::Box`: The box that you want to replicate
- `repfactors::Tuple{Int, Int, Int}`: The factor you want to replicate the box by

# Returns
- `box::Box`: Fully formed Box object
"""
function replicate(box::Box, repfactors::Tuple{Int, Int, Int})
    return Box(box.a * repfactors[1], box.b * repfactors[2], box.c * repfactors[3],
               box.α, box.β, box.γ)
end

"""
    write_vtk(box, filename; verbose=true, center_at_origin=false)
    write_vtk(framework)

Write a `Box` to a .vtk file for visualizing e.g. the unit cell boundary of a crystal.
If a `Framework` is passed, the `Box` of that framework is written to a file that is the
same as the crystal structure filename but with a .vtk extension.

Appends ".vtk" extension to `filename` automatically if not passed.

# Arguments
- `box::Box`: a Bravais lattice
- `filename::AbstractString`: filename of the .vtk file output (absolute path)
- `framework::Framework`: A framework containing the crystal structure information
- `center_at_origin::Bool`: center box at origin if true. if false, the origin is the corner of the box.
"""
function write_vtk(box::Box, filename::AbstractString; verbose::Bool=true, 
                   center_at_origin::Bool=false)
    if ! occursin(".vtk", filename)
        filename *= ".vtk"
    end
    vtk_file = open(filename, "w")

    @printf(vtk_file, "# vtk DataFile Version 2.0\nunit cell boundary\n
                       ASCII\nDATASET POLYDATA\nPOINTS 8 double\n")
    
    x_shift = zeros(3)
    if center_at_origin
        x_shift = box.f_to_c * [0.5, 0.5, 0.5]
    end

    # write points on boundary of unit cell
    for i = 0:1
        for j = 0:1
            for k = 0:1
                xf = [i, j, k] # fractional coordinates of corner
                cornerpoint = box.f_to_c * xf - x_shift
                @printf(vtk_file, "%.3f %.3f %.3f\n",
                        cornerpoint[1], cornerpoint[2], cornerpoint[3])
            end
        end
    end

    # define connections
    @printf(vtk_file, "LINES 12 36\n2 0 1\n2 0 2\n2 1 3\n2 2 3\n2 4 5\n2 4 6\n2 5 7\n2 6 7\n2 0 4\n2 1 5\n2 2 6\n2 3 7\n")
    close(vtk_file)
    if verbose
        println("See ", filename)
    end
end

"""
    inside_box = inside(x, box) # true or false

Determine whether a Cartesian vector `x` lays inside a `Box`. This works by computing the 
fractional coordinates of vector `x` and ensuring each lie within the interval `[0, 1]`.
"""
function inside(x::Array{Float64, 1}, box::Box)
    xf = box.c_to_f * x
    return all(xf .<= 1.0) && all(xf .>= 0.0)
end

function Base.show(io::IO, box::Box)
    println(io, "Bravais unit cell of a crystal.")
    @printf(io, "\tUnit cell angles α = %f deg. β = %f deg. γ = %f deg.\n",
        box.α * 180.0 / π, box.β * 180.0 / π, box.γ * 180.0 / π)
    @printf(io, "\tUnit cell dimensions a = %f Å. b = %f Å, c = %f Å\n",
        box.a, box.b, box.c)
    @printf(io, "\tVolume of unit cell: %f Å³\n", box.Ω)
end

function Base.isapprox(box1::Box, box2::Box; rtol::Real=sqrt(eps()))
    return (isapprox(box1.a, box2.a, rtol=rtol) &&
            isapprox(box1.b, box2.b, rtol=rtol) &&
            isapprox(box1.c, box2.c, rtol=rtol) &&
            isapprox(box1.α, box2.α, rtol=rtol) &&
            isapprox(box1.β, box2.β, rtol=rtol) &&
            isapprox(box1.γ, box2.γ, rtol=rtol) &&
            isapprox(box1.Ω, box2.Ω, rtol=rtol) &&
            isapprox(box1.f_to_c, box2.f_to_c, rtol=rtol) &&
            isapprox(box1.c_to_f, box2.c_to_f, rtol=rtol) &&
            isapprox(box1.reciprocal_lattice, box2.reciprocal_lattice, rtol=rtol))
end
