"""All things crystal structures"""
module Crystal 

export Framework, readcssr, replicate_to_xyz

global PATH_TO_STRUCTURE_FILES = homedir() * "/Box/School/Code/StructureFiles"

"""
    framework = Framework(a, b, c, α, β, γ, N, atoms, f_coords, f_to_c, c_to_f)

Data structure for a 3D crystal structure.

# Arguments
- `a,b,c::Float64`: unit cell dimensions (units: Angstroms)
- `α,β,γ::Float64`: unit cell angles (units: degrees)
- `n_atoms::Int64`: number of atoms in a unit cell
- `Ω::Float64`: the volume of the unit cell
- `atoms::Array{String,1}`: list of atoms composing crystal unit cell, in strict order
- `f_coords::Array{Float64,2}`: a 2D array of fractional coordinates of the atoms in order
- `f_to_c::Array{Float64,2}`: is a 3x3 matrix used to convert fractional coordinates to cartesian coordinates
- `c_to_f::Array{Float64,2}`: is a 3x3 matrix used to convert cartesian coordinates to fractional coordinates
"""
struct Framework
    a::Float64
    b::Float64
    c::Float64

    α::Float64
    β::Float64
    γ::Float64

    Ω::Float64

    n_atoms::Int64
    atoms::Array{String, 1}
    f_coords::Array{Float64, 2}

    f_to_c::Array{Float64, 2}
    c_to_f::Array{Float64, 2}
end

"""
    framework = readcssr("filename.cssr")

Read a .cssr file and construct a Framework object
"""
function readcssr(cssrfilename::String)
    f = open(cssrfilename,"r")
    lines = readlines(f)
    N = length(lines)-4
    a,b,c,α,β,γ = Array{Float64}(6)
    x = Array{Float64}(N)
    y = similar(x)
    z = similar(x)
    atoms = Array{String}(N)
    for (i,line) in enumerate(lines)
        str = split(line)
        if (i == 1)
            a,b,c = map(x->parse(Float64,x),str)
        elseif (i == 2)
            α,β,γ = map(x->parse(Float64,x),str[1:3])
        elseif (i > 4)
            atoms[i-4] = str[2]
            x[i-4],y[i-4],z[i-4] = map(x->parse(Float64,x),str[3:5])
        end
    end
    close(f)
    Ω = a*b*c*sqrt(1-cos(α)^2-cos(β)^2-cos(γ)^2+2*cos(α)*cos(β)*cos(γ))
    M = [[a, b*cos(γ), c*cos(β)] [0, b*sin(γ), c*(cos(α)-cos(β)*cos(γ))/sin(γ)] [0, 0, Ω/(a*b*sin(γ))]]
    return Framework(a,b,c,α,β,γ,N,Ω,atoms,([x y z]),M)
end

"""
    replicate_to_xyz(framework, xyzfilename, comment="")

Write a .xyz file from a Framework object. Write an optional comment to the .xyz file if desired.
Extend the structure in the x-,y- or z-direction by changing nx, ny or nz respectively. 
A value of 1 replicates the structure once in the desired direction
"""
function replicate_to_xyz(framework::Framework, xyzfilename::String; comment::String="", nx = 0, ny = 0, nz = 0)
    f = open(xyzfilename,"w")
    @printf(f,"%d\n%s\n", framework.n_atoms*(nx+1)*(ny+1)*(nz+1),comment)
    c_coords = similar(framework.f_coords)
    fcoord_temp = similar(framework.f_coords)
    
    for i = 0:nx, j = 0:ny, k = 0:nz
        fcoord_temp = framework.f_coords .+ [i j k]
        c_coords = framework.f_to_c*fcoord_temp'
        for ii = 1:size(c_coords,2)
            @printf(f,"%s\t%.4f\t%.4f\t%.4f\n",framework.atoms[ii],c_coords[1,ii],c_coords[2,ii],c_coords[3,ii])
        end
    end
    close(f)
end

end # end module
