"""
    Frame = Framework(a,b,c,α,β,γ,N,Ele,f_coords,convM)

Represents a 3D crystal structure read through a .cssr file.

#Arguments
* a,b,c::Float64 : are unit cell dimensions
* α,β,γ::Float64 : are unit cell angles
* N::Int64 : is number of atoms
* Ele::Array{String,1} : is a vector of elements in order
* f_coords::Array{Float64,2} : is a 2D array of fractional coordinates of the atoms in order
* convM::Array{Float64,2} : is a 3x3 matrix used to convert fractional coordinates to cartesian coordinates

"""
struct Framework
    a::Float64
    b::Float64
    c::Float64
    α::Float64
    β::Float64
    γ::Float64
    N::Int64
    Ele::Array{String,1}
    f_coords::Array{Float64,2}
    convM::Array{Float64,2}
end

"""
    Frame = readcssr("filename.cssr")

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
    Ele = Array{String}(N)
    for (i,line) in enumerate(lines)
        str = split(line)
        if (i == 1)
            a,b,c = map(x->parse(Float64,x),str)
        elseif (i == 2)
            α,β,γ = map(x->parse(Float64,x),str[1:3])
        elseif (i > 4)
            Ele[i-4] = str[2]
            x[i-4],y[i-4],z[i-4] = map(x->parse(Float64,x),str[3:5])
        end
    end
    close(f)
    Ω = a*b*c*sqrt(1-cos(α)^2-cos(β)^2-cos(γ)^2+2*cos(α)*cos(β)*cos(γ))
    M = [[a, b*cos(γ), c*cos(β)] [0, b*sin(γ), c*(cos(α)-cos(β)*cos(γ))/sin(γ)] [0, 0, Ω/(a*b*sin(γ))]]
    return Framework(a,b,c,α,β,γ,N,Ele,([x y z]),M)
end

HKUST = readcssr("HKUST-1.cssr");

"""
    replicate_to_xyz(framework, xyzfilename, comment="")

Write a .xyz file from a Framework object. Write an optional comment to the .xyz file if desired.
"""
function replicate_to_xyz(framework::Framework, xyzfilename::String, comment::String="")
    n_x = 1
    n_y = 2
    n_z = 1
    N = framework.N*n_x*n_y*n_z
    Ele = framework.Ele
    c_coords = framework.convM*framework.f_coords'
    f = open(xyzfilename,"w")
    write(f, "$N\n")
    write(f, "$comment\n")
    for i = 1:n_x, j = 1:n_y, k = 1:n_z
        for ii = 1:size(c_coords,2)
            str = "$(Ele[ii])\t$(c_coords[1,ii])\t$(c_coords[2,ii])\t$(c_coords[3,ii])\n"
            write(f,str)
        end
    end
    close(f)
end

replicate_to_xyz(HKUST,"testfile1.xyz","Test run")


