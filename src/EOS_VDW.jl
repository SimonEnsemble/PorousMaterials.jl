module EOS

export compute_properties, VDWGas
export density,vdw,rho,props
using Roots
using Polynomials
using CSV
using DataFrames

# [J/mol K] gas constant
const R = 0.083145

struct vdWMolecule
  a::Float64
  b::Float64
end

#Reads VDW CSV and returns the VDW constants for later use
#STILL NEED TO FIGURE OUT HOW TO READ FILE WITHOUT SPECIFIYING PATH DIRECTLY
#IE CSV.read("vdw_constants")
function VDWGas(gas::Symbol)

        gas = string(gas)
        vdwfile = CSV.read("C:\\Users\\Caleb\\Programming Work\\Julia\\Actual\\vdw_constants.csv")
        A = vdwfile[vdwfile[:,1] .== gas,2]
        B = vdwfile[vdwfile[:,1] .== gas,3]
        vdw = vdWMolecule(A[1],B[1])

end

#Function to solve for fugacity and rho using
function compute_properties(T,P,vdw)

        A = -P
        B = (P * vdw.b + R * T)
        C = -vdw.a
        D = vdw.a * vdw.b

        #Creates a polynomial for the vdw cubic function
        pol = Poly([A,B,C,D])
        #finds roots of that polynomial
        polroots = roots(pol)
        #assigns rho to be the real root and then makes it real to get rid of the 0im
        rho = real.(polroots[isreal.(polroots)])

        #Finds fugacity using the derivation from the vander waals
        fug = P .* exp. (- log. (((1 ./ rho) - vdw.b) * P./(R * T))+(vdw.b ./ ((1 ./ rho)-vdw.b) - 2*vdw.a*rho/(R*T)))
        #Dict('D' => rho, 'F' => fug)
end

end
