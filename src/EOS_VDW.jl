
#module EOS

#export compute_properties
#export density,vdw,rho,props
using Roots
using Polynomials

struct vdWMolecule
    a::Float64
    b::Float64
end


# [J/mol K] gas constant
const R = 0.083145
#Stand in for VDW constants
vdw = vdWMolecule(0.2476,0.02661)

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
    #looks at each root and determines if it is real or imag
    i = isreal(polroots[1])
    j = isreal(polroots[2])
    k = isreal(polroots[3])
    #assigns rho to be the real root and then makes it real to get rid of the 0im
    rho = real.(polroots[isreal.(polroots)])

    #Finds fugacity using the derivation from the vander waals
    fug = P .* exp. (- log. (((1 ./ rho) - vdw.b) * P./(R * T))+(vdw.b ./ ((1 ./ rho)-vdw.b) - 2*vdw.a*rho/(R*T)))
    #Dict('D' => rho, 'F' => fug)
end

#end

#using EOS
T = 10000
P = 1

props = compute_properties(T,P,vdw)
println(props)
return props
