<<<<<<< HEAD
# Begin Module EOS
module EOS

# calls and exports all packages and functions used in module
export compute_properties, VDWGas
export vdw, rho
=======
module EOS

export compute_properties, VDWGas
export density,vdw,rho,props
>>>>>>> 1fb0f528849d859b5fd78da7a9293de719b0f292
using Roots
using Polynomials
using CSV
using DataFrames

<<<<<<< HEAD
# Assigns R to be a constant
const R = 0.083145  # [L bar/ mol K]

=======
# [J/mol K] gas constant
const R = 0.083145
>>>>>>> 1fb0f528849d859b5fd78da7a9293de719b0f292

# Assigns type of variable used in calculations. Constructs vdWMolecule array.
struct vdWMolecule
  a::Float64
  b::Float64
end

<<<<<<< HEAD
"""

The VDWGas function is used to identify and gather the vanderwaals gas constants
specific to the molecule being analyzed. These values are stored in the function
to be used in later calculations.

"""

function VDWGas(gas::Symbol)

        if ! isfile("data//vdw_constants.csv")
        error(@sprintf("No file 'vdw_constants.csv' exists. This file is needed to determine specific Vanderwaals constants"))
        end

        gas = string(gas)
        vdwfile = CSV.read("data//vdw_constants.csv")
        A = vdwfile[vdwfile[:, 1] .== gas, 2]
        B = vdwfile[vdwfile[:, 1] .== gas, 3]
        vdw = vdWMolecule(A[1], B[1])

end

"""
The compute_properties function is used to calculate the density and fugacity of
a specified molecule based on a known temperature and pressure using the
Vanderwaals gas equation. The gas constants gathered above will be utilized in
this function.

"""

function compute_properties(T, P, vdw)

        # calculates coefficients used in polynomial equation derived from
        # Vanderwaals
        A = -P
        B = (P * vdw.b + R * T)
        C = -vdw.a
        D = vdw.a * vdw.b

        # Creates a third order polynomial equation from the coefficients
        # calculated above
        pol = Poly([A, B, C, D])

        # Determines the roots of the equation above
        polroots = roots(pol)

        # Returns the value of the real root(s) and assigns this value to be the
        # Density of the molecule. (Refer to compute_properties function for
        # Vanderwaals derivation results)
        rho = real.(polroots[isreal.(polroots)])

        # Calculates fugacity of desired molecule from calculations and
        # functions above.
        fug = P .* exp. (- log. (((1 ./ rho) - vdw.b) * P./(R * T))+(vdw.b ./ ((1 ./ rho)-vdw.b) - 2*vdw.a*rho/(R*T)))

        # Enters fugacity and density into a dictionary to be displayed once
        # the function is complete.
        Dict("Density" => rho, "Fugacity" => fug)
end

=======
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
        Dict('D' => rho, 'F' => fug)
end

>>>>>>> 1fb0f528849d859b5fd78da7a9293de719b0f292
end
