using Polynomials
using DataFrames
using CSV
using Roots

"""
Calculates a Peng-Robinson gas' compressibility factor and returns a dictionary of its properties.
Properties include density, molar volume, and fugacity.
Compressibility factor found using polynomial form of the Peng-Robinson Equation of State.
Compressibility factor used to find density, molar volume, and fugacity.
"""

# Universal gas constant (R). units: m³-bar/(K-mol)
const R = 8.3144598e-5

# Data structure stating characteristics of a Peng-Robinson gas
struct PengRobinsonGas
    gas::Symbol
    # Peng-Robinson Gas
    Tc::Float64
    # Critical temperature (units: Kelvin)
    Pc::Float64
    # Critical pressure (units: bar)
    ω::Float64
    # Acentric factor (units: unitless)
end

#Characteristics of a Van der Waals gas
struct vdWMolecule
  #VDW constant a (units: m⁶bar/mol)
  a::Float64
  #VDW constant b (unites: m³/mol)
  b::Float64
  gas::Symbol
end

"""
Evaluates the Peng-Robinson Equation of State to determine the compressibility factor(z)
of a real gas using its critical temperature and pressure, and acentric factor.

# Arguments
- `gas::PengRobinsonGas`: Peng-Robinson gas structure
- `κ::Float64`: Definition of α in relation to the acentric factor
- `Tr::Float64`: Reduced temperature
- `T::Float64`: Temperature given in Kelvin
- `P::Float64`: Pressure given in bar
"""

a(gas::PengRobinsonGas) = (0.457235 * R ^ 2 * gas.Tc ^ 2) / gas.Pc
b(gas::PengRobinsonGas) = (0.0777961 * R * gas.Tc) / gas.Pc
κ(gas::PengRobinsonGas) = 0.37464 + (1.54226 * gas.ω) - (0.26992 * gas.ω ^ 2)
α(κ::Float64, Tr::Float64) = (1 + κ * (1 - √(Tr))) ^ 2
A(T::Float64, P::Float64, gas::PengRobinsonGas) = α(κ(gas), T / gas.Tc) * a(gas) * P / (R ^ 2 * T ^ 2)
B(T::Float64, P::Float64, gas::PengRobinsonGas) = b(gas) * P / (R * T)

# Calculates three outputs for compressibility factor using the polynomial form of
# the Peng-Robinson Equation of State. Filters for only real roots and returns the
# root closest to unity.

function compressibility_factor(gas::PengRobinsonGas, T::Float64, P::Float64)
    # construct cubic polynomial in z
    p = Poly([-(A(T, P, gas) * B(T, P, gas) - B(T, P, gas) ^ 2 - B(T, P, gas) ^ 3),
              A(T, P, gas) - 2 * B(T, P, gas) - 3 * B(T, P, gas) ^ 2,
              -(1.0 - B(T, P, gas)),
              1.0])
    # solve for the roots of the cubic polynomial
    z_roots = roots(p)
    # select real roots only.
    z_factor = z_roots[isreal.(z_roots)]
    # find the index of the root that is closest to unity
    id_closest_to_unity = indmin(abs.(z_factor - 1.0))
    # return root closest to unity.
    return real(z_factor[id_closest_to_unity])
end

# Calculating for fugacity coefficient from an integration (bar).
function calculate_ϕ(gas::PengRobinsonGas, T::Float64, P::Float64)
    z = compressibility_factor(gas, T, P)
    log_ϕ = z - 1.0 - log(z - B(T, P, gas)) +
            - A(T, P, gas) / (√8 * B(T, P, gas)) * log(
            (z + (1 + √2) * B(T, P, gas)) / (z + (1 - √(2)) * B(T, P, gas)))
    return exp(log_ϕ)
end

"""
    calculate_properties = (gas, T, P)

Use equation of state to calculate density, fugacity, and molar volume of a real gas at a
given temperature and pressure.

# Arguments
- `gas::PengRobinsonGas`: Peng-Robinson gas structure
- `T::Float64`: Temperature given in Kelvin
- `P::Float64`: Pressure given in bar

# Returns
- `prop_dict:: `: Dictionary of Peng-Robinson gas properties
"""
function calculate_properties(gas::PengRobinsonGas, T::Float64, P::Float64; verbose::Bool=true)
    # Compressbility factor (unitless)
    z = compressibility_factor(gas, T, P)
    # Density (mol/m^3)
    ρ = P / (z * R * T)
    # Molar volume (L/mol)
    Vm = 1.0 / ρ * 1000.0
    # Fugacity (bar)
    ϕ = calculate_ϕ(gas, T, P)
    f = ϕ * P
    # Prints a dictionary holding values for compressibility factor, molar volume, density, and fugacity.
    prop_dict = Dict("compressibility factor" => z, "molar volume (L/mol)"=> Vm ,
                     "density (mol/m³)" => ρ, "fugacity (bar)" => f,
                     "fugacity coefficient" => ϕ)
    if verbose
        @printf("%s properties at T = %f K, P = %f bar:\n", gas.gas, T, P)
        for (property, value) in prop_dict
            println(property * ": ", value)
        end
    end
    return prop_dict
end


""" 
    gas = PengRobinsonGas(gas)

Reads in properties file in the directory `PorousMaterials.PATH_TO_DATA * "PengRobinsonGasProps.csv"`
with gas parameters using dataframes and queries values for gas.

# Returns
- `PengRobinsonGas::struct`: Data structure containing Peng-Robinson gas parameters.
"""

#Function to solve for fugacity and rho using
function calculate_properties(gas::vdWMolecule, T::Float64, P::Float64)

        A = -P
        B = (P * gas.b + R * T)
        C = -gas.a
        D = gas.a * gas.b

        #Creates a polynomial for the vdw cubic function
        pol = Poly([A, B, C, D])
        #finds roots of that polynomial
        polroots = roots(pol)
        #assigns rho to be the real root and then makes it real to get rid of the 0im
        rho = real.(polroots[isreal.(polroots)])

        #specifies that molar volume is the reciprocal of the density
        # In units of L/mol
        vm = (1./ rho) * 1000
        #specifies the compressibility factor
        z = (P * (1./ rho))./ (R * T)


        #Finds fugacity using the derivation from the vander waals
        fug = P .* exp. (- log. (((1 ./ rho) - gas.b) * P./(R * T))+(gas.b ./ ((1 ./ rho)-gas.b) - 2*gas.a*rho/(R*T)))
        #defines the fugacity coefficient as fugacity over pressure
        ϕ = fug ./ P

        prop_dict = Dict("Density (mol/m³)" => rho, "Fugacity (bar)" => fug, "Molar Volume (L/mol)" => vm, "Fugacity Coefficient" => ϕ, "Compressibility Factor" => z )
end

function VDWGas(gas::Symbol)

        gas = string(gas)
        vdwfile = CSV.read("C:\\Users\\Caleb\\Programming Work\\Julia\\Actual\\vdw_constants.csv")
        if ! (:gas in vdwfile(:molecule))
              error(@sprintf("Gas %s properties not found in %sPengRobinsonGasProps.csv", gas, PATH_TO_DATA))
        end
        A = vdwfile[vdwfile[:molecule].== gas, Symbol("a(m6bar/mol)")]
        B = vdwfile[vdwfile[:molecule].== gas, Symbol("b(m3/mol)")]
        return vdWMolecule(A[1], B[1], gas)

end


function PengRobinsonGas(gas::Symbol)
    df = CSV.read(PATH_TO_DATA * "PengRobinsonGasProps.csv")
    if ! (string(gas) in df[:gas])
        error(@sprintf("Gas %s properties not found in %sPengRobinsonGasProps.csv", gas, PATH_TO_DATA))
    end
    Tc = df[df[:gas].== string(gas), Symbol("Tc(K)")][1]
    Pc = df[df[:gas].== string(gas), Symbol("Pc(bar)")][1]
    ω = df[df[:gas].== string(gas), Symbol("acentric_factor")][1]
    return PengRobinsonGas(gas, Tc, Pc, ω)
end

# Prints resulting values for Peng-Robinson gas properties
function Base.show(io::IO, gas::PengRobinsonGas)
    println(io, "Gas species: ", gas.gas)
    println(io, "Critical temperature (K): ", gas.Tc)
    println(io, "Critical pressure (bar): ", gas.Pc)
    println(io, "Acenteric factor: ", gas.ω)
end

function Base.show(io::IO, gas::vdWMolecule)
    println(io, "Gas species: ", gas.gas)
    println(io, "Van der Waals constant a (m⁶bar/mol): ", gas.a)
    println(io, "Van der Waals constant b (m³/mol): ", gas.b)
end
