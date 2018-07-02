# EOS.jl
# ----------------------------------------------------------------------------------------
# Description

# This script calculates a real gas' compressibility factor and returns a dictionary of
# the real gas' properties such as fugacity, molar volume, and density.
# The script first calculates a real gas' compressibility factor using the polynomial
# form of the Peng-Robinson Equation of State, which requires the gas' critical
# temperature and pressure, and acentric factor, along with user input values for
# temperature and pressure in which the gas resides in, and uses this value to determine
# the real gas' fugacity. Next, it will calculate the gas' density, which is then
# used to find molar volume. These values are then returned to the user

using Polynomials
using DataFrames
using CSV

# Universal gas constant (R). units: m³-bar/(K-mol)
const R = 8.3144598e-5

# Characteristics of a Peng-Robinson gas
struct PengRobinsonGas
    gas::Symbol
    "critical temperature (units: K)"
    Tc::Float64
    "critical pressure (units: bar)"
    Pc::Float64 
    "acentric factor (unitless)"
    ω::Float64
end

# These functions evaluate the Peng-Robinson Equation of State to determine the the compressibility factor(z)
# of a real gas using its critical temperature and pressure, and acentric factor.
a(gas::PengRobinsonGas) = (0.457235 * R ^ 2 * gas.Tc ^ 2) / gas.Pc
b(gas::PengRobinsonGas) = (0.0777961 * R * gas.Tc) / gas.Pc
κ(gas::PengRobinsonGas) = 0.37464 + (1.54226 * gas.ω) - (0.26992 * gas.ω ^ 2)
α(κ::Float64, Tr::Float64) = (1 + κ * (1 - √(Tr))) ^ 2
A(T::Float64, P::Float64, gas::PengRobinsonGas) = α(κ(gas), T / gas.Tc) * a(gas) * P / (R ^ 2 * T ^ 2)
B(T::Float64, P::Float64, gas::PengRobinsonGas) = b(gas) * P / (R * T)

# PREOS is a third degree polynomial. Therefore, there will be three outputs for z.
# When z is equal to 1, the gas will behave as an ideal gas.
# Finding the value of the real number closest to 1 will be the compressibility factor of the real gas.
# By subtracting each number in the 3x1 array by 1 and taking the absolute value, the lowest number can be determined.
# That number is then returned and called out of the function and further calculated with
# user input values of temperature and pressure to determine z (unitless).
# TODO put in Julia documentation format
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

# Calculating for fugacity coefficient from an integration (bar)
function calculate_ϕ(gas::PengRobinsonGas, T::Float64, P::Float64)
    z = compressibility_factor(gas, T, P)
    log_ϕ = z - 1.0 - log(z - B(T, P, gas)) +
            - A(T, P, gas) / (√8 * B(T, P, gas)) * log(
            (z + (1 + √2) * B(T, P, gas)) / (z + (1 - √(2)) * B(T, P, gas)))
    return exp(log_ϕ) #fc_1 - fc_2 * fc_3)
end

# Returns a dictionary of gas properties based on user input for T and P
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

function PengRobinsonGas(gas::Symbol)
    # Read in file with gas parameters using dataframes
    # TODO change to Porousmaterials.PATH_TO_DATA
    df = CSV.read("data/PengRobinsonGasProps.csv")
    if ! (string(gas) in df[:gas])
        error(@sprintf("Gas %s properties not found in %sPengRobinsonGasProps.csv", gas, PATH_TO_DATA))
    end
    # Query values for gas
    Tc = df[df[:gas].== string(gas), Symbol("Tc(K)")][1]
    Pc = df[df[:gas].== string(gas), Symbol("Pc(bar)")][1]
    ω = df[df[:gas].== string(gas), Symbol("acentric_factor")][1]
    return PengRobinsonGas(gas, Tc, Pc, ω)
end

function Base.show(io::IO, gas::PengRobinsonGas)
    println(io, "Gas species: ", gas.gas)
    println(io, "Critical temperature (K): ", gas.Tc)
    println(io, "Critical pressure (bar): ", gas.Pc)
    println(io, "Acenteric factor: ", gas.ω)
end
