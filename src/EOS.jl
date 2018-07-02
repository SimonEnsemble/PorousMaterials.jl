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

using PorousMaterials
using Polynomials
using DataFrames
using CSV

# This specifies the value of the ideal gas constant (R).
const R = 8.314

# Characteristics of a real gas
struct PengRobinsonGas
    gas::Symbol
    Tc::Float64
    Pc::Float64
    ω::Float64
end

# These functions evaluate the Peng-Robinson Equation of State to determine the the compressibility factor(z)
# of a real gas using its critical temperature and pressure, and acentric factor.
a(gas::PengRobinsonGas) = (0.457235 * R ^ 2 * (gas.Tc) ^ 2) / (gas.Pc)
b(gas::PengRobinsonGas) = (0.0777961 * R * (gas.Tc)) / (gas.Pc)
κ(gas::PengRobinsonGas) = 0.37464 + (1.54226 * (gas.ω)) - (0.26992 * (gas.ω) ^ 2)
Tr(gas::PengRobinsonGas, T::Float64) = T / (gas.Tc)
α(κ::Float64, Tr::Float64) = (1 + κ * (1 - √(Tr))) ^ 2
A(T::Float64, P::Float64) = α(κ(gas), Tr(gas, T) * a(gas) * P) / (R ^ 2 * T ^ 2)
B(T::Float64, P::Float64) = (b(gas) * P)/(R * T)

# PREOS is a third degree polynomial. Therefore, there will be three outputs for z.
# When z is equal to 1, the gas will behave as an ideal gas.
# Finding the value of the real number closest to 1 will be the compressibility factor of the real gas.
# By subtracting each number in the 3x1 array by 1 and taking the absolute value, the lowest number can be determined.
# That number is then returned and called out of the function and further calculated with
# user input values of temperature and pressure to determine z (unitless).
function compressibility_factor(gas::PengRobinsonGas, T::Float64, P::Float64)
    z_roots = roots(Poly([-(A(T, P) * B(T, P) - B(T, P) ^ 2 - B(T, P) ^ 3), (A(T, P) - 2 * B(T, P) - 3 * B(T, P) ^ 2), -(1 - B(T, P)), 1]))
    z_factor = z_roots[isreal.(z_roots)]
    return minimum(abs.(z_factor - 1))
end

# Calculating for density with compressibility factor (mol/m^3)
calculate_ρ(T::Float64, P::Float64) = P / ( compressibility_factor(gas, T, P) * R * T)

# Calculating for molar volume with compressibility factor (L/mol)
molar_volume(T::Float64, P::Float64) = 1 / ρ(T, P) * 1000

# Calculating for fugacity coefficient from an integration (bar)
function calculate_ϕ(gas::PengRobinsonGas, T::Float64, P::Float64)
    fc_1 = compressibility_factor(gas, T, P) - 1 - log(e, compressibility_factor(gas, T, P) - B(T, P))
    fc_2 = A(T, P) / (sqrt(8) * B(T, P))
    fc_3 = log(e, (compressibility_factor(gas, T, P) + (1 + √(2)) * B(T, P)) / (compressibility_factor(gas, T, P) - (1 - √(2)) * B(T, P)))
    solved_ϕ = exp(fc_1 - fc_2 * fc_3)
    return solved_ϕ
end

# Returns a dictionary of gas properties based on user input for T and P
function calculate_properties(gas::PengRobinsonGas, T::Float64, P::Float64; verbose::Bool=true)
    # Compressbility factor (unitless)
    z = compressibility_factor(gas, T, P)
    # Molar volume (L/mol)
    Vm = molar_volume(T, P)
    # Density (mol/m^3)
    ρ = calculate_ρ(T, P)
    # Fugacity (bar)
    f = calculate_ϕ(gas, T, P) * P
    # Prints a dictionary holding values for compressibility factor, molar volume, density, and fugacity.
    prop_dict = Dict("compressibility factor"=>z, "molar volume (m^3/mol)"=>Vm, "density (mol/m^3)"=>ρ, "fugacity (bar)"=>f)
    if verbose
        println(prop_dict)
    end
    return prop_dict
end

function PengRobinsonGas(gas::Symbol)
    # Read in file with gas parameters using dataframes
    df = CSV.read("PengRobinsonGasProps.csv")
    if ! (string(gas) in df[:gas])
        error(@sprintf("Gas %s properties not found in PengRobinsonGasProps.csv", gas))
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
