using Polynomials
using DataFrames
using CSV

"""
Calculates a Peng-Robinson gas' compressibility factor and returns a dictionary of its properties.
Properties include density, molar volume, and fugacity.
Compressibility factor found using polynomial form of the Peng-Robinson Equation of State.
Compressibility factor used to find density, molar volume, and fugacity.
"""

# Universal gas constant (R). units: m³-bar/(K-mol)
const R = 8.3144598e-5

"""
Data structure stating characteristics of a Peng-Robinson gas
# Attributes
- `gas::Symbol`: Peng-Robinson Gas
- `Tc::Float64`: Peng-Robinson Gas critical temperature (units: Kelvin)
- `Pc::Float64`: Peng-Robinson Gas critical pressure (units: bar)
- `ω::Float64`: Peng-Robinson Gas acentric factor (units: unitless)
"""
struct PengRobinsonGas
    gas::Symbol
    Tc::Float64
    Pc::Float64
    ω::Float64
end

"""
Evaluates the Peng-Robinson Equation of State to determine the the compressibility factor(z)
of a real gas using its critical temperature and pressure, and acentric factor.

# Arguments
- `gas:: PengRobinsonGas`:
- `κ::Float64`:
- `Tr::Float64`:
- `T::Float64`:
- `P::Float64`:

"""
a(gas::PengRobinsonGas) = (0.457235 * R ^ 2 * gas.Tc ^ 2) / gas.Pc
b(gas::PengRobinsonGas) = (0.0777961 * R * gas.Tc) / gas.Pc
κ(gas::PengRobinsonGas) = 0.37464 + (1.54226 * gas.ω) - (0.26992 * gas.ω ^ 2)
α(κ::Float64, Tr::Float64) = (1 + κ * (1 - √(Tr))) ^ 2
A(T::Float64, P::Float64, gas::PengRobinsonGas) = α(κ(gas), T / gas.Tc) * a(gas) * P / (R ^ 2 * T ^ 2)
B(T::Float64, P::Float64, gas::PengRobinsonGas) = b(gas) * P / (R * T)

"""
Calculates three outputs for compressibility factor using the polynomial form of
the Peng-Robinson Equation of State.
Filters for only real roots, and then returns root closest to unity.
"""
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

"""
Calculating for fugacity coefficient from an integration (bar)
"""
function calculate_ϕ(gas::PengRobinsonGas, T::Float64, P::Float64)
    z = compressibility_factor(gas, T, P)
    log_ϕ = z - 1.0 - log(z - B(T, P, gas)) +
            - A(T, P, gas) / (√8 * B(T, P, gas)) * log(
            (z + (1 + √2) * B(T, P, gas)) / (z + (1 - √(2)) * B(T, P, gas)))
    return exp(log_ϕ)
end

"""
Returns a dictionary of Peng-Robinson gas properties based on user input for T and P
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
Reads in file with gas parameters using dataframes and queries values for gas
"""
function PengRobinsonGas(gas::Symbol)
    df = CSV.read(PorousMaterials.PATH_TO_DATA * "PengRobinsonGasProps.csv")
    if ! (string(gas) in df[:gas])
        error(@sprintf("Gas %s properties not found in %sPengRobinsonGasProps.csv", gas, PATH_TO_DATA))
    end
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
