# Calculates the properties of a real gas, such as the compressibility factor, fugacity,
#   and molar volume.

# Universal gas constant (R). units: m³-bar/(K-mol)
const R = 8.3144598e-5

# Data structure stating characteristics of a Peng-Robinson gas
struct PengRobinsonGas
    "Peng-Robinson Gas species. e.g. :CO2"
    gas::Symbol
    "Critical temperature (units: Kelvin)"
    Tc::Float64
    "Critical pressure (units: bar)"
    Pc::Float64
    "Acentric factor (units: unitless)"
    ω::Float64
end

# Parameters in the Peng-Robinson Equation of State
# T in Kelvin, P in bar
a(gas::PengRobinsonGas) = (0.457235 * R ^ 2 * gas.Tc ^ 2) / gas.Pc
b(gas::PengRobinsonGas) = (0.0777961 * R * gas.Tc) / gas.Pc
κ(gas::PengRobinsonGas) = 0.37464 + (1.54226 * gas.ω) - (0.26992 * gas.ω ^ 2)
α(κ::Float64, Tr::Float64) = (1 + κ * (1 - √Tr)) ^ 2
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
    id_closest_to_unity = argmin(abs.(z_factor .- 1.0))
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
    props = calculate_properties(gas, T, P, verbose=true)

Use equation of state to calculate density, fugacity, and molar volume of a real gas at a
given temperature and pressure.

# Arguments
- `gas::PengRobinsonGas`: Peng-Robinson gas data structure
- `T::Float64`: Temperature (units: Kelvin)
- `P::Float64`: Pressure (units: bar)
- `verbose::Bool`: will print results if `true`

# Returns
- `prop_dict::Dict`: Dictionary of Peng-Robinson gas properties
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
            println("\t" * property * ": ", value)
        end
    end
    return prop_dict
end


"""
    gas = PengRobinsonGas(gas)

Reads in critical temperature, critical pressure, and acentric factor of the `gas::Symbol`
from the properties .csv file `joinpath(PorousMaterials.PATH_TO_DATA, "PengRobinsonGasProps.csv")`
and returns a complete `PengRobinsonGas` data structure.
**NOTE: Do not delete the last three comment lines in PengRobinsonGasProps.csv

# Arguments
- `gas::Symbol`: The gas molecule you wish to construct a PengRobinsonGas struct for

# Returns
- `PengRobinsonGas::struct`: Data structure containing Peng-Robinson gas parameters.
"""
function PengRobinsonGas(gas::Symbol)
    df = CSV.read(joinpath(PATH_TO_DATA, "PengRobinsonGasProps.csv"); footerskip=3)
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
    println(io, "\tCritical temperature (K): ", gas.Tc)
    println(io, "\tCritical pressure (bar): ", gas.Pc)
    println(io, "\tAcenteric factor: ", gas.ω)
end
