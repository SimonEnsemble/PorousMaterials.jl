using Polynomials
using DataFrames
using CSV
using Roots

"""
Calculates the properties of a real gas, such as the compressibility factor, fugacity,
and molar volume.
"""

# Universal gas constant (R). units: m³-bar/(K-mol)
const R = 8.3144598e-5

# Data structure stating characteristics of a Peng-Robinson gas
struct PengRobinsonFluid
  "Peng-Robinson Gas species. e.g. :CO2"
  gas::Symbol
  "Critical temperature (units: Kelvin)"
  Tc::Float64
  "Critical pressure (units: bar)"
  Pc::Float64
  "Acentric factor (units: unitless)"
  ω::Float64
end

#Data structure stating characteristics of a Van der Waals gas
struct VDWFluid
    "VDW constant a (units: m⁶bar/mol²)"
    a::Float64
    "VDW constant b (units: m³/mol)"
    b::Float64
    "Van der Waals Gas species e.g. :H2"
    gas::Symbol
end

# Parameters in the Peng-Robinson Equation of State
# T in Kelvin, P in bar
a(gas::PengRobinsonFluid) = (0.457235 * R ^ 2 * gas.Tc ^ 2) / gas.Pc
b(gas::PengRobinsonFluid) = (0.0777961 * R * gas.Tc) / gas.Pc
κ(gas::PengRobinsonFluid) = 0.37464 + (1.54226 * gas.ω) - (0.26992 * gas.ω ^ 2)
α(κ::Float64, Tr::Float64) = (1 + κ * (1 - √(Tr))) ^ 2
A(T::Float64, P::Float64, gas::PengRobinsonFluid) = α(κ(gas), T / gas.Tc) * a(gas) * P / (R ^ 2 * T ^ 2)
B(T::Float64, P::Float64, gas::PengRobinsonFluid) = b(gas) * P / (R * T)

# Calculates three outputs for compressibility factor using the polynomial form of
# the Peng-Robinson Equation of State. Filters for only real roots and returns the
# root closest to unity.
function compressibility_factor(gas::PengRobinsonFluid, T::Float64, P::Float64)

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
function calculate_ϕ(gas::PengRobinsonFluid, T::Float64, P::Float64)
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
- `gas::Union[PengRobinsonFluid, VDWFluid]`:  Gas structure
- `T::Float64`: Temperature (units: Kelvin)
- `P::Float64`: Pressure (units: bar)
- `verbose::Bool`: print results

# Returns
- `prop_dict::Dict`: Dictionary of Peng-Robinson gas properties
"""
function calculate_properties(gas::PengRobinsonFluid, T::Float64, P::Float64; verbose::Bool=true)
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

function calculate_properties(gas::VDWFluid, T::Float64, P::Float64; verbose::Bool=true)

    A = -P
    B = (P * gas.b + R * T)
    C = -gas.a
    D = gas.a * gas.b

    #Creates a polynomial for the vdw cubic function
    pol = Poly([A, B, C, D])
    #finds roots of that polynomial
    polroots = roots(pol)
    #assigns rho to be the real root(s) and then makes it real to get rid of the 0im
    rho = real.(polroots[isreal.(polroots)])
    #disregards all roots except the lowest one, as the lowest real root
    #is the density corresponding to the gas phase
    rho = rho[argmin(rho)]
    #specifies that molar volume is the reciprocal of the density
    # In units of [L/mol]
    vm = (1 ./ rho) * 1000
    #specifies the compressibility factor
    z = (P * (1 ./ rho))./ (R * T)

    #Finds fugacity using the derivation from the Van der Waals
    fug = P .* exp. (- log. (((1 ./ rho) - gas.b) * P./(R * T))+(gas.b ./ ((1 ./ rho)-gas.b) - 2*gas.a*rho/(R*T)))
    #defines the fugacity coefficient as fugacity over pressure
    ϕ = fug ./ P

    prop_dict = Dict("Density (mol/m³)" => rho, "Fugacity (bar)" => fug,
        "Molar Volume (L/mol)" => vm, "Fugacity Coefficient" => ϕ,
        "Compressibility Factor" => z )

    if verbose
        @printf("%s properties at T = %f K, P = %f bar:\n", gas.gas, T, P)
        for (property, value) in prop_dict
            println(property * ": ", value)
        end
    end
    return prop_dict
end

"""
    gas = VDWFluid(gas)

Reads in vdw constants a and b of the `gas::Symbol`
from the properties .csv file `PorousMaterials.PATH_TO_DATA * "vdw_constants.csv"`
and returns a complete `VDWFluid` data structure.

# Returns
- `VDWFluid::struct`: Data structure containing Van der Waals gas parameters.
"""
function VDWFluid(gas::Symbol)
    vdwfile = CSV.read(joinpath(PATH_TO_DATA, "vdw_constants.csv"))
    if ! (string(gas) in vdwfile[:molecule])
          error(@sprintf("Gas %s properties not found in %sVDW_Constants.csv", gas, PATH_TO_DATA))
    end
    gas = string(gas)
    A = vdwfile[vdwfile[:molecule].== gas, Symbol("a(m6bar/mol2)")]
    B = vdwfile[vdwfile[:molecule].== gas, Symbol("b(m3/mol)")]
    return VDWFluid(A[1], B[1], Symbol(gas))
end

"""
    gas = PengRobinsonFluid(gas)

Reads in critical temperature, critical pressure, and acentric factor of the `gas::Symbol`
from the properties .csv file `PorousMaterials.PATH_TO_DATA * "PengRobinsonGasProps.csv"`
and returns a complete `PengRobinsonFluid` data structure.
**NOTE: Do not delete the last three comment lines in PengRobinsonGasProps.csv

# Returns
- `PengRobinsonFluid::struct`: Data structure containing Peng-Robinson gas parameters.
"""
function PengRobinsonFluid(gas::Symbol)
    df = CSV.read(joinpath(PATH_TO_DATA, "PengRobinsonGasProps.csv"); footerskip=3)
    if ! (string(gas) in df[:gas])
        error(@sprintf("Gas %s properties not found in %sPengRobinsonGasProps.csv", gas, PATH_TO_DATA))
    end
    Tc = df[df[:gas].== string(gas), Symbol("Tc(K)")][1]
    Pc = df[df[:gas].== string(gas), Symbol("Pc(bar)")][1]
    ω = df[df[:gas].== string(gas), Symbol("acentric_factor")][1]
    return PengRobinsonFluid(gas, Tc, Pc, ω)
end

# Prints resulting values for Peng-Robinson gas properties
function Base.show(io::IO, gas::PengRobinsonFluid)
    println(io, "Gas species: ", gas.gas)
    println(io, "Critical temperature (K): ", gas.Tc)
    println(io, "Critical pressure (bar): ", gas.Pc)
    println(io, "Acenteric factor: ", gas.ω)
end

# Prints resulting values for Van der Waals gas properties
function Base.show(io::IO, gas::VDWFluid)
    println(io, "\tGas species: ", gas.gas)
    println(io, "\tVan der Waals constant a (m⁶bar/mol²): ", gas.a)
    println(io, "\tVan der Waals constant b (m³/mol): ", gas.b)
end
