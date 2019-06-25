# Calculates the properties of a real fluid, such as the compressibility factor, fugacity,
#   and molar volume.

# Universal fluid constant (R). units: m³-bar/(K-mol)
const R = 8.3144598e-5

# Data structure stating characteristics of a Peng-Robinson fluid
struct PengRobinsonFluid
    "Peng-Robinson Fluid species. e.g. :CO2"
    fluid::Symbol
    "Critical temperature (units: Kelvin)"
    Tc::Float64
    "Critical pressure (units: bar)"
    Pc::Float64
    "Acentric factor (units: unitless)"
    ω::Float64
end

# Data structure stating

# Parameters in the Peng-Robinson Equation of State
# T in Kelvin, P in bar
a(fluid::PengRobinsonFluid) = (0.457235 * R ^ 2 * fluid.Tc ^ 2) / fluid.Pc
b(fluid::PengRobinsonFluid) = (0.0777961 * R * fluid.Tc) / fluid.Pc
κ(fluid::PengRobinsonFluid) = 0.37464 + (1.54226 * fluid.ω) - (0.26992 * fluid.ω ^ 2)
α(κ::Float64, Tr::Float64) = (1 + κ * (1 - √Tr)) ^ 2
A(T::Float64, P::Float64, fluid::PengRobinsonFluid) = α(κ(fluid), T / fluid.Tc) * a(fluid) * P / (R ^ 2 * T ^ 2)
B(T::Float64, P::Float64, fluid::PengRobinsonFluid) = b(fluid) * P / (R * T)

# Calculates three outputs for compressibility factor using the polynomial form of
# the Peng-Robinson Equation of State. Filters for only real roots and returns the
# root closest to unity.
function compressibility_factor(fluid::PengRobinsonFluid, T::Float64, P::Float64)
    # construct cubic polynomial in z
    p = Poly([-(A(T, P, fluid) * B(T, P, fluid) - B(T, P, fluid) ^ 2 - B(T, P, fluid) ^ 3),
              A(T, P, fluid) - 2 * B(T, P, fluid) - 3 * B(T, P, fluid) ^ 2,
              -(1.0 - B(T, P, fluid)),
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
function calculate_ϕ(fluid::PengRobinsonFluid, T::Float64, P::Float64)
    z = compressibility_factor(fluid, T, P)
    log_ϕ = z - 1.0 - log(z - B(T, P, fluid)) +
            - A(T, P, fluid) / (√8 * B(T, P, fluid)) * log(
            (z + (1 + √2) * B(T, P, fluid)) / (z + (1 - √(2)) * B(T, P, fluid)))
    return exp(log_ϕ)
end

"""
    props = calculate_properties(fluid, T, P, verbose=true)

Use equation of state to calculate density, fugacity, and molar volume of a real fluid at a
given temperature and pressure.

# Arguments
- `fluid::PengRobinsonFluid`: Peng-Robinson fluid data structure
- `T::Float64`: Temperature (units: Kelvin)
- `P::Float64`: Pressure (units: bar)
- `verbose::Bool`: will print results if `true`

# Returns
- `prop_dict::Dict`: Dictionary of Peng-Robinson fluid properties
"""
function calculate_properties(fluid::PengRobinsonFluid, T::Float64, P::Float64; verbose::Bool=true)
    # Compressbility factor (unitless)
    z = compressibility_factor(fluid, T, P)
    # Density (mol/m^3)
    ρ = P / (z * R * T)
    # Molar volume (L/mol)
    Vm = 1.0 / ρ * 1000.0
    # Fugacity (bar)
    ϕ = calculate_ϕ(fluid, T, P)
    f = ϕ * P
    # Prints a dictionary holding values for compressibility factor, molar volume, density, and fugacity.
    prop_dict = Dict("compressibility factor" => z, "molar volume (L/mol)"=> Vm ,
                     "density (mol/m³)" => ρ, "fugacity (bar)" => f,
                     "fugacity coefficient" => ϕ)
    if verbose
        @printf("%s properties at T = %f K, P = %f bar:\n", fluid.fluid, T, P)
        for (property, value) in prop_dict
            println("\t" * property * ": ", value)
        end
    end
    return prop_dict
end


"""
    fluid = PengRobinsonFluid(fluid)

Reads in critical temperature, critical pressure, and acentric factor of the `fluid::Symbol`
from the properties .csv file `joinpath(PorousMaterials.PATH_TO_DATA, "PengRobinsonFluidProps.csv")`
and returns a complete `PengRobinsonFluid` data structure.
**NOTE: Do not delete the last three comment lines in PengRobinsonFluidProps.csv

# Arguments
- `fluid::Symbol`: The fluid molecule you wish to construct a PengRobinsonFluid struct for

# Returns
- `PengRobinsonFluid::struct`: Data structure containing Peng-Robinson fluid parameters.
"""
function PengRobinsonFluid(fluid::Symbol)
    df = CSV.read(joinpath(PATH_TO_DATA, "PengRobinsonFluidProps.csv"); footerskip=3)
    if ! (string(fluid) in df[:fluid])
        error(@sprintf("fluid %s properties not found in %sPengRobinsonFluidProps.csv", fluid, PATH_TO_DATA))
    end
    Tc = df[df[:fluid].== string(fluid), Symbol("Tc(K)")][1]
    Pc = df[df[:fluid].== string(fluid), Symbol("Pc(bar)")][1]
    ω = df[df[:fluid].== string(fluid), Symbol("acentric_factor")][1]
    return PengRobinsonFluid(fluid, Tc, Pc, ω)
end

# Prints resulting values for Peng-Robinson fluid properties
function Base.show(io::IO, fluid::PengRobinsonFluid)
    println(io, "fluid species: ", fluid.fluid)
    println(io, "\tCritical temperature (K): ", fluid.Tc)
    println(io, "\tCritical pressure (bar): ", fluid.Pc)
    println(io, "\tAcenteric factor: ", fluid.ω)
end
