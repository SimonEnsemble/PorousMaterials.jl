# Equation of State

`PorousMaterials.jl` provides Peng-Robinson and van der Waals equation of state calculations to find the properties of a real fluid.

## Peng-Robinson equation of state
The Peng-Robinson equation of state can be written as:

![PREOS](https://latex.codecogs.com/gif.latex?P%20%3D%5Cfrac%7BRT%7D%7BV%7Bm%7D-b%7D-%5Cfrac%7Ba%5Calpha%7D%7BV%7B%7Bm%7D%7D%5E%7B2%7D&plus;2bV%7Bm%7D-b%5E%7B2%7D%7D)

where $V_{m}$ is the molar volume, $R$ is the gas constant, and $T$ is temperature.

Variables $a$ and $b$ can be calculated using the critical temperature, $T_{c}$ and pressure, $P_{c}$, of the fluid:

![PREOS_a](https://latex.codecogs.com/gif.latex?a%20%5Capprox%20%5Cfrac%7B0.457235%26space%3BR%5E%7B2%7DT%20_%7Bc%7D%5E%7B2%7D%7D%7BP_%7Bc%7D%7D)

![PREOS_b](https://latex.codecogs.com/gif.latex?b%26space%3B%5Capprox%26space%3B%5Cfrac%7B0.07780%26space%3BR%26space%3B%20T_%7Bc%7D%7D%7BP_%7Bc%7D%7D)

and $\alpha$ can be calculated using acentric factor $\omega$ and critical temperature:

![PREOS_alpha](https://latex.codecogs.com/gif.latex?%5Calpha%26space%3B%3D%26space%3B%281%26plus%3B%5Ckappa%281-T_%7Br%7D%5E%7B%5Cfrac%20%7B1%7D%7B2%7D%7D%29%29%5E%7B2%7D)

where

![PREOS_kappa](https://latex.codecogs.com/gif.latex?%5Ckappa%20%5Capprox%200.37464&plus;1.54226%5Comega-0.26992%5Comega%5E%7B2%7D)

and

![PREOS_Tr](https://latex.codecogs.com/gif.latex?T_%7Br%7D%26space%3B%3D%26space%3B%5Cfrac%7BT%7D%7BT_%7Bc%7D%7D)

## Van der Waals equation of state
The van der Waals equation can be written as:
![VDWEOS](https://latex.codecogs.com/gif.latex?%28P%26space%3B%26plus%3B%26space%3B%5Cfrac%7Ba%7D%7BV%7B_%7Bm%7D%7D%5E%7B2%7D%7D%29%28V_%7Bm%7D%20%26space%3B-%26space%3Bb%29%26space%3B%3D%26space%3BRT)

where $a$ and $b$ are the van der Waals constants. $a$ and $b$ can be calculated from fluid critical propertis, but `PorousMaterials.jl` reads them in as experimentally determined values.

## reading in fluid characteristics

`PengRobinsonFluid` and `VdWFluid` are structs defining the characteristics of a fluid of interest, depending on which equation of state is used.

For Peng-Robinson fluids, `PorousMaterials.jl` reads in the critical temperature, critical pressure, and acentric factor of `fluid::Symbol` from the properties .csv file `PorousMaterials.PATH_TO_DATA, "PengRobinson_fluid_props.csv")`.

For van der Waals fluids, van der Waals constants of `fluid::Symbol` are read in from the properties .csv file `joinpath(PorousMaterials.PATH_TO_DATA, "VdW_fluid_props.csv")`.

*** NOTE: DO NOT DELETE LAST THREE COMMENT LINES IN `PengRobinson_fluid_props.csv` AND `VdW_fluid_props.csv`

The characteristics can be read as:

For Peng-Robinson fluids

```julia
fluid = PengRobinsonFluid(:Xe)       # Input fluid as a symbol. The fluids reader stores the information in fluid as a struct
fluid.fluid                          # The name of the fluid
fluid.Pc                             # The critical pressure of the fluid
fluid.Tc                             # The critical temperature of the fluid
fluid.ω                              # The acentric factor of the fluid
```

For van der Waals fluids

```julia
fluid = VdWFluid(:Xe)                # Input fluid as a symbol. The fluids reader stores the information in fluid as a struct
fluid.fluid                          # The name of the fluid
fluid.a                              # The van der Waals constant a of the fluid
fluid.b                              # The van der Waals constant b of the fluid
```

## calculating density, fugacity, and molar volume
Using a given temperature and pressure, `PorousMaterials.jl` the equation of state can be used to calculate the dnesity, fugacity, and molar volume of a real fluid, stored as a dictionary.
```julia
T = 298.0 # K                        # The temperature in Kelvin of interest type Float64.
P = 1.0 # bar                        # The pressure in bar of interest type Float64.
props = calculate_properties(fluid, T, P, verbose=true) # verbose::Bool will print results if `true`
```

The output is a dictionary containing the following keys:
```julia
props["compressibility factor"]      # the compressibility factor
props["density [mol/$m^{3}$]"]       # fluid density in mol/$m^{3}$
props["fugacity [bar]"]              # the fugacity in bar
props["fugacity coefficient"]        # the fugacity coefficient
props["molar volume [L/mol]"]        # the molar volume in L/mol
```

# details
```@docs
    PengRobinsonFluid
    VdWFluid
    calculate_properties
```

