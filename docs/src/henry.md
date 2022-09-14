# Henry Coefficient

`PorousMaterials.jl` allows for Henry coefficient calculations using Widom insertions.

## Preparing the Henry coefficient simulation

The simulation requires the following `PorousMaterials.jl` objects:
* `Crystal` structure
* `Molecule` adsorbate
* `LJForceField` forcefield

In addition the the list above, one has to specify the temperature (in K) and the number of Widom insertions per unit volume (in Angstrom).

```julia
xtal = Crystal("SBMOF-1.cif")           # The crystal structure we are interested in
strip_numbers_from_atom_labels!(xtal)   # We have to make sure the atom species have no numbers appended to them
methane = Molecule("CH4")               # Here we choose to use methane as the adsorbate
ljff = LJForceField("UFF")              # We will use the Universal Force Field (UFF) to calculate the interatomic interactions
temp = 298.0                            # Standard temperature (K)
widom_insertions = 2000                 # Number of insertions per unit volume

results = henry_coefficient(xtal, methane, temp, ljff, insertions_per_volume=widom_insertions)
```

The results are also saved to `rc[:paths][:simulations]` as a `.jld2` file that can be read using the `JLD2` package.

The output (and saved file) is a dictionary:

```julia
results
# output
Dict{String, Any} with 16 entries:
  "⟨U⟩ (K)"                              => -3623.51
  "err Qst (kJ/mol)"                     => 0.0917643
  "⟨U, vdw⟩ (kJ/mol)"                    => -30.1275
  "⟨U, es⟩ (kJ/mol)"                     => 0.0
  "elapsed time (min)"                   => 0.0978277
  "Qst (kJ/mol)"                         => 32.6052
  "err henry coefficient [mmol/(g-bar)]" => 6.12256
  "xtal"                                 => "SBMOF-1.cif"
  "henry coefficient [mmol/(g-bar)]"     => 46.1068
  "err ⟨U, es⟩ (kJ/mol)"                 => 0.0
  "⟨U, vdw⟩ (K)"                         => -3623.51
  "err ⟨U, vdw⟩ (kJ/mol)"                => 0.0917643
  "⟨U, es⟩ (K)"                          => 0.0
  "⟨U⟩ (kJ/mol)"                         => -30.1275
  "henry coefficient [mol/(kg-Pa)]"      => 0.000461068
  "henry coefficient [mol/(m³-bar)]"     => 72406.2
```

## locating the saved results

The name of the result filenames follow a convention outlined in `henry_result_savename`.
```julia
using JLD2
# determine the canonical filename for the simulation
result_filename = henry_result_savename(xtal, methane, temp, ljff, widom_insertions)
# load the results dictionary
@load joinpath(rc[:paths][:simulations], result_filename) results
```

# detailed docs

```@docs
    henry_coefficient
    henry_result_savename
```
