![PorousMaterials.jl](assets/PMlogo.png)
A pure-[Julia](https://julialang.org/) package for classical molecular modeling of adsorption in porous crystals such as metal-organic frameworks (MOFs).

🔨 Compute the potential energy of a molecule at particular position and orientation inside of a porous crystal

🔨 Write a potential energy grid of a molecule inside a porous material to visualize binding sites

🔨 Compute the Henry coefficient of a gas in a porous crystal

🔨 Run grand-canonical Monte Carlo simulations of gas adsorption in a porous crystal

Designed for high-throughput computations to minimize input files and querying results from output files. User-friendly. Instructive error messages thrown when they should be. Well-documented (eventually). Easy to install (eventually).

*In development, please contribute, post issues 🐛, and improve!*

## Installation

1. Download and install the [Julia programming language](https://julialang.org/), v1.0.

2. In Julia, the following code will clone this repository and install Julia package dependencies in `REQUIRE`:

```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/SimonEnsemble/PorousMaterials.jl"))
Pkg.test("PorousMaterials") # run tests
```

3. In Julia, load all functions in `PorousMaterials.jl` into the namespace:

```julia
using PorousMaterials # that's it
```

Note: This package is in development. After stabilized and fully documented, installation will be as easy as `Pkg.add("PorousMaterials")`.

## Tests
Run the unit-ish tests in the script `tests/runtests.jl` manually or type `Pkg.test("PorousMaterials")` into Julia.

Direct tests for Henry coefficients and grand-canonical Monte Carlo simulations take much longer and are found in `tests/henry_test.jl` and `tests/gcmc_test.jl`.

## FAQ

**How do I type out the math symbols? e.g. `box.α`?**

Julia supports [unicode input](https://docs.julialang.org/en/release-0.4/manual/unicode-input/)! Type `box.\alpha`, then hit tab. Voilà. There is a vim extension for Julia [here](https://github.com/JuliaEditorSupport/julia-vim).


**How do I run as a script in the command line?**

It is instructive to first run an example in the Julia REPL so you can print out and interact with attributes of your `forcefield`, `framework`, and `molecule` to ensure they are correct. If you want to then run the Julia code in the command line, simply put the commands in a text file with a `.jl` extension and run in terminal as `julia my_script.jl`. For parallelization in `adsorption_isotherm` and `henry_coefficient`, call e.g. 4 cores with `julia -p 4 my_script.jl`.

**Can I use `PorousMaterials.jl` in Jupyter Notebook/ Jupyter Lab?**

Yes! See [here](https://github.com/JuliaLang/IJulia.jl).

**How can I convert my `.cif` into P1 symmetry for `PorousMaterials.jl`?**

We hope someone will contribute this feature to `PorousMaterials.jl` eventually. For now, you can use [OpenBabel](http://openbabel.org/wiki/Main_Page):

```
obabel -icif non-P1.cif -ocif -O P1.cif --fillUC strict
```
