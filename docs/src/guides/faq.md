## FAQ

**How do I type out the math symbols? e.g. `box.α`?**

Julia supports [unicode input](https://docs.julialang.org/en/release-0.4/manual/unicode-input/)! Type `box.\alpha`, then hit tab. Voilà. There is a vim extension for Julia [here](https://github.com/JuliaEditorSupport/julia-vim).


**How do I run as a script in the command line?**

It is instructive to first run an example in the Julia REPL so you can print out and interact with attributes of your `forcefield`, `framework`, and `molecule` to ensure they are correct. If you want to then run the Julia code in the command line, simply put the commands in a text file with a `.jl` extension and run in terminal as `julia my_script.jl`. For parallelization in `adsorption_isotherm` and `henry_coefficient`, call e.g. 4 cores with `julia -p 4 my_script.jl`.

**Can I use `PorousMaterials.jl` in Jupyter Notebook/ Jupyter Lab?**

Yes! See [here](https://github.com/JuliaLang/IJulia.jl).

**How can I convert my `.cif` into P1 symmetry for `PorousMaterials.jl`?**

`PorousMaterials.jl` will automatically do this for you! It looks for the
`_symmetry_equiv_pos_as_xyz` tag in the `.cif` file and uses those symmetry operations to replicate the structure in a lower symmetry into P1 symmetry.

It is important to note that `PorousMaterials.jl` will read in the space group
name, but it does **NOT** use this for converting your structure to P1.
