PorousMaterials relies on functions that fall outside of the sections discussed up to this point. They are useful, but operate on their own for the most part.

## Control over directories for input and output files
Call `print_file_paths()` to show where input/output files are read from/written to.

To change the directory where certain files are read from/output to, for example:

```
@eval PorousMaterials PATH_TO_CRYSTALS = joinpath(pwd(), "my_xtals")
@eval PorousMaterials PATH_TO_FORCEFIELDS = joinpath(pwd(), "my_ffs")
@eval PorousMaterials PATH_TO_GRIDS = joinpath(pwd(), "my_grids")
@eval PorousMaterials PATH_TO_MOLECULES = joinpath(pwd(), "my_molecules")
@eval PorousMaterials PATH_TO_DATA = pwd()
```

To use absolute file paths when reading in e.g. crystals via the `Framework` constructor, you can simply:
```
@eval PorousMaterials PATH_TO_CRYSTALS = ""
```

Call `set_default_file_paths()` to set input/output file paths back to default.


The PATH\_TO\_DATA is crucial for loading in data files. These functions allow the user to control this after they have done `using PorousMaterials`. The `set_tutorial_mode()` function has been discussed before, and it allows you to recreate our example and try PorousMaterials before loading your data. The other functions allow the user to reset the PATH\_TO\_DATA to the data folder in their current directory or to another directory on their machine if their files are stored in many places.

## Reading in Atomic Values

These functions are used to read in the `atomicmasses.csv`, `atom_properties.csv`, and `cpk_atom_colors.csv` files from the current `PATH_TO_DATA` directory. They contain information on the masses, radii, and cpk color scheme for atoms.

## Using .xyz files

These functions allow the user to load and save .xyz files describing where molecules appear in a given space. This can be used to save the location of molecules in the middle of a simulation or to visualize what is happening.

## Fitting data to adsorption models

PorousMaterials allows for a `DataFrame` to be read in and fitted to a single-site Langmuir model or to a Henry's law model.
```
using PorousMaterials
using DataFrames

adsorption_data = DataFrame(Dict("Pressure (bar)" => [0.008045, 0.042218, 0.078772, 0.108018, 0.156741, 0.312670, 0.414986, 0.517303, 0.619628, 
                                                      0.719519, 0.821872, 0.863296, 0.912055, 0.960878, 0.982918, 0.990353, 0.995361, 0.998043,
                                                      1.000610, 1.005600, 1.005720],
                                 "Adsorption (mmol/g)" => [4.062603, 4.462196, 4.560714, 4.659598, 4.707321, 4.950402, 5.045670, 5.140938,
                                                           5.286339, 5.431875, 5.727768, 5.826027, 6.074420, 6.673929, 7.324955, 8.026830,
                                                           8.778973, 10.133170, 10.835313, 11.487143, 12.189375]))
# We can fit the adsorption data to a Langmuir isotherm with the following call. Note that we need to enter the column names for the pressure and adsorption.
results = fit_adsorption_isotherm(adsorption_data, Symbol("Pressure (bar)"), Symbol("Adsorption (mmol/g)"), :langmuir)
M, K = results["M"], results["K"]

# We can also use the function to get the Henry Coefficient
results = fit_adsorption_isotherm(adsorption_data, Symbol("Pressure (bar)"), Symbol("Adsorption (mmol/g)"), :henry)
H = results["H"]
```

## PATH\_TO\_DATA Control
```@docs
    print_file_paths
    set_default_file_paths
    set_tutorial_mode
```

## Reading in Atomic Values
```@docs
    read_atomic_radii
    read_atomic_masses
    read_cpk_colors
```

## Using .xyz files
```@docs
    read_xyz
    write_xyz
```

## Generic Rotations
```@docs
    rotation_matrix
```

## Fitting data
```@docs
    fit_adsorption_isotherm
```
