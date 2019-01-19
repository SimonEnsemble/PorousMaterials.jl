PorousMaterials relies on functions that fall outside of the sections discussed up to this point. They are useful, but operate on their own for the most part.

## PATH\_TO\_DATA Control

The PATH\_TO\_DATA is crucial for loading in data files. These functions allow the user to control this after they have done `using PorusMaterials`. The `set_tutorial_mode()` function has been discussed before, and it allows you to recreate our example and try PorouMaterials before loading your data. The other functions allow the user to reset the PATH\_TO\_DATA to the data folder in their current directory or to another directory on their machine if their files are stored in many places.

## Reading in Atomic Values

These functions are used to read in the `atomicmasses.csv`, `atom_properties.csv`, and `cpk_atom_colors.csv` files from the current `PATH_TO_DATA` directory. They contain information on the masses, radii, and cpk color scheme for atoms.

## Using .xyz files

These functions allow the user to load and save .xyz files describing where molecules appear in a given space. This can be used to save the location of molecules in the middle of a simulation or to visualize what is happening.

## PATH\_TO\_DATA Control
```@docs
    set_path_to_data
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
