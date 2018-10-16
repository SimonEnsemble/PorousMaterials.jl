
PorousMaterials relies on functions that fall outside of the sections discussed up to this point. They are useful, but operate on their own for the most part.


<a id='PATH_TO_DATA-Control-1'></a>

## PATH_TO_DATA Control


The PATH_TO_DATA is crucial for loading in data files. These functions allow the user to control this after they have done `using PorusMaterials`. The `set_tutorial_mode()` function has been discussed before, and it allows you to recreate our example and try PorouMaterials before loading your data. The other functions allow the user to reset the PATH_TO_DATA to the data folder in their current directory or to another directory on their machine if their files are stored in many places.


<a id='Reading-in-Atomic-Values-1'></a>

## Reading in Atomic Values


These functions are used to read in the `atomicmasses.csv`, `atom_properties.csv`, and `cpk_atom_colors.csv` files from the current `PATH_TO_DATA` directory. They contain information on the masses, radii, and cpk color scheme for atoms.


<a id='Using-.xyz-files-1'></a>

## Using .xyz files


These functions allow the user to load and save .xyz files describing where molecules appear in a given space. This can be used to save the location of molecules in the middle of a simulation or to visualize what is happening.

