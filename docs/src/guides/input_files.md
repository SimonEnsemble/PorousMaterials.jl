## Input files to describe crystals, molecules, and forcefields

All input files are stored in the path `PorousMaterials.PATH_TO_DATA` (type into Julia). By default, this path is set to be in the present working directory (type `pwd()` into Julia) in a folder `data/`. Go inside `PorousMaterials.jl/test/data` to see example input files for each case below.

There will be example code snippets through the documentation showing how to load in various files. To get a feel for this we have included a Tutorial Mode in `PorousMaterials.jl` that sets the `PorousMaterial.PATH_TO_DATA` to the data folder in our testing directory. To follow along with the examples without downloading your own data simply do the following:

```

```

#### Molecules/Adsorbates

Molecule input files are stored in `data/molecules`. Each molecule possesses its own directory and contains two files: `point_charges.csv` and `lennard_jones_spheres.csv`, comma-separated-value files describing the point charges and Lennard Jones spheres, respectively, comprising the molecule. Only rigid molecules are currently supported. Units of length are in Angstrom; units of charges are electrons.

#### Force field parameters

Lennard-Jones forcefield parameters are stored in comma-separated-value format in `data/forcefields/`.

Interaction of an adsorbate with the framework is modeled as pair-wise additive and with Lennard-Jones potentials of the form:

`V(r) = 4 * epsilon * [ x ^ 12 - x ^ 6 ]`, where `x = sigma / r`

The Lennard Jones force field input files, e.g. `UFF.csv` contain a list of pure (i.e. X-X, where X is an atom) sigmas and epsilons in units Angstrom and Kelvin, respectively. Note that, e.g., in the UFF paper, the Lennard Jones potential is written in a different form and thus parameters need to be converted to correspond to the functional form used in `PorousMaterials.jl`.

#### Atomic masses

Add fancy pseudo-atoms to `data/atomic_masses.csv`.

#### Peng-Robinson gas parameters

Critical temperatures and pressures and acentric factors are stored in `data/PengRobinsonGasProps.csv`.
