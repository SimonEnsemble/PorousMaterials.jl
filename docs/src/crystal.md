# Crystals

`PorousMaterials.jl` maintains a data structure `Crystal` that stores information about a crystal structure file.

## reading in a crystal structure file

Currently, the crystal structure file reader accepts `.cif` and `.cssr` file formats. `PorousMaterials.jl` looks for the crystal structure files in `PorousMaterials.PATH_TO_CRYSTALS` which is by default `./data/crystals`. By typing `@eval PorousMaterials PATH_TO_CRYSTALS = "my_crystal_dir"`, `PorousMaterials.jl` now looks for the crystal structure file in `my_crystal_dir`.
The files can be read as:

```julia
xtal = Crystal("IRMOF-1.cif")       # The crystal reader stores the information in xtal
xtal.name                           # The name of the crystal structure file
xtal.box                            # The unit cell information
xtal.atoms                          # The atom coordinates (in fractional space) and the atom identities
xtal.charges                        # The charge magnitude and coordinates (in fractional space)
xtal.bonds                          # Bonding information in the structure. By default this is an empty graph,
                                    #  but use `read_bonds_from_file=true` argument in `Crystal` to read from crystal structure file
xtal.symmetry                       # Symmetry information of the crystal. By default converts the symmetry to P1 symmetry.
                                    #  Use `convert_to_p1=false` argument in `Crystal` to keep original symmetry
```

### fixing atom species

Often, the atoms species are appended by numbers. This messes with the internal workings of `PorousMaterials.jl`.
To circumvent this problem, the function `strip_numbers_from_atom_labels!(xtal)` removes the appending numbers.
It is important to use this function prior to GCMC or Henry coefficient calculations.
```julia
xtal.atoms.species              # [:C1, :C2, :O1, ...]
strip_numbers_from_atom_labels!(xtal)
xtal.atoms.species              # [:C, :C, :O, ...]
```

### converting the coordinates to cartesian space

The coordinates of the crystals are stored in fractional coordinates. If one needs to analyze the cartesian coordinates of the crystal,
that can be done by using the unit cell information of the crystal.
```julia
xtal.atoms.coords.xf                                    # array of fractional coordinates
cart_coords = xtal.box.f_to_c * xtal.atoms.coords.xf    # array of cartesian coordinates
```

## creating a super cell

For many simulations, one needs to replicate the unit cell multiple times to create a bigger super cell.
```
super_xtal = replicate(xtal, (2,2,2))       # Replicates the original unit cell once in each dimension
```

## finding other properties
```julia
rho = crystal_density(xtal)         # Crystal density of the crystal in kg/m^2
mw = molecular_weight(xtal)         # The molecular weight of the unit cell in amu
formula = chemical_formula(xtal)    # The irreducible chemical formula of the crystal
```

## assigning new charges
If the crystal structure file does not contains partial charges, we provide methods to assign new charges to the crystal
```julia
species_to_charges = Dict(:Ca => 2.0, :C => 1.0, :H => -1.0)                # This method assigns a static charge to atom species
charged_xtal = assign_charges(xtal, species_to_charge, 1e-5)                # This function creates a new charged `Crystal` object.
                                                                            #   The function checks for charge neutrality with a tolerance of 1e-5
new_charges = Charges([2.0, 1.0, -1.0, -1.0, ...], xtal.atoms.coords)
other_charged_xtal = Crystal(xtal.name, xtal.box, xtal.atoms,               # Here we create a new `Charges` object using an array of new charges.
                             new_charges, xtal.bonds, xtal.symmetry)        #   The number of charges in the array has to be equal to the number of atoms
                                                                            #   and finally a new `Crystal` object is manually created
```

## writing crystal files
We provide methods to write both `.xyz` and `.cif` files
```julia
write_cif(xtal, "my_new_cif_file.cif")      # Stored in the current directory
write_xyz(xtal, "my_new_xyz_file.xyz")      # stored in the current directory
```


# detailed docs

```@docs
    Crystal
    SymmetryInfo
    replicate
    molecular_weight
    crystal_density
    chemical_formula
    assign_charges
    write_cif
    write_xyz
```
