```@meta
DocTestSetup = quote
  using PorousMaterials
end
```

# Crystals

`PorousMaterials.jl` maintains a data structure `Crystal` that stores information about a crystal structure file.

## reading in a crystal structure file

Currently, the crystal structure file reader accepts `.cif` and `.cssr` file formats. `PorousMaterials.jl` looks for the crystal structure files in `rc[:paths][:crystals]` which is by default `./data/crystals`. By typing `rc[:paths][:crystals] = "my_crystal_dir"`, `PorousMaterials.jl` now looks for the crystal structure file in `my_crystal_dir`.
The files can be read as:

```jldoctest crystal; output=false
xtal = Crystal("SBMOF-1.cif")       # The crystal reader stores the information in xtal
xtal.name                           # The name of the crystal structure file
xtal.box                            # The unit cell information
xtal.atoms                          # The atom coordinates (in fractional space) and the atom identities
xtal.charges                        # The charge magnitude and coordinates (in fractional space)
xtal.bonds                          # Bonding information in the structure. By default this is an empty graph,
                                    #  but use `read_bonds_from_file=true` argument in `Crystal` to read from crystal structure file
xtal.symmetry                       # Symmetry information of the crystal. By default converts the symmetry to P1 symmetry.
                                    #  Use `convert_to_p1=false` argument in `Crystal` to keep original symmetry
# output
Xtals.SymmetryInfo(["x"; "y"; "z"], "P1", true)
```

## fixing atom species

Often, the atoms species are appended by numbers. This messes with the internal workings of `PorousMaterials.jl`.
To circumvent this problem, the function `strip_numbers_from_atom_labels!(xtal)` removes the appending numbers.
It is important to use this function prior to GCMC or Henry coefficient calculations.

```jldoctest crystal; output=false
strip_numbers_from_atom_labels!(xtal)
# output

```

## converting the coordinates to cartesian space

The coordinates of the crystals are stored in fractional coordinates. If one needs to analyze the cartesian coordinates of the crystal,
that can be done by using the unit cell information of the crystal.

```jldoctest crystal
xtal.atoms.coords.xf                                    # array of fractional coordinates
cart_coords = xtal.box.f_to_c * xtal.atoms.coords.xf    # array of cartesian coordinates
# output
3×120 Matrix{Float64}:
 4.59487  -0.95272   2.68943   8.23701  …  8.8164    0.839249  -1.53211
 1.43955   4.2229    4.12715   1.3438      1.35443   1.42892    4.21227
 5.89964   5.35922  16.6181   17.1585      6.27862  17.5375    16.2391
```

## creating a super cell

For many simulations, one needs to replicate the unit cell multiple times to create a bigger super cell.

```jldoctest crystal
super_xtal = replicate(xtal, (2,2,2))       # Replicates the original unit cell once in each dimension
# output
Name: SBMOF-1.cif
Bravais unit cell of a crystal.
	Unit cell angles α = 90.000000 deg. β = 100.897000 deg. γ = 90.000000 deg.
	Unit cell dimensions a = 23.238600 Å. b = 11.133400 Å, c = 45.862400 Å
	Volume of unit cell: 11651.776815 Å³

	# atoms = 960
	# charges = 960
	chemical formula: Dict(:H => 8, :S => 1, :Ca => 1, :O => 6, :C => 14)
	space Group: P1
	symmetry Operations:
		'x, y, z'
```

## finding other properties

```jldoctest crystal; output=false
rho = crystal_density(xtal)         # Crystal density of the crystal in kg/m^2
mw = molecular_weight(xtal)         # The molecular weight of the unit cell in amu
formula = chemical_formula(xtal)    # The irreducible chemical formula of the crystal
# output
Dict{Symbol, Int64} with 5 entries:
  :H  => 8
  :S  => 1
  :Ca => 1
  :O  => 6
  :C  => 14
```

## assigning new charges

If the crystal structure file does not contains partial charges, we provide methods to assign new charges to the crystal

```julia
species_to_charges = Dict(:Ca => 2.0, :C => 1.0, :H => -1.0)                # This method assigns a static charge to atom species
charged_xtal = assign_charges(xtal, species_to_charges, 1e-5)                # This function creates a new charged `Crystal` object.
                                                                            #   The function checks for charge neutrality with a tolerance of 1e-5
new_charges = Charges([2.0, 1.0, -1.0, -1.0, ...], xtal.atoms.coords)
other_charged_xtal = Crystal(xtal.name, xtal.box, xtal.atoms,               # Here we create a new `Charges` object using an array of new charges.
                             new_charges, xtal.bonds, xtal.symmetry)        #   The number of charges in the array has to be equal to the number of atoms
                                                                            #   and finally a new `Crystal` object is manually created
```

## writing crystal files

We provide methods to write both `.xyz` and `.cif` files

```jldoctest crystal; output=false
write_cif(xtal, "my_new_cif_file.cif")      # Stored in the current directory
write_xyz(xtal, "my_new_xyz_file.xyz")      # stored in the current directory
# output

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
