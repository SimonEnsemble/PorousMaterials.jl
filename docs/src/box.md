```@meta
DocTestSetup = quote
  using PorousMaterials
end
```

# The Spatial Box

Within `PorousMaterials.jl`, the 3D space in which all `Coords` are located is the
`Box`.  Each `Crystal` has its own `Box`, equivalent to the unit cell of a material,
containing as attributes the unit cell edge lengths (`a` `b` `c`), crystallographic
dihedral angles (`α` `β` `γ`), volume, conversion factors for translating
between `Frac`tional and `Cart`esian coordinates, and the reciprocal (Fourier
transform) vectors for the Bravais lattice.

## defining a box

A `Box` is most conveniently constructed from its basic spatial data (`a` `b` `c`
`α` `β` `γ`).  For example, given the unit cell of Co-MOF-74, we can define its `Box`:

```jldoctest
a = 26.13173 # Å
b = 26.13173
c = 6.722028
α = π/2 # radians
β = π/2
γ = 2*π/3
box = Box(a, b, c, α, β, γ)
# output
Bravais unit cell of a crystal.
    Unit cell angles α = 90.000000 deg. β = 90.000000 deg. γ = 120.000000 deg.
    Unit cell dimensions a = 26.131730 Å. b = 26.131730 Å, c = 6.722028 Å
    Volume of unit cell: 3975.275878 Å³
```

A `Box` may also be defined by providing only the `Frac`tional-to-`Cart`esian conversion
matrix:
```jldoctest box
box = Box([26.1317 -13.0659 0; 0 22.6307 0; 0 0 6.72203])
# output
Bravais unit cell of a crystal.
    Unit cell angles α = 90.000000 deg. β = 90.000000 deg. γ = 120.000113 deg.
    Unit cell dimensions a = 26.131700 Å. b = 26.131711 Å, c = 6.722030 Å
    Volume of unit cell: 3975.265115 Å³
```

To quickly get a simple unit-cubic `Box`, use the `unit_cube` function.
```jldoctest
unit_cube()
# output
Bravais unit cell of a crystal.
    Unit cell angles α = 90.000000 deg. β = 90.000000 deg. γ = 90.000000 deg.
    Unit cell dimensions a = 1.000000 Å. b = 1.000000 Å, c = 1.000000 Å
    Volume of unit cell: 1.000000 Å³
```


## transforming coordinates

Conversions are provided for switching between `Frac`tional and `Cart`esian `Coords`
using the `Box` (works for `Atoms` and `Charges`, too)

```jldoctest box
xtal = Crystal("SBMOF-1.cif")
Cart(xtal.atoms.coords, xtal.box)
# output
Cart([4.594867082350715 -0.952720283971488 … 0.8392490029633858 -1.5321086078257065; 1.4395486200000005 4.2228986200000005 … 1.4289162230000012 4.212266223; 5.89964228469024 5.359217037237699 … 17.537474811394276 16.239103154389543])
```


## replicating a box

For simulations in larger volumes than a single crystallograhic unit cell, the
`Box` may be replicated along each or any of the three crystallographic axes.

```jldoctest box
replicated_box = replicate(box, (2,2,2))
# output
Bravais unit cell of a crystal.
    Unit cell angles α = 90.000000 deg. β = 90.000000 deg. γ = 120.000113 deg.
    Unit cell dimensions a = 52.263400 Å. b = 52.263422 Å, c = 13.444060 Å
    Volume of unit cell: 31802.120923 Å³
```


## exporting a box

For visualization of the unit cell boundaries, the `Box` may be written out to a
`.vtk` file for use in [Visit](https://wci.llnl.gov/simulation/computer-codes/visit/)

```jldoctest box
write_vtk(box, "box.vtk")
# output

```


# detailed docs

```@docs
    Box
    unit_cube
    replicate
    write_vtk
    Frac
```
