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

```julia
a = 26.13173 # Å
b = 26.13173
c = 6.722028
α = π/2 # radians
β = π/2
γ = 2*π/3
box = Box(a, b, c, α, β, γ)
```

A `Box` may also be defined by providing only the `Frac`tional-to-`Cart`esian conversion
matrix:
```julia
box = Box([26.1317 -13.0659 0; 0 22.6307 0; 0 0 6.72203])
```

To quickly get a simple unit-cubic `Box`, use the `unit_cube` function.
```julia
@info unit_cube()
#┌ Info: Bravais unit cell of a crystal.
#│       Unit cell angles α = 90.000000 deg. β = 90.000000 deg. γ = 90.000000 deg.
#│       Unit cell dimensions a = 1.000000 Å. b = 1.000000 Å, c = 1.000000 Å
#└       Volume of unit cell: 1.000000 Å³
```


## transforming coordinates

Conversions are provided for switching between `Frac`tional and `Cart`esian `Coords`
using the `Box` (works for `Atoms` and `Charges`, too)

```julia
xtal = Crystal("Co-MOF-74.cif")
Cart(xtal.atoms.coords, xtal.box)
#Cart([-5.496156112249995 7.181391379950001 … 15.131970232450003 2.4686645331000063;
# 22.270234304380295 2.8331425940892103 … 0.7607701110682343 22.13256395706254;
# 1.231811631 0.32198514120000005 … 6.2082409932000004 2.2119953472])
```


## replicating a box

For simulations in larger volumes than a single crystallograhic unit cell, the
`Box` may be replicated along each or any of the three crystallographic axes.

```julia
replicated_box = replicate(box, (2,2,2))
```


## exporting a box

For visualization of the unit cell boundaries, the `Box` may be written out to a
`.vtk` file for use in [Visit](https://wci.llnl.gov/simulation/computer-codes/visit/)

```julia
write_vtk(box, "box.vtk")
```


# detailed docs

```@docs
    Box
    unit_cube
    replicate
    write_vtk
    Frac
```
