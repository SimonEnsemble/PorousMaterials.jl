## Boxes
```@docs
    Box
    replicate
    UnitCube
    write_vtk
```

## Crystals
```@docs
    Framework
    remove_overlapping_atoms_and_charges
    strip_numbers_from_atom_labels!
    chemical_formula
    molecular_weight
    crystal_density
    replicate(::Framework, ::Tuple{Int, Int, Int})
    charged(::Framework; ::Bool)
    write_cif
    assign_charges
```

## Grids
```@docs
    Grid
    apply_periodic_boundary_condition!
    write_cube
    read_cube
    energy_grid
```
