
<a id='Boxes-1'></a>

## Boxes


```
    Box
    replicate
    UnitCube
    write_vtk
```


<a id='Crystals-1'></a>

## Crystals


```
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


<a id='Grids-1'></a>

## Grids


```
    Grid
    apply_periodic_boundary_condition!
    write_cube
    read_cube
    energy_grid
```

