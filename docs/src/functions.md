#Functions
This page contains all of the functions exported by PorousMaterials. They are sorted by the .jl files they are found in.

##Box.jl
```@docs
    Box
    replicate
    UnitCube
    write_vtk
```

##Crystal.jl
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

##ElectrostaticsEnergetics.jl
```@docs
    electrostatic_potential_energy
    precompute_kvec_wts
    setup_Ewald_sum
    total_electrostatic_potential_energy
```

##Energetics_Util.jl
```@docs
    PotentialEnergy
    SystemPotentialEnergy
```

##EOS.jl
```@docs
    PengRobinsonsGas
    calculate_properties
```

##Forcefield.jl
```@docs
    LJForcefield
    replication_factors
    check_forcefield_coverage
```

##GCMC.jl
```@docs
    gcmc_simulation
    adsorption_isotherm
    stepwise_adsorption_isotherm
    gcmc_result_savename
```

##Grid.jl
```@docs
    Grid
    apply_periodic_boundary_condition
    write_cube
    read_cube
    energy_grid
```

##Henry.jl
```@docs
    henry_coefficient
    henry_result_savename
```

##Matter.jl
```@docs
    LJSphere
    PtCharge
```

##MChelpers.jl
```@docs
    insert_molecule!
    delete_molecule!
    translate_molecule!
    reinsert_molecule!
    rotatable
```

##Misc.jl
```@docs
    read_xyz
    read_cpk_colors
    read_atomic_radii
    read_atomic_masses
    write_to_xyz
```

##Molecules.jl
```@docs
    translate_to!
    rotate!
    rotation_matrix
    rand_point_on_unit_sphere
    charged(::Molecule, ::Bool)
```

##NearestImage.jl
```@docs
    nearest_image!
    nearest_r²
    nearest_r
```

##VdWEnergetics.jl
```@docs
    lennard_jones
    vdw_energy
    vdw_energy_no_PBC
```
