```@meta
DocTestSetup = quote
  using PorousMaterials
end
```

# Lennard-Jones Force Fields and Potential Energy

Lennard-Jones force field parameters are stored in comma-separated-value format in `rc[:paths][:forcefields]`.

Interaction of an adsorbate with the crystal is modeled as pair-wise additive and with Lennard-Jones potentials of the form:

`V(r) = 4 * ϵ * [ x ^ 12 - x ^ 6 ]`, where `x = σ / r`

The Lennard-Jones force field input files, e.g. `UFF.csv` contain a list of pure (i.e. X-X, where X is an atom) sigmas (σ) and epsilons (ϵ) with units Angstrom (Å) and Kelvin (K), respectively. Note that, e.g., in the [UFF paper](https://doi.org/10.1021/ja00051a040), the Lennard-Jones potential is written in a different form; thus, parameters need to be converted to correspond to the functional form used in `PorousMaterials.jl`.

## Building Blocks of PorousMaterials: Lennard-Jones Force Fields

### Loading Force Field Files and Accessing Attributes

Reading in Lennard-Jones force field parameters is made easy with the [`LJForceField`](@ref) function. Let's load in the parameters from the Universal Force Field file (`UFF.csv`):

```jldoctest force_field
# read in Lennard-Jones force field parameters from the Universal Force Field
ljforcefield = LJForceField("UFF"; r_cutoff=14.0, mixing_rules="Lorentz-Berthelot")

# output

Force field: UFF
Number of atoms included: 108
Cut-off radius (Å) = 14.0
   Ra-   Ra ϵ =  203.30076 K, σ =    3.27583 Å
   Cl-   Cl ϵ =  114.23087 K, σ =    3.51638 Å
   Al-   Al ϵ =  254.12595 K, σ =    4.00815 Å
   Be-   Be ϵ =   42.77367 K, σ =    2.44552 Å
   Re-   Re ϵ =   33.21250 K, σ =    2.63171 Å
   Cr-   Cr ϵ =    7.54830 K, σ =    2.69319 Å
   Na-   Na ϵ =   15.09659 K, σ =    2.65755 Å
   Sb-   Sb ϵ =  225.94565 K, σ =    3.93777 Å
   Cf-   Cf ϵ =    6.54186 K, σ =    2.95155 Å
   Kr-   Kr ϵ =  110.70833 K, σ =    3.68921 Å
   Ni-   Ni ϵ =    7.54830 K, σ =    2.52481 Å
    S-    S ϵ =  137.88220 K, σ =    3.59478 Å
  CH4-  CH4 ϵ =  148.00000 K, σ =    3.73000 Å
   Fm-   Fm ϵ =    6.03864 K, σ =    2.92749 Å
   Ru-   Ru ϵ =   28.18030 K, σ =    2.63973 Å
   Tl-   Tl ϵ =  342.18940 K, σ =    3.87274 Å
   re-   re ϵ =    0.00010 K, σ =    5.00000 Å
   Tm-   Tm ϵ =    3.01932 K, σ =    3.00589 Å
C_CO2-C_CO2 ϵ =   27.00000 K, σ =    2.80000 Å
    W-    W ϵ =   33.71572 K, σ =    2.73417 Å
    O-    O ϵ =   30.19318 K, σ =    3.11815 Å
   Nd-   Nd ϵ =    5.03220 K, σ =    3.18496 Å
   Tb-   Tb ϵ =    3.52254 K, σ =    3.07449 Å
   Th-   Th ϵ =   13.08371 K, σ =    3.02549 Å
   Zr-   Zr ϵ =   34.72216 K, σ =    2.78317 Å
    F-    F ϵ =   25.16099 K, σ =    2.99698 Å
   Co-   Co ϵ =    7.04508 K, σ =    2.55866 Å
   Fr-   Fr ϵ =   25.16099 K, σ =    4.36540 Å
   Gd-   Gd ϵ =    4.52898 K, σ =    3.00055 Å
   Rh-   Rh ϵ =   26.67064 K, σ =    2.60944 Å
   Pu-   Pu ϵ =    8.05152 K, σ =    3.05044 Å
   Lw-   Lw ϵ =    5.53542 K, σ =    2.88295 Å
   Ar-   Ar ϵ =   93.09564 K, σ =    3.44600 Å
   Ca-   Ca ϵ =  119.76629 K, σ =    3.02816 Å
   Cm-   Cm ϵ =    6.54186 K, σ =    2.96313 Å
    N-    N ϵ =   34.72216 K, σ =    3.26069 Å
   As-   As ϵ =  155.49489 K, σ =    3.76850 Å
   Yb-   Yb ϵ =  114.73409 K, σ =    2.98897 Å
   Se-   Se ϵ =  146.43693 K, σ =    3.74623 Å
    Y-    Y ϵ =   36.23182 K, σ =    2.98006 Å
   Am-   Am ϵ =    7.04508 K, σ =    3.01213 Å
   Pt-   Pt ϵ =   40.25758 K, σ =    2.45354 Å
    I-    I ϵ =  170.59148 K, σ =    4.00904 Å
   Fe-   Fe ϵ =    6.54186 K, σ =    2.59430 Å
   Ba-   Ba ϵ =  183.17197 K, σ =    3.29900 Å
   Hf-   Hf ϵ =   36.23182 K, σ =    2.79831 Å
   Es-   Es ϵ =    6.03864 K, σ =    2.93907 Å
   Po-   Po ϵ =  163.54640 K, σ =    4.19524 Å
   Eu-   Eu ϵ =    4.02576 K, σ =    3.11191 Å
    C-    C ϵ =   52.83807 K, σ =    3.43085 Å
   Zn-   Zn ϵ =   62.39924 K, σ =    2.46155 Å
   Cs-   Cs ϵ =   22.64489 K, σ =    4.02419 Å
   Mn-   Mn ϵ =    6.54186 K, σ =    2.63795 Å
   Rn-   Rn ϵ =  124.79849 K, σ =    4.24513 Å
   Bk-   Bk ϵ =    6.54186 K, σ =    2.97471 Å
   Ir-   Ir ϵ =   36.73504 K, σ =    2.53015 Å
   Rb-   Rb ϵ =   20.12879 K, σ =    3.66516 Å
   In-   In ϵ =  301.42860 K, σ =    3.97608 Å
   Hg-   Hg ϵ =  193.73958 K, σ =    2.40988 Å
   Te-   Te ϵ =  200.28144 K, σ =    3.98232 Å
   At-   At ϵ =  142.91439 K, σ =    4.23177 Å
   Bi-   Bi ϵ =  260.66780 K, σ =    3.89323 Å
   Cu-   Cu ϵ =    2.51610 K, σ =    3.11369 Å
   Tc-   Tc ϵ =   24.15455 K, σ =    2.67091 Å
   Sn-   Sn ϵ =  285.32557 K, σ =    3.91283 Å
   Pa-   Pa ϵ =   11.07083 K, σ =    3.05044 Å
   Lu-   Lu ϵ =   20.63201 K, σ =    3.24287 Å
   Mo-   Mo ϵ =   28.18030 K, σ =    2.71902 Å
   Ac-   Ac ϵ =   16.60625 K, σ =    3.09855 Å
    U-    U ϵ =   11.07083 K, σ =    3.02460 Å
   Li-   Li ϵ =   12.58049 K, σ =    2.18359 Å
   Er-   Er ϵ =    3.52254 K, σ =    3.02104 Å
   Ta-   Ta ϵ =   40.76080 K, σ =    2.82415 Å
S_H2S-S_H2S ϵ =  122.00000 K, σ =    3.60000 Å
   Cd-   Cd ϵ =  114.73409 K, σ =    2.53728 Å
   Os-   Os ϵ =   18.61913 K, σ =    2.77960 Å
   Ti-   Ti ϵ =    8.55473 K, σ =    2.82860 Å
    B-    B ϵ =   90.57955 K, σ =    3.63754 Å
    V-    V ϵ =    8.05152 K, σ =    2.80099 Å
H_H2S-H_H2S ϵ =   50.00000 K, σ =    2.50000 Å
   Si-   Si ϵ =  202.29432 K, σ =    3.82641 Å
   Ga-   Ga ϵ =  208.83618 K, σ =    3.90481 Å
   Au-   Au ϵ =   19.62557 K, σ =    2.93373 Å
   Mg-   Mg ϵ =   55.85739 K, σ =    2.69141 Å
    K-    K ϵ =   17.61269 K, σ =    3.39611 Å
   Ag-   Ag ϵ =   18.11591 K, σ =    2.80455 Å
   Sc-   Sc ϵ =    9.56117 K, σ =    2.93551 Å
   Ge-   Ge ϵ =  190.72027 K, σ =    3.81305 Å
   Nb-   Nb ϵ =   29.68996 K, σ =    2.81969 Å
   Ce-   Ce ϵ =    6.54186 K, σ =    3.16804 Å
   Pm-   Pm ϵ =    4.52898 K, σ =    3.16002 Å
   Pd-   Pd ϵ =   24.15455 K, σ =    2.58272 Å
   Dy-   Dy ϵ =    3.52254 K, σ =    3.05400 Å
   Sr-   Sr ϵ =  118.25663 K, σ =    3.24376 Å
   Ho-   Ho ϵ =    3.52254 K, σ =    3.03707 Å
   No-   No ϵ =    5.53542 K, σ =    2.89364 Å
O_CO2-O_CO2 ϵ =   79.00000 K, σ =    3.05000 Å
   Sm-   Sm ϵ =    4.02576 K, σ =    3.13596 Å
   Br-   Br ϵ =  126.30814 K, σ =    3.73197 Å
   Pb-   Pb ϵ =  333.63466 K, σ =    3.82819 Å
   Xe-   Xe ϵ =  167.06894 K, σ =    3.92352 Å
   Np-   Np ϵ =    9.56117 K, σ =    3.05044 Å
    P-    P ϵ =  153.48201 K, σ =    3.69456 Å
   La-   La ϵ =    8.55473 K, σ =    3.13775 Å
   Md-   Md ϵ =    5.53542 K, σ =    2.91680 Å
    H-    H ϵ =   22.14167 K, σ =    2.57113 Å
   Pr-   Pr ϵ =    5.03220 K, σ =    3.21258 Å
   He-   He ϵ =   28.18319 K, σ =    2.10430 Å
```

This also prints all of the atoms included in the loaded forcefield with their given ϵ and σ. This was excluded because it would use too much space on this page.

We can access attributes `LJForceField` such as `pure_σ`, `pure_ϵ`, and interaction values:

```jldoctest force_field; output=false
# access the Lennard-Jones epsilon & sigma for Xe 
ljforcefield.pure_ϵ[:Xe] # K
ljforcefield.pure_σ[:Xe] # Å

# access the Lennard-Jones epsilon & sigma for Xe-C interactions
ljforcefield.ϵ[:Xe][:C]  # K
ljforcefield.σ²[:Xe][:C] # Å (store σ² for faster computation)

# output

13.521685546424905
```

### Checking Force Field Coverage

When running simulations, it is necessary to have the force field terms for all of the atoms. This can be checked using [`forcefield_coverage`](@ref):

```jldoctest force_field; output=false
# check is the atoms in a crystal are covered
xtal = Crystal("SBMOF-1.cif")
forcefield_coverage(xtal, ljforcefield)

# check if the atoms in a molecule are covered
molecule = Molecule("CO2")
forcefield_coverage(molecule, ljforcefield)

# output

true
```

### Simulation Box and the Cutoff Radius

Find the replication factors needed to make a supercell big enough to fit a sphere with the specified cutoff radius.In PorousMaterials.jl, rather than replicating the atoms in the home unit cell to build the supercell that serves as a simulation box, we replicate the home unit cell to form the supercell (simulation box) in a for loop.The [`replication_factors`](@ref) function ensures enough replication factors such that the nearest image convention can be applied.

```jldoctest force_field
r_cutoff = 14.0 # Å
repfactors = replication_factors(xtal.box, r_cutoff)

# output

(3, 6, 2)
```

### Potential Energies: Van der Waals

What is the van der Waals potential energy of a Xe adsorbate inside SBMOF-1 at Cartesian coordinates `[0.0, 1.0, 3.0]` using the UFF as a molecular model?

```jldoctest force_field
# load molecule and convert it to fractional
molecule = Molecule("Xe")
molecule = Frac(molecule, xtal.box)
translate_to!(molecule, Cart([0.0, 1.0, 3.0]), xtal.box) # need box b/c we're in Cartesian
energy = vdw_energy(xtal, molecule, ljforcefield) # K

# output

5.73882798944654e6
```

# detailed docs

## Forcefields

```@docs
    LJForceField
    replication_factors
    forcefield_coverage
```

## Potential Energy

```@docs
    PotentialEnergy
    SystemPotentialEnergy
```

## Nearest Image Conventions

```@docs
    nearest_image!
```

## Electrostatics Energy

```@docs
    Eikr
    total
    electrostatic_potential
    electrostatic_potential_energy
    precompute_kvec_wts
    setup_Ewald_sum
    total_electrostatic_potential_energy
```

## Van der Waals Energy

```@docs
    lennard_jones
    vdw_energy
    vdw_energy_no_PBC
```
