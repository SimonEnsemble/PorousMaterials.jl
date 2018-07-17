using PorousMaterials
using Base.Test
    
framework = read_crystal_structure_file("SBMOF-1.cif")
ljff = read_forcefield_file("Dreiding.csv", cutoffradius=12.5)
molecule = read_molecule_file("Xe")
temperature = 298.0

result = henry_coefficient(framework, molecule, temperature, ljff, nb_insertions=10000000, verbose=true)
@test isapprox(result["henry coefficient [mol/(kg-Pa)]"], 0.00985348, atol=0.000135146)
, # 0.00985348 +/- 0.000135146 mol/kg/Pa
