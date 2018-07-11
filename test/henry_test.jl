using PorousMaterials
using Base.Test
 
###
#  Henry test 1: Xe in SBMOF-1
###
framework = read_crystal_structure_file("SBMOF-1.cif")
ljff = read_forcefield_file("Dreiding.csv", cutoffradius=12.5)
molecule = read_molecule_file("Xe")
temperature = 298.0

result = henry_coefficient(framework, molecule, temperature, ljff, nb_insertions=1000000, verbose=true)
@test isapprox(result["henry coefficient [mol/(kg-Pa)]"], 0.00985348, atol=0.000135146)
@test isapprox(result["⟨U⟩ (kJ/mol)"], -39.6811, atol=0.03)

###
#  Henry test 2: CO2 in CAXVII_clean.cif
###
framework = read_crystal_structure_file("CAXVII_clean.cif")
ljff = read_forcefield_file("Dreiding.csv", cutoffradius=12.5)
molecule = read_molecule_file("CO2")
temperature = 298.0

result = henry_coefficient(framework, molecule, temperature, ljff, nb_insertions=100000, verbose=true)
warn("waiting for CO2 test data")
 # @test isapprox(result["henry coefficient [mol/(kg-Pa)]"], 0.00985348, atol=0.000135146)
 # @test isapprox(result["Qst (kJ/mol)"], -39.6811, atol=0.03)
