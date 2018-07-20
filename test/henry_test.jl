using PorousMaterials
using Base.Test

insertions_per_volume = 1

@testset "Henry coefficient tests" begin
###
#  Henry test 1: Xe in SBMOF-1
###
framework = read_crystal_structure_file("SBMOF-1.cif")
ljff = read_forcefield_file("Dreiding.csv", cutoffradius=12.5)
molecule = read_molecule_file("Xe")
temperature = 298.0

result = henry_coefficient(framework, molecule, temperature, ljff, 
                           insertions_per_volume=insertions_per_volume, verbose=true)
@test isapprox(result["henry coefficient [mol/(kg-Pa)]"], 0.00985348, atol=0.0002)
@test isapprox(result["⟨U⟩ (kJ/mol)"], -39.6811, atol=0.1)

###
#  Henry test 2: CO2 in CAXVII_clean.cif
###
framework = read_crystal_structure_file("CAXVII_clean.cif")
ljff = read_forcefield_file("Dreiding.csv", cutoffradius=12.5)
molecule = read_molecule_file("CO2")
temperature = 298.0

result = henry_coefficient(framework, molecule, temperature, ljff, 
                           insertions_per_volume=insertions_per_volume, verbose=true)
@test isapprox(result["henry coefficient [mol/(kg-Pa)]"], 2.88317e-05, atol=1.5e-7)
@test isapprox(result["⟨U⟩ (kJ/mol)"], -18.69582223, atol=0.1)

end
