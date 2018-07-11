using PorousMaterials
molecule = read_molecule_file("CO2")
framework = read_crystal_structure_file("SBMOF-1.cif")
ljff = read_forcefield_file("UFF.csv")
temperature = 298.0

henry_coeff = henry_coefficient(framework, molecule, temperature, ljff, nb_insertions=100)
