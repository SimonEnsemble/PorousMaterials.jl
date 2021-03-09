# module Grid_Test

using PorousMaterials
using OffsetArrays
using LinearAlgebra
using Test
using JLD2
using Statistics
using Random
    

###
#  test robustness of find energy minimum
###
# grid params
xtal = replicate(Crystal("SBMOF-1.cif"), (1, 1, 1))
strip_numbers_from_atom_labels!(xtal)
mol  = Molecule("Xe")
ljff = LJForceField("UFF")
mesh = (100, 100, 100)
temp = 298.0
rot  = 750
# run energy grid calculation
grid = energy_grid(xtal, mol, ljff; n_pts=mesh, temperature=temp, n_rotations=rot)
# find location of the minumum energy on the grid
grid_min  = findmin(grid.data)
xyz_1     = Tuple([grid_min[2][i] for i in 1:3])
xf_minE_1 = id_to_xf(xyz_1, grid.n_pts) # fractional coords of min
# find the minimum energy using grid search optimization
res_1 = find_energy_minimum(xtal, Frac(mol, xtal.box), ljff, Frac(xf_minE_1))
# perturb molecule and re-apply grid search 



########
# 1. perform a coarse energy_grid calculation 
# 2. find the minumim energy location on that grid
# 3. apply grid search optimization to find "actual" minimum
# 4. move molecule and re-apply grid search optimization
# 5. test 1: is the energy minimum at the same location?
#    test 2: is the minimum energy the same value?
#
# How much should I perturb the probe molecule before passing it into the grid search again?
########
