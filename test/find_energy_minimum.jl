using PorousMaterials
using OffsetArrays
using LinearAlgebra
using Test
using JLD2
using Statistics
using Random

# TODO:
# 1. hard code a test value
# 2. do a very fine energy_grid calc and test res.minimum < findmin(grid)[1]

###
#  test find_energy_minimum
###

#    # grid params
#    xtal = replicate(Crystal("SBMOF-1.cif"), (1, 1, 1)) 
#    strip_numbers_from_atom_labels!(xtal)
#    mol  = Molecule("Xe")
#    ljff = LJForceField("UFF")
#    mesh = (50, 50, 50) 
#    temp = 298.0
#    rot  = 750 
#
#    # run energy grid calculation
#    grid = energy_grid(xtal, mol, ljff; n_pts=mesh, temperature=temp, n_rotations=rot)
#
#    # find location of the minumum energy on the grid
#    grid_min  = findmin(grid.data)
#    xyz_1     = Tuple([grid_min[2][i] for i in 1:3])
#    xf_minE_1 = id_to_xf(xyz_1, grid.n_pts) # fractional coords of min
#
#    # find the minimum energy using grid search optimization
#    res_1 = find_energy_minimum(xtal, Frac(mol, xtal.box), ljff, Frac(xf_minE_1))
#
#    # perturb molecule and re-apply grid search 
#    xyz_2 = (xyz_1[1], xyz_1[2] + 2, xyz_1[3]) # move two voxels along y-axis
#    xf_minE_2 = id_to_xf(xyz_2, grid.n_pts)
#
#    # find the minimum energy using grid search optimization
#    res_2 = find_energy_minimum(xtal, Frac(mol, xtal.box), ljff, Frac(xf_minE_2))
#
#    # test the the value and location of minimum is the same
#    @test isapprox(res_1.minimum, res_2.minimum)
#    @test isapprox(res_1.minimizer, res_2.minimizer)

@info "running grid test"

###
#  remakeing test
###
xtal = Crystal("SBMOF-1.cif")
strip_numbers_from_atom_labels!(xtal)
mol  = Molecule("Xe")
ljff = LJForceField("UFF")
mesh = (50, 50, 50)
temp = 298.0

grid = energy_grid(xtal, mol, ljff; n_pts=mesh, temperature=temp)
grid_min = findmin(grid.data)

# results of energy_grid calc with n_pts=(250, 250, 250)
test_energy = grid_min[1] #-44.45608710623562 # kJ/mol
text_vox_id = Tuple([grid_min[2][i] for i in 1:3])  #(130, 141, 128) # Cartesian Index of voxel
test_xyz    = id_to_xf(text_vox_id, mesh) #id_to_xf(text_vox_id, (250, 250, 250))

@info "test_energy = $test_energy"
@info "test_xyz = $test_xyz"

# run grid search without passing xf₀ will automatically run a coarse energy_grid 
#   calc to get initial estimate with n_pts=(50, 50, 50)
res_1 = find_energy_minimum(xtal, Frac(mol, xtal.box), ljff) 

@info "res_1.minimum = $(res_1.minimum)"
@info "res_1.minimizer = $(res_1.minimizer)"

@test res_1.minimum < test_energy
@test all(i -> isapprox(res_1.minimizer[i], test_xyz[i], atol=1/mesh[i]), 1:length(mesh))

# choose a different starting location and recover the energy and location
xyz_2 = id_to_xf((text_vox_id[1], text_vox_id[2] + 3, text_vox_id[3]), (250, 250, 250))
res_2 = find_energy_minimum(xtal, Frac(mol, xtal.box), ljff; xf₀=Frac(xyz_2))

@test isapprox(res_1.minimum, res_2.minimum)
@test isapprox(res_1.minimizer, res_2.minimizer)

