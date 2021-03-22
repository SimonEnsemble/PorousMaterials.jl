using PorousMaterials
using OffsetArrays
using LinearAlgebra
using Test
using JLD2
using Statistics
using Random

run_tests = Dict("energy_min" => true, 
                 "grid_calc"  => false
                )

###
#  test find_energy_minimum
###
if run_tests["energy_min"]
    @info "running grid test"

    xtal = Crystal("SBMOF-1.cif")
    mol  = Molecule("Xe")
    ljff = LJForceField("UFF")
    mesh = (50, 50, 50)

    rep_factors = replication_factors(xtal, ljff)
    if run_tests["grid_calc"]
        grid = energy_grid(xtal, mol, ljff)
        grid_min = findmin(grid.data)

        # results of energy_grid calc 
        test_energy = grid_min[1] # kJ/mol
        text_vox_id = Tuple([grid_min[2][i] for i in 1:3]) # Cartesian Index of voxel
        test_xyz    = id_to_xf(text_vox_id, mesh) # fractional coords of voxel
    else
        # results from energy_grid calc with n_pts=(50, 50, 50)
        test_energy = -37.65575654162971 # kJ/mol
        text_vox_id = Tuple([26, 28, 26])
        test_xyz = [0.5102040816326531, 0.5510204081632653, 0.5102040816326531]
    end

    # run grid search without passing xf₀ will automatically run a coarse energy_grid 
    #   calc to get initial estimate with n_pts=(50, 50, 50)
    res_1 = find_energy_minimum(xtal, Frac(mol, xtal.box), ljff) 

    # rescale to original xtal
    res_1.minimizer .*= rep_factors

    @test res_1.minimum <= test_energy
    @test all(j -> (isapprox(res_1.minimizer[j], test_xyz[j], atol=1/mesh[j])), [1, 2, 3])

    # choose a different starting location and recover the energy and location
    xyz_2 = id_to_xf((text_vox_id[1], text_vox_id[2] + 2, text_vox_id[3]), mesh)
    res_2 = find_energy_minimum(xtal, Frac(mol, xtal.box), ljff; xf₀=Frac(xyz_2))

    # rescale to original xtal box
    res_2.minimizer .*= rep_factors

    @test isapprox(res_1.minimum, res_2.minimum)
    @test all(k -> (isapprox(res_1.minimizer[k], res_2.minimizer[k], atol=1/mesh[k])), [1, 2, 3])
end
