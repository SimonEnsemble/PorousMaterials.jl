module Path_Test

using PorousMaterials
using Test

@testset "Path Tests" begin
    set_tutorial_mode()
    @test PorousMaterials.PATH_TO_DATA == joinpath(dirname(pathof(PorousMaterials)),"..","test","data")
    # recommended way to change path to files
    @eval PorousMaterials PATH_TO_CRYSTALS = "blah"
    @test PorousMaterials.PATH_TO_CRYSTALS == "blah"

    # set back to default for testing below
    set_default_file_paths()

    # test changing paths and reading information from new paths
    # testing files were copied over directly from the original data folder, so
    #   they should match
    # MOLECULES
    old_molecule_path = Molecule("Xe")
    @eval PorousMaterials PATH_TO_MOLECULES = joinpath(pwd(), "other_data", "other_molecules")
    new_molecule_path = Molecule("Xe")
    @test isapprox(old_molecule_path, new_molecule_path)
    # CRYSTALS
    old_crystal_path = Framework("SBMOF-1.cif")
    @eval PorousMaterials PATH_TO_CRYSTALS = joinpath(pwd(), "other_data", "other_crystals")
    new_crystal_path = Framework("SBMOF-1.cif")
    @test isapprox(old_crystal_path, new_crystal_path)
    # FORCEFIELDS
    old_forcefield_path = LJForceField("Dreiding.csv")
    @eval PorousMaterials PATH_TO_FORCEFIELDS = joinpath(pwd(), "other_data", "other_forcefields")
    new_forcefield_path = LJForceField("Dreiding.csv")
    # TODO make an approximation function or find something like `@test_throws NoError`
    # GRIDS
    old_grid_path = read_cube("test_grid.cube")
    if ! isdir(joinpath("other_data", "other_grids"))
        mkdir(joinpath("other_data", "other_grids"))
    end
    run(`cp data/grids/test_grid.cube other_data/other_grids/test_grid.cube`)
    @eval PorousMaterials PATH_TO_GRIDS = joinpath(pwd(), "other_data", "other_grids")
    new_grid_path = read_cube("test_grid.cube")
    @test isapprox(old_grid_path, new_grid_path, atol=1e-5)
end
end
