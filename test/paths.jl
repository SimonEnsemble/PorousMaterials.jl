module Path_Test

using PorousMaterials
using Test

@testset "Path Tests" begin
    
    set_path_to_crystals(joinpath(pwd(), "other_data", "other_crystals"))
    @eval PorousMaterials PATH_TO_MOLECULES = joinpath(pwd(), "other_data", "other_molecules")
    @eval PorousMaterials PATH_TO_FORCEFIELDS = joinpath(pwd(), "other_data", "other_forcefields")
    @eval PorousMaterials PATH_TO_GRIDS = joinpath(pwd(), "my_grids")
    Crystal("other_SBMOF-1.cif")
    Molecule("other_Xe")
    LJForceField("other_Dreiding")
    if ! isdir("my_grids")
        mkdir("my_grids")
    end
    cp("data/grids/test_grid.cube", "my_grids/other_test_grid.cube", force=true)
    read_cube("other_test_grid.cube")

    @test true # if made it this far :)
end
end
