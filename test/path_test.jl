module Path_Test

using PorousMaterials
using Test

@testset "Path Tests" begin
    set_tutorial_mode()
    @test PorousMaterials.PATH_TO_DATA == joinpath(dirname(pathof(PorousMaterials)),"..","test","data")
    # recommended way to change path to files
    @eval PorousMaterials PATH_TO_CRYSTALS = "blah"
    @test PorousMaterials.PATH_TO_CRYSTALS == "blah"
end
end
