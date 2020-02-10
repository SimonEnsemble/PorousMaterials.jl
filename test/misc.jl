module Misc_Test

using PorousMaterials
using Test
using DataFrames
using CSV
using Optim
using Random

@testset "Misc Tests" begin
    am = read_atomic_masses()
    @test isapprox(am[:H], 1.00794, atol=0.001)
    @test isapprox(am[:Co], 58.9332, atol=0.001)
    
    test_xyz_filename = "atoms_test"
    c = Cart([1.0 4.0;
              2.0 5.0;
              3.0 6.0]
             )
    s = [:C, :H]
    atoms = Atoms(s, c)
    write_xyz(atoms, test_xyz_filename)
    atoms_read = read_xyz(test_xyz_filename * ".xyz")
    @test isapprox(atoms, atoms_read)
    rm(test_xyz_filename * ".xyz")

    @test read_cpk_colors()[:Li] == (204,128,255)
end
end
