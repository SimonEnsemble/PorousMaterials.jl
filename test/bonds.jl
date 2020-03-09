module Bonds_test

using PorousMaterials
using LightGraphs
using Test

@testset "Bond Tests" begin
    c = Crystal("NiPyC2.cif")
    strip_numbers_from_atom_labels!(c)
    c = replicate(c, (4, 4, 4))
    bonding_rules = [BondingRule(:H, :*, 0.4, 1.2),
                     BondingRule(:Ni, :O, 0.4, 2.5),
                     BondingRule(:Ni, :N, 0.4, 2.5),
                     BondingRule(:*, :*, 0.4, 1.9)]
    infer_bonds!(c, true, bonding_rules)
    conn_comps = connected_components(c.bonds)
    @test length(conn_comps) == 2
    c_red = getindex(c, conn_comps[1])
    c_blue = getindex(c, conn_comps[2])
    @test ne(c_red.bonds) + ne(c_blue.bonds) == ne(c.bonds)

end
end
