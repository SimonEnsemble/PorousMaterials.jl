module Bonds_test

using PorousMaterials
using LightGraphs
using Test

function visual_check(xtal::String)
    c = Crystal(xtal)
    c = replicate(c, (2, 2, 2))
    strip_numbers_from_atom_labels!(c)
    infer_geometry_based_bonds!(c, true) # must
    write_xyz(c)
    write_bond_information(c)
    @info c.name * " see .vtk and .xyz to visually check bonds"
end

@testset "Bond Tests" begin
    # NiPyC2 bonding
    c = Crystal("NiPyC2_relax.cif")
    strip_numbers_from_atom_labels!(c)
    c = replicate(c, (4, 4, 4))
    bonding_rules = [BondingRule(:H, :*, 0.4, 1.2),
                     BondingRule(:Ni, :O, 0.4, 2.5),
                     BondingRule(:Ni, :N, 0.4, 2.5),
                     BondingRule(:*, :*, 0.4, 1.9)]
    infer_bonds!(c, true, bonding_rules)
    g1 = deepcopy(c.bonds)
    conn_comps = connected_components(c.bonds)
    @test length(conn_comps) == 2 # interpenetrated
    c_red = getindex(c, conn_comps[1])
    c_blue = getindex(c, conn_comps[2])
    @test ne(c_red.bonds) + ne(c_blue.bonds) == ne(c.bonds)
    remove_bonds!(c)
    infer_geometry_based_bonds!(c, true)
    @test length(conn_comps) == 2 # interpenetrated
    @test g1 == c.bonds # consistency between two different bonding schemes
    
    # FIQCEN bonding
    c = Crystal("FIQCEN_clean.cif")
    strip_numbers_from_atom_labels!(c)
    infer_geometry_based_bonds!(c, true)
    @test length(connected_components(c.bonds)) == 1 # not interpenetrated
    @test c.atoms.species[neighbors(c.bonds, 1)] == [:Cu, :O, :O, :O, :O]
    visual_check("FIQCEN_clean.cif")
    # reduce covalant radii to see Cu-O bond disappear
    cov_radii = cordero_covalent_atomic_radii()
    cov_radii[:Cu] = 1.125
    remove_bonds!(c)
    @test ne(c.bonds) == 0
    infer_geometry_based_bonds!(c, true, covalent_radii=cov_radii)
    @test c.atoms.species[neighbors(c.bonds, 1)] == [:O, :O, :O, :O]
end
end
