using PorousMaterials
using LinearAlgebra
using LightGraphs
using Test

@testset "Subtraction Tests" begin
    f = Framework("SBMOF-1.cif")
    bonding_rules = [BondingRule(:*, :H, 0.4, 1.2),
                     BondingRule(:*, :*, 0.4, 1.9),
                     BondingRule(:Ca, :O, 0.4, 2.5)]
    infer_bonds!(f, false, bonding_rules)
    num_bonds_original = ne(f.bonds)
    all_ids = collect(1:f.atoms.n_atoms)
    # get ids for all S atoms so we can remove them
    ids_to_remove = filter(id -> f.atoms.species[id] == :S, all_ids)
    num_sulfur_bonds = 0
    for edge in collect(edges(f.bonds))
        if edge.src in ids_to_remove || edge.dst in ids_to_remove
            num_sulfur_bonds += 1
        end
    end

    # remove the sulfur from the framework
    new_f = subtract_atoms(f, ids_to_remove)

    @test num_bonds_original == num_sulfur_bonds + ne(new_f.bonds)

    # write out the original xyz and bonds plus the new xyz and bonds
    write_xyz(f)
    write_bond_information(f)

    write_xyz(new_f)
    write_bond_information(new_f)

 #      # test framework partitioning
 #      sbmof1 = Framework("SBMOF-1.cif")
 #      # create 3 partitions to test splitting the framework
 #      #   there are 120 atoms in sbmof1, so each partition will have 40 atoms
 #      lo = collect(1:40)
 #      mid = collect(41:80)
 #      hi = collect(81:120)
 #      sbmof1_partitions = partition_framework(sbmof1, [lo, mid, hi])
 #      @test length(lo) == sbmof1_partitions[1].atoms.n_atoms
 #      @test length(mid) == sbmof1_partitions[2].atoms.n_atoms
 #      @test length(hi) == sbmof1_partitions[3].atoms.n_atoms
 #      # add the framework partitions back together
 #      sbmof1_rebuilt = +(sbmof1_partitions...)
 #      @test isapprox(sbmof1, sbmof1_rebuilt)
 #  
 #      # test framework partitioning with the errors it should throw
 #      lo_hi = vcat(collect(1:20), collect(101:120))
 #      @test_throws ErrorException partition_framework(sbmof1, [lo, hi])
end
