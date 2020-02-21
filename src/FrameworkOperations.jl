# file containing Framework operations
function Base.isapprox(f1::Framework, f2::Framework)
    box_flag = isapprox(f1.box, f2.box)
    if f1.charges.n_charges != f2.charges.n_charges
        return false
    end
    if f1.atoms.n_atoms != f2.atoms.n_atoms
        return false
    end
    charges_flag = isapprox(f1.charges, f2.charges)
    atoms_flag = isapprox(f1.atoms, f2.atoms)
    symmetry_flag = is_symmetry_equal(f1.symmetry, f2.symmetry)
    return box_flag && charges_flag && atoms_flag && symmetry_flag
end

function Base.:+(frameworks::Framework...; check_overlap::Bool=true)
    new_framework = deepcopy(frameworks[1])
    for (i, f) in enumerate(frameworks)
        if i == 1
            continue
        end
        @assert isapprox(new_framework.box, f.box) @sprintf("Framework %s has a different box\n", f.name)
        @assert is_symmetry_equal(new_framework.symmetry, f.symmetry) @sprintf("Framework %s has different symmetry rules\n", f.name)
        @assert new_framework.space_group == f.space_group

        new_atoms = new_framework.atoms + f.atoms
        new_charges = new_framework.charges + f.charges

        nf_n_atoms = new_framework.atoms.n_atoms
        new_bonds = SimpleGraph(nf_n_atoms + f.atoms.n_atoms)
        for edge in collect(edges(new_framework.bonds))
            if !add_edge!(new_bonds, edge.src, edge.dst)
                @warn @sprintf("Edge %d->%d in framework %s is not being added correctly\nThe new edge should be %d->%d",
                               edge.src, edge.dst, f.name, nf_n_atoms + edge.src, nf_n_atoms + edge.dst)
            end
        end
	for edge in collect(edges(f.bonds))
            if !add_edge!(new_bonds, nf_n_atoms + edge.src, nf_n_atoms + edge.dst)
                @warn @sprintf("Edge %d->%d in framework %s is not being added correctly\nThe new edge should be %d->%d",
                               edge.src, edge.dst, f.name, nf_n_atoms + edge.src, nf_n_atoms + edge.dst)
            end
        end

        new_framework = Framework(split(new_framework.name, ".")[1] * "_" * split(f.name, ".")[1],
                                 new_framework.box, new_atoms, new_charges,
                                 symmetry=new_framework.symmetry,space_group=new_framework.space_group,
                                 is_p1=new_framework.is_p1, bonds=new_bonds)
    end
    if check_overlap
        if atom_overlap(new_framework)
            @warn "This new framework has overlapping atoms, use:\n`remove_overlapping_atoms_and_charges(framework)`\nto remove them"
        end

        if charge_overlap(new_framework)
            @warn "This new framework has overlapping charges, use:\n`remove_overlapping_atoms_and_charges(framework)`\nto remove them"
        end
    end

    return new_framework
end

"""
    partitions = partition_framework(framework, [[<indices>], [<indices>], ...])

Partitions a framework based on atom ids passed in with arrays. Each array of
indices will be create an additional partition being created. 

# Arguments
- `framework::Framework`: The framework to be partitioned
- `partition_atom_ids::Array{Array{Int, 1}, 1}`: Each element in this array is
    an array of atom indices to be included in a given partition. The number of
    arrays in this element is the number of partitions that will be returned
- `complete_partition::Bool`: Whether or not every atom in the framework must be
    included in a partition. i.e. if some atoms can be excluded from the
    resulting partitions.
- `atoms_in_multiple_partitions::Bool`: Determines if a specific atom can be
    included in more than one partition
"""
function partition_framework(framework::Framework,
                             all_partition_ids::Array{Array{Int, 1}, 1};
                             complete_partition::Bool=true,
                             atoms_in_multiple_partitions::Bool=false)
    if complete_partition && !issetequal(collect(1:framework.atoms.n_atoms),
                                         unique(vcat(all_partition_ids...)))
        error("Not all atoms in %s are accounted for in the partition ids\n\trun again with `complete_partition=false` to only partition some atoms")
    end

    if !atoms_in_multiple_partitions && !allunique(vcat(all_partition_ids...))
        error("Atoms appear in multiple partitions\n\trun again with `atoms_in_multiple_partitions=false` to allow atoms to appear in more than one partition")
    end

    # charges are not paired with atoms, so will only partition when there are
    #   no charges present
    @assert !charged(framework) "Cannot partition a charged framework"
    framework_partitions = []
    for (i, partition_ids) in enumerate(all_partition_ids)
        atoms_xf = framework.atoms.xf[:, partition_ids]
        atoms_species = framework.atoms.species[partition_ids]
        atoms = Atoms(atoms_species, atoms_xf)
        charges = Charges(Array{Float64, 1}(undef, 0),
                          Array{Float64, 2}(undef, 3, 0))
        old_to_new = Dict(partition_ids[i] => i for i in 1:length(partition_ids))
        bonds = SimpleGraph(length(partition_ids))
        for edge in collect(edges(framework.bonds))
            if edge.src in partition_ids && edge.dst in partition_ids
                add_edge!(bonds, old_to_new[edge.src], old_to_new[edge.dst])
            end
        end
        partitioned_framework = Framework(framework.name * "_partition_" * string(i),
                                        framework.box, atoms, charges; bonds=bonds)
        push!(framework_partitions, partitioned_framework)
    end
    #if setdiff(collect(edges(framework.bonds)), [collect(edges(framework_partitions[i].bonds)) for i in 1:length(all_partition_ids)]...) != []
    if ne(framework.bonds) != sum(ne.([framework_partitions[i].bonds for i in 1:length(framework_partitions)])) && 
        complete_partition && !atoms_in_multiple_partitions
        @warn @sprintf("Some bonds present in %s are between partitions and have been lost\n", framework.name)
    end

    if !complete_partition || atoms_in_multiple_partitions
        @warn "Some bonds may be missing from the original structure"
    end
    return framework_partitions
end

"""
    updated_framework = subtract_atoms(framework, ids_to_remove)

removes the atoms given from the framework. It then updates the bond graph
accordingly. This does not work for charged frameworks
"""
function subtract_atoms(framework::Framework, ids_to_remove::Array{Int, 1})
    # does not allow charged frameworks
    @assert !charged(framework)
    # get array of ids to keep
    ids_to_keep = setdiff(collect(1:framework.atoms.n_atoms), ids_to_remove)
    # create dictionary to get new atom id from old atom id
    old_to_new = Dict(ids_to_keep[i] => i for i in 1:length(ids_to_keep))
    # get the atom locations
    atoms_xf = framework.atoms.xf[:, ids_to_keep]
    atoms_species = framework.atoms.species[ids_to_keep]
    atoms = Atoms(atoms_species, atoms_xf)
    # create emtpy charges group
    charges = Charges(Array{Float64, 1}(undef, 0),
                      Array{Float64, 2}(undef, 3, 0))
    # create the new bonding graph
    bonds = SimpleGraph(length(ids_to_keep))
    for edge in collect(edges(framework.bonds))
        if edge.src in ids_to_keep && edge.dst in ids_to_keep
            add_edge!(bonds, old_to_new[edge.src], old_to_new[edge.dst])
        end
    end

    return Framework("removed_atoms_" * framework.name, framework.box, atoms,
                     charges, bonds=bonds, symmetry=framework.symmetry,
		     space_group=framework.space_group, is_p1=framework.is_p1)
end
