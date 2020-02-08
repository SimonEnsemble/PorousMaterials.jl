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
        add_vertices!(new_framework.bonds, nf_n_atoms)
        for edge in collect(edges(f.bonds))
            add_edge!(new_framework.bonds, nf_n_atoms + edge.src, nf_n_atoms + edge.dst)
        end

        new_framework = Framework(split(new_framework.name, ".")[1] * "_" * split(f.name, ".")[1],
                                 new_framework.box, new_atoms, new_charges,
                                 symmetry=new_framework.symmetry,space_group=new_framework.space_group,
                                 is_p1=new_framework.is_p1, bonds=new_framework.bonds)
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
    # charges are not paired with atoms, so will only partition when there are
    #   no charges present
    @assert !charged(framework)
    framework_partitions = []
    for (i, partition_ids) in enumerate(all_partition_ids)
        atoms_xf = framework.atoms.xf[:, partition_ids]
        atoms_species = framework.atoms.species[partition_ids]
        atoms = Atoms(atoms_species, atoms_xf)
        charges = Charges(Array{Float64, 1}(undef, 0),
                          Array{Float64, 2}(undef, 3, 0))
        partitioned_framework = Framework(framework.name * "_partition_" * string(i),
                                        framework.box, atoms, charges)
        push!(framework_partitions, partitioned_framework)
    end
    return framework_partitions
end
