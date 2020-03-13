"""
    bonding_rule = BondingRule(:Ca, :O, 0.4, 2.0)
    bonding_rules = [BondingRule(:H, :*, 0.4, 1.2),
                     BondingRule(:*, :*, 0.4, 1.9)]

A rule for determining if two atoms within a crystal are bonded. 

# Attributes
-`species_i::Symbol`: One of the atoms types for this bond rule
-`species_j::Symbol`: The other atom type for this bond rule
-`min_dist`: The minimum distance between the atoms for bonding to occur
-`max_dist`: The maximum distance between the atoms for bonding to occur
"""
struct BondingRule
    species_i::Symbol
    species_j::Symbol
    min_dist::Float64
    max_dist::Float64
end

"""
    default_bondingrules = default_bondingrules()

Returns the default bonding rules. Using `append!` and/or `prepend!` to add to the default bonding rules:

# Example
```
bond_rules = default_bondingrules()
prepend!(bond_rules, BondingRule(:Cu, :*, 0.1, 2.6))
```

# Returns
-`default_bondingrules::Array{BondingRule, 1}`: The default bonding rules: `[BondingRule(:*, :*, 0.4, 1.2), BondingRule(:*, :*, 0.4, 1.9)]`
"""
default_bondingrules() = [BondingRule(:H, :*, 0.4, 1.2), BondingRule(:*, :*, 0.4, 1.9)]

"""
    are_atoms_bonded = is_bonded(crystal, i, j, bonding_rules=[BondingRule(:H, :*, 0.4, 1.2), BondingRule(:*, :*, 0.4, 1.9)],
                                 include_bonds_across_periodic_boundaries=true)

Checks to see if atoms `i` and `j` in `crystal` are bonded according to the `bonding_rules`.

# Arguments
-`crystal::Crystal`: The crystal that bonds will be added to
-`i::Int`: Index of the first atom
-`j::Int`: Index of the second atom
-`bonding_rules::Array{BondingRule, 1}`: The array of bonding rules that will
    be used to fill the bonding information. They are applied in the order that
    they appear.
-`include_bonds_across_periodic_boundaries::Bool`: Whether to check across the
    periodic boundary when calculating bonds

# Returns
-`are_atoms_bonded::Bool`: Whether atoms `i` and `j` are bonded according to `bonding_rules`

"""
function is_bonded(crystal::Crystal, i::Int64, j::Int64,
                   bonding_rules::Array{BondingRule, 1}=[BondingRule(:H, :*, 0.4, 1.2), BondingRule(:*, :*, 0.4, 1.9)];
                   include_bonds_across_periodic_boundaries::Bool=true)
    species_i = crystal.atoms.species[i]
    species_j = crystal.atoms.species[j]
    
    r = distance(crystal.atoms, crystal.box, i, j, include_bonds_across_periodic_boundaries)

    # loop over possible bonding rules
    for br in bonding_rules
        # determine if the atom species correspond to the species in `bonding_rules`
        species_match = false
        if br.species_i == :* && br.species_j == :*
            species_match = true
        elseif br.species_i == :* && (species_i == br.species_j || species_j == br.species_j)
            species_match = true
        elseif br.species_j == :* && (species_i == br.species_i || species_j == br.species_j)
            species_match = true
        elseif (species_i == br.species_i && species_j == br.species_j) || (species_j == br.species_i && species_i == br.species_j)
            species_match = true
        end

        if species_match
            # determine if the atoms are close enough to bond
            if br.min_dist < r && br.max_dist > r
                return true
            else
                return false # found relevant bonding rule, don't apply others
            end
        end
    end
    return false # no bonding rule applied
end

"""
    remove_bonds!(crystal)

Remove all bonds from a crystal structure, `crystal::Crystal`.
"""
function remove_bonds!(crystal::Crystal)
    while ne(crystal.bonds) > 0
        rem_edge!(crystal.bonds, collect(edges(crystal.bonds))[1].src, collect(edges(crystal.bonds))[1].dst)
    end
end

"""
    infer_bonds!(crystal, include_bonds_across_periodic_boundaries, 
                    bonding_rules=[BondingRule(:H, :*, 0.4, 1.2), BondingRule(:*, :*, 0.4, 1.9)])

Populate the bonds in the crystal object based on the bonding rules. If a
pair doesn't have a suitable rule then they will not be considered bonded. 

`:*` is considered a wildcard and can be substituted for any species. It is a
good idea to include a bonding rule between two `:*` to allow any atoms to bond
as long as they are close enough.

The bonding rules are hierarchical, i.e. the first bonding rule takes precedence over the latter ones.

# Arguments
-`crystal::Crystal`: The crystal that bonds will be added to
-`include_bonds_across_periodic_boundaries::Bool`: Whether to check across the
    periodic boundary when calculating bonds
-`bonding_rules::Array{BondingRule, 1}`: The array of bonding rules that will
    be used to fill the bonding information. They are applied in the order that
    they appear.
"""
function infer_bonds!(crystal::Crystal, include_bonds_across_periodic_boundaries::Bool,
                      bonding_rules::Array{BondingRule, 1}=
                      [BondingRule(:H, :*, 0.4, 1.2), BondingRule(:*, :*, 0.4, 1.9)])
    @assert ne(crystal.bonds) == 0 @sprintf("The crystal %s already has bonds. Remove them with the `remove_bonds!` function before inferring new ones.", crystal.name)

    # loop over every atom
    for i in 1:crystal.atoms.n
        # loop over every unique pair of atoms
        for j in i+1:crystal.atoms.n
            if is_bonded(crystal, i, j, bonding_rules; include_bonds_across_periodic_boundaries=include_bonds_across_periodic_boundaries)
                add_edge!(crystal.bonds, i, j)
            end
        end
    end
end

"""
    sane_bonds = bond_sanity_check(crystal)

Run sanity checks on `crystal.bonds`.
* is the bond graph fully connected? i.e. does every vertex (=atom) in the bond graph have at least one edge?
* each hydrogen can have only one bond
* each carbon can have a maximum of four bonds

if sanity checks fail, refer to [`write_bond_information`](@ref) to write a .vtk to visualize the bonds.

Print warnings when sanity checks fail.
Return `true` if sanity checks pass, `false` otherwise.
"""
function bond_sanity_check(crystal::Crystal)
    sane_bonds = true
    for a = 1:crystal.atoms.n
        ns = neighbors(crystal.bonds, a)
        # is the graph fully connected?
        if length(ns) == 0
            @warn "atom $a = $(crystal.atoms.species[a]) in $(crystal.name) is not bonded to any other atom."
            sane_bonds = false
        end
        # does hydrogen have only one bond?
        if (crystal.atoms.species[a] == :H) && (length(ns) > 1)
            @warn "hydrogen atom $a in $(crystal.name) is bonded to more than one atom!"
            sane_bonds = false
        end
        # does carbon have greater than four bonds?
        if (crystal.atoms.species[a] == :C) && (length(ns) > 4)
            @warn "carbon atom $a in $(crystal.name) is bonded to more than four atoms!"
            sane_bonds = false
        end
    end
    return sane_bonds
end

# TODO remove? why is this needed?
"""
    bonds_equal = compare_bonds_in_crystal(crystal1, crystal2, atol=0.0)

Returns whether the bonds defined in crystal1 are the same as the bonds
defined in crystal2. It checks whether the atoms in the same positions
have the same bonds.

# Arguments
-`crystal1::Crystal`: The first crystal
-`crystal2::Crystal`: The second crystal
-`atol::Float64`: absolute tolerance for the comparison of coordinates in the crystal

# Returns
-`bonds_equal::Bool`: Wether the bonds in crystal1 and crystal2 are equal
"""
function compare_bonds_in_crystal(fi::Crystal, fj::Crystal; atol::Float64=0.0)
    if ne(fi.bonds) != ne(fj.bonds)
        return false
    end

    num_in_common = 0
    for edge_i in collect(edges(fi.bonds))
        for edge_j in collect(edges(fj.bonds))
            # either the bond matches going src-src dst-dst
            if  (fi.atoms.species[edge_i.src] == fj.atoms.species[edge_j.src] &&
                 fi.atoms.species[edge_i.dst] == fj.atoms.species[edge_j.dst] &&
                 isapprox(fi.atoms.xf[:, edge_i.src], fj.atoms.xf[:, edge_j.src]; atol=atol) &&
                 isapprox(fi.atoms.xf[:, edge_i.dst], fj.atoms.xf[:, edge_j.dst]; atol=atol)) ||
                # or the bond matches going src-dst dst-src
                (fi.atoms.species[edge_i.src] == fj.atoms.species[edge_j.dst] &&
                 fi.atoms.species[edge_i.dst] == fj.atoms.species[edge_j.src] &&
                 isapprox(fi.atoms.xf[:, edge_i.src], fj.atoms.xf[:, edge_j.dst]; atol=atol) &&
                 isapprox(fi.atoms.xf[:, edge_i.dst], fj.atoms.xf[:, edge_j.src]; atol=atol))
                num_in_common += 1
                break
            end
        end
    end
    return num_in_common == ne(fi.bonds) && num_in_common == ne(fj.bonds)
end
 
"""
    write_bond_information(crystal, filename)
    write_bond_information(crystal, center_at_origin=false)

Writes the bond information from a crystal to the selected filename.

# Arguments
-`crystal::Crystal`: The crystal to have its bonds written to a vtk file
-`filename::String`: The filename the bond information will be saved to. If left out, will default to crystal name.
- `center_at_origin::Bool`: center the coordinates at the origin of the crystal
"""
function write_bond_information(crystal::Crystal, filename::String; center_at_origin::Bool=false)
    if ne(crystal.bonds) == 0
        @warn("Crystal %s has no bonds present. To get bonding information for this crystal run `infer_bonds!` with an array of bonding rules\n", crystal.name)
    end
    if ! occursin(".vtk", filename)
        filename *= ".vtk"
    end

    vtk_file = open(filename, "w")

    @printf(vtk_file, "# vtk DataFile Version 2.0\n%s bond information\nASCII\nDATASET POLYDATA\nPOINTS %d double\n", crystal.name, nv(crystal.bonds))

    for i = 1:crystal.atoms.n
        if center_at_origin
            @printf(vtk_file, "%0.5f\t%0.5f\t%0.5f\n", (crystal.box.f_to_c * (crystal.atoms.coords.xf[:, i] - [0.5, 0.5, 0.5]))...)
        else
            @printf(vtk_file, "%0.5f\t%0.5f\t%0.5f\n", (crystal.box.f_to_c * crystal.atoms.coords.xf[:, i])...)
        end
    end
    @printf(vtk_file, "\nLINES %d %d\n", ne(crystal.bonds), 3 * ne(crystal.bonds))
    for edge in collect(edges(crystal.bonds))
        @printf(vtk_file, "2\t%d\t%d\n", edge.src - 1, edge.dst - 1)
    end
    close(vtk_file)
    @printf("Saving bond information for crystal %s to %s.\n", crystal.name, joinpath(pwd(), filename))
end

write_bond_information(crystal::Crystal; center_at_origin::Bool=false) = write_bond_information(crystal, split(crystal.name, ".")[1] * "_bonds.vtk", center_at_origin=center_at_origin)

# TODO remove bonds with atom i?
