```@meta
DocTestSetup = quote
  using PorousMaterials
end
```

# Data paths

`PorousMaterials.jl` is built on [`Xtals.jl`](https://SimonEnsemble.github.io/Xtals.jl/dev).  All bindings are re-exported.  As such, the global
variable dictionary `rc` is part of the `PorousMaterials` namespace.  This is where data such as atomic masses and covalent radii are stored.

```jldoctest
keys(rc)
# output
KeySet for a Dict{Symbol, Any} with 7 entries. Keys:
  :covalent_radii
  :bonding_rules
  :scipy
  :atomic_masses
  :cpk_colors
  :paths
  :pymatgen
```

`rc[:paths]` gives a listing of the read/write paths for various kinds of information.  By default, `rc[:paths][:data]` will be set to `./data` based on
the present working directory when the module is imported.  All the other paths are set relative to `rc[:paths][:data]`:

```jldoctest
keys(rc[:paths])
# output
KeySet for a Dict{Symbol, String} with 6 entries. Keys:
  :forcefields
  :grids
  :molecules
  :crystals
  :data
  :simulations
```

So, for example, if `PorousMaterials` is loaded in `/my_project`, then `rc[:paths][:data]` will by default be `/my_project/data` and the other paths will be subdirectories, e.g. `/my_project/data/simulations`.

To move the root of the data tree and change all the paths together, use `set_paths`:

```julia
set_paths("other_project")
# now rc[:paths][:data] is "other_project/data"
# rc[:paths][:crystals] is "other_project/data/crystals"
# etc.
```
