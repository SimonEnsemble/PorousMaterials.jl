testfiles = [
    "molecule.jl",
    "gcmc_quick.jl",
    "grid.jl",
    "misc.jl",
    "forcefield.jl",
    "vdw_energetics.jl",
    "energy_utilities.jl",
    "electrostatics.jl",
    "mc_helpers.jl",
    "eos.jl"
    ]

@info "\n\n\tPorousMaterials\n\n\n"

for testfile ∈ testfiles
    @info "Running test/$testfile"
    include(testfile)
end

@info "Done!"
