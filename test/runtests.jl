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

using Test, FIGlet

FIGlet.render("PorousMaterials.jl", FIGlet.availablefonts()[6])

for testfile ∈ testfiles
    @info "Running test/$testfile"
    include(testfile)
end

@info "Tests complete!"
