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

using Test, Documenter, PorousMaterials

PorousMaterials.banner()

for testfile in testfiles
    @info "Running test/$testfile"
    @time include(testfile)
end

doctest(PorousMaterials)

@info "Tests complete!"
