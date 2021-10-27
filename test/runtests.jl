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

using Test, FIGlet, Documenter

font_num = 579
FIGlet.render("Porous", FIGlet.availablefonts()[font_num])
FIGlet.render("Materials", FIGlet.availablefonts()[font_num])

using PorousMaterials

for testfile ∈ testfiles
    @info "Running test/$testfile"
    include(testfile)
end

doctest(PorousMaterials)

@info "Tests complete!"
