# Details from http://www.stochasticlifestyle.com/finalizing-julia-package-documentation-testing-coverage-publishing/

function runtest(testfile::String)
    @info "Testing $(testfile)"
    include(testfile)
end

testfiles = ["box.jl",
             "matter.jl",
             "crystal.jl",
             "distance.jl",
             "misc.jl",
             "forcefield.jl",
             "molecule.jl",
             "vdw_energetics.jl",
             "energy_utilities.jl",
             "electrostatics.jl",
             "mc_helpers.jl",
             "eos.jl",
             "assert_p1_symmetry.jl",
             "grid.jl",
             #"paths.jl" ## NOTE run(`cp ...`) incompatible with Windows
             ]

[runtest(testfile) for testfile in testfiles]
