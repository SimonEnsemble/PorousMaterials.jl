# Details from http://www.stochasticlifestyle.com/finalizing-julia-package-documentation-testing-coverage-publishing/

function runtest(testfile::String)
    @info "Testing $(testfile)"
    try
        include(testfile)
    catch exception
        @error "Exception in $(testfile)" exception
    end
end

testfiles = [
    "crystal.jl",
    "box.jl",
    "matter.jl",
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
    "gcmc_quick.jl",
    "paths.jl" # must be last so as to not interfere with reading in files
    ]

runtest.(testfiles)
