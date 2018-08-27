#TODO split these into different files for different tests so e.g. framework in one suite is not used in another suite so it's clear what's goign on.

# Details from http://www.stochasticlifestyle.com/finalizing-julia-package-documentation-testing-coverage-publishing/
# Start Test Script

include("grid_test.jl")

include("box_test.jl")

include("crystal_test.jl")

include("molecule_test.jl")

include("nearest_image_test.jl")

include("forcefield_test.jl")

include("vdw_energetics_test.jl")

include("energetics_util_test.jl")

include("electrostatics_energetics_test.jl")

include("mchelpers_test.jl")

include("guest_guest_energetics_test.jl")

include("eos_test.jl")

include("gcmc_checkpoints_test.jl")
