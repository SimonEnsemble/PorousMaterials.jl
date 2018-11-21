# Details from http://www.stochasticlifestyle.com/finalizing-julia-package-documentation-testing-coverage-publishing/
# Start Test Script

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
include("henry_checkpoint_test.jl")
include("grid_test.jl")
include("path_test.jl")
include("generic_rotation_test.jl")
