#!/bin/env bash

julia -p 24 -e '
    using Distributed
    @everywhere begin
        import Pkg
        Pkg.activate(".")
    end
    #Pkg.test("PorousMaterials")
    @info "Running test/gcmc_long.jl on $(length(workers())) workers."
    @warn "This will take HOURS!"
    cd("test")
    include("gcmc_long.jl")
    @info "Long GCMC tests complete!"
'
