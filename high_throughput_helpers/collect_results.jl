using PorousMaterials
using JLD
 # using PyPlot

include("jobs_to_run.jl")

for (structure, conditions) in jobs_to_run
    ljff = jobs_to_run[structure]["forcefield"]
    n_burn_cycles = jobs_to_run[structure]["n_burn_cycles"]
    n_sample_cycles = jobs_to_run[structure]["n_sample_cycles"]

    for gas in conditions["gases"]
        for temperature in conditions["temperatures"]
            # store results per adsorption isotherm
            results = []
            for pressure in conditions["pressures"]

                save_results_filename = joinpath(PorousMaterials.PATH_TO_DATA, "gcmc_sims", PorousMaterials.root_save_filename(
                    structure, Symbol(gas), ljff.name, temperature, pressure, 
                    n_burn_cycles, n_sample_cycles) * ".jld")

                result = JLD.load(save_results_filename)["results"]
                push!(results, result)
            end
        end
    end
end
