include("muvt.jl")
include("nvt.jl")

using .MuVT, .NVT

import Base: +, /

# count proposed/accepted for each subtype
mutable struct MarkovCounts
    n_proposed::Array{Int, 1}
    n_accepted::Array{Int, 1}
end

###
#   collecting statistics
###
mutable struct GCMCstats
    n_samples::Int

    n::Array{Int, 1}
    n²::Array{Int, 1}

    U::SystemPotentialEnergy
    U²::SystemPotentialEnergy

    Un::Float64 # ⟨U n⟩
end

function GCMCstats(nb_species::Int)
    return GCMCstats(
        0,
        zeros(Int, nb_species),
        zeros(Int, nb_species),
        SystemPotentialEnergy(),
        SystemPotentialEnergy(),
        0.0
    )
end

function +(s1::GCMCstats, s2::GCMCstats)
    return GCMCstats(
        s1.n_samples + s2.n_samples,
        s1.n .+ s2.n,
        s1.n² .+ s2.n²,
        s1.U + s2.U,
        s1.U² + s2.U²,
        s1.Un + s2.Un
    )
end

function Base.sum(gcmc_stats::Array{GCMCstats, 1})
    nb_species = length(gcmc_stats[1].n)
    sum_stats = GCMCstats(nb_species)
    for gs in gcmc_stats
        sum_stats += gs
    end
    return sum_stats
end

function Base.print(gcmc_stats::GCMCstats)
    println("\t# samples: ", gcmc_stats.n_samples)
    println("\t⟨N⟩ (molecules) = ", gcmc_stats.n / gcmc_stats.n_samples)

    println("\t⟨U_gh, vdw⟩ (K) = ", gcmc_stats.U.gh.vdw / gcmc_stats.n_samples)
    println("\t⟨U_gh, Coulomb⟩ (K) = ", gcmc_stats.U.gh.es / gcmc_stats.n_samples)
    println("\t⟨U_gg, vdw⟩ (K) = ", gcmc_stats.U.gg.vdw / gcmc_stats.n_samples)
    println("\t⟨U_gg, Coulomb⟩ (K) = ", gcmc_stats.U.gg.es / gcmc_stats.n_samples)

    return println("\t⟨U⟩ (K) = ", sum(gcmc_stats.U) / gcmc_stats.n_samples)
end

# Compute average and standard error of the number of molecules and potential
# energy from an array of `GCMCstats`, each corresponding to statitics from an
# independent block/bin during the simulation. The average from each bin is
# treated as an independent sample and used to estimate the error in the
# estimate as a confidence interval.
function mean_stderr_n_U(gcmc_stats::Array{GCMCstats, 1})
    sum_gcmc_stats = sum(gcmc_stats)     # over entire simulation
    nb_blocks = length(gcmc_stats)       # number of blocks
    nb_species = length(gcmc_stats[1].n) # number of species

    # ⟨N⟩, ⟨U⟩ over entire simulation
    avg_n = sum_gcmc_stats.n / sum_gcmc_stats.n_samples
    avg_U = sum_gcmc_stats.U / sum_gcmc_stats.n_samples

    # ⟨N⟩, ⟨U⟩ over each block
    avg_n_blocks = [gs.n / gs.n_samples for gs in gcmc_stats]
    avg_U_blocks = [gs.U / gs.n_samples for gs in gcmc_stats]

    # std over the blocks
    std_n = std(avg_n_blocks)
    std_U = std(avg_U_blocks)

    # errors, computed treating each block average as an independent sample.
    err_n = 1.96 * std_n / sqrt(nb_blocks)
    err_U = 1.96 * std_U / sqrt(nb_blocks)

    return avg_n, err_n, avg_U, err_U
end

function print_results(results::Dict; print_title::Bool=true)
    if print_title
        # already print in GCMC tests...
        @printf(
            "GCMC simulation of %s in %s at %f K and %f bar pressure, %f bar fugacity using %s forcefield.\n\n",
            results["adsorbate"],
            results["xtal"],
            results["temperature (K)"],
            results["pressure (bar)"],
            results["fugacity (bar)"] / 100000.0,
            results["forcefield"]
        )
    end

    @printf(
        "\nUnit cell replication factors: %d %d %d\n\n",
        results["repfactors"][1],
        results["repfactors"][2],
        results["repfactors"][3]
    )
    # Markov stats
    println("")
    for key in ["# sample cycles", "# burn cycles", "# samples"]
        println(key * ": ", results[key])
    end

    for (_, proposal_description) in PROPOSAL_ENCODINGS
        total_proposals = results[@sprintf("Total # %s proposals", proposal_description)]
        fraction_accepted =
            results[@sprintf("Fraction of %s proposals accepted", proposal_description)]
        if total_proposals > 0
            printstyled(proposal_description; color=:yellow)
            @printf("\t%d total proposals.\n", total_proposals)
            @printf("\t%f %% proposals accepted.\n", 100.0 * fraction_accepted)
        end
    end

    println("")
    for s in 1:length(results["adsorbate"])
        println(results["adsorbate"][s])
        for key in ["⟨N⟩ (molecules)", "⟨N⟩ (molecules/unit cell)", "⟨N⟩ (mmol/g)"]
            @printf("%s: %f +/- %f\n", key, results[key][s], results["err " * key][s])
            if key == "⟨N⟩ (mmol/g)"
                println("")
            end
        end
    end

    for key in [
        "⟨U_gg, vdw⟩ (K)",
        "⟨U_gh, vdw⟩ (K)",
        "⟨U_gg, electro⟩ (K)",
        "⟨U_gh, electro⟩ (K)",
        "⟨U⟩ (K)"
    ]
        @printf("%s: %f +/- %f\n", key, results[key], results["err " * key])
        if key == "⟨N⟩ (mmol/g)"
            println("")
        end
    end

    for s in 1:length(results["adsorbate"])
        @printf(
            "\nQ_st (K) = %f = %f kJ/mol\n\n",
            results["Q_st (K)"][s],
            results["Q_st (K)"][s] * 8.314 / 1000.0
        )
    end
    return
end

function pretty_print(
    xtal::Crystal,
    molecule_templates::Array{Molecule{Cart}, 1},
    temperature::Float64,
    pressures::Array{Float64, 1},
    ljff::LJForceField
)
    printstyled("μVT simulation\n"; color=:yellow)

    print("crystal: ")
    printstyled(xtal.name; color=:green)

    print("\ntemperature: ")
    printstyled(@sprintf("%f K\n", temperature); color=:green)

    println("\nadsorbates:")
    for m in 1:length(molecule_templates)
        printstyled(molecule_templates[m].species; color=:yellow)
        print("\tpartial pressure = ")
        printstyled(@sprintf("%f bar\n", pressures[m]); color=:green)
    end

    println(" force field.")
    return printstyled(split(ljff.name, ".")[1]; color=:green)
end

"""
    results = stepwise_adsorption_isotherm(xtal, molecule_templates, temperature, pressures,
                                  ljff; n_sample_cycles=5000,
                                  n_burn_cycles=5000, sample_frequency=1,
                                  verbose=true, ewald_precision=1e-6, eos=:ideal,
                                  write_adsorbate_snapshots=false,
                                  snapshot_frequency=1, calculate_density_grid=false,
                                  density_grid_dx=1.0, density_grid_species=nothing,
                                  density_grid_sim_box::Bool=true,
                                  results_filename_comment="", show_progress_bar=false)

Run a set of grand-canonical (μVT) Monte Carlo simulations in series. Arguments are the
same as [`μVT_sim`](@ref), as this is the function run behind the scenes. An
exception is that we pass an array of pressures. The adsorption isotherm is computed step-
wise, where the ending configuration from the previous simulation (array of molecules) is
passed into the next simulation as a starting point. The ordering of `pressures` is
honored. By giving each simulation a good starting point, (if the next pressure does not
differ significantly from the previous pressure), we can reduce the number of burn cycles
required to reach equilibrium in the Monte Carlo simulation. Also see
[`adsorption_isotherm`](@ref) which runs the μVT simulation at each pressure in parallel.
"""
function stepwise_adsorption_isotherm(
    xtal::Crystal,
    molecule_templates::Array{Molecule{Cart}, 1},
    temperature::Float64,
    pressures::Array{Array{Float64, 1}, 1},
    ljff::LJForceField;
    kwargs...
)

    # simulation only works if xtal is in P1
    assert_P1_symmetry(xtal)

    # we only consider a pure gas isotherm
    if length(molecule_templates) != 1
        error("We only calculate pure gas isotherms")
    end

    results = Dict{String, Any}[] # push results to this array
    molecules = nothing # initialize with empty xtal
    for pressure in pressures
        result, molecules = μVT_sim(
            xtal,
            molecule_templates,
            temperature,
            pressure,
            ljff;
            molecules=molecules, # essential step here
            kwargs...
        )
        push!(results, result)
    end

    return results
end

# overload function for backward compatibility
function stepwise_adsorption_isotherm(
    xtal::Crystal,
    molecule_templates::Molecule{Cart},
    temperature::Float64,
    pressures::Array{Float64, 1},
    ljff::LJForceField;
    kwargs...
)
    return stepwise_adsorption_isotherm(
        xtal,
        [molecule_templates],
        temperature,
        [[p] for p in pressures],
        ljff;
        kwargs...
    )
end

"""
    results = adsorption_isotherm(xtal, molecule_templates, temperature, 
                                  partial_pressures,
                                  ljff; n_sample_cycles=5000,
                                  n_burn_cycles=5000, sample_frequency=1,
                                  verbose=true, ewald_precision=1e-6, eos=:ideal,
                                  write_adsorbate_snapshots=false,
                                  snapshot_frequency=1, calculate_density_grid=false,
                                  density_grid_dx=1.0, density_grid_species=nothing,
                                  density_grid_sim_box=true,
                                  results_filename_comment="", show_progress_bar=false)

Run a set of grand-canonical (μVT) Monte Carlo simulations in parallel. Arguments are the
same as [`μVT_sim`](@ref), as this is the function run in parallel behind the scenes.
The only exception is that we pass an array of pressures, and we only consider a single species.
To give Julia access to multiple cores, run your script as `julia -p 4 mysim.jl` to allocate e.g. four cores. See
[Parallel Computing](https://docs.julialang.org/en/stable/manual/parallel-computing/#Parallel-Computing-1).
"""
function adsorption_isotherm(
    xtal::Crystal,
    molecule_templates::Array{Molecule{Cart}, 1},
    temperature::Float64,
    pressures::Array{Array{Float64, 1}, 1},
    ljff::LJForceField;
    kwargs...
)

    # simulation only works if xtal is in P1
    assert_P1_symmetry(xtal)

    # we only consider a pure gas isotherm
    if length(molecule_templates) != 1
        error("We only calculate pure gas isotherms")
    end

    # make a function of pressure only to facilitate uses of `pmap`
    function run_pressure(partial_pressures::Array{Float64, 1})
        return μVT_sim(
            xtal,
            molecule_templates,
            temperature,
            partial_pressures,
            ljff;
            kwargs...
        )[1] # only return results
    end # only return results

    # for load balancing, larger pressures with longer computation time goes first
    ids = sortperm(pressures; rev=true)

    # run gcmc simulations in parallel using Julia's pmap parallel computing function
    results = pmap(run_pressure, pressures[ids])

    # return results in same order as original pressure even though we permuted them for
    #  better load balancing.
    return results[[findall(x -> x == i, ids)[1] for i in eachindex(ids)]]
end

# overload function for backward compatibility
function adsorption_isotherm(
    xtal::Crystal,
    molecule_templates::Molecule{Cart},
    temperature::Float64,
    pressures::Array{Float64, 1},
    ljff::LJForceField;
    kwargs...
)
    return adsorption_isotherm(
        xtal,
        [molecule_templates],
        temperature,
        [[p] for p in pressures],
        ljff;
        kwargs...
    )
end

"""
    dataframe = isotherm_sim_results_to_dataframe(desired_props, xtal,
                                                  molecule, temperature,
                                                  pressures, ljff,
                                                  n_burn_cycles, n_sample_cycles;
                                                  where_are_jld_files=nothing)`

convert the `.jld2` results output files auto-saved from [`μVT_sim`](@ref) into a `DataFrame`.
each row of the `DataFrame` corresponds to a pressure in the adsorption isotherm.
`desired_props` is an array of desired properties from the results.
to locate the requested output files, this function calls [`μVT_output_filename`](@ref) to retrieve the file names.

# Arguments

  - `desired_props::Array{String, 1}`: An array containing names of properties to be retrieved from the `.jld2` file
  - `xtal::Crystal`: The porous crystal
  - `molecule::Array{Molecule{Cart}, 1}`: The adsorbate molecule
  - `temperature::Float64`: The temperature in the simulation, units: Kelvin (K)
  - `pressures::Array{Array{Float64, 1}, 1}`: The pressures in the adsorption isotherm simulation(s), units: bar
  - `ljff::LJForceField`: The molecular model being used in the simulation to describe intermolecular Van
    der Waals interactions
  - `n_burn_cycles::Int64`: The number of burn cycles used in this simulation
  - `n_sample_cycles::Int64`: The number of sample cycles used in the simulations
  - `where_are_jld_files::Union{String, Nothing}=nothing`: The location to the simulation files. defaults to
    `PorousMaterials.rc[:paths][:simulations]`.
  - `comment::String=""`: comment appended to outputfilename

# Returns

  - `dataframe::DataFrame`: A `DataFrame` containing the simulated data along the adsorption isotherm, whose
    columns are for the specified properties

# Note

A range of pressures can be used to select a batch of simulation files to be included in the `DataFrame`.

# Example

```julia
xtal = Crystal("SBMOF-1.cif")
molecule = Molecule("Xe")
ljff = LJForceField("UFF"; mixing_rules="Lorentz-Berthelot")
temperature = 298.0 # K
pressures = 10 .^ range(-2; stop=2, length=15)

dataframe = isotherm_sim_results_to_dataframe(
    ["pressure (bar)", "⟨N⟩ (mmol/g)"],
    xtal,
    molecule,
    temperature,
    pressures,
    ljff,
    10000,
    10000
)
dataframe[Symbol("pressure (bar)")] # pressures
dataframe[Symbol("⟨N⟩ (mmol/g)")] # adsorption at corresponding pressures
```
"""
function isotherm_sim_results_to_dataframe(
    desired_props::Array{String, 1},
    xtal::Crystal,
    molecule_templates::Array{Molecule{Cart}, 1},
    temperature::Float64,
    pressures::Array{Array{Float64, 1}, 1},
    ljff::LJForceField,
    n_burn_cycles::Int64,
    n_sample_cycles::Int64;
    comment::String="",
    where_are_jld_files::Union{String, Nothing}=nothing
)
    # determine the location of the data files
    if isnothing(where_are_jld_files)
        where_are_jld_files = rc[:paths][:simulations]
    end
    # prepare dataframe to populate
    df = DataFrame()
    # loop over pressures and populate dataframe
    for (i, pressure) in enumerate(pressures)
        jld2_filename = μVT_output_filename(
            xtal,
            molecule_templates,
            temperature,
            pressure,
            ljff,
            n_burn_cycles,
            n_sample_cycles;
            comment=comment
        )
        # load in the results as a dictionary
        @load joinpath(where_are_jld_files, jld2_filename) results

        # population column names, taking into account types
        if i == 1
            for col in desired_props
                insertcols!(
                    df,
                    length(names(df)) + 1,
                    Symbol(col) => typeof(results[col])[]
                )
            end
        end
        push!(df, [results[prop] for prop in desired_props])
    end
    return df
end

# overload function for backward compatibility
function isotherm_sim_results_to_dataframe(
    desired_props::Array{String, 1},
    xtal::Crystal,
    molecule_templates::Molecule{Cart},
    temperature::Float64,
    pressures::Array{Float64, 1},
    ljff::LJForceField,
    n_burn_cycles::Int64,
    n_sample_cycles::Int64;
    kwargs...
)
    return isotherm_sim_results_to_dataframe(
        desired_props,
        xtal,
        [molecule_templates],
        temperature,
        [[p] for p in pressures],
        ljff,
        n_burn_cycles,
        n_sample_cycles;
        kwargs...
    )
end
