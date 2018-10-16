
<a id='Henry-Coefficient-Calculations-1'></a>

## Henry Coefficient Calculations

<a id='PorousMaterials.henry_coefficient' href='#PorousMaterials.henry_coefficient'>#</a>
**`PorousMaterials.henry_coefficient`** &mdash; *Function*.



```
result = henry_coefficient(framework, molecule, temperature, ljforcefield,
                            nb_insertions=1e6, verbose=true, ewald_precision=1e-6,
                            autosave=true)
```

Conduct particle insertions to compute the Henry coefficient Kₕ of a molecule in a framework. Also, for free, the heat of adsorption and ensemble average energy of adsorption is computed. The Henry coefficient is a model for adsorption at infinite dilution (low coverage): ⟨N⟩ = Kₕ P, where P is pressure and Kₕ is the Henry coefficient.

Kₕ = β ⟨e^{-β U}⟩, where the average is over positions and orientations of the molecule in the framework.

**Arguments**

  * `framework::Framework`: the porous crystal in which we seek to simulate adsorption
  * `molecule::Molecule`: the adsorbate molecule
  * `temperature::Float64`: temperature of bulk gas phase in equilibrium with adsorbed phase   in the porous material. units: Kelvin (K)
  * `ljforcefield::LJForceField`: the molecular model used to describe the   energetics of the adsorbate-adsorbate and adsorbate-host van der Waals interactions.
  * `insertions_per_volume::Int`: number of Widom insertions to perform for computing the

average, per unit cell volume (Å³)

  * `verbose::Bool`: whether or not to print off information during the simulation.
  * `ewald_precision::Float64`: desired precision for Ewald summations; used to determine

the replication factors in reciprocal space.

  * `autosave::Bool`: save results file as a .jld in PATH*TO*DATA * `sims`
  * `filename_comment::AbstractString`: An optional comment that will be appended to the name of the saved file.
  * `write_checkpoint::Bool`: Will periodically save checkpoints to start the job from a previous state.
  * `load_checkpoint::Bool`: Instructs the program to look for a checkpoint file in `data/henry_checkpoints`

and start the simulation from that point.

  * `checkpoint_frequency::Int`: The frequency at which we will save a checkpoint file. Is only used if `write_checkpoint=true`

**Returns**

  * `result::Dict{String, Float64}`: A dictionary containing all the results from the Henry coefficient simulation


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/Henry.jl#L4-L41' class='documenter-source'>source</a><br>

<a id='PorousMaterials.henry_result_savename' href='#PorousMaterials.henry_result_savename'>#</a>
**`PorousMaterials.henry_result_savename`** &mdash; *Function*.



```
save_name = henry_result_savename(framework, molecule, temperature,
                               ljforcefield, insertions_per_volume;
                               comment="")
```

Determine the name of files saved while calculating the henry coefficient. It uses many pieces of information from the simulation to ensure the file name accurately describes what it holds.

**Arguments**

  * `framework::Framework`: The porous crystal being tested
  * `molecule::Molecule`: The molecule being tested inside the crystal
  * `temperature::Float64`: The temperature used in the simulation units: Kelvin (K)
  * `ljforcefield::LJForceField`: The molecular model being used in the simulation   to describe the intermolecular Van der Waals forces
  * `insertions_per_volume::Union{Int, Float64}`: The number of widom insertions per unit volume.   Will be scaled according to the framework we're working with
  * `comment::AbstractString`: An optional comment that will be appended to the filename


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/Henry.jl#L327-L345' class='documenter-source'>source</a><br>


<a id='Grand-Canonical-Monte-Carlo-Simulations-1'></a>

## Grand-Canonical Monte Carlo Simulations

<a id='PorousMaterials.gcmc_simulation' href='#PorousMaterials.gcmc_simulation'>#</a>
**`PorousMaterials.gcmc_simulation`** &mdash; *Function*.



```
results, molecules = gcmc_simulation(framework, molecule, temperature, pressure,
                                     ljforcefield; n_sample_cycles=5000,
                                     n_burn_cycles=5000, sample_frequency=1,
                                     verbose=false, molecules=Molecule[],
                                     eos=:ideal, ewald_precision=1e-6,
                                     load_checkpoint_file=false,
                                     show_progress_bar=false, checkpoint=Dict(),
                                     write_checkpoints=false, checkpoint_frequency=100,
                                     filename_comment="")
```

Runs a grand-canonical (μVT) Monte Carlo simulation of the adsorption of a molecule in a framework at a particular temperature and pressure using a Lennard Jones force field.

A cycle is defined as max(20, number of adsorbates currently in the system) Markov chain proposals. Current Markov chain moves implemented are particle insertion/deletion and translation.

**Arguments**

  * `framework::Framework`: the porous crystal in which we seek to simulate adsorption
  * `molecule::Molecule`: a template of the adsorbate molecule of which we seek to simulate
  * `temperature::Float64`: temperature of bulk gas phase in equilibrium with adsorbed phase   in the porous material. units: Kelvin (K)
  * `pressure::Float64`: pressure of bulk gas phase in equilibrium with adsorbed phase in the   porous material. units: bar   the adsorption
  * `ljforcefield::LJForceField`: the molecular model used to describe the   energetics of the adsorbate-adsorbate and adsorbate-host van der Waals interactions.
  * `n_burn_cycles::Int`: number of cycles to allow the system to reach equilibrium before   sampling.
  * `n_sample_cycles::Int`: number of cycles used for sampling
  * `sample_frequency::Int`: during the sampling cycles, sample e.g. the number of adsorbed   gas molecules every this number of Markov proposals.
  * `verbose::Bool`: whether or not to print off information during the simulation.
  * `molecules::Array{Molecule, 1}`: a starting configuration of molecules in the framework.

Note that we assume these coordinates are Cartesian, i.e. corresponding to a unit box.

  * `ewald_precision::Float64`: The desired precision for the long range Ewald summation
  * `eos::Symbol`: equation of state to use for calculation of fugacity from pressure. Default

is ideal gas, where fugacity = pressure.

  * `load_checkpoint_file::Bool`: Will find a checkpoint file corresponding to the [`gcmc_result_savename`](docs.md#PorousMaterials.gcmc_result_savename) if true.   If that file is not found, function will throw an error.
  * `checkpoint::Dict`: A checkpoint dictionary that will work as a starting configuration for the run.   The dictionary has to have the following keys: `outer_cycle`, `molecules`, `system_energy`, `current_block`, `gcmc_stats`, `markov_counts`, `markov_chain_time` and `time`. If this argument is used, keep `load_checkpoint_file=false`.
  * `write_checkpoints::Bool`: Will save checkpoints in data/gcmc_checkpoints if this is true.
  * `checkpoint_frequency::Int`: Will save checkpoint files every `checkpoint_frequency` cycles.
  * `filename_comment::AbstractString`: An optional comment that will be appended to the name of the saved file (if autosaved)


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/GCMC.jl#L219-L266' class='documenter-source'>source</a><br>

<a id='PorousMaterials.adsorption_isotherm' href='#PorousMaterials.adsorption_isotherm'>#</a>
**`PorousMaterials.adsorption_isotherm`** &mdash; *Function*.



```
results = adsorption_isotherm(framework, molecule, temperature, pressures,
                              ljforcefield; n_sample_cycles=5000,
                              n_burn_cycles=5000, sample_frequency=1,
                              verbose=true, ewald_precision=1e-6, eos=:ideal, 
                              load_checkpoint_file=false, checkpoint=Dict(), 
                              write_checkpoints=false, checkpoint_frequency=50,
                              filename_comment="", show_progress_bar=false)
```

Run a set of grand-canonical (μVT) Monte Carlo simulations in parallel. Arguments are the same as [`gcmc_simulation`](docs.md#PorousMaterials.gcmc_simulation), as this is the function run in parallel behind the scenes. The only exception is that we pass an array of pressures. To give Julia access to multiple cores, run your script as `julia -p 4 mysim.jl` to allocate e.g. four cores. See [Parallel Computing](https://docs.julialang.org/en/stable/manual/parallel-computing/#Parallel-Computing-1).


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/GCMC.jl#L172-L186' class='documenter-source'>source</a><br>

<a id='PorousMaterials.stepwise_adsorption_isotherm' href='#PorousMaterials.stepwise_adsorption_isotherm'>#</a>
**`PorousMaterials.stepwise_adsorption_isotherm`** &mdash; *Function*.



```
results = stepwise_adsorption_isotherm(framework, molecule, temperature, pressures,
                              ljforcefield; n_sample_cycles=5000,
                              n_burn_cycles=5000, sample_frequency=1,
                              verbose=true, ewald_precision=1e-6, eos=:ideal,
                              load_checkpoint_file=false, checkpoint=Dict(),
                              write_checkpoints=false, checkpoint_frequency=50,
                              filename_comment="", show_progress_bar=false)
```

Run a set of grand-canonical (μVT) Monte Carlo simulations in series. Arguments are the same as [`gcmc_simulation`](docs.md#PorousMaterials.gcmc_simulation), as this is the function run behind the scenes. An exception is that we pass an array of pressures. The adsorption isotherm is computed step- wise, where the ending configuration from the previous simulation (array of molecules) is passed into the next simulation as a starting point. The ordering of `pressures` is honored. By giving each simulation a good starting point, (if the next pressure does not differ significantly from the previous pressure), we can reduce the number of burn cycles required to reach equilibrium in the Monte Carlo simulation. Also see [`adsorption_isotherm`](docs.md#PorousMaterials.adsorption_isotherm) which runs the μVT simulation at each pressure in parallel.


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/GCMC.jl#L125-L143' class='documenter-source'>source</a><br>

<a id='PorousMaterials.gcmc_result_savename' href='#PorousMaterials.gcmc_result_savename'>#</a>
**`PorousMaterials.gcmc_result_savename`** &mdash; *Function*.



```
file_save_name = gcmc_result_savename(framework_name, molecule_species,
                                    ljforcefield_name, temperature, pressure,
                                    n_burn_cycles, n_sample_cycles; comment="")
```

Determine the name of files saved during the GCMC simulation, be molecule positions or results. It uses many pieces of information from the simulation to ensure the file name accurately describes what it holds.

**Arguments**

  * `framework_name::AbstractString`: The porous crystal being tested
  * `molecule_species::Symbol`: The molecule being tested inside the porous crystal
  * `ljforcefield_name::AbstractString`: The molecular model being used in this   simulation to describe intermolecular Van der Waals interactions
  * `temperature::Float64`: The temperature used in the simulation units: Kelvin (K)
  * `pressure::Float64`: The pressure used in the simulation units: bar
  * `n_burn_cycles::Int`: The number of burn cycles used in this simulation
  * `n_sample_cycles::Int`: The number of sample cycles used in this simulation
  * `comment::AbstractString`: An optional comment that will be appended to the end of the filename


<a target='_blank' href='https://github.com/SimonEnsemble/PorousMaterials.jl/blob/0b314a309d738169927abfd7afcb30c6f7d7a651/src/GCMC.jl#L761-L780' class='documenter-source'>source</a><br>

