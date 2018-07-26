using PorousMaterials
using JLD

jobs_to_run = Dict("zif71_bogus_charges" => Dict("gases" => ["CO2EPM2"],
                                                 "temperatures" => [298.0],
                                                 "pressures" => vcat(collect(linspace(0.0, 1.0, 11))[2:end],
                                                                     collect(linspace(1.0, 20.0, 20))[2:end]),
                                                 # this applies to all pressures/temperatures
                                                 "forcefield" => LJForceField("Greg_bogus_ZIF71.csv", cutoffradius=12.8),
                                                 "n_burn_cycles" => 5,
                                                 "n_sample_cycles" => 5
                                                ),
                  )
