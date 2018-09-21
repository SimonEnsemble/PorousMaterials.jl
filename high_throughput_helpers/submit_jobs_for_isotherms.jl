"""
Write a job submission script to submit for a μVT simulation of `gas` in `structurename`
at `temperature` K and `pressure` bar.
"""
function write_gcmc_submit_script(structurename::AbstractString, gas::AbstractString, 
                                  temperature::Float64, pressure::Float64)
    jobscriptdir = "jobz"
    if ~ isdir(jobscriptdir)
        mkdir(jobscriptdir)
    end

    @printf("Writing submission script for μVT sim of %s in %s at %f bar and %f K.\n",
        gas, structurename, temperature, pressure)

    ## build gcmc_submit.sh script
    gcmc_submit = open("gcmc_submit.sh", "w")
    @printf(gcmc_submit,
    """
    #!/bin/bash                                                                                 
                                                                                                
    # use current working directory for input and output                                        
    # default is to use the users home directory                                                
    #\$ -cwd                                                                                     
                                                                                                
    # name this job                                                                             
    #\$ -N %s
                                                                                                
    # send stdout and stderror to this file                                                     
    #\$ -o %s.o                                                                          
    #\$ -e %s.e                                                                          
                                                                                                
    ## select queue - if needed; mime5 is SimonEnsemble priority queue but is restrictive.       
    ##\$ -q mime5                                                                                 
                                                                                                
    # print date and time
    date

    julia run_gcmc_simulation.jl %s %s %f %f
    """,
	structurename, 
    joinpath(jobscriptdir, structurename * "_" * string(temperature) * "_" * string(pressure)), 
    joinpath(jobscriptdir, structurename * "_" * string(temperature) * "_" * string(pressure)),
    structurename, gas, temperature, pressure)

    close(gcmc_submit)
end

# import list of materials
include("jobs_to_run.jl")

for (structure, conditions) in jobs_to_run
    for gas in conditions["gases"]
        for temperature in conditions["temperatures"]
            for pressure in conditions["pressures"]
                write_gcmc_submit_script(structure, gas, temperature, pressure)
                run(`qsub gcmc_submit.sh`)
                sleep(1)
                run(`rm gcmc_submit.sh`)
            end
        end
    end
end
