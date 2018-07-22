"""
Write a job submission script to submit for a μVT simulation of `gas` in `structurename`
at `temperature` K and `pressure` bar.
"""
function write_gcmc_submit_script(structurename::AbstractString, gas::AbstractString, 
                                  temperature::Float64, pressure::Float64;
                                  email_address::AbstractString="")
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
                                                                                                
    #the list of users who will recieve mail about this job                                     
    #\$ -M %s                                                           
    # options for when mail is sent out, this will send mail when the job begins,                
    #       ends, or is aborted                                                                 
    #\$ -m bea                                                                                   
                                                                                                
    ## select queue - if needed; mime5 is SimonEnsemble priority queue but is restrictive.       
    ##\$ -q mime5                                                                                 
                                                                                                
    # print date and time
    date

    julia run_gcmc_simulation.jl %s %s %f %f
    """,
	structurename, 
    jobscriptdir * "/" * structurename * "_" * string(temperature) * "_" * string(pressure), 
    jobscriptdir * "/" * structurename * "_" * string(temperature) * "_" * string(pressure),
    email_address,
    structurename, gas, temperature, pressure)

    close(gcmc_submit)
end

# import list of materials
include("jobs_to_run.jl")

for (structurename, conditions) in keys(jobs_to_run)
    for gasname in conditions["gases"]
        for temperature in conditions["temperatures"]
            for pressure in conditions["pressures"]
                write_gcmc_submit_script(structurename, gas, temperature, pressure)
                run(`qsub gcmc_submit.sh`)
                sleep(1)
            end
        end
    end
end
		
