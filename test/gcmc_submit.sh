#!/bin/bash

# use current working directory for input and output
# default is to use the users home directory
#$ -cwd

# name this job
#$ -N gcmc_simulation

# send stdout and stderror to this file
#$ -o gcmc_simulation.o
#$ -e gcmc_simulation.e
#$ -j y

#the list of users who will recieve mail about this job
#$ -M yourONID@oregonstate.edu
#options for when mail is sent out, this will send mail when the job begins,
#       ends, or is aborted
#$ -m bea

# select queue - if needed; mime5 is SimonEnsemble priority queue but is restrictive.
#$ -q mime5

#set up a parallel environment
#$ -pe openmpi 4

# print date and time
date

julia -p 4 gcmc_test.jl 
