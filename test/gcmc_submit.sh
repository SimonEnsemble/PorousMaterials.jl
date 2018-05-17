#!/bin/bash

# use current working directory for input and output
# default is to use the users home directory
#$ -cwd

# name this job
#$ -N GCMC-tester

# send stdout and stderror to this file
#$ -o gcmc_results
#$ -e gcmc_errors
#$ -j y

#the list of users who will recieve mail about this job
#$ -M yorkar@oregonstate.edu
#options for when mail is sent out, this will send mail when the job begins,
#       ends, or is aborted
#$ -m bea

# select queue - if needed 
#$ -q mime5,share,share2,share3,share4

#set up a parallel environment
#$ -pe OpenMP 4 

# print date and time
date

julia -p 4 gcmc_test.jl 
