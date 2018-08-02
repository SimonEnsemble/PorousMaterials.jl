#!/bin/bash
cd test
if [ julia runtests.jl ];
then
    cd ..
    exit 0
else
    cd ..
    exit 1
fi
