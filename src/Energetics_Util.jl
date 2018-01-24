module Energetics_Util

export nearest_image!

"""
    nearest_image!(fractional_distance,repfactors)

runs the nearest image convention given an array of fractional coordinates and
the repfactors of the supercell. This code was pulled from Arni's vdw_energy
function in his Energetics module

See comments in vdw_energy for more description
"""
#TODO put in Energetics_Utils.jl
function nearest_image!(fractional_distance::Array{Float64}, repfactors::Tuple{Int64, Int64, Int64})
    for coords = 1:3
        if abs(fractional_distance[coords]) > repfactors[coords] / 2
            fractional_distance[coords] -= sign(fractional_distance[coords]) * repfactors[coords]
        end #if statement
    end #for loop for going over
end #nearest_image

end #module
