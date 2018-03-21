"""
    nearest_image!(fractional_distance,repfactors)

runs the nearest image convention given an array of fractional coordinates and
the repfactors of the supercell. This code was pulled from Arni's vdw_energy
function in his Energetics module

See comments in vdw_energy for more description
"""
function nearest_image!(fractional_distance::Array{Float64}, repfactors::Tuple{Int64, Int64, Int64})
    for coords = 1:3
        if abs(fractional_distance[coords]) > repfactors[coords] / 2.0
            fractional_distance[coords] -= sign(fractional_distance[coords]) * repfactors[coords]
        end
    end
end

"""
    outside_box = completely_outside_box(molecule, box)

returns true if each atom of a given molecule is completely outside of a given box and false otherwise
"""
function completely_outside_box(molecule::Molecule, box::Box)
    xf = box.c_to_f * molecule.x
    for xyz = 1:3 # loop over x, y, z coordinate
        # if none of the coords are less than 1 it must be outside of the box
        if sum(xf[xyz, :] .<= 1.0) == 0
            return true
        # if none of the coords are greater than 0 it must be outside of the box
        elseif sum(xf[xyz, :] .>= 0.0) == 0
            return true
        end
    end
    return false
end
