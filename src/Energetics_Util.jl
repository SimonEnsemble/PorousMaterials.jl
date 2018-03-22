"""
    nearest_image!(fractional_distance, repfactors)

applies the nearest image convention on a vector `fractional_distance` between two atoms
in fractional space; modifies `fractional_distance` for nearest image convention.

See comments in vdw_energy for more description
"""
function nearest_image!(fractional_distance::Array{Float64}, repfactors::Tuple{Int64, Int64, Int64})
    for xyz = 1:3 # xf, yf, or zf coordinate
        if abs(fractional_distance[xyz]) > repfactors[xyz] / 2.0
            fractional_distance[xyz] -= sign(fractional_distance[xyz]) * repfactors[xyz]
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
