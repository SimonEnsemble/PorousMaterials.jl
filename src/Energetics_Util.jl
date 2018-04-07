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
    outside_box = completely_outside_box(molecule::Molecule, box::Box)

returns true if the center of mass of a molecule is outside of the box.
"""
function outside_box(molecule::Molecule, box::Box)
    xf = box.c_to_f * molecule.center_of_mass
    for k = 1:3
        if (xf[k] > 1.0) | (xf[k] < 0.0)
            return true
        end
    end
    return false
end
