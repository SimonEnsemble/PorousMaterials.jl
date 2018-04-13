"""
    nearest_image!(dxf, repfactors)

Applies the nearest image convention on a vector `dxf` between two atoms
in fractional space; modifies `dxf` for nearest image convention.

See comments in vdw_energy for more description.
"""
function nearest_image!(dxf::Array{Float64, 1}, repfactors::Tuple{Int, Int, Int})
    for k = 1:3 # loop over components
        if abs(dxf[k]) > repfactors[k] / 2.0
            dxf[k] -= sign(dxf[k]) * repfactors[k]
        end
    end
end

function nearest_image!(dxf::Array{Float64, 2}, repfactors::Tuple{Int, Int, Int})
    for a = 1:size(dxf)[2] # loop over atoms
        for k = 1:3 # loop over components
            if abs(dxf[k, a]) > repfactors[k] / 2.0
                dxf[k, a] -= sign(dxf[k, a]) * repfactors[k]
            end
        end
    end
end
