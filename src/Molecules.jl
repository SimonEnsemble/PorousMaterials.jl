# Molecules are composed of Lennard-Jones spheres and point charges.
"A Lennard-Jones sphere"
struct LJsphere
    "atom name corresponding to force field"
    atom::Symbol
    "Cartesian coordinates (units: A)"
    x::Array{Float64, 1}
end

"A point charge"
struct Charge
    "charge magnitude (electrons)"
    q::Float64
    "Cartesian coordinates (units: A)"
    x::Array{Float64, 1}
end

"""
Data structure for a molecule/adsorbate.

# Attributes
- `species::Symbol`: Species of molecule, e.g. `CO2` or `ethane`
- `ljspheres::Array{LJspheres, 1}`: array of Lennard-Jones spheres comprising the molecule
- `charges::Array{Charges, 1}`: array of point charges comprising the molecule
"""
struct Molecule
    species::Symbol
    ljspheres::Array{LJsphere, 1}
    charges::Array{Charge, 1}
end

function Base.show(io::IO, molecule::Molecule)
    println(io, "Molecule species: ", molecule.species)
    if length(molecule.ljspheres) > 0
        println(io, "Lennard-Jones spheres: ")
        for ljsphere in molecule.ljspheres
            @printf(io, "\tatom = %s, x = [%.3f, %.3f, %.3f]\n", ljsphere.atom,
                    ljsphere.x[1], ljsphere.x[2], ljsphere.x[3])
        end
    end
    if length(molecule.charges) > 0
        println(io, "Point charges: ")
        for charge in molecule.charges
            @printf(io, "\tq = %f, x = [%.3f, %.3f, %.3f]\n", charge.q, 
                    charge.x[1], charge.x[2], charge.x[3])
        end
    end
end
