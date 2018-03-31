# Molecules are composed of Lennard-Jones spheres and point charges.
struct LennardJonesSphere
    "atom name corresponding to force field"
    atom::Symbol
    "Cartesian coordinates (units: A)"
    x::Array{Float64, 1}
end

"A point charge"
struct PointCharge
    "charge magnitude (electrons)"
    q::Float64
    "Cartesian coordinates (units: A)"
    x::Array{Float64, 1}
end

"""
Data structure for a molecule/adsorbate.

# Attributes
- `species::Symbol`: Species of molecule, e.g. `CO2` or `ethane`
- `ljspheres::Array{LennardJonesSphere, 1}`: array of Lennard-Jones spheres comprising the molecule
- `charges::Array{PointCharges, 1}`: array of point charges comprising the molecule
"""
struct Molecule
    species::Symbol
    ljspheres::Array{LennardJonesSphere, 1}
    charges::Array{PointCharge, 1}
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

"""
    molecule = read_molecule_file("Xe")

Reads molecule files in the directory `PorousMaterials.PATH_TO_DATA * "/molecule/" * species * "/"`.
"""
function read_molecule_file(species::String)
    if ! isdir(PATH_TO_DATA * "molecules/" * species)
        error(@sprintf("No directory created for %s in %s\n", species,
                       PATH_TO_DATA * "molecules/"))
    end
    
    # Read in Lennard Jones spheres
    ljspheresfilename = PATH_TO_DATA * "molecules/" * species * "/lennard_jones_spheres.csv"
    if ! isfile(ljspheresfilename)
        error(@sprintf("No file %s exists. Even if there are no Lennard Jones spheres in your molecule, include a .csv file with the proper headers but no rows.",
                       ljspheresfilename))
    end
    df_lj = CSV.read(ljspheresfilename)
    
    ljspheres = LennardJonesSphere[]
    for row in eachrow(df_lj)
        x = [row[:x], row[:y], row[:z]]
        push!(ljspheres, LennardJonesSphere(Symbol(row[:atom]), x))
    end

    # Read in point charges
    chargesfilename = PATH_TO_DATA * "molecules/" * species * "/point_charges.csv"
    if ! isfile(chargesfilename)
        error(@sprintf("No file %s exists. Even if there are no point charges in your molecule, include a .csv file with the proper headers but no rows.",
                       chargesfilename))
    end
    df_c = CSV.read(chargesfilename)

    charges = PointCharge[]
    for row in eachrow(df_c)
        x = [row[:x], row[:y], row[:z]]
        push!(charges, PointCharge(row[:q], x))
    end
    
    return Molecule(Symbol(species), ljspheres, charges)
end
