using Documenter, PorousMaterials

    makedocs(
        modules = [PorousMaterials],
        sitename = "PorousMaterials.jl",
        pages = [
                "Home" => "index.md"
                "Functions" => "functions.md"
#=
                ["functions.md"
                "Box.jl" => "box.md"
                "Crystal.jl" => "crystal.md"
                "ElectrostaticsEnergetics.jl" => "electrostatic_energetics.md"
                "Energetics_Util.jl" => "energetics_util.md"
                "EOS.jl" => "eos.md"
                "Forcefield.jl" => "forcefield.md"
                "GCMC.jl" => "gcmc.md"
                "Grid.jl" => "grid.md"
                "Henry.jl" => "henry.md"
                "Matter.jl" => "matter.md"
                "MChelpers.jl" => "mchelpers.md"
                "Misc.jl" => "misc.md"
                "Molecules.jl" => "molecules.md"
                "NearestImage.jl" => "nearest_image.md"
                "VdWEnergetics.jl" => "vdw_energetics.md"
                ]
=#
                # Any other pages go here
                ]
    )

    deploydocs(
        repo = "https://github.com/SimonEnsemble/PorousMaterials.jl",
        # This is a link to the main repo and the master branch
        target = "build",
        julia = "0.6",
        osname = "linux",
        deps = nothing, # These are dependencies for the site, not the package
        make = nothing
    )
