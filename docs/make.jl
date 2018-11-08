using Documenter
using PorousMaterials

    makedocs(
        root = joinpath(dirname(pathof(PorousMaterials)), "..", "docs"),
        modules = [PorousMaterials],
        sitename = "PorousMaterials.jl",
        clean = true,
        pages = [
                "PorousMaterials" => "index.md",
                "Loading Data" => "guides/input_files.md",

                "Manual" => [
                    "Matter" =>
                        "manual/matter.md",

                    "Boxes, Crystals, and Grids" => 
                        "manual/boxes_crystals_grids.md",

                    "Molecules" => 
                        "manual/molecules.md",

                    "Atomic Interactions" => 
                        "manual/atomic_interactions.md",

                    "Molecular Simulations" => 
                        "manual/mof_simulations.md",

                    "Other" => 
                        "manual/other.md",
                    ],

                "FAQ" => "guides/faq.md",
                "Help Wanted" => "guides/help_wanted.md",
                ]
    )

    deploydocs(
        repo = "github.com/SimonEnsemble/PorousMaterials.jl.git",
        # This is a link to the main repo and the master branch
        # target = "build",
        julia = "1.0",
        osname = "linux",
        deps = Deps.pip("mkdocs", "mkdocs-windmill", "pymdown-extensions") # These are dependencies for the site, not the package
    )
