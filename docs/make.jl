using Documenter
import PorousMaterials

    makedocs(
        root = joinpath(dirname(pathof(PorousMaterials)), "..", "docs"),
        modules = [PorousMaterials],
        sitename = "PorousMaterials.jl",
        clean = true,
        pages = [
                "PorousMaterials" => "index.md",
                "Loading Data" => "guides/input_files.md",

                "Manual" => [
                    "Matter" => [
                        "Matter" => "manual/matter/docs.md",
                        "Matter Documentation" => "manual/matter/about.md"
                        ],

                    "Boxes, Crystals, and Grids" => [
                        "Boxes, Crystals, and Grids" => "manual/boxes_crystals_grids/docs.md",
                        "Boxes, Crystals, and Grids Documentation" => "manual/boxes_crystals_grids/about.md"
                        ],

                    "Molecules" => [
                        "Molecules" => "manual/molecules/docs.md",
                        "Molecules Documentation" => "manual/molecules/about.md"
                        ],

                    "Atomic Interactions" => [
                        "Atomic Interactions" => "manual/atomic_interactions/docs.md",
                        "Atomic Interactions" => "manual/atomic_interactions/about.md"
                        ],

                    "Molecular Simulations" => [ 
                        "Molecular Simulations" => "manual/mof_simulations/docs.md",
                        "Molecular Simulations Documentation" => "manual/mof_simulations/about.md"
                        ],

                    "Other" => [
                        "Other" => "manual/other/docs.md",
                        "Other Documentation" => "manual/other/about.md"
                        ]
                    ],

                "FAQ" => "guides/faq.md",
                "Help Wanted" => "guides/help_wanted.md",
                ]
    )

    deploydocs(
        repo = "github.com/SimonEnsemble/PorousMaterials.jl.git",
        # This is a link to the main repo and the master branch
        target = "build",
        julia = "1.0",
        osname = "linux",
        deps = Deps.pip("mkdocs", "mkdocs-windmill", "pymdown-extensions") # These are dependencies for the site, not the package
    )
