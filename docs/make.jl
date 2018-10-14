using Documenter, PorousMaterials

    makedocs(
        root = joinpath(dirname(pathof(PorousMaterials)), "..", "docs"),
        modules = [PorousMaterials],
        sitename = "PorousMaterials.jl",
        clean = true,
        pages = [
                "index.md"
                "guides/input_files.md"

                # Matter
                "manual/matter/docs.md"
                "manual/matter/about.md"

                # Boxes, Crystals, and Grids
                "manual/boxes_crystals_grids/docs.md"
                "manual/boxes_crystals_grids/about.md"

                # Molecules
                "manual/molecules/docs.md"
                "manual/molecules/about.md"

                # Atomic Interactions
                "manual/atomic_interactions/docs.md"
                "manual/atomic_interactions/about.md"

                # MOF Simulations
                "manual/mof_simulations/docs.md"
                "manual/mof_simulations/about.md"

                # Other
                "manual/other/docs.md"
                "manual/other/about.md"

                "guides/faq.md"
                "guides/help_wanted.md"
                ]
    )

    deploydocs(
        repo = "github.com/SimonEnsemble/PorousMaterials.jl.git",
        # This is a link to the main repo and the master branch
        target = "build",
        julia = "1.0",
        osname = "linux",
        deps = Deps.pip3("mkdocs", "mkdocs-windmill", "pymdown-extensions") # These are dependencies for the site, not the package
    )
