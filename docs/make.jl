using Documenter, PorousMaterials

    makedocs(
        modules = [PorousMaterials],
        sitename = "PorousMaterials.jl",
        pages = [
                "Home" => "index.md"
                "Functions" => "functions.md"
                # Any other pages go here
                ]
    )

    deploydocs(
        repo = "github.com/SimonEnsemble/PorousMaterials.jl.git",
        # This is a link to the main repo and the master branch
        target = "build",
        julia = "1.0",
        osname = "linux",
        deps = Deps.pip("mkdocs"), # These are dependencies for the site, not the package
        make = nothing
    )
