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
        repo = "https://github.com/SimonEnsemble/PorousMaterials.jl",
        # This is a link to the main repo and the master branch
        target = "build",
        julia = "0.6",
        osname = "linux",
        deps = nothing, # These are dependencies for the site, not the package
        make = nothing
    )
