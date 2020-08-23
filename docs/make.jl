using Documenter
using PorousMaterials

makedocs(
    root = joinpath(dirname(pathof(PorousMaterials)), "..", "docs"),
    modules = [PorousMaterials],
    sitename = "PorousMaterials.jl",
    clean = true,
    pages = [
            "PorousMaterials" => "index.md",
            "Matter" => "matter.md"
            ],
    format = Documenter.HTML(assets = ["assets/flux.css"])
)

 # deploydocs(
 #     repo = "github.com/SimonEnsemble/PorousMaterials.jl.git",
 #     push_preview=false,
 #     deps = Deps.pip("mkdocs", "mkdocs-material", "pymdown-extensions") # These are dependencies for the site, not the package
 # )
