using Documenter
using PorousMaterials

makedocs(
    root = joinpath(dirname(pathof(PorousMaterials)), "..", "docs"),
    modules = [PorousMaterials],
    sitename = "PorousMaterials.jl",
    clean = true,
    pages = [
            "PorousMaterials" => "index.md",
            "paths" => "paths.md",
            "matter" => "matter.md",
            "boxes" => "box.md",
            "crystals" => "crystal.md",
            "molecules" => "molecule.md",
            "computing distances" => "distance.md",
            "forcefields" => "force_field.md",
            "equations of state" => "eos.md",
            "Henry coefficients" => "henry.md",
            "grand-canonical Monte Carlo simulations" => "gcmc.md",
            "energy minimum" => "energy_min.md",
            "grids" => "grid.md"
            ],
    format = Documenter.HTML(assets = ["assets/flux.css"])
)

deploydocs(repo="github.com/SimonEnsemble/PorousMaterials.jl.git")
