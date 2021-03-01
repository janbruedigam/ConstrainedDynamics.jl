using Documenter
using ConstrainedDynamics

makedocs(;
    modules = [ConstrainedDynamics],
    format = Documenter.HTML(
        canonical = "https://www.ei.tum.de/itr/software/ConstrainedDynamics.jl/stable/",
        assets = ["assets/favicon.ico"],
    ),
    pages = [
        "Home" => "index.md",
        "Library" => "pages/library.md",
        "Body" => "pages/body.md",
        "Constraint" => "pages/constraint.md",
        "Mechanism" => "pages/mechanism.md",
        "Interface" => "pages/interface.md",
        "Simulation" => "pages/simulation.md",
        "State" => "pages/state.md",
        ],
    sitename = "ConstrainedDynamics.jl",
)

deploydocs(; repo = "github.com/janbruedigam/ConstrainedDynamics.jl")