using Documenter
using ConstrainedDynamics

makedocs(;
    modules = [ConstrainedDynamics],
    format = Documenter.HTML(
        canonical = "https://janbruedigam.github.io/ConstrainedDynamics.jl/latest/",
        assets = ["assets/favicon.ico"],
    ),
    pages = [
        "Home" => "index.md",
        "Library" => "library/library.md",
        "Body" => "library/body.md",
        "Constraint" => "library/constraint.md",
        "Mechanism" => "library/mechanism.md",
        "Interface" => "library/interface.md",
        "Simulation" => "library/simulation.md",
        "State" => "library/state.md",
        ],
    sitename = "ConstrainedDynamics.jl",
)

deploydocs(; repo = "github.com/janbruedigam/ConstrainedDynamics.jl")